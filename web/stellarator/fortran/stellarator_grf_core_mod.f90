! Fortran twin of stellarator_grf.m: full-surface GRF on the stellarator.
!
! Same structure as the .m: analytic interior-harmonic density from
! torus_surf_den, FMM far field over all surface nodes (self term dropped by
! the FMM), then per close panel the rrq correction minus the naive block,
! added onto the FMM value.  lap3dsdlpmat_r64 returns zero for a coincident
! pair, which matches the self-masking the .m does by index.  Geometry comes
! from utils/stellarator_geo_mex.c through bind(C) -- nothing read from disk,
! MATLAB not involved.
!
! The far field needs the FMM3D Fortran library (integer*8 build):
!   FMM3D_DIR ?= $(HOME)/git/FMM3D      (lib-static/libfmm3d.a)
!
! Case table as in the .m:
!   case:    1     2     3     4     5     6     7
!   order:   4     6     8    10    12    14    16
!   mp:     24    24    36    36    48    60    72
!   np:     72    72   108   108   144   180   216
!   tol:  1e-4  1e-6  1e-8  1e-9 1e-11 1e-14 1e-14
!
! Reference table from the .m (pre edge-ladder work):
!   case 1: order= 4  GRF max rel err = 1.793e-04
!   case 2: order= 6                    2.084e-05
!   case 3: order= 8                    2.711e-07
!   case 4: order=10                    3.704e-08
!   case 5: order=12                    8.510e-10
!   case 6: order=14                    8.679e-12
!   case 7: order=16                    2.416e-12
!
! Build:  make sgrf   (from test/)
! Run:    ./build/stellarator_grf [ncases] [isimd] [ichart] [iadap]

module stellarator_grf_core_mod
  use iso_c_binding,         only: c_long_long, c_double, c_int
  use quatapproximation_mod, only: r64
  use linequaaadrature_mod,  only: gauss_r64
  use koorn_geom_mod,        only: get_vioreanu_nodes, get_vioreanu_wts, &
       koorn_vals2coefs_coefs2vals, koorn_pols_batch_r64
#ifdef BIESOLVER_WASM_SCALAR_ONLY
  use w7x_modes_mod,         only: W7X_NMODE, W7X_NFP, load_w7x_modes
#else
  use stellarator_mesh_mod,  only: W7X_NMODE, W7X_NFP, &
       get_w7x_modes_r64, stellarator_mesh_init_r64, &
       create_stellarator_tri_mesh_r64, stellarator_tri_uv2x_r64
#endif
  use lap3d_close_mod,       only: rrq_r64
  use lap3d_mod,             only: lap3dsdlpmat_r64
  use patch_refine_mod,      only: PR_ADAPTIVE, PR_SQUARE
#ifdef BIESOLVER_WASM_PROGRESS
  use stellarator_progress_mod, only: stellarator_progress_state, &
       stellarator_progress_step, stellarator_progress_reset
#endif
#ifndef BIESOLVER_WASM_SCALAR_ONLY
  use omp_lib,               only: omp_get_max_threads, omp_get_wtime
#endif
  implicit none
  private
  integer(8), parameter, public :: STELLARATOR_TARGET_BLOCK = 64_8

#ifdef BIESOLVER_WASM_PROGRESS
  ! Stage codes are a stable part of the browser progress contract.
  integer(c_int), parameter :: PROGRESS_GEOMETRY_BEGIN   = 1_c_int
  integer(c_int), parameter :: PROGRESS_GEOMETRY_READY   = 2_c_int
  integer(c_int), parameter :: PROGRESS_FAR_BEGIN        = 3_c_int
  integer(c_int), parameter :: PROGRESS_FAR_STEP         = 4_c_int
  integer(c_int), parameter :: PROGRESS_COUNT_BEGIN      = 5_c_int
  integer(c_int), parameter :: PROGRESS_COUNT_STEP       = 6_c_int
  integer(c_int), parameter :: PROGRESS_RRQ_BEGIN        = 7_c_int
  integer(c_int), parameter :: PROGRESS_RRQ_STEP         = 8_c_int
  integer(c_int), parameter :: PROGRESS_SCATTER          = 9_c_int
  integer(c_int), parameter :: PROGRESS_RENDER_BEGIN     = 10_c_int
  integer(c_int), parameter :: PROGRESS_RESULT           = 11_c_int
#endif
  real(r64), parameter :: STELLARATOR_LAM = 0.5_r64

  type, public :: stellarator_case_config
    integer(8) :: mp = 24_8, np = 72_8, order = 4_8
    integer(8) :: isimd = 0_8, ichart = 1_8
    logical :: use_fmm = .true., use_w7x = .false.
    ! W7-X curvature refinement tolerance; the mesh criterion runs at the
    ! solve order.  Ignored for the built-in surface.
    real(r64) :: restol = 1.0e-1_r64
    real(r64) :: fmm_eps = 1.0e-3_r64
  end type stellarator_case_config

  type, public :: stellarator_case_result
    integer(8) :: ntri = 0_8, nsrc = 0_8, nrender = 0_8
    real(r64) :: grf_error = 0.0_r64
    real(r64), allocatable :: sx(:,:), snx(:,:), sw(:), ub(:), ubn(:), u(:)
    real(r64), allocatable :: render_xyz(:,:), render_log_error(:)
    integer(8), allocatable :: render_triangles(:,:)
  end type stellarator_case_result

  public :: stellarator_run_case, stellarator_result_clear

  interface
    function stellarator_geo_ntri(mp, np, p, x0) &
             bind(C, name='stellarator_geo_ntri')
      import :: c_long_long, c_double
      integer(c_long_long)             :: stellarator_geo_ntri
      integer(c_long_long), intent(in) :: mp, np, p
      real(c_double),       intent(in) :: x0(*)
    end function stellarator_geo_ntri

    subroutine stellarator_geo_uv2x(mp, np, p, x0, itri, nuv, uv, x) &
               bind(C, name='stellarator_geo_uv2x')
      import :: c_long_long, c_double
      integer(c_long_long), intent(in)  :: mp, np, p, itri, nuv
      real(c_double),       intent(in)  :: x0(*), uv(*)
      real(c_double),       intent(out) :: x(*)
    end subroutine stellarator_geo_uv2x

    subroutine stellarator_geo(mp, np, p, x0, D, uvs, wts, &
                               sx, snx, sw, rts, rps) &
               bind(C, name='stellarator_geo')
      import :: c_long_long, c_double
      integer(c_long_long), intent(in)  :: mp, np, p
      real(c_double),       intent(in)  :: x0(*), D(*), uvs(*), wts(*)
      real(c_double),       intent(out) :: sx(*), snx(*), sw(*), rts(*), rps(*)
    end subroutine stellarator_geo

#ifdef BIESOLVER_C_BACKEND_ROW_MAJOR
    function biesolver_charts_rowmajor(mp, np, p, x0, D, nfp, nmode, mn, rc, zs, &
                                       restol, cap, cl) &
             bind(C, name='biesolver_charts_rowmajor')
      import :: c_long_long, c_double
      integer(c_long_long)             :: biesolver_charts_rowmajor
      integer(c_long_long), intent(in) :: mp, np, p, nfp, nmode, cap
      real(c_double), intent(in)        :: x0(*), D(*), mn(*), rc(*), zs(*), restol
      real(c_double), intent(out)       :: cl(*)
    end function biesolver_charts_rowmajor

    function biesolver_geo_charts_rowmajor(cl, nchart, mp, np, p, x0, D, uvs, wts, &
                                            nfp, nmode, mn, rc, zs, &
                                            sx, snx, sw, rts, rps) result(ierr) &
             bind(C, name='biesolver_geo_charts_rowmajor')
      import :: c_long_long, c_double, c_int
      integer(c_int)                    :: ierr
      integer(c_long_long), intent(in)  :: nchart, mp, np, p, nfp, nmode
      real(c_double), intent(in)         :: cl(*), x0(*), D(*), uvs(*), wts(*)
      real(c_double), intent(in)         :: mn(*), rc(*), zs(*)
      real(c_double), intent(out)        :: sx(*), snx(*), sw(*), rts(*), rps(*)
    end function biesolver_geo_charts_rowmajor

    function biesolver_uv2x_charts_rowmajor(cl, nchart, mp, np, p, x0, itri, nuv, uv, &
                                             nfp, nmode, mn, rc, zs, x) result(ierr) &
             bind(C, name='biesolver_uv2x_charts_rowmajor')
      import :: c_long_long, c_double, c_int
      integer(c_int)                    :: ierr
      integer(c_long_long), intent(in)  :: nchart, mp, np, p, itri, nuv, nfp, nmode
      real(c_double), intent(in)         :: cl(*), x0(*), uv(*), mn(*), rc(*), zs(*)
      real(c_double), intent(out)        :: x(*)
    end function biesolver_uv2x_charts_rowmajor

    function biesolver_geo_ntri_rowmajor(mp, np, p, x0) &
             bind(C, name='biesolver_geo_ntri_rowmajor')
      import :: c_long_long, c_double
      integer(c_long_long)             :: biesolver_geo_ntri_rowmajor
      integer(c_long_long), intent(in) :: mp, np, p
      real(c_double), intent(in)        :: x0(*)
    end function biesolver_geo_ntri_rowmajor

    function biesolver_geo_rowmajor(ntri, mp, np, p, x0, D, uvs, wts, &
                                    sx, snx, sw, rts, rps) result(ierr) &
             bind(C, name='biesolver_geo_rowmajor')
      import :: c_long_long, c_double, c_int
      integer(c_int)                    :: ierr
      integer(c_long_long), intent(in)  :: ntri, mp, np, p
      real(c_double), intent(in)         :: x0(*), D(*), uvs(*), wts(*)
      real(c_double), intent(out)        :: sx(*), snx(*), sw(*), rts(*), rps(*)
    end function biesolver_geo_rowmajor

    function biesolver_geo_uv2x_rowmajor(mp, np, p, x0, itri, nuv, uv, x) &
             result(ierr) &
             bind(C, name='biesolver_geo_uv2x_rowmajor')
      import :: c_long_long, c_double, c_int
      integer(c_int)                    :: ierr
      integer(c_long_long), intent(in)  :: mp, np, p, itri, nuv
      real(c_double), intent(in)         :: x0(*), uv(*)
      real(c_double), intent(out)        :: x(*)
    end function biesolver_geo_uv2x_rowmajor
#else
    function stellarator_charts(mp, np, p, x0, D, nfp, nmode, mn, rc, zs, &
                                restol, cap, cl) &
             bind(C, name='stellarator_charts')
      import :: c_long_long, c_double
      integer(c_long_long)             :: stellarator_charts
      integer(c_long_long), intent(in) :: mp, np, p, nfp, nmode, cap
      real(c_double), intent(in)        :: x0(*), D(*), mn(*), rc(*), zs(*), restol
      real(c_double), intent(out)       :: cl(*)
    end function stellarator_charts

    subroutine stellarator_geo_charts(cl, nchart, mp, np, p, x0, D, uvs, wts, &
                                      nfp, nmode, mn, rc, zs, &
                                      sx, snx, sw, rts, rps) &
               bind(C, name='stellarator_geo_charts')
      import :: c_long_long, c_double
      integer(c_long_long), intent(in)  :: nchart, mp, np, p, nfp, nmode
      real(c_double), intent(in)        :: cl(*), x0(*), D(*), uvs(*), wts(*)
      real(c_double), intent(in)        :: mn(*), rc(*), zs(*)
      real(c_double), intent(out)       :: sx(*), snx(*), sw(*), rts(*), rps(*)
    end subroutine stellarator_geo_charts

    subroutine stellarator_uv2x_charts(cl, nchart, mp, np, p, x0, itri, nuv, uv, &
                                       nfp, nmode, mn, rc, zs, x) &
               bind(C, name='stellarator_uv2x_charts')
      import :: c_long_long, c_double
      integer(c_long_long), intent(in)  :: nchart, mp, np, p, itri, nuv, nfp, nmode
      real(c_double), intent(in)        :: cl(*), x0(*), uv(*), mn(*), rc(*), zs(*)
      real(c_double), intent(out)       :: x(*)
    end subroutine stellarator_uv2x_charts
#endif

#ifndef BIESOLVER_WASM_SCALAR_ONLY
    subroutine lfmm3d_t_cd_p(eps, nsource, source, charge, dipvec, &
                             ntarg, targ, pottarg, ier) &
        bind(C, name='lfmm3d_t_cd_p_')
      integer(8), intent(in) :: nsource, ntarg
      integer(8), intent(out) :: ier
      real(8), intent(in) :: eps, source(3,*), charge(*), dipvec(3,*)
      real(8), intent(in) :: targ(3,*)
      real(8), intent(out) :: pottarg(*)
    end subroutine lfmm3d_t_cd_p
#endif

#ifdef BIESOLVER_WASM_FMM3D
    subroutine biesolver_lfmm3d_t_cd_p_rowmajor(eps, nsource, source, &
        charge, dipvec, ntarg, targ, pottarg, ier, adapter_status, elapsed) &
        bind(C, name='biesolver_lfmm3d_t_cd_p_rowmajor')
      import :: c_double, c_long_long
      real(c_double), value :: eps
      integer(c_long_long), value :: nsource, ntarg
      real(c_double), intent(in) :: source(3,nsource), charge(nsource)
      real(c_double), intent(in) :: dipvec(3,nsource), targ(3,ntarg)
      real(c_double), intent(out) :: pottarg(ntarg)
      integer(c_long_long), intent(out) :: ier, adapter_status
      real(c_double), intent(out) :: elapsed
    end subroutine biesolver_lfmm3d_t_cd_p_rowmajor
#endif

#ifdef BIESOLVER_WASM_SCALAR_ONLY
    ! By-value C ABI: unlike the geometry entry points above, this takes its
    ! argument by value, matching the plain C declaration in
    ! native/wasm_memory_limits.h.
    function biesolver_wasm_memory_preflight(nsrc) result(ok) &
        bind(C, name='biesolver_wasm_memory_preflight')
      import :: c_int, c_long_long
      integer(c_long_long), value :: nsrc
      integer(c_int) :: ok
    end function biesolver_wasm_memory_preflight
#endif

#ifdef BIESOLVER_WASM_PROGRESS
    subroutine biesolver_progress(stage, current, total, aux0, aux1, value) &
        bind(C, name='biesolver_progress')
      import :: c_int, c_long_long, c_double
      integer(c_int), value :: stage
      integer(c_long_long), value :: current, total, aux0, aux1
      real(c_double), value :: value
    end subroutine biesolver_progress
#endif
  end interface

contains

  subroutine stellarator_run_case(cfg, tgl, wgl, Dgl, w_bclag, Legmat, &
      umatr, vmatr, result, status, t_fmm_out, t_close_out, timeinfo_out)
  type(stellarator_case_config), intent(in) :: cfg
  type(stellarator_case_result), intent(inout) :: result
  integer(8), intent(out) :: status
  real(r64), intent(out), optional :: t_fmm_out, t_close_out, timeinfo_out(20)
  real(r64), intent(in) :: tgl(cfg%order+2_8), wgl(cfg%order+2_8)
  real(r64), intent(in) :: Dgl(cfg%order+2_8,cfg%order+2_8)
  real(r64), intent(in) :: w_bclag(cfg%order+2_8)
  real(r64), intent(in) :: Legmat(cfg%order+2_8,cfg%order+2_8)
  real(r64), intent(in) :: umatr(cfg%order*(cfg%order+1_8)/2_8, &
                                 cfg%order*(cfg%order+1_8)/2_8)
  real(r64), intent(in) :: vmatr(cfg%order*(cfg%order+1_8)/2_8, &
                                 cfg%order*(cfg%order+1_8)/2_8)

  integer(8) :: mp, np, nterms, hdim, nquad, orderff, ntri, nsrc
  integer(8) :: isimd, ichart
  integer(8) :: k, i, j, ier, adapter_status, nbd, sbdnp, nnz, mtcmax
  integer(8) :: cap, nchart, ib, mb
  integer(c_int) :: geometry_ierr
  real(r64)  :: distff, timeinfo(20), err, umax, ubnmax
  real(r64)  :: pi, tt, th1, th2, thi1
  real(r64)  :: t_fmm, t_close, t_geo, w0, w1p, w2p, w3p
  integer(8) :: c0, c1, crate
  logical(kind(.true.)) :: exterior

  real(r64), allocatable :: tgl_geo(:), wgl_geo(:), Dgl_geo(:,:)
  real(r64), allocatable :: uvs(:,:), wts(:)
  real(r64), allocatable :: uvbd(:,:), tpan(:)
  real(r64), allocatable :: sx(:,:), snx(:,:), sw(:), rts(:,:), rps(:,:)
  real(r64), allocatable :: ub(:), ubn(:), u(:), charge(:), dipvec(:,:)
#ifdef BIESOLVER_WASM_SCALAR_ONLY
  real(r64), allocatable :: cl(:), mn(:), rc(:), zs(:)
#else
  real(r64), allocatable :: charts(:,:), rc(:), zs(:)
  integer(8), allocatable :: mn(:,:)
  integer(8) :: geometry_nfp, geometry_nmode
#endif
  real(r64), allocatable :: csr_val(:)
  integer(8),allocatable :: mtcs(:), offs(:), csr_idx(:)
#ifdef BIESOLVER_WASM_PROGRESS
  type(stellarator_progress_state) :: direct_progress, count_progress
  type(stellarator_progress_state) :: rrq_progress
  logical    :: emit_progress
  integer(8) :: completed, target_blocks
#endif

  call stellarator_result_clear(result)
  status = 0_8
  if (present(t_fmm_out)) t_fmm_out = 0.0_r64
  if (present(t_close_out)) t_close_out = 0.0_r64
  if (present(timeinfo_out)) timeinfo_out = 0.0_r64
  exterior = .true.
  distff   = 1.4_r64
  pi       = 4.0_r64*atan(1.0_r64)
#ifndef BIESOLVER_WASM_SCALAR_ONLY
  write(*,'(a,i0,a)') 'close loop on ', omp_get_max_threads(), &
      ' threads (timeinfo is cpu-summed across threads)'
#endif

    nterms  = cfg%order
    mp      = cfg%mp
    np      = cfg%np
    isimd   = cfg%isimd
    ichart  = cfg%ichart
    hdim    = nterms*(nterms+1_8)/2_8
    nquad   = nterms + 2_8
    orderff = nterms + 2_8
    timeinfo = 0.0_r64

    ! ---- geometry ----
#ifdef BIESOLVER_WASM_PROGRESS
    call biesolver_progress(PROGRESS_GEOMETRY_BEGIN, 0_c_long_long, &
      0_c_long_long, int(merge(1,0,cfg%use_w7x),c_long_long), &
      0_c_long_long, 0.0_c_double)
#endif
#ifndef BIESOLVER_WASM_SCALAR_ONLY
    call system_clock(c0, crate)
#else
    t_geo = 0.0_r64
#endif
    allocate(tgl_geo(nterms), wgl_geo(nterms), Dgl_geo(nterms,nterms))
    allocate(uvs(2,hdim), wts(hdim))
    call gauss_r64(nterms, tgl_geo, wgl_geo, Dgl_geo)
    call get_vioreanu_nodes(nterms-1_8, hdim, uvs)
    call get_vioreanu_wts(nterms-1_8, hdim, wts)
#ifdef BIESOLVER_WASM_SCALAR_ONLY
    if (cfg%use_w7x) then
      if (cfg%restol <= 0.0_r64 .or. cfg%restol /= cfg%restol) then
        status = 14_8
        return
      end if
      cap = 200000_8
      allocate(cl(6*cap), mn(2*W7X_NMODE), rc(W7X_NMODE), zs(W7X_NMODE))
      call load_w7x_modes(mn, rc, zs)
#ifdef BIESOLVER_C_BACKEND_ROW_MAJOR
      nchart = biesolver_charts_rowmajor(int(mp,c_long_long), int(np,c_long_long), &
                   int(nterms,c_long_long), tgl_geo, Dgl_geo, int(W7X_NFP,c_long_long), &
                   int(W7X_NMODE,c_long_long), mn, rc, zs, &
                   cfg%restol, int(cap,c_long_long), cl)
#else
      nchart = stellarator_charts(int(mp,c_long_long), int(np,c_long_long), &
                   int(nterms,c_long_long), tgl_geo, Dgl_geo, int(W7X_NFP,c_long_long), &
                   int(W7X_NMODE,c_long_long), mn, rc, zs, &
                   cfg%restol, int(cap,c_long_long), cl)
#endif
      ! The builder returns -1 when the refinement would exceed cap; a clean
      ! error tells the caller to loosen restol rather than truncating.
      if (nchart <= 0_8 .or. nchart > cap) then
        status = 10_8
        return
      end if
      ntri = 2_8*nchart
      nsrc = ntri*hdim
#ifdef BIESOLVER_WASM_SCALAR_ONLY
      ! A 64-by-nsrc direct block above the allocator's single-request limit
      ! would abort inside the allocator; reject it before allocating instead.
      if (.not. cfg%use_fmm) then
        if (biesolver_wasm_memory_preflight(int(nsrc,c_long_long)) == 0_c_int) then
          status = 15_8
          return
        end if
      end if
#endif
      allocate(sx(3,nsrc), snx(3,nsrc), sw(nsrc), rts(3,nsrc), rps(3,nsrc))
#ifdef BIESOLVER_C_BACKEND_ROW_MAJOR
      geometry_ierr = biesolver_geo_charts_rowmajor(cl, int(nchart,c_long_long), &
                   int(mp,c_long_long), int(np,c_long_long), int(nterms,c_long_long), &
                   tgl_geo, Dgl_geo, uvs, wts, int(W7X_NFP,c_long_long), &
                   int(W7X_NMODE,c_long_long), mn, rc, zs, &
                   sx, snx, sw, rts, rps)
      if (geometry_ierr /= 0_c_int) then
        status = 12_8
        return
      end if
#else
      call stellarator_geo_charts(cl, int(nchart,c_long_long), &
                   int(mp,c_long_long), int(np,c_long_long), int(nterms,c_long_long), &
                   tgl_geo, Dgl_geo, uvs, wts, int(W7X_NFP,c_long_long), &
                   int(W7X_NMODE,c_long_long), mn, rc, zs, &
                   sx, snx, sw, rts, rps)
#endif
    else
      nchart = 0_8
#ifdef BIESOLVER_C_BACKEND_ROW_MAJOR
      ntri = biesolver_geo_ntri_rowmajor(int(mp,c_long_long), &
                 int(np,c_long_long), int(nterms,c_long_long), tgl_geo)
      if (ntri <= 0_8) then
        status = 11_8
        return
      end if
      nsrc = ntri*hdim
#ifdef BIESOLVER_WASM_SCALAR_ONLY
      ! A 64-by-nsrc direct block above the allocator's single-request limit
      ! would abort inside the allocator; reject it before allocating instead.
      if (.not. cfg%use_fmm) then
        if (biesolver_wasm_memory_preflight(int(nsrc,c_long_long)) == 0_c_int) then
          status = 15_8
          return
        end if
      end if
#endif
      allocate(sx(3,nsrc), snx(3,nsrc), sw(nsrc), rts(3,nsrc), rps(3,nsrc))
      geometry_ierr = biesolver_geo_rowmajor(int(ntri,c_long_long), &
                   int(mp,c_long_long), int(np,c_long_long), &
                   int(nterms,c_long_long), tgl_geo, Dgl_geo, uvs, wts, &
                   sx, snx, sw, rts, rps)
      if (geometry_ierr /= 0_c_int) then
        status = 12_8
        return
      end if
#else
      ntri = stellarator_geo_ntri(int(mp,c_long_long), int(np,c_long_long), &
                                  int(nterms,c_long_long), tgl_geo)
      nsrc = ntri*hdim
#ifdef BIESOLVER_WASM_SCALAR_ONLY
      ! A 64-by-nsrc direct block above the allocator's single-request limit
      ! would abort inside the allocator; reject it before allocating instead.
      if (.not. cfg%use_fmm) then
        if (biesolver_wasm_memory_preflight(int(nsrc,c_long_long)) == 0_c_int) then
          status = 15_8
          return
        end if
      end if
#endif
      allocate(sx(3,nsrc), snx(3,nsrc), sw(nsrc), rts(3,nsrc), rps(3,nsrc))
      call stellarator_geo(int(mp,c_long_long), int(np,c_long_long), &
                           int(nterms,c_long_long), tgl_geo, Dgl_geo, uvs, wts, &
                           sx, snx, sw, rts, rps)
#endif
    end if
#else
    if (cfg%use_w7x) then
      if (cfg%restol <= 0.0_r64 .or. cfg%restol /= cfg%restol) then
        status = 14_8
        return
      end if
      cap = 200000_8
      allocate(charts(6,cap), mn(2,W7X_NMODE), rc(W7X_NMODE), zs(W7X_NMODE))
      call get_w7x_modes_r64(geometry_nfp, geometry_nmode, mn, rc, zs)
    else
      cap = 4_8*mp*np
      geometry_nfp = 1_8
      geometry_nmode = 0_8
      allocate(charts(6,cap), mn(2,0), rc(0), zs(0))
    end if
    call stellarator_mesh_init_r64(mp, np, nterms, geometry_nfp, &
         geometry_nmode, mn, rc, zs, merge(cfg%restol,0.0_r64,cfg%use_w7x), &
         cap, charts, nchart, ntri, ier)
    if (ier /= 0_8) then
      status = merge(10_8,11_8,cfg%use_w7x)
      return
    end if
    nsrc = ntri*hdim
    allocate(sx(3,nsrc), snx(3,nsrc), sw(nsrc), rts(3,nsrc), rps(3,nsrc))
    call create_stellarator_tri_mesh_r64(mp, np, nterms, geometry_nfp, &
         geometry_nmode, mn, rc, zs, nchart, charts(:,1:nchart), ntri, &
         sx, snx, sw, rts, rps, ier)
    if (ier /= 0_8) then
      status = 12_8
      return
    end if
#endif
#ifndef BIESOLVER_WASM_SCALAR_ONLY
    call system_clock(c1);  t_geo = real(c1-c0, r64)/real(crate, r64)
#endif

#ifdef BIESOLVER_WASM_PROGRESS
    call biesolver_progress(PROGRESS_GEOMETRY_READY, int(ntri,c_long_long), &
      int(nsrc,c_long_long), int(nchart,c_long_long), &
      int(nterms,c_long_long), 0.0_c_double)
#endif

    ! ---- GRF density: interior harmonic u, and du/dn ----
    allocate(ub(nsrc), ubn(nsrc))
    do i = 1, nsrc
      ub(i)  = f_harm(sx(:,i))
      ubn(i) = dot_product(snx(:,i), gradf_harm(sx(:,i)))
    end do

    ! ---- far field: SLP(ubn) - DLP(ub) in one FMM, targets = sources ----
    ! kernel 1/(4 pi r), self term dropped by the FMM
    allocate(u(nsrc))
#ifdef BIESOLVER_WASM_SCALAR_ONLY
#ifndef BIESOLVER_WASM_FMM3D
    if (cfg%use_fmm) then
      status = 1_8
      return
    end if
#endif
#endif
    if (cfg%use_fmm) then
      allocate(charge(nsrc), dipvec(3,nsrc))
      charge = sw*ubn
      do i = 1, nsrc
        dipvec(:,i) = -sw(i)*ub(i)*snx(:,i)
      end do
      ier = 0_8
      adapter_status = 0_8
#ifdef BIESOLVER_WASM_PROGRESS
      call biesolver_progress(PROGRESS_FAR_BEGIN, 0_c_long_long, &
        1_c_long_long, 1_c_long_long, 0_c_long_long, 0.0_c_double)
#endif
#ifdef BIESOLVER_WASM_FMM3D
      call biesolver_lfmm3d_t_cd_p_rowmajor(cfg%fmm_eps, nsrc, sx, &
        charge, dipvec, nsrc, sx, u, ier, adapter_status, t_fmm)
#else
#ifndef BIESOLVER_WASM_SCALAR_ONLY
      call system_clock(c0)
      call lfmm3d_t_cd_p(cfg%fmm_eps, nsrc, sx, charge, dipvec, nsrc, sx, u, ier)
      call system_clock(c1);  t_fmm = real(c1-c0, r64)/real(crate, r64)
#endif
#endif
      deallocate(charge, dipvec)
      if (adapter_status /= 0_8) then
        status = 16_8
      else if (ier == 4_8) then
        status = 17_8
      else if (ier == 8_8) then
        status = 18_8
      else if (ier /= 0_8) then
        status = 19_8
      end if
      if (status /= 0_8) return
      do i = 1, nsrc
        if (u(i) /= u(i) .or. abs(u(i)) > 1.0e300_r64) then
          status = 19_8
          return
        end if
      end do
#ifdef BIESOLVER_WASM_PROGRESS
      call biesolver_progress(PROGRESS_FAR_STEP, 1_c_long_long, &
        1_c_long_long, 1_c_long_long, 0_c_long_long, real(t_fmm,c_double))
#endif
    else
#ifndef BIESOLVER_WASM_SCALAR_ONLY
      call system_clock(c0, crate)
#else
      t_fmm = 0.0_r64
#endif
      block
        real(r64), allocatable :: targets(:,:), as_block(:,:), ad_block(:,:)
        allocate(targets(3,STELLARATOR_TARGET_BLOCK))
        allocate(as_block(STELLARATOR_TARGET_BLOCK,nsrc))
        allocate(ad_block(STELLARATOR_TARGET_BLOCK,nsrc))
#ifdef BIESOLVER_WASM_PROGRESS
        target_blocks = (nsrc + STELLARATOR_TARGET_BLOCK - 1_8) / &
                        STELLARATOR_TARGET_BLOCK
        completed = 0_8
        call stellarator_progress_reset(direct_progress)
        call biesolver_progress(PROGRESS_FAR_BEGIN, 0_c_long_long, &
          int(target_blocks,c_long_long), 0_c_long_long, 0_c_long_long, &
          0.0_c_double)
#endif
        do ib = 1, nsrc, STELLARATOR_TARGET_BLOCK
          mb = min(STELLARATOR_TARGET_BLOCK, nsrc-ib+1_8)
          do i = 1, STELLARATOR_TARGET_BLOCK
            do k = 1, 3_8
              targets(k,i) = sx(k,min(ib+i-1_8,nsrc))
            end do
          end do
          call lap3dsdlpmat_r64(STELLARATOR_TARGET_BLOCK, targets, nsrc, &
            sx, snx, sw, as_block, ad_block)
          do i = 1, mb
            u(ib+i-1_8) = 0.0_r64
            do j = 1, nsrc
              u(ib+i-1_8) = u(ib+i-1_8) + &
                as_block(i,j)*ubn(j) - ad_block(i,j)*ub(j)
            end do
          end do
#ifdef BIESOLVER_WASM_PROGRESS
          completed = completed + 1_8
          call stellarator_progress_step(direct_progress, completed, &
            target_blocks, emit_progress)
          if (emit_progress) call biesolver_progress(PROGRESS_FAR_STEP, &
            int(completed,c_long_long), int(target_blocks,c_long_long), &
            0_c_long_long, 0_c_long_long, 0.0_c_double)
#endif
        end do
        deallocate(targets, as_block, ad_block)
      end block
#ifndef BIESOLVER_WASM_SCALAR_ONLY
      call system_clock(c1)
      t_fmm = real(c1-c0, r64)/real(crate, r64)
#endif
    end if

    ! ---- chart boundary parameterization table ----
    sbdnp = 3_8;  nbd = sbdnp*nquad
    allocate(tpan(sbdnp+1), uvbd(2,nbd+3))
    do k = 1, sbdnp+1_8
      tpan(k) = real(k-1, r64)*2.0_r64*pi/real(sbdnp, r64)
    end do
    th1 = 2.0_r64*pi/3.0_r64;  th2 = 4.0_r64*pi/3.0_r64;  thi1 = 1.0_r64/th1
    do k = 1, sbdnp
      do i = 1, nquad
        j  = (k-1_8)*nquad + i
        tt = tpan(k) + 0.5_r64*(1.0_r64 + tgl(i))*(tpan(k+1) - tpan(k))
        uvbd(:,j) = tpar2uv(tt, th1, th2, thi1)
      end do
    end do
    uvbd(:,nbd+1) = [0.0_r64, 0.0_r64]
    uvbd(:,nbd+2) = [1.0_r64, 0.0_r64]
    uvbd(:,nbd+3) = [0.0_r64, 1.0_r64]

    ! ---- close eval: per panel, correct the near targets ----
    ! OpenMP with buffered scatter: (0) parallel near-target count per panel,
    ! (1) parallel fill of each panel's CSR slice, (2) serial scatter in panel
    ! order.  The scatter reproduces the serial loop's summation order into u,
    ! so the printed digits are independent of the thread count.  One serial
    ! warm-up call precedes the parallel region so patch-refinement's shared
    ! per-order reference state is initialized before threaded evaluation.
    allocate(mtcs(ntri), offs(ntri+1_8))
#ifndef BIESOLVER_WASM_SCALAR_ONLY
    call system_clock(c0)
#else
    t_close = 0.0_r64
#endif
    w0 = stellarator_wall_time()

#ifdef BIESOLVER_WASM_PROGRESS
    call stellarator_progress_reset(count_progress)
    call biesolver_progress(PROGRESS_COUNT_BEGIN, 0_c_long_long, &
      int(ntri,c_long_long), 0_c_long_long, 0_c_long_long, 0.0_c_double)
#endif
#ifndef BIESOLVER_WASM_SCALAR_ONLY
    !$omp parallel do schedule(static)
#endif
    do k = 1, ntri
      block
        real(r64) :: qr, qp(3), swk(hdim), sxk(3,hdim)
        integer(8) :: i2, j2, cnt
        qp(1) = 0.0_r64
        qp(2) = 0.0_r64
        qp(3) = 0.0_r64
        do i2 = 1, hdim
          sxk(:,i2) = sx(:, (k-1_8)*hdim + i2)
          swk(i2)   = sw((k-1_8)*hdim + i2)
          qp(1) = qp(1) + sxk(1,i2)
          qp(2) = qp(2) + sxk(2,i2)
          qp(3) = qp(3) + sxk(3,i2)
        end do
        qr    = 1.75_r64*sqrt(sum(swk))
        qp(1) = qp(1)/real(hdim, r64)
        qp(2) = qp(2)/real(hdim, r64)
        qp(3) = qp(3)/real(hdim, r64)
        cnt = 0_8
        do j2 = 1, nsrc
          if ((sx(1,j2)-qp(1))**2 + (sx(2,j2)-qp(2))**2 + &
              (sx(3,j2)-qp(3))**2 < qr**2) cnt = cnt + 1_8
        end do
        mtcs(k) = cnt
      end block
#ifdef BIESOLVER_WASM_PROGRESS
      call stellarator_progress_step(count_progress, k, ntri, emit_progress)
      if (emit_progress) call biesolver_progress(PROGRESS_COUNT_STEP, &
        int(k,c_long_long), int(ntri,c_long_long), 0_c_long_long, &
        0_c_long_long, 0.0_c_double)
#endif
    end do
#ifndef BIESOLVER_WASM_SCALAR_ONLY
    !$omp end parallel do
#endif
    offs(1) = 0_8
    do k = 1, ntri
      offs(k+1_8) = offs(k) + mtcs(k)
    end do
    nnz = offs(ntri+1_8)
    allocate(csr_idx(nnz), csr_val(nnz))
    w1p = stellarator_wall_time()

    ! serial warm-up: build the lazy per-order caches before threading
    block
      real(r64) :: sxk(3,hdim), snxk(3,hdim), swk(hdim)
      real(r64) :: rtsk(3,hdim), rpsk(3,hdim)
      real(r64), allocatable :: xbuf(:,:)
      real(r64) :: tw(3,1), S1w(1,hdim), K1w(1,hdim), I1w(1)
      real(r64), allocatable :: O1w(:,:)
      real(r64) :: tinfo0(20)
      integer(8) :: i2
      allocate(xbuf(3,nbd+3), O1w(4*hdim,1))
      do i2 = 1, hdim
        sxk(:,i2)  = sx(:,i2);   snxk(:,i2) = snx(:,i2);  swk(i2) = sw(i2)
        rtsk(:,i2) = rts(:,i2);  rpsk(:,i2) = rps(:,i2)
      end do
      tw(:,1) = sx(:,1) + 0.01_r64*snx(:,1)
      xbuf = 0.0_r64
#ifdef BIESOLVER_WASM_SCALAR_ONLY
      if (ichart == 1_8) call map_panel_uv_bridge(cfg%use_w7x, cl, nchart, &
        mp, np, nterms, tgl_geo, 1_8, nbd+3_8, uvbd, mn, rc, zs, xbuf, geometry_ierr)
#else
      if (ichart == 1_8) call map_panel_uv_bridge(charts, nchart, mp, np, &
        nterms, 1_8, nbd+3_8, uvbd, mn, rc, zs, xbuf, geometry_ierr)
#endif
#ifdef BIESOLVER_C_BACKEND_ROW_MAJOR
      if (geometry_ierr /= 0_c_int) then
        status = 13_8
        return
      end if
#endif
      do i2 = 1, hdim
        S1w(1,i2) = 0.0_r64
        K1w(1,i2) = 0.0_r64
      end do
      do i2 = 1, 4_8*hdim
        O1w(i2,1) = 0.0_r64
      end do
      I1w(1) = 0.0_r64
      tinfo0 = 0.0_r64
      call rrq_r64(1_8, tw, hdim, sxk, snxk, swk, rtsk, rpsk, &
                   nterms, nquad, orderff, distff, exterior, isimd, &
                   ichart, xbuf(:,1:nbd), xbuf(:,nbd+1:nbd+3), &
                   tgl, wgl, Dgl, w_bclag, Legmat, umatr, vmatr, &
                   S1w, K1w, O1w, I1w, tinfo0)
      deallocate(xbuf, O1w)
    end block

    ! per-thread work buffers, allocated once outside the do loop.  A
    ! per-panel allocate/deallocate here puts every thread on the macOS
    ! large-allocation mutex (visible in `sample` as _pthread_mutex_firstfit
    ! under _omp_fn, with madvise on the free side), which is what kept the
    ! close loop from scaling once the arrays got big.  The buffers are 1-D
    ! and sequence-associated to rrq_r64's explicit-shape dummies: a 2-D
    ! buffer sliced to (1:mtc,1:hdim) would be non-contiguous and copied.
    mtcmax = 0_8
    do k = 1, ntri
      mtcmax = max(mtcmax, mtcs(k))
    end do
#ifdef BIESOLVER_WASM_PROGRESS
    call stellarator_progress_reset(rrq_progress)
    call biesolver_progress(PROGRESS_RRQ_BEGIN, 0_c_long_long, &
      int(ntri,c_long_long), 0_c_long_long, 0_c_long_long, 0.0_c_double)
#endif
#ifndef BIESOLVER_WASM_SCALAR_ONLY
    !$omp parallel default(shared) reduction(+:timeinfo)
#endif
    block
      real(r64), allocatable :: wtcx(:), wAs(:), wAd(:)
      real(r64), allocatable :: wS(:), wK(:), wOm(:), wIa(:)
      allocate(wtcx(3*mtcmax), wAs(mtcmax*hdim), wAd(mtcmax*hdim))
      allocate(wS(mtcmax*hdim), wK(mtcmax*hdim))
      allocate(wOm(4*hdim*mtcmax), wIa(mtcmax))

#ifndef BIESOLVER_WASM_SCALAR_ONLY
    !$omp do schedule(dynamic,4)
#endif
    do k = 1, ntri
      if (mtcs(k) /= 0_8) then
        block
          real(r64) :: qr, qp(3)
          real(r64), allocatable :: xbuf(:,:)
          real(r64) :: sxk(3,hdim), snxk(3,hdim), swk(hdim)
          real(r64) :: rtsk(3,hdim), rpsk(3,hdim)
          integer(8) :: i2, j2, mtc, o0, i0
          allocate(xbuf(3,nbd+3))
          i0 = (k-1_8)*hdim
          qp(1) = 0.0_r64
          qp(2) = 0.0_r64
          qp(3) = 0.0_r64
          do i2 = 1, hdim
            sxk(:,i2)  = sx(:,i0+i2);   snxk(:,i2) = snx(:,i0+i2)
            swk(i2)    = sw(i0+i2)
            rtsk(:,i2) = rts(:,i0+i2);  rpsk(:,i2) = rps(:,i0+i2)
            qp(1) = qp(1) + sxk(1,i2)
            qp(2) = qp(2) + sxk(2,i2)
            qp(3) = qp(3) + sxk(3,i2)
          end do
          qr    = 1.75_r64*sqrt(sum(swk))
          qp(1) = qp(1)/real(hdim, r64)
          qp(2) = qp(2)/real(hdim, r64)
          qp(3) = qp(3)/real(hdim, r64)
          mtc = mtcs(k);  o0 = offs(k)
          j2 = 0_8
          do i2 = 1, nsrc
            if ((sx(1,i2)-qp(1))**2 + (sx(2,i2)-qp(2))**2 + &
                (sx(3,i2)-qp(3))**2 < qr**2) then
              j2 = j2 + 1_8
              if (j2 > mtc) exit
              csr_idx(o0+j2) = i2
            end if
          end do
          if (j2 /= mtc) then
            write(*,'(a,i0,a,i0,a,i0)') 'COUNT MISMATCH panel ', k, &
                ': phase0 ', mtc, '  phase1 ', j2
            error stop 'csr count mismatch'
          end if

          do j2 = 1, mtc
#ifdef BIESOLVER_C_BACKEND_ROW_MAJOR
            do i2 = 1, 3
              wtcx((i2-1_8)*mtc+j2) = sx(i2,csr_idx(o0+j2))
            end do
#else
            wtcx(3*j2-2_8:3*j2) = sx(:, csr_idx(o0+j2))
#endif
          end do
          wS(1:mtc*hdim)     = 0.0_r64
          wK(1:mtc*hdim)     = 0.0_r64
          wOm(1:4*hdim*mtc)  = 0.0_r64
          wIa(1:mtc)         = 0.0_r64
          xbuf = 0.0_r64
#ifdef BIESOLVER_WASM_SCALAR_ONLY
          if (ichart == 1_8) call map_panel_uv_bridge(cfg%use_w7x, cl, nchart, &
            mp, np, nterms, tgl_geo, k, nbd+3_8, uvbd, mn, rc, zs, xbuf, geometry_ierr)
#else
          if (ichart == 1_8) call map_panel_uv_bridge(charts, nchart, mp, np, &
            nterms, k, nbd+3_8, uvbd, mn, rc, zs, xbuf, geometry_ierr)
#endif
#ifdef BIESOLVER_C_BACKEND_ROW_MAJOR
          if (geometry_ierr /= 0_c_int) then
            status = 13_8
            return
          end if
#endif
          call rrq_r64(mtc, wtcx, hdim, sxk, snxk, swk, rtsk, rpsk, &
                       nterms, nquad, orderff, distff, exterior, isimd, &
                       ichart, xbuf(:,1:nbd), xbuf(:,nbd+1:nbd+3), &
                       tgl, wgl, Dgl, w_bclag, Legmat, umatr, vmatr, &
                       wS, wK, wOm, wIa, timeinfo)


          ! naive block; zero at coincident pairs, matching the FMM's dropped
          ! self term and the .m's index masking
          call lap3dsdlpmat_r64(mtc, wtcx, hdim, sxk, snxk, swk, wAs, wAd)


          ! The 1-D buffers are sequence-associated with (mtc,hdim) dummies.
          ! gfortran uses column-major storage; LFortran's C backend lowers
          ! the same explicit-shape arrays in row-major order.
          do j2 = 1, mtc
            csr_val(o0+j2) = 0.0_r64
          end do
          do i2 = 1, hdim
            do j2 = 1, mtc
#ifdef BIESOLVER_C_BACKEND_ROW_MAJOR
              csr_val(o0+j2) = csr_val(o0+j2) &
                + (wS((j2-1_8)*hdim+i2) - wAs((j2-1_8)*hdim+i2))*ubn(i0+i2) &
                - (wK((j2-1_8)*hdim+i2) - wAd((j2-1_8)*hdim+i2))*ub(i0+i2)
#else
              csr_val(o0+j2) = csr_val(o0+j2) &
                + (wS((i2-1_8)*mtc+j2) - wAs((i2-1_8)*mtc+j2))*ubn(i0+i2) &
                - (wK((i2-1_8)*mtc+j2) - wAd((i2-1_8)*mtc+j2))*ub(i0+i2)
#endif
            end do
          end do
          deallocate(xbuf)
        end block
      end if
#ifdef BIESOLVER_WASM_PROGRESS
      call stellarator_progress_step(rrq_progress, k, ntri, emit_progress)
      if (emit_progress) call biesolver_progress(PROGRESS_RRQ_STEP, &
        int(k,c_long_long), int(ntri,c_long_long), 0_c_long_long, &
        0_c_long_long, 0.0_c_double)
#endif
    end do
#ifndef BIESOLVER_WASM_SCALAR_ONLY
    !$omp end do
#endif
      deallocate(wtcx, wAs, wAd, wS, wK, wOm, wIa)
    end block
#ifndef BIESOLVER_WASM_SCALAR_ONLY
    !$omp end parallel
#endif
    w2p = stellarator_wall_time()

#ifdef BIESOLVER_WASM_PROGRESS
    call biesolver_progress(PROGRESS_SCATTER, int(ntri,c_long_long), &
      int(ntri,c_long_long), 0_c_long_long, 0_c_long_long, 0.0_c_double)
#endif
    do k = 1, ntri
      do j = offs(k)+1_8, offs(k+1_8)
        u(csr_idx(j)) = u(csr_idx(j)) + csr_val(j)
      end do
    end do
    w3p = stellarator_wall_time()
#ifndef BIESOLVER_WASM_SCALAR_ONLY
    call system_clock(c1);  t_close = real(c1-c0, r64)/real(crate, r64)
#endif
#ifndef BIESOLVER_WASM_SCALAR_ONLY
    write(*,'(a,f8.3,a,f8.3,a,f8.3,a)') '  phases: ballcount ', w1p-w0, &
        ' s | rrq+fill ', w2p-w1p, ' s | scatter ', w3p-w2p, ' s'
#endif

    umax = 0.0_r64
    ubnmax = 0.0_r64
    do i = 1, nsrc
      if (u(i) /= u(i) .or. ubn(i) /= ubn(i) .or. &
          abs(u(i)) > 1.0e300_r64 .or. abs(ubn(i)) > 1.0e300_r64) then
        status = 30_8
        return
      end if
      umax = max(umax, abs(u(i)))
      ubnmax = max(ubnmax, abs(ubn(i)))
    end do
    if (umax /= umax .or. ubnmax /= ubnmax .or. umax > 1.0e300_r64 .or. &
        ubnmax > 1.0e300_r64 .or. ubnmax <= 0.0_r64) then
      status = 31_8
      return
    end if
    err = umax/ubnmax
    result%ntri = ntri
    result%nsrc = nsrc
    result%grf_error = err
#ifdef BIESOLVER_WASM_PROGRESS
    call biesolver_progress(PROGRESS_RENDER_BEGIN, int(ntri,c_long_long), &
      0_c_long_long, 0_c_long_long, 0_c_long_long, 0.0_c_double)
#endif
#ifdef BIESOLVER_WASM_SCALAR_ONLY
    call build_render_lattice(ubnmax, cfg%use_w7x, result, status, &
      nterms, hdim, ntri, nchart, mp, np, tgl_geo, cl, mn, rc, zs, u, &
      geometry_ierr)
#else
    call build_render_lattice(ubnmax, result, status, nterms, hdim, ntri, &
      nchart, mp, np, charts, mn, rc, zs, u, geometry_ierr)
#endif
    if (status /= 0_8) return
    ! Keep the generated C ABI portable: LFortran's current C backend does
    ! not transfer allocatable descriptors correctly for MOVE_ALLOC.
    allocate(result%sx(3,nsrc), result%snx(3,nsrc), result%sw(nsrc))
    allocate(result%ub(nsrc), result%ubn(nsrc), result%u(nsrc))
    result%sx = sx
    result%snx = snx
    result%sw = sw
    result%ub = ub
    result%ubn = ubn
    result%u = u
    deallocate(sx, snx, sw, ub, ubn, u)
#ifdef BIESOLVER_WASM_PROGRESS
    call biesolver_progress(PROGRESS_RESULT, int(result%nsrc,c_long_long), &
      int(result%nrender,c_long_long), 0_c_long_long, 0_c_long_long, &
      real(result%grf_error,c_double))
#endif
    if (present(t_fmm_out)) t_fmm_out = t_fmm
    if (present(t_close_out)) t_close_out = t_close
    if (present(timeinfo_out)) timeinfo_out = timeinfo

    deallocate(tgl_geo, wgl_geo, Dgl_geo, uvs, wts)
    deallocate(tpan, uvbd)
    deallocate(rts, rps)
#ifdef BIESOLVER_WASM_SCALAR_ONLY
    deallocate(cl, stat=ier)
#else
    deallocate(charts, stat=ier)
#endif
    deallocate(mn, stat=ier)
    deallocate(mtcs, offs, csr_idx, csr_val)

  end subroutine stellarator_run_case

#ifdef BIESOLVER_WASM_SCALAR_ONLY
  subroutine build_render_lattice(scale, use_w7x, result, status, &
      nterms, hdim, ntri, nchart, mp, np, tgl, cl, mn, rc, zs, u, &
      geometry_ierr)
    real(r64), intent(in) :: scale
    logical, intent(in) :: use_w7x
#else
  subroutine build_render_lattice(scale, result, status, nterms, hdim, &
      ntri, nchart, mp, np, charts, mn, rc, zs, u, geometry_ierr)
    real(r64), intent(in) :: scale
#endif
    type(stellarator_case_result), intent(inout) :: result
    integer(8), intent(inout) :: status
    integer(8), intent(in) :: nterms, hdim, ntri, nchart, mp, np
#ifdef BIESOLVER_WASM_SCALAR_ONLY
    real(r64), intent(in) :: tgl(nterms), u(:)
    real(r64), allocatable, intent(in) :: cl(:), mn(:), rc(:), zs(:)
#else
    real(r64), intent(in) :: charts(6,nchart), rc(:), zs(:), u(:)
    integer(8), intent(in) :: mn(:,:)
#endif
    integer(c_int), intent(out) :: geometry_ierr
    integer(8), parameter :: KRENDER = 4_8
    integer(8), parameter :: NVR = 15_8
    integer(8), parameter :: NFR = 16_8
    real(r64) :: render_uv(2,NVR), pols(NVR,hdim)
    real(r64) :: umatr(hdim,hdim), vmatr(hdim,hdim)
    real(r64) :: interp_matrix(NVR,hdim), render_u(NVR)
    integer(8) :: ii, jj, q, face, panel, voff, foff

    q = 0_8
    do ii = 0, KRENDER
      do jj = 0, KRENDER-ii
        q = q + 1_8
        render_uv(1,q) = real(ii,r64)/real(KRENDER,r64)
        render_uv(2,q) = real(jj,r64)/real(KRENDER,r64)
      end do
    end do
    call koorn_vals2coefs_coefs2vals(nterms-1_8, hdim, umatr, vmatr)
    call koorn_pols_batch_r64(NVR, render_uv, nterms-1_8, hdim, pols)
    interp_matrix = matmul(pols, umatr)

    result%nrender = ntri*NVR
    allocate(result%render_xyz(3,result%nrender))
    allocate(result%render_log_error(result%nrender))
    allocate(result%render_triangles(3,ntri*NFR))
    do panel = 1, ntri
      voff = (panel-1_8)*NVR
      foff = (panel-1_8)*NFR
#ifdef BIESOLVER_WASM_SCALAR_ONLY
      call map_panel_uv_bridge(use_w7x, cl, nchart, mp, np, nterms, tgl, &
        panel, NVR, render_uv, mn, rc, zs, &
        result%render_xyz(:,voff+1_8:voff+NVR), &
        geometry_ierr)
#else
      call map_panel_uv_bridge(charts, nchart, mp, np, nterms, panel, NVR, &
        render_uv, mn, rc, zs, result%render_xyz(:,voff+1_8:voff+NVR), &
        geometry_ierr)
#endif
#ifdef BIESOLVER_C_BACKEND_ROW_MAJOR
      if (geometry_ierr /= 0_c_int) then
        status = 13_8
        return
      end if
#endif
      render_u = matmul(interp_matrix, u((panel-1_8)*hdim+1_8:panel*hdim))
      do q = 1, NVR
        if (render_u(q) /= render_u(q) .or. abs(render_u(q)) > 1.0e300_r64) then
          status = 32_8
          return
        end if
        result%render_log_error(voff+q) = &
          log10(max(abs(render_u(q))/scale, 1.0e-16_r64))
      end do
      face = 0_8
      do ii = 0, KRENDER-1_8
        do jj = 0, KRENDER-1_8-ii
          face = face + 1_8
          result%render_triangles(:,foff+face) = [ &
            voff + render_vertex(ii,jj)-1_8, &
            voff + render_vertex(ii+1_8,jj)-1_8, &
            voff + render_vertex(ii,jj+1_8)-1_8 ]
          if (jj < KRENDER-1_8-ii) then
            face = face + 1_8
            result%render_triangles(:,foff+face) = [ &
              voff + render_vertex(ii+1_8,jj)-1_8, &
              voff + render_vertex(ii+1_8,jj+1_8)-1_8, &
              voff + render_vertex(ii,jj+1_8)-1_8 ]
          end if
        end do
      end do
    end do
  end subroutine build_render_lattice

  pure function render_vertex(ii, jj) result(index)
    integer(8), intent(in) :: ii, jj
    integer(8) :: index
    index = 1_8 + ii*5_8 - ii*(ii-1_8)/2_8 + jj
  end function render_vertex

  pure function tpar2uv(t, thold1, thold2, tholdinv1) result(uv)
    real(r64), intent(in) :: t, thold1, thold2, tholdinv1
    real(r64) :: uv(2)
    uv = 0.0_r64
    if (t < thold1) then
      uv(2) = 1.0_r64 - tholdinv1*t
    else if (t < thold2) then
      uv(1) = tholdinv1*(t - thold1)
    else
      uv(1) = 1.0_r64 - tholdinv1*(t - thold2)
      uv(2) = tholdinv1*(t - thold2)
    end if
  end function tpar2uv

  pure function f_harm(x) result(v)
    real(r64), intent(in) :: x(3)
    real(r64) :: v
    v = exp(STELLARATOR_LAM*x(1))*cos(STELLARATOR_LAM*x(2))
  end function f_harm

  pure function gradf_harm(x) result(g)
    real(r64), intent(in) :: x(3)
    real(r64) :: g(3)
    g(1) =  STELLARATOR_LAM*exp(STELLARATOR_LAM*x(1))* &
      cos(STELLARATOR_LAM*x(2))
    g(2) = -STELLARATOR_LAM*exp(STELLARATOR_LAM*x(1))* &
      sin(STELLARATOR_LAM*x(2))
    g(3) =  0.0_r64
  end function gradf_harm

#ifdef BIESOLVER_WASM_SCALAR_ONLY
  subroutine map_panel_uv_bridge(use_w7x, cl, nchart, mp, np, nterms, tgl, &
                                 itri, nuv, uv, mn, rc, zs, x, ierr)
    logical, intent(in) :: use_w7x
    real(r64), allocatable, intent(in) :: cl(:), mn(:), rc(:), zs(:)
#else
  subroutine map_panel_uv_bridge(charts, nchart, mp, np, nterms, itri, &
                                 nuv, uv, mn, rc, zs, x, ierr)
    real(r64), intent(in) :: charts(6,nchart), rc(:), zs(:)
    integer(8), intent(in) :: mn(:,:)
#endif
    integer(8), intent(in) :: nchart, mp, np, nterms, itri, nuv
#ifdef BIESOLVER_WASM_SCALAR_ONLY
    real(r64), intent(in) :: tgl(nterms), uv(2,nuv)
#else
    real(r64), intent(in) :: uv(2,nuv)
    integer(8) :: mesh_ier
#endif
    real(r64), intent(out) :: x(3,nuv)
    integer(c_int), intent(out) :: ierr
    ierr = 0_c_int
#ifdef BIESOLVER_WASM_SCALAR_ONLY
    if (use_w7x) then
#ifdef BIESOLVER_C_BACKEND_ROW_MAJOR
      ierr = biesolver_uv2x_charts_rowmajor(cl, int(nchart,c_long_long), &
             int(mp,c_long_long), int(np,c_long_long), int(nterms,c_long_long), &
             tgl, int(itri,c_long_long), int(nuv,c_long_long), uv, &
             int(W7X_NFP,c_long_long), int(W7X_NMODE,c_long_long), &
             mn, rc, zs, x)
#else
      call stellarator_uv2x_charts(cl, int(nchart,c_long_long), &
             int(mp,c_long_long), int(np,c_long_long), int(nterms,c_long_long), &
             tgl, int(itri,c_long_long), int(nuv,c_long_long), uv, &
             int(W7X_NFP,c_long_long), int(W7X_NMODE,c_long_long), &
             mn, rc, zs, x)
#endif
    else
#ifdef BIESOLVER_C_BACKEND_ROW_MAJOR
      ierr = biesolver_geo_uv2x_rowmajor(int(mp,c_long_long), &
             int(np,c_long_long), int(nterms,c_long_long), tgl, &
             int(itri,c_long_long), int(nuv,c_long_long), uv, x)
#else
      call stellarator_geo_uv2x(int(mp,c_long_long), int(np,c_long_long), &
             int(nterms,c_long_long), tgl, int(itri,c_long_long), &
             int(nuv,c_long_long), uv, x)
#endif
    end if
#else
    call stellarator_tri_uv2x_r64(mp, np, nterms, merge(W7X_NFP,1_8, &
         size(mn,2)>0), int(size(mn,2),8), mn, rc, zs, nchart, charts, &
         itri, nuv, uv, x, mesh_ier)
    ierr = int(mesh_ier,c_int)
#endif
  end subroutine map_panel_uv_bridge

  subroutine stellarator_result_clear(result)
    type(stellarator_case_result), intent(inout) :: result
    integer :: stat
    deallocate(result%sx, stat=stat)
    deallocate(result%snx, stat=stat)
    deallocate(result%sw, stat=stat)
    deallocate(result%ub, stat=stat)
    deallocate(result%ubn, stat=stat)
    deallocate(result%u, stat=stat)
    deallocate(result%render_xyz, stat=stat)
    deallocate(result%render_log_error, stat=stat)
    deallocate(result%render_triangles, stat=stat)
    result%ntri = 0_8
    result%nsrc = 0_8
    result%nrender = 0_8
    result%grf_error = 0.0_r64
  end subroutine stellarator_result_clear

  function stellarator_wall_time() result(t)
    real(r64) :: t
#ifndef BIESOLVER_WASM_SCALAR_ONLY
    t = omp_get_wtime()
#else
    t = 0.0_r64
#endif
  end function stellarator_wall_time

end module stellarator_grf_core_mod
