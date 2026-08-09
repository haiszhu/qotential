program stellarator_fmm_main
  use stellarator_grf_core_mod, only: stellarator_case_config, &
       stellarator_case_result, stellarator_run_case, stellarator_result_clear
  use lap3d_close_mod, only: simplex_precomp_r64
  use patch_refine_mod, only: PR_ADAPTIVE, PR_SQUARE
  implicit none

  type(stellarator_case_config) :: cfg
  type(stellarator_case_result) :: result
  integer(8) :: status, nfaces, nquad, hdim
  real(8), allocatable :: tgl(:), wgl(:), Dgl(:,:), w_bclag(:)
  real(8), allocatable :: Legmat(:,:), umatr(:,:), vmatr(:,:)
  integer :: unit, ios
  character(len=512) :: output
  character(len=64) :: argument
  character(len=8), parameter :: MAGIC = 'STGRF001'

  if (command_argument_count() /= 7) then
    error stop &
      'usage: stellarator_fmm <output.bin> <mp> <np> <order> <surface> <restol> <eps>'
  end if
  call get_command_argument(1, output)
  call get_command_argument(2, argument)
  read(argument, *, iostat=ios) cfg%mp
  if (ios /= 0 .or. cfg%mp <= 0_8) error stop 'mp must be a positive integer'
  call get_command_argument(3, argument)
  read(argument, *, iostat=ios) cfg%np
  if (ios /= 0 .or. cfg%np <= 0_8) error stop 'np must be a positive integer'
  call get_command_argument(4, argument)
  read(argument, *, iostat=ios) cfg%order
  if (ios /= 0 .or. cfg%order < 4_8 .or. cfg%order > 16_8 .or. &
      mod(cfg%order, 2_8) /= 0_8) then
    error stop 'order must be one of 4, 6, 8, 10, 12, 14, 16'
  end if
  call get_command_argument(5, argument)
  select case (trim(argument))
  case ('builtin')
    cfg%use_w7x = .false.
  case ('w7x')
    cfg%use_w7x = .true.
  case default
    error stop 'surface must be builtin or w7x'
  end select
  call get_command_argument(6, argument)
  read(argument, *, iostat=ios) cfg%restol
  if (ios /= 0 .or. cfg%restol <= 0.0_8) then
    error stop 'restol must be a positive real'
  end if
  call get_command_argument(7, argument)
  read(argument, *, iostat=ios) cfg%fmm_eps
  if (ios /= 0 .or. .not. (cfg%fmm_eps == 1.0e-3_8 .or. &
      cfg%fmm_eps == 1.0e-6_8 .or. cfg%fmm_eps == 1.0e-9_8 .or. &
      cfg%fmm_eps == 1.0e-12_8 .or. cfg%fmm_eps == 1.0e-15_8)) then
    error stop 'eps must be one of 1e-3, 1e-6, 1e-9, 1e-12, 1e-15'
  end if

  PR_ADAPTIVE = .false.
  PR_SQUARE = .true.
  cfg%isimd = 0_8
  cfg%ichart = 1_8
  cfg%use_fmm = .true.
  nquad = cfg%order + 2_8
  hdim = cfg%order*(cfg%order+1_8)/2_8
  allocate(tgl(nquad), wgl(nquad), Dgl(nquad,nquad), w_bclag(nquad))
  allocate(Legmat(nquad,nquad), umatr(hdim,hdim), vmatr(hdim,hdim))
  call simplex_precomp_r64(nquad, cfg%order-1_8, hdim, tgl, wgl, Dgl, &
                           w_bclag, Legmat, umatr, vmatr)
  call stellarator_run_case(cfg, tgl, wgl, Dgl, w_bclag, Legmat, &
      umatr, vmatr, result, status)
  if (status /= 0_8) then
    write(*,'(a,i0)') 'FMM solve status=', status
    error stop 'FMM solve failed'
  end if

  nfaces = size(result%render_triangles, 2, kind=8)
  open(newunit=unit, file=trim(output), access='stream', form='unformatted', &
       status='replace', action='write', convert='little_endian')
  write(unit) MAGIC, result%ntri, result%nsrc, result%nrender, nfaces, result%grf_error
  write(unit) result%sx, result%snx, result%sw, result%ub, result%ubn, result%u
  write(unit) result%render_xyz, result%render_log_error, result%render_triangles
  close(unit)

  write(*,'(a,i0,a,i0,a,i0,a,es24.16,a,es10.2)') 'FMM ntri=', result%ntri, &
       ' nsrc=', result%nsrc, ' nrender=', result%nrender, &
       ' GRF max rel err=', result%grf_error, ' eps=', cfg%fmm_eps
  call stellarator_result_clear(result)
  deallocate(tgl, wgl, Dgl, w_bclag, Legmat, umatr, vmatr)
end program stellarator_fmm_main
