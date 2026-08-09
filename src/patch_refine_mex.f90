module qol_patch_levels_mex_state_mod
  use iso_c_binding, only: c_funloc, c_funptr, c_int
  use patch_refine_mod, only: patch_levels_t, build_patch_levels_r64, &
      lap3dsdlpmat_levels_r64, clear_patch_levels_r64
  implicit none
  private

  integer(8), parameter :: max_handles = 1_8

  type :: patch_levels_slot_t
    type(patch_levels_t) :: levels
    logical :: live = .false.
    integer(8) :: generation = 1_8
    integer(8) :: npols = 0_8
  end type patch_levels_slot_t

  type(patch_levels_slot_t), save :: slots(max_handles)
  logical, save :: registry_is_locked = .false.
  logical, save :: registry_has_exit_hook = .false.

  public :: patch_levels_create, patch_levels_evaluate
  public :: patch_levels_free, patch_levels_free_all

  interface
    subroutine mex_lock() bind(C, name='mexLock')
    end subroutine mex_lock

    subroutine mex_unlock() bind(C, name='mexUnlock')
    end subroutine mex_unlock

    integer(c_int) function mex_at_exit(callback) bind(C, name='mexAtExit')
      import :: c_funptr, c_int
      type(c_funptr), value :: callback
    end function mex_at_exit
  end interface

contains

  subroutine patch_levels_create(korder, npols, sx, rts, rps, order, &
                                 dist, nlevel, handle, ier)
    integer(8), intent(in) :: korder, npols, order, nlevel
    real(8), intent(in) :: sx(3,npols), rts(3,npols), rps(3,npols), dist
    integer(8), intent(out) :: handle, ier
    integer(8) :: slot, lifecycle_ier

    handle = 0_8
    ier = 1_8
    do slot = 1, max_handles
      if (.not. slots(slot)%live) then
        call build_patch_levels_r64(korder, npols, sx, rts, rps, order, &
                                    dist, nlevel, slots(slot)%levels)
        slots(slot)%live = .true.
        slots(slot)%npols = npols
        call patch_registry_lock(lifecycle_ier)
        if (lifecycle_ier /= 0_8) then
          call clear_patch_levels_r64(slots(slot)%levels)
          slots(slot)%live = .false.
          slots(slot)%npols = 0_8
          ier = lifecycle_ier
          return
        end if
        handle = slots(slot)%generation*max_handles + slot
        ier = 0_8
        return
      end if
    end do
  end subroutine patch_levels_create

  subroutine patch_levels_evaluate(m, tx, npols, handle, isimd, As, Ad, &
                                   idxs, ms, ier)
    integer(8), intent(in) :: m, npols, handle, isimd
    real(8), intent(in) :: tx(3,m)
    real(8), intent(inout) :: As(m,npols), Ad(m,npols)
    integer(8), intent(inout) :: idxs(m), ms
    integer(8), intent(out) :: ier
    integer(8) :: slot
    logical :: valid

    call decode_handle(handle, slot, valid)
    if (.not. valid) then
      ms = 0_8
      ier = 2_8
      return
    end if
    if (npols /= slots(slot)%npols) then
      ms = 0_8
      ier = 3_8
      return
    end if

    call lap3dsdlpmat_levels_r64(m, tx, npols, slots(slot)%levels, &
                                 isimd, As, Ad, idxs, ms)
    ier = 0_8
  end subroutine patch_levels_evaluate

  subroutine patch_levels_free(handle, ier)
    integer(8), intent(in) :: handle
    integer(8), intent(out) :: ier
    integer(8) :: slot
    logical :: valid

    call decode_handle(handle, slot, valid)
    if (.not. valid) then
      ier = 2_8
      return
    end if

    call clear_patch_levels_r64(slots(slot)%levels)
    slots(slot)%live = .false.
    slots(slot)%npols = 0_8
    slots(slot)%generation = slots(slot)%generation + 1_8
    call patch_registry_unlock()
    ier = 0_8
  end subroutine patch_levels_free

  subroutine patch_levels_free_all(ier)
    integer(8), intent(out) :: ier
    integer(8) :: slot
    logical :: had_live

    had_live = .false.
    do slot = 1, max_handles
      if (slots(slot)%live) then
        had_live = .true.
        call clear_patch_levels_r64(slots(slot)%levels)
        slots(slot)%live = .false.
        slots(slot)%npols = 0_8
        slots(slot)%generation = slots(slot)%generation + 1_8
      end if
    end do
    if (had_live) call patch_registry_unlock()
    ier = 0_8
  end subroutine patch_levels_free_all

  subroutine patch_levels_free_all_at_exit() &
      bind(C, name='qol_patch_levels_free_all_at_exit')
    integer(8) :: slot

    do slot = 1, max_handles
      if (slots(slot)%live) then
        call clear_patch_levels_r64(slots(slot)%levels)
        slots(slot)%live = .false.
        slots(slot)%npols = 0_8
        slots(slot)%generation = slots(slot)%generation + 1_8
      end if
    end do
    registry_is_locked = .false.
    registry_has_exit_hook = .false.
  end subroutine patch_levels_free_all_at_exit

  subroutine patch_registry_lock(ier)
    integer(8), intent(out) :: ier
    integer(c_int) :: mex_ier

    ier = 0_8
    if (.not. registry_has_exit_hook) then
      mex_ier = mex_at_exit(c_funloc(patch_levels_free_all_at_exit))
      if (mex_ier /= 0_c_int) then
        ier = 4_8
        return
      end if
      registry_has_exit_hook = .true.
    end if
    if (.not. registry_is_locked) then
      call mex_lock()
      registry_is_locked = .true.
    end if
  end subroutine patch_registry_lock

  subroutine patch_registry_unlock()
    if (registry_is_locked) then
      registry_is_locked = .false.
      call mex_unlock()
    end if
  end subroutine patch_registry_unlock

  subroutine decode_handle(handle, slot, valid)
    integer(8), intent(in) :: handle
    integer(8), intent(out) :: slot
    logical, intent(out) :: valid
    integer(8) :: generation

    if (handle <= 0_8) then
      slot = 0_8
      valid = .false.
      return
    end if

    slot = mod(handle-1_8, max_handles) + 1_8
    generation = (handle-slot)/max_handles
    valid = slots(slot)%live .and. &
            slots(slot)%generation == generation
  end subroutine decode_handle

end module qol_patch_levels_mex_state_mod


subroutine qol_build_patch_levels_mex(korder, npols, sx, rts, rps, order, &
                                      dist, nlevel, lv, ier)
  use qol_patch_levels_mex_state_mod, only: patch_levels_create
  implicit none
  integer(8), intent(in) :: korder, npols, order, nlevel
  real(8), intent(in) :: sx(3,npols), rts(3,npols), rps(3,npols), dist
  integer(8), intent(out) :: lv, ier

  call patch_levels_create(korder, npols, sx, rts, rps, order, dist, &
                           nlevel, lv, ier)
end subroutine qol_build_patch_levels_mex


subroutine qol_lap3dsdlpmat_levels_mex(m, tx, npols, lv, isimd, &
                                      As, Ad, idxs, ms, ier)
  use qol_patch_levels_mex_state_mod, only: patch_levels_evaluate
  implicit none
  integer(8), intent(in) :: m, npols, lv, isimd
  real(8), intent(in) :: tx(3,m)
  real(8), intent(inout) :: As(m,npols), Ad(m,npols)
  integer(8), intent(inout) :: idxs(m), ms
  integer(8), intent(out) :: ier

  call patch_levels_evaluate(m, tx, npols, lv, isimd, As, Ad, idxs, &
                             ms, ier)
end subroutine qol_lap3dsdlpmat_levels_mex


subroutine qol_clear_patch_levels_mex(lv, ier)
  use qol_patch_levels_mex_state_mod, only: patch_levels_free
  implicit none
  integer(8), intent(in) :: lv
  integer(8), intent(out) :: ier

  call patch_levels_free(lv, ier)
end subroutine qol_clear_patch_levels_mex


subroutine qol_clear_patch_levels_all_mex(ier)
  use qol_patch_levels_mex_state_mod, only: patch_levels_free_all
  implicit none
  integer(8), intent(out) :: ier

  call patch_levels_free_all(ier)
end subroutine qol_clear_patch_levels_all_mex
