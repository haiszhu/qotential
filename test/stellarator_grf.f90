! Thin command-line driver for the shared stellarator GRF orchestration.
program stellarator_grf
  use stellarator_grf_core_mod, only: stellarator_case_config, &
       stellarator_case_result, stellarator_run_case, stellarator_result_clear
  use patch_refine_mod, only: PR_ADAPTIVE, PR_SQUARE
  implicit none

  integer(8), parameter :: NC = 7_8
  integer(8), parameter :: ORDER_L(NC) = [ 4_8,  6_8,  8_8, 10_8, 12_8, 14_8, 16_8]
  integer(8), parameter :: MP_L(NC)    = [24_8, 24_8, 36_8, 36_8, 48_8, 60_8, 72_8]
  integer(8), parameter :: NP_L(NC)    = [72_8, 72_8,108_8,108_8,144_8,180_8,216_8]

  type(stellarator_case_config) :: cfg
  type(stellarator_case_result) :: result
  integer(8) :: ncases, isimd, ichart, iadap, icase, status
  real(8) :: t_fmm, t_close, timeinfo(20)
  character(len=32) :: arg

  ncases = 3_8
  isimd = 0_8
  ichart = 1_8
  iadap = 0_8
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg)
    read(arg,*) ncases
  end if
  if (command_argument_count() >= 2) then
    call get_command_argument(2, arg)
    read(arg,*) isimd
  end if
  if (command_argument_count() >= 3) then
    call get_command_argument(3, arg)
    read(arg,*) ichart
  end if
  if (command_argument_count() >= 4) then
    call get_command_argument(4, arg)
    read(arg,*) iadap
  end if

  PR_ADAPTIVE = (iadap == 1_8)
  if (iadap == 1_8 .or. iadap == 3_8) PR_SQUARE = .false.

  do icase = 1, min(ncases, NC)
    cfg%order = ORDER_L(icase)
    cfg%mp = MP_L(icase)
    cfg%np = NP_L(icase)
    cfg%isimd = isimd
    cfg%ichart = ichart
    cfg%use_fmm = .true.
    cfg%use_w7x = .false.
    call stellarator_run_case(cfg, result, status, t_fmm, t_close, timeinfo)
    if (status /= 0_8) then
      write(*,'(a,i0,a,i0)') 'ERROR: case ', icase, ' failed with status ', status
      error stop 'stellarator solve failed'
    end if
    write(*,'(a,i0,a,i2,a,i3,a,i3,a,i0,a,es10.3,a,f8.2,a,f8.2,a)') &
        'case ', icase, ': order=', cfg%order, '  mp=', cfg%mp, '  np=', cfg%np, &
        '  N=', result%nsrc, '  GRF max rel err = ', result%grf_error, &
        '  (fmm ', t_fmm, ' s, close ', t_close, ' s)'
    write(*,'(a,8es10.2)') '  timeinfo ', timeinfo(1:8)
    call stellarator_result_clear(result)
  end do
end program stellarator_grf
