! Milestone sequences must be deterministic, monotone, bounded, and end
! exactly once at 100 percent.
program progress_milestone_test
  use stellarator_progress_mod, only: stellarator_progress_state, &
                                      stellarator_progress_step, &
                                      stellarator_progress_reset
  implicit none

  integer(8), allocatable :: emitted(:)
  integer(8) :: n

  ! total smaller than the milestone spacing: every step emits
  call collect(3_8, emitted, n)
  call require(n == 3_8, 'total=3 must emit three milestones')
  call require(emitted(1) == 1_8 .and. emitted(2) == 2_8 .and. &
               emitted(3) == 3_8, 'total=3 must emit 1,2,3')

  ! total exactly divisible into 2-percent steps
  call collect(50_8, emitted, n)
  call require(n == 50_8, 'total=50 must emit fifty milestones')
  call require(emitted(1) == 1_8 .and. emitted(50) == 50_8, &
               'total=50 must span 1..50')

  ! non-divisible total: monotone, unique, bounded, final present once
  call collect(73_8, emitted, n)
  call check_sequence(emitted, n, 73_8)

  call require_final_once()

  call require_reset_reuse()

  print '(A)', 'PROGRESS_MILESTONE_OK'

contains

  subroutine collect(total, values, count)
    integer(8), intent(in) :: total
    integer(8), allocatable, intent(out) :: values(:)
    integer(8), intent(out) :: count
    type(stellarator_progress_state) :: state
    logical :: emit
    integer(8) :: i

    ! Match production: reset explicitly rather than relying on derived-type
    ! default initialization, which the LFortran C backend does not lower.
    call stellarator_progress_reset(state)
    allocate(values(total))
    count = 0_8
    do i = 1_8, total
      call stellarator_progress_step(state, i, total, emit)
      if (emit) then
        count = count + 1_8
        values(count) = i
      end if
    end do
  end subroutine collect

  subroutine check_sequence(values, count, total)
    integer(8), intent(in) :: values(:), count, total
    integer(8) :: i, ordinary, finals

    call require(count >= 1_8, 'sequence must not be empty')
    do i = 2_8, count
      call require(values(i) > values(i-1_8), &
                   'milestones must strictly increase')
    end do
    finals = 0_8
    ordinary = 0_8
    do i = 1_8, count
      if (values(i) == total) then
        finals = finals + 1_8
      else
        ordinary = ordinary + 1_8
      end if
    end do
    call require(finals == 1_8, 'final milestone must appear exactly once')
    call require(ordinary <= 50_8, 'at most fifty ordinary milestones')
  end subroutine check_sequence

  ! Reaching the end twice must not produce a second final milestone.
  subroutine require_final_once()
    type(stellarator_progress_state) :: state
    logical :: first, second

    call stellarator_progress_reset(state)
    call stellarator_progress_step(state, 73_8, 73_8, first)
    call stellarator_progress_step(state, 73_8, 73_8, second)
    call require(first, 'first call at current=total must emit')
    call require(.not. second, 'second call at current=total must not emit')
  end subroutine require_final_once

  ! A state reused across loops must be reset: reset restores next_pct to the
  ! first milestone and clears final_sent, so the final 100% event is emitted
  ! again.  This locks the fix for the LFortran C-backend default-init gap; a
  ! regression that drops the reset would leave final_sent latched.
  subroutine require_reset_reuse()
    type(stellarator_progress_state) :: state
    logical :: a, b, c, d

    call stellarator_progress_reset(state)
    call stellarator_progress_step(state, 1_8, 10_8, a)    ! first 2% milestone
    call stellarator_progress_step(state, 10_8, 10_8, b)   ! final, emits
    call stellarator_progress_step(state, 10_8, 10_8, c)   ! no repeat final
    call stellarator_progress_reset(state)
    call stellarator_progress_step(state, 10_8, 10_8, d)   ! after reset, emits
    call require(a, 'reset: first milestone must emit from current=1')
    call require(b, 'reset: final must emit at current=total')
    call require(.not. c, 'reset: repeated final must not emit before reset')
    call require(d, 'reset: after reset the final must emit again')
  end subroutine require_reset_reuse

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      print '(A,A)', 'PROGRESS_MILESTONE_FAIL: ', message
      stop 1
    end if
  end subroutine require

end program progress_milestone_test
