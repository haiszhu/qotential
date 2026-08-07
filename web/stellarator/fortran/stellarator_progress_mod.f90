! Deterministic 2-percent milestone state for the browser progress log.
!
! The decision is a pure function of integer counters, so the emitted sequence
! depends only on the problem size, never on wall-clock timing.  A loop of
! `total` iterations therefore produces at most 50 ordinary milestones plus one
! final milestone, and the final one is emitted exactly once even if the caller
! reaches `current == total` more than once.
!
! Deliberately dependency-free: it must lower cleanly through the LFortran C
! backend ahead of the solver core.

module stellarator_progress_mod
  implicit none
  private

  type, public :: stellarator_progress_state
    integer(8) :: next_pct = 2_8
    logical :: final_sent = .false.
  end type stellarator_progress_state

  public :: stellarator_progress_step
  public :: stellarator_progress_reset

contains

  ! Explicitly (re)initialize a milestone state before a loop.
  !
  ! Derived-type default component initialization is NOT relied upon here: the
  ! LFortran C backend that lowers the solver core does not emit a C initializer
  ! for a local `type(stellarator_progress_state)`, leaving next_pct/final_sent
  ! indeterminate on the stack.  Every caller must reset a state through this
  ! routine before its loop so the emitted milestone sequence is deterministic
  ! and the single final 100% event is always produced.
  subroutine stellarator_progress_reset(state)
    type(stellarator_progress_state), intent(inout) :: state
    state%next_pct = 2_8
    state%final_sent = .false.
  end subroutine stellarator_progress_reset

  ! Report whether `current` of `total` crosses the next 2-percent milestone.
  ! Real solver totals have already passed the ABI allocation guards and are
  ! far below huge(0_8)/100, so the percentage product cannot overflow.
  subroutine stellarator_progress_step(state, current, total, emit)
    type(stellarator_progress_state), intent(inout) :: state
    integer(8), intent(in) :: current, total
    logical, intent(out) :: emit
    integer(8) :: pct

    emit = .false.
    if (total <= 0_8 .or. current <= 0_8) return
    if (current >= total) then
      if (.not. state%final_sent) then
        emit = .true.
        state%final_sent = .true.
      end if
      return
    end if
    pct = (100_8*current)/total
    if (pct < state%next_pct) return
    emit = .true.
    do while (state%next_pct <= pct)
      state%next_pct = state%next_pct + 2_8
    end do
  end subroutine stellarator_progress_step

end module stellarator_progress_mod
