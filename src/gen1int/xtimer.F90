!!  gen1int: compute one-electron integrals using rotational London atomic-orbitals
!!  Copyright 2009-2012 Bin Gao, and Andreas Thorvaldsen
!! 
!!  gen1int is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!! 
!!  gen1int is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU Lesser General Public License for more details.
!! 
!!  You should have received a copy of the GNU Lesser General Public License
!!  along with gen1int. If not, see <http://www.gnu.org/licenses/>.
!!
!!  This file returns the time used by calling \fn(cpu_time).
!!
!!  2009-10-23, Bin Gao:
!!  * first version

  !> \brief sets current CPU time
  !> \author Bin Gao
  !> \date 2009-10-24
  !> \param curr_time is current time from \fn(cpu_time)
  subroutine xtimer_set(curr_time)
    use xkind
    implicit none
    real(REALK), intent(out) :: curr_time
!f2py intent(out) :: curr_time
    ! gets the CPU elapsed time
    call cpu_time(curr_time)
    return
  end subroutine xtimer_set

  !> \brief prints the CPU elapsed time
  !> \author Bin Gao
  !> \date 2009-10-24
  !> \param prev_time is the previous time
  !> \param msg_timer contains message to print
  !> \param io_viewer is the logical unit number of the viewer
  subroutine xtimer_view(prev_time, msg_timer, io_viewer)
    use xkind
    implicit none
    real(REALK), intent(in) :: prev_time
    character*(*), intent(in) :: msg_timer
    integer, intent(in) :: io_viewer
!f2py intent(inout) :: prev_time
!f2py intent(in) :: msg_timer
!f2py intent(in) :: io_viewer
    real(REALK) curr_time  !current CPU elapsed time
    ! gets the CPU elapsed time
    call cpu_time(curr_time)
    write(io_viewer,100) msg_timer, curr_time-prev_time
    !-write(io_viewer,'()')
    !-! sets current CPU time
    !-call cpu_time(prev_time)
    return
100 format('$$ CPU TIME of ',A,1X,F12.6,' seconds')
  end subroutine xtimer_view
