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
!!  This file handles error stop.
!!
!!  2011-12-13, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief handles error stop
  !> \author Bin Gao
  !> \date 2011-12-13
  !> \param prev_time is the previous time
  !> \param msg_timer contains message to print
  subroutine error_stop(sub_name, error_msg, error_code)
    implicit none
    character*(*), intent(in) :: sub_name
    character*(*), intent(in) :: error_msg
    integer, intent(in) :: error_code
!f2py intent(in) :: sub_name
!f2py intent(in) :: error_msg
!f2py intent(in) :: error_code
    write(STDOUT,100) sub_name, error_msg, error_code
    stop
100 format("ERROR occurred in sub. ",A," with message '",A,"' and '",I12,"'")
  end subroutine error_stop
