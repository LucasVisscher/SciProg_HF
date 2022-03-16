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
!!  This file defines the information of London atomic orbitals.
!!
!!  2013-06-11, Bin Gao:
!!  * first version

#include "stdout.h"

module london_ao
  use xkind
  implicit none

  ! types of Gaussian type orbitals (GTOs)
  integer, parameter, public :: NON_LAO = 0
  integer, parameter, public :: LONDON = 1
  integer, parameter, public :: ROT_LAO = 2

  ! London atomic orbital (LAO)
  type, public :: london_ao_t
    private
    ! type of GTOs
    integer :: gto_type = NON_LAO
    ! gauge origin of the magnetic vector potential
    real(REALK) :: gauge_origin(3) = (/0.0,0.0,0.0/)
    ! origin of the London phase factor
    real(REALK) :: origin_London_PF(3) = (/0.0,0.0,0.0/)
  end type london_ao_t

  public :: LondonAOSet
  public :: LondonAOUsed
  public :: LondonAOGetGaugeOrigin
  public :: LondonAOGetLPFOrigin
  public :: LondonAOView

  contains

  !> \brief sets the information of London atomic orbital
  !> \author Bin Gao
  !> \date 2013-01-31
  !> \param info_LAO contains the information of LAO
  !> \param gauge_origin contains the coordinates of gauge origin of the magnetic vector potential
  !> \param origin_London_PF contains the coordinates of origin of the London phase factor
  subroutine LondonAOSet(info_LAO, gauge_origin, origin_London_PF)
    type(london_ao_t), intent(inout) :: info_LAO
    real(REALK), optional, intent(in) :: gauge_origin(3)
    real(REALK), optional, intent(in) :: origin_London_PF(3)
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
!FIXME: implement rotational LAO
    if (present(gauge_origin)) then
      info_LAO%gauge_origin = gauge_origin
    else
      info_LAO%gauge_origin = 0.0
    end if
    if (present(origin_London_PF)) then
      info_LAO%origin_London_PF = origin_London_PF
      info_LAO%gto_type = LONDON
    else
      info_LAO%origin_London_PF = 0.0
      info_LAO%gto_type = NON_LAO
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "LondonAOSet", STDOUT)
#endif
  end subroutine LondonAOSet

  !> \brief returns if London atomic orbital used or not
  !> \author Bin Gao
  !> \date 2013-06-11
  !> \param info_LAO contains the information of LAO
  function LondonAOUsed(info_LAO)
    logical LondonAOUsed
    type(london_ao_t), intent(in) :: info_LAO
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    LondonAOUsed = (info_LAO%gto_type/=NON_LAO)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "LondonAOUsed", STDOUT)
#endif
  end function LondonAOUsed

  !> \brief gets the gauge origin
  !> \author Bin Gao
  !> \date 2013-06-11
  !> \param info_LAO contains the information of LAO
  !> \return gauge_origin contains the coordinates of gauge origin of the magnetic vector potential
  subroutine LondonAOGetGaugeOrigin(info_LAO, gauge_origin)
    type(london_ao_t), intent(in) :: info_LAO
    real(REALK), intent(out) :: gauge_origin(3)
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    gauge_origin = info_LAO%gauge_origin
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "LondonAOGetGaugeOrigin", STDOUT)
#endif
  end subroutine LondonAOGetGaugeOrigin

  !> \brief gets the origin of the London phase factor
  !> \author Bin Gao
  !> \date 2013-06-11
  !> \param info_LAO contains the information of LAO
  !> \return origin_London_PF contains the coordinates of origin of the London phase factor
  subroutine LondonAOGetLPFOrigin(info_LAO, origin_London_PF)
    type(london_ao_t), intent(in) :: info_LAO
    real(REALK), intent(out) :: origin_London_PF(3)
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    origin_London_PF = info_LAO%origin_London_PF
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "LondonAOGetLPFOrigin", STDOUT)
#endif
  end subroutine LondonAOGetLPFOrigin

  !> \brief visualizes the information of LAO
  !> \author Bin Gao
  !> \date 2013-06-11
  !> \param info_LAO contains the information of LAO
  !> \param io_viewer is the logical unit number of the viewer
  subroutine LondonAOView(info_LAO, io_viewer)
    type(london_ao_t), intent(in) :: info_LAO
    integer, intent(in) :: io_viewer
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    select case(info_LAO%gto_type)
    case(NON_LAO)
      write(io_viewer,100) "non-LAO"
    case(LONDON)
      write(io_viewer,100) "London AO"
      write(io_viewer,100) "gauge origin of the magnetic vector potential", info_LAO%gauge_origin
      write(io_viewer,100) "origin of the London phase factor", info_LAO%origin_London_PF
    case(ROT_LAO)
      write(io_viewer,100) "rotational LAO"
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "LondonAOView", STDOUT)
#endif
100 format("LondonAOView>> ",A,3F16.8)
  end subroutine LondonAOView

end module london_ao
