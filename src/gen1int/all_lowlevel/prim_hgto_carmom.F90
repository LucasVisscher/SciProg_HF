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
!!  This file contains the recurrence relations of Cartesian multipole moment
!!  integrals using primitive Hermite Gaussians.
!!
!!  2012-02-22, Bin Gao:
!!  * rewrites using new recurrence relations
!!
!!  2010-10-10, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief recurrence relations of Cartesian multipole moment integrals using
  !>        primitive Hermite Gaussians
  !> \author Bin Gao
  !> \date 2010-10-10
  !> \param orders_hgto_bra is the range of orders of Hermite Gaussians on bra center
  !> \param coord_bra is the coordinates of bra center
  !> \param exponent_bra is the exponent of primitive Gaussian of bra center
  !> \param orders_hgto_ket is the range of orders of Hermite Gaussians on ket center
  !> \param coord_ket is the coordinates of ket center
  !> \param exponent_ket is the exponent of primitive Gaussian of ket center
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scale constant for Cartesian multipole moments
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_elec is the order of electronic derivatives
  !> \param dim_hgto_bra is the dimension of Hermite Gaussians on bra center
  !> \param dim_hgto_ket is the dimension of Hermite Gaussians on ket center
  !> \param num_elec is the number of xyz components of electronic derivatives,
  !>        equals to \f$(\var(order_elec)+1)(\var(order_elec)+2)/2\f$
  !> \param num_mom is the number of xyz components of Cartesian multipole moment,
  !>        equals to \f$(\var(order_mom)+1)(\var(order_mom)+2)/2\f$
  !> \return hgto_pints contains the primitive HGTO integrals
  subroutine prim_hgto_carmom(orders_hgto_bra, coord_bra, exponent_bra, &
                              orders_hgto_ket, coord_ket, exponent_ket, &
                              order_elec, dipole_origin, scal_const,    &
                              order_mom, dim_hgto_bra, dim_hgto_ket,    &
                              num_elec, num_mom, hgto_pints)
    use xkind
    implicit none
    integer, intent(in) :: orders_hgto_bra(2)
    real(REALK), intent(in) :: coord_bra(3)
    real(REALK), intent(in) :: exponent_bra
    integer, intent(in) :: orders_hgto_ket(2)
    real(REALK), intent(in) :: coord_ket(3)
    real(REALK), intent(in) :: exponent_ket
    integer, intent(in) :: order_elec
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: order_mom
    integer, intent(in) :: dim_hgto_bra
    integer, intent(in) :: dim_hgto_ket
    integer, intent(in) :: num_elec
    integer, intent(in) :: num_mom
    real(REALK), intent(out) :: hgto_pints(dim_hgto_bra,dim_hgto_ket,num_elec,num_mom)
!f2py intent(in) :: orders_hgto_bra
!f2py intent(in) :: coord_bra
!f2py intent(in) :: exponent_bra
!f2py intent(in) :: orders_hgto_ket
!f2py intent(in) :: coord_ket
!f2py intent(in) :: exponent_ket
!f2py intent(in) :: order_elec
!f2py intent(in) :: dipole_origin
!f2py intent(in) :: scal_const
!f2py intent(in) :: order_mom
!f2py intent(in) :: dim_hgto_bra
!f2py intent(in) :: dim_hgto_ket
!f2py intent(in) :: num_elec
!f2py intent(in) :: num_mom
!f2py intent(out) :: hgto_pints
!f2py depend(dim_hgto_bra) :: hgto_pints
!f2py depend(dim_hgto_ket) :: hgto_pints
!f2py depend(num_elec) :: hgto_pints
!f2py depend(num_mom) :: hgto_pints
    integer orders_hket_elec(2)                    !orders of HGTOs on ket center with electronic derivatives
    integer orders_hbra(2)                         !orders of HGTOs on bra center for \fn(carmom_hbra)
    integer dim_hbra                               !dimension of integrals from \fn(carmom_hbra)
    real(REALK), allocatable :: hbra_pints(:)      !integrals from \fn(carmom_hbra)
    integer orders_hbra_zero(2)                    !orders of HGTOs on bra center for \fn(carmom_hrr_ket)
    integer orders_hket_zero(2)                    !orders of HGTOs on ket center for \fn(carmom_hrr_ket)
    integer dim_hbra_zero                          !dimensions of integrals from \fn(carmom_hrr_ket)
    integer dim_hket_zero
    real(REALK), allocatable :: hket_pints(:,:)    !integrals from \fn(carmom_hrr_ket)
    integer dim_hket_mom                           !dimension of integrals from \fn(carmom_moment)
    real(REALK), allocatable :: hmom_pints(:,:,:)  !integrals from \fn(carmom_moment)
    integer ierr                                   !error information
#if defined(XTIME)
    real(REALK) curr_time                          !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! sets the orders of HGTOs on ket center with electronic derivatives
    orders_hket_elec(1) = orders_hgto_ket(1)+order_elec
    orders_hket_elec(2) = orders_hgto_ket(2)+order_elec
    ! sets the orders of HGTOs for \fn(carmom_hrr_ket)
    orders_hbra_zero(1) = max(0,orders_hgto_bra(1)-order_mom)
    orders_hbra_zero(2) = orders_hgto_bra(2)
    orders_hket_zero(1) = max(0,orders_hket_elec(1)-order_mom)
    orders_hket_zero(2) = orders_hket_elec(2)
    ! sets the orders of HGTOs on bra center for \fn(carmom_hbra)
    orders_hbra = orders_hbra_zero+orders_hket_zero
    ! allocates memory for the integrals from \fn(carmom_hbra)
    if (orders_hbra(1)>0) then
      dim_hbra = ((orders_hbra(2)+1)*(orders_hbra(2)+2)*(orders_hbra(2)+3) &
               -  orders_hbra(1)*(orders_hbra(1)+1)*(orders_hbra(1)+2))/6
    else
      dim_hbra = (orders_hbra(2)+1)*(orders_hbra(2)+2)*(orders_hbra(2)+3)/6
    end if
    allocate(hbra_pints(dim_hbra), stat=ierr)
    if (ierr/=0) &
      call error_stop("prim_hgto_carmom", "failed to allocate hbra_pints", dim_hbra)
    ! calculates the integrals using \fn(carmom_hbra)
    call carmom_hbra(orders_hbra, coord_bra, exponent_bra, &
                     coord_ket, exponent_ket, order_elec,  &
                     scal_const, dim_hbra, hbra_pints)
    ! allocates memory for the integrals from \fn(carmom_hrr_ket)
    if (orders_hbra_zero(1)>0) then
      dim_hbra_zero = ((orders_hbra_zero(2)+1)*(orders_hbra_zero(2)+2) &
                       *(orders_hbra_zero(2)+3)                        &
                    -  orders_hbra_zero(1)*(orders_hbra_zero(1)+1)     &
                       *(orders_hbra_zero(1)+2))/6
    else
      dim_hbra_zero = (orders_hbra_zero(2)+1)*(orders_hbra_zero(2)+2) &
                    * (orders_hbra_zero(2)+3)/6
    end if
    if (orders_hket_zero(1)>0) then
      dim_hket_zero = ((orders_hket_zero(2)+1)*(orders_hket_zero(2)+2) &
                       *(orders_hket_zero(2)+3)                        &
                    -  orders_hket_zero(1)*(orders_hket_zero(1)+1)     &
                       *(orders_hket_zero(1)+2))/6
    else
      dim_hket_zero = (orders_hket_zero(2)+1)*(orders_hket_zero(2)+2) &
                    * (orders_hket_zero(2)+3)/6
    end if
    allocate(hket_pints(dim_hbra_zero,dim_hket_zero), stat=ierr)
    if (ierr/=0)                          &
      call error_stop("prim_hgto_carmom", &
                      "failed to allocate hket_pints", dim_hbra_zero*dim_hket_zero)
    ! calculates the integrals using \fn(carmom_hrr_ket)
    call carmom_hrr_ket(orders_hbra_zero, exponent_bra, &
                        orders_hket_zero, exponent_ket, &
                        dim_hbra, hbra_pints,           &
                        dim_hbra_zero, dim_hket_zero, hket_pints)
    deallocate(hbra_pints)
    if (order_elec==0) then
      ! calculates the integrals using \fn(carmom_moment)
      call carmom_moment(orders_hgto_bra, coord_bra, exponent_bra, &
                         orders_hgto_ket, coord_ket, exponent_ket, &
                         dipole_origin, order_mom, dim_hbra_zero,  &
                         dim_hket_zero, hket_pints, dim_hgto_bra,  &
                         dim_hgto_ket, hgto_pints)
      deallocate(hket_pints)
    else
      ! allocates memory for the integrals from \fn(carmom_moment)
      if (orders_hket_elec(1)>0) then
        dim_hket_mom = ((orders_hket_elec(2)+1)*(orders_hket_elec(2)+2) &
                        *(orders_hket_elec(2)+3)                        &
                     -  orders_hket_elec(1)*(orders_hket_elec(1)+1)     &
                        *(orders_hket_elec(1)+2))/6
      else
        dim_hket_mom = (orders_hket_elec(2)+1)*(orders_hket_elec(2)+2) &
                     * (orders_hket_elec(2)+3)/6
      end if
      allocate(hmom_pints(dim_hgto_bra,dim_hket_mom,num_mom), stat=ierr)
      ! calculates the integrals using \fn(carmom_moment)
      call carmom_moment(orders_hgto_bra, coord_bra, exponent_bra,  &
                         orders_hket_elec, coord_ket, exponent_ket, &
                         dipole_origin, order_mom, dim_hbra_zero,   &
                         dim_hket_zero, hket_pints, dim_hgto_bra,   &
                         dim_hket_mom, hmom_pints)
      deallocate(hket_pints)
      ! recovers the electronic derivatives
      call scatter_multi_inner(dim_hgto_bra, dim_hket_mom, 1, num_mom,  &
                               hmom_pints, orders_hgto_ket, order_elec, &
                               dim_hgto_ket, num_elec, hgto_pints)
      deallocate(hmom_pints)
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "prim_hgto_carmom", STDOUT)
#endif
    return
  end subroutine prim_hgto_carmom
