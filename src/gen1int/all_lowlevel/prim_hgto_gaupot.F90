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
!!  This file contains the recurrence relations of Gaussian charge potential
!!  (multiplied by Cartesian multipole moment) integrals using primitive Hermite
!!  Gaussians.
!!
!!  2011-12-08, Bin Gao
!!  * first version

#include "stdout.h"

  !> \brief recurrence relations of Gaussian charge potential (multiplied by
  !>        Cartesian multipole moment) integrals using primitive Hermite Gaussians
  !> \author Bin Gao
  !> \date 2011-12-08
  !> \param orders_hgto_bra is the range of orders of Hermite Gaussians on bra center
  !> \param coord_bra contains the coordinates of bra center
  !> \param exponent_bra is the exponent of primitive Gaussian of bra center
  !> \param orders_hgto_ket is the range of orders of Hermite Gaussians on ket center
  !> \param coord_ket contains the coordinates of ket center
  !> \param exponent_ket is the exponent of primitive Gaussian of ket center
  !> \param gaupot_origin contains the coordinates of Gaussian charge potential origin
  !> \param gaupot_expt is the exponent used in the Gaussian broadening function of the charge
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scale constant for Gaussian charge potential
  !> \param order_geo_pot is the order of geometric derivatives on potential origin
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_elec is the order of electronic derivatives
  !> \param dim_hgto_bra is the dimension of Hermite Gaussians on bra center
  !> \param dim_hgto_ket is the dimension of Hermite Gaussians on ket center
  !> \param num_elec is the number of xyz components of electronic derivatives,
  !>        equals to \f$(\var(order_elec)+1)(\var(order_elec)+2)/2\f$
  !> \param num_mom is the number of xyz components of Cartesian multipole moment,
  !>        equals to \f$(\var(order_mom)+1)(\var(order_mom)+2)/2\f$
  !> \param num_geo_pot is the number of geometric derivatives on potential origin,
  !>        equals to \f$(\var(order_geo_pot)+1)(\var(order_geo_pot)+2)/2\f$
  !> \return hgto_pints contains the primitive HGTO integrals
  subroutine prim_hgto_gaupot(orders_hgto_bra, coord_bra, exponent_bra,  &
                              orders_hgto_ket, coord_ket, exponent_ket,  &
                              order_elec, gaupot_origin, gaupot_expt,    &
                              dipole_origin, scal_const, order_mom,      &
                              order_geo_pot, dim_hgto_bra, dim_hgto_ket, &
                              num_elec, num_mom, num_geo_pot, hgto_pints)
    use xkind
    implicit none
    integer, intent(in) :: orders_hgto_bra(2)
    real(REALK), intent(in) :: coord_bra(3)
    real(REALK), intent(in) :: exponent_bra
    integer, intent(in) :: orders_hgto_ket(2)
    real(REALK), intent(in) :: coord_ket(3)
    real(REALK), intent(in) :: exponent_ket
    integer, intent(in) :: order_elec
    real(REALK), intent(in) :: gaupot_origin(3)
    real(REALK), intent(in) :: gaupot_expt
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_geo_pot
    integer, intent(in) :: dim_hgto_bra
    integer, intent(in) :: dim_hgto_ket
    integer, intent(in) :: num_elec
    integer, intent(in) :: num_mom
    integer, intent(in) :: num_geo_pot
    real(REALK), intent(out) :: hgto_pints(dim_hgto_bra,dim_hgto_ket, &
                                           num_elec,num_mom,num_geo_pot)
!f2py intent(in) :: orders_hgto_bra
!f2py intent(in) :: coord_bra
!f2py intent(in) :: exponent_bra
!f2py intent(in) :: orders_hgto_ket
!f2py intent(in) :: coord_ket
!f2py intent(in) :: exponent_ket
!f2py intent(in) :: order_elec
!f2py intent(in) :: gaupot_origin
!f2py intent(in) :: gaupot_expt
!f2py intent(in) :: dipole_origin
!f2py intent(in) :: scal_const
!f2py intent(in) :: order_mom
!f2py intent(in) :: order_geo_pot
!f2py intent(in) :: dim_hgto_bra
!f2py intent(in) :: dim_hgto_ket
!f2py intent(in) :: num_elec
!f2py intent(in) :: num_mom
!f2py intent(in) :: num_geo_pot
!f2py intent(out) :: hgto_pints
!f2py depend(dim_hgto_bra) :: hgto_pints
!f2py depend(dim_hgto_ket) :: hgto_pints
!f2py depend(num_elec) :: hgto_pints
!f2py depend(num_mom) :: hgto_pints
!f2py depend(num_geo_pot) :: hgto_pints
    integer orders_hket_elec(2)  !orders of HGTOs on ket center with electronic derivatives
    integer orders_hbra_mom(2)   !orders of HGTOs on bra center with Cartesian multipole moments
    integer max_order_gaupot     !maximum order of geometric derivatives for \fn(gaupot_geom)
    integer dim_gaupot           !dimension of geometric derivatives of potential origin
                                 !with zeroth order HGTOs
    real(REALK), allocatable :: gaupot_pints(:)
                                 !geometric derivatives of potential origin with zeroth order HGTOs
    integer min_hgto_ket         !minimum order of HGTOs on ket center after recovering HGTOs on ket center
    integer max_geo_hket         !maximum order of geometric derivatives after recovering HGTOs on ket center
    integer dim_hket             !dimension of HGTOs on ket center after recovering HGTOs on ket center
    integer dim_geo_hket         !dimension of geometric derivatives after recovering HGTOs on ket center
    real(REALK), allocatable :: hket_pints(:,:)
                                 !integrals after recovering HGTOs on ket center
    integer dim_hbra_zero        !dimension of HGTOs on bra center with zeroth order Cartesian multipole moment
    integer dim_hket_zero        !dimension of HGTOs on ket center with zeroth order Cartesian multipole moment
    real(REALK), allocatable :: hbra_pints(:,:,:)
                                 !integrals after recovering HGTOs on bra center
    real(REALK), allocatable :: hmom_pints(:,:,:,:)
                                 !integrals after recovering Cartesian multipole moments
    integer ierr                 !error information
#if defined(XTIME)
    real(REALK) curr_time        !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! sets the orders of HGTOs on ket center with electronic derivatives
    orders_hket_elec(1) = orders_hgto_ket(1)+order_elec
    orders_hket_elec(2) = orders_hgto_ket(2)+order_elec
    ! sets the orders of HGTOs on bra center with Cartesian multipole moments
    orders_hbra_mom(1) = max(0,orders_hgto_bra(1)-order_mom)
    orders_hbra_mom(2) = orders_hgto_bra(2)+order_mom
    ! sets the maximum order and dimension of geometric derivatives for \fn(gaupot_geom)
    max_order_gaupot = orders_hbra_mom(2)+orders_hket_elec(2)+order_geo_pot
    if (order_geo_pot>0) then
      dim_gaupot = ((max_order_gaupot+1)*(max_order_gaupot+2)*(max_order_gaupot+3) &
                 -  order_geo_pot*(order_geo_pot+1)*(order_geo_pot+2))/6
    else
      dim_gaupot = (max_order_gaupot+1)*(max_order_gaupot+2)*(max_order_gaupot+3)/6
    end if
    ! allocates memory for the integrals from \fn(gaupot_geom)
    allocate(gaupot_pints(dim_gaupot), stat=ierr)
    if (ierr/=0)                                                             &
      call error_stop("prim_hgto_gaupot", "failed to allocate gaupot_pints", &
                      dim_gaupot)
    ! transfers Boys functions to geometric derivatives of Gaussian charge potential
    ! origin with zeroth order HGTOs
    call gaupot_geom(coord_bra, exponent_bra, coord_ket, exponent_ket, &
                     gaupot_origin, gaupot_expt, scal_const,           &
                     (/order_geo_pot,max_order_gaupot/), order_elec,   &
                     dim_gaupot, gaupot_pints)
    ! allocates memory for the integrals from \fn(nucpot_hket)
    min_hgto_ket = max(orders_hket_elec(1)-orders_hbra_mom(2),0)
    if (min_hgto_ket>0) then
      dim_hket = ((orders_hket_elec(2)+1)*(orders_hket_elec(2)+2) &
                  *(orders_hket_elec(2)+3)                        &
               -  min_hgto_ket*(min_hgto_ket+1)*(min_hgto_ket+2))/6
    else
      dim_hket = (orders_hket_elec(2)+1)*(orders_hket_elec(2)+2) &
               * (orders_hket_elec(2)+3)/6
    end if
    max_geo_hket = orders_hbra_mom(2)+order_geo_pot
    if (order_geo_pot>0) then
      dim_geo_hket = ((max_geo_hket+1)*(max_geo_hket+2)*(max_geo_hket+3) &
                   -  order_geo_pot*(order_geo_pot+1)*(order_geo_pot+2))/6
    else
      dim_geo_hket = (max_geo_hket+1)*(max_geo_hket+2)*(max_geo_hket+3)/6
    end if
    allocate(hket_pints(dim_hket,dim_geo_hket), stat=ierr)
    if (ierr/=0)                                                           &
      call error_stop("prim_hgto_gaupot", "failed to allocate hket_pints", &
                      dim_hket*dim_geo_hket)
    ! recovers the HGTOs on ket center
    call nucpot_hket((/min_hgto_ket,orders_hket_elec(2)/), &
                     (/order_geo_pot,max_geo_hket/),       &
                     coord_bra, exponent_bra,              &
                     coord_ket, exponent_ket,              &
                     dim_gaupot, gaupot_pints,             &
                     dim_hket, dim_geo_hket, hket_pints)
    deallocate(gaupot_pints)
    if (order_mom==0) then
      ! pure Gaussian charge potential
      if (order_elec==0) then
        ! recovers the HGTOs on bra center
        call nucpot_hbra(orders_hbra_mom, orders_hket_elec,  &
                         (/order_geo_pot,order_geo_pot/),    &
                         coord_bra, exponent_bra,            &
                         coord_ket, exponent_ket,            &
                         dim_hket, dim_geo_hket, hket_pints, &
                         dim_hgto_bra, dim_hgto_ket, num_geo_pot, hgto_pints)
        deallocate(hket_pints)
      ! Gaussian charge potential + electronic derivatives
      else
        ! allocates the intermediate integrals from \fn(nucpot_hbra)
        if (orders_hket_elec(1)==0) then
          dim_hket_zero = (orders_hket_elec(2)+1)*(orders_hket_elec(2)+2) &
                        * (orders_hket_elec(2)+3)/6
        else
          dim_hket_zero = ((orders_hket_elec(2)+1)*(orders_hket_elec(2)+2) &
                           *(orders_hket_elec(2)+3)                        &
                        -  orders_hket_elec(1)*(orders_hket_elec(1)+1)     &
                           *(orders_hket_elec(1)+2))/6
        end if
        allocate(hbra_pints(dim_hgto_bra,dim_hket_zero,num_geo_pot), stat=ierr)
        if (ierr/=0)                                                              &
          call error_stop("prim_hgto_gaupot", "failed to allocate hbra_pints/ne", &
                          dim_hgto_bra*dim_hket_zero*num_geo_pot)
        ! recovers the HGTOs on bra center
        call nucpot_hbra(orders_hbra_mom, orders_hket_elec,  &
                         (/order_geo_pot,order_geo_pot/),    &
                         coord_bra, exponent_bra,            &
                         coord_ket, exponent_ket,            &
                         dim_hket, dim_geo_hket, hket_pints, &
                         dim_hgto_bra, dim_hket_zero, num_geo_pot, hbra_pints)
        deallocate(hket_pints)
        ! recovers the electronic derivatives
        call scatter_multi_inner(dim_hgto_bra, dim_hket_zero, 1, num_geo_pot, &
                                 hbra_pints, orders_hgto_ket, order_elec,     &
                                 dim_hgto_ket, num_elec, hgto_pints)
        deallocate(hbra_pints)
      end if
    else
      ! Gaussian charge potential + Cartesian multipole moments
      if (order_elec==0) then
        ! allocates the intermediate integrals from \fn(nucpot_hbra)
        if (orders_hbra_mom(1)==0) then
          dim_hbra_zero = (orders_hbra_mom(2)+1)*(orders_hbra_mom(2)+2) &
                        * (orders_hbra_mom(2)+3)/6
        else
          dim_hbra_zero = ((orders_hbra_mom(2)+1)*(orders_hbra_mom(2)+2) &
                           *(orders_hbra_mom(2)+3)                       &
                        -  orders_hbra_mom(1)*(orders_hbra_mom(1)+1)     &
                           *(orders_hbra_mom(1)+2))/6
        end if
        allocate(hbra_pints(dim_hbra_zero,dim_hgto_ket,num_geo_pot), stat=ierr)
        if (ierr/=0)                                                              &
          call error_stop("prim_hgto_gaupot", "failed to allocate hbra_pints/nm", &
                          dim_hbra_zero*dim_hgto_ket*num_geo_pot)
        ! recovers the HGTOs on bra center
        call nucpot_hbra(orders_hbra_mom, orders_hket_elec,  &
                         (/order_geo_pot,order_geo_pot/),    &
                         coord_bra, exponent_bra,            &
                         coord_ket, exponent_ket,            &
                         dim_hket, dim_geo_hket, hket_pints, &
                         dim_hbra_zero, dim_hgto_ket, num_geo_pot, hbra_pints)
        deallocate(hket_pints)
        ! recovers the Cartesian multipole moments
        call london_mom_hgto(orders_hgto_bra, (/order_mom,order_mom/),    &
                             coord_bra, exponent_bra, dipole_origin,      &
                             1, dim_hbra_zero, dim_hgto_ket, num_geo_pot, &
                             hbra_pints, dim_hgto_bra, num_mom, hgto_pints)
        deallocate(hbra_pints)
      ! Gaussian charge potential + Cartesian multipole moments + electronic derivatives
      else
        ! allocates the intermediate integrals from \fn(nucpot_hbra)
        if (orders_hbra_mom(1)==0) then
          dim_hbra_zero = (orders_hbra_mom(2)+1)*(orders_hbra_mom(2)+2) &
                        * (orders_hbra_mom(2)+3)/6
        else
          dim_hbra_zero = ((orders_hbra_mom(2)+1)*(orders_hbra_mom(2)+2) &
                           *(orders_hbra_mom(2)+3)                       &
                        -  orders_hbra_mom(1)*(orders_hbra_mom(1)+1)     &
                           *(orders_hbra_mom(1)+2))/6
        end if
        if (orders_hket_elec(1)==0) then
          dim_hket_zero = (orders_hket_elec(2)+1)*(orders_hket_elec(2)+2) &
                        * (orders_hket_elec(2)+3)/6
        else
          dim_hket_zero = ((orders_hket_elec(2)+1)*(orders_hket_elec(2)+2) &
                           *(orders_hket_elec(2)+3)                        &
                        -  orders_hket_elec(1)*(orders_hket_elec(1)+1)     &
                           *(orders_hket_elec(1)+2))/6
        end if
        allocate(hbra_pints(dim_hbra_zero,dim_hket_zero,num_geo_pot), stat=ierr)
        if (ierr/=0)                                                               &
          call error_stop("prim_hgto_gaupot", "failed to allocate hbra_pints/nme", &
                          dim_hbra_zero*dim_hket_zero*num_geo_pot)
        ! recovers the HGTOs on bra center
        call nucpot_hbra(orders_hbra_mom, orders_hket_elec,  &
                         (/order_geo_pot,order_geo_pot/),    &
                         coord_bra, exponent_bra,            &
                         coord_ket, exponent_ket,            &
                         dim_hket, dim_geo_hket, hket_pints, &
                         dim_hbra_zero, dim_hket_zero, num_geo_pot, hbra_pints)
        deallocate(hket_pints)
        ! recovers the Cartesian multipole moments
        allocate(hmom_pints(dim_hgto_bra,dim_hket_zero,num_mom,num_geo_pot), stat=ierr)
        if (ierr/=0)                                                           &
          call error_stop("prim_hgto_gaupot", "failed to allocate hmom_pints", &
                          dim_hgto_bra*dim_hket_zero*num_mom*num_geo_pot)
        call london_mom_hgto(orders_hgto_bra, (/order_mom,order_mom/),     &
                             coord_bra, exponent_bra, dipole_origin,       &
                             1, dim_hbra_zero, dim_hket_zero, num_geo_pot, &
                             hbra_pints, dim_hgto_bra, num_mom, hmom_pints)
        deallocate(hbra_pints)
        ! recovers the electronic derivatives
        call scatter_multi_inner(dim_hgto_bra, dim_hket_zero, 1, num_mom*num_geo_pot, &
                                 hmom_pints, orders_hgto_ket, order_elec ,            &
                                 dim_hgto_ket, num_elec, hgto_pints)
        deallocate(hmom_pints)
      end if
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "prim_hgto_gaupot", STDOUT)
#endif
    return
  end subroutine prim_hgto_gaupot
