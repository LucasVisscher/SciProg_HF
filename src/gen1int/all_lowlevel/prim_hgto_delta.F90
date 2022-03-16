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
!!  This file contains the recurrence relations of Dirac delta function (multiplied
!!  by Cartesian multipole moment) integrals using primitive Hermite Gaussians.
!!
!!  2012-03-17, Bin Gao:
!!  * rewrites with new recurrence relations
!!
!!  2010-10-10, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief recurrence relations of Dirac delta function (multiplied by Cartesian
  !>        multipole moment) integrals using primitive Hermite Gaussians
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param orders_hgto_bra is the range of orders of Hermite Gaussians on bra center
  !> \param coord_bra contains the coordinates of bra center
  !> \param exponent_bra is the exponent of primitive Gaussian of bra center
  !> \param orders_hgto_ket is the range of orders of Hermite Gaussians on ket center
  !> \param coord_ket contains the coordinates of ket center
  !> \param exponent_ket is the exponent of primitive Gaussian of ket center
  !> \param delta_origin contains the coordinates of Dirac delta function origin
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scale constant for Dirac delta function
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
  subroutine prim_hgto_delta(orders_hgto_bra, coord_bra, exponent_bra, &
                             orders_hgto_ket, coord_ket, exponent_ket, &
                             order_elec, delta_origin, dipole_origin,  &
                             scal_const, order_mom, order_geo_pot,     &
                             dim_hgto_bra, dim_hgto_ket, num_elec,     &
                             num_mom, num_geo_pot, hgto_pints)
    use xkind
    implicit none
    integer, intent(in) :: orders_hgto_bra(2)
    real(REALK), intent(in) :: coord_bra(3)
    real(REALK), intent(in) :: exponent_bra
    integer, intent(in) :: orders_hgto_ket(2)
    real(REALK), intent(in) :: coord_ket(3)
    real(REALK), intent(in) :: exponent_ket
    integer, intent(in) :: order_elec
    real(REALK), intent(in) :: delta_origin(3)
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
!f2py intent(in) :: delta_origin
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
    integer max_dim_geo          !maximum dimension of geometric derivatives on Dirac delta function
    integer min_delta_geom       !minimum order of geometric derivatives on Dirac delta function ...
                                 !from \fn(delta_geom)
    integer dim_geo_pot          !dimension of integrals from \fn(delta_geom)
    real(REALK), allocatable :: geo_pot_pints(:)
                                 !integrals from \fn(delta_geom)
    integer min_geo_hbra         !minimum order of geometric derivatives on Dirac delta function ...
                                 !after recovering HGTOs on bra center
    integer dim_geo_hbra         !dimension of geometric derivatives on Dirac delta function ...
                                 !after recovering HGTOs on bra center
    real(REALK), allocatable :: hbra_pints(:,:)
                                 !integrals after recovering HGTOs on bra center
    integer min_geo_hket         !minimum order of geometric derivatives on Dirac delta function ...
                                 !after recovering HGTOs on ket center
    integer dim_hket             !dimension of HGTOs on ket center after recovering HGTOs on ket center
    integer dim_geo_hket         !dimension of geometric derivatives on Dirac delta function ...
                                 !after recovering HGTOs on ket center
    real(REALK), allocatable :: hket_pints(:,:,:)
                                 !integrals after recovering HGTOs on ket center
    real(REALK), allocatable :: elec_pints(:,:,:,:)
                                 !integrals after recovering electronic derivatives
    integer ierr                 !error information
#if defined(XTIME)
    real(REALK) curr_time        !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! sets the orders of HGTOs on ket center with electronic derivatives
    orders_hket_elec(1) = orders_hgto_ket(1)+order_elec
    orders_hket_elec(2) = orders_hgto_ket(2)+order_elec
    ! sets the minimum order of geometric derivatives on Dirac delta function
    ! after recovering HGTOs on ket center
    min_geo_hket = max(0,order_geo_pot-order_mom)
    ! sets the minimum order of geometric derivatives on Dirac delta function
    ! after recovering HGTOs on bra center
    min_geo_hbra = max(0,min_geo_hket-orders_hket_elec(2))
    ! sets the minimum order of geometric derivatives on Dirac delta function from \fn(delta_geom)
    min_delta_geom = max(0,min_geo_hbra-orders_hgto_bra(2))
    ! sets the maximum dimension of geometric derivatives on Dirac delta function
    max_dim_geo = (order_geo_pot+1)*(order_geo_pot+2)*(order_geo_pot+3)/6
    ! allocates memory for integrals from \fn(delta_geom)
    if (min_delta_geom==0) then
      dim_geo_pot = max_dim_geo
    else
      dim_geo_pot = max_dim_geo &
                  - min_delta_geom*(min_delta_geom+1)*(min_delta_geom+2)/6
    end if
    allocate(geo_pot_pints(dim_geo_pot), stat=ierr)
    if (ierr/=0)                                                             &
      call error_stop("prim_hgto_delta", "failed to allocate geo_pot_pints", &
                      dim_geo_pot)
    ! recovers the geometric derivatives on Dirac delta function using \fn(delta_geom)
    call delta_geom(coord_bra, exponent_bra, coord_ket, exponent_ket,           &
                    delta_origin, scal_const, (/min_delta_geom,order_geo_pot/), &
                    order_elec, dim_geo_pot, geo_pot_pints)
    ! allocates memory for integrals from \fn(delta_hket) on bra center
    if (min_geo_hbra==0) then
      dim_geo_hbra = max_dim_geo
    else
      dim_geo_hbra = max_dim_geo &
                   - min_geo_hbra*(min_geo_hbra+1)*(min_geo_hbra+2)/6
    end if
    allocate(hbra_pints(dim_hgto_bra,dim_geo_hbra), stat=ierr)
    if (ierr/=0)                                                          &
      call error_stop("prim_hgto_delta", "failed to allocate hbra_pints", &
                      dim_hgto_bra*dim_geo_hbra)
    ! recovers the HGTOs on bra center using \fn(delta_hket)
    call delta_hket(orders_hgto_bra, (/min_geo_hbra,order_geo_pot/),       &
                    coord_bra, exponent_bra, delta_origin, 1, dim_geo_pot, &
                    geo_pot_pints, dim_hgto_bra, dim_geo_hbra, hbra_pints)
    deallocate(geo_pot_pints)
    if (order_mom==0) then
      ! pure Dirac delta function integrals
      if (order_elec==0) then
        ! recovers the HGTOs on ket center using \fn(delta_hket)
        call delta_hket(orders_hgto_ket, (/order_geo_pot,order_geo_pot/),    &
                        coord_ket, exponent_ket, delta_origin, dim_hgto_bra, &
                        dim_geo_hbra, hbra_pints, dim_hgto_ket, num_geo_pot, &
                        hgto_pints)
        deallocate(hbra_pints)
      ! Dirac delta function + electronic derivatives
      else
        ! allocates memory for integrals from \fn(delta_hket) on ket center
        if (orders_hket_elec(1)==0) then
          dim_hket = (orders_hket_elec(2)+1)*(orders_hket_elec(2)+2) &
                   * (orders_hket_elec(2)+3)/6
        else
          dim_hket = ((orders_hket_elec(2)+1)*(orders_hket_elec(2)+2) &
                   *  (orders_hket_elec(2)+3)                         &
                   -  orders_hket_elec(1)*(orders_hket_elec(1)+1)     &
                   *  (orders_hket_elec(1)+2))/6
        end if
        allocate(hket_pints(dim_hgto_bra,dim_hket,num_geo_pot), stat=ierr)
        if (ierr/=0)                                                          &
          call error_stop("prim_hgto_delta", "failed to allocate hket_pints", &
                          dim_hgto_bra*dim_hket*num_geo_pot)
        ! recovers the HGTOs on ket center using \fn(delta_hket)
        call delta_hket(orders_hket_elec, (/order_geo_pot,order_geo_pot/),   &
                        coord_ket, exponent_ket, delta_origin, dim_hgto_bra, &
                        dim_geo_hbra, hbra_pints, dim_hket, num_geo_pot,     &
                        hket_pints)
        deallocate(hbra_pints)
        ! recovers the electronic derivatives
        call scatter_multi_inner(dim_hgto_bra, dim_hket, 1, num_geo_pot,  &
                                 hket_pints, orders_hgto_ket, order_elec, &
                                 dim_hgto_ket, num_elec, hgto_pints)
        deallocate(hket_pints)
      end if
    else
      ! Dirac delta function + Cartesian multipole moments
      if (order_elec==0) then
        ! allocates memory for integrals from \fn(delta_hket) on ket center
        if (min_geo_hket==0) then
          dim_geo_hket = max_dim_geo
        else
          dim_geo_hket = max_dim_geo &
                       - min_geo_hket*(min_geo_hket+1)*(min_geo_hket+2)/6
        end if
        allocate(hket_pints(dim_hgto_bra,dim_hgto_ket,dim_geo_hket), stat=ierr)
        if (ierr/=0)                                                          &
          call error_stop("prim_hgto_delta", "failed to allocate hket_pints", &
                          dim_hgto_bra*dim_hgto_ket*dim_geo_hket)
        ! recovers the HGTOs on ket center using \fn(delta_hket)
        call delta_hket(orders_hgto_ket, (/min_geo_hket,order_geo_pot/),      &
                        coord_ket, exponent_ket, delta_origin, dim_hgto_bra,  &
                        dim_geo_hbra, hbra_pints, dim_hgto_ket, dim_geo_hket, &
                        hket_pints)
        deallocate(hbra_pints)
        ! recovers the Cartesian multipole moments
        call delta_moment(order_mom, order_geo_pot, dipole_origin,  &
                          delta_origin, dim_hgto_bra, dim_hgto_ket, &
                          dim_geo_hket, hket_pints, num_mom,        &
                          num_geo_pot, hgto_pints)
        deallocate(hket_pints)
      ! Dirac delta function + Cartesian multipole moments + electronic derivatives
      else
        ! allocates memory for integrals from \fn(delta_hket) on ket center
        if (orders_hket_elec(1)==0) then
          dim_hket = (orders_hket_elec(2)+1)*(orders_hket_elec(2)+2) &
                   * (orders_hket_elec(2)+3)/6
        else
          dim_hket = ((orders_hket_elec(2)+1)*(orders_hket_elec(2)+2) &
                   *  (orders_hket_elec(2)+3)                         &
                   -  orders_hket_elec(1)*(orders_hket_elec(1)+1)     &
                   *  (orders_hket_elec(1)+2))/6
        end if
        if (min_geo_hket==0) then
          dim_geo_hket = max_dim_geo
        else
          dim_geo_hket = max_dim_geo &
                       - min_geo_hket*(min_geo_hket+1)*(min_geo_hket+2)/6
        end if
        allocate(hket_pints(dim_hgto_bra,dim_hket,dim_geo_hket), stat=ierr)
        if (ierr/=0)                                                          &
          call error_stop("prim_hgto_delta", "failed to allocate hket_pints", &
                          dim_hgto_bra*dim_hket*dim_geo_hket)
        ! recovers the HGTOs on ket center using \fn(delta_hket)
        call delta_hket(orders_hket_elec, (/min_geo_hket,order_geo_pot/),     &
                        coord_ket, exponent_ket, delta_origin, dim_hgto_bra,  &
                        dim_geo_hbra, hbra_pints, dim_hket, dim_geo_hket,     &
                        hket_pints)
        deallocate(hbra_pints)
        ! recovers the electronic derivatives
        allocate(elec_pints(dim_hgto_bra,dim_hgto_ket,num_elec,dim_geo_hket), stat=ierr)
        if (ierr/=0)                                                          &
          call error_stop("prim_hgto_delta", "failed to allocate elec_pints", &
                          dim_hgto_bra*dim_hgto_ket*num_elec*dim_geo_hket)
        call scatter_multi_inner(dim_hgto_bra, dim_hket, 1, dim_geo_hket, &
                                 hket_pints, orders_hgto_ket, order_elec, &
                                 dim_hgto_ket, num_elec, elec_pints)
        deallocate(hket_pints)
        ! recovers the Cartesian multipole moments
        call delta_moment(order_mom, order_geo_pot, dipole_origin, delta_origin, &
                          dim_hgto_bra, dim_hgto_ket*num_elec, dim_geo_hket,     &
                          elec_pints, num_mom, num_geo_pot, hgto_pints)
        deallocate(elec_pints)
      end if
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "prim_hgto_delta", STDOUT)
#endif
    return
  end subroutine prim_hgto_delta
