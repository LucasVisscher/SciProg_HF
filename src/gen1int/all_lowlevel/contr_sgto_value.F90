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
!!  This file calculates the values of contracted spheical Gaussians.
!!
!!  2012-03-10, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief calculates the values of contracted spherical Gaussians
  !> \author Bin Gao
  !> \date 2012-03-10
  !> \param coord_bra contains the coordinates of bra center
  !> \param angular_bra is the angular number of bra center (s=0, p=1, d=2, ...)
  !> \param num_prim_bra is the number of primitive Gaussians of bra center
  !> \param exponent_bra contains the exponents of primitive Gaussians of bra center
  !> \param num_contr_bra is the number of contractions of bra center
  !> \param contr_coef_bra contains the contraction coefficients of bra center
  !> \param order_geo_bra is the order of geometric derivatives with respect to bra center
  !> \param num_points is the number of grid points
  !> \param grid_points contains the coordinates of grid points
  !> \param num_sgto_bra is the number of spherical GTOs on bra center,
  !>        equals to \f$2\var(angular_bra)+1\f$
  !> \param num_derv is the number of geometric derivatives
  !> \return contr_value contains the contracted SGTOs
  subroutine contr_sgto_value(coord_bra, angular_bra, num_prim_bra,        &
                              exponent_bra, num_contr_bra, contr_coef_bra, &
                              order_geo_bra, num_points, grid_points,      &
                              num_sgto_bra, num_derv, contr_value)
    use xkind
    implicit none
    real(REALK), intent(in) :: coord_bra(3)
    integer, intent(in) :: angular_bra
    integer, intent(in) :: num_prim_bra
    real(REALK), intent(in) :: exponent_bra(num_prim_bra)
    integer, intent(in) :: num_contr_bra
    real(REALK), intent(in) :: contr_coef_bra(num_contr_bra,num_prim_bra)
    integer, intent(in) :: order_geo_bra
    integer, intent(in) :: num_points
    real(REALK), intent(in) :: grid_points(3,num_points)
    integer, intent(in) :: num_sgto_bra
    integer, intent(in) :: num_derv
    real(REALK), intent(out) :: contr_value(num_sgto_bra,num_contr_bra, &
                                            num_points, num_derv)
!f2py intent(in) :: coord_bra
!f2py intent(in) :: angular_bra
!f2py intent(hide) :: num_prim_bra
!f2py intent(in) :: exponent_bra
!f2py intent(hide) :: num_contr_bra
!f2py intent(in) :: contr_coef_bra
!f2py depend(num_prim_bra) :: contr_coef_bra
!f2py intent(in) :: order_geo_bra
!f2py intent(hide) :: num_points
!f2py intent(in) :: grid_points
!f2py intent(in) :: num_sgto_bra
!f2py intent(in) :: num_derv
!f2py intent(out) :: contr_value
!f2py depend(num_sgto_bra) :: contr_value
!f2py depend(num_contr_bra) :: contr_value
!f2py depend(num_points) :: contr_value
!f2py depend(num_derv) :: contr_value
    integer orders_hgto_bra(2)  !range of orders of Hermite Gaussians on bra center
    integer dim_hgto_bra        !dimension of Hermite Gaussians on bra center including geometric derivatives
    integer num_hgto_bra        !number of Hermite Gaussians on bra center
    integer num_geo_bra         !number of geometric derivatives on bra center
    integer num_opt_geo         !number of operators and geometric derivatives
    real(REALK), allocatable :: hgto_value(:,:,:)       !primitive or contracted HGTOs
    real(REALK), allocatable :: geom_hgto_value(:,:,:)  !geometric derivatives of primitive HGTOs
    real(REALK) scal_geo_bra    !scale constant when recovering geometric derivatives on bra center
    integer iprim               !incremental recorder over primitives
    integer ierr                !error information
#if defined(XTIME)
    real(REALK) curr_time       !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
#if defined(DEBUG)
    ! dumps the contracted GTOs and derivatives to check
    call dump_gto_pd("SGTO-bra", 0, coord_bra, angular_bra,     &
                     num_prim_bra, exponent_bra, num_contr_bra, &
                     contr_coef_bra, 0, 0, order_geo_bra, STDOUT)
#endif
    ! sets the orders of Hermite Gaussians on bra center
    orders_hgto_bra(1) = order_geo_bra+angular_bra
    orders_hgto_bra(2) = orders_hgto_bra(1)
    ! sets the dimension of primitive Hermite Gaussians including geometric derivatives
    dim_hgto_bra = (orders_hgto_bra(1)+1)*(orders_hgto_bra(1)+2)/2
    ! computes the number of primitive Hermite Gaussians
    num_hgto_bra = (angular_bra+1)*(angular_bra+2)/2
    ! computes the number of geometric derivatives
    num_geo_bra = (order_geo_bra+1)*(order_geo_bra+2)/2
    num_opt_geo = num_points*num_geo_bra
    ! allocates the memory for primitive HGTOs
    allocate(geom_hgto_value(num_hgto_bra,num_opt_geo,num_prim_bra), stat=ierr)
    if (ierr/=0)                                            &
      call error_stop("contr_sgto_value",                   &
                      "failed to allocate geom_hgto_value", &
                      num_hgto_bra*num_opt_geo*num_prim_bra)
#if defined(DEBUG)
    write(STDOUT,100) "size of geometric derivatives of primitive HGTOs:", &
                      size(geom_hgto_value)
#endif
    if (order_geo_bra>0) then
      allocate(hgto_value(dim_hgto_bra,num_points,1), stat=ierr)
      if (ierr/=0)                                       &
        call error_stop("contr_sgto_value",              &
                        "failed to allocate hgto_value", &
                        dim_hgto_bra*num_points)
#if defined(DEBUG)
      write(STDOUT,100) "size of primitive HGTOs:", size(hgto_value)
#endif
      do iprim = 1, num_prim_bra
        ! sets the scale constant when recovering geometric derivatives on bra center
        scal_geo_bra = (exponent_bra(iprim)+exponent_bra(iprim))**order_geo_bra
        ! computes the primitive HGTOs
        call prim_hgto_value(orders_hgto_bra, coord_bra, exponent_bra(iprim), &
                             num_points, grid_points, dim_hgto_bra, hgto_value)
        ! scales the primitive HGTOs
        hgto_value = scal_geo_bra*hgto_value
        ! recovers the geometric derivatives on bra center
        call scatter_single(1, dim_hgto_bra, num_points, 1, hgto_value, &
                            angular_bra, order_geo_bra, num_hgto_bra,   &
                            num_geo_bra, geom_hgto_value(:,:,iprim))
      end do
      ! cleans
      deallocate(hgto_value)
    else
      do iprim = 1, num_prim_bra
        ! computes the primitive HGTOs
        call prim_hgto_value(orders_hgto_bra, coord_bra, exponent_bra(iprim), &
                             num_points, grid_points, dim_hgto_bra,           &
                             geom_hgto_value(:,:,iprim))
      end do
    end if
    ! allocates the memory for geometric derivatives of contracted HGTOs
    allocate(hgto_value(num_hgto_bra,num_contr_bra,num_opt_geo), stat=ierr)
    if (ierr/=0)                                       &
      call error_stop("contr_sgto_value",              &
                      "failed to allocate hgto_value", &
                      num_hgto_bra*num_contr_bra*num_opt_geo)
#if defined(DEBUG)
    write(STDOUT,100) "size of geometric derivatives of contracted HGTOs:", &
                      size(hgto_value)
#endif
    ! constructs contracted HGTOs in the order of
    ! hgto_value(num_hgto_bra,num_contr_bra,num_points,num_geo_bra)
    call const_contr_gto(num_contr_bra, num_prim_bra, contr_coef_bra, &
                         num_hgto_bra, num_opt_geo, geom_hgto_value, hgto_value)
    deallocate(geom_hgto_value)
    ! transforms the HGTOs to SGTOs on bra center
    call hgto_to_sgto_oop(angular_bra, 1, num_hgto_bra, num_contr_bra*num_opt_geo, &
                          hgto_value, num_sgto_bra, contr_value)
    deallocate(hgto_value)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "contr_sgto_value", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("contr_sgto_value>> ",A,I8,1X,A)
#endif
  end subroutine contr_sgto_value
