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
!!  This file calculates the values of contracted Cartesian Gaussians.
!!
!!  2012-03-10, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief calculates the values of contracted Cartesian Gaussians
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
  !> \param num_cgto_bra is the number of Cartesian GTOs on bra center,
  !>        equals to \f$(\var(angular_bra)+1)(\var(angular_bra)+2)/2\f$
  !> \param num_derv is the number of geometric derivatives
  !> \return contr_value contains the contracted CGTOs
  subroutine contr_cgto_value(coord_bra, angular_bra, num_prim_bra,        &
                              exponent_bra, num_contr_bra, contr_coef_bra, &
                              order_geo_bra, num_points, grid_points,      &
                              num_cgto_bra,num_derv, contr_value)
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
    integer, intent(in) :: num_cgto_bra
    integer, intent(in) :: num_derv
    real(REALK), intent(out) :: contr_value(num_cgto_bra,num_contr_bra, &
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
!f2py intent(in) :: num_cgto_bra
!f2py intent(in) :: num_derv
!f2py intent(out) :: contr_value
!f2py depend(num_cgto_bra) :: contr_value
!f2py depend(num_contr_bra) :: contr_value
!f2py depend(num_points) :: contr_value
!f2py depend(num_derv) :: contr_value
    integer orders_hgto_bra(2)  !range of orders of Hermite Gaussians on bra center
    integer dim_hgto_bra        !dimension of Hermite Gaussians on bra center including geometric derivatives
    integer num_geo_bra         !number of geometric derivatives on bra center
    integer num_opt_geo         !number of operators and geometric derivatives
    real(REALK), allocatable :: hgto_value(:,:)    !primitive HGTOs
    real(REALK), allocatable :: cgto_value(:,:,:)  !primitive CGTOs
    integer iprim               !incremental recorder over primitives
    integer ierr                !error information
#if defined(XTIME)
    real(REALK) curr_time       !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
#if defined(DEBUG)
    ! dumps the contracted GTOs and derivatives to check
    call dump_gto_pd("CGTO-bra", 0, coord_bra, angular_bra,     &
                     num_prim_bra, exponent_bra, num_contr_bra, &
                     contr_coef_bra, 0, 0, order_geo_bra, STDOUT)
#endif
    ! computes the minimum and maximum orders of Hermite Gaussians on bra center
    orders_hgto_bra(1) = order_geo_bra+mod(angular_bra,2)
    orders_hgto_bra(2) = order_geo_bra+angular_bra
    ! computes the number of primitive Hermite Gaussians used in recurrence relations
    if (orders_hgto_bra(1)==0) then
      dim_hgto_bra = (orders_hgto_bra(2)+1)*(orders_hgto_bra(2)+2) &
                   * (orders_hgto_bra(2)+3)/6
    else
      dim_hgto_bra = ((orders_hgto_bra(2)+1)*(orders_hgto_bra(2)+2) &
                      *(orders_hgto_bra(2)+3)                       &
                   -  orders_hgto_bra(1)*(orders_hgto_bra(1)+1)     &
                      *(orders_hgto_bra(1)+2))/6
    end if
    ! computes the number of geometric derivatives
    num_geo_bra = (order_geo_bra+1)*(order_geo_bra+2)/2
    num_opt_geo = num_points*num_geo_bra
    ! allocates the memory for primitive HGTOs
    allocate(hgto_value(dim_hgto_bra,num_points), stat=ierr)
    if (ierr/=0)                                       &
      call error_stop("contr_cgto_value",              &
                      "failed to allocate hgto_value", &
                      dim_hgto_bra*num_points)
    allocate(cgto_value(num_cgto_bra,num_opt_geo,num_prim_bra), stat=ierr)
    if (ierr/=0)                                       &
      call error_stop("contr_cgto_value",              &
                      "failed to allocate cgto_value", &
                      num_cgto_bra*num_opt_geo*num_prim_bra)
#if defined(DEBUG)
    write(STDOUT,100) "size of primitive HGTOs:", size(hgto_value)
    write(STDOUT,100) "size of primitive CGTOs:", size(cgto_value)
#endif
    ! computes the values of primitive CGTOs
    do iprim = 1, num_prim_bra
      call prim_hgto_value(orders_hgto_bra, coord_bra, exponent_bra(iprim), &
                           num_points, grid_points, dim_hgto_bra, hgto_value)
      ! transforms Hermite Gaussians on bra center
      call hgto_to_cgto(angular_bra, order_geo_bra, exponent_bra(iprim), &
                        1, dim_hgto_bra, num_points, 1, hgto_value,      &
                        num_cgto_bra, num_geo_bra, cgto_value(:,:,iprim))
    end do
    ! cleans
    deallocate(hgto_value)
    ! constructs contracted CGTOs in the order of
    ! contr_value(num_cgto_bra,num_contr_bra,num_points,num_geo_bra)
    call const_contr_gto(num_contr_bra, num_prim_bra, contr_coef_bra, &
                         num_cgto_bra, num_opt_geo, cgto_value, contr_value)
    ! cleans
    deallocate(cgto_value)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "contr_cgto_value", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("contr_cgto_value>> ",A,I8,1X,A)
#endif
  end subroutine contr_cgto_value
