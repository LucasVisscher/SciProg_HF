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
!!  This file calculates the overlap distribution of two contracted spherical Gaussians.
!!
!!  2012-01-24, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief calculates the overlap distribution of two contracted spherical Gaussians
  !> \author Bin Gao
  !> \date 2012-01-24
  !> \param idx_bra is the atomic index of bra center
  !> \param coord_bra contains the coordinates of bra center
  !> \param angular_bra is the angular number of bra center (s=0, p=1, d=2, ...)
  !> \param num_prim_bra is the number of primitive Gaussians of bra center
  !> \param exponent_bra contains the exponents of primitive Gaussians of bra center
  !> \param num_contr_bra is the number of contractions of bra center
  !> \param contr_coef_bra contains the contraction coefficients of bra center
  !> \param idx_ket is the atomic index of ket center
  !> \param coord_ket contains the coordinates of ket center
  !> \param angular_ket is the angular number of ket center
  !> \param num_prim_ket is the number of primitive Gaussians of ket center
  !> \param exponent_ket contains the exponents of primitive Gaussians of ket center
  !> \param num_contr_ket is the number of contractions of ket center
  !> \param contr_coef_ket contains the contraction coefficients of ket center
  !> \param order_geo_bra is the order of geometric derivatives with respect to bra center
  !> \param order_geo_ket is the order of geometric derivatives with respect to ket center
  !> \param num_cents is the number of differentiated centers
  !> \param idx_cent contains the indices of differentiated centers
  !> \param order_cent contains the order of derivatives of the differentiated centers
  !> \param num_points is the number of grid points
  !> \param grid_points contains the coordinates of grid points
  !> \param num_sgto_bra is the number of spherical GTOs on bra center,
  !>        equals to \f$2\var(angular_bra)+1\f$
  !> \param num_sgto_ket is the number of spherical GTOs on ket center,
  !>        equals to \f$2\var(angular_ket)+1\f$
  !> \param num_derv is the number of geometric derivatives
  !> \return contr_odist contains the overlap distribution
  subroutine contr_sgto_odist(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                              exponent_bra, num_contr_bra, contr_coef_bra,   &
                              idx_ket, coord_ket, angular_ket, num_prim_ket, &
                              exponent_ket, num_contr_ket, contr_coef_ket,   &
                              order_geo_bra, order_geo_ket, num_cents,       &
                              idx_cent, order_cent, num_points, grid_points, &
                              num_sgto_bra, num_sgto_ket, num_derv, contr_odist)
    use xkind
    implicit none
    integer, intent(in) :: idx_bra
    real(REALK), intent(in) :: coord_bra(3)
    integer, intent(in) :: angular_bra
    integer, intent(in) :: num_prim_bra
    real(REALK), intent(in) :: exponent_bra(num_prim_bra)
    integer, intent(in) :: num_contr_bra
    real(REALK), intent(in) :: contr_coef_bra(num_contr_bra,num_prim_bra)
    integer, intent(in) :: idx_ket
    real(REALK), intent(in) :: coord_ket(3)
    integer, intent(in) :: angular_ket
    integer, intent(in) :: num_prim_ket
    real(REALK), intent(in) :: exponent_ket(num_prim_ket)
    integer, intent(in) :: num_contr_ket
    real(REALK), intent(in) :: contr_coef_ket(num_contr_ket,num_prim_ket)
    integer, intent(in) :: order_geo_bra
    integer, intent(in) :: order_geo_ket
    integer, intent(in) :: num_cents
    integer, intent(in) :: idx_cent(num_cents)
    integer, intent(in) :: order_cent(num_cents)
    integer, intent(in) :: num_points
    real(REALK), intent(in) :: grid_points(3,num_points)
    integer, intent(in) :: num_sgto_bra
    integer, intent(in) :: num_sgto_ket
    integer, intent(in) :: num_derv
    real(REALK), intent(out) :: contr_odist(num_sgto_bra,num_contr_bra, &
                                            num_sgto_ket,num_contr_ket, &
                                            num_points,num_derv)
!f2py intent(in) :: idx_bra
!f2py intent(in) :: coord_bra
!f2py intent(in) :: angular_bra
!f2py intent(hide) :: num_prim_bra
!f2py intent(in) :: exponent_bra
!f2py intent(hide) :: num_contr_bra
!f2py intent(in) :: contr_coef_bra
!f2py depend(num_prim_bra) :: contr_coef_bra
!f2py intent(in) :: idx_ket
!f2py intent(in) :: coord_ket
!f2py intent(in) :: angular_ket
!f2py intent(hide) :: num_prim_ket
!f2py intent(in) :: exponent_ket
!f2py intent(hide) :: num_contr_ket
!f2py intent(in) :: contr_coef_ket
!f2py depend(num_prim_ket) :: contr_coef_ket
!f2py intent(in) :: order_geo_bra
!f2py intent(in) :: order_geo_ket
!f2py intent(hide) :: num_cents
!f2py intent(in) :: idx_cent
!f2py intent(in) :: order_cent
!f2py depend(num_cents) :: order_cent
!f2py intent(hide) :: num_points
!f2py intent(in) :: grid_points
!f2py intent(in) :: num_sgto_bra
!f2py intent(in) :: num_sgto_ket
!f2py intent(in) :: num_derv
!f2py intent(out) :: contr_odist
!f2py depend(num_sgto_bra) :: contr_odist
!f2py depend(num_contr_bra) :: contr_odist
!f2py depend(num_sgto_ket) :: contr_odist
!f2py depend(num_contr_ket) :: contr_odist
!f2py depend(num_points) :: contr_odist
!f2py depend(num_derv) :: contr_odist
    integer order_part_cent(2)  !orders of partial geometric derivatives of bra, ket and dipole origin centers
    logical zero_odist          !indicates if the overlap distribution is zero
    logical scatter_deriv       !indicates if scattering the geometric derivatives
    integer seq_part_geo(2)     !sequence of operator, bra and ket centers for partial geometric derivatives
    integer dim_odist           !dimension of overlap distribution using contracted GTOs
    integer dim_geo_bra         !dimension of geometric derivatives on bra center before scattering
    integer dim_geo_ket         !dimension of geometric derivatives on ket center before scattering
    integer num_part_derv       !number of derivatives before scattering
    integer num_geo_bra         !number of partial geometric derivatives on bra center after scattering
    integer num_geo_ket         !number of partial geometric derivatives on ket center after scattering
    integer num_tot_geo         !number of total geometric derivatives after scattering
    real(REALK), allocatable :: part_odist(:,:,:)
                                !contains the derivatives of overlap distribution before scattering
    integer ierr                !error information
#if defined(XTIME)
    real(REALK) curr_time       !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
#if defined(DEBUG)
    ! dumps the contracted GTOs and derivatives to check
    call dump_gto_pd("SGTO-bra", idx_bra, coord_bra, angular_bra, &
                     num_prim_bra, exponent_bra, num_contr_bra,   &
                     contr_coef_bra, 0, 0, order_geo_bra, STDOUT)
    call dump_gto_pd("SGTO-ket", idx_ket, coord_ket, angular_ket, &
                     num_prim_ket, exponent_ket, num_contr_ket,   &
                     contr_coef_ket, 0, 0, order_geo_ket, STDOUT)
    call dump_total_derv(0, 0, num_cents, idx_cent, order_cent, STDOUT)
#endif
    ! no total geometric derivatives
    if (num_cents==0) then
      call contr_sgto_odist_recurr(coord_bra, angular_bra, num_prim_bra,        &
                                   exponent_bra, num_contr_bra, contr_coef_bra, &
                                   coord_ket, angular_ket, num_prim_ket,        &
                                   exponent_ket, num_contr_ket, contr_coef_ket, &
                                   order_geo_bra, order_geo_ket, num_points,    &
                                   grid_points, num_sgto_bra, num_sgto_ket,     &
                                   num_derv, contr_odist)
    else
      ! sets the orders of partial geometric derivatives
      order_part_cent(1) = order_geo_bra
      order_part_cent(2) = order_geo_ket
      call geom_part_zero_param(num_cents, idx_cent, order_cent,      &
                                (/idx_bra,idx_ket/), order_part_cent, &
                                zero_odist, scatter_deriv, seq_part_geo)
      if (zero_odist) then
        contr_odist = 0.0_REALK
      ! we need to scatter the partial and total geometric derivatives into appropriate places
      else if (scatter_deriv) then
        ! allocates memory for the derivatives of overlap distribution before scattering
        dim_odist = size(contr_odist)/num_derv
        dim_geo_bra = (order_part_cent(1)+1)*(order_part_cent(1)+2)/2
        dim_geo_ket = (order_part_cent(2)+1)*(order_part_cent(2)+2)/2
        allocate(part_odist(dim_odist,dim_geo_bra,dim_geo_ket), stat=ierr)
        if (ierr/=0)                                       &
          call error_stop("contr_sgto_odist",              &
                          "failed to allocate part_odist", &
                          dim_odist*dim_geo_bra*dim_geo_ket)
        num_part_derv = dim_geo_bra*dim_geo_ket
        call contr_sgto_odist_recurr(coord_bra, angular_bra, num_prim_bra,        &
                                     exponent_bra, num_contr_bra, contr_coef_bra, &
                                     coord_ket, angular_ket, num_prim_ket,        &
                                     exponent_ket, num_contr_ket, contr_coef_ket, &
                                     order_part_cent(1), order_part_cent(2),      &
                                     num_points, grid_points, num_sgto_bra,       &
                                     num_sgto_ket, num_part_derv, part_odist)
        ! scatters the partial and total geometric derivatives
        num_geo_bra = (order_geo_bra+1)*(order_geo_bra+2)/2
        num_geo_ket = (order_geo_ket+1)*(order_geo_ket+2)/2
        num_tot_geo = num_derv/(num_geo_bra*num_geo_ket)
        call geom_part_zero_scatter(num_cents, order_cent, seq_part_geo(1:num_cents), &
                                    order_geo_bra, order_geo_ket, dim_odist,          &
                                    dim_geo_bra, dim_geo_ket, 1, part_odist,          &
                                    num_geo_bra, num_geo_ket, num_tot_geo, contr_odist)
        deallocate(part_odist)
      else
        call contr_sgto_odist_recurr(coord_bra, angular_bra, num_prim_bra,        &
                                     exponent_bra, num_contr_bra, contr_coef_bra, &
                                     coord_ket, angular_ket, num_prim_ket,        &
                                     exponent_ket, num_contr_ket, contr_coef_ket, &
                                     order_part_cent(1), order_part_cent(2),      &
                                     num_points, grid_points, num_sgto_bra,       &
                                     num_sgto_ket, num_derv, contr_odist)
      end if
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "contr_sgto_odist", STDOUT)
#endif
    return
  end subroutine contr_sgto_odist

  !> \brief recurrence relations of overlap distribution between two contracted spherical Gaussians
  !> \author Bin Gao
  !> \date 2012-01-24
  !> \param coord_bra contains the coordinates of bra center
  !> \param angular_bra is the angular number of bra center (s=0, p=1, d=2, ...)
  !> \param num_prim_bra is the number of primitive Gaussians of bra center
  !> \param exponent_bra contains the exponents of primitive Gaussians of bra center
  !> \param num_contr_bra is the number of contractions of bra center
  !> \param contr_coef_bra contains the contraction coefficients of bra center
  !> \param coord_ket contains the coordinates of ket center
  !> \param angular_ket is the angular number of ket center
  !> \param num_prim_ket is the number of primitive Gaussians of ket center
  !> \param exponent_ket contains the exponents of primitive Gaussians of ket center
  !> \param num_contr_ket is the number of contractions of ket center
  !> \param contr_coef_ket contains the contraction coefficients of ket center
  !> \param order_geo_bra is the order of geometric derivatives with respect to bra center
  !> \param order_geo_ket is the order of geometric derivatives with respect to ket center
  !> \param num_points is the number of grid points
  !> \param grid_points contains the coordinates of grid points
  !> \param num_sgto_bra is the number of spherical GTOs on bra center,
  !>        equals to \f$2\var(angular_bra)+1\f$
  !> \param num_sgto_ket is the number of spherical GTOs on ket center,
  !>        equals to \f$2\var(angular_ket)+1\f$
  !> \param num_derv is the number of geometric derivatives
  !> \return contr_odist contains the overlap distribution
  subroutine contr_sgto_odist_recurr(coord_bra, angular_bra, num_prim_bra,        &
                                     exponent_bra, num_contr_bra, contr_coef_bra, &
                                     coord_ket, angular_ket, num_prim_ket,        &
                                     exponent_ket, num_contr_ket, contr_coef_ket, &
                                     order_geo_bra, order_geo_ket, num_points,    &
                                     grid_points, num_sgto_bra, num_sgto_ket,     &
                                     num_derv, contr_odist)
    use xkind
    implicit none
    real(REALK), intent(in) :: coord_bra(3)
    integer, intent(in) :: angular_bra
    integer, intent(in) :: num_prim_bra
    real(REALK), intent(in) :: exponent_bra(num_prim_bra)
    integer, intent(in) :: num_contr_bra
    real(REALK), intent(in) :: contr_coef_bra(num_contr_bra,num_prim_bra)
    real(REALK), intent(in) :: coord_ket(3)
    integer, intent(in) :: angular_ket
    integer, intent(in) :: num_prim_ket
    real(REALK), intent(in) :: exponent_ket(num_prim_ket)
    integer, intent(in) :: num_contr_ket
    real(REALK), intent(in) :: contr_coef_ket(num_contr_ket,num_prim_ket)
    integer, intent(in) :: order_geo_bra
    integer, intent(in) :: order_geo_ket
    integer, intent(in) :: num_points
    real(REALK), intent(in) :: grid_points(3,num_points)
    integer, intent(in) :: num_sgto_bra
    integer, intent(in) :: num_sgto_ket
    integer, intent(in) :: num_derv
    real(REALK), intent(out) :: contr_odist(num_sgto_bra,num_contr_bra, &
                                            num_sgto_ket,num_contr_ket, &
                                            num_points, num_derv)
!f2py intent(in) :: coord_bra
!f2py intent(in) :: angular_bra
!f2py intent(hide) :: num_prim_bra
!f2py intent(in) :: exponent_bra
!f2py intent(hide) :: num_contr_bra
!f2py intent(in) :: contr_coef_bra
!f2py depend(num_prim_bra) :: contr_coef_bra
!f2py intent(in) :: coord_ket
!f2py intent(in) :: angular_ket
!f2py intent(hide) :: num_prim_ket
!f2py intent(in) :: exponent_ket
!f2py intent(hide) :: num_contr_ket
!f2py intent(in) :: contr_coef_ket
!f2py depend(num_prim_ket) :: contr_coef_ket
!f2py intent(in) :: order_geo_bra
!f2py intent(in) :: order_geo_ket
!f2py intent(hide) :: num_points
!f2py intent(in) :: grid_points
!f2py intent(in) :: num_sgto_bra
!f2py intent(in) :: num_sgto_ket
!f2py intent(in) :: num_derv
!f2py intent(out) :: contr_odist
!f2py depend(num_sgto_bra) :: contr_odist
!f2py depend(num_contr_bra) :: contr_odist
!f2py depend(num_sgto_ket) :: contr_odist
!f2py depend(num_contr_ket) :: contr_odist
!f2py depend(num_points) :: contr_odist
!f2py depend(num_derv) :: contr_odist
    integer orders_hgto_bra(2)  !range of orders of Hermite Gaussians on bra center
    integer orders_hgto_ket(2)  !range of orders of Hermite Gaussians on ket center
    integer dim_hgto_bra        !dimension of Hermite Gaussians on bra center including geometric derivatives
    integer dim_hgto_ket        !dimension of Hermite Gaussians on ket center including geometric derivatives
    integer num_hgto_bra        !number of Hermite Gaussians on bra center
    integer num_hgto_ket        !number of Hermite Gaussians on ket center
    integer num_geo_bra         !number of geometric derivatives on bra center
    integer num_geo_ket         !number of geometric derivatives on ket center
    integer num_opt_geo         !number of operators and geometric derivatives
    real(REALK), allocatable :: hgto_odist(:,:,:,:,:)       !overlap distribution with HGTOs
    real(REALK), allocatable :: kgeo_hgto_odist(:,:,:,:)    !geometric derivatives (on ket center) of
                                                            !overlap distribution with HGTOs
    real(REALK), allocatable :: geom_hgto_odist(:,:,:,:,:)  !geometric derivatives of overlap distribution
                                                            !with HGTOs
    real(REALK), allocatable :: sgto_ket_odist(:,:,:,:,:)   !overlap distribution with contracted SGTOs
                                                            !on ket center
    real(REALK), allocatable :: scal_geo_bra(:)             !scale constants when recovering
                                                            !geometric derivatives on bra center
    real(REALK) scal_geo_ket                                !scale constant when recovering
                                                            !geometric derivatives on ket center
    integer iprim, jprim        !incremental recorders over primitives
    integer ierr                !error information
#if defined(XTIME)
    real(REALK) curr_time       !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! sets the orders of Hermite Gaussians on bra and ket centers
    orders_hgto_bra(1) = order_geo_bra+angular_bra
    orders_hgto_bra(2) = orders_hgto_bra(1)
    orders_hgto_ket(1) = order_geo_ket+angular_ket
    orders_hgto_ket(2) = orders_hgto_ket(1)
    ! sets the dimension of primitive Hermite Gaussians including geometric derivatives
    dim_hgto_bra = (orders_hgto_bra(1)+1)*(orders_hgto_bra(1)+2)/2
    dim_hgto_ket = (orders_hgto_ket(1)+1)*(orders_hgto_ket(1)+2)/2
    ! computes the number of primitive Hermite Gaussians
    num_hgto_bra = (angular_bra+1)*(angular_bra+2)/2
    num_hgto_ket = (angular_ket+1)*(angular_ket+2)/2
    ! computes the number of geometric derivatives
    num_geo_bra = (order_geo_bra+1)*(order_geo_bra+2)/2
    num_geo_ket = (order_geo_ket+1)*(order_geo_ket+2)/2
    num_opt_geo = num_points*num_geo_bra*num_geo_ket
    ! allocates the memory for overlap distribution between primitive HGTOs
    allocate(hgto_odist(dim_hgto_bra,dim_hgto_ket,num_points,1,1), stat=ierr)
    if (ierr/=0)                                       &
      call error_stop("contr_sgto_odist_recurr",       &
                      "failed to allocate hgto_odist", &
                      dim_hgto_bra*dim_hgto_ket*num_points)
    allocate(kgeo_hgto_odist(dim_hgto_bra,num_hgto_ket, &
                             num_points,num_geo_ket), stat=ierr)
    if (ierr/=0)                                            &
      call error_stop("contr_sgto_odist_recurr",            &
                      "failed to allocate kgeo_hgto_odist", &
                      dim_hgto_bra*num_hgto_ket*num_points*num_geo_ket)
    allocate(geom_hgto_odist(num_hgto_bra,num_hgto_ket,num_opt_geo, &
                             num_prim_bra,num_prim_ket), stat=ierr)
    if (ierr/=0)                                            &
      call error_stop("contr_sgto_odist_recurr",            &
                      "failed to allocate geom_hgto_odist", &
                      num_hgto_bra*num_hgto_ket*num_opt_geo &
                      *num_prim_bra*num_prim_ket)
#if defined(DEBUG)
    write(STDOUT,100) "size of overlap distribution with primitive HGTOs:", &
                      size(hgto_odist)
    write(STDOUT,100) "size of geometric derivatives of overlap distribution "// &
                      "with primitive HGTOs (ket):", size(kgeo_hgto_odist)
    write(STDOUT,100) "size of geometric derivatives of overlap distribution "// &
                      "with primitive HGTOs:", size(geom_hgto_odist)
#endif
    ! sets the scale constant when recovering geometric derivatives on bra center
    if (order_geo_bra>0) then
      allocate(scal_geo_bra(num_prim_bra), stat=ierr)
      if (ierr/=0)                                         &
        call error_stop("contr_sgto_odist_recurr",         &
                        "failed to allocate scal_geo_bra", &
                        num_prim_bra)
      do iprim = 1, num_prim_bra
        scal_geo_bra(iprim) = (exponent_bra(iprim)+exponent_bra(iprim))**order_geo_bra
      end do
      if (order_geo_ket>0) then
        ! computes the overlap distribution with primitive HGTOs
        do jprim = 1, num_prim_ket
          ! sets the scale constant when recovering geometric derivatives on ket center
          scal_geo_ket = (exponent_ket(jprim)+exponent_ket(jprim))**order_geo_ket
          do iprim = 1, num_prim_bra
            ! computes the overlap distribution with primitive HGTOs
            call prim_hgto_odist(orders_hgto_bra, coord_bra, exponent_bra(iprim), &
                                 orders_hgto_ket, coord_ket, exponent_ket(jprim), &
                                 num_points, grid_points, dim_hgto_bra,           &
                                 dim_hgto_ket, hgto_odist)
            ! scales the primitive Hermite integrals
            hgto_odist = scal_geo_ket*hgto_odist
            hgto_odist = scal_geo_bra(iprim)*hgto_odist
            ! recovers the geometric derivatives on ket center
            call scatter_single(dim_hgto_bra, dim_hgto_ket, num_points, 1, hgto_odist, &
                                angular_ket, order_geo_ket, num_hgto_ket, num_geo_ket, &
                                kgeo_hgto_odist)
            ! recovers the geometric derivatives on bra center
            call scatter_single(1, dim_hgto_bra, num_hgto_ket*num_points, num_geo_ket,     &
                                kgeo_hgto_odist, angular_bra, order_geo_bra, num_hgto_bra, &
                                num_geo_bra, geom_hgto_odist(:,:,:,iprim,jprim))
          end do
        end do
      else
        ! computes the overlap distribution with primitive HGTOs
        do jprim = 1, num_prim_ket
          do iprim = 1, num_prim_bra
            ! computes the overlap distribution with primitive HGTOs
            call prim_hgto_odist(orders_hgto_bra, coord_bra, exponent_bra(iprim), &
                                 orders_hgto_ket, coord_ket, exponent_ket(jprim), &
                                 num_points, grid_points, dim_hgto_bra,           &
                                 dim_hgto_ket, hgto_odist)
            ! scales the primitive Hermite integrals
            hgto_odist = scal_geo_bra(iprim)*hgto_odist
            ! recovers the geometric derivatives on ket center
            call scatter_single(dim_hgto_bra, dim_hgto_ket, num_points, 1, hgto_odist, &
                                angular_ket, order_geo_ket, num_hgto_ket, num_geo_ket, &
                                kgeo_hgto_odist)
            ! recovers the geometric derivatives on bra center
            call scatter_single(1, dim_hgto_bra, num_hgto_ket*num_points, num_geo_ket,     &
                                kgeo_hgto_odist, angular_bra, order_geo_bra, num_hgto_bra, &
                                num_geo_bra, geom_hgto_odist(:,:,:,iprim,jprim))
          end do
        end do
      end if
      deallocate(scal_geo_bra)
    else if (order_geo_ket>0) then
      ! computes the overlap distribution with primitive HGTOs
      do jprim = 1, num_prim_ket
        ! sets the scale constant when recovering geometric derivatives on ket center
        scal_geo_ket = (exponent_ket(jprim)+exponent_ket(jprim))**order_geo_ket
        do iprim = 1, num_prim_bra
          ! computes the overlap distribution with primitive HGTOs
          call prim_hgto_odist(orders_hgto_bra, coord_bra, exponent_bra(iprim), &
                               orders_hgto_ket, coord_ket, exponent_ket(jprim), &
                               num_points, grid_points, dim_hgto_bra,           &
                               dim_hgto_ket, hgto_odist)
          ! scales the primitive Hermite integrals
          hgto_odist = scal_geo_ket*hgto_odist
          ! recovers the geometric derivatives on ket center
          call scatter_single(dim_hgto_bra, dim_hgto_ket, num_points, 1, hgto_odist, &
                              angular_ket, order_geo_ket, num_hgto_ket, num_geo_ket, &
                              kgeo_hgto_odist)
          ! recovers the geometric derivatives on bra center
          call scatter_single(1, dim_hgto_bra, num_hgto_ket*num_points, num_geo_ket,     &
                              kgeo_hgto_odist, angular_bra, order_geo_bra, num_hgto_bra, &
                              num_geo_bra, geom_hgto_odist(:,:,:,iprim,jprim))
        end do
      end do
    else
      ! computes the overlap distribution with primitive HGTOs
      do jprim = 1, num_prim_ket
        do iprim = 1, num_prim_bra
          ! computes the overlap distribution with primitive HGTOs
          call prim_hgto_odist(orders_hgto_bra, coord_bra, exponent_bra(iprim), &
                               orders_hgto_ket, coord_ket, exponent_ket(jprim), &
                               num_points, grid_points, dim_hgto_bra,           &
                               dim_hgto_ket, hgto_odist)
          ! recovers the geometric derivatives on ket center
          call scatter_single(dim_hgto_bra, dim_hgto_ket, num_points, 1, hgto_odist, &
                              angular_ket, order_geo_ket, num_hgto_ket, num_geo_ket, &
                              kgeo_hgto_odist)
          ! recovers the geometric derivatives on bra center
          call scatter_single(1, dim_hgto_bra, num_hgto_ket*num_points, num_geo_ket,     &
                              kgeo_hgto_odist, angular_bra, order_geo_bra, num_hgto_bra, &
                              num_geo_bra, geom_hgto_odist(:,:,:,iprim,jprim))
        end do
      end do
    end if
    ! cleans
    deallocate(hgto_odist)
    deallocate(kgeo_hgto_odist)
    ! allocates the memory for geometric derivatives of overlap distribition
    ! with contracted HGTOs
    allocate(hgto_odist(num_hgto_bra,num_contr_bra, &
                        num_hgto_ket,num_contr_ket,num_opt_geo), stat=ierr)
    if (ierr/=0)                                       &
      call error_stop("contr_sgto_odist_recurr",       &
                      "failed to allocate hgto_odist", &
                      num_hgto_bra*num_contr_bra       &
                      *num_hgto_ket*num_contr_ket*num_opt_geo)
    ! allocates the memory for geometric derivatives of overlap distribution
    ! with contracted SGTOs on ket center
    allocate(sgto_ket_odist(num_hgto_bra,num_contr_bra, &
                            num_sgto_ket,num_contr_ket,num_opt_geo), stat=ierr)
    if (ierr/=0)                                           &
      call error_stop("contr_sgto_odist_recurr",           &
                      "failed to allocate sgto_ket_odist", &
                      num_hgto_bra*num_contr_bra           &
                      *num_sgto_ket*num_contr_ket*num_opt_geo)
#if defined(DEBUG)
    write(STDOUT,100) "size of geometric derivatives of overlap distribution "// &
                      "with contracted HGTOs:", size(hgto_odist)
    write(STDOUT,100) "size of overlap distribution with contracted SGTOs "// &
                      "on ket center:", size(sgto_ket_odist)
#endif
    ! constructs overlap distribution with contracted HGTOs
    ! in the order of hgto_odist(num_hgto_bra,num_contr_bra,num_hgto_ket,num_contr_ket, &
    !                            num_points,num_geo_bra,num_geo_ket)
    call const_contr_ints(num_contr_bra, num_prim_bra, contr_coef_bra, &
                          num_contr_ket, num_prim_ket, contr_coef_ket, &
                          num_hgto_bra, num_hgto_ket, num_opt_geo,     &
                          geom_hgto_odist, hgto_odist)
    deallocate(geom_hgto_odist)
    ! transforms the HGTOs to SGTOs on ket center
    call hgto_to_sgto_oop(angular_ket, num_hgto_bra*num_contr_bra, &
                          num_hgto_ket, num_contr_ket*num_opt_geo, &
                          hgto_odist, num_sgto_ket, sgto_ket_odist)
    deallocate(hgto_odist)
    ! transforms the HGTOs to SGTOs on bra center
    call hgto_to_sgto_oop(angular_bra, 1, num_hgto_bra,                         &
                          num_contr_bra*num_sgto_ket*num_contr_ket*num_opt_geo, &
                          sgto_ket_odist, num_sgto_bra, contr_odist)
    deallocate(sgto_ket_odist)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "contr_sgto_odist_recurr", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("contr_sgto_odist_recurr>> ",A,I8,1X,A)
#endif
  end subroutine contr_sgto_odist_recurr
