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
!!  This file recovers the HGTOs on ket center in nuclear attraction potential integrals.
!!
!!  2012-03-04, Bin Gao
!!  * rewrites to improve efficiency
!!
!!  2011-10-18, Bin Gao
!!  * first version

#include "stdout.h"

  !> \brief recovers the HGTOs on ket center in nuclear attraction potential integrals
  !> \author Bin Gao
  !> \date 2011-10-18
  !> \param orders_hgto_ket contains the minimum and maximum orders of HGTOs to return
  !> \param orders_geo_pot contains the minimum and maximum orders of geometric derivatives to return
  !> \param coord_ket contains the coordinates of bra center
  !> \param exponent_ket is the exponent of HGTOs of bra center
  !> \param coord_ket contains the coordinates of ket center
  !> \param exponent_ket is the exponent of HGTOs of ket center
  !> \param dim_geo_pot is the dimension of geometric derivatives
  !> \param geo_pot_pints contains the nuclear attraction integrals with zeroth order
  !>        HGTOs on ket center
  !> \param dim_hgto_ket is the dimension of HGTOs of ket center afterwards
  !> \param dim_geo_hket is the dimension of geometric derivatives afterwards
  !> \return hket_pints contains the integrals with required HGTOs on ket center
  subroutine nucpot_hket(orders_hgto_ket, orders_geo_pot, coord_bra, exponent_bra, &
                         coord_ket, exponent_ket, dim_geo_pot, geo_pot_pints,      &
                         dim_hgto_ket, dim_geo_hket, hket_pints)
    use xkind
    implicit none
    integer, intent(in) :: orders_hgto_ket(2)
    integer, intent(in) :: orders_geo_pot(2)
    real(REALK), intent(in) :: coord_bra(3)
    real(REALK), intent(in) :: exponent_bra
    real(REALK), intent(in) :: coord_ket(3)
    real(REALK), intent(in) :: exponent_ket
    integer, intent(in) :: dim_geo_pot
    real(REALK), intent(in) :: geo_pot_pints(dim_geo_pot)
    integer, intent(in) :: dim_hgto_ket
    integer, intent(in) :: dim_geo_hket
    real(REALK), intent(out) :: hket_pints(dim_hgto_ket,dim_geo_hket)
!f2py intent(in) :: orders_hgto_ket
!f2py intent(in) :: orders_geo_pot
!f2py intent(in) :: coord_bra
!f2py intent(in) :: exponent_bra
!f2py intent(in) :: coord_ket
!f2py intent(in) :: exponent_ket
!f2py intent(hide) :: dim_geo_pot
!f2py intent(in) :: geo_pot_pints
!f2py intent(in) :: dim_hgto_ket
!f2py intent(in) :: dim_geo_hket
!f2py intent(out) :: hket_pints
!f2py depend(dim_hgto_ket) :: hket_pints
!f2py depend(dim_geo_hket) :: hket_pints
    logical zero_hket          !if returning zeroth order HGTO on ket center
    integer max_low_geo        !maximum of lower order of geometric derivatives
    integer max_cur_hket       !maximum of current order of HGTOs
    real(REALK) half_neg_rp    !half of the negative reciprocal of total exponent \f$p_{ij}\f$
    real(REALK) cc_wrt_ket(3)  !relative coordinates of center-of-charge w.r.t. ket center
    real(REALK) bra_to_ket     !ratio of exponent on bra center to that on ket center
    integer order_geo          !incremental recorder over orders of geometric derivatives
    integer num_up_geo         !number of upper order geometric derivatives
    integer num_low_geo        !number of lower order geometric derivatives
    integer start_up_geo       !start address of upper order geometric derivatives
    integer end_up_geo         !end address of upper order geometric derivatives
    integer start_low_geo      !start address of lower order geometric derivatives
    integer end_low_geo        !end address of lower order geometric derivatives
    integer dim_low_hket       !dimension of lower order HGTOs on ket center
    integer dim_up_hket        !dimension of upper order HGTOs on ket center
    integer size_low_hket      !size of temporary integrals of lower order HGTOs on ket center
    integer size_up_hket       !size of temporary integrals of upper order HGTOs on ket center
    integer low_hket_int       !pointer to temporary integrals of lower order HGTOs on ket center
    integer up_hket_int        !pointer to temporary integrals of upper order HGTOs on ket center
    integer dim_tmp            !dimension of temporary integrals
    real(REALK), allocatable :: tmp_ints(:,:)
                               !temporary integrals
    integer offset_hket_geo    !offset of geometric derivatives in returned integrals
    integer ierr               !error information
#if defined(XTIME)
    real(REALK) curr_time      !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    select case(orders_hgto_ket(2))
    ! returns s-shell HGTO
    case(0)
      hket_pints(1,:) = geo_pot_pints
    ! maximum returned HGTOs are p-shell
    case(1)
      end_low_geo = dim_geo_pot                  !end address of lower order geometric derivatives
      num_low_geo = (orders_geo_pot(2)+2) &      !number of lower order geometric derivatives
                  * (orders_geo_pot(2)+3)/2
      start_low_geo = end_low_geo-num_low_geo+1  !start address of lower order derivatives
      ! allocates memory for temporary integrals
      allocate(tmp_ints(3*(num_low_geo-(orders_geo_pot(2)+2)),1), stat=ierr)
      if (ierr/=0)                                                    &
        call error_stop("nucpot_hket", "failed to allocate tmp_ints", &
                        3*(num_low_geo-orders_geo_pot(2)-2))
      ! if returing zero order HGTO
      zero_hket = orders_hgto_ket(1)==0
      ! sets the offset of geometric derivatives in returned integrals
      offset_hket_geo = dim_geo_hket
      ! calculates the relative coordinates of center-of-charge w.r.t. ket center,
      ! and the half of negative reciprocal of total exponent
      half_neg_rp = 1.0_REALK/(exponent_bra+exponent_ket)
      do ierr = 1, 3
        cc_wrt_ket(ierr) = half_neg_rp*(exponent_bra*coord_bra(ierr) &
                         + exponent_ket*coord_ket(ierr))-coord_ket(ierr)
      end do
      half_neg_rp = -0.5_REALK*half_neg_rp
      ! loops over orders of geometric derivatives
      do order_geo = orders_geo_pot(2), orders_geo_pot(1), -1
        ! updates the numbers of lower and upper order geometric derivatives
        num_up_geo = num_low_geo
        num_low_geo = num_low_geo-(order_geo+2)  !=(order_geo+1)*(order_geo+2)/2
        ! updates the start and end addresses of lower and upper order geometric derivatives
        end_up_geo = end_low_geo
        start_up_geo = start_low_geo
        end_low_geo = start_low_geo-1
        start_low_geo = end_low_geo-num_low_geo+1
#if defined(DEBUG)
        write(STDOUT,100) "GEO_HGTO/p/loop/upper/start/end:", &
                          order_geo+1, start_up_geo, end_up_geo
        write(STDOUT,100) "GEO-HGTO/p/loop/lower/start/end:", &
                          order_geo, start_low_geo, end_low_geo
#endif
        ! sets the size of temporary integrals
        size_up_hket = 3*num_low_geo
        ! gets the p-shell HGTO and order \var(orders_geo_pot(2)) geometric derivatives integrals
        call nucpot_hket_p(order_geo, cc_wrt_ket, half_neg_rp,        &
               num_up_geo, geo_pot_pints(start_up_geo:end_up_geo),    &
               num_low_geo, geo_pot_pints(start_low_geo:end_low_geo), &
               3, tmp_ints(1:size_up_hket,1))
        ! updates the offset of geometric derivatives
        offset_hket_geo = offset_hket_geo-num_low_geo
        ! assigns the returned integrals
        call nucpot_hket_assign(start_low_geo-1, dim_geo_pot, geo_pot_pints, &
                                3, num_low_geo, tmp_ints(1:size_up_hket,1),  &
                                zero_hket, offset_hket_geo, dim_hgto_ket,    &
                                dim_geo_hket, hket_pints)
      end do
      deallocate(tmp_ints)
    ! maximum returned HGTOs are d-shell
    case(2)
      ! (1) \var(max_low_geo)-1 order geometric derivatives
      max_low_geo = orders_geo_pot(2)+2               !maximum of lower order of geometric derivatives
      end_up_geo = dim_geo_pot                        !end address of upper order geometric derivatives
      num_up_geo = (max_low_geo+1)*(max_low_geo+2)/2  !number of upper order geometric derivatives
      end_low_geo = end_up_geo-num_up_geo             !end address of lower order geometric derivatives
      start_up_geo = end_low_geo+1                    !start address of upper order geometric derivatives
      max_low_geo = max_low_geo-1
      num_low_geo = num_up_geo-(max_low_geo+2)        !number of lower order geometric derivatives
      start_low_geo = end_low_geo-num_low_geo+1       !end address of \var(orders_geo_pot(2)) order derivatives
      ! allocates memory for temporary integrals
      allocate(tmp_ints(9*num_low_geo,2), stat=ierr)
      if (ierr/=0)                                                    &
        call error_stop("nucpot_hket", "failed to allocate tmp_ints", &
                        9*num_low_geo*2)
#if defined(DEBUG)
      write(STDOUT,100) "GEO-s-HGTO/d/upper/start/end:", &
                        max_low_geo+1, start_up_geo, end_up_geo
      write(STDOUT,100) "GEO-s-HGTO/d/upper/start/end:", &
                        max_low_geo, start_low_geo, end_low_geo
#endif
      ! calculates the relative coordinates of center-of-charge w.r.t. ket center,
      ! and the half of negative reciprocal of total exponent
      half_neg_rp = 1.0_REALK/(exponent_bra+exponent_ket)
      do ierr = 1, 3
        cc_wrt_ket(ierr) = half_neg_rp*(exponent_bra*coord_bra(ierr) &
                         + exponent_ket*coord_ket(ierr))-coord_ket(ierr)
      end do
      half_neg_rp = -0.5_REALK*half_neg_rp
      ! sets the dimension of lower order HGTOs and size of temporary integrals, p-shell
      dim_low_hket = 3
      size_low_hket = dim_low_hket*num_low_geo
      ! gets the temporary p-shell HGTO integrals
      call nucpot_hket_p(max_low_geo, cc_wrt_ket, half_neg_rp,      &
             num_up_geo, geo_pot_pints(start_up_geo:end_up_geo),    &
             num_low_geo, geo_pot_pints(start_low_geo:end_low_geo), &
             dim_low_hket, tmp_ints(1:size_low_hket,1))
      ! (2) \var(max_low_geo)-2 order geometric derivatives
      max_low_geo = max_low_geo-1
      num_up_geo = num_low_geo
      num_low_geo = num_low_geo-(max_low_geo+2)
      end_up_geo = end_low_geo
      start_up_geo = start_low_geo
      end_low_geo = start_low_geo-1
      start_low_geo = end_low_geo-num_low_geo+1
#if defined(DEBUG)
      write(STDOUT,100) "GEO-p-HGTO/d/upper/start/end:", &
                        max_low_geo+1, start_up_geo, end_up_geo
      write(STDOUT,100) "GEO-p-HGTO/d/lower/start/end:", &
                        max_low_geo, start_low_geo, end_low_geo
#endif
      ! sets the dimension of upper order HGTOs and size of temporary integrals, p- and d-shell
      dim_up_hket = 9
      size_up_hket = dim_up_hket*num_low_geo
      ! gets the temporary p-shell HGTO integrals
      call nucpot_hket_p(max_low_geo, cc_wrt_ket, half_neg_rp,      &
             num_up_geo, geo_pot_pints(start_up_geo:end_up_geo),    &
             num_low_geo, geo_pot_pints(start_low_geo:end_low_geo), &
             dim_up_hket, tmp_ints(1:size_up_hket,2))
      ! sets the ratio of exponent on bra center to that on ket center
      bra_to_ket = exponent_bra/exponent_ket
      ! recovers d-shell HGTOs
      call nucpot_hket_d(max_low_geo, cc_wrt_ket, half_neg_rp, bra_to_ket,      &
                         num_low_geo, geo_pot_pints(start_low_geo:end_low_geo), &
                         dim_low_hket, num_up_geo, tmp_ints(1:size_low_hket,1), &
                         dim_up_hket, tmp_ints(1:size_up_hket,2))
      ! if returing zero order HGTO
      zero_hket = orders_hgto_ket(1)==0
      ! sets the offset of geometric derivatives
      offset_hket_geo = dim_geo_hket-num_low_geo
      ! assigns the returned integrals
      call nucpot_hket_assign(start_low_geo-1, dim_geo_pot, geo_pot_pints, &
                              dim_up_hket, num_low_geo,                    &
                              tmp_ints(1:size_up_hket,2), zero_hket,       &
                              offset_hket_geo, dim_hgto_ket, dim_geo_hket, &
                              hket_pints)
      ! initializes the pointers of HGTOs on ket center
      low_hket_int = 1
      up_hket_int = 2
      ! updates the dimension of lower order HGTOs
      dim_low_hket = dim_up_hket
      ! loops over other returned orders of geometric derivatives
      do order_geo = orders_geo_pot(2)-1, orders_geo_pot(1), -1
        ! updates the numbers of lower and upper order geometric derivatives
        num_up_geo = num_low_geo
        num_low_geo = num_low_geo-(order_geo+2)  !=(order_geo+1)*(order_geo+2)/2
        ! updates the start and end addresses of lower and upper order geometric derivatives
        end_up_geo = end_low_geo
        start_up_geo = start_low_geo
        end_low_geo = start_low_geo-1
        start_low_geo = end_low_geo-num_low_geo+1
#if defined(DEBUG)
        write(STDOUT,100) "GEO-HGTO/d/loop/upper/start/end:", &
                          order_geo+1, start_up_geo, end_up_geo
        write(STDOUT,100) "GEO-HGTO/d/loop/lower/start/end:", &
                          order_geo, start_low_geo, end_low_geo
#endif
        ! updates the sizes of temporary integrals
        size_low_hket = size_up_hket
        size_up_hket = dim_up_hket*num_low_geo
        ! switches the pointers
        low_hket_int = 3-low_hket_int
        up_hket_int = 3-up_hket_int
        ! gets the temporary p-shell HGTO integrals
        call nucpot_hket_p(order_geo, cc_wrt_ket, half_neg_rp,        &
               num_up_geo, geo_pot_pints(start_up_geo:end_up_geo),    &
               num_low_geo, geo_pot_pints(start_low_geo:end_low_geo), &
               dim_up_hket, tmp_ints(1:size_up_hket,up_hket_int))
        ! gets the temporary d-shell HGTO integrals
        call nucpot_hket_d(order_geo, cc_wrt_ket, half_neg_rp,       &
                           bra_to_ket, num_low_geo,                  &
                           geo_pot_pints(start_low_geo:end_low_geo), &
                           dim_low_hket, num_up_geo,                 &
                           tmp_ints(1:size_low_hket,low_hket_int),   &
                           dim_up_hket, tmp_ints(1:size_up_hket,up_hket_int))
        ! updates the offset of geometric derivatives
        offset_hket_geo = offset_hket_geo-num_low_geo
        ! assigns the returned integrals
        call nucpot_hket_assign(start_low_geo-1, dim_geo_pot, geo_pot_pints, &
                                dim_up_hket, num_low_geo,                    &
                                tmp_ints(1:size_up_hket,up_hket_int),        &
                                zero_hket, offset_hket_geo, dim_hgto_ket,    &
                                dim_geo_hket, hket_pints)
      end do
      deallocate(tmp_ints)
    ! maximum returned HGTOs are other shells, at least f-shell
    case default
      ! (1) \var(max_low_geo)-1 order geometric derivatives
      max_low_geo = orders_geo_pot(2) &               !maximum of lower order of geometric derivatives
                  + orders_hgto_ket(2)
      ! allocates memory for temporary integrals
      call dim_nucpot_hket(max_low_geo, orders_geo_pot(2), dim_tmp)
      allocate(tmp_ints(dim_tmp,2), stat=ierr)
      if (ierr/=0)                                                    &
        call error_stop("nucpot_hket", "failed to allocate tmp_ints", &
                        dim_tmp*2)
      end_up_geo = dim_geo_pot                        !end address of upper order geometric derivatives
      num_up_geo = (max_low_geo+1)*(max_low_geo+2)/2  !number of upper order geometric derivatives
      end_low_geo = end_up_geo-num_up_geo             !end address of lower order geometric derivatives
      start_up_geo = end_low_geo+1                    !start address of upper order geometric derivatives
      max_low_geo = max_low_geo-1
      num_low_geo = num_up_geo-(max_low_geo+2)        !number of lower order geometric derivatives
      start_low_geo = end_low_geo-num_low_geo+1       !end address of \var(orders_geo_pot(2)) order derivatives
#if defined(DEBUG)
      write(STDOUT,100) "GEO-s-HGTO/upper/start/end:", &
                        max_low_geo+1, start_up_geo, end_up_geo
      write(STDOUT,100) "GEO-s-HGTO/upper/start/end:", &
                        max_low_geo, start_low_geo, end_low_geo
#endif
      ! calculates the relative coordinates of center-of-charge w.r.t. ket center,
      ! and the half of negative reciprocal of total exponent
      half_neg_rp = 1.0_REALK/(exponent_bra+exponent_ket)
      do ierr = 1, 3
        cc_wrt_ket(ierr) = half_neg_rp*(exponent_bra*coord_bra(ierr) &
                         + exponent_ket*coord_ket(ierr))-coord_ket(ierr)
      end do
      half_neg_rp = -0.5_REALK*half_neg_rp
      ! sets the dimension of lower order HGTOs and size of temporary integrals, p-shell
      dim_low_hket = 3
      size_low_hket = dim_low_hket*num_low_geo
      ! gets the temporary p-shell HGTO integrals
      call nucpot_hket_p(max_low_geo, cc_wrt_ket, half_neg_rp,      &
             num_up_geo, geo_pot_pints(start_up_geo:end_up_geo),    &
             num_low_geo, geo_pot_pints(start_low_geo:end_low_geo), &
             dim_low_hket, tmp_ints(1:size_low_hket,1))
      ! (2) \var(max_low_geo)-2 order geometric derivatives
      max_low_geo = max_low_geo-1
      num_up_geo = num_low_geo
      num_low_geo = num_low_geo-(max_low_geo+2)
      end_up_geo = end_low_geo
      start_up_geo = start_low_geo
      end_low_geo = start_low_geo-1
      start_low_geo = end_low_geo-num_low_geo+1
#if defined(DEBUG)
      write(STDOUT,100) "GEO-p-HGTO/upper/start/end:", &
                        max_low_geo+1, start_up_geo, end_up_geo
      write(STDOUT,100) "GEO-p-HGTO/lower/start/end:", &
                        max_low_geo, start_low_geo, end_low_geo
#endif
      ! sets the dimension of upper order HGTOs and size of temporary integrals, p- and d-shell
      dim_up_hket = 9
      size_up_hket = dim_up_hket*num_low_geo
      ! gets the temporary p-shell HGTO integrals
      call nucpot_hket_p(max_low_geo, cc_wrt_ket, half_neg_rp,      &
             num_up_geo, geo_pot_pints(start_up_geo:end_up_geo),    &
             num_low_geo, geo_pot_pints(start_low_geo:end_low_geo), &
             dim_up_hket, tmp_ints(1:size_up_hket,2))
      ! sets the ratio of exponent on bra center to that on ket center
      bra_to_ket = exponent_bra/exponent_ket
      ! recovers d-shell HGTOs
      call nucpot_hket_d(max_low_geo, cc_wrt_ket, half_neg_rp, bra_to_ket,      &
                         num_low_geo, geo_pot_pints(start_low_geo:end_low_geo), &
                         dim_low_hket, num_up_geo, tmp_ints(1:size_low_hket,1), &
                         dim_up_hket, tmp_ints(1:size_up_hket,2))
      ! initializes the maximum of current order of HGTOs on ket center
      max_cur_hket = 1
      ! initializes the pointers of HGTOs on ket center
      low_hket_int = 1
      up_hket_int = 2
      ! (3) loops over the orders of geometric derivatives till the maximum returned
      ! order, the maximum of current order of HGTOs \var(max_cur_hket) needs to
      ! update each iteration
       do order_geo = max_low_geo-1, orders_geo_pot(2), -1
        ! updates the numbers of lower and upper order geometric derivatives
        num_up_geo = num_low_geo
        num_low_geo = num_low_geo-(order_geo+2)  !=(order_geo+1)*(order_geo+2)/2
        ! updates the start and end addresses of lower and upper order geometric derivatives
        end_up_geo = end_low_geo
        start_up_geo = start_low_geo
        end_low_geo = start_low_geo-1
        start_low_geo = end_low_geo-num_low_geo+1
#if defined(DEBUG)
        write(STDOUT,100) "GEO-HGTO/loop/1/upper/start/end:", &
                          order_geo+1, start_up_geo, end_up_geo
        write(STDOUT,100) "GEO-HGTO/loop/1/lower/start/end:", &
                          order_geo, start_low_geo, end_low_geo
#endif
        ! updates the maximum of current order of HGTOs
        max_cur_hket = max_cur_hket+1
        ! updates the dimensions of HGTOs on ket center
        dim_low_hket = dim_up_hket
        dim_up_hket = dim_low_hket+(max_cur_hket+2)*(max_cur_hket+3)/2
        ! updates the sizes of temporary integrals
        size_low_hket = size_up_hket
        size_up_hket = dim_up_hket*num_low_geo
        ! switches the pointers
        low_hket_int = 3-low_hket_int
        up_hket_int = 3-up_hket_int
        ! gets the temporary p-shell HGTO integrals
        call nucpot_hket_p(order_geo, cc_wrt_ket, half_neg_rp,        &
               num_up_geo, geo_pot_pints(start_up_geo:end_up_geo),    &
               num_low_geo, geo_pot_pints(start_low_geo:end_low_geo), &
               dim_up_hket, tmp_ints(1:size_up_hket,up_hket_int))
        ! gets the temporary d-shell HGTO integrals
        call nucpot_hket_d(order_geo, cc_wrt_ket, half_neg_rp,       &
                           bra_to_ket, num_low_geo,                  &
                           geo_pot_pints(start_low_geo:end_low_geo), &
                           dim_low_hket, num_up_geo,                 &
                           tmp_ints(1:size_low_hket,low_hket_int),   &
                           dim_up_hket, tmp_ints(1:size_up_hket,up_hket_int))
        ! gets the temporary integrals for other HGTOs (f, g, ...)
        ! and \var(order_geo) order geometric derivatives
        call sub_nucpot_hket(order_geo, max_cur_hket, cc_wrt_ket, half_neg_rp,    &
                             bra_to_ket, dim_low_hket, num_up_geo,                &
                             tmp_ints(1:size_low_hket,low_hket_int), dim_up_hket, &
                             num_low_geo, tmp_ints(1:size_up_hket,up_hket_int))
      end do
      ! if returing zero order HGTO
      zero_hket = orders_hgto_ket(1)==0
      ! sets the offset of geometric derivatives
      offset_hket_geo = dim_geo_hket-num_low_geo
      ! assigns the returned integrals
      call nucpot_hket_assign(start_low_geo-1, dim_geo_pot, geo_pot_pints, &
                              dim_up_hket, num_low_geo,                    &
                              tmp_ints(1:size_up_hket,up_hket_int),        &
                              zero_hket, offset_hket_geo, dim_hgto_ket,    &
                              dim_geo_hket, hket_pints)
      ! updates the dimension of lower order HGTOs
      dim_low_hket = dim_up_hket
      ! loops over other returned orders of geometric derivatives, the maximum
      ! of current order of HGTOs \var(max_cur_hket) does not need to update
      do order_geo = orders_geo_pot(2)-1, orders_geo_pot(1), -1
        ! updates the numbers of lower and upper order geometric derivatives
        num_up_geo = num_low_geo
        num_low_geo = num_low_geo-(order_geo+2)  !=(order_geo+1)*(order_geo+2)/2
        ! updates the start and end addresses of lower and upper order geometric derivatives
        end_up_geo = end_low_geo
        start_up_geo = start_low_geo
        end_low_geo = start_low_geo-1
        start_low_geo = end_low_geo-num_low_geo+1
#if defined(DEBUG)
        write(STDOUT,100) "GEO-HGTO/loop/2/upper/start/end:", &
                          order_geo+1, start_up_geo, end_up_geo
        write(STDOUT,100) "GEO-HGTO/loop/2/lower/start/end:", &
                          order_geo, start_low_geo, end_low_geo
#endif
        ! updates the sizes of temporary integrals
        size_low_hket = size_up_hket
        size_up_hket = dim_up_hket*num_low_geo
        ! switches the pointers
        low_hket_int = 3-low_hket_int
        up_hket_int = 3-up_hket_int
        ! gets the temporary p-shell HGTO integrals
        call nucpot_hket_p(order_geo, cc_wrt_ket, half_neg_rp,        &
               num_up_geo, geo_pot_pints(start_up_geo:end_up_geo),    &
               num_low_geo, geo_pot_pints(start_low_geo:end_low_geo), &
               dim_up_hket, tmp_ints(1:size_up_hket,up_hket_int))
        ! gets the temporary d-shell HGTO integrals
        call nucpot_hket_d(order_geo, cc_wrt_ket, half_neg_rp,       &
                           bra_to_ket, num_low_geo,                  &
                           geo_pot_pints(start_low_geo:end_low_geo), &
                           dim_low_hket, num_up_geo,                 &
                           tmp_ints(1:size_low_hket,low_hket_int),   &
                           dim_up_hket, tmp_ints(1:size_up_hket,up_hket_int))
        ! gets the temporary integrals for other HGTOs (f, g, ...)
        ! and \var(order_geo) order geometric derivatives
        call sub_nucpot_hket(order_geo, max_cur_hket, cc_wrt_ket, half_neg_rp,    &
                             bra_to_ket, dim_low_hket, num_up_geo,                &
                             tmp_ints(1:size_low_hket,low_hket_int), dim_up_hket, &
                             num_low_geo, tmp_ints(1:size_up_hket,up_hket_int))
        ! updates the offset of geometric derivatives
        offset_hket_geo = offset_hket_geo-num_low_geo
        ! assigns the returned integrals
        call nucpot_hket_assign(start_low_geo-1, dim_geo_pot, geo_pot_pints, &
                                dim_up_hket, num_low_geo,                    &
                                tmp_ints(1:size_up_hket,up_hket_int),        &
                                zero_hket, offset_hket_geo, dim_hgto_ket,    &
                                dim_geo_hket, hket_pints)
      end do
      deallocate(tmp_ints)
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "nucpot_hket", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("nucpot_hket>> ",A,I6,2I8)
#endif
  end subroutine nucpot_hket

  !> \brief gets the maximum dimension of temporary integrals used in recurrence relations
  !> \author Bin Gao
  !> \date 2012-03-05
  !> \param max_order_geo is the maximum order of geometric derivatives
  !> \param min_order_geo is the minimum order of geometric derivatives
  !> \return dim_ints is the maximum dimension of temporary integrals
  subroutine dim_nucpot_hket(max_order_geo, min_order_geo, dim_ints)
    use xkind
    implicit none
    integer, intent(in) :: max_order_geo
    integer, intent(in) :: min_order_geo
    integer, intent(out) :: dim_ints
!f2py intent(in) :: max_order_geo
!f2py intent(in) :: min_order_geo
!f2py intent(out) :: dim_ints
    integer num_geo_pot    !number of upper order geometric derivatives
    integer max_hgto_ket   !maximum order of HGTOs on ket center
    integer dim_hgto_ket   !dimension of HGTOs on ket center
    integer order_geo      !incremental recorder over orders of geometric derivatives
    integer dim_tmp        !temporary result of dimension
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! initializes the return value
    dim_ints = 0
    ! initializes the number of geometric derivatives
    num_geo_pot = (max_order_geo+1)*(max_order_geo+2)/2
    ! initializes the maximum order and dimension of HGTOs on ket center
    max_hgto_ket = 0
    dim_hgto_ket = 0
    ! loops over the orders of geometric derivatives till the maximum order
    ! of returned geometric derivatives, the maximum order of HGTOs needs to
    ! update each iteration
    do order_geo = max_order_geo-1, min_order_geo, -1
      ! updates the numbers of geometric derivatives
      num_geo_pot = num_geo_pot-(order_geo+2)  !=(order_geo+1)*(order_geo+2)/2
      ! updates the maximum order of HGTOs
      max_hgto_ket = max_hgto_ket+1
      ! updates the dimension of HGTOs on ket center
      dim_hgto_ket = dim_hgto_ket+(max_hgto_ket+1)*(max_hgto_ket+2)/2
      ! updates the maximum dimension
      dim_tmp = dim_hgto_ket*num_geo_pot
      if (dim_tmp>dim_ints) dim_ints = dim_tmp
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "dim_nucpot_hket", STDOUT)
#endif
    return
  end subroutine dim_nucpot_hket

  !> \brief recovers p-shell HGTOs
  !> \author Bin Gao
  !> \date 2011-10-18
  !> \param cur_order_geo is current order of geometric derivatives
  !> \param cc_wrt_ket contains the relative coordinates of center-of-charge w.r.t. ket center
  !> \param half_neg_rp is the half of the negative reciprocal of total exponent \f$p_{ij}\f$
  !> \param num_up_geo is the number of upper order geometric derivatives
  !> \param up_geo_pints contains the integrals with s-shell HGTO and upper order
  !>        geometric derivatives
  !> \param num_low_geo is the number lower order geometric derivatives
  !> \param low_geo_pints contains the integrals with s-shell HGTO and lower order
  !>        geometric derivatives
  !> \param dim_up_hket is the dimension of upper order HGTOs
  !> \return up_hket_pints contains the integrals of p-shell HGTOs and lower order
  !>         geometric derivatives
  subroutine nucpot_hket_p(cur_order_geo, cc_wrt_ket, half_neg_rp, &
                           num_up_geo, up_geo_pints, num_low_geo,  &
                           low_geo_pints, dim_up_hket, up_hket_pints)
    use xkind
    implicit none
    integer, intent(in) :: cur_order_geo
    real(REALK), intent(in) :: cc_wrt_ket(3)
    real(REALK), intent(in) :: half_neg_rp
    integer, intent(in) :: num_up_geo
    real(REALK), intent(in) :: up_geo_pints(num_up_geo)
    integer, intent(in) :: num_low_geo
    real(REALK), intent(in) :: low_geo_pints(num_low_geo)
    integer, intent(in) :: dim_up_hket
    real(REALK), intent(inout) :: up_hket_pints(dim_up_hket,num_low_geo)
!f2py intent(in) :: cur_order_geo
!f2py intent(in) :: cc_wrt_ket
!f2py intent(in) :: half_neg_rp
!f2py intent(hide) :: num_up_geo
!f2py intent(in) :: up_geo_pints
!f2py intent(hide) :: num_low_geo
!f2py intent(in) :: low_geo_pints
!f2py intent(hide) :: dim_up_hket
!f2py intent(inout) :: up_hket_pints
!f2py depend(num_low_geo) :: up_hket_pints
    integer addr_up_geo    !address of upper order geometric derivatives
    integer addr_cur_geo   !address of current order geometric derivatives
    integer igeo, jgeo     !incremental recorders over geometric derivatives
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    addr_up_geo = 0
    addr_cur_geo = 0
    do igeo = cur_order_geo, 0, -1
      do jgeo = 0, igeo
        addr_up_geo = addr_up_geo+1
        addr_cur_geo = addr_cur_geo+1
        ! px
        up_hket_pints(1,addr_cur_geo)                 &
          = cc_wrt_ket(1)*low_geo_pints(addr_cur_geo) &
          + half_neg_rp*up_geo_pints(addr_up_geo)
        ! py
        up_hket_pints(2,addr_cur_geo)                 &
          = cc_wrt_ket(2)*low_geo_pints(addr_cur_geo) &
          + half_neg_rp*up_geo_pints(addr_up_geo+1)
        ! pz
        up_hket_pints(3,addr_cur_geo)                 &
          = cc_wrt_ket(3)*low_geo_pints(addr_cur_geo) &
          + half_neg_rp*up_geo_pints(addr_up_geo+igeo+2)
      end do
      addr_up_geo = addr_up_geo+1
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "nucpot_hket_p", STDOUT)
#endif
    return
  end subroutine nucpot_hket_p

  !> \brief recovers d-shell HGTOs
  !> \author Bin Gao
  !> \date 2011-10-18
  !> \param cur_order_geo is current order of geometric derivatives
  !> \param cc_wrt_ket contains the relative coordinates of center-of-charge w.r.t. ket center
  !> \param half_neg_rp is the half of the negative reciprocal of total exponent \f$p_{ij}\f$
  !> \param bra_to_ket is the ratio of exponent on bra center to that on ket center
  !> \param low_geo_pints contains the integrals with s-shell HGTO and lower order
  !>        geometric derivatives
  !> \param dim_low_hket is the dimension of lower order HGTOs
  !> \param num_up_geo is the number of upper order geometric derivatives
  !> \param low_hket_pints contains the integrals with p-shell HGTOs and upper order
  !>        geometric derivatives
  !> \param dim_up_hket is the dimension of upper order HGTOs
  !> \return up_hket_pints contains the integrals of d-shell HGTOs and lower order
  !>         geometric derivatives
  subroutine nucpot_hket_d(cur_order_geo, cc_wrt_ket, half_neg_rp,   &
                           bra_to_ket, num_low_geo, low_geo_pints,   &
                           dim_low_hket, num_up_geo, low_hket_pints, &
                           dim_up_hket, up_hket_pints)
    use xkind
    implicit none
    integer, intent(in) :: cur_order_geo
    real(REALK), intent(in) :: cc_wrt_ket(3)
    real(REALK), intent(in) :: half_neg_rp
    real(REALK), intent(in) :: bra_to_ket
    integer, intent(in) :: num_low_geo
    real(REALK), intent(in) :: low_geo_pints(num_low_geo)
    integer, intent(in) :: dim_low_hket
    integer, intent(in) :: num_up_geo
    real(REALK), intent(in) :: low_hket_pints(dim_low_hket,num_up_geo)
    integer, intent(in) :: dim_up_hket
    real(REALK), intent(inout) :: up_hket_pints(dim_up_hket,num_low_geo)
!f2py intent(in) :: cur_order_geo
!f2py intent(in) :: cc_wrt_ket
!f2py intent(in) :: half_neg_rp
!f2py intent(in) :: bra_to_ket
!f2py intent(hide) :: num_low_geo
!f2py intent(in) :: low_geo_pints
!f2py intent(hide) :: dim_low_hket
!f2py intent(hide) :: num_up_geo
!f2py intent(in) :: low_hket_pints
!f2py intent(hide) :: dim_up_hket
!f2py intent(inout) :: up_hket_pints
!f2py depend(num_low_geo) :: up_hket_pints
    integer addr_up_geo    !addresses of upper order geometric derivatives
    integer addr_up_geo_y
    integer addr_up_geo_z
    integer addr_cur_geo   !address of current order geometric derivatives
    integer igeo, jgeo     !incremental recorders over geometric derivatives
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    addr_up_geo = 0
    addr_cur_geo = 0
    do igeo = cur_order_geo, 0, -1
      do jgeo = 0, igeo
        addr_up_geo = addr_up_geo+1
        addr_cur_geo = addr_cur_geo+1
        ! dxx
        up_hket_pints(4,addr_cur_geo)                           &
          = cc_wrt_ket(1)*up_hket_pints(1,addr_cur_geo)         &
          + half_neg_rp*(bra_to_ket*low_geo_pints(addr_cur_geo) &
          + low_hket_pints(1,addr_up_geo))
        addr_up_geo_y = addr_up_geo+1
        ! dxy
        up_hket_pints(5,addr_cur_geo)                   &
          = cc_wrt_ket(2)*up_hket_pints(1,addr_cur_geo) &
          + half_neg_rp*low_hket_pints(1,addr_up_geo_y)
        ! dyy
        up_hket_pints(6,addr_cur_geo)                           &
          = cc_wrt_ket(2)*up_hket_pints(2,addr_cur_geo)         &
          + half_neg_rp*(bra_to_ket*low_geo_pints(addr_cur_geo) &
          + low_hket_pints(2,addr_up_geo_y))
        addr_up_geo_z = addr_up_geo+igeo+2
        ! dxz
        up_hket_pints(7,addr_cur_geo)                   &
          = cc_wrt_ket(3)*up_hket_pints(1,addr_cur_geo) &
          + half_neg_rp*low_hket_pints(1,addr_up_geo_z)
        ! dyz
        up_hket_pints(8,addr_cur_geo)                   &
          = cc_wrt_ket(3)*up_hket_pints(2,addr_cur_geo) &
          + half_neg_rp*low_hket_pints(2,addr_up_geo_z)
        ! dzz
        up_hket_pints(9,addr_cur_geo)                           &
          = cc_wrt_ket(3)*up_hket_pints(3,addr_cur_geo)         &
          + half_neg_rp*(bra_to_ket*low_geo_pints(addr_cur_geo) &
          + low_hket_pints(3,addr_up_geo_z))
      end do
      addr_up_geo = addr_up_geo+1
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "nucpot_hket_d", STDOUT)
#endif
    return
  end subroutine nucpot_hket_d

  !> \brief sub-recurrence relations by recovering upper order HGTOs on ket center
  !> \author Bin Gao
  !> \date 2011-10-18
  !> \param cur_order_geo is current order of geometric derivatives
  !> \param max_cur_hket is maximum of current order of HGTOs on ket center
  !> \param cc_wrt_ket contains the relative coordinates of center-of-charge w.r.t. ket center
  !> \param half_neg_rp is the half of the negative reciprocal of total exponent \f$p_{ij}\f$
  !> \param bra_to_ket is the ratio of exponent on bra center to that on ket center
  !> \param dim_low_hket is the dimension of lower order HGTOs
  !> \param num_up_geo is the number of upper order geometric derivatives
  !> \param low_hket_pints contains the integrals with lower order HGTOs and upper order
  !>        geometric derivatives
  !> \param dim_up_hket is the dimension of upper order HGTOs
  !> \param num_low_geo is the number lower order geometric derivatives
  !> \return up_hket_pints contains the integrals of upper order HGTOs and lower order
  !>         geometric derivatives
  subroutine sub_nucpot_hket(cur_order_geo, max_cur_hket, cc_wrt_ket, &
                             half_neg_rp, bra_to_ket, dim_low_hket,   &
                             num_up_geo, low_hket_pints, dim_up_hket, &
                             num_low_geo, up_hket_pints)
    use xkind
    implicit none
    integer, intent(in) :: cur_order_geo
    integer, intent(in) :: max_cur_hket
    real(REALK), intent(in) :: cc_wrt_ket(3)
    real(REALK), intent(in) :: half_neg_rp
    real(REALK), intent(in) :: bra_to_ket
    integer, intent(in) :: dim_low_hket
    integer, intent(in) :: num_up_geo
    real(REALK), intent(in) :: low_hket_pints(dim_low_hket,num_up_geo)
    integer, intent(in) :: dim_up_hket
    integer, intent(in) :: num_low_geo
    real(REALK), intent(inout) :: up_hket_pints(dim_up_hket,num_low_geo)
!f2py intent(in) :: cur_order_geo
!f2py intent(in) :: max_cur_hket
!f2py intent(in) :: cc_wrt_ket
!f2py intent(in) :: half_neg_rp
!f2py intent(in) :: bra_to_ket
!f2py intent(hide) :: dim_low_hket
!f2py intent(hide) :: num_up_geo
!f2py intent(in) :: low_hket_pints
!f2py intent(hide) :: dim_up_hket
!f2py intent(hide) :: num_low_geo
!f2py intent(inout) :: up_hket_pints
    integer addr_up_geo    !address of upper order geometric derivatives
    integer addr_cur_geo   !address of current order geometric derivatives
    integer igeo, jgeo     !incremental recorders over geometric derivatives
    integer order_hket     !order of HGTOs on ket center
    integer addr_up_hket   !address of upper order HGTOs
    integer addr_cur_hket  !address of current order HGTOs
    integer addr_low_hket  !address of lower order HGTOs
    integer iket, jket     !incremental recorder over HGTOs
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    addr_up_geo = 0
    addr_cur_geo = 0
    ! loops over xyz components of geometric derivatives
    do igeo = cur_order_geo, 0, -1
      do jgeo = 0, igeo
        addr_up_geo = addr_up_geo+1
        addr_cur_geo = addr_cur_geo+1
        addr_up_hket = 9   !base address of the f-shell HGTOs
        addr_cur_hket = 3  !base address of the d-shell HGTOs
        addr_low_hket = 0  !base address of the p-shell HGTOs
        ! loops over current order of HGTOs, starting from d-shell
        do order_hket = 2, max_cur_hket
          addr_up_hket = addr_up_hket+1
          addr_cur_hket = addr_cur_hket+1
#if defined(DEBUG)
          write(STDOUT,100) "orders:", order_hket, cur_order_geo
          write(STDOUT,100) "x-direction"
          write(STDOUT,110) addr_up_hket, addr_cur_geo
          write(STDOUT,110) addr_cur_hket, addr_cur_geo
          write(STDOUT,110) addr_low_hket+1, addr_cur_geo
          write(STDOUT,110) addr_cur_hket, addr_up_geo
          write(STDOUT,100) "------------------------------------"
#endif
          ! x-direction
          up_hket_pints(addr_up_hket,addr_cur_geo)                    &
            = cc_wrt_ket(1)*up_hket_pints(addr_cur_hket,addr_cur_geo) &
            + half_neg_rp*(real(order_hket,REALK)*bra_to_ket          &
            * up_hket_pints(addr_low_hket+1,addr_cur_geo)             &
            + low_hket_pints(addr_cur_hket,addr_up_geo))
          ! y-direction
          addr_up_hket = addr_up_hket+1
#if defined(DEBUG)
          write(STDOUT,100) "y-direction"
          write(STDOUT,110) addr_up_hket, addr_cur_geo
          write(STDOUT,110) addr_cur_hket, addr_cur_geo
          write(STDOUT,110) addr_cur_hket, addr_up_geo+1
          write(STDOUT,100) "------------------------------------"
#endif
          up_hket_pints(addr_up_hket,addr_cur_geo)                    &
            = cc_wrt_ket(2)*up_hket_pints(addr_cur_hket,addr_cur_geo) &
            + half_neg_rp*low_hket_pints(addr_cur_hket,addr_up_geo+1)
          do iket = 1, order_hket
            addr_up_hket = addr_up_hket+1
#if defined(DEBUG)
            write(STDOUT,110) addr_up_hket, addr_cur_geo
            write(STDOUT,110) addr_cur_hket+iket, addr_cur_geo
            write(STDOUT,110) addr_low_hket+iket, addr_cur_geo
            write(STDOUT,110) addr_cur_hket+iket, addr_up_geo+1
            write(STDOUT,100) "------------------------------------"
#endif
            up_hket_pints(addr_up_hket,addr_cur_geo)                         &
              = cc_wrt_ket(2)*up_hket_pints(addr_cur_hket+iket,addr_cur_geo) &
              + half_neg_rp*(real(iket,REALK)*bra_to_ket                     &
              * up_hket_pints(addr_low_hket+iket,addr_cur_geo)               &
              + low_hket_pints(addr_cur_hket+iket,addr_up_geo+1))
          end do
          ! z-direction
#if defined(DEBUG)
          write(STDOUT,100) "z-direction"
#endif
          addr_cur_hket = addr_cur_hket-1
          do jket = 0, order_hket
            addr_up_hket = addr_up_hket+1
            addr_cur_hket = addr_cur_hket+1
#if defined(DEBUG)
            write(STDOUT,110) addr_up_hket, addr_cur_geo
            write(STDOUT,110) addr_cur_hket, addr_cur_geo
            write(STDOUT,110) addr_cur_hket, addr_up_geo+igeo+2
            write(STDOUT,100) "------------------------------------"
#endif
            up_hket_pints(addr_up_hket,addr_cur_geo)                    &
              = cc_wrt_ket(3)*up_hket_pints(addr_cur_hket,addr_cur_geo) &
              + half_neg_rp*low_hket_pints(addr_cur_hket,addr_up_geo+igeo+2)
          end do
          do iket = 1, order_hket
            do jket = 0, order_hket-iket
              addr_up_hket = addr_up_hket+1
              addr_cur_hket = addr_cur_hket+1
              addr_low_hket = addr_low_hket+1
#if defined(DEBUG)
              write(STDOUT,110) addr_up_hket, addr_cur_geo
              write(STDOUT,110) addr_cur_hket, addr_cur_geo
              write(STDOUT,110) addr_low_hket, addr_cur_geo
              write(STDOUT,110) addr_cur_hket, addr_up_geo+igeo+2
              write(STDOUT,100) "------------------------------------"
#endif
              up_hket_pints(addr_up_hket,addr_cur_geo)                    &
                = cc_wrt_ket(3)*up_hket_pints(addr_cur_hket,addr_cur_geo) &
                + half_neg_rp*(real(iket,REALK)*bra_to_ket                &
                * up_hket_pints(addr_low_hket,addr_cur_geo)               &
                + low_hket_pints(addr_cur_hket,addr_up_geo+igeo+2))
            end do
          end do
        end do
      end do
      addr_up_geo = addr_up_geo+1
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "sub_nucpot_hket", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("sub_nucpot_hket>> ",A,2I6)
110 format("sub_nucpot_hket>> ","HGTO",I8,4X,"GEO",I8)
#endif
  end subroutine sub_nucpot_hket

  !> \brief assigns the integrals with required HGTOs on ket center
  !> \author Bin Gao
  !> \date 2012-03-04
  !> \param offset_geo_pot is the offset of geometric derivatives on nuclear potential origin
  !> \param dim_geo_pot is the dimension of geometric derivatives
  !> \param geo_pot_pints contains the nuclear attraction integrals with zeroth order
  !>        HGTOs on ket center
  !> \param dim_up_hket is the dimension of upper order HGTOs in temporary integrals
  !> \param num_low_geo is the number of lower order geometric derivatives in temporary integrals
  !> \param recur_pints contains the temporary integrals from recurrence relations
  !> \param zero_hket indicates if zeroth order HGTO returned
  !> \param offset_hket_geo is the offset of geometric derivatives on nuclear potential origin
  !>        in returned integrals
  !> \param dim_hgto_ket is the dimension of HGTOs of ket center afterwards
  !> \param dim_geo_hket is the dimension of geometric derivatives afterwards
  !> \return hket_pints contains the integrals with required HGTOs on ket center
  subroutine nucpot_hket_assign(offset_geo_pot, dim_geo_pot, geo_pot_pints, &
                                dim_up_hket, num_low_geo, recur_pints,      &
                                zero_hket, offset_hket_geo, dim_hgto_ket,   &
                                dim_geo_hket, hket_pints)
    use xkind
    implicit none
    integer, intent(in) :: offset_geo_pot
    integer, intent(in) :: dim_geo_pot
    real(REALK), intent(in) :: geo_pot_pints(dim_geo_pot)
    integer, intent(in) :: dim_up_hket
    integer, intent(in) :: num_low_geo
    real(REALK), intent(in) :: recur_pints(dim_up_hket,num_low_geo)
    logical, intent(in) :: zero_hket
    integer, intent(in) :: offset_hket_geo
    integer, intent(in) :: dim_hgto_ket
    integer, intent(in) :: dim_geo_hket
    real(REALK), intent(inout) :: hket_pints(dim_hgto_ket,dim_geo_hket)
!f2py intent(in) :: offset_geo_pot
!f2py intent(hide) :: dim_geo_pot
!f2py intent(in) :: geo_pot_pints
!f2py intent(hide) :: dim_up_hket
!f2py intent(hide) :: num_low_geo
!f2py intent(in) :: recur_pints
!f2py intent(in) :: zero_hket
!f2py intent(in) :: offset_hket_geo
!f2py intent(hide) :: dim_hgto_ket
!f2py intent(hide) :: dim_geo_hket
!f2py intent(inout) :: hket_pints
    integer start_up_hket  !start address of upper order HGTOs in temporary integrals
    integer addr_geo_pot   !address of geometric derivatives on nuclear potential origin
    integer addr_hket_geo  !address of geometric derivatives on nuclear potential origin in returned integrals
    integer igeo           !incremental recorder over geometric derivatives
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! s-shell HGTO returned
    if (zero_hket) then
      start_up_hket = dim_up_hket-dim_hgto_ket+2
      do igeo = 1, num_low_geo
        addr_geo_pot = offset_geo_pot+igeo
        addr_hket_geo = offset_hket_geo+igeo
        hket_pints(1,addr_hket_geo) = geo_pot_pints(addr_geo_pot)
        hket_pints(2:dim_hgto_ket,addr_hket_geo) &
          = recur_pints(start_up_hket:dim_up_hket,igeo)
      end do
    else
      start_up_hket = dim_up_hket-dim_hgto_ket+1
      do igeo = 1, num_low_geo
        addr_hket_geo = offset_hket_geo+igeo
        hket_pints(:,addr_hket_geo) &
          = recur_pints(start_up_hket:dim_up_hket,igeo)
      end do
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "nucpot_hket_assign", STDOUT)
#endif
    return
  end subroutine nucpot_hket_assign
