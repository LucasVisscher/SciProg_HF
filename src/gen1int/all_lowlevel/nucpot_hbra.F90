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
!!  This file recovers the HGTOs on bra center in nuclear attraction potential integrals.
!!
!!  2012-03-04, Bin Gao
!!  * rewrites to improve efficiency
!!
!!  2011-10-18, Bin Gao
!!  * first version

#include "stdout.h"

  !> \brief recovers the HGTOs on bra center in nuclear attraction potential integrals
  !> \author Bin Gao
  !> \date 2011-10-18
  !> \param orders_hgto_bra contains the minimum and maximum orders of HGTOs on bra center to return
  !> \param orders_hgto_ket contains the minimum and maximum orders of HGTOs on ket center to return
  !> \param orders_geo_pot contains the minimum and maximum orders of geometric derivatives to return
  !> \param coord_ket contains the coordinates of bra center
  !> \param exponent_ket is the exponent of HGTOs of bra center
  !> \param coord_ket contains the coordinates of bra center
  !> \param exponent_ket is the exponent of HGTOs of bra center
  !> \param dim_hket is the dimension of HGTOs on ket center
  !> \param dim_geo_hket is the dimension of geometric derivatives on nuclear potential origin
  !> \param hket_pints contains the nuclear attraction integrals with zeroth order
  !>        HGTOs on bra center
  !> \param dim_hgto_bra is the dimension of HGTOs of bra center afterwards
  !> \param dim_hgto_ket is the dimension of HGTOs of ket center afterwards
  !> \param dim_geo_hbra is the dimension of geometric derivatives afterwards
  !> \return hbra_pints contains the integrals with required HGTOs on bra center
  subroutine nucpot_hbra(orders_hgto_bra, orders_hgto_ket, orders_geo_pot, &
                         coord_bra, exponent_bra, coord_ket, exponent_ket, &
                         dim_hket, dim_geo_hket, hket_pints, dim_hgto_bra, &
                         dim_hgto_ket, dim_geo_hbra, hbra_pints)
    use xkind
    implicit none
    integer, intent(in) :: orders_hgto_bra(2)
    integer, intent(in) :: orders_hgto_ket(2)
    integer, intent(in) :: orders_geo_pot(2)
    real(REALK), intent(in) :: coord_bra(3)
    real(REALK), intent(in) :: exponent_bra
    real(REALK), intent(in) :: coord_ket(3)
    real(REALK), intent(in) :: exponent_ket
    integer, intent(in) :: dim_hket
    integer, intent(in) :: dim_geo_hket
    real(REALK), intent(in) :: hket_pints(dim_hket,dim_geo_hket)
    integer, intent(in) :: dim_hgto_bra
    integer, intent(in) :: dim_hgto_ket
    integer, intent(in) :: dim_geo_hbra
    real(REALK), intent(out) :: hbra_pints(dim_hgto_bra,dim_hgto_ket,dim_geo_hbra)
!f2py intent(in) :: orders_hgto_bra
!f2py intent(in) :: orders_hgto_ket
!f2py intent(in) :: orders_geo_pot
!f2py intent(in) :: coord_bra
!f2py intent(in) :: exponent_bra
!f2py intent(in) :: coord_ket
!f2py intent(in) :: exponent_ket
!f2py intent(hide) :: dim_hket
!f2py intent(hide) :: dim_geo_hket
!f2py intent(in) :: hket_pints
!f2py intent(in) :: dim_hgto_bra
!f2py intent(in) :: dim_hgto_ket
!f2py intent(in) :: dim_geo_hbra
!f2py intent(out) :: hbra_pints
!f2py depend(dim_hgto_bra) :: hbra_pints
!f2py depend(dim_hgto_ket) :: hbra_pints
!f2py depend(dim_geo_hbra) :: hbra_pints
    real(REALK) half_neg_rp    !half of the negative reciprocal of total exponent \f$p_{ij}\f$
    real(REALK) ket_to_bra     !ratio of exponent on ket center to that on bra center
    real(REALK) cc_wrt_bra(3)  !relative coordinates of center-of-charge w.r.t. bra center
    integer max_low_geo        !maximum of lower order of geometric derivatives
    integer num_up_geo         !number of higher order geometric derivatives
    integer num_low_geo        !number of lower order geometric derivatives
    integer start_up_geo       !start address of upper order geometric derivatives
    integer end_up_geo         !end address of upper order geometric derivatives
    integer start_low_geo      !start address of lower order geometric derivatives
    integer end_low_geo        !end address of lower order geometric derivatives
    integer min_order_hket     !minimum of current order of HGTOs on ket center
    integer max_order_hbra     !maximum of current order of HGTOs on bra center
    integer dim_cur_hket       !dimension of current order HGTOs on ket cetner of temporary integrals
    integer dim_hket_zero      !dimension of HGTOs on ket center with zeroth order HGTO on bra center
    integer dim_hket_low       !dimension of HGTOs on ket center with lower order HGTOs on bra center
    integer dim_low_hbra       !dimension of lower order HGTOs on bra center of temporary integrals
    integer dim_up_hbra        !dimension of upper order HGTOs on bra center of temporary integrals
    integer order_geo          !incremental recorder over orders of geometric derivatives
    integer max_dim_hket       !maximum dimension of HGTOs on ket center
    integer size_low_hbra      !size of temporary integrals of lower order HGTOs on bra center
    integer size_up_hbra       !size of temporary integrals of upper order HGTOs on bra center
    integer size_bra_ket       !size of HGTOs on bra and ket centers
    integer low_hbra_int       !pointer to temporary integrals of lower order HGTOs on bra center
    integer up_hbra_int        !pointer to temporary integrals of upper order HGTOs on bra center
    integer dim_tmp            !dimension of temporary integrals
    real(REALK), allocatable :: tmp_ints(:,:)
                               !temporary integrals
    integer offset_hbra_geo    !offset of geometric derivatives in returned integrals
    integer offset_hgto_ket    !offset of HGTOs on ket center in integrals with zeroth
                               !order HGTO on bra center
    logical zero_hbra          !if returning zeroth order HGTO on bra center
    integer offset_hket_low    !offset of HGTOs on ket center in temporary integrals with
                               !lower order HGTOs on bra center
    integer start_hket_up      !start addresses of HGTOs on ket center in integrals with zeroth order
    integer start_hket_low     !HGTO on bra center, and upper/lower order geometric derivatives
    integer ierr               !error information
#if defined(XTIME)
    real(REALK) curr_time      !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    select case(orders_hgto_bra(2))
    ! returns s-shell HGTO
    case(0)
      hbra_pints(1,:,:) = hket_pints
    ! maximum returned HGTOs are other shells, at least p-shell
    case default
      ! calculates the relative coordinates of center-of-charge w.r.t. bra center,
      ! and the half of negative reciprocal of total exponent
      half_neg_rp = 1.0_REALK/(exponent_bra+exponent_ket)
      do ierr = 1, 3
        cc_wrt_bra(ierr) = half_neg_rp*(exponent_bra*coord_bra(ierr) &
                         + exponent_ket*coord_ket(ierr))-coord_bra(ierr)
      end do
      half_neg_rp = -0.5_REALK*half_neg_rp
      ! sets the ratio of exponent on ket center to that on bra center
      ket_to_bra = exponent_ket/exponent_bra
      ! maximum dimension of HGTOs on ket center
      max_dim_hket = (orders_hgto_ket(2)+1)*(orders_hgto_ket(2)+2) &
                   * (orders_hgto_ket(2)+3)/6
      ! sets the maximum of lower order of geometric derivatives
      max_low_geo = orders_geo_pot(2)+orders_hgto_bra(2)
      ! allocates memory for temporary integrals
      call dim_nucpot_hbra(orders_hgto_ket(1), max_low_geo, orders_geo_pot(2), &
                           max_dim_hket, dim_tmp)
      allocate(tmp_ints(dim_tmp,2), stat=ierr)
      if (ierr/=0)                                                    &
        call error_stop("nucpot_hbra", "failed to allocate tmp_ints", &
                        dim_tmp*2)
      ! \var(max_low_geo)-1 order geometric derivatives
      end_up_geo = dim_geo_hket                       !end address of upper order geometric derivatives
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
      ! sets the minimum of current order of HGTOs on ket center
      min_order_hket = orders_hgto_ket(1)
      ! sets the dimensions of HGTOs on ket center
      if (min_order_hket==0) then
        dim_cur_hket = max_dim_hket
        dim_hket_zero = dim_cur_hket
      else
        ierr = min_order_hket*(min_order_hket+1)/2
        dim_cur_hket = max_dim_hket-ierr*(min_order_hket+2)/3
        dim_hket_zero = dim_cur_hket+ierr
      end if
      ! initializes the maximum of current order of HGTOs on bra center
      max_order_hbra = 0
      ! sets the dimensions and size of temporary integrals with lower order
      ! HGTOs on bra center (not used in recurrence relations)
      dim_low_hbra = 1
      dim_hket_low = dim_cur_hket
      size_low_hbra = dim_hket_low*num_up_geo
      ! sets the dimension and size of temporary integrals with upper order HGTOs on bra center
      dim_up_hbra = 3
      size_bra_ket = dim_up_hbra*dim_cur_hket
      size_up_hbra = size_bra_ket*num_low_geo
      ! sets the start addresses of HGTOs on ket center in integrals with zeroth
      ! order HGTO on bra center, and upper/lower order geometric derivatives
      start_hket_up = dim_hket-dim_cur_hket+1
      start_hket_low = dim_hket-dim_hket_zero+1
      ! gets the temporary integrals for order of HGTOs on bra center up to 
      ! \var(max_order_hbra)+1 and order of geometric derivatives as \var(max_low_geo)
      call sub_nucpot_hbra(max_low_geo, orders_hgto_ket, min_order_hket,    &
             max_order_hbra, cc_wrt_bra, half_neg_rp,                       &
             ket_to_bra, dim_cur_hket, num_up_geo,                          &
             hket_pints(start_hket_up:dim_hket,start_up_geo:end_up_geo),    &
             dim_hket_zero, num_low_geo,                                    &
             hket_pints(start_hket_low:dim_hket,start_low_geo:end_low_geo), &
             0, dim_low_hbra, dim_hket_low, tmp_ints(1:size_low_hbra,1),    &
             dim_up_hbra, tmp_ints(1:size_up_hbra,2))
      ! initializes the pointers of HGTOs on bra center
      low_hbra_int = 1
      up_hbra_int = 2
      ! loops over the orders of geometric derivatives not returned, the maximum
      ! of current order of HGTOs \var(max_order_hbra) needs to update in the cycle
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
        ! updates maximum of current order of HGTOs on bra center
        max_order_hbra = max_order_hbra+1
        ! updates the dimensions of HGTOs on bra center
        dim_low_hbra = dim_up_hbra
        dim_up_hbra = dim_low_hbra+(max_order_hbra+2)*(max_order_hbra+3)/2
        ! sets minimum of current order of HGTOs on ket center
        min_order_hket = max(orders_hgto_ket(1)-max_order_hbra,0)
        ! sets the dimensions of HGTOs on ket center
        dim_hket_low = dim_cur_hket
        if (min_order_hket==0) then
          dim_cur_hket = max_dim_hket
          dim_hket_zero = dim_cur_hket
        else
          ierr = min_order_hket*(min_order_hket+1)/2
          dim_cur_hket = max_dim_hket-ierr*(min_order_hket+2)/3
          dim_hket_zero = dim_cur_hket+ierr
        end if
        ! updates the sizes of temporary integrals
        size_low_hbra = size_up_hbra
        size_bra_ket = dim_up_hbra*dim_cur_hket
        size_up_hbra = size_bra_ket*num_low_geo
        ! switches the pointers
        low_hbra_int = 3-low_hbra_int
        up_hbra_int = 3-up_hbra_int
        ! sets the start addresses of HGTOs on ket center in integrals with zeroth
        ! order HGTO on bra center, and upper/lower order geometric derivatives
        start_hket_up = dim_hket-dim_cur_hket+1
        start_hket_low = dim_hket-dim_hket_zero+1
        ! gets the temporary integrals for order of HGTOs on bra center up to 
        ! \var(max_order_hbra)+1 and order of geometric derivatives as \var(order_geo)
        call sub_nucpot_hbra(order_geo, orders_hgto_ket, min_order_hket,      &
               max_order_hbra, cc_wrt_bra, half_neg_rp,                       &
               ket_to_bra, dim_cur_hket, num_up_geo,                          &
               hket_pints(start_hket_up:dim_hket,start_up_geo:end_up_geo),    &
               dim_hket_zero, num_low_geo,                                    &
               hket_pints(start_hket_low:dim_hket,start_low_geo:end_low_geo), &
               0, dim_low_hbra, dim_hket_low,                                 &
               tmp_ints(1:size_low_hbra,low_hbra_int),                        &
               dim_up_hbra, tmp_ints(1:size_up_hbra,up_hbra_int))
      end do
      ! if returing zero order HGTO
      zero_hbra = orders_hgto_bra(1)==0
      ! sets the offsets of geometric derivatives and HGTOs
      offset_hbra_geo = dim_geo_hbra-num_low_geo
      offset_hgto_ket = dim_hket-dim_hgto_ket
      ! assigns the returned integrals
      call nucpot_hbra_assign(offset_hgto_ket, start_low_geo-1, dim_hket, &
                              dim_geo_hket, hket_pints, dim_up_hbra,      &
                              dim_cur_hket, num_low_geo,                  &
                              tmp_ints(1:size_up_hbra,up_hbra_int),       &
                              zero_hbra, offset_hbra_geo, dim_hgto_bra,   &
                              dim_hgto_ket, dim_geo_hbra, hbra_pints)
      ! sets offset of HGTOs on ket center in temporary integrals with lower order HGTOs on bra center
      offset_hket_low = dim_cur_hket-dim_hket_low
      ! loops over other returned orders of geometric derivatives, maximum of current
      ! order of HGTOs on bra center \var(max_order_hbra) does not need to update
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
        size_low_hbra = size_up_hbra
        size_up_hbra = size_bra_ket*num_low_geo
        ! switches the pointers
        low_hbra_int = 3-low_hbra_int
        up_hbra_int = 3-up_hbra_int
        ! gets the temporary integrals for order of HGTOs on bra center up to 
        ! \var(max_order_hbra)+1 and order of geometric derivatives as \var(order_geo)
        call sub_nucpot_hbra(order_geo, orders_hgto_ket, min_order_hket,      &
               max_order_hbra, cc_wrt_bra, half_neg_rp,                       &
               ket_to_bra, dim_cur_hket, num_up_geo,                          &
               hket_pints(start_hket_up:dim_hket,start_up_geo:end_up_geo),    &
               dim_hket_zero, num_low_geo,                                    &
               hket_pints(start_hket_low:dim_hket,start_low_geo:end_low_geo), &
               offset_hket_low, dim_up_hbra, dim_cur_hket,                    &
               tmp_ints(1:size_low_hbra,low_hbra_int),                        &
               dim_up_hbra, tmp_ints(1:size_up_hbra,up_hbra_int))
        ! updates the offset of geometric derivatives
        offset_hbra_geo = offset_hbra_geo-num_low_geo
        ! assigns the returned integrals
        call nucpot_hbra_assign(offset_hgto_ket, start_low_geo-1, dim_hket, &
                                dim_geo_hket, hket_pints, dim_up_hbra,      &
                                dim_cur_hket, num_low_geo,                  &
                                tmp_ints(1:size_up_hbra,up_hbra_int),       &
                                zero_hbra, offset_hbra_geo, dim_hgto_bra,   &
                                dim_hgto_ket, dim_geo_hbra, hbra_pints)
      end do
      deallocate(tmp_ints)
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "nucpot_hbra", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("nucpot_hbra>> ",A,I6,2I8)
#endif
  end subroutine nucpot_hbra

  !> \brief gets the maximum dimension of temporary integrals used in recurrence relations
  !> \author Bin Gao
  !> \date 2012-03-05
  !> \param min_order_hket is the minimum order of HGTOs on ket center
  !> \param max_order_geo is the maximum order of geometric derivatives
  !> \param min_order_geo is the minimum order of geometric derivatives
  !> \param max_dim_hket is the maximum dimension of HGTOs on ket center
  !> \return dim_ints is the maximum dimension of temporary integrals
  subroutine dim_nucpot_hbra(min_order_hket, max_order_geo, min_order_geo, &
                             max_dim_hket, dim_ints)
    use xkind
    implicit none
    integer, intent(in) :: min_order_hket
    integer, intent(in) :: max_order_geo
    integer, intent(in) :: min_order_geo
    integer, intent(in) :: max_dim_hket
    integer, intent(out) :: dim_ints
!f2py intent(in) :: min_order_hket
!f2py intent(in) :: max_order_geo
!f2py intent(in) :: min_order_geo
!f2py intent(in) :: max_dim_hket
!f2py intent(out) :: dim_ints
    integer num_geo_pot    !number of upper order geometric derivatives
    integer max_hgto_bra   !maximum order of HGTOs on bra center
    integer dim_hgto_bra   !dimension of HGTOs on bra center
    integer order_geo      !incremental recorder over orders of geometric derivatives
    integer min_hgto_ket   !minimum order of HGTOs on ket center in the recurrence relations
    integer dim_hgto_ket   !dimension of HGTOs on ket center
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
    ! initializes the maximum order and dimension of HGTOs on bra center
    max_hgto_bra = 0
    dim_hgto_bra = 0
    ! loops over the orders of geometric derivatives not returned, the maximum
    ! order of HGTOs \var(max_hgto_bra) needs to update in the cycle
    do order_geo = max_order_geo-1, min_order_geo, -1
      ! updates the number of geometric derivatives
      num_geo_pot = num_geo_pot-(order_geo+2)  !=(order_geo+1)*(order_geo+2)/2
      ! updates maximum order of HGTOs on bra center
      max_hgto_bra = max_hgto_bra+1
      ! updates the dimension of HGTOs on bra center
      dim_hgto_bra = dim_hgto_bra+(max_hgto_bra+1)*(max_hgto_bra+2)/2
      ! sets minimum order of HGTOs on ket center
      min_hgto_ket = max(min_order_hket-max_hgto_bra,0)
      ! sets the dimension of HGTOs on ket center
      if (min_hgto_ket==0) then
        dim_hgto_ket = max_dim_hket
      else
        dim_hgto_ket = max_dim_hket-min_hgto_ket*(min_hgto_ket+1)*(min_hgto_ket+2)/6
      end if
      ! updates the maximum dimension
      dim_tmp = dim_hgto_bra*dim_hgto_ket*num_geo_pot
      if (dim_tmp>dim_ints) dim_ints = dim_tmp
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "dim_nucpot_hbra", STDOUT)
#endif
    return
  end subroutine dim_nucpot_hbra

  !> \brief sub-recurrence relations by recovering upper order HGTOs on bra center
  !> \author Bin Gao
  !> \date 2011-10-18
  !> \param cur_order_geo is current order of geometric derivatives
  !> \param orders_hgto_ket contains the minimum and maximum orders HGTOs on ket center
  !> \param min_order_hket is minimum of current order of HGTOs on ket center
  !> \param max_order_hbra is maximum of current order of HGTOs on bra center
  !> \param cc_wrt_bra contains the relative coordinates of center-of-charge w.r.t. bra center
  !> \param half_neg_rp is the half of the negative reciprocal of total exponent \f$p_{ij}\f$
  !> \param ket_to_bra is the ratio of exponent on ket center to that on bra center
  !> \param dim_cur_hket is the dimension of HGTOs on ket center for integrals with
  !>        upper order geometric derivatives and zeroth order HGTOs on bra center
  !> \param num_up_geo is the number of upper order geometric derivatives
  !> \param up_geo_zero contains the integrals with zeroth order HGTO on bra center and
  !>        upper order geometric derivatives
  !> \param dim_hket_zero is the dimension of HGTOs on ket center for integrals with
  !>        lower order geometric derivatives and zeroth order HGTO on bra center
  !> \param num_low_geo is the number lower order geometric derivatives
  !> \param low_geo_zero contains the integrals with zeroth order HGTO on bra center and
  !>        lower order geometric derivatives
  !> \param offset_hket_low is the offset of HGTOs on ket center in temporary integrals with
  !>        lower order HGTOs on bra center
  !> \param dim_low_hbra is the dimension of lower order HGTOs on bra center
  !> \param dim_hket_low is the dimension of HGTOs on ket center with upper order
  !>        geometric derivatives and lower order HGTOs on bra center
  !> \param low_hbra_pints contains the integrals with lower order HGTOs on bra center and
  !>        upper order geometric derivatives
  !> \param dim_up_hbra is the dimension of upper order HGTOs on bra center
  !> \return up_hbra_pints contains the integrals with upper order HGTOs on bra center and
  !>         lower order geometric derivatives
  subroutine sub_nucpot_hbra(cur_order_geo, orders_hgto_ket, min_order_hket, &
                             max_order_hbra, cc_wrt_bra, half_neg_rp,        &
                             ket_to_bra, dim_cur_hket, num_up_geo,           &
                             up_geo_zero, dim_hket_zero, num_low_geo,        &
                             low_geo_zero, offset_hket_low, dim_low_hbra,    &
                             dim_hket_low, low_hbra_pints, dim_up_hbra, up_hbra_pints)
    use xkind
    implicit none
    integer, intent(in) :: cur_order_geo
    integer, intent(in) :: orders_hgto_ket(2)
    integer, intent(in) :: min_order_hket
    integer, intent(in) :: max_order_hbra
    real(REALK), intent(in) :: cc_wrt_bra(3)
    real(REALK), intent(in) :: half_neg_rp
    real(REALK), intent(in) :: ket_to_bra
    integer, intent(in) :: dim_cur_hket
    integer, intent(in) :: num_up_geo
    real(REALK), intent(in) :: up_geo_zero(dim_cur_hket,num_up_geo)
    integer, intent(in) :: dim_hket_zero
    integer, intent(in) :: num_low_geo
    real(REALK), intent(in) :: low_geo_zero(dim_hket_zero,num_low_geo)
    integer, intent(in) :: offset_hket_low
    integer, intent(in) :: dim_low_hbra
    integer, intent(in) :: dim_hket_low
    real(REALK), intent(in) :: low_hbra_pints(dim_low_hbra,dim_hket_low,num_up_geo)
    integer, intent(in) :: dim_up_hbra
    real(REALK), intent(out) :: up_hbra_pints(dim_up_hbra,dim_cur_hket,num_low_geo)
!f2py intent(in) :: cur_order_geo
!f2py intent(in) :: orders_hgto_ket
!f2py intent(in) :: min_order_hket
!f2py intent(in) :: max_order_hbra
!f2py intent(in) :: cc_wrt_bra
!f2py intent(in) :: half_neg_rp
!f2py intent(in) :: ket_to_bra
!f2py intent(hide) :: dim_cur_hket
!f2py intent(hide) :: num_up_geo
!f2py intent(in) :: up_geo_zero
!f2py intent(hide) :: dim_hket_zero
!f2py intent(hide) :: num_low_geo
!f2py intent(in) :: low_geo_zero
!f2py intent(in) :: offset_hket_low
!f2py intent(hide) :: dim_low_hbra
!f2py intent(hide) :: dim_hket_low
!f2py intent(in) :: low_hbra_pints
!f2py depend(num_up_geo) :: low_hbra_pints
!f2py intent(in) :: dim_up_hbra
!f2py intent(out) :: up_hbra_pints
!f2py depend(dim_up_hbra) :: up_hbra_pints
!f2py depend(dim_cur_hket) :: up_hbra_pints
!f2py depend(num_low_geo) :: up_hbra_pints
    integer addr_up_geo     !addresses of upper order geometric derivatives
    integer addr_up_geo_y
    integer addr_up_geo_z
    integer addr_cur_geo    !address of current order geometric derivatives
    integer igeo, jgeo      !incremental recorders over geometric derivatives
    integer min_cur_hket    !minimum of current order of HGTOs on ket center
    logical zero_cur_hket   !if the minimum order of HGTOs on ket center is zeroth order
    integer max_cur_hbra    !maximum of current order of HGTOs on bra center for given order HGTOs on ket center
    integer base_low_hket   !base address of lower order HGTOs on ket center
    integer base_hket_zero  !base address of lower order HGTOs on ket center in integrals with zeroth
                            !order HGTO on bra center and lower order geometric derivatives
    integer addr_cur_hket   !address of current order HGTOs on ket center
    integer addr_hket_zero  !address of current order HGTOs on ket center in integrals with zeroth
                            !order HGTO on bra center and lower order geometric derivatives
    integer addr_hket_low   !address of current order HGTOs on ket center in integrals with lower
                            !order HGTOs on bra center and lower order geometric derivatives
    integer addr_low_xket   !addresses of lower order HGTOs on ket center
    integer addr_low_yket
    integer addr_low_zket
    integer addr_xket_zero  !addresses of lower order HGTOs on ket center in integrals with zeroth
    integer addr_yket_zero  !order HGTO on bra center and lower order geometric derivatives
    integer addr_zket_zero
    integer order_hket      !order of HGTOs on ket center
    integer order_xket      !orders of xyz components of HGTOs on ket center
    integer order_yket
    integer order_zket
    integer addr_up_hbra    !address of upper order HGTOs on bra center
    integer addr_cur_hbra   !address of current order HGTOs on bra center
    integer addr_low_hbra   !address of lower order HGTOs on bra center
    integer order_hbra      !order of HGTOs on bra center
    integer ibra, jbra      !incremental recorder over HGTOs on bra center
#if defined(XTIME)
    real(REALK) curr_time   !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (min_order_hket==0) then
      zero_cur_hket = .true.
      min_cur_hket = 1
    else
      zero_cur_hket = .false.
      min_cur_hket = min_order_hket
    end if
    addr_up_geo = 0
    addr_cur_geo = 0
    ! loops over xyz components of geometric derivatives
    do igeo = cur_order_geo, 0, -1
      do jgeo = 0, igeo
        addr_up_geo = addr_up_geo+1
        addr_cur_geo = addr_cur_geo+1
        addr_cur_hket = 0
        addr_hket_low = offset_hket_low
        ! zeroth order HGTOs on ket center
        if (zero_cur_hket) then
          addr_cur_hket = addr_cur_hket+1
          addr_hket_zero = 1
          ! sets maximum of current order of HGTOs on bra center for zeroth order HGTOs on ket center
          max_cur_hbra = max_order_hbra-orders_hgto_ket(1)
          ! px HGTO on bra center
          up_hbra_pints(1,addr_cur_hket,addr_cur_geo)                 &
            = cc_wrt_bra(1)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
            + half_neg_rp*up_geo_zero(addr_cur_hket,addr_up_geo)
          ! py HGTO on bra center
          addr_up_geo_y = addr_up_geo+1
          up_hbra_pints(2,addr_cur_hket,addr_cur_geo)                 &
            = cc_wrt_bra(2)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
            + half_neg_rp*up_geo_zero(addr_cur_hket,addr_up_geo_y)
          ! pz HGTO on bra center
          addr_up_geo_z = addr_up_geo+igeo+2
          up_hbra_pints(3,addr_cur_hket,addr_cur_geo)                 &
            = cc_wrt_bra(3)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
            + half_neg_rp*up_geo_zero(addr_cur_hket,addr_up_geo_z)
          ! d-shell on bra center
          if (max_cur_hbra>0) then
            addr_hket_low = addr_hket_low+1
            ! dxx
            up_hbra_pints(4,addr_cur_hket,addr_cur_geo)                           &
              = cc_wrt_bra(1)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo)         &
              + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + low_hbra_pints(1,addr_hket_low,addr_up_geo))
            ! dxy
            up_hbra_pints(5,addr_cur_hket,addr_cur_geo)                   &
              = cc_wrt_bra(2)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo) &
              + half_neg_rp*low_hbra_pints(1,addr_hket_low,addr_up_geo_y)
            ! dyy
            up_hbra_pints(6,addr_cur_hket,addr_cur_geo)                           &
              = cc_wrt_bra(2)*up_hbra_pints(2,addr_cur_hket,addr_cur_geo)         &
              + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + low_hbra_pints(2,addr_hket_low,addr_up_geo_y))
            ! dxz
            up_hbra_pints(7,addr_cur_hket,addr_cur_geo)                   &
              = cc_wrt_bra(3)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo) &
              + half_neg_rp*low_hbra_pints(1,addr_hket_low,addr_up_geo_z)
            ! dyz
            up_hbra_pints(8,addr_cur_hket,addr_cur_geo)                   &
              = cc_wrt_bra(3)*up_hbra_pints(2,addr_cur_hket,addr_cur_geo) &
              + half_neg_rp*low_hbra_pints(2,addr_hket_low,addr_up_geo_z)
            ! dzz
            up_hbra_pints(9,addr_cur_hket,addr_cur_geo)                           &
              = cc_wrt_bra(3)*up_hbra_pints(3,addr_cur_hket,addr_cur_geo)         &
              + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + low_hbra_pints(3,addr_hket_low,addr_up_geo_z))
          end if
          if (max_cur_hbra>1) then
            addr_up_hbra = 9   !base address of the f-shell HGTOs on bra center
            addr_cur_hbra = 3  !base address of the d-shell HGTOs on bra center
            addr_low_hbra = 0  !base address of the p-shell HGTOs on bra center
            ! loops over other current order of HGTOs on bra center, starting from d-shell
            do order_hbra = 2, max_cur_hbra
              addr_up_hbra = addr_up_hbra+1
              addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
              write(STDOUT,100) "orders:", order_hbra, order_hket, cur_order_geo
              write(STDOUT,100) "x-direction"
              write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
              write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
              write(STDOUT,110) addr_low_hbra+1, addr_cur_hket, addr_cur_geo
              write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo
              write(STDOUT,100) "------------------------------------"
#endif
              ! x-direction
              up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)        &
                = cc_wrt_bra(1)                                             &
                * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo)   &
                + half_neg_rp*(real(order_hbra,REALK)*ket_to_bra            &
                * up_hbra_pints(addr_low_hbra+1,addr_cur_hket,addr_cur_geo) &
                + low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo))
              ! y-direction
              addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
              write(STDOUT,100) "y-direction"
              write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
              write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
              write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+1
              write(STDOUT,100) "------------------------------------"
#endif
              up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)      &
                = cc_wrt_bra(2)                                           &
                * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo) &
                + half_neg_rp                                             &
                * low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo+1)
              do ibra = 1, order_hbra
                addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_low_hbra+ibra, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_up_geo+1
                write(STDOUT,100) "------------------------------------"
#endif
                up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)           &
                  = cc_wrt_bra(2)                                                &
                  * up_hbra_pints(addr_cur_hbra+ibra,addr_cur_hket,addr_cur_geo) &
                  + half_neg_rp*(real(ibra,REALK)*ket_to_bra                     &
                  * up_hbra_pints(addr_low_hbra+ibra,addr_cur_hket,addr_cur_geo) &
                  + low_hbra_pints(addr_cur_hbra+ibra,addr_hket_low,addr_up_geo+1))
              end do
              ! z-direction
#if defined(DEBUG)
              write(STDOUT,100) "z-direction"
#endif
              addr_cur_hbra = addr_cur_hbra-1
              do jbra = 0, order_hbra
                addr_up_hbra = addr_up_hbra+1
                addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+igeo+2
                write(STDOUT,100) "------------------------------------"
#endif
                up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)      &
                  = cc_wrt_bra(3)                                           &
                  * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo) &
                  + half_neg_rp                                             &
                  * low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo+igeo+2)
              end do
              do ibra = 1, order_hbra
                do jbra = 0, order_hbra-ibra
                  addr_up_hbra = addr_up_hbra+1
                  addr_cur_hbra = addr_cur_hbra+1
                  addr_low_hbra = addr_low_hbra+1
#if defined(DEBUG)
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_low_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+igeo+2
                  write(STDOUT,100) "------------------------------------"
#endif
                  up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)      &
                    = cc_wrt_bra(3)                                           &
                    * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo) &
                    + half_neg_rp*(real(ibra,REALK)*ket_to_bra                &
                    * up_hbra_pints(addr_low_hbra,addr_cur_hket,addr_cur_geo) &
                    + low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo+igeo+2))
                end do
              end do
            end do
          end if
        else
          addr_hket_zero = min_order_hket*(min_order_hket+1)/2
        end if
        base_low_hket = 0
        base_hket_zero = 0
        ! loops over other current order HGTOs on ket center
        do order_hket = min_cur_hket, orders_hgto_ket(2)
          addr_cur_hket = addr_cur_hket+1
          addr_hket_zero = addr_hket_zero+1
          addr_xket_zero = base_hket_zero+1
          ! sets maximum of current order of HGTOs on bra center for order \var(order_hket)
          ! HGTOs on ket center
          max_cur_hbra = min(max_order_hbra,max_order_hbra-orders_hgto_ket(1)+order_hket)
          ! (1) x...x component of upper order HGTOs on ket center
          !
          ! px HGTO on bra center
          up_hbra_pints(1,addr_cur_hket,addr_cur_geo)                 &
            = cc_wrt_bra(1)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
            + half_neg_rp*(up_geo_zero(addr_cur_hket,addr_up_geo)     &
            - real(order_hket,REALK)*low_geo_zero(addr_xket_zero,addr_cur_geo))
          ! py HGTO on bra center
          addr_up_geo_y = addr_up_geo+1
          up_hbra_pints(2,addr_cur_hket,addr_cur_geo)                 &
            = cc_wrt_bra(2)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
            + half_neg_rp*up_geo_zero(addr_cur_hket,addr_up_geo_y)
          ! pz HGTO on bra center
          addr_up_geo_z = addr_up_geo+igeo+2
          up_hbra_pints(3,addr_cur_hket,addr_cur_geo)                 &
            = cc_wrt_bra(3)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
            + half_neg_rp*up_geo_zero(addr_cur_hket,addr_up_geo_z)
          ! d-shell on bra center
          addr_low_xket = base_low_hket+1  !should put after if statement, but to avoid warning from compiler
          if (max_cur_hbra>0) then
            addr_hket_low = addr_hket_low+1
            ! dxx
            up_hbra_pints(4,addr_cur_hket,addr_cur_geo)                           &
              = cc_wrt_bra(1)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo)         &
              + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + low_hbra_pints(1,addr_hket_low,addr_up_geo)                       &
              - real(order_hket,REALK)*up_hbra_pints(1,addr_low_xket,addr_cur_geo))
            ! dxy
            up_hbra_pints(5,addr_cur_hket,addr_cur_geo)                   &
              = cc_wrt_bra(2)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo) &
              + half_neg_rp*low_hbra_pints(1,addr_hket_low,addr_up_geo_y)
            ! dyy
            up_hbra_pints(6,addr_cur_hket,addr_cur_geo)                           &
              = cc_wrt_bra(2)*up_hbra_pints(2,addr_cur_hket,addr_cur_geo)         &
              + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + low_hbra_pints(2,addr_hket_low,addr_up_geo_y))
            ! dxz
            up_hbra_pints(7,addr_cur_hket,addr_cur_geo)                   &
              = cc_wrt_bra(3)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo) &
              + half_neg_rp*low_hbra_pints(1,addr_hket_low,addr_up_geo_z)
            ! dyz
            up_hbra_pints(8,addr_cur_hket,addr_cur_geo)                   &
              = cc_wrt_bra(3)*up_hbra_pints(2,addr_cur_hket,addr_cur_geo) &
              + half_neg_rp*low_hbra_pints(2,addr_hket_low,addr_up_geo_z)
            ! dzz
            up_hbra_pints(9,addr_cur_hket,addr_cur_geo)                           &
              = cc_wrt_bra(3)*up_hbra_pints(3,addr_cur_hket,addr_cur_geo)         &
              + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + low_hbra_pints(3,addr_hket_low,addr_up_geo_z))
            if (max_cur_hbra>1) then
              addr_up_hbra = 9   !base address of the f-shell HGTOs on bra center
              addr_cur_hbra = 3  !base address of the d-shell HGTOs on bra center
              addr_low_hbra = 0  !base address of the p-shell HGTOs on bra center
              ! loops over other current order of HGTOs on bra center, starting from d-shell
              do order_hbra = 2, max_cur_hbra
                addr_up_hbra = addr_up_hbra+1
                addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                write(STDOUT,100) "orders:", order_hbra, order_hket, cur_order_geo
                write(STDOUT,100) "x-direction"
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_low_hbra+1, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo
                write(STDOUT,110) addr_up_hbra, addr_low_xket, addr_cur_geo
                write(STDOUT,100) "------------------------------------"
#endif
                ! x-direction
                up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)        &
                  = cc_wrt_bra(1)                                             &
                  * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo)   &
                  + half_neg_rp*(real(order_hbra,REALK)*ket_to_bra            &
                  * up_hbra_pints(addr_low_hbra+1,addr_cur_hket,addr_cur_geo) &
                  + low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo)   &
                  - real(order_hket,REALK)                                    &
                  * up_hbra_pints(addr_cur_hbra,addr_low_xket,addr_cur_geo))
                ! y-direction
                addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                write(STDOUT,100) "y-direction"
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+1
                write(STDOUT,100) "------------------------------------"
#endif
                up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)      &
                  = cc_wrt_bra(2)                                           &
                  * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo) &
                  + half_neg_rp                                             &
                  * low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo+1)
                do ibra = 1, order_hbra
                  addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_low_hbra+ibra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_up_geo+1
                  write(STDOUT,100) "------------------------------------"
#endif
                  up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)           &
                    = cc_wrt_bra(2)                                                &
                    * up_hbra_pints(addr_cur_hbra+ibra,addr_cur_hket,addr_cur_geo) &
                    + half_neg_rp*(real(ibra,REALK)*ket_to_bra                     &
                    * up_hbra_pints(addr_low_hbra+ibra,addr_cur_hket,addr_cur_geo) &
                    + low_hbra_pints(addr_cur_hbra+ibra,addr_hket_low,addr_up_geo+1))
                end do
                ! z-direction
#if defined(DEBUG)
                write(STDOUT,100) "z-direction"
#endif
                addr_cur_hbra = addr_cur_hbra-1
                do jbra = 0, order_hbra
                  addr_up_hbra = addr_up_hbra+1
                  addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+igeo+2
                  write(STDOUT,100) "------------------------------------"
#endif
                  up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)      &
                    = cc_wrt_bra(3)                                           &
                    * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo) &
                    + half_neg_rp                                             &
                    * low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo+igeo+2)
                end do
                do ibra = 1, order_hbra
                  do jbra = 0, order_hbra-ibra
                    addr_up_hbra = addr_up_hbra+1
                    addr_cur_hbra = addr_cur_hbra+1
                    addr_low_hbra = addr_low_hbra+1
#if defined(DEBUG)
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_low_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+igeo+2
                    write(STDOUT,100) "------------------------------------"
#endif
                    up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)      &
                      = cc_wrt_bra(3)                                           &
                      * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo) &
                      + half_neg_rp*(real(ibra,REALK)*ket_to_bra                &
                      * up_hbra_pints(addr_low_hbra,addr_cur_hket,addr_cur_geo) &
                      + low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo+igeo+2))
                  end do
                end do
              end do
            end if
          end if
          ! (2) x...xy to xy...y components of upper order HGTOs on ket center
          addr_low_yket = base_low_hket
          addr_yket_zero = base_hket_zero
          do order_yket = 1, order_hket-1
            addr_cur_hket = addr_cur_hket+1
            addr_hket_zero = addr_hket_zero+1
            addr_xket_zero = addr_xket_zero+1
            addr_yket_zero = addr_yket_zero+1
            ! sets the order along x-direction
            order_xket = order_hket-order_yket
            ! px HGTO on bra center
            up_hbra_pints(1,addr_cur_hket,addr_cur_geo)                 &
              = cc_wrt_bra(1)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + half_neg_rp*(up_geo_zero(addr_cur_hket,addr_up_geo)     &
              - real(order_xket,REALK)*low_geo_zero(addr_xket_zero,addr_cur_geo))
            ! py HGTO on bra center
            addr_up_geo_y = addr_up_geo+1
            up_hbra_pints(2,addr_cur_hket,addr_cur_geo)                 &
              = cc_wrt_bra(2)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + half_neg_rp*(up_geo_zero(addr_cur_hket,addr_up_geo_y)   &
              - real(order_yket,REALK)*low_geo_zero(addr_yket_zero,addr_cur_geo))
            ! pz HGTO on bra center
            addr_up_geo_z = addr_up_geo+igeo+2
            up_hbra_pints(3,addr_cur_hket,addr_cur_geo)                 &
              = cc_wrt_bra(3)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + half_neg_rp*up_geo_zero(addr_cur_hket,addr_up_geo_z)
            ! d-shell on bra center
            if (max_cur_hbra>0) then
              addr_low_xket = addr_low_xket+1
              addr_low_yket = addr_low_yket+1
              addr_hket_low = addr_hket_low+1
              ! dxx
              up_hbra_pints(4,addr_cur_hket,addr_cur_geo)                           &
                = cc_wrt_bra(1)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo)         &
                + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
                + low_hbra_pints(1,addr_hket_low,addr_up_geo)                       &
                - real(order_xket,REALK)*up_hbra_pints(1,addr_low_xket,addr_cur_geo))
              ! dxy
              up_hbra_pints(5,addr_cur_hket,addr_cur_geo)                    &
                = cc_wrt_bra(2)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo)  &
                + half_neg_rp*(low_hbra_pints(1,addr_hket_low,addr_up_geo_y) &
                - real(order_yket,REALK)*up_hbra_pints(1,addr_low_yket,addr_cur_geo))
              ! dyy
              up_hbra_pints(6,addr_cur_hket,addr_cur_geo)                           &
                = cc_wrt_bra(2)*up_hbra_pints(2,addr_cur_hket,addr_cur_geo)         &
                + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
                + low_hbra_pints(2,addr_hket_low,addr_up_geo_y)                     &
                - real(order_yket,REALK)*up_hbra_pints(2,addr_low_yket,addr_cur_geo))
              ! dxz
              up_hbra_pints(7,addr_cur_hket,addr_cur_geo)                   &
                = cc_wrt_bra(3)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo) &
                + half_neg_rp*low_hbra_pints(1,addr_hket_low,addr_up_geo_z)
              ! dyz
              up_hbra_pints(8,addr_cur_hket,addr_cur_geo)                   &
                = cc_wrt_bra(3)*up_hbra_pints(2,addr_cur_hket,addr_cur_geo) &
                + half_neg_rp*low_hbra_pints(2,addr_hket_low,addr_up_geo_z)
              ! dzz
              up_hbra_pints(9,addr_cur_hket,addr_cur_geo)                           &
                = cc_wrt_bra(3)*up_hbra_pints(3,addr_cur_hket,addr_cur_geo)         &
                + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
                + low_hbra_pints(3,addr_hket_low,addr_up_geo_z))
              if (max_cur_hbra>1) then
                addr_up_hbra = 9   !base address of the f-shell HGTOs on bra center
                addr_cur_hbra = 3  !base address of the d-shell HGTOs on bra center
                addr_low_hbra = 0  !base address of the p-shell HGTOs on bra center
                ! loops over other current order of HGTOs on bra center, starting from d-shell
                do order_hbra = 2, max_cur_hbra
                  addr_up_hbra = addr_up_hbra+1
                  addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                  write(STDOUT,100) "orders:", order_hbra, order_xket, cur_order_geo
                  write(STDOUT,100) "x-direction"
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_low_hbra+1, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo
                  write(STDOUT,110) addr_cur_hbra, addr_low_xket, addr_cur_geo
                  write(STDOUT,100) "------------------------------------"
#endif
                  ! x-direction
                  up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)        &
                    = cc_wrt_bra(1)                                             &
                    * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo)   &
                    + half_neg_rp*(real(order_hbra,REALK)*ket_to_bra            &
                    * up_hbra_pints(addr_low_hbra+1,addr_cur_hket,addr_cur_geo) &
                    + low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo)   &
                    - real(order_xket,REALK)                                    &
                    * up_hbra_pints(addr_cur_hbra,addr_low_xket,addr_cur_geo))
                  ! y-direction
                  addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                  write(STDOUT,100) "y-direction"
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+1
                  write(STDOUT,110) addr_up_hbra, addr_low_yket, addr_cur_geo
                  write(STDOUT,100) "------------------------------------"
#endif
                  up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)         &
                    = cc_wrt_bra(2)                                              &
                    * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo)    &
                    + half_neg_rp                                                &
                    * (low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo+1) &
                    - real(order_yket,REALK)                                     &
                    * up_hbra_pints(addr_cur_hbra,addr_low_yket,addr_cur_geo))
                  do ibra = 1, order_hbra
                    addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_low_hbra+ibra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_up_geo+1
                    write(STDOUT,110) addr_cur_hbra+ibra, addr_low_yket, addr_cur_geo
                    write(STDOUT,100) "------------------------------------"
#endif
                    up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)             &
                      = cc_wrt_bra(2)                                                  &
                      * up_hbra_pints(addr_cur_hbra+ibra,addr_cur_hket,addr_cur_geo)   &
                      + half_neg_rp*(real(ibra,REALK)*ket_to_bra                       &
                      * up_hbra_pints(addr_low_hbra+ibra,addr_cur_hket,addr_cur_geo)   &
                      + low_hbra_pints(addr_cur_hbra+ibra,addr_hket_low,addr_up_geo+1) &
                      - real(order_yket,REALK)                                         &
                      * up_hbra_pints(addr_cur_hbra+ibra,addr_low_yket,addr_cur_geo))
                  end do
                  ! z-direction
#if defined(DEBUG)
                  write(STDOUT,100) "z-direction"
#endif
                  addr_cur_hbra = addr_cur_hbra-1
                  do jbra = 0, order_hbra
                    addr_up_hbra = addr_up_hbra+1
                    addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+igeo+2
                    write(STDOUT,100) "------------------------------------"
#endif
                    up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)      &
                      = cc_wrt_bra(3)                                           &
                      * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo) &
                      + half_neg_rp                                             &
                      * low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo+igeo+2)
                  end do
                  do ibra = 1, order_hbra
                    do jbra = 0, order_hbra-ibra
                      addr_up_hbra = addr_up_hbra+1
                      addr_cur_hbra = addr_cur_hbra+1
                      addr_low_hbra = addr_low_hbra+1
#if defined(DEBUG)
                      write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                      write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                      write(STDOUT,110) addr_low_hbra, addr_cur_hket, addr_cur_geo
                      write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+igeo+2
                      write(STDOUT,100) "------------------------------------"
#endif
                      up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)      &
                        = cc_wrt_bra(3)                                           &
                        * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo) &
                        + half_neg_rp*(real(ibra,REALK)*ket_to_bra                &
                        * up_hbra_pints(addr_low_hbra,addr_cur_hket,addr_cur_geo) &
                        + low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo+igeo+2))
                    end do
                  end do
                end do
              end if
            end if
          end do
          ! (3) y...y component of upper order HGTOs on ket center
          addr_cur_hket = addr_cur_hket+1
          addr_hket_zero = addr_hket_zero+1
          addr_yket_zero = addr_yket_zero+1
          ! px HGTO on bra center
          up_hbra_pints(1,addr_cur_hket,addr_cur_geo)                 &
            = cc_wrt_bra(1)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
            + half_neg_rp*up_geo_zero(addr_cur_hket,addr_up_geo)
          ! py HGTO on bra center
          addr_up_geo_y = addr_up_geo+1
          up_hbra_pints(2,addr_cur_hket,addr_cur_geo)                 &
            = cc_wrt_bra(2)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
            + half_neg_rp*(up_geo_zero(addr_cur_hket,addr_up_geo_y)   &
            - real(order_hket,REALK)*low_geo_zero(addr_yket_zero,addr_cur_geo))
          ! pz HGTO on bra center
          addr_up_geo_z = addr_up_geo+igeo+2
          up_hbra_pints(3,addr_cur_hket,addr_cur_geo)                 &
            = cc_wrt_bra(3)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
            + half_neg_rp*up_geo_zero(addr_cur_hket,addr_up_geo_z)
          ! d-shell on bra center
          if (max_cur_hbra>0) then
            addr_low_yket = addr_low_yket+1
            addr_hket_low = addr_hket_low+1
            ! dxx
            up_hbra_pints(4,addr_cur_hket,addr_cur_geo)                           &
              = cc_wrt_bra(1)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo)         &
              + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + low_hbra_pints(1,addr_hket_low,addr_up_geo))
            ! dxy
            up_hbra_pints(5,addr_cur_hket,addr_cur_geo)                    &
              = cc_wrt_bra(2)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo)  &
              + half_neg_rp*(low_hbra_pints(1,addr_hket_low,addr_up_geo_y) &
              - real(order_hket,REALK)*up_hbra_pints(1,addr_low_yket,addr_cur_geo))
            ! dyy
            up_hbra_pints(6,addr_cur_hket,addr_cur_geo)                           &
              = cc_wrt_bra(2)*up_hbra_pints(2,addr_cur_hket,addr_cur_geo)         &
              + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + low_hbra_pints(2,addr_hket_low,addr_up_geo_y)                     &
              - real(order_hket,REALK)*up_hbra_pints(2,addr_low_yket,addr_cur_geo))
            ! dxz
            up_hbra_pints(7,addr_cur_hket,addr_cur_geo)                   &
              = cc_wrt_bra(3)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo) &
              + half_neg_rp*low_hbra_pints(1,addr_hket_low,addr_up_geo_z)
            ! dyz
            up_hbra_pints(8,addr_cur_hket,addr_cur_geo)                   &
              = cc_wrt_bra(3)*up_hbra_pints(2,addr_cur_hket,addr_cur_geo) &
              + half_neg_rp*low_hbra_pints(2,addr_hket_low,addr_up_geo_z)
            ! dzz
            up_hbra_pints(9,addr_cur_hket,addr_cur_geo)                           &
              = cc_wrt_bra(3)*up_hbra_pints(3,addr_cur_hket,addr_cur_geo)         &
              + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + low_hbra_pints(3,addr_hket_low,addr_up_geo_z))
            if (max_cur_hbra>1) then
              addr_up_hbra = 9   !base address of the f-shell HGTOs on bra center
              addr_cur_hbra = 3  !base address of the d-shell HGTOs on bra center
              addr_low_hbra = 0  !base address of the p-shell HGTOs on bra center
              ! loops over other current order of HGTOs on bra center, starting from d-shell
              do order_hbra = 2, max_cur_hbra
                addr_up_hbra = addr_up_hbra+1
                addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                write(STDOUT,100) "orders:", order_hbra, order_xket, cur_order_geo
                write(STDOUT,100) "x-direction"
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_low_hbra+1, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo
                write(STDOUT,100) "------------------------------------"
#endif
                ! x-direction
                up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)        &
                  = cc_wrt_bra(1)                                             &
                  * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo)   &
                  + half_neg_rp*(real(order_hbra,REALK)*ket_to_bra            &
                  * up_hbra_pints(addr_low_hbra+1,addr_cur_hket,addr_cur_geo) &
                  + low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo))
                ! y-direction
                addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                write(STDOUT,100) "y-direction"
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+1
                write(STDOUT,110) addr_cur_hbra, addr_low_yket, addr_cur_geo
                write(STDOUT,100) "------------------------------------"
#endif
                up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)         &
                  = cc_wrt_bra(2)                                              &
                  * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo)    &
                  + half_neg_rp                                                &
                  * (low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo+1) &
                  - real(order_hket,REALK)                                     &
                  * up_hbra_pints(addr_cur_hbra,addr_low_yket,addr_cur_geo))
                do ibra = 1, order_hbra
                  addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_low_hbra+ibra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_up_geo+1
                  write(STDOUT,110) addr_cur_hbra+ibra, addr_low_yket, addr_cur_geo
                  write(STDOUT,100) "------------------------------------"
#endif
                  up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)             &
                    = cc_wrt_bra(2)                                                  &
                    * up_hbra_pints(addr_cur_hbra+ibra,addr_cur_hket,addr_cur_geo)   &
                    + half_neg_rp*(real(ibra,REALK)*ket_to_bra                       &
                    * up_hbra_pints(addr_low_hbra+ibra,addr_cur_hket,addr_cur_geo)   &
                    + low_hbra_pints(addr_cur_hbra+ibra,addr_hket_low,addr_up_geo+1) &
                    - real(order_hket,REALK)                                         &
                    * up_hbra_pints(addr_cur_hbra+ibra,addr_low_yket,addr_cur_geo))
                end do
                ! z-direction
#if defined(DEBUG)
                write(STDOUT,100) "z-direction"
#endif
                addr_cur_hbra = addr_cur_hbra-1
                do jbra = 0, order_hbra
                  addr_up_hbra = addr_up_hbra+1
                  addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+igeo+2
                  write(STDOUT,100) "------------------------------------"
#endif
                  up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)      &
                    = cc_wrt_bra(3)                                           &
                    * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo) &
                    + half_neg_rp                                             &
                    * low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo+igeo+2)
                end do
                do ibra = 1, order_hbra
                  do jbra = 0, order_hbra-ibra
                    addr_up_hbra = addr_up_hbra+1
                    addr_cur_hbra = addr_cur_hbra+1
                    addr_low_hbra = addr_low_hbra+1
#if defined(DEBUG)
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_low_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+igeo+2
                    write(STDOUT,100) "------------------------------------"
#endif
                    up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)      &
                      = cc_wrt_bra(3)                                           &
                      * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo) &
                      + half_neg_rp*(real(ibra,REALK)*ket_to_bra                &
                      * up_hbra_pints(addr_low_hbra,addr_cur_hket,addr_cur_geo) &
                      + low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo+igeo+2))
                  end do
                end do
              end do
            end if
          end if
          ! (4) x...xz to yz...z components of upper order HGTOs on ket center
          addr_low_xket = base_low_hket+order_hket
          addr_low_zket = base_low_hket
          addr_xket_zero = base_hket_zero+order_hket
          addr_zket_zero = base_hket_zero
          do order_zket = 1, order_hket-1
            addr_cur_hket = addr_cur_hket+1
            addr_hket_zero = addr_hket_zero+1
            addr_xket_zero = addr_xket_zero+1
            addr_zket_zero = addr_zket_zero+1
            ! (4.1) x...xz...z component of upper order HGTOs on ket center
            order_xket = order_hket-order_zket
            ! px HGTO on bra center
            up_hbra_pints(1,addr_cur_hket,addr_cur_geo)                 &
              = cc_wrt_bra(1)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + half_neg_rp*(up_geo_zero(addr_cur_hket,addr_up_geo)     &
              - real(order_xket,REALK)*low_geo_zero(addr_xket_zero,addr_cur_geo))
            ! py HGTO on bra center
            addr_up_geo_y = addr_up_geo+1
            up_hbra_pints(2,addr_cur_hket,addr_cur_geo)                 &
              = cc_wrt_bra(2)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + half_neg_rp*up_geo_zero(addr_cur_hket,addr_up_geo_y)
            ! pz HGTO on bra center
            addr_up_geo_z = addr_up_geo+igeo+2
            up_hbra_pints(3,addr_cur_hket,addr_cur_geo)                 &
              = cc_wrt_bra(3)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + half_neg_rp*(up_geo_zero(addr_cur_hket,addr_up_geo_z)   &
              - real(order_zket,REALK)*low_geo_zero(addr_zket_zero,addr_cur_geo))
            ! d-shell on bra center
            if (max_cur_hbra>0) then
              addr_low_xket = addr_low_xket+1
              addr_low_zket = addr_low_zket+1
              addr_hket_low = addr_hket_low+1
              ! dxx
              up_hbra_pints(4,addr_cur_hket,addr_cur_geo)                           &
                = cc_wrt_bra(1)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo)         &
                + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
                + low_hbra_pints(1,addr_hket_low,addr_up_geo)                       &
                - real(order_xket,REALK)*up_hbra_pints(1,addr_low_xket,addr_cur_geo))
              ! dxy
              up_hbra_pints(5,addr_cur_hket,addr_cur_geo)                   &
                = cc_wrt_bra(2)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo) &
                + half_neg_rp*low_hbra_pints(1,addr_hket_low,addr_up_geo_y)
              ! dyy
              up_hbra_pints(6,addr_cur_hket,addr_cur_geo)                           &
                = cc_wrt_bra(2)*up_hbra_pints(2,addr_cur_hket,addr_cur_geo)         &
                + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
                + low_hbra_pints(2,addr_hket_low,addr_up_geo_y))
              ! dxz
              up_hbra_pints(7,addr_cur_hket,addr_cur_geo)                    &
                = cc_wrt_bra(3)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo)  &
                + half_neg_rp*(low_hbra_pints(1,addr_hket_low,addr_up_geo_z) &
                - real(order_zket,REALK)*up_hbra_pints(1,addr_low_zket,addr_cur_geo))
              ! dyz
              up_hbra_pints(8,addr_cur_hket,addr_cur_geo)                    &
                = cc_wrt_bra(3)*up_hbra_pints(2,addr_cur_hket,addr_cur_geo)  &
                + half_neg_rp*(low_hbra_pints(2,addr_hket_low,addr_up_geo_z) &
                - real(order_zket,REALK)*up_hbra_pints(2,addr_low_zket,addr_cur_geo))
              ! dzz
              up_hbra_pints(9,addr_cur_hket,addr_cur_geo)                           &
                = cc_wrt_bra(3)*up_hbra_pints(3,addr_cur_hket,addr_cur_geo)         &
                + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
                + low_hbra_pints(3,addr_hket_low,addr_up_geo_z)                     &
                - real(order_zket,REALK)*up_hbra_pints(3,addr_low_zket,addr_cur_geo))
              if (max_cur_hbra>1) then
                addr_up_hbra = 9   !base address of the f-shell HGTOs on bra center
                addr_cur_hbra = 3  !base address of the d-shell HGTOs on bra center
                addr_low_hbra = 0  !base address of the p-shell HGTOs on bra center
                ! loops over other current order of HGTOs on bra center, starting from d-shell
                do order_hbra = 2, max_cur_hbra
                  addr_up_hbra = addr_up_hbra+1
                  addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                  write(STDOUT,100) "orders:", order_hbra, order_xket, cur_order_geo
                  write(STDOUT,100) "x-direction"
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_low_hbra+1, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra, addr_low_xket, addr_up_geo
                  write(STDOUT,100) "------------------------------------"
#endif
                  ! x-direction
                  up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)        &
                    = cc_wrt_bra(1)                                             &
                    * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo)   &
                    + half_neg_rp*(real(order_hbra,REALK)*ket_to_bra            &
                    * up_hbra_pints(addr_low_hbra+1,addr_cur_hket,addr_cur_geo) &
                    + low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo)   &
                    - real(order_xket,REALK)                                    &
                    * up_hbra_pints(addr_cur_hbra,addr_low_xket,addr_cur_geo))
                  ! y-direction
                  addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                  write(STDOUT,100) "y-direction"
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+1
                  write(STDOUT,100) "------------------------------------"
#endif
                  up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)      &
                    = cc_wrt_bra(2)                                           &
                    * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo) &
                    + half_neg_rp                                             &
                    * (low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo+1))
                  do ibra = 1, order_hbra
                    addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_low_hbra+ibra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_up_geo+1
                    write(STDOUT,100) "------------------------------------"
#endif
                    up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)           &
                      = cc_wrt_bra(2)                                                &
                      * up_hbra_pints(addr_cur_hbra+ibra,addr_cur_hket,addr_cur_geo) &
                      + half_neg_rp*(real(ibra,REALK)*ket_to_bra                     &
                      * up_hbra_pints(addr_low_hbra+ibra,addr_cur_hket,addr_cur_geo) &
                      + low_hbra_pints(addr_cur_hbra+ibra,addr_hket_low,addr_up_geo+1))
                  end do
                  ! z-direction
#if defined(DEBUG)
                  write(STDOUT,100) "z-direction"
#endif
                  addr_cur_hbra = addr_cur_hbra-1
                  do jbra = 0, order_hbra
                    addr_up_hbra = addr_up_hbra+1
                    addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+igeo+2
                    write(STDOUT,110) addr_cur_hbra, addr_low_zket, addr_cur_geo
                    write(STDOUT,100) "------------------------------------"
#endif
                    up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)              &
                      = cc_wrt_bra(3)                                                   &
                      * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo)         &
                      + half_neg_rp                                                     &
                      * (low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo+igeo+2) &
                      - real(order_zket,REALK)                                          &
                      * up_hbra_pints(addr_cur_hbra,addr_low_zket,addr_cur_geo))
                  end do
                  do ibra = 1, order_hbra
                    do jbra = 0, order_hbra-ibra
                      addr_up_hbra = addr_up_hbra+1
                      addr_cur_hbra = addr_cur_hbra+1
                      addr_low_hbra = addr_low_hbra+1
#if defined(DEBUG)
                      write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                      write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                      write(STDOUT,110) addr_low_hbra, addr_cur_hket, addr_cur_geo
                      write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+igeo+2
                      write(STDOUT,110) addr_cur_hbra, addr_low_zket, addr_cur_geo
                      write(STDOUT,100) "------------------------------------"
#endif
                      up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)      &
                        = cc_wrt_bra(3)                                           &
                        * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo) &
                        + half_neg_rp*(real(ibra,REALK)*ket_to_bra                &
                        * up_hbra_pints(addr_low_hbra,addr_cur_hket,addr_cur_geo) &
                        + low_hbra_pints(addr_cur_hbra,addr_hket_low,             &
                                         addr_up_geo+igeo+2)                      &
                        - real(order_zket,REALK)                                  &
                        * up_hbra_pints(addr_cur_hbra,addr_low_zket,addr_cur_geo))
                    end do
                  end do
                end do
              end if
            end if
            ! (4.2) x...xyz...z to xy...yz...z components of upper order HGTOs on ket center
            do order_yket = 1, order_hket-(order_zket+1)
              addr_cur_hket = addr_cur_hket+1
              addr_hket_zero = addr_hket_zero+1
              addr_xket_zero = addr_xket_zero+1
              addr_yket_zero = addr_yket_zero+1
              addr_zket_zero = addr_zket_zero+1
              order_xket = order_xket-1
              ! px HGTO on bra center
              up_hbra_pints(1,addr_cur_hket,addr_cur_geo)                 &
                = cc_wrt_bra(1)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
                + half_neg_rp*(up_geo_zero(addr_cur_hket,addr_up_geo)     &
                - real(order_xket,REALK)*low_geo_zero(addr_xket_zero,addr_cur_geo))
              ! py HGTO on bra center
              addr_up_geo_y = addr_up_geo+1
              up_hbra_pints(2,addr_cur_hket,addr_cur_geo)                 &
                = cc_wrt_bra(2)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
                + half_neg_rp*(up_geo_zero(addr_cur_hket,addr_up_geo_y)   &
                - real(order_yket,REALK)*low_geo_zero(addr_yket_zero,addr_cur_geo))
              ! pz HGTO on bra center
              addr_up_geo_z = addr_up_geo+igeo+2
              up_hbra_pints(3,addr_cur_hket,addr_cur_geo)                 &
                = cc_wrt_bra(3)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
                + half_neg_rp*(up_geo_zero(addr_cur_hket,addr_up_geo_z)   &
                - real(order_zket,REALK)*low_geo_zero(addr_zket_zero,addr_cur_geo))
              ! d-shell on bra center
              if (max_cur_hbra>0) then
                addr_low_xket = addr_low_xket+1
                addr_low_yket = addr_low_yket+1
                addr_low_zket = addr_low_zket+1
                addr_hket_low = addr_hket_low+1
                ! dxx
                up_hbra_pints(4,addr_cur_hket,addr_cur_geo)                           &
                  = cc_wrt_bra(1)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo)         &
                  + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
                  + low_hbra_pints(1,addr_hket_low,addr_up_geo)                       &
                  - real(order_xket,REALK)*up_hbra_pints(1,addr_low_xket,addr_cur_geo))
                ! dxy
                up_hbra_pints(5,addr_cur_hket,addr_cur_geo)                    &
                  = cc_wrt_bra(2)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo)  &
                  + half_neg_rp*(low_hbra_pints(1,addr_hket_low,addr_up_geo_y) &
                  - real(order_yket,REALK)*up_hbra_pints(1,addr_low_yket,addr_cur_geo))
                ! dyy
                up_hbra_pints(6,addr_cur_hket,addr_cur_geo)                   &
                  = cc_wrt_bra(2)*up_hbra_pints(2,addr_cur_hket,addr_cur_geo) &
                  + half_neg_rp*(ket_to_bra                                   &
                  * low_geo_zero(addr_hket_zero,addr_cur_geo)                 &
                  + low_hbra_pints(2,addr_hket_low,addr_up_geo_y)             &
                  - real(order_yket,REALK)*up_hbra_pints(2,addr_low_yket,addr_cur_geo))
                ! dxz
                up_hbra_pints(7,addr_cur_hket,addr_cur_geo)                    &
                  = cc_wrt_bra(3)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo)  &
                  + half_neg_rp*(low_hbra_pints(1,addr_hket_low,addr_up_geo_z) &
                  - real(order_zket,REALK)*up_hbra_pints(1,addr_low_zket,addr_cur_geo))
                ! dyz
                up_hbra_pints(8,addr_cur_hket,addr_cur_geo)                    &
                  = cc_wrt_bra(3)*up_hbra_pints(2,addr_cur_hket,addr_cur_geo)  &
                  + half_neg_rp*(low_hbra_pints(2,addr_hket_low,addr_up_geo_z) &
                  - real(order_zket,REALK)*up_hbra_pints(2,addr_low_zket,addr_cur_geo))
                ! dzz
                up_hbra_pints(9,addr_cur_hket,addr_cur_geo)                   &
                  = cc_wrt_bra(3)*up_hbra_pints(3,addr_cur_hket,addr_cur_geo) &
                  + half_neg_rp*(ket_to_bra                                   &
                  * low_geo_zero(addr_hket_zero,addr_cur_geo)                 &
                  + low_hbra_pints(3,addr_hket_low,addr_up_geo_z)             &
                  - real(order_zket,REALK)*up_hbra_pints(3,addr_low_zket,addr_cur_geo))
                if (max_cur_hbra>1) then
                  addr_up_hbra = 9   !base address of the f-shell HGTOs on bra center
                  addr_cur_hbra = 3  !base address of the d-shell HGTOs on bra center
                  addr_low_hbra = 0  !base address of the p-shell HGTOs on bra center
                  ! loops over other current order of HGTOs on bra center, starting from d-shell
                  do order_hbra = 2, max_cur_hbra
                    addr_up_hbra = addr_up_hbra+1
                    addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                    write(STDOUT,100) "orders:", order_hbra, order_xket, cur_order_geo
                    write(STDOUT,100) "x-direction"
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_low_hbra+1, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra, addr_low_xket, addr_up_geo
                    write(STDOUT,100) "------------------------------------"
#endif
                    ! x-direction
                    up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)        &
                      = cc_wrt_bra(1)                                             &
                      * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo)   &
                      + half_neg_rp*(real(order_hbra,REALK)*ket_to_bra            &
                      * up_hbra_pints(addr_low_hbra+1,addr_cur_hket,addr_cur_geo) &
                      + low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo)   &
                      - real(order_xket,REALK)                                    &
                      * up_hbra_pints(addr_cur_hbra,addr_low_xket,addr_cur_geo))
                    ! y-direction
                    addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                    write(STDOUT,100) "y-direction"
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+1
                    write(STDOUT,110) addr_cur_hbra, addr_low_yket, addr_cur_geo
                    write(STDOUT,100) "------------------------------------"
#endif
                    up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)         &
                      = cc_wrt_bra(2)                                              &
                      * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo)    &
                      + half_neg_rp                                                &
                      * (low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo+1) &
                      - real(order_yket,REALK)                                     &
                      * up_hbra_pints(addr_cur_hbra,addr_low_yket,addr_cur_geo))
                    do ibra = 1, order_hbra
                      addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                      write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                      write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_cur_geo
                      write(STDOUT,110) addr_low_hbra+ibra, addr_cur_hket, addr_cur_geo
                      write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_up_geo+1
                      write(STDOUT,110) addr_cur_hbra+ibra, addr_low_yket, addr_cur_geo
                      write(STDOUT,100) "------------------------------------"
#endif
                      up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo) &
                        = cc_wrt_bra(2)                                      &
                        * up_hbra_pints(addr_cur_hbra+ibra,addr_cur_hket,    &
                                        addr_cur_geo)                        &
                        + half_neg_rp*(real(ibra,REALK)*ket_to_bra           &
                        * up_hbra_pints(addr_low_hbra+ibra,addr_cur_hket,    &
                                        addr_cur_geo)                        &
                        + low_hbra_pints(addr_cur_hbra+ibra,addr_hket_low,   &
                                         addr_up_geo+1)                      &
                        - real(order_yket,REALK)                             &
                        * up_hbra_pints(addr_cur_hbra+ibra,addr_low_yket,addr_cur_geo))
                    end do
                    ! z-direction
#if defined(DEBUG)
                    write(STDOUT,100) "z-direction"
#endif
                    addr_cur_hbra = addr_cur_hbra-1
                    do jbra = 0, order_hbra
                      addr_up_hbra = addr_up_hbra+1
                      addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                      write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                      write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                      write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+igeo+2
                      write(STDOUT,110) addr_cur_hbra, addr_low_zket, addr_cur_geo
                      write(STDOUT,100) "------------------------------------"
#endif
                      up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)       &
                        = cc_wrt_bra(3)                                            &
                        * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo)  &
                        + half_neg_rp*(low_hbra_pints(addr_cur_hbra,addr_hket_low, &
                                                        addr_up_geo+igeo+2)        &
                        - real(order_zket,REALK)                                   &
                        * up_hbra_pints(addr_cur_hbra,addr_low_zket,addr_cur_geo))
                    end do
                    do ibra = 1, order_hbra
                      do jbra = 0, order_hbra-ibra
                        addr_up_hbra = addr_up_hbra+1
                        addr_cur_hbra = addr_cur_hbra+1
                        addr_low_hbra = addr_low_hbra+1
#if defined(DEBUG)
                        write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                        write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                        write(STDOUT,110) addr_low_hbra, addr_cur_hket, addr_cur_geo
                        write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+igeo+2
                        write(STDOUT,110) addr_cur_hbra, addr_low_zket, addr_cur_geo
                        write(STDOUT,100) "------------------------------------"
#endif
                        up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)      &
                          = cc_wrt_bra(3)                                           &
                          * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo) &
                          + half_neg_rp*(real(ibra,REALK)*ket_to_bra                &
                          * up_hbra_pints(addr_low_hbra,addr_cur_hket,addr_cur_geo) &
                          + low_hbra_pints(addr_cur_hbra,addr_hket_low,             &
                                           addr_up_geo+igeo+2)                      &
                          - real(order_zket,REALK)                                  &
                          * up_hbra_pints(addr_cur_hbra,addr_low_zket,addr_cur_geo))
                      end do
                    end do
                  end do
                end if
              end if
            end do
            ! (4.3) y...yz...z component of upper order HGTOs on ket center
            addr_cur_hket = addr_cur_hket+1
            addr_hket_zero = addr_hket_zero+1
            addr_yket_zero = addr_yket_zero+1
            addr_zket_zero = addr_zket_zero+1
            order_yket = order_hket-order_zket
            ! px HGTO on bra center
            up_hbra_pints(1,addr_cur_hket,addr_cur_geo)                 &
              = cc_wrt_bra(1)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + half_neg_rp*up_geo_zero(addr_cur_hket,addr_up_geo)
            ! py HGTO on bra center
            addr_up_geo_y = addr_up_geo+1
            up_hbra_pints(2,addr_cur_hket,addr_cur_geo)                 &
              = cc_wrt_bra(2)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + half_neg_rp*(up_geo_zero(addr_cur_hket,addr_up_geo_y)   &
              - real(order_yket,REALK)*low_geo_zero(addr_yket_zero,addr_cur_geo))
            ! pz HGTO on bra center
            addr_up_geo_z = addr_up_geo+igeo+2
            up_hbra_pints(3,addr_cur_hket,addr_cur_geo)                 &
              = cc_wrt_bra(3)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + half_neg_rp*(up_geo_zero(addr_cur_hket,addr_up_geo_z)   &
              - real(order_zket,REALK)*low_geo_zero(addr_zket_zero,addr_cur_geo))
            ! d-shell on bra center
            if (max_cur_hbra>0) then
              addr_low_yket = addr_low_yket+1
              addr_low_zket = addr_low_zket+1
              addr_hket_low = addr_hket_low+1
              ! dxx
              up_hbra_pints(4,addr_cur_hket,addr_cur_geo)                           &
                = cc_wrt_bra(1)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo)         &
                + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
                + low_hbra_pints(1,addr_hket_low,addr_up_geo))
              ! dxy
              up_hbra_pints(5,addr_cur_hket,addr_cur_geo)                    &
                = cc_wrt_bra(2)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo)  &
                + half_neg_rp*(low_hbra_pints(1,addr_hket_low,addr_up_geo_y) &
                - real(order_yket,REALK)*up_hbra_pints(1,addr_low_yket,addr_cur_geo))
              ! dyy
              up_hbra_pints(6,addr_cur_hket,addr_cur_geo)                           &
                = cc_wrt_bra(2)*up_hbra_pints(2,addr_cur_hket,addr_cur_geo)         &
                + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
                + low_hbra_pints(2,addr_hket_low,addr_up_geo_y)                     &
                - real(order_yket,REALK)*up_hbra_pints(2,addr_low_yket,addr_cur_geo))
              ! dxz
              up_hbra_pints(7,addr_cur_hket,addr_cur_geo)                    &
                = cc_wrt_bra(3)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo)  &
                + half_neg_rp*(low_hbra_pints(1,addr_hket_low,addr_up_geo_z) &
                - real(order_zket,REALK)*up_hbra_pints(1,addr_low_zket,addr_cur_geo))
              ! dyz
              up_hbra_pints(8,addr_cur_hket,addr_cur_geo)                    &
                = cc_wrt_bra(3)*up_hbra_pints(2,addr_cur_hket,addr_cur_geo)  &
                + half_neg_rp*(low_hbra_pints(2,addr_hket_low,addr_up_geo_z) &
                - real(order_zket,REALK)*up_hbra_pints(2,addr_low_zket,addr_cur_geo))
              ! dzz
              up_hbra_pints(9,addr_cur_hket,addr_cur_geo)                           &
                = cc_wrt_bra(3)*up_hbra_pints(3,addr_cur_hket,addr_cur_geo)         &
                + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
                + low_hbra_pints(3,addr_hket_low,addr_up_geo_z)                     &
                - real(order_zket,REALK)*up_hbra_pints(3,addr_low_zket,addr_cur_geo))
              if (max_cur_hbra>1) then
                addr_up_hbra = 9   !base address of the f-shell HGTOs on bra center
                addr_cur_hbra = 3  !base address of the d-shell HGTOs on bra center
                addr_low_hbra = 0  !base address of the p-shell HGTOs on bra center
                ! loops over other current order of HGTOs on bra center, starting from d-shell
                do order_hbra = 2, max_cur_hbra
                  addr_up_hbra = addr_up_hbra+1
                  addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                  write(STDOUT,100) "orders:", order_hbra, order_xket, cur_order_geo
                  write(STDOUT,100) "x-direction"
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_low_hbra+1, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,100) "------------------------------------"
#endif
                  ! x-direction
                  up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)        &
                    = cc_wrt_bra(1)                                             &
                    * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo)   &
                    + half_neg_rp*(real(order_hbra,REALK)*ket_to_bra            &
                    * up_hbra_pints(addr_low_hbra+1,addr_cur_hket,addr_cur_geo) &
                    + low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo))
                  ! y-direction
                  addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                  write(STDOUT,100) "y-direction"
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+1
                  write(STDOUT,110) addr_cur_hbra, addr_low_yket, addr_cur_geo
                  write(STDOUT,100) "------------------------------------"
#endif
                  up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)         &
                    = cc_wrt_bra(2)                                              &
                    * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo)    &
                    + half_neg_rp                                                &
                    * (low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo+1) &
                    - real(order_yket,REALK)                                     &
                    * up_hbra_pints(addr_cur_hbra,addr_low_yket,addr_cur_geo))
                  do ibra = 1, order_hbra
                    addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_low_hbra+ibra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_up_geo+1
                    write(STDOUT,110) addr_cur_hbra+ibra, addr_low_yket, addr_cur_geo
                    write(STDOUT,100) "------------------------------------"
#endif
                    up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)             &
                      = cc_wrt_bra(2)                                                  &
                      * up_hbra_pints(addr_cur_hbra+ibra,addr_cur_hket,addr_cur_geo)   &
                      + half_neg_rp*(real(ibra,REALK)*ket_to_bra                       &
                      * up_hbra_pints(addr_low_hbra+ibra,addr_cur_hket,addr_cur_geo)   &
                      + low_hbra_pints(addr_cur_hbra+ibra,addr_hket_low,addr_up_geo+1) &
                      - real(order_yket,REALK)                                         &
                      * up_hbra_pints(addr_cur_hbra+ibra,addr_low_yket,addr_cur_geo))
                  end do
                  ! z-direction
#if defined(DEBUG)
                  write(STDOUT,100) "z-direction"
#endif
                  addr_cur_hbra = addr_cur_hbra-1
                  do jbra = 0, order_hbra
                    addr_up_hbra = addr_up_hbra+1
                    addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+igeo+2
                    write(STDOUT,110) addr_cur_hbra, addr_low_zket, addr_cur_geo
                    write(STDOUT,100) "------------------------------------"
#endif
                    up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)      &
                      = cc_wrt_bra(3)                                           &
                      * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo) &
                      + half_neg_rp                                             &
                      * (low_hbra_pints(addr_cur_hbra,addr_hket_low,            &
                                        addr_up_geo+igeo+2)                     &
                      - real(order_zket,REALK)                                  &
                      * up_hbra_pints(addr_cur_hbra,addr_low_zket,addr_cur_geo))
                  end do
                  do ibra = 1, order_hbra
                    do jbra = 0, order_hbra-ibra
                      addr_up_hbra = addr_up_hbra+1
                      addr_cur_hbra = addr_cur_hbra+1
                      addr_low_hbra = addr_low_hbra+1
#if defined(DEBUG)
                      write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                      write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                      write(STDOUT,110) addr_low_hbra, addr_cur_hket, addr_cur_geo
                      write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+igeo+2
                      write(STDOUT,110) addr_cur_hbra, addr_low_zket, addr_cur_geo
                      write(STDOUT,100) "------------------------------------"
#endif
                      up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)      &
                        = cc_wrt_bra(3)                                           &
                        * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo) &
                        + half_neg_rp*(real(ibra,REALK)*ket_to_bra                &
                        * up_hbra_pints(addr_low_hbra,addr_cur_hket,addr_cur_geo) &
                        + low_hbra_pints(addr_cur_hbra,addr_hket_low,             &
                                         addr_up_geo+igeo+2)                      &
                        - real(order_zket,REALK)                                  &
                        * up_hbra_pints(addr_cur_hbra,addr_low_zket,addr_cur_geo))
                    end do
                  end do
                end do
              end if
            end if
          end do
          ! (5) z...z component of upper order HGTOs on ket center
          addr_cur_hket = addr_cur_hket+1
          addr_hket_zero = addr_hket_zero+1
          addr_zket_zero = addr_zket_zero+1
          ! px HGTO on bra center
          up_hbra_pints(1,addr_cur_hket,addr_cur_geo)                 &
            = cc_wrt_bra(1)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
            + half_neg_rp*up_geo_zero(addr_cur_hket,addr_up_geo)
          ! py HGTO on bra center
          addr_up_geo_y = addr_up_geo+1
          up_hbra_pints(2,addr_cur_hket,addr_cur_geo)                 &
            = cc_wrt_bra(2)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
            + half_neg_rp*up_geo_zero(addr_cur_hket,addr_up_geo_y)
          ! pz HGTO on bra center
          addr_up_geo_z = addr_up_geo+igeo+2
          up_hbra_pints(3,addr_cur_hket,addr_cur_geo)                 &
            = cc_wrt_bra(3)*low_geo_zero(addr_hket_zero,addr_cur_geo) &
            + half_neg_rp*(up_geo_zero(addr_cur_hket,addr_up_geo_z)   &
            - real(order_hket,REALK)*low_geo_zero(addr_zket_zero,addr_cur_geo))
          ! d-shell on bra center
          if (max_cur_hbra>0) then
            addr_low_zket = addr_low_zket+1
            addr_hket_low = addr_hket_low+1
            ! dxx
            up_hbra_pints(4,addr_cur_hket,addr_cur_geo)                           &
              = cc_wrt_bra(1)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo)         &
              + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + low_hbra_pints(1,addr_hket_low,addr_up_geo))
            ! dxy
            up_hbra_pints(5,addr_cur_hket,addr_cur_geo)                   &
              = cc_wrt_bra(2)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo) &
              + half_neg_rp*low_hbra_pints(1,addr_hket_low,addr_up_geo_y)
            ! dyy
            up_hbra_pints(6,addr_cur_hket,addr_cur_geo)                           &
              = cc_wrt_bra(2)*up_hbra_pints(2,addr_cur_hket,addr_cur_geo)         &
              + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + low_hbra_pints(2,addr_hket_low,addr_up_geo_y))
            ! dxz
            up_hbra_pints(7,addr_cur_hket,addr_cur_geo)                    &
              = cc_wrt_bra(3)*up_hbra_pints(1,addr_cur_hket,addr_cur_geo)  &
              + half_neg_rp*(low_hbra_pints(1,addr_hket_low,addr_up_geo_z) &
              - real(order_hket,REALK)*up_hbra_pints(1,addr_low_zket,addr_cur_geo))
            ! dyz
            up_hbra_pints(8,addr_cur_hket,addr_cur_geo)                    &
              = cc_wrt_bra(3)*up_hbra_pints(2,addr_cur_hket,addr_cur_geo)  &
              + half_neg_rp*(low_hbra_pints(2,addr_hket_low,addr_up_geo_z) &
              - real(order_hket,REALK)*up_hbra_pints(2,addr_low_zket,addr_cur_geo))
            ! dzz
            up_hbra_pints(9,addr_cur_hket,addr_cur_geo)                           &
              = cc_wrt_bra(3)*up_hbra_pints(3,addr_cur_hket,addr_cur_geo)         &
              + half_neg_rp*(ket_to_bra*low_geo_zero(addr_hket_zero,addr_cur_geo) &
              + low_hbra_pints(3,addr_hket_low,addr_up_geo_z)                     &
              - real(order_hket,REALK)*up_hbra_pints(3,addr_low_zket,addr_cur_geo))
            if (max_cur_hbra>1) then
              addr_up_hbra = 9   !base address of the f-shell HGTOs on bra center
              addr_cur_hbra = 3  !base address of the d-shell HGTOs on bra center
              addr_low_hbra = 0  !base address of the p-shell HGTOs on bra center
              ! loops over other current order of HGTOs on bra center, starting from d-shell
              do order_hbra = 2, max_cur_hbra
                addr_up_hbra = addr_up_hbra+1
                addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                write(STDOUT,100) "orders:", order_hbra, order_xket, cur_order_geo
                write(STDOUT,100) "x-direction"
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_low_hbra+1, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                write(STDOUT,100) "------------------------------------"
#endif
                ! x-direction
                up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)        &
                  = cc_wrt_bra(1)                                             &
                  * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo)   &
                  + half_neg_rp*(real(order_hbra,REALK)*ket_to_bra            &
                  * up_hbra_pints(addr_low_hbra+1,addr_cur_hket,addr_cur_geo) &
                  + low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo))
                ! y-direction
                addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                write(STDOUT,100) "y-direction"
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+1
                write(STDOUT,100) "------------------------------------"
#endif
                up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)      &
                  = cc_wrt_bra(2)                                           &
                  * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo) &
                  + half_neg_rp                                             &
                  * (low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo+1))
                do ibra = 1, order_hbra
                  addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_low_hbra+ibra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_up_geo+1
                  write(STDOUT,100) "------------------------------------"
#endif
                  up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)           &
                    = cc_wrt_bra(2)                                                &
                    * up_hbra_pints(addr_cur_hbra+ibra,addr_cur_hket,addr_cur_geo) &
                    + half_neg_rp*(real(ibra,REALK)*ket_to_bra                     &
                    * up_hbra_pints(addr_low_hbra+ibra,addr_cur_hket,addr_cur_geo) &
                    + low_hbra_pints(addr_cur_hbra+ibra,addr_hket_low,addr_up_geo+1))
                end do
                ! z-direction
#if defined(DEBUG)
                write(STDOUT,100) "z-direction"
#endif
                addr_cur_hbra = addr_cur_hbra-1
                do jbra = 0, order_hbra
                  addr_up_hbra = addr_up_hbra+1
                  addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+igeo+2
                  write(STDOUT,110) addr_cur_hbra, addr_low_zket, addr_cur_geo
                  write(STDOUT,100) "------------------------------------"
#endif
                  up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)              &
                    = cc_wrt_bra(3)                                                   &
                    * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo)         &
                    + half_neg_rp                                                     &
                    * (low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo+igeo+2) &
                    - real(order_hket,REALK)                                          &
                    * up_hbra_pints(addr_cur_hbra,addr_low_zket,addr_cur_geo))
                end do
                do ibra = 1, order_hbra
                  do jbra = 0, order_hbra-ibra
                    addr_up_hbra = addr_up_hbra+1
                    addr_cur_hbra = addr_cur_hbra+1
                    addr_low_hbra = addr_low_hbra+1
#if defined(DEBUG)
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_low_hbra, addr_cur_hket, addr_cur_geo
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_geo+igeo+2
                    write(STDOUT,110) addr_cur_hbra, addr_low_zket, addr_cur_geo
                    write(STDOUT,100) "------------------------------------"
#endif
                    up_hbra_pints(addr_up_hbra,addr_cur_hket,addr_cur_geo)             &
                      = cc_wrt_bra(3)                                                  &
                      * up_hbra_pints(addr_cur_hbra,addr_cur_hket,addr_cur_geo)        &
                      + half_neg_rp*(real(ibra,REALK)*ket_to_bra                       &
                      * up_hbra_pints(addr_low_hbra,addr_cur_hket,addr_cur_geo)        &
                      + low_hbra_pints(addr_cur_hbra,addr_hket_low,addr_up_geo+igeo+2) &
                      - real(order_hket,REALK)                                         &
                      * up_hbra_pints(addr_cur_hbra,addr_low_zket,addr_cur_geo))
                  end do
                end do
              end do
            end if
          end if
          ! updates the base addresses of lower order HGTOs on ket center
          base_low_hket = addr_low_zket
          base_hket_zero = addr_zket_zero
        end do
      end do
      addr_up_geo = addr_up_geo+1
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "sub_nucpot_hbra", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("sub_nucpot_hbra>> ",A,3I6)
110 format("sub_nucpot_hbra>> ","HBRA",I8,4X,"HKET",I8,4X,"DERV",I8)
#endif
  end subroutine sub_nucpot_hbra

  !> \brief assigns the integrals with required HGTOs on bra center
  !> \author Bin Gao
  !> \date 2012-03-04
  !> \param offset_hgto_ket is the offset of HGTOs on ket center
  !> \param offset_hket_geo is the offset of geometric derivatives on nuclear potential origin
  !> \param dim_hket is the dimension of HGTOs on ket center
  !> \param dim_geo_hket is the dimension of geometric derivatives on nuclear potential origin
  !> \param hket_pints contains the nuclear attraction integrals with zeroth order
  !>        HGTOs on bra center
  !> \param dim_up_hbra is the dimension of upper order HGTOs on bra center in temporary integrals
  !> \param dim_cur_hket is the dimension of current order HGTOs on ket center in temporary integrals
  !> \param num_low_geo is the number of lower order geometric derivatives in temporary integrals
  !> \param recur_pints contains the temporary integrals from recurrence relations
  !> \param zero_hbra indicates if zeroth order HGTO returned
  !> \param offset_hbra_geo is the offset of geometric derivatives on nuclear potential origin
  !>        in returned integrals
  !> \param dim_hgto_bra is the dimension of HGTOs of bra center afterwards
  !> \param dim_hgto_ket is the dimension of HGTOs of ket center afterwards
  !> \param dim_geo_hbra is the dimension of geometric derivatives afterwards
  !> \return hbra_pints contains the integrals with required HGTOs on bra center
  subroutine nucpot_hbra_assign(offset_hgto_ket, offset_hket_geo, dim_hket, &
                                dim_geo_hket, hket_pints, dim_up_hbra,      &
                                dim_cur_hket, num_low_geo, recur_pints,     &
                                zero_hbra, offset_hbra_geo, dim_hgto_bra,   &
                                dim_hgto_ket, dim_geo_hbra, hbra_pints)
    use xkind
    implicit none
    integer, intent(in) :: offset_hgto_ket
    integer, intent(in) :: offset_hket_geo
    integer, intent(in) :: dim_hket
    integer, intent(in) :: dim_geo_hket
    real(REALK), intent(in) :: hket_pints(dim_hket,dim_geo_hket)
    integer, intent(in) :: dim_up_hbra
    integer, intent(in) :: dim_cur_hket
    integer, intent(in) :: num_low_geo
    real(REALK), intent(in) :: recur_pints(dim_up_hbra,dim_cur_hket,num_low_geo)
    logical, intent(in) :: zero_hbra
    integer, intent(in) :: offset_hbra_geo
    integer, intent(in) :: dim_hgto_bra
    integer, intent(in) :: dim_hgto_ket
    integer, intent(in) :: dim_geo_hbra
    real(REALK), intent(inout) :: hbra_pints(dim_hgto_bra,dim_hgto_ket,dim_geo_hbra)
!f2py intent(in) :: offset_hgto_ket
!f2py intent(in) :: offset_hket_geo
!f2py intent(hide) :: dim_hket
!f2py intent(hide) :: dim_geo_hket
!f2py intent(in) :: hket_pints
!f2py intent(hide) :: dim_up_hbra
!f2py intent(hide) :: dim_cur_hket
!f2py intent(hide) :: num_low_geo
!f2py intent(in) :: recur_pints
!f2py intent(in) :: zero_hbra
!f2py intent(in) :: offset_hbra_geo
!f2py intent(hide) :: dim_hgto_bra
!f2py intent(hide) :: dim_hgto_ket
!f2py intent(hide) :: dim_geo_hbra
!f2py intent(inout) :: hbra_pints
    integer offset_cur_hket  !offset of current order HGTOs on ket center in temporary integrals
    integer start_up_hbra    !start address of upper order HGTOs on bra center in temporary integrals
    integer addr_hket_geo    !address of geometric derivatives on nuclear potential origin
    integer addr_hbra_geo    !address of geometric derivatives on nuclear potential origin in returned integrals
    integer igeo             !incremental recorder over geometric derivatives
    integer iket             !incremental recorder over HGTOs on ket center
#if defined(XTIME)
    real(REALK) curr_time    !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! sets the offset of current order HGTOs on ket center in temporary integrals
    offset_cur_hket = dim_cur_hket-dim_hgto_ket
    ! s-shell HGTO returned
    if (zero_hbra) then
      start_up_hbra = dim_up_hbra-dim_hgto_bra+2
      do igeo = 1, num_low_geo
        addr_hket_geo = offset_hket_geo+igeo
        addr_hbra_geo = offset_hbra_geo+igeo
        do iket = 1, dim_hgto_ket
          hbra_pints(1,iket,addr_hbra_geo) &
            = hket_pints(offset_hgto_ket+iket,addr_hket_geo)
          hbra_pints(2:dim_hgto_bra,iket,addr_hbra_geo) &
            = recur_pints(start_up_hbra:dim_up_hbra,offset_cur_hket+iket,igeo)
        end do
      end do
    else
      start_up_hbra = dim_up_hbra-dim_hgto_bra+1
      do igeo = 1, num_low_geo
        addr_hbra_geo = offset_hbra_geo+igeo
        do iket = 1, dim_hgto_ket
          hbra_pints(:,iket,addr_hbra_geo) &
            = recur_pints(start_up_hbra:dim_up_hbra,offset_cur_hket+iket,igeo)
        end do
      end do
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "nucpot_hbra_assign", STDOUT)
#endif
    return
  end subroutine nucpot_hbra_assign
