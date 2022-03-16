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
!!  This file contains the horizontal recurrence relations of recovering
!!  the HGTOs on ket center from integrals from \fn(carmom_hbra).
!!
!!  2012-03-04, Bin Gao
!!  * rewrites to improve efficiency
!!
!!  2012-02-12, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief horizontal recurrence relations (HRR) of recovering the HGTOs
  !>        on ket center from integrals from \fn(carmom_hbra)
  !> \author Bin Gao
  !> \date 2012-02-12
  !> \param orders_hgto_bra contains the range of orders of HGTOs on ket center to return
  !> \param exponent_bra is the exponent of primitive Gaussians of bra center
  !> \param orders_hgto_ket contains the range of orders of HGTOs on bra center to return
  !> \param exponent_ket is the exponent of primitive Gaussians of ket center
  !> \param dim_hbra is the dimension of HGTOs on bra center from \fn(carmom_hbra)
  !> \param hbra_pints contains the primitive Hermite integrals from \fn(carmom_hbra)
  !> \param dim_hgto_bra is the dimension of HGTOs on bra center
  !> \param dim_hgto_ket is the dimension of HGTOs on ket center
  !> \return hket_pints contains the integrals with specific orders of HGTOs on bra and ket centers
  subroutine carmom_hrr_ket(orders_hgto_bra, exponent_bra, &
                            orders_hgto_ket, exponent_ket, &
                            dim_hbra, hbra_pints,          &
                            dim_hgto_bra, dim_hgto_ket, hket_pints)
    use xkind
    implicit none
    integer, intent(in) :: orders_hgto_bra(2)
    real(REALK), intent(in) :: exponent_bra
    integer, intent(in) :: orders_hgto_ket(2)
    real(REALK), intent(in) :: exponent_ket
    integer, intent(in) :: dim_hbra
    real(REALK), intent(in) :: hbra_pints(dim_hbra)
    integer, intent(in) :: dim_hgto_bra
    integer, intent(in) :: dim_hgto_ket
    real(REALK), intent(out) :: hket_pints(dim_hgto_bra,dim_hgto_ket)
!f2py intent(in) :: orders_hgto_bra
!f2py intent(in) :: exponent_bra
!f2py intent(in) :: orders_hgto_ket
!f2py intent(in) :: exponent_ket
!f2py intent(hide) :: dim_hbra
!f2py intent(in) :: hbra_pints
!f2py intent(in) :: dim_hgto_bra
!f2py intent(in) :: dim_hgto_ket
!f2py intent(out) :: hket_pints
!f2py depend(dim_hgto_bra) :: hket_pints
!f2py depend(dim_hgto_ket) :: hket_pints
    real(REALK) neg_ratio_braket  !negative ratio between exponents of bra and ket centers
    integer offset_up_hbra        !offset of upper order HGTOs on bra center in recurrence relations
    integer low_order_hket        !lower order HGTOs on ket center
    integer low_orders_hbra(2)    !range of lower orders of HGTOs on bra center
    integer dim_up_hbra           !dimension of upper order HGTOs on bra center
    integer num_low_hket          !number of lower order HGTOs on ket center
    integer low_hket_int          !pointer to temporary integrals of lower order HGTOs on ket center
    integer dim_low_hbra          !dimension of lower order HGTOs on bra center
    integer num_up_hket           !number of upper order HGTOs on ket center
    integer up_hket_int           !pointer of upper order HGTOs on ket center
    integer dim_tmp               !dimension of temporary integrals
    real(REALK), allocatable :: tmp_ints(:,:)
                                  !temporary integrals
    integer offset_hket           !offset of returned HGTOs on ket center
    integer ierr                  !error information
#if defined(XTIME)
    real(REALK) curr_time         !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    select case(orders_hgto_ket(2))
    ! only zeroth order of HGTO on ket center returned
    case(0)
      hket_pints(:,1) = hbra_pints(1:dim_hgto_bra)
    ! the maximum order of HGTOs on ket center returned is the first
    case(1)
      neg_ratio_braket = -exponent_bra/exponent_ket
      select case(orders_hgto_ket(1))
      ! the minimum order of HGTOs on ket center returned is zeroth
      case(0)
        ! assigns the zeroth order HGTO on ket center
        hket_pints(:,1) = hbra_pints(1:dim_hgto_bra)
        ! the first few integrals do not need in the following recurrence relations
        offset_up_hbra = (orders_hgto_bra(1)+1)*(orders_hgto_bra(1)+2)/2
        call sub_carmom_hrr_ket(orders_hgto_bra, 0, neg_ratio_braket,    &
                                dim_hbra, 1, offset_up_hbra, hbra_pints, &
                                dim_hgto_bra, 3, hket_pints(:,2:4))
      ! only first order of HGTOs on ket center returned
      case default
        call sub_carmom_hrr_ket(orders_hgto_bra, 0, neg_ratio_braket, &
                                dim_hbra, 1, 0, hbra_pints,           &
                                dim_hgto_bra, 3, hket_pints(:,1:3))
      end select
    ! the maximum order of HGTOs on ket center returned is, at least, the second
    case default
      neg_ratio_braket = -exponent_bra/exponent_ket
      ! sets the offset of upper order HGTOs on bra center in recurrence relations
      offset_up_hbra = (orders_hgto_bra(1)+1)*(orders_hgto_bra(1)+2)/2
      ! initializes the range of lower orders of HGTOs on bra center
      low_orders_hbra(1) = orders_hgto_bra(1)+orders_hgto_ket(1)
      low_orders_hbra(2) = orders_hgto_bra(2)+orders_hgto_ket(2)
      ! allocates memory for temporary integrals
      dim_tmp = dim_hbra
      call dim_carmom_hrr_ket(orders_hgto_bra, orders_hgto_ket, dim_tmp)
      allocate(tmp_ints(dim_tmp,2), stat=ierr)
      if (ierr/=0)                        &
        call error_stop("carmom_hrr_ket", &
                        "failed to allocate tmp_ints", dim_tmp*2)
      select case(orders_hgto_ket(1))
      ! the minimum order of HGTOs on ket center returned is zeroth
      case(0)
        ! sets the dimension of low order HGTOs on bra center
        dim_low_hbra = dim_hbra-(low_orders_hbra(2)+1)*(low_orders_hbra(2)+2)/2
        ! updates the maximum lower order of HGTOs on bra center, \var(low_orders_hbra(1)) does not change
        low_orders_hbra(2) = low_orders_hbra(2)-1
        ! assigns the zeroth order HGTO on ket center
        hket_pints(:,1) = hbra_pints(1:dim_hgto_bra)
        ! gets the first order HGTOs on ket center
        call sub_carmom_hrr_ket(low_orders_hbra, 0, neg_ratio_braket,    &
                                dim_hbra, 1, offset_up_hbra, hbra_pints, &
                                dim_low_hbra, 3, tmp_ints(:,1))
        ! assigns the first order HGTO on ket center
        call carmom_hrr_ket_assign(dim_low_hbra, 3, tmp_ints(:,1), 1, &
                                   dim_hgto_bra, dim_hgto_ket, hket_pints)
        ! initializes the number of upper order HGTOs on ket center
        num_up_hket = 3
        ! initializes the pointers of HGTOs on ket center
        low_hket_int = 2
        up_hket_int = 1
        ! sets the offset of returned HGTOs on ket center
        offset_hket = 4
        ! other returned orders with constant \var(low_orders_hbra(1))
        do low_order_hket = 1, orders_hgto_ket(2)-1
          ! updates the dimensions of HGTOs on bra center
          dim_up_hbra = dim_low_hbra
          dim_low_hbra = dim_low_hbra-(low_orders_hbra(2)+1)*(low_orders_hbra(2)+2)/2
          ! updates the maximum lower order of HGTOs on bra center
          low_orders_hbra(2) = low_orders_hbra(2)-1
          ! updates the number of HGTOs on ket center
          num_low_hket = num_up_hket
          num_up_hket = num_low_hket+low_order_hket+2
          ! switches the pointers
          low_hket_int = 3-low_hket_int
          up_hket_int = 3-up_hket_int
          ! gets the temporary integrals
          call sub_carmom_hrr_ket(low_orders_hbra, low_order_hket, &
                                  neg_ratio_braket, dim_up_hbra,   &
                                  num_low_hket, offset_up_hbra,    &
                                  tmp_ints(:,low_hket_int),        &
                                  dim_low_hbra, num_up_hket,       &
                                  tmp_ints(:,up_hket_int))
          ! assigns the returned integrals
          call carmom_hrr_ket_assign(dim_low_hbra, num_up_hket, &
                                     tmp_ints(:,up_hket_int),   &
                                     offset_hket, dim_hgto_bra, &
                                     dim_hgto_ket, hket_pints)
          ! updates the offset of returned HGTOs on ket center
          offset_hket = offset_hket+num_up_hket
        end do
      ! the minimum order of HGTOs on ket center returned is first
      case(1)
        ! sets the dimension of low order HGTOs on bra center
        dim_low_hbra = dim_hbra                                   &
                     + (low_orders_hbra(1)*(low_orders_hbra(1)+1) &
                     -  (low_orders_hbra(2)+1)*(low_orders_hbra(2)+2))/2
        ! updates the range of lower orders of HGTOs on bra center
        low_orders_hbra = low_orders_hbra-1
        ! gets the first order HGTOs on ket center
        call sub_carmom_hrr_ket(low_orders_hbra, 0, neg_ratio_braket, &
                                dim_hbra, 1, 0, hbra_pints,           &
                                dim_low_hbra, 3, tmp_ints(:,1))
        ! assigns the first order HGTO on ket center
        call carmom_hrr_ket_assign(dim_low_hbra, 3, tmp_ints(:,1), 0, &
                                   dim_hgto_bra, dim_hgto_ket, hket_pints)
        ! initializes the pointers of HGTOs on ket center
        low_hket_int = 1
        up_hket_int = 2
        ! initializes the number of upper order HGTOs on ket center
        num_up_hket = 3
        ! sets the offset of returned HGTOs on ket center
        offset_hket = 3
        ! other returned orders with constant \var(low_orders_hbra(1))
        do low_order_hket = 1, orders_hgto_ket(2)-1
          ! updates the dimensions of HGTOs on bra center
          dim_up_hbra = dim_low_hbra
          dim_low_hbra = dim_low_hbra-(low_orders_hbra(2)+1)*(low_orders_hbra(2)+2)/2
          ! updates the maximum lower order of HGTOs on bra center
          low_orders_hbra(2) = low_orders_hbra(2)-1
          ! updates the number of HGTOs on ket center
          num_low_hket = num_up_hket
          num_up_hket = num_low_hket+low_order_hket+2
          call sub_carmom_hrr_ket(low_orders_hbra, low_order_hket, &
                                  neg_ratio_braket, dim_up_hbra,   &
                                  num_low_hket, offset_up_hbra,    &
                                  tmp_ints(:,low_hket_int),        &
                                  dim_low_hbra, num_up_hket,       &
                                  tmp_ints(:,up_hket_int))
          ! assigns the returned integrals
          call carmom_hrr_ket_assign(dim_low_hbra, num_up_hket, &
                                     tmp_ints(:,up_hket_int),   &
                                     offset_hket, dim_hgto_bra, &
                                     dim_hgto_ket, hket_pints)
          ! updates the offset of returned HGTOs on ket center
          offset_hket = offset_hket+num_up_hket
          ! switches the pointers
          low_hket_int = 3-low_hket_int
          up_hket_int = 3-up_hket_int
        end do
      ! the minimum order of HGTOs on ket center returned is, at least, the second
      case default
        ! sets the dimension of low order HGTOs on bra center
        dim_low_hbra = dim_hbra                                   &
                     + (low_orders_hbra(1)*(low_orders_hbra(1)+1) &
                     -  (low_orders_hbra(2)+1)*(low_orders_hbra(2)+2))/2
        ! updates the range of lower orders of HGTOs on bra center
        low_orders_hbra = low_orders_hbra-1
        ! gets the first order HGTOs on ket center
        call sub_carmom_hrr_ket(low_orders_hbra, 0, neg_ratio_braket, &
                                dim_hbra, 1, 0, hbra_pints,           &
                                dim_low_hbra, 3, tmp_ints(:,1))
        ! initializes the number of upper order HGTOs on ket center
        num_up_hket = 3
        ! initializes the pointers of HGTOs on ket center
        low_hket_int = 2
        up_hket_int = 1
        ! orders of HGTOs on ket center which are not returned, with decreased \var(low_orders_hbra)
        do low_order_hket = 1, orders_hgto_ket(1)-1
          ! updates the dimensions of HGTOs on bra center
          dim_up_hbra = dim_low_hbra
          dim_low_hbra = dim_low_hbra                               &
                       + (low_orders_hbra(1)*(low_orders_hbra(1)+1) &
                       -  (low_orders_hbra(2)+1)*(low_orders_hbra(2)+2))/2
          ! updates the minimum and maximum lower orders of HGTOs on bra center
          low_orders_hbra = low_orders_hbra-1
          ! updates the number of HGTOs on ket center
          num_low_hket = num_up_hket
          num_up_hket = num_low_hket+low_order_hket+2
          ! switches the pointers
          low_hket_int = 3-low_hket_int
          up_hket_int = 3-up_hket_int
          ! gets the temporary integrals
          call sub_carmom_hrr_ket(low_orders_hbra, low_order_hket,           &
                                  neg_ratio_braket, dim_up_hbra,             &
                                  num_low_hket, 0, tmp_ints(:,low_hket_int), &
                                  dim_low_hbra, num_up_hket, tmp_ints(:,up_hket_int))
        end do
        ! assigns the returned integrals
        call carmom_hrr_ket_assign(dim_low_hbra, num_up_hket,  &
                                   tmp_ints(:,up_hket_int), 0, &
                                   dim_hgto_bra, dim_hgto_ket, hket_pints)
        ! sets the offset of returned HGTOs on ket center
        offset_hket = num_up_hket
        ! other returned orders with constant \var(low_orders_hbra(1))
        do low_order_hket = orders_hgto_ket(1), orders_hgto_ket(2)-1
          ! updates the dimensions of HGTOs on bra center
          dim_up_hbra = dim_low_hbra
          dim_low_hbra = dim_low_hbra-(low_orders_hbra(2)+1)*(low_orders_hbra(2)+2)/2
          ! updates the maximum lower order of HGTOs on bra center
          low_orders_hbra(2) = low_orders_hbra(2)-1
          ! updates the number of HGTOs on ket center
          num_low_hket = num_up_hket
          num_up_hket = num_low_hket+low_order_hket+2
          ! switches the pointers
          low_hket_int = 3-low_hket_int
          up_hket_int = 3-up_hket_int
          ! gets the temporary integrals
          call sub_carmom_hrr_ket(low_orders_hbra, low_order_hket, &
                                  neg_ratio_braket, dim_up_hbra,   &
                                  num_low_hket, offset_up_hbra,    &
                                  tmp_ints(:,low_hket_int),        &
                                  dim_low_hbra, num_up_hket,       &
                                  tmp_ints(:,up_hket_int))
          ! assigns the returned integrals
          call carmom_hrr_ket_assign(dim_low_hbra, num_up_hket, &
                                     tmp_ints(:,up_hket_int),   &
                                     offset_hket, dim_hgto_bra, &
                                     dim_hgto_ket, hket_pints)
          ! updates the offset of returned HGTOs on ket center
          offset_hket = offset_hket+num_up_hket
        end do
      end select
      deallocate(tmp_ints)
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "carmom_hrr_ket", STDOUT)
#endif
    return
  end subroutine carmom_hrr_ket

  !> \brief gets the maximum dimension of temporary integrals used in recurrence relations
  !> \author Bin Gao
  !> \date 2012-03-05
  !> \param orders_hgto_bra contains the range of orders of HGTOs on ket center to return
  !> \param orders_hgto_ket contains the range of orders of HGTOs on bra center to return
  !> \return dim_ints is the maximum dimension of temporary integrals
  subroutine dim_carmom_hrr_ket(orders_hgto_bra, orders_hgto_ket, dim_ints)
    use xkind
    implicit none
    integer, intent(in) :: orders_hgto_bra(2)
    integer, intent(in) :: orders_hgto_ket(2)
    integer, intent(inout) :: dim_ints
!f2py intent(in) :: orders_hgto_bra
!f2py intent(in) :: orders_hgto_ket
!f2py intent(inout) :: dim_ints
    integer range_hgto_bra(2)  !range of orders of HGTOs on bra center
    integer dim_hgto_bra       !dimension of HGTOs on bra center
    integer num_hgto_ket       !number of HGTOs on ket center
    integer order_hket         !incremental recorder over orders of HGTOs on ket center
    integer dim_tmp            !temporary result of dimension
#if defined(XTIME)
    real(REALK) curr_time      !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! initializes the range of orders of HGTOs on bra center
    range_hgto_bra(1) = orders_hgto_bra(1)+orders_hgto_ket(1)
    range_hgto_bra(2) = orders_hgto_bra(2)+orders_hgto_ket(2)
    ! initializes the dimension of HGTOs on bra center
    dim_hgto_bra = dim_ints
    ! resets the return value
    dim_ints = 0
    ! initializes the number of HGTOs on ket center
    num_hgto_ket = 1
    ! loops over the orders of HGTOs on ket center till to the minimum returned order
    ! the minimum and maximum of orders of HGTOs on bra center decrease
    do order_hket = 1, orders_hgto_ket(1)
      ! updates the dimension of HGTOs on bra center
      dim_hgto_bra = dim_hgto_bra                             &
                   + (range_hgto_bra(1)*(range_hgto_bra(1)+1) & 
                   -  (range_hgto_bra(2)+1)*(range_hgto_bra(2)+2))/2
      ! updates the minimum and maximum of orders of HGTOs on bra center
      range_hgto_bra = range_hgto_bra-1
      ! updates the number of HGTOs on ket center
      num_hgto_ket = num_hgto_ket+order_hket+1
      ! updates the maximum dimension
      dim_tmp = dim_hgto_bra*num_hgto_ket
      if (dim_tmp>dim_ints) dim_ints = dim_tmp
    end do
    ! other returned orders with constant \var(range_hgto_bra(1))
    do order_hket = orders_hgto_ket(1)+1, orders_hgto_ket(2)
      ! updates the dimension of HGTOs on bra center
      dim_hgto_bra = dim_hgto_bra-(range_hgto_bra(2)+1)*(range_hgto_bra(2)+2)/2
      ! updates the maximum of orders of HGTOs on bra center
      range_hgto_bra(2) = range_hgto_bra(2)-1
      ! updates the number of HGTOs on ket center
      num_hgto_ket = num_hgto_ket+order_hket+1
      ! updates the maximum dimension
      dim_tmp = dim_hgto_bra*num_hgto_ket
      if (dim_tmp>dim_ints) dim_ints = dim_tmp
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "dim_carmom_hrr_ket", STDOUT)
#endif
    return
  end subroutine dim_carmom_hrr_ket

  !> \brief sub-recurrence relations by recovering HGTOs on ket center from
  !>        those on bra center
  !> \author Bin Gao
  !> \date 2012-02-12
  !> \param low_orders_hbra contains the range of orders of HGTOs on bra center
  !> \param low_order_hket is the order of HGTOs on ket center
  !> \param neg_ratio_braket is the negative ratio between exponents of bra and ket centers
  !> \param dim_up_hbra is the dimension of upper order HGTOs on bra center
  !> \param num_low_hket is the number of lower order HGTOs on ket center
  !> \param offset_up_hbra is the offset of upper order HGTOs on bra center in recurrence relations
  !> \param low_hket_pints contains the integrals with upper order HGTOs on bra center,
  !>        and lower order HGTOs on ket center
  !> \param dim_low_hbra is the dimension of lower order HGTOs on bra center
  !> \param num_up_hket is the number of upper order HGTOs on ket center
  !> \return up_hket_pints contains the integrals with lower order HGTOs on bra center,
  !>         and upper order HGTOs on ket center
  subroutine sub_carmom_hrr_ket(low_orders_hbra, low_order_hket, &
                                neg_ratio_braket, dim_up_hbra,   &
                                num_low_hket, offset_up_hbra,    &
                                low_hket_pints, dim_low_hbra,    &
                                num_up_hket, up_hket_pints)
    use xkind
    implicit none
    integer, intent(in) :: low_orders_hbra(2)
    integer, intent(in) :: low_order_hket
    real(REALK), intent(in) :: neg_ratio_braket
    integer, intent(in) :: dim_up_hbra
    integer, intent(in) :: num_low_hket
    integer, intent(in) :: offset_up_hbra
    real(REALK), intent(in) :: low_hket_pints(dim_up_hbra,num_low_hket)
    integer, intent(in) :: dim_low_hbra
    integer, intent(in) :: num_up_hket
    real(REALK), intent(out) :: up_hket_pints(dim_low_hbra,num_up_hket)
!f2py intent(in) :: low_orders_hbra
!f2py intent(in) :: low_order_hket
!f2py intent(in) :: neg_ratio_braket
!f2py intent(hide) :: dim_up_hbra
!f2py intent(hide) :: num_low_hket
!f2py intent(in) :: offset_up_hbra
!f2py intent(in) :: low_hket_pints
!f2py intent(in) :: dim_low_hbra
!f2py intent(in) :: num_up_hket
!f2py intent(out) :: up_hket_pints
!f2py depend(dim_low_hbra) :: up_hket_pints
!f2py depend(num_up_hket) :: up_hket_pints
    integer addr_up_hbra   !address of upper order HGTOs on bra center
    integer addr_low_hbra  !address of lower order HGTOs on bra center
    integer addr_up_hket   !address of upper order HGTOs on ket center
    integer addr_low_hket  !address of lower order HGTOs on ket center
    integer order_hbra     !incremental recorder over orders of HGTOs on bra center
    integer ihbra, jhbra   !incremental recorders over xyz components of HGTOs on bra center
    integer ihket, jhket   !incremental recorders over xyz components of HGTOs on ket center
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    addr_up_hket = 1
    addr_low_hket = 1
    addr_up_hbra = offset_up_hbra
    addr_low_hbra = 0
    ! recurrence relations along x-direction
    do order_hbra = low_orders_hbra(1), low_orders_hbra(2)
      do ihbra = order_hbra, 0, -1
        do jhbra = ihbra, 0, -1
          addr_up_hbra = addr_up_hbra+1
          addr_low_hbra = addr_low_hbra+1
          up_hket_pints(addr_low_hbra,addr_up_hket) &
            = neg_ratio_braket*low_hket_pints(addr_up_hbra,addr_low_hket)
        end do
        addr_up_hbra = addr_up_hbra+1
      end do
      addr_up_hbra = addr_up_hbra+1
    end do
    ! recurrence relations along y-direction
    do ihket = 0, low_order_hket
      addr_up_hket = addr_up_hket+1
      addr_up_hbra = offset_up_hbra+1
      addr_low_hbra = 0
      do order_hbra = low_orders_hbra(1), low_orders_hbra(2)
        do ihbra = order_hbra, 0, -1
          do jhbra = 0, ihbra
            addr_up_hbra = addr_up_hbra+1
            addr_low_hbra = addr_low_hbra+1
            up_hket_pints(addr_low_hbra,addr_up_hket) &
              = neg_ratio_braket*low_hket_pints(addr_up_hbra,addr_low_hket+ihket)
          end do
          addr_up_hbra = addr_up_hbra+1
        end do
        addr_up_hbra = addr_up_hbra+1
      end do
    end do
    ! recurrence relations along z-direction
    addr_low_hket = addr_low_hket-1
    do ihket = 0, low_order_hket
      do jhket = 0, low_order_hket-ihket
        addr_low_hket = addr_low_hket+1
        addr_up_hket = addr_up_hket+1
        addr_up_hbra = offset_up_hbra+low_orders_hbra(1)+2
        addr_low_hbra = 0
        do order_hbra = low_orders_hbra(1), low_orders_hbra(2)
          do ihbra = 0, order_hbra
            do jhbra = ihbra, order_hbra
              addr_up_hbra = addr_up_hbra+1
              addr_low_hbra = addr_low_hbra+1
              up_hket_pints(addr_low_hbra,addr_up_hket) &
                = neg_ratio_braket*low_hket_pints(addr_up_hbra,addr_low_hket)
            end do
          end do
          addr_up_hbra = addr_up_hbra+order_hbra+3
        end do
      end do
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "sub_carmom_hrr_ket", STDOUT)
#endif
    return
  end subroutine sub_carmom_hrr_ket

  !> \brief assigns the integrals with required HGTOs on ket center
  !> \author Bin Gao
  !> \date 2012-03-04
  !> \param dim_low_hbra is the dimension of lower order HGTOs on bra center in temporary integrals
  !> \param num_up_hket is the number of upper order HGTOs on ket center in temporary integrals
  !> \param recur_pints contains the temporary integrals from recurrence relations
  !> \param offset_hket is the offset of HGTOs on ket center in returned integrals
  !> \param dim_hgto_bra is the dimension of HGTOs on bra center
  !> \param dim_hgto_ket is the dimension of HGTOs on ket center
  !> \return hket_pints contains the integrals with specific orders of HGTOs on bra and ket centers
  subroutine carmom_hrr_ket_assign(dim_low_hbra, num_up_hket, recur_pints,  &
                                   offset_hket, dim_hgto_bra, dim_hgto_ket, &
                                   hket_pints)
    use xkind
    implicit none
    integer, intent(in) :: dim_low_hbra
    integer, intent(in) :: num_up_hket
    real(REALK), intent(in) :: recur_pints(dim_low_hbra,num_up_hket)
    integer, intent(in) :: offset_hket
    integer, intent(in) :: dim_hgto_bra
    integer, intent(in) :: dim_hgto_ket
    real(REALK), intent(inout) :: hket_pints(dim_hgto_bra,dim_hgto_ket)
!f2py intent(hide) :: dim_low_hbra
!f2py intent(hide) :: num_up_hket
!f2py intent(in) :: recur_pints
!f2py intent(in) :: offset_hket
!f2py intent(hide) :: dim_hgto_bra
!f2py intent(hide) :: dim_hgto_ket
!f2py intent(inout) :: hket_pints
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    hket_pints(:,offset_hket+1:offset_hket+num_up_hket) = recur_pints(1:dim_hgto_bra,:)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "carmom_hrr_ket_assign", STDOUT)
#endif
    return
  end subroutine carmom_hrr_ket_assign
