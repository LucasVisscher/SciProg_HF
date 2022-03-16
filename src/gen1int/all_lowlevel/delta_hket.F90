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
!!  This file recovers HGTOs on ket center in Dirac delta function integrals.
!!
!!  2012-03-16, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief recovers HGTOs on ket center in Dirac delta function integrals
  !> \author Bin Gao
  !> \date 2012-03-16
  !> \param orders_hgto_ket is the range of orders of Hermite Gaussians on ket center
  !> \param orders_geo_pot is the range of orders of geometric derivatives on Dirac delta function
  !> \param coord_ket contains the coordinates of ket center
  !> \param exponent_ket is the exponent of primitive Gaussian of ket center
  !> \param delta_origin contains the coordinates of Dirac delta function origin
  !> \param dim_hgto_bra is the dimension of HGTOs on bra center
  !> \param dim_geo_hbra is the dimension of geometric derivatives on Dirac delta function
  !> \param hbra_pints contains the primitive Hermite integrals with zeroth order HGTO on ket center
  !> \param dim_hgto_ket is the dimension of HGTOs on ket center
  !> \param dim_geo_hket is dimension of geometric derivatives on Dirac delta function afterwards
  !> \return hket_pints contains the primitive Hermite integrals with specified
  !>         orders of HGTOs on ket center and geometric derivatives on Dirac delta function
  subroutine delta_hket(orders_hgto_ket, orders_geo_pot, coord_ket, exponent_ket, &
                        delta_origin, dim_hgto_bra, dim_geo_hbra, hbra_pints,     &
                        dim_hgto_ket, dim_geo_hket, hket_pints)
    use xkind
    implicit none
    integer, intent(in) :: orders_hgto_ket(2)
    integer, intent(in) :: orders_geo_pot(2)
    real(REALK), intent(in) :: coord_ket(3)
    real(REALK), intent(in) :: exponent_ket
    real(REALK), intent(in) :: delta_origin(3)
    integer, intent(in) :: dim_hgto_bra
    integer, intent(in) :: dim_geo_hbra
    real(REALK), intent(in) :: hbra_pints(dim_hgto_bra,dim_geo_hbra)
    integer, intent(in) :: dim_hgto_ket
    integer, intent(in) :: dim_geo_hket
    real(REALK), intent(out) :: hket_pints(dim_hgto_bra,dim_hgto_ket,dim_geo_hket)
!f2py intent(in) :: orders_hgto_ket
!f2py intent(in) :: orders_geo_pot
!f2py intent(in) :: coord_ket
!f2py intent(in) :: exponent_ket
!f2py intent(in) :: delta_origin
!f2py intent(hide) :: dim_hgto_bra
!f2py intent(hide) :: dim_geo_hbra
!f2py intent(in) :: hbra_pints
!f2py intent(in) :: dim_hgto_ket
!f2py intent(in) :: dim_geo_hket
!f2py intent(out) :: hket_pints
!f2py depend(dim_hgto_bra) :: hket_pints
!f2py depend(dim_hgto_ket) :: hket_pints
!f2py depend(dim_geo_hket) :: hket_pints
    logical zero_hket             !if returning zeroth order HGTO on ket center
    integer min_geo_pot           !minimum geometric derivatives on Dirac delta function
    integer max_cur_hket          !maximum of current order HGTOs on ket center
    real(REALK) half_nr_expnt     !half of the negative reciprocal of exponent on ket center
    real(REALK) delta_wrt_ket(3)  !relative coordinates of Dirac delta function origin w.r.t. ket center
    integer dim_cur_hket          !dimension of current order HGTOs on ket center
    integer dim_up_hket           !dimension of upper order HGTOs on ket center
    integer num_low_geo           !number of xyz components of lower order geometric derivatives
    integer num_cur_geo           !number of xyz components of current order geometric derivatives
    integer size_low_geo          !size of temporary integrals of lower order geometric derivatives
    integer size_cur_geo          !size of temporary integrals of current order geometric derivatives
    integer low_geo_int           !pointer of lower order geometric derivatives
    integer cur_geo_int           !pointer of current order geometric derivatives
    real(REALK), allocatable :: cur_zero_pints(:,:)
                                  !temporary integrals with zeroth order HGTOs on ket center
    real(REALK), allocatable :: tmp_ints(:,:,:)
                                  !temporary integrals
    integer start_low_hbra        !start address of lower order geometric derivatives in \var(hbra_pints)
    integer end_low_hbra          !end address of lower order geometric derivatives in \var(hbra_pints)
    integer start_cur_hbra        !start address of current order geometric derivatives in \var(hbra_pints)
    integer end_cur_hbra          !end address of current order geometric derivatives in \var(hbra_pints)
    integer start_cur_hket        !start address of current order geometric derivatives in var(hket_pints)
    integer end_cur_hket          !end address of current order geometric derivatives in var(hket_pints)
    integer cur_order_geo         !incremental recorder over current order of geometric derivatives
    integer ierr                  !error information
#if defined(XTIME)
    real(REALK) curr_time         !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! only zeroth order HGTO on ket center returns
    if (orders_hgto_ket(2)==0) then
      hket_pints(:,1,:) = hbra_pints
    else
      ! if returning zeroth order HGTO on ket center
      zero_hket = orders_hgto_ket(1)==0
      ! sets the minimum geometric derivatives on Dirac delta function
      min_geo_pot = orders_geo_pot(1)-orders_hgto_ket(2)+1
      ! sets the maximum of current order HGTOs on ket center
      max_cur_hket = -min_geo_pot
      ! sets the half of the negative reciprocal of exponent on ket center
      half_nr_expnt = -0.5_REALK/exponent_ket
      ! sets the relative coordinates of Dirac delta function origin w.r.t. ket center
      delta_wrt_ket = delta_origin-coord_ket
      ! allocates memory for temporary integrals
      dim_up_hket = (orders_hgto_ket(2)+1)*(orders_hgto_ket(2)+2) &
                  * (orders_hgto_ket(2)+3)/6-1
      num_cur_geo = (orders_geo_pot(2)+1)*(orders_geo_pot(2)+2)/2
      allocate(cur_zero_pints(dim_hgto_bra,num_cur_geo), stat=ierr)
      if (ierr/=0)                                                &
        call error_stop("delta_hket", "failed to cur_zero_pints", &
                        dim_hgto_bra*num_cur_geo)
      allocate(tmp_ints(dim_hgto_bra,dim_up_hket*num_cur_geo,2), stat=ierr)
      if (ierr/=0)                                          &
        call error_stop("delta_hket", "failed to tmp_ints", &
                        dim_hgto_bra*dim_up_hket*num_cur_geo*2)
      ! gets temporary integrals with zeroth order geometric derivative on Dirac delta function
      if (min_geo_pot<=0) then
        ! sets the dimension of upper order HGTOs on ket center
        dim_up_hket = (max_cur_hket+2)*(max_cur_hket+3)*(max_cur_hket+4)/6-1
        ! sets the temporary integrals with zeroth order HGTOs on ket center,
        ! which will be changed during recurrence relations
        cur_zero_pints(:,1) = hbra_pints(:,1)
        call zero_delta_hket(max_cur_hket, half_nr_expnt, delta_wrt_ket,     &
                             dim_hgto_bra, cur_zero_pints(:,1), dim_up_hket, &
                             tmp_ints(:,1:dim_up_hket,1))
        ! sets the minimum geometric derivatives on Dirac delta function
        min_geo_pot = 1
        ! sets number of xyz components of current order geometric derivatives
        num_cur_geo = 1
        ! sets the start and end addresses of current order geometric derivatives in \var(hbra_pints)
        start_cur_hbra = 1
        end_cur_hbra = 1
        ! sets the size of temporary integrals of current order geometric derivatives
        size_cur_geo = dim_up_hket
        ! initializes the pointers
        low_geo_int = 2
        cur_geo_int = 1
      else
        max_cur_hket = -1
        dim_up_hket = 0
        num_cur_geo = min_geo_pot*(min_geo_pot+1)/2
        start_cur_hbra = 1
        end_cur_hbra = num_cur_geo
        size_cur_geo = 0
        low_geo_int = 1
        cur_geo_int = 2
      end if
      ! initializes the start and end addresses of HGTOs on ket center in returned integrals
      start_cur_hket = 0
      end_cur_hket = 0
      ! loops over the orders of geometric derivatives on Dirac delta function
      ! till the minimum returned order, in which the maximum of current order
      ! HGTOs on ket center increases
      do cur_order_geo = min_geo_pot, orders_geo_pot(1)
        ! updates the maximum of current order HGTOs on ket center
        max_cur_hket = max_cur_hket+1
        ! updates the dimensions of HGTOs on ket center
        dim_cur_hket = dim_up_hket
        dim_up_hket = dim_up_hket+(max_cur_hket+2)*(max_cur_hket+3)/2
        ! updates the number of geometric derivatives on Dirac delta function
        num_low_geo = num_cur_geo
        num_cur_geo = num_cur_geo+cur_order_geo+1
        ! updates the start and end addresses of geometric derivatives in \var(hbra_pints)
        start_low_hbra = start_cur_hbra
        end_low_hbra = end_cur_hbra
        start_cur_hbra = end_cur_hbra+1
        end_cur_hbra = end_cur_hbra+num_cur_geo
        ! sets the size of temporary integrals
        size_low_geo = size_cur_geo
        size_cur_geo = dim_up_hket*num_cur_geo
        ! switch the pointers
        low_geo_int = 3-low_geo_int
        cur_geo_int = 3-cur_geo_int
        ! sets the temporary integrals with zeroth order HGTOs on ket center,
        ! which will be changed during recurrence relations
        cur_zero_pints(:,1:num_cur_geo) = hbra_pints(:,start_cur_hbra:end_cur_hbra)
        ! gets the temporary integrals with \var(cur_order_geo) order geometric derivatives
        ! on Dirac delta function
        call sub_delta_hket(max_cur_hket, cur_order_geo, half_nr_expnt,           &
                            delta_wrt_ket, dim_hgto_bra, num_low_geo,             &
                            hbra_pints(:,start_low_hbra:end_low_hbra),            &
                            dim_cur_hket, tmp_ints(:,1:size_low_geo,low_geo_int), &
                            num_cur_geo, cur_zero_pints(:,1:num_cur_geo),         &
                            dim_up_hket, tmp_ints(:,1:size_cur_geo,cur_geo_int))
      end do
      ! assigns integrals of minimum order of geometric derivatives on Dirac delta function
      start_cur_hket = end_cur_hket+1
      end_cur_hket = end_cur_hket+num_cur_geo
      call delta_hket_assign(dim_hgto_bra, num_cur_geo,              &
             hbra_pints(:,start_cur_hbra:end_cur_hbra), dim_up_hket, &
             tmp_ints(:,1:size_cur_geo,cur_geo_int), zero_hket,      &
             dim_hgto_ket, hket_pints(:,:,start_cur_hket:end_cur_hket))
      ! updates the dimension of current order HGTOs on ket center
      dim_cur_hket = dim_up_hket
      ! loops over other orders of geometric derivatives on Dirac delta function,
      ! in which the maximum of current order HGTOs on ket center does not change
      do cur_order_geo = orders_geo_pot(1)+1, orders_geo_pot(2)
        ! updates the number of geometric derivatives on Dirac delta function
        num_low_geo = num_cur_geo
        num_cur_geo = num_cur_geo+cur_order_geo+1
        ! updates the start and end addresses of geometric derivatives in \var(hbra_pints)
        start_low_hbra = start_cur_hbra
        end_low_hbra = end_cur_hbra
        start_cur_hbra = end_cur_hbra+1
        end_cur_hbra = end_cur_hbra+num_cur_geo
        ! sets the size of temporary integrals
        size_low_geo = size_cur_geo
        size_cur_geo = dim_up_hket*num_cur_geo
        ! switch the pointers
        low_geo_int = 3-low_geo_int
        cur_geo_int = 3-cur_geo_int
        ! sets the temporary integrals with zeroth order HGTOs on ket center,
        ! which will be changed during recurrence relations
        cur_zero_pints(:,1:num_cur_geo) = hbra_pints(:,start_cur_hbra:end_cur_hbra)
        ! gets the temporary integrals with \var(cur_order_geo) order geometric derivatives
        ! on Dirac delta function
        call sub_delta_hket(max_cur_hket, cur_order_geo, half_nr_expnt,           &
                            delta_wrt_ket, dim_hgto_bra, num_low_geo,             &
                            hbra_pints(:,start_low_hbra:end_low_hbra),            &
                            dim_cur_hket, tmp_ints(:,1:size_low_geo,low_geo_int), &
                            num_cur_geo, cur_zero_pints(:,1:num_cur_geo),         &
                            dim_up_hket, tmp_ints(:,1:size_cur_geo,cur_geo_int))
        ! assigns integrals back
        start_cur_hket = end_cur_hket+1
        end_cur_hket = end_cur_hket+num_cur_geo
        call delta_hket_assign(dim_hgto_bra, num_cur_geo,              &
               hbra_pints(:,start_cur_hbra:end_cur_hbra), dim_up_hket, &
               tmp_ints(:,1:size_cur_geo,cur_geo_int), zero_hket,      &
               dim_hgto_ket, hket_pints(:,:,start_cur_hket:end_cur_hket))
      end do
      ! cleans
      deallocate(cur_zero_pints)
      deallocate(tmp_ints)
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "delta_hket", STDOUT)
#endif
    return
  end subroutine delta_hket

  !> \brief recovers HGTOs on ket center for the zeroth order geometric derivatives
  !>        on Dirac delta function
  !> \author Bin Gao
  !> \date 2012-03-16
  !> \param max_cur_hket is the maximum of current order HGTOs on ket center
  !> \param half_nr_expnt is the half of the negative reciprocal of exponent on ket center
  !> \param delta_wrt_ket contains the relative coordinates of Dirac delta function
  !>        w.r.t. ket center
  !> \param dim_hgto_bra is the dimension of HGTOs on bra center
  !> \param zero_geo_pints contains integrals of zeroth order geometric derivatives on
  !>        Dirac delta function and zeroth order HGTOs on ket center
  !> \param dim_up_hket is the dimension of upper order HGTOs on ket center
  !> \return cur_geo_pints contains the integrals with orders of HGTOs on ket center from
  !>         1 to \var(max_cur_hket)+1
  subroutine zero_delta_hket(max_cur_hket, half_nr_expnt, delta_wrt_ket, &
                             dim_hgto_bra, zero_geo_pints, dim_up_hket,  &
                             cur_geo_pints)
    use xkind
    implicit none
    integer, intent(in) :: max_cur_hket
    real(REALK), intent(in) :: half_nr_expnt
    real(REALK), intent(in) :: delta_wrt_ket(3)
    integer, intent(in) :: dim_hgto_bra
    real(REALK), intent(inout) :: zero_geo_pints(dim_hgto_bra)
    integer, intent(in) :: dim_up_hket
    real(REALK), intent(out) :: cur_geo_pints(dim_hgto_bra,dim_up_hket)
!f2py intent(in) :: max_cur_hket
!f2py intent(in) :: half_nr_expnt
!f2py intent(in) :: delta_wrt_ket
!f2py intent(hide) :: dim_hgto_bra
!f2py intent(inout) :: zero_geo_pints
!f2py intent(in) :: dim_up_hket
!f2py intent(out) :: cur_geo_pints
!f2py depend(dim_hgto_bra) :: cur_geo_pints
!f2py depend(dim_up_hket) :: cur_geo_pints
    integer base_low_hket  !base address of lower order HGTOs on ket center
    integer base_cur_hket  !base address of current order HGTOs on ket center
    integer addr_low_hket  !address of lower order HGTOs on ket center
    integer addr_cur_hket  !address of current order HGTOs on ket center
    integer addr_up_hket   !address of upper order HGTOs on ket center
    integer order_hket     !incremental recorder over the orders of HGTOs on ket center
    integer iket, jket     !incremental recorders over xyz components of HGTOs on ket center
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time   
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    select case(max_cur_hket)
    ! only the first order HGTOs on ket center required
    case(0)
#if defined(DEBUG)
      write(STDOUT,100) "1st order HGTOs return"
#endif
      cur_geo_pints(:,1) = delta_wrt_ket(1)*zero_geo_pints  !px
      cur_geo_pints(:,2) = delta_wrt_ket(2)*zero_geo_pints  !py
      cur_geo_pints(:,3) = delta_wrt_ket(3)*zero_geo_pints  !pz
    ! the first and second order HGTOs on ket center
    case(1)
#if defined(DEBUG)
      write(STDOUT,100) "1st and 2nd order HGTOs return"
#endif
      cur_geo_pints(:,1) = delta_wrt_ket(1)*zero_geo_pints  !px
      cur_geo_pints(:,2) = delta_wrt_ket(2)*zero_geo_pints  !py
      cur_geo_pints(:,3) = delta_wrt_ket(3)*zero_geo_pints  !pz
      zero_geo_pints = half_nr_expnt*zero_geo_pints
      cur_geo_pints(:,4) = delta_wrt_ket(1)*cur_geo_pints(:,1) &  !dxx
                         + zero_geo_pints
      cur_geo_pints(:,5) = delta_wrt_ket(2)*cur_geo_pints(:,1)    !dxy
      cur_geo_pints(:,6) = delta_wrt_ket(2)*cur_geo_pints(:,2) &  !dyy
                         + zero_geo_pints
      cur_geo_pints(:,7) = delta_wrt_ket(3)*cur_geo_pints(:,1)    !dxz
      cur_geo_pints(:,8) = delta_wrt_ket(3)*cur_geo_pints(:,2)    !dyz
      cur_geo_pints(:,9) = delta_wrt_ket(3)*cur_geo_pints(:,3) &  !dzz
                         + zero_geo_pints
    ! the maximum order of HGTOs on ket center required is > 2
    case default
#if defined(DEBUG)
      write(STDOUT,100) "higher order HGTOs return", max_cur_hket
#endif
      cur_geo_pints(:,1) = delta_wrt_ket(1)*zero_geo_pints  !px
      cur_geo_pints(:,2) = delta_wrt_ket(2)*zero_geo_pints  !py
      cur_geo_pints(:,3) = delta_wrt_ket(3)*zero_geo_pints  !pz
      zero_geo_pints = half_nr_expnt*zero_geo_pints
      cur_geo_pints(:,4) = delta_wrt_ket(1)*cur_geo_pints(:,1) &  !dxx
                         + zero_geo_pints
      cur_geo_pints(:,5) = delta_wrt_ket(2)*cur_geo_pints(:,1)    !dxy
      cur_geo_pints(:,6) = delta_wrt_ket(2)*cur_geo_pints(:,2) &  !dyy
                         + zero_geo_pints
      cur_geo_pints(:,7) = delta_wrt_ket(3)*cur_geo_pints(:,1)    !dxz
      cur_geo_pints(:,8) = delta_wrt_ket(3)*cur_geo_pints(:,2)    !dyz
      cur_geo_pints(:,9) = delta_wrt_ket(3)*cur_geo_pints(:,3) &  !dzz
                         + zero_geo_pints
      ! initializes the (base) addresses of HGTOs on ket center
      base_low_hket = 0
      base_cur_hket = 3
      addr_up_hket = 9
      ! other order (>2) HGTOs on ket center
      do order_hket = 2, max_cur_hket
        ! recurrence relation along x-direction
        addr_up_hket = addr_up_hket+1
        addr_cur_hket = base_cur_hket+1
        cur_geo_pints(:,addr_up_hket)                       &
          = delta_wrt_ket(1)*cur_geo_pints(:,addr_cur_hket) &
          + half_nr_expnt*real(order_hket,REALK)*cur_geo_pints(:,base_low_hket+1)
        ! recurrence relation along y-direction
        addr_up_hket = addr_up_hket+1
        cur_geo_pints(:,addr_up_hket) &
          = delta_wrt_ket(2)*cur_geo_pints(:,addr_cur_hket)
        addr_low_hket = base_low_hket
        do iket = 1, order_hket
          addr_low_hket = addr_low_hket+1
          addr_cur_hket = addr_cur_hket+1
          addr_up_hket = addr_up_hket+1
          cur_geo_pints(:,addr_up_hket)                       &
            = delta_wrt_ket(2)*cur_geo_pints(:,addr_cur_hket) &
            + half_nr_expnt*real(iket,REALK)*cur_geo_pints(:,addr_low_hket)
        end do
        ! recurrence relation along z-direction
        addr_cur_hket = base_cur_hket
        do jket = 0, order_hket
          addr_cur_hket = addr_cur_hket+1
          addr_up_hket = addr_up_hket+1
          cur_geo_pints(:,addr_up_hket) &
            = delta_wrt_ket(3)*cur_geo_pints(:,addr_cur_hket)
        end do
        addr_low_hket = base_low_hket
        do iket = 1, order_hket
          do jket = 0, order_hket-iket
            addr_low_hket = addr_low_hket+1
            addr_cur_hket = addr_cur_hket+1
            addr_up_hket = addr_up_hket+1
            cur_geo_pints(:,addr_up_hket)                       &
              = delta_wrt_ket(3)*cur_geo_pints(:,addr_cur_hket) &
              + half_nr_expnt*real(iket,REALK)*cur_geo_pints(:,addr_low_hket)
          end do
        end do
        ! updates the base addresses
        base_low_hket = addr_low_hket
        base_cur_hket = addr_cur_hket
      end do
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "zero_delta_hket", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("zero_delta_hket>> ",A,I6)
#endif
  end subroutine zero_delta_hket

  !> \brief recovers HGTOs on ket center for a given order (>0) of geometric derivatives
  !>        on Dirac delta function
  !> \author Bin Gao
  !> \date 2012-03-16
  !> \param max_cur_hket is the maximum of current order HGTOs on ket center
  !> \param cur_order_geo is the current order of geometric derivatives on Dirac delta function
  !> \param half_nr_expnt is the half of the negative reciprocal of exponent on ket center
  !> \param delta_wrt_ket contains the relative coordinates of Dirac delta function
  !>        w.r.t. ket center
  !> \param dim_hgto_bra is the dimension of HGTOs on bra center
  !> \param num_low_geo is the number of lower order geometric derivatives on Dirac delta function
  !> \param low_zero_pints contains integrals of lower order geometric derivatives on
  !>        Dirac delta function and zeroth order HGTOs on ket center
  !> \param dim_cur_hket is the dimension of current order HGTOs on ket center
  !> \param low_geo_pints contains integrals of lower order geometric derivatives on
  !>        Dirac delta function and lower order HGTOs on ket center
  !> \param num_cur_geo is the number of current order geometric derivatives on Dirac delta function
  !> \param cur_zero_pints contains integrals of current order geometric derivatives on
  !>        Dirac delta function and zeroth order HGTOs on ket center
  !> \param dim_up_hket is the dimension of upper order HGTOs on ket center
  !> \return cur_geo_pints contains the integrals of current order geometric derivatives on
  !>         Dirac delta function and orders of HGTOs on ket center from 1 to \var(max_cur_hket)+1
  subroutine sub_delta_hket(max_cur_hket, cur_order_geo, half_nr_expnt, delta_wrt_ket, &
                            dim_hgto_bra, num_low_geo, low_zero_pints, dim_cur_hket,   &
                            low_geo_pints, num_cur_geo, cur_zero_pints, dim_up_hket,   &
                            cur_geo_pints)
    use xkind
    implicit none
    integer, intent(in) :: max_cur_hket
    integer, intent(in) :: cur_order_geo
    real(REALK), intent(in) :: half_nr_expnt
    real(REALK), intent(in) :: delta_wrt_ket(3)
    integer, intent(in) :: dim_hgto_bra
    integer, intent(in) :: num_low_geo
    real(REALK), intent(in) :: low_zero_pints(dim_hgto_bra,num_low_geo)
    integer, intent(in) :: dim_cur_hket
    real(REALK), intent(in) :: low_geo_pints(dim_hgto_bra,dim_cur_hket,num_low_geo)
    integer, intent(in) :: num_cur_geo
    real(REALK), intent(inout) :: cur_zero_pints(dim_hgto_bra,num_cur_geo)
    integer, intent(in) :: dim_up_hket
    real(REALK), intent(out) :: cur_geo_pints(dim_hgto_bra,dim_up_hket,num_cur_geo)
!f2py intent(in) :: max_cur_hket
!f2py intent(in) :: cur_order_geo
!f2py intent(in) :: half_nr_expnt
!f2py intent(in) :: delta_wrt_ket
!f2py intent(hide) :: dim_hgto_bra
!f2py intent(hide) :: num_low_geo
!f2py intent(in) :: low_zero_pints
!f2py intent(hide) :: dim_cur_hket
!f2py intent(in) :: low_geo_pints
!f2py depend(dim_hgto_bra) :: low_geo_pints
!f2py depend(num_low_geo) :: low_geo_pints
!f2py intent(hide) :: num_cur_geo
!f2py intent(inout) :: cur_zero_pints
!f2py depend(dim_hgto_bra) :: cur_zero_pints
!f2py intent(in) :: dim_up_hket
!f2py intent(out) :: cur_geo_pints
!f2py depend(dim_hgto_bra) :: cur_geo_pints
!f2py depend(dim_up_hket) :: cur_geo_pints
!f2py depend(num_cur_geo) :: cur_geo_pints
    integer addr_low_xgeo  !address of lower order geometric derivatives along x-direction
    integer addr_low_ygeo  !address of lower order geometric derivatives along y-direction
    integer addr_low_zgeo  !address of lower order geometric derivatives along z-direction
    integer addr_cur_geo   !address of current order geometric derivatives
    integer igeo, jgeo     !incremental recorders of address of geometric derivatives
    integer order_geo_xy   !order of geometric derivatives along x- and y-direction
    integer order_geo_x    !order of geometric derivatives along x-direction
    integer base_low_hket  !base address of lower order HGTOs on ket center
    integer base_cur_hket  !base address of current order HGTOs on ket center
    integer addr_low_hket  !address of lower order HGTOs on ket center
    integer addr_cur_hket  !address of current order HGTOs on ket center
    integer addr_up_hket   !address of upper order HGTOs on ket center
    integer order_hket     !incremental recorder over the orders of HGTOs on ket center
    integer iket, jket     !incremental recorders over xyz components of HGTOs on ket center
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time   
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    select case(max_cur_hket)
    ! only the first order HGTOs on ket center required
    case(0)
#if defined(DEBUG)
      write(STDOUT,100) "1st order HGTOs return"
      write(STDOUT,100) "order of geometric derivatives", cur_order_geo
#endif
      ! (1) x...x component of geometric derivatives
      cur_geo_pints(:,1,1) = delta_wrt_ket(1)*cur_zero_pints(:,1) &         !px
                           + real(cur_order_geo,REALK)*low_zero_pints(:,1)
      cur_geo_pints(:,2,1) = delta_wrt_ket(2)*cur_zero_pints(:,1)           !py
      cur_geo_pints(:,3,1) = delta_wrt_ket(3)*cur_zero_pints(:,1)           !pz
      ! (2) x...xy to xy...y components of geometric derivatives
      addr_low_xgeo = 1
      addr_low_ygeo = 0
      addr_cur_geo = 1
      do igeo = 1, cur_order_geo-1
        addr_low_xgeo = addr_low_xgeo+1
        addr_low_ygeo = addr_low_ygeo+1
        addr_cur_geo = addr_cur_geo+1
        cur_geo_pints(:,1,addr_cur_geo)                     &  !px
          = delta_wrt_ket(1)*cur_zero_pints(:,addr_cur_geo) &
          + real(cur_order_geo-igeo,REALK)*low_zero_pints(:,addr_low_xgeo)
        cur_geo_pints(:,2,addr_cur_geo)                     &  !py
          = delta_wrt_ket(2)*cur_zero_pints(:,addr_cur_geo) &
          + real(igeo,REALK)*low_zero_pints(:,addr_low_ygeo)
        cur_geo_pints(:,3,addr_cur_geo) &                      !pz
          = delta_wrt_ket(3)*cur_zero_pints(:,addr_cur_geo)
      end do
      ! (3) y...y component of geometric derivatives
      addr_low_ygeo = addr_low_ygeo+1
      addr_cur_geo = addr_cur_geo+1
      cur_geo_pints(:,1,addr_cur_geo) &                      !px
        = delta_wrt_ket(1)*cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,2,addr_cur_geo)                     &  !py
        = delta_wrt_ket(2)*cur_zero_pints(:,addr_cur_geo) &  
        + real(cur_order_geo,REALK)*low_zero_pints(:,addr_low_ygeo)
      cur_geo_pints(:,3,addr_cur_geo) &                      !pz
        = delta_wrt_ket(3)*cur_zero_pints(:,addr_cur_geo)    
      ! (4) x...xz...z to y...yz...z components of geometric derivatives
      addr_low_zgeo = 0
      do igeo = 1, cur_order_geo-1
        ! (4.1) x...xz...z component of geometric derivatives
        addr_low_xgeo = addr_low_xgeo+1
        addr_low_zgeo = addr_low_zgeo+1
        addr_cur_geo = addr_cur_geo+1
        order_geo_xy = cur_order_geo-igeo
        cur_geo_pints(:,1,addr_cur_geo)                     &  !px
          = delta_wrt_ket(1)*cur_zero_pints(:,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_zero_pints(:,addr_low_xgeo)
        cur_geo_pints(:,2,addr_cur_geo) &                      !py
          = delta_wrt_ket(2)*cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,3,addr_cur_geo)                     &  !pz
          = delta_wrt_ket(3)*cur_zero_pints(:,addr_cur_geo) &
          + real(igeo,REALK)*low_zero_pints(:,addr_low_zgeo)
        ! (4.2) x...xyz...z to xy...yz...z components of geometric derivatives
        do jgeo = 1, cur_order_geo-(igeo+1)
          addr_low_xgeo = addr_low_xgeo+1
          addr_low_ygeo = addr_low_ygeo+1
          addr_low_zgeo = addr_low_zgeo+1
          addr_cur_geo = addr_cur_geo+1
          cur_geo_pints(:,1,addr_cur_geo)                     &  !px
            = delta_wrt_ket(1)*cur_zero_pints(:,addr_cur_geo) &
            + real(order_geo_xy-jgeo,REALK)*low_zero_pints(:,addr_low_xgeo)
          cur_geo_pints(:,2,addr_cur_geo)                     &  !py
            = delta_wrt_ket(2)*cur_zero_pints(:,addr_cur_geo) &
            + real(jgeo,REALK)*low_zero_pints(:,addr_low_ygeo)
          cur_geo_pints(:,3,addr_cur_geo)                     &  !pz
            = delta_wrt_ket(3)*cur_zero_pints(:,addr_cur_geo) &
            + real(igeo,REALK)*low_zero_pints(:,addr_low_zgeo)
        end do
        ! (4.3) y...yz...z component of geometric derivatives
        addr_low_ygeo = addr_low_ygeo+1
        addr_low_zgeo = addr_low_zgeo+1
        addr_cur_geo = addr_cur_geo+1
        cur_geo_pints(:,1,addr_cur_geo) &                      !px
          = delta_wrt_ket(1)*cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,2,addr_cur_geo)                     &  !py
          = delta_wrt_ket(2)*cur_zero_pints(:,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_zero_pints(:,addr_low_ygeo)
        cur_geo_pints(:,3,addr_cur_geo)                     &  !pz
          = delta_wrt_ket(3)*cur_zero_pints(:,addr_cur_geo) &
          + real(igeo,REALK)*low_zero_pints(:,addr_low_zgeo)
      end do
      ! (5) z...z component of geometric derivatives
      addr_low_zgeo = addr_low_zgeo+1
      addr_cur_geo = addr_cur_geo+1
      cur_geo_pints(:,1,addr_cur_geo) &                      !px
        = delta_wrt_ket(1)*cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,2,addr_cur_geo) &                      !py
        = delta_wrt_ket(2)*cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,3,addr_cur_geo)                     &  !pz
        = delta_wrt_ket(3)*cur_zero_pints(:,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_zero_pints(:,addr_low_zgeo)
    ! the first and second order HGTOs on ket center
    case(1)
#if defined(DEBUG)
      write(STDOUT,100) "1st and 2nd order HGTOs return"
      write(STDOUT,100) "order of geometric derivatives", cur_order_geo
#endif
      ! (1) x...x component of geometric derivatives
      cur_geo_pints(:,1,1) = delta_wrt_ket(1)*cur_zero_pints(:,1) &         !px
                           + real(cur_order_geo,REALK)*low_zero_pints(:,1)
      cur_geo_pints(:,2,1) = delta_wrt_ket(2)*cur_zero_pints(:,1)           !py
      cur_geo_pints(:,3,1) = delta_wrt_ket(3)*cur_zero_pints(:,1)           !pz
      cur_zero_pints(:,1) = half_nr_expnt*cur_zero_pints(:,1)
      cur_geo_pints(:,4,1) = delta_wrt_ket(1)*cur_geo_pints(:,1,1)          &  !dxx
                           + real(cur_order_geo,REALK)*low_geo_pints(:,1,1) &
                           + cur_zero_pints(:,1)
      cur_geo_pints(:,5,1) = delta_wrt_ket(2)*cur_geo_pints(:,1,1)             !dxy
      cur_geo_pints(:,6,1) = delta_wrt_ket(2)*cur_geo_pints(:,2,1) &           !dyy
                           + cur_zero_pints(:,1)
      cur_geo_pints(:,7,1) = delta_wrt_ket(3)*cur_geo_pints(:,1,1)             !dxz
      cur_geo_pints(:,8,1) = delta_wrt_ket(3)*cur_geo_pints(:,2,1)             !dyz
      cur_geo_pints(:,9,1) = delta_wrt_ket(3)*cur_geo_pints(:,3,1) &           !dzz
                           + cur_zero_pints(:,1)
      ! (2) x...xy to xy...y components of geometric derivatives
      addr_low_xgeo = 1
      addr_low_ygeo = 0
      addr_cur_geo = 1
      do igeo = 1, cur_order_geo-1
        addr_low_xgeo = addr_low_xgeo+1
        addr_low_ygeo = addr_low_ygeo+1
        addr_cur_geo = addr_cur_geo+1
        order_geo_x = cur_order_geo-igeo
        cur_geo_pints(:,1,addr_cur_geo)                     &  !px
          = delta_wrt_ket(1)*cur_zero_pints(:,addr_cur_geo) &
          + real(order_geo_x,REALK)*low_zero_pints(:,addr_low_xgeo)
        cur_geo_pints(:,2,addr_cur_geo)                     &  !py
          = delta_wrt_ket(2)*cur_zero_pints(:,addr_cur_geo) &
          + real(igeo,REALK)*low_zero_pints(:,addr_low_ygeo)
        cur_geo_pints(:,3,addr_cur_geo) &                      !pz
          = delta_wrt_ket(3)*cur_zero_pints(:,addr_cur_geo)
        cur_zero_pints(:,addr_cur_geo) = half_nr_expnt*cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,4,addr_cur_geo)                              &  !dxx
          = delta_wrt_ket(1)*cur_geo_pints(:,1,addr_cur_geo)         &
          + real(order_geo_x,REALK)*low_geo_pints(:,1,addr_low_xgeo) &
          + cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,5,addr_cur_geo)                      &          !dxy
          = delta_wrt_ket(2)*cur_geo_pints(:,1,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,1,addr_low_ygeo)
        cur_geo_pints(:,6,addr_cur_geo)                       &         !dyy
          = delta_wrt_ket(2)*cur_geo_pints(:,2,addr_cur_geo)  &
          + real(igeo,REALK)*low_geo_pints(:,2,addr_low_ygeo) &
          + cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,7,addr_cur_geo) &                               !dxz
          = delta_wrt_ket(3)*cur_geo_pints(:,1,addr_cur_geo)
        cur_geo_pints(:,8,addr_cur_geo) &                               !dyz
          = delta_wrt_ket(3)*cur_geo_pints(:,2,addr_cur_geo)
        cur_geo_pints(:,9,addr_cur_geo)                      &          !dzz
          = delta_wrt_ket(3)*cur_geo_pints(:,3,addr_cur_geo) &
          + cur_zero_pints(:,addr_cur_geo)
      end do
      ! (3) y...y component of geometric derivatives
      addr_low_ygeo = addr_low_ygeo+1
      addr_cur_geo = addr_cur_geo+1
      cur_geo_pints(:,1,addr_cur_geo) &                      !px
        = delta_wrt_ket(1)*cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,2,addr_cur_geo)                     &  !py
        = delta_wrt_ket(2)*cur_zero_pints(:,addr_cur_geo) &  
        + real(cur_order_geo,REALK)*low_zero_pints(:,addr_low_ygeo)
      cur_geo_pints(:,3,addr_cur_geo) &                      !pz
        = delta_wrt_ket(3)*cur_zero_pints(:,addr_cur_geo)
      cur_zero_pints(:,addr_cur_geo) = half_nr_expnt*cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,4,addr_cur_geo)                      &            !dxx
        = delta_wrt_ket(1)*cur_geo_pints(:,1,addr_cur_geo) &
        + cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,5,addr_cur_geo)                      &            !dxy
        = delta_wrt_ket(2)*cur_geo_pints(:,1,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_geo_pints(:,1,addr_low_ygeo)
      cur_geo_pints(:,6,addr_cur_geo)                                &  !dyy
        = delta_wrt_ket(2)*cur_geo_pints(:,2,addr_cur_geo)           &
        + real(cur_order_geo,REALK)*low_geo_pints(:,2,addr_low_ygeo) &
        + cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,7,addr_cur_geo) &                                 !dxz
        = delta_wrt_ket(3)*cur_geo_pints(:,1,addr_cur_geo)
      cur_geo_pints(:,8,addr_cur_geo) &                                 !dyz
        = delta_wrt_ket(3)*cur_geo_pints(:,2,addr_cur_geo)
      cur_geo_pints(:,9,addr_cur_geo)                      &            !dzz
        = delta_wrt_ket(3)*cur_geo_pints(:,3,addr_cur_geo) &
        + cur_zero_pints(:,addr_cur_geo)
      ! (4) x...xz...z to y...yz...z components of geometric derivatives
      addr_low_zgeo = 0
      do igeo = 1, cur_order_geo-1
        ! (4.1) x...xz...z component of geometric derivatives
        addr_low_xgeo = addr_low_xgeo+1
        addr_low_zgeo = addr_low_zgeo+1
        addr_cur_geo = addr_cur_geo+1
        order_geo_xy = cur_order_geo-igeo
        cur_geo_pints(:,1,addr_cur_geo)                     &  !px
          = delta_wrt_ket(1)*cur_zero_pints(:,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_zero_pints(:,addr_low_xgeo)
        cur_geo_pints(:,2,addr_cur_geo) &                      !py
          = delta_wrt_ket(2)*cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,3,addr_cur_geo)                     &  !pz
          = delta_wrt_ket(3)*cur_zero_pints(:,addr_cur_geo) &
          + real(igeo,REALK)*low_zero_pints(:,addr_low_zgeo)
        cur_zero_pints(:,addr_cur_geo) = half_nr_expnt*cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,4,addr_cur_geo)                               &  !dxx
          = delta_wrt_ket(1)*cur_geo_pints(:,1,addr_cur_geo)          &
          + real(order_geo_xy,REALK)*low_geo_pints(:,1,addr_low_xgeo) &
          + cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,5,addr_cur_geo) &                                !dxy
          = delta_wrt_ket(2)*cur_geo_pints(:,1,addr_cur_geo)
        cur_geo_pints(:,6,addr_cur_geo)                      &           !dyy
          = delta_wrt_ket(2)*cur_geo_pints(:,2,addr_cur_geo) &
          + cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,7,addr_cur_geo)                      &           !dxz
          = delta_wrt_ket(3)*cur_geo_pints(:,1,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,1,addr_low_zgeo)
        cur_geo_pints(:,8,addr_cur_geo)                      &           !dyz
          = delta_wrt_ket(3)*cur_geo_pints(:,2,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,2,addr_low_zgeo)
        cur_geo_pints(:,9,addr_cur_geo)                       &          !dzz
          = delta_wrt_ket(3)*cur_geo_pints(:,3,addr_cur_geo)  &
          + real(igeo,REALK)*low_geo_pints(:,3,addr_low_zgeo) &
          + cur_zero_pints(:,addr_cur_geo)
        ! (4.2) x...xyz...z to xy...yz...z components of geometric derivatives
        do jgeo = 1, cur_order_geo-(igeo+1)
          addr_low_xgeo = addr_low_xgeo+1
          addr_low_ygeo = addr_low_ygeo+1
          addr_low_zgeo = addr_low_zgeo+1
          addr_cur_geo = addr_cur_geo+1
          cur_geo_pints(:,1,addr_cur_geo)                     &  !px
            = delta_wrt_ket(1)*cur_zero_pints(:,addr_cur_geo) &
            + real(order_geo_xy-jgeo,REALK)*low_zero_pints(:,addr_low_xgeo)
          cur_geo_pints(:,2,addr_cur_geo)                     &  !py
            = delta_wrt_ket(2)*cur_zero_pints(:,addr_cur_geo) &
            + real(jgeo,REALK)*low_zero_pints(:,addr_low_ygeo)
          cur_geo_pints(:,3,addr_cur_geo)                     &  !pz
            = delta_wrt_ket(3)*cur_zero_pints(:,addr_cur_geo) &
            + real(igeo,REALK)*low_zero_pints(:,addr_low_zgeo)
          cur_zero_pints(:,addr_cur_geo) = half_nr_expnt*cur_zero_pints(:,addr_cur_geo)
          cur_geo_pints(:,4,addr_cur_geo)                                    &  !dxx
            = delta_wrt_ket(1)*cur_geo_pints(:,1,addr_cur_geo)               &
            + real(order_geo_xy-jgeo,REALK)*low_geo_pints(:,1,addr_low_xgeo) &
            + cur_zero_pints(:,addr_cur_geo)
          cur_geo_pints(:,5,addr_cur_geo)                      &                !dxy
            = delta_wrt_ket(2)*cur_geo_pints(:,1,addr_cur_geo) &
            + real(jgeo,REALK)*low_geo_pints(:,1,addr_low_ygeo)
          cur_geo_pints(:,6,addr_cur_geo)                       &               !dyy
            = delta_wrt_ket(2)*cur_geo_pints(:,2,addr_cur_geo)  &
            + real(jgeo,REALK)*low_geo_pints(:,2,addr_low_ygeo) &
            + cur_zero_pints(:,addr_cur_geo)
          cur_geo_pints(:,7,addr_cur_geo)                      &                !dxz
            = delta_wrt_ket(3)*cur_geo_pints(:,1,addr_cur_geo) &
            + real(igeo,REALK)*low_geo_pints(:,1,addr_low_zgeo)
          cur_geo_pints(:,8,addr_cur_geo)                      &                !dyz
            = delta_wrt_ket(3)*cur_geo_pints(:,2,addr_cur_geo) &
            + real(igeo,REALK)*low_geo_pints(:,2,addr_low_zgeo)
          cur_geo_pints(:,9,addr_cur_geo)                       &               !dzz
            = delta_wrt_ket(3)*cur_geo_pints(:,3,addr_cur_geo)  &
            + real(igeo,REALK)*low_geo_pints(:,3,addr_low_zgeo) &
            + cur_zero_pints(:,addr_cur_geo)
        end do
        ! (4.3) y...yz...z component of geometric derivatives
        addr_low_ygeo = addr_low_ygeo+1
        addr_low_zgeo = addr_low_zgeo+1
        addr_cur_geo = addr_cur_geo+1
        cur_geo_pints(:,1,addr_cur_geo) &                      !px
          = delta_wrt_ket(1)*cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,2,addr_cur_geo)                     &  !py
          = delta_wrt_ket(2)*cur_zero_pints(:,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_zero_pints(:,addr_low_ygeo)
        cur_geo_pints(:,3,addr_cur_geo)                     &  !pz
          = delta_wrt_ket(3)*cur_zero_pints(:,addr_cur_geo) &
          + real(igeo,REALK)*low_zero_pints(:,addr_low_zgeo)
        cur_zero_pints(:,addr_cur_geo) = half_nr_expnt*cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,4,addr_cur_geo)                      &           !dxx
          = delta_wrt_ket(1)*cur_geo_pints(:,1,addr_cur_geo) &
          + cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,5,addr_cur_geo)                      &           !dxy
          = delta_wrt_ket(2)*cur_geo_pints(:,1,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_geo_pints(:,1,addr_low_ygeo)
        cur_geo_pints(:,6,addr_cur_geo)                               &  !dyy
          = delta_wrt_ket(2)*cur_geo_pints(:,2,addr_cur_geo)          &
          + real(order_geo_xy,REALK)*low_geo_pints(:,2,addr_low_ygeo) &
          + cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,7,addr_cur_geo)                      &           !dxz
          = delta_wrt_ket(3)*cur_geo_pints(:,1,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,1,addr_low_zgeo)
        cur_geo_pints(:,8,addr_cur_geo)                      &           !dyz
          = delta_wrt_ket(3)*cur_geo_pints(:,2,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,2,addr_low_zgeo)
        cur_geo_pints(:,9,addr_cur_geo)                       &          !dzz
          = delta_wrt_ket(3)*cur_geo_pints(:,3,addr_cur_geo)  &
          + real(igeo,REALK)*low_geo_pints(:,3,addr_low_zgeo) &
          + cur_zero_pints(:,addr_cur_geo)
      end do
      ! (5) z...z component of geometric derivatives
      addr_low_zgeo = addr_low_zgeo+1
      addr_cur_geo = addr_cur_geo+1
      cur_geo_pints(:,1,addr_cur_geo) &                      !px
        = delta_wrt_ket(1)*cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,2,addr_cur_geo) &                      !py
        = delta_wrt_ket(2)*cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,3,addr_cur_geo)                     &  !pz
        = delta_wrt_ket(3)*cur_zero_pints(:,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_zero_pints(:,addr_low_zgeo)
      cur_zero_pints(:,addr_cur_geo) = half_nr_expnt*cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,4,addr_cur_geo)                      &            !dxx
        = delta_wrt_ket(1)*cur_geo_pints(:,1,addr_cur_geo) &
        + cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,5,addr_cur_geo) &                                 !dxy
        = delta_wrt_ket(2)*cur_geo_pints(:,1,addr_cur_geo)
      cur_geo_pints(:,6,addr_cur_geo)                      &            !dyy
        = delta_wrt_ket(2)*cur_geo_pints(:,2,addr_cur_geo) &
        + cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,7,addr_cur_geo)                      &            !dxz
        = delta_wrt_ket(3)*cur_geo_pints(:,1,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_geo_pints(:,1,addr_low_zgeo)
      cur_geo_pints(:,8,addr_cur_geo)                      &            !dyz
        = delta_wrt_ket(3)*cur_geo_pints(:,2,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_geo_pints(:,2,addr_low_zgeo)
      cur_geo_pints(:,9,addr_cur_geo)                                &  !dzz
        = delta_wrt_ket(3)*cur_geo_pints(:,3,addr_cur_geo)           &
        + real(cur_order_geo,REALK)*low_geo_pints(:,3,addr_low_zgeo) &
        + cur_zero_pints(:,addr_cur_geo)
    ! the maximum order of HGTOs on ket center required is > 2
    case default
#if defined(DEBUG)
      write(STDOUT,100) "higher order HGTOs return", max_cur_hket
      write(STDOUT,100) "order of geometric derivatives", cur_order_geo
#endif
      ! (1) x...x component of geometric derivatives
      cur_geo_pints(:,1,1) = delta_wrt_ket(1)*cur_zero_pints(:,1) &         !px
                           + real(cur_order_geo,REALK)*low_zero_pints(:,1)
      cur_geo_pints(:,2,1) = delta_wrt_ket(2)*cur_zero_pints(:,1)           !py
      cur_geo_pints(:,3,1) = delta_wrt_ket(3)*cur_zero_pints(:,1)           !pz
      cur_zero_pints(:,1) = half_nr_expnt*cur_zero_pints(:,1)
      cur_geo_pints(:,4,1) = delta_wrt_ket(1)*cur_geo_pints(:,1,1)          &  !dxx
                           + real(cur_order_geo,REALK)*low_geo_pints(:,1,1) &
                           + cur_zero_pints(:,1)
      cur_geo_pints(:,5,1) = delta_wrt_ket(2)*cur_geo_pints(:,1,1)             !dxy
      cur_geo_pints(:,6,1) = delta_wrt_ket(2)*cur_geo_pints(:,2,1) &           !dyy
                           + cur_zero_pints(:,1)
      cur_geo_pints(:,7,1) = delta_wrt_ket(3)*cur_geo_pints(:,1,1)             !dxz
      cur_geo_pints(:,8,1) = delta_wrt_ket(3)*cur_geo_pints(:,2,1)             !dyz
      cur_geo_pints(:,9,1) = delta_wrt_ket(3)*cur_geo_pints(:,3,1) &           !dzz
                           + cur_zero_pints(:,1)
      ! initializes the (base) addresses of HGTOs on ket center
      base_low_hket = 0
      base_cur_hket = 3
      addr_up_hket = 9
      ! other order (>2) HGTOs on ket center
      do order_hket = 2, max_cur_hket
        ! recurrence relation along x-direction
        addr_up_hket = addr_up_hket+1
        addr_cur_hket = base_cur_hket+1
        cur_geo_pints(:,addr_up_hket,1)                                &
          = delta_wrt_ket(1)*cur_geo_pints(:,addr_cur_hket,1)          &
          + real(cur_order_geo,REALK)*low_geo_pints(:,addr_cur_hket,1) &
          + half_nr_expnt*real(order_hket,REALK)*cur_geo_pints(:,base_low_hket+1,1)
        ! recurrence relation along y-direction
        addr_up_hket = addr_up_hket+1
        cur_geo_pints(:,addr_up_hket,1) &
          = delta_wrt_ket(2)*cur_geo_pints(:,addr_cur_hket,1)
        addr_low_hket = base_low_hket
        do iket = 1, order_hket
          addr_low_hket = addr_low_hket+1
          addr_cur_hket = addr_cur_hket+1
          addr_up_hket = addr_up_hket+1
          cur_geo_pints(:,addr_up_hket,1)                       &
            = delta_wrt_ket(2)*cur_geo_pints(:,addr_cur_hket,1) &
            + half_nr_expnt*real(iket,REALK)*cur_geo_pints(:,addr_low_hket,1)
        end do
        ! recurrence relation along z-direction
        addr_cur_hket = base_cur_hket
        do jket = 0, order_hket
          addr_cur_hket = addr_cur_hket+1
          addr_up_hket = addr_up_hket+1
          cur_geo_pints(:,addr_up_hket,1) &
            = delta_wrt_ket(3)*cur_geo_pints(:,addr_cur_hket,1)
        end do
        addr_low_hket = base_low_hket
        do iket = 1, order_hket
          do jket = 0, order_hket-iket
            addr_low_hket = addr_low_hket+1
            addr_cur_hket = addr_cur_hket+1
            addr_up_hket = addr_up_hket+1
            cur_geo_pints(:,addr_up_hket,1)                       &
              = delta_wrt_ket(3)*cur_geo_pints(:,addr_cur_hket,1) &
              + half_nr_expnt*real(iket,REALK)*cur_geo_pints(:,addr_low_hket,1)
          end do
        end do
        ! updates the base addresses
        base_low_hket = addr_low_hket
        base_cur_hket = addr_cur_hket
      end do
      ! (2) x...xy to xy...y components of geometric derivatives
      addr_low_xgeo = 1
      addr_low_ygeo = 0
      addr_cur_geo = 1
      do igeo = 1, cur_order_geo-1
        addr_low_xgeo = addr_low_xgeo+1
        addr_low_ygeo = addr_low_ygeo+1
        addr_cur_geo = addr_cur_geo+1
        order_geo_x = cur_order_geo-igeo
        cur_geo_pints(:,1,addr_cur_geo)                     &  !px
          = delta_wrt_ket(1)*cur_zero_pints(:,addr_cur_geo) &
          + real(order_geo_x,REALK)*low_zero_pints(:,addr_low_xgeo)
        cur_geo_pints(:,2,addr_cur_geo)                     &  !py
          = delta_wrt_ket(2)*cur_zero_pints(:,addr_cur_geo) &
          + real(igeo,REALK)*low_zero_pints(:,addr_low_ygeo)
        cur_geo_pints(:,3,addr_cur_geo) &                      !pz
          = delta_wrt_ket(3)*cur_zero_pints(:,addr_cur_geo)
        cur_zero_pints(:,addr_cur_geo) = half_nr_expnt*cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,4,addr_cur_geo)                              &  !dxx
          = delta_wrt_ket(1)*cur_geo_pints(:,1,addr_cur_geo)         &
          + real(order_geo_x,REALK)*low_geo_pints(:,1,addr_low_xgeo) &
          + cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,5,addr_cur_geo)                      &          !dxy
          = delta_wrt_ket(2)*cur_geo_pints(:,1,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,1,addr_low_ygeo)
        cur_geo_pints(:,6,addr_cur_geo)                       &         !dyy
          = delta_wrt_ket(2)*cur_geo_pints(:,2,addr_cur_geo)  &
          + real(igeo,REALK)*low_geo_pints(:,2,addr_low_ygeo) &
          + cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,7,addr_cur_geo) &                               !dxz
          = delta_wrt_ket(3)*cur_geo_pints(:,1,addr_cur_geo)
        cur_geo_pints(:,8,addr_cur_geo) &                               !dyz
          = delta_wrt_ket(3)*cur_geo_pints(:,2,addr_cur_geo)
        cur_geo_pints(:,9,addr_cur_geo)                      &          !dzz
          = delta_wrt_ket(3)*cur_geo_pints(:,3,addr_cur_geo) &
          + cur_zero_pints(:,addr_cur_geo)
        ! initializes the (base) addresses of HGTOs on ket center
        base_low_hket = 0
        base_cur_hket = 3
        addr_up_hket = 9
        ! other order (>2) HGTOs on ket center
        do order_hket = 2, max_cur_hket
          ! recurrence relation along x-direction
          addr_up_hket = addr_up_hket+1
          addr_cur_hket = base_cur_hket+1
          cur_geo_pints(:,addr_up_hket,addr_cur_geo)                                &
            = delta_wrt_ket(1)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo)          &
            + real(order_geo_x,REALK)*low_geo_pints(:,addr_cur_hket,addr_low_xgeo)  &
            + half_nr_expnt*real(order_hket,REALK)                                  &
            * cur_geo_pints(:,base_low_hket+1,addr_cur_geo)
          ! recurrence relation along y-direction
          addr_up_hket = addr_up_hket+1
          cur_geo_pints(:,addr_up_hket,addr_cur_geo)                       &
            = delta_wrt_ket(2)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo) &
            + real(igeo,REALK)*low_geo_pints(:,addr_cur_hket,addr_low_ygeo)
          addr_low_hket = base_low_hket
          do iket = 1, order_hket
            addr_low_hket = addr_low_hket+1
            addr_cur_hket = addr_cur_hket+1
            addr_up_hket = addr_up_hket+1
            cur_geo_pints(:,addr_up_hket,addr_cur_geo)                        &
              = delta_wrt_ket(2)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo)  &
              + real(igeo,REALK)*low_geo_pints(:,addr_cur_hket,addr_low_ygeo) &
              + half_nr_expnt*real(iket,REALK)                                &
              * cur_geo_pints(:,addr_low_hket,addr_cur_geo)
          end do
          ! recurrence relation along z-direction
          addr_cur_hket = base_cur_hket
          do jket = 0, order_hket
            addr_cur_hket = addr_cur_hket+1
            addr_up_hket = addr_up_hket+1
            cur_geo_pints(:,addr_up_hket,addr_cur_geo) &
              = delta_wrt_ket(3)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo)
          end do
          addr_low_hket = base_low_hket
          do iket = 1, order_hket
            do jket = 0, order_hket-iket
              addr_low_hket = addr_low_hket+1
              addr_cur_hket = addr_cur_hket+1
              addr_up_hket = addr_up_hket+1
              cur_geo_pints(:,addr_up_hket,addr_cur_geo)                       &
                = delta_wrt_ket(3)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo) &
                + half_nr_expnt*real(iket,REALK)                               &
                * cur_geo_pints(:,addr_low_hket,addr_cur_geo)
            end do
          end do
          ! updates the base addresses
          base_low_hket = addr_low_hket
          base_cur_hket = addr_cur_hket
        end do
      end do
      ! (3) y...y component of geometric derivatives
      addr_low_ygeo = addr_low_ygeo+1
      addr_cur_geo = addr_cur_geo+1
      cur_geo_pints(:,1,addr_cur_geo) &                      !px
        = delta_wrt_ket(1)*cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,2,addr_cur_geo)                     &  !py
        = delta_wrt_ket(2)*cur_zero_pints(:,addr_cur_geo) &  
        + real(cur_order_geo,REALK)*low_zero_pints(:,addr_low_ygeo)
      cur_geo_pints(:,3,addr_cur_geo) &                      !pz
        = delta_wrt_ket(3)*cur_zero_pints(:,addr_cur_geo)
      cur_zero_pints(:,addr_cur_geo) = half_nr_expnt*cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,4,addr_cur_geo)                      &            !dxx
        = delta_wrt_ket(1)*cur_geo_pints(:,1,addr_cur_geo) &
        + cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,5,addr_cur_geo)                      &            !dxy
        = delta_wrt_ket(2)*cur_geo_pints(:,1,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_geo_pints(:,1,addr_low_ygeo)
      cur_geo_pints(:,6,addr_cur_geo)                                &  !dyy
        = delta_wrt_ket(2)*cur_geo_pints(:,2,addr_cur_geo)           &
        + real(cur_order_geo,REALK)*low_geo_pints(:,2,addr_low_ygeo) &
        + cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,7,addr_cur_geo) &                                 !dxz
        = delta_wrt_ket(3)*cur_geo_pints(:,1,addr_cur_geo)
      cur_geo_pints(:,8,addr_cur_geo) &                                 !dyz
        = delta_wrt_ket(3)*cur_geo_pints(:,2,addr_cur_geo)
      cur_geo_pints(:,9,addr_cur_geo)                      &            !dzz
        = delta_wrt_ket(3)*cur_geo_pints(:,3,addr_cur_geo) &
        + cur_zero_pints(:,addr_cur_geo)
      ! initializes the (base) addresses of HGTOs on ket center
      base_low_hket = 0
      base_cur_hket = 3
      addr_up_hket = 9
      ! other order (>2) HGTOs on ket center
      do order_hket = 2, max_cur_hket
        ! recurrence relation along x-direction
        addr_up_hket = addr_up_hket+1
        addr_cur_hket = base_cur_hket+1
        cur_geo_pints(:,addr_up_hket,addr_cur_geo)                       &
          = delta_wrt_ket(1)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo) &
          + half_nr_expnt*real(order_hket,REALK)                         &
          * cur_geo_pints(:,base_low_hket+1,addr_cur_geo)
        ! recurrence relation along y-direction
        addr_up_hket = addr_up_hket+1
        cur_geo_pints(:,addr_up_hket,addr_cur_geo)                       &
          = delta_wrt_ket(2)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo) &
          + real(cur_order_geo,REALK)*low_geo_pints(:,addr_cur_hket,addr_low_ygeo)
        addr_low_hket = base_low_hket
        do iket = 1, order_hket
          addr_low_hket = addr_low_hket+1
          addr_cur_hket = addr_cur_hket+1
          addr_up_hket = addr_up_hket+1
          cur_geo_pints(:,addr_up_hket,addr_cur_geo)                                 &
            = delta_wrt_ket(2)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo)           &
            + real(cur_order_geo,REALK)*low_geo_pints(:,addr_cur_hket,addr_low_ygeo) &
            + half_nr_expnt*real(iket,REALK)                                         &
            * cur_geo_pints(:,addr_low_hket,addr_cur_geo)
        end do
        ! recurrence relation along z-direction
        addr_cur_hket = base_cur_hket
        do jket = 0, order_hket
          addr_cur_hket = addr_cur_hket+1
          addr_up_hket = addr_up_hket+1
          cur_geo_pints(:,addr_up_hket,addr_cur_geo) &
            = delta_wrt_ket(3)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo)
        end do
        addr_low_hket = base_low_hket
        do iket = 1, order_hket
          do jket = 0, order_hket-iket
            addr_low_hket = addr_low_hket+1
            addr_cur_hket = addr_cur_hket+1
            addr_up_hket = addr_up_hket+1
            cur_geo_pints(:,addr_up_hket,addr_cur_geo)                       &
              = delta_wrt_ket(3)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo) &
              + half_nr_expnt*real(iket,REALK)                               &
              * cur_geo_pints(:,addr_low_hket,addr_cur_geo)
          end do
        end do
        ! updates the base addresses
        base_low_hket = addr_low_hket
        base_cur_hket = addr_cur_hket
      end do
      ! (4) x...xz...z to y...yz...z components of geometric derivatives
      addr_low_zgeo = 0
      do igeo = 1, cur_order_geo-1
        ! (4.1) x...xz...z component of geometric derivatives
        addr_low_xgeo = addr_low_xgeo+1
        addr_low_zgeo = addr_low_zgeo+1
        addr_cur_geo = addr_cur_geo+1
        order_geo_xy = cur_order_geo-igeo
        cur_geo_pints(:,1,addr_cur_geo)                     &  !px
          = delta_wrt_ket(1)*cur_zero_pints(:,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_zero_pints(:,addr_low_xgeo)
        cur_geo_pints(:,2,addr_cur_geo) &                      !py
          = delta_wrt_ket(2)*cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,3,addr_cur_geo)                     &  !pz
          = delta_wrt_ket(3)*cur_zero_pints(:,addr_cur_geo) &
          + real(igeo,REALK)*low_zero_pints(:,addr_low_zgeo)
        cur_zero_pints(:,addr_cur_geo) = half_nr_expnt*cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,4,addr_cur_geo)                               &  !dxx
          = delta_wrt_ket(1)*cur_geo_pints(:,1,addr_cur_geo)          &
          + real(order_geo_xy,REALK)*low_geo_pints(:,1,addr_low_xgeo) &
          + cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,5,addr_cur_geo) &                                !dxy
          = delta_wrt_ket(2)*cur_geo_pints(:,1,addr_cur_geo)
        cur_geo_pints(:,6,addr_cur_geo)                      &           !dyy
          = delta_wrt_ket(2)*cur_geo_pints(:,2,addr_cur_geo) &
          + cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,7,addr_cur_geo)                      &           !dxz
          = delta_wrt_ket(3)*cur_geo_pints(:,1,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,1,addr_low_zgeo)
        cur_geo_pints(:,8,addr_cur_geo)                      &           !dyz
          = delta_wrt_ket(3)*cur_geo_pints(:,2,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,2,addr_low_zgeo)
        cur_geo_pints(:,9,addr_cur_geo)                       &          !dzz
          = delta_wrt_ket(3)*cur_geo_pints(:,3,addr_cur_geo)  &
          + real(igeo,REALK)*low_geo_pints(:,3,addr_low_zgeo) &
          + cur_zero_pints(:,addr_cur_geo)
        ! initializes the (base) addresses of HGTOs on ket center
        base_low_hket = 0
        base_cur_hket = 3
        addr_up_hket = 9
        ! other order (>2) HGTOs on ket center
        do order_hket = 2, max_cur_hket
          ! recurrence relation along x-direction
          addr_up_hket = addr_up_hket+1
          addr_cur_hket = base_cur_hket+1
          cur_geo_pints(:,addr_up_hket,addr_cur_geo)                                &
            = delta_wrt_ket(1)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo)          &
            + real(order_geo_xy,REALK)*low_geo_pints(:,addr_cur_hket,addr_low_xgeo) &
            + half_nr_expnt*real(order_hket,REALK)                                  &
            * cur_geo_pints(:,base_low_hket+1,addr_cur_geo)
          ! recurrence relation along y-direction
          addr_up_hket = addr_up_hket+1
          cur_geo_pints(:,addr_up_hket,addr_cur_geo) &
            = delta_wrt_ket(2)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo)
          addr_low_hket = base_low_hket
          do iket = 1, order_hket
            addr_low_hket = addr_low_hket+1
            addr_cur_hket = addr_cur_hket+1
            addr_up_hket = addr_up_hket+1
            cur_geo_pints(:,addr_up_hket,addr_cur_geo)                       &
              = delta_wrt_ket(2)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo) &
              + half_nr_expnt*real(iket,REALK)                               &
              * cur_geo_pints(:,addr_low_hket,addr_cur_geo)
          end do
          ! recurrence relation along z-direction
          addr_cur_hket = base_cur_hket
          do jket = 0, order_hket
            addr_cur_hket = addr_cur_hket+1
            addr_up_hket = addr_up_hket+1
            cur_geo_pints(:,addr_up_hket,addr_cur_geo)                       &
              = delta_wrt_ket(3)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo) &
              + real(igeo,REALK)*low_geo_pints(:,addr_cur_hket,addr_low_zgeo)
          end do
          addr_low_hket = base_low_hket
          do iket = 1, order_hket
            do jket = 0, order_hket-iket
              addr_low_hket = addr_low_hket+1
              addr_cur_hket = addr_cur_hket+1
              addr_up_hket = addr_up_hket+1
              cur_geo_pints(:,addr_up_hket,addr_cur_geo)                        &
                = delta_wrt_ket(3)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo)  &
                + real(igeo,REALK)*low_geo_pints(:,addr_cur_hket,addr_low_zgeo) &
                + half_nr_expnt*real(iket,REALK)                                &
                * cur_geo_pints(:,addr_low_hket,addr_cur_geo)
            end do
          end do
          ! updates the base addresses
          base_low_hket = addr_low_hket
          base_cur_hket = addr_cur_hket
        end do
        ! (4.2) x...xyz...z to xy...yz...z components of geometric derivatives
        do jgeo = 1, cur_order_geo-(igeo+1)
          addr_low_xgeo = addr_low_xgeo+1
          addr_low_ygeo = addr_low_ygeo+1
          addr_low_zgeo = addr_low_zgeo+1
          addr_cur_geo = addr_cur_geo+1
          cur_geo_pints(:,1,addr_cur_geo)                     &  !px
            = delta_wrt_ket(1)*cur_zero_pints(:,addr_cur_geo) &
            + real(order_geo_xy-jgeo,REALK)*low_zero_pints(:,addr_low_xgeo)
          cur_geo_pints(:,2,addr_cur_geo)                     &  !py
            = delta_wrt_ket(2)*cur_zero_pints(:,addr_cur_geo) &
            + real(jgeo,REALK)*low_zero_pints(:,addr_low_ygeo)
          cur_geo_pints(:,3,addr_cur_geo)                     &  !pz
            = delta_wrt_ket(3)*cur_zero_pints(:,addr_cur_geo) &
            + real(igeo,REALK)*low_zero_pints(:,addr_low_zgeo)
          cur_zero_pints(:,addr_cur_geo) = half_nr_expnt*cur_zero_pints(:,addr_cur_geo)
          cur_geo_pints(:,4,addr_cur_geo)                                    &  !dxx
            = delta_wrt_ket(1)*cur_geo_pints(:,1,addr_cur_geo)               &
            + real(order_geo_xy-jgeo,REALK)*low_geo_pints(:,1,addr_low_xgeo) &
            + cur_zero_pints(:,addr_cur_geo)
          cur_geo_pints(:,5,addr_cur_geo)                      &                !dxy
            = delta_wrt_ket(2)*cur_geo_pints(:,1,addr_cur_geo) &
            + real(jgeo,REALK)*low_geo_pints(:,1,addr_low_ygeo)
          cur_geo_pints(:,6,addr_cur_geo)                       &               !dyy
            = delta_wrt_ket(2)*cur_geo_pints(:,2,addr_cur_geo)  &
            + real(jgeo,REALK)*low_geo_pints(:,2,addr_low_ygeo) &
            + cur_zero_pints(:,addr_cur_geo)
          cur_geo_pints(:,7,addr_cur_geo)                      &                !dxz
            = delta_wrt_ket(3)*cur_geo_pints(:,1,addr_cur_geo) &
            + real(igeo,REALK)*low_geo_pints(:,1,addr_low_zgeo)
          cur_geo_pints(:,8,addr_cur_geo)                      &                !dyz
            = delta_wrt_ket(3)*cur_geo_pints(:,2,addr_cur_geo) &
            + real(igeo,REALK)*low_geo_pints(:,2,addr_low_zgeo)
          cur_geo_pints(:,9,addr_cur_geo)                       &               !dzz
            = delta_wrt_ket(3)*cur_geo_pints(:,3,addr_cur_geo)  &
            + real(igeo,REALK)*low_geo_pints(:,3,addr_low_zgeo) &
            + cur_zero_pints(:,addr_cur_geo)
          ! initializes the (base) addresses of HGTOs on ket center
          base_low_hket = 0
          base_cur_hket = 3
          addr_up_hket = 9
          ! other order (>2) HGTOs on ket center
          do order_hket = 2, max_cur_hket
            ! recurrence relation along x-direction
            addr_up_hket = addr_up_hket+1
            addr_cur_hket = base_cur_hket+1
            cur_geo_pints(:,addr_up_hket,addr_cur_geo)                       &
              = delta_wrt_ket(1)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo) &
              + real(order_geo_xy-jgeo,REALK)                                &
              * low_geo_pints(:,addr_cur_hket,addr_low_xgeo)                 &
              + half_nr_expnt*real(order_hket,REALK)                         &
              * cur_geo_pints(:,base_low_hket+1,addr_cur_geo)
            ! recurrence relation along y-direction
            addr_up_hket = addr_up_hket+1
            cur_geo_pints(:,addr_up_hket,addr_cur_geo)                       &
              = delta_wrt_ket(2)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo) &
              + real(jgeo,REALK)*low_geo_pints(:,addr_cur_hket,addr_low_ygeo)
            addr_low_hket = base_low_hket
            do iket = 1, order_hket
              addr_low_hket = addr_low_hket+1
              addr_cur_hket = addr_cur_hket+1
              addr_up_hket = addr_up_hket+1
              cur_geo_pints(:,addr_up_hket,addr_cur_geo)                        &
                = delta_wrt_ket(2)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo)  &
                + real(jgeo,REALK)*low_geo_pints(:,addr_cur_hket,addr_low_ygeo) &
                + half_nr_expnt*real(iket,REALK)                                &
                * cur_geo_pints(:,addr_low_hket,addr_cur_geo)
            end do
            ! recurrence relation along z-direction
            addr_cur_hket = base_cur_hket
            do jket = 0, order_hket
              addr_cur_hket = addr_cur_hket+1
              addr_up_hket = addr_up_hket+1
              cur_geo_pints(:,addr_up_hket,addr_cur_geo)                       &
                = delta_wrt_ket(3)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo) &
                + real(igeo,REALK)*low_geo_pints(:,addr_cur_hket,addr_low_zgeo)
            end do
            addr_low_hket = base_low_hket
            do iket = 1, order_hket
              do jket = 0, order_hket-iket
                addr_low_hket = addr_low_hket+1
                addr_cur_hket = addr_cur_hket+1
                addr_up_hket = addr_up_hket+1
                cur_geo_pints(:,addr_up_hket,addr_cur_geo)                        &
                  = delta_wrt_ket(3)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo)  &
                  + real(igeo,REALK)*low_geo_pints(:,addr_cur_hket,addr_low_zgeo) &
                  + half_nr_expnt*real(iket,REALK)                                &
                  * cur_geo_pints(:,addr_low_hket,addr_cur_geo)
              end do
            end do
            ! updates the base addresses
            base_low_hket = addr_low_hket
            base_cur_hket = addr_cur_hket
          end do
        end do
        ! (4.3) y...yz...z component of geometric derivatives
        addr_low_ygeo = addr_low_ygeo+1
        addr_low_zgeo = addr_low_zgeo+1
        addr_cur_geo = addr_cur_geo+1
        cur_geo_pints(:,1,addr_cur_geo) &                      !px
          = delta_wrt_ket(1)*cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,2,addr_cur_geo)                     &  !py
          = delta_wrt_ket(2)*cur_zero_pints(:,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_zero_pints(:,addr_low_ygeo)
        cur_geo_pints(:,3,addr_cur_geo)                     &  !pz
          = delta_wrt_ket(3)*cur_zero_pints(:,addr_cur_geo) &
          + real(igeo,REALK)*low_zero_pints(:,addr_low_zgeo)
        cur_zero_pints(:,addr_cur_geo) = half_nr_expnt*cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,4,addr_cur_geo)                      &           !dxx
          = delta_wrt_ket(1)*cur_geo_pints(:,1,addr_cur_geo) &
          + cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,5,addr_cur_geo)                      &           !dxy
          = delta_wrt_ket(2)*cur_geo_pints(:,1,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_geo_pints(:,1,addr_low_ygeo)
        cur_geo_pints(:,6,addr_cur_geo)                               &  !dyy
          = delta_wrt_ket(2)*cur_geo_pints(:,2,addr_cur_geo)          &
          + real(order_geo_xy,REALK)*low_geo_pints(:,2,addr_low_ygeo) &
          + cur_zero_pints(:,addr_cur_geo)
        cur_geo_pints(:,7,addr_cur_geo)                      &           !dxz
          = delta_wrt_ket(3)*cur_geo_pints(:,1,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,1,addr_low_zgeo)
        cur_geo_pints(:,8,addr_cur_geo)                      &           !dyz
          = delta_wrt_ket(3)*cur_geo_pints(:,2,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,2,addr_low_zgeo)
        cur_geo_pints(:,9,addr_cur_geo)                       &          !dzz
          = delta_wrt_ket(3)*cur_geo_pints(:,3,addr_cur_geo)  &
          + real(igeo,REALK)*low_geo_pints(:,3,addr_low_zgeo) &
          + cur_zero_pints(:,addr_cur_geo)
        ! initializes the (base) addresses of HGTOs on ket center
        base_low_hket = 0
        base_cur_hket = 3
        addr_up_hket = 9
        ! other order (>2) HGTOs on ket center
        do order_hket = 2, max_cur_hket
          ! recurrence relation along x-direction
          addr_up_hket = addr_up_hket+1
          addr_cur_hket = base_cur_hket+1
          cur_geo_pints(:,addr_up_hket,addr_cur_geo)                       &
            = delta_wrt_ket(1)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo) &
            + half_nr_expnt*real(order_hket,REALK)                         &
            * cur_geo_pints(:,base_low_hket+1,addr_cur_geo)
          ! recurrence relation along y-direction
          addr_up_hket = addr_up_hket+1
          cur_geo_pints(:,addr_up_hket,addr_cur_geo)                       &
            = delta_wrt_ket(2)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo) &
            + real(order_geo_xy,REALK)*low_geo_pints(:,addr_cur_hket,addr_low_ygeo)
          addr_low_hket = base_low_hket
          do iket = 1, order_hket
            addr_low_hket = addr_low_hket+1
            addr_cur_hket = addr_cur_hket+1
            addr_up_hket = addr_up_hket+1
            cur_geo_pints(:,addr_up_hket,addr_cur_geo)                                &
              = delta_wrt_ket(2)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo)          &
              + real(order_geo_xy,REALK)*low_geo_pints(:,addr_cur_hket,addr_low_ygeo) &
              + half_nr_expnt*real(iket,REALK)                                        &
              * cur_geo_pints(:,addr_low_hket,addr_cur_geo)
          end do
          ! recurrence relation along z-direction
          addr_cur_hket = base_cur_hket
          do jket = 0, order_hket
            addr_cur_hket = addr_cur_hket+1
            addr_up_hket = addr_up_hket+1
            cur_geo_pints(:,addr_up_hket,addr_cur_geo)                       &
              = delta_wrt_ket(3)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo) &
              + real(igeo,REALK)*low_geo_pints(:,addr_cur_hket,addr_low_zgeo)
          end do
          addr_low_hket = base_low_hket
          do iket = 1, order_hket
            do jket = 0, order_hket-iket
              addr_low_hket = addr_low_hket+1
              addr_cur_hket = addr_cur_hket+1
              addr_up_hket = addr_up_hket+1
              cur_geo_pints(:,addr_up_hket,addr_cur_geo)                        &
                = delta_wrt_ket(3)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo)  &
                + real(igeo,REALK)*low_geo_pints(:,addr_cur_hket,addr_low_zgeo) &
                + half_nr_expnt*real(iket,REALK)                                &
                * cur_geo_pints(:,addr_low_hket,addr_cur_geo)
            end do
          end do
          ! updates the base addresses
          base_low_hket = addr_low_hket
          base_cur_hket = addr_cur_hket
        end do
      end do
      ! (5) z...z component of geometric derivatives
      addr_low_zgeo = addr_low_zgeo+1
      addr_cur_geo = addr_cur_geo+1
      cur_geo_pints(:,1,addr_cur_geo) &                      !px
        = delta_wrt_ket(1)*cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,2,addr_cur_geo) &                      !py
        = delta_wrt_ket(2)*cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,3,addr_cur_geo)                     &  !pz
        = delta_wrt_ket(3)*cur_zero_pints(:,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_zero_pints(:,addr_low_zgeo)
      cur_zero_pints(:,addr_cur_geo) = half_nr_expnt*cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,4,addr_cur_geo)                      &            !dxx
        = delta_wrt_ket(1)*cur_geo_pints(:,1,addr_cur_geo) &
        + cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,5,addr_cur_geo) &                                 !dxy
        = delta_wrt_ket(2)*cur_geo_pints(:,1,addr_cur_geo)
      cur_geo_pints(:,6,addr_cur_geo)                      &            !dyy
        = delta_wrt_ket(2)*cur_geo_pints(:,2,addr_cur_geo) &
        + cur_zero_pints(:,addr_cur_geo)
      cur_geo_pints(:,7,addr_cur_geo)                      &            !dxz
        = delta_wrt_ket(3)*cur_geo_pints(:,1,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_geo_pints(:,1,addr_low_zgeo)
      cur_geo_pints(:,8,addr_cur_geo)                      &            !dyz
        = delta_wrt_ket(3)*cur_geo_pints(:,2,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_geo_pints(:,2,addr_low_zgeo)
      cur_geo_pints(:,9,addr_cur_geo)                                &  !dzz
        = delta_wrt_ket(3)*cur_geo_pints(:,3,addr_cur_geo)           &
        + real(cur_order_geo,REALK)*low_geo_pints(:,3,addr_low_zgeo) &
        + cur_zero_pints(:,addr_cur_geo)
      ! initializes the (base) addresses of HGTOs on ket center
      base_low_hket = 0
      base_cur_hket = 3
      addr_up_hket = 9
      ! other order (>2) HGTOs on ket center
      do order_hket = 2, max_cur_hket
        ! recurrence relation along x-direction
        addr_up_hket = addr_up_hket+1
        addr_cur_hket = base_cur_hket+1
        cur_geo_pints(:,addr_up_hket,addr_cur_geo)                       &
          = delta_wrt_ket(1)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo) &
          + half_nr_expnt*real(order_hket,REALK)                         &
          * cur_geo_pints(:,base_low_hket+1,addr_cur_geo)
        ! recurrence relation along y-direction
        addr_up_hket = addr_up_hket+1
        cur_geo_pints(:,addr_up_hket,addr_cur_geo) &
          = delta_wrt_ket(2)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo)
        addr_low_hket = base_low_hket
        do iket = 1, order_hket
          addr_low_hket = addr_low_hket+1
          addr_cur_hket = addr_cur_hket+1
          addr_up_hket = addr_up_hket+1
          cur_geo_pints(:,addr_up_hket,addr_cur_geo)                       &
            = delta_wrt_ket(2)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo) &
            + half_nr_expnt*real(iket,REALK)                               &
            * cur_geo_pints(:,addr_low_hket,addr_cur_geo)
        end do
        ! recurrence relation along z-direction
        addr_cur_hket = base_cur_hket
        do jket = 0, order_hket
          addr_cur_hket = addr_cur_hket+1
          addr_up_hket = addr_up_hket+1
          cur_geo_pints(:,addr_up_hket,addr_cur_geo)                       &
            = delta_wrt_ket(3)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo) &
            + real(cur_order_geo,REALK)*low_geo_pints(:,addr_cur_hket,addr_low_zgeo)
        end do
        addr_low_hket = base_low_hket
        do iket = 1, order_hket
          do jket = 0, order_hket-iket
            addr_low_hket = addr_low_hket+1
            addr_cur_hket = addr_cur_hket+1
            addr_up_hket = addr_up_hket+1
            cur_geo_pints(:,addr_up_hket,addr_cur_geo)                                 &
              = delta_wrt_ket(3)*cur_geo_pints(:,addr_cur_hket,addr_cur_geo)           &
              + real(cur_order_geo,REALK)*low_geo_pints(:,addr_cur_hket,addr_low_zgeo) &
              + half_nr_expnt*real(iket,REALK)                                         &
              * cur_geo_pints(:,addr_low_hket,addr_cur_geo)
          end do
        end do
        ! updates the base addresses
        base_low_hket = addr_low_hket
        base_cur_hket = addr_cur_hket
      end do
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "sub_delta_hket", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("sub_delta_hket>> ",A,I6)
#endif
  end subroutine sub_delta_hket

  !> \brief assigns the Dirac delta function integrals
  !> \author Bin Gao
  !> \date 2012-03-16
  !> \param dim_hgto_bra is the dimension of HGTOs on bra center
  !> \param num_cur_geo is the number of xyz components of current order geometric derivatives
  !> \param hbra_pints contains the primitive Hermite integrals with zeroth order HGTO on ket center
  !> \param dim_up_hket is the dimension of upper order HGTOs on ket center
  !> \param cur_geo_pints contains the integrals with upper order HGTOs on ket center
  !>        and current order geometric derivatives on Dirac delta function
  !> \param zero_hket indicates if returning zeroth order HGTO on ket center
  !> \param dim_hgto_ket is the dimension of HGTOs on ket center
  !> \return hket_pints contains the primitive Hermite integrals with the specified orders
  !>         of HGTOs on ket center and geometric derivatives on Dirac delta function
  subroutine delta_hket_assign(dim_hgto_bra, num_cur_geo, hbra_pints, &
                               dim_up_hket, cur_geo_pints, zero_hket, &
                               dim_hgto_ket, hket_pints)
    use xkind
    implicit none
    integer, intent(in) :: dim_hgto_bra
    integer, intent(in) :: num_cur_geo
    real(REALK), intent(in) :: hbra_pints(dim_hgto_bra,num_cur_geo)
    integer, intent(in) :: dim_up_hket
    real(REALK), intent(in) :: cur_geo_pints(dim_hgto_bra,dim_up_hket,num_cur_geo)
    logical, intent(in) :: zero_hket
    integer, intent(in) :: dim_hgto_ket
    real(REALK), intent(out) :: hket_pints(dim_hgto_bra,dim_hgto_ket,num_cur_geo)
!f2py intent(hide) :: dim_hgto_bra
!f2py intent(hide) :: num_cur_geo
!f2py intent(in) :: hbra_pints
!f2py intent(hide) :: dim_up_hket
!f2py intent(in) :: cur_geo_pints
!f2py depend(dim_hgto_bra) :: cur_geo_pints
!f2py depend(num_cur_geo) :: cur_geo_pints
!f2py intent(in) :: zero_hket
!f2py intent(in) :: dim_hgto_ket
!f2py intent(out) :: hket_pints
!f2py depend(dim_hgto_bra) :: hket_pints
!f2py depend(dim_hgto_ket) :: hket_pints
!f2py depend(num_cur_geo) :: hket_pints
    integer start_hket     !start address of HGTOs on ket center in temporary integrals
    integer igeo           !incremental recorder over geometric derivatives
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time   
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (zero_hket) then
      do igeo = 1, num_cur_geo
        hket_pints(:,1,igeo) = hbra_pints(:,igeo)
        hket_pints(:,2:dim_hgto_ket,igeo) = cur_geo_pints(:,:,igeo)
      end do
    else
      start_hket = dim_up_hket-dim_hgto_ket+1
      do igeo = 1, num_cur_geo
        hket_pints(:,:,igeo) = cur_geo_pints(:,start_hket:dim_up_hket,igeo)
      end do
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "delta_hket_assign", STDOUT)
#endif
    return
  end subroutine delta_hket_assign
