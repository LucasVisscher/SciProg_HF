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
!!  This file recovers Cartesian multipole moments (using HGTOs) at the origin of
!!  London phase factor.
!!
!!  2012-03-05, Bin Gao
!!  * rewrites to improve efficiency
!!
!!  2011-10-24, Bin Gao
!!  * first version

#include "stdout.h"

  !> \brief recovers Cartesian multipole moments (using HGTOs) at the origin of
  !>        London phase factor
  !> \author Bin Gao
  !> \date 2011-10-24
  !> \param orders_hgto_ket contains the minimum and maximum orders of HGTOs to return
  !> \param orders_mom contains the minimum and maximum orders of Cartesian multipole moments to return
  !> \param coord_ket contains the coordinates of ket center
  !> \param exponent_ket is the exponent of HGTOs of ket center
  !> \param london_origin contains the coordinates of the origin of London phase factor
  !> \param dim_hgto_bra is the dimension of HGTOs of bra center
  !> \param dim_hket_zero is the dimension of HGTOs on ket center with zeroth order
  !>        Cartesian multipole moment
  !> \param num_opt is the number of operators
  !> \param num_geo is the number of geometric derivatives
  !> \param hgto_pints contains the primitive HGTO integrals
  !> \param dim_hgto_ket is the dimension of HGTOs on ket center in returned integrals
  !> \param dim_mom is the dimension of Cartesian multipole moments in returned integrals
  !> \return lmom_pints contains the integrals with Cartesian multipole moments at
  !>         the origin of London phase factor
  subroutine london_mom_hgto(orders_hgto_ket, orders_mom, coord_ket, exponent_ket, &
                             london_origin, dim_hgto_bra, dim_hket_zero, num_opt,  &
                             num_geo, hgto_pints, dim_hgto_ket, dim_mom, lmom_pints)
    use xkind
    implicit none
    integer, intent(in) :: orders_hgto_ket(2)
    integer, intent(in) :: orders_mom(2)
    real(REALK), intent(in) :: coord_ket(3)
    real(REALK), intent(in) :: exponent_ket
    real(REALK), intent(in) :: london_origin(3)
    integer, intent(in) :: dim_hgto_bra
    integer, intent(in) :: dim_hket_zero
    integer, intent(in) :: num_opt
    integer, intent(in) :: num_geo
    real(REALK), intent(in) :: hgto_pints(dim_hgto_bra,dim_hket_zero, &
                                          num_opt,num_geo)
    integer, intent(in) :: dim_hgto_ket
    integer, intent(in) :: dim_mom
    real(REALK), intent(out) :: lmom_pints(dim_hgto_bra,dim_hgto_ket, &
                                           num_opt,dim_mom,num_geo)
!f2py intent(in) :: orders_hgto_ket
!f2py intent(in) :: orders_mom
!f2py intent(in) :: coord_ket
!f2py intent(in) :: exponent_ket
!f2py intent(in) :: london_origin
!f2py intent(hide) :: dim_hgto_bra
!f2py intent(hide) :: dim_hket_zero
!f2py intent(hide) :: num_opt
!f2py intent(hide) :: num_geo
!f2py intent(in) :: hgto_pints
!f2py intent(in) :: dim_hgto_ket
!f2py intent(in) :: dim_mom
!f2py intent(out) :: lmom_pints
!f2py depend(dim_hgto_bra) :: lmom_pints
!f2py depend(dim_hgto_ket) :: lmom_pints
!f2py depend(num_opt) :: lmom_pints
!f2py depend(dim_mom) :: lmom_pints
!f2py depend(num_geo) :: lmom_pints
    real(REALK) half_recip_expnt   !half of the reciprocal of exponent
    real(REALK) ket_wrt_london(3)  !relative coordinates of ket center w.r.t. the origin of London phase factor
    integer min_order_hket         !minimum of current order HGTOs (might be <0)
    integer cur_order_hket(2)      !minimum and maximum of current order HGTOs in recurrence relations
    integer dim_low_hket           !dimension of lower order HGTOs of temporary integrals
    integer dim_up_hket            !dimension of upper order HGTOs of temporary integrals
    integer num_low_mom            !number of lower order Cartesian multipole moments
    integer num_up_mom             !number of upper order Cartesian multipole moments
    integer size_opt_geo           !size of operators and geometric derivatives
    integer size_low_mom           !size of temporary integrals of lower order Cartesian multipole moments
    integer size_up_mom            !size of temporary integrals of upper order Cartesian multipole moments
    integer low_mom_int            !pointer to temporary integrals of lower order Cartesian multipole moments
    integer up_mom_int             !pointer to temporary integrals of upper order Cartesian multipole moments
    integer dim_tmp                !dimension of temporary integrals
    real(REALK), allocatable :: tmp_ints(:,:,:)
                                   !temporary integrals
    integer order_mom              !incremental recorder over orders of Cartesian multipole moments
    integer start_addr_mom         !start address of Cartesian multipole moments in returned integrals
    integer end_addr_mom           !end address of Cartesian multipole moments in returned integrals
    integer start_min_hket         !start address of minimum order HGTOs on ket center
    integer dim_max_hket           !dimension of HGTOs on ket center from zeroth order to maximum order
    integer offset_low_hket        !offset of lower order HGTOs on ket center in temporary integrals
    integer start_low_hket         !start address of HGTOs in temporary integrals
    integer end_low_hket           !end address of HGTOs in temporary integrals
    integer ierr                   !error information
#if defined(XTIME)
    real(REALK) curr_time          !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    select case(orders_mom(2))
    ! returns zeroth order Cartesian multipole moment
    case(0)
      lmom_pints(:,:,:,1,:) = hgto_pints
    ! other order Cartesian multipole moments
    case default
      ! sets the start address of minimum order HGTOs on ket center
      start_min_hket = orders_hgto_ket(1)*(orders_hgto_ket(1)+1) &
                     * (orders_hgto_ket(1)+2)/6+1
      ! sets the dimension of HGTOs on ket center from zeroth order to maximum order
      dim_max_hket = (orders_hgto_ket(2)+1)*(orders_hgto_ket(2)+2) &
                   * (orders_hgto_ket(2)+3)/6
      ! sets the minimum order of HGTOs with zeroth order Cartesian multipole moment
      min_order_hket = orders_hgto_ket(1)-orders_mom(2)
      cur_order_hket(1) = max(0,min_order_hket)
      ! assigns the zeroth order Cartesian multipole moment
      if (orders_mom(1)==0) then
        ! sets the start and end addresses of lower order HGTOs on ket center to return
        if (cur_order_hket(1)==0) then
          start_low_hket = start_min_hket
          end_low_hket = dim_max_hket
        else
          offset_low_hket = cur_order_hket(1)*(cur_order_hket(1)+1) &
                          * (cur_order_hket(1)+2)/6
          start_low_hket = start_min_hket-offset_low_hket
          end_low_hket = dim_max_hket-offset_low_hket
        end if
        lmom_pints(:,:,:,1,:) &
          = hgto_pints(:,start_low_hket:end_low_hket,:,:)
        ! initializes the start address of Cartesian multipole moments in returned integrals
        start_addr_mom = 1
      else
        start_addr_mom = 0
      end if
      ! half of the reciprocal of exponent
      half_recip_expnt = 0.5_REALK/exponent_ket
      ! relative coordinates of ket center w.r.t. the origin of London phase factor
      ket_wrt_london = coord_ket-london_origin
      ! sets the size of operators and geometric derivatives
      size_opt_geo = num_opt*num_geo
      ! allocates memory for temporary integrals
      call dim_london_mom_hgto(orders_hgto_ket, orders_mom(2), dim_tmp)
      allocate(tmp_ints(dim_hgto_bra,dim_tmp*size_opt_geo,2), stat=ierr)
      if (ierr/=0)                                                        &
        call error_stop("london_mom_hgto", "failed to allocate tmp_ints", &
                        dim_hgto_bra*dim_tmp*size_opt_geo*2)
      ! sets the minimum and maximum orders of HGTOs with first order Cartesian multipole moments
      min_order_hket = min_order_hket+1
      cur_order_hket(1) = max(0,min_order_hket)
      cur_order_hket(2) = orders_hgto_ket(2)+orders_mom(2)-1
      ! sets the dimensions of HGTOs and Cartesian multipole moments
      dim_up_hket = dim_hket_zero
      if (cur_order_hket(1)==0) then
        offset_low_hket = 0
        dim_low_hket = (cur_order_hket(2)+1)*(cur_order_hket(2)+2) &
                     * (cur_order_hket(2)+3)/6
      else
        offset_low_hket = cur_order_hket(1)*(cur_order_hket(1)+1) &
                        * (cur_order_hket(1)+2)/6
        dim_low_hket = (cur_order_hket(2)+1)*(cur_order_hket(2)+2) &
                     * (cur_order_hket(2)+3)/6-offset_low_hket
      end if
      num_low_mom = 1
      num_up_mom = 3
      ! sets the size of temporary integrals of first order Cartesian multipole moments
      size_up_mom = dim_low_hket*num_up_mom*size_opt_geo
      ! gets the first order Cartesian multipole moments integrals
      call sub_london_mom_hgto(0, cur_order_hket, ket_wrt_london, half_recip_expnt, &
                               dim_hgto_bra, dim_up_hket, num_opt, num_low_mom,     &
                               num_geo, hgto_pints, dim_low_hket, num_up_mom,       &
                               tmp_ints(:,1:size_up_mom,1))
      ! initializes the pointers of Cartesian multipole moments
      low_mom_int = 2
      up_mom_int = 1
      ! loops over the orders of Cartesian multipole moments till the minimum returned order
      do order_mom = 1, orders_mom(1)-1, 1
        ! sets the minimum and maximum orders of HGTOs on ket center
        min_order_hket = min_order_hket+1
        cur_order_hket(1) = max(0,min_order_hket)
        cur_order_hket(2) = cur_order_hket(2)-1
        ! sets the dimensions of HGTOs and Cartesian multipole moments
        dim_up_hket = dim_low_hket
        if (cur_order_hket(1)==0) then
          offset_low_hket = 0
          dim_low_hket = (cur_order_hket(2)+1)*(cur_order_hket(2)+2) &
                       * (cur_order_hket(2)+3)/6
        else
          offset_low_hket = cur_order_hket(1)*(cur_order_hket(1)+1) &
                          * (cur_order_hket(1)+2)/6
          dim_low_hket = (cur_order_hket(2)+1)*(cur_order_hket(2)+2) &
                       * (cur_order_hket(2)+3)/6-offset_low_hket
        end if
        num_low_mom = num_up_mom
        num_up_mom = num_up_mom+(order_mom+2)  !=(order_mom+2)*(order_mom+3)/2
        ! updates the sizes of temporary integrals
        size_low_mom = size_up_mom
        size_up_mom = dim_low_hket*num_up_mom*size_opt_geo
        ! switches the pointers
        low_mom_int = 3-low_mom_int
        up_mom_int = 3-up_mom_int
        ! gets the order \var(order_mom)+1 Cartesian multipole moments integrals
        call sub_london_mom_hgto(order_mom, cur_order_hket, ket_wrt_london,   &
                                 half_recip_expnt, dim_hgto_bra, dim_up_hket, &
                                 num_opt, num_low_mom, num_geo,               &
                                 tmp_ints(:,1:size_low_mom,low_mom_int),      &
                                 dim_low_hket, num_up_mom,                    &
                                 tmp_ints(:,1:size_up_mom,up_mom_int))
      end do
      ! assigns the integrals
      if (cur_order_hket(1)==0) then
        start_low_hket = start_min_hket
        end_low_hket = dim_max_hket
      else
        start_low_hket = start_min_hket-offset_low_hket
        end_low_hket = dim_max_hket-offset_low_hket
      end if
      end_addr_mom = start_addr_mom+num_up_mom
      start_addr_mom = start_addr_mom+1
      call london_mom_hgto_assign(start_low_hket, end_low_hket, dim_hgto_bra, &
                                  dim_low_hket, num_opt, num_up_mom, num_geo, &
                                  tmp_ints(:,1:size_up_mom,up_mom_int),       &
                                  start_addr_mom, end_addr_mom, dim_hgto_ket, &
                                  dim_mom, lmom_pints)
      ! updates the start address of Cartesian multipole moments in returned integrals
      start_addr_mom = end_addr_mom
      ! loops over other returned orders of Cartesian multipole moments
      do order_mom = max(orders_mom(1),1), orders_mom(2)-1, 1
        ! sets the minimum and maximum orders of HGTOs on ket center
        min_order_hket = min_order_hket+1
        cur_order_hket(1) = max(0,min_order_hket)
        cur_order_hket(2) = cur_order_hket(2)-1
        ! sets the dimensions of HGTOs and Cartesian multipole moments
        dim_up_hket = dim_low_hket
        if (cur_order_hket(1)==0) then
          offset_low_hket = 0
          dim_low_hket = (cur_order_hket(2)+1)*(cur_order_hket(2)+2) &
                       * (cur_order_hket(2)+3)/6
        else
          offset_low_hket = cur_order_hket(1)*(cur_order_hket(1)+1) &
                          * (cur_order_hket(1)+2)/6
          dim_low_hket = (cur_order_hket(2)+1)*(cur_order_hket(2)+2) &
                       * (cur_order_hket(2)+3)/6-offset_low_hket
        end if
        num_low_mom = num_up_mom
        num_up_mom = num_up_mom+(order_mom+2)  !=(order_mom+2)*(order_mom+3)/2
        ! updates the sizes of temporary integrals
        size_low_mom = size_up_mom
        size_up_mom = dim_low_hket*num_up_mom*size_opt_geo
        ! switches the pointers
        low_mom_int = 3-low_mom_int
        up_mom_int = 3-up_mom_int
        ! gets the order \var(order_mom)+1 Cartesian multipole moments integrals
        call sub_london_mom_hgto(order_mom, cur_order_hket, ket_wrt_london,   &
                                 half_recip_expnt, dim_hgto_bra, dim_up_hket, &
                                 num_opt, num_low_mom, num_geo,               &
                                 tmp_ints(:,1:size_low_mom,low_mom_int),      &
                                 dim_low_hket, num_up_mom,                    &
                                 tmp_ints(:,1:size_up_mom,up_mom_int))
        ! assigns the integrals
        start_low_hket = start_min_hket-offset_low_hket
        end_low_hket = dim_max_hket-offset_low_hket
        end_addr_mom = start_addr_mom+num_up_mom
        start_addr_mom = start_addr_mom+1
        call london_mom_hgto_assign(start_low_hket, end_low_hket, dim_hgto_bra, &
                                    dim_low_hket, num_opt, num_up_mom, num_geo, &
                                    tmp_ints(:,1:size_up_mom,up_mom_int),       &
                                    start_addr_mom, end_addr_mom, dim_hgto_ket, &
                                    dim_mom, lmom_pints)
        ! updates the start address of Cartesian multipole moments in returned integrals
        start_addr_mom = end_addr_mom
      end do
      deallocate(tmp_ints)
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "london_mom_hgto", STDOUT)
#endif
    return
  end subroutine london_mom_hgto

  !> \brief gets the maximum dimension of temporary integrals used in recurrence relations
  !> \author Bin Gao
  !> \date 2012-03-05
  !> \param orders_hgto_ket contains the minimum and maximum orders of HGTOs to return
  !> \param max_order_mom is the maximum order of Cartesian multipole moments
  !> \return dim_ints is the maximum dimension of temporary integrals
  subroutine dim_london_mom_hgto(orders_hgto_ket, max_order_mom, dim_ints)
    use xkind
    implicit none
    integer, intent(in) :: orders_hgto_ket(2)
    integer, intent(in) :: max_order_mom
    integer, intent(out) :: dim_ints
!f2py intent(in) :: orders_hgto_ket
!f2py intent(in) :: max_order_mom
!f2py intent(out) :: dim_ints
    integer min_hgto_ket       !minimum order of HGTOs on ket center, might be <0
    integer range_hgto_ket(2)  !range of orders of HGTOs on ket center in recurrence relations
    integer order_mom          !incremental recorder over orders of Cartesian multipole moments
    integer dim_hgto_ket       !dimension of HGTOs on ket center
    integer num_mom            !number of Cartesian multipole moments
    integer dim_tmp            !temporary result of dimension
#if defined(XTIME)
    real(REALK) curr_time      !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! initializes the return value
    dim_ints = 0
    ! initializes the minimum and maximum order of HGTOs
    min_hgto_ket = orders_hgto_ket(1)-max_order_mom
    range_hgto_ket(2) = orders_hgto_ket(2)+max_order_mom
    ! initializes the number of Cartesian multipole moments
    num_mom = 1
    ! loops over the orders of Cartesian multipole moments
    do order_mom = 1, max_order_mom
      num_mom = num_mom+order_mom+1  !=(order_mom+1)*(order_mom+2)/2
      ! updates the minimum and maximum orders of HGTOs
      min_hgto_ket = min_hgto_ket+1
      range_hgto_ket(1) = max(0,min_hgto_ket)
      range_hgto_ket(2) = range_hgto_ket(2)-1
      ! sets the dimension of HGTOs and Cartesian multipole moments
      if (range_hgto_ket(1)==0) then
        dim_hgto_ket = (range_hgto_ket(2)+1)*(range_hgto_ket(2)+2) &
                     * (range_hgto_ket(2)+3)/6
      else
        dim_hgto_ket = ((range_hgto_ket(2)+1)*(range_hgto_ket(2)+2) &
                     *  (range_hgto_ket(2)+3)                       &
                     -  range_hgto_ket(1)*(range_hgto_ket(1)+1)     &
                     *  (range_hgto_ket(1)+2))/6
      end if
      ! updates the maximum dimension
      dim_tmp = dim_hgto_ket*num_mom
      if (dim_tmp>dim_ints) dim_ints = dim_tmp 
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "dim_london_mom_hgto", STDOUT)
#endif
    return
  end subroutine dim_london_mom_hgto

  !> \brief sub-recurrence relations by recovering upper order Cartesian multipole moments
  !> \author Bin Gao
  !> \date 2011-10-24
  !> \param cur_order_mom is the current order of Cartesian multipole moments
  !> \param cur_order_hket contains the minimum and maximum current orders of HGTOs
  !> \param ket_wrt_london contains relative coordinates of ket center w.r.t. the origin
  !>        of London phase factor
  !> \param half_recip_expnt is the half of the reciprocal of exponent on ket center
  !> \param dim_hgto_bra is the dimension of HGTOs on bra center
  !> \param dim_up_hket is the dimension of HGTOs with lower order Cartesian multipole moments
  !> \param num_low_mom is the number of lower order Cartesian multipole moments
  !> \param num_geo is the number of geometric derivatives
  !> \param low_lmom_pints contains the integrals with lower order Cartesian multipole moments
  !> \param dim_low_hket is the dimension of HGTOs with upper order Cartesian multipole moments
  !> \param num_up_mom is the number upper order Cartesian multipole moments
  !> \return up_lmom_pints contains the integrals with upper order Cartesian multipole moments
  subroutine sub_london_mom_hgto(cur_order_mom, cur_order_hket, ket_wrt_london, &
                                 half_recip_expnt, dim_hgto_bra, dim_up_hket,   &
                                 num_opt, num_low_mom, num_geo, low_lmom_pints, &
                                 dim_low_hket, num_up_mom, up_lmom_pints)
    use xkind
    implicit none
    integer, intent(in) :: cur_order_mom
    integer, intent(in) :: cur_order_hket(2)
    real(REALK), intent(in) :: ket_wrt_london(3)
    real(REALK), intent(in) :: half_recip_expnt
    integer, intent(in) :: dim_hgto_bra
    integer, intent(in) :: dim_up_hket
    integer, intent(in) :: num_opt
    integer, intent(in) :: num_low_mom
    integer, intent(in) :: num_geo
    real(REALK), intent(in) :: low_lmom_pints(dim_hgto_bra,dim_up_hket, &
                                              num_opt,num_low_mom,num_geo)
    integer, intent(in) :: dim_low_hket
    integer, intent(in) :: num_up_mom
    real(REALK), intent(out) :: up_lmom_pints(dim_hgto_bra,dim_low_hket, &
                                              num_opt,num_up_mom,num_geo)
!f2py intent(in) :: cur_order_mom
!f2py intent(in) :: cur_order_hket
!f2py intent(in) :: ket_wrt_london
!f2py intent(in) :: half_recip_expnt
!f2py intent(hide) :: dim_hgto_bra
!f2py intent(hide) :: dim_up_hket
!f2py intent(hide) :: num_opt
!f2py intent(hide) :: num_low_mom
!f2py intent(hide) :: num_geo
!f2py intent(in) :: low_lmom_pints
!f2py intent(in) :: dim_low_hket
!f2py intent(in) :: num_up_mom
!f2py intent(out) :: up_lmom_pints
!f2py depend(dim_hgto_bra) :: up_lmom_pints
!f2py depend(dim_low_hket) :: up_lmom_pints
!f2py depend(num_opt) :: up_lmom_pints
!f2py depend(num_up_mom) :: up_lmom_pints
!f2py depend(num_geo) :: up_lmom_pints
    integer min_cur_hket   !current minimum order of HGTOs on ket center
    logical zero_cur_hket  !if the minimum order of HGTOs on ket center is zeroth order
    integer igeo           !incremental recorder over geometric derivatives
    integer imom, jmom     !incremental recorders over Cartesian multipole moments
    integer addr_up_mom    !address of upper order Cartesian multipole moments
    integer addr_low_mom   !addresses of lower order Cartesian multipole moments
    integer addr_low_ymom
    integer addr_cur_hket  !address of current order HGTOs in integrals with upper order
                           !Cartesian multipole moments
    integer addr_hket      !address of current order HGTOs in integrals with lower order
                           !Cartesian multipole moments
    integer addr_low_hket  !address of lower order HGTOs in integrals with lower order
                           !Cartesian multipole moments
    integer addr_up_hket   !address of upper order HGTOs in integrals with lower order
                           !Cartesian multipole moments
    integer iket, jket     !incremental recorder over HGTOs
    integer order_hket     !order of HGTOs
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (cur_order_hket(1)==0) then
      zero_cur_hket = .true.
      min_cur_hket = 1
    else
      zero_cur_hket = .false.
      min_cur_hket = cur_order_hket(1)
    end if
    ! loops over geometric derivatives
    do igeo = 1, num_geo
      addr_up_mom = 1
      addr_low_mom = 1
      ! (1) x-direction along recurrence relation of Cartesian multipole moments
#if defined(DEBUG)
      write(STDOUT,100) "x-direction"
#endif
      if (zero_cur_hket) then
#if defined(DEBUG)
        write(STDOUT,100) "orders:", 0, cur_order_mom
        write(STDOUT,110) 1, addr_up_mom
        write(STDOUT,110) 1, addr_low_mom
        write(STDOUT,110) 2, addr_low_mom
        write(STDOUT,100) "------------------------------------"
#endif
        up_lmom_pints(:,1,:,addr_up_mom,igeo)       &
          = ket_wrt_london(1)                       &
          * low_lmom_pints(:,1,:,addr_low_mom,igeo) &
          + low_lmom_pints(:,2,:,addr_low_mom,igeo)
        addr_cur_hket = 1  !base address of p-shell HGTOs
        addr_hket = 1      !base address of p-shell HGTOs
        addr_up_hket = 4   !base address of d-shell HGTOs
      else
        addr_cur_hket = 0
        addr_hket = min_cur_hket*(min_cur_hket+1)/2
        addr_up_hket = addr_hket+(min_cur_hket+1)*(min_cur_hket+2)/2
      end if
      addr_low_hket = 0
      ! loops over other order of HGTOs, starting from d-shell
      do order_hket = min_cur_hket, cur_order_hket(2)
#if defined(DEBUG)
        write(STDOUT,100) "orders:", order_hket, cur_order_mom
#endif
        do iket = order_hket, 0, -1
          do jket = iket, 1, -1
            addr_cur_hket = addr_cur_hket+1
            addr_hket = addr_hket+1
            addr_low_hket = addr_low_hket+1
            addr_up_hket = addr_up_hket+1
#if defined(DEBUG)
            write(STDOUT,110) addr_cur_hket, addr_up_mom
            write(STDOUT,110) addr_hket, addr_low_mom
            write(STDOUT,110) addr_low_hket, addr_low_mom
            write(STDOUT,110) addr_up_hket, addr_low_mom
            write(STDOUT,100) "------------------------------------"
#endif
            up_lmom_pints(:,addr_cur_hket,:,addr_up_mom,igeo)       &
              = ket_wrt_london(1)                                   &
              * low_lmom_pints(:,addr_hket,:,addr_low_mom,igeo)     &
              + half_recip_expnt*real(jket,REALK)                   &
              * low_lmom_pints(:,addr_low_hket,:,addr_low_mom,igeo) &
              + low_lmom_pints(:,addr_up_hket,:,addr_low_mom,igeo)
          end do
          addr_cur_hket = addr_cur_hket+1
          addr_hket = addr_hket+1
          addr_up_hket = addr_up_hket+1
#if defined(DEBUG)
          write(STDOUT,110) addr_cur_hket, addr_up_mom
          write(STDOUT,110) addr_hket, addr_low_mom
          write(STDOUT,110) addr_up_hket, addr_low_mom
          write(STDOUT,100) "------------------------------------"
#endif
          up_lmom_pints(:,addr_cur_hket,:,addr_up_mom,igeo)       &
            = ket_wrt_london(1)                                   &
            * low_lmom_pints(:,addr_hket,:,addr_low_mom,igeo)     &
            + low_lmom_pints(:,addr_up_hket,:,addr_low_mom,igeo)
          addr_up_hket = addr_up_hket+1
        end do
        addr_up_hket = addr_up_hket+1
      end do
      ! (2) y-direction along recurrence relation of Cartesian multipole moments
#if defined(DEBUG)
      write(STDOUT,100) "y-direction"
#endif
      do imom = 0, cur_order_mom
        addr_up_mom = addr_up_mom+1
        addr_low_ymom = addr_low_mom+imom
        if (zero_cur_hket) then
#if defined(DEBUG)
          write(STDOUT,100) "orders:", 0, cur_order_mom
          write(STDOUT,110) 1, addr_up_mom
          write(STDOUT,110) 1, addr_low_ymom
          write(STDOUT,110) 3, addr_low_ymom
          write(STDOUT,100) "------------------------------------"
#endif
          up_lmom_pints(:,1,:,addr_up_mom,igeo)        &
            = ket_wrt_london(2)                        &
            * low_lmom_pints(:,1,:,addr_low_ymom,igeo) &
            + low_lmom_pints(:,3,:,addr_low_ymom,igeo)
          addr_cur_hket = 1  !base address of p-shell HGTOs
          addr_hket = 1      !base address of p-shell HGTOs
          addr_up_hket = 4   !base address of d-shell HGTOs
        else
          addr_cur_hket = 0
          addr_hket = min_cur_hket*(min_cur_hket+1)/2
          addr_up_hket = addr_hket+(min_cur_hket+1)*(min_cur_hket+2)/2
        end if
        addr_low_hket = 0
        ! loops over other order of HGTOs, starting from d-shell
        do order_hket = min_cur_hket, cur_order_hket(2)
#if defined(DEBUG)
          write(STDOUT,100) "orders:", order_hket, cur_order_mom
#endif
          addr_up_hket = addr_up_hket+1
          do iket = order_hket, 0, -1
            addr_cur_hket = addr_cur_hket+1
            addr_hket = addr_hket+1
            addr_up_hket = addr_up_hket+1
#if defined(DEBUG)
            write(STDOUT,110) addr_cur_hket, addr_up_mom
            write(STDOUT,110) addr_hket, addr_low_ymom
            write(STDOUT,110) addr_up_hket, addr_low_ymom
            write(STDOUT,100) "------------------------------------"
#endif
            up_lmom_pints(:,addr_cur_hket,:,addr_up_mom,igeo)        &
              = ket_wrt_london(2)                                    &
              * low_lmom_pints(:,addr_hket,:,addr_low_ymom,igeo)     &
              + low_lmom_pints(:,addr_up_hket,:,addr_low_ymom,igeo)
            do jket = 1, iket
              addr_cur_hket = addr_cur_hket+1
              addr_hket = addr_hket+1
              addr_low_hket = addr_low_hket+1
              addr_up_hket = addr_up_hket+1
#if defined(DEBUG)
              write(STDOUT,110) addr_cur_hket, addr_up_mom
              write(STDOUT,110) addr_hket, addr_low_ymom
              write(STDOUT,110) addr_low_hket, addr_low_ymom
              write(STDOUT,110) addr_up_hket, addr_low_ymom
              write(STDOUT,100) "------------------------------------"
#endif
              up_lmom_pints(:,addr_cur_hket,:,addr_up_mom,igeo)        &
                = ket_wrt_london(2)                                    &
                * low_lmom_pints(:,addr_hket,:,addr_low_ymom,igeo)     &
                + half_recip_expnt*real(jket,REALK)                    &
                * low_lmom_pints(:,addr_low_hket,:,addr_low_ymom,igeo) &
                + low_lmom_pints(:,addr_up_hket,:,addr_low_ymom,igeo)
            end do
            addr_up_hket = addr_up_hket+1
          end do
        end do
      end do
      ! (3) z-direction along recurrence relation of Cartesian multipole moments
#if defined(DEBUG)
      write(STDOUT,100) "z-direction"
#endif
      addr_low_mom = addr_low_mom-1
      do imom = 0, cur_order_mom
        do jmom = 0, cur_order_mom-imom
          addr_up_mom = addr_up_mom+1
          addr_low_mom = addr_low_mom+1
          if (zero_cur_hket) then
#if defined(DEBUG)
            write(STDOUT,100) "orders:", 0, cur_order_mom
            write(STDOUT,110) 1, addr_up_mom
            write(STDOUT,110) 1, addr_low_mom
            write(STDOUT,110) 4, addr_low_mom
            write(STDOUT,100) "------------------------------------"
#endif
            up_lmom_pints(:,1,:,addr_up_mom,igeo)       &
              = ket_wrt_london(3)                       &
              * low_lmom_pints(:,1,:,addr_low_mom,igeo) &
              + low_lmom_pints(:,4,:,addr_low_mom,igeo)
            addr_cur_hket = 1  !base address of p-shell HGTOs
            addr_hket = 1      !base address of p-shell HGTOs
            addr_up_hket = 4   !base address of d-shell HGTOs
          else
            addr_cur_hket = 0
            addr_hket = min_cur_hket*(min_cur_hket+1)/2
            addr_up_hket = addr_hket+(min_cur_hket+1)*(min_cur_hket+2)/2
          end if
          addr_low_hket = 0
          ! loops over other order of HGTOs, starting from d-shell
          do order_hket = min_cur_hket, cur_order_hket(2)
#if defined(DEBUG)
            write(STDOUT,100) "orders:", order_hket, cur_order_mom
#endif
            addr_up_hket = addr_up_hket+order_hket+2
            do jket = 0, order_hket
              addr_cur_hket = addr_cur_hket+1
              addr_hket = addr_hket+1
              addr_up_hket = addr_up_hket+1
#if defined(DEBUG)
              write(STDOUT,110) addr_cur_hket, addr_up_mom
              write(STDOUT,110) addr_hket, addr_low_mom
              write(STDOUT,110) addr_up_hket, addr_low_mom
              write(STDOUT,100) "------------------------------------"
#endif
              up_lmom_pints(:,addr_cur_hket,:,addr_up_mom,igeo)       &
                = ket_wrt_london(3)                                   &
                * low_lmom_pints(:,addr_hket,:,addr_low_mom,igeo)     &
                + low_lmom_pints(:,addr_up_hket,:,addr_low_mom,igeo)
            end do
            do iket = 1, order_hket
              do jket = iket, order_hket
                addr_cur_hket = addr_cur_hket+1
                addr_hket = addr_hket+1
                addr_low_hket = addr_low_hket+1
                addr_up_hket = addr_up_hket+1
#if defined(DEBUG)
                write(STDOUT,110) addr_cur_hket, addr_up_mom
                write(STDOUT,110) addr_hket, addr_low_mom
                write(STDOUT,110) addr_low_hket, addr_low_mom
                write(STDOUT,110) addr_up_hket, addr_low_mom
                write(STDOUT,100) "------------------------------------"
#endif
                up_lmom_pints(:,addr_cur_hket,:,addr_up_mom,igeo)       &
                  = ket_wrt_london(3)                                   &
                  * low_lmom_pints(:,addr_hket,:,addr_low_mom,igeo)     &
                  + half_recip_expnt*real(iket,REALK)                   &
                  * low_lmom_pints(:,addr_low_hket,:,addr_low_mom,igeo) &
                  + low_lmom_pints(:,addr_up_hket,:,addr_low_mom,igeo)
              end do
            end do
          end do
        end do
      end do
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "sub_london_mom_hgto", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("sub_london_mom_hgto>> ",A,2I6)
110 format("sub_london_mom_hgto>> ","HGTO",I8,4X,"MOM",I8)
#endif
  end subroutine sub_london_mom_hgto

  !> \brief assigns the integrals with Cartesian multipole moments at the origin
  !>        of London phase factor
  !> \author Bin Gao
  !> \date 2012-03-05
  !> \param start_low_hket is the start address of HGTOs in temporary integrals
  !> \param end_low_hket is the end address of HGTOs in temporary integrals
  !> \param dim_hgto_bra is the dimension of HGTOs of bra center
  !> \param dim_low_hket is the dimension of lower order HGTOs of temporary integrals
  !> \param num_opt is the number of operators
  !> \param num_up_mom is the number of upper order Cartesian multipole moments
  !> \param num_geo is the number of geometric derivatives
  !> \param recur_pints contains the temporary integrals from recurrence relations
  !> \param start_addr_mom is the start address of Cartesian multipole moments in returned integrals
  !> \param end_addr_mom is the end address of Cartesian multipole moments in returned integrals
  !> \param dim_hgto_ket is the dimension of HGTOs on ket center in returned integrals
  !> \param dim_mom is the dimension of Cartesian multipole moments in returned integrals
  !> \return lmom_pints contains the integrals with Cartesian multipole moments at
  !>         the origin of London phase factor
  subroutine london_mom_hgto_assign(start_low_hket, end_low_hket, dim_hgto_bra, &
                                    dim_low_hket, num_opt, num_up_mom, num_geo, &
                                    recur_pints, start_addr_mom, end_addr_mom,  &
                                    dim_hgto_ket, dim_mom, lmom_pints)
    use xkind
    implicit none
    integer, intent(in) :: start_low_hket
    integer, intent(in) :: end_low_hket
    integer, intent(in) :: dim_hgto_bra
    integer, intent(in) :: dim_low_hket
    integer, intent(in) :: num_opt
    integer, intent(in) :: num_up_mom
    integer, intent(in) :: num_geo
    real(REALK), intent(in) :: recur_pints(dim_hgto_bra,dim_low_hket, &
                                           num_opt,num_up_mom,num_geo)
    integer, intent(in) :: start_addr_mom
    integer, intent(in) :: end_addr_mom
    integer, intent(in) :: dim_hgto_ket
    integer, intent(in) :: dim_mom
    real(REALK), intent(inout) :: lmom_pints(dim_hgto_bra,dim_hgto_ket, &
                                             num_opt,dim_mom,num_geo)
!f2py intent(in) :: start_low_hket
!f2py intent(in) :: end_low_hket
!f2py intent(hide) :: dim_hgto_bra
!f2py intent(hide) :: dim_low_hket
!f2py intent(hide) :: num_opt
!f2py intent(hide) :: num_up_mom
!f2py intent(hide) :: num_geo
!f2py intent(in) :: recur_pints
!f2py intent(in) :: start_addr_mom
!f2py intent(in) :: end_addr_mom
!f2py intent(hide) :: dim_hgto_ket
!f2py intent(hide) :: dim_mom
!f2py intent(inout) :: lmom_pints
!f2py depend(dim_hgto_bra) :: lmom_pints
!f2py depend(num_opt) :: lmom_pints
!f2py depend(num_geo) :: lmom_pints
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    lmom_pints(:,:,:,start_addr_mom:end_addr_mom,:) &
      = recur_pints(:,start_low_hket:end_low_hket,:,:,:)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "london_mom_hgto_assign", STDOUT)
#endif
    return
  end subroutine london_mom_hgto_assign
