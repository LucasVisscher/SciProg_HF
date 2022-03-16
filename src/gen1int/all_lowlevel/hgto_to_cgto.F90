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
!!  This file transforms primitive HGTOs to CGTOs on bra or ket center.
!!
!!  2012-03-04, Bin Gao
!!  * rewrites to improve the efficiency
!!
!!  2011-06-22, Bin Gao
!!  * rewrites for efficiency, only for integrals with non-London atomic orbitals
!!
!!  2011-03-29, Bin Gao
!!  * first version

#include "stdout.h"

  !> \brief transforms primitive HGTOs to CGTOs on ket center
  !> \author Bin Gao
  !> \date 2011-03-29
  !> \param angular_ket is the angular number to return
  !> \param order_geo_ket is the order of geometric derivatives to return
  !> \param exponent_ket is the exponent of primitive GTOs of ket center
  !> \param num_cgto_bra is the number of CGTOs of bra center
  !> \param dim_hgto_ket is the dimension of HGTOs of ket center
  !> \param num_geo_bra is the number of geometric derivatives of bra center
  !> \param num_opt is the number of operators
  !> \param hgto_pints contains the primitive HGTO integrals
  !> \param num_cgto_ket is the number of CGTOs of ket center,
  !>        equals to \f$(\var(angular_ket)+1)(\var(angular_ket)+2)/2\f$
  !> \param num_geo_ket is the number of geometric derivatives of ket center,
  !>        equals to \f$(\var(order_geo_ket)+1)(\var(order_geo_ket)+2)/2\f$
  !> \return cgto_pints contains the primitive CGTO integrals
  subroutine hgto_to_cgto(angular_ket, order_geo_ket, exponent_ket, num_cgto_bra, &
                          dim_hgto_ket, num_geo_bra, num_opt, hgto_pints,         &
                          num_cgto_ket, num_geo_ket, cgto_pints)
    use xkind
    implicit none
    integer, intent(in) :: angular_ket
    integer, intent(in) :: order_geo_ket
    real(REALK), intent(in) :: exponent_ket
    integer, intent(in) :: num_cgto_bra
    integer, intent(in) :: dim_hgto_ket
    integer, intent(in) :: num_geo_bra
    integer, intent(in) :: num_opt
    real(REALK), intent(inout) :: hgto_pints(num_cgto_bra,dim_hgto_ket, &
                                             num_geo_bra,num_opt)
    integer, intent(in) :: num_cgto_ket
    integer, intent(in) :: num_geo_ket
    real(REALK), intent(out) :: cgto_pints(num_cgto_bra,num_cgto_ket, &
                                           num_geo_bra,num_geo_ket,num_opt)
!f2py intent(in) :: angular_ket
!f2py intent(in) :: order_geo_ket
!f2py intent(in) :: exponent_ket
!f2py intent(hide) :: num_cgto_bra
!f2py intent(hide) :: dim_hgto_ket
!f2py intent(hide) :: num_geo_bra
!f2py intent(hide) :: num_opt
!f2py intent(inout) :: hgto_pints
!f2py intent(in) :: num_cgto_ket
!f2py intent(in) :: num_geo_ket
!f2py intent(out) :: cgto_pints
!f2py depend(num_cgto_bra) :: cgto_pints
!f2py depend(num_cgto_ket) :: cgto_pints
!f2py depend(num_geo_bra) :: cgto_pints
!f2py depend(num_geo_ket) :: cgto_pints
!f2py depend(num_opt) :: cgto_pints
    integer max_low_hgto            !maximum order of lower order HGTOs
    real(REALK) twice_expnt         !twice of the exponent
    real(REALK) half_recip_expnt    !half of the reciprocal of exponent
    real(REALK) square_half_rexpnt  !square of \var(half_recip_expnt)
    real(REALK) prefact_hgto        !prefactor when transforming HGTOs to CGTOs
    integer max_low_cgto            !maximum order of lower order CGTOs
    integer num_low_hgto            !number of lower order HGTOs
    integer num_up_hgto             !number of higher order HGTOs
    integer start_up_hgto           !start address of upper or current order HGTOs
    integer end_up_hgto             !end address of upper or current order HGTOs
    integer dim_cgto(2)             !dimensions of CGTOs of temporary integrals
    integer size_geo_opt            !size of geometric derivatives on bra center and other operators
    integer size_low_cgto           !size of temporary integrals of lower order CGTOs
    integer size_up_cgto            !size of temporary integrals of upper order CGTOs
    integer low_cgto_int            !pointer to temporary integrals of lower order CGTOs
    integer up_cgto_int             !pointer to temporary integrals of upper order CGTOs
    integer dim_tmp                 !dimension of temporary integrals
    real(REALK), allocatable :: tmp_ints(:,:,:)
                                    !temporary integrals
    logical recur_p_cgto            !if recurrence relation starts from p-shell
    integer offset_cgto(2)          !offsets of CGTOs on ket center
    integer low_order_cgto(2)       !minimum and maximum orders of lower order CGTOs
    integer order_hgto              !incremental recorder over orders of HGTOs
    integer iopt                    !incremental recorder over operators
    integer igeo                    !incremental recorder over geometric derivatives
    integer ierr                    !error information
#if defined(XTIME)
    real(REALK) curr_time           !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    max_low_hgto = order_geo_ket+angular_ket  !maximum order of lower order HGTOs
    twice_expnt = exponent_ket+exponent_ket   !twice of the exponent
    prefact_hgto = twice_expnt**max_low_hgto  !prefactor for the maximum order HGTOs
    select case(angular_ket)
    ! s-shell CGTO returned, only multiplies the prefactor
    case(0)
      do iopt = 1, num_opt
        do igeo = 1, num_geo_bra
          cgto_pints(:,1,igeo,:,iopt) = prefact_hgto*hgto_pints(:,:,igeo,iopt)
        end do
      end do
    ! p-shell CGTOs returned
    case(1)
      hgto_pints = prefact_hgto*hgto_pints      !multiplies the prefactor
      max_low_hgto = max_low_hgto-1             !updates the maximum order of HGTOs
      half_recip_expnt = 1.0_REALK/twice_expnt  !1/(2*\var(exponent_ket)) for recurrence relations
      ! recovers p-shell CGTOs
      call hgto_to_cgto_p(max_low_hgto, half_recip_expnt, num_cgto_bra,   &
                          dim_hgto_ket, num_geo_bra, num_opt, hgto_pints, &
                          num_cgto_ket, num_geo_ket, cgto_pints)
    ! d-shell CGTOs returned
    case(2)
      prefact_hgto = twice_expnt**max_low_hgto           !prefactor for the maximum order HGTOs
      end_up_hgto = dim_hgto_ket                         !end address of upper order HGTOs
      num_up_hgto = (max_low_hgto+1)*(max_low_hgto+2)/2  !number of upper order HGTOs
      start_up_hgto = end_up_hgto-num_up_hgto+1          !start address of upper order HGTOs
#if defined(DEBUG)
      write(STDOUT,100) "HGTO-s-CGTO/order/start/end:", &
                        max_low_hgto, start_up_hgto, end_up_hgto
#endif
      ! multiplies the prefactor for order \var(max_low_hgto) HGTOs
      hgto_pints(:,start_up_hgto:end_up_hgto,:,:) &
        = prefact_hgto*hgto_pints(:,start_up_hgto:end_up_hgto,:,:)
      max_low_hgto = max_low_hgto-1                       !updates the maximum order of HGTOs
      num_low_hgto = num_up_hgto-(max_low_hgto+2)         !number of lower order HGTOs
      ! temporary p-shell CGTO integrals
      allocate(tmp_ints(num_cgto_bra,3*num_geo_bra*num_low_hgto,num_opt), stat=ierr)
      if (ierr/=0)                                                     &
        call error_stop("hgto_to_cgto", "failed to allocate tmp_ints", &
                        num_cgto_bra*3*num_geo_bra*num_low_hgto*num_opt)
      ! 1/(2*\var(exponent_ket)) for recurrence relations
      half_recip_expnt = 1.0_REALK/twice_expnt
      ! gets the temporary p-shell CGTO integrals
      call hgto_to_cgto_p(max_low_hgto, half_recip_expnt, num_cgto_bra, &
                          num_up_hgto, num_geo_bra, num_opt,            &
                          hgto_pints(:,start_up_hgto:end_up_hgto,:,:),  &
                          3, num_low_hgto, tmp_ints)
      ! updates the maximum order of HGTOs
      max_low_hgto = max_low_hgto-1
      ! updates the numbers of lower and upper order HGTOs
      num_up_hgto = num_low_hgto
      num_low_hgto = num_low_hgto-(max_low_hgto+2)
      ! updates the prefactor for order \var(max_low_hgto)-1 HGTOs
      prefact_hgto = half_recip_expnt*half_recip_expnt*prefact_hgto
      ! sets the start and end addresses of HGTOs for order \var(max_low_hgto)-1
      end_up_hgto = start_up_hgto-num_up_hgto
      start_up_hgto = end_up_hgto-num_low_hgto
      end_up_hgto = end_up_hgto-1
#if defined(DEBUG)
      write(STDOUT,100) "HGTO-d-CGTO/order/start/end:", &
                        max_low_hgto, start_up_hgto, end_up_hgto
#endif
      ! multiplies the prefactor for order \var(max_low_hgto)-1 HGTOs
      hgto_pints(:,start_up_hgto:end_up_hgto,:,:) &
        = prefact_hgto*hgto_pints(:,start_up_hgto:end_up_hgto,:,:)
      ! recovers d-shell CGTOs
      call hgto_to_cgto_d(max_low_hgto, half_recip_expnt, num_cgto_bra, &
                          num_low_hgto, num_geo_bra, num_opt,           &
                          hgto_pints(:,start_up_hgto:end_up_hgto,:,:),  &
                          3, num_up_hgto, tmp_ints, 0, 6, cgto_pints)
      deallocate(tmp_ints)
    ! other shell CGTOs returned, at least f-shell
    case default
      ! (1) \var(max_low_hgto)-1 order HGTOs
      prefact_hgto = twice_expnt**max_low_hgto           !prefactor for the maximum order HGTOs
      end_up_hgto = dim_hgto_ket                         !end address of upper order HGTOs
      num_up_hgto = (max_low_hgto+1)*(max_low_hgto+2)/2  !number of upper order HGTOs
      start_up_hgto = end_up_hgto-num_up_hgto+1          !start address of upper order HGTOs
#if defined(DEBUG)
      write(STDOUT,100) "HGTO-s-CGTO/order/start/end:", &
                        max_low_hgto, start_up_hgto, end_up_hgto
#endif
      ! multiplies the prefactor for order \var(max_low_hgto) HGTOs
      hgto_pints(:,start_up_hgto:end_up_hgto,:,:) &
        = prefact_hgto*hgto_pints(:,start_up_hgto:end_up_hgto,:,:)
      max_low_hgto = max_low_hgto-1                      !updates the maximum order of HGTOs
      num_low_hgto = num_up_hgto-(max_low_hgto+2)        !number of lower order HGTOs
      ! allocates memory for temporary integrals
      call dim_hgto_to_cgto(angular_ket, order_geo_ket, dim_tmp)
      size_geo_opt = num_geo_bra*num_opt
      allocate(tmp_ints(num_cgto_bra,dim_tmp*size_geo_opt,2), stat=ierr)
      if (ierr/=0)                                                     &
        call error_stop("hgto_to_cgto", "failed to allocate tmp_ints", &
                        num_cgto_bra*dim_tmp*size_geo_opt*2)
      ! sets the dimension of CGTOs and size of temporary odd order shell integrals, starting from p-shell
      dim_cgto(1) = 3
      size_low_cgto = dim_cgto(1)*num_low_hgto*size_geo_opt
      ! 1/(2*\var(exponent_ket)) for recurrence relations
      half_recip_expnt = 1.0_REALK/twice_expnt
      ! gets the temporary p-shell CGTO integrals
      call hgto_to_cgto_p(max_low_hgto, half_recip_expnt, num_cgto_bra, &
                          num_up_hgto, num_geo_bra, num_opt,            &
                          hgto_pints(:,start_up_hgto:end_up_hgto,:,:),  &
                          dim_cgto(1), num_low_hgto, tmp_ints(:,1:size_low_cgto,1))
      ! (2) \var(max_low_hgto)-2 order HGTOs
      max_low_hgto = max_low_hgto-1
      num_up_hgto = num_low_hgto
      num_low_hgto = num_low_hgto-(max_low_hgto+2)
      ! updates the prefactor for order \var(max_low_hgto)-2 HGTOs
      square_half_rexpnt = half_recip_expnt*half_recip_expnt
      prefact_hgto = square_half_rexpnt*prefact_hgto
      ! sets the start and end addresses of HGTOs for order \var(max_low_hgto)-2
      end_up_hgto = start_up_hgto-num_up_hgto
      start_up_hgto = end_up_hgto-num_low_hgto
      end_up_hgto = end_up_hgto-1
#if defined(DEBUG)
      write(STDOUT,100) "HGTO-d-CGTO/order/start/end:", &
                        max_low_hgto, start_up_hgto, end_up_hgto
#endif
      ! multiplies the prefactor for order \var(max_low_hgto)-2 HGTOs
      hgto_pints(:,start_up_hgto:end_up_hgto,:,:) &
        = prefact_hgto*hgto_pints(:,start_up_hgto:end_up_hgto,:,:)
      ! sets the dimension of CGTOs and size of temporary even order shell integrals, starting from d-shell
      dim_cgto(2) = 6
      size_up_cgto = dim_cgto(2)*num_low_hgto*size_geo_opt
      ! gets the temporary d-shell CGTO integrals
      call hgto_to_cgto_d(max_low_hgto, half_recip_expnt, num_cgto_bra,  &
                          num_low_hgto, num_geo_bra, num_opt,            &
                          hgto_pints(:,start_up_hgto:end_up_hgto,:,:),   &
                          dim_cgto(1), num_up_hgto,                      &
                          tmp_ints(:,1:size_low_cgto,1), 0, dim_cgto(2), &
                          tmp_ints(:,1:size_up_cgto,2))
      ! initializes the maximum order of lower order CGTOs
      max_low_cgto = 1
      ! initializes the pointers of CGTOs
      low_cgto_int = 1
      up_cgto_int = 2
      ! first recurrence relation starts from p-shell
      recur_p_cgto = .true.
      ! loops over the orders of HGTOs on ket center
      do order_hgto = max_low_hgto-1, order_geo_ket, -1
        ! updates the numbers of lower and upper order HGTOs
        num_up_hgto = num_low_hgto
        num_low_hgto = num_low_hgto-(order_hgto+2)  !=(order_hgto+1)*(order_hgto+2)/2
        ! updates the maximum order of CGTOs
        max_low_cgto = max_low_cgto+1
        ! switches the pointers
        low_cgto_int = 3-low_cgto_int
        up_cgto_int = 3-up_cgto_int
        ! updates the dimension of CGTOs
        dim_cgto(up_cgto_int) = dim_cgto(up_cgto_int)+(max_low_cgto+2)*(max_low_cgto+3)/2
        ! updates the sizes of temporary integrals
        size_low_cgto = size_up_cgto
        size_up_cgto = dim_cgto(up_cgto_int)*num_low_hgto*size_geo_opt
        ! gets the temporary p-shell CGTO integrals
        if (recur_p_cgto) then
          call hgto_to_cgto_p(order_hgto, half_recip_expnt, num_cgto_bra,  &
                              num_up_hgto, num_geo_bra, num_opt,           &
                              hgto_pints(:,start_up_hgto:end_up_hgto,:,:), &
                              dim_cgto(up_cgto_int), num_low_hgto,         &
                              tmp_ints(:,1:size_up_cgto,up_cgto_int))
          ! recurrence relation in next cycle starts from d-shell
          recur_p_cgto = .false.
          ! sets the offset and orders of lower order CGTOs for following recurrence relations
          offset_cgto(1) = 0     !offset for d-shell
          offset_cgto(2) = 3     !offset for f-shell
          low_order_cgto(1) = 2  !d-shell
        else
          ! updates the prefactor for order \var(order_hgto) HGTOs
          prefact_hgto = square_half_rexpnt*prefact_hgto
          ! sets the start and end addresses of HGTOs for order \var(order_hgto)
          end_up_hgto = start_up_hgto-num_up_hgto
          start_up_hgto = end_up_hgto-num_low_hgto
          end_up_hgto = end_up_hgto-1
#if defined(DEBUG)
          write(STDOUT,100) "HGTO-CGTO/loop/start/end:", &
                            order_hgto, start_up_hgto, end_up_hgto
#endif
          ! multiplies the prefactor for order \var(order_hgto) HGTOs
          hgto_pints(:,start_up_hgto:end_up_hgto,:,:) &
            = prefact_hgto*hgto_pints(:,start_up_hgto:end_up_hgto,:,:)
          ! gets the temporary d-shell CGTO integrals
          call hgto_to_cgto_d(order_hgto, half_recip_expnt, num_cgto_bra,  &
                              num_low_hgto, num_geo_bra, num_opt,          &
                              hgto_pints(:,start_up_hgto:end_up_hgto,:,:), &
                              dim_cgto(low_cgto_int), num_up_hgto,         &
                              tmp_ints(:,1:size_low_cgto,low_cgto_int),    &
                              0, dim_cgto(up_cgto_int),                    &
                              tmp_ints(:,1:size_up_cgto,up_cgto_int))
          ! recurrence relation in next cycle starts from p-shell
          recur_p_cgto = .true.
          ! sets the offset and orders of lower order CGTOs for following recurrence relations
          offset_cgto(1) = 3     !offset for f-shell
          offset_cgto(2) = 6     !offset for g-shell
          low_order_cgto(1) = 3  !f-shell
        end if
        low_order_cgto(2) = max_low_cgto
        ! gets the temporary integrals for other CGTOs (f, h, ...; or g, i, ...)
        ! and \var(order_hgto) order HGTOs
        call sub_hgto_to_cgto(order_hgto, low_order_cgto, half_recip_expnt, &
                              num_cgto_bra, dim_cgto(low_cgto_int),         &
                              num_geo_bra, num_up_hgto, num_opt,            &
                              tmp_ints(:,1:size_low_cgto,low_cgto_int),     &
                              2, offset_cgto, dim_cgto(up_cgto_int),        &
                              num_low_hgto, tmp_ints(:,1:size_up_cgto,up_cgto_int))
      end do
#if defined(DEBUG)
      write(STDOUT,110) start_up_hgto, "==", 1
      write(STDOUT,110) dim_cgto(up_cgto_int)-num_cgto_ket+1, "->", &
                        dim_cgto(up_cgto_int)
      write(STDOUT,110) num_low_hgto, "==", num_geo_ket
#endif
      ! assigns integrals
      call hgto_to_cgto_assign(num_cgto_bra, dim_cgto(up_cgto_int),    &
                               num_geo_bra, num_geo_ket, num_opt,      &
                               tmp_ints(:,1:size_up_cgto,up_cgto_int), &
                               num_cgto_ket, cgto_pints)
      deallocate(tmp_ints)
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "hgto_to_cgto", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("hgto_to_cgto>> ",A,I6,2I8)
110 format("hgto_to_cgto>> ",I8,2X,A,I8)
#endif
  end subroutine hgto_to_cgto

  !> \brief gets the maximum dimension of temporary integrals used in recurrence relations
  !> \author Bin Gao
  !> \date 2012-03-05
  !> \param angular_ket is the angular number to return
  !> \param order_geo_ket is the order of geometric derivatives to return
  !> \return dim_ints is the maximum dimension of temporary integrals
  subroutine dim_hgto_to_cgto(angular_ket, order_geo_ket, dim_ints)
    use xkind
    implicit none
    integer, intent(in) :: angular_ket
    integer, intent(in) :: order_geo_ket
    integer, intent(out) :: dim_ints
!f2py intent(in) :: angular_ket
!f2py intent(in) :: order_geo_ket
!f2py intent(out) :: dim_ints
    integer max_order_hgto   !maximum order of HGTOs on ket center
    integer num_hgto_ket     !number of HGTOs on ket center
    integer max_order_cgto   !maximum order of CGTOs on ket center
    integer dim_cgto_ket(2)  !dimensions of CGTOs on ket center
    integer which_cgto       !pointer of CGTOs
    integer order_hgto       !incremental recorder over orders of HGTOs
    integer dim_tmp          !temporary result of dimension
#if defined(XTIME)
    real(REALK) curr_time    !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! sets the maximum order and number of HGTOs on ket center
    max_order_hgto = order_geo_ket+angular_ket
    num_hgto_ket = (max_order_hgto+1)*(max_order_hgto+2)/2
    ! initializes the maximum order and dimensions of CGTOs on ket center
    max_order_cgto = 0
    dim_cgto_ket = 0
    ! initializes the pointer of CGTOs
    which_cgto = 2
    ! initializes the return value
    dim_ints = 0
    ! loops over the orders of HGTOs on ket center
    do order_hgto = max_order_hgto-1, order_geo_ket, -1
      ! updates the number of HGTOs
      num_hgto_ket = num_hgto_ket-(order_hgto+2)  !=(order_hgto+1)*(order_hgto+2)/2
      ! updates the maximum order of CGTOs
      max_order_cgto = max_order_cgto+1
      ! resets the pointer of CGTOs
      which_cgto = 3-which_cgto
      ! updates the dimension of CGTOs
      dim_cgto_ket(which_cgto) = dim_cgto_ket(which_cgto) &
                               + (max_order_cgto+1)*(max_order_cgto+2)/2
      ! updates the maximum dimension
      dim_tmp = dim_cgto_ket(which_cgto)*num_hgto_ket
      if (dim_tmp>dim_ints) dim_ints = dim_tmp
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "dim_hgto_to_cgto", STDOUT)
#endif
    return
  end subroutine dim_hgto_to_cgto

  !> \brief transforms p-shell HGTOs to CGTOs on ket center
  !> \author Bin Gao
  !> \date 2011-03-29
  !> \param order_hgto is the order of HGTOs
  !> \param half_recip_expnt is the half reciprocal of the exponent on ket center
  !> \param num_cgto_bra is the number of CGTOs of bra center
  !> \param num_up_hgto is the number of upper order HGTOs
  !> \param num_geo_bra is the number of geometric derivatives on bra center
  !> \param num_opt is the number of operators
  !> \param hgto_pints contains the integrals with s-shell CGTO and upper order HGTOs
  !> \param dim_up_cgto is the dimension of upper order CGTOs
  !> \param num_low_hgto is the number lower order HGTOs
  !> \return up_cgto_pints contains the integrals of p-shell CGTOs and lower order HGTOs
  subroutine hgto_to_cgto_p(order_hgto, half_recip_expnt, num_cgto_bra,    &
                            num_up_hgto, num_geo_bra, num_opt, hgto_pints, &
                            dim_up_cgto, num_low_hgto, up_cgto_pints)
    use xkind
    implicit none
    integer, intent(in) :: order_hgto
    real(REALK), intent(in) :: half_recip_expnt
    integer, intent(in) :: num_cgto_bra
    integer, intent(in) :: num_up_hgto
    integer, intent(in) :: num_geo_bra
    integer, intent(in) :: num_opt
    real(REALK), intent(in) :: hgto_pints(num_cgto_bra,num_up_hgto,num_geo_bra,num_opt)
    integer, intent(in) :: dim_up_cgto
    integer, intent(in) :: num_low_hgto
    real(REALK), intent(inout) :: up_cgto_pints(num_cgto_bra,dim_up_cgto,num_geo_bra, &
                                                num_low_hgto,num_opt)
!f2py intent(in) :: order_hgto
!f2py intent(in) :: half_recip_expnt
!f2py intent(hide) :: num_cgto_bra
!f2py intent(hide) :: num_up_hgto
!f2py intent(hide) :: num_geo_bra
!f2py intent(hide) :: num_opt
!f2py intent(in) :: hgto_pints
!f2py intent(hide) :: dim_up_cgto
!f2py intent(hide) :: num_low_hgto
!f2py intent(inout) :: up_cgto_pints
!f2py depend(num_cgto_bra) :: up_cgto_pints
!f2py depend(num_geo_bra) :: up_cgto_pints
!f2py depend(num_opt) :: up_cgto_pints
    integer iopt           !incremental recorder over operators
    integer iherm, jherm   !incremental recorder over HGTOs
    integer addr_up_hgto   !address of upper order HGTOs
    integer addr_low_hgto  !address of lower order HGTOs
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    do iopt = 1, num_opt
      addr_up_hgto = 0
      addr_low_hgto = 0
      do iherm = order_hgto, 0, -1
        do jherm = 0, iherm
          addr_up_hgto = addr_up_hgto+1
          addr_low_hgto = addr_low_hgto+1
          ! px
          up_cgto_pints(:,1,:,addr_low_hgto,iopt) &
            = half_recip_expnt*hgto_pints(:,addr_up_hgto,:,iopt)
          ! py
          up_cgto_pints(:,2,:,addr_low_hgto,iopt) &
            = half_recip_expnt*hgto_pints(:,addr_up_hgto+1,:,iopt)
          ! pz
          up_cgto_pints(:,3,:,addr_low_hgto,iopt) &
            = half_recip_expnt*hgto_pints(:,addr_up_hgto+iherm+2,:,iopt)
        end do
        addr_up_hgto = addr_up_hgto+1
      end do
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "hgto_to_cgto_p", STDOUT)
#endif
    return
  end subroutine hgto_to_cgto_p

  !> \brief transforms d-shell HGTOs to CGTOs on ket center
  !> \author Bin Gao
  !> \date 2011-03-29
  !> \param order_hgto is the order of HGTOs
  !> \param half_recip_expnt is the half reciprocal of the exponent on ket center
  !> \param num_cgto_bra is the number of CGTOs of bra center
  !> \param num_low_hgto is the number lower order HGTOs
  !> \param num_geo_bra is the number of geometric derivatives on bra center
  !> \param num_opt is the number of operators
  !> \param hgto_pints contains the integrals with s-shell CGTO and lower order HGTOs
  !> \param dim_low_cgto is the dimension of lower order CGTOs
  !> \param num_up_hgto is the number of upper order HGTOs
  !> \param low_cgto_pints contains the integrals with p-shell CGTOs and upper order HGTOs
  !> \param offset_cgto is the offset of upper order CGTOs
  !> \param dim_up_cgto is the dimension of upper order CGTOs
  !> \return up_cgto_pints contains the integrals of d-shell CGTOs and lower order HGTOs
  subroutine hgto_to_cgto_d(order_hgto, half_recip_expnt, num_cgto_bra,     &
                            num_low_hgto, num_geo_bra, num_opt, hgto_pints, &
                            dim_low_cgto, num_up_hgto, low_cgto_pints,      &
                            offset_cgto, dim_up_cgto, up_cgto_pints)
    use xkind
    implicit none
    integer, intent(in) :: order_hgto
    real(REALK), intent(in) :: half_recip_expnt
    integer, intent(in) :: num_cgto_bra
    integer, intent(in) :: num_low_hgto
    integer, intent(in) :: num_geo_bra
    integer, intent(in) :: num_opt
    real(REALK), intent(in) :: hgto_pints(num_cgto_bra,num_low_hgto,num_geo_bra,num_opt)
    integer, intent(in) :: dim_low_cgto
    integer, intent(in) :: num_up_hgto
    real(REALK), intent(in) :: low_cgto_pints(num_cgto_bra,dim_low_cgto,num_geo_bra, &
                                              num_up_hgto,num_opt)
    integer, intent(in) :: offset_cgto
    integer, intent(in) :: dim_up_cgto
    real(REALK), intent(inout) :: up_cgto_pints(num_cgto_bra,dim_up_cgto,num_geo_bra, &
                                                num_low_hgto,num_opt)
!f2py intent(in) :: order_hgto
!f2py intent(in) :: half_recip_expnt
!f2py intent(hide) :: num_cgto_bra
!f2py intent(hide) :: num_low_hgto
!f2py intent(hide) :: num_geo_bra
!f2py intent(hide) :: num_opt
!f2py intent(in) :: hgto_pints
!f2py intent(hide) :: dim_low_cgto
!f2py intent(hide) :: num_up_hgto
!f2py intent(in) :: low_cgto_pints
!f2py depend(num_cgto_bra) :: low_cgto_pints
!f2py depend(num_geo_bra) :: low_cgto_pints
!f2py depend(num_opt) :: low_cgto_pints
!f2py intent(in) :: offset_cgto
!f2py intent(hide) :: dim_up_cgto
!f2py intent(inout) :: up_cgto_pints
!f2py depend(num_cgto_bra) :: up_cgto_pints
!f2py depend(num_geo_bra) :: up_cgto_pints
!f2py depend(num_low_hgto) :: up_cgto_pints
!f2py depend(num_opt) :: up_cgto_pints
    integer iopt            !incremental recorder over operators
    integer iherm, jherm    !incremental recorder over HGTOs
    integer addr_up_hgto    !addresses of upper order HGTOs
    integer addr_up_hgto_y
    integer addr_up_hgto_z
    integer addr_low_hgto   !address of lower order HGTOs
#if defined(XTIME)
    real(REALK) curr_time   !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    do iopt = 1, num_opt
      addr_up_hgto = 0
      addr_low_hgto = 0
      do iherm = order_hgto, 0, -1
        do jherm = 0, iherm
          addr_up_hgto = addr_up_hgto+1
          addr_low_hgto = addr_low_hgto+1
          ! dxx
          up_cgto_pints(:,offset_cgto+1,:,addr_low_hgto,iopt)           &
            = half_recip_expnt*(low_cgto_pints(:,1,:,addr_up_hgto,iopt) &
            + hgto_pints(:,addr_low_hgto,:,iopt))
          addr_up_hgto_y = addr_up_hgto+1
          ! dxy
          up_cgto_pints(:,offset_cgto+2,:,addr_low_hgto,iopt) &
            = half_recip_expnt*low_cgto_pints(:,1,:,addr_up_hgto_y,iopt)
          ! dyy
          up_cgto_pints(:,offset_cgto+3,:,addr_low_hgto,iopt)             &
            = half_recip_expnt*(low_cgto_pints(:,2,:,addr_up_hgto_y,iopt) &
            + hgto_pints(:,addr_low_hgto,:,iopt))
          addr_up_hgto_z = addr_up_hgto+iherm+2
          ! dxz
          up_cgto_pints(:,offset_cgto+4,:,addr_low_hgto,iopt) &
            = half_recip_expnt*low_cgto_pints(:,1,:,addr_up_hgto_z,iopt)
          ! dyz
          up_cgto_pints(:,offset_cgto+5,:,addr_low_hgto,iopt) &
            = half_recip_expnt*low_cgto_pints(:,2,:,addr_up_hgto_z,iopt)
          ! dzz
          up_cgto_pints(:,offset_cgto+6,:,addr_low_hgto,iopt)             &
            = half_recip_expnt*(low_cgto_pints(:,3,:,addr_up_hgto_z,iopt) &
            + hgto_pints(:,addr_low_hgto,:,iopt))
        end do
        addr_up_hgto = addr_up_hgto+1
      end do
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "hgto_to_cgto_d", STDOUT)
#endif
    return
  end subroutine hgto_to_cgto_d

  !> \brief sub-recurrence relations by transforming lower order CGTOs to upper order ones,
  !>        starting from at least f-shell upper order CGTOs
  !> \author Bin Gao
  !> \date 2011-03-29
  !> \param order_hgto is the order of HGTOs
  !> \param low_order_cgto contains the minimum and maximum orders of lower order CGTOs,
  !>        the valid minimum order might be:
  !>        -# 2 means the lower order CGTOs involved in recurrence relations start
  !>           from d-shell, so that the upper order CGTOs start from f-shell
  !>        -# 3 means the lower order CGTOs involved in recurrence relations start
  !>           from f-shell, so that the upper order CGTOs start from g-shell
  !> \param half_recip_expnt is the half reciprocal of the exponent on ket center
  !> \param num_cgto_bra is the number of CGTOs of bra center
  !> \param dim_low_cgto is the dimension of lower order CGTOs
  !> \param num_geo_bra is the number of geometric derivatives on bra center
  !> \param num_up_hgto is the number of upper order HGTOs
  !> \param num_opt is the number of operators
  !> \param low_cgto_pints contains the integrals with lower order CGTOs and upper order HGTOs
  !> \param step_cgto is the step size used when looping different order of CGTOs,
  !>        it should be chosen as 1 or 2, usually for London CGTOs or non-London
  !>        CGTOs respectively
  !> \param offset_cgto contains the offsets of CGTOs in the lower and upper order CGTO
  !>        integrals, the valid options might be:
  !>        -# (/0,3/) all the CGTOs in lower order integrals will be addressed during
  !>           recurrence relations, while the first three CGTOs in upper order integrals
  !>           are p-shell, which have been calculated before, this option should be used
  !>           with \var(low_order_cgto(1))=2 and \var(step_cgto)=2 for getting odd order CGTOs
  !>        -# (/3,6/) the first three CGTOs in lower order integrals are p-shell which
  !>           are not used in the recurrence relations, and the first six CGTOs in upper
  !>           order integrals are d-shell, which have been calculated before, this option
  !>           should be used with \var(low_order_cgto(1))=3 and \var(step_cgto)=2 for
  !>           getting even order CGTOs
  !>        -# (/3,9/) the first three GTOs in lower order integrals are p-shell which
  !>           are not used in the recurrence relations, and the first nine CGTOs in upper
  !>           order integrals are p- and d-shell, which have been calculated before, this
  !>           option should be used with \var(low_order_cgto(1))=2 and \var(step_cgto)=1
  !>           for getting both odd and even orders CGTOs
  !> \param dim_up_cgto is the dimension of upper order CGTOs
  !> \param num_low_hgto is the number lower order HGTOs
  !> \return up_cgto_pints contains the integrals of upper order CGTOs and lower order HGTOs
  subroutine sub_hgto_to_cgto(order_hgto, low_order_cgto, half_recip_expnt,    &
                              num_cgto_bra, dim_low_cgto, num_geo_bra,         &
                              num_up_hgto, num_opt, low_cgto_pints, step_cgto, &
                              offset_cgto, dim_up_cgto, num_low_hgto, up_cgto_pints)
    use xkind
    implicit none
    integer, intent(in) :: order_hgto
    integer, intent(in) :: low_order_cgto(2)
    real(REALK), intent(in) :: half_recip_expnt
    integer, intent(in) :: num_cgto_bra
    integer, intent(in) :: dim_low_cgto
    integer, intent(in) :: num_geo_bra
    integer, intent(in) :: num_up_hgto
    integer, intent(in) :: num_opt
    real(REALK), intent(in) :: low_cgto_pints(num_cgto_bra,dim_low_cgto, &
                                              num_geo_bra,num_up_hgto,num_opt)
    integer, intent(in) :: step_cgto
    integer, intent(in) :: offset_cgto(2)
    integer, intent(in) :: dim_up_cgto
    integer, intent(in) :: num_low_hgto
    real(REALK), intent(inout) :: up_cgto_pints(num_cgto_bra,dim_up_cgto, &
                                                num_geo_bra,num_low_hgto,num_opt)
!f2py intent(in) :: order_hgto
!f2py intent(in) :: low_order_cgto
!f2py intent(in) :: half_recip_expnt
!f2py intent(hide) :: num_cgto_bra
!f2py intent(hide) :: dim_low_cgto
!f2py intent(hide) :: num_geo_bra
!f2py intent(hide) :: num_up_hgto
!f2py intent(hide) :: num_opt
!f2py intent(in) :: low_cgto_pints
!f2py intent(in) :: step_cgto
!f2py intent(in) :: offset_cgto
!f2py intent(hide) :: dim_up_cgto
!f2py intent(hide) :: num_low_hgto
!f2py intent(inout) :: up_cgto_pints
!f2py depend(num_cgto_bra) :: up_cgto_pints
!f2py depend(num_geo_bra) :: up_cgto_pints
!f2py depend(num_opt) :: up_cgto_pints
    integer iopt           !incremental recorder over operators
    integer iherm, jherm   !incremental recorder over HGTOs
    integer addr_up_hgto   !address of upper order HGTOs
    integer addr_low_hgto  !address of lower order HGTOs
    integer addr_up_cgto   !address of upper order CGTOs
    integer addr_cur_cgto  !address of current order CGTOs
    integer addr_low_cgto  !address of lower order CGTOs
    integer icart, jcart   !incremental recorder over CGTOs
    integer order_cgto     !order of CGTOs
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! loops over different operators
    do iopt = 1, num_opt
      addr_up_hgto = 0
      addr_low_hgto = 0
      ! loops over xyz components of HGTOs
      do iherm = order_hgto, 0, -1
        do jherm = 0, iherm
          addr_up_hgto = addr_up_hgto+1
          addr_low_hgto = addr_low_hgto+1
          ! base addresses of CGTOs in the lower and upper order CGTO integrals
          addr_up_cgto = offset_cgto(2)
          addr_cur_cgto = offset_cgto(1)
          addr_low_cgto = 0
          ! loops over the lower order of CGTOs, starting from d- or f-shell
          do order_cgto = low_order_cgto(1), low_order_cgto(2), step_cgto
            addr_up_cgto = addr_up_cgto+1
            addr_cur_cgto = addr_cur_cgto+1
#if defined(DEBUG)
            write(STDOUT,100) "orders:", order_cgto, order_hgto
            write(STDOUT,100) "x-direction"
            write(STDOUT,110) addr_up_cgto, addr_low_hgto
            write(STDOUT,110) addr_cur_cgto, addr_up_hgto
            write(STDOUT,110) addr_low_cgto+1, addr_low_hgto
            write(STDOUT,100) "------------------------------------"
#endif
            ! x-direction
            up_cgto_pints(:,addr_up_cgto,:,addr_low_hgto,iopt)       &
              = half_recip_expnt                                     &
              * (low_cgto_pints(:,addr_cur_cgto,:,addr_up_hgto,iopt) &
              + real(order_cgto,REALK)                               &
              * up_cgto_pints(:,addr_low_cgto+1,:,addr_low_hgto,iopt))
            ! y-direction
            addr_up_cgto = addr_up_cgto+1
#if defined(DEBUG)
            write(STDOUT,100) "y-direction"
            write(STDOUT,110) addr_up_cgto, addr_low_hgto
            write(STDOUT,110) addr_cur_cgto, addr_up_hgto+1
            write(STDOUT,100) "------------------------------------"
#endif
            up_cgto_pints(:,addr_up_cgto,:,addr_low_hgto,iopt) &
              = half_recip_expnt*low_cgto_pints(:,addr_cur_cgto,:,addr_up_hgto+1,iopt)
            do icart = 1, order_cgto
              addr_up_cgto = addr_up_cgto+1
#if defined(DEBUG)
              write(STDOUT,110) addr_up_cgto, addr_low_hgto
              write(STDOUT,110) addr_cur_cgto+icart, addr_up_hgto+1
              write(STDOUT,110) addr_low_cgto+icart, addr_low_hgto
              write(STDOUT,100) "------------------------------------"
#endif
              up_cgto_pints(:,addr_up_cgto,:,addr_low_hgto,iopt)               &
                = half_recip_expnt                                             &
                * (low_cgto_pints(:,addr_cur_cgto+icart,:,addr_up_hgto+1,iopt) &
                + real(icart,REALK)                                            &
                * up_cgto_pints(:,addr_low_cgto+icart,:,addr_low_hgto,iopt))
            end do
            ! z-direction
#if defined(DEBUG)
            write(STDOUT,100) "z-direction"
#endif
            addr_cur_cgto = addr_cur_cgto-1
            do jcart = 0, order_cgto
              addr_up_cgto = addr_up_cgto+1
              addr_cur_cgto = addr_cur_cgto+1
#if defined(DEBUG)
              write(STDOUT,110) addr_up_cgto, addr_low_hgto
              write(STDOUT,110) addr_cur_cgto, addr_up_hgto+iherm+2
              write(STDOUT,100) "------------------------------------"
#endif
              up_cgto_pints(:,addr_up_cgto,:,addr_low_hgto,iopt) &
                = half_recip_expnt                               &
                * low_cgto_pints(:,addr_cur_cgto,:,addr_up_hgto+iherm+2,iopt)
            end do
            do icart = 1, order_cgto
              do jcart = 0, order_cgto-icart
                addr_up_cgto = addr_up_cgto+1
                addr_cur_cgto = addr_cur_cgto+1
                addr_low_cgto = addr_low_cgto+1
#if defined(DEBUG)
                write(STDOUT,110) addr_up_cgto, addr_low_hgto
                write(STDOUT,110) addr_cur_cgto, addr_up_hgto+iherm+2
                write(STDOUT,110) addr_low_cgto, addr_low_hgto
                write(STDOUT,100) "------------------------------------"
#endif
                up_cgto_pints(:,addr_up_cgto,:,addr_low_hgto,iopt)               &
                  = half_recip_expnt                                             &
                  * (low_cgto_pints(:,addr_cur_cgto,:,addr_up_hgto+iherm+2,iopt) &
                  + real(icart,REALK)                                            &
                  * up_cgto_pints(:,addr_low_cgto,:,addr_low_hgto,iopt))
              end do
            end do
          end do
        end do
        addr_up_hgto = addr_up_hgto+1
      end do
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "sub_hgto_to_cgto", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("sub_hgto_to_cgto>> ",A,2I6)
110 format("sub_hgto_to_cgto>> ","CGTO",I8,4X,"HGTO",I8)
#endif
  end subroutine sub_hgto_to_cgto

  !> \brief assigns the primitive CGTO integrals
  !> \author Bin Gao
  !> \date 2012-03-04
  !> \param num_cgto_bra is the number of CGTOs of bra center
  !> \param dim_up_cgto is the dimension of upper order CGTOs in temporary integrals
  !> \param num_geo_bra is the number of geometric derivatives on bra center
  !> \param num_geo_ket is the number of geometric derivatives of ket center
  !> \param num_opt is the number of operators
  !> \param recur_pints contains the temporary integrals from recurrence relations
  !> \param num_cgto_ket is the number of CGTOs of ket center
  !> \return cgto_pints contains the primitive CGTO integrals
  subroutine hgto_to_cgto_assign(num_cgto_bra, dim_up_cgto, num_geo_bra, &
                                 num_geo_ket, num_opt, recur_pints,      &
                                 num_cgto_ket, cgto_pints)
    use xkind
    implicit none
    integer, intent(in) :: num_cgto_bra
    integer, intent(in) :: dim_up_cgto
    integer, intent(in) :: num_geo_bra
    integer, intent(in) :: num_geo_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(in) :: recur_pints(num_cgto_bra,dim_up_cgto, &
                                           num_geo_bra,num_geo_ket,num_opt)
    integer, intent(in) :: num_cgto_ket
    real(REALK), intent(inout) :: cgto_pints(num_cgto_bra,num_cgto_ket, &
                                             num_geo_bra,num_geo_ket,num_opt)
!f2py intent(hide) :: num_cgto_bra
!f2py intent(hide) :: dim_up_cgto
!f2py intent(hide) :: num_geo_bra
!f2py intent(hide) :: num_geo_ket
!f2py intent(hide) :: num_opt
!f2py intent(in) :: recur_pints
!f2py intent(hide) :: num_cgto_ket
!f2py intent(inout) :: cgto_pints
!f2py depend(num_cgto_bra) :: cgto_pints
!f2py depend(num_geo_bra) :: cgto_pints
!f2py depend(num_geo_ket) :: cgto_pints
!f2py depend(num_opt) :: cgto_pints
    integer start_up_cgto  !start address of upper order CGTOs in temporary integrals
    integer iopt           !incremental recorder over operators
    integer iket           !incremental recorder over HGTOs on ket center
    integer igeo           !incremental recorder over geometric derivatives
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    start_up_cgto = dim_up_cgto-num_cgto_ket+1
    do iopt = 1, num_opt
      do iket = 1, num_geo_ket
        do igeo = 1, num_geo_bra
          cgto_pints(:,:,igeo,iket,iopt) &
            = recur_pints(:,start_up_cgto:dim_up_cgto,igeo,iket,iopt)
        end do
      end do
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "hgto_to_cgto_assign", STDOUT)
#endif
    return
  end subroutine hgto_to_cgto_assign
