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
!!  This file transforms several primitive HGTOs to CGTOs on bra or ket center.
!!
!!  2012-03-04, Bin Gao
!!  * rewrites to improve the efficiency
!!
!!  2011-03-29, Bin Gao
!!  * first version

#include "stdout.h"

  !> \brief transforms several primitive HGTOs to CGTOs on ket center
  !> \author Bin Gao
  !> \date 2011-03-29
  !> \note you may refer to \fn(hgto_to_cgto) for non-London atomic orbitals for efficiency
  !> \param angular_ket contains the minimum and maximum angular numbers to return
  !> \param orders_geo_ket contains the minimum and maximum orders of geometric derivatives to return
  !> \param exponent_ket is the exponent of primitive Gaussians of ket center
  !> \param dim_cgto_bra is the dimension of CGTOs of bra center
  !> \param dim_hgto_ket is the dimension of HGTOs of ket center
  !> \param dim_geo_bra is the dimension of geometric derivatives on bra center
  !> \param num_opt is the number of operators
  !> \param hgto_pints contains the primitive HGTO integrals
  !> \param dim_cgto_ket is the dimension of CGTOs of ket center
  !> \param dim_geo_ket is the dimension of geometric derivatives of ket center
  !> \return cgto_pints contains the primitive CGTO integrals
  subroutine hgto_to_lcgto(angular_ket, orders_geo_ket, exponent_ket, &
                           dim_cgto_bra, dim_hgto_ket, dim_geo_bra,   &
                           num_opt, hgto_pints, dim_cgto_ket,         &
                           dim_geo_ket, cgto_pints)
    use xkind
    implicit none
    integer, intent(in) :: angular_ket(2)
    integer, intent(in) :: orders_geo_ket(2)
    real(REALK), intent(in) :: exponent_ket
    integer, intent(in) :: dim_cgto_bra
    integer, intent(in) :: dim_hgto_ket
    integer, intent(in) :: dim_geo_bra
    integer, intent(in) :: num_opt
    real(REALK), intent(inout) :: hgto_pints(dim_cgto_bra,dim_hgto_ket, &
                                             dim_geo_bra,num_opt)
    integer, intent(in) :: dim_cgto_ket
    integer, intent(in) :: dim_geo_ket
    real(REALK), intent(out) :: cgto_pints(dim_cgto_bra,dim_cgto_ket, &
                                           dim_geo_bra,dim_geo_ket,num_opt)
!f2py intent(in) :: angular_ket
!f2py intent(in) :: orders_geo_ket
!f2py intent(in) :: exponent_ket
!f2py intent(hide) :: dim_cgto_bra
!f2py intent(hide) :: dim_hgto_ket
!f2py intent(hide) :: dim_geo_bra
!f2py intent(hide) :: num_opt
!f2py intent(inout) :: hgto_pints
!f2py intent(in) :: dim_cgto_ket
!f2py intent(in) :: dim_geo_ket
!f2py intent(out) :: cgto_pints
!f2py depend(dim_cgto_bra) :: cgto_pints
!f2py depend(dim_cgto_ket) :: cgto_pints
!f2py depend(dim_geo_bra) :: cgto_pints
!f2py depend(dim_geo_ket) :: cgto_pints
!f2py depend(num_opt) :: cgto_pints
    real(REALK) twice_expnt       !twice of the exponent
    real(REALK) min_prefactor     !prefactor for the minimum order HGTOs
    real(REALK) prefact_hgto      !prefactor when transforming HGTOs to CGTOs
    integer start_up_hgto         !start address of upper order HGTOs
    integer end_up_hgto           !end address of upper order HGTOs
    integer num_up_hgto           !number of upper order HGTOs
    integer num_low_hgto          !number of lower order HGTOs
    real(REALK) half_recip_expnt  !half of the reciprocal of exponent
    logical zero_cgto_ket         !if returning zeroth order CGTO
    integer max_low_hgto          !maximum order of lower order HGTOs
    integer max_low_cgto          !maximum order of lower order CGTOs
    integer dim_low_cgto          !dimension of lower order CGTOs of temporary integrals
    integer dim_up_cgto           !dimension of upper order CGTOs of temporary integrals
    integer size_geo_opt          !size of geometric derivatives on bra center and other operators
    integer size_low_cgto         !size of temporary integrals of lower order CGTOs
    integer size_up_cgto          !size of temporary integrals of upper order CGTOs
    integer low_cgto_int          !pointer to temporary integrals of lower order CGTOs
    integer up_cgto_int           !pointer to temporary integrals of upper order CGTOs
    integer dim_tmp               !dimension of temporary integrals
    real(REALK), allocatable :: tmp_ints(:,:,:)
                                  !temporary integrals
    integer offset_geo_ket        !offset of geometric derivatives on ket center in returned integrals
    integer order_hgto            !incremental recorder over orders of HGTOs
    integer iopt                  !incremental recorder over operators
    integer igeo                  !incremental recorder over geometric derivatives
    integer ierr                  !error information
#if defined(XTIME)
    real(REALK) curr_time         !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! twice of the exponent
    twice_expnt = exponent_ket+exponent_ket
    select case(angular_ket(2))
    ! s-shell CGTO returned, only multiplies the prefactor
    case(0)
      ! prefactor for the minimum order HGTOs
      min_prefactor = twice_expnt**orders_geo_ket(1)
      ! loops over different operators
      do iopt = 1, num_opt
        prefact_hgto = min_prefactor
        start_up_hgto = 0
        ! loops over HGTOs on ket center
        do order_hgto = orders_geo_ket(1), orders_geo_ket(2)
          end_up_hgto = start_up_hgto+(order_hgto+1)*(order_hgto+2)/2
          start_up_hgto = start_up_hgto+1
          do igeo = 1, dim_geo_bra
            cgto_pints(:,1,igeo,start_up_hgto:end_up_hgto,iopt) &
              = prefact_hgto*hgto_pints(:,start_up_hgto:end_up_hgto,igeo,iopt)
          end do
          prefact_hgto = twice_expnt*prefact_hgto
          start_up_hgto = end_up_hgto
        end do
      end do
    ! maximum returned CGTOs are p-shell
    case(1)
      prefact_hgto = twice_expnt**(orders_geo_ket(2)+1)            !prefactor for maximum order HGTOs
      end_up_hgto = dim_hgto_ket                                   !end address of upper order HGTOs
      num_up_hgto = (orders_geo_ket(2)+2)*(orders_geo_ket(2)+3)/2  !number of upper order HGTOs
      start_up_hgto = end_up_hgto-num_up_hgto+1                    !start address of upper order HGTOs
#if defined(DEBUG)
      write(STDOUT,100) "HGTO-s-CGTO/p/order/start/end:", &
                        orders_geo_ket(2)+1, start_up_hgto, end_up_hgto
#endif
      ! multiplies the prefactor for order \var(orders_geo_ket(2))+1 HGTOs
      hgto_pints(:,start_up_hgto:end_up_hgto,:,:) &
        = prefact_hgto*hgto_pints(:,start_up_hgto:end_up_hgto,:,:)
      ! number of \var(orders_geo_ket(2)) order HGTOs
      num_low_hgto = num_up_hgto-(orders_geo_ket(2)+2)
      ! sets the size of geometric derivatives on bra center and other operators
      size_geo_opt = 3*dim_geo_bra*num_opt
      ! allocates memory for temporary integrals
      allocate(tmp_ints(dim_cgto_bra,size_geo_opt*num_low_hgto,1), stat=ierr)
      if (ierr/=0)                                                      &
        call error_stop("hgto_to_lcgto", "failed to allocate tmp_ints", &
                        dim_cgto_bra*size_geo_opt*num_low_hgto)
      ! 1/(2*\var(exponent_ket)) for recurrence relations
      half_recip_expnt = 1.0_REALK/twice_expnt
      ! initializes the offset of geometric derivatives on ket center in returned integrals
      offset_geo_ket = dim_geo_ket
      ! if returning zeroth order CGTO on ket center
      zero_cgto_ket = angular_ket(1)==0
      ! loops over returned orders of HGTOs
      do order_hgto = orders_geo_ket(2), orders_geo_ket(1), -1
        ! sets the size of temporary integrals of upper order CGTOs
        size_up_cgto = num_low_hgto*size_geo_opt
        ! gets the p-shell CGTO and order \var(order_hgto) HGTO integrals
        call hgto_to_cgto_p(order_hgto, half_recip_expnt, dim_cgto_bra,  &
                            num_up_hgto, dim_geo_bra, num_opt,           &
                            hgto_pints(:,start_up_hgto:end_up_hgto,:,:), &
                            3, num_low_hgto, tmp_ints(:,1:size_up_cgto,1))
        ! updates the prefactor for order \var(order_hgto)+1 HGTOs
        prefact_hgto = half_recip_expnt*prefact_hgto
        ! sets the start and end addresses of HGTOs for order \var(order_hgto)
        end_up_hgto = start_up_hgto-1
        start_up_hgto = end_up_hgto-num_low_hgto+1
#if defined(DEBUG)
        write(STDOUT,100) "HGTO-CGTO/p/loop/order/start/end:", &
                          order_hgto, start_up_hgto, end_up_hgto
#endif
        ! multiplies the prefactor for HGTOs
        hgto_pints(:,start_up_hgto:end_up_hgto,:,:) &
          = prefact_hgto*hgto_pints(:,start_up_hgto:end_up_hgto,:,:)
        ! updates the offset of geometric derivatives on ket center in returned integrals
        offset_geo_ket = offset_geo_ket-num_low_hgto
        ! assigns the integrals
        call hgto_to_lcgto_assign(start_up_hgto-1, dim_cgto_bra, dim_hgto_ket, &
                                  dim_geo_bra, num_opt, hgto_pints, 3,         &
                                  num_low_hgto, tmp_ints, zero_cgto_ket,       &
                                  offset_geo_ket, dim_cgto_ket, dim_geo_ket,   &
                                  cgto_pints)
        ! updates the numbers of lower and upper order HGTOs
        num_up_hgto = num_low_hgto
        num_low_hgto = num_low_hgto-(order_hgto+1)  !=order_hgto*(order_hgto+1)/2
      end do
      deallocate(tmp_ints)
    ! maximum returned CGTOs are d-shell
    case(2)
      ! (1) \var(max_low_hgto)-1 order HGTOs
      max_low_hgto = orders_geo_ket(2)+2                 !maximum order of HGTOs
      prefact_hgto = twice_expnt**max_low_hgto           !prefactor for maximum order HGTOs
      end_up_hgto = dim_hgto_ket                         !end address of upper order HGTOs
      num_up_hgto = (max_low_hgto+1)*(max_low_hgto+2)/2  !number of upper order HGTOs
      start_up_hgto = end_up_hgto-num_up_hgto+1          !start address of upper order HGTOs
#if defined(DEBUG)
      write(STDOUT,100) "HGTO-s-CGTO/d/order/start/end:", &
                        max_low_hgto, start_up_hgto, end_up_hgto
#endif
      ! multiplies the prefactor for order \var(max_low_hgto) HGTOs
      hgto_pints(:,start_up_hgto:end_up_hgto,:,:) &
        = prefact_hgto*hgto_pints(:,start_up_hgto:end_up_hgto,:,:)
      max_low_hgto = max_low_hgto-1                      !updates the maximum order of HGTOs
      num_low_hgto = num_up_hgto-(max_low_hgto+2)        !number of lower order HGTOs
      ! sets the size of geometric derivatives on bra center and other operators
      size_geo_opt = dim_geo_bra*num_opt
      ! allocates memory for temporary integrals
      allocate(tmp_ints(dim_cgto_bra,9*num_low_hgto*size_geo_opt,2), stat=ierr)
      if (ierr/=0)                                                      &
        call error_stop("hgto_to_lcgto", "failed to allocate tmp_ints", &
                        dim_cgto_bra*9*num_low_hgto*size_geo_opt*2)
      ! sets the size of temporary integrals of lower order CGTOs
      size_low_cgto = 3*num_low_hgto*size_geo_opt
      ! 1/(2*\var(exponent_ket)) for recurrence relations
      half_recip_expnt = 1.0_REALK/twice_expnt
      ! gets the temporary p-shell CGTO integrals
      call hgto_to_cgto_p(max_low_hgto, half_recip_expnt, dim_cgto_bra, &
                          num_up_hgto, dim_geo_bra, num_opt,            &
                          hgto_pints(:,start_up_hgto:end_up_hgto,:,:),  &
                          3, num_low_hgto, tmp_ints(:,1:size_low_cgto,1))
      ! sets the start and end addresses of HGTOs for order \var(max_low_hgto)-1
      end_up_hgto = start_up_hgto-1
      start_up_hgto = end_up_hgto-num_low_hgto+1
      ! updates the prefactor for order \var(max_low_hgto)-1 HGTOs
      prefact_hgto = half_recip_expnt*prefact_hgto
#if defined(DEBUG)
      write(STDOUT,100) "HGTO-p-CGTO/d/order/start/end:", &
                        max_low_hgto, start_up_hgto, end_up_hgto
#endif
      ! multiplies the prefactor for order \var(max_low_hgto)-1 HGTOs
      hgto_pints(:,start_up_hgto:end_up_hgto,:,:) &
        = prefact_hgto*hgto_pints(:,start_up_hgto:end_up_hgto,:,:)
      ! (2) \var(max_low_hgto)-2 order HGTOs
      max_low_hgto = max_low_hgto-1                 !updates the maximum order of HGTOs
      num_up_hgto = num_low_hgto                    !updates the number of upper order HGTOs
      num_low_hgto = num_low_hgto-(max_low_hgto+2)  !updates the number of lower order HGTOs
      ! multiplies the number of CGTOs on ket center, constant 9
      size_geo_opt = 9*size_geo_opt
      ! sets the size of temporary integrals of upper order CGTOs
      size_up_cgto = num_low_hgto*size_geo_opt
      ! gets the temporary p-shell CGTO integrals
      call hgto_to_cgto_p(max_low_hgto, half_recip_expnt, dim_cgto_bra, &
                          num_up_hgto, dim_geo_bra, num_opt,            &
                          hgto_pints(:,start_up_hgto:end_up_hgto,:,:),  &
                          9, num_low_hgto, tmp_ints(:,1:size_up_cgto,2))
      ! updates the prefactor for order \var(max_low_hgto)-2 HGTOs
      prefact_hgto = half_recip_expnt*prefact_hgto
      ! sets the start and end addresses of HGTOs for order \var(max_low_hgto)-2
      end_up_hgto = start_up_hgto-1
      start_up_hgto = end_up_hgto-num_low_hgto+1
#if defined(DEBUG)
      write(STDOUT,100) "HGTO-d-CGTO/d/order/start/end:", &
                        max_low_hgto, start_up_hgto, end_up_hgto
#endif
      ! multiplies the prefactor for order \var(max_low_hgto)-2 HGTOs
      hgto_pints(:,start_up_hgto:end_up_hgto,:,:) &
        = prefact_hgto*hgto_pints(:,start_up_hgto:end_up_hgto,:,:)
      ! recovers d-shell CGTOs
      call hgto_to_cgto_d(max_low_hgto, half_recip_expnt, dim_cgto_bra,  &
                          num_low_hgto, dim_geo_bra, num_opt,            &
                          hgto_pints(:,start_up_hgto:end_up_hgto,:,:),   &
                          3, num_up_hgto, tmp_ints(:,1:size_low_cgto,1), &
                          3, 9, tmp_ints(:,1:size_up_cgto,2))
      ! initializes the offset of geometric derivatives on ket center in returned integrals
      offset_geo_ket = dim_geo_ket-num_low_hgto
      ! if returning zeroth order CGTO on ket center
      zero_cgto_ket = angular_ket(1)==0
      ! assigns the integrals
      call hgto_to_lcgto_assign(start_up_hgto-1, dim_cgto_bra, dim_hgto_ket, &
                                dim_geo_bra, num_opt, hgto_pints, 9,         &
                                num_low_hgto, tmp_ints(:,1:size_up_cgto,2),  &
                                zero_cgto_ket, offset_geo_ket, dim_cgto_ket, &
                                dim_geo_ket, cgto_pints)
      ! initializes the pointers of CGTOs
      low_cgto_int = 1
      up_cgto_int = 2
      ! (3) loops over other returned orders of HGTOs
      do order_hgto = orders_geo_ket(2)-1, orders_geo_ket(1), -1
        ! updates the numbers of lower and upper order HGTOs
        num_up_hgto = num_low_hgto
        num_low_hgto = num_low_hgto-(order_hgto+2)  !=(order_hgto+1)*(order_hgto+2)/2
        ! switches the pointers
        low_cgto_int = 3-low_cgto_int
        up_cgto_int = 3-up_cgto_int
        ! updates the sizes of temporary integrals
        size_low_cgto = size_up_cgto
        size_up_cgto = num_low_hgto*size_geo_opt
        ! gets the temporary p-shell CGTO integrals
        call hgto_to_cgto_p(order_hgto, half_recip_expnt, dim_cgto_bra,  &
                            num_up_hgto, dim_geo_bra, num_opt,           &
                            hgto_pints(:,start_up_hgto:end_up_hgto,:,:), &
                            9, num_low_hgto, tmp_ints(:,1:size_up_cgto,up_cgto_int))
        ! updates the prefactor for order \var(order_hgto) HGTOs
        prefact_hgto = half_recip_expnt*prefact_hgto
        ! sets the start and end addresses of HGTOs for order \var(order_hgto)
        end_up_hgto = start_up_hgto-1
        start_up_hgto = end_up_hgto-num_low_hgto+1
#if defined(DEBUG)
        write(STDOUT,100) "HGTO-CGTO/d/loop/order/start/end:", &
                          order_hgto, start_up_hgto, end_up_hgto
#endif
        ! multiplies the prefactor for order \var(order_hgto) HGTOs
        hgto_pints(:,start_up_hgto:end_up_hgto,:,:) &
          = prefact_hgto*hgto_pints(:,start_up_hgto:end_up_hgto,:,:)
        ! gets the temporary d-shell CGTO integrals
        call hgto_to_cgto_d(order_hgto, half_recip_expnt, dim_cgto_bra,  &
                            num_low_hgto, dim_geo_bra, num_opt,          &
                            hgto_pints(:,start_up_hgto:end_up_hgto,:,:), &
                            9, num_up_hgto,                              &
                            tmp_ints(:,1:size_low_cgto,low_cgto_int),    &
                            3, 9, tmp_ints(:,1:size_up_cgto,up_cgto_int))
        ! updates the offset of geometric derivatives on ket center in returned integrals
        offset_geo_ket = offset_geo_ket-num_low_hgto
        ! assigns the integrals
        call hgto_to_lcgto_assign(start_up_hgto-1, dim_cgto_bra, dim_hgto_ket,       &
                                  dim_geo_bra, num_opt, hgto_pints, 9, num_low_hgto, &
                                  tmp_ints(:,1:size_up_cgto,up_cgto_int),            &
                                  zero_cgto_ket, offset_geo_ket, dim_cgto_ket,       &
                                  dim_geo_ket, cgto_pints)
      end do
      deallocate(tmp_ints)
    ! maximum returned CGTOs are other shells, at least f-shell
    case default
      ! (1) \var(max_low_hgto)-1 order HGTOs
      max_low_hgto = orders_geo_ket(2)+angular_ket(2)    !maximum order of HGTOs
      prefact_hgto = twice_expnt**max_low_hgto           !prefactor for maximum order HGTOs
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
      ! sets the size of geometric derivatives on bra center and other operators
      size_geo_opt = dim_geo_bra*num_opt
      ! allocates memory for temporary integrals
      call dim_nucpot_hket(max_low_hgto, orders_geo_ket(2), dim_tmp)
      allocate(tmp_ints(dim_cgto_bra,dim_tmp*size_geo_opt,2), stat=ierr)
      if (ierr/=0)                                                      &
        call error_stop("hgto_to_lcgto", "failed to allocate tmp_ints", &
                        dim_cgto_bra*dim_tmp*size_geo_opt*2)
      max_low_hgto = max_low_hgto-1                      !updates the maximum order of HGTOs
      num_low_hgto = num_up_hgto-(max_low_hgto+2)        !number of lower order HGTOs
      ! dimension of CGTOs and size of temporary integrals, p-shell
      dim_low_cgto = 3
      size_low_cgto = 3*num_low_hgto*size_geo_opt
      ! 1/(2*\var(exponent_ket)) for recurrence relations
      half_recip_expnt = 1.0_REALK/twice_expnt
      ! gets the temporary p-shell CGTO integrals
      call hgto_to_cgto_p(max_low_hgto, half_recip_expnt, dim_cgto_bra, &
                          num_up_hgto, dim_geo_bra, num_opt,            &
                          hgto_pints(:,start_up_hgto:end_up_hgto,:,:),  &
                          3, num_low_hgto, tmp_ints(:,1:size_low_cgto,1))
      ! sets the start and end addresses of HGTOs for order \var(max_low_hgto)-1
      end_up_hgto = start_up_hgto-1
      start_up_hgto = end_up_hgto-num_low_hgto+1
      ! updates the prefactor for order \var(max_low_hgto)-1 HGTOs
      prefact_hgto = half_recip_expnt*prefact_hgto
#if defined(DEBUG)
      write(STDOUT,100) "HGTO-p-CGTO/order/start/end:", &
                        max_low_hgto, start_up_hgto, end_up_hgto
#endif
      ! multiplies the prefactor for order \var(max_low_hgto)-1 HGTOs
      hgto_pints(:,start_up_hgto:end_up_hgto,:,:) &
        = prefact_hgto*hgto_pints(:,start_up_hgto:end_up_hgto,:,:)
      ! (2) \var(max_low_hgto)-2 order HGTOs
      max_low_hgto = max_low_hgto-1                 !updates the maximum order of HGTOs
      num_up_hgto = num_low_hgto                    !updates the number of upper order HGTOs
      num_low_hgto = num_low_hgto-(max_low_hgto+2)  !updates the number of lower order HGTOs
      ! dimension of CGTOs and size of temporary integrals, p- and d-shell
      dim_up_cgto = 9
      size_up_cgto = 9*num_low_hgto*size_geo_opt
      ! gets the temporary p-shell CGTO integrals
      call hgto_to_cgto_p(max_low_hgto, half_recip_expnt, dim_cgto_bra, &
                          num_up_hgto, dim_geo_bra, num_opt,            &
                          hgto_pints(:,start_up_hgto:end_up_hgto,:,:),  &
                          9, num_low_hgto, tmp_ints(:,1:size_up_cgto,2))
      ! updates the prefactor for order \var(max_low_hgto)-2 HGTOs
      prefact_hgto = half_recip_expnt*prefact_hgto
      ! sets the start and end addresses of HGTOs for order \var(max_low_hgto)-2
      end_up_hgto = start_up_hgto-1
      start_up_hgto = end_up_hgto-num_low_hgto+1
#if defined(DEBUG)
      write(STDOUT,100) "HGTO-d-CGTO/order/start/end:", &
                        max_low_hgto, start_up_hgto, end_up_hgto
#endif
      ! multiplies the prefactor for order \var(max_low_hgto)-2 HGTOs
      hgto_pints(:,start_up_hgto:end_up_hgto,:,:) &
        = prefact_hgto*hgto_pints(:,start_up_hgto:end_up_hgto,:,:)
      ! recovers d-shell CGTOs
      call hgto_to_cgto_d(max_low_hgto, half_recip_expnt, dim_cgto_bra,  &
                          num_low_hgto, dim_geo_bra, num_opt,            &
                          hgto_pints(:,start_up_hgto:end_up_hgto,:,:),   &
                          3, num_up_hgto, tmp_ints(:,1:size_low_cgto,1), &
                          3, 9, tmp_ints(:,1:size_up_cgto,2))
      ! initializes the maximum order of lower order CGTOs
      max_low_cgto = 1
      ! initializes the pointers of CGTOs
      low_cgto_int = 1
      up_cgto_int = 2
      ! (3) loops over the orders of HGTOs not returned and the maximum returned order,
      ! the maximum of lower order CGTOs \var(max_low_cgto) needs to update each iteration
      do order_hgto = max_low_hgto-1, orders_geo_ket(2), -1
        ! updates the numbers of lower and upper order HGTOs
        num_up_hgto = num_low_hgto
        num_low_hgto = num_low_hgto-(order_hgto+2)  !=(order_hgto+1)*(order_hgto+2)/2
        ! updates the maximum order of CGTOs
        max_low_cgto = max_low_cgto+1
        ! updates the dimensions of CGTOs
        dim_low_cgto = dim_up_cgto
        dim_up_cgto = dim_low_cgto+(max_low_cgto+2)*(max_low_cgto+3)/2
        ! updates the sizes of temporary integrals
        size_low_cgto = size_up_cgto
        size_up_cgto = dim_up_cgto*num_low_hgto*size_geo_opt
        ! switches the pointers
        low_cgto_int = 3-low_cgto_int
        up_cgto_int = 3-up_cgto_int
        ! gets the temporary p-shell CGTO integrals
        call hgto_to_cgto_p(order_hgto, half_recip_expnt, dim_cgto_bra,  &
                            num_up_hgto, dim_geo_bra, num_opt,           &
                            hgto_pints(:,start_up_hgto:end_up_hgto,:,:), &
                            dim_up_cgto, num_low_hgto,                   &
                            tmp_ints(:,1:size_up_cgto,up_cgto_int))
        ! updates the prefactor for order \var(order_hgto) HGTOs
        prefact_hgto = half_recip_expnt*prefact_hgto
        ! sets the start and end addresses of HGTOs for order \var(order_hgto)
        end_up_hgto = start_up_hgto-1
        start_up_hgto = end_up_hgto-num_low_hgto+1
#if defined(DEBUG)
        write(STDOUT,100) "HGTO-CGTO/loop/1/order/start/end:", &
                          order_hgto, start_up_hgto, end_up_hgto
#endif
        ! multiplies the prefactor for order \var(order_hgto) HGTOs
        hgto_pints(:,start_up_hgto:end_up_hgto,:,:) &
          = prefact_hgto*hgto_pints(:,start_up_hgto:end_up_hgto,:,:)
        ! gets the temporary d-shell CGTO integrals
        call hgto_to_cgto_d(order_hgto, half_recip_expnt, dim_cgto_bra,  &
                            num_low_hgto, dim_geo_bra, num_opt,          &
                            hgto_pints(:,start_up_hgto:end_up_hgto,:,:), &
                            dim_low_cgto, num_up_hgto,                   &
                            tmp_ints(:,1:size_low_cgto,low_cgto_int), 3, &
                            dim_up_cgto, tmp_ints(:,1:size_up_cgto,up_cgto_int))
        ! gets the temporary integrals for other CGTOs (f, g, ...) and \var(order_hgto) order HGTOs
        call sub_hgto_to_cgto(order_hgto, (/2,max_low_cgto/), half_recip_expnt,     &
                              dim_cgto_bra, dim_low_cgto, dim_geo_bra, num_up_hgto, &
                              num_opt, tmp_ints(:,1:size_low_cgto,low_cgto_int),    &
                              1, (/3,9/), dim_up_cgto, num_low_hgto,                &
                              tmp_ints(:,1:size_up_cgto,up_cgto_int))
      end do
      ! initializes the offset of geometric derivatives on ket center in returned integrals
      offset_geo_ket = dim_geo_ket-num_low_hgto
      ! if returning zeroth order CGTO on ket center
      zero_cgto_ket = angular_ket(1)==0
      ! assigns the integrals
      call hgto_to_lcgto_assign(start_up_hgto-1, dim_cgto_bra, dim_hgto_ket,          &
                                dim_geo_bra, num_opt, hgto_pints, dim_up_cgto,        &
                                num_low_hgto, tmp_ints(:,1:size_up_cgto,up_cgto_int), &
                                zero_cgto_ket, offset_geo_ket, dim_cgto_ket,          &
                                dim_geo_ket, cgto_pints)
      ! updates the dimension of lower order CGTOs
      dim_low_cgto = dim_up_cgto
      ! dimensions of lower and upper order CGTOs do not change from now on
      size_geo_opt = dim_up_cgto*size_geo_opt
      ! (4) loops over other returned orders of HGTOs, the maximum of lower order CGTOs
      ! \var(max_low_cgto) does not need to update
      do order_hgto = orders_geo_ket(2)-1, orders_geo_ket(1), -1
        ! updates the numbers of lower and upper order HGTOs
        num_up_hgto = num_low_hgto
        num_low_hgto = num_low_hgto-(order_hgto+2)  !=(order_hgto+1)*(order_hgto+2)/2
        ! updates the sizes of temporary integrals
        size_low_cgto = size_up_cgto
        size_up_cgto = num_low_hgto*size_geo_opt
        ! switches the pointers
        low_cgto_int = 3-low_cgto_int
        up_cgto_int = 3-up_cgto_int
        ! gets the temporary p-shell CGTO integrals
        call hgto_to_cgto_p(order_hgto, half_recip_expnt, dim_cgto_bra,  &
                            num_up_hgto, dim_geo_bra, num_opt,           &
                            hgto_pints(:,start_up_hgto:end_up_hgto,:,:), &
                            dim_up_cgto, num_low_hgto,                   &
                            tmp_ints(:,1:size_up_cgto,up_cgto_int))
        ! updates the prefactor for order \var(order_hgto) HGTOs
        prefact_hgto = half_recip_expnt*prefact_hgto
        ! sets the start and end addresses of HGTOs for order \var(order_hgto)
        end_up_hgto = start_up_hgto-1
        start_up_hgto = end_up_hgto-num_low_hgto+1
#if defined(DEBUG)
        write(STDOUT,100) "HGTO-CGTO/loop/2/order/start/end:", &
                          order_hgto, start_up_hgto, end_up_hgto
#endif
        ! multiplies the prefactor for order \var(order_hgto) HGTOs
        hgto_pints(:,start_up_hgto:end_up_hgto,:,:) &
          = prefact_hgto*hgto_pints(:,start_up_hgto:end_up_hgto,:,:)
        ! gets the temporary d-shell CGTO integrals
        call hgto_to_cgto_d(order_hgto, half_recip_expnt, dim_cgto_bra,  &
                            num_low_hgto, dim_geo_bra, num_opt,          &
                            hgto_pints(:,start_up_hgto:end_up_hgto,:,:), &
                            dim_low_cgto, num_up_hgto,                   &
                            tmp_ints(:,1:size_low_cgto,low_cgto_int), 3, &
                            dim_up_cgto, tmp_ints(:,1:size_up_cgto,up_cgto_int))
        ! gets the temporary integrals for other CGTOs (f, g, ...) and \var(order_hgto) order HGTOs
        call sub_hgto_to_cgto(order_hgto, (/2,max_low_cgto/), half_recip_expnt,     &
                              dim_cgto_bra, dim_low_cgto, dim_geo_bra, num_up_hgto, &
                              num_opt, tmp_ints(:,1:size_low_cgto,low_cgto_int),    &
                              1, (/3,9/), dim_up_cgto, num_low_hgto,                &
                              tmp_ints(:,1:size_up_cgto,up_cgto_int))
        ! updates the offset of geometric derivatives on ket center in returned integrals
        offset_geo_ket = offset_geo_ket-num_low_hgto
        ! assigns the integrals
        call hgto_to_lcgto_assign(start_up_hgto-1, dim_cgto_bra, dim_hgto_ket,          &
                                  dim_geo_bra, num_opt, hgto_pints, dim_up_cgto,        &
                                  num_low_hgto, tmp_ints(:,1:size_up_cgto,up_cgto_int), &
                                  zero_cgto_ket, offset_geo_ket, dim_cgto_ket,          &
                                  dim_geo_ket, cgto_pints)
      end do
      deallocate(tmp_ints)
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "hgto_to_lcgto", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("hgto_to_lcgto>> ",A,I6,2I8)
#endif
  end subroutine hgto_to_lcgto

  !> \brief assigns the primitive CGTO integrals
  !> \author Bin Gao
  !> \date 2012-03-04
  !> \param offset_hgto_ket is the offset of HGTOs on ket center in the primitive
  !>        HGTO integrals
  !> \param dim_cgto_bra is the dimension of CGTOs of bra center
  !> \param dim_hgto_ket is the dimension of HGTOs of ket center
  !> \param dim_geo_bra is the dimension of geometric derivatives on bra center
  !> \param num_opt is the number of operators
  !> \param hgto_pints contains the primitive HGTO integrals
  !> \param dim_up_cgto is the dimension of upper order CGTOs in temporary integrals
  !> \param num_low_hgto is the number of lower order HGTOs in temporary integrals
  !> \param recur_pints contains the temporary integrals from recurrence relations
  !> \param zero_cgto_ket indicates if zeroth order CGTO returned
  !> \param offset_geo_ket is the offset of geometric derivatives on ket center
  !>        in primitive CGTO integrals
  !> \param dim_cgto_ket is the dimension of CGTOs of ket center
  !> \param dim_geo_ket is the dimension of geometric derivatives of ket center
  !> \return cgto_pints contains the primitive CGTO integrals
  subroutine hgto_to_lcgto_assign(offset_hgto_ket, dim_cgto_bra, dim_hgto_ket,   &
                                  dim_geo_bra, num_opt, hgto_pints, dim_up_cgto, &
                                  num_low_hgto, recur_pints, zero_cgto_ket,      &
                                  offset_geo_ket, dim_cgto_ket, dim_geo_ket,     &
                                  cgto_pints)
    use xkind
    implicit none
    integer, intent(in) :: offset_hgto_ket
    integer, intent(in) :: dim_cgto_bra
    integer, intent(in) :: dim_hgto_ket
    integer, intent(in) :: dim_geo_bra
    integer, intent(in) :: num_opt
    real(REALK), intent(in) :: hgto_pints(dim_cgto_bra,dim_hgto_ket, &
                                          dim_geo_bra,num_opt)
    integer, intent(in) :: dim_up_cgto
    integer, intent(in) :: num_low_hgto
    real(REALK), intent(in) :: recur_pints(dim_cgto_bra,dim_up_cgto,dim_geo_bra, &
                                           num_low_hgto,num_opt)
    logical, intent(in) :: zero_cgto_ket
    integer, intent(in) :: offset_geo_ket
    integer, intent(in) :: dim_cgto_ket
    integer, intent(in) :: dim_geo_ket
    real(REALK), intent(inout) :: cgto_pints(dim_cgto_bra,dim_cgto_ket, &
                                             dim_geo_bra,dim_geo_ket,num_opt)
!f2py intent(in) :: offset_hgto_ket
!f2py intent(hide) :: dim_cgto_bra
!f2py intent(hide) :: dim_hgto_ket
!f2py intent(hide) :: dim_geo_bra
!f2py intent(hide) :: num_opt
!f2py intent(in) :: hgto_pints
!f2py intent(hide) :: dim_up_cgto
!f2py intent(hide) :: num_low_hgto
!f2py intent(in) :: recur_pints
!f2py depend(dim_cgto_bra) :: recur_pints
!f2py depend(dim_geo_bra) :: recur_pints
!f2py depend(num_opt) :: recur_pints
!f2py intent(in) :: zero_cgto_ket
!f2py intent(in) :: offset_geo_ket
!f2py intent(hide) :: dim_cgto_ket
!f2py intent(hide) :: dim_geo_ket
!f2py intent(inout) :: cgto_pints
!f2py depend(dim_cgto_bra) :: cgto_pints
!f2py depend(dim_geo_bra) :: cgto_pints
!f2py depend(num_opt) :: cgto_pints
    integer start_up_cgto  !start address of upper order CGTOs in temporary integrals
    integer addr_hgto_ket  !address of HGTOs on ket center in primitive HGTO integrals
    integer addr_geo_ket   !address of geometric derivatives on ket center in primitive CGTO integrals
    integer iopt           !incremental recorder over operators
    integer iket           !incremental recorder over HGTOs on ket center
    integer igeo           !incremental recorder over geometric derivatives
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! s-shell CGTO returned
    if (zero_cgto_ket) then
      start_up_cgto = dim_up_cgto-dim_cgto_ket+2
      do iopt = 1, num_opt
        do iket = 1, num_low_hgto
          addr_hgto_ket = offset_hgto_ket+iket
          addr_geo_ket = offset_geo_ket+iket
          do igeo = 1, dim_geo_bra
            cgto_pints(:,1,igeo,addr_geo_ket,iopt) &
              = hgto_pints(:,addr_hgto_ket,igeo,iopt)
            cgto_pints(:,2:dim_cgto_ket,igeo,addr_geo_ket,iopt) &
              = recur_pints(:,start_up_cgto:dim_up_cgto,igeo,iket,iopt)
          end do
        end do
      end do
    else
      start_up_cgto = dim_up_cgto-dim_cgto_ket+1
      do iopt = 1, num_opt
        do iket = 1, num_low_hgto
          addr_geo_ket = offset_geo_ket+iket
          do igeo = 1, dim_geo_bra
            cgto_pints(:,:,igeo,addr_geo_ket,iopt) &
              = recur_pints(:,start_up_cgto:dim_up_cgto,igeo,iket,iopt)
          end do
        end do
      end do
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "hgto_to_lcgto_assign", STDOUT)
#endif
    return
  end subroutine hgto_to_lcgto_assign
