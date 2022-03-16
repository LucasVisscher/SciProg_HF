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
!!  This file recovers the order of Cartesian multipole moment.
!!
!!  2012-02-16, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief recovers the order of Cartesian multipole moment
  !> \author Bin Gao
  !> \date 2012-02-16
  !> \param orders_hgto_bra is the range of orders of Hermite Gaussians on bra center
  !> \param coord_bra contains the coordinates of bra center
  !> \param exponent_bra is the exponent of primitive Gaussian of bra center
  !> \param orders_hgto_ket is the range of orders of Hermite Gaussians on ket center
  !> \param coord_ket contains the coordinates of ket center
  !> \param exponent_ket is the exponent of primitive Gaussian of ket center
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param dim_hbra is the dimension of HGTOs on bra center from \fn(carmom_hrr_ket)
  !> \param dim_hket is the dimension of HGTOs on ket center from \fn(carmom_hrr_ket)
  !> \param hket_pints contains the primitive Hermite integrals from \fn(carmom_hrr_ket)
  !> \param dim_hgto_bra is the dimension of HGTOs on bra center
  !> \param dim_hgto_ket is the dimension of HGTOs on ket center
  !> \return hmom_pints contains the primitive Hermite integrals with specified
  !>         orders of HGTOs and Cartesian multipole moment
  subroutine carmom_moment(orders_hgto_bra, coord_bra, exponent_bra, &
                           orders_hgto_ket, coord_ket, exponent_ket, &
                           dipole_origin, order_mom,                 &
                           dim_hbra, dim_hket, hket_pints,           &
                           dim_hgto_bra, dim_hgto_ket, hmom_pints)
    use xkind
    implicit none
    integer, intent(in) :: orders_hgto_bra(2)
    real(REALK), intent(in) :: coord_bra(3)
    real(REALK), intent(in) :: exponent_bra
    integer, intent(in) :: orders_hgto_ket(2)
    real(REALK), intent(in) :: coord_ket(3)
    real(REALK), intent(in) :: exponent_ket
    real(REALK), intent(in) :: dipole_origin(3)
    integer, intent(in) :: order_mom
    integer, intent(in) :: dim_hbra
    integer, intent(in) :: dim_hket
    real(REALK), intent(in) :: hket_pints(dim_hbra,dim_hket)
    integer, intent(in) :: dim_hgto_bra
    integer, intent(in) :: dim_hgto_ket
    real(REALK), intent(out) :: hmom_pints(dim_hgto_bra,dim_hgto_ket, &
                                           (order_mom+1)*(order_mom+2)/2)
!f2py intent(in) :: orders_hgto_bra
!f2py intent(in) :: coord_bra
!f2py intent(in) :: exponent_bra
!f2py intent(in) :: orders_hgto_ket
!f2py intent(in) :: coord_ket
!f2py intent(in) :: exponent_ket
!f2py intent(in) :: dipole_origin
!f2py intent(in) :: order_mom
!f2py intent(hide) :: dim_hbra
!f2py intent(hide) :: dim_hket
!f2py intent(in) :: hket_pints
!f2py intent(in) :: dim_hgto_bra
!f2py intent(in) :: dim_hgto_ket
!f2py intent(out) :: hmom_pints
!f2py depend(dim_hgto_bra) :: hmom_pints
!f2py depend(dim_hgto_ket) :: hmom_pints
!f2py depend(order_mom) :: hmom_pints
    real(REALK) hrp_total_expnt   !half reciprocal of total exponent
    real(REALK) cc_wrt_diporg(3)  !relative coordinates of center-of-charge w.r.t. dipole origin
    integer min_mom_hbra          !minimum order of Cartesian multipole moment from which the minimum order
                                  !of HGTOs on bra center increases
    integer min_mom_hket          !minimum order of Cartesian multipole moment from which the minimum order
                                  !of HGTOs on ket center increases
    integer orders_hbra_up(2)     !orders of HGTOs on bra center in upper order Cartesian multipole moments
    integer orders_hket_up(2)     !orders of HGTOs on ket center in upper order Cartesian multipole moments
    integer dim_hbra_tmp          !dimension of HGTOs on bra center in temporary integrals
    integer dim_hket_tmp          !dimension of HGTOs on ket center in temporary integrals
    integer low_mom_int           !pointer of lower order Cartesian multipole moments
    integer cur_mom_int           !pointer of current order Cartesian multipole moments
    integer up_mom_int            !pointer of upper order Cartesian multipole moments
    real(REALK), allocatable :: tmp_ints(:,:,:,:)
                                  !temporary integrals
    integer offset_hbra_low       !offset of HGTOs on bra center in lower order Cartesian multipole moments
    integer offset_hket_low       !offset of HGTOs on ket center in lower order Cartesian multipole moments
    integer num_up_mom            !number of xyz components of upper order Cartesian multipole moments
    integer num_cur_mom           !number of xyz components of current order Cartesian multipole moments
    integer num_low_mom           !number of xyz components of lower order Cartesian multipole moments
    integer iorder                !incremental recorder over orders of Cartesian multipole moments
    integer ierr                  !error information
#if defined(XTIME)
    real(REALK) curr_time         !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    select case(order_mom)
    ! zeroth order Cartesian multipole moment returned
    case(0)
      hmom_pints(:,:,1) = hket_pints
    ! first order Cartesian multipole moment returned
    case(1)
      ! computes the relative coordinates of center-of-charge w.r.t. dipole origin
      hrp_total_expnt = 1.0_REALK/(exponent_bra+exponent_ket)
      do iorder = 1, 3
        cc_wrt_diporg(iorder) = hrp_total_expnt*(exponent_bra*coord_bra(iorder) &
                              + exponent_ket*coord_ket(iorder))-dipole_origin(iorder)
      end do
      ! computes the half reciprocal of total exponent
      hrp_total_expnt = 0.5_REALK*hrp_total_expnt
      ! sets the minimum orders of Cartesian multipole moment from which the minimum orders
      ! of HGTOs increase
      min_mom_hbra = order_mom-orders_hgto_bra(1)
      min_mom_hket = order_mom-orders_hgto_ket(1)
      ! sets the orders of HGTOs in upper order Cartesian multipole moments
      orders_hbra_up(1) = max(0,1-min_mom_hbra)
      orders_hbra_up(2) = orders_hgto_bra(2)
      orders_hket_up(1) = max(0,1-min_mom_hket)
      orders_hket_up(2) = orders_hgto_ket(2)
      ! gets the first order Cartesian multipole moments
      call carmom_moment_p(orders_hbra_up, orders_hket_up, hrp_total_expnt, &
                           cc_wrt_diporg, dim_hbra, dim_hket, hket_pints,   &
                           dim_hgto_bra, dim_hgto_ket, hmom_pints)
    ! second order Cartesian multipole moment returned
    case(2)
      ! computes the relative coordinates of center-of-charge w.r.t. dipole origin
      hrp_total_expnt = 1.0_REALK/(exponent_bra+exponent_ket)
      do iorder = 1, 3
        cc_wrt_diporg(iorder) = hrp_total_expnt*(exponent_bra*coord_bra(iorder) &
                              + exponent_ket*coord_ket(iorder))-dipole_origin(iorder)
      end do
      ! computes the half reciprocal of total exponent
      hrp_total_expnt = 0.5_REALK*hrp_total_expnt
      ! sets the minimum orders of Cartesian multipole moment from which the minimum orders
      ! of HGTOs increase
      min_mom_hbra = order_mom-orders_hgto_bra(1)
      min_mom_hket = order_mom-orders_hgto_ket(1)
      ! sets the orders of HGTOs in upper order Cartesian multipole moments
      orders_hbra_up(1) = max(0,1-min_mom_hbra)
      orders_hbra_up(2) = orders_hgto_bra(2)
      orders_hket_up(1) = max(0,1-min_mom_hket)
      orders_hket_up(2) = orders_hgto_ket(2)
      ! allocates memory for temporary integrals
      if (orders_hbra_up(1)==0) then
        dim_hbra_tmp = dim_hbra
      else
        dim_hbra_tmp = dim_hbra-orders_hbra_up(1)*(orders_hbra_up(1)+1)/2
      end if
      if (orders_hket_up(1)==0) then
        dim_hket_tmp = dim_hket
      else
        dim_hket_tmp = dim_hket-orders_hket_up(1)*(orders_hket_up(1)+1)/2
      end if
      allocate(tmp_ints(dim_hbra_tmp,dim_hket_tmp,3,1), stat=ierr)
      if (ierr/=0)                                     &
        call error_stop("carmom_moment",               &
                        "failed to allocate tmp_ints", &
                        dim_hbra_tmp*dim_hket_tmp*3)
      ! gets the first order Cartesian multipole moments
      call carmom_moment_p(orders_hbra_up, orders_hket_up, hrp_total_expnt, &
                           cc_wrt_diporg, dim_hbra, dim_hket, hket_pints,   &
                           dim_hbra_tmp, dim_hket_tmp, tmp_ints(:,:,:,1))
      ! sets the minimum order of HGTOs on bra center in upper order Cartesian multipole moments,
      ! and the offset of HGTOs on bra center in lower order Cartesian multipole moments
      if (1<min_mom_hbra) then
        offset_hbra_low = 0
      else if (1==min_mom_hbra) then
        orders_hbra_up(1) = orders_hbra_up(1)+1
        ! skips order \var(orders_hbra_up(1))-1
        offset_hbra_low = orders_hbra_up(1)*(orders_hbra_up(1)+1)/2
      else
        orders_hbra_up(1) = orders_hbra_up(1)+1
        ! skips orders \var(orders_hbra_up(1))-2 and \var(orders_hbra_up(1))-1
        offset_hbra_low = orders_hbra_up(1)*orders_hbra_up(1)
      end if
      ! sets the minimum order of HGTOs on ket center in upper order Cartesian multipole moments,
      ! and the offset of HGTOs on ket center in lower order Cartesian multipole moments
      if (1<min_mom_hket) then
        offset_hket_low = 0
      else if (1==min_mom_hket) then
        orders_hket_up(1) = orders_hket_up(1)+1
        offset_hket_low = orders_hket_up(1)*(orders_hket_up(1)+1)/2
      else
        orders_hket_up(1) = orders_hket_up(1)+1
        offset_hket_low = orders_hket_up(1)*orders_hket_up(1)
      end if
      ! gets the second order Cartesian multipole moments
      call sub_carmom_moment(orders_hbra_up, orders_hket_up, 1, &
                             hrp_total_expnt, cc_wrt_diporg,    &
                             dim_hbra_tmp, dim_hket_tmp, 3,     &
                             tmp_ints(:,:,:,1),                 &
                             offset_hbra_low, offset_hket_low,  &
                             dim_hbra, dim_hket, 1, hket_pints, &
                             dim_hgto_bra, dim_hgto_ket, 6, hmom_pints)
      deallocate(tmp_ints)
    ! higher order (>2) Cartesian multipole moment returned
    case default
      ! computes the relative coordinates of center-of-charge w.r.t. dipole origin
      hrp_total_expnt = 1.0_REALK/(exponent_bra+exponent_ket)
      do iorder = 1, 3
        cc_wrt_diporg(iorder) = hrp_total_expnt*(exponent_bra*coord_bra(iorder) &
                              + exponent_ket*coord_ket(iorder))-dipole_origin(iorder)
      end do
      ! computes the half reciprocal of total exponent
      hrp_total_expnt = 0.5_REALK*hrp_total_expnt
      ! sets the minimum orders of Cartesian multipole moment from which the minimum orders
      ! of HGTOs increase
      min_mom_hbra = order_mom-orders_hgto_bra(1)
      min_mom_hket = order_mom-orders_hgto_ket(1)
      ! sets the orders of HGTOs in upper order Cartesian multipole moments
      orders_hbra_up(1) = max(0,1-min_mom_hbra)
      orders_hbra_up(2) = orders_hgto_bra(2)
      orders_hket_up(1) = max(0,1-min_mom_hket)
      orders_hket_up(2) = orders_hgto_ket(2)
      ! allocates memory for temporary integrals
      if (orders_hbra_up(1)==0) then
        dim_hbra_tmp = dim_hbra
      else
        dim_hbra_tmp = dim_hbra-orders_hbra_up(1)*(orders_hbra_up(1)+1)/2
      end if
      if (orders_hket_up(1)==0) then
        dim_hket_tmp = dim_hket
      else
        dim_hket_tmp = dim_hket-orders_hket_up(1)*(orders_hket_up(1)+1)/2
      end if
      allocate(tmp_ints(dim_hbra_tmp,dim_hket_tmp, &
                        order_mom*(order_mom+1)/2,3), stat=ierr)
      if (ierr/=0)                                     &
        call error_stop("carmom_moment",               &
                        "failed to allocate tmp_ints", &
                        dim_hbra_tmp*dim_hket_tmp*order_mom*(order_mom+1)/2*3)
      ! gets the first order Cartesian multipole moments
      call carmom_moment_p(orders_hbra_up, orders_hket_up, hrp_total_expnt, &
                           cc_wrt_diporg, dim_hbra, dim_hket, hket_pints,   &
                           dim_hbra_tmp, dim_hket_tmp, tmp_ints(:,:,1:3,1))
      ! sets the minimum order of HGTOs on bra center in upper order Cartesian multipole moments,
      ! and the offset of HGTOs on bra center in lower order Cartesian multipole moments
      if (1<min_mom_hbra) then
        offset_hbra_low = 0
      else if (1==min_mom_hbra) then
        orders_hbra_up(1) = orders_hbra_up(1)+1
        ! skips order \var(orders_hbra_up(1))-1
        offset_hbra_low = orders_hbra_up(1)*(orders_hbra_up(1)+1)/2
      else
        orders_hbra_up(1) = orders_hbra_up(1)+1
        ! skips orders \var(orders_hbra_up(1))-2 and \var(orders_hbra_up(1))-1
        offset_hbra_low = orders_hbra_up(1)*orders_hbra_up(1)
      end if
      ! sets the minimum order of HGTOs on ket center in upper order Cartesian multipole moments,
      ! and the offset of HGTOs on ket center in lower order Cartesian multipole moments
      if (1<min_mom_hket) then
        offset_hket_low = 0
      else if (1==min_mom_hket) then
        orders_hket_up(1) = orders_hket_up(1)+1
        offset_hket_low = orders_hket_up(1)*(orders_hket_up(1)+1)/2
      else
        orders_hket_up(1) = orders_hket_up(1)+1
        offset_hket_low = orders_hket_up(1)*orders_hket_up(1)
      end if
      ! gets the second order Cartesian multipole moments
      call sub_carmom_moment(orders_hbra_up, orders_hket_up, 1, &
                             hrp_total_expnt, cc_wrt_diporg,    &
                             dim_hbra_tmp, dim_hket_tmp, 3,     &
                             tmp_ints(:,:,1:3,1),               &
                             offset_hbra_low, offset_hket_low,  &
                             dim_hbra, dim_hket, 1, hket_pints, &
                             dim_hbra_tmp, dim_hket_tmp, 6,     &
                             tmp_ints(:,:,1:6,2))
      ! initializes the number of xyz components of Cartesian multipole moments
      num_cur_mom = 3
      num_up_mom = 6
      ! initializes the pointers of Cartesian multipole moments
      low_mom_int = 3
      cur_mom_int = 1
      up_mom_int = 2
      ! loops over other lower orders of Cartesian multipole moments
      do iorder = 2, order_mom-2
        ! sets the minimum order of HGTOs on bra center in upper order Cartesian multipole moments,
        ! and the offset of HGTOs on bra center in lower order Cartesian multipole moments
        if (iorder<min_mom_hbra) then
          offset_hbra_low = 0
        else if (iorder==min_mom_hbra) then
          orders_hbra_up(1) = orders_hbra_up(1)+1
          offset_hbra_low = orders_hbra_up(1)*(orders_hbra_up(1)+1)/2
        else
          orders_hbra_up(1) = orders_hbra_up(1)+1
          offset_hbra_low = orders_hbra_up(1)*orders_hbra_up(1)
        end if
        ! sets the minimum order of HGTOs on ket center in upper order Cartesian multipole moments,
        ! and the offset of HGTOs on ket center in lower order Cartesian multipole moments
        if (iorder<min_mom_hket) then
          offset_hket_low = 0
        else if (iorder==min_mom_hket) then
          orders_hket_up(1) = orders_hket_up(1)+1
          offset_hket_low = orders_hket_up(1)*(orders_hket_up(1)+1)/2
        else
          orders_hket_up(1) = orders_hket_up(1)+1
          offset_hket_low = orders_hket_up(1)*orders_hket_up(1)
        end if
        ! updates the number of xyz components of Cartesian multipole moments
        num_low_mom = num_cur_mom
        num_cur_mom = num_up_mom
        num_up_mom = num_up_mom+iorder+2  !=(iorder+2)*(iorder+3)/2
        ! updates the pointers of Cartesian multipole moments
        ierr = low_mom_int
        low_mom_int = cur_mom_int
        cur_mom_int = up_mom_int
        up_mom_int = ierr
        ! gets the \var(iorder)+1 order Cartesian multipole moments
        call sub_carmom_moment(orders_hbra_up, orders_hket_up, iorder,  &
                               hrp_total_expnt, cc_wrt_diporg,          &
                               dim_hbra_tmp, dim_hket_tmp, num_cur_mom, &
                               tmp_ints(:,:,1:num_cur_mom,cur_mom_int), &
                               offset_hbra_low, offset_hket_low,        &
                               dim_hbra_tmp, dim_hket_tmp, num_low_mom, &
                               tmp_ints(:,:,1:num_low_mom,low_mom_int), &
                               dim_hbra_tmp, dim_hket_tmp, num_up_mom,  &
                               tmp_ints(:,:,1:num_up_mom,up_mom_int))
      end do
      iorder = order_mom-1
      ! sets the minimum order of HGTOs on bra center in upper order Cartesian multipole moments,
      ! and the offset of HGTOs on bra center in lower order Cartesian multipole moments
      if (iorder<min_mom_hbra) then
        offset_hbra_low = 0
      else if (iorder==min_mom_hbra) then
        orders_hbra_up(1) = orders_hbra_up(1)+1
        offset_hbra_low = orders_hbra_up(1)*(orders_hbra_up(1)+1)/2
      else
        orders_hbra_up(1) = orders_hbra_up(1)+1
        offset_hbra_low = orders_hbra_up(1)*orders_hbra_up(1)
      end if
      ! sets the minimum order of HGTOs on ket center in upper order Cartesian multipole moments,
      ! and the offset of HGTOs on ket center in lower order Cartesian multipole moments
      if (iorder<min_mom_hket) then
        offset_hket_low = 0
      else if (iorder==min_mom_hket) then
        orders_hket_up(1) = orders_hket_up(1)+1
        offset_hket_low = orders_hket_up(1)*(orders_hket_up(1)+1)/2
      else
        orders_hket_up(1) = orders_hket_up(1)+1
        offset_hket_low = orders_hket_up(1)*orders_hket_up(1)
      end if
      ! updates the number of xyz components of Cartesian multipole moments
      num_low_mom = num_cur_mom
      num_cur_mom = num_up_mom
      num_up_mom = num_up_mom+iorder+2  !=(iorder+2)*(iorder+3)/2
      ! updates the pointers of Cartesian multipole moments
      low_mom_int = cur_mom_int
      cur_mom_int = up_mom_int
      ! gets the \var(order_mom) order Cartesian multipole moments
      call sub_carmom_moment(orders_hbra_up, orders_hket_up, iorder,  &
                             hrp_total_expnt, cc_wrt_diporg,          &
                             dim_hbra_tmp, dim_hket_tmp, num_cur_mom, &
                             tmp_ints(:,:,1:num_cur_mom,cur_mom_int), &
                             offset_hbra_low, offset_hket_low,        &
                             dim_hbra_tmp, dim_hket_tmp, num_low_mom, &
                             tmp_ints(:,:,1:num_low_mom,low_mom_int), &
                             dim_hgto_bra, dim_hgto_ket, num_up_mom,  &
                             hmom_pints)
      deallocate(tmp_ints)
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "carmom_moment", STDOUT)
#endif
    return
  end subroutine carmom_moment

  !> \brief recovers the first order Cartesian multipole moments
  !> \author Bin Gao
  !> \date 2012-02-16
  !> \param orders_hbra_up contains the orders of HGTOs on bra center in upper order
  !>        Cartesian multipole moments
  !> \param orders_hket_up contains the orders of HGTOs on ket center in upper order
  !>        Cartesian multipole moments
  !> \param hrp_total_expnt is the half reciprocal of total exponent
  !> \param cc_wrt_diporg contains the relative coordinates of center-of-charge
  !>        w.r.t. dipole origin
  !> \param dim_hbra_cur is the dimension of HGTOs on bra center in current order
  !>        Cartesian multipole moments
  !> \param dim_hket_cur is the dimension of HGTOs on ket center in current order
  !>        Cartesian multipole moments
  !> \param cur_mom_pints contains the integrals of zeroth order Cartesian multipole moment
  !> \param dim_hbra_up is the dimension of HGTOs on bra center in upper order
  !>        Cartesian multipole moments
  !> \param dim_hket_up is the dimension of HGTOs on ket center in upper order
  !>        Cartesian multipole moments
  !> \return up_mom_pints contains the integrals of first order Cartesian multipole moments
  subroutine carmom_moment_p(orders_hbra_up, orders_hket_up, hrp_total_expnt, &
                             cc_wrt_diporg, dim_hbra_cur, dim_hket_cur,       &
                             cur_mom_pints, dim_hbra_up, dim_hket_up, up_mom_pints)
    use xkind
    implicit none
    integer, intent(in) :: orders_hbra_up(2)
    integer, intent(in) :: orders_hket_up(2)
    real(REALK), intent(in) :: hrp_total_expnt
    real(REALK), intent(in) :: cc_wrt_diporg(3)
    integer, intent(in) :: dim_hbra_cur
    integer, intent(in) :: dim_hket_cur
    real(REALK), intent(in) :: cur_mom_pints(dim_hbra_cur,dim_hket_cur)
    integer, intent(in) :: dim_hbra_up
    integer, intent(in) :: dim_hket_up
    real(REALK), intent(out) :: up_mom_pints(dim_hbra_up,dim_hket_up,3)
!f2py intent(in) :: orders_hbra_up
!f2py intent(in) :: orders_hket_up
!f2py intent(in) :: hrp_total_expnt
!f2py intent(in) :: cc_wrt_diporg
!f2py intent(hide) :: dim_hbra_cur
!f2py intent(hide) :: dim_hket_cur
!f2py intent(in) :: cur_mom_pints
!f2py intent(in) :: dim_hbra_up
!f2py intent(in) :: dim_hket_up
!f2py intent(out) :: up_mom_pints
!f2py depend(dim_hbra_up) :: up_mom_pints
!f2py depend(dim_hket_up) :: up_mom_pints
    logical zero_hbra       !if calculating zeroth order HGTO on bra center
    logical zero_hket       !if calculating zeroth order HGTO on ket center
    integer min_order_hbra  !minimum order of HGTOs on bra center
    integer min_order_hket  !minimum order of HGTOs on ket center
    integer base_hbra_cur   !base address of HGTOs on bra center in current order Cartesian multipole moments
    integer base_hket_cur   !base address of HGTOs on ket center in current order Cartesian multipole moments
    integer addr_hbra_up    !address of HGTOs on bra center in upper order Cartesian multipole moments
    integer addr_hbra_cur   !address of HGTOs on bra center in current order Cartesian multipole moments
    integer addr_hket_up    !address of HGTOs on ket center in upper order Cartesian multipole moments
    integer addr_hket_cur   !address of HGTOs on ket center in current order Cartesian multipole moments
    integer addr_low_hbra   !address of lower order HGTOs on bra center
    integer addr_low_hket   !address of lower order HGTOs on ket center
    integer order_hbra      !incremental recorder over the orders of HGTOs on bra center
    integer ibra, jbra      !incremental recorders over xyz components of HGTOs on bra center
    integer order_hket      !incremental recorder over the orders of HGTOs on ket center
    integer iket, jket      !incremental recorders over xyz components of HGTOs on ket center
#if defined(XTIME)
    real(REALK) curr_time   !current CPU time   
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! (1) x component of Cartesian multipole moment
    ! (1.1) zeroth order HGTO on ket center if required
    if (orders_hket_up(1)==0) then
      zero_hket = .true.
      ! sets the minimum order of HGTOs on ket center
      min_order_hket = 1
      ! initializes the address of HGTOs on ket center in upper order Cartesian multipole moments
      addr_hket_up = 1
      ! sets the base address of HGTOs on ket center in current order Cartesian multipole moments
      base_hket_cur = 1
      if (orders_hbra_up(1)==0) then
        zero_hbra = .true.
        ! sets the minimum order of HGTOs on bra center
        min_order_hbra = 1
        ! initializes the address of HGTOs on bra center in upper order Cartesian multipole moments
        addr_hbra_up = 1
        ! sets the base address of HGTOs on bra center in current order Cartesian multipole moments
        base_hbra_cur = 1
        ! recovers the zeroth order HGTO on bra center
        up_mom_pints(1,1,1) = cc_wrt_diporg(1)*cur_mom_pints(1,1)
      else
        zero_hbra = .false.
        min_order_hbra = orders_hbra_up(1)
        addr_hbra_up = 0
        base_hbra_cur = min_order_hbra*(min_order_hbra+1)/2
      end if
      ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
      addr_hbra_cur = base_hbra_cur
      ! initializes the address of lower order HGTOs on bra center
      addr_low_hbra = 0
      ! other order (>0) HGTOs on bra center
      do order_hbra = min_order_hbra, orders_hbra_up(2)
        do ibra = order_hbra, 1, -1
          ! components x...xz...z to xy...yz..z on bra center
          do jbra = ibra, 1, -1
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            addr_low_hbra = addr_low_hbra+1
            up_mom_pints(addr_hbra_up,1,1)                      &
              = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,1) &
              + hrp_total_expnt*real(jbra,REALK)*cur_mom_pints(addr_low_hbra,1)
          end do
          ! component y...yz...z on bra center
          addr_hbra_up = addr_hbra_up+1
          addr_hbra_cur = addr_hbra_cur+1
          up_mom_pints(addr_hbra_up,1,1) &
            = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,1)
        end do
        ! component z...z on bra center
        addr_hbra_up = addr_hbra_up+1
        addr_hbra_cur = addr_hbra_cur+1
        up_mom_pints(addr_hbra_up,1,1) &
          = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,1)
      end do
    else
      zero_hket = .false.
      min_order_hket = orders_hket_up(1)
      addr_hket_up = 0
      base_hket_cur = min_order_hket*(min_order_hket+1)/2
      if (orders_hbra_up(1)==0) then
        zero_hbra = .true.
        min_order_hbra = 1
        base_hbra_cur = 1
      else
        zero_hbra = .false.
        min_order_hbra = orders_hbra_up(1)
        base_hbra_cur = min_order_hbra*(min_order_hbra+1)/2
      end if
    end if
    ! initializes the address of HGTOs on ket center in current order Cartesian multipole moments
    addr_hket_cur = base_hket_cur
    ! initializes the address of lower order HGTOs on ket center
    addr_low_hket = 0
    ! (1.2) other order (>0) HGTOs on ket center
    do order_hket = min_order_hket, orders_hket_up(2)
      do iket = order_hket, 1, -1
        ! components x...xz...z to xy...yz..z on ket center
        do jket = iket, 1, -1
          addr_hket_up = addr_hket_up+1
          addr_hket_cur = addr_hket_cur+1
          addr_low_hket = addr_low_hket+1
          ! zeroth order HGTo on bra center
          if (zero_hbra) then
            ! initializes the address of HGTOs on bra center in upper order Cartesian multipole moments
            addr_hbra_up = 1
            ! recovers the zeroth order HGTO on bra center
            up_mom_pints(1,addr_hket_up,1)                      &
              = cc_wrt_diporg(1)*cur_mom_pints(1,addr_hket_cur) &
              + hrp_total_expnt*real(jket,REALK)*cur_mom_pints(1,addr_low_hket)
          else
            addr_hbra_up = 0
          end if
          ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
          addr_hbra_cur = base_hbra_cur
          ! initializes the address of lower order HGTOs on bra center
          addr_low_hbra = 0
          ! other order (>0) HGTOs on bra center
          do order_hbra = min_order_hbra, orders_hbra_up(2)
            do ibra = order_hbra, 1, -1
              ! components x...xz...z to xy...yz..z on bra center
              do jbra = ibra, 1, -1
                addr_hbra_up = addr_hbra_up+1
                addr_hbra_cur = addr_hbra_cur+1
                addr_low_hbra = addr_low_hbra+1
                up_mom_pints(addr_hbra_up,addr_hket_up,1)                       &
                  = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,addr_hket_cur) &
                  + hrp_total_expnt*(real(jbra,REALK)                           &
                  * cur_mom_pints(addr_low_hbra,addr_hket_cur)                  &
                  + real(jket,REALK)*cur_mom_pints(addr_hbra_cur,addr_low_hket))
              end do
              ! component y...yz...z on bra center
              addr_hbra_up = addr_hbra_up+1
              addr_hbra_cur = addr_hbra_cur+1
              up_mom_pints(addr_hbra_up,addr_hket_up,1)                       &
                = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,addr_hket_cur) &
                + hrp_total_expnt*real(jket,REALK)                            &
                * cur_mom_pints(addr_hbra_cur,addr_low_hket)
            end do
            ! component z...z on bra center
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            up_mom_pints(addr_hbra_up,addr_hket_up,1)                       &
              = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,addr_hket_cur) &
              + hrp_total_expnt*real(jket,REALK)                            &
              * cur_mom_pints(addr_hbra_cur,addr_low_hket)
          end do
        end do
        ! component y...yz...z on ket center
        addr_hket_up = addr_hket_up+1
        addr_hket_cur = addr_hket_cur+1
        ! zeroth order HGTo on bra center
        if (zero_hbra) then
          ! initializes the address of HGTOs on bra center in upper order Cartesian multipole moments
          addr_hbra_up = 1
          ! recovers the zeroth order HGTO on bra center
          up_mom_pints(1,addr_hket_up,1) &
            = cc_wrt_diporg(1)*cur_mom_pints(1,addr_hket_cur)
        else
          addr_hbra_up = 0
        end if
        ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
        addr_hbra_cur = base_hbra_cur
        ! initializes the address of lower order HGTOs on bra center
        addr_low_hbra = 0
        ! other order (>0) HGTOs on bra center
        do order_hbra = min_order_hbra, orders_hbra_up(2)
          do ibra = order_hbra, 1, -1
            ! components x...xz...z to xy...yz..z on bra center
            do jbra = ibra, 1, -1
              addr_hbra_up = addr_hbra_up+1
              addr_hbra_cur = addr_hbra_cur+1
              addr_low_hbra = addr_low_hbra+1
              up_mom_pints(addr_hbra_up,addr_hket_up,1)                       &
                = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,addr_hket_cur) &
                + hrp_total_expnt*real(jbra,REALK)                            &
                * cur_mom_pints(addr_low_hbra,addr_hket_cur)
            end do
            ! component y...yz...z on bra center
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            up_mom_pints(addr_hbra_up,addr_hket_up,1) &
              = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,addr_hket_cur)
          end do
          ! component z...z on bra center
          addr_hbra_up = addr_hbra_up+1
          addr_hbra_cur = addr_hbra_cur+1
          up_mom_pints(addr_hbra_up,addr_hket_up,1) &
            = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,addr_hket_cur)
        end do
      end do
      ! component z...z on ket center
      addr_hket_up = addr_hket_up+1
      addr_hket_cur = addr_hket_cur+1
      ! zeroth order HGTo on bra center
      if (zero_hbra) then
        ! initializes the address of HGTOs on bra center in upper order Cartesian multipole moments
        addr_hbra_up = 1
        ! recovers the zeroth order HGTO on bra center
        up_mom_pints(1,addr_hket_up,1) &
          = cc_wrt_diporg(1)*cur_mom_pints(1,addr_hket_cur)
      else
        addr_hbra_up = 0
      end if
      ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
      addr_hbra_cur = base_hbra_cur
      ! initializes the address of lower order HGTOs on bra center
      addr_low_hbra = 0
      ! other order (>0) HGTOs on bra center
      do order_hbra = min_order_hbra, orders_hbra_up(2)
        do ibra = order_hbra, 1, -1
          ! components x...xz...z to xy...yz..z on bra center
          do jbra = ibra, 1, -1
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            addr_low_hbra = addr_low_hbra+1
            up_mom_pints(addr_hbra_up,addr_hket_up,1)                       &
              = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,addr_hket_cur) &
              + hrp_total_expnt*real(jbra,REALK)                            &
              * cur_mom_pints(addr_low_hbra,addr_hket_cur)
          end do
          ! component y...yz...z on bra center
          addr_hbra_up = addr_hbra_up+1
          addr_hbra_cur = addr_hbra_cur+1
          up_mom_pints(addr_hbra_up,addr_hket_up,1) &
            = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,addr_hket_cur)
        end do
        ! component z...z on bra center
        addr_hbra_up = addr_hbra_up+1
        addr_hbra_cur = addr_hbra_cur+1
        up_mom_pints(addr_hbra_up,addr_hket_up,1) &
          = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,addr_hket_cur)
      end do
    end do
    ! (2) y component of Cartesian multipole moment
    ! (2.1) zeroth order HGTO on ket center if required
    if (zero_hket) then
      addr_hket_up = 1
      ! recovers the zeroth order HGTO on bra center
      if (zero_hbra) then
        addr_hbra_up = 1
        up_mom_pints(1,1,2) = cc_wrt_diporg(2)*cur_mom_pints(1,1)
      else
        addr_hbra_up = 0
      end if
      ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
      addr_hbra_cur = base_hbra_cur
      ! initializes the address of lower order HGTOs on bra center
      addr_low_hbra = 0
      ! other order (>0) HGTOs on bra center
      do order_hbra = min_order_hbra, orders_hbra_up(2)
        do ibra = order_hbra, 1, -1
          ! component x...xz...z on bra center
          addr_hbra_up = addr_hbra_up+1
          addr_hbra_cur = addr_hbra_cur+1
          up_mom_pints(addr_hbra_up,1,2) &
            = cc_wrt_diporg(2)*cur_mom_pints(addr_hbra_cur,1)
          ! components x...xyz...z to y...yz..z on bra center
          do jbra = 1, ibra
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            addr_low_hbra = addr_low_hbra+1
            up_mom_pints(addr_hbra_up,1,2)                      &
              = cc_wrt_diporg(2)*cur_mom_pints(addr_hbra_cur,1) &
              + hrp_total_expnt*real(jbra,REALK)*cur_mom_pints(addr_low_hbra,1)
          end do
        end do
        ! component z...z on bra center
        addr_hbra_up = addr_hbra_up+1
        addr_hbra_cur = addr_hbra_cur+1
        up_mom_pints(addr_hbra_up,1,2) &
          = cc_wrt_diporg(2)*cur_mom_pints(addr_hbra_cur,1)
      end do
    else
      addr_hket_up = 0
    end if
    ! initializes the address of HGTOs on ket center in current order Cartesian multipole moments
    addr_hket_cur = base_hket_cur
    ! initializes the address of lower order HGTOs on ket center
    addr_low_hket = 0
    ! (2.2) other order (>0) HGTOs on ket center
    do order_hket = min_order_hket, orders_hket_up(2)
      do iket = order_hket, 1, -1
        ! component x...xz...z on ket center
        addr_hket_up = addr_hket_up+1
        addr_hket_cur = addr_hket_cur+1
        ! recovers the zeroth order HGTO on bra center
        if (zero_hbra) then
          addr_hbra_up = 1
          up_mom_pints(1,addr_hket_up,2) &
            = cc_wrt_diporg(2)*cur_mom_pints(1,addr_hket_cur)
        else
          addr_hbra_up = 0
        end if
        ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
        addr_hbra_cur = base_hbra_cur
        ! initializes the address of lower order HGTOs on bra center
        addr_low_hbra = 0
        ! other order (>0) HGTOs on bra center
        do order_hbra = min_order_hbra, orders_hbra_up(2)
          do ibra = order_hbra, 1, -1
            ! component x...xz...z on bra center
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            up_mom_pints(addr_hbra_up,addr_hket_up,2) &
              = cc_wrt_diporg(2)*cur_mom_pints(addr_hbra_cur,addr_hket_cur)
            ! components x...xyz...z to y...yz..z on bra center
            do jbra = 1, ibra
              addr_hbra_up = addr_hbra_up+1
              addr_hbra_cur = addr_hbra_cur+1
              addr_low_hbra = addr_low_hbra+1
              up_mom_pints(addr_hbra_up,addr_hket_up,2)                       &
                = cc_wrt_diporg(2)*cur_mom_pints(addr_hbra_cur,addr_hket_cur) &
                + hrp_total_expnt*real(jbra,REALK)                            &
                * cur_mom_pints(addr_low_hbra,addr_hket_cur)
            end do
          end do
          ! component z...z on bra center
          addr_hbra_up = addr_hbra_up+1
          addr_hbra_cur = addr_hbra_cur+1
          up_mom_pints(addr_hbra_up,addr_hket_up,2) &
            = cc_wrt_diporg(2)*cur_mom_pints(addr_hbra_cur,addr_hket_cur)
        end do
        ! components x...xyz...z to y...yz..z on ket center
        do jket = 1, iket
          addr_hket_up = addr_hket_up+1
          addr_hket_cur = addr_hket_cur+1
          addr_low_hket = addr_low_hket+1
          ! recovers the zeroth order HGTO on bra center
          if (zero_hbra) then
            addr_hbra_up = 1
            up_mom_pints(1,addr_hket_up,2)                      &
              = cc_wrt_diporg(2)*cur_mom_pints(1,addr_hket_cur) &
              + hrp_total_expnt*real(jket,REALK)*cur_mom_pints(1,addr_low_hket)
          else
            addr_hbra_up = 0
          end if
          ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
          addr_hbra_cur = base_hbra_cur
          ! initializes the address of lower order HGTOs on bra center
          addr_low_hbra = 0
          ! other order (>0) HGTOs on bra center
          do order_hbra = min_order_hbra, orders_hbra_up(2)
            do ibra = order_hbra, 1, -1
              ! component x...xz...z on bra center
              addr_hbra_up = addr_hbra_up+1
              addr_hbra_cur = addr_hbra_cur+1
              up_mom_pints(addr_hbra_up,addr_hket_up,2)                       &
                = cc_wrt_diporg(2)*cur_mom_pints(addr_hbra_cur,addr_hket_cur) &
                + hrp_total_expnt*real(jket,REALK)                            &
                * cur_mom_pints(addr_hbra_cur,addr_low_hket)
              ! components x...xyz...z to y...yz..z on bra center
              do jbra = 1, ibra
                addr_hbra_up = addr_hbra_up+1
                addr_hbra_cur = addr_hbra_cur+1
                addr_low_hbra = addr_low_hbra+1
                up_mom_pints(addr_hbra_up,addr_hket_up,2)                       &
                  = cc_wrt_diporg(2)*cur_mom_pints(addr_hbra_cur,addr_hket_cur) &
                  + hrp_total_expnt*(real(jbra,REALK)                           &
                  * cur_mom_pints(addr_low_hbra,addr_hket_cur)                  &
                  + real(jket,REALK)*cur_mom_pints(addr_hbra_cur,addr_low_hket))
              end do
            end do
            ! component z...z on bra center
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            up_mom_pints(addr_hbra_up,addr_hket_up,2)                       &
              = cc_wrt_diporg(2)*cur_mom_pints(addr_hbra_cur,addr_hket_cur) &
              + hrp_total_expnt*real(jket,REALK)                            &
              * cur_mom_pints(addr_hbra_cur,addr_low_hket)
          end do
        end do
      end do
      ! component z...z on ket center
      addr_hket_up = addr_hket_up+1
      addr_hket_cur = addr_hket_cur+1
      ! recovers the zeroth order HGTO on bra center
      if (zero_hbra) then
        addr_hbra_up = 1
        up_mom_pints(1,addr_hket_up,2) &
          = cc_wrt_diporg(2)*cur_mom_pints(1,addr_hket_cur)
      else
        addr_hbra_up = 0
      end if
      ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
      addr_hbra_cur = base_hbra_cur
      ! initializes the address of lower order HGTOs on bra center
      addr_low_hbra = 0
      ! other order (>0) HGTOs on bra center
      do order_hbra = min_order_hbra, orders_hbra_up(2)
        do ibra = order_hbra, 1, -1
          ! component x...xz...z on bra center
          addr_hbra_up = addr_hbra_up+1
          addr_hbra_cur = addr_hbra_cur+1
          up_mom_pints(addr_hbra_up,addr_hket_up,2) &
            = cc_wrt_diporg(2)*cur_mom_pints(addr_hbra_cur,addr_hket_cur)
          ! components x...xyz...z to y...yz..z on bra center
          do jbra = 1, ibra
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            addr_low_hbra = addr_low_hbra+1
            up_mom_pints(addr_hbra_up,addr_hket_up,2)                       &
              = cc_wrt_diporg(2)*cur_mom_pints(addr_hbra_cur,addr_hket_cur) &
              + hrp_total_expnt*real(jbra,REALK)                            &
              * cur_mom_pints(addr_low_hbra,addr_hket_cur)
          end do
        end do
        ! component z...z on bra center
        addr_hbra_up = addr_hbra_up+1
        addr_hbra_cur = addr_hbra_cur+1
        up_mom_pints(addr_hbra_up,addr_hket_up,2) &
          = cc_wrt_diporg(2)*cur_mom_pints(addr_hbra_cur,addr_hket_cur)
      end do
    end do
    ! (3) z component of Cartesian multipole moment
    ! (3.1) zeroth order HGTO on ket center if required
    if (zero_hket) then
      addr_hket_up = 1
      ! recovers the zeroth order HGTO on bra center
      if (zero_hbra) then
        addr_hbra_up = 1
        up_mom_pints(1,1,3) = cc_wrt_diporg(3)*cur_mom_pints(1,1)
      else
        addr_hbra_up = 0
      end if
      ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
      addr_hbra_cur = base_hbra_cur
      ! initializes the address of lower order HGTOs on bra center
      addr_low_hbra = 0
      ! other order (>0) HGTOs on bra center
      do order_hbra = min_order_hbra, orders_hbra_up(2)
        ! components x...x to y...y on bra center
        do ibra = 0, order_hbra
          addr_hbra_up = addr_hbra_up+1
          addr_hbra_cur = addr_hbra_cur+1
          up_mom_pints(addr_hbra_up,1,3) &
            = cc_wrt_diporg(3)*cur_mom_pints(addr_hbra_cur,1)
        end do
        ! components x...xz to z...z on bra center
        do ibra = 1, order_hbra
          do jbra = 0, order_hbra-ibra
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            addr_low_hbra = addr_low_hbra+1
            up_mom_pints(addr_hbra_up,1,3)                      &
              = cc_wrt_diporg(3)*cur_mom_pints(addr_hbra_cur,1) &
              + hrp_total_expnt*real(ibra,REALK)*cur_mom_pints(addr_low_hbra,1)
          end do
        end do
      end do
    else
      addr_hket_up = 0
    end if
    ! initializes the address of HGTOs on ket center in current order Cartesian multipole moments
    addr_hket_cur = base_hket_cur
    ! initializes the address of lower order HGTOs on ket center
    addr_low_hket = 0
    ! (3.2) other order (>0) HGTOs on ket center
    do order_hket = min_order_hket, orders_hket_up(2)
      ! components x...x to y...y on ket center
      do iket = 0, order_hket
        addr_hket_up = addr_hket_up+1
        addr_hket_cur = addr_hket_cur+1
        ! recovers the zeroth order HGTO on bra center
        if (zero_hbra) then
          addr_hbra_up = 1
          up_mom_pints(1,addr_hket_up,3) &
            = cc_wrt_diporg(3)*cur_mom_pints(1,addr_hket_cur)
        else
          addr_hbra_up = 0
        end if
        ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
        addr_hbra_cur = base_hbra_cur
        ! initializes the address of lower order HGTOs on bra center
        addr_low_hbra = 0
        ! other order (>0) HGTOs on bra center
        do order_hbra = min_order_hbra, orders_hbra_up(2)
          ! components x...x to y...y on bra center
          do ibra = 0, order_hbra
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            up_mom_pints(addr_hbra_up,addr_hket_up,3) &
              = cc_wrt_diporg(3)*cur_mom_pints(addr_hbra_cur,addr_hket_cur)
          end do
          ! components x...xz to z...z on bra center
          do ibra = 1, order_hbra
            do jbra = 0, order_hbra-ibra
              addr_hbra_up = addr_hbra_up+1
              addr_hbra_cur = addr_hbra_cur+1
              addr_low_hbra = addr_low_hbra+1
              up_mom_pints(addr_hbra_up,addr_hket_up,3)                       &
                = cc_wrt_diporg(3)*cur_mom_pints(addr_hbra_cur,addr_hket_cur) &
                + hrp_total_expnt*real(ibra,REALK)                            &
                * cur_mom_pints(addr_low_hbra,addr_hket_cur)
            end do
          end do
        end do
      end do
      ! components x...xz to z...z on ket center
      do iket = 1, order_hket
        do jket = 0, order_hket-iket
          addr_hket_up = addr_hket_up+1
          addr_hket_cur = addr_hket_cur+1
          addr_low_hket = addr_low_hket+1
          ! recovers the zeroth order HGTO on bra center
          if (zero_hbra) then
            addr_hbra_up = 1
            up_mom_pints(1,addr_hket_up,3)                      &
              = cc_wrt_diporg(3)*cur_mom_pints(1,addr_hket_cur) &
              + hrp_total_expnt*real(iket,REALK)*cur_mom_pints(1,addr_low_hket)
          else
            addr_hbra_up = 0
          end if
          ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
          addr_hbra_cur = base_hbra_cur
          ! initializes the address of lower order HGTOs on bra center
          addr_low_hbra = 0
          ! other order (>0) HGTOs on bra center
          do order_hbra = min_order_hbra, orders_hbra_up(2)
            ! components x...x to y...y on bra center
            do ibra = 0, order_hbra
              addr_hbra_up = addr_hbra_up+1
              addr_hbra_cur = addr_hbra_cur+1
              up_mom_pints(addr_hbra_up,addr_hket_up,3)                       &
                = cc_wrt_diporg(3)*cur_mom_pints(addr_hbra_cur,addr_hket_cur) &
                + hrp_total_expnt*real(iket,REALK)                            &
                * cur_mom_pints(addr_hbra_cur,addr_low_hket)
            end do
            ! components x...xz to z...z on bra center
            do ibra = 1, order_hbra
              do jbra = 0, order_hbra-ibra
                addr_hbra_up = addr_hbra_up+1
                addr_hbra_cur = addr_hbra_cur+1
                addr_low_hbra = addr_low_hbra+1
                up_mom_pints(addr_hbra_up,addr_hket_up,3)                       &
                  = cc_wrt_diporg(3)*cur_mom_pints(addr_hbra_cur,addr_hket_cur) &
                  + hrp_total_expnt*(real(ibra,REALK)                           &
                  * cur_mom_pints(addr_low_hbra,addr_hket_cur)                  &
                  + real(iket,REALK)*cur_mom_pints(addr_hbra_cur,addr_low_hket))
              end do
            end do
          end do
        end do
      end do
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "carmom_moment_p", STDOUT)
#endif
    return
  end subroutine carmom_moment_p

  !> \brief recovers the specific order (>1) Cartesian multipole moments
  !> \author Bin Gao
  !> \date 2012-02-16
  !> \param orders_hbra_up contains the orders of HGTOs on bra center in upper order
  !>        Cartesian multipole moments
  !> \param orders_hket_up contains the orders of HGTOs on ket center in upper order
  !>        Cartesian multipole moments
  !> \param cur_order_mom is the current order of Cartesian multipole moments
  !> \param hrp_total_expnt is the half reciprocal of total exponent
  !> \param cc_wrt_diporg contains the relative coordinates of center-of-charge
  !>        w.r.t. dipole origin
  !> \param dim_hbra_cur is the dimension of HGTOs on bra center in current order
  !>        Cartesian multipole moments
  !> \param dim_hket_cur is the dimension of HGTOs on ket center in current order
  !>        Cartesian multipole moments
  !> \param num_cur_mom is the number of xyz components of current order Cartesian multipole moments
  !> \param cur_mom_pints contains the integrals of current order Cartesian multipole moments
  !> \param offset_hbra_low is the offset of HGTOs on bra center in lower order
  !>        Cartesian multipole moments
  !> \param offset_hket_low is the offset of HGTOs on ket center in lower order
  !>        Cartesian multipole moments
  !> \param dim_hbra_low is the dimension of HGTOs on bra center in lower order
  !>        Cartesian multipole moments
  !> \param dim_hket_low is the dimension of HGTOs on ket center in lower order
  !>        Cartesian multipole moments
  !> \param num_low_mom number of xyz components of lower order Cartesian multipole moments
  !> \param low_mom_pints contains the integrals of lower order Cartesian multipole moments
  !> \param dim_hbra_up is the dimension of HGTOs on bra center in upper order
  !>        Cartesian multipole moments
  !> \param dim_hket_up is the dimension of HGTOs on ket center in upper order
  !>        Cartesian multipole moments
  !> \param num_up_mom is the number of xyz components of upper order Cartesian multipole moments
  !> \return up_mom_pints contains the integrals of upper order Cartesian multipole moments
  subroutine sub_carmom_moment(orders_hbra_up, orders_hket_up, cur_order_mom,  &
                               hrp_total_expnt, cc_wrt_diporg, dim_hbra_cur,   &
                               dim_hket_cur, num_cur_mom, cur_mom_pints,       &
                               offset_hbra_low, offset_hket_low, dim_hbra_low, &
                               dim_hket_low, num_low_mom, low_mom_pints,       &
                               dim_hbra_up, dim_hket_up, num_up_mom, up_mom_pints)
    use xkind
    implicit none
    integer, intent(in) :: orders_hbra_up(2)
    integer, intent(in) :: orders_hket_up(2)
    integer, intent(in) :: cur_order_mom
    real(REALK), intent(in) :: hrp_total_expnt
    real(REALK), intent(in) :: cc_wrt_diporg(3)
    integer, intent(in) :: dim_hbra_cur
    integer, intent(in) :: dim_hket_cur
    integer, intent(in) :: num_cur_mom
    real(REALK), intent(in) :: cur_mom_pints(dim_hbra_cur,dim_hket_cur,num_cur_mom)
    integer, intent(in) :: offset_hbra_low
    integer, intent(in) :: offset_hket_low
    integer, intent(in) :: dim_hbra_low
    integer, intent(in) :: dim_hket_low
    integer, intent(in) :: num_low_mom
    real(REALK), intent(in) :: low_mom_pints(dim_hbra_low,dim_hket_low,num_low_mom)
    integer, intent(in) :: dim_hbra_up
    integer, intent(in) :: dim_hket_up
    integer, intent(in) :: num_up_mom
    real(REALK), intent(out) :: up_mom_pints(dim_hbra_up,dim_hket_up,num_up_mom)
!f2py intent(in) :: orders_hbra_up
!f2py intent(in) :: orders_hket_up
!f2py intent(in) :: cur_order_mom
!f2py intent(in) :: hrp_total_expnt
!f2py intent(in) :: cc_wrt_diporg
!f2py intent(hide) :: dim_hbra_cur
!f2py intent(hide) :: dim_hket_cur
!f2py intent(hide) :: num_cur_mom
!f2py intent(in) :: cur_mom_pints
!f2py intent(in) :: offset_hbra_low
!f2py intent(in) :: offset_hket_low
!f2py intent(hide) :: dim_hbra_low
!f2py intent(hide) :: dim_hket_low
!f2py intent(hide) :: num_low_mom
!f2py intent(in) :: low_mom_pints
!f2py intent(in) :: dim_hbra_up
!f2py intent(in) :: dim_hket_up
!f2py intent(in) :: num_up_mom
!f2py intent(out) :: up_mom_pints
!f2py depend(dim_hbra_up) :: up_mom_pints
!f2py depend(dim_hket_up) :: up_mom_pints
!f2py depend(num_up_mom) :: up_mom_pints
    logical zero_hbra       !if calculating zeroth order HGTO on bra center
    logical zero_hket       !if calculating zeroth order HGTO on ket center
    integer min_order_hbra  !minimum order of HGTOs on bra center
    integer min_order_hket  !minimum order of HGTOs on ket center
    integer base_hbra_cur   !base address of HGTOs on bra center in current order Cartesian multipole moments
    integer base_hket_cur   !base address of HGTOs on ket center in current order Cartesian multipole moments
    integer addr_hbra_up    !address of HGTOs on bra center in upper order Cartesian multipole moments
    integer addr_hbra_cur   !address of HGTOs on bra center in current order Cartesian multipole moments
    integer addr_hbra_low   !address of HGTOs on bra center in lower order Cartesian multipole moments
    integer addr_hket_up    !address of HGTOs on ket center in upper order Cartesian multipole moments
    integer addr_hket_cur   !address of HGTOs on ket center in current order Cartesian multipole moments
    integer addr_hket_low   !address of HGTOs on ket center in lower order Cartesian multipole moments
    integer addr_low_hbra   !address of lower order HGTOs on bra center
    integer addr_low_hket   !address of lower order HGTOs on ket center
    integer order_hbra      !incremental recorder over the orders of HGTOs on bra center
    integer ibra, jbra      !incremental recorders over xyz components of HGTOs on bra center
    integer order_hket      !incremental recorder over the orders of HGTOs on ket center
    integer iket, jket      !incremental recorders over xyz components of HGTOs on ket center
    integer addr_up_mom     !address of upper order Cartesian multipole moment
    integer addr_cur_mom    !address of current order Cartesian multipole moment
    integer addr_low_mom    !address of lower order Cartesian multipole moment
    integer imom, jmom      !incremental recorders over xyz components of Cartesian multipole moments
#if defined(XTIME)
    real(REALK) curr_time   !current CPU time   
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! (1) x...x component of Cartesian multipole moment
    ! (1.1) zeroth order HGTO on ket center if required
    if (orders_hket_up(1)==0) then
      zero_hket = .true.
      ! sets the minimum order of HGTOs on ket center
      min_order_hket = 1
      ! initializes the address of HGTOs on ket center in upper order Cartesian multipole moments
      addr_hket_up = 1
      ! sets the base address of HGTOs on ket center in current order Cartesian multipole moments
      base_hket_cur = 1
      ! sets the address of HGTOs on ket center in lower order Cartesian multipole moments
      addr_hket_low = offset_hket_low+1
      if (orders_hbra_up(1)==0) then
        zero_hbra = .true.
        ! sets the minimum order of HGTOs on bra center
        min_order_hbra = 1
        ! initializes the address of HGTOs on bra center in upper order Cartesian multipole moments
        addr_hbra_up = 1
        ! sets the base address of HGTOs on bra center in current order Cartesian multipole moments
        base_hbra_cur = 1
        ! sets the address of HGTOs on bra center in lower order Cartesian multipole moments
        addr_hbra_low = offset_hbra_low+1
        ! recovers the zeroth order HGTO on bra center
        up_mom_pints(1,1,1)                           &
          = cc_wrt_diporg(1)*cur_mom_pints(1,1,1)     &
          + hrp_total_expnt*real(cur_order_mom,REALK) &
          * low_mom_pints(addr_hbra_low,addr_hket_low,1)
      else
        zero_hbra = .false.
        min_order_hbra = orders_hbra_up(1)
        addr_hbra_up = 0
        base_hbra_cur = min_order_hbra*(min_order_hbra+1)/2
        addr_hbra_low = offset_hbra_low
      end if
      ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
      addr_hbra_cur = base_hbra_cur
      ! initializes the address of lower order HGTOs on bra center
      addr_low_hbra = 0
      ! other order (>0) HGTOs on bra center
      do order_hbra = min_order_hbra, orders_hbra_up(2)
        do ibra = order_hbra, 1, -1
          ! components x...xz...z to xy...yz..z on bra center
          do jbra = ibra, 1, -1
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            addr_hbra_low = addr_hbra_low+1
            addr_low_hbra = addr_low_hbra+1
            up_mom_pints(addr_hbra_up,1,1)                                         &
              = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,1,1)                  &
              + hrp_total_expnt*(real(jbra,REALK)*cur_mom_pints(addr_low_hbra,1,1) &
              + real(cur_order_mom,REALK)*low_mom_pints(addr_hbra_low,addr_hket_low,1))
          end do
          ! component y...yz...z on bra center
          addr_hbra_up = addr_hbra_up+1
          addr_hbra_cur = addr_hbra_cur+1
          addr_hbra_low = addr_hbra_low+1
          up_mom_pints(addr_hbra_up,1,1)                        &
            = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,1,1) &
            + hrp_total_expnt*real(cur_order_mom,REALK)         &
            * low_mom_pints(addr_hbra_low,addr_hket_low,1)
        end do
        ! component z...z on bra center
        addr_hbra_up = addr_hbra_up+1
        addr_hbra_cur = addr_hbra_cur+1
        addr_hbra_low = addr_hbra_low+1
        up_mom_pints(addr_hbra_up,1,1)                        &
          = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,1,1) &
          + hrp_total_expnt*real(cur_order_mom,REALK)         &
          * low_mom_pints(addr_hbra_low,addr_hket_low,1)
      end do
    else
      zero_hket = .false.
      min_order_hket = orders_hket_up(1)
      addr_hket_up = 0
      base_hket_cur = min_order_hket*(min_order_hket+1)/2
      addr_hket_low = offset_hket_low
      if (orders_hbra_up(1)==0) then
        zero_hbra = .true.
        min_order_hbra = 1
        base_hbra_cur = 1
      else
        zero_hbra = .false.
        min_order_hbra = orders_hbra_up(1)
        base_hbra_cur = min_order_hbra*(min_order_hbra+1)/2
      end if
    end if
    ! initializes the address of HGTOs on ket center in current order Cartesian multipole moments
    addr_hket_cur = base_hket_cur
    ! initializes the address of lower order HGTOs on ket center
    addr_low_hket = 0
    ! (1.2) other order (>0) HGTOs on ket center
    do order_hket = min_order_hket, orders_hket_up(2)
      do iket = order_hket, 1, -1
        ! components x...xz...z to xy...yz..z on ket center
        do jket = iket, 1, -1
          addr_hket_up = addr_hket_up+1
          addr_hket_cur = addr_hket_cur+1
          addr_hket_low = addr_hket_low+1
          addr_low_hket = addr_low_hket+1
          ! recovers the zeroth order HGTo on bra center
          if (zero_hbra) then
            addr_hbra_up = 1
            addr_hbra_low = offset_hbra_low+1
            up_mom_pints(1,addr_hket_up,1)                        &
              = cc_wrt_diporg(1)*cur_mom_pints(1,addr_hket_cur,1) &
              + hrp_total_expnt*(real(jket,REALK)                 &
              * cur_mom_pints(1,addr_low_hket,1)                  &
              + real(cur_order_mom,REALK)                         &
              * low_mom_pints(addr_hbra_low,addr_hket_low,1))
          else
            addr_hbra_up = 0
            addr_hbra_low = offset_hbra_low
          end if
          ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
          addr_hbra_cur = base_hbra_cur
          ! initializes the address of lower order HGTOs on bra center
          addr_low_hbra = 0
          ! other order (>0) HGTOs on bra center
          do order_hbra = min_order_hbra, orders_hbra_up(2)
            do ibra = order_hbra, 1, -1
              ! components x...xz...z to xy...yz..z on bra center
              do jbra = ibra, 1, -1
                addr_hbra_up = addr_hbra_up+1
                addr_hbra_cur = addr_hbra_cur+1
                addr_hbra_low = addr_hbra_low+1
                addr_low_hbra = addr_low_hbra+1
                up_mom_pints(addr_hbra_up,addr_hket_up,1)                         &
                  = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,addr_hket_cur,1) &
                  + hrp_total_expnt*(real(jbra,REALK)                             &
                  * cur_mom_pints(addr_low_hbra,addr_hket_cur,1)                  &
                  + real(jket,REALK)*cur_mom_pints(addr_hbra_cur,addr_low_hket,1) &
                  + real(cur_order_mom,REALK)                                     &
                  * low_mom_pints(addr_hbra_low,addr_hket_low,1))
              end do
              ! component y...yz...z on bra center
              addr_hbra_up = addr_hbra_up+1
              addr_hbra_cur = addr_hbra_cur+1
              addr_hbra_low = addr_hbra_low+1
              up_mom_pints(addr_hbra_up,addr_hket_up,1)                         &
                = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,addr_hket_cur,1) &
                + hrp_total_expnt*(real(jket,REALK)                             &
                * cur_mom_pints(addr_hbra_cur,addr_low_hket,1)                  &
                + real(cur_order_mom,REALK)                                     &
                * low_mom_pints(addr_hbra_low,addr_hket_low,1))
            end do
            ! component z...z on bra center
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            addr_hbra_low = addr_hbra_low+1
            up_mom_pints(addr_hbra_up,addr_hket_up,1)                         &
              = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,addr_hket_cur,1) &
              + hrp_total_expnt*(real(jket,REALK)                             &
              * cur_mom_pints(addr_hbra_cur,addr_low_hket,1)                  &
              + real(cur_order_mom,REALK)                                     &
              * low_mom_pints(addr_hbra_low,addr_hket_low,1))
          end do
        end do
        ! component y...yz...z on ket center
        addr_hket_up = addr_hket_up+1
        addr_hket_cur = addr_hket_cur+1
        addr_hket_low = addr_hket_low+1
        ! recovers the zeroth order HGTo on bra center
        if (zero_hbra) then
          addr_hbra_up = 1
          addr_hbra_low = offset_hbra_low+1
          up_mom_pints(1,addr_hket_up,1)                        &
            = cc_wrt_diporg(1)*cur_mom_pints(1,addr_hket_cur,1) &
            + hrp_total_expnt*real(cur_order_mom,REALK)         &
            * low_mom_pints(addr_hbra_low,addr_hket_low,1)
        else
          addr_hbra_up = 0
          addr_hbra_low = offset_hbra_low
        end if
        ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
        addr_hbra_cur = base_hbra_cur
        ! initializes the address of lower order HGTOs on bra center
        addr_low_hbra = 0
        ! other order (>0) HGTOs on bra center
        do order_hbra = min_order_hbra, orders_hbra_up(2)
          do ibra = order_hbra, 1, -1
            ! components x...xz...z to xy...yz..z on bra center
            do jbra = ibra, 1, -1
              addr_hbra_up = addr_hbra_up+1
              addr_hbra_cur = addr_hbra_cur+1
              addr_hbra_low = addr_hbra_low+1
              addr_low_hbra = addr_low_hbra+1
              up_mom_pints(addr_hbra_up,addr_hket_up,1)                         &
                = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,addr_hket_cur,1) &
                + hrp_total_expnt*(real(jbra,REALK)                             &
                * cur_mom_pints(addr_low_hbra,addr_hket_cur,1)                  &
                + real(cur_order_mom,REALK)                                     &
                * low_mom_pints(addr_hbra_low,addr_hket_low,1))
            end do
            ! component y...yz...z on bra center
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            addr_hbra_low = addr_hbra_low+1
            up_mom_pints(addr_hbra_up,addr_hket_up,1)                         &
              = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,addr_hket_cur,1) &
              + hrp_total_expnt*real(cur_order_mom,REALK)                     &
              * low_mom_pints(addr_hbra_low,addr_hket_low,1)
          end do
          ! component z...z on bra center
          addr_hbra_up = addr_hbra_up+1
          addr_hbra_cur = addr_hbra_cur+1
          addr_hbra_low = addr_hbra_low+1
          up_mom_pints(addr_hbra_up,addr_hket_up,1)                         &
            = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,addr_hket_cur,1) &
            + hrp_total_expnt*real(cur_order_mom,REALK)                     &
            * low_mom_pints(addr_hbra_low,addr_hket_low,1)
        end do
      end do
      ! component z...z on ket center
      addr_hket_up = addr_hket_up+1
      addr_hket_cur = addr_hket_cur+1
      addr_hket_low = addr_hket_low+1
      ! recovers the zeroth order HGTo on bra center
      if (zero_hbra) then
        addr_hbra_up = 1
        addr_hbra_low = offset_hbra_low+1
        up_mom_pints(1,addr_hket_up,1)                        &
          = cc_wrt_diporg(1)*cur_mom_pints(1,addr_hket_cur,1) &
          + hrp_total_expnt*real(cur_order_mom,REALK)         &
          * low_mom_pints(addr_hbra_low,addr_hket_low,1)
      else
        addr_hbra_up = 0
        addr_hbra_low = offset_hbra_low
      end if
      ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
      addr_hbra_cur = base_hbra_cur
      ! initializes the address of lower order HGTOs on bra center
      addr_low_hbra = 0
      ! other order (>0) HGTOs on bra center
      do order_hbra = min_order_hbra, orders_hbra_up(2)
        do ibra = order_hbra, 1, -1
          ! components x...xz...z to xy...yz..z on bra center
          do jbra = ibra, 1, -1
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            addr_hbra_low = addr_hbra_low+1
            addr_low_hbra = addr_low_hbra+1
            up_mom_pints(addr_hbra_up,addr_hket_up,1)                         &
              = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,addr_hket_cur,1) &
              + hrp_total_expnt*(real(jbra,REALK)                             &
              * cur_mom_pints(addr_low_hbra,addr_hket_cur,1)                  &
              + real(cur_order_mom,REALK)                                     &
              * low_mom_pints(addr_hbra_low,addr_hket_low,1))
          end do
          ! component y...yz...z on bra center
          addr_hbra_up = addr_hbra_up+1
          addr_hbra_cur = addr_hbra_cur+1
          addr_hbra_low = addr_hbra_low+1
          up_mom_pints(addr_hbra_up,addr_hket_up,1)                         &
            = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,addr_hket_cur,1) &
            + hrp_total_expnt*real(cur_order_mom,REALK)                     &
            * low_mom_pints(addr_hbra_low,addr_hket_low,1)
        end do
        ! component z...z on bra center
        addr_hbra_up = addr_hbra_up+1
        addr_hbra_cur = addr_hbra_cur+1
        addr_hbra_low = addr_hbra_low+1
        up_mom_pints(addr_hbra_up,addr_hket_up,1)                         &
          = cc_wrt_diporg(1)*cur_mom_pints(addr_hbra_cur,addr_hket_cur,1) &
          + hrp_total_expnt*real(cur_order_mom,REALK)                     &
          * low_mom_pints(addr_hbra_low,addr_hket_low,1)
      end do
    end do
    ! (2) x...xy component of Cartesian multipole moment
    addr_up_mom = 2
    addr_cur_mom = 1
    ! (2.1) zeroth order HGTO on ket center if required
    if (zero_hket) then
      addr_hket_up = 1
      ! recovers the zeroth order HGTO on bra center
      if (zero_hbra) then
        addr_hbra_up = 1
        up_mom_pints(1,1,addr_up_mom) &
          = cc_wrt_diporg(2)*cur_mom_pints(1,1,addr_cur_mom)
      else
        addr_hbra_up = 0
      end if
      ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
      addr_hbra_cur = base_hbra_cur
      ! initializes the address of lower order HGTOs on bra center
      addr_low_hbra = 0
      ! other order (>0) HGTOs on bra center
      do order_hbra = min_order_hbra, orders_hbra_up(2)
        do ibra = order_hbra, 1, -1
          ! component x...xz...z on bra center
          addr_hbra_up = addr_hbra_up+1
          addr_hbra_cur = addr_hbra_cur+1
          up_mom_pints(addr_hbra_up,1,addr_up_mom) &
            = cc_wrt_diporg(2)*cur_mom_pints(addr_hbra_cur,1,addr_cur_mom)
          ! components x...xyz...z to y...yz..z on bra center
          do jbra = 1, ibra
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            addr_low_hbra = addr_low_hbra+1
            up_mom_pints(addr_hbra_up,1,addr_up_mom)                         &
              = cc_wrt_diporg(2)*cur_mom_pints(addr_hbra_cur,1,addr_cur_mom) &
              + hrp_total_expnt*real(jbra,REALK)                             &
              * cur_mom_pints(addr_low_hbra,1,addr_cur_mom)
          end do
        end do
        ! component z...z on bra center
        addr_hbra_up = addr_hbra_up+1
        addr_hbra_cur = addr_hbra_cur+1
        up_mom_pints(addr_hbra_up,1,addr_up_mom) &
          = cc_wrt_diporg(2)*cur_mom_pints(addr_hbra_cur,1,addr_cur_mom)
      end do
    else
      addr_hket_up = 0
    end if
    ! initializes the address of HGTOs on ket center in current order Cartesian multipole moments
    addr_hket_cur = base_hket_cur
    ! initializes the address of lower order HGTOs on ket center
    addr_low_hket = 0
    ! (2.2) other order (>0) HGTOs on ket center
    do order_hket = min_order_hket, orders_hket_up(2)
      do iket = order_hket, 1, -1
        ! component x...xz...z on ket center
        addr_hket_up = addr_hket_up+1
        addr_hket_cur = addr_hket_cur+1
        ! recovers the zeroth order HGTO on bra center
        if (zero_hbra) then
          addr_hbra_up = 1
          up_mom_pints(1,addr_hket_up,addr_up_mom) &
            = cc_wrt_diporg(2)*cur_mom_pints(1,addr_hket_cur,addr_cur_mom)
        else
          addr_hbra_up = 0
        end if
        ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
        addr_hbra_cur = base_hbra_cur
        ! initializes the address of lower order HGTOs on bra center
        addr_low_hbra = 0
        ! other order (>0) HGTOs on bra center
        do order_hbra = min_order_hbra, orders_hbra_up(2)
          do ibra = order_hbra, 1, -1
            ! component x...xz...z on bra center
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom) &
              = cc_wrt_diporg(2)                                &
              * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) 
            ! components x...xyz...z to y...yz..z on bra center
            do jbra = 1, ibra
              addr_hbra_up = addr_hbra_up+1
              addr_hbra_cur = addr_hbra_cur+1
              addr_low_hbra = addr_low_hbra+1
              up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
                = cc_wrt_diporg(2)                                        &
                * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
                + hrp_total_expnt*real(jbra,REALK)                        &
                * cur_mom_pints(addr_low_hbra,addr_hket_cur,addr_cur_mom)
            end do
          end do
          ! component z...z on bra center
          addr_hbra_up = addr_hbra_up+1
          addr_hbra_cur = addr_hbra_cur+1
          up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom) &
            = cc_wrt_diporg(2)                                &
            * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom)
        end do
        ! components x...xyz...z to y...yz..z on ket center
        do jket = 1, iket
          addr_hket_up = addr_hket_up+1
          addr_hket_cur = addr_hket_cur+1
          addr_low_hket = addr_low_hket+1
          ! recovers the zeroth order HGTO on bra center
          if (zero_hbra) then
            addr_hbra_up = 1
            up_mom_pints(1,addr_hket_up,addr_up_mom)                         &
              = cc_wrt_diporg(2)*cur_mom_pints(1,addr_hket_cur,addr_cur_mom) &
              + hrp_total_expnt*real(jket,REALK)                             &
              * cur_mom_pints(1,addr_low_hket,addr_cur_mom)
          else
            addr_hbra_up = 0
          end if
          ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
          addr_hbra_cur = base_hbra_cur
          ! initializes the address of lower order HGTOs on bra center
          addr_low_hbra = 0
          ! other order (>0) HGTOs on bra center
          do order_hbra = min_order_hbra, orders_hbra_up(2)
            do ibra = order_hbra, 1, -1
              ! component x...xz...z on bra center
              addr_hbra_up = addr_hbra_up+1
              addr_hbra_cur = addr_hbra_cur+1
              up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
                = cc_wrt_diporg(2)                                        &
                * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
                + hrp_total_expnt*real(jket,REALK)                        &
                * cur_mom_pints(addr_hbra_cur,addr_low_hket,addr_cur_mom)
              ! components x...xyz...z to y...yz..z on bra center
              do jbra = 1, ibra
                addr_hbra_up = addr_hbra_up+1
                addr_hbra_cur = addr_hbra_cur+1
                addr_low_hbra = addr_low_hbra+1
                up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
                  = cc_wrt_diporg(2)                                        &
                  * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
                  + hrp_total_expnt*(real(jbra,REALK)                       &
                  * cur_mom_pints(addr_low_hbra,addr_hket_cur,addr_cur_mom) &
                  + real(jket,REALK)                                        &
                  * cur_mom_pints(addr_hbra_cur,addr_low_hket,addr_cur_mom))
              end do
            end do
            ! component z...z on bra center
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
              = cc_wrt_diporg(2)                                        &
              * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
              + hrp_total_expnt*real(jket,REALK)                        &
              * cur_mom_pints(addr_hbra_cur,addr_low_hket,addr_cur_mom)
          end do
        end do
      end do
      ! component z...z on ket center
      addr_hket_up = addr_hket_up+1
      addr_hket_cur = addr_hket_cur+1
      ! recovers the zeroth order HGTO on bra center
      if (zero_hbra) then
        addr_hbra_up = 1
        up_mom_pints(1,addr_hket_up,addr_up_mom) &
          = cc_wrt_diporg(2)*cur_mom_pints(1,addr_hket_cur,addr_cur_mom)
      else
        addr_hbra_up = 0
      end if
      ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
      addr_hbra_cur = base_hbra_cur
      ! initializes the address of lower order HGTOs on bra center
      addr_low_hbra = 0
      ! other order (>0) HGTOs on bra center
      do order_hbra = min_order_hbra, orders_hbra_up(2)
        do ibra = order_hbra, 1, -1
          ! component x...xz...z on bra center
          addr_hbra_up = addr_hbra_up+1
          addr_hbra_cur = addr_hbra_cur+1
          up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom) &
            = cc_wrt_diporg(2)                                &
            * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom)
          ! components x...xyz...z to y...yz..z on bra center
          do jbra = 1, ibra
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            addr_low_hbra = addr_low_hbra+1
            up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
              = cc_wrt_diporg(2)                                        &
              * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
              + hrp_total_expnt*real(jbra,REALK)                        &
              * cur_mom_pints(addr_low_hbra,addr_hket_cur,addr_cur_mom)
          end do
        end do
        ! component z...z on bra center
        addr_hbra_up = addr_hbra_up+1
        addr_hbra_cur = addr_hbra_cur+1
        up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom) &
          = cc_wrt_diporg(2)                                &
          * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom)
      end do
    end do
    ! (3) x...xyy to y...y components of Cartesian multipole moment
    addr_low_mom = 0
    do imom = 1, cur_order_mom
      addr_up_mom = addr_up_mom+1
      addr_cur_mom = addr_cur_mom+1
      addr_low_mom = addr_low_mom+1
      ! (3.1) zeroth order HGTO on ket center if required
      if (zero_hket) then
        addr_hket_up = 1
        addr_hket_low = offset_hket_low+1
        ! recovers the zeroth order HGTO on bra center
        if (zero_hbra) then
          addr_hbra_up = 1
          addr_hbra_low = offset_hbra_low+1
          up_mom_pints(1,1,addr_up_mom)                        &
            = cc_wrt_diporg(2)*cur_mom_pints(1,1,addr_cur_mom) &
            + hrp_total_expnt*real(imom,REALK)                 &
            * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom)
        else
          addr_hbra_up = 0
          addr_hbra_low = offset_hbra_low
        end if
        ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
        addr_hbra_cur = base_hbra_cur
        ! initializes the address of lower order HGTOs on bra center
        addr_low_hbra = 0
        ! other order (>0) HGTOs on bra center
        do order_hbra = min_order_hbra, orders_hbra_up(2)
          do ibra = order_hbra, 1, -1
            ! component x...xz...z on bra center
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            addr_hbra_low = addr_hbra_low+1
            up_mom_pints(addr_hbra_up,1,addr_up_mom)                         &
              = cc_wrt_diporg(2)*cur_mom_pints(addr_hbra_cur,1,addr_cur_mom) &
              + hrp_total_expnt*real(imom,REALK)                             &
              * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom)
            ! components x...xyz...z to y...yz..z on bra center
            do jbra = 1, ibra
              addr_hbra_up = addr_hbra_up+1
              addr_hbra_cur = addr_hbra_cur+1
              addr_hbra_low = addr_hbra_low+1
              addr_low_hbra = addr_low_hbra+1
              up_mom_pints(addr_hbra_up,1,addr_up_mom)                         &
                = cc_wrt_diporg(2)*cur_mom_pints(addr_hbra_cur,1,addr_cur_mom) &
                + hrp_total_expnt*(real(jbra,REALK)                            &
                * cur_mom_pints(addr_low_hbra,1,addr_cur_mom)                  &
                + real(imom,REALK)                                             &
                * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom))
            end do
          end do
          ! component z...z on bra center
          addr_hbra_up = addr_hbra_up+1
          addr_hbra_cur = addr_hbra_cur+1
          addr_hbra_low = addr_hbra_low+1
          up_mom_pints(addr_hbra_up,1,addr_up_mom)                         &
            = cc_wrt_diporg(2)*cur_mom_pints(addr_hbra_cur,1,addr_cur_mom) &
            + hrp_total_expnt*real(imom,REALK)                             &
            * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom)
        end do
      else
        addr_hket_up = 0
        addr_hket_low = offset_hket_low
      end if
      ! initializes the address of HGTOs on ket center in current order Cartesian multipole moments
      addr_hket_cur = base_hket_cur
      ! initializes the address of lower order HGTOs on ket center
      addr_low_hket = 0
      ! (3.2) other order (>0) HGTOs on ket center
      do order_hket = min_order_hket, orders_hket_up(2)
        do iket = order_hket, 1, -1
          ! component x...xz...z on ket center
          addr_hket_up = addr_hket_up+1
          addr_hket_cur = addr_hket_cur+1
          addr_hket_low = addr_hket_low+1
          ! recovers the zeroth order HGTO on bra center
          if (zero_hbra) then
            addr_hbra_up = 1
            addr_hbra_low = offset_hbra_low+1
            up_mom_pints(1,addr_hket_up,addr_up_mom)                         &
              = cc_wrt_diporg(2)*cur_mom_pints(1,addr_hket_cur,addr_cur_mom) &
              + hrp_total_expnt*real(imom,REALK)                             &
              * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom)
          else
            addr_hbra_up = 0
            addr_hbra_low = offset_hbra_low
          end if
          ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
          addr_hbra_cur = base_hbra_cur
          ! initializes the address of lower order HGTOs on bra center
          addr_low_hbra = 0
          ! other order (>0) HGTOs on bra center
          do order_hbra = min_order_hbra, orders_hbra_up(2)
            do ibra = order_hbra, 1, -1
              ! component x...xz...z on bra center
              addr_hbra_up = addr_hbra_up+1
              addr_hbra_cur = addr_hbra_cur+1
              addr_hbra_low = addr_hbra_low+1
              up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
                = cc_wrt_diporg(2)                                        &
                * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
                + hrp_total_expnt*real(imom,REALK)                        &
                * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom)
              ! components x...xyz...z to y...yz..z on bra center
              do jbra = 1, ibra
                addr_hbra_up = addr_hbra_up+1
                addr_hbra_cur = addr_hbra_cur+1
                addr_hbra_low = addr_hbra_low+1
                addr_low_hbra = addr_low_hbra+1
                up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
                  = cc_wrt_diporg(2)                                        &
                  * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
                  + hrp_total_expnt*(real(jbra,REALK)                       &
                  * cur_mom_pints(addr_low_hbra,addr_hket_cur,addr_cur_mom) &
                  + real(imom,REALK)                                        &
                  * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom))
              end do
            end do
            ! component z...z on bra center
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            addr_hbra_low = addr_hbra_low+1
            up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
              = cc_wrt_diporg(2)                                        &
              * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
              + hrp_total_expnt*real(imom,REALK)                        &
              * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom)
          end do
          ! components x...xyz...z to y...yz..z on ket center
          do jket = 1, iket
            addr_hket_up = addr_hket_up+1
            addr_hket_cur = addr_hket_cur+1
            addr_hket_low = addr_hket_low+1
            addr_low_hket = addr_low_hket+1
            ! recovers the zeroth order HGTO on bra center
            if (zero_hbra) then
              addr_hbra_up = 1
              addr_hbra_low = offset_hbra_low+1
              up_mom_pints(1,addr_hket_up,addr_up_mom)                         &
                = cc_wrt_diporg(2)*cur_mom_pints(1,addr_hket_cur,addr_cur_mom) &
                + hrp_total_expnt*(real(jket,REALK)                            &
                * cur_mom_pints(1,addr_low_hket,addr_cur_mom)                  &
                + real(imom,REALK)                                             &
                * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom))
            else
              addr_hbra_up = 0
              addr_hbra_low = offset_hbra_low
            end if
            ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
            addr_hbra_cur = base_hbra_cur
            ! initializes the address of lower order HGTOs on bra center
            addr_low_hbra = 0
            ! other order (>0) HGTOs on bra center
            do order_hbra = min_order_hbra, orders_hbra_up(2)
              do ibra = order_hbra, 1, -1
                ! component x...xz...z on bra center
                addr_hbra_up = addr_hbra_up+1
                addr_hbra_cur = addr_hbra_cur+1
                addr_hbra_low = addr_hbra_low+1
                up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
                  = cc_wrt_diporg(2)                                        &
                  * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
                  + hrp_total_expnt*(real(jket,REALK)                       &
                  * cur_mom_pints(addr_hbra_cur,addr_low_hket,addr_cur_mom) &
                  + real(imom,REALK)                                        &
                  * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom))
                ! components x...xyz...z to y...yz..z on bra center
                do jbra = 1, ibra
                  addr_hbra_up = addr_hbra_up+1
                  addr_hbra_cur = addr_hbra_cur+1
                  addr_hbra_low = addr_hbra_low+1
                  addr_low_hbra = addr_low_hbra+1
                  up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
                    = cc_wrt_diporg(2)                                        &
                    * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
                    + hrp_total_expnt*(real(jbra,REALK)                       &
                    * cur_mom_pints(addr_low_hbra,addr_hket_cur,addr_cur_mom) &
                    + real(jket,REALK)                                        &
                    * cur_mom_pints(addr_hbra_cur,addr_low_hket,addr_cur_mom) &
                    + real(imom,REALK)                                        &
                    * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom))
                end do
              end do
              ! component z...z on bra center
              addr_hbra_up = addr_hbra_up+1
              addr_hbra_cur = addr_hbra_cur+1
              addr_hbra_low = addr_hbra_low+1
              up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
                = cc_wrt_diporg(2)                                        &
                * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
                + hrp_total_expnt*(real(jket,REALK)                       &
                * cur_mom_pints(addr_hbra_cur,addr_low_hket,addr_cur_mom) &
                + real(imom,REALK)                                        &
                * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom))
            end do
          end do
        end do
        ! component z...z on ket center
        addr_hket_up = addr_hket_up+1
        addr_hket_cur = addr_hket_cur+1
        addr_hket_low = addr_hket_low+1
        ! recovers the zeroth order HGTO on bra center
        if (zero_hbra) then
          addr_hbra_up = 1
          addr_hbra_low = offset_hbra_low+1
          up_mom_pints(1,addr_hket_up,addr_up_mom)                         &
            = cc_wrt_diporg(2)*cur_mom_pints(1,addr_hket_cur,addr_cur_mom) &
            + hrp_total_expnt*real(imom,REALK)                             &
            * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom)
        else
          addr_hbra_up = 0
          addr_hbra_low = offset_hbra_low
        end if
        ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
        addr_hbra_cur = base_hbra_cur
        ! initializes the address of lower order HGTOs on bra center
        addr_low_hbra = 0
        ! other order (>0) HGTOs on bra center
        do order_hbra = min_order_hbra, orders_hbra_up(2)
          do ibra = order_hbra, 1, -1
            ! component x...xz...z on bra center
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            addr_hbra_low = addr_hbra_low+1
            up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
              = cc_wrt_diporg(2)                                        &
              * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
              + hrp_total_expnt*real(imom,REALK)                        &
              * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom)
            ! components x...xyz...z to y...yz..z on bra center
            do jbra = 1, ibra
              addr_hbra_up = addr_hbra_up+1
              addr_hbra_cur = addr_hbra_cur+1
              addr_hbra_low = addr_hbra_low+1
              addr_low_hbra = addr_low_hbra+1
              up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
                = cc_wrt_diporg(2)                                        &
                * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
                + hrp_total_expnt*(real(jbra,REALK)                       &
                * cur_mom_pints(addr_low_hbra,addr_hket_cur,addr_cur_mom) &
                + real(imom,REALK)                                        &
                * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom))
            end do
          end do
          ! component z...z on bra center
          addr_hbra_up = addr_hbra_up+1
          addr_hbra_cur = addr_hbra_cur+1
          addr_hbra_low = addr_hbra_low+1
          up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
            = cc_wrt_diporg(2)                                        &
            * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
            + hrp_total_expnt*real(imom,REALK)                        &
            * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom)
        end do
      end do
    end do
    ! (4) x...xz to y...yz components of Cartesian multipole moment
    addr_cur_mom = 0
    do jmom = 0, cur_order_mom
      addr_up_mom = addr_up_mom+1
      addr_cur_mom = addr_cur_mom+1
      ! (4.1) zeroth order HGTO on ket center if required
      if (zero_hket) then
        addr_hket_up = 1
        ! recovers the zeroth order HGTO on bra center
        if (zero_hbra) then
          addr_hbra_up = 1
          up_mom_pints(1,1,addr_up_mom) &
            = cc_wrt_diporg(3)*cur_mom_pints(1,1,addr_cur_mom)
        else
          addr_hbra_up = 0
        end if
        ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
        addr_hbra_cur = base_hbra_cur
        ! initializes the address of lower order HGTOs on bra center
        addr_low_hbra = 0
        ! other order (>0) HGTOs on bra center
        do order_hbra = min_order_hbra, orders_hbra_up(2)
          ! components x...x to y...y on bra center
          do ibra = 0, order_hbra
            addr_hbra_up = addr_hbra_up+1
            addr_hbra_cur = addr_hbra_cur+1
            up_mom_pints(addr_hbra_up,1,addr_up_mom) &
              = cc_wrt_diporg(3)*cur_mom_pints(addr_hbra_cur,1,addr_cur_mom)
          end do
          ! components x...xz to z...z on bra center
          do ibra = 1, order_hbra
            do jbra = 0, order_hbra-ibra
              addr_hbra_up = addr_hbra_up+1
              addr_hbra_cur = addr_hbra_cur+1
              addr_low_hbra = addr_low_hbra+1
              up_mom_pints(addr_hbra_up,1,addr_up_mom)                         &
                = cc_wrt_diporg(3)*cur_mom_pints(addr_hbra_cur,1,addr_cur_mom) &
                + hrp_total_expnt*real(ibra,REALK)                             &
                * cur_mom_pints(addr_low_hbra,1,addr_cur_mom)
            end do
          end do
        end do
      else
        addr_hket_up = 0
      end if
      ! initializes the address of HGTOs on ket center in current order Cartesian multipole moments
      addr_hket_cur = base_hket_cur
      ! initializes the address of lower order HGTOs on ket center
      addr_low_hket = 0
      ! (4.2) other order (>0) HGTOs on ket center
      do order_hket = min_order_hket, orders_hket_up(2)
        ! components x...x to y...y on ket center
        do iket = 0, order_hket
          addr_hket_up = addr_hket_up+1
          addr_hket_cur = addr_hket_cur+1
          ! recovers the zeroth order HGTO on bra center
          if (zero_hbra) then
            addr_hbra_up = 1
            up_mom_pints(1,addr_hket_up,addr_up_mom) &
              = cc_wrt_diporg(3)*cur_mom_pints(1,addr_hket_cur,addr_cur_mom)
          else
            addr_hbra_up = 0
          end if
          ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
          addr_hbra_cur = base_hbra_cur
          ! initializes the address of lower order HGTOs on bra center
          addr_low_hbra = 0
          ! other order (>0) HGTOs on bra center
          do order_hbra = min_order_hbra, orders_hbra_up(2)
            ! components x...x to y...y on bra center
            do ibra = 0, order_hbra
              addr_hbra_up = addr_hbra_up+1
              addr_hbra_cur = addr_hbra_cur+1
              up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom) &
                = cc_wrt_diporg(3)                                &
                * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom)
            end do
            ! components x...xz to z...z on bra center
            do ibra = 1, order_hbra
              do jbra = 0, order_hbra-ibra
                addr_hbra_up = addr_hbra_up+1
                addr_hbra_cur = addr_hbra_cur+1
                addr_low_hbra = addr_low_hbra+1
                up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
                  = cc_wrt_diporg(3)                                        &
                  * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
                  + hrp_total_expnt*real(ibra,REALK)                        &
                  * cur_mom_pints(addr_low_hbra,addr_hket_cur,addr_cur_mom)
              end do
            end do
          end do
        end do
        ! components x...xz to z...z on ket center
        do iket = 1, order_hket
          do jket = 0, order_hket-iket
            addr_hket_up = addr_hket_up+1
            addr_hket_cur = addr_hket_cur+1
            addr_low_hket = addr_low_hket+1
            ! recovers the zeroth order HGTO on bra center
            if (zero_hbra) then
              addr_hbra_up = 1
              up_mom_pints(1,addr_hket_up,addr_up_mom)                         &
                = cc_wrt_diporg(3)*cur_mom_pints(1,addr_hket_cur,addr_cur_mom) &
                + hrp_total_expnt*real(iket,REALK)                             &
                * cur_mom_pints(1,addr_low_hket,addr_cur_mom)
            else
              addr_hbra_up = 0
            end if
            ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
            addr_hbra_cur = base_hbra_cur
            ! initializes the address of lower order HGTOs on bra center
            addr_low_hbra = 0
            ! other order (>0) HGTOs on bra center
            do order_hbra = min_order_hbra, orders_hbra_up(2)
              ! components x...x to y...y on bra center
              do ibra = 0, order_hbra
                addr_hbra_up = addr_hbra_up+1
                addr_hbra_cur = addr_hbra_cur+1
                up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
                  = cc_wrt_diporg(3)                                        &
                  * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
                  + hrp_total_expnt*real(iket,REALK)                        &
                  * cur_mom_pints(addr_hbra_cur,addr_low_hket,addr_cur_mom)
              end do
              ! components x...xz to z...z on bra center
              do ibra = 1, order_hbra
                do jbra = 0, order_hbra-ibra
                  addr_hbra_up = addr_hbra_up+1
                  addr_hbra_cur = addr_hbra_cur+1
                  addr_low_hbra = addr_low_hbra+1
                  up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
                    = cc_wrt_diporg(3)                                        &
                    * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
                    + hrp_total_expnt*(real(ibra,REALK)                       &
                    * cur_mom_pints(addr_low_hbra,addr_hket_cur,addr_cur_mom) &
                    + real(iket,REALK)                                        &
                    * cur_mom_pints(addr_hbra_cur,addr_low_hket,addr_cur_mom))
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    ! (5) x...xzz to z...z components of Cartesian multipole moment
    addr_low_mom = 0
    do imom = 1, cur_order_mom
      do jmom = 0, cur_order_mom-imom
        addr_up_mom = addr_up_mom+1
        addr_cur_mom = addr_cur_mom+1
        addr_low_mom = addr_low_mom+1
        ! (5.1) zeroth order HGTO on ket center if required
        if (zero_hket) then
          addr_hket_up = 1
          addr_hket_low = offset_hket_low+1
          ! recovers the zeroth order HGTO on bra center
          if (zero_hbra) then
            addr_hbra_up = 1
            addr_hbra_low = offset_hbra_low+1
            up_mom_pints(1,1,addr_up_mom)                        &
              = cc_wrt_diporg(3)*cur_mom_pints(1,1,addr_cur_mom) &
              + hrp_total_expnt*real(imom,REALK)                 &
              * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom)
          else
            addr_hbra_up = 0
            addr_hbra_low = offset_hbra_low
          end if
          ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
          addr_hbra_cur = base_hbra_cur
          ! initializes the address of lower order HGTOs on bra center
          addr_low_hbra = 0
          ! other order (>0) HGTOs on bra center
          do order_hbra = min_order_hbra, orders_hbra_up(2)
            ! components x...x to y...y on bra center
            do ibra = 0, order_hbra
              addr_hbra_up = addr_hbra_up+1
              addr_hbra_cur = addr_hbra_cur+1
              addr_hbra_low = addr_hbra_low+1
              up_mom_pints(addr_hbra_up,1,addr_up_mom)                         &
                = cc_wrt_diporg(3)*cur_mom_pints(addr_hbra_cur,1,addr_cur_mom) &
                + hrp_total_expnt*real(imom,REALK)                             &
                * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom)
            end do
            ! components x...xz to z...z on bra center
            do ibra = 1, order_hbra
              do jbra = 0, order_hbra-ibra
                addr_hbra_up = addr_hbra_up+1
                addr_hbra_cur = addr_hbra_cur+1
                addr_hbra_low = addr_hbra_low+1
                addr_low_hbra = addr_low_hbra+1
                up_mom_pints(addr_hbra_up,1,addr_up_mom)                         &
                  = cc_wrt_diporg(3)*cur_mom_pints(addr_hbra_cur,1,addr_cur_mom) &
                  + hrp_total_expnt*(real(ibra,REALK)                            &
                  * cur_mom_pints(addr_low_hbra,1,addr_cur_mom)                  &
                  + real(imom,REALK)                                             &
                  * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom))
              end do
            end do
          end do
        else
          addr_hket_up = 0
          addr_hket_low = offset_hket_low
        end if
        ! initializes the address of HGTOs on ket center in current order Cartesian multipole moments
        addr_hket_cur = base_hket_cur
        ! initializes the address of lower order HGTOs on ket center
        addr_low_hket = 0
        ! (5.2) other order (>0) HGTOs on ket center
        do order_hket = min_order_hket, orders_hket_up(2)
          ! components x...x to y...y on ket center
          do iket = 0, order_hket
            addr_hket_up = addr_hket_up+1
            addr_hket_cur = addr_hket_cur+1
            addr_hket_low = addr_hket_low+1
            ! recovers the zeroth order HGTO on bra center
            if (zero_hbra) then
              addr_hbra_up = 1
              addr_hbra_low = offset_hbra_low+1
              up_mom_pints(1,addr_hket_up,addr_up_mom)                         &
                = cc_wrt_diporg(3)*cur_mom_pints(1,addr_hket_cur,addr_cur_mom) &
                + hrp_total_expnt*real(imom,REALK)                             &
                * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom)
            else
              addr_hbra_up = 0
              addr_hbra_low = offset_hbra_low
            end if
            ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
            addr_hbra_cur = base_hbra_cur
            ! initializes the address of lower order HGTOs on bra center
            addr_low_hbra = 0
            ! other order (>0) HGTOs on bra center
            do order_hbra = min_order_hbra, orders_hbra_up(2)
              ! components x...x to y...y on bra center
              do ibra = 0, order_hbra
                addr_hbra_up = addr_hbra_up+1
                addr_hbra_cur = addr_hbra_cur+1
                addr_hbra_low = addr_hbra_low+1
                up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
                  = cc_wrt_diporg(3)                                        &
                  * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
                  + hrp_total_expnt*real(imom,REALK)                        &
                  * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom)
              end do
              ! components x...xz to z...z on bra center
              do ibra = 1, order_hbra
                do jbra = 0, order_hbra-ibra
                  addr_hbra_up = addr_hbra_up+1
                  addr_hbra_cur = addr_hbra_cur+1
                  addr_hbra_low = addr_hbra_low+1
                  addr_low_hbra = addr_low_hbra+1
                  up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
                    = cc_wrt_diporg(3)                                        &
                    * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
                    + hrp_total_expnt*(real(ibra,REALK)                       &
                    * cur_mom_pints(addr_low_hbra,addr_hket_cur,addr_cur_mom) &
                    + real(imom,REALK)                                        &
                    * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom))
                end do
              end do
            end do
          end do
          ! components x...xz to z...z on ket center
          do iket = 1, order_hket
            do jket = 0, order_hket-iket
              addr_hket_up = addr_hket_up+1
              addr_hket_cur = addr_hket_cur+1
              addr_hket_low = addr_hket_low+1
              addr_low_hket = addr_low_hket+1
              ! recovers the zeroth order HGTO on bra center
              if (zero_hbra) then
                addr_hbra_up = 1
                addr_hbra_low = offset_hbra_low+1
                up_mom_pints(1,addr_hket_up,addr_up_mom)                         &
                  = cc_wrt_diporg(3)*cur_mom_pints(1,addr_hket_cur,addr_cur_mom) &
                  + hrp_total_expnt*(real(iket,REALK)                            &
                  * cur_mom_pints(1,addr_low_hket,addr_cur_mom)                  &
                  + real(imom,REALK)                                             &
                  * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom))
              else
                addr_hbra_up = 0
                addr_hbra_low = offset_hbra_low
              end if
              ! initializes the address of HGTOs on bra center in current order Cartesian multipole moments
              addr_hbra_cur = base_hbra_cur
              ! initializes the address of lower order HGTOs on bra center
              addr_low_hbra = 0
              ! other order (>0) HGTOs on bra center
              do order_hbra = min_order_hbra, orders_hbra_up(2)
                ! components x...x to y...y on bra center
                do ibra = 0, order_hbra
                  addr_hbra_up = addr_hbra_up+1
                  addr_hbra_cur = addr_hbra_cur+1
                  addr_hbra_low = addr_hbra_low+1
                  up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
                    = cc_wrt_diporg(3)                                        &
                    * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
                    + hrp_total_expnt*(real(iket,REALK)                       &
                    * cur_mom_pints(addr_hbra_cur,addr_low_hket,addr_cur_mom) &
                    + real(imom,REALK)                                        &
                    * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom))
                end do
                ! components x...xz to z...z on bra center
                do ibra = 1, order_hbra
                  do jbra = 0, order_hbra-ibra
                    addr_hbra_up = addr_hbra_up+1
                    addr_hbra_cur = addr_hbra_cur+1
                    addr_hbra_low = addr_hbra_low+1
                    addr_low_hbra = addr_low_hbra+1
                    up_mom_pints(addr_hbra_up,addr_hket_up,addr_up_mom)         &
                      = cc_wrt_diporg(3)                                        &
                      * cur_mom_pints(addr_hbra_cur,addr_hket_cur,addr_cur_mom) &
                      + hrp_total_expnt*(real(ibra,REALK)                       &
                      * cur_mom_pints(addr_low_hbra,addr_hket_cur,addr_cur_mom) &
                      + real(iket,REALK)                                        &
                      * cur_mom_pints(addr_hbra_cur,addr_low_hket,addr_cur_mom) &
                      + real(imom,REALK)                                        &
                      * low_mom_pints(addr_hbra_low,addr_hket_low,addr_low_mom))
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "sub_carmom_moment", STDOUT)
#endif
    return
  end subroutine sub_carmom_moment
