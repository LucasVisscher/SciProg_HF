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
!!  This file recovers Cartesian multipole moments in Dirac delta function integrals.
!!
!!  2012-03-16, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief recovers Cartesian multipole moments in Dirac delta function integrals
  !> \author Bin Gao
  !> \date 2012-03-16
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_geo_pot is the order of geometric derivatives on Dirac delta function
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param delta_origin contains the coordinates of Dirac delta function origin
  !> \param dim_hgto_bra is the dimension of HGTOs on bra center
  !> \param dim_hgto_ket is the dimension of Hermite Gaussians on ket center
  !> \param dim_geo_hket is the dimension of geometric derivatives on Dirac delta function
  !> \param hket_pints contains the primitive Hermite integrals with zeroth order Cartesian
  !>        multipole moment
  !> \param num_mom is the number of xyz components of Cartesian multipole moment,
  !>        equals to \f$(\var(order_mom)+1)(\var(order_mom)+2)/2\f$
  !> \param num_geo_pot is the number of geometric derivatives on potential origin,
  !>        equals to \f$(\var(order_geo_pot)+1)(\var(order_geo_pot)+2)/2\f$
  !> \return hmom_pints contains the primitive Hermite integrals with specified orders of
  !>         Cartesian multipole moments and geometric derivatives on Dirac delta function
  subroutine delta_moment(order_mom, order_geo_pot, dipole_origin, delta_origin, &
                          dim_hgto_bra, dim_hgto_ket, dim_geo_hket, hket_pints,  &
                          num_mom, num_geo_pot, hmom_pints)
    use xkind
    implicit none
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_geo_pot
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: delta_origin(3)
    integer, intent(in) :: dim_hgto_bra
    integer, intent(in) :: dim_hgto_ket
    integer, intent(in) :: dim_geo_hket
    real(REALK), intent(in) :: hket_pints(dim_hgto_bra,dim_hgto_ket,dim_geo_hket)
    integer, intent(in) :: num_mom
    integer, intent(in) :: num_geo_pot
    real(REALK), intent(out) :: hmom_pints(dim_hgto_bra,dim_hgto_ket,num_mom,num_geo_pot)
!f2py intent(in) :: order_mom
!f2py intent(in) :: order_geo_pot
!f2py intent(in) :: dipole_origin
!f2py intent(in) :: delta_origin
!f2py intent(hide) :: dim_hgto_bra
!f2py intent(hide) :: dim_hgto_ket
!f2py intent(hide) :: dim_geo_hket
!f2py intent(in) :: hket_pints
!f2py intent(in) :: num_mom
!f2py intent(in) :: num_geo_pot
!f2py intent(out) :: hmom_pints
!f2py depend(dim_hgto_bra) :: hmom_pints
!f2py depend(dim_hgto_ket) :: hmom_pints
!f2py depend(num_mom) :: hmom_pints
!f2py depend(num_geo_pot) :: hmom_pints
    integer min_geo_pot              !minimum geometric derivatives on Dirac delta function
    integer max_cur_mom              !maximum of current order Cartesian multipole moments
    real(REALK) delta_wrt_diporg(3)  !relative coordinates of Dirac delta function origin w.r.t. dipole origin
    integer dim_cur_mom              !dimension of current order Cartesian multipole moments
    integer dim_up_mom               !dimension of upper order Cartesian multipole moments
    integer num_low_geo              !number of xyz components of lower order geometric derivatives
    integer num_cur_geo              !number of xyz components of current order geometric derivatives
    integer size_low_geo             !size of temporary integrals of lower order geometric derivatives
    integer size_cur_geo             !size of temporary integrals of current order geometric derivatives
    integer low_geo_int              !pointer of lower order geometric derivatives
    integer cur_geo_int              !pointer of current order geometric derivatives
    real(REALK), allocatable :: tmp_ints(:,:,:,:)
                                     !temporary integrals
    integer start_low_hket           !start address of lower order geometric derivatives in \var(hket_pints)
    integer end_low_hket             !end address of lower order geometric derivatives in \var(hket_pints)
    integer start_cur_hket           !start address of current order geometric derivatives in \var(hket_pints)
    integer end_cur_hket             !end address of current order geometric derivatives in \var(hket_pints)
    integer cur_order_geo            !incremental recorder over current order of geometric derivatives
    integer ierr                     !error information
#if defined(XTIME)
    real(REALK) curr_time            !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! only zeroth order Cartesian multipole moment returns
    if (order_mom==0) then
      hmom_pints(:,:,1,:) = hket_pints
    else
      ! sets the minimum geometric derivatives on Dirac delta function
      min_geo_pot = order_geo_pot-order_mom+1
      ! sets the maximum of current order Cartesian multipole moments
      max_cur_mom = -min_geo_pot
      ! sets the relative coordinates of Dirac delta function origin w.r.t. dipole origin
      delta_wrt_diporg = delta_origin-dipole_origin
      ! allocates memory for temporary integrals
      dim_up_mom = (order_mom+1)*(order_mom+2)*(order_mom+3)/6-1
      num_cur_geo = (order_geo_pot+1)*(order_geo_pot+2)/2
      allocate(tmp_ints(dim_hgto_bra,dim_hgto_ket,dim_up_mom*num_cur_geo,2), stat=ierr)
      if (ierr/=0)                                            &
        call error_stop("delta_moment", "failed to tmp_ints", &
                        dim_hgto_bra*dim_hgto_ket*dim_up_mom*num_cur_geo*2)
      ! gets temporary integrals with zeroth order geometric derivative on Dirac delta function
      if (min_geo_pot<=0) then
        ! sets the dimension of upper order Cartesian multipole moments
        dim_up_mom = (max_cur_mom+2)*(max_cur_mom+3)*(max_cur_mom+4)/6-1
        call zero_delta_moment(max_cur_mom, delta_wrt_diporg, dim_hgto_bra, &
                               dim_hgto_ket, hket_pints(:,:,1), dim_up_mom, &
                               tmp_ints(:,:,1:dim_up_mom,1))
        ! sets the minimum geometric derivatives on Dirac delta function
        min_geo_pot = 1
        ! sets number of xyz components of current order geometric derivatives
        num_cur_geo = 1
        ! sets the start and end addresses of current order geometric derivatives in \var(hket_pints)
        start_cur_hket = 1
        end_cur_hket = 1
        ! sets the size of temporary integrals of current order geometric derivatives
        size_cur_geo = dim_up_mom
        ! initializes the pointers
        low_geo_int = 2
        cur_geo_int = 1
      else
        max_cur_mom = -1
        dim_up_mom = 0
        num_cur_geo = min_geo_pot*(min_geo_pot+1)/2
        start_cur_hket = 1
        end_cur_hket = num_cur_geo
        size_cur_geo = 0
        low_geo_int = 1
        cur_geo_int = 2
      end if
      ! loops over the orders of geometric derivatives on Dirac delta function
      do cur_order_geo = min_geo_pot, order_geo_pot
        ! updates the maximum of current order Cartesian multipole moments
        max_cur_mom = max_cur_mom+1
        ! updates the dimensions of Cartesian multipole moments
        dim_cur_mom = dim_up_mom
        dim_up_mom = dim_up_mom+(max_cur_mom+2)*(max_cur_mom+3)/2
        ! updates the number of geometric derivatives on Dirac delta function
        num_low_geo = num_cur_geo
        num_cur_geo = num_cur_geo+cur_order_geo+1
        ! updates the start and end addresses of geometric derivatives in \var(hket_pints)
        start_low_hket = start_cur_hket
        end_low_hket = end_cur_hket
        start_cur_hket = end_cur_hket+1
        end_cur_hket = end_cur_hket+num_cur_geo
        ! sets the size of temporary integrals
        size_low_geo = size_cur_geo
        size_cur_geo = dim_up_mom*num_cur_geo
        ! switch the pointers
        low_geo_int = 3-low_geo_int
        cur_geo_int = 3-cur_geo_int
        ! gets the temporary integrals with \var(cur_order_geo) order geometric derivatives
        ! on Dirac delta function
        call sub_delta_moment(max_cur_mom, cur_order_geo, delta_wrt_diporg,             &
                              dim_hgto_bra, dim_hgto_ket, num_low_geo,                  &
                              hket_pints(:,:,start_low_hket:end_low_hket),              &
                              dim_cur_mom, tmp_ints(:,:,1:size_low_geo,low_geo_int),    &
                              num_cur_geo, hket_pints(:,:,start_cur_hket:end_cur_hket), &
                              dim_up_mom, tmp_ints(:,:,1:size_cur_geo,cur_geo_int))
      end do
      ! assigns integrals of minimum order of geometric derivatives on Dirac delta function
      call delta_moment_assign(dim_hgto_bra, dim_hgto_ket, dim_up_mom, num_geo_pot, &
                               tmp_ints(:,:,1:size_cur_geo,cur_geo_int), num_mom,   &
                               hmom_pints)
      ! cleans
      deallocate(tmp_ints)
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "delta_moment", STDOUT)
#endif
    return
  end subroutine delta_moment

  !> \brief recovers Cartesian multipole moments for the zeroth order geometric derivatives
  !>        on Dirac delta function
  !> \author Bin Gao
  !> \date 2012-03-16
  !> \param max_cur_mom is the maximum of current order Cartesian multipole moments
  !> \param delta_wrt_diporg contains the relative coordinates of Dirac delta function
  !>        w.r.t. dipole origin
  !> \param dim_hgto_bra is the dimension of HGTOs on bra center
  !> \param dim_hgto_ket is the dimension of HGTOs on ket center
  !> \param zero_geo_pints contains integrals of zeroth order geometric derivatives on
  !>        Dirac delta function and zeroth order Cartesian multipole moments
  !> \param dim_up_mom is the dimension of upper order Cartesian multipole moments
  !> \return cur_geo_pints contains the integrals with orders of Cartesian multipole
  !>         moments from 1 to \var(max_cur_mom)+1
  subroutine zero_delta_moment(max_cur_mom, delta_wrt_diporg, dim_hgto_bra, &
                               dim_hgto_ket, zero_geo_pints, dim_up_mom,    &
                               cur_geo_pints)
    use xkind
    implicit none
    integer, intent(in) :: max_cur_mom
    real(REALK), intent(in) :: delta_wrt_diporg(3)
    integer, intent(in) :: dim_hgto_bra
    integer, intent(in) :: dim_hgto_ket
    real(REALK), intent(in) :: zero_geo_pints(dim_hgto_bra,dim_hgto_ket)
    integer, intent(in) :: dim_up_mom
    real(REALK), intent(out) :: cur_geo_pints(dim_hgto_bra,dim_hgto_ket,dim_up_mom)
!f2py intent(in) :: max_cur_mom
!f2py intent(in) :: delta_wrt_diporg
!f2py intent(hide) :: dim_hgto_bra
!f2py intent(hide) :: dim_hgto_ket
!f2py intent(in) :: zero_geo_pints
!f2py intent(in) :: dim_up_mom
!f2py intent(out) :: cur_geo_pints
!f2py depend(dim_hgto_bra) :: cur_geo_pints
!f2py depend(dim_hgto_ket) :: cur_geo_pints
!f2py depend(dim_up_mom) :: cur_geo_pints
    integer base_cur_mom  !base address of current order Cartesian multipole moments
    integer addr_cur_mom  !address of current order Cartesian multipole moments
    integer addr_up_mom   !address of upper order Cartesian multipole moments
    integer order_mom     !incremental recorder over the orders of Cartesian multipole moments
    integer imom, jmom    !incremental recorders over xyz components of Cartesian multipole moments
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time   
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    select case(max_cur_mom)
    ! only the first order Cartesian multipole moments required
    case(0)
#if defined(DEBUG)
      write(STDOUT,100) "1st order Cartesian multipole moments return"
#endif
      cur_geo_pints(:,:,1) = delta_wrt_diporg(1)*zero_geo_pints  !px
      cur_geo_pints(:,:,2) = delta_wrt_diporg(2)*zero_geo_pints  !py
      cur_geo_pints(:,:,3) = delta_wrt_diporg(3)*zero_geo_pints  !pz
    ! the first and second order Cartesian multipole moments
    case(1)
#if defined(DEBUG)
      write(STDOUT,100) "1st and 2nd order Cartesian multipole moments return"
#endif
      cur_geo_pints(:,:,1) = delta_wrt_diporg(1)*zero_geo_pints  !px
      cur_geo_pints(:,:,2) = delta_wrt_diporg(2)*zero_geo_pints  !py
      cur_geo_pints(:,:,3) = delta_wrt_diporg(3)*zero_geo_pints  !pz
      cur_geo_pints(:,:,4) = delta_wrt_diporg(1)*cur_geo_pints(:,:,1)  !dxx
      cur_geo_pints(:,:,5) = delta_wrt_diporg(2)*cur_geo_pints(:,:,1)  !dxy
      cur_geo_pints(:,:,6) = delta_wrt_diporg(2)*cur_geo_pints(:,:,2)  !dyy
      cur_geo_pints(:,:,7) = delta_wrt_diporg(3)*cur_geo_pints(:,:,1)  !dxz
      cur_geo_pints(:,:,8) = delta_wrt_diporg(3)*cur_geo_pints(:,:,2)  !dyz
      cur_geo_pints(:,:,9) = delta_wrt_diporg(3)*cur_geo_pints(:,:,3)  !dzz
    ! the maximum order of Cartesian multipole moments required is > 2
    case default
#if defined(DEBUG)
      write(STDOUT,100) "higher order Cartesian multipole moments return", max_cur_mom
#endif
      cur_geo_pints(:,:,1) = delta_wrt_diporg(1)*zero_geo_pints  !px
      cur_geo_pints(:,:,2) = delta_wrt_diporg(2)*zero_geo_pints  !py
      cur_geo_pints(:,:,3) = delta_wrt_diporg(3)*zero_geo_pints  !pz
      cur_geo_pints(:,:,4) = delta_wrt_diporg(1)*cur_geo_pints(:,:,1)  !dxx
      cur_geo_pints(:,:,5) = delta_wrt_diporg(2)*cur_geo_pints(:,:,1)  !dxy
      cur_geo_pints(:,:,6) = delta_wrt_diporg(2)*cur_geo_pints(:,:,2)  !dyy
      cur_geo_pints(:,:,7) = delta_wrt_diporg(3)*cur_geo_pints(:,:,1)  !dxz
      cur_geo_pints(:,:,8) = delta_wrt_diporg(3)*cur_geo_pints(:,:,2)  !dyz
      cur_geo_pints(:,:,9) = delta_wrt_diporg(3)*cur_geo_pints(:,:,3)  !dzz
      ! initializes the (base) addresses of Cartesian multipole moments
      base_cur_mom = 3
      addr_up_mom = 9
      ! other order (>2) Cartesian multipole moments
      do order_mom = 2, max_cur_mom
        ! recurrence relation along x-direction
        addr_up_mom = addr_up_mom+1
        addr_cur_mom = base_cur_mom
        cur_geo_pints(:,:,addr_up_mom) &
          = delta_wrt_diporg(1)*cur_geo_pints(:,:,addr_cur_mom+1)
        ! recurrence relation along y-direction
        do imom = 0, order_mom
          addr_cur_mom = addr_cur_mom+1
          addr_up_mom = addr_up_mom+1
          cur_geo_pints(:,:,addr_up_mom) &
            = delta_wrt_diporg(2)*cur_geo_pints(:,:,addr_cur_mom)
        end do
        ! recurrence relation along z-direction
        addr_cur_mom = base_cur_mom
        do imom = 0, order_mom
          do jmom = 0, order_mom-imom
            addr_cur_mom = addr_cur_mom+1
            addr_up_mom = addr_up_mom+1
            cur_geo_pints(:,:,addr_up_mom) &
              = delta_wrt_diporg(3)*cur_geo_pints(:,:,addr_cur_mom)
          end do
        end do
        ! updates the base addresses
        base_cur_mom = addr_cur_mom
      end do
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "zero_delta_moment", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("zero_delta_moment>> ",A,I6)
#endif
  end subroutine zero_delta_moment

  !> \brief recovers Cartesian multipole moments for a given order (>0) of geometric derivatives
  !>        on Dirac delta function
  !> \author Bin Gao
  !> \date 2012-03-16
  !> \param max_cur_mom is the maximum of current order Cartesian multipole moments
  !> \param cur_order_geo is the current order of geometric derivatives on Dirac delta function
  !> \param delta_wrt_diporg contains the relative coordinates of Dirac delta function
  !>        w.r.t. dipole origin
  !> \param dim_hgto_bra is the dimension of HGTOs on bra center
  !> \param dim_hgto_ket is the dimension of HGTOs on ket center
  !> \param num_low_geo is the number of lower order geometric derivatives on Dirac delta function
  !> \param low_zero_pints contains integrals of lower order geometric derivatives on
  !>        Dirac delta function and zeroth order Cartesian multipole moments
  !> \param dim_cur_mom is the dimension of current order Cartesian multipole moments
  !> \param low_geo_pints contains integrals of lower order geometric derivatives on
  !>        Dirac delta function and lower order Cartesian multipole moments
  !> \param num_cur_geo is the number of current order geometric derivatives on Dirac delta function
  !> \param cur_zero_pints contains integrals of current order geometric derivatives on
  !>        Dirac delta function and zeroth order Cartesian multipole moments
  !> \param dim_up_mom is the dimension of upper order Cartesian multipole moments
  !> \return cur_geo_pints contains the integrals of current order geometric derivatives
  !>         on Dirac delta function and orders of Cartesian multipole moments from 1 to
  !>         \var(max_cur_mom)+1
  subroutine sub_delta_moment(max_cur_mom, cur_order_geo, delta_wrt_diporg, &
                              dim_hgto_bra, dim_hgto_ket, num_low_geo,      &
                              low_zero_pints, dim_cur_mom, low_geo_pints,   &
                              num_cur_geo, cur_zero_pints, dim_up_mom,      &
                              cur_geo_pints)
    use xkind
    implicit none
    integer, intent(in) :: max_cur_mom
    integer, intent(in) :: cur_order_geo
    real(REALK), intent(in) :: delta_wrt_diporg(3)
    integer, intent(in) :: dim_hgto_bra
    integer, intent(in) :: dim_hgto_ket
    integer, intent(in) :: num_low_geo
    real(REALK), intent(in) :: low_zero_pints(dim_hgto_bra,dim_hgto_ket,num_low_geo)
    integer, intent(in) :: dim_cur_mom
    real(REALK), intent(in) :: low_geo_pints(dim_hgto_bra,dim_hgto_ket, &
                                             dim_cur_mom,num_low_geo)
    integer, intent(in) :: num_cur_geo
    real(REALK), intent(in) :: cur_zero_pints(dim_hgto_bra,dim_hgto_ket,num_cur_geo)
    integer, intent(in) :: dim_up_mom
    real(REALK), intent(out) :: cur_geo_pints(dim_hgto_bra,dim_hgto_ket, &
                                              dim_up_mom,num_cur_geo)
!f2py intent(in) :: max_cur_mom
!f2py intent(in) :: cur_order_geo
!f2py intent(in) :: delta_wrt_diporg
!f2py intent(hide) :: dim_hgto_bra
!f2py intent(hide) :: dim_hgto_ket
!f2py intent(hide) :: num_low_geo
!f2py intent(in) :: low_zero_pints
!f2py intent(hide) :: dim_cur_mom
!f2py intent(in) :: low_geo_pints
!f2py depend(dim_hgto_bra) :: low_geo_pints
!f2py depend(dim_hgto_ket) :: low_geo_pints
!f2py depend(num_low_geo) :: low_geo_pints
!f2py intent(hide) :: num_cur_geo
!f2py intent(in) :: cur_zero_pints
!f2py depend(dim_hgto_bra) :: cur_zero_pints
!f2py depend(dim_hgto_ket) :: cur_zero_pints
!f2py intent(in) :: dim_up_mom
!f2py intent(out) :: cur_geo_pints
!f2py depend(dim_hgto_bra) :: cur_geo_pints
!f2py depend(dim_hgto_ket) :: cur_geo_pints
!f2py depend(dim_up_mom) :: cur_geo_pints
!f2py depend(num_cur_geo) :: cur_geo_pints
    integer addr_low_xgeo  !address of lower order geometric derivatives along x-direction
    integer addr_low_ygeo  !address of lower order geometric derivatives along y-direction
    integer addr_low_zgeo  !address of lower order geometric derivatives along z-direction
    integer addr_cur_geo   !address of current order geometric derivatives
    integer igeo, jgeo     !incremental recorders of address of geometric derivatives
    integer order_geo_xy   !order of geometric derivatives along x- and y-direction
    integer order_geo_x    !order of geometric derivatives along x-direction
    integer base_cur_mom   !base address of current order Cartesian multipole moments
    integer addr_cur_mom   !address of current order Cartesian multipole moments
    integer addr_up_mom    !address of upper order Cartesian multipole moments
    integer order_mom      !incremental recorder over the orders of Cartesian multipole moments
    integer imom, jmom     !incremental recorders over xyz components of Cartesian multipole moments
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time   
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    select case(max_cur_mom)
    ! only the first order Cartesian multipole moments required
    case(0)
#if defined(DEBUG)
      write(STDOUT,100) "1st order Cartesian multipole moments return"
      write(STDOUT,100) "order of geometric derivatives", cur_order_geo
#endif
      ! (1) x...x component of geometric derivatives
      cur_geo_pints(:,:,1,1) = delta_wrt_diporg(1)*cur_zero_pints(:,:,1) &      !px
                             + real(cur_order_geo,REALK)*low_zero_pints(:,:,1)
      cur_geo_pints(:,:,2,1) = delta_wrt_diporg(2)*cur_zero_pints(:,:,1)        !py
      cur_geo_pints(:,:,3,1) = delta_wrt_diporg(3)*cur_zero_pints(:,:,1)        !pz
      ! (2) x...xy to xy...y components of geometric derivatives
      addr_low_xgeo = 1
      addr_low_ygeo = 0
      addr_cur_geo = 1
      do igeo = 1, cur_order_geo-1
        addr_low_xgeo = addr_low_xgeo+1
        addr_low_ygeo = addr_low_ygeo+1
        addr_cur_geo = addr_cur_geo+1
        cur_geo_pints(:,:,1,addr_cur_geo)                        &  !px
          = delta_wrt_diporg(1)*cur_zero_pints(:,:,addr_cur_geo) &
          + real(cur_order_geo-igeo,REALK)*low_zero_pints(:,:,addr_low_xgeo)
        cur_geo_pints(:,:,2,addr_cur_geo)                        &  !py
          = delta_wrt_diporg(2)*cur_zero_pints(:,:,addr_cur_geo) &
          + real(igeo,REALK)*low_zero_pints(:,:,addr_low_ygeo)
        cur_geo_pints(:,:,3,addr_cur_geo) &                         !pz
          = delta_wrt_diporg(3)*cur_zero_pints(:,:,addr_cur_geo)
      end do
      ! (3) y...y component of geometric derivatives
      addr_low_ygeo = addr_low_ygeo+1
      addr_cur_geo = addr_cur_geo+1
      cur_geo_pints(:,:,1,addr_cur_geo) &                         !px
        = delta_wrt_diporg(1)*cur_zero_pints(:,:,addr_cur_geo)
      cur_geo_pints(:,:,2,addr_cur_geo)                        &  !py
        = delta_wrt_diporg(2)*cur_zero_pints(:,:,addr_cur_geo) &  
        + real(cur_order_geo,REALK)*low_zero_pints(:,:,addr_low_ygeo)
      cur_geo_pints(:,:,3,addr_cur_geo) &                         !pz
        = delta_wrt_diporg(3)*cur_zero_pints(:,:,addr_cur_geo)    
      ! (4) x...xz...z to y...yz...z components of geometric derivatives
      addr_low_zgeo = 0
      do igeo = 1, cur_order_geo-1
        ! (4.1) x...xz...z component of geometric derivatives
        addr_low_xgeo = addr_low_xgeo+1
        addr_low_zgeo = addr_low_zgeo+1
        addr_cur_geo = addr_cur_geo+1
        order_geo_xy = cur_order_geo-igeo
        cur_geo_pints(:,:,1,addr_cur_geo)                        &  !px
          = delta_wrt_diporg(1)*cur_zero_pints(:,:,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_zero_pints(:,:,addr_low_xgeo)
        cur_geo_pints(:,:,2,addr_cur_geo) &                         !py
          = delta_wrt_diporg(2)*cur_zero_pints(:,:,addr_cur_geo)
        cur_geo_pints(:,:,3,addr_cur_geo)                        &  !pz
          = delta_wrt_diporg(3)*cur_zero_pints(:,:,addr_cur_geo) &
          + real(igeo,REALK)*low_zero_pints(:,:,addr_low_zgeo)
        ! (4.2) x...xyz...z to xy...yz...z components of geometric derivatives
        do jgeo = 1, cur_order_geo-(igeo+1)
          addr_low_xgeo = addr_low_xgeo+1
          addr_low_ygeo = addr_low_ygeo+1
          addr_low_zgeo = addr_low_zgeo+1
          addr_cur_geo = addr_cur_geo+1
          cur_geo_pints(:,:,1,addr_cur_geo)                        &  !px
            = delta_wrt_diporg(1)*cur_zero_pints(:,:,addr_cur_geo) &
            + real(order_geo_xy-jgeo,REALK)*low_zero_pints(:,:,addr_low_xgeo)
          cur_geo_pints(:,:,2,addr_cur_geo)                        &  !py
            = delta_wrt_diporg(2)*cur_zero_pints(:,:,addr_cur_geo) &
            + real(jgeo,REALK)*low_zero_pints(:,:,addr_low_ygeo)
          cur_geo_pints(:,:,3,addr_cur_geo)                        &  !pz
            = delta_wrt_diporg(3)*cur_zero_pints(:,:,addr_cur_geo) &
            + real(igeo,REALK)*low_zero_pints(:,:,addr_low_zgeo)
        end do
        ! (4.3) y...yz...z component of geometric derivatives
        addr_low_ygeo = addr_low_ygeo+1
        addr_low_zgeo = addr_low_zgeo+1
        addr_cur_geo = addr_cur_geo+1
        cur_geo_pints(:,:,1,addr_cur_geo) &                         !px
          = delta_wrt_diporg(1)*cur_zero_pints(:,:,addr_cur_geo)
        cur_geo_pints(:,:,2,addr_cur_geo)                        &  !py
          = delta_wrt_diporg(2)*cur_zero_pints(:,:,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_zero_pints(:,:,addr_low_ygeo)
        cur_geo_pints(:,:,3,addr_cur_geo)                        &  !pz
          = delta_wrt_diporg(3)*cur_zero_pints(:,:,addr_cur_geo) &
          + real(igeo,REALK)*low_zero_pints(:,:,addr_low_zgeo)
      end do
      ! (5) z...z component of geometric derivatives
      addr_low_zgeo = addr_low_zgeo+1
      addr_cur_geo = addr_cur_geo+1
      cur_geo_pints(:,:,1,addr_cur_geo) &                         !px
        = delta_wrt_diporg(1)*cur_zero_pints(:,:,addr_cur_geo)
      cur_geo_pints(:,:,2,addr_cur_geo) &                         !py
        = delta_wrt_diporg(2)*cur_zero_pints(:,:,addr_cur_geo)
      cur_geo_pints(:,:,3,addr_cur_geo)                        &  !pz
        = delta_wrt_diporg(3)*cur_zero_pints(:,:,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_zero_pints(:,:,addr_low_zgeo)
    ! the first and second order Cartesian multipole moments
    case(1)
#if defined(DEBUG)
      write(STDOUT,100) "1st and 2nd order Cartesian multipole moments return"
      write(STDOUT,100) "order of geometric derivatives", cur_order_geo
#endif
      ! (1) x...x component of geometric derivatives
      cur_geo_pints(:,:,1,1) = delta_wrt_diporg(1)*cur_zero_pints(:,:,1) &       !px
                             + real(cur_order_geo,REALK)*low_zero_pints(:,:,1)
      cur_geo_pints(:,:,2,1) = delta_wrt_diporg(2)*cur_zero_pints(:,:,1)         !py
      cur_geo_pints(:,:,3,1) = delta_wrt_diporg(3)*cur_zero_pints(:,:,1)         !pz
      cur_geo_pints(:,:,4,1) = delta_wrt_diporg(1)*cur_geo_pints(:,:,1,1) &      !dxx
                             + real(cur_order_geo,REALK)*low_geo_pints(:,:,1,1)
      cur_geo_pints(:,:,5,1) = delta_wrt_diporg(2)*cur_geo_pints(:,:,1,1)        !dxy
      cur_geo_pints(:,:,6,1) = delta_wrt_diporg(2)*cur_geo_pints(:,:,2,1)        !dyy
      cur_geo_pints(:,:,7,1) = delta_wrt_diporg(3)*cur_geo_pints(:,:,1,1)        !dxz
      cur_geo_pints(:,:,8,1) = delta_wrt_diporg(3)*cur_geo_pints(:,:,2,1)        !dyz
      cur_geo_pints(:,:,9,1) = delta_wrt_diporg(3)*cur_geo_pints(:,:,3,1)        !dzz
      ! (2) x...xy to xy...y components of geometric derivatives
      addr_low_xgeo = 1
      addr_low_ygeo = 0
      addr_cur_geo = 1
      do igeo = 1, cur_order_geo-1
        addr_low_xgeo = addr_low_xgeo+1
        addr_low_ygeo = addr_low_ygeo+1
        addr_cur_geo = addr_cur_geo+1
        order_geo_x = cur_order_geo-igeo
        cur_geo_pints(:,:,1,addr_cur_geo)                        &  !px
          = delta_wrt_diporg(1)*cur_zero_pints(:,:,addr_cur_geo) &
          + real(order_geo_x,REALK)*low_zero_pints(:,:,addr_low_xgeo)
        cur_geo_pints(:,:,2,addr_cur_geo)                        &  !py
          = delta_wrt_diporg(2)*cur_zero_pints(:,:,addr_cur_geo) &
          + real(igeo,REALK)*low_zero_pints(:,:,addr_low_ygeo)
        cur_geo_pints(:,:,3,addr_cur_geo) &                         !pz
          = delta_wrt_diporg(3)*cur_zero_pints(:,:,addr_cur_geo)
        cur_geo_pints(:,:,4,addr_cur_geo)                         &  !dxx
          = delta_wrt_diporg(1)*cur_geo_pints(:,:,1,addr_cur_geo) &
          + real(order_geo_x,REALK)*low_geo_pints(:,:,1,addr_low_xgeo)
        cur_geo_pints(:,:,5,addr_cur_geo)                         &  !dxy
          = delta_wrt_diporg(2)*cur_geo_pints(:,:,1,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,:,1,addr_low_ygeo)
        cur_geo_pints(:,:,6,addr_cur_geo)                         &  !dyy
          = delta_wrt_diporg(2)*cur_geo_pints(:,:,2,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,:,2,addr_low_ygeo)
        cur_geo_pints(:,:,7,addr_cur_geo) &                          !dxz
          = delta_wrt_diporg(3)*cur_geo_pints(:,:,1,addr_cur_geo)
        cur_geo_pints(:,:,8,addr_cur_geo) &                          !dyz
          = delta_wrt_diporg(3)*cur_geo_pints(:,:,2,addr_cur_geo)
        cur_geo_pints(:,:,9,addr_cur_geo) &                          !dzz
          = delta_wrt_diporg(3)*cur_geo_pints(:,:,3,addr_cur_geo)
      end do
      ! (3) y...y component of geometric derivatives
      addr_low_ygeo = addr_low_ygeo+1
      addr_cur_geo = addr_cur_geo+1
      cur_geo_pints(:,:,1,addr_cur_geo) &                         !px
        = delta_wrt_diporg(1)*cur_zero_pints(:,:,addr_cur_geo)
      cur_geo_pints(:,:,2,addr_cur_geo)                        &  !py
        = delta_wrt_diporg(2)*cur_zero_pints(:,:,addr_cur_geo) &  
        + real(cur_order_geo,REALK)*low_zero_pints(:,:,addr_low_ygeo)
      cur_geo_pints(:,:,3,addr_cur_geo) &                         !pz
        = delta_wrt_diporg(3)*cur_zero_pints(:,:,addr_cur_geo)
      cur_geo_pints(:,:,4,addr_cur_geo) &                          !dxx
        = delta_wrt_diporg(1)*cur_geo_pints(:,:,1,addr_cur_geo)
      cur_geo_pints(:,:,5,addr_cur_geo)                         &  !dxy
        = delta_wrt_diporg(2)*cur_geo_pints(:,:,1,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_geo_pints(:,:,1,addr_low_ygeo)
      cur_geo_pints(:,:,6,addr_cur_geo)                         &  !dyy
        = delta_wrt_diporg(2)*cur_geo_pints(:,:,2,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_geo_pints(:,:,2,addr_low_ygeo)
      cur_geo_pints(:,:,7,addr_cur_geo) &                          !dxz
        = delta_wrt_diporg(3)*cur_geo_pints(:,:,1,addr_cur_geo)
      cur_geo_pints(:,:,8,addr_cur_geo) &                          !dyz
        = delta_wrt_diporg(3)*cur_geo_pints(:,:,2,addr_cur_geo)
      cur_geo_pints(:,:,9,addr_cur_geo) &                          !dzz
        = delta_wrt_diporg(3)*cur_geo_pints(:,:,3,addr_cur_geo)
      ! (4) x...xz...z to y...yz...z components of geometric derivatives
      addr_low_zgeo = 0
      do igeo = 1, cur_order_geo-1
        ! (4.1) x...xz...z component of geometric derivatives
        addr_low_xgeo = addr_low_xgeo+1
        addr_low_zgeo = addr_low_zgeo+1
        addr_cur_geo = addr_cur_geo+1
        order_geo_xy = cur_order_geo-igeo
        cur_geo_pints(:,:,1,addr_cur_geo)                        &  !px
          = delta_wrt_diporg(1)*cur_zero_pints(:,:,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_zero_pints(:,:,addr_low_xgeo)
        cur_geo_pints(:,:,2,addr_cur_geo) &                         !py
          = delta_wrt_diporg(2)*cur_zero_pints(:,:,addr_cur_geo)
        cur_geo_pints(:,:,3,addr_cur_geo)                        &  !pz
          = delta_wrt_diporg(3)*cur_zero_pints(:,:,addr_cur_geo) &
          + real(igeo,REALK)*low_zero_pints(:,:,addr_low_zgeo)
        cur_geo_pints(:,:,4,addr_cur_geo)                         &  !dxx
          = delta_wrt_diporg(1)*cur_geo_pints(:,:,1,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_geo_pints(:,:,1,addr_low_xgeo)
        cur_geo_pints(:,:,5,addr_cur_geo) &                          !dxy
          = delta_wrt_diporg(2)*cur_geo_pints(:,:,1,addr_cur_geo)
        cur_geo_pints(:,:,6,addr_cur_geo) &                          !dyy
          = delta_wrt_diporg(2)*cur_geo_pints(:,:,2,addr_cur_geo)
        cur_geo_pints(:,:,7,addr_cur_geo)                         &  !dxz
          = delta_wrt_diporg(3)*cur_geo_pints(:,:,1,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,:,1,addr_low_zgeo)
        cur_geo_pints(:,:,8,addr_cur_geo)                         &  !dyz
          = delta_wrt_diporg(3)*cur_geo_pints(:,:,2,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,:,2,addr_low_zgeo)
        cur_geo_pints(:,:,9,addr_cur_geo)                         &  !dzz
          = delta_wrt_diporg(3)*cur_geo_pints(:,:,3,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,:,3,addr_low_zgeo)
        ! (4.2) x...xyz...z to xy...yz...z components of geometric derivatives
        do jgeo = 1, cur_order_geo-(igeo+1)
          addr_low_xgeo = addr_low_xgeo+1
          addr_low_ygeo = addr_low_ygeo+1
          addr_low_zgeo = addr_low_zgeo+1
          addr_cur_geo = addr_cur_geo+1
          cur_geo_pints(:,:,1,addr_cur_geo)                        &  !px
            = delta_wrt_diporg(1)*cur_zero_pints(:,:,addr_cur_geo) &
            + real(order_geo_xy-jgeo,REALK)*low_zero_pints(:,:,addr_low_xgeo)
          cur_geo_pints(:,:,2,addr_cur_geo)                        &  !py
            = delta_wrt_diporg(2)*cur_zero_pints(:,:,addr_cur_geo) &
            + real(jgeo,REALK)*low_zero_pints(:,:,addr_low_ygeo)
          cur_geo_pints(:,:,3,addr_cur_geo)                        &  !pz
            = delta_wrt_diporg(3)*cur_zero_pints(:,:,addr_cur_geo) &
            + real(igeo,REALK)*low_zero_pints(:,:,addr_low_zgeo)
          cur_geo_pints(:,:,4,addr_cur_geo)                         &  !dxx
            = delta_wrt_diporg(1)*cur_geo_pints(:,:,1,addr_cur_geo) &
            + real(order_geo_xy-jgeo,REALK)*low_geo_pints(:,:,1,addr_low_xgeo)
          cur_geo_pints(:,:,5,addr_cur_geo)                         &  !dxy
            = delta_wrt_diporg(2)*cur_geo_pints(:,:,1,addr_cur_geo) &
            + real(jgeo,REALK)*low_geo_pints(:,:,1,addr_low_ygeo)
          cur_geo_pints(:,:,6,addr_cur_geo)                         &  !dyy
            = delta_wrt_diporg(2)*cur_geo_pints(:,:,2,addr_cur_geo) &
            + real(jgeo,REALK)*low_geo_pints(:,:,2,addr_low_ygeo)
          cur_geo_pints(:,:,7,addr_cur_geo)                         &  !dxz
            = delta_wrt_diporg(3)*cur_geo_pints(:,:,1,addr_cur_geo) &
            + real(igeo,REALK)*low_geo_pints(:,:,1,addr_low_zgeo)
          cur_geo_pints(:,:,8,addr_cur_geo)                         &  !dyz
            = delta_wrt_diporg(3)*cur_geo_pints(:,:,2,addr_cur_geo) &
            + real(igeo,REALK)*low_geo_pints(:,:,2,addr_low_zgeo)
          cur_geo_pints(:,:,9,addr_cur_geo)                         &  !dzz
            = delta_wrt_diporg(3)*cur_geo_pints(:,:,3,addr_cur_geo) &
            + real(igeo,REALK)*low_geo_pints(:,:,3,addr_low_zgeo)
        end do
        ! (4.3) y...yz...z component of geometric derivatives
        addr_low_ygeo = addr_low_ygeo+1
        addr_low_zgeo = addr_low_zgeo+1
        addr_cur_geo = addr_cur_geo+1
        cur_geo_pints(:,:,1,addr_cur_geo) &                         !px
          = delta_wrt_diporg(1)*cur_zero_pints(:,:,addr_cur_geo)
        cur_geo_pints(:,:,2,addr_cur_geo)                        &  !py
          = delta_wrt_diporg(2)*cur_zero_pints(:,:,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_zero_pints(:,:,addr_low_ygeo)
        cur_geo_pints(:,:,3,addr_cur_geo)                        &  !pz
          = delta_wrt_diporg(3)*cur_zero_pints(:,:,addr_cur_geo) &
          + real(igeo,REALK)*low_zero_pints(:,:,addr_low_zgeo)
        cur_geo_pints(:,:,4,addr_cur_geo) &                          !dxx
          = delta_wrt_diporg(1)*cur_geo_pints(:,:,1,addr_cur_geo)
        cur_geo_pints(:,:,5,addr_cur_geo)                         &  !dxy
          = delta_wrt_diporg(2)*cur_geo_pints(:,:,1,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_geo_pints(:,:,1,addr_low_ygeo)
        cur_geo_pints(:,:,6,addr_cur_geo)                         &  !dyy
          = delta_wrt_diporg(2)*cur_geo_pints(:,:,2,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_geo_pints(:,:,2,addr_low_ygeo)
        cur_geo_pints(:,:,7,addr_cur_geo)                         &  !dxz
          = delta_wrt_diporg(3)*cur_geo_pints(:,:,1,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,:,1,addr_low_zgeo)
        cur_geo_pints(:,:,8,addr_cur_geo)                         &  !dyz
          = delta_wrt_diporg(3)*cur_geo_pints(:,:,2,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,:,2,addr_low_zgeo)
        cur_geo_pints(:,:,9,addr_cur_geo)                         &  !dzz
          = delta_wrt_diporg(3)*cur_geo_pints(:,:,3,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,:,3,addr_low_zgeo)
      end do
      ! (5) z...z component of geometric derivatives
      addr_low_zgeo = addr_low_zgeo+1
      addr_cur_geo = addr_cur_geo+1
      cur_geo_pints(:,:,1,addr_cur_geo) &                         !px
        = delta_wrt_diporg(1)*cur_zero_pints(:,:,addr_cur_geo)
      cur_geo_pints(:,:,2,addr_cur_geo) &                         !py
        = delta_wrt_diporg(2)*cur_zero_pints(:,:,addr_cur_geo)
      cur_geo_pints(:,:,3,addr_cur_geo)                        &  !pz
        = delta_wrt_diporg(3)*cur_zero_pints(:,:,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_zero_pints(:,:,addr_low_zgeo)
      cur_geo_pints(:,:,4,addr_cur_geo) &                          !dxx
        = delta_wrt_diporg(1)*cur_geo_pints(:,:,1,addr_cur_geo)
      cur_geo_pints(:,:,5,addr_cur_geo) &                          !dxy
        = delta_wrt_diporg(2)*cur_geo_pints(:,:,1,addr_cur_geo)
      cur_geo_pints(:,:,6,addr_cur_geo) &                          !dyy
        = delta_wrt_diporg(2)*cur_geo_pints(:,:,2,addr_cur_geo)
      cur_geo_pints(:,:,7,addr_cur_geo)                         &  !dxz
        = delta_wrt_diporg(3)*cur_geo_pints(:,:,1,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_geo_pints(:,:,1,addr_low_zgeo)
      cur_geo_pints(:,:,8,addr_cur_geo)                         &  !dyz
        = delta_wrt_diporg(3)*cur_geo_pints(:,:,2,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_geo_pints(:,:,2,addr_low_zgeo)
      cur_geo_pints(:,:,9,addr_cur_geo)                         &  !dzz
        = delta_wrt_diporg(3)*cur_geo_pints(:,:,3,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_geo_pints(:,:,3,addr_low_zgeo)
    ! the maximum order of Cartesian multipole moments required is > 2
    case default
#if defined(DEBUG)
      write(STDOUT,100) "higher order Cartesian multipole moments return", max_cur_mom
      write(STDOUT,100) "order of geometric derivatives", cur_order_geo
#endif
      ! (1) x...x component of geometric derivatives
      cur_geo_pints(:,:,1,1) = delta_wrt_diporg(1)*cur_zero_pints(:,:,1) &   !px
                             + real(cur_order_geo,REALK)*low_zero_pints(:,:,1)
      cur_geo_pints(:,:,2,1) = delta_wrt_diporg(2)*cur_zero_pints(:,:,1)     !py
      cur_geo_pints(:,:,3,1) = delta_wrt_diporg(3)*cur_zero_pints(:,:,1)     !pz
      cur_geo_pints(:,:,4,1) = delta_wrt_diporg(1)*cur_geo_pints(:,:,1,1) &  !dxx
                             + real(cur_order_geo,REALK)*low_geo_pints(:,:,1,1)
      cur_geo_pints(:,:,5,1) = delta_wrt_diporg(2)*cur_geo_pints(:,:,1,1)    !dxy
      cur_geo_pints(:,:,6,1) = delta_wrt_diporg(2)*cur_geo_pints(:,:,2,1)    !dyy
      cur_geo_pints(:,:,7,1) = delta_wrt_diporg(3)*cur_geo_pints(:,:,1,1)    !dxz
      cur_geo_pints(:,:,8,1) = delta_wrt_diporg(3)*cur_geo_pints(:,:,2,1)    !dyz
      cur_geo_pints(:,:,9,1) = delta_wrt_diporg(3)*cur_geo_pints(:,:,3,1)    !dzz
      ! initializes the (base) addresses of Cartesian multipole moments
      base_cur_mom = 3
      addr_up_mom = 9
      ! other order (>2) Cartesian multipole moments
      do order_mom = 2, max_cur_mom
        ! recurrence relation along x-direction
        addr_up_mom = addr_up_mom+1
        addr_cur_mom = base_cur_mom+1
        cur_geo_pints(:,:,addr_up_mom,1)                          &
          = delta_wrt_diporg(1)*cur_geo_pints(:,:,addr_cur_mom,1) &
          + real(cur_order_geo,REALK)*low_geo_pints(:,:,addr_cur_mom,1)
        ! recurrence relation along y-direction
        addr_cur_mom = base_cur_mom
        do imom = 0, order_mom
          addr_cur_mom = addr_cur_mom+1
          addr_up_mom = addr_up_mom+1
          cur_geo_pints(:,:,addr_up_mom,1) &
            = delta_wrt_diporg(2)*cur_geo_pints(:,:,addr_cur_mom,1)
        end do
        ! recurrence relation along z-direction
        addr_cur_mom = base_cur_mom
        do imom = 0, order_mom
          do jmom = 0, order_mom-imom
            addr_cur_mom = addr_cur_mom+1
            addr_up_mom = addr_up_mom+1
            cur_geo_pints(:,:,addr_up_mom,1) &
              = delta_wrt_diporg(3)*cur_geo_pints(:,:,addr_cur_mom,1)
          end do
        end do
        ! updates the base addresses
        base_cur_mom = addr_cur_mom
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
        cur_geo_pints(:,:,1,addr_cur_geo)                        &  !px
          = delta_wrt_diporg(1)*cur_zero_pints(:,:,addr_cur_geo) &
          + real(order_geo_x,REALK)*low_zero_pints(:,:,addr_low_xgeo)
        cur_geo_pints(:,:,2,addr_cur_geo)                        &  !py
          = delta_wrt_diporg(2)*cur_zero_pints(:,:,addr_cur_geo) &
          + real(igeo,REALK)*low_zero_pints(:,:,addr_low_ygeo)
        cur_geo_pints(:,:,3,addr_cur_geo) &                         !pz
          = delta_wrt_diporg(3)*cur_zero_pints(:,:,addr_cur_geo)
        cur_geo_pints(:,:,4,addr_cur_geo)                         &  !dxx
          = delta_wrt_diporg(1)*cur_geo_pints(:,:,1,addr_cur_geo) &
          + real(order_geo_x,REALK)*low_geo_pints(:,:,1,addr_low_xgeo)
        cur_geo_pints(:,:,5,addr_cur_geo)                         &  !dxy
          = delta_wrt_diporg(2)*cur_geo_pints(:,:,1,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,:,1,addr_low_ygeo)
        cur_geo_pints(:,:,6,addr_cur_geo)                         &  !dyy
          = delta_wrt_diporg(2)*cur_geo_pints(:,:,2,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,:,2,addr_low_ygeo)
        cur_geo_pints(:,:,7,addr_cur_geo) &                          !dxz
          = delta_wrt_diporg(3)*cur_geo_pints(:,:,1,addr_cur_geo)
        cur_geo_pints(:,:,8,addr_cur_geo) &                          !dyz
          = delta_wrt_diporg(3)*cur_geo_pints(:,:,2,addr_cur_geo)
        cur_geo_pints(:,:,9,addr_cur_geo) &                          !dzz
          = delta_wrt_diporg(3)*cur_geo_pints(:,:,3,addr_cur_geo)
        ! initializes the (base) addresses of HGTOs on ket center
        base_cur_mom = 3
        addr_up_mom = 9
        ! other order (>2) Cartesian multipole moments
        do order_mom = 2, max_cur_mom
          ! recurrence relation along x-direction
          addr_up_mom = addr_up_mom+1
          addr_cur_mom = base_cur_mom+1
          cur_geo_pints(:,:,addr_up_mom,addr_cur_geo)                          &
            = delta_wrt_diporg(1)*cur_geo_pints(:,:,addr_cur_mom,addr_cur_geo) &
            + real(order_geo_x,REALK)*low_geo_pints(:,:,addr_cur_mom,addr_low_xgeo)
          ! recurrence relation along y-direction
          addr_cur_mom = base_cur_mom
          do imom = 0, order_mom
            addr_cur_mom = addr_cur_mom+1
            addr_up_mom = addr_up_mom+1
            cur_geo_pints(:,:,addr_up_mom,addr_cur_geo)                          &
              = delta_wrt_diporg(2)*cur_geo_pints(:,:,addr_cur_mom,addr_cur_geo) &
              + real(igeo,REALK)*low_geo_pints(:,:,addr_cur_mom,addr_low_ygeo)
          end do
          ! recurrence relation along z-direction
          addr_cur_mom = base_cur_mom
          do imom = 0, order_mom
            do jmom = 0, order_mom-imom
              addr_cur_mom = addr_cur_mom+1
              addr_up_mom = addr_up_mom+1
              cur_geo_pints(:,:,addr_up_mom,addr_cur_geo) &
                = delta_wrt_diporg(3)*cur_geo_pints(:,:,addr_cur_mom,addr_cur_geo)
            end do
          end do
          ! updates the base addresses
          base_cur_mom = addr_cur_mom
        end do
      end do
      ! (3) y...y component of geometric derivatives
      addr_low_ygeo = addr_low_ygeo+1
      addr_cur_geo = addr_cur_geo+1
      cur_geo_pints(:,:,1,addr_cur_geo) &                         !px
        = delta_wrt_diporg(1)*cur_zero_pints(:,:,addr_cur_geo)
      cur_geo_pints(:,:,2,addr_cur_geo)                        &  !py
        = delta_wrt_diporg(2)*cur_zero_pints(:,:,addr_cur_geo) &  
        + real(cur_order_geo,REALK)*low_zero_pints(:,:,addr_low_ygeo)
      cur_geo_pints(:,:,3,addr_cur_geo) &                         !pz
        = delta_wrt_diporg(3)*cur_zero_pints(:,:,addr_cur_geo)
      cur_geo_pints(:,:,4,addr_cur_geo) &                          !dxx
        = delta_wrt_diporg(1)*cur_geo_pints(:,:,1,addr_cur_geo)
      cur_geo_pints(:,:,5,addr_cur_geo)                         &  !dxy
        = delta_wrt_diporg(2)*cur_geo_pints(:,:,1,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_geo_pints(:,:,1,addr_low_ygeo)
      cur_geo_pints(:,:,6,addr_cur_geo)                         &  !dyy
        = delta_wrt_diporg(2)*cur_geo_pints(:,:,2,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_geo_pints(:,:,2,addr_low_ygeo)
      cur_geo_pints(:,:,7,addr_cur_geo) &                          !dxz
        = delta_wrt_diporg(3)*cur_geo_pints(:,:,1,addr_cur_geo)
      cur_geo_pints(:,:,8,addr_cur_geo) &                          !dyz
        = delta_wrt_diporg(3)*cur_geo_pints(:,:,2,addr_cur_geo)
      cur_geo_pints(:,:,9,addr_cur_geo) &                          !dzz
        = delta_wrt_diporg(3)*cur_geo_pints(:,:,3,addr_cur_geo)
      ! initializes the (base) addresses of Cartesian multipole moments
      base_cur_mom = 3
      addr_up_mom = 9
      ! other order (>2) Cartesian multipole moments
      do order_mom = 2, max_cur_mom
        ! recurrence relation along x-direction
        addr_up_mom = addr_up_mom+1
        addr_cur_mom = base_cur_mom
        cur_geo_pints(:,:,addr_up_mom,addr_cur_geo) &
          = delta_wrt_diporg(1)*cur_geo_pints(:,:,addr_cur_mom+1,addr_cur_geo)
        ! recurrence relation along y-direction
        do imom = 0, order_mom
          addr_cur_mom = addr_cur_mom+1
          addr_up_mom = addr_up_mom+1
          cur_geo_pints(:,:,addr_up_mom,addr_cur_geo)                          &
            = delta_wrt_diporg(2)*cur_geo_pints(:,:,addr_cur_mom,addr_cur_geo) &
            + real(cur_order_geo,REALK)*low_geo_pints(:,:,addr_cur_mom,addr_low_ygeo)
        end do
        ! recurrence relation along z-direction
        addr_cur_mom = base_cur_mom
        do imom = 0, order_mom
          do jmom = 0, order_mom-imom
            addr_cur_mom = addr_cur_mom+1
            addr_up_mom = addr_up_mom+1
            cur_geo_pints(:,:,addr_up_mom,addr_cur_geo) &
              = delta_wrt_diporg(3)*cur_geo_pints(:,:,addr_cur_mom,addr_cur_geo)
          end do
        end do
        ! updates the base addresses
        base_cur_mom = addr_cur_mom
      end do
      ! (4) x...xz...z to y...yz...z components of geometric derivatives
      addr_low_zgeo = 0
      do igeo = 1, cur_order_geo-1
        ! (4.1) x...xz...z component of geometric derivatives
        addr_low_xgeo = addr_low_xgeo+1
        addr_low_zgeo = addr_low_zgeo+1
        addr_cur_geo = addr_cur_geo+1
        order_geo_xy = cur_order_geo-igeo
        cur_geo_pints(:,:,1,addr_cur_geo)                        &  !px
          = delta_wrt_diporg(1)*cur_zero_pints(:,:,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_zero_pints(:,:,addr_low_xgeo)
        cur_geo_pints(:,:,2,addr_cur_geo) &                         !py
          = delta_wrt_diporg(2)*cur_zero_pints(:,:,addr_cur_geo)
        cur_geo_pints(:,:,3,addr_cur_geo)                        &  !pz
          = delta_wrt_diporg(3)*cur_zero_pints(:,:,addr_cur_geo) &
          + real(igeo,REALK)*low_zero_pints(:,:,addr_low_zgeo)
        cur_geo_pints(:,:,4,addr_cur_geo)                         &  !dxx
          = delta_wrt_diporg(1)*cur_geo_pints(:,:,1,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_geo_pints(:,:,1,addr_low_xgeo)
        cur_geo_pints(:,:,5,addr_cur_geo) &                          !dxy
          = delta_wrt_diporg(2)*cur_geo_pints(:,:,1,addr_cur_geo)
        cur_geo_pints(:,:,6,addr_cur_geo) &                          !dyy
          = delta_wrt_diporg(2)*cur_geo_pints(:,:,2,addr_cur_geo)
        cur_geo_pints(:,:,7,addr_cur_geo)                         &  !dxz
          = delta_wrt_diporg(3)*cur_geo_pints(:,:,1,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,:,1,addr_low_zgeo)
        cur_geo_pints(:,:,8,addr_cur_geo)                         &  !dyz
          = delta_wrt_diporg(3)*cur_geo_pints(:,:,2,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,:,2,addr_low_zgeo)
        cur_geo_pints(:,:,9,addr_cur_geo)                         &  !dzz
          = delta_wrt_diporg(3)*cur_geo_pints(:,:,3,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,:,3,addr_low_zgeo)
        ! initializes the (base) addresses of Cartesian multipole moments
        base_cur_mom = 3
        addr_up_mom = 9
        ! other order (>2) Cartesian multipole moments
        do order_mom = 2, max_cur_mom
          ! recurrence relation along x-direction
          addr_up_mom = addr_up_mom+1
          addr_cur_mom = base_cur_mom+1
          cur_geo_pints(:,:,addr_up_mom,addr_cur_geo)                          &
            = delta_wrt_diporg(1)*cur_geo_pints(:,:,addr_cur_mom,addr_cur_geo) &
            + real(order_geo_xy,REALK)*low_geo_pints(:,:,addr_cur_mom,addr_low_xgeo)
          ! recurrence relation along y-direction
          addr_cur_mom = base_cur_mom
          do imom = 0, order_mom
            addr_cur_mom = addr_cur_mom+1
            addr_up_mom = addr_up_mom+1
            cur_geo_pints(:,:,addr_up_mom,addr_cur_geo) &
              = delta_wrt_diporg(2)*cur_geo_pints(:,:,addr_cur_mom,addr_cur_geo)
          end do
          ! recurrence relation along z-direction
          addr_cur_mom = base_cur_mom
          do imom = 0, order_mom
            do jmom = 0, order_mom-imom
              addr_cur_mom = addr_cur_mom+1
              addr_up_mom = addr_up_mom+1
              cur_geo_pints(:,:,addr_up_mom,addr_cur_geo)                          &
                = delta_wrt_diporg(3)*cur_geo_pints(:,:,addr_cur_mom,addr_cur_geo) &
                + real(igeo,REALK)*low_geo_pints(:,:,addr_cur_mom,addr_low_zgeo)
            end do
          end do
          ! updates the base addresses
          base_cur_mom = addr_cur_mom
        end do
        ! (4.2) x...xyz...z to xy...yz...z components of geometric derivatives
        do jgeo = 1, cur_order_geo-(igeo+1)
          addr_low_xgeo = addr_low_xgeo+1
          addr_low_ygeo = addr_low_ygeo+1
          addr_low_zgeo = addr_low_zgeo+1
          addr_cur_geo = addr_cur_geo+1
          cur_geo_pints(:,:,1,addr_cur_geo)                        &  !px
            = delta_wrt_diporg(1)*cur_zero_pints(:,:,addr_cur_geo) &
            + real(order_geo_xy-jgeo,REALK)*low_zero_pints(:,:,addr_low_xgeo)
          cur_geo_pints(:,:,2,addr_cur_geo)                        &  !py
            = delta_wrt_diporg(2)*cur_zero_pints(:,:,addr_cur_geo) &
            + real(jgeo,REALK)*low_zero_pints(:,:,addr_low_ygeo)
          cur_geo_pints(:,:,3,addr_cur_geo)                        &  !pz
            = delta_wrt_diporg(3)*cur_zero_pints(:,:,addr_cur_geo) &
            + real(igeo,REALK)*low_zero_pints(:,:,addr_low_zgeo)
          cur_geo_pints(:,:,4,addr_cur_geo)                         &  !dxx
            = delta_wrt_diporg(1)*cur_geo_pints(:,:,1,addr_cur_geo) &
            + real(order_geo_xy-jgeo,REALK)*low_geo_pints(:,:,1,addr_low_xgeo)
          cur_geo_pints(:,:,5,addr_cur_geo)                         &  !dxy
            = delta_wrt_diporg(2)*cur_geo_pints(:,:,1,addr_cur_geo) &
            + real(jgeo,REALK)*low_geo_pints(:,:,1,addr_low_ygeo)
          cur_geo_pints(:,:,6,addr_cur_geo)                         &  !dyy
            = delta_wrt_diporg(2)*cur_geo_pints(:,:,2,addr_cur_geo) &
            + real(jgeo,REALK)*low_geo_pints(:,:,2,addr_low_ygeo)
          cur_geo_pints(:,:,7,addr_cur_geo)                         &  !dxz
            = delta_wrt_diporg(3)*cur_geo_pints(:,:,1,addr_cur_geo) &
            + real(igeo,REALK)*low_geo_pints(:,:,1,addr_low_zgeo)
          cur_geo_pints(:,:,8,addr_cur_geo)                         &  !dyz
            = delta_wrt_diporg(3)*cur_geo_pints(:,:,2,addr_cur_geo) &
            + real(igeo,REALK)*low_geo_pints(:,:,2,addr_low_zgeo)
          cur_geo_pints(:,:,9,addr_cur_geo)                         &  !dzz
            = delta_wrt_diporg(3)*cur_geo_pints(:,:,3,addr_cur_geo) &
            + real(igeo,REALK)*low_geo_pints(:,:,3,addr_low_zgeo)
          ! initializes the (base) addresses of Cartesian multipole moments
          base_cur_mom = 3
          addr_up_mom = 9
          ! other order (>2) Cartesian multipole moments
          do order_mom = 2, max_cur_mom
            ! recurrence relation along x-direction
            addr_up_mom = addr_up_mom+1
            addr_cur_mom = base_cur_mom+1
            cur_geo_pints(:,:,addr_up_mom,addr_cur_geo)                          &
              = delta_wrt_diporg(1)*cur_geo_pints(:,:,addr_cur_mom,addr_cur_geo) &
              + real(order_geo_xy-jgeo,REALK)                                    &
              * low_geo_pints(:,:,addr_cur_mom,addr_low_xgeo)
            ! recurrence relation along y-direction
            addr_cur_mom = base_cur_mom
            do imom = 0, order_mom
              addr_cur_mom = addr_cur_mom+1
              addr_up_mom = addr_up_mom+1
              cur_geo_pints(:,:,addr_up_mom,addr_cur_geo)                          &
                = delta_wrt_diporg(2)*cur_geo_pints(:,:,addr_cur_mom,addr_cur_geo) &
                + real(jgeo,REALK)*low_geo_pints(:,:,addr_cur_mom,addr_low_ygeo)
            end do
            ! recurrence relation along z-direction
            addr_cur_mom = base_cur_mom
            do imom = 0, order_mom
              do jmom = 0, order_mom-imom
                addr_cur_mom = addr_cur_mom+1
                addr_up_mom = addr_up_mom+1
                cur_geo_pints(:,:,addr_up_mom,addr_cur_geo)                          &
                  = delta_wrt_diporg(3)*cur_geo_pints(:,:,addr_cur_mom,addr_cur_geo) &
                  + real(igeo,REALK)*low_geo_pints(:,:,addr_cur_mom,addr_low_zgeo)
              end do
            end do
            ! updates the base addresses
            base_cur_mom = addr_cur_mom
          end do
        end do
        ! (4.3) y...yz...z component of geometric derivatives
        addr_low_ygeo = addr_low_ygeo+1
        addr_low_zgeo = addr_low_zgeo+1
        addr_cur_geo = addr_cur_geo+1
        cur_geo_pints(:,:,1,addr_cur_geo) &                         !px
          = delta_wrt_diporg(1)*cur_zero_pints(:,:,addr_cur_geo)
        cur_geo_pints(:,:,2,addr_cur_geo)                        &  !py
          = delta_wrt_diporg(2)*cur_zero_pints(:,:,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_zero_pints(:,:,addr_low_ygeo)
        cur_geo_pints(:,:,3,addr_cur_geo)                        &  !pz
          = delta_wrt_diporg(3)*cur_zero_pints(:,:,addr_cur_geo) &
          + real(igeo,REALK)*low_zero_pints(:,:,addr_low_zgeo)
        cur_geo_pints(:,:,4,addr_cur_geo) &                          !dxx
          = delta_wrt_diporg(1)*cur_geo_pints(:,:,1,addr_cur_geo)
        cur_geo_pints(:,:,5,addr_cur_geo)                         &  !dxy
          = delta_wrt_diporg(2)*cur_geo_pints(:,:,1,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_geo_pints(:,:,1,addr_low_ygeo)
        cur_geo_pints(:,:,6,addr_cur_geo)                         &  !dyy
          = delta_wrt_diporg(2)*cur_geo_pints(:,:,2,addr_cur_geo) &
          + real(order_geo_xy,REALK)*low_geo_pints(:,:,2,addr_low_ygeo)
        cur_geo_pints(:,:,7,addr_cur_geo)                         &  !dxz
          = delta_wrt_diporg(3)*cur_geo_pints(:,:,1,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,:,1,addr_low_zgeo)
        cur_geo_pints(:,:,8,addr_cur_geo)                         &  !dyz
          = delta_wrt_diporg(3)*cur_geo_pints(:,:,2,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,:,2,addr_low_zgeo)
        cur_geo_pints(:,:,9,addr_cur_geo)                         &  !dzz
          = delta_wrt_diporg(3)*cur_geo_pints(:,:,3,addr_cur_geo) &
          + real(igeo,REALK)*low_geo_pints(:,:,3,addr_low_zgeo)
        ! initializes the (base) addresses of Cartesian multipole moments
        base_cur_mom = 3
        addr_up_mom = 9
        ! other order (>2) Cartesian multipole moments
        do order_mom = 2, max_cur_mom
          ! recurrence relation along x-direction
          addr_up_mom = addr_up_mom+1
          addr_cur_mom = base_cur_mom
          cur_geo_pints(:,:,addr_up_mom,addr_cur_geo) &
            = delta_wrt_diporg(1)*cur_geo_pints(:,:,addr_cur_mom+1,addr_cur_geo)
          ! recurrence relation along y-direction
          do imom = 0, order_mom
            addr_cur_mom = addr_cur_mom+1
            addr_up_mom = addr_up_mom+1
            cur_geo_pints(:,:,addr_up_mom,addr_cur_geo)                          &
              = delta_wrt_diporg(2)*cur_geo_pints(:,:,addr_cur_mom,addr_cur_geo) &
              + real(order_geo_xy,REALK)*low_geo_pints(:,:,addr_cur_mom,addr_low_ygeo)
          end do
          ! recurrence relation along z-direction
          addr_cur_mom = base_cur_mom
          do imom = 0, order_mom
            do jmom = 0, order_mom-imom
              addr_cur_mom = addr_cur_mom+1
              addr_up_mom = addr_up_mom+1
              cur_geo_pints(:,:,addr_up_mom,addr_cur_geo)                          &
                = delta_wrt_diporg(3)*cur_geo_pints(:,:,addr_cur_mom,addr_cur_geo) &
                + real(igeo,REALK)*low_geo_pints(:,:,addr_cur_mom,addr_low_zgeo)
            end do
          end do
          ! updates the base addresses
          base_cur_mom = addr_cur_mom
        end do
      end do
      ! (5) z...z component of geometric derivatives
      addr_low_zgeo = addr_low_zgeo+1
      addr_cur_geo = addr_cur_geo+1
      cur_geo_pints(:,:,1,addr_cur_geo) &                         !px
        = delta_wrt_diporg(1)*cur_zero_pints(:,:,addr_cur_geo)
      cur_geo_pints(:,:,2,addr_cur_geo) &                         !py
        = delta_wrt_diporg(2)*cur_zero_pints(:,:,addr_cur_geo)
      cur_geo_pints(:,:,3,addr_cur_geo)                        &  !pz
        = delta_wrt_diporg(3)*cur_zero_pints(:,:,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_zero_pints(:,:,addr_low_zgeo)
      cur_geo_pints(:,:,4,addr_cur_geo) &                          !dxx
        = delta_wrt_diporg(1)*cur_geo_pints(:,:,1,addr_cur_geo)
      cur_geo_pints(:,:,5,addr_cur_geo) &                          !dxy
        = delta_wrt_diporg(2)*cur_geo_pints(:,:,1,addr_cur_geo)
      cur_geo_pints(:,:,6,addr_cur_geo) &                          !dyy
        = delta_wrt_diporg(2)*cur_geo_pints(:,:,2,addr_cur_geo)
      cur_geo_pints(:,:,7,addr_cur_geo)                         &  !dxz
        = delta_wrt_diporg(3)*cur_geo_pints(:,:,1,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_geo_pints(:,:,1,addr_low_zgeo)
      cur_geo_pints(:,:,8,addr_cur_geo)                         &  !dyz
        = delta_wrt_diporg(3)*cur_geo_pints(:,:,2,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_geo_pints(:,:,2,addr_low_zgeo)
      cur_geo_pints(:,:,9,addr_cur_geo)                         &  !dzz
        = delta_wrt_diporg(3)*cur_geo_pints(:,:,3,addr_cur_geo) &
        + real(cur_order_geo,REALK)*low_geo_pints(:,:,3,addr_low_zgeo)
      ! initializes the (base) addresses of Cartesian multipole moments
      base_cur_mom = 3
      addr_up_mom = 9
      ! other order (>2) Cartesian multipole moments
      do order_mom = 2, max_cur_mom
        ! recurrence relation along x-direction
        addr_up_mom = addr_up_mom+1
        addr_cur_mom = base_cur_mom
        cur_geo_pints(:,:,addr_up_mom,addr_cur_geo) &
          = delta_wrt_diporg(1)*cur_geo_pints(:,:,addr_cur_mom+1,addr_cur_geo)
        ! recurrence relation along y-direction
        do imom = 0, order_mom
          addr_cur_mom = addr_cur_mom+1
          addr_up_mom = addr_up_mom+1
          cur_geo_pints(:,:,addr_up_mom,addr_cur_geo) &
            = delta_wrt_diporg(2)*cur_geo_pints(:,:,addr_cur_mom,addr_cur_geo)
        end do
        ! recurrence relation along z-direction
        addr_cur_mom = base_cur_mom
        do imom = 0, order_mom
          do jmom = 0, order_mom-imom
            addr_cur_mom = addr_cur_mom+1
            addr_up_mom = addr_up_mom+1
            cur_geo_pints(:,:,addr_up_mom,addr_cur_geo)                          &
              = delta_wrt_diporg(3)*cur_geo_pints(:,:,addr_cur_mom,addr_cur_geo) &
              + real(cur_order_geo,REALK)*low_geo_pints(:,:,addr_cur_mom,addr_low_zgeo)
          end do
        end do
        ! updates the base addresses
        base_cur_mom = addr_cur_mom
      end do
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "sub_delta_moment", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("sub_delta_moment>> ",A,I6)
#endif
  end subroutine sub_delta_moment

  !> \brief assigns the Dirac delta function integrals
  !> \author Bin Gao
  !> \date 2012-03-16
  !> \param dim_hgto_bra is the dimension of HGTOs on bra center
  !> \param dim_hgto_ket is the dimension of HGTOs on ket center
  !> \param dim_up_mom is the dimension of upper order Cartesian multipole moments
  !> \param num_geo_pot is the number of geometric derivatives on potential origin,
  !>        equals to \f$(\var(order_geo_pot)+1)(\var(order_geo_pot)+2)/2\f$
  !> \param cur_geo_pints contains the integrals with upper order Cartesian multipole moments
  !>        and current order geometric derivatives on Dirac delta function
  !> \param num_mom is the number of xyz components of Cartesian multipole moment,
  !>        equals to \f$(\var(order_mom)+1)(\var(order_mom)+2)/2\f$
  !> \return hmom_pints contains the primitive Hermite integrals with specified orders of
  !>         Cartesian multipole moments and geometric derivatives on Dirac delta function
  subroutine delta_moment_assign(dim_hgto_bra, dim_hgto_ket, dim_up_mom, num_geo_pot, &
                                 cur_geo_pints, num_mom, hmom_pints)
    use xkind
    implicit none
    integer, intent(in) :: dim_hgto_bra
    integer, intent(in) :: dim_hgto_ket
    integer, intent(in) :: dim_up_mom
    integer, intent(in) :: num_geo_pot
    real(REALK), intent(in) :: cur_geo_pints(dim_hgto_bra,dim_hgto_ket, &
                                             dim_up_mom,num_geo_pot)
    integer, intent(in) :: num_mom
    real(REALK), intent(out) :: hmom_pints(dim_hgto_bra,dim_hgto_ket, &
                                           num_mom,num_geo_pot)
!f2py intent(hide) :: dim_hgto_bra
!f2py intent(hide) :: dim_hgto_ket
!f2py intent(hide) :: dim_up_mom
!f2py intent(hide) :: num_geo_pot
!f2py intent(in) :: cur_geo_pints
!f2py intent(in) :: num_mom
!f2py intent(out) :: hmom_pints
!f2py depend(dim_hgto_bra) :: hmom_pints
!f2py depend(dim_hgto_ket) :: hmom_pints
!f2py depend(num_mom) :: hmom_pints
!f2py depend(num_geo_pot) :: hmom_pints
    integer start_mom      !start address of Cartesian multipole moments in temporary integrals
    integer igeo           !incremental recorder over geometric derivatives
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time   
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    start_mom = dim_up_mom-num_mom+1
    do igeo = 1, num_geo_pot
      hmom_pints(:,:,:,igeo) = cur_geo_pints(:,:,start_mom:dim_up_mom,igeo)
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "delta_moment_assign", STDOUT)
#endif
    return
  end subroutine delta_moment_assign
