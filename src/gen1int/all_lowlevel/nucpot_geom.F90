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
!!  This file transfers Boys functions to geometric derivatives of nuclear
!!  potential origin with zeroth order HGTOs.
!!
!!  2012-03-05, Bin Gao
!!  * rewrites to use fixed size temporary integrals, the size is determined
!!    from \fn(dim_nucpot_geom)
!!
!!  2010-10-23, Bin Gao
!!  * first version

#include "stdout.h"

  !> \brief transfers Boys functions to geometric derivatives of nuclear
  !>        potential origin with zeroth order HGTOs
  !> \author Bin Gao
  !> \date 2010-10-23
  !> \param coord_bra contains the coordinates of bra center
  !> \param exponent_bra is the exponent of primitive HGTO of bra center
  !> \param coord_ket contains the coordinates of ket center
  !> \param exponent_ket is the exponent of primitive HGTO of ket center
  !> \param nucpot_origin contains the coordinates of nuclear potential origin
  !> \param scal_const is the scale constant for nuclear attraction potential
  !> \param orders_geo_pot contains the orders of geometric derivatives of
  !>        nuclear potential origin
  !> \param order_elec is the order of electronic derivatives
  !> \param dim_geo_pot is the dimension of integrals
  !> \return geo_pot_pints contains the primitive integrals with zeroth order HGTOs
  subroutine nucpot_geom(coord_bra, exponent_bra, coord_ket, exponent_ket,      &
                         nucpot_origin, scal_const, orders_geo_pot, order_elec, &
                         dim_geo_pot, geo_pot_pints)
    use xkind
    implicit none
    real(REALK), intent(in) :: coord_bra(3)
    real(REALK), intent(in) :: exponent_bra
    real(REALK), intent(in) :: coord_ket(3)
    real(REALK), intent(in) :: exponent_ket
    real(REALK), intent(in) :: nucpot_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: orders_geo_pot(2)
    integer, intent(in) :: order_elec
    integer, intent(in) :: dim_geo_pot
    real(REALK), intent(out) :: geo_pot_pints(dim_geo_pot)
!f2py intent(in) :: coord_bra
!f2py intent(in) :: exponent_bra
!f2py intent(in) :: coord_ket
!f2py intent(in) :: exponent_ket
!f2py intent(in) :: nucpot_origin
!f2py intent(in) :: scal_const
!f2py intent(in) :: orders_geo_pot
!f2py intent(in) :: order_elec
!f2py intent(in) :: dim_geo_pot
!f2py intent(out) :: geo_pot_pints
!f2py depend(dim_geo_pot) :: geo_pot_pints
#include "private/pi.h"
    real(REALK) total_expnt        !total exponent
    real(REALK) recip_total_expnt  !reciprocal of total exponent
    real(REALK) reduced_expnt      !reduced exponent
    real(REALK) nucorg_wrt_cc(3)   !relative coordinates of nuclear potential origin w.r.t. center-of-charge
    real(REALK) sd_bra_ket         !square of the relative distance between bra and ket
    integer min_order_boys         !minimum order of Boys functions
    logical odd_min_geo            !indicates if the minimum order of geometric derivatives is odd or even
    real(REALK) arg_boys           !argument of Boys functions
    real(REALK), allocatable :: aux_pints(:)
                                   !auxiliary integrals
    integer max_up_geo             !maximum upper order of geometric derivatives of potential origin
    integer min_up_geo             !minimum upper order of geometric derivatives of potential origin
    integer dim_cur_geo            !dimension of integrals with current order geometric derivatives
    integer dim_up_geo             !dimension of integrals with upper order geometric derivatives
    integer cur_boys               !pointer to current order Boys function in temporary integrals
    integer up_boys                !pointer of upper order Boys function in temporary integrals
    integer dim_tmp                !dimension of temporary integrals
    real(REALK), allocatable :: tmp_ints(:,:)
                                   !temporary integrals
    integer order_boys             !current order of Boys function
    integer ierr                   !error information
#if defined(XTIME)
    real(REALK) curr_time          !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! prepares the prerequisite quantities
    total_expnt = exponent_bra+exponent_ket
    recip_total_expnt = 1.0_REALK/total_expnt
    reduced_expnt = exponent_bra*exponent_ket*recip_total_expnt
    sd_bra_ket = 0.0_REALK
    arg_boys = 0.0_REALK
    do ierr = 1, 3
      nucorg_wrt_cc(ierr) = nucpot_origin(ierr)-recip_total_expnt &
                          * (exponent_bra*coord_bra(ierr)+exponent_ket*coord_ket(ierr))
      sd_bra_ket = sd_bra_ket+(coord_ket(ierr)-coord_bra(ierr))**2
      arg_boys = arg_boys+nucorg_wrt_cc(ierr)**2
    end do
    arg_boys = total_expnt*arg_boys
    ! computes the minimum order of Boys functions
    if (mod(orders_geo_pot(1),2)==1) then
      min_order_boys = (orders_geo_pot(1)+1)/2
      odd_min_geo = .true.
    else
      min_order_boys = orders_geo_pot(1)/2
      odd_min_geo = .false.
    end if
    ! allocates the memory for auxiliary integrals
    allocate(aux_pints(min_order_boys:orders_geo_pot(2)), stat=ierr)
    if (ierr/=0)                                                     &
      call error_stop("nucpot_geom", "failed to allocate aux_pints", &
                      orders_geo_pot(2)-min_order_boys+1)
    ! computes the Boys functions
    call aux_boys_vec(min_order_boys, orders_geo_pot(2), arg_boys, aux_pints)
    ! computes \f$\bar{C}(-2b_{j\lambda})^{|\boldsymbol{n}|}\frac{2\pi}{p_{ij}}%
    !             \mathrm{e}^{-u_{ij}R_{\kappa\lambda}^2}\f$
    arg_boys = scal_const*(-exponent_ket-exponent_ket)**order_elec &
             * exp(-reduced_expnt*sd_bra_ket)*(PI+PI)*recip_total_expnt
    select case(orders_geo_pot(2))
    ! returns s-shell integrals
    case(0)
      geo_pot_pints(1) = arg_boys*aux_pints(0)
      deallocate(aux_pints)
    ! returns integrals up to p-shell
    case(1)
      ! computes the first order auxiliary integral
      aux_pints(1) = -arg_boys*(total_expnt+total_expnt)*aux_pints(1)
      ! returns integrals of s- and p-shells
      if (orders_geo_pot(1)==0) then
        geo_pot_pints(1) = arg_boys*aux_pints(0)          !s
        geo_pot_pints(2) = nucorg_wrt_cc(1)*aux_pints(1)  !px
        geo_pot_pints(3) = nucorg_wrt_cc(2)*aux_pints(1)  !py
        geo_pot_pints(4) = nucorg_wrt_cc(3)*aux_pints(1)  !pz
      ! only returns integrals of p-shell
      else
        geo_pot_pints(1) = nucorg_wrt_cc(1)*aux_pints(1)  !px
        geo_pot_pints(2) = nucorg_wrt_cc(2)*aux_pints(1)  !py
        geo_pot_pints(3) = nucorg_wrt_cc(3)*aux_pints(1)  !pz
      end if
      deallocate(aux_pints)
    ! returns integrals whose maximum order is > 1
    case default
      ! computes \f$\bar{C}(-2b_{j\lambda})^{|\boldsymbol{n}|}\frac{2\pi}{p_{ij}}%
      !             \mathrm{e}^{-u_{ij}R_{\kappa\lambda}^2}%
      !             (-2p_{ij})^{|\boldsymbol{L}_{C}|_{\max}}\f$
      arg_boys = arg_boys*(-2.0_REALK*total_expnt)**orders_geo_pot(2)
      ! computes the last auxiliary integral
      aux_pints(orders_geo_pot(2)) = arg_boys*aux_pints(orders_geo_pot(2))
      ! \f$-\frac{1}{2p_{ij}}\f$
      recip_total_expnt = -0.5_REALK*recip_total_expnt
      ! allocates memory for temporary integrals
      call dim_nucpot_geom(odd_min_geo, orders_geo_pot(2), min_order_boys, dim_tmp)
      allocate(tmp_ints(dim_tmp,2), stat=ierr)
      if (ierr/=0) &
        call error_stop("nucpot_geom", "failed to allocate tmp_ints", dim_tmp*2)
      ! sets the temporary integrals of s- and p-shells
      arg_boys = recip_total_expnt*arg_boys                          !updates the scale factor
      tmp_ints(1,1) = arg_boys*aux_pints(orders_geo_pot(2)-1)        !s
      tmp_ints(2,1) = nucorg_wrt_cc(1)*aux_pints(orders_geo_pot(2))  !px
      tmp_ints(3,1) = nucorg_wrt_cc(2)*aux_pints(orders_geo_pot(2))  !py
      tmp_ints(4,1) = nucorg_wrt_cc(3)*aux_pints(orders_geo_pot(2))  !pz
      ! initializes the maximum upper order of geometric derivatives of potential origin
      max_up_geo = 1
      ! initializes the dimension of integrals with current order geometric derivatives
      dim_cur_geo = 4
      ! initializes the pointer of current order Boys function in temporary integrals
      cur_boys = 1
      up_boys = 2
      ! loops over the orders of Boys fucntions
      do order_boys = orders_geo_pot(2)-2, min_order_boys, -1
        ! updates the maximum upper order of geometric derivatives of potential origin
        max_up_geo = max_up_geo+1
        ! sets the dimension of integrals with upper order geometric derivatives
        dim_up_geo = dim_cur_geo+(max_up_geo+1)*(max_up_geo+2)/2
        ! switches the pointers
        cur_boys = 3-cur_boys
        up_boys = 3-up_boys
        ! updates the scale factor
        arg_boys = recip_total_expnt*arg_boys
        ! sets the temporary integrals of s- and p-shells
        tmp_ints(1,cur_boys) = arg_boys*aux_pints(order_boys)        !s
        tmp_ints(2,cur_boys) = nucorg_wrt_cc(1)*tmp_ints(1,up_boys)  !px
        tmp_ints(3,cur_boys) = nucorg_wrt_cc(2)*tmp_ints(1,up_boys)  !py
        tmp_ints(4,cur_boys) = nucorg_wrt_cc(3)*tmp_ints(1,up_boys)  !pz
        ! performs the recurrecen relation for other shells
        call sub_nucpot_geom_d(2, max_up_geo, nucorg_wrt_cc, dim_cur_geo,     &
                               tmp_ints(1:dim_cur_geo,up_boys), dim_up_geo-4, &
                               tmp_ints(5:dim_up_geo,cur_boys))
        ! updates the dimension of integrals with current order geometric derivatives
        dim_cur_geo = dim_up_geo
      end do
      ! cleans the auxiliary integrals
      deallocate(aux_pints)
      ! recurrence relations for geometric derivatives starting from p-shell,
      ! followed by f-, h-shells, and so on
      if (odd_min_geo) then
        ! updates the maximum upper order of geometric derivatives of potential origin
        max_up_geo = max_up_geo+1
        ! sets the dimension of integrals with upper order geometric derivatives,
        ! -1 due to removing s-shell
        dim_up_geo = dim_cur_geo+(max_up_geo+1)*(max_up_geo+2)/2-1
        ! switches the pointers
        cur_boys = 3-cur_boys
        up_boys = 3-up_boys
        ! sets the temporary integrals of s- and p-shells
        tmp_ints(1,cur_boys) = nucorg_wrt_cc(1)*tmp_ints(1,up_boys)  !px
        tmp_ints(2,cur_boys) = nucorg_wrt_cc(2)*tmp_ints(1,up_boys)  !py
        tmp_ints(3,cur_boys) = nucorg_wrt_cc(3)*tmp_ints(1,up_boys)  !pz
        ! performs the recurrecen relation for other shells
        call sub_nucpot_geom_d(2, max_up_geo, nucorg_wrt_cc, dim_cur_geo,     &
                               tmp_ints(1:dim_cur_geo,up_boys), dim_up_geo-3, &
                               tmp_ints(4:dim_up_geo,cur_boys))
        ! updates the dimension of integrals with current order geometric derivatives
        dim_cur_geo = dim_up_geo
        ! recurrence relations of other shells starting from f-shell
        min_up_geo = 1
      ! recurrence relations for geometric derivatives starting from d-shell,
      ! followed by g-, i-shells, and so on
      else
        min_up_geo = 0
      end if
      ! loops over the left orders of Boys functions till to the zeroth order
      do order_boys = min_order_boys-min_up_geo-1, 0, -1
        ! updates the upper orders of geometric derivatives of potential origin
        min_up_geo = min_up_geo+2
        max_up_geo = max_up_geo+1
        ! sets the dimension of integrals with upper order geometric derivatives,
        ! adds one maximum order shell, and removes two minimum order shells
        dim_up_geo = dim_cur_geo+(max_up_geo+1)*(max_up_geo+2)/2 &
                   - min_up_geo*min_up_geo
        ! switches the pointers
        cur_boys = 3-cur_boys
        up_boys = 3-up_boys
        ! performs the recurrence relations for different shells
        call sub_nucpot_geom_d(min_up_geo, max_up_geo, nucorg_wrt_cc,        &
                               dim_cur_geo, tmp_ints(1:dim_cur_geo,up_boys), &
                               dim_up_geo, tmp_ints(1:dim_up_geo,cur_boys))
        ! updates the dimension of integrals with current order geometric derivatives
        dim_cur_geo = dim_up_geo
      end do
      ! assigns the returned integrals
      geo_pot_pints = tmp_ints(1:dim_geo_pot,cur_boys)
      deallocate(tmp_ints)
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "nucpot_geom", STDOUT)
#endif
    return
  end subroutine nucpot_geom

  !> \brief gets the maximum dimension of temporary integrals used in recurrence relations
  !> \author Bin Gao
  !> \date 2012-03-05
  !> \param odd_min_geo indicates if the minimum order of geometric derivatives is odd or even
  !> \param max_order_boys is the maximum order of Boys functions
  !> \param min_order_boys is the minimum order of Boys functions
  !> \return dim_ints is the maximum dimension of temporary integrals
  subroutine dim_nucpot_geom(odd_min_geo, max_order_boys, min_order_boys, dim_ints)
    use xkind
    implicit none
    logical, intent(in) :: odd_min_geo
    integer, intent(in) :: max_order_boys
    integer, intent(in) :: min_order_boys
    integer, intent(out) :: dim_ints
!f2py intetn(in) :: odd_min_geo
!f2py intent(in) :: max_order_boys
!f2py intent(in) :: min_order_boys
!f2py intent(out) :: dim_ints
    integer max_geo_pot    !maximum order of geometric derivatives
    integer dim_geo_pot    !dimension of geometric derivatives
    integer order_boys     !incremental recorder over orders of Boys functions
    integer min_geo_pot    !minimum order of geometric derivatives
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! initializes the return value
    dim_ints = 0
    ! initializes the maximum order and dimension of geometric derivatives
    max_geo_pot = 1
    dim_geo_pot = 4
    ! loops over the orders of Boys fucntions
    do order_boys = max_order_boys-2, min_order_boys, -1
      ! updates the maximum order of geometric derivatives
      max_geo_pot = max_geo_pot+1
      ! sets the dimension of geometric derivatives
      dim_geo_pot = dim_geo_pot+(max_geo_pot+1)*(max_geo_pot+2)/2
      ! updates the maximum dimension
      if (dim_geo_pot>dim_ints) dim_ints = dim_geo_pot
    end do
    ! recurrence relations for geometric derivatives starting from p-shell,
    ! followed by f-, h-shells, and so on
    if (odd_min_geo) then
      ! updates the maximum order of geometric derivatives
      max_geo_pot = max_geo_pot+1
      ! sets the dimension of geometric derivatives, -1 due to removing s-shell
      dim_geo_pot = dim_geo_pot+(max_geo_pot+1)*(max_geo_pot+2)/2-1
      ! updates the maximum dimension
      if (dim_geo_pot>dim_ints) dim_ints = dim_geo_pot
      ! recurrence relations of other shells starting from f-shell
      min_geo_pot = 1
    ! recurrence relations for geometric derivatives starting from d-shell,
    ! followed by g-, i-shells, and so on
    else
      min_geo_pot = 0
    end if
    ! loops over the left orders of Boys functions till to the zeroth order
    do order_boys = min_order_boys-min_geo_pot-1, 0, -1
      ! updates the orders of geometric derivatives
      min_geo_pot = min_geo_pot+2
      max_geo_pot = max_geo_pot+1
      ! sets the dimension of geometric derivatives, adds one maximum order shell,
      ! and removes two minimum order shells
      dim_geo_pot = dim_geo_pot+(max_geo_pot+1)*(max_geo_pot+2)/2 &
                  - min_geo_pot*min_geo_pot
      ! updates the maximum dimension
      if (dim_geo_pot>dim_ints) dim_ints = dim_geo_pot
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "dim_nucpot_geom", STDOUT)
#endif
    return
  end subroutine dim_nucpot_geom

  !> \brief sub-recurrence relation by transferring Boys functions to geometric
  !>        derivatives (order\f$\ge2\f$) of nuclear potential origin with zeroth
  !>        order HGTOs
  !> \author Bin Gao
  !> \date 2010-10-23
  !> \param min_up_geo is minimum upper order of geometric derivatives of potential origin
  !> \param max_up_geo is maximum upper order of geometric derivatives of potential origin
  !> \param nucorg_wrt_cc contains the relative coordinates of nuclear potential
  !>        origin w.r.t. center-of-charge
  !> \param dim_cur_geo is the dimension of integrals of current order geometric derivatives
  !> \param cur_geo_pints contains the integrals of current order geometric derivatives
  !> \param dim_up_geo is the dimension of integrals of upper order geometric derivatives
  !> \return up_geo_pints contains the integrals of upper order geometric derivatives
  subroutine sub_nucpot_geom_d(min_up_geo, max_up_geo, nucorg_wrt_cc,  &
                               dim_cur_geo, cur_geo_pints, dim_up_geo, &
                               up_geo_pints)
    use xkind
    implicit none
    integer, intent(in) :: min_up_geo
    integer, intent(in) :: max_up_geo
    real(REALK), intent(in) :: nucorg_wrt_cc(3)
    integer, intent(in) :: dim_cur_geo
    real(REALK), intent(in) :: cur_geo_pints(dim_cur_geo)
    integer, intent(in) :: dim_up_geo
    real(REALK), intent(out) :: up_geo_pints(dim_up_geo)
!f2py intent(in) :: min_up_geo
!f2py intent(in) :: max_up_geo
!f2py intent(in) :: nucorg_wrt_cc
!f2py intent(hide) :: dim_cur_geo
!f2py intent(in) :: cur_geo_pints
!f2py intent(in) :: dim_up_geo
!f2py intent(out) :: up_geo_pints
!f2py depend(dim_up_geo) :: up_geo_pints
    integer min_cur_geo    !minimum current order of geometric derivatives
    integer addr_up_geo    !address of upper order geometric derivatives
    integer addr_cur_geo   !address of current order geometric derivatives
    integer addr_low_geo   !address of lower order geometric derivatives
    integer order_geo      !order of geometric derivatives
    integer igeo, jgeo     !incremental recorders
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time   
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    addr_up_geo = 0
    min_cur_geo = min_up_geo-1
    addr_cur_geo = min_cur_geo*min_up_geo/2
    addr_low_geo = 0
    ! loops over current orders of geometric derivatives of potential origin
    do order_geo = min_cur_geo, max_up_geo-1
      addr_up_geo = addr_up_geo+1
      addr_cur_geo = addr_cur_geo+1
#if defined(DEBUG)
      write(STDOUT,100) "order:", order_geo
      write(STDOUT,100) "x-direction"
      write(STDOUT,110) addr_up_geo, addr_cur_geo, addr_low_geo+1
      write(STDOUT,100) "------------------------------------"
#endif
      ! x-direction
      up_geo_pints(addr_up_geo)                        &
        = nucorg_wrt_cc(1)*cur_geo_pints(addr_cur_geo) &
        + real(order_geo,REALK)*cur_geo_pints(addr_low_geo+1)
      ! y-direction
      addr_up_geo = addr_up_geo+1
#if defined(DEBUG)
      write(STDOUT,100) "y-direction"
      write(STDOUT,120) addr_up_geo, addr_cur_geo
      write(STDOUT,100) "------------------------------------"
#endif
      up_geo_pints(addr_up_geo) = nucorg_wrt_cc(2)*cur_geo_pints(addr_cur_geo)
      do igeo = 1, order_geo
        addr_up_geo = addr_up_geo+1
#if defined(DEBUG)
        write(STDOUT,110) addr_up_geo, addr_cur_geo+igeo, addr_low_geo+igeo
        write(STDOUT,100) "------------------------------------"
#endif
        up_geo_pints(addr_up_geo)                             &
          = nucorg_wrt_cc(2)*cur_geo_pints(addr_cur_geo+igeo) &
          + real(igeo,REALK)*cur_geo_pints(addr_low_geo+igeo)
      end do
      ! z-direction
#if defined(DEBUG)
      write(STDOUT,100) "z-direction"
#endif
      addr_cur_geo = addr_cur_geo-1
      do jgeo = 0, order_geo
        addr_up_geo = addr_up_geo+1
        addr_cur_geo = addr_cur_geo+1
#if defined(DEBUG)
        write(STDOUT,120) addr_up_geo, addr_cur_geo
        write(STDOUT,100) "------------------------------------"
#endif
        up_geo_pints(addr_up_geo) = nucorg_wrt_cc(3)*cur_geo_pints(addr_cur_geo)
      end do
      do igeo = 1, order_geo
        do jgeo = 0, order_geo-igeo
          addr_up_geo = addr_up_geo+1
          addr_cur_geo = addr_cur_geo+1
          addr_low_geo = addr_low_geo+1
#if defined(DEBUG)
          write(STDOUT,110) addr_up_geo, addr_cur_geo, addr_low_geo
          write(STDOUT,100) "------------------------------------"
#endif
          up_geo_pints(addr_up_geo)                        &
            = nucorg_wrt_cc(3)*cur_geo_pints(addr_cur_geo) &
            + real(igeo,REALK)*cur_geo_pints(addr_low_geo)
        end do
      end do
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "sub_nucpot_geom_d", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("sub_nucpot_geom_d>> ",A,I6)
110 format("sub_nucpot_geom_d>> ","up",I6,", cur",I6,", low",I6)
120 format("sub_nucpot_geom_d>> ","up",I6,", cur",I6)
#endif
  end subroutine sub_nucpot_geom_d
