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
!!  This file transfers Boys functions to geometric derivatives of Gaussian charge
!!  potential origin with zeroth order HGTOs.
!!
!!  2012-03-06, Bin Gao
!!  * rewrites to use fixed size temporary integrals, the size is determined
!!    from \fn(dim_gaupot_geom)
!!
!!  2010-12-09, Bin Gao
!!  * first version

#include "stdout.h"

  !> \brief transfers Boys functions to geometric derivatives of Gaussian charge
  !>        potential origin with zeroth order HGTOs
  !> \author Bin Gao
  !> \date 2010-12-09
  !> \param coord_bra contains the coordinates of bra center
  !> \param exponent_bra is the exponent of primitive HGTO of bra center
  !> \param coord_ket contains the coordinates of ket center
  !> \param exponent_ket is the exponent of primitive HGTO of ket center
  !> \param gaupot_origin contains the coordinates of Gaussian charge potential origin
  !> \param gaupot_expt is the exponent used in the Gaussian broadening function of the charge
  !> \param scal_const is the scale constant for Gaussian charge potential
  !> \param orders_geo_pot contains the orders of geometric derivatives of
  !>        Gaussian charge potential origin
  !> \param order_elec is the order of electronic derivatives
  !> \param dim_geo_pot is the dimension of integrals
  !> \return geo_pot_pints contains the primitive integrals with zeroth order HGTOs
  subroutine gaupot_geom(coord_bra, exponent_bra, coord_ket, exponent_ket,       &
                         gaupot_origin, gaupot_expt, scal_const, orders_geo_pot, &
                         order_elec, dim_geo_pot, geo_pot_pints)
    use xkind
    implicit none
    real(REALK), intent(in) :: coord_bra(3)
    real(REALK), intent(in) :: exponent_bra
    real(REALK), intent(in) :: coord_ket(3)
    real(REALK), intent(in) :: exponent_ket
    real(REALK), intent(in) :: gaupot_origin(3)
    real(REALK), intent(in) :: gaupot_expt
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: orders_geo_pot(2)
    integer, intent(in) :: order_elec
    integer, intent(in) :: dim_geo_pot
    real(REALK), intent(out) :: geo_pot_pints(dim_geo_pot)
!f2py intent(in) :: coord_bra
!f2py intent(in) :: exponent_bra
!f2py intent(in) :: coord_ket
!f2py intent(in) :: exponent_ket
!f2py intent(in) :: gaupot_origin
!f2py intent(in) :: gaupot_expt
!f2py intent(in) :: scal_const
!f2py intent(in) :: orders_geo_pot
!f2py intent(in) :: order_elec
!f2py intent(in) :: dim_geo_pot
!f2py intent(out) :: geo_pot_pints
!f2py depend(dim_geo_pot) :: geo_pot_pints
#include "private/pi.h"
    real(REALK) scaled_total_expnt  !scaled total exponent
    real(REALK) scale_factor        !scale factor
    real(REALK) recip_total_expnt   !reciprocal of total exponent
    real(REALK) reduced_expnt       !reduced exponent
    real(REALK) gauorg_wrt_cc(3)    !relative coordinates of Gaussian charge potential origin
                                    !w.r.t. center-of-charge
    real(REALK) sd_bra_ket          !square of the relative distance between bra and ket
    integer min_order_boys          !minimum order of Boys functions
    logical odd_min_geo             !indicates if the minimum order of geometric derivatives is odd or even
    real(REALK) arg_boys            !argument of Boys functions
    real(REALK), allocatable :: aux_pints(:)
                                    !auxiliary integrals
    integer max_up_geo              !maximum upper order of geometric derivatives of potential origin
    integer min_up_geo              !minimum upper order of geometric derivatives of potential origin
    integer dim_cur_geo             !dimension of integrals with current order geometric derivatives
    integer dim_up_geo              !dimension of integrals with upper order geometric derivatives
    integer cur_boys                !pointer to current order Boys function in temporary integrals
    integer up_boys                 !pointer of upper order Boys function in temporary integrals
    integer dim_tmp                 !dimension of temporary integrals
    real(REALK), allocatable :: tmp_ints(:,:)
                                    !temporary integrals
    integer order_boys              !current order of Boys function
    integer ierr                    !error information
#if defined(XTIME)
    real(REALK) curr_time           !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! prepares the prerequisite quantities
    scaled_total_expnt = exponent_bra+exponent_ket
    scale_factor = gaupot_expt/(gaupot_expt+scaled_total_expnt)
    recip_total_expnt = 1.0_REALK/scaled_total_expnt
    reduced_expnt = exponent_bra*exponent_ket*recip_total_expnt
    sd_bra_ket = 0.0_REALK
    arg_boys = 0.0_REALK
    do ierr = 1, 3
      gauorg_wrt_cc(ierr) = gaupot_origin(ierr)-recip_total_expnt &
                          * (exponent_bra*coord_bra(ierr)+exponent_ket*coord_ket(ierr))
      sd_bra_ket = sd_bra_ket+(coord_ket(ierr)-coord_bra(ierr))**2
      arg_boys = arg_boys+gauorg_wrt_cc(ierr)**2
    end do
    scaled_total_expnt = scale_factor*scaled_total_expnt
    arg_boys = scaled_total_expnt*arg_boys
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
      call error_stop("gaupot_geom", "failed to allocate aux_pints", &
                      orders_geo_pot(2)-min_order_boys+1)
    ! computes the Boys functions
    call aux_boys_vec(min_order_boys, orders_geo_pot(2), arg_boys, aux_pints)
    ! computes \f$\bar{C}(-2b_{j\lambda})^{|\boldsymbol{n}|}\frac{2\pi}{p_{ij}}%
    !             \mathrm{e}^{-u_{ij}R_{\kappa\lambda}^2}\sqrt{\frac{\rho}{\rho+p_{ij}}}\f$
    arg_boys = scal_const*(-exponent_ket-exponent_ket)**order_elec &
             * exp(-reduced_expnt*sd_bra_ket)*(PI+PI)              &
             * recip_total_expnt*sqrt(scale_factor)
    select case(orders_geo_pot(2))
    ! returns s-shell integrals
    case(0)
      geo_pot_pints(1) = arg_boys*aux_pints(0)
      deallocate(aux_pints)
    ! returns integrals up to p-shell
    case(1)
      ! computes the first order auxiliary integral
      aux_pints(1) = -arg_boys*(scaled_total_expnt+scaled_total_expnt)*aux_pints(1)
      ! returns integrals of s- and p-shells
      if (orders_geo_pot(1)==0) then
        geo_pot_pints(1) = arg_boys*aux_pints(0)          !s
        geo_pot_pints(2) = gauorg_wrt_cc(1)*aux_pints(1)  !px
        geo_pot_pints(3) = gauorg_wrt_cc(2)*aux_pints(1)  !py
        geo_pot_pints(4) = gauorg_wrt_cc(3)*aux_pints(1)  !pz
      ! only returns integrals of p-shell
      else
        geo_pot_pints(1) = gauorg_wrt_cc(1)*aux_pints(1)  !px
        geo_pot_pints(2) = gauorg_wrt_cc(2)*aux_pints(1)  !py
        geo_pot_pints(3) = gauorg_wrt_cc(3)*aux_pints(1)  !pz
      end if
      deallocate(aux_pints)
    ! returns integrals whose maximum order is > 1
    case default
      ! computes \f$\bar{C}(-2b_{j\lambda})^{|\boldsymbol{n}|}\frac{2\pi}{p_{ij}}%
      !             \mathrm{e}^{-u_{ij}R_{\kappa\lambda}^2}\sqrt{\frac{\rho}{\rho+p_{ij}}}%
      !             \left(-2\frac{\rho p_{ij}}{\rho+p_{ij}}\right)^{|\boldsymbol{L}_{C}|_{\max}}\f$
      arg_boys = arg_boys*(-2.0_REALK*scaled_total_expnt)**orders_geo_pot(2)
      ! computes the last auxiliary integral
      aux_pints(orders_geo_pot(2)) = arg_boys*aux_pints(orders_geo_pot(2))
      ! \f$-\frac{1}{2p_{ij}}\f$
      recip_total_expnt = -0.5_REALK/scaled_total_expnt
      ! allocates memory for temporary integrals
      call dim_nucpot_geom(odd_min_geo, orders_geo_pot(2), min_order_boys, dim_tmp)
      allocate(tmp_ints(dim_tmp,2), stat=ierr)
      if (ierr/=0) &
        call error_stop("gaupot_geom", "failed to allocate tmp_ints", dim_tmp*2)
      ! sets the temporary integrals of s- and p-shells
      arg_boys = recip_total_expnt*arg_boys                          !updates the scale factor
      tmp_ints(1,1) = arg_boys*aux_pints(orders_geo_pot(2)-1)        !s
      tmp_ints(2,1) = gauorg_wrt_cc(1)*aux_pints(orders_geo_pot(2))  !px
      tmp_ints(3,1) = gauorg_wrt_cc(2)*aux_pints(orders_geo_pot(2))  !py
      tmp_ints(4,1) = gauorg_wrt_cc(3)*aux_pints(orders_geo_pot(2))  !pz
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
        tmp_ints(2,cur_boys) = gauorg_wrt_cc(1)*tmp_ints(1,up_boys)  !px
        tmp_ints(3,cur_boys) = gauorg_wrt_cc(2)*tmp_ints(1,up_boys)  !py
        tmp_ints(4,cur_boys) = gauorg_wrt_cc(3)*tmp_ints(1,up_boys)  !pz
        ! performs the recurrecen relation for other shells
        call sub_nucpot_geom_d(2, max_up_geo, gauorg_wrt_cc, dim_cur_geo,     &
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
        tmp_ints(1,cur_boys) = gauorg_wrt_cc(1)*tmp_ints(1,up_boys)  !px
        tmp_ints(2,cur_boys) = gauorg_wrt_cc(2)*tmp_ints(1,up_boys)  !py
        tmp_ints(3,cur_boys) = gauorg_wrt_cc(3)*tmp_ints(1,up_boys)  !pz
        ! performs the recurrecen relation for other shells
        call sub_nucpot_geom_d(2, max_up_geo, gauorg_wrt_cc, dim_cur_geo,     &
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
        call sub_nucpot_geom_d(min_up_geo, max_up_geo, gauorg_wrt_cc,        &
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
    call xtimer_view(curr_time, "gaupot_geom", STDOUT)
#endif
    return
  end subroutine gaupot_geom
