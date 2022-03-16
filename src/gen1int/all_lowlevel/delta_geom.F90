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
!!  This file calculates the geometric derivatives of Dirac delta function
!!  with zeroth order geometric derivatives.
!!
!!  2012-03-15, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief calculates the geometric derivatives of Dirac delta function with zeroth order geometric derivatives
  !> \author Bin Gao
  !> \date 2012-03-15
  !> \param coord_bra contains the coordinates of bra center
  !> \param exponent_bra is the exponent of primitive geometric derivative of bra center
  !> \param coord_ket contains the coordinates of ket center
  !> \param exponent_ket is the exponent of primitive geometric derivative of ket center
  !> \param delta_origin contains the coordinates of Dirac delta function origin
  !> \param scal_const is the scale constant for Dirac delta function
  !> \param orders_geo_pot contains the orders of geometric derivatives of
  !>        Dirac delta function origin
  !> \param order_elec is the order of electronic derivatives
  !> \param dim_geo_pot is the dimension of integrals
  !> \return geo_pot_pints contains the primitive integrals with zeroth order geometric derivatives
  subroutine delta_geom(coord_bra, exponent_bra, coord_ket, exponent_ket,     &
                        delta_origin, scal_const, orders_geo_pot, order_elec, &
                        dim_geo_pot, geo_pot_pints)
    use xkind
    implicit none
    real(REALK), intent(in) :: coord_bra(3)
    real(REALK), intent(in) :: exponent_bra
    real(REALK), intent(in) :: coord_ket(3)
    real(REALK), intent(in) :: exponent_ket
    real(REALK), intent(in) :: delta_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: orders_geo_pot(2)
    integer, intent(in) :: order_elec
    integer, intent(in) :: dim_geo_pot
    real(REALK), intent(out) :: geo_pot_pints(dim_geo_pot)
!f2py intent(in) :: coord_bra
!f2py intent(in) :: exponent_bra
!f2py intent(in) :: coord_ket
!f2py intent(in) :: exponent_ket
!f2py intent(in) :: delta_origin
!f2py intent(in) :: scal_const
!f2py intent(in) :: orders_geo_pot
!f2py intent(in) :: order_elec
!f2py intent(in) :: dim_geo_pot
!f2py intent(out) :: geo_pot_pints
!f2py depend(dim_geo_pot) :: geo_pot_pints
    real(REALK) neg_total_expnt              !negative total exponent
    real(REALK) recip_total_expnt            !reciprocal of total exponent
    real(REALK) neg_reduced_expnt            !negative reduced exponent
    real(REALK) sd_bra_ket                   !square of the relative distance between bra and ket centers
    real(REALK) delta_wrt_cc(3)              !relative coordinates of delta function w.r.t. center-of-charge
    real(REALK) sd_delta_cc                  !square of the relative distance between delta function ...
                                             !and center-of-charge
    real(REALK) recur_param                  !parameter \f$-2p_{ij}\f$ used in recurrence relations
    integer base_low_geo                     !base address of lower order geometric derivatives
    integer base_cur_geo                     !base address of current order geometric derivatives
    integer base_up_geo                      !base address of upper order geometric derivatives
    integer order_geo, iorder, jorder        !incremental recorders over orders of geometric derivatives
    integer addr_low_geo                     !address of lower order geometric derivatives
    integer addr_cur_geo                     !address of current order geometric derivatives
    integer addr_up_geo                      !address of upper order geometric derivatives
    integer num_min_geo                      !number of minimum returned order geometric derivatives
    integer dim_tmp                          !dimension of temporary integrals
    real(REALK), allocatable :: tmp_ints(:)  !temporary integrals
    real(REALK) zero_pint                    !integral with zeroth order geometric derivative
    integer ierr                             !error information
#if defined(XTIME)
    real(REALK) curr_time                    !current CPU time   
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! computes the negative total exponent, reciprocal of total exponent and negative reduced exponent
    neg_total_expnt = -(exponent_bra+exponent_ket)
    recip_total_expnt = 1.0_REALK/neg_total_expnt
    neg_reduced_expnt = exponent_bra*exponent_ket*recip_total_expnt
    recip_total_expnt = -recip_total_expnt
    ! computes the relative coordinates of delta function w.r.t. center-of-charge,
    ! and squares of the relative distance between bra and ket centers, delta function
    ! and center-of-charge
    sd_bra_ket = 0.0_REALK
    sd_delta_cc = 0.0_REALK
    do iorder = 1, 3
      delta_wrt_cc(iorder) = delta_origin(iorder)-(exponent_bra*coord_bra(iorder) &
                           + exponent_ket*coord_ket(iorder))*recip_total_expnt
      sd_bra_ket = sd_bra_ket+(coord_ket(iorder)-coord_bra(iorder))**2
      sd_delta_cc = sd_delta_cc+delta_wrt_cc(iorder)**2
    end do
    select case(orders_geo_pot(1))
    ! zeroth order geometric derivative returned
    case(0)
      select case(orders_geo_pot(2))
      ! only zeroth order geometric derivatives returned
      case(0)
        geo_pot_pints(1) = scal_const*(-exponent_ket-exponent_ket)**order_elec &
                         * exp(neg_reduced_expnt*sd_bra_ket+neg_total_expnt*sd_delta_cc)
      ! zeroth and first order geometric derivatives returned
      case(1)
        ! zeroth order
        geo_pot_pints(1) = scal_const*(-exponent_ket-exponent_ket)**order_elec &
                         * exp(neg_reduced_expnt*sd_bra_ket+neg_total_expnt*sd_delta_cc)
        ! first order
        recur_param = neg_total_expnt+neg_total_expnt
        geo_pot_pints(2) = recur_param*delta_wrt_cc(1)*geo_pot_pints(1)
        geo_pot_pints(3) = recur_param*delta_wrt_cc(2)*geo_pot_pints(1)
        geo_pot_pints(4) = recur_param*delta_wrt_cc(3)*geo_pot_pints(1)
      ! up to, at least, second order geometric derivatives returned
      case default
        ! zeroth order
        geo_pot_pints(1) = scal_const*(-exponent_ket-exponent_ket)**order_elec &
                         * exp(neg_reduced_expnt*sd_bra_ket+neg_total_expnt*sd_delta_cc)
        ! first order
        recur_param = neg_total_expnt+neg_total_expnt
        geo_pot_pints(2) = recur_param*delta_wrt_cc(1)*geo_pot_pints(1)
        geo_pot_pints(3) = recur_param*delta_wrt_cc(2)*geo_pot_pints(1)
        geo_pot_pints(4) = recur_param*delta_wrt_cc(3)*geo_pot_pints(1)
        ! initializes base addresses
        base_low_geo = 0
        base_cur_geo = 1
        base_up_geo = 4
        ! left returned orders up to \var(orders_geo_pot(2))
        do order_geo = 1, orders_geo_pot(2)-1
          ! recurrence relations along x-direction
          addr_up_geo = base_up_geo+1
          addr_cur_geo = base_cur_geo+1
          geo_pot_pints(addr_up_geo)                                   &
            = recur_param*(delta_wrt_cc(1)*geo_pot_pints(addr_cur_geo) &
            + real(order_geo,REALK)*geo_pot_pints(base_low_geo+1))
          ! recurrence relations along y-direction
          addr_up_geo = addr_up_geo+1
          geo_pot_pints(addr_up_geo) = recur_param*delta_wrt_cc(2) &
                                     * geo_pot_pints(addr_cur_geo)
          do iorder = 1, order_geo
            addr_up_geo = addr_up_geo+1
            addr_cur_geo = addr_cur_geo+1
            geo_pot_pints(addr_up_geo)                                   &
              = recur_param*(delta_wrt_cc(2)*geo_pot_pints(addr_cur_geo) &
              + real(iorder,REALK)*geo_pot_pints(base_low_geo+iorder))
          end do
          ! recurrence relations along z-direction
          addr_cur_geo = base_cur_geo
          do jorder = 0, order_geo
            addr_up_geo = addr_up_geo+1
            addr_cur_geo = addr_cur_geo+1
            geo_pot_pints(addr_up_geo) = recur_param*delta_wrt_cc(3) &
                                       * geo_pot_pints(addr_cur_geo)
          end do
          addr_low_geo = base_low_geo
          do iorder = 1, order_geo
            do jorder = 0, order_geo-iorder
              addr_up_geo = addr_up_geo+1
              addr_cur_geo = addr_cur_geo+1
              addr_low_geo = addr_low_geo+1
              geo_pot_pints(addr_up_geo)                                   &
                = recur_param*(delta_wrt_cc(3)*geo_pot_pints(addr_cur_geo) &
                + real(iorder,REALK)*geo_pot_pints(addr_low_geo))
            end do
          end do
          ! updates base addresses
          base_low_geo = base_cur_geo
          base_cur_geo = base_up_geo
          base_up_geo = addr_up_geo
        end do
      end select
    ! first order geometric derivatives returned
    case (1)
      ! parameter \f$-2p_{ij}\f$ used in recurrence relations
      recur_param = neg_total_expnt+neg_total_expnt
      select case(orders_geo_pot(2))
      ! only first order geometric derivatives returned
      case(1)
        ! zeroth order
        geo_pot_pints(1) = scal_const*(-exponent_ket-exponent_ket)**order_elec &
                         * exp(neg_reduced_expnt*sd_bra_ket+neg_total_expnt*sd_delta_cc)
        ! first order
        geo_pot_pints(2) = recur_param*delta_wrt_cc(2)*geo_pot_pints(1)
        geo_pot_pints(3) = recur_param*delta_wrt_cc(3)*geo_pot_pints(1)
        geo_pot_pints(1) = recur_param*delta_wrt_cc(1)*geo_pot_pints(1)
      ! first and second orders geometric derivatives returned
      case(2)
        ! computes the integral with zeroth order geometric derivative
        zero_pint = scal_const*(-exponent_ket-exponent_ket)**order_elec &
                  * exp(neg_reduced_expnt*sd_bra_ket+neg_total_expnt*sd_delta_cc)
        ! scales the integral with zeroth order geometric derivative by \f$-2p_{ij}\f$
        zero_pint = recur_param*zero_pint
        ! first order
        do iorder = 1, 3
          geo_pot_pints(iorder) = delta_wrt_cc(iorder)*zero_pint
          ! scales the relative coordinates of delta function w.r.t. center-of-charge by \f$-2p_{ij}\f$
          delta_wrt_cc(iorder) = recur_param*delta_wrt_cc(iorder)
        end do
        ! second order
        geo_pot_pints(4) = delta_wrt_cc(1)*geo_pot_pints(1)+zero_pint
        geo_pot_pints(5) = delta_wrt_cc(2)*geo_pot_pints(1)
        geo_pot_pints(6) = delta_wrt_cc(2)*geo_pot_pints(2)+zero_pint
        geo_pot_pints(7) = delta_wrt_cc(3)*geo_pot_pints(1)
        geo_pot_pints(8) = delta_wrt_cc(3)*geo_pot_pints(2)
        geo_pot_pints(9) = delta_wrt_cc(3)*geo_pot_pints(3)+zero_pint
      ! up to, at least, third order geometric derivatives returned
      case default
        ! computes the integral with zeroth order geometric derivative
        zero_pint = scal_const*(-exponent_ket-exponent_ket)**order_elec &
                  * exp(neg_reduced_expnt*sd_bra_ket+neg_total_expnt*sd_delta_cc)
        ! scales the integral with zeroth order geometric derivative by \f$-2p_{ij}\f$
        zero_pint = recur_param*zero_pint
        ! first order
        do iorder = 1, 3
          geo_pot_pints(iorder) = delta_wrt_cc(iorder)*zero_pint
        end do
        ! second order
        geo_pot_pints(4) = recur_param*delta_wrt_cc(1)*geo_pot_pints(1)+zero_pint
        geo_pot_pints(5) = recur_param*delta_wrt_cc(2)*geo_pot_pints(1)
        geo_pot_pints(6) = recur_param*delta_wrt_cc(2)*geo_pot_pints(2)+zero_pint
        geo_pot_pints(7) = recur_param*delta_wrt_cc(3)*geo_pot_pints(1)
        geo_pot_pints(8) = recur_param*delta_wrt_cc(3)*geo_pot_pints(2)
        geo_pot_pints(9) = recur_param*delta_wrt_cc(3)*geo_pot_pints(3)+zero_pint
        ! initializes base addresses
        base_low_geo = 0
        base_cur_geo = 3
        base_up_geo = 9
        ! left returned orders up to \var(orders_geo_pot(2))
        do order_geo = 2, orders_geo_pot(2)-1
          ! recurrence relations along x-direction
          addr_up_geo = base_up_geo+1
          addr_cur_geo = base_cur_geo+1
          geo_pot_pints(addr_up_geo)                                   &
            = recur_param*(delta_wrt_cc(1)*geo_pot_pints(addr_cur_geo) &
            + real(order_geo,REALK)*geo_pot_pints(base_low_geo+1))
          ! recurrence relations along y-direction
          addr_up_geo = addr_up_geo+1
          geo_pot_pints(addr_up_geo) = recur_param*delta_wrt_cc(2) &
                                     * geo_pot_pints(addr_cur_geo)
          do iorder = 1, order_geo
            addr_up_geo = addr_up_geo+1
            addr_cur_geo = addr_cur_geo+1
            geo_pot_pints(addr_up_geo)                                   &
              = recur_param*(delta_wrt_cc(2)*geo_pot_pints(addr_cur_geo) &
              + real(iorder,REALK)*geo_pot_pints(base_low_geo+iorder))
          end do
          ! recurrence relations along z-direction
          addr_cur_geo = base_cur_geo
          do jorder = 0, order_geo
            addr_up_geo = addr_up_geo+1
            addr_cur_geo = addr_cur_geo+1
            geo_pot_pints(addr_up_geo) = recur_param*delta_wrt_cc(3) &
                                       * geo_pot_pints(addr_cur_geo)
          end do
          addr_low_geo = base_low_geo
          do iorder = 1, order_geo
            do jorder = 0, order_geo-iorder
              addr_up_geo = addr_up_geo+1
              addr_cur_geo = addr_cur_geo+1
              addr_low_geo = addr_low_geo+1
              geo_pot_pints(addr_up_geo)                                   &
                = recur_param*(delta_wrt_cc(3)*geo_pot_pints(addr_cur_geo) &
                + real(iorder,REALK)*geo_pot_pints(addr_low_geo))
            end do
          end do
          ! updates base addresses
          base_low_geo = base_cur_geo
          base_cur_geo = base_up_geo
          base_up_geo = addr_up_geo
        end do
      end select
    ! high order geometric derivatives (>1) returned
    case default
      ! parameter \f$-2p_{ij}\f$ used in recurrence relations
      recur_param = neg_total_expnt+neg_total_expnt
      ! only order \var(orders_geo_pot(1)) geometric derivatives returned
      if (orders_geo_pot(1)==orders_geo_pot(2)) then
        dim_tmp = (orders_geo_pot(1)+1)*(orders_geo_pot(1)+2)/2
        allocate(tmp_ints(3*dim_tmp), stat=ierr)
        if (ierr/=0)                                     &
          call error_stop("delta_geom",                  &
                          "failed to allocate tmp_ints", &
                          3*dim_tmp)
        ! computes the integral with zeroth order geometric derivative
        tmp_ints(1) = scal_const*(-exponent_ket-exponent_ket)**order_elec &
                    * exp(neg_reduced_expnt*sd_bra_ket+neg_total_expnt*sd_delta_cc)
        ! first order
        do iorder = 1, 3
          tmp_ints(dim_tmp+iorder) = recur_param*delta_wrt_cc(iorder)*tmp_ints(1)
        end do
        ! initializes base addresses
        base_low_geo = 0
        base_cur_geo = dim_tmp
        base_up_geo = 2*dim_tmp
        do order_geo = 1, orders_geo_pot(1)-1
          ! recurrence relations along x-direction
          addr_up_geo = base_up_geo+1
          addr_cur_geo = base_cur_geo+1
          tmp_ints(addr_up_geo)                                   &
            = recur_param*(delta_wrt_cc(1)*tmp_ints(addr_cur_geo) &
            + real(order_geo,REALK)*tmp_ints(base_low_geo+1))
          ! recurrence relations along y-direction
          addr_up_geo = addr_up_geo+1
          tmp_ints(addr_up_geo) = recur_param*delta_wrt_cc(2)*tmp_ints(addr_cur_geo)
          do iorder = 1, order_geo
            addr_up_geo = addr_up_geo+1
            addr_cur_geo = addr_cur_geo+1
            tmp_ints(addr_up_geo)                                   &
              = recur_param*(delta_wrt_cc(2)*tmp_ints(addr_cur_geo) &
              + real(iorder,REALK)*tmp_ints(base_low_geo+iorder))
          end do
          ! recurrence relations along z-direction
          addr_cur_geo = base_cur_geo
          do jorder = 0, order_geo
            addr_up_geo = addr_up_geo+1
            addr_cur_geo = addr_cur_geo+1
            tmp_ints(addr_up_geo) = recur_param*delta_wrt_cc(3)*tmp_ints(addr_cur_geo)
          end do
          addr_low_geo = base_low_geo
          do iorder = 1, order_geo
            do jorder = 0, order_geo-iorder
              addr_up_geo = addr_up_geo+1
              addr_cur_geo = addr_cur_geo+1
              addr_low_geo = addr_low_geo+1
              tmp_ints(addr_up_geo)                                   &
                = recur_param*(delta_wrt_cc(3)*tmp_ints(addr_cur_geo) &
                + real(iorder,REALK)*tmp_ints(addr_low_geo))
            end do
          end do
          ! updates base addresses
          addr_up_geo = base_low_geo
          base_low_geo = base_cur_geo
          base_cur_geo = base_up_geo
          base_up_geo = addr_up_geo
        end do
        ! assigns order \var(orders_geo_pot(1))
        geo_pot_pints(:) = tmp_ints(base_cur_geo+1:base_cur_geo+dim_tmp)
        deallocate(tmp_ints)
      else
        num_min_geo = (orders_geo_pot(1)+1)*(orders_geo_pot(1)+2)/2
        dim_tmp = num_min_geo+orders_geo_pot(1)+2
        allocate(tmp_ints(3*dim_tmp), stat=ierr)
        if (ierr/=0)                                     &
          call error_stop("delta_geom",                  &
                          "failed to allocate tmp_ints", &
                          3*dim_tmp)
        ! computes the integral with zeroth order geometric derivative
        tmp_ints(1) = scal_const*(-exponent_ket-exponent_ket)**order_elec &
                    * exp(neg_reduced_expnt*sd_bra_ket+neg_total_expnt*sd_delta_cc)
        ! first order
        do iorder = 1, 3
          tmp_ints(dim_tmp+iorder) = recur_param*delta_wrt_cc(iorder)*tmp_ints(1)
        end do
        ! initializes base addresses
        base_low_geo = 0
        base_cur_geo = dim_tmp
        base_up_geo = 2*dim_tmp
        ! orders up to \var(orders_geo_pot(1))+1
        do order_geo = 1, orders_geo_pot(1)
          ! recurrence relations along x-direction
          addr_up_geo = base_up_geo+1
          addr_cur_geo = base_cur_geo+1
          tmp_ints(addr_up_geo)                                   &
            = recur_param*(delta_wrt_cc(1)*tmp_ints(addr_cur_geo) &
            + real(order_geo,REALK)*tmp_ints(base_low_geo+1))
          ! recurrence relations along y-direction
          addr_up_geo = addr_up_geo+1
          tmp_ints(addr_up_geo) = recur_param*delta_wrt_cc(2)*tmp_ints(addr_cur_geo)
          do iorder = 1, order_geo
            addr_up_geo = addr_up_geo+1
            addr_cur_geo = addr_cur_geo+1
            tmp_ints(addr_up_geo)                                   &
              = recur_param*(delta_wrt_cc(2)*tmp_ints(addr_cur_geo) &
              + real(iorder,REALK)*tmp_ints(base_low_geo+iorder))
          end do
          ! recurrence relations along z-direction
          addr_cur_geo = base_cur_geo
          do jorder = 0, order_geo
            addr_up_geo = addr_up_geo+1
            addr_cur_geo = addr_cur_geo+1
            tmp_ints(addr_up_geo) = recur_param*delta_wrt_cc(3)*tmp_ints(addr_cur_geo)
          end do
          addr_low_geo = base_low_geo
          do iorder = 1, order_geo
            do jorder = 0, order_geo-iorder
              addr_up_geo = addr_up_geo+1
              addr_cur_geo = addr_cur_geo+1
              addr_low_geo = addr_low_geo+1
              tmp_ints(addr_up_geo)                                   &
                = recur_param*(delta_wrt_cc(3)*tmp_ints(addr_cur_geo) &
                + real(iorder,REALK)*tmp_ints(addr_low_geo))
            end do
          end do
          ! updates base addresses
          addr_up_geo = base_low_geo
          base_low_geo = base_cur_geo
          base_cur_geo = base_up_geo
          base_up_geo = addr_up_geo
        end do
        ! assigns orders \var(orders_geo_pot(1)) and \var(orders_geo_pot(1))+1,
        ! and initializes base addresses
        geo_pot_pints(1:num_min_geo) &
          = tmp_ints(base_low_geo+1:base_low_geo+num_min_geo)
        base_up_geo = num_min_geo+dim_tmp
        geo_pot_pints(num_min_geo+1:base_up_geo) &
          = tmp_ints(base_cur_geo+1:base_cur_geo+dim_tmp)
        deallocate(tmp_ints)
        base_low_geo = 0
        base_cur_geo = num_min_geo
        ! left returned orders up to \var(orders_geo_pot(2))
        do order_geo = orders_geo_pot(1)+1, orders_geo_pot(2)-1
          ! recurrence relations along x-direction
          addr_up_geo = base_up_geo+1
          addr_cur_geo = base_cur_geo+1
          geo_pot_pints(addr_up_geo)                                   &
            = recur_param*(delta_wrt_cc(1)*geo_pot_pints(addr_cur_geo) &
            + real(order_geo,REALK)*geo_pot_pints(base_low_geo+1))
          ! recurrence relations along y-direction
          addr_up_geo = addr_up_geo+1
          geo_pot_pints(addr_up_geo) = recur_param*delta_wrt_cc(2) &
                                     * geo_pot_pints(addr_cur_geo)
          do iorder = 1, order_geo
            addr_up_geo = addr_up_geo+1
            addr_cur_geo = addr_cur_geo+1
            geo_pot_pints(addr_up_geo)                                   &
              = recur_param*(delta_wrt_cc(2)*geo_pot_pints(addr_cur_geo) &
              + real(iorder,REALK)*geo_pot_pints(base_low_geo+iorder))
          end do
          ! recurrence relations along z-direction
          addr_cur_geo = base_cur_geo
          do jorder = 0, order_geo
            addr_up_geo = addr_up_geo+1
            addr_cur_geo = addr_cur_geo+1
            geo_pot_pints(addr_up_geo) = recur_param*delta_wrt_cc(3) &
                                       * geo_pot_pints(addr_cur_geo)
          end do
          addr_low_geo = base_low_geo
          do iorder = 1, order_geo
            do jorder = 0, order_geo-iorder
              addr_up_geo = addr_up_geo+1
              addr_cur_geo = addr_cur_geo+1
              addr_low_geo = addr_low_geo+1
              geo_pot_pints(addr_up_geo)                                   &
                = recur_param*(delta_wrt_cc(3)*geo_pot_pints(addr_cur_geo) &
                + real(iorder,REALK)*geo_pot_pints(addr_low_geo))
            end do
          end do
          ! updates base addresses
          base_low_geo = base_cur_geo
          base_cur_geo = base_up_geo
          base_up_geo = addr_up_geo
        end do
      end if
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "delta_geom", STDOUT)
#endif
    return
  end subroutine delta_geom
