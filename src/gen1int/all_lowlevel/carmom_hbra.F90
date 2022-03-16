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
!!  This file returns the zeroth order Cartesian multipole moment integrals
!!  with specific orders of HGTOs on bra center.
!!
!!  2012-02-12, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief computes the zeroth order Cartesian multipole moment integrals
  !>        with specific orders of HGTOs on bra center
  !> \author Bin Gao
  !> \date 2012-02-12
  !> \param orders_hgto_bra contains the minimum and maximum orders of HGTOs on bra center
  !> \param coord_bra contains the coordinates of bra center
  !> \param exponent_bra is the exponent of primitive HGTO of bra center
  !> \param coord_ket contains the coordinates of ket center
  !> \param exponent_ket is the exponent of primitive HGTO of ket center
  !> \param order_elec is the order of electronic derivatives
  !> \param scal_const is the scale constant for Cartesian multipole moments
  !> \param dim_hgto_bra is the dimension of integrals
  !> \return hbra_pints contains the primitive HGTO integrals
  subroutine carmom_hbra(orders_hgto_bra, coord_bra, exponent_bra, &
                         coord_ket, exponent_ket, order_elec,      &
                         scal_const, dim_hgto_bra, hbra_pints)
    use xkind
    implicit none
    integer, intent(in) :: orders_hgto_bra(2)
    real(REALK), intent(in) :: coord_bra(3)
    real(REALK), intent(in) :: exponent_bra
    real(REALK), intent(in) :: coord_ket(3)
    real(REALK), intent(in) :: exponent_ket
    integer, intent(in) :: order_elec
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: dim_hgto_bra
    real(REALK), intent(out) :: hbra_pints(dim_hgto_bra)
!f2py intent(in) :: orders_hgto_bra
!f2py intent(in) :: coord_bra
!f2py intent(in) :: exponent_bra
!f2py intent(in) :: coord_ket
!f2py intent(in) :: exponent_ket
!f2py intent(in) :: order_elec
!f2py intent(in) :: scal_const
!f2py intent(in) :: dim_hgto_bra
!f2py intent(out) :: hbra_pints
!f2py depend(dim_hgto_bra) :: hbra_pints
#include "private/pi.h"
    real(REALK) recip_total_expnt            !reciprocal of total exponent
    real(REALK) reduced_expnt                !reduced exponent
    real(REALK) sd_bra_ket                   !square of the relative distance between bra and ket centers
    real(REALK) cc_wrt_bra(3)                !relative coordinates of center-of-charge w.r.t. bra center
    real(REALK) half_recip_nup               !\f$-\frac{1}{2p_{ij}}\frac{b_{j\lambda}}{a_{i\kappa}}\f$
    integer base_low_hgto                    !base address of lower order HGTOs
    integer base_cur_hgto                    !base address of current order HGTOs
    integer base_up_hgto                     !base address of upper order HGTOs
    integer order_hgto, iorder, jorder       !incremental recorders over orders of HGTOs
    integer addr_low_hgto                    !address of lower order HGTOs
    integer addr_cur_hgto                    !address of current order HGTOs
    integer addr_up_hgto                     !address of upper order HGTOs
    integer num_min_hgto                     !number of minimum returned order HGTOs
    integer dim_tmp                          !dimension of temporary integrals
    real(REALK), allocatable :: tmp_ints(:)  !temporary integrals
    real(REALK) zero_pint                    !integral with zeroth order HGTO
    integer ierr                             !error information
#if defined(XTIME)
    real(REALK) curr_time                    !current CPU time   
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! computes the reciprocal of total exponent
    recip_total_expnt = 1.0_REALK/(exponent_bra+exponent_ket)
    ! computes the reduced exponent
    reduced_expnt = exponent_bra*exponent_ket*recip_total_expnt
    ! computes the square of the relative distance between bra and ket centers
    sd_bra_ket = 0.0_REALK
    do iorder = 1, 3
      sd_bra_ket = sd_bra_ket+(coord_ket(iorder)-coord_bra(iorder))**2
    end do
    select case(orders_hgto_bra(1))
    ! zeroth order HGTOs returned
    case(0)
      select case(orders_hgto_bra(2))
      ! only zeroth order HGTO returned
      case(0)
        hbra_pints(1) = scal_const*(-exponent_ket-exponent_ket)**order_elec &
                      * exp(-reduced_expnt*sd_bra_ket)*sqrt((PI*recip_total_expnt)**3)
      ! zeroth and first order HGTOs returned
      case(1)
        ! zeroth order
        hbra_pints(1) = scal_const*(-exponent_ket-exponent_ket)**order_elec &
                      * exp(-reduced_expnt*sd_bra_ket)*sqrt((PI*recip_total_expnt)**3)
        ! computes parameter \f$-\frac{b_{j\lambda}}{p_{ij}}\f$
        half_recip_nup = -exponent_ket*recip_total_expnt
        ! computes the relative coordinates of center-of-charge w.r.t. bra center
        do iorder = 1, 3
          cc_wrt_bra(iorder) = half_recip_nup*(coord_bra(iorder)-coord_ket(iorder))
        end do
        ! first order
        hbra_pints(2) = cc_wrt_bra(1)*hbra_pints(1)
        hbra_pints(3) = cc_wrt_bra(2)*hbra_pints(1)
        hbra_pints(4) = cc_wrt_bra(3)*hbra_pints(1)
      ! up to, at least, second order HGTOs returned
      case default
        ! zeroth order
        hbra_pints(1) = scal_const*(-exponent_ket-exponent_ket)**order_elec &
                      * exp(-reduced_expnt*sd_bra_ket)*sqrt((PI*recip_total_expnt)**3)
        ! computes parameter \f$-\frac{b_{j\lambda}}{p_{ij}}\f$
        half_recip_nup = -exponent_ket*recip_total_expnt
        ! computes the relative coordinates of center-of-charge w.r.t. bra center
        do iorder = 1, 3
          cc_wrt_bra(iorder) = half_recip_nup*(coord_bra(iorder)-coord_ket(iorder))
        end do
        ! first order
        hbra_pints(2) = cc_wrt_bra(1)*hbra_pints(1)
        hbra_pints(3) = cc_wrt_bra(2)*hbra_pints(1)
        hbra_pints(4) = cc_wrt_bra(3)*hbra_pints(1)
        ! computes parameter \f$-\frac{1}{2p_{ij}}\frac{b_{j\lambda}}{a_{i\kappa}}\f$
        half_recip_nup = 0.5_REALK*half_recip_nup/exponent_bra
        ! initializes base addresses
        base_low_hgto = 0
        base_cur_hgto = 1
        base_up_hgto = 4
        ! left returned orders up to \var(orders_hgto_bra(2))
        do order_hgto = 1, orders_hgto_bra(2)-1
          ! recurrence relations along x-direction
          addr_up_hgto = base_up_hgto+1
          addr_cur_hgto = base_cur_hgto+1
          hbra_pints(addr_up_hgto)                    &
            = cc_wrt_bra(1)*hbra_pints(addr_cur_hgto) &
            + real(order_hgto,REALK)*half_recip_nup*hbra_pints(base_low_hgto+1)
          ! recurrence relations along y-direction
          addr_up_hgto = addr_up_hgto+1
          hbra_pints(addr_up_hgto) = cc_wrt_bra(2)*hbra_pints(addr_cur_hgto)
          do iorder = 1, order_hgto
            addr_up_hgto = addr_up_hgto+1
            addr_cur_hgto = addr_cur_hgto+1
            hbra_pints(addr_up_hgto)                    &
              = cc_wrt_bra(2)*hbra_pints(addr_cur_hgto) &
              + real(iorder,REALK)*half_recip_nup*hbra_pints(base_low_hgto+iorder)
          end do
          ! recurrence relations along z-direction
          addr_cur_hgto = base_cur_hgto
          do jorder = 0, order_hgto
            addr_up_hgto = addr_up_hgto+1
            addr_cur_hgto = addr_cur_hgto+1
            hbra_pints(addr_up_hgto) = cc_wrt_bra(3)*hbra_pints(addr_cur_hgto)
          end do
          addr_low_hgto = base_low_hgto
          do iorder = 1, order_hgto
            do jorder = 0, order_hgto-iorder
              addr_up_hgto = addr_up_hgto+1
              addr_cur_hgto = addr_cur_hgto+1
              addr_low_hgto = addr_low_hgto+1
              hbra_pints(addr_up_hgto)                    &
                = cc_wrt_bra(3)*hbra_pints(addr_cur_hgto) &
                + real(iorder,REALK)*half_recip_nup*hbra_pints(addr_low_hgto)
            end do
          end do
          ! updates base addresses
          base_low_hgto = base_cur_hgto
          base_cur_hgto = base_up_hgto
          base_up_hgto = addr_up_hgto
        end do
      end select
    ! first order HGTOs returned
    case (1)
      select case(orders_hgto_bra(2))
      ! only first order HGTOs returned
      case(1)
        ! zeroth order
        hbra_pints(1) = scal_const*(-exponent_ket-exponent_ket)**order_elec &
                      * exp(-reduced_expnt*sd_bra_ket)*sqrt((PI*recip_total_expnt)**3)
        ! computes parameter \f$-\frac{b_{j\lambda}}{p_{ij}}\f$
        half_recip_nup = -exponent_ket*recip_total_expnt
        ! computes the relative coordinates of center-of-charge w.r.t. bra center
        do iorder = 1, 3
          cc_wrt_bra(iorder) = half_recip_nup*(coord_bra(iorder)-coord_ket(iorder))
        end do
        ! first order
        hbra_pints(2) = cc_wrt_bra(2)*hbra_pints(1)
        hbra_pints(3) = cc_wrt_bra(3)*hbra_pints(1)
        hbra_pints(1) = cc_wrt_bra(1)*hbra_pints(1)
      ! first and second orders HGTOs returned
      case(2)
        ! computes the integral with zeroth order HGTO
        zero_pint = scal_const*(-exponent_ket-exponent_ket)**order_elec &
                  * exp(-reduced_expnt*sd_bra_ket)*sqrt((PI*recip_total_expnt)**3)
        ! computes parameter \f$-\frac{b_{j\lambda}}{p_{ij}}\f$
        half_recip_nup = -exponent_ket*recip_total_expnt
        ! computes the relative coordinates of center-of-charge w.r.t. bra center
        do iorder = 1, 3
          cc_wrt_bra(iorder) = half_recip_nup*(coord_bra(iorder)-coord_ket(iorder))
        end do
        ! computes parameter \f$-\frac{1}{2p_{ij}}\frac{b_{j\lambda}}{a_{i\kappa}}\f$
        half_recip_nup = 0.5_REALK*half_recip_nup/exponent_bra
        ! first order
        hbra_pints(1) = cc_wrt_bra(1)*zero_pint
        hbra_pints(2) = cc_wrt_bra(2)*zero_pint
        hbra_pints(3) = cc_wrt_bra(3)*zero_pint
        ! second order
        hbra_pints(4) = cc_wrt_bra(1)*hbra_pints(1)+half_recip_nup*zero_pint
        hbra_pints(5) = cc_wrt_bra(2)*hbra_pints(1)
        hbra_pints(6) = cc_wrt_bra(2)*hbra_pints(2)+half_recip_nup*zero_pint
        hbra_pints(7) = cc_wrt_bra(3)*hbra_pints(1)
        hbra_pints(8) = cc_wrt_bra(3)*hbra_pints(2)
        hbra_pints(9) = cc_wrt_bra(3)*hbra_pints(3)+half_recip_nup*zero_pint
      ! up to, at least, third order HGTOs returned
      case default
        ! computes the integral with zeroth order HGTO
        zero_pint = scal_const*(-exponent_ket-exponent_ket)**order_elec &
                  * exp(-reduced_expnt*sd_bra_ket)*sqrt((PI*recip_total_expnt)**3)
        ! computes parameter \f$-\frac{b_{j\lambda}}{p_{ij}}\f$
        half_recip_nup = -exponent_ket*recip_total_expnt
        ! computes the relative coordinates of center-of-charge w.r.t. bra center
        do iorder = 1, 3
          cc_wrt_bra(iorder) = half_recip_nup*(coord_bra(iorder)-coord_ket(iorder))
        end do
        ! computes parameter \f$-\frac{1}{2p_{ij}}\frac{b_{j\lambda}}{a_{i\kappa}}\f$
        half_recip_nup = 0.5_REALK*half_recip_nup/exponent_bra
        ! first order
        hbra_pints(1) = cc_wrt_bra(1)*zero_pint
        hbra_pints(2) = cc_wrt_bra(2)*zero_pint
        hbra_pints(3) = cc_wrt_bra(3)*zero_pint
        ! second order
        hbra_pints(4) = cc_wrt_bra(1)*hbra_pints(1)+half_recip_nup*zero_pint
        hbra_pints(5) = cc_wrt_bra(2)*hbra_pints(1)
        hbra_pints(6) = cc_wrt_bra(2)*hbra_pints(2)+half_recip_nup*zero_pint
        hbra_pints(7) = cc_wrt_bra(3)*hbra_pints(1)
        hbra_pints(8) = cc_wrt_bra(3)*hbra_pints(2)
        hbra_pints(9) = cc_wrt_bra(3)*hbra_pints(3)+half_recip_nup*zero_pint
        ! initializes base addresses
        base_low_hgto = 0
        base_cur_hgto = 3
        base_up_hgto = 9
        ! left returned orders up to \var(orders_hgto_bra(2))
        do order_hgto = 2, orders_hgto_bra(2)-1
          ! recurrence relations along x-direction
          addr_up_hgto = base_up_hgto+1
          addr_cur_hgto = base_cur_hgto+1
          hbra_pints(addr_up_hgto)                    &
            = cc_wrt_bra(1)*hbra_pints(addr_cur_hgto) &
            + real(order_hgto,REALK)*half_recip_nup*hbra_pints(base_low_hgto+1)
          ! recurrence relations along y-direction
          addr_up_hgto = addr_up_hgto+1
          hbra_pints(addr_up_hgto) = cc_wrt_bra(2)*hbra_pints(addr_cur_hgto)
          do iorder = 1, order_hgto
            addr_up_hgto = addr_up_hgto+1
            addr_cur_hgto = addr_cur_hgto+1
            hbra_pints(addr_up_hgto)                    &
              = cc_wrt_bra(2)*hbra_pints(addr_cur_hgto) &
              + real(iorder,REALK)*half_recip_nup*hbra_pints(base_low_hgto+iorder)
          end do
          ! recurrence relations along z-direction
          addr_cur_hgto = base_cur_hgto
          do jorder = 0, order_hgto
            addr_up_hgto = addr_up_hgto+1
            addr_cur_hgto = addr_cur_hgto+1
            hbra_pints(addr_up_hgto) = cc_wrt_bra(3)*hbra_pints(addr_cur_hgto)
          end do
          addr_low_hgto = base_low_hgto
          do iorder = 1, order_hgto
            do jorder = 0, order_hgto-iorder
              addr_up_hgto = addr_up_hgto+1
              addr_cur_hgto = addr_cur_hgto+1
              addr_low_hgto = addr_low_hgto+1
              hbra_pints(addr_up_hgto)                    &
                = cc_wrt_bra(3)*hbra_pints(addr_cur_hgto) &
                + real(iorder,REALK)*half_recip_nup*hbra_pints(addr_low_hgto)
            end do
          end do
          ! updates base addresses
          base_low_hgto = base_cur_hgto
          base_cur_hgto = base_up_hgto
          base_up_hgto = addr_up_hgto
        end do
      end select
    ! high order HGTOs (>1) returned
    case default
      ! only order \var(orders_hgto_bra(1)) HGTOs returned
      if (orders_hgto_bra(1)==orders_hgto_bra(2)) then
        dim_tmp = (orders_hgto_bra(1)+1)*(orders_hgto_bra(1)+2)/2
        allocate(tmp_ints(3*dim_tmp), stat=ierr)
        if (ierr/=0)                                     &
          call error_stop("carmom_hbra",                 &
                          "failed to allocate tmp_ints", &
                          3*dim_tmp)
        ! computes the integral with zeroth order HGTO
        tmp_ints(1) = scal_const*(-exponent_ket-exponent_ket)**order_elec &
                    * exp(-reduced_expnt*sd_bra_ket)*sqrt((PI*recip_total_expnt)**3)
        ! computes parameter \f$-\frac{b_{j\lambda}}{p_{ij}}\f$
        half_recip_nup = -exponent_ket*recip_total_expnt
        ! computes the relative coordinates of center-of-charge w.r.t. bra center
        do iorder = 1, 3
          cc_wrt_bra(iorder) = half_recip_nup*(coord_bra(iorder)-coord_ket(iorder))
        end do
        ! computes parameter \f$-\frac{1}{2p_{ij}}\frac{b_{j\lambda}}{a_{i\kappa}}\f$
        half_recip_nup = 0.5_REALK*half_recip_nup/exponent_bra
        ! first order
        tmp_ints(dim_tmp+1) = cc_wrt_bra(1)*tmp_ints(1)
        tmp_ints(dim_tmp+2) = cc_wrt_bra(2)*tmp_ints(1)
        tmp_ints(dim_tmp+3) = cc_wrt_bra(3)*tmp_ints(1)
        ! initializes base addresses
        base_low_hgto = 0
        base_cur_hgto = dim_tmp
        base_up_hgto = 2*dim_tmp
        do order_hgto = 1, orders_hgto_bra(1)-1
          ! recurrence relations along x-direction
          addr_up_hgto = base_up_hgto+1
          addr_cur_hgto = base_cur_hgto+1
          tmp_ints(addr_up_hgto)                    &
            = cc_wrt_bra(1)*tmp_ints(addr_cur_hgto) &
            + real(order_hgto,REALK)*half_recip_nup*tmp_ints(base_low_hgto+1)
          ! recurrence relations along y-direction
          addr_up_hgto = addr_up_hgto+1
          tmp_ints(addr_up_hgto) = cc_wrt_bra(2)*tmp_ints(addr_cur_hgto)
          do iorder = 1, order_hgto
            addr_up_hgto = addr_up_hgto+1
            addr_cur_hgto = addr_cur_hgto+1
            tmp_ints(addr_up_hgto)                    &
              = cc_wrt_bra(2)*tmp_ints(addr_cur_hgto) &
              + real(iorder,REALK)*half_recip_nup*tmp_ints(base_low_hgto+iorder)
          end do
          ! recurrence relations along z-direction
          addr_cur_hgto = base_cur_hgto
          do jorder = 0, order_hgto
            addr_up_hgto = addr_up_hgto+1
            addr_cur_hgto = addr_cur_hgto+1
            tmp_ints(addr_up_hgto) = cc_wrt_bra(3)*tmp_ints(addr_cur_hgto)
          end do
          addr_low_hgto = base_low_hgto
          do iorder = 1, order_hgto
            do jorder = 0, order_hgto-iorder
              addr_up_hgto = addr_up_hgto+1
              addr_cur_hgto = addr_cur_hgto+1
              addr_low_hgto = addr_low_hgto+1
              tmp_ints(addr_up_hgto)                    &
                = cc_wrt_bra(3)*tmp_ints(addr_cur_hgto) &
                + real(iorder,REALK)*half_recip_nup*tmp_ints(addr_low_hgto)
            end do
          end do
          ! updates base addresses
          addr_up_hgto = base_low_hgto
          base_low_hgto = base_cur_hgto
          base_cur_hgto = base_up_hgto
          base_up_hgto = addr_up_hgto
        end do
        ! assigns order \var(orders_hgto_bra(1))
        hbra_pints(:) = tmp_ints(base_cur_hgto+1:base_cur_hgto+dim_tmp)
        deallocate(tmp_ints)
      else
        num_min_hgto = (orders_hgto_bra(1)+1)*(orders_hgto_bra(1)+2)/2
        dim_tmp = num_min_hgto+orders_hgto_bra(1)+2
        allocate(tmp_ints(3*dim_tmp), stat=ierr)
        if (ierr/=0)                                     &
          call error_stop("carmom_hbra",                 &
                          "failed to allocate tmp_ints", &
                          3*dim_tmp)
        ! computes the integral with zeroth order HGTO
        tmp_ints(1) = scal_const*(-exponent_ket-exponent_ket)**order_elec &
                    * exp(-reduced_expnt*sd_bra_ket)*sqrt((PI*recip_total_expnt)**3)
        ! computes parameter \f$-\frac{b_{j\lambda}}{p_{ij}}\f$
        half_recip_nup = -exponent_ket*recip_total_expnt
        ! computes the relative coordinates of center-of-charge w.r.t. bra center
        do iorder = 1, 3
          cc_wrt_bra(iorder) = half_recip_nup*(coord_bra(iorder)-coord_ket(iorder))
        end do
        ! computes parameter \f$-\frac{1}{2p_{ij}}\frac{b_{j\lambda}}{a_{i\kappa}}\f$
        half_recip_nup = 0.5_REALK*half_recip_nup/exponent_bra
        ! first order
        tmp_ints(dim_tmp+1) = cc_wrt_bra(1)*tmp_ints(1)
        tmp_ints(dim_tmp+2) = cc_wrt_bra(2)*tmp_ints(1)
        tmp_ints(dim_tmp+3) = cc_wrt_bra(3)*tmp_ints(1)
        ! initializes base addresses
        base_low_hgto = 0
        base_cur_hgto = dim_tmp
        base_up_hgto = 2*dim_tmp
        ! orders up to \var(orders_hgto_bra(1))+1
        do order_hgto = 1, orders_hgto_bra(1)
          ! recurrence relations along x-direction
          addr_up_hgto = base_up_hgto+1
          addr_cur_hgto = base_cur_hgto+1
          tmp_ints(addr_up_hgto)                    &
            = cc_wrt_bra(1)*tmp_ints(addr_cur_hgto) &
            + real(order_hgto,REALK)*half_recip_nup*tmp_ints(base_low_hgto+1)
          ! recurrence relations along y-direction
          addr_up_hgto = addr_up_hgto+1
          tmp_ints(addr_up_hgto) = cc_wrt_bra(2)*tmp_ints(addr_cur_hgto)
          do iorder = 1, order_hgto
            addr_up_hgto = addr_up_hgto+1
            addr_cur_hgto = addr_cur_hgto+1
            tmp_ints(addr_up_hgto)                    &
              = cc_wrt_bra(2)*tmp_ints(addr_cur_hgto) &
              + real(iorder,REALK)*half_recip_nup*tmp_ints(base_low_hgto+iorder)
          end do
          ! recurrence relations along z-direction
          addr_cur_hgto = base_cur_hgto
          do jorder = 0, order_hgto
            addr_up_hgto = addr_up_hgto+1
            addr_cur_hgto = addr_cur_hgto+1
            tmp_ints(addr_up_hgto) = cc_wrt_bra(3)*tmp_ints(addr_cur_hgto)
          end do
          addr_low_hgto = base_low_hgto
          do iorder = 1, order_hgto
            do jorder = 0, order_hgto-iorder
              addr_up_hgto = addr_up_hgto+1
              addr_cur_hgto = addr_cur_hgto+1
              addr_low_hgto = addr_low_hgto+1
              tmp_ints(addr_up_hgto)                    &
                = cc_wrt_bra(3)*tmp_ints(addr_cur_hgto) &
                + real(iorder,REALK)*half_recip_nup*tmp_ints(addr_low_hgto)
            end do
          end do
          ! updates base addresses
          addr_up_hgto = base_low_hgto
          base_low_hgto = base_cur_hgto
          base_cur_hgto = base_up_hgto
          base_up_hgto = addr_up_hgto
        end do
        ! assigns orders \var(orders_hgto_bra(1)) and \var(orders_hgto_bra(1))+1,
        ! and initializes base addresses
        hbra_pints(1:num_min_hgto) &
          = tmp_ints(base_low_hgto+1:base_low_hgto+num_min_hgto)
        base_up_hgto = num_min_hgto+dim_tmp
        hbra_pints(num_min_hgto+1:base_up_hgto) &
          = tmp_ints(base_cur_hgto+1:base_cur_hgto+dim_tmp)
        deallocate(tmp_ints)
        base_low_hgto = 0
        base_cur_hgto = num_min_hgto
        ! left returned orders up to \var(orders_hgto_bra(2))
        do order_hgto = orders_hgto_bra(1)+1, orders_hgto_bra(2)-1
          ! recurrence relations along x-direction
          addr_up_hgto = base_up_hgto+1
          addr_cur_hgto = base_cur_hgto+1
          hbra_pints(addr_up_hgto)                    &
            = cc_wrt_bra(1)*hbra_pints(addr_cur_hgto) &
            + real(order_hgto,REALK)*half_recip_nup*hbra_pints(base_low_hgto+1)
          ! recurrence relations along y-direction
          addr_up_hgto = addr_up_hgto+1
          hbra_pints(addr_up_hgto) = cc_wrt_bra(2)*hbra_pints(addr_cur_hgto)
          do iorder = 1, order_hgto
            addr_up_hgto = addr_up_hgto+1
            addr_cur_hgto = addr_cur_hgto+1
            hbra_pints(addr_up_hgto)                    &
              = cc_wrt_bra(2)*hbra_pints(addr_cur_hgto) &
              + real(iorder,REALK)*half_recip_nup*hbra_pints(base_low_hgto+iorder)
          end do
          ! recurrence relations along z-direction
          addr_cur_hgto = base_cur_hgto
          do jorder = 0, order_hgto
            addr_up_hgto = addr_up_hgto+1
            addr_cur_hgto = addr_cur_hgto+1
            hbra_pints(addr_up_hgto) = cc_wrt_bra(3)*hbra_pints(addr_cur_hgto)
          end do
          addr_low_hgto = base_low_hgto
          do iorder = 1, order_hgto
            do jorder = 0, order_hgto-iorder
              addr_up_hgto = addr_up_hgto+1
              addr_cur_hgto = addr_cur_hgto+1
              addr_low_hgto = addr_low_hgto+1
              hbra_pints(addr_up_hgto)                    &
                = cc_wrt_bra(3)*hbra_pints(addr_cur_hgto) &
                + real(iorder,REALK)*half_recip_nup*hbra_pints(addr_low_hgto)
            end do
          end do
          ! updates base addresses
          base_low_hgto = base_cur_hgto
          base_cur_hgto = base_up_hgto
          base_up_hgto = addr_up_hgto
        end do
      end if
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "carmom_hbra", STDOUT)
#endif
    return
  end subroutine carmom_hbra
