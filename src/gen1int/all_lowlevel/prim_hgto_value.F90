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
!!  This file contains the recurrence relations of primitive Hermite Gaussians.
!!
!!  2012-03-09, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief recurrence relations of primitive Hermite Gaussians
  !> \author Bin Gao
  !> \date 2012-03-09
  !> \param orders_hgto_bra is the range of orders of Hermite Gaussians on bra center
  !> \param coord_bra is the coordinates of bra center
  !> \param exponent_bra is the exponent of primitive Gaussian of bra center
  !> \param num_points is the number of grid points
  !> \param grid_points contains the coordinates of grid points
  !> \param dim_hgto_bra is the number of Hermite Gaussians on bra center
  !> \return hgto_value contains the values of primitive Hermite Gaussians
  subroutine prim_hgto_value(orders_hgto_bra, coord_bra, exponent_bra, &
                             num_points, grid_points, dim_hgto_bra,    &
                             hgto_value)
    use xkind
    implicit none
    integer, intent(in) :: orders_hgto_bra(2)
    real(REALK), intent(in) :: coord_bra(3)
    real(REALK), intent(in) :: exponent_bra
    integer, intent(in) :: num_points
    real(REALK), intent(in) :: grid_points(3,num_points)
    integer, intent(in) :: dim_hgto_bra
    real(REALK), intent(out) :: hgto_value(dim_hgto_bra,num_points)
!f2py intent(in) :: orders_hgto_bra
!f2py intent(in) :: coord_bra
!f2py intent(in) :: exponent_bra
!f2py intent(hide) :: num_points
!f2py intent(in) :: grid_points
!f2py intent(in) :: dim_hgto_bra
!f2py intent(out) :: hgto_value
!f2py depend(dim_hgto_bra) :: hgto_value
!f2py depend(num_points) :: hgto_value
    real(REALK) half_nr_expnt                   !half of the negative reciprocal of exponent
    integer ipoint                              !incremental recorder over grid points
    real(REALK) elec_coord(3)                   !electron coordinates
    integer base_low_hgto                       !base address of lower order HGTOs
    integer base_cur_hgto                       !base address of current order HGTOs
    integer base_up_hgto                        !base address of upper order HGTOs
    integer order_hgto, iorder, jorder          !incremental recorders over orders of HGTOs
    integer addr_low_hgto                       !address of lower order HGTOs
    integer addr_cur_hgto                       !address of current order HGTOs
    integer addr_up_hgto                        !address of upper order HGTOs
    integer num_min_hgto                        !number of minimum returned order HGTOs
    integer dim_tmp_hgto                        !dimension of temporary values of HGTOs
    real(REALK), allocatable :: tmp_hgto(:,:)   !temporary values of HGTOs
    real(REALK) zero_hgto                       !value of zeroth order HGTO
    integer ierr                                !error information
#if defined(XTIME)
    real(REALK) curr_time                       !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! sets the half of the negative reciprocal of exponent on bra center
    half_nr_expnt = -0.5_REALK/exponent_bra
    ! recovers the Hermite Gaussians on bra center
    select case(orders_hgto_bra(1))
    ! zeroth order HGTOs returned
    case(0)
      select case(orders_hgto_bra(2))
      ! only zeroth order HGTOs returned
      case(0)
        do ipoint = 1, num_points
          elec_coord = grid_points(:,ipoint)-coord_bra
          hgto_value(1,ipoint) = &
            exp(-exponent_bra*(elec_coord(1)**2+elec_coord(2)**2+elec_coord(3)**2))
        end do
      ! up to, at least, first order HGTOs returned
      case default
        do ipoint = 1, num_points
          elec_coord = grid_points(:,ipoint)-coord_bra
          ! zeroth order
          hgto_value(1,ipoint) = &
            exp(-exponent_bra*(elec_coord(1)**2+elec_coord(2)**2+elec_coord(3)**2))
          ! first order
          hgto_value(2,ipoint) = elec_coord(1)*hgto_value(1,ipoint)
          hgto_value(3,ipoint) = elec_coord(2)*hgto_value(1,ipoint)
          hgto_value(4,ipoint) = elec_coord(3)*hgto_value(1,ipoint)
          ! initializes base addresses
          base_low_hgto = 0
          base_cur_hgto = 1
          base_up_hgto = 4
          ! left returned orders up to \var(orders_hgto_bra(2))
          do order_hgto = 1, orders_hgto_bra(2)-1
            ! recurrence relations along x-direction
            addr_up_hgto = base_up_hgto+1
            addr_cur_hgto = base_cur_hgto+1
            hgto_value(addr_up_hgto,ipoint)                    &
              = elec_coord(1)*hgto_value(addr_cur_hgto,ipoint) &
              + real(order_hgto,REALK)*half_nr_expnt           &
              * hgto_value(base_low_hgto+1,ipoint)
            ! recurrence relations along y-direction
            addr_up_hgto = addr_up_hgto+1
            hgto_value(addr_up_hgto,ipoint) &
              = elec_coord(2)*hgto_value(addr_cur_hgto,ipoint)
            do iorder = 1, order_hgto
              addr_up_hgto = addr_up_hgto+1
              addr_cur_hgto = addr_cur_hgto+1
              hgto_value(addr_up_hgto,ipoint)                    &
                = elec_coord(2)*hgto_value(addr_cur_hgto,ipoint) &
                + real(iorder,REALK)*half_nr_expnt               &
                * hgto_value(base_low_hgto+iorder,ipoint)
            end do
            ! recurrence relations along z-direction
            addr_cur_hgto = base_cur_hgto
            do jorder = 0, order_hgto
              addr_up_hgto = addr_up_hgto+1
              addr_cur_hgto = addr_cur_hgto+1
              hgto_value(addr_up_hgto,ipoint) &
                = elec_coord(3)*hgto_value(addr_cur_hgto,ipoint)
            end do
            addr_low_hgto = base_low_hgto
            do iorder = 1, order_hgto
              do jorder = 0, order_hgto-iorder
                addr_up_hgto = addr_up_hgto+1
                addr_cur_hgto = addr_cur_hgto+1
                addr_low_hgto = addr_low_hgto+1
                hgto_value(addr_up_hgto,ipoint)                    &
                  = elec_coord(3)*hgto_value(addr_cur_hgto,ipoint) &
                  + real(iorder,REALK)*half_nr_expnt               &
                  * hgto_value(addr_low_hgto,ipoint)
              end do
            end do
            ! updates base addresses
            base_low_hgto = base_cur_hgto
            base_cur_hgto = base_up_hgto
            base_up_hgto = addr_up_hgto
          end do
        end do
      end select
    ! first order HGTOs returned
    case (1)
      select case(orders_hgto_bra(2))
      ! only first order HGTOs returned
      case(1)
        do ipoint = 1, num_points
          elec_coord = grid_points(:,ipoint)-coord_bra
          ! zeroth order
          hgto_value(1,ipoint) = &
            exp(-exponent_bra*(elec_coord(1)**2+elec_coord(2)**2+elec_coord(3)**2))
          ! first order
          hgto_value(2,ipoint) = elec_coord(2)*hgto_value(1,ipoint)
          hgto_value(3,ipoint) = elec_coord(3)*hgto_value(1,ipoint)
          hgto_value(1,ipoint) = elec_coord(1)*hgto_value(1,ipoint)
        end do
      ! up to, at least, second order HGTOs returned
      case default
        do ipoint = 1, num_points
          elec_coord = grid_points(:,ipoint)-coord_bra
          ! zeroth order
          zero_hgto = &
            exp(-exponent_bra*(elec_coord(1)**2+elec_coord(2)**2+elec_coord(3)**2))
          ! first order
          hgto_value(1,ipoint) = elec_coord(1)*zero_hgto
          hgto_value(2,ipoint) = elec_coord(2)*zero_hgto
          hgto_value(3,ipoint) = elec_coord(3)*zero_hgto
          ! second order
          hgto_value(4,ipoint) = elec_coord(1)*hgto_value(1,ipoint) &
                               + half_nr_expnt*zero_hgto
          hgto_value(5,ipoint) = elec_coord(2)*hgto_value(1,ipoint)
          hgto_value(6,ipoint) = elec_coord(2)*hgto_value(2,ipoint) &
                               + half_nr_expnt*zero_hgto
          hgto_value(7,ipoint) = elec_coord(3)*hgto_value(1,ipoint) 
          hgto_value(8,ipoint) = elec_coord(3)*hgto_value(2,ipoint) 
          hgto_value(9,ipoint) = elec_coord(3)*hgto_value(3,ipoint) &
                               + half_nr_expnt*zero_hgto
          ! initializes base addresses
          base_low_hgto = 0
          base_cur_hgto = 3
          base_up_hgto = 9
          ! left returned orders up to \var(orders_hgto_bra(2))
          do order_hgto = 2, orders_hgto_bra(2)-1
            ! recurrence relations along x-direction
            addr_up_hgto = base_up_hgto+1
            addr_cur_hgto = base_cur_hgto+1
            hgto_value(addr_up_hgto,ipoint)                    &
              = elec_coord(1)*hgto_value(addr_cur_hgto,ipoint) &
              + real(order_hgto,REALK)*half_nr_expnt           &
              * hgto_value(base_low_hgto+1,ipoint)
            ! recurrence relations along y-direction
            addr_up_hgto = addr_up_hgto+1
            hgto_value(addr_up_hgto,ipoint) &
              = elec_coord(2)*hgto_value(addr_cur_hgto,ipoint)
            do iorder = 1, order_hgto
              addr_up_hgto = addr_up_hgto+1
              addr_cur_hgto = addr_cur_hgto+1
              hgto_value(addr_up_hgto,ipoint)                    &
                = elec_coord(2)*hgto_value(addr_cur_hgto,ipoint) &
                + real(iorder,REALK)*half_nr_expnt               &
                * hgto_value(base_low_hgto+iorder,ipoint)
            end do
            ! recurrence relations along z-direction
            addr_cur_hgto = base_cur_hgto
            do jorder = 0, order_hgto
              addr_up_hgto = addr_up_hgto+1
              addr_cur_hgto = addr_cur_hgto+1
              hgto_value(addr_up_hgto,ipoint) &
                = elec_coord(3)*hgto_value(addr_cur_hgto,ipoint)
            end do
            addr_low_hgto = base_low_hgto
            do iorder = 1, order_hgto
              do jorder = 0, order_hgto-iorder
                addr_up_hgto = addr_up_hgto+1
                addr_cur_hgto = addr_cur_hgto+1
                addr_low_hgto = addr_low_hgto+1
                hgto_value(addr_up_hgto,ipoint)                    &
                  = elec_coord(3)*hgto_value(addr_cur_hgto,ipoint) &
                  + real(iorder,REALK)*half_nr_expnt               &
                  * hgto_value(addr_low_hgto,ipoint)
              end do
            end do
            ! updates base addresses
            base_low_hgto = base_cur_hgto
            base_cur_hgto = base_up_hgto
            base_up_hgto = addr_up_hgto
          end do
        end do
      end select
    ! high order HGTOs (>1) returned
    case default
      ! only order \var(orders_hgto_bra(1)) HGTOs returned
      if (orders_hgto_bra(1)==orders_hgto_bra(2)) then
        dim_tmp_hgto = (orders_hgto_bra(1)+1)*(orders_hgto_bra(1)+2)/2
        allocate(tmp_hgto(3*dim_tmp_hgto,1), stat=ierr)
        if (ierr/=0)                                     &
          call error_stop("prim_hgto_value",             &
                          "failed to allocate tmp_hgto", &
                          3*dim_tmp_hgto)
        do ipoint = 1, num_points
          elec_coord = grid_points(:,ipoint)-coord_bra
          ! zeroth order
          tmp_hgto(1,1) = &
            exp(-exponent_bra*(elec_coord(1)**2+elec_coord(2)**2+elec_coord(3)**2))
          ! first order
          tmp_hgto(dim_tmp_hgto+1,1) = elec_coord(1)*tmp_hgto(1,1)
          tmp_hgto(dim_tmp_hgto+2,1) = elec_coord(2)*tmp_hgto(1,1)
          tmp_hgto(dim_tmp_hgto+3,1) = elec_coord(3)*tmp_hgto(1,1)
          ! initializes base addresses
          base_low_hgto = 0
          base_cur_hgto = dim_tmp_hgto
          base_up_hgto = 2*dim_tmp_hgto
          do order_hgto = 1, orders_hgto_bra(1)-1
            ! recurrence relations along x-direction
            addr_up_hgto = base_up_hgto+1
            addr_cur_hgto = base_cur_hgto+1
            tmp_hgto(addr_up_hgto,1)                    &
              = elec_coord(1)*tmp_hgto(addr_cur_hgto,1) &
              + real(order_hgto,REALK)*half_nr_expnt*tmp_hgto(base_low_hgto+1,1)
            ! recurrence relations along y-direction
            addr_up_hgto = addr_up_hgto+1
            tmp_hgto(addr_up_hgto,1) = elec_coord(2)*tmp_hgto(addr_cur_hgto,1)
            do iorder = 1, order_hgto
              addr_up_hgto = addr_up_hgto+1
              addr_cur_hgto = addr_cur_hgto+1
              tmp_hgto(addr_up_hgto,1)                    &
                = elec_coord(2)*tmp_hgto(addr_cur_hgto,1) &
                + real(iorder,REALK)*half_nr_expnt*tmp_hgto(base_low_hgto+iorder,1)
            end do
            ! recurrence relations along z-direction
            addr_cur_hgto = base_cur_hgto
            do jorder = 0, order_hgto
              addr_up_hgto = addr_up_hgto+1
              addr_cur_hgto = addr_cur_hgto+1
              tmp_hgto(addr_up_hgto,1) = elec_coord(3)*tmp_hgto(addr_cur_hgto,1)
            end do
            addr_low_hgto = base_low_hgto
            do iorder = 1, order_hgto
              do jorder = 0, order_hgto-iorder
                addr_up_hgto = addr_up_hgto+1
                addr_cur_hgto = addr_cur_hgto+1
                addr_low_hgto = addr_low_hgto+1
                tmp_hgto(addr_up_hgto,1)                    &
                  = elec_coord(3)*tmp_hgto(addr_cur_hgto,1) &
                  + real(iorder,REALK)*half_nr_expnt*tmp_hgto(addr_low_hgto,1)
              end do
            end do
            ! updates base addresses
            addr_up_hgto = base_low_hgto
            base_low_hgto = base_cur_hgto
            base_cur_hgto = base_up_hgto
            base_up_hgto = addr_up_hgto
          end do
          ! assigns order \var(orders_hgto_bra(1))
          hgto_value(:,ipoint) &
            = tmp_hgto(base_cur_hgto+1:base_cur_hgto+dim_tmp_hgto,1)
        end do
      else
        num_min_hgto = (orders_hgto_bra(1)+1)*(orders_hgto_bra(1)+2)/2
        dim_tmp_hgto = num_min_hgto+orders_hgto_bra(1)+2
        allocate(tmp_hgto(3*dim_tmp_hgto,1), stat=ierr)
        if (ierr/=0)                                     &
          call error_stop("prim_hgto_value",             &
                          "failed to allocate tmp_hgto", &
                          3*dim_tmp_hgto)
        do ipoint = 1, num_points
          elec_coord = grid_points(:,ipoint)-coord_bra
          ! zeroth order
          tmp_hgto(1,1) = &
            exp(-exponent_bra*(elec_coord(1)**2+elec_coord(2)**2+elec_coord(3)**2))
          ! first order
          tmp_hgto(dim_tmp_hgto+1,1) = elec_coord(1)*tmp_hgto(1,1)
          tmp_hgto(dim_tmp_hgto+2,1) = elec_coord(2)*tmp_hgto(1,1)
          tmp_hgto(dim_tmp_hgto+3,1) = elec_coord(3)*tmp_hgto(1,1)
          ! initializes base addresses
          base_low_hgto = 0
          base_cur_hgto = dim_tmp_hgto
          base_up_hgto = 2*dim_tmp_hgto
          ! orders up to \var(orders_hgto_bra(1))+1
          do order_hgto = 1, orders_hgto_bra(1)
            ! recurrence relations along x-direction
            addr_up_hgto = base_up_hgto+1
            addr_cur_hgto = base_cur_hgto+1
            tmp_hgto(addr_up_hgto,1)                    &
              = elec_coord(1)*tmp_hgto(addr_cur_hgto,1) &
              + real(order_hgto,REALK)*half_nr_expnt*tmp_hgto(base_low_hgto+1,1)
            ! recurrence relations along y-direction
            addr_up_hgto = addr_up_hgto+1
            tmp_hgto(addr_up_hgto,1) = elec_coord(2)*tmp_hgto(addr_cur_hgto,1)
            do iorder = 1, order_hgto
              addr_up_hgto = addr_up_hgto+1
              addr_cur_hgto = addr_cur_hgto+1
              tmp_hgto(addr_up_hgto,1)                    &
                = elec_coord(2)*tmp_hgto(addr_cur_hgto,1) &
                + real(iorder,REALK)*half_nr_expnt*tmp_hgto(base_low_hgto+iorder,1)
            end do
            ! recurrence relations along z-direction
            addr_cur_hgto = base_cur_hgto
            do jorder = 0, order_hgto
              addr_up_hgto = addr_up_hgto+1
              addr_cur_hgto = addr_cur_hgto+1
              tmp_hgto(addr_up_hgto,1) = elec_coord(3)*tmp_hgto(addr_cur_hgto,1)
            end do
            addr_low_hgto = base_low_hgto
            do iorder = 1, order_hgto
              do jorder = 0, order_hgto-iorder
                addr_up_hgto = addr_up_hgto+1
                addr_cur_hgto = addr_cur_hgto+1
                addr_low_hgto = addr_low_hgto+1
                tmp_hgto(addr_up_hgto,1)                    &
                  = elec_coord(3)*tmp_hgto(addr_cur_hgto,1) &
                  + real(iorder,REALK)*half_nr_expnt*tmp_hgto(addr_low_hgto,1)
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
          hgto_value(1:num_min_hgto,ipoint) &
            = tmp_hgto(base_low_hgto+1:base_low_hgto+num_min_hgto,1)
          base_up_hgto = num_min_hgto+dim_tmp_hgto
          hgto_value(num_min_hgto+1:base_up_hgto,ipoint) &
            = tmp_hgto(base_cur_hgto+1:base_cur_hgto+dim_tmp_hgto,1)
          base_low_hgto = 0
          base_cur_hgto = num_min_hgto
          ! left returned orders up to \var(orders_hgto_bra(2))
          do order_hgto = orders_hgto_bra(1)+1, orders_hgto_bra(2)-1
            ! recurrence relations along x-direction
            addr_up_hgto = base_up_hgto+1
            addr_cur_hgto = base_cur_hgto+1
            hgto_value(addr_up_hgto,ipoint)                    &
              = elec_coord(1)*hgto_value(addr_cur_hgto,ipoint) &
              + real(order_hgto,REALK)*half_nr_expnt           &
              * hgto_value(base_low_hgto+1,ipoint)
            ! recurrence relations along y-direction
            addr_up_hgto = addr_up_hgto+1
            hgto_value(addr_up_hgto,ipoint) &
              = elec_coord(2)*hgto_value(addr_cur_hgto,ipoint)
            do iorder = 1, order_hgto
              addr_up_hgto = addr_up_hgto+1
              addr_cur_hgto = addr_cur_hgto+1
              hgto_value(addr_up_hgto,ipoint)                    &
                = elec_coord(2)*hgto_value(addr_cur_hgto,ipoint) &
                + real(iorder,REALK)*half_nr_expnt               &
                * hgto_value(base_low_hgto+iorder,ipoint)
            end do
            ! recurrence relations along z-direction
            addr_cur_hgto = base_cur_hgto
            do jorder = 0, order_hgto
              addr_up_hgto = addr_up_hgto+1
              addr_cur_hgto = addr_cur_hgto+1
              hgto_value(addr_up_hgto,ipoint) &
                = elec_coord(3)*hgto_value(addr_cur_hgto,ipoint)
            end do
            addr_low_hgto = base_low_hgto
            do iorder = 1, order_hgto
              do jorder = 0, order_hgto-iorder
                addr_up_hgto = addr_up_hgto+1
                addr_cur_hgto = addr_cur_hgto+1
                addr_low_hgto = addr_low_hgto+1
                hgto_value(addr_up_hgto,ipoint)                    &
                  = elec_coord(3)*hgto_value(addr_cur_hgto,ipoint) &
                  + real(iorder,REALK)*half_nr_expnt               &
                  * hgto_value(addr_low_hgto,ipoint)
              end do
            end do
            ! updates base addresses
            base_low_hgto = base_cur_hgto
            base_cur_hgto = base_up_hgto
            base_up_hgto = addr_up_hgto
          end do
        end do
      end if
      deallocate(tmp_hgto)
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "prim_hgto_value", STDOUT)
#endif
    return
  end subroutine prim_hgto_value
