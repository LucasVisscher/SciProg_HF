!!  gen1int: compute derivatives of one-electron integrals using Hermite Gaussians
!!  Copyright 2009-2012 Bin Gao, and Andreas Thorvaldsen
!!
!!  This file is part of gen1int.
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
!!  This file contains subroutines related to partial geometric derivatives
!!  of integrals of zero-center operators.
!!
!!  2011-07-02, Bin Gao:
!!  * first version

#include "stdout.h"
#include "private/tag_cent.h"

  !> \brief sets the orders of partial derivative terms for operators without any center
  !> \author Bin Gao
  !> \date 2011-07-02
  !> \param num_cents is the number of differentiated centers, could be 1 or 2
  !> \param idx_cent contains the indices of different differentiated centers
  !> \param order_cent contains the order of total geometric derivatives of differentiated centers
  !> \param idx_part_cent contains the indices of bra and ket centers
  !> \return order_part_cent contains the orders of geometric derivatives with respect to
  !>         bra and ket centers on entry, while the orders of resulted partial geometric
  !>         derivatives are added on exit
  !> \return zero_ints indicates if the total geometric derivatives are zero
  !> \return scatter_deriv indicates if scattering the geometric derivatives later on
  !> \return seq_part_geo contains the sequence of bra and ket centers for partial derivative terms
  subroutine geom_part_zero_param(num_cents, idx_cent, order_cent, &
                                  idx_part_cent, order_part_cent,  &
                                  zero_ints, scatter_deriv, seq_part_geo)
    use xkind
    implicit none
    integer, intent(in) :: num_cents
    integer, intent(in) :: idx_cent(num_cents)
    integer, intent(in) :: order_cent(num_cents)
    integer, intent(in) :: idx_part_cent(2)
    integer, intent(inout) :: order_part_cent(2)
    logical, intent(out) :: zero_ints
    logical, intent(out) :: scatter_deriv
    integer, intent(out) :: seq_part_geo(num_cents)
!f2py intent(hide) :: num_cents
!f2py intent(in) :: idx_cent
!f2py intent(in) :: order_cent
!f2py depend(num_cents) :: order_cent
!f2py intent(in) :: idx_part_cent
!f2py intent(inout) :: order_part_cent
!f2py intent(out) :: zero_ints
!f2py intent(out) :: scatter_deriv
!f2py intent(out) :: seq_part_geo
!f2py depend(num_cents) :: seq_part_geo
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! centers of bra and ket are identical
    if (idx_part_cent(1)==idx_part_cent(2)) then
      zero_ints = .true.
    else
      if (num_cents>0) then
        select case(num_cents)
        ! one-center geometric derivatives
        case(1)
          if (idx_cent(1)==idx_part_cent(1)) then
            zero_ints = .false.
            ! needs to scatter derivatives when there are partial geometric derivatives
            ! either on bra or ket center
            scatter_deriv = any(order_part_cent/=0)
            order_part_cent(1) = order_part_cent(1)+order_cent(1)
            seq_part_geo(1) = 1
          else if (idx_cent(1)==idx_part_cent(2)) then
            zero_ints = .false.
            ! needs to scatter derivatives when there are partial geometric derivatives on ket center
            scatter_deriv = order_part_cent(2)/=0
            order_part_cent(2) = order_part_cent(2)+order_cent(1)
            seq_part_geo(1) = 2
          else
            zero_ints = .true.
          end if
        ! two-center geometric derivatives
        case(2)
          if (all(idx_cent==idx_part_cent)) then
            zero_ints = .false.
            ! needs to scatter derivatives when there are partial geometric derivatives
            ! either on bra or ket center
            scatter_deriv = any(order_part_cent/=0)
            order_part_cent = order_part_cent+order_cent
            seq_part_geo(1) = 1
            seq_part_geo(2) = 2
          else if (idx_cent(1)==idx_part_cent(2) .and. idx_cent(2)==idx_part_cent(1)) then
            zero_ints = .false.
            ! always needs to scatter derivatives
            scatter_deriv = .true.
            order_part_cent(1) = order_part_cent(1)+order_cent(2)
            order_part_cent(2) = order_part_cent(2)+order_cent(1)
            seq_part_geo(1) = 2
            seq_part_geo(2) = 1
          else
            zero_ints = .true.
          end if
        ! too many differentiated centers
        case default
          zero_ints = .true.
        end select
      else
        call error_stop("geom_part_zero_param", "invalid num_cents", num_cents)
      end if
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "geom_part_zero_param", STDOUT)
#endif
    return
  end subroutine geom_part_zero_param

  !> \brief scatters geometric derivatives for operators without any center
  !> \author Bin Gao
  !> \date 2011-07-02
  !> \param num_cents is the number of differentiated centers, could be 1 or 2
  !> \param order_cent contains the order of total geometric derivatives of differentiated centers
  !> \param seq_part_geo contains the sequence of bra and ket centers for
  !>        partial derivative terms, from \fn(geom_part_zero_param)
  !> \param order_geo_bra is the order of partial geometric derivatives on bra center
  !> \param order_geo_ket is the order of partial geometric derivatives on ket center
  !> \param dim_cints is the dimension of contracted integrals
  !> \param dim_geo_bra is the dimension of geometric derivatives on bra center before scattering
  !> \param dim_geo_ket is the dimension of geometric derivatives on ket center before scattering
  !> \param num_opt is the number of other operators
  !> \param part_cints contains the derivatives of integrals before scattering
  !> \param num_geo_bra is the number of partial geometric derivatives on bra center after scattering,
  !>        should equal to \f$(\var(order_geo_bra)+1)(\var(order_geo_bra)+2)/2\f$
  !> \param num_geo_ket is the number of partial geometric derivatives on ket center after scattering,
  !>        should equal to \f$(\var(order_geo_ket)+1)(\var(order_geo_ket)+2)/2\f$
  !> \param num_tot_geo is the number of total geometric derivatives after scattering
  !> \return total_cints contains the derivatives of integrals after scattering
  subroutine geom_part_zero_scatter(num_cents, order_cent, seq_part_geo,     &
                                    order_geo_bra, order_geo_ket, dim_cints, &
                                    dim_geo_bra, dim_geo_ket, num_opt,       &
                                    part_cints, num_geo_bra, num_geo_ket,    &
                                    num_tot_geo, total_cints)
    use xkind
    implicit none
    integer, intent(in) :: num_cents
    integer, intent(in) :: order_cent(num_cents)
    integer, intent(in) :: seq_part_geo(num_cents)
    integer, intent(in) :: order_geo_bra
    integer, intent(in) :: order_geo_ket
    integer, intent(in) :: dim_cints
    integer, intent(in) :: dim_geo_bra
    integer, intent(in) :: dim_geo_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(in) :: part_cints(dim_cints,dim_geo_bra,dim_geo_ket,num_opt)
    integer, intent(in) :: num_geo_bra
    integer, intent(in) :: num_geo_ket
    integer, intent(in) :: num_tot_geo
    real(REALK), intent(out) :: total_cints(dim_cints,num_geo_bra,num_geo_ket, &
                                            num_tot_geo,num_opt)
!f2py intent(hide) :: num_cents
!f2py intent(in) :: order_cent
!f2py intent(in) :: seq_part_geo
!f2py depend(num_cents) :: seq_part_geo
!f2py intent(in) :: order_geo_bra
!f2py intent(in) :: order_geo_ket
!f2py intent(hide) :: dim_cints
!f2py intent(hide) :: dim_geo_bra
!f2py intent(hide) :: dim_geo_ket
!f2py intent(hide) :: num_opt
!f2py intent(in) :: part_cints
!f2py intent(in) :: num_geo_bra
!f2py intent(in) :: num_geo_ket
!f2py intent(in) :: num_tot_geo
!f2py intent(out) :: total_cints
!f2py depend(dim_cints) :: total_cints
!f2py depend(num_geo_bra) :: total_cints
!f2py depend(num_geo_ket) :: total_cints
!f2py depend(num_tot_geo) :: total_cints
!f2py depend(num_opt) :: total_cints
    real(REALK), allocatable :: scatter_cints(:,:,:,:,:)
                              !integrals after scattering derivatives on bra or ket
    integer num_tot_bra       !number of total geometric derivatives on bra center
    integer num_tot_ket       !number of total geometric derivatives on ket center
    integer nopt              !incremental recorder over other operators
    integer ibra, iket, igeo  !incremental recorder over geometric derivatives
    integer ierr              !error information
#if defined(XTIME)
    real(REALK) curr_time     !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    select case(num_cents)
    ! one-center total geometric derivatives
    case(1)
      ! geometric derivatives on bra center
      if (seq_part_geo(1)==TAG_BRA) then
        ! \var(order_geo_ket)/=0
        if (order_geo_bra==0) then
          do nopt = 1, num_opt
            do ibra = 1, dim_geo_bra
              total_cints(:,1,:,ibra,nopt) = part_cints(:,ibra,:,nopt)
            end do
          end do
        else
          call scatter_single(dim_cints, dim_geo_bra, dim_geo_ket, num_opt, &
                              part_cints, order_geo_bra, order_cent(1),     &
                              num_geo_bra, num_tot_geo, total_cints)
        end if
      ! geometric derivatives on ket center
      else
        ! notice that \var(order_geo_ket)/=0
        call scatter_single(dim_cints*dim_geo_bra, dim_geo_ket, 1, num_opt, &
                            part_cints, order_geo_ket, order_cent(1),       &
                            num_geo_ket, num_tot_geo, total_cints)
      end if
    ! two-center total geometric derivatives
    case(2)
      ! the first differentiated center is bra center
      if (seq_part_geo(1)==TAG_BRA) then
        ! no partial geometric derivatives on ket center
        if (order_geo_ket==0) then
          ! sets the number of total geometric derivatives on bra center
          num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
          call scatter_single(dim_cints, dim_geo_bra, 1, dim_geo_ket*num_opt, &
                              part_cints, order_geo_bra, order_cent(1),       &
                              num_geo_bra, num_tot_bra, total_cints)
        else
          ! sets the number of total geometric derivatives on ket center
          num_tot_ket = (order_cent(2)+1)*(order_cent(2)+2)/2
          ! allocates the integrals after scattering geometric derivatives on ket
          allocate(scatter_cints(dim_cints,dim_geo_bra,num_geo_ket, &
                                 num_tot_ket,num_opt), stat=ierr)
          if (ierr/=0)                                               &
            call error_stop("geom_part_zero_scatter",                &
                            "failed to allocate scatter_cints/1/?K", &
                            dim_cints*dim_geo_bra*num_geo_ket*num_tot_ket*num_opt)
          ! scatters geometric derivatives on ket center
          call scatter_single(dim_cints*dim_geo_bra, dim_geo_ket, 1, num_opt, &
                              part_cints, order_geo_ket, order_cent(2),       &
                              num_geo_ket, num_tot_ket, scatter_cints)
          ! no partial geometric derivatives on bra center
          if (order_geo_bra==0) then
            ! scatters geometric derivatives on bra center
            do nopt = 1, num_opt
              igeo = 0
              do iket = 1, num_tot_ket
                do ibra = 1, dim_geo_bra
                  igeo = igeo+1
                  total_cints(:,1,:,igeo,nopt) = scatter_cints(:,ibra,:,iket,nopt)
                end do
              end do
            end do
          else
            ! sets the number of total geometric derivatives on bra center
            num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
            ! scatters geometric derivatives on bra center
            call scatter_single(dim_cints, dim_geo_bra, num_geo_ket, &
                                num_tot_ket*num_opt, scatter_cints,  &
                                order_geo_bra, order_cent(1),        &
                                num_geo_bra, num_tot_bra, total_cints)
          end if
          deallocate(scatter_cints)
        end if
      ! the first differentiated center is ket center
      else
        ! no partial geometric derivatives on bra center
        if (order_geo_bra==0) then
          if (order_geo_ket==0) then
            ! switches the derivatives on bra and ket centers
            do nopt = 1, num_opt
              igeo = 0
              do ibra = 1, dim_geo_bra
                do iket = 1, dim_geo_ket
                  igeo = igeo+1
                  total_cints(:,1,1,igeo,nopt) = part_cints(:,ibra,iket,nopt)
                end do
              end do
            end do
          else
            ! allocates the integrals after scattering geometric derivatives on bra
            allocate(scatter_cints(dim_cints,1,dim_geo_ket,dim_geo_bra,num_opt), stat=ierr)
            if (ierr/=0)                                               &
              call error_stop("geom_part_zero_scatter",                &
                              "failed to allocate scatter_cints/2/K0", &
                              dim_cints*dim_geo_ket*dim_geo_bra*num_opt)
            ! scatters geometric derivatives on bra center
            do nopt = 1, num_opt
              do ibra = 1, dim_geo_bra
                scatter_cints(:,1,:,ibra,nopt) = part_cints(:,ibra,:,nopt)
              end do
            end do
            ! sets the number of total geometric derivatives on ket center
            num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
            ! scatters geometric derivatives on ket center
            call scatter_single(dim_cints, dim_geo_ket, 1, dim_geo_bra*num_opt, &
                                scatter_cints, order_geo_ket, order_cent(1),    &
                                num_geo_ket, num_tot_ket, total_cints)
            deallocate(scatter_cints)
          end if
        ! no partial geometric derivatives on ket center
        else if (order_geo_ket==0) then
          ! sets the number of total geometric derivatives on bra center
          num_tot_bra = (order_cent(2)+1)*(order_cent(2)+2)/2
          call scatter_single(dim_cints, dim_geo_bra, dim_geo_ket, num_opt, &
                              part_cints, order_geo_bra, order_cent(2),     &
                              num_geo_bra, num_tot_bra, total_cints)
        else
          ! sets the number of total geometric derivatives on bra center
          num_tot_bra = (order_cent(2)+1)*(order_cent(2)+2)/2
          ! allocates the integrals after scattering geometric derivatives on bra
          allocate(scatter_cints(dim_cints,num_geo_bra,dim_geo_ket, &
                                 num_tot_bra,num_opt), stat=ierr)
          if (ierr/=0)                                               &
            call error_stop("geom_part_zero_scatter",                &
                            "failed to allocate scatter_cints/2/KB", &
                            dim_cints*num_geo_bra*dim_geo_ket*num_tot_bra*num_opt)
          ! scatters geometric derivatives on bra center
          call scatter_single(dim_cints, dim_geo_bra, dim_geo_ket, num_opt, &
                              part_cints, order_geo_bra, order_cent(2),     &
                              num_geo_bra, num_tot_bra, scatter_cints)
          ! sets the number of total geometric derivatives on ket center
          num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
          ! scatters geometric derivatives on ket center
          call scatter_single(dim_cints*num_geo_bra, dim_geo_ket, 1, &
                              num_tot_bra*num_opt, scatter_cints,    &
                              order_geo_ket, order_cent(1),          &
                              num_geo_ket, num_tot_ket, total_cints)
          deallocate(scatter_cints)
        end if
      end if
    case default
      call error_stop("geom_part_zero_scatter", "invalid num_cents", num_cents)
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "geom_part_zero_scatter", STDOUT)
#endif
    return
  end subroutine geom_part_zero_scatter
