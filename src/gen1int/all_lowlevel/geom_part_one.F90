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
!!  of integrals of one-center operator.
!!
!!  2011-07-02, Bin Gao:
!!  * first version

#include "stdout.h"
#include "private/tag_cent.h"

  !> \brief sets the orders of partial derivative terms and the parameters
  !>        for operators with one center
  !> \author Bin Gao
  !> \date 2011-07-02
  !> \param num_cents is the number of differentiated centers, could be 1, 2 or 3
  !> \param idx_cent contains the indices of different differentiated centers,
  !>        should be in ascending order
  !> \param order_cent contains the order of total geometric derivatives of differentiated centers
  !> \param idx_part_cent contains the indices of bra, ket and operator centers; it will
  !>        be sorted in ascending order on exit
  !> \return order_part_cent contains the orders of geometric derivatives with respect to
  !>         bra, ket and operator centers on entry, while the orders of resulted partial
  !>         geometric derivatives are added on exit
  !> \return zero_ints indicates if the total geometric derivatives are zero
  !> \return neg_one indicates if the integrals will be multiplied by -1
  !> \return scatter_deriv indicates if scattering the geometric derivatives later on
  !> \return seq_part_geo contains the sequence of bra, ket and operator centers
  !>         for partial derivative terms
  subroutine geom_part_one_param(num_cents, idx_cent, order_cent,   &
                                 idx_part_cent, order_part_cent,    &
                                 zero_ints, neg_one, scatter_deriv, &
                                 seq_part_geo)
    use xkind
    implicit none
    integer, intent(in) :: num_cents
    integer, intent(in) :: idx_cent(num_cents)
    integer, intent(in) :: order_cent(num_cents)
    integer, intent(inout) :: idx_part_cent(3)
    integer, intent(inout) :: order_part_cent(3)
    logical, intent(out) :: zero_ints
    logical, intent(out) :: neg_one
    logical, intent(out) :: scatter_deriv
    integer, intent(out) :: seq_part_geo(num_cents)
!f2py intent(hide) :: num_cents
!f2py intent(in) :: idx_cent
!f2py intent(in) :: order_cent
!f2py depend(num_cents) :: order_cent
!f2py intent(inout) :: idx_part_cent
!f2py intent(inout) :: order_part_cent
!f2py intent(out) :: zero_ints
!f2py intent(out) :: neg_one
!f2py intent(out) :: scatter_deriv
!f2py intent(out) :: seq_part_geo
!f2py depend(num_cents) :: seq_part_geo
    integer, parameter :: NUM_PART_CENTS = 3  !number of bra, ket and operator centers
    integer tag_part_cent(NUM_PART_CENTS)     !marking the indices before sorting
    integer num_non_cents                     !number of unique non-atomic centers
    integer id_non_cent(NUM_PART_CENTS)       !position of unique non-atomic centers
    integer num_non_iden(NUM_PART_CENTS)      !number of identical non-atomic centers
    integer num_atom_cents                    !number of unique atomic centers
    integer id_atom_cent(NUM_PART_CENTS)      !position of unique atomic centers
    integer num_atom_iden(NUM_PART_CENTS)     !number of identical atomic centers
#if defined(XTIME)
    real(REALK) curr_time                     !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (num_cents<=0) &
      call error_stop("geom_part_one_param", "invalid num_cents", num_cents)
    ! sorts the centers, and finds the unique ones and corresponding number of identical ones
    tag_part_cent = (/1,2,3/)
    call sort_gen_cents(NUM_PART_CENTS, idx_part_cent, tag_part_cent, &
                        num_non_cents, id_non_cent, num_non_iden,     &
                        num_atom_cents, id_atom_cent, num_atom_iden)
    ! too many differentiated centers
    if (num_cents>num_atom_cents) then
      zero_ints = .true.
      return
    else
      select case(num_non_cents)
      ! there does not exist any non-atomic center
      case(0)
        select case(num_atom_cents)
        ! zeros due to translational invariance, no need to check the differentiated centers
        case(1)
          zero_ints = .true.
          return
        ! there are two unique atomic centers
        case(2)
          select case(num_cents)
          ! one-center total geometric derivatives
          case(1)
            if (idx_cent(1)==idx_part_cent(id_atom_cent(1))) then
              select case(num_atom_iden(1))
              ! partial geometric derivative is performed only on this atomic center
              case(1)
                neg_one = .false.
                seq_part_geo(1) = tag_part_cent(id_atom_cent(1))
              ! translation invariance by performing partial geometric derivatives
              ! on another unique atomic center
              case(2)
                neg_one = mod(order_cent(1),2)==1
                seq_part_geo(1) = tag_part_cent(id_atom_cent(2))
              case default
                call error_stop("geom_part_one_param",              &
                                "invalid num_atom_iden(1) (0/2/1)", &
                                num_atom_iden(1))
              end select
            else if (idx_cent(1)==idx_part_cent(id_atom_cent(2))) then
              select case(num_atom_iden(2))
              ! partial geometric derivative is performed only on this atomic center
              case(1)
                neg_one = .false.
                seq_part_geo(1) = tag_part_cent(id_atom_cent(2))
              ! translation invariance by performing partial geometric derivatives
              ! on another unique atomic center
              case(2)
                neg_one = mod(order_cent(1),2)==1
                seq_part_geo(1) = tag_part_cent(id_atom_cent(1))
              case default
                call error_stop("geom_part_one_param",              &
                                "invalid num_atom_iden(2) (0/2/1)", &
                                num_atom_iden(2))
              end select
            else
              zero_ints = .true.
              return
            end if
            zero_ints = .false.
            ! needs to scatter derivatives when there are partial geometric derivatives
            ! on this center or its succedent centers
            scatter_deriv = any(order_part_cent(seq_part_geo(1):3)/=0)
            order_part_cent(seq_part_geo(1)) = order_part_cent(seq_part_geo(1)) &
                                             + order_cent(1)
          ! two-center total geometric derivatives
          case(2)
            ! notice that both center sequences are sorted in ascending order
            if (idx_cent(1)==idx_part_cent(id_atom_cent(1)) .and. &
                idx_cent(2)==idx_part_cent(id_atom_cent(2))) then
              select case(num_atom_iden(1))
              ! translational invairance by performing partial geometric derivatives
              ! on the first unique atomic center
              case(1)
                neg_one = mod(order_cent(2),2)==1
                seq_part_geo(1) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2) = seq_part_geo(1)
              ! translational invairance by performing partial geometric derivatives
              ! on the second unique atomic center
              case(2)
                neg_one = mod(order_cent(1),2)==1
                seq_part_geo(1) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(2) = seq_part_geo(1)
              case default
                call error_stop("geom_part_one_param",              &
                                "invalid num_atom_iden(1) (0/2/2)", &
                                num_atom_iden(1))
              end select
              zero_ints = .false.
              ! always needs to scatter derivatives
              scatter_deriv = .true.
              order_part_cent(seq_part_geo(1)) = order_part_cent(seq_part_geo(1)) &
                                               + order_cent(1)+order_cent(2)
            else
              zero_ints = .true.
              return
            end if
          ! too many differentiated centers
          case default
            zero_ints = .true.
            return
          end select
        ! there are three unique atomic centers
        case(3)
          select case(num_cents)
          ! one-center total geometric derivatives
          case(1)
            ! partial geometric derivative is performed on the first atomic center
            if (idx_cent(1)==idx_part_cent(id_atom_cent(1))) then
              seq_part_geo(1) = tag_part_cent(id_atom_cent(1))
            ! partial geometric derivative is performed on the second atomic center
            else if (idx_cent(1)==idx_part_cent(id_atom_cent(2))) then
              seq_part_geo(1) = tag_part_cent(id_atom_cent(2))
            ! partial geometric derivative is performed on the third atomic center
            else if (idx_cent(1)==idx_part_cent(id_atom_cent(3))) then
              seq_part_geo(1) = tag_part_cent(id_atom_cent(3))
            else
              zero_ints = .true.
              return
            end if
            zero_ints = .false.
            neg_one = .false.
            ! needs to scatter derivatives when there are partial geometric derivatives
            ! on this center or its succedent centers
            scatter_deriv = any(order_part_cent(seq_part_geo(1):3)/=0)
            order_part_cent(seq_part_geo(1)) = order_part_cent(seq_part_geo(1)) &
                                             + order_cent(1)
          ! two-center total geometric derivatives
          case(2)
            ! notice that both centers are sorted in ascending order, so that
            ! we divide them into two parts according to the third atomic center
            if (idx_cent(2)<idx_part_cent(id_atom_cent(3))) then
              ! partial geometric derivatives on the first and second atomic centers
              if (idx_cent(1)==idx_part_cent(id_atom_cent(1)) .and. &
                  idx_cent(2)==idx_part_cent(id_atom_cent(2))) then
                seq_part_geo(1) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2) = tag_part_cent(id_atom_cent(2))
              else
                zero_ints = .true.
                return
              end if
            else if (idx_cent(2)==idx_part_cent(id_atom_cent(3))) then
              ! partial geometric derivatives on the first and third atomic centers
              if (idx_cent(1)==idx_part_cent(id_atom_cent(1))) then
                seq_part_geo(1) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2) = tag_part_cent(id_atom_cent(3))
              ! partial geometric derivatives on the second and third atomic centers
              else if (idx_cent(1)==idx_part_cent(id_atom_cent(2))) then
                seq_part_geo(1) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(2) = tag_part_cent(id_atom_cent(3))
              else
                zero_ints = .true.
                return
              end if
            else
              zero_ints = .true.
              return
            end if
            zero_ints = .false.
            neg_one = .false.
            ! needs to scatter derivatives when the sequence of differentiated centers
            ! is different from the bra, ket and/or operators centers
            if (seq_part_geo(1)>seq_part_geo(2)) then
              scatter_deriv = .true.
            ! needs to scatter derivatives when there are partial geometric derivatives
            ! on the first differentiated center or its succedent centers
            else
              scatter_deriv = any(order_part_cent(seq_part_geo(1):3)/=0)
            end if
            order_part_cent(seq_part_geo(1)) = order_part_cent(seq_part_geo(1)) &
                                             + order_cent(1)
            order_part_cent(seq_part_geo(2)) = order_part_cent(seq_part_geo(2)) &
                                             + order_cent(2)
          ! three-center total geometric derivatives
          case(3)
            ! notice that both center sequences are sorted in ascending order
            if (idx_cent(1)==idx_part_cent(id_atom_cent(1)) .and. &
                idx_cent(2)==idx_part_cent(id_atom_cent(2)) .and. &
                idx_cent(3)==idx_part_cent(id_atom_cent(3))) then
              zero_ints = .false.
              neg_one = .false.
              seq_part_geo(1) = tag_part_cent(id_atom_cent(1))
              seq_part_geo(2) = tag_part_cent(id_atom_cent(2))
              seq_part_geo(3) = tag_part_cent(id_atom_cent(3))
              ! needs to scatter derivatives when the sequence of differentiated centers
              ! is different from the bra, ket and/or operators centers
              if (seq_part_geo(1)>seq_part_geo(2) .or. &
                  seq_part_geo(2)>seq_part_geo(3)) then
                scatter_deriv = .true.
              ! needs to scatter derivatives when there are partial geometric derivatives
              ! either on bra, ket or operator center
              else
                scatter_deriv = any(order_part_cent/=0)
              end if
              order_part_cent(seq_part_geo(1)) = order_part_cent(seq_part_geo(1)) &
                                               + order_cent(1)
              order_part_cent(seq_part_geo(2)) = order_part_cent(seq_part_geo(2)) &
                                               + order_cent(2)
              order_part_cent(seq_part_geo(3)) = order_part_cent(seq_part_geo(3)) &
                                               + order_cent(3)
            else
              zero_ints = .true.
              return
            end if
          ! too many differentiated centers
          case default
            zero_ints = .true.
            return
          end select
        case default
          call error_stop("geom_part_one_param", &
                          "invalid num_atom_cents (0/X)", num_atom_cents)
        end select
      ! one unique non-atomic center
      case(1)
        select case(num_atom_cents)
        ! there is one unique atomic center (two identical atomic centers)
        case(1)
          ! the only non-zero one-center total geometric derivatives
          if (idx_cent(1)==idx_part_cent(id_atom_cent(1))) then
            ! translation invariance by performing partial geometric derivatives
            ! on the unique non-atomic center
            zero_ints = .false.
            neg_one = mod(order_cent(1),2)==1
            seq_part_geo(1) = tag_part_cent(id_non_cent(1))
            ! needs to scatter derivatives when there are partial geometric derivatives
            ! on this center or its succedent centers
            scatter_deriv = any(order_part_cent(seq_part_geo(1):3)/=0)
            order_part_cent(seq_part_geo(1)) = order_part_cent(seq_part_geo(1)) &
                                             + order_cent(1)
          else
            zero_ints = .true.
            return
          end if
        ! there are two unique atomic centers
        case(2)
          select case(num_cents)
          ! one-center total geometric derivatives
          case(1)
            ! partial geometric derivative is performed on the first atomic center
            if (idx_cent(1)==idx_part_cent(id_atom_cent(1))) then
              seq_part_geo(1) = tag_part_cent(id_atom_cent(1))
            ! partial geometric derivative is performed on the second atomic center
            else if (idx_cent(1)==idx_part_cent(id_atom_cent(2))) then
              seq_part_geo(1) = tag_part_cent(id_atom_cent(2))
            else
              zero_ints = .true.
              return
            end if
            zero_ints = .false.
            neg_one = .false.
            ! needs to scatter derivatives when there are partial geometric derivatives
            ! on this center or its succedent centers
            scatter_deriv = any(order_part_cent(seq_part_geo(1):3)/=0)
            order_part_cent(seq_part_geo(1)) = order_part_cent(seq_part_geo(1)) &
                                             + order_cent(1)
          ! two-center total geometric derivatives
          case(2)
            ! notice that both center sequences are sorted in ascending order
            if (idx_cent(1)==idx_part_cent(id_atom_cent(1)) .and. &
                idx_cent(2)==idx_part_cent(id_atom_cent(2))) then
              zero_ints = .false.
              neg_one = .false.
              seq_part_geo(1) = tag_part_cent(id_atom_cent(1))
              seq_part_geo(2) = tag_part_cent(id_atom_cent(2))
              ! needs to scatter derivatives when the sequence of differentiated centers
              ! is different from the bra, ket and/or operators centers
              if (seq_part_geo(1)>seq_part_geo(2)) then
                scatter_deriv = .true.
              ! needs to scatter derivatives when there are partial geometric derivatives
              ! on the first differentiated center or its succedent centers
              else
                scatter_deriv = any(order_part_cent(seq_part_geo(1):3)/=0)
              end if
              order_part_cent(seq_part_geo(1)) = order_part_cent(seq_part_geo(1)) &
                                               + order_cent(1)
              order_part_cent(seq_part_geo(2)) = order_part_cent(seq_part_geo(2)) &
                                               + order_cent(2)
            else
              zero_ints = .true.
              return
            end if
          ! too many differentiated centers
          case default
            zero_ints = .true.
            return
          end select
        case default
          call error_stop("geom_part_one_param", &
                          "invalid num_atom_cents (1/X)", num_atom_cents)
        end select
      case default
        call error_stop("geom_part_one_param", &
                        "invalid num_non_cents", num_non_cents)
      end select
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "geom_part_one_param", STDOUT)
#endif
    return
  end subroutine geom_part_one_param

  !> \brief scatters geometric derivatives for operators with one center
  !> \author Bin Gao
  !> \date 2011-07-02
  !> \param num_cents is the number of differentiated centers, could be 1, 2 or 3
  !> \param order_cent contains the order of total geometric derivatives of differentiated centers
  !> \param seq_part_geo contains the sequence of bra, ket and operator centers for
  !>        partial derivative terms, from \fn(geom_part_one_param)
  !> \param order_geo_bra is the order of partial geometric derivatives on bra center
  !> \param order_geo_ket is the order of partial geometric derivatives on ket center
  !> \param order_geo_opt is the order of partial geometric derivatives on operator center
  !> \param dim_cints is the dimension of contracted integrals
  !> \param dim_geo_bra is the dimension of geometric derivatives on bra center before scattering
  !> \param dim_geo_ket is the dimension of geometric derivatives on ket center before scattering
  !> \param dim_geo_opt is the dimension of geometric derivatives on operator center before scattering
  !> \param num_opt is the number of other operators
  !> \param part_cints contains the derivatives of integrals before scattering
  !> \param num_geo_bra is the number of partial geometric derivatives on bra center after scattering,
  !>        should equal to \f$(\var(order_geo_bra)+1)(\var(order_geo_bra)+2)/2\f$
  !> \param num_geo_ket is the number of partial geometric derivatives on ket center after scattering,
  !>        should equal to \f$(\var(order_geo_ket)+1)(\var(order_geo_ket)+2)/2\f$
  !> \param num_geo_opt is the number of partial geometric derivatives on operator center after scattering,
  !>        should equal to \f$(\var(order_geo_opt)+1)(\var(order_geo_opt)+2)/2\f$
  !> \param num_tot_geo is the number of total geometric derivatives after scattering
  !> \return total_cints contains the derivatives of integrals after scattering
  subroutine geom_part_one_scatter(num_cents, order_cent, seq_part_geo,         &
                                   order_geo_bra, order_geo_ket, order_geo_opt, &
                                   dim_cints, dim_geo_bra, dim_geo_ket,         &
                                   dim_geo_opt, num_opt, part_cints,            &
                                   num_geo_bra, num_geo_ket, num_geo_opt,       &
                                   num_tot_geo, total_cints)
    use xkind
    implicit none
    integer, intent(in) :: num_cents
    integer, intent(in) :: order_cent(num_cents)
    integer, intent(in) :: seq_part_geo(num_cents)
    integer, intent(in) :: order_geo_bra
    integer, intent(in) :: order_geo_ket
    integer, intent(in) :: order_geo_opt
    integer, intent(in) :: dim_cints
    integer, intent(in) :: dim_geo_bra
    integer, intent(in) :: dim_geo_ket
    integer, intent(in) :: dim_geo_opt
    integer, intent(in) :: num_opt
    real(REALK), intent(in) :: part_cints(dim_cints,dim_geo_bra,dim_geo_ket, &
                                          dim_geo_opt,num_opt)
    integer, intent(in) :: num_geo_bra
    integer, intent(in) :: num_geo_ket
    integer, intent(in) :: num_geo_opt
    integer, intent(in) :: num_tot_geo
    real(REALK), intent(out) :: total_cints(dim_cints,num_geo_bra,num_geo_ket, &
                                            num_geo_opt,num_tot_geo,num_opt)
!f2py intent(hide) :: num_cents
!f2py intent(in) :: order_cent
!f2py intent(in) :: seq_part_geo
!f2py depend(num_cents) :: seq_part_geo
!f2py intent(in) :: order_geo_bra
!f2py intent(in) :: order_geo_ket
!f2py intent(in) :: order_geo_opt
!f2py intent(hide) :: dim_cints
!f2py intent(hide) :: dim_geo_bra
!f2py intent(hide) :: dim_geo_ket
!f2py intent(hide) :: dim_geo_opt
!f2py intent(hide) :: num_opt
!f2py intent(in) :: part_cints
!f2py intent(in) :: num_geo_bra
!f2py intent(in) :: num_geo_ket
!f2py intent(in) :: num_geo_opt
!f2py intent(in) :: num_tot_geo
!f2py intent(out) :: total_cints
!f2py depend(dim_cints) :: total_cints
!f2py depend(num_geo_bra) :: total_cints
!f2py depend(num_geo_ket) :: total_cints
!f2py depend(num_geo_opt) :: total_cints
!f2py depend(num_tot_geo) :: total_cints
!f2py depend(num_opt) :: total_cints
    real(REALK), allocatable :: scatt1_cints(:,:,:,:,:,:)
    real(REALK), allocatable :: scatt2_cints(:,:,:,:,:,:)
                                    !intermediate integrals after scattering derivatives
    integer num_tot_bra             !number of total geometric derivatives on bra center
    integer num_tot_ket             !number of total geometric derivatives on ket center
    integer num_tot_opt             !number of total geometric derivatives on operator center
    integer num_tot_cent2           !number of total geometric derivatives on two centers
    integer nopt                    !incremental recorder over other operators
    integer iopt, ibra, iket, igeo  !incremental recorder over geometric derivatives
    integer ierr                    !error information
#if defined(XTIME)
    real(REALK) curr_time           !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    select case(num_cents)
    ! one-center total geometric derivatives
    case(1)
      select case(seq_part_geo(1))
      ! geometric derivatives on bra center
      case(TAG_BRA)
        ! \var(order_geo_ket)/=0 and/or \var(order_geo_opt)/=0
        if (order_geo_bra==0) then
          do nopt = 1, num_opt
            do ibra = 1, dim_geo_bra
              total_cints(:,1,:,:,ibra,nopt) = part_cints(:,ibra,:,:,nopt)
            end do
          end do
        else
          call scatter_single(dim_cints, dim_geo_bra, dim_geo_ket*dim_geo_opt,   &
                              num_opt, part_cints, order_geo_bra, order_cent(1), &
                              num_geo_bra, num_tot_geo, total_cints)
        end if
      ! geometric derivatives on ket center
      case(TAG_KET)
        ! \var(order_geo_opt)/=0
        if (order_geo_ket==0) then
          do nopt = 1, num_opt
            do iket = 1, dim_geo_ket
              total_cints(:,:,1,:,iket,nopt) = part_cints(:,:,iket,:,nopt)
            end do
          end do
        else
          call scatter_single(dim_cints*dim_geo_bra, dim_geo_ket, dim_geo_opt,   &
                              num_opt, part_cints, order_geo_ket, order_cent(1), &
                              num_geo_ket, num_tot_geo, total_cints)
        end if
      ! geometric derivatives on operator center
      case default
        ! notice that \var(order_geo_opt)/=0
        call scatter_single(dim_cints*dim_geo_bra*dim_geo_ket, dim_geo_opt, 1, &
                            num_opt, part_cints, order_geo_opt, order_cent(1), &
                            num_geo_opt, num_tot_geo, total_cints)
      end select
    ! two-center total geometric derivatives
    case(2)
      select case(seq_part_geo(1))
      ! the first differentiated center is bra center
      case(TAG_BRA)
        select case(seq_part_geo(2))
        ! the second differentiated center is also bra center
        case(TAG_BRA)
          ! no partial geometric derivatives
          if (order_geo_opt==0 .and. order_geo_ket==0 .and. order_geo_bra==0) then
            ! sets the number of total geometric derivatives on bra center
            num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
            num_tot_ket = (order_cent(2)+1)*(order_cent(2)+2)/2
            call scatter_single(dim_cints, dim_geo_bra, 1, num_opt, part_cints, &
                                order_cent(1), order_cent(2), num_tot_bra,      &
                                num_tot_ket, total_cints)
          ! partial geometric derivatives on bra, ket or operator center
          else
            ! sets the number of total geometric derivatives on bra center
            num_tot_opt = (order_cent(2)+1)*(order_cent(2)+2)/2
            ! sets the number of partial geometric derivatives on bra center
            igeo = order_geo_bra+order_cent(1)
            num_tot_ket = (igeo+1)*(igeo+2)/2
            ! allocates the integrals after scattering geometric derivatives on bra center
            allocate(scatt1_cints(dim_cints,num_tot_ket,dim_geo_ket, &
                                  dim_geo_opt,num_tot_opt,num_opt), stat=ierr)
            if (ierr/=0)                                                    &
              call error_stop("geom_part_one_scatter",                      &
                              "failed to allocate scatt1_cints/2/BB/X??",   &
                              dim_cints*num_tot_ket*dim_geo_ket*dim_geo_opt &
                              *num_tot_opt*num_opt)
            ! scatters on bra center
            call scatter_single(dim_cints, dim_geo_bra, dim_geo_ket*dim_geo_opt, &
                                num_opt, part_cints, igeo, order_cent(2),        &
                                num_tot_ket, num_tot_opt, scatt1_cints)
            ! sets the number of total geometric derivatives on bra center
            num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
            ! scatters on bra center
            call scatter_single(dim_cints, num_tot_ket, dim_geo_ket*dim_geo_opt,  &
                                num_tot_opt*num_opt, scatt1_cints, order_geo_bra, &
                                order_cent(1), num_geo_bra, num_tot_bra, total_cints)
            deallocate(scatt1_cints)
          end if
        ! the second differentiated center is ket center
        case(TAG_KET)
          ! no partial geoemtric derivatives on operator center, so that \var(dim_geo_opt)==1
          if (order_geo_opt==0) then
            ! only \var(order_geo_bra)/=0
            if (order_geo_ket==0) then
              ! sets the number of total geometric derivatives on bra center
              num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
              call scatter_single(dim_cints, dim_geo_bra, 1, dim_geo_ket*num_opt, &
                                  part_cints, order_geo_bra, order_cent(1),       &
                                  num_geo_bra, num_tot_bra, total_cints)
            ! at least \var(order_geo_ket)/=0, but \var(order_geo_opt)==0
            else
              ! sets the number of total geometric derivatives on ket center
              num_tot_ket = (order_cent(2)+1)*(order_cent(2)+2)/2
              ! allocates the integrals after scattering geometric derivatives on ket center
              allocate(scatt1_cints(dim_cints,dim_geo_bra,num_geo_ket, &
                                    num_tot_ket,num_opt,1), stat=ierr)
              if (ierr/=0)                                                  &
                call error_stop("geom_part_one_scatter",                    &
                                "failed to allocate scatt1_cints/2/BK/?X0", &
                                dim_cints*dim_geo_bra*num_geo_ket*num_tot_ket*num_opt)
              ! scatters geometric derivatives on ket center
              call scatter_single(dim_cints*dim_geo_bra, dim_geo_ket, 1, num_opt, &
                                  part_cints, order_geo_ket, order_cent(2),       &
                                  num_geo_ket, num_tot_ket, scatt1_cints)
              ! only \var(order_geo_ket)/=0
              if (order_geo_bra==0) then
                ! scatters geometric derivatives on bra center
                do nopt = 1, num_opt
                  igeo = 0
                  do iket = 1, num_tot_ket
                    do ibra = 1, dim_geo_bra
                      igeo = igeo+1
                      total_cints(:,1,:,1,igeo,nopt) &
                        = scatt1_cints(:,ibra,:,iket,nopt,1)
                    end do
                  end do
                end do
              ! only \var(order_geo_opt)==0
              else
                ! sets the number of total geometric derivatives on bra center
                num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, num_geo_ket, &
                                    num_tot_ket*num_opt, scatt1_cints,   &
                                    order_geo_bra, order_cent(1),        &
                                    num_geo_bra, num_tot_bra, total_cints)
              end if
              deallocate(scatt1_cints)
            end if
          ! at least \var(order_geo_opt)/=0, and we should have \var(dim_geo_opt)==\var(num_geo_opt)
          else
            if (order_geo_ket==0) then
              ! only \var(order_geo_opt)/=0
              if (order_geo_bra==0) then
                ! scatters geometric derivatives on bra and ket centers
                do nopt = 1, num_opt
                  do iopt = 1, dim_geo_opt
                    igeo = 0
                    do iket = 1, dim_geo_ket
                      do ibra = 1, dim_geo_bra
                        igeo = igeo+1
                        total_cints(:,1,1,iopt,igeo,nopt) &
                          = part_cints(:,ibra,iket,iopt,nopt)
                      end do
                    end do
                  end do
                end do
              ! only \var(order_geo_ket)==0
              else
                ! allocates the integrals after scattering geometric derivatives on ket center
                allocate(scatt1_cints(dim_cints,dim_geo_bra,dim_geo_opt, &
                                      dim_geo_ket,num_opt,1), stat=ierr)
                if (ierr/=0)                                                  &
                  call error_stop("geom_part_one_scatter",                    &
                                  "failed to allocate scatt1_cints/2/BK/X0X", &
                                  dim_cints*dim_geo_bra*dim_geo_opt*dim_geo_ket*num_opt)
                ! scatters geometric derivatives on ket center
                do nopt = 1, num_opt
                  do iopt = 1, dim_geo_opt
                    do iket = 1, dim_geo_ket
                      scatt1_cints(:,:,iopt,iket,nopt,1) &
                        = part_cints(:,:,iket,iopt,nopt)
                    end do
                  end do
                end do
                ! sets the number of total geometric derivatives on bra center
                num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, dim_geo_opt, &
                                    dim_geo_ket*num_opt, scatt1_cints,   &
                                    order_geo_bra, order_cent(1),        &
                                    num_geo_bra, num_tot_bra, total_cints)
                deallocate(scatt1_cints)
              end if
            ! at least \var(order_geo_ket)/=0 and \var(order_geo_opt)/=0
            else 
              ! sets the number of total geometric derivatives on ket center
              num_tot_ket = (order_cent(2)+1)*(order_cent(2)+2)/2
              ! allocates the integrals after scattering geometric derivatives on ket center
              allocate(scatt1_cints(dim_cints,dim_geo_bra,num_geo_ket,dim_geo_opt, &
                                    num_tot_ket,num_opt), stat=ierr)
              if (ierr/=0)                                                  &
                call error_stop("geom_part_one_scatter",                    &
                                "failed to allocate scatt1_cints/2/BK/?XX", &
                                dim_cints*dim_geo_bra*num_geo_ket           &
                                *dim_geo_opt*num_tot_ket*num_opt)
              ! scatters geometric derivatives on ket center
              call scatter_single(dim_cints*dim_geo_bra, dim_geo_ket, dim_geo_opt,   &
                                  num_opt, part_cints, order_geo_ket, order_cent(2), &
                                  num_geo_ket, num_tot_ket, scatt1_cints)
              ! only \var(order_geo_bra)==0
              if (order_geo_bra==0) then
                ! scatters geometric derivatives on bra center
                do nopt = 1, num_opt
                  igeo = 0
                  do iket = 1, num_tot_ket
                    do ibra = 1, dim_geo_bra
                      igeo = igeo+1
                      total_cints(:,1,:,:,igeo,nopt) &
                        = scatt1_cints(:,ibra,:,:,iket,nopt)
                    end do
                  end do
                end do
              ! none of the orders is zero
              else
                ! sets the number of total geometric derivatives on bra center
                num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, num_geo_ket*dim_geo_opt,  &
                                    num_tot_ket*num_opt, scatt1_cints, order_geo_bra, &
                                    order_cent(1), num_geo_bra, num_tot_bra, total_cints)
              end if
              deallocate(scatt1_cints)
            end if
          end if
        ! the second differentiated center is operator center
        case default
          ! no partial geoemtric derivatives on ket center, so that \var(dim_geo_ket)==1
          if (order_geo_ket==0) then
            ! \var(order_geo_ket)==0 and \var(order_geo_opt)==0, only \var(order_geo_bra)/=0
            if (order_geo_opt==0) then
              ! sets the number of total geometric derivatives on bra center
              num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
              call scatter_single(dim_cints, dim_geo_bra, 1, dim_geo_opt*num_opt, &
                                  part_cints, order_geo_bra, order_cent(1),       &
                                  num_geo_bra, num_tot_bra, total_cints)
            ! \var(order_geo_ket)==0 but \var(order_geo_opt)/=0
            else
              ! sets the number of total geometric derivatives on operator center
              num_tot_opt = (order_cent(2)+1)*(order_cent(2)+2)/2
              ! allocates the integrals after scattering geometric derivatives on operator center
              allocate(scatt1_cints(dim_cints,dim_geo_bra,num_geo_opt, &
                                    num_tot_opt,num_opt,1), stat=ierr)
              if (ierr/=0)                                                  &
                call error_stop("geom_part_one_scatter",                    &
                                "failed to allocate scatt1_cints/2/BC/?X0", &
                                dim_cints*dim_geo_bra*num_geo_opt*num_tot_opt*num_opt)
              ! scatters geometric derivatives on operator center
              call scatter_single(dim_cints*dim_geo_bra, dim_geo_opt, 1, num_opt, &
                                  part_cints, order_geo_opt, order_cent(2),       &
                                  num_geo_opt, num_tot_opt, scatt1_cints)
              ! only \var(order_geo_opt)/=0
              if (order_geo_bra==0) then
                ! scatters geometric derivatives on bra center
                do nopt = 1, num_opt
                  igeo = 0
                  do iopt = 1, num_tot_opt
                    do ibra = 1, dim_geo_bra
                      igeo = igeo+1
                      total_cints(:,1,1,:,igeo,nopt) &
                        = scatt1_cints(:,ibra,:,iopt,nopt,1)
                    end do
                  end do
                end do
              ! only \var(order_geo_ket)==0
              else
                ! sets the number of total geometric derivatives on bra center
                num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, num_geo_opt, &
                                    num_tot_opt*num_opt, scatt1_cints,   &
                                    order_geo_bra, order_cent(1),        &
                                    num_geo_bra, num_tot_bra, total_cints)
              end if
              deallocate(scatt1_cints)
            end if
          ! at least \var(order_geo_ket)/=0
          else
            if (order_geo_opt==0) then
              ! only \var(order_geo_ket)/=0
              if (order_geo_bra==0) then
                ! switches the derivatives on bra and ket centers
                do nopt = 1, num_opt
                  igeo = 0
                  do iopt = 1, dim_geo_opt
                    do ibra = 1, dim_geo_bra
                      igeo = igeo+1
                      total_cints(:,1,:,1,igeo,nopt) = part_cints(:,ibra,:,iopt,nopt)
                    end do
                  end do
                end do
              ! only \var(order_geo_opt)==0
              else
                ! sets the number of total geometric derivatives on bra center
                num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, dim_geo_ket, &
                                    dim_geo_opt*num_opt, part_cints,     &
                                    order_geo_bra, order_cent(1),        &
                                    num_geo_bra, num_tot_bra, total_cints)
              end if
            ! at least \var(order_geo_ket)/=0 and \var(order_geo_opt)/=0
            else 
              ! sets the number of total geometric derivatives on operator center
              num_tot_opt = (order_cent(2)+1)*(order_cent(2)+2)/2
              ! allocates the integrals after scattering geometric derivatives on operator center
              allocate(scatt1_cints(dim_cints,dim_geo_bra,dim_geo_ket,num_geo_opt, &
                                    num_tot_opt,num_opt), stat=ierr)
              if (ierr/=0)                                                  &
                call error_stop("geom_part_one_scatter",                    &
                                "failed to allocate scatt1_cints/2/BC/?XX", &
                                dim_cints*dim_geo_bra*dim_geo_ket           &
                                *num_geo_opt*num_tot_opt*num_opt)
              ! scatters geometric derivatives on operator center
              call scatter_single(dim_cints*dim_geo_bra*dim_geo_ket,   &
                                  dim_geo_opt, 1, num_opt, part_cints, &
                                  order_geo_opt, order_cent(2),        &
                                  num_geo_opt, num_tot_opt, scatt1_cints)
              ! only \var(order_geo_bra)==0
              if (order_geo_bra==0) then
                ! scatters geometric derivatives on bra center
                do nopt = 1, num_opt
                  igeo = 0
                  do iopt = 1, num_tot_opt
                    do ibra = 1, dim_geo_bra
                      igeo = igeo+1
                      total_cints(:,1,:,:,igeo,nopt) &
                        = scatt1_cints(:,ibra,:,:,iopt,nopt)
                    end do
                  end do
                end do
              ! none of the orders is zero
              else
                ! sets the number of total geometric derivatives on bra center
                num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, dim_geo_ket*num_geo_opt,  &
                                    num_tot_opt*num_opt, scatt1_cints, order_geo_bra, &
                                    order_cent(1), num_geo_bra, num_tot_bra, total_cints)
              end if
              deallocate(scatt1_cints)
            end if
          end if
        end select
      ! the first differentiated center is ket center
      case(TAG_KET)
        select case(seq_part_geo(2))
        ! the second differentiated center is bra center
        case(TAG_BRA)
          ! no partial geoemtric derivatives on operator center, so that \var(dim_geo_opt)==1
          if (order_geo_opt==0) then
            ! \var(order_geo_opt)==0 and \var(order_geo_bra)==0
            if (order_geo_bra==0) then
              ! no partial geometric derivatives
              if (order_geo_ket==0) then
                ! switches the derivatives on bra and ket centers
                do nopt = 1, num_opt
                  igeo = 0
                  do ibra = 1, dim_geo_bra
                    do iket = 1, dim_geo_ket
                      igeo = igeo+1
                      total_cints(:,1,1,1,igeo,nopt) = part_cints(:,ibra,iket,1,nopt)
                    end do
                  end do
                end do
              ! only \var(order_geo_ket)/=0
              else
                ! switches the derivatives on bra and ket centers
                allocate(scatt1_cints(dim_cints,dim_geo_ket, &
                                      dim_geo_bra,1,num_opt,1), stat=ierr)
                if (ierr/=0)                                                  &
                  call error_stop("geom_part_one_scatter",                    &
                                  "failed to allocate scatt1_cints/2/KB/X00", &
                                  dim_cints*dim_geo_ket*dim_geo_bra*num_opt)
                do nopt = 1, num_opt
                  do iket = 1, dim_geo_ket
                    do ibra = 1, dim_geo_bra
                      scatt1_cints(:,iket,ibra,1,nopt,1) &
                        = part_cints(:,ibra,iket,1,nopt)
                    end do
                  end do
                end do
                ! sets the number of total geometric derivatives on ket center
                num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
                call scatter_single(dim_cints, dim_geo_ket, 1, dim_geo_bra*num_opt, &
                                    scatt1_cints, order_geo_ket, order_cent(1),     &
                                    num_geo_ket, num_tot_ket, total_cints)
                deallocate(scatt1_cints)
              end if
            ! \var(order_geo_opt)==0 but \var(order_geo_bra)/=0
            else
              ! only \var(order_geo_bra)/=0
              if (order_geo_ket==0) then
                ! sets the number of total geometric derivatives on bra center
                num_tot_bra = (order_cent(2)+1)*(order_cent(2)+2)/2
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, dim_geo_ket, num_opt, &
                                    part_cints, order_geo_bra, order_cent(2),     &
                                    num_geo_bra, num_tot_bra, total_cints)
              ! only \var(order_geo_opt)==0
              else
                ! sets the number of total geometric derivatives on bra center
                num_tot_bra = (order_cent(2)+1)*(order_cent(2)+2)/2
                ! allocates the integrals after scattering geometric derivatives on bra center
                allocate(scatt1_cints(dim_cints,num_geo_bra,dim_geo_ket, &
                                      num_tot_bra,num_opt,1), stat=ierr)
                if (ierr/=0)                                                  &
                  call error_stop("geom_part_one_scatter",                    &
                                  "failed to allocate scatt1_cints/2/KB/XX0", &
                                  dim_cints*num_geo_bra*dim_geo_ket*num_tot_bra*num_opt)
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, dim_geo_ket, num_opt, &
                                    part_cints, order_geo_bra, order_cent(2),     &
                                    num_geo_bra, num_tot_bra, scatt1_cints)
                ! sets the number of total geometric derivatives on ket center
                num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on ket center
                call scatter_single(dim_cints*num_geo_bra, dim_geo_ket, 1, &
                                    num_tot_bra*num_opt, scatt1_cints,     &
                                    order_geo_ket, order_cent(1),          &
                                    num_geo_ket, num_tot_ket, total_cints)
                deallocate(scatt1_cints)
              end if
            end if
          ! at least \var(order_geo_opt)/=0
          else
            if (order_geo_bra==0) then
              ! only \var(order_geo_opt)/=0
              if (order_geo_ket==0) then
                ! scatters geometric derivatives on bra and ket centers
                do nopt = 1, num_opt
                  do iopt = 1, dim_geo_opt
                    igeo = 0
                    do ibra = 1, dim_geo_bra
                      do iket = 1, dim_geo_ket
                        igeo = igeo+1
                        total_cints(:,1,1,iopt,igeo,nopt) &
                          = part_cints(:,ibra,iket,iopt,nopt)
                      end do
                    end do
                  end do
                end do
              ! only \var(order_geo_bra)==0
              else
                ! allocates the integrals after scattering geometric derivatives on bra center
                allocate(scatt1_cints(dim_cints,dim_geo_ket,dim_geo_opt, &
                                      dim_geo_bra,num_opt,1), stat=ierr)
                if (ierr/=0)                                                  &
                  call error_stop("geom_part_one_scatter",                    &
                                  "failed to allocate scatt1_cints/2/KB/X0X", &
                                  dim_cints*dim_geo_ket*dim_geo_opt*dim_geo_bra*num_opt)
                ! scatters geometric derivatives on bra center
                do nopt = 1, num_opt
                  do iopt = 1, dim_geo_opt
                    do iket = 1, dim_geo_ket
                      do ibra = 1, dim_geo_bra
                        scatt1_cints(:,iket,iopt,ibra,nopt,1) &
                          = part_cints(:,ibra,iket,iopt,nopt)
                      end do
                    end do
                  end do
                end do
                ! sets the number of total geometric derivatives on ket center
                num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on ket center
                call scatter_single(dim_cints, dim_geo_ket, dim_geo_opt, &
                                    dim_geo_bra*num_opt, scatt1_cints,   &
                                    order_geo_ket, order_cent(1),        &
                                    num_geo_ket, num_tot_ket, total_cints)
                deallocate(scatt1_cints)
              end if
            ! at least \var(order_geo_opt)/=0 and \var(order_geo_bra)/=0
            else 
              ! sets the number of total geometric derivatives on bra center
              num_tot_bra = (order_cent(2)+1)*(order_cent(2)+2)/2
              ! allocates the integrals after scattering geometric derivatives on bra center
              allocate(scatt1_cints(dim_cints,num_geo_bra,dim_geo_ket,dim_geo_opt, &
                                    num_tot_bra,num_opt), stat=ierr)
              if (ierr/=0)                                                    &
                call error_stop("geom_part_one_scatter",                      &
                                "failed to allocate scatt1_cints/2/KB/?XX",   &
                                dim_cints*num_geo_bra*dim_geo_ket*dim_geo_opt &
                                *num_tot_bra*num_opt)
              ! scatters geometric derivatives on bra center
              call scatter_single(dim_cints, dim_geo_bra, dim_geo_ket*dim_geo_opt,   &
                                  num_opt, part_cints, order_geo_bra, order_cent(2), &
                                  num_geo_bra, num_tot_bra, scatt1_cints)
              ! only \var(order_geo_ket)==0
              if (order_geo_ket==0) then
                ! scatters geometric derivatives on ket center
                do nopt = 1, num_opt
                  igeo = 0
                  do ibra = 1, num_tot_bra
                    do iket = 1, dim_geo_ket
                      igeo = igeo+1
                      total_cints(:,:,1,:,igeo,nopt) &
                        = scatt1_cints(:,:,iket,:,ibra,nopt)
                    end do
                  end do
                end do
              ! none of the orders is zero
              else
                ! sets the number of total geometric derivatives on ket center
                num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on ket center
                call scatter_single(dim_cints*num_geo_bra, dim_geo_ket, dim_geo_opt,  &
                                    num_tot_bra*num_opt, scatt1_cints, order_geo_ket, &
                                    order_cent(1), num_geo_ket, num_tot_ket, total_cints)
              end if
              deallocate(scatt1_cints)
            end if
          end if
        ! the second differentiated center is also ket center
        case(TAG_KET)
          ! no partial geometric derivatives on ket and operator centers
          if (order_geo_opt==0 .and. order_geo_ket==0) then
            ! sets the number of total geometric derivatives on ket center
            num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
            num_tot_ket = (order_cent(2)+1)*(order_cent(2)+2)/2
            call scatter_single(dim_cints*dim_geo_bra, dim_geo_ket, 1, num_opt, &
                                part_cints, order_cent(1), order_cent(2),       &
                                num_tot_bra, num_tot_ket, total_cints)
          ! partial geometric derivatives on ket or operator center
          else
            ! sets the number of total geometric derivatives on ket center
            num_tot_opt = (order_cent(2)+1)*(order_cent(2)+2)/2
            ! sets the number of partial geometric derivatives on ket center
            igeo = order_geo_ket+order_cent(1)
            num_tot_ket = (igeo+1)*(igeo+2)/2
            ! allocates the integrals after scattering geometric derivatives on ket center
            allocate(scatt1_cints(dim_cints,dim_geo_bra,num_tot_ket, &
                                  dim_geo_opt,num_tot_opt,num_opt), stat=ierr)
            if (ierr/=0)                                                    &
              call error_stop("geom_part_one_scatter",                      &
                              "failed to allocate scatt1_cints/2/KK/X??",   &
                              dim_cints*dim_geo_bra*num_tot_ket*dim_geo_opt &
                              *num_tot_opt*num_opt)
            ! scatters on ket center
            call scatter_single(dim_cints*dim_geo_bra, dim_geo_ket, dim_geo_opt, &
                                num_opt, part_cints, igeo, order_cent(2),        &
                                num_tot_ket, num_tot_opt, scatt1_cints)
            ! sets the number of total geometric derivatives on ket center
            num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
            ! scatters on ket center
            call scatter_single(dim_cints*dim_geo_bra, num_tot_ket, dim_geo_opt,  &
                                num_tot_opt*num_opt, scatt1_cints, order_geo_ket, &
                                order_cent(1), num_geo_ket, num_tot_bra, total_cints)
            deallocate(scatt1_cints)
          end if
        ! the second differentiated center is operator center
        case default
          ! \var(order_geo_opt)==0, only \var(order_geo_ket)/=0
          if (order_geo_opt==0) then
            ! sets the number of total geometric derivatives on ket center
            num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
            call scatter_single(dim_cints*dim_geo_bra, dim_geo_ket, 1, &
                                dim_geo_opt*num_opt, part_cints,       &
                                order_geo_ket, order_cent(1),          &
                                num_geo_ket, num_tot_ket, total_cints)
          ! at least \var(order_geo_opt)/=0
          else
            ! sets the number of total geometric derivatives on operator center
            num_tot_opt = (order_cent(2)+1)*(order_cent(2)+2)/2
            ! allocates the integrals after scattering geometric derivatives on operator center
            allocate(scatt1_cints(dim_cints,dim_geo_bra,dim_geo_ket, &
                                  num_geo_opt,num_tot_opt,num_opt), stat=ierr)
            if (ierr/=0)                                                    &
              call error_stop("geom_part_one_scatter",                      &
                              "failed to allocate scatt1_cints/2/KC/?X?",   &
                              dim_cints*dim_geo_bra*dim_geo_ket*num_geo_opt &
                              *num_tot_opt*num_opt)
            ! scatters geometric derivatives on operator center
            call scatter_single(dim_cints*dim_geo_bra*dim_geo_ket, dim_geo_opt, 1, &
                                num_opt, part_cints, order_geo_opt, order_cent(2), &
                                num_geo_opt, num_tot_opt, scatt1_cints)
            ! only \var(order_geo_opt)/=0
            if (order_geo_ket==0) then
              ! scatters geometric derivatives on ket center
              do nopt = 1, num_opt
                igeo = 0
                do iopt = 1, num_tot_opt
                  do iket = 1, dim_geo_ket
                    igeo = igeo+1
                    total_cints(:,:,1,:,igeo,nopt) = scatt1_cints(:,:,iket,:,iopt,nopt)
                  end do
                end do
              end do
            ! both \var(order_geo_ket)/=0 and \var(order_geo_opt)/=0
            else
              ! sets the number of total geometric derivatives on ket center
              num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
              ! scatters geometric derivatives on ket center
              call scatter_single(dim_cints*dim_geo_bra, dim_geo_ket, num_geo_opt,  &
                                  num_tot_opt*num_opt, scatt1_cints, order_geo_ket, &
                                  order_cent(1), num_geo_ket, num_tot_ket, total_cints)
            end if
            deallocate(scatt1_cints)
          end if
        end select
      ! the first differentiated center is operator center
      case default
        select case(seq_part_geo(2))
        ! the second differentiated center is bra center
        case(TAG_BRA)
          ! \var(order_geo_bra)==0
          if (order_geo_bra==0) then
            ! no partial geometric derivatives
            if (order_geo_opt==0) then
              ! switches the derivatives on bra and operator centers
              do nopt = 1, num_opt
                igeo = 0
                do ibra = 1, dim_geo_bra
                  do iopt = 1, dim_geo_opt
                    igeo = igeo+1
                    total_cints(:,1,:,1,igeo,nopt) = part_cints(:,ibra,:,iopt,nopt)
                  end do
                end do
              end do
            ! only \var(order_geo_opt)/=0
            else
              ! switches the derivatives on bra and ket, operator centers
              allocate(scatt1_cints(dim_cints,dim_geo_ket,dim_geo_opt, &
                                    dim_geo_bra,num_opt,1), stat=ierr)
              if (ierr/=0)                                                  &
                call error_stop("geom_part_one_scatter",                    &
                                "failed to allocate scatt1_cints/2/CB/X0?", &
                                dim_cints*dim_geo_ket*dim_geo_opt*dim_geo_bra*num_opt)
              do nopt = 1, num_opt
                do iopt = 1, dim_geo_opt
                  do iket = 1, dim_geo_ket
                    do ibra = 1, dim_geo_bra
                      scatt1_cints(:,iket,iopt,ibra,nopt,1) &
                        = part_cints(:,ibra,iket,iopt,nopt)
                    end do
                  end do
                end do
              end do
              ! sets the number of total geometric derivatives on operator center
              num_tot_opt = (order_cent(1)+1)*(order_cent(1)+2)/2
              call scatter_single(dim_cints*dim_geo_ket, dim_geo_opt, 1, &
                                  dim_geo_bra*num_opt, scatt1_cints,     &
                                  order_geo_opt, order_cent(1),          &
                                  num_geo_opt, num_tot_opt, total_cints)
              deallocate(scatt1_cints)
            end if
          ! \var(order_geo_bra)/=0
          else
            ! sets the number of total geometric derivatives on bra center
            num_tot_bra = (order_cent(2)+1)*(order_cent(2)+2)/2
            ! only \var(order_geo_bra)/=0
            if (order_geo_opt==0) then
              ! scatters geometric derivatives on bra center
              call scatter_single(dim_cints, dim_geo_bra, dim_geo_ket*dim_geo_opt,   &
                                  num_opt, part_cints, order_geo_bra, order_cent(2), &
                                  num_geo_bra, num_tot_bra, total_cints)
            ! both \var(order_geo_bra)/=0 and \var(order_geo_opt)/=0
            else
              ! allocates the integrals after scattering geometric derivatives on bra center
              allocate(scatt1_cints(dim_cints,num_geo_bra,dim_geo_ket, &
                                    dim_geo_opt,num_tot_bra,num_opt), stat=ierr)
              if (ierr/=0)                                                    &
                call error_stop("geom_part_one_scatter",                      &
                                "failed to allocate scatt1_cints/2/CB/XX?",   &
                                dim_cints*num_geo_bra*dim_geo_ket*dim_geo_opt &
                                *num_tot_bra*num_opt)
              ! scatters geometric derivatives on bra center
              call scatter_single(dim_cints, dim_geo_bra, dim_geo_ket*dim_geo_opt,   &
                                  num_opt, part_cints, order_geo_bra, order_cent(2), &
                                  num_geo_bra, num_tot_bra, scatt1_cints)
              ! sets the number of total geometric derivatives on operator center
              num_tot_opt = (order_cent(1)+1)*(order_cent(1)+2)/2
              ! scatters geometric derivatives on operator center
              call scatter_single(dim_cints*num_geo_bra*dim_geo_ket, dim_geo_opt, 1, &
                                  num_tot_bra*num_opt, scatt1_cints, order_geo_opt,  &
                                  order_cent(1), num_geo_opt, num_tot_opt, total_cints)
              deallocate(scatt1_cints)
            end if
          end if
        ! the second differentiated center is ket center
        case(TAG_KET)
          ! \var(order_geo_ket)==0
          if (order_geo_ket==0) then
            ! no partial geometric derivatives
            if (order_geo_opt==0) then
              ! switches the derivatives on ket and operator centers
              do nopt = 1, num_opt
                igeo = 0
                do iket = 1, dim_geo_ket
                  do iopt = 1, dim_geo_opt
                    igeo = igeo+1
                    total_cints(:,:,1,1,igeo,nopt) = part_cints(:,:,iket,iopt,nopt)
                  end do
                end do
              end do
            ! only \var(order_geo_opt)/=0
            else
              ! switches the derivatives on ket and operator centers
              allocate(scatt1_cints(dim_cints,dim_geo_bra,dim_geo_opt, &
                                    dim_geo_ket,num_opt,1), stat=ierr)
              if (ierr/=0)                                                  &
                call error_stop("geom_part_one_scatter",                    &
                                "failed to allocate scatt1_cints/2/CK/X0?", &
                                dim_cints*dim_geo_bra*dim_geo_opt*dim_geo_ket*num_opt)
              do nopt = 1, num_opt
                do iopt = 1, dim_geo_opt
                  do iket = 1, dim_geo_ket
                    scatt1_cints(:,:,iopt,iket,nopt,1) = part_cints(:,:,iket,iopt,nopt)
                  end do
                end do
              end do
              ! sets the number of total geometric derivatives on operator center
              num_tot_opt = (order_cent(1)+1)*(order_cent(1)+2)/2
              call scatter_single(dim_cints*dim_geo_bra, dim_geo_opt, 1, &
                                  dim_geo_ket*num_opt, scatt1_cints,     &
                                  order_geo_opt, order_cent(1),          &
                                  num_geo_opt, num_tot_opt, total_cints)
              deallocate(scatt1_cints)
            end if
          ! at least \var(order_geo_ket)/=0
          else
            ! sets the number of total geometric derivatives on ket center
            num_tot_ket = (order_cent(2)+1)*(order_cent(2)+2)/2
            ! only \var(order_geo_ket)/=0
            if (order_geo_opt==0) then
              ! scatters geometric derivatives on ket center
              call scatter_single(dim_cints*dim_geo_bra, dim_geo_ket, dim_geo_opt,   &
                                  num_opt, part_cints, order_geo_ket, order_cent(2), &
                                  num_geo_ket, num_tot_ket, total_cints)
            ! both \var(order_geo_ket)/=0 and \var(order_geo_opt)/=0
            else
              ! allocates the integrals after scattering geometric derivatives on ket center
              allocate(scatt1_cints(dim_cints,dim_geo_bra,num_geo_ket, &
                                    dim_geo_opt,num_tot_ket,num_opt), stat=ierr)
              if (ierr/=0)                                                    &
                call error_stop("geom_part_one_scatter",                      &
                                "failed to allocate scatt1_cints/2/CK/XX?",   &
                                dim_cints*dim_geo_bra*num_geo_ket*dim_geo_opt &
                                *num_tot_ket*num_opt)
              ! scatters geometric derivatives on ket center
              call scatter_single(dim_cints*dim_geo_bra, dim_geo_ket, dim_geo_opt,   &
                                  num_opt, part_cints, order_geo_ket, order_cent(2), &
                                  num_geo_ket, num_tot_ket, scatt1_cints)
              ! sets the number of total geometric derivatives on operator center
              num_tot_opt = (order_cent(1)+1)*(order_cent(1)+2)/2
              ! scatters geometric derivatives on operator center
              call scatter_single(dim_cints*dim_geo_bra*num_geo_ket, dim_geo_opt, 1, &
                                  num_tot_ket*num_opt, scatt1_cints, order_geo_opt,  &
                                  order_cent(1), num_geo_opt, num_tot_opt, total_cints)
              deallocate(scatt1_cints)
            end if
          end if
        ! the second differentiated center is also operator center
        case default
          ! no partial geometric derivatives on operator center
          if (order_geo_opt==0) then
            ! sets the number of total geometric derivatives on operator center
            num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
            num_tot_ket = (order_cent(2)+1)*(order_cent(2)+2)/2
            call scatter_single(dim_cints*dim_geo_bra*dim_geo_ket, dim_geo_opt, 1, &
                                num_opt, part_cints, order_cent(1), order_cent(2), &
                                num_tot_bra, num_tot_ket, total_cints)
          ! partial geometric derivatives on operator center
          else
            ! sets the number of total geometric derivatives on operator center
            num_tot_opt = (order_cent(2)+1)*(order_cent(2)+2)/2
            ! sets the number of partial geometric derivatives on operator center
            igeo = order_geo_opt+order_cent(1)
            num_tot_ket = (igeo+1)*(igeo+2)/2
            ! allocates the integrals after scattering geometric derivatives on operator center
            allocate(scatt1_cints(dim_cints,dim_geo_bra,dim_geo_ket, &
                                  num_tot_ket,num_tot_opt,num_opt), stat=ierr)
            if (ierr/=0)                                                    &
              call error_stop("geom_part_one_scatter",                      &
                              "failed to allocate scatt1_cints/2/CC/X??",   &
                              dim_cints*dim_geo_bra*dim_geo_ket*num_tot_ket &
                              *num_tot_opt*num_opt)
            ! scatters on operator center
            call scatter_single(dim_cints*dim_geo_bra*dim_geo_ket, dim_geo_opt, &
                                1, num_opt, part_cints, igeo, order_cent(2),    &
                                num_tot_ket, num_tot_opt, scatt1_cints)
            ! sets the number of total geometric derivatives on operator center
            num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
            ! scatters on operator center
            call scatter_single(dim_cints*dim_geo_bra*dim_geo_ket, num_tot_ket, 1, &
                                num_tot_opt*num_opt, scatt1_cints, order_geo_opt,  &
                                order_cent(1), num_geo_opt, num_tot_bra, total_cints)
            deallocate(scatt1_cints)
          end if
        end select
      end select
    ! three-center total geometric derivatives
    case(3)
      select case(seq_part_geo(1))
      ! the first differentiated center is bra center
      case(TAG_BRA)
        select case(seq_part_geo(2))
        ! the second differentiated center is ket center, and the third is operator center
        case(TAG_KET)
          if (order_geo_opt==0) then
            ! only \var(order_geo_bra)/=0
            if (order_geo_ket==0) then
              ! sets the number of total geometric derivatives on bra center
              num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
              call scatter_single(dim_cints, dim_geo_bra, 1,                &
                                  dim_geo_ket*dim_geo_opt*num_opt,          &
                                  part_cints, order_geo_bra, order_cent(1), &
                                  num_geo_bra, num_tot_bra, total_cints)
            ! at least \var(order_geo_ket)/=0, but \var(order_geo_opt)==0
            else
              ! sets the number of total geometric derivatives on ket center
              num_tot_ket = (order_cent(2)+1)*(order_cent(2)+2)/2
              ! allocates the integrals after scattering geometric derivatives on ket center
              allocate(scatt1_cints(dim_cints,dim_geo_bra,num_geo_ket, &
                                    num_tot_ket,dim_geo_opt,num_opt), stat=ierr)
              if (ierr/=0)                                                    &
                call error_stop("geom_part_one_scatter",                      &
                                "failed to allocate scatt1_cints/3/BKC/?X0",  &
                                dim_cints*dim_geo_bra*num_geo_ket*num_tot_ket &
                                *dim_geo_opt*num_opt)
              ! scatters geometric derivatives on ket center
              call scatter_single(dim_cints*dim_geo_bra, dim_geo_ket, &
                                  1, dim_geo_opt*num_opt, part_cints, &
                                  order_geo_ket, order_cent(2),       &
                                  num_geo_ket, num_tot_ket, scatt1_cints)
              ! only \var(order_geo_ket)/=0
              if (order_geo_bra==0) then
                ! scatters geometric derivatives on bra center
                do nopt = 1, num_opt
                  igeo = 0
                  do iopt = 1, dim_geo_opt
                    do iket = 1, num_tot_ket
                      do ibra = 1, dim_geo_bra
                        igeo = igeo+1
                        total_cints(:,1,:,1,igeo,nopt) &
                          = scatt1_cints(:,ibra,:,iket,iopt,nopt)
                      end do
                    end do
                  end do
                end do
              ! only \var(order_geo_opt)==0
              else
                ! sets the number of total geometric derivatives on bra center
                num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, num_geo_ket,        &
                                    num_tot_ket*dim_geo_opt*num_opt,            &
                                    scatt1_cints, order_geo_bra, order_cent(1), &
                                    num_geo_bra, num_tot_bra, total_cints)
              end if
              deallocate(scatt1_cints)
            end if
          ! at least \var(order_geo_opt)/=0
          else
            ! sets the number of total geometric derivatives on operator center
            num_tot_opt = (order_cent(3)+1)*(order_cent(3)+2)/2
            ! allocates the integrals after scattering geometric derivatives on operator center
            allocate(scatt1_cints(dim_cints,dim_geo_bra,dim_geo_ket, &
                                  num_geo_opt,num_tot_opt,num_opt), stat=ierr)
            if (ierr/=0)                                                    &
              call error_stop("geom_part_one_scatter",                      &
                              "failed to allocate scatt1_cints/3/BKC/??X",  &
                              dim_cints*dim_geo_bra*dim_geo_ket*num_geo_opt &
                              *num_tot_opt*num_opt)
            ! scatters geometric derivatives on operator center
            call scatter_single(dim_cints*dim_geo_bra*dim_geo_ket, dim_geo_opt, 1, &
                                num_opt, part_cints, order_geo_opt, order_cent(3), &
                                num_geo_opt, num_tot_opt, scatt1_cints)
            if (order_geo_ket==0) then
              ! only \var(order_geo_opt)/=0
              if (order_geo_bra==0) then
                ! scatters geometric derivatives on bra and ket centers
                do nopt = 1, num_opt
                  igeo = 0
                  do iopt = 1, num_tot_opt
                    do iket = 1, dim_geo_ket
                      do ibra = 1, dim_geo_bra
                        igeo = igeo+1
                        total_cints(:,1,1,:,igeo,nopt) &
                          = scatt1_cints(:,ibra,iket,:,iopt,nopt)
                      end do
                    end do
                  end do
                end do
                deallocate(scatt1_cints)
              ! only \var(order_geo_ket)==0
              else
                ! allocates the integrals after scattering geometric derivatives on ket center
                allocate(scatt2_cints(dim_cints,dim_geo_bra,num_geo_opt, &
                                      dim_geo_ket,num_tot_opt,num_opt), stat=ierr)
                if (ierr/=0)                                                    &
                  call error_stop("geom_part_one_scatter",                      &
                                  "failed to allocate scatt2_cints/3/BKC/X0X",  &
                                  dim_cints*dim_geo_bra*num_geo_opt*dim_geo_ket &
                                  *num_tot_opt*num_opt)
                ! scatters geometric derivatives on ket center
                do nopt = 1, num_opt
                  do iopt = 1, num_tot_opt
                    do iket = 1, dim_geo_ket
                      scatt2_cints(:,:,:,iket,iopt,nopt) &
                        = scatt1_cints(:,:,iket,:,iopt,nopt)
                    end do
                  end do
                end do
                deallocate(scatt1_cints)
                ! sets the number of total geometric derivatives on bra center
                num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, num_geo_opt,           &
                                    dim_geo_ket*num_tot_opt*num_opt, scatt2_cints, &
                                    order_geo_bra, order_cent(1),                  &
                                    num_geo_bra, num_tot_bra, total_cints)
                deallocate(scatt2_cints)
              end if
            ! at least \var(order_geo_ket)/=0 and \var(order_geo_opt)/=0
            else 
              ! sets the number of total geometric derivatives on ket center
              num_tot_ket = (order_cent(2)+1)*(order_cent(2)+2)/2
              num_tot_cent2 = num_tot_ket*num_tot_opt
              ! allocates the integrals after scattering geometric derivatives on ket center
              allocate(scatt2_cints(dim_cints,dim_geo_bra,num_geo_ket,num_geo_opt, &
                                    num_tot_cent2,num_opt), stat=ierr)
              if (ierr/=0)                                                    &
                call error_stop("geom_part_one_scatter",                      &
                                "failed to allocate scatt2_cints/3/BKC/?XX",  &
                                dim_cints*dim_geo_bra*num_geo_ket*num_geo_opt &
                                *num_tot_cent2*num_opt)
              ! scatters geometric derivatives on ket center
              call scatter_single(dim_cints*dim_geo_bra, dim_geo_ket, num_geo_opt,  &
                                  num_tot_opt*num_opt, scatt1_cints, order_geo_ket, &
                                  order_cent(2), num_geo_ket, num_tot_ket, scatt2_cints)
              deallocate(scatt1_cints)
              ! only \var(order_geo_bra)==0
              if (order_geo_bra==0) then
                ! scatters geometric derivatives on bra center
                do nopt = 1, num_opt
                  igeo = 0
                  do iket = 1, num_tot_cent2
                    do ibra = 1, dim_geo_bra
                      igeo = igeo+1
                      total_cints(:,1,:,:,igeo,nopt) &
                        = scatt2_cints(:,ibra,:,:,iket,nopt)
                    end do
                  end do
                end do
              ! none of the orders is zero
              else
                ! sets the number of total geometric derivatives on bra center
                num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, num_geo_ket*num_geo_opt,    &
                                    num_tot_cent2*num_opt, scatt2_cints, order_geo_bra, &
                                    order_cent(1), num_geo_bra, num_tot_bra, total_cints)
              end if
              deallocate(scatt2_cints)
            end if
          end if
        ! the second differentiated center is operator center, and the third is ket center
        case default
          if (order_geo_ket==0) then
            ! \var(order_geo_ket)==0 and \var(order_geo_opt)==0
            if (order_geo_opt==0) then
              ! no partial geometric derivatives
              if (order_geo_bra==0) then
                ! switches the derivatives on ket and operator centers
                do nopt = 1, num_opt
                  igeo = 0
                  do iket = 1, dim_geo_ket
                    do iopt = 1, dim_geo_opt
                      do ibra = 1, dim_geo_bra
                        igeo = igeo+1
                        total_cints(:,1,1,1,igeo,nopt) &
                          = part_cints(:,ibra,iket,iopt,nopt)
                      end do
                    end do
                  end do
                end do
              ! only \var(order_geo_bra)/=0
              else
                ! switches the derivatives on ket and operator centers
                allocate(scatt1_cints(dim_cints,dim_geo_bra,dim_geo_opt, &
                                      dim_geo_ket,num_opt,1), stat=ierr)
                if (ierr/=0)                                                   &
                  call error_stop("geom_part_one_scatter",                     &
                                  "failed to allocate scatt1_cints/3/BCK/X00", &
                                  dim_cints*dim_geo_bra*dim_geo_opt*dim_geo_ket*num_opt)
                do nopt = 1, num_opt
                  do iket = 1, dim_geo_ket
                    do iopt = 1, dim_geo_opt
                      scatt1_cints(:,:,iopt,iket,nopt,1) &
                        = part_cints(:,:,iket,iopt,nopt)
                    end do
                  end do
                end do
                ! sets the number of total geometric derivatives on bra center
                num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
                call scatter_single(dim_cints, dim_geo_bra, 1,                  &
                                    dim_geo_opt*dim_geo_ket*num_opt,            &
                                    scatt1_cints, order_geo_bra, order_cent(1), &
                                    num_geo_bra, num_tot_bra, total_cints)
                deallocate(scatt1_cints)
              end if
            ! \var(order_geo_ket)==0 but \var(order_geo_opt)/=0
            else
              ! switches the derivatives on ket and operator centers
              allocate(scatt1_cints(dim_cints,dim_geo_bra,dim_geo_opt, &
                                    dim_geo_ket,num_opt,1), stat=ierr)
              if (ierr/=0)                                                   &
                call error_stop("geom_part_one_scatter",                     &
                                "failed to allocate scatt1_cints/3/BCK/?X0", &
                                dim_cints*dim_geo_bra*dim_geo_opt*dim_geo_ket*num_opt)
              do nopt = 1, num_opt
                do iopt = 1, dim_geo_opt
                  do iket = 1, dim_geo_ket
                    scatt1_cints(:,:,iopt,iket,nopt,1) = part_cints(:,:,iket,iopt,nopt)
                  end do
                end do
              end do
              ! sets the number of total geometric derivatives on operator center
              num_tot_opt = (order_cent(2)+1)*(order_cent(2)+2)/2
              ! allocates the integrals after scattering geometric derivatives on operator center
              allocate(scatt2_cints(dim_cints,dim_geo_bra,num_geo_opt, &
                                    num_tot_opt,dim_geo_ket,num_opt), stat=ierr)
              if (ierr/=0)                                                    &
                call error_stop("geom_part_one_scatter",                      &
                                "failed to allocate scatt2_cints/3/BCK/?X0",  &
                                dim_cints*dim_geo_bra*num_geo_opt*num_tot_opt &
                                *dim_geo_ket*num_opt)
              ! scatters geometric derivatives on operator center
              call scatter_single(dim_cints*dim_geo_bra, dim_geo_opt,   &
                                  1, dim_geo_ket*num_opt, scatt1_cints, &
                                  order_geo_opt, order_cent(2),         &
                                  num_geo_opt, num_tot_opt, scatt2_cints)
              deallocate(scatt1_cints)
              ! only \var(order_geo_opt)/=0
              if (order_geo_bra==0) then
                ! scatters geometric derivatives on bra center
                do nopt = 1, num_opt
                  igeo = 0
                  do iket = 1, dim_geo_ket
                    do iopt = 1, num_tot_opt
                      do ibra = 1, dim_geo_bra
                        igeo = igeo+1
                        total_cints(:,1,1,:,igeo,nopt) &
                          = scatt2_cints(:,ibra,:,iopt,iket,nopt)
                      end do
                    end do
                  end do
                end do
              ! only \var(order_geo_ket)==0
              else
                ! sets the number of total geometric derivatives on bra center
                num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, num_geo_opt,        &
                                    num_tot_opt*dim_geo_ket*num_opt,            &
                                    scatt2_cints, order_geo_bra, order_cent(1), &
                                    num_geo_bra, num_tot_bra, total_cints)
              end if
              deallocate(scatt2_cints)
            end if
          ! at least \var(order_geo_ket)/=0
          else
            ! sets the number of total geometric derivatives on ket center
            num_tot_ket = (order_cent(3)+1)*(order_cent(3)+2)/2
            ! allocates the integrals after scattering geometric derivatives on ket center
            allocate(scatt1_cints(dim_cints,dim_geo_bra,num_geo_ket, &
                                  dim_geo_opt,num_tot_ket,num_opt), stat=ierr)
            if (ierr/=0)                                                    &
              call error_stop("geom_part_one_scatter",                      &
                              "failed to allocate scatt1_cints/3/BCK/??X",  &
                              dim_cints*dim_geo_bra*num_geo_ket*dim_geo_opt &
                              *num_tot_ket*num_opt)
            ! scatters geometric derivatives on ket center
            call scatter_single(dim_cints*dim_geo_bra, dim_geo_ket, dim_geo_opt,   &
                                num_opt, part_cints, order_geo_ket, order_cent(3), &
                                num_geo_ket, num_tot_ket, scatt1_cints)
            if (order_geo_opt==0) then
              ! only \var(order_geo_ket)/=0
              if (order_geo_bra==0) then
                ! scatters geometric derivatives on bra center
                do nopt = 1, num_opt
                  igeo = 0
                  do iket = 1, num_tot_ket
                    do iopt = 1, dim_geo_opt
                      do ibra = 1, dim_geo_bra
                        igeo = igeo+1
                        total_cints(:,1,:,1,igeo,nopt) &
                          = scatt1_cints(:,ibra,:,iopt,iket,nopt)
                      end do
                    end do
                  end do
                end do
              ! only \var(order_geo_opt)==0
              else
                ! sets the number of total geometric derivatives on bra center
                num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, num_geo_ket,        &
                                    dim_geo_opt*num_tot_ket*num_opt,            &
                                    scatt1_cints, order_geo_bra, order_cent(1), &
                                    num_geo_bra, num_tot_bra, total_cints)
              end if
              deallocate(scatt1_cints)
            ! at least \var(order_geo_ket)/=0 and \var(order_geo_opt)/=0
            else 
              ! sets the number of total geometric derivatives on operator center
              num_tot_opt = (order_cent(2)+1)*(order_cent(2)+2)/2
              num_tot_cent2 = num_tot_opt*num_tot_ket
              ! allocates the integrals after scattering geometric derivatives on operator center
              allocate(scatt2_cints(dim_cints,dim_geo_bra,num_geo_ket,num_geo_opt, &
                                    num_tot_cent2,num_opt), stat=ierr)
              if (ierr/=0)                                                    &
                call error_stop("geom_part_one_scatter",                      &
                                "failed to allocate scatt2_cints/3/BCK/?XX",  &
                                dim_cints*dim_geo_bra*num_geo_ket*num_geo_opt &
                                *num_tot_cent2*num_opt)
              ! scatters geometric derivatives on operator center
              call scatter_single(dim_cints*dim_geo_bra*num_geo_ket,          &
                                  dim_geo_opt, 1, num_tot_ket*num_opt,        &
                                  scatt1_cints, order_geo_opt, order_cent(2), &
                                  num_geo_opt, num_tot_opt, scatt2_cints)
              deallocate(scatt1_cints)
              ! only \var(order_geo_bra)==0
              if (order_geo_bra==0) then
                ! scatters geometric derivatives on bra center
                do nopt = 1, num_opt
                  igeo = 0
                  do iopt = 1, num_tot_cent2
                    do ibra = 1, dim_geo_bra
                      igeo = igeo+1
                      total_cints(:,1,:,:,igeo,nopt) &
                        = scatt2_cints(:,ibra,:,:,iopt,nopt)
                    end do
                  end do
                end do
              ! none of the orders is zero
              else
                ! sets the number of total geometric derivatives on bra center
                num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, num_geo_ket*num_geo_opt, &
                                    num_tot_cent2*num_opt, scatt2_cints,             &
                                    order_geo_bra, order_cent(1), num_geo_bra,       &
                                    num_tot_bra, total_cints)
              end if
              deallocate(scatt2_cints)
            end if
          end if
        end select
      ! the first differentiated center is ket center
      case(TAG_KET)
        select case(seq_part_geo(2))
        ! the second differentiated center is bra center, and the third is operator center
        case(TAG_BRA)
          if (order_geo_opt==0) then
            ! \var(order_geo_opt)==0 and \var(order_geo_bra)==0
            if (order_geo_bra==0) then
              ! no partial geometric derivatives
              if (order_geo_ket==0) then
                ! switches the derivatives on bra and ket centers
                do nopt = 1, num_opt
                  igeo = 0
                  do iopt = 1, dim_geo_opt
                    do ibra = 1, dim_geo_bra
                      do iket = 1, dim_geo_ket
                        igeo = igeo+1
                        total_cints(:,1,1,1,igeo,nopt) &
                          = part_cints(:,ibra,iket,iopt,nopt)
                      end do
                    end do
                  end do
                end do
              ! only \var(order_geo_ket)/=0
              else
                ! switches the derivatives on bra and ket centers
                allocate(scatt1_cints(dim_cints,dim_geo_ket,dim_geo_bra, &
                                      dim_geo_opt,num_opt,1), stat=ierr)
                if (ierr/=0)                                                   &
                  call error_stop("geom_part_one_scatter",                     &
                                  "failed to allocate scatt1_cints/3/KBC/X00", &
                                  dim_cints*dim_geo_ket*dim_geo_bra*dim_geo_opt*num_opt)
                do nopt = 1, num_opt
                  do iopt = 1, dim_geo_opt
                    do iket = 1, dim_geo_ket
                      do ibra = 1, dim_geo_bra
                        scatt1_cints(:,iket,ibra,iopt,nopt,1) &
                          = part_cints(:,ibra,iket,iopt,nopt)
                      end do
                    end do
                  end do
                end do
                ! sets the number of total geometric derivatives on ket center
                num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
                call scatter_single(dim_cints, dim_geo_ket, 1,                  &
                                    dim_geo_bra*dim_geo_opt*num_opt,            &
                                    scatt1_cints, order_geo_ket, order_cent(1), &
                                    num_geo_ket, num_tot_ket, total_cints)
                deallocate(scatt1_cints)
              end if
            ! \var(order_geo_opt)==0 but \var(order_geo_bra)/=0
            else
              ! only \var(order_geo_bra)/=0
              if (order_geo_ket==0) then
                ! sets the number of total geometric derivatives on bra center
                num_tot_bra = (order_cent(2)+1)*(order_cent(2)+2)/2
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, dim_geo_ket, &
                                    dim_geo_opt*num_opt, part_cints,     &
                                    order_geo_bra, order_cent(2),        &
                                    num_geo_bra, num_tot_bra, total_cints)
              ! only \var(order_geo_opt)==0
              else
                ! sets the number of total geometric derivatives on bra center
                num_tot_bra = (order_cent(2)+1)*(order_cent(2)+2)/2
                ! allocates the integrals after scattering geometric derivatives on bra center
                allocate(scatt1_cints(dim_cints,num_geo_bra,dim_geo_ket, &
                                      num_tot_bra,dim_geo_opt,num_opt), stat=ierr)
                if (ierr/=0)                                                    &
                  call error_stop("geom_part_one_scatter",                      &
                                  "failed to allocate scatt1_cints/3/KBC/XX0",  &
                                  dim_cints*num_geo_bra*dim_geo_ket*num_tot_bra &
                                  *dim_geo_opt*num_opt)
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, dim_geo_ket, &
                                    dim_geo_opt*num_opt, part_cints,     &
                                    order_geo_bra, order_cent(2),        &
                                    num_geo_bra, num_tot_bra, scatt1_cints)
                ! sets the number of total geometric derivatives on ket center
                num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on ket center
                call scatter_single(dim_cints*num_geo_bra, dim_geo_ket, 1,      &
                                    num_tot_bra*dim_geo_opt*num_opt,            &
                                    scatt1_cints, order_geo_ket, order_cent(1), &
                                    num_geo_ket, num_tot_ket, total_cints)
                deallocate(scatt1_cints)
              end if
            end if
          ! at least \var(order_geo_opt)/=0
          else
            ! sets the number of total geometric derivatives on operator center
            num_tot_opt = (order_cent(3)+1)*(order_cent(3)+2)/2
            ! allocates the integrals after scattering geometric derivatives on operator center
            allocate(scatt1_cints(dim_cints,dim_geo_bra,dim_geo_ket, &
                                  num_geo_opt,num_tot_opt,num_opt), stat=ierr)
            if (ierr/=0)                                                    &
              call error_stop("geom_part_one_scatter",                      &
                              "failed to allocate scatt1_cints/3/KBC/??X",  &
                              dim_cints*dim_geo_bra*dim_geo_ket*num_geo_opt &
                              *num_tot_opt*num_opt)
            ! scatters geometric derivatives on operator center
            call scatter_single(dim_cints*dim_geo_bra*dim_geo_ket, dim_geo_opt, 1, &
                                num_opt, part_cints, order_geo_opt, order_cent(3), &
                                num_geo_opt, num_tot_opt, scatt1_cints)
            if (order_geo_bra==0) then
              ! only \var(order_geo_opt)/=0
              if (order_geo_ket==0) then
                ! scatters geometric derivatives on bra and ket centers
                do nopt = 1, num_opt
                  igeo = 0
                  do iopt = 1, num_tot_opt
                    do ibra = 1, dim_geo_bra
                      do iket = 1, dim_geo_ket
                        igeo = igeo+1
                        total_cints(:,1,1,:,igeo,nopt) &
                          = scatt1_cints(:,ibra,iket,:,iopt,nopt)
                      end do
                    end do
                  end do
                end do
                deallocate(scatt1_cints)
              ! only \var(order_geo_bra)==0
              else
                ! allocates the integrals after scattering geometric derivatives on bra center
                allocate(scatt2_cints(dim_cints,dim_geo_ket,num_geo_opt, &
                                      dim_geo_bra,num_tot_opt,num_opt), stat=ierr)
                if (ierr/=0)                                                    &
                  call error_stop("geom_part_one_scatter",                      &
                                  "failed to allocate scatt2_cints/3/KBC/X0X",  &
                                  dim_cints*dim_geo_ket*num_geo_opt*dim_geo_bra &
                                  *num_tot_opt*num_opt)
                ! scatters geometric derivatives on bra center
                do nopt = 1, num_opt
                  do iopt = 1, num_tot_opt
                    do iket = 1, dim_geo_ket
                      do ibra = 1, dim_geo_bra
                        scatt2_cints(:,iket,:,ibra,iopt,nopt) &
                          = scatt1_cints(:,ibra,iket,:,iopt,nopt)
                      end do
                    end do
                  end do
                end do
                deallocate(scatt1_cints)
                ! sets the number of total geometric derivatives on ket center
                num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on ket center
                call scatter_single(dim_cints, dim_geo_ket, num_geo_opt,        &
                                    dim_geo_bra*num_tot_opt*num_opt,            &
                                    scatt2_cints, order_geo_ket, order_cent(1), &
                                    num_geo_ket, num_tot_ket, total_cints)
                deallocate(scatt2_cints)
              end if
            ! at least \var(order_geo_opt)/=0 and \var(order_geo_bra)/=0
            else 
              ! sets the number of total geometric derivatives on bra center
              num_tot_bra = (order_cent(2)+1)*(order_cent(2)+2)/2
              num_tot_cent2 = num_tot_bra*num_tot_opt
              ! allocates the integrals after scattering geometric derivatives on bra center
              allocate(scatt2_cints(dim_cints,num_geo_bra,dim_geo_ket,num_geo_opt, &
                                    num_tot_cent2,num_opt), stat=ierr)
              if (ierr/=0)                                                    &
                call error_stop("geom_part_one_scatter",                      &
                                "failed to allocate scatt2_cints/3/KBC/?XX",  &
                                dim_cints*num_geo_bra*dim_geo_ket*num_geo_opt &
                                *num_tot_cent2*num_opt)
              ! scatters geometric derivatives on bra center
              call scatter_single(dim_cints, dim_geo_bra, dim_geo_ket*num_geo_opt,  &
                                  num_tot_opt*num_opt, scatt1_cints, order_geo_bra, &
                                  order_cent(2), num_geo_bra, num_tot_bra, scatt2_cints)
              deallocate(scatt1_cints)
              ! only \var(order_geo_ket)==0
              if (order_geo_ket==0) then
                ! scatters geometric derivatives on ket center
                do nopt = 1, num_opt
                  igeo = 0
                  do ibra = 1, num_tot_cent2
                    do iket = 1, dim_geo_ket
                      igeo = igeo+1
                      total_cints(:,:,1,:,igeo,nopt) &
                        = scatt2_cints(:,:,iket,:,ibra,nopt)
                    end do
                  end do
                end do
              ! none of the orders is zero
              else
                ! sets the number of total geometric derivatives on ket center
                num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on ket center
                call scatter_single(dim_cints*num_geo_bra, dim_geo_ket,         &
                                    num_geo_opt, num_tot_cent2*num_opt,         &
                                    scatt2_cints, order_geo_ket, order_cent(1), &
                                    num_geo_ket, num_tot_ket, total_cints)
              end if
              deallocate(scatt2_cints)
            end if
          end if
        ! the second differentiated center is operator center, and the third is bra center
        case default
          if (order_geo_bra==0) then
            ! \var(order_geo_bra)==0 and \var(order_geo_opt)==0
            if (order_geo_opt==0) then
              ! no partial geometric derivatives
              if (order_geo_ket==0) then
                ! switches the derivatives on bra, ket and operators centers
                do nopt = 1, num_opt
                  igeo = 0
                  do ibra = 1, dim_geo_bra
                    do iopt = 1, dim_geo_opt
                      do iket = 1, dim_geo_ket
                        igeo = igeo+1
                        total_cints(:,1,1,1,igeo,nopt) &
                          = part_cints(:,ibra,iket,iopt,nopt)
                      end do
                    end do
                  end do
                end do
              ! only \var(order_geo_ket)/=0
              else
                ! switches the derivatives on bra, ket and operator centers
                allocate(scatt1_cints(dim_cints,dim_geo_ket,dim_geo_opt, &
                                      dim_geo_bra,num_opt,1), stat=ierr)
                if (ierr/=0)                                                   &
                  call error_stop("geom_part_one_scatter",                     &
                                  "failed to allocate scatt1_cints/3/KCB/X00", &
                                  dim_cints*dim_geo_ket*dim_geo_opt*dim_geo_bra*num_opt)
                do nopt = 1, num_opt
                  do iopt = 1, dim_geo_opt
                    do iket = 1, dim_geo_ket
                      do ibra = 1, dim_geo_bra
                        scatt1_cints(:,iket,iopt,ibra,nopt,1) &
                          = part_cints(:,ibra,iket,iopt,nopt)
                      end do
                    end do
                  end do
                end do
                ! sets the number of total geometric derivatives on ket center
                num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
                call scatter_single(dim_cints, dim_geo_ket, 1,                  &
                                    dim_geo_opt*dim_geo_bra*num_opt,            &
                                    scatt1_cints, order_geo_ket, order_cent(1), &
                                    num_geo_ket, num_tot_ket, total_cints)
                deallocate(scatt1_cints)
              end if
            ! \var(order_geo_bra)==0 but \var(order_geo_opt)/=0
            else
              ! switches the derivatives on bra, ket and operator centers
              allocate(scatt1_cints(dim_cints,dim_geo_ket,dim_geo_opt, &
                                    dim_geo_bra,num_opt,1), stat=ierr)
              if (ierr/=0)                                                   &
                call error_stop("geom_part_one_scatter",                     &
                                "failed to allocate scatt1_cints/3/KCB/?X0", &
                                dim_cints*dim_geo_ket*dim_geo_opt*dim_geo_bra*num_opt)
              do nopt = 1, num_opt
                do iopt = 1, dim_geo_opt
                  do iket = 1, dim_geo_ket
                    do ibra = 1, dim_geo_bra
                      scatt1_cints(:,iket,iopt,ibra,nopt,1) &
                        = part_cints(:,ibra,iket,iopt,nopt)
                    end do
                  end do
                end do
              end do
              ! sets the number of total geometric derivatives on operator center
              num_tot_opt = (order_cent(2)+1)*(order_cent(2)+2)/2
              ! allocates the integrals after scattering geometric derivatives on operator center
              allocate(scatt2_cints(dim_cints,dim_geo_ket,num_geo_opt, &
                                    num_tot_opt,dim_geo_bra,num_opt), stat=ierr)
              if (ierr/=0)                                                    &
                call error_stop("geom_part_one_scatter",                      &
                                "failed to allocate scatt2_cints/3/KCB/?X0",  &
                                dim_cints*dim_geo_ket*num_geo_opt*num_tot_opt &
                                *dim_geo_bra*num_opt)
              ! scatters geometric derivatives on operator center
              call scatter_single(dim_cints*dim_geo_ket, dim_geo_opt, 1,            &
                                  dim_geo_bra*num_opt, scatt1_cints, order_geo_opt, &
                                  order_cent(2), num_geo_opt, num_tot_opt, scatt2_cints)
              deallocate(scatt1_cints)
              ! only \var(order_geo_opt)/=0
              if (order_geo_ket==0) then
                ! scatters geometric derivatives on ket center
                do nopt = 1, num_opt
                  igeo = 0
                  do ibra = 1, dim_geo_bra
                    do iopt = 1, num_tot_opt
                      do iket = 1, dim_geo_ket
                        igeo = igeo+1
                        total_cints(:,1,1,:,igeo,nopt) &
                          = scatt2_cints(:,iket,:,iopt,ibra,nopt)
                      end do
                    end do
                  end do
                end do
              ! only \var(order_geo_bra)==0
              else
                ! sets the number of total geometric derivatives on ket center
                num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on ket center
                call scatter_single(dim_cints, dim_geo_ket, num_geo_opt,        &
                                    num_tot_opt*dim_geo_bra*num_opt,            &
                                    scatt2_cints, order_geo_ket, order_cent(1), &
                                    num_geo_ket, num_tot_ket, total_cints)
              end if
              deallocate(scatt2_cints)
            end if
          ! at least \var(order_geo_bra)/=0
          else
            ! sets the number of total geometric derivatives on bra center
            num_tot_bra = (order_cent(3)+1)*(order_cent(3)+2)/2
            if (order_geo_opt==0) then
              ! only \var(order_geo_bra)/=0
              if (order_geo_ket==0) then
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, dim_geo_ket*dim_geo_opt,   &
                                    num_opt, part_cints, order_geo_bra, order_cent(3), &
                                    num_geo_bra, num_tot_bra, total_cints)
              ! only \var(order_geo_opt)==0
              else
                ! allocates the integrals after scattering geometric derivatives on bra center
                allocate(scatt1_cints(dim_cints,num_geo_bra,dim_geo_ket, &
                                      dim_geo_opt,num_tot_bra,num_opt), stat=ierr)
                if (ierr/=0)                                                    &
                  call error_stop("geom_part_one_scatter",                      &
                                  "failed to allocate scatt1_cints/3/KCB/X0X",  &
                                  dim_cints*num_geo_bra*dim_geo_ket*dim_geo_opt &
                                  *num_tot_bra*num_opt)
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, dim_geo_ket*dim_geo_opt,   &
                                    num_opt, part_cints, order_geo_bra, order_cent(3), &
                                    num_geo_bra, num_tot_bra, scatt1_cints)
                ! sets the number of total geometric derivatives on ket center
                num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on ket center
                call scatter_single(dim_cints*num_geo_bra, dim_geo_ket, 1,      &
                                    dim_geo_opt*num_tot_bra*num_opt,            &
                                    scatt1_cints, order_geo_ket, order_cent(1), &
                                    num_geo_ket, num_tot_ket, total_cints)
                deallocate(scatt1_cints)
              end if
            ! at least \var(order_geo_bra)/=0 and \var(order_geo_opt)/=0
            else 
              ! allocates the integrals after scattering geometric derivatives on bra center
              allocate(scatt1_cints(dim_cints,num_geo_bra,dim_geo_ket, &
                                    dim_geo_opt,num_tot_bra,num_opt), stat=ierr)
              if (ierr/=0)                                                    &
                call error_stop("geom_part_one_scatter",                      &
                                "failed to allocate scatt1_cints/3/KCB/?XX",  &
                                dim_cints*num_geo_bra*dim_geo_ket*dim_geo_opt &
                                *num_tot_bra*num_opt)
              ! scatters geometric derivatives on bra center
              call scatter_single(dim_cints, dim_geo_bra, dim_geo_ket*dim_geo_opt,   &
                                  num_opt, part_cints, order_geo_bra, order_cent(3), &
                                  num_geo_bra, num_tot_bra, scatt1_cints)
              ! sets the number of total geometric derivatives on operator center
              num_tot_opt = (order_cent(2)+1)*(order_cent(2)+2)/2
              num_tot_cent2 = num_tot_opt*num_tot_bra
              ! allocates the integrals after scattering geometric derivatives on operator center
              allocate(scatt2_cints(dim_cints,num_geo_bra,dim_geo_ket,num_geo_opt, &
                                    num_tot_cent2,num_opt), stat=ierr)
              if (ierr/=0)                                                    &
                call error_stop("geom_part_one_scatter",                      &
                                "failed to allocate scatt2_cints/3/KCB/?XX",  &
                                dim_cints*num_geo_bra*dim_geo_ket*num_geo_opt &
                                *num_tot_cent2*num_opt)
              ! scatters geometric derivatives on operator center
              call scatter_single(dim_cints*num_geo_bra*dim_geo_ket,          &
                                  dim_geo_opt, 1, num_tot_bra*num_opt,        &
                                  scatt1_cints, order_geo_opt, order_cent(2), &
                                  num_geo_opt, num_tot_opt, scatt2_cints)
              deallocate(scatt1_cints)
              ! only \var(order_geo_ket)==0
              if (order_geo_ket==0) then
                ! scatters geometric derivatives on ket center
                do nopt = 1, num_opt
                  igeo = 0
                  do iopt = 1, num_tot_cent2
                    do iket = 1, dim_geo_ket
                      igeo = igeo+1
                      total_cints(:,:,1,:,igeo,nopt) &
                        = scatt2_cints(:,:,iket,:,iopt,nopt)
                    end do
                  end do
                end do
              ! none of the orders is zero
              else
                ! sets the number of total geometric derivatives on ket center
                num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on ket center
                call scatter_single(dim_cints*num_geo_bra, dim_geo_ket,         &
                                    num_geo_opt, num_tot_cent2*num_opt,         &
                                    scatt2_cints, order_geo_ket, order_cent(1), &
                                    num_geo_ket, num_tot_ket, total_cints)
              end if
              deallocate(scatt2_cints)
            end if
          end if
        end select
      ! the first differentiated center is operator center
      case default
        select case(seq_part_geo(2))
        ! the second differentiated center is bra center, and the third is ket center
        case(TAG_BRA)
          if (order_geo_ket==0) then
            ! \var(order_geo_ket)==0 and \var(order_geo_bra)==0
            if (order_geo_bra==0) then
              ! no partial geometric derivatives
              if (order_geo_opt==0) then
                ! switches the derivatives on bra and ket, operator centers
                do nopt = 1, num_opt
                  igeo = 0
                  do iket = 1, dim_geo_ket
                    do ibra = 1, dim_geo_bra
                      do iopt = 1, dim_geo_opt
                        igeo = igeo+1
                        total_cints(:,1,1,1,igeo,nopt) &
                          = part_cints(:,ibra,iket,iopt,nopt)
                      end do
                    end do
                  end do
                end do
              ! only \var(order_geo_opt)/=0
              else
                ! switches the derivatives on bra and ket, operator centers
                allocate(scatt1_cints(dim_cints,dim_geo_opt,dim_geo_bra, &
                                      dim_geo_ket,num_opt,1), stat=ierr)
                if (ierr/=0)                                                   &
                  call error_stop("geom_part_one_scatter",                     &
                                  "failed to allocate scatt1_cints/3/CBK/X00", &
                                  dim_cints*dim_geo_opt*dim_geo_bra*dim_geo_ket*num_opt)
                do nopt = 1, num_opt
                  do iopt = 1, dim_geo_opt
                    do iket = 1, dim_geo_ket
                      do ibra = 1, dim_geo_bra
                        scatt1_cints(:,iopt,ibra,iket,nopt,1) &
                          = part_cints(:,ibra,iket,iopt,nopt)
                      end do
                    end do
                  end do
                end do
                ! sets the number of total geometric derivatives on operator center
                num_tot_opt = (order_cent(1)+1)*(order_cent(1)+2)/2
                call scatter_single(dim_cints, dim_geo_opt, 1,                  &
                                    dim_geo_bra*dim_geo_ket*num_opt,            &
                                    scatt1_cints, order_geo_opt, order_cent(1), &
                                    num_geo_opt, num_tot_opt, total_cints)
                deallocate(scatt1_cints)
              end if
            ! \var(order_geo_ket)==0 but \var(order_geo_bra)/=0
            else
              ! switches the derivatives on ket and operator centers
              allocate(scatt1_cints(dim_cints,dim_geo_bra,dim_geo_opt, &
                                    dim_geo_ket,num_opt,1), stat=ierr)
              if (ierr/=0)                                                   &
                call error_stop("geom_part_one_scatter",                     &
                                "failed to allocate scatt1_cints/3/CBK/?X0", &
                                dim_cints*dim_geo_bra*dim_geo_opt*dim_geo_ket*num_opt)
              do nopt = 1, num_opt
                do iopt = 1, dim_geo_opt
                  do iket = 1, dim_geo_ket
                    do ibra = 1, dim_geo_bra
                      scatt1_cints(:,ibra,iopt,iket,nopt,1) &
                        = part_cints(:,ibra,iket,iopt,nopt)
                    end do
                  end do
                end do
              end do
              ! sets the number of total geometric derivatives on bra center
              num_tot_bra = (order_cent(2)+1)*(order_cent(2)+2)/2
              ! only \var(order_geo_bra)/=0
              if (order_geo_opt==0) then
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, dim_geo_opt, &
                                    dim_geo_ket*num_opt, scatt1_cints,   &
                                    order_geo_bra, order_cent(2),        &
                                    num_geo_bra, num_tot_bra, total_cints)
                deallocate(scatt1_cints)
              ! only \var(order_geo_ket)==0
              else
                ! allocates the integrals after scattering geometric derivatives on bra center
                allocate(scatt2_cints(dim_cints,num_geo_bra,dim_geo_opt, &
                                      num_tot_bra,dim_geo_ket,num_opt), stat=ierr)
                if (ierr/=0)                                                    &
                  call error_stop("geom_part_one_scatter",                      &
                                  "failed to allocate scatt2_cints/3/CBK/XX0",  &
                                  dim_cints*num_geo_bra*dim_geo_opt*num_tot_bra &
                                  *dim_geo_ket*num_opt)
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, dim_geo_opt, &
                                    dim_geo_ket*num_opt, scatt1_cints,   &
                                    order_geo_bra, order_cent(2),        &
                                    num_geo_bra, num_tot_bra, scatt2_cints)
                deallocate(scatt1_cints)
                ! sets the number of total geometric derivatives on operator center
                num_tot_opt = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on operator center
                call scatter_single(dim_cints*num_geo_bra, dim_geo_opt,         &
                                    1, num_tot_bra*dim_geo_ket*num_opt,         &
                                    scatt2_cints, order_geo_opt, order_cent(1), &
                                    num_geo_opt, num_tot_opt, total_cints)
                deallocate(scatt2_cints)
              end if
            end if
          ! at least \var(order_geo_ket)/=0
          else
            ! sets the number of total geometric derivatives on ket center
            num_tot_ket = (order_cent(3)+1)*(order_cent(3)+2)/2
            ! allocates the integrals after scattering geometric derivatives on ket center
            allocate(scatt1_cints(dim_cints,dim_geo_bra,num_geo_ket, &
                                  dim_geo_opt,num_tot_ket,num_opt), stat=ierr)
            if (ierr/=0)                                                    &
              call error_stop("geom_part_one_scatter",                      &
                              "failed to allocate scatt1_cints/3/CBK/??X",  &
                              dim_cints*dim_geo_bra*num_geo_ket*dim_geo_opt &
                              *num_tot_ket*num_opt)
            ! scatters geometric derivatives on ket center
            call scatter_single(dim_cints*dim_geo_bra, dim_geo_ket, dim_geo_opt,   &
                                num_opt, part_cints, order_geo_ket, order_cent(3), &
                                num_geo_ket, num_tot_ket, scatt1_cints)
            if (order_geo_bra==0) then
              ! only \var(order_geo_ket)/=0
              if (order_geo_opt==0) then
                ! scatters geometric derivatives on bra and operator centers
                do nopt = 1, num_opt
                  igeo = 0
                  do iket = 1, num_tot_ket
                    do ibra = 1, dim_geo_bra
                      do iopt = 1, dim_geo_opt
                        igeo = igeo+1
                        total_cints(:,1,:,1,igeo,nopt) &
                          = scatt1_cints(:,ibra,:,iopt,iket,nopt)
                      end do
                    end do
                  end do
                end do
                deallocate(scatt1_cints)
              ! only \var(order_geo_bra)==0
              else
                ! allocates the integrals after scattering geometric derivatives on bra center
                allocate(scatt2_cints(dim_cints,num_geo_ket,dim_geo_opt, &
                                      dim_geo_bra,num_tot_ket,num_opt), stat=ierr)
                if (ierr/=0)                                                    &
                  call error_stop("geom_part_one_scatter",                      &
                                  "failed to allocate scatt2_cints/3/CBK/X0X",  &
                                  dim_cints*num_geo_ket*dim_geo_opt*dim_geo_bra &
                                  *num_tot_ket*num_opt)
                ! scatters geometric derivatives on bra center
                do nopt = 1, num_opt
                  do iket = 1, num_tot_ket
                    do iopt = 1, dim_geo_opt
                      do ibra = 1, dim_geo_bra
                        scatt2_cints(:,:,iopt,ibra,iket,nopt) &
                          = scatt1_cints(:,ibra,:,iopt,iket,nopt)
                      end do
                    end do
                  end do
                end do
                deallocate(scatt1_cints)
                ! sets the number of total geometric derivatives on operator center
                num_tot_opt = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on operator center
                call scatter_single(dim_cints*num_geo_ket, dim_geo_opt,         &
                                    1, dim_geo_bra*num_tot_ket*num_opt,         &
                                    scatt2_cints, order_geo_opt, order_cent(1), &
                                    num_geo_opt, num_tot_opt, total_cints)
                deallocate(scatt2_cints)
              end if
            ! at least \var(order_geo_ket)/=0 and \var(order_geo_bra)/=0
            else 
              ! sets the number of total geometric derivatives on bra center
              num_tot_bra = (order_cent(2)+1)*(order_cent(2)+2)/2
              ! only \var(order_geo_opt)==0
              if (order_geo_opt==0) then
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, num_geo_ket*dim_geo_opt,  &
                                    num_tot_ket*num_opt, scatt1_cints, order_geo_bra, &
                                    order_cent(2), num_geo_bra, num_tot_bra, total_cints)
                deallocate(scatt1_cints)
              ! none of the orders is zero
              else
                num_tot_cent2 = num_tot_bra*num_tot_ket
                ! allocates the integrals after scattering geometric derivatives on bra center
                allocate(scatt2_cints(dim_cints,num_geo_bra,num_geo_ket,dim_geo_opt, &
                                      num_tot_cent2,num_opt), stat=ierr)
                if (ierr/=0)                                                    &
                  call error_stop("geom_part_one_scatter",                      &
                                  "failed to allocate scatt2_cints/3/CBK/XXX",  &
                                  dim_cints*num_geo_bra*num_geo_ket*dim_geo_opt &
                                  *num_tot_cent2*num_opt)
                ! scatters geometric derivatives on bra center
                call scatter_single(dim_cints, dim_geo_bra, num_geo_ket*dim_geo_opt,  &
                                    num_tot_ket*num_opt, scatt1_cints, order_geo_bra, &
                                    order_cent(2), num_geo_bra, num_tot_bra, scatt2_cints)
                deallocate(scatt1_cints)
                ! sets the number of total geometric derivatives on operator center
                num_tot_opt = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on operator center
                call scatter_single(dim_cints*num_geo_bra*num_geo_ket,          &
                                    dim_geo_opt, 1, num_tot_cent2*num_opt,      &
                                    scatt2_cints, order_geo_opt, order_cent(1), &
                                    num_geo_opt, num_tot_opt, total_cints)
                deallocate(scatt2_cints)
              end if
            end if
          end if
        ! the second differentiated center is ket center, and the third is bra center
        case default
          if (order_geo_bra==0) then
            ! \var(order_geo_bra)==0 and \var(order_geo_ket)==0
            if (order_geo_ket==0) then
              ! no partial geometric derivatives
              if (order_geo_opt==0) then
                ! switches the derivatives on bra and operator centers
                do nopt = 1, num_opt
                  igeo = 0
                  do ibra = 1, dim_geo_bra
                    do iket = 1, dim_geo_ket
                      do iopt = 1, dim_geo_opt
                        igeo = igeo+1
                        total_cints(:,1,1,1,igeo,nopt) &
                          = part_cints(:,ibra,iket,iopt,nopt)
                      end do
                    end do
                  end do
                end do
              ! only \var(order_geo_opt)/=0
              else
                ! switches the derivatives on bra and operator centers
                allocate(scatt1_cints(dim_cints,dim_geo_opt,dim_geo_ket, &
                                      dim_geo_bra,num_opt,1), stat=ierr)
                if (ierr/=0)                                                   &
                  call error_stop("geom_part_one_scatter",                     &
                                  "failed to allocate scatt1_cints/3/CKB/X00", &
                                  dim_cints*dim_geo_opt*dim_geo_ket*dim_geo_bra*num_opt)
                do nopt = 1, num_opt
                  do iopt = 1, dim_geo_opt
                    do iket = 1, dim_geo_ket
                      do ibra = 1, dim_geo_bra
                        scatt1_cints(:,iopt,iket,ibra,nopt,1) &
                          = part_cints(:,ibra,iket,iopt,nopt)
                      end do
                    end do
                  end do
                end do
                ! sets the number of total geometric derivatives on operator center
                num_tot_opt = (order_cent(1)+1)*(order_cent(1)+2)/2
                call scatter_single(dim_cints, dim_geo_opt, 1,                  &
                                    dim_geo_ket*dim_geo_bra*num_opt,            &
                                    scatt1_cints, order_geo_opt, order_cent(1), &
                                    num_geo_opt, num_tot_opt, total_cints)
                deallocate(scatt1_cints)
              end if
            ! \var(order_geo_bra)==0 but \var(order_geo_ket)/=0
            else
              ! switches the derivatives on bra, ket and operator centers
              allocate(scatt1_cints(dim_cints,dim_geo_ket,dim_geo_opt, &
                                    dim_geo_bra,num_opt,1), stat=ierr)
              if (ierr/=0)                                                   &
                call error_stop("geom_part_one_scatter",                     &
                                "failed to allocate scatt1_cints/3/CKB/?X0", &
                                dim_cints*dim_geo_ket*dim_geo_opt*dim_geo_bra*num_opt)
              do nopt = 1, num_opt
                do iopt = 1, dim_geo_opt
                  do iket = 1, dim_geo_ket
                    do ibra = 1, dim_geo_bra
                      scatt1_cints(:,iket,iopt,ibra,nopt,1) &
                        = part_cints(:,ibra,iket,iopt,nopt)
                    end do
                  end do
                end do
              end do
              ! sets the number of total geometric derivatives on ket center
              num_tot_ket = (order_cent(2)+1)*(order_cent(2)+2)/2
              ! only \var(order_geo_ket)/=0
              if (order_geo_opt==0) then
                ! scatters geometric derivatives on ket center
                call scatter_single(dim_cints, dim_geo_ket, dim_geo_opt, &
                                    dim_geo_bra*num_opt, scatt1_cints,   &
                                    order_geo_ket, order_cent(2),        &
                                    num_geo_ket, num_tot_ket, total_cints)
                deallocate(scatt1_cints)
              ! only \var(order_geo_bra)==0
              else
                ! allocates the integrals after scattering geometric derivatives on ket center
                allocate(scatt2_cints(dim_cints,num_geo_ket,dim_geo_opt, &
                                      num_tot_ket,dim_geo_bra,num_opt), stat=ierr)
                if (ierr/=0)                                                    &
                  call error_stop("geom_part_one_scatter",                      &
                                  "failed to allocate scatt2_cints/3/CKB/XX0",  &
                                  dim_cints*num_geo_ket*dim_geo_opt*num_tot_ket &
                                  *dim_geo_bra*num_opt)
                ! scatters geometric derivatives on ket center
                call scatter_single(dim_cints, dim_geo_ket, dim_geo_opt, &
                                    dim_geo_bra*num_opt, scatt1_cints,   &
                                    order_geo_ket, order_cent(2),        &
                                    num_geo_ket, num_tot_ket, scatt2_cints)
                deallocate(scatt1_cints)
                ! sets the number of total geometric derivatives on operator center
                num_tot_opt = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on operator center
                call scatter_single(dim_cints*num_geo_ket, dim_geo_opt,         &
                                    1, num_tot_ket*dim_geo_bra*num_opt,         &
                                    scatt2_cints, order_geo_opt, order_cent(1), &
                                    num_geo_opt, num_tot_opt, total_cints)
                deallocate(scatt2_cints)
              end if
            end if
          ! at least \var(order_geo_bra)/=0
          else
            ! sets the number of total geometric derivatives on bra center
            num_tot_bra = (order_cent(3)+1)*(order_cent(3)+2)/2
            ! allocates the integrals after scattering geometric derivatives on bra center
            allocate(scatt1_cints(dim_cints,num_geo_bra,dim_geo_ket, &
                                  dim_geo_opt,num_tot_bra,num_opt), stat=ierr)
            if (ierr/=0)                                                    &
              call error_stop("geom_part_one_scatter",                      &
                              "failed to allocate scatt1_cints/3/CKB/??X",  &
                              dim_cints*num_geo_bra*dim_geo_ket*dim_geo_opt &
                              *num_tot_bra*num_opt)
            ! scatters geometric derivatives on bra center
            call scatter_single(dim_cints, dim_geo_bra, dim_geo_ket*dim_geo_opt,   &
                                num_opt, part_cints, order_geo_bra, order_cent(3), &
                                num_geo_bra, num_tot_bra, scatt1_cints)
            if (order_geo_ket==0) then
              ! only \var(order_geo_bra)/=0
              if (order_geo_opt==0) then
                ! scatters geometric derivatives on ket and operator centers
                do nopt = 1, num_opt
                  igeo = 0
                  do ibra = 1, num_tot_bra
                    do iket = 1, dim_geo_ket
                      do iopt = 1, dim_geo_opt
                        igeo = igeo+1
                        total_cints(:,:,1,1,igeo,nopt) &
                          = scatt1_cints(:,:,iket,iopt,ibra,nopt)
                      end do
                    end do
                  end do
                end do
                deallocate(scatt1_cints)
              ! only \var(order_geo_ket)==0
              else
                ! allocates the integrals after scattering geometric derivatives on ket center
                allocate(scatt2_cints(dim_cints,num_geo_bra,dim_geo_opt, &
                                      dim_geo_ket,num_tot_bra,num_opt), stat=ierr)
                if (ierr/=0)                                                    &
                  call error_stop("geom_part_one_scatter",                      &
                                  "failed to allocate scatt2_cints/3/CKB/X0X",  &
                                  dim_cints*num_geo_bra*dim_geo_opt*dim_geo_ket &
                                  *num_tot_bra*num_opt)
                ! scatters geometric derivatives on ket center
                do nopt = 1, num_opt
                  do ibra = 1, num_tot_bra
                    do iopt = 1, dim_geo_opt
                      do iket = 1, dim_geo_ket
                        scatt2_cints(:,:,iopt,iket,ibra,nopt) &
                          = scatt1_cints(:,:,iket,iopt,ibra,nopt)
                      end do
                    end do
                  end do
                end do
                deallocate(scatt1_cints)
                ! sets the number of total geometric derivatives on operator center
                num_tot_opt = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on operator center
                call scatter_single(dim_cints*num_geo_bra, dim_geo_opt,         &
                                    1, dim_geo_ket*num_tot_bra*num_opt,         &
                                    scatt2_cints, order_geo_opt, order_cent(1), &
                                    num_geo_opt, num_tot_opt, total_cints)
                deallocate(scatt2_cints)
              end if
            ! at least \var(order_geo_bra)/=0 and \var(order_geo_ket)/=0
            else 
              ! sets the number of total geometric derivatives on ket center
              num_tot_ket = (order_cent(2)+1)*(order_cent(2)+2)/2
              ! only \var(order_geo_opt)==0
              if (order_geo_opt==0) then
                ! scatters geometric derivatives on ket center
                call scatter_single(dim_cints*num_geo_bra, dim_geo_ket,         &
                                    dim_geo_opt, num_tot_bra*num_opt,           &
                                    scatt1_cints, order_geo_ket, order_cent(2), &
                                    num_geo_ket, num_tot_ket, total_cints)
                deallocate(scatt1_cints)
              ! none of the orders is zero
              else
                num_tot_cent2 = num_tot_ket*num_tot_bra
                ! allocates the integrals after scattering geometric derivatives on ket center
                allocate(scatt2_cints(dim_cints,num_geo_bra,num_geo_ket,dim_geo_opt, &
                                      num_tot_cent2,num_opt), stat=ierr)
                if (ierr/=0)                                                    &
                  call error_stop("geom_part_one_scatter",                      &
                                  "failed to allocate scatt2_cints/3/CKB/XXX",  &
                                  dim_cints*num_geo_bra*num_geo_ket*dim_geo_opt &
                                  *num_tot_cent2*num_opt)
                ! scatters geometric derivatives on ket center
                call scatter_single(dim_cints*num_geo_bra, dim_geo_ket,         &
                                    dim_geo_opt, num_tot_bra*num_opt,           &
                                    scatt1_cints, order_geo_ket, order_cent(2), &
                                    num_geo_ket, num_tot_ket, scatt2_cints)
                deallocate(scatt1_cints)
                ! sets the number of total geometric derivatives on operator center
                num_tot_opt = (order_cent(1)+1)*(order_cent(1)+2)/2
                ! scatters geometric derivatives on operator center
                call scatter_single(dim_cints*num_geo_bra*num_geo_ket,          &
                                    dim_geo_opt, 1, num_tot_cent2*num_opt,      &
                                    scatt2_cints, order_geo_opt, order_cent(1), &
                                    num_geo_opt, num_tot_opt, total_cints)
                deallocate(scatt2_cints)
              end if
            end if
          end if
        end select
      end select
    case default
      call error_stop("geom_part_one_scatter", "invalid num_cents", num_cents)
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "geom_part_one_scatter", STDOUT)
#endif
    return
  end subroutine geom_part_one_scatter
