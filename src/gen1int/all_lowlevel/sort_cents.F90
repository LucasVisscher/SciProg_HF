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
!!  This file contains subroutines sorting the given centers.
!!
!!  2011-07-02, Bin Gao:
!!  * first version

#include "stdout.h"
#include "max_idx_non.h"

  !> \brief sorts the indices of centers using insertion sort algorithm, and
  !>        finds the unique centers and corresponding number of identical centers
  !> \author Bin Gao
  !> \date 2011-07-02
  !> \param num_cents is the number of centers
  !> \return idx_cent contains the indices of the centers, numbers less than 1 indicate
  !>         non-atomic centers, sorted in ascending order on exit
  !> \return tag_cent contains the tags marking the places of centers before sorting
  !> \return num_non_cents is the number of unique non-atomic centers
  !> \return id_non_cent contains the places of unique non-atomic centers
  !> \return num_non_iden contains the number of identical non-atomic centers
  !> \return num_atom_cents is the number of unique atomic centers
  !> \return id_atom_cent contains the places of unique atomic centers
  !> \return num_atom_iden contains the number of identical atomic centers
  subroutine sort_gen_cents(num_cents, idx_cent, tag_cent,            &
                            num_non_cents, id_non_cent, num_non_iden, &
                            num_atom_cents, id_atom_cent, num_atom_iden)
    use xkind
    implicit none
    integer, intent(in) :: num_cents
    integer, intent(inout) :: idx_cent(num_cents)
    integer, intent(inout) :: tag_cent(num_cents)
    integer, intent(out) :: num_non_cents
    integer, intent(out) :: id_non_cent(num_cents)
    integer, intent(out) :: num_non_iden(num_cents)
    integer, intent(out) :: num_atom_cents
    integer, intent(out) :: id_atom_cent(num_cents)
    integer, intent(out) :: num_atom_iden(num_cents)
!f2py intent(hide) :: num_cents
!f2py intent(inout) :: idx_cent
!f2py intent(inout) :: tag_cent
!f2py depend(num_cents) :: tag_cent
!f2py intent(out) :: num_non_cents
!f2py intent(out) :: id_non_cent
!f2py depend(num_cents) :: id_non_cent
!f2py intent(out) :: num_non_iden
!f2py depend(num_cents) :: num_non_iden
!f2py intent(out) :: num_atom_cents
!f2py intent(out) :: id_atom_cent
!f2py depend(num_cents) :: id_atom_cent
!f2py intent(out) :: num_atom_iden
!f2py depend(num_cents) :: num_atom_iden
    integer ircd, jrcd, krcd  !incremental recorders
    integer idx_ins           !inserting index
    integer tag_ins           !tag of inserting index
    logical not_ins           !not found right place to insert the index
#if defined(XTIME)
    real(REALK) curr_time     !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
#if defined(DEBUG)
    write(STDOUT,100) "indices before sorting:", idx_cent
#endif
    ! checks if the first one is atomic center
    if (idx_cent(1)<=MAX_IDX_NON) then
      num_non_cents = 1
      num_atom_cents = 0
    else
      num_non_cents = 0
      num_atom_cents = 1
    end if
    ! loops over the indices
    do ircd = 2, num_cents
      idx_ins = idx_cent(ircd)
      tag_ins = tag_cent(ircd)
      jrcd = ircd-1
      not_ins = .true.
      ! when an index is less than the previous one,
      ! loops down to insert it at its right place indicated by \var(jrcd)+1
      do while (not_ins)
        if (idx_ins<idx_cent(jrcd)) then
          idx_cent(jrcd+1) = idx_cent(jrcd)
          tag_cent(jrcd+1) = tag_cent(jrcd)
          jrcd = jrcd-1
          if (jrcd<1) then
            ! found another unique non-atomic center \var(idx_ins)
            if (idx_ins<=MAX_IDX_NON) then
              num_non_cents = num_non_cents+1
            ! this is an atomic center
            else
              num_atom_cents = num_atom_cents+1
            end if
            not_ins = .false.
          end if
        else
          ! found another unique center \var(idx_ins)
          if (idx_ins>idx_cent(jrcd)) then
            ! found another unique non-atomic center \var(idx_ins)
            if (idx_ins<=MAX_IDX_NON) then
              num_non_cents = num_non_cents+1
            ! this is an atomic center
            else
              num_atom_cents = num_atom_cents+1
            end if
          end if
          not_ins = .false.
        end if
      end do
      idx_cent(jrcd+1) = idx_ins
      tag_cent(jrcd+1) = tag_ins
    end do
    ! tries to find the positions of unique centers, and the number of identical centers
    ! non-atomic centers
    if (num_non_cents>0) then
      id_non_cent(1) = 1
      num_non_iden(1) = 1
      jrcd = 1
      do ircd = 2, num_cents
        ! exists for atomic center
        if (idx_cent(ircd)>MAX_IDX_NON) then
          exit
        else
          ! found another unique non-atomic center
          if (idx_cent(ircd)/=idx_cent(id_non_cent(jrcd))) then
            jrcd = jrcd+1
            id_non_cent(jrcd) = ircd
            num_non_iden(jrcd) = 1
          else
            num_non_iden(jrcd) = num_non_iden(jrcd)+1
          end if
        end if
      end do
    else
      ircd = 1
    end if
    ! atomic-centers
    id_atom_cent(1) = ircd
    num_atom_iden(1) = 1
    jrcd = 1
    do krcd = ircd+1, num_cents
      ! found another unique atomic center
      if (idx_cent(krcd)/=idx_cent(id_atom_cent(jrcd))) then
        jrcd = jrcd+1
        id_atom_cent(jrcd) = krcd
        num_atom_iden(jrcd) = 1
      else
        num_atom_iden(jrcd) = num_atom_iden(jrcd)+1
      end if
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "sort_gen_cents", STDOUT)
#endif
#if defined(DEBUG)
    write(STDOUT,100) "indices after sorting: ", idx_cent
    write(STDOUT,100) "previous places:       ", tag_cent
    write(STDOUT,100) "number of unique non-atomic centers:   ", num_non_cents
    if (num_non_cents>0) then
      write(STDOUT,100) "position of unique non-atomic centers: ", id_non_cent(1:num_non_cents)
      write(STDOUT,100) "number of identical non-atomic centers:", num_non_iden(1:num_non_cents)
    end if
    write(STDOUT,100) "number of unique atomic centers:   ", num_atom_cents
    write(STDOUT,100) "position of unique atomic centers: ", id_atom_cent(1:num_atom_cents)
    write(STDOUT,100) "number of identical atomic centers:", num_atom_iden(1:num_atom_cents)
100 format("sort_gen_cents>> ",A,100I4)
#endif
    return
  end subroutine sort_gen_cents

  !> \brief sorts the indices of atomic centers using insertion sort algorithm, and
  !>        finds the unique centers and corresponding number of identical centers
  !> \author Bin Gao
  !> \date 2011-07-02
  !> \param num_cents is the number of centers
  !> \return idx_cent contains the indices of the centers, numbers should be greater than 0,
  !>         sorted in ascending order on exit
  !> \return tag_cent contains the tags marking the places of centers before sorting
  !> \return num_atom_cents is the number of unique atomic centers
  !> \return id_atom_cent contains the places of unique atomic centers
  !> \return num_atom_iden contains the number of identical atomic centers
  subroutine sort_atom_cents(num_cents, idx_cent, tag_cent, &
                             num_atom_cents, id_atom_cent, num_atom_iden)
    use xkind
    implicit none
    integer, intent(in) :: num_cents
    integer, intent(inout) :: idx_cent(num_cents)
    integer, intent(inout) :: tag_cent(num_cents)
    integer, intent(out) :: num_atom_cents
    integer, intent(out) :: id_atom_cent(num_cents)
    integer, intent(out) :: num_atom_iden(num_cents)
!f2py intent(hide) :: num_cents
!f2py intent(inout) :: idx_cent
!f2py intent(inout) :: tag_cent
!f2py depend(num_cents) :: tag_cent
!f2py intent(out) :: num_atom_cents
!f2py intent(out) :: id_atom_cent
!f2py depend(num_cents) :: id_atom_cent
!f2py intent(out) :: num_atom_iden
!f2py depend(num_cents) :: num_atom_iden
    integer ircd, jrcd     !incremental recorders
    integer idx_ins        !inserting index
    integer tag_ins        !tag of inserting index
    logical not_ins        !not found right place to insert the index
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
#if defined(DEBUG)
    write(STDOUT,100) "indices before sorting:", idx_cent
#endif
    ! initializes the number of unique atomic centers
    num_atom_cents = 1
    ! loops over the indices
    do ircd = 2, num_cents
      idx_ins = idx_cent(ircd)
      tag_ins = tag_cent(ircd)
      jrcd = ircd-1
      not_ins = .true.
      ! when an index is less than the previous one,
      ! loops down to insert it at its right place indicated by \var(jrcd)+1
      do while (not_ins)
        if (idx_ins<idx_cent(jrcd)) then
          idx_cent(jrcd+1) = idx_cent(jrcd)
          tag_cent(jrcd+1) = tag_cent(jrcd)
          jrcd = jrcd-1
          if (jrcd<1) then
            ! found another unique atomic center \var(idx_ins)
            num_atom_cents = num_atom_cents+1
            not_ins = .false.
          end if
        else
          ! found another unique atomic center \var(idx_ins)
          if (idx_ins>idx_cent(jrcd)) num_atom_cents = num_atom_cents+1
          not_ins = .false.
        end if
      end do
      idx_cent(jrcd+1) = idx_ins
      tag_cent(jrcd+1) = tag_ins
    end do
    ! tries to find the positions of unique centers, and the number of identical centers
    id_atom_cent(1) = 1
    num_atom_iden(1) = 1
    jrcd = 1
    do ircd = 2, num_cents
      ! found another unique atomic center
      if (idx_cent(ircd)/=idx_cent(id_atom_cent(jrcd))) then
        jrcd = jrcd+1
        id_atom_cent(jrcd) = ircd
        num_atom_iden(jrcd) = 1
      else
        num_atom_iden(jrcd) = num_atom_iden(jrcd)+1
      end if
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "sort_atom_cents", STDOUT)
#endif
#if defined(DEBUG)
    write(STDOUT,100) "indices after sorting: ", idx_cent
    write(STDOUT,100) "previous places:       ", tag_cent
    write(STDOUT,100) "number of unique atomic centers:   ", num_atom_cents
    write(STDOUT,100) "position of unique atomic centers: ", id_atom_cent(1:num_atom_cents)
    write(STDOUT,100) "number of identical atomic centers:", num_atom_iden(1:num_atom_cents)
100 format("sort_atom_cents>> ",A,100I4)
#endif
    return
  end subroutine sort_atom_cents
