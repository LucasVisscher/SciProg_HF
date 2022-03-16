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
!!  This file gives the next permutation for a given integer sequence.
!!
!!  2011-11-07, Bin Gao:
!!  * first version

  !> \brief gives the next permutation for a given integer sequence,
  !>        modified from function \fn(next_permutation) in C++ library,
  !>        see for example, http://marknelson.us/2002/03/01/next-permutation/
  !> \author Bin Gao
  !> \date 2011-11-07
  !> \param first is the index of first number
  !> \param last is the index of last number
  !> \return perm_num contains the next permutation on exit
  !> \return do_perm indicates if \fn(next_permutation) still could be called
  !>         for next permutation
  subroutine next_permutation(first, last, perm_num, do_perm)
    implicit none
    integer, intent(in) :: first
    integer, intent(in) :: last
    integer, intent(inout) :: perm_num(last)
    logical, intent(out) :: do_perm
!f2py intent(in) :: first
!f2py intent(in) :: last
!f2py intent(inout) :: perm_num
!f2py depend(first) :: perm_num
!f2py depend(last) :: perm_num
!f2py intent(out) :: do_perm
    integer idx, jdx, kdx, ldx  !incremental recorders
    integer tmp_num             !temporary number for swapping
    if (first==last+1) then
      do_perm = .false.
      return
    end if
    idx = first+1
    if (idx==last+1) then
      do_perm = .false.
      return
    end if
    idx = last
    do while (.true.)
      jdx = idx
      idx = idx-1
      if (perm_num(idx)<perm_num(jdx)) then
        kdx = last
        do while (perm_num(idx)>=perm_num(kdx))
          kdx = kdx-1
        end do
        ! swaps numbers at positions \var(idx) and \var(kdx)
        tmp_num = perm_num(idx)
        perm_num(idx) = perm_num(kdx)
        perm_num(kdx) = tmp_num
        ! reverses numbers from \var(jdx) to \var(last)
        ldx = last+jdx
        do kdx = jdx, (ldx-mod(ldx,2))/2
          tmp_num = perm_num(kdx)
          perm_num(kdx) = perm_num(ldx-kdx)
          perm_num(ldx-kdx) = tmp_num
        end do
        do_perm = .true.
        exit
      end if
      if (idx==first) then
        ! reverses numbers from \var(first) to \var(last)
        ldx = last+first
        do kdx = first, (ldx-mod(ldx,2))/2
          tmp_num = perm_num(kdx)
          perm_num(kdx) = perm_num(ldx-kdx)
          perm_num(ldx-kdx) = tmp_num
        end do
        do_perm = .false.
        exit
      end if
    end do
    return
  end subroutine next_permutation
