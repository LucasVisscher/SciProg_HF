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
!!  This file provides the Fortran 90 module of geometric derivatives.
!!
!!  2012-03-20, Bin Gao:
!!  * first version

#include "stdout.h"
#include "err_info.h"

!> \brief Fortran 90 module of geometric derivatives
!> \author Bin Gao
!> \date 2012-03-20
module gen1int_geom

  use xkind
  use london_ao
  implicit none

  ! information of N-ary tree for geometric derivatives
  type, public :: nary_tree_t
    private
    ! number of atoms
    integer :: num_atoms = 0
    ! indices of atoms to be used as differentiated centers
    integer, allocatable :: idx_atoms(:)
    ! order of geometric derivatives
    integer :: order_geo = 0
    ! maximum number of differentiated centers
    integer :: max_num_cent = 0
    ! total number of different allowed paths in N-ary tree
    integer :: num_paths = 0
    ! number of unique geometric derivatives in the N-ary tree
    integer :: num_unique_geo = 0
    ! index of current path
    integer :: idx_path = 0
    ! height of atom to visit
    integer visit_height
    ! indices of the selected atom nodes
    integer, allocatable :: idx_node(:)
    ! weights of the selected atom nodes
    integer, allocatable :: wt_node(:)
    ! indices of the generated differentiated centers
    integer, allocatable :: idx_cent(:)
    ! orders of derivatives of the differentiated centers
    integer, allocatable :: order_cent(:)
    ! offset of unique derivatives of current path
    integer :: path_offset = 0
    ! number of unique derivatives of current path
    integer :: path_num_unique = 0
    ! number of redundant derivatives of current path
    integer :: path_num_redunt = 0
  end type nary_tree_t

  public :: NaryTreeCreate
  public :: NaryTreeSetAtoms
  public :: NaryTreeGetNumAtoms
  public :: NaryTreeGetOrder
  public :: NaryTreeGetMaxNumCenters
  public :: NaryTreeGetNumPaths
  public :: NaryTreeGetNumGeo
  public :: NaryTreeSearch
  public :: NaryTreePathGetIndex
  public :: NaryTreePathGetNumCenters
  public :: NaryTreePathGetIdxCent
  public :: NaryTreePathGetOrderCent
  public :: NaryTreePathGetOffset
  public :: NaryTreePathGetNumUnique
  public :: NaryTreePathGetNumRedunt
  public :: NaryTreePathGetReduntList
  public :: NaryTreePathSetReduntExpt
  !-public :: NaryTreePathReorderInts
  public :: NaryTreeView
  public :: NaryTreeDestroy

  public :: IntGetCARMOM
  public :: IntGetNUCPOT
  public :: IntGetGAUPOT
  public :: IntGetODST

  contains

  !> \brief initializes the information of N-ary tree for geometric derivatives, and
  !>        returns the total number of different paths, and generates the first path
  !> \author Bin Gao
  !> \date 2011-12-12
  !> \param num_atoms is the number of atoms
  !> \param order_geo is the order of geometric derivatives
  !> \param max_num_cent is the maximum number of differentiated centers
  !> \return nary_tree contains the information of N-ary tree for geometric derivatives
  !> \return info_geom (==ERR_INFO) indicates the N-ary tree is not successfully created
  !> \return num_paths is the total number of different paths
  !> \return path_num_unique is the number of unique derivatives of the first path
  !> \return path_num_redunt is the number of redundant geometric derivatives of the first path
  subroutine NaryTreeCreate(num_atoms, order_geo, max_num_cent, nary_tree, info_geom, &
                            num_paths, path_num_unique, path_num_redunt)
    integer, intent(in) :: num_atoms
    integer, intent(in) :: order_geo
    integer, intent(in) :: max_num_cent
    type(nary_tree_t), intent(inout) :: nary_tree
    integer, intent(out) :: info_geom
    integer, optional, intent(out) :: num_paths
    integer, optional, intent(out) :: path_num_unique
    integer, optional, intent(out) :: path_num_redunt
    integer num_cent       !number of differentiated centers
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (num_atoms<1 .or. order_geo<0) then
      info_geom = ERR_INFO
    else
      info_geom = 0
      nary_tree%num_atoms = num_atoms
      nary_tree%order_geo = order_geo
      ! no geometric derivatives
      if (nary_tree%order_geo==0) then
        nary_tree%max_num_cent = 0
        allocate(nary_tree%idx_node(0:0), stat=info_geom)
        if (info_geom/=0) return  
        allocate(nary_tree%wt_node(0:0), stat=info_geom)
        if (info_geom/=0) then
          deallocate(nary_tree%idx_node)
          return                   
        end if
        allocate(nary_tree%idx_cent(1), stat=info_geom)
        if (info_geom/=0) then
          deallocate(nary_tree%idx_node)
          deallocate(nary_tree%wt_node)
          return
        end if
        allocate(nary_tree%order_cent(1), stat=info_geom)
        if (info_geom/=0) then
          deallocate(nary_tree%idx_node)
          deallocate(nary_tree%wt_node)
          deallocate(nary_tree%idx_cent)
          return
        end if
        nary_tree%num_paths = 1
        if (present(num_paths)) num_paths = nary_tree%num_paths
        nary_tree%idx_path = 1
        nary_tree%visit_height = 0
        nary_tree%idx_node = 0
        nary_tree%wt_node = 0
        nary_tree%idx_cent = 0
        nary_tree%order_cent = 0
        nary_tree%path_offset = 0
        nary_tree%path_num_unique = 1
        nary_tree%path_num_redunt = 1
        nary_tree%num_unique_geo = 1
        if (present(path_num_unique)) path_num_unique = nary_tree%path_num_unique
        if (present(path_num_redunt)) path_num_redunt = nary_tree%path_num_redunt
      ! having geometric derivatives
      else
        nary_tree%max_num_cent = max_num_cent
        ! if the number of atoms is less than the number of differentiated centers
        nary_tree%max_num_cent = min(nary_tree%num_atoms,nary_tree%max_num_cent)
        if (nary_tree%max_num_cent<1) then
          info_geom = ERR_INFO
          return
        end if
        allocate(nary_tree%idx_node(order_geo), stat=info_geom)
        if (info_geom/=0) return
        allocate(nary_tree%wt_node(order_geo), stat=info_geom)
        if (info_geom/=0) then
          deallocate(nary_tree%idx_node)
          return
        end if
        allocate(nary_tree%idx_cent(nary_tree%max_num_cent), stat=info_geom)
        if (info_geom/=0) then
          deallocate(nary_tree%idx_node)
          deallocate(nary_tree%wt_node)
          return
        end if
        allocate(nary_tree%order_cent(nary_tree%max_num_cent), stat=info_geom)
        if (info_geom/=0) then
          deallocate(nary_tree%idx_node)
          deallocate(nary_tree%wt_node)
          deallocate(nary_tree%idx_cent)
          return
        end if
        ! returns the total number of different paths, and generates the first path
        call geom_total_tree_init(num_atoms, order_geo,   &
                                  nary_tree%max_num_cent, &
                                  nary_tree%num_paths,    &
                                  nary_tree%visit_height, &
                                  nary_tree%idx_node,     &
                                  nary_tree%wt_node,      &
                                  nary_tree%idx_cent,     &
                                  nary_tree%order_cent,   &
                                  nary_tree%path_num_unique)
        ! gets the number of redundant derivatives of current path
        num_cent = nary_tree%wt_node(order_geo)
        call geom_total_num_redunt(num_cent, nary_tree%order_cent(1:num_cent), &
                                   nary_tree%path_num_redunt)
        if (present(path_num_unique)) path_num_unique = nary_tree%path_num_unique
        if (present(path_num_redunt)) path_num_redunt = nary_tree%path_num_redunt
        if (present(num_paths)) num_paths = nary_tree%num_paths
        nary_tree%idx_path = 1
        ! gets the number of unique geometric derivatives in the N-ary tree
        call geom_total_num_derv(order_geo, nary_tree%max_num_cent, num_atoms, &
                                 nary_tree%num_unique_geo)
      end if
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "NaryTreeCreate", STDOUT)
#endif
  end subroutine NaryTreeCreate

  !> \brief sets the indices of atoms to be used as differentiated centers
  !> \author Bin Gao
  !> \date 2012-05-09
  !> \param num_atoms is the number of atoms
  !> \param idx_atoms contains the indices of atoms to be used as differentiated centers
  !> \return nary_tree contains the information of N-ary tree for geometric derivatives
  !> \return info_geom (==ERR_INFO) indicates the N-ary tree is not successfully created
  subroutine NaryTreeSetAtoms(num_atoms, idx_atoms, nary_tree, info_geom)
    integer, intent(in) :: num_atoms
    integer, intent(in) :: idx_atoms(num_atoms)
    type(nary_tree_t), intent(inout) :: nary_tree
    integer, intent(out) :: info_geom
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (num_atoms==nary_tree%num_atoms) then
      allocate(nary_tree%idx_atoms(num_atoms), stat=info_geom)
      if (info_geom==0) then
        nary_tree%idx_atoms = idx_atoms
      end if
    else
      info_geom = num_atoms
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "NaryTreeSetAtoms", STDOUT)
#endif
  end subroutine NaryTreeSetAtoms

  !> \brief gets the number of atoms for a given N-ary tree
  !> \author Bin Gao
  !> \date 2012-05-09
  !> \param nary_tree contains the information of N-ary tree for geometric derivatives
  !> \return num_atoms is the number of atoms
  subroutine NaryTreeGetNumAtoms(nary_tree, num_atoms)
    type(nary_tree_t), intent(in) :: nary_tree
    integer, intent(out) :: num_atoms
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    num_atoms = nary_tree%num_atoms
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "NaryTreeGetNumAtoms", STDOUT)
#endif
  end subroutine NaryTreeGetNumAtoms

  !> \brief gets the order of geometric derivatives for a given N-ary tree
  !> \author Bin Gao
  !> \date 2012-05-09
  !> \param nary_tree contains the information of N-ary tree for geometric derivatives
  !> \return order_geo is the order of geometric derivatives
  subroutine NaryTreeGetOrder(nary_tree, order_geo)
    type(nary_tree_t), intent(in) :: nary_tree
    integer, intent(out) :: order_geo
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    order_geo = nary_tree%order_geo
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "NaryTreeGetOrder", STDOUT)
#endif
  end subroutine NaryTreeGetOrder

  !> \brief gets the maximum number of differentiated centers for a given N-ary tree
  !> \author Bin Gao
  !> \date 2012-05-09
  !> \param nary_tree contains the information of N-ary tree for geometric derivatives
  !> \return max_num_cent is the maximum number of differentiated centers
  subroutine NaryTreeGetMaxNumCenters(nary_tree, max_num_cent)
    type(nary_tree_t), intent(in) :: nary_tree
    integer, intent(out) :: max_num_cent
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    max_num_cent = nary_tree%max_num_cent
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "NaryTreeGetMaxNumCenters", STDOUT)
#endif
  end subroutine NaryTreeGetMaxNumCenters

  !> \brief gets the total number of different allowed paths in N-ary tree
  !> \author Bin Gao
  !> \date 2012-03-09
  !> \param nary_tree contains the information of N-ary tree for geometric derivatives
  !> \return num_paths is the total number of different allowed paths in N-ary tree
  subroutine NaryTreeGetNumPaths(nary_tree, num_paths)
    type(nary_tree_t), intent(in) :: nary_tree
    integer, intent(out) :: num_paths
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! checks if the N-ary tree is created
    if (nary_tree%idx_path==0) then
      call error_stop("NaryTreeGetNumPaths", "N-ary tree is not created", 0)
    else
      num_paths = nary_tree%num_paths
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "NaryTreeGetNumPaths", STDOUT)
#endif
  end subroutine NaryTreeGetNumPaths

  !> \brief gets the number of unique geometric derivatives in the N-ary tree
  !> \author Bin Gao
  !> \date 2012-01-12
  !> \param nary_tree contains the information of N-ary tree for geometric derivatives
  !> \return num_unique_geo is the number of unique geometric derivatives in the N-ary tree
  subroutine NaryTreeGetNumGeo(nary_tree, num_unique_geo)
    type(nary_tree_t), intent(in) :: nary_tree
    integer, intent(out) :: num_unique_geo
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! checks if the N-ary tree is created
    if (nary_tree%idx_path==0) then
      call error_stop("NaryTreeGetNumGeo", "N-ary tree is not created", 0)
    else
      num_unique_geo = nary_tree%num_unique_geo
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "NaryTreeGetNumGeo", STDOUT)
#endif
  end subroutine NaryTreeGetNumGeo

  !> \brief searches for the next satisfied path from a given path, could be called recursively
  !> \author Bin Gao
  !> \date 2011-12-12
  !> \param nary_tree contains the information of N-ary tree for geometric derivatives
  !> \return path_num_unique is the number of unique derivatives of current path
  !> \return path_num_redunt is the number of redundant derivatives of current path
  subroutine NaryTreeSearch(nary_tree, path_num_unique, path_num_redunt)
    type(nary_tree_t), intent(inout) :: nary_tree
    integer, optional, intent(out) :: path_num_unique
    integer, optional, intent(out) :: path_num_redunt
    integer num_cent       !number of differentiated centers
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! checks if the N-ary tree is created
    if (nary_tree%idx_path==0) then
      call error_stop("NaryTreeSearch", "N-ary tree is not created", 0)
    else if (nary_tree%order_geo==0) then
      if (present(path_num_unique)) path_num_unique = nary_tree%path_num_unique
      if (present(path_num_redunt)) path_num_redunt = nary_tree%path_num_redunt
    else
      ! updates the offset of unique derivatives of current path
      if (nary_tree%idx_path==nary_tree%num_paths) then
        nary_tree%path_offset = 0
        nary_tree%idx_path = 0
      else
        nary_tree%path_offset = nary_tree%path_offset+nary_tree%path_num_unique
      end if
      ! searches for the next satisfied path from the given path
      call geom_total_tree_search(nary_tree%num_atoms, nary_tree%order_geo,       &
                                  nary_tree%max_num_cent, nary_tree%visit_height, &
                                  nary_tree%idx_node, nary_tree%wt_node,          &
                                  nary_tree%idx_cent, nary_tree%order_cent,       &
                                  nary_tree%path_num_unique)
      ! gets the number of redundant derivatives of current path
      num_cent = nary_tree%wt_node(nary_tree%order_geo)
      call geom_total_num_redunt(num_cent, nary_tree%order_cent(1:num_cent), &
                                 nary_tree%path_num_redunt)
      if (present(path_num_unique)) path_num_unique = nary_tree%path_num_unique
      if (present(path_num_redunt)) path_num_redunt = nary_tree%path_num_redunt
      nary_tree%idx_path = nary_tree%idx_path+1
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "NaryTreeSearch", STDOUT)
#endif
  end subroutine NaryTreeSearch

  !> \brief gets the index of current path
  !> \author Bin Gao
  !> \date 2012-03-09
  !> \param nary_tree contains the information of N-ary tree for geometric derivatives
  !> \return idx_path is the index of current path
  subroutine NaryTreePathGetIndex(nary_tree, idx_path)
    type(nary_tree_t), intent(in) :: nary_tree
    integer, intent(out) :: idx_path
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! checks if the N-ary tree is created
    if (nary_tree%idx_path==0) then
      call error_stop("NaryTreePathGetIndex", "N-ary tree is not created", 0)
    else
      idx_path = nary_tree%idx_path
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "NaryTreePathGetIndex", STDOUT)
#endif
  end subroutine NaryTreePathGetIndex

  !> \brief returns the number of differentiated centers of current path for given N-ary tree
  !> \author Bin Gao
  !> \date 2012-03-20
  !> \param nary_tree contains the information of N-ary tree for geometric derivatives
  !> \return num_centers is the number of differentiated centers of current path
  subroutine NaryTreePathGetNumCenters(nary_tree, num_centers)
    type(nary_tree_t), intent(in) :: nary_tree
    integer, intent(out) :: num_centers
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! checks if the N-ary tree is created
    if (nary_tree%idx_path==0) then
      call error_stop("NaryTreePathGetNumCenters", "N-ary tree is not created", 0)
    else
      num_centers = nary_tree%wt_node(nary_tree%order_geo)
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "NaryTreePathGetNumCenters", STDOUT)
#endif
  end subroutine NaryTreePathGetNumCenters

  !> \brief returns the indices of differentiated centers of current path for given N-ary tree
  !> \author Bin Gao
  !> \date 2013-05-07
  !> \param nary_tree contains the information of N-ary tree for geometric derivatives
  !> \return idx_cent contains the indices of differentiated centers of current path
  subroutine NaryTreePathGetIdxCent(nary_tree, idx_cent)
    type(nary_tree_t), intent(in) :: nary_tree
    integer, intent(out) :: idx_cent(:)
    integer num_centers    !number of differentiated centers of current path
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! checks if the N-ary tree is created
    if (nary_tree%idx_path==0) then
      call error_stop("NaryTreePathGetIdxCent", "N-ary tree is not created", 0)
    else
      num_centers = nary_tree%wt_node(nary_tree%order_geo)
      if (size(idx_cent)<num_centers) then
        call error_stop("NaryTreePathGetIdxCent",      &
                        "size of idx_cent not enough", &
                        num_centers)
      else
        idx_cent(1:num_centers) = nary_tree%idx_cent(1:num_centers)
      end if
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "NaryTreePathGetIdxCent", STDOUT)
#endif
  end subroutine NaryTreePathGetIdxCent

  !> \brief returns the orders of differentiated centers of current path for given N-ary tree
  !> \author Bin Gao
  !> \date 2013-05-07
  !> \param nary_tree contains the information of N-ary tree for geometric derivatives
  !> \return order_cent contains the orders of differentiated centers of current path
  subroutine NaryTreePathGetOrderCent(nary_tree, order_cent)
    type(nary_tree_t), intent(in) :: nary_tree
    integer, intent(out) :: order_cent(:)
    integer num_centers    !number of differentiated centers of current path
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! checks if the N-ary tree is created
    if (nary_tree%idx_path==0) then
      call error_stop("NaryTreePathGetOrderCent", "N-ary tree is not created", 0)
    else
      num_centers = nary_tree%wt_node(nary_tree%order_geo)
      if (size(order_cent)<num_centers) then
        call error_stop("NaryTreePathGetOrderCent",      &
                        "size of order_cent not enough", &
                        num_centers)
      else
        order_cent(1:num_centers) = nary_tree%order_cent(1:num_centers)
      end if
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "NaryTreePathGetOrderCent", STDOUT)
#endif
  end subroutine NaryTreePathGetOrderCent

  !> \brief returns the offset of unique derivatives of current path
  !> \author Bin Gao
  !> \date 2012-01-12
  !> \param nary_tree contains the information of N-ary tree for geometric derivatives
  !> \return path_offset is the offset of unique derivatives of current path
  subroutine NaryTreePathGetOffset(nary_tree, path_offset)
    type(nary_tree_t), intent(in) :: nary_tree
    integer, intent(out) :: path_offset
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! checks if the N-ary tree is created
    if (nary_tree%idx_path==0) then
      call error_stop("NaryTreePathGetOffset", "N-ary tree is not created", 0)
    else
      path_offset = nary_tree%path_offset
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "NaryTreePathGetOffset", STDOUT)
#endif
  end subroutine NaryTreePathGetOffset

  !> \brief returns the number of unique derivatives of current path
  !> \author Bin Gao
  !> \date 2012-01-12
  !> \param nary_tree contains the information of N-ary tree for geometric derivatives
  !> \return path_num_unique is the number of unique derivatives of current path
  subroutine NaryTreePathGetNumUnique(nary_tree, path_num_unique)
    type(nary_tree_t), intent(in) :: nary_tree
    integer, intent(out) :: path_num_unique
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! checks if the N-ary tree is created
    if (nary_tree%idx_path==0) then
      call error_stop("NaryTreePathGetNumUnique", "N-ary tree is not created", 0)
    else
      path_num_unique = nary_tree%path_num_unique
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "NaryTreePathGetNumUnique", STDOUT)
#endif
  end subroutine NaryTreePathGetNumUnique

  !> \brief returns the number of redundant derivatives of current path
  !> \author Bin Gao
  !> \date 2012-01-12
  !> \param nary_tree contains the information of N-ary tree for geometric derivatives
  !> \return path_num_redunt is the number of redundant derivatives of current path
  subroutine NaryTreePathGetNumRedunt(nary_tree, path_num_redunt)
    type(nary_tree_t), intent(in) :: nary_tree
    integer, intent(out) :: path_num_redunt
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! checks if the N-ary tree is created
    if (nary_tree%idx_path==0) then
      call error_stop("NaryTreePathGetNumRedunt", "N-ary tree is not created", 0)
    else
      path_num_redunt = nary_tree%path_num_redunt
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "NaryTreePathGetNumRedunt", STDOUT)
#endif
  end subroutine NaryTreePathGetNumRedunt

  !> \brief gets the list addresses of redundant derivatives of current path
  !> \author Bin Gao
  !> \date 2011-12-12
  !> \param nary_tree contains the information of N-ary tree for geometric derivatives
  !> \return redunt_list contains the list addresses of redundant derivatives of current path
  subroutine NaryTreePathGetReduntList(nary_tree, redunt_list)
    type(nary_tree_t), intent(in) :: nary_tree
    integer, intent(out) :: redunt_list(:,:)
    integer path_num_redunt  !number of redundant derivatives of current path
    integer num_cent         !number of differentiated centers
#if defined(XTIME)
    real(REALK) curr_time   !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! checks if the N-ary tree is created
    if (nary_tree%idx_path==0) then
      call error_stop("NaryTreePathGetReduntList", "N-ary tree is not created", 0)
    else if (nary_tree%order_geo==0) then
      redunt_list = 1
    else
      ! checks the sizes of the list addresses
      if (size(redunt_list,1)/=2)                                 &
        call error_stop("NaryTreePathGetReduntList",                  &
                        "incorrect size of variable redunt_list", &
                        size(redunt_list,1))
      path_num_redunt = size(redunt_list,2)
      if (path_num_redunt/=nary_tree%path_num_redunt)                &
        call error_stop("NaryTreePathGetReduntList",                     &
                        "incorrect number of redundant derivatives", &
                        path_num_redunt)
      ! gets the list addresses of redundant derivatives of current path
      num_cent = nary_tree%wt_node(nary_tree%order_geo)
      call geom_total_redunt_list(nary_tree%num_atoms, num_cent,    &
                                  nary_tree%idx_cent(1:num_cent),   &
                                  nary_tree%order_cent(1:num_cent), &
                                  .true., path_num_redunt, redunt_list)
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "NaryTreePathGetReduntList", STDOUT)
#endif
  end subroutine NaryTreePathGetReduntList

  !> \brief returns the expectation values with redundant derivatives of current path
  !> \author Bin Gao
  !> \date 2011-01-14
  !> \param nary_tree contains the information of N-ary tree for geometric derivatives
  !> \param num_opt is the number of operators
  !> \param path_num_unique is the number of unique derivatives of current path
  !> \param num_dens is the number of AO density matrices
  !> \param unique_expt contains the expectation values of unique geometric derivatives
  !> \param num_redunt_geo is the total number of redundant derivatives, equals to
  !>        \var(3*num_atoms)^sum(\var(order_cent))
  !> \return redunt_expt contains the updated expectation values of redundant derivatives on exit
  subroutine NaryTreePathSetReduntExpt(nary_tree, num_opt, path_num_unique, num_dens, &
                                   unique_expt, num_redunt_geo, redunt_expt)
    type(nary_tree_t), intent(in) :: nary_tree
    integer, intent(in) :: num_opt
    integer, intent(in) :: path_num_unique
    integer, intent(in) :: num_dens
    real(REALK), intent(in) :: unique_expt(num_opt,path_num_unique,num_dens)
    integer, intent(in) :: num_redunt_geo
    real(REALK), intent(inout) :: redunt_expt(num_opt,num_redunt_geo,num_dens)
    integer num_cent       !number of differentiated centers
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! checks if the N-ary tree is created
    if (nary_tree%idx_path==0) then
      call error_stop("NaryTreePathSetReduntExpt", "N-ary tree is not created", 0)
    else if (nary_tree%order_geo==0) then
      if (num_redunt_geo/=path_num_unique)                      &
        call error_stop("NaryTreePathSetReduntExpt",                &
                        "incorrect size of expectation values", &
                        path_num_unique)
      redunt_expt = redunt_expt+unique_expt
    else
      ! checks the sizes of the expectation values
      if (path_num_unique/=nary_tree%path_num_unique)         &
        call error_stop("NaryTreePathSetReduntExpt",              &
                        "incorrect variable path_num_unique", &
                        path_num_unique)
      if (num_redunt_geo/=(3*nary_tree%num_atoms)**nary_tree%order_geo) &
        call error_stop("NaryTreePathSetReduntExpt",                        &
                        "incorrect variable num_redunt_geo",            &
                        num_redunt_geo)
      ! gets the expectation values with redundant geometric derivatives for current path
      num_cent = nary_tree%wt_node(nary_tree%order_geo)
      call geom_total_redunt_expectation(nary_tree%num_atoms, num_cent,      &
                                         nary_tree%idx_cent(1:num_cent),     &
                                         nary_tree%order_cent(1:num_cent),   &
                                         num_opt, path_num_unique, num_dens, &
                                         unique_expt, num_redunt_geo, redunt_expt)
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "NaryTreePathSetReduntExpt", STDOUT)
#endif
  end subroutine NaryTreePathSetReduntExpt

!-  !> \brief reorders geometric derivatives
!-  !> \author Bin Gao
!-  !> \date 2012-08-30
!-  !> \param nary_tree contains the information of N-ary tree for geometric derivatives
!-  !> \param dim_ints is the dimension of integrals
!-  !> \param path_num_unique is the number of unique derivatives of current path
!-  !> \param dim_opt is the dimension of operators
!-  !> \return contr_ints contains the contracted integrals
!-  subroutine NaryTreePathReorderInts(nary_tree, dim_ints, path_num_unique, dim_opt, contr_ints)
!-    type(nary_tree_t), intent(in) :: nary_tree
!-    integer, intent(in) :: dim_ints
!-    integer, intent(in) :: path_num_unique
!-    integer, intent(in) :: dim_opt
!-    real(REALK), intent(inout) :: contr_ints(dim_ints,path_num_unique,dim_opt)
!-    integer icent     !incremental recorder over differentiated centers
!-    integer num_ints  !number of integrals
!-    integer num_geo   !number of geometric derivatives
!-    integer num_opt   !number of operators
!-#if defined(XTIME)
!-    real(REALK) curr_time  !current CPU time
!-    ! sets current CPU time
!-    call xtimer_set(curr_time)
!-#endif
!-    ! checks if the N-ary tree is created
!-    if (nary_tree%idx_path==0) then
!-      call error_stop("NaryTreePathReorderInts", "N-ary tree is not created", 0)
!-    else if (nary_tree%order_geo>1) then
!-      ! checks the sizes of the expectation values
!-      if (path_num_unique/=nary_tree%path_num_unique)         &
!-        call error_stop("NaryTreePathReorderInts",                &
!-                        "incorrect variable path_num_unique", &
!-                        path_num_unique)
!-      num_ints = dim_ints
!-      num_opt = path_num_unique*dim_opt
!-      do icent = 1, nary_tree%wt_node(nary_tree%order_geo)
!-        num_geo = (nary_tree%order_cent(icent)+1)*(nary_tree%order_cent(icent)+2)/2
!-        num_opt = num_opt/num_geo
!-        call geom_part_reorder(nary_tree%order_cent(icent), &
!-                               num_ints, num_opt, contr_ints) 
!-        num_ints = num_ints*num_geo
!-      end do
!-    end if
!-#if defined(XTIME)
!-    ! prints the CPU elapsed time
!-    call xtimer_view(curr_time, "NaryTreePathReorderInts", STDOUT)
!-#endif
!-  end subroutine NaryTreePathReorderInts

  !> \brief visualizes the information of N-ary tree for geometric derivatives
  !> \author Bin Gao
  !> \date 2011-12-12
  !> \param nary_tree contains the information of N-ary tree for geometric derivatives
  !> \param io_viewer is the logical unit number of the viewer
  subroutine NaryTreeView(nary_tree, io_viewer)
    type(nary_tree_t), intent(in) :: nary_tree
    integer, intent(in) :: io_viewer
    integer icent          !incremental recorder over centers
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! checks if the N-ary tree is created
    if (nary_tree%idx_path==0) then
      write(io_viewer,100) "N-ary tree is not created!"
      return
    else
      write(io_viewer,100) "number of atoms", nary_tree%num_atoms
      if (allocated(nary_tree%idx_atoms)) then
        do icent = 1, nary_tree%num_atoms
          write(io_viewer,100) "indices of atoms to be used as differentiated centers", &
                               nary_tree%idx_atoms(icent)
        end do
      end if
      write(io_viewer,100) "order of geometric derivatives", nary_tree%order_geo
      write(io_viewer,100) "maximum number of differentiated centers", nary_tree%max_num_cent
      write(io_viewer,100) "total number of different allowed paths in N-ary tree", &
                           nary_tree%num_paths
      write(io_viewer,100) "number of unique geometric derivatives in the N-ary tree", &
                           nary_tree%num_unique_geo
      write(io_viewer,100) "index of current path", nary_tree%idx_path
      if (nary_tree%order_geo>0) then
        write(io_viewer,110) (nary_tree%idx_node(icent),"(",nary_tree%wt_node(icent), &
                              ")",icent=1,nary_tree%order_geo)
        write(io_viewer,120) (nary_tree%idx_cent(icent),"(",nary_tree%order_cent(icent), &
                              ")",icent=1,nary_tree%wt_node(nary_tree%order_geo))
      end if
      write(io_viewer,100) "offset of unique derivatives of current path", &
        nary_tree%path_offset
      write(io_viewer,100) "number of unique derivatives of current path", &
        nary_tree%path_num_unique
      write(io_viewer,100) "number of redundant derivatives of current path", &
        nary_tree%path_num_redunt
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "NaryTreeView", STDOUT)
#endif
100 format("NaryTreeView>> ",A,I8)
110 format("NaryTreeView>> indices and weights of atom nodes",12(I4,A,I3,A))
120 format("NaryTreeView>> generated differentiated centers",16(I3,A,I2,A))
  end subroutine NaryTreeView

  !> \brief frees space taken by a N-ary tree for geometric derivatives
  !> \author Bin Gao
  !> \date 2011-12-12
  !> \param nary_tree contains the information of N-ary tree for geometric derivatives
  subroutine NaryTreeDestroy(nary_tree)
    type(nary_tree_t), intent(inout) :: nary_tree
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (allocated(nary_tree%idx_atoms)) deallocate(nary_tree%idx_atoms)
    if (allocated(nary_tree%idx_node)) deallocate(nary_tree%idx_node)
    if (allocated(nary_tree%wt_node)) deallocate(nary_tree%wt_node)
    if (allocated(nary_tree%idx_cent)) deallocate(nary_tree%idx_cent)
    if (allocated(nary_tree%order_cent)) deallocate(nary_tree%order_cent)
    nary_tree%num_atoms = 0
    nary_tree%order_geo = 0
    nary_tree%max_num_cent = 0
    nary_tree%num_paths = 0
    nary_tree%idx_path = 0
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "NaryTreeDestroy", STDOUT)
#endif
  end subroutine NaryTreeDestroy

  !> \brief calculates the Cartesian multipole moment integrals using contracted
  !>        Gaussian type orbitals (GTOs)
  !> \author Bin Gao
  !> \date 2011-12-12
  !> \param idx_bra is the atomic index of bra center
  !> \param coord_bra contains the coordinates of bra center
  !> \param angular_bra is the angular number of bra center (s=0, p=1, d=2, ...)
  !> \param num_prim_bra is the number of primitive Gaussians of bra center
  !> \param exponent_bra contains the exponents of primitive Gaussians of bra center
  !> \param num_contr_bra is the number of contractions of bra center
  !> \param contr_coef_bra contains the contraction coefficients of bra center
  !> \param idx_ket is the atomic index of ket center
  !> \param coord_ket contains the coordinates of ket center
  !> \param angular_ket is the angular number of ket center
  !> \param num_prim_ket is the number of primitive Gaussians of ket center
  !> \param exponent_ket contains the exponents of primitive Gaussians of ket center
  !> \param num_contr_ket is the number of contractions of ket center
  !> \param contr_coef_ket contains the contraction coefficients of ket center
  !> \param spher_gto indicates if using spherical GTOs, otherwise Cartesian GTOs
  !> \param info_LAO contains the information of London atomic orbital
  !> \param order_elec is the order of electronic derivatives
  !> \param idx_diporg is the atomic center of dipole origin (<1 for non-atomic center)
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scale constant for Cartesian multipole moments
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_mag_bra is the order of magnetic derivatives on bra center
  !> \param order_mag_ket is the order of magnetic derivatives on ket center
  !> \param order_mag_total is the order of total magnetic derivatives
  !> \param order_ram_bra is the order of derivatives w.r.t. total rotational angular momentum on bra center
  !> \param order_ram_ket is the order of derivatives w.r.t. total rotational angular momentum on ket center
  !> \param order_ram_total is the order of total derivatives w.r.t. total rotational angular momentum
  !> \param order_geo_bra is the order of geometric derivatives with respect to bra center
  !> \param order_geo_ket is the order of geometric derivatives with respect to ket center
  !> \param order_geo_mom is the order of geometric derivatives on dipole origin
  !> \param nary_tree_total contains the information of N-ary tree for total geometric derivatives
  !> \param num_gto_bra is the number of spherical/Cartesian GTOs on bra center
  !> \param num_gto_ket is the number of spherical/Cartesian GTOs on ket center
  !> \param num_opt is the number of operators including derivatives
  !> \param mag_num_bra contains the magnetic numbers of spherical GTOs on bra center
  !> \param mag_num_ket contains the magnetic numbers of spherical GTOs on ket center
  !> \param powers_bra contains the Cartesian powers of Cartesian GTOs on bra center
  !> \param powers_ket contains the Cartesian powers of Cartesian GTOs on ket center
  !> \return contr_ints contains the calculated contracted integrals
  subroutine IntGetCARMOM(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                          exponent_bra, num_contr_bra, contr_coef_bra,   &
                          idx_ket, coord_ket, angular_ket, num_prim_ket, &
                          exponent_ket, num_contr_ket, contr_coef_ket,   &
                          spher_gto, info_LAO, order_elec, idx_diporg,   &
                          dipole_origin, scal_const, order_mom,          &
                          order_mag_bra, order_mag_ket, order_mag_total, &
                          order_ram_bra, order_ram_ket, order_ram_total, &
                          order_geo_bra, order_geo_ket, order_geo_mom,   &
                          nary_tree_total, num_gto_bra, num_gto_ket,     &
                          num_opt, contr_ints, mag_num_bra, mag_num_ket, &
                          powers_bra, powers_ket)
    integer, intent(in) :: idx_bra
    real(REALK), intent(in) :: coord_bra(3)
    integer, intent(in) :: angular_bra
    integer, intent(in) :: num_prim_bra
    real(REALK), intent(in) :: exponent_bra(num_prim_bra)
    integer, intent(in) :: num_contr_bra
    real(REALK), intent(in) :: contr_coef_bra(num_contr_bra,num_prim_bra)
    integer, intent(in) :: idx_ket
    real(REALK), intent(in) :: coord_ket(3)
    integer, intent(in) :: angular_ket
    integer, intent(in) :: num_prim_ket
    real(REALK), intent(in) :: exponent_ket(num_prim_ket)
    integer, intent(in) :: num_contr_ket
    real(REALK), intent(in) :: contr_coef_ket(num_contr_ket,num_prim_ket)
    logical, optional, intent(in) :: spher_gto
    type(london_ao_t), intent(in) :: info_LAO
    integer, intent(in) :: order_elec
    integer, intent(in) :: idx_diporg
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_mag_bra
    integer, intent(in) :: order_mag_ket
    integer, intent(in) :: order_mag_total
    integer, intent(in) :: order_ram_bra
    integer, intent(in) :: order_ram_ket
    integer, intent(in) :: order_ram_total
    integer, optional, intent(in) :: order_geo_bra
    integer, optional, intent(in) :: order_geo_ket
    integer, optional, intent(in) :: order_geo_mom
    type(nary_tree_t), optional, intent(in) :: nary_tree_total
    integer, intent(in) :: num_gto_bra
    integer, intent(in) :: num_gto_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(out) :: contr_ints(num_gto_bra,num_contr_bra, &
                                           num_gto_ket,num_contr_ket,num_opt)
    integer, optional, intent(in) :: mag_num_bra(num_gto_bra)
    integer, optional, intent(in) :: mag_num_ket(num_gto_ket)
    integer, optional, intent(in) :: powers_bra(3,num_gto_bra)
    integer, optional, intent(in) :: powers_ket(3,num_gto_ket)
#include "max_idx_non.h"
    real(REALK) gauge_origin(3)                      !gauge origin of the magnetic vector potential
    real(REALK) origin_London_PF(3)                  !origin of the London phase factor
    integer iopt                                     !incremental recorder over operators
    logical p_spher_gto                              !arguments for Gen1Int (local)
    integer p_order_geo_bra
    integer p_order_geo_ket
    integer p_order_geo_mom
    integer p_num_cents
    integer, allocatable :: p_idx_cent(:)
    integer, allocatable :: p_order_cent(:)
    integer icent
    real(REALK), allocatable :: tmp_ints(:,:,:,:,:)  !contracted integrals from Gen1Int
    integer ierr                                     !error information
    ! sets the arguments for Gen1Int (local)
    if (present(spher_gto)) then
      p_spher_gto = spher_gto
    else
      p_spher_gto = .true.
    end if
    if (present(order_geo_bra)) then
      p_order_geo_bra = order_geo_bra
    else
      p_order_geo_bra = 0
    end if
    if (present(order_geo_ket)) then
      p_order_geo_ket = order_geo_ket
    else
      p_order_geo_ket = 0
    end if
    if (present(order_geo_mom)) then
      p_order_geo_mom = order_geo_mom
    else
      p_order_geo_mom = 0
    end if
    if (present(nary_tree_total)) then
      p_num_cents = nary_tree_total%wt_node(nary_tree_total%order_geo)
    else
      p_num_cents = 0
    end if
    if (p_num_cents>0) then
      allocate(p_idx_cent(p_num_cents), stat=ierr)
      if (ierr/=0) &
        call error_stop("IntGetCARMOM", "failed to allocate p_idx_cent", p_num_cents)
      allocate(p_order_cent(p_num_cents), stat=ierr)
      if (ierr/=0) &
        call error_stop("IntGetCARMOM", "failed to allocate p_order_cent", p_num_cents)
      if (allocated(nary_tree_total%idx_atoms)) then
        do icent = 1, p_num_cents
          p_idx_cent(icent) = nary_tree_total%idx_atoms(nary_tree_total%idx_cent(icent))
        end do
      else
        p_idx_cent = nary_tree_total%idx_cent(1:p_num_cents)
      end if
      p_order_cent = nary_tree_total%order_cent(1:p_num_cents)
    else
      allocate(p_idx_cent(1), stat=ierr)
      if (ierr/=0) call error_stop("IntGetCARMOM", "failed to allocate p_idx_cent", 1)
      allocate(p_order_cent(1), stat=ierr)
      if (ierr/=0) call error_stop("IntGetCARMOM", "failed to allocate p_order_cent", 1)
      p_idx_cent = 0
      p_order_cent = 0
    end if
    ! magnetic derivatives or derivatives with respect to rotational angular momentum
    if (order_mag_total>0 .or. order_mag_bra>0 .or. order_mag_ket>0 .or. &
        order_ram_total>0 .or. order_ram_bra>0 .or. order_ram_ket>0) then
      ! London atomic orbitals
      if (LondonAOUsed(info_LAO)) then
        if ((order_mag_total==0 .and. order_mom==0 .and.       &
             ((order_mag_bra==1 .and. order_mag_ket==0) .or.   &
              (order_mag_bra==0 .and. order_mag_ket==1)) .and. &
             order_elec==0 .and. p_order_geo_bra==0 .and.      &
             p_order_geo_ket==0) .or.                          &
            (order_mag_total==1 .and. order_mom==0 .and.       &
             order_mag_bra==0 .and. order_mag_ket==0 .and.     &
             order_elec==0 .and. p_order_geo_bra==0 .and.      &
             p_order_geo_ket==0)) then
          call LondonAOGetLPFOrigin(info_LAO, origin_London_PF)
          allocate(tmp_ints(num_gto_bra,num_contr_bra, &
                            num_gto_ket,num_contr_ket,num_opt), stat=ierr)
          if (ierr/=0)                                                     &
            call error_stop("IntGetCARMOM", "failed to allocate tmp_ints", &
                            num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket*num_opt)
          ! SGTOs
          if (p_spher_gto) then
            if ((angular_bra==0 .and. angular_ket==0) .or. &
                .not.(present(mag_num_bra) .or. present(mag_num_ket))) then
              call contr_sgto_carmom(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                     exponent_bra, num_contr_bra, contr_coef_bra,   &
                                     idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                     exponent_ket, num_contr_ket, contr_coef_ket,   &
                                     order_elec, MAX_IDX_NON, origin_London_PF,     &
                                     scal_const*0.5_REALK, 1, p_order_geo_bra,      &
                                     p_order_geo_ket, p_order_geo_mom,              &
                                     p_num_cents, p_idx_cent, p_order_cent,         &
                                     num_gto_bra, num_gto_ket, num_opt, tmp_ints)
            ! reorders integrals if required
            else
              call contr_sgto_carmom(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                     exponent_bra, num_contr_bra, contr_coef_bra,   &
                                     idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                     exponent_ket, num_contr_ket, contr_coef_ket,   &
                                     order_elec, MAX_IDX_NON, origin_London_PF,     &
                                     scal_const*0.5_REALK, 1, p_order_geo_bra,      &
                                     p_order_geo_ket, p_order_geo_mom,              &
                                     p_num_cents, p_idx_cent, p_order_cent,         &
                                     num_gto_bra, num_gto_ket, num_opt, contr_ints)
              if (angular_bra>0 .and. present(mag_num_bra) .and. &
                  angular_ket>0 .and. present(mag_num_ket)) then
                call reorder_sgto_ints(angular_bra, num_gto_bra, mag_num_bra, &
                                       angular_ket, num_gto_ket, mag_num_ket, &
                                       num_contr_bra, num_contr_ket, num_opt, &
                                       contr_ints, tmp_ints)
              else if (angular_bra>0 .and. present(mag_num_bra)) then
                call reorder_sgtos(angular_bra, num_gto_bra, mag_num_bra, 1, num_contr_bra, &
                                   num_gto_ket*num_contr_ket*num_opt, contr_ints, tmp_ints)
              else if (angular_ket>0 .and. present(mag_num_ket)) then
                call reorder_sgtos(angular_ket, num_gto_ket, mag_num_ket,    &
                                   num_gto_bra*num_contr_bra, num_contr_ket, &
                                   num_opt, contr_ints, tmp_ints)
              end if
            end if
          ! CGTOs
          else
            if ((angular_bra==0 .and. angular_ket==0) .or. &
                .not.(present(powers_bra) .or. present(powers_ket))) then
              call contr_cgto_carmom(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                     exponent_bra, num_contr_bra, contr_coef_bra,   &
                                     idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                     exponent_ket, num_contr_ket, contr_coef_ket,   &
                                     order_elec, MAX_IDX_NON, origin_London_PF,     &
                                     scal_const*0.5_REALK, 1, p_order_geo_bra,      &
                                     p_order_geo_ket, p_order_geo_mom,              &
                                     p_num_cents, p_idx_cent, p_order_cent,         &
                                     num_gto_bra, num_gto_ket, num_opt, tmp_ints)
            ! reorders integrals if required
            else
              call contr_cgto_carmom(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                     exponent_bra, num_contr_bra, contr_coef_bra,   &
                                     idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                     exponent_ket, num_contr_ket, contr_coef_ket,   &
                                     order_elec, MAX_IDX_NON, origin_London_PF,     &
                                     scal_const*0.5_REALK, 1, p_order_geo_bra,      &
                                     p_order_geo_ket, p_order_geo_mom,              &
                                     p_num_cents, p_idx_cent, p_order_cent,         &
                                     num_gto_bra, num_gto_ket, num_opt, contr_ints)
              if (angular_bra>0 .and. present(powers_bra) .and. &
                  angular_ket>0 .and. present(powers_ket)) then
                call reorder_cgto_ints(angular_bra, num_gto_bra, powers_bra,  &
                                       angular_ket, num_gto_ket, powers_ket,  &
                                       num_contr_bra, num_contr_ket, num_opt, &
                                       contr_ints, tmp_ints)
              else if (angular_bra>0 .and. present(powers_bra)) then
                call reorder_cgtos(angular_bra, num_gto_bra, powers_bra, 1, num_contr_bra, &
                                   num_gto_ket*num_contr_ket*num_opt, contr_ints, tmp_ints)
              else if (angular_ket>0 .and. present(powers_ket)) then
                call reorder_cgtos(angular_ket, num_gto_ket, powers_ket,     &
                                   num_gto_bra*num_contr_bra, num_contr_ket, &
                                   num_opt, contr_ints, tmp_ints)
              end if
            end if
          end if
          ! sets the first order partial magnetic derivatives
          call LondonAOGetGaugeOrigin(info_LAO, gauge_origin)
          if (order_mag_bra==1 .and. order_mag_ket==0) then
            do iopt = 1, num_opt-2, 3
              contr_ints(:,:,:,:,iopt) = (coord_bra(2)-gauge_origin(2))*tmp_ints(:,:,:,:,iopt+2) &
                                       - (coord_bra(3)-gauge_origin(3))*tmp_ints(:,:,:,:,iopt+1)
              contr_ints(:,:,:,:,iopt+1) = (coord_bra(3)-gauge_origin(3))*tmp_ints(:,:,:,:,iopt) &
                                         - (coord_bra(1)-gauge_origin(1))*tmp_ints(:,:,:,:,iopt+2)
              contr_ints(:,:,:,:,iopt+2) = (coord_bra(1)-gauge_origin(1))*tmp_ints(:,:,:,:,iopt+1) &
                                         - (coord_bra(2)-gauge_origin(2))*tmp_ints(:,:,:,:,iopt)
            end do
          else if (order_mag_bra==0 .and. order_mag_ket==1) then
            do iopt = 1, num_opt-2, 3
              contr_ints(:,:,:,:,iopt) = (coord_ket(3)-gauge_origin(3))*tmp_ints(:,:,:,:,iopt+1) &
                                       - (coord_ket(2)-gauge_origin(2))*tmp_ints(:,:,:,:,iopt+2)
              contr_ints(:,:,:,:,iopt+1) = (coord_ket(1)-gauge_origin(1))*tmp_ints(:,:,:,:,iopt+2) &
                                         - (coord_ket(3)-gauge_origin(3))*tmp_ints(:,:,:,:,iopt)
              contr_ints(:,:,:,:,iopt+2) = (coord_ket(2)-gauge_origin(2))*tmp_ints(:,:,:,:,iopt) &
                                         - (coord_ket(1)-gauge_origin(1))*tmp_ints(:,:,:,:,iopt+1)
            end do
          else
            do iopt = 1, num_opt-2, 3
              contr_ints(:,:,:,:,iopt) = (coord_bra(2)-coord_ket(2))*tmp_ints(:,:,:,:,iopt+2) &
                                       - (coord_bra(3)-coord_ket(3))*tmp_ints(:,:,:,:,iopt+1)
              contr_ints(:,:,:,:,iopt+1) = (coord_bra(3)-coord_ket(3))*tmp_ints(:,:,:,:,iopt) &
                                         - (coord_bra(1)-coord_ket(1))*tmp_ints(:,:,:,:,iopt+2)
              contr_ints(:,:,:,:,iopt+2) = (coord_bra(1)-coord_ket(1))*tmp_ints(:,:,:,:,iopt+1) &
                                         - (coord_bra(2)-coord_ket(2))*tmp_ints(:,:,:,:,iopt)
            end do
          end if
          deallocate(tmp_ints)
        !FIXME: quick implementation of CM-1 integrals
        else if ((order_mag_total==1 .and. order_mom==1 .and.   &
                  order_mag_bra==0 .and. order_mag_ket==0 .and. &
                  order_elec==0 .and. p_order_geo_bra==0 .and.  &
                  p_order_geo_ket==0)) then
          call LondonAOGetLPFOrigin(info_LAO, origin_London_PF)
          allocate(tmp_ints(num_gto_bra,num_contr_bra, &
                            num_gto_ket,num_contr_ket,9), stat=ierr)
          if (ierr/=0)                                                     &
            call error_stop("IntGetCARMOM", "failed to allocate tmp_ints", &
                            num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket*9)
          ! SGTOs
          if (p_spher_gto) then
            if ((angular_bra==0 .and. angular_ket==0) .or. &
                .not.(present(mag_num_bra) .or. present(mag_num_ket))) then
              call contr_sgto_carmom(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                     exponent_bra, num_contr_bra, contr_coef_bra,   &
                                     idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                     exponent_ket, num_contr_ket, contr_coef_ket,   &
                                     order_elec, MAX_IDX_NON, origin_London_PF,     &
                                     scal_const*0.5_REALK, 1, p_order_geo_bra,      &
                                     p_order_geo_ket, p_order_geo_mom,              &
                                     p_num_cents, p_idx_cent, p_order_cent,         &
                                     num_gto_bra, num_gto_ket, 3, tmp_ints)
              call contr_sgto_carmom(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                     exponent_bra, num_contr_bra, contr_coef_bra,   &
                                     idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                     exponent_ket, num_contr_ket, contr_coef_ket,   &
                                     order_elec, MAX_IDX_NON, origin_London_PF,     &
                                     scal_const*0.5_REALK, 2, p_order_geo_bra,      &
                                     p_order_geo_ket, p_order_geo_mom,              &
                                     p_num_cents, p_idx_cent, p_order_cent,         &
                                     num_gto_bra, num_gto_ket, 6, tmp_ints(:,:,:,:,4:9))
            ! reorders integrals if required
            else
              call contr_sgto_carmom(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                     exponent_bra, num_contr_bra, contr_coef_bra,   &
                                     idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                     exponent_ket, num_contr_ket, contr_coef_ket,   &
                                     order_elec, MAX_IDX_NON, origin_London_PF,     &
                                     scal_const*0.5_REALK, 1, p_order_geo_bra,      &
                                     p_order_geo_ket, p_order_geo_mom,              &
                                     p_num_cents, p_idx_cent, p_order_cent,         &
                                     num_gto_bra, num_gto_ket, 3, contr_ints)
              call contr_sgto_carmom(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                     exponent_bra, num_contr_bra, contr_coef_bra,   &
                                     idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                     exponent_ket, num_contr_ket, contr_coef_ket,   &
                                     order_elec, MAX_IDX_NON, origin_London_PF,     &
                                     scal_const*0.5_REALK, 2, p_order_geo_bra,      &
                                     p_order_geo_ket, p_order_geo_mom,              &
                                     p_num_cents, p_idx_cent, p_order_cent,         &
                                     num_gto_bra, num_gto_ket, 6, contr_ints(:,:,:,:,4:9))
              if (angular_bra>0 .and. present(mag_num_bra) .and. &
                  angular_ket>0 .and. present(mag_num_ket)) then
                call reorder_sgto_ints(angular_bra, num_gto_bra, mag_num_bra, &
                                       angular_ket, num_gto_ket, mag_num_ket, &
                                       num_contr_bra, num_contr_ket, 9,       &
                                       contr_ints, tmp_ints)
              else if (angular_bra>0 .and. present(mag_num_bra)) then
                call reorder_sgtos(angular_bra, num_gto_bra, mag_num_bra, 1, num_contr_bra, &
                                   num_gto_ket*num_contr_ket*9, contr_ints, tmp_ints)
              else if (angular_ket>0 .and. present(mag_num_ket)) then
                call reorder_sgtos(angular_ket, num_gto_ket, mag_num_ket,    &
                                   num_gto_bra*num_contr_bra, num_contr_ket, &
                                   9, contr_ints, tmp_ints)
              end if
            end if
          ! CGTOs
          else
            if ((angular_bra==0 .and. angular_ket==0) .or. &
                .not.(present(powers_bra) .or. present(powers_ket))) then
              call contr_cgto_carmom(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                     exponent_bra, num_contr_bra, contr_coef_bra,   &
                                     idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                     exponent_ket, num_contr_ket, contr_coef_ket,   &
                                     order_elec, MAX_IDX_NON, origin_London_PF,     &
                                     scal_const*0.5_REALK, 1, p_order_geo_bra,      &
                                     p_order_geo_ket, p_order_geo_mom,              &
                                     p_num_cents, p_idx_cent, p_order_cent,         &
                                     num_gto_bra, num_gto_ket, 3, tmp_ints)
              call contr_cgto_carmom(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                     exponent_bra, num_contr_bra, contr_coef_bra,   &
                                     idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                     exponent_ket, num_contr_ket, contr_coef_ket,   &
                                     order_elec, MAX_IDX_NON, origin_London_PF,     &
                                     scal_const*0.5_REALK, 2, p_order_geo_bra,      &
                                     p_order_geo_ket, p_order_geo_mom,              &
                                     p_num_cents, p_idx_cent, p_order_cent,         &
                                     num_gto_bra, num_gto_ket, 6, tmp_ints(:,:,:,:,4:9))
            ! reorders integrals if required
            else
              call contr_cgto_carmom(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                     exponent_bra, num_contr_bra, contr_coef_bra,   &
                                     idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                     exponent_ket, num_contr_ket, contr_coef_ket,   &
                                     order_elec, MAX_IDX_NON, origin_London_PF,     &
                                     scal_const*0.5_REALK, 1, p_order_geo_bra,      &
                                     p_order_geo_ket, p_order_geo_mom,              &
                                     p_num_cents, p_idx_cent, p_order_cent,         &
                                     num_gto_bra, num_gto_ket, 3, contr_ints)
              call contr_cgto_carmom(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                     exponent_bra, num_contr_bra, contr_coef_bra,   &
                                     idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                     exponent_ket, num_contr_ket, contr_coef_ket,   &
                                     order_elec, MAX_IDX_NON, origin_London_PF,     &
                                     scal_const*0.5_REALK, 2, p_order_geo_bra,      &
                                     p_order_geo_ket, p_order_geo_mom,              &
                                     p_num_cents, p_idx_cent, p_order_cent,         &
                                     num_gto_bra, num_gto_ket, 6, contr_ints(:,:,:,:,4:9))
              if (angular_bra>0 .and. present(powers_bra) .and. &
                  angular_ket>0 .and. present(powers_ket)) then
                call reorder_cgto_ints(angular_bra, num_gto_bra, powers_bra,  &
                                       angular_ket, num_gto_ket, powers_ket,  &
                                       num_contr_bra, num_contr_ket, 9,       &
                                       contr_ints, tmp_ints)
              else if (angular_bra>0 .and. present(powers_bra)) then
                call reorder_cgtos(angular_bra, num_gto_bra, powers_bra, 1, num_contr_bra, &
                                   num_gto_ket*num_contr_ket*9, contr_ints, tmp_ints)
              else if (angular_ket>0 .and. present(powers_ket)) then
                call reorder_cgtos(angular_ket, num_gto_ket, powers_ket,     &
                                   num_gto_bra*num_contr_bra, num_contr_ket, &
                                   9, contr_ints, tmp_ints)
              end if
            end if
          end if
          ! sets the first order magnetic derivatives
          origin_London_PF(1) = origin_London_PF(1)-dipole_origin(1)
          origin_London_PF(2) = origin_London_PF(2)-dipole_origin(2)
          origin_London_PF(3) = origin_London_PF(3)-dipole_origin(3)
          ! x,Bx
          contr_ints(:,:,:,:,1) = (coord_bra(2)-coord_ket(2)) &
              * (tmp_ints(:,:,:,:,7)+origin_London_PF(1)*tmp_ints(:,:,:,:,3)) &
              - (coord_bra(3)-coord_ket(3)) &
              * (tmp_ints(:,:,:,:,5)+origin_London_PF(1)*tmp_ints(:,:,:,:,2))
          ! y,Bx
          contr_ints(:,:,:,:,2) = (coord_bra(2)-coord_ket(2)) &
              * (tmp_ints(:,:,:,:,8)+origin_London_PF(2)*tmp_ints(:,:,:,:,3)) &
              - (coord_bra(3)-coord_ket(3)) &
              * (tmp_ints(:,:,:,:,6)+origin_London_PF(2)*tmp_ints(:,:,:,:,2))
          ! z,Bx
          contr_ints(:,:,:,:,3) = (coord_bra(2)-coord_ket(2)) &
              * (tmp_ints(:,:,:,:,9)+origin_London_PF(3)*tmp_ints(:,:,:,:,3)) &
              - (coord_bra(3)-coord_ket(3)) &
              * (tmp_ints(:,:,:,:,8)+origin_London_PF(3)*tmp_ints(:,:,:,:,2))
          ! x,By
          contr_ints(:,:,:,:,4) = (coord_bra(3)-coord_ket(3)) &
              * (tmp_ints(:,:,:,:,4)+origin_London_PF(1)*tmp_ints(:,:,:,:,1)) &
              - (coord_bra(1)-coord_ket(1)) &
              * (tmp_ints(:,:,:,:,7)+origin_London_PF(1)*tmp_ints(:,:,:,:,3))
          ! y,By
          contr_ints(:,:,:,:,5) = (coord_bra(3)-coord_ket(3)) &
              * (tmp_ints(:,:,:,:,5)+origin_London_PF(2)*tmp_ints(:,:,:,:,1)) &
              - (coord_bra(1)-coord_ket(1)) &
              * (tmp_ints(:,:,:,:,8)+origin_London_PF(2)*tmp_ints(:,:,:,:,3))
          ! z,By
          contr_ints(:,:,:,:,6) = (coord_bra(3)-coord_ket(3)) &
              * (tmp_ints(:,:,:,:,7)+origin_London_PF(3)*tmp_ints(:,:,:,:,1)) &
              - (coord_bra(1)-coord_ket(1)) &
              * (tmp_ints(:,:,:,:,9)+origin_London_PF(3)*tmp_ints(:,:,:,:,3))
          ! x,Bz
          contr_ints(:,:,:,:,7) = (coord_bra(1)-coord_ket(1)) &
              * (tmp_ints(:,:,:,:,5)+origin_London_PF(1)*tmp_ints(:,:,:,:,2)) &
              - (coord_bra(2)-coord_ket(2)) &
              * (tmp_ints(:,:,:,:,4)+origin_London_PF(1)*tmp_ints(:,:,:,:,1))
          ! y,Bz
          contr_ints(:,:,:,:,8) = (coord_bra(1)-coord_ket(1)) &
              * (tmp_ints(:,:,:,:,6)+origin_London_PF(2)*tmp_ints(:,:,:,:,2)) &
              - (coord_bra(2)-coord_ket(2)) &
              * (tmp_ints(:,:,:,:,5)+origin_London_PF(2)*tmp_ints(:,:,:,:,1))
          ! z,Bz
          contr_ints(:,:,:,:,9) = (coord_bra(1)-coord_ket(1)) &
              * (tmp_ints(:,:,:,:,8)+origin_London_PF(3)*tmp_ints(:,:,:,:,2)) &
              - (coord_bra(2)-coord_ket(2)) &
              * (tmp_ints(:,:,:,:,7)+origin_London_PF(3)*tmp_ints(:,:,:,:,1))
          deallocate(tmp_ints)
        else
          call error_stop("IntGetCARMOM", "LAO is not implemented", -1)
        end if
      ! non-LAOs
      else
        contr_ints = 0.0
      end if
    else
      ! SGTOs
      if (p_spher_gto) then
        if ((angular_bra==0 .and. angular_ket==0) .or. &
            .not.(present(mag_num_bra) .or. present(mag_num_ket))) then
          call contr_sgto_carmom(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                 exponent_bra, num_contr_bra, contr_coef_bra,   &
                                 idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                 exponent_ket, num_contr_ket, contr_coef_ket,   &
                                 order_elec, idx_diporg, dipole_origin,         &
                                 scal_const, order_mom, p_order_geo_bra,        &
                                 p_order_geo_ket, p_order_geo_mom,              &
                                 p_num_cents, p_idx_cent, p_order_cent,         &
                                 num_gto_bra, num_gto_ket, num_opt, contr_ints)
        ! reorders integrals if required
        else
          allocate(tmp_ints(num_gto_bra,num_contr_bra, &
                            num_gto_ket,num_contr_ket,num_opt), stat=ierr)
          if (ierr/=0)                                                     &
            call error_stop("IntGetCARMOM", "failed to allocate tmp_ints", &
                            num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket*num_opt)
          call contr_sgto_carmom(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                 exponent_bra, num_contr_bra, contr_coef_bra,   &
                                 idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                 exponent_ket, num_contr_ket, contr_coef_ket,   &
                                 order_elec, idx_diporg, dipole_origin,         &
                                 scal_const, order_mom, p_order_geo_bra,        &
                                 p_order_geo_ket, p_order_geo_mom,              &
                                 p_num_cents, p_idx_cent, p_order_cent,         &
                                 num_gto_bra, num_gto_ket, num_opt, tmp_ints)
          if (angular_bra>0 .and. present(mag_num_bra) .and. &
              angular_ket>0 .and. present(mag_num_ket)) then
            call reorder_sgto_ints(angular_bra, num_gto_bra, mag_num_bra, &
                                   angular_ket, num_gto_ket, mag_num_ket, &
                                   num_contr_bra, num_contr_ket, num_opt, &
                                   tmp_ints, contr_ints)
          else if (angular_bra>0 .and. present(mag_num_bra)) then
            call reorder_sgtos(angular_bra, num_gto_bra, mag_num_bra, 1, num_contr_bra, &
                               num_gto_ket*num_contr_ket*num_opt, tmp_ints, contr_ints)
          else if (angular_ket>0 .and. present(mag_num_ket)) then
            call reorder_sgtos(angular_ket, num_gto_ket, mag_num_ket,    &
                               num_gto_bra*num_contr_bra, num_contr_ket, &
                               num_opt, tmp_ints, contr_ints)
          end if
          deallocate(tmp_ints)
        end if
      ! CGTOs
      else
        if ((angular_bra==0 .and. angular_ket==0) .or. &
            .not.(present(powers_bra) .or. present(powers_ket))) then
          call contr_cgto_carmom(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                 exponent_bra, num_contr_bra, contr_coef_bra,   &
                                 idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                 exponent_ket, num_contr_ket, contr_coef_ket,   &
                                 order_elec, idx_diporg, dipole_origin,         &
                                 scal_const, order_mom, p_order_geo_bra,        &
                                 p_order_geo_ket, p_order_geo_mom,              &
                                 p_num_cents, p_idx_cent, p_order_cent,         &
                                 num_gto_bra, num_gto_ket, num_opt, contr_ints)
        ! reorders integrals if required
        else
          allocate(tmp_ints(num_gto_bra,num_contr_bra, &
                            num_gto_ket,num_contr_ket,num_opt), stat=ierr)
          if (ierr/=0)                                                     &
            call error_stop("IntGetCARMOM", "failed to allocate tmp_ints", &
                            num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket*num_opt)
          call contr_cgto_carmom(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                 exponent_bra, num_contr_bra, contr_coef_bra,   &
                                 idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                 exponent_ket, num_contr_ket, contr_coef_ket,   &
                                 order_elec, idx_diporg, dipole_origin,         &
                                 scal_const, order_mom, p_order_geo_bra,        &
                                 p_order_geo_ket, p_order_geo_mom,              &
                                 p_num_cents, p_idx_cent, p_order_cent,         &
                                 num_gto_bra, num_gto_ket, num_opt, tmp_ints)
          if (angular_bra>0 .and. present(powers_bra) .and. &
              angular_ket>0 .and. present(powers_ket)) then
            call reorder_cgto_ints(angular_bra, num_gto_bra, powers_bra,  &
                                   angular_ket, num_gto_ket, powers_ket,  &
                                   num_contr_bra, num_contr_ket, num_opt, &
                                   tmp_ints, contr_ints)
          else if (angular_bra>0 .and. present(powers_bra)) then
            call reorder_cgtos(angular_bra, num_gto_bra, powers_bra, 1, num_contr_bra, &
                               num_gto_ket*num_contr_ket*num_opt, tmp_ints, contr_ints)
          else if (angular_ket>0 .and. present(powers_ket)) then
            call reorder_cgtos(angular_ket, num_gto_ket, powers_ket,     &
                               num_gto_bra*num_contr_bra, num_contr_ket, &
                               num_opt, tmp_ints, contr_ints)
          end if
          deallocate(tmp_ints)
        end if
      end if
    end if
    deallocate(p_idx_cent)
    deallocate(p_order_cent)
  end subroutine IntGetCARMOM

  !> \brief calculates the nuclear attraction potential integrals using contracted
  !>        Gaussian type orbitals (GTOs)
  !> \author Bin Gao
  !> \date 2011-12-12
  !> \param idx_bra is the atomic index of bra center
  !> \param coord_bra contains the coordinates of bra center
  !> \param angular_bra is the angular number of bra center (s=0, p=1, d=2, ...)
  !> \param num_prim_bra is the number of primitive Gaussians of bra center
  !> \param exponent_bra contains the exponents of primitive Gaussians of bra center
  !> \param num_contr_bra is the number of contractions of bra center
  !> \param contr_coef_bra contains the contraction coefficients of bra center
  !> \param idx_ket is the atomic index of ket center
  !> \param coord_ket contains the coordinates of ket center
  !> \param angular_ket is the angular number of ket center
  !> \param num_prim_ket is the number of primitive Gaussians of ket center
  !> \param exponent_ket contains the exponents of primitive Gaussians of ket center
  !> \param num_contr_ket is the number of contractions of ket center
  !> \param contr_coef_ket contains the contraction coefficients of ket center
  !> \param spher_gto indicates if using spherical GTOs, otherwise Cartesian GTOs
  !> \param info_LAO contains the information of London atomic orbital
  !> \param order_elec is the order of electronic derivatives
  !> \param idx_nucorg is the atomic center of nuclear potential origin (<1 for non-atomic center)
  !> \param nucpot_origin contains the coordinates of nuclear potential origin
  !> \param idx_diporg is the atomic center of dipole origin (<1 for non-atomic center)
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scale constant for nuclear attraction potential operators
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_mag_bra is the order of magnetic derivatives on bra center
  !> \param order_mag_ket is the order of magnetic derivatives on ket center
  !> \param order_mag_total is the order of total magnetic derivatives
  !> \param order_ram_bra is the order of derivatives w.r.t. total rotational angular momentum on bra center
  !> \param order_ram_ket is the order of derivatives w.r.t. total rotational angular momentum on ket center
  !> \param order_ram_total is the order of total derivatives w.r.t. total rotational angular momentum
  !> \param order_geo_bra is the order of geometric derivatives with respect to bra center
  !> \param order_geo_ket is the order of geometric derivatives with respect to ket center
  !> \param order_geo_pot is the order of geometric derivatives on nuclear attraction potential origin
  !> \param order_geo_mom is the order of geometric derivatives on dipole origin
  !> \param nary_tree_total contains the information of N-ary tree for total geometric derivatives
  !> \param num_gto_bra is the number of spherical/Cartesian GTOs on bra center
  !> \param num_gto_ket is the number of spherical/Cartesian GTOs on ket center
  !> \param num_opt is the number of operators including derivatives
  !> \param mag_num_bra contains the magnetic numbers of spherical GTOs on bra center
  !> \param mag_num_ket contains the magnetic numbers of spherical GTOs on ket center
  !> \param powers_bra contains the Cartesian powers of Cartesian GTOs on bra center
  !> \param powers_ket contains the Cartesian powers of Cartesian GTOs on ket center
  !> \return contr_ints contains the calculated contracted integrals
  subroutine IntGetNUCPOT(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                          exponent_bra, num_contr_bra, contr_coef_bra,   &
                          idx_ket, coord_ket, angular_ket, num_prim_ket, &
                          exponent_ket, num_contr_ket, contr_coef_ket,   &
                          spher_gto, info_LAO, order_elec, idx_nucorg,   &
                          nucpot_origin, idx_diporg, dipole_origin,      &
                          scal_const, order_mom,                         &
                          order_mag_bra, order_mag_ket, order_mag_total, &
                          order_ram_bra, order_ram_ket, order_ram_total, &
                          order_geo_bra, order_geo_ket,                  &
                          order_geo_pot, order_geo_mom,                  &
                          nary_tree_total, num_gto_bra, num_gto_ket,     &
                          num_opt, contr_ints, mag_num_bra, mag_num_ket, &
                          powers_bra, powers_ket)
    integer, intent(in) :: idx_bra
    real(REALK), intent(in) :: coord_bra(3)
    integer, intent(in) :: angular_bra
    integer, intent(in) :: num_prim_bra
    real(REALK), intent(in) :: exponent_bra(num_prim_bra)
    integer, intent(in) :: num_contr_bra
    real(REALK), intent(in) :: contr_coef_bra(num_contr_bra,num_prim_bra)
    integer, intent(in) :: idx_ket
    real(REALK), intent(in) :: coord_ket(3)
    integer, intent(in) :: angular_ket
    integer, intent(in) :: num_prim_ket
    real(REALK), intent(in) :: exponent_ket(num_prim_ket)
    integer, intent(in) :: num_contr_ket
    real(REALK), intent(in) :: contr_coef_ket(num_contr_ket,num_prim_ket)
    logical, optional, intent(in) :: spher_gto
    type(london_ao_t), intent(in) :: info_LAO
    integer, intent(in) :: order_elec
    integer, intent(in) :: idx_nucorg
    real(REALK), intent(in) :: nucpot_origin(3)
    integer, intent(in) :: idx_diporg
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_mag_bra
    integer, intent(in) :: order_mag_ket
    integer, intent(in) :: order_mag_total
    integer, intent(in) :: order_ram_bra
    integer, intent(in) :: order_ram_ket
    integer, intent(in) :: order_ram_total
    integer, optional, intent(in) :: order_geo_bra
    integer, optional, intent(in) :: order_geo_ket
    integer, optional, intent(in) :: order_geo_pot
    integer, optional, intent(in) :: order_geo_mom
    type(nary_tree_t), optional, intent(in) :: nary_tree_total
    integer, intent(in) :: num_gto_bra
    integer, intent(in) :: num_gto_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(out) :: contr_ints(num_gto_bra,num_contr_bra, &
                                           num_gto_ket,num_contr_ket,num_opt)
    integer, optional, intent(in) :: mag_num_bra(num_gto_bra)
    integer, optional, intent(in) :: mag_num_ket(num_gto_ket)
    integer, optional, intent(in) :: powers_bra(3,num_gto_bra)
    integer, optional, intent(in) :: powers_ket(3,num_gto_ket)
#include "max_idx_non.h"
    real(REALK) gauge_origin(3)                      !gauge origin of the magnetic vector potential
    real(REALK) origin_London_PF(3)                  !origin of the London phase factor
    integer iopt                                     !incremental recorder over operators
    logical p_spher_gto                              !arguments for Gen1Int (local)
    integer p_order_geo_bra
    integer p_order_geo_ket
    integer p_order_geo_mom
    integer p_order_geo_pot
    integer p_num_cents
    integer, allocatable :: p_idx_cent(:)
    integer, allocatable :: p_order_cent(:)
    integer icent
    real(REALK), allocatable :: tmp_ints(:,:,:,:,:)  !contracted integrals from Gen1Int
    real(REALK), allocatable :: ro_ints(:,:,:,:,:)
    integer ierr                                     !error information
    ! sets the arguments for Gen1Int (local)
    if (present(spher_gto)) then
      p_spher_gto = spher_gto
    else
      p_spher_gto = .true.
    end if
    if (present(order_geo_bra)) then
      p_order_geo_bra = order_geo_bra
    else
      p_order_geo_bra = 0
    end if
    if (present(order_geo_ket)) then
      p_order_geo_ket = order_geo_ket
    else
      p_order_geo_ket = 0
    end if
    if (present(order_geo_pot)) then
      p_order_geo_pot = order_geo_pot
    else
      p_order_geo_pot = 0
    end if
    if (present(order_geo_mom)) then
      p_order_geo_mom = order_geo_mom
    else
      p_order_geo_mom = 0
    end if
    if (present(nary_tree_total)) then
      p_num_cents = nary_tree_total%wt_node(nary_tree_total%order_geo)
    else
      p_num_cents = 0
    end if
    if (p_num_cents>0) then
      allocate(p_idx_cent(p_num_cents), stat=ierr)
      if (ierr/=0) &
        call error_stop("IntGetNUCPOT", "failed to allocate p_idx_cent", p_num_cents)
      allocate(p_order_cent(p_num_cents), stat=ierr)
      if (ierr/=0) &
        call error_stop("IntGetNUCPOT", "failed to allocate p_order_cent", p_num_cents)
      if (allocated(nary_tree_total%idx_atoms)) then
        do icent = 1, p_num_cents
          p_idx_cent(icent) = nary_tree_total%idx_atoms(nary_tree_total%idx_cent(icent))
        end do
      else
        p_idx_cent = nary_tree_total%idx_cent(1:p_num_cents)
      end if
      p_order_cent = nary_tree_total%order_cent(1:p_num_cents)
    else
      allocate(p_idx_cent(1), stat=ierr)
      if (ierr/=0) call error_stop("IntGetNUCPOT", "failed to allocate p_idx_cent", 1)
      allocate(p_order_cent(1), stat=ierr)
      if (ierr/=0) call error_stop("IntGetNUCPOT", "failed to allocate p_order_cent", 1)
      p_idx_cent = 0
      p_order_cent = 0
    end if
    ! magnetic derivatives or derivatives with respect to rotational angular momentum
    if (order_mag_total>0 .or. order_mag_bra>0 .or. order_mag_ket>0 .or. &
        order_ram_total>0 .or. order_ram_bra>0 .or. order_ram_ket>0) then
      ! London atomic orbitals
      if (LondonAOUsed(info_LAO)) then
        if ((order_mag_total==0 .and. order_mom==0 .and.       &
             ((order_mag_bra==1 .and. order_mag_ket==0) .or.   &
              (order_mag_bra==0 .and. order_mag_ket==1)) .and. &
             order_elec==0 .and. p_order_geo_bra==0 .and.      &
             p_order_geo_ket==0 .and. p_num_cents<2) .or.      &
            (order_mag_total==1 .and. order_mom==0 .and.       &
             order_mag_bra==0 .and. order_mag_ket==0 .and.     &
             order_elec==0 .and. p_order_geo_bra==0 .and.      &
             p_order_geo_ket==0 .and. p_num_cents<2)) then
          if ((order_mag_total==1 .and. idx_bra/=idx_ket) .or. order_mag_total==0) then
            call LondonAOGetLPFOrigin(info_LAO, origin_London_PF)
            allocate(tmp_ints(num_gto_bra,num_contr_bra, &
                              num_gto_ket,num_contr_ket,num_opt), stat=ierr)
            if (ierr/=0)                                                     &
              call error_stop("IntGetNUCPOT", "failed to allocate tmp_ints", &
                              num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket*num_opt)
            ! SGTOs
            if (p_spher_gto) then
              if ((angular_bra==0 .and. angular_ket==0) .or. &
                  .not.(present(mag_num_bra) .or. present(mag_num_ket))) then
                call contr_sgto_nucpot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                       exponent_bra, num_contr_bra, contr_coef_bra,   &
                                       idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                       exponent_ket, num_contr_ket, contr_coef_ket,   &
                                       order_elec, idx_nucorg, nucpot_origin,         &
                                       MAX_IDX_NON, origin_London_PF,                 &
                                       scal_const*0.5_REALK, 1,                       &
                                       p_order_geo_bra, p_order_geo_ket,              &
                                       p_order_geo_pot, p_order_geo_mom,              &
                                       p_num_cents, p_idx_cent, p_order_cent,         &
                                       num_gto_bra, num_gto_ket, num_opt, tmp_ints)
              ! reorders integrals if required
              else
                call contr_sgto_nucpot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                       exponent_bra, num_contr_bra, contr_coef_bra,   &
                                       idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                       exponent_ket, num_contr_ket, contr_coef_ket,   &
                                       order_elec, idx_nucorg, nucpot_origin,         &
                                       MAX_IDX_NON, origin_London_PF,                 &
                                       scal_const*0.5_REALK, 1,                       &
                                       p_order_geo_bra, p_order_geo_ket,              &
                                       p_order_geo_pot, p_order_geo_mom,              &
                                       p_num_cents, p_idx_cent, p_order_cent,         &
                                       num_gto_bra, num_gto_ket, num_opt, contr_ints)
                if (angular_bra>0 .and. present(mag_num_bra) .and. &
                    angular_ket>0 .and. present(mag_num_ket)) then
                  call reorder_sgto_ints(angular_bra, num_gto_bra, mag_num_bra, &
                                         angular_ket, num_gto_ket, mag_num_ket, &
                                         num_contr_bra, num_contr_ket, num_opt, &
                                         contr_ints, tmp_ints)
                else if (angular_bra>0 .and. present(mag_num_bra)) then
                  call reorder_sgtos(angular_bra, num_gto_bra, mag_num_bra, 1, num_contr_bra, &
                                     num_gto_ket*num_contr_ket*num_opt, contr_ints, tmp_ints)
                else if (angular_ket>0 .and. present(mag_num_ket)) then
                  call reorder_sgtos(angular_ket, num_gto_ket, mag_num_ket,    &
                                     num_gto_bra*num_contr_bra, num_contr_ket, &
                                     num_opt, contr_ints, tmp_ints)
                end if
              end if
            ! CGTOs
            else
              if ((angular_bra==0 .and. angular_ket==0) .or. &
                  .not.(present(powers_bra) .or. present(powers_ket))) then
                call contr_cgto_nucpot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                       exponent_bra, num_contr_bra, contr_coef_bra,   &
                                       idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                       exponent_ket, num_contr_ket, contr_coef_ket,   &
                                       order_elec, idx_nucorg, nucpot_origin,         &
                                       MAX_IDX_NON, origin_London_PF,                 &
                                       scal_const*0.5_REALK, 1,                       &
                                       p_order_geo_bra, p_order_geo_ket,              &
                                       p_order_geo_pot, p_order_geo_mom,              &
                                       p_num_cents, p_idx_cent, p_order_cent,         &
                                       num_gto_bra, num_gto_ket, num_opt, tmp_ints)
              ! reorders integrals if required
              else
                call contr_cgto_nucpot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                       exponent_bra, num_contr_bra, contr_coef_bra,   &
                                       idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                       exponent_ket, num_contr_ket, contr_coef_ket,   &
                                       order_elec, idx_nucorg, nucpot_origin,         &
                                       MAX_IDX_NON, origin_London_PF,                 &
                                       scal_const*0.5_REALK, 1,                       &
                                       p_order_geo_bra, p_order_geo_ket,              &
                                       p_order_geo_pot, p_order_geo_mom,              &
                                       p_num_cents, p_idx_cent, p_order_cent,         &
                                       num_gto_bra, num_gto_ket, num_opt, contr_ints)
                if (angular_bra>0 .and. present(powers_bra) .and. &
                    angular_ket>0 .and. present(powers_ket)) then
                  call reorder_cgto_ints(angular_bra, num_gto_bra, powers_bra,  &
                                         angular_ket, num_gto_ket, powers_ket,  &
                                         num_contr_bra, num_contr_ket, num_opt, &
                                         contr_ints, tmp_ints)
                else if (angular_bra>0 .and. present(powers_bra)) then
                  call reorder_cgtos(angular_bra, num_gto_bra, powers_bra, 1, num_contr_bra, &
                                     num_gto_ket*num_contr_ket*num_opt, contr_ints, tmp_ints)
                else if (angular_ket>0 .and. present(powers_ket)) then
                  call reorder_cgtos(angular_ket, num_gto_ket, powers_ket,     &
                                     num_gto_bra*num_contr_bra, num_contr_ket, &
                                     num_opt, contr_ints, tmp_ints)
                end if
              end if
            end if
            ! sets the first order partial magnetic derivatives
            call LondonAOGetGaugeOrigin(info_LAO, gauge_origin)
            if (order_mag_bra==1) then
              do iopt = 1, num_opt-2, 3
                contr_ints(:,:,:,:,iopt) = (coord_bra(2)-gauge_origin(2))*tmp_ints(:,:,:,:,iopt+2) &
                                         - (coord_bra(3)-gauge_origin(3))*tmp_ints(:,:,:,:,iopt+1)
                contr_ints(:,:,:,:,iopt+1) = (coord_bra(3)-gauge_origin(3))*tmp_ints(:,:,:,:,iopt) &
                                           - (coord_bra(1)-gauge_origin(1))*tmp_ints(:,:,:,:,iopt+2)
                contr_ints(:,:,:,:,iopt+2) = (coord_bra(1)-gauge_origin(1))*tmp_ints(:,:,:,:,iopt+1) &
                                           - (coord_bra(2)-gauge_origin(2))*tmp_ints(:,:,:,:,iopt)
              end do
            else if (order_mag_ket==1) then
              do iopt = 1, num_opt-2, 3
                contr_ints(:,:,:,:,iopt) = (coord_ket(3)-gauge_origin(3))*tmp_ints(:,:,:,:,iopt+1) &
                                         - (coord_ket(2)-gauge_origin(2))*tmp_ints(:,:,:,:,iopt+2)
                contr_ints(:,:,:,:,iopt+1) = (coord_ket(1)-gauge_origin(1))*tmp_ints(:,:,:,:,iopt+2) &
                                           - (coord_ket(3)-gauge_origin(3))*tmp_ints(:,:,:,:,iopt)
                contr_ints(:,:,:,:,iopt+2) = (coord_ket(2)-gauge_origin(2))*tmp_ints(:,:,:,:,iopt) &
                                           - (coord_ket(1)-gauge_origin(1))*tmp_ints(:,:,:,:,iopt+1)
              end do
            else
              do iopt = 1, num_opt-2, 3
                contr_ints(:,:,:,:,iopt) = (coord_bra(2)-coord_ket(2))*tmp_ints(:,:,:,:,iopt+2) &
                                         - (coord_bra(3)-coord_ket(3))*tmp_ints(:,:,:,:,iopt+1)
                contr_ints(:,:,:,:,iopt+1) = (coord_bra(3)-coord_ket(3))*tmp_ints(:,:,:,:,iopt) &
                                           - (coord_bra(1)-coord_ket(1))*tmp_ints(:,:,:,:,iopt+2)
                contr_ints(:,:,:,:,iopt+2) = (coord_bra(1)-coord_ket(1))*tmp_ints(:,:,:,:,iopt+1) &
                                           - (coord_bra(2)-coord_ket(2))*tmp_ints(:,:,:,:,iopt)
              end do
            end if
            deallocate(tmp_ints)
            ! mixed total geometric derivatives and magnetic derivatives
            if (p_num_cents==1) then
              if (p_order_cent(1)==1) then
                if ((p_idx_cent(1)==idx_bra .and. (order_mag_bra==1 .or. order_mag_total==1)) .or. &
                    (p_idx_cent(1)==idx_ket .and. (order_mag_ket==1 .or. order_mag_total==1))) then
                  allocate(tmp_ints(num_gto_bra,num_contr_bra, &
                                    num_gto_ket,num_contr_ket,num_opt/3), stat=ierr)
                  if (ierr/=0)                                                     &
                    call error_stop("IntGetNUCPOT", "failed to allocate tmp_ints", &
                                    num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket*num_opt/3)
                  allocate(ro_ints(num_gto_bra,num_contr_bra, &
                                    num_gto_ket,num_contr_ket,num_opt/3), stat=ierr)
                  if (ierr/=0)                                                    &
                    call error_stop("IntGetNUCPOT", "failed to allocate ro_ints", &
                                    num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket*num_opt/3)
                  ! SGTOs
                  if (p_spher_gto) then
                    if ((angular_bra==0 .and. angular_ket==0) .or. &
                        .not.(present(mag_num_bra) .or. present(mag_num_ket))) then
                      call contr_sgto_nucpot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                             exponent_bra, num_contr_bra, contr_coef_bra,   &
                                             idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                             exponent_ket, num_contr_ket, contr_coef_ket,   &
                                             order_elec, idx_nucorg, nucpot_origin,         &
                                             MAX_IDX_NON, origin_London_PF,                 &
                                             scal_const*0.5_REALK, 1,                       &
                                             p_order_geo_bra, p_order_geo_ket,              &
                                             p_order_geo_pot, p_order_geo_mom,              &
                                             0, (/0/), (/0/),                               &
                                             num_gto_bra, num_gto_ket, num_opt/3, ro_ints)
                    ! reorders integrals if required
                    else
                      call contr_sgto_nucpot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                             exponent_bra, num_contr_bra, contr_coef_bra,   &
                                             idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                             exponent_ket, num_contr_ket, contr_coef_ket,   &
                                             order_elec, idx_nucorg, nucpot_origin,         &
                                             MAX_IDX_NON, origin_London_PF,                 &
                                             scal_const*0.5_REALK, 1,                       &
                                             p_order_geo_bra, p_order_geo_ket,              &
                                             p_order_geo_pot, p_order_geo_mom,              &
                                             0, (/0/), (/0/),                               &
                                             num_gto_bra, num_gto_ket, num_opt/3, tmp_ints)
                      if (angular_bra>0 .and. present(mag_num_bra) .and. &
                          angular_ket>0 .and. present(mag_num_ket)) then
                        call reorder_sgto_ints(angular_bra, num_gto_bra, mag_num_bra,   &
                                               angular_ket, num_gto_ket, mag_num_ket,   &
                                               num_contr_bra, num_contr_ket, num_opt/3, &
                                               tmp_ints, ro_ints)
                      else if (angular_bra>0 .and. present(mag_num_bra)) then
                        call reorder_sgtos(angular_bra, num_gto_bra, mag_num_bra, 1, num_contr_bra, &
                                           num_gto_ket*num_contr_ket*num_opt/3, tmp_ints, ro_ints)
                      else if (angular_ket>0 .and. present(mag_num_ket)) then
                        call reorder_sgtos(angular_ket, num_gto_ket, mag_num_ket,    &
                                           num_gto_bra*num_contr_bra, num_contr_ket, &
                                           num_opt/3, tmp_ints, ro_ints)
                      end if
                    end if
                  ! CGTOs
                  else
                    if ((angular_bra==0 .and. angular_ket==0) .or. &
                        .not.(present(powers_bra) .or. present(powers_ket))) then
                      call contr_cgto_nucpot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                             exponent_bra, num_contr_bra, contr_coef_bra,   &
                                             idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                             exponent_ket, num_contr_ket, contr_coef_ket,   &
                                             order_elec, idx_nucorg, nucpot_origin,         &
                                             MAX_IDX_NON, origin_London_PF,                 &
                                             scal_const*0.5_REALK, 1,                       &
                                             p_order_geo_bra, p_order_geo_ket,              &
                                             p_order_geo_pot, p_order_geo_mom,              &
                                             0, (/0/), (/0/),                               &
                                             num_gto_bra, num_gto_ket, num_opt/3, ro_ints)
                    ! reorders integrals if required
                    else
                      call contr_cgto_nucpot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                             exponent_bra, num_contr_bra, contr_coef_bra,   &
                                             idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                             exponent_ket, num_contr_ket, contr_coef_ket,   &
                                             order_elec, idx_nucorg, nucpot_origin,         &
                                             MAX_IDX_NON, origin_London_PF,                 &
                                             scal_const*0.5_REALK, 1,                       &
                                             p_order_geo_bra, p_order_geo_ket,              &
                                             p_order_geo_pot, p_order_geo_mom,              &
                                             0, (/0/), (/0/),                               &
                                             num_gto_bra, num_gto_ket, num_opt/3, tmp_ints)
                      if (angular_bra>0 .and. present(powers_bra) .and. &
                          angular_ket>0 .and. present(powers_ket)) then
                        call reorder_cgto_ints(angular_bra, num_gto_bra, powers_bra,    &
                                               angular_ket, num_gto_ket, powers_ket,    &
                                               num_contr_bra, num_contr_ket, num_opt/3, &
                                               tmp_ints, ro_ints)
                      else if (angular_bra>0 .and. present(powers_bra)) then
                        call reorder_cgtos(angular_bra, num_gto_bra, powers_bra, 1, num_contr_bra, &
                                           num_gto_ket*num_contr_ket*num_opt/3, tmp_ints, ro_ints)
                      else if (angular_ket>0 .and. present(powers_ket)) then
                        call reorder_cgtos(angular_ket, num_gto_ket, powers_ket,     &
                                           num_gto_bra*num_contr_bra, num_contr_ket, &
                                           num_opt/3, tmp_ints, ro_ints)
                      end if
                    end if
                  end if
                  deallocate(tmp_ints)
                  ! first order partial magnetic derivatives on the bra center or total magnetic derivatives
                  if (p_idx_cent(1)==idx_bra) then
                      ! Bx, Gx
                      ! By, Gx
                      contr_ints(:,:,:,:,2) = contr_ints(:,:,:,:,2)-ro_ints(:,:,:,:,3)
                      ! Bz, Gx
                      contr_ints(:,:,:,:,3) = contr_ints(:,:,:,:,3)+ro_ints(:,:,:,:,2)
                      ! Bx, Gy
                      contr_ints(:,:,:,:,4) = contr_ints(:,:,:,:,4)+ro_ints(:,:,:,:,3)
                      ! By, Gy
                      ! Bz, Gy
                      contr_ints(:,:,:,:,6) = contr_ints(:,:,:,:,6)-ro_ints(:,:,:,:,1)
                      ! Bx, Gz
                      contr_ints(:,:,:,:,7) = contr_ints(:,:,:,:,7)-ro_ints(:,:,:,:,2)
                      ! By, Gz
                      contr_ints(:,:,:,:,8) = contr_ints(:,:,:,:,8)+ro_ints(:,:,:,:,1)
                      ! Bz, Gz
                  ! first order partial magnetic derivatives on the ket center or total magnetic derivatives
                  else
                      ! Bx, Gx
                      ! By, Gx
                      contr_ints(:,:,:,:,2) = contr_ints(:,:,:,:,2)+ro_ints(:,:,:,:,3)
                      ! Bz, Gx
                      contr_ints(:,:,:,:,3) = contr_ints(:,:,:,:,3)-ro_ints(:,:,:,:,2)
                      ! Bx, Gy
                      contr_ints(:,:,:,:,4) = contr_ints(:,:,:,:,4)-ro_ints(:,:,:,:,3)
                      ! By, Gy
                      ! Bz, Gy
                      contr_ints(:,:,:,:,6) = contr_ints(:,:,:,:,6)+ro_ints(:,:,:,:,1)
                      ! Bx, Gz
                      contr_ints(:,:,:,:,7) = contr_ints(:,:,:,:,7)+ro_ints(:,:,:,:,2)
                      ! By, Gz
                      contr_ints(:,:,:,:,8) = contr_ints(:,:,:,:,8)-ro_ints(:,:,:,:,1)
                      ! Bz, Gz
                  end if
                  deallocate(ro_ints)
                end if
              else
                call error_stop("IntGetNUCPOT", "LAO is not implemented", -1)
              end if
            end if
          else
            contr_ints = 0.0
          end if
        else
          call error_stop("IntGetNUCPOT", "LAO is not implemented", -1)
        end if
      ! non-LAOs
      else
        contr_ints = 0.0
      end if
    else
      ! SGTOs
      if (p_spher_gto) then
        if ((angular_bra==0 .and. angular_ket==0) .or. &
            .not.(present(mag_num_bra) .or. present(mag_num_ket))) then
          call contr_sgto_nucpot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                 exponent_bra, num_contr_bra, contr_coef_bra,   &
                                 idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                 exponent_ket, num_contr_ket, contr_coef_ket,   &
                                 order_elec, idx_nucorg, nucpot_origin,         &
                                 idx_diporg, dipole_origin,                     &
                                 scal_const, order_mom,                         &
                                 p_order_geo_bra, p_order_geo_ket,              &
                                 p_order_geo_pot, p_order_geo_mom,              &
                                 p_num_cents, p_idx_cent, p_order_cent,         &
                                 num_gto_bra, num_gto_ket, num_opt, contr_ints)
        ! reorders integrals if required
        else
          allocate(tmp_ints(num_gto_bra,num_contr_bra, &
                            num_gto_ket,num_contr_ket,num_opt), stat=ierr)
          if (ierr/=0)                                                     &
            call error_stop("IntGetNUCPOT", "failed to allocate tmp_ints", &
                            num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket*num_opt)
          call contr_sgto_nucpot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                 exponent_bra, num_contr_bra, contr_coef_bra,   &
                                 idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                 exponent_ket, num_contr_ket, contr_coef_ket,   &
                                 order_elec, idx_nucorg, nucpot_origin,         &
                                 idx_diporg, dipole_origin,                     &
                                 scal_const, order_mom,                         &
                                 p_order_geo_bra, p_order_geo_ket,              &
                                 p_order_geo_pot, p_order_geo_mom,              &
                                 p_num_cents, p_idx_cent, p_order_cent,         &
                                 num_gto_bra, num_gto_ket, num_opt, tmp_ints)
          if (angular_bra>0 .and. present(mag_num_bra) .and. &
              angular_ket>0 .and. present(mag_num_ket)) then
            call reorder_sgto_ints(angular_bra, num_gto_bra, mag_num_bra, &
                                   angular_ket, num_gto_ket, mag_num_ket, &
                                   num_contr_bra, num_contr_ket, num_opt, &
                                   tmp_ints, contr_ints)
          else if (angular_bra>0 .and. present(mag_num_bra)) then
            call reorder_sgtos(angular_bra, num_gto_bra, mag_num_bra, 1, num_contr_bra, &
                               num_gto_ket*num_contr_ket*num_opt, tmp_ints, contr_ints)
          else if (angular_ket>0 .and. present(mag_num_ket)) then
            call reorder_sgtos(angular_ket, num_gto_ket, mag_num_ket,    &
                               num_gto_bra*num_contr_bra, num_contr_ket, &
                               num_opt, tmp_ints, contr_ints)
          end if
          deallocate(tmp_ints)
        end if
      ! CGTOs
      else
        if ((angular_bra==0 .and. angular_ket==0) .or. &
            .not.(present(powers_bra) .or. present(powers_ket))) then
          call contr_cgto_nucpot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                 exponent_bra, num_contr_bra, contr_coef_bra,   &
                                 idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                 exponent_ket, num_contr_ket, contr_coef_ket,   &
                                 order_elec, idx_nucorg, nucpot_origin,         &
                                 idx_diporg, dipole_origin,                     &
                                 scal_const, order_mom,                         &
                                 p_order_geo_bra, p_order_geo_ket,              &
                                 p_order_geo_pot, p_order_geo_mom,              &
                                 p_num_cents, p_idx_cent, p_order_cent,         &
                                 num_gto_bra, num_gto_ket, num_opt, contr_ints)
        ! reorders integrals if required
        else
          allocate(tmp_ints(num_gto_bra,num_contr_bra, &
                            num_gto_ket,num_contr_ket,num_opt), stat=ierr)
          if (ierr/=0)                                                     &
            call error_stop("IntGetNUCPOT", "failed to allocate tmp_ints", &
                            num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket*num_opt)
          call contr_cgto_nucpot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                 exponent_bra, num_contr_bra, contr_coef_bra,   &
                                 idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                 exponent_ket, num_contr_ket, contr_coef_ket,   &
                                 order_elec, idx_nucorg, nucpot_origin,         &
                                 idx_diporg, dipole_origin,                     &
                                 scal_const, order_mom,                         &
                                 p_order_geo_bra, p_order_geo_ket,              &
                                 p_order_geo_pot, p_order_geo_mom,              &
                                 p_num_cents, p_idx_cent, p_order_cent,         &
                                 num_gto_bra, num_gto_ket, num_opt, tmp_ints)
          if (angular_bra>0 .and. present(powers_bra) .and. &
              angular_ket>0 .and. present(powers_ket)) then
            call reorder_cgto_ints(angular_bra, num_gto_bra, powers_bra,  &
                                   angular_ket, num_gto_ket, powers_ket,  &
                                   num_contr_bra, num_contr_ket, num_opt, &
                                   tmp_ints, contr_ints)
          else if (angular_bra>0 .and. present(powers_bra)) then
            call reorder_cgtos(angular_bra, num_gto_bra, powers_bra, 1, num_contr_bra, &
                               num_gto_ket*num_contr_ket*num_opt, tmp_ints, contr_ints)
          else if (angular_ket>0 .and. present(powers_ket)) then
            call reorder_cgtos(angular_ket, num_gto_ket, powers_ket,     &
                               num_gto_bra*num_contr_bra, num_contr_ket, &
                               num_opt, tmp_ints, contr_ints)
          end if
          deallocate(tmp_ints)
        end if
      end if
    end if
    deallocate(p_idx_cent)
    deallocate(p_order_cent)
  end subroutine IntGetNUCPOT

  !> \brief calculates the Gaussian charge potential integrals using contracted
  !>        Gaussian type orbitals (GTOs)
  !> \author Bin Gao
  !> \date 2011-12-12
  !> \param idx_bra is the atomic index of bra center
  !> \param coord_bra contains the coordinates of bra center
  !> \param angular_bra is the angular number of bra center (s=0, p=1, d=2, ...)
  !> \param num_prim_bra is the number of primitive Gaussians of bra center
  !> \param exponent_bra contains the exponents of primitive Gaussians of bra center
  !> \param num_contr_bra is the number of contractions of bra center
  !> \param contr_coef_bra contains the contraction coefficients of bra center
  !> \param idx_ket is the atomic index of ket center
  !> \param coord_ket contains the coordinates of ket center
  !> \param angular_ket is the angular number of ket center
  !> \param num_prim_ket is the number of primitive Gaussians of ket center
  !> \param exponent_ket contains the exponents of primitive Gaussians of ket center
  !> \param num_contr_ket is the number of contractions of ket center
  !> \param contr_coef_ket contains the contraction coefficients of ket center
  !> \param spher_gto indicates if using spherical GTOs, otherwise Cartesian GTOs
  !> \param info_LAO contains the information of London atomic orbital
  !> \param order_elec is the order of electronic derivatives
  !> \param idx_gauorg is the atomic center of Gaussian charge potential origin (<1 for non-atomic center)
  !> \param gaupot_origin contains the coordinates of Gaussian charge potential origin
  !> \param gaupot_expt is the exponent used in the Gaussian broadening function of the charge
  !> \param idx_diporg is the atomic center of dipole origin (<1 for non-atomic center)
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scale constant for potential operators
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_mag_bra is the order of magnetic derivatives on bra center
  !> \param order_mag_ket is the order of magnetic derivatives on ket center
  !> \param order_mag_total is the order of total magnetic derivatives
  !> \param order_ram_bra is the order of derivatives w.r.t. total rotational angular momentum on bra center
  !> \param order_ram_ket is the order of derivatives w.r.t. total rotational angular momentum on ket center
  !> \param order_ram_total is the order of total derivatives w.r.t. total rotational angular momentum
  !> \param order_geo_bra is the order of geometric derivatives with respect to bra center
  !> \param order_geo_ket is the order of geometric derivatives with respect to ket center
  !> \param order_geo_pot is the order of geometric derivatives on potential origin
  !> \param order_geo_mom is the order of geometric derivatives on dipole origin
  !> \param nary_tree_total contains the information of N-ary tree for total geometric derivatives
  !> \param num_gto_bra is the number of spherical/Cartesian GTOs on bra center
  !> \param num_gto_ket is the number of spherical/Cartesian GTOs on ket center
  !> \param num_opt is the number of operators including derivatives
  !> \param mag_num_bra contains the magnetic numbers of spherical GTOs on bra center
  !> \param mag_num_ket contains the magnetic numbers of spherical GTOs on ket center
  !> \param powers_bra contains the Cartesian powers of Cartesian GTOs on bra center
  !> \param powers_ket contains the Cartesian powers of Cartesian GTOs on ket center
  !> \return contr_ints contains the calculated contracted integrals
  subroutine IntGetGAUPOT(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                          exponent_bra, num_contr_bra, contr_coef_bra,   &
                          idx_ket, coord_ket, angular_ket, num_prim_ket, &
                          exponent_ket, num_contr_ket, contr_coef_ket,   &
                          spher_gto, info_LAO, order_elec, idx_gauorg,   &
                          gaupot_origin, gaupot_expt, idx_diporg,        &
                          dipole_origin, scal_const, order_mom,          &
                          order_mag_bra, order_mag_ket, order_mag_total, &
                          order_ram_bra, order_ram_ket, order_ram_total, &
                          order_geo_bra, order_geo_ket,                  &
                          order_geo_pot, order_geo_mom,                  &
                          nary_tree_total, num_gto_bra, num_gto_ket,     &
                          num_opt, contr_ints, mag_num_bra, mag_num_ket, &
                          powers_bra, powers_ket)
    integer, intent(in) :: idx_bra
    real(REALK), intent(in) :: coord_bra(3)
    integer, intent(in) :: angular_bra
    integer, intent(in) :: num_prim_bra
    real(REALK), intent(in) :: exponent_bra(num_prim_bra)
    integer, intent(in) :: num_contr_bra
    real(REALK), intent(in) :: contr_coef_bra(num_contr_bra,num_prim_bra)
    integer, intent(in) :: idx_ket
    real(REALK), intent(in) :: coord_ket(3)
    integer, intent(in) :: angular_ket
    integer, intent(in) :: num_prim_ket
    real(REALK), intent(in) :: exponent_ket(num_prim_ket)
    integer, intent(in) :: num_contr_ket
    real(REALK), intent(in) :: contr_coef_ket(num_contr_ket,num_prim_ket)
    logical, optional, intent(in) :: spher_gto
    type(london_ao_t), intent(in) :: info_LAO
    integer, intent(in) :: order_elec
    integer, intent(in) :: idx_gauorg
    real(REALK), intent(in) :: gaupot_origin(3)
    real(REALK), intent(in) :: gaupot_expt
    integer, intent(in) :: idx_diporg
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_mag_bra
    integer, intent(in) :: order_mag_ket
    integer, intent(in) :: order_mag_total
    integer, intent(in) :: order_ram_bra
    integer, intent(in) :: order_ram_ket
    integer, intent(in) :: order_ram_total
    integer, optional, intent(in) :: order_geo_bra
    integer, optional, intent(in) :: order_geo_ket
    integer, optional, intent(in) :: order_geo_pot
    integer, optional, intent(in) :: order_geo_mom
    type(nary_tree_t), optional, intent(in) :: nary_tree_total
    integer, intent(in) :: num_gto_bra
    integer, intent(in) :: num_gto_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(out) :: contr_ints(num_gto_bra,num_contr_bra, &
                                           num_gto_ket,num_contr_ket,num_opt)
    integer, optional, intent(in) :: mag_num_bra(num_gto_bra)
    integer, optional, intent(in) :: mag_num_ket(num_gto_ket)
    integer, optional, intent(in) :: powers_bra(3,num_gto_bra)
    integer, optional, intent(in) :: powers_ket(3,num_gto_ket)
#include "max_idx_non.h"
    real(REALK) gauge_origin(3)                      !gauge origin of the magnetic vector potential
    real(REALK) origin_London_PF(3)                  !origin of the London phase factor
    integer iopt                                     !incremental recorder over operators
    logical p_spher_gto                              !arguments for Gen1Int (local)
    integer p_order_geo_bra
    integer p_order_geo_ket
    integer p_order_geo_mom
    integer p_order_geo_pot
    integer p_num_cents
    integer, allocatable :: p_idx_cent(:)
    integer, allocatable :: p_order_cent(:)
    integer icent
    real(REALK), allocatable :: tmp_ints(:,:,:,:,:)  !contracted integrals from Gen1Int
    integer ierr                                     !error information
    ! sets the arguments for Gen1Int (local)
    if (present(spher_gto)) then
      p_spher_gto = spher_gto
    else
      p_spher_gto = .true.
    end if
    if (present(order_geo_bra)) then
      p_order_geo_bra = order_geo_bra
    else
      p_order_geo_bra = 0
    end if
    if (present(order_geo_ket)) then
      p_order_geo_ket = order_geo_ket
    else
      p_order_geo_ket = 0
    end if
    if (present(order_geo_pot)) then
      p_order_geo_pot = order_geo_pot
    else
      p_order_geo_pot = 0
    end if
    if (present(order_geo_mom)) then
      p_order_geo_mom = order_geo_mom
    else
      p_order_geo_mom = 0
    end if
    if (present(nary_tree_total)) then
      p_num_cents = nary_tree_total%wt_node(nary_tree_total%order_geo)
    else
      p_num_cents = 0
    end if
    if (p_num_cents>0) then
      allocate(p_idx_cent(p_num_cents), stat=ierr)
      if (ierr/=0) &
        call error_stop("IntGetGAUPOT", "failed to allocate p_idx_cent", p_num_cents)
      allocate(p_order_cent(p_num_cents), stat=ierr)
      if (ierr/=0) &
        call error_stop("IntGetGAUPOT", "failed to allocate p_order_cent", p_num_cents)
      if (allocated(nary_tree_total%idx_atoms)) then
        do icent = 1, p_num_cents
          p_idx_cent(icent) = nary_tree_total%idx_atoms(nary_tree_total%idx_cent(icent))
        end do
      else
        p_idx_cent = nary_tree_total%idx_cent(1:p_num_cents)
      end if
      p_order_cent = nary_tree_total%order_cent(1:p_num_cents)
    else
      allocate(p_idx_cent(1), stat=ierr)
      if (ierr/=0) call error_stop("IntGetGAUPOT", "failed to allocate p_idx_cent", 1)
      allocate(p_order_cent(1), stat=ierr)
      if (ierr/=0) call error_stop("IntGetGAUPOT", "failed to allocate p_order_cent", 1)
      p_idx_cent = 0
      p_order_cent = 0
    end if
    ! magnetic derivatives or derivatives with respect to rotational angular momentum
    if (order_mag_total>0 .or. order_mag_bra>0 .or. order_mag_ket>0 .or. &
        order_ram_total>0 .or. order_ram_bra>0 .or. order_ram_ket>0) then
      ! London atomic orbitals
      if (LondonAOUsed(info_LAO)) then
        if ((order_mag_total==0 .and. order_mom==0 .and.       &
             ((order_mag_bra==1 .and. order_mag_ket==0) .or.   &
              (order_mag_bra==0 .and. order_mag_ket==1)) .and. &
             order_elec==0 .and. p_order_geo_bra==0 .and.      &
             p_order_geo_ket==0 .and. p_num_cents==0) .or.     &
            (order_mag_total==1 .and. order_mom==0 .and.       &
             order_mag_bra==0 .and. order_mag_ket==0 .and.     &
             order_elec==0 .and. p_order_geo_bra==0 .and.      &
             p_order_geo_ket==0 .and. p_num_cents==0)) then
          call LondonAOGetLPFOrigin(info_LAO, origin_London_PF)
          allocate(tmp_ints(num_gto_bra,num_contr_bra, &
                            num_gto_ket,num_contr_ket,num_opt), stat=ierr)
          if (ierr/=0)                                                     &
            call error_stop("IntGetGAUPOT", "failed to allocate tmp_ints", &
                            num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket*num_opt)
          ! SGTOs
          if (p_spher_gto) then
            if ((angular_bra==0 .and. angular_ket==0) .or. &
                .not.(present(mag_num_bra) .or. present(mag_num_ket))) then
              call contr_sgto_gaupot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                     exponent_bra, num_contr_bra, contr_coef_bra,   &
                                     idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                     exponent_ket, num_contr_ket, contr_coef_ket,   &
                                     order_elec, idx_gauorg, gaupot_origin,         &
                                     gaupot_expt, MAX_IDX_NON, origin_London_PF,    &
                                     scal_const*0.5_REALK, 1,                       &
                                     p_order_geo_bra, p_order_geo_ket,              &
                                     p_order_geo_pot, p_order_geo_mom,              &
                                     p_num_cents, p_idx_cent, p_order_cent,         &
                                     num_gto_bra, num_gto_ket, num_opt, tmp_ints)
            ! reorders integrals if required
            else
              call contr_sgto_gaupot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                     exponent_bra, num_contr_bra, contr_coef_bra,   &
                                     idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                     exponent_ket, num_contr_ket, contr_coef_ket,   &
                                     order_elec, idx_gauorg, gaupot_origin,         &
                                     gaupot_expt, MAX_IDX_NON, origin_London_PF,    &
                                     scal_const*0.5_REALK, 1,                       &
                                     p_order_geo_bra, p_order_geo_ket,              &
                                     p_order_geo_pot, p_order_geo_mom,              &
                                     p_num_cents, p_idx_cent, p_order_cent,         &
                                     num_gto_bra, num_gto_ket, num_opt, contr_ints)
              if (angular_bra>0 .and. present(mag_num_bra) .and. &
                  angular_ket>0 .and. present(mag_num_ket)) then
                call reorder_sgto_ints(angular_bra, num_gto_bra, mag_num_bra, &
                                       angular_ket, num_gto_ket, mag_num_ket, &
                                       num_contr_bra, num_contr_ket, num_opt, &
                                       contr_ints, tmp_ints)
              else if (angular_bra>0 .and. present(mag_num_bra)) then
                call reorder_sgtos(angular_bra, num_gto_bra, mag_num_bra, 1, num_contr_bra, &
                                   num_gto_ket*num_contr_ket*num_opt, contr_ints, tmp_ints)
              else if (angular_ket>0 .and. present(mag_num_ket)) then
                call reorder_sgtos(angular_ket, num_gto_ket, mag_num_ket,    &
                                   num_gto_bra*num_contr_bra, num_contr_ket, &
                                   num_opt, contr_ints, tmp_ints)
              end if
            end if
          ! CGTOs
          else
            if ((angular_bra==0 .and. angular_ket==0) .or. &
                .not.(present(powers_bra) .or. present(powers_ket))) then
              call contr_cgto_gaupot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                     exponent_bra, num_contr_bra, contr_coef_bra,   &
                                     idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                     exponent_ket, num_contr_ket, contr_coef_ket,   &
                                     order_elec, idx_gauorg, gaupot_origin,         &
                                     gaupot_expt, MAX_IDX_NON, origin_London_PF,    &
                                     scal_const*0.5_REALK, 1,                       &
                                     p_order_geo_bra, p_order_geo_ket,              &
                                     p_order_geo_pot, p_order_geo_mom,              &
                                     p_num_cents, p_idx_cent, p_order_cent,         &
                                     num_gto_bra, num_gto_ket, num_opt, tmp_ints)
            ! reorders integrals if required
            else
              call contr_cgto_gaupot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                     exponent_bra, num_contr_bra, contr_coef_bra,   &
                                     idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                     exponent_ket, num_contr_ket, contr_coef_ket,   &
                                     order_elec, idx_gauorg, gaupot_origin,         &
                                     gaupot_expt, MAX_IDX_NON, origin_London_PF,    &
                                     scal_const*0.5_REALK, 1,                       &
                                     p_order_geo_bra, p_order_geo_ket,              &
                                     p_order_geo_pot, p_order_geo_mom,              &
                                     p_num_cents, p_idx_cent, p_order_cent,         &
                                     num_gto_bra, num_gto_ket, num_opt, contr_ints)
              if (angular_bra>0 .and. present(powers_bra) .and. &
                  angular_ket>0 .and. present(powers_ket)) then
                call reorder_cgto_ints(angular_bra, num_gto_bra, powers_bra,  &
                                       angular_ket, num_gto_ket, powers_ket,  &
                                       num_contr_bra, num_contr_ket, num_opt, &
                                       contr_ints, tmp_ints)
              else if (angular_bra>0 .and. present(powers_bra)) then
                call reorder_cgtos(angular_bra, num_gto_bra, powers_bra, 1, num_contr_bra, &
                                   num_gto_ket*num_contr_ket*num_opt, contr_ints, tmp_ints)
              else if (angular_ket>0 .and. present(powers_ket)) then
                call reorder_cgtos(angular_ket, num_gto_ket, powers_ket,     &
                                   num_gto_bra*num_contr_bra, num_contr_ket, &
                                   num_opt, contr_ints, tmp_ints)
              end if
            end if
          end if
          ! sets the first order partial magnetic derivatives
          call LondonAOGetGaugeOrigin(info_LAO, gauge_origin)
          if (order_mag_bra==1 .and. order_mag_ket==0) then
            do iopt = 1, num_opt-2, 3
              contr_ints(:,:,:,:,iopt) = (coord_bra(2)-gauge_origin(2))*tmp_ints(:,:,:,:,iopt+2) &
                                       - (coord_bra(3)-gauge_origin(3))*tmp_ints(:,:,:,:,iopt+1)
              contr_ints(:,:,:,:,iopt+1) = (coord_bra(3)-gauge_origin(3))*tmp_ints(:,:,:,:,iopt) &
                                         - (coord_bra(1)-gauge_origin(1))*tmp_ints(:,:,:,:,iopt+2)
              contr_ints(:,:,:,:,iopt+2) = (coord_bra(1)-gauge_origin(1))*tmp_ints(:,:,:,:,iopt+1) &
                                         - (coord_bra(2)-gauge_origin(2))*tmp_ints(:,:,:,:,iopt)
            end do
          else if (order_mag_bra==0 .and. order_mag_ket==1) then
            do iopt = 1, num_opt-2, 3
              contr_ints(:,:,:,:,iopt) = (coord_ket(3)-gauge_origin(3))*tmp_ints(:,:,:,:,iopt+1) &
                                       - (coord_ket(2)-gauge_origin(2))*tmp_ints(:,:,:,:,iopt+2)
              contr_ints(:,:,:,:,iopt+1) = (coord_ket(1)-gauge_origin(1))*tmp_ints(:,:,:,:,iopt+2) &
                                         - (coord_ket(3)-gauge_origin(3))*tmp_ints(:,:,:,:,iopt)
              contr_ints(:,:,:,:,iopt+2) = (coord_ket(2)-gauge_origin(2))*tmp_ints(:,:,:,:,iopt) &
                                         - (coord_ket(1)-gauge_origin(1))*tmp_ints(:,:,:,:,iopt+1)
            end do
          else
            do iopt = 1, num_opt-2, 3
              contr_ints(:,:,:,:,iopt) = (coord_bra(2)-coord_ket(2))*tmp_ints(:,:,:,:,iopt+2) &
                                       - (coord_bra(3)-coord_ket(3))*tmp_ints(:,:,:,:,iopt+1)
              contr_ints(:,:,:,:,iopt+1) = (coord_bra(3)-coord_ket(3))*tmp_ints(:,:,:,:,iopt) &
                                         - (coord_bra(1)-coord_ket(1))*tmp_ints(:,:,:,:,iopt+2)
              contr_ints(:,:,:,:,iopt+2) = (coord_bra(1)-coord_ket(1))*tmp_ints(:,:,:,:,iopt+1) &
                                         - (coord_bra(2)-coord_ket(2))*tmp_ints(:,:,:,:,iopt)
            end do
          end if
          deallocate(tmp_ints)
        else
          call error_stop("IntGetGAUPOT", "LAO is not implemented", -1)
        end if
      ! non-LAOs
      else
        contr_ints = 0.0
      end if
    else
      ! SGTOs
      if (p_spher_gto) then
        if ((angular_bra==0 .and. angular_ket==0) .or. &
            .not.(present(mag_num_bra) .or. present(mag_num_ket))) then
          call contr_sgto_gaupot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                 exponent_bra, num_contr_bra, contr_coef_bra,   &
                                 idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                 exponent_ket, num_contr_ket, contr_coef_ket,   &
                                 order_elec, idx_gauorg, gaupot_origin,         &
                                 gaupot_expt, idx_diporg, dipole_origin,        &
                                 scal_const, order_mom,                         &
                                 p_order_geo_bra, p_order_geo_ket,              &
                                 p_order_geo_pot, p_order_geo_mom,              &
                                 p_num_cents, p_idx_cent, p_order_cent,         &
                                 num_gto_bra, num_gto_ket, num_opt, contr_ints)
        ! reorders integrals if required
        else
          allocate(tmp_ints(num_gto_bra,num_contr_bra, &
                            num_gto_ket,num_contr_ket,num_opt), stat=ierr)
          if (ierr/=0)                                                     &
            call error_stop("IntGetGAUPOT", "failed to allocate tmp_ints", &
                            num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket*num_opt)
          call contr_sgto_gaupot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                 exponent_bra, num_contr_bra, contr_coef_bra,   &
                                 idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                 exponent_ket, num_contr_ket, contr_coef_ket,   &
                                 order_elec, idx_gauorg, gaupot_origin,         &
                                 gaupot_expt, idx_diporg, dipole_origin,        &
                                 scal_const, order_mom,                         &
                                 p_order_geo_bra, p_order_geo_ket,              &
                                 p_order_geo_pot, p_order_geo_mom,              &
                                 p_num_cents, p_idx_cent, p_order_cent,         &
                                 num_gto_bra, num_gto_ket, num_opt, tmp_ints)
          if (angular_bra>0 .and. present(mag_num_bra) .and. &
              angular_ket>0 .and. present(mag_num_ket)) then
            call reorder_sgto_ints(angular_bra, num_gto_bra, mag_num_bra, &
                                   angular_ket, num_gto_ket, mag_num_ket, &
                                   num_contr_bra, num_contr_ket, num_opt, &
                                   tmp_ints, contr_ints)
          else if (angular_bra>0 .and. present(mag_num_bra)) then
            call reorder_sgtos(angular_bra, num_gto_bra, mag_num_bra, 1, num_contr_bra, &
                               num_gto_ket*num_contr_ket*num_opt, tmp_ints, contr_ints)
          else if (angular_ket>0 .and. present(mag_num_ket)) then
            call reorder_sgtos(angular_ket, num_gto_ket, mag_num_ket,    &
                               num_gto_bra*num_contr_bra, num_contr_ket, &
                               num_opt, tmp_ints, contr_ints)
          end if
          deallocate(tmp_ints)
        end if
      ! CGTOs
      else
        if ((angular_bra==0 .and. angular_ket==0) .or. &
            .not.(present(powers_bra) .or. present(powers_ket))) then
          call contr_cgto_gaupot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                 exponent_bra, num_contr_bra, contr_coef_bra,   &
                                 idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                 exponent_ket, num_contr_ket, contr_coef_ket,   &
                                 order_elec, idx_gauorg, gaupot_origin,         &
                                 gaupot_expt, idx_diporg, dipole_origin,        &
                                 scal_const, order_mom,                         &
                                 p_order_geo_bra, p_order_geo_ket,              &
                                 p_order_geo_pot, p_order_geo_mom,              &
                                 p_num_cents, p_idx_cent, p_order_cent,         &
                                 num_gto_bra, num_gto_ket, num_opt, contr_ints)
        ! reorders integrals if required
        else
          allocate(tmp_ints(num_gto_bra,num_contr_bra, &
                            num_gto_ket,num_contr_ket,num_opt), stat=ierr)
          if (ierr/=0)                                                     &
            call error_stop("IntGetGAUPOT", "failed to allocate tmp_ints", &
                            num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket*num_opt)
          call contr_cgto_gaupot(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                 exponent_bra, num_contr_bra, contr_coef_bra,   &
                                 idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                 exponent_ket, num_contr_ket, contr_coef_ket,   &
                                 order_elec, idx_gauorg, gaupot_origin,         &
                                 gaupot_expt, idx_diporg, dipole_origin,        &
                                 scal_const, order_mom,                         &
                                 p_order_geo_bra, p_order_geo_ket,              &
                                 p_order_geo_pot, p_order_geo_mom,              &
                                 p_num_cents, p_idx_cent, p_order_cent,         &
                                 num_gto_bra, num_gto_ket, num_opt, tmp_ints)
          if (angular_bra>0 .and. present(powers_bra) .and. &
              angular_ket>0 .and. present(powers_ket)) then
            call reorder_cgto_ints(angular_bra, num_gto_bra, powers_bra,  &
                                   angular_ket, num_gto_ket, powers_ket,  &
                                   num_contr_bra, num_contr_ket, num_opt, &
                                   tmp_ints, contr_ints)
          else if (angular_bra>0 .and. present(powers_bra)) then
            call reorder_cgtos(angular_bra, num_gto_bra, powers_bra, 1, num_contr_bra, &
                               num_gto_ket*num_contr_ket*num_opt, tmp_ints, contr_ints)
          else if (angular_ket>0 .and. present(powers_ket)) then
            call reorder_cgtos(angular_ket, num_gto_ket, powers_ket,     &
                               num_gto_bra*num_contr_bra, num_contr_ket, &
                               num_opt, tmp_ints, contr_ints)
          end if
          deallocate(tmp_ints)
        end if
      end if
    end if
    deallocate(p_idx_cent)
    deallocate(p_order_cent)
  end subroutine IntGetGAUPOT

  !> \brief calculates the overlap distribution using contracted Gaussian type orbitals (GTOs)
  !> \author Bin Gao
  !> \date 2012-02-10
  !> \param idx_bra is the atomic index of bra center
  !> \param coord_bra contains the coordinates of bra center
  !> \param angular_bra is the angular number of bra center (s=0, p=1, d=2, ...)
  !> \param num_prim_bra is the number of primitive Gaussians of bra center
  !> \param exponent_bra contains the exponents of primitive Gaussians of bra center
  !> \param num_contr_bra is the number of contractions of bra center
  !> \param contr_coef_bra contains the contraction coefficients of bra center
  !> \param idx_ket is the atomic index of ket center
  !> \param coord_ket contains the coordinates of ket center
  !> \param angular_ket is the angular number of ket center
  !> \param num_prim_ket is the number of primitive Gaussians of ket center
  !> \param exponent_ket contains the exponents of primitive Gaussians of ket center
  !> \param num_contr_ket is the number of contractions of ket center
  !> \param contr_coef_ket contains the contraction coefficients of ket center
  !> \param spher_gto indicates if using spherical GTOs, otherwise Cartesian GTOs
  !> \param info_LAO contains the information of London atomic orbital
  !> \param order_mag_bra is the order of magnetic derivatives on bra center
  !> \param order_mag_ket is the order of magnetic derivatives on ket center
  !> \param order_mag_total is the order of total magnetic derivatives
  !> \param order_ram_bra is the order of derivatives w.r.t. total rotational angular momentum on bra center
  !> \param order_ram_ket is the order of derivatives w.r.t. total rotational angular momentum on ket center
  !> \param order_ram_total is the order of total derivatives w.r.t. total rotational angular momentum
  !> \param order_geo_bra is the order of geometric derivatives with respect to bra center
  !> \param order_geo_ket is the order of geometric derivatives with respect to ket center
  !> \param nary_tree_total contains the information of N-ary tree for total geometric derivatives
  !> \param num_points is the number of grid points
  !> \param grid_points contains the coordinates of grid points
  !> \param num_gto_bra is the number of spherical/Cartesian GTOs on bra center
  !> \param num_gto_ket is the number of spherical/Cartesian GTOs on ket center
  !> \param num_derv is the number of derivatives
  !> \param mag_num_bra contains the magnetic numbers of spherical GTOs on bra center
  !> \param mag_num_ket contains the magnetic numbers of spherical GTOs on ket center
  !> \param powers_bra contains the Cartesian powers of Cartesian GTOs on bra center
  !> \param powers_ket contains the Cartesian powers of Cartesian GTOs on ket center
  !> \return contr_odist contains the calculated contracted overlap distribution
  subroutine IntGetODST(idx_bra, coord_bra, angular_bra, num_prim_bra,  &
                        exponent_bra, num_contr_bra, contr_coef_bra,    &
                        idx_ket, coord_ket, angular_ket, num_prim_ket,  &
                        exponent_ket, num_contr_ket, contr_coef_ket,    &
                        spher_gto, info_LAO,                            &
                        order_mag_bra, order_mag_ket, order_mag_total,  &
                        order_ram_bra, order_ram_ket, order_ram_total,  &
                        order_geo_bra, order_geo_ket, nary_tree_total,  &
                        num_points, grid_points, num_gto_bra,           &
                        num_gto_ket, num_derv, contr_odist,             &
                        mag_num_bra, mag_num_ket, powers_bra, powers_ket)
    integer, intent(in) :: idx_bra
    real(REALK), intent(in) :: coord_bra(3)
    integer, intent(in) :: angular_bra
    integer, intent(in) :: num_prim_bra
    real(REALK), intent(in) :: exponent_bra(num_prim_bra)
    integer, intent(in) :: num_contr_bra
    real(REALK), intent(in) :: contr_coef_bra(num_contr_bra,num_prim_bra)
    integer, intent(in) :: idx_ket
    real(REALK), intent(in) :: coord_ket(3)
    integer, intent(in) :: angular_ket
    integer, intent(in) :: num_prim_ket
    real(REALK), intent(in) :: exponent_ket(num_prim_ket)
    integer, intent(in) :: num_contr_ket
    real(REALK), intent(in) :: contr_coef_ket(num_contr_ket,num_prim_ket)
    logical, optional, intent(in) :: spher_gto
    type(london_ao_t), intent(in) :: info_LAO
    integer, intent(in) :: order_mag_bra
    integer, intent(in) :: order_mag_ket
    integer, intent(in) :: order_mag_total
    integer, intent(in) :: order_ram_bra
    integer, intent(in) :: order_ram_ket
    integer, intent(in) :: order_ram_total
    integer, optional, intent(in) :: order_geo_bra
    integer, optional, intent(in) :: order_geo_ket
    type(nary_tree_t), optional, intent(in) :: nary_tree_total
    integer, intent(in) :: num_points
    real(REALK), intent(in) :: grid_points(3,num_points)
    integer, intent(in) :: num_gto_bra
    integer, intent(in) :: num_gto_ket
    integer, intent(in) :: num_derv
    real(REALK), intent(out) :: contr_odist(num_gto_bra,num_contr_bra, &
                                            num_gto_ket,num_contr_ket, &
                                            num_points,num_derv)
    integer, optional, intent(in) :: mag_num_bra(num_gto_bra)
    integer, optional, intent(in) :: mag_num_ket(num_gto_ket)
    integer, optional, intent(in) :: powers_bra(3,num_gto_bra)
    integer, optional, intent(in) :: powers_ket(3,num_gto_ket)
    logical p_spher_gto                                 !arguments for Gen1Int (local)
    integer p_order_geo_bra
    integer p_order_geo_ket
    integer p_num_cents
    integer, allocatable :: p_idx_cent(:)
    integer, allocatable :: p_order_cent(:)
    integer icent
    real(REALK), allocatable :: tmp_odist(:,:,:,:,:,:)  !contracted overlap distribution from Gen1Int
    integer ierr                                        !error information
    ! sets the arguments for Gen1Int (local)
    if (present(spher_gto)) then
      p_spher_gto = spher_gto
    else
      p_spher_gto = .true.
    end if
    if (present(order_geo_bra)) then
      p_order_geo_bra = order_geo_bra
    else
      p_order_geo_bra = 0
    end if
    if (present(order_geo_ket)) then
      p_order_geo_ket = order_geo_ket
    else
      p_order_geo_ket = 0
    end if
    if (present(nary_tree_total)) then
      p_num_cents = nary_tree_total%wt_node(nary_tree_total%order_geo)
    else
      p_num_cents = 0
    end if
    if (p_num_cents>0) then
      allocate(p_idx_cent(p_num_cents), stat=ierr)
      if (ierr/=0) &
        call error_stop("IntGetODST", "failed to allocate p_idx_cent", p_num_cents)
      allocate(p_order_cent(p_num_cents), stat=ierr)
      if (ierr/=0) &
        call error_stop("IntGetODST", "failed to allocate p_order_cent", p_num_cents)
      if (allocated(nary_tree_total%idx_atoms)) then
        do icent = 1, p_num_cents
          p_idx_cent(icent) = nary_tree_total%idx_atoms(nary_tree_total%idx_cent(icent))
        end do
      else
        p_idx_cent = nary_tree_total%idx_cent(1:p_num_cents)
      end if
      p_order_cent = nary_tree_total%order_cent(1:p_num_cents)
    else
      allocate(p_idx_cent(1), stat=ierr)
      if (ierr/=0) call error_stop("IntGetODST", "failed to allocate p_idx_cent", 1)
      allocate(p_order_cent(1), stat=ierr)
      if (ierr/=0) call error_stop("IntGetODST", "failed to allocate p_order_cent", 1)
      p_idx_cent = 0
      p_order_cent = 0
    end if
    ! magnetic derivatives or derivatives with respect to rotational angular momentum
    if (order_mag_total>0 .or. order_mag_bra>0 .or. order_mag_ket>0 .or. &
        order_ram_total>0 .or. order_ram_bra>0 .or. order_ram_ket>0) then
      ! London atomic orbitals
      if (LondonAOUsed(info_LAO)) then
        call error_stop("IntGetODST", "LAO is not implemented", -1)
      ! non-LAOs
      else
        contr_odist = 0.0
      end if
    else
      ! SGTOs
      if (p_spher_gto) then
        if ((angular_bra==0 .and. angular_ket==0) .or. &
            .not.(present(mag_num_bra) .or. present(mag_num_ket))) then
          call contr_sgto_odist(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                exponent_bra, num_contr_bra, contr_coef_bra,   &
                                idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                exponent_ket, num_contr_ket, contr_coef_ket,   &
                                p_order_geo_bra, p_order_geo_ket,              &
                                p_num_cents, p_idx_cent, p_order_cent,         &
                                num_points, grid_points,                       &
                                num_gto_bra, num_gto_ket, num_derv, contr_odist)
        ! reorders integrals if required
        else
          allocate(tmp_odist(num_gto_bra,num_contr_bra, &
                             num_gto_ket,num_contr_ket, &
                             num_points,num_derv), stat=ierr)
          if (ierr/=0)                                                    &
            call error_stop("IntGetODST", "failed to allocate tmp_odist", &
                            num_gto_bra*num_contr_bra*num_gto_ket         &
                            *num_contr_ket*num_points*num_derv)
          call contr_sgto_odist(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                exponent_bra, num_contr_bra, contr_coef_bra,   &
                                idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                exponent_ket, num_contr_ket, contr_coef_ket,   &
                                p_order_geo_bra, p_order_geo_ket,              &
                                p_num_cents, p_idx_cent, p_order_cent,         &
                                num_points, grid_points,                       &
                                num_gto_bra, num_gto_ket, num_derv, tmp_odist)
          if (angular_bra>0 .and. present(mag_num_bra) .and. &
              angular_ket>0 .and. present(mag_num_ket)) then
            call reorder_sgto_ints(angular_bra, num_gto_bra, mag_num_bra, &
                                   angular_ket, num_gto_ket, mag_num_ket, &
                                   num_contr_bra, num_contr_ket,          &
                                   num_points*num_derv, tmp_odist, contr_odist)
          else if (angular_bra>0 .and. present(mag_num_bra)) then
            call reorder_sgtos(angular_bra, num_gto_bra, mag_num_bra, 1, num_contr_bra, &
                               num_gto_ket*num_contr_ket*num_points*num_derv,           &
                               tmp_odist, contr_odist)
          else if (angular_ket>0 .and. present(mag_num_ket)) then
            call reorder_sgtos(angular_ket, num_gto_ket, mag_num_ket,    &
                               num_gto_bra*num_contr_bra, num_contr_ket, &
                               num_points*num_derv, tmp_odist, contr_odist)
          end if
          deallocate(tmp_odist)
        end if
      ! CGTOs
      else
        if ((angular_bra==0 .and. angular_ket==0) .or. &
            .not.(present(powers_bra) .or. present(powers_ket))) then
          call contr_cgto_odist(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                exponent_bra, num_contr_bra, contr_coef_bra,   &
                                idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                exponent_ket, num_contr_ket, contr_coef_ket,   &
                                p_order_geo_bra, p_order_geo_ket,              &
                                p_num_cents, p_idx_cent, p_order_cent,         &
                                num_points, grid_points,                       &
                                num_gto_bra, num_gto_ket, num_derv, contr_odist)
        ! reorders integrals if required
        else
          allocate(tmp_odist(num_gto_bra,num_contr_bra, &
                             num_gto_ket,num_contr_ket, &
                             num_points,num_derv), stat=ierr)
          if (ierr/=0)                                                    &
            call error_stop("IntGetODST", "failed to allocate tmp_odist", &
                            num_gto_bra*num_contr_bra*num_gto_ket         &
                            *num_contr_ket*num_points*num_derv)
          call contr_cgto_odist(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                exponent_bra, num_contr_bra, contr_coef_bra,   &
                                idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                exponent_ket, num_contr_ket, contr_coef_ket,   &
                                p_order_geo_bra, p_order_geo_ket,              &
                                p_num_cents, p_idx_cent, p_order_cent,         &
                                num_points, grid_points,                       &
                                num_gto_bra, num_gto_ket, num_derv, tmp_odist)
          if (angular_bra>0 .and. present(powers_bra) .and. &
              angular_ket>0 .and. present(powers_ket)) then
            call reorder_cgto_ints(angular_bra, num_gto_bra, powers_bra,  &
                                   angular_ket, num_gto_ket, powers_ket,  &
                                   num_contr_bra, num_contr_ket,          &
                                   num_points*num_derv, tmp_odist, contr_odist)
          else if (angular_bra>0 .and. present(powers_bra)) then
            call reorder_cgtos(angular_bra, num_gto_bra, powers_bra, 1, num_contr_bra,   &
                               num_gto_ket*num_contr_ket*num_points*num_derv, tmp_odist, &
                               contr_odist)
          else if (angular_ket>0 .and. present(powers_ket)) then
            call reorder_cgtos(angular_ket, num_gto_ket, powers_ket,     &
                               num_gto_bra*num_contr_bra, num_contr_ket, &
                               num_points*num_derv, tmp_odist, contr_odist)
          end if
          deallocate(tmp_odist)
        end if
      end if
    end if
    deallocate(p_idx_cent)
    deallocate(p_order_cent)
  end subroutine IntGetODST

end module gen1int_geom
