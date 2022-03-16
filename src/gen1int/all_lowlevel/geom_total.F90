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
!!  This file contains subroutines related to total geometric derivatives.
!!
!!  2011-12-13, Bin Gao:
!!  * adds subroutine \fn(geom_total_num_redunt) to calculate the number of
!!    redundant geometric derivatives for a given path
!!  * adds subroutine \fn(geom_total_base_unique) to return the list addresses of
!!    unique total geometric derivatives for a given path
!!  * adds subroutine \fn(geom_total_redunt_list) to return the list addresses of
!!    redundant total geometric derivatives for a given path
!!
!!  2011-11-06, Bin Gao:
!!  * adds subroutine \fn(geom_total_redunt_expectation) to return the expectation
!!    values of redundant total geometric derivatives
!!
!!  2011-10-05, Bin Gao:
!!  * adds subroutine \fn(geom_total_num_derv) to calculate the total number
!!    of total geometric derivatives
!!
!!  2011-07-03, Bin Gao:
!!  * calling \fn(geom_total_num_paths) from \fn(geom_total_tree_init)
!!
!!  2010-02-06, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief computes the number of unique total geometric derivatives with the given
  !>        (1) order of geometric derivatives, (2) maximum number of differentiated centers,
  !>        and (3) number of atoms to select
  !> \author Bin Gao
  !> \date 2009-08-27
  !> \param order_geo is the order of total geometric derivatives
  !> \param max_num_cent is the maximum number of differentiated centers
  !> \param num_atoms is the number of atoms
  !> \return num_derv is the total number of total geometric derivatives
  subroutine geom_total_num_derv(order_geo, max_num_cent, num_atoms, num_derv)
    use xkind
    implicit none
    integer, intent(in) :: order_geo
    integer, intent(in) :: max_num_cent
    integer, intent(in) :: num_atoms
    integer, intent(out) :: num_derv
!f2py intent(in) :: order_geo
!f2py intent(in) :: max_num_cent
!f2py intent(in) :: num_atoms
!f2py intent(out) :: num_derv
    real(REALK) num_select                 !number of selections of atoms
    real(REALK) int_result                 !intermediate result
    real(REALK) real_order                 !order of geometric derivatives, real number
    real(REALK) real_num_atoms             !number of atoms, real number
    real(REALK) real_num_derv              !total number of geometric derivatives, real number
    integer num_paths                      !total number of different paths
    integer visit_height                   !height of atom to visit
    integer, allocatable :: idx_node(:)    !indices of the selected atom nodes
    integer, allocatable :: wt_node(:)     !weights of the selected atom nodes
    integer, allocatable :: idx_cent(:)    !indices of generated differentiated centers
    integer, allocatable :: order_cent(:)  !orders of derivatives of the differentiated centers
    integer num_unique_geo                 !number of unique total geometric derivatives for the generated path
    integer ipath                          !incremental recorder over the paths
    integer ierr                           !error information
#if defined(XTIME)
    real(REALK) curr_time       !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    select case(max_num_cent)
    ! one-center geometric derivatives
    case(1)
      num_derv = num_atoms*(order_geo+1)*(order_geo+2)/2
    ! one- and two-center geometric derivatives
    case(2)
      real_order = real(order_geo,REALK)
      real_num_atoms = real(num_atoms,REALK)
      int_result = real((order_geo+1)*(order_geo+2)/2,REALK)
      ! one-center
      real_num_derv = real_num_atoms*int_result
      ! two-center
      num_select = real((num_atoms-1)*num_atoms/2,REALK)
      num_derv = nint(real_num_derv+num_select*int_result*(real_order-1.0_REALK) &
                      *((real_order+13.0_REALK)*real_order+60.0_REALK)/60.0_REALK)
    ! one-, two- and three-center geometric derivatives
    case(3)
      real_order = real(order_geo,REALK)
      real_num_atoms = real(num_atoms,REALK)
      int_result = real((order_geo+1)*(order_geo+2)/2,REALK)
      ! one-center
      real_num_derv = real_num_atoms*int_result
      ! two-center
      num_select = real((num_atoms-1)*num_atoms/2,REALK)
      int_result = int_result*(real_order-1.0_REALK)/60.0_REALK
      real_num_derv = real_num_derv+num_select*int_result &
                    * ((real_order+13.0_REALK)*real_order+60.0_REALK)
      ! three-center
      num_select = real(((num_atoms-3)*num_atoms+2)*num_atoms/6,REALK)
      num_derv = nint(real_num_derv+num_select*int_result*(real_order-2)  &
                      *((((real_order+36.0_REALK)*real_order+551.0_REALK) &
                      *real_order+3708.0_REALK)*real_order+10080.0_REALK)/336.0_REALK)
    ! one-, two-, three- and four-center geometric derivatives
    case(4)
      real_order = real(order_geo,REALK)
      real_num_atoms = real(num_atoms,REALK)
      int_result = real((order_geo+1)*(order_geo+2)/2,REALK)
      ! one-center
      real_num_derv = real_num_atoms*int_result
      ! two-center
      num_select = real((num_atoms-1)*num_atoms/2,REALK)
      int_result = int_result*(real_order-1.0_REALK)/60.0_REALK
      real_num_derv = real_num_derv+num_select*int_result &
                    * ((real_order+13.0_REALK)*real_order+60.0_REALK)
      ! three-center
      num_select = real(((num_atoms-3)*num_atoms+2)*num_atoms/6,REALK)
      int_result = int_result*(real_order-2.0_REALK)/336.0_REALK
      real_num_derv = real_num_derv+num_select*int_result                &
                    * ((((real_order+36.0_REALK)*real_order+551.0_REALK) &
                    * real_order+3708.0_REALK)*real_order+10080.0_REALK)
      ! four-center
      num_select = real((((num_atoms-6)*num_atoms+11)*num_atoms-6)*num_atoms/24,REALK)
      num_derv = nint(real_num_derv+num_select*int_result*(real_order-3.0_REALK) &
                      *((((((real_order+69.0_REALK)*real_order+2137.0_REALK)     &
                      *real_order+35451.0_REALK)*real_order+330862.0_REALK)      &
                      *real_order+1612920.0_REALK)*real_order+3326400.0_REALK)/990.0_REALK)
    ! not implemented
    case default
      ! allocates memory for total geometric derivatives
      allocate(idx_node(order_geo), stat=ierr)
      if (ierr/=0)                                     &
        call error_stop("geom_total_num_derv",         &
                        "failed to allocate idx_node", &
                        order_geo)
      allocate(wt_node(order_geo), stat=ierr)
      if (ierr/=0)                                    &
        call error_stop("geom_total_num_derv",        &
                        "failed to allocate wt_node", &
                        order_geo)
      allocate(idx_cent(max_num_cent), stat=ierr)
      if (ierr/=0)                                     &
        call error_stop("geom_total_num_derv",         &
                        "failed to allocate idx_cent", &
                        max_num_cent)
      allocate(order_cent(max_num_cent), stat=ierr)
      if (ierr/=0)                                       &
        call error_stop("geom_total_num_derv",           &
                        "failed to allocate order_cent", &
                        max_num_cent)
      ! generates the first path
      call geom_total_tree_init(num_atoms, order_geo, max_num_cent, &
                                num_paths, visit_height, idx_node,  &
                                wt_node, idx_cent, order_cent, num_derv)
      ! searches for the next satisfied path from a given path
      do ipath = 2, num_paths
        call geom_total_tree_search(num_atoms, order_geo, max_num_cent, &
                                    visit_height, idx_node, wt_node,    &
                                    idx_cent, order_cent, num_unique_geo)
        num_derv = num_derv+num_unique_geo
      end do
      ! cleans
      deallocate(idx_node)
      deallocate(wt_node)
      deallocate(idx_cent)
      deallocate(order_cent)
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "geom_total_num_derv", STDOUT)
#endif
    return
  end subroutine geom_total_num_derv

  !> \brief returns the total number of different paths, and generates the first path
  !> \author Bin Gao
  !> \date 2011-04-04
  !> \param num_atoms is the number of atoms
  !> \param order_geo is the order of total geometric derivatives
  !> \param max_num_cent is the maximum number of differentiated centers
  !> \return num_paths is the total number of different paths
  !> \return visit_height is the height of atom to visit
  !> \return idx_node contains the indices of the selected atom nodes
  !> \return wt_node contains the weights of the selected atom nodes
  !> \return idx_cent contains the indices of generated differentiated centers
  !> \return order_cent contains the orders of derivatives of the differentiated centers
  !> \return num_unique_geo is the number of unique total geometric derivatives for the generated path
  subroutine geom_total_tree_init(num_atoms, order_geo, max_num_cent, &
                                  num_paths, visit_height, idx_node,  &
                                  wt_node, idx_cent, order_cent, num_unique_geo)
    use xkind
    implicit none
    integer, intent(in) :: num_atoms
    integer, intent(in) :: order_geo
    integer, intent(in) :: max_num_cent
    integer, intent(out) :: num_paths
    integer, intent(out) :: visit_height
    integer, intent(out) :: idx_node(order_geo)
    integer, intent(out) :: wt_node(order_geo)
    integer, intent(out) :: idx_cent(max_num_cent)
    integer, intent(out) :: order_cent(max_num_cent)
    integer, intent(out) :: num_unique_geo
!f2py intent(in) :: num_atoms
!f2py intent(in) :: order_geo
!f2py intent(in) :: max_num_cent
!f2py intent(out) :: num_paths
!f2py intent(out) :: visit_height
!f2py intent(out) :: idx_node
!f2py depend(order_geo) :: idx_node
!f2py intent(out) :: wt_node
!f2py depend(order_geo) :: wt_node
!f2py intent(out) :: idx_cent
!f2py depend(max_num_cent) :: idx_cent
!f2py intent(out) :: order_cent
!f2py depend(max_num_cent) :: order_cent
!f2py intent(out) :: num_unique_geo
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! computes the number of different paths
    call geom_total_num_paths(num_atoms, order_geo, max_num_cent, num_paths)
    ! the first path
    visit_height = order_geo
    idx_node = 1
    wt_node = 1
    idx_cent(1) = 1
    idx_cent(2:max_num_cent) = 0
    order_cent(1) = order_geo
    order_cent(2:max_num_cent) = 0
    ! computes the number of unique total geometric derivatives of the first path
    num_unique_geo = (order_geo+1)*(order_geo+2)/2
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "geom_total_tree_init", STDOUT)
#endif
    return
  end subroutine geom_total_tree_init

  !> \brief computes the total number of different paths, more explicitly, it equals to
  !>        \f$\sum_{k=1}^{\min(N_{cent},N_{order})}\binom{N_{atom}}{k}\binom{N_{order}-1}{k-1}\f$
  !> \author Bin Gao
  !> \date 2011-03-06
  !> \param num_atoms is the number of atoms
  !> \param order_geo is the order of total geometric derivatives
  !> \param max_num_cent is the maximum number of centers
  !> \return num_paths is the total number of different paths
  subroutine geom_total_num_paths(num_atoms, order_geo, max_num_cent, num_paths)
    use xkind
    implicit none
    integer, intent(in) :: num_atoms
    integer, intent(in) :: order_geo
    integer, intent(in) :: max_num_cent
    integer, intent(out) :: num_paths
!f2py intent(in) :: num_atoms
!f2py intent(in) :: order_geo
!f2py intent(in) :: max_num_cent
!f2py intent(out) :: num_paths
    integer icent          !incremental recorder of number of centers
    integer order_geo1     !\var(order_geo)-1
    integer nnz_comps      !number of different non-zero compositions of the order
    integer num_select     !number of selections of different atoms to be the differentiated centers
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    order_geo1 = order_geo-1
    num_paths = 0
    ! considers number of centers from 1 to \fn(min)(\var(max_num_cent),\var(order_geo))
    do icent = 1, min(max_num_cent, order_geo)
      ! computes the number of different non-zero compositions of the order
      ! puts \var(icent)-1 |'s into the left \var(order_geo)-1 spaces
      ! in |X X ... X X X| (\var(order_geo) X's)
      call dbinom_coeff(order_geo1, icent-1, nnz_comps)
      ! computes the number of selections of different atoms to be the differentiated centers
      call dbinom_coeff(num_atoms, icent, num_select)
      ! updates the number of paths
      num_paths = num_paths+num_select*nnz_comps
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "geom_total_num_paths", STDOUT)
#endif
    return
  end subroutine geom_total_num_paths

  !> \brief searches for the next satisfied path from a given path, could be called recursively
  !> \author Bin Gao
  !> \date 2010-02-06
  !> \param num_atoms is the number of atoms
  !> \param order_geo is the order of total geometric derivatives
  !> \param max_num_cent is the maximum number of differentiated centers
  !> \return visit_height is the height of atom to visit
  !> \return idx_node contains the indices of the selected atom nodes
  !> \return wt_node contains the weights of the selected atom nodes
  !> \return idx_cent contains the indices of the generated differentiated centers
  !> \return order_cent contains the orders of derivatives of the differentiated centers
  !> \return num_unique_geo is the number of unique total geometric derivatives for the generated path
  subroutine geom_total_tree_search(num_atoms, order_geo, max_num_cent, &
                                    visit_height, idx_node, wt_node,    &
                                    idx_cent, order_cent, num_unique_geo)
    use xkind
    implicit none
    integer, intent(in) :: num_atoms
    integer, intent(in) :: order_geo
    integer, intent(in) :: max_num_cent
    integer, intent(inout) :: visit_height
    integer, intent(inout) :: idx_node(order_geo)
    integer, intent(inout) :: wt_node(order_geo)
    integer, intent(inout) :: idx_cent(max_num_cent)
    integer, intent(inout) :: order_cent(max_num_cent)
    integer, intent(out) :: num_unique_geo
!f2py intent(in) :: num_atoms
!f2py intent(in) :: order_geo
!f2py intent(in) :: max_num_cent
!f2py intent(inout) :: visit_height
!f2py intent(inout) :: idx_node
!f2py depend(order_geo) :: idx_node
!f2py intent(inout) :: wt_node
!f2py depend(order_geo) :: wt_node
!f2py intent(inout) :: idx_cent
!f2py depend(max_num_cent) :: idx_cent
!f2py intent(inout) :: order_cent
!f2py depend(max_num_cent) :: order_cent
!f2py intent(out) :: num_unique_geo
    !-integer num_cents   !number of differentiated centers in the given path
    !-integer num_diff    !difference between the number of differentiated centers in the given and new paths
    integer inode          !incremental recorder
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    !-! saves the number of differentiated centers in the given path
    !-num_cents = wt_node(order_geo)
    ! generates new path of differentiated centers from the given path
    call geom_total_new_path(num_atoms, order_geo, visit_height, idx_node, wt_node)
    ! if the number of atomic centers exceeds, we then need to visit its parent node
    do while (wt_node(visit_height)>max_num_cent)
      visit_height = visit_height-1
      ! generates new path of differentiated centers from current path
      call geom_total_new_path(num_atoms, order_geo, visit_height, idx_node, wt_node)
    end do
    ! if previously idx_node(visit_height)==idx_node(visit_height-1), i.e.
    ! this node was the same as its parent node ago, they are different now
    ! after replaced by the sibling, so that we need to decrease the order of
    ! derivatives of the parent node
    if (visit_height>1) then
      if (idx_node(visit_height-1)+1==idx_node(visit_height)) then
        !-inode = num_cents-wt_node(visit_height-1)+1
        inode = wt_node(visit_height-1)
        order_cent(inode) = order_cent(inode)-1
      end if
    end if
    ! the index of the atomic center,
    idx_cent(wt_node(visit_height)) = idx_node(visit_height)
    ! and its order of geometric derivatives
    order_cent(wt_node(visit_height)) = order_geo-visit_height+1
    !-centers are in descending order
    !-! we move the differentiated centers idx_cent(num_cents-wt_node(visit_height)+2) ...
    !-! idx_cent(num_cents) to idx_cent(2) ... idx_cent(wt_node(visit_height))
    !-!
    !-! just replacing the first differentiated atomic center
    !-if (num_cents==wt_node(visit_height)) then
    !-  ! the index of the first differentiated atomic center,
    !-  ! and its order of geometric derivatives
    !-  idx_cent(1) = idx_node(visit_height)
    !-  order_cent(1) = order_geo-visit_height+1
    !-  ! initializes the number of unique total geometric derivatives
    !-  num_unique_geo = (order_cent(1)+1)*(order_cent(1)+2)/2
    !-  ! upates the number of unique total geometric derivatives
    !-  do inode = 2, num_cents
    !-    num_unique_geo = num_unique_geo*(order_cent(inode)+1)*(order_cent(inode)+2)/2
    !-  end do
    !-! moving centers forward
    !-else if (num_cents>wt_node(visit_height)) then
    !-  ! the index of the first differentiated atomic center,
    !-  ! and its order of geometric derivatives
    !-  idx_cent(1) = idx_node(visit_height)
    !-  order_cent(1) = order_geo-visit_height+1
    !-  ! initializes the number of unique total geometric derivatives
    !-  num_unique_geo = (order_cent(1)+1)*(order_cent(1)+2)/2
    !-  ! moving the centers and upates the number of unique total geometric derivatives
    !-  num_diff = num_cents-wt_node(visit_height)
    !-  do inode = 2, wt_node(visit_height)
    !-    idx_cent(inode) = idx_cent(inode+num_diff)
    !-    order_cent(inode) = order_cent(inode+num_diff)
    !-    num_unique_geo = num_unique_geo*(order_cent(inode)+1)*(order_cent(inode)+2)/2
    !-  end do
    !-! moving centers afterward
    !-else if (num_cents<wt_node(visit_height)) then
    !-  ! initializes the number of unique total geometric derivatives
    !-  num_unique_geo = 1
    !-  ! moving the centers and upates the number of unique total geometric derivatives
    !-  num_diff = num_cents-wt_node(visit_height)
    !-  do inode = wt_node(visit_height), 2, -1
    !-    idx_cent(inode) = idx_cent(inode+num_diff)
    !-    order_cent(inode) = order_cent(inode+num_diff)
    !-    num_unique_geo = num_unique_geo*(order_cent(inode)+1)*(order_cent(inode)+2)/2
    !-  end do
    !-  ! the index of the first differentiated atomic center,
    !-  ! and its order of geometric derivatives
    !-  idx_cent(1) = idx_node(visit_height)
    !-  order_cent(1) = order_geo-visit_height+1
    !-  ! upates the number of unique total geometric derivatives
    !-  num_unique_geo = num_unique_geo*(order_cent(1)+1)*(order_cent(1)+2)/2
    !-end if
    ! sets the child nodes as the same as current visiting node
    do inode = visit_height+1, order_geo
      idx_node(inode) = idx_node(visit_height)
      wt_node(inode) = wt_node(visit_height)
    end do
    ! gets the number of unique total geometric derivatives
    num_unique_geo = (order_cent(1)+1)*(order_cent(1)+2)/2
    do inode = 2, wt_node(order_geo)
      num_unique_geo = num_unique_geo*(order_cent(inode)+1)*(order_cent(inode)+2)/2
    end do
    ! we will visit the leaf node next time
    if (idx_node(visit_height)<num_atoms) then
      visit_height = order_geo
    ! we will visit the parent node next time
    else
      visit_height = visit_height-1
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "geom_total_tree_search", STDOUT)
#endif
    return
  end subroutine geom_total_tree_search

  !> \brief generates a new path of differentiated centers from a given path,
  !>        the first path will return if the last path is given
  !> \author Bin Gao
  !> \date 2010-02-06
  !> \param num_atoms is the number of atoms
  !> \param order_geo is the order of total geometric derivatives
  !> \return visit_height is the height of atom to visit
  !> \return idx_node contains the indices of the selected atom nodes
  !> \return wt_node contains the weights of the selected atom nodes
  subroutine geom_total_new_path(num_atoms, order_geo, visit_height, idx_node, wt_node)
    use xkind
    implicit none
    integer, intent(in) :: num_atoms
    integer, intent(in) :: order_geo
    integer, intent(inout) :: visit_height
    integer, intent(inout) :: idx_node(order_geo)
    integer, intent(inout) :: wt_node(order_geo)
!f2py intent(in) :: num_atoms
!f2py intent(in) :: order_geo
!f2py intent(in) :: max_num_cent
!f2py intent(inout) :: visit_height
!f2py intent(inout) :: idx_node
!f2py depend(order_geo) :: idx_node
!f2py intent(inout) :: wt_node
!f2py depend(order_geo) :: wt_node
    logical do_replace  !indicates that we found the highest node to replace with its sibling
    integer inode       !incremental recorder
    ! if the visiting node is the last atom, we then need to visit its parent node
    do_replace = .false.
    do inode = visit_height, 1, -1
      ! we found the highest node to replace with its sibling
      if (idx_node(inode)<num_atoms) then
        visit_height = inode
        do_replace = .true.
        exit
      end if
    end do
    if (do_replace) then
      ! replace this node with its sibling
      idx_node(visit_height) = idx_node(visit_height)+1
      ! this node (replaced by its sibling) is then different from its parent node,
      ! since previously idx_node(visit_height)>=idx_node(visit_height-1), it becomes
      ! definitely idx_node(visit_height)>idx_node(visit_height-1), so that we need to
      ! update the weight of this node
      if (visit_height>1) wt_node(visit_height) = wt_node(visit_height-1)+1
    ! we arrive at the final path since last call, go back to the first path
    else
      visit_height = 1
      idx_node(1) = 1
      wt_node(1) = 1
    end if
    return
  end subroutine geom_total_new_path

  !> \brief returns the base addresses of unique total geometric derivatives for a given path
  !> \author Bin Gao
  !> \date 2011-12-13
  !> \param num_cent is the number of differentiated centers
  !> \param idx_cent contains the indices of differentiated centers
  !> \param order_cent contains the orders of differentiated centers
  !> \param order_geo is the order of total geometric derivatives
  !> \param num_unique_geo is the number of unique total geometric derivatives
  !> \return base_unique contains the base addresses of unique total geometric derivatives
  subroutine geom_total_base_unique(num_cent, idx_cent, order_cent, &
                                    order_geo, num_unique_geo, base_unique)
    use xkind
    implicit none
    integer, intent(in) :: num_cent
    integer, intent(in) :: idx_cent(num_cent)
    integer, intent(in) :: order_cent(num_cent)
    integer, intent(in) :: order_geo
    integer, intent(in) :: num_unique_geo
    integer, intent(out) :: base_unique(order_geo,num_unique_geo)
!f2py intent(hide) :: num_cent
!f2py intent(in) :: idx_cent
!f2py intent(in) :: order_cent
!f2py depend(num_cent) :: order_cent
!f2py intent(in) :: order_geo
!f2py intent(in) :: num_unique_geo
!f2py intent(out) :: base_unique
!f2py depend(order_geo) :: base_unique
!f2py depend(num_unique_geo) :: base_unique
    integer pnum_unique_geo         !number of unique total geometric derivatives (private)
    integer last_prod_nderv         !product of number of derivatives up to last differentiated center
    integer curr_prod_nderv         !product of number of derivatives up to current differentiated center
    integer base_xyz_derv           !base address of derivatives along xyz directions in \var(base_unique)
    integer base_cent_derv(3)       !base addresses of the first order derivatives for a given center
    integer num_cent_derv           !number of derivatives for a given center
    integer strt_xderv, end_xderv   !start and end addresses of derivatives in \var(base_unique)
    integer strt_yderv, end_yderv
    integer strt_zderv, end_zderv
    integer icent                   !incremental recorders
    integer iderv, jderv
    integer x_derv, y_derv, z_derv
    integer ixyz
    integer p_unique_geo            !pointer for the unique total geometric derivatives
#if defined(XTIME)
    real(REALK) curr_time           !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! checks the order of total geometric derivatives
    if (order_geo/=sum(order_cent))                                     &
      call error_stop("geom_total_base_unique",                         &
                      "incorrect order of total geometric derivatives", &
                      order_geo)
    ! checks the number of unique total geometric derivatives
    pnum_unique_geo = (order_cent(1)+1)*(order_cent(1)+2)/2
    do icent = 2, num_cent
      pnum_unique_geo = pnum_unique_geo*(order_cent(icent)+1)*(order_cent(icent)+2)/2
    end do
    if (pnum_unique_geo/=num_unique_geo)                                        &
      call error_stop("geom_total_base_unique",                                 &
                      "incorrect number of unique total geometric derivatives", &
                      pnum_unique_geo)
    ! initializes quantities
    last_prod_nderv = 1
    base_xyz_derv = 0
    ! loops over differentiated centers \var(idx_cent) in ascending order
    do icent = 1, num_cent
      base_cent_derv(1) = 3*idx_cent(icent)-3  !px
      base_cent_derv(2) = base_cent_derv(1)+1  !py
      base_cent_derv(3) = base_cent_derv(2)+1  !pz
      num_cent_derv = (order_cent(icent)+1)*(order_cent(icent)+2)/2
      curr_prod_nderv = num_cent_derv*last_prod_nderv
      strt_xderv = base_xyz_derv+1
      end_zderv = base_xyz_derv+order_cent(icent)
      ! generates the base addresses of unique total geometric derivatives
      p_unique_geo = 0
      do iderv = 1, num_unique_geo/curr_prod_nderv
        ! xyz components for the current differentitated center
        do z_derv = 0, order_cent(icent)
          do y_derv = 0, order_cent(icent)-z_derv
            x_derv = order_cent(icent)-(y_derv+z_derv)
            end_xderv = base_xyz_derv+x_derv
            strt_yderv = end_xderv+1
            end_yderv = end_xderv+y_derv
            strt_zderv = end_yderv+1
            ! geometric derivatives of preceding differentitated centers
            do jderv = 1, last_prod_nderv
              p_unique_geo = p_unique_geo+1
              do ixyz = strt_xderv, end_xderv
                base_unique(ixyz,p_unique_geo) = base_cent_derv(1)
              end do
              do ixyz = strt_yderv, end_yderv
                base_unique(ixyz,p_unique_geo) = base_cent_derv(2)
              end do
              do ixyz = strt_zderv, end_zderv
                base_unique(ixyz,p_unique_geo) = base_cent_derv(3)
              end do
            end do
          end do
        end do
      end do
      last_prod_nderv = curr_prod_nderv
      base_xyz_derv = base_xyz_derv+order_cent(icent)
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "geom_total_base_unique", STDOUT)
#endif
    return
  end subroutine geom_total_base_unique

  !> \brief returns the number of redundant total geometric derivatives
  !> \author Bin Gao
  !> \date 2011-12-13
  !> \param num_cent is the number of differentiated centers
  !> \param order_cent contains the orders of differentiated centers
  !> \return num_redunt_geo is the number of redundant total geometric derivatives
  subroutine geom_total_num_redunt(num_cent, order_cent, num_redunt_geo)
    use xkind
    implicit none
    integer, intent(in) :: num_cent
    integer, intent(in) :: order_cent(num_cent)
    integer, intent(out) :: num_redunt_geo
!f2py intent(hide) :: num_cent
!f2py intent(in) :: order_cent
!f2py intent(out) :: num_redunt_geo
    integer order_geo      !current order of total geometric derivatives
    integer icent          !incremental recorder over centers
    integer num_subset     !number of subsets for current differentiated order
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! sets the order of total geometric derivatives
    order_geo = sum(order_cent)
!FIXME: checks if the numbers exceed bound
    num_redunt_geo = 3**order_geo
    do icent = 1, num_cent
      call dbinom_coeff(order_geo, order_cent(icent), num_subset)
      num_redunt_geo = num_subset*num_redunt_geo
      order_geo = order_geo-order_cent(icent)
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "geom_total_num_redunt", STDOUT)
#endif
    return
  end subroutine geom_total_num_redunt

  !> \brief returns the list addresses of redundant total geometric derivatives for a given path
  !> \author Bin Gao
  !> \date 2011-12-13
  !> \param num_atoms is the number of atoms
  !> \param num_cent is the number of differentiated centers
  !> \param idx_cent contains the indices of differentiated centers
  !> \param order_cent contains the orders of differentiated centers
  !> \param sort_redunt indicates if sorting the list addresses according to the redundant
  !>        total geometric derivatives
  !> \param num_redunt_geo is the number of redundant total geometric derivatives
  !> \return redunt_list contains the list addresses of redundant total geometric derivatives
  subroutine geom_total_redunt_list(num_atoms, num_cent, idx_cent, order_cent, &
                                    sort_redunt, num_redunt_geo, redunt_list)
    use xkind
    implicit none
    integer, intent(in) :: num_atoms
    integer, intent(in) :: num_cent
    integer, intent(in) :: idx_cent(num_cent)
    integer, intent(in) :: order_cent(num_cent)
    logical, intent(in) :: sort_redunt
    integer, intent(in) :: num_redunt_geo
    integer, intent(out) :: redunt_list(2,num_redunt_geo)
!f2py intent(in) :: num_atoms
!f2py intent(hide) :: num_cent
!f2py intent(in) :: idx_cent
!f2py intent(in) :: order_cent
!f2py depend(num_cent) :: order_cent
!f2py intent(in) :: sort_redunt
!f2py intent(in) :: num_redunt_geo
!f2py intent(out) :: redunt_list
!f2py depend(num_redunt_geo) :: redunt_list
    integer pnum_redunt_geo                   !number of redundant total geometric derivatives (private)
    integer order_geo                         !order of total geometric derivatives
    integer num_unique_geo                    !number of unique total geometric derivatives
    integer, allocatable :: base_unique(:,:)  !base addresses of unique total geometric derivatives
    integer num_coord                         !number of atomic coordinates
    integer, allocatable :: powers_nc(:)      !powers of \var(num_coord)
    integer icent                             !incremental recorders
    integer iderv, jderv
    integer ixyz
    logical do_perm                           !indicates if there exits next permutation
    integer addr_redunt                       !address of redundant total geometric derivatives
    integer ierr                              !error information
#if defined(XTIME)
    real(REALK) curr_time                     !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! checks the number of redundant total geometric derivatives
    call geom_total_num_redunt(num_cent, order_cent, pnum_redunt_geo)
    if (pnum_redunt_geo/=num_redunt_geo)                                           &
      call error_stop("geom_total_redunt_list",                                    &
                      "incorrect number of redundant total geometric derivatives", &
                      pnum_redunt_geo)
    ! gets the order of total geometric derivatives
    order_geo = sum(order_cent)
    ! gets the number of unique total geometric derivatives
    num_unique_geo = (order_cent(1)+1)*(order_cent(1)+2)/2
    do icent = 2, num_cent
      num_unique_geo = num_unique_geo*(order_cent(icent)+1)*(order_cent(icent)+2)/2
    end do
    ! gets the base addresses of unique total geometric derivatives for the path
    allocate(base_unique(order_geo,num_unique_geo), stat=ierr)
    if (ierr/=0)                                        &
      call error_stop("geom_total_redunt_list",         &
                      "failed to allocate base_unique", &
                      order_geo*num_unique_geo)
    call geom_total_base_unique(num_cent, idx_cent, order_cent, &
                                order_geo, num_unique_geo, base_unique)
    ! sets the number of atomic coordinates
    num_coord = 3*num_atoms
    ! sets the powers of \var(num_coord)
    allocate(powers_nc(order_geo), stat=ierr)
    if (ierr/=0)                                &
      call error_stop("geom_total_redunt_list", &
                      "failed to allocate powers_nc", order_geo)
    powers_nc(1) = 1
    do ixyz = 2, order_geo
      powers_nc(ixyz) = num_coord*powers_nc(ixyz-1)
    end do
    ! loops over the unique total geometric derivatives
    jderv = 0
    do iderv = 1, num_unique_geo
      do_perm = .true.
      do while (do_perm)
        ! gets the address of redundant total geometric derivatives
        addr_redunt = base_unique(1,iderv)+1
        do ixyz = 2, order_geo
          addr_redunt = addr_redunt+base_unique(ixyz,iderv)*powers_nc(ixyz)
        end do
        jderv = jderv+1
        redunt_list(1,jderv) = iderv
        redunt_list(2,jderv) = addr_redunt
        ! gets the next posible redundant total geometric derivatives
        call next_permutation(1, order_geo, base_unique(:,iderv), do_perm)
      end do
    end do
    deallocate(base_unique)
    deallocate(powers_nc)
    ! sorts the list addresses according to the redundant total geometric derivatives
    if (sort_redunt) then
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "geom_total_redunt_list", STDOUT)
#endif
    return
  end subroutine geom_total_redunt_list

  !> \brief returns the expectation values of redundant total geometric derivatives
  !> \author Bin Gao
  !> \date 2011-11-06
  !> \param num_atoms is the number of atoms
  !> \param num_cent is the number of differentiated centers
  !> \param idx_cent contains the indices of differentiated centers
  !> \param order_cent contains the orders of differentiated centers
  !> \param num_opt is the number of operators
  !> \param num_unique_geo is the number of unique total geometric derivatives
  !> \param num_dens is the number of AO density matrices
  !> \param unique_expt contains the expectation values of unique total geometric derivatives
  !> \param dim_redunt_geo is the size of redundant total geometric derivatives,
  !>        equals to \var(3*num_atoms)^sum(\var(order_cent))
  !> \return redunt_expt contains the updated expectation values of redundant total
  !>         geometric derivatives on exit
  subroutine geom_total_redunt_expectation(num_atoms, num_cent, idx_cent, order_cent, &
                                           num_opt, num_unique_geo, num_dens,         &
                                           unique_expt, dim_redunt_geo, redunt_expt)
    use xkind
    implicit none
    integer, intent(in) :: num_atoms
    integer, intent(in) :: num_cent
    integer, intent(in) :: idx_cent(num_cent)
    integer, intent(in) :: order_cent(num_cent)
    integer, intent(in) :: num_opt
    integer, intent(in) :: num_unique_geo
    integer, intent(in) :: num_dens
    real(REALK), intent(in) :: unique_expt(num_opt,num_unique_geo,num_dens)
    integer, intent(in) :: dim_redunt_geo
    real(REALK), intent(inout) :: redunt_expt(num_opt,dim_redunt_geo,num_dens)
!f2py intent(in) :: num_atoms
!f2py intent(hide) :: num_cent
!f2py intent(in) :: idx_cent
!f2py intent(in) :: order_cent
!f2py depend(num_cent) :: order_cent
!f2py intent(hide) :: num_opt
!f2py intent(hide) :: num_unique_geo
!f2py intent(hide) :: num_dens
!f2py intent(in) :: unique_expt
!f2py intent(hide) :: dim_redunt_geo
!f2py intent(inout) :: redunt_expt
!f2py depend(num_opt) :: redunt_expt
!f2py depend(num_dens) :: redunt_expt
    integer pdim_redunt_geo                   !size of redundant total geometric derivatives (private)
    integer order_geo                         !order of total geometric derivatives
    integer pnum_unique_geo                   !number of unique total geometric derivatives
    integer, allocatable :: base_unique(:,:)  !base addresses of unique total geometric derivatives
    integer num_coord                         !number of atomic coordinates
    integer, allocatable :: powers_nc(:)      !powers of \var(num_coord)
    integer icent                             !incremental recorders
    integer iderv
    integer ixyz
    logical do_perm                           !indicates if there exits next permutation
    integer addr_redunt                       !address of redundant total geometric derivatives
    integer ierr                              !error information
#if defined(XTIME)
    real(REALK) curr_time                     !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! gets the order of total geometric derivatives
    order_geo = sum(order_cent)
    ! checks the number of unique total geometric derivatives
    pnum_unique_geo = (order_cent(1)+1)*(order_cent(1)+2)/2
    do icent = 2, num_cent
      pnum_unique_geo = pnum_unique_geo*(order_cent(icent)+1)*(order_cent(icent)+2)/2
    end do
    if (pnum_unique_geo/=num_unique_geo)                                        &
      call error_stop("geom_total_redunt_expectation",                          &
                      "incorrect number of unique total geometric derivatives", &
                      pnum_unique_geo)
    ! sets the number of atomic coordinates
    num_coord = 3*num_atoms
    ! checks the size of redundant total geometric derivatives
    pdim_redunt_geo = num_coord**order_geo
    if (pdim_redunt_geo/=dim_redunt_geo)                                         &
      call error_stop("geom_total_redunt_expectation",                           &
                      "incorrect size of redundant total geometric derivatives", &
                      pdim_redunt_geo)
    ! gets the base addresses of unique total geometric derivatives for the path
    allocate(base_unique(order_geo,num_unique_geo), stat=ierr)
    if (ierr/=0)                                        &
      call error_stop("geom_total_redunt_expectation",  &
                      "failed to allocate base_unique", &
                      order_geo*num_unique_geo)
    call geom_total_base_unique(num_cent, idx_cent, order_cent, &
                                order_geo, num_unique_geo, base_unique)
    ! sets the powers of \var(num_coord)
    allocate(powers_nc(order_geo), stat=ierr)
    if (ierr/=0)                                       &
      call error_stop("geom_total_redunt_expectation", &
                      "failed to allocate powers_nc", order_geo)
    powers_nc(1) = 1
    do ixyz = 2, order_geo
      powers_nc(ixyz) = num_coord*powers_nc(ixyz-1)
    end do
    ! puts the expectation values into appropriate positions
    do iderv = 1, num_unique_geo
      do_perm = .true.
      do while (do_perm)
        ! gets the address of redundant total geometric derivatives
        addr_redunt = base_unique(1,iderv)+1
        do ixyz = 2, order_geo
          addr_redunt = addr_redunt+base_unique(ixyz,iderv)*powers_nc(ixyz)
        end do
        redunt_expt(:,addr_redunt,:) = redunt_expt(:,addr_redunt,:) &
                                     + unique_expt(:,iderv,:)
        ! gets the next posible redundant total geometric derivatives
        call next_permutation(1, order_geo, base_unique(:,iderv), do_perm)
      end do
    end do
    deallocate(base_unique)
    deallocate(powers_nc)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "geom_total_redunt_expectation", STDOUT)
#endif
    return
  end subroutine geom_total_redunt_expectation

  !-!> \brief reorders partial geometric derivatives
  !-!> \author Bin Gao
  !-!> \date 2012-08-30
  !-!> \param order_geo is the order of partial geometric derivatives
  !-!> \param dim_ints is the dimension of integrals
  !-!> \param dim_opt is the dimension of operators
  !-!> \return contr_ints contains the contracted integrals
  !-subroutine geom_part_reorder(order_geo, dim_ints, dim_opt, contr_ints)
  !-  use xkind
  !-  implicit none
  !-  integer, intent(in) :: order_geo
  !-  integer, intent(in) :: dim_ints
  !-  integer, intent(in) :: dim_opt
  !-  real(REALK), intent(inout) :: contr_ints(dim_ints,(order_geo+1)*(order_geo+2)/2,dim_opt)
  !-  real(REALK), allocatable :: tmp_ints(:,:)  !temporary integrals
  !-  integer twice_geo_tri                      !twice of the order of geometric derivatives +3
  !-  integer order_x, order_y, order_z          !orders of xyz components
  !-  integer iopt                               !incremental recorder over operators
  !-  integer igeo                               !incremental recorder over geometric derivatives
  !-  integer ierr                               !error information
  !-  allocate(tmp_ints(dim_ints,(order_geo+1)*(order_geo+2)/2), stat=ierr)
  !-  if (ierr/=0)                                     &
  !-    call error_stop("geom_part_reorder",           &
  !-                    "failed to allocate tmp_ints", &
  !-                    dim_ints*(order_geo+1)*(order_geo+2)/2)
  !-  twice_geo_tri = order_geo+order_geo+3
  !-  do iopt = 1, dim_opt
  !-    tmp_ints = contr_ints(:,:,iopt)
  !-    ! generates derivatives, for instance, in the order of xx, xy, xz, yy, yz, zz
  !-    igeo = 0
  !-    do order_x = order_geo, 0, -1
  !-      do order_y = order_geo-order_x, 0, -1
  !-        igeo = igeo+1
  !-        order_z = order_geo-(order_x+order_y)
  !-        ! index of x^{l}y^{m}z^{n} with l+m+n=angm is 1+m+(2*angm+3-n)*n/2
  !-        contr_ints(:,igeo,iopt) = tmp_ints(:,1+order_y+(twice_geo_tri-order_z)*order_z/2)
  !-      end do
  !-    end do
  !-  end do
  !-  deallocate(tmp_ints)
  !-  return
  !-end subroutine geom_part_reorder
