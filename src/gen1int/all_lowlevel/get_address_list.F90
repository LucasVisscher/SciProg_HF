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
!!  This file gets the addresses of xyz powers/derivatives.
!!
!!  2012-09-10, Bin Gao
!!  * first version

#include "stdout.h"

  !> \brief gets the addresses of xyz powers/derivatives
  !> \author Bin Gao
  !> \date 2012-09-10
  !> \param num_derv is the number of orders of derivatives
  !> \param order_derv contains the orders of derivatives
  !> \param num_atoms is the number of atoms
  !> \param order_total_geo is the order of total geometric derivatives
  !> \param max_num_cent is the maximum number of differentiated centers
  !> \param num_order is the number of orders
  !> \param num_addr is the number of addresses
  !> \return address_list contains the addresses of xyz powers/derivatives
!FIXME: to implement partial geometric derivatives
  subroutine get_address_list(num_derv, order_derv, num_atoms, order_total_geo, &
                              max_num_cent, num_order, num_addr, address_list)
    use xkind
    implicit none
    integer, intent(in) :: num_derv
    integer, intent(in) :: order_derv(num_derv)
    integer, intent(in) :: num_atoms
    integer, intent(in) :: order_total_geo
    integer, intent(in) :: max_num_cent
    integer, intent(in) :: num_order
    integer, intent(in) :: num_addr
    integer, intent(out) :: address_list(num_order,num_addr)
!f2py intent(hide) :: num_derv
!f2py intent(in) :: order_derv
!f2py intent(in) :: num_atoms
!f2py intent(in) :: order_total_geo
!f2py intent(in) :: max_num_cent
!f2py intent(in) :: num_order
!f2py intent(in) :: num_addr
!f2py intent(out) :: address_list
!f2py depend(num_order) :: address_list
!f2py depend(num_addr) :: address_list
    integer prod_num_geo                   !product of number of (geometric) derivatives
    integer num_geo                        !number of (geometric) derivatives
    integer num_paths                      !total number of different paths
    integer visit_height                   !height of atom to visit
    integer, allocatable :: idx_node(:)    !indices of the selected atom nodes
    integer, allocatable :: wt_node(:)     !weights of the selected atom nodes
    integer, allocatable :: idx_cent(:)    !indices of generated differentiated centers
    integer, allocatable :: order_cent(:)  !orders of derivatives of the differentiated centers
    integer num_unique_geo                 !number of unique total geometric derivatives for the generated path
    integer dim_in                         !dimension of inner geometric derivatives
    integer dim_out                        !dimension of outer geometric derivatives
    integer start_addr, end_addr           !start and end addresses
    integer offset_order                   !offset of order
    integer address_derv(3)                !addresses of (geometric) derivatives
    integer ipath                          !incremental recorder over the paths
    integer igeo                           !incremental recorder over geometric derivatives
    integer ierr                           !error information
#if defined(XTIME)
    real(REALK) curr_time                  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! gets the product of number of (geometric) derivatives
    prod_num_geo = 1
    do igeo = 1, num_derv
      num_geo = (order_derv(igeo)+1)*(order_derv(igeo)+2)/2
      prod_num_geo = prod_num_geo*num_geo
    end do
    if (order_total_geo/=0) then
      ! allocates memory for total geometric derivatives
      allocate(idx_node(order_total_geo), stat=ierr)
      if (ierr/=0)                                     &
        call error_stop("get_address_list",            &
                        "failed to allocate idx_node", &
                        order_total_geo)
      allocate(wt_node(order_total_geo), stat=ierr)
      if (ierr/=0)                                    &
        call error_stop("get_address_list",           &
                        "failed to allocate wt_node", &
                        order_total_geo)
      allocate(idx_cent(max_num_cent), stat=ierr)
      if (ierr/=0)                                     &
        call error_stop("get_address_list",            &
                        "failed to allocate idx_cent", &
                        max_num_cent)
      allocate(order_cent(max_num_cent), stat=ierr)
      if (ierr/=0)                                       &
        call error_stop("get_address_list",              &
                        "failed to allocate order_cent", &
                        max_num_cent)
      ! generates the first path
      call geom_total_tree_init(num_atoms, order_total_geo, max_num_cent, &
                                num_paths, visit_height, idx_node,        &
                                wt_node, idx_cent, order_cent, num_unique_geo)
      ! generates the addresses of the first path
      dim_in = 1
      dim_out = prod_num_geo*num_unique_geo
      start_addr = 1
      end_addr = dim_out
      if (end_addr>num_addr) then
        call error_stop("get_address_list",                  &
                        "more geometric derivatives needed", &
                        end_addr)
      end if
      offset_order = 0
      address_derv = (/1,2,3/)
      do igeo = 1, num_derv
        num_geo = (order_derv(igeo)+1)*(order_derv(igeo)+2)/2
        dim_out = dim_out/num_geo
        call get_address_derv(order_derv(igeo), address_derv, offset_order, num_order, &
                              dim_in, dim_out, address_list(:,start_addr:end_addr))
        dim_in = dim_in*num_geo
        offset_order = offset_order+order_derv(igeo)
      end do
      do igeo = 1, wt_node(order_total_geo)
        address_derv(1) = 3*idx_cent(igeo)-2
        address_derv(2) = address_derv(1)+1
        address_derv(3) = address_derv(2)+1
        num_geo = (order_cent(igeo)+1)*(order_cent(igeo)+2)/2
        dim_out = dim_out/num_geo
        call get_address_derv(order_cent(igeo), address_derv, offset_order, num_order, &
                              dim_in, dim_out, address_list(:,start_addr:end_addr))
        dim_in = dim_in*num_geo
        offset_order = offset_order+order_cent(igeo)
      end do
      ! searches for the next satisfied path from a given path
      do ipath = 2, num_paths
        call geom_total_tree_search(num_atoms, order_total_geo, max_num_cent, &
                                    visit_height, idx_node, wt_node,          &
                                    idx_cent, order_cent, num_unique_geo)
        ! generates the addresses of the path
        dim_in = 1
        dim_out = prod_num_geo*num_unique_geo
        start_addr = end_addr
        end_addr = start_addr+dim_out
        if (end_addr>num_addr) then
          call error_stop("get_address_list",                  &
                          "more geometric derivatives needed", &
                          end_addr)
        end if
        start_addr = start_addr+1
        offset_order = 0
        address_derv = (/1,2,3/)
        do igeo = 1, num_derv
          num_geo = (order_derv(igeo)+1)*(order_derv(igeo)+2)/2
          dim_out = dim_out/num_geo
          call get_address_derv(order_derv(igeo), address_derv, offset_order, num_order, &
                                dim_in, dim_out, address_list(:,start_addr:end_addr))
          dim_in = dim_in*num_geo
          offset_order = offset_order+order_derv(igeo)
        end do
        do igeo = 1, wt_node(order_total_geo)
          address_derv(1) = 3*idx_cent(igeo)-2
          address_derv(2) = address_derv(1)+1
          address_derv(3) = address_derv(2)+1
          num_geo = (order_cent(igeo)+1)*(order_cent(igeo)+2)/2
          dim_out = dim_out/num_geo
          call get_address_derv(order_cent(igeo), address_derv, offset_order, num_order, &
                                dim_in, dim_out, address_list(:,start_addr:end_addr))
          dim_in = dim_in*num_geo
          offset_order = offset_order+order_cent(igeo)
        end do
      end do
      ! cleans
      deallocate(idx_node)
      deallocate(wt_node)
      deallocate(idx_cent)
      deallocate(order_cent)
    else
      dim_in = 1
      dim_out = prod_num_geo
      offset_order = 0
      address_derv = (/1,2,3/)
      do igeo = 1, num_derv
        num_geo = (order_derv(igeo)+1)*(order_derv(igeo)+2)/2
        dim_out = dim_out/num_geo
        call get_address_derv(order_derv(igeo), address_derv, offset_order, &
                              num_order, dim_in, dim_out, address_list)
        dim_in = dim_in*num_geo
        offset_order = offset_order+order_derv(igeo)
      end do
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "get_address_list", STDOUT)
#endif
    return
  end subroutine get_address_list

  !> \brief gets the addresses of xyz powers/derivatives
  !> \author Bin Gao
  !> \date 2012-09-10
  !> \param order_derv contains the orders of derivatives
  !> \param address_derv contains the addresses of derivatives
  !> \param offset_order is the offset of order
  !> \param num_order is the number of orders
  !> \param dim_in is the dimension of inner geometric derivatives
  !> \param dim_out is the dimension of outer geometric derivatives
  !> \return address_list contains the addresses of xyz powers/derivatives
  subroutine get_address_derv(order_derv, address_derv, offset_order, num_order, &
                              dim_in, dim_out, address_list)
    use xkind
    implicit none
    integer, intent(in) :: order_derv
    integer, intent(in) :: address_derv(3)
    integer, intent(in) :: offset_order
    integer, intent(in) :: num_order
    integer, intent(in) :: dim_in
    integer, intent(in) :: dim_out
    integer, intent(inout) :: address_list(num_order,dim_in,                &
                                           (order_derv+1)*(order_derv+2)/2, &
                                           dim_out)
!f2py intent(in) :: order_derv
!f2py intent(in) :: address_derv
!f2py intent(in) :: offset_order
!f2py intent(in) :: num_order
!f2py intent(in) :: dim_in
!f2py intent(in) :: dim_out
!f2py intent(inout) :: address_list
!f2py depend(dim_in) :: address_list
!f2py depend(order_derv) :: address_list
!f2py depend(dim_out) :: address_list
    integer order_x, order_y, order_z  !orders along xyz direction
    integer iaddr                      !incremental recorder over addresses
    integer iorder                     !incremental recorder over orders
#if defined(XTIME)
    real(REALK) curr_time              !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    iaddr = 0
    do order_z = 0, order_derv
      do order_y = 0, order_derv-order_z
        order_x = order_derv-(order_y+order_z)
        iaddr = iaddr+1
        ! x-direction
        do iorder = offset_order+1, offset_order+order_x
          address_list(iorder,:,iaddr,:) = address_derv(1)
        end do
        ! y-direction
        do iorder = offset_order+order_x+1, offset_order+order_derv-order_z
          address_list(iorder,:,iaddr,:) = address_derv(2)
        end do
        ! z-direction
        do iorder = offset_order+order_x+order_y+1, offset_order+order_derv
          address_list(iorder,:,iaddr,:) = address_derv(3)
        end do
      end do
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "get_address_derv", STDOUT)
#endif
    return
  end subroutine get_address_derv
