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
!!  This file scatters higher order shells into two lower order shells.
!! 
!!  2011-06-09, Bin Gao
!!  * first version

#include "stdout.h"

  !> \brief scatters higher order shells into two lower order shells
  !> \author Bin Gao
  !> \date 2011-06-09
  !> \param dim_ints is the dimension of integrals
  !> \param dim_gather is the dimension of higher order shells
  !> \param dim_between is the dimension of shells between these two lower order shells
  !> \param dim_outmost is the dimension of other outmost shells
  !> \param gather_ints contains the integrals of higher order shells
  !> \param order_inner is the order of inner lower order shell
  !> \param order_outer is the order of outer lower order shell
  !> \param num_inner is the number of xyz components of inner lower order shell,
  !>        should equal to \f$(\var(order_inner)+1)(\var(order_inner)+2)/2\f$
  !> \param num_outer is the number of xyz components of outer lower order shell,
  !>        should equal to \f$(\var(order_outer)+1)(\var(order_outer)+2)/2\f$
  !> \return scatter_ints contains the integrals of two lower order shells
  subroutine scatter_single(dim_ints, dim_gather, dim_between, dim_outmost,   &
                            gather_ints, order_inner, order_outer, num_inner, &
                            num_outer, scatter_ints)
    use xkind
    implicit none
    integer, intent(in) :: dim_ints
    integer, intent(in) :: dim_gather
    integer, intent(in) :: dim_between
    integer, intent(in) :: dim_outmost
    real(REALK), intent(in) :: gather_ints(dim_ints,dim_gather,dim_between,dim_outmost)
    integer, intent(in) :: order_inner
    integer, intent(in) :: order_outer
    integer, intent(in) :: num_inner
    integer, intent(in) :: num_outer
    real(REALK), intent(out) :: scatter_ints(dim_ints,num_inner,dim_between, &
                                             num_outer,dim_outmost)
!f2py intent(hide) :: dim_ints
!f2py intent(hide) :: dim_gather
!f2py intent(hide) :: dim_between
!f2py intent(hide) :: dim_outmost
!f2py intent(in) :: gather_ints
!f2py intent(in) :: order_inner
!f2py intent(in) :: order_outer
!f2py intent(in) :: num_inner
!f2py intent(in) :: num_outer
!f2py intent(out) :: scatter_ints
!f2py depend(dim_ints) :: scatter_ints
!f2py depend(num_inner) :: scatter_ints
!f2py depend(dim_between) :: scatter_ints
!f2py depend(num_outer) :: scatter_ints
!f2py depend(dim_outmost) :: scatter_ints
    integer start_xyz_in    !start address of xyz components in current inner shell
    integer end_xyz_in      !end address of xyz components in current inner shell
    integer xyz_min_high    !address of current xyz component in the minimum higher order shell
    integer start_xyz_high  !start address of xyz components in current higher order shell
    integer end_xyz_high    !end address of xyz components in current higher order shell
    integer step_outer      !step size of xyz components of outer lower order shell
    integer xyz_outer       !address of xyz components of outer lower order shell
    integer jinner          !incremental recorder over the components of inner lower order shell
    integer iouter, jouter  !incremental recorders over xyz components of outer lower order shell
    integer imost           !incremental recorder over components of outmost shells
#if defined(XTIME)
    real(REALK) curr_time   !current CPU time   
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! loops over outmost shells
    do imost = 1, dim_outmost
      xyz_min_high = 0
      xyz_outer = 0
      ! loops over the xyz components of outer shell
      do iouter = order_outer, 0, -1
        step_outer = iouter+1
        do jouter = 0, iouter
          start_xyz_in = 1
          xyz_min_high = xyz_min_high+1
          start_xyz_high = xyz_min_high
          xyz_outer = xyz_outer+1
          ! loops over the xyz components of inner shell
          do jinner = order_inner, 0, -1
            end_xyz_in = start_xyz_in+jinner
            end_xyz_high = start_xyz_high+jinner
#if defined(DEBUG)
            write(STDOUT,100) "inner/outer address", start_xyz_in, end_xyz_in, xyz_outer
            write(STDOUT,100) "higher address", start_xyz_high, end_xyz_high
#endif
            scatter_ints(:,start_xyz_in:end_xyz_in,:,xyz_outer,imost) &
              = gather_ints(:,start_xyz_high:end_xyz_high,:,imost)
            start_xyz_in = end_xyz_in+1
            start_xyz_high = end_xyz_high+step_outer
          end do
        end do
        ! updates the address of current xyz component in the minimum higher order shell
        xyz_min_high = xyz_min_high+order_inner
      end do
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "scatter_single", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("scatter_single>> ",A,I8," -->",2I8)
#endif
  end subroutine scatter_single

  !> \brief scatters higher order shells into one lower order outer shell,
  !>        and multiple lower order inner shells
  !> \author Bin Gao
  !> \date 2011-06-09
  !> \param dim_ints is the dimension of integrals
  !> \param dim_gather is the dimension of higher order shells
  !> \param dim_between is the dimension of shells between these two lower order shells
  !> \param dim_outmost is the dimension of other outmost shells
  !> \param gather_ints contains the integrals of higher order shells
  !> \param orders_inner contains the minimum and maximum orders of inner lower order shells
  !> \param order_outer is the order of outer lower order shell
  !> \param dim_inner is the dimension of inner lower order shells
  !> \param num_outer is the number of xyz components of outer lower order shell,
  !>        should equal to \f$(\var(order_outer)+1)(\var(order_outer)+2)/2\f$
  !> \return scatter_ints contains the integrals of two lower order shells
  subroutine scatter_multi_inner(dim_ints, dim_gather, dim_between, dim_outmost,    &
                                 gather_ints, orders_inner, order_outer, dim_inner, &
                                 num_outer, scatter_ints)
    use xkind
    implicit none
    integer, intent(in) :: dim_ints
    integer, intent(in) :: dim_gather
    integer, intent(in) :: dim_between
    integer, intent(in) :: dim_outmost
    real(REALK), intent(in) :: gather_ints(dim_ints,dim_gather,dim_between,dim_outmost)
    integer, intent(in) :: orders_inner(2)
    integer, intent(in) :: order_outer
    integer, intent(in) :: dim_inner
    integer, intent(in) :: num_outer
    real(REALK), intent(out) :: scatter_ints(dim_ints,dim_inner,dim_between, &
                                             num_outer,dim_outmost)
!f2py intent(hide) :: dim_ints
!f2py intent(hide) :: dim_gather
!f2py intent(hide) :: dim_between
!f2py intent(hide) :: dim_outmost
!f2py intent(in) :: gather_ints
!f2py intent(in) :: orders_inner
!f2py intent(in) :: order_outer
!f2py intent(in) :: dim_inner
!f2py intent(in) :: num_outer
!f2py intent(out) :: scatter_ints
!f2py depend(dim_ints) :: scatter_ints
!f2py depend(dim_inner) :: scatter_ints
!f2py depend(dim_between) :: scatter_ints
!f2py depend(num_outer) :: scatter_ints
!f2py depend(dim_outmost) :: scatter_ints
    integer start_xyz_in    !start address of xyz components in current inner shell
    integer end_xyz_in      !end address of xyz components in current inner shell
    integer xyz_min_high    !address of current xyz component in the minimum higher order shell
    integer start_higher    !start address of current higher order shell
    integer start_xyz_high  !start address of xyz components in current higher order shell
    integer end_xyz_high    !end address of xyz components in current higher order shell
    integer step_inner      !step size between the addresses of different inner lower order shells
    integer step_outer      !step size of xyz components of outer lower order shell
    integer xyz_outer       !address of xyz components of outer lower order shell
    integer order_inner     !order of current inner lower order shell
    integer jinner          !incremental recorder over the components of inner lower order shell
    integer order_high      !order of current higher order shell
    integer iouter, jouter  !incremental recorders over xyz components of outer lower order shell
    integer imost           !incremental recorder over components of outmost shells
#if defined(XTIME)
    real(REALK) curr_time   !current CPU time   
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! loops over outmost shells
    do imost = 1, dim_outmost
      xyz_min_high = 0
      step_inner = 0
      xyz_outer = 0
      ! loops over the xyz components of outer shell
      do iouter = order_outer, 0, -1
        step_outer = iouter+1
        do jouter = 0, iouter
          ! for each component of outer shell, loops over different inner shells
          start_xyz_in = 1
          xyz_min_high = xyz_min_high+1
          start_higher = xyz_min_high
          xyz_outer = xyz_outer+1
          ! loops over different orders of inner shells
          do order_inner = orders_inner(1), orders_inner(2)
            order_high = order_inner+order_outer
#if defined(DEBUG)
            write(STDOUT,100) "higher/inner orders", order_high, order_inner
#endif
            start_xyz_high = start_higher
            do jinner = order_inner, 0, -1
              end_xyz_in = start_xyz_in+jinner
              end_xyz_high = start_xyz_high+jinner
#if defined(DEBUG)
              write(STDOUT,100) "inner/outer address", start_xyz_in, end_xyz_in, xyz_outer
              write(STDOUT,100) "higher address", start_xyz_high, end_xyz_high
#endif
              scatter_ints(:,start_xyz_in:end_xyz_in,:,xyz_outer,imost) &
                = gather_ints(:,start_xyz_high:end_xyz_high,:,imost)
              start_xyz_in = end_xyz_in+1
              start_xyz_high = end_xyz_high+step_outer
            end do
            ! gets the start address of next higher order shell
            start_higher = start_higher+(order_high+1)*(order_high+2)/2+step_inner
          end do
        end do
        ! updates the address of current xyz component in the minimum higher order shell
        xyz_min_high = xyz_min_high+orders_inner(1)
        step_inner = step_inner+1
      end do
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "scatter_multi_inner", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("scatter_multi_inner>> ",A,I8," -->",2I8)
#endif
  end subroutine scatter_multi_inner
