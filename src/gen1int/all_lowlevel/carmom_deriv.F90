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
!!  This file calculates the geometric derivatives of Cartesian multipole moment integrals.
!!
!!  2011-07-22, Bin Gao:
!!  * rewrites at least four times, ..., I should write down the times of rewriting ...
!!  * the sequence of dimensions of \var(geo_mom_cints) might be a problem, in particular,
!!    the dimensions \var(num_geo_mom) and \var(num_mom)
!!
!!  2011-06-17, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief computes the geometric derivatives of Cartesian multipole moment center
  !> \author Bin Gao
  !> \date 2011-06-17
  !> \param order_low_mom is the order of lower order Cartesian multipole moment
  !> \param dim_cints is the dimension of contracted integrals
  !> \param num_low_mom is the number of xyz components of lower order Cartesian multipole moment,
  !>        should be \f$(\var(order_low_mom)+1)(\var(order_low_mom)+2)/2\f$
  !> \param num_opt is the number of other operators
  !> \param low_mom_cints contains the contracted integrals of lower order artesian multipole moment
  !> \param order_mom is the order of Cartesian multipole moment, should be
  !>        \f$\var(order_low_mom)+\var(order_geo_mom)\f$
  !> \param order_geo_mom is the order of geometric derivatives
  !> \param num_mom is the number of xyz components of Cartesian multipole moment,
  !>        should be \f$(\var(order_mom)+1)(\var(order_mom)+2)/2\f$
  !> \param num_geo_mom is the number of xyz components of geometric derivatives,
  !>        should be \f$(\var(order_geo_mom)+1)(\var(order_geo_mom)+2)/2\f$
  !> \return geo_mom_cints contains the contracted integrals of geometric derivatives of
  !>         Cartesian multipole moment
  subroutine carmom_deriv(order_low_mom, dim_cints, num_low_mom, num_opt, low_mom_cints, &
                          order_mom, order_geo_mom, num_mom, num_geo_mom, geo_mom_cints)
    use xkind
    implicit none
    integer, intent(in) :: order_low_mom
    integer, intent(in) :: dim_cints
    integer, intent(in) :: num_low_mom
    integer, intent(in) :: num_opt
    real(REALK), intent(in) :: low_mom_cints(dim_cints,num_low_mom,num_opt)
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_geo_mom
    integer, intent(in) :: num_mom
    integer, intent(in) :: num_geo_mom
    real(REALK), intent(out) :: geo_mom_cints(dim_cints,num_mom,num_opt,num_geo_mom)
!f2py intent(in) :: order_low_mom
!f2py intent(hide) :: dim_cints
!f2py intent(hide) :: num_low_mom
!f2py intent(hide) :: num_opt
!f2py intent(in) :: low_mom_cints
!f2py intent(in) :: order_mom
!f2py intent(in) :: order_geo_mom
!f2py intent(in) :: num_mom
!f2py intent(in) :: num_geo_mom
!f2py intent(out) :: geo_mom_cints
!f2py depend(dim_cints) :: geo_mom_cints
!f2py depend(num_mom) :: geo_mom_cints
!f2py depend(num_opt) :: geo_mom_cints
!f2py depend(num_geo_mom) :: geo_mom_cints
    real(REALK), allocatable :: coeff_mom(:,:)  !coefficients of geometric derivatives of
                                                !Cartesian multipole moment
    integer x_mom, y_mom, z_mom  !order of Cartesian multipole moments
    integer x_geo, y_geo, z_geo  !order of geometric derivatives
    integer x_low, y_low, z_low  !order of lower order Cartesian multipole moments
    integer addr_mom             !address of Cartesian multipole moment
    integer addr_geo             !address of geometric derivatives
    integer addr_low             !address of lower order Cartesian multipole moment
    real(REALK) prod_coeff       !prodcut coefficient from xyz components
    integer iopt                 !incremental recorder over other operators
    integer ierr                 !error information
#if defined(XTIME)
    real(REALK) curr_time        !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! checks the validity of the orders
    if (order_mom/=order_low_mom+order_geo_mom) &
      call error_stop("carmom_deriv", "invalid order_mom", order_mom)
    ! considers the following geometric derivative of Cartesian multipole moment:
    ! \f$\frac{\partial^{L_{X}}}{\partial X_{M}^{L_{X}}}
    !    \frac{\partial^{L_{Y}}}{\partial Y_{M}^{L_{Y}}}
    !    \frac{\partial^{L_{Z}}}{\partial Z_{M}^{L_{Z}}}
    !    x_{M}^{m_{x}}y_{M}^{m_{y}}z_{M}^{m_{z}}\f$
    ! (where \f$m_{\xi}\ge L_{\xi}\f$, \f$\xi=x,y,z\f$), which gives the coefficient
    ! \f$(-1)^{L_{X}+L_{Y}+L_{Z}}\frac{m_{x}!}{(m_{x}-L_{X})!}
    !    \frac{m_{y}!}{(m_{y}-L_{Y})!}\frac{m_{z}!}{(m_{z}-L_{Z})!}\f$;
    ! as regards \f$x\f$ direction, we have the coefficient coeff_mom(m_x,L_X) as
    ! \f$(-1)^{L_{X}}\frac{m_{x}!}{(m_{x}-L_{X})!}\f$, and the following recurrence relation
    ! coeff_mom(m_x,L_X+1)=-(m_x-L_X)*coeff_mom(m_x,L_X)=(L_X-m_x)*coeff_mom(m_x,L_X),
    ! coeff_mom(:,0)=1, and coeff_mom(:,1)=-m_x
    allocate(coeff_mom(0:order_mom,0:order_geo_mom), stat=ierr)
    if (ierr/=0)                                                      &
      call error_stop("carmom_deriv", "failed to allocate coeff_mom", &
                      (order_mom+1)*(order_geo_mom+1))
    ! prepares all the possible coefficients according to the above relations
    coeff_mom(0,0) = 1.0_REALK
    do x_mom = 1, order_geo_mom
      coeff_mom(x_mom,0) = 1.0_REALK
      coeff_mom(x_mom,1) = -real(x_mom,REALK)
      do x_geo = 1, x_mom-1
        coeff_mom(x_mom,x_geo+1) = real(x_geo-x_mom,REALK)*coeff_mom(x_mom,x_geo)
      end do
    end do
    do x_mom = order_geo_mom+1, order_mom
      coeff_mom(x_mom,0) = 1.0_REALK
      coeff_mom(x_mom,1) = -real(x_mom,REALK)
      do x_geo = 1, order_geo_mom-1
        coeff_mom(x_mom,x_geo+1) = real(x_geo-x_mom,REALK)*coeff_mom(x_mom,x_geo)
      end do
    end do
    ! loops over the xyz componets of geometric derivative
    addr_geo = 0
    do z_geo = 0, order_geo_mom
      do y_geo = 0, order_geo_mom-z_geo
        x_geo = order_geo_mom-(z_geo+y_geo)
        addr_geo = addr_geo+1
        !FIXME: this may be more efficient?
        !addr_geo = 1+y_geo+(2*order_geo_mom+3-z_geo)*z_geo/2
        ! loops over other operators
        do iopt = 1, num_opt
          ! loops over the xyz components of Cartesian multipole moment
          addr_mom = 0
          do z_mom = 0, order_mom
            do y_mom = 0, order_mom-z_mom
              x_mom = order_mom-(z_mom+y_mom)
              addr_mom = addr_mom+1
              !FIXME: this may be more efficient?
              !addr_mom = 1+y_mom+(2*order_mom+3-z_mom)*z_mom/2
              ! the corresponding orders of lower order Cartesian multipole moment
              x_low = x_mom-x_geo
              y_low = y_mom-y_geo
              z_low = z_mom-z_geo
              if (x_low<0 .or. y_low<0 .or. z_low<0) then
                geo_mom_cints(:,addr_mom,iopt,addr_geo) = 0.0_REALK
              else
                ! sets the product coefficient
                prod_coeff = coeff_mom(x_mom,x_geo)*coeff_mom(y_mom,y_geo) &
                           * coeff_mom(z_mom,z_geo)
                ! index of x^{l}y^{m}z^{n} with l+m+n=angm is 1+m+(2*angm+3-n)*n/2
                addr_low = 1+y_low+(2*order_low_mom+3-z_low)*z_low/2
#if defined(DEBUG)
                ! dumps to check
                write(STDOUT,100) "place (MOM,GEO,LOW)", addr_mom, addr_geo, addr_low
                write(STDOUT,100) "order of MOM:      ", x_mom, y_mom, z_mom
                write(STDOUT,100) "order of GEO:      ", x_geo, y_geo, z_geo
                write(STDOUT,110) "coeffcient:        ", prod_coeff
#endif
                geo_mom_cints(:,addr_mom,iopt,addr_geo) &
                  = prod_coeff*low_mom_cints(:,addr_low,iopt)
              end if
            end do
          end do
        end do
      end do
    end do
    ! cleans
    deallocate(coeff_mom)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "carmom_deriv", STDOUT)
#endif
#if defined(DEBUG)
100 format("carmom_deriv>> ",A,3I8)
110 format("carmom_deriv>> ",A,Es16.8)
#endif
    return
  end subroutine carmom_deriv
