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
!!  This file computes the binomial coefficient.
!!
!!  2011-07-07, Bin Gao:
!!  * adds subroutine "pascal_triangle" to calculate the binomial coefficients
!!    in the form of Pascal's triangle
!!
!!  2009-09-23, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief computes the binomial coefficient \f$\binom{n}{k}=\frac{n!}{k!(n-k)!}\f$
  !> \author Bin Gao
  !> \date 2009-09-23
  !> \param num_n is the number \f$n\f$
  !> \param num_k is the number \f$k\f$
  !> \return bcoeff is the binomial coefficient
  subroutine binom_coeff(num_n, num_k, bcoeff)
    use xkind
    implicit none
    integer, intent(in) :: num_n
    integer, intent(in) :: num_k
    integer, intent(out) :: bcoeff
!f2py intent(in) :: num_n
!f2py intent(in) :: num_k
!f2py intent(out) :: bcoeff
    integer n_inc  !number \f$n+1\f$
    integer inum    !incremental recorder
!FIXME: how to check if the numbers do not exceed the bound
    bcoeff = 1
    n_inc = num_n+1
    if (num_k<num_n-num_k) then
      do inum = num_n, n_inc-num_k, -1
        bcoeff = bcoeff*inum/(n_inc-inum)
      end do
    else
      do inum = num_n, num_k+1, -1
        bcoeff = bcoeff*inum/(n_inc-inum)
      end do
    end if
    return
  end subroutine binom_coeff

  !> \brief computes the binomial coefficient \f$\binom{n}{k}=\frac{n!}{k!(n-k)!}\f$
  !> \author Bin Gao
  !> \date 2009-09-23
  !> \param num_n is the number \f$n\f$
  !> \param num_k is the number \f$k\f$
  !> \return bcoeff is the binomial coefficient
  subroutine dbinom_coeff(num_n, num_k, bcoeff)
    use xkind
    implicit none
    integer, intent(in) :: num_n
    integer, intent(in) :: num_k
    integer, intent(out) :: bcoeff
!f2py intent(in) :: num_n
!f2py intent(in) :: num_k
!f2py intent(out) :: bcoeff
    real(REALK) dn_inc               !number \f$n+1\f$ in double precision
    real(REALK) dnum_k                !number \f$n\f$ in double precision
    real(REALK) d_bcoeff              !binomial coefficient in double precision
    integer inum, dnum                !incremental recorders
!FIXME: how to check if the numbers do not exceed the bound
    d_bcoeff = 1.0_REALK
    dn_inc = real(num_n+1,REALK)
    dnum_k = real(num_k,REALK)
    if (num_k<num_n-num_k) then
      do inum = num_n, num_n-num_k+1, -1
        dnum = real(inum,REALK)
        d_bcoeff = d_bcoeff*dnum/(dn_inc-dnum)
      end do
    else
      do inum = num_n, num_k+1, -1
        dnum = real(inum,REALK)
        d_bcoeff = d_bcoeff*dnum/(dn_inc-dnum)
      end do
    end if
    bcoeff = d_bcoeff
    if (d_bcoeff/=real(bcoeff,REALK)) then
      write(STDOUT,100) "numbers", num_n, num_k
      call error_stop("dbinom_coeff", "numbers exceed the bound", -1)
    end if
    return
100 format("dbinom_coeff>> ",A,I6," and",I6)
  end subroutine dbinom_coeff

  !> \brief returns the Pascal's triangle using \f$\binom{n+1}{k+1}=\binom{n}{k}+\binom{n}{k+1}\f$
  !> \author Bin Gao
  !> \date 2011-07-07
  !> \param binom_order is the maximum order
  !> \return binom_coeff contains the binomial coefficients
  subroutine pascal_triangle(binom_order, binom_coeff)
    use xkind
    implicit none
    integer, intent(in) :: binom_order
    real(REALK), intent(out) :: binom_coeff(0:binom_order,0:binom_order)
!f2py intent(in) :: binom_order
!f2py intent(out) :: binom_coeff
!f2py depend(binom_order) :: binom_coeff
    integer num_n, num_k, n_inc, k_inc  !incremental recorders
    binom_coeff(0,0) = 1.0_REALK
    do num_n = 0, binom_order-1
      n_inc = num_n+1
      binom_coeff(0,n_inc) = 1.0_REALK
      do num_k = 0, num_n-1
        k_inc = num_k+1
        binom_coeff(k_inc,n_inc) = binom_coeff(num_k,num_n)+binom_coeff(k_inc,num_n)
      end do
      binom_coeff(n_inc,n_inc) = 1.0_REALK
    end do
    return
  end subroutine pascal_triangle
