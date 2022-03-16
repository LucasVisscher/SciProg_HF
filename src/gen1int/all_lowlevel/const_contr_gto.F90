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
!!  This file constructs the contracted GTOs.
!!
!!  2012-03-10, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief constructs the contracted Cartesian or Hermite Gaussians
  !> \author Bin Gao
  !> \date 2012-03-10
  !> \param num_contr_bra is the number of contractions of bra center
  !> \param num_prim_bra is the number of primitive Gaussians of bra center
  !> \param contr_coef_bra contains the contraction coefficients of bra center
  !> \param num_gto_bra is the number of Cartesian/Hermite GTOs on bra
  !> \param num_opt is the number of operators including derivatives
  !> \param prim_gto contains the primitive Cartesian/Hermite Gaussians
  !> \return contr_gto contains the contracted Cartesian/Hermite Gaussians
  subroutine const_contr_gto(num_contr_bra, num_prim_bra, contr_coef_bra, &
                             num_gto_bra, num_opt, prim_gto, contr_gto)
    use xkind
    implicit none
    integer, intent(in) :: num_contr_bra
    integer, intent(in) :: num_prim_bra
    real(REALK), intent(in) :: contr_coef_bra(num_contr_bra,num_prim_bra)
    integer, intent(in) :: num_gto_bra
    integer, intent(in) :: num_opt
    real(REALK), intent(in) :: prim_gto(num_gto_bra,num_opt,num_prim_bra)
    real(REALK), intent(out) :: contr_gto(num_gto_bra,num_contr_bra,num_opt)
!f2py intent(hide) :: num_contr_bra
!f2py intent(hide) :: num_prim_bra
!f2py intent(in) :: contr_coef_bra
!f2py intent(hide) :: num_gto_bra
!f2py intent(hide) :: num_opt
!f2py intent(in) :: prim_gto
!f2py depend(num_prim_bra) :: prim_gto
!f2py intent(out) :: contr_gto
!f2py depend(num_gto_bra) :: contr_gto
!f2py depend(num_contr_bra) :: contr_gto
!f2py depend(num_opt) :: contr_gto
    integer iprim          !incremental recorder over primitives
    integer icontr         !incremental recorder over contractions
    integer igto           !incremental recorder over Cartesian/Hermite Gaussians
    integer iopt           !incremental recorder over operators
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! zeros the contracted integrals
    contr_gto = 0.0_REALK
    ! contracted on bra
    do iprim = 1, num_prim_bra
      do iopt = 1, num_opt
        do icontr = 1, num_contr_bra
          do igto = 1, num_gto_bra
            contr_gto(igto,icontr,iopt) = contr_gto(igto,icontr,iopt) &
              + contr_coef_bra(icontr,iprim)*prim_gto(igto,iopt,iprim)
          end do
        end do
      end do
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "const_contr_gto", STDOUT)
#endif
    return
  end subroutine const_contr_gto
