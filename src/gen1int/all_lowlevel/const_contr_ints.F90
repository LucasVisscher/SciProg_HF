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
!!  This file constructs the contracted integrals.
!!
!!  2011-07-24, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief constructs the contracted Cartesian or Hermite integrals
  !> \author Bin Gao
  !> \date 2011-07-24
  !> \param num_contr_bra is the number of contractions of bra center
  !> \param num_prim_bra is the number of primitive Gaussians of bra center
  !> \param contr_coef_bra contains the contraction coefficients of bra center
  !> \param num_contr_ket is the number of contractions of ket center
  !> \param num_prim_ket is the number of primitive Gaussians of ket center
  !> \param contr_coef_ket contains the contraction coefficients of ket center
  !> \param num_gto_bra is the number of Cartesian/Hermite GTOs on bra
  !> \param num_gto_ket is the number of Cartesian/Hermite GTOs on ket
  !> \param num_opt is the number of operators including derivatives
  !> \param prim_ints contains the primitive Cartesian/Hermite integrals
  !> \return contr_ints contains the contracted Cartesian/Hermite integrals
  subroutine const_contr_ints(num_contr_bra, num_prim_bra, contr_coef_bra, &
                              num_contr_ket, num_prim_ket, contr_coef_ket, &
                              num_gto_bra, num_gto_ket, num_opt,           &
                              prim_ints, contr_ints)
    use xkind
    implicit none
    integer, intent(in) :: num_contr_bra
    integer, intent(in) :: num_prim_bra
    real(REALK), intent(in) :: contr_coef_bra(num_contr_bra,num_prim_bra)
    integer, intent(in) :: num_contr_ket
    integer, intent(in) :: num_prim_ket
    real(REALK), intent(in) :: contr_coef_ket(num_contr_ket,num_prim_ket)
    integer, intent(in) :: num_gto_bra
    integer, intent(in) :: num_gto_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(in) :: prim_ints(num_gto_bra,num_gto_ket,num_opt, &
                                         num_prim_bra,num_prim_ket)
    real(REALK), intent(out) :: contr_ints(num_gto_bra,num_contr_bra, &
                                           num_gto_ket,num_contr_ket,num_opt)
!f2py intent(hide) :: num_contr_bra
!f2py intent(hide) :: num_prim_bra
!f2py intent(in) :: contr_coef_bra
!f2py intent(hide) :: num_contr_ket
!f2py intent(hide) :: num_prim_ket
!f2py intent(in) :: contr_coef_ket
!f2py intent(hide) :: num_gto_bra
!f2py intent(hide) :: num_gto_ket
!f2py intent(hide) :: num_opt
!f2py intent(in) :: prim_ints
!f2py depend(num_prim_bra) :: prim_ints
!f2py depend(num_prim_ket) :: prim_ints
!f2py intent(out) :: contr_ints
!f2py depend(num_gto_bra) :: contr_ints
!f2py depend(num_contr_bra) :: contr_ints
!f2py depend(num_gto_ket) :: contr_ints
!f2py depend(num_contr_ket) :: contr_ints
!f2py depend(num_opt) :: contr_ints
    real(REALK), allocatable :: bra_cints(:,:,:,:,:)
                            !intermediate integrals by contracting on bra
    integer iprim, jprim    !incremental recorders over primitives
    integer icontr, jcontr  !incremental recorders over contractions
    integer igto, jgto      !incremental recorders over Cartesian/Hermite Gaussians
    integer iopt            !incremental recorder over operators
    integer ierr            !error information
#if defined(XTIME)
    real(REALK) curr_time   !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! zeros the contracted integrals
    contr_ints = 0.0_REALK
    ! allocates the memory for the intermediate integrals by contracting on bra
    allocate(bra_cints(num_gto_bra,num_contr_bra,num_gto_ket, &
                       num_opt,num_prim_ket), stat=ierr)
    if (ierr/=0)                                                                     &
      call error_stop("const_contr_ints", "failed to allocate memory for bra_cints", &
                      num_gto_bra*num_contr_bra*num_gto_ket*num_opt*num_prim_ket)
    bra_cints = 0.0_REALK
    ! contracted on bra
    do jprim = 1, num_prim_ket
      do iprim = 1, num_prim_bra
        do iopt = 1, num_opt
          do jgto = 1, num_gto_ket
            do icontr = 1, num_contr_bra
              do igto = 1, num_gto_bra
                bra_cints(igto,icontr,jgto,iopt,jprim)     &
                  = bra_cints(igto,icontr,jgto,iopt,jprim) &
                  + contr_coef_bra(icontr,iprim)           &
                  * prim_ints(igto,jgto,iopt,iprim,jprim)
              end do
            end do
          end do
        end do
      end do
    end do
    ! contracted on ket
    do jprim = 1, num_prim_ket
      do iopt = 1, num_opt
        do jcontr = 1, num_contr_ket
          do jgto = 1, num_gto_ket
            do icontr = 1, num_contr_bra
              do igto = 1, num_gto_bra
                contr_ints(igto,icontr,jgto,jcontr,iopt)     &
                  = contr_ints(igto,icontr,jgto,jcontr,iopt) &
                  + contr_coef_ket(jcontr,jprim)             &
                  * bra_cints(igto,icontr,jgto,iopt,jprim)
              end do
            end do
          end do
        end do
      end do
    end do
    ! cleans
    deallocate(bra_cints)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "const_contr_ints", STDOUT)
#endif
    return
  end subroutine const_contr_ints
