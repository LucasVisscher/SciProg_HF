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
!!  This file computes the trace of integrals with atomic-orbital density matrix.
!!
!!  2011-08-03, Bin Gao
!!  * revised for contracted GTOs
!! 
!!  2010-10-07, Bin Gao
!!  * first version

#include "stdout.h"

  !> \brief reorders integrals of contracted real solid-harmonic Gaussians, and computes
  !>        their traces with atomic-orbital density matrix at the same time
  !> \details \f$tr(AD) = \sum_{ij}A_{ij}D_{ji}\f$
  !> \author Bin Gao
  !> \date 2010-10-07
  !> \param ang_bra is the orbital quantum number (or angular number) of bra center
  !> \param num_sgto_bra is the number of basis functions on bra center
  !> \param mag_bra contains the magnetic numbers of basis functions on bra center
  !> \param ang_ket is the orbital quantum number (or angular number) of ket center
  !> \param num_sgto_ket is the number of basis functions on ket center
  !> \param mag_ket contains the magnetic numbers of basis functions on ket center
  !> \param num_contr_bra is the number of contractions of bra center
  !> \param num_contr_ket is the number of contractions of ket center
  !> \param num_opt is the number of operators
  !> \param gen_ints contains the contracted integrals from Gen1Int
  !> \param ao_dens contains the corresponding blocks of the atomic-orbital density matrix
  !> \return tr_ints contains the traces of integrals with atomic-orbital density matrix
  subroutine trace_sgto_ints(ang_bra, num_sgto_bra, mag_bra, &
                             ang_ket, num_sgto_ket, mag_ket, &
                             num_contr_bra, num_contr_ket,   &
                             num_opt, gen_ints, ao_dens, tr_ints)
    use xkind
    implicit none
    integer, intent(in) :: ang_bra
    integer, intent(in) :: num_sgto_bra
    integer, intent(in) :: mag_bra(num_sgto_bra)
    integer, intent(in) :: ang_ket
    integer, intent(in) :: num_sgto_ket
    integer, intent(in) :: mag_ket(num_sgto_ket)
    integer, intent(in) :: num_contr_bra
    integer, intent(in) :: num_contr_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(in) :: gen_ints(num_sgto_bra,num_contr_bra, &
                                        num_sgto_ket,num_contr_ket,num_opt)
    real(REALK), intent(in) :: ao_dens(num_sgto_ket,num_contr_ket, &
                                       num_sgto_bra,num_contr_bra)
    real(REALK), intent(inout) :: tr_ints(num_opt)
!f2py intent(in) :: ang_bra
!f2py intent(hide) :: num_sgto_bra
!f2py intent(in) :: mag_bra
!f2py intent(in) :: ang_ket
!f2py intent(hide) :: num_sgto_ket
!f2py intent(in) :: mag_ket
!f2py intent(hide) :: num_contr_bra
!f2py intent(hide) :: num_contr_ket
!f2py intent(hide) :: num_opt
!f2py intent(in) :: gen_ints
!f2py depend(num_sgto_bra) :: gen_ints
!f2py depend(num_sgto_ket) :: gen_ints
!f2py intent(in) :: ao_dens
!f2py depend(num_sgto_ket) :: ao_dens
!f2py depend(num_contr_ket) :: ao_dens
!f2py depend(num_sgto_bra) :: ao_dens
!f2py depend(num_contr_bra) :: ao_dens
!f2py intent(inout) :: tr_ints
!f2py depend(num_opt) :: tr_ints
    integer ang_bra_one, ang_ket_one          !orbital quantum numbers of bra and ket +1
    integer ibra, iket, icontr, jcontr, iopt  !increment recorders
    integer addr_bra, addr_ket                !addresses of SGTOs on bra and ket centers in Gen1Int
#if defined(XTIME)
    real(REALK) curr_time                     !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ang_bra_one = ang_bra+1
    ang_ket_one = ang_ket+1
    ! loops over operators
    do iopt = 1, num_opt
      ! loops over contractions on ket center
      do jcontr = 1, num_contr_ket
        ! loops over SGTOs on ket center
        do iket = 1, num_sgto_ket
          addr_ket = mag_ket(iket)+ang_ket_one
          ! loops over contractions on bra center
          do icontr = 1, num_contr_bra
            ! loops over SGTOs on bra center
            do ibra = 1, num_sgto_bra
              addr_bra = mag_bra(ibra)+ang_bra_one
              tr_ints(iopt) = tr_ints(iopt)                                  &
                            + gen_ints(addr_bra,icontr,addr_ket,jcontr,iopt) &
                            * ao_dens(iket,jcontr,ibra,icontr)
            end do
          end do
        end do
      end do
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "trace_sgto_ints", STDOUT)
#endif
    return
  end subroutine trace_sgto_ints

  !> \brief reorders integrals of contracted Cartesian Gaussians, and computes
  !>        their traces with atomic-orbital density matrix at the same time
  !> \details \f$tr(AD) = \sum_{ij}A_{ij}D_{ji}\f$
  !> \author Bin Gao
  !> \date 2010-10-07
  !> \param ang_bra is the orbital quantum number (or angular number) of bra center
  !> \param num_cgto_bra is the number of basis functions on bra center
  !> \param power_bra contains the Cartesian powers of basis functions on bra center
  !> \param ang_ket is the orbital quantum number (or angular number) of ket center
  !> \param num_cgto_ket is the number of basis functions on ket center
  !> \param power_ket contains the Cartesian powers of basis functions on ket center
  !> \param num_contr_bra is the number of contractions of bra center
  !> \param num_contr_ket is the number of contractions of ket center
  !> \param num_opt is the number of operators
  !> \param gen_ints contains the contracted integrals from Gen1Int
  !> \param ao_dens contains the corresponding blocks of the atomic-orbital density matrix
  !> \return tr_ints contains the traces of integrals with atomic-orbital density matrix
  subroutine trace_cgto_ints(ang_bra, num_cgto_bra, power_bra, &
                             ang_ket, num_cgto_ket, power_ket, &
                             num_contr_bra, num_contr_ket,     &
                             num_opt, gen_ints, ao_dens, tr_ints)
    use xkind
    implicit none
    integer, intent(in) :: ang_bra
    integer, intent(in) :: num_cgto_bra
    integer, intent(in) :: power_bra(3,num_cgto_bra)
    integer, intent(in) :: ang_ket
    integer, intent(in) :: num_cgto_ket
    integer, intent(in) :: power_ket(3,num_cgto_ket)
    integer, intent(in) :: num_contr_bra
    integer, intent(in) :: num_contr_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(in) :: gen_ints(num_cgto_bra,num_contr_bra, &
                                        num_cgto_ket,num_contr_ket,num_opt)
    real(REALK), intent(in) :: ao_dens(num_cgto_ket,num_contr_ket, &
                                       num_cgto_bra,num_contr_bra)
    real(REALK), intent(inout) :: tr_ints(num_opt)
!f2py intent(in) :: ang_bra
!f2py intent(hide) :: num_cgto_bra
!f2py intent(in) :: power_bra
!f2py intent(in) :: ang_ket
!f2py intent(hide) :: num_cgto_ket
!f2py intent(in) :: power_ket
!f2py intent(hide) :: num_contr_bra
!f2py intent(hide) :: num_contr_ket
!f2py intent(hide) :: num_opt
!f2py intent(in) :: gen_ints
!f2py depend(num_cgto_bra) :: gen_ints
!f2py depend(num_cgto_ket) :: gen_ints
!f2py depend(num_cgto_ket) :: ao_dens
!f2py depend(num_contr_ket) :: ao_dens
!f2py depend(num_cgto_bra) :: ao_dens
!f2py depend(num_contr_bra) :: ao_dens
!f2py intent(inout) :: tr_ints
!f2py depend(num_opt) :: tr_ints
    integer dang_bra_tri, dang_ket_tri        !twice of the angular numbers on bra and ket +3
    integer ibra, iket, icontr, jcontr, iopt  !increment recorders
    integer addr_bra, addr_ket                !addresses of CGTOs on bra and ket centers in Gen1Int
#if defined(XTIME)
    real(REALK) curr_time                     !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! index of x^{l}y^{m}z^{n} with l+m+n=angm is 1+m+(2*angm+3-n)*n/2
    dang_bra_tri = ang_bra+ang_bra+3
    dang_ket_tri = ang_ket+ang_ket+3
    ! loops over operators
    do iopt = 1, num_opt
      ! loops over contractions on ket center
      do jcontr = 1, num_contr_ket
        ! loops over CGTOs on ket center
        do iket = 1, num_cgto_ket
          addr_ket = 1+power_ket(2,iket)+(dang_ket_tri-power_ket(3,iket))*power_ket(3,iket)/2
          ! loops over contractions on bra center
          do icontr = 1, num_contr_bra
            ! loops over CGTOs on bra center
            do ibra = 1, num_cgto_bra
              addr_bra = 1+power_bra(2,ibra)+(dang_bra_tri-power_bra(3,ibra))*power_bra(3,ibra)/2
              tr_ints(iopt) = tr_ints(iopt)                                  &
                            + gen_ints(addr_bra,icontr,addr_ket,jcontr,iopt) &
                            * ao_dens(iket,jcontr,ibra,icontr)
            end do
          end do
        end do
      end do
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "trace_cgto_ints", STDOUT)
#endif
    return
  end subroutine trace_cgto_ints

  !> \brief computes the traces of integrals of contracted Gaussians with
  !>        atomic-orbital density matrix
  !> \details \f$tr(AD) = \sum_{ij}A_{ij}D_{ji}\f$
  !> \author Bin Gao
  !> \date 2010-10-07
  !> \param num_gto_bra is the number of basis functions on bra center
  !> \param num_contr_bra is the number of contractions of bra center
  !> \param num_gto_ket is the number of basis functions on ket center
  !> \param num_contr_ket is the number of contractions of ket center
  !> \param num_opt is the number of operators
  !> \param gen_ints contains the contracted integrals from Gen1Int
  !> \param ao_dens contains the corresponding blocks of the atomic-orbital density matrix
  !> \return tr_ints contains the traces of integrals with atomic-orbital density matrix
  subroutine trace_gto_ints(num_gto_bra, num_contr_bra, num_gto_ket, num_contr_ket, &
                            num_opt, gen_ints, ao_dens, tr_ints)
    use xkind
    implicit none
    integer, intent(in) :: num_gto_bra
    integer, intent(in) :: num_contr_bra
    integer, intent(in) :: num_gto_ket
    integer, intent(in) :: num_contr_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(in) :: gen_ints(num_gto_bra,num_contr_bra, &
                                        num_gto_ket,num_contr_ket,num_opt)
    real(REALK), intent(in) :: ao_dens(num_gto_ket,num_contr_ket, &
                                       num_gto_bra,num_contr_bra)
    real(REALK), intent(inout) :: tr_ints(num_opt)
!f2py intent(hide) :: num_gto_bra
!f2py intent(hide) :: num_contr_bra
!f2py intent(hide) :: num_gto_ket
!f2py intent(hide) :: num_contr_ket
!f2py intent(hide) :: num_opt
!f2py intent(in) :: gen_ints
!f2py intent(in) :: ao_dens
!f2py depend(num_gto_ket) :: ao_dens
!f2py depend(num_contr_ket) :: ao_dens
!f2py depend(num_gto_bra) :: ao_dens
!f2py depend(num_contr_bra) :: ao_dens
!f2py intent(inout) :: tr_ints
!f2py depend(num_opt) :: tr_ints
    integer ibra, iket, icontr, jcontr, iopt  !increment recorders
#if defined(XTIME)
    real(REALK) curr_time                     !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! loops over operators
    do iopt = 1, num_opt
      ! loops over contractions on ket center
      do jcontr = 1, num_contr_ket
        ! loops over GTOs on ket center
        do iket = 1, num_gto_ket
          ! loops over contractions on bra center
          do icontr = 1, num_contr_bra
            ! loops over GTOs on bra center
            do ibra = 1, num_gto_bra
              tr_ints(iopt) = tr_ints(iopt)                          &
                            + gen_ints(ibra,icontr,iket,jcontr,iopt) &
                            * ao_dens(iket,jcontr,ibra,icontr)
            end do
          end do
        end do
      end do
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "trace_gto_ints", STDOUT)
#endif
    return
  end subroutine trace_gto_ints
