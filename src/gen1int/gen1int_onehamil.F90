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
!!  This file provides the Fortran 90 module of one-electron Hamiltonian.
!!
!!  2012-03-20, Bin Gao:
!!  * first version

#include "stdout.h"
#include "kind_matrix.h"

!> \brief Fortran 90 module of one-electron Hamiltonian
!> \author Bin Gao
!> \date 2012-03-20
module gen1int_onehamil

  use xkind
  use london_ao
  use gen1int_geom
  use gen1int_carmom
  use gen1int_nucpot
  implicit none

  ! one-electron Hamiltonian
  type, public :: one_hamil_t
    private
    ! number of property integral matrices
    integer :: num_prop = 1
    ! symmetry of property integral matrices
    integer :: prop_sym = SYMMETRIC_MATRIX
    ! kinetic energy integrals
    type(kin_energy_t) kin_energy
    ! one-electron potential energy integrals
    type(pot_energy_t) pot_energy
  end type one_hamil_t

  public :: OneHamilCreate
  public :: OneHamilView
  public :: OneHamilGetNumProp
  public :: OneHamilGetSymmetry
  public :: OneHamilGetIntegral
  public :: OneHamilDestroy

  contains

  !> \brief initializes the information of one-electron Hamiltonian integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param idx_nuclei contains the atomic centers of nuclei (<1 for non-atomic center)
  !> \param coord_nuclei contains the coordinates of nuclei
  !> \param charge_nuclei contains the charges of nuclei
  !> \return one_hamil contains the information of one-electron Hamiltonian integrals
  !> \return info_prop (==ERR_INFO) indicates the one-electron Hamiltonian integrals
  !>         are not successfully created
  subroutine OneHamilCreate(one_hamil, info_prop, idx_nuclei, coord_nuclei, &
                            charge_nuclei)
    type(one_hamil_t), intent(inout) :: one_hamil
    integer, intent(out) :: info_prop
    integer, optional, intent(in) :: idx_nuclei(:)
    real(REALK), optional, intent(in) :: coord_nuclei(:,:)
    real(REALK), optional, intent(in) :: charge_nuclei(:)
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    call KinEnergyCreate(kin_energy=one_hamil%kin_energy, &
                         info_prop=info_prop)
    call PotEnergyCreate(pot_energy=one_hamil%pot_energy, &
                         info_prop=info_prop,             &
                         idx_nuclei=idx_nuclei,           &
                         coord_nuclei=coord_nuclei,       &
                         charge_nuclei=charge_nuclei)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OneHamilCreate", STDOUT)
#endif
  end subroutine OneHamilCreate

  !> \brief visualizes the information of one-electron Hamiltonian integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param one_hamil contains the information of one-electron Hamiltonian integrals
  !> \param io_viewer is the logical unit number of the viewer
  subroutine OneHamilView(one_hamil, io_viewer)
    type(one_hamil_t), intent(in) :: one_hamil
    integer, intent(in) :: io_viewer
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    write(io_viewer,100) "one-electron Hamiltonian integrals"
    write(io_viewer,100) "number of property integral matrices", one_hamil%num_prop
    write(io_viewer,100) "symmetry of property integral matrices", one_hamil%prop_sym
    call KinEnergyView(one_hamil%kin_energy, io_viewer)
    call PotEnergyView(one_hamil%pot_energy, io_viewer)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OneHamilView", STDOUT)
#endif
100 format("OneHamilView>> ",A,I6)
  end subroutine OneHamilView

  !> \brief returns the number of integral matrices for given one-electron Hamiltonian integrals
  !> \param one_hamil contains the information of one-electron Hamiltonian integrals
  !> \return num_prop is the number of property integral matrices
  subroutine OneHamilGetNumProp(one_hamil, num_prop)
    type(one_hamil_t), intent(in) :: one_hamil
    integer, intent(out) :: num_prop
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif 
    num_prop = one_hamil%num_prop
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OneHamilGetNumProp", STDOUT)
#endif
  end subroutine OneHamilGetNumProp

  !> \brief returns the symmetry of integral matrices for given one-electron Hamiltonian integrals
  !> \param one_hamil contains the information of one-electron Hamiltonian integrals
  !> \return prop_sym indicates the symmetry of property integral matrices
  subroutine OneHamilGetSymmetry(one_hamil, prop_sym)
    type(one_hamil_t), intent(in) :: one_hamil
    integer, intent(out) :: prop_sym
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif 
    prop_sym = one_hamil%prop_sym
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OneHamilGetSymmetry", STDOUT)
#endif
  end subroutine OneHamilGetSymmetry

  !> \brief evaluates the one-electron Hamiltonian integrals
  !> \author Bin Gao
  !> \date 2011-12-28
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
  !> \param one_hamil contains the information of one-electron Hamiltonian integrals
  !> \param order_mag_bra is the order of magnetic derivatives on bra center
  !> \param order_mag_ket is the order of magnetic derivatives on ket center
  !> \param order_mag_total is the order of total magnetic derivatives
  !> \param order_ram_bra is the order of derivatives w.r.t. total rotational angular momentum on bra center
  !> \param order_ram_ket is the order of derivatives w.r.t. total rotational angular momentum on ket center
  !> \param order_ram_total is the order of total derivatives w.r.t. total rotational angular momentum
  !> \param order_geo_bra is the order of geometric derivatives with respect to bra center
  !> \param order_geo_ket is the order of geometric derivatives with respect to ket center
  !> \param nary_tree_total contains the information of N-ary tree for total geometric derivatives
  !> \param num_gto_bra is the number of spherical/Cartesian GTOs on bra center
  !> \param num_gto_ket is the number of spherical/Cartesian GTOs on ket center
  !> \param num_opt is the number of operators including derivatives
  !> \param mag_num_bra contains the magnetic numbers of spherical GTOs on bra center
  !> \param mag_num_ket contains the magnetic numbers of spherical GTOs on ket center
  !> \param powers_bra contains the Cartesian powers of Cartesian GTOs on bra center
  !> \param powers_ket contains the Cartesian powers of Cartesian GTOs on ket center
  !> \return contr_ints contains the calculated contracted integrals
  subroutine OneHamilGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                 exponent_bra, num_contr_bra, contr_coef_bra,   &
                                 idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                 exponent_ket, num_contr_ket, contr_coef_ket,   &
                                 spher_gto, info_LAO, one_hamil,                &
                                 order_mag_bra, order_mag_ket, order_mag_total, &
                                 order_ram_bra, order_ram_ket, order_ram_total, &
                                 order_geo_bra, order_geo_ket, nary_tree_total, &
                                 num_gto_bra, num_gto_ket, num_opt, contr_ints, &
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
    type(one_hamil_t), intent(in) :: one_hamil
    integer, intent(in) :: order_mag_bra
    integer, intent(in) :: order_mag_ket
    integer, intent(in) :: order_mag_total
    integer, intent(in) :: order_ram_bra
    integer, intent(in) :: order_ram_ket
    integer, intent(in) :: order_ram_total
    integer, optional, intent(in) :: order_geo_bra
    integer, optional, intent(in) :: order_geo_ket
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
    integer num_diff_cent                            !number of differentiated centers of current path
    real(REALK), allocatable :: tmp_ints(:,:,:,:,:)  !temporary integrals
    integer ierr                                     !error information
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (present(nary_tree_total)) then
      call NaryTreePathGetNumCenters(nary_tree_total, num_diff_cent)
      if (num_diff_cent>3) then
        contr_ints = 0.0_REALK
        return
      else if (num_diff_cent==3) then
        call PotEnergyGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                  exponent_bra, num_contr_bra, contr_coef_bra,   &
                                  idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                  exponent_ket, num_contr_ket, contr_coef_ket,   &
                                  spher_gto, info_LAO, one_hamil%pot_energy,     &
                                  order_mag_bra, order_mag_ket, order_mag_total, &
                                  order_ram_bra, order_ram_ket, order_ram_total, &
                                  order_geo_bra, order_geo_ket, nary_tree_total, &
                                  num_gto_bra, num_gto_ket, num_opt, contr_ints, &
                                  mag_num_bra, mag_num_ket, powers_bra, powers_ket)
        return
      end if
    end if
    allocate(tmp_ints(num_gto_bra,num_contr_bra,num_gto_ket,num_contr_ket, &
                      num_opt), stat=ierr)
    if (ierr/=0)                                                            &
      call error_stop("OneHamilGetIntegral", "failed to allocate tmp_ints", &
                      num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket*num_opt)
    call KinEnergyGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                              exponent_bra, num_contr_bra, contr_coef_bra,   &
                              idx_ket, coord_ket, angular_ket, num_prim_ket, &
                              exponent_ket, num_contr_ket, contr_coef_ket,   &
                              spher_gto, info_LAO, one_hamil%kin_energy,     &
                              order_mag_bra, order_mag_ket, order_mag_total, &
                              order_ram_bra, order_ram_ket, order_ram_total, &
                              order_geo_bra, order_geo_ket, nary_tree_total, &
                              num_gto_bra, num_gto_ket, num_opt, tmp_ints,   &
                              mag_num_bra, mag_num_ket, powers_bra, powers_ket)
    call PotEnergyGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                              exponent_bra, num_contr_bra, contr_coef_bra,   &
                              idx_ket, coord_ket, angular_ket, num_prim_ket, &
                              exponent_ket, num_contr_ket, contr_coef_ket,   &
                              spher_gto, info_LAO, one_hamil%pot_energy,     &
                              order_mag_bra, order_mag_ket, order_mag_total, &
                              order_ram_bra, order_ram_ket, order_ram_total, &
                              order_geo_bra, order_geo_ket, nary_tree_total, &
                              num_gto_bra, num_gto_ket, num_opt, contr_ints, &
                              mag_num_bra, mag_num_ket, powers_bra, powers_ket)
    contr_ints = contr_ints+tmp_ints
    deallocate(tmp_ints)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OneHamilGetIntegral", STDOUT)
#endif
  end subroutine OneHamilGetIntegral

  !> \brief frees space taken by the one-electron Hamiltonian integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param one_hamil contains the information of one-electron Hamiltonian integrals
  subroutine OneHamilDestroy(one_hamil)
    type(one_hamil_t), intent(inout) :: one_hamil
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    call KinEnergyDestroy(one_hamil%kin_energy)
    call PotEnergyDestroy(one_hamil%pot_energy)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OneHamilDestroy", STDOUT)
#endif
  end subroutine OneHamilDestroy

end module gen1int_onehamil
