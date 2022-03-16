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
!!  This file provides the Fortran 90 module of nuclear attraction potential integrals.
!!
!!  2012-03-20, Bin Gao:
!!  * first version

#include "stdout.h"
#include "err_info.h"
#include "kind_matrix.h"
#include "max_idx_non.h"

!> \brief Fortran 90 module of nuclear attraction potential integrals
!> \author Bin Gao
!> \date 2012-03-20
module gen1int_nucpot

  use xkind
  use london_ao
  use gen1int_geom
  implicit none

  ! one-electron potential energy integrals
  type, public :: pot_energy_t
    private
    ! number of property integral matrices
    integer :: num_prop = 1
    ! symmetry of property integral matrices
    integer :: prop_sym = SYMMETRIC_MATRIX
    ! number of nuclei
    integer :: num_nuclei = 0
    ! atomic centers of nuclei (<1 for non-atomic center)
    integer, allocatable :: idx_nuclei(:)
    ! coordinates of nuclei
    real(REALK), allocatable :: coord_nuclei(:,:)
    ! charges of nuclei
    real(REALK), allocatable :: charge_nuclei(:)
    ! order of geometric derivatives with respect to the nuclei
    integer :: order_geo_pot = 0
    ! order of electronic derivatives
    integer :: order_elec = 0
  end type pot_energy_t

  ! PSO integrals
  type, public :: pso_t
    private
    ! number of property integral matrices
    integer :: num_prop = 0
    ! symmetry of property integral matrices
    integer :: prop_sym = ANTI_SYM_MATRIX
    ! number of nuclei
    integer :: num_nuclei = 0
    ! atomic centers of nuclei (<1 for non-atomic center)
    integer, allocatable :: idx_nuclei(:)
    ! coordinates of nuclei
    real(REALK), allocatable :: coord_nuclei(:,:)
  end type pso_t

  public :: PotEnergyCreate
  public :: PotEnergyView
  public :: PotEnergyGetNumProp
  public :: PotEnergyGetSymmetry
  public :: PotEnergyGetIntegral
  public :: PotEnergyDestroy

  public :: PSOCreate
  public :: PSOView
  public :: PSOGetNumProp
  public :: PSOGetSymmetry
  public :: PSOGetIntegral
  public :: PSODestroy

  contains

  !> \brief initializes the information of one-electron potential energy integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param idx_nuclei contains the atomic centers of nuclei (<1 for non-atomic center)
  !> \param coord_nuclei contains the coordinates of nuclei
  !> \param charge_nuclei contains the charges of nuclei
  !> \param order_geo_pot is the orders of geometric derivatives with respect to the nuclei
  !> \param order_elec is the electronic derivatives
  !> \return pot_energy contains the information of one-electron potential energy integrals
  !> \return info_prop (==ERR_INFO) indicates the one-electron potential energy integrals
  !>         are not successfully created
  subroutine PotEnergyCreate(pot_energy, info_prop, idx_nuclei, coord_nuclei, &
                             charge_nuclei, order_geo_pot, order_elec)
    type(pot_energy_t), intent(inout) :: pot_energy
    integer, intent(out) :: info_prop
    integer, optional, intent(in) :: idx_nuclei(:)
    real(REALK), optional, intent(in) :: coord_nuclei(:,:)
    real(REALK), optional, intent(in) :: charge_nuclei(:)
    integer, optional, intent(in) :: order_geo_pot
    integer, optional, intent(in) :: order_elec
    integer icent          !incremental recorder
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (present(order_geo_pot)) pot_energy%order_geo_pot = order_geo_pot
    if (present(order_elec)) pot_energy%order_elec = max(0,order_elec)
    if (mod(pot_energy%order_elec,2)==1) then
      pot_energy%prop_sym = ANTI_SYM_MATRIX
    end if
    ! negative order, returns with error code
    if (pot_energy%order_geo_pot<0) then
      info_prop = pot_energy%order_geo_pot
    else if (present(coord_nuclei) .and. present(charge_nuclei)) then
      pot_energy%num_nuclei = size(coord_nuclei,2)
      if (size(coord_nuclei,1)==3 .and. &
          pot_energy%num_nuclei==size(charge_nuclei)) then
        info_prop = 0
        allocate(pot_energy%idx_nuclei(pot_energy%num_nuclei), stat=info_prop)
        if (info_prop/=0) return
        allocate(pot_energy%coord_nuclei(3,pot_energy%num_nuclei), stat=info_prop)
        if (info_prop/=0) then
          deallocate(pot_energy%idx_nuclei)
          return
        end if
        allocate(pot_energy%charge_nuclei(pot_energy%num_nuclei), stat=info_prop)
        if (info_prop/=0) then
          deallocate(pot_energy%idx_nuclei)
          deallocate(pot_energy%coord_nuclei)
          return
        end if
        if (present(idx_nuclei)) then
          if (size(idx_nuclei)/=pot_energy%num_nuclei) then
            deallocate(pot_energy%idx_nuclei)
            deallocate(pot_energy%coord_nuclei)
            deallocate(pot_energy%charge_nuclei)
            info_prop = ERR_INFO
            return
          end if
          pot_energy%idx_nuclei = idx_nuclei
        else
          do icent = 1, pot_energy%num_nuclei
            pot_energy%idx_nuclei(icent) = icent
          end do
        end if
        pot_energy%coord_nuclei = coord_nuclei
        pot_energy%charge_nuclei = charge_nuclei
        pot_energy%num_prop = (pot_energy%order_geo_pot+1) &
                            * (pot_energy%order_geo_pot+2) &
                            * (pot_energy%order_elec+1)    &
                            * (pot_energy%order_elec+2)/4
      else
        info_prop = ERR_INFO
      end if
    else
      info_prop = ERR_INFO
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "PotEnergyCreate", STDOUT)
#endif
  end subroutine PotEnergyCreate

  !> \brief visualizes the information of one-electron potential energy integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param pot_energy contains the information of one-electron potential energy integrals
  !> \param io_viewer is the logical unit number of the viewer
  subroutine PotEnergyView(pot_energy, io_viewer)
    type(pot_energy_t), intent(in) :: pot_energy
    integer, intent(in) :: io_viewer
    integer icent          !incremental recorder over centers
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    write(io_viewer,100) "one-electron potential energy integrals"
    write(io_viewer,100) "number of property integral matrices", pot_energy%num_prop
    write(io_viewer,100) "symmetry of property integral matrices", pot_energy%prop_sym
    write(io_viewer,100) "number of nuclei", pot_energy%num_nuclei
    write(io_viewer,100) "order of geometric derivatives with respect to nuclei", &
                         pot_energy%order_geo_pot
    write(io_viewer,100) "order of electronic derivatives", pot_energy%order_elec
    write(io_viewer,110) "atomic center", "coordinates", "charges"
    do icent = 1, pot_energy%num_nuclei
      write(io_viewer,120) pot_energy%idx_nuclei(icent),     &
                           pot_energy%coord_nuclei(:,icent), &
                           pot_energy%charge_nuclei(icent)
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "PotEnergyView", STDOUT)
#endif
100 format("PotEnergyView>> ",A,I6)
110 format("PotEnergyView>> ",A,10X,A,15X,A)
120 format("PotEnergyView>> ",I6,3F14.8,Es14.6)
  end subroutine PotEnergyView

  !> \brief returns the number of integral matrices for given one-electron potential energy integrals
  !> \param pot_energy contains the information of one-electron potential integrals
  !> \return num_prop is the number of property integral matrices
  subroutine PotEnergyGetNumProp(pot_energy, num_prop)
    type(pot_energy_t), intent(in) :: pot_energy
    integer, intent(out) :: num_prop
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif 
    num_prop = pot_energy%num_prop
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "PotEnergyGetNumProp", STDOUT)
#endif
  end subroutine PotEnergyGetNumProp

  !> \brief returns the symmetry of integral matrices for given one-electron potential energy integrals
  !> \param pot_energy contains the information of one-electron potential integrals
  !> \return prop_sym indicates the symmetry of property integral matrices
  subroutine PotEnergyGetSymmetry(pot_energy,  prop_sym)
    type(pot_energy_t), intent(in) :: pot_energy
    integer, intent(out) :: prop_sym
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif 
    prop_sym = pot_energy%prop_sym
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "PotEnergyGetSymmetry", STDOUT)
#endif
  end subroutine PotEnergyGetSymmetry

  !> \brief evaluates the one-electron potential energy integrals
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
  !> \param pot_energy contains the information of one-electron potential energy integrals
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
  subroutine PotEnergyGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                  exponent_bra, num_contr_bra, contr_coef_bra,   &
                                  idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                  exponent_ket, num_contr_ket, contr_coef_ket,   &
                                  spher_gto, info_LAO, pot_energy,               &
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
    type(pot_energy_t), intent(in) :: pot_energy
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
    integer icent                                    !incremental recorder over centers
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
      end if
    end if
    ! the first atomic center
    call IntGetNUCPOT(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                      exponent_bra, num_contr_bra, contr_coef_bra,   &
                      idx_ket, coord_ket, angular_ket, num_prim_ket, &
                      exponent_ket, num_contr_ket, contr_coef_ket,   &
                      spher_gto, info_LAO, pot_energy%order_elec,    &
                      pot_energy%idx_nuclei(1),                      &
                      pot_energy%coord_nuclei(:,1), MAX_IDX_NON,     &
                      (/0.0_REALK,0.0_REALK,0.0_REALK/),             &
                      pot_energy%charge_nuclei(1), 0,                &
                      order_mag_bra, order_mag_ket, order_mag_total, &
                      order_ram_bra, order_ram_ket, order_ram_total, &
                      order_geo_bra, order_geo_ket,                  &
                      pot_energy%order_geo_pot, 0, nary_tree_total,  &
                      num_gto_bra, num_gto_ket, num_opt, contr_ints, &
                      mag_num_bra, mag_num_ket, powers_bra, powers_ket)
    allocate(tmp_ints(num_gto_bra,num_contr_bra,num_gto_ket,num_contr_ket, &
                      num_opt), stat=ierr)
    if (ierr/=0)                                                             &
      call error_stop("PotEnergyGetIntegral", "failed to allocate tmp_ints", &
                      num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket*num_opt)
    ! other atomic centers
    do icent = 2, pot_energy%num_nuclei
      call IntGetNUCPOT(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                        exponent_bra, num_contr_bra, contr_coef_bra,   &
                        idx_ket, coord_ket, angular_ket, num_prim_ket, &
                        exponent_ket, num_contr_ket, contr_coef_ket,   &
                        spher_gto, info_LAO, pot_energy%order_elec,    &
                        pot_energy%idx_nuclei(icent),                  &
                        pot_energy%coord_nuclei(:,icent), MAX_IDX_NON, &
                        (/0.0_REALK,0.0_REALK,0.0_REALK/),             &
                        pot_energy%charge_nuclei(icent), 0,            &
                        order_mag_bra, order_mag_ket, order_mag_total, &
                        order_ram_bra, order_ram_ket, order_ram_total, &
                        order_geo_bra, order_geo_ket,                  &
                        pot_energy%order_geo_pot, 0, nary_tree_total,  &
                        num_gto_bra, num_gto_ket, num_opt, tmp_ints,   &
                        mag_num_bra, mag_num_ket, powers_bra, powers_ket)
      contr_ints = contr_ints+tmp_ints
    end do
    deallocate(tmp_ints)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "PotEnergyGetIntegral", STDOUT)
#endif
  end subroutine PotEnergyGetIntegral

  !> \brief frees space taken by the one-electron potential energy integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param pot_energy contains the information of one-electron potential energy integrals
  subroutine PotEnergyDestroy(pot_energy)
    type(pot_energy_t), intent(inout) :: pot_energy
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    pot_energy%num_nuclei = 0
    deallocate(pot_energy%idx_nuclei)
    deallocate(pot_energy%coord_nuclei)
    deallocate(pot_energy%charge_nuclei)
    pot_energy%order_geo_pot = 0
    pot_energy%order_elec = 0
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "PotEnergyDestroy", STDOUT)
#endif
  end subroutine PotEnergyDestroy

  !> \brief initializes the information of PSO integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param idx_nuclei contains the atomic centers of nuclei (<1 for non-atomic center)
  !> \param coord_nuclei contains the coordinates of nuclei
  !> \return pso contains the information of PSO integrals
  !> \return info_prop (==ERR_INFO) indicates the PSO integrals
  !>         are not successfully created
  subroutine PSOCreate(pso, info_prop, idx_nuclei, coord_nuclei)
    type(pso_t), intent(inout) :: pso
    integer, intent(out) :: info_prop
    integer, intent(in) :: idx_nuclei(:)
    real(REALK), intent(in) :: coord_nuclei(:,:)
    integer icent          !incremental recorder
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    pso%num_nuclei = size(idx_nuclei)
    if (pso%num_nuclei/=size(coord_nuclei,2)) then
      info_prop = size(coord_nuclei,2)-pso%num_nuclei
      return
    else if (size(coord_nuclei,1)/=3) then
      info_prop = size(coord_nuclei,1)
      return
    end if
    allocate(pso%idx_nuclei(pso%num_nuclei), stat=info_prop)
    if (info_prop/=0) return
    allocate(pso%coord_nuclei(3,pso%num_nuclei), stat=info_prop)
    if (info_prop/=0) then
      deallocate(pso%idx_nuclei)
      return
    end if
    pso%idx_nuclei = idx_nuclei
    pso%coord_nuclei = coord_nuclei
    pso%num_prop = 3*pso%num_nuclei
    info_prop = 0 
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "PSOCreate", STDOUT)
#endif
  end subroutine PSOCreate

  !> \brief visualizes the information of PSO integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param pso contains the information of PSO integrals
  !> \param io_viewer is the logical unit number of the viewer
  subroutine PSOView(pso, io_viewer)
    type(pso_t), intent(in) :: pso
    integer, intent(in) :: io_viewer
    integer icent          !incremental recorder over centers
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    write(io_viewer,100) "PSO integrals"
    write(io_viewer,100) "number of property integral matrices", pso%num_prop
    write(io_viewer,100) "symmetry of property integral matrices", pso%prop_sym
    write(io_viewer,100) "number of nuclei", pso%num_nuclei
    write(io_viewer,110) "atomic center", "coordinates"
    do icent = 1, pso%num_nuclei
      write(io_viewer,120) pso%idx_nuclei(icent), &
                           pso%coord_nuclei(:,icent)
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "PSOView", STDOUT)
#endif
100 format("PSOView>> ",A,I6)
110 format("PSOView>> ",A,10X,A,15X)
120 format("PSOView>> ",I6,3F14.8)
  end subroutine PSOView

  !> \brief returns the number of integral matrices for given PSO integrals
  !> \param pso contains the information of one-electron potential integrals
  !> \return num_prop is the number of property integral matrices
  subroutine PSOGetNumProp(pso, num_prop)
    type(pso_t), intent(in) :: pso
    integer, intent(out) :: num_prop
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif 
    num_prop = pso%num_prop
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "PSOGetNumProp", STDOUT)
#endif
  end subroutine PSOGetNumProp

  !> \brief returns the symmetry of integral matrices for given PSO integrals
  !> \param pso contains the information of one-electron potential integrals
  !> \return prop_sym indicates the symmetry of property integral matrices
  subroutine PSOGetSymmetry(pso,  prop_sym)
    type(pso_t), intent(in) :: pso
    integer, intent(out) :: prop_sym
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif 
    prop_sym = pso%prop_sym
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "PSOGetSymmetry", STDOUT)
#endif
  end subroutine PSOGetSymmetry

  !> \brief evaluates the PSO integrals
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
  !> \param pso contains the information of PSO integrals
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
  subroutine PSOGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                            exponent_bra, num_contr_bra, contr_coef_bra,   &
                            idx_ket, coord_ket, angular_ket, num_prim_ket, &
                            exponent_ket, num_contr_ket, contr_coef_ket,   &
                            spher_gto, info_LAO, pso,                      &
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
    type(pso_t), intent(in) :: pso
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
    integer num_diff_cent                              !number of differentiated centers of current path
    real(REALK), allocatable :: tmp_ints(:,:,:,:,:,:)  !temporary integrals
    integer ierr                                       !error information
    integer icent, jcent                               !incremental recorder over centers
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
      end if
    end if
    allocate(tmp_ints(num_gto_bra,num_contr_bra,num_gto_ket,num_contr_ket,3,3), stat=ierr)
    if (ierr/=0)                                                       &
      call error_stop("PSOGetIntegral", "failed to allocate tmp_ints", &
                      num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket*9)
    jcent = 0
    do icent = 1, pso%num_nuclei
      call IntGetNUCPOT(idx_bra, coord_bra, angular_bra, num_prim_bra,       &
                        exponent_bra, num_contr_bra, contr_coef_bra,         &
                        idx_ket, coord_ket, angular_ket, num_prim_ket,       &
                        exponent_ket, num_contr_ket, contr_coef_ket,         &
                        spher_gto, info_LAO, 1,                              &
                        pso%idx_nuclei(icent), pso%coord_nuclei(:,icent),    &
                        MAX_IDX_NON, (/0.0_REALK,0.0_REALK,0.0_REALK/),      &
                        1.0_REALK, 0,                                        &
                        order_mag_bra, order_mag_ket, order_mag_total,       &
                        order_ram_bra, order_ram_ket, order_ram_total,       &
                        order_geo_bra, order_geo_ket, 1, 0, nary_tree_total, &
                        num_gto_bra, num_gto_ket, num_opt, tmp_ints,         &
                        mag_num_bra, mag_num_ket, powers_bra, powers_ket)
      ! operator (y_k*d/dz - z_k*d/dy)/r_k^3
      jcent = jcent+1
      contr_ints(:,:,:,:,jcent) = tmp_ints(:,:,:,:,3,2)-tmp_ints(:,:,:,:,2,3)
      ! operator (z_k*d/dx - x_k*d/dz)/r_k^3
      jcent = jcent+1
      contr_ints(:,:,:,:,jcent) = tmp_ints(:,:,:,:,1,3)-tmp_ints(:,:,:,:,3,1)
      ! operator (x_k*d/dy - y_k*d/dx)/r_k^3
      jcent = jcent+1
      contr_ints(:,:,:,:,jcent) = tmp_ints(:,:,:,:,2,1)-tmp_ints(:,:,:,:,1,2)
    end do
    deallocate(tmp_ints)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "PSOGetIntegral", STDOUT)
#endif
  end subroutine PSOGetIntegral

  !> \brief frees space taken by the PSO integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param pso contains the information of PSO integrals
  subroutine PSODestroy(pso)
    type(pso_t), intent(inout) :: pso
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    pso%num_nuclei = 0
    deallocate(pso%idx_nuclei)
    deallocate(pso%coord_nuclei)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "PSODestroy", STDOUT)
#endif
  end subroutine PSODestroy

end module gen1int_nucpot
