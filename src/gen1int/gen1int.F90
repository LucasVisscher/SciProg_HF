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
!!  This file provides the Fortran 90 module of using Gen1Int library.
!!
!!  2012-02-10, Bin Gao:
!!  * adds overlap distribution
!!
!!  2011-12-11, Bin Gao:
!!  * first version

#include "stdout.h"
#include "err_info.h"
#include "kind_matrix.h"

!> \brief Fortran 90 module of Gen1Int library
!> \author Bin Gao
!> \date 2011-12-11
module gen1int

  use xkind
  use london_ao
  use gen1int_geom
  use gen1int_carmom
  use gen1int_nucpot
  use gen1int_onehamil
  use gen1int_gaupot
  implicit none

  ! label of angular momentum integrals
  character*(*), parameter, public :: INT_ANGMOM = "INT_ANGMOM"
  ! label of overlap integrals
  character*(*), parameter, public :: INT_OVERLAP = "INT_OVERLAP"
  ! label of kinetic energy integrals
  character*(*), parameter, public :: INT_KIN_ENERGY = "INT_KIN_ENERGY"
  ! label of one-electron potential energy integrals
  character*(*), parameter, public :: INT_POT_ENERGY = "INT_POT_ENERGY"
  ! label of one-electron Hamiltonian
  character*(*), parameter, public :: INT_ONE_HAMIL = "INT_ONE_HAMIL"
  ! label of Cartesian multipole integrals
  character*(*), parameter, public :: INT_CART_MULTIPOLE = "INT_CART_MULTIPOLE"
  ! label of spherical multipole integrals
  character*(*), parameter, public :: INT_SPHER_MULTIPOLE = "INT_SPHER_MULTIPOLE"
  ! label of Gaussian charge potential integrals
  character*(*), parameter, public :: INT_GAUSSIAN_POT = "INT_GAUSSIAN_POT"
  ! label of PSO integrals
  character*(*), parameter, public :: INT_PSO = "INT_PSO"

  ! kind of property integral matrices, 1 for symmetric, -1 for anti-symmetric, 0 for square
  integer, parameter, public :: SYMM_INT_MAT = SYMMETRIC_MATRIX
  integer, parameter, public :: ANTI_INT_MAT = ANTI_SYM_MATRIX
  integer, parameter, public :: SQUARE_INT_MAT = SQUARE_MATRIX

  ! one-electron property integrals
  type, public :: one_prop_t
    private
    ! indicates if the property integrals are created
    logical :: prop_init = .false.
    ! angular momentum integrals
    type(ang_mom_t), allocatable :: ang_mom(:)
    ! overlap integrals
    type(overlap_t), allocatable :: overlap(:)
    ! kinetic energy integrals
    type(kin_energy_t), allocatable :: kin_energy(:)
    ! one-electron potential energy integrals
    type(pot_energy_t), allocatable :: pot_energy(:)
    ! one-electron Hamiltonian
    type(one_hamil_t), allocatable :: one_hamil(:)
    ! multipole integrals
    type(multipole_t), allocatable :: multipole(:)
    ! Gaussian charge potential integrals
    type(gaussian_pot_t), allocatable :: gaussian_pot(:)
    ! PSO integrals
    type(pso_t), allocatable :: pso(:)
    ! order of total magnetic derivatives
    integer :: order_mag = 0
    ! order of partial magnetic derivatives on bra center
    integer :: order_mag_bra = 0
    ! order of partial magnetic derivatives on ket center
    integer :: order_mag_ket = 0
    ! order of total derivatives w.r.t. total rotational angular momentum (RAM)
    integer :: order_ram = 0
    ! order of partial derivatives on bra center w.r.t. total RAM
    integer :: order_ram_bra = 0
    ! order of partial derivatives on ket center w.r.t. total RAM
    integer :: order_ram_ket = 0
    ! information of London atomic orbital
    type(london_ao_t) :: info_LAO
  end type one_prop_t

  public :: OnePropCreate
  public :: OnePropSetMag
  public :: OnePropSetLAO
  public :: OnePropSetRAM
  public :: OnePropView
  public :: OnePropGetNumProp
  public :: OnePropGetSymmetry
  public :: OnePropGetIntegral
  public :: OnePropGetFunction
  public :: OnePropDestroy

  contains

  !> \brief initializes the information of one-electron property integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param prop_name is the name of one-electron property integrals
  !> \param idx_nuclei contains the atomic centers of nuclei (<1 for non-atomic center)
  !> \param coord_nuclei contains the coordinates of nuclei
  !> \param charge_nuclei contains the charges of nuclei
  !> \param order_geo_pot is the order of geometric derivatives with respect to the potential center
  !> \param order_elec is the order of electronic derivatives
  !> \param idx_diporg is the atomic center of dipole origin (<1 for non-atomic center)
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param order_mom is the order of multipole integrals
  !> \param idx_gauorg contains the atomic centers of Gaussian charge potential origins
  !>        (<1 for non-atomic center)
  !> \param gaupot_origin contains the coordinates of Gaussian charge potential origins
  !> \param gaupot_charge contains the charges of Gaussian charge potential
  !> \param gaupot_expt contains the exponents used in the Gaussian broadening function of charges
  !> \return one_prop contains the information of one-electron property integrals
  !> \return info_prop (==ERR_INFO) indicates the one-electron property integrals
  !>         are not successfully created
  subroutine OnePropCreate(prop_name, one_prop, info_prop,           &
                           idx_nuclei, coord_nuclei,                 &
                           charge_nuclei, order_geo_pot, order_elec, &
                           idx_diporg, dipole_origin, order_mom,     &
                           idx_gauorg, gaupot_origin,                &
                           gaupot_charge, gaupot_expt)
    character*(*), intent(in) :: prop_name
    type(one_prop_t), intent(inout) :: one_prop
    integer, intent(out) :: info_prop
    integer, optional, intent(in) :: idx_nuclei(:)
    real(REALK), optional, intent(in) :: coord_nuclei(:,:)
    real(REALK), optional, intent(in) :: charge_nuclei(:)
    integer, optional, intent(in) :: order_geo_pot
    integer, optional, intent(in) :: order_elec
    integer, optional, intent(in) :: idx_diporg
    real(REALK), optional, intent(in) :: dipole_origin(3)
    integer, optional, intent(in) :: order_mom
    integer, optional, intent(in) :: idx_gauorg(:)
    real(REALK), optional, intent(in) :: gaupot_origin(:,:)
    real(REALK), optional, intent(in) :: gaupot_charge(:)
    real(REALK), optional, intent(in) :: gaupot_expt(:)
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! frees space if the property integrals are created ago
    if (one_prop%prop_init) call OnePropDestroy(one_prop)
    ! different one-electron property integrals
    select case(trim(prop_name))
    ! angular momentum integrals
    case(INT_ANGMOM)
      allocate(one_prop%ang_mom(1), stat=info_prop)
      if (info_prop==0) then
        call AngMomCreate(ang_mom=one_prop%ang_mom(1), &
                          info_prop=info_prop,         &
                          idx_diporg=idx_diporg,       &
                          dipole_origin=dipole_origin)
        if (info_prop==0) then
          one_prop%prop_init = .true.
        else
          one_prop%prop_init = .false.
          deallocate(one_prop%ang_mom)
        end if
      end if
    ! overlap integrals
    case(INT_OVERLAP)
      allocate(one_prop%overlap(1), stat=info_prop)
      if (info_prop==0) then
        call OverlapCreate(overlap=one_prop%overlap(1), info_prop=info_prop)
        if (info_prop==0) then
          one_prop%prop_init = .true.
        else
          one_prop%prop_init = .false.
          deallocate(one_prop%overlap)
        end if
      end if
    ! kinetic energy integrals
    case(INT_KIN_ENERGY)
      allocate(one_prop%kin_energy(1), stat=info_prop)
      if (info_prop==0) then
        call KinEnergyCreate(kin_energy=one_prop%kin_energy(1), info_prop=info_prop)
        if (info_prop==0) then
          one_prop%prop_init = .true.
        else
          one_prop%prop_init = .false.
          deallocate(one_prop%kin_energy)
        end if
      end if
    ! one-electron potential energy integrals
    case(INT_POT_ENERGY)
      allocate(one_prop%pot_energy(1), stat=info_prop)
      if (info_prop==0) then
        call PotEnergyCreate(pot_energy=one_prop%pot_energy(1), &
                             info_prop=info_prop,               &
                             idx_nuclei=idx_nuclei,             &
                             coord_nuclei=coord_nuclei,         &
                             charge_nuclei=charge_nuclei,       &
                             order_geo_pot=order_geo_pot,       &
                             order_elec=order_elec)
        if (info_prop==0) then
          one_prop%prop_init = .true.
        else
          one_prop%prop_init = .false.
          deallocate(one_prop%pot_energy)
        end if
      end if
    ! one-electron Hamiltonian
    case(INT_ONE_HAMIL)
      allocate(one_prop%one_hamil(1), stat=info_prop)
      if (info_prop==0) then
        call OneHamilCreate(one_hamil=one_prop%one_hamil(1), &
                            info_prop=info_prop,             &
                            idx_nuclei=idx_nuclei,           &
                            coord_nuclei=coord_nuclei,       &
                            charge_nuclei=charge_nuclei)
        if (info_prop==0) then
          one_prop%prop_init = .true.
        else
          one_prop%prop_init = .false.
          deallocate(one_prop%one_hamil)
        end if
      end if
    ! Cartesian multipole integrals
    case(INT_CART_MULTIPOLE)
      allocate(one_prop%multipole(1), stat=info_prop)
      if (info_prop==0) then
        call MultipoleCreate(multipole=one_prop%multipole(1), &
                             info_prop=info_prop,             &
                             idx_diporg=idx_diporg,           &
                             dipole_origin=dipole_origin,     &
                             order_mom=order_mom,             &
                             spher_mom=.false.,               &
                             order_elec=order_elec)
        if (info_prop==0) then
          one_prop%prop_init = .true.
        else
          one_prop%prop_init = .false.
          deallocate(one_prop%multipole)
        end if
      end if
    ! spherical multipole integrals
    case(INT_SPHER_MULTIPOLE)
      allocate(one_prop%multipole(1), stat=info_prop)
      if (info_prop==0) then
        call MultipoleCreate(multipole=one_prop%multipole(1), &
                             info_prop=info_prop,             &
                             idx_diporg=idx_diporg,           &
                             dipole_origin=dipole_origin,     &
                             order_mom=order_mom,             &
                             spher_mom=.true.,                &
                             order_elec=order_elec)
        if (info_prop==0) then
          one_prop%prop_init = .true.
        else
          one_prop%prop_init = .false.
          deallocate(one_prop%multipole)
        end if
      end if
    ! Gaussian charge potential integrals
    case(INT_GAUSSIAN_POT)
      allocate(one_prop%gaussian_pot(1), stat=info_prop)
      if (info_prop==0) then
        call GaussianPotCreate(gaussian_pot=one_prop%gaussian_pot(1), &
                               info_prop=info_prop,                   &
                               idx_gauorg=idx_gauorg,                 &
                               gaupot_origin=gaupot_origin,           &
                               gaupot_charge=gaupot_charge,           &
                               gaupot_expt=gaupot_expt,               &
                               order_geo_pot=order_geo_pot)
        if (info_prop==0) then
          one_prop%prop_init = .true.
        else
          one_prop%prop_init = .false.
          deallocate(one_prop%gaussian_pot)
        end if
      end if
    case(INT_PSO)
      allocate(one_prop%pso(1), stat=info_prop)
      if (info_prop==0) then
        call PSOCreate(pso=one_prop%pso(1),   &
                       info_prop=info_prop,   &
                       idx_nuclei=idx_nuclei, &
                       coord_nuclei=coord_nuclei)
        if (info_prop==0) then
          one_prop%prop_init = .true.
        else
          one_prop%prop_init = .false.
          deallocate(one_prop%pso)
        end if
      end if
    case default
      info_prop = ERR_INFO
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OnePropCreate", STDOUT)
#endif
  end subroutine OnePropCreate

  !> \brief sets the magnetic derivatives
  !> \author Bin Gao
  !> \date 2012-05-07
  !> \param one_prop contains the information of one-electron property integrals
  !> \param order_mag is the order of total magnetic derivatives
  !> \param order_mag_bra is the order of partial magnetic derivatives on bra center
  !> \param order_mag_ket is the order of partial magnetic derivatives on ket center
  subroutine OnePropSetMag(one_prop, order_mag, order_mag_bra, order_mag_ket)
    type(one_prop_t), intent(inout) :: one_prop
    integer, optional, intent(in) :: order_mag
    integer, optional, intent(in) :: order_mag_bra
    integer, optional, intent(in) :: order_mag_ket
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (present(order_mag)) one_prop%order_mag = max(0,order_mag)
    if (present(order_mag_bra)) one_prop%order_mag_bra = max(0,order_mag_bra)
    if (present(order_mag_ket)) one_prop%order_mag_ket = max(0,order_mag_ket)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OnePropSetMag", STDOUT)
#endif
  end subroutine OnePropSetMag

  !> \brief sets the information of London atomic orbitals (LAOs)
  !> \author Bin Gao
  !> \date 2013-01-31
  !> \param one_prop contains the information of one-electron property integrals
  !> \param gauge_origin contains the coordinates of gauge origin of the magnetic vector potential
  !> \param origin_London_PF contains the coordinates of origin of the London phase factor
  subroutine OnePropSetLAO(one_prop, gauge_origin, origin_London_PF)
    type(one_prop_t), intent(inout) :: one_prop
    real(REALK), optional, intent(in) :: gauge_origin(3)
    real(REALK), optional, intent(in) :: origin_London_PF(3)
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    call LondonAOSet(one_prop%info_LAO, gauge_origin, origin_London_PF)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OnePropSetLAO", STDOUT)
#endif
  end subroutine OnePropSetLAO

  !> \brief sets the derivatives w.r.t. total rotational angular momentum (RAM)
  !> \author Bin Gao
  !> \date 2012-05-07
  !> \param one_prop contains the information of one-electron property integrals
  !> \param order_ram is the order of total derivatives w.r.t. total RAM
  !> \param order_ram_bra is the order of partial derivatives on bra center w.r.t. total RAM
  !> \param order_ram_ket is the order of partial derivatives on ket center w.r.t. total RAM
  subroutine OnePropSetRAM(one_prop, order_ram, order_ram_bra, order_ram_ket)
    type(one_prop_t), intent(inout) :: one_prop
    integer, optional, intent(in) :: order_ram
    integer, optional, intent(in) :: order_ram_bra
    integer, optional, intent(in) :: order_ram_ket
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (present(order_ram)) one_prop%order_ram = max(0,order_ram)
    if (present(order_ram_bra)) one_prop%order_ram_bra = max(0,order_ram_bra)
    if (present(order_ram_ket)) one_prop%order_ram_ket = max(0,order_ram_ket)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OnePropSetRam", STDOUT)
#endif
  end subroutine OnePropSetRAM

  !> \brief visualizes the information of one-electron property integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param one_prop contains the information of one-electron property integrals
  !> \param io_viewer is the logical unit number of the viewer
  subroutine OnePropView(one_prop, io_viewer)
    type(one_prop_t), intent(in) :: one_prop
    integer, intent(in) :: io_viewer
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (allocated(one_prop%ang_mom)) then
      call AngMomView(one_prop%ang_mom(1), io_viewer)
    else if (allocated(one_prop%overlap)) then
      call OverlapView(one_prop%overlap(1), io_viewer)
    else if (allocated(one_prop%kin_energy)) then
      call KinEnergyView(one_prop%kin_energy(1), io_viewer)
    else if (allocated(one_prop%pot_energy)) then
      call PotEnergyView(one_prop%pot_energy(1), io_viewer)
    else if (allocated(one_prop%one_hamil)) then
      call OneHamilView(one_prop%one_hamil(1), io_viewer)
    else if (allocated(one_prop%multipole)) then
      call MultipoleView(one_prop%multipole(1), io_viewer)
    else if (allocated(one_prop%gaussian_pot)) then
      call GaussianPotView(one_prop%gaussian_pot(1), io_viewer)
    else if (allocated(one_prop%pso)) then
      call PSOView(one_prop%pso(1), io_viewer)
    else
      write(io_viewer,100) "unknown one-electron property integrals"
    end if
    write(io_viewer,100) "order of total magnetic derivatives", one_prop%order_mag
    write(io_viewer,100) "order of partial magnetic derivatives on bra center", &
                         one_prop%order_mag_bra
    write(io_viewer,100) "order of partial magnetic derivatives on ket center", &
                         one_prop%order_mag_ket
    write(io_viewer,100)                                                           &
      "order of total derivatives w.r.t. total rotational angular momentum (RAM)", &
      one_prop%order_ram
    write(io_viewer,100) "order of partial derivatives on bra center w.r.t. total RAM", &
                         one_prop%order_ram_bra
    write(io_viewer,100) "order of partial derivatives on ket center w.r.t. total RAM", &
                         one_prop%order_ram_ket
    call LondonAOView(one_prop%info_LAO, io_viewer)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OnePropView", STDOUT)
#endif
100 format("OnePropView>> ",A,I6)
  end subroutine OnePropView

  !> \brief returns the number of property integral matrices for given one-electron property integrals
  !> \author Bin Gao
  !> \date 2012-01-12
  !> \param one_prop contains the information of one-electron property integrals
  !> \return num_prop is the number of property integral matrices
  subroutine OnePropGetNumProp(one_prop, num_prop)
    type(one_prop_t), intent(in) :: one_prop
    integer, intent(out) :: num_prop
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (allocated(one_prop%ang_mom)) then
      call AngMomGetNumProp(one_prop%ang_mom(1), num_prop)
    else if (allocated(one_prop%overlap)) then
      call OverlapGetNumProp(one_prop%overlap(1), num_prop)
    else if (allocated(one_prop%kin_energy)) then
      call KinEnergyGetNumProp(one_prop%kin_energy(1), num_prop)
    else if (allocated(one_prop%pot_energy)) then
      call PotEnergyGetNumProp(one_prop%pot_energy(1), num_prop)
    else if (allocated(one_prop%one_hamil)) then
      call OneHamilGetNumProp(one_prop%one_hamil(1), num_prop)
    else if (allocated(one_prop%multipole)) then
      call MultipoleGetNumProp(one_prop%multipole(1), num_prop)
    else if (allocated(one_prop%gaussian_pot)) then
      call GaussianPotGetNumProp(one_prop%gaussian_pot(1), num_prop)
    else if (allocated(one_prop%pso)) then
      call PSOGetNumProp(one_prop%pso(1), num_prop)
    else
      call error_stop("OnePropGetNumProp", &
                      "unknown one-electron property integrals", -1)
    end if
    ! gets the number of derivatives
    if (one_prop%order_mag>0) &
       num_prop = num_prop*(one_prop%order_mag+1)*(one_prop%order_mag+2)/2
    if (one_prop%order_mag_bra>0) &
      num_prop = num_prop*(one_prop%order_mag_bra+1)*(one_prop%order_mag_bra+2)/2
    if (one_prop%order_mag_ket>0) &
       num_prop = num_prop*(one_prop%order_mag_ket+1)*(one_prop%order_mag_ket+2)/2
    if (one_prop%order_ram>0) &
      num_prop = num_prop*(one_prop%order_ram+1)*(one_prop%order_ram+2)/2
    if (one_prop%order_ram_bra>0) &
      num_prop = num_prop*(one_prop%order_ram_bra+1)*(one_prop%order_ram_bra+2)/2
    if (one_prop%order_ram_ket>0) &
      num_prop = num_prop*(one_prop%order_ram_ket+1)*(one_prop%order_ram_ket+2)/2
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OnePropGetNumProp", STDOUT)
#endif
  end subroutine OnePropGetNumProp

  !> \brief returns the symmetry of property integral matrices for given one-electron property integrals
  !> \author Bin Gao
  !> \date 2012-01-12
  !> \param one_prop contains the information of one-electron property integrals
  !> \return prop_sym indicates the symmetry of property integral matrices (SYMM_INT_MAT,
  !>         ANTI_INT_MAT, or SQUARE_INT_MAT)
  subroutine OnePropGetSymmetry(one_prop, prop_sym)
    type(one_prop_t), intent(in) :: one_prop
    integer, intent(out) :: prop_sym
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! non-symmetric if having partial derivatives
    if (one_prop%order_mag_bra>0 .or. &
        one_prop%order_mag_ket>0 .or. &
        one_prop%order_ram_bra>0 .or. &
        one_prop%order_ram_ket>0) then
      prop_sym = SQUARE_INT_MAT
    else
      if (allocated(one_prop%ang_mom)) then
        call AngMomGetSymmetry(one_prop%ang_mom(1), prop_sym)
      else if (allocated(one_prop%overlap)) then
        call OverlapGetSymmetry(one_prop%overlap(1), prop_sym)
      else if (allocated(one_prop%kin_energy)) then
        call KinEnergyGetSymmetry(one_prop%kin_energy(1), prop_sym)
      else if (allocated(one_prop%pot_energy)) then
        call PotEnergyGetSymmetry(one_prop%pot_energy(1), prop_sym)
      else if (allocated(one_prop%one_hamil)) then
        call OneHamilGetSymmetry(one_prop%one_hamil(1), prop_sym)
      else if (allocated(one_prop%multipole)) then
        call MultipoleGetSymmetry(one_prop%multipole(1), prop_sym)
      else if (allocated(one_prop%gaussian_pot)) then
        call GaussianPotGetSymmetry(one_prop%gaussian_pot(1), prop_sym)
      else if (allocated(one_prop%pso)) then
        call PSOGetSymmetry(one_prop%pso(1), prop_sym)
      else
        call error_stop("OnePropGetSymmetry", &
                        "unknown one-electron property integrals", -1)
      end if
      if (mod(one_prop%order_mag+one_prop%order_ram,2)==1) then
        prop_sym = ANTI_SYM_MATRIX*prop_sym
      end if
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OnePropGetSymmetry", STDOUT)
#endif
  end subroutine OnePropGetSymmetry

  !> \brief evaluates the one-electron property integrals
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
  !> \param one_prop contains the information of one-electron property integrals
  !> \param order_geo_bra is the order of partial geometric derivatives on bra center
  !> \param order_geo_ket is the order of partial geometric derivatives on ket center
  !> \param nary_tree_total contains the information of N-ary tree for total geometric derivatives
  !> \param num_gto_bra is the number of spherical/Cartesian GTOs on bra center
  !> \param num_gto_ket is the number of spherical/Cartesian GTOs on ket center
  !> \param num_opt is the number of operators including different derivatives
  !> \param mag_num_bra contains the magnetic numbers of spherical GTOs on bra center
  !> \param mag_num_ket contains the magnetic numbers of spherical GTOs on ket center
  !> \param powers_bra contains the Cartesian powers of Cartesian GTOs on bra center
  !> \param powers_ket contains the Cartesian powers of Cartesian GTOs on ket center
  !> \return contr_ints contains the calculated contracted integrals
  subroutine OnePropGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                exponent_bra, num_contr_bra, contr_coef_bra,   &
                                idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                exponent_ket, num_contr_ket, contr_coef_ket,   &
                                spher_gto, one_prop, order_geo_bra,            &
                                order_geo_ket, nary_tree_total,                &
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
    type(one_prop_t), intent(in) :: one_prop
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
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (allocated(one_prop%ang_mom)) then
      call AngMomGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra,     &
                             exponent_bra, num_contr_bra, contr_coef_bra,       &
                             idx_ket, coord_ket, angular_ket, num_prim_ket,     &
                             exponent_ket, num_contr_ket, contr_coef_ket,       &
                             spher_gto, one_prop%info_LAO, one_prop%ang_mom(1), &
                             one_prop%order_mag_bra, one_prop%order_mag_ket,    &
                             one_prop%order_mag, one_prop%order_ram_bra,        &
                             one_prop%order_ram_ket, one_prop%order_ram,        &
                             order_geo_bra, order_geo_ket, nary_tree_total,     &
                             num_gto_bra, num_gto_ket, num_opt, contr_ints,     &
                             mag_num_bra, mag_num_ket, powers_bra, powers_ket)
    else if (allocated(one_prop%overlap)) then
      call OverlapGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra,     &
                              exponent_bra, num_contr_bra, contr_coef_bra,       &
                              idx_ket, coord_ket, angular_ket, num_prim_ket,     &
                              exponent_ket, num_contr_ket, contr_coef_ket,       &
                              spher_gto, one_prop%info_LAO, one_prop%overlap(1), &
                              one_prop%order_mag_bra, one_prop%order_mag_ket,    &
                              one_prop%order_mag, one_prop%order_ram_bra,        &
                              one_prop%order_ram_ket, one_prop%order_ram,        &
                              order_geo_bra, order_geo_ket, nary_tree_total,     &
                              num_gto_bra, num_gto_ket, num_opt, contr_ints,     &
                              mag_num_bra, mag_num_ket, powers_bra, powers_ket)
    else if (allocated(one_prop%kin_energy)) then
      call KinEnergyGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra,        &
                                exponent_bra, num_contr_bra, contr_coef_bra,          & 
                                idx_ket, coord_ket, angular_ket, num_prim_ket,        & 
                                exponent_ket, num_contr_ket, contr_coef_ket,          & 
                                spher_gto, one_prop%info_LAO, one_prop%kin_energy(1), &
                                one_prop%order_mag_bra, one_prop%order_mag_ket,       &
                                one_prop%order_mag, one_prop%order_ram_bra,           &
                                one_prop%order_ram_ket, one_prop%order_ram,           &
                                order_geo_bra, order_geo_ket, nary_tree_total,        &
                                num_gto_bra, num_gto_ket, num_opt, contr_ints,        &
                                mag_num_bra, mag_num_ket, powers_bra, powers_ket)
    else if (allocated(one_prop%pot_energy)) then
      call PotEnergyGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra,        &
                                exponent_bra, num_contr_bra, contr_coef_bra,          &
                                idx_ket, coord_ket, angular_ket, num_prim_ket,        &
                                exponent_ket, num_contr_ket, contr_coef_ket,          &
                                spher_gto, one_prop%info_LAO, one_prop%pot_energy(1), &
                                one_prop%order_mag_bra, one_prop%order_mag_ket,       &
                                one_prop%order_mag, one_prop%order_ram_bra,           &
                                one_prop%order_ram_ket, one_prop%order_ram,           &
                                order_geo_bra, order_geo_ket, nary_tree_total,        &
                                num_gto_bra, num_gto_ket, num_opt, contr_ints,        &
                                mag_num_bra, mag_num_ket, powers_bra, powers_ket)
    else if (allocated(one_prop%one_hamil)) then
      call OneHamilGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra,       &
                               exponent_bra, num_contr_bra, contr_coef_bra,         &
                               idx_ket, coord_ket, angular_ket, num_prim_ket,       &
                               exponent_ket, num_contr_ket, contr_coef_ket,         &
                               spher_gto, one_prop%info_LAO, one_prop%one_hamil(1), &
                               one_prop%order_mag_bra, one_prop%order_mag_ket,      &
                               one_prop%order_mag, one_prop%order_ram_bra,          &
                               one_prop%order_ram_ket, one_prop%order_ram,          &
                               order_geo_bra, order_geo_ket, nary_tree_total,       &
                               num_gto_bra, num_gto_ket, num_opt, contr_ints,       &
                               mag_num_bra, mag_num_ket, powers_bra, powers_ket)
    else if (allocated(one_prop%multipole)) then
      call MultipoleGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra,       &
                                exponent_bra, num_contr_bra, contr_coef_bra,         &
                                idx_ket, coord_ket, angular_ket, num_prim_ket,       &
                                exponent_ket, num_contr_ket, contr_coef_ket,         &
                                spher_gto, one_prop%info_LAO, one_prop%multipole(1), &
                                one_prop%order_mag_bra, one_prop%order_mag_ket,      &
                                one_prop%order_mag, one_prop%order_ram_bra,          &
                                one_prop%order_ram_ket, one_prop%order_ram,          &
                                order_geo_bra, order_geo_ket, nary_tree_total,       &
                                num_gto_bra, num_gto_ket, num_opt, contr_ints,       &
                                mag_num_bra, mag_num_ket, powers_bra, powers_ket)
    else if (allocated(one_prop%gaussian_pot)) then
      call GaussianPotGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra,          &
                                  exponent_bra, num_contr_bra, contr_coef_bra,            &
                                  idx_ket, coord_ket, angular_ket, num_prim_ket,          &
                                  exponent_ket, num_contr_ket, contr_coef_ket,            &
                                  spher_gto, one_prop%info_LAO, one_prop%gaussian_pot(1), &
                                  one_prop%order_mag_bra, one_prop%order_mag_ket,         &
                                  one_prop%order_mag, one_prop%order_ram_bra,             &
                                  one_prop%order_ram_ket, one_prop%order_ram,             &
                                  order_geo_bra, order_geo_ket, nary_tree_total,          &
                                  num_gto_bra, num_gto_ket, num_opt, contr_ints,          &
                                  mag_num_bra, mag_num_ket, powers_bra, powers_ket)
    else if (allocated(one_prop%pso)) then
      call PSOGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra,  &
                          exponent_bra, num_contr_bra, contr_coef_bra,    &
                          idx_ket, coord_ket, angular_ket, num_prim_ket,  &
                          exponent_ket, num_contr_ket, contr_coef_ket,    &
                          spher_gto, one_prop%info_LAO, one_prop%pso(1),  &
                          one_prop%order_mag_bra, one_prop%order_mag_ket, &
                          one_prop%order_mag, one_prop%order_ram_bra,     &
                          one_prop%order_ram_ket, one_prop%order_ram,     &
                          order_geo_bra, order_geo_ket, nary_tree_total,  &
                          num_gto_bra, num_gto_ket, num_opt, contr_ints,  &
                          mag_num_bra, mag_num_ket, powers_bra, powers_ket)
    else
      call error_stop("OnePropGetIntegral", "unknown one-electron property integrals", -1)
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OnePropGetIntegral", STDOUT)
#endif
  end subroutine OnePropGetIntegral

  !> \brief evaluates the one-electron property integrand
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
  !> \param one_prop contains the information of one-electron property integrals
  !> \param order_geo_bra is the order of partial geometric derivatives on bra center
  !> \param order_geo_ket is the order of partial geometric derivatives on ket center
  !> \param nary_tree_total contains the information of N-ary tree for total geometric derivatives
  !> \param num_points is the number of grid points
  !> \param grid_points contains the coordinates of grid points
  !> \param num_gto_bra is the number of spherical/Cartesian GTOs on bra center
  !> \param num_gto_ket is the number of spherical/Cartesian GTOs on ket center
  !> \param num_opt is the number of operators including different derivatives
  !> \param mag_num_bra contains the magnetic numbers of spherical GTOs on bra center
  !> \param mag_num_ket contains the magnetic numbers of spherical GTOs on ket center
  !> \param powers_bra contains the Cartesian powers of Cartesian GTOs on bra center
  !> \param powers_ket contains the Cartesian powers of Cartesian GTOs on ket center
  !> \return contr_ints contains the calculated contracted integrals
  subroutine OnePropGetFunction(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                exponent_bra, num_contr_bra, contr_coef_bra,   &
                                idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                exponent_ket, num_contr_ket, contr_coef_ket,   &
                                spher_gto, one_prop, order_geo_bra,            &
                                order_geo_ket, nary_tree_total,                &
                                num_points, grid_points,                       &
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
    type(one_prop_t), intent(in) :: one_prop
    integer, optional, intent(in) :: order_geo_bra
    integer, optional, intent(in) :: order_geo_ket
    type(nary_tree_t), optional, intent(in) :: nary_tree_total
    integer, intent(in) :: num_points
    real(REALK), intent(in) :: grid_points(3,num_points)
    integer, intent(in) :: num_gto_bra
    integer, intent(in) :: num_gto_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(out) :: contr_ints(num_gto_bra,num_contr_bra, &
                                           num_gto_ket,num_contr_ket, &
                                           num_points,num_opt)
    integer, optional, intent(in) :: mag_num_bra(num_gto_bra)
    integer, optional, intent(in) :: mag_num_ket(num_gto_ket)
    integer, optional, intent(in) :: powers_bra(3,num_gto_bra)
    integer, optional, intent(in) :: powers_ket(3,num_gto_ket)
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (allocated(one_prop%ang_mom)) then
      call error_stop("OnePropGetFunction", "ang_mom not implemented", -1)
    else if (allocated(one_prop%overlap)) then
      call OverlapGetFunction(idx_bra, coord_bra, angular_bra, num_prim_bra,     &
                              exponent_bra, num_contr_bra, contr_coef_bra,       &
                              idx_ket, coord_ket, angular_ket, num_prim_ket,     &
                              exponent_ket, num_contr_ket, contr_coef_ket,       &
                              spher_gto, one_prop%info_LAO, one_prop%overlap(1), &
                              one_prop%order_mag_bra, one_prop%order_mag_ket,    &
                              one_prop%order_mag, one_prop%order_ram_bra,        &
                              one_prop%order_ram_ket, one_prop%order_ram,        &
                              order_geo_bra, order_geo_ket, nary_tree_total,     &
                              num_points, grid_points,                           &
                              num_gto_bra, num_gto_ket, num_opt, contr_ints,     &
                              mag_num_bra, mag_num_ket, powers_bra, powers_ket)
    else if (allocated(one_prop%kin_energy)) then
      call error_stop("OnePropGetFunction", "kin_energy not implemented", -1)
    else if (allocated(one_prop%pot_energy)) then
      call error_stop("OnePropGetFunction", "pot_energy not implemented", -1)
    else if (allocated(one_prop%one_hamil)) then
      call error_stop("OnePropGetFunction", "one_hamil not implemented", -1)
    else if (allocated(one_prop%multipole)) then
      call error_stop("OnePropGetFunction", "multipole not implemented", -1)
    else if (allocated(one_prop%gaussian_pot)) then
      call error_stop("OnePropGetFunction", "gaussian_pot not implemented", -1)
    else if (allocated(one_prop%pso)) then
      call error_stop("OnePropGetFunction", "pso not implemented", -1)
    else
      call error_stop("OnePropGetFunction", "unknown one-electron property integrand", -1)
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OnePropGetFunction", STDOUT)
#endif
  end subroutine OnePropGetFunction

  !> \brief frees space taken by the information of one-electron property integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param one_prop contains the information of one-electron property integrals
  subroutine OnePropDestroy(one_prop)
    type(one_prop_t), intent(inout) :: one_prop
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (allocated(one_prop%ang_mom)) then
      call AngMomDestroy(one_prop%ang_mom(1))
      deallocate(one_prop%ang_mom)
    else if (allocated(one_prop%overlap)) then
      call OverlapDestroy(one_prop%overlap(1))
      deallocate(one_prop%overlap)
    else if (allocated(one_prop%kin_energy)) then
      call KinEnergyDestroy(one_prop%kin_energy(1))
      deallocate(one_prop%kin_energy)
    else if (allocated(one_prop%pot_energy)) then
      call PotEnergyDestroy(one_prop%pot_energy(1))
      deallocate(one_prop%pot_energy)
    else if (allocated(one_prop%one_hamil)) then
      call OneHamilDestroy(one_prop%one_hamil(1))
      deallocate(one_prop%one_hamil)
    else if (allocated(one_prop%multipole)) then
      call MultipoleDestroy(one_prop%multipole(1))
      deallocate(one_prop%multipole)
    else if (allocated(one_prop%gaussian_pot)) then
      call GaussianPotDestroy(one_prop%gaussian_pot(1))
      deallocate(one_prop%gaussian_pot)
    else if (allocated(one_prop%pso)) then
      call PSODestroy(one_prop%pso(1))
      deallocate(one_prop%pso)
    end if
    ! resets the information of LAO
    call LondonAOSet(one_prop%info_LAO)
    one_prop%prop_init = .false.
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OnePropDestroy", STDOUT)
#endif
  end subroutine OnePropDestroy

end module gen1int
