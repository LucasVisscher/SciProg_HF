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
!!  This file provides the Fortran 90 module of Cartesian multipole moments.
!!
!!  2012-03-20, Bin Gao:
!!  * first version

#include "stdout.h"
#include "err_info.h"
#include "kind_matrix.h"
#include "max_idx_non.h"

!> \brief Fortran 90 module of Cartesian multipole moments
!> \author Bin Gao
!> \date 2012-03-20
module gen1int_carmom

  use xkind
  use london_ao
  use gen1int_geom
  implicit none

  ! angular momentum integrals
  type, public :: ang_mom_t
    private
    ! number of property integral matrices
    integer :: num_prop = 3
    ! symmetry of property integral matrices
    integer :: prop_sym = ANTI_SYM_MATRIX
    ! atomic center of dipole origin
    integer :: idx_diporg = MAX_IDX_NON
    ! coordinates of dipole origin
    real(REALK) :: dipole_origin(3) = (/0.0_REALK,0.0_REALK,0.0_REALK/)
  end type ang_mom_t

  ! overlap integrals
  type, public :: overlap_t
    private
    ! number of property integral matrices
    integer :: num_prop = 1
    ! symmetry of property integral matrices
    integer :: prop_sym = SYMMETRIC_MATRIX
    ! coordinates of dipole origin
    real(REALK) :: dipole_origin(3) = (/0.0_REALK,0.0_REALK,0.0_REALK/)
  end type overlap_t

  ! kinetic energy integrals
  type, public :: kin_energy_t
    private
    ! number of property integral matrices
    integer :: num_prop = 1
    ! symmetry of property integral matrices
    integer :: prop_sym = SYMMETRIC_MATRIX
    ! coordinates of dipole origin
    real(REALK) :: dipole_origin(3) = (/0.0_REALK,0.0_REALK,0.0_REALK/)
  end type kin_energy_t

  ! multipole integrals
  type, public :: multipole_t
    private
    ! number of property integral matrices
    integer :: num_prop = 1
    ! symmetry of property integral matrices
    integer :: prop_sym = SYMMETRIC_MATRIX
    ! order of multipole integrals
    integer :: order_mom = 0
    ! if calculating spherical multipole integrals
    logical :: spher_mom = .false.
    ! atomic center of dipole origin
    integer :: idx_diporg = MAX_IDX_NON
    ! coordinates of dipole origin
    real(REALK) :: dipole_origin(3) = (/0.0_REALK,0.0_REALK,0.0_REALK/)
    ! order of electronic derivatives
    integer :: order_elec = 0
  end type multipole_t

  public :: AngMomCreate
  public :: AngMomView
  public :: AngMomGetNumProp
  public :: AngMomGetSymmetry
  public :: AngMomGetIntegral
  public :: AngMomDestroy

  public :: OverlapCreate
  public :: OverlapView
  public :: OverlapGetNumProp
  public :: OverlapGetSymmetry
  public :: OverlapGetIntegral
  public :: OverlapGetFunction
  public :: OverlapDestroy

  public :: KinEnergyCreate
  public :: KinEnergyView
  public :: KinEnergyGetNumProp
  public :: KinEnergyGetSymmetry
  public :: KinEnergyGetIntegral
  public :: KinEnergyDestroy

  public :: MultipoleCreate
  public :: MultipoleView
  public :: MultipoleGetNumProp
  public :: MultipoleGetSymmetry
  public :: MultipoleGetIntegral
  public :: MultipoleDestroy

  contains

  !> \brief initializes the information of angular momentum integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param idx_diporg is the atomic center of dipole origin (<1 for non-atomic center)
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \return ang_mom contains the information of angular momentum integrals
  !> \return info_prop (==ERR_INFO) indicates the angular momentum integrals are not successfully created
  subroutine AngMomCreate(ang_mom, info_prop, idx_diporg, dipole_origin)
    type(ang_mom_t), intent(inout) :: ang_mom
    integer, intent(out) :: info_prop
    integer, optional, intent(in) :: idx_diporg
    real(REALK), optional, intent(in) :: dipole_origin(3)
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (present(idx_diporg)) then
      ang_mom%idx_diporg = idx_diporg
      if (ang_mom%idx_diporg>=1) ang_mom%prop_sym = SQUARE_MATRIX
    end if
    if (present(dipole_origin)) ang_mom%dipole_origin = dipole_origin
    info_prop = 0
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "AngMomCreate", STDOUT)
#endif
  end subroutine AngMomCreate

  !> \brief visualizes the information of angular momentum integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param ang_mom contains the information of angular momentum integrals
  !> \param io_viewer is the logical unit number of the viewer
  subroutine AngMomView(ang_mom, io_viewer)
    type(ang_mom_t), intent(in) :: ang_mom
    integer, intent(in) :: io_viewer
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    write(io_viewer,100) "angular momentum integrals"
    write(io_viewer,100) "number of property integral matrices", ang_mom%num_prop
    write(io_viewer,100) "symmetry of property integral matrices", ang_mom%prop_sym
    write(io_viewer,100) "atomic center of dipole origin", ang_mom%idx_diporg
    write(io_viewer,110) "coordinates of dipole origin", ang_mom%dipole_origin
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "AngMomView", STDOUT)
#endif
100 format("AngMomView>> ",A,I6)
110 format("AngMomView>> ",A,3F12.6)
  end subroutine AngMomView

  !> \brief returns the number of integral matrices for given angular momentum integrals
  !> \param ang_mom contains the information of angular momentum integrals
  !> \return num_prop is the number of property integral matrices
  subroutine AngMomGetNumProp(ang_mom, num_prop)
    type(ang_mom_t), intent(in) :: ang_mom
    integer, intent(out) :: num_prop
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    num_prop = ang_mom%num_prop
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "AngMomGetNumProp", STDOUT)
#endif
  end subroutine AngMomGetNumProp

  !> \brief returns the symmetry of integral matrices for given angular momentum integrals
  !> \param ang_mom contains the information of angular momentum integrals
  !> \return prop_sym indicates the symmetry of property integral matrices
  subroutine AngMomGetSymmetry(ang_mom, prop_sym)
    type(ang_mom_t), intent(in) :: ang_mom
    integer, intent(out) :: prop_sym
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    prop_sym = ang_mom%prop_sym
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "AngMomGetSymmetry", STDOUT)
#endif
  end subroutine AngMomGetSymmetry

  !> \brief evaluates the angular momentum integrals
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
  !> \param ang_mom contains the information of angular momentum integrals
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
  subroutine AngMomGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                               exponent_bra, num_contr_bra, contr_coef_bra,   &
                               idx_ket, coord_ket, angular_ket, num_prim_ket, &
                               exponent_ket, num_contr_ket, contr_coef_ket,   &
                               spher_gto, info_LAO, ang_mom,                  &
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
    type(ang_mom_t), intent(in) :: ang_mom
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
    integer iopt, jopt                                 !incremental recorder over operators
    integer ierr                                       !error information
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (present(nary_tree_total)) then
      call NaryTreePathGetNumCenters(nary_tree_total, num_diff_cent)
      if (num_diff_cent>2) then
        contr_ints = 0.0_REALK
        return
      end if
    end if
    allocate(tmp_ints(num_gto_bra,num_contr_bra,num_gto_ket,num_contr_ket, &
                      9,num_opt/3), stat=ierr)
    if (ierr/=0)                                                          &
      call error_stop("AngMomGetIntegral", "failed to allocate tmp_ints", &
                      num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket*6*num_opt)
    call IntGetCARMOM(idx_bra, coord_bra, angular_bra, num_prim_bra,    &
                      exponent_bra, num_contr_bra, contr_coef_bra,      &
                      idx_ket, coord_ket, angular_ket, num_prim_ket,    &
                      exponent_ket, num_contr_ket, contr_coef_ket,      &
                      spher_gto, info_LAO, 1, ang_mom%idx_diporg,       &
                      ang_mom%dipole_origin, 1.0_REALK, 1,              &
                      order_mag_bra, order_mag_ket, order_mag_total,    &
                      order_ram_bra, order_ram_ket, order_ram_total,    &
                      order_geo_bra, order_geo_ket, 0, nary_tree_total, &
                      num_gto_bra, num_gto_ket, 3*num_opt, tmp_ints,    &
                      mag_num_bra, mag_num_ket, powers_bra, powers_ket)
    ! we get: x*d/dx, x*d/dy, x*d/dz, y*d/dx, y*d/dy, y*d/dz, z*d/dx, z*d/dy, z*d/dz;
    ! angular momentum operators are: -(y*d/dz-z*d/dy), -(z*d/dx-x*d/dz), -(x*d/dy-y*d/dx)
    jopt = 0
    do iopt = 1, num_opt/3
      jopt = jopt+1
      contr_ints(:,:,:,:,jopt) = tmp_ints(:,:,:,:,8,iopt)-tmp_ints(:,:,:,:,6,iopt)
      jopt = jopt+1
      contr_ints(:,:,:,:,jopt) = tmp_ints(:,:,:,:,3,iopt)-tmp_ints(:,:,:,:,7,iopt)
      jopt = jopt+1
      contr_ints(:,:,:,:,jopt) = tmp_ints(:,:,:,:,4,iopt)-tmp_ints(:,:,:,:,2,iopt)
    end do
    deallocate(tmp_ints)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "AngMomGetIntegral", STDOUT)
#endif
  end subroutine AngMomGetIntegral

  !> \brief frees space taken by the angular momentum integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param ang_mom contains the information of angular momentum integrals
  subroutine AngMomDestroy(ang_mom)
    type(ang_mom_t), intent(inout) :: ang_mom
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ang_mom%prop_sym = ANTI_SYM_MATRIX
    ang_mom%idx_diporg = MAX_IDX_NON
    ang_mom%dipole_origin = 0.0_REALK
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "AngMomDestroy", STDOUT)
#endif
  end subroutine AngMomDestroy

  !> \brief initializes the information of overlap integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \return overlap contains the information of overlap integrals
  !> \return info_prop (==ERR_INFO) indicates the overlap integrals are not successfully created
  subroutine OverlapCreate(overlap, info_prop)
    type(overlap_t), intent(inout) :: overlap
    integer, intent(out) :: info_prop
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    info_prop = 0
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OverlapCreate", STDOUT)
#endif
  end subroutine OverlapCreate

  !> \brief visualizes the information of overlap integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param overlap contains the information of overlap integrals
  !> \param io_viewer is the logical unit number of the viewer
  subroutine OverlapView(overlap, io_viewer)
    type(overlap_t), intent(in) :: overlap
    integer, intent(in) :: io_viewer
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    write(io_viewer,100) "overlap integrals"
    write(io_viewer,100) "number of property integral matrices", overlap%num_prop
    write(io_viewer,100) "symmetry of property integral matrices", overlap%prop_sym
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OverlapView", STDOUT)
#endif
100 format("OverlapView>> ",A,I6)
  end subroutine OverlapView

  !> \brief returns the number of integral matrices for given overlap integrals
  !> \param overlap contains the information of overlap integrals
  !> \return num_prop is the number of property integral matrices
  subroutine OverlapGetNumProp(overlap, num_prop)
    type(overlap_t), intent(in) :: overlap
    integer, intent(out) :: num_prop
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    num_prop = overlap%num_prop
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OverlapGetNumProp", STDOUT)
#endif
  end subroutine OverlapGetNumProp

  !> \brief returns the symmetry of integral matrices for given overlap integrals
  !> \param overlap contains the information of overlap integrals
  !> \return prop_sym indicates the symmetry of property integral matrices
  subroutine OverlapGetSymmetry(overlap, prop_sym)
    type(overlap_t), intent(in) :: overlap
    integer, intent(out) :: prop_sym
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    prop_sym = overlap%prop_sym
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OverlapGetSymmetry", STDOUT)
#endif
  end subroutine OverlapGetSymmetry

  !> \brief evaluates the overlap integrals
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
  !> \param overlap contains the information of overlap integrals
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
  subroutine OverlapGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                exponent_bra, num_contr_bra, contr_coef_bra,   &
                                idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                exponent_ket, num_contr_ket, contr_coef_ket,   &
                                spher_gto, info_LAO, overlap,                  &
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
    type(overlap_t), intent(in) :: overlap
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
    integer num_diff_cent  !number of differentiated centers of current path
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (present(nary_tree_total)) then
      call NaryTreePathGetNumCenters(nary_tree_total, num_diff_cent)
      if (num_diff_cent>2) then
        contr_ints = 0.0_REALK
        return
      end if
    end if
    call IntGetCARMOM(idx_bra, coord_bra, angular_bra, num_prim_bra,    &
                      exponent_bra, num_contr_bra, contr_coef_bra,      &
                      idx_ket, coord_ket, angular_ket, num_prim_ket,    &
                      exponent_ket, num_contr_ket, contr_coef_ket,      &
                      spher_gto, info_LAO, 0, MAX_IDX_NON,              &
                      overlap%dipole_origin, 1.0_REALK, 0,              &
                      order_mag_bra, order_mag_ket, order_mag_total,    &
                      order_ram_bra, order_ram_ket, order_ram_total,    &
                      order_geo_bra, order_geo_ket, 0, nary_tree_total, &
                      num_gto_bra, num_gto_ket, num_opt, contr_ints,    &
                      mag_num_bra, mag_num_ket, powers_bra, powers_ket)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OverlapGetIntegral", STDOUT)
#endif
  end subroutine OverlapGetIntegral

  !> \brief evaluates the overlap distribution
  !> \author Bin Gao
  !> \date 2012-02-10
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
  !> \param overlap contains the information of overlap integrals
  !> \param order_mag_bra is the order of magnetic derivatives on bra center
  !> \param order_mag_ket is the order of magnetic derivatives on ket center
  !> \param order_mag_total is the order of total magnetic derivatives
  !> \param order_ram_bra is the order of derivatives w.r.t. total rotational angular momentum on bra center
  !> \param order_ram_ket is the order of derivatives w.r.t. total rotational angular momentum on ket center
  !> \param order_ram_total is the order of total derivatives w.r.t. total rotational angular momentum
  !> \param order_geo_bra is the order of geometric derivatives with respect to bra center
  !> \param order_geo_ket is the order of geometric derivatives with respect to ket center
  !> \param nary_tree_total contains the information of N-ary tree for total geometric derivatives
  !> \param num_points is the number of grid points
  !> \param grid_points contains the coordinates of grid points
  !> \param num_gto_bra is the number of spherical/Cartesian GTOs on bra center
  !> \param num_gto_ket is the number of spherical/Cartesian GTOs on ket center
  !> \param num_opt is the number of operators including derivatives
  !> \param mag_num_bra contains the magnetic numbers of spherical GTOs on bra center
  !> \param mag_num_ket contains the magnetic numbers of spherical GTOs on ket center
  !> \param powers_bra contains the Cartesian powers of Cartesian GTOs on bra center
  !> \param powers_ket contains the Cartesian powers of Cartesian GTOs on ket center
  !> \return contr_val contains the calculated contracted overlap distribution
  subroutine OverlapGetFunction(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                exponent_bra, num_contr_bra, contr_coef_bra,   &
                                idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                exponent_ket, num_contr_ket, contr_coef_ket,   &
                                spher_gto, info_LAO, overlap,                  &
                                order_mag_bra, order_mag_ket, order_mag_total, &
                                order_ram_bra, order_ram_ket, order_ram_total, &
                                order_geo_bra, order_geo_ket, nary_tree_total, &
                                num_points, grid_points,                       &
                                num_gto_bra, num_gto_ket, num_opt, contr_val,  &
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
    type(overlap_t), intent(in) :: overlap
    integer, intent(in) :: order_mag_bra
    integer, intent(in) :: order_mag_ket
    integer, intent(in) :: order_mag_total
    integer, intent(in) :: order_ram_bra
    integer, intent(in) :: order_ram_ket
    integer, intent(in) :: order_ram_total
    integer, optional, intent(in) :: order_geo_bra
    integer, optional, intent(in) :: order_geo_ket
    type(nary_tree_t), optional, intent(in) :: nary_tree_total
    integer, intent(in) :: num_points
    real(REALK), intent(in) :: grid_points(3,num_points)
    integer, intent(in) :: num_gto_bra
    integer, intent(in) :: num_gto_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(out) :: contr_val(num_gto_bra,num_contr_bra, &
                                          num_gto_ket,num_contr_ket, &
                                          num_points,num_opt)
    integer, optional, intent(in) :: mag_num_bra(num_gto_bra)
    integer, optional, intent(in) :: mag_num_ket(num_gto_ket)
    integer, optional, intent(in) :: powers_bra(3,num_gto_bra)
    integer, optional, intent(in) :: powers_ket(3,num_gto_ket)
    integer num_diff_cent  !number of differentiated centers of current path
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (present(nary_tree_total)) then
      call NaryTreePathGetNumCenters(nary_tree_total, num_diff_cent)
      if (num_diff_cent>2) then
        contr_val = 0.0_REALK
        return
      end if
    end if
    call IntGetODST(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                    exponent_bra, num_contr_bra, contr_coef_bra,   &
                    idx_ket, coord_ket, angular_ket, num_prim_ket, &
                    exponent_ket, num_contr_ket, contr_coef_ket,   &
                    spher_gto, info_LAO,                           &
                    order_mag_bra, order_mag_ket, order_mag_total, &
                    order_ram_bra, order_ram_ket, order_ram_total, &
                    order_geo_bra, order_geo_ket, nary_tree_total, &
                    num_points, grid_points,                       &
                    num_gto_bra, num_gto_ket, num_opt, contr_val,  &
                    mag_num_bra, mag_num_ket, powers_bra, powers_ket)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OverlapGetFunction", STDOUT)
#endif
  end subroutine OverlapGetFunction

  !> \brief frees space taken by the overlap integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param overlap contains the information of overlap integrals
  subroutine OverlapDestroy(overlap)
    type(overlap_t), intent(inout) :: overlap
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    overlap%dipole_origin = 0.0_REALK
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "OverlapDestroy", STDOUT)
#endif
  end subroutine OverlapDestroy

  !> \brief initializes the information of kinetic energy integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \return kin_energy contains the information of kinetic energy integrals
  !> \return info_prop (==ERR_INFO) indicates the kinetic energy integrals are not successfully created
  subroutine KinEnergyCreate(kin_energy, info_prop)
    type(kin_energy_t), intent(inout) :: kin_energy
    integer, intent(out) :: info_prop
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    info_prop = 0
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "KinEnergyCreate", STDOUT)
#endif
  end subroutine KinEnergyCreate

  !> \brief visualizes the information of kinetic energy integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param kin_energy contains the information of kinetic energy integrals
  !> \param io_viewer is the logical unit number of the viewer
  subroutine KinEnergyView(kin_energy, io_viewer)
    type(kin_energy_t), intent(in) :: kin_energy
    integer, intent(in) :: io_viewer
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    write(io_viewer,100) "kinetic energy integrals"
    write(io_viewer,100) "number of property integral matrices", kin_energy%num_prop
    write(io_viewer,100) "symmetry of property integral matrices", kin_energy%prop_sym
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "KinEnergyView", STDOUT)
#endif
100 format("KinEnergyView>> ",A,I6)
  end subroutine KinEnergyView

  !> \brief returns the number of integral matrices for given kinetic energy integrals
  !> \param kin_energy contains the information of kinetic energy integrals
  !> \return num_prop is the number of property integral matrices
  subroutine KinEnergyGetNumProp(kin_energy, num_prop)
    type(kin_energy_t), intent(in) :: kin_energy
    integer, intent(out) :: num_prop
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    num_prop = kin_energy%num_prop
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "KinEnergyGetNumProp", STDOUT)
#endif
  end subroutine KinEnergyGetNumProp

  !> \brief returns the symmetry of integral matrices for given kinetic energy integrals
  !> \param kin_energy contains the information of kinetic energy integrals
  !> \return prop_sym indicates the symmetry of property integral matrices
  subroutine KinEnergyGetSymmetry(kin_energy, prop_sym)
    type(kin_energy_t), intent(in) :: kin_energy
    integer, intent(out) :: prop_sym
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    prop_sym = kin_energy%prop_sym
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "KinEnergyGetSymmetry", STDOUT)
#endif
  end subroutine KinEnergyGetSymmetry

  !> \brief evaluates the kinetic energy integrals
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
  !> \param kin_energy contains the information of kinetic energy integrals
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
  subroutine KinEnergyGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                  exponent_bra, num_contr_bra, contr_coef_bra,   &
                                  idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                  exponent_ket, num_contr_ket, contr_coef_ket,   &
                                  spher_gto, info_LAO, kin_energy,               &
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
    type(kin_energy_t), intent(in) :: kin_energy
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
    type(ang_mom_t) ang_mom                            !angular momentum integrals
    real(REALK) gauge_origin(3)                        !gauge origin of the magnetic vector potential
    real(REALK) origin_London_PF(3)                    !origin of the London phase factor
    integer p_order_geo_bra
    integer p_order_geo_ket
    integer p_order_geo_total
    integer p_idx_cent_total(1)
    integer iopt                                       !incremental recorder over operators
    integer num_diff_cent                              !number of differentiated centers of current path
    real(REALK), allocatable :: tmp_ints(:,:,:,:,:,:)  !temporary integrals
    integer ierr                                       !error information
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (present(nary_tree_total)) then
      call NaryTreePathGetNumCenters(nary_tree_total, num_diff_cent)
      if (num_diff_cent>2) then
        contr_ints = 0.0_REALK
        return
      end if
    end if
    ! external magnetic field
    if (order_mag_total>0 .or. order_mag_bra>0 .or. order_mag_ket>0 .or. &
        order_ram_total>0 .or. order_ram_bra>0 .or. order_ram_ket>0) then
      ! only first order total magnetic derivatives implemented
      if (order_mag_total>1 .or. order_mag_bra>0 .or. order_mag_ket>0 .or. &
          order_ram_total>0 .or. order_ram_bra>0 .or. order_ram_ket>0) then
        call error_stop("KinEnergyGetIntegral", &
                        "only first order total magnetic derivatives implemented", 1)
      end if
      ! LAO used
      if (LondonAOUsed(info_LAO)) then
        if (present(order_geo_bra)) then
          p_order_geo_bra = order_geo_bra
        else
          p_order_geo_bra = 0
        end if
        if (present(order_geo_ket)) then
          p_order_geo_ket = order_geo_ket
        else
          p_order_geo_ket = 0
        end if
        if (present(nary_tree_total)) then
          call NaryTreeGetOrder(nary_tree=nary_tree_total, &
                                order_geo=p_order_geo_total)
        else
          p_order_geo_total = 0
        end if
        if (p_order_geo_bra>0 .or. p_order_geo_ket>0 .or. p_order_geo_total>1) then
          call error_stop("KinEnergyGetIntegral", &
                          "no mixed partial geometric and magnetic derivatives", 1)
        ! magnetic derivatives (and mixed total geometric derivatives)
        else
          ! L_{N} = -i*r_{N}x\nabla
          call AngMomCreate(ang_mom=ang_mom,    &
                            info_prop=ierr,     &
                            idx_diporg=idx_ket, &
                            dipole_origin=coord_ket)
          if (ierr/=0) then
            call error_stop("KinEnergyGetIntegral", "failed to call AngMomCreate", ierr)
          end if
          call AngMomGetIntegral(idx_bra=idx_bra, &
                                 coord_bra=coord_bra, &
                                 angular_bra=angular_bra, &
                                 num_prim_bra=num_prim_bra, &
                                 exponent_bra=exponent_bra, &
                                 num_contr_bra=num_contr_bra, &
                                 contr_coef_bra=contr_coef_bra, &
                                 idx_ket=idx_ket, &
                                 coord_ket=coord_ket, &
                                 angular_ket=angular_ket, &
                                 num_prim_ket=num_prim_ket, &
                                 exponent_ket=exponent_ket, &
                                 num_contr_ket=num_contr_ket, &
                                 contr_coef_ket=contr_coef_ket, &
                                 spher_gto=spher_gto, &
                                 info_LAO=info_LAO, &
                                 ang_mom=ang_mom, &
                                 order_mag_bra=0, &
                                 order_mag_ket=0, &
                                 order_mag_total=0, &
                                 order_ram_bra=0, &
                                 order_ram_ket=0, &
                                 order_ram_total=0, &
                                 order_geo_bra=order_geo_bra, &
                                 order_geo_ket=order_geo_ket, &
                                 nary_tree_total=nary_tree_total, &
                                 num_gto_bra=num_gto_bra, &
                                 num_gto_ket=num_gto_ket, &
                                 num_opt=num_opt, &
                                 contr_ints=contr_ints, &
                                 mag_num_bra=mag_num_bra, &
                                 mag_num_ket=mag_num_ket, &
                                 powers_bra=powers_bra, &
                                 powers_ket=powers_ket)
          call AngMomDestroy(ang_mom=ang_mom)
          if (idx_bra/=idx_ket) then
            ! -1/2*Q_{MN}*r_{P}*\nabla^{2} with total geometric derivatives
            allocate(tmp_ints(num_gto_bra,num_contr_bra,num_gto_ket,num_contr_ket, &
                              6,num_opt), stat=ierr)
            if (ierr/=0)                                                             &
              call error_stop("KinEnergyGetIntegral", "failed to allocate tmp_ints", &
                              num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket*6*num_opt)
            call LondonAOGetLPFOrigin(info_LAO=info_LAO, &
                                      origin_London_PF=origin_London_PF)
            call IntGetCARMOM(idx_bra=idx_bra, &
                              coord_bra=coord_bra, &
                              angular_bra=angular_bra, &
                              num_prim_bra=num_prim_bra, &
                              exponent_bra=exponent_bra, &
                              num_contr_bra=num_contr_bra, &
                              contr_coef_bra=contr_coef_bra, &
                              idx_ket=idx_ket, &
                              coord_ket=coord_ket, &
                              angular_ket=angular_ket, &
                              num_prim_ket=num_prim_ket, &
                              exponent_ket=exponent_ket, &
                              num_contr_ket=num_contr_ket, &
                              contr_coef_ket=contr_coef_ket, &
                              spher_gto=spher_gto, &
                              info_LAO=info_LAO, &
                              order_elec=2, &
                              idx_diporg=MAX_IDX_NON, &
                              dipole_origin=origin_London_PF, &
                              scal_const=-0.5_REALK, &
                              order_mom=1, &
                              order_mag_bra=0, &
                              order_mag_ket=0, &
                              order_mag_total=0, &
                              order_ram_bra=0, &
                              order_ram_ket=0, &
                              order_ram_total=0, &
                              order_geo_bra=order_geo_bra, &
                              order_geo_ket=order_geo_ket, &
                              order_geo_mom=0, &
                              nary_tree_total=nary_tree_total, &
                              num_gto_bra=num_gto_bra, &
                              num_gto_ket=num_gto_ket, &
                              num_opt=6*num_opt, &
                              contr_ints=tmp_ints, &
                              mag_num_bra=mag_num_bra, &
                              mag_num_ket=mag_num_ket, &
                              powers_bra=powers_bra, &
                              powers_ket=powers_ket)
            ! sums xx, yy and zz components
            tmp_ints(:,:,:,:,2,:) = tmp_ints(:,:,:,:,1,:)+tmp_ints(:,:,:,:,3,:)+tmp_ints(:,:,:,:,6,:)
            ! sums L_{N} and -1/2*Q_{MN}*r_{P}*\nabla^{2}
            do iopt = 1, num_opt-2, 3
              contr_ints(:,:,:,:,iopt) = contr_ints(:,:,:,:,iopt)                               &
                                       + (coord_bra(2)-coord_ket(2))*tmp_ints(:,:,:,:,2,iopt+2) &
                                       - (coord_bra(3)-coord_ket(3))*tmp_ints(:,:,:,:,2,iopt+1)
              contr_ints(:,:,:,:,iopt+1) = contr_ints(:,:,:,:,iopt+1)                           &
                                         + (coord_bra(3)-coord_ket(3))*tmp_ints(:,:,:,:,2,iopt) &
                                         - (coord_bra(1)-coord_ket(1))*tmp_ints(:,:,:,:,2,iopt+2)
              contr_ints(:,:,:,:,iopt+2) = contr_ints(:,:,:,:,iopt+2)                             &
                                         + (coord_bra(1)-coord_ket(1))*tmp_ints(:,:,:,:,2,iopt+1) &
                                         - (coord_bra(2)-coord_ket(2))*tmp_ints(:,:,:,:,2,iopt)
            end do

            ! mixed total geometric and magnetic derivatives
            if (p_order_geo_total==1) then
              call NaryTreePathGetIdxCent(nary_tree=nary_tree_total, &
                                          idx_cent=p_idx_cent_total)
              ! geometric derivatives on the bra or ket center
              if (p_idx_cent_total(1)==idx_bra .or. p_idx_cent_total(1)==idx_ket) then
                tmp_ints = 0.0_REALK
                ! -1/2*(Q_{MN})^{GEO}*r_{P}*\nabla^{2}
                call IntGetCARMOM(idx_bra=idx_bra, &
                                  coord_bra=coord_bra, &
                                  angular_bra=angular_bra, &
                                  num_prim_bra=num_prim_bra, &
                                  exponent_bra=exponent_bra, &
                                  num_contr_bra=num_contr_bra, &
                                  contr_coef_bra=contr_coef_bra, &
                                  idx_ket=idx_ket, &
                                  coord_ket=coord_ket, &
                                  angular_ket=angular_ket, &
                                  num_prim_ket=num_prim_ket, &
                                  exponent_ket=exponent_ket, &
                                  num_contr_ket=num_contr_ket, &
                                  contr_coef_ket=contr_coef_ket, &
                                  spher_gto=spher_gto, &
                                  info_LAO=info_LAO, &
                                  order_elec=2, &
                                  idx_diporg=MAX_IDX_NON, &
                                  dipole_origin=origin_London_PF, &
                                  scal_const=-0.5_REALK, &
                                  order_mom=1, &
                                  order_mag_bra=0, &
                                  order_mag_ket=0, &
                                  order_mag_total=0, &
                                  order_ram_bra=0, &
                                  order_ram_ket=0, &
                                  order_ram_total=0, &
                                  order_geo_bra=order_geo_bra, &
                                  order_geo_ket=order_geo_ket, &
                                  order_geo_mom=0, &
                                  num_gto_bra=num_gto_bra, &
                                  num_gto_ket=num_gto_ket, &
                                  num_opt=18, &
                                  contr_ints=tmp_ints(:,:,:,:,:,1:3), &
                                  mag_num_bra=mag_num_bra, &
                                  mag_num_ket=mag_num_ket, &
                                  powers_bra=powers_bra, &
                                  powers_ket=powers_ket)
                ! sums xx, yy and zz components
                tmp_ints(:,:,:,:,2,1:3) = tmp_ints(:,:,:,:,1,1:3) &
                                        + tmp_ints(:,:,:,:,3,1:3) &
                                        + tmp_ints(:,:,:,:,6,1:3)
                ! adds -1/2*Q_{MN}^{GEO}*r_{P}*\nabla^{2}
                if (p_idx_cent_total(1)==idx_bra) then
                    ! Bx, Gx
                    ! By, Gx
                    contr_ints(:,:,:,:,2) = contr_ints(:,:,:,:,2)-tmp_ints(:,:,:,:,2,3)
                    ! Bz, Gx
                    contr_ints(:,:,:,:,3) = contr_ints(:,:,:,:,3)+tmp_ints(:,:,:,:,2,2)
                    ! Bx, Gy
                    contr_ints(:,:,:,:,4) = contr_ints(:,:,:,:,4)+tmp_ints(:,:,:,:,2,3)
                    ! By, Gy
                    ! Bz, Gy
                    contr_ints(:,:,:,:,6) = contr_ints(:,:,:,:,6)-tmp_ints(:,:,:,:,2,1)
                    ! Bx, Gz
                    contr_ints(:,:,:,:,7) = contr_ints(:,:,:,:,7)-tmp_ints(:,:,:,:,2,2)
                    ! By, Gz
                    contr_ints(:,:,:,:,8) = contr_ints(:,:,:,:,8)+tmp_ints(:,:,:,:,2,1)
                    ! Bz, Gz
                else
                    ! Bx, Gx
                    ! By, Gx
                    contr_ints(:,:,:,:,2) = contr_ints(:,:,:,:,2)+tmp_ints(:,:,:,:,2,3)
                    ! Bz, Gx
                    contr_ints(:,:,:,:,3) = contr_ints(:,:,:,:,3)-tmp_ints(:,:,:,:,2,2)
                    ! Bx, Gy
                    contr_ints(:,:,:,:,4) = contr_ints(:,:,:,:,4)-tmp_ints(:,:,:,:,2,3)
                    ! By, Gy
                    ! Bz, Gy
                    contr_ints(:,:,:,:,6) = contr_ints(:,:,:,:,6)+tmp_ints(:,:,:,:,2,1)
                    ! Bx, Gz
                    contr_ints(:,:,:,:,7) = contr_ints(:,:,:,:,7)+tmp_ints(:,:,:,:,2,2)
                    ! By, Gz
                    contr_ints(:,:,:,:,8) = contr_ints(:,:,:,:,8)-tmp_ints(:,:,:,:,2,1)
                    ! Bz, Gz
                end if
              end if
            end if

            deallocate(tmp_ints)
          end if
          contr_ints = 0.5_REALK*contr_ints
        end if
      ! non-LAO
      else
        call LondonAOGetGaugeOrigin(info_LAO, gauge_origin)
        call AngMomCreate(ang_mom, ierr, MAX_IDX_NON, gauge_origin)
        if (ierr/=0) then
          call error_stop("KinEnergyGetIntegral", "failed to call AngMomCreate", ierr)
        end if
        call AngMomGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                               exponent_bra, num_contr_bra, contr_coef_bra,   &
                               idx_ket, coord_ket, angular_ket, num_prim_ket, &
                               exponent_ket, num_contr_ket, contr_coef_ket,   &
                               spher_gto, info_LAO, ang_mom,                  &
                               0, 0, 0,                                       &
                               0, 0, 0,                                       &
                               order_geo_bra, order_geo_ket, nary_tree_total, &
                               num_gto_bra, num_gto_ket, num_opt, contr_ints, &
                               mag_num_bra, mag_num_ket, powers_bra, powers_ket)
        contr_ints = 0.5_REALK*contr_ints
        call AngMomDestroy(ang_mom)
      end if
    ! no magnetic derivatives
    else
      allocate(tmp_ints(num_gto_bra,num_contr_bra,num_gto_ket,num_contr_ket, &
                        6,num_opt), stat=ierr)
      if (ierr/=0)                                                             &
        call error_stop("KinEnergyGetIntegral", "failed to allocate tmp_ints", &
                        num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket*6*num_opt)
      call IntGetCARMOM(idx_bra, coord_bra, angular_bra, num_prim_bra,    &
                        exponent_bra, num_contr_bra, contr_coef_bra,      &
                        idx_ket, coord_ket, angular_ket, num_prim_ket,    &
                        exponent_ket, num_contr_ket, contr_coef_ket,      &
                        spher_gto, info_LAO, 2, MAX_IDX_NON,              &
                        kin_energy%dipole_origin, -0.5_REALK, 0,          &
                        order_mag_bra, order_mag_ket, order_mag_total,    &
                        order_ram_bra, order_ram_ket, order_ram_total,    &
                        order_geo_bra, order_geo_ket, 0, nary_tree_total, &
                        num_gto_bra, num_gto_ket, 6*num_opt, tmp_ints,    &
                        mag_num_bra, mag_num_ket, powers_bra, powers_ket)
      ! sums xx, yy and zz components
      contr_ints = tmp_ints(:,:,:,:,1,:)+tmp_ints(:,:,:,:,3,:)+tmp_ints(:,:,:,:,6,:)
      deallocate(tmp_ints)
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "KinEnergyGetIntegral", STDOUT)
#endif
  end subroutine KinEnergyGetIntegral

  !> \brief frees space taken by the kinetic energy integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param kin_energy contains the information of kinetic energy integrals
  subroutine KinEnergyDestroy(kin_energy)
    type(kin_energy_t), intent(inout) :: kin_energy
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    kin_energy%dipole_origin = 0.0_REALK
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "KinEnergyDestroy", STDOUT)
#endif
  end subroutine KinEnergyDestroy

  !> \brief initializes the information of multipole integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param idx_diporg is the atomic center of dipole origin (<1 for non-atomic center)
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param order_mom is the order of multipole integrals
  !> \param spher_mom indicates if calculating spherical multipole integrals
  !> \param order_elec is the order of electronic derivatives
  !> \return multipole contains the information of multipole integrals
  !> \return info_prop (==ERR_INFO) indicates the Multipole integrals are not successfully created
  subroutine MultipoleCreate(multipole, info_prop, idx_diporg, &
                             dipole_origin, order_mom, spher_mom, order_elec)
    type(multipole_t), intent(inout) :: multipole
    integer, intent(out) :: info_prop
    integer, optional, intent(in) :: idx_diporg
    real(REALK), optional, intent(in) :: dipole_origin(3)
    integer, optional, intent(in) :: order_mom
    logical, optional, intent(in) :: spher_mom
    integer, optional, intent(in) :: order_elec
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! checks the validity of input order
    if (present(order_mom)) then
      if (order_mom>=0) then
        multipole%order_mom = order_mom
        info_prop = 0
      ! negative order, returns with error code
      else
        info_prop = order_mom
        return
      end if
    end if
    if (present(idx_diporg)) multipole%idx_diporg = idx_diporg
    if (present(dipole_origin)) multipole%dipole_origin = dipole_origin
    if (present(spher_mom)) multipole%spher_mom = spher_mom
    if (present(order_elec)) then
      multipole%order_elec = max(0,order_elec)
    end if
    if (multipole%spher_mom) then
      multipole%num_prop = (2*multipole%order_mom+1)*(multipole%order_elec+1) &
                         * (multipole%order_elec+2)/2
    else
      multipole%num_prop = (multipole%order_mom+1)*(multipole%order_mom+2) &
                         * (multipole%order_elec+1)*(multipole%order_elec+2)/4
    end if
    if (mod(multipole%order_elec,2)==1) then
      multipole%prop_sym = ANTI_SYM_MATRIX
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "MultipoleCreate", STDOUT)
#endif
  end subroutine MultipoleCreate

  !> \brief visualizes the information of multipole integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param multipole contains the information of multipole integrals
  !> \param io_viewer is the logical unit number of the viewer
  subroutine MultipoleView(multipole, io_viewer)
    type(multipole_t), intent(in) :: multipole
    integer, intent(in) :: io_viewer
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (multipole%spher_mom) then
      write(io_viewer,100) "spherical multipole integrals"
    else
      write(io_viewer,100) "Cartesian multipole integrals"
    end if
    write(io_viewer,100) "number of property integral matrices", multipole%num_prop
    write(io_viewer,100) "symmetry of property integral matrices", multipole%prop_sym
    write(io_viewer,100) "order of multipole integrals", multipole%order_mom
    write(io_viewer,100) "atomic center of dipole origin", multipole%idx_diporg
    write(io_viewer,110) "coordinates of dipole origin", multipole%dipole_origin
    write(io_viewer,100) "order of electronic derivatives", multipole%order_elec
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "MultipoleView", STDOUT)
#endif
100 format("MultipoleView>> ",A,I6)
110 format("MultipoleView>> ",A,3F12.6)
  end subroutine MultipoleView

  !> \brief returns the number of integral matrices for given multipole integrals
  !> \param multipole contains the information of multipole integrals
  !> \return num_prop is the number of property integral matrices
  subroutine MultipoleGetNumProp(multipole, num_prop)
    type(multipole_t), intent(in) :: multipole
    integer, intent(out) :: num_prop
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    num_prop = multipole%num_prop
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "MultipoleGetNumProp", STDOUT)
#endif
  end subroutine MultipoleGetNumProp

  !> \brief returns the symmetry of integral matrices for given multipole integrals
  !> \param multipole contains the information of multipole integrals
  !> \return prop_sym indicates the symmetry of property integral matrices
  subroutine MultipoleGetSymmetry(multipole, prop_sym)
    type(multipole_t), intent(in) :: multipole
    integer, intent(out) :: prop_sym
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    prop_sym = multipole%prop_sym
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "MultipoleGetSymmetry", STDOUT)
#endif
  end subroutine MultipoleGetSymmetry

  !> \brief evaluates the multipole integrals
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
  !> \param multipole contains the information of multipole integrals
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
  subroutine MultipoleGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                  exponent_bra, num_contr_bra, contr_coef_bra,   &
                                  idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                  exponent_ket, num_contr_ket, contr_coef_ket,   &
                                  spher_gto, info_LAO, multipole,                &
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
    type(multipole_t), intent(in) :: multipole
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
    integer num_xyz_mom                              !number of xyz components of Cartesian multipole moment
    integer num_derv                                 !number of derivatives
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
      end if
    end if
    ! spherical multipole integrals
    if (multipole%spher_mom) then
      ! calculates the Cartesian multipole integrals
      num_xyz_mom = (multipole%order_mom+1)*(multipole%order_mom+2)/2
      num_derv = num_opt/multipole%num_prop
      allocate(tmp_ints(num_gto_bra,num_contr_bra,num_gto_ket,num_contr_ket, &
                        num_xyz_mom*num_derv), stat=ierr)
      if (ierr/=0)                                                             &
        call error_stop("MultipoleGetIntegral", "failed to allocate tmp_ints", &
                        num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket    &
                        *num_xyz_mom*num_derv)
      call IntGetCARMOM(idx_bra, coord_bra, angular_bra, num_prim_bra,    &
                        exponent_bra, num_contr_bra, contr_coef_bra,      &
                        idx_ket, coord_ket, angular_ket, num_prim_ket,    &
                        exponent_ket, num_contr_ket, contr_coef_ket,      &
                        spher_gto, info_LAO, multipole%order_elec,        &
                        multipole%idx_diporg, multipole%dipole_origin,    &
                        1.0_REALK, multipole%order_mom,                   &
                        order_mag_bra, order_mag_ket, order_mag_total,    &
                        order_ram_bra, order_ram_ket, order_ram_total,    &
                        order_geo_bra, order_geo_ket, 0, nary_tree_total, &
                        num_gto_bra, num_gto_ket, num_opt, tmp_ints,      &
                        mag_num_bra, mag_num_ket, powers_bra, powers_ket)
      ! transforms to spherical multipole integrals
      call hgto_to_sgto_oop(multipole%order_mom, size(contr_ints)/num_opt, num_xyz_mom, &
                            num_derv, tmp_ints, multipole%num_prop, contr_ints)
      deallocate(tmp_ints)
    ! Cartesian multipole integrals
    else
      call IntGetCARMOM(idx_bra, coord_bra, angular_bra, num_prim_bra,    &
                        exponent_bra, num_contr_bra, contr_coef_bra,      &
                        idx_ket, coord_ket, angular_ket, num_prim_ket,    &
                        exponent_ket, num_contr_ket, contr_coef_ket,      &
                        spher_gto, info_LAO, multipole%order_elec,        &
                        multipole%idx_diporg, multipole%dipole_origin,    &
                        1.0_REALK, multipole%order_mom,                   &
                        order_mag_bra, order_mag_ket, order_mag_total,    &
                        order_ram_bra, order_ram_ket, order_ram_total,    &
                        order_geo_bra, order_geo_ket, 0, nary_tree_total, &
                        num_gto_bra, num_gto_ket, num_opt, contr_ints,    &
                        mag_num_bra, mag_num_ket, powers_bra, powers_ket)
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "MultipoleGetIntegral", STDOUT)
#endif
  end subroutine MultipoleGetIntegral

  !> \brief frees space taken by the multipole integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param multipole contains the information of multipole integrals
  subroutine MultipoleDestroy(multipole)
    type(multipole_t), intent(inout) :: multipole
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    multipole%num_prop = 1
    multipole%order_mom = 0
    multipole%spher_mom = .false.
    multipole%idx_diporg = MAX_IDX_NON
    multipole%dipole_origin = 0.0_REALK
    multipole%order_elec = 0
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "MultipoleDestroy", STDOUT)
#endif
  end subroutine MultipoleDestroy

end module gen1int_carmom
