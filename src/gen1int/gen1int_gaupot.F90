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
!!  This file provides the Fortran 90 module of Gaussian charge potential integrals.
!!
!!  2012-03-20, Bin Gao:
!!  * first version

#include "stdout.h"
#include "err_info.h"
#include "kind_matrix.h"
#include "max_idx_non.h"

!> \brief Fortran 90 module of Gaussian charge potential integrals
!> \author Bin Gao
!> \date 2012-03-20
module gen1int_gaupot

  use xkind
  use london_ao
  use gen1int_geom
  implicit none

  ! Gaussian charge potential integrals
  type, public :: gaussian_pot_t
    private
    ! number of property integral matrices
    integer :: num_prop = 1
    ! symmetry of property integral matrices
    integer :: prop_sym = SYMMETRIC_MATRIX
    ! number of Gaussian charge potential origins
    integer :: num_origin = 0
    ! atomic centers of Gaussian charge potential origins (<1 for non-atomic center)
    integer, allocatable :: idx_gauorg(:)
    ! coordinates of Gaussian charge potential origins
    real(REALK), allocatable :: gaupot_origin(:,:)
    ! charges of Gaussian charge potential
    real(REALK), allocatable :: gaupot_charge(:)
    ! exponets of Gaussian charge potential
    real(REALK), allocatable :: gaupot_expt(:)
    ! order of geometric derivatives with respect to the Gaussian charge potential origins
    integer :: order_geo_pot = 0
  end type gaussian_pot_t

  public :: GaussianPotCreate
  public :: GaussianPotView
  public :: GaussianPotGetNumProp
  public :: GaussianPotGetSymmetry
  public :: GaussianPotGetIntegral
  public :: GaussianPotDestroy

  contains

  !> \brief initializes the information of Gaussian charge potential integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param idx_gauorg contains the atomic centers of Gaussian charge potential origins
  !>        (<1 for non-atomic center)
  !> \param gaupot_origin contains the coordinates of Gaussian charge potential origins
  !> \param gaupot_charge contains the charges of Gaussian charge potential
  !> \param gaupot_expt contains the exponents used in the Gaussian broadening function of charges
  !> \param order_geo_pot is the orders of geometric derivatives with respect to
  !>         the Gaussian charge potential origins
  !> \return gaussian_pot contains the information of Gaussian charge potential integrals
  !> \return info_prop (==ERR_INFO) indicates the Gaussian charge potential integrals
  !>         are not successfully created
  subroutine GaussianPotCreate(gaussian_pot, info_prop, idx_gauorg, gaupot_origin, &
                               gaupot_charge, gaupot_expt, order_geo_pot)
    type(gaussian_pot_t), intent(inout) :: gaussian_pot
    integer, intent(out) :: info_prop
    integer, optional, intent(in) :: idx_gauorg(:)
    real(REALK), optional, intent(in) :: gaupot_origin(:,:)
    real(REALK), optional, intent(in) :: gaupot_charge(:)
    real(REALK), optional, intent(in) :: gaupot_expt(:)
    integer, optional, intent(in) :: order_geo_pot
    integer icent          !incremental recorder
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (present(order_geo_pot)) gaussian_pot%order_geo_pot = order_geo_pot
    ! negative order, returns with error code
    if (gaussian_pot%order_geo_pot<0) then
      info_prop = gaussian_pot%order_geo_pot
    else if (present(gaupot_origin) .and. present(gaupot_charge) .and. &
             present(gaupot_expt)) then
      gaussian_pot%num_origin = size(gaupot_origin,2)
      if (size(gaupot_origin,1)==3 .and.                     &
          gaussian_pot%num_origin==size(gaupot_charge) .and. &
          gaussian_pot%num_origin==size(gaupot_expt)) then
        info_prop = 0
        allocate(gaussian_pot%idx_gauorg(gaussian_pot%num_origin), stat=info_prop)
        if (info_prop/=0) return
        allocate(gaussian_pot%gaupot_origin(3,gaussian_pot%num_origin), stat=info_prop)
        if (info_prop/=0) then
          deallocate(gaussian_pot%idx_gauorg)
          return
        end if
        allocate(gaussian_pot%gaupot_charge(gaussian_pot%num_origin), stat=info_prop)
        if (info_prop/=0) then
          deallocate(gaussian_pot%idx_gauorg)
          deallocate(gaussian_pot%gaupot_origin)
          return
        end if
        allocate(gaussian_pot%gaupot_expt(gaussian_pot%num_origin), stat=info_prop)
        if (info_prop/=0) then
          deallocate(gaussian_pot%idx_gauorg)
          deallocate(gaussian_pot%gaupot_origin)
          deallocate(gaussian_pot%gaupot_charge)
          return
        end if
        if (present(idx_gauorg)) then
          if (size(idx_gauorg)/=gaussian_pot%num_origin) then
            deallocate(gaussian_pot%idx_gauorg)
            deallocate(gaussian_pot%gaupot_origin)
            deallocate(gaussian_pot%gaupot_charge)
            deallocate(gaussian_pot%gaupot_expt)
            info_prop = ERR_INFO
            return
          end if
          gaussian_pot%idx_gauorg = idx_gauorg
        else
          do icent = 1, gaussian_pot%num_origin
            gaussian_pot%idx_gauorg(icent) = icent
          end do
        end if
        gaussian_pot%gaupot_origin = gaupot_origin
        gaussian_pot%gaupot_charge = gaupot_charge
        gaussian_pot%gaupot_expt = gaupot_expt
        gaussian_pot%num_prop = (gaussian_pot%order_geo_pot+1) &
                              * (gaussian_pot%order_geo_pot+2)/2
      else
        info_prop = ERR_INFO
      end if
    else
      info_prop = ERR_INFO
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "GaussianPotCreate", STDOUT)
#endif
  end subroutine GaussianPotCreate

  !> \brief visualizes the information of Gaussian charge potential integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param gaussian_pot contains the information of Gaussian charge potential integrals
  !> \param io_viewer is the logical unit number of the viewer
  subroutine GaussianPotView(gaussian_pot, io_viewer)
    type(gaussian_pot_t), intent(in) :: gaussian_pot
    integer, intent(in) :: io_viewer
    integer icent          !incremental recorder over centers
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    write(io_viewer,100) "Gaussian charge potential integrals"
    write(io_viewer,100) "number of property integral matrices", gaussian_pot%num_prop
    write(io_viewer,100) "symmetry of property integral matrices", gaussian_pot%prop_sym
    write(io_viewer,100) "number of origins", gaussian_pot%num_origin
    write(io_viewer,100) "order of geometric derivatives with respect to origins", &
                         gaussian_pot%order_geo_pot
    write(io_viewer,110) "atomic center", "coordinates", "charge", "exponent"
    do icent = 1, gaussian_pot%num_origin
      write(io_viewer,120) gaussian_pot%idx_gauorg(icent),      &
                           gaussian_pot%gaupot_origin(:,icent), &
                           gaussian_pot%gaupot_charge(icent),   &
                           gaussian_pot%gaupot_expt(icent)
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "GaussianPotView", STDOUT)
#endif
100 format("GaussianPotView>> ",A,I6)
110 format("GaussianPotView>> ",A,10X,A,15X,A)
120 format("GaussianPotView>> ",I6,3F14.8,F8.3,F14.8)
  end subroutine GaussianPotView

  !> \brief returns the number of integral matrices for given Gaussian charge potential integrals
  !> \param gaussian_pot contains the information of one-electron potential integrals
  !> \return num_prop is the number of property integral matrices
  subroutine GaussianPotGetNumProp(gaussian_pot, num_prop)
    type(gaussian_pot_t), intent(in) :: gaussian_pot
    integer, intent(out) :: num_prop
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif 
    num_prop = gaussian_pot%num_prop
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "GaussianPotGetNumProp", STDOUT)
#endif
  end subroutine GaussianPotGetNumProp

  !> \brief returns the symmetry of integral matrices for given Gaussian charge potential integrals
  !> \param gaussian_pot contains the information of one-electron potential integrals
  !> \return prop_sym indicates the symmetry of property integral matrices
  subroutine GaussianPotGetSymmetry(gaussian_pot, prop_sym)
    type(gaussian_pot_t), intent(in) :: gaussian_pot
    integer, intent(out) :: prop_sym
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif 
    prop_sym = gaussian_pot%prop_sym
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "GaussianPotGetSymmetry", STDOUT)
#endif
  end subroutine GaussianPotGetSymmetry

  !> \brief evaluates the Gaussian charge potential integrals
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
  !> \param gaussian_pot contains the information of Gaussian charge potential integrals
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
  subroutine GaussianPotGetIntegral(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                                    exponent_bra, num_contr_bra, contr_coef_bra,   &
                                    idx_ket, coord_ket, angular_ket, num_prim_ket, &
                                    exponent_ket, num_contr_ket, contr_coef_ket,   &
                                    spher_gto, info_LAO, gaussian_pot,             &
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
    type(gaussian_pot_t), intent(in) :: gaussian_pot
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
    call IntGetGAUPOT(idx_bra, coord_bra, angular_bra, num_prim_bra,  &
                      exponent_bra, num_contr_bra, contr_coef_bra,    &
                      idx_ket, coord_ket, angular_ket, num_prim_ket,  &
                      exponent_ket, num_contr_ket, contr_coef_ket,    &
                      spher_gto, info_LAO, 0,                         &
                      gaussian_pot%idx_gauorg(1),                     &
                      gaussian_pot%gaupot_origin(:,1),                &
                      gaussian_pot%gaupot_expt(1), MAX_IDX_NON,       &
                      (/0.0_REALK,0.0_REALK,0.0_REALK/),              &
                      gaussian_pot%gaupot_charge(1), 0,               &
                      order_mag_bra, order_mag_ket, order_mag_total,  &
                      order_ram_bra, order_ram_ket, order_ram_total,  &
                      order_geo_bra, order_geo_ket,                   &
                      gaussian_pot%order_geo_pot, 0, nary_tree_total, &
                      num_gto_bra, num_gto_ket, num_opt, contr_ints,  &
                      mag_num_bra, mag_num_ket, powers_bra, powers_ket)
    allocate(tmp_ints(num_gto_bra,num_contr_bra,num_gto_ket,num_contr_ket, &
                      num_opt), stat=ierr)
    if (ierr/=0)                                                               &
      call error_stop("GaussianPotGetIntegral", "failed to allocate tmp_ints", &
                      num_gto_bra*num_contr_bra*num_gto_ket*num_contr_ket*num_opt)
    ! other atomic centers
    do icent = 2, gaussian_pot%num_origin
      call IntGetGAUPOT(idx_bra, coord_bra, angular_bra, num_prim_bra,  &
                        exponent_bra, num_contr_bra, contr_coef_bra,    &
                        idx_ket, coord_ket, angular_ket, num_prim_ket,  &
                        exponent_ket, num_contr_ket, contr_coef_ket,    &
                        spher_gto, info_LAO, 0,                         &
                        gaussian_pot%idx_gauorg(icent),                 &
                        gaussian_pot%gaupot_origin(:,icent),            &
                        gaussian_pot%gaupot_expt(icent), MAX_IDX_NON,   &
                        (/0.0_REALK,0.0_REALK,0.0_REALK/),              &
                        gaussian_pot%gaupot_charge(icent), 0,           &
                        order_mag_bra, order_mag_ket, order_mag_total,  &
                        order_ram_bra, order_ram_ket, order_ram_total,  &
                        order_geo_bra, order_geo_ket,                   &
                        gaussian_pot%order_geo_pot, 0, nary_tree_total, &
                        num_gto_bra, num_gto_ket, num_opt, tmp_ints,    &
                        mag_num_bra, mag_num_ket, powers_bra, powers_ket)
      contr_ints = contr_ints+tmp_ints
    end do
    deallocate(tmp_ints)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "GaussianPotGetIntegral", STDOUT)
#endif
  end subroutine GaussianPotGetIntegral

  !> \brief frees space taken by the Gaussian charge potential integrals
  !> \author Bin Gao
  !> \date 2011-12-28
  !> \param gaussian_pot contains the information of Gaussian charge potential integrals
  subroutine GaussianPotDestroy(gaussian_pot)
    type(gaussian_pot_t), intent(inout) :: gaussian_pot
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    gaussian_pot%num_origin = 0
    deallocate(gaussian_pot%idx_gauorg)
    deallocate(gaussian_pot%gaupot_origin)
    deallocate(gaussian_pot%gaupot_charge)
    deallocate(gaussian_pot%gaupot_expt)
    gaussian_pot%order_geo_pot = 0
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "GaussianPotDestroy", STDOUT)
#endif
  end subroutine GaussianPotDestroy

end module gen1int_gaupot
