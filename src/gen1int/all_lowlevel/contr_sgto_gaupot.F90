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
!!  This file calculates the Gaussian charge potential (with Cartesian
!!  multipole moment) integrals using contracted spherical Gaussians.
!!
!!  2011-12-08, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief calculates the Gaussian charge potential (with Cartesian multipole
  !>        moment) integrals using contracted spherical Gaussians
  !> \author Bin Gao
  !> \date 2011-12-08
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
  !> \param order_elec is the order of electronic derivatives
  !> \param idx_gauorg is the atomic center of Gaussian charge potential origin (<1 for non-atomic center)
  !> \param gaupot_origin contains the coordinates of Gaussian charge potential origin
  !> \param gaupot_expt is the exponent used in the Gaussian broadening function of the charge
  !> \param idx_diporg is the atomic center of dipole origin (<1 for non-atomic center)
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scale constant for Gaussian charge potential
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_geo_bra is the order of geometric derivatives with respect to bra center
  !> \param order_geo_ket is the order of geometric derivatives with respect to ket center
  !> \param order_geo_pot is the order of geometric derivatives on potential origin
  !> \param order_geo_mom is the order of geometric derivatives on dipole origin
  !> \param num_cents is the number of differentiated centers
  !> \param idx_cent contains the indices of differentiated centers
  !> \param order_cent contains the order of derivatives of the differentiated centers
  !> \param num_sgto_bra is the number of spherical GTOs on bra center,
  !>        equals to \f$2\var(angular_bra)+1\f$
  !> \param num_sgto_ket is the number of spherical GTOs on ket center,
  !>        equals to \f$2\var(angular_ket)+1\f$
  !> \param num_opt is the number of operators including derivatives
  !> \return contr_ints contains the contracted Cartesian multipole moment integrals
  subroutine contr_sgto_gaupot(idx_bra, coord_bra, angular_bra, num_prim_bra,      &
                               exponent_bra, num_contr_bra, contr_coef_bra,        &
                               idx_ket, coord_ket, angular_ket, num_prim_ket,      &
                               exponent_ket, num_contr_ket, contr_coef_ket,        &
                               order_elec, idx_gauorg, gaupot_origin, gaupot_expt, &
                               idx_diporg, dipole_origin, scal_const,              &
                               order_mom, order_geo_bra, order_geo_ket,            &
                               order_geo_pot, order_geo_mom, num_cents,            &
                               idx_cent, order_cent, num_sgto_bra,                 &
                               num_sgto_ket, num_opt, contr_ints)
    use xkind
    implicit none
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
    integer, intent(in) :: order_elec
    integer, intent(in) :: idx_gauorg
    real(REALK), intent(in) :: gaupot_origin(3)
    real(REALK), intent(in) :: gaupot_expt
    integer, intent(in) :: idx_diporg
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_geo_bra
    integer, intent(in) :: order_geo_ket
    integer, intent(in) :: order_geo_pot
    integer, intent(in) :: order_geo_mom
    integer, intent(in) :: num_cents
    integer, intent(in) :: idx_cent(num_cents)
    integer, intent(in) :: order_cent(num_cents)
    integer, intent(in) :: num_sgto_bra
    integer, intent(in) :: num_sgto_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(out) :: contr_ints(num_sgto_bra,num_contr_bra, &
                                           num_sgto_ket,num_contr_ket,num_opt)
!f2py intent(in) :: idx_bra
!f2py intent(in) :: coord_bra
!f2py intent(in) :: angular_bra
!f2py intent(hide) :: num_prim_bra
!f2py intent(in) :: exponent_bra
!f2py intent(hide) :: num_contr_bra
!f2py intent(in) :: contr_coef_bra
!f2py depend(num_prim_bra) :: contr_coef_bra
!f2py intent(in) :: idx_ket
!f2py intent(in) :: coord_ket
!f2py intent(in) :: angular_ket
!f2py intent(hide) :: num_prim_ket
!f2py intent(in) :: exponent_ket
!f2py intent(hide) :: num_contr_ket
!f2py intent(in) :: contr_coef_ket
!f2py depend(num_prim_ket) :: contr_coef_ket
!f2py intent(in) :: order_elec
!f2py intent(in) :: idx_gauorg
!f2py intent(in) :: gaupot_origin
!f2py intent(in) :: gaupot_expt
!f2py intent(in) :: idx_diporg
!f2py intent(in) :: dipole_origin
!f2py intent(in) :: scal_const
!f2py intent(in) :: order_mom
!f2py intent(in) :: order_geo_bra
!f2py intent(in) :: order_geo_ket
!f2py intent(in) :: order_geo_pot
!f2py intent(in) :: order_geo_mom
!f2py intent(hide) :: num_cents
!f2py intent(in) :: idx_cent
!f2py intent(in) :: order_cent
!f2py depend(num_cents) :: order_cent
!f2py intent(in) :: num_sgto_bra
!f2py intent(in) :: num_sgto_ket
!f2py intent(in) :: num_opt
!f2py intent(out) :: contr_ints
!f2py depend(num_sgto_bra) :: contr_ints
!f2py depend(num_contr_bra) :: contr_ints
!f2py depend(num_sgto_ket) :: contr_ints
!f2py depend(num_contr_ket) :: contr_ints
!f2py depend(num_opt) :: contr_ints
    integer idx_part_cent(4)    !indices of bra, ket and operator centers
    integer order_part_cent(4)  !orders of partial geometric derivatives of bra, ket and operator centers
    logical zero_ints           !indicates if the integrals are zero
    logical neg_one             !indicates if the integrals will be multiplied by -1
    logical scatter_deriv       !indicates if scattering the geometric derivatives
    real(REALK) neg_scal_const  !scale constant by considering \var(neg_one) for Cartesian multipole moments
    integer seq_part_geo(4)     !sequence of operator, bra and ket centers for partial geometric derivatives
    integer dim_cints           !dimension of contracted integrals
    integer num_elec            !number of xyz components of electronic derivatives
    integer num_mom             !number of xyz components of Cartesian multipole moments
    integer num_elec_mom        !number of electronic derivatives and Cartesian multipole moments
    integer dim_geo_bra         !dimension of geometric derivatives on bra center before scattering
    integer dim_geo_ket         !dimension of geometric derivatives on ket center before scattering
    integer dim_geo_pot         !dimension of geometric derivatives on operator centers before scattering
    integer dim_geo_mom
    integer num_part_opt        !number of operators and derivatives before scattering
    integer num_geo_bra         !number of partial geometric derivatives on bra center after scattering
    integer num_geo_ket         !number of partial geometric derivatives on ket center after scattering
    integer num_geo_pot         !number of partial geometric derivatives on operator centers after scattering
    integer num_geo_mom
    integer num_tot_geo         !number of total geometric derivatives after scattering
    real(REALK), allocatable :: part_cints(:,:,:,:,:)
                                !contains the derivatives of integrals before scattering
    integer ierr                !error information
#if defined(XTIME)
    real(REALK) curr_time       !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
#if defined(DEBUG)
    ! dumps the contracted GTOs and derivatives to check
    call dump_gto_pd("SGTO-bra", idx_bra, coord_bra, angular_bra, &
                     num_prim_bra, exponent_bra, num_contr_bra,   &
                     contr_coef_bra, 0, 0, order_geo_bra, STDOUT)
    call dump_gto_pd("SGTO-ket", idx_ket, coord_ket, angular_ket, &
                     num_prim_ket, exponent_ket, num_contr_ket,   &
                     contr_coef_ket, 0, 0, order_geo_ket, STDOUT)
    call dump_total_derv(0, 0, num_cents, idx_cent, order_cent, STDOUT)
    ! dumps the operator to check
    call dump_gaupot(order_elec, idx_gauorg, gaupot_origin, gaupot_expt, &
                     idx_diporg, dipole_origin, scal_const, order_mom,   &   
                     order_geo_pot, order_geo_mom, STDOUT)
#endif
    ! no geometric derivatives
    if (num_cents==0) then
      call contr_sgto_gaupot_recurr(coord_bra, angular_bra, num_prim_bra,        &
                                    exponent_bra, num_contr_bra, contr_coef_bra, &
                                    coord_ket, angular_ket, num_prim_ket,        &
                                    exponent_ket, num_contr_ket, contr_coef_ket, &
                                    order_elec, gaupot_origin, gaupot_expt,      &
                                    dipole_origin, scal_const, order_mom,        &
                                    order_geo_bra, order_geo_ket, order_geo_pot, &
                                    order_geo_mom, num_sgto_bra, num_sgto_ket,   &
                                    num_opt, contr_ints)
    else
      ! total geometric derivatives of Gaussian charge potential integrals
      if (order_mom==0) then
        ! sets the orders of partial geometric derivatives, the last is dipole origin
        idx_part_cent(1) = idx_bra
        idx_part_cent(2) = idx_ket
        idx_part_cent(3) = idx_gauorg
        order_part_cent(1) = order_geo_bra
        order_part_cent(2) = order_geo_ket
        order_part_cent(3) = order_geo_pot
        call geom_part_one_param(num_cents, idx_cent, order_cent,          &
                                 idx_part_cent(1:3), order_part_cent(1:3), &
                                 zero_ints, neg_one, scatter_deriv, seq_part_geo)
        ! zero integrals
        if (zero_ints) then
          contr_ints = 0.0_REALK
        ! we need to scatter the partial and total geometric derivatives into appropriate places
        else if (scatter_deriv) then
          ! allocates memory for the derivatives of integrals before scattering
          dim_cints = size(contr_ints)/num_opt
          num_elec = (order_elec+1)*(order_elec+2)/2
          dim_geo_bra = (order_part_cent(1)+1)*(order_part_cent(1)+2)/2
          dim_geo_ket = (order_part_cent(2)+1)*(order_part_cent(2)+2)/2
          dim_geo_pot = (order_part_cent(3)+1)*(order_part_cent(3)+2)/2
          allocate(part_cints(dim_cints,num_elec,dim_geo_bra, &
                              dim_geo_ket,dim_geo_pot), stat=ierr)
          if (ierr/=0)                                                              &
            call error_stop("contr_sgto_gaupot", "failed to allocate part_cints/1", &
                            dim_cints*num_elec*dim_geo_bra*dim_geo_ket*dim_geo_pot)
          num_part_opt = size(part_cints)/dim_cints
          ! sets the scale constant
          if (neg_one) then
            neg_scal_const = -scal_const
          else
            neg_scal_const = scal_const
          end if
          call contr_sgto_gaupot_recurr(coord_bra, angular_bra, num_prim_bra,        &
                                        exponent_bra, num_contr_bra, contr_coef_bra, &
                                        coord_ket, angular_ket, num_prim_ket,        &
                                        exponent_ket, num_contr_ket, contr_coef_ket, &
                                        order_elec, gaupot_origin, gaupot_expt,      &
                                        dipole_origin, neg_scal_const, 0,            &
                                        order_part_cent(1), order_part_cent(2),      &
                                        order_part_cent(3), 0, num_sgto_bra,         &
                                        num_sgto_ket, num_part_opt, part_cints)
          ! scatters the partial and total geometric derivatives
          num_geo_bra = (order_geo_bra+1)*(order_geo_bra+2)/2
          num_geo_ket = (order_geo_ket+1)*(order_geo_ket+2)/2
          num_geo_pot = (order_geo_pot+1)*(order_geo_pot+2)/2
          num_tot_geo = num_opt/(num_elec*num_geo_bra*num_geo_ket*num_geo_pot)
          call geom_part_one_scatter(num_cents, order_cent, seq_part_geo(1:num_cents), &
                                     order_geo_bra, order_geo_ket, order_geo_pot,      &
                                     dim_cints*num_elec, dim_geo_bra, dim_geo_ket,     &
                                     dim_geo_pot, 1, part_cints, num_geo_bra,          &
                                     num_geo_ket, num_geo_pot, num_tot_geo, contr_ints)
          deallocate(part_cints)
        else
          ! sets the scale constant
          if (neg_one) then
            neg_scal_const = -scal_const
          else
            neg_scal_const = scal_const
          end if
          call contr_sgto_gaupot_recurr(coord_bra, angular_bra, num_prim_bra,        &
                                        exponent_bra, num_contr_bra, contr_coef_bra, &
                                        coord_ket, angular_ket, num_prim_ket,        &
                                        exponent_ket, num_contr_ket, contr_coef_ket, &
                                        order_elec, gaupot_origin, gaupot_expt,      &
                                        dipole_origin, neg_scal_const, 0,            &
                                        order_part_cent(1), order_part_cent(2),      &
                                        order_part_cent(3), 0, num_sgto_bra,         &
                                        num_sgto_ket, num_part_opt, contr_ints)
        end if
      ! total geometric derivatives of Gaussian charge potential with Cartesian
      ! multipole moment integrals
      else
        call error_stop("contr_sgto_gaupot",                        &
                        "not implemented! please write to authors", &
                        order_mom)
        !! higher order geometric derivatives vanish for lower order Cartesian multipole moments
        !if (zero_ints .or. order_part_cent(4)>order_mom) then
        !  contr_ints = 0.0_REALK
        !else
        !end if
      end if
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "contr_sgto_gaupot", STDOUT)
#endif
    return
  end subroutine contr_sgto_gaupot

  !> \brief recurrence relations of Gaussian charge potential (with Cartesian multipole
  !>        moment) integrals using contracted spherical Gaussians
  !> \author Bin Gao
  !> \date 2011-09-24
  !> \param coord_bra contains the coordinates of bra center
  !> \param angular_bra is the angular number of bra center (s=0, p=1, d=2, ...)
  !> \param num_prim_bra is the number of primitive Gaussians of bra center
  !> \param exponent_bra contains the exponents of primitive Gaussians of bra center
  !> \param num_contr_bra is the number of contractions of bra center
  !> \param contr_coef_bra contains the contraction coefficients of bra center
  !> \param coord_ket contains the coordinates of ket center
  !> \param angular_ket is the angular number of ket center
  !> \param num_prim_ket is the number of primitive Gaussians of ket center
  !> \param exponent_ket contains the exponents of primitive Gaussians of ket center
  !> \param num_contr_ket is the number of contractions of ket center
  !> \param contr_coef_ket contains the contraction coefficients of ket center
  !> \param order_elec is the order of electronic derivatives
  !> \param gaupot_origin contains the coordinates of Gaussian charge potential origin
  !> \param gaupot_expt is the exponent used in the Gaussian broadening function of the charge
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scale constant for Gaussian charge potential
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_geo_bra is the order of geometric derivatives with respect to bra center
  !> \param order_geo_ket is the order of geometric derivatives with respect to ket center
  !> \param order_geo_pot is the order of geometric derivatives on potential origin
  !> \param order_geo_mom is the order of geometric derivatives on dipole origin
  !> \param num_sgto_bra is the number of spherical GTOs on bra center,
  !>        equals to \f$2\var(angular_bra)+1\f$
  !> \param num_sgto_ket is the number of spherical GTOs on ket center,
  !>        equals to \f$2\var(angular_ket)+1\f$
  !> \param num_opt is the number of operators including derivatives
  !> \return contr_ints contains the contracted Cartesian multipole moment integrals
  subroutine contr_sgto_gaupot_recurr(coord_bra, angular_bra, num_prim_bra,        &
                                      exponent_bra, num_contr_bra, contr_coef_bra, &
                                      coord_ket, angular_ket, num_prim_ket,        &
                                      exponent_ket, num_contr_ket, contr_coef_ket, &
                                      order_elec, gaupot_origin, gaupot_expt,      &
                                      dipole_origin, scal_const, order_mom,        &
                                      order_geo_bra, order_geo_ket, order_geo_pot, &
                                      order_geo_mom, num_sgto_bra, num_sgto_ket,   &
                                      num_opt, contr_ints)
    use xkind
    implicit none
    real(REALK), intent(in) :: coord_bra(3)
    integer, intent(in) :: angular_bra
    integer, intent(in) :: num_prim_bra
    real(REALK), intent(in) :: exponent_bra(num_prim_bra)
    integer, intent(in) :: num_contr_bra
    real(REALK), intent(in) :: contr_coef_bra(num_contr_bra,num_prim_bra)
    real(REALK), intent(in) :: coord_ket(3)
    integer, intent(in) :: angular_ket
    integer, intent(in) :: num_prim_ket
    real(REALK), intent(in) :: exponent_ket(num_prim_ket)
    integer, intent(in) :: num_contr_ket
    real(REALK), intent(in) :: contr_coef_ket(num_contr_ket,num_prim_ket)
    integer, intent(in) :: order_elec
    real(REALK), intent(in) :: gaupot_origin(3)
    real(REALK), intent(in) :: gaupot_expt
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_geo_bra
    integer, intent(in) :: order_geo_ket
    integer, intent(in) :: order_geo_pot
    integer, intent(in) :: order_geo_mom
    integer, intent(in) :: num_sgto_bra
    integer, intent(in) :: num_sgto_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(out) :: contr_ints(num_sgto_bra,num_contr_bra, &
                                           num_sgto_ket,num_contr_ket,num_opt)
!f2py intent(in) :: coord_bra
!f2py intent(in) :: angular_bra
!f2py intent(hide) :: num_prim_bra
!f2py intent(in) :: exponent_bra
!f2py intent(hide) :: num_contr_bra
!f2py intent(in) :: contr_coef_bra
!f2py depend(num_prim_bra) :: contr_coef_bra
!f2py intent(in) :: coord_ket
!f2py intent(in) :: angular_ket
!f2py intent(hide) :: num_prim_ket
!f2py intent(in) :: exponent_ket
!f2py intent(hide) :: num_contr_ket
!f2py intent(in) :: contr_coef_ket
!f2py depend(num_prim_ket) :: contr_coef_ket
!f2py intent(in) :: order_elec
!f2py intent(in) :: gaupot_origin
!f2py intent(in) :: gaupot_expt
!f2py intent(in) :: dipole_origin
!f2py intent(in) :: scal_const
!f2py intent(in) :: order_mom
!f2py intent(in) :: order_geo_bra
!f2py intent(in) :: order_geo_ket
!f2py intent(in) :: order_geo_pot
!f2py intent(in) :: order_geo_mom
!f2py intent(in) :: num_sgto_bra
!f2py intent(in) :: num_sgto_ket
!f2py intent(in) :: num_opt
!f2py intent(out) :: contr_ints
!f2py depend(num_sgto_bra) :: contr_ints
!f2py depend(num_contr_bra) :: contr_ints
!f2py depend(num_sgto_ket) :: contr_ints
!f2py depend(num_contr_ket) :: contr_ints
!f2py depend(num_opt) :: contr_ints
    integer orders_hgto_bra(2)  !range of orders of Hermite Gaussians on bra center
    integer orders_hgto_ket(2)  !range of orders of Hermite Gaussians on ket center
    integer order_low_mom       !order of Cartesian multipole moment by considering its derivative
    integer num_elec            !number of xyz components of electronic derivatives
    integer num_low_mom         !number of xyz components of lower order Cartesian multipole moment
    integer num_elec_mom        !number of xyz components of lower order Cartesian multipole moment
                                !and electronic derivatives
    integer dim_hgto_bra        !dimension of Hermite Gaussians on bra center including geometric derivatives
    integer dim_hgto_ket        !dimension of Hermite Gaussians on ket center including geometric derivatives
    integer num_hgto_bra        !number of Hermite Gaussians on bra center
    integer num_hgto_ket        !number of Hermite Gaussians on ket center
    integer num_geo_bra         !number of geometric derivatives on bra center
    integer num_geo_ket         !number of geometric derivatives on ket center
    integer num_geo_pot         !number of geometric derivatives on potential origin
    integer num_opt_geo         !number of operators and geometric derivatives
    real(REALK), allocatable :: hgto_pints(:,:,:,:)         !primitive Hermite integrals
    real(REALK), allocatable :: kgeo_hgto_pints(:,:,:,:,:)  !geometric derivatives (ket) of primitive
                                                            !Hermite integrals
    real(REALK), allocatable :: geom_hgto_pints(:,:,:,:,:)  !geometric derivatives of primitive
                                                            !Hermite integrals
    real(REALK), allocatable :: geom_hgto_cints(:,:,:,:,:)  !geometric derivatives of contracted
                                                            !Hermite integrals
    real(REALK), allocatable :: sgto_ket_cints(:,:,:,:,:)   !contracted SGTO (on ket center) integrals
    real(REALK), allocatable :: low_mom_cints(:,:,:,:,:)    !contracted integrals of lower order
                                                            !Cartesian multipole moment
    real(REALK), allocatable :: scal_geo_bra(:)             !scale constants when recovering
                                                            !geometric derivatives on bra center
    real(REALK) scal_geo_ket                                !scale constant when recovering
                                                            !geometric derivatives on ket center
    integer iprim, jprim        !incremental recorders over primitives
    integer ierr                !error information
#if defined(XTIME)
    real(REALK) curr_time       !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! sets the orders of Hermite Gaussians on bra and ket centers
    orders_hgto_bra(1) = order_geo_bra+angular_bra
    orders_hgto_bra(2) = orders_hgto_bra(1)
    orders_hgto_ket(1) = order_geo_ket+angular_ket
    orders_hgto_ket(2) = orders_hgto_ket(1)
    ! computes the order of Cartesian multipole moment after considering its geometric derivative
    order_low_mom = order_mom-order_geo_mom
    ! higher order geometric derivatives vanish for lower order Cartesian multipole moments
    if (order_low_mom<0) then
      contr_ints = 0.0_REALK
      return
    end if
    ! sets the dimension of primitive Hermite Gaussians including geometric derivatives
    dim_hgto_bra = (orders_hgto_bra(1)+1)*(orders_hgto_bra(1)+2)/2
    dim_hgto_ket = (orders_hgto_ket(1)+1)*(orders_hgto_ket(1)+2)/2
    ! computes the number of primitive Hermite Gaussians
    num_hgto_bra = (angular_bra+1)*(angular_bra+2)/2
    num_hgto_ket = (angular_ket+1)*(angular_ket+2)/2
    ! computes the number of geometric derivatives
    num_geo_bra = (order_geo_bra+1)*(order_geo_bra+2)/2
    num_geo_ket = (order_geo_ket+1)*(order_geo_ket+2)/2
    num_geo_pot = (order_geo_pot+1)*(order_geo_pot+2)/2
    ! computes the number of operators, which is more efficient than using \var(TRI_SIZE)
    num_elec = (order_elec+1)*(order_elec+2)/2
    num_low_mom = (order_low_mom+1)*(order_low_mom+2)/2
    num_elec_mom = num_elec*num_low_mom
    num_opt_geo = num_elec_mom*num_geo_bra*num_geo_ket*num_geo_pot
    ! allocates the memory for primitive Hermite integrals
    allocate(hgto_pints(dim_hgto_bra,dim_hgto_ket,num_elec_mom,num_geo_pot), stat=ierr)
    if (ierr/=0)                                       &
      call error_stop("contr_sgto_gaupot_recurr",      &
                      "failed to allocate hgto_pints", &
                      dim_hgto_bra*dim_hgto_ket*num_elec_mom*num_geo_pot)
    allocate(kgeo_hgto_pints(dim_hgto_bra,num_hgto_ket, &
                            num_elec_mom,num_geo_ket,num_geo_pot), stat=ierr)
    if (ierr/=0)                                            &
      call error_stop("contr_sgto_gaupot_recurr",           &
                      "failed to allocate kgeo_hgto_pints", &
                      dim_hgto_bra*num_hgto_ket*num_elec_mom*num_geo_ket*num_geo_pot)
    allocate(geom_hgto_pints(num_hgto_bra,num_hgto_ket,num_opt_geo, &
                        num_prim_bra,num_prim_ket), stat=ierr)
    if (ierr/=0)                                            &
      call error_stop("contr_sgto_gaupot_recurr",           &
                      "failed to allocate geom_hgto_pints", &
                      num_hgto_bra*num_hgto_ket*num_opt_geo*num_prim_bra*num_prim_ket)
#if defined(DEBUG)
    write(STDOUT,100) "size of primitive Hermite integrals:", size(hgto_pints)
    write(STDOUT,100) "size of geometric derivatives of primitive "// &
                      "Hermite (ket) integrals:", size(kgeo_hgto_pints)
    write(STDOUT,100) "size of geometric derivatives of primitive "// &
                      "Hermite integrals:", size(geom_hgto_pints)
#endif
    ! sets the scale constant when recovering geometric derivatives on bra center
    if (order_geo_bra>0) then
      allocate(scal_geo_bra(num_prim_bra), stat=ierr)
      if (ierr/=0)                                  &
        call error_stop("contr_sgto_gaupot_recurr", &
                        "failed to allocate scal_geo_bra", num_prim_bra)
      do iprim = 1, num_prim_bra
        scal_geo_bra(iprim) = (exponent_bra(iprim)+exponent_bra(iprim))**order_geo_bra
      end do
      if (order_geo_ket>0) then
        ! computes the primitive integrals
        do jprim = 1, num_prim_ket
          ! sets the scale constant when recovering geometric derivatives on ket center
          scal_geo_ket = (exponent_ket(jprim)+exponent_ket(jprim))**order_geo_ket
          do iprim = 1, num_prim_bra
            ! computes the primitive Hermite integrals
            call prim_hgto_gaupot(orders_hgto_bra, coord_bra, exponent_bra(iprim), &
                                  orders_hgto_ket, coord_ket, exponent_ket(jprim), &
                                  order_elec, gaupot_origin, gaupot_expt,          &
                                  dipole_origin, scal_const, order_low_mom,        &
                                  order_geo_pot, dim_hgto_bra, dim_hgto_ket,       &
                                  num_elec, num_low_mom, num_geo_pot, hgto_pints)
            ! scales the primitive Hermite integrals
            hgto_pints = scal_geo_ket*hgto_pints
            hgto_pints = scal_geo_bra(iprim)*hgto_pints
            ! recovers the geometric derivatives on ket center
            call scatter_single(dim_hgto_bra, dim_hgto_ket, num_elec_mom, num_geo_pot, &
                                hgto_pints, angular_ket, order_geo_ket, num_hgto_ket,  &
                                num_geo_ket, kgeo_hgto_pints)
            ! recovers the geometric derivatives on bra center
            call scatter_single(1, dim_hgto_bra, num_hgto_ket*num_elec_mom, &
                                num_geo_ket*num_geo_pot, kgeo_hgto_pints,   &
                                angular_bra, order_geo_bra, num_hgto_bra,   &
                                num_geo_bra, geom_hgto_pints(:,:,:,iprim,jprim))
          end do
        end do
      else
        ! computes the primitive integrals
        do jprim = 1, num_prim_ket
          do iprim = 1, num_prim_bra
            ! computes the primitive Hermite integrals
            call prim_hgto_gaupot(orders_hgto_bra, coord_bra, exponent_bra(iprim), &
                                  orders_hgto_ket, coord_ket, exponent_ket(jprim), &
                                  order_elec, gaupot_origin, gaupot_expt,          &
                                  dipole_origin, scal_const, order_low_mom,        &
                                  order_geo_pot, dim_hgto_bra, dim_hgto_ket,       &
                                  num_elec, num_low_mom, num_geo_pot, hgto_pints)
            ! scales the primitive Hermite integrals
            hgto_pints = scal_geo_bra(iprim)*hgto_pints
            ! recovers the geometric derivatives on ket center
            call scatter_single(dim_hgto_bra, dim_hgto_ket, num_elec_mom, num_geo_pot, &
                                hgto_pints, angular_ket, order_geo_ket, num_hgto_ket,  &
                                num_geo_ket, kgeo_hgto_pints)
            ! recovers the geometric derivatives on bra center
            call scatter_single(1, dim_hgto_bra, num_hgto_ket*num_elec_mom, &
                                num_geo_ket*num_geo_pot, kgeo_hgto_pints,   &
                                angular_bra, order_geo_bra, num_hgto_bra,   &
                                num_geo_bra, geom_hgto_pints(:,:,:,iprim,jprim))
          end do
        end do
      end if
      deallocate(scal_geo_bra)
    else if (order_geo_ket>0) then
      ! computes the primitive integrals
      do jprim = 1, num_prim_ket
        ! sets the scale constant when recovering geometric derivatives on ket center
        scal_geo_ket = (exponent_ket(jprim)+exponent_ket(jprim))**order_geo_ket
        do iprim = 1, num_prim_bra
          ! computes the primitive Hermite integrals
          call prim_hgto_gaupot(orders_hgto_bra, coord_bra, exponent_bra(iprim), &
                                orders_hgto_ket, coord_ket, exponent_ket(jprim), &
                                order_elec, gaupot_origin, gaupot_expt,          &
                                dipole_origin, scal_const, order_low_mom,        &
                                order_geo_pot, dim_hgto_bra, dim_hgto_ket,       &
                                num_elec, num_low_mom, num_geo_pot, hgto_pints)
          ! scales the primitive Hermite integrals
          hgto_pints = scal_geo_ket*hgto_pints
          ! recovers the geometric derivatives on ket center
          call scatter_single(dim_hgto_bra, dim_hgto_ket, num_elec_mom, num_geo_pot, &
                              hgto_pints, angular_ket, order_geo_ket, num_hgto_ket,  &
                              num_geo_ket, kgeo_hgto_pints)
          ! recovers the geometric derivatives on bra center
          call scatter_single(1, dim_hgto_bra, num_hgto_ket*num_elec_mom, &
                              num_geo_ket*num_geo_pot, kgeo_hgto_pints,   &
                              angular_bra, order_geo_bra, num_hgto_bra,   &
                              num_geo_bra, geom_hgto_pints(:,:,:,iprim,jprim))
        end do
      end do
    else
      ! computes the primitive integrals
      do jprim = 1, num_prim_ket
        do iprim = 1, num_prim_bra
          ! computes the primitive Hermite integrals
          call prim_hgto_gaupot(orders_hgto_bra, coord_bra, exponent_bra(iprim), &
                                orders_hgto_ket, coord_ket, exponent_ket(jprim), &
                                order_elec, gaupot_origin, gaupot_expt,          &
                                dipole_origin, scal_const, order_low_mom,        &
                                order_geo_pot, dim_hgto_bra, dim_hgto_ket,       &
                                num_elec, num_low_mom, num_geo_pot, hgto_pints)
          ! recovers the geometric derivatives on ket center
          call scatter_single(dim_hgto_bra, dim_hgto_ket, num_elec_mom, num_geo_pot, &
                              hgto_pints, angular_ket, order_geo_ket, num_hgto_ket,  &
                              num_geo_ket, kgeo_hgto_pints)
          ! recovers the geometric derivatives on bra center
          call scatter_single(1, dim_hgto_bra, num_hgto_ket*num_elec_mom, &
                              num_geo_ket*num_geo_pot, kgeo_hgto_pints,   &
                              angular_bra, order_geo_bra, num_hgto_bra,   &
                              num_geo_bra, geom_hgto_pints(:,:,:,iprim,jprim))
        end do
      end do
    end if
    ! cleans
    deallocate(hgto_pints)
    deallocate(kgeo_hgto_pints)
    ! allocates the memory for geometric derivatives of contracted Hermite integrals
    allocate(geom_hgto_cints(num_hgto_bra,num_contr_bra, &
                             num_hgto_ket,num_contr_ket,num_opt_geo), stat=ierr)
    if (ierr/=0)                                            &
      call error_stop("contr_sgto_gaupot_recurr",           &
                      "failed to allocate geom_hgto_cints", &
                      num_hgto_bra*num_contr_bra*num_hgto_ket*num_contr_ket*num_opt_geo)
    ! allocates the memory for contracted SGTO (on ket center) integrals
    allocate(sgto_ket_cints(num_hgto_bra,num_contr_bra, &
                             num_sgto_ket,num_contr_ket,num_opt_geo), stat=ierr)
    if (ierr/=0)                                           &
      call error_stop("contr_sgto_gaupot_recurr",          &
                      "failed to allocate sgto_ket_cints", &
                      num_hgto_bra*num_contr_bra*num_sgto_ket*num_contr_ket*num_opt_geo)
#if defined(DEBUG)
    write(STDOUT,100) "size of geometric derivatives of contracted "// &
                      "Hermite integrals:", size(geom_hgto_cints)
    write(STDOUT,100) "size of contracted SGTO (on ket center) integrals:", &
                      size(sgto_ket_cints)
#endif
    ! constructs contracted Hermite integrals
    ! in the order of geom_hgto_cints(num_hgto_bra,num_contr_bra,num_hgto_ket,num_contr_ket, &
    !                            num_elec,num_low_mom,num_geo_bra,num_geo_ket,num_geo_pot)
    call const_contr_ints(num_contr_bra, num_prim_bra, contr_coef_bra, &
                          num_contr_ket, num_prim_ket, contr_coef_ket, &
                          num_hgto_bra, num_hgto_ket, num_opt_geo,     &
                          geom_hgto_pints, geom_hgto_cints)
    deallocate(geom_hgto_pints)
    ! transforms the HGTOs to SGTOs on ket center
    call hgto_to_sgto_oop(angular_ket, num_hgto_bra*num_contr_bra, &
                          num_hgto_ket, num_contr_ket*num_opt_geo, &
                          geom_hgto_cints, num_sgto_ket, sgto_ket_cints)
    deallocate(geom_hgto_cints)
    ! transforms the HGTOs to SGTOs on bra center
    if (order_geo_mom==0) then
      call hgto_to_sgto_oop(angular_bra, 1, num_hgto_bra,                         &
                            num_contr_bra*num_sgto_ket*num_contr_ket*num_opt_geo, &
                            sgto_ket_cints, num_sgto_bra, contr_ints)
      deallocate(sgto_ket_cints)
    else
      ! transforms to contracted spherical integrals of lower order Cartesian multipole memonts,
      ! in the order of low_mom_cints(num_sgto_bra,num_contr_bra,num_sgto_ket,num_contr_ket, &
      !                               num_elec,num_low_mom,num_geo_bra,num_geo_ket,num_geo_pot)
      allocate(low_mom_cints(num_sgto_bra,num_contr_bra, &
                             num_sgto_ket,num_contr_ket,num_opt_geo), stat=ierr)
      if (ierr/=0)                                          &
        call error_stop("contr_sgto_gaupot_recurr",         &
                        "failed to allocate low_mom_cints", &
                        num_sgto_bra*num_contr_bra          &
                        *num_sgto_ket*num_contr_ket*num_opt_geo)
#if defined(DEBUG)
      write(STDOUT,100)                                                               &
        "number of contracted integrals of lower order Cartesian multipole memonts:", &
        size(low_mom_cints)
#endif
      ! transforms the HGTOs to SGTOs on bra center
      call hgto_to_sgto_oop(angular_bra, 1, num_hgto_bra,                         &
                            num_contr_bra*num_sgto_ket*num_contr_ket*num_opt_geo, &
                            sgto_ket_cints, num_sgto_bra, low_mom_cints)
      deallocate(sgto_ket_cints)
      ! constructs the geometric derivatives of Cartesian multipole moments
      ! in the order of contr_ints(num_sgto_bra,num_contr_bra,num_sgto_ket,num_contr_ket, &
      !                            num_elec,num_mom,num_geo_bra,num_geo_ket,num_geo_pot,num_geo_mom)
      call carmom_deriv(order_low_mom, size(low_mom_cints)/num_opt_geo*num_elec, &
                        num_low_mom, num_geo_bra*num_geo_ket*num_geo_pot,        &
                        low_mom_cints, order_mom, order_geo_mom,                 &
                        (order_mom+1)*(order_mom+2)/2,                           &
                        (order_geo_mom+1)*(order_geo_mom+2)/2, contr_ints)
      ! cleans
      deallocate(low_mom_cints)
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "contr_sgto_gaupot_recurr", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("contr_sgto_gaupot_recurr>> ",A,I8,1X,A)
#endif
  end subroutine contr_sgto_gaupot_recurr
