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
!!  This file calculates the Dirac delta function (with Cartesian multipole
!!  moment) integrals using contracted Cartesian Gaussians.
!!
!!  2011-09-24, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief calculates the Dirac delta function (with Cartesian multipole
  !>        moment) integrals using contracted Cartesian Gaussians
  !> \author Bin Gao
  !> \date 2011-09-24
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
  !> \param idx_delta is the atomic center of delta function origin (<1 for non-atomic center)
  !> \param delta_origin contains the coordinates of delta function origin
  !> \param idx_diporg is the atomic center of dipole origin (<1 for non-atomic center)
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scale constant for Dirac delta function
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_geo_bra is the order of geometric derivatives with respect to bra center
  !> \param order_geo_ket is the order of geometric derivatives with respect to ket center
  !> \param order_geo_pot is the order of geometric derivatives on delta function origin
  !> \param order_geo_mom is the order of geometric derivatives on dipole origin
  !> \param num_cents is the number of differentiated centers
  !> \param idx_cent contains the indices of differentiated centers
  !> \param order_cent contains the order of derivatives of the differentiated centers
  !> \param num_cgto_bra is the number of Cartesian GTOs on bra center,
  !>        equals to \f$(\var(angular_bra)+1)(\var(angular_bra)+2)/2\f$
  !> \param num_cgto_ket is the number of Cartesian GTOs on ket center,
  !>        equals to \f$(\var(angular_ket)+1)(\var(angular_ket)+2)/2\f$
  !> \param num_opt is the number of operators including derivatives
  !> \return contr_ints contains the contracted Cartesian multipole moment integrals
  subroutine contr_cgto_delta(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                              exponent_bra, num_contr_bra, contr_coef_bra,   &
                              idx_ket, coord_ket, angular_ket, num_prim_ket, &
                              exponent_ket, num_contr_ket, contr_coef_ket,   &
                              order_elec, idx_delta, delta_origin,           &
                              idx_diporg, dipole_origin, scal_const,         &
                              order_mom, order_geo_bra, order_geo_ket,       &
                              order_geo_pot, order_geo_mom, num_cents,       &
                              idx_cent, order_cent, num_cgto_bra,            &
                              num_cgto_ket, num_opt, contr_ints)
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
    integer, intent(in) :: idx_delta
    real(REALK), intent(in) :: delta_origin(3)
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
    integer, intent(in) :: num_cgto_bra
    integer, intent(in) :: num_cgto_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(out) :: contr_ints(num_cgto_bra,num_contr_bra, &
                                           num_cgto_ket,num_contr_ket,num_opt)
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
!f2py intent(in) :: idx_delta
!f2py intent(in) :: delta_origin
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
!f2py intent(in) :: num_cgto_bra
!f2py intent(in) :: num_cgto_ket
!f2py intent(in) :: num_opt
!f2py intent(out) :: contr_ints
!f2py depend(num_cgto_bra) :: contr_ints
!f2py depend(num_contr_bra) :: contr_ints
!f2py depend(num_cgto_ket) :: contr_ints
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
    call dump_gto_pd("CGTO-bra", idx_bra, coord_bra, angular_bra, &
                     num_prim_bra, exponent_bra, num_contr_bra,   &
                     contr_coef_bra, 0, 0, order_geo_bra, STDOUT)
    call dump_gto_pd("CGTO-ket", idx_ket, coord_ket, angular_ket, &
                     num_prim_ket, exponent_ket, num_contr_ket,   &
                     contr_coef_ket, 0, 0, order_geo_ket, STDOUT)
    call dump_total_derv(0, 0, num_cents, idx_cent, order_cent, STDOUT)
    ! dumps the operator to check
    call dump_delta(order_elec, idx_delta, delta_origin, idx_diporg,     &
                    dipole_origin, scal_const, order_mom, order_geo_pot, &
                    order_geo_mom, STDOUT)
#endif
    ! integrals will be zero if the origin of Dirac delta function is the same as
    ! (1) dipole origin, (2) bra or ket center with angular number >0 
    if (idx_delta==idx_diporg .or. (idx_delta==idx_bra .and. angular_bra>0) .or. &
        (idx_delta==idx_ket .and. angular_ket>0)) then
      contr_ints = 0.0_REALK
    else
      ! no geometric derivatives
      if (num_cents==0) then
        call contr_cgto_delta_recurr(coord_bra, angular_bra, num_prim_bra,        &
                                     exponent_bra, num_contr_bra, contr_coef_bra, &
                                     coord_ket, angular_ket, num_prim_ket,        &
                                     exponent_ket, num_contr_ket, contr_coef_ket, &
                                     order_elec, delta_origin, dipole_origin,     &
                                     scal_const, order_mom, order_geo_bra,        &
                                     order_geo_ket, order_geo_pot, order_geo_mom, &
                                     num_cgto_bra, num_cgto_ket, num_opt, contr_ints)
      else
        ! total geometric derivatives of Dirac delta function integrals
        if (order_mom==0) then
          ! sets the orders of partial geometric derivatives, the last is delta function origin
          idx_part_cent(1) = idx_bra
          idx_part_cent(2) = idx_ket
          idx_part_cent(3) = idx_delta
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
            if (ierr/=0)                                                             &
              call error_stop("contr_cgto_delta", "failed to allocate part_cints/1", &
                              dim_cints*num_elec*dim_geo_bra*dim_geo_ket*dim_geo_pot)
            num_part_opt = size(part_cints)/dim_cints
            ! sets the scale constant
            if (neg_one) then
              neg_scal_const = -scal_const
            else
              neg_scal_const = scal_const
            end if
            call contr_cgto_delta_recurr(coord_bra, angular_bra, num_prim_bra,        &
                                         exponent_bra, num_contr_bra, contr_coef_bra, &
                                         coord_ket, angular_ket, num_prim_ket,        &
                                         exponent_ket, num_contr_ket, contr_coef_ket, &
                                         order_elec, delta_origin, dipole_origin,     &
                                         neg_scal_const, 0, order_part_cent(1),       &
                                         order_part_cent(2), order_part_cent(3), 0,   &
                                         num_cgto_bra, num_cgto_ket, num_part_opt,    &
                                         part_cints)
            ! scatters the partial and total geometric derivatives
            num_geo_bra = (order_geo_bra+1)*(order_geo_bra+2)/2
            num_geo_ket = (order_geo_ket+1)*(order_geo_ket+2)/2
            num_geo_pot = (order_geo_pot+1)*(order_geo_pot+2)/2
            num_tot_geo = num_opt/(num_elec*num_geo_bra*num_geo_ket*num_geo_pot)
            call geom_part_one_scatter(num_cents, order_cent,             &
                                       seq_part_geo(1:num_cents),         &
                                       order_geo_bra, order_geo_ket,      &
                                       order_geo_pot, dim_cints*num_elec, &
                                       dim_geo_bra, dim_geo_ket,          &
                                       dim_geo_pot, 1, part_cints,        &
                                       num_geo_bra, num_geo_ket,          &
                                       num_geo_pot, num_tot_geo, contr_ints)
            deallocate(part_cints)
          else
            ! sets the scale constant
            if (neg_one) then
              neg_scal_const = -scal_const
            else
              neg_scal_const = scal_const
            end if
            call contr_cgto_delta_recurr(coord_bra, angular_bra, num_prim_bra,        &
                                         exponent_bra, num_contr_bra, contr_coef_bra, &
                                         coord_ket, angular_ket, num_prim_ket,        &
                                         exponent_ket, num_contr_ket, contr_coef_ket, &
                                         order_elec, delta_origin, dipole_origin,     &
                                         neg_scal_const, 0, order_part_cent(1),       &
                                         order_part_cent(2), order_part_cent(3), 0,   &
                                         num_cgto_bra, num_cgto_ket, num_part_opt,    &
                                         contr_ints)
          end if
        ! total geometric derivatives of Dirac delta function with Cartesian
        ! multipole moment integrals
        else
          call error_stop("contr_cgto_delta",                         &
                          "not implemented! please write to authors", &
                          order_mom)
          !! higher order geometric derivatives vanish for lower order Cartesian multipole moments
          !if (zero_ints .or. order_part_cent(4)>order_mom) then
          !  contr_ints = 0.0_REALK
          !else
          !end if
        end if
      end if
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "contr_cgto_delta", STDOUT)
#endif
    return
  end subroutine contr_cgto_delta

  !> \brief recurrence relations of Dirac delta function (with Cartesian multipole
  !>        moment) integrals using contracted Cartesian Gaussians
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
  !> \param delta_origin contains the coordinates of delta function origin
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scale constant for Dirac delta function
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_geo_bra is the order of geometric derivatives with respect to bra center
  !> \param order_geo_ket is the order of geometric derivatives with respect to ket center
  !> \param order_geo_pot is the order of geometric derivatives on delta function origin
  !> \param order_geo_mom is the order of geometric derivatives on dipole origin
  !> \param num_cgto_bra is the number of Cartesian GTOs on bra center,
  !>        equals to \f$(\var(angular_bra)+1)(\var(angular_bra)+2)/2\f$
  !> \param num_cgto_ket is the number of Cartesian GTOs on ket center,
  !>        equals to \f$(\var(angular_ket)+1)(\var(angular_ket)+2)/2\f$
  !> \param num_opt is the number of operators including derivatives
  !> \return contr_ints contains the contracted Cartesian multipole moment integrals
  subroutine contr_cgto_delta_recurr(coord_bra, angular_bra, num_prim_bra,        &
                                     exponent_bra, num_contr_bra, contr_coef_bra, &
                                     coord_ket, angular_ket, num_prim_ket,        &
                                     exponent_ket, num_contr_ket, contr_coef_ket, &
                                     order_elec, delta_origin, dipole_origin,     &
                                     scal_const, order_mom, order_geo_bra,        &
                                     order_geo_ket, order_geo_pot, order_geo_mom, &
                                     num_cgto_bra, num_cgto_ket, num_opt, contr_ints)
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
    real(REALK), intent(in) :: delta_origin(3)
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_geo_bra
    integer, intent(in) :: order_geo_ket
    integer, intent(in) :: order_geo_pot
    integer, intent(in) :: order_geo_mom
    integer, intent(in) :: num_cgto_bra
    integer, intent(in) :: num_cgto_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(out) :: contr_ints(num_cgto_bra,num_contr_bra, &
                                           num_cgto_ket,num_contr_ket,num_opt)
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
!f2py intent(in) :: delta_origin
!f2py intent(in) :: dipole_origin
!f2py intent(in) :: scal_const
!f2py intent(in) :: order_mom
!f2py intent(in) :: order_geo_bra
!f2py intent(in) :: order_geo_ket
!f2py intent(in) :: order_geo_pot
!f2py intent(in) :: order_geo_mom
!f2py intent(in) :: num_cgto_bra
!f2py intent(in) :: num_cgto_ket
!f2py intent(in) :: num_opt
!f2py intent(out) :: contr_ints
!f2py depend(num_cgto_bra) :: contr_ints
!f2py depend(num_contr_bra) :: contr_ints
!f2py depend(num_cgto_ket) :: contr_ints
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
    integer num_geo_bra         !number of geometric derivatives on bra center
    integer num_geo_ket         !number of geometric derivatives on ket center
    integer num_geo_pot         !number of geometric derivatives on delta function origin
    integer num_opt_geo         !number of operators and geometric derivatives
    real(REALK), allocatable :: hgto_pints(:,:,:,:)        !primitive Hermite integrals
    real(REALK), allocatable :: cgto_ket_pints(:,:,:,:,:)  !primitive Cartesian (ket) integrals
    real(REALK), allocatable :: cgto_pints(:,:,:,:,:)      !primitive Cartesian integrals
    real(REALK), allocatable :: low_mom_cints(:,:,:,:,:)   !contracted integrals of lower order
                                                           !Cartesian multipole moment
    integer iprim, jprim        !incremental recorders over primitives
    integer ierr                !error information
#if defined(XTIME)
    real(REALK) curr_time       !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! computes the minimum and maximum orders of Hermite Gaussians on bra and ket centers
    orders_hgto_bra(1) = order_geo_bra+mod(angular_bra,2)
    orders_hgto_bra(2) = order_geo_bra+angular_bra
    orders_hgto_ket(1) = order_geo_ket+mod(angular_ket,2)
    orders_hgto_ket(2) = order_geo_ket+angular_ket
    ! computes the order of Cartesian multipole moment after considering its geometric derivative
    order_low_mom = order_mom-order_geo_mom
    ! higher order geometric derivatives vanish for lower order Cartesian multipole moments
    if (order_low_mom<0) then
      contr_ints = 0.0_REALK
      return
    end if
    ! computes the number of primitive Hermite Gaussians used in recurrence relations
    if (orders_hgto_bra(1)==0) then
      dim_hgto_bra = (orders_hgto_bra(2)+1)*(orders_hgto_bra(2)+2) &
                   * (orders_hgto_bra(2)+3)/6
    else
      dim_hgto_bra = ((orders_hgto_bra(2)+1)*(orders_hgto_bra(2)+2) &
                      *(orders_hgto_bra(2)+3)                       &
                   -  orders_hgto_bra(1)*(orders_hgto_bra(1)+1)     &
                      *(orders_hgto_bra(1)+2))/6
    end if
    if (orders_hgto_ket(1)==0) then
      dim_hgto_ket = (orders_hgto_ket(2)+1)*(orders_hgto_ket(2)+2) &
                   * (orders_hgto_ket(2)+3)/6
    else
      dim_hgto_ket = ((orders_hgto_ket(2)+1)*(orders_hgto_ket(2)+2) &
                      *(orders_hgto_ket(2)+3)                       &
                   -  orders_hgto_ket(1)*(orders_hgto_ket(1)+1)     &
                      *(orders_hgto_ket(1)+2))/6
    end if
    ! computes the number of geometric derivatives
    num_geo_bra = (order_geo_bra+1)*(order_geo_bra+2)/2
    num_geo_ket = (order_geo_ket+1)*(order_geo_ket+2)/2
    num_geo_pot = (order_geo_pot+1)*(order_geo_pot+2)/2
    ! computes the number of operators, which is more efficient than using \var(TRI_SIZE)
    num_elec = (order_elec+1)*(order_elec+2)/2
    num_low_mom = (order_low_mom+1)*(order_low_mom+2)/2
    num_elec_mom = num_elec*num_low_mom
    num_opt_geo = num_elec_mom*num_geo_bra*num_geo_ket*num_geo_pot
    ! allocates the memory for primitive integrals
    allocate(hgto_pints(dim_hgto_bra,dim_hgto_ket,num_elec_mom,num_geo_pot), stat=ierr)
    if (ierr/=0)                                                                  &
      call error_stop("contr_cgto_delta_recurr", "failed to allocate hgto_pints", &
                      dim_hgto_bra*dim_hgto_ket*num_elec_mom*num_geo_pot)
    allocate(cgto_ket_pints(dim_hgto_bra,num_cgto_ket, &
                            num_elec_mom,num_geo_ket,num_geo_pot), stat=ierr)
    if (ierr/=0)                                                                      &
      call error_stop("contr_cgto_delta_recurr", "failed to allocate cgto_ket_pints", &
                      dim_hgto_bra*num_cgto_ket*num_elec_mom*num_geo_ket*num_geo_pot)
    allocate(cgto_pints(num_cgto_bra,num_cgto_ket,num_opt_geo, &
                        num_prim_bra,num_prim_ket), stat=ierr)
    if (ierr/=0)                                                                  &
      call error_stop("contr_cgto_delta_recurr", "failed to allocate cgto_pints", &
                      num_cgto_bra*num_cgto_ket*num_opt_geo*num_prim_bra*num_prim_ket)
#if defined(DEBUG)
    write(STDOUT,100) "size of primitive Hermite integrals:", size(hgto_pints)
    write(STDOUT,100) "size of primitive Cartesian (ket) integrals:", size(cgto_ket_pints)
    write(STDOUT,100) "size of primitive Cartesian integrals:", size(cgto_pints)
#endif
    ! computes the primitive integrals
    do jprim = 1, num_prim_ket
      do iprim = 1, num_prim_bra
        ! computes the primitive Hermite integrals
        call prim_hgto_delta(orders_hgto_bra, coord_bra, exponent_bra(iprim), &
                             orders_hgto_ket, coord_ket, exponent_ket(jprim), &
                             order_elec, delta_origin, dipole_origin,         &
                             scal_const, order_low_mom, order_geo_pot,        &
                             dim_hgto_bra, dim_hgto_ket, num_elec,            &
                             num_low_mom, num_geo_pot, hgto_pints)
        ! transforms Hermite Gaussians on ket center to Cartesian ones first
        call hgto_to_cgto(angular_ket, order_geo_ket, exponent_ket(jprim),       &
                          dim_hgto_bra, dim_hgto_ket, num_elec_mom, num_geo_pot, &
                          hgto_pints, num_cgto_ket, num_geo_ket, cgto_ket_pints)
        ! then transforms Hermite Gaussians on bra center
        call hgto_to_cgto(angular_bra, order_geo_bra, exponent_bra(iprim), &
                          1, dim_hgto_bra, num_cgto_ket*num_elec_mom,      &
                          num_geo_ket*num_geo_pot, cgto_ket_pints,         &
                          num_cgto_bra, num_geo_bra, cgto_pints(:,:,:,iprim,jprim))
      end do
    end do
    ! cleans
    deallocate(hgto_pints)
    deallocate(cgto_ket_pints)
    if (order_geo_mom==0) then
      ! constructs contracted integrals without geometric derivative on dipole origin,
      ! in the order of contr_ints(num_cgto_bra,num_contr_bra,num_cgto_ket,num_contr_ket, &
      !                            num_elec,num_mom,num_geo_bra,num_geo_ket,num_geo_pot)
      call const_contr_ints(num_contr_bra, num_prim_bra, contr_coef_bra, &
                            num_contr_ket, num_prim_ket, contr_coef_ket, &
                            num_cgto_bra, num_cgto_ket, num_opt_geo,     &
                            cgto_pints, contr_ints)
      ! cleans
      deallocate(cgto_pints)
    else
      ! constructs contracted integrals of lower order Cartesian multipole memonts,
      ! in the order of low_mom_cints(num_cgto_bra,num_contr_bra,num_cgto_ket,num_contr_ket, &
      !                               num_elec,num_low_mom,num_geo_bra,num_geo_ket,num_geo_pot)
      allocate(low_mom_cints(num_cgto_bra,num_contr_bra, &
                             num_cgto_ket,num_contr_ket,num_opt_geo), stat=ierr)
      if (ierr/=0)                                                                     &
        call error_stop("contr_cgto_delta_recurr", "failed to allocate low_mom_cints", &
                        num_cgto_bra*num_contr_bra*num_cgto_ket*num_contr_ket*num_opt_geo)
#if defined(DEBUG)
      write(STDOUT,100)                                                             &
        "size of contracted integrals of lower order Cartesian multipole memonts:", &
        size(low_mom_cints)
#endif
      call const_contr_ints(num_contr_bra, num_prim_bra, contr_coef_bra, &
                            num_contr_ket, num_prim_ket, contr_coef_ket, &
                            num_cgto_bra, num_cgto_ket, num_opt_geo,     &
                            cgto_pints, low_mom_cints)
      ! cleans
      deallocate(cgto_pints)
      ! constructs the geometric derivatives of Cartesian multipole moments
      ! in the order of contr_ints(num_cgto_bra,num_contr_bra,num_cgto_ket,num_contr_ket, &
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
    call xtimer_view(curr_time, "contr_cgto_delta_recurr", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("contr_cgto_delta_recurr>> ",A,I8,1X,A)
#endif
  end subroutine contr_cgto_delta_recurr
