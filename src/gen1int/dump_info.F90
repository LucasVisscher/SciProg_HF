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
!!  This file dumps the information of contracted GTOs, derivatives and operators.
!!
!!  2011-03-24, Bin Gao:
!!  * first version

  !> \brief dumps the information of contracted GTOs and partial derivatives
  !> \author Bin Gao
  !> \date 2011-03-24
  !> \param type_gto specifies the type of GTOs
  !> \param idx_gto is the atomic index of GTO center
  !> \param coord_gto contains the coordinates of GTO center
  !> \param angular_num is the angular number (s=0, p=1, d=2, ...)
  !> \param num_prim is the number of primitive Gaussians
  !> \param exponents contains the exponents of primitive Gaussians
  !> \param num_contr is the number of contractions
  !> \param contr_coef contains the contraction coefficients
  !> \param order_part_mag is the order of partial magnetic derivatives
  !> \param order_part_ram is the order of partial derivatives with respect to
  !>        total rotational angular momentum
  !> \param order_part_geo is the order of partial geometric derivatives
  !> \param io_viewer is the logical unit number of the viewer
  subroutine dump_gto_pd(type_gto, idx_gto, coord_gto, angular_num,      &
                         num_prim, exponents, num_contr, contr_coef,     &
                         order_part_mag, order_part_ram, order_part_geo, &
                         io_viewer)
    use xkind
    implicit none
    character*(*), intent(in) :: type_gto
    integer, intent(in) :: idx_gto
    real(REALK), intent(in) :: coord_gto(3)
    integer, intent(in) :: angular_num
    integer, intent(in) :: num_prim
    real(REALK), intent(in) :: exponents(num_prim)
    integer, intent(in) :: num_contr
    real(REALK), intent(in) :: contr_coef(num_contr,num_prim)
    integer, intent(in) :: order_part_mag
    integer, intent(in) :: order_part_ram
    integer, intent(in) :: order_part_geo
    integer, intent(in) :: io_viewer
!f2py intent(in) :: type_gto
!f2py intent(in) :: idx_gto
!f2py intent(in) :: coord_gto
!f2py intent(in) :: angular_num
!f2py intent(hide) :: num_prim
!f2py intent(in) :: exponents
!f2py intent(hide) :: num_contr
!f2py intent(in) :: contr_coef
!f2py depend(num_prim) :: contr_coef
!f2py intent(in) :: order_part_mag
!f2py intent(in) :: order_part_ram
!f2py intent(in) :: order_part_geo
!f2py intent(in) :: io_viewer
    integer iprim  !incremental recorders over primitives
    ! information of GTOs
    write(io_viewer,100) "contracted "//type_gto
    write(io_viewer,100) "atom center:", idx_gto
    write(io_viewer,110) "coordinates:", coord_gto
    write(io_viewer,100) "angular number:", angular_num
    write(io_viewer,100) "number of primitives:", num_prim
    write(io_viewer,100) "number of contractions:", num_contr
    write(io_viewer,100) "    exponents      coefficients"
    do iprim = 1, num_prim
      write(io_viewer,120) exponents(iprim), contr_coef(1:min(4,num_contr),iprim)
      write(io_viewer,130) contr_coef(5:num_contr,iprim)
    end do
    ! information of partial derivatives
    if (order_part_mag>0) &
      write(io_viewer,100) "order of partial magnetic derivatives:", order_part_mag
    if (order_part_ram>0)                                                       &
      write(io_viewer,100)                                                      &
      "order of partial derivatives w.r.t. total rotational angular momentum ", &
      order_part_ram
    if (order_part_geo>0) &
      write(io_viewer,100) "order of partial geometric derivatives", order_part_geo
    return
100 format("dump_gto_pd>> ",A,I8)
110 format("dump_gto_pd>> ",A,3F16.8)
120 format("dump_gto_pd>> ",5Es16.8)
130 format("dump_gto_pd>> ",16X,4Es16.8)
  end subroutine dump_gto_pd

  !> \brief dumps the information of total derivatives
  !> \author Bin Gao
  !> \date 2011-03-24
  !> \param order_mag is the order of total magnetic derivatives
  !> \param order_ram is the order of total derivatives with respect to
  !>        total rotational angular momentum
  !> \param num_cents is the number of differentiated centers
  !> \param idx_cent contains the indices of differentiated centers
  !> \param order_cent contains the order of total geometric derivatives of
  !>        the corresponding differentiated centers
  !> \param io_viewer is the logical unit number of the viewer
  subroutine dump_total_derv(order_mag, order_ram, num_cents, &
                             idx_cent, order_cent, io_viewer)
    use xkind
    implicit none
    integer, intent(in) :: order_mag
    integer, intent(in) :: order_ram
    integer, intent(in) :: num_cents
    integer, intent(in) :: idx_cent(num_cents)
    integer, intent(in) :: order_cent(num_cents)
    integer, intent(in) :: io_viewer
!f2py intent(in) :: order_mag
!f2py intent(in) :: order_ram
!f2py intent(hide) :: num_cents
!f2py intent(in) :: idx_cent
!f2py intent(in) :: order_cent
!f2py depend(num_cents) :: order_cent
!f2py intent(in) :: io_viewer
    integer icent  !incremental recorders over centers
    if (order_mag>0) &
      write(io_viewer,100) "order of total magnetic derivatives:", order_mag
    if (order_ram>0) write(io_viewer,100)                                     &
      "order of total derivatives w.r.t. total rotational angular momentum:", &
      order_ram
    if (num_cents>0) then
      write(io_viewer,100) "number of differentiated centers:", num_cents
      write(io_viewer,110) (idx_cent(icent), order_cent(icent), icent=1,num_cents)
    end if
    return
100 format("dump_total_derv>> ",A,I8)
110 format("dump_total_derv>> center(order):",5(I6,"(",I3,")"))
  end subroutine dump_total_derv

  !> \brief dumps the information of Cartesian multipole moments
  !> \author Bin Gao
  !> \date 2011-03-24
  !> \param order_elec is the order of electronic derivatives
  !> \param idx_diporg is the atomic center of dipole origin (<1 for non-atomic center)
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scale constant for Cartesian multipole moments
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_geo_mom is the order of geometric derivatives on dipole origin
  !> \param io_viewer is the logical unit number of the viewer
  subroutine dump_carmom(order_elec, idx_diporg, dipole_origin, scal_const, &
                         order_mom, order_geo_mom, io_viewer)
    use xkind
    implicit none
    integer, intent(in) :: order_elec
    integer, intent(in) :: idx_diporg
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_geo_mom
    integer, intent(in) :: io_viewer
!f2py intent(in) :: order_elec
!f2py intent(in) :: idx_diporg
!f2py intent(in) :: dipole_origin
!f2py intent(in) :: scal_const
!f2py intent(in) :: order_mom
!f2py intent(in) :: order_geo_mom
!f2py intent(in) :: io_viewer
    write(io_viewer,100) "atom center of dipole origin:", idx_diporg
    write(io_viewer,110) "coordinates of dipole origin:", dipole_origin
    write(io_viewer,110) "scale constant:", scal_const
    write(io_viewer,100) "order of multipole moments:", order_mom
    write(io_viewer,100) "order of partial geometric derivatives on dipole origin:", &
                         order_geo_mom
    write(io_viewer,100) "order of electronic derivatives:", order_elec
    return
100 format("dump_carmom>> ",A,I8)
110 format("dump_carmom>> ",A,3F16.8)
  end subroutine dump_carmom

  !> \brief dumps the information of Dirac delta function integrals
  !> \author Bin Gao
  !> \date 2011-03-24
  !> \param order_elec is the order of electronic derivatives
  !> \param idx_delta is the atomic center of delta function origin (<1 for non-atomic center)
  !> \param delta_origin contains the coordinates of delta function origin
  !> \param idx_diporg is the atomic center of dipole origin (<1 for non-atomic center)
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scale constant for Dirac delta function
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_geo_pot is the order of geometric derivatives on delta function origin
  !> \param order_geo_mom is the order of geometric derivatives on dipole origin
  !> \param io_viewer is the IO unit to dump the information
  subroutine dump_delta(order_elec, idx_delta, delta_origin, idx_diporg,     &
                        dipole_origin, scal_const, order_mom, order_geo_pot, &
                        order_geo_mom, io_viewer)
    use xkind
    implicit none
    integer, intent(in) :: order_elec
    integer, intent(in) :: idx_delta
    real(REALK), intent(in) :: delta_origin(3)
    integer, intent(in) :: idx_diporg
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_geo_pot
    integer, intent(in) :: order_geo_mom
    integer, intent(in) :: io_viewer
!f2py intent(in) :: order_elec
!f2py intent(in) :: idx_delta
!f2py intent(in) :: delta_origin
!f2py intent(in) :: idx_diporg
!f2py intent(in) :: dipole_origin
!f2py intent(in) :: scal_const
!f2py intent(in) :: order_mom
!f2py intent(in) :: order_geo_pot
!f2py intent(in) :: order_geo_mom
!f2py intent(in) :: io_viewer
    write(io_viewer,100) "atom center of Dirac delta function:", idx_delta
    write(io_viewer,110) "coordinates of center of Dirac delta function:", &
                         delta_origin
    write(io_viewer,110) "scale constant:", scal_const
    write(io_viewer,100)                                                 &
      "order of partial geometric derivatives on Dirac delta function:", &
      order_geo_pot
    write(io_viewer,100) "atom center of dipole origin:", idx_diporg
    write(io_viewer,110) "coordinates of dipole origin:", dipole_origin
    write(io_viewer,100) "order of multipole moments:", order_mom
    write(io_viewer,100) "order of partial geometric derivatives on dipole origin:", &
                         order_geo_mom
    write(io_viewer,100) "order of electronic derivatives:", order_elec
    return
100 format("dump_delta>> ",A,I8)
110 format("dump_delta>> ",A,3F16.8)
  end subroutine dump_delta

  !> \brief dumps the information of nuclear attraction potential integrals
  !> \author Bin Gao
  !> \date 2011-03-24
  !> \param order_elec is the order of electronic derivatives
  !> \param idx_nucorg is the atomic center of nuclear potential origin (<1 for non-atomic center)
  !> \param nucpot_origin contains the coordinates of nuclear potential origin
  !> \param idx_diporg is the atomic center of dipole origin (<1 for non-atomic center)
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scale constant for nuclear attraction potential
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_geo_pot is the order of geometric derivatives on nuclear potential origin
  !> \param order_geo_mom is the order of geometric derivatives on dipole origin
  !> \param io_viewer is the IO unit to dump the information
  subroutine dump_nucpot(order_elec, idx_nucorg, nucpot_origin, idx_diporg,   &
                         dipole_origin, scal_const, order_mom, order_geo_pot, &
                         order_geo_mom, io_viewer)
    use xkind
    implicit none
    integer, intent(in) :: order_elec
    integer, intent(in) :: idx_nucorg
    real(REALK), intent(in) :: nucpot_origin(3)
    integer, intent(in) :: idx_diporg
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_geo_pot
    integer, intent(in) :: order_geo_mom
    integer, intent(in) :: io_viewer
!f2py intent(in) :: order_elec
!f2py intent(in) :: idx_nucorg
!f2py intent(in) :: nucpot_origin
!f2py intent(in) :: idx_diporg
!f2py intent(in) :: dipole_origin
!f2py intent(in) :: scal_const
!f2py intent(in) :: order_mom
!f2py intent(in) :: order_geo_pot
!f2py intent(in) :: order_geo_mom
!f2py intent(in) :: io_viewer
    write(io_viewer,100) "atom center of nuclear attraction potential:", idx_nucorg
    write(io_viewer,110) "coordinates of center of nuclear attraction potential:", &
                         nucpot_origin
    write(io_viewer,110) "scale constant:", scal_const
    write(io_viewer,100)                                                         &
      "order of partial geometric derivatives on nuclear attraction potential:", &
      order_geo_pot
    write(io_viewer,100) "atom center of dipole origin:", idx_diporg
    write(io_viewer,110) "coordinates of dipole origin:", dipole_origin
    write(io_viewer,100) "order of multipole moments:", order_mom
    write(io_viewer,100) "order of partial geometric derivatives on dipole origin:", &
                         order_geo_mom
    write(io_viewer,100) "order of electronic derivatives:", order_elec
    return
100 format("dump_nucpot>> ",A,I8)
110 format("dump_nucpot>> ",A,3F16.8)
  end subroutine dump_nucpot

  !> \brief dumps the information of inverse square distance potential integrals
  !> \author Bin Gao
  !> \date 2011-03-24
  !> \param order_elec is the order of electronic derivatives
  !> \param idx_isdorg is the atomic center of inverse square distance potential
  !>        origin (<1 for non-atomic center)
  !> \param isdpot_origin contains the coordinates of inverse square distance potential origin
  !> \param idx_diporg is the atomic center of dipole origin (<1 for non-atomic center)
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scale constant for inverse square distance potential
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_geo_pot is the order of geometric derivatives on inverse square distance potential origin
  !> \param order_geo_mom is the order of geometric derivatives on dipole origin
  !> \param io_viewer is the IO unit to dump the information
  subroutine dump_isdpot(order_elec, idx_isdorg, isdpot_origin, idx_diporg,   &
                         dipole_origin, scal_const, order_mom, order_geo_pot, &
                         order_geo_mom, io_viewer)
    use xkind
    implicit none
    integer, intent(in) :: order_elec
    integer, intent(in) :: idx_isdorg
    real(REALK), intent(in) :: isdpot_origin(3)
    integer, intent(in) :: idx_diporg
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_geo_pot
    integer, intent(in) :: order_geo_mom
    integer, intent(in) :: io_viewer
!f2py intent(in) :: order_elec
!f2py intent(in) :: idx_isdorg
!f2py intent(in) :: isdpot_origin
!f2py intent(in) :: idx_diporg
!f2py intent(in) :: dipole_origin
!f2py intent(in) :: scal_const
!f2py intent(in) :: order_mom
!f2py intent(in) :: order_geo_pot
!f2py intent(in) :: order_geo_mom
!f2py intent(in) :: io_viewer
    write(io_viewer,100) "atom center of inverse square distance potential:", idx_isdorg
    write(io_viewer,110) "coordinates of center of inverse square distance potential:", &
                         isdpot_origin
    write(io_viewer,110) "scale constant:", scal_const
    write(io_viewer,100)                                                              &
      "order of partial geometric derivatives on inverse square distance potential:", &
      order_geo_pot
    write(io_viewer,100) "atom center of dipole origin:", idx_diporg
    write(io_viewer,110) "coordinates of dipole origin:", dipole_origin
    write(io_viewer,100) "order of multipole moments:", order_mom
    write(io_viewer,100) "order of partial geometric derivatives on dipole origin:", &
                         order_geo_mom
    write(io_viewer,100) "order of electronic derivatives:", order_elec
    return
100 format("dump_isdpot>> ",A,I8)
110 format("dump_isdpot>> ",A,3F16.8)
  end subroutine dump_isdpot

  !> \brief dumps the information of Gaussian charge potential integrals
  !> \author Bin Gao
  !> \date 2011-03-24
  !> \param order_elec is the order of electronic derivatives
  !> \param idx_gauorg is the atomic center of Gaussian charge potential origin (<1 for non-atomic center)
  !> \param gaupot_origin contains the coordinates of Gaussian charge potential origin
  !> \param gaupot_expt is the exponent used in the Gaussian broadening function of the charge
  !> \param idx_diporg is the atomic center of dipole origin (<1 for non-atomic center)
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scale constant for Gaussian charge potential
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_geo_pot is the order of geometric derivatives on gaupot function origin
  !> \param order_geo_mom is the order of geometric derivatives on dipole origin
  !> \param io_viewer is the IO unit to dump the information
  subroutine dump_gaupot(order_elec, idx_gauorg, gaupot_origin, gaupot_expt, &
                         idx_diporg, dipole_origin, scal_const, order_mom,   &
                         order_geo_pot, order_geo_mom, io_viewer)
    use xkind
    implicit none
    integer, intent(in) :: order_elec
    integer, intent(in) :: idx_gauorg
    real(REALK), intent(in) :: gaupot_origin(3)
    real(REALK), intent(in) :: gaupot_expt
    integer, intent(in) :: idx_diporg
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_geo_pot
    integer, intent(in) :: order_geo_mom
    integer, intent(in) :: io_viewer
!f2py intent(in) :: order_elec
!f2py intent(in) :: idx_gauorg
!f2py intent(in) :: gaupot_origin
!f2py intent(in) :: gaupot_expt
!f2py intent(in) :: idx_diporg
!f2py intent(in) :: dipole_origin
!f2py intent(in) :: scal_const
!f2py intent(in) :: order_mom
!f2py intent(in) :: order_geo_pot
!f2py intent(in) :: order_geo_mom
!f2py intent(in) :: io_viewer
    write(io_viewer,100) "atom center of Gaussian charge potential:", idx_gauorg
    write(io_viewer,110) "coordinates of center of Gaussian charge potential:", &
                         gaupot_origin
    write(io_viewer,110) "exponent in Gaussian broadening function of the charge", &
                         gaupot_expt
    write(io_viewer,110) "scale constant:", scal_const
    write(io_viewer,100)                                                      &
      "order of partial geometric derivatives on Gaussian charge potential:", &
      order_geo_pot
    write(io_viewer,100) "atom center of dipole origin:", idx_diporg
    write(io_viewer,110) "coordinates of dipole origin:", dipole_origin
    write(io_viewer,100) "order of multipole moments:", order_mom
    write(io_viewer,100) "order of partial geometric derivatives on dipole origin:", &
                         order_geo_mom
    write(io_viewer,100) "order of electronic derivatives:", order_elec
    return
100 format("dump_gaupot>> ",A,I8)
110 format("dump_gaupot>> ",A,3F16.8)
  end subroutine dump_gaupot

  !> \brief dumps the information of diamagnetic spin-orbit coupling integrals
  !> \author Bin Gao
  !> \date 2011-03-24
  !> \param order_elec is the order of electronic derivatives
  !> \param idx_dso contains the atomic centers of diamagnetic spin-orbit coupling
  !>        origins (<1 for non-atomic center)
  !> \param coord_dso contains the coordinates of diamagnetic spin-orbit coupling origins
  !> \param scal_const is the scale constant for diamagnetic spin-orbit coupling
  !> \param order_geo_dso contains the orders of geometric derivatives on diamagnetic
  !>        spin-orbit coupling origins
  !> \param io_viewer is the IO unit to dump the information
  subroutine dump_dso(order_elec, idx_dso, coord_dso, scal_const, &
                      order_geo_dso, io_viewer)
    use xkind
    implicit none
    integer, intent(in) :: order_elec
    integer, intent(in) :: idx_dso(2)
    real(REALK), intent(in) :: coord_dso(3,2)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: order_geo_dso(2)
    integer, intent(in) :: io_viewer
!f2py intent(in) :: order_elec
!f2py intent(in) :: idx_dso
!f2py intent(in) :: coord_dso
!f2py intent(in) :: scal_const
!f2py intent(in) :: order_geo_dso
!f2py intent(in) :: io_viewer
    write(io_viewer,100) "first atom center:", idx_dso(1)
    write(io_viewer,110) "coordinates of the first atom center:", coord_dso(:,1)
    write(io_viewer,100)                                                  &
      "order of partial geometric derivatives on the first atom center:", &
      order_geo_dso(1)
    write(io_viewer,100) "second atom center:", idx_dso(2)
    write(io_viewer,110) "coordinates of the second atom center:", coord_dso(:,2)
    write(io_viewer,100)                                                   &
      "order of partial geometric derivatives on the second atom center:", &
      order_geo_dso(2)
    write(io_viewer,110) "scale constant:", scal_const
    write(io_viewer,100) "order of electronic derivatives:", order_elec
    return
100 format("dump_dso>> ",A,I8)
110 format("dump_dso>> ",A,3F16.8)
  end subroutine dump_dso

  !> \brief dumps the information of effective core potential (ECP)
  !> \author Bin Gao
  !> \date 2011-03-24
  !> \param ecp_local indicates if it is the local or non-local part
  !> \param idx_ecp is the atomic index of ECP center
  !> \param coord_ecp is the coordinates of ECP center
  !> \param core_angl is the number of core electrons (local part) or angular number (non-local part)
  !> \param num_gauss is the number of Gaussians combining the radical function
  !> \param gau_rpower contains the powers of r in Gaussians
  !> \param gau_exponent contains the exponents of Gaussians
  !> \param gau_coeff contains the combination coefficients of Gaussians
  !> \param io_viewer is the logical unit number of the viewer
  subroutine dump_ecp(ecp_local, idx_ecp, coord_ecp, core_angl, num_gauss, &
                      gau_rpower, gau_exponent, gau_coeff, io_viewer)
    use xkind
    implicit none
    logical, intent(in) :: ecp_local
    integer, intent(in) :: idx_ecp
    real(REALK), intent(in) :: coord_ecp(3)
    integer, intent(in) :: core_angl
    integer, intent(in) :: num_gauss
    integer, intent(in) :: gau_rpower(num_gauss)
    real(REALK), intent(in) :: gau_exponent(num_gauss)
    real(REALK), intent(in) :: gau_coeff(num_gauss)
    integer, intent(in) :: io_viewer
!f2py intent(in) :: ecp_local
!f2py intent(in) :: idx_ecp
!f2py intent(in) :: coord_ecp
!f2py intent(in) :: core_angl
!f2py intent(hide) :: num_gauss
!f2py intent(in) :: gau_rpower
!f2py intent(in) :: gau_exponent
!f2py depend(num_gauss) :: gau_exponent
!f2py intent(in) :: gau_coeff
!f2py depend(num_gauss) :: gau_coeff
!f2py intent(in) :: io_viewer
    integer igau  !incremental recorder over Gaussians
    if (ecp_local) then
      write(io_viewer,100) "operator>> effective core potential (local part)"
      write(io_viewer,100) "number of core electrons:", core_angl
    else
      write(io_viewer,100) "operator>> effective core potential (non-local part)"
      write(io_viewer,100) "angular number:", core_angl
    end if
    write(io_viewer,100) "atom center of ECP:", idx_ecp
    write(io_viewer,110) "coordinates of ECP center:", coord_ecp
    write(io_viewer,100) "    powers      exponents      coefficients"
    do igau = 1, num_gauss
      write(io_viewer,120) -gau_rpower(igau), gau_exponent(igau), gau_coeff(igau)
    end do
    return
100 format("dump_ecp>> ",A,I8)
110 format("dump_ecp>> ",A,3F16.8)
120 format("dump_ecp>> ",I4,2F16.8)
  end subroutine dump_ecp

  !> \brief dumps the information of model core potential (Version 1) -- potential
  !> \author Bin Gao
  !> \date 2011-03-24
  !> \param idx_mcp is the atomic index of MCP1 center
  !> \param coord_mcp is the coordinates of MCP1 center
  !> \param charge_mcp is the effective core charge
  !> \param num_gauss is the number of Gaussians in the model potential
  !> \param gau_rpower contains the powers of r in Gaussians
  !> \param gau_exponent contains the exponents of Gaussians
  !> \param gau_coeff contains the combination coefficients of Gaussians
  !> \param io_viewer is the logical unit number of the viewer
  subroutine dump_mcp1_pot(idx_mcp, coord_mcp, charge_mcp, num_gauss, &
                           gau_rpower, gau_exponent, gau_coeff, io_viewer)
    use xkind
    implicit none
    integer, intent(in) :: idx_mcp
    real(REALK), intent(in) :: coord_mcp(3)
    real(REALK), intent(in) :: charge_mcp
    integer, intent(in) :: num_gauss
    integer, intent(in) :: gau_rpower(num_gauss)
    real(REALK), intent(in) :: gau_exponent(num_gauss)
    real(REALK), intent(in) :: gau_coeff(num_gauss)
    integer, intent(in) :: io_viewer
!f2py intent(in) :: idx_mcp
!f2py intent(in) :: coord_mcp
!f2py intent(in) :: charge_mcp
!f2py intent(hide) :: num_gauss
!f2py intent(in) :: gau_rpower
!f2py intent(in) :: gau_exponent
!f2py depend(num_gauss) :: gau_exponent
!f2py intent(in) :: gau_coeff
!f2py depend(num_gauss) :: gau_coeff
!f2py intent(in) :: io_viewer
    integer igau  !incremental recorder over Gaussians
    write(io_viewer,100) "operator>> model core potential (Version 1) -- potential"
    write(io_viewer,100) "atom center of MCP1:", idx_mcp
    write(io_viewer,110) "coordinates of MCP1 center:", coord_mcp
    write(io_viewer,110) "effective core charge:", charge_mcp
    write(io_viewer,100) "    powers      exponents      coefficients"
    do igau = 1, num_gauss
      write(io_viewer,120) gau_rpower(igau), gau_exponent(igau), gau_coeff(igau)
    end do
    return
100 format("dump_mcp1_pot>> ",A,I8)
110 format("dump_mcp1_pot>> ",A,3F16.8)
120 format("dump_mcp1_pot>> ",I4,2F16.8)
  end subroutine dump_mcp1_pot

  !> \brief dumps the information of model core potential (Version 1) -- core orbitals
  !> \author Bin Gao
  !> \date 2011-03-24
  !> \param idx_mcp is the atomic index of MCP1 center
  !> \param coord_mcp is the coordinates of MCP1 center
  !> \param core_angnum is the angular number of core orbitals
  !> \param param_eshfit is the parameter of energy shift operator
  !> \param core_nexp is the number of exponents of core orbtials
  !> \param core_exponent contains the orbital exponents
  !> \param core_nshell is the number of shells of core orbitals
  !> \param core_coeff contains the expansion coefficients
  !> \param io_viewer is the logical unit number of the viewer
  subroutine dump_mcp1_core(idx_mcp, coord_mcp, core_angnum, param_eshfit, &
                            core_nexp, core_exponent, core_nshell,         &
                            core_coeff, io_viewer)
    use xkind
    implicit none
    integer, intent(in) :: idx_mcp
    real(REALK), intent(in) :: coord_mcp(3)
    integer, intent(in) :: core_angnum
    real(REALK), intent(in) :: param_eshfit
    integer, intent(in) :: core_nexp
    real(REALK), intent(in) :: core_exponent(core_nexp)
    integer, intent(in) :: core_nshell
    real(REALK), intent(in) :: core_coeff(core_nshell,core_nexp)
    integer, intent(in) :: io_viewer
!f2py intent(in) :: idx_mcp
!f2py intent(in) :: coord_mcp
!f2py intent(in) :: core_angnum
!f2py intent(in) :: param_eshfit
!f2py intent(hide) :: core_nexp
!f2py intent(in) :: core_exponent
!f2py intent(hide) :: core_nshell
!f2py intent(in) :: core_coeff
!f2py depend(core_nexp) :: core_coeff
!f2py intent(in) :: io_viewer
    integer iexp  !incremental recorder over orbital exponents
    write(io_viewer,100) "operator>> model core potential (Version 1) -- core orbitals"
    write(io_viewer,100) "atom center of MCP1:", idx_mcp
    write(io_viewer,110) "coordinates of MCP1 center:", coord_mcp
    write(io_viewer,100) "angular number of core orbitals:", core_angnum
    write(io_viewer,110) "parameter of energy shift operator:", param_eshfit
    write(io_viewer,100) "number of exponents of core orbtials:", core_nexp
    write(io_viewer,100) "number of shells of core orbitals:", core_nshell
    write(io_viewer,100) "    exponents      coefficients"
    do iexp = 1, core_nexp
      write(io_viewer,120) core_exponent(iexp), core_coeff(1:min(4,core_nshell),iexp)
      write(io_viewer,130) core_coeff(5:core_nshell,iexp)
    end do
    return
100 format("dump_mcp1_core>> ",A,I8)
110 format("dump_mcp1_core>> ",A,3F16.8)
120 format("dump_dump_mcp1_core>> ",5F14.7)
130 format("dump_dump_mcp1_core>> ",14X,4F14.7)
  end subroutine dump_mcp1_core
