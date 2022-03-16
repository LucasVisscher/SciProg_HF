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
!!  This file normalizes the contracted spherical Gaussians.
!!
!!  2009-10-26, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief normalizes the contracted spherical Gaussians
  !> \details the normalized contraction coefficient \f$\bar{C}(\alpha_{nl})\f$ is,
  !> \f{eqnarray}{
  !>   \bar{C}(\alpha_{nl})
  !>   &=&\left(\frac{2}{\pi}\right)^{1/4}\frac{(4\alpha_{nl})^{l/2+3/4}}{\sqrt{(2l+1)!!}} %
  !>     \frac{C(\alpha_{nl})}{\sqrt{\sum_{\alpha_{i},\alpha_{j}}C(\alpha_{i})C(\alpha_{j})%
  !>     \left(\frac{\sqrt{4\alpha_{i}\alpha_{j}}}{\alpha_{i}+\alpha_{j}}\right)^{l+3/2}}},
  !> \f}
  !> where \f$\left(\frac{2}{\pi}\right)^{1/4}\frac{(4\alpha_{nl})^{l/2+3/4}}{\sqrt{(2l+1)!!}}\f$
  !> comes from the normalization on primitive Gaussians, and \f$\frac{1}{\sqrt{%
  !> \sum_{\alpha_{i},\alpha_{j}}C(\alpha_{i})C(\alpha_{j})%
  !> \left(\frac{\sqrt{4\alpha_{i}\alpha_{j}}}{\alpha_{i}+\alpha_{j}}\right)^{l+3/2}}}\f$
  !> comes from the radial overlap between contracted GTOs (see, for example, Eq. (6S.9.3),
  !> "Molecular Electronic Structure Theory", pp. 254)
  !> \author Bin Gao
  !> \date 2009-10-26
  !> \param angular_num is the angular number
  !> \param num_prim is the number of primitives
  !> \param exponents contains the exponents of primitives
  !> \param num_contr is the number of contractions
  !> \return contr_coef contains the contraction coefficients
  !> \note the normalization over the angles will be taken care by subroutine \fn(cgto_to_sgto)
  subroutine norm_contr_sgto(angular_num, num_prim, exponents, num_contr, contr_coef)
    use xkind
    implicit none
    integer, intent(in) :: angular_num
    integer, intent(in) :: num_prim
    real(REALK), intent(in) :: exponents(num_prim)
    integer, intent(in) :: num_contr
    real(REALK), intent(inout) :: contr_coef(num_contr,num_prim)
!f2py intent(in) :: angular_num
!f2py intent(hide) :: num_prim
!f2py intent(in) :: exponents
!f2py intent(hide) :: num_contr
!f2py intent(inout) :: contr_coef
!f2py depend(num_prim) :: contr_coef
#include "private/pi.h"
    real(REALK), parameter :: NORM_THRSH = epsilon(1.0_REALK)
                                      !threshold of normalization constant
!radovan: not portable
!   real(REALK), parameter :: PREFACT_NRM = sqrt(sqrt(2.0_REALK/PI))
    real(REALK), parameter :: PREFACT_NRM = 0.893243841738002_REALK
                                      !constant \f$\left(\frac{2}{\pi}\right)^{1/4}\f$
    real(REALK) exp_power             !\f$l/2+3/4\f$
    real(REALK) overlap_power         !\f$l+3/2\f$
    real(REALK) ang_prefact           !\f$\left(\frac{2}{\pi}\right)^{1/4}\frac{1}{\sqrt{(2l+1)!!}}\f$
    real(REALK) rad_overlap           !radial overlap between two GTOs
    real(REALK) norm_const            !normalization constant
    integer icontr                    !incremental recorder over contractions
    integer iprim, jprim              !incremental recorders over primitive Gaussians
    integer iang                      !incremental recorder over angular numbers
    ! computes \f$l/2+3/4\f$ and \f$l+3/2\f$
    overlap_power = real(angular_num,REALK)+1.5_REALK
    exp_power = 0.5_REALK*overlap_power
    ! computes \f$\left(\frac{2}{\pi}\right)^{1/4}\frac{1}{\sqrt{(2l+1)!!}}\f$
    ang_prefact = 1.0_REALK
    do iang = 3, angular_num+angular_num+1, 2
      ang_prefact = real(iang,REALK)*ang_prefact
    end do
    ang_prefact = sqrt(ang_prefact)
    ang_prefact = PREFACT_NRM/ang_prefact
    ! loops over contractions
    do icontr = 1, num_contr
      norm_const = 0.0_REALK
      ! computes \f$\frac{1}{\sqrt{\sum_{\alpha_{i},\alpha_{j}}C(\alpha_{i})C(\alpha_{j})%
      ! \left(\frac{\sqrt{4\alpha_{i}\alpha_{j}}}{\alpha_{i}+\alpha_{j}}\right)^{l+3/2}}}\f$
      do iprim = 1, num_prim
        do jprim = 1, num_prim
          !FIXME: checks something, like exponents^(-angular_num-3/2)?
          ! computes the radial overlap between two GTOs
          rad_overlap = sqrt(exponents(iprim)*exponents(jprim)) &
                      / (exponents(iprim)+exponents(jprim))
          rad_overlap = rad_overlap+rad_overlap
          rad_overlap = rad_overlap**overlap_power
          ! multiplied by the contraction coefficients
          norm_const = norm_const+contr_coef(icontr,iprim) &
                     * contr_coef(icontr,jprim)*rad_overlap
        end do
      end do
      norm_const = sqrt(norm_const)
      ! checks the normalization constant
      if (norm_const<NORM_THRSH)                                       &
        call error_stop("norm_contr_sgto",                             &
                        "contraction has zero normalization constant", &
                        icontr)
      ! multiplied by \f$\left(\frac{2}{\pi}\right)^{1/4}\frac{1}{\sqrt{(2l+1)!!}}\f$
      norm_const = ang_prefact/norm_const
      ! performs normalization
      do iprim = 1, num_prim
        ! multiplied by \f$(4\alpha_{nl})^{l/2+3/4}\f$
        contr_coef(icontr,iprim) = contr_coef(icontr,iprim)*norm_const &
                                 * (4.0_REALK*exponents(iprim))**exp_power
      end do
    end do
    return
  end subroutine norm_contr_sgto
