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
!!  This file transforms the integrals with Hermite/Cartesian GTOs to those
!!  with spherical GTOs.
!!
!!  2011-08-27, Bin Gao:
!!  * adds subroutine hgto_to_sgto_oop to transforms the integrals out of place
!!
!!  2009-07-16, Bin Gao:
!!  * modified from Andreas' subroutine cart_to_spher

#include "stdout.h"

  !> \brief transforms the integrals with Hermite/Cartesian GTOs to those
  !>        with spherical GTOs, in place
  !> \author Bin Gao
  !> \date 2009-07-16
  !> \param angular_ket is the angular number on ket center
  !> \param num_sgto_bra is the number of spherical Gaussians on bra center
  !> \param num_hgto_ket is the number of Hermite Gaussians on ket center
  !> \param num_opt is the number of other operators
  !> \return hgto_ints contains the Hermite integrals on entry, while it contains
  !>         the transformed spherical integrals on exit
  subroutine hgto_to_sgto_in_situ(angular_ket, num_sgto_bra, num_hgto_ket, &
                                  num_opt, hgto_ints)
    use xkind
    implicit none
    integer, intent(in) :: angular_ket
    integer, intent(in) :: num_sgto_bra
    integer, intent(in) :: num_hgto_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(inout) :: hgto_ints(num_sgto_bra,num_hgto_ket,num_opt)
!f2py intent(in) :: angular_ket
!f2py intent(hide) :: num_sgto_bra
!f2py intent(hide) :: num_hgto_ket
!f2py intent(hide) :: num_opt
!f2py intent(inout) :: hgto_ints
#include "private/pi.h"
#include "private/hgto_to_sgto.h"
    real(REALK), allocatable :: herm_px_ints(:)  !x-component of p-shell Hermite integrals
    real(REALK) legendre                         !Legendre coefficient
    real(REALK), allocatable :: invnorm(:)       !normalization constants
    integer, allocatable :: binomial(:)          !binomial coefficients for (x+iy)^?
    integer, allocatable :: binomixy(:)          !binomial coefficients for (xx+yy)^?
    integer iaddr, jaddr
    integer iang, jang                           !incremental recorders over angular number
    integer imag                                 !incremental recorder over magnetic quantum number
    integer iopt                                 !incremental recorder over operators
    integer ierr                                 !error information
#if defined(XTIME)
    real(REALK) curr_time                        !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    !-----hardcode for s-, p- and d-shells-----
    if (angular_ket==0) then
      hgto_ints(:,1,:) = NRM_S_SGTO*hgto_ints(:,1,:)
      return
    else if (angular_ket==1) then
      allocate(herm_px_ints(num_sgto_bra), stat=ierr)
      if (ierr/=0)                              &
        call error_stop("hgto_to_sgto_in_situ", &
                        "failed to allocate herm_px_ints", num_sgto_bra)
      ! in the order of -1(y), 0(z), 1(x) for p-shells
      do iopt = 1, num_opt
        herm_px_ints = hgto_ints(:,1,iopt)
        hgto_ints(:,1,iopt) = NRM_P_SGTO*hgto_ints(:,2,iopt)
        hgto_ints(:,2,iopt) = NRM_P_SGTO*hgto_ints(:,3,iopt)
        hgto_ints(:,3,iopt) = NRM_P_SGTO*herm_px_ints
      end do
      deallocate(herm_px_ints)
      return
    else if (angular_ket==2) then
      hgto_ints(:,6,:) = 2.0_REALK*hgto_ints(:,6,:) &       !2zz-xx-yy
                       - hgto_ints(:,1,:)-hgto_ints(:,3,:)
      hgto_ints(:,3,:) = hgto_ints(:,1,:)-hgto_ints(:,3,:)  !xx-yy
      hgto_ints(:,1,:) = NRM_D_SGTO(1)*hgto_ints(:,2,:)     !xy
      hgto_ints(:,2,:) = NRM_D_SGTO(1)*hgto_ints(:,5,:)     !yz
      hgto_ints(:,5,:) = NRM_D_SGTO(2)*hgto_ints(:,3,:)     !normalizes xx-yy
      hgto_ints(:,3,:) = NRM_D_SGTO(3)*hgto_ints(:,6,:)     !normalizes 2zz-xx-yy
      hgto_ints(:,4,:) = NRM_D_SGTO(1)*hgto_ints(:,4,:)     !xz
      return
    end if
!FIXME: adds more comments, and considers improving efficiency?
    !-----higher angular momentum than 2------
    allocate(binomial(angular_ket+1), stat=ierr)
    if (ierr/=0)                              &
      call error_stop("hgto_to_sgto_in_situ", &
                      "failed to allocate binomial", angular_ket+1)
    allocate(binomixy(angular_ket+1), stat=ierr)
    if (ierr/=0)                              &
      call error_stop("hgto_to_sgto_in_situ", &
                      "failed to allocate binomixy", angular_ket+1)
    allocate(invnorm(angular_ket+1), stat=ierr)
    if (ierr/=0)                              &
      call error_stop("hgto_to_sgto_in_situ", &
                      "failed to allocate invnorm", angular_ket+1)
    !--first m=0, which only has real part
    legendre = 1.0_REALK
    binomial(1) = 1
    binomial(2:(angular_ket+2)/2) = 0
    iaddr = (angular_ket+1)*(angular_ket+2)/2
    jaddr = iaddr
    do iang = 2, angular_ket, 2
      ! increments addition coefficient
      jang = angular_ket-iang+1
      legendre = -legendre*real(jang*(jang+1),REALK)/real(iang*iang,REALK)
      ! jumps two rows up the triangle
      jaddr = jaddr-2*iang-1
      ! increments binomial coefficients
      jang = iang/2
      binomial(2:jang+1) = binomial(2:jang+1)+binomial(1:jang)
      do jang = iang, 0, -2
        hgto_ints(:,iaddr,:) = hgto_ints(:,iaddr,:)      &
                             + hgto_ints(:,jaddr+jang,:) &
                             * real(binomial(jang/2+1),REALK)*legendre
      end do
    end do
    !--then 0<m<=angular_ket, which have both real and imaginary parts
    binomixy(1) = 1
    binomixy(2:) = 0
    ! for m=0
    invnorm(1) = sqrt(real(4*angular_ket+2,REALK))*NRM_S_SGTO
    do imag = 1, angular_ket
      ! jumps one row up the triangle
      iaddr = iaddr-imag-1
      ! increments starting binomial coefficients to the next m
      binomial = binomixy
      binomial(2:imag+1:2) = binomial(2:imag+1:2)+binomixy(1:imag:2)
      binomial(3:imag+1:2) = binomial(3:imag+1:2)-binomixy(2:imag:2)
      binomixy(1:imag+1) = binomial(1:imag+1)
      ! first row is special
      if (binomial(2)/=1)                               &
        hgto_ints(:,iaddr+1,:) = hgto_ints(:,iaddr+1,:) &
                               * real(binomial(2),REALK)
      do iang = imag, 2, -1
        hgto_ints(:,iaddr+mod(iang,2),:) = hgto_ints(:,iaddr+mod(iang,2),:) &
                                         + hgto_ints(:,iaddr+iang,:)        &
                                         * real(binomial(iang+1),REALK)
      end do
      ! loops over remaining rows
      legendre = 1.0_REALK
      jaddr = iaddr
      do iang = imag+2, angular_ket, 2
        ! increments addition coefficient
        jang = angular_ket-iang+1
        legendre = -legendre*real(jang*(jang+1),REALK)/real((iang-imag)*(iang+imag),REALK)
        ! jumps two rows up the triangle
        jaddr = jaddr-2*iang-1
        ! increments binomial coefs
        binomial(3:iang+1) = binomial(3:iang+1)+binomial(1:iang-1)
        do jang = iang, 0, -1
          hgto_ints(:,iaddr+mod(jang,2),:) = hgto_ints(:,iaddr+mod(jang,2),:) &
                                           + hgto_ints(:,jaddr+jang,:)        &
                                           * real(binomial(jang+1),REALK)*legendre
        end do
      end do
      ! calculates the normalization constant
      invnorm(imag+1) = invnorm(imag)                                             &
                      * sqrt(real((angular_ket+1-imag)*(angular_ket+imag),REALK)) &
                      / real(2*imag,REALK)
    end do
    deallocate(binomial)
    deallocate(binomixy)
    !-----normalizes and reorders-----
    hgto_ints(:,2*angular_ket+1,:) = hgto_ints(:,1,:)*invnorm(angular_ket+1)  !re Yll
    hgto_ints(:,1,:) = hgto_ints(:,2,:)*invnorm(angular_ket+1)                !im Yll
    iaddr = angular_ket+2
    do imag = angular_ket-1, 1, -1
      hgto_ints(:,angular_ket+1-imag,:) = hgto_ints(:,iaddr+1,:)*invnorm(imag+1)  !im Ylm
      hgto_ints(:,angular_ket+1+imag,:) = hgto_ints(:,iaddr,:)*invnorm(imag+1)    !re Ylm
      iaddr = iaddr+imag+1
    end do
    hgto_ints(:,angular_ket+1,:) = hgto_ints(:,iaddr,:)*invnorm(1)/sqrt(2.0_REALK)  !Yl0
    deallocate(invnorm)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "hgto_to_sgto_in_situ", STDOUT)
#endif
    return
  end subroutine hgto_to_sgto_in_situ

  !> \brief transforms the integrals with Hermite/Cartesian GTOs to those
  !>        with spherical GTOs, out of place
  !> \author Bin Gao
  !> \date 2011-08-27
  !> \param angular_ket is the angular number on ket center
  !> \param num_sgto_bra is the number of spherical Gaussians on bra center
  !> \param num_hgto_ket is the number of Hermite Gaussians on ket center
  !> \param num_opt is the number of other operators
  !> \param hgto_ints contains the Hermite integrals
  !> \param num_sgto_ket is the number of spherical Gaussians on ket center
  !> \return sgto_ints contains the transformed spherical integrals
  subroutine hgto_to_sgto_oop(angular_ket, num_sgto_bra, num_hgto_ket, &
                              num_opt, hgto_ints, num_sgto_ket, sgto_ints)
    use xkind
    implicit none
    integer, intent(in) :: angular_ket
    integer, intent(in) :: num_sgto_bra
    integer, intent(in) :: num_hgto_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(inout) :: hgto_ints(num_sgto_bra,num_hgto_ket,num_opt)
    integer, intent(in) :: num_sgto_ket
    real(REALK), intent(out) :: sgto_ints(num_sgto_bra,num_sgto_ket,num_opt)
!f2py intent(in) :: angular_ket
!f2py intent(hide) :: num_sgto_bra
!f2py intent(hide) :: num_hgto_ket
!f2py intent(hide) :: num_opt
!f2py intent(inout) :: hgto_ints
!f2py intent(in) :: num_sgto_ket
!f2py intent(out) :: sgto_ints
!f2py depend(num_sgto_bra) :: sgto_ints
!f2py depend(num_sgto_ket) :: sgto_ints
!f2py depend(num_opt) :: sgto_ints
#include "private/pi.h"
#include "private/hgto_to_sgto.h"
    real(REALK) legendre                         !Legendre coefficient
    real(REALK), allocatable :: invnorm(:)       !normalization constants
    integer, allocatable :: binomial(:)          !binomial coefficients for (x+iy)^?
    integer, allocatable :: binomixy(:)          !binomial coefficients for (xx+yy)^?
    integer iaddr, jaddr
    integer iang, jang                           !incremental recorders over angular number
    integer imag                                 !incremental recorder over magnetic quantum number
    integer iopt                                 !incremental recorder over operators
    integer ierr                                 !error information
#if defined(XTIME)
    real(REALK) curr_time                        !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    !-----hardcode for s-, p- and d-shells-----
    if (angular_ket==0) then
      sgto_ints(:,1,:) = NRM_S_SGTO*hgto_ints(:,1,:)
      return
    else if (angular_ket==1) then
      ! in the order of -1(y), 0(z), 1(x) for p-shells
      do iopt = 1, num_opt
        sgto_ints(:,1,iopt) = NRM_P_SGTO*hgto_ints(:,2,iopt)
        sgto_ints(:,2,iopt) = NRM_P_SGTO*hgto_ints(:,3,iopt)
        sgto_ints(:,3,iopt) = NRM_P_SGTO*hgto_ints(:,1,iopt)
      end do
      return
    else if (angular_ket==2) then
      hgto_ints(:,6,:) = 2.0_REALK*hgto_ints(:,6,:) &       !2zz-xx-yy
                       - hgto_ints(:,1,:)-hgto_ints(:,3,:)
      hgto_ints(:,3,:) = hgto_ints(:,1,:)-hgto_ints(:,3,:)  !xx-yy
      sgto_ints(:,1,:) = NRM_D_SGTO(1)*hgto_ints(:,2,:)    !xy
      sgto_ints(:,2,:) = NRM_D_SGTO(1)*hgto_ints(:,5,:)    !yz
      sgto_ints(:,5,:) = NRM_D_SGTO(2)*hgto_ints(:,3,:)    !normalizes xx-yy
      sgto_ints(:,3,:) = NRM_D_SGTO(3)*hgto_ints(:,6,:)    !normalizes 2zz-xx-yy
      sgto_ints(:,4,:) = NRM_D_SGTO(1)*hgto_ints(:,4,:)    !xz
      return
    end if
!FIXME: adds more comments, and considers improving efficiency?
    !-----higher angular momentum than 2------
    allocate(binomial(angular_ket+1), stat=ierr)
    if (ierr/=0)                          &
      call error_stop("hgto_to_sgto_oop", &
                      "failed to allocate binomial", angular_ket+1)
    allocate(binomixy(angular_ket+1), stat=ierr)
    if (ierr/=0)                          &
      call error_stop("hgto_to_sgto_oop", &
                      "failed to allocate binomixy", angular_ket+1)
    allocate(invnorm(angular_ket+1), stat=ierr)
    if (ierr/=0)                          &
      call error_stop("hgto_to_sgto_oop", &
                      "failed to allocate invnorm", angular_ket+1)
    !--first m=0, which only has real part
    legendre = 1.0_REALK
    binomial(1) = 1
    binomial(2:(angular_ket+2)/2) = 0
    iaddr = (angular_ket+1)*(angular_ket+2)/2
    jaddr = iaddr
    do iang = 2, angular_ket, 2
      ! increments addition coefficient
      jang = angular_ket-iang+1
      legendre = -legendre*real(jang*(jang+1),REALK)/real(iang*iang,REALK)
      ! jumps two rows up the triangle
      jaddr = jaddr-2*iang-1
      ! increments binomial coefficients
      jang = iang/2
      binomial(2:jang+1) = binomial(2:jang+1)+binomial(1:jang)
      do jang = iang, 0, -2
        hgto_ints(:,iaddr,:) = hgto_ints(:,iaddr,:)      &
                             + hgto_ints(:,jaddr+jang,:) &
                             * real(binomial(jang/2+1),REALK)*legendre
      end do
    end do
    !--then 0<m<=angular_ket, which have both real and imaginary parts
    binomixy(1) = 1
    binomixy(2:) = 0
    ! for m=0
    invnorm(1) = sqrt(real(4*angular_ket+2,REALK))*NRM_S_SGTO
    do imag = 1, angular_ket
      ! jumps one row up the triangle
      iaddr = iaddr-imag-1
      ! increments starting binomial coefficients to the next m
      binomial = binomixy
      binomial(2:imag+1:2) = binomial(2:imag+1:2)+binomixy(1:imag:2)
      binomial(3:imag+1:2) = binomial(3:imag+1:2)-binomixy(2:imag:2)
      binomixy(1:imag+1) = binomial(1:imag+1)
      ! first row is special
      if (binomial(2)/=1)                               &
        hgto_ints(:,iaddr+1,:) = hgto_ints(:,iaddr+1,:) &
                               * real(binomial(2),REALK)
      do iang = imag, 2, -1
        hgto_ints(:,iaddr+mod(iang,2),:) = hgto_ints(:,iaddr+mod(iang,2),:) &
                                         + hgto_ints(:,iaddr+iang,:)        &
                                         * real(binomial(iang+1),REALK)
      end do
      ! loops over remaining rows
      legendre = 1.0_REALK
      jaddr = iaddr
      do iang = imag+2, angular_ket, 2
        ! increments addition coefficient
        jang = angular_ket-iang+1
        legendre = -legendre*real(jang*(jang+1),REALK)/real((iang-imag)*(iang+imag),REALK)
        ! jumps two rows up the triangle
        jaddr = jaddr-2*iang-1
        ! increments binomial coefs
        binomial(3:iang+1) = binomial(3:iang+1)+binomial(1:iang-1)
        do jang = iang, 0, -1
          hgto_ints(:,iaddr+mod(jang,2),:) = hgto_ints(:,iaddr+mod(jang,2),:) &
                                           + hgto_ints(:,jaddr+jang,:)        &
                                           * real(binomial(jang+1),REALK)*legendre
        end do
      end do
      ! calculates the normalization constant
      invnorm(imag+1) = invnorm(imag)                                             &
                      * sqrt(real((angular_ket+1-imag)*(angular_ket+imag),REALK)) &
                      / real(2*imag,REALK)
    end do
    deallocate(binomial)
    deallocate(binomixy)
    !-----normalizes and reorders-----
    sgto_ints(:,2*angular_ket+1,:) = hgto_ints(:,1,:)*invnorm(angular_ket+1)  !re Yll
    sgto_ints(:,1,:) = hgto_ints(:,2,:)*invnorm(angular_ket+1)                !im Yll
    iaddr = angular_ket+2
    do imag = angular_ket-1, 1, -1
      sgto_ints(:,angular_ket+1-imag,:) = hgto_ints(:,iaddr+1,:)*invnorm(imag+1)  !im Ylm
      sgto_ints(:,angular_ket+1+imag,:) = hgto_ints(:,iaddr,:)*invnorm(imag+1)    !re Ylm
      iaddr = iaddr+imag+1
    end do
    sgto_ints(:,angular_ket+1,:) = hgto_ints(:,iaddr,:)*invnorm(1)/sqrt(2.0_REALK)  !Yl0
    deallocate(invnorm)
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "hgto_to_sgto_oop", STDOUT)
#endif
    return
  end subroutine hgto_to_sgto_oop
