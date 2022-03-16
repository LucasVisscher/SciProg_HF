!radovan: this is more elegant but is unfortunately not standard
!   ! inverse of constant Pi
!   real(REALK), parameter :: INV_PI = 1.0_REALK/PI
!   ! normalization constant of s-shell spherical Gaussians
!   real(REALK), parameter :: NRM_S_SGTO = sqrt(INV_PI)*0.5_REALK
!   ! normalization constant of p-shell spherical Gaussians
!   real(REALK), parameter :: NRM_P_SGTO = sqrt(3.0_REALK)*NRM_S_SGTO
!   ! normalization constant of d-shell spherical Gaussians
!   real(REALK), parameter :: NRM_D_SGTO(3) = (/                &
!     sqrt(15.0_REALK)*NRM_S_SGTO, sqrt(3.75_REALK)*NRM_S_SGTO, &
!     sqrt(1.25_REALK)*NRM_S_SGTO/)

    real(REALK), parameter :: INV_PI = 0.318309886183791_REALK
    real(REALK), parameter :: NRM_S_SGTO = 0.282094791773878_REALK
    real(REALK), parameter :: NRM_P_SGTO = 0.488602511902920_REALK
    real(REALK), parameter :: NRM_D_SGTO(3) = (/1.09254843059208_REALK,  &
                                                0.546274215296040_REALK, &
                                                0.315391565252520_REALK/)
