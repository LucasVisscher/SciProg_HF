! ------------------------------------------------------------------------------------
!
! Program:      Dirac 
!
! Library:      InteRest
!
! Module:       module_interest_one.f90 
!
! Description:  InteRest library routine 
!
! Contains:     
!
! Licensing:    Copyright (c) by the authors of DIRAC.
!
!               This program is free software; you can redistribute it and/or
!               modify it under the terms of the GNU Lesser General Public
!               License version 2.1 as published by the Free Software Foundation.
!
!               This program is distributed in the hope that it will be useful,
!               but WITHOUT ANY WARRANTY; without even the implied warranty of
!               MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!               Lesser General Public License for more details.
!
!               If a copy of the GNU LGPL v2.1 was not distributed with this
!               code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
! Author:       Michal Repisky (michal.repisky@uit.no)
!
! Revision:     2.0   
!
! ------------------------------------------------------------------------------------
  subroutine interest_initialize( info )

    use module_interest_osr
    implicit none
    !> input
    logical :: info

    if( is_interest_initialized )then
      return
     !write(6,*)
     !write(6,'(2x,a)')'Integral module InteRest was already initialized'
     !write(6,*)
    else
      if( info )then
        write(6,*)
        write(6,'(2x,a)')'Integral module InteRest is initialized'
        call print_interest_git_revision_info
        write(6,*)
      endif
      call interest_osr_initialize()
      is_interest_initialized = .true.
    endif

  end subroutine
