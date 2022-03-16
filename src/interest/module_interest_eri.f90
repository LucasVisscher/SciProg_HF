! ------------------------------------------------------------------------------------
!
! Program:      Dirac 
!
! Library:      InteRest
!
! Module:       module_interest_eri.f90 
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
! Comments:     
!                            
! ------------------------------------------------------------------------------------
MODULE module_interest_eri

  use module_interest_osr
  use module_interest_hrr
 
  implicit none

  public interest_eri

  private

  interface interest_eri
    module procedure interest_eri_4c_new
  end interface

CONTAINS

! -------------------------------------------------------------------------
  SUBROUTINE interest_eri_4c_new(class,factor,fint_inp,                          &
                                 la_inp,alpha_inp,ax_inp,ay_inp,az_inp,anorm_inp,&
                                 lb_inp,beta_inp, bx_inp,by_inp,bz_inp,bnorm_inp,&
                                 lc_inp,gamma_inp,cx_inp,cy_inp,cz_inp,cnorm_inp,&
                                 ld_inp,delta_inp,dx_inp,dy_inp,dz_inp,dnorm_inp )

    !output
    real(8), intent(out) :: fint_inp(*)
    !input      
    character(len=4), intent(in) :: class
    real(8),          intent(in) :: factor 
    integer,          intent(in) :: la_inp,lb_inp,lc_inp,ld_inp
    real(8),          intent(in) :: alpha_inp,ax_inp,ay_inp,az_inp,anorm_inp
    real(8),          intent(in) :: beta_inp, bx_inp,by_inp,bz_inp,bnorm_inp
    real(8),          intent(in) :: gamma_inp,cx_inp,cy_inp,cz_inp,cnorm_inp
    real(8),          intent(in) :: delta_inp,dx_inp,dy_inp,dz_inp,dnorm_inp
    !local
    integer :: nint,index
    integer :: la,lb,lc,ld
    real(8), save :: fint(nccx*nccx)
    !$OMP threadprivate(fint)
    real(8) :: alpha,ax,ay,az,anorm
    real(8) :: beta, bx,by,bz,bnorm
    real(8) :: gamma,cx,cy,cz,cnorm
    real(8) :: delta,dx,dy,dz,dnorm
    integer :: n1,n2,nab,ncd
    real(8) :: rxab,ryab,rzab,rrab
    real(8) :: rxcd,rycd,rzcd,rrcd
    real(8) :: pexp,rexp,rxpa,rypa,rzpa
    real(8) :: qexp,sexp,rxqc,ryqc,rzqc
    real(8) :: px,py,pz,qx,qy,qz
    real(8) :: psq,ppq,rho,rxpq,rypq,rzpq
    real(8) :: rxwp,rywp,rzwp,rxwq,rywq,rzwq
    real(8) :: cntr,cntab,cntcd,tval,fnbra,fdbra,fnket,fdket
    real(8), parameter :: pi52 = 2.d0*pi**(5.d0/2.d0) !2*pow(pi,5/2)


    if( la_inp >= lb_inp )then
      la = la_inp    
      lb = lb_inp   
      ax = ax_inp   
      ay = ay_inp   
      az = az_inp    
      bx = bx_inp   
      by = by_inp   
      bz = bz_inp   
      index = 1200
      alpha = alpha_inp
      anorm = anorm_inp
      beta  = beta_inp  
      bnorm = bnorm_inp
    else
      la = lb_inp    
      lb = la_inp   
      ax = bx_inp   
      ay = by_inp   
      az = bz_inp    
      bx = ax_inp   
      by = ay_inp   
      bz = az_inp   
      index = 2100
      alpha = beta_inp
      anorm = bnorm_inp
      beta  = alpha_inp  
      bnorm = anorm_inp
    endif

    if( lc_inp >= ld_inp )then
      lc = lc_inp    
      ld = ld_inp   
      cx = cx_inp   
      cy = cy_inp   
      cz = cz_inp    
      dx = dx_inp   
      dy = dy_inp   
      dz = dz_inp   
      gamma = gamma_inp
      cnorm = cnorm_inp
      delta = delta_inp  
      dnorm = dnorm_inp
      index = index +34 
    else
      lc = ld_inp    
      ld = lc_inp   
      cx = dx_inp   
      cy = dy_inp   
      cz = dz_inp    
      dx = cx_inp   
      dy = cy_inp   
      dz = cz_inp   
      gamma = delta_inp
      cnorm = dnorm_inp
      delta = gamma_inp  
      dnorm = cnorm_inp
      index = index +43 
    endif
      
    rxab = ax - bx
    ryab = ay - by
    rzab = az - bz
    rrab = rxab*rxab + ryab*ryab + rzab*rzab

    pexp = alpha + beta
    if( rrab.lt.1.d-12 )then
      px    = ax 
      py    = ay
      pz    = az
      rxpa  = 0.0d0 
      rypa  = 0.0d0
      rzpa  = 0.0d0
      rxab  = 0.0d0  
      ryab  = 0.0d0  
      rzab  = 0.0d0  
      cntab = anorm*bnorm 
    else
      px    = ( alpha*ax + beta*bx )/pexp
      py    = ( alpha*ay + beta*by )/pexp
      pz    = ( alpha*az + beta*bz )/pexp
      rxpa  = px - ax
      rypa  = py - ay
      rzpa  = pz - az
      rexp  = alpha*beta/pexp
      cntab = anorm*bnorm*dexp(-rexp*rrab)
    endif

    rxcd = cx - dx
    rycd = cy - dy
    rzcd = cz - dz
    rrcd = rxcd*rxcd + rycd*rycd + rzcd*rzcd

    qexp = gamma + delta
    if( rrcd.lt.1.d-12 )then
      qx    = cx 
      qy    = cy
      qz    = cz
      rxqc  = 0.0d0 
      ryqc  = 0.0d0
      rzqc  = 0.0d0
      rxcd  = 0.0d0  
      rycd  = 0.0d0  
      rzcd  = 0.0d0  
      cntcd = cnorm*dnorm 
    else
      qx    = ( gamma*cx + delta*dx )/qexp
      qy    = ( gamma*cy + delta*dy )/qexp
      qz    = ( gamma*cz + delta*dz )/qexp
      rxqc  = qx - cx
      ryqc  = qy - cy
      rzqc  = qz - cz
      sexp  = gamma*delta/qexp
      cntcd = cnorm*dnorm*dexp(-sexp*rrcd)
    endif


    psq = pexp+qexp
    ppq = pexp*qexp
    rho = ppq/psq

    rxpq = px - qx
    rypq = py - qy
    rzpq = pz - qz
    tval = rho*((rxpq*rxpq)+(rypq*rypq)+(rzpq*rzpq))
    cntr = pi52*(1.0d0/dsqrt(psq))*(1.0d0/ppq)*cntab*cntcd*factor

    fnbra = qexp/psq
    fnket = pexp/psq 
    fdbra = 1.0d0/(2.0d0*pexp)
    fdket = 1.0d0/(2.0d0*qexp)
    rxwp  = ( pexp*px + qexp*qx )/psq - px
    rywp  = ( pexp*py + qexp*qy )/psq - py
    rzwp  = ( pexp*pz + qexp*qz )/psq - pz
    rxwq  = ( pexp*px + qexp*qx )/psq - qx
    rywq  = ( pexp*py + qexp*qy )/psq - qy
    rzwq  = ( pexp*pz + qexp*qz )/psq - qz


    nab  = ncc(la)*ncc(lb)
    ncd  = ncc(lc)*ncc(ld)
    nint = nab*ncd 

    select case(class)
      case('llll')
          call interest_osr_class_nuclear(n1,n2,fint,cntr,tval,                                      &
                                          la,(la+lb-1),nab,fnbra,fdbra,rxpa,rypa,rzpa,rxwp,rywp,rzwp,&
                                          lc,(lc+ld-1),ncd,fnket,fdket,rxqc,ryqc,rzqc,rxwq,rywq,rzwq )
          call interest_hrr_ket(fint,ncd,nab,la,lb,rxab,ryab,rzab,n2)
          call interest_hrr_bra(fint,ncd,nab,lc,ld,rxcd,rycd,rzcd) 
          select case( index )
            case( 1234 ); fint_inp(1:nint)=fint(1:nint) 
            case( 2134 ); call transpose_ab( fint, fint_inp, ncc(la), ncc(lb), ncd ) 
            case( 1243 ); call transpose_cd( fint, fint_inp, ncc(lc), ncc(ld), nab ) 
            case( 2143 ); fint_inp(1:nint)=fint(1:nint)
                          call transpose_ab( fint_inp, fint,     ncc(la), ncc(lb), ncd )
                          call transpose_cd( fint,     fint_inp, ncc(lc), ncc(ld), nab )
          end select 
      case default
          write(6,*)'Error: unknown eri class inside InteRest => stop'; stop 
    end select

    CONTAINS

    !---------
    subroutine transpose_ab( A, B, na, nb, ncd ) 
      !input
      integer, intent(in) :: na, nb, ncd 
      real(8), intent(in) :: A(ncd,na,nb)
      !output
      real(8), intent(out) :: B(ncd,nb,na)
      !local
      integer :: i,j,k 
      !
      do j=1,nb
        do i=1,na
          do k=1,ncd
            B(k,j,i) = A(k,i,j) 
          enddo
        enddo
      enddo
    end subroutine
    !-------------
    subroutine transpose_cd( A, B, nc, nd, nab ) 
      !input
      integer, intent(in) :: nc, nd, nab 
      real(8), intent(in) :: A(nc,nd,nab)
      !output
      real(8), intent(out) :: B(nd,nc,nab)
      !local
      integer :: i,j,k 
      !
      do k=1,nab
        do j=1,nd
          do i=1,nc
            B(j,i,k) = A(i,j,k) 
          enddo
        enddo
      enddo
    end subroutine

  END SUBROUTINE

! -------------------------------------------------------------------------
  SUBROUTINE abcd(factor,                       &
                  alpha,ca,ax,ay,az,anorm,      &
                  beta, cb,bx,by,bz,bnorm,      &
                  gamma,cc,cx,cy,cz,cnorm,      &
                  delta,cd,dx,dy,dz,dnorm,      &
                  cntr,tval,                    &
                  fnbra,fdbra,fnket,fdket,      &
                  rxab,ryab,rzab,rxcd,rycd,rzcd,&
                  rxpa,rypa,rzpa,rxwp,rywp,rzwp,&
                  rxqc,ryqc,rzqc,rxwq,rywq,rzwq )

    !-- input --!
    real(8), intent(in) :: factor 
    integer, intent(in) :: ca,cb,cc,cd
    real(8), intent(in) :: alpha,ax,ay,az,anorm
    real(8), intent(in) :: beta, bx,by,bz,bnorm
    real(8), intent(in) :: gamma,cx,cy,cz,cnorm
    real(8), intent(in) :: delta,dx,dy,dz,dnorm
    !-- output --!
    real(8), intent(out) :: rxab,ryab,rzab
    real(8), intent(out) :: rxcd,rycd,rzcd
    real(8), intent(out) :: rxpa,rypa,rzpa
    real(8), intent(out) :: rxqc,ryqc,rzqc
    real(8), intent(out) :: rxwp,rywp,rzwp
    real(8), intent(out) :: rxwq,rywq,rzwq
    real(8), intent(out) :: cntr,tval,fnbra,fdbra,fnket,fdket
    !-- local --!
    logical :: AeqB
    logical :: CeqD
    logical :: ABeqCD
    real(8) :: pexp,rexp
    real(8) :: qexp,sexp
    real(8) :: cntab,cntcd
    real(8) :: px,py,pz,qx,qy,qz
    real(8) :: psq,ppq,rho,rxpq,rypq,rzpq
    real(8), parameter :: pi52 = 2.d0*pi**(5.d0/2.d0) !2*pow(pi,5/2)


    AeqB   = ( ca==cb )
    CeqD   = ( cc==cd )
    ABeqCD = ( AeqB .and. CeqD .and. ca==cc )

    pexp=alpha+beta
    if( AeqB )then
      px    = ax 
      py    = ay
      pz    = az
      rxpa  = 0.0d0 
      rypa  = 0.0d0
      rzpa  = 0.0d0
      rxab  = 0.0d0  
      ryab  = 0.0d0  
      rzab  = 0.0d0  
      cntab = anorm*bnorm 
      AeqB  = .true.
    else
      px    = ( alpha*ax + beta*bx )/pexp
      py    = ( alpha*ay + beta*by )/pexp
      pz    = ( alpha*az + beta*bz )/pexp
      rxpa  = px - ax
      rypa  = py - ay
      rzpa  = pz - az
      rxab  = ax - bx
      ryab  = ay - by
      rzab  = az - bz
      rexp  = alpha*beta/pexp
      cntab = anorm*bnorm*dexp(-rexp*(rxab*rxab+ryab*ryab+rzab*rzab))
      AeqB  = .true.
    endif


    qexp=gamma+delta
    if( CeqD )then
      qx    = cx 
      qy    = cy
      qz    = cz
      rxqc  = 0.0d0 
      ryqc  = 0.0d0
      rzqc  = 0.0d0
      rxcd  = 0.0d0  
      rycd  = 0.0d0  
      rzcd  = 0.0d0  
      cntcd = cnorm*dnorm 
    else
      qx    = ( gamma*cx + delta*dx )/qexp
      qy    = ( gamma*cy + delta*dy )/qexp
      qz    = ( gamma*cz + delta*dz )/qexp
      rxqc  = qx - cx
      ryqc  = qy - cy
      rzqc  = qz - cz
      rxcd  = cx - dx
      rycd  = cy - dy
      rzcd  = cz - dz
      sexp  = gamma*delta/qexp
      cntcd = cnorm*dnorm*dexp(-sexp*(rxcd*rxcd+rycd*rycd+rzcd*rzcd))
    endif


    psq   = pexp+qexp
    ppq   = pexp*qexp
    rho   = ppq/psq
    fnbra = qexp/psq
    fnket = pexp/psq 
    fdbra = 1.0d0/(2.0d0*pexp)
    fdket = 1.0d0/(2.0d0*qexp)
    if( ABeqCD )then
      tval = 0.0d0 
      rxwp = 0.0d0 
      rywp = 0.0d0
      rzwp = 0.0d0
      rxwq = 0.0d0
      rywq = 0.0d0
      rzwq = 0.0d0
      cntr = pi52*(1.0d0/dsqrt(psq))*(1.0d0/ppq)*cntab*cntcd*factor
    else
      rxpq = px - qx
      rypq = py - qy
      rzpq = pz - qz
      tval = rho*((rxpq*rxpq)+(rypq*rypq)+(rzpq*rzpq))
      cntr = pi52*(1.0d0/dsqrt(psq))*(1.0d0/ppq)*cntab*cntcd*factor
      rxwq = ( pexp*px + qexp*qx )/psq
      rywq = ( pexp*py + qexp*qy )/psq
      rzwq = ( pexp*pz + qexp*qz )/psq
      rxwp = rxwq - px
      rywp = rywq - py
      rzwp = rzwq - pz
      rxwq = rxwq - qx
      rywq = rywq - qy
      rzwq = rzwq - qz
    endif

  END SUBROUTINE
! -------------------------------------------------------------------------
END MODULE
