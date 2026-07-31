 !***********************************************************************
 !
 !    Copyright (c) 2016, Lawrence Livermore National Security, LLC.
 !                        Produced at the Lawrence Livermore National
 !                        Laboratory.
 !                        Written by Nicolas Schunck, schunck1@llnl.gov
 !
 !    LLNL-CODE-710577 All rights reserved.
 !    LLNL-CODE-470611 All rights reserved.
 !
 !    Copyright 2016, N. Schunck, J. Dobaczewski, W. Satula, P. Baczyk,
 !                    J. Dudek, Y. Gao, M. Konieczka, K. Sato, Y. Shi,
 !                    X.B. Wang, T.R. Werner
 !    Copyright 2012, N. Schunck, J. Dobaczewski, J. McDonnell,
 !                    W. Satula, J.A. Sheikh, A. Staszczak,
 !                    M. Stoitsov, P. Toivanen
 !    Copyright 2009, J. Dobaczewski, W. Satula, B.G. Carlsson, J. Engel,
 !                    P. Olbratowski, P. Powalowski, M. Sadziak,
 !                    J. Sarich, N. Schunck, A. Staszczak, M. Stoitsov,
 !                    M. Zalewski, H. Zdunczuk
 !    Copyright 2004, 2005, J. Dobaczewski, P. Olbratowski
 !    Copyright 1997, 2000, J. Dobaczewski, J. Dudek
 !
 !    This file is part of HFODD.
 !
 !    HFODD is free software: you can redistribute it and/or modify it
 !    under the terms of the GNU General Public License as published by
 !    the Free Software Foundation, either version 3 of the License, or
 !    (at your option) any later version.
 !
 !    HFODD is distributed in the hope that it will be useful, but
 !    WITHOUT ANY WARRANTY; without even the implied warranty of
 !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 !    GNU General Public License for more details.
 !
 !    You should have received a copy of the GNU General Public License
 !    along with HFODD. If not, see <http://www.gnu.org/licenses/>.
 !
 !    OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
 !    LICENSE
 !
 !    Our Preamble Notice
 !
 !      A. This notice is required to be provided under our contract
 !         with the U.S. Department of Energy (DOE). This work was
 !         produced at the Lawrence Livermore National Laboratory under
 !         Contract No. DE-AC52-07NA27344 with the DOE.
 !      B. Neither the United States Government nor Lawrence Livermore
 !         National Security, LLC nor any of their employees, makes any
 !         warranty, express or implied, or assumes any liability or
 !         responsibility for the accuracy, completeness, or usefulness
 !         of any information, apparatus, product, or process disclosed,
 !         or represents that its use would not infringe privately-owned
 !         rights.
 !      C. Also, reference herein to any specific commercial products,
 !         process, or services by trade name, trademark, manufacturer
 !         or otherwise does not necessarily constitute or imply its
 !         endorsement, recommendation, or favoring by the United States
 !         Government or Lawrence Livermore National Security, LLC. The
 !         views and opinions of authors expressed herein do not
 !         necessarily state or reflect those of the United States
 !         Government or Lawrence Livermore National Security, LLC, and
 !         shall not be used for advertising or product endorsement
 !         purposes.
 !
 !    The precise terms and conditions for copying, distribution and
 !    modification are contained in the file COPYING.
 !
 !***********************************************************************

MODULE PARS
 REAL::PI=ATAN(1.0D0)*4
END MODULE PARS
!=======================================================================
MODULE PFAFFI
!routines to calculate overlap with pfaffian method
USE PARS
USE BLOSAV
USE hfodd_sizes
INTEGER::LTIEXT,LUPEXT,LSTEXT,LTOEXT
COMPLEX::OVKTMP
INTEGER,DIMENSION(0:NDISOS)::IODDS
REAL   ,DIMENSION(0:NDISOS)::LV2AUX,PV2AUX
COMPLEX,ALLOCATABLE::VUAUX(:,:,:)
!COMPLEX,DIMENSION(1:2*NDSTAT,1:2*NDSTAT,0:NDISOS)::VUAUX

COMPLEX,ALLOCATABLE::WSINGL(:,:,:),WBAZER(:,:,:,:),&
                     WALTMP(:,:,:),WARTMP(:,:,:),&
                     SALTMP(:,:,:,:),SALODD(:,:,:,:),SARODD(:,:,:,:)
CONTAINS
!=======================================================================
!CALCULATE OVERLAP BY USING PFAFFIAN
!=======================================================================
 subroutine norm_n(iswitch,ihfb,lupper,ldbase,waleft,warigh,vu,logv2,phav2,norm)
 implicit none
!calculate norm in normal and canonical basis
!return logv2, or norm
 integer::iswitch,ihfb,mrever,isimpy,ltimup,lupper,ldbase
 complex,dimension(1:ndbase,1:4*ndstat,0:ndspin)::waleft,warigh
 complex,dimension(1:2*ndstat,1:2*ndstat)::vu
 real   ::logv2,phav2
 complex::norm
 !
 integer::signs,istate,jstate,ibase,ispin,iflag,ih,dims,npi
 real   ::logpf,phase,netphase
 complex::pf,c_i,c_zero,c_unit
 complex,dimension(1:lupper,1:lupper)::vrv
 complex,dimension(1:2*lupper,1:2*lupper)::pfmat,work
 complex::zdotu
 !
 c_i=cmplx(0.0d0,1.0d0)
 c_zero=cmplx(0.0d0,0.0d0)
 c_unit=cmplx(1.0d0,0.0d0)
 vrv(:,:)=c_zero
 pfmat(:,:)=c_zero
 work(:,:)=c_zero

 if(iswitch.eq.0)then
!(v^t.u)i>j. for odd nuclei, the last line of v^t.u is not skew
    vu(:,:)=c_zero
    if(ihfb.ge.1)then
       do istate=1,lupper
          vu(istate,istate)=c_zero
          do jstate=istate+1,lupper
             vu(istate,jstate)&
      =conjg(zdotu(ldbase,waleft(1:ldbase,istate,0),1,waleft(1:ldbase,jstate+lupper,1),1)&
            -zdotu(ldbase,waleft(1:ldbase,istate,1),1,waleft(1:ldbase,jstate+lupper,0),1))
             vu(jstate,istate)=-vu(istate,jstate)
          enddo
       enddo
    endif
 endif
!VRV
 CALL ZGEMM('C','N',LUPPER,LUPPER,LDBASE,C_UNIT,&
             WALEFT(1:LDBASE,1:LUPPER,0),LDBASE,&
             WARIGH(1:LDBASE,1:LUPPER,0),LDBASE,C_ZERO,VRV,LUPPER)
 CALL ZGEMM('C','N',LUPPER,LUPPER,LDBASE,C_UNIT,&
             WALEFT(1:LDBASE,1:LUPPER,1),LDBASE,&
             WARIGH(1:LDBASE,1:LUPPER,1),LDBASE,C_UNIT,VRV,LUPPER)

!construct the skew matrix
 pfmat(1:lupper,1:lupper)=vu(1:lupper,1:lupper)
 pfmat(1:lupper,lupper+1:2*lupper)=vrv(1:lupper,1:lupper)
 pfmat(lupper+1:2*lupper,1:lupper)=-transpose(vrv(1:lupper,1:lupper))
 pfmat(lupper+1:2*lupper,lupper+1:2*lupper)=transpose(conjg(vu(1:lupper,1:lupper)))
!calculate pfaffian
 dims=2*lupper
 ih=2
 signs=1
 iflag=1
 if(iswitch.eq.0)then!calculate Logv2
    logv2=0.0d0
    phav2=0.0d0
    call ZPfaffianH (pfmat,dims,dims,pf,logv2,phav2,signs,iflag,work,ih)
    if(signs.eq.(-1)) phav2=phav2+pi
 else!calculate Log(pfaffian) and phase factor
    logpf=0.0d0
    phase=0.0d0
    call ZPfaffianH (pfmat,dims,dims,pf,logpf,phase,signs,iflag,work,ih)
    if(signs.eq.(-1)) phase=phase+pi
!
    npi=int((phase-phav2)/pi)
    netphase=phase-phav2-npi*pi
    if(abs(netphase).le.1.0d-8) then
       netphase=0.0d0
    elseif(abs(abs(netphase)-pi).le.1.0d-8) then
       netphase=pi
    elseif(abs(netphase-pi/2).le.1.0d-8)then
       netphase=pi/2
    elseif(abs(netphase+pi/2).le.1.0d-8)then
       netphase=-pi/2
    endif
!
!    write(*,*)npi,(phase-phav2)/pi,npi*pi,netphase
    norm=dexp(logpf-logv2)*cdexp(c_i*(netphase))*(1-2*abs(mod(npi,2)))
!    norm=dexp(logpf-logv2)*cmplx(dcos(netphase),dsin(netphase))*(1-2*abs(mod(npi,2)))
    if(iflag.eq.0) norm=0.0d0
 endif
 return
 end subroutine norm_n
!=======================================================================
Subroutine ZPfaffianH (SK,LDS,N,Pf,Logpf,phase,SIGNS,IFLAG,WORK,IH)
Implicit none
Integer IS,M,i,j,k,IDEN,IH
Integer,intent(in)::LDS,N
Real epsln,xmod,u2
Complex ki,one,zero,al
Complex,intent(inout)::Pf
REAL   ::logpf,phase
!=====
integer,intent(inout)::signs,IFLAG
!=====
Complex,intent(inout),dimension(LDS,N):: SK
Complex,intent(inout),dimension(LDS,N):: WORK
Complex ZDOTC
one = dcmplx(1.0d+00,0.0d+00)
zero= dcmplx(0.0d+00,0.0d+00)
epsln = 1.0d-13   ! smallest number such that 1+epsln=1
if(IH.ne.1.and.IH.ne.2) then
   write(6,'(" Value of IH undefined in pfaffianH",I2)') IH
end if
if(mod(N,2)==1) then
!                                                   N odd   Pf=0
Pf = zero
!=====
Logpf=0.0d0
phase=0.0d0
signs=1
IFLAG=0
!=====
else
!                                                   N even
if (N==2) then
!                          The pfaffian of a 2x2 matrix is Pf=Sk(1,2)
Pf = Sk(1,2)
Logpf=dlog(cdabs(Sk(1,2)))
phase=datan2(aimag(sk(1,2)),real(sk(1,2)))
signs=1
IFLAG=1
else   ! N.gt.2
!
!                                      N even and > 2
Pf = one
LogPf=0.0d0
phase=0.0d0
signs=1
IFLAG=1
IDen=0
do IS=1,N-2,IH
M = N - IS
! Column N-IS+1 of Sk is copied into column IS of working vector
! This is vector x with dimension M of the Householder transformation
 call zcopy(M,Sk(1,N-IS+1),1,WORK(1,IS),1)
! Norm of x
xmod = dsqrt( abs( ZDOTC(M,WORK(1,IS),1,WORK(1,IS),1) ) ) ! |x|
! if Norm of x=0 the pfaffian is zero
if (xmod>epsln) then
if(abs(Sk(M,M+1))>epsln) then
!   u = x -+ sign(Sk(N-IS,N-IS+1))*xmod e (N-IS)
!   e (M) is the unit cartesian vector with 1 in position M and
!   zero elsewhere
ki    = xmod*Sk(M,M+1)/abs(Sk(M,M+1))
else
ki    = dcmplx(xmod,0.0d+00)
end if
!
!      u = x -+ e**[i*arg(Sk(N-IS,N-IS+1))]*xmod e (N-IS)
!
work(M,IS) = work(M,IS) - ki
!   Norm of u
u2 = abs(ZDOTC(M,WORK(1,IS),1,WORK(1,IS),1))  ! |u|**2
!  If the norm of u is too small the negative sign in the definition
!  of u is taken
if (u2 < epsln**2) then
!  ki changes sign
ki    = -ki
work(M,IS) = work(M,IS) - 2.0d+00*ki
u2 = abs(ZDOTC(M,WORK(1,IS),1,WORK(1,IS),1))  ! |u|**2
end if
!     *
!    u  --> work(1,N-1)
do i=1,M
   work( i , N-1 ) = dconjg(work( i , IS ))
end do
!                *
!    v = A(N-1) u --> work(1,N)
 call Zgemv('N',M,M,one,Sk,LDS,work(1,N-1),1,zero,work(1,N),1)
!      Update of A(N-IS)
al = dcmplx(2.0d+00/u2,0.0d+00)
!     Sk -> Sk + al*u v'
 call zgeru(M,M, al,work(1,IS),1,work(1,N),1,Sk,LDS)
!     Sk -> Sk - al*v u'
 call zgeru(M,M,-al,work(1,N),1,work(1,IS),1,Sk,LDS)
else ! xmod.gt.epsln
ki = 0.0d+00
Iden = Iden + 1
end if   ! xmod.gt.epsln
if(mod(IS,2)==1) then
   Pf = Pf * ki
   if(abs(ki).gt.1.0d-32)then
      Logpf=Logpf+dlog(cdabs(ki))
      phase=phase+datan2(aimag(ki),real(ki))
   else
      Logpf=0
      IFLAG=0!in case ki=0, log(ki) will be illegal
   endif
endif
end do ! IS
Pf = Pf * Sk(1,2) * dfloat(1-2*mod(Iden,2))
if(abs(Sk(1,2)).gt.1.0d-32)then
   Logpf=logpf+dlog(cdabs(Sk(1,2)))
   phase=phase+datan2(aimag(sk(1,2)),real(sk(1,2)))
   signs=(1-2*mod(Iden,2))
else
   Logpf=0
   IFLAG=0!in case ki=0, log(ki) will be illegal
endif
if (IH==2) then
   Pf = Pf * float(1-2*mod(N/2-1,2))
   signs=signs*(1-2*mod(N/2-1,2))
endif
end if !    N==2b
end if !    mod(N,2)==1
return
end subroutine ZPfaffianH
!=======================================================================
!MODIFY WF SO THAT IT CAN BE USED IN OVERLAP CALCULATION
!=======================================================================
 SUBROUTINE WADD_AZ(IRIGHT,ICHARG,ISIMPY,LDBASE,LTIMUP,LUPPER,NOTOCC,PHASPI,WALEFT)
 IMPLICIT NONE
!=======================================================================
!PUT ORIGIN WF TO FLIPPED COLOMN AND ADD FLIPPED COLOMN TO THE END
!NEED TO CALL RETU_AZ BEFORE THIS TO REMOVE POSSIBLE 0 COLOMN
!RETURN MODIFIED WALEFT AND UPDATE LUPPER
!=======================================================================
 INTEGER::IRIGHT,ICHARG,ISIMPY,LDBASE,LTIMUP,LUPPER
 REAL   ::NOTOCC
 COMPLEX,DIMENSION(1:NDBASE,0:NDREVE,0:NDSPIN)::PHASPI
 COMPLEX,DIMENSION(1:NDBASE,1:4*NDSTAT,0:NDSPIN)::WALEFT!,WSINGL
 INTEGER::BR,BSTATE,ISTATE,JSTATE,IBASE,ISPIN,WB,WA
 REAL   ::VSQ
 COMPLEX::C_ZERO
 COMPLEX,ALLOCATABLE::WLAUX(:,:)
!
 CALL INIADD
 C_ZERO=CMPLX(0.0D0,0.0D0)
!
 BR=BLOLAB(ICHARG)%IR
 WB=WLAB_B(ICHARG)
 WA=WLAB_A(ICHARG)
!
 IF(V2BLO(ICHARG).GT.NOTOCC)THEN !IF FLIPPED STATE IS NOT STORED, DO NOTHING
    IF(IRIGHT.EQ.0)THEN
    DO ISPIN=0,NDSPIN
       IF(ISIMPY.EQ.0)THEN
          WBAZER(1:LDBASE,1,ISPIN,ICHARG)&
        =(B_ZERO(1:LDBASE,0,ICHARG)*PHASPI(1:LDBASE,0,ISPIN)&
         +B_ZERO(1:LDBASE,1,ICHARG)*PHASPI(1:LDBASE,1,ISPIN))
!
          WBAZER(1:LDBASE,2,ISPIN,ICHARG)&
        =(A_ZERO(1:LDBASE,0,ICHARG)*(2*ISPIN-1)*CONJG(PHASPI(1:LDBASE,0,1-ISPIN))&
         +A_ZERO(1:LDBASE,1,ICHARG)*(2*ISPIN-1)*CONJG(PHASPI(1:LDBASE,1,1-ISPIN)))
       ELSE
!
          WBAZER(1:LDBASE,1,ISPIN,ICHARG)&
         =B_ZERO(1:LDBASE,BR,ICHARG)*PHASPI(1:LDBASE,BR,ISPIN)
!
          WBAZER(1:LDBASE,2,ISPIN,ICHARG)&
         =A_ZERO(1:LDBASE,1-BR,ICHARG)*(2*ISPIN-1)*CONJG(PHASPI(1:LDBASE,1-BR,1-ISPIN))
       ENDIF
    ENDDO
!
    WSINGL(1:LDBASE,1:2,:)=WBAZER(1:LDBASE,1:2,:,ICHARG)
    ENDIF

    ALLOCATE(WLAUX(1:2*LDBASE,0:NDSPIN))!add the flipped state to the end
    WLAUX(1       :  LDBASE,:)=WALEFT(1:LDBASE,WB,:)
    WLAUX(1+LDBASE:2*LDBASE,:)=WALEFT(1:LDBASE,WA,:)
!
    WALEFT(1:LDBASE,WB,0:NDSPIN)=WSINGL(1:LDBASE,1,:)
    WALEFT(1:LDBASE,WA,0:NDSPIN)=WSINGL(1:LDBASE,2,:)
!
    WALEFT(1:LDBASE,2+LUPPER:2*LUPPER+1,:)=WALEFT(1:LDBASE,LUPPER+1:2*LUPPER,:)
!
    WALEFT(1:LDBASE,  LUPPER+1,:)=WLAUX(1       :  LDBASE,:)
    WALEFT(1:LDBASE,2*LUPPER+2,:)=WLAUX(1+LDBASE:2*LDBASE,:)
!
    LUPPER=LUPPER+1
    DEALLOCATE(WLAUX)
 ENDIF
!
 RETURN
 END SUBROUTINE WADD_AZ
!=======================================================================
 SUBROUTINE RETU_AZ(ICHARG,ISIMPY,LDBASE,LTIMUP,LUPPER,WLAB_BOLD,NOTOCC,LTINEW,LUPNEW,WLAB_BNEW,WALEFT)
 IMPLICIT NONE
!=======================================================================
! WHEN ISIMPY.EQ.1, PUT FILPPED COLOMN BACK TO 0 COLOMN IF FILPPED STATE IS OCCUPIED
! REMOVE ALL REMAIN 0 COLUMN IN BOTH ISIMPY CASE
! CALL THIS BEFORE CALL WADD_AZ
!=======================================================================
 INTEGER::ICHARG,ISIMPY,LDBASE,LTIMUP,LUPPER,WLAB_BOLD,LTINEW,LUPNEW,WLAB_BNEW
 REAL   ::NOTOCC
 COMPLEX,DIMENSION(1:NDBASE,1:4*NDSTAT,0:NDSPIN)::WALEFT
!
 INTEGER::BR,WB,WA,ADDSTA
 COMPLEX,ALLOCATABLE::WAAUX(:,:,:)
 BR=BLOLAB(ICHARG)%IR
 WB=WLAB_B(ICHARG)!LABEL OF BLOCKED STATE
 ADDSTA=LUPPER-BR*(LUPPER-LTIMUP)! LABEL OF ADDED CLOLOMN
!
 LUPNEW=LUPPER!LUPPER,LTIMUP,WLAB_B MAY BE CHANGED WHEN REMOVE THE ADDED TERM
 LTINEW=LTIMUP
 WLAB_BNEW=WLAB_BOLD
!
 IF(V2BLO(ICHARG).GT.NOTOCC .AND. ISIMPY.EQ.1)THEN!IF FLIPPED STATE IS NOT STORED, DO NOTHING
    ALLOCATE(WAAUX(1:2*LDBASE,1:LUPPER,0:NDSPIN))
    WAAUX(1       :  LDBASE,1:LUPPER,:)=WALEFT(1:LDBASE,1       :  LUPPER,:)
    WAAUX(1+LDBASE:2*LDBASE,1:LUPPER,:)=WALEFT(1:LDBASE,1+LUPPER:2*LUPPER,:)
!   PUT ADDED COLOMN TO ZERO COLOMN
    WAAUX(1:2*LDBASE,WB,:)=WAAUX(1:2*LDBASE,ADDSTA,:)
!   REMOVE ADDED COLOMN
    IF(BR.EQ.1)THEN!IF THE ADDED COLOMN IS AT LTIMUP, LTIMUP-1, WLAB_B-1
    WAAUX(1:2*LDBASE,LTIMUP:LUPPER-1,:)=WAAUX(1:2*LDBASE,LTIMUP+1:LUPPER,:)
    LTINEW=LTIMUP-1
    WLAB_BNEW=WLAB_BOLD-1
    ENDIF
!
    LUPNEW=LUPPER-1! LUPPER WILL BECOME LUPPER-1
!    WAAUX(1:2*LDBASE,LUPPER,:)=CMPLX(0.0D0,0.0D0)
!
    WALEFT(1:LDBASE,1       :  LUPNEW,:)=WAAUX(1       :  LDBASE,1:LUPNEW,:)
    WALEFT(1:LDBASE,1+LUPNEW:2*LUPNEW,:)=WAAUX(1+LDBASE:2*LDBASE,1:LUPNEW,:)
    WALEFT(1:LDBASE,2*LUPNEW+1:2*LUPNEW+2,:)=CMPLX(0.0D0,0.0D0)
    DEALLOCATE(WAAUX)
 ENDIF
!REMOVE ALL REMAIN 0 COLOMN, THIS IS NEED WHEN PAIRING DISAPEARED
 CALL RMEMPT(LDBASE,LTINEW,LUPNEW,WALEFT)
 RETURN
 END SUBROUTINE RETU_AZ
!=======================================================================
 SUBROUTINE RMEMPT(LDBASE,LTIMUP,LUPPER,WALEFT)
 IMPLICIT NONE
 INTEGER::LDBASE,LTIMUP,LUPPER
 COMPLEX,DIMENSION(1:NDBASE,1:4*NDSTAT,0:NDSPIN)::WALEFT
!
 INTEGER::ISTATE,JSTATE,KSTATE
 REAL   ::VSQ
 COMPLEX::C_ZERO
 COMPLEX,ALLOCATABLE::WFAUX(:,:,:)
!REMOVE ALL EMPTY STATE
 ALLOCATE(WFAUX(1:2*LDBASE,1:LUPPER,0:NDSPIN))
 C_ZERO=CMPLX(0.0D0,0.0D0)
 JSTATE=0
 KSTATE=0
 DO ISTATE=1,LUPPER
    VSQ=SUM(ABS(WALEFT(1:LDBASE,ISTATE,0))**2)+SUM(ABS(WALEFT(1:LDBASE,ISTATE,1))**2)
    IF(VSQ.GT.1.0D-8)THEN
       JSTATE=JSTATE+1
       IF(ISTATE.LE.LTIMUP)KSTATE=KSTATE+1
       WFAUX(1       :  LDBASE,JSTATE,:)=WALEFT(1:LDBASE,ISTATE,:)
       WFAUX(1+LDBASE:2*LDBASE,JSTATE,:)=WALEFT(1:LDBASE,ISTATE+LUPPER,:)
    ELSE
    WRITE(*,*)'ISTATE REMOVED',ISTATE,VSQ
    ENDIF
 ENDDO
!
 LTIMUP=KSTATE
 LUPPER=JSTATE
!
 WALEFT=C_ZERO
 WALEFT(1:LDBASE,1       :  LUPPER,:)=WFAUX(1       :  LDBASE,1:LUPPER,:)
 WALEFT(1:LDBASE,1+LUPPER:2*LUPPER,:)=WFAUX(1+LDBASE:2*LDBASE,1:LUPPER,:)
 DEALLOCATE(WFAUX)
 RETURN
 END SUBROUTINE RMEMPT
!=======================================================================
 SUBROUTINE ALLPFF(IALLO,IBA,ISG,IWL,IWR,IVU,ISLT,ISLO,ISRO)
 IMPLICIT NONE
 INTEGER::IALLO,IBA,ISG,IWL,IWR,IVR,IUR,IVU,ISLT,ISLO,ISRO
 IF(IALLO.EQ.1)THEN!ALLOCATE SPACE
    IF(.NOT. ALLOCATED(WBAZER) .AND. IBA.GE.1)THEN
       ALLOCATE(WBAZER(1:NDBASE,1:2,0:NDSPIN,0:NDISOS))
    ENDIF
    IF(.NOT. ALLOCATED(WSINGL) .AND. ISG.GE.1)THEN
       ALLOCATE(WSINGL(1:NDBASE,1:4*NDSTAT,0:NDSPIN))
    ENDIF
    IF(.NOT. ALLOCATED(WALTMP) .AND. IWL.GE.1)THEN
       ALLOCATE(WALTMP(1:NDBASE,1:4*NDSTAT,0:NDSPIN))
    ENDIF
    IF(.NOT. ALLOCATED(WARTMP) .AND. IWR.GE.1)THEN
       ALLOCATE(WARTMP(1:NDBASE,1:4*NDSTAT,0:NDSPIN))
    ENDIF
    IF(.NOT. ALLOCATED(VUAUX) .AND. IVU.GE.1)THEN
       ALLOCATE (VUAUX(1:2*NDSTAT,1:2*NDSTAT,0:NDISOS))
    ENDIF
    IF(.NOT. ALLOCATED(SALTMP) .AND. ISLT.GE.1)THEN
       ALLOCATE (SALTMP(1:NDBASE,1:4*NDSTAT,0:NDSPIN,0:NDISOS))
    ENDIF
    IF(.NOT. ALLOCATED(SALODD) .AND. ISLO.GE.1)THEN
       ALLOCATE (SALODD(1:NDBASE,1:4*NDSTAT,0:NDSPIN,0:NDISOS))
    ENDIF
    IF(.NOT. ALLOCATED(SARODD) .AND. ISRO.GE.1)THEN
       ALLOCATE (SARODD(1:NDBASE,1:4*NDSTAT,0:NDSPIN,0:NDISOS))
    ENDIF

 ELSE!DEALLOCATE SPACE
     IF(IBA .GE.1) DEALLOCATE(WBAZER)
     IF(ISG .GE.1) DEALLOCATE(WSINGL)
     IF(IWL .GE.1) DEALLOCATE(WALTMP)
     IF(IWR .GE.1) DEALLOCATE(WARTMP)
     IF(IVU .GE.1) DEALLOCATE(VUAUX)
     IF(ISLT.GE.1) DEALLOCATE(SALTMP)
     IF(ISLO.GE.1) DEALLOCATE(SALODD)
     IF(ISRO.GE.1) DEALLOCATE(SARODD)
 ENDIF
 RETURN
 END SUBROUTINE ALLPFF
!=======================================================================
END MODULE PFAFFI
!=======================================================================
SUBROUTINE INIADD
USE BLOSAV
USE hfodd_sizes
 IMPLICIT NONE
!=======================================================================
!ALLOCATE SPACE TO ARRAY
!=======================================================================
 INTEGER::IALLOC
 COMPLEX::C_ZERO

 IALLOC=0
 C_ZERO=CMPLX(0.0D0,0.0D0)
!
! IF(.NOT.ALLOCATED(WBZERO)) THEN
!    ALLOCATE(WBZERO(1:NDBASE,0:NDREVE,0:NDSPIN))
!    IF (IALLOC.NE.0) CALL NOALLO('WBZERO','INIADD')
!    WBZERO(:,:,:)=C_ZERO
! ENDIF
! IF(.NOT.ALLOCATED(WAZERO)) THEN
!    ALLOCATE(WAZERO(1:NDBASE,0:NDREVE,0:NDSPIN))
!    IF (IALLOC.NE.0) CALL NOALLO('WAZERO','INIADD')
!    WAZERO(:,:,:)=C_ZERO
! ENDIF
!
 IF (.NOT.ALLOCATED(A_ZERO)) THEN
    ALLOCATE (A_ZERO(1:NDBASE,0:NDREVE,0:NDISOS),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('A_ZERO','INIADD')
    A_ZERO(:,:,:)=C_ZERO
 END IF
 IF (.NOT.ALLOCATED(B_ZERO)) THEN
    ALLOCATE (B_ZERO(1:NDBASE,0:NDREVE,0:NDISOS),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('B_ZERO','INIADD')
    B_ZERO(:,:,:)=C_ZERO
 END IF
!
 IF (.NOT.ALLOCATED(BLOLAB)) THEN
    ALLOCATE (BLOLAB(0:NDISOS),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('BLOLAB','INIADD')
    BLOLAB(:)%IR=0
    BLOLAB(:)%N=0
 END IF
!
 IF (.NOT.ALLOCATED(WLAB_B)) THEN
    ALLOCATE (WLAB_B(0:NDISOS),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('WLAB_B','INIADD')
    WLAB_B(:)=0
 END IF
!
 IF (.NOT.ALLOCATED(WLAB_A)) THEN
    ALLOCATE (WLAB_A(0:NDISOS),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('WLAB_A','INIADD')
    WLAB_A(:)=0
 END IF
!
 IF (.NOT.ALLOCATED(V2BLO)) THEN
    ALLOCATE (V2BLO(0:NDISOS),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('V2BLO','INIADD')
    V2BLO(:)=0.0D0
 END IF
!
 RETURN
END SUBROUTINE INIADD
!=======================================================================
SUBROUTINE TESTOVP(IBLOCK,ICHARG,MREVER,ISIMPY,LDBASE,LTIMUP,LUPPER,PHASPI,WALEFT,WARIGH,OVKERN)
USE PFAFFI
!USE SAVRIG
USE BLOSAV
USE hfodd_sizes
IMPLICIT NONE
 INTEGER::IBLOCK,ICHARG,MREVER,ISIMPY,LDBASE,LTIMUP,LUPPER,&
          LTINEW,LUPNEW,RTINEW,RUPNEW,LLAB_BNEW,RLAB_BNEW,&
          LLAB_BOLD,RLAB_BOLD
 COMPLEX::OVKERN
 COMPLEX,DIMENSION(1:NDBASE,1:4*NDSTAT,0:NDSPIN)::WALEFT,WARIGH
 REAL   ::NOTOCC
 COMPLEX,DIMENSION(1:NDBASE,0:NDREVE,0:NDSPIN)::PHASPI
!
 IF(.NOT.ALLOCATED(WALTMP))ALLOCATE(WALTMP(1:NDBASE,1:4*NDSTAT,0:NDSPIN))
 IF(.NOT.ALLOCATED(WARTMP))ALLOCATE(WARTMP(1:NDBASE,1:4*NDSTAT,0:NDSPIN))
 IF(.NOT.ALLOCATED(Wsingl))ALLOCATE(Wsingl(1:NDBASE,1:4*NDSTAT,0:NDSPIN))
 IF(.NOT.ALLOCATED(WBAZER))ALLOCATE(WBAZER(1:NDBASE,1:2,0:NDSPIN,0:NDISOS))
 WALTMP(:,:,:)=WALEFT(:,:,:)!SARIGH(:,:,:,ICHARG)
 WARTMP(:,:,:)=WARIGH(:,:,:)

!
 NOTOCC=1.0d-8
 LTINEW=LTIMUP
 LUPNEW=LUPPER
 RTINEW=LTIMUP
 RUPNEW=LUPPER
!
 IF(IBLOCK.EQ.1)THEN
 LLAB_BOLD=WLAB_B(ICHARG)
 RLAB_BOLD=WLAB_B(ICHARG)
!
 CALL RETU_AZ(ICHARG,ISIMPY,LDBASE,LTIMUP,LUPPER,LLAB_BOLD,NOTOCC,LTINEW,LUPNEW,LLAB_BNEW,WALTMP)
 CALL RETU_AZ(ICHARG,ISIMPY,LDBASE,LTIMUP,LUPPER,RLAB_BOLD,NOTOCC,RTINEW,RUPNEW,RLAB_BNEW,WARTMP)
!
 WLAB_B(ICHARG)=LLAB_BNEW
 WLAB_A(ICHARG)=LLAB_BNEW+LUPNEW

 CALL WADD_AZ(0,ICHARG,ISIMPY,LDBASE,LTINEW,LUPNEW,NOTOCC,PHASPI,WALTMP)
!
 WSINGL(1:LDBASE,1:2,:)=WBAZER(1:LDBASE,1:2,:,ICHARG)
 CALL ROTWAV(2,WSINGL)
!
 CALL WADD_AZ(1,ICHARG,ISIMPY,LDBASE,RTINEW,RUPNEW,NOTOCC,PHASPI,WARTMP)
 ENDIF

 IF(MREVER.EQ.0) STOP 'PFAFFIAN NEED MREVER.EQ.1'
 CALL NORM_N(0,1,LUPNEW,LDBASE,WALTMP,WALTMP,VUAUX(:,:,ICHARG),LV2AUX(ICHARG),PV2AUX(ICHARG),OVKERN)
 CALL NORM_N(1,1,LUPNEW,LDBASE,WALTMP,WARTMP,VUAUX(:,:,ICHARG),LV2AUX(ICHARG),PV2AUX(ICHARG),OVKERN)
 WRITE(*,*)'OVRLAP',OVKERN,ABS(OVKERN),ICHARG
! IF(ICHARG.EQ.1) STOP 'TEST FINISHIED'
RETURN
END SUBROUTINE TESTOVP
!=======================================================================
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!=======================================================================
MODULE O2AVRG
!routines to calculate <o^2>, where o is a one-body operator
CONTAINS
!=======================================================================
SUBROUTINE AVRXSD(HAUXIL,AOBSLI,AKBSLI,RHO_PP)
!=======================================================================
!Calculate Tr[x.rho], x=[x(0),0;0,x(1)]
!=======================================================================
USE hfodd_sizes
!=======================================================================
IMPLICIT NONE
!=======================================================================
 REAL   ::AOBSLI
 COMPLEX::AKBSLI
 COMPLEX,DIMENSION(1:NDBASE,1:NDBASE,0:NDREVE)::HAUXIL,RHO_PP
!
! INTEGER::INCX,INCY
! COMPLEX::RES1,RES2,ZDOTU
 INTEGER::IBRA,IKET
!
 INTEGER::LDBASE
 COMPLEX::RESULT,C_ZERO
!=======================================================================
 COMMON /DIMENS/ LDBASE
!=======================================================================
 C_ZERO=CMPLX(0.0D0,0.0D0)
 RESULT=C_ZERO
!
 DO IBRA=1,LDBASE
    DO IKET=1,LDBASE
       RESULT=RESULT+HAUXIL(IBRA,IKET,0)*RHO_PP(IKET,IBRA,0)&
                    +HAUXIL(IBRA,IKET,1)*RHO_PP(IKET,IBRA,1)

    END DO
 END DO
 AOBSLI=REAL(RESULT)
 AKBSLI=RESULT
 !INCX=1; INCY=1
 !RES1=ZDOTU(NDBASE*NDBASE,HAUXIL(:,:,0),INCX,TRANSPOSE(RHO_PP(:,:,0)),INCY)
 !RES2=ZDOTU(NDBASE*NDBASE,HAUXIL(:,:,1),INCX,TRANSPOSE(RHO_PP(:,:,0)),INCY)
 !AOBSLI=REAL(RES1+RES2)
 !AKBSLI=RES1+RES2
!=======================================================================
 RETURN
END SUBROUTINE AVRXSD
!=======================================================================
SUBROUTINE AVRXSO(HAUXIL,AOBSLI,AKBSLI,RHO_PM)
!=======================================================================
! Calculate Tr[x.rho], x=[0,x(0);x(1),0]
!=======================================================================
USE hfodd_sizes
!=======================================================================
IMPLICIT NONE
!=======================================================================
 REAL   ::AOBSLI
 COMPLEX::AKBSLI
 COMPLEX,DIMENSION(1:NDBASE,1:NDBASE,0:NDREVE)::HAUXIL,RHO_PM
!
! INTEGER::INCX,INCY
! COMPLEX::RES1,RES2,ZDOTU
 INTEGER::IBRA,IKET
 COMPLEX::RESULT,C_ZERO
!
 INTEGER::LDBASE
 COMMON /DIMENS/ LDBASE
!=======================================================================
 C_ZERO=CMPLX(0.0D0,0.0D0)
 RESULT=C_ZERO
!
 DO IBRA=1,LDBASE
    DO IKET=1,LDBASE
       RESULT=RESULT+HAUXIL(IBRA,IKET,0)*RHO_PM(IKET,IBRA,1)&
                    +HAUXIL(IBRA,IKET,1)*RHO_PM(IKET,IBRA,0)
    END DO
 END DO
 AOBSLI=REAL(RESULT)
 AKBSLI=RESULT
 !INCX=1; INCY=1
 !RES1=ZDOTU(NDBASE*NDBASE,HAUXIL(:,:,0),INCX,TRANSPOSE(RHO_PM(:,:,1)),INCY)
 !RES2=ZDOTU(NDBASE*NDBASE,HAUXIL(:,:,1),INCX,TRANSPOSE(RHO_PM(:,:,0)),INCY)
 !AOBSLI=REAL(RES1+RES2)
 !AKBSLI=RES1+RES2
!=======================================================================
 RETURN
END SUBROUTINE AVRXSO
!=======================================================================
SUBROUTINE EXCXSD(IMSIGN,MREVER,ISIMPY,JETACM,EOBSSQ,EKBSSQ,XMAT,RPP,RPM,APP,APM)
!=======================================================================
!Tr[(x.rho)^2], x=[x(0),0;0,x(1)]
!=======================================================================
USE hfodd_sizes
!=======================================================================
IMPLICIT NONE
!=======================================================================
 INTEGER::IMSIGN,MREVER,ISIMPY,JETACM
 REAL   ::EOBSSQ
 COMPLEX::EKBSSQ
 COMPLEX,DIMENSION(1:NDBASE,1:NDBASE,0:NDREVE)::XMAT,RPP
 COMPLEX,DIMENSION(1:NDBASE,1:NDBASE,0:NDREVE),OPTIONAL::RPM,APP,APM
!
 INTEGER::IREVER,INCX,INCY
 COMPLEX::C_ZERO,RESULT,C_UNIT,CTWO,RES1,RES2,ZDOTU
 COMPLEX,DIMENSION(1:NDBASE,1:NDBASE,0:NDREVE)::AUXILI
!
 INTEGER::LDBASE
 COMMON /DIMENS/ LDBASE
!=======================================================================
 C_ZERO=CMPLX(0.0D0,0.0D0)
 IF(JETACM.EQ.2) THEN
                    APP(1:NDBASE,1:NDBASE,0:NDREVE)=C_ZERO
    IF(ISIMPY.EQ.0) APM(1:NDBASE,1:NDBASE,0:NDREVE)=C_ZERO
 END IF
!=======================================================================
!X.RHO
!=======================================================================
 C_UNIT=CMPLX(1.0D0,0.0D0)
 DO IREVER=0,MREVER
    AUXILI(:,:,IREVER) = C_ZERO
    CALL ZGEMM('N','N',LDBASE,LDBASE,LDBASE, &
                       C_UNIT,XMAT(:,:,IREVER),NDBASE, &
                       RPP(:,:,IREVER),NDBASE,C_ZERO,AUXILI(:,:,IREVER),NDBASE)
 END DO
!=======================================================================
!TR[(X.RHO)^2] AND 2(X.RHO.X)_{IJ}
!======================================================================
 RESULT=C_ZERO
 INCX=1; INCY=1
 DO IREVER=0,MREVER
    RESULT=RESULT+ZDOTU(NDBASE*NDBASE,AUXILI(:,:,IREVER),INCX,TRANSPOSE(AUXILI(:,:,IREVER)),INCY)
    IF (JETACM.EQ.2) THEN
        CTWO=CMPLX(-2.0D0,0.0D0)
        CALL ZGEMM('N','N',LDBASE,LDBASE,LDBASE, &
                           CTWO,AUXILI(:,:,IREVER),NDBASE, &
                           XMAT(:,:,IREVER),NDBASE,C_ZERO,APP(:,:,IREVER),NDBASE)
    END IF
 END DO
 !
 IF (MREVER.EQ.0) THEN
    RESULT=RESULT+CONJG(RESULT)*(IMSIGN**2)
    IF(JETACM.EQ.2) APP(:,:,1)=CONJG(APP(:,:,0))*(IMSIGN**2)
 END IF
 !
 EOBSSQ = REAL(RESULT)
 EKBSSQ = RESULT
 !
!=======================================================================
!ADDITIONAL TERMS WHEN SIMPLEX IS NOT CONSERVED
!=======================================================================
 IF (ISIMPY.NE.1) THEN
!=======================================================================
!X.RHO
!=======================================================================
    C_UNIT=CMPLX(1.0D0,0.0D0); AUXILI = C_ZERO
    CALL ZGEMM('N','N',LDBASE,LDBASE,LDBASE, &
                       C_UNIT,XMAT(:,:,0),NDBASE, &
                       RPM(:,:,0),NDBASE,C_ZERO,AUXILI(:,:,0),NDBASE)
    CALL ZGEMM('N','N',LDBASE,LDBASE,LDBASE, &
                       C_UNIT,XMAT(:,:,1),NDBASE, &
                       RPM(:,:,1),NDBASE,C_ZERO,AUXILI(:,:,1),NDBASE)
!=======================================================================
!TR[(X.RHO)^2] AND (X.RHO.X)_{IJ}
!=======================================================================
    INCX=1; INCY=1
    RESULT=ZDOTU(NDBASE*NDBASE,AUXILI(:,:,0),INCX,TRANSPOSE(AUXILI(:,:,1)),INCY)
    IF (JETACM.EQ.2) THEN
        IF (MREVER.EQ.0) STOP ' WRONG MREVER IN EXCXSD'
        CTWO=CMPLX(-2.0D0,0.0D0)
        CALL ZGEMM('N','N',LDBASE,LDBASE,LDBASE, &
                           CTWO,AUXILI(:,:,0),NDBASE, &
                           XMAT(:,:,1),NDBASE,C_ZERO,APM(:,:,0),NDBASE)
        CALL ZGEMM('N','N',LDBASE,LDBASE,LDBASE, &
                           CTWO,AUXILI(:,:,1),NDBASE, &
                           XMAT(:,:,0),NDBASE,C_ZERO,APM(:,:,1),NDBASE)
    END IF
    !
    RESULT=2*RESULT
    !
    EOBSSQ = EOBSSQ + REAL(RESULT)
    EKBSSQ = EKBSSQ + RESULT
    !
 END IF
!======================================================================
 RETURN
 END SUBROUTINE EXCXSD
!=======================================================================
SUBROUTINE EXCXSO(IMSIGN,MREVER,ISIMPY,JETACM,EOBSSQ,EKBSSQ,XMAT,RPP,RPM,APP,APM)
!=======================================================================
! Tr[(x.rho)^2], x=[0,x(0);x(1),0]
!=======================================================================
USE hfodd_sizes
!=======================================================================
IMPLICIT NONE
!=======================================================================
 INTEGER::IMSIGN,MREVER,ISIMPY,JETACM
 REAL   ::EOBSSQ
 COMPLEX::EKBSSQ
 COMPLEX::C_ZERO,RESULT,C_UNIT,CTWO,RES1,RES2,ZDOTU
 COMPLEX,DIMENSION(1:NDBASE,1:NDBASE,0:NDREVE)::XMAT,RPP
 COMPLEX,DIMENSION(1:NDBASE,1:NDBASE,0:NDREVE),OPTIONAL::RPM,APP,APM
!
 INTEGER::IREVER,INCX,INCY
 COMPLEX,DIMENSION(1:NDBASE,1:NDBASE,0:NDREVE)::AUXILI
!
 INTEGER::LDBASE
 COMMON /DIMENS/ LDBASE
!=======================================================================
 C_ZERO=CMPLX(0.0D0,0.0D0)
  IF(JETACM.EQ.2)THEN
                    APP(1:NDBASE,1:NDBASE,0:NDREVE)=C_ZERO
    IF(ISIMPY.EQ.0) APM(1:NDBASE,1:NDBASE,0:NDREVE)=C_ZERO
 END IF
!=======================================================================
!X.RHO
!=======================================================================
 C_UNIT=CMPLX(1.0D0,0.0D0); AUXILI = C_ZERO
 DO IREVER=0,MREVER
    CALL ZGEMM('N','N',LDBASE,LDBASE,LDBASE, &
                       C_UNIT,XMAT(:,:,1-IREVER),NDBASE, &
                       RPP(:,:,IREVER),NDBASE,C_ZERO,AUXILI(:,:,IREVER),NDBASE)
 END DO
 IF (MREVER.EQ.0) THEN
    AUXILI(:,:,1)=CONJG(AUXILI(:,:,0))*IMSIGN
 END IF
!=======================================================================
!TR[(X.RHO)^2] AND 2(X.RHO.X)
!=======================================================================
 INCX=1; INCY=1
 RESULT=ZDOTU(NDBASE*NDBASE,AUXILI(:,:,1),INCX,TRANSPOSE(AUXILI(:,:,0)),INCY)
 IF (JETACM.EQ.2) THEN
     CTWO=CMPLX(-2.0D0,0.0D0)
     CALL ZGEMM('N','N',LDBASE,LDBASE,LDBASE, &
                        CTWO,AUXILI(:,:,1),NDBASE, &
                        XMAT(:,:,1),NDBASE,C_ZERO,APP(:,:,0),NDBASE)
     CALL ZGEMM('N','N',LDBASE,LDBASE,LDBASE, &
                        CTWO,AUXILI(:,:,0),NDBASE, &
                        XMAT(:,:,0),NDBASE,C_ZERO,APP(:,:,1),NDBASE)
 END IF
 !
 EOBSSQ = 2*REAL(RESULT)
 EKBSSQ = 2*RESULT
 !
!=======================================================================
!ADDITIONAL TERMS WHEN SIMPLEX IS NOT CONSERVED
!=======================================================================
 IF (ISIMPY.NE.1) THEN
!=======================================================================
!X.RHO
!=======================================================================
     C_UNIT=CMPLX(1.0D0,0.0D0); AUXILI = C_ZERO
     CALL ZGEMM('N','N',LDBASE,LDBASE,LDBASE, &
                        C_UNIT,XMAT(:,:,1),NDBASE, &
                        RPM(:,:,0),NDBASE,C_ZERO,AUXILI(:,:,0),NDBASE)
     CALL ZGEMM('N','N',LDBASE,LDBASE,LDBASE, &
                        C_UNIT,XMAT(:,:,0),NDBASE, &
                        RPM(:,:,1),NDBASE,C_ZERO,AUXILI(:,:,1),NDBASE)
!=======================================================================
!TR[(X.RHO)^2] AND (X.RHO.X)_{IJ}
!=======================================================================
     INCX=1; INCY=1
     RES1=ZDOTU(NDBASE*NDBASE,AUXILI(:,:,0),INCX,TRANSPOSE(AUXILI(:,:,0)),INCY)
     RES2=ZDOTU(NDBASE*NDBASE,AUXILI(:,:,1),INCX,TRANSPOSE(AUXILI(:,:,1)),INCY)
     RESULT=RES1+RES2
     IF (JETACM.EQ.2) THEN
         CTWO=CMPLX(-2.0D0,0.0D0)
         CALL ZGEMM('N','N',LDBASE,LDBASE,LDBASE, &
                            CTWO,AUXILI(:,:,1),NDBASE, &
                            XMAT(:,:,0),NDBASE,C_ZERO,APM(:,:,0),NDBASE)
         CALL ZGEMM('N','N',LDBASE,LDBASE,LDBASE, &
                            CTWO,AUXILI(:,:,0),NDBASE, &
                            XMAT(:,:,1),NDBASE,C_ZERO,APM(:,:,1),NDBASE)
     END IF
     !
     EOBSSQ = EOBSSQ + REAL(RESULT)
     EKBSSQ = EKBSSQ + RESULT
     !
 END IF
!=======================================================================
 RETURN
END SUBROUTINE EXCXSO
!=======================================================================
SUBROUTINE PAIXSD(MREVER,ISIMPY,JETACM,POBSSQ,PKBSSQ,XMAT,KPM,KPP1,KPP2,PPM,PPP)
!=======================================================================
! Tr[(X.kappa)(X.kappa')*] x=[x(0),0;0,x(1)]
!=======================================================================
USE hfodd_sizes
!=======================================================================
IMPLICIT NONE
!=======================================================================
 INTEGER::MREVER,ISIMPY,JETACM
 REAL   ::POBSSQ
 COMPLEX::PKBSSQ
 COMPLEX,DIMENSION(1:NDBASE,1:NDBASE,0:NDREVE)::XMAT,KPM
 COMPLEX,DIMENSION(1:NDBASE,1:NDBASE,0:NDREVE),OPTIONAL::KPP1,KPP2,PPM,PPP
!
 INTEGER::IREVER,INCX,INCY
 COMPLEX::C_ZERO,RESULT,C_UNIT,CTWO,RES1,RES2,ZDOTU
 COMPLEX,DIMENSION(1:NDBASE,1:NDBASE,0:NDREVE)::AUXILI
!
 INTEGER::LDBASE
 COMMON /DIMENS/ LDBASE
!=======================================================================
 C_ZERO=CMPLX(0.0D0,0.0D0)
 IF(JETACM.EQ.2)THEN
                    PPM(1:NDBASE,1:NDBASE,0:NDREVE)=C_ZERO
    IF(ISIMPY.EQ.0) PPP(1:NDBASE,1:NDBASE,0:NDREVE)=C_ZERO
 END IF
!=======================================================================
!X.KAPPA AND X.KAPPA'
!=======================================================================
 C_UNIT=CMPLX(1.0D0,0.0D0); AUXILI = C_ZERO
 CALL ZGEMM('N','N',LDBASE,LDBASE,LDBASE, &
                    C_UNIT,XMAT(:,:,0),NDBASE, &
                    KPM(:,:,0),NDBASE,C_ZERO,AUXILI(:,:,0),NDBASE)
 C_UNIT=CMPLX(-1.0D0,0.0D0)
 CALL ZGEMM('N','T',LDBASE,LDBASE,LDBASE, &
                    C_UNIT,XMAT(:,:,1),NDBASE, &
                    KPM(:,:,0),NDBASE,C_ZERO,AUXILI(:,:,1),NDBASE)
!=======================================================================
!TR[(X.KAPPA)(X.KAPPA)*] AND 2[X.KAPPA.X']
!=======================================================================
 INCX=1; INCY=1
 RES1=ZDOTU(NDBASE*NDBASE,AUXILI(:,:,0),INCX,TRANSPOSE(CONJG(AUXILI(:,:,1))),INCY)
 RES2=ZDOTU(NDBASE*NDBASE,AUXILI(:,:,1),INCX,TRANSPOSE(CONJG(AUXILI(:,:,0))),INCY)
 RESULT=RES1+RES2
 IF (JETACM.EQ.2) THEN
     CTWO=CMPLX(2.0D0,0.0D0)
     CALL ZGEMM('N','T',LDBASE,LDBASE,LDBASE, &
                        CTWO,AUXILI(:,:,1),NDBASE, &
                        XMAT(:,:,0),NDBASE,C_ZERO,PPM(:,:,1),NDBASE)
     CALL ZGEMM('N','T',LDBASE,LDBASE,LDBASE, &
                        CTWO,AUXILI(:,:,0),NDBASE, &
                        XMAT(:,:,1),NDBASE,C_ZERO,PPM(:,:,0),NDBASE)
 END IF
!=======================================================================
!    attention: before version 11 and 278h, the two IREVER components
!               above where exchanged, which did not conform  to  the
!               standard of defining arrays PLINPM and PROTPM.
!               This bug was corrected on 06/12/2016
!=======================================================================
 !
 POBSSQ = -REAL(RESULT)
 PKBSSQ = -RESULT
 !
!=======================================================================
!ADDITIONAL TERMS WHEN SIMPLEX IS NOT CONSERVED
!=======================================================================
 IF (ISIMPY.NE.1) THEN
!=======================================================================
!X.KAPPA
!=======================================================================
     C_UNIT=CMPLX(1.0D0,0.0D0); AUXILI = C_ZERO
     CALL ZGEMM('N','N',LDBASE,LDBASE,LDBASE, &
                        C_UNIT,XMAT(:,:,0),NDBASE, &
                        KPP1(:,:,0),NDBASE,C_ZERO,AUXILI(:,:,0),NDBASE)
     CALL ZGEMM('N','N',LDBASE,LDBASE,LDBASE, &
                        C_UNIT,XMAT(:,:,1),NDBASE, &
                        KPP1(:,:,1),NDBASE,C_ZERO,AUXILI(:,:,1),NDBASE)
!=======================================================================
!TR[(X.KAPPA)(X.KAPPA)*] AND 2[X.KAPPA.X]
!=======================================================================
     INCX=1; INCY=1
     RES1=ZDOTU(NDBASE*NDBASE,AUXILI(:,:,0),INCX,TRANSPOSE(CONJG(AUXILI(:,:,0))),INCY)
     RES2=ZDOTU(NDBASE*NDBASE,AUXILI(:,:,1),INCX,TRANSPOSE(CONJG(AUXILI(:,:,1))),INCY)
     RESULT=RES1+RES2
     IF (JETACM.EQ.2) THEN
         CTWO=CMPLX(2.0D0,0.0D0)
         CALL ZGEMM('N','T',LDBASE,LDBASE,LDBASE, &
                            CTWO,AUXILI(:,:,0),NDBASE, &
                            XMAT(:,:,0),NDBASE,C_ZERO,PPP(:,:,0),NDBASE)
         CALL ZGEMM('N','T',LDBASE,LDBASE,LDBASE, &
                            CTWO,AUXILI(:,:,1),NDBASE, &
                            XMAT(:,:,1),NDBASE,C_ZERO,PPP(:,:,1),NDBASE)
     END IF
     !
     POBSSQ = POBSSQ - REAL(RESULT)
     PKBSSQ = PKBSSQ - RESULT
     !
 END IF
!=======================================================================
 RETURN
END SUBROUTINE PAIXSD
!=======================================================================
SUBROUTINE PAIXSO(MREVER,ISIMPY,JETACM,POBSSQ,PKBSSQ,XMAT,KPM,KPP1,KPP2,PPM,PPP)
!=======================================================================
! Tr[(X.kappa)(X.kappa')*] x=[0,x(0);x(1),0]
!=======================================================================
USE hfodd_sizes
!=======================================================================
IMPLICIT NONE
!=======================================================================
 INTEGER::MREVER,ISIMPY,JETACM
 REAL   ::POBSSQ
 COMPLEX::PKBSSQ
 COMPLEX,DIMENSION(1:NDBASE,1:NDBASE,0:NDREVE)::XMAT,KPM
 COMPLEX,DIMENSION(1:NDBASE,1:NDBASE,0:NDREVE),OPTIONAL::KPP1,KPP2,PPM,PPP
!
 INTEGER::IREVER,INCX,INCY
 COMPLEX::C_ZERO,RESULT,C_UNIT,CTWO,RES1,RES2,ZDOTU
 COMPLEX,DIMENSION(1:NDBASE,1:NDBASE,0:NDREVE)::AUXILI
!
 INTEGER::LDBASE
 COMMON /DIMENS/ LDBASE
!=======================================================================
 C_ZERO=CMPLX(0.0D0,0.0D0)
 IF(JETACM.EQ.2)THEN
                    PPM(1:NDBASE,1:NDBASE,0:NDREVE)=C_ZERO
    IF(ISIMPY.EQ.0) PPP(1:NDBASE,1:NDBASE,0:NDREVE)=C_ZERO
 END IF
!=======================================================================
!X.KAPPA
!=======================================================================
 C_UNIT=CMPLX(1.0D0,0.0D0); AUXILI = C_ZERO
 CALL ZGEMM('N','N',LDBASE,LDBASE,LDBASE, &
                    C_UNIT,XMAT(:,:,0),NDBASE, &
                    KPM(:,:,1),NDBASE,C_ZERO,AUXILI(:,:,0),NDBASE)
 C_UNIT=CMPLX(-1.0D0,0.0D0)
 CALL ZGEMM('N','T',LDBASE,LDBASE,LDBASE, &
                    C_UNIT,XMAT(:,:,1),NDBASE, &
                    KPM(:,:,1),NDBASE,C_ZERO,AUXILI(:,:,1),NDBASE)
!=======================================================================
!TR[(X.KAPPA)(X.KAPPA)*] AND 2[X.KAPPA.X']
!=======================================================================
 INCX=1; INCY=1
 RES1=ZDOTU(NDBASE*NDBASE,AUXILI(:,:,0),INCX,TRANSPOSE(CONJG(AUXILI(:,:,0))),INCY)
 RES2=ZDOTU(NDBASE*NDBASE,AUXILI(:,:,1),INCX,TRANSPOSE(CONJG(AUXILI(:,:,1))),INCY)
 RESULT=RES1+RES2
 IF (JETACM.EQ.2) THEN
     CTWO=CMPLX(2.0D0,0.0D0)
     CALL ZGEMM('N','T',LDBASE,LDBASE,LDBASE, &
                        CTWO,AUXILI(:,:,1),NDBASE, &
                        XMAT(:,:,0),NDBASE,C_ZERO,PPM(:,:,1),NDBASE)
     CALL ZGEMM('N','T',LDBASE,LDBASE,LDBASE, &
                        CTWO,AUXILI(:,:,0),NDBASE, &
                        XMAT(:,:,1),NDBASE,C_ZERO,PPM(:,:,0),NDBASE)
 END IF
!=======================================================================
!    attention: before version 11 and 278h, the two IREVER components
!               above where exchanged, which did not conform  to  the
!               standard of defining arrays PLINPM and PROTPM.
!               This bug was corrected on 06/12/2016
!=======================================================================
 !
 POBSSQ = -REAL(RESULT)
 PKBSSQ = -RESULT
 !
!=======================================================================
!ADDITIONAL TERMS WHEN SIMPLEX IS NOT CONSERVED
!=======================================================================
 IF (ISIMPY.NE.1) THEN
!=======================================================================
!X.KAPPA
!=======================================================================
     C_UNIT=CMPLX(1.0D0,0.0D0); AUXILI = C_ZERO
     CALL ZGEMM('N','N',LDBASE,LDBASE,LDBASE, &
                        C_UNIT,XMAT(:,:,1),NDBASE, &
                        KPP1(:,:,0),NDBASE,C_ZERO,AUXILI(:,:,0),NDBASE)
     CALL ZGEMM('N','N',LDBASE,LDBASE,LDBASE, &
                        C_UNIT,XMAT(:,:,0),NDBASE, &
                        KPP1(:,:,1),NDBASE,C_ZERO,AUXILI(:,:,1),NDBASE)
!=======================================================================
!TR[(X.KAPPA)(X.KAPPA')*] AND 2[X.KAPPA.X]
!=======================================================================
     INCX=1; INCY=1
     RES1=ZDOTU(NDBASE*NDBASE,AUXILI(:,:,1),INCX,TRANSPOSE(CONJG(AUXILI(:,:,0))),INCY)
     RES2=ZDOTU(NDBASE*NDBASE,AUXILI(:,:,0),INCX,TRANSPOSE(CONJG(AUXILI(:,:,1))),INCY)
     RESULT=RES1+RES2
     IF (JETACM.EQ.2) THEN
         CTWO=CMPLX(2.0D0,0.0D0)
         CALL ZGEMM('N','T',LDBASE,LDBASE,LDBASE, &
                            CTWO,AUXILI(:,:,1),NDBASE, &
                            XMAT(:,:,0),NDBASE,C_ZERO,PPP(:,:,0),NDBASE)
         CALL ZGEMM('N','T',LDBASE,LDBASE,LDBASE, &
                            CTWO,AUXILI(:,:,0),NDBASE, &
                            XMAT(:,:,1),NDBASE,C_ZERO,PPP(:,:,1),NDBASE)
     END IF
     !
     POBSSQ = POBSSQ - REAL(RESULT)
     PKBSSQ = PKBSSQ - RESULT
     !
 END IF
!=======================================================================
 RETURN
END SUBROUTINE PAIXSO
!=======================================================================
END MODULE O2AVRG
!=======================================================================
MODULE OFFDIA
!routines to calculate transformation matrix, transformed wave function or energy kernal
USE hfodd_sizes
COMPLEX,ALLOCATABLE::VRVAUX(:,:),URUAUX(:,:)
!
CONTAINS
SUBROUTINE WLRMAT(IHFB,LDBASE,LUPPER,WALEFT,WARIGH,OVRV,OVRU)
!CALCULATE OVRLAP MATRIX FROM WL, WR
IMPLICIT NONE
 INTEGER::IHFB,LDBASE,LUPPER
 COMPLEX,DIMENSION(1:NDBASE,1:4*NDSTAT,0:NDSPIN)::WALEFT,WARIGH
 COMPLEX,DIMENSION(1:LUPPER,1:LUPPER)::OVRV
 COMPLEX,DIMENSION(1:LUPPER,1:LUPPER),OPTIONAL::OVRU
!
 COMPLEX::C_ZERO,C_UNIT
!
 C_ZERO=CMPLX(0.0D0,0.0D0)
 C_UNIT=CMPLX(1.0D0,0.0D0)
!
 CALL ZGEMM('C','N',LUPPER,LUPPER,LDBASE,C_UNIT,&
             WALEFT(1:LDBASE,1:LUPPER,0),LDBASE,&
             WARIGH(1:LDBASE,1:LUPPER,0),LDBASE,C_ZERO,OVRV,LUPPER)
!
 CALL ZGEMM('C','N',LUPPER,LUPPER,LDBASE,C_UNIT,&
             WALEFT(1:LDBASE,1:LUPPER,1),LDBASE,&
             WARIGH(1:LDBASE,1:LUPPER,1),LDBASE,C_UNIT,OVRV,LUPPER)
!
 IF(IHFB.GE.1)THEN
    CALL ZGEMM('C','N',LUPPER,LUPPER,LDBASE,C_UNIT,&
                WALEFT(1:LDBASE,1+LUPPER:2*LUPPER,0),LDBASE,&
                WARIGH(1:LDBASE,1+LUPPER:2*LUPPER,0),LDBASE,C_ZERO,OVRU,LUPPER)
!
    CALL ZGEMM('C','N',LUPPER,LUPPER,LDBASE,C_UNIT,&
                WALEFT(1:LDBASE,1+LUPPER:2*LUPPER,1),LDBASE,&
                WARIGH(1:LDBASE,1+LUPPER:2*LUPPER,1),LDBASE,C_UNIT,OVRU,LUPPER)
 ENDIF
!
 RETURN
END SUBROUTINE WLRMAT
!=======================================================================
SUBROUTINE MAT_WR(IHFB,LDBASE,LUPPER,OVRLAP,WARIGH)
IMPLICIT NONE
 INTEGER::IHFB,LDBASE,LUPPER
 COMPLEX,DIMENSION(1:NDBASE,1:4*NDSTAT,0:NDSPIN)::WARIGH
 COMPLEX,DIMENSION(1:LUPPER,1:LUPPER)::OVRLAP
!
 INTEGER::ISPIN
 COMPLEX::C_UNIT,C_ZERO
 COMPLEX,ALLOCATABLE::WARNEW(:,:)
!
 C_ZERO=CMPLX(0.0D0,0.0D0)
 C_UNIT=CMPLX(1.0D0,0.0D0)
 ALLOCATE(WARNEW(1:LDBASE,1:LUPPER))

 DO ISPIN=0,NDSPIN
    CALL ZGEMM('N','N',LDBASE,LUPPER,LUPPER,C_UNIT,&
             WARIGH(1:LDBASE,1:LUPPER,ISPIN),LDBASE,&
             OVRLAP(1:LUPPER,1:LUPPER),      LUPPER,&
             C_ZERO,WARNEW,LDBASE)
    WARIGH(1:LDBASE,1:LUPPER,ISPIN)=WARNEW(1:LDBASE,1:LUPPER)
!
    IF(IHFB.GE.1)THEN
       CALL ZGEMM('N','N',LDBASE,LUPPER,LUPPER,C_UNIT,&
                WARIGH(1:LDBASE,1+LUPPER:2*LUPPER,ISPIN),LDBASE,&
                OVRLAP(1:LUPPER,1       :  LUPPER),      LUPPER,&
                C_ZERO,WARNEW,LDBASE)
       WARIGH(1:LDBASE,1+LUPPER:2*LUPPER,ISPIN)=WARNEW(1:LDBASE,1:LUPPER)
    ENDIF
 ENDDO
 DEALLOCATE(WARNEW)
 RETURN
END SUBROUTINE MAT_WR
!=======================================================================
SUBROUTINE RMATCN(MREVER,ISIMPY,LDBASE,LTIMUP,LUPPER,LTOTAL,INROTA,ICHARG,&
                  RMATRX,RMAINV,WALEFT,WARIGH)
!=======================================================================
! RETURN ROTATION MATRIX FOR UPPER AND LOWER BLOCK IN CANONICAL
! BASIS, THE LATTER WILL BE INVERSED LATER
!=======================================================================
USE CANBAS
!=======================================================================
IMPLICIT NONE
!=======================================================================
 INTEGER::MREVER,ISIMPY,LDBASE,LTIMUP,LUPPER,LTOTAL,INROTA,ICHARG
 COMPLEX,DIMENSION(1:LUPPER,1:LUPPER)::RMATRX, RMAINV
 COMPLEX,DIMENSION(1:NDBASE,1:4*NDSTAT,0:NDSPIN)::WALEFT,WARIGH
!
 INTEGER:: JSTATE,ISTATE,IREVER,ISPIN,IBASE
 COMPLEX::C_ZERO!,AVRCMB,AVRCMA
!
 INTEGER::NUMBQP
 COMPLEX::PHASPI
 COMMON /SPIPHA/ PHASPI(1:NDBASE,0:NDREVE,0:NDSPIN)
 COMMON /QPNUMS/ NUMBQP(0:NDREVE,0:NDISOS)
!=======================================================================
 CALL CPUTIM('RMATCN',1)
!=======================================================================
 IF(ISIMPY.EQ.0) STOP 'NEED ISIMPY=1 IN RMATCN'
 C_ZERO=CMPLX(0.0D0,0.0D0)
 WALEFT(:,:,:)=C_ZERO
!=======================================================================
 JSTATE=0
 DO IREVER=0,MREVER
    DO ISTATE=1,LTIMUP
       JSTATE=JSTATE+1
       DO ISPIN=0,NDSPIN
!D FOR UPPER BLOCK
          WALEFT(1:LDBASE,JSTATE,ISPIN)=&
          WAVCAN(1:LDBASE,ISTATE,IREVER,ICHARG)*PHASPI(1:LDBASE,IREVER,ISPIN)
!D FOR LOWER BLOCK
          WALEFT(1:LDBASE,JSTATE+LUPPER,ISPIN)=(2*ISPIN-1)*&
          CONJG(WAVCAN(1:LDBASE,ISTATE,1-IREVER,ICHARG)*PHASPI(1:LDBASE,1-IREVER,1-ISPIN))
       END DO
    END DO
 END DO
! RD
 WARIGH(:,:,:)=WALEFT(:,:,:)
 IF(INROTA.NE.0) CALL ROTWAV(LTOTAL,WARIGH)
! DRD
 RMATRX(:,:)=C_ZERO
 RMAINV(:,:)=C_ZERO
!
 CALL WLRMAT(1,LDBASE,LUPPER,WALEFT,WARIGH,RMATRX,RMAINV)

!(rotated) wavefunction
 JSTATE=0
 DO IREVER=0,MREVER
    DO ISTATE=1,LTIMUP
       JSTATE=JSTATE+1
       DO ISPIN=0,NDSPIN
!
          WALEFT(1:LDBASE,JSTATE,ISPIN)&
         =V_CAN(ISTATE,IREVER,ICHARG)*WALEFT(1:LDBASE,JSTATE,ISPIN)
          WALEFT(1:LDBASE,JSTATE+LUPPER,ISPIN)&
         =U_CAN(ISTATE,1-IREVER,ICHARG)*WALEFT(1:LDBASE,JSTATE+LUPPER,ISPIN)
!
          WARIGH(1:LDBASE,JSTATE,ISPIN)&
         =V_CAN(ISTATE,IREVER,ICHARG)*WARIGH(1:LDBASE,JSTATE,ISPIN)
          WARIGH(1:LDBASE,JSTATE+LUPPER,ISPIN)&
         =U_CAN(ISTATE,1-IREVER,ICHARG)*WARIGH(1:LDBASE,JSTATE+LUPPER,ISPIN)
          END DO
    END DO
 END DO

!=======================================================================
 CALL CPUTIM('RMATCN',0)
!=======================================================================
 RETURN
END SUBROUTINE RMATCN
!=======================================================================
SUBROUTINE VRVMAT(MREVER,ISIMPY,LDBASE,LTIMUP,LUPPER,ICHARG, RMAT_V,RMAT_U,OVRLAP)
!=======================================================================
! calculate V*R_V*V+U*CONJG(TRANSPOSE(R_U))*U matrix in canonical basis.
!=======================================================================
USE CANBAS
!=======================================================================
IMPLICIT NONE
!=======================================================================
 INTEGER::MREVER,ISIMPY,LDBASE,LUPPER,LTIMUP,INROTA,ICHARG
 COMPLEX,DIMENSION(1:LUPPER,1:LUPPER)::OVRLAP,RMAT_V,RMAT_U
!
 INTEGER::JSTATE,ISTATE,JREVER,IREVER,IBRA,JKET,ISPIN,IBASE
 COMPLEX::C_ZERO,AVRCMP
!
 INTEGER::NUMBQP
 COMPLEX::PHASPI
 COMMON /SPIPHA/ PHASPI(1:NDBASE,0:NDREVE,0:NDSPIN)
 COMMON /QPNUMS/ NUMBQP(0:NDREVE,0:NDISOS)
!=======================================================================
 CALL CPUTIM('VRVMAT',1)
!=======================================================================
 IF(ISIMPY.EQ.0) STOP 'NEED ISIMPY=1 IN VRVMAT'
 C_ZERO=CMPLX(0.0D0,0.0D0)
!=======================================================================
 OVRLAP(:,:)=C_ZERO
! CALCULATE URU+VRV
 DO IREVER=0,MREVER!SIMPLEX BLOCK OF LEFT V,U
    DO JREVER=0,MREVER!SIMPLEX BLOCK OF RIGHT V,U
       DO ISTATE=1,LTIMUP
          DO JSTATE=1,LTIMUP
             IBRA=ISTATE+IREVER*LTIMUP!THE LABEL OF R_U,R_V, IT CAN HAVE OFF DIAGONAL PART SINCE ROTATION CAN BREAK SIMPLEX
             JKET=JSTATE+JREVER*LTIMUP
             AVRCMP=V_CAN(ISTATE,IREVER,ICHARG)*RMAT_V(IBRA,JKET)*V_CAN(JSTATE,JREVER,ICHARG)&
                             +U_CAN(ISTATE,1-IREVER,ICHARG)*CONJG(RMAT_U(JKET,IBRA))*U_CAN(JSTATE,1-JREVER,ICHARG)
             OVRLAP(IBRA,JKET)=AVRCMP
          ENDDO
       ENDDO
    ENDDO
 ENDDO
!=======================================================================
      CALL CPUTIM('VRVMAT',0)
!=======================================================================
 RETURN
END SUBROUTINE VRVMAT
!=======================================================================
SUBROUTINE ROTWLR(MREVER,ISIMPY,ICANON,IDOGOA,IPAIRI,INROTA,IEXTEN,LDBASE,LTIMUP,LUPPER,LTOTAL,ICHARG,OVKERN)
!=======================================================================
!CALCULATE B.OVERLAP^-1,A.OVERLAP^-1 IN ROTATIONAL CASE
!ICANON>1 ONLY WHEN ISIMPY=1
!=======================================================================
USE SAVLEF
USE SAVRIG
USE WAVR_L
USE PFAFFI, ONLY:NORM_N,VUAUX,LV2AUX,PV2AUX
!=======================================================================
IMPLICIT NONE
!=======================================================================
 INTEGER::MREVER,ISIMPY,ICANON,IDOGOA,IPAIRI,INROTA,IEXTEN,LDBASE,LTIMUP,LUPPER,LTOTAL,ICHARG
 COMPLEX::OVKERN
!
 INTEGER::IALLOC,ISTATE,JSTATE,KSTATE,ISPIN,IBASE,SIGNS,IFLAG,SIGV2,IFGVI2,LBOLD,LBNEW
 INTEGER::LOCCUD(0:NDREVE)
 INTEGER::IAUXDI(1:2*NDSTAT)
 REAL   ::NOTOCC,RICOND,LOGNORM,PHASE,LOGPF,LOGV2,phav2
 COMPLEX::C_ZERO,C_UNIT,C_I,AVRCMP,CPHASE
 COMPLEX::DETWRK(2),DETRMA(2),AUXDIA(1:2*NDSTAT)
 COMPLEX,DIMENSION(1:2*NDSTAT,1:2*NDSTAT)::VU
 COMPLEX,ALLOCATABLE::OVRLAP(:,:),RMATRX(:,:),RMAINV(:,:),EWAVEF(:,:),WARAUX(:,:)
!
 INTEGER::NFIPRI
 REAL   ::ROTATX,ROTATY,ROTATZ
 COMPLEX::PHASPI
 COMMON /CFIPRI/ NFIPRI
 COMMON /ROTATT/ ROTATX,ROTATY,ROTATZ
 COMMON /SPIPHA/ PHASPI(1:NDBASE,0:NDREVE,0:NDSPIN)
!=======================================================================
 IALLOC=1
!=======================================================================
 ALLOCATE (OVRLAP(1:2*NDSTAT,1:2*NDSTAT),STAT=IALLOC)
 IF (IALLOC.NE.0) CALL NOALLO('OVRLAP','ROTWLR')
!
 IF(ICANON.EQ.2)THEN
    ALLOCATE (RMATRX(1:2*NDSTAT,1:2*NDSTAT),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('RMATRX','ROTWLR')
    ALLOCATE (RMAINV(1:2*NDSTAT,1:2*NDSTAT),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('RMAINV','ROTWLR')
 ENDIF
!=======================================================================
 C_ZERO=CMPLX(0.0D0,0.0D0)
 C_UNIT=CMPLX(1.0D0,0.0D0)
 C_I=CMPLX(0.0D0,1.0D0)
 LOCCUD(:)=LTIMUP
!=======================================================================
!(BRB+ARA) IN HFB CASE
!=======================================================================
 IF(IPAIRI.EQ.1)THEN

    IF(ICANON.LT.2)THEN
!=======================================================================
!PUT WAVE FUNCTION INTO WALEFT AND WARIGH
!=======================================================================
       IF(ICANON.EQ.1)THEN
          CALL CANWAV(MREVER,ISIMPY,WALEFT,ICHARG,LDBASE,LOCCUD)
          WARIGH(:,:,:)=WALEFT(:,:,:)
       ELSE
          WALEFT(:,:,:)=SARIGH(:,:,:,ICHARG)
          WARIGH(:,:,:)=SARIGH(:,:,:,ICHARG)
       ENDIF
!=======================================================================
!ROTATE WARIGH
!=======================================================================
       IF(INROTA.NE.0) CALL ROTWAV(LTOTAL,WARIGH)
!=======================================================================
!/CALCULATE OVERLAP MATRIX
!=======================================================================
       OVRLAP(:,:)=C_ZERO
       IAUXDI(:)=0
       AUXDIA(:)=0.0D0
!
       ALLOCATE(VRVAUX(1:LUPPER,1:LUPPER))
       ALLOCATE(URUAUX(1:LUPPER,1:LUPPER))
!
       CALL WLRMAT(IPAIRI,LDBASE,LUPPER,WALEFT,WARIGH,VRVAUX,URUAUX)
       OVRLAP(1:LUPPER,1:LUPPER)=VRVAUX(1:LUPPER,1:LUPPER)+URUAUX(1:LUPPER,1:LUPPER)
!
       DEALLOCATE(VRVAUX,URUAUX)
!=======================================================================
! INVERTING THE OVERLAP MATRIX
!=======================================================================
       CALL ZGECO(OVRLAP,2*NDSTAT,LUPPER,IAUXDI,RICOND,AUXDIA)
       CALL ZGEDI(OVRLAP,2*NDSTAT,LUPPER,IAUXDI,DETWRK,AUXDIA,01)
       IF(IDOGOA.NE.0)THEN!GOA NEED THE KERNAL
          IF(MREVER.EQ.0) STOP 'CALCULATE KERNAL WITH PFAFFIAN IN HFB CASE NEED MREVER=1 IN ROTWLR'
          !USE PFAFFIAN METHOD TO CALCULATE THE KERNAL
           IF(INROTA.NE.0)THEN
              ALLOCATE(VUAUX(1:2*NDSTAT,1:2*NDSTAT,0:NDISOS))
              CALL NORM_N(0,IPAIRI,LUPPER,LDBASE,WALEFT,WALEFT,VUAUX(:,:,ICHARG),LV2AUX(ICHARG),PV2AUX(ICHARG),OVKERN)
              CALL NORM_N(1,IPAIRI,LUPPER,LDBASE,WALEFT,WARIGH,VUAUX(:,:,ICHARG),LV2AUX(ICHARG),PV2AUX(ICHARG),OVKERN)
!in hfb case, pfaffian is the overlap; in hf case, determinant of d^+ d is the overlap
              DEALLOCATE(VUAUX)
           ELSE!the kernal is 1 if there is no rotation
              OVKERN=1.0d0
           ENDIF
       ENDIF
!
       IF (ABS(RICOND).LT.1.0D-12) THEN
          WRITE(NFIPRI,*)'BE CAREFUL - R MATRIX ,MAY BE SINGULAR IN RENINE FOR:'
          WRITE(NFIPRI,*)'INROTA=',INROTA,' RICOND=',RICOND
       END IF
!
    ELSEIF(ICANON.EQ.2)THEN
!=======================================================================
! CONSTRUT R AND R^-1 MATRIX UNDER CANONICAL BASIS
! ALSO CALCULATE WAVE FUNCTION AND ROTATED WAVE FUNCTION
!=======================================================================
       CALL RMATCN(MREVER,ISIMPY,LDBASE,LTIMUP,LUPPER,LTOTAL,INROTA,ICHARG,&
                   RMATRX(1:LUPPER,1:LUPPER),RMAINV(1:LUPPER,1:LUPPER),&
                   WALEFT,WARIGH)
!
       IAUXDI(:)=0
       AUXDIA(:)=0.0D0
       CALL ZGECO(RMAINV,2*NDSTAT,LUPPER,IAUXDI,RICOND,AUXDIA)
       CALL ZGEDI(RMAINV,2*NDSTAT,LUPPER,IAUXDI,DETRMA,AUXDIA,11)
!
       IF (ABS(RICOND).LT.1.0D-12) THEN
          WRITE(NFIPRI,*)'BE CAREFUL - R MATRIX ,MAY BE SINGULAR IN RENINE FOR:'
          WRITE(NFIPRI,*)'INROTA=',INROTA,' RICOND=',RICOND
       END IF
!=======================================================================
!(BRB+AR^-1A) IN HFB CASE
!=======================================================================
       CALL VRVMAT(MREVER,ISIMPY,LDBASE,LTIMUP,LUPPER,ICHARG,&
                   RMATRX(1:LUPPER,1:LUPPER),RMAINV(1:LUPPER,1:LUPPER),&
                   OVRLAP(1:LUPPER,1:LUPPER))
!=======================================================================
! INVERTING THE OVERLAP MATRIX
!=======================================================================
       IAUXDI(:)=0
       AUXDIA(:)=0.0D0
       CALL ZGECO(OVRLAP,2*NDSTAT,LUPPER,IAUXDI,RICOND,AUXDIA)
       CALL ZGEDI(OVRLAP,2*NDSTAT,LUPPER,IAUXDI,DETWRK,AUXDIA,11)
!
       IF (ABS(RICOND).LT.1.0D-12) THEN
          WRITE(NFIPRI,*)'BE CAREFUL - R MATRIX ,MAY BE SINGULAR IN RENINE FOR:'
          WRITE(NFIPRI,*)'INROTA=',INROTA,' RICOND=',RICOND
       END IF
!
       OVKERN=DETWRK(1)*10.0D0**DETWRK(2)*(DETRMA(1)*10.0D0**DETRMA(2))

    ELSE
          STOP 'WRONG ICANO IN RENINE'

    ENDIF
!=======================================================================
!MULTIPLYING THE RIGHT WAVE FUNCTIONS BY THE INVERTED OVERLAP MATRIX
!=======================================================================
    CALL MAT_WR(IPAIRI,LDBASE,LUPPER,OVRLAP(1:LUPPER,1:LUPPER),WARIGH)
 ELSE
!=======================================================================
!CALCULATE B(BRB+ARA)^-1, A(BRB+ARA)^-1 IN HF CASE
!=======================================================================
    WALEFT(:,:,:)=SARIGH(:,:,:,ICHARG)
    WARIGH(:,:,:)=SARIGH(:,:,:,ICHARG)
    IF(INROTA.NE.0) CALL ROTWAV(LTOTAL,WARIGH)
!
    CALL WLRMAT(IPAIRI,LDBASE,LTOTAL,WALEFT,WARIGH,OVRLAP(1:LTOTAL,1:LTOTAL))
!
!=======================================================================
! INVERTING THE OVERLAP MATRIX
!=======================================================================
    CALL ZGECO(OVRLAP,2*NDSTAT,LTOTAL,IAUXDI,RICOND,AUXDIA)
!
    IF(IDOGOA.GE.1)THEN
       CALL ZGEDI(OVRLAP,2*NDSTAT,LTOTAL,IAUXDI,DETWRK,AUXDIA,11)
       OVKERN=DETWRK(1)*10.0D0**DETWRK(2)
    ELSE
       CALL ZGEDI(OVRLAP,2*NDSTAT,LTOTAL,IAUXDI,DETWRK,AUXDIA,01)
    ENDIF
!
    IF (ABS(RICOND).LT.1.0D-12) THEN
        WRITE(NFIPRI,*) 'BE CAREFULL - OVERLAP MATRIX ,MAY BE SINGULAR IN RENINE FOR:'
        WRITE(NFIPRI,*) 'INROTA=',INROTA,' RICOND=',RICOND
    END IF
!=======================================================================
!MULTIPLYING THE RIGHT WAVE FUNCTIONS BY THE INVERTED OVERLAP MATRIX
!=======================================================================
    CALL MAT_WR(IPAIRI,LDBASE,LTOTAL,OVRLAP(1:LTOTAL,1:LTOTAL),WARIGH)
 END IF ! End IF (IPAIRI.EQ.1)
!=======================================================================
 DEALLOCATE(OVRLAP)
 IF(ICANON.EQ.2)DEALLOCATE(RMATRX,RMAINV)
!=======================================================================
END SUBROUTINE ROTWLR
!=======================================================================
SUBROUTINE LMATCN(MREVER,ISIMPY,LDBASE,LTIMUP,LUPPER,LTOTAL,INSHIF,ICHARG,&
                  LMATRX,LMAINV,WALEFT,WARIGH)
!=======================================================================
! RETURN ROTATION MATRIX FOR UPPER AND LOWER BLOCK IN CANONICAL
! BASIS, THE LATTER WILL BE INVERSED LATER
!=======================================================================
USE CANBAS
!=======================================================================
IMPLICIT NONE
!=======================================================================
 INTEGER::MREVER,ISIMPY,LDBASE,LTIMUP,LUPPER,LTOTAL,INSHIF,ICHARG
 COMPLEX,DIMENSION(1:LUPPER,1:LUPPER)::LMATRX, LMAINV
 COMPLEX,DIMENSION(1:NDBASE,1:4*NDSTAT,0:NDSPIN)::WALEFT,WARIGH
!
 INTEGER:: JSTATE,ISTATE,IREVER,ISPIN,IBASE
 COMPLEX::C_ZERO,AVRCMB,AVRCMA
!
 INTEGER::NUMBQP
 COMPLEX::PHASPI
 COMMON /SPIPHA/ PHASPI(1:NDBASE,0:NDREVE,0:NDSPIN)
 COMMON /QPNUMS/ NUMBQP(0:NDREVE,0:NDISOS)
!=======================================================================
 CALL CPUTIM('LMATCN',1)
!=======================================================================
 IF(ISIMPY.EQ.0) STOP 'NEED ISIMPY=1 IN LMATCN'
 C_ZERO=CMPLX(0.0D0,0.0D0)
 WALEFT(:,:,:)=C_ZERO
!=======================================================================
 JSTATE=0
 DO IREVER=0,MREVER
    DO ISTATE=1,LTIMUP
       JSTATE=JSTATE+1
       DO ISPIN=0,NDSPIN
!D FOR UPPER BLOCK
          WALEFT(1:LDBASE,JSTATE,ISPIN)=&
          WAVCAN(1:LDBASE,ISTATE,IREVER,ICHARG)*PHASPI(1:LDBASE,IREVER,ISPIN)
!D FOR LOWER BLOCK
          WALEFT(1:LDBASE,JSTATE+LUPPER,ISPIN)=(2*ISPIN-1)*&
          CONJG(WAVCAN(1:LDBASE,ISTATE,1-IREVER,ICHARG)*PHASPI(1:LDBASE,1-IREVER,1-ISPIN))
       END DO
    END DO
 END DO
! RD
 WARIGH(:,:,:)=WALEFT(:,:,:)
 IF(INSHIF.NE.0) CALL SHIWAV(LTOTAL,WARIGH)
! DRD
 LMATRX(:,:)=C_ZERO
 LMAINV(:,:)=C_ZERO
!
 CALL WLRMAT(1,LDBASE,LUPPER,WALEFT,WARIGH,LMATRX,LMAINV)
!(shifted) wavefunction
 JSTATE=0
 DO IREVER=0,MREVER
    DO ISTATE=1,LTIMUP
       JSTATE=JSTATE+1
       DO ISPIN=0,NDSPIN
!
          WALEFT(1:LDBASE,JSTATE,ISPIN)&
         =V_CAN(ISTATE,IREVER,ICHARG)*WALEFT(1:LDBASE,JSTATE,ISPIN)
          WALEFT(1:LDBASE,JSTATE+LUPPER,ISPIN)&
         =U_CAN(ISTATE,1-IREVER,ICHARG)*WALEFT(1:LDBASE,JSTATE+LUPPER,ISPIN)
!
          WARIGH(1:LDBASE,JSTATE,ISPIN)&
         =V_CAN(ISTATE,IREVER,ICHARG)*WARIGH(1:LDBASE,JSTATE,ISPIN)
          WARIGH(1:LDBASE,JSTATE+LUPPER,ISPIN)&
         =U_CAN(ISTATE,1-IREVER,ICHARG)*WARIGH(1:LDBASE,JSTATE+LUPPER,ISPIN)
          END DO
    END DO
 END DO

!=======================================================================
 CALL CPUTIM('LMATCN',0)
!=======================================================================
 RETURN
END SUBROUTINE LMATCN
!=======================================================================
SUBROUTINE LINWLR(MREVER,ISIMPY,ICANON,IDOGOA,IPAIRI,INSHIF,IEXTEN,LDBASE,LTIMUP,LUPPER,LTOTAL,ICHARG,OVKERN)
!=======================================================================
!CALCULATE B.OVERLAP^-1,A.OVERLAP^-1 IN ROTATIONAL CASE
!ICANON>1 ONLY WHEN ISIMPY=1
!=======================================================================
USE SAVLEF
USE SAVRIG
USE WAVR_L
USE PFAFFI, ONLY:NORM_N,VUAUX,LV2AUX,PV2AUX
!=======================================================================
IMPLICIT NONE
!=======================================================================
 INTEGER::MREVER,ISIMPY,ICANON,IDOGOA,IPAIRI,INSHIF,IEXTEN,LDBASE,LTIMUP,LUPPER,LTOTAL,ICHARG
 COMPLEX::OVKERN
!
 INTEGER::IALLOC,ISTATE,JSTATE,KSTATE,ISPIN,IBASE
 INTEGER::LOCCUD(0:NDREVE)
 INTEGER::IAUXDI(1:2*NDSTAT)
 REAL   ::RICOND
 COMPLEX::C_ZERO,C_UNIT,AVRCMP
 COMPLEX::DETWRK(2),DETRMA(2),AUXDIA(1:2*NDSTAT)
 COMPLEX,ALLOCATABLE::OVRLAP(:,:),LMATRX(:,:),LMAINV(:,:),EWAVEF(:,:),WARAUX(:,:)
!
 INTEGER::NFIPRI
 COMMON /CFIPRI/ NFIPRI
!=======================================================================
 IALLOC=1
!=======================================================================
 ALLOCATE (OVRLAP(1:2*NDSTAT,1:2*NDSTAT),STAT=IALLOC)
 IF (IALLOC.NE.0) CALL NOALLO('OVRLAP','LINWLR')
!
 IF(ICANON.EQ.2)THEN
    ALLOCATE (LMATRX(1:2*NDSTAT,1:2*NDSTAT),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('LMATRX','LINWLR')
    ALLOCATE (LMAINV(1:2*NDSTAT,1:2*NDSTAT),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('LMAINV','LINWLR')
 ENDIF
!=======================================================================
 C_ZERO=CMPLX(0.0D0,0.0D0)
 C_UNIT=CMPLX(1.0D0,0.0D0)
 LOCCUD(:)=LTIMUP
!=======================================================================
!(BRB+ARA) IN HFB CASE
!=======================================================================
 IF(IPAIRI.EQ.1)THEN

    IF(ICANON.LT.2)THEN
!=======================================================================
!PUT WAVE FUNCTION INTO WALEFT AND WARIGH
!=======================================================================
       IF(ICANON.EQ.1)THEN
          CALL CANWAV(MREVER,ISIMPY,WALEFT,ICHARG,LDBASE,LOCCUD)
          WARIGH(:,:,:)=WALEFT(:,:,:)
       ELSE
          WALEFT(:,:,:)=SARIGH(:,:,:,ICHARG)
          WARIGH(:,:,:)=SARIGH(:,:,:,ICHARG)
       ENDIF
!=======================================================================
!ROTATE WARIGH
!=======================================================================
       IF(INSHIF.NE.0) CALL SHIWAV(LTOTAL,WARIGH)
!=======================================================================
!CALCULATE OVERLAP MATRIX
!=======================================================================
       OVRLAP(:,:)=C_ZERO
       IAUXDI(:)=0
       AUXDIA(:)=0.0D0
!
       DO JSTATE=1,LUPPER
          DO ISTATE=1,LUPPER
             AVRCMP=C_ZERO
             DO ISPIN=0,NDSPIN
                DO IBASE=1,LDBASE
                   AVRCMP=AVRCMP&
            +CONJG(WALEFT(IBASE,ISTATE,ISPIN))*WARIGH(IBASE,JSTATE,ISPIN)&
            +CONJG(WALEFT(IBASE,ISTATE+LUPPER,ISPIN))*WARIGH(IBASE,JSTATE+LUPPER,ISPIN)
                END DO
             END DO
             OVRLAP(ISTATE,JSTATE)=AVRCMP
          END DO
       END DO
!=======================================================================
! INVERTING THE OVERLAP MATRIX
!=======================================================================
       CALL ZGECO(OVRLAP,2*NDSTAT,LUPPER,IAUXDI,RICOND,AUXDIA)
       CALL ZGEDI(OVRLAP,2*NDSTAT,LUPPER,IAUXDI,DETWRK,AUXDIA,01)
       IF(IDOGOA.NE.0)THEN!GOA NEED THE KERNAL
          IF(MREVER.EQ.0) STOP 'CALCULATE KERNAL WITH PFAFFIAN IN HFB CASE NEED MREVER=1 IN LINWLR'
          !USE PFAFFIAN METHOD TO CALCULATE THE KERNAL
           IF(INSHIF.NE.0)THEN
              ALLOCATE(VUAUX(1:2*NDSTAT,1:2*NDSTAT,0:NDISOS))
              CALL NORM_N(0,IPAIRI,LUPPER,LDBASE,WALEFT,WALEFT,VUAUX(:,:,ICHARG),LV2AUX(ICHARG),PV2AUX(ICHARG),OVKERN)
              CALL NORM_N(1,IPAIRI,LUPPER,LDBASE,WALEFT,WARIGH,VUAUX(:,:,ICHARG),LV2AUX(ICHARG),PV2AUX(ICHARG),OVKERN)
!in hfb case, pfaffian is the overlap; in hf case, determinant of d^+ d is the overlap
              DEALLOCATE(VUAUX)
           ELSE
              OVKERN=1.0d0
           ENDIF
       ENDIF
!
       IF (ABS(RICOND).LT.1.0D-12) THEN
          WRITE(NFIPRI,*)'BE CAREFUL - R MATRIX ,MAY BE SINGULAR IN RENINE FOR:'
          WRITE(NFIPRI,*)'INSHIF=',INSHIF,' RICOND=',RICOND
       END IF
!
    ELSEIF(ICANON.EQ.2)THEN
!=======================================================================
! CONSTRUT R AND R^-1 MATRIX UNDER CANONICAL BASIS
! ALSO CALCULATE WAVE FUNCTION AND ROTATED WAVE FUNCTION
!=======================================================================
       CALL LMATCN(MREVER,ISIMPY,LDBASE,LTIMUP,LUPPER,LTOTAL,INSHIF,ICHARG,&
                   LMATRX(1:LUPPER,1:LUPPER),LMAINV(1:LUPPER,1:LUPPER),&
                   WALEFT,WARIGH)
!
       IAUXDI(:)=0
       AUXDIA(:)=0.0D0
       CALL ZGECO(LMAINV,2*NDSTAT,LUPPER,IAUXDI,RICOND,AUXDIA)
       CALL ZGEDI(LMAINV,2*NDSTAT,LUPPER,IAUXDI,DETRMA,AUXDIA,11)
!
       IF (ABS(RICOND).LT.1.0D-12) THEN
          WRITE(NFIPRI,*)'BE CAREFUL - R MATRIX ,MAY BE SINGULAR IN RENINE FOR:'
          WRITE(NFIPRI,*)'INSHIF=',INSHIF,' RICOND=',RICOND
       END IF
!=======================================================================
!(BRB+AR^-1A) IN HFB CASE
!=======================================================================
       CALL VRVMAT(MREVER,ISIMPY,LDBASE,LTIMUP,LUPPER,ICHARG,&
                   LMATRX(1:LUPPER,1:LUPPER),LMAINV(1:LUPPER,1:LUPPER),&
                   OVRLAP(1:LUPPER,1:LUPPER))
! write(*,*)ovrlap(1,1),ovrlap(1,1+ltimup),'ovrlap'
!=======================================================================
! INVERTING THE OVERLAP MATRIX
!=======================================================================
       IAUXDI(:)=0
       AUXDIA(:)=0.0D0
       CALL ZGECO(OVRLAP,2*NDSTAT,LUPPER,IAUXDI,RICOND,AUXDIA)
       CALL ZGEDI(OVRLAP,2*NDSTAT,LUPPER,IAUXDI,DETWRK,AUXDIA,11)
!
       IF (ABS(RICOND).LT.1.0D-12) THEN
          WRITE(NFIPRI,*)'BE CAREFUL - R MATRIX ,MAY BE SINGULAR IN RENINE FOR:'
          WRITE(NFIPRI,*)'INSHIF=',INSHIF,' RICOND=',RICOND
       END IF
!
       OVKERN=DETWRK(1)*10.0D0**DETWRK(2)*(DETRMA(1)*10.0D0**DETRMA(2))
    ELSE
          STOP 'WRONG ICANO IN RENINE'

    ENDIF
!=======================================================================
!MULTIPLYING THE RIGHT WAVE FUNCTIONS BY THE INVERTED OVERLAP MATRIX
!=======================================================================
    ALLOCATE(EWAVEF(1:NDBASE,1:2*NDSTAT))
    ALLOCATE(WARAUX(1:NDBASE,1:2*NDSTAT))
!
    DO ISPIN=0,NDSPIN
       EWAVEF(:,:)=C_ZERO
       WARAUX(:,:)=C_ZERO
!UPPER BLOCK
       DO ISTATE=1,LUPPER
          CALL ZCOPY(NDBASE,WARIGH(1,ISTATE,ISPIN),1,EWAVEF(1,ISTATE),1)
       END DO
       CALL ZGEMM('N','N',LDBASE,LUPPER,LUPPER,&
                   C_UNIT,EWAVEF,NDBASE,OVRLAP,2*NDSTAT,&
                   C_ZERO,WARAUX,NDBASE)
       DO ISTATE=1,LUPPER
          CALL ZCOPY(NDBASE,WARAUX(1,ISTATE),1,WARIGH(1,ISTATE,ISPIN),1)
       END DO
!LOWER BLOCK
       DO ISTATE=1,LUPPER
          CALL ZCOPY(NDBASE,WARIGH(1,ISTATE+LUPPER,ISPIN),1,EWAVEF(1,ISTATE),1)
       END DO
       CALL ZGEMM('N','N',LDBASE,LUPPER,LUPPER,&
                   C_UNIT,EWAVEF,NDBASE,OVRLAP,2*NDSTAT,&
                   C_ZERO,WARAUX,NDBASE)
       DO ISTATE=1,LUPPER
          CALL ZCOPY(NDBASE, WARAUX(1,ISTATE),1,WARIGH(1,ISTATE+LUPPER,ISPIN),1)
       END DO
    END DO
!
    DEALLOCATE(EWAVEF)
    DEALLOCATE(WARAUX)
 ELSE
!=======================================================================
!CALCULATE B(BRB+ARA)^-1, A(BRB+ARA)^-1 IN HF CASE
!=======================================================================
    WALEFT(:,:,:)=SARIGH(:,:,:,ICHARG)
    WARIGH(:,:,:)=SARIGH(:,:,:,ICHARG)
    IF(INSHIF.NE.0) CALL SHIWAV(LTOTAL,WARIGH)
!
    DO JSTATE=1,LTOTAL
       DO ISTATE=1,LTOTAL
          AVRCMP=C_ZERO
          DO ISPIN=0,NDSPIN
             DO IBASE=1,LDBASE
                AVRCMP=AVRCMP+CONJG(WALEFT(IBASE,ISTATE,ISPIN))*WARIGH(IBASE,JSTATE,ISPIN)
             END DO
          END DO
          OVRLAP(ISTATE,JSTATE)=AVRCMP
       END DO
    END DO
!=======================================================================
! INVERTING THE OVERLAP MATRIX
!=======================================================================
    CALL ZGECO(OVRLAP,2*NDSTAT,LTOTAL,IAUXDI,RICOND,AUXDIA)
    IF(IDOGOA.GE.1)THEN
       CALL ZGEDI(OVRLAP,2*NDSTAT,LTOTAL,IAUXDI,DETWRK,AUXDIA,11)
       OVKERN=DETWRK(1)*10.0D0**DETWRK(2)
    ELSE
       CALL ZGEDI(OVRLAP,2*NDSTAT,LTOTAL,IAUXDI,DETWRK,AUXDIA,01)
    ENDIF
    IF (ABS(RICOND).LT.1.0D-12) THEN
        WRITE(NFIPRI,*) 'BE CAREFULL - OVERLAP MATRIX ,MAY BE SINGULAR IN RENINE FOR:'
        WRITE(NFIPRI,*) 'INSHIF=',INSHIF,' RICOND=',RICOND
    END IF
!=======================================================================
!MULTIPLYING THE RIGHT WAVE FUNCTIONS BY THE INVERTED OVERLAP MATRIX
!=======================================================================
    ALLOCATE (WARAUX(1:2*NDSTAT,0:NDSPIN),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('WARAUX','LINWRL')

    DO JSTATE=1,LDBASE
       DO ISTATE=1,LTOTAL
          DO ISPIN=0,NDSPIN
             WARAUX(ISTATE,ISPIN)=C_ZERO
             DO KSTATE=1,LTOTAL
                WARAUX(ISTATE,ISPIN)=WARAUX(ISTATE,ISPIN)&
          +WARIGH(JSTATE,KSTATE,ISPIN)*OVRLAP(KSTATE,ISTATE)
             END DO
          END DO
       END DO
!
       DO ISTATE=1,LTOTAL
          DO ISPIN=0,NDSPIN
              WARIGH(JSTATE,ISTATE,ISPIN)=WARAUX(ISTATE,ISPIN)
          END DO
       END DO
    END DO
    DEALLOCATE(WARAUX)
 END IF ! End IF (IPAIRI.EQ.1)
!=======================================================================
 DEALLOCATE(OVRLAP)
 IF(ICANON.EQ.2)DEALLOCATE(LMATRX,LMAINV)
!=======================================================================
END SUBROUTINE LINWLR
!=======================================================================
SUBROUTINE H_KERN(NXHERM,NYHERM,NZHERM,NXMAXX,NYMAXX,NZMAXX,&
                  ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY,IPAHFB,MREVER,&
                  ICHARG,MIN_QP,IPNMIX,ITIREP,&
                  NAMEPN,PRINIT,IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR,ISAWAV,IKERNE,&
                  ISHIFY,JMULMO,NUMCOU,BOUCOU,IN_FIX,IZ_FIX,&
                  ISKYRM,EKEKIN,EPAIKI,DROTSQ,DKOTSQ,EKECOD,EKECOX,EKESKY)
!=======================================================================
!CALCULATE RHO, KINETIC, PAIRING AND COULOMB ENERGY
!=======================================================================
IMPLICIT NONE
!=======================================================================
 CHARACTER::NAMEPN*8
 LOGICAL::PRINIT
 INTEGER::NXHERM,NYHERM,NZHERM,NXMAXX,NYMAXX,NZMAXX,&
          ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY,IPAHFB,MREVER,&
          ICHARG,MIN_QP,IPNMIX,ITPNMX,ITIREP,&
          IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR,ISAWAV,IKERNE,&
          ISHIFY,JMULMO,NUMCOU,IN_FIX,IZ_FIX,ISKYRM
 REAL   ::BOUCOU,ENEKIN,DROTSQ
 COMPLEX::EKEKIN,DKOTSQ,EPAIKI,EKECOD,EKECOX,EKESKY
!
 INTEGER::LDPNMX
 REAL   ::COMULT,EKIN_N,EKIN_P,EKIN_T,EPOT_N,EPOT_P,EPOT_T,ESUM_N,ESUM_P,ESUM_T,EPAI_N,EPAI_P,EPAI_T,&
          EREA_N,EREA_P,EREA_T,ELIP_N,ELIP_P,ELIP_T,ECOULD,ECOULE,ECOULT,ECOULS,ECOULV, EMULCO,&
          EMUSLO,EMUREA,ESIFCO,ESISLO,ESIREA,ESPICO,ESPSLO,ESPREA,ENREAR,ECORCM,ECOR_R,EEVEW0,&
          EODDW0,ENE_W0,ENEVEN,ENEODD,ENESKY,ESTABN,ETOTSP,ETOTFU
 COMPLEX::QMUT_N,QMUT_P,QMUT_T,QMUL_N,QMUL_P,QMUL_T
 COMMON /QMUTTI/ QMUT_N(0:NDMULT,-NDMULT:NDMULT),QMUT_P(0:NDMULT,-NDMULT:NDMULT),QMUT_T(0:NDMULT,-NDMULT:NDMULT)
 COMMON /QMULTI/ QMUL_N(0:NDMULT,-NDMULT:NDMULT),QMUL_P(0:NDMULT,-NDMULT:NDMULT),QMUL_T(0:NDMULT,-NDMULT:NDMULT)
 COMMON /COEMUL/ COMULT(0:NDMULT,0:NDOSCI,0:NDOSCI,1:NDKART)
 COMMON /ALLENE/ EKIN_N,EKIN_P,EKIN_T,EPOT_N,EPOT_P,EPOT_T,ESUM_N,ESUM_P,ESUM_T,&
                 EPAI_N,EPAI_P,EPAI_T,EREA_N,EREA_P,EREA_T,ELIP_N,ELIP_P,ELIP_T,&
                 ECOULD,ECOULE,ECOULT,ECOULS,ECOULV,&
                 EMULCO,EMUSLO,EMUREA,ESIFCO,ESISLO,ESIREA,ESPICO,ESPSLO,ESPREA,&
                 ENREAR,ECORCM,ECOR_R,&
                 EEVEW0,EODDW0,ENE_W0,ENEVEN,ENEODD,ENESKY,ESTABN,ETOTSP,ETOTFU
!
 COMPLEX::C_ZERO
!=======================================================================
 ITPNMX=ICHARG
 C_ZERO=CMPLX(0.0D0,0.0D0)
!=======================================================================
!              DENSITY
!=======================================================================
 CALL ZEDENS(ITPNMX)
 CALL DENSHF(NXHERM,NYHERM,NZHERM,ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY,IPAHFB,MREVER,&
             ICHARG,MIN_QP,IPNMIX,ITPNMX,ITIREP,NXMAXX,&
             NAMEPN,PRINIT,IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR,ISAWAV,IKERNE)
 IF(MREVER.EQ.0) CALL DBLING(ITPNMX,NXHERM,NYHERM,NZHERM)
!=======================================================================
!              KINETIC ENERGY
!=======================================================================
! CALL DENMAC(MREVER,ICHARG,ISIMPY,IPAHFB,WALEFT,WARIGH)!CALL BEFORE THIS SUBROUTINE
 CALL EKINET(ENEKIN,DROTSQ,EKEKIN,DKOTSQ)
!=======================================================================
!              PAIRING ENERGY
!=======================================================================
 IF (IPAHFB.EQ.1) THEN
    CALL EPAIRI(NXHERM,NYHERM,NZHERM,IN_FIX,IZ_FIX,ITPNMX,EPAIKI)
    IF (MREVER.EQ.0) EPAIKI=EPAIKI*4.0d0
 ELSE
    EPAIKI=C_ZERO
 END IF
!
 IF(ICHARG.EQ.1)THEN
! if(.false.)then
!=======================================================================
!              CALCULATING THE COULOMB ENERGY
!=======================================================================
    CALL MOMETS(ISIMPY,ISIGNY,ISIQTY,QMUL_P,QMUT_P,COMULT,JMULMO,ISHIFY)
    CALL COUMAT(NXHERM,NYHERM,NZHERM,NXMAXX,NYMAXX,NZMAXX,NUMCOU,BOUCOU,ISIMPY,IKERNE)
    CALL COULOD(ISIMPY,EKECOD)
    CALL TRUCHD(NXHERM,NYHERM,NZHERM)
    CALL COULOE(NXHERM,NYHERM,NZHERM,EKECOX)
 ENDIF
!
 IF(ISKYRM.EQ.1)THEN
!=======================================================================
!              CALCULATING THE SKYRME ENERGY
!=======================================================================
       LDPNMX=1
       CALL TRUTOD(NXHERM,NYHERM,NZHERM)
       CALL ESKYRM(NXHERM,NYHERM,NZHERM,ENESKY,ENEVEN,ENEODD,&
                   ENREAR,ENE_W0,EEVEW0,EODDW0,EKESKY,LDPNMX)
 ENDIF
 RETURN
END SUBROUTINE H_KERN
!=======================================================================
END MODULE OFFDIA
!=======================================================================
!ROUTINES CALLED IN PROGRAM MAIN
SUBROUTINE CANWAV(MREVER,ISIMPY,WALEFT,ICHARG,LDBASE,LOCCUD)
!=======================================================================
! use u_can, v_can ,wavcan to build A=DU, B=DV
!=======================================================================
USE CANBAS
!=======================================================================
USE hfodd_sizes
!=======================================================================
IMPLICIT NONE
!=======================================================================
 INTEGER::MREVER,ISIMPY,ICHARG,LDBASE,LTIMUP,LUPPER,LSTATE,LTOTAL
 INTEGER::LOCCUD(0:NDREVE)
!
 INTEGER::JSTATE,IREVER,ISTATE,ISPIN
 COMPLEX::C_ZERO
 COMPLEX::WALEFT(1:NDBASE,1:4*NDSTAT,0:NDSPIN)
!
 INTEGER::NUMBQP
 COMPLEX::PHASPI
 COMMON /SPIPHA/ PHASPI(1:NDBASE,0:NDREVE,0:NDSPIN)
 COMMON /QPNUMS/ NUMBQP(0:NDREVE,0:NDISOS)
!=======================================================================
 CALL CPUTIM('CANWAV',1)
!=======================================================================
 IF(ISIMPY.EQ.0) STOP 'NOT IMPLEMENT YET'
!=======================================================================
 C_ZERO=CMPLX(0.0D0,0.0D0)
 WALEFT(:,:,:)=C_ZERO
!=======================================================================
 JSTATE=0
! B
! do istate=1,loccud(0)
! write(*,*)v_can(istate,0,icharg)**2+u_can(istate,1,icharg)**2,&
!           v_can(istate,0,icharg)**2+u_can(istate,0,icharg)**2,istate
! enddo
! stop
! V_CAN(:,1,:)=-V_CAN(:,0,:)
! U_CAN(:,1,:)=U_CAN(:,0,:)
 !WAVCAN(:,:,1,:)=CONJG(WAVCAN(:,:,0,:))
 DO IREVER=0,MREVER
    DO ISTATE=1,LOCCUD(IREVER)
          JSTATE=JSTATE+1
          DO ISPIN=0,NDSPIN
             WALEFT(1:LDBASE,JSTATE,ISPIN)=V_CAN(ISTATE,IREVER,ICHARG)&
           * WAVCAN(1:LDBASE,ISTATE,IREVER,ICHARG)* PHASPI(1:LDBASE,IREVER,ISPIN)
          END DO
    END DO
 END DO
! A
 DO IREVER=1,(1-MREVER),-1
     DO ISTATE=1,LOCCUD(IREVER)
          JSTATE=JSTATE+1
          DO ISPIN=0,NDSPIN
             WALEFT(1:LDBASE,JSTATE,ISPIN)=U_CAN(ISTATE,1-IREVER,ICHARG) &
          * (2*ISPIN-1)*CONJG(WAVCAN(1:LDBASE,ISTATE,IREVER,ICHARG)* PHASPI(1:LDBASE,IREVER,1-ISPIN))
          END DO
     END DO
 END DO
!=======================================================================
      CALL CPUTIM('CANWAV',0)
!=======================================================================
 RETURN
END SUBROUTINE CANWAV
!=======================================================================
SUBROUTINE ROTAVR(NXMAXX,NYMAXX,NZMAXX,ISIMPY,ISIGNY,ISIQTY,&
                  IPAIRI,ICHARG,KETA_R,IROTAT,MREVER,CORROT,&
                  DROTSQ,EROTSQ,PROTSQ,TROTSQ,AROTLI,PROTLI,PROTKI,&
                  DKOTSQ,EKOTSQ,PKOTSQ,TKOTSQ,AKOTLI,PKOTLI,PKOTKI)
!=======================================================================
USE MAT_PP
USE MAT_PM
USE MAD_PP
USE MAD_PM
USE KAP_PP
USE KAP2PP
USE KAP_PM
USE MAP_PP
USE MAP_PM
USE AROSTO
USE PROSTO
USE O2AVRG
!=======================================================================
USE hfodd_sizes
!=======================================================================
IMPLICIT NONE
!=======================================================================
 LOGICAL::CORROT
 INTEGER::NXMAXX,NYMAXX,NZMAXX,ISIMPY,ISIGNY,ISIQTY,IPAIRI,ICHARG,KETA_R,IROTAT,MREVER
 REAL   ::DROTSQ(0:NDKART),EROTSQ(0:NDKART),PROTSQ(0:NDKART),TROTSQ(0:NDKART),&
          AROTLI(0:NDKART),PROTLI(0:NDKART),PROTKI(0:NDKART)
 COMPLEX::DKOTSQ(0:NDKART),EKOTSQ(0:NDKART),PKOTSQ(0:NDKART),TKOTSQ(0:NDKART),&
          AKOTLI(0:NDKART),PKOTLI(0:NDKART),PKOTKI(0:NDKART)
!
 INTEGER::IALLOC,IBRA,KARTEZ,IMSIGN
 REAL   ::FACTOR
 COMPLEX:: RESULT,C_ZERO,UNIT_I
 REAL   ,ALLOCATABLE:: HAUXIL(:,:)
 COMPLEX,ALLOCATABLE:: HAUXIJ(:,:,:)
 COMPLEX,ALLOCATABLE:: RHO_PP(:,:,:),RHO_PM(:,:,:),KAI_PM(:,:,:),KAI_PP(:,:,:),KAI2PP(:,:,:)
!
 INTEGER::IPHAPP,IPHAPM,IPHAMP,IPHAMM,LDBASE,ITESTA,ITESTO,NOITER,NUMITE
 COMMON /T_PHAS/ IPHAPP(0:NDYMAX,0:NDYMAX,0:NDKART),IPHAPM(0:NDYMAX,0:NDYMAX,0:NDKART),&
                 IPHAMP(0:NDYMAX,0:NDYMAX,0:NDKART),IPHAMM(0:NDYMAX,0:NDYMAX,0:NDKART)
 COMMON /DIMENS/ LDBASE
 COMMON /ITERTS/ ITESTA,ITESTO,NOITER,NUMITE
!=======================================================================
 CALL CPUTIM('ROTAVR',1)
 IALLOC=0
 !
 C_ZERO=CMPLX(0.0D0,0.0D0)
 UNIT_I=CMPLX(0.0D0,1.0D0)
 FACTOR=1.0D0
!=======================================================================
!allocate space for 2[J.rho.J] and 2[J.kAPPa.J] matrix
!=======================================================================
 IF(KETA_R.EQ.2)THEN
    IF(.NOT.ALLOCATED(AROTPP))THEN
       ALLOCATE(AROTPP(NDBASE,NDBASE,1:NDKART,0:NDREVE,0:NDISOS),STAT=IALLOC)
       IF (IALLOC.NE.0) CALL NOALLO('AROTPP','ROTAVR')
    ENDIF
    IF(.NOT.ALLOCATED(PROTPM) .AND. IPAIRI.EQ.1)THEN
       ALLOCATE(PROTPM(NDBASE,NDBASE,1:NDKART,0:NDREVE,0:NDISOS),STAT=IALLOC)
       IF (IALLOC.NE.0) CALL NOALLO('PROTPM','ROTAVR')
    ENDIF
!
    IF(ISIMPY.NE.1)THEN
       IF(.NOT.ALLOCATED(AROTPM))THEN
          ALLOCATE(AROTPM(NDBASE,NDBASE,1:NDKART,0:NDREVE,0:NDISOS),STAT=IALLOC)
          IF (IALLOC.NE.0) CALL NOALLO('AROTPM','ROTAVR')
       ENDIF
       IF(.NOT.ALLOCATED(PROTPP) .AND. IPAIRI.EQ.1)THEN
          ALLOCATE(PROTPP(NDBASE,NDBASE,1:NDKART,0:NDREVE,0:NDISOS),STAT=IALLOC)
          IF (IALLOC.NE.0) CALL NOALLO('PROTPP','ROTAVR')
       ENDIF
    ENDIF
 ENDIF
!=======================================================================
!allocate space for rho and kappa in case they are not allocate yet
!=======================================================================
 IF(.NOT.ALLOCATED(BIG_PP)) THEN
    ALLOCATE(BIG_PP(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('BIG_PP','ROTAVR')
 ENDIF
 IF(.NOT.ALLOCATED(DEN_PP)) THEN
    ALLOCATE(DEN_PP(1:NDBASE,1:NDBASE,0:NDREVE,0:NDISOS),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('DEN_PP','SAVDEN')
 ENDIF
!
 IF(.NOT.ALLOCATED(PAI_PM) .AND. IPAIRI.EQ.1)THEN
    ALLOCATE(PAI_PM(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('PAI_PM','ROTAVR')
 ENDIF
 IF(.NOT.ALLOCATED(CHI_PM) .AND. IPAIRI.EQ.1)THEN
    ALLOCATE(CHI_PM(1:NDBASE,1:NDBASE,0:NDREVE,0:NDISOS), STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('CHI_PM','ROTAVR')
 ENDIF
!
 IF(ISIMPY.NE.1)THEN
    IF(.NOT.ALLOCATED(BIG_PM))THEN
        ALLOCATE (BIG_PM(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
        IF (IALLOC.NE.0) CALL NOALLO('BIG_PM','ROTAVR')
    END IF
    IF(.NOT.ALLOCATED(DEN_PM))THEN
        ALLOCATE (DEN_PM(1:NDBASE,1:NDBASE,0:NDREVE,0:NDISOS),STAT=IALLOC)
        IF (IALLOC.NE.0) CALL NOALLO('DEN_PM','ROTAVR')
    ENDIF
!
    IF (.NOT.ALLOCATED(PAI_PP) .AND. IPAIRI.EQ.1)THEN
       ALLOCATE (PAI_PP(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
       IF (IALLOC.NE.0) CALL NOALLO('PAI_PP','ROTAVR')
    END IF
    IF (.NOT.ALLOCATED(PAI2PP) .AND. IPAIRI.EQ.1)THEN
       ALLOCATE (PAI2PP(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
       IF (IALLOC.NE.0) CALL NOALLO('PAI2PP','ROTAVR')
    END IF
    IF (.NOT.ALLOCATED(CHI_PP) .AND. IPAIRI.EQ.1)THEN
       ALLOCATE (CHI_PP(1:NDBASE,1:NDBASE,0:NDREVE,0:NDISOS),STAT=IALLOC)
       IF (IALLOC.NE.0) CALL NOALLO('CHI_PP','ROTAVR')
    END IF
 ENDIF
!=======================================================================
!allocate work space
!=======================================================================
 ALLOCATE (HAUXIL(1:NDBASE,1:NDBASE),STAT=IALLOC)
 IF (IALLOC.NE.0) CALL NOALLO('HAUXIL','ROTAVR')
 ALLOCATE (HAUXIJ(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
 IF (IALLOC.NE.0) CALL NOALLO('HAUXIJ','ROTAVR')
!
 ALLOCATE (RHO_PP(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
 IF (IALLOC.NE.0) CALL NOALLO('RHO_PP','ROTAVR')
!
 IF(ISIMPY.NE.1)THEN
    ALLOCATE (RHO_PM(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('RHO_PM','ROTAVR')
 ENDIF
!
 IF (IPAIRI.EQ.1) THEN
    ALLOCATE (KAI_PM(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('KAI_PM','ROTAVR')
!
    IF(ISIMPY.NE.1)THEN
       ALLOCATE (KAI_PP(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
       IF (IALLOC.NE.0) CALL NOALLO('KAI_PP','ROTAVR')
       ALLOCATE (KAI2PP(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
       IF (IALLOC.NE.0) CALL NOALLO('KAI2PP','ROTAVR')
    ENDIF
 ENDIF
!=======================================================================
!    ZEROING ARRAYS OF AVERAGE VALUES
!=======================================================================
 DROTSQ(:)=0.0D0
 EROTSQ(:)=0.0D0
 PROTSQ(:)=0.0D0
 TROTSQ(:)=0.0D0
 AROTLI(:)=0.0D0
 PROTLI(:)=0.0D0
 PROTKI(:)=0.0D0
!
 DKOTSQ(:)=C_ZERO
 EKOTSQ(:)=C_ZERO
 PKOTSQ(:)=C_ZERO
 TKOTSQ(:)=C_ZERO
 AKOTLI(:)=C_ZERO
 PKOTLI(:)=C_ZERO
 PKOTKI(:)=C_ZERO
!
 IF(KETA_R.EQ.2)THEN
                    AROTPP(:,:,:,:,ICHARG)=C_ZERO
    IF(ISIMPY.NE.1) AROTPM(:,:,:,:,ICHARG)=C_ZERO
!
    IF(IPAIRI.EQ.1)THEN
                        PROTPM(:,:,:,:,ICHARG)=C_ZERO
       IF (ISIMPY.NE.1) PROTPP(:,:,:,:,ICHARG)=C_ZERO
    ENDIF
 ENDIF
!
                 RHO_PP(:,:,:)=C_ZERO
 IF(ISIMPY.NE.1) RHO_PM(:,:,:)=C_ZERO
 IF (IPAIRI.EQ.1) THEN
    KAI_PM(:,:,:)=C_ZERO
     IF(ISIMPY.NE.1) THEN
        KAI_PP(:,:,:)=C_ZERO
        KAI2PP(:,:,:)=C_ZERO
     ENDIF
 ENDIF
!=======================================================================
!    store rho and kAPPa into matrix
!=======================================================================
 IF (KETA_R.EQ.2) THEN !when calculating matrix element,use mixed rho and kAPPa
    IF(ITESTA.EQ.NUMITE) THEN!initial den, chi at the first iteration
                        DEN_PP(:,:,:,ICHARG)= BIG_PP(:,:,:)
       IF (ISIMPY.NE.1) DEN_PM(:,:,:,ICHARG)= BIG_PM(:,:,:)
       IF (IPAIRI.EQ.1) THEN
                           CHI_PM(:,:,:,ICHARG)=PAI_PM(:,:,:)
          IF (ISIMPY.NE.1) CHI_PP(:,:,:,ICHARG)=PAI_PP(:,:,:)
       ENDIF
    ENDIF
!
                     RHO_PP(:,:,:)=DEN_PP(:,:,:,ICHARG)
    IF (ISIMPY.NE.1) RHO_PM(:,:,:)=DEN_PM(:,:,:,ICHARG)
    IF (IPAIRI.EQ.1) THEN
       KAI_PM(:,:,:)=CHI_PM(:,:,:,ICHARG)
       IF (ISIMPY.NE.1)THEN
          KAI_PP(:,:,:)=CHI_PP(:,:,:,ICHARG)
          KAI2PP(:,:,:)=CHI_PP(:,:,:,ICHARG)
       ENDIF
    ENDIF
!
 ELSE!when there is no correction or inside RENINE, use rho and kAPPa from DENMAC and PAIMAC
                     RHO_PP(:,:,:)=BIG_PP(:,:,:)
    IF (ISIMPY.NE.1) RHO_PM(:,:,:)=BIG_PM(:,:,:)
    IF (IPAIRI.EQ.1) THEN
        KAI_PM(:,:,:)=PAI_PM(:,:,:)
        IF (ISIMPY.NE.1) THEN
           KAI_PP(:,:,:)=PAI_PP(:,:,:)
           KAI2PP(:,:,:)=PAI2PP(:,:,:)
        ENDIF
    ENDIF
 ENDIF
!=======================================================================
!                         X-DIRECTION
!=======================================================================
!
!=======================================================================
![Jx]
!=======================================================================
 CALL INT_LX(NXMAXX,NYMAXX,NZMAXX,IPHAPM(0,0,0),HAUXIL,FACTOR)
 HAUXIJ(1:NDBASE,1:NDBASE,0)= UNIT_I*HAUXIL(1:NDBASE,1:NDBASE)
 HAUXIJ(1:NDBASE,1:NDBASE,1)=-UNIT_I*HAUXIL(1:NDBASE,1:NDBASE)
!
 CALL INT_SX(NXMAXX,NYMAXX,NZMAXX,HAUXIL,FACTOR)
 HAUXIJ(1:NDBASE,1:NDBASE,0)=HAUXIJ(1:NDBASE,1:NDBASE,0)+HAUXIL(1:NDBASE,1:NDBASE)
 HAUXIJ(1:NDBASE,1:NDBASE,1)=HAUXIJ(1:NDBASE,1:NDBASE,1)+HAUXIL(1:NDBASE,1:NDBASE)
!=======================================================================
!Tr[Jx.rho]
!=======================================================================
 IF (ISIMPY.NE.1) THEN
! IF (IROTAT.EQ.1 .AND. ISIMPY.NE.1 .AND. ISIQTY.NE.1) THEN
     CALL AVRXSO(HAUXIJ,AROTLI(1),AKOTLI(1),RHO_PM)
 END IF
!
 IF (CORROT) THEN
    IMSIGN=1
!=======================================================================
!Tr[(Jx.rho)^2] and 2[Jx.rho.Jx]
!=======================================================================

    IF(ISIMPY.EQ.1 .AND. KETA_R.NE.2)then
       CALL EXCXSO(IMSIGN,MREVER,ISIMPY,KETA_R,EROTSQ(1),EKOTSQ(1),HAUXIJ,RHO_PP)
    ELSEIF(ISIMPY.EQ.0 .AND. KETA_R.NE.2)then
       CALL EXCXSO(IMSIGN,MREVER,ISIMPY,KETA_R,EROTSQ(1),EKOTSQ(1),HAUXIJ,RHO_PP,RPM=RHO_PM)
    ELSEIF(ISIMPY.EQ.1 .AND.  KETA_R.EQ.2)then
       CALL EXCXSO(IMSIGN,MREVER,ISIMPY,KETA_R,EROTSQ(1),EKOTSQ(1),HAUXIJ,RHO_PP,APP=AROTPP(:,:,1,:,ICHARG))
    ELSEIF(ISIMPY.EQ.0 .AND.  KETA_R.EQ.2)then
       CALL EXCXSO(IMSIGN,MREVER,ISIMPY,KETA_R,EROTSQ(1),EKOTSQ(1),HAUXIJ,RHO_PP,RPM=RHO_PM,&
                   APP=AROTPP(:,:,1,:,ICHARG),APM=AROTPM(:,:,1,:,ICHARG))
    ELSE
       stop 'wrong ISIMPY or KETA_R in rotavr'
    ENDIF
!=======================================================================
!Tr[(Jx.kAPPa)(Jx.kAPPa')*] and 2[Jx.kAPPa.Jx]
!=======================================================================
    IF (IPAIRI.EQ.1) THEN
       IF(ISIMPY.EQ.1 .AND. KETA_R.NE.2)then
          CALL PAIXSO(MREVER,ISIMPY,KETA_R,PROTSQ(1),PKOTSQ(1),HAUXIJ,KAI_PM)
       ELSEIF(ISIMPY.EQ.0 .AND. KETA_R.NE.2)then
          CALL PAIXSO(MREVER,ISIMPY,KETA_R,PROTSQ(1),PKOTSQ(1),HAUXIJ,KAI_PM,KPP1=KAI_PP,KPP2=KAI2PP)
       ELSEIF(ISIMPY.EQ.1 .AND.  KETA_R.EQ.2)then
          CALL PAIXSO(MREVER,ISIMPY,KETA_R,PROTSQ(1),PKOTSQ(1),HAUXIJ,KAI_PM,PPM=PROTPM(:,:,1,:,ICHARG))
       ELSEIF(ISIMPY.EQ.0 .AND.  KETA_R.EQ.2)then
          CALL PAIXSO(MREVER,ISIMPY,KETA_R,PROTSQ(1),PKOTSQ(1),HAUXIJ,KAI_PM,KPP1=KAI_PP,KPP2=KAI2PP,&
                      PPM=PROTPM(:,:,1,:,ICHARG),PPP=PROTPP(:,:,1,:,ICHARG))
       ELSE
         stop 'wrong ISIMPY or KETA_R in rotavr'
       ENDIF
    END IF

!=======================================================================
![Jx^2]=[Lx^2]+2[Lx.Sx]+[Sx^2]
!=======================================================================
    CALL INTLX2(HAUXIL,FACTOR)
    HAUXIJ(1:NDBASE,1:NDBASE,0)=HAUXIL(1:NDBASE,1:NDBASE)
    HAUXIJ(1:NDBASE,1:NDBASE,1)=HAUXIL(1:NDBASE,1:NDBASE)
!
    CALL INT_LX(NXMAXX,NYMAXX,NZMAXX,IPHAPP(0,0,1),HAUXIL,FACTOR)
    HAUXIJ(1:NDBASE,1:NDBASE,0)=HAUXIJ(1:NDBASE,1:NDBASE,0)+UNIT_I*HAUXIL(1:NDBASE,1:NDBASE)
!
    HAUXIJ(1:NDBASE,1:NDBASE,1)=HAUXIJ(1:NDBASE,1:NDBASE,1)-UNIT_I*HAUXIL(1:NDBASE,1:NDBASE)
!
    DO IBRA=1,LDBASE
       HAUXIJ(IBRA,IBRA,0)=HAUXIJ(IBRA,IBRA,0)+CMPLX(0.25D0,0.0D0)
       HAUXIJ(IBRA,IBRA,1)=HAUXIJ(IBRA,IBRA,1)+CMPLX(0.25D0,0.0D0)
    END DO
!=======================================================================
!Tr[Jx^2.rho]
!=======================================================================
    CALL AVRXSD(HAUXIJ,DROTSQ(1),DKOTSQ(1),RHO_PP)
 END IF
!=======================================================================
!                         Z-DIRECTION
!=======================================================================
!
!=======================================================================
![Jz]
!=======================================================================
 CALL INT_LZ(NXMAXX,NYMAXX,NZMAXX,IPHAPM(0,0,0),HAUXIL,FACTOR)
 HAUXIJ(1:NDBASE,1:NDBASE,0)=UNIT_I*HAUXIL(1:NDBASE,1:NDBASE)
 HAUXIJ(1:NDBASE,1:NDBASE,1)=-UNIT_I*HAUXIL(1:NDBASE,1:NDBASE)
!
 CALL INT_SZ(NXMAXX,NYMAXX,NZMAXX,HAUXIL,FACTOR)
 HAUXIJ(1:NDBASE,1:NDBASE,0)= HAUXIJ(1:NDBASE,1:NDBASE,0)+UNIT_I*HAUXIL(1:NDBASE,1:NDBASE)
 HAUXIJ(1:NDBASE,1:NDBASE,1)= HAUXIJ(1:NDBASE,1:NDBASE,1)-UNIT_I*HAUXIL(1:NDBASE,1:NDBASE)
!=======================================================================
!Tr[Jz.rho]
!=======================================================================
 IF (ISIMPY.NE.1) THEN
! IF (IROTAT.EQ.1 .AND. ISIMPY.NE.1 .AND. ISIQTY.NE.1) THEN
     CALL AVRXSO(HAUXIJ,AROTLI(3),AKOTLI(3),RHO_PM)
 END IF
!
 IF (CORROT) THEN
    IMSIGN=1
!=======================================================================
!Tr[(Jz.rho)^2] and 2[Jz.rho.Jz]
!=======================================================================

    IF(ISIMPY.EQ.1 .AND. KETA_R.NE.2)then
       CALL EXCXSO(IMSIGN,MREVER,ISIMPY,KETA_R,EROTSQ(3),EKOTSQ(3),HAUXIJ,RHO_PP)
    ELSEIF(ISIMPY.EQ.0 .AND. KETA_R.NE.2)then
       CALL EXCXSO(IMSIGN,MREVER,ISIMPY,KETA_R,EROTSQ(3),EKOTSQ(3),HAUXIJ,RHO_PP,RPM=RHO_PM)
    ELSEIF(ISIMPY.EQ.1 .AND.  KETA_R.EQ.2)then
       CALL EXCXSO(IMSIGN,MREVER,ISIMPY,KETA_R,EROTSQ(3),EKOTSQ(3),HAUXIJ,RHO_PP,APP=AROTPP(:,:,3,:,ICHARG))
    ELSEIF(ISIMPY.EQ.0 .AND.  KETA_R.EQ.2)then
       CALL EXCXSO(IMSIGN,MREVER,ISIMPY,KETA_R,EROTSQ(3),EKOTSQ(3),HAUXIJ,RHO_PP,RPM=RHO_PM,&
                   APP=AROTPP(:,:,3,:,ICHARG), APM=AROTPM(:,:,3,:,ICHARG))
    ELSE
       stop 'wrong ISIMPY or KETA_R in rotavr'
    ENDIF
!=======================================================================
!Tr[(Jz.kAPPa)(Jz.kAPPa')*] and 2[Jz.kAPPa.Jz]
!=======================================================================
    IF (IPAIRI.EQ.1) THEN
       IF(ISIMPY.EQ.1 .AND. KETA_R.NE.2)then
          CALL PAIXSO(MREVER,ISIMPY,KETA_R,PROTSQ(3),PKOTSQ(3),HAUXIJ,KAI_PM)
       ELSEIF(ISIMPY.EQ.0 .AND. KETA_R.NE.2)then
          CALL PAIXSO(MREVER,ISIMPY,KETA_R,PROTSQ(3),PKOTSQ(3),HAUXIJ,KAI_PM,KPP1=KAI_PP,KPP2=KAI2PP)
       ELSEIF(ISIMPY.EQ.1 .AND.  KETA_R.EQ.2)then
          CALL PAIXSO(MREVER,ISIMPY,KETA_R,PROTSQ(3),PKOTSQ(3),HAUXIJ,KAI_PM,PPM=PROTPM(:,:,3,:,ICHARG))
       ELSEIF(ISIMPY.EQ.0 .AND.  KETA_R.EQ.2)then
          CALL PAIXSO(MREVER,ISIMPY,KETA_R,PROTSQ(3),PKOTSQ(3),HAUXIJ,KAI_PM,KPP1=KAI_PP,KPP2=KAI2PP,&
                      PPM=PROTPM(:,:,3,:,ICHARG),PPP=PROTPP(:,:,3,:,ICHARG))
       ELSE
         stop 'wrong ISIMPY or KETA_R in rotavr'
       ENDIF
    END IF
!=======================================================================
![Jz^2]
!=======================================================================
    CALL INTLZ2(HAUXIL,FACTOR)
    HAUXIJ(1:NDBASE,1:NDBASE,0)=HAUXIL(1:NDBASE,1:NDBASE)
    HAUXIJ(1:NDBASE,1:NDBASE,1)=HAUXIL(1:NDBASE,1:NDBASE)
!
    CALL INT_LZ(NXMAXX,NYMAXX,NZMAXX,IPHAPP(0,0,3),HAUXIL,FACTOR)
    HAUXIJ(1:NDBASE,1:NDBASE,0)=HAUXIJ(1:NDBASE,1:NDBASE,0)-HAUXIL(1:NDBASE,1:NDBASE)
    HAUXIJ(1:NDBASE,1:NDBASE,1)=HAUXIJ(1:NDBASE,1:NDBASE,1)-HAUXIL(1:NDBASE,1:NDBASE)
!
    DO IBRA=1,LDBASE
       HAUXIJ(IBRA,IBRA,0)=HAUXIJ(IBRA,IBRA,0)+CMPLX(0.25D0,0.0D0)
       HAUXIJ(IBRA,IBRA,1)=HAUXIJ(IBRA,IBRA,1)+CMPLX(0.25D0,0.0D0)
    END DO
!=======================================================================
!Tr[(Jz^2.rho]
!=======================================================================
    CALL AVRXSD(HAUXIJ,DROTSQ(3),DKOTSQ(3),RHO_PP)
 END IF
!=======================================================================
!                         Y-DIRECTION
!=======================================================================
!
!=======================================================================
![Jy]
!=======================================================================
 CALL INT_LY(NXMAXX,NYMAXX,NZMAXX,IPHAPP(0,0,0),HAUXIL,FACTOR)
 HAUXIJ(1:NDBASE,1:NDBASE,0)=UNIT_I*HAUXIL(1:NDBASE,1:NDBASE)
 HAUXIJ(1:NDBASE,1:NDBASE,1)=UNIT_I*HAUXIL(1:NDBASE,1:NDBASE)
!
 CALL INT_SY(NXMAXX,NYMAXX,NZMAXX,HAUXIL,FACTOR)
 HAUXIJ(1:NDBASE,1:NDBASE,0)=HAUXIJ(1:NDBASE,1:NDBASE,0)+HAUXIL(1:NDBASE,1:NDBASE)
 HAUXIJ(1:NDBASE,1:NDBASE,1)=HAUXIJ(1:NDBASE,1:NDBASE,1)-HAUXIL(1:NDBASE,1:NDBASE)
!=======================================================================
!Tr[Jy.rho]
!=======================================================================
!  IF (IROTAT.EQ.1 .AND. ((ISIMPY.EQ.1 .AND. ISIGNY.NE.1) .OR. (ISIMPY.NE.1 .AND. ISIQTY.NE.1))) THEN
 CALL AVRXSD(HAUXIJ,AROTLI(2),AKOTLI(2),RHO_PP)
!  END IF
!
 IF (CORROT) THEN
    IMSIGN=-1
!=======================================================================
!Tr[(Jy.rho)^2] and 2[Jy.rho.Jy]
!=======================================================================
    IF(ISIMPY.EQ.1 .AND. KETA_R.NE.2)then
       CALL EXCXSD(IMSIGN,MREVER,ISIMPY,KETA_R,EROTSQ(2),EKOTSQ(2),HAUXIJ,RHO_PP)
    ELSEIF(ISIMPY.EQ.0 .AND. KETA_R.NE.2)then
       CALL EXCXSD(IMSIGN,MREVER,ISIMPY,KETA_R,EROTSQ(2),EKOTSQ(2),HAUXIJ,RHO_PP,RPM=RHO_PM)
    ELSEIF(ISIMPY.EQ.1 .AND.  KETA_R.EQ.2)then
       CALL EXCXSD(IMSIGN,MREVER,ISIMPY,KETA_R,EROTSQ(2),EKOTSQ(2),HAUXIJ,RHO_PP,APP=AROTPP(:,:,2,:,ICHARG))
    ELSEIF(ISIMPY.EQ.0 .AND.  KETA_R.EQ.2)then
       CALL EXCXSD(IMSIGN,MREVER,ISIMPY,KETA_R,EROTSQ(2),EKOTSQ(2),HAUXIJ,RHO_PP,RPM=RHO_PM,&
                   APP=AROTPP(:,:,2,:,ICHARG), APM=AROTPM(:,:,2,:,ICHARG))
    ELSE
       stop 'wrong ISIMPY or KETA_R in rotavr'
    ENDIF
!=======================================================================
!Tr[(Jy.kAPPa)(Jy.kAPPa')*] and 2[Jy.kAPPa.Jy]
!=======================================================================
    IF (IPAIRI.EQ.1) THEN
       IF(ISIMPY.EQ.1 .AND. KETA_R.NE.2)then
          CALL PAIXSD(MREVER,ISIMPY,KETA_R,PROTSQ(2),PKOTSQ(2),HAUXIJ,KAI_PM)
       ELSEIF(ISIMPY.EQ.0 .AND. KETA_R.NE.2)then
          CALL PAIXSD(MREVER,ISIMPY,KETA_R,PROTSQ(2),PKOTSQ(2),HAUXIJ,KAI_PM,KPP1=KAI_PP,KPP2=KAI2PP)
       ELSEIF(ISIMPY.EQ.1 .AND.  KETA_R.EQ.2)then
          CALL PAIXSD(MREVER,ISIMPY,KETA_R,PROTSQ(2),PKOTSQ(2), HAUXIJ,KAI_PM,PPM=PROTPM(:,:,2,:,ICHARG))
       ELSEIF(ISIMPY.EQ.0 .AND.  KETA_R.EQ.2)then
          CALL PAIXSD(MREVER,ISIMPY,KETA_R,PROTSQ(2),PKOTSQ(2),HAUXIJ,KAI_PM,KPP1=KAI_PP,KPP2=KAI2PP,&
                      PPM=PROTPM(:,:,2,:,ICHARG),PPP=PROTPP(:,:,2,:,ICHARG))
       ELSE
         stop 'wrong ISIMPY or KETA_R in rotavr'
       ENDIF
    END IF
!=======================================================================
![Jy^2]
!=======================================================================
    CALL INTLY2(HAUXIL,FACTOR)
    HAUXIJ(1:NDBASE,1:NDBASE,0)=HAUXIL(1:NDBASE,1:NDBASE)
    HAUXIJ(1:NDBASE,1:NDBASE,1)=HAUXIL(1:NDBASE,1:NDBASE)
!
    CALL INT_LY(NXMAXX,NYMAXX,NZMAXX,IPHAPP(0,0,2),HAUXIL,FACTOR)
    HAUXIJ(1:NDBASE,1:NDBASE,0)=HAUXIJ(1:NDBASE,1:NDBASE,0)+UNIT_I*HAUXIL(1:NDBASE,1:NDBASE)

    HAUXIJ(1:NDBASE,1:NDBASE,1)=HAUXIJ(1:NDBASE,1:NDBASE,1)-UNIT_I*HAUXIL(1:NDBASE,1:NDBASE)
!
    DO IBRA=1,LDBASE
       HAUXIJ(IBRA,IBRA,0)=HAUXIJ(IBRA,IBRA,0)+CMPLX(0.25D0,0.0D0)
       HAUXIJ(IBRA,IBRA,1)=HAUXIJ(IBRA,IBRA,1)+CMPLX(0.25D0,0.0D0)
    END DO
!=======================================================================
!Tr[(Jy^2.rho]
!=======================================================================
    CALL AVRXSD(HAUXIJ,DROTSQ(2),DKOTSQ(2),RHO_PP)
 END IF
!
!=======================================================================
!    HERE ALL TERMS ARE COMBINED INTO THE AVERAGE OF THE SQUARE OF
!    THE TOTAL ANGULAR MOMENTUM.
!=======================================================================
 DO KARTEZ=1,NDKART
    TROTSQ(KARTEZ)=DROTSQ(KARTEZ)-EROTSQ(KARTEZ)&
                  +PROTSQ(KARTEZ)&
                  +AROTLI(KARTEZ)**2
    TKOTSQ(KARTEZ)=DKOTSQ(KARTEZ)-EKOTSQ(KARTEZ)&
                  +PKOTSQ(KARTEZ)&
                  +AKOTLI(KARTEZ)**2
 END DO
! write(*,*)'drs',drotsq
! write(*,*)'ers',EROTSQ
! write(*,*)'prs',protsq
! write(*,*)'trs',trotsq
!  write(*,*)'a',akotli

 DO KARTEZ=1,NDKART
    DROTSQ(0)=DROTSQ(0)+DROTSQ(KARTEZ)
    EROTSQ(0)=EROTSQ(0)+EROTSQ(KARTEZ)
    PROTSQ(0)=PROTSQ(0)+PROTSQ(KARTEZ)
    TROTSQ(0)=TROTSQ(0)+TROTSQ(KARTEZ)
    AROTLI(0)=AROTLI(0)+AROTLI(KARTEZ)**2
    PROTLI(0)=PROTLI(0)+PROTLI(KARTEZ)**2
    PROTKI(0)=PROTKI(0)+PROTKI(KARTEZ)**2
!
    DKOTSQ(0)=DKOTSQ(0)+DKOTSQ(KARTEZ)
    EKOTSQ(0)=EKOTSQ(0)+EKOTSQ(KARTEZ)
    PKOTSQ(0)=PKOTSQ(0)+PKOTSQ(KARTEZ)
    TKOTSQ(0)=TKOTSQ(0)+TKOTSQ(KARTEZ)
    AKOTLI(0)=AKOTLI(0)+AKOTLI(KARTEZ)**2
    PKOTLI(0)=PKOTLI(0)+PKOTLI(KARTEZ)**2
    PKOTKI(0)=PKOTKI(0)+PKOTKI(KARTEZ)**2
 END DO
!=======================================================================
!Deallocate work space
!=======================================================================
 DEALLOCATE (HAUXIL)
 DEALLOCATE (HAUXIJ)
 DEALLOCATE (RHO_PP)
 IF(ISIMPY.NE.1) DEALLOCATE (RHO_PM)
 IF(IPAIRI.EQ.1) DEALLOCATE (KAI_PM)
 IF(IPAIRI.EQ.1 .AND. ISIMPY.NE.1) DEALLOCATE (KAI_PP,KAI2PP)
!=======================================================================
 CALL CPUTIM('ROTAVR',0)
 RETURN
END SUBROUTINE ROTAVR
!=======================================================================
SUBROUTINE RENINE(MIN_QP,IRENIN,IDOGOA,&
                  NXHERM,NYHERM,NZHERM,NXMAXX,NYMAXX,NZMAXX,&
                  KSIMTX,JSIMTY,KSIMTZ,KSIGNY,KSIMPY,KSIQTY,&
                  IPAIRI,IPAHFB,IROTAT,ITIREP,KREVER,IPNMIX,&
                  IDEVAR,ITERUN,NMUCOU,ISHIFY,NUMCOU,BOUCOU,&
                  ISYMDE,INIROT,INIINV,INIKAR,ISAWAV,IKERNE,&
                  NOSCIL,NUMETA,NFICOU,FILCOU,ICOULI,ICOULO,&
                                IN_FIX,IZ_FIX,ICANRP,SLOWRP)
!=======================================================================
 USE SAVLEF
 USE SAVRIG
 USE PD_DEN
 USE WAVR_L
 USE HCOULO
 USE OFFDIA, ONLY:RMATCN,ROTWLR,VRVMAT,H_KERN
!=======================================================================
 USE hfodd_sizes
!=======================================================================
IMPLICIT NONE
!=======================================================================
 CHARACTER::FILCOU*68
 INTEGER::MIN_QP,IRENIN,IDOGOA,&
          NXHERM,NYHERM,NZHERM,NXMAXX,NYMAXX,NZMAXX,&
          KSIMTX,JSIMTY,KSIMTZ,KSIGNY,KSIMPY,KSIQTY,&
          IPAIRI,IPAHFB,IROTAT,ITIREP,KREVER,IPNMIX,&
          IDEVAR,ITERUN,NMUCOU,ISHIFY,NUMCOU,ISYMDE,&
          INIROT,INIINV,INIKAR,ISAWAV,IKERNE,&
          NOSCIL,NUMETA,NFICOU,ICOULI,ICOULO,&
          IN_FIX,IZ_FIX,ICANRP,I_SLOW
 REAL   ::BOUCOU,SLOWRP
!=======================================================================
 CHARACTER::NAMEPN*8,FCOUOF*68,CLAB*10
 LOGICAL::PRINIT,CORROT
 INTEGER::IALLOC,JMULMO,JDKART,ICANON,IEXTEN,IFULLD,INROTA,ICHARG,IDKART,KETA_R,&
          LTIMUP,LUPPER,LTOTAL,LDPNMX,MREVER,ISIMPY,ISIGNY,ISIQTY,INPCOU,I_SUCC,&
          ITRUNC,IDIGTA,IMASKN,IREINE,ISKYRM,ISIMTX,ISIMTZ,NLAB
 REAL   ::TSHIFT,H2PARA,A_PARA,REINER,GOAMAS,BAREKE
 COMPLEX::C_ZERO,C_UNIT,EKESKY,BAREKK
 INTEGER,DIMENSION(0:NDISOS)::LSTEMP,LUTEMP,LTTEMP,LDTEMP
 REAL   ,DIMENSION(0:NDKART)::DROTSQ,EROTSQ,PROTSQ,TROTSQ,AROTLI,PROTLI,PROTKI
 COMPLEX,DIMENSION(0:NDSHIF,1:NDKART,0:NDISOS)::AKKERN
 COMPLEX,DIMENSION(0:NDSHIF,0:NDISOS)::OVKERN,EKKERN,EPKERN,DKKERN
 COMPLEX,DIMENSION(0:NDSHIF)::SKKERN,CDKERN,CXKERN,ETKERN,JTKERN
 COMPLEX,DIMENSION(0:NDKART)::DKOTSQ,EKOTSQ,PKOTSQ,TKOTSQ,AKOTLI,PKOTLI,PKOTKI
 COMPLEX,ALLOCATABLE::COUHPP(:,:),COUHPM(:,:)
!=======================================================================
 INTEGER::LDBASE,LDTOTS,LDSTAT,LDUPPE,LDTIMU,NFIPRI,ITESTA,ITESTO,NOITER,NUMITE
 REAL   ::ALPROT,BETROT,GAMROT,ROTATX,ROTATY,ROTATZ,HBMASS,HBMRPA,HBMINP,ROTREN,&
          GOAREN,GOAREI,ROTRIN
!
 REAL   ::COMULT
 COMPLEX::QMUT_N,QMUT_P,QMUT_T,QMUL_N,QMUL_P,QMUL_T
 REAL   ,DIMENSION(0:NDMULT,0:NDOSCI,0:NDOSCI,1:NDKART)::COM_TMP
 COMPLEX,DIMENSION(0:NDMULT,-NDMULT:NDMULT)::QTN_TMP,QTP_TMP,QTT_TMP,QLN_TMP,QLP_TMP,QLT_TMP
!
 COMMON /QMUTTI/ QMUT_N(0:NDMULT,-NDMULT:NDMULT),QMUT_P(0:NDMULT,-NDMULT:NDMULT),QMUT_T(0:NDMULT,-NDMULT:NDMULT)
 COMMON /QMULTI/ QMUL_N(0:NDMULT,-NDMULT:NDMULT),QMUL_P(0:NDMULT,-NDMULT:NDMULT),QMUL_T(0:NDMULT,-NDMULT:NDMULT)
 COMMON /COEMUL/ COMULT(0:NDMULT,0:NDOSCI,0:NDOSCI,1:NDKART)
!
! REAL   ::COEYLM,QUNITS
! INTEGER::NREYLM,KREYLM,NIMYLM,KIMYLM
! REAL   ,DIMENSION(0:NDMULT,0:NDMULT)::COETMP,QUNTMP
! INTEGER,DIMENSION(0:NDMULT,0:NDMULT)::NRETMP
! INTEGER,DIMENSION(0:NDMULT,0:NDMULT,1:NDTERM,0:NDKART)::KRETMP
! INTEGER,DIMENSION(0:NDMULT,1:NDMULT,1:NDTERM,0:NDKART)::KIMTMP
! INTEGER,DIMENSION(0:NDMULT,1:NDMULT)::NIMTMP
! COMMON /DATYLM/ COEYLM(0:NDMULT,0:NDMULT),NREYLM(0:NDMULT,0:NDMULT),KREYLM(0:NDMULT,0:NDMULT,1:NDTERM,0:NDKART),&
!                 NIMYLM(0:NDMULT,1:NDMULT),KIMYLM(0:NDMULT,1:NDMULT,1:NDTERM,0:NDKART)
! COMMON /OURUNI/ QUNITS(0:NDMULT,0:NDMULT)
!
 COMMON /RANGLE/ ALPROT,BETROT,GAMROT
 COMMON /DIMENS/ LDBASE
 COMMON /DIMSTA/ LDTOTS(0:NDISOS),LDSTAT(0:NDISOS),LDUPPE(0:NDISOS),LDTIMU(0:NDISOS)
 COMMON /ROTATT/ ROTATX,ROTATY,ROTATZ
 COMMON /PLANCK/ HBMASS,HBMRPA,HBMINP
 COMMON /RENROT/ ROTREN(NDKART)
 COMMON /RENGOA/ GOAREN(NDKART),GOAREI(NDKART)
 COMMON /CFIPRI/ NFIPRI
 COMMON /ITERTS/ ITESTA,ITESTO,NOITER,NUMITE
 COMMON /RINROT/ ROTRIN(NDKART)
!=======================================================================
 CALL CPUTIM('RENINE',1)
!=======================================================================
 IALLOC=0
 I_SLOW=0
 C_ZERO=CMPLX(0.0D0,0.0D0)
 C_UNIT=CMPLX(1.0D0,0.0D0)
!=======================================================================
 IF (.NOT.ALLOCATED(WARIGH)) THEN
     ALLOCATE (WARIGH(1:NDBASE,1:4*NDSTAT,0:NDSPIN),STAT=IALLOC)
     IF (IALLOC.NE.0) CALL NOALLO('WARIGH','RENINE')
     WARIGH(:,:,:)=C_ZERO
 END IF
 IF (.NOT.ALLOCATED(WALEFT)) THEN
     ALLOCATE (WALEFT(1:NDBASE,1:4*NDSTAT,0:NDSPIN),STAT=IALLOC)
     IF (IALLOC.NE.0) CALL NOALLO('WALEFT','RENINE')
     WALEFT(:,:,:)=C_ZERO
 END IF
!=======================================================================
 ALLOCATE (COUHPP(1:NDBASE,1:NDBASE),STAT=IALLOC)
 IF (IALLOC.NE.0) CALL NOALLO('COUHPP','RENINE')
 ALLOCATE (COUHPM(1:NDBASE,1:NDBASE),STAT=IALLOC)
 IF (IALLOC.NE.0) CALL NOALLO('COUHPM','RENINE')
!=======================================================================
!
!=======================================================================
!         THIS SUBROUTINE CALCULATES THE INERTIA RENORMALIZATION FACTORS
!=======================================================================
 NAMEPN='DUMMY   '
 PRINIT=.FALSE.
 JMULMO=NMUCOU
 CORROT=.TRUE.
 KETA_R=0
!=======================================================================
!         STORING THE COULOMB MATRIX ELEMENTS AND SO ON.
!         THESE WILL BE RETURNED AFTER THE INERTIA RENORMALIZATION.
!=======================================================================
                  COUHPP(:,:)=HPPCOU(:,:)
 IF (KSIMPY.NE.1) COUHPM(:,:)=HPMCOU(:,:)
!
 COM_TMP=COMULT
 QTN_TMP=QMUT_N
 QTP_TMP=QMUT_P
 QTT_TMP=QMUT_T
 QLN_TMP=QMUL_N
 QLP_TMP=QMUL_P
 QLT_TMP=QMUL_T
! COETMP=COEYLM
! NRETMP=NREYLM
! KRETMP=KREYLM
! NIMTMP=NIMYLM
! KIMTMP=KIMYLM
! QUNTMP=QUNITS

 LSTEMP(:)=LDTOTS(:)
 LUTEMP(:)=LDUPPE(:)
 LTTEMP(:)=LDTIMU(:)
 LDTEMP(:)=LDSTAT(:)
!go to 1
!=======================================================================
!STORE SYMMETRY PARAMETERS,
!THESE VALUES MAY BE CHANGED DEPENDING ON ROTATION AXIS
!=======================================================================
 ISIMPY=KSIMPY
 ISIGNY=KSIGNY
 ISIQTY=KSIQTY
 ISIMTX=KSIMTX
 ISIMTZ=KSIMTZ
 MREVER=KREVER
!
 JDKART = 0
 ICANON = ABS(ICANRP)!calculate under canonical basis
 IFULLD = 0
 IEXTEN = 0!extend mrever=0 wave function to mrever=1
 IF(ICANRP .EQ. -1) IFULLD=1!use full canonical basis when ICANRP=-1
 IF(ABS(ROTATX).GT. 1.0D-8)THEN!WHEN ROTATING ALONG X OR Z AXIX, SYMEMTRY WOULD BE BROKEN
    ISIMPY = 0
    ISIGNY = 0
    ISIMTZ = 0
!
    IF(KREVER.EQ.0)THEN!After rotate along x or z axis ,wave function would have off diagonal part, both simplex block is need.
       IEXTEN=1
    ENDIF
 ENDIF
 IF(ABS(ROTATZ).GT. 1.0D-8)THEN!WHEN ROTATING ALONG X OR Z AXIX, SYMEMTRY WOULD BE BROKEN
    ISIMPY = 0
    ISIGNY = 0
    ISIMTX = 0
!
    IF(KREVER.EQ.0)THEN!After rotate along x or z axis ,wave function would have off diagonal part, both simplex block is need.
       IEXTEN=1
    ENDIF
 ENDIF
 IF(IEXTEN.EQ.1) MREVER=1!use both simpliex block
 IF(IPAHFB.NE.1) ICANON=0!in HF and HFB+canonical case one dosen't need to transform into canonical basis
!=======================================================================
!CHECK INPUT PARAMTER
!=======================================================================
 IF (KREVER.EQ.0 .AND. ICANON.EQ.0 .AND.(ABS(ROTATX) .GT. 1.0d-8 .OR. ABS(ROTATZ) .GT. 1.0d-8))&
    STOP ' ROTATION ALONG X OR Z WITH MREVER=0 IN RENINE'
 IF ((ABS(ROTATX).GT. 1.0d-8).AND.(ABS(ROTATY).GT. 1.0d-8).OR.&
     (ABS(ROTATY).GT. 1.0d-8).AND.(ABS(ROTATZ).GT. 1.0d-8).OR.&
     (ABS(ROTATZ).GT. 1.0d-8).AND.(ABS(ROTATX).GT. 1.0d-8))&
    STOP ' ROTATION ALONG MULTIPLE DIRECTION NOT IMPLEMENTED YET in RENINE'
 IF ((ABS(ROTATX).LE. 1.0d-8).AND.(ABS(ROTATY).LE. 1.0d-8).AND.(ABS(ROTATZ).LE. 1.0d-8))&
    STOP ' NO ROTATION ANGLE DETECTED'
! IF (ICANON.GT.0 .AND. KSIMPY.NE.1) STOP 'ICANON>0 WITH ISIMPY =0 NOT IMPLINENTED YET in RENINE'
! IF (ICANON.EQ.0 .AND. IFULLD.EQ.1) STOP 'ONLY FULL CANONICAL BASIS ARE STORED'
 IF (IRENIN.GT.NDSHIF) STOP ' IRENIN.GT.NDSHIF IN RENINE'
! IF (IEXTEN.EQ.1 .AND. (KREVER.EQ.1 .OR. ICANON.EQ.0)) STOP 'IEXTEN=1 needs KREVER=0 and ICANON>0 in RENINE'
!=======================================================================
!    TREAT COULOMB PART WHEN SYMMETRY IS BROKEN INSIDE THIS ROUTINE
!=======================================================================
 IF(ISIMPY.NE.KSIMPY .OR. ISIGNY.NE.KSIGNY .OR. ISIQTY.NE.KSIQTY )THEN
    CALL MULCOU(ISIMPY,ISIGNY,ISIQTY,NMUCOU)
    FCOUOF=TRIM(FILCOU)//'OFD'
    IF(ITESTA.EQ.NUMITE)THEN!DO CALCULATION ONLY AT THE FIRST ITERATION
      ICOULI=0
      ICOULO=1
      CALL PRECOU(NOSCIL,NUMCOU,NUMETA,BOUCOU,NFICOU,FCOUOF,ICOULI,ICOULO)
    ELSE!READ FROM FILE
      INPCOU=1
      I_SUCC=1
      OPEN(NFICOU,FILE=FCOUOF,STATUS='OLD',FORM='UNFORMATTED')
      CALL RECOUL(NFICOU,INPCOU,I_SUCC,NUMCOU,NUMETA,BOUCOU)
      CLOSE(NFICOU)
    ENDIF
 ENDIF
!=======================================================================
!         HERE BEGIN THE LOOP THE ROTATION ANGLES.
!=======================================================================
 DO INROTA=IRENIN,0,-1
    IF (ABS(ROTATY) .GT. 1.0d-8) THEN
     ALPROT=0.0d0
     BETROT=INROTA*ROTATY
     GAMROT=0.0d0
    ELSEIF (ABS(ROTATX) .GT. 1.0d-8) THEN
     ALPROT=ATAN(1.0d0)*2.0d0
     BETROT=INROTA*ROTATX
     GAMROT=-ALPROT
    ELSE
     ALPROT=0.0d0
     BETROT=0.0d0
     GAMROT=INROTA*ROTATZ
    ENDIF
!
    DO ICHARG=0,NDISOS
       IF(IFULLD.EQ.1)THEN!USE FULL CANONICAL BASIS
          LDTIMU(ICHARG)=LDBASE
          LDUPPE(ICHARG)=LDBASE*(1+MREVER)
          LDSTAT(ICHARG)=LDBASE*(2+MREVER)
          LDTOTS(ICHARG)=LDBASE*2*(1+MREVER)
       ENDIF
!
       IF(IEXTEN.EQ.1)THEN!modify orbits number countings
          LDUPPE(ICHARG)=LDTIMU(ICHARG)*2
          LDSTAT(ICHARG)=LDTIMU(ICHARG)*3
          LDTOTS(ICHARG)=LDTIMU(ICHARG)*4
       ENDIF
!
       LTOTAL=LDTOTS(ICHARG)
       LUPPER=LDUPPE(ICHARG)
       LTIMUP=LDTIMU(ICHARG)
       IF (LTOTAL.GT.4*NDSTAT)STOP ' LTOTAL.GT.4*NDSTAT IN RENINE'
!=======================================================================
!FETCHING THE LEFT WAVE FUNCTIONS FOR THE GIVEN CHARGE
!AND CALCULATE ROTATED RIGHT WAVE FUNCTIONS MULTIPLYING THE INVERTED OVERLAP MATRIX
!USE 'KSIMPY' ,since symmetry in wave function is broken after this routine is called
!USE 'mrever' since both blocks are already stored, by change lduppe and so on, the extension is done
!=======================================================================
       CALL ROTWLR(MREVER,KSIMPY,ICANON,IDOGOA,IPAIRI,INROTA,IEXTEN,&
                   LDBASE,LTIMUP,LUPPER,LTOTAL,ICHARG,OVKERN(INROTA,ICHARG))
!=======================================================================
!        HERE BEGIN CALCULATIONS OF DENSITIES AND OBSERVABLES
!=======================================================================
!=======================================================================
!        RHO MATRIX
!=======================================================================
       CALL DENMAC(MREVER,ICHARG,ISIMPY,IPAHFB,WALEFT,WARIGH)
!=======================================================================
!        KAPPA MATRIX
!=======================================================================
       IF (IDOGOA.EQ.0) THEN
          IF (IPAIRI.EQ.1 .AND. CORROT) CALL PAIMAC(MREVER,ICHARG,ISIMPY,WALEFT,WARIGH,IKERNE)
!=======================================================================
!        TERMS FOR <|J^2|BETA>
!=======================================================================
          CALL ROTAVR(NXMAXX,NYMAXX,NZMAXX,ISIMPY,ISIGNY,ISIQTY,&
                      IPAIRI,ICHARG,KETA_R,IROTAT,MREVER,CORROT,&
                      DROTSQ,EROTSQ,PROTSQ,TROTSQ,AROTLI,PROTLI,PROTKI,&
                      DKOTSQ,EKOTSQ,PKOTSQ,TKOTSQ,AKOTLI,PKOTLI,PKOTKI)
!
          DKKERN(INROTA,ICHARG)=TKOTSQ(JDKART)
          AKKERN(INROTA,1,ICHARG)=AKOTLI(1)
          AKKERN(INROTA,2,ICHARG)=AKOTLI(2)
          AKKERN(INROTA,3,ICHARG)=AKOTLI(3)
      ELSE
          DKKERN(INROTA,ICHARG)=C_ZERO
          AKKERN(INROTA,1,ICHARG)=C_ZERO
          AKKERN(INROTA,2,ICHARG)=C_ZERO
          AKKERN(INROTA,3,ICHARG)=C_ZERO
      END IF
!=======================================================================
!        <|H|BETA>
!=======================================================================
!goto 1
      ISKYRM=ICHARG
      CALL H_KERN(NXHERM,NYHERM,NZHERM,NXMAXX,NYMAXX,NZMAXX,&
                  ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY,&
                  IPAHFB,MREVER,ICHARG,MIN_QP,IPNMIX,ITIREP,NAMEPN,&
                  PRINIT,IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,&
                  INIKAR,ISAWAV,IKERNE,ISHIFY,JMULMO,NUMCOU,&
                  BOUCOU,IN_FIX,IZ_FIX,ISKYRM,&
                  EKKERN(INROTA,ICHARG),EPKERN(INROTA,ICHARG),&
                  BAREKE,BAREKK,CDKERN(INROTA),CXKERN(INROTA),SKKERN(INROTA))
!goto 1
    END DO!END ICHARG=0,NDISOS
!goto 1
!=======================================================================
!        SKYREM PART OF <|H|Ri>
!=======================================================================
!       CALL TRUTOD(NXHERM,NYHERM,NZHERM)
!       CALL ESKYRM(NXHERM,NYHERM,NZHERM,ENESKY,ENEVEN,ENEODD,&
!                   ENREAR,ENE_W0,EEVEW0,EODDW0,EKESKY,LDPNMX)
!       SKKERN(INSHIF)=EKESKY
!=======================================================================
!        <|H|BETA> and <|J^2|BETA>
!=======================================================================
    ETKERN(INROTA)=EKKERN(INROTA,0)+EKKERN(INROTA,1)&
                  +EPKERN(INROTA,0)+EPKERN(INROTA,1)&
                  +CDKERN(INROTA)+CXKERN(INROTA)+SKKERN(INROTA)
    JTKERN(INROTA)=DKKERN(INROTA,0) +DKKERN(INROTA,1)&
                +2*AKKERN(INROTA,1,0)*AKKERN(INROTA,1,1)&
                +2*AKKERN(INROTA,2,0)*AKKERN(INROTA,2,1)
!
! write(*,*)abs(etkern(inrota)),abs(jtkern(inrota)),abs(OVKERN(INROTA,0)*OVKERN(INROTA,1)),inrota
! NLAB=INT((ROTATX+ROTATY+ROTATZ)/3.14*50)
! WRITE(CLAB,'(I3)')NLAB
! OPEN(1,FILE='pf'//trim(adjustl(clab))//'.dat')
! write(1,*)etkern(inrota)-(EPKERN(INROTA,0)+EPKERN(INROTA,1)),&
!           EPKERN(INROTA,0)+EPKERN(INROTA,1),OVKERN(INROTA,0)*OVKERN(INROTA,1)
! close(1)
! stop
 END DO!END INROTA=IRENIN,0
!=======================================================================
!        BELOW WE DEFINE THE LIPKIN PARAMETERS, WHICH ARE POSITIVE
!=======================================================================
 IF (IDOGOA.EQ.0) THEN
    REINER=(ETKERN(IRENIN)-ETKERN(0))/(JTKERN(IRENIN)-JTKERN(0))
 ELSE
    TSHIFT=IRENIN**2*(ROTATX**2+ROTATY**2+ROTATZ**2)
!
    H2PARA=-2*   (ETKERN(IRENIN)  -ETKERN(0)       )/TSHIFT
    A_PARA=-2*LOG(OVKERN(IRENIN,0)*OVKERN(IRENIN,1))/TSHIFT
!
    GOAMAS=H2PARA/(2*A_PARA**2)
    IF (MREVER.EQ.0) GOAMAS=GOAMAS/4
! MREVER.EQ.0 with shifting or rotating is not allowed yet
! But this option is still kept inside RENMAS and RENINE.
    REINER=GOAMAS
 END IF
!=======================================================================
!   Only keep a given digital of lipkin parameter to accelate converge
!=======================================================================
 ITRUNC=1
 IDIGTA=8
 IF(ITRUNC.EQ.1)THEN
    IMASKN=10**IDIGTA
    IREINE=REINER*IMASKN
    REINER=IREINE/(IMASKN*1.0d0)
 ENDIF
! IF(REINER.GE.0.1d0) REINER=5.0d-2
! IF(REINER.LE.0.0d0) REINER=0.0d0
!=======================================================================
!        BELOW WE DEFINE THE COEFFICIENTS OF <J^2>, WHICH ARE NEGATIVE
!=======================================================================
 DO IDKART = 1,NDKART
    IF (ROTRIN(IDKART).LE. 0.0d0)THEN
        ROTREN(IDKART)=-ROTRIN(IDKART)
    ELSE
        IF (I_SLOW.NE.1) THEN
            ROTREN(IDKART)=ROTREN(IDKART)*SLOWRP-REINER*(1-SLOWRP)
!           ROTREN(IDKART)=-ROTRIN(IDKART)
        ELSE
            ROTREN(IDKART)=-REINER
        ENDIF
    ENDIF
 ENDDO
!
 GOAREI(1)=-GOAMAS
 GOAREI(2)=-GOAMAS
 GOAREI(3)=-GOAMAS
!=======================================================================
! Rstore quantities temporarily changed inside the routine
!=======================================================================
 LDTOTS(:)=LSTEMP(:)
 LDUPPE(:)=LUTEMP(:)
 LDTIMU(:)=LTTEMP(:)
 LDSTAT(:)=LDTEMP(:)
!
!goto 1
 IF(ISIMPY.NE.KSIMPY .OR. ISIGNY.NE.KSIGNY .OR. ISIQTY.NE.KSIQTY )THEN
    CALL MULCOU(KSIMPY,KSIGNY,KSIQTY,NMUCOU)
    OPEN(NFICOU,FILE=FILCOU,STATUS='OLD',FORM='UNFORMATTED')
    INPCOU=1
    I_SUCC=1
    CALL RECOUL(NFICOU,INPCOU,I_SUCC,NUMCOU,NUMETA,BOUCOU)
    CLOSE(NFICOU)
 ENDIF
! CALL MOMETS(ISIMPY,ISIGNY,ISIQTY,QMUL_P,QMUT_P,COMULT,JMULMO,ISHIFY)
                 HPPCOU(:,:)=COUHPP(:,:)
 IF (ISIMPY.NE.1)HPMCOU(:,:)=COUHPM(:,:)
! COEYLM=COETMP
! NREYLM=NRETMP
! KREYLM=KRETMP
! NIMYLM=NIMTMP
! KIMYLM=KIMTMP
! QUNITS=QUNTMP
 COMULT=COM_TMP
 QMUT_N=QTN_TMP
 QMUT_P=QTP_TMP
 QMUT_T=QTT_TMP
 QMUL_N=QLN_TMP
 QMUL_P=QLP_TMP
 QMUL_T=QLT_TMP
1 continue
!=======================================================================
! WRITE(*,*)'REINER',REINER,GOAMAS,IPAHFB
!=======================================================================
 CALL CPUTIM('RENINE',0)
 RETURN
 END SUBROUTINE RENINE
!=======================================================================
SUBROUTINE LINAVR(NXMAXX,NYMAXX,NZMAXX,ISIMPY,ISIGNY,ISIQTY,&
                  IPAIRI,ICHARG,JETACM,IROTAT,MREVER,COR_CM,&
                  DLINSQ,ELINSQ,PLINSQ,TLINSQ,ALINLI,PLINLI,PLINKI,&
                  DKINSQ,EKINSQ,PKINSQ,TKINSQ,AKINLI,PKINLI,PKINKI)
!=======================================================================
USE MAT_PP
USE MAT_PM
USE MAD_PP
USE MAD_PM
USE KAP_PM
USE KAP_PP
USE KAP2PP
USE MAP_PP
USE MAP_PM
USE ALISTO
USE PLISTO
USE O2AVRG
! USE OSQUARE
!=======================================================================
USE hfodd_sizes
!=======================================================================
IMPLICIT NONE
!=======================================================================
 LOGICAL::COR_CM
 INTEGER::NXMAXX,NYMAXX,NZMAXX,ISIMPY,ISIGNY,ISIQTY,IPAIRI,ICHARG,JETACM,IROTAT,MREVER
 REAL   ::DLINSQ(0:NDKART),ELINSQ(0:NDKART),PLINSQ(0:NDKART),TLINSQ(0:NDKART),&
          ALINLI(0:NDKART),PLINLI(0:NDKART),PLINKI(0:NDKART)
 COMPLEX::DKINSQ(0:NDKART),EKINSQ(0:NDKART),PKINSQ(0:NDKART),TKINSQ(0:NDKART),&
          AKINLI(0:NDKART),PKINLI(0:NDKART),PKINKI(0:NDKART)
!
 INTEGER::IALLOC,IBRA,KARTEZ,IMSIGN
 REAL   ::FACTOR
 COMPLEX:: RESULT,C_ZERO,UNIT_I
 REAL   ,ALLOCATABLE:: HAUXIL(:,:)
 COMPLEX,ALLOCATABLE:: HAUXIJ(:,:,:)
 COMPLEX,ALLOCATABLE:: RHO_PP(:,:,:),RHO_PM(:,:,:),KAI_PM(:,:,:),KAI_PP(:,:,:),KAI2PP(:,:,:)
!
 INTEGER::IPHAPP,IPHAPM,IPHAMP,IPHAMM,LDBASE,ITESTA,ITESTO,NOITER,NUMITE
 COMMON /T_PHAS/ IPHAPP(0:NDYMAX,0:NDYMAX,0:NDKART),IPHAPM(0:NDYMAX,0:NDYMAX,0:NDKART),&
                 IPHAMP(0:NDYMAX,0:NDYMAX,0:NDKART),IPHAMM(0:NDYMAX,0:NDYMAX,0:NDKART)
 COMMON /DIMENS/ LDBASE
 COMMON /ITERTS/ ITESTA,ITESTO,NOITER,NUMITE
!=======================================================================
 CALL CPUTIM('LINAVR',1)
 IALLOC=0
 !
 C_ZERO=CMPLX(0.0D0,0.0D0)
 UNIT_I=CMPLX(0.0D0,1.0D0)
 FACTOR=1.0D0
!=======================================================================
!allocate space for 2[J.rho.J] and 2[J.kappa.J] matrix
!=======================================================================
 IF(JETACM.EQ.2)THEN
    IF(.NOT.ALLOCATED(ALINPP))THEN
       ALLOCATE(ALINPP(NDBASE,NDBASE,1:NDKART,0:NDREVE,0:NDISOS),STAT=IALLOC)
       IF (IALLOC.NE.0) CALL NOALLO('ALINPP','LINAVR')
    ENDIF
    IF(.NOT.ALLOCATED(PLINPM) .AND. IPAIRI.EQ.1)THEN
       ALLOCATE(PLINPM(NDBASE,NDBASE,1:NDKART,0:NDREVE,0:NDISOS),STAT=IALLOC)
       IF (IALLOC.NE.0) CALL NOALLO('PLINPM','LINAVR')
    ENDIF
!
    IF(ISIMPY.NE.1)THEN
       IF(.NOT.ALLOCATED(ALINPM))THEN
          ALLOCATE(ALINPM(NDBASE,NDBASE,1:NDKART,0:NDREVE,0:NDISOS),STAT=IALLOC)
          IF (IALLOC.NE.0) CALL NOALLO('ALINPM','LINAVR')
       ENDIF
       IF(.NOT.ALLOCATED(PLINPP) .AND. IPAIRI.EQ.1)THEN
          ALLOCATE(PLINPP(NDBASE,NDBASE,1:NDKART,0:NDREVE,0:NDISOS),STAT=IALLOC)
          IF (IALLOC.NE.0) CALL NOALLO('PLINPP','LINAVR')
       ENDIF
    ENDIF
 ENDIF
!=======================================================================
!allocate space for rho and kappa in case they are not allocated yet
!=======================================================================
 IF(.NOT.ALLOCATED(BIG_PP)) THEN
    ALLOCATE(BIG_PP(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('BIG_PP','LINAVR')
 ENDIF
 IF(.NOT.ALLOCATED(DEN_PP)) THEN
    ALLOCATE(DEN_PP(1:NDBASE,1:NDBASE,0:NDREVE,0:NDISOS),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('DEN_PP','SAVDEN')
 ENDIF
!
 IF(.NOT.ALLOCATED(PAI_PM) .AND. IPAIRI.EQ.1)THEN
    ALLOCATE(PAI_PM(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('PAI_PM','LINAVR')
 ENDIF
 IF(.NOT.ALLOCATED(CHI_PM) .AND. IPAIRI.EQ.1)THEN
    ALLOCATE(CHI_PM(1:NDBASE,1:NDBASE,0:NDREVE,0:NDISOS), STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('CHI_PM','LINAVR')
 ENDIF
!
 IF(ISIMPY.NE.1)THEN
    IF(.NOT.ALLOCATED(BIG_PM))THEN
        ALLOCATE (BIG_PM(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
        IF (IALLOC.NE.0) CALL NOALLO('BIG_PM','LINAVR')
    END IF
    IF(.NOT.ALLOCATED(DEN_PM))THEN
        ALLOCATE (DEN_PM(1:NDBASE,1:NDBASE,0:NDREVE,0:NDISOS),STAT=IALLOC)
        IF (IALLOC.NE.0) CALL NOALLO('DEN_PM','LINAVR')
    ENDIF
!
    IF (.NOT.ALLOCATED(PAI_PP) .AND. IPAIRI.EQ.1)THEN
       ALLOCATE (PAI_PP(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
       IF (IALLOC.NE.0) CALL NOALLO('PAI_PP','LINAVR')
    END IF
    IF (.NOT.ALLOCATED(PAI2PP) .AND. IPAIRI.EQ.1)THEN
       ALLOCATE (PAI2PP(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
       IF (IALLOC.NE.0) CALL NOALLO('PAI2PP','LINAVR')
    END IF
    IF (.NOT.ALLOCATED(CHI_PP) .AND. IPAIRI.EQ.1)THEN
       ALLOCATE (CHI_PP(1:NDBASE,1:NDBASE,0:NDREVE,0:NDISOS),STAT=IALLOC)
       IF (IALLOC.NE.0) CALL NOALLO('CHI_PP','LINAVR')
    END IF
 ENDIF
!=======================================================================
!allocate work space
!=======================================================================
 ALLOCATE (HAUXIL(1:NDBASE,1:NDBASE),STAT=IALLOC)
 IF (IALLOC.NE.0) CALL NOALLO('HAUXIL','LINAVR')
 ALLOCATE (HAUXIJ(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
 IF (IALLOC.NE.0) CALL NOALLO('HAUXIJ','LINAVR')
!
 ALLOCATE (RHO_PP(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
 IF (IALLOC.NE.0) CALL NOALLO('RHO_PP','LINAVR')
!
 IF(ISIMPY.NE.1)THEN
    ALLOCATE (RHO_PM(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('RHO_PM','LINAVR')
 ENDIF
!
 IF (IPAIRI.EQ.1) THEN
    ALLOCATE (KAI_PM(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
    IF (IALLOC.NE.0) CALL NOALLO('KAI_PM','LINAVR')
!
    IF(ISIMPY.NE.1)THEN
       ALLOCATE (KAI_PP(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
       IF (IALLOC.NE.0) CALL NOALLO('KAI_PP','LINAVR')
       ALLOCATE (KAI2PP(1:NDBASE,1:NDBASE,0:NDREVE),STAT=IALLOC)
       IF (IALLOC.NE.0) CALL NOALLO('KAI2PP','LINAVR')
    ENDIF
 ENDIF
!=======================================================================
!    ZEROING ARRAYS OF AVERAGE VALUES
!=======================================================================
 DLINSQ(:)=0.0D0
 ELINSQ(:)=0.0D0
 PLINSQ(:)=0.0D0
 TLINSQ(:)=0.0D0
 ALINLI(:)=0.0D0
 PLINLI(:)=0.0D0
 PLINKI(:)=0.0D0
!
 DKINSQ(:)=C_ZERO
 EKINSQ(:)=C_ZERO
 PKINSQ(:)=C_ZERO
 TKINSQ(:)=C_ZERO
 AKINLI(:)=C_ZERO
 PKINLI(:)=C_ZERO
 PKINKI(:)=C_ZERO
!
 IF(JETACM.EQ.2)THEN
                    ALINPP(:,:,:,:,ICHARG)=C_ZERO
    IF(ISIMPY.NE.1) ALINPM(:,:,:,:,ICHARG)=C_ZERO
!
    IF(IPAIRI.EQ.1)THEN
                        PLINPM(:,:,:,:,ICHARG)=C_ZERO
       IF (ISIMPY.NE.1) PLINPP(:,:,:,:,ICHARG)=C_ZERO
    ENDIF
 ENDIF
!
 IF (JETACM.EQ.2) THEN !when calculating matrix element,use mixed rho and kappa
    IF(ITESTA.EQ.NUMITE) THEN!initial den, chi at the first iteration
                        DEN_PP(:,:,:,ICHARG)= BIG_PP(:,:,:)
       IF (ISIMPY.NE.1) DEN_PM(:,:,:,ICHARG)= BIG_PM(:,:,:)
       IF (IPAIRI.EQ.1) THEN
                           CHI_PM(:,:,:,ICHARG)=PAI_PM(:,:,:)
          IF (ISIMPY.NE.1) CHI_PP(:,:,:,ICHARG)=PAI_PP(:,:,:)
       ENDIF
    ENDIF
!
                     RHO_PP(:,:,:)=DEN_PP(:,:,:,ICHARG)
    IF (ISIMPY.NE.1) RHO_PM(:,:,:)=DEN_PM(:,:,:,ICHARG)
    IF (IPAIRI.EQ.1) THEN
       KAI_PM(:,:,:)=CHI_PM(:,:,:,ICHARG)
       IF (ISIMPY.NE.1)THEN
          KAI_PP(:,:,:)=CHI_PP(:,:,:,ICHARG)
          KAI2PP(:,:,:)=CHI_PP(:,:,:,ICHARG)
       ENDIF
    ENDIF
!
 ELSE!when there is no correction or inside RENMAS, use rho and kappa from DENMAC and PAIMAC
                     RHO_PP(:,:,:)=BIG_PP(:,:,:)
    IF (ISIMPY.NE.1) RHO_PM(:,:,:)=BIG_PM(:,:,:)
    IF (IPAIRI.EQ.1) THEN
        KAI_PM(:,:,:)=PAI_PM(:,:,:)
        IF (ISIMPY.NE.1) THEN
           KAI_PP(:,:,:)=PAI_PP(:,:,:)
           KAI2PP(:,:,:)=PAI2PP(:,:,:)
        ENDIF
    ENDIF
 ENDIF
!=======================================================================
!                         X-DIRECTION
!=======================================================================
!
!=======================================================================
![Px]
!=======================================================================
 CALL INT_PX(NXMAXX,NYMAXX,NZMAXX,IPHAPP(0,0,0),HAUXIL,FACTOR)
 HAUXIJ(1:NDBASE,1:NDBASE,0)= UNIT_I*HAUXIL(1:NDBASE,1:NDBASE)
 HAUXIJ(1:NDBASE,1:NDBASE,1)= UNIT_I*HAUXIL(1:NDBASE,1:NDBASE)
!=======================================================================
!Tr[Px.rho]
!=======================================================================
 CALL AVRXSD(HAUXIJ,ALINLI(1),AKINLI(1),RHO_PP)
!
 IF (COR_CM) THEN
    IMSIGN=-1
!=======================================================================
!Tr[(Px.rho)^2] and 2[Px.rho.Px]
!=======================================================================

    IF(ISIMPY.EQ.1 .AND. JETACM.NE.2)then
       CALL EXCXSD(IMSIGN,MREVER,ISIMPY,JETACM,ELINSQ(1),EKINSQ(1),HAUXIJ,RHO_PP)
    ELSEIF(ISIMPY.EQ.0 .AND. JETACM.NE.2)then
       CALL EXCXSD(IMSIGN,MREVER,ISIMPY,JETACM,ELINSQ(1),EKINSQ(1),HAUXIJ,RHO_PP,RPM=RHO_PM)
    ELSEIF(ISIMPY.EQ.1 .AND.  JETACM.EQ.2)then
       CALL EXCXSD(IMSIGN,MREVER,ISIMPY,JETACM,ELINSQ(1),EKINSQ(1),HAUXIJ,RHO_PP,APP=ALINPP(:,:,1,:,ICHARG))
    ELSEIF(ISIMPY.EQ.0 .AND.  JETACM.EQ.2)then
       CALL EXCXSD(IMSIGN,MREVER,ISIMPY,JETACM,ELINSQ(1),EKINSQ(1),HAUXIJ,RHO_PP,RPM=RHO_PM,&
                   APP=ALINPP(:,:,1,:,ICHARG),APM=ALINPM(:,:,1,:,ICHARG))
    ELSE
       stop 'wrong ISIMPY or JETACM in LINAVR'
    ENDIF
!=======================================================================
!Tr[(Px.kappa)(Px.kappa')*] and 2[Px.kappa.Px]
!=======================================================================
    IF (IPAIRI.EQ.1) THEN
       IF(ISIMPY.EQ.1 .AND. JETACM.NE.2)then
          CALL PAIXSD(MREVER,ISIMPY,JETACM,PLINSQ(1),PKINSQ(1),HAUXIJ,KAI_PM)
       ELSEIF(ISIMPY.EQ.0 .AND. JETACM.NE.2)then
          CALL PAIXSD(MREVER,ISIMPY,JETACM,PLINSQ(1),PKINSQ(1),HAUXIJ,KAI_PM,KPP1=KAI_PP,KPP2=KAI2PP)
       ELSEIF(ISIMPY.EQ.1 .AND.  JETACM.EQ.2)then
          CALL PAIXSD(MREVER,ISIMPY,JETACM,PLINSQ(1),PKINSQ(1),HAUXIJ,KAI_PM,PPM=PLINPM(:,:,1,:,ICHARG))
       ELSEIF(ISIMPY.EQ.0 .AND.  JETACM.EQ.2)then
          CALL PAIXSD(MREVER,ISIMPY,JETACM,PLINSQ(1),PKINSQ(1),HAUXIJ,KAI_PM,KPP1=KAI_PP,KPP2=KAI2PP,&
                      PPM=PLINPM(:,:,1,:,ICHARG),PPP=PLINPP(:,:,1,:,ICHARG))
       ELSE
         stop 'wrong ISIMPY or JETACM in LINAVR'
       ENDIF
    END IF

!=======================================================================
![Px^2]
!=======================================================================
    CALL INTPX2(NXMAXX,NYMAXX,NZMAXX,IPHAPP(0,0,0),HAUXIL,FACTOR)
    HAUXIJ(1:NDBASE,1:NDBASE,0)=HAUXIL(1:NDBASE,1:NDBASE)
    HAUXIJ(1:NDBASE,1:NDBASE,1)=HAUXIL(1:NDBASE,1:NDBASE)
!=======================================================================
!Tr[Px^2.rho]
!=======================================================================
    CALL AVRXSD(HAUXIJ,DLINSQ(1),DKINSQ(1),RHO_PP)
 END IF
!=======================================================================
!                         Z-DIRECTION
!=======================================================================
!
!=======================================================================
![Pz]
!=======================================================================
 CALL INT_PZ(NXMAXX,NYMAXX,NZMAXX,IPHAPP(0,0,0),HAUXIL,FACTOR)
 HAUXIJ(1:NDBASE,1:NDBASE,0)=UNIT_I*HAUXIL(1:NDBASE,1:NDBASE)
 HAUXIJ(1:NDBASE,1:NDBASE,1)=UNIT_I*HAUXIL(1:NDBASE,1:NDBASE)
!=======================================================================
!Tr[Pz.rho]
!=======================================================================
 CALL AVRXSD(HAUXIJ,ALINLI(3),AKINLI(3),RHO_PP)
!
 IF (COR_CM) THEN
    IMSIGN=-1
!=======================================================================
!Tr[(Pz.rho)^2] and 2[Pz.rho.Pz]
!=======================================================================

    IF(ISIMPY.EQ.1 .AND. JETACM.NE.2)then
       CALL EXCXSD(IMSIGN,MREVER,ISIMPY,JETACM,ELINSQ(3),EKINSQ(3),HAUXIJ,RHO_PP)
    ELSEIF(ISIMPY.EQ.0 .AND. JETACM.NE.2)then
       CALL EXCXSD(IMSIGN,MREVER,ISIMPY,JETACM,ELINSQ(3),EKINSQ(3),HAUXIJ,RHO_PP,RPM=RHO_PM)
    ELSEIF(ISIMPY.EQ.1 .AND.  JETACM.EQ.2)then
       CALL EXCXSD(IMSIGN,MREVER,ISIMPY,JETACM,ELINSQ(3),EKINSQ(3),HAUXIJ,RHO_PP,APP=ALINPP(:,:,3,:,ICHARG))
    ELSEIF(ISIMPY.EQ.0 .AND.  JETACM.EQ.2)then
       CALL EXCXSD(IMSIGN,MREVER,ISIMPY,JETACM,ELINSQ(3),EKINSQ(3),HAUXIJ,RHO_PP,RPM=RHO_PM,&
                   APP=ALINPP(:,:,3,:,ICHARG), APM=ALINPM(:,:,3,:,ICHARG))
    ELSE
       stop 'wrong ISIMPY or JETACM in LINAVR'
    ENDIF
!=======================================================================
!Tr[(Pz.kappa)(Pz.kappa')*] and 2[Pz.kappa.Pz]
!=======================================================================
    IF (IPAIRI.EQ.1) THEN
       IF(ISIMPY.EQ.1 .AND. JETACM.NE.2)then
          CALL PAIXSD(MREVER,ISIMPY,JETACM,PLINSQ(3),PKINSQ(3),HAUXIJ,KAI_PM)
       ELSEIF(ISIMPY.EQ.0 .AND. JETACM.NE.2)then
          CALL PAIXSD(MREVER,ISIMPY,JETACM,PLINSQ(3),PKINSQ(3),HAUXIJ,KAI_PM,KPP1=KAI_PP,KPP2=KAI2PP)
       ELSEIF(ISIMPY.EQ.1 .AND.  JETACM.EQ.2)then
          CALL PAIXSD(MREVER,ISIMPY,JETACM,PLINSQ(3),PKINSQ(3),HAUXIJ,KAI_PM,PPM=PLINPM(:,:,3,:,ICHARG))
       ELSEIF(ISIMPY.EQ.0 .AND.  JETACM.EQ.2)then
          CALL PAIXSD(MREVER,ISIMPY,JETACM,PLINSQ(3),PKINSQ(3),HAUXIJ,KAI_PM,KPP1=KAI_PP,KPP2=KAI2PP,&
                      PPM=PLINPM(:,:,3,:,ICHARG),PPP=PLINPP(:,:,3,:,ICHARG))
       ELSE
         stop 'wrong ISIMPY or JETACM in LINAVR'
       ENDIF
    END IF
!=======================================================================
![Pz^2]
!=======================================================================
    CALL INTPZ2(NXMAXX,NYMAXX,NZMAXX,IPHAPP(0,0,0),HAUXIL,FACTOR)
    HAUXIJ(1:NDBASE,1:NDBASE,0)=HAUXIL(1:NDBASE,1:NDBASE)
    HAUXIJ(1:NDBASE,1:NDBASE,1)=HAUXIL(1:NDBASE,1:NDBASE)
!=======================================================================
!Tr[(Pz^2.rho]
!=======================================================================
    CALL AVRXSD(HAUXIJ,DLINSQ(3),DKINSQ(3),RHO_PP)
 END IF
!=======================================================================
!                         Y-DIRECTION
!=======================================================================
!
!=======================================================================
![Py]
!=======================================================================
 CALL INT_PY(NXMAXX,NYMAXX,NZMAXX,IPHAPM(0,0,0),HAUXIL,FACTOR)
 HAUXIJ(1:NDBASE,1:NDBASE,0)= UNIT_I*HAUXIL(1:NDBASE,1:NDBASE)
 HAUXIJ(1:NDBASE,1:NDBASE,1)=-UNIT_I*HAUXIL(1:NDBASE,1:NDBASE)
!=======================================================================
!Tr[Py.rho]
!=======================================================================
 IF (ISIMPY.NE.1 ) THEN
 CALL AVRXSO(HAUXIJ,ALINLI(2),AKINLI(2),RHO_PM)
  END IF
!
 IF (COR_CM) THEN
    IMSIGN=1
!=======================================================================
!Tr[(Py.rho)^2] and 2[Py.rho.Py]
!=======================================================================
    IF(ISIMPY.EQ.1 .AND. JETACM.NE.2)then
       CALL EXCXSO(IMSIGN,MREVER,ISIMPY,JETACM,ELINSQ(2),EKINSQ(2),HAUXIJ,RHO_PP)
    ELSEIF(ISIMPY.EQ.0 .AND. JETACM.NE.2)then
       CALL EXCXSO(IMSIGN,MREVER,ISIMPY,JETACM,ELINSQ(2),EKINSQ(2),HAUXIJ,RHO_PP,RPM=RHO_PM)
    ELSEIF(ISIMPY.EQ.1 .AND.  JETACM.EQ.2)then
       CALL EXCXSO(IMSIGN,MREVER,ISIMPY,JETACM,ELINSQ(2),EKINSQ(2),HAUXIJ,RHO_PP,APP=ALINPP(:,:,2,:,ICHARG))
    ELSEIF(ISIMPY.EQ.0 .AND.  JETACM.EQ.2)then
       CALL EXCXSO(IMSIGN,MREVER,ISIMPY,JETACM,ELINSQ(2),EKINSQ(2),HAUXIJ,RHO_PP,RPM=RHO_PM,&
                   APP=ALINPP(:,:,2,:,ICHARG), APM=ALINPM(:,:,2,:,ICHARG))
    ELSE
       stop 'wrong ISIMPY or JETACM in LINAVR'
    ENDIF
!=======================================================================
!Tr[(Py.kappa)(Py.kappa')*] and 2[Py.kappa.Py]
!=======================================================================
    IF (IPAIRI.EQ.1) THEN
       IF(ISIMPY.EQ.1 .AND. JETACM.NE.2)then
          CALL PAIXSO(MREVER,ISIMPY,JETACM,PLINSQ(2),PKINSQ(2),HAUXIJ,KAI_PM)
       ELSEIF(ISIMPY.EQ.0 .AND. JETACM.NE.2)then
          CALL PAIXSO(MREVER,ISIMPY,JETACM,PLINSQ(2),PKINSQ(2),HAUXIJ,KAI_PM,KPP1=KAI_PP,KPP2=KAI2PP)
       ELSEIF(ISIMPY.EQ.1 .AND.  JETACM.EQ.2)then
          CALL PAIXSO(MREVER,ISIMPY,JETACM,PLINSQ(2),PKINSQ(2), HAUXIJ,KAI_PM,PPM=PLINPM(:,:,2,:,ICHARG))
       ELSEIF(ISIMPY.EQ.0 .AND.  JETACM.EQ.2)then
          CALL PAIXSO(MREVER,ISIMPY,JETACM,PLINSQ(2),PKINSQ(2),HAUXIJ,KAI_PM,KPP1=KAI_PP,KPP2=KAI2PP,&
                      PPM=PLINPM(:,:,2,:,ICHARG),PPP=PLINPP(:,:,2,:,ICHARG))
       ELSE
         stop 'wrong ISIMPY or JETACM in LINAVR'
       ENDIF
    END IF
!=======================================================================
![Py^2]
!=======================================================================
    CALL INTPY2(NXMAXX,NYMAXX,NZMAXX,IPHAPP(0,0,0),HAUXIL,FACTOR)
    HAUXIJ(1:NDBASE,1:NDBASE,0)=HAUXIL(1:NDBASE,1:NDBASE)
    HAUXIJ(1:NDBASE,1:NDBASE,1)=HAUXIL(1:NDBASE,1:NDBASE)
!=======================================================================
!Tr[(Py^2.rho]
!=======================================================================
    CALL AVRXSD(HAUXIJ,DLINSQ(2),DKINSQ(2),RHO_PP)
 END IF
!
!=======================================================================
!    HERE ALL TERMS ARE COMBINED INTO THE AVERAGE OF THE SQUARE OF
!    THE TOTAL LINEAR MOMENTUM.
!=======================================================================
 DO KARTEZ=1,NDKART
    TLINSQ(KARTEZ)=DLINSQ(KARTEZ)-ELINSQ(KARTEZ)+PLINSQ(KARTEZ) +ALINLI(KARTEZ)**2
    TKINSQ(KARTEZ)=DKINSQ(KARTEZ)-EKINSQ(KARTEZ)+PKINSQ(KARTEZ) +AKINLI(KARTEZ)**2
 END DO

 DO KARTEZ=1,NDKART
    DLINSQ(0)=DLINSQ(0)+DLINSQ(KARTEZ)
    ELINSQ(0)=ELINSQ(0)+ELINSQ(KARTEZ)
    PLINSQ(0)=PLINSQ(0)+PLINSQ(KARTEZ)
    TLINSQ(0)=TLINSQ(0)+TLINSQ(KARTEZ)
    ALINLI(0)=ALINLI(0)+ALINLI(KARTEZ)**2
    PLINLI(0)=PLINLI(0)+PLINLI(KARTEZ)**2
    PLINKI(0)=PLINKI(0)+PLINKI(KARTEZ)**2
!
    DKINSQ(0)=DKINSQ(0)+DKINSQ(KARTEZ)
    EKINSQ(0)=EKINSQ(0)+EKINSQ(KARTEZ)
    PKINSQ(0)=PKINSQ(0)+PKINSQ(KARTEZ)
    TKINSQ(0)=TKINSQ(0)+TKINSQ(KARTEZ)
    AKINLI(0)=AKINLI(0)+AKINLI(KARTEZ)**2
    PKINLI(0)=PKINLI(0)+PKINLI(KARTEZ)**2
    PKINKI(0)=PKINKI(0)+PKINKI(KARTEZ)**2
 END DO
!=======================================================================
!Deallocate work space
!=======================================================================
 DEALLOCATE (HAUXIL)
 DEALLOCATE (HAUXIJ)
 DEALLOCATE (RHO_PP)
 IF(ISIMPY.NE.1) DEALLOCATE (RHO_PM)
 IF(IPAIRI.EQ.1) DEALLOCATE (KAI_PM)
 IF(IPAIRI.EQ.1 .AND. ISIMPY.NE.1) DEALLOCATE (KAI_PP,KAI2PP)
!=======================================================================
 CALL CPUTIM('LINAVR',0)
 RETURN
END SUBROUTINE LINAVR
!=======================================================================
SUBROUTINE RENMAS(MIN_QP,IRENMA,IDOGOA,&
                  NXHERM,NYHERM,NZHERM,NXMAXX,NYMAXX,NZMAXX,&
                  KSIMTX,JSIMTY,KSIMTZ,KSIGNY,KSIMPY,KSIQTY,&
                  IPAIRI,IPAHFB,IROTAT,ITIREP,KREVER,IPNMIX,&
                  IDEVAR,ITERUN,NMUCOU,ISHIFY,NUMCOU,BOUCOU,&
                  ISYMDE,INIROT,INIINV,INIKAR,ISAWAV,IKERNE,&
                  NOSCIL,NUMETA,NFICOU,FILCOU,ICOULI,ICOULO,&
                                IN_FIX,IZ_FIX,ICANTP,SLOWTP)
!=======================================================================
USE SAVLEF
USE SAVRIG
USE PD_DEN
USE WAVR_L
USE HCOULO
USE OFFDIA, ONLY:LMATCN,LINWLR,VRVMAT,H_KERN
!=======================================================================
USE hfodd_sizes
!=======================================================================
 CHARACTER::FILCOU*68
 CHARACTER::NAMEPN*8,FCOUOF*68
 LOGICAL::PRINIT,COR_CM
 COMPLEX::EKESKY,EKEKIN,EKECOD,EKECOX
!=======================================================================
 COMPLEX::C_ZERO,C_UNIT
 COMPLEX::AVRCMP,BUFCMP
 COMPLEX::EPAIKN,EPAIKP,BAREKK
!=======================================================================
 COMPLEX,ALLOCATABLE::COUHPP(:,:),COUHPM(:,:)
!=======================================================================
 INTEGER,DIMENSION(NDKART)::IREMAS
 COMPLEX,DIMENSION(2)::DETWRK
 COMPLEX,DIMENSION(0:NDSHIF,1:NDKART,0:NDISOS)::AKKERN
 COMPLEX,DIMENSION(0:NDSHIF,0:NDISOS)::OVKERN,EKKERN,EPKERN,DKKERN
 COMPLEX,DIMENSION(0:NDSHIF)::SKKERN,CDKERN,CXKERN,ETKERN
 REAL   ,DIMENSION(0:NDKART)::DLINSQ,ELINSQ,PLINSQ,TLINSQ,ALINLI,PLINLI,PLINKI
 COMPLEX,DIMENSION(0:NDKART)::DKINSQ,EKINSQ,PKINSQ,TKINSQ,AKINLI,PKINLI,PKINKI
 COMPLEX,DIMENSION(0:NDSHIF,1:NDKART)::ETKERS,DKKERS,OVKERS
 REAL   ,DIMENSION(NDKART)::REMASS,GOAMAS
!
 COMPLEX::QMUT_N,QMUT_P,QMUT_T,QMUL_N,QMUL_P,QMUL_T
 REAL   ,DIMENSION(0:NDMULT,0:NDOSCI,0:NDOSCI,1:NDKART)::COM_TMP
 COMPLEX,DIMENSION(0:NDMULT,-NDMULT:NDMULT)::QTN_TMP,QTP_TMP,QTT_TMP,QLN_TMP,QLP_TMP,QLT_TMP
!
 COMMON /QMUTTI/ QMUT_N(0:NDMULT,-NDMULT:NDMULT),QMUT_P(0:NDMULT,-NDMULT:NDMULT),QMUT_T(0:NDMULT,-NDMULT:NDMULT)
 COMMON /QMULTI/ QMUL_N(0:NDMULT,-NDMULT:NDMULT),QMUL_P(0:NDMULT,-NDMULT:NDMULT),QMUL_T(0:NDMULT,-NDMULT:NDMULT)
 COMMON /COEMUL/ COMULT(0:NDMULT,0:NDOSCI,0:NDOSCI,1:NDKART)
!
 COMMON /SVECTO/ BSHIFT(NDKART)
 COMMON /DIMENS/ LDBASE
 COMMON /DIMSTA/ LDTOTS(0:NDISOS),LDSTAT(0:NDISOS),LDUPPE(0:NDISOS),LDTIMU(0:NDISOS)
 COMMON /DISTAN/ DISTAX,DISTAY,DISTAZ
 COMMON /PLANCK/ HBMASS,HBMRPA,HBMINP
 COMMON /RENORM/ HBMREN(NDKART)
 COMMON /RENGOA/ GOAREN(NDKART),GOAREI(NDKART)
 COMMON /CFIPRI/ NFIPRI
 COMMON /ITERTS/ ITESTA,ITESTO,NOITER,NUMITE
!=======================================================================
 CALL CPUTIM('RENMAS',1)
!=======================================================================
 IALLOC=0
 I_SLOW=0
 !IF(MREVER.EQ.0)STOP 'MREVER=0 in RENMAS'
!=======================================================================
 IF (.NOT.ALLOCATED(WARIGH)) THEN
     ALLOCATE (WARIGH(1:NDBASE,1:4*NDSTAT,0:NDSPIN),STAT=IALLOC)
     IF (IALLOC.NE.0) CALL NOALLO('WARIGH','RENMAS')
     WARIGH(:,:,:)=CMPLX(0.0D0,0.0D0)
 END IF
 IF (.NOT.ALLOCATED(WALEFT)) THEN
     ALLOCATE (WALEFT(1:NDBASE,1:4*NDSTAT,0:NDSPIN),STAT=IALLOC)
     IF (IALLOC.NE.0) CALL NOALLO('WALEFT','RENMAS')
     WALEFT(:,:,:)=CMPLX(0.0D0,0.0D0)
 END IF
!=======================================================================
 ALLOCATE (COUHPP(1:NDBASE,1:NDBASE),STAT=IALLOC)
 IF (IALLOC.NE.0) CALL NOALLO('COUHPP','RENMAS')
 ALLOCATE (COUHPM(1:NDBASE,1:NDBASE),STAT=IALLOC)
 IF (IALLOC.NE.0) CALL NOALLO('COUHPM','RENMAS')
!=======================================================================
!         THIS SUBROUTINE CALCULATES THE MASS RENORMALIZATION FACTORS
!=======================================================================
 C_ZERO=CMPLX(0.0D0,0.0D0)
 C_UNIT=CMPLX(1.0D0,0.0D0)
 SIF=1.0D-8
!
 NAMEPN='DUMMY   '
 PRINIT=.FALSE.
 JMULMO=NMUCOU
 COR_CM=.TRUE.
 JETACM=0
!=======================================================================
!         STORING THE COULOMB MATRIX ELEMENTS. THESE WILL BE RETURNED
!         AFTER THE MASS RENORMALIZATION.
!=======================================================================
                  COUHPP(:,:)=HPPCOU(:,:)
 IF (KSIMPY.NE.1) COUHPM(:,:)=HPMCOU(:,:)
 COM_TMP=COMULT
 QTN_TMP=QMUT_N
 QTP_TMP=QMUT_P
 QTT_TMP=QMUT_T
 QLN_TMP=QMUL_N
 QLP_TMP=QMUL_P
 QLT_TMP=QMUL_T
!
 ICANON=ABS(ICANTP)
 IEXTEN=0
!=======================================================================
!STORE SYMMETRY PARAMETERS,
!THESE VALUES MAY BE CHANGED DEPENDING ON SHIFTING DIRECTION
!=======================================================================
 ISIMPY=KSIMPY
 ISIGNY=KSIGNY
 ISIQTY=KSIQTY
 ISIMTX=KSIMTX
 ISIMTZ=KSIMTZ
 MREVER=KREVER
!
 ISIQTY = 0
 IF(ABS(DISTAX).GT. SIF) THEN
    ISIGNY = 0
    ISIMTZ = 0
 ENDIF
!
 IF(ABS(DISTAY).GT. SIF) THEN
    ISIMPY = 0
    ISIMTX = 0
    ISIMTZ = 0
    IF(KREVER.EQ.0) IEXTEN=1!when breaking time reversal, the stored wave function should be modified to include all simplex blocks
 ENDIF
!
 IF(ABS(DISTAZ).GT. SIF) THEN
    ISIGNY = 0
    ISIMTX = 0
 ENDIF
!
 IF(IEXTEN.EQ.1) MREVER=1
 IF(IPAHFB.NE.1) ICANON=0!
!
!=======================================================================
!CHECK INPUT PARAMTER
!=======================================================================
 IF (KREVER.EQ.0 .AND. ICANON.EQ.0 .AND. ABS(DISTAY).GT. 1.0d-8)&
    STOP 'SHIFTING ALONG Y WITH MREVER=0 IN RENMAS'
! IF (ICANON.GT.0 .AND. KSIMPY.NE.1) STOP 'ICANON>0 WITH ISIMPY =0 NOT IMPLINENTED YET in RENMAS'
 IF (IRENMA.GT.NDSHIF) STOP ' IRENMA.GT.NDSHIF IN RENMAS'
! IF (IEXTEN.EQ.1 .AND. (KREVER.EQ.1 .OR. ICANON.EQ.0)) STOP 'IEXTEN=1 needs KREVER=0 and ICANON>0 in RENMAS'
! IF(ICANON.GT.1)STOP 'ICANON>1 in RENMAS'
!=======================================================================
!    TREAT COULOMB PART WHEN SYMMETRY IS BROKEN INSIDE THIS ROUTINE
!=======================================================================
 IF(ISIMPY.NE.KSIMPY .OR. ISIGNY.NE.KSIGNY .OR. ISIQTY.NE.KSIQTY )THEN
    CALL MULCOU(ISIMPY,ISIGNY,ISIQTY,NMUCOU)
    FCOUOF=TRIM(FILCOU)//'OFD'
    IF(ITESTA.EQ.NUMITE)THEN!DO CALCULATION ONLY AT THE FIRST ITERATION
      ICOULI=0
      ICOULO=1
      CALL PRECOU(NOSCIL,NUMCOU,NUMETA,BOUCOU,NFICOU,FCOUOF,ICOULI,ICOULO)
    ELSE!READ FROM FILE
      INPCOU=1
      I_SUCC=1
      OPEN(NFICOU,FILE=FCOUOF,STATUS='OLD',FORM='UNFORMATTED')
      CALL RECOUL(NFICOU,INPCOU,I_SUCC,NUMCOU,NUMETA,BOUCOU)
      CLOSE(NFICOU)
    ENDIF
 ENDIF
!
 ETKERS(:,:)=C_ZERO
 DKKERS(:,:)=C_ZERO
 OVKERS(:,:)=C_ZERO
!=======================================================================
!         HERE BEGIN THE LOOP THE SHIFT VECTORS.
!=======================================================================!
 IF (IRENMA.GT.NDSHIF) STOP ' IRENMA.GT.NDSHIF IN RENMAS'
 MDKART=NDKART!UPPER LIMIT OF DIRECTION LOOP
 DO INSHIF=IRENMA,0,-1!SHIFT DISTANCE
    IF (INSHIF.EQ.0) MDKART=1!WHEN THERE IS NO SHIFT,ONLY DO 1 TIME
    DO KARTEZ=1,MDKART!X,Y,Z DIRECTION
     BSHIFT(:)=0.0d0
     IF (KARTEZ.EQ.1) BSHIFT(1)=INSHIF*DISTAX
     IF (KARTEZ.EQ.2) BSHIFT(2)=INSHIF*DISTAY
     IF (KARTEZ.EQ.3) BSHIFT(3)=INSHIF*DISTAZ
     IF ((INSHIF.EQ.0).OR.(ABS(BSHIFT(KARTEZ)).GT.SIF)) THEN
       DO ICHARG=0,NDISOS
          IF(IEXTEN.EQ.1)THEN
             LDUPPE(ICHARG)=LDTIMU(ICHARG)*2
             LDSTAT(ICHARG)=LDTIMU(ICHARG)*3
             LDTOTS(ICHARG)=LDTIMU(ICHARG)*4
          ENDIF
          LTIMUP=LDTIMU(ICHARG)
          LUPPER=LDUPPE(ICHARG)
          LTOTAL=LDTOTS(ICHARG)
!=======================================================================
!FETCHING THE LEFT WAVE FUNCTIONS FOR THE GIVEN CHARGE
!AND CALCULATE ROTATED RIGHT WAVE FUNCTIONS MULTIPLYING THE INVERTED OVERLAP MATRIX
!USE 'KSIMPY' HERE,SINCE CHANGE IN WAVEFUNCTION HAPPENS AFTER THIS ROUTINE IS CALLED
!USE 'mrever' since both blocks are already stored, by change lduppe and so on, the extension is done
!=======================================================================
          CALL LINWLR(MREVER,KSIMPY,ICANON,IDOGOA,IPAIRI,INSHIF,IEXTEN,&
                      LDBASE,LTIMUP,LUPPER,LTOTAL,ICHARG,OVKERN(INSHIF,ICHARG))
!=======================================================================
!        RHO MATRIX
!=======================================================================
          CALL DENMAC(MREVER,ICHARG,ISIMPY,IPAHFB,WALEFT,WARIGH)
!=======================================================================
!        KAPPA MATRIX
!=======================================================================
          IF (IDOGOA.EQ.0) THEN
             IF (IPAIRI.EQ.1 .AND. COR_CM) CALL PAIMAC(MREVER,ICHARG,ISIMPY,WALEFT,WARIGH,IKERNE)
!=======================================================================
!        TERMS FOR <|J^2|Ri>
!=======================================================================
             CALL LINAVR(NXMAXX,NYMAXX,NZMAXX,ISIMPY,ISIGNY,ISIQTY,&
                         IPAIRI,ICHARG,JETACM,IROTAT,MREVER,COR_CM,&
                         DLINSQ,ELINSQ,PLINSQ,TLINSQ,ALINLI,PLINLI,PLINKI,&
                         DKINSQ,EKINSQ,PKINSQ,TKINSQ,AKINLI,PKINLI,PKINKI)
!
             DKKERN(INSHIF,ICHARG)=TKINSQ(0)
             AKKERN(INSHIF,1,ICHARG)=AKINLI(1)
             AKKERN(INSHIF,2,ICHARG)=AKINLI(2)
             AKKERN(INSHIF,3,ICHARG)=AKINLI(3)
          ELSE
             DKKERN(INSHIF,ICHARG)=C_ZERO
             AKKERN(INSHIF,1,ICHARG)=C_ZERO
             AKKERN(INSHIF,2,ICHARG)=C_ZERO
             AKKERN(INSHIF,3,ICHARG)=C_ZERO
          END IF

!=======================================================================
!        <|H|Ri>
!=======================================================================
          ISKYRM=ICHARG
          CALL H_KERN(NXHERM,NYHERM,NZHERM,NXMAXX,NYMAXX,NZMAXX,&
                      ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY,&
                      IPAHFB,MREVER,ICHARG,MIN_QP,IPNMIX,ITIREP,NAMEPN,&
                      PRINIT,IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,&
                      INIKAR,ISAWAV,IKERNE,ISHIFY,JMULMO,NUMCOU,&
                      BOUCOU,IN_FIX,IZ_FIX,ISKYRM,&
                      EKKERN(INSHIF,ICHARG),EPKERN(INSHIF,ICHARG),&
                      BAREKE,BAREKK,CDKERN(INSHIF),CXKERN(INSHIF),SKKERN(INSHIF))
       ENDDO!END ICHARG
!=======================================================================
!        SKYREM PART OF <|H|Ri>
!=======================================================================
!       CALL TRUTOD(NXHERM,NYHERM,NZHERM)
!       CALL ESKYRM(NXHERM,NYHERM,NZHERM,ENESKY,ENEVEN,ENEODD,&
!                   ENREAR,ENE_W0,EEVEW0,EODDW0,EKESKY,LDPNMX)
!       SKKERN(INSHIF)=EKESKY
!=======================================================================
!        <|H|Ri> and <|J^2|Ri>
!=======================================================================
       ETKERS(INSHIF,KARTEZ)=EKKERN(INSHIF,0)+EKKERN(INSHIF,1)&
                            +EPKERN(INSHIF,0)+EPKERN(INSHIF,1)&
                            +CDKERN(INSHIF)+CXKERN(INSHIF)+SKKERN(INSHIF)
       DKKERS(INSHIF,KARTEZ)= DKKERN(INSHIF,0)+DKKERN(INSHIF,1)&
               +2*AKKERN(INSHIF,1,0)*AKKERN(INSHIF,1,1)&
               +2*AKKERN(INSHIF,2,0)*AKKERN(INSHIF,2,1)&
               +2*AKKERN(INSHIF,3,0)*AKKERN(INSHIF,3,1)
       OVKERS(INSHIF,KARTEZ)= OVKERN(INSHIF,0)*OVKERN(INSHIF,1)
!
!     write(*,*)'e',etkers(inshif,kartez),inshif,kartez
!     write(*,*)'p',dkkers(inshif,kartez),inshif,kartez
     ENDIF! IF ((BSHIFT(KARTEZ).NE.0.0D0).OR.(INSHIF.EQ.0))
    ENDDO! End KARTEZ
 END DO ! End INSHIF
!=======================================================================
!        BELOW WE DEFINE THE LIPKIN PARAMETERS, WHICH ARE POSITIVE
!=======================================================================
 REMASS(:)=0.0d0
 GOAMAS(:)=0.0d0
 ABX=ABS(DISTAX)
 ABY=ABS(DISTAY)
 ABZ=ABS(DISTAZ)

 IF(ABX.GT.SIF)THEN!SHIFT ALONG X AXIS EXISTS
!
    IF (IDOGOA.EQ.0) THEN
       REMASS(1)=REAL((ETKERS(IRENMA,1)-ETKERS(0,1))/(DKKERS(IRENMA,1)-DKKERS(0,1)))
    ELSE
       TSHIFT=REAL(IRENMA**2,KIND=8)*DISTAX**2
       H2PARA=-2*REAL(ETKERS(IRENMA,1)-ETKERS(0,1)) /TSHIFT
       A_PARA=-2*LOG(REAL(OVKERS(IRENMA,1)))/TSHIFT
       GOAMAS(1)=H2PARA/(2*A_PARA**2)
!       REMASS(1)=GOAMAS(1)
    END IF
 ENDIF
!
 IF(ABY.GT.SIF)THEN!SHIFT ALONG Y AXIS EXISTS
!
    IF (IDOGOA.EQ.0) THEN
       REMASS(2)=REAL((ETKERS(IRENMA,2)-ETKERS(0,1))/(DKKERS(IRENMA,2)-DKKERS(0,1)))
    ELSE
       TSHIFT=REAL(IRENMA**2,KIND=8)*DISTAY**2
       H2PARA=-2*REAL(ETKERS(IRENMA,2)-ETKERS(0,1)) /TSHIFT
       A_PARA=-2*LOG(REAL(OVKERS(IRENMA,2)))/TSHIFT
       GOAMAS(2)=H2PARA/(2*A_PARA**2)
!       REMASS(2)=GOAMAS(2)
    END IF
 ENDIF
!
 IF(ABZ.GT.SIF)THEN!SHIFT ALONG Y AXIS EXISTS

!
    IF (IDOGOA.EQ.0) THEN
       REMASS(3)=REAL((ETKERS(IRENMA,3)-ETKERS(0,1))/(DKKERS(IRENMA,3)-DKKERS(0,1)))
    ELSE
       TSHIFT=REAL(IRENMA**2,KIND=8)*DISTAZ**2
       H2PARA=-2*REAL(ETKERS(IRENMA,3)-ETKERS(0,1)) /TSHIFT
       A_PARA=-2*LOG(REAL(OVKERS(IRENMA,3)))/TSHIFT
       GOAMAS(3)=H2PARA/(2*A_PARA**2)
!       REMASS(3)=GOAMAS(3)
    END IF
 ENDIF

 IF (MREVER.EQ.0) GOAMAS(:)=GOAMAS(:)/4
! MREVER.EQ.0 with shifting or rotating is not allowed yet
! But this option is still kept inside RENMAS and RENINE.
 IF (IDOGOA.EQ.1) REMASS(:)=GOAMAS(:)
!
 IF(ABX.LT.SIF) THEN!NO SHIFT ALONG X AXIS
    IF(DISTAY.GT.SIF) THEN
       REMASS(1)=REMASS(2)
       GOAMAS(1)=GOAMAS(2)
    ELSE
       REMASS(1)=REMASS(3)
       GOAMAS(1)=GOAMAS(3)
    ENDIF
 ENDIF
!
 IF(ABY.LT.SIF) THEN!NO SHIFT ALONG Y AXIS
    IF(DISTAX.GT.SIF) THEN
       REMASS(2)=REMASS(1)
       GOAMAS(2)=GOAMAS(1)
     ELSE
       REMASS(2)=REMASS(3)
       GOAMAS(2)=GOAMAS(3)
     ENDIF
 ENDIF
!
 IF(ABZ.LT.SIF) THEN!NO SHIFT ALONG Z AXIS
    IF(DISTAX.GT.SIF) THEN
       REMASS(3)=REMASS(1)
       GOAMAS(3)=GOAMAS(1)
     ELSE
       REMASS(3)=REMASS(2)
       GOAMAS(3)=GOAMAS(2)
     ENDIF
 ENDIF
!=======================================================================
!   Only keep a given digital of lipkin parameter to accelate converge
!=======================================================================
 ITRUNC=1
 IDIGTA=8
 IF(ITRUNC.EQ.1)THEN
    IMASK=10**IDIGTA
    IREMAS(:)=NINT(REMASS(:)*IMASK)
    REMASS(:)=IREMAS(:)/(IMASK*1.0d0)
 ENDIF
!
!=======================================================================
!       BELOW WE DEFINE THE COEFFICIENTS OF <P^2>, WHICH ARE NEGATIVE
!=======================================================================
 IF (I_SLOW.NE.1) THEN
     HBMREN(1)=HBMREN(1)*SLOWTP-REMASS(1)*(1-SLOWTP)
     HBMREN(2)=HBMREN(2)*SLOWTP-REMASS(2)*(1-SLOWTP)
     HBMREN(3)=HBMREN(3)*SLOWTP-REMASS(3)*(1-SLOWTP)
 ELSE
     HBMREN(1)=-REMASS(1)
     HBMREN(2)=-REMASS(2)
     HBMREN(3)=-REMASS(3)
 ENDIF
 !
 GOAREN(1)=-GOAMAS(1)
 GOAREN(2)=-GOAMAS(2)
 GOAREN(3)=-GOAMAS(3)
 !write(*,*)'renmas',remass(1),goamas(1)
!=======================================================================
!    RESTORE COULOMB PART WHEN SYMMETRY IS BROKEN INSIDE THIS ROUTINE
!=======================================================================
 IF(ISIMPY.NE.KSIMPY .OR. ISIGNY.NE.KSIGNY .OR. ISIQTY.NE.KSIQTY )THEN
    CALL MULCOU(KSIMPY,KSIGNY,KSIQTY,NMUCOU)
    OPEN(NFICOU,FILE=FILCOU,STATUS='OLD',FORM='UNFORMATTED')
    INPCOU=1
    I_SUCC=1
    CALL RECOUL(NFICOU,INPCOU,I_SUCC,NUMCOU,NUMETA,BOUCOU)
    CLOSE(NFICOU)
 ENDIF
                  HPPCOU(:,:)=COUHPP(:,:)
 IF (ISIMPY.NE.1) HPMCOU(:,:)=COUHPM(:,:)
 COMULT=COM_TMP
 QMUT_N=QTN_TMP
 QMUT_P=QTP_TMP
 QMUT_T=QTT_TMP
 QMUL_N=QLN_TMP
 QMUL_P=QLP_TMP
 QMUL_T=QLT_TMP
!=======================================================================
 CALL CPUTIM('RENMAS',0)
 RETURN
END SUBROUTINE RENMAS
!=======================================================================
