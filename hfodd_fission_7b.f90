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
 !         with the U.S. Department of Energy (DoE). This work was
 !         produced at the Lawrence Livermore National Laboratory under
 !         Contract No. DE-AC52-07NA27344 with the DoE.
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

 !----------------------------------------------------------------------!
 !  This module sets the precision                                      !
 !----------------------------------------------------------------------!
 Module hfodd_fission_precision
    Implicit None
    Integer, Parameter :: ipr=Kind(1), pr=Kind(1.0)
    Integer(ipr) :: debug = 0
 End Module hfodd_fission_precision

 ! ==================================================================== !
 !                                                                      !
 !                  GLOBAL FRAGMENT PROPERTIES PACKAGE                  !
 !                                                                      !
 ! ==================================================================== !

 !----------------------------------------------------------------------!
 !                                                                      !
 !  This module contains a number of routines and functions dealing     !
 !  with the properties of the fission fragments. It contains all       !
 !  routines related to the quantum mechanical operators used to study  !
 !  fragment poperties. In particular, it contains the definition of    !
 !  the Gaussian neck, fragment distance and the mass asymmetry opera-  !
 !  tors, their expectation value and matrix elements in the simplex HO !
 !  basis. It also contains the routines to compute global fragment     !
 !  properties such as charge, mass and multipole moments.              !
 !                                                                      !
 !  Inputs:                                                             !
 !    - DENSIT ...: total isoscalar density for the compound nucleus    !
 !    - DENCHA ...: charge density for the compound nucleus             !
 !                                                                      !
 !  Outputs:                                                            !
 !    - Z_NECK .........: position of the neck (in QNFIND)              !
 !    - CENLEF,CENRIG ..: position of the center of mass for the left   !
 !                        and right fragments (in center_of_mass)       !
 !    - RIGMAS,RIGCHA ..: mass and charge of the right fragment (in     !
 !                        RFRAGM)                                       !
 !    - QLMLEF,QLMRIG ..: multipole moments in the left and right       !
 !                        fragments                                     !
 !                                                                      !
 !----------------------------------------------------------------------!

 Module hfodd_fission_fragments

    Use hfodd_sizes
    Use hfodd_fission_precision

    Implicit None

    Integer(ipr), PRIVATE, SAVE :: NPOINT=1000
    Integer(ipr), PRIVATE, SAVE :: MXHERM=1,MYHERM=1,MZHERM=1
    Real(pr), PRIVATE, SAVE :: prec_heaviside = 0.001_pr

    ! New constraints
    Real(pr), PUBLIC, SAVE :: distance_prefragments=0.0_pr, asymmetry_prefragments=0.0_pr
    Real(pr), PUBLIC, SAVE :: mass_left=0.0_pr, mass_right=0.0_pr
    ! Flags to switch on/off calculations
    Logical, PUBLIC, SAVE :: distance_is_calculated = .False.
    Logical, PUBLIC, SAVE :: asymmetry_is_calculated = .False.
    ! Constraints on the neck
    Real(pr), PUBLIC, SAVE :: AN_VAL=1.0_pr, Q_NECK=0.0_pr
    Real(pr), PUBLIC, SAVE :: X_NECK=0.0_pr,Y_NECK=0.0_pr,Z_NECK=0.0_pr
    Real(pr), PUBLIC :: ZNMINI=0.0_pr, ZNMAXI=0.0_pr

 Contains

    !---------------------------------------------------------------------!
    ! Function computes the expectation value of the multipole moment     !
    ! operators (l,m) = (LAMACT, MIUACT) in the left and right fragments, !
    ! assuming the position of the neck is Z_POSI                         !
    !---------------------------------------------------------------------!
    Subroutine QLMFRA(NXHERM,NYHERM,NZHERM,Z_POSI,LAMACT,MIUACT,QLMLEF,QLMRIG,CENLEF,CENRIG,I_TYPE)

      Integer(ipr), INTENT(IN) :: NXHERM,NYHERM,NZHERM,LAMACT,MIUACT,I_TYPE
      Real(pr), INTENT(IN) :: CENLEF,CENRIG
      Real(pr), INTENT(INOUT) :: Z_POSI,QLMLEF,QLMRIG

      Real(pr) :: DENSIT,DENCHA,EXPAUX,FOURWG,FOURPT,HOMSCA

      Integer(ipr) :: IX,IY,IZ,IERROR,i
      Real(pr) :: DERIV1,DERIVN,XBEGIN,XFINIS,XARGUM,X_STEP,pi,dh
      Real(pr) :: XPOINT,YPOINT,ZPOINT,DENLOC,W_HERM,QLMSUM,QLMVAL,RESULT
      Real(pr) :: RADIUS,COSTHE,ANGPHI,epsilon,one
      Real(pr), Dimension(1:2) :: SPHHAR
      Real(pr), Allocatable :: FUNCTI(:),XVALUE(:),ZINTER(:),AUXSTO(:)
      Real(pr), Allocatable :: DINTER(:,:,:)

      COMMON                                        &
             /DENTOT/ DENSIT(NDXHRM,NDYHRM,NDZHRM), &
                      DENCHA(NDXHRM,NDYHRM,NDZHRM)
      COMMON                                        &
             /DENEXP/ EXPAUX(NDXHRM,NDYHRM,NDZHRM)
      COMMON                                        &
             /INTSTO/ FOURWG(1:NDGAUS,1:NDKART),    &
                      FOURPT(1:NDGAUS,1:NDKART)
      COMMON                                        &
             /SCALNG/ HOMSCA(NDKART)

      pi = 4.0_pr*Atan(1.0_pr)

      Allocate(FUNCTI(1:NZHERM),XVALUE(1:NZHERM))
      Allocate(AUXSTO(1:NZHERM),ZINTER(1:NPOINT))
      Allocate(DINTER(1:NZHERM,1:NZHERM,1:NPOINT))

      ! Define boundaries for the integration over z, see comments in subroutine center_of_mass()
      XBEGIN=FOURPT(NZHERM,3)/HOMSCA(3) + 1.D-14
      XFINIS=Z_POSI/HOMSCA(3)
      ! Boundary test: make sure there is a valid integration interval
      If(Z_POSI.LT.XBEGIN) Then
         Write(6,'("ATTENTION! In QLMFRA(): changing neck position from ",f20.14," to ",f20.14)') &
                    Z_POSI,XBEGIN + 1.D-12
         Z_POSI = XBEGIN + 1.D-12
      End If

      Do IY=1,NYHERM
         Do IX=1,NXHERM

            ! See comments in subroutine center_of_mass()
            ! Total density
            If(I_TYPE.Eq.1) Then
               Do IZ=1,NZHERM
                  XVALUE(IZ)=FOURPT(NZHERM-IZ+1,3)/HOMSCA(3)
                  FUNCTI(IZ)=DENSIT(IX,IY,NZHERM-IZ+1) &
                            /Exp(-2.0_pr*FOURPT(IX,1)**2)/Exp(-2.0_pr*FOURPT(IY,2)**2)
               End Do
            End If
            ! Proton (=charge in HFODD) density
            If(I_TYPE.Eq.2) Then
               Do IZ=1,NZHERM
                  XVALUE(IZ)=FOURPT(NZHERM-IZ+1,3)/HOMSCA(3)
                  FUNCTI(IZ)=DENCHA(IX,IY,NZHERM-IZ+1) &
                            /Exp(-2.0_pr*FOURPT(IX,1)**2)/Exp(-2.0_pr*FOURPT(IY,2)**2)
               End Do
            End If

            ! Approximation of derivatives at the boundaries
            dh=XVALUE(2)-XVALUE(1) ! >0
            DERIV1=+0.5_pr*(-3.0_pr*FUNCTI(1)     +4.0_pr*FUNCTI(2)       -FUNCTI(3))       /dh
            DERIVN=-0.5_pr*(-3.0_pr*FUNCTI(NZHERM)+4.0_pr*FUNCTI(NZHERM-1)-FUNCTI(NZHERM-2))/dh

            ! Define spline coefficients
            IERROR=0; Call SPLINE(XVALUE,FUNCTI,NZHERM,DERIV1,DERIVN,AUXSTO,IERROR)

            ! Define interpolated mesh and value of the function on that mesh
            If(IERROR.Eq.0) Then
               Do i=1,NPOINT
                  XARGUM=XBEGIN+(XFINIS-XBEGIN)*Real(i-1,Kind=pr)/Real(NPOINT-1,Kind=pr)
                  ZINTER(i)=XARGUM
                  DINTER(IX,IY,i)=SPLINT(XVALUE,FUNCTI,AUXSTO,NZHERM,XARGUM)
               End Do
            Else
               Write(6,'("In QLMFRA() - Error in interpolating the function (right fragment)!")')
            End If

         End Do
      End Do

      X_STEP=(ZINTER(NPOINT)-ZINTER(1))/Real(NPOINT-1,Kind=pr)

      Deallocate(FUNCTI); Allocate(FUNCTI(1:NPOINT))

      ! Computing the expectation values of multipole moment for the left fragment
      epsilon=1.D-14; one=1.0_pr
      Do i=1,NPOINT
         FUNCTI(i)=0.0_pr; QLMSUM=0.0_pr
         Do IY=1,NYHERM
            Do IX=1,NXHERM

               DENLOC=DINTER(IX,IY,i)
               W_HERM=FOURWG(IX,1)*FOURWG(IY,2)

               ! Coordinates (x,y,z) for the multipole moments must be in fermis
               XPOINT=FOURPT(IX,1)/HOMSCA(1)
               YPOINT=FOURPT(IY,2)/HOMSCA(2)
               ZPOINT=ZINTER(i)-CENLEF ! shift with respect to the c.o.m. of left fragment

               RADIUS=Sqrt(XPOINT**2+YPOINT**2+ZPOINT**2)

               ! Angle theta
               If(RADIUS.Le.epsilon) Then
                  COSTHE=0.0_pr
               Else
                  COSTHE=ZPOINT/RADIUS
               End If
               If(Abs(COSTHE-one).Le.epsilon) Then
                  COSTHE=one-2.0_pr*epsilon
               End If
               If(Abs(COSTHE+one).Le.epsilon) Then
                  COSTHE=2.0_pr*epsilon-one
               End If

               ! Angle phi
               If(Abs(XPOINT).Le.epsilon) Then
                  If(YPOINT*XPOINT.Lt.0.0_pr) ANGPHI=-0.5_pr*pi
                  If(YPOINT*XPOINT.Gt.0.0_pr) ANGPHI=+0.5_pr*pi
               Else
                  ANGPHI=Atan(YPOINT/XPOINT)
               End If

               Call DEFSPH(LAMACT,MIUACT,COSTHE,ANGPHI,SPHHAR)

               If(LAMACT.Eq.0) Then
                  QLMVAL=1.0_pr
               Else
                  If(RADIUS.Ge.epsilon) Then
                     QLMVAL=RADIUS**(LAMACT)*SPHHAR(1)
                  Else
                     QLMVAL=0.0_pr
                  End If
               End If
               QLMSUM=QLMSUM+W_HERM*DENLOC*QLMVAL

            End Do
         End Do
         FUNCTI(i)=QLMSUM
      End Do

      ! Integrating over z, from -infty to zN
      Call SIMP38(FUNCTI,NPOINT,X_STEP,RESULT)

      QLMLEF=RESULT

      Deallocate(FUNCTI,XVALUE,DINTER,AUXSTO,ZINTER)

      ! Repeating the procedure for the right fragment

      Allocate(FUNCTI(1:NZHERM),XVALUE(1:NZHERM))
      Allocate(DINTER(1:NZHERM,1:NZHERM,1:NPOINT))
      Allocate(AUXSTO(1:NZHERM),ZINTER(1:NPOINT))

      XBEGIN=Z_POSI/HOMSCA(3)
      XFINIS=FOURPT(1,3)/HOMSCA(3) - 1.D-14
      If(Z_POSI.Gt.XFINIS) Then
         Write(6,'("ATTENTION! In QLMFRA(): changing neck position from ",f20.14," to ",f20.14)') &
                    Z_POSI,XFINIS - 1.D-12
         Z_POSI = XFINIS - 1.D-12
      End If

      Do IY=1,NYHERM
         Do IX=1,NXHERM

            If(I_TYPE.Eq.1) Then
               Do IZ=1,NZHERM
                  XVALUE(IZ)=FOURPT(NZHERM-IZ+1,3)/HOMSCA(3)
                  FUNCTI(IZ)=DENSIT(IX,IY,NZHERM-IZ+1) &
                            /Exp(-2.0_pr*FOURPT(IX,1)**2)/Exp(-2.0_pr*FOURPT(IY,2)**2)
               End Do
            End If
            If(I_TYPE.Eq.2) Then
               Do IZ=1,NZHERM
                  XVALUE(IZ)=FOURPT(NZHERM-IZ+1,3)/HOMSCA(3)
                  FUNCTI(IZ)=DENCHA(IX,IY,NZHERM-IZ+1) &
                            /Exp(-2.0_pr*FOURPT(IX,1)**2)/Exp(-2.0_pr*FOURPT(IY,2)**2)
               End Do
            End If

            dh=XVALUE(2)-XVALUE(1) ! >0
            DERIV1=+0.5_pr*(-3.0_pr*FUNCTI(1)     +4.0_pr*FUNCTI(2)       -FUNCTI(3))       /dh
            DERIVN=-0.5_pr*(-3.0_pr*FUNCTI(NZHERM)+4.0_pr*FUNCTI(NZHERM-1)-FUNCTI(NZHERM-2))/dh

            IERROR=0; Call SPLINE(XVALUE,FUNCTI,NZHERM,DERIV1,DERIVN,AUXSTO,IERROR)

            If(IERROR.Eq.0) Then
               Do i=1,NPOINT
                  XARGUM=XBEGIN+(XFINIS-XBEGIN)*Real(i-1,Kind=pr)/Real(NPOINT-1,Kind=pr)
                  ZINTER(i)=XARGUM
                  DINTER(IX,IY,i)=SPLINT(XVALUE,FUNCTI,AUXSTO,NZHERM,XARGUM)
               End Do
            Else
               Write(6,'("In QLMFRA() - Error in interpolating the function (left fragment)!")')
            End If

         End Do
      End Do

      X_STEP=(ZINTER(NPOINT)-ZINTER(1))/Real(NPOINT-1,Kind=pr)

      Deallocate(FUNCTI); Allocate(FUNCTI(1:NPOINT))

      epsilon=1.D-14; one=1.0_pr
      Do i=1,NPOINT
         FUNCTI(i)=0.0_pr; QLMSUM=0.0_pr
         Do IY=1,NYHERM
            Do IX=1,NXHERM

               DENLOC=DINTER(IX,IY,i)
               W_HERM=FOURWG(IX,1)*FOURWG(IY,2)

               ! Coordinates (x,y,z) in fermis
               XPOINT=FOURPT(IX,1)/HOMSCA(1)
               YPOINT=FOURPT(IY,2)/HOMSCA(2)
               ZPOINT=ZINTER(i)-CENRIG ! shift with respect to the c.o.m. of right fragment

               RADIUS=Sqrt(XPOINT**2+YPOINT**2+ZPOINT**2)

               ! Angle theta
               If(RADIUS.Le.epsilon) Then
                  COSTHE=0.0_pr
               Else
                  COSTHE=ZPOINT/RADIUS
               End If
               If(Abs(COSTHE-one).Le.epsilon) Then
                  COSTHE=one-2.0_pr*epsilon
               End If
               If(Abs(COSTHE+one).Le.epsilon) Then
                  COSTHE=2.0_pr*epsilon-one
               End If

               ! Angle phi
               If(Abs(XPOINT).Le.epsilon) Then
                  If(YPOINT*XPOINT.Lt.0.0_pr) ANGPHI=-0.5_pr*pi
                  If(YPOINT*XPOINT.Gt.0.0_pr) ANGPHI=+0.5_pr*pi
               Else
                  ANGPHI=Atan(YPOINT/XPOINT)
               End If

               Call DEFSPH(LAMACT,MIUACT,COSTHE,ANGPHI,SPHHAR)

               If(LAMACT.Eq.0) Then
                  QLMVAL=1.0_pr
               Else
                  If(RADIUS.Ge.epsilon) Then
                     QLMVAL=RADIUS**(LAMACT)*SPHHAR(1)
                  Else
                     QLMVAL=0.0_pr
                  End If
               End If
               QLMSUM=QLMSUM+W_HERM*DENLOC*QLMVAL

            End Do
         End Do
         FUNCTI(i)=QLMSUM
      End Do

      Call SIMP38(FUNCTI,NPOINT,X_STEP,RESULT)

      QLMRIG=RESULT

      Deallocate(FUNCTI,XVALUE,DINTER,AUXSTO,ZINTER)

    End Subroutine QLMFRA

    !---------------------------------------------------------------------!
    ! The routine computes the value of the multipole moment operators    !
    ! (l,m) = (LAMACT, MIUACT) at angles theta and phi such that          !
    !        cos(theta) = COSTHE                                          !
    !        phi        = ANGPHI                                          !
    ! and returns the real and imaginary part in array SPHHAR             !
    !---------------------------------------------------------------------!
    Subroutine DEFSPH(LAMACT,MIUACT,COSTHE,ANGPHI,SPHHAR)

      Integer(ipr), INTENT(IN) :: LAMACT,MIUACT
      Real(pr), INTENT(IN) :: COSTHE,ANGPHI
      Real(pr), Dimension(1:2), INTENT(INOUT) :: SPHHAR

      Real(pr) :: QUNITS

      Integer(ipr) :: i
      Real(pr) :: ZLEGPO,PIARGU,FACMUL
      Real(pr), Allocatable :: FACTOR(:)

      COMMON                                    &
             /OURUNI/ QUNITS(0:NDMULT,0:NDMULT)

      ! Computing P_{l,m}(cos(theta))
      ZLEGPO=DEFLEG(LAMACT,MIUACT,COSTHE)

      ! Defining the factorials
      Allocate(FACTOR(0:LAMACT+MIUACT))
      FACTOR(0)=1.0_pr
      Do i=1,LAMACT+MIUACT
         FACTOR(I)=FACTOR(i-1)*Real(i,Kind=pr)
      End Do

      ! Computing the spherical harmonics
      PIARGU=4.0_pr*Atan(1.0_pr)

      FACMUL = Sqrt(FACTOR(LAMACT-MIUACT)/FACTOR(LAMACT+MIUACT)) &
             * Sqrt(0.25_pr*Real(2*LAMACT+1,Kind=pr)/PIARGU)     &
             * QUNITS(LAMACT,MIUACT)

      SPHHAR(1) = FACMUL*Cos(Real(MIUACT,Kind=pr)*ANGPHI)*ZLEGPO
      SPHHAR(2) = FACMUL*Sin(Real(MIUACT,Kind=pr)*ANGPHI)*ZLEGPO

      Deallocate(FACTOR)

    End Subroutine DEFSPH

    !---------------------------------------------------------------------!
    !    This subroutine computes:                                        !
    !                                                                     !
    !       SFACTO(MZ,NZ) = SUM_{K}^{NZ+MZ} C_{NZ,MZ}^{K}(00)             !
    !                                                                     !
    !                \INT_{ZPOINT}^{\INFTY} D\XI H_{K}(\XI) EXP(-\XI^2)   !
    !                                                                     !
    !    This is needed to evaluate the  expectation value of the         !
    !    density in the fission fragments. Inputs are NZMAXX, the         !
    !    maximum value for either nz or mz, and ZPOINT, the lower         !
    !    bound of the generalized integral                                !
    !---------------------------------------------------------------------!
    Subroutine DEFMAS(NZHERM,NZMAXX,ZPOINT,SFACTO)

      Integer(ipr), INTENT(IN) :: NZHERM,NZMAXX
      Real(pr), INTENT(IN) :: ZPOINT
      Real(pr), Allocatable, INTENT(INOUT) :: SFACTO(:,:)

      Real(pr) :: HERFAC,COEF00,COEF01,COEF11,COEF02,HOMSCA

      Integer(ipr) :: KARTEZ,NORDER,NZ,MZ,K,I
      Real(pr) :: PIARGU,VALINT,XPOINT,SUMVAL,F_INTE
      Real(pr), Dimension(1:ND2MAX+1) :: PHERMI,DHERMI,ZHERMI

      COMMON                                      &
             /FACHER/ HERFAC(0:ND2MAX)
      COMMON                                                       &
             /COEXYZ/ COEF00(0:ND2MAX,0:NDOSCI,0:NDOSCI,1:NDKART), &
                      COEF01(0:ND2MAX,0:NDOSCI,0:NDOSCI,1:NDKART), &
                      COEF11(0:ND2MAX,0:NDOSCI,0:NDOSCI,1:NDKART), &
                      COEF02(0:ND2MAX,0:NDOSCI,0:NDOSCI,1:NDKART)
      COMMON                                      &
             /SCALNG/ HOMSCA(NDKART)

      PIARGU=4.0_pr*Atan(1.0_pr); KARTEZ=3

      ! Computing the value of the probability distribution function
      VALINT=PROINT(ZPOINT)

      ! Compute value of Hermite polynomials up to order Max(Nz+Mz) at z=0
      XPOINT=0.0_pr; NORDER=2*NZMAXX
      Call D_HERM(XPOINT,NORDER,PHERMI,DHERMI,ND2MAX+1)

      ! Compute value of Hermite polynomials up to order Max(Nz+Mz) at z=Abs(ZPOINT)
      XPOINT=Abs(ZPOINT); NORDER=2*NZMAXX
      Call D_HERM(XPOINT,NORDER,ZHERMI,DHERMI,ND2MAX+1)

      ! Loop over bra and ket z quantum number
      Do NZ=0,NZMAXX
         Do MZ=0,NZMAXX

            ! Initial value of the sum for K = 0. Coefficients COEF00 are scaled by
            ! SQRT(HOMSCA(KARTEZ)) in DEVHER, we need to remove this scaling here
            SUMVAL=COEF00(0,MZ,NZ,KARTEZ)/Sqrt(HOMSCA(KARTEZ)) * 0.5_pr*Sqrt(Sqrt(PIARGU))*(1.0_pr - VALINT)
            If(ZPOINT.Lt.0.0_pr) SUMVAL=COEF00(0,MZ,NZ,KARTEZ)/Sqrt(HOMSCA(KARTEZ))* 0.5_pr*Sqrt(Sqrt(PIARGU))*(1.0_pr + VALINT)

            Do K=1,NZ+MZ

               ! Hermite polynomial of order k is stored in PHERMI(k+1), but this shift does
               ! *not* apply to its normalization coefficient, which is stored in HERFAC(k)
               ! as expected. Below, we need Hermite polynomials of order k-1 and normalization
               ! coefficients of order k
               If(Mod(K,2).Eq.1) Then
                  F_INTE = ZHERMI(K)*Exp(-ZPOINT*ZPOINT)/HERFAC(K)
               Else
                  F_INTE =(ZHERMI(K)*Exp(-ZPOINT*ZPOINT) - PHERMI(K)) / HERFAC(K)
                  If(ZPOINT.Lt.0.0_pr) F_INTE = - F_INTE
               End If

               SUMVAL=SUMVAL+COEF00(K,MZ,NZ,KARTEZ)/Sqrt(HOMSCA(KARTEZ)) * F_INTE

            End Do

            SFACTO(MZ,NZ)=SUMVAL

         End Do
      End Do

    End Subroutine DEFMAS

    !---------------------------------------------------------------------!
    !      This subroutine computes the expectation value of the total    !
    !      iso-scalar and charge density in the 'right fragment'.  The    !
    !      integration domain for the right fragment is defined by:       !
    !             X = [-\infty, +\infty]	                          !
    !             Y = [-\infty, +\infty]                                  !
    !             Z = [ ZPOINT, +\infty]                                  !
    !---------------------------------------------------------------------!
    Subroutine RFRAGM(NZHERM,NZMAXX,ZPOINT,RIGMAS,RIGCHA,SFACTO)

      Use MAD_PP

      Integer(ipr), INTENT(IN) :: NZHERM,NZMAXX
      Real(pr), INTENT(IN) :: ZPOINT
      Real(pr), Allocatable, INTENT(IN) :: SFACTO(:,:)
      Real(pr), INTENT(INOUT) :: RIGMAS,RIGCHA

      Integer(ipr) :: NXVECT,NYVECT,NZVECT,LDBASE,IPHAPP,IPHAPM,IPHAMP,IPHAMM

      Logical :: flip, asymmetry
      Integer(ipr) :: I,J,NX,NY,NZ,MX,MY,MZ,KARTEZ
      Real(pr), Allocatable :: D_MATR(:,:)
      Complex(pr) :: CRIGMAS,CRIGCHA,CLEFMAS,CLEFCHA

      COMMON                                              &
             /T_PHAS/ IPHAPP(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAPM(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAMP(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAMM(0:NDYMAX,0:NDYMAX,0:NDKART)
      COMMON                                                            &
             /BASISO/ NXVECT(1:NDBASE),NYVECT(1:NDBASE),NZVECT(1:NDBASE)
      COMMON                 &
             /DIMENS/ LDBASE

      CRIGMAS=Cmplx(0.0_pr,0.0_pr); CRIGCHA=Cmplx(0.0_pr,0.0_pr)

      ! Actual expectation value
      Do J=1,LDBASE
         NX=NXVECT(J); NY=NYVECT(J); NZ=NZVECT(J)
         Do I=1,LDBASE
            MX=NXVECT(I); MY=NYVECT(I); MZ=NZVECT(I)

            If(NX.Eq.MX.And.NY.Eq.MY) Then
               CRIGMAS=CRIGMAS + (DEN_PP(I,J,0,0)+DEN_PP(I,J,0,1)  &
                               +  DEN_PP(I,J,1,0)+DEN_PP(I,J,1,1)) &
                                 *SFACTO(MZ,NZ)*IPHAPP(MY,NY,0)
               CRIGCHA=CRIGCHA + (DEN_PP(I,J,0,1)+DEN_PP(I,J,1,1)) &
                                 *SFACTO(MZ,NZ)*IPHAPP(MY,NY,0)
            End If

         End Do
      End Do
      RIGMAS=Real(CRIGMAS); RIGCHA=Real(CRIGCHA)

      ! Debugging: Checking that the two alternative ways of computing the number of
      !            particles in the left and right fragments give comparable numbers.
      !            The direct way is by integrating the density from zN to +\infty and
      !            is the default mode of this routine; The indirect way is by computing
      !            the matrix of the corresponding operator in the HO basis (using
      !            Heaviside functions) and computing the trace with the density matrix
      !            in configuration space: this option is tested below.
      If(debug.Ge.2) Then
         Allocate(D_MATR(0:NZMAXX,0:NZMAXX))
         KARTEZ=3; flip = .False.; asymmetry = .False.
         Call define_asymmetry(NZHERM,NZMAXX,KARTEZ,D_MATR,ZPOINT,flip,asymmetry)
         CRIGMAS=Cmplx(0.0_pr,0.0_pr); CRIGCHA=Cmplx(0.0_pr,0.0_pr)
         Do J=1,LDBASE
            NX=NXVECT(J); NY=NYVECT(J); NZ=NZVECT(J)
            Do I=1,LDBASE
               MX=NXVECT(I); MY=NYVECT(I); MZ=NZVECT(I)
               If(NX.Eq.MX.And.NY.Eq.MY) Then
                  CRIGMAS=CRIGMAS + (DEN_PP(I,J,0,0)+DEN_PP(I,J,0,1)  &
                                  +  DEN_PP(I,J,1,0)+DEN_PP(I,J,1,1)) &
                                    *D_MATR(MZ,NZ)*IPHAPP(MY,NY,0)
                  CRIGCHA=CRIGCHA + (DEN_PP(I,J,0,1)+DEN_PP(I,J,1,1)) &
                                    *D_MATR(MZ,NZ)*IPHAPP(MY,NY,0)
               End If
            End Do
         End Do
         Write(6,'("In RFRAGM - RIGMAS = ",f20.14," RIGMAS_ = ",f20.14)') RIGMAS,Real(CRIGMAS)
         Write(6,'("In RFRAGM - RIGCHA = ",f20.14," RIGCHA_ = ",f20.14)') RIGCHA,Real(CRIGCHA)
         KARTEZ=3; flip = .True.; asymmetry = .False.
         Call define_asymmetry(NZHERM,NZMAXX,KARTEZ,D_MATR,ZPOINT,flip,asymmetry)
         CLEFMAS=Cmplx(0.0_pr,0.0_pr); CLEFCHA=Cmplx(0.0_pr,0.0_pr)
         Do J=1,LDBASE
            NX=NXVECT(J); NY=NYVECT(J); NZ=NZVECT(J)
            Do I=1,LDBASE
               MX=NXVECT(I); MY=NYVECT(I); MZ=NZVECT(I)
               If(NX.Eq.MX.And.NY.Eq.MY) Then
                  CLEFMAS=CLEFMAS + (DEN_PP(I,J,0,0)+DEN_PP(I,J,0,1)  &
                                  +  DEN_PP(I,J,1,0)+DEN_PP(I,J,1,1)) &
                                    *D_MATR(MZ,NZ)*IPHAPP(MY,NY,0)
                  CLEFCHA=CLEFCHA + (DEN_PP(I,J,0,1)+DEN_PP(I,J,1,1)) &
                                    *D_MATR(MZ,NZ)*IPHAPP(MY,NY,0)
               End If
            End Do
         End Do
         Write(6,'("In RFRAGM - LEFMAS_ = ",f20.14," Total = ",f20.14)') Real(CLEFMAS),Real(CLEFMAS+CRIGMAS)
         Write(6,'("In RFRAGM - LEFCHA_ = ",f20.14," Total = ",f20.14)') Real(CLEFCHA),Real(CLEFCHA+CRIGCHA)
         Deallocate(D_MATR)
      End If

    End Subroutine RFRAGM

    !---------------------------------------------------------------------!
    ! This function computes position of the center of mass of the left   !
    ! and right fragments. The resulting quantity is in fermis.           !
    !---------------------------------------------------------------------!
    Subroutine center_of_mass(NXHERM,NYHERM,NZHERM,Z_POSI,CENLEF,CENRIG)

      Integer(ipr), INTENT(IN) :: NXHERM,NYHERM,NZHERM
      Real(pr), INTENT(INOUT) :: Z_POSI,CENLEF,CENRIG

      Real(pr) :: DENSIT,DENCHA,EXPAUX,FOURWG,FOURPT,HOMSCA

      Integer(ipr) :: IX,IY,IZ,IERROR,i
      Real(pr) :: DERIV1,DERIVN,XBEGIN,XFINIS,XARGUM,X_STEP,RHOTMP,ZPOTMP,dh
      Real(pr) :: XPOINT,YPOINT,ZPOINT,DENLOC,W_HERM,RHONOR,Z_NORM
      Real(pr), Allocatable :: FUNCTI(:),GUNCTI(:),XVALUE(:),ZINTER(:),AUXSTO(:)
      Real(pr), Allocatable :: DINTER(:,:,:)

      COMMON                                        &
             /DENTOT/ DENSIT(NDXHRM,NDYHRM,NDZHRM), &
                      DENCHA(NDXHRM,NDYHRM,NDZHRM)
      COMMON                                        &
             /DENEXP/ EXPAUX(NDXHRM,NDYHRM,NDZHRM)
      COMMON                                        &
             /INTSTO/ FOURWG(1:NDGAUS,1:NDKART),    &
                      FOURPT(1:NDGAUS,1:NDKART)
      COMMON                                        &
             /SCALNG/ HOMSCA(NDKART)

      Allocate(FUNCTI(1:NZHERM),XVALUE(1:NZHERM))
      Allocate(AUXSTO(1:NZHERM),ZINTER(1:NPOINT))
      Allocate(DINTER(1:NZHERM,1:NZHERM,1:NPOINT))

      ! Integrations over the transverse coordinates are done using the dimensionless units
      ! xi_x and xi_y. In the longitudinal direction, the density is defined as a function
      ! of the dimensionless variable xi_z = bz * z. The operator O(z) = z must thus be
      ! defined on the grid z = xi_z/bz with integration limits given below.
      XBEGIN=FOURPT(NZHERM,3)/HOMSCA(3) + 1.D-14
      XFINIS=Z_POSI/HOMSCA(3)
      ! Boundary test: make sure there is a valid integration interval
      If(Z_POSI.LT.XBEGIN) Then
         Write(6,'("ATTENTION! In center_of_mass(): changing neck position from ",f20.14," to ",f20.14)') &
                    Z_POSI,XBEGIN + 1.D-12
         Z_POSI = XBEGIN + 1.D-12
      End If

      Do IY=1,NYHERM
         Do IX=1,NXHERM

            ! We define below the function at points z (=xi_z/bz) where the density is given by DENSIT.
            ! Dividing by the exponential factors is needed when using Gauss integrations, see routine
            ! test_density() for example.
            Do IZ=1,NZHERM
               XVALUE(IZ)=FOURPT(NZHERM-IZ+1,3)/HOMSCA(3)
               FUNCTI(IZ)=DENSIT(IX,IY,NZHERM-IZ+1)&
                         /Exp(-2.0_pr*FOURPT(IX,1)**2)/Exp(-2.0_pr*FOURPT(IY,2)**2)
            End Do

            ! Approximation of derivatives at the boundaries
            dh=XVALUE(2)-XVALUE(1) ! >0
            DERIV1=+0.5_pr*(-3.0_pr*FUNCTI(1)     +4.0_pr*FUNCTI(2)       -FUNCTI(3))       /dh
            DERIVN=-0.5_pr*(-3.0_pr*FUNCTI(NZHERM)+4.0_pr*FUNCTI(NZHERM-1)-FUNCTI(NZHERM-2))/dh

            ! Define spline coefficients
            IERROR=0; Call SPLINE(XVALUE,FUNCTI,NZHERM,DERIV1,DERIVN,AUXSTO,IERROR)

            ! Define interpolated mesh and value of the function of that mesh
            If(IERROR.Eq.0) Then
               Do i=1,NPOINT
                  XARGUM=XBEGIN+(XFINIS-XBEGIN)*Real(i-1,Kind=pr)/Real(NPOINT-1,Kind=pr)
                  ZINTER(i)=XARGUM
                  DINTER(IX,IY,i)=SPLINT(XVALUE,FUNCTI,AUXSTO,NZHERM,XARGUM)
               End Do
            Else
               Write(6,'("In center_of_mass() - Error in interpolating the function (right fragment)!")')
            End If

         End Do
      End Do

      X_STEP=(ZINTER(NPOINT)-ZINTER(1))/Real(NPOINT-1,Kind=pr)

      Deallocate(FUNCTI); Allocate(FUNCTI(1:NPOINT),GUNCTI(1:NPOINT))

      ! Computing the integrals \int \rho(z) and \int z\rho(z)
      !          -\infty < x < +\infty
      !          -\infty < y < +\infty
      !          -\infty < z < +z_{N}
      ! Integrations over x and y are performed "exactly" by Gauss-Hermite
      ! quadratures, integration over z is performed numerically using the
      ! Simpson 3/8 rule
      Do i=1,NPOINT
         FUNCTI(i)=0.0_pr; RHOTMP=0.0_pr
         GUNCTI(i)=0.0_pr; ZPOTMP=0.0_pr
         Do IY=1,NYHERM
            Do IX=1,NXHERM

               ZPOINT=ZINTER(i)

               DENLOC=DINTER(IX,IY,i)
               W_HERM=FOURWG(IX,1)*FOURWG(IY,2)

               RHOTMP=RHOTMP+W_HERM*DENLOC
               ZPOTMP=ZPOTMP+W_HERM*DENLOC*ZPOINT

            End Do
         End Do
         FUNCTI(i)=RHOTMP; GUNCTI(i)=ZPOTMP
      End Do

      ! Integrating over z, from -infty to zN
      RHONOR=0.0_pr; Call SIMP38(FUNCTI,NPOINT,X_STEP,RHONOR)
      Z_NORM=0.0_pr; Call SIMP38(GUNCTI,NPOINT,X_STEP,Z_NORM)

      CENLEF=Z_NORM/RHONOR

      Deallocate(FUNCTI,GUNCTI,XVALUE,DINTER,AUXSTO,ZINTER)

      ! Repeating this whole procedure for the right fragment
      Allocate(FUNCTI(1:NZHERM),XVALUE(1:NZHERM))
      Allocate(AUXSTO(1:NZHERM),ZINTER(1:NPOINT))
      Allocate(DINTER(1:NZHERM,1:NZHERM,1:NPOINT))

      XBEGIN=Z_POSI/HOMSCA(3)
      XFINIS=FOURPT(1,3)/HOMSCA(3) - 1.D-14
      If(Z_POSI.Gt.XFINIS) Then
         Write(6,'("ATTENTION! In center_of_mass(): changing neck position from ",f20.14," to ",f20.14)') &
                    Z_POSI,XFINIS - 1.D-12
         Z_POSI = XFINIS - 1.D-12
      End If

      Do IY=1,NYHERM
         Do IX=1,NXHERM

            Do IZ=1,NZHERM
               XVALUE(IZ)=FOURPT(NZHERM-IZ+1,3)/HOMSCA(3)
               FUNCTI(IZ)=DENSIT(IX,IY,NZHERM-IZ+1)&
                         /Exp(-2.0_pr*FOURPT(IX,1)**2)/Exp(-2.0_pr*FOURPT(IY,2)**2)
            End Do

            dh=XVALUE(2)-XVALUE(1)
            DERIV1=+0.5_pr*(-3.0_pr*FUNCTI(1)     +4.0_pr*FUNCTI(2)       -FUNCTI(3))       /dh
            DERIVN=-0.5_pr*(-3.0_pr*FUNCTI(NZHERM)+4.0_pr*FUNCTI(NZHERM-1)-FUNCTI(NZHERM-2))/dh

            IERROR=0; Call SPLINE(XVALUE,FUNCTI,NZHERM,DERIV1,DERIVN,AUXSTO,IERROR)

            If(IERROR.Eq.0) Then
               Do i=1,NPOINT
                  XARGUM=XBEGIN+(XFINIS-XBEGIN)*Real(i-1,Kind=pr)/Real(NPOINT-1,Kind=pr)
                  ZINTER(i)=XARGUM
                  DINTER(IX,IY,i)=SPLINT(XVALUE,FUNCTI,AUXSTO,NZHERM,XARGUM)
               End Do
            Else
               Write(6,'("In center_of_mass() - Error in interpolating the function (left fragment)!")')
            End If

         End Do
      End Do

      X_STEP=(ZINTER(NPOINT)-ZINTER(1))/Real(NPOINT-1,Kind=pr)

      Deallocate(FUNCTI); Allocate(FUNCTI(1:NPOINT),GUNCTI(1:NPOINT))

      Do i=1,NPOINT
         FUNCTI(i)=0.0_pr; RHOTMP=0.0_pr
         GUNCTI(i)=0.0_pr; ZPOTMP=0.0_pr
         Do IY=1,NYHERM
            Do IX=1,NXHERM

               ZPOINT=ZINTER(i)

               DENLOC=DINTER(IX,IY,i)
               W_HERM=FOURWG(IX,1)*FOURWG(IY,2)

               RHOTMP=RHOTMP+W_HERM*DENLOC
               ZPOTMP=ZPOTMP+W_HERM*DENLOC*ZPOINT

            End Do
         End Do
         FUNCTI(i)=RHOTMP; GUNCTI(i)=ZPOTMP
      End Do

      RHONOR=0.0_pr; Call SIMP38(FUNCTI,NPOINT,X_STEP,RHONOR)
      Z_NORM=0.0_pr; Call SIMP38(GUNCTI,NPOINT,X_STEP,Z_NORM)

      CENRIG=Z_NORM/RHONOR

      Deallocate(FUNCTI,GUNCTI,AUXSTO,ZINTER,DINTER,XVALUE)

    End Subroutine center_of_mass

    !---------------------------------------------------------------------!
    !  The routine computes the matrix elements of the Gaussian neck      !
    !  operator qn = exp(-(z - zn)^2/an^2), with an = 1 fm and zn is      !
    !  the  position of the  neck (defined as the  point  inside the      !
    !  nucleus where the  density  is the lowest). The neck operator      !
    !  is  real,  time-invariant and  simplex-invariant. It has only      !
    !  non-zero matrix elements in the (t+,t+) block.                     !
    !---------------------------------------------------------------------!
    Subroutine INT_QN(NZMAXX,HAUXPP)

      Real(pr), Allocatable, INTENT(INOUT) :: HAUXPP(:,:)
      Integer(ipr), INTENT(IN) :: NZMAXX

      Integer(ipr) :: NXVECT,NYVECT,NZVECT,LDBASE,IPHAPP,IPHAPM,IPHAMP,IPHAMM

      Integer(ipr) :: nx,ny,nz,mx,my,mz,i,j,KARTEZ
      Real(pr), Allocatable :: QNMATR(:,:)

      COMMON                                              &
             /T_PHAS/ IPHAPP(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAPM(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAMP(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAMM(0:NDYMAX,0:NDYMAX,0:NDKART)
      COMMON                                                             &
             /BASISO/ NXVECT(1:NDBASE),NYVECT(1:NDBASE),NZVECT(1:NDBASE)
      COMMON                 &
             /DIMENS/ LDBASE

      ! Computing the integral of a Gaussian multiplied by a Hermite
      ! polynomials, see:   Table of Integrals, Series, and Products
      !                     I. S. Gradshteyn / I. M. Ryzhik
      !                     Eq. (8), Sec. 7.374, Ch. 7, P. 837
      Allocate(QNMATR(0:NZMAXX,0:NZMAXX)); KARTEZ=3
      Call DEF_QN(NZMAXX,KARTEZ,QNMATR)

      ! The matrix elements of QN are simply the product of such integrals
      HAUXPP(:,:)=0.0_pr
      Do j=1,LDBASE
         nx=NXVECT(j); ny=NYVECT(j); nz=NZVECT(j)
         Do i=1,LDBASE
            mx=NXVECT(i); my=NYVECT(i); mz=NZVECT(i)

            If(mx.Eq.nx.And.my.Eq.ny) Then
               HAUXPP(i,j)=QNMATR(mz,nz)*IPHAPP(my,ny,0)
            End If

         End Do
      End Do

      Deallocate(QNMATR)

    End Subroutine INT_QN

    !---------------------------------------------------------------------!
    ! The  routine  computes the integral of a  Gaussian multiplied by    !
    ! a Hermite polynomial. It needs to evaluate the  HO  wavefunction    !
    ! at the neck, the coordinates of which are X_NECK, Y_NECK, Z_NECK    !
    !                                                                     !
    ! Ref.: Table of integrals, series, and products                      !
    !       I. S. Gradshteyn / I. M. Ryzhik                               !
    !       Eq. (8), sec. 7.374, ch. 7, p. 837                            !
    !                                                                     !
    ! Attention: coefficients  COEF00()  are scaled by the factors        !
    !            SQRT(HOMSCA(KARTEZ)) in subroutine DEVHER(). This        !
    !            scaling is removed here.                                 !
    !---------------------------------------------------------------------!
    Subroutine DEF_QN(N_MAXX,KARTEZ,QNMATR)

      Integer(ipr), INTENT(IN) :: N_MAXX,KARTEZ
      Real(pr), Allocatable, INTENT(INOUT) :: QNMATR(:,:)

      Real(pr) :: COEF00,COEF01,COEF11,COEF02,HERFAC,HOMSCA

      Integer(ipr) :: k,n,m
      Real(pr) :: PIARGU,G_KART,X_KART,COORDI,ADDEXP,SUM_QN
      Real(pr), Dimension(1:ND2MAX+1) :: HERVAL,DHERVA
      Real(pr), Allocatable :: WAVFUN(:),GFACTO(:)

      COMMON                                                       &
             /COEXYZ/ COEF00(0:ND2MAX,0:NDOSCI,0:NDOSCI,1:NDKART), &
                      COEF01(0:ND2MAX,0:NDOSCI,0:NDOSCI,1:NDKART), &
                      COEF11(0:ND2MAX,0:NDOSCI,0:NDOSCI,1:NDKART), &
                      COEF02(0:ND2MAX,0:NDOSCI,0:NDOSCI,1:NDKART)
      COMMON                           &
             /FACHER/ HERFAC(0:ND2MAX)
      COMMON                         &
             /SCALNG/ HOMSCA(NDKART)

      PIARGU=4.0_pr*Atan(1.0_pr); G_KART=1.0_pr+(HOMSCA(KARTEZ)*AN_VAL)**2
      X_KART=AN_VAL*Sqrt(PIARGU*HOMSCA(KARTEZ)/G_KART)

                      COORDI=0.0_pr
      If(KARTEZ.Eq.1) COORDI=X_NECK/Sqrt(G_KART)
      If(KARTEZ.Eq.2) COORDI=Y_NECK/Sqrt(G_KART)
      If(KARTEZ.Eq.3) COORDI=Z_NECK/Sqrt(G_KART)

      ADDEXP=Exp(-0.5_pr*COORDI**2)

      Call D_HERM(COORDI,2*N_MAXX,HERVAL,DHERVA,ND2MAX+1)

      Allocate(WAVFUN(0:2*N_MAXX)); Allocate(GFACTO(0:2*N_MAXX))
      Do K=0,2*N_MAXX
         WAVFUN(K)=Sqrt(HOMSCA(KARTEZ))*ADDEXP*HERVAL(K+1)/HERFAC(K)
         GFACTO(K)=G_KART**(K/2.0_pr)
      End Do

      Do n=0,N_MAXX
         Do m=0,N_MAXX
            QNMATR(m,n)=0.0_pr; SUM_QN=0.0_pr
            Do k=0,n+m
               SUM_QN=SUM_QN+COEF00(k,m,n,KARTEZ)*WAVFUN(k)/GFACTO(k)
            End Do
            QNMATR(m,n)=SUM_QN*X_KART*ADDEXP/Sqrt(HOMSCA(KARTEZ))
         End Do
      End Do

      Deallocate(WAVFUN,GFACTO)

    End Subroutine DEF_QN

    !---------------------------------------------------------------------!
    !  The routine finds the position of the neck, where <QN> is the      !
    !  lowest.  At the minimum, all partial derivatives of <QN> with      !
    !  respect to  xn, yn, zn (coordinates of the neck) should be 0.      !
    !  We find these 0 by  three  successive  Newton methods for the      !
    !  three coordinates.  The starting point for all Newton methods      !
    !  is the position of the neck at the previous HF iteration.          !
    !---------------------------------------------------------------------!
    Subroutine QNFIND(NXHERM,NYHERM,NZHERM,SLOWOD,NUMITE,ISIMPY,ISIGNY,IPARTY)

      Integer(ipr), INTENT(IN) :: NXHERM,NYHERM,NZHERM,NUMITE,ISIMPY,ISIGNY,IPARTY
      Real(pr), INTENT(IN) :: SLOWOD

      Integer(ipr) :: IFNECK,NFIPRI

      Logical :: symmetrical,interpolate,enforced_symmetry,any_odd_multipole
      Integer(ipr) :: i,l,m,Ninter,IERROR,IFOUND,N_SCAN,KARTEZ,KORDER,found
      Real(pr) :: QN_VAL,QBEFOR,ZBEFOR,DQNPRE,DQNCUR,XLOWER,XUPPER
      Real(pr) :: dq,z_optim,Z_INIT_back,Qtotal
      Real(pr) :: Z_MINI,Z_MAXI,X_INIT,Y_INIT,Z_INIT,Q_INIT,Z_POSI,TOLERA
      Real(pr) :: XARGUM,XBEGIN,XFINIS,DERIV1,DERIVN,dh
      Real(pr) :: QMUL_N,QMUL_P,QMUL_T
      Real(pr), Allocatable :: FUNCTI(:),XVALUE(:),AUXSTO(:),ZINTER(:),DINTER(:)

      COMMON                                           &
             /QMULTI/ QMUL_N(0:NDMULT,-NDMULT:NDMULT), &
                      QMUL_P(0:NDMULT,-NDMULT:NDMULT), &
                      QMUL_T(0:NDMULT,-NDMULT:NDMULT)
      COMMON                 &
             /NCKFLA/ IFNECK &
             /CFIPRI/ NFIPRI

      ! Initialization of the neck coordinates from previously found solutions.
      ! At the very first iteration, (X_NECK,Y_NECK,Z_NECK)=(0,0,0)
      X_INIT=X_NECK; Y_INIT=Y_NECK; Z_INIT=Z_NECK; Q_INIT=0.0_pr
      MXHERM=NXHERM; MYHERM=NYHERM; MZHERM=NZHERM

      symmetrical = .True.

      ! Test if symmetry flags enforce symmetric nucleus
         enforced_symmetry = .True.
      If(ISIMPY.Eq.1.And.(ISIGNY.Eq.0.Or.IPARTY.Eq.0)) Then
         enforced_symmetry = .False.
      End If
      ! Test if all multipole moments are actually zero
      any_odd_multipole = .False.
      Do l=1,NDMULT,2
         Do m=0,l
            Qtotal = QMUL_N(l,m) + QMUL_P(l,m)
            If(Abs(Qtotal).Gt.1.e-4) Then
               any_odd_multipole =.True.
            End If
         End Do
      End Do
      ! System is symmetric either if symmetry flags enforce it, or if self-consistency makes it so
      If(.Not.enforced_symmetry) Then
         symmetrical = .Not.any_odd_multipole
      Else
         symmetrical = enforced_symmetry
      End If

      KARTEZ=3; KORDER=0

      ! Initial scanning of the values of QN: we need this to get a somewhat
      ! reliable initial estimate of the position of the neck (only if
      ! reflection symmetry is broken). We search around the last known
      ! position of the neck.
      If(.Not.symmetrical) Then

         DQNPRE=0.1_pr; N_SCAN=60
         Z_MINI=-5.0_pr; Z_MAXI=+5.0_pr
         dq = (Z_MAXI-Z_MINI)/Real(N_SCAN-1,Kind=pr)
         ZBEFOR=Z_MINI; QBEFOR=QMOM_K(KARTEZ,ZBEFOR-dq,KORDER)

            interpolate=.True.
         If(interpolate) Then

            Ninter=3*N_SCAN

            Allocate(FUNCTI(1:N_SCAN),XVALUE(1:N_SCAN),AUXSTO(1:N_SCAN))
            Allocate(ZINTER(1:Ninter),DINTER(1:Ninter))
            Do i=1,N_SCAN
               XVALUE(i)=dq*Real(i-1,Kind=pr)
               FUNCTI(i)=QMOM_K(KARTEZ,XVALUE(i),KORDER)
            End Do
            ! Approximation of derivatives at the boundaries
            DERIV1=+0.5_pr*(-3.0_pr*FUNCTI(1)     +4.0_pr*FUNCTI(2)       -FUNCTI(3))       /dq
            DERIVN=-0.5_pr*(-3.0_pr*FUNCTI(N_SCAN)+4.0_pr*FUNCTI(N_SCAN-1)-FUNCTI(N_SCAN-2))/dq

            ! Define spline coefficients
             IERROR=0; Call SPLINE(XVALUE,FUNCTI,N_SCAN,DERIV1,DERIVN,AUXSTO,IERROR)

            ! Define interpolated mesh and value of QN and its the derivative
            If(IERROR.Eq.0) Then
               Do i=1,Ninter
                  XARGUM=Z_MINI+(Z_MAXI-Z_MINI)*Real(i-1,Kind=pr)/Real(Ninter-1,Kind=pr)
                  ZINTER(i)=SPLINT(XVALUE,FUNCTI,AUXSTO,N_SCAN,XARGUM)
                  DINTER(i)=SPLIND(XVALUE,FUNCTI,AUXSTO,N_SCAN,XARGUM)
                  If(debug.Ge.2) Then
                     Write(6,'("SPLINE - zint=",f20.14," QN=",f20.14," dQN=",f20.14)') &
                                XARGUM,ZINTER(i),DINTER(i)
                  End If
               End Do
            End If
         End If

         found = 0
         Do i=1,N_SCAN
            Z_POSI=Z_MINI+dq*Real(i-1,Kind=pr)
            QN_VAL=QMOM_K(KARTEZ,Z_POSI,KORDER)
            DQNCUR=DERIVE(Z_POSI)
            If(debug.Ge.2) Then
               Write(6,'("z=",f20.14," QN=",f20.14," dQN(analytic)=",f20.14," dQN(numerical)=",f20.14)') &
                          Z_POSI,QN_VAL,DQNCUR,(QN_VAL-QBEFOR)/dq
            End If
            ! Found the minimum
            If(DQNCUR.Gt.0.0_pr.And.DQNPRE.Lt.0.0_pr) Then
               Z_INIT=0.5_pr*(Z_POSI+ZBEFOR)
               ZNMINI=Z_INIT-0.5_pr*dq; ZNMAXI=Z_INIT+0.5_pr*dq
               Q_INIT=QN_VAL
               found =found + 1
            End If
            ZBEFOR=Z_POSI
            QBEFOR=QN_VAL
            DQNPRE=DQNCUR
         End Do

         If(debug.Ge.1) Then
            Write(6,'("NUMITE=",i4," Z_NECK=",f20.14)') NUMITE,Z_NECK
            Write(6,'("ZNMINI=",f20.14," ZNMAXI=",f20.14)') ZNMINI,ZNMAXI
            Write(6,'("found=",i4," Z_POSI=",f20.14," QN=",f20.14)') found,Z_INIT,Q_INIT
         End If

         ! Search for the z-coordinate of the neck (most likely to be
         ! different from 0)
         If(found.Gt.0) Then
            Z_INIT_back = Z_INIT ! backup of rough estimate
            XLOWER=ZNMINI; XUPPER=ZNMAXI; TOLERA=1.D-5; IFOUND=0
            Z_INIT=ZBRENT(DERIVE,XLOWER,XUPPER,TOLERA,IFOUND)
            If(debug.Ge.1) Then
               Write(6,'("found (Newton)=",i4," Z_INIT=",f20.14)') ifound,Z_INIT
            End If
            ! If ZBRENT gives a value too different from the rough estimate, there is a problem...
            If(Abs(Z_INIT-Z_INIT_back).Gt.0.5_pr*Abs(dq)) Then
               If(debug.Ge.1) Then
                  Write(6,'("Initial guess: ",f20.14," ZBRENT found: ", F20.14)') Z_INIT_back, Z_INIT
               End If
               Z_INIT = Z_INIT_back
            End If
            If(IFOUND.EQ.0) Z_INIT = Z_INIT_back
         Else
            ! If no candidate for the neck position was found, we default to 0
            Z_INIT=0.0_pr
         End If

         ! Updating the neck coordinates
         If(NUMITE.Eq.0.Or.IFNECK.Eq.0) Then
            X_NECK=X_INIT; Y_NECK=Y_INIT; Z_NECK=Z_INIT
         Else
            X_NECK=SLOWOD*X_NECK+(1.0_pr-SLOWOD)*X_INIT
            Y_NECK=SLOWOD*Y_NECK+(1.0_pr-SLOWOD)*Y_INIT
            Z_NECK=SLOWOD*Z_NECK+(1.0_pr-SLOWOD)*Z_INIT
         End If

         Q_NECK=QMOM_K(KARTEZ,Z_NECK,KORDER)

      Else
         ! Reflection symmetry is conserved, the neck is at z=0, we compute QN
         X_NECK=0.0_pr; Y_NECK=0.0_pr; Z_NECK=0.0_pr
         Q_NECK=QMOM_K(KARTEZ,Z_NECK,KORDER)
      End If

    End Subroutine QNFIND

    !---------------------------------------------------------------------!
    !  Function computes the first derivative of the gaussian neck        !
    !  expectation value with respect to the z-coordinate of the neck.    !
    !  Note that the coordinate Z_POSI is dimensionless here (equivalent  !
    !  to \chi_N in the notes).                                           !
    !---------------------------------------------------------------------!
    Real(pr) Function DERIVE(Z_POSI)

      Real(pr), INTENT(IN) :: Z_POSI

      Real(pr) :: HOMSCA

      Integer(ipr) :: KORDER,KARTEZ
      Real(pr) :: DERIV1,PREFAC

      COMMON                                &
             /SCALNG/ HOMSCA(NDKART)

      DERIV1=0.0_pr; KARTEZ=3

      PREFAC=2.0_pr/(HOMSCA(KARTEZ)*AN_VAL)**2

      KORDER=1; DERIV1=QMOM_K(KARTEZ,Z_POSI,KORDER)
      KORDER=0; DERIV1=DERIV1-Z_POSI*QMOM_K(KARTEZ,Z_POSI,KORDER)

      DERIVE=PREFAC*DERIV1

    End Function DERIVE

    !---------------------------------------------------------------------!
    !  Function computes moments i_k of order k of the gaussian neck      !
    !---------------------------------------------------------------------!
    Real(pr) Function QMOM_K(KARTEZ,Z_POSI,KORDER)

      Integer(ipr), INTENT(IN) :: KARTEZ,KORDER
      Real(pr), INTENT(IN) :: Z_POSI

      Real(pr) :: DENSIT,DENCHA,EXPAUX,FOURWG,FOURPT,HOMSCA

      Integer(ipr) :: kx,ky,kz

      Real(pr) :: QN_VAL,W_HERM,DENLOC,ZNECKZ,COORDO,PREFAC

      COMMON                                        &
             /DENTOT/ DENSIT(NDXHRM,NDYHRM,NDZHRM), &
                      DENCHA(NDXHRM,NDYHRM,NDZHRM)
      COMMON                                        &
             /DENEXP/ EXPAUX(NDXHRM,NDYHRM,NDZHRM)
      COMMON                                        &
             /INTSTO/ FOURWG(1:NDGAUS,1:NDKART),    &
                      FOURPT(1:NDGAUS,1:NDKART)
      COMMON                                        &
             /SCALNG/ HOMSCA(NDKART)

      QN_VAL=0.0_pr
      Do kz=1,MZHERM
         Do ky=1,MYHERM
            Do kx=1,MXHERM
               DENLOC=DENSIT(kx,ky,kz)/EXPAUX(kx,ky,kz)/EXPAUX(kx,ky,kz)
               W_HERM=FOURWG(kx,1)*FOURWG(ky,2)*FOURWG(kz,3)
               ZNECKZ=Exp(-((FOURPT(kz,KARTEZ)-Z_POSI)/(HOMSCA(KARTEZ)*AN_VAL))**2)
               COORDO=FOURPT(kz,KARTEZ)
               If(KORDER.Eq.0) Then
                  PREFAC=1.0_pr
               Else
                  PREFAC=COORDO**KORDER
               End If
               QN_VAL=QN_VAL+W_HERM*DENLOC*ZNECKZ*PREFAC
            End Do
         End Do
      End Do

      QMOM_K=QN_VAL

    End Function QMOM_K

    !---------------------------------------------------------------------!
    !  Function computes moments i_k of order k of the gaussian neck      !
    !---------------------------------------------------------------------!
    Real(pr) Function protons_in_neck(KARTEZ,Z_POSI,KORDER)

      Integer(ipr), INTENT(IN) :: KARTEZ,KORDER
      Real(pr), INTENT(IN) :: Z_POSI

      Real(pr) :: DENSIT,DENCHA,EXPAUX,FOURWG,FOURPT,HOMSCA

      Integer(ipr) :: kx,ky,kz

      Real(pr) :: QN_VAL,W_HERM,DENLOC,ZNECKZ,COORDO,PREFAC

      COMMON                                        &
             /DENTOT/ DENSIT(NDXHRM,NDYHRM,NDZHRM), &
                      DENCHA(NDXHRM,NDYHRM,NDZHRM)
      COMMON                                        &
             /DENEXP/ EXPAUX(NDXHRM,NDYHRM,NDZHRM)
      COMMON                                        &
             /INTSTO/ FOURWG(1:NDGAUS,1:NDKART),    &
                      FOURPT(1:NDGAUS,1:NDKART)
      COMMON                                        &
             /SCALNG/ HOMSCA(NDKART)

      QN_VAL=0.0_pr
      Do kz=1,MZHERM
         Do ky=1,MYHERM
            Do kx=1,MXHERM
               DENLOC=DENCHA(kx,ky,kz)/EXPAUX(kx,ky,kz)/EXPAUX(kx,ky,kz)
               W_HERM=FOURWG(kx,1)*FOURWG(ky,2)*FOURWG(kz,3)
               ZNECKZ=Exp(-((FOURPT(kz,KARTEZ)-Z_POSI)/(HOMSCA(KARTEZ)*AN_VAL))**2)
               COORDO=FOURPT(kz,KARTEZ)
               If(KORDER.Eq.0) Then
                  PREFAC=1.0_pr
               Else
                  PREFAC=COORDO**KORDER
               End If
               QN_VAL=QN_VAL+W_HERM*DENLOC*ZNECKZ*PREFAC
            End Do
         End Do
      End Do

      protons_in_neck = QN_VAL

    End Function protons_in_neck

    !---------------------------------------------------------------------!
    !  The function computes the expectation value of the operator for    !
    !  the distance between the two pre-fragments. The form of the opera- !
    !  tor is:                                                            !
    !       D(z) = | 1/A_R * z * H(z-zN) - 1/A_L * z* [1 - H(z-zN)] |     !
    !                                                                     !
    !  where H(x) is the Heaviside step function and zN the z-coordinate  !
    !  of the neck.                                                       !
    !---------------------------------------------------------------------!
    Subroutine distance_value()

      Use MAD_PP
      Use NEWCNT

      Integer(ipr) :: i,j,LDBASE
      Complex(pr) :: distance

      COMMON &
             /DIMENS/ LDBASE

      distance = Cmplx(0.0_pr,0.0_pr)
      Do j=1,LDBASE
         Do i=1,LDBASE
            distance = distance + Cmplx(ARDIST(i,j),0.0_pr)        &
                                *(DEN_PP(j,i,0,0)+DEN_PP(j,i,1,0)  &
                                + DEN_PP(j,i,0,1)+DEN_PP(j,i,1,1))
         End Do
      End Do

      distance_prefragments = Real(distance,Kind=pr)

    End Subroutine distance_value

    !---------------------------------------------------------------------!
    !  The routine computes the matrix elements of the operator for the   !
    !  distance between the two pre-fragments. The form of the operator   !
    !  is:                                                                !
    !       D(z) = | 1/A_R * z * H(z-zN) - 1/A_L * z* [1 - H(z-zN)] |     !
    !                                                                     !
    !  where H(x) is the Heaviside step function and zN the z-coordinate  !
    !  of the neck. The distance operator is real, time-invariant and     !
    !  simplex-invariant. It has only non-zero matrix elements in the     !
    !  (t+,t+) block.                                                     !
    !---------------------------------------------------------------------!
    Subroutine distance_operator(IN_FIX,IZ_FIX,NXHERM,NYHERM,NZHERM,NZMAXX, &
                                 ITERUN,ISIMPY,ISIGNY,IPARTY,SLOWOD)

      Use NEWCNT

      Integer(ipr), INTENT(IN) :: IN_FIX,IZ_FIX,NXHERM,NYHERM,NZHERM,NZMAXX,ITERUN,ISIMPY,ISIGNY,IPARTY
      Real(pr), INTENT(IN) :: SLOWOD

      Integer(ipr) :: NXVECT,NYVECT,NZVECT,LDBASE,IPHAPP,IPHAPM,IPHAMP,IPHAMM

      Integer(ipr) :: KARTEZ,i,j,nx,ny,nz,mx,my,mz,I_TYPE,LAMBDA,MIU
      Real(pr) :: A_LEFT,A_RIGH,CEFLEF,RIGCHA,mixing
      Real(pr), Allocatable :: HAUXPP(:,:),D_MATR(:,:)

      COMMON                                              &
             /T_PHAS/ IPHAPP(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAPM(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAMP(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAMM(0:NDYMAX,0:NDYMAX,0:NDKART)
      COMMON                                                             &
             /BASISO/ NXVECT(1:NDBASE),NYVECT(1:NDBASE),NZVECT(1:NDBASE)
      COMMON                 &
             /DIMENS/ LDBASE

      If(.Not.Allocated(ARDIST)) Then
         Allocate(ARDIST(1:NDBASE,1:NDBASE))
         ARDIST(:,:) = 0.0_pr
      End If

      If(distance_is_calculated) Then
         mixing = SLOWOD
      Else
         mixing = 0.0_pr
      End If

     ! Computing the mass and charge of the fragments. Routine QLMFRA relies on
      ! density on Gauss-Hermite mesh, while RFRAGM requires the density in the
      ! HO basis
      I_TYPE=1; LAMBDA=0; MIU=0; CEFLEF=0.0_pr; A_LEFT=0.0_pr; A_RIGH=0.0_pr
      CALL QLMFRA(NXHERM,NYHERM,NZHERM,Z_NECK,LAMBDA,MIU,A_LEFT,A_RIGH,CEFLEF,CEFLEF,I_TYPE)

      mass_left =(1.0_pr-mixing)*A_LEFT + mixing*mass_left
      mass_right=(1.0_pr-mixing)*A_RIGH + mixing*mass_right

      ! Defining the distance between fragments
      Allocate(D_MATR(0:NZMAXX,0:NZMAXX)); D_MATR(:,:)=0.0_pr
      KARTEZ=3
      Call define_distance(NZHERM,NZMAXX,KARTEZ,D_MATR,Z_NECK)

      Allocate(HAUXPP(1:NDBASE,1:NDBASE))
      HAUXPP(:,:)=0.0_pr
      Do j=1,LDBASE
         nx=NXVECT(j); ny=NYVECT(j); nz=NZVECT(j)
         Do i=1,LDBASE
            mx=NXVECT(i); my=NYVECT(i); mz=NZVECT(i)

            If(mx.Eq.nx.And.my.Eq.ny) Then
               HAUXPP(i,j)=D_MATR(mz,nz)*IPHAPP(my,ny,0)
            End If

         End Do
      End Do

      Deallocate(D_MATR)

      ! Linear mixing of matrix elements
      Do j=1,LDBASE
         Do i=1,LDBASE
            ARDIST(i,j)=(1.0_pr-mixing)*HAUXPP(i,j)+mixing*ARDIST(i,j)
         End Do
      End Do

      Deallocate(HAUXPP)

    End Subroutine distance_operator

    !---------------------------------------------------------------------!
    !  The routine computes the integral of the distance operator along   !
    !  the z-direction:                                                   !
    !                                                                     !
    !      D_mn = \sum_k C_mn^k(00) \int dxi H_k(xi) e^(-xi^2) D(z)       !
    !                                                                     !
    !  where xi = b_z * z, H_k(xi) is the (normalized) Hermite polynomial !
    !  used throughout in HFODD, and C_mn^k(00) relates the product of 2  !
    !  Hermite polynomials to a sum of HErmite polynomials. The output is !
    !  the matrix D(m_z, n_z) along the z-direction only. Integrations    !
    !  are performed by Gauss-Hermite quadrature or numericall by using   !
    !  Simpson's 3/8 rule.                                                !
    !---------------------------------------------------------------------!
    Subroutine define_distance(NZHERM,N_MAXX,KARTEZ,D_MATR,Z_POSI)

      Integer(ipr), INTENT(IN) :: NZHERM,N_MAXX,KARTEZ
      Real(pr), INTENT(IN) :: Z_POSI
      Real(pr), Allocatable, INTENT(INOUT) :: D_MATR(:,:)

      Real(pr) :: COEF00,COEF01,COEF11,COEF02,HOMSCA,TWOWGT,TWOPNT,HERONE,A_LEFT,A_RIGH

      Integer(ipr) :: k,n,m,IERROR,i,iz,Nx
      Real(pr) :: SUM_DN,COORDI,operator
      Real(pr) :: DERIV1,DERIVN,XBEGIN,XFINIS,XARGUM,X_STEP,dh
      Real(pr), Allocatable :: WAVFUN(:,:),value(:)
      Real(pr), Allocatable :: XVALUE(:),FUNCTI(:),ZINTER(:),AUXSTO(:),DINTER(:)

      COMMON                                            &
             /HERMEM/ HERONE(0:ND2MAX,1:NDGAUS,1:NDKART)
      COMMON                                                       &
             /COEXYZ/ COEF00(0:ND2MAX,0:NDOSCI,0:NDOSCI,1:NDKART), &
                      COEF01(0:ND2MAX,0:NDOSCI,0:NDOSCI,1:NDKART), &
                      COEF11(0:ND2MAX,0:NDOSCI,0:NDOSCI,1:NDKART), &
                      COEF02(0:ND2MAX,0:NDOSCI,0:NDOSCI,1:NDKART)
      COMMON                         &
             /SCALNG/ HOMSCA(NDKART)
      COMMON                                     &
             /INTMEM/ TWOWGT(1:NDGAUS,1:NDKART), &
                      TWOPNT(1:NDGAUS,1:NDKART)

      Allocate(value(0:2*N_MAXX))
      Allocate(WAVFUN(1:NZHERM,0:2*N_MAXX)); WAVFUN(:,:)=0.0_pr

      A_LEFT = mass_left
      A_RIGH = mass_right

      Nx = 1000000

      ! See comments in routine define_asymmetry()
      Do iz=1,NZHERM
         COORDI = TWOPNT(NZHERM-iz+1,KARTEZ)
         Do k=0,2*N_MAXX
         WAVFUN(iz,k) = HERONE(k,NZHERM-iz+1,KARTEZ)/Sqrt(HOMSCA(KARTEZ)) * Exp(-COORDI*COORDI)
         End Do
      End Do

      ! Define boundaries
      XBEGIN=TWOPNT(NZHERM,KARTEZ) + 1.0D-14
      XFINIS=TWOPNT(1,KARTEZ) - 1.0D-14

!$OMP PARALLEL DO &
!$OMP DEFAULT(NONE) &
!$OMP SCHEDULE(STATIC) &
!$OMP SHARED(N_MAXX,NZHERM,TWOPNT,KARTEZ,WAVFUN,Nx,XBEGIN,XFINIS,A_LEFT,A_RIGH,&
!$OMP        HOMSCA,value,Z_POSI,prec_heaviside) &
!$OMP PRIVATE(k,iz,XVALUE,FUNCTI,dh,DERIV1,DERIVN,IERROR,AUXSTO,i,XARGUM,&
!$OMP         operator,ZINTER,DINTER,X_STEP,SUM_DN)
      Do k=0,2*N_MAXX

         Allocate(FUNCTI(1:NZHERM),XVALUE(1:NZHERM),AUXSTO(1:NZHERM))
         Allocate(DINTER(1:Nx),ZINTER(1:Nx))

         ! Define function to interpolate
         Do iz=1,NZHERM
            XVALUE(iz)=TWOPNT(NZHERM-iz+1,KARTEZ)
            FUNCTI(iz)=WAVFUN(iz,k)
         End Do

         ! Approximation of derivatives at the boundaries
         dh=XVALUE(2)-XVALUE(1) ! >0
         DERIV1=+0.5_pr*(-3.0_pr*FUNCTI(1)     +4.0_pr*FUNCTI(2)       -FUNCTI(3))       /dh
         DERIVN=-0.5_pr*(-3.0_pr*FUNCTI(NZHERM)+4.0_pr*FUNCTI(NZHERM-1)-FUNCTI(NZHERM-2))/dh

         ! Define spline coefficients
         IERROR=0; Call SPLINE(XVALUE,FUNCTI,NZHERM,DERIV1,DERIVN,AUXSTO,IERROR)

         ! Define interpolated mesh and value of the function of that mesh
         If(IERROR.Eq.0) Then
            Do i=1,Nx
               XARGUM=XBEGIN+(XFINIS-XBEGIN)*Real(i-1,Kind=pr)/Real(Nx-1,Kind=pr)
               operator =XARGUM * (Heaviside(XARGUM-Z_POSI,prec_heaviside)/A_RIGH &
                         - (1.0_pr-Heaviside(XARGUM-Z_POSI,prec_heaviside))/A_LEFT) / HOMSCA(KARTEZ)
               ZINTER(i)=XARGUM
               DINTER(i)=operator*SPLINT(XVALUE,FUNCTI,AUXSTO,NZHERM,XARGUM)
            End Do
         Else
            Write(6,'("In define_distance() - Error in interpolating the function (left fragment)!")')
         End If

         ! Define step for integration
         X_STEP=(ZINTER(Nx)-ZINTER(1))/Real(Nx-1,Kind=pr)

         ! Integrate
         SUM_DN=0.0_pr; Call SIMP38(DINTER,Nx,X_STEP,SUM_DN)
         value(k)=SUM_DN

         Deallocate(FUNCTI,XVALUE,AUXSTO,DINTER,ZINTER)

      End Do
!$OMP END PARALLEL DO

      ! Filling the matrix. Attention: Coefficients COEF00 are scaled by
      ! SQRT(HOMSCA(KARTEZ)) in DEVHER, we need to remove this scaling here
      Do n=0,N_MAXX
         Do m=0,N_MAXX
            SUM_DN=0.0_pr
            Do k=0,n+m
               SUM_DN=SUM_DN+COEF00(k,m,n,KARTEZ)/Sqrt(HOMSCA(KARTEZ))*value(k)
            End Do
            D_MATR(m,n)=SUM_DN
         End Do
      End Do

      Deallocate(WAVFUN)
      Deallocate(value)

    End Subroutine define_distance

    !---------------------------------------------------------------------!
    !  The function computes the expectation value of the operator for    !
    !  the mass asymmetry between the two pre-fragments. The form of the  !
    !  operator is:                                                       !
    !                   A(z) = | 2 H(z-zN) - 1 | / A                      !
    !                                                                     !
    !  where H(x) is the Heaviside step function and zN the z-coordinate  !
    !  of the neck.                                                       !
    !---------------------------------------------------------------------!
    Subroutine asymmetry_value()

      Use MAD_PP
      Use NEWCNT

      Integer(ipr) :: i,j,LDBASE
      Complex(pr) :: asymmetry

      COMMON &
             /DIMENS/ LDBASE

      asymmetry = Cmplx(0.0_pr,0.0_pr)
      Do j=1,LDBASE
         Do i=1,LDBASE
            asymmetry = asymmetry + Cmplx(ARASYM(i,j),0.0_pr)        &
                                *(DEN_PP(j,i,0,0)+DEN_PP(j,i,1,0)  &
                                + DEN_PP(j,i,0,1)+DEN_PP(j,i,1,1))
         End Do
      End Do

      asymmetry_prefragments = Real(asymmetry,Kind=pr)

    End Subroutine asymmetry_value

    !---------------------------------------------------------------------!
    !  The routine computes the matrix elements of the operator for the   !
    !  mass-asymmetry between the two pre-fragments. The mass asymmetry   !
    !  is defined by |A_R - A_L|/A, and the corresponding operator reads  !
    !                                                                     !
    !                   A(z) = | 2 H(z-zN) - 1 | / A                      !
    !                                                                     !
    !  where H(x) is the Heaviside step function and zN the z-coordinate  !
    !  of the neck. The asymmetry operator is real, time-invariant and    !
    !  simplex-invariant. It has only non-zero matrix elements in the     !
    !  (t+,t+) block.                                                     !
    !---------------------------------------------------------------------!
    Subroutine asymmetry_operator(IN_FIX,IZ_FIX,NXHERM,NYHERM,NZHERM,NZMAXX, &
                                  ITERUN,ISIMPY,ISIGNY,IPARTY,SLOWOD)

      Use NEWCNT

      Integer(ipr), INTENT(IN) :: IN_FIX,IZ_FIX,NXHERM,NYHERM,NZHERM,NZMAXX,ITERUN,ISIMPY,ISIGNY,IPARTY
      Real(pr), INTENT(IN) :: SLOWOD

      Integer(ipr) :: NXVECT,NYVECT,NZVECT,LDBASE,IPHAPP,IPHAPM,IPHAMP,IPHAMM
      Logical :: flip,asymmetry
      Integer(ipr) :: KARTEZ,i,j,nx,ny,nz,mx,my,mz
      Real(pr) :: mixing
      Real(pr), Allocatable :: HAUXPP(:,:),D_MATR(:,:)

      COMMON                                              &
             /T_PHAS/ IPHAPP(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAPM(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAMP(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAMM(0:NDYMAX,0:NDYMAX,0:NDKART)
      COMMON                                                             &
             /BASISO/ NXVECT(1:NDBASE),NYVECT(1:NDBASE),NZVECT(1:NDBASE)
      COMMON                 &
             /DIMENS/ LDBASE

      If(.Not.Allocated(ARASYM)) Then
         Allocate(ARASYM(1:NDBASE,1:NDBASE))
         ARASYM(:,:) = 0.0_pr
      End If

      ! Computing the mass asymmetry between fragments
      Allocate(D_MATR(0:NZMAXX,0:NZMAXX)); D_MATR(:,:)=0.0_pr
      KARTEZ=3; flip=.False.; asymmetry = .True.
      Call define_asymmetry(NZHERM,NZMAXX,KARTEZ,D_MATR,Z_NECK,flip,asymmetry)

      Allocate(HAUXPP(1:NDBASE,1:NDBASE))
      HAUXPP(:,:)=0.0_pr
      Do j=1,LDBASE
         nx=NXVECT(j); ny=NYVECT(j); nz=NZVECT(j)
         Do i=1,LDBASE
            mx=NXVECT(i); my=NYVECT(i); mz=NZVECT(i)

            If(mx.Eq.nx.And.my.Eq.ny) Then
               HAUXPP(i,j)=D_MATR(mz,nz)*IPHAPP(my,ny,0)/Real(IN_FIX+IZ_FIX,Kind=pr)
            End If

         End Do
      End Do

      Deallocate(D_MATR)

      ! Linear mixing of matrix elements
      If(asymmetry_is_calculated) Then
         mixing = SLOWOD
      Else
         mixing = 0.0_pr
      End If

      Do j=1,LDBASE
         Do i=1,LDBASE
            ARASYM(i,j)=(1.0_pr-mixing)*HAUXPP(i,j)+mixing*ARASYM(i,j)
         End Do
      End Do

      Deallocate(HAUXPP)

    End Subroutine asymmetry_operator

    !---------------------------------------------------------------------!
    !  The routine computes the integral of the asymmetry operator along  !
    !  the z-direction:                                                   !
    !                                                                     !
    !      D_mn = \sum_k C_mn^k(00) \int dxi H_k(xi) e^(-xi^2) A(z)       !
    !                                                                     !
    !  where xi = b_z * z, H_k(xi) is the (normalized) Hermite polynomial !
    !  used throughout in HFODD, and C_mn^k(00) relates the product of 2  !
    !  Hermite polynomials to a sum of Hermite polynomials. The output is !
    !  the matrix D(m_z, n_z) along the z-direction only. Integrations    !
    !  are performed by Gauss-Hermite quadrature or numericall by using   !
    !  Simpson's 3/8 rule.                                                !
    !---------------------------------------------------------------------!
    Subroutine define_asymmetry(NZHERM,N_MAXX,KARTEZ,D_MATR,Z_POSI,flip,asymmetry)

      Logical :: flip,asymmetry
      Integer(ipr), INTENT(IN) :: NZHERM,N_MAXX,KARTEZ
      Real(pr), INTENT(IN) :: Z_POSI
      Real(pr), Allocatable, INTENT(INOUT) :: D_MATR(:,:)

      Real(pr) :: COEF00,COEF01,COEF11,COEF02,HOMSCA,TWOWGT,TWOPNT,HERONE

      Logical :: quadrature
      Integer(ipr) :: k,n,m,IERROR,i,iz,Nx
      Real(pr) :: SUM_DN,COORDI,operator
      Real(pr) :: DERIV1,DERIVN,XBEGIN,XFINIS,XARGUM,X_STEP,dh
      Real(pr), Allocatable :: WAVFUN(:,:),value(:)
      Real(pr), Allocatable :: XVALUE(:),FUNCTI(:),ZINTER(:),AUXSTO(:),DINTER(:)

      COMMON                                            &
             /HERMEM/ HERONE(0:ND2MAX,1:NDGAUS,1:NDKART)
      COMMON                                                       &
             /COEXYZ/ COEF00(0:ND2MAX,0:NDOSCI,0:NDOSCI,1:NDKART), &
                      COEF01(0:ND2MAX,0:NDOSCI,0:NDOSCI,1:NDKART), &
                      COEF11(0:ND2MAX,0:NDOSCI,0:NDOSCI,1:NDKART), &
                      COEF02(0:ND2MAX,0:NDOSCI,0:NDOSCI,1:NDKART)
      COMMON                         &
             /SCALNG/ HOMSCA(NDKART)
      COMMON                                     &
             /INTMEM/ TWOWGT(1:NDGAUS,1:NDKART), &
                      TWOPNT(1:NDGAUS,1:NDKART)

      Allocate(value(0:2*N_MAXX))
      Allocate(WAVFUN(1:NZHERM,0:2*N_MAXX)); WAVFUN(:,:)=0.0_pr

      Nx = 1000000

      ! The function to integrate is H_k^(0)(xi) e^(-xi^2) and Hermite polynomials are
      ! multiplied by the factor \sqrt(HOMSCA(KARTEZ)) in routine DEFINT of HFODD. We
      ! need to remove this scaling here.
      Do iz=1,NZHERM
         COORDI = TWOPNT(NZHERM-iz+1,KARTEZ)
         Do k=0,2*N_MAXX
         WAVFUN(iz,k) = HERONE(k,NZHERM-iz+1,KARTEZ)/Sqrt(HOMSCA(KARTEZ)) * Exp(-COORDI*COORDI)
      End DO
      End Do
      Allocate(FUNCTI(1:NZHERM),XVALUE(1:NZHERM),AUXSTO(1:NZHERM))
      Allocate(DINTER(1:Nx),ZINTER(1:Nx))

      ! Define boundaries for integration
      XBEGIN=TWOPNT(NZHERM,KARTEZ) + 1.0D-14
      XFINIS=TWOPNT(1,KARTEZ) - 1.0D-14

!$OMP PARALLEL DO &
!$OMP DEFAULT(NONE) &
!$OMP SCHEDULE(STATIC) &
!$OMP SHARED(N_MAXX,NZHERM,TWOPNT,KARTEZ,WAVFUN,Nx,XBEGIN,XFINIS,HOMSCA,value,Z_POSI,asymmetry,flip,prec_heaviside) &
!$OMP PRIVATE(k,iz,XVALUE,FUNCTI,dh,DERIV1,DERIVN,IERROR,AUXSTO,i,XARGUM,&
!$OMP         operator,ZINTER,DINTER,X_STEP,SUM_DN)
      Do k=0,2*N_MAXX

         ! Define function to interpolate
         Do iz=1,NZHERM
            XVALUE(iz)=TWOPNT(NZHERM-iz+1,KARTEZ)
            FUNCTI(iz)=WAVFUN(iz,k)
         End Do

         ! Approximation of derivatives at the boundaries
         dh=XVALUE(2)-XVALUE(1) ! >0
         DERIV1=+0.5_pr*(-3.0_pr*FUNCTI(1)     +4.0_pr*FUNCTI(2)       -FUNCTI(3))       /dh
         DERIVN=-0.5_pr*(-3.0_pr*FUNCTI(NZHERM)+4.0_pr*FUNCTI(NZHERM-1)-FUNCTI(NZHERM-2))/dh

         ! Define spline coefficients
         IERROR=0; Call SPLINE(XVALUE,FUNCTI,NZHERM,DERIV1,DERIVN,AUXSTO,IERROR)

         ! Define interpolated mesh and value of the function of that mesh
         If(IERROR.Eq.0) Then
            Do i=1,Nx
               XARGUM=XBEGIN+(XFINIS-XBEGIN)*Real(i-1,Kind=pr)/Real(Nx-1,Kind=pr)
                             operator = Heaviside(XARGUM-Z_POSI,prec_heaviside)
               If(flip)      operator = 1.0_pr-Heaviside(XARGUM-Z_POSI,prec_heaviside)
               If(asymmetry) operator = 2.0_pr*Heaviside(XARGUM-Z_POSI,prec_heaviside) - 1.0_pr
               ZINTER(i)=XARGUM
               DINTER(i)=operator*SPLINT(XVALUE,FUNCTI,AUXSTO,NZHERM,XARGUM)
            End Do
         Else
            Write(6,'("In define_asymmetry() - Error in interpolating the function!")')
         End If

         ! Define step for integration
         X_STEP=(ZINTER(Nx)-ZINTER(1))/Real(Nx-1,Kind=pr)

         ! Integrate
         Call SIMP38(DINTER,Nx,X_STEP,SUM_DN)
         value(k)=SUM_DN

      End Do
!$OMP END PARALLEL DO

      Deallocate(FUNCTI,XVALUE,AUXSTO,DINTER,ZINTER)

      ! Filling the matrix. Attention: Coefficients COEF00 are scaled by
      ! SQRT(HOMSCA(KARTEZ)) in DEVHER, we need to remove this scaling here
      Do n=0,N_MAXX
         Do m=0,N_MAXX
            SUM_DN=0.0_pr
            Do k=0,n+m
               SUM_DN=SUM_DN+COEF00(k,m,n,KARTEZ)/Sqrt(HOMSCA(KARTEZ))*value(k)
            End Do
            D_MATR(m,n)=SUM_DN
         End Do
      End Do

      Deallocate(WAVFUN)
      Deallocate(value)

    End Subroutine define_asymmetry

    !---------------------------------------------------------------------!
    !  This function computes the Heaviside step ffunction by using the
    !  approximation
    !      H(x) ~ 1/2 lim_{t->0} [ 1 + tanh(x)]
    !---------------------------------------------------------------------!
    Real(pr) Function Heaviside(x,t)

      Real(pr), INTENT(IN) :: x,t

      Heaviside = 0.5_pr*(1.0_pr + tanh(x/t))

    End Function Heaviside

    !---------------------------------------------------------------------!
    !  This function computes the probability integral at point z > 0     !
    !      PHI(Z) = FRAC{2}{SQRT{PI}} INT_{0}^{Z} EXP(-T^2)DT             !
    !  by using the Simpson 3/8 integration rule                          !
    !---------------------------------------------------------------------!
    Real(pr) Function PROINT(ZPOINT)

      Real(pr), INTENT(IN) :: ZPOINT

      Integer(ipr) :: i, Nx
      Real(pr) :: dx,x,PIARGU,VALRES
      Real(pr), Allocatable :: FUNCTI(:)

      Nx = 10000

      Allocate(FUNCTI(1:Nx))

      dx = Abs(ZPOINT/Real(Nx-1,Kind=pr))
      Do i=1,Nx
         x=dx*Real(i-1,Kind=pr); FUNCTI(i)=Exp(-x**2)
      End Do

      VALRES=0.0_pr; Call SIMP38(FUNCTI,Nx,dx,VALRES)

      Deallocate(FUNCTI)

      PIARGU=4.0_pr*Atan(1.0_pr)
      PROINT=2.0_pr*VALRES/Sqrt(PIARGU)

    End Function PROINT

    !---------------------------------------------------------------------!
    !  This routine performs the composite simpson's rule integration     !
    !  of a function F defined by a table of N equispaced values.         !
    !                                                                     !
    !                       See: KOONIN, Computational physics, P.9       !
    !                                                                     !
    !   The parameters are:                                               !
    !    F = array of values of the function F(X)                         !
    !    N = number of points X_K                                         !
    !    H = the uniform spacing between x values: H = X_K+1 - X_K        !
    !    RESULT = estimate of the integral that is returned to caller.    !
    !---------------------------------------------------------------------!
    Subroutine SIMP38(FUNCTI,NINTEG,X_STEP,RESULT)

      Integer(ipr), INTENT(IN) :: NINTEG
      Real(pr), INTENT(INOUT) :: RESULT
      Real(pr), INTENT(IN) :: X_STEP
      Real(pr), Dimension(1:NINTEG), INTENT(IN) :: FUNCTI

      Integer(ipr) :: NBEGIN,NPANEL,N_HALF,NENDIN,i
      Real(pr) :: VALSUM

      ! Check to see if number of panels is even. Number of panels is n - 1.
      NBEGIN=1; NPANEL=NINTEG-NBEGIN; N_HALF=NPANEL/2; RESULT=0.0_pr

      ! Number of panels is odd.  Use Simpson's 3/8 rule on first three
      ! panels, 1/3 rule on rest of them.
      If((NPANEL-2*N_HALF).Ne.0) Then

         RESULT = 3.0_pr*X_STEP*(FUNCTI(NBEGIN)               &
                + 3.0_pr*(FUNCTI(NBEGIN+1)+FUNCTI(NBEGIN+2))  &
                        + FUNCTI(NBEGIN+3))/8.0_pr

         If((NINTEG-NBEGIN).Eq.3) Return

         NBEGIN=NBEGIN+3

      End If

      ! Apply 1/3 rule - add in first, second, last values
      RESULT = RESULT + X_STEP*(FUNCTI(NBEGIN)        &
                       + 4.0_pr*FUNCTI(NBEGIN+1)      &
                              + FUNCTI(NINTEG))/3.0_pr
      NBEGIN = NBEGIN+2

      If(NBEGIN.Eq.NINTEG) Then
         Return
      Else

         VALSUM=0.0_pr; NENDIN=NINTEG - 1
         Do i=NBEGIN,NENDIN,2
            VALSUM=VALSUM+FUNCTI(i)+2.0_pr*FUNCTI(i+1)
         End Do

         RESULT=RESULT+2.0_pr*X_STEP*VALSUM/3.0_pr

         Return

      End If

    End Subroutine SIMP38

    !---------------------------------------------------------------------!
    !     This subroutine computes the value of the associated Legendre   !
    !     polynomial Plm(x) at value x = ZVALUE and with order and        !
    !     degree l = LAMACT and m = MIUACT                                !
    !---------------------------------------------------------------------!
    Real(pr) Function DEFLEG(LAMACT,MIUACT,ZVALUE)

      Integer(ipr), INTENT(IN) :: LAMACT,MIUACT
      Real(pr), INTENT(IN) :: ZVALUE

      Integer(ipr) :: L,M
      Real(pr) :: ZVALPP,ZVAL_P,ZPOLYN,TMPVAL,epsilon,one

      ! Argument z of polynomial P_{l,m}(z) must be lower than 1

      epsilon=1.D-14; one=1.0_pr
      If(Abs(Abs(ZVALUE)-one).Le.epsilon) Then
         Write(6,'("ZVALUE=",f20.16," one=",f20.16," epsilon=",f20.16)') ZVALUE,one,epsilon
         Stop 'Error in DEFLEG - ARGUMENT |z| = 1'
      End If

      ! Initialization
      ZVALPP = one; ZVAL_P = ZVALUE

      ! P{0,0}(z) = 1
      If (LAMACT.Eq.0) Then
          DEFLEG=one; Return
      End If

      ! P{1,0}(z) = 1
      If (LAMACT.Eq.1.And.MIUACT.Eq.0) Then
          DEFLEG=ZVALUE; Return
      End If

      ! P{1,1}(z) = 1
      If (LAMACT.Eq.1.And.MIUACT.Eq.1) Then
          DEFLEG=Sqrt(one-ZVALUE**2); Return
      End If

      ! Obtaining P_{l,0}(z) by recurrence
      If (LAMACT.Gt.1) Then

          Do L = 1,LAMACT-1
             ZPOLYN = (Real(2*L+1,Kind=pr)*ZVALUE*ZVAL_P - Real(L,Kind=pr)*ZVALPP) /Real(L+1,Kind=pr)
             TMPVAL=ZVAL_P; ZVAL_P=ZPOLYN; ZVALPP=TMPVAL
          End Do

          If (MIUACT.Eq.0) Then
              DEFLEG=ZVAL_P
          Else

              ! Obtaining P_{l,1} from P_{l,0}(z) and P_{l-1,0}(z)
              !   - ZVAL_P contains P_{l,0}(z)
              !   - ZVALPP contains P_{l-1,0}(z)
              ZPOLYN =( Real(LAMACT,Kind=pr)*ZVALUE*ZVAL_P                 &
                      - Real(LAMACT,Kind=pr)*ZVALPP)/Sqrt(one - ZVALUE**2)
              TMPVAL=ZVAL_P; ZVAL_P=ZPOLYN; ZVALPP=TMPVAL

              ! Obtaining P_{l,m+1} from P_{l,m}(z) and P_{l,m-1}(z)
              !   - ZVAL_P contains P_{l,1}(z)
              !   - ZVALPP contains P_{l,0}(z)
              Do M = 1,MIUACT-1
                 ZPOLYN =-2.0_pr*Real(M,Kind=pr)*ZVALUE*ZVAL_P/Sqrt(one - ZVALUE**2) &
                               - Real((LAMACT+M)*(LAMACT-M+1),Kind=pr)*ZVALPP
                 TMPVAL=ZVAL_P; ZVAL_P=ZPOLYN; ZVALPP=TMPVAL
              End Do

              DEFLEG=(-1)**MIUACT * ZVAL_P

          End If

      End If

    End Function DEFLEG

    !---------------------------------------------------------------------!
    !  Using Brent's method, find the root of a function 'FONCTI' known   !
    !  to lie between 'XLOWER' and 'XUPPER'. The root, returned as        !
    !  'ZBRENT', will be refined until its accuracy is 'TOLERA'.          !
    !  Parameters: - maximum allowed number of iterations, and            !
    !              - machine floating-point precision.                    !
    !---------------------------------------------------------------------!
    Real(pr) Function ZBRENT(FONCTI,XLOWER,XUPPER,TOLERA,IFOUND)

      Integer(ipr), INTENT(INOUT) :: IFOUND
      Real(pr), INTENT(IN) :: XLOWER,XUPPER,TOLERA

      Integer(ipr), Parameter :: ITRMAX=100
      Integer(ipr) :: ITERAT
      Real(pr), Parameter :: EPSCPU=3.D-8
      Real(pr) :: A,B,C,D,E,FA,FB,FC,P,Q,R,S,EPSILO,XM

      Interface
         Real(pr) Function FONCTI(X)
           Use hfodd_fission_precision
           Implicit None
           Real(pr), INTENT(IN) :: X
         End Function FONCTI
      End Interface

      IFOUND=1; A=XLOWER; B=XUPPER; FA=FONCTI(A); FB=FONCTI(B)

      If((FA.Gt.0.0_pr.And.FB.Gt.0.0_pr).Or.(FA.Lt.0.0_pr.And.FB.Lt.0.0_pr)) Then
         IFOUND=0
         If(debug.Ge.2) Then
            Write(6,'("A = ",F20.14," FA = ",f20.14)') A,FA
            Write(6,'("B = ",F20.14," FB = ",f20.14)') B,FB
            Write(6,'("ROOT MUST BE BRACKETED FOR ZBRENT")')
         End if
         Return
      End If

      C=B; FC=FB

      Do ITERAT=1,ITRMAX

         If(debug.Ge.3) Then
            Write(6,'("ITERAT=",i6)') ITERAT
            Write(6,'("A = ",f20.14," B = ",f20.14)') A,B
            Write(6,'("FA= ",f20.14," FB= ",f20.14)') FA,FB
         End If
         If((FB.Gt.0.0_pr.And.FC.Gt.0.0_pr).Or.(FB.Lt.0.0_pr.And.FC.Lt.0.0_pr)) Then
           C=A; FC=FA; D=B-A; E=D
         End If

         If(Abs(FC).LT.Abs(FB)) Then
           A=B; B=C; C=A; FA=FB; FB=FC; FC=FA
         End If

         EPSILO=2.0_pr*EPSCPU*Abs(B) + 0.5_pr*TOLERA; XM=0.5_pr*(C - B)

         If(Abs(XM).Le.EPSILO .Or. Abs(FB).Le.EPSCPU) Then
            ZBRENT=B
            If(debug.Ge.2) Then
               Write(6,'("Returning: ITERAT=",i6)') ITERAT
               Write(6,'("             A = ",f20.14," B = ",f20.14)') A,B
               Write(6,'("             FA= ",f20.14," FB= ",f20.14)') FA,FB
            End If
            Return
         End If

         If(Abs(E).Ge.EPSILO .And. Abs(FA).Ge.Abs(FB)) Then

             S = FB/FA

             If(Abs(A-C).Le.EPSCPU) Then
                P=2.0_pr*XM*S
                Q=1.0_pr - S
             Else
                ! Attempt inverse quadratic interpolation
                Q=FA/FC; R=FB/FC
                P=S*(2.0_pr*XM*Q*(Q - R) - (B - A)*(R - 1.0_pr))
                Q=(Q - 1.0_pr)*(R - 1.0_pr)*(S - 1.0_pr)
             End If

             ! Check whether in bounds
             If (P.Gt.0.0_pr) Q = -Q

             P=Abs(P)

             ! Test quality of interpolation. If too bad, switch to bisection method
             If(2.0_pr*P .Lt. Min(3.0_pr*XM*Q-Abs(EPSILO*Q),Abs(E*Q))) Then
                E=D; D=P/Q
             Else
                D=XM; E=D
             End If

         Else ! Bounds decreasing too slowly, use bisection

            D=XM; E=D

         End If

         ! Move last best guess to A
         A=B; FA=FB

         ! Evaluate new trial root
         If(Abs(D) .Gt. EPSILO) Then
            B=B+D
         Else
            B=B+Sign(EPSILO,XM)
         End If

         FB=FONCTI(B)

      End Do

      IFOUND=2
      ZBRENT=B

    End Function ZBRENT

    !---------------------------------------------------------------------!
    !  Cubic spline interpolation. Given the arrays XENTRY and YENTRY of  !
    !  length NPOINT  (of x- and y- values of a function;  XENTRY must be !
    !  in ascending order).  For the array of second derivatives,         !
    !  Y2A, the routine gives the value of the function in point          !
    !  XARGUM. The vector Y2A can be obtained by using the                !
    !  routine SPLINE given below (SPLINE is called only once)            !
    !---------------------------------------------------------------------!
    Real(pr) Function SPLINT(XENTRY,YENTRY,Y2A,NPOINT,XARGUM)

      Integer(ipr), INTENT(IN) :: NPOINT
      Real(pr), INTENT(IN) :: XARGUM
      Real(pr), Dimension(1:NPOINT), INTENT(IN) :: XENTRY,YENTRY,Y2A

      Integer(ipr) :: KLOWER, KHIGHR, K
      Real(pr) :: HLNGTH,A,B

      KLOWER=1; KHIGHR=NPOINT
      Do While((KHIGHR-KLOWER).Gt.1)
         K=(KHIGHR+KLOWER)/2
         If(XENTRY(K).Gt.XARGUM) Then
            KHIGHR=K
         Else
            KLOWER=K
         End If
      End Do

      HLNGTH=XENTRY(KHIGHR)-XENTRY(KLOWER)

      If(HLNGTH.EQ.0) Stop 'SPLINT01'

      A=(XENTRY(KHIGHR)-XARGUM)/HLNGTH
      B=(XARGUM-XENTRY(KLOWER))/HLNGTH

      SPLINT=A*YENTRY(KLOWER)+B*YENTRY(KHIGHR) &
                        +((A**3-A)*Y2A(KLOWER) &
                        + (B**3-B)*Y2A(KHIGHR))*(HLNGTH**2)/6.0_pr

    End Function SPLINT

    !---------------------------------------------------------------------!
    !  Cubic spline interpolation. Given the arrays XENTRY and YENTRY of  !
    !  length NPOINT  (of x- and y- values of a function;  XENTRY must be !
    !  in ascending order).  For the array of second derivatives,         !
    !  Y2A, the routine gives the value of the function in point          !
    !  XARGUM. The vector Y2A can be obtained by using the                !
    !  routine SPLINE given below (SPLINE is called only once)            !
    !---------------------------------------------------------------------!
    Real(pr) Function SPLIND(XENTRY,YENTRY,Y2A,NPOINT,XARGUM)

      Integer(ipr), INTENT(IN) :: NPOINT
      Real(pr), INTENT(IN) :: XARGUM
      Real(pr), Dimension(1:NPOINT), INTENT(IN) :: XENTRY,YENTRY,Y2A

      Integer(ipr) :: KLOWER, KHIGHR, K
      Real(pr) :: HLNGTH,A,B

      KLOWER=1; KHIGHR=NPOINT
      Do While((KHIGHR-KLOWER).Gt.1)
         K=(KHIGHR+KLOWER)/2
         If(XENTRY(K).Gt.XARGUM) Then
            KHIGHR=K
         Else
            KLOWER=K
         End If
      End Do

      HLNGTH=XENTRY(KHIGHR)-XENTRY(KLOWER)

      If(HLNGTH.EQ.0) Stop 'SPLIND01'

      A=(XENTRY(KHIGHR)-XARGUM)/HLNGTH
      B=(XARGUM-XENTRY(KLOWER))/HLNGTH

      SPLIND=(YENTRY(KHIGHR)-YENTRY(KLOWER))/HLNGTH        &
                        -((2.0_pr*A**2-1.0_pr)*Y2A(KLOWER) &
                        - (3.0_pr*B**2-1.0_pr)*Y2A(KHIGHR))*HLNGTH/6.0_pr

    End Function SPLIND

    !---------------------------------------------------------------------!
    !  Given X and Y (x's and y's of a function) the routine gives        !
    !  the array Y2 of length N containing the second derivatives of      !
    !  the function. YP1 and YPN are the values of the first deriva-      !
    !  tive at first and last points. The value .GT.1.D+30 means that     !
    !  natural spline is to be used at this point.                        !
    !---------------------------------------------------------------------!
    Subroutine SPLINE(X,Y,N,YP1,YPN,Y2,IERROR)

      Integer(ipr), INTENT(IN) :: N
      Integer(ipr), INTENT(INOUT) :: IERROR
      Real(pr), INTENT(IN) :: YP1,YPN
      Real(pr), Dimension(1:N), INTENT(IN) :: X,Y
      Real(pr), Dimension(1:N), INTENT(OUT) :: Y2

      Integer(ipr), Parameter :: NDSPLN=100
      Integer(ipr) :: I,K
      Real(pr) :: SIG,P,QN,UN
      Real(pr), Dimension(1:NDSPLN) :: U

                                IERROR=0
      If(N.Lt.4.Or.N.Gt.NDSPLN) IERROR=1

      If(YP1.Gt.0.99D+30) Then
         Y2(1)=0.0_pr
         U(1)=0.0_pr
      Else
         Y2(1)=-0.5_pr
         U(1)=(3.0_pr/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      End If

      Do I=2,N-1
         SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
         P=SIG*Y2(I-1)+2.0_pr
         Y2(I)=(SIG-1.0_pr)/P
         U(I)=(6.0_pr*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      End Do

      If(YPN.Gt.0.99D+30) Then
         QN=0.0_pr
         UN=0.0_pr
      Else
         QN=0.5_pr
         UN=(3.0_pr/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      End If

      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1)

      Do K=N-1,1,-1
         Y2(K)=Y2(K)*Y2(K+1)+U(K)
      End Do

    End Subroutine SPLINE

 End Module hfodd_fission_fragments

 ! ==================================================================== !
 !                                                                      !
 !               QUANTUM LOCALIZATION PACKAGE                           !
 !                                                                      !
 ! ==================================================================== !

 !----------------------------------------------------------------------!
 !                                                                      !
 !  This module contains all routines needed to identify each fragment  !
 !  based on the occupation of their quasiparticles, as well as to      !
 !  rotate said quasiparticles to minimize the tails between fragments. !
 !                                                                      !
 !  Inputs (direct or indirect):                                        !
 !    - FERMFN ......: Fermi-Dirac occupations of quasiparticles (=0 at !
 !                     zero temperature)                                !
 !    - EQPISO ......: Quasiparticle energies                           !
 !    - ASVQUA,BSVQUA: U and V (U=A*, B=V*) of the Bogoliubov transfor- !
 !                     mation                                           !
 !    - Z_NECK ......: position of the neck                             !
 !                                                                      !
 !      With the U and the V, Z_NECK, and FERMFN, the subroutine        !
 !      wave_localization() can define occupations in the left and      !
 !      right fragments, as well as the localization indicator. This    !
 !      is done at T>=0 always before rotation of quasiparticles. The   !
 !      localization method also requires quasiparticle energies, the   !
 !      equivalent single-particle energies and the Fermi level. This   !
 !      allows defining a set of candidate pairs for the rotation.      !
 !                                                                      !
 !  Outputs:                                                            !
 !    - FERMFN ......: At T>0, when quasiparticles are rotated, the     !
 !                     generalized density matrix is not diagonal any   !
 !                     longer. FERMFN contains the diagonal terms.      !
 !    - EQPISO ......: When quasiparticle are rotated, the HFB matrix   !
 !                     not diagonal anymore. EQPISO contains the        !
 !                     diagonal term (by convention).                   !
 !    - F_FLAG ......: the value of element number 'mu' indicates if    !
 !                     the corresponding quasiparticle is occupied      !
 !                     (=1.0) or not (=0.0). It is used in DENSHF to    !
 !                     compute densities in r-space for each fragment   !
 !                     based on qp occupations (in filter_density).     !
 !    - flag_left ...: Same as F_FLAG                                   !
 !    - flag_right ..: Convenience array. It is defined as 1-f_left     !
 !                                                                      !
 !      F_FLAG is used in DENSHF to filter out given qp from the calcu- !
 !      tion of the densities depending on whether they belong to a     !
 !      given fragment. flag_left and flag_right are used in the        !
 !      hfodd_fission_interaction module (subroutine filter_density()   !
 !      to set the value of F_FLAG.                                     !
 !                                                                      !
 !      Localization package for DENSHF                                 !
 !       - number_of_pairs ....: Number of rotated qp pairs             !
 !       - table_singles ......: Array keeping track of which qp are    !
 !                               rotated                                !
 !       - table_pair_mu ......: Index (tail size ordering) of q.p. mu  !
 !       - table_pair_nu       : and nu                                 !
 !       - table_pair_order_qp : Array giving the mapping between tail  !
 !                               size ordering and regular ordering     !
 !       - table_pair_fTmu ....: Diagonal terms of the Fermi-Dirac      !
 !         table_pair_fTnu       occupation matrix for the qp pair      !
 !                               (mu,nu)                                !
 !       - table_pair_fTmunu ..: Off-diagonal germ of the FD occupation !
 !                               matrix (which is symmetric)            !
 !       - table_pair_gTmu ....: By definition, g = 1-f                 !
 !         table_pair_gTnu                                              !
 !         table_pair_gTmunu                                            !
 !       - table_pair_Emu .....: Diagonal term of the HFB matrix after  !
 !         table_pair_Enu        rotation for the qp pair (mu,nu)       !
 !       - table_pair_Emunu ...: Off-diagonal term of the HFB matrix    !
 !                                                                      !
 !         All these arrays are used in DENSHF at T>0 when localization !
 !         is applied. At T=0, or at T>0 without localization, these    !
 !         arrays are not needed since the generalized density is still !
 !         diagonal.                                                    !
 !                                                                      !
 !----------------------------------------------------------------------!

 Module hfodd_fission_rotated_qp

    Use hfodd_sizes
    Use hfodd_fission_precision
    Use hfodd_fission_fragments
    Use tilings_utilities
#if(USE_MPI==1)
    Use hfodd_mpiio
#endif

    Implicit None

    ! PRIVATE VARIABLES

    ! Quantities pertaining to the localization package
    Integer(ipr), Parameter, PRIVATE :: number_of_pairs_max = 100
    Integer(ipr), PRIVATE, SAVE :: iteration_max, number_angles, Nstate
    Real(pr), PRIVATE, SAVE :: energy_window, loc_max, occ_min

    ! Array giving the correspondence between tail and regular indexing: regular = table(tail)
    Integer(ipr), Allocatable, PRIVATE, SAVE :: table(:)

    ! Arrays giving quantities with 'regular' indexing
    Integer(ipr), Allocatable, PRIVATE, SAVE :: table_partner(:)
    Real(pr), Allocatable, PRIVATE, SAVE :: localisation(:,:), v_left(:,:), vright(:,:)
    Real(pr), Allocatable, PRIVATE, SAVE :: diagonal_off_f(:,:), diagonal_off_g(:,:), diagonal_off_E(:,:)
    Complex(pr), Allocatable, PRIVATE, SAVE :: densityR_lower(:,:), densityR_upper(:,:)

    ! Arrays giving quantities with 'pair' indexing
    Integer(ipr), Allocatable, PRIVATE, SAVE :: final_mu(:), final_nu(:)

    ! Arrays giving quantities with 'tail' indexing
    Real(pr), Allocatable, PRIVATE, SAVE :: table_theta(:),table_Esp(:), table_occ(:)
    Real(pr), Allocatable, PRIVATE, SAVE :: table_loc(:), table_vlef(:), table_vrig(:), table_flag(:)
    Real(pr), Allocatable, PRIVATE, SAVE :: table_Eqp(:),off_diagonal_E(:)
    Real(pr), Allocatable, PRIVATE, SAVE :: table_fT(:), off_diagonal_f(:)
    Real(pr), Allocatable, PRIVATE, SAVE :: table_gT(:), off_diagonal_g(:)

    ! Internal arrays needed by the localization method
    Integer(ipr), Allocatable, PRIVATE, SAVE :: table_iter(:,:), table_ipair(:,:)
    Integer(ipr), Allocatable, PRIVATE, SAVE :: temp_order_qp(:,:),temp_pair_mu(:,:), temp_pair_nu(:,:)
    Real(pr), Allocatable, PRIVATE, SAVE :: temp_pair_mu_vlef(:,:), temp_pair_mu_vrig(:,:), temp_pair_mu_loc(:,:), temp_pair_mu_flag(:,:)
    Real(pr), Allocatable, PRIVATE, SAVE :: temp_pair_nu_vlef(:,:), temp_pair_nu_vrig(:,:), temp_pair_nu_loc(:,:), temp_pair_nu_flag(:,:)
    Real(pr), Allocatable, PRIVATE, SAVE :: temp_pair_fTmu(:,:), temp_pair_fTnu(:,:), temp_pair_fTmunu(:,:)
    Real(pr), Allocatable, PRIVATE, SAVE :: temp_pair_gTmu(:,:), temp_pair_gTnu(:,:), temp_pair_gTmunu(:,:)
    Real(pr), Allocatable, PRIVATE, SAVE :: temp_pair_Emu(:,:), temp_pair_Enu(:,:), temp_pair_Emunu(:,:)

    Type :: liste_int
        Integer(ipr) :: numero
        Integer(ipr) :: value_mu
        Integer(ipr) :: value_nu
        Real(pr):: score
        TYPE(liste_int), pointer :: next
    End type liste_int

    ! PUBLIC VARIABLES

    ! Arrays needed by DENSHF to compute densities at T>0 after localization
    Integer(ipr), Allocatable, PUBLIC, SAVE :: number_of_pairs(:), table_singles(:,:)
    Integer(ipr), Allocatable, PUBLIC, SAVE :: table_pair_mu(:,:), table_pair_nu(:,:), table_pair_order_qp(:,:)
    Real(pr), Allocatable, PUBLIC, SAVE :: table_pair_mu_vlef(:,:), table_pair_mu_vrig(:,:), table_pair_mu_loc(:,:), table_pair_mu_flag(:,:)
    Real(pr), Allocatable, PUBLIC, SAVE :: table_pair_nu_vlef(:,:), table_pair_nu_vrig(:,:), table_pair_nu_loc(:,:), table_pair_nu_flag(:,:)
    Real(pr), Allocatable, PUBLIC, SAVE :: table_pair_fTmu(:,:), table_pair_fTnu(:,:), table_pair_fTmunu(:,:)
    Real(pr), Allocatable, PUBLIC, SAVE :: table_pair_gTmu(:,:), table_pair_gTnu(:,:), table_pair_gTmunu(:,:)
    Real(pr), Allocatable, PUBLIC, SAVE :: table_pair_Emu(:,:), table_pair_Enu(:,:), table_pair_Emunu(:,:)

    ! Arrays needed by filter_density() to define what is in the left/right fragment
    Real(pr), Allocatable, PUBLIC, SAVE :: flag_left(:,:,:), flag_right(:,:,:)

 Contains

    !---------------------------------------------------------------------!
    ! This set of subroutines allocates (and initializes) and deallocates !
    ! arrays used in wave_localization() to characterize properties of    !
    ! fission fragments.                                                  !
    !---------------------------------------------------------------------!
    Subroutine allocate_loc(MIN_QP)
      Integer, INTENT(IN) :: MIN_QP

      Allocate(vright(1:NDSTAT,0:NDISOS),v_left(1:NDSTAT,0:NDISOS))
      Allocate(localisation(1:NDSTAT,0:NDISOS))
      Allocate(flag_left(1:NDSTAT,0:NDREVE,0:NDISOS),flag_right(1:NDSTAT,0:NDREVE,0:NDISOS))
      vright(:,:)=0.0_pr
      v_left(:,:)=0.0_pr
      localisation(:,:)=1.0_pr
      flag_right(:,:,:)=1.0_pr
      flag_left(:,:,:)=0.0_pr
      Allocate(diagonal_off_E(1:NDSTAT,0:NDISOS))
      diagonal_off_E(:,:)= 0.0_pr

      If(MIN_QP.EQ.1) THEN
         Allocate(diagonal_off_f(1:NDSTAT,0:NDISOS),diagonal_off_g(1:NDSTAT,0:NDISOS))
         diagonal_off_f(:,:)= 0.0_pr
         diagonal_off_g(:,:)= 0.0_pr
         If(debug.Ge.2) Then
            Allocate(densityR_lower(1:NDSTAT,1:NDSTAT),densityR_upper(1:NDSTAT,1:NDSTAT))
            densityR_lower(:,:)= Cmplx(0.0_pr,0.0_pr)
            densityR_upper(:,:)= Cmplx(0.0_pr,0.0_pr)
         End If
      End If

    End Subroutine allocate_loc

    Subroutine deallocate_loc(MIN_QP)
      Integer, INTENT(IN) :: MIN_QP
      Deallocate(vright,v_left)
      Deallocate(localisation)
      Deallocate(flag_left,flag_right)
      Deallocate(diagonal_off_E)
      If(MIN_QP.EQ.1) Then
         Deallocate(diagonal_off_f,diagonal_off_g)
         If(debug.Ge.2) Deallocate(densityR_lower,densityR_upper)
      End If
    End Subroutine deallocate_loc

    !---------------------------------------------------------------------!
    ! This set of subroutines allocates (and initializes) and deallocates !
    ! arrays used in DENSHF() to compute densities after rotation of qp   !
    ! at T>0.                                                             !
    !---------------------------------------------------------------------!
    Subroutine allocate_qp(n)
      Integer, INTENT(IN) :: n

      Allocate(table_pair_order_qp(1:n,0:NDISOS))
      Allocate(table_iter(1:n,0:NDISOS))
      Allocate(table_ipair(1:n,0:NDISOS))

      Allocate(table_pair_mu_vlef(1:number_of_pairs_max,0:NDISOS))
      Allocate(table_pair_mu_vrig(1:number_of_pairs_max,0:NDISOS))
      Allocate(table_pair_mu_loc(1:number_of_pairs_max,0:NDISOS))
      Allocate(table_pair_mu_flag(1:number_of_pairs_max,0:NDISOS))
      Allocate(table_pair_nu_vlef(1:number_of_pairs_max,0:NDISOS))
      Allocate(table_pair_nu_vrig(1:number_of_pairs_max,0:NDISOS))
      Allocate(table_pair_nu_loc(1:number_of_pairs_max,0:NDISOS))
      Allocate(table_pair_nu_flag(1:number_of_pairs_max,0:NDISOS))

      Allocate(table_pair_mu(1:number_of_pairs_max,0:NDISOS))
      Allocate(table_pair_nu(1:number_of_pairs_max,0:NDISOS))
      Allocate(table_pair_fTmu(1:number_of_pairs_max,0:NDISOS))
      Allocate(table_pair_fTnu(1:number_of_pairs_max,0:NDISOS))
      Allocate(table_pair_fTmunu(1:number_of_pairs_max,0:NDISOS))
      Allocate(table_pair_gTmu(1:number_of_pairs_max,0:NDISOS))
      Allocate(table_pair_gTnu(1:number_of_pairs_max,0:NDISOS))
      Allocate(table_pair_gTmunu(1:number_of_pairs_max,0:NDISOS))
      Allocate(table_pair_Emu(1:number_of_pairs_max,0:NDISOS))
      Allocate(table_pair_Enu(1:number_of_pairs_max,0:NDISOS))
      Allocate(table_pair_Emunu(1:number_of_pairs_max,0:NDISOS))

    End Subroutine allocate_qp

    Subroutine deallocate_qp()
      Deallocate(table_pair_order_qp,table_iter,table_ipair)
      Deallocate(table_pair_mu,table_pair_nu)
      Deallocate(table_pair_mu_vlef,table_pair_mu_vrig,table_pair_mu_loc,table_pair_mu_flag)
      Deallocate(table_pair_nu_vlef,table_pair_nu_vrig,table_pair_nu_loc,table_pair_nu_flag)
      Deallocate(table_pair_fTmu,table_pair_fTnu,table_pair_fTmunu)
      Deallocate(table_pair_gTmu,table_pair_gTnu,table_pair_gTmunu)
      Deallocate(table_pair_Emu,table_pair_Enu,table_pair_Emunu)
    End Subroutine deallocate_qp

    !---------------------------------------------------------------------!
    !  This routine allocates local (to this module) arrays that are used !
    !  during the localization method                                     !
    !---------------------------------------------------------------------!
    Subroutine allocate_arrays(ICHARG)

      Integer(ipr), INTENT(IN) :: ICHARG

      Integer(ipr) :: NUMBQP,IREVER

      COMMON                                   &
             /QPNUMS/ NUMBQP(0:NDREVE,0:NDISOS)

      IREVER=0
      Nstate=NUMBQP(IREVER,ICHARG)

      Allocate(table(1:Nstate))

      Allocate(table_Eqp(1:Nstate),table_Esp(1:Nstate))
      Allocate(table_occ(1:Nstate),table_loc(1:Nstate))
      Allocate(table_vlef(1:Nstate),table_vrig(1:Nstate))
      Allocate(table_fT(1:Nstate),table_gT(1:Nstate),table_flag(1:Nstate))
      ! The rotation of qp makes the HFB matrix (and generalized density
      ! matrix at T>0) no diagonal any more. We store the off-diagonal terms
      ! here
      Allocate(table_partner(1:Nstate))
      table_partner(:) = 0
      Allocate(off_diagonal_E(1:Nstate))
      off_diagonal_E(:) = 0.0_pr
      Allocate(off_diagonal_f(1:Nstate),off_diagonal_g(1:Nstate))
      off_diagonal_f(:) = 0.0_pr
      off_diagonal_g(:) = 0.0_pr

      Allocate(table_theta(1:Nstate))
      table_theta(:) = 0.0_pr

      ! Arrays needed at T>0 to keep track of pairs of rotated qp at each
      ! iteration (for a given isospin)
      Allocate(temp_order_qp(1:Nstate,1:iteration_max))
      Allocate(temp_pair_mu(1:number_of_pairs_max,1:iteration_max))
      Allocate(temp_pair_nu(1:number_of_pairs_max,1:iteration_max))
      Allocate(temp_pair_mu_vlef(1:number_of_pairs_max,1:iteration_max))
      Allocate(temp_pair_mu_vrig(1:number_of_pairs_max,1:iteration_max))
      Allocate(temp_pair_mu_loc(1:number_of_pairs_max,1:iteration_max))
      Allocate(temp_pair_mu_flag(1:number_of_pairs_max,1:iteration_max))
      Allocate(temp_pair_nu_vlef(1:number_of_pairs_max,1:iteration_max))
      Allocate(temp_pair_nu_vrig(1:number_of_pairs_max,1:iteration_max))
      Allocate(temp_pair_nu_loc(1:number_of_pairs_max,1:iteration_max))
      Allocate(temp_pair_nu_flag(1:number_of_pairs_max,1:iteration_max))
      Allocate(temp_pair_fTmu(1:number_of_pairs_max,1:iteration_max))
      Allocate(temp_pair_fTnu(1:number_of_pairs_max,1:iteration_max))
      Allocate(temp_pair_fTmunu(1:number_of_pairs_max,1:iteration_max))
      Allocate(temp_pair_gTmu(1:number_of_pairs_max,1:iteration_max))
      Allocate(temp_pair_gTnu(1:number_of_pairs_max,1:iteration_max))
      Allocate(temp_pair_gTmunu(1:number_of_pairs_max,1:iteration_max))
      Allocate(temp_pair_Emu(1:number_of_pairs_max,1:iteration_max))
      Allocate(temp_pair_Enu(1:number_of_pairs_max,1:iteration_max))
      Allocate(temp_pair_Emunu(1:number_of_pairs_max,1:iteration_max))

    End Subroutine allocate_arrays

    Subroutine deallocate_arrays()
      Deallocate(table)
      Deallocate(table_Eqp,table_Esp)
      Deallocate(table_occ,table_loc)
      Deallocate(table_vlef,table_vrig)
      Deallocate(table_fT,table_gT,table_flag)
      Deallocate(table_partner)
      Deallocate(off_diagonal_E)
      Deallocate(off_diagonal_f,off_diagonal_g)
      Deallocate(table_theta)
      Deallocate(temp_order_qp)
      Deallocate(temp_pair_mu,temp_pair_nu)
      Deallocate(temp_pair_mu_vlef,temp_pair_mu_vrig,temp_pair_mu_loc,temp_pair_mu_flag)
      Deallocate(temp_pair_nu_vlef,temp_pair_nu_vrig,temp_pair_nu_loc,temp_pair_nu_flag)
      Deallocate(temp_pair_fTmu,temp_pair_fTnu,temp_pair_fTmunu)
      Deallocate(temp_pair_gTmu,temp_pair_gTnu,temp_pair_gTmunu)
      Deallocate(temp_pair_Emu,temp_pair_Enu,temp_pair_Emunu)
    End Subroutine deallocate_arrays

    !---------------------------------------------------------------------!
    ! This subroutine passes input variables read in HFODD into module    !
    ! variables for local routines.                                       !
    !---------------------------------------------------------------------!
    Subroutine fission_wrapper(DELTAE,XLOCMX,V2_MIN,ITRMAX,NTHETA)

      Integer(ipr), INTENT(IN) :: ITRMAX,NTHETA
      Real(pr), INTENT(IN) :: DELTAE,XLOCMX,V2_MIN

      iteration_max = ITRMAX
      number_angles = NTHETA

      energy_window = DELTAE
      loc_max = XLOCMX
      occ_min = V2_MIN

    End Subroutine fission_wrapper

    !---------------------------------------------------------------------!
    !  This routine copies arrays coming from HFODD into local arrays. In !
    !  addition, it sorts q.p. states by the occupation of their tail: if !
    !  the q.p. \mu is localized to the right (left), the tail is the     !
    !  extent of the wave-function to the left (right) of the neck.       !
    !---------------------------------------------------------------------!
    Subroutine copy_fromHFODD(ICHARG,iter)

      Integer(ipr), INTENT(IN) :: ICHARG,iter

      Real(pr) :: EQPISO,VQPISO,ESPEQU,DELEQU,FERMFN

      Integer(ipr) :: IREVER
      Integer(ipr) :: mu,ii

      Real(pr) :: tail_mu, Nmu
      Real(pr), Allocatable :: table_score(:)

      COMMON                                              &
             /FER_IS/ FERMFN(1:NDBASE,0:NDISOS)
      COMMON                                              &
             /EQUISO/ ESPEQU(1:NDBASE,0:NDREVE,0:NDISOS), &
                      DELEQU(1:NDBASE,0:NDREVE,0:NDISOS)
      COMMON                                              &
             /QUAISO/ EQPISO(1:NDBASE,0:NDREVE,0:NDISOS), &
                      VQPISO(1:NDBASE,0:NDREVE,0:NDISOS)

      IREVER=0

      ! Rank q.p. states by tail size. Routine sort_vector sorts in ascending
      ! order, and we want to make sure that states with the largest occupation
      ! and lowest tail size (or localization) come first. Vector 'table' keeps
      ! track of original ordering
      Allocate(table_score(1:Nstate))
      ii = 0
      Do mu=1,Nstate
         ii = ii + 1
         table(ii) = ii
         ! occupation at T>=0 (0 <= vmu <= 1)
         Nmu = v_left(mu,ICHARG) + vright(mu,ICHARG)
         ! size of the tail at T>=0 (0 <= tail_mu <= 1)
         tail_mu = localisation(mu,ICHARG) * (1.0_pr - Nmu)
         table_score(ii) = tail_mu
      End Do
      Call sort_vector(table_score, table, ii)
      Deallocate(table_score)

      ! Storing the q.p. ordered by tail size
      Do mu=1,Nstate
         table_Eqp(mu)      = EQPISO(table(mu),IREVER,ICHARG)
         table_Esp(mu)      = ESPEQU(table(mu),IREVER,ICHARG)
         table_loc(mu)      = localisation(table(mu),ICHARG)
         table_vlef(mu)     = v_left(table(mu),ICHARG)
         table_vrig(mu)     = vright(table(mu),ICHARG)
         table_fT(mu)       = FERMFN(table(mu),ICHARG)
         table_gT(mu)       = 1.0_pr - FERMFN(table(mu),ICHARG)
         table_flag(mu)     = flag_left(table(mu),IREVER,ICHARG)
         table_occ(mu)      = v_left(table(mu),ICHARG) + vright(table(mu),ICHARG)
         off_diagonal_f(mu) = diagonal_off_f(table(mu),ICHARG)
         off_diagonal_g(mu) = diagonal_off_g(table(mu),ICHARG)
         off_diagonal_E(mu) = diagonal_off_E(table(mu),ICHARG)
         If(debug.Ge.3) Then
            If(Abs(table_occ(mu)).Gt.1.D-5) Then
               Write(6,'("mu=",i4," vleft=",f10.7," vright=",f10.7," lmu=",f10.7," fT=",f10.7)') &
                          mu,table_vlef(mu),table_vrig(mu),table_loc(mu),table_fT(mu)
            End If
         End If
      End Do

    End Subroutine copy_fromHFODD

    !---------------------------------------------------------------------!
    !  This subroutine copies back from this module the q.p. energies,    !
    !  s.p. energies and left/right occupation into the proper common     !
    !  blocks. Array VQPISO is not updated, since (i) it is not used      !
    !  after the calls to these routines (ii) it does not contain the     !
    !  temperature-dependent Fermi-Dirac factors. If the occupation of    !
    !  q.p. are needed beyond this routine, the arrays v_left and v_right !
    !  should be used preferably.                                         !
    !---------------------------------------------------------------------!
    Subroutine copy_toHFODD(ICHARG)

      Integer(ipr), INTENT(IN) :: ICHARG

      Real(pr) :: EQPISO,VQPISO,ESPEQU,DELEQU,FERMFN

      Integer(ipr) :: mu

      COMMON                                              &
             /FER_IS/ FERMFN(1:NDBASE,0:NDISOS)
      COMMON                                              &
             /EQUISO/ ESPEQU(1:NDBASE,0:NDREVE,0:NDISOS), &
                      DELEQU(1:NDBASE,0:NDREVE,0:NDISOS)
      COMMON                                              &
             /QUAISO/ EQPISO(1:NDBASE,0:NDREVE,0:NDISOS), &
                      VQPISO(1:NDBASE,0:NDREVE,0:NDISOS)

      Do mu=1,Nstate
         EQPISO(table(mu),0,ICHARG)=table_Eqp(mu)
         EQPISO(table(mu),1,ICHARG)=table_Eqp(mu)
         ESPEQU(table(mu),0,ICHARG)=table_Esp(mu)
         ESPEQU(table(mu),1,ICHARG)=table_Esp(mu)
         localisation(table(mu),ICHARG)=table_loc(mu)
         v_left(table(mu),ICHARG)=table_vlef(mu)
         vright(table(mu),ICHARG)=table_vrig(mu)
         FERMFN(table(mu),ICHARG)=table_fT(mu)
         diagonal_off_f(table(mu),ICHARG)=off_diagonal_f(mu)
         diagonal_off_g(table(mu),ICHARG)=off_diagonal_g(mu)
         diagonal_off_E(table(mu),ICHARG)=off_diagonal_E(mu)
         flag_left(table(mu),0,ICHARG)=table_flag(mu)
         flag_left(table(mu),1,ICHARG)=table_flag(mu)
         flag_right(table(mu),0,ICHARG)=1.0_pr-table_flag(mu)
         flag_right(table(mu),1,ICHARG)=1.0_pr-table_flag(mu)
      End Do

    End Subroutine copy_toHFODD

    !---------------------------------------------------------------------!
    !  This subroutine finds all the pairs of quasi-particles that fit a  !
    !  certain set of criteria, and returns the list of indexes mu and nu !
    !  of these pairs. Original list of q.p. was ordered by tail size.    !
    !---------------------------------------------------------------------!
    Subroutine find_pairs_new(EFERMI,n_pair,iter)

      Integer(ipr), INTENT(IN) :: iter
      Integer(ipr), INTENT(INOUT) :: n_pair
      Real(pr), INTENT(IN) :: EFERMI

      Logical :: valid_partners,old_pair,available
      Integer(ipr) :: mu,nu,phnu,phmu,i_pair,ii,mu_act,nu_act
      Integer(ipr), Allocatable :: taken_pair(:,:)
      Integer(ipr), Allocatable :: index_mu(:),index_nu(:)
      Integer(ipr), Allocatable :: temp_mu(:),temp_nu(:)

      Real(pr) :: Enu,espn,xlnu,Nnu,fnumu
      Real(pr) :: Emu,espm,xlmu,Nmu,fmunu

      Real(pr), Allocatable :: table_score(:)

      Type(liste_int), pointer :: liste_pair, current_node, previous_node

      ! Nullify all points used in this routine
      NULLIFY(liste_pair)
      NULLIFY(current_node)
      NULLIFY(previous_node)

      i_pair = 1
      Allocate(liste_pair)
      liste_pair%numero   = i_pair
      liste_pair%value_mu = 1
      liste_pair%value_nu = 1
      liste_pair%score    = 0.0_pr
      NULLIFY(liste_pair%next)

      Allocate(taken_pair(1:Nstate,1:Nstate)); taken_pair(:,:) = 0

      ! Loop over q.p. states
      Do mu = 2, Nstate
         Do nu = 1, mu-1

            Enu  = table_Eqp(nu) ! Eqp (nu)
            espn = table_Esp(nu) ! Esp (nu)
            xlnu = table_loc(nu) ! l_nu
            Nnu  = table_occ(nu) ! N_nu
            fnumu= off_diagonal_f(nu)

                               phnu=0
            If(espn.Lt.EFERMI) phnu=1

            Emu  = table_Eqp(mu) ! Eqp (mu)
            espm = table_Esp(mu) ! Esp (mu)
            xlmu = table_loc(mu) ! l_mu
            Nmu  = table_occ(mu) ! N_mu
            fmunu= off_diagonal_f(mu)

                               phmu=0
            If(espm.Lt.EFERMI) phmu=1

            mu_act = table(mu); nu_act = table(nu)

            ! Considering only pairs of q.p. that (i) are close in energy, (ii) are occupied, (iii) are
            ! delocalized, (iv) have not yet been listed, and (v) have not been rotated yet at a previous
            ! iteration. This last condition is necessary to deal with the case T>0, where the generalized
            ! density is no longer diagonal after rotation. In this implementation, pairs at iteration n+1
            ! are either new q.p. partners, or old pairs. Therefore, the number of pairs can only become
            ! larger. This guarantees that the arrays table_pair_[] filled out in minimize_tails() will
            ! contain all necessary information to leave the full density matrix invariant.
            !
            ! Note: This is an arbitrary limitation that ensures that the generalized density matrix is only
            ! ----  block-diagonal with at most 2x2 blocks in the diagonal. If we allowed schemes such that
            !       (n1,n2) and (n3,n4) at iteration k become (n1,n3), (n2,n4) at iteration k+1, the bookkeeping
            !       would be a lot more complicated, since the full density matrix would have a 4x4 block
            !       along the diagonal at iteration k+1.
            valid_partners = Abs(Emu-Enu).Le.energy_window .And. (Nmu.Gt.occ_min .And. Nnu.Gt.occ_min) .And. &
                            (xlmu.Lt.loc_max .And. xlnu.Lt.loc_max) .And.( (phnu+phmu).Eq.0 .Or. (phnu+phmu).Eq.2 ) .And. &
                            taken_pair(mu,nu).Eq.0 .And. ( Abs(fmunu).Lt.1.E-14_pr .And. Abs(fnumu).Le.1.0E-14_pr )
            old_pair = (nu_act.Eq.table_partner(mu_act) .And. mu_act.Eq.table_partner(nu_act))
            available = .Not. (table_partner(mu_act).Gt.0 .Or. table_partner(nu_act).Gt.0)

            ! Considering only pairs with similar energy where both q.p. are (i) occupied and
            ! (ii) delocalized
            If( (valid_partners .And. available) .Or. old_pair) Then

               ! Make sure we won't consider this pair again
               taken_pair(mu,nu)=1

               ! Record the information on the new pair, including its score for delocalization
               i_pair = i_pair + 1

               Allocate(current_node)
               current_node%numero   = i_pair
               current_node%value_mu = mu
               current_node%value_nu = nu
               current_node%score    = (xlmu**2+xlnu**2)**2
               current_node%next     => liste_pair
               liste_pair => current_node

               If(debug.Ge.5) Then
                  Write(6,'("In find_pairs_new() - Pair ",i10)') i_pair-1
                  If(valid_partners) Write(6,'("   valid partners")')
                  If(old_pair) Write(6,'("   old_pair")')
                  Write(6,'("   mu=",i4," nu=",i4," mu_act=",i4," nu_act=",i4)') mu,nu,mu_act,nu_act
                  Write(6,'("   Emu=",f10.5," Enu=",f10.5)') Emu,Enu
                  Write(6,'("   Nmu=",f10.5," Nnu=",f10.5)') Nmu,Nnu
                  Write(6,'("   lmu=",f10.5," lnu=",f10.5)') xlmu,xlnu
               End If

            End If

         End Do ! End nu
      End Do ! End mu
      If(debug.Ge.1) Write(6,'("Total number of possible pairs: ",i10)') i_pair-1

      ! Now that the number of pairs is known, we store information in vectors
      Allocate(index_mu(1:i_pair),index_nu(1:i_pair),table_score(1:i_pair))

      Allocate(current_node); current_node => liste_pair

      ii = 0

      Do While ( associated(current_node) )

         ii = ii + 1
         index_mu(ii) = current_node%value_mu
         index_nu(ii) = current_node%value_nu
         table_score(ii) = current_node%score

         Allocate(previous_node)
         previous_node => current_node
         current_node  => current_node%next
         Deallocate(previous_node)

      End Do

      ! Remove spurious first point of the linked list
      n_pair=i_pair-1
      Allocate(temp_mu(1:n_pair),temp_nu(1:n_pair))
      Do ii=1,n_pair
         temp_mu(ii)=index_mu(ii)
         temp_nu(ii)=index_nu(ii)
      End Do
      Deallocate(index_mu,index_nu,table_score)
      If(Allocated(final_mu)) Deallocate(final_mu,final_nu); Allocate(final_mu(1:n_pair),final_nu(1:n_pair))
      Do ii=1,n_pair
         final_mu(ii)=temp_mu(ii)
         final_nu(ii)=temp_nu(ii)
      End Do

      Deallocate(temp_mu,temp_nu)
      Deallocate(taken_pair)

    End Subroutine find_pairs_new

    !---------------------------------------------------------------------!
    !  This subroutine minimizes the tails between two fragments by       !
    !  trying to maximize the localization of highly-delocalized pairs    !
    !  of quasi-particles.                                                !
    !---------------------------------------------------------------------!
    Subroutine minimize_tails(EFERMI,NZHERM,NZMAXX,ZPOINT,ICHARG,SFACTO)

      Integer(ipr), INTENT(IN) :: NZHERM,NZMAXX,ICHARG
      Real(pr), INTENT(IN) :: EFERMI,ZPOINT
      Real(pr), Allocatable, INTENT(IN) :: SFACTO(:,:)

      Integer(ipr) :: IFTEMP

      Integer(ipr) :: iter,i_pair,j_pair,n_pair,Ntiles,tiling,Ntilings,different,tiling_optimal
      Integer(ipr) :: size_tiling,back_rotate,finalize,mu,nu,mu_act,nu_act,i
#if(USE_MPI==1)
      Integer(ipr) :: mpi_rank,mpi_err
#endif

      Integer(ipr), Allocatable :: partner(:)
      Real(pr) :: tail_left,tail_right,tail_size
      Real(pr), Allocatable :: tail_vector(:)

      Type(PairSet) :: new_pair_vector
      Type(ListePairSet) :: tilings,tilings_unique

      COMMON &
             /T_FLAG/ IFTEMP

      ! Allocate and initialize necessary arrays for the localization
      Call allocate_arrays(ICHARG)

      Do iter=1,iteration_max

         ! Sort qp by tail size and initializes/updates arrays table_Eqp,
         ! table_Esp, table_loc, table_vlef, table_vrig, table_fT, table_gT,
         ! table_flag, table_occ, table
         Call copy_fromHFODD(ICHARG,iter)

         ! Computing all possible pairs
         n_pair=0; Call find_pairs_new(EFERMI,n_pair,iter)

         If(n_pair.Le.0) Then
#if(USE_MPI==1)
            Call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, mpi_err)
            Write(6,'("WARNING! mpi_rank= ",i8," ICHARG=",i1," iter=",i1," n_pair= ",i10," - Exiting ...")') &
                                mpi_rank,ICHARG,iter,n_pair
#else
            Write(6,'("WARNING! ICHARG=",i1," iter=",i1," n_pair= ",i10," - Exiting loop...")') &
                                ICHARG,iter,n_pair
#endif
            size_tiling = 0
            number_of_pairs(ICHARG) = size_tiling
            Call copy_toHFODD(ICHARG)
            ! Deallocate necessary arrays
            Call deallocate_arrays()
            Return

         Else

            If(debug.Ge.5) Write(6,'("ICHARG=",i1," iter=",i1," - Creating a PairSet...")') ICHARG,iter

            ! Create a new object 'list of pairs'
            Call PairSet_new(n_pair, new_pair_vector)
            Do i_pair=1,n_pair
               new_pair_vector%pair_vector(i_pair)%mu = final_mu(i_pair)
               new_pair_vector%pair_vector(i_pair)%nu = final_nu(i_pair)
            End Do
            If(debug.Ge.5) Then
               Write(6,'("ICHARG=",i1," iter=",i1," - List of all pairs")') ICHARG,iter
               Call PairSet_print(n_pair, new_pair_vector)
            End If

            If(n_pair.Gt.1) Then
               Ntilings=n_pair*(n_pair-1) ! maximum number of possible tilings
               Ntiles=n_pair ! maximum number of tiles in a tiling
            Else
               Ntilings=1 ! maximum number of possible tilings
               Ntiles=1 ! maximum number of tiles in a tiling
            End If
            If(debug.Ge.1) Write(6,'("ICHARG=",i1," iter=",i1," Ntilings=",i10," Ntiles=",i10," (max)")') &
                                      ICHARG,iter,Ntilings,Ntiles

            ! Create a new object 'list of lists of pairs' with default sizes
            Call ListePairSet_new(Ntilings, Ntiles, tilings)

            ! Identify all the possible tilings: the number of tilings is updated after this call
            Call tilings_find_all(n_pair, Ntilings, tilings, new_pair_vector)
            If(debug.Ge.7) Then
               Write(6,'("ICHARG=",i1," iter=",i1," Total number of mappings: ",i4)') ICHARG,iter,Ntilings
               Write(6,'(" ---- ")')
               Write(6,'("List of all mappings")')
               Call ListePairSet_print(Ntilings, tilings)
               Write(6,'(" ---- ")')
            End If

            ! Sanity check: there should be tilings to deal with
            If(Ntilings.Le.0) Then
#if(USE_MPI==1)
               Call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, mpi_err)
               Write(6,'("mpi_rank= ",i8," ICHARG=",i1," iter=",i1," Ntilings= ",i10)') &
                          mpi_rank,ICHARG,iter,Ntilings
#else
               Write(6,'("ICHARG=",i1," iter=",i1," Ntilings= ",i10)') &
                          ICHARG,iter,Ntilings
#endif
               Stop 'Error in minimize_pairs - no tiling left after tilings_find_all()'
            End If

            ! Remove all the duplicate tilings
            Call tilings_remove_duplicates(Ntiles, Ntilings, different, tilings, tilings_unique)
            If(debug.Ge.3) Then
               Write(6,'("ICHARG=",i1," iter=",i1," Total number of UNIQUE mappings: ",i4)') ICHARG,iter,different
               Write(6,'(" ---- ")')
               Write(6,'("List of UNIQUE mappings")')
               Call ListePairSet_print(different, tilings_unique)
               Write(6,'(" ---- ")')
            End If

            ! Sanity check: at least one tiling remains after removing duplicates
            If(different.Le.0) Then
#if(USE_MPI==1)
               Call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, mpi_err)
               Write(6,'("mpi_rank= ",i8," ICHARG=",i1," iter=",i1," different= ",i10)') &
                          mpi_rank,ICHARG,iter,different
#else
               Write(6,'("ICHARG= ",i2," different= ",i10)') ICHARG,different
#endif
               Stop 'Error in minimize_pairs - no tiling left after tilings_remove_duplicates()'
            End If

            ! Reallocate space for list of pairs based on tilings
            If(Allocated(final_mu)) Deallocate(final_mu,final_nu)
            Allocate(final_mu(1:Ntiles),final_nu(1:Ntiles))

            ! Rotate q.p. for all tilings
            Allocate(tail_vector(1:different))
            Do tiling=1,different

               size_tiling = tilings_unique%size_tiling(tiling)
               Do i_pair=1,size_tiling
                  final_mu(i_pair) = tilings_unique%liste(tiling)%pair_vector(i_pair)%mu
                  final_nu(i_pair) = tilings_unique%liste(tiling)%pair_vector(i_pair)%nu
               End Do

               If(size_tiling.Le.0) Then
#if(USE_MPI==1)
                  Call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, mpi_err)
                  Write(6,'("WARNING! mpi_rank= ",i8," ICHARG=",i1," iter=",i1," size_tiling= ",i10)') &
                                      mpi_rank,ICHARG,iter,size_tiling
#else
                  Write(6,'("WARNING! ICHARG=",i1," iter=",i1," size_tiling= ",i10)') &
                                      ICHARG,iter,size_tiling
#endif
               Else

                  If(debug.Ge.3) Then
                     Write(6,'("ICHARG=",i1," iter=",i1," tiling:",i6," number of pairs:",i6," (max:",i6,")")') &
                                ICHARG,iter,tiling,size_tiling,Ntiles
                  End If
                  tail_left = 0.0_pr; tail_right = 0.0_pr; back_rotate=1
                  Call rotate_qp(ICHARG,SFACTO,size_tiling,tail_left,tail_right,back_rotate)
                  If(debug.Ge.3) Then
                     Write(6,'("  Size of the Tails")')
                     Write(6,'("     Left: ",f20.14," right: ",f20.14)') tail_left,tail_right
                  End If

                  tail_vector(tiling) = tail_left+tail_right

               End If

            End Do ! End tiling

            ! Determine the optimal tiling
            tail_size = 100000.0_pr; tiling_optimal = 1
            Do tiling=1,different
               If(tail_vector(tiling).Lt.tail_size) Then
                  tiling_optimal = tiling
                  tail_size = tail_vector(tiling)
               End If
            End Do
            Deallocate(tail_vector)

            ! Re-do q.p. rotation for optimal tiling and keep rotated q.p.
            size_tiling = tilings_unique%size_tiling(tiling_optimal)
            Do i_pair=1,size_tiling
               final_mu(i_pair) = tilings_unique%liste(tiling_optimal)%pair_vector(i_pair)%mu
               final_nu(i_pair) = tilings_unique%liste(tiling_optimal)%pair_vector(i_pair)%nu
               ! Correspondence between partners should be established based on real indices, not
               ! indices mu and nu corresponding to qp ranked by tail size (since ranking may change
               ! from iteration to iteration)
               mu_act = table(final_mu(i_pair))
               nu_act = table(final_nu(i_pair))
               table_partner(mu_act) = nu_act
               table_partner(nu_act) = mu_act
            End Do

            If(debug.Ge.1) Then
               Write(6,'("------------------------------")')
               Write(6,'("QP ROTATION FOR OPTIMAL TILING")')
               Write(6,'("  ICHARG=",i1," iter=",i1)') ICHARG,iter
               Write(6,'("  Optimal tiling definition - Numero:",i6,"/",i6," Size:",i6," (max:",i6,")")') &
                                                        tiling_optimal,different,size_tiling,Ntiles
            End If
            tail_left = 0.0_pr; tail_right = 0.0_pr; back_rotate=0
            Call rotate_qp(ICHARG,SFACTO,size_tiling,tail_left,tail_right,back_rotate)
            If(debug.Ge.1) Then
               Write(6,'("  Left tail  : ",f20.14)') tail_left
               Write(6,'("  Right tail : ",f20.14)') tail_right
               Write(6,'("------------------------------")')
            End If

            ! Close everything
            If(associated(tilings%liste)) Call ListePairSet_del(Ntilings, tilings)
            If(associated(tilings_unique%liste)) Call ListePairSet_del(different, tilings_unique)
            If(associated(new_pair_vector%pair_vector)) Call PairSet_del(new_pair_vector)

         End If ! end n_pair=0

         ! Computing all localizations for this iteration, including left and
         ! right occupations. Updates arrays 'flag_left', 'flag_right',
         ! 'localization', 'v_left' and 'vright'
         If(IFTEMP.Eq.0) Then
            Call wave_localization(NZHERM,NZMAXX,ZPOINT,ICHARG,SFACTO)
         Else
            finalize = 1
            Call wave_localization_T(NZHERM,NZMAXX,ZPOINT,ICHARG,SFACTO,size_tiling,iter,finalize)
         End If

         ! Saving information about rotated qp pairs after each iteration of
         ! localization method
         !   - table_iter ........: at which iteration the pair was rotated
         !   - table_ipair .......: mapping between regular qp indexing and pair indexing at given iteration
         !   - temp_order_qp .....: mapping between tailsize indexing and regular indexing
         !   - table_singles_X ...: which qp are rotated (regular indexing)
         !   - temp_pair_mu ......: index of rotated qp (tailsize indexing)
         !   - temp_pair_fTmu ....: diagonal term of the statistical occupation matrix (upper part of the generalized
         !                          density matrix) for member mu of the pair (mu,nu)
         !   - temp_pair_fTmunu ..: off-diagonal term of statistical occupation matrix
         !   - temp_pair_gTmu ....: g is defined as g=1-f (lower part of the generalized density matrix)
         !   - temp_pair_gTmunu ..: off-diagonal term of g
         !   - temp_pair_Emu .....: diagonal term of the HFB matrix for member mu of the pair (mu,nu)
         !   - temp_pair_Emunu ...: off-diagonal term of the HFB matrix
         temp_order_qp(1:Nstate,iter) = table(1:Nstate)
         Do i_pair=1,size_tiling
            mu = final_mu(i_pair)
            nu = final_nu(i_pair)
            mu_act = table(mu)
            nu_act = table(nu)
            ! Remember at which iteration this pair was found
            table_iter(mu_act,ICHARG) = iter
            table_iter(nu_act,ICHARG) = iter
            table_ipair(mu_act,ICHARG) = i_pair
            table_ipair(nu_act,ICHARG) = i_pair
            table_singles(mu_act,ICHARG) = 0
            table_singles(nu_act,ICHARG) = 0
            ! Characteristics of the pair
            temp_pair_mu(i_pair,iter) = mu
            temp_pair_nu(i_pair,iter) = nu
            temp_pair_mu_vlef(i_pair,iter) = table_vlef(mu)
            temp_pair_mu_vrig(i_pair,iter) = table_vrig(mu)
            temp_pair_mu_loc(i_pair,iter)  = table_loc(mu)
            temp_pair_mu_flag(i_pair,iter) = table_flag(mu)
            temp_pair_nu_vlef(i_pair,iter) = table_vlef(nu)
            temp_pair_nu_vrig(i_pair,iter) = table_vrig(nu)
            temp_pair_nu_loc(i_pair,iter)  = table_loc(nu)
            temp_pair_nu_flag(i_pair,iter) = table_flag(nu)
            temp_pair_fTmu(i_pair,iter)   = table_fT(mu)
            temp_pair_fTnu(i_pair,iter)   = table_fT(nu)
            temp_pair_fTmunu(i_pair,iter) = off_diagonal_f(mu)
            temp_pair_gTmu(i_pair,iter)   = table_gT(mu)
            temp_pair_gTnu(i_pair,iter)   = table_gT(nu)
            temp_pair_gTmunu(i_pair,iter) = off_diagonal_g(mu)
            temp_pair_Emu(i_pair,iter)    = table_Eqp(mu)
            temp_pair_Enu(i_pair,iter)    = table_Eqp(nu)
            temp_pair_Emunu(i_pair,iter)  = off_diagonal_E(mu)
         End Do

      End Do ! End iteration

      i_pair=0
      Allocate(partner(1:Nstate))
      partner(:)=1
      ! Loop over all qp to identify if it belongs to a pair, and if yes, at
      ! which iteration it was included
      Do i=1,Nstate
         ! Current qp belongs to a pair
         If(table_singles(i,ICHARG).Eq.0.And.partner(i).Eq.1) Then
            i_pair = i_pair + 1
            ! get the iteration number
            iter = table_iter(i,ICHARG)
            ! get the pair number at this iteration
            j_pair = table_ipair(i,ICHARG)
            ! get the id of the partner and remove it from list
            mu = temp_pair_mu(j_pair,iter)
            nu = temp_pair_nu(j_pair,iter)
            mu_act = temp_order_qp(mu,iter)
            nu_act = temp_order_qp(nu,iter)
            If(debug.Ge.3) Then
               Write(6,'("Found pair number",i4)') i_pair
               Write(6,'("     mu=",i4," nu=",i4," mu_act=",i4," nu_act=",i4," iter=",i2)') mu,nu,mu_act,nu_act,iter
            End If
            If(mu_act.Ne.i .And. nu_act.Ne.i) Then
               Write(6,'("Error in minimize_pairs() - Inconsistency in pair indices")')
               Write(6,'("i=",i4," mu_act",i4," nu_act",i4)') i,mu_act,nu_act
            End If
            partner(mu_act) = 0
            partner(nu_act) = 0
            ! fill out the pair characteristics
            table_pair_order_qp(mu,ICHARG)   = mu_act
            table_pair_order_qp(nu,ICHARG)   = nu_act
            table_pair_mu(i_pair,ICHARG)     = mu
            table_pair_nu(i_pair,ICHARG)     = nu
            table_pair_mu_vlef(i_pair,ICHARG) = temp_pair_mu_vlef(j_pair,iter)
            table_pair_mu_vrig(i_pair,ICHARG) = temp_pair_mu_vrig(j_pair,iter)
            table_pair_mu_loc(i_pair,ICHARG)  = temp_pair_mu_loc(j_pair,iter)
            table_pair_mu_flag(i_pair,ICHARG) = temp_pair_mu_flag(j_pair,iter)
            table_pair_nu_vlef(i_pair,ICHARG) = temp_pair_nu_vlef(j_pair,iter)
            table_pair_nu_vrig(i_pair,ICHARG) = temp_pair_nu_vrig(j_pair,iter)
            table_pair_nu_loc(i_pair,ICHARG)  = temp_pair_nu_loc(j_pair,iter)
            table_pair_nu_flag(i_pair,ICHARG) = temp_pair_nu_flag(j_pair,iter)
            table_pair_fTmu(i_pair,ICHARG)   = temp_pair_fTmu(j_pair,iter)
            table_pair_fTnu(i_pair,ICHARG)   = temp_pair_fTnu(j_pair,iter)
            table_pair_fTmunu(i_pair,ICHARG) = temp_pair_fTmunu(j_pair,iter)
            table_pair_gTmu(i_pair,ICHARG)   = temp_pair_gTmu(j_pair,iter)
            table_pair_gTnu(i_pair,ICHARG)   = temp_pair_gTnu(j_pair,iter)
            table_pair_gTmunu(i_pair,ICHARG) = temp_pair_gTmunu(j_pair,iter)
            table_pair_Emu(i_pair,ICHARG)    = temp_pair_Emu(j_pair,iter)
            table_pair_Enu(i_pair,ICHARG)    = temp_pair_Enu(j_pair,iter)
            table_pair_Emunu(i_pair,ICHARG)  = temp_pair_Emunu(j_pair,iter)
         End If
      End Do
      number_of_pairs(ICHARG) = i_pair
      If(debug.Ge.1) Then
         Write(6,'("TOTAL NUMBER OF PAIRS FOR ICHARG=",i1," number_of_pairs=",i4,/,"---")') ICHARG,i_pair
      End If
      Deallocate(partner)

      ! Deallocate necessary arrays
      Call deallocate_arrays()

    End Subroutine minimize_tails

    !---------------------------------------------------------------------!
    !  This subroutine rotates the q.p. wave-functions for all pairs of   !
    !  q.p. in the given tiling. This routine is called within a loop     !
    !  over each tiling to determine the tiling that yields the smallest  !
    !  tails; it is called a second time for the optimal tiling.          !
    !---------------------------------------------------------------------!
    Subroutine rotate_qp(ICHARG,SFACTO,size_tiling,tail_left,tail_right,back_rotate)

      Use SAVQUA

      Integer(ipr), INTENT(IN) :: ICHARG,size_tiling,back_rotate
      Real(pr), Allocatable, INTENT(IN) :: SFACTO(:,:)
      Real(pr), INTENT(INOUT) :: tail_left,tail_right

      Integer(ipr) :: IFTEMP,LDBASE

      Integer(ipr) :: IREVER,LREVER
      Integer(ipr) :: itheta,ithetaMax,i_pair,i,i1,i2,iprint,full_pair
      Integer(ipr) :: mu,nu,mu_act,nu_act,icount
      Integer(ipr), Allocatable :: taken(:)

      Real(pr) :: xlnu,vnu,v1nu,v2nu,Emu,fmu,gmu
      Real(pr) :: xlmu,vmu,v1mu,v2mu,Enu,fnu,gnu
      Real(pr) :: vnew1mu,vnew2mu,xlnew_mu,vnew1nu,vnew2nu,xlnew_nu,Emunu,fmunu,gmunu
      Real(pr) :: omega_p_B,omega_t_B,t1,t2,max_dev
      Real(pr) :: Pi,DeltaMax,ct,st,diff,N1,N2
      Real(pr), Dimension(1:number_angles) :: theta,delta

      Complex(pr) :: C_UNIT,C_ZERO,c1,c2
      Complex(pr), Allocatable :: vec_mu(:),vec_nu(:)
      Complex(pr), Allocatable :: Umat(:,:),Vmat(:,:),dmat_low(:,:),dmat_upp(:,:),temp1(:,:),temp2(:,:)
      Complex(pr), Allocatable :: DensityMatrix_bef(:,:),DensityMatrix_aft(:,:)

      COMMON                 &
             /T_FLAG/ IFTEMP
      COMMON                 &
             /DIMENS/ LDBASE

      IREVER=0; LREVER=1

      C_UNIT=cmplx(1.0_pr,0.0_pr)
      C_ZERO=cmplx(0.0_pr,0.0_pr)

      Allocate(vec_mu(1:LDBASE),vec_nu(1:LDBASE))

      ! Defining the rotation angle
      Pi = 4.0_pr*Atan(1.0_pr)
      Do itheta=1,number_angles
         theta(itheta) = Real(itheta-1)*2.0_pr*Pi/Real(number_angles-1)
      End Do

      ! Test that the rotation followed by back-rotation has left the total
      ! U and V unchanged, as it should have - Storing U and V before the action
      If(debug.Ge.2) Then
         Allocate(Umat(1:NDBASE,1:NDSTAT),Vmat(1:NDBASE,1:NDSTAT))
         Umat(1:NDBASE,1:NDSTAT)=ASVQUA(1:NDBASE,1:NDSTAT,LREVER,ICHARG)
         Vmat(1:NDBASE,1:NDSTAT)=BSVQUA(1:NDBASE,1:NDSTAT,IREVER,ICHARG)
         Allocate(dmat_low(1:NDSTAT,1:NDSTAT),dmat_upp(1:NDSTAT,1:NDSTAT))
         Call rotated_R(ICHARG,size_tiling)
         dmat_low(1:NDSTAT,1:NDSTAT) = densityR_lower(1:NDSTAT,1:NDSTAT)
         dmat_upp(1:NDSTAT,1:NDSTAT) = densityR_upper(1:NDSTAT,1:NDSTAT)
      End If

      ! Computing the bulk occupation on the left and right, i.e. *not* including
      ! any of the q.p. that are in the various pairs. For each tiling, the
      ! number of pairs 'size_tiling' is different, as are the vectors
      ! 'final_mu' and 'final_nu'
      Allocate(taken(1:Nstate)); taken(:)=0
      Do i_pair=1,size_tiling
         mu=final_mu(i_pair); nu=final_nu(i_pair)
         taken(mu)=1; taken(nu)=1
      End Do
      Do mu=1,Nstate
         If(taken(mu).Eq.0) Then
            If(table_flag(mu).Eq.1.0_pr) Then
               tail_left = tail_left + table_vrig(mu)
            Else
               tail_right = tail_right + table_vlef(mu)
            End If
         End If
      End Do
      ! Sanity check: computing the initial size of the tail
      If(debug.Ge.2) Then
         t1 = tail_left
         t2 = tail_right
         Do mu=1,Nstate
            If(taken(mu).Eq.1) Then
               If(table_flag(mu).Eq.1.0_pr) Then
                  t1 = t1 + table_vrig(mu)
               Else
                  t2 = t2 + table_vlef(mu)
               End If
            End If
         End Do
         Write(6,'("  Initial Conditions")')
         Write(6,'("     NUMBER OF PAIRS: ",i8)') size_tiling
         Write(6,'("     INITIAL SIZE OF THE TAILS - Left: ",f20.14," Right: ",f20.14)') &
                    t1,t2
      End If
      Deallocate(taken)

      ! No point wasting time in this routine if there are no tilings
      If(size_tiling.Le.0) Return

      ! Loop over the number of pairs in a given tiling. For each pair, we
      ! search the angle theta that gives the maximum localization of the pair.
      ! We then compute the size of the tails for this angle and update the 2
      ! quantities 'tail_left' and 'tail_right'.
      Do i_pair=1,size_tiling

         mu = final_mu(i_pair)
         nu = final_nu(i_pair)

         ! Correspondence between the old indices (before sorting q.p. based on
         ! their tail size), mu_act, nu_act, and the new indices (after sorting),
         ! mu, nu. The old indices are needed because the U and V matrices have
         ! not been copied, so we must rotate the proper q.p. pairs.
         mu_act = table(mu)
         nu_act = table(nu)

         xlnu = table_loc(nu)  ! l_nu
         vnu  = table_occ(nu)  ! v^2
         v1nu = table_vlef(nu) ! v^2 for ]-infty, zn]
         v2nu = table_vrig(nu) ! v^2 for [zn, +infty[

         xlmu = table_loc(mu)  ! l_mu
         vmu  = table_occ(mu)  ! v^2
         v1mu = table_vlef(mu) ! v^2 for ]-infty, zn]
         v2mu = table_vrig(mu) ! v^2 for [zn, +infty[

         ! HFB matrix for this pair
         Enu  = table_Eqp(nu)  ! qp energy
         Emu  = table_Eqp(mu)  ! qp energy
         Emunu= off_diagonal_E(mu)
         ! Fermi-Dirac occupations for this pair
         fnu  = table_fT(nu)
         fmu  = table_fT(mu)
         fmunu= off_diagonal_f(mu)
         gnu  = table_gT(nu)
         gmu  = table_gT(mu)
         gmunu= off_diagonal_g(mu)

         ! Computing auxiliary summations. These quantities, \omega_mn(z) and
         ! \omega_mn(-\infty), are defined for each pair from the (U,V) before
         ! rotation, and should thus be constant for each value of the rotation
         ! angle. However, they should be recalculated after each iteration.
         If(IFTEMP.Eq.0) Then
            omega_p_B = omega_partial_B(ICHARG,SFACTO,mu_act,nu_act)
            omega_t_B = omega_full_B(ICHARG,mu_act,nu_act)
         End If

         If(debug.Ge.3) Then
            Write(6,'("     Pair ",i4,"/",i4)') i_pair,size_tiling
            Write(6,'("       mu=",i4," Emu=",f10.5," vmu=",f10.7," xlmu=",f10.7," fmu=",f10.7," fmunu=",f10.7)') &
                              mu,table_Eqp(mu),vmu,xlmu,fmu,fmunu
            Write(6,'("       nu=",i4," Enu=",f10.5," vnu=",f10.7," xlnu=",f10.7," fnu=",f10.7)') &
                              nu,table_Eqp(nu),vnu,xlnu,fnu
         End If

         If(debug.Ge.5) Write(6,'(" TESTS BEFORE: v=v1+v2 - mu: ",2f20.14," nu=",2f20.14)') &
                                    vmu,v1mu+v2mu,vnu,v1nu+v2nu

         ! Searching for the angle that maximizes the localization of the pair
         ! of q.p. After each q.p. has been rotated and the localization has
         ! been computed, we back-rotate the q.p. to restore the original U and
         ! V matrices for this tiling.
         Do itheta=1,number_angles

            ct=Cos(theta(itheta))
            st=Sin(theta(itheta))

            ! Rotation of the HFB eigenvectors
            vec_mu(1:LDBASE) = BSVQUA(1:LDBASE,mu_act,IREVER,ICHARG)
            vec_nu(1:LDBASE) = BSVQUA(1:LDBASE,nu_act,IREVER,ICHARG)
            Call wave_rotation(vec_mu,vec_nu,ct,st)
            BSVQUA(1:LDBASE,mu_act,IREVER,ICHARG) = vec_mu(1:LDBASE)
            BSVQUA(1:LDBASE,nu_act,IREVER,ICHARG) = vec_nu(1:LDBASE)
            vec_mu(1:LDBASE) = ASVQUA(1:LDBASE,mu_act,LREVER,ICHARG)
            vec_nu(1:LDBASE) = ASVQUA(1:LDBASE,nu_act,LREVER,ICHARG)
            Call wave_rotation(vec_mu,vec_nu,ct,st)
            ASVQUA(1:LDBASE,mu_act,LREVER,ICHARG) = vec_mu(1:LDBASE)
            ASVQUA(1:LDBASE,nu_act,LREVER,ICHARG) = vec_nu(1:LDBASE)

            ! Left and right occupations, and localization, for q.p. mu and nu after
            ! rotation
            If(IFTEMP.Eq.0) Then
               vnew1mu = ct**2 * v1mu + st**2 * v1nu + st*ct*(omega_t_B - omega_p_B)
               vnew2mu = ct**2 * v2mu + st**2 * v2nu + st*ct* omega_p_B
               vnew1nu = ct**2 * v1nu + st**2 * v1mu - st*ct*(omega_t_B - omega_p_B)
               vnew2nu = ct**2 * v2nu + st**2 * v2mu - st*ct* omega_p_B
            Else
               full_pair=1; iprint=0
               If(debug.Ge.4) iprint=1
               Call occup_pair_T(ICHARG,SFACTO,mu,nu,mu_act,nu_act,fmu,fnu,fmunu,gmu,gnu,gmunu, &
                                               vnew1mu,vnew2mu,vnew1nu,vnew2nu,full_pair,iprint)
            End If
            xlnew_mu = Abs(vnew1mu - vnew2mu)/(vnew1mu + vnew2mu)
            xlnew_nu = Abs(vnew1nu - vnew2nu)/(vnew1nu + vnew2nu)

            ! Norm of the localization vector
            delta(itheta) = (vnew1mu - vnew2mu)**2 + (vnew1nu - vnew2nu)**2

            If(debug.Ge.7) Then
               Write(6,'("        theta =",f10.5," before: vm + vn = ",f14.9," - after: vm + vn = ",f14.9)') &
                            theta(itheta)*180.0_pr/Pi, v1mu+v2mu+v1nu+v2nu,vnew1mu+vnew2mu+vnew1nu+vnew2nu
               Write(6,'("        omega_mu_nu(z) =",f14.8," omega_mu_nu(infty) = ",f14.8," Delta = ",f20.14)') &
                               omega_p_B,omega_t_B,delta(itheta)
            End If

            ! Back-rotation of the HFB eigenvectors to return to original
            ! vectors
            vec_mu(1:LDBASE) = BSVQUA(1:LDBASE,mu_act,IREVER,ICHARG)
            vec_nu(1:LDBASE) = BSVQUA(1:LDBASE,nu_act,IREVER,ICHARG)
            Call wave_rotation(vec_mu,vec_nu,ct,-st)
            BSVQUA(1:LDBASE,mu_act,IREVER,ICHARG) = vec_mu(1:LDBASE)
            BSVQUA(1:LDBASE,nu_act,IREVER,ICHARG) = vec_nu(1:LDBASE)
            vec_mu(1:LDBASE) = ASVQUA(1:LDBASE,mu_act,LREVER,ICHARG)
            vec_nu(1:LDBASE) = ASVQUA(1:LDBASE,nu_act,LREVER,ICHARG)
            Call wave_rotation(vec_mu,vec_nu,ct,-st)
            ASVQUA(1:LDBASE,mu_act,LREVER,ICHARG) = vec_mu(1:LDBASE)
            ASVQUA(1:LDBASE,nu_act,LREVER,ICHARG) = vec_nu(1:LDBASE)

         End Do ! End of loop over itheta

         ! Determining the angle that maximizes the localization of the pairs
         DeltaMax = 0.0_pr
         ithetaMax = 1
         Do itheta=1,number_angles
            If(delta(itheta).Gt.DeltaMax) Then
               DeltaMax=delta(itheta)
               ithetaMax=itheta
            End If
         End Do

         ! Repeat the rotation of the HFB eigenvectors at the right angle
         ! (index 'ithetaMax') to update the size of the tails, 'tail_left'
         ! and 'tail_right'.
         ct=Cos(theta(ithetaMax))
         st=Sin(theta(ithetaMax))

         ! Rotation of the HFB eigenvectors
         vec_mu(1:LDBASE) = BSVQUA(1:LDBASE,mu_act,IREVER,ICHARG)
         vec_nu(1:LDBASE) = BSVQUA(1:LDBASE,nu_act,IREVER,ICHARG)
         Call wave_rotation(vec_mu,vec_nu,ct,st)
         BSVQUA(1:LDBASE,mu_act,IREVER,ICHARG) = vec_mu(1:LDBASE)
         BSVQUA(1:LDBASE,nu_act,IREVER,ICHARG) = vec_nu(1:LDBASE)
         vec_mu(1:LDBASE) = ASVQUA(1:LDBASE,mu_act,LREVER,ICHARG)
         vec_nu(1:LDBASE) = ASVQUA(1:LDBASE,nu_act,LREVER,ICHARG)
         Call wave_rotation(vec_mu,vec_nu,ct,st)
         ASVQUA(1:LDBASE,mu_act,LREVER,ICHARG) = vec_mu(1:LDBASE)
         ASVQUA(1:LDBASE,nu_act,LREVER,ICHARG) = vec_nu(1:LDBASE)

         ! Left and right occupations, and localization, for q.p. mu and nu
         ! after rotation
         If(IFTEMP.Eq.0) Then
            vnew1mu = ct**2 * v1mu + st**2 * v1nu + st*ct*(omega_t_B - omega_p_B)
            vnew2mu = ct**2 * v2mu + st**2 * v2nu + st*ct* omega_p_B
            vnew1nu = ct**2 * v1nu + st**2 * v1mu - st*ct*(omega_t_B - omega_p_B)
            vnew2nu = ct**2 * v2nu + st**2 * v2mu - st*ct* omega_p_B
         Else
            full_pair=1; iprint=0
            Call occup_pair_T(ICHARG,SFACTO,mu,nu,mu_act,nu_act,fmu,fnu,fmunu,gmu,gnu,gmunu, &
                                            vnew1mu,vnew2mu,vnew1nu,vnew2nu,full_pair,iprint)
         End If
         xlnew_mu = Abs(vnew1mu - vnew2mu)/(vnew1mu + vnew2mu)
         xlnew_nu = Abs(vnew1nu - vnew2nu)/(vnew1nu + vnew2nu)

         ! Update of table_[] vector only after final rotation (needed to print)
         If(back_rotate.Eq.0) Then
            ! For the current pair of q.p. mu, nu, store the proper rotation
            ! angle (same for mu and nu)
            table_theta(mu)=theta(ithetaMax)
            ! HFB matrix
            table_Eqp(mu) = Emu*ct**2 + Enu*st**2 + 2.0_pr*st*ct*Emunu
            table_Eqp(nu) = Enu*ct**2 + Emu*st**2 - 2.0_pr*st*ct*Emunu
            off_diagonal_E(mu) = -st*ct*(Emu-Enu) + (ct**2-st**2)*Emunu
            off_diagonal_E(nu) = off_diagonal_E(mu)
            ! Fermi-Dirac occupation matrix
            table_fT(mu)  = fmu*ct**2 + fnu*st**2 + 2.0_pr*st*ct*fmunu
            table_fT(nu)  = fnu*ct**2 + fmu*st**2 - 2.0_pr*st*ct*fmunu
            off_diagonal_f(mu) = -st*ct*(fmu-fnu) + (ct**2-st**2)*fmunu
            off_diagonal_f(nu) = off_diagonal_f(mu)
            table_gT(mu)  = gmu*ct**2 + gnu*st**2 + 2.0_pr*st*ct*gmunu
            table_gT(nu)  = gnu*ct**2 + gmu*st**2 - 2.0_pr*st*ct*gmunu
            off_diagonal_g(mu) = -st*ct*(gmu-gnu) + (ct**2-st**2)*gmunu
            off_diagonal_g(nu) = off_diagonal_g(mu)
            ! Occupations and localization for q.p. mu
            table_loc(mu) = xlnew_mu
            table_vlef(mu)= vnew1mu
            table_vrig(mu)= vnew2mu
            ! Occupations and localization for q.p. nu
            table_loc(nu) = xlnew_nu
            table_vlef(nu)= vnew1nu
            table_vrig(nu)= vnew2nu
         End If

         If(vnew1mu.Gt.0.5_pr*(vnew1mu+vnew2mu)) Then
            tail_left = tail_left + vnew2mu
         Else
            tail_right = tail_right + vnew1mu
         End If
         If(vnew1nu.Gt.0.5_pr*(vnew1nu+vnew2nu)) Then
            tail_left = tail_left + vnew2nu
         Else
            tail_right = tail_right + vnew1nu
         End If

         If(debug.Ge.4.And.back_rotate.Eq.1) Then
            Write(6,'("       OPTIMAL ANGLE THETA: ",f10.5)') theta(ithetaMax)*180.0_pr/Pi
            Write(6,'("       FINAL SIZE OF THE TAILS - Left: ",f20.14," Right: ",f20.14)') &
                                tail_left, tail_right
            Write(6,'("       CONSERVATION OF OCCUPATION Nmu - old:",f17.14," new:",f17.14)') &
                                table_occ(mu),vnew1mu+vnew2mu
            Write(6,'("                                  Nnu - old:",f17.14," new:",f17.14)') &
                                table_occ(nu),vnew1nu+vnew2nu
            Write(6,'("       LOCALIZATION lmu - old:",f17.14," new:",f17.14)') xlmu,xlnew_mu
            Write(6,'("                    lnu - old:",f17.14," new:",f17.14)') xlnu,xlnew_nu
            Write(6,'("       ENERGIES - DIAGONAL -> Emu - old:",f17.14," new:",f17.14)') Emu,Emu*ct**2+Enu*st**2
            Write(6,'("                              Enu - old:",f17.14," new:",f17.14)') Enu,Enu*ct**2+Emu*st**2
            Write(6,'("                  OFF-DIAGONAL -> Emunu:",f17.14)') -st*ct*(Emu-Enu)
            Write(6,'("       STATISTICS DIAGONAL -> fT_mu - old:",f17.14," new:",f17.14)') &
                                                     fmu,fmu*ct**2 + fnu*st**2 + 2.0_pr*st*ct*fmunu
            Write(6,'("                              fT_nu - old:",f17.14," new:",f17.14)') &
                                                     fnu,fnu*ct**2 + fmu*st**2 - 2.0_pr*st*ct*fmunu
            Write(6,'("                  OFF-DIAGONAL -> fmunu:",f17.14)') -st*ct*(fmu-fnu) + (ct**2-st**2)*fmunu
            Write(6,'("                  DIAGONAL -> gT_mu - old:",f17.14," new:",f17.14)')  &
                                                     gmu,gmu*ct**2 + gnu*st**2 + 2.0_pr*st*ct*gmunu
            Write(6,'("                              gT_nu - old:",f17.14," new:",f17.14)') &
                                                     gnu,gnu*ct**2 + gmu*st**2 - 2.0_pr*st*ct*gmunu
            Write(6,'("                  OFF-DIAGONAL -> gmunu:",f17.14)') -st*ct*(gmu-gnu) + (ct**2-st**2)*gmunu

         End If

         If(debug.Ge.1.And.back_rotate.Eq.0) Then
            Write(6,'("     Final Results (no back-rotation)")')
            Write(6,'("       OPTIMAL ANGLE THETA: ",f10.5)') theta(ithetaMax)*180.0_pr/Pi
            Write(6,'("       FINAL SIZE OF THE TAILS - Left: ",f20.14," Right: ",f20.14)') &
                                tail_left, tail_right
            Write(6,'("       CONSERVATION OF OCCUPATION Nmu - old:",f17.14," new:",f17.14)') &
                                table_occ(mu),vnew1mu+vnew2mu
            Write(6,'("                                  Nnu - old:",f17.14," new:",f17.14)') &
                                table_occ(nu),vnew1nu+vnew2nu
            Write(6,'("       LOCALIZATION lmu - old:",f17.14," new:",f17.14)') xlmu,xlnew_mu
            Write(6,'("                    lnu - old:",f17.14," new:",f17.14)') xlnu,xlnew_nu
            Write(6,'("       ENERGIES - DIAGONAL -> Emu - old:",f17.14," new:",f17.14)') Emu,table_Eqp(mu)
            Write(6,'("                              Enu - old:",f17.14," new:",f17.14)') Enu,table_Eqp(nu)
            Write(6,'("                  OFF-DIAGONAL -> Emunu:",f17.14)') off_diagonal_E(mu)
            Write(6,'("       STATISTICS DIAGONAL -> fT_mu - old:",f17.14," new:",f17.14)') fmu,table_fT(mu)
            Write(6,'("                              fT_nu - old:",f17.14," new:",f17.14)') fnu,table_fT(nu)
            Write(6,'("                  OFF-DIAGONAL -> fmunu:",f17.14)') off_diagonal_f(mu)
            Write(6,'("                  DIAGONAL -> gT_mu - old:",f17.14," new:",f17.14)') gmu,table_gT(mu)
            Write(6,'("                              gT_nu - old:",f17.14," new:",f17.14)') gnu,table_gT(nu)
            Write(6,'("                  OFF-DIAGONAL -> gmunu:",f17.14)') off_diagonal_g(mu)
         End If

         ! Back-rotation of the HFB eigenvectors to return to original vectors.
         ! This step is executed for each pair of a given tiling to restore the
         ! original U and V matrices; when the optimal tiling has been found,
         ! the routine is called one more time and this step is skipped so that
         ! the U and V matrices can be updated with the result of the rotation.
         If(back_rotate.Eq.1) Then
            vec_mu(1:LDBASE) = BSVQUA(1:LDBASE,mu_act,IREVER,ICHARG)
            vec_nu(1:LDBASE) = BSVQUA(1:LDBASE,nu_act,IREVER,ICHARG)
            Call wave_rotation(vec_mu,vec_nu,ct,-st)
            BSVQUA(1:LDBASE,mu_act,IREVER,ICHARG) = vec_mu(1:LDBASE)
            BSVQUA(1:LDBASE,nu_act,IREVER,ICHARG) = vec_nu(1:LDBASE)
            vec_mu(1:LDBASE) = ASVQUA(1:LDBASE,mu_act,LREVER,ICHARG)
            vec_nu(1:LDBASE) = ASVQUA(1:LDBASE,nu_act,LREVER,ICHARG)
            Call wave_rotation(vec_mu,vec_nu,ct,-st)
            ASVQUA(1:LDBASE,mu_act,LREVER,ICHARG) = vec_mu(1:LDBASE)
            ASVQUA(1:LDBASE,nu_act,LREVER,ICHARG) = vec_nu(1:LDBASE)
         End If

      End Do ! End pair

      ! Printing the off-diagonal term of the HFB matrix after rotation
      If(debug.Ge.1) Then
         If(back_rotate.Eq.0) Then
            icount=0
            Do i_pair=1,size_tiling
               mu = final_mu(i_pair)
               If(Abs(off_diagonal_E(mu)).Gt.1.D-14) Then
                  icount=icount+1
               End If
            End Do
            Write(6,'("Total number of non-zero off-diagonal terms: ",i4)') icount
         End If
      End If
      ! Test that the rotation followed by back-rotation has left the total
      ! U and V unchanged (within machine precision), as it should have
      If(debug.Ge.2) Then
         i1=0; i2=0
         Do nu=1,Nstate
            Do mu=1,LDBASE
               If(Abs(Vmat(mu,nu)-BSVQUA(mu,nu,IREVER,ICHARG)).Ge.1.D-14) i1=i1+1
               If(Abs(Umat(mu,nu)-ASVQUA(mu,nu,LREVER,ICHARG)).Ge.1.D-14) i2=i2+1
            End Do
         End Do
         If(back_rotate.Eq.1) Then
            Write(6,'("Numbers of non-equal elements after full cycle (back-rotation, should be 0): ",2i12)') i1,i2
         Else
            Write(6,'("Numbers of non-equal elements after full cycle (no back_rotation, should be non-zero): ",2i12)') i1,i2
            ! Count the number non-zero matrix elements in U and V
            i1=0; i2=0
            Do nu=1,Nstate
               Do mu=1,LDBASE
                  If(Abs(BSVQUA(mu,nu,IREVER,ICHARG)).Ge.1.D-14) i1=i1+1
                  If(Abs(ASVQUA(mu,nu,LREVER,ICHARG)).Ge.1.D-14) i2=i2+1
               End Do
            End Do
            Write(6,'("Number of non-zero elements in V: ",i12)') i1
            Write(6,'("Number of non-zero elements in U: ",i12)') i2
            ! Compute density matrices
            Allocate(DensityMatrix_bef(1:NDBASE,1:NDBASE),DensityMatrix_aft(1:NDBASE,1:NDBASE))
            DensityMatrix_bef(:,:)=C_ZERO
            DensityMatrix_aft(:,:)=C_ZERO
            If(IFTEMP.Eq.0) Then
               ! Compute density matrix in configuration space before rotation
               ! Remember: B = V*, A = U*
               CALL ZGEMM('N','C',LDBASE,LDBASE,NState, &
                                    C_UNIT,Vmat,NDBASE, &
                                           Vmat,NDBASE, &
                       C_ZERO,DensityMatrix_bef,NDBASE)
               ! Compute density matrix in configuration space after rotation
               CALL ZGEMM('N','C',LDBASE,LDBASE,NState, &
                           C_UNIT,BSVQUA(1,1,IREVER,ICHARG),NDBASE, &
                                  BSVQUA(1,1,IREVER,ICHARG),NDBASE, &
                                   C_ZERO,DensityMatrix_aft,NDBASE)
            Else
               ! Compute density matrix in configuration space before rotation
               Allocate(temp1(1:NDBASE,1:NDSTAT))
               Allocate(temp2(1:NDSTAT,1:NDBASE))
               temp1(:,:)=C_ZERO; temp2(:,:)=C_ZERO
               ! Part V* E V^T = B E B^+
               CALL ZGEMM('N','N',LDBASE,NState,NState, &
                                    C_UNIT,Vmat,NDBASE, &
                                       dmat_low,NDSTAT, &
                                   C_ZERO,temp1,NDBASE)
               CALL ZGEMM('N','C',LDBASE,LDBASE,NState, &
                                   C_UNIT,temp1,NDBASE, &
                                           Vmat,NDBASE, &
                       C_ZERO,DensityMatrix_bef,NDBASE)
               ! Part U F U^+ = A* F A^T = (F^+ A^T)^+ A^T
               CALL ZGEMM('C','T',Nstate,LDBASE,NState, &
                                C_UNIT,dmat_upp,NDSTAT, &
                                           Umat,NDBASE, &
                                   C_ZERO,temp2,NDSTAT)
               CALL ZGEMM('C','T',LDBASE,LDBASE,NState, &
                                   C_UNIT,temp2,NDSTAT, &
                                           Umat,NDBASE, &
                       C_UNIT,DensityMatrix_bef,NDBASE)
               ! Compute density matrix in configuration space after rotation
               Call rotated_R(ICHARG,size_tiling)
               i1=0;i2=0
               Do nu=1,Nstate
                  Do mu=1,Nstate
                     diff = Abs(densityR_lower(mu,nu)-dmat_low(mu,nu))
                     If(diff.Ge.1.D-14) Then
                        i1=i1+1
                        !Write(6,'("i=",i4," j=",i4," Emn(bef)=",2f20.14," Emn(aft)=",2f20.14)') &
                        !           mu,nu,dmat_low(mu,nu),densityR_lower(mu,nu)
                     End If
                     diff = Abs(densityR_upper(mu,nu)-dmat_upp(mu,nu))
                     If(diff.Ge.1.D-14) Then
                        i2=i2+1
                        !Write(6,'("i=",i4," j=",i4," Fmn(bef)=",2f20.14," Fmn(aft)=",2f20.14)') &
                        !           mu,nu,dmat_upp(mu,nu),densityR_upper(mu,nu)
                     End If
                  End Do
               End Do
               Write(6,'("Number of non equal elements in Emn before and after: ",i12)') i1
               Write(6,'("Number of non equal elements in Fmn before and after: ",i12)') i2
               ! Part V* E V^T = B E B^+
               Umat(1:NDBASE,1:NDSTAT)=ASVQUA(1:NDBASE,1:NDSTAT,LREVER,ICHARG)
               Vmat(1:NDBASE,1:NDSTAT)=BSVQUA(1:NDBASE,1:NDSTAT,IREVER,ICHARG)
               temp1(:,:)=C_ZERO; temp2(:,:)=C_ZERO
               CALL ZGEMM('N','N',LDBASE,NState,NState, &
                                    C_UNIT,Vmat,NDBASE, &
                                 densityR_lower,NDSTAT, &
                                   C_ZERO,temp1,NDBASE)
               CALL ZGEMM('N','C',LDBASE,LDBASE,NState, &
                                   C_UNIT,temp1,NDBASE, &
                                           Vmat,NDBASE, &
                       C_ZERO,DensityMatrix_aft,NDBASE)
               ! Part U F U^+ = A* F A^T = (F^+ A^T)^+ A^T
               CALL ZGEMM('C','T',Nstate,LDBASE,NState, &
                          C_UNIT,densityR_upper,NDSTAT, &
                                           Umat,NDBASE, &
                                   C_ZERO,temp2,NDSTAT)
               CALL ZGEMM('C','T',LDBASE,LDBASE,NState, &
                                   C_UNIT,temp2,NDSTAT, &
                                           Umat,NDBASE, &
                       C_UNIT,DensityMatrix_aft,NDBASE)
               Deallocate(temp1,temp2)
            End If ! IFTEMP=0
            ! Count the number of matrix elements that differ
            icount=0; max_dev=1.D-30; i1=0; i2=0
            c1=Cmplx(0.0_pr,0.0_pr); c2=Cmplx(0.0_pr,0.0_pr)
            Do nu=1,LDBASE
               c1=c1+DensityMatrix_bef(nu,nu)
               c2=c2+DensityMatrix_aft(nu,nu)
               Do mu=1,LDBASE
                  If(Abs(DensityMatrix_bef(mu,nu)).Ge.1.D-14) i1=i1+1
                  If(Abs(DensityMatrix_aft(mu,nu)).Ge.1.D-14) i2=i2+1
                  diff = Abs(DensityMatrix_aft(mu,nu)-DensityMatrix_bef(mu,nu))
                  If(diff.Ge.1.D-12) Then
                     icount=icount+1
                     If(diff.Ge.max_dev) Then
                        max_dev = diff
                     End If
                  End If
               End Do
            End Do
            Write(6,'("Trace(\rho_bef)=",2f20.14," Trace(\rho_aft)=",2f20.14)') c1,c2
            Write(6,'("Number of non-zero elements in \rho_bef: ",i12)') i1
            Write(6,'("Number of non-zero elements in \rho_aft: ",i12)') i2
            Write(6,'("Number of non equal elements in \rho before and after: ",i12)') icount
            Write(6,'("   --> Maximum deviation ............................: ",e21.14)') max_dev
            Deallocate(DensityMatrix_bef,DensityMatrix_aft)
         End If ! back_rotate=0
      End If ! debug=2

      If(debug.Ge.2) Deallocate(Umat,Vmat,dmat_low,dmat_upp)

      Deallocate(vec_mu,vec_nu)

    End Subroutine rotate_qp

    !---------------------------------------------------------------------!
    ! This function returns the quantity omega_\mu\nu(z) for q.p. states  !
    ! \mu and \nu and argument z = ZPOINT. This function is independent   !
    ! of the temperature T by definition.                                 !
    !---------------------------------------------------------------------!
    Real(pr) Function omega_partial_B(ICHARG,SFACTO,mu,nu)

      Use SAVQUA

      Integer(ipr), INTENT(IN) :: ICHARG,mu,nu
      Real(pr), Allocatable, INTENT(IN) :: SFACTO(:,:)

      Integer(ipr) :: NXVECT,NYVECT,NZVECT,LDBASE,IPHAPP,IPHAPM,IPHAMP,IPHAMM

      Integer(ipr) :: i,j,IREVER
      Integer(ipr) :: NX,NY,NZ,MX,MY,MZ
      Real(pr) :: cN

      COMMON                                              &
             /T_PHAS/ IPHAPP(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAPM(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAMP(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAMM(0:NDYMAX,0:NDYMAX,0:NDKART)
      COMMON                                                            &
             /BASISO/ NXVECT(1:NDBASE),NYVECT(1:NDBASE),NZVECT(1:NDBASE)
      COMMON                 &
             /DIMENS/ LDBASE

      IREVER=0
      cN = 0.0_pr
      Do j=1,LDBASE
         NX=NXVECT(j); NY=NYVECT(j); NZ=NZVECT(j)
         Do i=1,LDBASE
            MX=NXVECT(i); MY=NYVECT(i); MZ=NZVECT(i)
            ! Attention: BSVQUA contains V*
            If(NX.Eq.MX.And.NY.Eq.MY) Then
               cN = cN + Real((BSVQUA(i,mu,IREVER,ICHARG)*Conjg(BSVQUA(j,nu,IREVER,ICHARG)) &
                              +BSVQUA(i,nu,IREVER,ICHARG)*Conjg(BSVQUA(j,mu,IREVER,ICHARG)))&
                              *SFACTO(MZ,NZ)*IPHAPP(MY,NY,0))
            End If
         End Do
      End Do
      omega_partial_B = cN

    End Function omega_partial_B

    !---------------------------------------------------------------------!
    ! This function returns the quantity Omega_\mu\nu(-infty) for q.p.    !
    ! states \mu and \nu. This function is independent of the temperature !
    ! T by definition.                                                         !
    !---------------------------------------------------------------------!
    Real(pr) Function omega_full_B(ICHARG,mu,nu)

      Use SAVQUA

      Integer(ipr), INTENT(IN) :: ICHARG,mu,nu

      Integer(ipr) :: LDBASE

      Integer(ipr) :: i,IREVER,LREVER
      Complex(pr) :: cN

      COMMON                 &
             /DIMENS/ LDBASE

      IREVER=0; LREVER=1
      cN=Cmplx(0.0_pr,0.0_pr)
      Do i=1,LDBASE
         cN = cN + Conjg(BSVQUA(i,mu,IREVER,ICHARG))*BSVQUA(i,nu,IREVER,ICHARG)
      End Do

      omega_full_B = 2.0_pr*Real(cN)

      If(debug.Ge.2) Then
         Do i=1,LDBASE
            cN = cN + Conjg(ASVQUA(i,mu,LREVER,ICHARG))*ASVQUA(i,nu,LREVER,ICHARG)
         End Do
         Write(6,'("     ORTHONORMALIZATION: ",2f20.14," (should be zero)")') cN
      End If

    End Function omega_full_B

    !---------------------------------------------------------------------!
    ! This routine rotates the vectors vec_mu and vec_nu based on the     !
    ! value of cost(theta) and sin(theta) stored in ct and st.            !
    !         vec_mu -> +cos(theta)*vec_mu + sin(theta)*vec_nu            !
    !         vec_nu -> -sin(theta)*vec_mu + cos(theta)*vec_nu            !
    !---------------------------------------------------------------------!
    Subroutine wave_rotation(vec_mu,vec_nu,ct,st)

      Real(pr), INTENT(IN) :: ct,st
      Complex(pr), Allocatable, INTENT(INOUT) :: vec_mu(:),vec_nu(:)

      Integer(ipr) :: LDBASE

      Complex(pr) :: ct_c, st_c
      Complex(pr), Allocatable :: vec1(:),vec2(:)

      COMMON                 &
             /DIMENS/ LDBASE

      ct_c = Cmplx(ct,0.0_pr); st_c = Cmplx(st,0.0_pr)

      Allocate(vec1(1:LDBASE),vec2(1:LDBASE))
      vec1(1:LDBASE) = vec_mu(1:LDBASE); vec2(1:LDBASE) = vec_nu(1:LDBASE)

      Call ZSCAL(LDBASE,ct_c,vec_mu,1)
      Call ZAXPY(LDBASE,st_c,vec_nu,1,vec_mu,1)

      vec_nu(1:LDBASE) = vec2(1:LDBASE)
      Call ZSCAL(LDBASE, ct_c,vec_nu,1)
      Call ZAXPY(LDBASE,-st_c,vec1,1,vec_nu,1)

      Deallocate(vec1,vec2)

    End Subroutine wave_rotation

    !---------------------------------------------------------------------!
    ! This routine constructs the rotated generalized density matrix at   !
    ! finite temperature. The generalized density contains diagonal terms !
    ! corresponding to q.p. that have not been rotated, and 2x2 blocks    !
    ! for pairs of rotated q.p.                                           !
    !---------------------------------------------------------------------!
    Subroutine rotated_R(ICHARG,size_tiling)

      Integer(ipr), INTENT(IN) :: ICHARG,size_tiling

      Integer(ipr) :: mu,nu,i_pair,mu_act,nu_act,count

      densityR_upper(:,:)=Cmplx(0.0_pr,0.0_pr)
      densityR_lower(:,:)=Cmplx(0.0_pr,0.0_pr)
      Do mu=1,Nstate
         mu_act = table(mu)
         densityR_upper(mu_act,mu_act) = table_fT(mu)
         densityR_lower(mu_act,mu_act) = table_gT(mu)
      End Do
      Do i_pair=1,size_tiling
         mu = final_mu(i_pair)
         nu = final_nu(i_pair)
         mu_act = table(mu)
         nu_act = table(nu)
         densityR_upper(mu_act,nu_act) = off_diagonal_f(mu)
         densityR_upper(nu_act,mu_act) = off_diagonal_f(mu)
         densityR_lower(mu_act,nu_act) = off_diagonal_g(mu)
         densityR_lower(nu_act,mu_act) = off_diagonal_g(mu)
         If(debug.Ge.2) Then
            Write(6,'("In rotated_R()")')
            Write(6,'("   i_pair",i2," (of ",i2,") - mu_act=",i4," nu_act=",i4)') &
                       i_pair,size_tiling,mu_act,nu_act
            Write(6,'(21x,"fmu=",f10.7," fnu=",f10.7," fmunu=",f10.7)') &
                           table_fT(mu),table_fT(nu),off_diagonal_f(mu)
            Write(6,'(21x,"gmu=",f10.7," gnu=",f10.7," gmunu=",f10.7)') &
                           table_gT(mu),table_gT(nu),off_diagonal_g(mu)
         End If
      End Do
      If(debug.Ge.2) Then
         count=0
         Do nu=1,Nstate
            Do mu=1,Nstate
               If(Abs(densityR_lower(mu,nu)).Gt.1.D-14) count=count+1
            End Do
         End Do
         Write(6,'("Total number of non-zero elements in rotated R: ",i12," (size:",i6,")")') count,Nstate
      End If

    End Subroutine rotated_R
    !---------------------------------------------------------------------!
    !                                                                     !
    !   This subroutine computes the localization indicator of every q.p. !
    !   and returns arrays containing occupations in the left and right   !
    !   fragments, the localization indicator, and a tag left/right for   !
    !   each fragment. Left (right) fragment is denoted '1'  ('2') by     !
    !   convention.                                                       !
    !                                                                     !
    !   Outputs:                                                          !
    !     - flag_left ...: if 1, .q. is localized in the left fragment,   !
    !                      if 0, it is localized in the right fragment,   !
    !     - flag_right ..: same as flag_left, for the right fragment.     !
    !     - v_left ......: occupation of the q.p. in the left fragment at !
    !                      T>=0 (hence includes possible thermal effects).!
    !     - vright ......: same as v_left, for the right fragment.        !
    !     - localization : localization indicator of q.p.                 !
    !                                                                     !
    !   The sum v_left + v_right gives the total occupation of q.p. at    !
    !   T>=0; it is equivalent to VQPISO only at T=0.                     !
    !                                                                     !
    !---------------------------------------------------------------------!
    Subroutine wave_localization(NZHERM,NZMAXX,ZPOINT,ICHARG,SFACTO)

      Use SAVQUA

      Integer(ipr), INTENT(IN) :: NZHERM,NZMAXX,ICHARG
      Real(pr), Allocatable, INTENT(IN) :: SFACTO(:,:)
      Real(pr), INTENT(IN) :: ZPOINT

      Integer(ipr) :: NXVECT,NYVECT,NZVECT,LDBASE,NUMBQP,IFTEMP,IPHAPP,IPHAPM,IPHAMP,IPHAMM
      Real(pr) :: EQPISO,VQPISO,FERMFN

      Logical :: flip,asymmetry
      Integer(ipr) :: I,J,NX,NY,NZ,MX,MY,MZ,state,IREVER,KARTEZ,LREVER,i1,i2
      Real(pr) :: occup,occup1,norme_total,norme_left,norme_right,vtotal
      Real(pr), Allocatable :: D_MATR(:,:),D1MATR(:,:)
      Complex(pr) :: cN,cM

      COMMON &
             /T_PHAS/ IPHAPP(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAPM(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAMP(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAMM(0:NDYMAX,0:NDYMAX,0:NDKART)
      COMMON &
             /BASISO/ NXVECT(1:NDBASE),NYVECT(1:NDBASE),NZVECT(1:NDBASE)
      COMMON &
             /DIMENS/ LDBASE
      COMMON &
             /QPNUMS/ NUMBQP(0:NDREVE,0:NDISOS)
      COMMON &
             /T_FLAG/ IFTEMP
      COMMON &
             /FER_IS/ FERMFN(1:NDBASE,0:NDISOS)
      COMMON &
             /QUAISO/ EQPISO(1:NDBASE,0:NDREVE,0:NDISOS), &
                      VQPISO(1:NDBASE,0:NDREVE,0:NDISOS)

      IREVER=0; LREVER=1

      ! Count the number non-zero matrix elements in U and V
      If(debug.Ge.2) Then
         i1=0; i2=0
         Do j=1,NUMBQP(IREVER,ICHARG)
            Do i=1,LDBASE
               If(Abs(BSVQUA(i,j,IREVER,ICHARG)).Ge.1.D-14) i1=i1+1
               If(Abs(ASVQUA(i,j,LREVER,ICHARG)).Ge.1.D-14) i2=i2+1
            End Do
         End Do
         Write(6,'("Number of non-zero elements in V (wave_localization): ",i12)') i1
         Write(6,'("Number of non-zero elements in U (wave_localization): ",i12)') i2
      End If

!$OMP PARALLEL DO &
!$OMP DEFAULT(NONE) &
!$OMP SCHEDULE(STATIC) &
!$OMP SHARED(NUMBQP,IREVER,LREVER,ICHARG,LDBASE,BSVQUA,FERMFN,ASVQUA,IPHAPP,SFACTO,&
!$OMP        debug,NXVECT,NYVECT,NZVECT,flag_left,flag_right,vright,v_left,localisation) &
!$OMP PRIVATE(state,cN,J,vtotal,NX,NY,NZ,I,MX,MY,MZ,occup)
      Do state=1,NUMBQP(IREVER,ICHARG)

         ! Total occupation of q.p. \mu, N_mu, at finite temperature T>=0
         cN=Cmplx(0.0_pr,0.0_pr)
         Do J=1,LDBASE
            NY=NYVECT(J)
            cN = cN + ( Abs(BSVQUA(J,state,IREVER,ICHARG))**2 &
                            *(1.0_pr-FERMFN(state,ICHARG))    &
                       +Abs(ASVQUA(J,state,LREVER,ICHARG))**2 &
                                    *FERMFN(state,ICHARG) )*IPHAPP(NY,NY,0)
         End Do
         vtotal = Real(cN)

         ! Trace of the density matrix \rho_{ij} * d_ij(zN) for fragment '2'
         ! (z >= zN) gives occup = N_mu(fragment '2'); occupation in fragment
         ! '1' is deduced by subtraction from the total occupation vtotal
         cN=Cmplx(0.0_pr,0.0_pr)
         Do J=1,LDBASE
            NX=NXVECT(J); NY=NYVECT(J); NZ=NZVECT(J)
            Do I=1,LDBASE
               MX=NXVECT(I); MY=NYVECT(I); MZ=NZVECT(I)
               ! Attention: BSVQUA contains V*, and ASVQUA, U*.
               If(NX.Eq.MX.And.NY.Eq.MY) Then
                  cN = cN + (BSVQUA(I,state,IREVER,ICHARG)  &
                      *Conjg(BSVQUA(J,state,IREVER,ICHARG)) &
                             *(1.0_pr-FERMFN(state,ICHARG)) &
                      +Conjg(ASVQUA(I,state,LREVER,ICHARG)) &
                            *ASVQUA(J,state,LREVER,ICHARG)  &
                                     *FERMFN(state,ICHARG)) &
                             *SFACTO(MZ,NZ)*IPHAPP(MY,NY,0)
               End If
            End Do
         End Do
         occup = Real(cN)

         ! Displays occupations
         If(debug.Ge.5) Then
            Write(6,'("Total occupation: ",f20.14," partial occupation:",f20.14," fn=",f20.14)') &
                       vtotal,occup,FERMFN(state,ICHARG)
         End If

         ! Sanity checks: total occupation of q.p. \mu should be between 0 and 1
         If(occup.Lt.0.0_pr) Then
            occup=0.0_pr
         End If
         If(occup.Gt.vtotal) Then
            occup=vtotal
         End If

         ! 'left'  or '1' is by definition between -\infty and zN,
         ! 'right' or '2' is by definition between zN and +\infty
         flag_left(state,0,ICHARG)=0.0_pr; flag_right(state,0,ICHARG)=0.0_pr
         flag_left(state,1,ICHARG)=0.0_pr; flag_right(state,1,ICHARG)=0.0_pr
         If(occup.Lt.0.5_pr*vtotal) Then
            flag_left(state,0,ICHARG)=1.0_pr
            flag_left(state,1,ICHARG)=1.0_pr
         Else
            flag_right(state,0,ICHARG)=1.0_pr
            flag_right(state,1,ICHARG)=1.0_pr
         End If

         ! Occupations of right and left fragments. Attention: v_left and
         ! vright contain the contribution from the Fermi-Dira distributions by
         ! opposition to VQPISO, which does not
         vright(state,ICHARG)=occup
         v_left(state,ICHARG)=vtotal-occup

         ! By definition, localization is l = |N_1 - N_2| / N, and 0 <= l <= 1
         localisation(state,ICHARG)=Abs(vtotal - 2.0_pr*occup)/vtotal

      End Do
!$OMP END PARALLEL DO

      ! Debugging: Checking the precision of computing the left and right
      !            occupation separately versus deducing the left occupation
      !            from the total and the right (or vice-versa). Separate
      !            calculations require computing two "numerical" integrals
      !            while the full occupation of q.p. is analytical.
      If(debug.Ge.6) Then

         norme_total=0.0_pr; norme_left=0.0_pr; norme_right=0.0_pr

         Allocate(D_MATR(0:NZMAXX,0:NZMAXX))
         KARTEZ=3; flip = .False.; asymmetry = .False.
         Call define_asymmetry(NZHERM,NZMAXX,KARTEZ,D_MATR,ZPOINT,flip,asymmetry)
         Allocate(D1MATR(0:NZMAXX,0:NZMAXX))
         KARTEZ=3; flip = .True.; asymmetry = .False.
         Call define_asymmetry(NZHERM,NZMAXX,KARTEZ,D1MATR,ZPOINT,flip,asymmetry)

         ! Loop over q.p. states. quantity 'occup' measures the occupation in
         ! the right fragment ('2'), i.e. for z >= zN
         Do state=1,NUMBQP(IREVER,ICHARG)

            ! Total occupation of q.p. \mu, N_mu, at finite temperature T>=0
            cN=Cmplx(0.0_pr,0.0_pr)
            Do J=1,LDBASE
               NY=NYVECT(J)
               cN = cN + ( Abs(BSVQUA(J,state,IREVER,ICHARG))**2 &
                               *(1.0_pr-FERMFN(state,ICHARG))    &
                          +Abs(ASVQUA(J,state,LREVER,ICHARG))**2 &
                                       *FERMFN(state,ICHARG) )*IPHAPP(NY,NY,0)
            End Do
            vtotal= Real(cN)

            ! Left and right occupation of q.p. \mu at finite temperature T>=0
            cN=Cmplx(0.0_pr,0.0_pr); cM=Cmplx(0.0_pr,0.0_pr)
            Do J=1,LDBASE
               NX=NXVECT(J); NY=NYVECT(J); NZ=NZVECT(J)
               Do I=1,LDBASE
                  MX=NXVECT(I); MY=NYVECT(I); MZ=NZVECT(I)
               ! Attention: BSVQUA contains V*, and ASVQUA, U*.
                  If(NX.Eq.MX.And.NY.Eq.MY) Then
                     cN = cN + (BSVQUA(I,state,IREVER,ICHARG)  &
                         *Conjg(BSVQUA(J,state,IREVER,ICHARG)) &
                                *(1.0_pr-FERMFN(state,ICHARG)) &
                         +Conjg(ASVQUA(I,state,LREVER,ICHARG)) &
                               *ASVQUA(J,state,LREVER,ICHARG)  &
                                        *FERMFN(state,ICHARG) )&
                                        *D_MATR(MZ,NZ)*IPHAPP(MY,NY,0)
                     cM = cM + (BSVQUA(I,state,IREVER,ICHARG)  &
                         *Conjg(BSVQUA(J,state,IREVER,ICHARG)) &
                                *(1.0_pr-FERMFN(state,ICHARG)) &
                         +Conjg(ASVQUA(I,state,LREVER,ICHARG)) &
                               *ASVQUA(J,state,LREVER,ICHARG)  &
                                        *FERMFN(state,ICHARG) )&
                                        *D1MATR(MZ,NZ)*IPHAPP(MY,NY,0)
                  End If
               End Do
            End Do
            occup=Real(cN); occup1=Real(cM)

            ! Sanity check: numerical precision of integrations
            If(Abs(occup+occup1-vtotal).Gt.1.D-10) Then
               Write(6,'("wave_localization: Precision test - mu= ",i4," N_mu= ",f20.14," N_mu(left)+N_mu(right)= ",f20.14," error =",f20.14)') &
                          state,vtotal,occup+occup1,occup+occup1-vtotal
            End If

            ! Sanity check: norme_total gives the total number of particles, norme_left the
            !               total number of particles in the left fragment ]-\infty, zN[,
            !               norme_right the total number of particles in the right fragment
            !               [zN, +\infty[
            norme_total=norme_total+2.0_pr*vtotal
            norme_right=norme_right+2.0_pr*occup
            norme_left =norme_left +2.0_pr*occup1

         End Do

         Write(6,'("wave_localization: ICHARG = ",i2," Norme: ",f20.14," Left: ",f20.14," Right: ",f20.14)') &
                    ICHARG, norme_total, norme_left, norme_right

         Deallocate(D_MATR,D1MATR)

      End If

    End Subroutine wave_localization

    !---------------------------------------------------------------------!
    !                                                                     !
    !   This subroutine is the analog of wave_localization() above at T>0 !
    !   when quasiparticles are rotated. Based on a set of rotated qp     !
    !   with given statistical occupations, it computes the spatial       !
    !   occupations of all quasiparticles at T>0, their localization, and !
    !   defines the arrays flag_left, flag_right which identify to which  !
    !   fragment the qp belong. Compared to wave_localization(), this     !
    !   routine needs specific information about the qp pairs that have   !
    !   been rotated, and for which the generalized density is not        !
    !   diagonal anymore.                                                 !
    !                                                                     !
    !   Inputs:                                                           !
    !     - size_tiling ...: current number of pairs                      !
    !     - final_mu ......: index of quasiparticles mu and nu in the     !
    !       final_nu         tailsize indexing scheme                     !
    !     - table .........: mapping between tailsize and regular         !
    !                        indexing                                     !
    !     - table_fT ......: statistical occupation f (diagonal term)     !
    !     - off_diagonal_f : statistical occupation f (off-diagonal term) !
    !     - table_gT ......: g = 1-f (diagonal term)                      !
    !     - off_diagonal_g : localization indicator of q.p.               !
    !                                                                     !
    !   Outputs:                                                          !
    !     - flag_left ...: if 1, .q. is localized in the left fragment,   !
    !                      if 0, it is localized in the right fragment,   !
    !     - flag_right ..: same as flag_left, for the right fragment.     !
    !     - v_left ......: occupation of the q.p. in the left fragment at !
    !                      T>=0 (hence includes possible thermal effects).!
    !     - vright ......: same as v_left, for the right fragment.        !
    !     - localization : localization indicator of q.p.                 !
    !                                                                     !
    !       These outputs are produced thanks to the call to copy_toHFODD !
    !       at the end of the routine.                                    !
    !                                                                     !
    !---------------------------------------------------------------------!
    Subroutine wave_localization_T(NZHERM,NZMAXX,ZPOINT,ICHARG,SFACTO,size_tiling,iter,finalize)

      Use SAVQUA

      Integer(ipr), INTENT(IN) :: NZHERM,NZMAXX,ICHARG,size_tiling,iter,finalize
      Real(pr), Allocatable, INTENT(IN) :: SFACTO(:,:)
      Real(pr), INTENT(IN) :: ZPOINT

      Integer(ipr) :: i_pair,IREVER,LREVER
      Integer(ipr) :: mu,nu,mu_act,nu_act,full_pair,iprint
      Integer(ipr), Allocatable :: taken(:)
      Real(pr) :: fmu,fnu,fmunu,gmu,gnu,gmunu,occup,vtotal,vnew1mu,vnew2mu,vnew1nu,vnew2nu

      IREVER=0; LREVER=1

      If(debug.Ge.1) Write(6,'("Number of pairs at iteration iter=",i4," - size_tiling =",i4)') &
                                                             iter,size_tiling

      Allocate(taken(1:Nstate)); taken(:)=0
      Do i_pair=1,size_tiling
         mu=final_mu(i_pair); nu=final_nu(i_pair)
         taken(mu)=1; taken(nu)=1
      End Do

      ! Computing left, right and total occupation, as well as localization
      ! indicator of each rotated pair of q.p.
      full_pair=1; iprint=1
      Do i_pair=1,size_tiling
         mu = final_mu(i_pair)
         nu = final_nu(i_pair)
         mu_act = table(mu)
         nu_act = table(nu)
         ! Statistical occupations
         fmu   = table_fT(mu)
         fnu   = table_fT(nu)
         fmunu = off_diagonal_f(mu)
         gmu   = table_gT(mu)
         gnu   = table_gT(nu)
         gmunu = off_diagonal_g(mu)
         Call occup_pair_T(ICHARG,SFACTO,mu,nu,mu_act,nu_act,fmu,fnu,fmunu,gmu,gnu,gmunu,&
                                         vnew1mu,vnew2mu,vnew1nu,vnew2nu,full_pair,iprint)
         If(debug.Ge.3.And.iprint.Eq.1) Then
            Write(6,'("In wave_localization_T()")')
            Write(6,'("   i_pair",i2," (of ",i2,") - mu_act=",i4," nu_act=",i4)') &
                       i_pair,size_tiling,mu_act,nu_act
            Write(6,'(21x,"fmu=",f10.7," fnu=",f10.7," fmunu=",f10.7)') &
                           fmu,fnu,fmunu
            Write(6,'(21x,"gmu=",f10.7," gnu=",f10.7," gmunu=",f10.7)') &
                           gmu,gnu,gmunu
            Write(6,'(21x,"N1mu=",f10.7," N2mu=",f10.7," N1nu=",f10.7," N2nu=",f10.7)') vnew1mu,vnew2mu,vnew1nu,vnew2nu
         End If
         ! Updating the table_flag, table_vlef, table_vrig and table_loc arrays
         ! This should only be done at the very end of the localization process,
         ! before getting ready to pass this info back to HFODD
         If(finalize.Eq.1) Then
            ! qp mu
            occup=vnew2mu; vtotal=vnew1mu+vnew2mu
            table_flag(mu)=0.0_pr
            table_flag(mu)=0.0_pr
            If(occup.Lt.0.5_pr*vtotal) Then
               table_flag(mu)=1.0_pr
               table_flag(mu)=1.0_pr
            End If
            table_vlef(mu)=vnew1mu
            table_vrig(mu)=vnew2mu
            table_loc(mu)=Abs(vtotal - 2.0_pr*occup)/vtotal
            ! qp nu
            occup=vnew2nu; vtotal=vnew1nu+vnew2nu
            table_flag(nu)=0.0_pr
            table_flag(nu)=0.0_pr
            If(occup.Lt.0.5_pr*vtotal) Then
               table_flag(nu)=1.0_pr
               table_flag(nu)=1.0_pr
            End If
            table_vlef(nu)=vnew1nu
            table_vrig(nu)=vnew2nu
            table_loc(nu)=Abs(vtotal - 2.0_pr*occup)/vtotal
         End If
      End Do

      ! Computing left, right and total occupation, as well as localization
      ! indicator of qp not belonging to a rotated pair
      full_pair=0; iprint=1
      Do mu=1,Nstate
         If(taken(mu).Eq.0) Then
            mu_act = table(mu)
            fmu=table_fT(mu)
            fmunu=0.0_pr
            gmu=table_gT(mu)
            ! Quasiparticle \nu is not needed below but we still set its
            ! characteristics identical to mu
            nu_act = mu_act
            fnu=fmu
            gnu=gmu
            gmunu=0.0_pr
            Call occup_pair_T(ICHARG,SFACTO,mu,nu,mu_act,nu_act,fmu,fnu,fmunu,gmu,gnu,gmunu,&
                                            vnew1mu,vnew2mu,vnew1nu,vnew2nu,full_pair,iprint)
            If(debug.Ge.6.And.iprint.Eq.1) Then
               Write(6,'("State ",i4,"/",i6)') mu,Nstate
               Write(6,'("     mu=",i4," mu_act=",i4)') mu,mu_act
               Write(6,'("     fmu=",f10.7)') fmu
               Write(6,'("     N1mu=",f10.7," N2mu=",f10.7," N1nu=",f10.7," N2nu=",f10.7)') vnew1mu,vnew2mu,vnew1nu,vnew2nu
            End If
            If(finalize.Eq.1) Then
               ! qp mu
               occup=vnew2mu; vtotal=vnew1mu+vnew2mu
               table_flag(mu)=0.0_pr
               table_flag(mu)=0.0_pr
               If(occup.Lt.0.5_pr*vtotal) Then
                  table_flag(mu)=1.0_pr
                  table_flag(mu)=1.0_pr
               End If
               table_vlef(mu)=vnew1mu
               table_vrig(mu)=vnew2mu
               table_loc(mu)=Abs(vtotal - 2.0_pr*occup)/vtotal
            End If
         End If
      End Do

      ! We need to call this routine because wave_localization_T()
      ! deals with qp ordered by tail size, not normal order
      Call copy_toHFODD(ICHARG)

      Deallocate(taken)

    End Subroutine wave_localization_T

    !---------------------------------------------------------------------!
    !                                                                     !
    !   This subroutine computes the total occupation and fragment        !
    !   occupations of quasiparticles \mu and \nu at T>0.                 !
    !                                                                     !
    !     Inputs:                                                         !
    !       - mu, nu ........: pair indexes (tailsize indexing)           !
    !       - mu_act, nu_act : pair indexes (regular indexing)            !
    !       - fmu, fnu ......: diagonal terms of the (mu,nu) block of the !
    !                          generalized density (upper part)           !
    !       - gmu, gnu ......: diagonal terms of the (mu,nu) block of the !
    !                          generalized density (lower part, g=1-f)    !
    !       - fmunu, gmunu ..: off-diagonal terms of the generalized      !
    !                          density !                                  !
    !                                                                     !
    !     Outputs:                                                        !
    !       - vnew1mu, vnew2mu: occupations of q.p. mu in the left and    !
    !                           right fragments                           !
    !       - vnew1nu, vnew2nu: occupations of q.p. nu in the left and    !
    !                           right fragments                           !
    !                                                                     !
    !---------------------------------------------------------------------!
    Subroutine occup_pair_T(ICHARG,SFACTO,mu,nu,mu_act,nu_act,fmu,fnu,fmunu,gmu,gnu,gmunu, &
                                          vnew1mu,vnew2mu,vnew1nu,vnew2nu,full_pair,iprint)

      Use SAVQUA

      Integer(ipr), INTENT(IN) :: ICHARG,mu,nu,mu_act,nu_act,full_pair,iprint
      Real(pr), INTENT(IN) :: fmu,fnu,fmunu,gmu,gnu,gmunu
      Real(pr), Allocatable, INTENT(IN) :: SFACTO(:,:)
      Real(pr), INTENT(INOUT) :: vnew1mu,vnew2mu,vnew1nu,vnew2nu

      Integer(ipr) :: NXVECT,NYVECT,NZVECT,LDBASE,IPHAPP,IPHAPM,IPHAMP,IPHAMM
      Integer(ipr) :: I,J,NX,NY,NZ,MX,MY,MZ,IREVER,LREVER
      Real(pr) :: occup,vtotal
      Real(pr) :: cN

      COMMON &
             /T_PHAS/ IPHAPP(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAPM(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAMP(0:NDYMAX,0:NDYMAX,0:NDKART), &
                      IPHAMM(0:NDYMAX,0:NDYMAX,0:NDKART)
      COMMON &
             /BASISO/ NXVECT(1:NDBASE),NYVECT(1:NDBASE),NZVECT(1:NDBASE)
      COMMON &
             /DIMENS/ LDBASE

      IREVER=0; LREVER=1

      ! Computing the total occupation of q.p. m at T>0. It also depends on its
      ! pair partner n, and is defined as
      !   N_m = \sum_{j} Fmm Ajm Amj^{+} + Gmm Bjm^{*}Bmj
      !                + Fmn Ajm Anj^{+} + Gmn Bjm^{*}Bnj
      ! At the very first iteration of the localization method, Fmn = Gmn = 0
      cN=0.0_pr
      Do j=1,LDBASE
         NY=NYVECT(j)
         cN = cN + Real( &
                   ( gmu*Abs(BSVQUA(j,mu_act,IREVER,ICHARG))**2   &
                 +   fmu*Abs(ASVQUA(j,mu_act,LREVER,ICHARG))**2 ) &
                 *IPHAPP(NY,NY,0)                                      &
                 + ( gmunu*BSVQUA(j,mu_act,IREVER,ICHARG) *Conjg(BSVQUA(j,nu_act,IREVER,ICHARG))   &
                    +fmunu*Conjg(ASVQUA(j,mu_act,LREVER,ICHARG))*ASVQUA(j,nu_act,LREVER,ICHARG) ) &
                 *Real(full_pair)*IPHAPP(NY,NY,0) &
                       )
      End Do
      vtotal = cN
      ! Computing the occupation in the left fragment of qp number m
      cN=0.0_pr
      Do j=1,LDBASE
         NX=NXVECT(j); NY=NYVECT(j); NZ=NZVECT(j)
         Do i=1,LDBASE
            MX=NXVECT(i); MY=NYVECT(i); MZ=NZVECT(i)
            If(NX.Eq.MX.And.NY.Eq.MY) Then
               cN = cN + Real( &
                         ( gmu*BSVQUA(i,mu_act,IREVER,ICHARG) *Conjg(BSVQUA(j,mu_act,IREVER,ICHARG))  &
                          +fmu*Conjg(ASVQUA(i,mu_act,LREVER,ICHARG))*ASVQUA(j,mu_act,LREVER,ICHARG) ) &
                       *SFACTO(MZ,NZ)*IPHAPP(NY,NY,0)                                      &
                       + ( gmunu*BSVQUA(i,mu_act,IREVER,ICHARG) *Conjg(BSVQUA(j,nu_act,IREVER,ICHARG))  &
                          +fmunu*Conjg(ASVQUA(i,mu_act,LREVER,ICHARG))*ASVQUA(j,nu_act,LREVER,ICHARG) ) &
                       *Real(full_pair)*SFACTO(MZ,NZ)*IPHAPP(NY,NY,0) &
                              )
            End If
         End Do
      End Do
      occup = cN

      ! Sanity checks: total occupation of q.p. \mu should be between 0 and 1
      If(occup.Lt.0.0_pr) Then
         occup=0.0_pr
      End If
      If(occup.Gt.vtotal) Then
         occup=vtotal
      End If

      vnew2mu=occup
      vnew1mu=vtotal-occup

      ! If full_pair=1, we also compute quantities related to q.p. \nu
      If(full_pair.Eq.1) Then

         cN=0.0_pr
         Do j=1,LDBASE
            NY=NYVECT(j)
            cN = cN + Real( &
                      ( gnu*Abs(BSVQUA(j,nu_act,IREVER,ICHARG))**2   &
                       +fnu*Abs(ASVQUA(j,nu_act,LREVER,ICHARG))**2 ) &
                    *IPHAPP(NY,NY,0)                                      &
                    + ( gmunu*BSVQUA(j,nu_act,IREVER,ICHARG) *Conjg(BSVQUA(j,mu_act,IREVER,ICHARG))  &
                       +fmunu*Conjg(ASVQUA(j,nu_act,LREVER,ICHARG))*ASVQUA(j,mu_act,LREVER,ICHARG) ) &
                    *IPHAPP(NY,NY,0) &
                          )
         End Do
         vtotal = cN

         cN=0.0_pr
         Do j=1,LDBASE
            NX=NXVECT(j); NY=NYVECT(j); NZ=NZVECT(j)
            Do i=1,LDBASE
               MX=NXVECT(i); MY=NYVECT(i); MZ=NZVECT(i)
               If(NX.Eq.MX.And.NY.Eq.MY) Then
                  cN = cN + Real( &
                            ( gnu*BSVQUA(i,nu_act,IREVER,ICHARG) *Conjg(BSVQUA(j,nu_act,IREVER,ICHARG))  &
                             +fnu*Conjg(ASVQUA(i,nu_act,LREVER,ICHARG))*ASVQUA(j,nu_act,LREVER,ICHARG) ) &
                          *SFACTO(MZ,NZ)*IPHAPP(NY,NY,0)                                      &
                          + ( gmunu*BSVQUA(i,nu_act,IREVER,ICHARG) *Conjg(BSVQUA(j,mu_act,IREVER,ICHARG))  &
                             +fmunu*Conjg(ASVQUA(i,nu_act,LREVER,ICHARG))*ASVQUA(j,mu_act,LREVER,ICHARG) ) &
                          *SFACTO(MZ,NZ)*IPHAPP(NY,NY,0) &
                                )
               End If
            End Do
         End Do
         occup = cN

         If(occup.Lt.0.0_pr) Then
            occup=0.0_pr
         End If
         If(occup.Gt.vtotal) Then
            occup=vtotal
         End If

         vnew2nu=occup
         vnew1nu=vtotal-occup

      End If

    End Subroutine occup_pair_T

    !---------------------------------------------------------------------!
    !  This subroutine prints information on localization of the q.p.     !
    !---------------------------------------------------------------------!
    Subroutine print_localization(EFERMI,ICHARG)

      Integer(ipr), INTENT(IN) :: ICHARG
      Real(pr), INTENT(IN) :: EFERMI

      Integer(ipr) :: NFIPRI,NUMBQP
      Real(pr) :: ESPEQU,DELEQU,EQPISO,VQPISO,FERMFN

      Integer(ipr) :: i
      Real(pr) :: DUMMY
      Character(Len=8) :: NAMEPN(0:1)

      COMMON                                               &
             /EQUISO/ ESPEQU(1:NDBASE,0:NDREVE,0:NDISOS),  &
                      DELEQU(1:NDBASE,0:NDREVE,0:NDISOS)
      COMMON                                               &
             /QUAISO/ EQPISO(1:NDBASE,0:NDREVE,0:NDISOS),  &
                      VQPISO(1:NDBASE,0:NDREVE,0:NDISOS)
      COMMON                                   &
             /QPNUMS/ NUMBQP(0:NDREVE,0:NDISOS)
      COMMON                                   &
             /FER_IS/ FERMFN(1:NDBASE,0:NDISOS)
      COMMON                 &
             /CFIPRI/ NFIPRI

      NAMEPN(0)='NEUTRONS'
      NAMEPN(1)='PROTONS '
      DUMMY=0.0_pr

      Write(NFIPRI,'(79("*"),/,"*",77X,"*")')
      Write(NFIPRI,'("*",2X,"LOCALIZATION - FERMI LEVEL: ",F8.4,28X,A8,"   *")') EFERMI,NAMEPN(ICHARG)
      Write(NFIPRI,'("*",77X,"*",/,79("*"),/,"*",77X,"*")')
      Write(NFIPRI,'("*  NO)   EQP",7X,"ESP     DELTA   ENERGY   V2(LEFT)  V2OCCUP   LOCALIZATION  *")')
      Write(NFIPRI,'("*",77X,"*")')
      Do i=1,NUMBQP(0,ICHARG)
         If(VQPISO(i,0,ICHARG)*(1.0_pr-FERMFN(i,ICHARG)).Gt.0.01_pr) Then
            Write(NFIPRI,'("*",i4,")",4F9.3,3F10.6,6X,"*")') &
              i,EQPISO(i,0,ICHARG),ESPEQU(i,0,ICHARG),DELEQU(i,0,ICHARG),DUMMY,&
                           v_left(i,ICHARG), v_left(i,ICHARG)+vright(i,ICHARG),   &
                           localisation(i,ICHARG)
         End If
      End Do
      Write(NFIPRI,'("*",77X,"*")')

    End Subroutine print_localization

    !---------------------------------------------------------------------!
    !  This subroutine prints the localization of the q.p. after rotation !
    !  on the standard HFODD output (characterized by Fortran logical     !
    !  unit NFIPRI.                                                       !
    !---------------------------------------------------------------------!
    Subroutine print_rotation(EFERMI,ICHARG)

      Integer(ipr), INTENT(IN) :: ICHARG
      Real(pr), INTENT(IN) :: EFERMI

      Integer(ipr) :: NFIPRI

      COMMON                 &
             /CFIPRI/ NFIPRI

      Write(NFIPRI,'(79("*"),/,"*",77X,"*")')
      Write(NFIPRI,'("*",23X,"ROTATION OF Q.P. WAVE-FUNCTIONS",23X,"*")')
      Write(NFIPRI,'("*",77X,"*")')
      Write(NFIPRI,'("*  ENERGY WINDOW: ",F5.2,"    MAXIMUM LOCALIZATION: ",F6.4,"   MINIMUM V2: ",F6.4,"  *")') &
                         energy_window, loc_max, occ_min
      Write(NFIPRI,'("*",77X,"*")')

      CALL print_localization(EFERMI,ICHARG)

    End Subroutine print_rotation

 End Module hfodd_fission_rotated_qp


 ! ==================================================================== !
 !                                                                      !
 !         FRAGMENT ENERGIES AND INTERACTION ENERGY PACKAGE             !
 !                                                                      !
 ! ==================================================================== !

 !----------------------------------------------------------------------!
 !                                                                      !
 !  This module contains a number the routines that compute the energy  !
 !  of each fragment and the interaction energy between fragments. This !
 !  includes the Skyrme, Coulomb and pairing energies.                  !
 !                                                                      !
 !  Inputs:                                                             !
 !    - flag_left ...: the value of element number 'mu' of this array   !
 !                     indicates if the corresponding quasiparticle     !
 !                     belongs to the left fragment (=1) or the right   !
 !                     fragment (=0).                                   !
 !    - flag_right ..: analog of flag_left for the right fragment       !
 !                                                                      !
 !  Outputs:                                                            !
 !    - F_FLAG ......: the value of element number 'mu' indicates if    !
 !                     the corresponding quasiparticle is occupied      !
 !                     (=1.0) or not (=0.0). It is used in DENSHF to    !
 !                     compute densities in r-space for each fragment   !
 !                     based on qp occupations (in filter_density).     !
                                                                        !
 !  Inputs (implicit in DENSHF):                                        !
 !   - number_of_pairs ....: Number of rotated qp pairs                 !
 !   - table_singles ......: Array keeping track of which qp are        !
 !                           rotated                                    !
 !   - table_pair_mu ......: Index (tail size ordering) of q.p. mu      !
 !   - table_pair_nu       : and nu                                     !
 !   - table_pair_order_qp : Array giving the mapping between tailsize  !
 !                           ordering and regular ordering              !
 !   - table_pair_fTmu ....: Diagonal terms of the Fermi-Dirac          !
 !     table_pair_fTnu       occupation matrix for the qp pair          !
 !                           (mu,nu)                                    !
 !   - table_pair_fTmunu ..: Off-diagonal germ of the FD occupatio      !
 !                           matrix (which is symmetric)                !
 !   - table_pair_gTmu ....: By definition, g = 1-f                     !
 !     table_pair_gTnu                                                  !
 !     table_pair_gTmunu                                                !
 !   - table_pair_Emu .....: Diagonal term of the HFB matrix after      !
 !     table_pair_Enu        rotation for the qp pair (mu,nu)           !
 !   - table_pair_Emunu ...: Off-diagonal term of the HFB matrix        !
 !  These inputs are only needed when implementing qp rotation at       !
 !  finite temperature T>0                                              !
 !                                                                      !
 !----------------------------------------------------------------------!

 Module hfodd_fission_interaction

    Use hfodd_sizes
    Use hfodd_fission_precision
    Use hfodd_fission_rotated_qp

    Implicit None

    ! Interaction energies between the two fragments
    Real(pr), PUBLIC, SAVE :: ESKINT=0.0_pr,EP_INT=0.0_pr,EC_INT=0.0_pr,EC1INT=0.0_pr,EC2INT=0.0_pr, &
                              ECDINT_ex=0.0_pr,ECXINT_ex=0.0_pr
    ! Energies for fragment 1
    Real(pr), PUBLIC, SAVE :: ESKFR1=0.0_pr,EPAIR1=0.0_pr,ECOUL1=0.0_pr,ECOEX1=0.0_pr,EKIN1N=0.0_pr,EKIN1P=0.0_pr
    Real(pr), Dimension(0:NDKART), PUBLIC, SAVE :: TPSQA1
    ! Energies for fragment 2
    Real(pr), PUBLIC, SAVE :: ESKFR2=0.0_pr,EPAIR2=0.0_pr,ECOUL2=0.0_pr,ECOEX2=0.0_pr,EKIN2N=0.0_pr,EKIN2P=0.0_pr
    Real(pr), Dimension(0:NDKART), PUBLIC, SAVE :: TPSQA2

    ! Breakdown of Skyrme energy for each fragment
    Real(pr), Dimension(NDCOUP), PUBLIC, SAVE :: fis_CCPALL
    Real(pr), Dimension(NDCOUP), PUBLIC, SAVE :: fis_DCPALL,   fis_ECPALL
    Real(pr), Dimension(NDCOUP), PUBLIC, SAVE :: fis_DCPALL_1, fis_ECPALL_1,  fis_EISALL_1
    Real(pr), Dimension(NDCOUP), PUBLIC, SAVE :: fis_DCPALL_2, fis_ECPALL_2,  fis_EISALL_2
    Real(pr), Dimension(NDCOUP), PUBLIC, SAVE :: fis_DCPALL_12,fis_ECPALL_12, fis_EISALL_12
    Complex(pr), Dimension(NDCOUP), PUBLIC, SAVE :: fis_ECCALL,fis_ECCALL_1,fis_ECCALL_2,fis_ECCALL_12

    Logical, PUBLIC :: compute_all = .True., compute_CoulExc = .False., compute_pairing = .True.

    ! Backup densities
    Complex(pr), Allocatable, PUBLIC :: backup_DENRHO(:,:,:,:)
    Complex(pr), Allocatable, PUBLIC :: backup_DENTAU(:,:,:,:)
    Complex(pr), Allocatable, PUBLIC :: backup_DENLPR(:,:,:,:)
    Complex(pr), Allocatable, PUBLIC :: backup_DENDIV(:,:,:,:)
    Complex(pr), Allocatable, PUBLIC :: backup_PENRHO(:,:,:,:)
    Complex(pr), Allocatable, PUBLIC :: backup_DENSPI(:,:,:,:,:)
    Complex(pr), Allocatable, PUBLIC :: backup_DENKIS(:,:,:,:,:)
    Complex(pr), Allocatable, PUBLIC :: backup_DENGRR(:,:,:,:,:)
    Complex(pr), Allocatable, PUBLIC :: backup_DENLPS(:,:,:,:,:)
    Complex(pr), Allocatable, PUBLIC :: backup_DENROS(:,:,:,:,:)
    Complex(pr), Allocatable, PUBLIC :: backup_DENROC(:,:,:,:,:)
    Complex(pr), Allocatable, PUBLIC :: backup_DENCUR(:,:,:,:,:)
    Complex(pr), Allocatable, PUBLIC :: backup_DENSCU(:,:,:,:,:,:)
    Complex(pr), Allocatable, PUBLIC :: backup_DENDES(:,:,:,:,:,:)

 Contains

    !---------------------------------------------------------------------!
    !  Filtering the density matrix based on conditions on q.p. states    !
    !   - I_TYPE=0: keep all q.p.                                         !
    !   - I_TYPE=1: keep all q.p. from the left fragment only             !
    !   - I_TYPE=2: keep all q.p. from the right fragment only            !
    !   - I_TYPE=3: keep all q.p. with +lambda <= E <= -lambda            !
    !   - I_TYPE=4: keep all q.p. with E < +lambda or E > -lambda         !
    !   - I_TYPE=5: keep all q.p. from the left fragment only with        !
    !                    +lambda <= E <= -lambda                          !
    !   - I_TYPE=6: keep all q.p. from the left fragment only with        !
    !                   E < +lambda or E > -lambda                        !
    !   - I_TYPE=7: keep all q.p. from the right fragment only with       !
    !                   +lambda <= E <= -lambda                           !
    !   - I_TYPE=8: keep all q.p. from the right fragment only with       !
    !                   E < +lambda or E > -lambda                        !
    !---------------------------------------------------------------------!
    Subroutine filter_density(EFERMI,ICHARG,I_TYPE)

      Use FRAGFL

      Integer(ipr), INTENT(IN) :: ICHARG,I_TYPE
      Real(pr), INTENT(IN) :: EFERMI

      Integer(ipr) :: NUMBQP
      Real(pr) :: EQPISO,VQPISO

      Integer(ipr) :: IREVER,state

      COMMON                                               &
             /QUAISO/ EQPISO(1:NDBASE,0:NDREVE,0:NDISOS),  &
                      VQPISO(1:NDBASE,0:NDREVE,0:NDISOS)
      COMMON                                   &
             /QPNUMS/ NUMBQP(0:NDREVE,0:NDISOS)

      If(.Not.Allocated(F_FLAG)) Allocate(F_FLAG(1:NDSTAT,0:NDREVE,0:NDISOS))
      F_FLAG(:,:,ICHARG)=1.0_pr

      IREVER=0

      Select Case (I_TYPE)

      ! Filter 0: no q.p. is filtered out
      Case (0)
         F_FLAG(:,:,ICHARG)=1.0_pr

      ! Filter 1: all q.p. from the right fragment only are filtered out
      Case (1)
         F_FLAG(:,:,ICHARG)=flag_left(:,:,ICHARG)

      ! Filter 2: all q.p. from the left fragment only are filtered out
      Case (2)
         F_FLAG(:,:,ICHARG)=flag_right(:,:,ICHARG)

      ! Filter 3: all q.p. with energy outside the [-lambda,+lambda] window are filtered out
      Case (3)
         Do state=1,NUMBQP(IREVER,ICHARG)
            If(EQPISO(state,IREVER,ICHARG).Gt.Abs(EFERMI)) F_FLAG(state,IREVER,ICHARG)=0.0_pr
         End Do

      ! Filter 4: all q.p. with energy inside the [-lambda,+lambda] window are filtered out
      Case(4)
         Do state=1,NUMBQP(IREVER,ICHARG)
            If(EQPISO(state,IREVER,ICHARG).Le.Abs(EFERMI)) F_FLAG(state,IREVER,ICHARG)=0.0_pr
         End Do

      ! Filter 5: all q.p. from left fragment with energy outside the [-lambda,+lambda] window are filtered out
      Case(5)
         F_FLAG(:,:,ICHARG)=flag_left(:,:,ICHARG)
         Do state=1,NUMBQP(IREVER,ICHARG)
            If(EQPISO(state,IREVER,ICHARG).Gt.Abs(EFERMI)) F_FLAG(state,IREVER,ICHARG)=0.0_pr
         End Do

      ! Filter 6: all q.p. from left with energy inside the [-lambda,+lambda] window are filtered out
      Case(6)
         F_FLAG(:,:,ICHARG)=flag_left(:,:,ICHARG)
         Do state=1,NUMBQP(IREVER,ICHARG)
            If(EQPISO(state,IREVER,ICHARG).Le.Abs(EFERMI)) F_FLAG(state,IREVER,ICHARG)=0.0_pr
         End Do

      ! Filter 7: all q.p. from right fragment with energy outside the [-lambda,+lambda] window are filtered out
      Case(7)
         F_FLAG(:,:,ICHARG)=flag_right(:,:,ICHARG)
         Do state=1,NUMBQP(IREVER,ICHARG)
            If(EQPISO(state,IREVER,ICHARG).Gt.Abs(EFERMI)) F_FLAG(state,IREVER,ICHARG)=0.0_pr
         End Do

      ! Filter 8: all q.p. from right with energy inside the [-lambda,+lambda] window are filtered out
      Case(8)
         F_FLAG(:,:,ICHARG)=flag_right(:,:,ICHARG)
         Do state=1,NUMBQP(IREVER,ICHARG)
            If(EQPISO(state,IREVER,ICHARG).Le.Abs(EFERMI)) F_FLAG(state,IREVER,ICHARG)=0.0_pr
         End Do

      End Select

    End Subroutine filter_density

    !---------------------------------------------------------------------!
    !  Writing filtered densities on disk                                 !
    !---------------------------------------------------------------------!
    Subroutine open_density(filename, fortran_unit, error)

      Character(Len=68), INTENT(IN) :: filename
      Integer(ipr), INTENT(IN) :: fortran_unit
      Integer(ipr), INTENT(INOUT) :: error

      Open(Unit=fortran_unit,File=filename,Status='Unknown',Form='Formatted',Iostat=error)

    End Subroutine open_density

    ! type_density = [nucleus] + [localization] + [type_qp]
    !   nucleus ......: 'n', 'p'
    !   localization .: 'all', 'lef', 'rig'
    !   type_qp ......: 'continum', 'discrete', 'spectrum'
    Subroutine write_density(fortran_unit,density,continuum,NXHERM,NYHERM,NZHERM,type_density)

      Integer(ipr), INTENT(IN) :: fortran_unit,NXHERM,NYHERM,NZHERM
      Complex(pr), Allocatable, INTENT(IN) :: density(:,:,:),continuum(:,:,:)
      Character(Len=18), INTENT(IN) :: type_density

      Real(pr) :: EXPAUX,HOMSCA

      Integer(ipr) :: IX,IY,IZ

      COMMON                                        &
             /SCALNG/ HOMSCA(NDKART)
      COMMON                                        &
             /DENEXP/ EXPAUX(NDXHRM,NDYHRM,NDZHRM)

      Write(fortran_unit,'(a18)') type_density
      Write(fortran_unit,'(6X,3PE18.10)') HOMSCA(1),HOMSCA(2),HOMSCA(3)
      Write(fortran_unit,'(6X,3i10)') NXHERM, NYHERM, NZHERM
      Write(fortran_unit,'((6X,10(1PE18.10)))') &
           (((Real(density(IX,IY,IZ)*EXPAUX(IX,IY,IZ)),IX=1,NXHERM),IY=1,NYHERM),IZ=1,NZHERM)
      Write(fortran_unit,'((6X,10(1PE18.10)))') &
           (((Real(continuum(IX,IY,IZ)*EXPAUX(IX,IY,IZ)),IX=1,NXHERM),IY=1,NYHERM),IZ=1,NZHERM)

    End Subroutine write_density

    Subroutine close_density(fortran_unit, error)

      Integer(ipr), INTENT(IN) :: fortran_unit
      Integer(ipr), INTENT(OUT) :: error

      error=0
      Close(Unit=fortran_unit,Iostat=error)

    End Subroutine close_density

    !---------------------------------------------------------------------!
    !  This subroutine computes the Skyrme interaction energy between two !
    !  sub-systems characterized by densities '1' and '2'. The density-   !
    !  dependent part of the interaction is computed from the total       !
    !  density of the two subsystems.                                     !
    !---------------------------------------------------------------------!
    Subroutine interaction_energy_Skyrme(NXHERM,NYHERM,NZHERM,ENESKY,ENEVEN,ENEODD, &
                                                ENREAR,ENE_W0,EEVEW0,EODDW0,EKESKY, &
                                                DN1RHO,DN1TAU,DN1LPR,DN1DIV,PN1RHO, &
                                  DN1SPI,DN1KIS,DN1GRR,DN1LPS,DN1ROS,DN1ROC,DN1CUR, &
                                                                     DN1SCU,DN1DES, &
                                                DP1RHO,DP1TAU,DP1LPR,DP1DIV,PP1RHO, &
                                  DP1SPI,DP1KIS,DP1GRR,DP1LPS,DP1ROS,DP1ROC,DP1CUR, &
                                                                     DP1SCU,DP1DES, &
                                                DN2RHO,DN2TAU,DN2LPR,DN2DIV,PN2RHO, &
                                  DN2SPI,DN2KIS,DN2GRR,DN2LPS,DN2ROS,DN2ROC,DN2CUR, &
                                                                     DN2SCU,DN2DES, &
                                                DP2RHO,DP2TAU,DP2LPR,DP2DIV,PP2RHO, &
                                  DP2SPI,DP2KIS,DP2GRR,DP2LPS,DP2ROS,DP2ROC,DP2CUR, &
                                                                     DP2SCU,DP2DES)
      Integer(ipr), INTENT(IN) :: NXHERM,NYHERM,NZHERM
      Real(pr), INTENT(INOUT) :: ENESKY,ENEVEN,ENEODD,ENREAR,ENE_W0,EEVEW0,EODDW0
      Complex(pr), INTENT(INOUT) ::  EKESKY

      Complex(pr), Allocatable, INTENT(INOUT) :: DN1RHO(:,:,:),DN1TAU(:,:,:),DN1LPR(:,:,:),&
                                                 DN1DIV(:,:,:),PN1RHO(:,:,:)
      Complex(pr), Allocatable, INTENT(INOUT) :: DN1SPI(:,:,:,:),DN1KIS(:,:,:,:),&
                                                 DN1GRR(:,:,:,:),DN1LPS(:,:,:,:)
      Complex(pr), Allocatable, INTENT(INOUT) :: DN1ROS(:,:,:,:),DN1ROC(:,:,:,:),DN1CUR(:,:,:,:)
      Complex(pr), Allocatable, INTENT(INOUT) :: DN1SCU(:,:,:,:,:),DN1DES(:,:,:,:,:)

      Complex(pr), Allocatable, INTENT(INOUT) :: DP1RHO(:,:,:),DP1TAU(:,:,:),DP1LPR(:,:,:),&
                                                 DP1DIV(:,:,:),PP1RHO(:,:,:)
      Complex(pr), Allocatable, INTENT(INOUT) :: DP1SPI(:,:,:,:),DP1KIS(:,:,:,:),&
                                                 DP1GRR(:,:,:,:),DP1LPS(:,:,:,:)
      Complex(pr), Allocatable, INTENT(INOUT) :: DP1ROS(:,:,:,:),DP1ROC(:,:,:,:),DP1CUR(:,:,:,:)
      Complex(pr), Allocatable, INTENT(INOUT) :: DP1SCU(:,:,:,:,:),DP1DES(:,:,:,:,:)

      Complex(pr), Allocatable, INTENT(INOUT) :: DN2RHO(:,:,:),DN2TAU(:,:,:),DN2LPR(:,:,:),&
                                                 DN2DIV(:,:,:),PN2RHO(:,:,:)
      Complex(pr), Allocatable, INTENT(INOUT) :: DN2SPI(:,:,:,:),DN2KIS(:,:,:,:),&
                                                 DN2GRR(:,:,:,:),DN2LPS(:,:,:,:)
      Complex(pr), Allocatable, INTENT(INOUT) :: DN2ROS(:,:,:,:),DN2ROC(:,:,:,:),DN2CUR(:,:,:,:)
      Complex(pr), Allocatable, INTENT(INOUT) :: DN2SCU(:,:,:,:,:),DN2DES(:,:,:,:,:)

      Complex(pr), Allocatable, INTENT(INOUT) :: DP2RHO(:,:,:),DP2TAU(:,:,:),DP2LPR(:,:,:),&
                                                 DP2DIV(:,:,:),PP2RHO(:,:,:)
      Complex(pr), Allocatable, INTENT(INOUT) :: DP2SPI(:,:,:,:),DP2KIS(:,:,:,:),&
                                                 DP2GRR(:,:,:,:),DP2LPS(:,:,:,:)
      Complex(pr), Allocatable, INTENT(INOUT) :: DP2ROS(:,:,:,:),DP2ROC(:,:,:,:),DP2CUR(:,:,:,:)
      Complex(pr), Allocatable, INTENT(INOUT) :: DP2SCU(:,:,:,:,:),DP2DES(:,:,:,:,:)

      Integer(ipr) :: NFIPRI
      Real(pr) :: CRHO_T,CRHO_S,CRHODT,CRHODS,CLPR_T,CLPR_S, &
                  CTAU_T,CTAU_S,CSCU_T,CSCU_S,CDIV_T,CDIV_S, &
                  CSPI_T,CSPI_S,CSPIDT,CSPIDS,CLPS_T,CLPS_S, &
                  CCUR_T,CCUR_S,CKIS_T,CKIS_S,CROT_T,CROT_S
      Real(pr) :: T0,X0,T1,X1,T2,X2,T3,X3,W0,POWER
      Real(pr) :: FOURWG,FOURPT
      Complex(pr) :: DENSIC,DENCHC

      Integer(ipr) :: IX,IY,IZ,K,L,NUCOUP

      Real(pr) :: QPPSCR,QMMSCR,QPMSCR,QMPSCR,QPPDIS,QMMDIS,QPMDIS,QMPDIS, &
                  QPPDSR,QMMDSR,QPMDSR,QMPDSR,QPPSCS,QMMSCS,QPMSCS,QMPSCS, &
                  QPPSCJ,QMMSCJ,QPMSCJ,QMPSCJ,QPPCUS,QMMCUS,QPMCUS,QMPCUS, &
                  QPPROS,QMMROS,QPMROS,QMPROS,QPPDES,QMMDES,QPMDES,QMPDES, &
                  W_HERM
      Real(pr) :: QUNDIS_1,QUPDIS_1,QUNDIS_2,QUPDIS_2

      Complex(pr) :: C_ZERO,DENPOW
      Complex(pr) :: QUNRHO,QUPRHO,QUTRHO,QUNLPR,QUPLPR,QUTLPR, &
                     QUNTAU,QUPTAU,QUTTAU,QUNSCU,QUPSCU,QUTSCU, &
                     QUNDIV,QUPDIV,QUTDIV,QUNSPI,QUPSPI,QUTSPI, &
                     QUNLPS,QUPLPS,QUTLPS,QUNCUR,QUPCUR,QUTCUR, &
                     QUNKIS,QUPKIS,QUTKIS,QUNROT,QUPROT,QUTROT
      Complex(pr), Dimension(NDCOUP) :: ZCPALL
      Complex(pr), Allocatable :: DN1ROJ(:,:,:,:),DP1ROJ(:,:,:,:)

      COMMON                                       &
             /DENTOC/ DENSIC(NDXHRM,NDYHRM,NDZHRM),&
                      DENCHC(NDXHRM,NDYHRM,NDZHRM)
      COMMON                                       &
             /CCPSKY/ CRHO_T,CRHO_S,CRHODT,CRHODS, &
                      CLPR_T,CLPR_S,               &
                      CTAU_T,CTAU_S,               &
                      CSCU_T,CSCU_S,               &
                      CDIV_T,CDIV_S,               &
                      CSPI_T,CSPI_S,CSPIDT,CSPIDS, &
                      CLPS_T,CLPS_S,               &
                      CCUR_T,CCUR_S,               &
                      CKIS_T,CKIS_S,               &
                      CROT_T,CROT_S
      COMMON                                           &
             /INTSKY/ T0,X0,T1,X1,T2,X2,T3,X3,W0,POWER
      COMMON                                     &
             /INTSTO/ FOURWG(1:NDGAUS,1:NDKART), &
                      FOURPT(1:NDGAUS,1:NDKART)
      COMMON                                     &
             /CFIPRI/ NFIPRI

      C_ZERO=Cmplx(0.0_pr,0.0_pr); EKESKY=C_ZERO; ZCPALL(:)=C_ZERO
      ENESKY=0.0_pr; ENREAR=0.0_pr

      fis_CCPALL( 1)=CRHO_T; fis_CCPALL( 2)=CRHO_S; fis_CCPALL( 3)=CRHODT; fis_CCPALL( 4)=CRHODS
      fis_CCPALL( 5)=CLPR_T; fis_CCPALL( 6)=CLPR_S; fis_CCPALL( 7)=CTAU_T; fis_CCPALL( 8)=CTAU_S
      fis_CCPALL( 9)=CSCU_T; fis_CCPALL(10)=CSCU_S; fis_CCPALL(11)=CDIV_T; fis_CCPALL(12)=CDIV_S
      fis_CCPALL(13)=CSPI_T; fis_CCPALL(14)=CSPI_S; fis_CCPALL(15)=CSPIDT; fis_CCPALL(16)=CSPIDS
      fis_CCPALL(17)=CLPS_T; fis_CCPALL(18)=CLPS_S; fis_CCPALL(19)=CCUR_T; fis_CCPALL(20)=CCUR_S
      fis_CCPALL(21)=CKIS_T; fis_CCPALL(22)=CKIS_S; fis_CCPALL(23)=CROT_T; fis_CCPALL(24)=CROT_S

      Allocate(DN1ROJ(NXHERM,NYHERM,NZHERM,NDKART),DP1ROJ(NXHERM,NYHERM,NZHERM,NDKART))

      Do IZ=1,NZHERM
         Do IY=1,NYHERM
            Do IX=1,NXHERM

               ! Calculating the vector part of the spin-current tensor

               DN1ROJ(IX,IY,IZ,1)=DN1SCU(IX,IY,IZ,2,3)-DN1SCU(IX,IY,IZ,3,2)
               DN1ROJ(IX,IY,IZ,2)=DN1SCU(IX,IY,IZ,3,1)-DN1SCU(IX,IY,IZ,1,3)
               DN1ROJ(IX,IY,IZ,3)=DN1SCU(IX,IY,IZ,1,2)-DN1SCU(IX,IY,IZ,2,1)

               DP1ROJ(IX,IY,IZ,1)=DP1SCU(IX,IY,IZ,2,3)-DP1SCU(IX,IY,IZ,3,2)
               DP1ROJ(IX,IY,IZ,2)=DP1SCU(IX,IY,IZ,3,1)-DP1SCU(IX,IY,IZ,1,3)
               DP1ROJ(IX,IY,IZ,3)=DP1SCU(IX,IY,IZ,1,2)-DP1SCU(IX,IY,IZ,2,1)

               ! Calculating  the  a u x i l i a r y   c o n t r a c t i o n s

               QUNRHO = DN1RHO(IX,IY,IZ)*DN2RHO(IX,IY,IZ)
               QUPRHO = DP1RHO(IX,IY,IZ)*DP2RHO(IX,IY,IZ)
               QUTRHO =(DN1RHO(IX,IY,IZ)+DP1RHO(IX,IY,IZ))*(DN2RHO(IX,IY,IZ)+DP2RHO(IX,IY,IZ))

               QUNLPR = DN1RHO(IX,IY,IZ)*DN2LPR(IX,IY,IZ)
               QUPLPR = DP1RHO(IX,IY,IZ)*DP2LPR(IX,IY,IZ)
               QUTLPR =(DN1RHO(IX,IY,IZ)+DP1RHO(IX,IY,IZ))*(DN2LPR(IX,IY,IZ)+DP2LPR(IX,IY,IZ))

               QUNTAU = DN1RHO(IX,IY,IZ)*DN2TAU(IX,IY,IZ)
               QUPTAU = DP1RHO(IX,IY,IZ)*DP2TAU(IX,IY,IZ)
               QUTTAU =(DN1RHO(IX,IY,IZ)+DP1RHO(IX,IY,IZ))*(DN2TAU(IX,IY,IZ)+DP2TAU(IX,IY,IZ))

               QUNDIV = DN1RHO(IX,IY,IZ)*DN2DIV(IX,IY,IZ)
               QUPDIV = DP1RHO(IX,IY,IZ)*DP2DIV(IX,IY,IZ)
               QUTDIV =(DN1RHO(IX,IY,IZ)+DP1RHO(IX,IY,IZ))*(DN2DIV(IX,IY,IZ)+DP2DIV(IX,IY,IZ))

               QUNSPI=C_ZERO; QUPSPI=C_ZERO; QUTSPI=C_ZERO
               QUNCUR=C_ZERO; QUPCUR=C_ZERO; QUTCUR=C_ZERO
               QUNLPS=C_ZERO; QUPLPS=C_ZERO; QUTLPS=C_ZERO
               QUNKIS=C_ZERO; QUPKIS=C_ZERO; QUTKIS=C_ZERO
               QUNROT=C_ZERO; QUPROT=C_ZERO; QUTROT=C_ZERO
               QUNSCU=C_ZERO; QUPSCU=C_ZERO; QUTSCU=C_ZERO

               QPPSCR=0.0_pr; QMMSCR=0.0_pr; QPMSCR=0.0_pr; QMPSCR=0.0_pr
               QPPDSR=0.0_pr; QMMDSR=0.0_pr; QPMDSR=0.0_pr; QMPDSR=0.0_pr
               QPPSCS=0.0_pr; QMMSCS=0.0_pr; QPMSCS=0.0_pr; QMPSCS=0.0_pr
               QPPSCJ=0.0_pr; QMMSCJ=0.0_pr; QPMSCJ=0.0_pr; QMPSCJ=0.0_pr
               QPPCUS=0.0_pr; QMMCUS=0.0_pr; QPMCUS=0.0_pr; QMPCUS=0.0_pr
               QPPROS=0.0_pr; QMMROS=0.0_pr; QPMROS=0.0_pr; QMPROS=0.0_pr

               QPPDES=0.0_pr; QMMDES=0.0_pr; QPMDES=0.0_pr; QMPDES=0.0_pr
               QPPDIS=0.0_pr; QMMDIS=0.0_pr; QPMDIS=0.0_pr; QMPDIS=0.0_pr

               QUNDIS_1=0.0_pr; QUPDIS_1=0.0_pr; QUNDIS_2=0.0_pr; QUPDIS_2=0.0_pr

               Do K=1,NDKART

                  QUNSPI=QUNSPI + DN1SPI(IX,IY,IZ,K)*DN2SPI(IX,IY,IZ,K)
                  QUPSPI=QUPSPI + DP1SPI(IX,IY,IZ,K)*DP2SPI(IX,IY,IZ,K)
                  QUTSPI=QUTSPI +(DN1SPI(IX,IY,IZ,K)+DP1SPI(IX,IY,IZ,K))*(DN2SPI(IX,IY,IZ,K)+DP2SPI(IX,IY,IZ,K))

                  QUNCUR=QUNCUR + DN1CUR(IX,IY,IZ,K)*DN2CUR(IX,IY,IZ,K)
                  QUPCUR=QUPCUR + DP1CUR(IX,IY,IZ,K)*DP2CUR(IX,IY,IZ,K)
                  QUTCUR=QUTCUR +(DN1CUR(IX,IY,IZ,K)+DP1CUR(IX,IY,IZ,K))*(DN2CUR(IX,IY,IZ,K)+DP2CUR(IX,IY,IZ,K))

                  QUNLPS=QUNLPS + DN1SPI(IX,IY,IZ,K)*DN2LPS(IX,IY,IZ,K)
                  QUPLPS=QUPLPS + DP1SPI(IX,IY,IZ,K)*DP2LPS(IX,IY,IZ,K)
                  QUTLPS=QUTLPS +(DN1SPI(IX,IY,IZ,K)+DP1SPI(IX,IY,IZ,K))*(DN2LPS(IX,IY,IZ,K)+DP2LPS(IX,IY,IZ,K))

                  QUNKIS=QUNKIS + DN1SPI(IX,IY,IZ,K)*DN2KIS(IX,IY,IZ,K)
                  QUPKIS=QUPKIS + DP1SPI(IX,IY,IZ,K)*DP2KIS(IX,IY,IZ,K)
                  QUTKIS=QUTKIS +(DN1SPI(IX,IY,IZ,K)+DP1SPI(IX,IY,IZ,K))*(DN2KIS(IX,IY,IZ,K)+DP2KIS(IX,IY,IZ,K))

                  QUNROT=QUNROT + DN1SPI(IX,IY,IZ,K)*DN2ROC(IX,IY,IZ,K)
                  QUPROT=QUPROT + DP1SPI(IX,IY,IZ,K)*DP2ROC(IX,IY,IZ,K)
                  QUTROT=QUTROT +(DN1SPI(IX,IY,IZ,K)+DP1SPI(IX,IY,IZ,K))*(DN2ROC(IX,IY,IZ,K)+DP2ROC(IX,IY,IZ,K))

                  ! Symmetry-violating terms

                  QPPSCR=QPPSCR+Real((DN1SCU(IX,IY,IZ,K,K)+DP1SCU(IX,IY,IZ,K,K)) &
                                    *(DN2RHO(IX,IY,IZ    )+DP2RHO(IX,IY,IZ    )))
                  QMMSCR=QMMSCR+Real((DN1SCU(IX,IY,IZ,K,K)-DP1SCU(IX,IY,IZ,K,K)) &
                                    *(DN2RHO(IX,IY,IZ    )-DP2RHO(IX,IY,IZ    )))
                  QPMSCR=QPMSCR+Real((DN1SCU(IX,IY,IZ,K,K)+DP1SCU(IX,IY,IZ,K,K)) &
                                    *(DN2RHO(IX,IY,IZ    )-DP2RHO(IX,IY,IZ    )))
                  QMPSCR=QMPSCR+Real((DN1SCU(IX,IY,IZ,K,K)-DP1SCU(IX,IY,IZ,K,K)) &
                                    *(DN2RHO(IX,IY,IZ    )+DP2RHO(IX,IY,IZ    )))

                  QUNDIS_1=QUNDIS_1+Real(DN1DES(IX,IY,IZ,K,K))
                  QUPDIS_1=QUPDIS_1+Real(DP1DES(IX,IY,IZ,K,K))
                  QUNDIS_2=QUNDIS_2+Real(DN2DES(IX,IY,IZ,K,K))
                  QUPDIS_2=QUPDIS_2+Real(DP2DES(IX,IY,IZ,K,K))

                  QPPDSR=QPPDSR+Real((DN1DES(IX,IY,IZ,K,K)+DP1DES(IX,IY,IZ,K,K)) &
                                    *(DN2RHO(IX,IY,IZ    )+DP2RHO(IX,IY,IZ    )))
                  QMMDSR=QMMDSR+Real((DN1DES(IX,IY,IZ,K,K)-DP1DES(IX,IY,IZ,K,K)) &
                                    *(DN2RHO(IX,IY,IZ    )-DP2RHO(IX,IY,IZ    )))
                  QPMDSR=QPMDSR+Real((DN1DES(IX,IY,IZ,K,K)+DP1DES(IX,IY,IZ,K,K)) &
                                    *(DN2RHO(IX,IY,IZ    )-DP2RHO(IX,IY,IZ    )))
                  QMPDSR=QMPDSR+Real((DN1DES(IX,IY,IZ,K,K)-DP1DES(IX,IY,IZ,K,K)) &
                                    *(DN2RHO(IX,IY,IZ    )+DP2RHO(IX,IY,IZ    )))

                  QPPSCS=QPPSCS + Real((DN1ROJ(IX,IY,IZ,K)+DP1ROJ(IX,IY,IZ,K))*(DN2SPI(IX,IY,IZ,K)+DP2SPI(IX,IY,IZ,K)))
                  QMMSCS=QMMSCS + Real((DN1ROJ(IX,IY,IZ,K)-DP1ROJ(IX,IY,IZ,K))*(DN2SPI(IX,IY,IZ,K)-DP2SPI(IX,IY,IZ,K)))
                  QPMSCS=QPMSCS + Real((DN1ROJ(IX,IY,IZ,K)+DP1ROJ(IX,IY,IZ,K))*(DN2SPI(IX,IY,IZ,K)-DP2SPI(IX,IY,IZ,K)))
                  QMPSCS=QMPSCS + Real((DN1ROJ(IX,IY,IZ,K)-DP1ROJ(IX,IY,IZ,K))*(DN2SPI(IX,IY,IZ,K)+DP2SPI(IX,IY,IZ,K)))

                  QPPSCJ=QPPSCJ + Real((DN1ROJ(IX,IY,IZ,K)+DP1ROJ(IX,IY,IZ,K))*(DN2CUR(IX,IY,IZ,K)+DP2CUR(IX,IY,IZ,K)))
                  QMMSCJ=QMMSCJ + Real((DN1ROJ(IX,IY,IZ,K)-DP1ROJ(IX,IY,IZ,K))*(DN2CUR(IX,IY,IZ,K)-DP2CUR(IX,IY,IZ,K)))
                  QPMSCJ=QPMSCJ + Real((DN1ROJ(IX,IY,IZ,K)+DP1ROJ(IX,IY,IZ,K))*(DN2CUR(IX,IY,IZ,K)-DP2CUR(IX,IY,IZ,K)))
                  QMPSCJ=QMPSCJ + Real((DN1ROJ(IX,IY,IZ,K)-DP1ROJ(IX,IY,IZ,K))*(DN2CUR(IX,IY,IZ,K)+DP2CUR(IX,IY,IZ,K)))

                  QPPCUS=QPPCUS + Real((DN1CUR(IX,IY,IZ,K)+DP1CUR(IX,IY,IZ,K))*(DN2SPI(IX,IY,IZ,K)+DP2SPI(IX,IY,IZ,K)))
                  QMMCUS=QMMCUS + Real((DN1CUR(IX,IY,IZ,K)-DP1CUR(IX,IY,IZ,K))*(DN2SPI(IX,IY,IZ,K)-DP2SPI(IX,IY,IZ,K)))
                  QPMCUS=QPMCUS + Real((DN1CUR(IX,IY,IZ,K)+DP1CUR(IX,IY,IZ,K))*(DN2SPI(IX,IY,IZ,K)-DP2SPI(IX,IY,IZ,K)))
                  QMPCUS=QMPCUS + Real((DN1CUR(IX,IY,IZ,K)-DP1CUR(IX,IY,IZ,K))*(DN2SPI(IX,IY,IZ,K)+DP2SPI(IX,IY,IZ,K)))

                  QPPROS=QPPROS + Real((DN1ROS(IX,IY,IZ,K)+DP1ROS(IX,IY,IZ,K))*(DN2SPI(IX,IY,IZ,K)+DP2SPI(IX,IY,IZ,K)))
                  QMMROS=QMMROS + Real((DN1ROS(IX,IY,IZ,K)-DP1ROS(IX,IY,IZ,K))*(DN2SPI(IX,IY,IZ,K)-DP2SPI(IX,IY,IZ,K)))
                  QPMROS=QPMROS + Real((DN1ROS(IX,IY,IZ,K)+DP1ROS(IX,IY,IZ,K))*(DN2SPI(IX,IY,IZ,K)-DP2SPI(IX,IY,IZ,K)))
                  QMPROS=QMPROS + Real((DN1ROS(IX,IY,IZ,K)-DP1ROS(IX,IY,IZ,K))*(DN2SPI(IX,IY,IZ,K)+DP2SPI(IX,IY,IZ,K)))

                  Do L=1,NDKART

                     QUNSCU=QUNSCU + DN1SCU(IX,IY,IZ,K,L) * DN2SCU(IX,IY,IZ,K,L)
                     QUPSCU=QUPSCU + DP1SCU(IX,IY,IZ,K,L) * DP2SCU(IX,IY,IZ,K,L)
                     QUTSCU=QUTSCU +(DN1SCU(IX,IY,IZ,K,L) + DP1SCU(IX,IY,IZ,K,L)) &
                                   *(DN2SCU(IX,IY,IZ,K,L) + DP2SCU(IX,IY,IZ,K,L))

                    ! Symmetry-violating terms

                     QPPDES=QPPDES + Real((DN1DES(IX,IY,IZ,K,L) + DP1DES(IX,IY,IZ,K,L)) &
                                         *(DN2SCU(IX,IY,IZ,K,L) + DP2SCU(IX,IY,IZ,K,L)))

                     QMMDES=QMMDES + Real((DN1DES(IX,IY,IZ,K,L) - DP1DES(IX,IY,IZ,K,L)) &
                                         *(DN2SCU(IX,IY,IZ,K,L) - DP2SCU(IX,IY,IZ,K,L)))

                     QPMDES=QPMDES + Real((DN1DES(IX,IY,IZ,K,L) + DP1DES(IX,IY,IZ,K,L)) &
                                         *(DN2SCU(IX,IY,IZ,K,L) - DP2SCU(IX,IY,IZ,K,L)))

                     QMPDES=QMPDES + Real((DN1DES(IX,IY,IZ,K,L) - DP1DES(IX,IY,IZ,K,L)) &
                                         *(DN2SCU(IX,IY,IZ,K,L) + DP2SCU(IX,IY,IZ,K,L)))

                  End Do ! End L

               End Do ! End K

               QPPDIS=Real((QUNDIS_1+QUPDIS_1)*(QUNDIS_2+QUPDIS_2))
               QMMDIS=Real((QUNDIS_1-QUPDIS_1)*(QUNDIS_2-QUPDIS_2))
               QPMDIS=Real((QUNDIS_1+QUPDIS_1)*(QUNDIS_2-QUPDIS_2))
               QMPDIS=Real((QUNDIS_1-QUPDIS_1)*(QUNDIS_2+QUPDIS_2))

               W_HERM=FOURWG(IX,1)*FOURWG(IY,2)*FOURWG(IZ,3)

               DENPOW=DENSIC(IX,IY,IZ) ** POWER

               ! Calculating the t0-terms of the skyrme energy density

               ZCPALL(1)  = ZCPALL(1)  + W_HERM *  QUTRHO
               ZCPALL(2)  = ZCPALL(2)  + W_HERM * (QUNRHO + QUPRHO)

               ZCPALL(13) = ZCPALL(13) + W_HERM *  QUTSPI
               ZCPALL(14) = ZCPALL(14) + W_HERM * (QUNSPI + QUPSPI)

               ! Calculating the t1-t2-terms of the skyrme energy density

               ZCPALL(5)  = ZCPALL(5)  + W_HERM *  QUTLPR
               ZCPALL(6)  = ZCPALL(6)  + W_HERM * (QUNLPR + QUPLPR)

               ZCPALL(7)  = ZCPALL(7)  + W_HERM *  QUTTAU
               ZCPALL(8)  = ZCPALL(8)  + W_HERM * (QUNTAU + QUPTAU)

               ZCPALL(19) = ZCPALL(19) + W_HERM *  QUTCUR
               ZCPALL(20) = ZCPALL(20) + W_HERM * (QUNCUR + QUPCUR)

               ZCPALL(17) = ZCPALL(17) + W_HERM *  QUTLPS
               ZCPALL(18) = ZCPALL(18) + W_HERM * (QUNLPS + QUPLPS)

               ZCPALL(21) = ZCPALL(21) + W_HERM *  QUTKIS
               ZCPALL(22) = ZCPALL(22) + W_HERM * (QUNKIS + QUPKIS)

               ZCPALL(9)  = ZCPALL(9)  + W_HERM *  QUTSCU
               ZCPALL(10) = ZCPALL(10) + W_HERM * (QUNSCU + QUPSCU)

               ! Calculating the t3-terms of the skyrme energy density

               ZCPALL(3)  = ZCPALL(3)  + W_HERM * QUTRHO           * DENPOW
               ZCPALL(4)  = ZCPALL(4)  + W_HERM *(QUNRHO + QUPRHO) * DENPOW
               ZCPALL(15) = ZCPALL(15) + W_HERM * QUTSPI           * DENPOW
               ZCPALL(16) = ZCPALL(16) + W_HERM *(QUNSPI + QUPSPI) * DENPOW

               ENREAR = ENREAR + W_HERM * Real( (fis_CCPALL(3 )*QUTRHO + fis_CCPALL(4 )*(QUNRHO + QUPRHO)  &
                                               + fis_CCPALL(15)*QUTSPI + fis_CCPALL(16)*(QUNSPI + QUPSPI)) &
                                               * DENPOW*(POWER/2) )

               ! Calculating the W0-terms of the skyrme energy density

               ZCPALL(11) = ZCPALL(11) + W_HERM *  QUTDIV
               ZCPALL(12) = ZCPALL(12) + W_HERM * (QUNDIV + QUPDIV)

               ZCPALL(23) = ZCPALL(23) + W_HERM *  QUTROT
               ZCPALL(24) = ZCPALL(24) + W_HERM * (QUNROT + QUPROT)

            End Do
         End Do
      End Do

      Deallocate(DN1ROJ,DP1ROJ)

      ! Partial summations of terms in energy
      Do NUCOUP=1,NDCOUP
         fis_DCPALL(NUCOUP)=Real(ZCPALL(NUCOUP))
         fis_ECPALL(NUCOUP)=fis_CCPALL(NUCOUP)*fis_DCPALL(NUCOUP)
         ENESKY            =ENESKY            +fis_ECPALL(NUCOUP)
         EKESKY            =EKESKY            +fis_CCPALL(NUCOUP)*ZCPALL(NUCOUP)
         fis_ECCALL(NUCOUP)=                   fis_CCPALL(NUCOUP)*ZCPALL(NUCOUP)
      End Do

      ENEVEN=0.0_pr; ENEODD=0.0_pr
      Do NUCOUP=1,NDCOUP/2
         ENEVEN=ENEVEN+fis_ECPALL(NUCOUP)
         ENEODD=ENEODD+fis_ECPALL(NUCOUP+NDCOUP/2)
      End Do

      EEVEW0=fis_ECPALL(11)+fis_ECPALL(12)
      EODDW0=fis_ECPALL(23)+fis_ECPALL(24)
      ENE_W0=EEVEW0+EODDW0

    End Subroutine interaction_energy_Skyrme

    !---------------------------------------------------------------------!
    !  This subroutine computes the full interaction energy between two   !
    !  sub-systems '1' and '2'. Based on a partition of the q.p. levels   !
    !  into '1' and '2' sub-systems, it constructs the Skyrme densities   !
    !  of each subsystem in r-space and in configuration space (density   !
    !  matrix), and computes the Coulomb and Skyrme interaction energy    !
    !  between the two sub-systems                                        !
    !---------------------------------------------------------------------!
    Subroutine interaction_energy_total(ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                                        ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                                        IPAHFB,MREVER,I_GOGA,IGOGPA,NXMAXX,NYMAXX, &
                                        NZMAXX,NAMEPN,PRINIT,IDEVAR,ITERUN,ISYMDE, &
                                        INIROT,INIINV,INIKAR,ISAWAV,IKERNE,NUMCOU, &
                                        BOUCOU,IN_FIX,IZ_FIX,IPAIRI,JETACM,IROTAT, &
                                        EFERMN,EFERMP,INDJOB,IDEALL,IDELOC,IDECON, &
                                        IPNMIX,ITIREP,MIN_QP,I_REGA)
      Use SAVQUA
      Use ALLQUA
      Use WAVR_L
      Use FRAGFL
      Use MAD_PP

      Logical,INTENT(IN) :: PRINIT
      Integer, INTENT(IN) :: ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                             ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                             IPAHFB,MREVER,I_GOGA,IGOGPA,NXMAXX,NYMAXX, &
                             NZMAXX,IN_FIX,IZ_FIX,IDEVAR,ITERUN,ISYMDE, &
                             INIROT,INIINV,INIKAR,ISAWAV,IKERNE,NUMCOU, &
                             IPAIRI,JETACM,IROTAT,INDJOB,IDEALL,IDELOC, &
                             IDECON,IPNMIX,ITIREP,MIN_QP,I_REGA
      Real(pr), INTENT(IN) :: BOUCOU,EFERMN,EFERMP

      Character(Len=8) :: NAMEPN
      Character(Len=18) :: type_density
      Character(Len=68) :: filename

      Logical :: COR_CM

      Integer(ipr) :: NOTOCC,NOCISO,NFIPRI,NUMBQP,ICHARG,ITPNMX
      Real(pr) :: PRHO_N,PRHODN,PRHOSN,POWERN,PRHO_P,PRHODP,PRHOSP,POWERP, &
                  EXPAUX,VCOISO,V_CORR,HBMASS,HBMRPA,HBMINP,FOURWG,FOURPT

      Integer(ipr) :: IX,IY,IZ,KFRAGM,ierr,i
      Integer(ipr) :: fortran_unit,error,I_TYPE,ICOTYP,ISEXAC,JCOUDI,JCOUEX
#if(USE_MPI==1)
      Integer(ipr) :: mpi_rank,mpi_err
#endif

      Real(pr) :: GPAIRN,DELTAN,EPAI_N,EREA_N,GPAIRP,DELTAP,EPAI_P,EREA_P, &
                  ENESKY,ENEVEN,ENEODD,ENREAR,ENE_W0,EEVEW0,EODDW0,ESUMN1, &
                  ESUMN2,ESUMP1,ESUMP2,W_HERM,SLOWEV

      Real(pr), Dimension(0:NDKART) :: DLINSN,ELINSN,PLINSN,TLINSN,ALINLN,PLINLN,PLINKN, &
                                       DLINSP,ELINSP,PLINSP,TLINSP,ALINLP,PLINLP,PLINKP, &
                                       ALINL1,TLINT1,ALINL2,TLINT2

      Complex(pr), Dimension(0:NDKART) :: DKINSQ,EKINSQ,PKINSQ,TKINSQ,AKINLI,PKINLI,PKINKI

      Complex(pr) :: EKESKY,EKECOD,EKECOE,EKESCA,EKEVEC,C_ZERO

      Complex(pr), Allocatable :: density_pp(:,:,:,:)

      Complex(pr), Allocatable :: WN_CEN(:,:,:),WP_CEN(:,:,:)

      Complex(pr), Allocatable :: DN1RHO(:,:,:),DN1TAU(:,:,:),DN1LPR(:,:,:),&
                                  DN1DIV(:,:,:),PN1RHO(:,:,:)
      Complex(pr), Allocatable :: DN1SPI(:,:,:,:),DN1KIS(:,:,:,:),DN1GRR(:,:,:,:),DN1LPS(:,:,:,:)
      Complex(pr), Allocatable :: DN1ROS(:,:,:,:),DN1ROC(:,:,:,:),DN1CUR(:,:,:,:)
      Complex(pr), Allocatable :: DN1SCU(:,:,:,:,:),DN1DES(:,:,:,:,:)

      Complex(pr), Allocatable :: DP1RHO(:,:,:),DP1TAU(:,:,:),DP1LPR(:,:,:),&
                                  DP1DIV(:,:,:),PP1RHO(:,:,:)
      Complex(pr), Allocatable :: DP1SPI(:,:,:,:),DP1KIS(:,:,:,:),DP1GRR(:,:,:,:),DP1LPS(:,:,:,:)
      Complex(pr), Allocatable :: DP1ROS(:,:,:,:),DP1ROC(:,:,:,:),DP1CUR(:,:,:,:)
      Complex(pr), Allocatable :: DP1SCU(:,:,:,:,:),DP1DES(:,:,:,:,:)

      Complex(pr), Allocatable :: DN2RHO(:,:,:),DN2TAU(:,:,:),DN2LPR(:,:,:),&
                                  DN2DIV(:,:,:),PN2RHO(:,:,:)
      Complex(pr), Allocatable :: DN2SPI(:,:,:,:),DN2KIS(:,:,:,:),DN2GRR(:,:,:,:),DN2LPS(:,:,:,:)
      Complex(pr), Allocatable :: DN2ROS(:,:,:,:),DN2ROC(:,:,:,:),DN2CUR(:,:,:,:)
      Complex(pr), Allocatable :: DN2SCU(:,:,:,:,:),DN2DES(:,:,:,:,:)

      Complex(pr), Allocatable :: DP2RHO(:,:,:),DP2TAU(:,:,:),DP2LPR(:,:,:),&
                                  DP2DIV(:,:,:),PP2RHO(:,:,:)
      Complex(pr), Allocatable :: DP2SPI(:,:,:,:),DP2KIS(:,:,:,:),DP2GRR(:,:,:,:),DP2LPS(:,:,:,:)
      Complex(pr), Allocatable :: DP2ROS(:,:,:,:),DP2ROC(:,:,:,:),DP2CUR(:,:,:,:)
      Complex(pr), Allocatable :: DP2SCU(:,:,:,:,:),DP2DES(:,:,:,:,:)

      COMMON                                     &
             /INTSTO/ FOURWG(1:NDGAUS,1:NDKART), &
                      FOURPT(1:NDGAUS,1:NDKART)
      COMMON                               &
             /PLANCK/ HBMASS,HBMRPA,HBMINP
      COMMON                                       &
             /CCPPAI/ PRHO_N,PRHODN,PRHOSN,POWERN, &
                      PRHO_P,PRHODP,PRHOSP,POWERP
      COMMON                                      &
             /DENEXP/ EXPAUX(NDXHRM,NDYHRM,NDZHRM)
      COMMON                                             &
             /LIMOCC/ NOTOCC(1:NDSTAT,0:NDREVE)          &
             /OCCISO/ NOCISO(1:NDSTAT,0:NDREVE,0:NDISOS)
      COMMON                                             &
             /CORPAR/ V_CORR(1:NDSTAT,0:NDREVE)          &
             /CORISO/ VCOISO(1:NDSTAT,0:NDREVE,0:NDISOS)
      COMMON                                    &
             /QPNUMS/ NUMBQP(0:NDREVE,0:NDISOS)
      COMMON                                             &
             /CFIPRI/ NFIPRI

      C_ZERO=Cmplx(0.0_pr,0.0_pr)

      ! Reinitialization of variables saved as public variables in module hfodd_fission_variables
      ESKINT=0.0_pr; EP_INT=0.0_pr; EC_INT=0.0_pr; EC1INT=0.0_pr; EC2INT=0.0_pr; ECDINT_ex=0.0_pr; ECXINT_ex=0.0_pr
      ESKFR1=0.0_pr; EPAIR1=0.0_pr; ECOUL1=0.0_pr; ECOEX1=0.0_pr; EKIN1N=0.0_pr; EKIN1P=0.0_pr
      ESKFR2=0.0_pr; EPAIR2=0.0_pr; ECOUL2=0.0_pr; ECOEX2=0.0_pr; EKIN2N=0.0_pr; EKIN2P=0.0_pr
      TPSQA1(:)=0.0_pr; TPSQA2(:)=0.0_pr

      ! Initialization of internal variables
      ENESKY=0.0_pr; ENEVEN=0.0_pr; ENEODD=0.0_pr; ENREAR=0.0_pr; ENE_W0=0.0_pr; EEVEW0=0.0_pr
      ESUMN1=0.0_pr; ESUMN2=0.0_pr; ESUMP1=0.0_pr; ESUMP2=0.0_pr
      EKESKY=Cmplx(0.0_pr,0.0_pr); EKECOD=Cmplx(0.0_pr,0.0_pr); EKECOE=Cmplx(0.0_pr,0.0_pr)

      If(compute_all) Then

         Write(filename,'("filtered_densities_",i8.8,".dat")') INDJOB
         fortran_unit=41
         error=0; Call open_density(filename, fortran_unit, error)

         !-----------------------------------------!
         !  NEUTRONS - LEFT AND RIGHT FRAGMENTS    !
         !-----------------------------------------!

         ICHARG=0; ITPNMX=0

         BWAQUA(:,:,:)=BSVQUA(:,:,:,ICHARG)
         AWAQUA(:,:,:)=ASVQUA(:,:,:,ICHARG)
         V_CORR(:,:)=1.0_pr
         NOTOCC(:,:)=0

         ! Left fragment
         ierr=0
         Allocate(DN1RHO(NDXHRM,NDYHRM,NDZHRM),DN1TAU(NDXHRM,NDYHRM,NDZHRM),&
                  DN1LPR(NDXHRM,NDYHRM,NDZHRM),DN1DIV(NDXHRM,NDYHRM,NDZHRM),&
                  PN1RHO(NDXHRM,NDYHRM,NDZHRM),STAT=ierr)
         If(ierr.Ne.0) Write(NFIPRI,'("Error 01a")')
         ierr=0
         Allocate(DN1SPI(NDXHRM,NDYHRM,NDZHRM,NDKART),DN1KIS(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                  DN1GRR(NDXHRM,NDYHRM,NDZHRM,NDKART),DN1LPS(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                  DN1ROS(NDXHRM,NDYHRM,NDZHRM,NDKART),DN1ROC(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                  DN1CUR(NDXHRM,NDYHRM,NDZHRM,NDKART),STAT=ierr)
         If(ierr.Ne.0) Write(NFIPRI,'("Error 01b")')
         ierr=0
         Allocate(DN1SCU(NDXHRM,NDYHRM,NDZHRM,NDKART,NDKART),DN1DES(NDXHRM,NDYHRM,NDZHRM,NDKART,NDKART),STAT=ierr)
         If(ierr.Ne.0) Write(NFIPRI,'("Error 01c")')

         If(IDEALL.Eq.1) Then
            I_TYPE=0; type_density="n_all_spectrum.dat"
            Call compute_density(DN1RHO,DN1TAU,DN1LPR,DN1DIV,DN1SPI,DN1KIS, &
                                 DN1GRR,DN1LPS,DN1ROS,DN1ROC,DN1CUR,DN1SCU, &
                                 DN1DES,PN1RHO,                             &
                                 ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                                 ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                                 IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                                 IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                                 ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMN, &
                                 ESUMN1,IPNMIX,ITIREP,ITPNMX,MIN_QP,I_REGA, &
                                 fortran_unit,type_density)
             If(debug.Ge.3) Write(6,'("Writing density n_all_spectrum")')
         End If

         If(IDELOC.Eq.1) Then

            I_TYPE=3; type_density="n_all_discrete.dat"
            Call compute_density(DN1RHO,DN1TAU,DN1LPR,DN1DIV,DN1SPI,DN1KIS, &
                                 DN1GRR,DN1LPS,DN1ROS,DN1ROC,DN1CUR,DN1SCU, &
                                 DN1DES,PN1RHO,                             &
                                 ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                                 ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                                 IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                                 IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                                 ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMN, &
                                 ESUMN1,IPNMIX,ITIREP,ITPNMX,MIN_QP,I_REGA, &
                                 fortran_unit,type_density)
             If(debug.Ge.3) Write(6,'("Writing density n_all_discrete")')

            I_TYPE=5; type_density="n_lef_discrete.dat"
            Call compute_density(DN1RHO,DN1TAU,DN1LPR,DN1DIV,DN1SPI,DN1KIS, &
                                 DN1GRR,DN1LPS,DN1ROS,DN1ROC,DN1CUR,DN1SCU, &
                                 DN1DES,PN1RHO,                             &
                                 ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                                 ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                                 IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                                 IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                                 ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMN, &
                                 ESUMN1,IPNMIX,ITIREP,ITPNMX,MIN_QP,I_REGA, &
                                 fortran_unit,type_density)
             If(debug.Ge.3) Write(6,'("Writing density n_lef_discrete")')

         End If

         If(IDECON.Eq.1) Then

            I_TYPE=4; type_density="n_all_continum.dat"
            Call compute_density(DN1RHO,DN1TAU,DN1LPR,DN1DIV,DN1SPI,DN1KIS, &
                                 DN1GRR,DN1LPS,DN1ROS,DN1ROC,DN1CUR,DN1SCU, &
                                 DN1DES,PN1RHO,                             &
                                 ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                                 ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                                 IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                                 IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                                 ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMN, &
                                 ESUMN1,IPNMIX,ITIREP,ITPNMX,MIN_QP,I_REGA, &
                                 fortran_unit,type_density)
             If(debug.Ge.3) Write(6,'("Writing density n_all_continum")')

            I_TYPE=6; type_density="n_lef_continum.dat"
            Call compute_density(DN1RHO,DN1TAU,DN1LPR,DN1DIV,DN1SPI,DN1KIS, &
                                 DN1GRR,DN1LPS,DN1ROS,DN1ROC,DN1CUR,DN1SCU, &
                                 DN1DES,PN1RHO,                             &
                                 ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                                 ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                                 IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                                 IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                                 ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMN, &
                                 ESUMN1,IPNMIX,ITIREP,ITPNMX,MIN_QP,I_REGA, &
                                 fortran_unit,type_density)
             If(debug.Ge.3) Write(6,'("Writing density n_lef_continum")')

         End If

         I_TYPE=1; type_density="n_lef_spectrum.dat"
         Call compute_density(DN1RHO,DN1TAU,DN1LPR,DN1DIV,DN1SPI,DN1KIS, &
                              DN1GRR,DN1LPS,DN1ROS,DN1ROC,DN1CUR,DN1SCU, &
                              DN1DES,PN1RHO,                             &
                              ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                              ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                              IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                              IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                              ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMN, &
                              ESUMN1,IPNMIX,ITIREP,ITPNMX,MIN_QP,I_REGA, &
                              fortran_unit,type_density)
         If(debug.Ge.3) Write(6,'("Writing density n_lef_spectrum")')

         ! Calculate neutron density matrix for fragment '1'
         Call DENMAC(MREVER,ICHARG,ISIMPY,IPAHFB,WARIGH,WARIGH)

         ! Calculate center of mass correction for neutrons (fragment '1')
         COR_CM = .True.

         Call LINAVR(NXMAXX,NYMAXX,NZMAXX,ISIMPY,ISIGNY,ISIQTY,        &
                     IPAIRI,ICHARG,JETACM,IROTAT,MREVER,COR_CM,        &
                     DLINSN,ELINSN,PLINSN,TLINSN,ALINLN,PLINLN,PLINKN, &
                     DKINSQ,EKINSQ,PKINSQ,TKINSQ,AKINLI,PKINLI,PKINKI)

         ! At this time, the ALINT1, TLINT1 arrays store the c.o.m correction for neutrons in fragment '1'
         ALINL1(:)=ALINLN(:)
         TLINT1(:)=TLINSN(:)

         ! Right fragment
         ierr=0
         Allocate(DN2RHO(NDXHRM,NDYHRM,NDZHRM),DN2TAU(NDXHRM,NDYHRM,NDZHRM),&
                  DN2LPR(NDXHRM,NDYHRM,NDZHRM),DN2DIV(NDXHRM,NDYHRM,NDZHRM),&
                  PN2RHO(NDXHRM,NDYHRM,NDZHRM),STAT=ierr)
         If(ierr.Ne.0) Write(NFIPRI,'("Error 02a")')
         ierr=0
         Allocate(DN2SPI(NDXHRM,NDYHRM,NDZHRM,NDKART),DN2KIS(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                  DN2GRR(NDXHRM,NDYHRM,NDZHRM,NDKART),DN2LPS(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                  DN2ROS(NDXHRM,NDYHRM,NDZHRM,NDKART),DN2ROC(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                  DN2CUR(NDXHRM,NDYHRM,NDZHRM,NDKART),STAT=ierr)
         If(ierr.Ne.0) Write(NFIPRI,'("Error 02b")')
         ierr=0
         Allocate(DN2SCU(NDXHRM,NDYHRM,NDZHRM,NDKART,NDKART),DN2DES(NDXHRM,NDYHRM,NDZHRM,NDKART,NDKART),STAT=ierr)
         If(ierr.Ne.0) Write(NFIPRI,'("Error 02c")')

         If(IDELOC.Eq.1) Then
            I_TYPE=7; type_density="n_rig_discrete.dat"
            Call compute_density(DN2RHO,DN2TAU,DN2LPR,DN2DIV,DN2SPI,DN2KIS, &
                                 DN2GRR,DN2LPS,DN2ROS,DN2ROC,DN2CUR,DN2SCU, &
                                 DN2DES,PN2RHO,                             &
                                 ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                                 ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                                 IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                                 IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                                 ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMN, &
                                 ESUMN2,IPNMIX,ITIREP,ITPNMX,MIN_QP,I_REGA, &
                                 fortran_unit,type_density)
            If(debug.Ge.3) Write(6,'("Writing density n_rig_discrete")')
         End If

         If(IDECON.Eq.1) Then
            I_TYPE=8; type_density="n_rig_continum.dat"
            Call compute_density(DN2RHO,DN2TAU,DN2LPR,DN2DIV,DN2SPI,DN2KIS, &
                                 DN2GRR,DN2LPS,DN2ROS,DN2ROC,DN2CUR,DN2SCU, &
                                 DN2DES,PN2RHO,                             &
                                 ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                                 ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                                 IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                                 IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                                 ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMN, &
                                 ESUMN2,IPNMIX,ITIREP,ITPNMX,MIN_QP,I_REGA, &
                                 fortran_unit,type_density)
            If(debug.Ge.3) Write(6,'("Writing density n_rig_continum")')
         End If

         I_TYPE=2; type_density="n_rig_spectrum.dat"
         Call compute_density(DN2RHO,DN2TAU,DN2LPR,DN2DIV,DN2SPI,DN2KIS, &
                              DN2GRR,DN2LPS,DN2ROS,DN2ROC,DN2CUR,DN2SCU, &
                              DN2DES,PN2RHO,                             &
                              ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                              ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                              IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                              IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                              ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMN, &
                              ESUMN2,IPNMIX,ITIREP,ITPNMX,MIN_QP,I_REGA, &
                              fortran_unit,type_density)
         If(debug.Ge.3) Write(6,'("Writing density n_rig_spectrum")')

         ! Calculate neutron density matrix for fragment '1'
         Call DENMAC(MREVER,ICHARG,ISIMPY,IPAHFB,WARIGH,WARIGH)

         ! Calculate center of mass correction for neutrons (fragment '2')
         COR_CM = .True.

         Call LINAVR(NXMAXX,NYMAXX,NZMAXX,ISIMPY,ISIGNY,ISIQTY,        &
                     IPAIRI,ICHARG,JETACM,IROTAT,MREVER,COR_CM,        &
                     DLINSN,ELINSN,PLINSN,TLINSN,ALINLN,PLINLN,PLINKN, &
                     DKINSQ,EKINSQ,PKINSQ,TKINSQ,AKINLI,PKINLI,PKINKI)

         ! At this time, the ALINT2, TLINT2 arrays store the c.o.m correction for neutrons in fragment '2'
         ALINL2(:)=ALINLN(:)
         TLINT2(:)=TLINSN(:)

         !-----------------------------------------!
         !  PROTONS - LEFT AND RIGHT FRAGMENTS     !
         !-----------------------------------------!

         ICHARG=1; ITPNMX=1

         BWAQUA(:,:,:)=BSVQUA(:,:,:,ICHARG)
         AWAQUA(:,:,:)=ASVQUA(:,:,:,ICHARG)
         V_CORR(:,:)=1.0_pr
         NOTOCC(:,:)=0

         ! Left fragment
         ierr=0
         Allocate(DP1RHO(NDXHRM,NDYHRM,NDZHRM),DP1TAU(NDXHRM,NDYHRM,NDZHRM),&
                  DP1LPR(NDXHRM,NDYHRM,NDZHRM),DP1DIV(NDXHRM,NDYHRM,NDZHRM),&
                  PP1RHO(NDXHRM,NDYHRM,NDZHRM),STAT=ierr)
         If(ierr.Ne.0) Write(NFIPRI,'("Error 03a")')
         ierr=0
         Allocate(DP1SPI(NDXHRM,NDYHRM,NDZHRM,NDKART),DP1KIS(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                  DP1GRR(NDXHRM,NDYHRM,NDZHRM,NDKART),DP1LPS(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                  DP1ROS(NDXHRM,NDYHRM,NDZHRM,NDKART),DP1ROC(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                  DP1CUR(NDXHRM,NDYHRM,NDZHRM,NDKART),STAT=ierr)
         If(ierr.Ne.0) Write(NFIPRI,'("Error 03b")')
         ierr=0
         Allocate(DP1SCU(NDXHRM,NDYHRM,NDZHRM,NDKART,NDKART),DP1DES(NDXHRM,NDYHRM,NDZHRM,NDKART,NDKART),STAT=ierr)
         If(ierr.Ne.0) Write(NFIPRI,'("Error 03c")')

         If(IDEALL.Eq.1) Then
            I_TYPE=0; type_density="p_all_spectrum.dat"
            Call compute_density(DP1RHO,DP1TAU,DP1LPR,DP1DIV,DP1SPI,DP1KIS, &
                                 DP1GRR,DP1LPS,DP1ROS,DP1ROC,DP1CUR,DP1SCU, &
                                 DP1DES,PP1RHO,                             &
                                 ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                                 ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                                 IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                                 IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                                 ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMP, &
                                 ESUMP1,IPNMIX,ITIREP,ITPNMX,MIN_QP,I_REGA, &
                                 fortran_unit,type_density)
            If(debug.Ge.3) Write(6,'("Writing density p_all_spectrum")')
         End If

         If(IDELOC.Eq.1) Then

            I_TYPE=3; type_density="p_all_discrete.dat"
            Call compute_density(DP1RHO,DP1TAU,DP1LPR,DP1DIV,DP1SPI,DP1KIS, &
                                 DP1GRR,DP1LPS,DP1ROS,DP1ROC,DP1CUR,DP1SCU, &
                                 DP1DES,PP1RHO,                             &
                                 ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                                 ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                                 IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                                 IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                                 ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMP, &
                                 ESUMP1,IPNMIX,ITIREP,ITPNMX,MIN_QP,I_REGA, &
                                 fortran_unit,type_density)
            If(debug.Ge.3) Write(6,'("Writing density p_all_discrete")')

            I_TYPE=5; type_density="p_lef_discrete.dat"
            Call compute_density(DP1RHO,DP1TAU,DP1LPR,DP1DIV,DP1SPI,DP1KIS, &
                                 DP1GRR,DP1LPS,DP1ROS,DP1ROC,DP1CUR,DP1SCU, &
                                 DP1DES,PP1RHO,                             &
                                 ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                                 ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                                 IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                                 IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                                 ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMP, &
                                 ESUMP1,IPNMIX,ITIREP,ITPNMX,MIN_QP,I_REGA, &
                                 fortran_unit,type_density)
            If(debug.Ge.3) Write(6,'("Writing density p_lef_discrete")')

         End If

         If(IDECON.Eq.1) Then

            I_TYPE=4; type_density="p_all_continum.dat"
            Call compute_density(DP1RHO,DP1TAU,DP1LPR,DP1DIV,DP1SPI,DP1KIS, &
                                 DP1GRR,DP1LPS,DP1ROS,DP1ROC,DP1CUR,DP1SCU, &
                                 DP1DES,PP1RHO,                             &
                                 ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                                 ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                                 IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                                 IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                                 ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMP, &
                                 ESUMP1,IPNMIX,ITIREP,ITPNMX,MIN_QP,I_REGA, &
                                 fortran_unit,type_density)
            If(debug.Ge.3) Write(6,'("Writing density p_all_continum")')

            I_TYPE=6; type_density="p_lef_continum.dat"
            Call compute_density(DP1RHO,DP1TAU,DP1LPR,DP1DIV,DP1SPI,DP1KIS, &
                                 DP1GRR,DP1LPS,DP1ROS,DP1ROC,DP1CUR,DP1SCU, &
                                 DP1DES,PP1RHO,                             &
                                 ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                                 ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                                 IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                                 IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                                 ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMP, &
                                 ESUMP1,IPNMIX,ITIREP,ITPNMX,MIN_QP,I_REGA, &
                                 fortran_unit,type_density)
            If(debug.Ge.3) Write(6,'("Writing density p_lef_continum")')

         End If

         I_TYPE=1; type_density="p_lef_spectrum.dat"
         Call compute_density(DP1RHO,DP1TAU,DP1LPR,DP1DIV,DP1SPI,DP1KIS, &
                              DP1GRR,DP1LPS,DP1ROS,DP1ROC,DP1CUR,DP1SCU, &
                              DP1DES,PP1RHO,                             &
                              ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                              ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                              IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                              IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                              ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMP, &
                              ESUMP1,IPNMIX,ITIREP,ITPNMX,MIN_QP,I_REGA, &
                              fortran_unit,type_density)
         If(debug.Ge.3) Write(6,'("Writing density p_lef_spectrum")')

         ! Calculate proton density matrix for fragment '1'
         Call DENMAC(MREVER,ICHARG,ISIMPY,IPAHFB,WARIGH,WARIGH)

         ! Saving the proton density of fragment '1' for calculation of exact Coulomb exchange
         If(compute_CoulExc) Then
            ! Backup existing density
            Allocate(density_pp(1:NDBASE,1:NDBASE,0:NDREVE,0:NDISOS))
            density_pp(1:NDBASE,1:NDBASE,0:NDREVE,0:NDISOS) = DEN_PP(1:NDBASE,1:NDBASE,0:NDREVE,0:NDISOS)
            ! Update array DEN_PP with proton density of fragment '1'
            Call SAVDEN(ISIMPY,ICHARG)
         End If

         ! Calculate center of mass correction for protons (fragment '1')
         COR_CM = .True.

         Call LINAVR(NXMAXX,NYMAXX,NZMAXX,ISIMPY,ISIGNY,ISIQTY,        &
                     IPAIRI,ICHARG,JETACM,IROTAT,MREVER,COR_CM,        &
                     DLINSP,ELINSP,PLINSP,TLINSP,ALINLP,PLINLP,PLINKP, &
                     DKINSQ,EKINSQ,PKINSQ,TKINSQ,AKINLI,PKINLI,PKINKI)

         Do i=0,NDKART
            TPSQA1(i)=TLINT1(i)+TLINSP(i)+2.0_pr*ALINL1(i)*ALINLP(i)
         End Do

         !-----------------------------------------!
         !     Coulomb energies for fragment 1     !
         !-----------------------------------------!

         ! Source `DP1RHO' in Coulomb matrix elements is '1', density matrix also comes from '1'
         KFRAGM=1
         Call interaction_energy_Coulomb(NXHERM,NYHERM,NZHERM,NXMAXX,NYMAXX,NZMAXX, &
                                         NUMCOU,BOUCOU,ISIMPY,IKERNE,KFRAGM,DP1RHO, &
                                                                            EKECOD)
         ECOUL1 = Real(EKECOD)

         Call exchange_Coulomb(NXHERM,NYHERM,NZHERM,DP1RHO,EKECOE)
         ECOEX1 = Real(EKECOE)

         ! Right fragment
         ierr=0
         Allocate(DP2RHO(NDXHRM,NDYHRM,NDZHRM),DP2TAU(NDXHRM,NDYHRM,NDZHRM),&
                  DP2LPR(NDXHRM,NDYHRM,NDZHRM),DP2DIV(NDXHRM,NDYHRM,NDZHRM),&
                  PP2RHO(NDXHRM,NDYHRM,NDZHRM),STAT=ierr)
         If(ierr.Ne.0) Write(NFIPRI,'("Error 04a")')
         ierr=0
         Allocate(DP2SPI(NDXHRM,NDYHRM,NDZHRM,NDKART),DP2KIS(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                  DP2GRR(NDXHRM,NDYHRM,NDZHRM,NDKART),DP2LPS(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                  DP2ROS(NDXHRM,NDYHRM,NDZHRM,NDKART),DP2ROC(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                  DP2CUR(NDXHRM,NDYHRM,NDZHRM,NDKART),STAT=ierr)
         If(ierr.Ne.0) Write(NFIPRI,'("Error 04b")')
         ierr=0
         Allocate(DP2SCU(NDXHRM,NDYHRM,NDZHRM,NDKART,NDKART),DP2DES(NDXHRM,NDYHRM,NDZHRM,NDKART,NDKART),STAT=ierr)
         If(ierr.Ne.0) Write(NFIPRI,'("Error 04c")')

         If(IDELOC.Eq.1) Then

            I_TYPE=7; type_density="p_rig_discrete.dat"
            Call compute_density(DP2RHO,DP2TAU,DP2LPR,DP2DIV,DP2SPI,DP2KIS, &
                                 DP2GRR,DP2LPS,DP2ROS,DP2ROC,DP2CUR,DP2SCU, &
                                 DP2DES,PP2RHO,                             &
                                 ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                                 ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                                 IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                                 IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                                 ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMP, &
                                 ESUMP2,IPNMIX,ITIREP,ITPNMX,MIN_QP,I_REGA, &
                                 fortran_unit,type_density)
            If(debug.Ge.3) Write(6,'("Writing density p_rig_discrete")')

         End If

         If(IDECON.Eq.1) Then

            I_TYPE=8; type_density="p_rig_continum.dat"
            Call compute_density(DP2RHO,DP2TAU,DP2LPR,DP2DIV,DP2SPI,DP2KIS, &
                                 DP2GRR,DP2LPS,DP2ROS,DP2ROC,DP2CUR,DP2SCU, &
                                 DP2DES,PP2RHO,                             &
                                 ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                                 ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                                 IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                                 IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                                 ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMP, &
                                 ESUMP2,IPNMIX,ITIREP,ITPNMX,MIN_QP,I_REGA, &
                                 fortran_unit,type_density)
            If(debug.Ge.3) Write(6,'("Writing density p_rig_continum")')

         End If

         I_TYPE=2; type_density="p_rig_spectrum.dat"
         Call compute_density(DP2RHO,DP2TAU,DP2LPR,DP2DIV,DP2SPI,DP2KIS, &
                              DP2GRR,DP2LPS,DP2ROS,DP2ROC,DP2CUR,DP2SCU, &
                              DP2DES,PP2RHO,                             &
                              ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                              ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                              IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                              IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                              ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMP, &
                              ESUMP2,IPNMIX,ITIREP,ITPNMX,MIN_QP,I_REGA, &
                              fortran_unit,type_density)
         If(debug.Ge.3) Write(6,'("Writing density p_rig_spectrum")')

         !-----------------------------------------!
         !   Closing the filtered densities file   !
         !-----------------------------------------!

         error=0; Call close_density(fortran_unit, error)

         !-----------------------------------------!
         !     Coulomb energies for fragment 2     !
         !-----------------------------------------!

         If(debug.Ge.2) Then
            Write(6,'("Testing normalization of neutrons")')
            Call TSTDEN(NXHERM,NYHERM,NZHERM,DN1RHO,DN2RHO)
            Write(6,'("Testing normalization of protons")')
            Call TSTDEN(NXHERM,NYHERM,NZHERM,DP1RHO,DP2RHO)
         End If

         ! Source `DP2RHO' in Coulomb matrix elements is '2', density matrix still comes from '1'
         KFRAGM=2
         Call interaction_energy_Coulomb(NXHERM,NYHERM,NZHERM,NXMAXX,NYMAXX,NZMAXX, &
                                         NUMCOU,BOUCOU,ISIMPY,IKERNE,KFRAGM,DP2RHO, &
                                                                            EKECOD)
         EC1INT = Real(EKECOD)

         ! Calculate proton density matrix for fragment '2'
         Call DENMAC(MREVER,ICHARG,ISIMPY,IPAHFB,WARIGH,WARIGH)

         If(compute_CoulExc) Then

            ! Calculate exact exchange contribution to Coulomb interaction energy
            ICOTYP=7; ISEXAC=1
            JCOUDI=2; JCOUEX=2
            SLOWEV=0.0_pr
            Call COUENE(NXHERM,NYHERM,NZHERM,NXMAXX,NYMAXX,NZMAXX, &
                                                           ISIMPY, &
                                             ICOTYP,JCOUDI,JCOUEX, &
                                             EKECOD,EKESCA,EKEVEC, &
                                                    SLOWEV,ISEXAC)
            ECDINT_ex = 2.0_pr* Real(EKECOD)
            ECXINT_ex = 2.0_pr*(Real(EKESCA)+Real(EKEVEC))

            ! Restore total density
            DEN_PP(1:NDBASE,1:NDBASE,0:NDREVE,0:NDISOS) = density_pp(1:NDBASE,1:NDBASE,0:NDREVE,0:NDISOS)
            Deallocate(density_pp)

         End If

         ! Calculate center of mass correction for protons (fragment '2')
         COR_CM = .True.

         Call LINAVR(NXMAXX,NYMAXX,NZMAXX,ISIMPY,ISIGNY,ISIQTY,        &
                     IPAIRI,ICHARG,JETACM,IROTAT,MREVER,COR_CM,        &
                     DLINSP,ELINSP,PLINSP,TLINSP,ALINLP,PLINLP,PLINKP, &
                     DKINSQ,EKINSQ,PKINSQ,TKINSQ,AKINLI,PKINLI,PKINKI)

         Do i=0,NDKART
            TPSQA2(i)=TLINT2(i)+TLINSP(i)+2.0_pr*ALINL2(i)*ALINLP(i)
         End Do

         ! Source `DP2RHO' in Coulomb matrix elements is '2', density matrix now comes from '2'
         KFRAGM=2
         Call interaction_energy_Coulomb(NXHERM,NYHERM,NZHERM,NXMAXX,NYMAXX,NZMAXX, &
                                         NUMCOU,BOUCOU,ISIMPY,IKERNE,KFRAGM,DP2RHO, &
                                                                            EKECOD)
         ECOUL2 = Real(EKECOD)

         Call exchange_Coulomb(NXHERM,NYHERM,NZHERM,DP2RHO,EKECOE)
         ECOEX2 = Real(EKECOE)

         ! Source `DP1RHO' in Coulomb matrix elements is back to '1', density matrix still comes from '2'
         KFRAGM=1
         Call interaction_energy_Coulomb(NXHERM,NYHERM,NZHERM,NXMAXX,NYMAXX,NZMAXX, &
                                         NUMCOU,BOUCOU,ISIMPY,IKERNE,KFRAGM,DP1RHO, &
                                                                            EKECOD)
         EC2INT = Real(EKECOD)

         ! Total kinetic energy available is the total interaction Coulomb energy
         EC_INT = EC1INT+EC2INT

         !---------------------------------------------------------------------!
         !                           PAIRING ENERGIES                          !
         !---------------------------------------------------------------------!

         If(compute_pairing) Then

            Allocate(WN_CEN(NDXHRM,NDYHRM,NDZHRM)); WN_CEN(:,:,:)=C_ZERO
            Allocate(WP_CEN(NDXHRM,NDYHRM,NDZHRM)); WP_CEN(:,:,:)=C_ZERO

            ! Compute pairing energy for fragment 1
            Call pairing_field(NXHERM,NYHERM,NZHERM,PN1RHO,PP1RHO,WN_CEN,WP_CEN)

            Call interaction_energy_pairing(NXHERM,NYHERM,NZHERM,PRHO_N,PRHODN,POWERN, &
                                            DN1RHO,PN1RHO,WN_CEN,GPAIRN,DELTAN,EPAI_N, &
                                                                        EREA_N,IN_FIX)

            Call interaction_energy_pairing(NXHERM,NYHERM,NZHERM,PRHO_P,PRHODP,POWERP, &
                                            DP1RHO,PP1RHO,WP_CEN,GPAIRP,DELTAP,EPAI_P, &
                                                                        EREA_P,IZ_FIX)

            EPAIR1=EPAI_N+EPAI_P

            ! Compute pairing energy for fragment 2
            Call pairing_field(NXHERM,NYHERM,NZHERM,PN2RHO,PP2RHO,WN_CEN,WP_CEN)

            Call interaction_energy_pairing(NXHERM,NYHERM,NZHERM,PRHO_N,PRHODN,POWERN, &
                                            DN2RHO,PN2RHO,WN_CEN,GPAIRN,DELTAN,EPAI_N, &
                                                                        EREA_N,IN_FIX)

            Call interaction_energy_pairing(NXHERM,NYHERM,NZHERM,PRHO_P,PRHODP,POWERP, &
                                            DP2RHO,PP2RHO,WP_CEN,GPAIRP,DELTAP,EPAI_P, &
                                                                        EREA_P,IZ_FIX)

            EPAIR2=EPAI_N+EPAI_P

            ! Compute interaction pairing energy: pairing field is created by fragment 1,
            ! energy is computed by taking the trace with pairing density from fragment 2
            Call pairing_field(NXHERM,NYHERM,NZHERM,PN1RHO,PP1RHO,WN_CEN,WP_CEN)

            Call interaction_energy_pairing(NXHERM,NYHERM,NZHERM,PRHO_N,PRHODN,POWERN, &
                                            DN2RHO,PN2RHO,WN_CEN,GPAIRN,DELTAN,EPAI_N, &
                                                                        EREA_N,IN_FIX)

            Call interaction_energy_pairing(NXHERM,NYHERM,NZHERM,PRHO_P,PRHODP,POWERP, &
                                            DP2RHO,PP2RHO,WP_CEN,GPAIRP,DELTAP,EPAI_P, &
                                                                        EREA_P,IZ_FIX)

            EP_INT=EPAI_N+EPAI_P

            Deallocate(WN_CEN,WP_CEN)

         End If

         !---------------------------------------------------------------------!
         !                           SKYRME ENERGIES                           !
         !---------------------------------------------------------------------!

         ! Compute Skyrme energy for fragment 1
         Call interaction_energy_Skyrme(NXHERM,NYHERM,NZHERM,ENESKY,ENEVEN,ENEODD, &
                                               ENREAR,ENE_W0,EEVEW0,EODDW0,EKESKY, &
                                               DN1RHO,DN1TAU,DN1LPR,DN1DIV,PN1RHO, &
                                 DN1SPI,DN1KIS,DN1GRR,DN1LPS,DN1ROS,DN1ROC,DN1CUR, &
                                                                    DN1SCU,DN1DES, &
                                               DP1RHO,DP1TAU,DP1LPR,DP1DIV,PP1RHO, &
                                 DP1SPI,DP1KIS,DP1GRR,DP1LPS,DP1ROS,DP1ROC,DP1CUR, &
                                                                    DP1SCU,DP1DES, &
                                               DN1RHO,DN1TAU,DN1LPR,DN1DIV,PN1RHO, &
                                 DN1SPI,DN1KIS,DN1GRR,DN1LPS,DN1ROS,DN1ROC,DN1CUR, &
                                                                    DN1SCU,DN1DES, &
                                               DP1RHO,DP1TAU,DP1LPR,DP1DIV,PP1RHO, &
                                 DP1SPI,DP1KIS,DP1GRR,DP1LPS,DP1ROS,DP1ROC,DP1CUR, &
                                                                    DP1SCU,DP1DES)
         fis_DCPALL_1(1:NDCOUP)=fis_DCPALL(1:NDCOUP)
         fis_ECPALL_1(1:NDCOUP)=fis_ECPALL(1:NDCOUP)
         fis_ECCALL_1(1:NDCOUP)=fis_ECCALL(1:NDCOUP)
         ESKFR1=Real(EKESKY)

         ! Compute Skyrme interaction energy between the two fragments: 1->2
         Call interaction_energy_Skyrme(NXHERM,NYHERM,NZHERM,ENESKY,ENEVEN,ENEODD, &
                                               ENREAR,ENE_W0,EEVEW0,EODDW0,EKESKY, &
                                               DN1RHO,DN1TAU,DN1LPR,DN1DIV,PN1RHO, &
                                 DN1SPI,DN1KIS,DN1GRR,DN1LPS,DN1ROS,DN1ROC,DN1CUR, &
                                                                    DN1SCU,DN1DES, &
                                               DP1RHO,DP1TAU,DP1LPR,DP1DIV,PP1RHO, &
                                 DP1SPI,DP1KIS,DP1GRR,DP1LPS,DP1ROS,DP1ROC,DP1CUR, &
                                                                    DP1SCU,DP1DES, &
                                               DN2RHO,DN2TAU,DN2LPR,DN2DIV,PN2RHO, &
                                 DN2SPI,DN2KIS,DN2GRR,DN2LPS,DN2ROS,DN2ROC,DN2CUR, &
                                                                    DN2SCU,DN2DES, &
                                               DP2RHO,DP2TAU,DP2LPR,DP2DIV,PP2RHO, &
                                 DP2SPI,DP2KIS,DP2GRR,DP2LPS,DP2ROS,DP2ROC,DP2CUR, &
                                                                    DP2SCU,DP2DES)
         fis_DCPALL_12(1:NDCOUP)=fis_DCPALL(1:NDCOUP)
         fis_ECPALL_12(1:NDCOUP)=fis_ECPALL(1:NDCOUP)
         fis_ECCALL_12(1:NDCOUP)=fis_ECCALL(1:NDCOUP)
         ESKINT=Real(EKESKY)

         ! Compute Skyrme interaction energy between the two fragments: 2->1
         Call interaction_energy_Skyrme(NXHERM,NYHERM,NZHERM,ENESKY,ENEVEN,ENEODD, &
                                               ENREAR,ENE_W0,EEVEW0,EODDW0,EKESKY, &
                                               DN2RHO,DN2TAU,DN2LPR,DN2DIV,PN2RHO, &
                                 DN2SPI,DN2KIS,DN2GRR,DN2LPS,DN2ROS,DN2ROC,DN2CUR, &
                                                                    DN2SCU,DN2DES, &
                                               DP2RHO,DP2TAU,DP2LPR,DP2DIV,PP2RHO, &
                                 DP2SPI,DP2KIS,DP2GRR,DP2LPS,DP2ROS,DP2ROC,DP2CUR, &
                                                                    DP2SCU,DP2DES, &
                                               DN1RHO,DN1TAU,DN1LPR,DN1DIV,PN1RHO, &
                                 DN1SPI,DN1KIS,DN1GRR,DN1LPS,DN1ROS,DN1ROC,DN1CUR, &
                                                                    DN1SCU,DN1DES, &
                                               DP1RHO,DP1TAU,DP1LPR,DP1DIV,PP1RHO, &
                                 DP1SPI,DP1KIS,DP1GRR,DP1LPS,DP1ROS,DP1ROC,DP1CUR, &
                                                                    DP1SCU,DP1DES)
         Do i=1,NDCOUP
            fis_DCPALL_12(i)=fis_DCPALL_12(i)+fis_DCPALL(i)
            fis_ECPALL_12(i)=fis_ECPALL_12(i)+fis_ECPALL(i)
            fis_ECCALL_12(i)=fis_ECCALL_12(i)+fis_ECCALL(i)
         End Do
         ESKINT=ESKINT+Real(EKESKY)

         ! Compute Skyrme energy for fragment 2
         Call interaction_energy_Skyrme(NXHERM,NYHERM,NZHERM,ENESKY,ENEVEN,ENEODD, &
                                               ENREAR,ENE_W0,EEVEW0,EODDW0,EKESKY, &
                                               DN2RHO,DN2TAU,DN2LPR,DN2DIV,PN2RHO, &
                                 DN2SPI,DN2KIS,DN2GRR,DN2LPS,DN2ROS,DN2ROC,DN2CUR, &
                                                                    DN2SCU,DN2DES, &
                                               DP2RHO,DP2TAU,DP2LPR,DP2DIV,PP2RHO, &
                                 DP2SPI,DP2KIS,DP2GRR,DP2LPS,DP2ROS,DP2ROC,DP2CUR, &
                                                                    DP2SCU,DP2DES, &
                                               DN2RHO,DN2TAU,DN2LPR,DN2DIV,PN2RHO, &
                                 DN2SPI,DN2KIS,DN2GRR,DN2LPS,DN2ROS,DN2ROC,DN2CUR, &
                                                                    DN2SCU,DN2DES, &
                                               DP2RHO,DP2TAU,DP2LPR,DP2DIV,PP2RHO, &
                                 DP2SPI,DP2KIS,DP2GRR,DP2LPS,DP2ROS,DP2ROC,DP2CUR, &
                                                                    DP2SCU,DP2DES)
         fis_DCPALL_2(1:NDCOUP)=fis_DCPALL(1:NDCOUP)
         fis_ECPALL_2(1:NDCOUP)=fis_ECPALL(1:NDCOUP)
         fis_ECCALL_2(1:NDCOUP)=fis_ECCALL(1:NDCOUP)
         ESKFR2=Real(EKESKY)

         !-----------------------------------------!
         !     Neutron kinetic energies        !
         !-----------------------------------------!

         EKIN1N=0.0_pr; EKIN2N=0.0_pr
         EKIN1P=0.0_pr; EKIN2P=0.0_pr
         Do IZ=1,NZHERM
            Do IY=1,NYHERM
               Do IX=1,NXHERM
                  W_HERM=FOURWG(IX,1)*FOURWG(IY,2)*FOURWG(IZ,3)
                  EKIN1N=EKIN1N+W_HERM*Real(DN1TAU(IX,IY,IZ))/EXPAUX(IX,IY,IZ)*HBMASS
                  EKIN2N=EKIN2N+W_HERM*Real(DN2TAU(IX,IY,IZ))/EXPAUX(IX,IY,IZ)*HBMASS
                  EKIN1P=EKIN1P+W_HERM*Real(DP1TAU(IX,IY,IZ))/EXPAUX(IX,IY,IZ)*HBMASS
                  EKIN2P=EKIN2P+W_HERM*Real(DP2TAU(IX,IY,IZ))/EXPAUX(IX,IY,IZ)*HBMASS
               End Do
            End Do
         End Do

         Deallocate(DN1RHO,DN1TAU,DN1LPR)
         Deallocate(DN1DIV,DN1SPI,DN1KIS)
         Deallocate(DN1GRR,DN1LPS,DN1ROS)
         Deallocate(DN1ROC,DN1CUR)
         Deallocate(DN1SCU,DN1DES)

         Deallocate(DP1RHO,DP1TAU,DP1LPR)
         Deallocate(DP1DIV,DP1SPI,DP1KIS)
         Deallocate(DP1GRR,DP1LPS,DP1ROS)
         Deallocate(DP1ROC,DP1CUR)
         Deallocate(DP1SCU,DP1DES)

         Deallocate(DN2RHO,DN2TAU,DN2LPR)
         Deallocate(DN2DIV,DN2SPI,DN2KIS)
         Deallocate(DN2GRR,DN2LPS,DN2ROS)
         Deallocate(DN2ROC,DN2CUR)
         Deallocate(DN2SCU,DN2DES)

         Deallocate(DP2RHO,DP2TAU,DP2LPR)
         Deallocate(DP2DIV,DP2SPI,DP2KIS)
         Deallocate(DP2GRR,DP2LPS,DP2ROS)
         Deallocate(DP2ROC,DP2CUR)
         Deallocate(DP2SCU,DP2DES)

      End If ! End compute_densities = True

    End Subroutine interaction_energy_total

    !---------------------------------------------------------------------!
    !  This subroutine computes the direct contribution to the Coulomb    !
    !  energy.                                                            !
    !---------------------------------------------------------------------!
    Subroutine compute_density(DENRHO,DENTAU,DENLPR,DENDIV,DENSPI,DENKIS, &
                               DENGRR,DENLPS,DENROS,DENROC,DENCUR,DENSCU, &
                               DENDES,PENRHO,                             &
                               ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                               ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                               IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                               IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                               ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMI, &
                               ESUMSP,IPNMIX,ITIREP,ITPNMX,MIN_QP,I_REGA, &
                               fortran_unit,type_density)

      Use HE_DEN
      Use PD_DEN
      Use DENCON

      Logical, INTENT(IN) :: PRINIT
      Integer, INTENT(IN) :: ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                             ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                             IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                             IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                             ISAWAV,IKERNE,I_TYPE,IPNMIX,ITIREP,ITPNMX, &
                             MIN_QP,I_REGA,fortran_unit
      Real(pr), INTENT(IN) :: EFERMI
      Real(pr), INTENT(INOUT) :: ESUMSP
      Character(Len=8), INTENT(IN) :: NAMEPN
      Character(Len=18), INTENT(IN) :: type_density
      Complex(pr), Allocatable, INTENT(INOUT) :: DENRHO(:,:,:),DENTAU(:,:,:),DENLPR(:,:,:),DENDIV(:,:,:),PENRHO(:,:,:)
      Complex(pr), Allocatable, INTENT(INOUT) :: DENSPI(:,:,:,:),DENKIS(:,:,:,:),DENGRR(:,:,:,:),DENLPS(:,:,:,:)
      Complex(pr), Allocatable, INTENT(INOUT) :: DENROS(:,:,:,:),DENROC(:,:,:,:),DENCUR(:,:,:,:)
      Complex(pr), Allocatable, INTENT(INOUT) :: DENSCU(:,:,:,:,:),DENDES(:,:,:,:,:)

      Integer(ipr) :: ix,iy,iz
      Complex(pr) :: C_ZERO

      C_ZERO=Cmplx(0.0_pr,0.0_pr)

      ! Get and write proton densities in r-space for left fragment and discrete states only
      Call filter_density(EFERMI,ICHARG,I_TYPE)

      DE_RHO(:,:,:,ICHARG)   = C_ZERO; DE_TAU(:,:,:,ICHARG) = C_ZERO
      DE_LPR(:,:,:,ICHARG)   = C_ZERO; DE_DIV(:,:,:,ICHARG) = C_ZERO
      PD_RHO(:,:,:,ICHARG)   = C_ZERO
      DE_SPI(:,:,:,:,ICHARG) = C_ZERO; DE_KIS(:,:,:,:,ICHARG) = C_ZERO
      DE_GRR(:,:,:,:,ICHARG) = C_ZERO; DE_LPS(:,:,:,:,ICHARG) = C_ZERO
      DE_ROS(:,:,:,:,ICHARG) = C_ZERO; DE_ROC(:,:,:,:,ICHARG) = C_ZERO
      DE_CUR(:,:,:,:,ICHARG) = C_ZERO
      DE_SCU(:,:,:,:,:,ICHARG) = C_ZERO; DE_DES(:,:,:,:,:,ICHARG) = C_ZERO

      Call DENSHF(NXHERM,NYHERM,NZHERM,ISIMTX,JSIMTY,ISIMTZ, &
                  ISIGNY,ISIMPY,ISIQTY,IPAHFB,MREVER,ICHARG, &
                                       MIN_QP,IPNMIX,ITPNMX, &
                  ITIREP,NXMAXX,NAMEPN,PRINIT,IDEVAR,ITERUN, &
                  ISYMDE,INIROT,INIINV,INIKAR,ISAWAV,IKERNE)

      Call DBLING(ITPNMX,NXHERM,NYHERM,NZHERM)

      DENRHO(:,:,:)=DE_RHO(:,:,:,ICHARG); DENTAU(:,:,:)=DE_TAU(:,:,:,ICHARG)
      DENLPR(:,:,:)=DE_LPR(:,:,:,ICHARG); DENDIV(:,:,:)=DE_DIV(:,:,:,ICHARG)
      PENRHO(:,:,:)=PD_RHO(:,:,:,ICHARG)
      DENSPI(:,:,:,:)=DE_SPI(:,:,:,:,ICHARG); DENKIS(:,:,:,:)=DE_KIS(:,:,:,:,ICHARG)
      DENGRR(:,:,:,:)=DE_GRR(:,:,:,:,ICHARG); DENLPS(:,:,:,:)=DE_LPS(:,:,:,:,ICHARG)
      DENROS(:,:,:,:)=DE_ROS(:,:,:,:,ICHARG); DENROC(:,:,:,:)=DE_ROC(:,:,:,:,ICHARG)
      DENCUR(:,:,:,:)=DE_CUR(:,:,:,:,ICHARG)
      DENSCU(:,:,:,:,:)=DE_SCU(:,:,:,:,:,ICHARG)
      DENDES(:,:,:,:,:)=DE_DES(:,:,:,:,:,ICHARG)

      Allocate(DCONTI(1:NXHERM,1:NYHERM,1:NZHERM))
      Do iz=1,NZHERM
         Do iy=1,NYHERM
            Do ix=1,NXHERM
               DCONTI(ix,iy,iz)=0.0_pr
            End Do
         End Do
      End Do
      Call write_density(fortran_unit, DENRHO, DCONTI, NXHERM, NYHERM, NZHERM, type_density)
      Deallocate(DCONTI)

    End Subroutine compute_density

    !---------------------------------------------------------------------!
    !  This subroutine computes the direct contribution to the Coulomb    !
    !  energy.                                                            !
    !---------------------------------------------------------------------!
    Subroutine save_density(ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                            ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                            IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                            IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                            ISAWAV,IKERNE,I_TYPE,NAMEPN,PRINIT,EFERMI, &
                            ESUMSP,IPNMIX,ITIREP,MIN_QP,I_REGA,fichier)

      Use SAVQUA
      Use ALLQUA
      Use HE_DEN
      Use PD_DEN
      Use DENCON

      Logical, INTENT(IN) :: PRINIT
      Integer, INTENT(IN) :: ICOUDI,ICOUEX,I_YUKA,NXHERM,NYHERM,NZHERM, &
                             ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
                             IPAHFB,MREVER,ICHARG,I_GOGA,IGOGPA,NXMAXX, &
                             IDEVAR,ITERUN,ISYMDE,INIROT,INIINV,INIKAR, &
                             ISAWAV,IKERNE,I_TYPE,IPNMIX,ITIREP,I_REGA, &
                             MIN_QP
      Real(pr), INTENT(IN) :: EFERMI
      Real(pr), INTENT(INOUT) :: ESUMSP
      Character(Len=8), INTENT(IN) :: NAMEPN
      Character(Len=18), INTENT(IN) :: fichier

      Integer(ipr) :: ix,iy,iz,ITPNMX,error,NOTOCC,NOCISO
      Real(pr) :: VCOISO,V_CORR
      Complex(pr) :: C_ZERO
      Complex(pr), Allocatable :: DENRHO(:,:,:)

      COMMON                                             &
             /LIMOCC/ NOTOCC(1:NDSTAT,0:NDREVE)          &
             /OCCISO/ NOCISO(1:NDSTAT,0:NDREVE,0:NDISOS)
      COMMON                                             &
             /CORPAR/ V_CORR(1:NDSTAT,0:NDREVE)          &
             /CORISO/ VCOISO(1:NDSTAT,0:NDREVE,0:NDISOS)

      C_ZERO=Cmplx(0.0_pr,0.0_pr)

      ! Get and write proton densities in r-space for left fragment and discrete states only
      Call filter_density(EFERMI,ICHARG,I_TYPE)

      DE_RHO(:,:,:,ICHARG)   = C_ZERO; DE_TAU(:,:,:,ICHARG) = C_ZERO
      DE_LPR(:,:,:,ICHARG)   = C_ZERO; DE_DIV(:,:,:,ICHARG) = C_ZERO
      PD_RHO(:,:,:,ICHARG)   = C_ZERO
      DE_SPI(:,:,:,:,ICHARG) = C_ZERO; DE_KIS(:,:,:,:,ICHARG) = C_ZERO
      DE_GRR(:,:,:,:,ICHARG) = C_ZERO; DE_LPS(:,:,:,:,ICHARG) = C_ZERO
      DE_ROS(:,:,:,:,ICHARG) = C_ZERO; DE_ROC(:,:,:,:,ICHARG) = C_ZERO
      DE_CUR(:,:,:,:,ICHARG) = C_ZERO
      DE_SCU(:,:,:,:,:,ICHARG) = C_ZERO; DE_DES(:,:,:,:,:,ICHARG) = C_ZERO

      ITPNMX=ICHARG

      BWAQUA(:,:,:)=BSVQUA(:,:,:,ICHARG)
      AWAQUA(:,:,:)=ASVQUA(:,:,:,ICHARG)
      V_CORR(:,:)=1.0_pr
      NOTOCC(:,:)=0

      Allocate(DENRHO(NDXHRM,NDYHRM,NDZHRM))

      Call DENSHF(NXHERM,NYHERM,NZHERM,ISIMTX,JSIMTY,ISIMTZ, &
                  ISIGNY,ISIMPY,ISIQTY,IPAHFB,MREVER,ICHARG, &
                                       MIN_QP,IPNMIX,ITPNMX, &
                  ITIREP,NXMAXX,NAMEPN,PRINIT,IDEVAR,ITERUN, &
                  ISYMDE,INIROT,INIINV,INIKAR,ISAWAV,IKERNE)

      Call DBLING(ITPNMX,NXHERM,NYHERM,NZHERM)

      DENRHO(:,:,:)=DE_RHO(:,:,:,ICHARG)

      Allocate(DCONTI(1:NXHERM,1:NYHERM,1:NZHERM))
      Do iz=1,NZHERM
         Do iy=1,NYHERM
            Do ix=1,NXHERM
               DCONTI(ix,iy,iz)=0.0_pr
            End Do
         End Do
      End Do

      Open(Unit=54,File=fichier,Status='Unknown',Form='Formatted',Iostat=error)
      Call write_density(54, DENRHO, DCONTI, NXHERM, NYHERM, NZHERM, 'n_all_spectrum.dat')
      Close(54)

      Deallocate(DCONTI)
      Deallocate(DENRHO)

    End Subroutine save_density

    !---------------------------------------------------------------------!
    !  This subroutine backups all densities coming from HFODD            !
    !---------------------------------------------------------------------!
    Subroutine backup_density()

      Use PD_DEN
      Use HE_DEN

      Integer(ipr) :: ialloc

      Allocate(backup_DENRHO(1:NDXHRM,1:NDYHRM,1:NDZHRM,0:NDISOS),STAT=IALLOC)
      If(IALLOC.Ne.0) Write(6,'("Error 00")')
      Allocate(backup_DENTAU(1:NDXHRM,1:NDYHRM,1:NDZHRM,0:NDISOS),STAT=IALLOC)
      If(IALLOC.Ne.0) Write(6,'("Error 01")')
      Allocate(backup_DENLPR(1:NDXHRM,1:NDYHRM,1:NDZHRM,0:NDISOS),STAT=IALLOC)
      If(IALLOC.Ne.0) Write(6,'("Error 02")')
      Allocate(backup_DENDIV(1:NDXHRM,1:NDYHRM,1:NDZHRM,0:NDISOS),STAT=IALLOC)
      If(IALLOC.Ne.0) Write(6,'("Error 03")')
      Allocate(backup_PENRHO(1:NDXHRM,1:NDYHRM,1:NDZHRM,0:NDISOS),STAT=IALLOC)
      If(IALLOC.Ne.0) Write(6,'("Error 04")')
      Allocate(backup_DENSPI(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,0:NDISOS),STAT=IALLOC)
      If(IALLOC.Ne.0) Write(6,'("Error 05")')
      Allocate(backup_DENKIS(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,0:NDISOS),STAT=IALLOC)
      If(IALLOC.Ne.0) Write(6,'("Error 06")')
      Allocate(backup_DENGRR(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,0:NDISOS),STAT=IALLOC)
      If(IALLOC.Ne.0) Write(6,'("Error 00")')
      Allocate(backup_DENLPS(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,0:NDISOS),STAT=IALLOC)
      If(IALLOC.Ne.0) Write(6,'("Error 07")')
      Allocate(backup_DENROS(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,0:NDISOS),STAT=IALLOC)
      If(IALLOC.Ne.0) Write(6,'("Error 08")')
      Allocate(backup_DENROC(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,0:NDISOS),STAT=IALLOC)
      If(IALLOC.Ne.0) Write(6,'("Error 09")')
      Allocate(backup_DENCUR(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,0:NDISOS),STAT=IALLOC)
      If(IALLOC.Ne.0) Write(6,'("Error 10")')
      Allocate(backup_DENSCU(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,1:NDKART,0:NDISOS),STAT=IALLOC)
      If(IALLOC.Ne.0) Write(6,'("Error 11")')
      Allocate(backup_DENDES(1:NDXHRM,1:NDYHRM,1:NDZHRM,1:NDKART,1:NDKART,0:NDISOS),STAT=IALLOC)
      If(IALLOC.Ne.0) Write(6,'("Error 12")')

      backup_DENRHO(:,:,:,:)    =DE_RHO(:,:,:,:)
      backup_DENTAU(:,:,:,:)    =DE_TAU(:,:,:,:)
      backup_DENLPR(:,:,:,:)    =DE_LPR(:,:,:,:)
      backup_DENDIV(:,:,:,:)    =DE_DIV(:,:,:,:)
      backup_PENRHO(:,:,:,:)    =PD_RHO(:,:,:,:)
      backup_DENSPI(:,:,:,:,:)  =DE_SPI(:,:,:,:,:)
      backup_DENKIS(:,:,:,:,:)  =DE_KIS(:,:,:,:,:)
      backup_DENGRR(:,:,:,:,:)  =DE_GRR(:,:,:,:,:)
      backup_DENLPS(:,:,:,:,:)  =DE_LPS(:,:,:,:,:)
      backup_DENROS(:,:,:,:,:)  =DE_ROS(:,:,:,:,:)
      backup_DENROC(:,:,:,:,:)  =DE_ROC(:,:,:,:,:)
      backup_DENCUR(:,:,:,:,:)  =DE_CUR(:,:,:,:,:)
      backup_DENSCU(:,:,:,:,:,:)=DE_SCU(:,:,:,:,:,:)
      backup_DENDES(:,:,:,:,:,:)=DE_DES(:,:,:,:,:,:)

    End Subroutine backup_density

    !---------------------------------------------------------------------!
    !  This subroutine restores all backup densities to their regular     !
    ! storage arrays
    !---------------------------------------------------------------------!
    Subroutine restore_density()

      Use PD_DEN
      Use HE_DEN

      DE_RHO(:,:,:,:)     = backup_DENRHO(:,:,:,:)
      DE_TAU(:,:,:,:)     = backup_DENTAU(:,:,:,:)
      DE_LPR(:,:,:,:)     = backup_DENLPR(:,:,:,:)
      DE_DIV(:,:,:,:)     = backup_DENDIV(:,:,:,:)
      PD_RHO(:,:,:,:)     = backup_PENRHO(:,:,:,:)
      DE_SPI(:,:,:,:,:)   = backup_DENSPI(:,:,:,:,:)
      DE_KIS(:,:,:,:,:)   = backup_DENKIS(:,:,:,:,:)
      DE_GRR(:,:,:,:,:)   = backup_DENGRR(:,:,:,:,:)
      DE_LPS(:,:,:,:,:)   = backup_DENLPS(:,:,:,:,:)
      DE_ROS(:,:,:,:,:)   = backup_DENROS(:,:,:,:,:)
      DE_ROC(:,:,:,:,:)   = backup_DENROC(:,:,:,:,:)
      DE_CUR(:,:,:,:,:)   = backup_DENCUR(:,:,:,:,:)
      DE_SCU(:,:,:,:,:,:) = backup_DENSCU(:,:,:,:,:,:)
      DE_DES(:,:,:,:,:,:) = backup_DENDES(:,:,:,:,:,:)

      table_singles(:,:) = 1
      number_of_pairs = 0

      Deallocate(backup_DENRHO,backup_DENTAU,backup_DENLPR,backup_DENDIV,backup_PENRHO, &
                 backup_DENSPI,backup_DENKIS,backup_DENGRR,backup_DENLPS,backup_DENROS, &
                 backup_DENROC,backup_DENCUR,backup_DENSCU,backup_DENDES)

    End Subroutine restore_density

    !---------------------------------------------------------------------!
    !  This subroutine computes the direct contribution to the Coulomb    !
    !  energy.                                                            !
    !---------------------------------------------------------------------!
    Subroutine interaction_energy_Coulomb(NXHERM,NYHERM,NZHERM,NXMAXX,NYMAXX,NZMAXX, &
                                          NUMCOU,BOUCOU,ISIMPY,IKERNE,KFRAGM,density, &
                                                                             EKECOD)

      Use HCOULO
      Use MAT_PP

      Integer, INTENT(IN) :: NXHERM,NYHERM,NZHERM,NXMAXX,NYMAXX,NZMAXX,NUMCOU, &
                             ISIMPY,IKERNE,KFRAGM
      Real(pr), INTENT(IN) :: BOUCOU
      Complex(pr), Allocatable, INTENT(INOUT) :: density(:,:,:)
      Complex(pr), INTENT(INOUT) :: EKECOD

      Integer(ipr) :: LDBASE,NFIPRI,NDEGRE,ISURFA
      Real(pr) :: HERONE,TWOWGT,TWOPNT,ECHAR2,PIOMAS
      Complex, Dimension(0:ND2MAX,0:ND2MAX,0:ND2MAX) :: SONBAS

      Integer(ipr) :: i,j
      Complex(pr) :: result

      COMMON                 &
             /DIMENS/ LDBASE
      COMMON                                            &
             /HERMEM/ HERONE(0:ND2MAX,1:NDGAUS,1:NDKART)
      COMMON                                    &
             /INTMEM/ TWOWGT(1:NDGAUS,1:NDKART),&
                      TWOPNT(1:NDGAUS,1:NDKART)
      COMMON                 &
             /CHARGE/ ECHAR2
      COMMON                                             &
             /CFIPRI/ NFIPRI

      If(.NOT.Allocated(HPPCOU)) Then
         Allocate(HPPCOU(1:NDBASE,1:NDBASE))
      End If
      If(.NOT.Allocated(HPMCOU)) Then
         Allocate(HPMCOU(1:NDBASE,1:NDBASE))
      End If

      NDEGRE=0
      Call BEGINC(density,SONBAS,NDEGRE,NXHERM,NYHERM,NZHERM,HERONE,TWOWGT, &
                                                      NXMAXX,NYMAXX,NZMAXX)

      PIOMAS=0.0_pr; ISURFA=1
      Call INTCOU(NXMAXX,NYMAXX,NZMAXX,NUMCOU,BOUCOU,SONBAS,ECHAR2,PIOMAS, &
                                ISURFA,IKERNE,ISIMPY,KFRAGM,HPPCOU,HPMCOU)

      result=Cmplx(0.0_pr,0.0_pr)
      Do i=1,LDBASE
         Do j=1,LDBASE
            result=result+HPPCOU(i,j) * (BIG_PP(j,i,0)+BIG_PP(j,i,1))
         End Do
      End Do

      EKECOD=0.5_pr*result

    End Subroutine interaction_energy_Coulomb

    !---------------------------------------------------------------------!
    !  This subroutine computes the exchange contribution to the Coulomb  !
    !  energy by using the Slater approximation. It is a variant of the   !
    !  routine TRUTOD of HFODD where the density is passed as argument to !
    !  the routine.
    !---------------------------------------------------------------------!
    Subroutine exchange_Coulomb(NXHERM,NYHERM,NZHERM,density,EKECOE)

      Integer, INTENT(IN) :: NXHERM,NYHERM,NZHERM
      Complex(pr), Allocatable, INTENT(IN) :: density(:,:,:)
      Complex(pr), INTENT(INOUT) :: EKECOE

      Real(pr) :: EXPAUX,FOURWG,FOURPT,COULEX

      Integer(ipr) :: IX,IY,IZ
      Real(pr) :: W_HERM
      Complex(pr) :: DENHAM,DENTST

      COMMON                                      &
             /DENEXP/ EXPAUX(NDXHRM,NDYHRM,NDZHRM)
      COMMON                                     &
             /INTSTO/ FOURWG(1:NDGAUS,1:NDKART), &
                      FOURPT(1:NDGAUS,1:NDKART)
      COMMON                 &
             /EXCHAN/ COULEX

      EKECOE=Cmplx(0.0_pr,0.0_pr)
      Do IZ=1,NZHERM
         Do IY=1,NYHERM
            Do IX=1,NXHERM
               DENTST = density(IX,IY,IZ)*EXPAUX(IX,IY,IZ)
               DENHAM = DENTST ** (4.0_pr/3.0_pr) / EXPAUX(IX,IY,IZ) / EXPAUX(IX,IY,IZ)
               W_HERM=FOURWG(IX,1)*FOURWG(IY,2)*FOURWG(IZ,3)
               EKECOE=EKECOE+W_HERM*DENHAM
            End Do
         End Do
      End Do

      EKECOE = COULEX * EKECOE

    End Subroutine exchange_Coulomb

    !---------------------------------------------------------------------!
    ! This subroutine calculates the pairing energy of a fragment, or the !
    ! pairing interaction energy between two fragments. The pairing field !
    ! is stored in PAICEN, the local density for particle t \rho_t(r) in  !
    ! DENRHO (without Gaussian factors), the total local density \rho(r)  !
    ! in density, and the pairing density for particle t in PAIRHO.       !
    !---------------------------------------------------------------------!
    Subroutine interaction_energy_pairing(NXHERM,NYHERM,NZHERM,PRHO_X,PRHODX,POWERX, &
                                          DENRHO,PAIRHO,PAICEN,G_PAIR,PDELTA,E_PAIR, &
                                                                      E_REAR,NOPART)
      Integer, INTENT(IN) :: NXHERM,NYHERM,NZHERM,NOPART
      Real(pr), INTENT(IN) :: PRHO_X,PRHODX,POWERX
      Real(pr), INTENT(INOUT) :: G_PAIR,PDELTA,E_PAIR,E_REAR
      Complex(pr), Allocatable, INTENT(IN) :: DENRHO(:,:,:),PAIRHO(:,:,:),PAICEN(:,:,:)

      Real(pr) :: FOURWG,FOURPT,EXPAUX
      Complex(pr) :: DENSIC,DENCHC

      Integer(ipr) :: IX,IY,IZ
      Real(pr) :: W_HERM,QRHO_X,QRHODX,QCEN_X,DENPOW,QUXRHO,QUXCEN,RRHO_X,RRHODX

      COMMON                                       &
             /DENTOC/ DENSIC(NDXHRM,NDYHRM,NDZHRM),&
                      DENCHC(NDXHRM,NDYHRM,NDZHRM)
      COMMON                                     &
             /INTSTO/ FOURWG(1:NDGAUS,1:NDKART), &
                      FOURPT(1:NDGAUS,1:NDKART)
      COMMON                                       &
             /DENEXP/ EXPAUX(NDXHRM,NDYHRM,NDZHRM)

      QRHO_X=0.0_pr; QRHODX=0.0_pr; QCEN_X=0.0_pr

      Do IZ=1,NZHERM
         Do IY=1,NYHERM
            Do IX=1,NXHERM

               ! Calculating the auxiliary contractions

               QUXRHO = Real(PAIRHO(IX,IY,IZ)*Conjg(PAIRHO(IX,IY,IZ)))
               QUXCEN = Real(DENRHO(IX,IY,IZ)*Conjg(PAICEN(IX,IY,IZ)))

               W_HERM = FOURWG(IX,1)*FOURWG(IY,2)*FOURWG(IZ,3)

               DENPOW = Real(DENSIC(IX,IY,IZ)) ** POWERX

               ! Calculating the t0-terms of the pairing energy density

               QRHO_X = QRHO_X + W_HERM * QUXRHO
               QCEN_X = QCEN_X + W_HERM * QUXCEN

               ! Calculating the t3-terms of the pairing energy density

               QRHODX = QRHODX + W_HERM * QUXRHO * DENPOW

               E_REAR = E_REAR + W_HERM * PRHODX * QUXRHO * 0.25_pr * DENPOW * (POWERX/2)

            End Do
         End Do
      End Do

      ! Summations of terms in energy.
      ! The division by 4 below corresponds to the factor of 1/4 which
      ! appears in the expressions for the pairing energy in terms of
      ! the t0*(1-x0) coupling constant

      RRHO_X = 0.25_pr * PRHO_X * QRHO_X
      RRHODX = 0.25_pr * PRHODX * QRHODX

      E_PAIR = RRHO_X+RRHODX
      PDELTA =-QCEN_X/NOPART

     ! This very crude estimate of the pairing strength is based on the simple BCS formula
      G_PAIR = 0.0_pr; If(E_PAIR.Ne.0.0_pr) G_PAIR = -PDELTA**2/E_PAIR

    End Subroutine interaction_energy_pairing

    !---------------------------------------------------------------------!
    ! This subroutine calculates the pairing fields.                      !
    ! Attention: at present, only the central fields coming from the t0   !
    ! and t3 terms are programed. The effective coupling constants are    !
    ! coded as:                                                           !
    !                  'P'//'XXX'//'N'//'TYP'                             !
    !                                                                     !
    ! where the symbols:  'P', 'XXX', 'N', and 'TYP' are explained in the !
    ! subroutine "ESKYRM"                                                 !
    !                                                                     !
    ! Ref. Nucl. Phys. A422 (1984) 103                                    !
    !---------------------------------------------------------------------!
    Subroutine pairing_field(NXHERM,NYHERM,NZHERM,PN_RHO,PP_RHO,WN_CEN,WP_CEN)

      Integer, INTENT(IN) :: NXHERM,NYHERM,NZHERM
      Complex(pr), Allocatable, INTENT(IN) :: PN_RHO(:,:,:),PP_RHO(:,:,:)
      Complex(pr), Allocatable, INTENT(INOUT) :: WN_CEN(:,:,:),WP_CEN(:,:,:)

      Real(pr) :: PRHO_N,PRHODN,PRHOSN,POWERN,PRHO_P,PRHODP,PRHOSP,POWERP
      Complex(pr) :: DENSIC,DENCHC

      Integer(ipr) :: IX,IY,IZ
      Complex(pr) :: PNEUTR,PPROTO,XNEUTR,XPROTO

      COMMON                                       &
             /DENTOC/ DENSIC(NDXHRM,NDYHRM,NDZHRM),&
                      DENCHC(NDXHRM,NDYHRM,NDZHRM)
      COMMON                                       &
             /CCPPAI/ PRHO_N,PRHODN,PRHOSN,POWERN, &
                      PRHO_P,PRHODP,PRHOSP,POWERP

      Do IZ=1,NZHERM
         Do IY=1,NYHERM
            Do IX=1,NXHERM

               ! Central term (<= t0) - (See ref Eq.(A5.b))
               PNEUTR = PRHO_N * PN_RHO(IX,IY,IZ) * 0.5_pr
               PPROTO = PRHO_P * PP_RHO(IX,IY,IZ) * 0.5_pr

               ! Central term (<= t3)
               XNEUTR = PRHODN * PN_RHO(IX,IY,IZ) * 0.5_pr * DENSIC(IX,IY,IZ)**POWERN
               XPROTO = PRHODP * PP_RHO(IX,IY,IZ) * 0.5_pr * DENSIC(IX,IY,IZ)**POWERP

               WN_CEN(IX,IY,IZ) = PNEUTR+XNEUTR
               WP_CEN(IX,IY,IZ) = PPROTO+XPROTO

            End Do
         End Do
      End Do

    End Subroutine pairing_field

    !---------------------------------------------------------------------!
    ! Subroutine used to test the normalization of the fragment densities !
    ! in coordinate space (on the Gauss-Hermite mesh)                     !
    !---------------------------------------------------------------------!
    Subroutine TSTDEN(NXHERM,NYHERM,NZHERM,D1_RHO,D2_RHO)

      Integer(ipr), INTENT(IN) :: NXHERM,NYHERM,NZHERM
      Complex(pr), Allocatable, INTENT(INOUT) :: D1_RHO(:,:,:),D2_RHO(:,:,:)

      Real(pr) :: EXPAUX,FOURWG,FOURPT

      Integer(ipr) :: IX,IY,IZ
      Real(pr) :: DENS_1,DENS_2,DENS_T,W_HERM,FRAGM1,FRAGM2,DTOTAL

      COMMON                                        &
             /DENEXP/ EXPAUX(NDXHRM,NDYHRM,NDZHRM)
      COMMON                                        &
             /INTSTO/ FOURWG(1:NDGAUS,1:NDKART),    &
                      FOURPT(1:NDGAUS,1:NDKART)

      FRAGM1=0.0_pr; FRAGM2=0.0_pr; DTOTAL=0.0_pr
      Do IZ=1,NZHERM
         Do IY=1,NYHERM
            Do IX=1,NXHERM

               DENS_1=Real(D1_RHO(IX,IY,IZ))/EXPAUX(IX,IY,IZ)
               DENS_2=Real(D2_RHO(IX,IY,IZ))/EXPAUX(IX,IY,IZ)
               DENS_T=DENS_1+DENS_2

               W_HERM=FOURWG(IX,1)*FOURWG(IY,2)*FOURWG(IZ,3)

               FRAGM1=FRAGM1+W_HERM*DENS_1
               FRAGM2=FRAGM2+W_HERM*DENS_2
               DTOTAL=DTOTAL+W_HERM*DENS_T

            End Do
         End Do
      End Do
      Write(6,'("Fragment 1: ", f20.14," - Fragment 2: ", f20.14," Total = ",f20.14)') &
                 FRAGM1,FRAGM2,DTOTAL

    End Subroutine TSTDEN

    !---------------------------------------------------------------------!
    ! Subroutine used to test the normalization of the fragment densities !
    ! in coordinate space (on the Gauss-Hermite mesh)                     !
    !---------------------------------------------------------------------!
    Subroutine test_density(NXHERM,NYHERM,NZHERM)

      Integer(ipr), INTENT(IN) :: NXHERM,NYHERM,NZHERM

      Real(pr) :: DENSIT,DENCHA,EXPAUX,FOURWG,FOURPT

      Integer(ipr) :: kx,ky,kz

      Real(pr) :: TOTDEN,CHADEN,W_HERM,DENLOC

      COMMON                                        &
             /DENTOT/ DENSIT(NDXHRM,NDYHRM,NDZHRM), &
                      DENCHA(NDXHRM,NDYHRM,NDZHRM)
      COMMON                                        &
             /DENEXP/ EXPAUX(NDXHRM,NDYHRM,NDZHRM)
      COMMON                                        &
             /INTSTO/ FOURWG(1:NDGAUS,1:NDKART),    &
                      FOURPT(1:NDGAUS,1:NDKART)

      TOTDEN=0.0_pr; CHADEN=0.0_pr
      Do kz=1,NZHERM
         Do ky=1,NYHERM
            Do kx=1,NXHERM
               ! Total density
               DENLOC=DENSIT(kx,ky,kz)/EXPAUX(kx,ky,kz)/EXPAUX(kx,ky,kz)
               W_HERM=FOURWG(kx,1)*FOURWG(ky,2)*FOURWG(kz,3)
               TOTDEN=TOTDEN+W_HERM*DENLOC
               ! Proton density
               DENLOC=DENCHA(kx,ky,kz)/EXPAUX(kx,ky,kz)/EXPAUX(kx,ky,kz)
               W_HERM=FOURWG(kx,1)*FOURWG(ky,2)*FOURWG(kz,3)
               CHADEN=CHADEN+W_HERM*DENLOC
            End Do
         End Do
      End Do
      Write(6,'("Total density DENSIT: ",f20.14)') TOTDEN
      Write(6,'("Total charge DENCHA: ",f20.14)') CHADEN

    End Subroutine test_density

 End Module hfodd_fission_interaction

 !----------------------------------------------------------------------!
 !  This module only prints some stuff on screen                        !
 !----------------------------------------------------------------------!
 Module hfodd_fission_print

    Use hfodd_sizes
    Use hfodd_fission_precision
    Use hfodd_fission_fragments
    Use hfodd_fission_interaction

    Implicit None

 Contains

    !---------------------------------------------------------------------!
    !    This subroutine prints the properties of the fragments at        !
    !    convergence. Left fragment has index 0, right fragment index 1   !
    !---------------------------------------------------------------------!
    Subroutine PRIFRA(IN_FIX,IZ_FIX,ZPOINT,QPOINT,CENLEF,CENRIG,RIGMAS,RIGCHA,ALLQLM,ORIQLM,DISFRA,NMUPRI)

      Use hfodd_sizes

      Integer(ipr), INTENT(IN) :: IZ_FIX,IN_FIX,NMUPRI
      Real(pr), INTENT(IN) :: ZPOINT,CENLEF,CENRIG,RIGMAS,RIGCHA,DISFRA,QPOINT
      Real(pr), Dimension(0:NDMULT,-NDMULT:NDMULT,0:1), INTENT(IN) :: ALLQLM,ORIQLM

      Integer(ipr) :: NFIPRI

      Character(Len=6), Dimension(NDCOUP/2) :: NCPALL
      Integer(ipr) :: LAMBDA,MIU,i,KARTEZ,KORDER
      Real(pr) :: A_MASS,DUMMY,ETOT_1,ETOT_2,protons_neck
      Real(pr) :: TPSQU1,TPSQU2,HBMREN
      Real(pr) :: sumEven_scalar_1,sumEven_vector_1,sumEven_scalar_2,sumEven_vector_2,sumEven_scalar_12,sumEven_vector_12
      Real(pr) :: sumEven_total_12,sumEven_sum_12
      Real(pr) :: sumOdd_scalar_1,sumOdd_vector_1,sumOdd_scalar_2,sumOdd_vector_2,sumOdd_scalar_12,sumOdd_vector_12
      Real(pr) :: sumOdd_total_12,sumOdd_sum_12

      COMMON                 &
             /CFIPRI/ NFIPRI

      COMMON                 &
             /RENORM/ HBMREN(NDKART)

      A_MASS=Real(IN_FIX+IZ_FIX,Kind=pr)

      EKIN1N=EKIN1N*(1.0_pr - 1.0_pr/(A_MASS-RIGMAS))
      EKIN1P=EKIN1P*(1.0_pr - 1.0_pr/(A_MASS-RIGMAS))
      EKIN2N=EKIN2N*(1.0_pr - 1.0_pr/RIGMAS)
      EKIN2P=EKIN2P*(1.0_pr - 1.0_pr/RIGMAS)

      ETOT_1=EKIN1N+EKIN1P+ESKFR1+EPAIR1+ECOUL1+ECOEX1
      ETOT_2=EKIN2N+EKIN2P+ESKFR2+EPAIR2+ECOUL2+ECOEX2

      TPSQU1=HBMREN(1)*TPSQA1(1)+HBMREN(2)*TPSQA1(2)+HBMREN(3)*TPSQA1(3)
      TPSQU2=HBMREN(1)*TPSQA2(1)+HBMREN(2)*TPSQA2(2)+HBMREN(3)*TPSQA2(3)

      NCPALL(1) ='ERHO_ '; NCPALL(2) ='ERHOD_'; NCPALL(3) ='ELPR_ '
      NCPALL(4) ='ETAU_ '; NCPALL(5) ='ESCU_ '; NCPALL(6) ='EDIV_ '
      NCPALL(7) ='ESPI_ '; NCPALL(8) ='ESPID_'; NCPALL(9) ='ELPS_ '
      NCPALL(10)='ECUR_ '; NCPALL(11)='EKIS_ '; NCPALL(12)='EROT_ '

      Do i=1,NDCOUP/2
         fis_EISALL_1(2*i-1)=   fis_DCPALL_1(2*i-1)* (fis_CCPALL(2*i-1)   + fis_CCPALL(2*i)/2.0_pr)
         fis_EISALL_1(2*i  )=(- fis_DCPALL_1(2*i-1)+  fis_DCPALL_1(2*i)*2)* fis_CCPALL(2*i)/2.0_pr
         fis_EISALL_2(2*i-1)=   fis_DCPALL_2(2*i-1)* (fis_CCPALL(2*i-1)   + fis_CCPALL(2*i)/2.0_pr)
         fis_EISALL_2(2*i  )=(- fis_DCPALL_2(2*i-1)+  fis_DCPALL_2(2*i)*2)* fis_CCPALL(2*i)/2.0_pr
         fis_EISALL_12(2*i-1)=  fis_DCPALL_12(2*i-1)*(fis_CCPALL(2*i-1)    +fis_CCPALL(2*i)/2.0_pr)
         fis_EISALL_12(2*i  )=(-fis_DCPALL_12(2*i-1)+ fis_DCPALL_12(2*i)*2)*fis_CCPALL(2*i)/2.0_pr
      End Do

      sumEven_scalar_1=0.0_pr
      sumEven_vector_1=0.0_pr
      sumEven_scalar_2=0.0_pr
      sumEven_vector_2=0.0_pr
      sumEven_scalar_12=0.0_pr
      sumEven_vector_12=0.0_pr
      sumEven_total_12=0.0_pr
      sumEven_sum_12=0.0_pr
      Do i=1,NDCOUP/4
         sumEven_scalar_1  = sumEven_scalar_1  + fis_EISALL_1(2*i-1)
         sumEven_vector_1  = sumEven_vector_1  + fis_EISALL_1(2*i  )
         sumEven_scalar_2  = sumEven_scalar_2  + fis_EISALL_2(2*i-1)
         sumEven_vector_2  = sumEven_vector_2  + fis_EISALL_2(2*i  )
         sumEven_scalar_12 = sumEven_scalar_12 + fis_EISALL_12(2*i-1)
         sumEven_vector_12 = sumEven_vector_12 + fis_EISALL_12(2*i  )
         sumEven_total_12  = sumEven_total_12  + fis_ECPALL_12(2*i-1)
         sumEven_sum_12    = sumEven_sum_12    + fis_ECPALL_12(2*i  )
      End Do

      sumOdd_scalar_1=0.0_pr
      sumOdd_vector_1=0.0_pr
      sumOdd_scalar_2=0.0_pr
      sumOdd_vector_2=0.0_pr
      sumOdd_scalar_12=0.0_pr
      sumOdd_vector_12=0.0_pr
      sumOdd_total_12=0.0_pr
      sumOdd_sum_12=0.0_pr
      Do i=NDCOUP/4+1,NDCOUP/2
         sumOdd_scalar_1  = sumOdd_scalar_1  + fis_EISALL_1(2*i-1)
         sumOdd_vector_1  = sumOdd_vector_1  + fis_EISALL_1(2*i  )
         sumOdd_scalar_2  = sumOdd_scalar_2  + fis_EISALL_2(2*i-1)
         sumOdd_vector_2  = sumOdd_vector_2  + fis_EISALL_2(2*i  )
         sumOdd_scalar_12 = sumOdd_scalar_12 + fis_EISALL_12(2*i-1)
         sumOdd_vector_12 = sumOdd_vector_12 + fis_EISALL_12(2*i  )
         sumOdd_total_12  = sumOdd_total_12  + fis_ECPALL_12(2*i-1)
         sumOdd_sum_12    = sumOdd_sum_12    + fis_ECPALL_12(2*i  )
      End Do

      KARTEZ=3;KORDER=0
      protons_neck = protons_in_neck(KARTEZ,ZPOINT,KORDER)

      Write(NFIPRI,'(79("*"),/,"*",77X,"*")')
      Write(NFIPRI,'("*",12X,"FISSION FRAGMENT PROPERTIES (SHARP ADIABATIC SCISSION)",11X,"*")')
      Write(NFIPRI,'("*",77X,"*")')
      Write(NFIPRI,'("*   POSITION OF THE NECK ALONG THE Z-AXIS..... Zn : ",F8.4," (DimensionLESS)  *")') &
                          ZPOINT
      Write(NFIPRI,'("*   SIZE OF THE NECK ALONG THE Z-AXIS......... Qn : ",F8.4," (DimensionLESS)  *")') &
                          QPOINT
      Write(NFIPRI,'("*   NUMBER OF PROTONS IN NECK ................ Nz : ",F8.4," (DimensionLESS)  *")') &
                          protons_neck
      Write(NFIPRI,'("*   POSITION OF THE CENTERS OF MASS .........z_CM : ",2F8.4," [FERMIS] *")') &
                          CENLEF,CENRIG
      Write(NFIPRI,'("*   DISTANCE BETWEEN THE TWO FRAGMENTS........ D  : ",F8.4," [FERMIS]         *")') &
                          DISFRA
      Write(NFIPRI,'("*   SKYRME INTERACTION ENERGY................. ESk:",F11.3,"  [MEV]         *")') &
                          ESKINT
      Write(NFIPRI,'("*   TOTAL KINETIC ENERGY (COULOMB DIRECT) .... TKE:",F11.3,"  [MEV]         *")') &
                          EC_INT
      Write(NFIPRI,'("*   EXACT COULOMB ENERGY (DIRECT) ............ ECd:",F11.3,"  [MEV]         *")') &
                          ECDINT_ex
      Write(NFIPRI,'("*   EXACT COULOMB ENERGY (EXCHANGE) .......... ECx:",F11.3,"  [MEV]         *")') &
                          ECXINT_ex

      Write(NFIPRI,'("*",77X,"*")')
      Write(NFIPRI,'("*                       BREAK-DOWN OF SKYRME INTERACTION ENERGY               *")')
      Write(NFIPRI,'("*                 TOTAL(T)        SUM(S)        ISOSCALAR(P)  ISOVECTOR(M)    *")')
      Write(NFIPRI,'("*                 --------        ------        ------------  ------------    *")')
      Do i=1,NDCOUP/4
         Write(NFIPRI,'("*",3X,A6,4X,"|",2F14.6," |",2F14.6,5X,"*")') &
                        NCPALL(i),fis_ECPALL_12(2*i-1),fis_ECPALL_12(2*i),fis_EISALL_12(2*i-1),fis_EISALL_12(2*i)
      End Do
      Write(NFIPRI,'("*             |  ============  ============ |  ============  ============     *")')
      Write(NFIPRI,'("*   SUM EVEN  |",2F14.6," |",2F14.6,5X,"*")') &
                          sumEven_total_12,sumEven_sum_12,sumEven_scalar_12,sumEven_vector_12
      Write(NFIPRI,'("*",77X,"*")')
      Do i=NDCOUP/4+1,NDCOUP/2
         Write(NFIPRI,'("*",3X,A6,4X,"|",2F14.6," |",2F14.6,5X,"*")') &
                        NCPALL(i),fis_ECPALL_12(2*i-1),fis_ECPALL_12(2*i),fis_EISALL_12(2*i-1),fis_EISALL_12(2*i)
      End Do
      Write(NFIPRI,'("*             |  ============  ============ |  ============  ============     *")')
      Write(NFIPRI,'("*   SUM ODD   |",2F14.6," |",2F14.6,5X,"*")') &
                          sumOdd_total_12,sumOdd_sum_12,sumOdd_scalar_12,sumOdd_vector_12


      Write(NFIPRI,'("*",77X,"*")')
      Write(NFIPRI,'("*   ",71("-"),"   *")')
      Write(NFIPRI,'("*",13X,"|   LEFT  FRAGMENT (z < zN)   |   RIGHT FRAGMENT (z > zN)",7x,"*")')
      Write(NFIPRI,'("*   ",71("-"),"   *")')
      Write(NFIPRI,'("*",7X,"<A>   |",6X,F11.4,12X,"|",6X,F11.4,16X,"*")') A_MASS-RIGMAS,RIGMAS
      Write(NFIPRI,'("*",7X,"<Z>   |",6X,F11.4,12X,"|",6X,F11.4,16X,"*")') Real(IZ_FIX,Kind=pr)-RIGCHA,RIGCHA
      Write(NFIPRI,'("*",7X,"Etot  |",6X,F11.4,12X,"|",6X,F11.4,16X,"*")') ETOT_1,ETOT_2
      Write(NFIPRI,'("*",7X,"Ekin  |",6X,F11.4,12X,"|",6X,F11.4,16X,"*")') EKIN1N+EKIN1P,EKIN2N+EKIN2P
      Write(NFIPRI,'("*",7X,"Enuc  |",6X,F11.4,12X,"|",6X,F11.4,16X,"*")') ESKFR1,ESKFR2
      Write(NFIPRI,'("*",7X,"Epair |",6X,F11.4,12X,"|",6X,F11.4,16X,"*")') EPAIR1,EPAIR2
      Write(NFIPRI,'("*",7X,"EcouD |",6X,F11.4,12X,"|",6X,F11.4,16X,"*")') ECOUL1,ECOUL2
      Write(NFIPRI,'("*",7X,"EcouX |",6X,F11.4,12X,"|",6X,F11.4,16X,"*")') ECOEX1,ECOEX2
      Write(NFIPRI,'("*",7X,"<P^2> |",6X,F11.4,12X,"|",6X,F11.4,16X,"*")') TPSQU1,TPSQU2

      Write(NFIPRI,'("*   ",71("-"),"   *")')
      Write(NFIPRI,'("*                     BREAK-DOWN OF FRAGMENT SKYRME ENERGY                    *")')
      Write(NFIPRI,'("*               ISOSCALAR(P)  ISOVECTOR(M)  |  ISOSCALAR(P)  ISOVECTOR(M)     *")')
      Write(NFIPRI,'("*   ",71("-"),"   *")')
      Do i=1,NDCOUP/4
         Write(NFIPRI,'("*",3X,A6,4X,"|",2F14.6," |",2F14.6,5X,"*")') &
                        NCPALL(i),fis_EISALL_1(2*i-1),fis_EISALL_1(2*i),fis_EISALL_2(2*i-1),fis_EISALL_2(2*i)
      End Do
      Write(NFIPRI,'("*             |  ============  ============ |  ============  ============     *")')
      Write(NFIPRI,'("*   SUM EVEN  |",2F14.6," |",2F14.6,5X,"*")') &
                          sumEven_scalar_1,sumEven_vector_1,sumEven_scalar_2,sumEven_vector_2
      Write(NFIPRI,'("*   ",71("-"),"   *")')
      Do i=NDCOUP/4+1,NDCOUP/2
         Write(NFIPRI,'("*",3X,A6,4X,"|",2F14.6," |",2F14.6,5X,"*")') &
                        NCPALL(i),fis_EISALL_1(2*i-1),fis_EISALL_1(2*i),fis_EISALL_2(2*i-1),fis_EISALL_2(2*i)
      End Do
      Write(NFIPRI,'("*             |  ============  ============ |  ============  ============     *")')
      Write(NFIPRI,'("*   SUM ODD   |",2F14.6," |",2F14.6,5X,"*")') &
                          sumOdd_scalar_1,sumOdd_vector_1,sumOdd_scalar_2,sumOdd_vector_2

      Write(NFIPRI,'("*   ",71("-"),"   *")')
      Write(NFIPRI,'("*                 MULTIPOLE MOMENTS IN FRAGMENT INTRINSIC FRAME               *")')
      Write(NFIPRI,'("*   ",71("-"),"   *")')
      Do LAMBDA=0,NMUPRI
         Do MIU=0,LAMBDA
            If(MIU.Eq.0.Or.(LAMBDA.Eq.2.And.MIU.Eq.2)) Then
               Write(NFIPRI,'("*",6X,"<Q",2I1,">  |",6X,F11.4,12X,"|",6X,F11.4,16X,"*")') &
                           LAMBDA,MIU, ALLQLM(LAMBDA,MIU,0),                                &
                                       ALLQLM(LAMBDA,MIU,1)
            End If
         End Do
      End Do

      Write(NFIPRI,'("*   ",71("-"),"   *")')
      Write(NFIPRI,'("*             MULTIPOLE MOMENTS IN COMPOUND NUCLEUS INTRINSIC FRAME           *")')
      Write(NFIPRI,'("*   ",71("-"),"   *")')
      Do LAMBDA=0,NMUPRI
         Do MIU=0,LAMBDA
            If(MIU.Eq.0.Or.(LAMBDA.Eq.2.And.MIU.Eq.2)) Then
               Write(NFIPRI,'("*",6X,"<Q",2I1,">  |",6X,F11.4,12X,"|",6X,F11.4,16X,"*")') &
                           LAMBDA,MIU, ORIQLM(LAMBDA,MIU,0),                                &
                                       ORIQLM(LAMBDA,MIU,1)
            End If
         End Do
      End Do

      DUMMY=0.0_pr

      Write(NFIPRI,'("*   ",71("-"),"   *")')
      Write(NFIPRI,'("*   <dN^2>_QF |",8X,F9.4,12X,"|",8X,F9.4,16X,"*")') DUMMY,DUMMY
      Write(NFIPRI,'("*   <dN^2>_SF |",8X,F9.4,12X,"|",8X,F9.4,16X,"*")') DUMMY,DUMMY
      Write(NFIPRI,'("*",77X,"*",/,79("*"),/)')

    End Subroutine PRIFRA

 End Module hfodd_fission_print

