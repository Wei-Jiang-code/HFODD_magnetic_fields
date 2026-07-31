!***********************************************************************
!
!  Written by T. Lesinski 04/2013
!
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

module hfodd_pnp
  implicit none

!***********************************************************************

! No globals here

!***********************************************************************

contains

  !=====================================================================
  !
  !     This subroutine calculates the energy after projection on
  !     neutron and proton particle numbers. It uses the canonical
  !     representation (canonical orbitals taken from MODULE CANBAS,
  !     populated by SUBROUTINE CANQUA). Canonical quasiparticles are
  !     built from canonical orbitals and u- and v- factors then passed
  !     to various routines calculating densities and energies.
  !
  !     TODO:
  !      * Final sanity check on the input parameters
  !      * Optimize:
  !        + Do left wf definition once
  !        + Store n/p densities for each value of n/p gauge angle
  !      * PN-projected kernel calculation ?
  !
  !=====================================================================

  SUBROUTINE PNPROJ(IN_FIX,IZ_FIX,NPPNPN,NPPNPP,               &
       &            NXHERM,NYHERM,NZHERM,NXMAXX,NYMAXX,NZMAXX, &
       &            ISIMTX,JSIMTY,ISIMTZ,ISIGNY,ISIMPY,ISIQTY, &
       &                                        IPAHFB,MREVER, &
       &                                        MIN_QP,IDOTHC, &
       &            IPNMIX,NUMCOU,BOUCOU,ICOTYP,ICOUDI,ICOUEX, &
       &            COUSCA,NAMEPN,PRINIT,IDEVAR,ITERUN,        &
       &            ISYMDE,INIROT,INIINV,INIKAR,ISAWAV,EPNPRJ)
      USE CANBAS ! Canonical orbitals and occupations from CANQUA/Z
      USE WAVR_L ! Left-right wave functions to be passed to DENSHF
      USE HFODD_SIZES

      !-----------------------------------------------------------------

      IMPLICIT NONE

      !-----------------------------------------------------------------
      ! Input

      ! Number of Points for Particle-Number Projection of Neutrons/Protons
      INTEGER, INTENT(IN) :: NPPNPN,NPPNPP

      ! Standard parameters -- see subroutine NAMELI
      INTEGER, INTENT(IN) :: IN_FIX,IZ_FIX
      INTEGER, INTENT(IN) :: NXHERM,NYHERM,NZHERM
      INTEGER, INTENT(IN) :: NXMAXX,NYMAXX,NZMAXX
      INTEGER, INTENT(IN) :: ISIMTX,JSIMTY,ISIMTZ
      INTEGER, INTENT(IN) :: ISIGNY,ISIMPY,ISIQTY
      INTEGER, INTENT(IN) :: IPAHFB
      INTEGER, INTENT(IN) :: MREVER
      INTEGER, INTENT(IN) :: IDOTHC
      INTEGER, INTENT(IN) :: MIN_QP
      INTEGER, INTENT(IN) :: IPNMIX
      INTEGER, INTENT(IN) :: NUMCOU,ICOTYP,ICOUDI,ICOUEX
      REAL,    INTENT(IN) :: BOUCOU,COUSCA
      LOGICAL, INTENT(IN) :: PRINIT
      INTEGER, INTENT(IN) :: IDEVAR
      INTEGER, INTENT(IN) :: ITERUN
      INTEGER, INTENT(IN) :: ISYMDE
      INTEGER, INTENT(IN) :: INIROT
      INTEGER, INTENT(IN) :: INIINV
      INTEGER, INTENT(IN) :: INIKAR
      INTEGER, INTENT(IN) :: ISAWAV

      CHARACTER(LEN=8), INTENT(IN) :: NAMEPN

      !-----------------------------------------------------------------
      ! Output: projected energy

      REAL, INTENT(INOUT) :: EPNPRJ

      !=================================================================
      ! Globals

      INTEGER :: LDBASE, LDTOTS, LDSTAT, LDUPPE, LDTIMU, NUMBQP
      COMPLEX :: PHASPI, ECCALL

      COMMON &
           /DIMENS/ LDBASE
      COMMON &
           /DIMSTA/ LDTOTS(0:NDISOS),LDSTAT(0:NDISOS), &
           &        LDUPPE(0:NDISOS),LDTIMU(0:NDISOS)
      COMMON &
           /SPIPHA/ PHASPI(1:NDBASE,0:NDREVE,0:NDSPIN)
      COMMON &
           /QPNUMS/ NUMBQP(0:NDREVE,0:NDISOS)
      COMMON &
           /ECCSKY/ ECCALL(NDCOUP)

     !=================================================================

     REAL, PARAMETER :: PI = 4.0D0*ATAN(1.0D0)

     !=================================================================

     INTEGER :: IALLOC, IBASE,  ISTATE, JSTATE, ISPIN,  ICHARG, IN, IP
     INTEGER :: IREVER, INN, IZZ

     REAL    :: THETAN, THETAP, ENEDUM, EKESCA, EKEVEC, XITFAC
     REAL    :: PARTICLENUMBER

     COMPLEX :: ZSHIFN, ZSHIFP, ZSHIFT, CFUPPR, CFLOWR
     COMPLEX :: OVRLAP, EKERNL, EDUMCX
     COMPLEX :: EKEKIN(0:NDISOS), EKESKY, &
          EKESVO, EKESSU, EKESSO, EKESTE, &
          EKECOU, EKECOD, EKECOE, EKEPAI(0:NDISOS)

     COMPLEX, DIMENSION(IN_FIX,IZ_FIX) ::                   &
          EKNITG,EKPITG,ESKITG,EVOITG,ESUITG,ESOITG,ETEITG, &
          ECOITG,ECDITG,ECEITG,                             &
          EPNITG,EPPITG,ENEITG,OVRITG

     COMPLEX, DIMENSION(NDCOUP,IN_FIX,IZ_FIX) :: ECCITG

     !=================================================================
     ! Time execution
     CALL CPUTIM('PNPROJ',1)


     PRINT '("*",77("*"),"*")'
     PRINT '("*",77x,"*")'
     PRINT '("*",24x,"PROJECTING ON PARTICLE NUMBER",24x,"*")'
     PRINT '("*",77x,"*")'
     PRINT '("*",77("*"),"*")'

     !=================================================================
     ! Prepare arrays for storing shifted wave functions

     IF (.NOT.ALLOCATED(WALEFT)) THEN
        ALLOCATE (WALEFT(1:NDBASE,1:4*NDSTAT,0:NDSPIN),STAT=IALLOC)
        IF (IALLOC.NE.0) CALL NOALLO('WALEFT','PNPROJ')
     END IF
     WALEFT(:,:,:) = 0.0d0

     IF (.NOT.ALLOCATED(WARIGH)) THEN
        ALLOCATE (WARIGH(1:NDBASE,1:4*NDSTAT,0:NDSPIN),STAT=IALLOC)
        IF (IALLOC.NE.0) CALL NOALLO('WARIGH','PNPROJ')
     END IF
     WARIGH(:,:,:) = 0.0d0

     !=================================================================

     ! TODO: Final sanity check on the parameters

     !=================================================================
     ! Get canonical WFs, v^2, uv
     ! ==> We assume CANQUA was called, populating the CANBAS module.

     !-----------------------------------------------------------------
     ! Init arrays for accumulating integrands

      ENEITG(:,:)=0.0d0  ! Total energy
      OVRITG(:,:)=0.0d0  ! Overlap

      EKNITG(:,:)=0.0d0  ! Kinetic n
      EKPITG(:,:)=0.0d0  ! Kinetic p
      ESKITG(:,:)=0.0d0  ! Skyrme total
      ECOITG(:,:)=0.0d0  ! Coulomb total
      ECDITG(:,:)=0.0d0  ! Coulomb direct
      ECEITG(:,:)=0.0d0  ! Coulomb exchange
      EPNITG(:,:)=0.0d0  ! Pairing n
      EPPITG(:,:)=0.0d0  ! Pairing p

      ECCITG(:,:,:)=0.0d0 ! Skyrme, term-by-term

      ! Integral prefactor
      XITFAC = 1.0d0/(NPPNPN*NPPNPP)

      !-----------------------------------------------------------------
      ! Begin loop on gauge-angle integration points

      DO IN = 1, NPPNPN
         DO IP = 1, NPPNPP

            ! Init energy kernel. We'll add up contributions there.
            EKERNL=0.0d0

            ! Gauge angles
            THETAN= (IN-1.0d0) * PI / NPPNPN
            THETAP= (IP-1.0d0) * PI / NPPNPP

            ! Complex shifts
            ZSHIFN=CMPLX(COS(THETAN), SIN(THETAN))
            ZSHIFP=CMPLX(COS(THETAP), SIN(THETAP))

            !...........................................................
            ! Rebuild quasiparticles with U(1) rotation built-in
            ! (we convert from simplex to spin representation on the fly)
            ! then calculate transition densities

            DO ICHARG=0,NDISOS

               ! Shift (z)
               IF ( ICHARG == 0 ) THEN
                  ! Neutrons
                  ZSHIFT = ZSHIFN
               ELSE
                  ! Protons
                  ZSHIFT = ZSHIFP
               END IF

               JSTATE=0 ! Counter for wfs we put in WALEFT/WARIGH

               ! . . . . . . . . . . . . . . . . . . . . . . . . . . . .

               ! First, upper components time+ then time-
               DO IREVER=0,1

                  ! NB: We iterate on both time-reversal blocks since we
                  ! need both for non-diagonal energy kernels (V_n and
                  ! V_(-n) cannot be recovered fom each other for the
                  ! right state because of the additional complex shift
                  ! factor).

                  ! Loop on canonical states
                  DO ISTATE=1,NUMBQP(IREVER,ICHARG)

                     JSTATE=JSTATE+1

                     ! Left wave function
                     DO ISPIN=0,1
                        WALEFT(1:LDBASE,JSTATE,ISPIN) =              &
                             V_CAN(ISTATE,IREVER,ICHARG)             &
                             * WAVCAN(1:LDBASE,ISTATE,IREVER,ICHARG) &
                             * PHASPI(1:LDBASE,IREVER,ISPIN)
                     END DO

                     ! Right (shifted and normalized) wave function...
                     CFUPPR = ZSHIFT**2 * V_CAN(ISTATE,IREVER,ICHARG)  &
                          / ( U_CAN(ISTATE,IREVER,ICHARG)**2           &
                          +   (ZSHIFT*V_CAN(ISTATE,IREVER,ICHARG))**2 )

                     DO ISPIN=0,1
                        WARIGH(1:LDBASE,JSTATE,ISPIN) = CFUPPR         &
                             * WAVCAN(1:LDBASE,ISTATE,IREVER,ICHARG)   &
                             * PHASPI(1:LDBASE,IREVER,ISPIN)
                     END DO
                  END DO ! ISTATE
                  ! Number of time-up states
                  IF (IREVER==0) LDTIMU(ICHARG) = JSTATE

               END DO ! IREVER
               ! Number of upper HFB wave functions
               LDUPPE(ICHARG) = JSTATE

               ! . . . . . . . . . . . . . . . . . . . . . . . . . . . .

               ! Then, lower components, time- then time+ We
               ! apply time-reversal on the wfs like in DENSHF to
               ! follow the Dobaczewski-Flocard-Treiner
               ! convention
               DO IREVER=1,0,-1

                  DO ISTATE=1,NUMBQP(IREVER,ICHARG)

                     JSTATE=JSTATE+1

                     ! Left wave function
                     DO ISPIN=0,1
                        WALEFT(1:LDBASE,JSTATE,ISPIN) =                &
                             U_CAN(ISTATE,IREVER,ICHARG) * (2*ISPIN-1) &
                             * CONJG(                                  &
                             WAVCAN(1:LDBASE,ISTATE,IREVER,ICHARG)     &
                             * PHASPI(1:LDBASE,IREVER,1-ISPIN) )
                     END DO

!                    Right (normalized) wave function

                     CFLOWR = U_CAN(ISTATE,IREVER,ICHARG)              &
                          / ( U_CAN(ISTATE,IREVER,ICHARG)**2           &
                          +   (ZSHIFT*V_CAN(ISTATE,IREVER,ICHARG))**2 )

                     DO ISPIN=0,1
                        WARIGH(1:LDBASE,JSTATE,ISPIN) =           &
                            CFLOWR * (2*ISPIN-1)                  &
                            * CONJG(                              &
                            WAVCAN(1:LDBASE,ISTATE,IREVER,ICHARG) &
                            * PHASPI(1:LDBASE,IREVER,1-ISPIN) )
                     END DO

                  END DO        ! ISTATE

                  ! Number of states up to now....
                  IF(IREVER.EQ.1) LDSTAT(ICHARG) = JSTATE

               END DO           ! IREVER
               ! Total number of states
               LDTOTS(ICHARG) = JSTATE

               ! . . . . . . . . . . . . . . . . . . . . . . . . . . . .

               ! Calculating transition density matrix for K.E. computation
               ! (setting MREVER=1, see NB: above, and IPAHFB=1)
               CALL DENMAC(1,ICHARG,ISIMPY,1,WALEFT,WARIGH)
               ! Calculating kinetic energy kernel for current isospin
               CALL EKINET(ENEDUM,ENEDUM,EKEKIN(ICHARG),EDUMCX)

               ! . . . . . . . . . . . . . . . . . . . . . . . . . . . .

               ! Calculating transition densities for current isospin
               ! (setting IPAHFB=1, MREVER=1,
               ! ITPNMX=ICHARG, ITIREP=0, IKERNE=1)
               CALL ZEDENS(ICHARG) ! Zero densities first
               CALL DENSHF(NXHERM,NYHERM,NZHERM,ISIMTX,JSIMTY,ISIMTZ, &
                    &      ISIGNY,ISIMPY,ISIQTY,1     ,1     ,ICHARG, &
                    &                           MIN_QP,IPNMIX,ICHARG, &
                    &      0     ,NXMAXX,NAMEPN,PRINIT,IDEVAR,ITERUN, &
                    &      ISYMDE,INIROT,INIINV,INIKAR,ISAWAV,1     )

            END DO ! ICHARG

            !...........................................................
            ! Calculate Coulomb energy kernel with previously calculated
            ! densities

            ! Default value if no Coulomb
            EKECOD = 0.0d0
            EKECOE = 0.0d0

            ! Calculate charge density
            CALL TRUCHD(NXHERM,NYHERM,NZHERM)

            ! Direct Coulomb field and energy by solution of Poisson's
            ! equation
            IF (ICOUDI.EQ.1) THEN
               ! (setting IKERNE=1)
               CALL COUMAT(NXHERM,NYHERM,NZHERM,NXMAXX,NYMAXX,NZMAXX, &
               &           NUMCOU,BOUCOU,ISIMPY,1)
               CALL COULOD(ISIMPY,EKECOD)
            END IF

            ! Slater approximation for Coulomb exchange energy
            IF (ICOUEX.EQ.1) THEN
               CALL COULOE(NXHERM,NYHERM,NZHERM,EKECOE)
            END IF

            ! NB: After exiting the loop on ICHARG, BIG_PP contains the
            ! proton density matrix.
            ! ***
            ! WE RELY ON THAT IN THE CALL TO COUENE !!!
            ! ***

            ! Direct and exact exchange through matrix elements
            ! (Gaussian expansion)
            IF (ICOUDI.EQ.2.OR.ICOUEX.EQ.2) THEN
               CALL COUENE(NXHERM,NYHERM,NZHERM,NXMAXX,NYMAXX,NZMAXX, &
                    &      ISIMPY,               &
                    &      ICOTYP,ICOUDI,ICOUEX, &
                    &      EKECOD,EKESCA,EKEVEC)
               !  Replace previous value if appropriate
               IF (ICOUEX.EQ.2) EKECOE=EKESCA+EKEVEC
            END IF

            ! To match HFBTHO
            EKECOD = CONJG(EKECOD)
            EKECOE = CONJG(EKECOE)

            ! Direct+exchange coulomb energy
            EKECOU=EKECOD + EKECOE*COUSCA

            !...........................................................
            ! Calculate Skyrme energy kernel with previously calculated
            ! densities

            ! Calculate density for density-dependent terms
            CALL TRUTOD(NXHERM,NYHERM,NZHERM)

            ! Skyrme energy (setting LDPNMX=1)
            CALL ESKYRM(NXHERM,NYHERM,NZHERM,        &
                 &      ENEDUM,ENEDUM,ENEDUM,        &
                 &      ENEDUM,ENEDUM,ENEDUM,ENEDUM, &
                 &      EKESKY,1)
            ! (Components of Skyrme energy now stored in ECCALL)

            ! Same for Pairing energy.
            DO ICHARG=0,NDISOS
               ! (setting ITPNMX=ICHARG)
               CALL EPAIRI(NXHERM,NYHERM,NZHERM,IN_FIX,IZ_FIX,ICHARG, &
                    &      EKEPAI(ICHARG))
            END DO

            ! To match HFBTHO
            EKESKY    = CONJG(EKESKY)
            ECCALL(:) = CONJG(ECCALL(:))
            EKEPAI(:) = CONJG(EKEPAI(:))

            !...........................................................
            ! Compute overlap between left and right wave functions
            ! (Onishi formula for canonical basis)

            OVRLAP=1.0d0
            DO ICHARG=0,NDISOS

               ! Shift (z)
               IF ( ICHARG .EQ. 0 ) THEN
                  ! Neutrons
                  ZSHIFT = ZSHIFN
               ELSE
                  ! Protons
                  ZSHIFT = ZSHIFP
               END IF

               DO ISTATE=1, NUMBQP(0,ICHARG)
                  OVRLAP=OVRLAP                      &
                       * ( U_CAN(ISTATE,0,ICHARG)**2 &
                       +   (ZSHIFT*V_CAN(ISTATE,0,ICHARG))**2 )
               END DO

            END DO

            !...........................................................
            ! Accumulate contributions

            WRITE(75, '(6es16.8)') THETAN, THETAP, OVRLAP, EKERNL

            EKERNL=EKEKIN(0)+EKEKIN(1)+EKESKY+EKECOU &
                 &    + EKEPAI(0)+EKEPAI(1)

            DO INN=1, IN_FIX
               DO IZZ=1, IZ_FIX
                  ! Components of the energy-integral

                  ! Energy-integral by coupling constants
                  ECCITG(:,INN,IZZ)=ECCITG(:,INN,IZZ)       &
                      + ZSHIFN**(-2*INN) * ZSHIFP**(-2*IZZ) &
                      * XITFAC*OVRLAP*ECCALL(:)
                  ! Skyrme-energy integral
                  ESKITG(INN,IZZ)=ESKITG(INN,IZZ)            &
                       + ZSHIFN**(-2*INN) * ZSHIFP**(-2*IZZ) &
                       * XITFAC*OVRLAP*EKESKY

                  ! Kinetic-energy integrals
                  EKNITG(INN,IZZ)=EKNITG(INN,IZZ)           &
                       + ZSHIFN**(-2*INN) * ZSHIFP**(-2*IZZ) &
                       * XITFAC*OVRLAP*EKEKIN(0)
                  EKPITG(INN,IZZ)=EKPITG(INN,IZZ)           &
                       + ZSHIFN**(-2*INN) * ZSHIFP**(-2*IZZ) &
                       * XITFAC*OVRLAP*EKEKIN(1)

                  ! Coulomb-energy integrals
                  ECOITG(INN,IZZ)=ECOITG(INN,IZZ)            &
                       + ZSHIFN**(-2*INN) * ZSHIFP**(-2*IZZ) &
                       * XITFAC*OVRLAP*EKECOU
                  ECDITG(INN,IZZ)=ECDITG(INN,IZZ)            &
                       + ZSHIFN**(-2*INN) * ZSHIFP**(-2*IZZ) &
                       * XITFAC*OVRLAP*EKECOD
                  ECEITG(INN,IZZ)=ECEITG(INN,IZZ)            &
                       + ZSHIFN**(-2*INN) * ZSHIFP**(-2*IZZ) &
                       * XITFAC*OVRLAP*EKECOE

                  ! Pairing-energy integrals
                  EPNITG(INN,IZZ)=EPNITG(INN,IZZ)            &
                       + ZSHIFN**(-2*INN) * ZSHIFP**(-2*IZZ) &
                       * XITFAC*OVRLAP*EKEPAI(0)
                  EPPITG(INN,IZZ)=EPPITG(INN,IZZ)            &
                       + ZSHIFN**(-2*INN) * ZSHIFP**(-2*IZZ) &
                       * XITFAC*OVRLAP*EKEPAI(1)

                  ! Energy-integral
                  ENEITG(INN,IZZ)=ENEITG(INN,IZZ)            &
                       + ZSHIFN**(-2*INN) * ZSHIFP**(-2*IZZ) &
                      * XITFAC*OVRLAP*EKERNL
                  ! Overlap-integral
                  OVRITG(INN,IZZ)=OVRITG(INN,IZZ)            &
                       + ZSHIFN**(-2*INN) * ZSHIFP**(-2*IZZ) &
                       * XITFAC*OVRLAP
               END DO
            END DO

            !...........................................................
            ! Print summary

            PRINT '("*",77x,"*")'
            PRINT '("*  z_n        ",2f12.5,"   z_p        ",2f12.5,"  *")', &
                 ZSHIFN,ZSHIFP
            PRINT '("*  Skyrme   E ",2f12.5,"   n kin. E   ",2f12.5,"  *")', &
                 EKESKY,EKEKIN(0)
            PRINT '("*   volume    ",2f12.5,"   p kin. E   ",2f12.5,"  *")', &
                   ECCALL(1) + ECCALL(2) + ECCALL(3) + ECCALL(4)      &
                 + ECCALL(7) + ECCALL(8), EKEKIN(1)
            PRINT '("*    contact  ",2f12.5,"   Coulomb E  ",2f12.5,"  *")', &
                 ECCALL(1) + ECCALL(2), EKECOU
            PRINT '("*    dens-dep ",2f12.5,"    direct    ",2f12.5,"  *")', &
                 ECCALL(3) + ECCALL(4), EKECOD
            PRINT '("*    eff mass ",2f12.5,"    exchange  ",2f12.5,"  *")', &
                 ECCALL(7) + ECCALL(8), EKECOE
            PRINT '("*   surface   ",2f12.5,"   n pair. E  ",2f12.5,"  *")', &
                 ECCALL(5) + ECCALL(6), EKEPAI(0)
            PRINT '("*   spin-o.   ",2f12.5,"   p pair. E  ",2f12.5,"  *")', &
                 ECCALL(11) + ECCALL(12), EKEPAI(1)
            PRINT '("*   tensor    ",2f12.5,38(" "),"  *")',                          &
                  ECCALL(9) + ECCALL(10)

            PRINT '("* Overlap     ",2f12.5,"  E kernel    ",2f12.5,"  *")',  &
                OVRLAP,   EKERNL

         END DO ! IP
      END DO ! IN

      ! (end loop on gauge-angle integration points)
      !-----------------------------------------------------------------

      ! Projected energy = ratio of integrals

      EPNPRJ = REAL(ENEITG(IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2))

      !     DEBUG DEBUG DEBUG
      !      WRITE(76, '(2a4,8(2x,2a15))')
      !     *           "# N", "Z",
      !     *           "ENEITG", "OVRITG",
      !     *           "EKNITG", "EKPITG",
      !     *           "ESKITG", "ECOITG",
      !     *           "EPNITG", "EPPITG"
      !      DO INN=1,IN_FIX
      !         DO IZZ=1,IZ_FIX
      !            WRITE(76, '(2i4,8(2x,2es15.6))')
      !     *           2*INN, 2*IZZ,
      !     *           ENEITG(INN,IZZ), OVRITG(INN,IZZ),
      !     *           EKNITG(INN,IZZ), EKPITG(INN,IZZ),
      !     *           ESKITG(INN,IZZ), ECOITG(INN,IZZ),
      !     *           EPNITG(INN,IZZ), EPPITG(INN,IZZ)
      !         END DO
      !         WRITE(76, *)
      !      END DO

      !------------------------------------------------------------------
      ! Print final results

      !PRINT '(79("*"))'
      !PRINT *, "* Projected term-by-term Skyrme energies"
      !PRINT '("*  1 ERHO_T = ",(2f13.6))', &
      !    ECCITG( 1,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("*  2 ERHO_S = ",(2f13.6))', &
      !    ECCITG( 2,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("*  3 ERHODT = ",(2f13.6))', &
      !    ECCITG( 3,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("*  4 ERHODS = ",(2f13.6))', &
      !    ECCITG( 4,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("*  5 ELPR_T = ",(2f13.6))', &
      !    ECCITG( 5,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("*  6 ELPR_S = ",(2f13.6))', &
      !    ECCITG( 6,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("*  7 ETAU_T = ",(2f13.6))', &
      !    ECCITG( 7,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("*  8 ETAU_S = ",(2f13.6))', &
      !    ECCITG( 8,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("*  9 ESCU_T = ",(2f13.6))', &
      !    ECCITG( 9,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("* 10 ESCU_S = ",(2f13.6))', &
      !    ECCITG(10,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("* 11 EDIV_T = ",(2f13.6))', &
      !    ECCITG(11,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("* 12 EDIV_S = ",(2f13.6))', &
      !    ECCITG(12,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("* 13 ESPI_T = ",(2f13.6))', &
      !    ECCITG(13,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("* 14 ESPI_S = ",(2f13.6))', &
      !    ECCITG(14,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("* 15 ESPIDT = ",(2f13.6))', &
      !    ECCITG(15,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("* 16 ESPIDS = ",(2f13.6))', &
      !    ECCITG(16,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("* 17 ELPS_T = ",(2f13.6))', &
      !    ECCITG(17,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("* 18 ELPS_S = ",(2f13.6))', &
      !    ECCITG(18,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("* 19 ECUR_T = ",(2f13.6))', &
      !    ECCITG(19,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("* 20 ECUR_S = ",(2f13.6))', &
      !    ECCITG(20,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("* 21 EKIS_T = ",(2f13.6))', &
      !    ECCITG(21,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("* 22 EKIS_S = ",(2f13.6))', &
      !    ECCITG(22,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("* 23 EROT_T = ",(2f13.6))', &
      !    ECCITG(23,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)
      !PRINT '("* 24 EROT_S = ",(2f13.6))', &
      !    ECCITG(24,IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/y2)

      PRINT '("*",77x,"*")'
      PRINT '("*",77("*"),"*")'
      PRINT '("*",77x,"*")'
      PRINT '("* E Skyrme     ",f15.6,"  E kinetic, n ",f15.6,18x,"*")', &
          REAL(ESKITG(IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2)), &
           REAL(EKNITG(IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2))

      PRINT '("*  volume      ",f15.6,"  E kinetic, p ",f15.6,18x,"*")', &
           REAL(                            &
           ( ECCITG( 1,IN_FIX/2,IZ_FIX/2) + ECCITG( 2,IN_FIX/2,IZ_FIX/2)  &
           + ECCITG( 3,IN_FIX/2,IZ_FIX/2) + ECCITG( 4,IN_FIX/2,IZ_FIX/2)  &
           + ECCITG( 7,IN_FIX/2,IZ_FIX/2) + ECCITG( 8,IN_FIX/2,IZ_FIX/2)) &
           /OVRITG(IN_FIX/2,IZ_FIX/2)),                                   &
           REAL(EKPITG(IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2))

      PRINT '("*   contact    ",f15.6,"  E Coulomb    ",f15.6,18x,"*")', &
           REAL(                            &
          ( ECCITG( 1,IN_FIX/2,IZ_FIX/2) + ECCITG( 2,IN_FIX/2,IZ_FIX/2))  &
          /OVRITG(IN_FIX/2,IZ_FIX/2)),                                    &
           REAL(ECOITG(IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2))

      PRINT '("*   dens. dep. ",f15.6,"    direct     ",f15.6,18x,"*")', &
           REAL(                            &
          ( ECCITG( 3,IN_FIX/2,IZ_FIX/2) + ECCITG( 4,IN_FIX/2,IZ_FIX/2))  &
          /OVRITG(IN_FIX/2,IZ_FIX/2)),                                    &
           REAL(ECDITG(IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2))

      PRINT '("*  eff. mass   ",f15.6,"    exchange   ",f15.6,18x,"*")', &
           REAL(                            &
          ( ECCITG( 7,IN_FIX/2,IZ_FIX/2) + ECCITG( 8,IN_FIX/2,IZ_FIX/2))  &
          /OVRITG(IN_FIX/2,IZ_FIX/2)),                                   &
           REAL(ECEITG(IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2))

      PRINT '("*  surface     ",f15.6,"  E pairing, n ",f15.6,18x,"*")', &
           REAL(                             &
           ( ECCITG( 5,IN_FIX/2,IZ_FIX/2) + ECCITG( 6,IN_FIX/2,IZ_FIX/2))  &
           /OVRITG(IN_FIX/2,IZ_FIX/2)),                                    &
           REAL(EPNITG(IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2))

      PRINT '("*  spin-o.     ",f15.6,"  E pairing, p ",f15.6,18x,"*")', &
           REAL(                             &
           ( ECCITG(11,IN_FIX/2,IZ_FIX/2) + ECCITG(12,IN_FIX/2,IZ_FIX/2))  &
           /OVRITG(IN_FIX/2,IZ_FIX/2)),                                   &
           REAL(EPPITG(IN_FIX/2,IZ_FIX/2)/OVRITG(IN_FIX/2,IZ_FIX/2))

      PRINT '("*  tensor      ",f15.6,48x,"*")', REAL(                     &
           ( ECCITG( 9,IN_FIX/2,IZ_FIX/2) + ECCITG(10,IN_FIX/2,IZ_FIX/2))  &
           /OVRITG(IN_FIX/2,IZ_FIX/2))
      PRINT '("*",77x,"*")'

      !------------------------------------------------------------------
      ! Norm of mean-field state from overlaps
      PRINT '("* PNP ENERGY                  ",f15.6,33x,"*")', EPNPRJ
      PRINT '("* N overlap                   ",2f15.6,18x,"*")', &
           OVRITG(IN_FIX/2,IZ_FIX/2)
      PRINT '("* Sum PNP overlap components  ",2f15.6,18x,"*")',      &
           SUM(OVRITG((IN_FIX-NPPNPN)/2+1:(IN_FIX-NPPNPN)/2+NPPNPN,   &
           (IZ_FIX-NPPNPP)/2+1:(IZ_FIX-NPPNPP)/2+NPPNPP))
      ! PNP sum rule
      PRINT '("* PNP sum Rule (Skyrme E)     ",2f15.6,18x,"*")',      &
           SUM(ESKITG((IN_FIX-NPPNPN)/2+1:(IN_FIX-NPPNPN)/2+NPPNPN,   &
           (IZ_FIX-NPPNPP)/2+1:(IZ_FIX-NPPNPP)/2+NPPNPP))
      PRINT '("* PNP sum Rule (total E)      ",2f15.6,18x,"*")',     &
           SUM(ENEITG((IN_FIX-NPPNPN)/2+1:(IN_FIX-NPPNPN)/2+NPPNPN,   &
           (IZ_FIX-NPPNPP)/2+1:(IZ_FIX-NPPNPP)/2+NPPNPP))

      PRINT '("*",77x,"*")'
      PRINT '("*",77("*"),"*")'

      !------------------------------------------------------------------
      ! Time execution
      CALL CPUTIM('PNPROJ',0)

    END SUBROUTINE PNPROJ

    !====================================================================

end module hfodd_pnp

!************************************************************************
