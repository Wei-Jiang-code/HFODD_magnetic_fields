C***********************************************************************
C
C    Copyright (c) 2016, Lawrence Livermore National Security, LLC.
c                        Produced at the Lawrence Livermore National
C                        Laboratory.
C                        Written by Nicolas Schunck, schunck1@llnl.gov
C
C    LLNL-CODE-710577 All rights reserved.
C    LLNL-CODE-470611 All rights reserved.
C
C    Copyright 2016, N. Schunck, J. Dobaczewski, W. Satula, P. Baczyk,
C                    J. Dudek, Y. Gao, M. Konieczka, K. Sato, Y. Shi,
C                    X.B. Wang, T.R. Werner
C    Copyright 2012, N. Schunck, J. Dobaczewski, J. McDonnell,
C                    W. Satula, J.A. Sheikh, A. Staszczak,
C                    M. Stoitsov, P. Toivanen
C    Copyright 2009, J. Dobaczewski, W. Satula, B.G. Carlsson, J. Engel,
C                    P. Olbratowski, P. Powalowski, M. Sadziak,
C                    J. Sarich, N. Schunck, A. Staszczak, M. Stoitsov,
C                    M. Zalewski, H. Zdunczuk
C    Copyright 2004, 2005, J. Dobaczewski, P. Olbratowski
C    Copyright 1997, 2000, J. Dobaczewski, J. Dudek
C
C    This file is part of HFODD.
C
C    HFODD is free software: you can redistribute it and/or modify it
C    under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    HFODD is distributed in the hope that it will be useful, but
C    WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with HFODD. If not, see <http://www.gnu.org/licenses/>.
C
C    OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
C    LICENSE
C
C    Our Preamble Notice
C
C      A. This notice is required to be provided under our contract
C         with the U.S. Department of Energy (DOE). This work was
c         produced at the Lawrence Livermore National Laboratory under
c         Contract No. DE-AC52-07NA27344 with the DOE.
C      B. Neither the United States Government nor Lawrence Livermore
C         National Security, LLC nor any of their employees, makes any
c         warranty, express or implied, or assumes any liability or
c         responsibility for the accuracy, completeness, or usefulness
c         of any information, apparatus, product, or process disclosed,
C         or represents that its use would not infringe privately-owned
c         rights.
C      C. Also, reference herein to any specific commercial products,
C         process, or services by trade name, trademark, manufacturer
C         or otherwise does not necessarily constitute or imply its
C         endorsement, recommendation, or favoring by the United States
C         Government or Lawrence Livermore National Security, LLC. The
c         views and opinions of authors expressed herein do not
C         necessarily state or reflect those of the United States
c         Government or Lawrence Livermore National Security, LLC, and
c         shall not be used for advertising or product endorsement
c         purposes.
C
C    The precise terms and conditions for copying, distribution and
C    modification are contained in the file COPYING.
C
C***********************************************************************

C----------------------------------------------------------------------!
C                                                                      !
C                Shell correction code for HFODD                       !
C                                                                      !
C Authors: N. Schunck, ORNL                                            !
C                                                                      !
C This module computes the shell correction based on the s.p. spectrum !
C generated in HFODD.                                                  !
C                                                                      !
C                          First included in official release v249h    !
C                                                                      !
C----------------------------------------------------------------------!

C=======================================================================
C     WSHELL  WSHELL  WSHELL  WSHELL   WSHELL  WSHELL  WSHELL  WSHELL
C=======================================================================
C                                    N. Schunck, ORNL,       2008-2010
C                                    J. Dudek,   Strasbourg,      2000
C=======================================================================
C
      SUBROUTINE HFINPU(IFSHEL,MREVER,IPAHFB,IMFHFB,NDACTU,ESINPP,
     *                  ESINPN,EGAS_P,EGAS_N,OCCU_P,OCCU_N,EFERMP,
     *                                                     EFERMN)
C=======================================================================
      USE hfodd_sizes
C=======================================================================
      DIMENSION
     *          NSTATE(1:2*NDACTU)
      DIMENSION
     *          ESINPP(1:2*NDACTU),ESINPN(1:2*NDACTU),
     *          EGAS_P(1:2*NDACTU),EGAS_N(1:2*NDACTU),
     *          OCCU_P(1:2*NDACTU),OCCU_N(1:2*NDACTU)
C
C@@@ TEMP - TEMP - TEMP - TEMP - TEMP - TEMP - TEMP - TEMP - TEMP
      COMMON
     *       /FER_IS/ FERMFN(1:NDBASE,0:NDISOS)
      COMMON
     *       /T_FLAG/ IFTEMP
      COMMON
     *       /GCTEMP/ TEMP_T
C@@@ TEMP - TEMP - TEMP - TEMP - TEMP - TEMP - TEMP - TEMP - TEMP
      COMMON
     *       /EQUISO/ ESPEQU(1:NDBASE,0:NDREVE,0:NDISOS),
     *                DELEQU(1:NDBASE,0:NDREVE,0:NDISOS)
      COMMON
     *       /SPLEVS/ SPENER(1:2*NDSTAT,0:NDISOS),
     *                SPOCCU(1:2*NDSTAT,0:NDISOS),
     *                SPALIG(1:2*NDSTAT,0:NDISOS),
     *                SPSPAL(1:2*NDSTAT,0:NDISOS),
     *                SPPARI(1:2*NDSTAT,0:NDISOS),
     *                SPSIMP(1:2*NDSTAT,0:NDISOS)
      COMMON
     *       /KIN_SP/ SPGASP(1:2*NDBASE),SPGASN(1:2*NDBASE)
      COMMON
     *       /QPNUMS/ NUMBQP(0:NDREVE,0:NDISOS)
      COMMON
     *       /SPNUMS/ NUMBSP(0:NDREVE,0:NDISOS)
C
      COMMON
     *       /WSHELL_KRAMER/ KRAMER
      COMMON
     *       /WSHELL_BASSIZ/ LDBASN,LDBASP
C
C=======================================================================
C
C      Interface between HFODD and the shell correction routines.
C      Single-particle energies
C
C      HFODD inputs
C         MREVER ..........: Flag for time-reversal symmetry
C         IPAHFB ..........: Flag for HFB versus HF calculations
C         IFTEMP ..........: Flag for temperature
C         NDACTU ..........: Static size of the s.p. arrays (should be NDBASE)
C         SPENER ..........: HF single-particle energies or HFB canonical
C                            spectrum
C         ESPEQU ..........: HFB equivalent spectrum
C         SPGASP, SPGASN ..: Eigenvalues of the Fermi free gas (proton
C                            and neutrons differ by the contribution
C                            from Coulomb)
C
C      Outputs
C         ESINPP, ESINPN ..: Protons and neutrons single-particle
C                            energies
C         EGAS_P, EGAS_N ..: Eigenvalues of the free gas
C         OCCU_P, OCCU_N ..: Statistical occupation of s.p. states
C         KRAMER ..........: Equivalent of MREVER, contains info on
C                            time-reversal symmetry
C
C      Conventions
C       - MREVER=0: time-reversal symmetry conserved, Kramers degeneracy
C                   applies, only first half of all HF s.p. needed
C       - MREVER=1: time-reversal symmetry broken, no Kramers degeneracy
C     	            full eigenvector is needed
C
C=======================================================================
C
      IF (NDBASE.NE.NDACTU) THEN
          WRITE(6,'(''NDBASE = '',I6,'' NDACTU = '',I6)')
     *                NDBASE,NDACTU
          STOP 'NDBASE .NE. NDACTU in HFINPU - WSHELL'
      END IF
C
C      KRAMER=1 for conserved time-reversal symmetry (MREVER=0)
C      KRAMER=0 for  broken   time-reversal symmetry (MREVER=1)
      KRAMER=1-MREVER
C
C      Sizes based on HF spectrum or HFB canonical spectrum
      IF (IPAHFB.EQ.0.OR.(IPAHFB.EQ.1.AND.IMFHFB.EQ.1)) THEN
      LDBASN=NUMBSP(0,0)+MREVER*NUMBSP(1,0)
      LDBASP=NUMBSP(0,1)+MREVER*NUMBSP(1,1)
      END IF
C
C      Sizes based on HFB equivalent spectrum
      IF (IPAHFB.EQ.1.AND.IMFHFB.EQ.0) THEN
          LDBASN=NUMBQP(0,0)+MREVER*NUMBQP(1,0)
          LDBASP=NUMBQP(0,1)+MREVER*NUMBQP(1,1)
      END IF
C
      IF (LDBASN.GT.NDACTU.OR.LDBASP.GT.NDACTU) THEN
          WRITE(6,'(''LDBASN = '',I6,'' LDBASP = '',I6,
     *              '' NDACTU = '',i6)')
     *                LDBASN,LDBASP,NDACTU
          STOP 'NUMBSP VERSUS NDACTU in HFINPU - WSHELL'
      END IF
C
C      Repackaging HF s.p. states into vectors ESINPP and ESINPN
C
      IF (IPAHFB.EQ.1.AND.IMFHFB.EQ.0) THEN
C
          ISTATE=0
          DO IREVER=0,MREVER
             DO I=1,NUMBQP(IREVER,0)
                ISTATE=ISTATE+1
                NSTATE(ISTATE)=ISTATE
                ESINPN(ISTATE)=ESPEQU(I,IREVER,0)
             END DO
          END DO
C
          CALL ORDER1(ESINPN,NSTATE,ISTATE)
C
          DO I=1,ISTATE
             EEQUIV=ESINPN(I)-EFERMN
             OCCU_N(I)=0.5D0*(1.0D0-TANH(EEQUIV/TEMP_T/2.0D0))
          END DO
C
          ISTATE=0
          DO IREVER=0,MREVER
             DO I=1,NUMBQP(IREVER,1)
                ISTATE=ISTATE+1
                NSTATE(ISTATE)=ISTATE
                ESINPP(ISTATE)=ESPEQU(I,IREVER,1)
             END DO
          END DO
C
          CALL ORDER1(ESINPP,NSTATE,ISTATE)
C
          DO I=1,ISTATE
             EEQUIV=ESINPP(I)-EFERMP
             OCCU_P(I)=0.5D0*(1.0D0-TANH(EEQUIV/TEMP_T/2.0D0))
          END DO
C
      END IF
C
      IF (IPAHFB.EQ.0.OR.(IPAHFB.EQ.1.AND.IMFHFB.EQ.1)) THEN
C
          ISTATE=0
          DO IREVER=0,MREVER
             DO I=1,NUMBSP(IREVER,0)
                ISTATE=ISTATE+1
                NSTATE(ISTATE)=ISTATE
                ESINPN(ISTATE)=SPENER(I+IREVER*NUMBSP(0,0),0)
             END DO
          END DO
C
          CALL ORDER1(ESINPN,NSTATE,ISTATE)
C
          DO I=1,ISTATE
             EEQUIV=ESINPN(I)-EFERMN
             OCCU_N(I)=0.5D0*(1.0D0-TANH(EEQUIV/TEMP_T/2.0D0))
          END DO
C
          ISTATE=0
          DO IREVER=0,MREVER
             DO I=1,NUMBSP(IREVER,1)
                ISTATE=ISTATE+1
                NSTATE(ISTATE)=ISTATE
                ESINPP(ISTATE)=SPENER(I+IREVER*NUMBSP(0,1),1)
             END DO
          END DO
C
          CALL ORDER1(ESINPP,NSTATE,ISTATE)
C
          DO I=1,ISTATE
             EEQUIV=ESINPP(I)-EFERMP
             OCCU_P(I)=0.5D0*(1.0D0-TANH(EEQUIV/TEMP_T/2.0D0))
          END DO
C
      END IF
C
C      Repackaging eigenvalues of the free gas into vectors EGAS_P and
C      EGAS_N.
C
      IF (IFSHEL.EQ.2) THEN
C
          ISTATE=0
          DO IREVER=0,MREVER
             DO I=1,NUMBSP(IREVER,0)
                ISTATE=ISTATE+1
                NSTATE(ISTATE)=ISTATE
                EGAS_N(ISTATE)=SPGASN(ISTATE)
             END DO
          END DO
C
          CALL ORDER1(EGAS_N,NSTATE,ISTATE)
C
          ISTATE=0
          DO IREVER=0,MREVER
             DO I=1,NUMBSP(IREVER,1)
                ISTATE=ISTATE+1
                NSTATE(ISTATE)=ISTATE
                EGAS_P(ISTATE)=SPGASP(ISTATE)
             END DO
          END DO
C
          CALL ORDER1(EGAS_P,NSTATE,ISTATE)
C
      END IF
C
C=======================================================================

      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE WSHELL(IN_FIX,IZ_FIX,IFSHEL,MREVER,IPAHFB,IMFHFB,
     *                  IZ_LOW,NOSTPZ,IZSTEP,IN_LOW,NOSTPN,INSTEP,
     *                  LD_BAS,IZNUCL,INNUCL,LDNUCL,DSHELP,DSHELN,
     *                                       EFERMP,EFERMN,NDACTU)
C=======================================================================
      USE hfodd_sizes
C=======================================================================
      DIMENSION
     *          ESINPP(1:2*NDACTU),ESINPN(1:2*NDACTU),
     *          EGAS_P(1:2*NDACTU),EGAS_N(1:2*NDACTU),
     *          OCCU_P(1:2*NDACTU),OCCU_N(1:2*NDACTU)
      DIMENSION
     *          ENCUTP(1:NDNUCL),ENCUTN(1:NDNUCL),
     *          NBCUTP(1:NDNUCL),NBCUTN(1:NDNUCL)
      DIMENSION
     *          IZNUCL(1:NDNUCL),INNUCL(1:NDNUCL)
      DIMENSION
     *          DSHELP(1:NDNUCL),DSHELN(1:NDNUCL),
     *          STRUTP(1:NDNUCL),STRUTN(1:NDNUCL)
C
      COMMON
     *       /STRUTI_PARAMS/ GSTRUN,GSTRUP,HOMFAC,EPSRTN
     *       /STRUTI_POLYNO/ NPOLYN
C
      COMMON
     *       /WSHELL_RES_WS/ ESINGP(1:2*NDBASE),
     *                       EGASSP(1:2*NDBASE),
     *                       OCCUPA(1:2*NDBASE)
      COMMON
     *       /WSHELL_BASSIZ/ LDBASN,LDBASP
      COMMON
     *       /WSHELL_DIMENS/ LDBASE
     *       /WSHELL_KRAMER/ KRAMER
      COMMON
     *       /WSHELL_STRUTP/ MPOLYN
C
C=======================================================================
C
C      IZ_FIX,IN_FIX .......: number of protons and neutrons
C      IFSHEL ..............: traditional (=1) or extended (=2) shell
C                             correction
C      MREVER ..............: contains info on time-reversal symmetry
C                             (1: broken, 0: conserved)
C      ISIMPY ..............: contains info on simplex symmetry
C                             (1: conserved, 0: broken)
C      IMFHFB ..............: activates HFB canonical mean-field
C      IZ_LOW,NOSTPZ,IZSTEP : mesh of the proton numbers for which shell
C                             correction will be computed (lower bound,
C                             number of points, stepsize)
C      IN_LOW,NOSTPN,INSTEP : same as above for neutrons
C      IZNUCL,INNUCL,LDNUCL : arrays storing the proton and neutron numbers
C                             on their respective mesh. LDNUCL is the size
C                             of the mesh.
C      DSHELP,DSHELN  ......: arrays (of size NDNUCL) containing proton
C                             and neutron shell correction
C      NDACTU ..............: alias for HFODD static size NDBASE
C
C=======================================================================
C
      CALL CPUTIM('WSHELL',1)
C
C=======================================================================
C
      KRAMER=1
      MPOLYN=NPOLYN
C
C      Converting HFODD inputs into usable arrays
C
      CALL HFINPU(IFSHEL,MREVER,IPAHFB,IMFHFB,NDACTU,ESINPP,ESINPN,
     *                   EGAS_P,EGAS_N,OCCU_P,OCCU_N,EFERMP,EFERMN)
C
C      Setting up a number of numerical parameters
C
      EPSRTN=0.0001D0
C
      A_MASS=REAL(IN_FIX+IZ_FIX,KIND=KIND(1.0))
C
      XLSTRN=0.0D0
      DNSHEL=0.0D0
      ESHARN=0.0D0
      ESMEAN=0.0D0
C
      XLSTRP=0.0D0
      DNSHEP=0.0D0
      ESHARP=0.0D0
      ESMEAP=0.0D0
C
C=======================================================================
C     Beginning the nucleus if-loop within the deformation loop
C=======================================================================
C
      IZ_UPP=IZ_LOW+IZSTEP*(NOSTPZ-1)
      IN_UPP=IN_LOW+INSTEP*(NOSTPN-1)
C
      AMASS0=REAL(IZ_FIX+IN_FIX,KIND=KIND(1.0))
C
C=======================================================================
C        Strutinsky energy at zero spin and zero frequency - NEUTRONS
C=======================================================================
C
C     I_NUCL refers to the nucleus: calculations are performed for a
C     mesh of nuclei between IZ_LOW and IZ_UPP and IN_LOW and IN_UPP
C     with steps of respectively IZSTEP and INSTEP.
C
      I_NUCL=0
C
      ISOSPI=0
      GSTRUT=GSTRUN
      LDBASE=LDBASN
C
      DO IZ_VAR=IZ_LOW,IZ_UPP,IZSTEP
         DO IN_VAR=IN_LOW,IN_UPP,INSTEP
C
            A_MASS=REAL(IZ_VAR+IN_VAR,KIND=KIND(1.0))
C
            I_NUCL=I_NUCL+1
C
            DO I=1,LDBASE
               ESINGP(I)=ESINPN(I)*(A_MASS/AMASS0)**(1.0D0/3.0D0)
               EGASSP(I)=EGAS_N(I)*(A_MASS/AMASS0)**(1.0D0/3.0D0)
               OCCUPA(I)=OCCU_N(I)
            END DO
C
            CALL SHCORR(IZ_VAR,IN_VAR,ISOSPI,HOMFAC,GSTRUT,NPOLYN,
     *                  EPSRTN,XLSTRN,DNSHEL,ESHARN,ESMEAN,ENCUTS,
     *                                              NBCUTS,IFSHEL)
C
            IZNUCL(I_NUCL)=IZ_VAR
            INNUCL(I_NUCL)=IN_VAR
C
            STRUTN(I_NUCL)=DNSHEL
            ENCUTN(I_NUCL)=ENCUTS
            NBCUTN(I_NUCL)=NBCUTS
            DSHELN(I_NUCL)=DNSHEL
C
         END DO
      END DO
C
C=======================================================================
C        Strutinsky energy at zero spin and zero frequency - PROTONS
C=======================================================================
C
      I_NUCL=I_NUCL-NOSTPZ*NOSTPN
C
      ISOSPI=1
      GSTRUT=GSTRUP
      LDBASE=LDBASP
C
      DO IZ_VAR=IZ_LOW,IZ_UPP,IZSTEP
         DO IN_VAR=IN_LOW,IN_UPP,INSTEP
C
            A_MASS=REAL(IZ_VAR+IN_VAR,KIND=KIND(1.0))
            I_NUCL=I_NUCL+1
C
            DO I=1,LDBASE
               ESINGP(I)=ESINPP(I)*(A_MASS/AMASS0)**(1.0D0/3.0D0)
               EGASSP(I)=EGAS_P(I)*(A_MASS/AMASS0)**(1.0D0/3.0D0)
               OCCUPA(I)=OCCU_P(I)
            END DO
C
            CALL SHCORR(IZ_VAR,IN_VAR,ISOSPI,HOMFAC,GSTRUT,NPOLYN,
     *                  EPSRTN,XLSTRP,DPSHEL,ESHARP,ESMEAP,EPCUTS,
     *                                              NBCUTS,IFSHEL)
C
            STRUTP(I_NUCL)=DPSHEL
            ENCUTP(I_NUCL)=EPCUTS
            NBCUTP(I_NUCL)=NBCUTS
            DSHELP(I_NUCL)=DPSHEL
C
         END DO
      END DO
C
      LDNUCL=I_NUCL
C
C=======================================================================
C
      CALL CPUTIM('WSHELL',0)
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE SHCORR(IZ_FIX,IN_FIX,ISOSPI,HOMFAC,GSTRUT,NPOLYN,
     *                  EPSRTN,XLMBAR,D_SHEL,ESHARP,ESMEAR,E_CUTS,
     *                                              NBCUTS,IFSHEL)
C=======================================================================
      USE hfodd_sizes
C=======================================================================
C
C@@@ TEMP - TEMP - TEMP - TEMP - TEMP - TEMP - TEMP - TEMP - TEMP
      COMMON
     *       /T_FLAG/ IFTEMP
C@@@ TEMP - TEMP - TEMP - TEMP - TEMP - TEMP - TEMP - TEMP - TEMP
      COMMON
     *       /WSHELL_RES_WS/ ESINGP(1:2*NDBASE),
     *                       EGASSP(1:2*NDBASE),
     *                       OCCUPA(1:2*NDBASE)
      COMMON
     *       /WSHELL_STRUTG/ GAMMAS,COEFFM(1:NDSMTH)
      COMMON
     *       /WSHELL_INTSTR/ LIMCUT,N_PART
     *       /WSHELL_DIMENS/ LDBASE
     *       /WSHELL_KRAMER/ KRAMER
C
C=======================================================================
C        Subroutine SHCORR calculates the Strutinsky shell correction
C        to the classical liquid drop model formula.  Comments on the
C        use of the parameters are given in detail in routine  DSHELL
C=======================================================================
C
      A_MASS=REAL(IN_FIX+IZ_FIX,KIND=KIND(1.0))
      HBAR_0=41.0D0/A_MASS**(1.0D0/3.0D0)
C
      IF (ISOSPI.EQ.0) THEN
          N_PART=IN_FIX
      END IF
C
      IF (ISOSPI.EQ.1) THEN
          N_PART=IZ_FIX
      END IF
C
C      Cutting the lowest "LAMBDA+HOMFAC * 41/[A^(1/3)]"  levels;
C      Here the single particle levels are supposed to be doubly
C      degenerate; HOMFAC is of the order of 3.5;
C
                       INUPOL=N_PART/2
      IF (KRAMER.EQ.0) INUPOL=N_PART
C
      E_CUTS=ESINGP(INUPOL)+HOMFAC*HBAR_0   ! HOMFAC is ~ 3.5
C
      LIMCUT=INUPOL-1
      DO WHILE (LIMCUT.LT.LDBASE)
         LIMCUT=LIMCUT+1
         IF (ESINGP(LIMCUT).GT.E_CUTS) EXIT
      END DO
C
      NBCUTS=LIMCUT
C
C      Calculating the Strutinsky shell correction
C
      GAMMAS=GSTRUT*HBAR_0      ! Here GSTRUT is of the order of 1.2
      CALL DSHELL(XLMBAR,EPSRTN,D_SHEL,ESHARP,ESMEAR,NPOLYN,
     *                                               IFSHEL)
C
C      Resetting gammas to its proper value for pairing applications
C
      GAMMAS=GSTRUT*HBAR_0
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE DSHELL(XLMBAR,EPSRTN,D_SHEL,ESHARP,ESMEAR,NPOLYN,
     *                                                     IFSHEL)
C=======================================================================
      USE hfodd_sizes
C=======================================================================
      EXTERNAL
     *          STRUTN
C
      DIMENSION
     *          HERPOL(0:NDSMTH),DHRPOL(0:NDSMTH)
C
C@@@ TEMP - TEMP - TEMP - TEMP - TEMP - TEMP - TEMP - TEMP - TEMP
      COMMON
     *       /T_FLAG/ IFTEMP
      COMMON
     *       /GCTEMP/ TEMP_T
C@@@ TEMP - TEMP - TEMP - TEMP - TEMP - TEMP - TEMP - TEMP - TEMP
      COMMON
     *       /WSHELL_RES_WS/ ESINGP(1:2*NDBASE),
     *                       EGASSP(1:2*NDBASE),
     *                       OCCUPA(1:2*NDBASE)
      COMMON
     *       /WSHELL_STRUTG/ GAMMAS,COEFFM(1:NDSMTH)
      COMMON
     *       /WSHELL_INTSTR/ LIMCUT,N_PART
     *       /WSHELL_DIMENS/ LDBASE
     *       /WSHELL_KRAMER/ KRAMER
C
      DATA
     *        ITERNS /  50  /
      DATA
     *        EXPLIM / 30.0 /
C
C=======================================================================
C
C      Subroutine DSHELL calculates Strutinsky shell correction
C
C      M e a n i n g    o f    s o m e    p a r a m e t e r s :
C
C         GAMMAS - Strutinsky's smearing constant ( in MeV ) [~3.5*hw0]
C         NPOLYN - Maximum order of the smoothing polynomial, usually 6
C         LIMCUT - The number of sing. partcl  states in the summations
C
C         XLMBAR - Resulting Strutinsky Fermi energy  in MeV
C         EPSRTN - Precision in the RTNI type algorithm when searching
C                  for the XLMBAR solution
C
C         D_SHEL - The resulting Strutinsky shell correction. Here and
C                  in the following D_SHEL is the same for an even and
C                  for the neighboring odd particle number
C
C         ESHARP - is the usual sum of all the lowest, occupied single
C                  particle energies
C
C         ESMEAR - is the corresponding Strutinsky-smeared version
C
C         ESINGP - are the single particle  (doubly degenerate) levels
C                  in MeV <--> as all the other energy-like quantities
C
C=======================================================================
C
      PINUMB=4.0D0*ATAN(1.0D0)
      SQRTPI=SQRT(PINUMB)
C
      FACTOR=1.0D0
C
      DO I=1,NDSMTH
         COEFFM(I)=0.0D0
      END DO
C
      DO KFACTO=2,NPOLYN,2
         FACTOR=FACTOR*KFACTO/2
         COEFFM(KFACTO)=(-1)**(KFACTO/2)/(2.0**KFACTO)/FACTOR
      END DO
C
      IF (KRAMER.EQ.0) THEN
          XLAMB0=0.5D0*(ESINGP(N_PART)+ESINGP(N_PART+1))
      ELSE
          XLAMB0=0.5D0*(ESINGP(N_PART/2)+ESINGP(N_PART/2+1))
      END IF
C
      CALL NONLIN(XLMBAR,FVALUE,FDERIV,STRUTN,XLAMB0,EPSRTN,
     *                                 NEEDED,ITERNS,IERROR)
C
      ESHARP=0.0D0
      ESMEAR=0.0D0
C
C=======================================================================
C      Contribution to smooth energy from bound states: e*g_bound(e)
C=======================================================================
C
      DO I_PART=1,LIMCUT
C
         UNSTRT=(XLMBAR-ESINGP(I_PART))/GAMMAS
C
         CALL DF_HER(UNSTRT,NPOLYN,HERPOL,DHRPOL,NDSMTH)
C
         EXPUN2=0.0D0
         UNSTR2=UNSTRT*UNSTRT
C
         SUMPOL=0.0D0
         OCCSTR=0.0D0
C
         IF (UNSTR2.LT.EXPLIM) THEN
C
             EXPUN2=EXP(-UNSTR2)
C
             DO MPOLYN=2,NPOLYN,2
C
                SUMPOL=SUMPOL+COEFFM(MPOLYN)
     *                 *(0.5D0*GAMMAS*HERPOL(MPOLYN)
     *                +ESINGP(I_PART)*HERPOL(MPOLYN-1)
     *                 +MPOLYN*GAMMAS*HERPOL(MPOLYN-2))
C
                OCCSTR=OCCSTR+COEFFM(MPOLYN)*HERPOL(MPOLYN-1)
C
             END DO
C
         END IF
C
         ERROFN=ERRORF(UNSTRT)
C
         ENESUM=+0.5D0*ESINGP(I_PART)*(1.0D0+ERROFN)
     *          -0.5D0*GAMMAS/SQRTPI*EXPUN2
     *          -EXPUN2*SUMPOL/SQRTPI
C
         ESMEAR=ESMEAR+ENESUM
C
      END DO
C
C=======================================================================
C      Contribution to smooth energy from continuum states: e*g_bound(e)
C=======================================================================
C
      IF (IFSHEL.EQ.2) THEN
C
          DO I_PART=1,LIMCUT
C
             UNSTRT=(XLMBAR-EGASSP(I_PART))/GAMMAS
C
             CALL DF_HER(UNSTRT,NPOLYN,HERPOL,DHRPOL,NDSMTH)
C
             EXPUN2=0.0D0
             UNSTR2=UNSTRT*UNSTRT
C
             SUMPOL=0.0D0
             OCCSTR=0.0D0
C
             IF (UNSTR2.LT.EXPLIM) THEN
C
                 EXPUN2=EXP(-UNSTR2)
C
                 DO MPOLYN=2,NPOLYN,2
C
                    SUMPOL=SUMPOL+COEFFM(MPOLYN)
     *                     *(0.5D0*GAMMAS*HERPOL(MPOLYN)
     *                    +EGASSP(I_PART)*HERPOL(MPOLYN-1)
     *                     +MPOLYN*GAMMAS*HERPOL(MPOLYN-2))
C
                    OCCSTR=OCCSTR+COEFFM(MPOLYN)*HERPOL(MPOLYN-1)
C
                 END DO
C
             END IF
C
             ERROFN=ERRORF(UNSTRT)
C
             ENESUM=+0.5D0*EGASSP(I_PART)*(1.+ERROFN)
     *              -0.5D0*GAMMAS/SQRTPI*EXPUN2
     *              -EXPUN2*SUMPOL/SQRTPI
C
             ESMEAR=ESMEAR-ENESUM
C
          END DO
C
      END IF
C
C=======================================================================
C      Computation of the actual shell correction: sum of single-
C      particle nergies minus smooth energy
C=======================================================================
C
      IF (KRAMER.EQ.1) THEN
C
          IMODUL=MOD(N_PART,2) ! For odd nuclei
C
C          Sum of single-particle energies up to the Fermi level
          ESHARP=0.0D0
          DO I_PART=1,N_PART/2
             ESHARP=ESHARP+ESINGP(I_PART)
          END DO
C          To obtain the thermal excitation energy
          ETHERM=0.0D0
          DO I_PART=1,LDBASE
             ETHERM=ETHERM+OCCUPA(I_PART)*ESINGP(I_PART)
          END DO
C
C          Doubling of all quantities to account for Kramers degeneracy
          ESHARP=ESHARP+ESHARP
          ESMEAR=ESMEAR+ESMEAR
          ETHERM=ETHERM+ETHERM
C
C          For odd nuclei, the sum of single-particle energies must also
C          contain the poor lonely particle on level N_PART/2+1
          IMODUL=MOD(N_PART,2)
C
          IF (IMODUL.EQ.1) THEN
              ESHARP=ESHARP+ESINGP(N_PART/2+1)
          END IF
C
      ELSE
C
C          Sum of single-particle energies up to the Fermi level
          ESHARP=0.0D0
          DO I_PART=1,N_PART
             ESHARP=ESHARP+ESINGP(I_PART)
          END DO
C          To obtain the thermal excitation energy
          ETHERM=0.0D0
          DO I_PART=1,LDBASE
             ETHERM=ETHERM+OCCUPA(I_PART)*ESINGP(I_PART)
          END DO
C
      END IF
C
C      Shell correction: sum of s.p. energies minus smooth energy plus
C                        thermal excitation energy
      IF (IFTEMP.EQ.0) THEN
          D_SHEL = ESHARP - ESMEAR
      ELSE
      	  D_SHEL = ETHERM - ESMEAR
      END IF
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE STRUTN(ENEACT,FUNSTR,DFNSTR,IFSHEL)
C=======================================================================
      USE hfodd_sizes
C=======================================================================
      DIMENSION
     *          HERPOL(0:NDSMTH),DHRPOL(0:NDSMTH)
C
      COMMON
     *       /WSHELL_RES_WS/ ESINGP(1:2*NDBASE),
     *                       EGASSP(1:2*NDBASE),
     *                       OCCUPA(1:2*NDBASE)
      COMMON
     *       /WSHELL_INTSTR/ LIMCUT,N_PART
      COMMON
     *       /WSHELL_STRUTG/ GAMMAS,COEFFM(1:NDSMTH)
     *       /WSHELL_STRUTP/ NPOLYN
      COMMON
     *       /WSHELL_KRAMER/ KRAMER
C
      DATA
     *     EXPLIM / 30.0 /
C
C=======================================================================
C
C      This subroutine calculates an auxiliary function  FUNSTR
C      and its derivative DFNSTR,  to find the Strutinsky Fermi
C      level the latter necessery in order to define Strutinsky
C      smearing function g(e).
C
C      Parameters
C
C         ESINGP - Single particle spectrum with single or double
C                  degeneracy
C         EGASSP - Eigenvalues of the free gas
C         HERPOL - Hermite polynomials for the Strutinsky smearing
C         DHRPOL - ... their derivatives used for the same purpose
C         NPOLYN - The order of smoothing in the Strutinski method
C         LIMCUT - The number of proton  or neutron states  summed
C                  in the  'smeared' version of related summations
C         COEFFM - Auxiliary coefficients (calculated outside)
C
C=======================================================================
C
      PINUMB=4.0D0*ATAN(1.0D0)
      SQRTPI=SQRT(PINUMB)
C
      SUMENE=0.0D0
      SUMDER=0.0D0
C
C      Contribution to level density from bound states
C
      DO N=1,LIMCUT
C
         UNSTRT=(ENEACT-ESINGP(N))/GAMMAS
         UNSTR2=UNSTRT*UNSTRT
C
         EXPUN2=0.0D0
C
         SUMDHE=0.0D0
         SUMPOL=0.0D0
C
         IF (UNSTR2.LT.EXPLIM) THEN
C
             EXPUN2=EXP(-UNSTR2)
C
             CALL DF_HER(UNSTRT,NPOLYN,HERPOL,DHRPOL,NDSMTH)
C
             DO MPOLYN=2,NPOLYN,2
C
                SUMPOL=SUMPOL+COEFFM(MPOLYN)*HERPOL(MPOLYN-1)
                SUMDHE=SUMDHE+COEFFM(MPOLYN)*DHRPOL(MPOLYN-1)
C
             END DO
C
         END IF
C
         ENESUM=+0.5D0*(1.+ERRORF(UNSTRT))-EXPUN2*SUMPOL/SQRTPI
         TEXPUN=EXPUN2*(1.+2.00D0*UNSTRT*SUMPOL-SUMDHE)
C
         SUMENE=SUMENE+ENESUM
         SUMDER=SUMDER+TEXPUN
C
      END DO
C
C      Contribution to level density from continuum states
C
      IF (IFSHEL.EQ.2) THEN
C
          DO N=1,LIMCUT
C
             UNSTRT=(ENEACT-EGASSP(N))/GAMMAS
             UNSTR2=UNSTRT*UNSTRT
C
             EXPUN2=0.0D0
C
             SUMDHE=0.0D0
             SUMPOL=0.0D0
C
             IF (UNSTR2.LT.EXPLIM) THEN
C
                 EXPUN2=EXP(-UNSTR2)
C
                 CALL DF_HER(UNSTRT,NPOLYN,HERPOL,DHRPOL,NDSMTH)
C
                 DO MPOLYN=2,NPOLYN,2
C
                    SUMPOL=SUMPOL+COEFFM(MPOLYN)*HERPOL(MPOLYN-1)
                    SUMDHE=SUMDHE+COEFFM(MPOLYN)*DHRPOL(MPOLYN-1)
C
                 END DO
C
             END IF
C
             ENESUM=+0.5D0*(1.+ERRORF(UNSTRT))-EXPUN2*SUMPOL/SQRTPI
             TEXPUN=EXPUN2*(1.+2.00D0*UNSTRT*SUMPOL-SUMDHE)
C
             SUMENE=SUMENE-ENESUM
             SUMDER=SUMDER-TEXPUN
C
          END DO
C
      END IF
C
      IF (KRAMER.EQ.1) THEN
          SUMENE=SUMENE+SUMENE
          SUMDER=SUMDER+SUMDER
      END IF
C
      FUNSTR=+N_PART-SUMENE
      DFNSTR=-SUMDER/GAMMAS/SQRTPI
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      FUNCTION ERRORF(ARGMNT)
      DATA
     *      A1COEF,A2COEF,A3COEF,A4COEF,A5COEF,A6COEF
     *
     *    /.0705230784D0, .0422820123D0, .0092705272D0,
     *     .0001520143D0, .0002765672D0, .0000430638D0 /
C
C=======================================================================
C        This program calculates the standard error function
C=======================================================================
C
      XAUXIL=ABS(ARGMNT)
      ZAUXIL=(XAUXIL*(A1COEF
     *       +XAUXIL*(A2COEF
     *       +XAUXIL*(A3COEF
     *       +XAUXIL*(A4COEF
     *       +XAUXIL*(A5COEF+XAUXIL*A6COEF)))))+1.0D0)
C
      IF (ABS(ZAUXIL).LT.100.0D0) GO TO 1
C
      ERRABS=1.0D0
C
      GO TO 2
C
   1  ERRABS=1.0D0-1.0D0/ZAUXIL**16
C
   2  CONTINUE
C
      ERRORF=SIGN(REAL(ERRABS,KIND=KIND(1.0)),
     *            REAL(ARGMNT,KIND=KIND(1.0)))
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE NONLIN(X_ROOT,FVALUE,FDERIV,FUNCTN,XFIRST,EPSRTN,
     *                                       NEEDED,ITERNS,IERROR)
      EXTERNAL
     *         FUNCTN
C
C=======================================================================
C
C        This subroutine is known in circles of connessairs under
C        the fantasy name RTNI. It solves a nonlinear equation of
C        the form
C                               FUNCTN(XARGUM) = 0
C        by Newton's method.
C
C=======================================================================
C
C        X_ROOT - On exit: root found of the equation in question
C
C        FVALUE - On exit: Function value at the argument x= root
C        FDERIV - On exit: The coreresponding derivative of FUNCT
C
C        FUNCTN - External subroutine
C
C                          FUNCTN(XARGUM,FUNVAL,FDERIV)
C
C                 It computes the function value FUNVAL, and its
C                 derivative  FDERIV  at  XARGUM
C
C        XFIRST - On entry: an initial guess of  XARGUM
C
C        EPSRTN - If |FUNVAL| < 100 * EPSRTN the iterations stop
C
C        ITERNS - Maximum allowed number  of Newton's iterations
C        NEEDED - The number  of iterations  actually  performed
C
C        IERROR - Resultant error parameter. The possible values
C                 are:
C                      IERROR=0 - solution fouind within
C                                 prescribed accuracy
C                      IERROR=1 - no convergence after  ITERNS
C                                 iterations
C                      IERROR=2 - at some iteration derivative
C                                 FDERIV was equal to 0,  thus
C                                 the next Newton iteration is
C                                 impossible
C
C=======================================================================
C
      IERROR=0
C
C=======================================================================
C        Prepare the first iteration
C=======================================================================
C
      X_ROOT=XFIRST
      SEARCH=X_ROOT
C
      CALL FUNCTN(SEARCH,FVALUE,FDERIV,IFSHEL)
C
      EPS100=100.0D0*EPSRTN
C
C=======================================================================
C        Start iteration loop
C=======================================================================
C
      DO I=1,ITERNS
C
         NEEDED=I
C
         IF (ABS(FVALUE).EQ.0.0) THEN
C
C            Equation was already satisfied by the argument
C
             RETURN
C
         ELSE
C
C            ... equation is not yet satisfied by X
C
             IF (ABS(FDERIV).EQ.0.0) THEN
C
C                Here the next iteration is NOT possible
C                =>>ERROR RETURN in case of zero divisor
C
                 IERROR=2
                 RETURN
C
             ELSE
C
C                ... and here the next iteration is possible
C
                 DXVALU=FVALUE/FDERIV
                 X_ROOT=X_ROOT-DXVALU
                 SEARCH=X_ROOT
C
                 CALL FUNCTN(SEARCH,FVALUE,FDERIV,IFSHEL)
C
C                Test on possibly satisfactory accuracy
C
                 SEARCH=EPSRTN
                 ABSARG=ABS(X_ROOT)
C
                 IF ((ABSARG-1.).GT.0.0) THEN
                     SEARCH=SEARCH*ABSARG
                 END IF
C
                 IF ((ABS(DXVALU)-SEARCH).GT.0.0) THEN
                     GO TO 1
                 ELSE
C
                     IF ((ABS(FVALUE)-EPS100).GT.0.0) THEN
                         GO TO 1
                     ELSE
                         RETURN
                     END IF
                 END IF
             END IF
         END IF
C
    1    CONTINUE
C
      END DO !  End of iteration loop
C
C     No convergence after iterns iteration steps, error return.
C
      IERROR=1
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE DF_HER(ARGHER,NORHER,HERPOL,DHERPO,NDHERM)
      DIMENSION
     *          HERPOL(0:NDHERM),DHERPO(0:NDHERM)
C
C=======================================================================
C
C         DFHERM calculates Hermite polynomials at point X=ARGHER up
C         to N=NORHER-th degree.  It also calculates the derivatives
C         DHERPO of these polynomials
C
C                            DHERPO(I)=2*I*HERPOL(I-1)
C
C=======================================================================
C
      HERPOL(0)=1.0D0
      DHERPO(0)=0.0D0
C
      IF (NORHER.EQ.0) RETURN
C
      HERPOL(1)=ARGHER+ARGHER
      DHERPO(1)=2.0000D0
C
      IF (NORHER.EQ.1) RETURN
C
      DO I=1,NORHER-1
C
         POLYNM=ARGHER*HERPOL(I)-REAL(I,KIND=KIND(1.0))*HERPOL(I-1)
         HERPOL(I+1)=POLYNM+POLYNM
C
         DPOLYN=REAL(I+1,KIND=KIND(1.0))*HERPOL(I)
         DHERPO(I+1)=DPOLYN+DPOLYN
C
      END DO
C
C=======================================================================
C
      RETURN
      END
C
