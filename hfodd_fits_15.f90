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

module hfodd_fits
      contains

      SUBROUTINE FITJAD(JADRO_NAME, DAT_SUFFIX)
      USE HFODD_SIZES
      CHARACTER JADRO_NAME*5
      CHARACTER DAT_SUFFIX*7
      CHARACTER DAT_FILE_NAME*64

      COMMON /CFIREA/ NFIREA
      COMMON /CFIPRI/ NFIPRI

      NFIREA=87
      WRITE (DAT_FILE_NAME, '(''data/'',A5,''-'',A7,''.dat'')') JADRO_NAME, DAT_SUFFIX
      OPEN(NFIREA,FILE=DAT_FILE_NAME,STATUS='UNKNOWN',FORM='FORMATTED')

      NFIPRI=88
      OPEN(NFIPRI,FILE='/dev/null',STATUS='UNKNOWN',FORM='FORMATTED')
      END SUBROUTINE FITJAD

  !=====================================================================
  !
  !     This subroutine calculates
  !=====================================================================

      SUBROUTINE FITINP(PARINP,NUMINP,I_TYPE,N3LORD)
      USE HFODD_SIZES
      INTEGER I_TYPE
      DIMENSION PARINP(NUMINP)
      COMMON /SPMSKY/ SRHO_P,SRHO_M,SRHODP,SRHODM, \
                      SLPR_P,SLPR_M, \
                      STAU_P,STAU_M, \
                      SSCU_P,SSCU_M, \
                      SDIV_P,SDIV_M, \
                      SSPI_P,SSPI_M,SSPIDP,SSPIDM, \
                      SLPS_P,SLPS_M, \
                      SCUR_P,SCUR_M, \
                      SKIS_P,SKIS_M, \
                      SROT_P,SROT_M
      COMMON /REGBIN/ REGBCC(NDREGA,NDFORC),IREREG
      COMMON /REGVIN/ REGVCC(NDREGA,NDFORC)
      COMMON /GOGFOR/ GOGSTR(NDGOGA,NDFORC),GOGWID(NDGOGA),NUGOGA

      ISTYPE=0

      IF (I_TYPE.EQ.1) THEN

          IF (NUMINP.NE.12) STOP ' INCREASE NUMINP 1'

          SRHO_P=SRHO_P*PARINP(1)
          SRHODP=SRHODP*PARINP(2)
          SLPR_P=SLPR_P*PARINP(3)
          STAU_P=STAU_P*PARINP(4)
          SSCU_P=SSCU_P*PARINP(5)
          SDIV_P=SDIV_P*PARINP(6)

          SRHO_M=SRHO_M*PARINP(7)
          SRHODM=SRHODM*PARINP(8)
          SLPR_M=SLPR_M*PARINP(9)
          STAU_M=STAU_M*PARINP(10)
          SSCU_M=SSCU_M*PARINP(11)
          SDIV_M=SDIV_M*PARINP(12)

          ISTYPE=1

      END IF

      IF (I_TYPE.EQ.2) THEN

          IF (NUMINP.NE.IREREG*NDFORC+2) STOP ' INCREASE NUMINP 2'

          IND=0
          DO INREGA=1,IREREG
             DO INFORC=1,NDFORC
                IND=IND+1
                IF (N3LORD.GE.0) THEN
                REGBCC(INREGA,INFORC)=REGBCC(INREGA,INFORC)*PARINP(IND)
                ELSE
                REGVCC(INREGA,INFORC)=REGVCC(INREGA,INFORC)*PARINP(IND)
                END IF
             END DO
          END DO

          SDIV_P=SDIV_P*PARINP(IND+1)
          SDIV_M=SDIV_M*PARINP(IND+2)

          ISTYPE=1

      END IF

      IF (I_TYPE.EQ.3) THEN

          IF (NUMINP.NE.IREREG*NDFORC+1) STOP ' INCREASE NUMINP 3'

          IND=0
          DO INREGA=1,IREREG
             DO INFORC=1,NDFORC
                IND=IND+1
                IF (N3LORD.GE.0) THEN
                REGBCC(INREGA,INFORC)=REGBCC(INREGA,INFORC)*PARINP(IND)
                ELSE
                REGVCC(INREGA,INFORC)=REGVCC(INREGA,INFORC)*PARINP(IND)
                END IF
             END DO
          END DO

          SDIV_P=SDIV_P*PARINP(IND+1)
          SDIV_M=SDIV_M*PARINP(IND+1)

          ISTYPE=1

      END IF

     IF (I_TYPE.EQ.4) THEN

          IF (NUMINP.NE.IREREG*NDFORC+2) STOP ' INCREASE NUMINP 4'

          IND=0
          DO INREGA=1,IREREG
             DO INFORC=1,NDFORC
                IND=IND+1
                IF (N3LORD.GE.0) THEN
                REGBCC(INREGA,INFORC)=REGBCC(INREGA,INFORC)*PARINP(IND)
                ELSE
                REGVCC(INREGA,INFORC)=REGVCC(INREGA,INFORC)*PARINP(IND)
                END IF
             END DO
          END DO

          SDIV_P=SDIV_P*PARINP(IND+1)
          SDIV_M=SDIV_M*PARINP(IND+1)
          SRHODP=SRHODP*PARINP(IND+2)
          SRHODM=SRHODM*PARINP(IND+2)

          ISTYPE=1

      END IF

      IF (I_TYPE.EQ.5) THEN

          IF (NUMINP.NE.8) STOP ' INCREASE NUMINP 5'

          SRHO_P=SRHO_P*PARINP(1)
          SLPR_P=SLPR_P*PARINP(2)
          STAU_P=STAU_P*PARINP(3)

          SRHO_M=SRHO_M*PARINP(4)
          SLPR_M=SLPR_M*PARINP(5)
          STAU_M=STAU_M*PARINP(6)

          SDIV_P=SDIV_P*PARINP(7)
          SDIV_M=SDIV_M*PARINP(7)
          SRHODP=SRHODP*PARINP(8)
          SRHODM=SRHODM*PARINP(8)

          ISTYPE=1

      END IF

      IF (I_TYPE.EQ.6) THEN

          IF (NUMINP.NE.10) STOP ' INCREASE NUMINP 6'

          GOGSTR(1,1)=GOGSTR(1,1)*PARINP(1)
          GOGSTR(1,2)=GOGSTR(1,2)*PARINP(2)
          GOGSTR(1,3)=GOGSTR(1,3)*PARINP(3)
          GOGSTR(1,4)=GOGSTR(1,4)*PARINP(4)

          GOGSTR(2,1)=GOGSTR(2,1)*PARINP(5)
          GOGSTR(2,2)=GOGSTR(2,2)*PARINP(6)
          GOGSTR(2,3)=GOGSTR(2,3)*PARINP(7)
          GOGSTR(2,4)=GOGSTR(2,4)*PARINP(8)

          SDIV_P=SDIV_P*PARINP(9)
          SDIV_M=SDIV_M*PARINP(9)
          SRHODP=SRHODP*PARINP(10)
          SRHODM=SRHODM*PARINP(10)

          ISTYPE=1

      END IF

      IF (I_TYPE.EQ.7) THEN

          IF (NUMINP.NE.9) STOP ' INCREASE NUMINP 1'

          SRHO_P=SRHO_P*PARINP(1)
          SLPR_P=SLPR_P*PARINP(2)
          STAU_P=STAU_P*PARINP(3)
          SSCU_P=SSCU_P*PARINP(4)

          SRHO_M=SRHO_M*PARINP(5)
          SLPR_M=SLPR_M*PARINP(6)
          STAU_M=STAU_M*PARINP(7)
          SSCU_M=SSCU_M*PARINP(8)

          SDIV_P=SDIV_P*PARINP(9)
          SDIV_M=SDIV_M*PARINP(9)

          ISTYPE=1

      END IF

      IF (ISTYPE.EQ.0) STOP ' INCORRECT I_TYPE IN FITINP'

      END SUBROUTINE FITINP


  !=====================================================================
  !
  !     This subroutine calculates
  !=====================================================================

      SUBROUTINE FITOUT(PAROUT,NUMOUT)
      USE HFODD_SIZES
      DIMENSION PAROUT(NUMOUT)
      COMMON \
            /ALLENE/ EKIN_N,EKIN_P,EKIN_T, \
                     EPOT_N,EPOT_P,EPOT_T, \
                     ESUM_N,ESUM_P,ESUM_T, \
                     EPAI_N,EPAI_P,EPAI_T, \
                     EREA_N,EREA_P,EREA_T, \
                     ELIP_N,ELIP_P,ELIP_T, \
     \
                     ECOULD,ECOULE,ECOULT, \
                            ECOULS,ECOULV, \
     \
                     EMULCO,EMUSLO,EMUREA, \
                     ESIFCO,ESISLO,ESIREA, \
                     ESPICO,ESPSLO,ESPREA, \
     \
                     ENREAR,ECORCM,ECOR_R, \
     \
                     EEVEW0,EODDW0,ENE_W0, \
                     ENEVEN,ENEODD,ENESKY, \
                     ESTABN,ETOTSP,ETOTFU
      COMMON  \
             /RMSRAD/ RADI_N(0:NDKART), \
                      RADI_P(0:NDKART), \
                      RADI_T(0:NDKART)
      COMMON  \
             /QMULTI/ QMUL_N(0:NDMULT,-NDMULT:NDMULT),  \
                      QMUL_P(0:NDMULT,-NDMULT:NDMULT),  \
                      QMUL_T(0:NDMULT,-NDMULT:NDMULT)
      COMMON /CFIREA/ NFIREA

      IF (NUMOUT.NE.4) STOP ' INCREASE NUMOUT'

      CLOSE(NFIREA)
      PAROUT(1)=ETOTFU
      IF (ABS(ESTABN).LT.0.001) THEN
          PAROUT(2)= 1.
      ELSE
          PAROUT(2)=-1.
      ENDIF
      PAROUT(3)=RADI_P(0)
      PAROUT(4)=QMUL_P(2,0)
      END SUBROUTINE FITOUT

    !====================================================================
end module hfodd_fits
