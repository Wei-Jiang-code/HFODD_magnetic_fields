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
C
C=======================================================================
      MODULE      HE_DEN
      SAVE        DE_RHO,DE_TAU,DE_LPR,DE_DIV,
     *            DE_SPI,DE_KIS,DE_GRR,DE_LPS,DE_ROS,DE_ROC,DE_CUR,
     *            DE_SCU,DE_DES
      COMPLEX     DE_RHO,DE_TAU,DE_LPR,DE_DIV,
     *            DE_SPI,DE_KIS,DE_GRR,DE_LPS,DE_ROS,DE_ROC,DE_CUR,
     *            DE_SCU,DE_DES
      ALLOCATABLE DE_RHO(:,:,:,:),DE_TAU(:,:,:,:),
     *            DE_LPR(:,:,:,:),DE_DIV(:,:,:,:)
      ALLOCATABLE DE_SPI(:,:,:,:,:),DE_KIS(:,:,:,:,:),
     *            DE_GRR(:,:,:,:,:),DE_LPS(:,:,:,:,:),
     *            DE_ROS(:,:,:,:,:),DE_ROC(:,:,:,:,:),
     *            DE_CUR(:,:,:,:,:)
      ALLOCATABLE DE_SCU(:,:,:,:,:,:),DE_DES(:,:,:,:,:,:)
      END MODULE  HE_DEN
C=======================================================================
      MODULE      PD_DEN
      SAVE        PD_RHO,PP_RHO,PD_TAU,PP_TAU,PD_LPR,PP_LPR,
     *            PD_GRR,PP_GRR,PD_SCU,PP_SCU
      COMPLEX     PD_RHO,PP_RHO,PD_TAU,PP_TAU,PD_LPR,PP_LPR,
     *            PD_GRR,PP_GRR,PD_SCU,PP_SCU
      ALLOCATABLE PD_RHO(:,:,:,:),PP_RHO(:,:,:,:),
     *            PD_TAU(:,:,:,:),PP_TAU(:,:,:,:),
     *            PD_LPR(:,:,:,:),PP_LPR(:,:,:,:)
      ALLOCATABLE PD_GRR(:,:,:,:,:),PP_GRR(:,:,:,:,:)
      ALLOCATABLE PD_SCU(:,:,:,:,:,:),PP_SCU(:,:,:,:,:,:)
      END MODULE  PD_DEN
C=======================================================================
      MODULE      VE_FLD
      SAVE        VE_MAS,VE_CEN,VE_KIS,VE_SPI,VE_CUR,VE_SOR
      COMPLEX     VE_MAS,VE_CEN,VE_KIS,VE_SPI,VE_CUR,VE_SOR
      ALLOCATABLE VE_MAS(:,:,:,:),VE_CEN(:,:,:,:),
     *            VE_KIS(:,:,:,:,:),VE_SPI(:,:,:,:,:),
     *            VE_CUR(:,:,:,:,:),
     *            VE_SOR(:,:,:,:,:,:)
      END MODULE  VE_FLD
C=======================================================================
      MODULE      WD_FLD
      SAVE        WD_CEN,WD_TAU,WD_SOR
      COMPLEX     WD_CEN,WD_TAU,WD_SOR
      ALLOCATABLE WD_CEN(:,:,:,:),WD_TAU(:,:,:,:),WD_SOR(:,:,:,:,:,:)
      END MODULE  WD_FLD
C=======================================================================
      MODULE      BROYIN
      SAVE        BRO_IN
      ALLOCATABLE BRO_IN(:)
      END MODULE  BROYIN
C=======================================================================
      MODULE      BROYOU
      SAVE        BRO_OU
      ALLOCATABLE BRO_OU(:)
      END MODULE  BROYOU
C=======================================================================
      MODULE      BROYIN_NEUTRON
      SAVE        BRONIN
      ALLOCATABLE BRONIN(:)
      END MODULE  BROYIN_NEUTRON
      MODULE      BROYOU_NEUTRON
      SAVE        BRONOU
      ALLOCATABLE BRONOU(:)
      END MODULE  BROYOU_NEUTRON
C=======================================================================
      MODULE      BROYIN_PROTON
      SAVE        BROPIN
      ALLOCATABLE BROPIN(:)
      END MODULE  BROYIN_PROTON
      MODULE      BROYOU_PROTON
      SAVE        BROPOU
      ALLOCATABLE BROPOU(:)
      END MODULE  BROYOU_PROTON
C=======================================================================
      MODULE      MAT_PP
      SAVE        BIG_PP
      COMPLEX     BIG_PP
      ALLOCATABLE BIG_PP(:,:,:)
      END MODULE  MAT_PP
C=======================================================================
      MODULE      MAT_PM
      SAVE        BIG_PM
      COMPLEX     BIG_PM
      ALLOCATABLE BIG_PM(:,:,:)
      END MODULE  MAT_PM
C=======================================================================
      MODULE      MAR_PP
      SAVE        RHO_PP
      COMPLEX     RHO_PP
      ALLOCATABLE RHO_PP(:,:,:)
      END MODULE  MAR_PP
C=======================================================================
      MODULE      MAR_PM
      SAVE        RHO_PM
      COMPLEX     RHO_PM
      ALLOCATABLE RHO_PM(:,:,:)
      END MODULE  MAR_PM
C=======================================================================
      MODULE      KAP_PP
      SAVE        PAI_PP
      COMPLEX     PAI_PP
      ALLOCATABLE PAI_PP(:,:,:)
      END MODULE  KAP_PP
C=======================================================================
      MODULE      KAP_PM
      SAVE        PAI_PM
      COMPLEX     PAI_PM
      ALLOCATABLE PAI_PM(:,:,:)
      END MODULE  KAP_PM
C=======================================================================
      MODULE      KAP2PP
      SAVE        PAI2PP
      COMPLEX     PAI2PP
      ALLOCATABLE PAI2PP(:,:,:)
      END MODULE  KAP2PP
C=======================================================================
      MODULE      MAD_PP
      SAVE        DEN_PP
      COMPLEX     DEN_PP
      ALLOCATABLE DEN_PP(:,:,:,:)
      END MODULE  MAD_PP
C=======================================================================
      MODULE      MAD_PM
      SAVE        DEN_PM
      COMPLEX     DEN_PM
      ALLOCATABLE DEN_PM(:,:,:,:)
      END MODULE  MAD_PM
C=======================================================================
      MODULE      MLI_PP
      SAVE        DLI_PP
      COMPLEX     DLI_PP
      ALLOCATABLE DLI_PP(:,:,:,:)
      END MODULE  MLI_PP
C=======================================================================
      MODULE      MLI_PM
      SAVE        DLI_PM
      COMPLEX     DLI_PM
      ALLOCATABLE DLI_PM(:,:,:,:)
      END MODULE  MLI_PM
C=======================================================================
      MODULE      MAP_PP
      SAVE        CHI_PP
      COMPLEX     CHI_PP
      ALLOCATABLE CHI_PP(:,:,:,:)
      END MODULE  MAP_PP
C=======================================================================
      MODULE      MAP_PM
      SAVE        CHI_PM
      COMPLEX     CHI_PM
      ALLOCATABLE CHI_PM(:,:,:,:)
      END MODULE  MAP_PM
C=======================================================================
      MODULE      MAF_PP
      SAVE        FIL_PP
      COMPLEX     FIL_PP
      ALLOCATABLE FIL_PP(:,:,:,:)
      END MODULE  MAF_PP
C=======================================================================
      MODULE      MAF_PM
      SAVE        FIL_PM
      COMPLEX     FIL_PM
      ALLOCATABLE FIL_PM(:,:,:,:)
      END MODULE  MAF_PM
C=======================================================================
      MODULE      MAPIDE
      SAVE        FIPIDE
      SAVE        FIP2DE
      COMPLEX     FIPIDE
      COMPLEX     FIP2DE
      ALLOCATABLE FIPIDE(:,:,:,:)
      ALLOCATABLE FIP2DE(:,:,:,:)
      END MODULE  MAPIDE
C=======================================================================
      MODULE      HCOULO
      SAVE        HPPCOU
      SAVE        HPMCOU
      COMPLEX     HPPCOU
      COMPLEX     HPMCOU
      ALLOCATABLE HPPCOU(:,:)
      ALLOCATABLE HPMCOU(:,:)
      END MODULE  HCOULO
C=======================================================================
      MODULE      HYUKAW
      SAVE        HPPYUK
      SAVE        HPMYUK
      COMPLEX     HPPYUK
      COMPLEX     HPMYUK
      ALLOCATABLE HPPYUK(:,:,:)
      ALLOCATABLE HPMYUK(:,:,:)
      END MODULE  HYUKAW
C=======================================================================
      MODULE      PAIDEL
      SAVE        HAMIDE
      SAVE        HAM2DE
      COMPLEX     HAMIDE
      COMPLEX     HAM2DE
      ALLOCATABLE HAMIDE(:,:,:)
      ALLOCATABLE HAM2DE(:,:,:)
      END MODULE  PAIDEL
C=======================================================================
      MODULE      SPINOR
      SAVE        WSWAVE
      COMPLEX     WSWAVE
      ALLOCATABLE WSWAVE(:,:,:)
      END MODULE  SPINOR
C=======================================================================
      MODULE      WAVSQU
      SAVE        SQUWAV
      ALLOCATABLE SQUWAV(:,:,:,:)
      END MODULE  WAVSQU
C=======================================================================
      MODULE      WAVR_L
      SAVE        WARIGH
      SAVE        WALEFT
      COMPLEX     WARIGH
      COMPLEX     WALEFT
      ALLOCATABLE WARIGH(:,:,:)
      ALLOCATABLE WALEFT(:,:,:)
      END MODULE  WAVR_L
C=======================================================================
      MODULE      SAVRIG
      SAVE        SARIGH
      COMPLEX     SARIGH
      ALLOCATABLE SARIGH(:,:,:,:)
      END MODULE  SAVRIG
C=======================================================================
      MODULE      SAVLEF
      SAVE        SALEFT
      COMPLEX     SALEFT
      ALLOCATABLE SALEFT(:,:,:,:)
      END MODULE  SAVLEF
C=======================================================================
      MODULE      ALLWAV
      SAVE        WAVOCC
      COMPLEX     WAVOCC
      ALLOCATABLE WAVOCC(:,:,:)
      END MODULE  ALLWAV
C=======================================================================
      MODULE      CANBAS
c     Canonical wave functions and u- and v- amplitudes
c     WAVCAN(IBASE, ISTATE, IREVER, ICHARG)
c     V_CAN (       ISTATE, IREVER, ICHARG)
c     U_CAN (       ISTATE, IREVER, ICHARG)
c     Where IBASE  = index of basis state
c     ISTATE = index of can.  state
c           IREVER = index of time-reversal block
c           ICHARG = index of charge (n/p) block
      COMPLEX, SAVE, ALLOCATABLE :: WAVCAN(:,:,:,:)
      REAL,    SAVE, ALLOCATABLE :: V_CAN(:,:,:), U_CAN(:,:,:)
      END MODULE  CANBAS
C=======================================================================
      MODULE      PNMWAV
      SAVE        PNMOCC
      COMPLEX     PNMOCC
      ALLOCATABLE PNMOCC(:,:,:)
      END MODULE  PNMWAV
C=======================================================================
      MODULE      ALLSIG
      SAVE        EWASIG,SPESVG,LDMEFG
      INTEGER     LDMEFG
      COMPLEX     EWASIG,SPESVG
      ALLOCATABLE EWASIG(:,:,:,:),SPESVG(:,:,:),LDMEFG(:,:)
      END MODULE  ALLSIG
C=======================================================================
      MODULE      ALLSIM
      SAVE        EWASIM,SPESVM,LDMEFM
      INTEGER     LDMEFM
      COMPLEX     EWASIM,SPESVM
      ALLOCATABLE EWASIM(:,:,:),SPESVM(:,:),LDMEFM(:)
      END MODULE  ALLSIM
C=======================================================================
      MODULE      ALLSIQ
      SAVE        EWASIQ,SPESVQ,LDMEFQ
      INTEGER     LDMEFQ
      COMPLEX     EWASIQ,SPESVQ
      ALLOCATABLE EWASIQ(:,:,:),SPESVQ(:,:),LDMEFQ(:)
      END MODULE  ALLSIQ
C=======================================================================
      MODULE      ALLSIZ
      SAVE        EWASIZ
      COMPLEX     EWASIZ
      ALLOCATABLE EWASIZ(:,:)
      END MODULE  ALLSIZ
C=======================================================================
      MODULE      ALLMIZ
      SAVE        EWAMIZ
      COMPLEX     EWAMIZ
      ALLOCATABLE EWAMIZ(:,:)
      END MODULE  ALLMIZ
C=======================================================================
      MODULE      SAVQUA
      SAVE        ASVQUA         ,BSVQUA
      COMPLEX     ASVQUA         ,BSVQUA
      ALLOCATABLE ASVQUA(:,:,:,:),BSVQUA(:,:,:,:)
      END MODULE  SAVQUA
C
      MODULE      ALLQUA
      SAVE        AWAQUA       ,BWAQUA
      COMPLEX     AWAQUA       ,BWAQUA
      ALLOCATABLE AWAQUA(:,:,:),BWAQUA(:,:,:)
      END MODULE  ALLQUA
C=======================================================================
      MODULE      ALLQUZ
      SAVE        AWAQUZ       ,BWAQUZ
      COMPLEX     AWAQUZ       ,BWAQUZ
      ALLOCATABLE AWAQUZ(:,:,:),BWAQUZ(:,:,:)
      END MODULE  ALLQUZ
C=======================================================================
      MODULE      BLOSAV
      TYPE BLAB
      INTEGER::N,IR
      END TYPE
      SAVE::BLOLAB,AZERO,BZERO,WAZERO,WBZERO
      TYPE(BLAB),ALLOCATABLE::BLOLAB(:)
      INTEGER,ALLOCATABLE::WLAB_B(:),WLAB_A(:)
      REAL,ALLOCATABLE::V2BLO(:)
      COMPLEX,ALLOCATABLE::A_ZERO(:,:,:),B_ZERO(:,:,:)
      END MODULE  BLOSAV
C=======================================================================
      MODULE      REDMOM
      SAVE        REDQ_P   ,REDM_P   ,REDS_P   ,INDRED
      COMPLEX     REDQ_P   ,REDM_P   ,REDS_P
      ALLOCATABLE REDQ_P(:),REDM_P(:),REDS_P(:),INDRED(:,:,:)
      END MODULE  REDMOM
C=======================================================================
      MODULE      REEMOM
      SAVE        REEQ_P   ,REEM_P   ,REES_P   ,INDREE
      COMPLEX     REEQ_P   ,REEM_P   ,REES_P
      ALLOCATABLE REEQ_P(:),REEM_P(:),REES_P(:),INDREE(:,:,:)
      END MODULE  REEMOM
C=======================================================================
      MODULE      GAUSTO
      SAVE        GAUFIL
      COMPLEX     GAUFIL
      ALLOCATABLE GAUFIL(:,:,:,:,:)
      SAVE        GAUAVR
      COMPLEX     GAUAVR
      ALLOCATABLE GAUAVR(:,:,:,:,:)
      END MODULE  GAUSTO
C=======================================================================
      MODULE      YUKSTO
      SAVE        YUKASV
      COMPLEX     YUKASV
      ALLOCATABLE YUKASV(:,:,:,:)
      END MODULE  YUKSTO
C=======================================================================
      MODULE      GOGSTO
      SAVE        GOGASV
      COMPLEX     GOGASV
      ALLOCATABLE GOGASV(:,:,:,:)
      END MODULE  GOGSTO
C=======================================================================
      MODULE      GOGKAP
      SAVE        GOGCHI
      COMPLEX     GOGCHI
      ALLOCATABLE GOGCHI(:,:,:,:)
      END MODULE  GOGKAP
C=======================================================================
      MODULE      REGSTO
      SAVE        REGASV
      COMPLEX     REGASV
      ALLOCATABLE REGASV(:,:,:,:)
      END MODULE  REGSTO
C=======================================================================
      MODULE      PAISTO
      SAVE        DEL_PP
      SAVE        DEL_PM
      COMPLEX     DEL_PP
      COMPLEX     DEL_PM
      ALLOCATABLE DEL_PP(:,:,:,:)
      ALLOCATABLE DEL_PM(:,:,:,:)
      END MODULE  PAISTO
C=======================================================================
      MODULE      ALISTO
      SAVE        ALINPP
      SAVE        ALINPM
      COMPLEX     ALINPP
      COMPLEX     ALINPM
      ALLOCATABLE ALINPP(:,:,:,:,:)
      ALLOCATABLE ALINPM(:,:,:,:,:)
      END MODULE  ALISTO
C=======================================================================
      MODULE      PLISTO
      SAVE        PLINPP
      SAVE        PLINPM
      COMPLEX     PLINPP
      COMPLEX     PLINPM
      ALLOCATABLE PLINPP(:,:,:,:,:)
      ALLOCATABLE PLINPM(:,:,:,:,:)
      END MODULE  PLISTO
C=======================================================================
      MODULE      AROSTO
      SAVE        AROTPP
      SAVE        AROTPM
      COMPLEX     AROTPP
      COMPLEX     AROTPM
      ALLOCATABLE AROTPP(:,:,:,:,:)
      ALLOCATABLE AROTPM(:,:,:,:,:)
      END MODULE  AROSTO
C=======================================================================
      MODULE      PROSTO
      SAVE        PROTPP
      SAVE        PROTPM
      COMPLEX     PROTPP
      COMPLEX     PROTPM
      ALLOCATABLE PROTPP(:,:,:,:,:)
      ALLOCATABLE PROTPM(:,:,:,:,:)
      END MODULE  PROSTO
C=======================================================================
      MODULE      COUSTO
      SAVE        COULSV
      COMPLEX     COULSV
      ALLOCATABLE COULSV(:,:,:)
      END MODULE  COUSTO
C=======================================================================
      MODULE      CYLTRA
      SAVE        AUXCYL
      COMPLEX     AUXCYL
      ALLOCATABLE AUXCYL(:,:)
      END MODULE  CYLTRA
C=======================================================================
      MODULE      PNMXPP
      SAVE        PNM_PP
      COMPLEX     PNM_PP
      ALLOCATABLE PNM_PP(:,:,:)
      END MODULE  PNMXPP
C=======================================================================
      MODULE      PNMXPM
      SAVE        PNM_PM
      COMPLEX     PNM_PM
      ALLOCATABLE PNM_PM(:,:,:)
      END MODULE  PNMXPM
C=======================================================================
      MODULE      FRAGFL
      SAVE        F_FLAG
      ALLOCATABLE F_FLAG(:,:,:)
      END MODULE  FRAGFL
C=======================================================================
      MODULE      RFACTO
#if(SWITCH_QUAD==1)
      REAL*16
     *            FACTOR,    FACTOE,    FACTOO
#endif
      SAVE        FACTOR,    FACTOE,    FACTOO
      ALLOCATABLE FACTOR(:), FACTOE(:), FACTOO(:)
      END MODULE  RFACTO
C=======================================================================
C=======================================================================
C=======================================================================
C      MODULE      CANWFV    !SIMIMAR AS CANBAS: REMOVE IT in lipnp_2
C     Canonical wave functions and u- and v- amplitudes
C     WAVCAN(IBASE, ISTATE, IREVER, ICHARG)
C     V_CAN (       ISTATE, ICHARG)
C     U_CAN (       ISTATE, ICHARG)
C     Where IBASE  = index of basis state
C           ISTATE = index of can.  state
C           IREVER = index of time-reversal block
C           ICHARG = index of charge (n/p) block
C      SAVE        WAVCAN,          V_CAN,        U_CAN
C      COMPLEX     WAVCAN
C      REAL                         V_CAN,        U_CAN
C      ALLOCATABLE WAVCAN(:,:,:,:), V_CAN(:,:), U_CAN(:,:)
C      END MODULE  CANWFV
C=======================================================================
C=======================================================================
      MODULE      MAD2PP
      SAVE        DEN2PP
      COMPLEX     DEN2PP
      ALLOCATABLE DEN2PP(:,:,:,:)
      END MODULE  MAD2PP
C=======================================================================
C=======================================================================
      MODULE      MAD3PP
      SAVE        DEN3PP
      COMPLEX     DEN3PP
      ALLOCATABLE DEN3PP(:,:,:,:)
      END MODULE  MAD3PP
C=======================================================================
C=======================================================================
      MODULE      MAD4PP
      SAVE        DEN4PP
      COMPLEX     DEN4PP
      ALLOCATABLE DEN4PP(:,:,:,:)
      END MODULE  MAD4PP
C=======================================================================
C=======================================================================
      MODULE      MAD5PP
      SAVE        DEN5PP
      COMPLEX     DEN5PP
      ALLOCATABLE DEN5PP(:,:,:,:)
      END MODULE  MAD5PP
C=======================================================================
      MODULE      DENCON
      SAVE        DCONTI
      COMPLEX     DCONTI
      ALLOCATABLE DCONTI(:,:,:)
      END MODULE  DENCON
C=======================================================================
      MODULE      NEWCNT
      SAVE        ARASYM,ARDIST
      ALLOCATABLE ARASYM(:,:),ARDIST(:,:)
      END MODULE  NEWCNT
C=======================================================================
