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

 MODULE hfodd_sizes

      IMPLICIT NONE

      !   Linear Algebra
      INTEGER, PARAMETER :: I_CRAY=0,IPARAL=0,IPARAT=0
      !
      !   General
      INTEGER, PARAMETER :: NDKART=3
      INTEGER, PARAMETER :: NDPARI=1, NDSIMP=1, NDSIGN=1, NDREVE=1
      INTEGER, PARAMETER :: NDISOS=1, NDSPIN=1, NDRHFB=1
      INTEGER, PARAMETER :: NDPAHO=5
      INTEGER, PARAMETER :: ND04=2, ND08=4, ND12=6, ND16=8, ND32=16, ND64=32, NDMAX=2*ND64
      INTEGER, PARAMETER :: NDENES=43
      !
      !   HO basis: number of states, integration mesh
      INTEGER, PARAMETER :: NDMAIN=22, NDXMAX=NDMAIN, NDYMAX=NDMAIN, NDZMAX=NDMAIN
      INTEGER, PARAMETER :: NDOSCI=NDMAIN, ND0MAX=2*NDOSCI, ND1MAX=NDOSCI+2
      INTEGER, PARAMETER :: ND2MAX=2*NDOSCI+2, NDBTOT=(NDOSCI+1)**3
      INTEGER, PARAMETER :: NDINTE=2*NDOSCI+3, NDINT1=NDOSCI+1, NDINT2=NDOSCI+2
      INTEGER, PARAMETER :: NDPERP=NDMAIN
      INTEGER, PARAMETER :: NDTRAN=((NDMAIN+1)*(NDMAIN+2)*(NDMAIN+3))/6
      INTEGER, PARAMETER :: NDBASE=680, NDBASD=4*NDBASE
      INTEGER, PARAMETER :: NDBASP=NDBASE/2+MAX(20,NDBASE/12)
      INTEGER, PARAMETER :: NDSTAT=680
      INTEGER, PARAMETER :: NDDATA=2*NDSTAT*(NDISOS+1)
      INTEGER, PARAMETER :: NDHALF=2, NDREIM=2
      INTEGER, PARAMETER :: NDXHRM=34, NDYHRM=34, NDZHRM=46
      INTEGER, PARAMETER :: NDGAUS=MAX(NDXHRM,NDYHRM,NDZHRM)
      !
      !   Densities
      INTEGER, PARAMETER :: NDALLD=4+7*NDKART+2*NDKART*NDKART, NDALLP=1
      !INTEGER, PARAMETER :: NDALLP=4+NDKART+NDKART*NDKART
      INTEGER, PARAMETER :: NDTYDE=2
      INTEGER, PARAMETER :: NDPNMX=1
      !
      !   Various: multipoles, Coulomb, iterations, binomial coefficients,
      !            functional
      INTEGER, PARAMETER :: NDMULT=9, NDLAMB=9, NDFACT=2*NDLAMB, NDCONS=10
      INTEGER, PARAMETER :: NDMULR=4, NDMULM=(NDMULR+1)**2
      INTEGER, PARAMETER :: NDDERI=2, NDSCAL=1
      INTEGER, PARAMETER :: NDCOUL=80, NDPOLS=25, NDNETA=NDCOUL
      INTEGER, PARAMETER :: ND_COU=20, NDCOTY=7
      INTEGER, PARAMETER :: NDYUKA=6, NDYUTY=2
      INTEGER, PARAMETER :: NDGOGA=2, NDFORC=4, NDREGA=4
      INTEGER, PARAMETER :: NDITER=5000, NDTERM=15
      INTEGER, PARAMETER :: NDSHIF=100, NDSMTH=20
      INTEGER, PARAMETER :: NDBINO=2*NDOSCI
      INTEGER, PARAMETER :: NDFACR=MAX(2*NDOSCI+7,NDGAUS)
      INTEGER, PARAMETER :: NDCOUP=24, NDVIOL=32
      INTEGER, PARAMETER :: NDJMAX=70, NDJMXD=2*NDJMAX
      INTEGER, PARAMETER :: NFIRST=1, NLAST=40000
      INTEGER, PARAMETER :: NDSUBR=100, NDCOLU=4
      INTEGER, PARAMETER :: MV25=150
      INTEGER, PARAMETER :: NDNUCL=1
      INTEGER, PARAMETER :: NDSEAR=10
      INTEGER, PARAMETER :: ND_PRO=140, ND_NEU=250
      !
      !    Projection
      INTEGER, PARAMETER :: NDAKNO=1, NDBKNO=1
      INTEGER, PARAMETER :: NDASAV=(1-IPARAL)*NDAKNO+IPARAL
      INTEGER, PARAMETER :: NDPROI=20, NDPROT=10, NDPROD=(NDPROI/2+1)*(NDPROT/2+1)+1
      INTEGER, PARAMETER :: NDPROK=(NDPROI/2+1)**2
      INTEGER, PARAMETER :: NDPROM=((NDPROI+1)*(NDPROI+2)*(NDPROI+3))/6
      INTEGER, PARAMETER :: NUMCOL=70
      INTEGER, PARAMETER :: NDATKN=1,NDBTKN=10
      INTEGER, PARAMETER :: NDASAT=(1-IPARAT)*NDATKN+IPARAT
      !INTEGER, PARAMETER :: NDISOM=1771)
      INTEGER, PARAMETER :: NDISOM=((NDPROT+1)*(NDPROT+2)*(NDPROT+3))/6
      !=========================================================================
      ! INTEGER, PARAMETER NDREDU BELOW IS EQUAL TO THE MINIMUM NECESSARY VALUE
      !=========================================================================
      ! INTEGER, PARAMETER
      !*          (NDREDU=NDPROK*(NDMULR+1)*NDPROK)
      !=========================================================================
      ! INTEGER, PARAMETER NDREDU BELOW IS EQUAL TO THE MINIMUM SUFFICIENT VALUE
      ! WHICH HAS BEEN OBTAINED BY FITTING THE  VALUES  CALCULATED  FOR
      ! NDMULR=4 WITH A 3RD-ORDER POLYNOMIAL, SEE THE ANALYSIS  IN  THE
      ! FILE HF222K.XLS OF 10/02/06
      !=========================================================================
      INTEGER, PARAMETER :: NDREDU=5*((((NDPROI+2)*NDPROI-5)*NDPROI+32))/2

 END MODULE hfodd_sizes
