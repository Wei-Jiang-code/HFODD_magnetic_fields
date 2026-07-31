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

 !----------------------------------------------------------------------!
 !                                                                      !
 !                Energy Density Functional for HFODD                   !
 !                                                                      !
 ! Authors: M. Stoitsov, M. Kortelainen, N. Schunck, ORNL               !
 !                                                                      !
 ! This module is/will be in charge of:                                 !
 !    - defining the coupling constants of a given EDF                  !
 !    - converting back/forth from nuclear matter properties to (t,x)   !
 !      to coupling constant representation.                            !
 !                                                                      !
 !                          First included in official release v249m    !
 !                                                                      !
 !----------------------------------------------------------------------!

Module hfodd_functional
  Use HFBTHO_utilities
  Use UNEDF
  Implicit None

Contains

   ! returns (t,x) filled out from definition of EDF
   Subroutine paramsFunctional(SKYRME,IETA_J,IETA_W,IETA_M,IETACM,&
                               ISKYRM,HBMASS,T0_EFF,T1_EFF,T2_EFF,&
                               T3_EFF,X0_EFF,X1_EFF,X2_EFF,X3_EFF,&
                               B4_EFF,B4PEFF,te_EFF,to_EFF,POWEFF,&
                                             V0NEFF,V0PEFF,COUSCA)

     Integer(ipr), INTENT(INOUT) :: IETA_J,IETA_W,IETA_M,IETACM,ISKYRM
     Real(pr), INTENT(INOUT) :: HBMASS,T0_EFF,T1_EFF,T2_EFF,T3_EFF, &
                                       X0_EFF,X1_EFF,X2_EFF,X3_EFF, &
                                       B4_EFF,B4PEFF,te_EFF,to_EFF, &
                                       POWEFF,V0NEFF,V0PEFF,COUSCA
     Character(Len=4), INTENT(IN) :: SKYRME
     Character(Len=30) :: fname
     Integer(ipr) :: noForces

     fname    = SKYRME
     noForces = 0

     ! Set coupling constants and nuclear matter properties
     Call read_UNEDF_NAMELIST(fname,noForces)

     ! Computes bulk coupling constants
     Call set_functional_parameters(fname,.True.)

     ! Extract (t,x) parameterization from coupling constants
     Call t_from_C()

     ! Export (t,x) parameterization into HFODD
     T0_EFF = t0_pub
     T1_EFF = t1_pub
     T2_EFF = t2_pub
     T3_EFF = t3_pub

     X0_EFF = x0_pub
     X1_EFF = x1_pub
     X2_EFF = x2_pub
     X3_EFF = x3_pub

     B4_EFF = b4_pub
     B4PEFF = b4p_pub

     te_EFF = te_pub
     to_EFF = to_pub

     POWEFF = sigma

     V0NEFF = CpV0(0)
     V0PEFF = CpV0(1)

     COUSCA = CExPar

     HBMASS = hbzero

                     IETA_J=0
     If(use_j2terms) IETA_J=1

     IETA_M=1
     IETA_W=1
                         IETACM=3
     If(use_cm_cor)      IETACM=0
     If(use_full_cm_cor) IETACM=1

  End Subroutine paramsFunctional

End Module hfodd_functional
