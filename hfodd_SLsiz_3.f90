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
 !                  HFODD, ScaLAPACK inclusion module                   !
 !                                                                      !
 ! Authors: J.D. McDonnell, N. Schunck, ORNL                            !
 !                                                                      !
 ! This module incorporates ScaLAPACK capabilities into HFODD.          !
 !                                                                      !
 !                          First included in official release v249g    !
 !----------------------------------------------------------------------!

 MODULE hfodd_SLsiz

   COMPLEX, ALLOCATABLE, DIMENSION(:,:), PUBLIC:: Hsub, Evec

   REAL, ALLOCATABLE, DIMENSION(:), PUBLIC:: Eval

   ! submatrix stored on the processor (HFB, eigenvectors... eigenvalues),
   ! with the local/absolute array's dimensions
   INTEGER, PUBLIC:: larow, lacol, absrow, abscol

   ! block rows, block cols (local version of MBLOCK, NBLOCK)
   INTEGER, PUBLIC:: MB, NB

   ! BLACS context, shape of process gridx2, row/col of processx2,
   ! communication subroutine (error) info
   INTEGER, PUBLIC:: ICTXT, NPROW, NPCOL, MYROW, MYCOL, icom_err

   ! The ScaLAPACK descriptor for HFB matrix
   INTEGER, DIMENSION(9), PUBLIC:: DESCHFB, DESCEIG

 CONTAINS

   !---------------------------------------------------------------------!
   !    Initialization of the ScaLAPACK environment and processor grid   !
   !---------------------------------------------------------------------!
   SUBROUTINE SL_INIT(inrow, incol, color, tribeRank, numberHFODDproc, rankTribe)

      INTEGER, INTENT(IN) :: inrow, incol, color, tribeRank, numberHFODDproc
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: rankTribe

      INTEGER, DIMENSION(:,:), ALLOCATABLE :: gridProc

      INTEGER :: i_row, j_col, i_tribe

      ALLOCATE(gridProc(inrow,incol))

      IF ( incol*inrow > numberHFODDproc ) STOP 'Error in SL_INIT'

      ! Advanced : the BLACS processor grid is defined within a MPI group of processors.
      ! feature    This requires to send to every processor the id of all other processors
      !            part of the same group (done at the beginning of HFODD). By default, a
      !            processor only has access to its own id number, via a call to mpi_comm_rank()
      i_tribe = 0
      DO j_col = 1, incol
         DO i_row = 1, inrow
            gridProc(i_row, j_col) = rankTribe(i_tribe, color)
            i_tribe = i_tribe + 1
         END DO
      END DO

      CALL BLACS_GET(0, 0, ICTXT)

      CALL BLACS_GRIDMAP(ICTXT, gridProc, inrow, inrow, incol)

      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)

   END SUBROUTINE SL_INIT

   !---------------------------------------------------------------------!
   !    Initialization of the ScaLAPACK context and local matrix sizes   !
   !---------------------------------------------------------------------!
   SUBROUTINE initSL_siz(ar,ac,lr,lc,inMB,inNB)

      INTEGER, INTENT(IN) :: ar, ac, lr, lc, inMB, inNB

      MB = inMB
      NB = inNB

      CALL DESCINIT(DESCHFB, ar, ac, MB, NB, 0, 0, ICTXT, lr, icom_err)
      CALL DESCINIT(DESCEIG, ar, ac, MB, NB, 0, 0, ICTXT, lr, icom_err)

      absrow = ar
      abscol = ac
      larow = lr
      lacol = lc

   END SUBROUTINE initSL_siz

   !---------------------------------------------------------------------!
   !  Parallel diagonalization of a complex Hermitian matrix Ham  using  !
   !  ScaLAPACK routines                                                 !
   !---------------------------------------------------------------------!
   SUBROUTINE diaSL_HFB_siz(Ham,KZHPEV)

     INTEGER, INTENT(IN) :: KZHPEV
     COMPLEX, POINTER, DIMENSION(:,:) :: Ham

     COMPLEX, ALLOCATABLE, DIMENSION(:) :: WORK, RWORK

     REAL :: realTime

     INTEGER :: LWORK, LRWORK, LIWORK
     INTEGER :: NP0, NQ0
     INTEGER :: IH, JH, IEv, JEv

     INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK

     NP0 = NUMROC(absrow, NB, MYROW, 0, NPROW)
     NQ0 = NUMROC(abscol, NB, MYCOL, 0, NPCOL)

        ! Definitions of sizes for work arrays

     ! Simple driver
     IF (KZHPEV.EQ.1) THEN
         LRWORK = 2*absrow + 2*absrow - 2
         LWORK  = (NP0 + NQ0 + NB)*NB + 3*absrow + (absrow**2)
     END IF

        ! Divide-and-Conquer driver
     IF (KZHPEV.EQ.4) THEN
         MP0 = NUMROC(absrow, NB, 0, 0, NPROW)
         MQ0 = NUMROC(abscol, NB, 0, 0, NPCOL)
         LRWORK = 1 + 9*absrow + 3*NP0*NQ0
         LWORK  = (MP0 + MQ0 + NB)*NB + absrow
     END IF

     IF (.NOT. ALLOCATED(Evec)) ALLOCATE(Evec(larow,lacol))
     IF (.NOT. ALLOCATED(Eval)) ALLOCATE(Eval(abscol))
     IF (.NOT. ALLOCATED(WORK)) ALLOCATE(WORK(LWORK), RWORK(LRWORK))

     IH  = 1
     JH  = 1
     IEv = 1
     JEv = 1

     CALL cpu_time(realTime)

     ! Diagonalization

     ! Simple driver
     IF (KZHPEV.EQ.1) THEN
         CALL PZHEEV('V', 'L', absrow, Ham, IH, JH, DESCHFB, Eval, Evec, IEv, JEv, &
                             DESCEIG, WORK, LWORK, RWORK, LRWORK, icom_err)
         CALL cpu_time(realTime)
     END IF

        ! Divide-and-Conquer driver
     IF (KZHPEV.EQ.4) THEN
         LIWORK=7*absrow + 8*NPCOL + 2
         ALLOCATE(IWORK(LIWORK))

         CALL PZHEEVD('V', 'L', absrow, Ham, IH, JH, DESCHFB, Eval, Evec, IEv, JEv, &
                      DESCEIG, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, icom_err)
         DEALLOCATE(IWORK)
         call cpu_time(realTime)
     END IF

     IF (ALLOCATED(WORK)) DEALLOCATE(WORK)
     IF (ALLOCATED(RWORK)) DEALLOCATE(RWORK)

   END SUBROUTINE diaSL_HFB_siz

   !---------------------------------------------------------------------!
   !  Reconstruction of the total eigenvectors and eigenvalues by using  !
   !  BLACS communication routines. EWAHFB and EQUHFB are global arrays, !
   !  i.e. identical arrays in all processes, which contain all eigen-   !
   !  vectors and values of the HFB matrix. Evec contains the local part !
   !  of the eigenvectors (local=in one process only), Eval the local    !
   !  part of the eigenvalues.                                           !
   !---------------------------------------------------------------------!
   SUBROUTINE diaSL_collapse_HFB_siz(EWAHFB,EQUHFB,bigs)

     INTEGER, INTENT(IN) :: bigs

     COMPLEX, DIMENSION(1:bigs, 1:bigs) :: EWAHFB
     COMPLEX, DIMENSION(:,:), ALLOCATABLE :: tEv

     REAL, DIMENSION(1:bigs) :: EQUHFB

     INTEGER, DIMENSION(9):: DESC_Ve

     IF (.NOT. ALLOCATED(tEv)) ALLOCATE(tEv(absrow,abscol))
     tEv = CMPLX(0.0d0,0.0d0)

     CALL DESCINIT(DESC_Ve, absrow, abscol, absrow, abscol, 0, 0, ICTXT, absrow, icom_err)

     CALL PZGEMR2D(absrow, abscol, Evec, 1, 1, DESCEIG, tEv, 1, 1, DESC_Ve, ICTXT)

     IF (MYROW .eq. 0 .and. MYCOL .eq. 0) THEN
        CALL ZGEBS2D(ICTXT, 'A', ' ', absrow, abscol, tEv, absrow)
     ELSE
        CALL ZGEBR2D(ICTXT, 'A', ' ', absrow, abscol, tEv, absrow, 0, 0)
     ENDIF

     EWAHFB(1:absrow,1:abscol) = tEv
     EQUHFB(1:absrow) = Eval(1:absrow)

     ! Synchronization at the end, just for security
     CALL BLACS_BARRIER( ICTXT, 'A' )

     IF (ALLOCATED(Eval)) DEALLOCATE(Eval)
     IF (ALLOCATED(Evec)) DEALLOCATE(Evec)
     IF (ALLOCATED(tEv))  DEALLOCATE(tEv)

   END SUBROUTINE diaSL_collapse_HFB_siz

   !---------------------------------------------------------------------!
   ! Closing down the ScaLAPACK environment and releasing processor grid !
   !---------------------------------------------------------------------!
   SUBROUTINE closeSL_all()
     CALL BLACS_GRIDEXIT(ICTXT)
     CALL BLACS_EXIT(0)
   END SUBROUTINE closeSL_all

END MODULE hfodd_SLsiz


