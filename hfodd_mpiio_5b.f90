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
 !                MPI I/O Environment for HFODD                         !
 !                                                                      !
 ! Authors: N. Schunck, ORNL                                            !
 !                                                                      !
 ! This module is/will be in charge of:                                 !
 !    - reading sequential and parallel input data from the hfodd.d     !
 !      and hfodd_mpiio.d files                                         !
 !    - broadcasting all input data to all processors                   !
 !    - converting data to the proper variables used in HFODD           !
 !                                                                      !
 !                          First included in official release v249e    !
 !                                                                      !
 !----------------------------------------------------------------------!

 Module linkedLists

    Implicit None

    ! Define derived types Node_int and Node_dbl, which are used
    ! to construct linked list of integers and real numbers

    Type :: Node_int
        Integer           :: value
        Integer           :: numero
        Character(Len=10) :: name
        TYPE(Node_int), pointer  :: next
    End type Node_int

    Type :: Node_dbl
        Real              :: value
        Integer           :: numero
        Character(Len=10) :: name
        TYPE(Node_dbl), pointer  :: next
    End type Node_dbl

    Type :: Node_str
        Character(Len=10) :: value
        Integer           :: numero
        Character(Len=10) :: name
        TYPE(Node_str), pointer  :: next
    End type Node_str

 Contains

    ! A bunch of useful subroutines to manipulate linked lists

    ! Routines dealing with Node_int type
    Recursive Subroutine printEntireList_int(ptr)
       Type(Node_int), pointer :: ptr
       If(associated(ptr%next)) Call printEntireList_int(ptr%next)
       print *, ptr%value, ptr%numero, ptr%name
    End Subroutine printEntireList_int

    Recursive Subroutine freeEntireList_int(ptr)
       Type(Node_int), pointer :: ptr
       If(associated(ptr%next)) Call freeEntireList_int(ptr%next)
       Deallocate(ptr)
    End Subroutine freeEntireList_int

    ! Routines dealing with Node_dbl type
    Recursive Subroutine printEntireList_dbl(ptr)
       Type(Node_dbl), pointer :: ptr
       If(associated(ptr%next)) Call printEntireList_dbl(ptr%next)
       Print *, ptr%value, ptr%numero, ptr%name
    End Subroutine printEntireList_dbl

    Recursive Subroutine freeEntireList_dbl(ptr)
       Type(Node_dbl), pointer :: ptr
       If(associated(ptr%next)) Call freeEntireList_dbl(ptr%next)
       Deallocate(ptr)
    End Subroutine freeEntireList_dbl

    ! Routines dealing with Node_str type
    Recursive Subroutine printEntireList_str(ptr)
       Type(Node_str), pointer :: ptr
       If(associated(ptr%next)) Call printEntireList_str(ptr%next)
       Print *, ptr%value, ptr%numero, ptr%name
    End Subroutine printEntireList_str

    Recursive Subroutine freeEntireList_str(ptr)
       Type(Node_str), pointer :: ptr
       If(associated(ptr%next)) Call freeEntireList_str(ptr%next)
       Deallocate(ptr)
    End Subroutine freeEntireList_str

 End Module linkedLists

 Module hfodd_mpiio

    Use linkedLists

    Implicit None

    Include "mpif.h"

    ! Private
    Integer, PRIVATE :: stringLength = 10
    Integer, PRIVATE :: counterInt, counterDbl, counterStr, counterTot
    Integer, PRIVATE :: counterVector(1:3)
    Integer, Allocatable, PRIVATE :: vectorInt_value(:)
    Real, Allocatable, PRIVATE :: vectorDbl_value(:)
    Character(Len=10), Allocatable, PRIVATE :: vectorStr_value(:)
    Character(Len=10), Allocatable, PRIVATE :: vectorAll_names(:)

    Type(Node_int), pointer :: listeInt
    Type(Node_dbl), pointer :: listeDbl
    Type(Node_str), pointer :: listeStr

    ! Public
    Integer, PUBLIC :: modeMPI_deformation,modeMPI_basis
    Integer, PUBLIC :: numberConstraints,numberConstraints_restart,counterForce, &
                       batch_mode,number_batch,number_constraint_batch
    Integer, PUBLIC :: neutron_stepNumber,neutron_stepValue,neutron_initValue,                &
                       proton_stepNumber,proton_stepValue,proton_initValue,                   &
                       basis_frequency_stepNumber,basis_deform_stepNumber,                    &
                       basis_size_stepNumber,basis_size_stepValue,basis_size_initValue,       &
                       basis_states_stepNumber,basis_states_stepValue,basis_states_initValue, &
                       turnConstraintOn,optimizePath
    Integer, Allocatable, PUBLIC :: qLambda(:),qMiu(:),qStepNumber(:),qLambda_restart(:),  &
                                    qMiu_restart(:),qStepNumber_restart(:),lambda_batch(:),&
                                    miu_batch(:)

    Real, PUBLIC :: x_neutron_initValue,x_neutron_stepValue, &
                    x_proton_initValue,x_proton_stepValue
    Real, PUBLIC :: basis_frequency_stepValue,basis_frequency_initValue, &
                    basis_deform_stepValue,basis_deform_initValue
    Real, Allocatable, PUBLIC :: qinitValue(:),qFinValue(:),qinitValue_restart(:), &
                                 qFinValue_restart(:)

    Character(Len=4), Allocatable, PUBLIC :: forceVector(:)

 Contains

    !------------------------------------------------------------------!
    ! This subroutine reads the standard HFODD input file and puts     !
    ! all encountered data into two linked lists of type Node_int      !
    ! (for data of type integer) and Node_dbl (for data of type real). !
    ! These derived types also contain a character string that is      !
    ! used to identify the keyword where this data come from. This     !
    ! will serve for the reconstruction of the original variables.     !
    !------------------------------------------------------------------!

    Subroutine mpi_getSequentialData()

      Integer :: NFIREA

      Integer, Parameter :: NDITER = 5000

      Integer :: ii, jj, counterDef, counterWS, counterColl
      Integer :: mpi_rank, mpi_err, mpi_string_10

      Integer :: ICOTYP, ICOUDI, ICOUEX, I_GOGA, NOITER, IBROYD, N_DUMM, &
                 NO_RUN, NULAST, NUPING, NUCHAO, ISIMPY, ISIGNY, IPARTY, &
                 IROTAT, ISIMTX, ISIMTY, ISIMTZ, IPAIRI, IPAHFB, IPABCS, &
                 IMFHFB, LIPKIN, LIPKIP, IOPTGS, NXHERM, NYHERM, NZHERM, &
                 NOSCIL, NLIMIT, I1LINE, NILXYZ, IREVIE, ICONTI, IPCONT, &
                 IGCONT, IGPCON, ILCONT, IFCONT, IACONT, IFSHEL, NPOLYN, &
                 LAMBDA, ITE_WS, IREAWS, ISTAND, KETA_J, KETA_W, KETACM, &
                 KETA_M, IF_RPA, IF_EDF, IF_THO, IBASIS, IHFPOT, IFRAGM, &
                 L_PUSH, M_PUSH, N_PUSH, DOPUSH, NPLAST, MIU,    NMUCON, &
                 NMUCOU, NMUPRI, NEXBET, IPRIBE, IPRIBL, IWRIFI, INNUMB, &
                 IZNUMB, IGOGPA, ITRMAX, NTHETA, MIN_QP, IDEALL, IDELOC, &
                 IDECON, IF2FIN, IF2FIP

      Real :: A_DUMM, B_DUMM, ENECUT, SLOWEV, SLOWOD, SLOWPA, &
              SLOWLI, EPSITE, EPSPNG, ECUTOF, PRHO_T, PRHOST, &
              PRHO_N, PRHOSN, POWERN, PRHO_P, PRHOSP, POWERP, &
              POWERT, OMEGAY, GSTRUN, GSTRUP, HOMFAC, ALPHAR, &
              STP_WS, CBETHO, TEMP_T, Q_PUSH, FCHOM0, R0PARM, &
              DELTAE, XLOCMX, V2_MIN, FE2FIN, FE2FIP

      Real :: PI, DEGRAD

      Character(Len=10) :: KEYWOR

      Type(Node_int), Pointer :: tempNodeInt, currentNodeInt, previousNodeInt
      Type(Node_dbl), Pointer :: tempNodeDbl, currentNodeDbl, previousNodeDbl
      Type(Node_str), Pointer :: tempNodeStr, currentNodeStr, previousNodeStr

      ! Initialize local variables

      ICOTYP=7
      ICOUDI=1;ICOUEX=1
      I_GOGA=0;IGOGPA=0
      NOITER=10
      IBROYD=1
      N_DUMM=8
      NO_RUN=0
      NULAST=3;NUPING=3;NUCHAO=3
      ISIMPY=1;ISIGNY=1;IPARTY=1
      IROTAT=0
      ISIMTX=-1;ISIMTY=-1;ISIMTZ=-1
      IPAIRI=1;IPAHFB=1
      IPABCS=0
      IMFHFB=0;LIPKIN=0;LIPKIP=0;IF2FIN=0;IF2FIP=0
      IOPTGS=1
      NXHERM=30;NYHERM=30;NZHERM=30
      NOSCIL=12
      NLIMIT=455
      INNUMB=86;IZNUMB=66
      LAMBDA=0; MIU=0
      I1LINE=1
      NILXYZ=3
      IREVIE=0;ICONTI=0;IPCONT=0;IGCONT=0;IGPCON=0;ILCONT=0;IFCONT=0;IACONT=0
      IFSHEL=0;NPOLYN=6
      ITE_WS=10;IREAWS=0
      ISTAND=1;KETA_J=0;KETA_W=1;KETACM=1;KETA_M=0
      IF_RPA=1;IF_EDF=0;IF_THO=0;CBETHO=0.0D0
      IBASIS=0
      IFRAGM=0
      ITRMAX=31;NTHETA=181
      IDEALL=0;IDELOC=0;IDECON=0
      IHFPOT=0
      L_PUSH=1;M_PUSH=0;N_PUSH=0;DOPUSH=0;NPLAST=1
      NMUCON=4;NMUCOU=4;NMUPRI=4
      NEXBET=4;IPRIBE=1;IPRIBL=0
      IWRIFI=0

      A_DUMM=0.0D0; B_DUMM=0.0D0
      ENECUT=800.0D0
      R0PARM=1.2D0
      ALPHAR=0.0D0
      SLOWEV=0.5D0; SLOWOD=0.5D0; SLOWPA=0.5D0; SLOWLI=0.5D0
      EPSITE=1.0D-5;EPSPNG=1.0d-5
      ECUTOF=60.0D0
      PRHO_T=-000.0D0; PRHOST=0.32D0; POWERT=1.0D0
      PRHO_N=-200.0D0; PRHOSN=0.32D0; POWERN=1.0D0
      PRHO_P=-200.0D0; PRHOSP=0.32D0; POWERP=1.0D0
      OMEGAY=0.0D0
      GSTRUN=1.2; GSTRUP=1.2;HOMFAC=3.5D0
      STP_WS=0.01D0
      Q_PUSH=0.00D0
      FCHOM0=1.20D0
      MIN_QP=0;DELTAE=2.50D0;XLOCMX=0.75D0;V2_MIN=0.005D0
      FE2FIN=0.1D0;FE2FIP=0.1D0
      TEMP_T=0.0D0

      ! Getting the rank of the current process
      Call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, mpi_err)

      ! Define MPI type for a character string of length stringLength and commit
      Call mpi_type_contiguous(stringLength, MPI_CHARACTER, mpi_string_10, mpi_err)
      Call mpi_type_commit(mpi_string_10, mpi_err)

      NO_RUN = 0

      ! Only the root process reads the data
      If(mpi_rank.Eq.0) Then

         PI=4.0D0*Atan(1.0D0)
         DEGRAD=180.0D0/PI

         counterInt = 0
         counterDbl = 0
         counterStr = 0

         ! Initializing the end of the linked list
         Write(KEYWOR,'(10x)')

         counterInt = counterInt + 1
         Allocate(listeInt)
         listeInt%value  = 0
         listeInt%numero = counterInt
         listeInt%name   = KEYWOR
         NULLIFY(listeInt%next)

         counterDbl = counterDbl + 1
         Allocate(listeDbl)
         listeDbl%value  = 0.0D0
         listeDbl%numero = counterDbl
         listeDbl%name   = KEYWOR
         NULLIFY(listeDbl%next)

         counterStr = counterStr + 1
         Allocate(listeStr)
         listeStr%value  = "INITIALIZE"
         listeStr%numero = counterStr
         listeStr%name   = KEYWOR
         NULLIFY(listeStr%next)

         NFIREA=87
         Open(NFIREA,FILE='hfodd.d',STATUS='UNKNOWN',FORM='FORMATTED')

         ! Main loop over the keyword value
         Do While (trim(KEYWOR) .NE. 'ALL_DONE' .AND. trim(KEYWOR) .NE. 'ALL-DONE')

            Read(NFIREA,'(A10)') KEYWOR

            !==============================!
            !   ITERATIONS CONTROL FLAGS   !
            !==============================!

            ! Number of iterations
            If(trim(KEYWOR).Eq.'ITERATIONS') Then

               Read(NFIREA,*) NOITER

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = NOITER
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               If(NOITER.GT.NDITER) STOP 'NOITER.GT.NDITER IN NAMELI'

            End If

            ! Switches on/off the Broyden method
            If(trim(KEYWOR).Eq.'BROYDEN') Then

               Read(NFIREA,*) IBROYD,N_DUMM,A_DUMM,B_DUMM

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IBROYD
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = N_DUMM
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = A_DUMM
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = B_DUMM
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

            End If

            ! Slow-down parameters for the mean-field
            If(trim(KEYWOR).Eq.'SLOW_DOWN'.OR.trim(KEYWOR).Eq.'SLOW-DOWN') Then

               Read(NFIREA,*) SLOWEV,SLOWOD

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = SLOWEV
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = SLOWOD
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

            End If

            ! Slow-down parameters for the pairing field
            If(trim(KEYWOR).Eq.'SLOW_PAIR'.OR.trim(KEYWOR).Eq.'SLOW-PAIR') Then

               Read(NFIREA,*) SLOWPA

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = SLOWPA
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

            End If

            ! Slow-down parameters for the Lipkin-Nogami correction
            If(trim(KEYWOR).Eq.'SLOWLIPKIN') Then

               Read(NFIREA,*) SLOWLI

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = SLOWLI
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

            End If

            ! Defines the precision of convergence
            If(trim(KEYWOR).Eq.'ITERAT_EPS'.OR.trim(KEYWOR).Eq.'ITERAT-EPS') Then

               Read(NFIREA,*) EPSITE

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = EPSITE
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

            End If

            ! Defines divergence
            If(trim(KEYWOR).Eq.'MAXANTIOSC') Then

               Read(NFIREA,*) NULAST

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = NULAST
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Control ping-pong effects
            If(trim(KEYWOR).Eq.'PING_PONG'.OR.trim(KEYWOR).Eq.'PING-PONG') Then

               Read(NFIREA,*) EPSPNG,NUPING

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = NUPING
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = EPSPNG
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

            End If

            ! Control chaotic divergence
            If(trim(KEYWOR).Eq.'CHAOTIC') Then

               Read(NFIREA,*) NUCHAO

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = NUCHAO
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            !=====================================!
            !   SPECIFIC FEATURES CONTROL FLAGS   !
            !=====================================!

            ! Temperature
            If(trim(KEYWOR).Eq.'FINITETEMP') Then

               Read(NFIREA,*) TEMP_T

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = TEMP_T
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

            End If

            ! Definition of the collective space for calculation of collective inertia
            If(trim(KEYWOR).Eq.'COLL_SPACE'.OR.trim(KEYWOR).Eq.'COLL-SPACE') Then

               counterColl = 0
               LAMBDA=-1

               Do While (LAMBDA < 0)

                  Read(NFIREA,*) LAMBDA,MIU

                  counterColl = counterColl + 1

                  counterInt = counterInt + 1
                  Allocate(tempNodeInt)
                  tempNodeInt%value  = ABS(LAMBDA)
                  tempNodeInt%numero = counterInt
                  tempNodeInt%name   = KEYWOR
                  tempNodeInt%next   => listeInt
                  listeInt => tempNodeInt

                  counterInt = counterInt + 1
                  Allocate(tempNodeInt)
                  tempNodeInt%value  = MIU
                  tempNodeInt%numero = counterInt
                  tempNodeInt%name   = KEYWOR
                  tempNodeInt%next   => listeInt
                  listeInt => tempNodeInt

               End Do

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = counterColl
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Type of shell correction
            If(trim(KEYWOR).Eq.'SHELLCORCT') Then

               Read(NFIREA,*) IFSHEL

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IFSHEL
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Shell correction parameters
            If(trim(KEYWOR).Eq.'SHELLPARAM') Then

               Read(NFIREA,*) GSTRUN,GSTRUP,HOMFAC,NPOLYN

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = GSTRUN
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = GSTRUP
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = HOMFAC
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = NPOLYN
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Switching to interface with HFBTHO
            If(trim(KEYWOR).Eq.'HFBTHOISON') Then

               Read(NFIREA,*) IF_THO,CBETHO

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IF_THO
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = CBETHO
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

            End If

            ! Compute mass and charge of fission fragments
            If(trim(KEYWOR).Eq.'MASSFRAGME') Then

               Read(NFIREA,*) IFRAGM

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IFRAGM
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Compute and record densities of the fission fragments
            If(trim(KEYWOR).Eq.'FRAGDENSIT') Then

               Read(NFIREA,*) IDEALL,IDELOC,IDECON

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IDEALL
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IDELOC
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IDECON
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Rotate selected pairs of q.p. wave-functions
            If(trim(KEYWOR).Eq.'QPROTATION') Then

               Read(NFIREA,*) MIN_QP,DELTAE,XLOCMX,V2_MIN,ITRMAX,NTHETA

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = MIN_QP
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = DELTAE
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = XLOCMX
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = V2_MIN
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = ITRMAX
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = NTHETA
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            !=====================================!
            !      INTERACTIONS CONTROL FLAGS     !
            !=====================================!

            ! Switching to Energy Density Functional formalism
            If(trim(KEYWOR).Eq.'UNEDF_PROJ'.OR.trim(KEYWOR).Eq.'UNEDF-PROJ') Then

               Read(NFIREA,*) IF_EDF

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IF_EDF
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Skyrme force flags
            If(trim(KEYWOR).Eq.'SKYRME_STD'.OR.trim(KEYWOR).Eq.'SKYRME-STD') Then

               Read(NFIREA,*) ISTAND,KETA_J,KETA_W,KETACM,KETA_M

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = ISTAND
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = KETA_J
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = KETA_W
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = KETACM
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = KETA_M
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Type of Coulomb potential
            If(trim(KEYWOR).Eq.'COULOMBPAR') Then

               Read(NFIREA,*) ICOTYP,ICOUDI,ICOUEX

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = ICOTYP
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = ICOUDI
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = ICOUEX
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Gogny force
            If(trim(KEYWOR).Eq.'GOGNY') Then

               Read(NFIREA,*) I_GOGA

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = I_GOGA
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            !==============================!
            !     PAIRING CONTROL FLAGS    !
            !==============================!

            ! Gogny pairing
            If(trim(KEYWOR).Eq.'GOGNY_PAIR') Then

               Read(NFIREA,*) IGOGPA

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IGOGPA
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Switches on/off pairing
            If(trim(KEYWOR).Eq.'PAIRING') Then

               Read(NFIREA,*) IPAIRI

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IPAIRI
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! HFB pairing
            If(trim(KEYWOR).Eq.'HFB') Then

               Read(NFIREA,*) IPAHFB

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IPAHFB
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Defines the pairing window cut-off
            If(trim(KEYWOR).Eq.'CUTOFF') Then

               Read(NFIREA,*) ECUTOF

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = ECUTOF
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

            End If

            ! Diagonalizes the mean-field (if HFB)
            If(trim(KEYWOR).Eq.'HFBMEANFLD') Then

               Read(NFIREA,*) IMFHFB

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IMFHFB
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Pairing interaction (neutrons only)
            If(trim(KEYWOR).Eq.'PAIRNINTER') Then

               Read(NFIREA,*) PRHO_N,PRHOSN,POWERN

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = PRHO_N
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = PRHOSN
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = POWERN
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

            End If

            ! Pairing interaction (protons only)
            If(trim(KEYWOR).Eq.'PAIRPINTER') Then

               Read(NFIREA,*) PRHO_P,PRHOSP,POWERP

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = PRHO_P
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = PRHOSP
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = POWERP
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

            End If

            ! Pairing interaction
            If(trim(KEYWOR).Eq.'PAIR_INTER') Then

               Read(NFIREA,*) PRHO_T,PRHOST,POWERT

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = PRHO_T
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = PRHOST
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = POWERT
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

            End If

            ! Lipkin-Nogami correction
            If(trim(KEYWOR).Eq.'LIPKIN') Then

               Read(NFIREA,*) LIPKIN,LIPKIP

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = LIPKIN
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = LIPKIP
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Reading the fixed Lipkin-Nogami lambda2 parameters
            If(trim(KEYWOR).Eq.'FIXLAMB2_N'.OR.trim(KEYWOR).EQ.'FIXLAMB2-N') Then

               Read(NFIREA,*) FE2FIN,IF2FIN

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = FE2FIN
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IF2FIN
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If
            If(trim(KEYWOR).Eq.'FIXLAMB2_P'.OR.trim(KEYWOR).EQ.'FIXLAMB2-P') Then

               Read(NFIREA,*) FE2FIP,IF2FIP

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = FE2FIP
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IF2FIP
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            !==============================!
            !     SYMMETRY CONTROL FLAGS   !
            !==============================!

            ! time-reversal
            If(trim(KEYWOR).Eq.'ROTATION') Then

               Read(NFIREA,*) IROTAT

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IROTAT
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! y-simplex
            If(trim(KEYWOR).Eq.'SIMPLEXY') Then

               Read(NFIREA,*) ISIMPY

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = ISIMPY
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! y-signature
            If(trim(KEYWOR).Eq.'SIGNATUREY') Then

               Read(NFIREA,*) ISIGNY

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = ISIGNY
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! parity
            If(trim(KEYWOR).Eq.'PARITY') Then

               Read(NFIREA,*) IPARTY

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IPARTY
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! All 3 simplexes at once
            If(trim(KEYWOR).Eq.'TSIMPLEX3D') Then

               Read(NFIREA,*) ISIMTX,ISIMTY,ISIMTZ

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = ISIMTX
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = ISIMTY
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = ISIMTZ
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            !=====================================!
            !       HO BASIS CONTROL FLAGS        !
            !=====================================!

            ! Automatic adjustment of basis deformation and frequency
            If(trim(KEYWOR).Eq.'BASISAUTOM') Then

               Read(NFIREA,*) IBASIS

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IBASIS
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Flag controling the basis frequency
            If(trim(KEYWOR).Eq.'HOMEGAZERO') Then

               Read(NFIREA,*) FCHOM0

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = FCHOM0
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

            End If

            ! Size of the HO basis
            If(trim(KEYWOR).Eq.'BASIS_SIZE'.OR.trim(KEYWOR).Eq.'BASIS-SIZE') Then

               Read(NFIREA,*) NOSCIL,NLIMIT,ENECUT

               If(NLIMIT.Eq.0) NLIMIT=((NOSCIL+1)*(NOSCIL+2)*(NOSCIL+3))/6

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = NOSCIL
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = NLIMIT
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = ENECUT
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

            End If

            ! Deformation of the HO basis
            If(trim(KEYWOR).Eq.'SURFAC_PAR'.OR.trim(KEYWOR).Eq.'SURFAC-PAR') Then

               Read(NFIREA,*) INNUMB,IZNUMB,R0PARM

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = INNUMB
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IZNUMB
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = R0PARM
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

            End If

            ! Deformation of the HO basis
            If(trim(KEYWOR).Eq.'SURFAC_DEF'.OR.trim(KEYWOR).Eq.'SURFAC-DEF') Then

               counterDef = 0
               LAMBDA=-1

               Do While (LAMBDA < 0)

                  Read(NFIREA,*) LAMBDA,MIU,ALPHAR

                  counterDef = counterDef + 1

                  counterInt = counterInt + 1
                  Allocate(tempNodeInt)
                  tempNodeInt%value  = ABS(LAMBDA)
                  tempNodeInt%numero = counterInt
                  tempNodeInt%name   = KEYWOR
                  tempNodeInt%next   => listeInt
                  listeInt => tempNodeInt

                  counterInt = counterInt + 1
                  Allocate(tempNodeInt)
                  tempNodeInt%value  = MIU
                  tempNodeInt%numero = counterInt
                  tempNodeInt%name   = KEYWOR
                  tempNodeInt%next   => listeInt
                  listeInt => tempNodeInt

                  counterDbl = counterDbl + 1
                  Allocate(tempNodeDbl)
                  tempNodeDbl%value  = ALPHAR
                  tempNodeDbl%numero = counterDbl
                  tempNodeDbl%name   = KEYWOR
                  tempNodeDbl%next   => listeDbl
                  listeDbl => tempNodeDbl

               End Do

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = counterDef
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Automatic number of Gauss-Hermite points
            If(trim(KEYWOR).Eq.'OPTI_GAUSS'.OR.trim(KEYWOR).Eq.'OPTI-GAUSS') Then

               Read(NFIREA,*) IOPTGS

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IOPTGS
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Number of Gauss-Hermite points
            If(trim(KEYWOR).Eq.'GAUSHERMIT') Then

               Read(NFIREA,*) NXHERM,NYHERM,NZHERM

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = NXHERM
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = NYHERM
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = NZHERM
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            !=====================================!
            !   MULTIPOLE MOMENTS CONTROL FLAGS   !
            !=====================================!

            ! Linear constraints based on the RPA matrix
            If(trim(KEYWOR).Eq.'RPA_CONSTR'.OR.trim(KEYWOR).Eq.'RPA-CONSTR') Then

               Read(NFIREA,*) IF_RPA

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IF_RPA
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Flags controlling the number of constraints
            If(trim(KEYWOR).Eq.'MAX_MULTIP'.OR.trim(KEYWOR).Eq.'MAX-MULTIP') Then

               Read(NFIREA,*) NMUCON,NMUCOU,NMUPRI

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = NMUCON
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = NMUCOU
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = NMUPRI
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Deformations of the WS potential used to provide initial solution
            If(trim(KEYWOR).Eq.'DEFORMS_WS'.OR.trim(KEYWOR).Eq.'DEFORMS-WS') Then

               counterWS = 0
               LAMBDA=-1

               Do While (LAMBDA < 0)

                  Read(NFIREA,*) LAMBDA,MIU

                  counterWS = counterWS + 1

                  counterInt = counterInt + 1
                  Allocate(tempNodeInt)
                  tempNodeInt%value  = ABS(LAMBDA)
                  tempNodeInt%numero = counterInt
                  tempNodeInt%name   = KEYWOR
                  tempNodeInt%next   => listeInt
                  listeInt => tempNodeInt

                  counterInt = counterInt + 1
                  Allocate(tempNodeInt)
                  tempNodeInt%value  = MIU
                  tempNodeInt%numero = counterInt
                  tempNodeInt%name   = KEYWOR
                  tempNodeInt%next   => listeInt
                  listeInt => tempNodeInt

               End Do

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = counterWS
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Additional deformations that will be used to optimize
            ! the WS deformation parameters (extra-push)
            If(trim(keywor).Eq.'PUSH_DEFOR'.OR.trim(keywor).Eq.'PUSH-DEFOR') Then

               L_PUSH = -1
               N_PUSH = 0

               Do While ( L_PUSH < 0 )

                  Read(NFIREA,*) L_PUSH, M_PUSH, Q_PUSH

                  N_PUSH = N_PUSH + 1

                  counterInt = counterInt + 1
                  Allocate(tempNodeInt)
                  tempNodeInt%value  = L_PUSH
                  tempNodeInt%numero = counterInt
                  tempNodeInt%name   = keywor
                  tempNodeInt%next   => listeInt
                  listeInt => tempNodeInt

                  counterInt = counterInt + 1
                  Allocate(tempNodeInt)
                  tempNodeInt%value  = M_PUSH
                  tempNodeInt%numero = counterInt
                  tempNodeInt%name   = keywor
                  tempNodeInt%next   => listeInt
                  listeInt => tempNodeInt

                  counterDbl = counterDbl + 1
                  Allocate(tempNodeDbl)
                  tempNodeDbl%value  = Q_PUSH
                  tempNodeDbl%numero = counterDbl
                  tempNodeDbl%name   = keywor
                  tempNodeDbl%next   => listeDbl
                  listeDbl => tempNodeDbl

               End Do

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = N_PUSH
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = keywor
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Number of iterations to keep the extra constraints on
            If(trim(KEYWOR).Eq.'PUSH_ITERS'.OR.trim(KEYWOR).Eq.'PUSH-ITERS') Then

               Read(NFIREA,*) DOPUSH,NPLAST

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = DOPUSH
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = NPLAST
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Step of the minimization over WS deformations
            If(trim(KEYWOR).Eq.'STEP_MINIM'.OR.trim(KEYWOR).Eq.'STEP-MINIM') Then

               Read(NFIREA,*) STP_WS,ITE_WS

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = STP_WS
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = ITE_WS
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Choosing between a WS or HF-based mean-field
            If(trim(KEYWOR).Eq.'MEAN-FIELD'.OR.trim(KEYWOR).Eq.'MEAN_FIELD') Then

               Read(NFIREA,*) IHFPOT

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IHFPOT
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            !=====================================!
            !    ANGULAR MOMENTUM CONTROL FLAGS   !
            !=====================================!

             ! Rotational frequency
            If(trim(KEYWOR).Eq.'OMEGAY') Then

               Read(NFIREA,*) OMEGAY

               counterDbl = counterDbl + 1
               Allocate(tempNodeDbl)
               tempNodeDbl%value  = OMEGAY
               tempNodeDbl%numero = counterDbl
               tempNodeDbl%name   = KEYWOR
               tempNodeDbl%next   => listeDbl
               listeDbl => tempNodeDbl

            End If

            !================================!
            !   OUTPUT FILE CONTROL FLAGS    !
            !================================!

            ! Type of print-out
            If(trim(KEYWOR).Eq.'ONE_LINE'.OR.trim(KEYWOR).Eq.'ONE-LINE') Then

               Read(NFIREA,*) I1LINE

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = I1LINE
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Type of Nilsson label
            If(trim(KEYWOR).Eq.'NILSSONLAB') Then

               Read(NFIREA,*) NILXYZ

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = NILXYZ
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Flags controlling the deformations extracted from multipoles
            If(trim(KEYWOR).Eq.'BOHR_BETAS'.OR.trim(KEYWOR).Eq.'BOHR-BETAS') Then

               Read(NFIREA,*) NEXBET,IPRIBE,IPRIBL

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = NEXBET
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IPRIBE
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IPRIBL
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            !================================!
            !     I/O FLAGS CONTROL FLAGS    !
            !================================!

            ! Optionally recording ascii file
            If(trim(KEYWOR).Eq.'REVIEW') Then

               Read(NFIREA,*) IREVIE

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IREVIE
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Restart on densities in Gauss-Hermite mesh
            If(trim(KEYWOR).Eq.'FIELD_SAVE'.OR.trim(KEYWOR).Eq.'FIELD-SAVE') Then

               Read(NFIREA,*) IWRIFI

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IWRIFI
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            !==================================!
            !    CHECKPOINTING CONTROL FLAGS   !
            !==================================!

            ! Restart on densities in Gauss-Hermite mesh
            If(trim(KEYWOR).Eq.'RESTART') Then

               Read(NFIREA,*) ICONTI

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = ICONTI
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Restart on pairing
            If(trim(KEYWOR).Eq.'CONT_PAIRI'.OR.trim(KEYWOR).Eq.'CONT-PAIRI') Then

               Read(NFIREA,*) IPCONT

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IPCONT
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Restart on Lipkin-Nogami
            If(trim(KEYWOR).Eq.'CONTLIPKIN') Then

               Read(NFIREA,*) ILCONT

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = ILCONT
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Restart on HFB mean-fields
            If(trim(KEYWOR).Eq.'CONTFIELDS') Then

               Read(NFIREA,*) IFCONT

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IFCONT
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Restart on linear constraint
            If(trim(KEYWOR).Eq.'CONTAUGMEN') Then

               Read(NFIREA,*) IACONT

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IACONT
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Restart on Gogny mean-field
            If(trim(KEYWOR).Eq.'CONTGOGNY') Then

               Read(NFIREA,*) IGCONT

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IGCONT
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Restart on Gogny pairing-field
            If(trim(KEYWOR).Eq.'CONTGOGPAI') Then

               Read(NFIREA,*) IGPCON

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IGPCON
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            ! Restart on WS binary file
            If(trim(KEYWOR).Eq.'READ_WOODS'.OR.trim(KEYWOR).Eq.'READ-WOODS') Then

               Read(NFIREA,*) IREAWS

               counterInt = counterInt + 1
               Allocate(tempNodeInt)
               tempNodeInt%value  = IREAWS
               tempNodeInt%numero = counterInt
               tempNodeInt%name   = KEYWOR
               tempNodeInt%next   => listeInt
               listeInt => tempNodeInt

            End If

            !===================================!
            !  BACKWARD COMPATIBILTY WITH HFODD !
            !===================================!

            If(trim(KEYWOR).Eq.'EXECUTE') Then

               NO_RUN=NO_RUN+1

            End If

         End Do ! end of main DO WHILE loop

         Close(NFIREA)

         counterVector(1) = counterInt
         counterVector(2) = counterDbl
         counterVector(3) = counterStr

      End If ! myrank.eq.0

      ! Broadcasting all the sizes (of integer, real and string data)
      Call MPI_Bcast(counterVector, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

      ! Reassigning local variables from broadcast values
      counterInt = counterVector(1)
      counterDbl = counterVector(2)
      counterStr = counterVector(3)

      Allocate(vectorInt_value(1:counterInt))
      Allocate(vectorDbl_value(1:counterDbl))
      Allocate(vectorStr_value(1:counterStr))

      counterTot = counterInt + counterDbl + counterStr
      Allocate(vectorAll_names(1:counterTot))

      ! Only processor 0 can fill out the vectors
      If(mpi_rank.Eq.0) Then

         ! Now that the final size of these vectors is known, we can
         ! allocate an array and broadcast to everyone else


         ! First, vector of integer data, and the relevant keywords

         Allocate(currentNodeInt)

         currentNodeInt => listeInt

         ii = 0

         Do While ( associated(currentNodeInt) )

            ii = ii + 1

            vectorInt_value(ii) = currentNodeInt%value
            vectorAll_names(ii) = currentNodeInt%name

            Allocate(previousNodeInt)

            previousNodeInt => currentNodeInt
            currentNodeInt  => currentNodeInt%next

            Deallocate(previousNodeInt)

         End Do

         ! Second, vector of real data, and the relevant keywords

         Allocate(currentNodeDbl)

         currentNodeDbl => listeDbl

         ii = 0

         Do While ( associated(currentNodeDbl) )

            ii = ii + 1
            jj = ii + counterInt

            vectorDbl_value(ii) = currentNodeDbl%value
            vectorAll_names(jj) = currentNodeDbl%name

            Allocate(previousNodeDbl)

            previousNodeDbl => currentNodeDbl
            currentNodeDbl  => currentNodeDbl%next

            Deallocate(previousNodeDbl)

         End Do

         ! Third: vector of string data, and the relevant keywords

         Allocate(currentNodeStr)

         currentNodeStr => listeStr

         ii = 0

         Do While ( associated(currentNodeStr) )

            ii = ii + 1
            jj = ii + counterInt + counterDbl

            vectorStr_value(ii) = currentNodeStr%value
            vectorAll_names(jj) = currentNodeStr%name

            Allocate(previousNodeStr)

            previousNodeStr => currentNodeStr
            currentNodeStr  => currentNodeStr%next

            Deallocate(previousNodeStr)

         End Do

      End If ! end of proc = 0 if

      ! Broadcasting all data of all types
      Call MPI_Bcast(vectorInt_value, counterInt, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
      Call MPI_Bcast(vectorDbl_value, counterDbl, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_err)
      Call MPI_Bcast(vectorStr_value, counterStr, mpi_string_10, 0, MPI_COMM_WORLD, mpi_err)

      ! Broadcasting all names (keywords needed to reconstruct data locally)
      Call MPI_Bcast(vectorAll_names, counterTot,  mpi_string_10, 0, MPI_COMM_WORLD, mpi_err)

      ! Derived MPI type is not needed anymore, we release it
      Call mpi_type_free(mpi_string_10, mpi_err)

   End Subroutine mpi_getSequentialData

   !-------------------------------------------------------------------!
   ! This subroutine fills out the common blocks and variables that    !
   ! HFODD needs based on the data read in mpi_nameli() and broadcast  !
   ! by it. The subroutine is executed by all processors               !
   !-------------------------------------------------------------------!

   Subroutine mpi_setSequentialData(IOPTGS,NXHERM,NYHERM,NZHERM, &
                                    IPAIRI,IPAHFB,IPABCS,IMFHFB, &
                                                  LIPKIN,LIPKIP, &
                                    FE2FIN,IF2FIN,FE2FIP,IF2FIP, &
                                                  IROTAT,IREVIE, &
                                    ICONTI,IPCONT,IGCONT,IGPCON, &
                                           ILCONT,IFCONT,IACONT, &
                                    NUCHAO,NUPING,EPSPNG,EPSITE, &
                                           ISIMPY,ISIGNY,IPARTY, &
                                           ISIMTX,ISIMTY,ISIMTZ, &
                                                  NOITER,NULAST, &
                                           NOSCIL,NLIMIT,ENECUT, &
                                           ICOTYP,ICOUDI,ICOUEX, &
                                                  I_GOGA,IGOGPA, &
                                                         I1LINE, &
                                    SLOWEV,SLOWOD,SLOWPA,SLOWLI, &
                                                         ECUTOF, &
                                                         NILXYZ, &
                                    IBROYD,N_ITER,BROTRI,ALPHAM, &
                                           PRHO_T,PRHOST,POWERT, &
                                                         OMEGAY, &
                             IFSHEL,MPOLYN,G_STUN,G_STUP,H_OMFA, &
                                    LAMBAS,MIUBAS,BETBAS,NUMDEF, &
                             LAM_WS,MIU_WS,NUM_WS,ITE_WS,STP_WS, &
                             LAMPUS,MIUPUS,NUMPUS,QMOPUS,IFPUSH, &
                                                  NITPUS,IREAWS, &
                             ISTAND,KETA_J,KETA_W,KETACM,KETA_M, &
                                    IF_RPA,IF_EDF,IF_THO,CBETHO, &
                                           IBASIS,IHFPOT,IFRAGM, &
                                           NMUCON,NMUCOU,NMUPRI, &
                                           NEXBET,IPRIBE,IPRIBL, &
                                                         FCHOM0, &
                                           L_COLL,M_COLL,N_COLL, &
                                                         IWRIFI, &
                                           INNUMB,IZNUMB,R0PARM, &
                             DELTAE,XLOCMX,V2_MIN,ITRMAX,NTHETA, &
                                    MIN_QP,IDEALL,IDELOC,IDECON, &
                                                         TEMP_T)

      Integer, intent(inout) :: ICOTYP, ICOUDI, ICOUEX, I_GOGA, IGOGPA, NOITER, IBROYD, &
                                N_ITER, NULAST, NUPING, NUCHAO, ISIMPY, ISIGNY, IPARTY, &
                                IROTAT, ISIMTX, ISIMTY, ISIMTZ, IPAIRI, IPAHFB, IPABCS, &
                                IMFHFB, LIPKIN, LIPKIP, IOPTGS, NXHERM, NYHERM, NZHERM, &
                                NOSCIL, NLIMIT, I1LINE, NILXYZ, IREVIE, ICONTI, IPCONT, &
                                IGCONT, IGPCON, ILCONT, IFCONT, IACONT, IFSHEL, MPOLYN, &
                                NUMDEF, NUM_WS, ITE_WS, IREAWS, ISTAND, KETA_J, KETA_W, &
                                KETACM, KETA_M, IF_RPA, IF_EDF, IF_THO, IBASIS, IHFPOT, &
                                IFRAGM, NUMPUS, IFPUSH, NITPUS, NMUCON, NMUCOU, NMUPRI, &
                                NEXBET, IPRIBE, IPRIBL, N_COLL, IWRIFI, INNUMB, IZNUMB, &
                                ITRMAX, NTHETA, MIN_QP, IDEALL, IDELOC, IDECON, IF2FIN, &
                                                                                IF2FIP

      Real, intent(inout) :: BROTRI, SLOWEV, SLOWOD, SLOWPA, SLOWLI, EPSITE, &
                             EPSPNG, ECUTOF, ENECUT, PRHO_T, PRHOST, POWERT, &
                             OMEGAY, G_STUN, G_STUP, H_OMFA, STP_WS, ALPHAM, &
                             CBETHO, FCHOM0, R0PARM, DELTAE, XLOCMX, V2_MIN, &
                                                     FE2FIN, FE2FIP, TEMP_T

      Integer :: ii, jj, N_DUMM, counterDef, counterWS, counterColl, counter_push

      Integer, Dimension(:), pointer :: LAMBAS, MIUBAS, LAM_WS, MIU_WS, LAMPUS, MIUPUS
      !Integer, Dimension(:), pointer :: L_COLL, M_COLL
      Integer, Allocatable :: L_COLL(:), M_COLL(:)

      Real :: B_DUMM, PRHO_N, PRHOSN, PRHODN, POWERN, PRHO_P, PRHOSP, PRHODP, POWERP

      Real, Dimension(:), pointer :: BETBAS, QMOPUS

      Character(Len=4)  :: SKYRME
      Character(Len=10) :: KEYWOR, string_10

      ! common blocks for HFODD

      COMMON  /CCPPAI/ PRHO_N,PRHODN,PRHOSN,POWERN, &
                       PRHO_P,PRHODP,PRHOSP,POWERP

      ! Initialize local variables which are not initialized in PREDEF()

      MPOLYN=6; G_STUN=1.2; G_STUP=1.2; H_OMFA=3.5

      NULLIFY(LAM_WS)
      NULLIFY(MIU_WS)

      !NULLIFY(L_COLL)
      !NULLIFY(M_COLL)

      NULLIFY(LAMPUS)
      NULLIFY(MIUPUS)
      NULLIFY(QMOPUS)

      NUM_WS=1
      STP_WS=0.01

      NULLIFY(LAMBAS)
      NULLIFY(MIUBAS)

      NULLIFY(BETBAS)

      NUMDEF=1
      NUMPUS=0

      !=================================================!
      !                                                 !
      !   RECONSTRUCTION OF THE INPUT OF TYPE INTEGER   !
      !                                                 !
      !=================================================!

      ii = 1

      Do While ( ii .le. counterInt )

         KEYWOR = vectorAll_names(ii)

         !==============================================!
         !   INTEGER: ITERATIONS CONTROL FLAGS          !
         !==============================================!

         ! Number of iterations
         If(trim(KEYWOR).Eq.'ITERATIONS') Then
             NOITER = vectorInt_value(ii)
         End If

         ! Switches on/off the Broyden method
         If(trim(KEYWOR).Eq.'BROYDEN') Then

               N_DUMM = vectorInt_value(ii)
            If(N_DUMM.GT.0) N_ITER=N_DUMM

            ii = ii +1
            IBROYD = vectorInt_value(ii)

         End If

         ! Defines divergence
         If(trim(KEYWOR).Eq.'MAXANTIOSC') Then
            NULAST = vectorInt_value(ii)
         End If

         ! Control ping-pong effects
         If(trim(KEYWOR).Eq.'PING_PONG'.OR.trim(KEYWOR).Eq.'PING-PONG') Then
            NUPING = vectorInt_value(ii)
         End If

         ! Control chaotic divergence
         If(trim(KEYWOR).Eq.'CHAOTIC') Then
            NUCHAO = vectorInt_value(ii)
         End If

         !==============================================!
         !   INTEGER: SPECIFIC FEATURES CONTROL FLAGS   !
         !==============================================!

         ! Type of shell correction
         If(trim(KEYWOR).Eq.'SHELLCORCT') Then
            IFSHEL = vectorInt_value(ii)
         End If

         ! Number of polynomials for the shell correction
         If(trim(KEYWOR).Eq.'SHELLPARAM') Then
            MPOLYN = vectorInt_value(ii)
         End If

         ! Switching to interface with HFBTHO
         If(trim(KEYWOR).Eq.'HFBTHOISON') Then
            IF_THO = vectorInt_value(ii)
         End If

         ! Compute mass and charge of fission fragments
         If(trim(KEYWOR).Eq.'MASSFRAGME') Then
            IFRAGM = vectorInt_value(ii)
         End If

         ! Compute and record densities of the fission fragments
         If(trim(KEYWOR).Eq.'FRAGDENSIT') Then
            IDECON = vectorInt_value(ii)

            ii = ii + 1
            IDELOC = vectorInt_value(ii)

            ii = ii + 1
            IDEALL = vectorInt_value(ii)
         End If

         ! Rotate selected pairs of q.p. wave-functions
         If(trim(KEYWOR).Eq.'QPROTATION') Then
            NTHETA = vectorInt_value(ii)

            ii = ii + 1
            ITRMAX = vectorInt_value(ii)

            ii = ii + 1
            MIN_QP = vectorInt_value(ii)
         End If

         ! Definition of the collective space used for calculation of the
         ! collective inertia
         If(trim(KEYWOR).Eq.'COLL_SPACE'.OR.trim(KEYWOR).Eq.'COLL-SPACE') Then

            counterColl = vectorInt_value(ii)

            N_COLL = counterColl

            Allocate(L_COLL(1:counterColl))
            Allocate(M_COLL(1:counterColl))

            Do jj = 1, counterColl

               ii = ii + 1
               M_COLL(jj) = vectorInt_value(ii)

               ii = ii + 1
               L_COLL(jj) = vectorInt_value(ii)

            End Do

         End If

         !==============================================!
         !   INTEGER: INTERACTIONS CONTROL FLAGS        !
         !==============================================!

         ! Switching to Energy Density Functional formalism
         If(trim(KEYWOR).Eq.'UNEDF_PROJ'.OR.trim(KEYWOR).Eq.'UNEDF-PROJ') Then
            IF_EDF = vectorInt_value(ii)
         End If

         ! Skyrme force flags
         If(trim(KEYWOR).Eq.'SKYRME_STD'.OR.trim(KEYWOR).Eq.'SKYRME-STD') Then

            KETA_M = vectorInt_value(ii)

            ii = ii + 1
            KETACM = vectorInt_value(ii)

            ii = ii + 1
            KETA_W = vectorInt_value(ii)

            ii = ii + 1
            KETA_J = vectorInt_value(ii)

            ii = ii + 1
            ISTAND = vectorInt_value(ii)

         End If

         ! Type of Coulomb potential
         If(trim(KEYWOR).Eq.'COULOMBPAR') Then

            ICOUEX = vectorInt_value(ii)

            ii = ii + 1
            ICOUDI = vectorInt_value(ii)

            ii = ii + 1
            ICOTYP = vectorInt_value(ii)

         End If

         ! Gogny force
         If(trim(KEYWOR).Eq.'GOGNY') Then
            I_GOGA = vectorInt_value(ii)
         End If

         !==============================================!
         !   INTEGER: PAIRING CONTROL FLAGS             !
         !==============================================!

         ! Gogny pairing
         If(trim(KEYWOR).Eq.'GOGNY_PAIR'.OR.trim(KEYWOR).Eq.'GOGNY-PAIR') Then
            IGOGPA = vectorInt_value(ii)
         End If

         ! Switches on/off pairing
         If(trim(KEYWOR).Eq.'PAIRING') Then
            IPAIRI = vectorInt_value(ii)
         End If

         ! HFB pairing
         If(trim(KEYWOR).Eq.'HFB') Then
            IPAHFB = vectorInt_value(ii)
         End If

         ! BCS pairing
         If(trim(KEYWOR).Eq.'BCS') Then
            IPABCS = vectorInt_value(ii)
         End If

         ! Diagonalizes the mean-field (if HFB)
         If(trim(KEYWOR).Eq.'HFBMEANFLD') Then
            IMFHFB = vectorInt_value(ii)
         End If

         ! Lipkin-Nogami correction
         If(trim(KEYWOR).Eq.'LIPKIN') Then

            LIPKIP = vectorInt_value(ii)

            ii = ii + 1
            LIPKIN = vectorInt_value(ii)

         End If

         ! Reading the fixed Lipkin-Nogami lambda2 parameters
         If(trim(KEYWOR).Eq.'FIXLAMB2_N'.OR.trim(KEYWOR).EQ.'FIXLAMB2-N') Then
            IF2FIN = vectorInt_value(ii)
         End If

         If(trim(KEYWOR).Eq.'FIXLAMB2_P'.OR.trim(KEYWOR).EQ.'FIXLAMB2-P') Then
            IF2FIP = vectorInt_value(ii)
         End If

         !==============================================!
         !   INTEGER: SYMMETRY CONTROL FLAGS            !
         !==============================================!

         ! time-reversal
         If(trim(KEYWOR).Eq.'ROTATION') Then
            IROTAT = vectorInt_value(ii)
         End If

         ! y-simplex
         If(trim(KEYWOR).Eq.'SIMPLEXY') Then
            ISIMPY = vectorInt_value(ii)
         End If

         ! y-signature
         If(trim(KEYWOR).Eq.'SIGNATUREY') Then
            ISIGNY = vectorInt_value(ii)
         End If

         ! parity
         If(trim(KEYWOR).Eq.'PARITY') Then
            IPARTY = vectorInt_value(ii)
         End If

         ! All 3 simplexes at once
         If(trim(KEYWOR).Eq.'TSIMPLEX3D') Then

            ISIMTZ = vectorInt_value(ii)

            ii = ii + 1
            ISIMTY = vectorInt_value(ii)

            ii = ii + 1
            ISIMTX = vectorInt_value(ii)

         End If

         !==============================================!
         !   INTEGER: HO BASIS CONTROL FLAGS            !
         !==============================================!

         ! Automatic adjustment of HO basis parameters
         If(trim(KEYWOR).Eq.'BASISAUTOM') Then
            IBASIS = vectorInt_value(ii)
         End If

         ! Size of the HO basis
         If(trim(KEYWOR).Eq.'BASIS_SIZE'.OR.trim(KEYWOR).Eq.'BASIS-SIZE') Then

            NLIMIT = vectorInt_value(ii)

            ii = ii + 1
            NOSCIL = vectorInt_value(ii)

         End If

         ! Particle numbers of the HO basis
         If(trim(KEYWOR).Eq.'SURFAC_PAR'.OR.trim(KEYWOR).Eq.'SURFAC-PAR') Then

            IZNUMB = vectorInt_value(ii)

            ii = ii + 1
            INNUMB = vectorInt_value(ii)

         End If

         ! Deformation of the HO basis
         If(trim(KEYWOR).Eq.'SURFAC_DEF'.OR.trim(KEYWOR).Eq.'SURFAC-DEF') Then

            counterDef = vectorInt_value(ii)

            NUMDEF = counterDef

            Allocate(LAMBAS(1:counterDef))
            Allocate(MIUBAS(1:counterDef))

            Do jj = 1, counterDef

               ii = ii + 1
               MIUBAS(jj) = vectorInt_value(ii)

               ii = ii + 1
               LAMBAS(jj) = vectorInt_value(ii)

             End Do

         End If

         ! Automatic number of Gauss-Hermite points
         If(trim(KEYWOR).Eq.'OPTI_GAUSS'.OR.trim(KEYWOR).Eq.'OPTI-GAUSS') Then
            IOPTGS = vectorInt_value(ii)
         End If

         ! Number of Gauss-Hermite points
         If(trim(KEYWOR).Eq.'GAUSHERMIT') Then

            NZHERM = vectorInt_value(ii)

            ii = ii + 1
            NYHERM = vectorInt_value(ii)

            ii = ii + 1
            NXHERM = vectorInt_value(ii)

         End If

         !==============================================!
         !   INTEGER: MULTIPOLE MOMENT CONTROL FLAGS    !
         !==============================================!

         If(trim(KEYWOR).Eq.'RPA_CONSTR'.OR.trim(KEYWOR).Eq.'RPA-CONSTR') Then
            IF_RPA = vectorInt_value(ii)
         End If

         ! Number of iterations to keep the extra constraints on
         If(trim(KEYWOR).Eq.'MAX_MULTIP'.OR.trim(KEYWOR).Eq.'MAX-MULTIP') Then

            NMUPRI = vectorInt_value(ii)

            ii = ii + 1
            NMUCOU = vectorInt_value(ii)

            ii = ii + 1
            NMUCON = vectorInt_value(ii)

         End If

         ! Deformations of the WS potential used to provide initial solution
         If(trim(KEYWOR).Eq.'DEFORMS_WS'.OR.trim(KEYWOR).Eq.'DEFORMS-WS') Then

            counterWS = vectorInt_value(ii)

            NUM_WS = counterWS

            Allocate(LAM_WS(1:counterWS))
            Allocate(MIU_WS(1:counterWS))

            Do jj = 1, counterWS

               ii = ii + 1
               MIU_WS(jj) = vectorInt_value(ii)

               ii = ii + 1
               LAM_WS(jj) = vectorInt_value(ii)

            End Do

         End If

         ! Additional deformations that will be used to optimize
         ! the WS deformation parameters (extra-push)
         If(trim(KEYWOR).Eq.'PUSH_DEFOR'.OR.trim(KEYWOR).Eq.'PUSH-DEFOR') Then

            counter_push = vectorInt_value(ii)

            NUMPUS = counter_push

            Allocate(LAMPUS(1:counter_push))
            Allocate(MIUPUS(1:counter_push))

            Do jj = 1, counter_push

               ii = ii + 1
               MIUPUS(jj) = vectorInt_value(ii)

               ii = ii + 1
               LAMPUS(jj) = vectorInt_value(ii)

            End Do

         End If

         ! Number of iterations to keep the extra constraints on
         If(trim(KEYWOR).Eq.'PUSH_ITERS'.OR.trim(KEYWOR).Eq.'PUSH-ITERS') Then

            NITPUS = vectorInt_value(ii)

            ii = ii + 1
            IFPUSH = vectorInt_value(ii)

         End If

         ! Number of iterations of the minimization
         If(trim(KEYWOR).Eq.'STEP_MINIM'.OR.trim(KEYWOR).Eq.'STEP-MINIM') Then
            ITE_WS = vectorInt_value(ii)
         End If

         ! Choosing between a WS or HF-based mean-field
         If(trim(KEYWOR).Eq.'MEAN-FIELD'.OR.trim(KEYWOR).Eq.'MEAN_FIELD') Then
            IHFPOT = vectorInt_value(ii)
         End If

         !==============================================!
         !   INTEGER: OUTPUT FILE CONTROL FLAGS         !
         !==============================================!

         ! Type of print-out
         If(trim(KEYWOR).Eq.'ONE_LINE'.OR.trim(KEYWOR).Eq.'ONE-LINE') Then
            I1LINE = vectorInt_value(ii)
         End If

         ! Type of Nilsson label
         If(trim(KEYWOR).Eq.'NILSSONLAB') Then
            NILXYZ = vectorInt_value(ii)
         End If

         ! Number of iterations to keep the extra constraints on
         If(trim(KEYWOR).Eq.'BOHR_BETAS'.OR.trim(KEYWOR).Eq.'BOHR-BETAS') Then

            IPRIBL = vectorInt_value(ii)

            ii = ii + 1
            IPRIBE = vectorInt_value(ii)

            ii = ii + 1
            NEXBET = vectorInt_value(ii)

         End If

         !==============================================!
         !   INTEGER: I/O FLAGS CONTROL FLAGS           !
         !==============================================!

         ! Optionally recording ascii file
         If(trim(KEYWOR).Eq.'REVIEW') Then
            IREVIE = vectorInt_value(ii)
         End If

         ! Restart on densities in Gauss-Hermite mesh
         If(trim(KEYWOR).Eq.'FIELD_SAVE'.OR.trim(KEYWOR).Eq.'FIELD-SAVE') Then
            IWRIFI = vectorInt_value(ii)
         End If

         !==============================================!
         !   INTEGER: CHECKPOINTING CONTROL FLAGS       !
         !==============================================!

         ! Restart on densities
         If(trim(KEYWOR).Eq.'RESTART') Then
             ICONTI = vectorInt_value(ii)
         End If

         ! Restart on pairing
         If(trim(KEYWOR).Eq.'CONT_PAIRI'.OR.trim(KEYWOR).Eq.'CONT-PAIRI') Then
            IPCONT = vectorInt_value(ii)
         End If

         ! Restart on Lipkin-Nogami
         If(trim(KEYWOR).Eq.'CONTLIPKIN') Then
            ILCONT = vectorInt_value(ii)
         End If

         ! Restart on HFB mean-fields
         If(trim(KEYWOR).Eq.'CONTFIELDS') Then
            IFCONT = vectorInt_value(ii)
         End If

         ! Restart on linear constraint
         If(trim(KEYWOR).Eq.'CONTAUGMEN') Then
            IACONT = vectorInt_value(ii)
         End If

         ! Restart on Gogny mean-field
         If(trim(KEYWOR).Eq.'CONTGOGNY') Then
            IGCONT = vectorInt_value(ii)
         End If

         ! Restart on Gogny pairing-field
         If(trim(KEYWOR).Eq.'CONTGOGPAI') Then
            IGPCON = vectorInt_value(ii)
         End If

         ! Restart on WS binary files
         If(trim(KEYWOR).Eq.'READ_WOODS'.OR.trim(KEYWOR).Eq.'READ-WOODS') Then
            IREAWS = vectorInt_value(ii)
         End If

         ! Go to the next item on the list
         ii = ii + 1

      End Do

      Deallocate(vectorInt_value)

      !==============================================!
      !                                              !
      !   RECONSTRUCTION OF THE INPUT OF TYPE REAL   !
      !                                              !
      !==============================================!

      ii = 1

      Do While ( ii .le. counterDbl )

         KEYWOR = vectorAll_names(ii + counterInt)

         !===========================================!
         !   REAL: ITERATIONS CONTROL FLAGS          !
         !===========================================!

         ! Switches on/off the Broyden method
         If(trim(KEYWOR).Eq.'BROYDEN') Then

               B_DUMM = vectorDbl_value(ii)
            If(B_DUMM.GT.0.0D0) BROTRI=B_DUMM

            ii = ii + 1
            ALPHAM = vectorDbl_value(ii)

         End If

         ! Slow-down parameters for the mean-field
         If(trim(KEYWOR).Eq.'SLOW_DOWN'.OR.trim(KEYWOR).Eq.'SLOW-DOWN') Then

            SLOWOD = vectorDbl_value(ii)

            ii = ii + 1
            SLOWEV = vectorDbl_value(ii)

         End If

         ! Slow-down parameters for the pairing field
         If(trim(KEYWOR).Eq.'SLOW_PAIR'.OR.trim(KEYWOR).Eq.'SLOW-PAIR') Then
            SLOWPA = vectorDbl_value(ii)
         End If

         ! Slow-down parameters for the Lipkin-Nogami correction
         If(trim(KEYWOR).Eq.'SLOWLIPKIN') Then
            SLOWLI = vectorDbl_value(ii)
         End If

         ! Defines the precision of convergence
         If(trim(KEYWOR).Eq.'ITERAT_EPS'.OR.trim(KEYWOR).Eq.'ITERAT-EPS') Then
            EPSITE = vectorDbl_value(ii)
         End If

         ! Control ping-pong effects
         If(trim(KEYWOR).Eq.'PING_PONG'.OR.trim(KEYWOR).Eq.'PING-PONG') Then
            EPSPNG = vectorDbl_value(ii)
         End If

         !===========================================!
         !   REAL: SPECIFIC FEATURES CONTROL FLAGS   !
         !===========================================!

         ! Temperature
         If(trim(KEYWOR).Eq.'FINITETEMP') Then
            TEMP_T = vectorDbl_value(ii)
         End If

         ! Parameters for the shell correction
         If(trim(KEYWOR).Eq.'SHELLPARAM') Then

            H_OMFA = vectorDbl_value(ii)

            ii = ii + 1
            G_STUP = vectorDbl_value(ii)

            ii = ii + 1
            G_STUN = vectorDbl_value(ii)

         End If

         ! Switching to interface with HFBTHO
         If(trim(KEYWOR).Eq.'HFBTHOISON') Then
            CBETHO = vectorDbl_value(ii)
         End If

         ! Rotate selected pairs of q.p. wave-functions
         If(trim(KEYWOR).Eq.'QPROTATION') Then

            V2_MIN = vectorDbl_value(ii)

            ii = ii + 1
            XLOCMX = vectorDbl_value(ii)

            ii = ii + 1
            DELTAE = vectorDbl_value(ii)

         End If

         !===========================================!
         !   REAL: PAIRING CONTROL FLAGS             !
         !===========================================!

         ! Control chaotic divergence
         If(trim(KEYWOR).Eq.'CUTOFF') Then
            ECUTOF = vectorDbl_value(ii)
         End If

         ! Pairing Interaction
         If(trim(KEYWOR).Eq.'PAIR_INTER'.OR.trim(KEYWOR).Eq.'PAIR-INTER') Then

            POWERT = vectorDbl_value(ii)

            ii = ii + 1
            PRHOST = vectorDbl_value(ii)

            ii = ii + 1
            PRHO_T = vectorDbl_value(ii)

         End If

         ! The fixed Lipkin-Nogami lambda2 parameters
         If(trim(KEYWOR).Eq.'FIXLAMB2_N'.OR.trim(KEYWOR).Eq.'FIXLAMB2-N') Then
            FE2FIN = vectorDbl_value(ii)
         End If

         If(trim(KEYWOR).Eq.'FIXLAMB2_P'.OR.trim(KEYWOR).Eq.'FIXLAMB2-P') Then
            FE2FIP = vectorDbl_value(ii)
         End If

         ! Pairing Interaction (neutrons only)
         If(trim(KEYWOR).Eq.'PAIRNINTER') Then

            POWERN = vectorDbl_value(ii)

            ii = ii + 1
            PRHOSN = vectorDbl_value(ii)

            ii = ii + 1
            PRHO_N = vectorDbl_value(ii)

         End If

         ! Pairing Interaction (protons only)
         If(trim(KEYWOR).Eq.'PAIRPINTER') Then

            POWERP = vectorDbl_value(ii)

            ii = ii + 1
            PRHOSP = vectorDbl_value(ii)

            ii = ii + 1
            PRHO_P = vectorDbl_value(ii)

         End If

         !===========================================!
         !   REAL: HO BASIS CONTROL FLAGS            !
         !===========================================!

         ! Flag controling the basis frequency
         If(trim(KEYWOR).Eq.'HOMEGAZERO') Then
            FCHOM0 = vectorDbl_value(ii)
         End If

         ! Size of the HO basis
         If(trim(KEYWOR).Eq.'BASIS_SIZE'.OR.trim(KEYWOR).Eq.'BASIS-SIZE') Then
            ENECUT = vectorDbl_value(ii)
         End If

         ! Particle numbers of the HO basis
         If(trim(KEYWOR).Eq.'SURFAC_PAR'.OR.trim(KEYWOR).Eq.'SURFAC-PAR') Then
            R0PARM = vectorDbl_value(ii)
         End If

         ! Deformation of the HO basis
         If(trim(KEYWOR).Eq.'SURFAC_DEF'.OR.trim(KEYWOR).Eq.'SURFAC-DEF') Then

            Allocate(BETBAS(1:counterDef))

            Do jj = 1, counterDef

               BETBAS(jj) = vectorDbl_value(ii)

               ii = ii + 1

             End Do

         End If

         !===========================================!
         !   REAL: MULTIPOLE MOMENT CONTROL FLAGS    !
         !===========================================!

         ! Deformations of the WS potential
         If(trim(KEYWOR).Eq.'PUSH_DEFOR'.OR.trim(KEYWOR).Eq.'PUSH-DEFOR') Then

            Allocate(QMOPUS(1:counter_push))

            DO jj = 1, counter_push

               QMOPUS(jj) = vectorDbl_value(ii)
               ii = ii + 1

            End Do
            ii = ii - 1

         End If

         ! Step of the minimization over WS deformations
         If(trim(KEYWOR).Eq.'STEP_MINIM'.OR.trim(KEYWOR).Eq.'STEP-MINIM') Then
            STP_WS = vectorDbl_value(ii)
         End If

         !===========================================!
         !   REAL: ANGULAR MOMENTUM CONTROL FLAGS    !
         !===========================================!

         ! Rotational frequency
         If(trim(KEYWOR).Eq.'OMEGAY') Then
            OMEGAY = vectorDbl_value(ii)
         End If

         ! Go to the next item on the list
         ii = ii + 1

      End Do

      Deallocate(vectorDbl_value)

      !================================================!
      !                                                !
      ! RECONSTRUCTION OF THE INPUT OF TYPE STRING*10  !
      !                                                !
      !================================================!

      ii = 1

      Do While ( ii .le. counterStr )

         KEYWOR = vectorAll_names(ii + counterInt + counterDbl)

         If(trim(KEYWOR).Eq.'SKYRME-SET') Then
            string_10 = vectorStr_value(ii)
            SKYRME = string_10(1:4)
         End If

         ! Go to the next item on the list
         ii = ii + 1

      End Do

      Deallocate(vectorStr_value)

      Deallocate(vectorAll_names)

   End Subroutine mpi_setSequentialData

   !---------------------------------------------------------------!
   ! This SUBROUTINE reads data from input file. These data will   !
   ! then be used to construct a collection of jobs for HFODD.     !
   !---------------------------------------------------------------!

   Subroutine mpi_getParallelData()

      Integer, Parameter :: ndatin = 83

      Character(Len=4) :: SKYRME
      Character(Len=80) :: filein
      Character(Len=10) :: keywor

      Integer :: ii,jj,ierror,numero,mpi_rank,mpi_err,mpi_string_10,mpidef,mpibas
      Integer :: lambda,miu,numberQ,iq,lambda_res,miu_res,numberQ_res,iq_res,iq_batch
      Integer :: izstrt,izstep,nstepz,instrt,instep,nstepn,n_mini,n_step,noffre,nofshl, &
                 nofb20,nsmini,nsstep,nofsta,ifcons,ibatch,nbatch,if_opt,read_next

      Real :: xzstrt,xzstep,xnstrt,xnstep,o_mini,o_step,b20min,b20stp, &
              qBegin,qFin,qBegin_res,qFin_res

      Type(Node_int), pointer :: tempNodeInt, currentNodeInt, previousNodeInt
      Type(Node_dbl), pointer :: tempNodeDbl, currentNodeDbl, previousNodeDbl
      Type(Node_str), pointer :: tempNodeStr, currentNodeStr, previousNodeStr

      ! Getting the rank of the current process
      Call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, mpi_err)

      ! Define MPI type for a character string of length stringLength and commit
      Call mpi_type_contiguous(stringLength, MPI_CHARACTER, mpi_string_10, mpi_err)
      Call mpi_type_commit(mpi_string_10, mpi_err)

      ! Only the root process reads the data
      If(mpi_rank.Eq.0) Then

         filein = "hfodd_mpiio.d"
         Open(ndatin,file=filein,status='old',form='formatted',iostat=ierror)

         If(ierror .NE. 0) Then
            Write(6,'("error in opening file : ",a80)') filein
            write(6,'("error code            : ",i12)') ierror
            Stop 'error in readin - i/o'
         End If

         keywor="Initialize"

         counterInt = 0
         counterDbl = 0
         counterStr = 0

         counterInt = counterInt + 1
         Allocate(listeInt)
         listeInt%value  = 0
         listeInt%numero = counterInt
         listeInt%name   = keywor
         NULLIFY(listeInt%next)

         counterDbl = counterDbl + 1
         Allocate(listeDbl)
         listeDbl%value  = 0.0
         listeDbl%numero = counterDbl
         listeDbl%name   = keywor
         NULLIFY(listeDbl%next)

         counterStr = counterStr + 1
         Allocate(listeStr)
         listeStr%value  = "INITIALIZE"
         listeStr%numero = counterStr
         listeStr%name   = KEYWOR
         NULLIFY(listeStr%next)

         Do While (trim(keywor) .NE. 'END_DATA')

             Read(ndatin, FMT="(A)") keywor

             ! Type of calculation (deformation, deformation with restart,
             ! basis)
             If(trim(keywor) == 'CALCULMODE') Then

                Read(ndatin,*) mpidef, mpibas

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = mpidef
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = mpibas
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

             End If

             ! Setting up constrained or unconstrained calculations
             If(trim(keywor) == 'CONSTR_DEF'.OR.trim(keywor) == 'CONSTR-DEF') Then

                Read(ndatin,*) ifcons

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = ifcons
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

             End If

             ! Setting up batch mode for constrained calculations
             If(trim(keywor) == 'BATCH_MODE'.OR.trim(keywor) == 'BATCH-MODE') Then

                Read(ndatin,*) ibatch, nbatch

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = ibatch
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = nbatch
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

             End If

             ! Specifying which multipoles will be 'batched'
             If(trim(keywor) == 'BATCH_SPEC'.OR.trim(keywor) == 'BATCH-SPEC') Then

                read_next = 1
                iq_batch = 0

                Do While ( read_next == 1 )

                   Read(ndatin,*) lambda, miu, read_next

                   iq_batch = iq_batch + 1

                   counterInt = counterInt + 1
                   Allocate(tempNodeInt)
                   tempNodeInt%value  = lambda
                   tempNodeInt%numero = counterInt
                   tempNodeInt%name   = keywor
                   tempNodeInt%next   => listeInt
                   listeInt => tempNodeInt

                   counterInt = counterInt + 1
                   Allocate(tempNodeInt)
                   tempNodeInt%value  = miu
                   tempNodeInt%numero = counterInt
                   tempNodeInt%name   = keywor
                   tempNodeInt%next   => listeInt
                   listeInt => tempNodeInt

                End Do

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = iq_batch
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

             End If

             ! Defining nuclear grid
             If(trim(keywor) == 'ALL_NUCLEI'.OR.trim(keywor) == 'ALL-NUCLEI') Then

                Read(ndatin,*) izstrt, izstep, nstepz, instrt, instep, nstepn

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = izstrt
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = izstep
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = nstepz
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = instrt
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = instep
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = nstepn
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

             End If

             ! Defining nuclear grid
             If(trim(keywor) == 'REALNUCLEI') Then

                Read(ndatin,*) xzstrt, xzstep, nstepz, xnstrt, xnstep, nstepn

                counterDbl = counterDbl + 1
                Allocate(tempNodeDbl)
                tempNodeDbl%value  = xzstrt
                tempNodeDbl%numero = counterDbl
                tempNodeDbl%name   = keywor
                tempNodeDbl%next   => listeDbl
                listeDbl => tempNodeDbl

                counterDbl = counterDbl + 1
                Allocate(tempNodeDbl)
                tempNodeDbl%value  = xzstep
                tempNodeDbl%numero = counterDbl
                tempNodeDbl%name   = keywor
                tempNodeDbl%next   => listeDbl
                listeDbl => tempNodeDbl

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = nstepz
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

                counterDbl = counterDbl + 1
                Allocate(tempNodeDbl)
                tempNodeDbl%value  = xnstrt
                tempNodeDbl%numero = counterDbl
                tempNodeDbl%name   = keywor
                tempNodeDbl%next   => listeDbl
                listeDbl => tempNodeDbl

                counterDbl = counterDbl + 1
                Allocate(tempNodeDbl)
                tempNodeDbl%value  = xnstep
                tempNodeDbl%numero = counterDbl
                tempNodeDbl%name   = keywor
                tempNodeDbl%next   => listeDbl
                listeDbl => tempNodeDbl

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = nstepn
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

             End If

             ! Defining default multipole grid
             If(trim(keywor) == 'MULTICONST') Then

                lambda = -1
                iq = 0

                Do While ( lambda < 0 )

                   iq = iq + 1

                   Read(ndatin,*) lambda, miu, qBegin, qFin, numberQ

                   counterInt = counterInt + 1
                   Allocate(tempNodeInt)
                   tempNodeInt%value  = abs(lambda)
                   tempNodeInt%numero = counterInt
                   tempNodeInt%name   = keywor
                   tempNodeInt%next   => listeInt
                   listeInt => tempNodeInt

                   counterInt = counterInt + 1
                   Allocate(tempNodeInt)
                   tempNodeInt%value  = miu
                   tempNodeInt%numero = counterInt
                   tempNodeInt%name   = keywor
                   tempNodeInt%next   => listeInt
                   listeInt => tempNodeInt

                   counterInt = counterInt + 1
                   Allocate(tempNodeInt)
                   tempNodeInt%value  = numberQ
                   tempNodeInt%numero = counterInt
                   tempNodeInt%name   = keywor
                   tempNodeInt%next   => listeInt
                   listeInt => tempNodeInt

                   counterDbl = counterDbl + 1
                   Allocate(tempNodeDbl)
                   tempNodeDbl%value  = qBegin
                   tempNodeDbl%numero = counterDbl
                   tempNodeDbl%name   = keywor
                   tempNodeDbl%next   => listeDbl
                   listeDbl => tempNodeDbl

                   counterDbl = counterDbl + 1
                   Allocate(tempNodeDbl)
                   tempNodeDbl%value  = qFin
                   tempNodeDbl%numero = counterDbl
                   tempNodeDbl%name   = keywor
                   tempNodeDbl%next   => listeDbl
                   listeDbl => tempNodeDbl

                End Do

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = iq
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

             End If

             ! Defining new multipole grid (active only if mpidef=2). If
             ! active, information under keyword MULTICONST will be used
             ! to identify the restart points
             If(trim(keywor) == 'MULTIRESTA') Then

                lambda_res = -1
                iq_res = 0

                Do While ( lambda_res < 0 )

                   iq_res = iq_res + 1

                   Read(ndatin,*) lambda_res, miu_res, qBegin_res, qFin_res, numberQ_res

                   counterInt = counterInt + 1
                   Allocate(tempNodeInt)
                   tempNodeInt%value  = abs(lambda_res)
                   tempNodeInt%numero = counterInt
                   tempNodeInt%name   = keywor
                   tempNodeInt%next   => listeInt
                   listeInt => tempNodeInt

                   counterInt = counterInt + 1
                   Allocate(tempNodeInt)
                   tempNodeInt%value  = miu_res
                   tempNodeInt%numero = counterInt
                   tempNodeInt%name   = keywor
                   tempNodeInt%next   => listeInt
                   listeInt => tempNodeInt

                   counterInt = counterInt + 1
                   Allocate(tempNodeInt)
                   tempNodeInt%value  = numberQ_res
                   tempNodeInt%numero = counterInt
                   tempNodeInt%name   = keywor
                   tempNodeInt%next   => listeInt
                   listeInt => tempNodeInt

                   counterDbl = counterDbl + 1
                   Allocate(tempNodeDbl)
                   tempNodeDbl%value  = qBegin_res
                   tempNodeDbl%numero = counterDbl
                   tempNodeDbl%name   = keywor
                   tempNodeDbl%next   => listeDbl
                   listeDbl => tempNodeDbl

                   counterDbl = counterDbl + 1
                   Allocate(tempNodeDbl)
                   tempNodeDbl%value  = qFin_res
                   tempNodeDbl%numero = counterDbl
                   tempNodeDbl%name   = keywor
                   tempNodeDbl%next   => listeDbl
                   listeDbl => tempNodeDbl

                End Do

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = iq_res
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

             End If

             ! Choose to optimize the deformation grid by centering new
             ! multipolarity constraints on the restart equilibrium value
             If(trim(keywor) == 'OPTIM_GRID'.OR.trim(keywor) == 'OPTIM-GRID') Then

                Read(ndatin,*) if_opt

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = if_opt
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

             End If

             ! Skyrme force name
             If(trim(KEYWOR).Eq.'ALL_FORCES'.OR.trim(KEYWOR).Eq.'ALL-FORCES') Then

                counterForce = 0
                numero = -1

                Do While (numero < 0)

                   Read(ndatin,*) numero, SKYRME

                   counterForce = counterForce + 1

                   counterStr = counterStr + 1
                   Allocate(tempNodeStr)
                   tempNodeStr%value  = trim(SKYRME)//"      "
                   tempNodeStr%numero = counterStr
                   tempNodeStr%name   = KEYWOR
                   tempNodeStr%next   => listeStr
                   listeStr => tempNodeStr

                End Do

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = counterForce
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = KEYWOR
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

             End If

             ! For mpibas=1, grid for the basis frequency
             If(trim(keywor) == 'BASIS_FREQ'.OR.trim(keywor) == 'BASIS-FREQ') Then

                Read(ndatin,*) o_mini, o_step, noffre

                counterDbl = counterDbl + 1
                Allocate(tempNodeDbl)
                tempNodeDbl%value  = o_mini
                tempNodeDbl%numero = counterDbl
                tempNodeDbl%name   = keywor
                tempNodeDbl%next   => listeDbl
                listeDbl => tempNodeDbl

                counterDbl = counterDbl + 1
                Allocate(tempNodeDbl)
                tempNodeDbl%value  = o_step
                tempNodeDbl%numero = counterDbl
                tempNodeDbl%name   = keywor
                tempNodeDbl%next   => listeDbl
                listeDbl => tempNodeDbl

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = noffre
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

             End If

             ! For mpibas=1, grid for the basis number of shells
             If(trim(keywor) == 'BASIS_NSHL'.OR.trim(keywor) == 'BASIS-NSHL') Then

                Read(ndatin,*) n_mini, n_step, nofshl

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = n_mini
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = n_step
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = nofshl
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

             End If

             ! For mpibas=1, grid for the basis number of states
             If(trim(keywor) == 'BASIS_NSTA'.OR.trim(keywor) == 'BASIS-NSTA') Then

                Read(ndatin,*) nsmini, nsstep, nofsta

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = nsmini
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = nsstep
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = nofsta
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

             End If

             ! For mpibas=1, grid for the basis deformation
             If(trim(keywor) == 'BASIS_DEFS'.OR.trim(keywor) == 'BASIS-DEFS') Then

                Read(ndatin,*) b20min, b20stp, nofb20

                counterDbl = counterDbl + 1
                Allocate(tempNodeDbl)
                tempNodeDbl%value  = b20min
                tempNodeDbl%numero = counterDbl
                tempNodeDbl%name   = keywor
                tempNodeDbl%next   => listeDbl
                listeDbl => tempNodeDbl

                counterDbl = counterDbl + 1
                Allocate(tempNodeDbl)
                tempNodeDbl%value  = b20stp
                tempNodeDbl%numero = counterDbl
                tempNodeDbl%name   = keywor
                tempNodeDbl%next   => listeDbl
                listeDbl => tempNodeDbl

                counterInt = counterInt + 1
                Allocate(tempNodeInt)
                tempNodeInt%value  = nofb20
                tempNodeInt%numero = counterInt
                tempNodeInt%name   = keywor
                tempNodeInt%next   => listeInt
                listeInt => tempNodeInt

             End If

             If(trim(keywor) == 'END_DATA') Then
                Exit
             End If

          End Do

         Close(ndatin)

          counterVector(1) = counterInt
          counterVector(2) = counterDbl
          counterVector(3) = counterStr

      End If

      ! Broadcasting all sizes of the different types of data
      Call MPI_Bcast(counterVector, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

      ! Reconstructing the counters locally and allocating local arrays
      counterInt = counterVector(1)
      counterDbl = counterVector(2)
      counterStr = counterVector(3)

      Allocate(vectorInt_value(1:counterInt))
      Allocate(vectorDbl_value(1:counterDbl))
      Allocate(vectorStr_value(1:counterStr))

      ! Allocating a local array for data names
      counterTot = counterInt + counterDbl + counterStr
      Allocate(vectorAll_names(1:counterTot))

      ! Only processor 0 can fill out the vectors
      If(mpi_rank.Eq.0) Then

         ! Now that the final size of these vectors is known, we can
         ! allocate an array and broadcast to everyone else


         ! First, vector of integer data, and the relevant keywords

         Allocate(currentNodeInt)

         currentNodeInt => listeInt

         ii = 0

         Do While ( associated(currentNodeInt) )

            ii = ii + 1

            vectorInt_value(ii) = currentNodeInt%value
            vectorAll_names(ii) = currentNodeInt%name

            Allocate(previousNodeInt)

            previousNodeInt => currentNodeInt
            currentNodeInt  => currentNodeInt%next

            Deallocate(previousNodeInt)

         End Do

         ! Second, vector of real data, and the relevant keywords

         Allocate(currentNodeDbl)

         currentNodeDbl => listeDbl

         ii = 0

         Do While ( associated(currentNodeDbl) )

            ii = ii + 1
            jj = ii + counterInt

            vectorDbl_value(ii) = currentNodeDbl%value
            vectorAll_names(jj) = currentNodeDbl%name

            Allocate(previousNodeDbl)

            previousNodeDbl => currentNodeDbl
            currentNodeDbl  => currentNodeDbl%next

            Deallocate(previousNodeDbl)

         End Do

         ! Third: vector of string data, and the relevant keywords

         Allocate(currentNodeStr)

         currentNodeStr => listeStr

         ii = 0

         Do While ( associated(currentNodeStr) )

            ii = ii + 1
            jj = ii + counterInt + counterDbl

            vectorStr_value(ii) = currentNodeStr%value
            vectorAll_names(jj) = currentNodeStr%name

            Allocate(previousNodeStr)

            previousNodeStr => currentNodeStr
            currentNodeStr  => currentNodeStr%next

            Deallocate(previousNodeStr)

         End Do

      End If ! end of proc = 0 if

      ! Broadcasting all data
      Call MPI_Bcast(vectorInt_value, counterInt, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
      Call MPI_Bcast(vectorDbl_value, counterDbl, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_err)
      Call MPI_Bcast(vectorStr_value, counterStr, mpi_string_10, 0, MPI_COMM_WORLD, mpi_err)

      ! Broadcasting all data names
      Call MPI_Bcast(vectorAll_names, counterTot,  mpi_string_10, 0, MPI_COMM_WORLD, mpi_err)

      ! Releasing derived MPI datatype
      Call mpi_type_free(mpi_string_10, mpi_err)

   End Subroutine mpi_getParallelData

   !-------------------------------------------------------------------!
   ! This subroutine initializes processor-dependent variables used to !
   ! construct the list of jobs. The subroutine is executed by all     !
   ! processors                                                        !
   !-------------------------------------------------------------------!

   Subroutine mpi_initializeParallelData()

      Integer :: iq

      modeMPI_deformation = 1
      modeMPI_basis       = 0

      turnConstraintOn        = 1
      optimizePath            = 0
      batch_mode              = 0
      number_batch            = 1
      number_constraint_batch = 1

      counterForce = 1
      Allocate(forceVector(1:counterForce))
      forceVector(:) = 'SKM*'

      Allocate(lambda_batch(1:number_constraint_batch))
      Allocate(   miu_batch(1:number_constraint_batch))
      lambda_batch(:) = 3
      miu_batch(:)    = 0

      proton_initValue  = 66; proton_stepValue  = 2; proton_stepNumber  = 1
      neutron_initValue = 86; neutron_stepValue = 2; neutron_stepNumber = 2

      x_proton_initValue  = -1.0D0; x_proton_stepValue  = 0.1
      x_neutron_initValue = -1.0D0; x_neutron_stepValue = 0.1

      numberConstraints = 1

      Allocate(qStepNumber(1:numberConstraints))
      Allocate(    qLambda(1:numberConstraints))
      Allocate(       qMiu(1:numberConstraints))
      Allocate( qinitValue(1:numberConstraints))
      Allocate(  qFinValue(1:numberConstraints))

      Do iq = 1, numberConstraints
         qStepNumber(iq) = 4
             qLambda(iq) = 2
                qMiu(iq) = 0
          qinitValue(iq) = 10.0D0
           qFinValue(iq) = 10.0D0
      End Do

      numberConstraints_restart = 1

      Allocate(qStepNumber_restart(1:numberConstraints_restart))
      Allocate(    qLambda_restart(1:numberConstraints_restart))
      Allocate(       qMiu_restart(1:numberConstraints_restart))
      Allocate( qinitValue_restart(1:numberConstraints_restart))
      Allocate(  qFinValue_restart(1:numberConstraints_restart))

      Do iq = 1, numberConstraints_restart
         qStepNumber_restart(iq) = 4
             qLambda_restart(iq) = 2
                qMiu_restart(iq) = 0
          qinitValue_restart(iq) = 10.0D0
           qFinValue_restart(iq) = 10.0D0
      End Do

      basis_frequency_initValue  = 8.0D0
      basis_frequency_stepValue  = 0.1D0
      basis_frequency_stepNumber = 1

      basis_size_initValue  = 8
      basis_size_stepValue  = 2
      basis_size_stepNumber = 1

      basis_states_initValue  = 165
      basis_states_stepValue  = 2
      basis_states_stepNumber = 1

      basis_deform_initValue  = 0.0D0
      basis_deform_stepValue  = 0.1D0
      basis_deform_stepNumber = 1

   End Subroutine mpi_initializeParallelData

   !-------------------------------------------------------------------!
   ! This subroutine sets variables that are used to construct the     !
   ! list of jobs. The subroutine is executed by all processors        !
   !-------------------------------------------------------------------!
   Subroutine mpi_setParallelData()

      Integer :: ii,jj,iq,iq_batch

      Character(Len=10) :: KEYWOR, string_10

      !=================================================!
      !                                                 !
      !   RECONSTRUCTION OF THE INPUT OF TYPE INTEGER   !
      !                                                 !
      !=================================================!

      ii = 1

      Do While ( ii .le. counterInt )

         KEYWOR = vectorAll_names(ii)

         If(trim(KEYWOR).Eq.'CALCULMODE') Then

            modeMPI_basis = vectorInt_value(ii)

            ii = ii + 1
            modeMPI_deformation = vectorInt_value(ii)

         End If

         If(trim(KEYWOR).Eq.'CONSTR_DEF'.OR.trim(KEYWOR).Eq.'CONSTR-DEF') Then
            turnConstraintOn = vectorInt_value(ii)
         End If

         If(trim(KEYWOR).Eq.'BATCH_MODE'.OR.trim(KEYWOR).Eq.'BATCH-MODE') Then

            number_batch = vectorInt_value(ii)

            ii = ii + 1
            batch_mode = vectorInt_value(ii)

            If(batch_mode.Eq.0) number_batch=1

         End If

         If(trim(KEYWOR).Eq.'BATCH_SPEC'.OR.trim(KEYWOR).Eq.'BATCH-SPEC') Then

            number_constraint_batch = vectorInt_value(ii)

            If(Allocated(lambda_batch)) Then
               Deallocate(lambda_batch,miu_batch)
            End If

            Allocate(lambda_batch(1:number_constraint_batch))
            Allocate(   miu_batch(1:number_constraint_batch))

            Do iq_batch = 1, number_constraint_batch

               ii = ii + 1
               miu_batch(iq_batch) = vectorInt_value(ii)

               ii = ii + 1
               lambda_batch(iq_batch) = vectorInt_value(ii)

             End Do

         End If

         If(trim(KEYWOR).Eq.'ALL_NUCLEI'.OR.trim(KEYWOR).Eq.'ALL-NUCLEI') Then

            neutron_stepNumber = vectorInt_value(ii)

            ii = ii + 1
            neutron_stepValue = vectorInt_value(ii)

            ii = ii + 1
            neutron_initValue = vectorInt_value(ii)

            ii = ii + 1
            proton_stepNumber = vectorInt_value(ii)

            ii = ii + 1
            proton_stepValue = vectorInt_value(ii)

            ii = ii + 1
            proton_initValue = vectorInt_value(ii)

         End If

         If(trim(KEYWOR).Eq.'REALNUCLEI') Then

            neutron_stepNumber = vectorInt_value(ii)

            ii = ii + 1
            proton_stepNumber = vectorInt_value(ii)

         End If

         If(trim(KEYWOR).Eq.'MULTICONST') Then

            numberConstraints = vectorInt_value(ii)

            If(Allocated(qLambda)) Then
               Deallocate(qLambda,qMiu,qStepNumber)
            End If

            Allocate(qLambda(1:numberConstraints))
            Allocate(qMiu(1:numberConstraints))
            Allocate(qStepNumber(1:numberConstraints))

            Do iq = 1, numberConstraints

               ii = ii + 1
               qStepNumber(iq) = vectorInt_value(ii)

               ii = ii + 1
               qMiu(iq) = vectorInt_value(ii)

               ii = ii + 1
               qLambda(iq) = vectorInt_value(ii)

            End Do

         End If

         If(trim(KEYWOR).Eq.'MULTIRESTA') Then

            numberConstraints_restart = vectorInt_value(ii)

            If(Allocated(qLambda_restart)) Then
            Deallocate(qLambda_restart,qMiu_restart,qStepNumber_restart)
            End If

            Allocate(    qLambda_restart(1:numberConstraints_restart))
            Allocate(       qMiu_restart(1:numberConstraints_restart))
            Allocate(qStepNumber_restart(1:numberConstraints_restart))

            Do iq = 1, numberConstraints_restart

                ii = ii + 1
                qStepNumber_restart(iq) = vectorInt_value(ii)

                ii = ii + 1
                qMiu_restart(iq) = vectorInt_value(ii)

                ii = ii + 1
                qLambda_restart(iq) = vectorInt_value(ii)

             End Do

         End If

         If(trim(KEYWOR).Eq.'OPTIM_GRID'.OR.trim(KEYWOR).Eq.'OPTIM-GRID') Then
            optimizePath = vectorInt_value(ii)
         End If

         If(trim(KEYWOR).Eq.'ALL_FORCES'.OR.trim(KEYWOR).Eq.'ALL-FORCES') Then
            counterForce = vectorInt_value(ii)
         End If

         If(trim(KEYWOR).Eq.'BASIS_FREQ'.OR.trim(KEYWOR).Eq.'BASIS-FREQ') Then
            basis_frequency_stepNumber = vectorInt_value(ii)
         End If

         If(trim(KEYWOR).Eq.'BASIS_NSHL'.OR.trim(KEYWOR).Eq.'BASIS-NSHL') Then

            basis_size_stepNumber = vectorInt_value(ii)

            ii = ii + 1
            basis_size_stepValue = vectorInt_value(ii)

            ii = ii + 1
            basis_size_initValue = vectorInt_value(ii)

         End If

         If(trim(KEYWOR).Eq.'BASIS_NSTA'.OR.trim(KEYWOR).Eq.'BASIS-NSTA') Then

            basis_states_stepNumber = vectorInt_value(ii)

            ii = ii + 1
            basis_states_stepValue = vectorInt_value(ii)

            ii = ii + 1
            basis_states_initValue = vectorInt_value(ii)

         End If

         If(trim(KEYWOR).Eq.'BASIS_DEFS'.OR.trim(KEYWOR).Eq.'BASIS-DEFS') Then
            basis_deform_stepNumber = vectorInt_value(ii)
         End If


         ! Go to the next item on the list
         ii = ii + 1

      End Do

      Deallocate(vectorInt_value)

      !==============================================!
      !                                              !
      !   RECONSTRUCTION OF THE INPUT OF TYPE REAL   !
      !                                              !
      !==============================================!

      ii = 1

      Do While ( ii .le. counterDbl )

         KEYWOR = vectorAll_names(ii + counterInt)

         If(trim(KEYWOR).Eq.'REALNUCLEI') Then

            x_neutron_stepValue = vectorDbl_value(ii)

            ii = ii + 1
            x_neutron_initValue = vectorDbl_value(ii)

            ii = ii + 1
            x_proton_stepValue = vectorDbl_value(ii)

            ii = ii + 1
            x_proton_initValue = vectorDbl_value(ii)

         End If

         If(trim(KEYWOR).Eq.'MULTICONST') Then

            If(Allocated(qinitValue)) Then
               Deallocate(qinitValue,qFinValue)
            End If

            Allocate(qinitValue(1:numberConstraints))
            Allocate(qFinValue(1:numberConstraints))

            Do iq = 1, numberConstraints

               qFinValue(iq) = vectorDbl_value(ii)
               ii = ii + 1

               qinitValue(iq) = vectorDbl_value(ii)
               ii = ii + 1

            End Do
            ii = ii - 1

         End If

         If(trim(KEYWOR).Eq.'MULTIRESTA') Then

            If(Allocated(qinitValue_restart)) Then
               Deallocate(qinitValue_restart,qFinValue_restart)
            End If

            Allocate(qinitValue_restart(1:numberConstraints_restart))
            Allocate( qFinValue_restart(1:numberConstraints_restart))

            Do iq = 1, numberConstraints_restart

               qFinValue_restart(iq) = vectorDbl_value(ii)
               ii = ii + 1

               qinitValue_restart(iq) = vectorDbl_value(ii)
               ii = ii + 1

            End Do
            ii = ii - 1

         End If

         If(trim(KEYWOR).Eq.'BASIS_FREQ'.OR.trim(KEYWOR).Eq.'BASIS-FREQ') Then

            basis_frequency_stepValue = vectorDbl_value(ii)

            ii = ii + 1
            basis_frequency_initValue = vectorDbl_value(ii)

         End If

         If(trim(KEYWOR).Eq.'BASIS_DEFS'.OR.trim(KEYWOR).Eq.'BASIS-DEFS') Then

            basis_deform_stepValue = vectorDbl_value(ii)

            ii = ii + 1
            basis_deform_initValue = vectorDbl_value(ii)

         End If


         ! Go to the next item on the list
         ii = ii + 1

      End Do

      Deallocate(vectorDbl_value)

      !================================================!
      !                                                !
      ! RECONSTRUCTION OF THE INPUT OF TYPE STRING*10  !
      !                                                !
      !================================================!

      ii = 1

      Do While ( ii .le. counterStr )

         KEYWOR = vectorAll_names(ii + counterInt + counterDbl)

         If(trim(KEYWOR).Eq.'ALL_FORCES'.OR.trim(KEYWOR).Eq.'ALL-FORCES') Then

            If(Allocated(forceVector)) Then
               Deallocate(forceVector)
            End If

            Allocate(forceVector(1:counterForce))

            Do jj=1, counterForce
               string_10 = vectorStr_value(jj)
               forceVector(jj) = string_10(1:4)
            End Do

         End If

         ! Go to the next item on the list
         ii = ii + counterForce

      End Do

      Deallocate(vectorStr_value)

      Deallocate(vectorAll_names)

   End Subroutine mpi_setParallelData

   !-------------------------------------------------------------------!
   ! This subroutine prints variables that are common to all processes.!
   ! The subroutine is executed by processor 0 only                    !
   !-------------------------------------------------------------------!
   Subroutine mpi_printSequentialData(IOPTGS,NXHERM,NYHERM,NZHERM, &
                                      IPAIRI,IPAHFB,IPABCS,IMFHFB, &
                                                    LIPKIN,LIPKIP, &
                                      FE2FIN,IF2FIN,FE2FIP,IF2FIP, &
                                                    IROTAT,IREVIE, &
                                      ICONTI,IPCONT,IGCONT,IGPCON, &
                                             ILCONT,IFCONT,IACONT, &
                                      NUCHAO,NUPING,EPSPNG,EPSITE, &
                                             ISIMPY,ISIGNY,IPARTY, &
                                             ISIMTX,ISIMTY,ISIMTZ, &
                                                    NOITER,NULAST, &
                                             NOSCIL,NLIMIT,ENECUT, &
                                             ICOTYP,ICOUDI,ICOUEX, &
                                                    I_GOGA,IGOGPA, &
                                                           I1LINE, &
                                      SLOWEV,SLOWOD,SLOWPA,SLOWLI, &
                                                           ECUTOF, &
                                                           NILXYZ, &
                                      IBROYD,N_ITER,BROTRI,ALPHAM, &
                                             PRHO_T,PRHOST,POWERT, &
                                      PRHO_N,PRHOSN,PRHODN,POWERN, &
                                      PRHO_P,PRHOSP,PRHODP,POWERP, &
                                                           OMEGAY, &
                               IFSHEL,NPOLYN,GSTRUN,GSTRUP,HOMFAC, &
                                      LAMBAS,MIUBAS,BETBAS,NUMDEF, &
                               LAM_WS,MIU_WS,NUM_WS,ITE_WS,STP_WS, &
                                             LMXINP,IFPUSH,IREAWS, &
                               ISTAND,KETA_J,KETA_W,KETACM,KETA_M, &
                                      IF_RPA,IF_EDF,IF_THO,CBETHO, &
                                                    IBASIS,IFRAGM, &
                                             NMUCON,NMUCOU,NMUPRI, &
                                             NEXBET,IPRIBE,IPRIBL, &
                                                           FCHOM0, &
                                             IWRIFI,IFTEMP,IFNECK, &
                                             INNUMB,IZNUMB,R0PARM, &
                               DELTAE,XLOCMX,V2_MIN,NTHETA,MIN_QP, &
                                             IDEALL,IDELOC,IDECON, &
                                                           TEMP_T, &
                                                       batch_init)

      Integer, intent(inout) :: ICOTYP, ICOUDI, ICOUEX, I_GOGA, IGOGPA, NOITER, IBROYD, &
                                N_ITER, NULAST, NUPING, NUCHAO, ISIMPY, ISIGNY, IPARTY, &
                                IROTAT, ISIMTX, ISIMTY, ISIMTZ, IPAIRI, IPAHFB, IPABCS, &
                                IMFHFB, LIPKIN, LIPKIP, IOPTGS, NXHERM, NYHERM, NZHERM, &
                                NOSCIL, NLIMIT, I1LINE, NILXYZ, IREVIE, ICONTI, IPCONT, &
                                IGCONT, IGPCON, ILCONT, IFCONT, IACONT, IFSHEL, NPOLYN, &
                                NUMDEF, NUM_WS, ITE_WS, IREAWS, ISTAND, KETA_J, KETA_W, &
                                KETACM, KETA_M, IF_RPA, IF_EDF, IF_THO, IBASIS, IFRAGM, &
                                IFPUSH, NMUCON, NMUCOU, NMUPRI, NEXBET, IPRIBE, IPRIBL, &
                                IWRIFI, IFTEMP, IFNECK, LMXINP, INNUMB, IZNUMB, NTHETA, &
                                MIN_QP, IDEALL, IDELOC, IDECON, batch_init,IF2FIN,IF2FIP

      Real, intent(inout) :: BROTRI, SLOWEV, SLOWOD, SLOWPA, SLOWLI, EPSITE, &
                             EPSPNG, ECUTOF, ENECUT, PRHO_T, PRHOST, POWERT, &
                             OMEGAY, GSTRUN, GSTRUP, HOMFAC, STP_WS, ALPHAM, &
                             CBETHO, FCHOM0, PRHO_N, PRHOSN, PRHODN, POWERN, &
                             PRHO_P, PRHOSP, PRHODP, POWERP, R0PARM, DELTAE, &
                             XLOCMX, V2_MIN, FE2FIN, FE2FIP, TEMP_T

      Integer :: ii,nfisum,LAMACT,MIUACT
      Integer :: mpi_rank,mpi_err

      Integer, Dimension(:), pointer :: LAMBAS,MIUBAS,LAM_WS,MIU_WS

      Real, Dimension(:), pointer :: BETBAS

      COMMON /OUTPUT/ nfisum

      Call mpi_comm_rank(MPI_COMM_WORLD,mpi_rank,mpi_err)

      ! Printing out the characteristics of the parallel run
      If(mpi_rank.Eq.0 .And. batch_init.Eq.1) Then

         Write(nfisum,'(/,14X,"/",42("-"),"/")')
         Write(nfisum,'(14X,"/ CHARACTERISTICS OF THE RUN: SERIAL PART  /")')
         Write(nfisum,'(14X,"/",42("-"),"/",/)')

         Write(nfisum,'("/",77("-"),"/",/,"/",77X,"/")')
         Write(nfisum,'("/   Iterations",64X,"/")')
         Write(nfisum,'("/   ----------",64X,"/")')
         Write(nfisum,'("/     IBROYD = ",i1," N_ITER = ",i3," BROTRI = ",f9.3," ALPHAM = ",f7.4,13x, &
        &               "/",/,"/",77X,"/")') IBROYD,N_ITER,BROTRI,ALPHAM
         Write(nfisum,'("/     NOITER = ",i5," NULAST = ",i5," NUCHAO = ",i5," NUPING = ",i5,13X,"/")') &
                          NOITER,NULAST,NUCHAO,NUPING
         Write(nfisum,'("/     SLOWEV = ",f7.4," SLOWOD = ",f7.4," SLOWPA = ",f7.4," SLOWLI = ",f7.4,5x,"/")') &
                          SLOWEV,SLOWOD,SLOWPA,SLOWLI
         Write(nfisum,'("/     EPSPNG = ",e20.12," EPSITE = ",e20.12,13X,"/")') &
                          EPSPNG,EPSITE
         Write(nfisum,'("/",77X,"/")')

         Write(nfisum,'("/   Characteristics of the HO basis",43X,"/")')
         Write(nfisum,'("/   -------------------------------",43X,"/")')
         Write(nfisum,'("/     IBASIS = ",i3," IOPTGS = ",i1,49X,"/")') IBASIS,IOPTGS
         Write(nfisum,'("/     NOSCIL = ",i3," NLIMIT = ",i5," ENECUT = ",f9.3,26X,"/")') &
                          NOSCIL,NLIMIT,ENECUT
         Write(nfisum,'("/     NXHERM = ",i3," NYHERM = ",i3," NZHERM = ",i3,34X,"/")') &
                          NXHERM,NYHERM,NZHERM
         Write(nfisum,'("/     INNUMB = ",i3," IZNUMB = ",i5," R0PARM = ",f9.3,26X,"/")') &
                          INNUMB,IZNUMB,R0PARM
         Write(nfisum,'("/",77X,"/")')

         If(Associated(LAMBAS)) Then
            Write(nfisum,'("/     Deformations of the nuclear surface",37X,"/")')
            Do ii=1,NUMDEF
               LAMACT=Abs(LAMBAS(ii)); MIUACT=MIUBAS(ii)
               Write(nfisum,'("/     l = ",i2," m = ",i2," beta_lm = ",f7.4,41X,"/")') &
                          LAMACT,MIUACT,BETBAS(ii)
            End Do
            Write(nfisum,'("/",77X,"/")')
         End If

         Write(nfisum,'("/   Type of calculation",55X,"/")')
         Write(nfisum,'("/   -------------------",55X,"/")')
         Write(nfisum,'("/     IPAIRI = ",i1," IPAHFB = ",i1," IPABCS =",i2, &
        &                    " IMFHFB = ",i1," LIPKIN = ",i1," LIPKIP = ",i1,7X,"/")') &
                          IPAIRI,IPAHFB,IPABCS,IMFHFB,LIPKIN,LIPKIP
         Write(nfisum,'("/     FE2FIN = ",f9.3,"    IF2FIN = ",i1," FE2FIP = ",f9.3,"    IF2FIP = ",i1,7X,"/")') &
                               FE2FIN,IF2FIN,FE2FIP,IF2FIP
         Write(nfisum,'("/     ICOUDI = ",i1," ICOUEX = ",i1," I_GOGA = ",i1, &
        &                    " IGOGPA = ",i1," IFTEMP = ",i1," IFSHEL = ",i1,7X,"/")') &
                          ICOUDI,ICOUEX,I_GOGA,IGOGPA,IFTEMP,IFSHEL
         Write(nfisum,'("/     IFRAGM = ",i1," IF_THO = ",i1," TEMP_T=",f5.2,38X,"/")') &
                               IFRAGM,IF_THO,TEMP_T
         Write(nfisum,'("/     IF_EDF = ",i1," ISTAND = ",i1," KETA_J = ",i1, &
        &                    " KETA_W = ",i1," KETACM = ",i1," KETA_M = ",i1,7X,"/")') &
                          IF_EDF,ISTAND,KETA_J,KETA_W,KETACM,KETA_M
         Write(nfisum,'("/",77X,"/")')

         Write(nfisum,'("/   Pairing options",59X,"/")')
         Write(nfisum,'("/   ---------------",59X,"/")')
         Write(nfisum,'("/     ECUTOF = ",f9.3,54X,"/")') ECUTOF
         Write(nfisum,'("/     PRHO_T = ",f9.3," PRHOST = ",f9.3," POWERT = ",f9.3,16X,"/")') &
                          PRHO_T,PRHOST,POWERT
         Write(nfisum,'("/     PRHO_N = ",f9.3," PRHOSN = ",f9.3," POWERN = ",f9.3,16X,"/")') &
                          PRHO_N,PRHOSN,POWERN
         Write(nfisum,'("/     PRHO_P = ",f9.3," PRHOSP = ",f9.3," POWERP = ",f9.3,16X,"/")') &
                          PRHO_P,PRHOSP,POWERP
         Write(nfisum,'("/",77X,"/")')

         Write(nfisum,'("/   Constraints",63X,"/")')
         Write(nfisum,'("/   -----------",63X,"/")')
         Write(nfisum,'("/     IF_RPA = ",i2," IFPUSH = ",i2," IFNECK = ",i2,37X,"/")') &
                               IF_RPA,IFPUSH,IFNECK
         Write(nfisum,'("/     NMUCON = ",i2," NMUCOU = ",i2," NMUPRI = ",i2,37X,"/")') &
                          NMUCON,NMUCOU,NMUPRI
         Write(nfisum,'("/     OMEGAY = ",f7.4,56X,"/")') &
                          OMEGAY
         Write(nfisum,'("/",77X,"/")')

         Write(nfisum,'("/   Restart options",59X,"/")')
         Write(nfisum,'("/   ---------------",59X,"/")')
         Write(nfisum,'("/     ICONTI = ",i1," IPCONT = ",i1," IGCONT = ",i1, &
        &                    " IGPCON = ",i1," ILCONT = ",i1," IFCONT = ",i1,7X,"/")') &
                          ICONTI,IPCONT,IGCONT,IGPCON,ILCONT,IFCONT
         Write(nfisum,'("/     IREAWS = ",i1,62X,"/")') IREAWS
         If(IREAWS.Eq.1) Then
            Write(nfisum,'("/     Minimization over deformations",42X,"/")')
            Write(nfisum,'("/        NUMINP = ",i3," ITEINP = ",i3," LMXINP = ",i2," STPINP = ",f7.4,15X,"/")') &
                                     NUM_WS,ITE_WS,LMXINP,STP_WS
            Do ii=1,NUM_WS
               LAMACT=Abs(LAM_WS(ii)); MIUACT=MIU_WS(ii)
               Write(nfisum,'("/     l = ",i2," m = ",i2,59X,"/")') LAMACT,MIUACT
            End Do
         End If
         Write(nfisum,'("/",77X,"/")')

         Write(nfisum,'("/   Symmetry flags",60X,"/")')
         Write(nfisum,'("/   --------------",60X,"/")')
         Write(nfisum,'("/     IROTAT = ",i2," ISIMPY = ",i2," ISIGNY = ",i2," IPARTY = ",i2,25X,"/")') &
                          IROTAT,ISIMPY,ISIGNY,IPARTY
         Write(nfisum,'("/                 ISIMTX = ",i2," ISIMTY = ",i2," ISIMTZ = ",i2,25X,"/")') &
                          ISIMTX,ISIMTY,ISIMTZ
         Write(nfisum,'("/",77X,"/")')

         Write(nfisum,'("/   Various",67X,"/")')
         Write(nfisum,'("/   -------",67X,"/")')
         Write(nfisum,'("/     I1LINE = ",i2," NILXYZ = ",i2," IREVIE = ",i2," IWRIFI = ",i2,25X,"/")') &
                          I1LINE,NILXYZ,IREVIE,IWRIFI
         Write(nfisum,'("/     NEXBET = ",i2," IPRIBE = ",i2," IPRIBL = ",i2,37X,"/")') &
                          NEXBET,IPRIBE,IPRIBL
         If(IFSHEL.Gt.0) then
            Write(nfisum,'("/",77X,"/")')
            Write(nfisum,'("/     NPOLYN = ",i2," HOMFAC = ",f7.4," GSTRUN = ",f7.4," GSTRUP = ",f7.4,10X,"/")') &
                             NPOLYN,HOMFAC,GSTRUN,GSTRUP
         End if
         Write(nfisum,'("/     MIN_QP = ",i2," DELTAE = ",f7.3," XLOCMX = ",f7.5," V2_MIN = ",f7.5,10X,"/")') &
                               MIN_QP,DELTAE,XLOCMX,V2_MIN
         Write(nfisum,'("/     IDEALL = ",i2," IDELOC = ",i2," IDECON = ",i2,37X,"/")') &
                               IDEALL,IDELOC,IDECON
         Write(nfisum,'("/",77X,"/",/,"/",77("-"),"/")')

      End If

      Call mpi_barrier(MPI_COMM_WORLD, mpi_err)

   End Subroutine mpi_printSequentialData

   !-------------------------------------------------------------------!
   ! This subroutine prints variables that are used to construct the   !
   ! list of jobs. The subroutine is executed by processor 0 only      !
   !-------------------------------------------------------------------!
   Subroutine mpi_printParallelData()

      Integer :: ii,nfisum,mpi_rank,mpi_err

      COMMON /OUTPUT/ nfisum

      Call mpi_comm_rank(MPI_COMM_WORLD,mpi_rank,mpi_err)

      ! Printing out the characteristics of the parallel run
      If(mpi_rank.Eq.0) Then

         Write(nfisum,'(/,14X,"/",44("-"),"/")')
         Write(nfisum,'(14X,"/ CHARACTERISTICS OF THE RUN: PARALLEL PART  /")')
         Write(nfisum,'(14X,"/",44("-"),"/",/)')

         Write(nfisum,'("/",77("-"),"/",/,"/",77X,"/")')
         Write(nfisum,'("/   Calculation mode",58X,"/")')
         Write(nfisum,'("/   ----------------",58X,"/")')
         Write(nfisum,'("/     Deformation .: ",i2,55X,"/")') modeMPI_deformation
         Write(nfisum,'("/     Batch .......: ",i2,55X,"/")') batch_mode
         Write(nfisum,'("/     Basis .......: ",i2,55X,"/")') modeMPI_basis
         Write(nfisum,'("/",77X,"/")')

         Write(nfisum,'("/   Nucleus grid",62X,"/")')
         Write(nfisum,'("/   ------------",62X,"/")')
         Write(nfisum,'("/     Particle   Initial    Step   Number of steps",28X,"/")')
         If(x_neutron_initValue.Gt.0.0D0) Then
            Write(nfisum,'("/      neutron   ",F8.4,"  ",F7.4,"      ",i4,34X,"/")') &
                         x_neutron_initValue, x_neutron_stepValue, neutron_stepNumber
         Else
            Write(nfisum,'("/      neutron     ",i4,"      ",i2,"         ",i4,34X,"/")') &
                         neutron_initValue, neutron_stepValue, neutron_stepNumber
         End If
         If(x_proton_initValue.Gt.0.0D0) Then
            Write(nfisum,'("/      proton    ",F8.4,"  ",F7.4,"      ",i4,34X,"/")') &
                         x_proton_initValue, x_proton_stepValue, proton_stepNumber
         Else
            Write(nfisum,'("/      proton      ",i4,"      ",i2,"         ",i4,34X,"/")') &
                         proton_initValue, proton_stepValue, proton_stepNumber
         End If
         Write(nfisum,'("/",77X,"/")')

         ! Displaying characteristics of the default deformation grid (regular)
         If(batch_mode.Eq.1) Then
            Write(nfisum,'("/   Batch mode grid",59X,"/")')
            Write(nfisum,'("/   ---------------",59X,"/")')
            Write(nfisum,'("/     Number of batches ......: ",i4,42X,"/")') number_batch
            Write(nfisum,'("/     Number of constraints ..: ",i4,42X,"/")') number_constraint_batch
            Write(nfisum,'("/     Constraint #    Lambda   Miu",44X,"/")')
            Do ii=1,number_constraint_batch
               Write(nfisum,'("/     ",6X,i2,9X,i2,6X,i2,45X,"/")') &
                       ii, lambda_batch(ii), miu_batch(ii)
            End Do
            Write(nfisum,'("/",77X,"/")')
         End If

         ! Displaying characteristics of the default deformation grid (regular)
         If(modeMPI_deformation.Eq.1.Or.modeMPI_basis.Eq.1) Then
            Write(nfisum,'("/   Deformation grid",58X,"/")')
            Write(nfisum,'("/   ----------------",58X,"/")')
            Write(nfisum,'("/     Constrained calculations: ",i4,42X,"/")') turnConstraintOn
            Write(nfisum,'("/     Number of constraints ..: ",i4,42X,"/")') numberConstraints
            Write(nfisum,'("/     Constraint #    Lambda   Miu   Initial    Final   Number of steps",7X,"/")')
            Do ii=1,numberConstraints
               Write(nfisum,'("/     ",6X,i2,9X,i2,6X,i2,2X,f8.2,2X,f8.2,6X,i5,14X,"/")') &
                       ii, qLambda(ii), qMiu(ii), qinitValue(ii), qFinValue(ii), qStepNumber(ii)
            End Do
            Write(nfisum,'("/",77X,"/")')
         End If

         ! Displaying characteristics of the initial grid (regular) and restart deformation grid (regular)
         If(modeMPI_deformation.Eq.2) Then
            Write(nfisum,'("/   Original deformation grid",49X,"/")')
            Write(nfisum,'("/   -------------------------",49X,"/")')
            Write(nfisum,'("/     Constrained calculations: ",i4,42X,"/")') turnConstraintOn
            Write(nfisum,'("/     Number of constraints ..: ",i4,42X,"/")') numberConstraints
            Write(nfisum,'("/     Constraint #    Lambda   Miu   Initial    Final   Number of steps",7X,"/")')
            Do ii=1,numberConstraints
               Write(nfisum,'("/     ",6X,i2,9X,i2,6X,i2,2X,f8.2,2X,f8.2,6X,i5,14X,"/")') &
                       ii, qLambda(ii), qMiu(ii), qinitValue(ii), qFinValue(ii), qStepNumber(ii)
            End Do
            Write(nfisum,'("/   New deformation grid",54X,"/")')
            Write(nfisum,'("/   --------------------",54X,"/")')
            Write(nfisum,'("/     Optimize path ..........: ",i4,42X,"/")') optimizePath
            Write(nfisum,'("/     Number of constraints ..: ",i4,42X,"/")') numberConstraints_restart
            Write(nfisum,'("/     Constraint #    Lambda   Miu   Initial    Final   Number of steps",7X,"/")')
            Do ii=1,numberConstraints_restart
               Write(nfisum,'("/     ",6X,i2,9X,i2,6X,i2,2X,f8.2,2X,f8.2,6X,i5,14X,"/")') &
                       ii, qLambda_restart(ii), qMiu_restart(ii), &
                           qinitValue_restart(ii), qFinValue_restart(ii), qStepNumber_restart(ii)
            End Do
            Write(nfisum,'("/",77X,"/")')
         End If

         ! Displaying characteristics of the initial grid (only what constraints were used) and of the
         ! restart grid (regular)
         If(modeMPI_deformation.Eq.3) Then
            Write(nfisum,'("/   Original deformation grid",49X,"/")')
            Write(nfisum,'("/   -------------------------",49X,"/")')
            Write(nfisum,'("/     Number of constraints ..: ",i4,42X,"/")') numberConstraints
            Write(nfisum,'("/     Constraint #    Lambda   Miu",44X,"/")')
            Do ii=1,numberConstraints
               Write(nfisum,'("/     ",6X,i2,9X,i2,6X,i2,45X,"/")') &
                       ii, qLambda(ii), qMiu(ii)
            End Do
            Write(nfisum,'("/   New deformation grid",54X,"/")')
            Write(nfisum,'("/   --------------------",54X,"/")')
            Write(nfisum,'("/     Optimize path ..........: ",i4,42X,"/")') optimizePath
            Write(nfisum,'("/     Number of constraints ..: ",i4,42X,"/")') numberConstraints_restart
            Write(nfisum,'("/     Constraint #    Lambda   Miu   Initial    Final   Number of steps",7X,"/")')
            Do ii=1,numberConstraints_restart
               Write(nfisum,'("/     ",6X,i2,9X,i2,6X,i2,2X,f8.2,2X,f8.2,6X,i5,14X,"/")') &
                       ii, qLambda_restart(ii), qMiu_restart(ii), &
                           qinitValue_restart(ii), qFinValue_restart(ii), qStepNumber_restart(ii)
            End Do
            write(nfisum,'("/",77X,"/")')
         End If

         ! Displaying characteristics of the origina deformation grid (grid). New calculation will not
         ! yield a regular grid but an irregular path defined by what is in the file hfodd_path.d
         If(modeMPI_deformation.Eq.4) Then
            Write(nfisum,'("/   Original deformation grid",49X,"/")')
            Write(nfisum,'("/   -------------------------",49X,"/")')
            Write(nfisum,'("/     Number of constraints ..: ",i4,42X,"/")') numberConstraints
            Write(nfisum,'("/     Constraint #    Lambda   Miu",44X,"/")')
            Do ii=1,numberConstraints
               Write(nfisum,'("/     ",6X,i2,9X,i2,6X,i2,45X,"/")') &
                       ii, qLambda(ii), qMiu(ii)
            End Do
            write(nfisum,'("/",77X,"/")')
         End If

         Write(nfisum,'("/   Forces",68X,"/")')
         Write(nfisum,'("/   ------",68X,"/")')
         Do ii=1,counterForce
            Write(nfisum,'("/     ",A4,68X,"/")') forceVector(ii)
         End Do
         write(nfisum,'("/",77X,"/")')

         ! Displaying characteristics of the basis grid
         If(modeMPI_basis.Eq.1) Then
            Write(nfisum,'("/   Basis grid",64X,"/")')
            Write(nfisum,'("/   ----------",64X,"/")')
            Write(nfisum,'("/        Type        Initial    Step   Number of steps",24X,"/")')
            Write(nfisum,'("/       Shells        ",i2,"        ",i2,"          ",i2,32X,"/")') &
                      basis_size_initValue,basis_size_stepValue,basis_size_stepNumber
            Write(nfisum,'("/       States        ",i4,"      ",i4,"        ",i2,32X,"/")') &
                      basis_states_initValue,basis_states_stepValue,basis_states_stepNumber
            Write(nfisum,'("/     HO frequency  ",f6.3,"    ",f6.3,"        ",i2,32X,"/")') &
                      basis_frequency_initValue,basis_frequency_stepValue,basis_frequency_stepNumber
            Write(nfisum,'("/    HO deformation ",f6.3,"    ",f6.3,"        ",i2,32X,"/")') &
                      basis_deform_initValue,basis_deform_stepValue,basis_deform_stepNumber
            Write(nfisum,'("/",77X,"/")')
         End If

         Write(nfisum,'("/",77("-"),"/")')

      End If

      Call mpi_barrier(MPI_COMM_WORLD, mpi_err)

   End Subroutine mpi_printParallelData

 End Module hfodd_mpiio

