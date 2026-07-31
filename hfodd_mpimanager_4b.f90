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

 !----------------------------------------------------------------------
 !>
 !>             MPI Management Environment for HFODD
 !>
 !> Role of this module
 !>    - defining the list of tasks to be run based on sequential and
 !>      parallel inputs read in hfodd_mpiio.f90
 !>    - managing the load balancing of the tasks
 !>
 !> Values of modeMPI_deformation
 !>   - 1: computes a regularly-spaced, rectangular grid of multipole
 !>        moments only
 !>   - 2: computes a regularly-spaced, rectangular grid of multipole
 !>        moments only by restarting from a regularly-spaced, rectangular
 !>        grid of multipole moments assuming the same number of
 !>        collective variables or more, but not less
 !>   - 3: computes an irregularly spaced non-rectangular grid of
 !>        arbitrary collective variables that can include multipole
 !>        moments, neck constraints, LN parameters, distance between
 !>        prefragments, etc. This grid is defined in the file
 !>        hfodd_path.d. It has the same format as described above. The
 !>        following conventions apply:
 !>            l=0,  m=0 ..: constraint on the neck <Qn>
 !>            l>=2, m>=0 .: usual multipole moments <Qlm> of compound
 !>                          nucleus
 !>            l=-1, m=-1 .: constraint on lambda_2 for protons
 !>            l=-1, m=+1 .: constraint on lambda_2 for neutrons
 !>            l=-2, m=0 ..: distance between the two centers of mass of
 !>                          the fragments (normalized), \xi
 !>            l=-3, m=0 ..: |A1-A2|/A, with A1, A2 the masses of the
 !>                          prefragments, A the mass of the compound
 !>                          nucleus
 !>   - 4: computes an irregularly spaced non-rectangular grid of
 !>        arbitrary collective variables as for modeMPI_deformation=5,
 !>        only:
 !>         + the configurations are defined in the file hfodd_path_new.d,
 !>         + the calculations are restarted from the configurations
 !>           listed in hfodd_path.d, with the convention that each line
 !>           in hfodd_path.d corresponds to point n and HFODD record file
 !>           named HFODD_[n].REC (n encoded in 8 digits padded with 0s)
 !----------------------------------------------------------------------
 Module hfodd_mpimanager

   Use hfodd_mpiio

   Implicit None

   ! Public
   Integer :: NOFJOB

   ! Private
   Logical :: filled

   Integer :: debug = 0

   Integer :: numberQ=4,numberQ_restart=4,numberQ_orig=4,numberConstraints_orig=1
   Integer, Allocatable :: IZ_JOB(:),IN_JOB(:),indexQjob(:),indexForce(:),indexRestartJob(:)
   Integer, Allocatable :: indexSizeJob(:),indexStatesJob(:),indexLengthJob(:),idxDeformJob(:),    &
                           basis_size(:),basis_states(:)
   Integer, Allocatable :: flex_lambda(:),flex_miu(:)
   Integer, Allocatable :: mapping_qlm(:)
   Integer, Allocatable :: lambdaJob(:,:),miuJob(:,:),lambdaJob_restart(:,:),miuJob_restart(:,:)
   Integer, Allocatable :: table_fichiers(:)

   Real, Allocatable :: XIZJOB(:),XINJOB(:)
   Real, Allocatable :: flex_qlm(:),basis_freq(:),basis_defo(:)
   Real, Allocatable :: qValueJob(:,:),qValueJob_restart(:,:)

   Character(Len=68) :: fichier_path='hfodd_path.d',fichier_path_new='hfodd_path_new.d'

   Character(Len=68), Allocatable :: fichier_rec_old(:),fichier_rec_new(:),&
                                     fichier_tho_old(:),fichier_tho_new(:),&
                                     fichier_lic_old(:),fichier_lic_new(:),&
                                     fichier_fic_old(:),fichier_fic_new(:)
   Character(Len=68), Allocatable :: fichier_rec_tmp(:),&
                                     fichier_tho_tmp(:),&
                                     fichier_lic_tmp(:),&
                                     fichier_fic_tmp(:)

 Contains

   !----------------------------------------------------------------------
   !>  Read the file 'fichier'
   !----------------------------------------------------------------------
   Subroutine read_path(fichier, restart_grid)

     Character(Len=68), Intent(IN) :: fichier
     Logical, Intent(IN) :: restart_grid

     Integer, Parameter :: ndatin=84
     Integer :: i,i_point,n_cons,n_points,ierr,mpi_rank,mpi_err
     Integer, Allocatable :: l_vector(:,:),m_vector(:,:)
     Real, Allocatable :: q_vector(:,:)

     ! Getting the rank of the current process
     Call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, mpi_err)

     If(mpi_rank.Eq.0) Then

        Open(ndatin,file=fichier,status='old',form='formatted',iostat=ierr)
        If(ierr.Ne.0) Then
           Write(6,'("Error in opening file : ",a80)') fichier
           Stop 'Error in read_path - I/O'
        End if

        ! Read the number of constraints used to define a trajectory
        Read(ndatin,*) n_cons, n_points
        If(debug.Ge.1) Write(6,'("In read_path() - Number of constraints: ",i10," Number of lines: ",i10)') n_cons, n_points
        ! Read trajectory in N-d space (N=number_cons)
        Allocate(q_vector(1:n_points,1:n_cons))
        Allocate(l_vector(1:n_points,1:n_cons),m_vector(1:n_points,1:n_cons))
        i_point=0
        Do While(i_point < n_points)
           i_point=i_point+1
           Read(ndatin,*) (l_vector(i_point,i),m_vector(i_point,i),q_vector(i_point,i),i=1,n_cons)
        End Do

        ! Close file
        Close(ndatin)

    End If ! mpi_rank=0

    ! Broadcast number of constraints and number of lines
    Call MPI_Bcast(n_cons,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
    Call MPI_Bcast(n_points,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
    ! Allocate space for vectors
    If(mpi_rank.Ne.0) Then
       Allocate(q_vector(1:n_points,1:n_cons))
       Allocate(l_vector(1:n_points,1:n_cons),m_vector(1:n_points,1:n_cons))
    End If

    ! Broadcast lambda values
    Call MPI_Bcast(l_vector,n_cons*n_points,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
    ! Broadcast miu values
    Call MPI_Bcast(m_vector,n_cons*n_points,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
    ! Broadcast qlm values
    Call MPI_Bcast(q_vector,n_cons*n_points,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpi_err)

    If(restart_grid) Then
       If(debug.Ge.5) Then
          Write(6,'("In read_path() - Setting the following quantities")')
          Write(6,'("  - numberConstraints_restart")')
          Write(6,'("  - numberQ_restart")')
          Write(6,'("  - lambdaJob_restart")')
          Write(6,'("  - miuJob_restart")')
          Write(6,'("  - qValueJob_restart")')
       End If
       numberConstraints_restart=n_cons
       numberQ_restart=n_points
       Allocate(qValueJob_restart(1:numberQ_restart,1:numberConstraints_restart))
       Allocate(lambdaJob_restart(1:numberQ_restart,1:numberConstraints_restart),&
                miuJob_restart(1:numberQ_restart,1:numberConstraints_restart))
       qValueJob_restart(:,:) = q_vector(:,:)
       lambdaJob_restart(:,:) = l_vector(:,:)
       miuJob_restart(:,:)    = m_vector(:,:)
    Else
       If(debug.Ge.5) Then
          Write(6,'("In read_path() - Setting the following quantities")')
          Write(6,'("  - numberConstraints")')
          Write(6,'("  - numberQ")')
          Write(6,'("  - lambdaJob")')
          Write(6,'("  - miuJob")')
          Write(6,'("  - qValueJob")')
          Write(6,'("mpi_rank = ",i4," n_cons = ",i12," n_points = ",i12)') mpi_rank,n_cons,n_points
       End If
       numberConstraints=n_cons
       numberQ=n_points
       Allocate(qValueJob(1:numberQ,1:numberConstraints))
       Allocate(lambdaJob(1:numberQ,1:numberConstraints),&
                miuJob(1:numberQ,1:numberConstraints))
       qValueJob(:,:) = q_vector(:,:)
       lambdaJob(:,:) = l_vector(:,:)
       miuJob(:,:)    = m_vector(:,:)
    End If

    Deallocate(q_vector,l_vector,m_vector)

   End Subroutine read_path

   !---------------------------------------------------------------------
   !>  Defines default deformation grid for modeMPI_deformation = 1, 3
   !>
   !>  Outputs:
   !>    - numberQ, total number of points on the mesh
   !>    - lambdaJob(:), miuJob(:), qValueJob(:), vectors storing
   !>      the lambda, miu and q(lambda,miu) values at each point
   !>      of the mesh
   !>    - numberQ_orig, numberConstraints_orig store the sizes of
   !>      the mesh (total size and number of constraints, resp.)
   !>      for use when modeMPI_deformation = 2, 4
   !---------------------------------------------------------------------
   Subroutine define_gridQlm()

     Integer :: iq,nQtotal,nQfactor,ii,lambda,miu,nQ,mpi_rank,mpi_err
     Integer, Allocatable :: nQmodul(:)
     Real, Allocatable :: deltaQ(:)

     Call mpi_comm_rank(MPI_COMM_WORLD,mpi_rank,mpi_err)

     If(modeMPI_deformation.Le.2) Then

        ! Determine total number of points on the mesh: number of collective
        ! variables 'numberConstraints' times the number of points for each
        numberQ = 1
        Do iq=1,numberConstraints
           numberQ = numberQ*qStepNumber(iq)
        End Do
        If(debug.Ge.1 .And. mpi_rank.Eq.0) &
           Write(6,'("In define_gridQlm() - numberQ = ",i10," numberConstraints=",i10)') numberQ,numberConstraints

        ! Storing characteristics of the default mesh for future use (if
        ! modeMPI_deformation = 2)
        numberQ_orig = numberQ
        numberConstraints_orig = numberConstraints

        ! Getting the meshsize for each collective variable
        Allocate(deltaQ(1:numberConstraints),nQmodul(0:numberConstraints))
        Do iq=1,numberConstraints
           If(qstepNumber(iq) > 1 ) Then
               deltaQ(iq) = (qFinValue(iq) - qinitValue(iq)) / Real(qstepNumber(iq) - 1)
           Else
               deltaQ(iq) = 0.0d0
           End If
        End Do

        ! Defining the mesh: lambdaJob(ii,iq), miuJob(ii,iq) and qValueJob(ii,iq) give,
        ! respectively, the value of lambda, miu and q(lambda,miu) at point ii of the
        ! total mesh and for the collective variable iq
        Allocate(qValueJob(1:numberQ,1:numberConstraints))
        Allocate(lambdaJob(1:numberQ,1:numberConstraints),miuJob(1:numberQ,1:numberConstraints))

        nQtotal=1; nQmodul(0)=1; nQfactor=0
        Do iq=1,numberConstraints

           nQ=qStepNumber(iq)
           lambda=Abs(qLambda(iq)); miu=qMiu(iq)
           nQtotal=nQtotal*nQ; nQmodul(iq)=numberQ/nQtotal

           Do ii=1,numberQ
              lambdaJob(ii,iq)=lambda
              miuJob(ii,iq)=miu
              qValueJob(ii,iq)=qinitValue(iq) &
                              + deltaQ(iq) * ((ii-1)/nQmodul(iq) - nQfactor*((ii-1)/nQmodul(iq-1))*nQ)
           End Do

           nQfactor = 1

        End Do

        Deallocate(deltaQ,nQmodul)

     Else

        Call read_path(fichier_path, .False.)
        ! Storing characteristics of the default mesh for future use (if
        ! modeMPI_deformation = 4)
        numberQ_orig = numberQ
        numberConstraints_orig = numberConstraints

     End If

   End Subroutine define_gridQlm

   !---------------------------------------------------------------------
   !> Defines restart deformation grid for modeMPI_deformation = 2, 4
   !>
   !> Outputs:
   !>   - numberQ_restart, total number of points on the restart
   !>     mesh
   !>   - lambdaJob_restart(:), miuJob_restart(:), and
   !>     qValueJob_restart(:), vectors storing the lambda, miu
   !>     and q(lambda,miu) values at each point of the restart
   !>     mesh
   !---------------------------------------------------------------------
   Subroutine define_gridQlm_new()

     Integer :: iq,nQtotal,nQfactor,ii,lambda,miu,nQ,mpi_rank,mpi_err
     Integer, Allocatable :: nQmodul(:)
     Real, Allocatable :: deltaQ(:)

     Call mpi_comm_rank(MPI_COMM_WORLD,mpi_rank,mpi_err)

     If(modeMPI_deformation.Eq.2) Then

        ! Determine total number of points on the restart mesh
        numberQ_restart = 1
        Do iq=1,numberConstraints_restart
           numberQ_restart = numberQ_restart * qStepNumber_restart(iq)
        End Do
        If(debug.Ge.1 .And. mpi_rank.Eq.0) &
           Write(6,'("In define_gridQlm_new() - numberQ_restart = ",i10," numberConstraints_restart=",i10)') &
                      numberQ_restart,numberConstraints_restart

        ! Getting the meshsize for each collective variable on the restart mesh
        Allocate(deltaQ(1:numberConstraints_restart),nQmodul(0:numberConstraints_restart))
        Do iq=1,numberConstraints_restart
           If(qstepNumber_restart(iq)>1) Then
              deltaQ(iq) = (qFinValue_restart(iq) - qinitValue_restart(iq)) / Real(qstepNumber_restart(iq) - 1)
           Else
              deltaQ(iq) = 0.0D0
           End If
        End Do

        ! Defining the restart mesh, see subroutine define_gridQlm for explanations
        Allocate(qValueJob_restart(1:numberQ_restart,1:numberConstraints_restart))
        Allocate(lambdaJob_restart(1:numberQ_restart,1:numberConstraints_restart),&
                 miuJob_restart(1:numberQ_restart,1:numberConstraints_restart))

        nQtotal=1; nQmodul(0)=1; nQfactor=0
        Do iq=1,numberConstraints_restart

           nQ=qStepNumber_restart(iq)
           lambda=Abs(qLambda_restart(iq)); miu=qMiu_restart(iq)
           nQtotal=nQtotal*nQ; nQmodul(iq)=numberQ_restart / nQtotal

           Do ii=1,numberQ_restart
              lambdaJob_restart(ii,iq)=lambda
              miuJob_restart(ii,iq)=miu
              qValueJob_restart(ii,iq)=qinitValue_restart(iq) &
                                      + deltaQ(iq) * ( (ii-1)/nQmodul(iq) - nQfactor*((ii-1)/nQmodul(iq-1))*nQ )
           End Do

           nQfactor = 1

        End Do

        Deallocate(deltaQ,nQmodul)

     Else

        Call read_path(fichier_path_new, .True.)

     End If

   End Subroutine define_gridQlm_new

   !---------------------------------------------------------------------
   !>  Subroutine mapping the default deformation grid to the restart
   !>  deformation grid.
   !>
   !>  Outputs:
   !>    - ii=mapping_qlm(jj) gives as before the index ii of
   !>      the original grid used to restart the calculation on
   !>      point jj of the restart grid
   !>    - numberQ, numberConstraints are overwritten and contain
   !>      the characteristics of the restart collective space
   !>    - qValueJob(:), lambdaJob(:) and miuJob(:) are similarly
   !>      associated with the restart grid
   !---------------------------------------------------------------------
   Subroutine mapping_gridQlm(IF_THO,LIPKIN,LIPKIP,ILCONT,IFCONT)

     Integer, INTENT(IN) :: IF_THO,LIPKIN,LIPKIP,ILCONT,IFCONT
     Integer :: ii,jj,index_optimal,iq,jq,lambda_restart,miu_restart, &
                lambda,miu,mpi_rank,mpi_err,mpi_world,exitstat,&
                pass_number
     Integer, Allocatable :: index_opt_vec(:)
     Real :: distance_min,distance,qlm_restart,qlm
     Real, Allocatable :: qmin(:),qmax(:),delta_q(:)
     Character(Len=132) :: directory_name
     Character(Len=68) :: fichier_courant,fichier_optimal
     Character(Len=150) :: commande

#if(USE_MANYCORES==1)
     Integer :: tribeID
     Integer :: numberMasters, numberHFODDproc, color,       &
                tribeRank, slaveRank, masterRank, worldRank, &
                tribeSize, slaveSize, masterSize, worldSize
     Integer :: worldGroup, groupMasters, mastersCOMM, &
                groupSlaves, slavesCOMM, tribeCOMM
     COMMON                                                 &
            /MPICOM/ worldGroup, groupMasters, groupSlaves, &
                     mastersCOMM, slavesCOMM, tribeCOMM
     COMMON                                                       &
            /MPIPRO/ numberMasters, numberHFODDproc, color,       &
                     tribeRank, slaveRank, masterRank, worldRank, &
                     tribeSize, slaveSize, masterSize, worldSize
#endif

     Call mpi_comm_rank(MPI_COMM_WORLD,mpi_rank,mpi_err)
     Call mpi_comm_size(MPI_COMM_WORLD,mpi_world,mpi_err)

     ! Default mapping: all new calculations restart from HFODD_00000001.REC
     Allocate(mapping_qlm(1:numberQ_restart))
     Do jj=1,numberQ_restart
        mapping_qlm(jj)=1
     End Do

     ! Move all restart files into current directory
     If(mpi_rank.Eq.0) Then
        commande = 'mv restart restart_old'; exitstat=0
        Call system(commande)
        If(debug.Ge.2) Write(6,'("mpi_rank = ",i10," executed command: ",a)') &
                                  mpi_rank,commande
        commande = 'mkdir restart'; exitstat=0
        Call system(commande)
        If(debug.Ge.2) Write(6,'("mpi_rank = ",i10," executed command: ",a)') &
                                  mpi_rank,commande
        ! Densities (for Lipkin method)
        If((LIPKIN.Ge.1.Or.LIPKIP.Ge.1).And.ILCONT.Eq.1) Then
           commande = 'mv lic lic_old'; exitstat=0
           Call system(commande)
           If(debug.Ge.2) Write(6,'("mpi_rank = ",i10," executed command: ",a)') &
                                     mpi_rank,commande
           commande = 'mkdir lic'; exitstat=0
           Call system(commande)
           If(debug.Ge.2) Write(6,'("mpi_rank = ",i10," executed command: ",a)') &
                                     mpi_rank,commande
        End If
        ! Fields (for Gogny)
        If(IFCONT.Eq.1) Then
           commande = 'mv fic fic_old'; exitstat=0
           Call system(commande)
           If(debug.Ge.2) Write(6,'("mpi_rank = ",i10," executed command: ",a)') &
                                     mpi_rank,commande
           commande = 'mkdir fic'; exitstat=0
           Call system(commande)
           If(debug.Ge.2) Write(6,'("mpi_rank = ",i10," executed command: ",a)') &
                                     mpi_rank,commande
        End If
        ! In batch mode, restart files will later be overwritten with record files. We make an
        ! additional copy here
        If(batch_mode.Eq.1) Then
           commande = 'mkdir restart_tmp'; exitstat=0
           Call system(commande)
           If(debug.Ge.2) Write(6,'("mpi_rank = ",i10," executed command: ",a)') &
                                                  mpi_rank,commande
           If((LIPKIN.Ge.1.Or.LIPKIP.Ge.1).And.ILCONT.Eq.1) Then
              commande = 'mkdir lic_tmp'; exitstat=0
              Call system(commande)
              If(debug.Ge.2) Write(6,'("mpi_rank = ",i10," executed command: ",a)') &
                                                     mpi_rank,commande
           End If
           If(IFCONT.Eq.1) Then
              commande = 'mkdir fic_tmp'; exitstat=0
              Call system(commande)
              If(debug.Ge.2) Write(6,'("mpi_rank = ",i10," executed command: ",a)') &
                                                     mpi_rank,commande
           End If
        End If
     End If

     ! Wait here for all processes
     Call mpi_barrier(MPI_COMM_WORLD,mpi_err)

     ! Allocate filenames
     Allocate(fichier_rec_old(1:numberQ_restart),fichier_rec_new(1:numberQ_restart))
     If(IF_THO.Ge.1) Allocate(fichier_tho_old(1:numberQ_restart),fichier_tho_new(1:numberQ_restart))
     If((LIPKIN.Ge.1.Or.LIPKIP.Ge.1).And.ILCONT.Eq.1) Allocate(fichier_lic_old(1:numberQ_restart),fichier_lic_new(1:numberQ_restart))
     If(IFCONT.Eq.1) Allocate(fichier_fic_old(1:numberQ_restart),fichier_fic_new(1:numberQ_restart))
     If(batch_mode.Eq.1) Then
        Allocate(fichier_rec_tmp(1:numberQ_restart))
        If(IF_THO.Ge.1) Allocate(fichier_tho_tmp(1:numberQ_restart))
        If((LIPKIN.Ge.1.Or.LIPKIP.Ge.1).And.ILCONT.Eq.1) Allocate(fichier_lic_tmp(1:numberQ_restart))
        If(IFCONT.Eq.1) Allocate(fichier_fic_tmp(1:numberQ_restart))
     End If

     ! Determine the metric for the new mesh
     Allocate(qmin(1:numberConstraints_restart),qmax(1:numberConstraints_restart),delta_q(1:numberConstraints_restart))
     qmin(1:numberConstraints_restart)=+1.0e9
     qmax(1:numberConstraints_restart)=-1.0e9
     Do jq=1,numberConstraints_restart
        Do jj=1,numberQ_restart
           If(qValueJob_restart(jj,jq).Le.qmin(jq)) qmin(jq)=qValueJob_restart(jj,jq)
           If(qValueJob_restart(jj,jq).Ge.qmax(jq)) qmax(jq)=qValueJob_restart(jj,jq)
        End Do
        If(Abs(qmax(jq)-qmin(jq)).Le.1.e-12) Then
           delta_q(jq) = 1.0
        Else
           delta_q(jq) = qmax(jq)-qmin(jq)
        End If
     End Do
     If(debug.Ge.1 .And. mpi_rank.Eq.0) Then
        Write(6,'("In mapping_gridQlm_new() - New grid")')
        Do jq=1,numberConstraints_restart
           Write(6,'("jq=",i2," qmin=",f20.10," qmax=",f20.10)') jq,qmin(jq),qmax(jq)
        End Do
     End If

     ! Determine index of the optimal restart file by using the normalized metric
     Allocate(index_opt_vec(1:numberQ_restart))
!$OMP  PARALLEL DO &
!$OMP  DEFAULT(NONE) &
!$OMP  SCHEDULE(STATIC) &
!$OMP  SHARED(numberQ_restart,numberQ,numberConstraints_restart,numberConstraints,lambdaJob_restart, &
!$OMP         miuJob_restart,qValueJob_restart,lambdaJob,miuJob,qValueJob,index_opt_vec,delta_q,qmin) &
!$OMP  PRIVATE(jj,index_optimal,distance_min,ii,distance,jq,iq,lambda_restart,miu_restart, &
!$OMP          qlm_restart,lambda,miu,qlm)
     ! TODO: ADD A METRIC TO COMPUTE THE DISTANCES
     Do jj=1,numberQ_restart
        index_optimal=1; distance_min=1.D9
        Do ii=1,numberQ
           distance = 0.0D0
           Do jq=1,numberConstraints_restart
              Do iq=1,numberConstraints
                 ! new constraints
                 lambda_restart=lambdaJob_restart(jj,jq)
                 miu_restart = miuJob_restart(jj,jq)
                 qlm_restart = (qValueJob_restart(jj,jq)-qmin(jq))/delta_q(jq) ! 0 <= qlm_rest. <= 1
                 ! old constraints
                 lambda=lambdaJob(ii,iq)
                 miu=miuJob(ii,iq)
                 qlm=(qValueJob(ii,iq)-qmin(iq))/delta_q(iq)
                 If(lambda.Eq.lambda_restart .And. miu.Eq.miu_restart) Then
                    distance = distance + (qlm - qlm_restart)**2
                 End If
              End Do
           End Do
           If(distance.Le.distance_min) Then
              distance_min=distance; index_optimal = ii
           End If
        End Do
        index_opt_vec(jj) = index_optimal
     End Do
!$OMP END PARALLEL DO
     Deallocate(qmin,qmax,delta_q)

     ! Create a list of optimal restart files, and define the mapping
     ! NB: In multicore mode, filenames are *different* for each process within the tribe, since they contain
     !     an additional label related to their tribe ID.
     Do jj=1,numberQ_restart
        index_optimal=index_opt_vec(jj)
        ! move proper file into restart directory with proper name
        Write(fichier_optimal,'("./restart_old/HFODD_",I8.8,".REC")') index_optimal
        Write(fichier_courant,'("./restart/HFODD_",I8.8,".REC")') jj
#if(USE_MANYCORES==1)
        Write(fichier_optimal,'("./restart_old/proc_",i3.3,"_HFODD_",I8.8,".REC")') tribeRank,index_optimal
        Write(fichier_courant,'("./restart/proc_",i3.3,"_HFODD_",I8.8,".REC")') tribeRank,jj
#endif
        fichier_rec_old(jj)=Trim(fichier_optimal)
        fichier_rec_new(jj)=Trim(fichier_courant)
        If(IF_THO.Ge.1) Then
           Write(fichier_optimal,'("./tho/t",I6.6,".hel")') index_optimal
           Write(fichier_courant,'("./t",I6.6,".hel")') jj
           fichier_tho_old(jj)=Trim(fichier_optimal)
           fichier_tho_new(jj)=Trim(fichier_courant)
        End If
        ! Densities
        If((LIPKIN.Ge.1.Or.LIPKIP.Ge.1).And.ILCONT.Eq.1) Then
           Write(fichier_optimal,'("./lic_old/HFODD_",I8.8,".LIC")') index_optimal
           Write(fichier_courant,'("./lic/HFODD_",I8.8,".LIC")') jj
#if(USE_MANYCORES==1)
           Write(fichier_optimal,'("./lic_old/proc_",i3.3,"_HFODD_",I8.8,".LIC")') tribeRank,index_optimal
           Write(fichier_courant,'("./lic/proc_",i3.3,"_HFODD_",I8.8,".LIC")') tribeRank,jj
#endif
           fichier_lic_old(jj)=Trim(fichier_optimal)
           fichier_lic_new(jj)=Trim(fichier_courant)
        End If
        ! Fields
        If(IFCONT.Eq.1) Then
           Write(fichier_optimal,'("./fic_old/HFODD_",I8.8,".FIC")') index_optimal
           Write(fichier_courant,'("./fic/HFODD_",I8.8,".FIC")') jj
#if(USE_MANYCORES==1)
           Write(fichier_optimal,'("./fic_old/proc_",i3.3,"_HFODD_",I8.8,".FIC")') tribeRank,index_optimal
           Write(fichier_courant,'("./fic/proc_",i3.3,"_HFODD_",I8.8,".FIC")') tribeRank,jj
#endif
           fichier_fic_old(jj)=Trim(fichier_optimal)
           fichier_fic_new(jj)=Trim(fichier_courant)
        End If
        mapping_qlm(jj)=jj
        ! In batch mode, we need to make an extra copy of all these restart files, since they will
        ! be overwritten later
        If(batch_mode.Eq.1) Then
           Write(fichier_courant,'("./restart_tmp/HFODD_",I8.8,".REC")') jj
#if(USE_MANYCORES==1)
           Write(fichier_courant,'("./restart_tmp/proc_",i3.3,"_HFODD_",I8.8,".REC")') tribeRank,jj
#endif
           fichier_rec_tmp(jj)=Trim(fichier_courant)
           If(IF_THO.Eq.1) Then
              Write(fichier_courant,'("./t",I6.6,".hel")') jj
              fichier_tho_tmp(jj)=Trim(fichier_courant)
           End If
           If((LIPKIN.Ge.1.Or.LIPKIP.Ge.1).And.ILCONT.Eq.1) Then
              Write(fichier_courant,'("./lic_tmp/HFODD_",I8.8,".LIC")') jj
#if(USE_MANYCORES==1)
              Write(fichier_courant,'("./lic_tmp/proc_",i3.3,"_HFODD_",I8.8,".LIC")') tribeRank,jj
#endif
              fichier_lic_tmp(jj)=Trim(fichier_courant)
           End If
           If(IFCONT.Eq.1) Then
              Write(fichier_courant,'("./fic_tmp/HFODD_",I8.8,".FIC")') jj
#if(USE_MANYCORES==1)
              Write(fichier_courant,'("./fic_tmp/proc_",i3.3,"_HFODD_",I8.8,".FIC")') tribeRank,jj
#endif
              fichier_fic_tmp(jj)=Trim(fichier_courant)
           End If
        End If

        If(mpi_rank.Eq.0.And.debug.Ge.4) Then
           Write(6,'("mpi_rank = ",i10," mapping_qlm - jj = ",i10," index_optimal = ",i10)') &
                      mpi_rank,jj,index_optimal
           Do jq=1,numberConstraints_restart
              Write(6,'("RESTART  - mpi_rank = ",i10," l=",i2,"m=",i2," Qlm=",f10.5)') &
                                    mpi_rank,lambdaJob_restart(jj,jq),                 &
                                             miuJob_restart(jj,jq),                    &
                                             qValueJob_restart(jj,jq)
           End do
           Do iq=1,numberConstraints
              Write(6,'("ORIGINAL - mpi_rank = ",i10," l=",i2,"m=",i2," Qlm=",f10.5)') &
                                    mpi_rank,lambdaJob(index_optimal,iq),              &
                                              miuJob(index_optimal,iq),                &
                                             qValueJob(index_optimal,iq)
           End do
        End If

     End Do
     Deallocate(index_opt_vec)

     ! System copy of the relevant files into their new names
     ! NB: As mentioned above, in multicore mode, filenames are different for each tribe member. This
     !     is the reason for the loop over tribe member and the IF/THEN condition below
     Do jj=1,numberQ_restart
#if(USE_MANYCORES==1)
        Do tribeID=0,tribeSize
           If(Mod(tribeRank,tribeSize).Eq.Mod(tribeID,tribeSize).And.color.Eq.Mod(jj-1,numberMasters)) Then
#else
        If(Mod(mpi_rank,mpi_world).Eq.Mod(jj-1,mpi_world)) Then
#endif
           commande = 'cp '//fichier_rec_old(jj)//' '//fichier_rec_new(jj); exitstat=0
           Call system(commande)
           If(debug.Ge.2) Write(6,'("mpi_rank = ",i10," jj = ",i10," - Executed command: ",a)') &
                                     mpi_rank,jj,commande
           If(IF_THO.Ge.1) Then
              commande = 'cp '//fichier_tho_old(jj)//' '//fichier_tho_new(jj); exitstat=0
              Call system(commande)
              If(debug.Ge.2) Write(6,'("mpi_rank = ",i10," jj = ",i10," - Executed command: ",a)') &
                                        mpi_rank,jj,commande
           End If
           ! Densities
           If((LIPKIN.Ge.1.Or.LIPKIP.Ge.1).And.ILCONT.Eq.1) Then
              commande = 'cp '//fichier_lic_old(jj)//' '//fichier_lic_new(jj); exitstat=0
              Call system(commande)
              If(debug.Ge.2) Write(6,'("mpi_rank = ",i10," jj = ",i10," - Executed command: ",a)') &
                                        mpi_rank,jj,commande
           End If
           ! Fields
           If(IFCONT.Eq.1) Then
              commande = 'cp '//fichier_fic_old(jj)//' '//fichier_fic_new(jj); exitstat=0
              Call system(commande)
              If(debug.Ge.2) Write(6,'("mpi_rank = ",i10," jj = ",i10," - Executed command: ",a)') &
                                        mpi_rank,jj,commande
           End If
           ! Make extra copy in batch mode
           If(batch_mode.Eq.1) Then
              commande = 'cp '//fichier_rec_old(jj)//' '//fichier_rec_tmp(jj); exitstat=0
              Call system(commande)
              If(debug.Ge.2) Write(6,'("mpi_rank = ",i10," jj = ",i10," - Executed command: ",a)') &
                                                    mpi_rank,jj,commande
              If(IF_THO.Ge.1) Then
                 commande = 'cp '//fichier_tho_old(jj)//' '//fichier_tho_tmp(jj); exitstat=0
                 Call system(commande)
                 If(debug.Ge.2) Write(6,'("mpi_rank = ",i10," jj = ",i10," - Executed command: ",a)') &
                                                       mpi_rank,jj,commande
              End If
              If((LIPKIN.Ge.1.Or.LIPKIP.Ge.1).And.ILCONT.Eq.1) Then
                 commande = 'cp '//fichier_lic_old(jj)//' '//fichier_lic_tmp(jj); exitstat=0
                 Call system(commande)
                 If(debug.Ge.2) Write(6,'("mpi_rank = ",i10," jj = ",i10," - Executed command: ",a)') &
                                                        mpi_rank,jj,commande
              End If
              If(IFCONT.Eq.1) Then
                 commande = 'cp '//fichier_fic_old(jj)//' '//fichier_fic_tmp(jj); exitstat=0
                 Call system(commande)
                 If(debug.Ge.2) Write(6,'("mpi_rank = ",i10," jj = ",i10," - Executed command: ",a)') &
                                                        mpi_rank,jj,commande
              End If
           End If
#if(USE_MANYCORES==1)
           End If ! tribeRank vs. tribe ID
        End Do ! loop over tribe ID
#else
        End If ! mpi_rank versus deformation number
#endif
     End Do ! loop over deformation number

     ! Wait here for all processes
     Call mpi_barrier(MPI_COMM_WORLD,mpi_err)

     ! Check that all files in restart/ have truly been copied
     Allocate(table_fichiers(1:numberQ_restart)); table_fichiers(:)=0
     directory_name='restart'

     pass_number=1; filled=.False.
     Do While (.Not.filled)

     If(mpi_rank.Eq.0.And.debug.Ge.2) Write(6,'("Pass =",i5)') pass_number
     Call check_directory(directory_name,numberQ_restart)

        Call fill_holes(numberQ_restart,IF_THO,LIPKIN,LIPKIP,ILCONT,IFCONT)

        pass_number=pass_number+1
        If(pass_number.Gt.4) Exit ! avoid infinite loops

     End Do
     Deallocate(table_fichiers)

     ! Check that all files in restart_tmp/ have truly been copied (batch mode)
     If(batch_mode.Eq.1) Then

        Allocate(table_fichiers(1:numberQ_restart)); table_fichiers(:)=0
        directory_name='restart_tmp'

        pass_number=1; filled=.False.
        Do While (.Not.filled)

           If(mpi_rank.Eq.0.And.debug.Ge.2) Write(6,'("Pass =",i5)') pass_number
           Call check_directory(directory_name,numberQ_restart)

           Call fill_holes(numberQ_restart,IF_THO,LIPKIN,LIPKIP,ILCONT,IFCONT)

           pass_number=pass_number+1
           If(pass_number.Gt.4) Exit ! avoid infinite loops

        End Do
        Deallocate(table_fichiers)
     End If

     Deallocate(fichier_rec_old,fichier_rec_new)
     If(IF_THO.Ge.1) Deallocate(fichier_tho_old,fichier_tho_new)
     If((LIPKIN.Ge.1.Or.LIPKIP.Ge.1).And.ILCONT.Eq.1) Deallocate(fichier_lic_old,fichier_lic_new)
     If(IFCONT.Eq.1) Deallocate(fichier_fic_old,fichier_fic_new)
     If(batch_mode.Eq.1) Then
        Deallocate(fichier_rec_tmp)
        If(IF_THO.Ge.1) Deallocate(fichier_tho_tmp)
        If((LIPKIN.Ge.1.Or.LIPKIP.Ge.1).And.ILCONT.Eq.1) Deallocate(fichier_lic_tmp)
        If(IFCONT.Eq.1) Deallocate(fichier_fic_tmp)
     End If

     Deallocate(qValueJob,lambdaJob,miuJob)

     numberQ=numberQ_restart
     numberConstraints=numberConstraints_restart
     If(mpi_rank.Eq.0.And.debug.Ge.2) &
        Write(6,'("In mapping_gridQlm_new - numberQ =",i10," numberConstraints = ",i5)') numberQ, numberConstraints

     Allocate(qValueJob(1:numberQ,1:numberConstraints))
     Allocate(lambdaJob(1:numberQ,1:numberConstraints),miuJob(1:numberQ,1:numberConstraints))
     qValueJob(:,:)=qValueJob_restart(:,:)
     lambdaJob(:,:)=lambdaJob_restart(:,:); miuJob(:,:)=miuJob_restart(:,:)
     Deallocate(qValueJob_restart,lambdaJob_restart,miuJob_restart)

   End Subroutine mapping_gridQlm

   !---------------------------------------------------------------------
   !>  Subroutine defining the tasks
   !---------------------------------------------------------------------
   Subroutine mpi_constructJobList(NUMITE,NXHERM,NYHERM,NZHERM,IFCONT, &
                                   IPCONT,ILCONT,IACONT,IMCONT,IRCONT, &
                                                        IRENMA,IRENIN, &
                                                        LIPKIN,LIPKIP, &
                                                 REFERN,REFERP,REDELN, &
                                                 REDELP,REFE2N,REFE2P, &
                                          IFIBLN,INIBLN,IFIBLP,INIBLP, &
                                          NMUMAX,NSIMAX,NMUCON,ISHIFT, &
                                                               IF_THO, &
                                                               IERROR)

     Integer, INTENT(IN) :: NUMITE,NXHERM,NYHERM,NZHERM,IPCONT,ILCONT,IACONT, &
                            IMCONT,IRCONT,IRENMA,IRENIN,LIPKIN,LIPKIP,IFIBLN, &
                            INIBLN,IFIBLP,INIBLP,NMUMAX,NSIMAX,NMUCON,ISHIFT, &
                            IF_THO,IFCONT
     Integer, INTENT(INOUT) :: IERROR
     Real, INTENT(INOUT) :: REFERN,REFERP,REDELN,REDELP,REFE2N,REFE2P

     Integer :: mpi_rank, mpi_err
     Integer :: iz,in,iq,nCommon,iCommon,force,indexJob,zNumber,nNumber,ii
     Integer :: i_shel,i_stat,i_freq,i_defo
     Integer, Allocatable :: zVector(:),nVector(:)

     Real :: xzNumber, xnNumber
     Real, Allocatable :: xzVector(:),xnVector(:)

     Call mpi_comm_rank(MPI_COMM_WORLD,mpi_rank,mpi_err)

     If(debug.Ge.1 .And. mpi_rank .Eq.0) &
        Write(6,'("modeMPI_deformation=",i3," modeMPI_basis=",i3)') &
                   modeMPI_deformation,modeMPI_basis

     ! Define deformation grid
     Call define_gridQlm()

     ! Define restart grid and mapping indexes
     If(modeMPI_deformation.Eq.2.Or.modeMPI_deformation.Eq.4) Then
        Call define_gridQlm_new()
        Call mapping_gridQlm(IF_THO,LIPKIN,LIPKIP,ILCONT,IFCONT)
     Else
        ! default allocations of mapping arrays: if no restart is required, these mappings simply
        ! provide identity vectors (v(i) = i)
        Allocate(mapping_qlm(1:numberQ))
        Do ii=1,numberQ
           mapping_qlm(ii)=ii
           If(debug.Ge.4 .And. mpi_rank .Eq.0) Write(6,'("ii = ",i10," mapping = ",i10)') ii,mapping_qlm(ii)
        End Do
     End If

     ! Define nucleus grid: integers
     Allocate(zVector(1:proton_stepNumber),nVector(1:neutron_stepNumber))
     Do iz=1,proton_stepNumber
        zVector(iz) = proton_initValue + (iz-1)*proton_stepValue
     End Do
     Do in=1, neutron_stepNumber
        nVector(in) = neutron_initValue + (in-1)*neutron_stepValue
     End Do

     ! Define nucleus grid: float
     If(x_proton_initValue.Gt.0.0D0) Then
        Allocate(xzVector(1:proton_stepNumber))
        Do iz=1,proton_stepNumber
           xzVector(iz) = x_proton_initValue + (iz-1)*x_proton_stepValue
        End Do
     End If
     If(x_neutron_initValue.Gt.0.0D0) Then
        Allocate(xnVector(1:neutron_stepNumber))
        Do in=1, neutron_stepNumber
           xnVector(in) = x_neutron_initValue + (in-1)*x_neutron_stepValue
        End Do
     End If

     ! Default deformation mode
     NOFJOB = proton_stepNumber * neutron_stepNumber * counterForce * numberQ

     Allocate(IZ_JOB(1:NOFJOB),IN_JOB(1:NOFJOB),indexQjob(1:NOFJOB),indexRestartJob(1:NOFJOB),indexForce(1:NOFJOB))
     Allocate(XIZJOB(1:NOFJOB),XINJOB(1:NOFJOB)); XIZJOB(:)=-1.0D0; XINJOB(:)=-1.0D0

     nCommon = proton_stepNumber * neutron_stepNumber * counterForce

     If(debug.Ge.1 .And. mpi_rank .Eq.0) Then
        Write(6,'("Sanity checks:")')
        Write(6,'("proton_stepNumber=",i10," neutron_stepNumber=",i10," numberConstraints=",i10)') &
                   proton_stepNumber,neutron_stepNumber,numberConstraints
        Write(6,'("counterForce=",i10," numberQ=",i10)') &
                   counterForce, numberQ
     End If

     iCommon=0; indexJob=0

     Do iz=1,proton_stepNumber
        zNumber = zVector(iz)
        If(x_proton_initValue.Gt.0.0D0) xzNumber=xzVector(iz)
        Do in=1,neutron_stepNumber
           nNumber = nVector(in)
           If(x_neutron_initValue.Gt.0.0D0) xnNumber=xnVector(in)
           Do force=1,counterForce

              iCommon = iCommon + 1

              Do ii=1,numberQ

                 indexJob=indexJob+1
                 indexQjob(indexJob)=ii

                 indexForce(indexJob)=force

                 IZ_JOB(indexJob)=zNumber
                 IN_JOB(indexJob)=nNumber

                 If(x_proton_initValue.Gt.0.0D0) Then
                    XIZJOB(indexJob)=xzNumber
                 Else
                    XIZJOB(indexJob)=Real(zNumber,Kind=Kind(1.0))
                 End If
                 If(x_neutron_initValue.Gt.0.0D0) Then
                    XINJOB(indexJob)=xnNumber
                 Else
                    XINJOB(indexJob)=Real(nNumber,Kind=Kind(1.0))
                 End If

                 indexRestartJob(indexJob)=indexJob

                  If(debug.Ge.4 .And. mpi_rank .Eq.0) Then
                     Do iq=1,numberConstraints
                        Write(6,'("Z=",i3," N=",i3," force=",a4,              &
                     &            " iq=",i2," l=",i2," m=",i2," Qlm=",f10.5)')&
                                zNumber,nNumber,                              &
                                forceVector(indexForce(indexJob)),            &
                                iq,lambdaJob(indexQjob(indexJob),iq),         &
                                miuJob(indexQjob(indexJob),iq),               &
                                qValueJob(indexQjob(indexJob),iq)
                     End Do
                  End If

              End Do

           End Do
        End Do
     End Do

     Deallocate(mapping_qlm)

     If(modeMPI_basis==1) Then

        Allocate(basis_size(1:basis_size_stepNumber))
        Allocate(basis_states(1:basis_states_stepNumber))
        Allocate(basis_freq(1:basis_frequency_stepNumber))
        Allocate(basis_defo(1:basis_deform_stepNumber))

        Do i_shel=1,basis_size_stepNumber
           basis_size(i_shel) = basis_size_initValue + (i_shel-1)*basis_size_stepValue
        End do
        Do i_stat=1,basis_states_stepNumber
           basis_states(i_stat) = basis_states_initValue + (i_stat-1)*basis_states_stepValue
        End do
        do i_freq=1,basis_frequency_stepNumber
           basis_freq(i_freq) = basis_frequency_initValue + Real(i_freq-1)*basis_frequency_stepValue
        End do
        do i_defo=1,basis_deform_stepNumber
           basis_defo(i_defo) = basis_deform_initValue + Real(i_defo-1)*basis_deform_stepValue
        End do

        NOFJOB = proton_stepNumber * neutron_stepNumber * numberQ   &
               * basis_frequency_stepNumber * basis_deform_stepNumber  &
               * basis_size_stepNumber * basis_states_stepNumber * counterForce

        Allocate(IZ_JOB(1:NOFJOB),IN_JOB(1:NOFJOB),indexQjob(1:NOFJOB),indexForce(1:NOFJOB),&
                 indexSizeJob(1:NOFJOB),indexStatesJob(1:NOFJOB),indexLengthJob(1:NOFJOB),  &
                 idxDeformJob(1:NOFJOB))
        Allocate(XIZJOB(1:NOFJOB),XINJOB(1:NOFJOB)); XIZJOB(:)=-1.0D0; XINJOB(:)=-1.0D0

        indexJob = 0

        Do iz=1,proton_stepNumber
           zNumber = zVector(iz)
           If(x_proton_initValue.Gt.0.0D0) xzNumber=xzVector(iz)
           Do in=1,neutron_stepNumber
              nNumber = nVector(in)
              If(x_neutron_initValue.Gt.0.0D0) xnNumber=xnVector(in)
              Do force=1,counterForce
                  Do ii=1,numberQ
                     Do i_shel=1,basis_size_stepNumber
                        Do i_stat=1,basis_states_stepNumber
                           Do i_freq=1,basis_frequency_stepNumber
                              Do i_defo=1,basis_deform_stepNumber

                                 indexJob = indexJob + 1

                                 IZ_JOB(indexJob)=zNumber
                                 IN_JOB(indexJob)=nNumber

                                 If(x_proton_initValue.Gt.0.0D0) Then
                                    XIZJOB(indexJob)=xzNumber
                                 Else
                                    XIZJOB(indexJob)=Real(zNumber,Kind=Kind(1.0))
                                 End If
                                 If(x_neutron_initValue.Gt.0.0D0) Then
                                    XINJOB(indexJob)=xnNumber
                                 Else
                                    XINJOB(indexJob)=Real(nNumber,Kind=Kind(1.0))
                                 End If

                                 indexQjob(indexJob)=ii
                                 indexForce(indexJob)=force

                                 indexSizeJob(indexJob)=i_shel
                                 indexStatesJob(indexJob)=i_stat
                                 indexLengthJob(indexJob)=i_freq
                                 idxDeformJob(indexJob)=i_defo

                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
              End Do
           End Do
        End Do

     End If

  End Subroutine mpi_constructJobList

  !---------------------------------------------------------------------
  !>  Subroutine setting local (process-wise) HFODD variables based on
  !>  existing list of tasks. Routine executed by each process.
  !---------------------------------------------------------------------
  Subroutine mpi_getCurrentData(INDJOB,IN_FIX,IZ_FIX,SKYRME,GOGNAM, &
                                              INSIQN,IPSIQN,IDSIQN, &
                                              INSIQP,IPSIQP,IDSIQP, &
                                       INSIGN,IPSIGN,ISSIGN,IDSIGN, &
                                       INSIGP,IPSIGP,ISSIGP,IDSIGP, &
                                       ICONTI,IPCONT,ILCONT,IFCONT, &
                                                            I_GOGA, &
                                       SLOWEV,SLOWOD,SLOWPA,SLOWLI, &
                                                            IROTAT, &
                                                            NILXYZ, &
                                                            N_EXIT, &
                                                            NOITER, &
                                              FCHOM0,NOSCIL,NLIMIT, &
                                                     LIPKIN,LIPKIP, &
                                       FE2FIN,IF2FIN,FE2FIP,IF2FIP, &
                                       NXHERM,NYHERM,NZHERM,IOPTGS, &
                                              IPAIRI,IPAHFB,IMFHFB, &
                                                            IWRILI, &
                                                            EPSITE, &
                                                            INIROT, &
                                                            IGAMMA, &
                                                            LANODD, &
                                       IFIBLN,INIBLN,IFIBLP,INIBLP, &
                                                            NMUCON, &
                                                            IFSHEL, &
                                                            IF_RPA, &
                                                            INDFIL, &
                                              L_COLL,M_COLL,N_COLL)

     Use hfodd_sizes

     Integer :: INDJOB,IN_FIX,IZ_FIX,INSIQN,IPSIQN,IDSIQN,INSIQP,IPSIQP,IDSIQP,&
                INSIGN,IPSIGN,ISSIGN,IDSIGN,INSIGP,IPSIGP,ISSIGP,IDSIGP,ICONTI,&
                IPCONT,ILCONT,IROTAT,NILXYZ,N_EXIT,NOITER,NOSCIL,NLIMIT,LIPKIN,&
                LIPKIP,NXHERM,NYHERM,NZHERM,IOPTGS,IPAIRI,IPAHFB,IMFHFB,IWRILI,&
                IFCONT,INIROT,IGAMMA,LANODD,IFIBLN,INIBLN,IFIBLP,INIBLP,NMUCON,&
                IFSHEL,IF_RPA,INDFIL,I_GOGA,IF2FIN,IF2FIP,N_COLL

     Integer :: IFLAGQ,IFLALQ,indexQ,ii,LAMBDA,MIU,IFNECK,IFASYM,IFDIST
     Integer :: i_shel,i_stat,i_freq,i_defo,force,is_deformed
     Integer :: mpi_rank,mpi_err

     !Integer, Dimension(:), POINTER :: L_COLL, M_COLL
     Integer, Allocatable :: L_COLL(:), M_COLL(:)

     Character(Len=4) :: SKYRME,GOGNAM

     Real :: SLOWEV,SLOWOD,SLOWPA,SLOWLI,EPSITE,FCHOM0,STIFFQ, &
             QASKED,GALMUQ,QLINEA,Q0NECK,G_NECK,HBMASS,HBMRPA, &
             HBMINP,H_BARC,XMASSP,XMASSN,HBCOE2,HBAROX,HBAROY, &
             HBAROZ,HOMSCA,ALPHAR,PIARGU,OVALUE,A_MASS,G_ASYM, &
             HOMEGA,XINFIX,XIZFIX,FE2FIN,FE2FIP,G_DIST,XIASYM, &
             D_COMS

     COMMON                                           &
            /FLOPAR/ XIZFIX,XINFIX
     COMMON                                           &
            /NCKVAL/ Q0NECK,G_NECK                    &
            /MASVAL/ XIASYM,G_ASYM                    &
            /DISVAL/ D_COMS,G_DIST
     COMMON                                           &
            /NCKFLA/ IFNECK                           &
            /MASFLA/ IFASYM                           &
            /DISFLA/ IFDIST
     COMMON                                           &
            /REALPH/ ALPHAR(0:NDLAMB,0:NDLAMB)
     COMMON                                           &
            /QCNSTR/ STIFFQ(0:NDMULT,-NDMULT:NDMULT), &
                     QASKED(0:NDMULT,-NDMULT:NDMULT), &
                     IFLAGQ(0:NDMULT,-NDMULT:NDMULT)
     COMMON                                           &
            /QLASTR/ GALMUQ(0:NDMULT,-NDMULT:NDMULT), &
                     QLINEA(0:NDMULT,-NDMULT:NDMULT), &
                     IFLALQ(0:NDMULT,-NDMULT:NDMULT)
     COMMON                                           &
            /PLANCK/ HBMASS,HBMRPA,HBMINP
     COMMON                                           &
            /PHYCON/ H_BARC,XMASSP,XMASSN,HBCOE2
     COMMON                                           &
            /BASPAR/ HBAROX,HBAROY,HBAROZ             &
            /SCALNG/ HOMSCA(NDKART)

     Call mpi_comm_rank(MPI_COMM_WORLD,mpi_rank,mpi_err)

     ! Set proton, neutron number and Skyrme functional
     IN_FIX=IN_JOB(INDJOB); IZ_FIX=IZ_JOB(INDJOB)
     XINFIX=XINJOB(INDJOB); XIZFIX=XIZJOB(INDJOB)
     If(debug.Ge.1) Write(6,'("In mpi_getCurrentData - mpi_rank = ",i10," XN = ",f20.14," XZ = ",f20.14)') &
                                                       mpi_rank, XINFIX, XIZFIX
     If(I_GOGA.Gt.0) Then
        force=indexForce(INDJOB)
        GOGNAM=forceVector(force)
        SKYRME=GOGNAM
     Else
        force=indexForce(INDJOB)
        SKYRME=forceVector(force)
     End If

     ! Set index of the restart file (assumed to be in restart/ directory,
     ! and to have the generic form HFODD_XXXXXXXX.REC, where XXXXXXXX is a 8 digit integer
     If(modeMPI_deformation.Ge.1) Then
        INDFIL=indexRestartJob(INDJOB)
        If(debug.Ge.1) Write(6,'("In mpi_getCurrentData - mpi_rank = ",i10," INDFIL = ",i10)') &
                                                          mpi_rank,INDFIL
     End If

     ! Set constraints on collective moments, and activate constraints
     indexQ=indexQjob(INDJOB)
     If(debug.Ge.1) Write(6,'("In mpi_getCurrentData - mpi_rank = ",i10," INDJOB = ",i10," indexQ = ",i10)') &
                                                       mpi_rank,INDJOB,indexQ
     Do ii=1,numberConstraints
        LAMBDA=lambdaJob(indexQ,ii)
        MIU = miuJob(indexQ, ii)
        If(LAMBDA.Ge.NMUCON) NMUCON=LAMBDA
        ! Regular multipole moments
        If(LAMBDA.GE.1.And.Miu.Ge.0) Then
           QASKED(LAMBDA,MIU) = qValueJob(indexQ, ii)
           STIFFQ(LAMBDA,MIU) = 0.5D0
           IFLAGQ(LAMBDA,MIU) = 0
           IFLALQ(LAMBDA,MIU) = 0
           If(turnConstraintOn.Eq.1) Then
              IFLAGQ(LAMBDA,MIU) = 1
                              IFLALQ(LAMBDA,MIU) = 1 ! ALM
              If(IF_RPA.Eq.1) IFLALQ(LAMBDA,MIU) =-1 ! RPA
           End If
        End If
        ! Constraint on the size of the neck
        If(LAMBDA.Eq.0.And.MIU.Eq.0.And.modeMPI_deformation.Ge.5) Then
           IFNECK=2
           Q0NECK=qValueJob(indexQ, ii)
        End If
        ! Constraint on the LN parameters for protons and neutrons
        If(LAMBDA.Eq.-1.And.MIU.Eq.-1.And.modeMPI_deformation.Ge.5) Then
           IF2FIP=1
           FE2FIP=qValueJob(indexQ, ii)
        End If
        If(LAMBDA.Eq.-1.And.MIU.Eq.1.And.modeMPI_deformation.Ge.5) Then
           IF2FIN=1
           FE2FIN=qValueJob(indexQ, ii)
        End If
        ! Constraint on distance between pre-fragments
        If(LAMBDA.Eq.-2.And.MIU.Eq.0.And.modeMPI_deformation.Ge.5) Then
           IFDIST=2
           D_COMS=qValueJob(indexQ, ii)
        End If
        ! Constraint on asymmetry between pre-fragment masses
        If(LAMBDA.Eq.-3.And.MIU.Eq.0.And.modeMPI_deformation.Ge.5) Then
           IFASYM=2
           XIASYM=qValueJob(indexQ, ii)
        End If
     End Do

     ! Set up the characteristics of collective inertia
     If(modeMPI_deformation.Ge.5) Then
        If(Allocated(L_COLL)) Deallocate(L_COLL,M_COLL)
        Allocate(L_COLL(1:NDMULT))
        Allocate(M_COLL(1:NDMULT))
        N_COLL=0
        Do ii=1,numberConstraints
           ! Do not include Q10 as a collective variable
           If(Abs(lambdaJob(indexQ,ii)).Ne.1.Or.miuJob(indexQ,ii).Ne.0) Then
              N_COLL=N_COLL+1
              L_COLL(N_COLL)=lambdaJob(indexQ,ii)
              M_COLL(N_COLL)=miuJob(indexQ,ii)
           End If
        End Do
     End If

     ! Set characteristics of the basis for this process
     If(modeMPI_basis.Eq.1) Then

        i_shel=indexSizeJob(INDJOB)
        i_stat=indexStatesJob(INDJOB)
        i_freq=indexLengthJob(INDJOB)
        i_defo=idxDeformJob(INDJOB)

        ! Number of shells and number of states
        NOSCIL = basis_size(i_shel)
        NLIMIT = basis_states(i_stat)

        ! Check if the requested configuration is deformed
        is_deformed = 0
        Do lambda=0,NDMULT
           Do miu=0,NDMULT
              If(QASKED(LAMBDA,MIU).Gt.1.D-12) Then
                 If(debug.Ge.4) Write(6,'("LAMBDA = ",i2," MIU = ",i2," QASKED = ",f20.14)') &
                                           LAMBDA,MIU,QASKED(LAMBDA,MIU)
                 is_deformed = 1
              End If
           End Do
        End Do
        ! Set frequency and deformation of the basis
        A_MASS=Real(IZ_FIX+IN_FIX); HOMEGA=41.0D0/A_MASS**(1.0D0/3.0D0)
        OVALUE=basis_freq(i_freq); FCHOM0=OVALUE/HOMEGA
        ALPHAR(2,0)=basis_defo(i_defo)
        If(debug.Ge.1) Write(6,'("NOSCIL=",i4," NLIMIT=",i4," FCHOM0=",f10.5," ALPHAR=",f10.5)') &
                                  NOSCIL,NLIMIT,FCHOM0,ALPHAR(2,0)

     End If

   End Subroutine mpi_getCurrentData

   !---------------------------------------------------------------------
   !>  Subroutine that allows to redefine the requested values of the
   !>  constraints based on the values of Qlm read from disk. This is
   !>  used in restart mode only to effectively constrain previously free
   !>  collective variables.
   !>
   !>  Example:
   !>    - Run 1: Q20 grid, Q40 unconstrained
   !>    - Run 2: Same Q20 grid, new Q40 grid between 0 and 50.
   !>             Calling mpi_pathExploration will recenter the Q40
   !>             grid around the equilibrium values of Q40 from
   !>             Run 1.
   !---------------------------------------------------------------------
   Subroutine mpi_pathExploration(IERROR,NMUCON)

     Use hfodd_sizes

     Integer :: IERROR,NMUCON,ITITLE,NFIPRI,IFLAGQ,LAMACT,MIUACT,i,lambda,miu

     Real :: QOLDIE,QSHIFT,STIFFQ,QASKED,QMUL_T,QMUL_P,QMUL_N

     COMMON                                           &
            /CFIPRI/ NFIPRI
     COMMON                                           &
            /QCNSTR/ STIFFQ(0:NDMULT,-NDMULT:NDMULT), &
                     QASKED(0:NDMULT,-NDMULT:NDMULT), &
                     IFLAGQ(0:NDMULT,-NDMULT:NDMULT)
     COMMON                                           &
            /QMULTI/ QMUL_N(0:NDMULT,-NDMULT:NDMULT), &
                     QMUL_P(0:NDMULT,-NDMULT:NDMULT), &
                     QMUL_T(0:NDMULT,-NDMULT:NDMULT)

     If(modeMPI_deformation.Eq.2.And.optimizePath.Eq.1.And.IERROR.Eq.0) Then

        ! Readjustement only for extra constraints not included in the original run
        Do i=1,numberConstraints
           LAMACT=qLambda_restart(i)
           MIUACT=qMiu_restart(i)
           QOLDIE=QASKED(LAMACT,MIUACT) ! Original value requested
           QSHIFT=QMUL_T(LAMACT,MIUACT) - qinitValue_restart(i) &
                                   -0.5D0*(qFinValue_restart(i) - qinitValue_restart(i))
           QASKED(LAMACT,MIUACT)=QSHIFT+QASKED(LAMACT,MIUACT) ! new value centered on unconstrained result
        End Do

        ! Displaying a message with the new values of the constraints
        Write(NFIPRI,'("*",77X,"*")')
        ITITLE=0
        Do lambda=1,NMUCON
           Do miu=-lambda,lambda
              If(IFLAGQ(lambda,miu).Eq.1) Then
                 If(ITITLE.Eq.0) Then
                     Write(NFIPRI,'(79("*"),/,"*",77X,"*")')
                     Write(NFIPRI,'("*",1X,"UPDATED MULTIPOLE CONSTRAINTS: ",&
                    &     "LAMBDA=",I2,"  MIU=",I2,2X,"MOMENT=",F7.3,12X,"*")') &
                           lambda,miu,QASKED(lambda,miu)
                     ITITLE=1
                 Else
                     Write(NFIPRI,'("*",32X,"LAMBDA=",I2,"  MIU=",I2,2X,"MOMENT=",F7.3,12X,"*")') &
                           lambda,miu,QASKED(lambda,miu)
                 End If
              End If
           End Do
        End Do

     End If

   End Subroutine mpi_pathExploration

   !---------------------------------------------------------------------
   !>  Subroutine reading the HFODD record file. It is a modified version
   !>  of routine RECORD(), where storage of big arrays is made through
   !>  allocatable Fortran arrays.
   !---------------------------------------------------------------------
   Subroutine read_hfodd(NFIREP,FILREP,NUMITE,NXHERM,NYHERM,NZHERM,IPCONT, &
                         ILCONT,IACONT,IMCONT,IRCONT,IRENMA,IRENIN,LIPKIN, &
                         LIPKIP,REFERN,REFERP,REDELN,REDELP,REFE2N,REFE2P, &
                         IFIBLN,INIBLN,IFIBLP,INIBLP,NMUMAX,NSIMAX,NMUCON, &
                                                            ISHIFT,IERROR)

     Use hfodd_sizes

     Integer :: currentVersion, minimalVersion, mpi_rank, mpi_err, iq
     Integer :: IX,IY,IZ,K,L,LDBASE,NFIPRI,IFNECK,IFIRED,LAMBDA,MIU, &
                NFIREP,NUMITE,NXHERM,NYHERM,NZHERM,MDMULT,IPCONT,    &
                ILCONT,IACONT,IMCONT,IRCONT,IRENMA,IRENIN,LIPKIN,    &
                LIPKIP,IFIBLN,INIBLN,IFIBLP,INIBLP,NMUMAX,NSIMAX,    &
                NMUCON,ISHIFT,IERROR,IVERIN,MXHERM,MYHERM,MZHERM,    &
                MMUMAX,MSIMAX,MDBASE,ICHARG,IREVER,KARTEZ,IBASE,ii
     Integer, Allocatable :: IFLALQ(:,:),IFLALS(:,:),IDEBLO(:)

     character(Len=68) :: FILREP

     Real :: REFERN,REFERP,REDELN,REDELP,REFE2N,REFE2P,Q0NECK,G_NECK,DUMMY
     Real, Allocatable :: ANGU_N(:),ANGU_P(:),ANGU_T(:),SPIN_N(:),SPIN_P(:),SPIN_T(:)
     Real, Allocatable :: QMUL_N(:,:),QMUL_P(:,:),QMUL_T(:,:),SMUL_N(:,:),SMUL_P(:,:), &
                          SMUL_T(:,:),GALMUQ(:,:),QLINEA(:,:),GALMUS(:,:),SLINEA(:,:), &
                          RALMUQ(:,:),RALMUS(:,:)
     Real, Allocatable :: FERINI(:),DELINI(:),FE2INI(:),RSHIFT(:),HBMREN(:),ROTREN(:), &
                          HBMREA(:),ROTREA(:)
     Real, Allocatable :: DLINSN(:),DLINSP(:),DLINST(:),ELINSN(:),ELINSP(:),ELINST(:), &
                          TLINSN(:),TLINSP(:),TLINST(:),ALINLN(:),ALINLP(:),ALINLT(:), &
                          DROTSN(:),DROTSP(:),DROTST(:),EROTSN(:),EROTSP(:),EROTST(:), &
                          TROTSN(:),TROTSP(:),TROTST(:),AROTLN(:),AROTLP(:),AROTLT(:)
     Real, Allocatable :: VN_MAS(:,:,:),VN_CEN(:,:,:),VP_MAS(:,:,:),VP_CEN(:,:,:), &
                          DPRRHO(:,:,:),DNRRHO(:,:,:)
     Real, Allocatable :: VN_KIS(:,:,:,:),VN_SPI(:,:,:,:),VN_CUR(:,:,:,:), &
                          VP_KIS(:,:,:,:),VP_SPI(:,:,:,:),VP_CUR(:,:,:,:)
     Real, Allocatable :: VN_SOR(:,:,:,:,:),VP_SOR(:,:,:,:,:)

     Complex :: C_ZERO
     Complex, Allocatable :: DN_RHO(:,:,:),DN_TAU(:,:,:),DN_LPR(:,:,:),DN_DIV(:,:,:), &
                             DP_RHO(:,:,:),DP_TAU(:,:,:),DP_LPR(:,:,:),DP_DIV(:,:,:),&
                             WN_CEN(:,:,:),WP_CEN(:,:,:),WAVBLO(:,:,:)
     Complex, Allocatable :: DN_SPI(:,:,:,:),DN_KIS(:,:,:,:),DN_GRR(:,:,:,:),DN_LPS(:,:,:,:),&
                             DN_ROS(:,:,:,:),DN_ROC(:,:,:,:),DN_CUR(:,:,:,:),                &
                             DP_SPI(:,:,:,:),DP_KIS(:,:,:,:),DP_GRR(:,:,:,:),DP_LPS(:,:,:,:),&
                             DP_ROS(:,:,:,:),DP_ROC(:,:,:,:),DP_CUR(:,:,:,:)
     Complex, Allocatable :: DN_SCU(:,:,:,:,:),DN_DES(:,:,:,:,:),DP_SCU(:,:,:,:,:),DP_DES(:,:,:,:,:)

     COMMON                 &
            /DIMENS/ LDBASE
     COMMON                 &
            /CFIPRI/ NFIPRI
     COMMON                 &
            /NCKFLA/ IFNECK

     C_ZERO=Cmplx(0.0D0,0.0D0)

     currentVersion = 16; minimalVersion =  4

     Call mpi_comm_rank(MPI_COMM_WORLD,mpi_rank,mpi_err)

     Open(UNIT=NFIREP,FILE=FILREP,STATUS='OLD',FORM='UNFORMATTED',IOSTAT=IERROR)
     If(IERROR.Ne.0) Then
         Write(NFIPRI,'("read_hfodd() - Error in opening file ",i5)') NFIREP
         Return
     End If

     IFIRED=0

     Read(NFIREP,IOSTAT=IERROR) IVERIN
     If(IERROR.Ne.0) Then
         Write(NFIPRI,'("read_hfodd() - error in reading IVERIN")')
         Return
     End If

     If(.Not.Allocated(WN_CEN)) &
        Allocate(WN_CEN(NDXHRM,NDYHRM,NDZHRM),VN_MAS(NDXHRM,NDYHRM,NDZHRM),VN_CEN(NDXHRM,NDYHRM,NDZHRM),&
                 WP_CEN(NDXHRM,NDYHRM,NDZHRM),VP_MAS(NDXHRM,NDYHRM,NDZHRM),VP_CEN(NDXHRM,NDYHRM,NDZHRM))
     If(.Not.Allocated(VN_KIS)) &
        Allocate(VN_KIS(NDXHRM,NDYHRM,NDZHRM,NDKART),VN_SPI(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                 VN_CUR(NDXHRM,NDYHRM,NDZHRM,NDKART),VP_KIS(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                 VP_SPI(NDXHRM,NDYHRM,NDZHRM,NDKART),VP_CUR(NDXHRM,NDYHRM,NDZHRM,NDKART))
     If(.Not.Allocated(VN_SOR)) &
        Allocate(VN_SOR(NDXHRM,NDYHRM,NDZHRM,NDKART,NDKART),VP_SOR(NDXHRM,NDYHRM,NDZHRM,NDKART,NDKART))

     If(.Not.Allocated(DN_RHO)) &
        Allocate(DN_RHO(NDXHRM,NDYHRM,NDZHRM),DN_TAU(NDXHRM,NDYHRM,NDZHRM),DN_LPR(NDXHRM,NDYHRM,NDZHRM),&
                 DN_DIV(NDXHRM,NDYHRM,NDZHRM),DP_RHO(NDXHRM,NDYHRM,NDZHRM),DP_TAU(NDXHRM,NDYHRM,NDZHRM),&
                 DP_LPR(NDXHRM,NDYHRM,NDZHRM),DP_DIV(NDXHRM,NDYHRM,NDZHRM))
     If(.Not.Allocated(DN_SPI)) &
        Allocate(DN_SPI(NDXHRM,NDYHRM,NDZHRM,NDKART),DN_KIS(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                 DN_GRR(NDXHRM,NDYHRM,NDZHRM,NDKART),DN_LPS(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                 DN_ROS(NDXHRM,NDYHRM,NDZHRM,NDKART),DN_ROC(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                 DN_CUR(NDXHRM,NDYHRM,NDZHRM,NDKART))
     If(.Not.Allocated(DP_SPI)) &
        Allocate(DP_SPI(NDXHRM,NDYHRM,NDZHRM,NDKART),DP_KIS(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                 DP_GRR(NDXHRM,NDYHRM,NDZHRM,NDKART),DP_LPS(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                 DP_ROS(NDXHRM,NDYHRM,NDZHRM,NDKART),DP_ROC(NDXHRM,NDYHRM,NDZHRM,NDKART),&
                 DP_CUR(NDXHRM,NDYHRM,NDZHRM,NDKART))
     If(.Not.Allocated(DN_SCU)) &
        Allocate(DN_SCU(NDXHRM,NDYHRM,NDZHRM,NDKART,NDKART),DN_DES(NDXHRM,NDYHRM,NDZHRM,NDKART,NDKART),&
                 DP_SCU(NDXHRM,NDYHRM,NDZHRM,NDKART,NDKART),DP_DES(NDXHRM,NDYHRM,NDZHRM,NDKART,NDKART))

     If(.Not.Allocated(ANGU_N)) &
        Allocate(ANGU_N(1:NDKART),ANGU_P(1:NDKART),ANGU_T(1:NDKART),&
                 SPIN_N(1:NDKART),SPIN_P(1:NDKART),SPIN_T(1:NDKART))

     If(.Not.Allocated(QMUL_N)) &
        Allocate(QMUL_N(0:NDMULT,-NDMULT:NDMULT),QMUL_P(0:NDMULT,-NDMULT:NDMULT),QMUL_T(0:NDMULT,-NDMULT:NDMULT),&
                 SMUL_N(0:NDMULT,-NDMULT:NDMULT),SMUL_P(0:NDMULT,-NDMULT:NDMULT),SMUL_T(0:NDMULT,-NDMULT:NDMULT),&
                 GALMUQ(0:NDMULT,-NDMULT:NDMULT),QLINEA(0:NDMULT,-NDMULT:NDMULT),RALMUQ(0:NDMULT,-NDMULT:NDMULT),&
                 GALMUS(0:NDMULT,-NDMULT:NDMULT),SLINEA(0:NDMULT,-NDMULT:NDMULT),RALMUS(0:NDMULT,-NDMULT:NDMULT))
     If(.Not.Allocated(IFLALQ)) &
        Allocate(IFLALQ(0:NDMULT,-NDMULT:NDMULT),IFLALS(0:NDMULT,-NDMULT:NDMULT))

     If(.Not.Allocated(FERINI)) &
        Allocate(FERINI(0:NDISOS),DELINI(0:NDISOS),FE2INI(0:NDISOS),RSHIFT(1:NDKART),&
                 HBMREN(1:NDKART),ROTREN(1:NDKART),HBMREA(1:NDKART),ROTREA(1:NDKART))

     If(.Not.Allocated(DLINSN)) &
        Allocate(DLINSN(0:NDKART),DLINSP(0:NDKART),DLINST(0:NDKART), &
                 ELINSN(0:NDKART),ELINSP(0:NDKART),ELINST(0:NDKART), &
                 TLINSN(0:NDKART),TLINSP(0:NDKART),TLINST(0:NDKART), &
                 ALINLN(0:NDKART),ALINLP(0:NDKART),ALINLT(0:NDKART), &
                 DROTSN(0:NDKART),DROTSP(0:NDKART),DROTST(0:NDKART), &
                 EROTSN(0:NDKART),EROTSP(0:NDKART),EROTST(0:NDKART), &
                 TROTSN(0:NDKART),TROTSP(0:NDKART),TROTST(0:NDKART), &
                 AROTLN(0:NDKART),AROTLP(0:NDKART),AROTLT(0:NDKART))

     If(.Not.Allocated(WAVBLO)) &
        Allocate(WAVBLO(1:NDBASE,0:NDISOS,0:NDREVE))
     If(.Not.Allocated(IDEBLO)) &
        Allocate(IDEBLO(0:NDISOS))

     If(IVERIN.Ge.minimalVersion.AND.IVERIN.Le.currentVersion) Then

        Read(NFIREP,IOSTAT=IERROR) NUMITE
        If(IERROR.Ne.0) Then
           Write(NFIPRI,'("read_hfodd() - Error in reading NUMITE")')
           Return
        End If

        Read(NFIREP,IOSTAT=IERROR) MXHERM,MYHERM,MZHERM
        If(IERROR.Ne.0) Then
           Write(NFIPRI,'("read_hfodd() - Error in reading MXHERM")')
           Return
        End If

        If(MXHERM.Ne.NXHERM.OR.MYHERM.Ne.NYHERM.OR.MZHERM.Ne.NZHERM) Then
           Write(NFIPRI,'(1X,19("/"),"  WRONG NXHERM,NYHERM,NZHERM:",3I3,2X,19("/"))') MXHERM,MYHERM,MZHERM
           Stop '  WRONG NXHERM,NYHERM,NZHERM IN RECORD'
        End If

        Read(NFIREP,IOSTAT=IERROR) (((VN_MAS(IX,IY,IZ),VN_CEN(IX,IY,IZ),IX=1,NXHERM),IY=1,NYHERM),IZ=1,NZHERM)
        If(IERROR.Ne.0) Then
           Write(NFIPRI,'("read_hfodd() - Error in reading VN_MAS")')
           Return
        End If

        Read(NFIREP,IOSTAT=IERROR) ((((VN_KIS(IX,IY,IZ,L),VN_SPI(IX,IY,IZ,L),VN_CUR(IX,IY,IZ,L), &
                                               IX=1,NXHERM),IY=1,NYHERM),IZ=1,NZHERM),L=1,NDKART)
        If(IERROR.Ne.0) Then
           Write(NFIPRI,'("read_hfodd() - Error in reading VN_KIS")')
           Return
        End If

        Read(NFIREP,IOSTAT=IERROR) (((((VN_SOR(IX,IY,IZ,L,K),IX=1,NXHERM),IY=1,NYHERM),IZ=1,NZHERM),L=1,NDKART),&
                                                                                                    K=1,NDKART)
        If(IERROR.Ne.0) Then
           Write(NFIPRI,'("read_hfodd() - Error in reading VN_SOR")')
           Return
        End If

        Read(NFIREP,IOSTAT=IERROR) (((VP_MAS(IX,IY,IZ),VP_CEN(IX,IY,IZ),IX=1,NXHERM),IY=1,NYHERM),IZ=1,NZHERM)
        If(IERROR.Ne.0) Then
           Write(NFIPRI,'("read_hfodd() - Error in reading VP_MAS")')
           Return
        End If

        Read(NFIREP,IOSTAT=IERROR) ((((VP_KIS(IX,IY,IZ,L),VP_SPI(IX,IY,IZ,L),VP_CUR(IX,IY,IZ,L), &
                                               IX=1,NXHERM),IY=1,NYHERM),IZ=1,NZHERM),L=1,NDKART)
        If(IERROR.Ne.0) Then
           Write(NFIPRI,'("read_hfodd() - Error in reading VP_KIS")')
           Return
        End If

        Read(NFIREP,IOSTAT=IERROR) (((((VP_SOR(IX,IY,IZ,L,K),IX=1,NXHERM),IY=1,NYHERM),IZ=1,NZHERM),L=1,NDKART),&
                                                                                                    K=1,NDKART)
        If(IERROR.Ne.0) Then
           Write(NFIPRI,'("read_hfodd() - Error in reading VP_SOR")')
           Return
        End If

        If(IVERIN.Le.6) Then

           If(IVERIN.Gt.5) Then
               Read(NFIREP,IOSTAT=IERROR) (((DNRRHO(IX,IY,IZ),IX=1,NXHERM),IY=1,NYHERM),IZ=1,NZHERM)
               If(IERROR.Ne.0) Then
                   Write(NFIPRI,'("read_hfodd() - Error in reading DNRRHO")')
                   Return
               End If
           Else
               DNRRHO(:,:,:)=0.0D0
           End If

           Read(NFIREP,IOSTAT=IERROR) (((DPRRHO(IX,IY,IZ),IX=1,NXHERM),IY=1,NYHERM),IZ=1,NZHERM)
           If(IERROR.Ne.0) Then
              Write(NFIPRI,'("read_hfodd() - Error in reading DPRRHO")')
              Return
           End If
           DN_RHO(:,:,:)=DNRRHO(:,:,:)*cmplx(1.0d0,0.0d0)
           DP_RHO(:,:,:)=DPRRHO(:,:,:)*cmplx(1.0d0,0.0d0)

        Else

           Read(NFIREP,IOSTAT=IERROR) (((DN_RHO(IX,IY,IZ),IX=1,NXHERM),IY=1,NYHERM),IZ=1,NZHERM)
           If(IERROR.Ne.0) Then
              Write(NFIPRI,'("read_hfodd() - Error in reading DN_RHO")')
              Return
           End If

           Read(NFIREP,IOSTAT=IERROR)(((DP_RHO(IX,IY,IZ),IX=1,NXHERM),IY=1,NYHERM),IZ=1,NZHERM)
           If(IERROR.Ne.0) Then
              Write(NFIPRI,'("read_hfodd() - Error in reading DP_RHO")')
              Return
           End If

        End If

        Read(NFIREP,IOSTAT=IERROR) MMUMAX
        If(IERROR.Ne.0) Then
           Write(NFIPRI,'("read_hfodd() - Error in reading NMUMAX")')
           Return
        End If

        If(MMUMAX.Lt.NMUMAX) Then

           Write(NFIPRI,'(1X,20("/"),"  WARNING:  TOO FEW MULTIPOLE MOMENTS",1X,20("/"))')
           Write(NFIPRI,'(1X,20("/"),"  FOUND ON THE REPLAY FILE.          ",1X,20("/"))')
           Write(NFIPRI,'(1X,20("/"),"  REQUESTED MAXIMUM  MULTIPOLARITY=",I2,1X,20("/"))') NMUMAX
           Write(NFIPRI,'(1X,20("/"),"  FOUND     MAXIMUM  MULTIPOLARITY=",I2,1X,20("/"))') MMUMAX
           Write(NFIPRI,'(1X,20("/"),"  A CONSTRAINT ON A MISSING MULTIPOLE",1X,20("/"))')
           Write(NFIPRI,'(1X,20("/"),"  WILL BE  INCLUDED  ONLY  AFTER  THE",1X,20("/"))')
           Write(NFIPRI,'(1X,20("/"),"  FIRST DIAGONALIZATION IS  PERFORMED",1X,20("/"))')

           Do LAMBDA=MMUMAX+1,NMUMAX
              Do MIU=-LAMBDA,LAMBDA
                 QMUL_N(LAMBDA,MIU)=0.0D0
                 QMUL_P(LAMBDA,MIU)=0.0D0
                 QMUL_T(LAMBDA,MIU)=0.0D0
              End Do
           End Do

        End If

        Read(NFIREP,IOSTAT=IERROR) ((QMUL_N(LAMBDA,MIU),QMUL_P(LAMBDA,MIU),QMUL_T(LAMBDA,MIU), &
                                                    MIU=-LAMBDA,LAMBDA),LAMBDA=0,MMUMAX)
        If(IERROR.Ne.0) Then
           Write(NFIPRI,'("read_hfodd() - Error in reading QMUL_N")')
           Return
        End If

        Do LAMBDA=0,NSIMAX
           Do MIU=-LAMBDA,LAMBDA
              SMUL_N(LAMBDA,MIU)=0.0D0
              SMUL_P(LAMBDA,MIU)=0.0D0
              SMUL_T(LAMBDA,MIU)=0.0D0
           End Do
        End Do

        If(IVERIN.Gt.4) Then

           Read(NFIREP,IOSTAT=IERROR) MSIMAX
           If(IERROR.Ne.0) Then
              Write(NFIPRI,'("read_hfodd() - Error in reading MSIMAX")')
              Return
           End If

           If(MSIMAX.Lt.NSIMAX) Then

              Write(NFIPRI,'(1X,20(1H/),"  WARNING:  TOO FEW SURFACE   MOMENTS",1X,20(1H/))')
              Write(NFIPRI,'(1X,20(1H/),"  FOUND ON THE REPLAY FILE.          ",1X,20(1H/))')
              Write(NFIPRI,'(1X,20(1H/),"  REQUESTED MAXIMUM  MULTIPOLARITY=",I2,1X,20(1H/))') NSIMAX
              Write(NFIPRI,'(1X,20(1H/),"  FOUND     MAXIMUM  MULTIPOLARITY=",I2,1X,20(1H/))') MSIMAX
              Write(NFIPRI,'(1X,20(1H/),"  A CONSTRAINT ON A MISSING MULTIPOLE",1X,20(1H/))')
              Write(NFIPRI,'(1X,20(1H/),"  WILL BE  INCLUDED  ONLY  AFTER  THE",1X,20(1H/))')
              Write(NFIPRI,'(1X,20(1H/),"  FIRST DIAGONALIZATION IS  PERFORMED",1X,20(1H/))')

              Do LAMBDA=MSIMAX+1,NSIMAX
                 Do MIU=-LAMBDA,LAMBDA
                    SMUL_N(LAMBDA,MIU)=0.0D0
                    SMUL_P(LAMBDA,MIU)=0.0D0
                    SMUL_T(LAMBDA,MIU)=0.0D0
                 End Do
              End Do

           End If

           Read(NFIREP,IOSTAT=IERROR) ((SMUL_N(LAMBDA,MIU),SMUL_P(LAMBDA,MIU),SMUL_T(LAMBDA,MIU), &
                                              MIU=-LAMBDA,LAMBDA),LAMBDA=0,MSIMAX)
           If(IERROR.Ne.0) Then
              Write(NFIPRI,'("read_hfodd() - Error in reading SMUL_N")')
              Return
           End If

        End If

        Read(NFIREP,IOSTAT=IERROR) ISHIFT,(RSHIFT(L),L=1,NDKART)
        If(IERROR.Ne.0) Then
           Write(NFIPRI,'("read_hfodd() - Error in reading ISHIFT")')
           Return
        End If

        Read(NFIREP,IOSTAT=IERROR) (ANGU_N(L),ANGU_P(L),ANGU_T(L),SPIN_N(L),SPIN_P(L),SPIN_T(L),L=1,NDKART)
        If(IERROR.Ne.0) Then
           Write(NFIPRI,'("read_hfodd() - Error in reading ANGU_N")')
           Return
        End If

        Read(NFIREP,IOSTAT=IERROR) (((WN_CEN(IX,IY,IZ),IX=1,NXHERM),IY=1,NYHERM),IZ=1,NZHERM)
        If(IERROR.Ne.0) Then
           Write(NFIPRI,'("read_hfodd() - Error in reading WN_CEN")')
           Return
        End If

        Read(NFIREP,IOSTAT=IERROR) (((WP_CEN(IX,IY,IZ),IX=1,NXHERM),IY=1,NYHERM),IZ=1,NZHERM)
        If(IERROR.Ne.0) Then
          Write(NFIPRI,'("read_hfodd() - Error in reading WP_CEN")')
          Return
        End If

        Read(NFIREP,IOSTAT=IERROR) REFERN,REDELN,REFERP,REDELP
        If(IERROR.Ne.0) Then
           Write(NFIPRI,'("read_hfodd() - Error in reading REFERN")')
           Return
        End If

        ! LN lambda_2 (IVERIN=8)
        If(IVERIN.Gt.7) Then
           Read(NFIREP,IOSTAT=IERROR) REFE2N,REFE2P
           If(IERROR.Ne.0) Then
              Write(NFIPRI,'("read_hfodd() - Error in reading REFE2N")')
              Return
           End If
        Else
           REFE2N=FE2INI(0); REFE2P=FE2INI(1)
           If((LIPKIN.EQ.1.OR.LIPKIP.EQ.1).AND.ILCONT.EQ.1) Then
               Write(NFIPRI,'(1X,13("/")," WARNING: THE LIPKIN-NOGAMI LAMBDA2S HAVE NOT BEEN",1X,13("/"))')
               Write(NFIPRI,'(1X,13("/")," FOUND ON THE  OLD  REPLAY  FILE  VERSION  NO. =",I2,1X,13("/"))') IVERIN
               Write(NFIPRI,'(1X,13("/")," THE  LN  CORRECTIONS CANNOT BE SMOOTHLY CONTINUED",1X,13("/"))')
           End If
        End If

        ! Blocked state wave-function (IVERIN=9)
        IDEBLO(0)=0; IDEBLO(1)=0
        If(IVERIN.Gt.8) Then

           Read(NFIREP,IOSTAT=IERROR) MDBASE,(IDEBLO(ICHARG),ICHARG=0,NDISOS)
           If(IERROR.Ne.0) Then
              Write(NFIPRI,'("read_hfodd() - Error in reading MDBASE")')
              Return
           End If

           If(MDBASE.Le.NDBASE) Then
              Read(NFIREP,IOSTAT=IERROR) (((WAVBLO(IBASE,ICHARG,IREVER),IBASE=1,MDBASE),ICHARG=0,NDISOS),&
                                                                IREVER=0,NDREVE)
              If(IERROR.Ne.0) Then
                 Write(NFIPRI,'("read_hfodd() - Error in reading WAVBLO")')
                 Return
              End If
           Else
              Read(NFIREP,IOSTAT=IERROR)
              If(IERROR.Ne.0) Then
                 Write(NFIPRI,'("read_hfodd() - Error in reading ...nothing...")')
                 Return
              End If
           End If

           If(((IFIBLN.EQ.1.AND.INIBLN.EQ.0).OR.(IFIBLP.EQ.1.AND.INIBLP.EQ.0)).AND. &
                (MDBASE.Ne.LDBASE.OR.MDBASE.Gt.NDBASE)) Then
                 Write(NFIPRI,'(1X,10("/")," INCOMPATIBLE DIMENSIONS FOUND ON THE RECORD FILE:",1X,10("/"))')
                 Write(NFIPRI,'(1X,10("/")," LDBASE=",I4," MDBASE=",I4," NDBASE=",I4,1X,10("/"))') LDBASE,MDBASE,NDBASE
                 Stop 'MDBASE<>LDBASE OR MDBASE>NDBASE IN RECORD'
           End If

        End If

        If(IVERIN.Ge.10) Then
           Read(NFIREP,IOSTAT=IERROR) (ALINLT(KARTEZ),KARTEZ=0,NDKART)
           If(IERROR.Ne.0) Then
              Write(NFIPRI,'("read_hfodd() - Error in reading ALINLT")')
              Return
           End If
        Else
           ALINLT(:)=0.0D0
           If(IMCONT.EQ.1) Then
              Write(NFIPRI,'(1X,13(1H/)," WARNING: THE AVERAGE LINEAR MOMENTA HAVE NOT BEEN",1X,13(1H/))')
              Write(NFIPRI,'(1X,13(1H/)," FOUND ON THE  OLD  REPLAY  FILE  VERSION  NO. =",I2,1X,13(1H/))') IVERIN
              Write(NFIPRI,'(1X,13(1H/)," THE EXACT CM CORREC. CANNOT BE SMOOTHLY CONTINUED",1X,13(1H/))')
           End If
        End If

        If(IVERIN.Ge.11) Then
           Read(NFIREP,IOSTAT=IERROR) (AROTLT(KARTEZ),KARTEZ=0,NDKART)
           If(IERROR.Ne.0) Then
              Write(NFIPRI,'("read_hfodd() - Error in reading AROTLT")')
              Return
           End If
        Else
           AROTLT(:)=0.0D0
           If(IRCONT.Eq.1) Then
              Write(NFIPRI,'(1X,13(1H/)," WARNING: THE AVERAGE ANGULAR MOMENTA HAVE NOT BEEN",1X,13(1H/))')
              Write(NFIPRI,'(1X,13(1H/)," FOUND ON  THE  OLD  REPLAY  FILE  VERSION  NO. =",I2,1X,13(1H/))') IVERIN
              Write(NFIPRI,'(1X,13(1H/)," THE EXACT ROT CORREC. CANNOT BE SMOOTHLY CONTINUED",1X,13(1H/))')
           End If
        End If

        If(IVERIN.Ge.12) Then

           Read(NFIREP,IOSTAT=IERROR) MDMULT
           If(IERROR.Ne.0) Then
              Write(NFIPRI,'("read_hfodd() - Error in reading MDMULT")')
              Return
           End If

           RALMUQ(:,:)=0.0D0; RALMUS(:,:)=0.0D0

           If(MDMULT.Ne.NDMULT) Then

              Write(NFIPRI,'(1X,20(1H/),"  WARNING: NOT ALL LINEAR CONSTRAINTS",1X,20(1H/))')
              Write(NFIPRI,'(1X,20(1H/),"  FOUND ON THE REPLAY FILE.          ",1X,20(1H/))')
              Write(NFIPRI,'(1X,20(1H/),"  REQUESTED MAXIMUM  MULTIPOLARITY=",I2,1X,20(1H/))') NDMULT
              Write(NFIPRI,'(1X,20(1H/),"  FOUND     MAXIMUM  MULTIPOLARITY=",I2,1X,20(1H/))') MDMULT
              Write(NFIPRI,'(1X,20(1H/),"  A CONSTRAINT ON A MISSING MULTIPOLE",1X,20(1H/))')
              Write(NFIPRI,'(1X,20(1H/),"  WILL BE  INCLUDED  ONLY  AFTER  THE",1X,20(1H/))')
              Write(NFIPRI,'(1X,20(1H/),"  FIRST DIAGONALIZATION IS  PERFORMED",1X,20(1H/))')

              RALMUQ(:,:)=QLINEA(:,:); RALMUS(:,:)=SLINEA(:,:)

              Read(NFIREP,IOSTAT=IERROR) ((DUMMY,DUMMY,MIU=-LAMBDA,LAMBDA),LAMBDA=0,MDMULT)
              If(IERROR.Ne.0) Then
                 Write(NFIPRI,'("read_hfodd() - Error in reading DUMMY")')
                 Return
              End If

           Else

               Read(NFIREP,IOSTAT=IERROR) ((RALMUQ(LAMBDA,MIU),RALMUS(LAMBDA,MIU),MIU=-LAMBDA,LAMBDA),&
                                                               LAMBDA=0,MDMULT)
               If(IERROR.Ne.0) Then
                  Write(NFIPRI,'("read_hfodd() - Error in reading RALMUQ")')
                  Return
               End If

               ! If not allocated, allocate; if already allocated, deallocate first and re-allocate (
               ! (enables vayring number of constraints)
               If(Allocated(flex_lambda)) Deallocate(flex_lambda,flex_miu,flex_qlm)
               Allocate(flex_lambda(numberConstraints),flex_miu(numberConstraints),flex_qlm(numberConstraints))
               Do iq=1,numberConstraints
                  lambda = qLambda(iq); miu = qMiu(iq)
                  flex_lambda(iq)=lambda
                  flex_miu(iq)=miu
                  flex_qlm(iq)=QMUL_T(lambda,miu)
               End Do

           End If

        Else

            RALMUQ(:,:)=QLINEA(:,:); RALMUS(:,:)=SLINEA(:,:)

            If(IACONT.EQ.1) Then
               Write(NFIPRI,'(1X,13(1H/)," WARNING: THE LINEAR CONSTRAINTS  H A V E  NOT BEEN",1X,13(1H/))')
               Write(NFIPRI,'(1X,13(1H/)," FOUND ON  THE  OLD  REPLAY  FILE  VERSION  NO. =",I2,1X,13(1H/))') IVERIN
               Write(NFIPRI,'(1X,13(1H/)," THE EXACT AUGM.LAGRA. CANNOT BE SMOOTHLY CONTINUED",1X,13(1H/))')
            End If

        End If

        If(IVERIN.Ge.13) Then

            Read(NFIREP,IOSTAT=IERROR) (HBMREA(KARTEZ),KARTEZ=1,NDKART)
            If(IERROR.Ne.0) Then
               Write(NFIPRI,'("read_hfodd() - Error in reading HBMREA")')
               Return
            End If

            Read(NFIREP,IOSTAT=IERROR) (ROTREA(KARTEZ),KARTEZ=1,NDKART)
            If(IERROR.Ne.0) Then
               Write(NFIPRI,'("read_hfodd() - Error in reading ROTREA")')
               Return
            End If

            If(IMCONT.EQ.1.AND.IRENMA.Ge.1) Then
               HBMREN(:)=HBMREA(:)
            End If

            If(IRCONT.EQ.1.AND.IRENIN.Ge.1) Then
               ROTREN(:)=ROTREA(:)
            End If

        Else

            If(IMCONT.EQ.1.AND.IRENMA.Ge.1) Then
               Write(NFIPRI,'(1X,13(1H/)," WARNING:  THE  RENORMALIZED MASSES  HAVE NOT BEEN",1X,13(1H/))')
               Write(NFIPRI,'(1X,13(1H/)," FOUND ON  THE  OLD  REPLAY FILE  VERSION  NO. =",I2,1X,13(1H/))') IVERIN
               Write(NFIPRI,'(1X,13(1H/)," THE EXACT CM CORREC. CANNOT BE SMOOTHLY CONTINUED",1X,13(1H/))')
            End If

            If(IRCONT.EQ.1.AND.IRENIN.Ge.1) Then
               Write(NFIPRI,'(1X,13(1H/)," WARNING:  THE   RENORMALIZED INERTIA HAVE NOT BEEN",1X,13(1H/))')
               Write(NFIPRI,'(1X,13(1H/)," FOUND  ON  THE  OLD  REPLAY FILE  VERSION  NO. =",I2,1X,13(1H/))') IVERIN
               Write(NFIPRI,'(1X,13(1H/)," THE EXACT ROT CORREC. CANNOT BE SMOOTHLY CONTINUED",1X,13(1H/))')
            End If

        End If

        ! Linear Constraints for the neck (IVERIN=14)
        If(IVERIN.Ge.14) Then
            Read(NFIREP,IOSTAT=IERROR) G_NECK
            If(IERROR.Ne.0) Then
               Write(NFIPRI,'("read_hfodd() - Error in reading MDMULT")')
               Return
            End If
        Else
            G_NECK=0.0D0
            If(IFNECK.Ge.2) Then
               Write(NFIPRI,'(1X,13(1H/)," WARNING : THE  NECK  CONSTRAINTS  H A S  NOT  BEEN",1X,13(1H/))')
               Write(NFIPRI,'(1X,13(1H/)," FOUND ON  THE  OLD  REPLAY  FILE  VERSION  NO. =",I2,1X,13(1H/))') IVERIN
               Write(NFIPRI,'(1X,13(1H/)," THE NECK CONSTRAINT CAN NOT BE  SMOOTHLY CONTINUED",1X,13(1H/))')
            End If
        End If

        IFIRED=1

     End If

     If(IFIRED.EQ.0) Then
        Write(NFIPRI,'(/,1X,10(1H/),"  CANNOT read THE REPLAY FILE VERSION:",I3,2X,10(1H/),/)') IVERIN
        Stop '  WRONG IVERIN IN RECORD'
     End If

     If(IFIBLN.EQ.1.AND.INIBLN.EQ.0.AND.IDEBLO(0).Ne.1) Then
        Write(NFIPRI,'(1X,10(1H/)," NEUTRON SINGLE-PARTICLE WF NOT FOUND ON THE RECORD FILE",1X,10(1H/))')
        Write(NFIPRI,'(1X,10(1H/)," RERUN THE BLOCKING CALCULATION WITH INIBLN=1           ",1X,10(1H/))')
        Stop ' NEUTRON SINGLE-PARTICLE WF NOT FOUND IN RECORD'
     End If

     If(IFIBLP.EQ.1.AND.INIBLP.EQ.0.AND.IDEBLO(1).Ne.1) Then
        Write(NFIPRI,'(1X,10(1H/)," PROTON  SINGLE-PARTICLE WF NOT FOUND ON THE RECORD FILE",1X,10(1H/))')
        Write(NFIPRI,'(1X,10(1H/)," RERUN THE BLOCKING CALCULATION WITH INIBLP=1           ",1X,10(1H/))')
        Stop ' PROTON SINGLE-PARTICLE WF NOT FOUND IN RECORD'
     End If

   End Subroutine read_hfodd

   !---------------------------------------------------------------------
   !>  Subroutine reads lists all the files in 'directory_name', extracts
   !>  all files of the form 'HFODD_XXXXXXXX.REC' (with XXXXXXXX a number)
   !>  and checks that all files with the number between 1 and N are
   !>  present. It then broadcasts a boolean flag indicating if all files
   !>  are there, and a vector of integers giving the breakdown.
   !---------------------------------------------------------------------
   Subroutine check_directory(directory_name,N)

      Integer, INTENT(IN) :: N
      Character(Len=132), INTENT(IN) :: directory_name

      Integer, Parameter :: ndatin=85
      Integer :: i,exitstat,iocheck,ierr,mpi_rank,mpi_err,numero,count
      Character(Len=80) :: filein
      Character(Len=150) :: commande
      Character(Len=132) :: string

      mpi_rank=0; mpi_err=0

      ! Getting the rank of the current process
      Call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, mpi_err)

      If(mpi_rank.Eq.0) Then

         commande = 'ls '//Trim(directory_name)// ' > toto'; exitstat=0
         Call system(commande)
         If(debug.Ge.2) Write(6,'("mpi_rank = ",i10," - Executed command: ",a)') &
                                   mpi_rank,commande

         ! open file
         Open(ndatin,file='toto',status='old',form='formatted',iostat=ierr)
         If(ierr.Ne.0) Then
            Write(6,'("Error in opening file : toto")')
            Stop 'Error in check_directory - I/O'
         End if

         ! Read file, extract number of file if relevant, and fill table
         iocheck=0
         Do While(iocheck.Eq.0)
            Read(ndatin,*,Iostat=iocheck) string
            numero = 0
            Call get_numero(string,numero)
            If(debug.Ge.2) Write(6,'("numero = ",i10)') numero
            If(numero.Gt.0) table_fichiers(numero)=1
         End Do

         Close(ndatin)
         If(debug.Ge.1) Write(6,'("Done reading with the file")')

         filled=.True.; count=0
         Do i=1,N
            If(table_fichiers(i).Eq.0) Then
               count=count+1
            End If
         End Do
         If(count.Gt.0) filled=.False.

         If(debug.Ge.2) Then
            If(filled) Then
               Write(6,'("All files are there")')
            Else
               Write(6,'(i10," files are missing")') count
            End If
         End If

      End If

      ! Broadcast flag telling if array is filled and table of missing files
      Call MPI_Bcast(filled,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_err)
      Call MPI_Bcast(table_fichiers,N,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)

      ! Wait for all threads
      Call mpi_barrier(MPI_COMM_WORLD,mpi_err)

    End Subroutine check_directory

   !---------------------------------------------------------------------
   !>  Subroutine returns in 'numero' the file number from a string of
   !>  the form 'HFODD_XXXXXXXX.REC' (with XXXXXXXX the number).
   !---------------------------------------------------------------------
    Subroutine get_numero(string,numero)

      Character(Len=132), INTENT(IN) :: string
      Integer, INTENT(INOUT) :: numero
      Integer :: start
      Character(Len=4) :: dum_4
      Character(Len=6) :: dum_6

      start = Index(string, 'HFODD_', .False.)
      If (start.Ge.0) Then
          Read(string(start:start+18), '(A6,I8.8,A4)') dum_6,numero,dum_4
      End If

    End Subroutine get_numero

   !---------------------------------------------------------------------
   !>  Subroutine copies all files that have not been copied yet.
   !---------------------------------------------------------------------
    Subroutine fill_holes(N,IF_THO,LIPKIN,LIPKIP,ILCONT,IFCONT)

      Integer, INTENT(IN) :: N,IF_THO,LIPKIN,LIPKIP,ILCONT,IFCONT

      Integer :: jj,exitstat,mpi_rank,mpi_world,mpi_err
      Character(Len=150) :: commande

#if(USE_MANYCORES==1)
     Integer :: tribeID
     Integer :: numberMasters, numberHFODDproc, color,       &
                tribeRank, slaveRank, masterRank, worldRank, &
                tribeSize, slaveSize, masterSize, worldSize
     Integer :: worldGroup, groupMasters, mastersCOMM, &
                groupSlaves, slavesCOMM, tribeCOMM
     COMMON                                                 &
            /MPICOM/ worldGroup, groupMasters, groupSlaves, &
                     mastersCOMM, slavesCOMM, tribeCOMM
     COMMON                                                       &
            /MPIPRO/ numberMasters, numberHFODDproc, color,       &
                     tribeRank, slaveRank, masterRank, worldRank, &
                     tribeSize, slaveSize, masterSize, worldSize
#endif

      mpi_rank=0; mpi_world=1; mpi_err=0

      ! Getting the rank of the current process
      Call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, mpi_err)
      Call mpi_comm_size(MPI_COMM_WORLD, mpi_world,mpi_err)

      Do jj=1,N

#if(USE_MANYCORES==1)
         Do tribeID=0,tribeSize
            If(Mod(tribeRank,tribeSize).Eq.Mod(tribeID,tribeSize).And.color.Eq.Mod(jj-1,numberMasters)) Then
#else
         If(Mod(mpi_rank,mpi_world).Eq.Mod(jj-1,mpi_world)) Then
#endif

            If(table_fichiers(jj).Eq.0) Then

               commande = 'cp '//fichier_rec_old(jj)//' '//fichier_rec_new(jj); exitstat=0
               Call system(commande)
               If(debug.Ge.2) Write(6,'("fill_holes: mpi_rank = ",i10," jj = ",i10," - Executed command: ",a)') &
                                                     mpi_rank,jj,commande
               ! HFBTHO files
               If(IF_THO.Ge.1) Then
                  commande = 'cp '//fichier_tho_old(jj)//' '//fichier_tho_new(jj); exitstat=0
                  Call system(commande)
                  If(debug.Ge.2) Write(6,'("fill_holes: mpi_rank = ",i10," jj = ",i10," - Executed command: ",a)') &
                                                        mpi_rank,jj,commande
               End If
               ! Densities
               If((LIPKIN.Ge.1.Or.LIPKIP.Ge.1).And.ILCONT.Eq.1) Then
                  commande = 'cp '//fichier_lic_old(jj)//' '//fichier_lic_new(jj); exitstat=0
                  Call system(commande)
                  If(debug.Ge.2) Write(6,'("fill_holes: mpi_rank = ",i10," jj = ",i10," - Executed command: ",a)') &
                                                        mpi_rank,jj,commande
               End If
               ! Fields
               If(IFCONT.Eq.1) Then
                  commande = 'cp '//fichier_fic_old(jj)//' '//fichier_fic_new(jj); exitstat=0
                  Call system(commande)
                  If(debug.Ge.2) Write(6,'("fill_holes: mpi_rank = ",i10," jj = ",i10," - Executed command: ",a)') &
                                                        mpi_rank,jj,commande
               End If

            End If

#if(USE_MANYCORES==1)
            End If ! tribeRank vs. tribe ID
         End Do ! loop over tribe ID
#else
         End If ! mpi_rank versus deformation number
#endif
      End Do ! loop over deformation number

      ! Wait for all threads
      Call mpi_barrier(MPI_COMM_WORLD,mpi_err)

   End Subroutine fill_holes

 End Module hfodd_mpimanager



