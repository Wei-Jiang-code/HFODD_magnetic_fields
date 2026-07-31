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
 !                HFBTHO interface for HFODD                            !
 !                                                                      !
 ! Authors: N., Schunck, ORNL, LLNL                                     !
 !                                                                      !
 ! This module reads the binary file produced by HFBTHO and passes the  !
 ! matrices of the p.h. and p.p. fields to HFODD (main code). Basis     !
 ! transformations (from cylindrical to cartesian coordinates, simplex  !
 ! basis) are performed in the main HFODD. 	                        !
 !                                                                      !
 !                          First included in official release v249k    !
 !                          Last modification: April 30, 2012           !
 !                                                                      !
 !----------------------------------------------------------------------!

 Module hfodd_interface

   Implicit None

   Integer, Parameter :: inik=Kind(1), dnik=Kind(1.0)

   ! PUBLIC VARIABLES (PASSED TO/FROM, AND USED IN HFODD)
   Integer(inik), PUBLIC :: NmxTHO,NsizeHFBTHO,NsizeHFODD,maxLambda,ntx_tho,nnm_tho,nt_tho
   Integer(inik), PUBLIC :: nhhdim_tho,nhhdim2_tho,nhhdim3_tho,nhhdim4_tho
   Integer(inik), Allocatable, PUBLIC :: id_tho(:),nl_tho(:),ns_tho(:),Sigma2(:)
   Integer(inik), Allocatable, PUBLIC :: Nm_tho(:),nz_tho(:),nr_tho(:),ml_tho(:)
   Integer(inik), Allocatable, PUBLIC :: Nm_odd(:),nz_odd(:),nr_odd(:),ml_odd(:)
   Integer(inik), Allocatable, PUBLIC :: validate_HFBTHO(:,:,:,:),index_HFBTHO(:,:,:,:,:),index_HFODD(:,:,:,:,:)
   !
   Real(dnik), Allocatable, PUBLIC :: xmultLag(:),brin_tho(:)
   !
   Complex(dnik), Allocatable, PUBLIC :: HF_PHN(:, :, :, :),HF_PHP(:, :, :, :),HF_PPN(:, :, :, :),HF_PPP(:, :, :, :)
   !
   Real(dnik), Allocatable, PUBLIC :: THOPHN(:, :),THOPHP(:, :),THOPPN(:, :),THOPPP(:, :)

   ! PRIVATE VARIABLES
   Real(dnik) :: zero = 0.0_dnik, one=1.0_dnik, two = 2.0_dnik, half = 0.5_dnik
   Integer(inik), PRIVATE :: nbx_tho,nzx_tho,nrx_tho,nlx_tho,ndx_tho,nb_tho,nom_tho
   Integer(inik), Allocatable, PRIVATE :: npar_tho(:),ia_tho(:),ikb_tho(:),ipb_tho(:)

 Contains
   !==========================================================================================================
   !
   !==========================================================================================================
   Subroutine ReadHFBTHO(EFERMN,EFERMP,EFER2N,EFER2P,DELTAN,DELTAP,&
                         PRHO_N,PRHODN,PRHOSN,POWERN,PRHO_P,PRHODP,&
                         PRHOSP,POWERP,BETTHO,N00THO,NLIMIT,FILTHO,&
                         B0_THO,BZ_THO,BP_THO)
      ! Inputs/outputs
      Character(12) :: FILTHO
      Integer :: N00THO,NLIMIT
      Real(dnik) :: EFERMN,EFERMP,EFER2N,EFER2P,DELTAN,DELTAP
      Real(dnik) :: PRHO_N,PRHODN,PRHOSN,POWERN,PRHO_P,PRHODP,PRHOSP,POWERP
      Real(dnik) :: BETTHO,B0_THO,BZ_THO,BP_THO
      ! From HFBTHO
      Integer(inik) :: it,bloall1,kindhfb,itass,iqqmax,ilst1,l,ierr,counter
      Integer(inik) :: npr11,npr1,ngh1,ngl1,n001,nb1,nt1
      Integer(inik) :: i,ibasis,ibro,ib,nd,n1,n2,n1n2,n2n1,N_rz
      Integer(inik) :: nb,nbx,ntx,ndx,nqx,nb2x,nhhdim,nhhdim2,nhhdim3,nhhdim4
      Integer(inik) :: npr(2)      ! number of neutrons (1) protons (2)
      Integer(inik) :: lwin = 41   ! number of the unformated input unit
      Character(2)  :: nucname     ! name of the nucleus
      Character(30) :: welfile     ! name of the file e.g. 's070_050.hel'
      Real(dnik) :: b0,bz,bp,beta0,bet,v0(2),v1,dalf,pwi,amas,pi
      Real(dnik) :: si,etot,xmix,varmas,skass,tz(2)
      Real(dnik) :: decay,rmm3,amm3,bmm3,cmm3
      Real(dnik) :: varmasnz(2),pjmassnz(2),ass(2)
      Real(dnik) :: rms(3),ept(3),del(2),ala(2),ala2(2),alast(2)
      Real(dnik), Allocatable :: brin(:)
      Integer(inik) :: iblocking,bloall; Parameter(bloall = 200)
      Integer(inik), Dimension(0:bloall,2) :: bloblo,blo123=0,blok1k2=0
      Real(dnik),    Dimension(0:bloall,2) :: bloqpdif
      Integer(inik), Dimension(3) :: keyblo
      Integer(inik), Dimension(2) :: blocross,blomax,blo123d,blok1k2d,blocanon
      ! Local variables
      Logical :: file_exists,file_opened
      Integer(inik), Allocatable :: id(:),Omega2(:),Parity(:)
      Integer(inik) :: N, nr, nz, ml, ms, mj, lambda_max, mpi_rank, mpi_err
      Integer(inik) :: IndexBra, IndexKet, IndexStart
      Integer(inik) :: Indexx, Nperpe, i_spin, j_spin, icount, j, ms_actual, i2
      !---------------------------------------------------------------------!
      !  This subroutine reads .hel files generated by HFBTHO. Basically,   !
      !  it is a variation of subroutine inout() of HFBTHO. Differences are !
      !  that a number of variables are passed in the header of the routine !
      !  and overwrite previous HFODD definitions.                          !
      !---------------------------------------------------------------------!
      welfile=TRIM(FILTHO)
      file_exists=.False.; inquire(file=welfile, exist=file_exists); ierr=0
      If(file_exists) Then
         file_opened=.False.; inquire(unit=lwin, opened=file_opened)
         If(file_opened) Then
            Close(lwin)
         Else
         End If
         Open(lwin,file=welfile,status='old',form='unformatted',IOSTAT=ierr)
         If(ierr.NE.0) Then
            Write(6,'("The file ",a30," could not be opened!")') welfile
            Stop 'welfile could not be opened!'
         End If
      Else
         Stop 'welfile does not exist!'
      End If
      !---------------------------------------------------------------------
      ! open and read input file 'welfile' if exists
      !---------------------------------------------------------------------
      Read(lwin,Err=100,End=100)  npr11,npr1,ngh1,ngl1,n001,nb1,nt1
      Read(lwin,Err=100,End=100)  b0,bz,bp,beta0,si,etot,rms,bet,xmix,v0,v1,dalf,pwi,del,ept, &
                                  ala,ala2,alast,tz,varmas,varmasNZ,pjmassNZ,ass,skass
      Read(lwin,Err=100,End=100)  ntx,nb,nhhdim
      NLIMIT=ntx; nbx=nb; nhhdim2=2*nhhdim; nhhdim3=3*nhhdim; nhhdim4=4*nhhdim
      Read(lwin,Err=100,End=100) lambda_max; maxLambda=lambda_max
      !
      Allocate(xmultLag(maxLambda)); Allocate(id(nb),brin(nhhdim4+maxLambda))
      !
      Read(lwin,Err=100,End=100) xmultLag
      Read(lwin,Err=100,End=100) id
      Read(lwin,Err=100,End=100) brin
      !---------------------------------------------------------------------
      ! First batch of variable assignment
      !---------------------------------------------------------------------
      BETTHO=beta0; B0_THO=b0; BZ_THO=bz; BP_THO=bp
      ! pairing forces
      PRHO_N=v0(1); PRHODN=zero; PRHOSN=0.32d0; POWERN=one
      PRHO_P=v0(2); PRHODP=zero; PRHOSP=0.32d0; POWERP=one
      ! pairing properties
      EFERMN=ala(1); EFERMP=ala(2); EFER2N=ala2(1); EFER2P=ala2(2); DELTAN=del(1); DELTAP=del(2)
      ! Sizes
      ndx=MaxVal(id); nqx=ndx*ndx; nb2x=nbx+nbx
      !
      ! Array allocations
      If(Allocated(Omega2)) Deallocate(Omega2,Parity); Allocate(Omega2(ntx),Parity(ntx))
      If(Allocated(Nm_tho)) Deallocate(Nm_tho,nz_tho,nr_tho)
      Allocate(Nm_tho(ntx),nz_tho(ntx),nr_tho(ntx),ml_tho(ntx),Sigma2(ntx))
      If(Allocated(THOPHN)) Deallocate(THOPHN,THOPHP,THOPPN,THOPPP)
      Allocate(THOPHN(1:ntx,1:ntx)); THOPHN = zero; Allocate(THOPHP(1:ntx,1:ntx)); THOPHP = zero
      Allocate(THOPPN(1:ntx,1:ntx)); THOPPN = zero; Allocate(THOPPP(1:ntx,1:ntx)); THOPPP = zero
      !---------------------------------------------------------------------
      ! Quantum Number: {nr, nz, Lambda, Omega}
      !
      ! N = 2nr + Lambda + nz
      ! Omega = Lambda + Sigma
      !
      ! Restrictions: - Omega > 0 only
      !               - ms = +/- 1 (both spin projections)
      !---------------------------------------------------------------------
      ibasis=0; ibro=0; IndexStart = 1; NmxTHO=-1; counter = 0
      Do ib=1,NB
         ND=id(ib); i=ibro
         Do N1=1,ND
            ibasis=ibasis+1; IndexBra = ibasis
            ! read and assign quantum numbers for HO basis
            Read(lwin,Err=100,End=100) Omega2(ibasis), Parity(ibasis), Nm_tho(ibasis), &
                                       nz_tho(ibasis), ml_tho(ibasis)
            If (Nm_tho(ibasis).GT. NmxTHO) NmxTHO=Nm_tho(ibasis)
            nr_tho(ibasis) = (Nm_tho(ibasis) - nz_tho(ibasis) - ml_tho(ibasis))/2
            Sigma2(ibasis) =  Omega2(ibasis) - 2*ml_tho(ibasis)
            Do N2=1,N1
               IndexKet = IndexStart + N2 - 1
               i=i+1; n1n2 = n1+(n2-1)*nd; n2n1 = n2+(n1-1)*nd
               ! Assign mean-field and pairing fields for protons and neutrons
               THOPHN(IndexBra,IndexKet) = brin(i)
               THOPHP(IndexBra,IndexKet) = brin(i+nhhdim)
               THOPPN(IndexBra,IndexKet) = brin(i+nhhdim2)
               THOPPP(IndexBra,IndexKet) = brin(i+nhhdim3)
               ! Both fields are symmetric
               THOPHN(IndexKet,IndexBra) = brin(i)
               THOPHP(IndexKet,IndexBra) = brin(i+nhhdim)
               THOPPN(IndexKet,IndexBra) = brin(i+nhhdim2)
               THOPPP(IndexKet,IndexBra) = brin(i+nhhdim3)
            End Do !N2
         End Do !N1
         ibro=i
         !  Starting index of the current block
         IndexStart = IndexStart + ND
      End Do !IB
      !---------------------------------------------------------------------
      ! Test compatibility of sizes
      !---------------------------------------------------------------------
      NsizeHFBTHO = ibasis
      ! If the basis is oblate, numbers of quanta in the x- and y-direction in HFODD
      ! must be updated
      If(beta0 .Le. 0.0d0) N00THO=NmxTHO
      If(NsizeHFBTHO .Ne. ntx) THEN
         Write(6,'(''NsizeHFBTHO = '',i6,'' ntx = '',i6)') NsizeHFBTHO, ntx
         Stop 'Error (3) in HFBTHO_HFBODD'
      End If
      !---------------------------------------------------------------------
      ! Read and set blocking quantities
      !---------------------------------------------------------------------
      Read(lwin,Err=100,End=100)  bloall1
      Read(lwin,Err=100,End=100)  bloblo,blo123,blok1k2,blomax,bloqpdif
      If(npr11.Eq.2*(npr11/2).And.npr1.Eq.2*(npr1/2)) blomax=0
      If(bloall1.Ne.bloall) go to 100
      !---------------------------------------------------------------------
      ! Close file
      !---------------------------------------------------------------------
      Close(lwin)
      Go To 200
      !---------------------------------------------------------------------
      ! If file is corrupted, arrive here, and close file
      !---------------------------------------------------------------------
100   Continue
      Close(lwin)
      Write(6,'("The file ",a30," is corrupted, missing, or wrong!")') welfile
      Write(6,'("STARTING FROM SCRATCH WITH ININ=IABS(ININ)!")')
      If(Allocated(xmultLag)) Deallocate(xmultLag)
      If(Allocated(id))       Deallocate(id)
      If(Allocated(brin))     Deallocate(brin)
      If(Allocated(Omega2))   Deallocate(Omega2)
      If(Allocated(Parity))   Deallocate(Parity)
      Return
      !---------------------------------------------------------------------
      ! If file is OK, arrive here, deallocate and exit
      !---------------------------------------------------------------------
200   Continue
      Deallocate(id,brin); Deallocate(Omega2,Parity)
      !
   End Subroutine ReadHFBTHO
   !==========================================================================================================
   !
   !==========================================================================================================
   Subroutine FromHFBTHO(NDTRAN,NDSPIN)
      Integer :: NDTRAN,NDSPIN
      Integer(inik) :: N, nr, nr_hfbtho, nz, ml, ml_hfbtho, ms, i, AllocateStatus = 0
      Integer(inik) :: IndexBra_HFODD, IndexKet_HFODD, IndexBra_HFBTHO, IndexKet_HFBTHO
      Integer(inik) :: Nm_bra, nz_bra, nr_bra, la_bra, Nm_ket, nz_ket, nr_ket, la_ket, ms_bra, ms_ket
      Integer(inik) :: Indexx, Nperpe, counter
      !--------------------------------------------------------------------!
      !  This subroutine defines the HFB matrix in cylindrical coordinates !
      !  according to HFODD ordering convention, taking as input the HFB   !
      !  matrix read in ReadHFBTHO(). Since HFBTHO conserves time-reversal !
      !  symmetry while HFODD does not, cylindrical quantum numbers Lambda !
      !  and nrho can take different values in the two codes: States with  !
      !  -|Lambda| and +|Lambda| are degenerate in HFBTHO, and at restart  !
      !  they will also be degenerate in HFODD. Later in the iterations,   !
      !  time-reversal symmetry may be broken. The routine can handle      !
      !  deformed, or stretched, bases, provided the deformation is along  !
      !  the z-axis (only possible convention in HFBTHO).                  !
      !--------------------------------------------------------------------!
      !--------------------------------------------------------------------------
      ! Setting the correspondence between HFBTHO quantum numbers and global index
      !            N = 2*nr + nz + Lambda, Lambda = ml, Omega  = mj
      !               ms = 0 <=> SPIN DOWN ms = -1  (mj = ml + ms)
      !               ms = 1 <=> SPIN UP   ms = +1  (mj = ml + ms)
      !    HFB-THO ........: 0 <= nr <= Nperpe/2, 0 <= Lambda <= +Nperpe, step 1
      !    Quantum Number .: {N, nr, nz, Lambda}
      !    Restrictions ...: - All Omega
      !                   - NO spin projection at this time
      !--------------------------------------------------------------------------
      Allocate(index_HFBTHO(0:NmxTHO,0:NmxTHO,0:NmxTHO,0:NmxTHO,0:NDSPIN)); index_HFBTHO=0
      If(Allocated(validate_HFBTHO)) Deallocate(validate_HFBTHO)
      Allocate(validate_HFBTHO(0:NmxTHO,0:NmxTHO,0:NmxTHO,0:NmxTHO)); validate_HFBTHO=0
      Do i = 1, NsizeHFBTHO
         N  = Nm_tho(i)
         nz = nz_tho(i)
         nr = nr_tho(i)
         ml = ml_tho(i)
         ms = (Sigma2(i)+1)/2
         index_HFBTHO(N, nz, nr, ml, ms) = i
         validate_HFBTHO(N, nz, nr, ml) = 1
      End Do
      !--------------------------------------------------------------------------
      ! Setting the correspondence between HFODD quantum numbers and global index
      !            N = 2*nr + nz + Lambda, Lambda = ml, Omega  = mj
      !               ms = 0 <=> SPIN DOWN ms = -1  (mj = ml + ms)
      !               ms = 1 <=> SPIN UP   ms = +1  (mj = ml + ms)
      !    HFODD ..........: 0 <= nr <= Nperpe,  -Nperpe <= Lambda <= +Nperpe, step 2
      !    Quantum Number .: {N, nr, nz, Lambda}
      !    Restrictions ...: - All Omega
      !                      - NO spin projection at this time
      !--------------------------------------------------------------------------
      ! Computing the size of the HFODD matrices
      Indexx = 0
      Do N = 0, NmxTHO
         Do nz = 0, NmxTHO
            Nperpe = N - nz
            Do nr = 0, Nperpe
                                    nr_hfbtho = nr
               If(nr .GT. Nperpe/2) nr_hfbtho = Nperpe - nr
               ml_hfbtho = Abs(Nperpe - 2*nr)
               ! selecting only those states which have been effectively used in HFBTHO
               If(validate_HFBTHO(N, nz, nr_hfbtho, ml_hfbtho) .Eq. 1) Then
                  Indexx = Indexx + 1
               End If
            End Do
         End Do
      End Do
      NsizeHFODD = Indexx
      ! Allocating memory for HFODD vectors of cylindrical quantum numbers
      Allocate(Nm_odd(1:NsizeHFODD),nz_odd(1:NsizeHFODD),nr_odd(1:NsizeHFODD),ml_odd(1:NsizeHFODD))
      ! Defining cylindrical quantum numbers according to HFODD convention
      Indexx = 0
      Do N = 0, NmxTHO
         Do nz = 0, NmxTHO
            Nperpe = N - nz
            Do nr = 0, Nperpe
               ml = Nperpe - 2*nr
                                    nr_hfbtho = nr
               If(nr .Gt. Nperpe/2) nr_hfbtho = Nperpe - nr
               ml_hfbtho = Abs(ml)
               If(validate_HFBTHO(N, nz, nr_hfbtho, ml_hfbtho) .Eq. 1) Then
                  Indexx = Indexx + 1
                  Nm_odd(Indexx) = N
                  nz_odd(Indexx) = nz
                  nr_odd(Indexx) = nr
                  ml_odd(Indexx) = ml
               End If
            End Do
         End Do
      End Do
      !--------------------------------------------------------------------------
      ! Allocating memory for the HFB matrix (mean-field: HF_PHN, HF_PHP,
      ! pairing: HF_PPN, HF_PPP) coming from HFBTHO
      !--------------------------------------------------------------------------
      Allocate(HF_PHN(1:NDTRAN,1:NDTRAN,0:NDSPIN,0:NDSPIN),STAT=AllocateStatus)
      If(AllocateStatus /= 0) Then
          Write(6,'("Error in Allocating HF_PHN")')
          Stop 'Error in Allocating HF_PHN - Routine FromHFBTHO in hfodd_interface.f90'
      End If
      HF_PHN = Cmplx(zero,zero)
      Do ms_ket = 0, 1
         Do ms_bra = 0, 1
            !
            Do IndexKet_HFODD = 1, NsizeHFODD
               Nm_ket = Nm_odd(IndexKet_HFODD)
               nz_ket = nz_odd(IndexKet_HFODD)
               Nperpe = Nm_ket - nz_ket
               nr_ket = nr_odd(IndexKet_HFODD)
               la_ket = ml_odd(IndexKet_HFODD)
               ! reset nr_ket to HFBTHO convention
               If(nr_ket .GT. Nperpe/2) Then
                  nr_ket =  Nperpe - nr_ket
                  ms = 1 - ms_ket
               Else
                  ms = ms_ket
               End If
               If(la_ket .EQ. 0 .AND. ms .EQ. 0) Then
                  IndexKet_HFBTHO = index_HFBTHO(Nm_ket, nz_ket, nr_ket, Abs(la_ket), 1)
               Else
                  IndexKet_HFBTHO = index_HFBTHO(Nm_ket, nz_ket, nr_ket, Abs(la_ket), ms)
               End If
               !
               Do IndexBra_HFODD = 1, NsizeHFODD
                  Nm_bra = Nm_odd(IndexBra_HFODD)
                  nz_bra = nz_odd(IndexBra_HFODD)
                  Nperpe = Nm_bra - nz_bra
                  nr_bra = nr_odd(IndexBra_HFODD)
                  la_bra = ml_odd(IndexBra_HFODD)
                  ! reset nr_bra to HFBTHO convention
                  If(nr_bra .GT. Nperpe/2) Then
                     nr_bra =  Nperpe - nr_bra
                     ms = 1 - ms_bra
                  Else
                     ms = ms_bra
                  End If
                  If(la_bra .EQ. 0 .AND. ms .EQ. 0) Then
                     IndexBra_HFBTHO = index_HFBTHO(Nm_bra, nz_bra, nr_bra, Abs(la_bra), 1)
                  Else
                     IndexBra_HFBTHO = index_HFBTHO(Nm_bra, nz_bra, nr_bra, Abs(la_bra), ms)
                  End If
                  !
                  If(IndexBra_HFBTHO .GT. 0 .And. IndexKet_HFBTHO .GT. 0) Then
                     HF_PHN(IndexBra_HFODD, IndexKet_HFODD, ms_bra, ms_ket) = Cmplx(THOPHN(IndexBra_HFBTHO, IndexKet_HFBTHO),0.0_dnik)
                 End If
                  !
               End Do
            End Do
            !
         End Do
      End Do
      Deallocate(THOPHN)
      !--------------------------------------------------------------------------
      ! Transfering the matrix elements of the HFB matrix from HFBTHO to HFODD
      ! This implies: - different ordering of quantum numbers, i.e. different
      !                 index <-> q.n. convention
      !               - different cylindrical quantum numbers: Lambda >= 0 in
      !                 HFBTHO
      !--------------------------------------------------------------------------
      Allocate(HF_PHP(1:NDTRAN,1:NDTRAN,0:NDSPIN,0:NDSPIN),STAT=AllocateStatus)
      If(AllocateStatus /= 0) Then
          Write(6,'("Error in Allocating HF_PHP")')
          Stop 'Error in Allocating HF_PHP - Routine FromHFBTHO in hfodd_interface.f90'
      End If
      HF_PHP = Cmplx(zero,zero)
      Do ms_ket = 0, 1
         Do ms_bra = 0, 1
            !
            Do IndexKet_HFODD = 1, NsizeHFODD
               Nm_ket = Nm_odd(IndexKet_HFODD)
               nz_ket = nz_odd(IndexKet_HFODD)
               Nperpe = Nm_ket - nz_ket
               nr_ket = nr_odd(IndexKet_HFODD)
               la_ket = ml_odd(IndexKet_HFODD)
               ! reset nr_ket to HFBTHO convention
               If(nr_ket .GT. Nperpe/2) Then
                  nr_ket =  Nperpe - nr_ket
                  ms = 1 - ms_ket
               Else
                  ms = ms_ket
               End If
               If(la_ket .EQ. 0 .AND. ms .EQ. 0) Then
                  IndexKet_HFBTHO = index_HFBTHO(Nm_ket, nz_ket, nr_ket, Abs(la_ket), 1)
               Else
                  IndexKet_HFBTHO = index_HFBTHO(Nm_ket, nz_ket, nr_ket, Abs(la_ket), ms)
               End If
               !
               Do IndexBra_HFODD = 1, NsizeHFODD
                  Nm_bra = Nm_odd(IndexBra_HFODD)
                  nz_bra = nz_odd(IndexBra_HFODD)
                  Nperpe = Nm_bra - nz_bra
                  nr_bra = nr_odd(IndexBra_HFODD)
                  la_bra = ml_odd(IndexBra_HFODD)
                  ! reset nr_bra to HFBTHO convention
                  If(nr_bra .GT. Nperpe/2) Then
                     nr_bra =  Nperpe - nr_bra
                     ms = 1 - ms_bra
                  Else
                     ms = ms_bra
                  End If
                  If(la_bra .EQ. 0 .AND. ms .EQ. 0) Then
                     IndexBra_HFBTHO = index_HFBTHO(Nm_bra, nz_bra, nr_bra, Abs(la_bra), 1)
                  Else
                     IndexBra_HFBTHO = index_HFBTHO(Nm_bra, nz_bra, nr_bra, Abs(la_bra), ms)
                  End If
                  !
                  If(IndexBra_HFBTHO .GT. 0 .And. IndexKet_HFBTHO .GT. 0) Then
                     HF_PHP(IndexBra_HFODD, IndexKet_HFODD, ms_bra, ms_ket) = Cmplx(THOPHP(IndexBra_HFBTHO, IndexKet_HFBTHO),0.0_dnik)
                 End If
                  !
               End Do
            End Do
            !
         End Do
      End Do
      Deallocate(THOPHP)
      !--------------------------------------------------------------------------
      ! Transfering the matrix elements of the HFB matrix from HFBTHO to HFODD
      ! This implies: - different ordering of quantum numbers, i.e. different
      !                 index <-> q.n. convention
      !               - different cylindrical quantum numbers: Lambda >= 0 in
      !                 HFBTHO
      !--------------------------------------------------------------------------
      Allocate(HF_PPN(1:NDTRAN,1:NDTRAN,0:NDSPIN,0:NDSPIN),STAT=AllocateStatus)
      If(AllocateStatus /= 0) Then
          Write(6,'("Error in Allocating HF_PPN")')
          Stop 'Error in Allocating HF_PPN - Routine FromHFBTHO in hfodd_interface.f90'
      End If
      HF_PPN = Cmplx(zero,zero)
      Do ms_ket = 0, 1
         Do ms_bra = 0, 1
            !
            Do IndexKet_HFODD = 1, NsizeHFODD
               Nm_ket = Nm_odd(IndexKet_HFODD)
               nz_ket = nz_odd(IndexKet_HFODD)
               Nperpe = Nm_ket - nz_ket
               nr_ket = nr_odd(IndexKet_HFODD)
               la_ket = ml_odd(IndexKet_HFODD)
               ! reset nr_ket to HFBTHO convention
               If(nr_ket .GT. Nperpe/2) Then
                  nr_ket =  Nperpe - nr_ket
                  ms = 1 - ms_ket
               Else
                  ms = ms_ket
               End If
               If(la_ket .EQ. 0 .AND. ms .EQ. 0) Then
                  IndexKet_HFBTHO = index_HFBTHO(Nm_ket, nz_ket, nr_ket, Abs(la_ket), 1)
               Else
                  IndexKet_HFBTHO = index_HFBTHO(Nm_ket, nz_ket, nr_ket, Abs(la_ket), ms)
               End If
               !
               Do IndexBra_HFODD = 1, NsizeHFODD
                  Nm_bra = Nm_odd(IndexBra_HFODD)
                  nz_bra = nz_odd(IndexBra_HFODD)
                  Nperpe = Nm_bra - nz_bra
                  nr_bra = nr_odd(IndexBra_HFODD)
                  la_bra = ml_odd(IndexBra_HFODD)
                  ! reset nr_bra to HFBTHO convention
                  If(nr_bra .GT. Nperpe/2) Then
                     nr_bra =  Nperpe - nr_bra
                     ms = 1 - ms_bra
                  Else
                     ms = ms_bra
                  End If
                  If(la_bra .EQ. 0 .AND. ms .EQ. 0) Then
                     IndexBra_HFBTHO = index_HFBTHO(Nm_bra, nz_bra, nr_bra, Abs(la_bra), 1)
                  Else
                     IndexBra_HFBTHO = index_HFBTHO(Nm_bra, nz_bra, nr_bra, Abs(la_bra), ms)
                  End If
                  !
                  If(IndexBra_HFBTHO .GT. 0 .And. IndexKet_HFBTHO .GT. 0) Then
                     HF_PPN(IndexBra_HFODD, IndexKet_HFODD, ms_bra, ms_ket) = Cmplx(THOPPN(IndexBra_HFBTHO, IndexKet_HFBTHO),0.0_dnik)
                 End If
                  !
               End Do
            End Do
            !
         End Do
      End Do
      Deallocate(THOPPN)
      !--------------------------------------------------------------------------
      ! Transfering the matrix elements of the HFB matrix from HFBTHO to HFODD
      ! This implies: - different ordering of quantum numbers, i.e. different
      !                 index <-> q.n. convention
      !               - different cylindrical quantum numbers: Lambda >= 0 in
      !                 HFBTHO
      !--------------------------------------------------------------------------
      Allocate(HF_PPP(1:NDTRAN,1:NDTRAN,0:NDSPIN,0:NDSPIN),STAT=AllocateStatus)
      If(AllocateStatus /= 0) Then
          Write(6,'("Error in Allocating HF_PPP")')
          Stop 'Error in Allocating HF_PPP - Routine FromHFBTHO in hfodd_interface.f90'
      End If
      HF_PPP = Cmplx(zero,zero)
      Do ms_ket = 0, 1
         Do ms_bra = 0, 1
            !
            Do IndexKet_HFODD = 1, NsizeHFODD
               Nm_ket = Nm_odd(IndexKet_HFODD)
               nz_ket = nz_odd(IndexKet_HFODD)
               Nperpe = Nm_ket - nz_ket
               nr_ket = nr_odd(IndexKet_HFODD)
               la_ket = ml_odd(IndexKet_HFODD)
               ! reset nr_ket to HFBTHO convention
               If(nr_ket .GT. Nperpe/2) Then
                  nr_ket =  Nperpe - nr_ket
                  ms = 1 - ms_ket
               Else
                  ms = ms_ket
               End If
               If(la_ket .EQ. 0 .AND. ms .EQ. 0) Then
                  IndexKet_HFBTHO = index_HFBTHO(Nm_ket, nz_ket, nr_ket, Abs(la_ket), 1)
               Else
                  IndexKet_HFBTHO = index_HFBTHO(Nm_ket, nz_ket, nr_ket, Abs(la_ket), ms)
               End If
               !
               Do IndexBra_HFODD = 1, NsizeHFODD
                  Nm_bra = Nm_odd(IndexBra_HFODD)
                  nz_bra = nz_odd(IndexBra_HFODD)
                  Nperpe = Nm_bra - nz_bra
                  nr_bra = nr_odd(IndexBra_HFODD)
                  la_bra = ml_odd(IndexBra_HFODD)
                  ! reset nr_bra to HFBTHO convention
                  If(nr_bra .GT. Nperpe/2) Then
                     nr_bra =  Nperpe - nr_bra
                     ms = 1 - ms_bra
                  Else
                     ms = ms_bra
                  End If
                  If(la_bra .EQ. 0 .AND. ms .EQ. 0) Then
                     IndexBra_HFBTHO = index_HFBTHO(Nm_bra, nz_bra, nr_bra, Abs(la_bra), 1)
                  Else
                     IndexBra_HFBTHO = index_HFBTHO(Nm_bra, nz_bra, nr_bra, Abs(la_bra), ms)
                  End If
                  !
                  If(IndexBra_HFBTHO .GT. 0 .And. IndexKet_HFBTHO .GT. 0) Then
                     HF_PPP(IndexBra_HFODD, IndexKet_HFODD, ms_bra, ms_ket) = Cmplx(THOPPP(IndexBra_HFBTHO, IndexKet_HFBTHO),0.0_dnik)
                 End If
                  !
               End Do
            End Do
            !
         End Do
      End Do
      !--------------------------------------------------------------------------
      ! Deallocations
      !--------------------------------------------------------------------------
      Deallocate(index_HFBTHO)
      Deallocate(THOPPP)
      Deallocate(Nm_odd,nr_odd,nz_odd,ml_odd)
      Deallocate(Nm_tho,nr_tho,nz_tho,ml_tho,Sigma2)
   End Subroutine FromHFBTHO
   !==========================================================================================================
   !
   !==========================================================================================================
   Subroutine Basis_THO(N_UNIT,NDTRAN,NDSPIN,NOSCIL,BP_THO,BZ_THO,HPERPE,HPARAL,Parity,EBASECUT,hbzero)
      ! Inputs with type defined in HFODD
      Logical :: Parity
      Integer :: N_UNIT,NDTRAN,NDSPIN,NOSCIL
      Real(dnik) :: BP_THO,BZ_THO,HPERPE,HPARAL,EBASECUT,hbzero
      ! Local with type defined in module (in practice, it will be the
      ! same Real(dnik))
      Real(dnik) :: hbz,hbp,bz,bp,ee
      Integer(inik) :: nze,nre,ib,ilauf,jlauf,k,la,le,ip,ir,iz,il,is,nn,nd,IBX,n1,n2,iw, &
                       n00,n00max,ierror_flag,nlm,nzm,nrn,nqp,nuv,counter
      !--------------------------------------------------------------------!
      ! This subroutine defines the basis needed by HFBTHO based on input  !
      ! data from HFODD. It is needed to perform the inverse tranformation !
      ! HFODD -> HFBTHO and enable restarts with HFBTHO from HFODD record  !
      ! files.
      !--------------------------------------------------------------------!
      n00=NOSCIL; n00max=50; bp=BP_THO; bz=BZ_THO
      !-------------------------------------------------------------------------------
      ! Determine total number of K-blocks nbx, actual total number of states ntx
      !-------------------------------------------------------------------------------
      !hbz=two*hbzero/bz**2; hbp=two*hbzero/bp**2;
      hbz=HPARAL; hbp=HPERPE
      nze=n00; nre=n00max/2; ib=0; ilauf=0; nzm=0; nlm=0; nqp=0; nuv=0
      nom_tho=0; nzx_tho=0; nrx_tho=0; nlx_tho=0; nnm_tho=0
      ! loop over k-quantum number
      Do k=1,n00max+1
         la=k-1; le=min0(n00max,k)
         ! loop over parity
         If(.Not.Parity) jlauf=ilauf !Nop
         Do ip=1,2
            If(Parity) jlauf=ilauf !Yesp
            Do ir=0,nre
               Do iz=0,nze
                  Do il=la,le
                     Do is=+1,-1,-2
                        If (iz+2*ir+il.Gt.n00max)    Cycle
                        If (il+(is+1)/2.Ne.k)        Cycle
                        If (Mod(iz+il,2).Ne.ip-1)    Cycle
                        ee=hbz*(Real(iz)+half)+hbp*(two*Real(ir)+Real(il)+one)
                        If(ee.Lt.EBASECUT) Then
                           ilauf=ilauf+1
                           nzx_tho=Max(nzx_tho,iz);    nrx_tho=Max(nrx_tho,ir); nlx_tho=Max(nlx_tho,il)
                           nom_tho=Max(nom_tho,2*k-1); nnm_tho=Max(nnm_tho,iz+2*ir+il)
                        End If
                     End Do
                  End Do
               End Do
            End Do
            If(Parity) Then                !Yesp
               If (ilauf.Gt.jlauf) Then
                  ib=ib+1
                  nd=ilauf-jlauf
                  ndx_tho=Max(ndx_tho,nd)
                  nqp=nqp+nd; nuv=nuv+nd*nd
               End If
            End If
         End Do
         If(.Not.Parity) Then              !Nop
            If(ilauf.Gt.jlauf) Then
               ib=ib+1
               nd=ilauf-jlauf
               ndx_tho=Max(ndx_tho,nd)
               nqp=nqp+nd; nuv=nuv+nd*nd
            End If
         End If
      End Do
      nbx_tho=ib; ntx_tho=ilauf
      !-----------------------------------------------
      ! loop over k-quantum number
      !-----------------------------------------------
      Allocate(nz_tho(ntx_tho),nr_tho(ntx_tho),nl_tho(ntx_tho),ns_tho(ntx_tho),id_tho(nbx_tho),Nm_tho(ntx_tho))
      Allocate(npar_tho(ntx_tho),ia_tho(nbx_tho),ikb_tho(nbx_tho),ipb_tho(nbx_tho))
      ierror_flag=0; ib=0; ilauf=0
      Do k=1,n00max+1
         la=k-1; le=min0(n00max,k)
         ! loop over parity
         If(.Not.Parity) jlauf=ilauf !Nop
         Do ip=1,2
            If(Parity) jlauf=ilauf   !Yesp
            Do ir=0,nre
               Do iz=0,nze
                  Do il=la,le
                     Do is=+1,-1,-2
                        If(iz+2*ir+il.Gt.n00max) Cycle
                        If(il+(is+1)/2.Ne.k)     Cycle
                        If(Mod(iz+il,2).Ne.ip-1) Cycle
                        ee=hbz*(Real(iz)+half)+hbp*(two*Real(ir)+Real(il)+one)
                        If(ee.Lt.EBASECUT) Then
                           ilauf=ilauf+1
                           If (ilauf.Gt.ntx_tho) Then
                              ierror_flag=ierror_flag+1
                              !ierror_info(ierror_flag)='STOP: in base: ntx too small'
                              STOP 'in base: ntx too small'
                              Return
                           End If
                           nz_tho(ilauf)=iz; nr_tho(ilauf)=ir; nl_tho(ilauf)=il
                           ns_tho(ilauf)=is; npar_tho(ilauf)=ip
                           nn =iz+2*ir+il; Nm_tho(ilauf)=nn
                        End If
                     End Do
                  End Do
               End Do
            End Do
            !-----------------------------------------------
            ! Block memory
            !-----------------------------------------------
            If(Parity) Then                !Yesp
               If (ilauf.Gt.jlauf) Then
                  ib=ib+1
                  ia_tho(ib)=jlauf; id_tho(ib)=ilauf-jlauf
                  ikb_tho(ib)=k; ipb_tho(ib)=ip
               End If
               If(id_tho(ib).Eq.0) Then
                  ierror_flag=ierror_flag+1
                  !ierror_info(ierror_flag)='STOP: in base Block Memory(1)'
                  Stop 'in base Block Memory(1)'
                  Return
               End If
            End If
         End Do ! end of ip
         !-----------------------------------------------
         ! Block memory
         !-----------------------------------------------
         If(.Not.Parity) Then               !Nop
            If (ilauf.Gt.jlauf) Then
               ib=ib+1
               ia_tho(ib)=jlauf; id_tho(ib)=ilauf-jlauf
               nn = nz_tho(ilauf)+2*nr_tho(ilauf)+nl_tho(ilauf); ip = 2 - Mod(nn,2)
               ikb_tho(ib)=k; ipb_tho(ib)=ip
            End If
            If(id_tho(ib).Eq.0) Then
               ierror_flag=ierror_flag+1
               !ierror_info(ierror_flag)='STOP: in base Block Memory(2)'
               Return
            End If
         End If
      End Do ! end k
      nb_tho=ib;  nt_tho=ilauf
      !-----------------------------------------------
      ! broyden/linear mixing (storage)
      !-----------------------------------------------
      nhhdim_tho=0
      Do ib=1,nb_tho
         nd=id_tho(ib)
         Do n1=1,nd
            Do n2=1,n1
               nhhdim_tho=nhhdim_tho+1
            End Do
         End Do
      End Do
      nhhdim2_tho=2*nhhdim_tho; nhhdim3_tho=3*nhhdim_tho; nhhdim4_tho=4*nhhdim_tho
      !-----------------------------------------------
      ! print statistics
      !-----------------------------------------------
      Write(N_UNIT,'(a)')    '  Comparison with bookkeeping spherical basis:'
      Write(N_UNIT,'(a,i4)') '  maximal number of shells n00 ...: ',n00
      Write(N_UNIT,'(a,i4)') '  number of blocks: nb ...........: ',nb_tho
      Write(N_UNIT,'(a,i4)') '  number of levels: nt ...........: ',nt_tho
      Write(N_UNIT,'(a,i4)') '  maximal 2*omega : nom ..........: ',nom_tho
      Write(N_UNIT,'(a,i4)') '  maximal nz:       nzm ..........: ',nzx_tho
      Write(N_UNIT,'(a,i4)') '  maximal nr:       nrm ..........: ',nrx_tho
      Write(N_UNIT,'(a,i4)') '  maximal ml:       nlm ..........: ',nlx_tho
      Write(N_UNIT,'(a,i4)') '  2 x biggest block dim. .........: ',ndx_tho
      Write(N_UNIT,'(a,i6)') '  Nonzero elements of HH .........: ',nhhdim_tho
      Write(N_UNIT,'(a,i6)') '  Number Broyden elements ........: ',nhhdim4_tho
      Write(N_UNIT,'(a)')
      If(nzx_tho.Ge.n00max.Or.(nom_tho-1)/2.Eq.n00max) Then
         Write(*,*) 'nzm=',nzx_tho,'  (nom-1)/2=',(nom_tho-1)/2,'  n00max=',n00max
         ierror_flag=ierror_flag+1
         !ierror_info(ierror_flag)='STOP: Please increase n00max to have correct basis'
      End If
      Deallocate(npar_tho,ia_tho,ikb_tho,ipb_tho)
      !
   End Subroutine Basis_THO
   !==========================================================================================================
   !
   !==========================================================================================================
   Subroutine ToHFBTHO(NDTRAN,NDSPIN,NOSCIL,LDBASE,NOFCNS)
      ! Inputs with type defined in HFODD
      Integer :: NDTRAN,NDSPIN,NOSCIL,LDBASE,NOFCNS
      ! Local with type defined in module (in practice, it will be the
      ! same Real(dnik))
      Integer(inik) :: N, nr, nr_hfbtho, nz, ml, ml_hfbtho, ms, i, AllocateStatus = 0
      Integer(inik) :: ib,n1,n2,nd,ibro,ibasis,indexStart
      Integer(inik) :: IndexBra_HFODD, IndexKet_HFODD, IndexBra_HFBTHO, IndexKet_HFBTHO
      Integer(inik) :: Nm_bra, nz_bra, nr_bra, la_bra, Nm_ket, nz_ket, nr_ket, la_ket, ms_bra, ms_ket
      Integer(inik) :: Indexx, Nperpe,counter
      !--------------------------------------------------------------------!
      !  This subroutine defines the HFB matrix in cylindrical coordinates !
      !  according to HFBTHO ordering convention, taking as input the HFB  !
      !  matrix read in ReadHFBTHO(). Since HFBTHO conserves time-reversal !
      !  symmetry while HFODD does not, cylindrical quantum numbers Lambda !
      !  and nrho can take different values in the two codes: States with  !
      !  -|Lambda| and +|Lambda| are degenerate in HFBTHO, and at restart  !
      !  they will also be degenerate in HFODD. Later in the iterations,   !
      !  time-reversal symmetry may be broken. The routine can handle      !
      !  deformed, or stretched, bases, provided the deformation is along  !
      !  the z-axis (only possible convention in HFBTHO).                  !
      !--------------------------------------------------------------------!
      NsizeHFODD=LDBASE
      !--------------------------------------------------------------------------
      ! Setting the correspondence between HFODD quantum numbers and global index
      !            N = 2*nr + nz + Lambda, Lambda = ml, Omega  = mj
      !               ms = 0 <=> SPIN DOWN ms = -1  (mj = ml + ms)
      !               ms = 1 <=> SPIN UP   ms = +1  (mj = ml + ms)
      !    HFODD ..........: 0 <= nr <= Nperpe,  -Nperpe <= Lambda <= +Nperpe, step 2
      !    Quantum Number .: {N, nr, nz, Lambda}
      !    Restrictions ...: - All Omega
      !                      - NO spin projection at this time
      !--------------------------------------------------------------------------
      ! Allocating memory for HFODD vectors of cylindrical quantum numbers
      If(.NOT. Allocated(index_HFODD)) Allocate(index_HFODD(0:nnm_tho,0:nnm_tho,0:nnm_tho,0:nnm_tho,0:1))
      index_HFODD = 0
      ! Defining cylindrical quantum numbers according to HFODD convention
      Do ms_ket = 0, 1
         Indexx = 0
         Do N = 0, nnm_tho
            Do nz = 0, nnm_tho
               Nperpe = N - nz
               Do nr_ket = 0, Nperpe
                  !
                  If(nr_ket .GT. Nperpe/2) Then
                     nr = Nperpe - nr_ket
                     ms = ms_ket
                  Else
                     nr = nr_ket
                     ms = ms_ket
                  End If
                  ml = Nperpe - 2*nr
                  !
                  If(validate_HFBTHO(N, nz, nr, iabs(ml)).Eq.1) Then
                     Indexx = Indexx + 1
                     If(ml .EQ. 0 .AND. ms .EQ. 0) Then
                        index_HFODD(N, nz, nr, iabs(ml), 1) = Indexx
                     Else
                        index_HFODD(N, nz, nr, iabs(ml), ms) = Indexx
                     End If
                  End If
                  !
               End Do
            End Do
         End Do
      End Do
      !--------------------------------------------------------------------------
      ! Allocating memory for the HFB matrix (mean-field: THOPHN, THOPHP,
      ! pairing: THOPPN, THOPPP) coming from HFBTHO
      !--------------------------------------------------------------------------
      Allocate(THOPHN(1:ntx_tho,1:ntx_tho)); THOPHN=zero
      Allocate(THOPHP(1:ntx_tho,1:ntx_tho)); THOPHP=zero
      Allocate(THOPPN(1:ntx_tho,1:ntx_tho)); THOPPN=zero
      Allocate(THOPPP(1:ntx_tho,1:ntx_tho)); THOPPP=zero
      !--------------------------------------------------------------------------
      ! Transfering the matrix elements of the HFB matrix from HFODD to HFBTHO
      ! This implies: - different ordering of quantum numbers, i.e. different
      !                 index <-> q.n. convention
      !               - different cylindrical quantum numbers: Lambda >= 0 in
      !                 HFBTHO
      !--------------------------------------------------------------------------
      Do IndexKet_HFBTHO = 1, NsizeHFBTHO
         !
         Nm_ket = Nm_tho(IndexKet_HFBTHO)
         nz_ket = nz_tho(IndexKet_HFBTHO)
         nr_ket = nr_tho(IndexKet_HFBTHO)
         la_ket = nl_tho(IndexKet_HFBTHO)
         ms_ket = (1+ns_tho(IndexKet_HFBTHO))/2
         !
         IndexKet_HFODD = index_HFODD(Nm_ket, nz_ket, nr_ket, la_ket, ms_ket)
         !
         Do IndexBra_HFBTHO = 1, NsizeHFBTHO
            !
            Nm_bra = Nm_tho(IndexBra_HFBTHO)
            nz_bra = nz_tho(IndexBra_HFBTHO)
            nr_bra = nr_tho(IndexBra_HFBTHO)
            la_bra = nl_tho(IndexBra_HFBTHO)
            ms_bra = (1+ns_tho(IndexBra_HFBTHO))/2
            !
            IndexBra_HFODD = index_HFODD(Nm_bra, nz_bra, nr_bra, la_bra, ms_bra)
            !
            If(IndexBra_HFODD .GT. 0 .And. IndexKet_HFODD .GT. 0) Then
               THOPHN(IndexBra_HFBTHO, IndexKet_HFBTHO) = Real(HF_PHN(IndexBra_HFODD, IndexKet_HFODD, ms_bra, ms_ket))
               THOPHP(IndexBra_HFBTHO, IndexKet_HFBTHO) = Real(HF_PHP(IndexBra_HFODD, IndexKet_HFODD, ms_bra, ms_ket))
               THOPPN(IndexBra_HFBTHO, IndexKet_HFBTHO) = Real(HF_PPN(IndexBra_HFODD, IndexKet_HFODD, ms_bra, ms_ket))
               THOPPP(IndexBra_HFBTHO, IndexKet_HFBTHO) = Real(HF_PPP(IndexBra_HFODD, IndexKet_HFODD, ms_bra, ms_ket))
            End If
            !
         End Do
      End Do
      Deallocate(HF_PHN,HF_PHP,HF_PPN,HF_PPP)
      !--------------------------------------------------------------------------
      ! Fill out brin_tho vector
      !--------------------------------------------------------------------------
      If(.NOT. Allocated(brin_tho)) Allocate(brin_tho(1:nhhdim4_tho+8))
      ibasis=0; ibro=0; IndexStart = 1
      Do ib=1,nb_tho
         nd=id_tho(ib); i=ibro
         Do n1=1,nd
            ibasis=ibasis+1; IndexBra_HFBTHO = ibasis
            Do n2=1,n1
               IndexKet_HFBTHO = IndexStart + n2 - 1; i=i+1
               ! Assign mean-field and pairing fields for protons and neutrons
               brin_tho(i)             = THOPHN(IndexBra_HFBTHO,IndexKet_HFBTHO)
               brin_tho(i+nhhdim_tho)  = THOPHP(IndexBra_HFBTHO,IndexKet_HFBTHO)
               brin_tho(i+nhhdim2_tho) = THOPPN(IndexBra_HFBTHO,IndexKet_HFBTHO)
               brin_tho(i+nhhdim3_tho) = THOPPP(IndexBra_HFBTHO,IndexKet_HFBTHO)
            End Do !N2
         End Do !N1
         ibro=i
         IndexStart = IndexStart + nd
      End Do !IB
      Deallocate(THOPHN,THOPHP,THOPPN,THOPPP)
      !--------------------------------------------------------------------------
      ! Deallocations
      !--------------------------------------------------------------------------
      Deallocate(index_HFODD)
      ! Exit
   End Subroutine ToHFBTHO
   !=======================================================================
   !
   !=======================================================================
   Subroutine WriteHFBTHO(unit_tho,file_tho,n_tho,z_tho,NOSCIL,REFERN,REFERP, &
                         REDELN,REDELP,REFE2N,REFE2P,PRHO_N,PRHO_P,BZ_THO,BP_THO, &
                         beta0_THO,multLag)
      Implicit None
      ! From HFODD
      Integer :: unit_tho,z_tho,n_tho,NOSCIL
      Real(dnik) :: REFERN,REFERP,REDELN,REDELP,REFE2N,REFE2P,PRHO_N,PRHO_P,BZ_THO,BP_THO,beta0_THO
      Real(dnik) :: multLag(8)
      Character(12) :: file_tho
      ! Local
      Integer(inik)  :: l,n1,n2,nd,ib,nza,nra,nla,nsa,nlansa,ibasis,lambdaMax=8
      Integer(inik)  :: bloall; Parameter(bloall=200)
      Integer(inik)  :: bloblo(0:bloall,2),blo123(0:bloall,2),blok1k2(0:bloall,2),blomax(2)
      Real(dnik)     :: bloqpdif(0:bloall,2)
      Integer(inik)  :: npr1,npr11,ngh1,ngl1,n001,nt1,ntx1,nb1,nhhdim1
      Real(dnik)     :: b0,bz,bp,beta0,si,etot,rms(3),bet,xmix,CpV0(2),CpV1(2),pwi, &
                        del(2),ept(3),ala(2),ala2(2),alast(2),tz(2),varmas,   &
                        varmasNZ(2),pjmassNZ(2),ass(2),skass
      !---------------------------------------------------------------------
      ! Open 'FILTHO' file
      !---------------------------------------------------------------------
      Open(unit_tho,file=file_tho,status='unknown',form='unformatted')
      !---------------------------------------------------------------------
      ! Set everything to default value
      !---------------------------------------------------------------------
      npr11=n_tho; npr1=z_tho; ngh1=0; ngl1=0; n001=NOSCIL
      b0=zero; bz=BZ_THO; bp=BP_THO; beta0=beta0_THO; beta0=beta0_THO; si=one; etot=zero; rms=zero; bet=zero; xmix=0.3
      CpV0(1)=PRHO_N; CpV0(2)=PRHO_N; CpV1=(/0.5, 0.5/)
      pwi=2.0; del(1)=REDELN; del(2)=REDELP; ept=zero
      ala(1)=REFERN; ala(2)=REFERP; ala2(1)=REFE2N; ala2(2)=REFE2P; alast(1)=REFERN; alast(2)=REFERP
      tz=zero; varmas=zero;varmasNZ=zero; pjmassNZ=zero; ass=zero; skass=zero
      Do l=1,8
         brin_tho(nhhdim4_tho+l)=multLag(l)
      End Do
      !---------------------------------------------------------------------
      ! Write matrix elements to 'FILTHO' file
      !---------------------------------------------------------------------
      Write(unit_tho) npr11,npr1,ngh1,ngl1,n001,nbx_tho,ntx_tho
      Write(unit_tho) b0,bz,bp,beta0,si,etot,rms,bet,xmix,CpV0,CpV1,pwi,del,ept,  &
                      ala,ala2,alast,tz,varmas,varmasNZ,pjmassNZ,ass,skass
      Write(unit_tho) ntx_tho,nb_tho,nhhdim_tho
      Write(unit_tho) lambdaMax
      Write(unit_tho) multLag
      Write(unit_tho) id_tho
      Write(unit_tho) brin_tho
      ibasis=0
      Do ib=1,nb_tho
         nd=id_tho(ib)
         Do n1=1,nd
            ibasis=ibasis+1
            nla=nl_tho(ibasis); nra=nr_tho(ibasis); nza=nz_tho(ibasis); nsa=ns_tho(ibasis)
            nlansa=(-1)**(nza+nla)
            Write(unit_tho) 2*nla+nsa,nlansa,nza+2*nra+nla,nza,nla
         End Do
      End Do
      !---------------------------------------------------------------------
      ! blocking: sort blocking candidates first
      !---------------------------------------------------------------------
      bloblo=0; blo123=0; blok1k2=0; bloqpdif=0.0_dnik; blomax=0
      Write(unit_tho) bloall
      Write(unit_tho) bloblo,blo123,blok1k2,blomax,bloqpdif
      !---------------------------------------------------------------------
      ! Close 'FILTHO' file
      !---------------------------------------------------------------------
      Close(unit_tho)
      !--------------------------------------------------------------------------
      ! Deallocations
      !--------------------------------------------------------------------------
      Deallocate(id_tho,brin_tho)
      Deallocate(Nm_tho,nz_tho,nr_tho,nl_tho,ns_tho)
   End Subroutine WriteHFBTHO

 End Module hfodd_interface
