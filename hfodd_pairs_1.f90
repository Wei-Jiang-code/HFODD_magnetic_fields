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
 !   Module defining various derived types related to pairs, lists of   !
 !   pairs and lists of lists of pairs, as well as the routines needed  !
 !   to manipulate these objects.                                       !
 !                                                                      !
 !----------------------------------------------------------------------!

 Module pair_definition

     Implicit None

     Integer, Parameter, Public :: iprec=Kind(1), prec=Kind(1.0) ! precision

     ! New type for pair of indexes
     Type :: Pair
         Integer(iprec) :: mu
         Integer(iprec) :: nu
     End Type Pair

     ! New type for list of pairs
     Type :: PairSet
         Type(Pair), Dimension(:), pointer :: pair_vector
     End Type PairSet

     ! New type for list of set of pairs
     Type :: ListePairSet
         Integer, Dimension(:), pointer :: size_tiling
         Type(PairSet), Dimension(:), pointer :: liste
     End Type ListePairSet

 Contains

     !---------------------------------------------------------------------!
     !  This subroutine sorts a vector EL(1:N) by ascending order and      !
     !  returns the vector of indices in KL(1:N)                           !
     !---------------------------------------------------------------------!
     Subroutine sort_vector(EL,KL,N)

       Integer(iprec), INTENT(IN) :: N
       Integer(iprec), Dimension(1:N), INTENT(INOUT) :: KL
       Real(prec), Dimension(1:N), INTENT(INOUT) :: EL

       Integer(iprec) :: i,j,IEQUAL,IFILLL
       Integer(iprec), Allocatable :: KM(:),KWHER(:),KFROM(:)
       Real(prec), Allocatable :: EM(:)

       Allocate(KM(1:N),KWHER(1:N),KFROM(1:N))
       Allocate(EM(1:N))

       ! Rewriting input arrays to auxiliary ones
       KWHER(1:N)=1
       EM(1:N)=EL(1:N)
       KM(1:N)=KL(1:N)

       ! How many elements are smaller than the i-th one? In other words
       ! where should it go with respect to the offset position 1?
       Do J=1,N
          Do I=1,N
             If(EM(j).Lt.EM(i)) KWHER(i)=KWHER(i)+1
          End Do
       End Do

       ! Now we only have to resolve the case of identical elements
       ! From where the i-th element should come?
       KFROM(1:N)=0
       Do i=1,N
          KFROM(KWHER(i))=i
       End Do

       ! If some elements are equal there will be unfilled places
       ! in the KFROM array: Fill them up now.

     1 Continue

       Do i=1,N
          If(KFROM(i).Eq.0) Then
             IEQUAL=i-1
             IFILLL=i-1
             GO TO 2
          End IF
       End Do

       GO TO 3

     2 Do i=1,N

          If(KWHER(i).Eq.IEQUAL) Then
             KFROM(IFILLL)=i
             IFILLL=IFILLL+1
          End If

       End Do

       GO TO 1

     3 Continue

       ! Finally we may rewrite the input arrays in the right order
       Do i=1,N
          EL(i)=EM(KFROM(i))
          KL(i)=KM(KFROM(i))
       End Do

     End Subroutine sort_vector

     !---------------------------------------------------------------------!
     !  This subroutine sorts a vector EL(1:N) by ascending order and      !
     !  returns the vector of indices in KL(1:N)                           !
     !---------------------------------------------------------------------!
     Subroutine isort_vector(EL,KL,N)

       Integer(iprec), INTENT(IN) :: N
       Integer(iprec), Dimension(1:N), INTENT(INOUT) :: KL
       Integer(iprec), Dimension(1:N), INTENT(INOUT) :: EL

       Integer(iprec) :: i,j,IEQUAL,IFILLL
       Integer(iprec), Allocatable :: KM(:),KWHER(:),KFROM(:)
       Integer(iprec), Allocatable :: EM(:)

       Allocate(KM(1:N),KWHER(1:N),KFROM(1:N))
       Allocate(EM(1:N))

       ! Rewriting input arrays to auxiliary ones
       KWHER(1:N)=1
       EM(1:N)=EL(1:N)
       KM(1:N)=KL(1:N)

       ! How many elements are smaller than the i-th one? In other words
       ! where should it go with respect to the offset position 1?
       Do J=1,N
          Do I=1,N
             If(EM(j).Lt.EM(i)) KWHER(i)=KWHER(i)+1
          End Do
       End Do

       ! Now we only have to resolve the case of identical elements
       ! From where the i-th element should come?
       KFROM(1:N)=0
       Do i=1,N
          KFROM(KWHER(i))=i
       End Do

       ! If some elements are equal there will be unfilled places
       ! in the KFROM array: Fill them up now.

     1 Continue

       Do i=1,N
          If(KFROM(i).Eq.0) Then
             IEQUAL=i-1
             IFILLL=i-1
             GO TO 2
          End IF
       End Do

       GO TO 3

     2 Do i=1,N

          If(KWHER(i).Eq.IEQUAL) Then
             KFROM(IFILLL)=i
             IFILLL=IFILLL+1
          End If

       End Do

       GO TO 1

     3 Continue

       ! Finally we may rewrite the input arrays in the right order
       Do i=1,N
          EL(i)=EM(KFROM(i))
          KL(i)=KM(KFROM(i))
       End Do

     End Subroutine isort_vector

     !------------------------------------------!
     ! ROUTINES FOR MANIPULATING LISTS OF PAIRS !
     !------------------------------------------!

     ! Creates a new list of pairs 'new_pair_vector' of size 'N' with all pair
     ! indices initialized to 0
     Subroutine PairSet_new(N, new_pair_vector)

       Integer(iprec), INTENT(IN) :: N
       Type(PairSet), INTENT(INOUT) :: new_pair_vector
       Type(Pair) :: pair_null

       Nullify(new_pair_vector%pair_vector)

       pair_null%mu = 0
       pair_null%nu = 0

       Allocate(new_pair_vector%pair_vector(N))
       new_pair_vector%pair_vector(:) = pair_null

     End Subroutine PairSet_new

     ! Print the members of the list of pairs 'new_pair_vector' of size 'N'
     Subroutine PairSet_print(N, new_pair_vector)

       Integer(iprec), INTENT(IN) :: N
       Type(PairSet), INTENT(In) :: new_pair_vector
       Integer(iprec) :: i

       Do i=1,N
          Write(6,'(i4," - (",i4,",",i4,")")') i,new_pair_vector%pair_vector(i)%mu,&
                                                 new_pair_vector%pair_vector(i)%nu
       End Do

     End Subroutine PairSet_print

     ! Computes the actual size of the list of pairs 'new_pair_vector' of
     ! original size 'N'
     Integer(iprec) Function PairSet_size(N, new_pair_vector)

       Integer(iprec), INTENT(IN) :: N
       Type(PairSet), INTENT(INOUT) :: new_pair_vector
       Integer(iprec) :: i,count

       count=0
       Do i=1,N
          If(new_pair_vector%pair_vector(i)%mu.Ne.0.And.new_pair_vector%pair_vector(i)%nu.Ne.0) count=count+1
       End Do
       PairSet_size=count

     End Function PairSet_size

     ! Sorts the list of pairs 'new_pair_vector' of size 'N' by ascending order
     ! over the first then second index of the pair
     Subroutine PairSet_sort(N, new_pair_vector)

       Integer(iprec), INTENT(IN) :: N
       Type(PairSet), INTENT(INOUT) :: new_pair_vector
       Integer(iprec) :: i,j,count
       Integer(iprec), Allocatable :: vec_mu(:),vec_nu(:),buf_nu(:),save_nu(:)
       Integer(iprec), Allocatable :: table1(:), table2(:)

       Allocate(vec_mu(1:N),vec_nu(1:N))
       Allocate(table1(1:N),save_nu(1:N))
       Do i=1,N
          vec_mu(i)=new_pair_vector%pair_vector(i)%mu
          vec_nu(i)=new_pair_vector%pair_vector(i)%nu
          table1(i)=i
       End Do
       ! Sort first index by ascending order
       Call isort_vector(vec_mu,table1,N)
       Do i=1,N
          save_nu(i)=vec_nu(table1(i))
       End Do

       Allocate(table2(1:N),buf_nu(1:N))
       i=0; count=0
       Do While (i < N-1 )
          i=i+1
          If(vec_mu(i).Eq.vec_mu(i+1)) Then
             count=count+1
             buf_nu(count)=save_nu(i)
             table2(count)=count
          Else
             count=count+1
             buf_nu(count)=save_nu(i)
             table2(count)=count
             Call isort_vector(buf_nu,table2,count)
             Do j=1,count
                vec_nu(i-count+j)=buf_nu(j)
             End Do
             count=0
          End If
       End Do

       Do i=1,N
          new_pair_vector%pair_vector(i)%mu=vec_mu(i)
          new_pair_vector%pair_vector(i)%nu=vec_nu(i)
       End Do

       Deallocate(vec_mu,vec_nu,buf_nu,save_nu,table1,table2)

     End Subroutine PairSet_sort

     ! Remove all pairs with members equal to either 'mu' or 'nu' in a list
     ! 'new_pair_vector' of size 'N'. The size of the resulting list is decreased
     ! automatically, so that N(output) <= N(input)
     Subroutine PairSet_remove(N, new_pair_vector, mu, nu)

       Integer(iprec), INTENT(IN) :: mu, nu
       Integer(iprec), INTENT(INOUT) :: N
       Type(PairSet), INTENT(INOUT) :: new_pair_vector
       Type(PairSet) :: buffer
       Integer(iprec) :: i,count

       ! If current pointer to the list is not associated, the list is already
       ! empty, we return
       If(.not.associated(new_pair_vector%pair_vector).Or.N.Eq.0) Then
          N=0
          Return
       End if

       count=0
       Do i=1,N
          If(new_pair_vector%pair_vector(i)%mu.Ne.mu .And. new_pair_vector%pair_vector(i)%nu.Ne.nu .And. &
             new_pair_vector%pair_vector(i)%mu.Ne.nu .And. new_pair_vector%pair_vector(i)%nu.Ne.mu) &
          count=count+1
       End do

       ! If the list is not empty after removing the necessary members, we
       ! proceed
       If(count>0) Then
       	 ! Create a temporary new list with the proper size
          Call PairSet_new(count, buffer)
          ! Copy relevant members of original list into temporary new list
          count=0
          Do i=1,N
             If(new_pair_vector%pair_vector(i)%mu.Ne.mu .And. new_pair_vector%pair_vector(i)%nu.Ne.nu .And. &
                new_pair_vector%pair_vector(i)%mu.Ne.nu .And. new_pair_vector%pair_vector(i)%nu.Ne.mu) Then
                count=count+1
                buffer%pair_vector(count)%mu=new_pair_vector%pair_vector(i)%mu
                buffer%pair_vector(count)%nu=new_pair_vector%pair_vector(i)%nu
             End If
          End do
          ! Delete original list
          Call PairSet_del(new_pair_vector)
          ! Update the size of the new list
          N=count
          ! If the new list is not empty, we copy the temporary list into the
          ! new list
          If(N>0) Then
             Call PairSet_new(N, new_pair_vector)
             Do i=1,N
                new_pair_vector%pair_vector(i)%mu = buffer%pair_vector(i)%mu
                new_pair_vector%pair_vector(i)%nu = buffer%pair_vector(i)%nu
             End Do
             Call PairSet_del(buffer)
          Else
             ! If the new list is empty, we update the size
             N=0
          End If
       Else
          ! If the new list is empty, we update the size
          N=0
       End If

     End Subroutine PairSet_remove

     ! Test whether a list of pairs 'new_pair_vector_1' of length 'N1' is the
     ! same as the list of pairs 'new_pair_vector_2' of length 'N2' (lists are
     ! unordered)
     Logical Function PairSet_equal(N1, new_pair_vector_1, N2, new_pair_vector_2)

       Integer(iprec), INTENT(IN) :: N1,N2
       Type(PairSet), INTENT(IN) :: new_pair_vector_1,new_pair_vector_2
       Integer(iprec) :: i,j,sum,m,n
       Integer(iprec), Allocatable :: mu1(:),nu1(:),mu2(:),nu2(:)

       PairSet_equal = .False.

       If(N1.Ne.N2) Return

       Allocate(mu1(1:N1),nu1(1:N1),mu2(1:N1),nu2(1:N1))
       Do i=1,N1
          mu1(i)=new_pair_vector_1%pair_vector(i)%mu
          nu1(i)=new_pair_vector_1%pair_vector(i)%nu
          mu2(i)=new_pair_vector_2%pair_vector(i)%mu
          nu2(i)=new_pair_vector_2%pair_vector(i)%nu
       End Do

       sum=0
       Do i=1,N1
          m=mu1(i); n=nu1(i)
          Do j=1,N1
             If(mu2(j).Eq.m.And.nu2(j).Eq.n) Then
                sum=sum+1
                Exit
             End If
          End Do
       End Do
       Deallocate(mu1,nu1,mu2,nu2)
       If(sum.eq.N1) PairSet_equal = .True.

     End Function PairSet_equal

     ! Removes list of pairs 'new_pair_vector' by deallocating the space and
     ! nullifying the pointer
     Subroutine PairSet_del(new_pair_vector)

       Type(PairSet), INTENT(INOUT) :: new_pair_vector

       Deallocate(new_pair_vector%pair_vector)
       Nullify(new_pair_vector%pair_vector)

     End Subroutine PairSet_del

     !---------------------------------------------------!
     ! ROUTINES FOR MANIPULATING LISTS OF LISTS OF PAIRS !
     !---------------------------------------------------!

     ! Create a new list 'new_liste' of 'M' lists, each of size 'N'
     Subroutine ListePairSet_new(M, N, new_liste)

       Integer(iprec), INTENT(IN) :: M,N
       Type(ListePairSet), INTENT(INOUT) :: new_liste
       Integer(iprec) :: i

       Nullify(new_liste%liste)
       Allocate(new_liste%liste(M),new_liste%size_tiling(M))
       Do i=1,M
          new_liste%size_tiling(i)=N
          Call PairSet_new(N, new_liste%liste(i))
       End Do

     End Subroutine ListePairSet_new

     ! Add a pair 'new_pair' to a list of lists 'new_liste' at the position
     ! 'pair_index' of list number 'tiling_number'
     Subroutine ListePairSet_add_to_tilings(new_pair, new_liste, tiling_number, pair_index)

       Integer(iprec), INTENT(IN) :: tiling_number,pair_index
       Type(ListePairSet), INTENT(INOUT) :: new_liste
       Type(Pair), INTENT(IN) :: new_pair

       new_liste%liste(tiling_number)%pair_vector(pair_index) = new_pair

     End Subroutine ListePairSet_add_to_tilings

     ! Test whether a pair 'new_pair' is anywhere in the list number
     ! 'tiling_number' (of size 'N') of the list of list 'new_liste'
     Logical Function ListePairSet_find_pair(N, new_pair, new_liste, tiling_number)

       Integer(iprec), INTENT(IN) :: N,tiling_number
       Type(ListePairSet), INTENT(IN) :: new_liste
       Type(Pair), INTENT(IN) :: new_pair
       Integer(iprec) :: i,mu0,nu0
       Integer(iprec), Allocatable :: mu(:),nu(:)

       ListePairSet_find_pair = .False.

       Allocate(mu(1:N),nu(1:N))
       Do i=1,N
          mu(i)=new_liste%liste(tiling_number)%pair_vector(i)%mu
          nu(i)=new_liste%liste(tiling_number)%pair_vector(i)%nu
       End Do
       mu0=new_pair%mu; nu0=new_pair%nu

       Do i=1,N
          If((mu0.Eq.mu(i).And.nu0.Eq.nu(i)).Or.(nu0.Eq.mu(i).And.mu0.Eq.nu(i))) Then
              ListePairSet_find_pair = .True.
              Exit
          End If
       End do
       Deallocate(mu,nu)

     End Function ListePairSet_find_pair

     ! Print the list (of size 'M') of lists of pairs (of size 'N' each)
     ! 'new_liste'
     Subroutine ListePairSet_print(M, new_liste)

       Integer(iprec), INTENT(IN) :: M
       Type(ListePairSet), INTENT(INOUT) :: new_liste
       Integer(iprec) :: i,Nact

       Do i=1,M
          Nact = new_liste%size_tiling(i)
          Write(6,'("Tiling number ",i10," of size ",i10)') i,Nact
          Call PairSet_print(Nact, new_liste%liste(i))
       End Do

     End Subroutine ListePairSet_print

     ! Remove a list of 'M' list of pairs 'new_liste'
     Subroutine ListePairSet_del(M, new_liste)

       Integer(iprec), INTENT(IN) :: M
       Type(ListePairSet), INTENT(INOUT) :: new_liste
       Integer(iprec) :: i

       Do i=1,M
          Call PairSet_del(new_liste%liste(i))
       End Do
       Nullify(new_liste%liste)
       Nullify(new_liste%size_tiling)

     End Subroutine ListePairSet_del

 End Module pair_definition

 !------------------------------------------------------------------------------!
 !                                                                              !
 !   Module defining various derived types related to pairs, lists of pairs     !
 !   and lists of lists of pairs, as well as the routines needed to manipulate  !
 !   these objects.                                                             !
 !                                                                              !
 !------------------------------------------------------------------------------!

 Module tilings_utilities

     Use pair_definition

     Implicit None

     Logical :: debug_pair = .False.

 Contains

     !--------------------------------------------------------------------------!
     !  This subroutine takes as input a list 'tilings' of 'M' lists of pairs   !
     !  of N' elements and remove all duplicates. It returns the new list in    !
     !  the variable 'tilings_unique'.                                          !
     !--------------------------------------------------------------------------!
     Subroutine tilings_remove_duplicates(N, M, different, tilings, tilings_unique)

       Integer(iprec), INTENT(INOUT) :: N,M,different
       Type(ListePairSet), INTENT(INOUT) :: tilings,tilings_unique

       Integer(iprec) :: i,j,equals,unique,Nact
       Integer(iprec), Allocatable :: duplicate(:)

       Type(PairSet) :: new_pair_vector_1,new_pair_vector_2

       ! Counting the number of duplicate tilings
       Allocate(duplicate(1:M)); duplicate(:)=0

       Call PairSet_new(N, new_pair_vector_1)
       Call PairSet_new(N, new_pair_vector_2)

       equals=0
       Do i=1,M
          new_pair_vector_1 = tilings%liste(i)
          Do j=i+1,M
             new_pair_vector_2 = tilings%liste(j)
             If(PairSet_equal(N, new_pair_vector_1, N, new_pair_vector_2).And.duplicate(j).Eq.0) Then
                duplicate(j)=1
                equals=equals+1
                Exit
             End If
          End Do
       End Do
       If(debug_pair) Then
          Write(6,'("In tilings_remove_duplicates() - Total number of tilings: ",i12)') M
          Write(6,'("In tilings_remove_duplicates() - Number of equal tilings: ",i12)') equals
       End If

       ! Retaining only the unique mappings
       unique = M - equals

       Call ListePairSet_new(unique, N, tilings_unique)

       different=0
       Do i=1,M
          If(duplicate(i).eq.0) Then
             different=different+1
             Nact = PairSet_size(N, tilings%liste(i))
             tilings_unique%liste(different)%pair_vector = tilings%liste(i)%pair_vector
             tilings_unique%size_tiling(different) = Nact
          End If
       End Do

       ! Cleanly deallocating space
       Nullify(new_pair_vector_1%pair_vector)

       Deallocate(duplicate)

     End Subroutine tilings_remove_duplicates

     !--------------------------------------------------------------------------!
     !  This subroutine creates a list of list of pairs 'mapping' based on the  !
     !  list of pairs 'new_pair_vector' of size 'N'                             !
     !--------------------------------------------------------------------------!
     Subroutine tilings_find_all(N, M, tilings, new_pair_vector)

        Integer(iprec), INTENT(INOUT) :: N, M
        Type(PairSet), INTENT(INOUT) :: new_pair_vector
        Type(ListePairSet), INTENT(INOUT) :: tilings

        Logical :: duplicates
        Integer(iprec) :: i,j,k,q,tiling_number,mu,nu,Ncur,Norig,pair_index,Nact,a,b,Nmax
        Integer(iprec), Allocatable :: vec1_mu(:),vec1_nu(:),vec_mu(:),vec_nu(:),itaken(:)
        Type(Pair) :: p

        tiling_number=0; Norig=N

        Allocate(vec_mu(Norig),vec_nu(Norig))
        Nmax=0 ! Maximum integer in the list
        Do q=1,Norig
           vec_mu(q) = new_pair_vector%pair_vector(q)%mu
           vec_nu(q) = new_pair_vector%pair_vector(q)%nu
           If(vec_mu(q).Ge.Nmax) Nmax=vec_mu(q)
           If(vec_nu(q).Ge.Nmax) Nmax=vec_nu(q)
        End Do

        Do i=1,Norig

           ! Reset original list of pairs
           Call PairSet_del(new_pair_vector)
           Call PairSet_new(Norig, new_pair_vector)
           Do q=1,Norig
              new_pair_vector%pair_vector(q)%mu = vec_mu(q)
              new_pair_vector%pair_vector(q)%nu = vec_nu(q)
           End Do

           ! Consider pair number i in the list
           p=new_pair_vector%pair_vector(i)
           mu=p%mu; nu=p%nu

           ! This will define a new tiling: we increment the tiling index and
           ! reset the pair index within this tiling
           tiling_number=tiling_number+1
           pair_index=1

           ! Add current pair to the new tiling
           Call ListePairSet_add_to_tilings(p, tilings, tiling_number, pair_index)

           ! Remove all indices involved in current pair from current list
           Ncur=Norig
           Call PairSet_remove(Ncur, new_pair_vector, mu, nu)

           ! Save resulting list of pairs
           Allocate(vec1_mu(Ncur),vec1_nu(Ncur))
           Do q=1,Ncur
              vec1_mu(q)=new_pair_vector%pair_vector(q)%mu
              vec1_nu(q)=new_pair_vector%pair_vector(q)%nu
           End Do

           ! Check if there are duplicates in the new list, i.e. indices that
           ! are repeated in more than one pair. Ex.: In [(2,3),(2,4),(5,7)],
           ! index 2 is repeated in 2 different pairs: there are duplicates.
           Allocate(itaken(1:Nmax)); itaken(:)=0
           duplicates = .False.
           Do j=1,Ncur
              a=new_pair_vector%pair_vector(j)%mu
              b=new_pair_vector%pair_vector(j)%nu
              itaken(a)=itaken(a)+1
              itaken(b)=itaken(b)+1
              if(itaken(a).Gt.1.or.itaken(b).gt.1) duplicates = .True.
           End Do
           Deallocate(itaken)

           ! The presence of duplicates means that there are more than one
           ! possible tiling associated with the current list: we will scan the
           ! list to identify all such tilings
           If(duplicates) Then
              Nact=Ncur
           Else
              Nact=1
           End If

           j=1
           Do While (j<=Nact)

              ! Move to position j of current list
              k=j; N=Ncur

              ! We identify all unique pairs in the current list, eliminating
              ! them as we proceed
              Do While (k<=N)

                 ! Take pair at position j (k=j at the beginning of this loop)
                 a=new_pair_vector%pair_vector(k)%mu
                 b=new_pair_vector%pair_vector(k)%nu

                 ! If the current pair is not already in the current tiling, we
                 ! add it to it and remove the corresponding indices from the
                 ! list
                 If(.Not.ListePairSet_find_pair(N,new_pair_vector%pair_vector(k),tilings,tiling_number)) Then
                    pair_index=pair_index+1
                    Call ListePairSet_add_to_tilings(new_pair_vector%pair_vector(k), tilings, tiling_number, pair_index)
                    Call PairSet_remove(N, new_pair_vector, a, b)
                 Else
                 ! If the current pair has already been included, we move to the
                 ! next element of the list
                    k=k+1
                 End If

              End Do
              tilings%size_tiling(tiling_number) = pair_index

              ! We restore the list (which has been emptied just before) with
              ! its original content
              Call PairSet_del(new_pair_vector)
              Call PairSet_new(Ncur, new_pair_vector)
              Do q=1,Ncur
                 new_pair_vector%pair_vector(q)%mu = vec1_mu(q)
                 new_pair_vector%pair_vector(q)%nu = vec1_nu(q)
              End Do

              ! Go one step further down the list
              j=j+1

              ! If there are duplicates, we will try to find the next tiling: we
              ! increment the tiling index and reset the pair index
              If(duplicates.and.j<=Nact) Then
                 tiling_number=tiling_number+1
                 pair_index=1
                 Call ListePairSet_add_to_tilings(p, tilings, tiling_number, pair_index)
              End If

           End Do ! end j

           Deallocate(vec1_mu,vec1_nu)

        End Do ! end i

        Deallocate(vec_mu,vec_nu)

        N = Norig
        M = tiling_number

     End Subroutine tilings_find_all

 End Module tilings_utilities

