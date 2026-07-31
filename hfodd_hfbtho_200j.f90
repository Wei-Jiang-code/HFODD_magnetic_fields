!***********************************************************************
!
!    Copyright (c) 2016, Lawrence Livermore National Security, LLC.
!                        Produced at the Lawrence Livermore National
!                        Laboratory.
!                        Written by Nicolas Schunck, schunck1@llnl.gov
!
!    LLNL-CODE-710577 All rights reserved.
!    LLNL-CODE-573953 All rights reserved.
!
!    Copyright 2016, N. Schunck, J. Dobaczewski, W. Satula, P. Baczyk,
!                    J. Dudek, Y. Gao, M. Konieczka, K. Sato, Y. Shi,
!                    X.B. Wang, T.R. Werner
!    Copyright 2012, M.V. Stoitsov, N. Schunck, M. Kortelainen, H.A. Nam,
!                    N. Michel, J. Sarich, S. Wild
!    Copyright 2005, M.V. Stoitsov, J. Dobaczewski, W. Nazarewicz, P.Ring
!
!    This file is part of HFBTHO.
!
!    HFBTHO is free software: you can redistribute it and/or modify it
!    under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    HFBTHO is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with HFBTHO. If not, see <http://www.gnu.org/licenses/>.
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
!==================================================================================================================================
!#START HFBTHO_VERSION MODULE
!==================================================================================================================================
Module HFBTHO_VERSION
  Implicit None
  Character(6) :: Version='200i'
  !--------------------------------------------------------------------------------------
  ! Version History
  !--------------------------------------------------------------------------------------
  ! ver#200e: added preprocessor option SWITCH_ESSL to use IBM libraries
  ! ver#200d: fixed bug in gfv; improved legibility and accuracy of coulom and coulom1
  ! ver#200c: added LLNL release number
  ! ver#200b: added module linear_algebra, analyzing THO, formatted output, fixed bug
  !           in calculation of entropy
  ! ver#200a: Restored option to compute all blocking configurations within given energy
  !           window, removed spurious output printing
  ! ver#199: added a few input options and a compatibility mode with HFODD. Release
  !          candidate before publication, last bug to fix is OpenMP in hfbdiag()
  ! ver#142: removed module pairing, module UNEDF, spurious preprocessing options
  !          for publication purposes; used the Lahey compiler to identify a few bugs,
  !          added routine check_consistency(), module HFBTHO_gauss, improved file
  !          handling in inout() and fixed bug in multipole moments
  ! ver#141: Reinstated the THO module
  ! ver#140: added a namelist for debugging purposes, added OpenMP in coulom()
  ! ver#139a: fixed bug in readjustment of constraint at finite T, cleaned up output
  !           system, fixed bug in calculation of Coulomb in parity breaking mode
  ! ver#139: added temperature
  ! ver#138: added new module HFBTHO_utilities to improve portability
  ! ver#137: new system of inputs based on several namelists contained in one unique
  !          file called hfbtho_NAMELIST.dat, multiple constraints up and running
  ! ver#136: ANL OpenMP optimizations included, and code clean-up
  ! ver#135: Tested with all previous versions of HFBTHO down to 101, full compatibility
  !          achieved
  ! ver#134: Begining of work toward publication
  ! ver#133: Single-file HFBTHO
  ! ver#130: tho.dat mdifications due to blocking, error indicator introduced
  ! ver#129: Even-even tested and equivalent with ptho101spt15sp.f90 used in ANL fit
  ! ver#128: EQP,U,V and their dimentions NUV,NEQ required for qrpa incoporated
  !          permanently in HFBTHO substituting old arrays eqp and uv
  ! ver#127: For easy development the module split in different F90 files which are
  !          invoked using INCLUDE statements (remove that when developement is over)
  !          Optimized qrpa_DENSIT_PLUS and qrpa_GAMDEL to be twice faster
  ! ver#126: Cleaning, optimizing, and isolating THO stuff
  !                parity good:   Time per iteration: 3.841 seconds
  !                parity broken: Time per iteration: 9.933 seconds
  !          HFB+THO tested and works in both parity regimes. iserial removed and
  !          substituted with Print_Screen, i.e, record results only when Nsh>0
  ! ver#125: Implemented and tested reflection symmetry as option. If parity is broken,
  !          computer time per iteration is almost 5 times bigger:
  !                parity good:   Time per iteration:  4.5276 seconds
  !                parity broken: Time per iteration: 19.4646 seconds
  !                difference in total energy:  0.001  keV
  ! ver#124: Rewrite to prepare for breaking reflection symmetry
  !          Preprocessing directives included: #ifndef hide_qrpa, hide_tho, hide_dme
  !          For preprocessing one needs -Dhide_qrpa -Dhide_tho -Dhide_dme
  !          If no preprocessing or -Uhide_qrpa -Uhide_tho -Uhide_dme then
  !          all modules are included
  ! ver#123: Playground for QRPA calculations [16/12/2010]:
  !          The changes are:
  !          - all variables in Module HFBTHO are Public
  !          - include_qrpa=0 is added to Module HFBTHO with asssociated
  !            declarations eventualy used by qrpa
  !          - Subroutine ByNucleus moved to PTHO_PROGRAM where is its place
  !            and it should be done long ago. Call Do_QRPA() is used only there.
  !          - So if the program is compiled with -Uhide_qrpa one can use Do_QRPA()
  !            to do qrpa calculations.
  ! ver#123: Fixed crash after iterations limit.
  !          Tested against anl version hfbtho101spt15.f90 - itterations go differently
  !          but the final results are identical.
  ! ver#122: (MK) Added CExPar for coulomb exchange. Parameter read from UNEDF module
  ! ver#121: (MK) Added possibility to use zero particle number for droplet calculations
  ! ver#120: (MK) Added external field, and all channels to direct Hartree. e^2 for Coulomb
  !          now read from the UNEDF module. Direct Hartree now always calculated based on
  !          module function regardless of the value of DMEorder parameter
  ! ver#117: Direct Hartree added when DME_order>-1
  ! ver#115: (MK) added use_cm_cor variable to hb0 calculation and (nabla rho)^2 terms to
  !          the calculate_U_parameters function calls
  ! ver#114: Name list, new tho.dat file, proton/neutron fields, confirms all results of
  !          recent published version after ANL optimization ptho101b_last_tested.f90,
  !          public/Public variables
  ! ver#113: Cleaning
  ! ver#112: No parameter functions
  ! ver#111: Main program detached from the file as PTHO_MAIN_PROGRAM.f90 which will not be
  !          versioned. ptho becomes jus a HFBTHO module.  Pairing constants V0(2),V1(2)
  !          replaced by CpV0(0:1), CpV1(0:1) coming as public from defined in UNEDF module
  !          Removed dalf and ippforce form the pairing. For compatibility, ippforce
  !          stays in the input file by now but the kind of pairing is given by CpV1 only
  !          Dropped corrections 'ecmcpavpj', 'erotcorrection' which should be added later.
  !          For compatibility inputfile stays the same. Added IDEUB.
  !          THO part in 'densit' (not densitpj), 'gamdel' commented HO/THO for speed
  ! ver#109: All public variables, expectpj works with a jump: not clear how UNEDF can work
  !          with complex numbers, just skip this part by now bu write results data
  !          New thodefh(iw1)
  ! ver#108: Removed all programs not used in ver#107
  !          expect contains a key DO_FITT:
  !            =0 calculare energy, delta, def  & rms only
  !            =1 the same+all integrals for the regression optimization
  !          V0,V1 pairing constant separated for neutrons and protons: v0(2),v1(2)
  !          HFBTHO collected in MODULE HFBTHO
  !          KOP3 removed
  ! ver#107: Towards UNEDF: complete rewrite based on Marcus to include N2LO
  !          LN for ZR110 at prolate solution with SLY4, mixed pairing and tensor terms:
  !          -agreement with previouse Skyrme implemetation to the last significant digit
  !          -agreement with previouse LO+LDA implemetation to the last significant digit
  !          -agreement with previouse LO+CB  implemetation to the last significant digit
  ! ver#106: Towards UNEDF: the standard functional rewritten in terms of UNEDF U-amplitudes
  !          The assumption U=U(tau_0,Delta rho_0,rho_0,rho_1) becomes possible after
  !          adding Nabla rho_ij terms  (STANDARD FUNCTIONAL ONLY)
  ! ver#105: Towards UNEDF: the standard functional rewritten in terns of UNEDF U-amplitudes
  !          The assumption U=U(rho_0,rho_1) becomes possible after adding Delta rho_ij terms
  ! ver#104: Broyden improved with linear search at negative curvature
  !          Implemented Agumented Lagrangian method for constraint calculations
  !          Manual blocking included and tested, key: manualBlocking
  ! ver#103: From this version on-no more support for VAP (VAP completely removed)
  !          The whole program in terms of C-parameters (including tenzor terms)
  ! ver#102: The whole program in terms of C-parameters (without tenzor terms)
  ! ver#101: Optimization in terms of nuclear matter: 'FITS' regime
  ! ver#100: Toward isovector pairing following Sagawa and Yamagami
  ! ver# 99: Subroutine HFBiterations. The isotopic line in tho.dat removed.
  !          Subroutines byNucleus, byConstraint, FitPairing, HFBTHO_HFODD isolated
  !          at the end and could be ported if necessary. skyrme='FITS' assumes the skyrme
  !          parameters as explicitely given. -N00 supresses completely the output and only
  !          hodef.dat and thodef.dat are charged (if iserial=0 even these files are supressed)
  !          HFBTHO_HFODD updated (think further about a constraint in Q2 terms)
  ! ver# 98: INOUT modified and added interface subroutine HFBTHO_HFODD
  !          Pairing fitted with MIX/(LN-NOLN) for SLY4,SKP,SKM* forces
  ! ver# 97: Corrected blocking candidates criteria
  ! ver# 96: Extended so term W0,W1 and SKLY4T forces
  ! ver# 95: Pairing regularization removed, linear HFBDIAG mixing for Lambda
  !          when blocking, DSYEV replaced by the faster DSYEVD, hfbdiag caculates
  !          canonical basis only at the last hfbdiag iteration, expect optimized
  !          new subroutines HFBiterations, FitPairing, byConstraint
  !          Work around a bug in LAPAK related to DSYEVD
  ! ver# 94: Misprints
  ! ver# 93: Removed hh and de matrices and related manipulations
  !          Broyden_min now escape maximum and inflex points
  ! ver# 92: Bug in blocking while no pairing cleaned
  ! ver# 91: BLAST  & LAPACK diagonalization
  ! ver# 90: If applying Broyden method to matrix elements then
  !          at 20 shells the total number of matrix elements
  !          is 2x2x65307=261228 or about 2.1 Mb and if one keeps
  !          8 iterations it will be  about 17 Mb-not too much
  !          This is the only way one can mix Lipkin-Nogami
  !
  !          If one uses the potentials at 30x30 grid points
  !          the numbers are 8x2x900=14400 or 115 Kb and if one keeps
  !          8 iterations it will be about 1 Mb-much better
  !          but Lipkin-Nogami is out of this scheme (?!)
  !          If one uses densities at 30x30 grid points they
  !          are 9x2x900-almost like the potential case.
  !          (sent to George)
  ! ver# 89: Reduced printout (no anymore lprinter)
  !          LN del+ala2 printed during the iterations
  !          Strength in the initial constraint calculations reduced to 0.3, requested
  !          deformation+/-0.3, untill si<1.1
  !          If too slow convergence (1000 iterations) and Lambda>0 interrupt iterations
  !          Odd nucleus right away from the even-even (even) solution
  !          When even solution missing/corrupt (even at inin<0) calculate it first
  !          and then odd one
  ! ver# 88: Synchronization for the parallel run
  !          FileLabel subroutine added.
  !          Modified inin control
  !            inin<0: Always start from a file if it exists, not corrupted and correct
  !                    otherwise inin=iabs(inin) and start from scratch
  !            inin>0: Always start from scratch
  !          SCRATCH calculations start with initial 20 constrained iteration if constrain
  !          is not requested (icstr=0). When icstr#0 standard constraint calculations.
  !          BY CHAIN calculations temporary removed due to blocking complications
  ! ver# 87: Approximate Blocking keeping time-reversal symmetry to PNP PAV
  ! ver# 86: Approximate Blocking keeping time-reversal symmetry and tested agings HFODD
  ! ver# 85: LN in canonical basis. Benchmark to HFODD
  ! ver# 84: Testing HFODD LN again HFB-HO
  ! ver# 83: Cleaning, SKM* mixed volume LN pairing fitted
  ! ver# 82: As 81 but prepared for jaguar
  ! ver# 81: Pairing regularization/renormalization. PAV done with unprojected v_k
  !          V0(Nsh=20,pwi=50) fitted for SLY4,SKP,Renormalized,Regularized,Mixed,Volume
  !          Removed delta and gamdel0 completely
  ! ver# 80: Accuracy for large number of shells increased by the number of gauss points
  !          Gaussian points now calculated
  !          Initial guest now from deformed Wood-Saxon
  !          Initial run now starts with requested shell number n00
  ! ver# 79: Cranking rotational correction implemented:
  !          Printed to screen, thoout.dat, hodef.dat and thodef.dat
  !          but not added to the energy
  ! ver# 78: Full CM correction implemented in HFB  & HFB(PAV)
  !          Printed to screen, thoout.dat, hodef.dat and thodef.dat
  !          but not added to the energy
  !          NB : ilpjnp(2) removed
  ! ver# 77: Automated Blocking:
  !          First is calculated N,Z without blocking, remembered *.hel *.tel files and
  !           determined blocking candidates according to pwablo criteria.
  !           if N(Z) is odd we have neutron(proton) blocking candidates.
  !           if both, N and Z, are odd we have both, neutron and proton, block.candidates
  !          Then we block state after state among the blocking candidates and calculate
  !           starting from the recorded unblocked (N,Z) solution
  !           if only N ( or Z) is odd then all neutron (or proton) blocking candidates
  !           are calculated
  !           if both N and Z are odd all pairs of proton and neutron blocking candidates
  !           are calculated                                               (PAV  &LN unclear)
  ! ver# 76: Manual Blocking for a particular state in a particular block
  !          overlap criteria used to avoid the level crossing problem    (PAV  &LN unclear)
  ! ver# 75: Manual Blocking for the minimal qpe within a given block.
  ! ver# 74: Bulgac procedure .. not done
  ! ver# 73: Cleanup, introducing the cpc notations, beyond unit circle removed (MARK 1)
  ! ver# 72: PNP: still valid version for integration over the unit circle
  ! ver# 71: PNP: towards beyond unit circle integration
  ! ver# 70: PNP: detached neutron from proton projection
  ! ver# 69: PNP: quadrupole constraint fixed to converge
  ! ver# 68: equivalent to var#67
  ! ver# 67: PNP: corrections to the tensor term and initial dumping factor
  ! ver# 66: PNP: V0 fitted to PLN energies at Sn126 for SKP mixed and volume at Nsh=20,HO
  ! ver# 65: deformed HO basis implemented and tested
  ! ver# 64: fixed byChain to go not 2 beyond the forced break
  ! ver# 63: iasswrong(3) fixed for correct multiprocessor run
  ! ver# 61/62: pthotop for Cheetah added at the end
  ! ver# 60: Thodef.dat header line fixed (added U:). Fixed P/N in ByChain calculations
  !          The 'Stop' is removed from  mishmatch conditions with alternative to use old one.
  !          Consistent pairing for SLY4 and SKP forces
  ! ver# 59: LST modified to accept negative aa-values. SLY4 And SKP with consistent
  !          (high densities regime) pairing for all cases. Old asymptotic prescription
  !          is used in the case of Mishmatch asymptotic. Temporary,still new SLY4
  !          pairing constants are commented.   ! ver# 58: Proton line in byChain
  !          calculations goes vertically. In case of wrong
  !          asymptotic parameter 'kindhfb' is recorded as 'kindhfb+100' in thodef.dat
  !          file where the results for this nucleus are substituted with HO results
  !          Only Nsh=20 pairing constants are already fitted to the higher density
  !          asymptotic prescription which is already enforced. (Temporary, still old
  !          SLY4 pairing constants are in the code).
  ! ver# 57: Partially refitted pairing constants according to the new asymptotic
  !          prescription. NB! old constants are still for the SLY4 force and not all of
  !          the cases with SKP are fitted. Old ass. regime is temporary enforced.
  ! ver# 56: Code optimization and checks
  ! ver# 55: Tensor term J.J implemented and tested
  ! ver# 54: Back to *.hel *.tel files; Including new hodef.dat like thodef.dat file
  ! ver# 53: LST is choosing the higher density in the asymptotic region
  !--------------------------------------------------------------------------------------
End Module HFBTHO_VERSION
!==================================================================================================================================
!#END HFBTHO_VERSION MODULE
!==================================================================================================================================
!#START HFBTHO_utilities
!==================================================================================================================================
Module HFBTHO_utilities

  Implicit None

  Integer, Parameter, Public :: ipr=Kind(1)     ! to set the precision of the DFT solver
  Integer, Parameter, Public :: pr =Kind(1.0)   ! to set the precision of the DFT solver

  ! I/O
  Integer, Public :: lout = 6, lfile = 7

  ! Global numbers
  Real(pr), Parameter :: zero=0.0_pr,half= 0.5_pr,one=1.0_pr,two  =2.0_pr,three=3.0_pr, &
                         four=4.0_pr,five= 5.0_pr,six=6.0_pr,seven=7.0_pr,eight=8.0_pr, &
                         nine=9.0_pr,ten =10.0_pr
  ! Whole global numbers pp#
  Real(pr), Parameter :: pp12=12.0_pr,pp16=16.0_pr,pp15=15.0_pr,pp20=20.0_pr, &
                         pp24=24.0_pr,pp27=27.0_pr,pp32=32.0_pr,pp64=64.0_pr, &
                         pp40=40.0_pr
  ! Fractional global numbers p#
  Real(pr), Parameter :: p12=one/two,   p13=one/three,  p14=one/four,  p23=two/three, &
                         p43=four/three,p32=three/two,  p34=three/four,p53=five/three,&
                         p18=one/eight, p38=three/eight,p59=five/nine, p52=five/two,  &
                         p54=five/four, p74=seven/four

 Contains
  !===============================================================================================
  !
  !===============================================================================================
  Subroutine get_CPU_time (subname,is)
    Implicit None
    Integer, Intent(in)   :: is
    Character*(*), Intent(in)  :: subname
    Real(pr), Save :: time1,time2
    !
    If(is.Eq.0) Then
       Call cpu_time(time1)
    Else
       Call cpu_time(time2)
       Write(lout,'(a,a,a,G16.6)') '  Time in seconds -> ',Trim(subname),':',time2-time1
    End If
    !
  End Subroutine get_CPU_time
  !
End Module HFBTHO_utilities
!==================================================================================================================================
!#END HFBTHO_utilities
!==================================================================================================================================
!#START linear_algebra MODULE
!==================================================================================================================================
Module linear_algebra

  Use HFBTHO_utilities
  Implicit None

Contains
  !=======================================================================
  Subroutine lingd(ma,mx,n,m,a,x,d,Ifl)
    !---------------------------------------------------------------------
    ! Solves the system of linear equations A*X = B
    ! At the beginning the matrix B is stored in X
    ! During the calculation it will be overwritten
    ! D is the determinant of A
    !---------------------------------------------------------------------
    Integer(ipr) :: ma,mx,n,m,Ifl
    Integer(ipr), Save :: i,j,k,l,k1,n1
    Real(pr) ::  a(ma,m),x(mx,m),d
    Real(pr), Save :: tollim,one,zero,p,q,tol,cp,cq
    Data tollim/1.d-10/,one/1._pr/,zero/0._pr/
    Ifl = 1; p = zero
    Do i=1,n
       q = zero
       Do j=1,n
          q = q + Abs(a(i,j))
       End Do
       If(q.Gt.p) p = q
    End Do
    tol = tollim*p; d   = one
    Do k=1,n
       p = zero
       Do j=k,n
          q = Abs(a(j,k))
          If(q.Lt.p) Cycle
          p = q; i = j
       End Do
       If (p.Lt.tol) Then
          Write (6,200) ('-',j=1,80),tol,i,k,a(i,k),('-',j=1,80)
200       Format (/1x,80a1/' *****  ERROR IN LINGD , TOLERANZ =',e14.4, &
               ' VALUE OF A(',i3,',',i3,') IS ',e14.4/1x,80a1)
          Ifl = -1
          Return
       End If
       cp = one/a(i,k)
       If(i.Ne.k) Then
          d = -d
          Do l=1,m
             cq = x(i,l); x(i,l) = x(k,l); x(k,l) = cq
          End Do
          Do l=k,n
             cq = a(i,l); a(i,l) = a(k,l); a(k,l) = cq
          End Do
       End If
       d = d*a(k,k)
       If(k.Eq.n) Exit
       k1 = k + 1
       Do i=k1,n
          cq=a(i,k)*cp
          Do l=1,m
             x(i,l)=x(i,l)-cq*x(k,l)
          End Do
          Do l=k1,n
             a(i,l)=a(i,l)-cq*a(k,l)
          End Do
       End Do
    End Do
    Do l=1,m
       x(n,l)=x(n,l)*cp
    End Do
    If(n.Eq.1) Return
    n1=n-1
    Do k=1,n1
       cp = one/a(n-k,n-k)
       Do l=1,m
          cq = x(n-k,l)
          Do i=1,k
             cq = cq-a(n-k,n+1-i)*x(n+1-i,l)
          End Do
          x(n-k,l) = cq*cp
       End Do
    End Do
    Return
  End Subroutine lingd
  !=======================================================================
  !
  !=======================================================================
  Subroutine csplin(n, x, y, b, c, d)
    !---------------------------------------------------------------------
    ! file: csplin.for  (from slac)
    ! The coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
    ! for a cubic interpolating spline
    ! s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
    ! for  x(i) <= x <= x(i+1)
    ! input..
    ! n = the number of data points or knots (n.ge.2)
    ! x = the abscissas of the knots in strictly increasing order
    ! y = the ordinates of the knots
    ! output..
    ! b, c, d  = arrays of spline coefficients as defined above.
    ! using  p  to denote dIfferentiation,
    ! y(i) = s(x(i))
    ! b(i) = sp(x(i))
    ! c(i) = spp(x(i))/2
    ! d(i) = sppp(x(i))/6  (derivative from the right)
    ! the accompanying function subprogram  cseval  can be used
    ! to evaluate the spline, its derivative or even its 2nd derivative.
    !---------------------------------------------------------------------
    Integer(ipr), Save :: nm1,i,ib
    Integer(ipr) :: n
    Real(pr) :: x(n), y(n), b(n), c(n), d(n)
    Real(pr), Save :: t,zero=0.0_pr,two=2.0_pr,tr=3.0_pr
    ! check input for consistency
    If(n.Lt.2) Stop '-n < 2 in csplin call--'
    nm1 = n-1
    Do i = 1, nm1
       If(x(i).Ge.x(i+1)) Stop 'x not strictly ascending in csplin call'
    End Do
    If (n.Ne.2) Then
       ! set up tridiagonal system
       ! b = diagonal, d = offdiagonal, c = right hand side.
       d(1) = x(2) - x(1); c(2) = (y(2) - y(1))/d(1)
       Do i = 2, nm1
          d(i) = x(i+1) - x(i); b(i) = two*(d(i-1) + d(i))
          c(i+1) = (y(i+1) - y(i))/d(i); c(i) = c(i+1) - c(i)
       End Do
       ! end conditions.  third derivatives at  x(1)  and  x(n)
       ! obtained from divided dIfferences
       b(1) = -d(1); b(n) = -d(n-1); c(1) = zero; c(n) = zero
       If (n.Ne.3) Then
          c(1) =  c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
          c(n) =  c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
          c(1) =  c(1)*d(1)**2/(x(4)-x(1))
          c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
          ! forward elimination
       Else
          Do i = 2, n
             t = d(i-1)/b(i-1); b(i) = b(i) - t*d(i-1); c(i) = c(i) - t*c(i-1)
          End Do
       End If
       ! back substitution
       c(n) = c(n)/b(n)
       Do ib = 1, nm1
          i = n-ib
          c(i) = (c(i) - d(i)*c(i+1))/b(i)
       End Do
       ! compute polynomial coefficients
       b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + two*c(n))
       Do i = 1, nm1
          b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + two*c(i))
          d(i) = (c(i+1) - c(i))/d(i); c(i) = tr*c(i)
       End Do
       c(n) = tr*c(n); d(n) = d(n-1)
       Return
    Else
       b(1) = (y(2)-y(1))/(x(2)-x(1)); c(1) = zero; d(1) = zero
       Return
    End If
  End Subroutine csplin
  !=======================================================================
  !
  !=======================================================================
  Subroutine cseval(n,u,x,y,b,c,d,splf0)
    !---------------------------------------------------------------------
    ! This subroutine is a copy of 'cseva' but only for the function
    !---------------------------------------------------------------------
    Integer(ipr) :: n
    Integer(ipr), Save :: i=1,j,k
    Real(pr) :: x(n),y(n),b(n),c(n),d(n),u,splf0
    Real(pr), Save :: dx
    If(i.Ge.n)      i = 1
    If(u.Lt.x(i))   Go To 10
    If(u.Le.x(i+1)) Go To 30
    ! binary search
10  i = 1
    j = n + 1
20  k = (i+j)/2
    If(u.Lt.x(k)) j = k
    If(u.Ge.x(k)) i = k
    If(j.Gt.i+1) Go To 20
    ! evaluate splf0
30 dx = u - x(i)
    splf0 = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
    Return
  End Subroutine cseval
  !=======================================================================
  !
  !=======================================================================
  Subroutine deri(h,n,f1,dunl)
    !---------------------------------------------------------------------
    ! First derivative of 'f1' if the step is 'h'
    !---------------------------------------------------------------------
    Integer(ipr) :: n
    Integer(ipr), Save :: k
    Real(pr) :: h,f1(n),dunl(n)
    Real(pr), Save :: t60,t12
    Real(pr), Save :: t8=8.0_pr,t45=45.0_pr,t9=9.0_pr
    t60 =1.0_pr/(h*60.0_pr); t12 =1.0_pr/(h*12.0_pr)
    !
    dunl(1)  =(t8*f1(2)-f1(3)+f1(1))*t12
    dunl(2)  =(t45*(f1(3)-f1(1))-t9*f1(4)+f1(5)-f1(1))*t60
    dunl(3)  =(t45*(f1(4)-f1(2))-t9*(f1(5)-f1(1))+f1(6))*t60
    dunl(n)  =(-t8*f1(n-1)+f1(n)+f1(n-2))*t12
    dunl(n-1)=(t45*(f1(n)-f1(n-2))+t9*f1(n-3)-f1(n)-f1(n-4))*t60
    dunl(n-2)=(t45*(f1(n-1)-f1(n-3))-t9*(f1(n)-f1(n-4))-f1(n-5))*t60
    Do k=4,n-3
       dunl(k) =(t45*(f1(k+1)-f1(k-1))-t9*(f1(k+2)-f1(k-2))+f1(k+3)-f1(k-3))*t60
    End Do
    Return
  End Subroutine deri
  !
End Module linear_algebra
!==================================================================================================================================
!#END linear_algebra MODULE
!==================================================================================================================================
!==================================================================================================================================
!#START UNEDF MODULE
!==================================================================================================================================
Module UNEDF
  !--------------------------------------------------------------------------------------
  ! M.Kortelainen & M.Stoitsov, 2009-2011
  ! UNEDF interface for Skyrme, DME(LO,NLO,N2LO) and other DFT solvers
  !--------------------------------------------------------------------------------------
  Use HFBTHO_utilities
  Implicit None
  !
  Character(16), Private :: Version='17'
  !
  ! Version History
  !--------------------------------------------------------------------------------------
  ! ver#17:(Mario)   use_TMR_pairing=0/1 standard/TMR pairing added
  !                  to Namelist. Using:
  !                  CpV0(0)=G,    CpV0(1)=a
  !                  CpV1(0)=vfacn,CpV1(1)=vfacp
  ! ver#16:(Mario)   #ifndef hide_dme preprocessing directive included
  ! ver#15:(Markus)  Added parameter CExPar, used in Coul. excange term.
  !                  Also, all the channels included in direct Hartree
  ! ver#14:(Markus)  Added function Vexternal for the external field,
  !                  and use_j2terms to switch off tensor terms.
  !                  Direct Hartree set to zero.
  ! Ver#13:(Mario)   Added ac2,ac3,acoord
  ! ver#12:(Mario)   hartree term temprorary dropped. rDr NNN terms taken
  !                  with a factor of 1/2
  ! ver#11:(Mario)   Gaussian approximation to the Hartree term added,
  ! [3/10/2010]      hatree_rc removed. NB! Function HartreeV is an
  !                  elemental function with possible array arguments
  ! ver#10: (Markus) Added e2charg (e^2 for Coulomb) to the public variables
  ! ver#9: (Mario)   Hartree 'CHrho' calculated in INM with rc='hatree_rc'
  ! [2/2/2010]       is subtracted from Crho(0)at DMEorder >= 0.
  !                  CHrho added to the public list, 'hatree_rc' added
  !                  to interaction parameters and the namelist.
  !                  In the case DMEorder=-1 (standard Skyrme)
  !                  both, 'CHrho' and 'hatree_rc', do not play.
  !                  New function HartreeV(u) defines Hatree energy as
  !                  E(Hartree)=(1/2)*Int[rho_0(r)*V(|r-r'|)*rho_0(r')]
  !                  HartreeV(u) is zero for u=<'hatree_rc'
  ! ver#8: (Markus)  Hartree DME terms dropped out.
  ! ver#7: (Markus)  Added switch to turn off the 3N terms.
  !        (Mario)   Added Abs to density and gradient dependent LDA
  !                  Public :: DMEorder,DMElda,use_DME3N_terms
  ! ver#6: (Mario)   Skyrme transformation added.
  ! ver#5: (Mario)   Print_Namelist=T/F added to the namelist
  ! ver#4: (Markus)  Added natural units to the module. Used only for printing.
  ! ver#3: (Mario)   Uamplitudes(0:3,0:7) in normal order
  !
  ! t for Uamplitudes(t,*)
  ! 0 -> 0,0
  ! 1 -> 1,1
  ! 2 -> 0,1
  ! 3 -> 1,0
  ! n for Uamplitudes(*,n)
  ! 0 -> U
  ! 1 -> dU/dRHO_0
  ! 2 -> dU/dRHO_1
  ! 3 -> d2U/(dRHO_0*dRHO_0)
  ! 4 -> d2U/(dRHO_1*dRHO_1)
  ! 5 -> d2U/(dRHO_0*dRHO_1)
  ! 6 -> dU/d(TAU_0)
  ! 7 -> dU/d(Delta RHO_0)
  !
  ! TESTED MATTHEMETICA<=>BIRUC & SCOTT; MATTHEMETICA<=>Module UNEDF (energy amplitudes only)
  ! ver#2: (Mario) Pairing included
  !  - set_functional_parameters(fname,lpr)
  !  - pairing incorporated into CpV0(0:1),CpV1(0:1)
  !    as public variables also serving two public amplitudes
  !     Urhorhopr(0:1,0)=CpV0(0:1)+CpV1(0:1)*rho(0)
  !     Urhorhopr(0:1,1)=CpV1(0:1)
  !     so, they can be used with appropriate values by the DME solver
  !  -need improvement later,
  !      currently HFBTHO uses CpV0(0:1), CpV0(0:1)  as before
  !      just substituting V0,V1 in pn-representation
  !      CpV0*(1-CpV1/0.16*rho_0)and this defines
  !      the default values in the module CpV0=V0,CpV1=1/2)
  !  -NAMELIST and input/output modified. RESERVED NAMES ARE:
  !      -namelist forbiden:
  !          'UNRDF'  - best UNEDF
  !          'SKYRME' - best SKYRME
  !      -namelist inforced but not for C-parameters (use_INM=F)
  !       or NM-parameters (use_INM=T) defined by the solver
  !          'FITS'
  !      -namelist inforced (one can overwrite all):
  !          'ANY OTHER NAME'
  !       i.e., the solver defines C-/NM- only using 'FITS'
  ! ver#1: (Mario) Complete rewrite consistent with HFBTHO
  !  -CB-LDA added
  !  -INM added
  !  -HFBTHO BENCHMARK: LN, ZR(110) prolate solution with SLY4,
  !   mixed pairing and tensor terms. Agreement with previouse
  !   implemetation to the last significant digit in the cases:
  !      - Standard Skyrme
  !      - LO+LDA
  !      - LO+CB-LDA
  !      - (NrNr=0,rDj=0), (rDr=0,jDr=0), 0.5(NrNr=-rDr,jDr=-rDj)
  !   -use_j2terms removed, i.e., in the SKYRME case CJ=0 removes all
  !    tensor terms, while in DME tensor terms are always present
  ! ver#0: (Marcus) Basic coding from scratch
  !   -DME(u) consistent with Mathematica numbers
  !   -including small 'u' approximation
  !--------------------------------------------------------------------------------------
  !
  ! === PUBLIC VARIABLES ===
  !
  ! Use pointers to prevent conflicts with UNEDF public variabes
  ! Example: Use UNEDF, pr=>my_pr, ipr=>my_ipr, Crho=>my_Crho ...
  !
  !--------------------------------------------------------------------------------------
  !Integer, Parameter, Public :: ipr=Kind(1)                                   ! to set the precision of the DFT solver
  !Integer, Parameter, Public :: pr=Kind(1.000D0)                              ! to set the precision of the DFT solver
  Logical, Public :: use_charge_density, use_cm_cor,use_DME3N_terms,   &
                     use_j2terms,use_full_cm_cor,use_INM,use_Namelist, &
                     Print_Namelist
  Integer(ipr), Public :: DMEorder,DMElda,use_TMR_pairing
  Real(pr), Public, Dimension(0:3,0:7) :: Urhorho,Urhotau,UrhoDrho,Unablarho  ! ph DME amplitudes
  Real(pr), Public, Dimension(0:3,0:7) :: UJnablarho,UrhonablaJ,UJJ
  Real(pr), Public, Dimension(0:3,0:7) :: Urhorhopr                           ! pp amplitudes
  Real(pr), Public, Dimension(0:1) :: UEnonstdr,UFnonstdr,URnonstdr           ! Other amplitudes
  Real(pr), Public :: hbzero,sigma,e2charg,CExPar                             ! hbr^2/2m, DD sigma, e^2 charge, coul.exch.
  Real(pr), Public, Dimension(0:1) :: Crho,Cdrho,Ctau,CrDr,CrdJ,CJ,CpV0,CpV1  ! basic coupling constants
  Real(pr), Public :: E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM,P_NM,KA_NM
  Real(pr), Public :: CHrho                                                  ! Crho(0) from the Hartree term in NM
  Real(pr), Public :: mpi,gA,fpi,c1,c3,c4,cd,ce,LambdaX
  Real(pr), PUBLIC :: t0s,t0a,drs,dra,ts,ta,t3alp,t3al0,t3alm,t324,alp,alm,wla0, &
                      wla1,TA7,TA8,TB7,TB8,tv1,tv2,tv3,tv4,tv5,tv6,ts1,ts2,t4o3
  Real(pr), PUBLIC :: t0_pub,t1_pub,t2_pub,t3_pub,x0_pub,x1_pub,x2_pub,x3_pub,&
                      b4_pub,b4p_pub,te_pub,to_pub
  !
  ! === PRIVATE VARIABLES ===
  !
  Real(pr), Private, Dimension(0:1) :: nuCrho,nuCdrho,nuCtau,nuCrDr  ! basic coupling constants in natural units
  Real(pr), Private, Dimension(0:1) :: nuCrdJ,nuCJ,nuCpV0,nuCpV1     !
  Real(pr), Private :: t0,t1,t2,t3,x0,x1,x2,x3,b4,b4p,te,to
  Real(pr), Private :: nuLambda,nufpi                                ! parameters associated to natural units
  Real(pr), Private, Dimension(0:1) :: Cnrho,CJdr                    ! hidden and always zero
  Integer(ipr), Private :: i_cut                                     ! dmeorder: -1=Standard Skyrme, 0=LO, 1=NLO, 2=N2LO
  Real(pr), Private :: Pi,eps                                        ! dmelda: 0=Kf-LDA, 1=CB-LDA
  Real(pr), Private :: kfconst,CK                                    ! (3Pi^2/2)^(1/3)
  Real(pr), Parameter, Private :: mevfm=197.30_pr;
  Real(pr), Private :: rho(0:1),tau(0:1),nrho2(0:1),lrho(0:1)
  Real(pr), Private :: mpi2,fpi2,fpi4,gA2,gA4,gA6,CHartree
  Real(pr), Private :: arhorho,brhorho,arhodrho,brhodrho,arhotau,brhotau,ajj,bjj,adrdr,bdrdr
  Real(pr), Private :: darhorho,dbrhorho,darhodrho,dbrhodrho,darhotau,dbrhotau,dajj,dbjj,dadrdr,dbdrdr
  Real(pr), Private :: ddarhodrho,ddbrhodrho,ddarhotau,ddbrhotau,ddarhorho,ddbrhorho
  Real(pr), Private :: hrho0rho0,hrho1rho1,hdr0dr0,hdr1dr1,hrho0Drho0,hrho1Drho0, &
       hrho1Drho1,hrho0tau0,hrho1tau0,hrho1tau1,hJ0dr0,hrho0DJ0,hJ1dr1,hrho1DJ1, &
       hJ0dr1,hrho1DJ0,hJ1dr0,hJ0J0,hJ0J1,hJ1J1
  Real(pr), Private :: dhrho0rho0,dhrho1rho1,dhdr0dr0,dhdr1dr1,dhrho0Drho0, &
       dhrho1Drho0,dhrho1Drho1,dhrho0tau0,dhrho1tau0,dhrho1tau1,dhJ0dr0,dhrho0DJ0, &
       dhJ1dr1,dhrho1DJ1,dhJ0dr1,dhrho1DJ0,dhJ1dr0,dhJ0J0,dhJ0J1,dhJ1J1
  Real(pr), Private :: ddhrho0rho0,ddhrho1rho1,ddhrho0Drho0,ddhrho1Drho0, &
       ddhrho1Drho1,ddhrho0tau0,ddhrho1tau0,ddhrho1tau1
  Real(pr), Private, Dimension(3,3,33) :: ctr0r0,ctr1r1,ctdr0dr0,ctdr1dr1, & ! coefficients for 3N part
       ctr0Dr0,ctr1Dr0,ctr1Dr1,ctr0t0,ctr1t0,ctr1t1,ctJ0dr0,ctr0dJ0,ctJ1dr1, &
       ctr1dJ1,ctJ0dr1,ctr1dJ0,ctJ1dr0,ctJ0J0,ctJ0J1,ctJ1J1
  Real(pr), Private :: u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12
  Real(pr), Private :: ual,lual,atu,asqu,asqu4
  Real(pr), Private :: ac2,ac3,acoord
  Parameter(acoord=0.50_pr,ac2=4.0_pr*(acoord**2-acoord+0.50_pr),ac3=2.0_pr*(acoord**2-acoord+0.50_pr))
  Character (30) :: FunctionalName
  !
  Real(pr), Private :: A1_1,A1_2,A1_3,A1_4,A1_5,b1_1,b1_2,b1_3,b1_4,b1_5
  Real(pr), Private :: A3_1,A3_2,A3_3,A3_4,A3_5,b3_1,b3_2,b3_3,b3_4,b3_5
  Real(pr), Private :: h0mpi6,h0mpi6c1,h0mpi6c3,h0mpi6NM,h0mpi6c1NM,h0mpi6c3NM
  !
  Namelist /UNEDF_NAMELIST/ FunctionalName, DMEorder, DMElda, use_INM, hbzero, use_TMR_pairing,   &
                            Crho, Cdrho, Ctau, CrDr, CrdJ, CJ, sigma, CpV0, CpV1, e2charg,        &
                            E_NM, K_NM, SMASS_NM, RHO_NM, ASS_NM, LASS_NM, VMASS_NM,              &
                            mpi, gA, fpi, c1, c3, c4, cd, ce, LambdaX,                            &
                            use_cm_cor, use_charge_density, use_DME3N_terms, use_j2terms, CExPar, &
                            Print_Namelist
 Contains
  !
  !===============================================================================================
  Subroutine calculate_U_parameters(rho0_in,rho1_in,tau0_in,tau1_in,laprho0,laprho1,nablarho0s,nablarho1s)
    Implicit None
    Real(pr), Intent(in) :: rho0_in,rho1_in,tau0_in,tau1_in
    Real(pr), Intent(in), Optional :: &
         nablarho0s,nablarho1s,laprho0,laprho1
    Integer(ipr) :: t,i,j,k,l
    Real(pr) :: u,du,ddu,dtu,dlu
    Real(pr) :: ph,aux,daux,ddaux
    Real(pr) :: y,dy,ddy,marc,dmarc,ddmarc,mlog,dmlog,ddmlog
    Real(pr) :: ucut,ucut3n
    !
    ucut=0.1_pr; ucut3n=0.6_pr
    !
    rho(0)=rho0_in; rho(1)=rho1_in;
    tau(0)=tau0_in; tau(1)=tau1_in;
    !
    lrho=0.0_pr; nrho2=0.0_pr;
    If (Present(laprho0)) lrho(0)=laprho0
    If (Present(laprho1)) lrho(1)=laprho1
    If (Present(nablarho0s)) nrho2(0)=nablarho0s
    If (Present(nablarho1s)) nrho2(1)=nablarho1s
    !
    arhorho=0.0_pr; darhorho=0.0_pr; ddarhorho=0.0_pr
    brhorho=0.0_pr; dbrhorho=0.0_pr; ddbrhorho=0.0_pr
    arhodrho=0.0_pr; darhodrho=0.0_pr; ddarhodrho=0.0_pr
    brhodrho=0.0_pr; dbrhodrho=0.0_pr; ddbrhodrho=0.0_pr
    arhotau=0.0_pr; darhotau=0.0_pr; ddarhotau=0.0_pr
    brhotau=0.0_pr; dbrhotau=0.0_pr; ddbrhotau=0.0_pr
    adrdr=0.0_pr; dadrdr=0.0_pr
    bdrdr=0.0_pr; dbdrdr=0.0_pr
    ajj=0.0_pr; dajj=0.0_pr
    bjj=0.0_pr; dbjj=0.0_pr
    !
    hrho0rho0=0.0_pr; hrho1rho1=0.0_pr; hdr0dr0=0.0_pr; hdr1dr1=0.0_pr
    hrho0Drho0=0.0_pr; hrho1Drho0=0.0_pr; hrho1Drho1=0.0_pr
    hrho0tau0=0.0_pr; hrho1tau0=0.0_pr; hrho1tau1=0.0_pr
    hJ0dr0=0.0_pr; hrho0DJ0=0.0_pr; hJ1dr1=0.0_pr; hrho1DJ1=0.0_pr
    hJ0dr1=0.0_pr; hrho1DJ0=0.0_pr; hJ1dr0=0.0_pr
    hJ0J0=0.0_pr; hJ0J1=0.0_pr; hJ1J1=0.0_pr
    dhrho0rho0=0.0_pr; dhrho1rho1=0.0_pr; dhdr0dr0=0.0_pr; dhdr1dr1=0.0_pr
    dhrho0Drho0=0.0_pr; dhrho1Drho0=0.0_pr; dhrho1Drho1=0.0_pr
    dhrho0tau0=0.0_pr; dhrho1tau0=0.0_pr; dhrho1tau1=0.0_pr
    dhJ0dr0=0.0_pr; dhrho0DJ0=0.0_pr; dhJ1dr1=0.0_pr; dhrho1DJ1=0.0_pr
    dhJ0dr1=0.0_pr; dhrho1DJ0=0.0_pr; dhJ1dr0=0.0_pr
    dhJ0J0=0.0_pr; dhJ0J1=0.0_pr; dhJ1J1=0.0_pr
    ddhrho0rho0=0.0_pr; ddhrho1rho1=0.0_pr
    ddhrho0Drho0=0.0_pr; ddhrho1Drho0=0.0_pr; ddhrho1Drho1=0.0_pr
    ddhrho0tau0=0.0_pr; ddhrho1tau0=0.0_pr; ddhrho1tau1=0.0_pr
    !
    u=0.0_pr; du=0.0_pr; ddu=0.0_pr; dtu=0.0_pr; dlu=0.0_pr
    !
    Urhorho=0.0_pr   ; Urhotau=0.0_pr
    UrhoDrho=0.0_pr  ; Unablarho=0.0_pr
    UJnablarho=0.0_pr; UrhonablaJ=0.0_pr
    Urhorhopr=0.0_pr ; UJJ=0.0_pr
    UEnonstdr=0.0_pr ; UFnonstdr=0.0_pr ; URnonstdr=0.0_pr
    !
    ! Notations for Uamplitudes(0:3,0:7)
    ! t for Uamplitudes(t,*)
    ! 0 -> 0,0
    ! 1 -> 1,1
    ! 2 -> 0,1
    ! 3 -> 1,0
    ! n for Uamplitudes(*,n)
    ! 0 -> U
    ! 1 -> dU/dRHO_0
    ! 2 -> dU/dRHO_1
    ! 3 -> d2U/(dRHO_0*dRHO_0)
    ! 4 -> d2U/(dRHO_1*dRHO_1)
    ! 5 -> d2U/(dRHO_0*dRHO_1)
    ! 6 -> dU/d(TAU_0)
    ! 7 -> dU/d(Delta RHO_0)
    !
    !! 2N terms
    Do t=0,1
       ph=1.0_pr
       If(t.Eq.1) ph=-1.0_pr
       Urhorho(t,0)=Crho(t)+Cdrho(t)*rho(0)**sigma &
            +0.50_pr*(arhorho+ph*brhorho)*mevfm
       Urhotau(t,0)=Ctau(t)+0.50_pr*(arhotau+ph*brhotau)*mevfm
       UrhoDrho(t,0)=Crdr(t)+ac2*0.50_pr*(arhoDrho+ph*brhoDrho)*mevfm
       UJJ(t,0)=CJ(t)+0.50_pr*(ajj+ph*bjj)*mevfm
       Unablarho(t,0)=Cnrho(t)+0.50_pr*(adrdr+ph*bdrdr)*mevfm
       UrhonablaJ(t,0)=Crdj(t)
       UJnablarho(t,0)=Cjdr(t)

       Urhorho(t,1)=sigma*Cdrho(t)*(rho(0)**sigma)/(rho(0)+eps) &
            +0.50_pr*(darhorho+ph*dbrhorho)*du*mevfm
       Urhotau(t,1)=0.50_pr*(darhotau+ph*dbrhotau)*du*mevfm
       UrhoDrho(t,1)=ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*du*mevfm
       UJJ(t,1)=0.50_pr*(dajj+ph*dbjj)*du*mevfm
       Unablarho(t,1)=0.50_pr*(dadrdr+ph*dbdrdr)*du*mevfm

       Urhorho(t,6)=0.50_pr*(darhorho+ph*dbrhorho)*dtu*mevfm
       Urhotau(t,6)=0.50_pr*(darhotau+ph*dbrhotau)*dtu*mevfm
       UrhoDrho(t,6)=ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*dtu*mevfm
       UJJ(t,6)=0.50_pr*(dajj+ph*dbjj)*dtu*mevfm
       Unablarho(t,6)=0.50_pr*(dadrdr+ph*dbdrdr)*dtu*mevfm

       Urhorho(t,7)=0.50_pr*(darhorho+ph*dbrhorho)*dlu*mevfm
       Urhotau(t,7)=0.50_pr*(darhotau+ph*dbrhotau)*dlu*mevfm
       UrhoDrho(t,7)=ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*dlu*mevfm
       UJJ(t,7)=0.50_pr*(dajj+ph*dbjj)*dlu*mevfm
       Unablarho(t,7)=0.50_pr*(dadrdr+ph*dbdrdr)*dlu*mevfm

       Urhorho(t,3)=sigma*(sigma-1.0_pr)*Cdrho(t)*(rho(0)**sigma)/(rho(0)**2+eps) &
            +0.50_pr*(darhorho+ph*dbrhorho)*ddu*mevfm &
            +0.50_pr*(ddarhorho+ph*ddbrhorho)*du*du*mevfm
       Urhotau(t,3)=0.50_pr*(darhotau+ph*dbrhotau)*ddu*mevfm &
            +0.50_pr*(ddarhotau+ph*ddbrhotau)*du*du*mevfm
       UrhoDrho(t,3)=ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*ddu*mevfm &
            +ac2*0.50_pr*(ddarhoDrho+ph*ddbrhoDrho)*du*du*mevfm

    End Do
    Urhorhopr(0,0) = CpV0(0)*(1.0_pr-CpV1(0)*rho(0)/0.16_pr)          &
         +CpV0(1)*(1.0_pr-CpV1(1)*rho(0)/0.16_pr)
    Urhorhopr(1,0) = CpV0(0)*(1.0_pr-CpV1(0)*rho(0)/0.16_pr)          &
         +CpV0(1)*(1.0_pr-CpV1(1)*rho(0)/0.16_pr)
    Urhorhopr(2,0) = (CpV0(0)*(1.0_pr-CpV1(0)*rho(0)/0.16_pr)         &
         -CpV0(1)*(1.0_pr-CpV1(1)*rho(0)/0.16_pr))*2.0_pr
    Urhorhopr(0,1) = (-CpV0(0)*CpV1(0)-CpV0(1)*CpV1(1))/0.16_pr
    Urhorhopr(1,1) = (-CpV0(0)*CpV1(0)-CpV0(1)*CpV1(1))/0.16_pr
    Urhorhopr(2,1) = 2.0_pr*(-CpV0(0)*CpV1(0)+CpV0(1)*CpV1(1))/0.16_pr
    Urhorhopr=Urhorhopr/16.0_pr
    !
    If (.Not.use_j2terms) Then
       UJJ=0.0_pr
    End If
    !
  End Subroutine calculate_U_parameters
  !===============================================================================================
  !
  !===============================================================================================
  Subroutine read_UNEDF_NAMELIST(fname,noForces)
    Use HFBTHO_utilities, Only: lout
    !--------------------------------------------------------------------------------
    ! RESERVED NAMES ARE:
    !  -namelist forbiden:
    !          'UNEDF'  - best UNEDF
    !          'SKYRME' - best SKYRME
    !  -namelist inforced but not for C-parameters (use_INM=F)
    !   or NM-parameters (use_INM=T) defined by the solver
    !          'FITS'
    !  -namelist inforced (one can overwrite all):
    !          'ANY OTHER NAME'
    ! i.e., the DME solver defines C-/NM- only using 'FITS'
    !--------------------------------------------------------------------------------
    Implicit None
    Character (30), Intent(inout) :: fname
    Character (30) :: inforcedname
    Logical        :: regularization
    Integer(ipr)   :: ios,lnamelist=16,noForces
    !
    ! parameters
    eps     = Spacing(1.0_pr)
    Pi      = 4.0_pr*Atan(1.0_pr)
    kfconst =(1.50_pr*Pi**2)**(1.0_pr/3.0_pr)    ! (3Pi^2/2)^(1/3)
    CK      = 3.0_pr/5.0_pr*kfconst**2
    !
    use_Namelist=.True.
    Do
       !---------------------------------------------------------------------
       ! Some default values for all cases
       !---------------------------------------------------------------------
       Print_Namelist=.False.
       FunctionalName="Bla-Bla"
       ! kind of the functional
       use_INM            = .False.
       use_DME3N_terms    = .False.
       use_charge_density = .False.
       regularization     = .False.
       use_cm_cor         = .False.
       use_full_cm_cor    = .False.
       use_j2terms        = .False.
       use_TMR_pairing    =  0
       DMEorder           = -1
       DMElda             =  0
       ! Coupling constants: ph channel
       Crho(0)  = -727.0933239596374733_pr; Crho(1)  =  474.8709969984467989_pr
       CDrho(0) =  612.1037411660222460_pr; CDrho(1) = -705.7204872069220301_pr
       Ctau(0)  =   33.8846741217252401_pr; Ctau(1)  =   32.4047409594248919_pr
       CrDr(0)  =  -76.9962031249999939_pr; CrDr(1)  =   15.6571351249999999_pr
       CrdJ(0)  =  -92.2500000000000000_pr; CrdJ(1)  =  -30.7500000000000000_pr
       CJ(0)    =   17.2096115000000012_pr; CJ(1)    =   64.5758124999999978_pr
       Cnrho    =    0.0000000000000000_pr; CJdr     =    0.0000000000000000_pr
       ! Coupling constants: pp channel
       CpV0     = -258.2000000000000000_pr; CpV1     =    0.5000000000000000_pr
       ! Various
       sigma    =    0.3062227576210547_pr;
       hbzero   =   20.7355300000000007_pr;
       e2charg  =    1.4399784085965135_pr ; CExPar = 1.0_pr
       ! DME
       mpi=   138.03_pr/197.3_pr; fpi = 92.4_pr/197.3_pr; gA = 1.29_pr
       c1 =    -0.81_pr/1000.0_pr*197.3_pr
       c3 =    -3.40_pr/1000.0_pr*197.3_pr
       c4 =     3.40_pr/1000.0_pr*197.3_pr
       cd = -2062.00_pr/1000.0_pr
       ce =  -625.00_pr/1000.0_pr
       ! Natural units
       LambdaX  = 700.0_pr/197.3_pr; nuLambda = 700.0_pr; nufpi = 93.0_pr
       ! Nuclear matter
       E_NM     = -15.972149141444596410_pr; RHO_NM   =  0.159538756711733398_pr
       K_NM     = 229.900964482603626493_pr; SMASS_NM =  1.439546988976078357_pr
       ASS_NM   =  32.004302815052007247_pr; LASS_NM  = 45.961751480461579433_pr
       VMASS_NM =   1.249838547196253424_pr
       !---------------------------------------------------------------------
       ! Select the functional: start with interaction
       !---------------------------------------------------------------------
       noForces=0 ! No forces to start with
       Call skforce(fname,noForces)
       !
       If (noForces.Eq.1) Then
           inforcedname='FORCE'
           use_Namelist=.False.
       Else
           FUNCTIONAL: Select Case (Trim(fname))
           Case ('UNE0')
              inforcedname='UNE0'
              use_Namelist=.False.
              ! kind of the functional
              use_INM    = .True.
              use_cm_cor = .True.
              ! Surface coefficients
              CrDr(0)  =  -55.260600000000000_pr
              CrDr(1)  =  -55.622600000000000_pr
              CpV0(0)  = -170.374000000000000_pr
              CpV0(1)  = -199.202000000000000_pr
              CrdJ(0)  =  -79.530800000000000_pr
              CrdJ(1)  =   45.630200000000000_pr
              CJ(0)    =    0.000000000000000_pr
              CJ(1)    =    0.000000000000000_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.160526000000000_pr
              E_NM     =  -16.055900000000000_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  230.000000000000000_pr
              ASS_NM   =   30.542900000000000_pr
              LASS_NM  =   45.080400000000000_pr
              SMASS_NM =    0.900000000000000_pr
              VMASS_NM =    1.249838000000000_pr
           Case ('UNE1')
              inforcedname='UNE1'
              use_Namelist=.False.
              ! kind of the functional
              use_INM  = .True.
              ! Surface coefficients
              CrDr(0)  =  -45.135131022237300_pr
              CrDr(1)  = -145.382167908057000_pr
              CpV0(0)  = -186.065399575124000_pr
              CpV0(1)  = -206.579593890106000_pr
              CrdJ(0)  =  -74.026333176459900_pr
              CrdJ(1)  =  -35.658261114791700_pr
              CJ(0)    =    0.000000000000000_pr
              CJ(1)    =    0.000000000000000_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.158706769332587_pr
              E_NM     =  -15.800000000000000_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  220.000000000000000_pr
              ASS_NM   =   28.986789057772100_pr
              LASS_NM  =   40.004790480413600_pr
              SMASS_NM =    0.992423332283364_pr
              VMASS_NM =    1.249838574232270_pr
           Case ('UNE2')
              inforcedname='UNE2'
              use_Namelist=.False.
              ! kind of the functional
              use_INM     = .True.
              use_j2terms = .True.
              ! Surface coefficients
              CrDr(0)  =  -46.831409147060600_pr
              CrDr(1)  = -113.163790795259000_pr
              CpV0(0)  = -208.889001962571000_pr
              CpV0(1)  = -230.329984038628000_pr
              CrdJ(0)  =  -64.308862415783800_pr
              CrdJ(1)  =  -38.650194685135500_pr
              CJ(0)    =  -54.433363597372100_pr
              CJ(1)    =  -65.903031044593800_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.156310622197074_pr
              E_NM     =  -15.800000000000000_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  239.929568022437000_pr
              ASS_NM   =   29.131006470773700_pr
              LASS_NM  =   40.000000000000000_pr
              SMASS_NM =    1.073763804147980_pr
              VMASS_NM =    1.249838574232270_pr
           Case ('HFB1')
              inforcedname='HFB1'
              use_Namelist=.False.
              ! kind of the functional
              use_INM  = .True.
              ! Surface coefficients
              CrDr(0)  =  -45.599763_pr
              CrDr(1)  = -143.935074_pr
              CpV0(0)  = -234.380010_pr
              CpV0(1)  = -260.437001_pr
              CrdJ(0)  =  -73.946388_pr
              CrdJ(1)  =  -51.912548_pr
              CJ(0)    =    0.000000_pr
              CJ(1)    =    0.000000_pr
              CExPar   =    1.000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.156247_pr
              E_NM     =  -15.800000_pr
              P_NM     =    0.000000_pr
              K_NM     =  244.668204_pr
              ASS_NM   =   28.668204_pr
              LASS_NM  =   40.109081_pr
              SMASS_NM =    1.067970_pr
              VMASS_NM =    1.249838574232270_pr
           Case ('N0LO')
              inforcedname='N0LO'
              use_Namelist=.False.
              ! kind of the functional
              use_INM         = .True.
              use_j2terms     = .False.
              use_DME3N_terms = .False.
              DMEorder        = 0
              ! Surface coefficients
              CrDr(0)  =  -67.437_pr
              CrDr(1)  =   21.551_pr
              CpV0(0)  = -241.203_pr
              CpV0(1)  = -252.818_pr
              CrdJ(0)  =  -95.451_pr
              CrdJ(1)  =  -65.906_pr
              CJ(0)    =    0.000_pr
              CJ(1)    =    0.000_pr
              CExPar   =    1.000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.1595_pr
              E_NM     =  -15.9700_pr
              P_NM     =    0.0000_pr
              K_NM     =  229.9000_pr
              ASS_NM   =   32.0000_pr
              LASS_NM  =   45.9600_pr
              SMASS_NM =    1.4400_pr
              VMASS_NM =    1.2500_pr
           Case ('N1LO')
              inforcedname='N1LO'
              use_Namelist=.False.
              ! kind of the functional
              use_INM         = .True.
              use_j2terms     = .False.
              use_DME3N_terms = .False.
              DMEorder        = 1
              ! Surface coefficients
              CrDr(0)  =  -63.996_pr
              CrDr(1)  =    -9.276_pr
              CpV0(0)  = -241.484_pr
              CpV0(1)  = -252.222
              CrdJ(0)  =  -95.463_pr
              CrdJ(1)  =  -60.800_pr
              CJ(0)    =    0.000_pr
              CJ(1)    =    0.000_pr
              CExPar   =    1.000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.1595_pr
              E_NM     =  -15.9700_pr
              P_NM     =    0.0000_pr
              K_NM     =  229.9000_pr
              ASS_NM   =   32.0000_pr
              LASS_NM  =   45.9600_pr
              SMASS_NM =    1.4400_pr
              VMASS_NM =    1.2500_pr
           Case ('N2LO')
              inforcedname='N2LO'
              use_Namelist=.False.
              ! kind of the functional
              use_INM         = .True.
              use_j2terms     = .False.
              use_DME3N_terms = .True.
              DMEorder        = 2
              ! Surface coefficients
              CrDr(0)  = -197.132_pr
              CrDr(1)  =  -12.503_pr
              CpV0(0)  = -272.164_pr
              CpV0(1)  = -193.188_pr
              CrdJ(0)  = -193.188_pr
              CrdJ(1)  =   37.790_pr
              CJ(0)    =    0.000_pr
              CJ(1)    =    0.000_pr
              CExPar   =    1.000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.1595_pr
              E_NM     =  -15.9700_pr
              P_NM     =    0.0000_pr
              K_NM     =  229.9000_pr
              ASS_NM   =   32.0000_pr
              LASS_NM  =   45.9600_pr
              SMASS_NM =    1.4400_pr
              VMASS_NM =    1.2500_pr
           Case default
              inforcedname=fname
              use_Namelist=.True.
           End Select FUNCTIONAL
       End If
       !---------------------------------------------------------------------
       ! Exit loop condition
       !---------------------------------------------------------------------
       If(.Not.use_Namelist) Exit
       !---------------------------------------------------------------------
       ! Read namelists
       !---------------------------------------------------------------------
       Open(lnamelist,file='UNEDF_NAMELIST.DAT',DELIM='APOSTROPHE') ! 'QUOTE'
       Read(UNIT=lnamelist,NML=UNEDF_NAMELIST,iostat=ios)
       If(ios.Ne.0) Then
          ! WRong entry within UNEDF_NAMELIST.DAT file
          Write(lout,'(1X,/,A)') 'ATTENTION: WRONG INPUT!'
          Write(lout,*) 'THE INPUT DATA WITH LABEL FUNCTIONALNAME=''',Trim(INFORCEDNAME),''''
          Write(lout,*) 'INSIDE THE UNEDF_NAMELIST.DAT FILE IS WRONG.'
          Write(lout,*) 'PLESE CORECT AND TRY AGAIN!'
          Stop 'PROGRAM STOP IN read_UNEDF_NAMELIST'
       End If
       Close(lnamelist)
       If(Trim(FunctionalName).Eq.Trim(inforcedname)) Exit
    End Do
    !---------------------------------------------------------------------
    ! See what the namelists modified
    !---------------------------------------------------------------------
    INFORCED_FUNCTIONAL: Select Case (Trim(inforcedname))
    Case ("FORCE")
       FunctionalName='FORCE'
    Case ("UNE0")
       FunctionalName='UNE0'
    Case ("UNE1")
       FunctionalName='UNE1'
    Case ("UNE2")
       FunctionalName='UNE2'
    Case ("HFB1")
       FunctionalName='HFB1'
    Case ("N0LO")
       FunctionalName='N0LO'
    Case ("N1LO")
       FunctionalName='N1LO'
    Case ("N2LO")
       FunctionalName='N2LO'
    Case default
       ! Missing entry within hfbtho_NAMELIST.dat file
       If(Trim(FunctionalName).Ne.Trim(inforcedname)) Then
          Write(lout,'(1X,/,A)') 'ATTENTION: MISSING INPUT!'
          Write(lout,*) 'THE INPUT DATA WITH LABEL FUNCTIONALNAME=''',Trim(INFORCEDNAME),''''
          Write(lout,*) 'IS MISSING INSIDE THE UNEDF_NAMELIST.DAT FILE.'
          Write(lout,*) 'PLESE CORECT AND TRY AGAIN!'
          Stop 'PROGRAM STOP IN SET_FUNCTIONAL_PARAMETERS'
       End If
    End Select INFORCED_FUNCTIONAL
    !
  End Subroutine read_UNEDF_NAMELIST
  !===============================================================================================
  !
  !===============================================================================================
  Subroutine skforce(fname,noForces)
    !---------------------------------------------------------------------
    ! Set up Pairing & Skyrme force parameters and their combinations
    !---------------------------------------------------------------------
    Use HFBTHO_utilities, Only: lout
    Implicit None
    Integer(ipr) :: noForces
    Real(pr) :: A,wls
    Real(pr) :: zero,one,two,three,four,five,six,seven,eight,nine
    Real(pr) :: half,pp16,pp24
    Character (30), Intent(inout) :: fname
    !
    zero = 0.0_pr; one = 1.0_pr; two = 2.0_pr; three = 3.0_pr; four = 4.0_pr
    five = 5.0_pr; six = 6.0_pr; seven = 7.0_pr; eight = 8.0_pr; nine = 9.0_pr
    half = 0.5_pr; pp16 = 16.0_pr; pp24 = 24.0_pr
    !
    ! Default for all forces if not modified
    hbzero = 1.0_pr/0.04823_pr ! DMSHB0=1/hbzero
    sigma = one
    t0 = zero; x0 = zero
    t1 = zero; x1 = zero
    t2 = zero; x2 = one
    t3 = zero; x3 = one
    wls= zero; b4 = wls/two; b4p=wls/two
    te = zero; to = zero
    CExPar=1.0_pr
    !
    noForces=0 ! No forces at all
    !
    INTERACTION: Select Case (Trim(fname))
    !---------------------------------------------------------------------
    ! SIII force, Beiner et al., NPA238 (1975) 29
    !---------------------------------------------------------------------
    Case ('SIII')
        ! ph-Force
        noForces=1
        use_cm_cor = .True.
        hbzero = 20.73533_pr
        t0 = -.1128750d+04; x0 = +0.4500000_pr
        t1 = +.3950000d+03; x1 = +0.0000000_pr
        t2 = -.9500000d+02; x2 = +0.0000000_pr
        t3 = +.1400000d+05; x3 = +1.0000000_pr
        wls= +.1200000d+03; sigma = one
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -265.2500_pr, -340.0625_pr /)
    !---------------------------------------------------------------------
    ! SKM* forces
    !---------------------------------------------------------------------
    Case ('SKM*')
        ! ph-Force
        noForces=1
        use_cm_cor = .True.
        hbzero = 20.73_pr
        t0 = -.2645000d+04; x0 = +.0900000_pr
        t1 = +.4100000d+03; x1 = +.0000000_pr
        t2 = -.1350000d+03; x2 = +.0000000_pr
        t3 = +.1559500d+05; x3 = +.0000000_pr
        wls= +.1300000d+03; sigma = one/six
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -265.2500_pr, -340.0625_pr /)
    !---------------------------------------------------------------------
    ! SKP force, Dobaczewski et al., NPA422 (1984) 103
    !---------------------------------------------------------------------
    Case ('SKP')
        ! ph-Force
        noForces=1
        use_cm_cor  = .True.
        use_j2terms = .True.
        hbzero = 20.730_pr
        t0 =-0.2931696d+04; x0 = 0.2921515_pr
        t1 = 0.3206182d+03; x1 = 0.6531765_pr
        t2 =-0.3374091d+03; x2 =-0.5373230_pr
        t3 = 0.1870896d+05; x3 = 0.1810269_pr
        wls= 0.1000000d+03; sigma=one/six
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -265.2500_pr, -340.0625_pr /)
    !---------------------------------------------------------------------
    ! SLY4 force
    !---------------------------------------------------------------------
    Case ('SLY4')
        ! ph-Force
        noForces=1
        use_cm_cor = .True.
        hbzero = 20.735530_pr
        t0 =-0.2488913d+04; x0 = 0.8340000_pr
        t1 = 0.4868180d+03; x1 =-0.3440000_pr
        t2 =-0.5463950d+03; x2 =-1.0000000_pr
        t3 = 0.1377700d+05; x3 = 1.3540000_pr
        wls= 0.1230000d+03; sigma=one/six
        b4 = wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -325.2500_pr, -340.0625_pr /) ! HFB
    !---------------------------------------------------------------------
    ! SLY5 force
    !---------------------------------------------------------------------
    Case ('SLY5')
        ! ph-Force
        noForces=1
        use_cm_cor  = .True.
        use_j2terms = .True.
        hbzero = 20.73553_pr
        t0 =-0.2483450d+04; x0 = 0.7760000_pr
        t1 = 0.4842300d+03; x1 =-0.3170000_pr
        t2 =-0.5566900d+03; x2 =-1.0000000_pr
        t3 = 0.1375700d+05; x3 = 1.2630000_pr
        wls= 0.1250000d+03; sigma=one/six
        b4 = wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -291.5000_pr, -297.7402_pr /) ! HFB
    !---------------------------------------------------------------------
    ! SLY6 forces
    !---------------------------------------------------------------------
    Case ('SLY6')
        ! ph-Force
        noForces=1
        use_cm_cor      = .True.
        use_full_cm_cor = .True.
        hbzero = 20.73553_pr
        t0 =-0.2479500d+04; x0 = 0.8250000_pr
        t1 = 0.4621800d+03; x1 =-0.4650000_pr
        t2 =-0.4486100d+03; x2 =-1.0000000_pr
        t3 = 0.1367300d+05; x3 = 1.3550000_pr
        wls= 0.1220000d+03; sigma=one/six
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -291.5000_pr, -297.7402_pr /) ! HFB
    !---------------------------------------------------------------------
    ! SLY6 forces
    !---------------------------------------------------------------------
    Case ('SLY7')
        ! ph-Force
        noForces=1
        use_cm_cor      = .True.
        use_j2terms     = .True.
        use_full_cm_cor = .True.
        hbzero = 20.73553_pr
        t0 =-0.2480800d+04; x0 = 0.8480000_pr
        t1 = 0.4612900d+03; x1 =-0.4920000_pr
        t2 =-0.4339300d+03; x2 =-1.0000000_pr
        t3 = 0.1366900d+05; x3 = 1.3930000_pr
        wls= 0.1250000d+03; sigma=one/six
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -291.5000_pr, -297.7402_pr /) ! HFB
    !!!!
    !!!! THE PARAMETERIZATIONS OF THE LYON FAMILY OF SKYRME FORCES COMMENTED
    !!!! OUT BELOW IS COMPATIBLE WITH VERSIONS OF HFODD >= 2.52I. FOR COMPA-
    !!!! RISONS WITH  OLDER VERSIONS OF HFODD PLEASE USE THE DEFAULT (UNCOM-
    !!!! MENTED) PARAMETERIZATIONS ABOVE,  WHICH  WAS USED  FOR ALL RELEVANT
    !!!! BENCHMARKS IN THE CPC ARTICLE DESCRIBING HFBTHO.
    !!!!
    !!!!---------------------------------------------------------------------
    !!!! SLY4 force
    !!!!---------------------------------------------------------------------
    !!!Case ('SLY4')
    !!!    ! ph-Force
    !!!    noForces=1
    !!!    use_cm_cor = .True.
    !!!    hbzero = 20.73553_pr
    !!!    t0 =-0.248891d+04; x0 = 0.834_pr
    !!!    t1 = 0.486820d+03; x1 =-0.344_pr
    !!!    t2 =-0.546390d+03; x2 =-1.000_pr
    !!!    t3 = 0.137770d+05; x3 = 1.354_pr
    !!!    wls= 0.123000d+03; sigma=one/six
    !!!    b4 = wls/two; b4p=wls/two
    !!!    ! pp-Forces
    !!!    CpV1=0.50_pr
    !!!    CpV0=(/ -325.2500_pr, -340.0625_pr /) ! HFB
    !!!!---------------------------------------------------------------------
    !!!! SLY5 force
    !!!!---------------------------------------------------------------------
    !!!Case ('SLY5')
    !!!    ! ph-Force
    !!!    noForces=1
    !!!    use_cm_cor  = .True.
    !!!    use_j2terms = .True.
    !!!    hbzero = 20.73553_pr
    !!!    t0 =-0.248488d+04; x0 = 0.778_pr
    !!!    t1 = 0.483130d+03; x1 =-0.328_pr
    !!!    t2 =-0.549400d+03; x2 =-1.000_pr
    !!!    t3 = 0.137630d+05; x3 = 1.267_pr
    !!!    wls= 0.126000d+03; sigma=one/six
    !!!    b4 = wls/two; b4p=wls/two
    !!!    ! pp-Forces
    !!!    CpV1=0.50_pr
    !!!    CpV0=(/ -291.5000_pr, -297.7402_pr /) ! HFB
    !!!!---------------------------------------------------------------------
    !!!! SLY6 forces
    !!!!---------------------------------------------------------------------
    !!!Case ('SLY6')
    !!!    ! ph-Force
    !!!    noForces=1
    !!!    use_cm_cor      = .True.
    !!!    use_full_cm_cor = .True.
    !!!    hbzero = 20.73553_pr
    !!!    t0 =-0.24795d+04; x0 = 0.825_pr
    !!!    t1 = 0.46218d+03; x1 =-0.465_pr
    !!!    t2 =-0.44861d+03; x2 =-1.000_pr
    !!!    t3 = 0.13673d+05; x3 = 1.355_pr
    !!!    wls= 0.12200d+03; sigma=one/six
    !!!    b4=wls/two; b4p=wls/two
    !!!    ! pp-Forces
    !!!    CpV1=0.50_pr
    !!!    CpV0=(/ -291.5000_pr, -297.7402_pr /) ! HFB
    !!!!---------------------------------------------------------------------
    !!!! SLY6 forces
    !!!!---------------------------------------------------------------------
    !!!Case ('SLY7')
    !!!    ! ph-Force
    !!!    noForces=1
    !!!    use_cm_cor      = .True.
    !!!    use_j2terms     = .True.
    !!!    use_full_cm_cor = .True.
    !!!    hbzero = 20.73553_pr
    !!!    t0 =-0.248241d+04; x0 = 0.846_pr
    !!!    t1 = 0.457970d+03; x1 =-0.511_pr
    !!!    t2 =-0.419850d+03; x2 =-1.000_pr
    !!!    t3 = 0.136770d+05; x3 = 1.391_pr
    !!!    wls= 0.126000d+03; sigma=one/six
    !!!    b4=wls/two; b4p=wls/two
    !!!    ! pp-Forces
    !!!    CpV1=0.50_pr
    !!!    CpV0=(/ -291.5000_pr, -297.7402_pr /) ! HFB
    !---------------------------------------------------------------------
    ! SKI3 force, P.G.-Reinhard et al. Nucl. Phys. A584 (1995) 467-488
    !---------------------------------------------------------------------
    Case ('SKI3')
        ! ph-Force
        noForces=1
        use_cm_cor      = .True.
        use_full_cm_cor = .True.
        hbzero = 20.7525_pr
        t0 =-0.176288d+04; x0 = 0.30830_pr
        t1 = 0.561608d+03; x1 =-1.17220_pr
        t2 =-0.227090d+03; x2 =-1.09070_pr
        t3 = 0.810620d+04; x3 = 1.29260_pr
        sigma=one/four
        b4 = 94.254_pr; b4p=zero
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -357.2324_pr, -388.5625_pr /)
    !---------------------------------------------------------------------
    ! SKO forces
    !---------------------------------------------------------------------
    Case ('SKO')
        ! ph-Force
        noForces=1
        use_cm_cor      = .True.
        use_full_cm_cor = .True.
        hbzero = 20.735530_pr
        t0 =-0.21036530d+04; x0 = -0.2107010_pr
        t1 = 0.30335200d+03; x1 = -2.8107520_pr
        t2 = 0.79167400d+03; x2 = -1.4615950_pr
        t3 = 0.13553252d+05; x3 = -0.4298810_pr
        wls= 0.12300000d+03; sigma=one/four
        b4 = 0.17657800d+03; b4p=-0.1987490d+03
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -259.0391_pr, -274.8433_pr /)
    !---------------------------------------------------------------------
    ! SKX forces, A.Brown; Phys.Rev. C58 (1998) 220
    !---------------------------------------------------------------------
    Case ('SKX')
        ! ph-Force
        noForces=1
        use_cm_cor  = .True.
        use_j2terms = .True.
        hbzero = 20.73_pr
        t0 = -1445.300_pr; x0 = 0.340_pr
        t1 =   246.900_pr; x1 = 0.580_pr
        t2 =  -131.800_pr; x2 = 0.127_pr
        t3 = 12103.900_pr; x3 = 0.030_pr
        sigma=one/two
        b4 = 0.0743d+03; b4p=zero
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -259.0391_pr, -274.8433_pr /)
    !---------------------------------------------------------------------
    ! HFB9 forces
    !---------------------------------------------------------------------
    Case ('HFB9')
        ! ph-Force
        noForces=1
        use_cm_cor  = .True.
        use_j2terms = .True.
        hbzero = 20.73553_pr
        t0 =-0.20439180d+04; x0 = 0.5149210_pr
        t1 = 0.41159870d+03; x1 =-0.9537990_pr
        t2 =-0.19418860d+03; x2 =-0.3322490_pr
        t3 = 0.12497170d+05; x3 = 0.8994350_pr
        wls= 0.14990000d+03; sigma=one/four
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -263.5000_pr, -274.9668_pr /)
    !---------------------------------------------------------------------
    ! Default
    !---------------------------------------------------------------------
    Case default
        Write(lout,'("No Skyrme interaction defined in routine skforce()")')
    End Select INTERACTION
    !
    If (noForces.Eq.1) Then
        ! obtain coupling constants
        Call C_from_t()
        ! Frequent combinations entering the energy
        tv1   =  t0*(one+half*x0)*half;    tv2 = t0*(x0+half)*half
        tv3   =  t3*(one+half*x3)/12.0_pr; tv4 = t3*(x3+half)/12.0_pr
        tv5   = (t1*(one+half*x1)+t2*(one+half*x2))/four
        tv6   = (t2*(half+x2)-t1*(half+x1))/four
        ts1   = (t2*(one+half*x2)-three*t1*(one+half*x1))/pp16
        ts2   = (t1*(half+x1)*three+t2*(half+x2))/pp16
        t4o3  =  four/three; t324 = t3/pp24
        ! Frequent combinations entering the potential
        t0s   =  t0*(one-x0)*half; t0a = t0*(one+x0*half)
        drs   = (t2*(one+x2)-t1*(one-x1))*three/pp16
        dra   = (t2*(one+half*x2)-three*t1*(one+half*x1))/eight
        ts    = (t1*(one-x1) + three*t2*(one+x2))/eight
        ta    = (t1*(one+half*x1) + t2*(one+half*x2))/four
        t3alp = t3*(two+sigma)*(two+x3)/pp24
        t3al0 = t3*(x3+half)/six; t3alm = t3*sigma*(one+two*x3)/pp24
        alp   = one + sigma; alm = sigma - one
        wla0  = CrdJ(0)+CrdJ(1); wla1  = CrdJ(0)-CrdJ(1);
        TA7   = zero; TA8 = zero
        If(use_j2terms) Then
          TA7=(T1*(ONE-X1)-T2*(ONE+X2))/eight + five*to/four
          TA8=-(T1*X1+T2*X2)/four             + five*(te+to)/four
        End If
        TB7 = TA7; TB8 = TA8*half
    End If
    !
    Return
  End Subroutine skforce
  !===============================================================================================
  !
  !===============================================================================================
  Subroutine set_functional_parameters(fname,lpr)
    !--------------------------------------------------------------------------------
    ! set functional parameters
    !--------------------------------------------------------------------------------
    Implicit None
    Logical, Intent(in) :: lpr
    Character (30), Intent(inout) :: fname
    Logical :: regularization
    Integer(ipr), Parameter :: lin=15
    !
    ! parameters
    FunctionalName=fname
    eps=Spacing(1.0_pr)
    Pi=4.0_pr*Atan(1.0_pr)
    kfconst=(1.50_pr*Pi**2)**(1.0_pr/3.0_pr)    ! (3Pi^2/2)^(1/3)
    CK=3.0_pr/5.0_pr*kfconst**2
    nuLambda=700.0_pr ; nufpi = 93.0_pr
    !
    Call Make_Parameter_Free_Useful_Combinations()
    !
    ! exact Hartree CHrho from INM
    CHrho=0.0_pr; !!!!If (dmeorder.eq.3) Call CHrho_from_NM()
    !
    If(use_INM) Then
       Call calculate_C_from_NM(E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM)
    Else
       Crho(0)=Crho(0)+CHrho
    End If
    Call calculate_NM_properties()
    !
    Crho(0)=Crho(0)-CHrho
    !
    Call calculate_natural_units()
    !
    ! Print output
    !If(lpr) Then
    !   Call print_functional_parameters()
    !End If
    !
  End Subroutine set_functional_parameters
  !===============================================================================================
  !
  !===============================================================================================
  Subroutine print_functional_parameters()
    Use HFBTHO_utilities, Only: lout
    Implicit None
    !
    Write(lout,'(a)')   '  ---------------------------------------'
    Write(lout,'(a,a)') '           UNEDF Module Version:', &
                                         Trim(Version)
    Write(lout,'(a)')   '         M.Kortelainen & M.Stoitsov '
    Write(lout,'(a)')   '  ---------------------------------------'
    !
    Write(lout,'(a)')
    Write(lout,'(100(2x,a,a,f15.8))') Trim(FunctionalName),' functional'
    Write(lout,'(100(2x,a,f15.8))') '----------------------------------'
    Write(lout,'("  Crho(0)= ",g26.18,"; Crho(1)= ",g26.18)') Crho
    Write(lout,'("  CDrho(0)=",g26.18,"; CDrho(1)=",g26.18)') CDrho
    Write(lout,'("  Ctau(0)= ",g26.18,"; Ctau(1)= ",g26.18)') Ctau
    Write(lout,'("  CrDr(0)= ",g26.18,"; CrDr(1)= ",g26.18)') Crdr
    Write(lout,'("  CrdJ(0)= ",g26.18,"; CrdJ(1)= ",g26.18)') CrdJ
    Write(lout,'("  CJ(0)=   ",g26.18,"; CJ(1)=   ",g26.18)') CJ
    Write(lout,'("  CpV0(0)= ",g26.18,"; CpV0(1)= ",g26.18)') CpV0
    Write(lout,'("  CpV1(0)= ",g26.18,"; CpV1(1)= ",g26.18)') CpV1
    Write(lout,'("  sigma=   ",g26.18,"; hbzero=  ",g26.18)') sigma,hbzero
    Write(lout,'("  e^2 chrg=",g26.18,"; CExPar=  ",g26.18)') e2charg,CExPar
    Write(lout,'("  c.m. correction: ",L1,", chr. density in direct Coul: ",L1)') use_cm_cor,use_charge_density
    Write(lout,'("  use tensor terms: ",L1)') use_j2terms
    Write(lout,'("  use TMR pairing:  ",I1)') use_TMR_pairing
    !
    Write(lout,'(100(2x,a,f15.8))')
    Write(lout,'(100(2x,a,f15.8))') 'Coupling constants in natural units'
    Write(lout,'(100(2x,a,f15.8))') '-----------------------------------'
    Write(lout,'("  Crho_nu(0)= ",g26.18,"; Crho_nu(1)= ",g26.18)') nuCrho
    Write(lout,'("  CDrho_nu(0)=",g26.18,"; CDrho_nu(1)=",g26.18)') nuCDrho
    Write(lout,'("  Ctau_nu(0)= ",g26.18,"; Ctau_nu(1)= ",g26.18)') nuCtau
    Write(lout,'("  CrDr_nu(0)= ",g26.18,"; CrDr_nu(1)= ",g26.18)') nuCrdr
    Write(lout,'("  CrdJ_nu(0)= ",g26.18,"; CrdJ_nu(1)= ",g26.18)') nuCrdJ
    Write(lout,'("  CJ_nu(0)=   ",g26.18,"; CJ_nu(1)=   ",g26.18)') nuCJ
    Write(lout,'("  CpV0_nu(0)= ",g26.18,"; CpV0_nu(1)= ",g26.18)') nuCpV0
    Write(lout,'("  CpV1_nu(0)= ",g26.18,"; CpV1_nu(1)= ",g26.18)') nuCpV1
    Write(lout,'("  fpi_nu=     ",g26.18,"; Lambda_nu=  ",g26.18)') nufpi,nuLambda
    !
    If(dmeorder.Ge.0) Then
       Write(lout,'(100(2x,a,f15.8))')
       Write(lout,'(100(2x,a,f15.8))') 'DME parameters'
       Write(lout,'(100(2x,a,f15.8))') '----------------------------------'
       Write(lout,'("       gA=",f12.6," mpi [1/fm]=",f12.6," fpi [1/fm]=",f12.6)') gA,mpi,fpi
       Write(lout,'("  c1 [fm]=",f12.6,"    c3 [fm]=",f12.6,"    c4 [fm]=",f12.6)') c1,c3,c4
       Write(lout,'("       cd=",f12.6,"         ce=",f12.6," LamX[1/fm]=",f12.6)') cd,ce,LambdaX
       Write(lout,'("  ->CHrho=",f12.6)') CHrho
       If(dmeorder.Ge.2) Write(lout,'("  use 3N terms: ",L1)') use_DME3N_terms
    End If
    !
    Write(lout,'(100(2x,a,f15.8))')
    Write(lout,'(100(2x,a,f15.8))') 'Nuclear matter properties'
    Write(lout,'(100(2x,a,f15.8))') '----------------------------------'
    Write(lout,'(100(2x,a9,f25.16))') 'E_NM=',E_NM,'K_NM=',K_NM
    Write(lout,'(100(2x,a9,f25.16))') 'P_NM=',P_NM,'RHO_NM=',RHO_NM
    Write(lout,'(100(2x,a9,f25.16))') 'ASS_NM=',ASS_NM,'LASS_NM=',LASS_NM
    Write(lout,'(100(2x,a9,f25.16))') 'SMASS_NM=',SMASS_NM,'VMASS_NM=',VMASS_NM
    !
    Call t_from_C()
    !
    Write(lout,'(100(2x,a,f15.8))')
    Write(lout,'(100(2x,a,f15.8))') 'Associated (t,x)-coupling constants'
    Write(lout,'(100(2x,a,f15.8))') '-----------------------------------'
    Write(lout,'("  t0=    ",g26.18,"; x0=     ",g26.18)') t0,x0
    Write(lout,'("  t1=    ",g26.18,"; x1=     ",g26.18)') t1,x1
    Write(lout,'("  t2=    ",g26.18,"; x2=     ",g26.18)') t2,x2
    Write(lout,'("  t3=    ",g26.18,"; x3=     ",g26.18)') t3,x3
    Write(lout,'("  b4=    ",g26.18,"; b4p=    ",g26.18)') b4,b4p
    Write(lout,'("  te=    ",g26.18,"; to=     ",g26.18)') te,to
    Write(lout,'("  sigma= ",g26.18,"; hbzero= ",g26.18)') sigma,hbzero
    !
    If(Print_Namelist) Then
       Write(lout,'(100(2x,a,f15.8))')
       SELECTED_FUNCTIONAL: Select Case (Trim(FunctionalName))
       Case ("UNEDF","SKYRME")
             Write(lout,'(100(2x,a,f15.8))') 'NAMELIST CONTENT (cannot be modified for functional names UNEDF,SKYRME)'
             Write(lout,'(100(2x,a,f15.8))') '-----------------------------------------------------------------------'
       Case ("FITS")
             Write(lout,'(100(2x,a,f15.8))') 'NAMELIST CONTENT (Advanced usage: modify all but not C-, NM-, and more...)'
             Write(lout,'(100(2x,a,f15.8))') '-----------------------------------------------------------------------'
       Case default
             Write(lout,'(100(2x,a,f15.8))') 'NAMELIST CONTENT (copy/past to UNEDF_NAMELIST.DAT and modify)'
             Write(lout,'(100(2x,a,f15.8))') '-------------------------------------------------------------'
       End Select SELECTED_FUNCTIONAL
       Write(lout,'(100(a,f15.8))')    ' !NB: FUNCTIONALNAME should be always in quotations'
       Write(lout,UNEDF_NAMELIST)
    End If
  End Subroutine print_functional_parameters
  !===============================================================================================
  !
  !===============================================================================================
  Subroutine calculate_natural_units
    !--------------------------------------------------------------------------------
    ! Calculates coupling constants in natural units
    !--------------------------------------------------------------------------------
    Implicit None
    nuCrho = Crho*(nufpi**2)/(mevfm**3)
    nuCdrho = Cdrho*(nufpi**2)*((nuLambda*nufpi*nufpi)**sigma)/(mevfm**(3.0_pr*(1.0_pr+sigma)))
    nuCtau = Ctau*((nufpi*nuLambda)**2)/(mevfm**5)
    nuCrDr = CrDr*((nufpi*nuLambda)**2)/(mevfm**5)
    nuCrdJ = CrdJ*((nufpi*nuLambda)**2)/(mevfm**5)
    nuCJ = CJ*((nufpi*nuLambda)**2)/(mevfm**5)
    nuCpV0 = CpV0*(nufpi**2)/(mevfm**3)
    nuCpV1 = CpV1*(nufpi**4)*nuLambda/(mevfm**6)
  End Subroutine calculate_natural_units
  !===============================================================================================
  !
  !===============================================================================================
  Subroutine calculate_C_from_NM(E,K,SMASS,RHO,ASS,LASS,VMASS,sigma_NM)
    !--------------------------------------------------------------------------------
    ! Calculates volume C-constants (and sigma) form NM properties
    ! Interface usage:
    !  hbzero,CK,kfconst,mpi,sigma
    !  aRhoRho,bRhoRho...
    !  hRho0Rho0,dhRho0Rho0...
    !  Crho(0),Crho(1),Cdrho(0),Cdrho(1),Ctau(0),Ctau(0)
    !  subroutine calculate_U_parameters
    !
    !  input: E,K,SMASS,RHO,ASS,LASS,VMASS,sigma_NM(optional)
    ! output: Crho(0),Crho(1),Cdrho(0),Cdrho(1),Ctau(0),Ctau(0),sigma(optional)
    !
    ! Options:
    !  When sigma_NM exists then 'sigma'=sigma_NM
    !  When sigma_NM does not exist then 'sigma' is defined from NM
    !--------------------------------------------------------------------------------
    Implicit None
    Real(pr), Intent(in) :: E,K,SMASS,RHO,ASS,LASS,VMASS
    Real(pr), Intent(in), Optional :: sigma_NM
    Real(pr) :: aRho0Rho0,daRho0Rho0,ddaRho0Rho0,aRho1Rho1,daRho1Rho1,ddaRho1Rho1
    Real(pr) :: aRho0Tau0,daRho0Tau0,ddaRho0Tau0,aRho1Tau1,daRho1Tau1,ddaRho1Tau1
    Real(pr) :: u,tauc,rho2
    Real(pr),Parameter :: c13=1.0_pr/3.0_pr,c23=2.0_pr/3.0_pr
    !
    tauc=CK*RHO**c23; u=(kfconst/mpi)*RHO**c13; rho2=rho**2
    !
    Call calculate_U_parameters(RHO,RHO,tauc*RHO,tauc*RHO,0.0_pr,0.0_pr)
    !
    aRho0Rho0=0.50_pr*(aRhoRho+bRhoRho)*mevfm
    aRho1Rho1=0.50_pr*(aRhoRho-bRhoRho)*mevfm
    aRho0Tau0=0.50_pr*(aRhoTau+bRhoTau)*mevfm
    aRho1Tau1=0.50_pr*(aRhoTau-bRhoTau)*mevfm
    daRho0Rho0=0.50_pr*(daRhoRho+dbRhoRho)*mevfm
    daRho1Rho1=0.50_pr*(daRhoRho-dbRhoRho)*mevfm
    daRho0Tau0=0.50_pr*(daRhoTau+dbRhoTau)*mevfm
    daRho1Tau1=0.50_pr*(daRhoTau-dbRhoTau)*mevfm
    ddaRho0Rho0=0.50_pr*(ddaRhoRho+ddbRhoRho)*mevfm
    ddaRho1Rho1=0.50_pr*(ddaRhoRho-ddbRhoRho)*mevfm
    ddaRho0Tau0=0.50_pr*(ddaRhoTau+ddbRhoTau)*mevfm
    ddaRho1Tau1=0.50_pr*(ddaRhoTau-ddbRhoTau)*mevfm
    !
    ! set/calculate sigma
    If (Present(sigma_NM)) Then
        sigma=sigma_NM
    Else
        sigma=((1.0_pr/3.0_pr)*(-K+tauc*hbzero*(-3.0_pr+4.0_pr*SMASS)-9.0_pr*E+9.0_pr*RHO2*hRho0Rho0 &
             +21.0_pr*tauc*RHO2*hRho0Tau0+u*RHO*(daRho0Rho0+5.0_pr*tauc*daRho0Tau0 &
             +7.0_pr*RHO*dhRho0Rho0+11.0_pr*tauc*RHO*dhRho0Tau0+u*ddaRho0Rho0 &
             +u*tauc*ddaRho0Tau0+u*RHO*ddhRho0Rho0+u*tauc*RHO*ddhRho0Tau0))) &
             /(tauc*hbzero*(-3.0_pr+2.0_pr*SMASS)+3.0_pr*E+3.0_pr*RHO2*hRho0Rho0 &
             +3.0_pr*tauc*RHO2*hRho0Tau0+u*RHO*(daRho0Rho0+tauc*daRho0Tau0 &
             + RHO*dhRho0Rho0+tauc*RHO*dhRho0Tau0))
    End If
    !
    Crho(0)=(c13*(tauc*hbzero*(-3.0_pr+(2.0_pr-3.0_pr*sigma)*SMASS) &
        +3.0_pr*(1.0_pr+sigma)*E-3.0_pr*sigma*RHO*aRho0Rho0 &
        +3.0_pr*(1.0_pr-sigma)*RHO2*hRho0Rho0+3.0_pr*tauc*RHO2*hRho0Tau0 &
        +u*RHO*(daRho0Rho0+tauc*daRho0Tau0+RHO*dhRho0Rho0 &
        +tauc*RHO*dhRho0Tau0)))/(sigma*RHO)
    Cdrho(0)=(c13*RHO**(-1.0_pr-sigma)*(tauc*hbzero*(3.0_pr-2.0_pr*SMASS)&
        -3.0_pr*E-3.0_pr*RHO**2*hRho0Rho0-3.0_pr*tauc*RHO2*hRho0Tau0&
        -u*RHO*(daRho0Rho0+tauc*daRho0Tau0+RHO*dhRho0Rho0 &
        +tauc*RHO*dhRho0Tau0)))/sigma
    Ctau(0)=(hbzero*(SMASS-1.0_pr)-RHO*(aRho0Tau0+RHO*hRho0Tau0))/RHO
    !
    Crho(1)=(27.0_pr*ASS*(1.0_pr+sigma)-9.0_pr*LASS &
        +5.0_pr*tauc*hbzero*(5.0_pr-6.0_pr*VMASS+3.0_pr*sigma*(-4.0_pr+3.0_pr*VMASS)) &
        +20.0_pr*tauc*(2.0_pr-3.0_pr*sigma)*RHO*aRho0Tau0 &
        +RHO*(-27.0_pr*sigma*aRho1Rho1+5.0_pr*tauc*(11.0_pr-12.0_pr*sigma)*RHO*hRho0Tau0 &
        -27.0_pr*(-1.0_pr+sigma)*RHO*hRho1Rho1+9.0_pr*tauc*(5.0_pr-3.0_pr*sigma)*RHO*hRho1Tau0 &
        +45.0_pr*tauc*RHO*hRho1Tau1+40.0_pr*tauc*Ctau(0)-60.0_pr*tauc*sigma*Ctau(0) &
        +5.0_pr*u*tauc*daRho0Tau0+9.0_pr*u*daRho1Rho1+15.0_pr*u*tauc*daRho1Tau1 &
        +5.0_pr*u*tauc*RHO*dhRho0Tau0+9.0_pr*u*RHO*dhRho1Rho1+9.0_pr*u*tauc*RHO*dhRho1Tau0 &
        +15.0_pr*u*tauc*RHO*dhRho1Tau1))/(27.0_pr*sigma*RHO)
    Cdrho(1)=-(RHO**(-1.0_pr-sigma)*(27.0_pr*ASS-9.0_pr*LASS &
        +5.0_pr*tauc*hbzero*(5.0_pr-6.0_pr*VMASS)+40.0_pr*tauc*RHO*aRho0Tau0 &
        +55.0_pr*tauc*RHO2*hRho0Tau0+27.0_pr*RHO**2*hRho1Rho1+45.0_pr*tauc*RHO2*hRho1Tau0 &
        +45.0_pr*tauc*RHO2*hRho1Tau1+40.0_pr*tauc*RHO*Ctau(0) +5.0_pr*u*tauc*RHO*daRho0Tau0 &
        +9.0_pr*u*RHO*daRho1Rho1+15.0_pr*u*tauc*RHO*daRho1Tau1 &
        +5.0_pr*u*tauc*RHO2*dhRho0Tau0+9.0_pr*u*RHO2*dhRho1Rho1 &
        +9.0_pr*u*tauc*RHO2*dhRho1Tau0 +15.0_pr*u*tauc*RHO2*dhRho1Tau1))/(27.0_pr*sigma)
    Ctau(1)=(hbzero-hbzero*VMASS+RHO*(aRho0Tau0-aRho1Tau1+RHO*hRho0Tau0-RHO*hRho1Tau1+Ctau(0)))/RHO
    !
  End Subroutine calculate_C_from_NM
  !===============================================================================================
  !
  !===============================================================================================
  Subroutine calculate_NM_properties()
    !--------------------------------------------------------------------------------
    ! Calculates INM properties
    ! Interface usage:
    !  hbzero,CK,kfconst,mpi,sigma
    !  aRhoRho,bRhoRho...
    !  hRho0Rho0,dhRho0Rho0...
    !  Crho(0),Crho(1),Cdrho(0),Cdrho(1),Ctau(0),Ctau(0)
    !  E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM,sigma,P_NM,KA_NM
    !  function find_NM_RHOC()
    ! input:  Crho(0),Crho(1),Cdrho(0),Cdrho(1),Ctau(0),Ctau(0),sigma
    ! output: E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM,sigma,P_NM,KA_NM
    ! Using:
    !  RHO_NM=find_NM_RHOC()
    !--------------------------------------------------------------------------------
    Implicit None
    Real(pr) :: aRho0Rho0,daRho0Rho0,ddaRho0Rho0,aRho1Rho1,daRho1Rho1,ddaRho1Rho1
    Real(pr) :: aRho0Tau0,daRho0Tau0,ddaRho0Tau0,aRho1Tau1,daRho1Tau1,ddaRho1Tau1
    Real(pr) :: u,tauc,rho_NM2
    Real(pr), Parameter :: c13=1.0_pr/3.0_pr,c23=2.0_pr/3.0_pr
    !
    RHO_NM=find_NM_RHOC()
    !
    aRho0Rho0=0.50_pr*(aRhoRho+bRhoRho)*mevfm
    aRho1Rho1=0.50_pr*(aRhoRho-bRhoRho)*mevfm
    aRho0Tau0=0.50_pr*(aRhoTau+bRhoTau)*mevfm
    aRho1Tau1=0.50_pr*(aRhoTau-bRhoTau)*mevfm
    daRho0Rho0=0.50_pr*(daRhoRho+dbRhoRho)*mevfm
    daRho1Rho1=0.50_pr*(daRhoRho-dbRhoRho)*mevfm
    daRho0Tau0=0.50_pr*(daRhoTau+dbRhoTau)*mevfm
    daRho1Tau1=0.50_pr*(daRhoTau-dbRhoTau)*mevfm
    ddaRho0Rho0=0.50_pr*(ddaRhoRho+ddbRhoRho)*mevfm
    ddaRho1Rho1=0.50_pr*(ddaRhoRho-ddbRhoRho)*mevfm
    ddaRho0Tau0=0.50_pr*(ddaRhoTau+ddbRhoTau)*mevfm
    ddaRho1Tau1=0.50_pr*(ddaRhoTau-ddbRhoTau)*mevfm
    tauc=CK*RHO_NM**c23; u=(kfconst/mpi)*RHO_NM**c13; rho_NM2=rho_NM**2
    !
    ! Symmetric Nuclear Matter
    E_NM=tauc*hbzero+RHO_NM*(aRho0Rho0+RHO_NM*hRho0Rho0+Crho(0)+RHO_NM**sigma*Cdrho(0)) &
      +tauc*RHO_NM*(aRho0Tau0+RHO_NM*hRho0Tau0+Ctau(0))
    P_NM=c13*RHO_NM**2*((2.0_pr*tauc*hbzero)/RHO_NM+3.0_pr*aRho0Rho0+5.0_pr*tauc*aRho0Tau0 &
      +6.0_pr*RHO_NM*hRho0Rho0+8.0_pr*tauc*RHO_NM*hRho0Tau0+3.0_pr*Crho(0) &
      +3.0_pr*(1+sigma)*RHO_NM**sigma*Cdrho(0)+5.0_pr*tauc*Ctau(0)+u*daRho0Rho0 &
      +u*tauc*daRho0Tau0+u*RHO_NM*dhRho0Rho0+u*tauc*RHO_NM*dhRho0Tau0)
    SMASS_NM=1.0_pr+(RHO_NM*(aRho0Tau0+RHO_NM*hRho0Tau0+Ctau(0)))/hbzero
    K_NM=9.0_pr*sigma*(1+sigma)*RHO_NM**(1+sigma)*Cdrho(0) &
      +(-2.0_pr*tauc*hbzero+10.0_pr*tauc*RHO_NM*aRho0Tau0+18.0_pr*RHO_NM2*hRho0Rho0 &
      +40.0_pr*tauc*RHO_NM**2*hRho0Tau0+4.0_pr*u*RHO_NM*daRho0Rho0 &
      +RHO_NM*(10.0_pr*tauc*Ctau(0)+u*(8.0_pr*tauc*daRho0Tau0+u*ddaRho0Rho0 &
      +(10.0_pr*RHO_NM*dhRho0Rho0+14.0_pr*tauc*RHO_NM*dhRho0Tau0 &
      +(u*tauc*ddaRho0Tau0+u*RHO_NM*ddhRho0Rho0+u*tauc*RHO_NM*ddhRho0Tau0)))))
    !
    ! Asymmetric Nuclear Matter
    ASS_NM=RHO_NM2*hRho1Rho1+RHO_NM*(aRho1Rho1+Crho(1)+RHO_NM**sigma*Cdrho(1)) &
       +(tauc*(5.0_pr*hbzero+RHO_NM*(5.0_pr*aRho0Tau0+15.0_pr*aRho1Tau1+5.0_pr*RHO_NM*hRho0Tau0 &
       +9.0_pr*RHO_NM*hRho1Tau0+5.0_pr*(3.0_pr*RHO_NM*hRho1Tau1+Ctau(0)+3.0_pr*Ctau(1)))))/9.0_pr
    VMASS_NM=(hbzero+RHO_NM*(aRho0Tau0-aRho1Tau1+RHO_NM*hRho0Tau0-RHO_NM*hRho1Tau1+Ctau(0)-Ctau(1)))/hbzero
    LASS_NM=6.0_pr*RHO_NM2*hRho1Rho1+3.0_pr*RHO_NM*(aRho1Rho1+Crho(1)+(1.0_pr+sigma)*RHO_NM**sigma*Cdrho(1)) &
       +u*RHO_NM*daRho1Rho1 +u*RHO_NM2*dhRho1Rho1 &
       +(tauc*(10.0_pr*hbzero+8.0_pr*RHO_NM2*(5.0_pr*hRho0Tau0+9.0_pr*hRho1Tau0+15.0_pr*hRho1Tau1) &
       +25.0_pr*RHO_NM*(aRho0Tau0+3.0_pr*aRho1Tau1+Ctau(0)+3*Ctau(1)) &
       +5.0_pr*u*RHO_NM*(daRho0Tau0+3.0_pr*daRho1Tau1) &
       +u*RHO_NM2*(5.0_pr*dhRho0Tau0+9.0_pr*dhRho1Tau0+15.0_pr*dhRho1Tau1)))/9.0_pr
    KA_NM=18.0_pr*RHO_NM2*hRho1Rho1+9.0_pr*sigma*(1.0_pr+sigma)*RHO_NM**(1.0_pr+sigma)*Cdrho(1) &
       +4.0_pr*u*RHO_NM*daRho1Rho1 +10.0_pr*u*RHO_NM2*dhRho1Rho1 &
       + u**2*RHO_NM*ddaRho1Rho1+u**2*RHO_NM2*ddhRho1Rho1 &
       +(tauc*(-10.0_pr*hbzero+40.0_pr*RHO_NM2*(5.0_pr*hRho0Tau0+9.0_pr*hRho1Tau0+15.0_pr*hRho1Tau1) &
       +50.0_pr*RHO_NM*(aRho0Tau0+3.0_pr*aRho1Tau1+Ctau(0)+3*Ctau(1)) &
       +40.0_pr*u*RHO_NM*(daRho0Tau0+3.0_pr*daRho1Tau1) &
       +14.0_pr*u*RHO_NM2*(5.0_pr*dhRho0Tau0+9.0_pr*dhRho1Tau0 &
       +15.0_pr*dhRho1Tau1)+5.0_pr*u**2*RHO_NM*(ddaRho0Tau0 &
       +3.0_pr*ddaRho1Tau1)+u**2*RHO_NM2*(5.0_pr*ddhRho0Tau0+9*ddhRho1Tau0+15*ddhRho1Tau1)))/9.
     !
  End Subroutine calculate_NM_properties
  !===============================================================================================
  !
  !===============================================================================================
  Real(pr) Function find_NM_RHOC()
    !--------------------------------------------------------------------------------
    ! Search for the INM saturation density RHO_NM using the Secant Method
    !--------------------------------------------------------------------------------
    Implicit None
    !Integer(pr) intent(out) :: ierr
    Integer(ipr) :: iter
    Real(pr) :: aRho0Rho0,daRho0Rho0,ddaRho0Rho0,aRho1Rho1,daRho1Rho1,ddaRho1Rho1
    Real(pr) :: aRho0Tau0,daRho0Tau0,ddaRho0Tau0,aRho1Tau1,daRho1Tau1,ddaRho1Tau1
    Real(pr) :: kfconstmpi,u,tauc
    Real(pr) :: rhom0,rhom,rhom2,w,w0,step,energy
    Real(pr),Parameter :: c13=1.0_pr/3.0_pr,c23=2.0_pr/3.0_pr
    !
    kfconstmpi=kfconst/mpi; step=-0.0010_pr; iter=0
    ! initial point
    rhom=0.170_pr; tauc=CK*rhom**c23; u=kfconstmpi*rhom**c13; rhom2=rhom**2
    !
    Call calculate_U_parameters(rhom,rhom,tauc*rhom,tauc*rhom,0.0_pr,0.0_pr)
    !
    aRho0Rho0=0.50_pr*(aRhoRho+bRhoRho)*mevfm; daRho0Rho0=0.50_pr*(daRhoRho+dbRhoRho)*mevfm
    aRho0Tau0=0.50_pr*(aRhoTau+bRhoTau)*mevfm; daRho0Tau0=0.50_pr*(daRhoTau+dbRhoTau)*mevfm
    w0=c13*rhom2*((2.0_pr*tauc*hbzero)/rhom+3.0_pr*aRho0Rho0+5.0_pr*tauc*aRho0Tau0 &
      +6.0_pr*rhom*hRho0Rho0+8.0_pr*tauc*rhom*hRho0Tau0+3.0_pr*Crho(0) &
      +3.0_pr*(1.0_pr+sigma)*rhom**sigma*Cdrho(0)+5.0_pr*tauc*Ctau(0)+u*daRho0Rho0 &
      +u*tauc*daRho0Tau0+u*rhom*dhRho0Rho0+u*tauc*rhom*dhRho0Tau0)
    rhom0=rhom; rhom=rhom+step
    !
    ! secant method
    Do While(Abs(step).Ge.eps*100.0_pr)
       iter=iter+1
       tauc=CK*rhom**c23; u=kfconstmpi*rhom**c13; rhom2=rhom**2
       !
       Call calculate_U_parameters(rhom,rhom,tauc*rhom,tauc*rhom,0.0_pr,0.0_pr)
       !
       aRho0Rho0=0.50_pr*(aRhoRho+bRhoRho)*mevfm; daRho0Rho0=0.50_pr*(daRhoRho+dbRhoRho)*mevfm
       aRho0Tau0=0.50_pr*(aRhoTau+bRhoTau)*mevfm; daRho0Tau0=0.50_pr*(daRhoTau+dbRhoTau)*mevfm
       w=c13*rhom2*((2.0_pr*tauc*hbzero)/rhom+3.0_pr*aRho0Rho0+5.0_pr*tauc*aRho0Tau0 &
         +6.0_pr*rhom*hRho0Rho0+8.0_pr*tauc*rhom*hRho0Tau0+3.0_pr*Crho(0) &
         +3.0_pr*(1.0_pr+sigma)*rhom**sigma*Cdrho(0)+5.0_pr*tauc*Ctau(0)+u*daRho0Rho0 &
         +u*tauc*daRho0Tau0+u*rhom*dhRho0Rho0+u*tauc*rhom*dhRho0Tau0)
       step=-w*(rhom-rhom0)/(w-w0)
       rhom0=rhom; w0=w; rhom=rhom+step
       If(iter.Gt.100) Stop 'STOP(In find_NM_RHOC)'
       !energy=tauc*hbzero+rhom*(aRho0Rho0+rhom*hRho0Rho0+Crho(0)+rhom**sigma*Cdrho(0)) &
       ! +tauc*rhom*(aRho0Tau0+rhom*hRho0Tau0+Ctau(0))
       !Write(6,'(a,15(1pg12.4))') ' rhom0,rhom,step,e,w=',rhom0,rhom,step,energy,w
    End Do
    find_NM_RHOC=rhom
  End Function find_NM_RHOC
  !===============================================================================================
  !
  !===============================================================================================
  Subroutine C_from_t()
    !--------------------------------------------------------------------------------
    ! C- from (t,x)-
    !--------------------------------------------------------------------------------
    Implicit None
    Crho(0)  =   3.0_pr/8.0_pr  * t0
    Cdrho(0) =  (1.0_pr/16.0_pr)* t3
    Crho(1)  = -(1.0_pr/4.0_pr) * t0*(0.50_pr+x0)
    Cdrho(1) = -(1.0_pr/24.0_pr)* t3*(0.50_pr+x3)
    Ctau(0)  =  (3.0_pr/16.0_pr)* t1+(1.0_pr/4.0_pr)*t2*(5.0_pr/4.0_pr+x2)
    Ctau(1)  = -(1.0_pr/8.0_pr) * t1*(0.5+x1)+(1.0_pr/8.0_pr)*t2*(0.50_pr+x2)
    CrDr(0)  =  (1.0_pr/16.0_pr)* t2*(5.0_pr/4.0_pr+x2)-(9.0_pr/64.0_pr)*t1
    CrDr(1)  =  (3.0_pr/32.0_pr)* t1*(0.5+x1)+(1.0_pr/32.0_pr)*t2*(0.50_pr+x2)
    CJ(0)    = -(1.0_pr/16.0_pr)*(t1*(2.0_pr*x1-1.0_pr)+t2*(2.0_pr*x2+1)-5*te-15*to)
    CJ(1)    = -(1.0_pr/16.0_pr)*(t2 -t1 + 5.0_pr*te -5.0_pr*to )
    CrdJ(0)  = -b4-(0.50_pr)*b4p
    CrdJ(1)  = -0.50_pr*b4p
  End Subroutine C_from_t
  !===============================================================================================
  !
  !===============================================================================================
  Subroutine t_from_C()
    !--------------------------------------------------------------------------------
    ! (t,x)- from C-
    !--------------------------------------------------------------------------------
    Implicit None
    t0  =  (8.0_pr/3)*Crho(0)
    t1  =  4.0_pr/3.0_pr*(Ctau(0)-4.0_pr*CrDr(0))
    t2  =  4.0_pr/3.0_pr*(3.0_pr*Ctau(0)-6.0_pr*Ctau(1)+4.0_pr*CrDr(0)-8.0_pr*CrDr(1))
    t3  =  16.0_pr*Cdrho(0)
    x0  = -0.50_pr*(3.0_pr*Crho(1)/Crho(0)+1.0_pr)
    x1  =  2.0_pr*(-Ctau(0)-3.0_pr*Ctau(1)+4.0_pr*CrDr(0)+12.0_pr*CrDr(1))/t1/3.0_pr
    x2  = -2.0_pr*(3.0_pr*Ctau(0)-15.0_pr*Ctau(1)+4.0_pr*CrDr(0)-20.0_pr*CrDr(1))/t2/3.0_pr
    x3  = -0.50_pr*(3.0_pr*Cdrho(1)/Cdrho(0)+1.0_pr)
    b4  =  CrdJ(1)-CrdJ(0)
    b4p = -2.0_pr*CrdJ(1)
    te  = (4.0_pr/15.0_pr)*(3.0_pr*CJ(0)-9.0_pr*CJ(1)-4.0_pr*CrDr(0)+12.0_pr*CrDr(1)-2.0_pr*Ctau(0)+6.0_pr*Ctau(1))
    to  = (4.0_pr/15.0_pr)*(3.0_pr*CJ(0)+3.0_pr*CJ(1)+4.0_pr*CrDr(0)+4.0_pr*CrDr(1))
    ! Set public variables
    t0_pub = t0; t1_pub  = t1;  t2_pub = t2; t3_pub = t3
    x0_pub = x0; x1_pub  = x1;  x2_pub = x2; x3_pub = x3
    b4_pub = b4; b4p_pub = b4p; te_pub = te; to_pub = to
  End Subroutine t_from_C
  !===============================================================================================
  !
  !===============================================================================================
  Subroutine CHrho_from_NM()
    !--------------------------------------------------------------------------------
    ! CHrho from NM, E_NM(Hartree)=CHrho*RHO_NM
    !--------------------------------------------------------------------------------
    Implicit None
    Real(pr) :: z3=1.50_pr
    !
    !!CHrho= &
    !!+h0mpi6c3NM*(A3_1/b3_1**z3+A3_2/b3_2**z3+A3_3/b3_3**z3+A3_4/b3_4**z3+A3_5/b3_5**z3) &
    !!+h0mpi6c1NM*(A1_1/b1_1**z3+A1_2/b1_2**z3+A1_3/b1_3**z3+A1_4/b1_4**z3+A1_5/b1_5**z3)
    CHrho = 0.0_pr
    !
  End Subroutine CHrho_from_NM
  !===============================================================================================
  !
  !===============================================================================================
  Elemental Function HartreeV00(u)
    !--------------------------------------------------------------------------------
    ! HartreeV(u), E(Hartree)=(1/2)*Int[rho_0(r)*V(|r-r'|)*rho_0(r')]
    !--------------------------------------------------------------------------------
    Implicit None
    Real(pr), Intent(in) :: u
    Real(pr)             :: x2,HartreeV00
    !
    !!x2=(u*mpi)**2
    !
    !!HartreeV=h0mpi6c1*(Exp(-x2*b1_1)*A1_1+Exp(-x2*b1_2)*A1_2+Exp(-x2*b1_3)*A1_3+Exp(-x2*b1_4)*A1_4+Exp(-x2*b1_5)*A1_5)+&
    !!h0mpi6c3*(Exp(-x2*b3_1)*A3_1+Exp(-x2*b3_2)*A3_2+Exp(-x2*b3_3)*A3_3+Exp(-x2*b3_4)*A3_4+Exp(-x2*b3_5)*A3_5)
    !
    HartreeV00=0.0_pr
    !
  End Function HartreeV00
  !
  Elemental Function HartreeV01(u)
    !--------------------------------------------------------------------------------
    ! HartreeV(u), E(Hartree)=(1/2)*Int[rho_0(r)*V(|r-r'|)*rho_1(r')]
    !--------------------------------------------------------------------------------
    Implicit None
    Real(pr), Intent(in) :: u
    Real(pr)             :: x2,HartreeV01
    !
    HartreeV01=0.0_pr
    !
  End Function HartreeV01
  !
  Elemental Function HartreeV11(u)
    !--------------------------------------------------------------------------------
    ! HartreeV(u), E(Hartree)=(1/2)*Int[rho_1(r)*V(|r-r'|)*rho_1(r')]
    !--------------------------------------------------------------------------------
    Implicit None
    Real(pr), Intent(in) :: u
    Real(pr)             :: x2,HartreeV11
    !
    HartreeV11=0.0_pr
    !
  End Function HartreeV11
  !===============================================================================================
  !
  !===============================================================================================
  Elemental Function ThetaFunction2(u)
    !--------------------------------------------------------------------------------
    ! ThetaFunction2(u)=0 or 1  when x2<2  or x2>2
    !--------------------------------------------------------------------------------
    Implicit None
    Real(pr), Intent(in) :: u
    Real(pr)             :: x2,ThetaFunction2
    !
    x2=(u*mpi)
    !
    ThetaFunction2=0.0_pr
    If(x2.Gt.2.0_pr) ThetaFunction2=1.0_pr
    !
  End Function ThetaFunction2
  !===============================================================================================
  !
  !===============================================================================================
  Subroutine Make_Parameter_Free_Useful_Combinations()
    !--------------------------------------------------------------------------------
    ! Make Useful combinations
    !--------------------------------------------------------------------------------
    Implicit None
    !
    If(dmeorder.Ge.0) Then
        !
       mpi2=mpi**2
       gA2=gA**2; gA4=gA2**2; gA6=gA2**3;
       fpi2=fpi**2; fpi4=fpi2**2;
       CHartree =mevfm*(3.0_pr*gA2)/(32.0_pr*fpi4*Pi**2)
       h0mpi6=197.30_pr*(mpi**6)*(3.0_pr*gA*gA)/(32.0_pr*fpi**4*Pi**2)
       h0mpi6c1=h0mpi6*c1;         h0mpi6c3=h0mpi6*c3
       !
       h0mpi6NM=197.30_pr*(3.0_pr*(mpi**3)*gA2)/(64.0_pr*fpi**4*Sqrt(Pi))
       h0mpi6c1NM=h0mpi6NM*c1;     h0mpi6c3NM=h0mpi6NM*c3
       !
       A3_1=42.7132145164590_pr;   A3_2=0.670441422115440_pr; A3_3=0.0525713896514650_pr;
       A3_4=0.0012545731701320_pr; A3_5=5.81008627207380_pr
       b3_1=3.0809379008590_pr;    b3_2=0.905186811964580_pr; b3_3=0.474514509597610_pr;
       b3_4=0.228138177966090_pr;  b3_5=1.66931540698090_pr;
       !
       A1_1=2.5000830618386_pr;    A1_2=0.619542286897850_pr; A1_3=0.169682589033730_pr;
       A1_4=0.0276112113725470_pr; A1_5=0.00108164458809540_pr
       b1_1=1.75854210706510_pr;   b1_2=0.88882524524657_pr;  b1_3=0.46377235143756_pr;
       b1_4=0.247665887704790_pr;  b1_5=0.132222413002680_pr
       !
    End If
    !
  End Subroutine Make_Parameter_Free_Useful_Combinations
  !===============================================================================================
  !
  !===============================================================================================
  Elemental Function Vexternal(t,x,y,z)
    !
    Implicit None
    Integer(ipr), Intent(in) :: t  !! isospin index: 0=isoscalar, 1=isovector
    Real(pr), Intent(in) :: x,y,z  !! position in cartesian basis
    Real(pr) :: Vexternal
    !
    Vexternal = 0.0_pr
    !
  End Function Vexternal
  !
End Module UNEDF
!==================================================================================================================================
!#END UNEDF MODULE
!==================================================================================================================================
!#START HFBTHO MODULE
!==================================================================================================================================
Module HFBTHO
  Use HFBTHO_VERSION
  Use UNEDF
  !
  ! Input for HFBiterations
  Integer(ipr)  :: n00_INI,iLST_INI,inin_INI,icou_INI
  Integer(ipr)  :: npr_INI(3),kindhfb_INI
  Integer(ipr)  :: keyblo1_INI,keyblo2_INI,IDEBUG_INI
  Integer(ipr)  :: ngh_INI,ngl_INI,nleg_INI,nstate_INI
  Real(pr)      :: b0_INI,bz_INI,bp_INI,q_INI
  Character(30) :: skyrme_INI
  Real(pr)      :: pwi_INI,V0n_INI,V0p_INI,cpv1_INI,epsi_INI
  Logical       :: basis_HFODD_INI,Add_Pairing_INI,Print_HFBTHO_Namelist_INI,DO_FITT_INI
  ! Output for regression optimization
  Real(pr)  :: efit_0
  Real(pr), Dimension(0:1) :: efit_rhorho,efit_rhorhoD,efit_rhotau,efit_rhoDrho
  Real(pr), Dimension(0:1) :: efit_rhonablaJ,efit_JJ,efitV0,dfitV0,efV_0
  ! serial output (1:on/0:off)
  Integer(ipr)      :: IDEBUG
  Logical           :: DO_FITT
  ! For loop over used particle types. For normal nuclei min=1, max=2. For n droplets min=max=1.
  Integer(ipr) :: itmin,itmax,irestart
  ! Frequent Constants
  Real(pr) :: PI,ffdef3,ffdef4,ffdef5,ffdef6,ffdef7
  ! Single constants
  Real(pr) :: bet,beta0,q,bp,bpp,bz,hom,hb0,b0,etot,coex,cex,ty20,vin,rin,ain,  &
              qin,pwi,si,siold,epsi,xmix,xmix0,xmax,alst,clst,sklst,alphi,amas, &
              skass,varmas,v0ws,akv,hqc,amu,amn,r0,r00,r02,r04,decay,rmm3,amm3, &
              bmm3,cmm3,chargee2,EBASECUT,rho_c,cdef,b2_0,b4_0
  Integer(ipr) :: lin,lwin,lwou,lplo,lwel,lres,icstr,icou,ncut,iLST1,iLST, &
                  maxi,iiter,inin,nzm,nrm,icacou,iqrpa,icacoupj,icahartree,nlm,  &
                  nb,nt,n00,itass,kindhfb,iappend,iError_in_HO,iError_in_THO,    &
                  ierest,esu,nstate
  Integer(ipr), Parameter :: n00max=50
  ! Results
  Integer(ipr), Parameter :: ieresu=50,ieresl=20,ieresj=50,ieresbl=6
  Integer(ipr), Parameter :: ieres=ieresu+ieresl+ieresj+ieresbl
  Real(pr) :: eres(ieres)
  Character(13) :: ereslbl(2)
  Character(2) :: nucname
  Real(pr)    :: eresu(ieresu),eresl(ieresl),eresbl(ieresbl),eresj(ieresj)
  Character(13) :: hlabels(ieres+3)
  ! Common small arrays
  Real(pr) :: alast(2),ala(2),ala1(2),tz(2),ass(2),drhoi(2),del(2),vso(2),r0v(2)  &
       ,av(2),rso(2),aso(2),Sumnz(2),Dispersion(2),v2min(2),v2minv(2),rms(3),ept(3),q2(3)  &
       ,Dnfactor(3),varmasnz(2),pjmassnz(2)
  Integer(ipr) :: npr(3),inz(2),ldel(2),nk(2),itbl(2),kbl(2),tpar(2),ipbl(2),nbl(2),ibbl(2)  &
       ,klmax(2),inner(2),iasswrong(3),lcc
  ! Lipkin-Nogami
  Real(pr) :: ala2(2),etr(3),ssln(3,2),Geff(2)
  ! Blocking
  Real(pr)     :: pwiblo=2.0_pr, eqpmin(2)=0.0_pr
  Integer(ipr) :: bloall; Parameter(bloall=200)
  Integer(ipr), Dimension(0:bloall,2) :: bloblo,blo123=0,blok1k2=0
  Real(pr),     Dimension(0:bloall,2) :: bloqpdif
  Integer(ipr) :: iparenti(2),keyblo(3),nkblo_INI(2,5),nkblo(2,5)=0
  Integer(ipr) :: blocross(2),blomax(2),blo123d(2),blok1k2d(2),blocanon(2)
  ! manualBlocking
  Integer(ipr):: manualBlocking=0
  ! Logical and character variables
  Character(1) :: tq,tp(2),tl(0:20),tis(2)
  Character(30) :: skyrme
  Character(8) :: tit(2)
  Character(7) :: protn(2)
  Data  protn/'neutron','proton '/
  ! Allocatable arrays
  Character(13), Allocatable  :: tb(:)
  Character(25), Allocatable  :: txb(:)
  Real(pr), Allocatable, Target  :: rk(:,:),ak(:,:),hh0(:,:),de0(:,:)  &
       ,ddc(:,:,:),ddc1(:,:,:),qh(:,:),qh1(:,:),ql(:,:,:),ql1(:,:,:)  &
       ,ek(:,:),dk(:,:),vk(:,:),vk1(:,:),uk(:,:),hfb1(:,:),vkmax(:,:)
  Real(pr), Allocatable  :: fdsx(:),fdsy(:),fdsy1(:),fdsy2(:),fdsy3(:),fspb0(:)  &
       ,fspc0(:),fspd0(:),fspb1(:),fspc1(:),fspd1(:),fspb2(:),fspc2(:),fspd2(:),fspb3(:)  &
       ,fspc3(:),fspd3(:),fak(:),fi(:),sq(:),sqi(:),wf(:),wfi(:),rkass(:,:)
  Integer(ipr), Allocatable  :: id(:),ia(:),ikb(:),ipb(:),nz(:),nr(:),nl(:),ns(:),npar(:)  &
       ,ka(:,:),kd(:,:),numax(:,:),iv(:), lcanon(:,:)
  Real(pr), Allocatable  :: AN(:),ANK(:),PFIU(:),PFID(:)
  Real(pr), Allocatable  :: FIU(:),FID(:),FIUR(:),FIDR(:)
  Real(pr), Allocatable  :: FIUD2N(:),FIDD2N(:),FIUZ(:),FIDZ(:)
  ! constraints
  Integer(ipr), Parameter :: lambdaMax=8
  Integer(ipr) :: numberCons
  Integer(ipr), Allocatable :: multLambda(:)
  Real(pr), Dimension(0:8,1:3) :: qmoment
  Real(pr), Allocatable :: q_units(:),multLag(:),multRequested(:)
  Real(pr), Allocatable :: multMatElems(:)
  ! Temperature
  Logical :: switch_on_temperature
  Real(pr) :: temper
  Real(pr), Dimension(3) :: entropy
  Real(pr), Allocatable, Target :: fn_T(:),fp_T(:)
  ! optimization arrays
  Real(pr), Allocatable  :: QHLA_opt(:,:),FI1R_opt(:,:),FI1Z_opt(:,:),FI2D_opt(:,:),y_opt(:)
  ! Arrays depending on mesh points
  Integer(ipr)  :: ngh,ngl,nleg,nghl,nbx,ntx,nzx,nrx,nlx,ndx,ndx2,ndxs,nqx
  Integer(ipr)  :: nhfbqx,nb2x,nhfbx,nkx,nzrlx,iqqmax
  Real(pr), Allocatable :: xh(:),wh(:),xl(:),sxl(:),wl(:),xleg(:),wleg(:),vc(:,:)
  Real(pr), Allocatable :: vhbn(:),vn(:),vrn(:),vzn(:),vdn(:),vsn(:),dvn(:)
  Real(pr), Allocatable :: vhbp(:),vp(:),vrp(:),vzp(:),vdp(:),vsp(:),dvp(:)
  Real(pr), Allocatable :: vSZFIn(:),vSFIZn(:),vSRFIn(:),vSFIRn(:)
  Real(pr), Allocatable :: vSZFIp(:),vSFIZp(:),vSRFIp(:),vSFIRp(:)
  Real(pr), Allocatable, Target :: aka(:,:),ro(:,:),tau(:,:),dro(:,:),dj(:,:) &
       ,SZFI(:,:),SFIZ(:,:),SRFI(:,:),SFIR(:,:),NABLAR(:,:),NABLAZ(:,:)
  Real(pr), Allocatable  :: fl(:),fli(:),fh(:),fd(:),fp1(:),fp2(:),fp3(:),fp4(:),fp5(:),fp6(:)  &
       ,fs1(:),fs2(:),fs3(:),fs4(:),fs5(:),fs6(:),wdcor(:),wdcori(:),cou(:)
  Real(pr), Allocatable  :: vDHartree(:,:),vhart00(:,:),vhart01(:,:),vhart11(:,:)
  ! PAV Projection
  Integer(ipr) :: keypj,ilpj,ilpj2,ilnqx,ilnghl
  Real(pr)     :: rehfbcan,ehfb,retotpj,depnp,iproj,npr1pj,npr2pj
  Complex(pr)  :: onei=(0.0_pr,1.0_pr)
  Complex(pr), Allocatable, Target  :: phypj(:),sinphy(:),exp1iphy(:)  &
       ,exp1iphym(:),exp2iphy(:),exp2iphym(:),coupj(:,:),ropj(:,:,:)  &
       ,taupj(:,:,:),dropj(:,:,:),djpj(:,:,:),akapj(:,:,:),pjk(:,:)  &
       ,SZFIpj(:,:,:),SFIZpj(:,:,:),SRFIpj(:,:,:),SFIRpj(:,:,:),epj(:,:)  &
       ,ddepj(:,:,:),cpj(:,:,:),ypj(:,:,:),rpj(:,:,:)
  Real(pr)  :: polem(2),polep(2)
  ! CMC
  Integer(ipr)  :: ICMinput
  Real(pr)     :: ECMHFB(3),ECMPAV(3)
  ! CRC
  Integer(ipr)  :: ICRinput
  Real(pr) :: DEROT(3),SQUJ(3),CRAN(3),ERIGHFB(3)
  ! hfbdiagonal
  Real(pr), Allocatable :: erhfb(:),drhfb(:),erhfb1(:),drhfb1(:)
  Real(pr), Allocatable :: hfb(:,:),zhfb(:),evvk(:),hfbcan(:,:),evvkcan(:)
  ! Jason: def derived types
  Type :: ptr_to_2darray
     Real(pr),    Dimension(:,:),Allocatable :: arr
  End Type ptr_to_2darray
  Type :: ptr_to_array
     Real(pr),    Dimension(:),Allocatable :: arr
  End Type ptr_to_array
  Type :: ptr_to_iarray
     Integer(ipr), Dimension(:),Allocatable :: arr
  End Type ptr_to_iarray
  ! Jason: use derived types
  Type(ptr_to_2darray), Allocatable :: allhfb(:)
  Type(ptr_to_array),   Allocatable :: allevvk(:),allalwork(:)
  Type(ptr_to_iarray),  Allocatable :: alllwork(:),allISUPPZ(:)
  Integer(ipr),        Allocatable :: allibro(:),allIALWORK(:),allILWORK(:)
  Integer(ipr) :: oldnb
  Real(pr) :: cutoff_tol=1.d-6
  ! Broyden
  Character(1)   ::  bbroyden
  Integer(ipr)  :: nbroyden=7
  Real(pr)     :: alphamix=0.70_pr
  Integer(ipr)  :: nhhdim,nhhdim2,nhhdim3,nhhdim4,ialwork,ilwork
  Real(pr),    Allocatable, Target  :: brout(:),brin(:)
  Real(pr),     Allocatable :: alwork(:)
  Integer(ipr), Allocatable :: lwork(:)
  ! cm
  Real(pr)  :: facECM=1.0_pr
  ! new keys
  Logical :: set_pairing,basis_HFODD,Parity,Parity_INI
  Logical :: Print_Screen=.False.
  Logical :: Add_Pairing,Print_HFBTHO_Namelist
  Integer(ipr)  :: MAX_ITER_INI,keypj_INI,iproj_INI,npr1pj_INI,npr2pj_INI
  ! Eqp U,V
  Integer(ipr) :: nuv,nqp
  Real(pr), Allocatable, Target :: RVqpN(:),RVqpP(:),RUqpN(:),RUqpP(:),REqpN(:),REqpP(:)
  Integer(ipr), Allocatable, Target :: KpwiN(:),KpwiP(:),KqpN(:),KqpP(:)
  ! error indicator
  Integer(ipr)  :: ierror_flag=0
  Character(60) :: ierror_info(0:10)
  Character(30) :: welfile_INI
  ! Namelists
  Logical      :: add_initial_pairing, set_temperature, compatibility_HFODD, force_parity, &
                  user_pairing
  Integer(ipr) :: number_of_shells, proton_number, neutron_number, type_of_calculation, &
                  number_iterations, type_of_coulomb, restart_file, projection_is_on,   &
                  gauge_points, delta_Z, delta_N, switch_to_THO, number_Gauss,          &
                  number_Laguerre, number_Legendre, number_states, print_time
  Integer(ipr) :: proton_blocking(1:5), neutron_blocking(1:5), lambda_values(1:lambdaMax), &
                                                               lambda_active(1:lambdaMax)
  Real(pr)     :: oscillator_length, basis_deformation, accuracy, temperature, &
                  vpair_n, vpair_p, pairing_cutoff, pairing_feature
  Real(pr)     :: expectation_values(1:lambdaMax)
  Character(Len=30) :: functional
  !Namelist /HFBTHO_GENERAL/ number_of_shells,oscillator_length, basis_deformation, &
  !                          proton_number,neutron_number,type_of_calculation
  !Namelist /HFBTHO_ITERATIONS/ number_iterations, accuracy, restart_file
  !Namelist /HFBTHO_FUNCTIONAL/ functional, add_initial_pairing, type_of_coulomb
  !Namelist /HFBTHO_PAIRING/ user_pairing, vpair_n, vpair_p, pairing_cutoff, pairing_feature
  !Namelist /HFBTHO_CONSTRAINTS/ lambda_values, lambda_active, expectation_values
  !Namelist /HFBTHO_BLOCKING/ proton_blocking, neutron_blocking
  !Namelist /HFBTHO_PROJECTION/ switch_to_THO,projection_is_on,gauge_points,delta_Z,delta_N
  !Namelist /HFBTHO_TEMPERATURE/ set_temperature, temperature
  !Namelist /HFBTHO_DEBUG/ number_Gauss, number_Laguerre, number_Legendre, &
  !                        compatibility_HFODD, number_states, force_parity, print_time
  !
End Module HFBTHO
!==================================================================================================================================
!#END HFBTHO MODULE
!==================================================================================================================================
!#START HFBTHO_gauss
!==================================================================================================================================
Module HFBTHO_gauss

  Use HFBTHO_utilities
  Use HFBTHO

  Implicit None

Contains
  Subroutine gausspoints
    !---------------------------------------------------------------------
    ! The routine determines the points and weights for Gauss quadratures
    ! in the cases of Gauss-Legendre, -Laguerre and -Hermite formulas.
    !---------------------------------------------------------------------
    Implicit None
    Real(pr):: al,be,sparity
    Real(pr), Allocatable :: endpts(:),b(:),t(:),w(:)
    Integer(ipr) :: N,i,j,KINDI,kpts,nparity
    !
    al=0.0_pr; be=0.0_pr; kpts=0
    !
    !--------------------------------------------------------------------
    !------------------>> Gauss-Hermite (positive nodes) <<--------------
    !--------------------------------------------------------------------
    If(Parity) Then
       KINDI=4; N=2*ngh ! Parity conserved
       nparity=ngh; sparity=two
    Else
       KINDI=4; N=ngh   ! Parity not conserved
       nparity=0;   sparity=one
    End If
    Allocate(endpts(2)); Allocate(b(N),t(N),w(N))
    Call Gaussq(KINDI,N,al,be,kpts,endpts,b,t,w)
    If(ierror_flag.Ne.0) Return
    Do i=nparity+1,N
       j=i-nparity
       xh(j)=t(i)
       ! Build in the Gaussian weight function into the weights wh
       wh(j)=sparity*Exp(xh(j)*xh(j)+Log(w(i)))
    End Do
    Deallocate(endpts,b,t,w)
    !--------------------------------------------------------------------
    !---------------------------->> Gauss-Laguerre <<--------------------|
    !--------------------------------------------------------------------
    KINDI=6; N=ngl
    Allocate(endpts(2)); Allocate(b(N),t(N),w(N))
    Call Gaussq(KINDI,N,al,be,kpts,endpts,b,t,w)
    If(ierror_flag.Ne.0) Return
    Do j=1,ngl
       xl(j)=t(j)
       ! Build in the exponential weight function into the weights wl
       wl(j)=Exp(xl(j)+Log(w(j)))
       sxl(j)=Sqrt(xl(j))
    End Do
    Deallocate(endpts,b,t,w)
    !--------------------------------------------------------------------
    !----------------->> Gauss-Legendre (positive nodes) <<--------------
    !--------------------------------------------------------------------
    If(nleg.Gt.0) Then
       KINDI=1; N=2*nleg
       Allocate(endpts(2)); Allocate(b(N),t(N),w(N))
       Call Gaussq(KINDI,N,al,be,kpts,endpts,b,t,w)
    If(ierror_flag.Ne.0) Return
       Do j=1,nleg
          i=nleg+j
          xleg(j)=t(i); wleg(j)=w(i)
    End Do
       Deallocate(endpts,b,t,w)
    End If
    !
  End Subroutine gausspoints
  !=======================================================================
  !
  !=======================================================================
  Subroutine Gaussq(kindi,n,alpha,beta,kpts,endpts,b,t,w)
    Implicit None
    Integer(ipr) :: N,kindi
    Real(pr):: alpha,beta,MUZERO,GAM,T1
    Integer(ipr) :: j1,J2,kpts,ierr
    Real(pr):: b(n),t(n),w(n),endpts(2)
    !--------------------------------------------------------------------
    ! This set of routines computes the nodes t(j) and weights w(j) for
    ! Gaussian-type quadrature rules with pre-assigned nodes. These are
    ! used when one wishes to approximate
    !
    !          integral (from a to b)  f(x) w(x) dx
    !
    !                       n
    ! by                   sum w  f(t )
    !                      j=1  j    j
    !
    ! (note w(x) and w(j) have no connection with each other). Here w(x)
    ! is one of six possible non-negative weight functions (listed below),
    ! and f(x) is the function to be integrated. Gaussian quadrature is
    ! particularly useful on infinite intervals (with appropriate weight
    ! functions), since then other techniques often fail. Associated with
    ! each weight function w(x) is a set of orthogonal polynomials. The
    ! nodes t(j) are just the zeroes of the proper n-th degree polynomial.
    !
    ! inputs (all real numbers are in double precision)
    !    kindi ..: an integer between 1 and 6 giving the type of
    !              quadrature rule:
    !                1:  Legendre quadrature, w(x) = 1 on [-1, 1]
    !                2:  Chebyshev quadrature of the first kind
    !                    w(x) = 1/sqrt(1 - x*x) on [-1, +1]
    !                3:  Chebyshev quadrature of the second kind
    !                    w(x) = sqrt(1 - x*x) on p-1, 1]
    !                4:  Hermite quadrature, w(x) = exp(-x*x) on
    !                    ]-infinity, +infinity[
    !                5:  Jacobi quadrature, w(x) = (1-x)**alpha * (1+x)**
    !                    beta on [-1, 1], alpha, beta > -1.
    !                    note: kind=2 and 3 are a special case of this.
    !                6:  generalized Laguerre quadrature, w(x) = exp(-x)*
    !                x**alpha on [0, +infinity[, alpha > -1
    !    n .....: the number of points used for the quadrature rule
    !    alpha .: real parameter used only for Gauss-Jacobi and Gauss-
    !             Laguerre quadrature (otherwise use 0.d0).
    !    beta ..: real parameter used only for Gauss-Jacobi quadrature
    !             (otherwise use 0.d0)
    !    kpts ..: (integer) normally 0, unless the left or right end-
    !             point (or both) of the interval is required to be a
    !             node (this is called Gauss-Radau or Gauss-Lobatto
    !             quadrature). Then kpts is the number of fixed
    !             endpoints (1 or 2).
    !    endpts : real array of length 2. Contains the values of any
    !             fixed endpoints, if kpts = 1 or 2.
    !    b .....: real scratch array of length n
    !
    ! outputs (both double precision arrays of length n)
    !    t .....: the desired nodes.
    !    w .....: the desired weights w(j).
    !
    ! NOTE: Underflow may sometimes occur, but is harmless.
    !
    ! Adapted from GO library at www.netlib.org
    !--------------------------------------------------------------------
    !
    Call Class(kindi,n,alpha,beta,b,t,muzero)
    !
    If(KPTS.Eq.0) Then
       W=0.0_pr; W(1)=1._pr
       Call GBTQL2(N,T,B,W,IERR)
       W=MUZERO*W*W
       Return
    End If
    If(KPTS.Eq.2) Then
       GAM=GBSLVE(ENDPTS(1),N,T,B)
       T1=((ENDPTS(1)-ENDPTS(2))/(GBSLVE(ENDPTS(2),N,T,B)-GAM))
       B(N-1)=Sqrt(T1)
       T(N)=ENDPTS(1)+GAM*T1
       W=0.0_pr; W(1)=1._pr
       Call GBTQL2(N,T,B,W,IERR)
       W=MUZERO*W*W
       Return
    End If
    T(N)=GBSLVE(ENDPTS(1),N,T,B)*B(N-1)**2+ENDPTS(1)
    W=0.0_pr; W(1)=1._pr
    Call GBTQL2(N,T,B,W,IERR)
    W=MUZERO*W*W
  End Subroutine Gaussq
  !=======================================================================
  !
  !=======================================================================
  Real(pr) Function GBSLVE(SHIFT,N,A,B)
    Implicit None
    Integer(ipr) :: N,NM1,i
    Real(pr) :: ALPHA,SHIFT,A(N),B(N)
    ALPHA=A(1)-SHIFT
    NM1=N-1
    Do I=2,NM1
       ALPHA=A(I)-SHIFT-B(I-1)**2/ALPHA
    End Do
    GBSLVE=1.0_pr/ALPHA
  End Function GBSLVE
  !=======================================================================
  !
  !=======================================================================
  Subroutine Class(kindi,N,ALPHA,BETA,B,A,MUZERO)
    Implicit None
    Integer(ipr) :: N,kindi,i,NM1
    Real(pr) :: MUZERO,ALPHA,BETA,A(N),B(N)
    Real(pr) :: PI,ABI,DI20,AB,A2B2,FI
    Data PI / 3.1415926535897930_pr /
    !--------------------------------------------------------------------
    ! This procedure supplies the coefficients a(j), b(j) of the
    ! recurrence relation
    !
    !      b p (x) = (x - a ) p   (x) - b   p   (x)
    !       j j            j   j-1       j-1 j-2
    !
    ! for the various classical (normalized) orthogonal polynomials,
    ! and the zero-th moment
    !
    !      muzero = integral w(x) dx
    !
    ! of the given polynomial's weight function w(x). Since the
    ! polynomials are orthonormalized, the tridiagonal matrix is
    ! guaranteed to be symmetric.
    !
    ! The input parameter alpha is used only for Laguerre and Jacobi
    ! polynomials, and the parameter beta is used only for Jacobi
    ! polynomials. The Laguerre and Jacobi polynomials require the Gamma
    ! function.
    !
    ! Adapted from GO library at www.netlib.org
    !--------------------------------------------------------------------
    NM1=N-1
    Select Case (kindi)
    Case (1) ! Legendre polynomials
       MUZERO=2.0_pr
       Do I=1,NM1
          A(I)=0.0_pr
          ABI=Real(I,Kind=pr)
          B(I)=ABI/Sqrt(4.0_pr*ABI*ABI-1.0_pr)
       End Do
       A(N)=0.0_pr
    Case (2) ! Chebyshev polynomials of the first kind
       MUZERO=PI
       Do I=1,NM1
          A(I)=0.0_pr
          B(I)=0.50_pr
       End Do
       B(1)=Sqrt(0.50_pr)
       A(N)=0.0_pr
    Case (3) ! Chebyshev polynomials of the second kind
       MUZERO=PI/2.0_pr
       Do I=1,NM1
          A(I)=0.0_pr
          B(I)=0.50_pr
       End Do
       A(N)=0.0_pr
    Case (4) ! Hermite polynomials
       MUZERO=Sqrt(PI)
       Do I=1,NM1
          A(I)=0.0_pr
          DI20=I/2.0_pr
          B(I)=Sqrt(DI20)
       End Do
       A(N)=0.0_pr
    Case (5) ! Jacobi polynomials
       AB=ALPHA+BETA
       ABI=2.0_pr+AB
       MUZERO=2.0_pr**(AB+1.0_pr)*pr_gamma(ALPHA+1.0_pr)*pr_gamma(BETA+1.0_pr)/pr_gamma(ABI)
       A(1)=(BETA-ALPHA)/ABI
       B(1)=Sqrt(4.0_pr*(1.0_pr+ALPHA)*(1.0_pr+BETA)/((ABI+1.0_pr)*ABI*ABI))
       A2B2=BETA*BETA-ALPHA*ALPHA
       Do I=2,NM1
          ABI=2.0_pr*I+AB
          A(I)=A2B2/((ABI-2.0_pr)*ABI)
          FI=I
          B(I)=Sqrt(4.0_pr*FI*(FI+ALPHA)*(FI+BETA)*(FI+AB)/((ABI*ABI-1.0_pr)*ABI*ABI))
       End Do
       ABI=2.0_pr*N+AB
       A(N)=A2B2/((ABI-2.0_pr)*ABI)
    Case (6) ! Laguerre polynomials
       MUZERO=pr_gamma(ALPHA+1.0_pr)
       Do I=1,NM1
          FI=I
          A(I)=2.0_pr*FI-1.0_pr+ALPHA
          B(I)=Sqrt(FI*(FI+ALPHA))
       End Do
       A(N)=2.0_pr*N-1.0_pr+ALPHA
    Case default
    End Select
  End Subroutine Class
  !=======================================================================
  !
  !=======================================================================
  Subroutine GBTQL2(N,D,E,Z,IERR)
    Implicit None
    Integer(ipr) :: N,IERR
    Real(pr) :: D(N),E(N),Z(N)
    Integer(ipr) :: I,J,K,L,M,II,MML
    Real(pr) :: MACHEP,P,G,R,S,C,F,B
    !MACHEP=16.0_pr**(-14)
    MACHEP=epsilon(1.0_pr)
    IERR=0
    If(N.Eq.1) Return
    E(N)=0.0_pr
    Do L= 1,N
       J=0
       Do
          Do M=L,N
             If(M .Eq. N) Exit
             If(Abs(E(M)) .Le. MACHEP*(Abs(D(M))+Abs(D(M+1)))) Exit
             Continue
          End Do
          P=D(L)
          If(M .Eq. L) Exit
          If(J .Eq. 30) Then
             IERR=L
             Return
          End If
          J=J+1
          G=(D(L+1)-P) / (2.0_pr*E(L))
          R=Sqrt(G*G+1.0_pr)
          G=D(M) - P + E(L)/(G+Sign(R,G))
          S=1.0_pr
          C=1.0_pr
          P=0.0_pr
          MML=M-L
          Do II=1, MML
             I=M-II
             F=S*E(I)
             B=C*E(I)
             If(Abs(F).Ge.Abs(G)) Then
                C=G/F
                R=Sqrt(C*C+1.0_pr)
                E(I+1)=F*R
                S=1.0_pr/R
                C=C*S
             Else
                S=F/G
                R=Sqrt(S*S+1.0_pr)
                E(I+1)=G*R
                C=1.0_pr/R
                S=S*C
             End If
             G=D(I+1)-P
             R=(D(I)-G)*S + 2.0_pr*C*B
             P=S*R
             D(I+1)=G+P
             G=C*R - B
             F=Z(I+1)
             Z(I+1)=S*Z(I) + C*F
             Z(I)=C*Z(I) - S*F
          End Do
          D(L)=D(L)-P
          E(L)=G
          E(M)=0.0_pr
       End Do
    End Do
    Do II=2, N
       I=II-1
       K=I
       P=D(I)
       Do J=II,N
          If(D(J) .Ge. P) Cycle
          K=J
          P=D(J)
       End Do
       If(K .Eq. I) Cycle
       D(K)=D(I)
       D(I)=P
       P=Z(I)
       Z(I)=Z(K)
       Z(K)=P
    End Do
  End Subroutine GBTQL2
  !=======================================================================
  !
  !=======================================================================
  Real(pr) Function pr_gamma(x)
    !---------------------------------------------------------------------
    !  pr_gamma evaluates Gamma(X) for a real argument.
    !
    !  Discussion:
    !    This function was originally named DGAMMA. However, a number of
    !    Fortran compilers now include a library function of this name. To
    !    avoid conflicts, this function was renamed pr_gamma.
    !
    !    This routine calculates the GAMMA function for a real argument X.
    !    Computation is based on an algorithm outlined in reference 1.
    !    The program uses rational functions that approximate the GAMMA
    !    function to at least 20 significant decimal digits. Coefficients
    !    for the approximation over the interval (1,2) are unpublished.
    !    Those for the approximation for 12 <= X are from reference 2.
    !
    !  Licensing:
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !    18 January 2008
    !
    !  Author:
    !    Original FORTRAN77 version by William Cody, Laura Stoltz.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !    William Cody,
    !      An Overview of Software Development for Special Functions,
    !      in Numerical Analysis Dundee, 1975,
    !      edited by GA Watson,
    !      Lecture Notes in Mathematics 506,
    !      Springer, 1976.
    !    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
    !    Charles Mesztenyi, John Rice, Henry Thatcher, Christoph Witzgall,
    !      Computer Approximations,
    !      Wiley, 1968,
    !      LC: QA297.C64.
    !
    !  Parameters:
    !    Input, real ( kind = 8 ) X, the argument of the function.
    !    Output, real ( kind = 8 ) R8_GAMMA, the value of the function.
    !---------------------------------------------------------------------
    Implicit None
    !
    !  Coefficients for minimax approximation over (12, INF).
    Logical :: parity
    Integer(ipr) :: i,n
    Real(pr), Dimension(7) :: c = (/ -1.910444077728000000000D-03, &
      8.417138778129500000000000D-04,-5.952379913043012000000D-04, &
      7.936507935003502480000000D-04,-2.777777777777681622553D-03, &
      8.333333333333333331554247D-02, 5.708383526100000000000D-03 /)
    Real(pr) :: eps,fact,pi,res,sqrtpi,sum,x,xbig,xden,xinf,xminin,xnum,y,y1,ysq,z
    Real(pr) :: p(8),q(8)
    !
    ! Mathematical constants
    Data sqrtpi /0.9189385332046727417803297D+00/
    Data pi /3.1415926535897932384626434D+00/
    !
    ! Numerator and denominator coefficients for rational minimax
    ! approximation over (1,2).
    Data p / -1.71618513886549492533811D+00,  2.47656508055759199108314D+01, &
             -3.79804256470945635097577D+02,  6.29331155312818442661052D+02, &
              8.66966202790413211295064D+02, -3.14512729688483675254357D+04, &
             -3.61444134186911729807069D+04,  6.64561438202405440627855D+04 /

    Data q / -3.08402300119738975254353D+01,  3.15350626979604161529144D+02, &
             -1.01515636749021914166146D+03, -3.10777167157231109440444D+03, &
              2.25381184209801510330112D+04,  4.75584627752788110767815D+03, &
             -1.34659959864969306392456D+05, -1.15132259675553483497211D+05 /

    parity = .False.; fact = one; n = 0; y = x
    xbig = 171.624D+00
    xminin = Tiny(1.0_pr); eps = Epsilon(1.0_pr) ; xinf = Huge(1.0_pr)
    !
    ! Argument is negative.
    If(y<=zero) Then

       y = - x
       y1 = Aint ( y )
       res = y - y1

       If(res/=zero) Then
          If(y1/=Aint(y1*half)*two) Then
             parity = .True.
          End If
          fact = - pi / Sin(pi*res)
          y = y + one
       Else
          res = xinf
          pr_gamma = res
          Return
       End If

    End If
    !
    ! Argument is positive.
    if (y < eps) Then
       !
       ! Argument < EPS.
       If(xminin <= y) Then
          res = one / y
       Else
          res = xinf
          pr_gamma = res
          Return
       End If

    Else If(y<pp12) Then

       y1 = y

       ! 0.0 < argument < 1.0.
       If(y<one) Then
          z = y
          y = y + one
       Else
          ! 1.0 < argument < 12.0.
          ! Reduce argument if necessary.
          n = Int(y) - 1
          y = y - Real(n,Kind=pr)
          z = y - one
       End If
       !
       ! Evaluate approximation for 1.0 < argument < 2.0.
       xnum = zero
       xden = one
       Do i = 1, 8
          xnum = ( xnum + p(i) ) * z
          xden = xden * z + q(i)
       End Do

       res = xnum / xden + one
       !
       ! Adjust result for case  0.0 < argument < 1.0.
       If(y1<y) Then
          res = res / y1
       Else If(y<y1) Then
          ! Adjust result for case 2.0 < argument < 12.0.
          Do i = 1, n
             res = res * y
             y = y + one
          End Do
       End If

    Else
       !
       !  Evaluate for 12.0 <= argument.
       If(y<=xbig) Then
          ysq = y * y
          sum = c(7)
          Do i = 1, 6
             sum = sum / ysq + c(i)
          End Do
          sum = sum / y - y + sqrtpi
          sum = sum + ( y - half ) * Log(y)
          res = Exp(sum)
       Else
          res = huge ( res )
          pr_gamma = res
          Return
       End If

    End If
    !
    !  Final adjustments and Return.
    If(parity) Then
       res = - res
    End If

    If(fact/=one) Then
       res = fact / res
    End If

    pr_gamma = res

    Return
  End Function pr_gamma  ! Gamma function in double precision
End Module HFBTHO_gauss
!==================================================================================================================================
!#END HFBTHO_gauss
!==================================================================================================================================
!#START THO MODULE
!==================================================================================================================================
Module HFBTHO_THO
  Use HFBTHO_utilities
  Use linear_algebra
  Implicit None
  !
  Character(6), Private :: THO_version='1'
  !
Contains
  !
  Subroutine f01234(lpr)
    !--------------------------------------------------------------------------
    ! Calculates LST-function 'f' its derivatives 'f1,f2,f3,f4', the J
    ! Jacobian 'fd=(f^2 f1/r^2)^(1/2)' and its first and second derivatives
    ! 'fd1' and 'fd2' at the point 'r'. The LST function is defined as in
    ! routine 'thofun()'.
    !
    ! All integrations are performed in the dimensionless variables
    ! 'u=xh(ih)' and 'v=xl(il)' being the dimensionless Gauss mesh
    ! points in cylindrical coordinates. Thus 'r' is found as the inverse
    ! function of the input variable 'qq=Sqrt[xh(ih)^2+xl(il)]', i.e.,
    ! solving f(r)=qq. The original 'z' and 'rho' variables in 'fm' are
    ! defined as 'fh(ihli) = zz*bz' and 'fl(ihli) = rr*bp'.
    !
    ! All required quantities for evaluating w.f., its first derivative
    ! and the associated Laplacian are calculated
    !
    ! Definitions used:
    !  - xh(ih),xl(il) : Gauss mesh points
    !  - bp, bz        : oscillator lengths in fm
    !  - ilst1         : 0->ho, #0->tho
    !
    ! In this way:
    ! (1) \delta \rho divided by j^2 is:
    !        dro = (cz*fs1+cr*fs2+fs3)*qhab*qlab
    !            +  fs1*qha1b1*qlab
    !            +  fs2*qhab*qla1b1/(4.*v**2)
    !            +  fs4*(qha1b+qhab1)*(qla1b+qlab1)/(2.*v)
    !            +  fs5*(qha1b+qhab1)*qlab
    !            +  fs6*qhab*(qla1b+qlab1)/(2.*v)
    !     with
    !        cr = 1./4. - (nr_a+nr_b+m+1)/(2.*v) + (m/(2.*v))**2
    !        cz = u*u - (nz_a+nz_b+1)
    ! (2) the first (r,z) derivatives are:
    !        fi'_z = fp1*(qh*ql) + fp2*(qh1*ql) + fp3*(qh*ql1)/(2.*v)
    !              -> HO -> qh1*ql/bpz
    !        fi'_r = fp4*(qh*ql) + fp5*(qh1*ql) + fp6*(qh*ql1)/(2.*v)
    !              -> HO -> qh*ql1/Sqrt(v)/bp
    !--------------------------------------------------------------------------
     Use HFBTHO
     Implicit None
     Logical :: lpr
     Integer(ipr), Save :: i,il,ih,ihli,iw,key0=0,key1=1
    Real(pr), Save :: bri,bri2,bzi,bzi2,wv1,u,qq,f,f1,f2,f3,rhoi,fd1,fd2,         &
                      fd12,fdd,r,r2,r3,r4,r5,r6,rr,rr2,rr4,zz,zz2,drr1,drr2,drz1, &
                      drz2,drr12,drz12,g,g1,g2,g3,gg,gg1,gg2,g1g1,uz1,ur1,vz1,vr1,&
                      uz2,ur2,vz2,vr2,ur12,vr12,uz12,vz12
     !
     bri= one/bp; bri2= bri*bri; bzi= one/bz; bzi2= bzi*bzi
     Do il=1,ngl
        wv1   = xl(il)
        Do ih=1,ngh
           ihli = ih+(il-1)*ngh
           u    = xh(ih); qq   = Sqrt(u*u+wv1)
           If(ilst1.Eq.0) Then
              ! ho-case
              r = qq; f = r; f1 = one; f2 = zero; f3 = zero
           Else
              ! tho-case: initial run
              If(ih*il.Eq.1.And.ilst.Lt.0) &
                 Call thofun(key0,g,f,f1,f2,f3,g1,.True.,.False.)
              If(iasswrong(3).ne.0) then
                 ! reinforce ho results
                 r = qq; f = r; f1 = one; f2 = zero; f3 = zero
              else
                 ! tho-case: f(r)=qq,f'(r),f''(r),f'''(r), r=Invers_f(r)
                 Call thofun(key1,qq,f,f1,f2,f3,r,.False.,.False.)
              End If
           End If
           ! Jacobian calculations
           r2= r*r; r3= r2*r; r4= r2*r2 ;r5= r3*r2 ;r6= r3*r3
           ! fdd=(f(r)^2 f'(r)/r^2)^(1/2),fd1=fdd'/fdd, fd2=fdd''/fdd
           fd1  = f1/f - one/r + half*f2/f1
           fdd  = f*Sqrt(f1)/r
           fd2  = two*(f2/f-fd1/r) - (f2/f1)**2/four + half*f3/f1
           fd12 = fd1**2
           ! g=(f/r)-derivatives
           g    = f/r; g1   =-(f - f1*r)/r2
           g2   = (two*(f - f1*r) + f2*r2)/r3
           g3   = (six*(f1*r - f) - three*f2*r2 + f3*r3)/r4
           ! g4   = (24.0d0*(f-f1*r+half*f2*r2)-four*f3*r3+f4*r4)/r5
           gg   = g*g; gg1= g*g1; gg2= g*g2; g1g1= g1*g1
           ! (rr,zz)-definitions
           rr   = Sqrt(wv1)/g; rhoi = bri/rr; zz   = u/g
           rr2  = rr*rr; rr4= rr2*rr2 ;zz2= zz*zz
           ! (r,z)-derivatives
           drr1  = bri*rr/r;           drz1  = bzi*zz/r
           drr2  = (bri2 - drr1**2)/r; drz2  = (bzi2 - drz1**2)/r
           drr12 = drr1**2;            drz12 = drz1**2
           ! (u,v)-derivatives
           uz1  = bzi*(g+g1*zz2/r);    ur1  = bri*g1*rr*zz/r
           vz1  = bzi*two*gg1*rr2*zz/r
           vr1  = bri*two*rr*(gg+gg1*rr2/r)
           uz2  = bzi2*zz*(three*g1*r2-g1*zz2+g2*r*zz2)/r3
           ur2  = bri2*(g1*r2-g1*rr2+g2*r*rr2)*zz/r3
           vz2  = bzi2*two*rr2*(gg1*r2-gg1*zz2+g1g1*r*zz2+gg2*r*zz2)/r3
           vr2  = gg*r3+five*gg1*r2*rr2-gg1*rr4+g1g1*r*rr4+gg2*r*rr4
           vr2  = vr2*bri2*two/r3; ur12=ur1**2; vr12=vr1**2; uz12=uz1**2
           vz12 = vz1**2
           ! storage
           fh(ihli)= zz*bz; fl(ihli)= rr*bp; fli(ihli)=one/fl(ihli); fd(ihli)= fdd*fdd
           ! for first derivatives
           fp1(ihli)= fd1*drz1; fp2(ihli)= uz1; fp3(ihli)= vz1
           fp4(ihli)= fd1*drr1; fp5(ihli)= ur1; fp6(ihli)= vr1
           ! for the Laplacian
           fs1(ihli) = two*(ur12 + uz12); fs2(ihli) = two*(vr12 + vz12)
           fs3(ihli) = two*(fd1*(drr1*rhoi+drr2+drz2)+ &
                (fd12+fd2)*(drr12+drz12))
           fs4(ihli) = two*(ur1*vr1 + uz1*vz1)
           fs5(ihli) = four*fd1*(drr1*ur1 + drz1*uz1) + &
                ur1*rhoi + ur2 + uz2
           fs6(ihli) = four*fd1*(drr1*vr1 + drz1*vz1) + &
                vr1*rhoi + vr2 + vz2 - (vr12 +vz12)/wv1
        End Do   !ihs
     End Do   !il
     !
     ! Associated (z,r)-weights
     Do il = 1,ngl
        Do ih = 1,ngh
           i = ih + (il-1)*ngh
           wdcor(i) = pi*wh(ih)*wl(il)*bz*bp*bp/fd(i)
           wdcori(i)=one/wdcor(i)
        End Do
     End Do
     !
     If(lpr) Then
        Do iw=lout,lfile
           Write(iw,*)
           If(ilst1.Eq.0) Then
              Write(iw,*) ' ### HO case: wdcor charged'
           Else
              Write(iw,*) ' ### THO case: wdcor charged'
           End If
           Write(iw,*)
        End Do
     End If
     Return
   End Subroutine f01234
   !=======================================================================
   !
   !=======================================================================
  Subroutine thofun(key,r,f,f1,f2,f3,fj,lpr,units)
     !---------------------------------------------------------------------
     ! Calculates LST-function 'f' its derivatives 'f1,f2,f3'
     ! at the point 'r' (all dimensionless).
     !---------------------------------------------------------------------
     Use HFBTHO
     Implicit None
     Logical :: lpr,units
     Integer(ipr) :: key,msw
     Integer(ipr) :: it,iter,ir,iqq,irmax,irmsit,immho,imm1,&
          imm2,immm,immmax,imm3
     Real(pr) ::  r,f,f1,f2,f3,fj,toto
     Real(pr), Allocatable :: dsx(:),dsy(:),dsyT(:),dsyi(:),dsyii(:),dsy1(:),         &
                              dsy1i(:),spb0(:),spc0(:),spd0(:),spbi(:),spci(:),spdi(:)
     Real(pr) :: h,hhb,pihhb,c00,snorm,snorm1,assm,asm1,asm2,asm3,               &
                 rmsit,rmmho,z1,z10,aaa,bex,rend,fj1,fj2,fj3,                    &
                 s,s1,sN,sP,sT,qq,qqup,qqdn,zqq,zqqi,df,zfj1,zfj1i,fjb,aa,bb,yyy,&
                 sqsq,rmmm,rmmmb0,z1mmm,rmm1,rmm1b0,z1mm1,rmm2,rmm2b0,z1mm2,     &
                 rmmmax,z1mmmax,rmmmaxb0,rmmx,z1mmx,z1mmxx,alaex,aldsy1,decay2
     Real(pr) :: denm1(2),denm2(2),rdenm(2)
     Real(pr) :: epsf=1.0d-14,epsdsy=1.0d-16,epsnorm
     Complex(pr) :: yyy1,bbb1,aac !for the e3rd order equation
     !

     ! Adjustable parameters
     ! Rend=40.0_pr
     epsnorm=0.01_pr

     ! ===========================
     ! KEY=0 INITIAL CALCULATIONS
     ! ===========================
     qq = r
     If(key.Eq.0.And.ilst.Le.0) Then
        write(*,*)
        write(*,*) ' LST transformation...'
        !
        ! steps
        h = 0.01_pr;
        hhb = b0*h; If(units) hhb = h
        pihhb = 4.0_pr*pi*hhb
        !
        ! correct density asymptotic
        !
        !================================
        ! Test neutron/proton asymptotic
        !================================
        Do it=1,2
           ! Neutron/proton density decay constant
           itass=it; decay=ass(itass); decay2=decay**2;
           rmsit=rms(itass)+one; irmsit=Int(rmsit/hhb)
           bb=Real((itass-1)*npr(2),Kind=pr)*Sqrt(1.440_pr)/hb0
           ! correct density Rend
           Rend=40.0_pr; irmax=Int(Rend/hhb)
           ! Deallocate/Allocate
           if(Allocated(dsx)) Deallocate(dsx,dsy,dsyT,dsyi,dsyii,dsy1,dsy1i,&
                                         spb0,spc0,spd0,spbi,spci,spdi)
           Allocate(dsx(irmax),dsy(irmax),dsyT(irmax),dsyi(irmax),dsyii(irmax),  &
                    dsy1(irmax),dsy1i(irmax),spb0(irmax),spc0(irmax),spd0(irmax),&
                    spbi(irmax),spci(irmax),spdi(irmax))
           ! HFB+HO_{L=0} density 'dsy' and its normalization
           ! integral 'dsyi' at points 'dsx' with step 'hhb=h*b0'
           ! up to the point where 'dsy*dsx*dsx < epsdsy'
           msw=0; snorm=zero
           Do While(Abs(snorm-Real(npr(itass),Kind=pr)).Gt.epsnorm.And.msw.Lt.25)
              msw=msw+6 ! increase for good norm of HFB+HO_{L=0}
              itass = - itass; s1 = zero
              Do ir=1,irmax
                 rmmho=hhb*Real(ir,Kind=pr); immho=ir
                ! L=0 component of density for isospin it (s) and isospin 1-it (sT)
                 Call densitr(itass,rmmho,sN,sP,msw)
                 if(itass.eq.1) then
                    s=sN; sT=sP
                 else
                    s=sP; sT=sN
                 End If
                 s=s*Dnfactor(itass)
                 z1=s*rmmho**2; s1=s1+z1
                 ! density dsy(ir) at point ir and its integral over r dsyi(ir)
                 ! up to that point
                 dsyT(ir)=sT; dsy(ir)=s; dsyi(ir)=hhb*s1
                 immho=ir !up to the point immho
                 If(z1.Lt.epsdsy) Exit
              End Do
              snorm=pihhb*s1  ! HFB+HO_{L=0} norm
           End Do
           ! dsy: density, spb0: first derivative with respect to r
           Call deri(hhb,immho,dsy,spb0)
           ! MIN: Find 'rmm1', the first minimum of Ln(HFB+HO_{L=0})'
           z10 = 1.0d10
           Do ir=irmsit,immho-5
              denm1(it)=dsy(ir); denm2(it)=dsyT(ir)
              z1 = spb0(ir)/dsy(ir)
              If(z1.Le.z10) Then
                 imm1 = ir; rmm1 = hhb*Real(ir,Kind=pr);  z1mm1= z1; Else; Exit
              End If
              z10 = z1
           End Do
           rdenm(itass) = rmm1/b0
           If(units) rdenm(itass) = rmm1
           ! no minimum of Ln(HFB+HO_{L=0})'
           If(rmm1.Ge.hhb*Real(immho-5,Kind=pr)) Then
              Write(*,*)
              Write(*,*) '#####################################'
              Write(*,*) 'Please increase Nsh NB!!!(NO THO RUN)'
              Write(*,*) '#####################################'
              Write(*,*)
              iasswrong(itass)=-1
              Stop
              !If(lpr) Then
              !   Open(1110,file='dat0.dat')
              !   Write(1110,*) ' r rhoh Log(rhoh)'
              !   Do ir=5,immho-5
              !      Write(1110,'(14(4x,e13.6))') hhb*Real(ir),dsy(ir),spb0(ir)/dsy(ir)
              !   End Do
              !   Close(1110)
              !End If
              !Return
           End If
        End do
        ! Asymptotics - denm1: density for isospin it, denm2: density for isospin 1-it
        If((denm1(1)-denm2(1))*(denm2(2)-denm1(2)).Le.zero) Then
           itass=1; if(ass(1).gt.ass(2)) itass=2      ! mismatch: use old asymptotic (lower decay)
        Else
           itass=2; If(denm1(1).gt.denm2(1)) itass=1  ! use new asymptotic (higher density)
        End If
        !
        iasswrong(3)=0
        If(iasswrong(itass).ne.0) iasswrong(3)=iasswrong(itass)  ! wrong assymptotic => reinforce HO results
        !
        Write(*,*) '                          min.point         neutron density       proton density'
        Write(*,*) '  1. Neutron min.point ',rdenm(1),denm1(1),denm2(1)
        Write(*,*) '  2. Protons min.point ',rdenm(2),denm2(2),denm1(2)
        Write(*,*) '     Neutron/Proton decay',ass(1),'/',ass(2)
        Write(*,*) '     Chosen Case=',itass
        !
        !===================
        ! Actual asymptotic
        !===================
        ! neutron/proton density decay constant
        decay = ass(itass); decay2=decay**2;
        rmsit= rms(itass)+one;irmsit=Int(rmsit/hhb)
        bb = Real((itass-1)*npr(2),Kind=pr)*Sqrt(1.440_pr)/hb0
        ! correct density Rend
        Rend  = 40.0_pr; irmax = Int(Rend/hhb)
        ! Deallocate/Allocate
        If(Allocated(dsx)) Deallocate(dsx,dsy,dsyT,dsyi,dsyii,dsy1,dsy1i,&
                                      spb0,spc0,spd0,spbi,spci,spdi)
        Allocate(dsx(irmax),dsy(irmax),dsyi(irmax),dsyii(irmax),   &
                 dsy1(irmax),dsy1i(irmax),spb0(irmax),spc0(irmax), &
                 spd0(irmax),spbi(irmax),spci(irmax),spdi(irmax))
        ! HFB+HO_{L=0} density 'dsy' and its normalization
        ! integral 'dsyi' at points 'dsx' with step 'hhb=h*b0'
        ! up to the point where 'dsy*dsx*dsx < epsdsy'
        msw = 0; snorm=zero
        Do While(Abs(snorm-Real(npr(itass),Kind=pr)).Gt.0.01.And.msw.Lt.25)
           msw = msw + 6 ! increase for good norm of HFB+HO_{L=0}
           itass = - itass; s1 = zero
           Do ir=1,irmax
              rmmho = hhb*Real(ir,Kind=pr); immho = ir
              Call densitr(itass,rmmho,sN,sP,msw)
              if(itass.eq.1) then
                 s=sN; sT=sP
              else
                 s=sP; sT=sN
              End If
              s=s*Dnfactor(itass)
              z1 = s*rmmho**2; s1 = s1 + z1
              dsy(ir)= s; dsyi(ir)= hhb*s1 !p-ho density and its integral
              immho = ir                   !up to the point immho
              If(z1.Lt.epsdsy) Exit
           End Do
           snorm = pihhb*s1  ! HFB+HO_{L=0} norm
        End Do
        Call deri(hhb,immho,dsy,spb0)
        ! MIN: Find 'rmm1', the first minimun of Ln(HFB+HO_{L=0})'
        z10 = 1.0d10
        Do ir=irmsit,immho-5
           z1 = spb0(ir)/dsy(ir)
           If(z1.Le.z10) Then
              imm1 = ir; rmm1 = hhb*Real(ir,Kind=pr);  z1mm1= z1; Else; Exit
           End If
           z10 = z1
        End Do
        rmm1b0 = rmm1/b0
        If(units) rmm1b0 = rmm1
        ! MAX: Find 'rmmmax', the first maximum of ln(HFB+HO_{L=0})'
        z10 = z1mm1
        Do ir=imm1,immho-5
           z1 = spb0(ir)/dsy(ir)
           If(z1.Ge.z10) Then
              immmax = ir; rmmmax = hhb*Real(ir,Kind=pr);  z1mmmax= z1; Else; Exit
           End If
           z10 = z1
        End Do
        rmmmaxb0 = rmmmax/b0
        If(units) rmmmaxb0 = rmmmax
        !
        ! END: Find 'rmm2', the last point of ln(HFB+HO_{L=0}) at the level of the first minimum
        z10 = z1mm1
        Do ir=immmax,immho-5
           z1 = spb0(ir)/dsy(ir)
           If(z1.Ge.z10) Then
              imm2 = ir; rmm2 = hhb*Real(ir,Kind=pr);  z1mm2= z1
           End If
        End Do
        rmm2b0 = rmm2/b0
        If(units) rmm2b0 = rmm2
        !
        ! MID: Find 'rmmm', the MinMaX mid point
        z10 = half*(z1mmmax+z1mm1)
        Do ir=imm1,immho
           z1 = spb0(ir)/dsy(ir)
           If(z1.Ge.z10) Exit
           immm = ir; rmmm = hhb*Real(ir,Kind=pr);  z1mmm= z1
        End Do
        rmmmb0 = rmmm/b0
        If(units) rmmmb0 = rmmm
        ! -------------------------------------------------------------
        ! Important points required:
        ! Minimum point         'rmm1'   and its log.density 'z1mm1'
        ! First maximum         'rmmmax' and its log.density 'z1mmmax'
        ! Last acceptable point 'rmm2'   and its log.density 'z1mm2'
        ! Mid point             'rmmm'   and its log.density 'z1mmm'
        ! -------------------------------------------------------------
        ! fit 'aa' from the mid match point 'rmmm','z1mmm'
        If(z1mmm.Ge.-decay) z1mmm= (-decay+z1mm1)/two   !just in case
        ! the 3rd order equation
        sqsq = rmmm*(two*bb+decay2*rmmm)
        bbb1 = one + rmmm*z1mmm
        yyy1 =(2.0_pr*bbb1**6 + 18.0_pr*bbb1**3*sqsq + 27.0_pr*sqsq**2 + &
             3.0_pr*Sqrt(3.0_pr)*sqsq**1.5_pr*Sqrt(4.0_pr*bbb1**3 + &
             27.0_pr*sqsq))**(1.0_pr/3.0_pr)
        aa = Real((2.*bbb1**2*yyy1 + 2.**(2.0_pr/3.0_pr)*yyy1**2 - &
             2.*2.**(1.0_pr/3.0_pr)*(-bbb1**4 - &
             6.0_pr*bbb1*sqsq))/(6.*yyy1),Kind=pr)
        aa = (-4.0_pr*bb*rmmm - decay2*rmmm**2 + aa)*0.250_pr
        !write(*,*)  ' aa= ',aa !If(aa.Le.zero) aa   = zero  ! in this case take l=0
        !
        ! matching logder at Rmin='rmm1' and 'rmmx'
        ! log density for rmm1<r<rmmx is z1mm1 + aaa*(r-rmm1)**2/r**assm
        ! rmmx  = rmm2     !rmm2 is last acceptable point
        rmmx  = rmmmax      !rmmmax is the firts maximum point
        sqsq  = Sqrt(decay2 + (4.0_pr*(aa + bb*rmmx))/rmmx**2)
        z1mmx=-one/rmmx - sqsq - (two*bb + decay2*rmmx)/&
             (4.0_pr*aa + rmmx*(4.0_pr*bb + decay2*rmmx))
        z1mmxx=(two*(8.0_pr*aa**2 + 16.0_pr*aa*bb*rmmx + &
             12.0_pr*bb**2*rmmx**2 + two*aa*decay2*rmmx**2 + &
             6.0_pr*bb*decay2*rmmx**3 + decay2**2*rmmx**4 + &
             rmmx*(two*aa + bb*rmmx)*sqsq*(4.0_pr*aa + &
             rmmx*(4.0_pr*bb + decay2*rmmx))))/(rmmx**2*(4.0_pr*aa + &
             rmmx*(4.0_pr*bb + decay2*rmmx))**2)
        assm  = rmmx*(two/(rmmx-rmm1)-z1mmxx/(z1mmx-z1mm1))
        asm1  =  one-assm; asm2= two-assm; asm3= three-assm
        aaa   = -rmmx*rmmx**assm*z1mmxx/ &
             ((rmm1-rmmx)*(assm*(rmm1-rmmx)+two*rmmx))
        alaex = Log(dsy(imm1))-&
             (z1mm1*rmm1 +two*aaa*rmm1**asm3/(asm3*asm2*asm1))
        ! correct density 'dsy1'
        Do ir=1,irmax
           qq = hhb*Real(ir,Kind=pr); dsx(ir) = qq
           If(ir.Le.imm1) Then
              ! inner region (r < rmm1)
              dsy1(ir) = dsy(ir)
           Else
              sqsq=Sqrt(decay2 + 4.0_pr*(aa/qq**2 + bb/qq))
              yyy=-qq*sqsq - Log(qq*qq*sqsq) - &
                   (two*bb/decay)*Log((two*bb + decay2*qq+decay*qq*sqsq)/decay)
              ! take complex in the case of negative aa
              aac=aa
              toto=two*Sqrt(aac)*&
                   Log((two*aac+bb*qq+Sqrt(aac)*qq*sqsq)/(two*aac**1.50_pr*qq))
              yyy=yyy+Cmplx(toto,0.0_pr)
              If(qq.Le.rmmx) Then
                 ! region (rmm1 < r < rmmx)
                 aldsy1 = alaex+z1mm1*qq+aaa*qq**asm1*&
                      (asm3*asm2*rmm1**2-two*asm3*asm1*rmm1*qq+asm2*asm1*qq*qq)/ &
                      (asm3*asm2*asm1)
                 dsy1(ir)=Exp(aldsy1)
                 bex = aldsy1-yyy
              Else
                 ! region (r > rmmx)
                 dsy1(ir) = Exp(bex + yyy)
              End If
           End If
        End Do
        ! correct density norm
        snorm1 = pihhb*Sum(dsy1*dsx*dsx)
        ! normalized correct density 'dsy1'
        dsy1   = snorm*dsy1/snorm1
        ! zero constant
        c00 = (dsy1(1)/dsy(1))**(1.0_pr/3.0_pr)
        ! splining correct density and its integral
        s1 = zero
        Do ir=1,irmax
           s1 = s1 +  dsy1(ir)*dsx(ir)**2
           dsy1i(ir) = hhb*s1
        End Do
        !
        ! correct density dsy1 and its integral dsy1i known up to irmax
        Call csplin(irmax,dsx,dsy1 ,spb0,spc0,spd0)
        Call csplin(irmax,dsx,dsy1i,spbi,spci,spdi)
        !
        ! print 'dat1.dat' with HFB+HO and 'correct' densities
        ! and their Log derivatives at lpr=.true.
        If(lpr) Then
           Open(1110,file='density.dat')
           !Write(1110,*) ' r rhoh rhoc Log(rhoh)'' Log(rhoc)'' '
           !Do ir=5,immho-5
           Do ir=5,irmax-5
              ! ho density derivative
              s =(45.0_pr*( dsy(ir+1)-dsy(ir-1))-9.0_pr*&
                   (dsy(ir+2)-dsy(ir-2))+dsy(ir+3)-dsy(ir-3))/(60.0_pr*hhb)/dsy(ir)
              !correct density derivative
              s1=(45.0_pr*(dsy1(ir+1)-dsy1(ir-1))-9.0_pr*(dsy1(ir+2)-&
                   dsy1(ir-2))+dsy1(ir+3)-dsy1(ir-3))/(60.0_pr*hhb)/dsy1(ir)
              Write(1110,'(2(1x,e13.6))') dsx(ir),dsy1(ir)
           End Do
           Close(1110)
           !Open(1111,file='dat1.dat')
           !Write(1111,*) ' Dimensionless_qq  Invers_f Invers_f1 Invers_f2 Invers_f3 '
        End If
        !
        ! =======================================
        ! Calculations at given dimensionless 'qq'
        ! =======================================
        ! f(R->0) = c00*R, therefore Invers_f(R)=R/c00
        fj = h/c00
        bmm3= zero; z1 = zero; ir=0
        Do iqq=1,immho
           qq = Real(iqq,Kind=pr)*h
           ! HFB+HO density and integral at 'b0*qq' or 'qq'
           zqq = dsy(iqq); zqqi = dsyi(iqq)
           ! Iterations to find 'fj = Invers_f(qq)'
           ! NB! f[Invers_f(qq)]=qq
           iter= 0; df= 0.00010_pr
           Do While(Abs(df).Ge.epsf.And.iter.Le.500)
              iter = iter + 1
              fjb = fj*b0; If(units) fjb=fj
              Call cseval(irmax,fjb,dsx,dsy1 ,spb0,spc0,spd0,zfj1 )
              Call cseval(irmax,fjb,dsx,dsy1i,spbi,spci,spdi,zfj1i)
              qqup = (Log(zfj1i/zqqi))*zfj1i;
              qqdn = zfj1*b0*fjb**2; If(units) qqdn=zfj1*fjb**2
              ! Secant & Newton
              If(zfj1i.Le.zqqi.And.df.Le.0.0_pr) df=-half*df
              If(zfj1i.Gt.zqqi.And.df.Gt.0.0_pr) df=-half*df
              If(Abs(qqdn).Gt.Abs(qqup).And.iter.Le.20) df= - qqup/qqdn
              fj = fj + df
           End Do
           fj1 = (zqq*qq*qq)/(zfj1*fj*fj)
           fj2 = (fj1 - z1)/h; z1 =fj1
           If(qq.Gt.rmm2b0) Then
              If(fj1.Ge.bmm3) Then
                 bmm3 = fj1; Else; Exit
              End If
           End If
           dsy(iqq)   = fj
           dsyi(iqq)  = fj1
           dsyii(iqq) = fj2
           iqqmax     = iqq
        End Do
        imm3 = iqqmax-50
        rmm3 = Real(imm3,Kind=pr)*h
        amm3 = dsy(imm3)
        bmm3 = dsyi(imm3)
        cmm3 = dsyii(imm3)
        ! second and third derivatives up to 'iqqmax''
        Call deri(h,iqqmax,dsyi,spb0)
        Call deri(h,iqqmax,spb0,spc0)
        !
        If(Allocated(fdsx)) Deallocate(fdsx,fdsy,fdsy1,fdsy2,fdsy3,fspb0,fspc0,  &
                                       fspd0,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,&
                                       fspb3,fspc3,fspd3)
        Allocate(fdsx(iqqmax),fdsy(iqqmax),fdsy1(iqqmax),fdsy2(iqqmax),  &
                 fdsy3(iqqmax),fspb0(iqqmax),fspc0(iqqmax),fspd0(iqqmax),&
                 fspb1(iqqmax),fspc1(iqqmax),fspd1(iqqmax),fspb2(iqqmax),&
                 fspc2(iqqmax),fspd2(iqqmax),fspb3(iqqmax),fspc3(iqqmax),&
                 fspd3(iqqmax))
        Do iqq=1,iqqmax
           fdsx(iqq)  = Real(iqq,Kind=pr)*h
           fdsy(iqq)  = dsy(iqq)
           fdsy1(iqq) = dsyi(iqq)
           fdsy2(iqq) = spb0(iqq)
           fdsy3(iqq) = spc0(iqq)
           !
           ! print 'dat1.dat' with fj=Inverse_f(qq) and its derivatives
           ! fj1..3 at no smoothing when lpr=.true.
           !If(lpr) Then
           !   Write(1111,'(14(4x,e13.6))') fdsx(iqq)*b0,fdsy(iqq),&
           !        fdsy1(iqq),fdsy2(iqq),fdsy3(iqq)
           !End If
        End Do
        !
        If(Allocated(dsx)) Deallocate(dsx,dsy,dsyi,dsyii,dsy1,dsy1i,&
             spb0,spc0,spd0,spbi,spci,spdi)
        !
        Call csplin(iqqmax,fdsx,fdsy ,fspb0,fspc0,fspd0)
        Call csplin(iqqmax,fdsx,fdsy1,fspb1,fspc1,fspd1)
        Call csplin(iqqmax,fdsx,fdsy2,fspb2,fspc2,fspd2)
        Call csplin(iqqmax,fdsx,fdsy3,fspb3,fspc3,fspd3)
        !
        Do ir=lout,lfile
           Write(ir,*)
           Write(ir,*) ' Legendre points = ',msw
           Write(ir,*) ' b0, decay=        ',b0,decay
           Write(ir,*) ' h, hhb=           ',h,hhb
           Write(ir,*) ' rms, rmsit=       ',rms(itass),rmsit
           Write(ir,*) ' Rend,   irmax=    ',Rend,irmax
           Write(ir,*) ' HORend, immho=    ',rmmho,immho
           Write(ir,*) ' snorm,snorm1=     ',snorm,snorm1
           Write(ir,*) ' snorm/snorm1=     ',snorm/snorm1
           Write(ir,*) ' min:rmm1,/b0=     ',rmm1,rmm1b0
           Write(ir,*) ' max:rmmmax,/b0=   ',rmmmax,rmmmaxb0
           Write(ir,*) ' last:rmm2,/b0=    ',rmm2,rmm2b0
           Write(ir,*) ' num:rmm3*b0,rmm3= ',rmm3*b0,rmm3
           Write(ir,*) ' rmmho, rmmho/b0=  ',rmmho,rmmho/b0
           Write(ir,*) ' alaex,bex=        ',alaex,bex
           Write(ir,*) ' aa,bb=            ',aa,bb
           Write(ir,*) ' L_eff=            ',(Sqrt(one + 4.0_pr*aac)-one)/two
           Write(ir,*) ' amm3, bmm3=       ',amm3,bmm3
           Write(ir,*) ' cmm3, one/c00=    ',cmm3,one/c00
           Write(ir,*) ' z1mm1,z1mmm=      ',z1mm1,z1mmm
           Write(ir,*)
        End Do
        !! print 'dat2..4.dat' after when lpr=.true.
        !! dat2.dat: 'r=b0*qq',correct density
        !! dat3.dat   qq, Invers_f,Invers_f1...3
        !! dat4.dat   qq, f,f1,f2,f3; 'qq' are the Gauss points
        !If(lpr) Then
        !   Close(1111)
        !   Open(1112,file='dat2.dat')
        !   Open(1113,file='dat3.dat')
        !   Open(1114,file='dat4.dat')
        !   Write(1112,*) ' r  den_correct'
        !   Write(1113,*) ' qq Invers_f Invers_f1 Invers_f2 Invers_f3'
        !   Write(1114,*) ' qq f f1 f2 f3'
        !   Do iqq=1,ngh
        !      Do ir=1,ngl
        !         qq   = Sqrt(xh(iqq)**2+xl(ir))
        !         Call densitr(itass,b0*qq,sN,sP,msw)
        !         if(itass.eq.1) then
        !            s=sN; sT=sP
        !         else
        !            s=sP; sT=sN
        !         End If
        !         s=s*Dnfactor(itass)
        !         If(qq.Le.rmm3) Then
        !            Call cseval(iqqmax,qq,fdsx,fdsy ,fspb0,fspc0,fspd0,fj)
        !            Call cseval(iqqmax,qq,fdsx,fdsy1,fspb1,fspc1,fspd1,fj1)
        !            Call cseval(iqqmax,qq,fdsx,fdsy2,fspb2,fspc2,fspd2,fj2)
        !            Call cseval(iqqmax,qq,fdsx,fdsy3,fspb3,fspc3,fspd3,fj3)
        !         Else
        !            fj  = amm3+bmm3*(qq-rmm3)+cmm3*(qq-rmm3)**2/two
        !            fj1 = bmm3+cmm3*(qq-rmm3)
        !            fj2 = cmm3; fj3 = zero
        !         End If
        !         f  = qq; f1 = one/fj1; f2 =-fj2*f1**3
        !         f3 = three*fj2**2*f1**5-fj3*f1**4
        !         s = s*qq*qq/(fj*fj*fj1)
        !         Write(1112,'(14(4x,e13.6))') b0*fj,s
        !         Write(1113,'(14(4x,e13.6))') qq,fj,fj1,fj2,fj3
        !         Write(1114,'(14(4x,e13.6))') qq,f,f1,f2,f3
        !      End Do
        !   End Do
        !   Close(1112); Close(1113);
        !   Close(1114)
        !End If
        ! =========================
     Else  !KEY=1 CALCULATIONS
        ! =========================
        !
        ! Calculations of Invers_f(qq),_f'(qq),_f''(qq),_f'''(qq)
        If(qq.Le.rmm3) Then
           Call cseval(iqqmax,qq,fdsx,fdsy ,fspb0,fspc0,fspd0,fj)
           Call cseval(iqqmax,qq,fdsx,fdsy1,fspb1,fspc1,fspd1,fj1)
           Call cseval(iqqmax,qq,fdsx,fdsy2,fspb2,fspc2,fspd2,fj2)
           Call cseval(iqqmax,qq,fdsx,fdsy3,fspb3,fspc3,fspd3,fj3)
        Else
           fj  = amm3+bmm3*(qq-rmm3)+cmm3*(qq-rmm3)**2/two
           fj1 = bmm3+cmm3*(qq-rmm3); fj2 = cmm3; fj3 = zero
        End If
        ! Calculations of f(fj),f'(fj),f''(fj),f'''(fj)
        f  = qq; f1 = one/fj1; f2 =-fj2*f1**3
        f3 = three*fj2**2*f1**5-fj3*f1**4
     End If
     Return
   End Subroutine thofun
   !=======================================================================
   !
   !=======================================================================
   Subroutine densitr(it,xr,yr,yrP,msw)
     !---------------------------------------------------------------------
     ! Calculates Legendre decomposition of neutron(proton) 'it=1(2)'
     ! HFB+HO_{L=0}(r) density 'yr' at point 'xr' (in fm)
     !---------------------------------------------------------------------
     Use HFBTHO
     Implicit None
     Integer(ipr) :: it,msw
     Integer(ipr), Save :: iw,ik,il,i0,i02,jk,nsa,nsb,nrb,ny,nyy,ib,nd,&
          n1,n2,n1n2nd,ibit,ibitnb,nzb,mlb,ngh1,ngl1
     ! msw=20 test for protons at the crazy case of U 212 120 92
     ! msw=3  s,a1= 107.665062 80.82396  s/s1= 1.33209
     ! msw=6  s,s1=  88.973061 80.99903  s/s1= 1.09844
     ! msw=12 s,s1=  92.025592 80.99595  s/s1= 1.13617
     ! msw=18 s,s1=  92.000017 80.99595  s/s1= 1.13585
     ! msw=24 s,s1=  92.000010 80.99595  s/s1= 1.13585
     ! Taken up to msw=24
     Real(pr), Allocatable, Save :: xmw(:),yi(:,:)
     Real(pr) :: phy(msw,nzrlx),anl(msw,msw),yl(msw,msw)
     Real(pr), Save :: sl,w,hw,ct2,s,frit,fritP,wdcorin,bzi,bri,ct,st,z,t
     Real(pr) :: xr,yr,yrp
     !
     ngh1=ngh+1; ngl1=ngl+1
     If(it.Lt.0) Then
        it = -it
        If(Allocated(xmw)) Deallocate(xmw,yi)
        Allocate(xmw(msw),yi(msw,msw))
        wdcorin=one/Sqrt(pi*bz*bpp); bzi=one/bz; bri=one/bp
        ! 'msw' mesh-points in 'angle' space
        xmw(1)=zero; hw=half*pi/Real(msw-1,Kind=pr)
        Do il=2,msw
           xmw(il)=hw*Real(il-1,Kind=pr)
        End Do
        ! coefficients for the L-decomposition
        sl=4.0_pr
        Do il=1,msw
           sl=0.250_pr*sl
           Do ny=1,il
              anl(ny,il) = iv(il-ny)*fak(2*(il+ny-2))*&
                   fi(il-ny)*fi(il+ny-2)*fi(2*ny-2)*sl
           End Do
        End Do
        Do iw=1,msw
           w = xmw(iw); ct2 = Cos(w)**2
           Do il=1,msw
              yi(iw,il) = zero; s = zero
              Do nyy=1,il
                 ny = il + 1 - nyy; s = s*ct2 + anl(ny,il)
              End Do
              yl(iw,il) = s*sq(4*il-3)
           End Do
           yi(iw,iw) = one
        End Do
        Call lingd(msw,msw,msw,msw,yl,yi,s,il)
     End If
     ! 'xr/yr' calculations
     Do iw=1,msw
        w = xmw(iw); ct = Cos(w); st = Sin(w)
        z = ct*xr*bzi; t = (st*xr*bri)**2
        Call gaupolr(z,t)
        ik = 0
        Do ib = 1,nb
           nd= id(ib); i0= ia(ib)
           Do n2 = 1,nd
              ik = ik + 1; i02= i0 + n2
              nzb= nz(i02);  nrb= nr(i02); mlb = nl(i02)
              phy(iw,ik) = qh(nzb,ngh1)*ql(nrb,mlb,ngl1)*wdcorin
           End Do
        End Do
     End Do
     ! 'yr' over the blocks
     yr=zero; yrP=zero; ik = 0
     Do ib = 1,nb
        nd= id(ib); i0= ia(ib)
        Do n2 = 1,nd
           jk= ik; ik= ik + 1; i02=i0 + n2; nsb= ns(i02)
           Do n1 = n2,nd
              jk= jk + 1; i02=i0 + n1; nsa= ns(i02)
              If(nsa.Eq.nsb) Then
                 ibit=ib; ibitnb= ib+nbx; n1n2nd= n1+(n2-1)*nd
                 frit = rk(n1n2nd,ibit); fritP = rk(n1n2nd,ibitnb)
                 If(n1.Ne.n2) then
                    frit = two*frit; fritP = two*fritP
                 End If
                 s = zero
                 Do iw=1,msw
                    s = s + yi(1,iw)*phy(iw,ik)*phy(iw,jk)
                 End Do
                 yr=yr+frit*s; yrP=yrP+fritP*s
              End If
           End Do !n2
        End Do !n1
     End Do !ib
     !
     Return
   End Subroutine densitr
   !=======================================================================
   !
   !=======================================================================
   Subroutine gaupolr(z,x)
     !---------------------------------------------------------------------
     ! see 'gaupol'
     !---------------------------------------------------------------------
     Use HFBTHO
     Implicit None
     Real(pr) :: z,x
     Real(pr) :: w0,w00,w4pii,dsq,d1,d2
     Integer(ipr) :: N,L,NGH1,NGL1
     !
     NGH1=NGH+1; NGL1=NGL+1
     W4PII = PI**(-0.250_pr); W0 = W4PII*Exp(-HALF*Z*Z)
     ! W0 = W0*SQRT(Z) NOT MULTIPLIED BY WDCOR
     QH(0,NGH1) = W0; QH(1,NGH1)= SQ(2)*W0*Z
     Do N = 2,NZM
        QH(N,NGH1) = SQI(N)*(SQ(2)*Z*QH(N-1,NGH1)-SQ(N-1)*QH(N-2,NGH1))
     End Do
     W00 = SQ(2)*Exp(-HALF*X)
     Do L = 0,NLM
        If(L.Eq.0) Then
           W0 = W00*Sqrt(HALF)
        Else
           W0 = W00*Sqrt(HALF*X**L)
        End If
        QL(0,L,NGL1) = WFI(L)*W0; QL(1,L,NGL1) = (Real(L+1,Kind=pr)-X)*WFI(L+1)*W0
        Do N = 2,NRM
           DSQ = SQ(N)*SQ(N+L); D1= Real(2*N + L - 1,Kind=pr) - X
           D2  = SQ(N-1)*SQ(N-1+L)
           QL(N,L,NGL1) = (D1*QL(N-1,L,NGL1)-D2*QL(N-2,L,NGL1))/DSQ
        End Do
     End Do
     Return
   End Subroutine gaupolr
End Module HFBTHO_THO
!==================================================================================================================================
!#END THO MODULE
!==================================================================================================================================
!#START EllipticIntegral MODULE
!==================================================================================================================================
Module EllipticIntegral
  !--------------------------------------------------------------------------------------
  ! This module provides a routine to compute the complete elliptic integral of
  ! the second kind.
  !
  ! Reference:
  !    Fukushima, T.,
  !    Fast Computation of Complete Elliptic Integrals and Jacobian Elliptic Functions,
  !    Celest. Mech. Dyn. Astron., 105, 305-328 (2009b)
  !--------------------------------------------------------------------------------------

  Use HFBTHO_utilities

  Implicit None

Contains

  Real(pr) Function CompleteEllipticFunction_2nd(x)
    Use HFBTHO_utilities
    Implicit None
    Real(pr), INTENT(IN) :: x

    Real(pr) :: Emp,Kmp,Em,Km,qp,pi,x_eff

    pi = four*Atan(one)

    If(x.Lt.zero.Or.x.Gt.one) Stop 'Error in CompleteEllipticFunction_2nd'

    If(x.Lt.0.9_pr) Then
       Em = elliptic_small_m(x)
    Else
       Call auxiliary(x,Emp,Kmp)
       x_eff = one-x; qp = nome(x_eff)
       Km = -Log(qp)*Kmp/pi
       Em = Km + (half*pi - Emp*Km)/Kmp
    End If

    CompleteEllipticFunction_2nd = Em

  End Function CompleteEllipticFunction_2nd
  !=======================================================================
  !
  !=======================================================================
  Subroutine auxiliary(x,Emp,Kmp)
    Use HFBTHO_utilities
    Implicit None
    Real(pr), INTENT(IN) :: x
    Real(pr), INTENT(INOUT) :: Emp,Kmp

    Integer(ipr) :: JE, JK
    Parameter (JE=16, JK=16)
    Integer(ipr) :: i
    Real(pr) :: x0, x_eff
    Real(pr), Dimension(0:JE) :: Coeff_E
    Real(pr), Dimension(0:JK) :: Coeff_K

    x0 = 0.05_pr; x_eff = one - x

    Coeff_E(0) = +1.550973351780472328_pr
    Coeff_E(1) = -0.400301020103198524_pr
    Coeff_E(2) = -0.078498619442941939_pr
    Coeff_E(3) = -0.034318853117591992_pr
    Coeff_E(4) = -0.019718043317365499_pr
    Coeff_E(5) = -0.013059507731993309_pr
    Coeff_E(6) = -0.009442372874146547_pr
    Coeff_E(7) = -0.007246728512402157_pr
    Coeff_E(8) = -0.005807424012956090_pr
    Coeff_E(9) = -0.004809187786009338_pr
    Coeff_E(10)=  0.000000000000000000_pr
    Coeff_E(11)=  0.000000000000000000_pr
    Coeff_E(12)=  0.000000000000000000_pr
    Coeff_E(13)=  0.000000000000000000_pr
    Coeff_E(14)=  0.000000000000000000_pr
    Coeff_E(15)=  0.000000000000000000_pr
    Coeff_E(16)=  0.000000000000000000_pr

    Emp = 0.0_pr
    Do i=0,JE
       Emp = Emp + Coeff_E(i)*(x_eff - x0)**i
    End Do

    Coeff_K(0) = 1.591003453790792180_pr
    Coeff_K(1) = 0.416000743991786912_pr
    Coeff_K(2) = 0.245791514264103415_pr
    Coeff_K(3) = 0.179481482914906162_pr
    Coeff_K(4) = 0.144556057087555150_pr
    Coeff_K(5) = 0.123200993312427711_pr
    Coeff_K(6) = 0.108938811574293531_pr
    Coeff_K(7) = 0.098853409871592910_pr
    Coeff_K(8) = 0.091439629201749751_pr
    Coeff_K(9) = 0.085842591595413900_pr
    Coeff_K(10)= 0.081541118718303215_pr
    Coeff_K(11)= 0.000000000000000000_pr
    Coeff_K(12)= 0.000000000000000000_pr
    Coeff_K(13)= 0.000000000000000000_pr
    Coeff_K(14)= 0.000000000000000000_pr
    Coeff_K(15)= 0.000000000000000000_pr
    Coeff_K(16)= 0.000000000000000000_pr

    Kmp = 0.0_pr
    Do i=0,JK
       Kmp = Kmp + Coeff_K(i)*(x_eff - x0)**i
    End Do

  End Subroutine auxiliary
  !=======================================================================
  !
  !=======================================================================
  Real(pr) Function elliptic_small_m(x)
    Use HFBTHO_utilities
    Implicit None
    Real(pr), INTENT(IN) :: x

    Integer(ipr) :: JE
    Parameter (JE=16)
    Integer(ipr) :: i
    Real(pr) :: x0, Em
    Real(pr), Dimension(0:JE) :: Coeff_E

    If(x.Lt.0.1_pr) Then
       x0 = 0.05_pr
       Coeff_E(0) = +1.550973351780472328_pr
       Coeff_E(1) = -0.400301020103198524_pr
       Coeff_E(2) = -0.078498619442941939_pr
       Coeff_E(3) = -0.034318853117591992_pr
       Coeff_E(4) = -0.019718043317365499_pr
       Coeff_E(5) = -0.013059507731993309_pr
       Coeff_E(6) = -0.009442372874146547_pr
       Coeff_E(7) = -0.007246728512402157_pr
       Coeff_E(8) = -0.005807424012956090_pr
       Coeff_E(9) = -0.004809187786009338_pr
       Coeff_E(10)=  0.000000000000000000_pr
       Coeff_E(11)=  0.000000000000000000_pr
       Coeff_E(12)=  0.000000000000000000_pr
       Coeff_E(13)=  0.000000000000000000_pr
       Coeff_E(14)=  0.000000000000000000_pr
       Coeff_E(15)=  0.000000000000000000_pr
       Coeff_E(16)=  0.000000000000000000_pr
    End If

    If(x.Ge.0.1_pr.And.x.Lt.0.2_pr) Then
       x0 = 0.15_pr
       Coeff_E(0) = +1.510121832092819728_pr
       Coeff_E(1) = -0.417116333905867549_pr
       Coeff_E(2) = -0.090123820404774569_pr
       Coeff_E(3) = -0.043729944019084312_pr
       Coeff_E(4) = -0.027965493064761785_pr
       Coeff_E(5) = -0.020644781177568105_pr
       Coeff_E(6) = -0.016650786739707238_pr
       Coeff_E(7) = -0.014261960828842520_pr
       Coeff_E(8) = -0.012759847429264803_pr
       Coeff_E(9) = -0.011799303775587354_pr
       Coeff_E(10)= -0.011197445703074968_pr
       Coeff_E(11)=  0.000000000000000000_pr
       Coeff_E(12)=  0.000000000000000000_pr
       Coeff_E(13)=  0.000000000000000000_pr
       Coeff_E(14)=  0.000000000000000000_pr
       Coeff_E(15)=  0.000000000000000000_pr
       Coeff_E(16)=  0.000000000000000000_pr
    End If

    If(x.Ge.0.2_pr.And.x.Lt.0.3_pr) Then
       x0 = 0.25_pr
       Coeff_E(0) = +1.467462209339427155_pr
       Coeff_E(1) = -0.436576290946337775_pr
       Coeff_E(2) = -0.105155557666942554_pr
       Coeff_E(3) = -0.057371843593241730_pr
       Coeff_E(4) = -0.041391627727340220_pr
       Coeff_E(5) = -0.034527728505280841_pr
       Coeff_E(6) = -0.031495443512532783_pr
       Coeff_E(7) = -0.030527000890325277_pr
       Coeff_E(8) = -0.030916984019238900_pr
       Coeff_E(9) = -0.032371395314758122_pr
       Coeff_E(10)= -0.034789960386404158_pr
       Coeff_E(11)=  0.000000000000000000_pr
       Coeff_E(12)=  0.000000000000000000_pr
       Coeff_E(13)=  0.000000000000000000_pr
       Coeff_E(14)=  0.000000000000000000_pr
       Coeff_E(15)=  0.000000000000000000_pr
       Coeff_E(16)=  0.000000000000000000_pr
    End If

    If(x.Ge.0.3_pr.And.x.Lt.0.4_pr) Then
       x0 = 0.35_pr
       Coeff_E(0) = +1.422691133490879171_pr
       Coeff_E(1) = -0.459513519621048674_pr
       Coeff_E(2) = -0.125250539822061878_pr
       Coeff_E(3) = -0.078138545094409477_pr
       Coeff_E(4) = -0.064714278472050002_pr
       Coeff_E(5) = -0.062084339131730311_pr
       Coeff_E(6) = -0.065197032815572477_pr
       Coeff_E(7) = -0.072793895362578779_pr
       Coeff_E(8) = -0.084959075171781003_pr
       Coeff_E(9) = -0.102539850131045997_pr
       Coeff_E(10)= -0.127053585157696036_pr
       Coeff_E(11)= -0.160791120691274606_pr
       Coeff_E(12)=  0.000000000000000000_pr
       Coeff_E(13)=  0.000000000000000000_pr
       Coeff_E(14)=  0.000000000000000000_pr
       Coeff_E(15)=  0.000000000000000000_pr
       Coeff_E(16)=  0.000000000000000000_pr
    End If

    If(x.Ge.0.4_pr.And.x.Lt.0.5_pr) Then
       x0 = 0.45_pr
       Coeff_E(0) = +1.375401971871116291_pr
       Coeff_E(1) = -0.487202183273184837_pr
       Coeff_E(2) = -0.153311701348540228_pr
       Coeff_E(3) = -0.111849444917027833_pr
       Coeff_E(4) = -0.108840952523135768_pr
       Coeff_E(5) = -0.122954223120269076_pr
       Coeff_E(6) = -0.152217163962035047_pr
       Coeff_E(7) = -0.200495323642697339_pr
       Coeff_E(8) = -0.276174333067751758_pr
       Coeff_E(9) = -0.393513114304375851_pr
       Coeff_E(10)= -0.575754406027879147_pr
       Coeff_E(11)= -0.860523235727239756_pr
       Coeff_E(12)= -1.308833205758540162_pr
       Coeff_E(13)=  0.000000000000000000_pr
       Coeff_E(14)=  0.000000000000000000_pr
       Coeff_E(15)=  0.000000000000000000_pr
       Coeff_E(16)=  0.000000000000000000_pr
    End If

    If(x.Ge.0.5_pr.And.x.Lt.0.6_pr) Then
       x0 = 0.55_pr
       Coeff_E(0) = +1.325024497958230082_pr
       Coeff_E(1) = -0.521727647557566767_pr
       Coeff_E(2) = -0.194906430482126213_pr
       Coeff_E(3) = -0.171623726822011264_pr
       Coeff_E(4) = -0.202754652926419141_pr
       Coeff_E(5) = -0.278798953118534762_pr
       Coeff_E(6) = -0.420698457281005762_pr
       Coeff_E(7) = -0.675948400853106021_pr
       Coeff_E(8) = -1.136343121839229244_pr
       Coeff_E(9) = -1.976721143954398261_pr
       Coeff_E(10)= -3.531696773095722506_pr
       Coeff_E(11)= -6.446753640156048150_pr
       Coeff_E(12)= -11.97703130208884026_pr
       Coeff_E(13)=  0.000000000000000000_pr
       Coeff_E(14)=  0.000000000000000000_pr
       Coeff_E(15)=  0.000000000000000000_pr
       Coeff_E(16)=  0.000000000000000000_pr
    End If

    If(x.Ge.0.6_pr.And.x.Lt.0.7_pr) Then
       x0 = 0.65_pr
       Coeff_E(0) = +1.270707479650149744_pr
       Coeff_E(1) = -0.566839168287866583_pr
       Coeff_E(2) = -0.262160793432492598_pr
       Coeff_E(3) = -0.292244173533077419_pr
       Coeff_E(4) = -0.440397840850423189_pr
       Coeff_E(5) = -0.774947641381397458_pr
       Coeff_E(6) = -1.498870837987561088_pr
       Coeff_E(7) = -3.089708310445186667_pr
       Coeff_E(8) = -6.667595903381001064_pr
       Coeff_E(9) = -14.89436036517319078_pr
       Coeff_E(10)= -34.18120574251449024_pr
       Coeff_E(11)= -80.15895841905397306_pr
       Coeff_E(12)= -191.3489480762984920_pr
       Coeff_E(13)= -463.5938853480342030_pr
       Coeff_E(14)= -1137.380822169360061_pr
       Coeff_E(15)=  0.000000000000000000_pr
       Coeff_E(16)=  0.000000000000000000_pr
    End If

    If(x.Ge.0.7_pr.And.x.Lt.0.8_pr) Then
       x0 = 0.75_pr
       Coeff_E(0) = +1.211056027568459525_pr
       Coeff_E(1) = -0.630306413287455807_pr
       Coeff_E(2) = -0.387166409520669145_pr
       Coeff_E(3) = -0.592278235311934603_pr
       Coeff_E(4) = -1.237555584513049844_pr
       Coeff_E(5) = -3.032056661745247199_pr
       Coeff_E(6) = -8.181688221573590762_pr
       Coeff_E(7) = -23.55507217389693250_pr
       Coeff_E(8) = -71.04099935893064956_pr
       Coeff_E(9) = -221.8796853192349888_pr
       Coeff_E(10)= -712.1364793277635425_pr
       Coeff_E(11)= -2336.125331440396407_pr
       Coeff_E(12)= -7801.945954775964673_pr
       Coeff_E(13)= -26448.19586059191933_pr
       Coeff_E(14)= -90799.48341621365251_pr
       Coeff_E(15)= -315126.0406449163424_pr
       Coeff_E(16)= -1104011.344311591159_pr
    End If

    If(x.Ge.0.8_pr.And.x.Lt.0.85_pr) Then
       x0 = 0.825_pr
       Coeff_E(0) = +1.161307152196282836_pr
       Coeff_E(1) = -0.701100284555289548_pr
       Coeff_E(2) = -0.580551474465437362_pr
       Coeff_E(3) = -1.243693061077786614_pr
       Coeff_E(4) = -3.679383613496634879_pr
       Coeff_E(5) = -12.81590924337895775_pr
       Coeff_E(6) = -49.25672530759985272_pr
       Coeff_E(7) = -202.1818735434090269_pr
       Coeff_E(8) = -869.8602699308701437_pr
       Coeff_E(9) = -3877.005847313289571_pr
       Coeff_E(10)= -17761.70710170939814_pr
       Coeff_E(11)= -83182.69029154232061_pr
       Coeff_E(12)= -396650.4505013548170_pr
       Coeff_E(13)= -1920033.413682634405_pr
       Coeff_E(14)=  0.000000000000000000_pr
       Coeff_E(15)=  0.000000000000000000_pr
       Coeff_E(16)=  0.000000000000000000_pr
    End If

    If(x.Ge.0.85_pr.And.x.Lt.0.9_pr) Then
       x0 = 0.875_pr
       Coeff_E(0) = +1.124617325119752213_pr
       Coeff_E(1) = -0.770845056360909542_pr
       Coeff_E(2) = -0.844794053644911362_pr
       Coeff_E(3) = -2.490097309450394453_pr
       Coeff_E(4) = -10.23971741154384360_pr
       Coeff_E(5) = -49.74900546551479866_pr
       Coeff_E(6) = -267.0986675195705196_pr
       Coeff_E(7) = -1532.665883825229947_pr
       Coeff_E(8) = -9222.313478526091951_pr
       Coeff_E(9) = -57502.51612140314030_pr
       Coeff_E(10)= -368596.1167416106063_pr
       Coeff_E(11)= -2415611.088701091428_pr
       Coeff_E(12)= -16120097.81581656797_pr
       Coeff_E(13)= -109209938.5203089915_pr
       Coeff_E(14)= -749380758.1942496220_pr
       Coeff_E(15)= -5198725846.725541393_pr
       Coeff_E(16)= -36409256888.12139973_pr
    End If

    Em = 0.0_pr
    Do i=0,JE
       Em = Em + Coeff_E(i)*(x - x0)**i
    End Do

    elliptic_small_m = Em

  End Function elliptic_small_m
  !=======================================================================
  !
  !=======================================================================
  Real(pr) Function nome(x)
    Use HFBTHO_utilities
    Implicit None
    Real(pr), INTENT(IN) :: x
    Integer(ipr) :: Jq
    Parameter (Jq = 14)
    Integer(ipr) :: i
    Real(pr) :: qp,epsilon
    Real(pr), Dimension(1:Jq) :: Coeff_q

    epsilon = 1.D-14

    Coeff_q(1)  = 1.0_pr/16.0_pr
    Coeff_q(2)  = 1.0_pr/32.0_pr
    Coeff_q(3)  = 21.0_pr/1024.0_pr
    Coeff_q(4)  = 31.0_pr/2048.0_pr
    Coeff_q(5)  = 6257.0_pr/524288.0_pr
    Coeff_q(6)  = 10293.0_pr/1048576.0_pr
    Coeff_q(7)  = 279025.0_pr/33554432.0_pr
    Coeff_q(8)  = 483127.0_pr/67108864.0_pr
    Coeff_q(9)  = 435506703.0_pr/68719476736.0_pr
    Coeff_q(10) = 776957575.0_pr/137438953472.0_pr
    Coeff_q(11) = 22417045555.0_pr/4398046511104.0_pr
    Coeff_q(12) = 40784671953.0_pr/8796093022208.0_pr
    Coeff_q(13) = 9569130097211.0_pr/2251799813685248.0_pr
    Coeff_q(14) = 17652604545791.0_pr/4503599627370496.0_pr

    qp = 0.0_pr
    If(x.Gt.epsilon) Then
       Do i=1,Jq
          qp = qp + Coeff_q(i) * (x**i)
       End Do
    Else
       qp = epsilon
    End If

    nome = qp

  End Function nome

End Module EllipticIntegral
!==================================================================================================================================
!#END EllipticIntegral
!==================================================================================================================================
!#START bessik
!==================================================================================================================================
Module bessik

  Use HFBTHO_utilities

  Implicit None

Contains

  Function besei0(x)

    !---------------------------------------------------------------------
    !  BESEI0 evaluates the exponentially scaled Bessel I0(X) function.
    !
    !  Discussion:
    !
    !    This routine computes approximate values for the modified Bessel
    !    function of the first kind of order zero multiplied by EXP(-ABS(X)).
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    03 April 2007
    !
    !  Author:
    !
    !    Original FORTRAN77 version by William Cody.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument of the function.
    !
    !    Output, real ( kind = 8 ) BESEI0, the value of the function.
    !---------------------------------------------------------------------
    Implicit None

    Real(Kind=pr) :: besei0
    Integer(Kind=ipr) :: jint
    Real(Kind=pr) :: result,x

    jint = 2
    Call calci0(x,result,jint)
    besei0 = result

    Return
  End Function besei0
  !=======================================================================
  !
  !=======================================================================
  Function besei1(x)
    !---------------------------------------------------------------------
    ! BESEI1 evaluates the exponentially scaled Bessel I1(X) function.
    !
    !  Discussion:
    !
    !    This routine computes approximate values for the
    !    modified Bessel function of the first kind of order one
    !    multiplied by EXP(-ABS(X)).
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    03 April 2007
    !
    !  Author:
    !
    !    Original FORTRAN77 version by William Cody.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument of the function.
    !
    !    Output, real ( kind = 8 ) BESEI1, the value of the function.
    !---------------------------------------------------------------------
    Implicit none

    Real(Kind=pr) :: besei1
    Integer(Kind=ipr) :: jint
    Real(Kind=pr) :: result,x

    jint = 2
    Call calci1 ( x, result, jint )
    besei1 = result

    Return
  End Function besei1
  !=======================================================================
  !
  !=======================================================================
  subroutine calci0 ( arg, result, jint )
    !---------------------------------------------------------------------
    !
    !! CALCI0 computes various I0 Bessel functions.
    !
    !  Discussion:
    !
    !    This routine computes modified Bessel functions of the first kind
    !    and order zero, I0(X) and EXP(-ABS(X))*I0(X), for real
    !    arguments X.
    !
    !    The main computation evaluates slightly modified forms of
    !    minimax approximations generated by Blair and Edwards, Chalk
    !    River (Atomic Energy of Canada Limited) Report AECL-4928,
    !    October, 1974.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    03 April 2007
    !
    !  Author:
    !
    !    Original FORTRAN77 version by William Cody, Laura Stoltz.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) ARG, the argument.  If JINT = 1, then
    !    the argument must be less than XMAX.
    !
    !    Output, real ( kind = 8 ) RESULT, the value of the function,
    !    which depends on the input value of JINT:
    !    1, RESULT = I0(x);
    !    2, RESULT = exp(-x) * I0(x);
    !
    !    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
    !    1, I0(x);
    !    2, exp(-x) * I0(x);
    !---------------------------------------------------------------------
    Implicit None

    Real(Kind=pr) :: a,arg,b,exp40,forty
    Integer(Kind=ipr) :: i,jint
    Real(Kind=pr) :: one5,p(15),pp(8),q(5),qq(7),result,rec15
    Real(Kind=pr) :: sump,sumq,two25,x,xinf,xmax,xsmall,xx
    !  Mathematical constants
    Data one5 /15.0_pr/
    Data exp40 /2.353852668370199854d17/
    Data forty /40.0_pr/
    Data rec15 /6.6666666666666666666d-2/
    Data two25 /225.0_pr/
    !  Machine-dependent constants
    Data xsmall /5.55d-17/
    !Data xinf /1.79d308/
    Data xmax /713.986_pr/
    !  Coefficients for XSMALL <= ABS(ARG) < 15.0
    Data  p/-5.2487866627945699800d-18,-1.5982226675653184646d-14, &
            -2.6843448573468483278d-11,-3.0517226450451067446d-08, &
            -2.5172644670688975051d-05,-1.5453977791786851041d-02, &
            -7.0935347449210549190d+00,-2.4125195876041896775d+03, &
            -5.9545626019847898221d+05,-1.0313066708737980747d+08, &
            -1.1912746104985237192d+10,-8.4925101247114157499d+11, &
            -3.2940087627407749166d+13,-5.5050369673018427753d+14, &
            -2.2335582639474375249d+15/
    Data  q/-3.7277560179962773046d+03, 6.5158506418655165707d+06, &
            -6.5626560740833869295d+09, 3.7604188704092954661d+12, &
            -9.7087946179594019126d+14/
    !  Coefficients for 15.0 <= ABS(ARG)
    Data pp/-3.9843750000000000000d-01, 2.9205384596336793945d+00, &
            -2.4708469169133954315d+00, 4.7914889422856814203d-01, &
            -3.7384991926068969150d-03,-2.6801520353328635310d-03, &
             9.9168777670983678974d-05,-2.1877128189032726730d-06/
    Data qq/-3.1446690275135491500d+01, 8.5539563258012929600d+01, &
            -6.0228002066743340583d+01, 1.3982595353892851542d+01, &
            -1.1151759188741312645d+00, 3.2547697594819615062d-02, &
            -5.5194330231005480228d-04/

    x = Abs(arg)

    If(x <xsmall) Then

       result = one
    !
    !  XSMALL <= ABS(ARG) < 15.0.
    !
    Else If (x<one5) Then

       xx = x * x
       sump = p(1)
       Do i = 2, 15
          sump = sump * xx + p(i)
       End Do
       xx = xx - two25

       sumq = (((( xx + q(1) ) * xx + q(2) ) * xx + q(3) ) * xx + q(4) ) * xx + q(5)

       result = sump / sumq

       If(jint == 2 ) Then
          result = result * Exp(-x)
       End If

    Else If (one5<=x) Then

       If(jint==1 .and. xmax<x) Then
          result = Huge(1.0_pr)
       Else
          !
          !  15.0 <= ABS(ARG).
          !
          xx = one / x - rec15

          sump = (((((( pp(1) * xx + pp(2) ) * xx + pp(3) ) * xx + pp(4) ) &
               * xx + pp(5) ) * xx + pp(6) ) * xx + pp(7) ) * xx + pp(8)

          sumq = (((((( xx + qq(1) ) * xx + qq(2) ) * xx + qq(3) ) * xx + qq(4) ) &
                      * xx + qq(5) ) * xx + qq(6) ) * xx + qq(7)

          result = sump / sumq

          If(jint==2) Then
             result = (result - pp(1)) / Sqrt(x)
          Else
             !
             !  Calculation reformulated to avoid premature overflow.
             !
             If(x.le.(xmax-one5)) Then
                a = Exp(x)
                b = one
             Else
                a = Exp(x-forty)
                b = exp40
             End If

             result = ( (result*a - pp(1)*a) / Sqrt(x) ) * b

          End If

       End If

    End If

    Return
  End Subroutine calci0
  !===============================================================================================
  !
  !===============================================================================================
  Subroutine calci1(arg,result,jint)
    !---------------------------------------------------------------------
    !
    !! CALCI1 computes various I1 Bessel functions.
    !
    !  Discussion:
    !
    !    This routine computes modified Bessel functioons of the first kind
    !    and order one, I1(X) and EXP(-ABS(X))*I1(X), for real
    !    arguments X.
    !
    !    The main computation evaluates slightly modified forms of
    !    minimax approximations generated by Blair and Edwards, Chalk
    !    River (Atomic Energy of Canada Limited) Report AECL-4928,
    !    October, 1974.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    03 April 2007
    !
    !  Author:
    !
    !    Original FORTRAN77 version by William Cody, Laura Stoltz.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) ARG, the argument.  If JINT = 1, then
    !    the argument must be less than XMAX.
    !
    !    Output, real ( kind = 8 ) RESULT, the value of the function,
    !    which depends on the input value of JINT:
    !    1, RESULT = I1(x);
    !    2, RESULT = exp(-x) * I1(x);
    !
    !    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
    !    1, I1(x);
    !    2, exp(-x) * I1(x);
    !---------------------------------------------------------------------
    Implicit None

    Integer(Kind=ipr) :: j,jint
    Real(Kind=pr) :: a,arg,b,exp40,forty
    Real(Kind=pr) :: one5,p(15),pbar,pp(8),q(5),qq(6),rec15,result,sump
    Real(Kind=pr) :: sumq,two25,x,xinf,xmax,xsmall,xx
    !  Mathematical constants
    Data one5 /15.0_pr/
    Data exp40 /2.353852668370199854d17/
    Data forty /40.0_pr/
    Data rec15 /6.6666666666666666666d-2/
    Data two25 /225.0_pr/
    !  Machine-dependent constants
    Data xsmall /5.55d-17/
    !Data xinf /1.79d308/
    Data xmax /713.987_pr/
    !  Coefficients for XSMALL <= ABS(ARG) < 15.0
    Data p/-1.9705291802535139930d-19,-6.5245515583151902910d-16, &
           -1.1928788903603238754d-12,-1.4831904935994647675d-09, &
           -1.3466829827635152875d-06,-9.1746443287817501309d-04, &
           -4.7207090827310162436d-01,-1.8225946631657315931d+02, &
           -5.1894091982308017540d+04,-1.0588550724769347106d+07, &
           -1.4828267606612366099d+09,-1.3357437682275493024d+11, &
           -6.9876779648010090070d+12,-1.7732037840791591320d+14, &
           -1.4577180278143463643d+15/
    Data q/-4.0076864679904189921d+03, 7.4810580356655069138d+06, &
           -8.0059518998619764991d+09, 4.8544714258273622913d+12, &
           -1.3218168307321442305d+15/
    !  Coefficients for 15.0 <= ABS(ARG)
    Data pp/-6.0437159056137600000d-02, 4.5748122901933459000d-01, &
            -4.2843766903304806403d-01, 9.7356000150886612134d-02, &
            -3.2457723974465568321d-03,-3.6395264712121795296d-04, &
             1.6258661867440836395d-05,-3.6347578404608223492d-07/
    Data qq/-3.8806586721556593450d+00, 3.2593714889036996297d+00, &
            -8.5017476463217924408d-01, 7.4212010813186530069d-02, &
            -2.2835624489492512649d-03, 3.7510433111922824643d-05/
    Data pbar/3.98437500d-01/

    x = Abs(arg)
    !
    !  Return for ABS(ARG) < XSMALL.
    !
    If(x<xsmall) Then

       result = half * x
    !
    !  XSMALL <= ABS(ARG) < 15.0.
    !
    Else If (x<one5) Then

       xx = x * x
       sump = p(1)
       Do j = 2, 15
          sump = sump * xx + p(j)
       End Do
       xx = xx - two25

       sumq = (((( xx + q(1) ) * xx + q(2) ) * xx + q(3) ) * xx + q(4) ) * xx + q(5)

       result = ( sump / sumq ) * x

       If(jint==2) Then
          result = result *Exp(-x)
       End If

    Else If (jint==1 .and. xmax<x) Then

       result = Huge(1.0_pr)

    Else
       !
       !  15.0 <= ABS(ARG).
       !
       xx = one / x - rec15

       sump = (((((( pp(1) * xx + pp(2) ) * xx + pp(3) ) * xx + pp(4) ) * xx + pp(5) ) &
                           * xx + pp(6) ) * xx + pp(7) ) * xx + pp(8)

       sumq = ((((( xx + qq(1) ) * xx + qq(2) ) * xx + qq(3) ) &
                  * xx + qq(4) ) * xx + qq(5) ) * xx + qq(6)

       result = sump / sumq

       If(jint/=1) Then
          result = (result + pbar) / Sqrt(x)
       Else
          !
          !  Calculation reformulated to avoid premature overflow.
          !
          If((xmax-one5)<x) Then
             a = Exp(x-forty)
             b = exp40
          Else
             a = Exp(x)
             b = one
          End If

          result = ( (result*a + pbar*a) / Sqrt(x) ) * b

       End If
    End If

    If(arg<zero) Then
       result = -result
    End If

    Return
  End Subroutine calci1

End Module bessik
!==================================================================================================================================
!#END bessik
!==================================================================================================================================
!================================================================================================================================
!#START hfodd_hfbtho MODULE
!================================================================================================================================
Module hfodd_hfbtho
  !
  Use HFBTHO_utilities
  Use linear_algebra
  !
  Implicit None
  !
 Contains
  Subroutine Main_Program(Nshell_THO, b0_THO, bz_THO, bp_THO, beta0_THO, b2_THO, b4_THO, neutron_THO, &
                          proton_THO, real_neutron_THO, real_proton_THO, type_THO, niter_THO,         &
                          accuracy_THO, restart_THO, functional_THO,                                  &
                          initial_pair_THO, coulomb_THO, VpairN_THO, VpairP_THO, Ecut_THO, V1_THO,    &
                          n_blocking_THO, p_blocking_THO, lambda_THO, active_THO, expect_THO,         &
                          useTHO_THO, gauge_THO, projection_THO, dN_THO, dZ_THO, set_temper_THO,      &
                          temper_THO, numGauss_THO, numLaguerre_THO, numLegendre_THO, HFODD_to_HFBTHO,&
                          nstate_THO, force_parity_THO, file_THO, output_THO)
    Use HFBTHO
    Implicit None
    Integer(ipr) :: iblocase(2),nkblocase(2,5)
    Integer(ipr) :: i,it,icount,l,noForce
    ! Inputs
    Logical      :: initial_pair_THO, set_temper_THO, force_parity_THO, HFODD_to_HFBTHO
    Integer(ipr) :: Nshell_THO, neutron_THO, proton_THO, type_THO, niter_THO, restart_THO, &
                    coulomb_THO, useTHO_THO, gauge_THO, projection_THO, dN_THO, dZ_THO,    &
                    output_THO, numGauss_THO, numLaguerre_THO, numLegendre_THO, nstate_THO
    Integer(ipr) :: p_blocking_THO(1:5), n_blocking_THO(1:5), lambda_THO(1:8), active_THO(1:8)
    Real(pr)     :: b0_THO, bz_THO, bp_THO, beta0_THO, accuracy_THO, VpairN_THO, VpairP_THO, &
                    temper_THO, b2_THO, b4_THO, V0n_THO, V0p_THO, Ecut_THO, V1_THO,          &
                    real_neutron_THO, real_proton_THO
    Real(pr)     :: expect_THO(1:8)
    Character(Len=30) :: functional_THO
    Character(Len=30) :: file_THO
    !
    !-------------------------------------------------------------------
    ! Read input and namelist data for the requested nucleus
    !-------------------------------------------------------------------
    Call initialize_HFBTHO_NAMELIST
    !
    !Call read_HFBTHO_NAMELIST
    !memo: Namelist /HFBTHO_GENERAL/ number_of_shells,oscillator_length,&
    !                                proton_number,neutron_number,type_of_calculation
    !      Namelist /HFBTHO_ITERATIONS/ number_iterations, accuracy
    !      Namelist /HFBTHO_FUNCTIONAL/ functional, add_initial_pairing, type_of_coulomb
    !      Namelist /HFBTHO_CONSTRAINTS/ lambda_values, lambda_active, expectation_values
    !      Namelist /HFBTHO_BLOCKING/ proton_blocking(1:5), neutron_blocking(1:5)
    !      Namelist /HFBTHO_BLOCKING/ switch_to_THO, projection_is_on, gauge_points, delta_N, delta_P
    !      Namelist /HFBTHO_TEMPERATURE/ set_temperature, temperature
    !      Namelist /HFBTHO_DEBUG/ number_Gauss, number_Laguerre, number_Legendre, &
    !                              force_parity, print_time
    welfile_INI               = file_THO             ! File containing HFBTHO binary results
    lout                      = output_THO           ! standard output redirected to file
    lfile                     = lout                 ! avoid creating spurious file hodef.dat
    !
    n00_INI                   = Nshell_THO           ! number of shells
    b0_INI                    = b0_THO               ! oscillator length (spherical)
    bz_INI                    = bz_THO               ! oscillator length (parallel)
    bp_INI                    = bp_THO               ! oscillator length (perpendicular)
    q_INI                     = beta0_THO            ! deformation beta_2 of the basis
    b2_0                      = b2_THO               ! initial deformation of the WS potential
    b4_0                      = b4_THO               ! initial deformation of the WS potential
    npr_INI(1)                = neutron_THO          ! N (as integer)
    npr_INI(2)                = proton_THO           ! Z (as integer)
    tz(1)                     = real_neutron_THO     ! N (as real number)
    tz(2)                     = real_proton_THO      ! Z (as real number)
    kindhfb_INI               = type_THO             ! 1: HFB, -1: HFB+LN
    !
    MAX_ITER_INI              = niter_THO            ! max number of iterations
    epsi_INI                  = accuracy_THO         ! convergence of iterations
    inin_INI                  = restart_THO          ! restart from file
    !
    skyrme_INI                = TRIM(functional_THO) ! functional
    Add_Pairing_INI           = initial_pair_THO     ! add pairing starting from file
    icou_INI                  = coulomb_THO          ! coul: no-(0), dir.only-(1), plus exchange-(2)
    !
    set_pairing               = .True.               ! pairing is defined by user if .True.
    V0n_INI                   = VpairN_THO           ! pairing strength for neutrons
    V0p_INI                   = VpairP_THO           ! pairing strength for protons
    pwi_INI                   = Ecut_THO             ! pairing q.p. cutoff
    cpv1_INI                  = V1_THO               ! Type of pairing: volume, surface, mixed
    !
    nkblocase(1,:)            = n_blocking_THO       ! config. of neutron blocked state
    nkblocase(2,:)            = p_blocking_THO       ! config. of proton blocked state
    !
    lambda_values             = lambda_THO           ! Values of lambda
    lambda_active             = active_THO           ! which lambda are active constraints
    expectation_values        = expect_THO           ! Requested values for the constraints
    !
    iLST_INI                  = useTHO_THO           ! 0:HO, -1:HO->THO, 1:THO
    keypj_INI                 = gauge_THO            ! PNP: number of gauge points
    iproj_INI                 = projection_THO       ! projecting on different nucleus
    npr1pj_INI                = dN_THO               ! its neutron number
    npr2pj_INI                = dZ_THO               ! its proton number
    !
    switch_on_temperature     = set_temper_THO       ! switches on temperature mode
    temper                    = temper_THO           ! value of the temperature
    !
    ngh_INI                   = numGauss_THO         ! number of Gauss-Hermite points for z-direction
    ngl_INI                   = numLaguerre_THO      ! number of Gauss-Laguerre points for rho-direction
    nleg_INI                  = numLegendre_THO      ! number of Gauss-Legendre points for Coulomb
    basis_HFODD_INI           = HFODD_to_HFBTHO      ! flag to enforce same basis as HFODD
    nstate_INI                = nstate_THO           ! total number of states in basis
    Parity_INI                = force_parity_THO     ! forces reflection symmetry
    IDEBUG_INI                = print_time           ! debug
    DO_FITT_INI               = .False.              ! calculates quantities for reg.optimization
    Print_HFBTHO_Namelist_INI = .True.               ! Print Namelist
    !
    If(SUM(lambda_active).Gt.0) Then
       icount=0
       Do l=1,lambdaMax,2
          If(lambda_active(l).Gt.0) icount=icount+1
       End Do
       If(icount.Gt.0) Parity_INI=.False.
    End If
    !
    Call read_UNEDF_NAMELIST(skyrme_INI,noForce)
    ! If functional is used, projection automaticaly switched off
    If(noForce.Eq.0) iproj_INI=0
    !---------------------------------------------------------------------------
    ! GROUND STATE BLOCKING WALKER: blocking candidates are predefined
    ! by the parent nucleus and we block them one by one
    !---------------------------------------------------------------------------
    iblocase=0; bloqpdif=zero !blomax will be charged from the previous solution
    Do it=1,2
       If(nkblocase(it,1).Ne.0.And.nkblocase(it,2).Eq.0) Then
          If(it.Eq.1) Then
             iblocase(1)=iblocase(1)+1
             If(iblocase(1).Gt.blomax(1)) iblocase(1)=1
          Else
             If(iblocase(1).Le.1) iblocase(2)=iblocase(2)+1
          End If
          nkblo_INI(it,1)=Sign(iblocase(it),nkblocase(it,1))
          nkblo_INI(it,2)=0
       Else
          ! case of external blocking
          nkblo_INI(it,:)=nkblocase(it,:)
       Endif
    End Do
    !---------------------------------------------------------------------------
    ! MANUAL BLOCKING: manualBlocking=1 in the module (ocasionaly used)
    ! One types which level to be blocked referencing the parent nucleus
    !---------------------------------------------------------------------------
    !If(manualBlocking.Ne.0) Then
    !   Write(*,'(a,5(1pg12.4))') 'Please print the number of the neutron level to block, num='
    !   Read(*,*) nkblo_INI(1,1)
    !   Write(*,'(a,5(1pg12.4))') 'Neutron blocked level num=',nkblo_INI(1,1)
    !   nkblo_INI(1,2)=0
    !   Write(*,'(a,5(1pg12.4))') 'Please print the number of the proton level to block, num='
    !   Read(*,*) nkblo_INI(2,1)
    !   Write(*,'(a,5(1pg12.4))') 'Proton blocked level num=',nkblo_INI(2,1)
    !   nkblo_INI(2,2)=0
    !End If
    !--------------------------------------------------------------------
    ! Calculations for 'FITS' functional (Modifies some values if needed)
    !--------------------------------------------------------------------
    If (Trim(skyrme_INI).Eq.'FITS') Then
       !
       DMEORDER=-1; DMELDA=0; use_cm_cor=.False.; use_TMR_pairing=0
       !
       HBZERO =  20.7355300000000_pr; E2CHARG =   1.4399784085965_pr
       CRHO(0)=-731.2227858295098_pr; CDRHO(0)= 855.6900515849785_pr; CTAU(0) = -0.5439888609059821_pr
       CRHO(1)= 263.7103055246761_pr; CDRHO(1)=-176.8641956040411_pr; CTAU(1) =-33.3618818665213400_pr
       CRDR(0)= -43.2900897553154_pr; CJ(0)   =   0.0000000000000_pr; CRDJ(0) =-75.2608700482894100_pr
       CRDR(1)=-164.1379857135440_pr; CJ(1)   =   0.0000000000000_pr; CRDJ(1) =-22.6528199648713000_pr
       CPV0(0)=-186.1922465962490_pr; CPV1(0) =   0.5000000000000_pr; SIGMA   =  0.2987839827782357_pr
       CPV0(1)=-206.7464168983860_pr; CPV1(1) =   0.5000000000000_pr; CEXPAR  =  0.6391295237623640_pr
       !
       RHO_NM=0.15732716296394680_pr;
       E_NM=-15.80000048487058000_pr; K_NM    =225.94339488338960_pr; SMASS_NM=  0.9958725808228520_pr
       ASS_NM=28.3483385569865900_pr; LASS_NM =40.001962979089330_pr; VMASS_NM=  1.2489999532699580_pr
       !
    End If
    !--------------------------------------------------------------------
    ! Run the solver in all cases EVEN/ODDS, FITS/NO-FITS
    !--------------------------------------------------------------------
    Call HFBTHO_SOLVER
    !--------------------------------------------------------------------
    ! Display error messages in case of problems
    !--------------------------------------------------------------------
    If (ierror_flag.Ne.0) Then
        Write(lout,*)
        Write(lout,'(a)') ' ERRORS IN HFBTHO_SOLVER'
        Do i=1,ierror_flag
           Write(lout,'(a,i2,2x,a)') ' error_flag=',i,ierror_info(i)
        End Do
        Write(lout,*)
    Else
        Write(lout,*)
        Write(lout,'(a)') ' HFBTHO_SOLVER ended without errors'
        Write(lout,*)
    End If
    !
    If (lout.Lt.lfile) Close(lfile)                 ! close the output
    !
    !help
    !=============================================================================================
    ! Output variables nucname,ereslbl(1:2),eres(1:ierest) contain:
    !  LBL,BLKN,BLKZ,Jsi,JININ,A,N',Z,Efn,Efp,JEtot,Jbett,Jbetn,Jbetp,JQt,JQn,JQp
    !  JpEn,JpEp,JpDn,JpDp,JAsn,JAsp,Jrt,Jrn,Jrp,Jrc,Jht,Jhn,Jhp,Jqht,Jqhn,Jqhp,
    !  JKINt,JKINn,JKINp,JSO,JCDIR,JCEX,JDisn,JDisp,JV2Mn,JV2Mp,JILST,JKIND,JL,
    !  JECMPAV1,JECMPAV2,JECMPAV3,JA,JN,JZ,ITER,UEtot,Ubett,Ubetn,Ubetp,UQt,UQn,
    !  UQp,Uln,Ulp,UpEn,UpEp,UpDn,UpDp,UAsn,UAsp,Urt,Urn,Urp,Urc,Uht,Uhn,Uhp,
    !  Uqht,Uqhn,Uqhp,UKINT,UKINN,UKINP,USO,UCDIR,UCEX,UDisn,UDisp,UV2Mn,UV2Mp,
    !  UECMT,UECMN,UECMP,UROTT,UROTN,UROTP,USQUJT,USQUJN,USQUJP,UCRANT,UCRANN,
    !  UCRANP,UERIGT,UERIGN,UERIGP,EHFBLN,EHFB,LNbet,LNben,LNbep,LNQt,LNQn,LNQp,
    !  LNpEn,LNpEp,LNpDn,LNpDp,LNrt,LNrn,LNrp,LNrC,LNam2n,LNam2p,LNe2n,LNe2p,
    !  BlEqpN,BlDEqpN,BlOvrN,BlEqpZ,BlDEqpZ,BlOvrZ
    !=============================================================================================
    !
  End Subroutine Main_Program
  !==================================================================================================================================
  !#END MAINPROGRAM
  !==================================================================================================================================
  !#START HFBTHO_SOLVER
  !==================================================================================================================================
  Subroutine HFBTHO_SOLVER
    !-----------------------------------------------------------------------------------------------
    ! Universal HFBTHO_SOLVER
    !
    !   Axially-deformed configurational constrained and/or unconstrained Hartree-Fock-Bogoliubov
    !   calculations with Skyrme-like functionals and delta pairing using the Harmonic-Oscillator
    !   (HO), and/or Transformed HO (THO) basis with or without reflection symmetry imposed, with
    !   or without the Lipkin-Nogami procedure. The solver can handle all Skyrme-like functionals,
    !   DME-functionals, Fayans-functionals, calculate infinite nuclear matter properties, finite
    !   nuclei (even-even, odd-even, odd-odd), and neutron drops, isoscalar and isovector monopo-
    !   le FAM QRPA calculations for spherical and deformed nuclei.
    !
    !   All necessary input variables contain the suffix _INI. Below, the complete list of these
    !   variables with some example values:
    !
    !   ======== hfbtho_NAMELIST.dat
    !   n00_INI=20; npr1_INI=70; npr2_INI=50;  kindhfb_INI=-1; inin_INI=-1
    !   b0_INI=2.234776; q_INI=0.0; cdef_INI=0.0; cqad_INI=0.5; skyrme_INI='SLY4'; nkblo_INI=0
    !   ILST_INI=0; keypj_INI=1; iproj_INI=0; npr1pj_INI=0;
    !   icou_INI=2; IDEBUG_INI=0; npr2pj_INI=0;
    !   Parity_INI=.False.; epsi_INI=0.00001_pr; MAX_ITER_INI=101
    !   Add_Pairing_INI=.False.; DO_FITT_INI=.False.; Print_PTHO_Namelist_INI=.True.
    !
    !   ======== from read_UNEDF_NAMELIST
    !   DMEORDER=-1; DMELDA=0; use_TMR_pairing=0
    !   HBZERO=20.73553000000000;    E2CHARG=1.4399784085965135; CRHO(0)=-933.3423749999999;  CRHO(1)=830.0524855000001;
    !   CDRHO(0)=861.0625000000000;  CDRHO(1)=-1064.2732500000;  CTAU(0)=57.12868750000000;   CTAU(1)=24.65673650000000;
    !   CRDR(0)=-76.99620312499999;  CRDR(1)=15.65713512500000;  CRDJ(0)=-92.25000000000000;  CRDJ(1)=-30.7500000000000;
    !   CJ(0)=17.20961150000000;     CJ(1)=64.57581250000000;    CPV0(0)=-258.2000000000000;  CPV0(1)=-258.2000000000000;
    !   CPV1(0)=0.5000000000000000;  CPV1(1)=0.500000000000000;  SIGMA=0.1666666666666667;    CEXPAR=1.000000000000000;
    !   E_NM=-15.97214914144462;     K_NM=229.9009644826037;     SMASS_NM =1.439546988976078; RHO_NM =0.1595387567117334;
    !   ASS_NM =32.00430281505202;   LASS_NM=45.96175148046161;  VMASS_NM =1.249838547196253;
    !   MPI=0.6995945261023822;      GA=1.290000000000000;       FPI=0.4683223517486062;      C1=-0.1598130000000000;
    !   C3 =-0.6708200000000;        C4 =0.6708200000000000;     CD =-2.062000000000000;      CE=-0.6250000000000;
    !   LAMBDAX =3.547896604156107;  USE_INM=.false.;            USE_CM_COR =.true.;          USE_DME3N_TERMS=.true.;
    !   USE_J2TERMS =.true.;         USE_CHARGE_DENSITY=.false.; PRINT_NAMELIST=.true.;
    !
    ! Memo:
    !  -  inin_INI switches scratch unconstrained (inin=1,2,3) or constrained (inin=100,200,300)
    !     calculations. Unconstrained mode begins with a small number of constrained iterations.
    !  -  inin_INI switches unconstrained (inin=-1,-2,-3) or constrained (inin=-100,-200,-300)
    !     calculations from a previous solution if the latter exists. If not, the solver sets
    !     inin=Abs(inin) and resumes from scratch.
    !  -  The same holds for odd nuclei. If even-even solution for the odd nucleus does not exists
    !     it is calculated first.
    !  -  Print_Screen=T/F for n00_INI=+/-: output is generated and written in thoout.dat file
    !	only if n00_INI>0. For n00_INI<0, the number of shells is set to abs(n00_INI) but all
    !     output is supressed.
    !
    !  -  At the end of the solution, the solver provides all results in the arrays
    !
    !                  nucname,ereslbl(1:2),eres(1:ierest)
    !
    !     which contain:
    !
    !        LBL,BLKN,BLKZ,Jsi,JININ,A,N',Z,Efn,Efp,JEtot,Jbett,Jbetn,Jbetp,JQt,JQn,JQp
    !        JpEn,JpEp,JpDn,JpDp,JAsn,JAsp,Jrt,Jrn,Jrp,Jrc,Jht,Jhn,Jhp,Jqht,Jqhn,Jqhp,
    !        JKINt,JKINn,JKINp,JSO,JCDIR,JCEX,JDisn,JDisp,JV2Mn,JV2Mp,JILST,JKIND,JL,
    !        JECMPAV1,JECMPAV2,JECMPAV3,JA,JN,JZ,ITER,UEtot,Ubett,Ubetn,Ubetp,UQt,UQn,
    !        UQp,Uln,Ulp,UpEn,UpEp,UpDn,UpDp,UAsn,UAsp,Urt,Urn,Urp,Urc,Uht,Uhn,Uhp,
    !        Uqht,Uqhn,Uqhp,UKINT,UKINN,UKINP,USO,UCDIR,UCEX,UDisn,UDisp,UV2Mn,UV2Mp,
    !        UECMT,UECMN,UECMP,UROTT,UROTN,UROTP,USQUJT,USQUJN,USQUJP,UCRANT,UCRANN,
    !        UCRANP,UERIGT,UERIGN,UERIGP,EHFBLN,EHFB,LNbet,LNben,LNbep,LNQt,LNQn,LNQp,
    !        LNpEn,LNpEp,LNpDn,LNpDp,LNrt,LNrn,LNrp,LNrC,LNam2n,LNam2p,LNe2n,LNe2p,
    !        BlEqpN,BlDEqpN,BlOvrN,BlEqpZ,BlDEqpZ,BlOvrZ
    !
    !  -  Standard inputs from               'hfbtho_NAMELIST.dat'
    !  -  USer-defined functional read from  'UNEDF_NAMELIST.DAT'
    !  -  Final solution stored to files     '*.hel' and/or '*.tel'
    !  -  Output is written to files         'thoout.dat', 'thodef.dat' and 'hodef.dat'
    !  -  Output files *.dat may exist as    'thoout','thores','hodenp',
    !                                      'thodenp','thodef','thoene',
    !                                      'thoprc','dat0.1.2.3.4'
    !  -  External accuracy pr/ipr always set in UNEDF module
    !
    !  -  n00    Number of oscillator shells
    !             n00>0 prints to thoout.dat & screen
    !             n00<0 no print at all
    !             n00=0 program stops (NB!)
    !  -  b0     Oscillator Basis parameter b0>0 (If b0<0 it takes a default value)
    !  -  beta0  Value of Basis deformation parameter
    !  -  AN     Number of neutrons N
    !  -  AZ     Number of protons Z
    !  -  FTST'  Fayance forces label
    !  -  kind   Kind of calculations 1: noLN, -1:LN
    !  -  inin   Unconstraint Iterations from peviouse solution
    !             -1: (from spherical *.hel or *.tel file)
    !             -2: (from prolate   *.hel or *.tel file)
    !             -3: (from oblate    *.hel or *.tel file)
    !            Unconstraint Iterations fom scratch with
    !            preliminary constraint at deformation Cbeta, (i):
    !              1: Spherical scratch
    !              2: Prolate scratch
    !              3: Oblate scratch
    !            Constrained calculation (icstr) at Cbeta, see (i)
    !              100, 200, 300 fom scratch
    !             -100,-200,-300 fom previouse solution
    !  -  blNeutrons: a group responsible for blocking a particular neutron level
    !     The group consists of 5 numbers, e.g., for 7-[ 3, 0, 3]: 7 -1 3 0 3
    !      k1  2 \times \Omega
    !         =0: the whole group (k) is disregarded (n0 blocking)
    !         >0: blocking in N+1 nucleus
    !         <0: blocking in N-1 nucleus
    !      k2  parity (+1 or -1); NB! when k2=0, the ground state walker is applied
    !      k3,k4,k5  Nilson quantum numbers
    !  -  blProtons: exactly the same as (k) but for protons
    !
    !-----------------------------------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Use HFBTHO_THO
    Implicit None
    Integer(ipr) :: iw,ib,j,i,it,l,maxi0,icstr0,iterMax,icons,il,kickoff
    Real(pr)     :: epsi0,qq,f,f1,f2,f3,r,g,g1
    !-------------------------------------------------------------
    ! Checking consistency of input values
    !-------------------------------------------------------------
    Call check_consistency
    If(ierror_flag.Ne.0) Return
    !-------------------------------------------------------------
    ! Initializing all according to *_INI values
    !-------------------------------------------------------------
    Call initialize_HFBTHO_SOLVER
    If(ierror_flag.Ne.0) Return
    If(lout.Lt.lfile) Open(lfile,file='thoout.dat',status='unknown')
    Call Constraint_or_not(inin_INI,inin,icstr)
    If(ierror_flag.Ne.0) Return
    !-------------------------------------------------------------------------
    ! Loop recalculating eventually the even-even solution for an odd nucleus
    !-------------------------------------------------------------------------
    irestart=0
    Do
       n00=Abs(n00_INI);  b0=b0_INI;           q=q_INI;           iLST=iLST_INI;
       maxi=MAX_ITER_INI; npr(1)=npr_INI(1);   npr(2)=npr_INI(2); npr(3)=npr(1)+npr(2);
       skyrme=skyrme_INI; kindhfb=kindhfb_INI
       keypj=keypj_INI;   iproj=iproj_INI;     npr1pj=npr1pj_INI; npr2pj=npr2pj_INI;
       nkblo=nkblo_INI;
       basis_HFODD=basis_HFODD_INI
       !-------------------------------------------------------------
       ! Define the set of constraints
       !-------------------------------------------------------------
       numberCons=0; kickoff=0
       Do l=1,lambdaMax
        If(Abs(lambda_active(l)).Gt.0) numberCons = numberCons + 1
        If(lambda_active(l).Lt.0) kickoff = kickoff + 1
       End Do
       !
       If(.Not.Allocated(multLag)) Allocate(multLag(1:lambdaMax)); multLag=zero
       If(.Not.Allocated(multLambda)) Allocate(multLambda(1:numberCons)); multLambda=0
       If(.Not.Allocated(multRequested)) Allocate(multRequested(0:lambdaMax)); multRequested=zero
       !
       icons=0
       Do l=1,lambdaMax
          If(Abs(lambda_active(l)).Gt.0) Then
              icons=icons+1
              multLambda(icons)=lambda_values(l)
          End If
          multRequested(l) = expectation_values(l)
       End Do
       !-------------------------------------------------------------
       ! Blocking
       !-------------------------------------------------------------
       Do it=1,2
          If(nkblo(it,1).Ne.0) Then
               If(nkblo(it,1).Gt.0) Then
                    ! particle state
                 npr(it)=npr(it)+1
                 iparenti(it)=-1
               Else
                    ! hole state
                 npr(it)=npr(it)-1
                 iparenti(it)=+1
               End If
               nkblo(it,1)=Abs(nkblo(it,1))
               If(nkblo(it,2).Eq.0) Then
                 ! ground state walker
                 keyblo(it)=keyblo(it)+1 !nkblo(it,1)
                 If(keyblo(it).Eq.blomax(it)) irestart=0
               Else
                 irestart=0
               End If
          End If
       End Do
       !-------------------------------------------------------------
       ! HFB+HO calculations
       !-------------------------------------------------------------
       If(ILST.Le.0) Then
          icacou=0; icahartree=0
          Call preparer(.True.)
          If(ierror_flag.Ne.0) Return
          Call inout(1)
          If(ierror_flag.Ne.0) Return
          !-------------------------------------------------------------
          ! Preliminary constrained calculations
          !-------------------------------------------------------------
          If(kickoff.Gt.0.And.icstr.Eq.0) Then
             icstr0=icstr; epsi0=epsi; ! remember accuracy
             icstr=1                   ! constraint true
             epsi=1.0_pr               ! small accuracy
             iterMax = maxi; maxi = 10
             numberCons=0
             Do l=1,lambdaMax
                If(Abs(multRequested(l)).Gt.1.D-14) Then
                   numberCons=numberCons+1
                   multLambda(numberCons)=lambda_values(l)
                End If
             End Do
             Do iw=lout,lfile
                If(Parity) Then
                   Write(iw,'(/,a,i3,a,i2,a,/)') '  ### INITIAL STAGE(constrained calculations, reflection symmetry used)'
                Else
                   Write(iw,'(/,a,i3,a,i2,a,/)') '  ### INITIAL STAGE(constrained calculations, no reflection symmetry used)'
                End If
             End Do
             Call iter(.True.)     ! small constraint iterations
             If(ierror_flag.Ne.0) Return
             icstr=icstr0; epsi=epsi0
             maxi = iterMax
             numberCons=0
          End If
          !-------------------------------------------------------------
          ! REGULAR HFB+HO ITERATIONS
          !-------------------------------------------------------------
          Do iw=lout,lfile
             If(Parity) Then
                Write(iw,'(/,a,i3,a,i2,a,/)')    '  ### REGULAR STAGE (reflection symmetry imposed)'
             Else
                Write(iw,'(/,a,i3,a,i2,a,/)')    '  ### REGULAR STAGE (no reflection symmetry imposed)'
             End If
          End Do
          Call iter(.True.)
          If(ierror_flag.Ne.0) Return
          Call resu(1)
          If(ierror_flag.Ne.0) Return
       End If
       !! write LST function on disk
       !Open(unit=66,file='LST.dat',status='unknown')
       !Write(66,'("#",15X,"R",20X,"f(R)",20X,"f^(1)",20X,"f^(2)",20X,"f^(3)")')
       !Do il=1,170
       !   qq=Real(il-1)/10.0_pr*1.0_pr
       !   If(il.Eq.1) Call thofun(0,g,f,f1,f2,f3,g1,.True.,.True.)
       !   Call thofun(1,qq,f,f1,f2,f3,r,.False.,.True.)
       !   Write(66,'(6E24.10)') r,qq,f1,f2,f3
       !End Do
       !Close(66)
       !-------------------------------------------------------------
       ! HFB+THO calculations from HFB+HO
       !-------------------------------------------------------------
       If(ILST.Lt.0) Then
          ILST1=1; icacou=0; icahartree=0
          Call coordinateLST(.False.) ! THO basis
          If(ierror_flag.Ne.0) Return
          Call densit                 ! THO densities
          If(ierror_flag.Ne.0) Return
          Call field                  ! Nuclear fields
          If(ierror_flag.Ne.0) Return
          Call iter(.True.)           ! HFB+THO iterations
          If(ierror_flag.Ne.0) Return
          Call resu(1)                ! print/record results
          If(ierror_flag.Ne.0) Return
       End If
       !-------------------------------------------------------------
       ! HFB+THO calculations from *.tel
       !-------------------------------------------------------------
       If(ILST.Gt.0) Then
          If(inin.Gt.0) Then
             ierror_flag=ierror_flag+1
             ierror_info(ierror_flag)= ' Forbidden iLST>0, inin>0 '
             Return
          End If
          icacou=0; icahartree=0
          Call preparer(.True.)
          If(ierror_flag.Ne.0) Return
          Call inout(1)               ! reading HFB matrices
          If(ierror_flag.Ne.0) Return
          Call iter(.True.)           ! HFB+THO iterations
          If(ierror_flag.Ne.0) Return
          Call resu(1)                ! print/record results
          If(ierror_flag.Ne.0) Return
       End If
       !-------------------------------------------------------------
       ! Go for the requested blocking state in a case of odd nuclei
       ! if restarted due to corrupted/missing previous solution
       !-------------------------------------------------------------
       inin=-Abs(inin)
       If(irestart.Eq.0) Exit
    End Do
    !
  End Subroutine HFBTHO_SOLVER
  !=======================================================================
  !
  !=======================================================================
  Subroutine heading
    !---------------------------------------------------------------------
    ! print heading to screen 'lout' and to tape thoout.dat 'lfile'
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
#if(USE_OPENMP==1)
    Use omp_lib
#endif
    Implicit None
    Integer(ipr) :: iw,idt(8),numThreads,idThread
    Character(len=12) rcl(3)
    Character(len=50) today
    !
#if(USE_OPENMP==1)
!$OMP PARALLEL PRIVATE(numThreads,idThread)
    numThreads = omp_get_num_threads()
    idThread = omp_get_thread_num()
    If (idThread .Eq. 0) Then
        Write(lout,'("Multi-threading framework with OpenMP:",i2," threads/task")') numThreads
    End If
!$OMP END PARALLEL
#endif
    Call Date_and_time(rcl(1),rcl(2),rcl(3),idt)
    Write(today,'(a,i2,a,i2,a,i4,a,i2,a,i2,a)')'(',idt(2),'/',idt(3),'/',idt(1),', ',idt(5),':',idt(6),')'
    Do iw=lout,lfile
       Write(iw,'(a)')
       Write(iw,'(a)')      '  ======================================='
       Write(iw,'(a,i2,a)') '           FORTRAN 95 CODE (KIND=',pr,') '
       Write(iw,'(a,a)')    '               Version: ',Version
       Write(iw,'(a)')      '  ======================================='
       Write(iw,'(a)')      '       AXIALLY DEFORMED CONFIGURATIONAL  '
       Write(iw,'(a)')      '     HARTREE-FOCK-BOGOLIUBOV CALCULATIONS'
       Write(iw,'(a)')      '                     WITH                '
       Write(iw,'(a)')      '            UNEDF AND DELTA PAIRING      '
       Write(iw,'(a)')      '                     USING               '
       Write(iw,'(a)')      '             HARMONIC-OSCILLATOR         '
       Write(iw,'(a)')      '                   AND/OR                '
       Write(iw,'(a)')      '       TRANSFORMED HARMONIC-OSCILLATOR   '
       Write(iw,'(a)')      '                    BASIS                '
       Write(iw,'(a)')      '                     ---                 '
       Write(iw,'(a)')      '               v1.66  (2005):,           '
       Write(iw,'(a)')      '   Stoitsov,Dobaczewski,Nazarewicz,Ring, '
       Write(iw,'(a)')      '               v2.00d (2012):,           '
       Write(iw,'(a)')      '         Stoitsov,Schunck,Kortelainen    '
       Write(iw,'(a)')      '  ======================================='
       Write(iw,'(a,a,a,i4,a,i3,a,i3,a)')'    Nucleus: ',nucname,' (A=',npr(1)+npr(2),&
                                                          ', N=',npr(1),', Z=',npr(2),')'
       If(Parity) Then
          Write(iw,'(a)')     '       Reflection Symmetry Imposed       '
       Else
          Write(iw,'(a)')     '      No Reflection Symmetry Imposed     '
       End If
       Write(iw,'(a,a)')    '            ',today
       Write(iw,'(a)')      '  ======================================='
       Write(iw,'(a)')
       Write(iw,'(a)')
    End Do
  End Subroutine heading
  !=======================================================================
  !
  !=======================================================================
  Subroutine thodefh(iw1)
    !---------------------------------------------------------------------
    ! print labels to hodef.dat or/and thodef.dat files
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: iw1
    hlabels(1)='LBL';       hlabels(11)='JEtot';    hlabels(21)='JpDp';
    hlabels(2)='BLKN';      hlabels(12)='Jbett';    hlabels(22)='JAsn';
    hlabels(3)='BLKZ';      hlabels(13)='Jbetn';    hlabels(23)='JAsp';
    hlabels(4)='Jsi';       hlabels(14)='Jbetp';    hlabels(24)='Jrt';
    hlabels(5)='JININ';     hlabels(15)='JQt';      hlabels(25)='Jrn';
    hlabels(6)='A';         hlabels(16)='JQn';      hlabels(26)='Jrp';
    hlabels(7)='N';         hlabels(17)='JQp';      hlabels(27)='Jrc';
    hlabels(8)='Z';         hlabels(18)='JpEn';     hlabels(28)='Jht';
    hlabels(9)='Efn';       hlabels(19)='JpEp';     hlabels(29)='Jhn';
    hlabels(10)='Efp';      hlabels(20)='JpDn';     hlabels(30)='Jhp';
    !
    hlabels(31)='Jqht';    hlabels(41)='JDisp';    hlabels(51)='JN';
    hlabels(32)='Jqhn';    hlabels(42)='JV2Mn';    hlabels(52)='JZ';
    hlabels(33)='Jqhp';    hlabels(43)='JV2Mp';    hlabels(53)='ITER';
    hlabels(34)='JKINt';   hlabels(44)='JILST';    hlabels(54)='UEtot';
    hlabels(35)='JKINn';   hlabels(45)='JKIND';    hlabels(55)='Ubett';
    hlabels(36)='JKINp';   hlabels(46)='JL';       hlabels(56)='Ubetn';
    hlabels(37)='JSO';     hlabels(47)='JECMPAV1'; hlabels(57)='Ubetp';
    hlabels(38)='JCDIR';   hlabels(48)='JECMPAV2'; hlabels(58)='UQt';
    hlabels(39)='JCEX';    hlabels(49)='JECMPAV3'; hlabels(59)='UQn';
    hlabels(40)='JDisn';   hlabels(50)='JA';       hlabels(60)='UQp';
    !
    hlabels(61)='Uln';     hlabels(71)='Urp';      hlabels(81)='UKINP';
    hlabels(62)='Ulp';     hlabels(72)='Urc';      hlabels(82)='USO';
    hlabels(63)='UpEn';    hlabels(73)='Uht';      hlabels(83)='UCDIR';
    hlabels(64)='UpEp';    hlabels(74)='Uhn';      hlabels(84)='UCEX';
    hlabels(65)='UpDn';    hlabels(75)='Uhp';      hlabels(85)='UDisn';
    hlabels(66)='UpDp';    hlabels(76)='Uqht';     hlabels(86)='UDisp';
    hlabels(67)='UAsn';    hlabels(77)='Uqhn';     hlabels(87)='UV2Mn';
    hlabels(68)='UAsp';    hlabels(78)='Uqhp';     hlabels(88)='UV2Mp';
    hlabels(69)='Urt';     hlabels(79)='UKINT';    hlabels(89)='UECMT';
    hlabels(70)='Urn';     hlabels(80)='UKINN';    hlabels(90)='UECMN';
    !
    hlabels(91)='UECMP';   hlabels(101)='UERIGT';  hlabels(111)='LNQp';
    hlabels(92)='UROTT';   hlabels(102)='UERIGN';  hlabels(112)='LNpEn';
    hlabels(93)='UROTN';   hlabels(103)='UERIGP';  hlabels(113)='LNpEp';
    hlabels(94)='UROTP';   hlabels(104)='EHFBLN';  hlabels(114)='LNpDn';
    hlabels(95)='USQUJT';  hlabels(105)='EHFB';    hlabels(115)='LNpDp';
    hlabels(96)='USQUJN';  hlabels(106)='LNbet';   hlabels(116)='LNrt';
    hlabels(97)='USQUJP';  hlabels(107)='LNben';   hlabels(117)='LNrn';
    hlabels(98)='UCRANT';  hlabels(108)='LNbep';   hlabels(118)='LNrp';
    hlabels(99)='UCRANN';  hlabels(109)='LNQt';    hlabels(119)='LNrC';
    hlabels(100)='UCRANP'; hlabels(110)='LNQn';    hlabels(120)='LNam2n';
    !
    hlabels(121)='LNam2p';
    hlabels(122)='LNe2n';
    hlabels(123)='LNe2p';
    hlabels(124)='BlEqpN';
    hlabels(125)='BlDEqpN';
    hlabels(126)='BlOvrN';
    hlabels(127)='BlEqpZ';
    hlabels(128)='BlDEqpZ';
    hlabels(129)='BlOvrZ';
    !
    Write(iw1,'((1x,a,2x),6x,660(a,2x))') hlabels
    !
    ! HELP
    !Do i=1,129
    ! Write(iw1,'(1x,i3,a,a)',advance='NO') i,':',trim(hlabels(i))
    !End Do
    ! 1:LBL  2:BLKN  3:BLKZ  4:Jsi  5:JININ  6:A  7:N  8:Z  9:Efn 10:Efp
    ! 11:JEtot 12:Jbett 13:Jbetn 14:Jbetp 15:JQt 16:JQn 17:JQp 18:JpEn 19:JpEp 20:JpDn
    ! 21:JpDp 22:JAsn 23:JAsp 24:Jrt 25:Jrn 26:Jrp 27:Jrc 28:Jht 29:Jhn 30:Jhp
    ! 31:Jqht 32:Jqhn 33:Jqhp 34:JKINt 35:JKINn 36:JKINp 37:JSO 38:JCDIR 39:JCEX 40:JDisn
    ! 41:JDisp 42:JV2Mn 43:JV2Mp 44:JILST 45:JKIND 46:JL 47:JECMPAV1 48:JECMPAV2 49:JECMPAV3 50:JA
    ! 51:JN 52:JZ 53:ITER 54:UEtot 55:Ubett 56:Ubetn 57:Ubetp 58:UQt 59:UQn 60:UQp
    ! 61:Uln 62:Ulp 63:UpEn 64:UpEp 65:UpDn 66:UpDp 67:UAsn 68:UAsp 69:Urt 70:Urn
    ! 71:Urp 72:Urc 73:Uht 74:Uhn 75:Uhp 76:Uqht 77:Uqhn 78:Uqhp 79:UKINT 80:UKINN
    ! 81:UKINP 82:USO 83:UCDIR 84:UCEX 85:UDisn 86:UDisp 87:UV2Mn 88:UV2Mp 89:UECMT 90:UECMN
    ! 91:UECMP 92:UROTT 93:UROTN 94:UROTP 95:USQUJT 96:USQUJN 97:USQUJP 98:UCRANT 99:UCRANN 100:UCRANP
    ! 101:UERIGT 102:UERIGN 103:UERIGP 104:EHFBLN 105:EHFB 106:LNbet 107:LNben 108:LNbep 109:LNQt 110:LNQn
    ! 111:LNQp 112:LNpEn 113:LNpEp 114:LNpDn 115:LNpDp 116:LNrt 117:LNrn 118:LNrp 119:LNrC 120:LNam2n
    ! 121:LNam2p 122:LNe2n 123:LNe2p 124:BlEqpN 125:BlDEqpN 126:BlOvrN 127:BlEqpZ 128:BlDEqpZ 129:BlOvrZ
  End Subroutine thodefh
  !=======================================================================
  !
  !=======================================================================
  Subroutine thoalloc
    !---------------------------------------------------------------------
    ! Allocates arrays at given number of oscillator shells 'n00'
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer :: ier,ib,ND
    !
    ! number of int.points
    If(Parity) Then
       ngh=ngh_INI; ngl=ngl_INI; nleg=nleg_INI     !Yesp
    Else
       ngh=2*ngh_INI; ngl=ngl_INI; nleg=nleg_INI   !Nop
    End If
    !
    !nbx=2*n00+1                   ! maximal number of k-blocks
    !ntx=(n00+1)*(n00+2)*(n00+3)/6 ! max.num. p/n levels
    !nzx=n00                       ! maximal nz-quantum number
    !nrx=n00/2+1                   ! maximal nr-quantum number
    !nlx=n00                       ! maximal ml-quantum number
    !ndx=(n00+2)*(n00+2)/4         ! maximal dim. of one k-block
    !nhhdim=number of nonzero HH matrix elements
    !
    nzrlx=(nzx+1)*(nrx+1)*(nlx+1)   ! phy(:,:,nzrlx)
    nghl=ngh*ngl                    ! nghl=ngh*ngl
    nqx=ndx*ndx; nb2x=nbx+nbx; ndx2=ndx+ndx
    ilnqx=ilpj*nqx; ilnghl=ilpj*nghl
    nhfbx=ndx+ndx; nhfbqx=nhfbx*nhfbx; nkx=ntx; ndxs=ndx*(ndx+1)/2
    !-----------------------------------------
    !Arrays depending on gauss points
    !-----------------------------------------
    If(Allocated(xleg)) Deallocate(xleg,wleg)
    If(nleg.Gt.0) Allocate(xleg(nleg),wleg(nleg))
    If(Allocated(xh)) Deallocate(xh,wh,xl,sxl,wl,vc &
         ,vhbn,vn,vrn,vzn,vdn,vsn,dvn,vhbp,vp,vrp,vzp,vdp,vsp,dvp  &
         ,vSZFIn,vSFIZn,vSRFIn,vSFIRn,vSZFIp,vSFIZp,vSRFIp,vSFIRp &
         ,fl,fli,fh,fd,fp1,fp2,fp3,fp4,fp5,fp6  &
         ,fs1,fs2,fs3,fs4,fs5,fs6,wdcor,wdcori,cou,vDHartree,vhart00,vhart01,vhart11)
    Allocate(xh(ngh),wh(ngh),xl(ngl),sxl(ngl),wl(ngl),vc(nghl,nghl))
    Allocate(vhbn(nghl),vn(nghl),vrn(nghl),vzn(nghl),vdn(nghl),vsn(nghl),dvn(nghl)  &
         ,vhbp(nghl),vp(nghl),vrp(nghl),vzp(nghl),vdp(nghl),vsp(nghl),dvp(nghl)  &
         ,vSZFIn(nghl),vSFIZn(nghl),vSRFIn(nghl),vSFIRn(nghl)  &
         ,vSZFIp(nghl),vSFIZp(nghl),vSRFIp(nghl),vSFIRp(nghl))
    Allocate(fl(nghl),fli(nghl),fh(nghl),fd(nghl),fp1(nghl),fp2(nghl),fp3(nghl)  &
         ,fp4(nghl),fp5(nghl),fp6(nghl),fs1(nghl),fs2(nghl),fs3(nghl),fs4(nghl)  &
         ,fs5(nghl),fs6(nghl),wdcor(nghl),wdcori(nghl),cou(nghl),vDHartree(nghl,2) &
         ,vhart00(nghl,nghl),vhart01(nghl,nghl),vhart11(nghl,nghl))
    If(Allocated(aka)) Deallocate(aka,ro,tau,dro,dj,NABLAR,NABLAZ,SZFI,SFIZ,SRFI,SFIR)
    Allocate(aka(nghl,2),ro(nghl,2),tau(nghl,2),dro(nghl,2),dj(nghl,2)  &
         ,SZFI(nghl,2),SFIZ(nghl,2),SRFI(nghl,2),SFIR(nghl,2)  &
         ,NABLAR(nghl,2),NABLAZ(nghl,2))
    !-----------------------------------------
    ! Arrays depending on configurations
    !-----------------------------------------
    If(Allocated(rk)) Deallocate(rk,ak,qh,qh1,ql,ql1,nz,nr,nl,ns,npar,id  &
         ,ia,ikb,ipb,ka,kd,tb,txb,numax,ek,dk,vk,vk1,uk,vkmax,ddc,ddc1,hfb1,lcanon)
    Allocate(rk(nqx,nb2x),ak(nqx,nb2x),qh(0:nzx,1:ngh+1)  &
         ,qh1(0:nzx,1:ngh+1),ql(0:nrx,0:nlx,1:ngl+1),ql1(0:nrx,0:nlx,1:ngl+1)  &
         ,nz(ntx),nr(ntx),nl(ntx),ns(ntx),npar(ntx),id(nbx),ia(nbx),ikb(nbx),lcanon(0:nbx,2)  &
         ,ipb(nbx),ka(nbx,2),kd(nbx,2),tb(ntx),txb(nbx),numax(0:nkx,2)  &
         ,ek(nkx,2),dk(nkx,2),vk(nkx,2),vk1(nkx,2),uk(nkx,2),vkmax(nkx,2)  &
         ,ddc(ndx,nkx,2),ddc1(ndx,nkx,2),hfb1(nhfbx,2))
    !-----------------------------------------
    ! HFB Arrays
    !-----------------------------------------
    If(Allocated(erhfb)) Deallocate(erhfb,drhfb,erhfb1,drhfb1)
    Allocate(erhfb(nkx),drhfb(nkx),erhfb1(nkx),drhfb1(nkx))
    If(Allocated(hfb)) Deallocate(hfb,zhfb,evvk,hfbcan,evvkcan)
    Allocate(hfb(ndx2,ndx2),zhfb(ndx2),evvk(ndx2),hfbcan(ndx,ndx),evvkcan(ndx))
    If(Allocated(AN)) Deallocate(AN,ANk,PFIU,PFID,FIU,FID,FIUR,FIDR,FIUD2N,FIDD2N,FIUZ,FIDZ)
    Allocate(AN(nqx),ANk(nqx),PFIU(ndx),PFID(ndx),FIU(ndx),FID(ndx)  &
         ,FIUR(ndx),FIDR(ndx),FIUD2N(ndx),FIDD2N(ndx),FIUZ(ndx),FIDZ(ndx))
    !-----------------------------------------
    ! Optimal LAPACK storage
    !-----------------------------------------
    If(Allocated(alwork)) Deallocate(alwork,lwork)
#if(SWITCH_ESSL==0)
    ialwork=1+6*ndx+2*ndx**2; ilwork=3+5*ndx;
    Allocate(alwork(ialwork),lwork(ilwork));alwork = 0.0; lwork = 1
#else
    ialwork=0; ilwork=5*ndx;
    Allocate(alwork(1),lwork(ilwork));alwork = 0.0; lwork = 0
#endif
    !ialwork=1; ilwork=1;
    !If(Allocated(alwork)) Deallocate(alwork,lwork)
    !Allocate(alwork(ialwork),lwork(ilwork))
    !ier=0; Call DSYEVD('V','L',ndx2,hfb,ndx2,evvk,ALWORK,-1,LWORK,-1,ier)
    !If(ier.Ne.0) Then
    !   ierror_flag=ierror_flag+1
    !   ierror_info(ierror_flag)='STOP: FATAL ERROR CONDITION IN DSYEVD'
    !   Return
    !End If
    !ialwork=Int(alwork(1)); ilwork=lwork(1)
    !If(Allocated(alwork)) Deallocate(alwork,lwork)
    !Allocate(alwork(ialwork),lwork(ilwork))
    !-----------------------------------------
    ! Eqp, U,V
    !-----------------------------------------
    If(Allocated(RVqpN)) Deallocate(RVqpN,RVqpP,RUqpN,RUqpP,REqpN,REqpP)
    Allocate(RVqpN(nuv),RVqpP(nuv),RUqpN(nuv),RUqpP(nuv),REqpN(nqp),REqpP(nqp))
    If(Allocated(KpwiP)) Deallocate(KpwiP,KpwiN,KqpN,KqpP)
    Allocate(KpwiN(nqp),KpwiP(nqp),KqpN(nqp),KqpP(nqp))
    If(Allocated(fn_T)) Deallocate(fn_T,fp_T)
    Allocate(fn_T(nqp),fp_T(nqp))
    !-----------------------------------------
    ! PNP ARRAYS: CONF. AND GAUGE ANGLE
    !-----------------------------------------
    If(Allocated(exp1iphy))Deallocate(ropj,taupj,dropj,djpj,akapj,coupj,pjk  &
         ,SZFIpj,SFIZpj,SRFIpj,SFIRpj,epj,cpj,ypj,rpj,ddepj,phypj,sinphy  &
         ,exp1iphy,exp2iphy,exp1iphym,exp2iphym)
    Allocate(ropj(nghl,ilpj,2),taupj(nghl,ilpj,2),dropj(nghl,ilpj,2)  &
         ,djpj(nghl,ilpj,2),akapj(nghl,ilpj,2),coupj(nghl,ilpj),pjk(ilpj,2)  &
         ,SZFIpj(nghl,ilpj,2),SFIZpj(nghl,ilpj,2),SRFIpj(nghl,ilpj,2)  &
         ,SFIRpj(nghl,ilpj,2),epj(ilpj,2),cpj(nkx,ilpj,2),ypj(nkx,ilpj,2)  &
         ,rpj(nkx,ilpj,2),ddepj(nqx,ilpj,nb2x),phypj(ilpj),sinphy(ilpj),  &
         exp1iphy(ilpj),exp2iphy(ilpj),exp1iphym(ilpj),exp2iphym(ilpj))
    !-----------------------------------------
    ! FIELDS INITIALIZATION (NB! optimize)
    !-----------------------------------------
    ro=zero;     tau=zero;    dro=zero;    dj=zero;  aka=zero; rk=zero
    vn=zero;     vsn=zero;    vhbn=zero;   vrn=zero; vzn=zero; vdn=zero;
    vp=zero;     vsp=zero;    vhbp=zero;   vrp=zero; vzp=zero; vdp=zero;
    dvn=zero;    dvp=zero;
    vSFIZn=zero; vSZFIn=zero; vSFIRn=zero; vSRFIn=zero;  vDHartree=zero;
    vSFIZp=zero; vSZFIp=zero; vSFIRp=zero; vSRFIp=zero;
    ! Jason
    If(Allocated(allhfb)) Then
       Do ib=1,oldnb
          Deallocate(allhfb(ib)%arr,allevvk(ib)%arr,allalwork(ib)%arr,alllwork(ib)%arr)
       End Do
       Deallocate (allhfb,allevvk,allalwork,alllwork)
       Deallocate (allIALWORK,allILWORK,allISUPPZ)
    End If
    If (Allocated(allibro)) Deallocate(allibro)
    !
  End Subroutine thoalloc
  !=======================================================================
  !
  !=======================================================================
  Subroutine preparer(lpr)
    !---------------------------------------------------------------------
    ! setup routine
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Use HFBTHO_gauss
    Implicit None
    Logical      :: lpr
    Integer(ipr) :: iw,l,icount
    !
    If(n00.Eq.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)=' STOP: No more nuclei pass to the solver'
       Return
    End If
    !-----------------------------------------
    ! select the symbol of the nucleus
    !-----------------------------------------
    Call nucleus(1,npr(2),nucname)
    If(ierror_flag.Ne.0) Return
    !-----------------------------------------
    ! print headings to screen/'thoout.dat'
    !-----------------------------------------
    If(lpr) Then
       Call heading
       Call print_functional_parameters()
       Do iw=lout,lfile
          If(ierror_flag.Ne.0) Return
          !If(Print_HFBTHO_Namelist) Then
          !   Write(iw,'(100(2x,a,f15.8))')
          !   Write(iw,'(100(2x,a,f15.8))') 'NAMELIST CONTENT (copy/past to hfbtho_NAMELIST.dat and modify)'
          !   Write(iw,'(100(2x,a,f15.8))') '-------------------------------------------------------------'
          !   Write(iw,HFBTHO_GENERAL)
          !   Write(iw,HFBTHO_ITERATIONS)
          !   Write(iw,HFBTHO_FUNCTIONAL)
          !   Write(iw,HFBTHO_CONSTRAINTS)
          !   Write(iw,HFBTHO_BLOCKING)
          !   Write(iw,HFBTHO_PROJECTION)
          !End If
       End Do
    End If
    !-----------------------------------------
    ! pairing parameters (NB! modify later)
    !-----------------------------------------
    rho_c=0.160_pr; pwi=60.0_pr;
    !-----------------------------------------
    ! particle number as real variable
    !-----------------------------------------
    !tz(1)=Real(npr(1),Kind=pr); tz(2)=Real(npr(2),Kind=pr)
    amas=tz(1)+tz(2)
    drhoi=zero
    !-----------------------------------------
    ! default combinations
    !-----------------------------------------
    chargee2=e2charg
    coex=-chargee2*(three/pi)**p13; cex=-0.750_pr*coex
    !-----------------------------------------
    ! hbzero from forces [hqc**2/(two*amu)]
    !-----------------------------------------
    hb0=hbzero; If (use_cm_cor) hb0=hb0*(one-one/amas)
    !-----------------------------------------
    ! basis parameter q
    !-----------------------------------------
    beta0=q; q=Exp((3.0_pr*Sqrt(5.0_pr/(16.0_pr*pi)))*beta0)
    !-----------------------------------------
    ! basis parameters b0,bp,bz
    !-----------------------------------------
    If(b0.Le.zero) Then
       ! define oscillator frequency from default with empirical factor 1.2,
       ! and set length accordingly
       r00=r0*amas**p13; r02=r00**2; r04=r02**2
       hom=41.0_pr*amas**(-p13)*r0
       b0=Sqrt(two*hbzero/hom)
    Else
       ! define oscillator frequency from user-defined length, and set default
       ! empirical factor accordingly
       hom=hqc**2/(amn*b0**2)
       r0=(hom/41.0_pr)*amas**(p13)
       r00=r0*amas**p13; r02=r00**2; r04=r02**2
    End If
    bp=b0*q**(-one/6.0_pr); bz=b0*q**(one/3.0_pr); bpp=bp*bp
    ! overwrites with input from HFODD
    bp=bp_INI;bz=bz_INI;b0=b0_INI; bpp=bp*bp
    !-----------------------------------------
    ! constraint in terms of beta
    !-----------------------------------------
    ty20=Sqrt(5.0_pr/pi)*hom/b0**2/two
    !-----------------------------------------
    ! projection: number of grid points
    !-----------------------------------------
    keypj=Max(1,keypj); ilpj=keypj; ilpj2=ilpj**2;
    If(iproj.Eq.0) Then
       npr1pj=npr(1); npr2pj=npr(2)
    Else
       npr1pj=npr(1)+npr1pj; npr2pj=npr(2)+npr2pj
    End If
    !-----------------------------------------
    ! blocking window
    !-----------------------------------------
    pwiblo=Min(Max(25.0_pr/Sqrt(Real(npr(1)+npr(2),Kind=pr)),2.0_pr),8.0_pr)
    !-----------------------------------------
    ! THO
    !-----------------------------------------
    ass=zero; iasswrong=0
    !-----------------------------------------
    ! iterations
    !-----------------------------------------
    etot=zero; varmas=zero; rms=zero; ept=-two; del=one; alast=-seven; siold=one
    varmasNZ=zero; pjmassNZ=zero; ass=zero; skass=zero
    !---------------------------------------------------------
    ! statistics to screen('lout')/file('lfile')
    !---------------------------------------------------------
    If(lpr) Then
       Do iw=lout,lfile
          Write(iw,*)
          Write(iw,'(a)')             '  ---------------------------------------'
          Write(iw,'(a)')             '        Characteristics of the run       '
          Write(iw,'(a)')             '  ---------------------------------------'
          Write(iw,'(a,i5)')          '  Output file ................: ',lfile
          Write(iw,'(a,2x,a2,i4)')    '  Nucleus ....................: ',nucname,npr(1)+npr(2)
          Write(iw,'(a,i5)')          '  Number of HO shells ........: ',n00
          Write(iw,'(a,f20.14)')      '  HO length b0 (fm) ..........: ',b0
          Write(iw,'(a,f8.3,a,f8.3)') '  Basis deformation ..........:  beta0=',beta0,' q=',q
          Write(iw,'(a,5(1x,e15.8))') '  HO: b0,1/b0,bp,bz,q ........: ',b0,one/b0,bp,bz,q
          Write(iw,'(a,3(1x,e15.8))') '  h**2/(2m), cmc, e**2 .......: ',hbzero,hb0,chargee2
          Write(iw,'(a,2(1X,e15.8))') '  hom=f*41.0_pr*A^{-1/3}, f...: ',hom,r0
          If(iLST.Eq.0)  Then         ! HFB+HO case only
             iLST1=0
             Write(iw,'(a)')          '  THO basis is ...............:  OFF'
          Else                        ! HFB+THO case
             Write(iw,'(a)')          '  THO basis is ...............:   ON'
             If(iLST.Gt.0) Then       ! HFB+THO  only
                iLST1=1
                If(inin.Gt.0) Then
                   ierror_flag=ierror_flag+1
                   ierror_info(ierror_flag)=' Stop: Forbidden iLST>0, inin>0 combination.'
                   Return
                End If
                Write(iw,'(a)')       '    THO parameters from tholst.wel'
             Else                     ! HFB+THO after HFB+HO
                iLST1=0
                Write(iw,'(a)')       '    HFB+THO after a HFB+HO run    '
             End If
          End If
          Write(iw,'(a,i5)')          '  Maximal number of iterations: ',maxi
          Write(iw,'(a,f6.3)')        '  Initial mixing parameter ...: ',xmix
          If(inin.Eq.1)  Then
             Write(iw,'(a)')          '  Initial w.f. ...............:  from spherical scratch'
          End If
          If(inin.Eq.2)  Then
             Write(iw,'(a)')          '  Initial w.f. ...............:  from prolate scratch'
          End If
          If(inin.Eq.3)  Then
             Write(iw,'(a)')          '  Initial w.f. ...............:  from oblate scratch'
          End If
          If(inin.Lt.0) Then
             Write(iw,'(a)')          '  Initial wave functions from :  tape'
          End If
          Write(iw,'(a,3x,a)')        '  Skyrme functional ..........: ',skyrme
          If(icou.Eq.0) Write(iw,'(a)')    '    without Coulomb forces'
          If(icou.Eq.1) Write(iw,'(a)')    '    with direct Coulomb force only'
          If(icou.Eq.2) Write(iw,'(a)')    '    with direct and exchange Coulomb'
          If(kindhfb.Lt.0) Then
             Write(iw,'(a)')          '  Lipkin-Nogami procedure is .:   ON'
          Else
             Write(iw,'(a)')          '  Lipkin-Nogami procedure is .:  OFF'
          End If
          If(ilpj-1.Eq.0) Then
             Write(iw,'(a)')          '  PAV procedure is ...........:  OFF'
          Else
             Write(iw,'(a)')          '  PAV procedure is ...........:   ON'
             Write(iw,'(a,i5)')       '    Number of gauge points....: ',keypj
          End If
          If(icstr.Eq.0) Then
             Write(iw,'(a)')          '  Constraint calculation is ..:  OFF'
          Else
             Write(iw,'(a)')          '  Constraint calculation is ..:   ON'
             icount=0
             Do l=1,8
                If(Abs(lambda_active(l)).Gt.0) Then
                   icount=icount+1
                   Write(iw,'(a,i1,a,i1,a,f8.3)') '    Constraint ',icount,' .............: lambda=',l, &
                                                  ' Ql=',multRequested(l)
                End If
             End Do
          End If
          If(keyblo(1).Ne.0) Then
             Write(iw,'(a)')          '  Neutron blocking is ........:   ON'
          End If
          If(keyblo(2).Ne.0) Then
             Write(iw,'(a)')          '  Proton blocking is .........:   ON'
          End If
          If(switch_on_temperature) Then
             Write(iw,'(a,f6.2,a)')   '  Temperature T ..............: ',temper,' MeV'
          Else
             Write(iw,'(a,f6.2)')     '  Temperature T ..............:   0.00 MeV'
          End If
          Write(iw,'(a,i3)')          '  Restart indicator ..........: ',inin
          If(nbroyden.Eq.0) Then
             Write(iw,'(a,i3)')       '  Linear mixing ..............: ',nbroyden
          Else
             Write(iw,'(a,i3)')       '  Broyden mixing (#iterations): ',nbroyden
          End If
       End Do
    End If
    !-----------------------------------------
    ! BASIS, GAUSS POINTS, HOWF
    !-----------------------------------------
    Call gfv                      ! factorials
    If(ierror_flag.Ne.0) Return
    Call base0(lpr)               ! basis space (calculate configurational space)
    If(ierror_flag.Ne.0) Return
    Call thoalloc                 ! global allocation
    If(ierror_flag.Ne.0) Return
    Call gausspoints              ! GAUSS mesh points
    If(ierror_flag.Ne.0) Return
    Call base(lpr)                ! oscillator configurations (set up quantum numbers)
    If(ierror_flag.Ne.0) Return
    Call gaupol(lpr)              ! basis wf at gauss mesh points
    If(ierror_flag.Ne.0) Return
    !
  End Subroutine preparer
  !====================================================================
  !
  !====================================================================
  Subroutine coordinateLST(lpr)
    !------------------------------------------------------------------
    ! HO/THO
    !------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Use HFBTHO_THO, Only: f01234
    Implicit None
    Logical :: lpr
    Integer(ipr) :: i,il,ih
    If(iLST1.Eq.0) Then
       ! HO-basis
       Do il=1,ngl
          Do ih=1,ngh
             i=ih+(il-1)*ngh
             fh(i)=bz*xh(ih)
             fl(i)=bp*Sqrt(xl(il))
             wdcor(i)=pi*wh(ih)*wl(il)*bz*bp*bp
             wdcori(i)=one/wdcor(i)
          End Do
       End Do
    Else
       ! THO basis
       Call f01234(.False.)
       If(ierror_flag.Ne.0) Return
    End If
    !
    Call optHFBTHO                ! optimal HO/THO combinations
    If(ierror_flag.Ne.0) Return
    !
  End Subroutine coordinateLST
  !====================================================================
  !
  !====================================================================
  Subroutine iter(lpr)
    !------------------------------------------------------------------
    ! Iterations through successive diagonalisation
    !------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Logical :: lpr
    Real(pr)       :: assprn,delln(2)
    Real(pr), Save :: time
    Integer(ipr)   :: iw,it,ite
    Real(pr)       :: time1,time2,time3,time4,time5
    !---------------------------------------------------
    ! print to screen('lout')/thoout.dat('lfile')
    !---------------------------------------------------
    Do iw=lout,lfile
       If(iLST.Eq.0) Then
          Write(iw,'(a,f7.3,4(a,i3),a)')  &
               '  |HFB+HO> iterations(b0=',b0,', Nsh=',n00,  &
               ', inin=',inin,', N=',npr(1),', Z=',npr(2),')...'
       Else
          If(iLST1.Eq.0.Or.iasswrong(3).Ne.0) Then
             If(iasswrong(3).Ne.0) Then
                Write(iw,'(a,f7.3,a,i3,a)')  &
                     '  |HFB+THO substituted by HFB+HO> iterations (b0=',  &
                     b0,', Nsh=',n00,')...'
             Else
                Write(iw,'(a,f7.3,a)')'  towards |hfb+tho> iterations...'
                Write(iw,'(a,f7.3,a)')
                Write(iw,'(a,f7.3,a,i3,a)')  &
                     '  |Preliminary HFB+HO> iterations (b0=',b0,', Nsh=',n00,')...'
             End If
          Else
             If(itass.Eq.1) Then
                Write(iw,'(2(a,f7.3),a,i3,a)')  &
                     '  |HFB+THO> iterations(b0=',b0,', neutron density decay=',  &
                     decay,', Nsh=',n00,')...'
             Else
                Write(iw,'(2(a,f7.3),a,i3,a)')  &
                     '  |HFB+THO> iterations(b0=',b0,', proton density decay=',  &
                     decay,', Nsh=',n00,')...'
             End If
          End If
       End If
       Write(iw,1)
       Write(iw,'(20(a))') '  i','          si ','    mix ','  beta',    &
      &    '      Etot ','      A ','      rn','      rp ','        En', &
      &    '      Dn','      Ep','      Dp','        Ln  ','    Lp ',    &
      &    '    time '
       Write(iw,1)
1      Format(2x,130('-'))
    End Do
    !---------------------------------------------------------------------
    ! main hfb iteration loop
    !---------------------------------------------------------------------
    iError_in_HO=0; iError_in_THO=0; time=0.0_pr; time5=0.0_pr
    Do ite=1,maxi
       Call Cpu_time(time1)
       !
       iiter=ite
       !
       If (lpr.Or.iiter.Eq.1) Then
          assprn=ass(1); If(assprn.Gt.ass(2)) assprn=-ass(2) ! protons come with '-'
          delLN=del; If(kindhfb.Lt.0) delLN=del+ala2         ! LN case
          ! during iterations print
          Do iw=lout,lfile
             If(Max(Abs(drhoi(1)),Abs(drhoi(2))).Gt.1.0e-10_pr) Then
              !Write(*,*) '  WARNING! Int(Dro)=',Max(Abs(drhoi(1)),Abs(drhoi(2)))
             End If
             Write(iw,2) iiter,bbroyden,si,xmix,bet,etot,varmas,rms(1),rms(2),ept(1),delLN(1), &
                  ept(2),delLN(2),alast(1),alast(2),time
          End Do
       End If
       !-------------------------------------------------
       ! HFBDIAG
       !-------------------------------------------------
       If(IDEBUG.Gt.0) Call Cpu_time (time3)
       Do it=itmin,itmax
          Call hfbdiag(it,0)   ! hfb diagonalization with minimal canonical
          If(ierror_flag.Ne.0) Return
       End Do
       If(Print_Screen.And.IDEBUG.Gt.0) Then
          Call Cpu_time (time4)
          Write(*,*) '  Time in hfbdiag:',time4-time3,' seconds'
       End If
       !-------------------------------------------------
       ! EXPECT, DENSIT, COULOMB, FIELD, GAMDEL
       !-------------------------------------------------
       Call expect(.False.)    ! expectation values
       If (numberCons.Gt.0) Call getLagrange(ite)   ! new Lagrange parameters for constraints
       If(ierror_flag.Ne.0) Return
       Call field              ! new fields
       If(ierror_flag.Ne.0) Return
       Call gamdel             ! hf-matrix
       If(ierror_flag.Ne.0) Return
       !-------------------------------------------------
       ! Dumping control (old linear mixing)
       !-------------------------------------------------
       xmix0=0.1 !original 0.1
       If(si.Lt.siold) Then
          xmix=Min(xmax,xmix * 1.130_pr);  !old value 1.13
       Else
          xmix=xmix0
       End If
       siold=si
       !-------------------------------------------------
       ! time per iteration
       !-------------------------------------------------
       Call Cpu_time(time2)
       time=time2-time1; time5=time5+time
       !-------------------------------------------------
       ! Solution is OK within the iteration limit
       !-------------------------------------------------
       If(iiter.Ge.2.And.si.Lt.epsi) Then
          If(iLST1.Eq.0) Then
             iError_in_HO=0
          Else
             iError_in_THO=0
          End If
          ! iteration interrupted print
          If(.Not.lpr) Then
             delLN=del; If(kindhfb.Lt.0) delLN=del+ala2
             Do iw=lout,lfile
                Write(iw,3) iiter,bbroyden,si,xmix,bet,etot,varmas,rms(1),rms(2),ept(1),delLN(1), &
                     ept(2),delLN(2),alast(1),alast(2),time !Max(Abs(drhoi(1)),Abs(drhoi(2)))
                Write(iw,'(a,f8.3,a)') '  Total CPU time=',time5/60.0_pr,' minutes'
             End Do
          End If
          ! converged print
          Do iw=lout,lfile
             Write(iw,4) iiter,si,iError_in_HO,iError_in_THO
             Write(iw,'(a,f8.3,a)') '  Total CPU time=',time5/60.0_pr,' minutes'
          End Do
          iiter=iiter+1
          Return
       End If
       !-------------------------------------------------
       ! Slow convergence and lambda >0 (stop iterations)
       !-------------------------------------------------
       If(iiter.Ge.1000.And.(alast(1).Gt.zero.Or.alast(2).Gt.zero)) Exit
       !
    End Do    ! ite
    iiter=iiter+1
    !-------------------------------------------------
    ! Solution interrupted due to iterations limit
    !-------------------------------------------------
    If(iLST1.Eq.0) Then
       iError_in_HO=-1
    Else
       iError_in_THO=-1
    End If
    delLN=del; If(kindhfb.Lt.0) delLN=del+ala2
    ! iterations limit print
    Do iw=lout,lfile
       Write(iw,2) iiter,bbroyden,si,xmix,bet,etot,varmas,rms(1),rms(2),ept(1),delLN(1), &
            ept(2),delLN(2),alast(1),alast(2),Max(Abs(drhoi(1)),Abs(drhoi(2)))
       Write(iw,5) iiter,si,iError_in_HO,iError_in_THO
       Write(iw,'(a,f8.3,a)') '  Total CPU time=',time5/60.0_pr,' minutes'
    End Do
    !-------------------------------------------------
2   Format(i4,a,1x,f12.8,f5.2,f7.3,f13.6,1x,f6.1,2(f8.3),' | ',4(f8.3),' | ',20(f8.3))
3   Format(2x,130('-'),/,'  *   iteration interrupted after',i4,' steps   si=',f17.10,' ho=',i3,' tho=',i3,/,2x,130('-'))
4   Format(2x,130('-'),/,'  *   iteration converged   after',i4,' steps   si=',f17.10,' ho=',i3,' tho=',i3,/,2x,130('-'))
5   Format(2x,130('-'),/,'  *   iterations limit interrupt after',i4,' steps   si=',f17.10,' ho=',i3,' tho=',i3,/,2x,130('-'))
    !-------------------------------------------------
  End Subroutine iter
  !====================================================================
  !
  !====================================================================
  Subroutine hfbdiag(it,icanon)
    !------------------------------------------------------------------
    ! Skyrme-HFB diagonalization in axial HO/THO basis
    !------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
#if(USE_OPENMP==1)
    Use omp_lib
#endif
    Implicit None
    Logical :: lpr_pwi,norm_to_improve
    Character(Len=1) :: char1,char2,char3
    Integer(ipr) :: iw,it,i0,icanon,ibiblo,ier,i,j,k,k0,kl,lc,ib,nd,  &
                    nhfb,n1,n2,kaib,m,ndk,nd1,nd2,kdib,k1,k2,id1,id2, &
                    n12,n21,ntz,nhhph,nhhpp,ibro,ibroib,i_uv,i_eqp,jj,&
                    tid,IL,IU,NUMFOU,jlwork,jalwork,ldw,ldi,ii
    Real(pr) :: al,al2,emin,hla,dla,pn,eqpe,ela,enb,enb1,ekb, &
                s1,s2,s3,alnorm,sitest,fac1,fac2,fT,exponent, &
                VL,VU,ABSTOL,buffer
    Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
    Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:),f_T(:)
    Integer(ipr), Allocatable :: ISUPPZ(:),lwork_p(:)
    Real(pr), Allocatable :: alwork_p(:),eigenv(:),eigenf(:,:),hfbmat(:,:)
    Real(pr), External :: DLAMCH
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('hfbdiag',0)
    !
    If(it.Eq.1) Then
       EqpPo=>REqpN; VqpPo=>RVqpN; UqpPo=>RUqpN; KpwiPo=>KpwiN; KqpPo=>KqpN; f_T=>fn_T
    Else
       EqpPo=>REqpP; VqpPo=>RVqpP; UqpPo=>RUqpP; KpwiPo=>KpwiP; KqpPo=>KqpP; f_T=>fp_T
    End If
    KpwiPo=0; KqpPo=0; f_T=zero
    !
    nhhph=(it-1)*nhhdim; nhhpp=(it+1)*nhhdim
    If(.Not. Allocated(allhfb)) Then
       oldnb = nb ! for destroying data structures in next computation
       Allocate(allhfb((nb)),allevvk((nb)),allalwork((nb)),alllwork((nb)),allIALWORK(nb),allILWORK(nb))
       Allocate(allISUPPZ(nb))
       Do ib=1,nb
          nhfb=2*id(ib)
          Allocate(allhfb(ib)%arr(1:nhfb,1:nhfb))
          Allocate(allevvk(ib)%arr(1:nhfb))
          !jalwork=1+6*nhfb+2*nhfb**2; allIALWORK(ib)=jalwork  ! DSYEVD
          !jlwork=3+5*nhfb; allILWORK(ib)=jlwork               ! DSYEVD
#if(SWITCH_ESSL==0)
          jalwork=26*nhfb; allIALWORK(ib)=jalwork  ! DSYEVR
          jlwork=10*nhfb; allILWORK(ib)=jlwork     ! DSYEVR
#else
          jalwork=8*nhfb; allIALWORK(ib)=jalwork   ! DSYEVX
          jlwork=5*nhfb; allILWORK(ib)=jlwork      ! DSYEVX
#endif
          Allocate(allalwork(ib)%arr(1:jalwork))
          Allocate(alllwork(ib)%arr(1:jlwork))
          Allocate(allISUPPZ(ib)%arr(1:2*nhfb)) ! DSYEVR
       End Do
    End If
    !
    ABSTOL=2.0_pr*DLAMCH('S')
    !
    If (.Not. Allocated(allibro)) Then
       Allocate(allibro(1:NB))
       allibro(1)=0
       Do ib=2,NB
          allibro(ib) = allibro(ib-1) + (ID(ib-1)*(ID(ib-1)+1)/2)
       End Do
    End If
    !
    !------------------------------------------------------------------
    ! Loop the internal normalization
    !------------------------------------------------------------------
    !sitest=Max(Min(0.10_pr,si*0.010_pr),0.000010_pr)
    sitest=Min(0.10_pr,si*0.010_pr)
    norm_to_improve=.True.; inner(it)=-1; sumnz(it)=one
    Do While(norm_to_improve)
       !
       inner(it)=inner(it)+1
       !
       If(Abs(sumnz(it)).Lt.sitest.Or.inner(it).Eq.20) norm_to_improve=.False.
       !
       sumnz(it)=zero; entropy(it)=zero; v2min(it)=one; Dispersion(it)=zero
       !
       kl=0; emin=1000.0_pr; al=ala(it)
       !
       ! blocking
       If(iparenti(it).Eq.0) blomax(it)=0
       blo123d(it)=0; blok1k2d(it)=0; blocanon(it)=0;
       ibiblo=bloblo(keyblo(it),it)
       !------------------------------------------------------------------
       ! Runs over blocks
       !------------------------------------------------------------------
       i_uv=0; i_eqp=0
       lc=0; lcanon(0,it)=0; klmax=0
!$OMP Parallel Default(None) &
!$OMP& SHARED(nb,id,ia,it,nbx,allibro,brin,allhfb,allevvk, &
!$OMP&        allALWORK,allLWORK,allIALWORK,allILWORK,nhhph,nhhpp,al, &
!$OMP&        zhfb,ndx2,allISUPPZ,ABSTOL) &
!$OMP& PRIVATE(ib,nd,nhfb,i0,m,ibro,n1,nd1,n2,nd2,hla,dla,ier,tid,char1,char2, &
!$OMP&         NUMFOU,IL,IU,VL,VU,eigenf,eigenv,hfbmat,ISUPPZ,alwork_p,lwork_p, &
!$OMP&         ldw,ldi,char3,jalwork,jlwork)
#if(USE_OPENMP==1)
       tid = OMP_GET_THREAD_NUM()
#endif
!$OMP DO SCHEDULE(DYNAMIC)
       Do ib=1,nb
          nd=id(ib); nhfb=nd+nd; i0=ia(ib); m=ib+(it-1)*nbx; ibro=allibro(ib)
#if(SWITCH_ESSL==0)
          jalwork=26*nhfb; jlwork=10*nhfb
#else
          jalwork=8*nhfb; jlwork=5*nhfb
#endif
          allhfb(ib)%arr(1:nhfb,1:nhfb)=0.0_pr; allevvk(ib)%arr(1:nhfb)=0.0_pr
          allALWORK(ib)%arr(1:jalwork)=0.0_pr; allLWORK(ib)%arr(1:jlwork)=0;
          allISUPPZ(ib)%arr(1:2*nhfb)=0
          !------------------------------------------------------------------
          !  hfb-matrix
          !------------------------------------------------------------------
          Allocate(hfbmat(nhfb,nhfb))
          Do n1=1,nd
             nd1=n1+nd
             Do n2=1,n1
                nd2=n2+nd; ibro=ibro+1
                hla=brin(nhhph+ibro); dla=brin(nhhpp+ibro)
                hfbmat(n1,n2)=hla;    hfbmat(nd2,n1)=dla
                hfbmat(nd1,n2)=dla;   hfbmat(nd1,nd2)=-hla
             End Do
             hfbmat(n1,n1)  =hfbmat(n1,n1)  -al
             hfbmat(nd1,nd1)=hfbmat(nd1,nd1)+al
          End Do
          char1='V'; char2='I'; char3='L'; NUMFOU=0
          VL=0.0_pr ;VU=0.0_pr; IL=1; IU=nhfb; ldw=allIALWORK(ib); ldi=allILWORK(ib)
          Allocate(eigenv(nhfb)); eigenv(1:nhfb)=0.0_pr
          Allocate(eigenf(nhfb,nhfb)); eigenf(1:nhfb,1:nhfb)=0.0_pr
          Allocate(alwork_p(ldw)); Allocate(lwork_p(ldi))
#if(SWITCH_ESSL==0)
          ier=0; Allocate(ISUPPZ(2*nhfb))
          Call DSYEVR(char1,char2,char3,nhfb,hfbmat,nhfb,VL,VU,IL,IU,ABSTOL,NUMFOU,    &
                      eigenv,eigenf,nhfb,ISUPPZ,alwork_p,ldw,lwork_p,ldi,ier)
#else
          ier=0; Allocate(ISUPPZ(nhfb))
          Call DSYEVX('V','A','L',nhfb,hfbmat,nhfb,VL,VU,IL,IU,ABSTOL,NUMFOU, &
                      eigenv,eigenf,nhfb,alwork_p,ldw,lwork_p,ISUPPZ,ier)
#endif
          allevvk(ib)%arr(1:nhfb) = eigenv(1:nhfb)
          allhfb(ib)%arr(1:nhfb,1:nhfb) = eigenf(1:nhfb,1:nhfb)
          If(ier.NE.0) Then
             Write(6,*)'The algorithm failed to compute eigenvalues.'
#if(USE_OPENMP==1)
             Write(6,*)'I am',tid,' and I am working on array ',ib,ier
#endif
          End If
          Deallocate(eigenf,eigenv,hfbmat,ISUPPZ,alwork_p,lwork_p)
       End Do ! ib
!$OMP End Do
!$OMP End Parallel
       Do ib=1,NB
          nd=id(ib); nhfb=nd+nd; i0=ia(ib); m=ib+(it-1)*nbx; ibro=allibro(ib)
          !------------------------------------------------------------------
          ! Blocking
          !------------------------------------------------------------------
          ! external blocking
          If(iiter.Eq.1.And.inner(it).Eq.0) Then
             If(iparenti(it).Ne.0.And.keyblo(it).Eq.0) Then
                ! eventually charging
                !   keyblo(it)=1
                !   bloblo(keyblo(it),it)=ib
                !   blo123(keyblo(it),it)=requested level (k0)
                Call requested_blocked_level(ib,it)
                If(ierror_flag.Ne.0) Return
                ibiblo=bloblo(keyblo(it),it)
             End If
          End If
          ! general blocking
          k0=0
          If(ibiblo.Eq.ib) Then
             If(iiter.Eq.1.And.inner(it).Eq.0) Then
                ! blocked level as in the even-even nucleus
                k0=blo123(keyblo(it),it); ndk=k0+nd
                Do n2=1,nd
                   nd2=n2+nd
                   hfb1(n2,it)=allhfb(ib)%arr(n2,ndk)    !U
                   hfb1(nd2,it)=allhfb(ib)%arr(nd2,ndk)  !V
                End Do
                ! number of states in the block to be tested
                blocross(it)=Min(blomax(it)+10,nd)
             End If
             ! overlap between new and old blocked levels
             s3=zero
             Do n1=1,blocross(it)
                ndk=n1+nd; s1=zero
                Do n2=1,nd
                   nd2=n2+nd
                   s1=s1+Abs(hfb1(nd2,it)*allhfb(ib)%arr(nd2,ndk)) !VV
                   s1=s1+Abs(hfb1(n2,it)*allhfb(ib)%arr(n2,ndk))   !UU
                End Do
                If(s1.Gt.s3) Then
                   s3=s1; k0=n1
                End If
             End Do
             blo123d(it)=k0
             If(.Not.norm_to_improve) Then
                ! find maximal HO component
                ndk=k0+nd
                s1=zero
                Do n1=1,nd
                   nd1=n1+nd
                   hfb1(n1,it)=allhfb(ib)%arr(n1,ndk); hfb1(nd1,it)=allhfb(ib)%arr(nd1,ndk)
                   s2=Max(s1,Abs(allhfb(ib)%arr(n1,ndk)),Abs(allhfb(ib)%arr(nd1,ndk)))
                   If(s2.Gt.s1) Then
                      s1=s2; i=n1+i0  ! labels in k[k1,k2] numbering
                   End If
                End Do
                ! print blocked state
                Do iw=lout,lfile
                   Write(iw,'(4x,a,2(a,i3),2x,3(a,1x,f12.8,1x),(i3,a,i3,1x),a)')  &
                        protn(it),' Blocking: block=',ib,  &
                        ' state=',k0,  &
                        ' Eqp=',allevvk(ib)%arr(k0+nd),  &
                        ' Dqpe=',allevvk(ib)%arr(k0+nd)-eqpmin(it),  &
                        ' Ovlp=',s3  &
                        , keyblo(it),'/',blomax(it)  &
                        , tb(i)
                End Do
                ! ieresbl=6, 'BLKN','BLKZ'
                ereslbl(it)=tb(i)
                If(it.Eq.1) Then
                   ! 'BlEqpN','BlDEqpN','BlOvrN'
                   eresbl(1)=allevvk(ib)%arr(k0+nd); eresbl(2)=allevvk(ib)%arr(k0+nd)-eqpmin(it); eresbl(3)=s1
                Else
                   ! 'BlEqpZ','BlDEqpZ','BlOvrZ'
                   eresbl(4)=allevvk(ib)%arr(k0+nd); eresbl(5)=allevvk(ib)%arr(k0+nd)-eqpmin(it); eresbl(6)=s1
                End If
             End If
          End If
          !------------------------------------------------------------------
          ! Run over all qp states k in the block
          !------------------------------------------------------------------
          kaib=kl
          Do k=1,nd
             ndk=k+nd
             ! referent spectra
             pn=zero
             Do i=1,nd
                hla=allhfb(ib)%arr(i+nd,ndk)**2; pn=pn+hla
             End Do
             ! Blocking
             If(k.Eq.k0) Then
                n1=k0+nd
                Do i=1,nd
                   hla=allhfb(ib)%arr(i+nd,n1)**2; dla=allhfb(ib)%arr(i,n1)**2; pn=pn-half*(hla-dla)
                End Do
             End If
             eqpe=allevvk(ib)%arr(nd+k); ela=eqpe*(one-two*pn)
             enb=ela+al;                 ekb=Sqrt(Abs(eqpe**2-ela**2))
             i_eqp=i_eqp+1
             !------------------------------------------------------------------
             ! cut-off condition: energy pwi + Fermi cut-off function
             !------------------------------------------------------------------
                                                              exponent=Huge(1.0_pr)
             If(Abs(100.0_pr*(enb-pwi)).Lt.Log(Huge(1.0_pr))) exponent=Exp(100.0_pr*(enb-pwi))
             If(basis_HFODD) Then
                lpr_pwi=enb.Le.pwi !jacek sharp cut off for hfodd
             Else
                lpr_pwi=enb.Le.pwi.Or.Abs(one/(one+exponent)).Gt.cutoff_tol
             End If
             !------------------------------------------------------------------
             ! Remember the whole qp solution
             !------------------------------------------------------------------
             If(.Not.norm_to_improve) Then
                EqpPo(i_eqp)=eqpe                       ! Eqp_k
                If(lpr_pwi) KqpPo(kl+1)=i_eqp           ! below pwi otherwise zero
                If(lpr_pwi) KpwiPo(kl+1)=i_uv           ! below pwi otherwise zero
                Do n2=1,nd
                   nd2=n2+nd; i_uv=i_uv+1
                   UqpPo(i_uv)=allhfb(ib)%arr(n2,ndk)   ! U_ak
                   VqpPo(i_uv)=allhfb(ib)%arr(nd2,ndk)  ! V_ak
                End Do
             End If
             !------------------------------------------------------------------
             ! Define Fermi-Dirac occupations
             !------------------------------------------------------------------
             fT=zero
             If(switch_on_temperature.And.temper.Gt.1.D-12) Then
                fT = half*(one-Tanh(half*eqpe/temper))
                ! factor two comes from K>0 states only
                buffer = zero
                If(fT.Gt.zero.And.fT.Lt.one) Then
                   buffer = two*fT*Log(fT) + two*(one-fT)*Log(one-fT)
                End If
                entropy(it) = entropy(it) - buffer
                f_T(i_eqp) = fT
             End If
             !------------------------------------------------------------------
             ! Pairing window
             !------------------------------------------------------------------
             If(lpr_pwi) Then
                kl=kl+1                                      !number of active states
                If(k0.Eq.k) blok1k2d(it)=kl                  !blocking: dynamic #: k[k1,k2] numbering
                If((eqpe.Le.emin).And.(pn.Gt.0.0001_pr)) Then   !to avoid unocc at magic numbers
                   emin=eqpe; alnorm=pn                      !min qpe and its occupation
                End If
                erhfb(kl)=enb; drhfb(kl)=ekb; uk(kl,it)=pn     !ref.s.p. energies, deltas, occupancies
                sumnz(it)=sumnz(it)+two*pn+two*(one-two*pn)*fT !internal normalization
             End If
          End Do ! End k
          !
          If(.Not.norm_to_improve) Then
             !
             !------------------------------------------------------------------
             !  Density matrices (computed only when norm_to_improve = False,
             !  i.e., when the particle number is conserved within 10^-5)
             !------------------------------------------------------------------
             kdib=kl-kaib; ka(ib,it)=kaib; kd(ib,it)=kdib
             k1=kaib+1; k2=kaib+kdib
             eqpe=0.0_pr
             Do n2=1,nd
                Do n1=n2,nd
                   s1=zero; s2=zero
                   If(k1.Le.k2) Then
                      Do k=k1,k2
                         ! temperature
                         fac1 = one; fac2 = zero
                         If(switch_on_temperature) Then
                            !i_eqp=KqpPo(k); fac1=one-f_T(i_eqp); fac2=f_T(i_eqp)
                            ii=KqpPo(k); fac1=one-f_T(ii); fac2=f_T(ii)
                         End If
                         nd1=KpwiPo(k)+n1; nd2=KpwiPo(k)+n2
                         s1=s1+VqpPo(nd1)*fac1*VqpPo(nd2)+UqpPo(nd1)*fac2*UqpPo(nd2)
                         s2=s2+UqpPo(nd1)*fac1*VqpPo(nd2)+VqpPo(nd1)*fac2*UqpPo(nd2) &
                              +VqpPo(nd2)*fac1*UqpPo(nd1)+UqpPo(nd2)*fac2*VqpPo(nd1)
                      End Do
                      s1=two*s1; s2=half*s2               ! two:due to m-projection, half:due to symmetrization
                      ! blocking
                      If(ibiblo.Eq.ib) Then
                         i=blok1k2d(it); id1=KpwiPo(i)+n1; id2=KpwiPo(i)+n2
                         s1=s1-VqpPo(id1)*VqpPo(id2)+UqpPo(id1)*UqpPo(id2)
                         s2=s2-half*(UqpPo(id1)*VqpPo(id2)+VqpPo(id1)*UqpPo(id2))
                      End If
                   End If
                   n12=n1+(n2-1)*nd; n21=n2+(n1-1)*nd
                   rk(n12,m)=s1; rk(n21,m)=s1              !  V V'
                   ak(n12,m)=-s2; ak(n21,m)=-s2            !- U V', ak=half*(pairing density)
                   hfbcan(n1,n2)=s1; allhfb(ib)%arr(n1,n2)=s1
                End Do !n1
             End Do !n2
             !------------------------------------------------------------------
             ! Canonical basis
             !------------------------------------------------------------------
             If(k1.Le.k2) Then
                Call Canonical(it,icanon,k2,k1,nd,i0,lc,ib,ibiblo,m,ibro)
                If(ierror_flag.Ne.0) Return
             End If
             lcanon(ib,it)=lc
             !
          End If
          !
       End Do !ib
       !
       If(kl.Eq.0) Then
          ierror_flag=ierror_flag+1
          ierror_info(ierror_flag)=' STOP: kl=zero, no states below pwi!!!'
          Return
       End If
       If(iparenti(it).Ne.0.And.ibiblo.Eq.0) Then
          ierror_flag=ierror_flag+1
          ierror_info(ierror_flag)='STOP: No blocking candidate found!!!'
          Return
       End If
       eqpmin(it)=emin; klmax(it)=kl; sumnz(it)=sumnz(it)-tz(it)
       !------------------------------------------------------------------
       ! Lambda search
       !------------------------------------------------------------------
       Call ALambda(al,it,kl)
       If(ierror_flag.Ne.0) Return
       If(keyblo(it).Eq.0) Then
          ala(it)=al
       Else
          ala(it)=ala(it)+0.50_pr*(al-ala(it))
       End If
       ! NB! 'alast' instead of 'al' at small pairing
       alast(it)=al
       If(Abs(ept(it)).Lt.0.0001_pr.And.(.Not.switch_on_temperature)) Then
          ntz=Int(tz(it)+0.1_pr); ntz=ntz/2
          Do k=1,kl
             drhfb(k)=erhfb(k)
          End Do
          Call ord(kl,drhfb)
          alast(it)=drhfb(ntz)  !last bound s.p. energy
       End If
       !------------------------------------------------------------------
       ! THO asymptotic decay
       !------------------------------------------------------------------
       ! density asymptotic decay \rho(r)->Exp(-ass(it)*r)
       ! ass(it)=2*Sqrt((E_min-\lambda)/((A-1)/A)*hbar**2/(2*m)))
       al2=zero
       If(kindhfb.Lt.0) Then
          al2=al+two*ala2(it)*(one-two*alnorm) ! al=al+two*ala2(it)
       End If
       al2=(emin-al2)/hb0
       ! wrong asymptotic
       iasswrong(it)=0; If(al2.Le.zero) iasswrong(it)=1; ass(it)=two*Sqrt(Abs(al2))
       !
    End Do ! While(norm_to_improve)
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('hfbdiag',1)
    !
  End Subroutine hfbdiag
  !=======================================================================
  !
  !=======================================================================
  Subroutine ALambda(al,it,kl)
    !---------------------------------------------------------------------
    ! Adjusting Fermi energy
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: it,i,k,kl,icze,lit,ntz,iw
    Real(pr), Save :: fm7=1.0d-7,fm10=1.0d-10
    Real(pr) :: al,vh,xinf,xsup,esup,ez,dez,dfz,dvh,y,a,b,einf,absez,sn
    Real(pr) :: fT,dfT
    Real(pr), Pointer :: f_T(:)
    !-------------------------------------------------
    ! Fermi-Dirac occupations
    !-------------------------------------------------
    If(switch_on_temperature) Then
       If(it.Eq.1) Then
          f_T=>fn_T
       Else
          f_T=>fp_T
       End If
    End If
    !-------------------------------------------------
    ! Chemical potential without pairing
    !-------------------------------------------------
    If(CpV0(it-1).Eq.zero) Then
       ntz=Int(tz(it)+0.1_pr); ntz=ntz/2
       Do k=1,kl
          drhfb(k)=erhfb(k)
       End Do
       Call ord(kl,drhfb)
       If (ntz.Lt.kl) Then
          al=half*(drhfb(ntz)+drhfb(ntz+1))
       Else
          al=drhfb(ntz)+0.001_pr
       End If
       Return
    End If
    !-------------------------------------------------
    ! Chemical potential with pairing
    !-------------------------------------------------
    xinf=-1000.0_pr; xsup=1000.0_pr; esup=one; icze=0
    Do lit=1,500
       sn=zero;dez=zero;dfz=zero
       Do i=1,kl
          vh=zero; dvh=zero; fT=zero; dfT=zero
          y=erhfb(i)-al; a=y*y+drhfb(i)**2; b=Sqrt(a)
          !
          If(switch_on_temperature.And.temper.Gt.1.D-12) Then
             fT =half*(one-Tanh(half*b/temper))
             dfT=y/b/temper*fT*(one-fT)
             f_T(i)=fT
          Else
             fT =zero
             dfT=zero
          End If
          !
          If(b.Gt.zero)  vh=half*(one-y/b)
          !
          !If(b.Lt.fm7.And.icze.Eq.1) vh=-einf/(esup-einf) !no pairing
          If(vh.Lt.1.D-12) vh = zero
          If((vh-one).Gt.1.D-12)  vh = one
          If(b.Gt.zero) dvh=half*drhfb(i)**2/(a*b)         ! D[ez,al](i)
          ! blocking
          If(i.Eq.blok1k2d(it)) Then
             vh=half; dvh=zero
          End If
          sn=sn+two*vh+two*(one-two*vh)*fT
          dez=dez+two*(one-two*fT)*dvh
          dfz=dfz+two*(one-two*vh)*dfT   ! D[ez,al]
       End Do
       ez=sn-tz(it); absez=Abs(ez)/tz(it)
       dez=dez+dfz
       !-------------------------------------------------
       ! Correcting bounds
       !-------------------------------------------------
       If(ez.Lt.zero) Then
          xinf=Max(xinf,al); einf=ez
       Else
          xsup=Min(xsup,al); esup=ez
       End If
       If(lit.Eq.1) Then
          If(absez.Le.0.10_pr) Then
             al=al-ez
          Else
             al=al-0.10_pr*Sign(one,ez)
          End If
       Else
          al=al-ez/(dez+1.d-20)                         ! newton method
       End If
       If(xsup-xinf.Lt.fm7) icze=1                      ! low/upp close
       If(al.Lt.xinf.Or.al.Gt.xsup) al=half*(xinf+xsup) ! mean upp/low
       If(absez.Le.fm10) Return
    End Do
    !-------------------------------------------------
    ! Low accuracy warning
    !-------------------------------------------------
    Do iw=lout,lfile
      Write(iw,'(a,2(e12.5,2x),a,2(2x,f8.4),a,i2)') ' Low accuracy=',sn,ez,' for N,Z=',tz,' it=',it
    End Do
  End Subroutine Alambda
  !=========================================================================================
  !
  !=========================================================================================
  Subroutine Canonical(it,icanon,k2,k1,nd,i0,lc,ib,ibiblo,m,ibroib)
    !---------------------------------------------------------------------------------------
    ! Canonical diagonalization
    !---------------------------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: it,i0,icanon,ibiblo,i,iw,k,kk,lc,ib,nd,n1,n2,m,nd1,k1,k2,n12,ier
    Integer(ipr) :: nhhph,nhhpp,ibro,ibroib,il,iu,NUMFOU
    Integer(ipr), Allocatable :: ifail(:)
    Real(pr) :: s1,s2,vx,h1,d1,h2,d2,ddn1,ddn2,vl,vu,ABSTOL
    Real(pr), Allocatable :: hh(:,:),de(:,:),ewavef(:,:)
    Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:)
    Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
    Integer(ipr), Allocatable :: ISUPPZ(:)
    Real(pr), Allocatable :: eigenv(:),eigenf(:,:)
    Real(pr), External :: DLAMCH
    !
    If(it.Eq.1) Then
       EqpPo=>REqpN; VqpPo=>RVqpN; UqpPo=>RUqpN; KpwiPo=>KpwiN; KqpPo=>KqpN
    Else
       EqpPo=>REqpP; VqpPo=>RVqpP; UqpPo=>RUqpP; KpwiPo=>KpwiP; KqpPo=>KqpP
    End If
    !
#if(SWITCH_ESSL==0)
    If(Allocated(alwork)) Deallocate(alwork,lwork)
    ialwork=26*nd; ilwork=10*nd
    Allocate(ALWORK(ialwork),LWORK(ilwork)); ALWORK = 0.0_pr; LWORK = 1
#else
    If(Allocated(alwork)) Deallocate(alwork,lwork)
    ialwork=8*nd; ilwork=5*nd;
    Allocate(ALWORK(ialwork),LWORK(ilwork)); ALWORK = 0.0_pr; LWORK = 0
#endif
    !
    ABSTOL=2.0_pr*DLAMCH('S')
    !
    If(Abs(ept(it)).Lt.0.0001_pr.And.(.Not.switch_on_temperature)) Then
       !------------------------------------------------------
       ! No pairing => just taking the HF states
       !------------------------------------------------------
       Do k=1,nd
          kk=k1+k-1; lc=lc+1                              ! total number of the canonical states
          ddc(1:nd,lc,it)=zero; vk(lc,it)=zero            ! zeros: nd could be larger then k2-k1+1
          If(kk.Gt.k2) Cycle
          vx=zero
          Do i=1,nd
             h1=VqpPo(KpwiPo(kk)+i)**2; vx=vx+h1
          End Do
          If (vx.Le.zero) vx=zero                         ! roundoff errors
          If (vx.Ge.one ) vx=one
          Do i=1,nd
             If(vx.Ge.half) Then
                ddc(i,lc,it)=VqpPo(KpwiPo(kk)+i)          ! (ph) s.p. orbitals in conf.space
             Else
                ddc(i,lc,it)=UqpPo(KpwiPo(kk)+i)          ! (ph) s.p. orbitals in conf.space
             End If
          End Do
          Dispersion(it)=Dispersion(it)+four*vx*(one-vx)  ! internal P/N Dispersion
          If(Abs(vx-half).Le.v2min(it)) Then
             v2min(it)=Abs(vx-half); v2minv(it)=vx        ! divergent condition
             lcc=lc
          End If
          vk(lc,it)=vx                                    ! (ph) s.p. occupations v^2
          !------------------------------------------------------
          ! RESU only
          !------------------------------------------------------
          If(icanon.Ne.0) Then
             ek(lc,it)=EqpPo(KqpPo(kk))*(one-two*vx)+ala(it)      ! (ph) s.p. energies
             dk(lc,it)=zero                                       ! (ph) s.p. deltas
          End If
       End Do !k
    Else
       !------------------------------------------------------
       ! Pairing => calculate canonical basis
       !------------------------------------------------------
#if(SWITCH_ESSL==0)
       VL=0.0_pr; VU=0.0_pr; IL=1; IU=nd; NUMFOU=0
       Allocate(ISUPPZ(2*nd))
       Allocate(eigenv(nd)); eigenv(1:nd)=0.0_pr
       Allocate(eigenf(nd,nd)); eigenf(1:nd,1:nd)=0.0_pr
       ier=0; Call DSYEVR('V','A','L',nd,hfbcan,ndx,VL,VU,IL,IU,ABSTOL,NUMFOU, &
                          eigenv,eigenf,nd,ISUPPZ,ALWORK,ialwork,LWORK,ilwork,ier)
       evvkcan(1:nd) = eigenv(1:nd)
       hfbcan(1:nd,1:nd) = eigenf(1:nd,1:nd)
       Deallocate(eigenv,eigenf,ISUPPZ)
       !ier=0; Call DSYEVD('V','L',nd,hfbcan,ndx,evvkcan,ALWORK,ialwork,LWORK,ilwork,ier)
#else
       vl=0.0_pr; vu=0.0_pr; il=1; iu=1; m=0; abstol=0.0_pr
       Allocate(ifail(nd),ewavef(ndx,ndx))
       ier=0; Call DSYEVX('V','A','L',nd,hfbcan,ndx,vl,vu,il,iu,abstol,m, &
                          evvkcan,ewavef,ndx,ALWORK,ialwork,LWORK,ifail,ier)
       hfbcan(1:ndx,1:ndx) = ewavef(1:ndx,1:ndx)
       Deallocate(ifail,ewavef)
#endif
       ! bug in LAPACK
       If(ier.Gt.0) Then
          Do iw=lout,lfile
             Write(iw,*) 'FATAL ERROR CONDITION IN CANONICAL DSYEVR, ier=',ier,'(RECOVERED)'
          End Do
          Do n2=1,nd
             Do n1=n2,nd
                vx=allhfb(ib)%arr(n1,n2)
                hfbcan(n2,n1)=vx; hfbcan(n1,n2)=vx
             End Do
          End Do
          Call sdiag(ndx,nd,hfbcan,evvkcan,hfbcan,zhfb,+1)
       End If
       !------------------------------------------------------
       ! Eigenvalues and wavefunctions
       !------------------------------------------------------
       Do k=1,nd
          lc=lc+1                                        ! total number of the canonical states
          Do i=1,nd
             ddc(i,lc,it)=hfbcan(i,k)                    ! (ph) canon orbitals in conf.space
          End Do
          vx=evvkcan(k)*half
          If (vx.Le.zero) vx=zero                        ! roundoff errors
          If (vx.Ge.one ) vx=one
          ! blocking
          If(ibiblo.Eq.ib.And.vx.Gt.0.49_pr.And.vx.Le.0.51_pr) blocanon(it)=lc
          Dispersion(it)=Dispersion(it)+four*vx*(one-vx) ! internal P/N Dispersion
          If(Abs(vx-half).Le.v2min(it)) Then
             v2min(it)=Abs(vx-half); v2minv(it)=vx       ! divergent condition
             lcc=lc
          End If
          vk(lc,it)=vx                                   ! (ph) canon occupations v^2
          !------------------------------------------------------
          ! RESU only
          !------------------------------------------------------
          If(icanon.Ne.0) Then
             ! canonical energies and deltas (no physical meaning in PNP)
             nhhph=(it-1)*nhhdim; nhhpp=(it+1)*nhhdim
             Allocate(hh(nd,nd),de(nd,nd))
             ibro=ibroib
             Do n1=1,nd
                Do n2=1,n1
                   ibro=ibro+1
                   vx=brin(nhhph+ibro); hh(n2,n1)=vx; hh(n1,n2)=vx
                   vx=brin(nhhpp+ibro); de(n2,n1)=vx; de(n1,n2)=vx
                End Do
             End Do
             h1=zero; d1=zero
             Do n2=1,nd
                h2=zero; d2=zero
                Do n1=1,nd
                   ddn1=hfbcan(n1,k)
                   h2=h2+ddn1*hh(n1,n2)
                   d2=d2+ddn1*de(n1,n2)
                End Do
                ddn2=hfbcan(n2,k)
                h1=h1+h2*ddn2
                d1=d1+d2*ddn2
             End Do
             ek(lc,it)=h1                              ! (ph) canon s.p. energies
             dk(lc,it)=d1                              ! (ph) canon s.p. deltas
             Deallocate(hh,de)
          End If
          !
       End Do !k
    End If
    !------------------------------------------------------
    ! RESU only
    !------------------------------------------------------
    If(icanon.Ne.0) Then
       !------------------------------------------------------
       ! Find maximal HO components of all qp states
       !------------------------------------------------------
       Do k=k1,k2
          s1=zero
          Do n1=1,nd
             nd1=nd+n1
             s2=Max(s1,Abs(VqpPo(KpwiPo(k)+n1)),Abs(UqpPo(KpwiPo(k)+n1)))
             If(s2.Gt.s1) Then
                s1=s2
                vkmax(k,it)=s1                       ! maximal overlap
                numax(k,it)=n1+i0                    ! its number in k[k1,k2] numbering
             End If
          End Do
       End Do
       !------------------------------------------------------
       ! Searching for possible blocking candidates
       !------------------------------------------------------
       If(iparenti(it).Eq.0) Then
          n1=0
          Do k=k1,k2
             n1=n1+1
            ! Search within |(1-2*N)*Eqpe| lover than 'pwiblo'
            ! The levels number n1 is 1,2,3,... for the given block ([123] numbering)
             If(Abs(EqpPo(KqpPo(k))-eqpmin(it)).Le.pwiblo) Then
                blomax(it)=blomax(it)+1                         ! blocked state #, maximel # of block candidates
                If(blomax(it).Gt.bloall) Then
                   ierror_flag=ierror_flag+1
                   ierror_info(ierror_flag)='Too many blocking candidates! Increase bloall and run again'
                   Return
                End If
                bloblo(blomax(it),it)=ib                        ! block where to block
                blo123(blomax(it),it)=n1                        ! state # [123] numbering
                blok1k2(blomax(it),it)=k                        ! state # k[k1,k2] numbering
                bloqpdif(blomax(it),it)=Abs(EqpPo(KqpPo(k))-eqpmin(it))
             End If
          End Do
       End If
    End If
    !
  End Subroutine Canonical
  !=======================================================================
  !
  !=======================================================================
  Subroutine resu(irecord)
    !---------------------------------------------------------------------
    ! prints results: single particle energies, densities, fields
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: it,iw,ib,im,m,nd,k,k0,k1,k2,j,n,imax,nhfb,irecord
    Real(pr)     :: sum,eqpe,pn,ela,enb,ek0,vk0,ekk,delb,ovmax,s,uuvv,  &
         dk0,skk,summ(4),vvs,vvc,enjacek
    Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:)
    Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
    !
    !--------------------------------------------
    ! last HFB run for full canon.calculations
    !--------------------------------------------
    Do it=itmin,itmax
       Call hfbdiag(it,1)       ! hfb with maximal canonical
       If(ierror_flag.Ne.0) Return
    End Do
    !PAV
    Call expect(.False.)        ! expectation values
    !
    !Call coulom_test
    !
    If(ierror_flag.Ne.0) Return
    Call field                  ! new fields
    If(ierror_flag.Ne.0) Return
    Call gamdel                 ! hf-matrix
    If(ierror_flag.Ne.0) Return
    ! inout(2): HFB matrices if nucleus is even-even
    !If(npr(1).Eq.2*(npr(1)/2).And.npr(2).Eq.2*(npr(2)/2)) Then
       Call inout(2)
       If(ierror_flag.Ne.0) Return
    !End If
    !--------------------------------------------
    ! Printing densities and fields
    !--------------------------------------------
    ! Call printLST(ro(1,1),ro(1,2))      !need fix
    ! Call printLST(tau(1,1),tau(1,2))
    ! Call printLST(dro(1,1),dro(1,2))
    ! Call printLST(dj (1,1),dj (1,2))
    ! printing of fields
    ! Call printLST(vhb(1,1),vhb(1,2))
    ! Call printLST(v (1,1),v (1,2))
    ! Call printLST(vs(1,1),vs(1,2))
    !--------------------------------------------
    ! Printing quasiparticle states
    !--------------------------------------------
    Do it=itmin,itmax
       If(it.Eq.1) Then
          EqpPo=>REqpN; VqpPo=>RVqpN; UqpPo=>RUqpN; KpwiPo=>KpwiN; KqpPo=>KqpN
       Else
          EqpPo=>REqpP; VqpPo=>RVqpP; UqpPo=>RUqpP; KpwiPo=>KpwiP; KqpPo=>KqpP
       End If
       !
       If(Print_Screen) Then
          iw=lfile
          Write(iw,200) tit(it)
          Write(iw,*) ' eqp(k) -> q.p. energy '
          Write(iw,*) ' e(k)   -> referent s.p. energy '
          Write(iw,*) ' p(k)   -> occ.probability '
          Write(iw,*) ' del(k) -> referent s.p. gap '
          Write(iw,*) ' fermi energy alast=',alast(it)
          Write(iw,'(a,a)')  &
               '  #k  block#    eqp(k)     e(k)       (1-2N)E      decay        p(k)',  &
               '        del(k)    overl      labels'
200     Format(//,' #quasiparticle energies ',a,/,1x,32('-'))
       End If
       sum=zero
       Do ib=1,nb
          nd=id(ib); im=ia(ib); m=ib+(it-1)*nbx; nhfb=nd+nd
          k1=ka(ib,it)+1
          k2=ka(ib,it)+kd(ib,it)
          If(k1.Le.k2) Then
             Do k=k1,k2                                 ! print active states only
                pn=uk(k,it)                             ! qp probabilities
                j=k
                If(pn.Gt.-1.d-14) Then                   ! print If signIficant pn
                   ! main oscillator component
                   ovmax=vkmax(k,it)                    ! maximal overlap
                   imax=numax(k,it)                     ! its number
                   ! printing
                   eqpe=EqpPo(KqpPo(k))                 ! qp energies
                   skk=two*Sqrt(Abs(eqpe-ala(it))/hb0)  ! qp decay
                   ela=eqpe*(one-two*pn)
                   enb=ela+ala(it)                      ! ref. s.p. energies
                   delb=Sqrt(Abs(eqpe**2-ela**2))       ! ref. s.p. delta
                   sum=sum+two*pn                       ! particle number
                   If(Print_Screen) Then
                      iw=lfile
                      Write(iw,201) k,ib,eqpe,enb,(one-two*pn)*eqpe,skk,pn,delb,ovmax,tb(imax)
201                 Format(i4,2x,i3,1x,f12.6,f12.6,f12.6,f12.6,2x,f12.8,  &
                           2(2x,f7.4),' ',a13)
                   End If
                End If
             End Do
          End If
       End Do !ib
       !--------------------------------------------
       ! Printing canonical single particle states
       !--------------------------------------------
       If(Print_Screen) Then
          iw=lfile
          Write(iw,'(a,i4,a,i4)')  &
               '#all active are ',j,' q.p. states out of ',nt
          Write(iw,'(a,f6.1)') '#since the cut off is pwi=',pwi
          Write(iw,'(3a,f6.1)')'#check: number of ',tit(it),'=',sum
          Write(iw,100) tit(it)
          Write(iw,*) ' labels -> {2*omega}{parity}[nn=nz+2*nr+nl,nz,nl]'
          Write(iw,*) ' cqpe   -> canonical q.p. energies'
          Write(iw,*) ' ce     -> canonical s.p. energies'
          Write(iw,*) ' fermi energy=',alast(it)
          Write(iw,*) ' average cdelt=',del(it)
          Write(iw,'(a,a)')'  k0      ceqp        ce         v*v',  &
               '       u*v        cdel     overl      labels'
100     Format(//,' #canonical s.p. energies ',a,/,1x,33('-'),//)
       End If
       k0=0
       summ=zero; enjacek=zero
       Do ib=1,nb
          nd=id(ib); im=ia(ib)
          k1=ka(ib,it)+1; k2=ka(ib,it)+kd(ib,it)
          If(k1.Le.k2) Then
             Do k=1,nd
                k0=k0+1
                ! for Lipkin Nogami
                vvs=two*Sqrt(vk(k0,it))*Sqrt(one-vk(k0,it))    !2vu
                vvc=two*vk(k0,it)-one                          !2v^2-1
                summ(1)=summ(1)+vvs**2
                summ(2)=summ(2)+vvs**2*vvc
                summ(3)=summ(3)+vvs**4
                summ(4)=summ(4)+(vvs*vvc)**2
                ! search for main oscillator component
                ovmax=zero
                Do n=1,nd
                   s=Abs(ddc(n,k0,it))                         !canon orbitals in conf.space
                   If (s.Ge.ovmax) Then
                      ovmax=s; imax=n
                   End If
                End Do
                ! printing
                ek0=ek(k0,it)                                  !canon s.p. energies
                enjacek=enjacek+ek0*vk(k0,it)
                If(ek0.Lt.pwi) Then                            !print up to 'pwi'
                   vk0=vk(k0,it)                               !canon occupations v^2
                   If(vk0.Gt.-1.d-4) Then                       !print If signIficant v^2
                      dk0=-dk(k0,it)                           !canon s.p. deltas
                      ekk=Sqrt((ek0-ala(it))**2+dk(k0,it)**2)  !resulting cqpe
                      uuvv=Sqrt(Abs(vk0*(one-vk0)))            !resulting u*v
                      If(Print_Screen) Then
                         iw=lfile
                         Write(iw,101) k0,ekk+ala(it),ek0,vk0,uuvv,dk0,ovmax,tb(im+imax)
101                    Format(i4,2f12.6,2(1x,f12.8),2(2x,f7.4),' ',a13)
                      End If
                   End If
                End If
             End Do !k0
          End If
       End Do !ib
       !--------------------------------------------
       ! Lipkin-Nogami
       !--------------------------------------------
       ssln(1,it)=summ(1)
       ssln(2,it)=summ(2)
       ssln(3,it)=summ(4)*summ(1)-summ(2)**2+summ(1)**3/4.0_pr-half*summ(3)*summ(1)
       If(Print_Screen) Then
          iw=lfile
          Write(iw,*) ' Sum canonical e_v*V^2_k=',two*enjacek
       End If
    End Do !it
    !--------------------------------------------
    ! To thoout.dat, thodef.dat and hodef.dat
    !--------------------------------------------
    If(irecord.Ne.0) Then
       iappend=1
       Call expect(.True.)    !print  & record HFB+PAV results
       If(ierror_flag.Ne.0) Return
       iappend=0
    Else
       Call expect(.True.)    !print HFB+PAV results
       If(ierror_flag.Ne.0) Return
    End If
    !
  End Subroutine resu
  !=======================================================================
  !
  !=======================================================================
  Subroutine initialize_HFBTHO_NAMELIST
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    !------------------------------------
    ! Namelist (default values)
    !------------------------------------
    ! HFBTHO_GENERAL
    number_of_shells    = 10
    oscillator_length   =-one
    basis_deformation   = zero
    proton_number       = 24
    neutron_number      = 26
    type_of_calculation = 1
    ! HFBTHO_ITERATIONS
    number_iterations   = 100
    accuracy            = 1.D-5
    restart_file        = -1
    ! HFBTHO_FUNCTIONAL
    functional          = 'SLY4'
    add_initial_pairing = .False.
    type_of_coulomb     = 2
    ! HFBTHO_PAIRING
    user_pairing        = .False.
    vpair_n             = -300.0_pr
    vpair_p             = -300.0_pr
    pairing_cutoff      =   60.0_pr
    pairing_feature     =    0.5_pr
    ! HFBTHO_CONSTRAINTS
    lambda_values       = (/ 0, 0, 0, 0, 0, 0, 0, 0 /)
    lambda_active       = (/ 0, 0, 0, 0, 0, 0, 0, 0 /)
    expectation_values  = (/ 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr /)
    ! HFBTHO_BLOCKING
    proton_blocking     = (/ 0, 0, 0, 0, 0 /)
    neutron_blocking    = (/ 0, 0, 0, 0, 0 /)
    ! HFBTHO_PROJECTION
    switch_to_THO       = 0
    projection_is_on    = 0
    gauge_points        = 1
    delta_Z             = 0
    delta_N             = 0
    ! HFBTHO_TEMPERATURE
    set_temperature     = .False.
    temperature         = zero
    ! HFBTHO_DEBUG
    number_Gauss        =  40
    number_Laguerre     =  40
    number_Legendre     =  80
    compatibility_HFODD = .False.
    number_states       = 500
    force_parity        = .True.
    print_time          = 0
    !
  End Subroutine initialize_HFBTHO_NAMELIST
  !=======================================================================
  !
  !=======================================================================
  !Subroutine read_HFBTHO_NAMELIST
  !  Use HFBTHO_utilities
  !  Use HFBTHO
  !  Implicit None
  !  Integer(ipr) :: ios,lnamelist=16
  !  !------------------------------------
  !  ! Namelist (handling)
  !  !------------------------------------
  !  Open(lnamelist,file='hfbtho_NAMELIST.dat',DELIM='APOSTROPHE') ! 'QUOTE'
  !  !
  !  ierror_flag = 0
  !  !
  !  ! General input data
  !  Read(UNIT=lnamelist,NML=HFBTHO_GENERAL,iostat=ios)
  !  If (ios.Ne.0) Then
  !     ierror_flag=ierror_flag+1
  !     ierror_info(ierror_flag)='Error in HFBTHO_GENERAL read'
  !     Return
  !  End If
  !  !
  !  ! Iterations
  !  Read(UNIT=lnamelist,NML=HFBTHO_ITERATIONS,iostat=ios)
  !  If (ios.Ne.0) Then
  !     ierror_flag=ierror_flag+1
  !     ierror_info(ierror_flag)='Error in HFBTHO_ITERATIONS read'
  !     Return
  !  End If
  !  !
  !  ! Type of functional
  !  Read(UNIT=lnamelist,NML=HFBTHO_FUNCTIONAL,iostat=ios)
  !  If (ios.Ne.0) Then
  !     ierror_flag=ierror_flag+1
  !     ierror_info(ierror_flag)='Error in HFBTHO_FUNCTIONAL read'
  !     Return
  !  End If
  !  !
  !  ! Characteristics of pairing
  !  Read(UNIT=lnamelist,NML=HFBTHO_PAIRING,iostat=ios)
  !  If (ios.Ne.0) Then
  !     ierror_flag=ierror_flag+1
  !     ierror_info(ierror_flag)='Error in HFBTHO_PAIRING read'
  !     Return
  !  End If
  !  !
  !  ! Constraints
  !  Read(UNIT=lnamelist,NML=HFBTHO_CONSTRAINTS,iostat=ios)
  !  If (ios.Ne.0) Then
  !     ierror_flag=ierror_flag+1
  !     ierror_info(ierror_flag)='Error in HFBTHO_CONSTRAINTS read'
  !     Return
  !  End If
  !  !
  !  ! Blocking
  !  Read(UNIT=lnamelist,NML=HFBTHO_BLOCKING,iostat=ios)
  !  If (ios.Ne.0) Then
  !     ierror_flag=ierror_flag+1
  !     ierror_info(ierror_flag)='Error in HFBTHO_BLOCKING read'
  !     Return
  !  End If
  !  !
  !  ! Particle number projection
  !  Read(UNIT=lnamelist,NML=HFBTHO_PROJECTION,iostat=ios)
  !  If (ios.Ne.0) Then
  !     ierror_flag=ierror_flag+1
  !     ierror_info(ierror_flag)='Error in HFBTHO_PROJECTION read'
  !     Return
  !  End If
  !  !
  !  ! Finite temperature
  !  Read(UNIT=lnamelist,NML=HFBTHO_TEMPERATURE,iostat=ios)
  !  If (ios.Ne.0) Then
  !     ierror_flag=ierror_flag+1
  !     ierror_info(ierror_flag)='Error in HFBTHO_TEMPERATURE read'
  !     Return
  !  End If
  !  !
  !  ! debug
  !  Read(UNIT=lnamelist,NML=HFBTHO_DEBUG,iostat=ios)
  !  If (ios.Ne.0) Then
  !     ierror_flag=ierror_flag+1
  !     ierror_info(ierror_flag)='Error in HFBTHO_DEBUG read'
  !     Return
  !  End If
  !  !
  !  Close(lnamelist)
  !  !
  !End Subroutine read_HFBTHO_NAMELIST
  !=======================================================================
  Subroutine check_consistency
    !---------------------------------------------------------------------
    ! Check consistency of input data
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: counter, i
    Real(pr) :: A, preset_inin(3)
    Character(30), Dimension(:) :: preset_forces(18)
    !
    If((n00_INI.Lt.1).Or.(n00_INI.GT.50)) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("number_of_shells = ",i6," out-of-bounds: [1,50]")') &
             n00_INI
       Return
    End If
    !
    If((npr_INI(1).Lt.1).Or.(npr_INI(2).Lt.1)) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("Z = ",i6," N = ",i6," out-of-bounds: (Z,N)>1")') &
             npr_INI(2),npr_INI(1)
       Return
    End If
    !
    If(Abs(kindhfb_INI).Ne.1) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("type_of_calculation = ",i6," unrecognized: (-1,1)")') &
             kindhfb_INI
       Return
    End If
    !
    If(epsi_INI.Lt.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("accuracy = ",e24.12," out-of-bounds: >0")') &
             epsi_INI
       Return
    End If
    !
    preset_inin( 1) = 1
    preset_inin( 2) = 2
    preset_inin( 3) = 3
    !
    counter=0
    Do i=1, 3
       If(Abs(inin_INI).Eq.preset_inin(i)) Then
          counter=1
          Exit
       End If
    End Do
    !
    If(counter.Eq.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("restart_file = ",i6," unrecognized: see list in publi")') &
             inin_INI
       Return
    End If
    !
    preset_forces( 1) = 'SIII'
    preset_forces( 2) = 'SKM*'
    preset_forces( 3) = 'SKP'
    preset_forces( 4) = 'SLY4'
    preset_forces( 5) = 'SLY5'
    preset_forces( 6) = 'SLY6'
    preset_forces( 7) = 'SLY7'
    preset_forces( 8) = 'SKI3'
    preset_forces( 9) = 'SKO'
    preset_forces(10) = 'SKX'
    preset_forces(11) = 'HFB9'
    preset_forces(12) = 'UNE0'
    preset_forces(13) = 'UNE1'
    preset_forces(14) = 'UNE2'
    preset_forces(15) = 'HFB1'
    preset_forces(16) = 'N0LO'
    preset_forces(17) = 'N1LO'
    preset_forces(18) = 'N2LO'
    !
    counter=0
    Do i=1, 18
       If(Trim(skyrme_INI).Eq.Trim(preset_forces(i))) Then
          counter=1
          Exit
       End If
    End Do
    ! Functional must be in preset list
    If(counter.Eq.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("functional = ",a30," unrecognized: see list in publi")') &
             skyrme_INI
       Return
    End If
    ! Pairing cut-off must be positive
    If(pwi_INI.Lt.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("pairing_cutoff = ",i4," out-of-bounds: >=0")') &
             pwi_INI
       Return
    End If
    ! Pairing cut-off must be positive
    If(cpv1_INI.Lt.0.0.Or.cpv1_INI.Gt.1.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("pairing_feature = ",i4," out-of-bounds: [0.0,1.0]")') &
             cpv1_INI
       Return
    End If
    ! Options for Coulomb: 0, 1, 2
    If(icou_INI.Lt.0.Or.icou_INI.Gt.2) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("type_of_coulomb = ",i4," unrecognized: (0,1,2)")') &
             icou_INI
       Return
    End If
    ! Choices of basis (HO or THO): -1, 0, 1
    If(Abs(iLST_INI).Gt.1) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("switch_to_THO = ",i4," unrecognized: (-1,0,1)")') &
             iLST_INI
       Return
    End If
    ! At least one gauge point if projection is required
    If(keypj_INI.Le.0.And.iproj_INI.Ne.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("gauge_points = ",i4," out-of-bounds: >=0")') &
             keypj_INI
       Return
    End If
    ! Number of protons must be greater than 0 for projection
    If((npr_INI(1)+npr1pj_INI).Lt.1.And.iproj_INI.Ne.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("delta_N = ",i4," out-of-bounds: N+dN>=1")') &
             npr1pj_INI
       Return
    End If
    ! Number of neutrons must be greater than 0 for projection
    If((npr_INI(2)+npr2pj_INI).Lt.1.And.iproj_INI.Ne.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("delta_Z = ",i4," out-of-bounds: Z+dZ>=1")') &
             npr2pj_INI
       Return
    End If
    ! Temperature must be positive
    If(temper.Lt.zero) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("temperature = ",i4," out-of-bounds: T>=0")') &
             temper
       Return
    End If
    ! Number of Gauss-Laguerre integration points between 0 and 100
    If(ngh_INI.Lt.1.Or.ngh_INI.Gt.100) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("number_Gauss = ",i4," out-of-bounds: [1,100]")') &
             ngh_INI
       Return
    End If
    ! Number of Gauss-Hermite integration points between 0 and 100
    If(ngl_INI.Lt.1.Or.ngl_INI.Gt.100) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("number_Laguerre = ",i4," out-of-bounds: [1,100]")') &
             ngl_INI
       Return
    End If
    ! Number of Gauss-Legendre integration points lower than 100
    If(nleg_INI.Gt.100) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("number_Legendre = ",i4," out-of-bounds: [-infty,100]")') &
             nleg_INI
       Return
    End If
    ! Number of Gauss-Legendre integration points between 1 and 100 for PNP
    If((nleg_INI.Lt.1.Or.nleg_INI.Gt.100).And.iproj_INI.Ne.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("number_Legendre = ",i4," out-of-bounds: [1,100]")') &
             nleg_INI
       Return
    End If
    ! Number of basis states must be greater than 0
    If(nstate_INI.Lt.1.And.basis_HFODD_INI) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("number_states = ",i4," out-of-bounds: >0")') &
             nstate_INI
       Return
    End If
    !
  End Subroutine check_consistency
  !=======================================================================
  !
  !=======================================================================
  Subroutine initialize_HFBTHO_SOLVER
    !---------------------------------------------------------------------
    ! default parameters
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Real(pr) :: A
    !------------------------------------
    ! tapes
    !------------------------------------
    lwin=41; lwou=42;  lwel=52; lres=57; lin=3
    !------------------------------------
    ! From Namelist or default values
    !------------------------------------
    nstate                = nstate_INI
    epsi                  = epsi_INI                  ! stop criteria
    Add_Pairing           = Add_Pairing_INI           ! add pairing starting from file
    icou                  = icou_INI                  ! coul: no-(0), dir.only-(1), plus exchange-(2)
    DO_FITT               = DO_FITT_INI               ! calculates quantities for reg.optimization
    IDEBUG                = IDEBUG_INI                ! debug
    Parity                = Parity_INI                ! reflection symmetry
    Print_HFBTHO_Namelist = Print_HFBTHO_Namelist_INI ! Print Namelist
    !---------------------------------------------------------------------
    ! Pairing set by user
    !---------------------------------------------------------------------
    If(set_pairing) Then
       CpV0(0)=V0n_INI
       CpV0(1)=V0p_INI
       CpV1(0)=cpv1_INI
       CpV1(1)=cpv1_INI
       pwi=pwi_INI
    End If
    !------------------------------------
    ! output control
    !------------------------------------
    !If(n00_INI.Gt.0) Then
    !   Print_Screen=.True.
    !   lfile=lout+1                                 ! output to screen & thoout.dat
    !Else
       Print_Screen=.False.
       lfile=lout                                 ! no output to screen & thoout.dat
    !End If
    !------------------------------------
    ! Pi
    !------------------------------------
    PI=four*Atan(one)
    !------------------------------------
    ! blocking
    !------------------------------------
    bloblo=0; blo123=0; blok1k2=0;  keyblo=0
    blomax=0; nkblo=0;  iparenti=0; irestart=0
    blocanon=0;         eqpmin=zero
    !------------------------------------
    ! buffers
    !------------------------------------
    eres=zero;  eresu=zero;  eresl=zero;
    eresj=zero; eresbl=zero; ereslbl=' 00[00,00,00]'
    !------------------------------------
    ! def parameters
    !------------------------------------
    ffdef3=Sqrt(five/(four*pi))/two
    ffdef4=Sqrt(117.0_pr)/(four*pi)
    ffdef5=Sqrt(nine/(four*pi))/eight
    ffdef6=Sqrt(five*pi)/three
    ffdef7=Sqrt(pi)/four
    !------------------------------------
    ! former linear mixing
    !------------------------------------
    xmix0=0.3_pr             ! lowest mixing parameter  (redefined later)
    xmix =0.3_pr             ! initial mixing parameter (changes every iteration)
    xmax =1.0_pr             ! mario
    !------------------------------------
    ! misc (redefined later)
    !------------------------------------
    rehfbcan=0.0_pr; depnp=0.0_pr; ala2=0.00_pr
    ept=-2.0_pr; del=1.0_pr; ala=-7.0_pr; bet=0.0_pr
    ala1(1)=-14.6851; ala1(2)=-3.7522; si=1.0_pr
    iqrpa=0; icacou=0;  icacoupj=0; icahartree=0; iasswrong=0
    iError_in_HO=0;  iError_in_THO=0
    ECMHFB=0.0_pr; ECMPAV=0.0_pr
    If(use_full_cm_cor) Then
       !A = npr_INI(1) + npr_INI(2)
       A = tz(1) + tz(2)
       facECM = A/(A-1.0_pr)
    End If
    entropy(:)=zero
    !------------------------------------
    ! Saxon-Woods: von koepf und ring, z.phys. (1991)
    !------------------------------------
    v0ws=-71.28_pr; akv=0.4616_pr; r0v=1.2334_pr; av=0.6150_pr
    vso=11.1175_pr; rso=1.1443_pr; aso=0.6476_pr
    !------------------------------------
    ! fixed text
    !------------------------------------
    tp(1)='+'; tp(2)='-'; tis(1)='n'; tis(2)='p';
    tit(1)='neutrons'; tit(2)='protons '
    tl(0)='s'; tl(1)='p'; tl(2)='d'; tl(3)='f'; tl(4)='g'
    tl(5)='h'; tl(6)='i'; tl(7)='j'; tl(8)='k'; tl(9)='l'
    tl(10)='m'; tl(11)='n'; tl(12)='o'; tl(13)='p'; tl(14)='q'
    tl(15)='r'; tl(16)='s'; tl(17)='t'; tl(18)='u'; tl(19)='v'; tl(20)='w'
    !------------------------------------
    ! fixed parity sign
    !------------------------------------
    tpar(1)=+1; tpar(2)=-1;
    !------------------------------------
    ! physical constants
    !------------------------------------
    amn=938.90590_pr
    amu=931.4940130_pr; r0=1.20_pr
    alphi=137.036020_pr; hqc=197.328910_pr
    !------------------------------------
    ! e2 for protons (set now in elsewhere)
    !------------------------------------
    !chargee2=hqc/alphi
    !chargee2=1.43997841_pr
    !-----------------------------------
    ! set the loops over particle types
    !-----------------------------------
    itmin=1 ; itmax = 2;
    If(npr_INI(1).Eq.0) itmin = 2
    If(npr_INI(2).Eq.0) itmax = 1
    !-----------------------------------
    ! error flag and info
    !-----------------------------------
    ierror_flag=0
    ierror_info(ierror_flag)='No errors in the solver!'
    !
    Call set_functional_parameters(skyrme_INI,.False.)
    !-----------------------------------
    ! set multipole moments units
    !-----------------------------------
    Call moments_setUnits
    !
  End Subroutine initialize_HFBTHO_SOLVER
  !=======================================================================
  !
  !=======================================================================
  Subroutine base0(lpr)
    !---------------------------------------------------------------------
    ! selects HO basis configurations in cylindrical coordinates
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Logical :: lpr
    Integer(ipr) :: iw,k,nre,nze,ke,la,le,ip,ir,iz,il,is,Iall,ilauf,jlauf,ib,nd
    Integer(ipr) :: NOSCIL
    Real(pr), Allocatable :: e(:)
    Real(pr) :: hbz,hbp,ee
    !
    If(n00.Gt.n00max) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='STOP: too large n00 versus n00max'
       Return
    End If
    !-----------------------------------------------
    ! MAXIMUM NUMBER OF THE HO SHELLS (n00,NOSCIL)
    ! (7,120),(8,165),(9,220),(10,286),(11,364)
    ! (12,455),(14,680),(16,969),(18,1330),(20,1771)
    !-----------------------------------------------
    NOSCIL=(n00+1)*(n00+2)*(n00+3)/6
    !-----------------------------------------------
    ! count all states for n00max
    !-----------------------------------------------
    nze=n00max; nre=n00max/2; ke=n00max
    If(basis_HFODD) Then
       nze=n00; nre=n00/2; ke=n00
    End If
    Iall=0;
    Do k=1,ke+1
       la=k-1; le=min0(ke,k)
       Do ip=1,2
          Do ir=0,nre
             Do iz=0,nze
                Do il=la,le
                   Do is=+1,-1,-2
                      If (iz+2*ir+il.Gt.n00max)    Cycle
                      If (il+(is+1)/2.Ne.k)        Cycle
                      If (Mod(iz+il,2).Ne.ip-1)    Cycle
                      Iall=Iall+1
                   End Do
                End Do
             End Do
          End Do
       End Do
    End Do
    !-----------------------------------------------
    ! charge all energies for n00max
    !-----------------------------------------------
    Allocate(e(Iall))
    hbz=two*hbzero/bz**2; hbp=two*hbzero/bp**2;
    Iall=0;
    Do k=1,ke+1
       la=k-1; le=min0(ke,k)
       Do ip=1,2
          Do ir=0,nre
             Do iz=0,nze
                Do il=la,le
                   Do is=+1,-1,-2
                      If (iz+2*ir+il.Gt.n00max)    Cycle
                      If (il+(is+1)/2.Ne.k)        Cycle
                      If (Mod(iz+il,2).Ne.ip-1)    Cycle
                      Iall=Iall+1
                      e(Iall)=hbz*(Real(iz,Kind=pr)+half) &
                         +hbp*(two*Real(ir,Kind=pr)+Real(il,Kind=pr)+one)
                   End Do
                End Do
             End Do
          End Do
       End Do
    End Do
    !-----------------------------------------------
    ! sort energies and derive base cut-off energy
    !-----------------------------------------------
    Call ord(Iall,e);
    If(Iall.Gt.NOSCIL) Then
       EBASECUT=E(NOSCIL)+1.0D-5
    Else
       EBASECUT=E(Iall)+1.0D-5
    End If
    If(basis_HFODD) EBASECUT=E(nstate)+1.0D-5
    Deallocate(e)
    !-----------------------------------------------
    ! calculate the actual states
    !-----------------------------------------------
    nze=n00max; nre=n00max/2; ke=n00max
    If(basis_HFODD) Then
       nze=n00; nre=n00/2; ke=n00
    EndIf
    ib=0; ilauf=0; ndx=0; nzx=0; nrx=0; nlx=0; nqp=0; nuv=0
    ! loop over k-quantum number
    Do k=1,ke+1
       la=k-1; le=min0(ke,k)
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
                      If (Mod(iz+il,2).Ne.(ip-1))  Cycle
                      ee=hbz*(Real(iz,Kind=pr)+half)&
                    +hbp*(two*Real(ir,Kind=pr)+Real(il,Kind=pr)+one)
                      If(ee.Lt.EBASECUT) Then
                         ilauf=ilauf+1
                         nzx=Max(nzx,iz); nrx=Max(nrx,ir); nlx=Max(nlx,il)
                      End If
                   End Do
                End Do
             End Do
          End Do
          If(Parity) Then                !Yesp
             If (ilauf.Gt.jlauf) Then
                ib=ib+1
                nd=ilauf-jlauf
                ndx=Max(ndx,nd)
                nqp=nqp+nd; nuv=nuv+nd*nd
             End If
          End If
       End Do
       If(.Not.Parity) Then              !Nop
          If(ilauf.Gt.jlauf) Then
             ib=ib+1
             nd=ilauf-jlauf
             ndx=Max(ndx,nd)
             nqp=nqp+nd; nuv=nuv+nd*nd
          End If
       End If
    End Do
    nbx=ib; ntx=ilauf
    !-----------------------------------------------
    ! print statistics
    !-----------------------------------------------
    If(lpr) Then
       Do iw=lout,lfile
          Write(iw,*)
          Write(iw,'(a)')  '  ---------------------------------------'
          Write(iw,'(a)')  '        Harmonic Oscillator Basis        '
          Write(iw,'(a)')  '  ---------------------------------------'
          Write(iw,'(a,2(i6,2x),a)') '  NUV, NQP:                      ',nuv,nqp
          Write(iw,'(a,2(i6,2x),a)') '  Comparison with bookkeeping spherical basis:'
          Write(iw,'(a,2(i6,2x),a)') '  n00:                           ',n00,n00,  &
             'Maximal number of shells'
          Write(iw,'(a,2(i6,2x),a)') '  nbx, 2*n00+1:                  ',nbx,2*n00+1,  &
             'Maximal number of K-blocks'
          Write(iw,'(a,2(i6,2x),a)') '  ntx, (n00+1)*(n00+2)*(n00+3)/6 ',ntx,(n00+1)*(n00+2)*(n00+3)/6,  &
             'Max.num. p/n levels'
          Write(iw,'(a,2(i6,2x),a)') '  nzx, n00:                      ',nzx,n00,  &
             'Maximal nz-quantum number'
          Write(iw,'(a,2(i6,2x),a)') '  nrx, n00/2  :                  ',nrx,n00/2,  &
             'Maximal nr-quantum number'
          Write(iw,'(a,2(i6,2x),a)') '  nlx, n00:                      ',nlx,n00,  &
             'Maximal ml-quantum number'
          Write(iw,'(a,2(i6,2x),a)') '  ndx, (n00+2)*(n00+2)/4:        ',ndx,(n00+2)*(n00+2)/4,  &
             'Maximal dim. of one k-block'
          Write(iw,*)
       End Do
    End If
    !
  End Subroutine base0
  !=======================================================================
  !
  !=======================================================================
  Subroutine base(lpr)
    !---------------------------------------------------------------------
    ! set HO basis configurations in cylindrical coordinates
    !---------------------------------------------------------------------
    Use HFBTHO
    Implicit None
    Logical :: lpr
    Integer(ipr) :: nze,nre,ke,ib,ilauf,jlauf,nom,nnm,  &
         k,la,le,ip,ir,iz,il,is,nn,ND,IBX,N1,N2,iw
    Real(pr) :: hbz,hbp,ee
    !
    hbz=two*hbzero/bz**2; hbp=two*hbzero/bp**2;
    !
    nze=n00max; nre=n00max/2; ke=n00max
    If(basis_HFODD) Then
       nze=n00; nre=n00/2; ke=n00
    End If
    ib=0; ilauf=0; nzm=0; nrm=0; nlm=0; nom=0; nnm=0
    !-----------------------------------------------
    ! loop over k-quantum number
    !-----------------------------------------------
    Do k=1,ke+1
       la=k-1; le=min0(ke,k)
       ! loop over parity
       If(.Not.Parity) jlauf=ilauf !Nop
       Do ip=1,2
          If(Parity) jlauf=ilauf   !Yesp
          Do ir=0,nre
             Do iz=0,nze
                Do il=la,le
                   Do is=+1,-1,-2
                      If (iz+2*ir+il.Gt.n00max) Cycle
                      If (il+(is+1)/2.Ne.k)     Cycle
                      If (Mod(iz+il,2).Ne.ip-1) Cycle
                      ee=hbz*(Real(iz,Kind=pr)+half)&
                    +hbp*(two*Real(ir,Kind=pr)+Real(il,Kind=pr)+one)
                      If(ee.Lt.EBASECUT) Then
                         ilauf=ilauf+1
                         If (ilauf.Gt.ntx) Then
                            ierror_flag=ierror_flag+1
                            ierror_info(ierror_flag)='STOP: in base: ntx too small'
                            Return
                         End If
                         nz(ilauf)=iz; nr(ilauf)=ir; nl(ilauf)=il; ns(ilauf)=is; npar(ilauf)=ip
                         nn =iz+2*ir+il
                         Write(tb(ilauf),100) 2*k-1,tp(ip),nn,iz,il
100                      Format(i2,a1,'[',i2,',',i2,',',i2,']')
                         Do iw=lout,lfile
                          If(lpr.And.IDEBUG.Gt.10) &
                             Write(iw,'(i4,a,i2,a,i2,a,i2,a,i2,a,i2,a,2x,a,1x,a,f14.8)')  &
                                 ilauf,'   nn=',nn,'   nz=',iz,'   nr=',ir,  &
                                 '   ml=',il,'  ms=',is,' /2',tb(ilauf),'e=',ee
                         End Do
                         nzm=Max(nzm,iz); nrm=Max(nrm,ir); nlm=Max(nlm,il)
                         nom=Max(nom,2*k-1); nnm=Max(nnm,iz+2*ir+il)
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
                ia(ib)=jlauf; id(ib)=ilauf-jlauf
                ikb(ib)=k; ipb(ib)=ip
                Write(txb(ib),'(i3,a,i2,a,a1)') ib,'. block:  k=',k+k-1,'/2',tp(ip)
                !ir=(ib+1)/2
                !Write(*,*)  ib,2*k-1,'2*Omega=',2*ir - 1
                Do iw=lout,lfile
                   If(lpr.And.IDEBUG.Gt.10) Write(iw,'(/,a,i3,a,a1)')'  For the above block:  k=',k+k-1,'/2',tp(ip)
                End Do
             End If
             If(id(ib).Eq.0) Then
                ierror_flag=ierror_flag+1
                ierror_info(ierror_flag)='STOP: in base Block Memory(1)'
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
             ia(ib)=jlauf; id(ib)=ilauf-jlauf
             nn = nz(ilauf)+2*nr(ilauf)+nl(ilauf); ip = 2 - Mod(nn,2)
             ikb(ib)=k; ipb(ib)=ip
             Write(txb(ib),'(i3,a,i2,a,a1)') ib,'. block:  k=',k+k-1,'/2',tp(ip)
             Do iw=lout,lfile
                If(lpr.And.IDEBUG.Gt.10) Write(iw,'(/,a,i3,a,a1)')'  For the above block:  k=',k+k-1,'/2',tp(ip)
             End Do
          End If
          If(id(ib).Eq.0) Then
             ierror_flag=ierror_flag+1
             ierror_info(ierror_flag)='STOP: in base Block Memory(2)'
             Return
          End If
       End If
    End Do ! end k
    nb=ib;  nt=ilauf
    !-----------------------------------------------
    ! broyden/linear mixing (storage)
    !-----------------------------------------------
    nhhdim=0
    Do ib=1,NB
       ND=ID(ib)
       Do N1=1,ND
          Do N2=1,N1
             nhhdim=nhhdim+1
          End Do
       End Do
    End Do
    nhhdim2=2*nhhdim; nhhdim3=3*nhhdim; nhhdim4=4*nhhdim
    If(Allocated(brin)) Deallocate(brin,brout)
    Allocate(brin(nhhdim4+lambdaMax),brout(nhhdim4+lambdaMax))
    !-----------------------------------------------
    ! Print statistics
    !-----------------------------------------------
    If(lpr) Then
       Do iw=lout,lfile
          Write(iw,'(a,i4)')   '  Actual basis used'
          Write(iw,'(a,i4)')   '  Number of blocks: nb .......: ',nb
          Write(iw,'(a,i4)')   '  Number of levels: nt .......: ',nt
          Write(iw,'(a,i4)')   '  Maximal 2*omega : nom ......: ',nom
          Write(iw,'(a,i4)')   '  Maximal nz:       nzm ......: ',nzm
          Write(iw,'(a,i4)')   '  Maximal nr:       nrm ......: ',nrm
          Write(iw,'(a,i4)')   '  Maximal ml:       nlm ......: ',nlm
          Write(iw,'(a,i4)')   '  Maximal N=nz+2*nr+nl .......: ',nnm
          Write(iw,'(a,i4)')   '  2 x biggest block dim. .....: ',ndx2
          Write(iw,'(a,i8)')   '  Non-zero elements of h .....: ',nhhdim
          Write(iw,'(a,i8)')   '  Number of Broyden elements .: ',nhhdim4
          Write(iw,'(a,i4)')
       End Do
    End If
    If(nzm.Ge.n00max.Or.(nom-1)/2.Eq.n00max) Then
       Write(*,*) 'nzm=',nzm,'  (nom-1)/2=',(nom-1)/2,'  n00max=',n00max
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='STOP: Please increase n00max to have correct basis'
    End If
  End Subroutine base
  !=======================================================================
  !
  !=======================================================================
  Subroutine gaupol(lpr)
    !---------------------------------------------------------------------
    ! HO wave functions in cylindrical coordinates
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Logical :: lpr
    Real(pr) :: w0,z,x,s,s0,s1,w00,w4pii,dsq,d1,d2,d3,d4,hs0,hs1
    Integer(ipr) :: ih,il,iw,ix,n,l,n1,n2
    !-----------------------------------------------
    ! Hermite
    !-----------------------------------------------
    w4pii=pi**(-0.250_pr)
    Do ih=1,ngh
       z=xh(ih); w0=w4pii*Exp(-half*z*z)
       ! functions qh, qh1; norm: \sum_{ih} qh(n1,ih)*qh(n2,ih)=\delta_{n1,n2}
       w0 =w0*Sqrt(wh(ih))
       qh(0,ih)=w0;       qh(1,ih)=sq(2)*w0*z
       qh1(0,ih)=-w0*z;   qh1(1,ih)=sq(2)*w0*(one-z*z)
       Do n=2,nzm
          qh(n,ih)=sqi(n)*(sq(2)*z*qh(n-1,ih)-sq(n-1)*qh(n-2,ih))
          qh1(n,ih)=sq(n+n)*qh(n-1,ih)-z*qh(n,ih)
       End Do
    End Do
    !-----------------------------------------------
    ! Laguerre
    !-----------------------------------------------
    Do il=1,ngl
       x=xl(il); w00=sq(2)*Exp(-half*x)
       Do l=0,nlm
          ! functions ql, ql1; norm: \sum_{il} ql(n1,l,il)*ql(n2,l,il)=\delta_{n1,n2}
          w0=w00*Sqrt(half*wl(il)*x**l)
          ql(0,l,il)=wfi(l)*w0;         ql(1,l,il)=(l+1-x)*wfi(l+1)*w0
          ql1(0,l,il)=(l-x)*wfi(l)*w0;  ql1(1,l,il)=(Real(l*l+l,Kind=pr) &
                                                  -x*Real(l+l+3,Kind=pr)+x*x)*wfi(l+1)*w0
          Do n=2,nrm
             dsq=sq(n)*sq(n+l); d1=Real(n+n+l-1,Kind=pr)-x
             d2=sq(n-1)*sq(n-1+l); d3=n+n+l-x; d4=two*dsq
             ql(n,l,il)=(d1*ql(n-1,l,il)-d2*ql(n-2,l,il))/dsq
             ql1(n,l,il)=d3*ql(n,l,il)-d4*ql(n-1,l,il)
          End Do
       End Do
    End Do
    !-----------------------------------------------
    ! Test accuracy for Hermite orthonormalization
    !-----------------------------------------------
    hs0=zero; hs1=two
    Do n1=0,nzm
       Do n2=0,n1
          If (Mod(n1-n2,2).Eq.0) Then
             s=zero
             Do ih=1,ngh
                s=s+qh(n1,ih)*qh(n2,ih)
             End Do
             If(n1.Ne.n2) Then
                hs0=Max(s,hs0)
             Else
                hs1=Min(s,hs1)
             End If
          End If
       End Do
    End Do
    !-----------------------------------------------
    ! Test accuracy for Laguerre orthonormalization
    !-----------------------------------------------
    s0=zero; s1=two
    Do l=0,nlm
       Do n1=0,nrm
          Do n2=0,n1
             s=zero
             Do il=1,ngl
                s=s+ql(n1,l,il)*ql(n2,l,il)
             End Do
             If(n1.Ne.n2) Then
                s0=Max(s,s0)
             Else
                s1=Min(s,s1)
             End If
          End Do
       End Do
    End Do
    !-----------------------------------------------
    ! print accuracy
    !-----------------------------------------------
    If(lpr) Then
       Do iw=lout,lfile
          Write(iw,'(a)')  '  ---------------------------------------'
          Write(iw,'(a)')  '            Integration Meshes           '
          Write(iw,'(a)')  '  ---------------------------------------'
          Write(iw,'(a,i3)')  &
               '  Number of Gauss-Hermite mesh points ngh ....: ',ngh
          Write(iw,'(a,i3)')  &
               '  Number of Gauss-Laguerre mesh points ngl ...: ',ngl
          Write(iw,'(a,i3)')  &
               '  Number of Gauss-Legendre mesh points nleg ..: ',nleg
          Write(iw,'(a)') &
               '  Integration boundaries'
          Write(iw,'(2(a,f12.8))')  &
               '    Hermite  - from xh(1)  =',xh(1),  ' to xh(ngh)   =',xh(ngh)
          Write(iw,'(2(a,f12.8))')  &
               '    Laguerre - From xl(1)  =',xl(1),  ' to xl(ngl)   =',xl(ngl)
          If(nleg.Gt.0) Then
            Write(iw,'(2(a,f12.8))')  &
                 '    Legendre - From xleg(1)=',xleg(1),' to xleg(nleg)=',xleg(nleg)
          End If
          Write(iw,*)  &
               ' Max.dev.in:     Orthogonality            Normalization'
          Write(iw,*) ' Hermite  ',hs0,Abs(one-hs1)
          Write(iw,*) ' Laguerre ',s0,Abs(one-s1)
       End Do
    End If
    !-----------------------------------------------
    ! debug
    !-----------------------------------------------
    If (lpr.And.IDEBUG.Gt.20) Then
       ix=3
       Do iw=lout,lfile
          Write(iw,*) ' nz    qh(nz,ih=1,...)'
          Do n=0,nzm
             Write(iw,'(i4,3f15.8)') n,(qh(n,ih),ih=1,ix)
             Write(iw,'(i4,3f15.8)') n,(qh1(n,ih),ih=1,ix)
             Write(iw,*) ' '
          End Do
          Do l=0,nlm
             Write(iw,*) ' nr ml    ql(nr,l,il=1,...)'
             Do n=0,nrm
                Write(iw,'(i4,i3,3f15.8)') n,l,(ql(n,l,il),il=1,ix)
                Write(iw,'(i4,i3,3f15.8)') n,l,(ql1(n,l,il),il=1,ix)
                Write(iw,*) ' '
             End Do
          End Do
       End Do
       !-----------------------------------------------
       ! Test for Hermite polynomials normalization
       !-----------------------------------------------
       Do n1=0,nzm
          Do n2=0,n1
             If (Mod(n1-n2,2).Eq.0) Then
                s=zero
                Do ih=1,ngh
                   s=s+qh(n1,ih)*qh(n2,ih)
                End Do
                Do iw=lout,lfile
                   Write(iw,100) n1,n2,s
                End Do
100             Format(' Gauss-Hermite: n1=',i3,'  n2=',i3,f20.8)
             End If
          End Do
       End Do
       !-----------------------------------------------
       ! Test for Laguerre polynomials normalization
       !-----------------------------------------------
       Do l=0,nlm
          Do n1=0,nrm
             Do n2=0,n1
                s=zero
                Do il=1,ngl
                   s=s+ql(n1,l,il)*ql(n2,l,il)
                End Do
                Do iw=lout,lfile
                   Write(iw,101) l,n1,n2,s
101                Format(' Gauss Laguerre: l='  &
                        ,i2,' n1=',i3,'  n2=',i3,f20.8)
                End Do
             End Do
          End Do
       End Do
    End If
    !
    Call coordinateLST(.False.)  ! coordinate LST
    !
  End Subroutine gaupol
  !=======================================================================
  !
  !=======================================================================
  Subroutine FileLabels(NPRI,ININL,FILELABEL)
    !---------------------------------------------------------------------
    ! file labels, e.g., filelabel='s070_040'
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: it,ininabs,ininl,nprt(2),NPRI(2)
    Character(1)  :: sinin
    Character(3)  :: snpr(2)
    Character(8)  :: filelabel
    !
    ininabs=iabs(ininl)
    filelabel='        ';sinin='t';snpr(1)='020'; snpr(2)='020'
    If(ininabs.Eq.4.or.ininabs.Eq.400) sinin='t'
    If(ininabs.Eq.3.or.ininabs.Eq.300) sinin='o'
    If(ininabs.Eq.2.or.ininabs.Eq.200) sinin='p'
    If(ininabs.Eq.1.or.ininabs.Eq.100) sinin='s'
    !
    nprt=npri
    Do it=itmin,itmax
       If(npri(it).Ne.2*(npri(it)/2)) nprt(it)=nprt(it)+iparenti(it)  !iparent=-/+ means particles/holes
       If(nprt(it).Lt.10  ) Then
          Write(snpr(it),'(a2,i1)') '00',nprt(it)
       Else
          If(nprt(it).Lt.100 ) Then
             Write(snpr(it),'(a1,i2)') '0',nprt(it)
          Else
             Write(snpr(it),'(i3)') nprt(it)
          End If
       End If
    End Do
    !
    Write(filelabel,'(a1,a3,a1,a3)')  sinin,snpr(1),'_',snpr(2)
    !
  End Subroutine FileLabels
  !=======================================================================
  !
  !=======================================================================
  Subroutine inout(is)
    !---------------------------------------------------------------------
    ! is=1: reads matrix elements from tape and exit
    ! is=2: writes matrix elements to tape and exit
    ! NB! if the welfile is missing or corrupt call start
    !     to restart calculations from scratch
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Logical       :: file_exists,file_opened
    Integer(ipr)  :: is,iw,n1,n2,nd,ib,bloall1,lambdaMax1,ierr,counterLine
    Character(8)  :: filelabel
    Character(30) :: welfile
    Real(pr)      :: tz1(2),b01,bz1,bp1,beta1,v0r(2),v1r(2),pwir
    Integer(ipr)  :: npr1,npr11,ngh1,ngl1,n001,nt1
    Integer(ipr)  :: ntx1,nb1,nhhdim1,NLANSA0,NLANSA1,NZA2NRA,NZA1,NLA1
    Integer(ipr)  :: ID1(nbx)
    !==== HFBODD interface
    Integer(ipr)  :: i,nza,nra,nla,nsa,ibasis
    Integer(ipr)  :: ibro
    !====
    ! label organization
    filelabel='        '
    Call FileLabels(NPR,ININ,filelabel)
    If(ierror_flag.Ne.0) Return
    If(iLST1.Le.0) Write(welfile,'(a8,a4)')  filelabel,'.hel'
    If(iLST1.Gt.0) Write(welfile,'(a8,a4)')  filelabel,'.tel'
    !
    welfile = welfile_INI
    If (is.Eq.1) Then
       !---------------------------------------------------------------------
       ! read matrix elements from 'welfile' file or start from scratch
       !---------------------------------------------------------------------
       If(inin.Gt.0) Then
          Call start
          Return
       End If
       !---------------------------------------------------------------------
       ! Check status of file on disk
       !---------------------------------------------------------------------
       file_exists=.False.; inquire(file=welfile, exist=file_exists); ierr=0
       If(file_exists) Then
          file_opened=.False.; inquire(unit=lwin, opened=file_opened)
          If(file_opened) Then
             Close(lwin)
          Else
          End If
          Open(lwin,file=welfile,status='old',form='unformatted',IOSTAT=ierr)
          If(ierr.NE.0) Then
             Do iw=lout,lfile
                Write(iw,'(1x,a,a,a)')
                Write(iw,'(1x,a,a,a)') ' The file ',welfile,' could not be opened!'
                Write(iw,'(1x,a,a,a)') ' STARTING FROM SCRATCH WITH ININ=IABS(ININ)!'
                Write(iw,'(1x,a,a,a)')
             End Do
             Call start
             Return
          End If
       Else
          Do iw=lout,lfile
             Write(iw,'(1x,a,a,a)')
             Write(iw,'(1x,a,a,a)') ' The file ',welfile,' is missing!'
             Write(iw,'(1x,a,a,a)') ' STARTING FROM SCRATCH WITH ININ=IABS(ININ)!'
             Write(iw,'(1x,a,a,a)')
          End Do
          Call start
          Return
       End If
       !---------------------------------------------------------------------
       ! Read data
       !---------------------------------------------------------------------
       counterLine = 0
       Read(lwin,ERR=100,End=100) npr11,npr1,ngh1,ngl1,n001,nb1,nt1
       counterLine = counterLine+1
       If(Abs(n001).Ne.Abs(n00).And.nb1.Ne.nb) go to 100
       Read(lwin,ERR=100,End=100) b01,bz1,bp1,beta1,si,etot,rms,bet,xmix,v0r,v1r,pwir, &
                                  del,ept,ala,ala2,alast,tz1,varmas,varmasNZ,pjmassNZ, &
                                  ass,skass
       brin=zero; bbroyden='L'; !si=one;
       counterLine = counterLine+1
       Read(lwin,ERR=100,End=100) ntx1,nb1,nhhdim1
       counterLine = counterLine+1
       Read(lwin,Err=100,End=100) lambdaMax1
       counterLine = counterLine+1
       Read(lwin,Err=100,End=100) multLag
       counterLine = counterLine+1
       Read(lwin,ERR=100,End=100) id1
       counterLine = counterLine+1
       Read(lwin,ERR=100,End=100) brin
       counterLine = counterLine+1
       !
       ! Add small pairing de=de+0.1 in the no-LN
       ! case to prevent pairing collaps
       If(kindhfb.Eq.1.And.Add_Pairing) Then
          ibro=0
          Do ib=1,NB
             ND=ID1(ib)
             I=ibro
             Do N1=1,ND
                Do N2=1,N1
                   I=I+1
                   brin(i+nhhdim2)=brin(i+nhhdim2)+0.10_pr
                   brin(i+nhhdim3)=brin(i+nhhdim3)+0.10_pr
                End Do !N2
             End Do !N1
             ibro=i
          End Do !IB
       End If
       Do ib=1,NB
          ND=ID1(ib)
          Do N1=1,ND
             Read(lwin,ERR=100,End=100) NLANSA0,NLANSA1,NZA2NRA,NZA1,NLA1
          End Do
       End Do
       counterLine = counterLine+1
       ! blocking
       Read(lwin,ERR=100,End=100)  bloall1
       counterLine = counterLine+1
       Read(lwin,ERR=100,End=100)  bloblo,blo123,blok1k2,blomax,bloqpdif
       counterLine = counterLine+1
       If(bloall1.Ne.bloall) go to 100
       !tel
       If(iLST.Gt.0) Then
          Read(lwin,ERR=100,End=100) decay,rmm3,cmm3,amm3,bmm3,itass,iqqmax
          If(Allocated(fdsx)) Deallocate(fdsx,fdsy,fdsy1,fdsy2,fdsy3,  &
               fspb0,fspc0,fspd0,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,  &
               fspb3,fspc3,fspd3)
          Allocate(fdsx(iqqmax),fdsy(iqqmax),fdsy1(iqqmax),  &
               fdsy2(iqqmax),fdsy3(iqqmax),fspb0(iqqmax),fspc0(iqqmax),  &
               fspd0(iqqmax),fspb1(iqqmax),fspc1(iqqmax),fspd1(iqqmax),  &
               fspb2(iqqmax),fspc2(iqqmax),fspd2(iqqmax),fspb3(iqqmax),  &
               fspc3(iqqmax),fspd3(iqqmax))
          Read(lwin,ERR=100,End=100) fdsx,fdsy,fdsy1,fdsy2,fdsy3,fspb0,fspc0,fspd0  &
               ,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,fspb3,fspc3,fspd3
       End If
       Do iw=lout,lfile
          Write(iw,*)
          Write(iw,*) ' Reading from wel_file: ',welfile
          Write(iw,*)
       End Do
       Close(lwin)
       Return
       !
100    Continue
       !---------------------------------------------------------------------
       ! missing or corrupt 'welfile' file
       !---------------------------------------------------------------------
       Close(lwin)
       Do iw=lout,lfile
          Write(iw,'(1x,a,a,a)')
          Write(iw,'(1x,a,a,a)')   ' The file ',welfile,' is corrupted!'
          Write(iw,'(1x,a,i2,a)')  ' Problem occurs at line ',counterLine,'        '
          Write(iw,'(1x,a,a,a)')   ' STARTING FROM SCRATCH WITH ININ=IABS(ININ)!'
          Write(iw,'(1x,a,a,a)')
       End Do
       Call start
       Return
    End If
    !---------------------------------------------------------------------
    ! write matrix elements to 'welfile' file
    !---------------------------------------------------------------------
    If (is.Eq.2.And.iasswrong(3).Eq.0) Then
       !---------------------------------------------------------------------
       ! Check status of file on disk
       !---------------------------------------------------------------------
       file_exists=.False.; inquire(file=welfile, exist=file_exists); ierr=0
       If(file_exists) Then
          file_opened=.False.; inquire(unit=lwou, opened=file_opened)
          If(file_opened) Then
             Close(lwou)
          Else
          End If
          Open(lwou,file=welfile,status='old',form='unformatted',IOSTAT=ierr)
          If(ierr.NE.0) Then
             Write(lout,'("Error in opening the old file, error code is ierr = ",i12)') ierr
             Return
          End If
       Else
          Open(lwou,file=welfile,status='new',form='unformatted',IOSTAT=ierr)
          If(ierr.NE.0) Then
             Write(lout,'("Error in opening the new file, error code is ierr = ",i12)') ierr
             Return
          End If
       End If
       !---------------------------------------------------------------------
       ! Write data
       !---------------------------------------------------------------------
       npr11=npr(1); npr1=npr(2)
       Write(lwou) npr11,npr1,ngh,ngl,n00,nb,nt
       Write(lwou) b0,bz,bp,beta0,si,etot,rms,bet,xmix,CpV0,CpV1,pwi,  &
                   del,ept,ala,ala2,alast,tz,varmas,varmasNZ,pjmassNZ, &
                   ass,skass
       Write(lwou) ntx,nb,nhhdim
       Write(lwou) lambdaMax
       Write(lwou) multLag
       Write(lwou) id
       Write(lwou) brin
       ibasis=0
       Do ib=1,NB
          ND=ID(ib)
          Do N1=1,ND
             ibasis=ibasis+1
             NLA=NL(ibasis); NRA=NR(ibasis); NZA=NZ(ibasis); NSA=NS(ibasis); NLANSA1=(-1)**(NZA+NLA)
             Write(lwou) 2*NLA+NSA,NLANSA1,NZA+2*NRA+NLA,NZA,NLA
          End Do
       End Do
       !---------------------------------------------------------------------
       ! blocking: sort blocking candidates first
       !---------------------------------------------------------------------
       Do ib=1,2
          Call blosort(ib,blomax(ib))
       End Do
       Write(lwou) bloall
       Write(lwou) bloblo,blo123,blok1k2,blomax,bloqpdif
       !tel
       If(iLST.Gt.0) Then
          If(Allocated(fdsx)) Then
             Write(lwou) decay,rmm3,cmm3,amm3,bmm3,itass,iqqmax
             Write(lwou) fdsx,fdsy,fdsy1,fdsy2,fdsy3,fspb0,fspc0,fspd0  &
                  ,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,fspb3,fspc3,fspd3
          End If
       End If
       Close(lwou)
       Do iw=lout,lfile
          Write(iw,'(a,a,a)')
          Write(iw,'(a,a,a)') '  Writing to wel_file: ',welfile
          Write(iw,'(a,a,a)') ' __________________________________  '
          Write(iw,'(a,a,a)') '  The tape ',welfile,' recorded:     '
          Write(iw,'(a,a,a)') '  nucname,npr,ngh,ngl,n00,nb,nt      '
          Write(iw,'(a,a,a)') '  b0,beta0,si,etot,rms,bet,xmix      '
          Write(iw,'(a,a,a)') '  pairing:     CpV0,CpV1,pwi         '
          Write(iw,'(a,a,a)') '  delta:       del,ept               '
          Write(iw,'(a,a,a)') '  lambda:      ala,ala2,alast,tz     '
          Write(iw,'(a,a,a)') '  asymptotic:  varmas,ass,skass      '
          Write(iw,'(a,a,a)') '  ntx,nb,nhhdim,id,N_rz,n_r,n_z      '
          Write(iw,'(a,a,a)') '  Omega2,Sigma2,Parity,Lambda        '
          Write(iw,'(a,a,a)') '  matrices(inbro):    hh,de          '
          Write(iw,'(a,a,a)') '  *all blocking candidates           '
          If(Allocated(fdsx)) Write(iw,'(a,a,a)') '  *all THO arrays                    '
          Write(iw,'(a,a,a)') ' __________________________________  '
          Write(iw,'(a,a,a)')
       End Do
    End If
    !
  End Subroutine inout
  !=======================================================================
  !
  !=======================================================================
  Subroutine start
    !---------------------------------------------------------------------
    ! initializes scratch Saxon-Woods potentials
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: iw,i,ih,il,ihl,it,ita
    Real(pr) :: zb(ngh),rrb(ngl),rb(ngl),rav,rao,vpws,vls,betas,gamma,fac, &
                facb,zz,rr,r,ctet,cphi,p2,p20,p22,s,u,w,f,rc,c,beta00,     &
                b2_ws,b4_ws,pleg2,pleg4
    !----------------------------------------------------------------------------
    ! Re-initializing all again since scratch calculation
    !----------------------------------------------------------------------------
    Call initialize_HFBTHO_SOLVER
    If(ierror_flag.Ne.0) Return
    Call Constraint_or_not(inin_INI,inin,icstr)
    If(ierror_flag.Ne.0) Return
    ! Setting restart to 1 for odd nuclei (blocking prescription)
    !Do it=itmin,itmax
    !   If(npr(it).Ne.2*(npr(it)/2)) Then
    !      irestart=irestart+1; npr(it)=npr_INI(it) ! Switch this off for HFODD
    !   End If
    !End Do
    !npr(3)=npr(1)+npr(2)
    amas=tz(1)+tz(2)
    If(irestart.Ne.0) Then
       ! odd nucleus requested but no even-even solution, recalculate the even-even nucleus from scratch
       Do iw=lout,lfile
          Write(iw,'(1x,a,2i4)')
          Write(iw,'(1x,a,2i4)') ' Initialization for the even-even core (N,Z)=: ',npr(1:2)
       End Do
    Else
       ! scratch for the even-even nucleus requested
       Do iw=lout,lfile
          Write(iw,'(1x,a,2i4)')
          Write(iw,'(a,a,3i4)')    '  Scratch initialization for the nucleus: ',nucname,npr(1:2)
          Write(iw,'(1x,a,2i4)')
       End Do
    End If
    n00=Abs(n00_INI);  b0=b0_INI;           q=q_INI; iLST=iLST_INI
    maxi=MAX_ITER_INI; inin=inin_INI;
    skyrme=skyrme_INI; kindhfb=kindhfb_INI
    iproj=iproj_INI;   npr1pj=npr1pj_INI;   npr2pj=npr2pj_INI;
    icacou=0; icahartree=0
    !
    Call preparer(.False.)
    !
    If(ierror_flag.Ne.0) Return
    inin=Abs(inin)          ! positive even if inin_INI is not
    If(Abs(b2_0).Gt.1.5_pr) b2_0=1.5_pr ! Avoid crazy initial points (quadrupole deformation)
    If(Abs(b4_0).Gt.1.0_pr) b4_0=1.0_pr ! Avoid crazy initial points (hexadecapole deformation)
    !-----------------------------------
    ! Saxon-Woods potentials
    !-----------------------------------
    Do iw=lout,lfile
       Write(iw,'(/,a)') '  Initial potentials of Saxon-Woods shape '
    End Do
    beta00=bet     ! wf to requested deformation
    Do iw=lout,lfile
       Write(iw,'(a,2f14.8)') '  v0ws   =',v0ws
       Write(iw,'(a,2f14.8)') '  kappa  =',akv
       Write(iw,'(a,2f14.8)') '  vs0    =',vso
       Write(iw,'(a,2f14.8)') '  r0     =',r0v
       Write(iw,'(a,2f14.8)') '  a      =',av
       Write(iw,'(a,2f14.8)') '  r0-so  =',rso
       Write(iw,'(a,2f14.8)') '  a-so   =',aso
       Write(iw,'(a,f14.8)')  '  b2_ws  =',b2_0
       Write(iw,'(a,f14.8)')  '  b4_ws  =',b4_0
    End Do
    !-----------------------------------
    ! Densities
    !-----------------------------------
    Do it=itmin,itmax
       ita=3-it; rav=r0v(it)*amas**p13; rao=rso(it)*amas**p13
       vpws=v0ws*(one-akv*(npr(it)-npr(ita))/amas)
       vls=half*(hqc/amu)**2*vpws*vso(it)
       ! Deformations of the surface
       b2_ws = b2_0 * Sqrt(5.0_pr/(4.0_pr*pi))
       b4_ws = b4_0 * Sqrt(9.0_pr/(4.0_pr*pi))
       ! Volume conservation condition
       !gamma=zero
       !fac= one+betas*Cos( gamma*pi/180.0_pr)
       !fac=(one+betas*Cos((gamma+120.0_pr)*pi/180.0_pr))*fac
       !fac=(one+betas*Cos((gamma-120.0_pr)*pi/180.0_pr))*fac
       !fac=fac**(-p13)
       fac = two + (143.0_pr*Sqrt(five)*b2_ws**3 + 1287.0_pr*b2_ws**2*b4_ws + 390.0_pr*Sqrt(five)*b2_ws*b4_ws**2 &
            + 243.0_pr*b4_ws**3)/(2002.0_pr*Pi**1.5_pr) + (three*(b2_ws**2 + b4_ws**2))/(two*Pi)
       fac=(two/fac)**(p13)
       ! z,r-coordinates in fm
       zb=xh*bz; rrb=xl*bp*bp; rb=Sqrt(rrb)
       Do ih=1,ngh
          zz=zb(ih)**2
          Do il=1,ngl
             rr=rrb(il)+zz; r=Sqrt(rr); ctet=zz/rr
             ! Deformed surface
             !p20=3.0_pr*ctet-one; p22=Sqrt(3.0_pr)*cphi
             !p2=p20*Cos(gamma*pi/180.0_pr)+p22*Sin(gamma*pi/180.0_pr)
             pleg2 = half*(three*ctet - one)
             pleg4 = (35.0_pr*ctet**2 - 30.0_pr*ctet + three)/eight
             facb=fac*(one + b2_ws*pleg2 + b4_ws*pleg4)
             ! Woods-Saxon potential
             u= vpws/( one+Exp( (r-rav*facb) / av(it) ))
             w=-vls /( one+Exp( (r-rao*facb) / aso(it)))
             ihl=ih+(il-1)*ngh
             If(it.Eq.1) Then
                vhbn(ihl)=hb0; vn(ihl)=u; vsn(ihl)=w;
                vrn(ihl)=zero; vzn(ihl)=zero; vdn(ihl)=zero;
                vSFIZn(ihl)=zero; vSZFIn(ihl)=zero;
                vSFIRn(ihl)=zero; vSRFIn(ihl)=zero;
             Else
                vhbp(ihl)=hb0; vp(ihl)=u; vsp(ihl)=w;
                vrp(ihl)=zero; vzp(ihl)=zero; vdp(ihl)=zero;
                vSFIZp(ihl)=zero; vSZFIp(ihl)=zero;
                vSFIRp(ihl)=zero; vSRFIp(ihl)=zero;
             End If
             ro(ihl,it)=u
             aka(ihl,it)=5.0d-3*Exp((r-rav*facb)/2.0_pr)
          End Do
       End Do
       !s=npr(it)/Sum(ro(:,it))
       s=tz(it)/Sum(ro(:,it))
       Do il=1,ngl
          Do ih=1,ngh
             ihl=ih+(il-1)*ngh
             f=s/(pi*wh(ih)*wl(il)* bz*bp*bp); ro(ihl,it)=f*ro(ihl,it)
          End Do
       End Do
       !-----------------------------------
       ! pairing
       !-----------------------------------
       Do il=1,nghl
          If(it.Eq.1) Then
             dvn(il)=-100.0_pr*aka(il,it)
          Else
             dvp(il)=-100.0_pr*aka(il,it)
          End If
       End Do
    End Do
    !-----------------------------------
    ! coulomb
    !-----------------------------------
    If(icou.Eq.0) Then
       cou=zero
    Else
       rc=r0v(2)*amas**p13
       Do il=1,ngl
          Do ih=1,ngh
             r=Sqrt(zb(ih)**2+rrb(il))
             If (r.Lt.rc) Then
                c=half*(3/rc-r*r/(rc**3))
             Else
                c=one/r
             End If
             !cou(ih+(il-1)*ngh)=c*npr(2)/alphi
             cou(ih+(il-1)*ngh)=c*tz(2)/alphi
          End Do
       End Do
    End If
    !-----------------------------------
    ! initial ph+pp matrix elements
    !-----------------------------------
    ak=0.1_pr; rk=0.1_pr ! initial density matrix elements (improve later)
    brin=zero             ! initial matrix elements to zero
    iiter=0               ! iteration number iiter to zero
    Call gamdel
    !
  End Subroutine start
  !=======================================================================
  !
  !=======================================================================
  Subroutine printRHO
    !---------------------------------------------------------------------
    ! prints rho, aka
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr), Save :: ihli,Ifle
    Ifle=76+iLST1
    If(Ifle.Eq.76) Open(Ifle,file='ho_den.dat',status='unknown')
    If(Ifle.Eq.77) Open(Ifle,file='tho_den.dat',status='unknown')
    Write(Ifle,*) 'r  denN  denP  akaN  akaP '
    Do ihli=1,nghl
       Write(Ifle,'(12(1x,e16.8))') Sqrt(fh(ihli)**2+fl(ihli)**2)  &
            ,ro(ihli,1),ro(ihli,2),aka(ihli,1),aka(ihli,2)
    End Do
    Close(Ifle)
  End Subroutine printRHO
  !=======================================================================
  !
  !=======================================================================
  Subroutine gfv
    !---------------------------------------------------------------------
    ! Calculates sign, Sqrt, factorials, etc. of integers and half int.
    ! iv(n)=(-1)**n, sq(n)=Sqrt(n), sqi(n)=1/Sqrt(n)
    ! fak(n)=n!; wf(n)=Sqrt(n!); wfi(n)=1/Sqrt(n!)
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: i,igfv
    Parameter(igfv=170)               !maximal number for GFV
    If(Allocated(iv)) Deallocate(iv,fak,fi,sq,sqi,wf,wfi)
    Allocate(iv(-igfv:igfv),fak(0:igfv),fi(0:igfv),sq(0:igfv),sqi(0:igfv))
    Allocate(wf(0:igfv),wfi(0:igfv))
    iv(0)=1; sq(0)=zero; sqi(0)=1.0d30
    fak(0)=one; fi(0)=one; wf(0)=one; wfi(0)=one
    Do i=1,igfv
       iv(i)=-iv(i-1)
       iv(-i) = iv(i)
       sq(i)=Sqrt(Real(i,Kind=pr)); sqi(i)=one/sq(i)
       fak(i)=Real(i,Kind=pr)*fak(i-1); fi(i)=one/fak(i)
       wf(i)=sq(i)*wf(i-1); wfi(i)=one/wf(i)
    End Do
  End Subroutine gfv
  !=======================================================================
  !
  !=======================================================================
  Subroutine sdiag(nmax,n,a,d,x,e,is)
    !---------------------------------------------------------------------
    ! A   matrix to be diagonalized
    ! D   eigenvalues,  X   eigenvectors, E   auxiliary field
    ! IS=1  eigenvalues are ordered (major component of X is positive)
    ! 0  eigenvalues are not ordered
    !---------------------------------------------------------------------
    Use HFBTHO_utilities, Only: pr,ipr
    Implicit None
    Integer(ipr), Save :: i,j,j1,k,l,im
    Integer(ipr)       :: n,nmax,is
    Real(pr), Save :: f,g,h,hi,s,p,b,r,pra,c
    Real(pr) :: a(nmax,nmax),x(nmax,nmax),e(n),d(n)
    Real(pr), Save :: tol=1.0D-32,eps=9.0D-12,one=1.0_pr,zero=0.0_pr
    !
    If (n.Le.1) Then
       d(1)=a(1,1); x(1,1)=one
       Return
    End If
    Do i=1,n
       Do j=1,i
          x(i,j)=a(i,j)
       End Do
    End Do
    ! householder-reduktion
    i=n
15 Continue
    If (i.Ge.2) Then
       l=i-2
       f=x(i,i-1); g=f; h=zero
       If (l.Gt.0) Then
          Do k=1,l
             h=h+x(i,k)*x(i,k)
          End Do
       End If
       s=h+f*f
       If (s.Lt.tol) Then
          h=zero
          Go To 100
       End If
       If (h.Gt.zero) Then
          l=l+1; g=Sqrt(s)
          If (f.Ge.zero) g=-g
          h=s-f*g; hi=one/h; x(i,i-1)=f-g; f=zero
          If (l.Gt.0) Then
             Do j=1,l
                x(j,i)=x(i,j)*hi
                s=zero
                Do k=1,j
                   s=s+x(j,k)*x(i,k)
                End Do
                j1=j+1
                If (l.Ge.j1) Then
                   Do k=j1,l
                      s=s+x(k,j)*x(i,k)
                   End Do
                End If
                e(j)=s*hi; f=f+s*x(j,i)
             End Do
          End If
          f=f*hi*0.50_pr
          If (l.Gt.0) Then
             Do j=1,l
                s=x(i,j); e(j)=e(j)-f*s; p=e(j)
                Do  k=1,j
                   x(j,k)=x(j,k)-s*e(k)-x(i,k)*p
                End Do
             End Do
          End If
       End If
100  Continue
       d(i)=h; e(i-1)=g; i=i-1
       Go To 15
       ! Bereitstellen der Transformationmatrix
    End If
    d(1)=zero; e(n)=zero; b=zero; f=zero
    Do i=1,n
       l=i-1
       If (d(i).Eq.0.) Go To 221
       If (l.Gt.0) Then
          Do J=1,L
             s=zero
             Do k=1,l
                s=s+x(i,k)*x(k,j)
             End Do
             Do k=1,l
                x(k,j)=x(k,j)-s*x(k,i)
             End Do
          End Do
       End If
221  Continue
       d(i)=x(i,i)
       x(i,i)=one
       If (l.Gt.0) Then
          Do j=1,l
             x(i,j)=zero; x(j,i)=zero
          End Do
       End If
    End Do
    ! Diagonalisieren der Tri-Diagonal-Matrix
    Do l=1,n
       h=eps*(Abs(d(l))+ Abs(e(l)))
       If (h.Gt.b) b=h
       ! Test fuer Splitting
       Do  j=l,n
          If (Abs(e(j)).Le.b) Exit
       End Do
       ! test fuer konvergenz
       If (j.Eq.l) Go To 300
340  p=(d(l+1)-d(l))/(2.0_pr*e(l))
       r=Sqrt(p*p+one); pra=p+r
       If (p.Lt.zero) pra=p-r
       h=d(l)-e(l)/pra
       Do i=l,n
          d(i)=d(i)-h
       End Do
       f=f+h
       ! QR-transformation
       p=d(j); c=one; s=zero; i=j
360  i=i-1
       If (i.Lt.l) Go To 362
       g=c*e(i); h=c*p
       If ( Abs(p).Ge.Abs(e(i))) Then
          c=e(i)/p
          r=Sqrt(c*c+one); e(i+1)=s*p*r; s=c/r; c=one/r
          Go To 365
       End If
       c=p/e(i)
       r=Sqrt(c*c+one); e(i+1)=s*e(i)*r; s=one/r; c=c/r
365  p=c*d(i)-s*g
       d(i+1)=h+s*(c*g+s*d(i))
       Do k=1,n
          h=x(k,i+1); x(k,i+1)=x(k,i)*s+h*c
          x(k,i)=x(k,i)*c-h*s
       End Do
       Go To 360
362  e(l)=s*p
       d(l)=c*p
       If ( Abs(e(l)).Gt.b) Go To 340
       ! konvergenz
300  d(l)=d(l)+f
    End Do
    If (is.Eq.0) Return
    ! ordnen der eigenwerte
    Do i=1,n
       k=i; p=d(i); j1=i+1
       If (j1.Le.n) Then
          Do j=j1,n
             If (d(j).Ge.p) Cycle
             k=j; p=d(j)
          End Do
          If (k.Eq.i) Cycle
          d(k)=d(i); d(i)=p
          Do j=1,n
             p=x(j,i); x(j,i)=x(j,k)
             x(j,k)=p
          End Do
       End If
    End Do
    ! signum
    Do  k=1,n
       s=zero
       Do i=1,n
          h=Abs(x(i,k))
          If (h.Gt.s) Then
             s=h; im=i
          End If
       End Do
       If (x(im,k).Lt.zero) Then
          Do i=1,n
             x(i,k)=-x(i,k)
          End Do
       End If
    End Do
  End Subroutine sdiag
  !=======================================================================
  !
  !=======================================================================
  Subroutine nucleus(is,npr2,te)
    !---------------------------------------------------------------------
    ! is=1 determines the symbol for a given proton number npr2
    ! 2 determines the proton number for a given symbol te
    !---------------------------------------------------------------------
    Use HFBTHO_utilities, Only: pr,ipr
    Use HFBTHO, Only: ierror_flag,ierror_info
    Implicit None
    Integer(ipr) :: is,npr2,np
    Integer(ipr) :: maxz
    Parameter (maxz=133)
    Character(2) te
    Character(2*maxz+2) t
    T(  1: 40)=' n HHELIBE B C N O FNENAMGALSI P SCLAR K'
    T( 41: 80)='CASCTI VCRMNFECONICUZNGAGEASSEBRKRRBSR Y'
    T( 81:120)='ZRNBMOTCRORHPDAGCDINSNSBTE IXECSBALACEPR'
    T(121:160)='NDPMSMEUGDTBDYHOERTMYBLUHFTA WREOSIRPTAU'
    T(161:200)='HGTLPBBIPOATRNFRRAACTHPA UNPPUAMCMBKCFES'
    T(201:220)='FMMDNOLR040506070809'
    T(221:265)='101112131415161718192021222324252627282930313233'
    If (is.Eq.1) Then
       If (npr2.Lt.0.Or.npr2.Gt.maxz) Then
          ierror_flag=ierror_flag+1
          ierror_info(ierror_flag)='STOP: in nucleus npr2 is wrong:'
          Return
       End If
       te=t(2*npr2+1:2*npr2+2)
       Return
    Else
       Do np=0,maxz
          If (te.Eq.t(2*np+1:2*np+2)) Then
             npr2=np
             Return
          End If
       End Do
    End If
    ierror_flag=ierror_flag+1
    ierror_info(ierror_flag)='STOP: in nucleus the nucleus is unknown!'
  End Subroutine nucleus
  !=======================================================================
  !
  !=======================================================================
  Subroutine stab(npr2,npr3)
    !---------------------------------------------------------------------
    ! given 'Z' returns mass number 'A' on the stability line
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: npr2,npr3
    Real(pr), Save :: sn,sz,dsn,c,c5
    c=0.0060_pr; c5=5.0_pr*c/3.0_pr; sz=npr2; sn=npr2
    Do While(Abs(dsn).Lt.1.0d-5)
       dsn=sz-sn+c*(sn+sz)**(5.0_pr/3.0_pr)
       dsn=dsn/(-1.0_pr+c5*(sn+sz)**(2.0_pr/3.0_pr))
       sn=sn-dsn
    End Do
    npr3=Int(sn)
  End Subroutine stab
  !=======================================================================
  !
  !=======================================================================
  Subroutine ord(n,e)
    !---------------------------------------------------------------------
    ! orders a set of numbers according to their size
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: n,i,k,j
    Real(pr), Save :: p
    Real(pr) :: e(n)
    Do i=1,n
       k=i; p=e(i)
       If (i.Lt.n) Then
          Do j=i+1,n
             If (e(j).Lt.p) Then
                k=j; p=e(j)
             End If
          End Do
          If (k.Ne.i) Then
             e(k)=e(i); e(i)=p
          End If
       End If
    End Do
  End Subroutine ord
  !=======================================================================
  !
  !=======================================================================
  Subroutine blosort(it,n)
    !---------------------------------------------------------------------
    ! sorting blocking candidates
    ! Integer(ipr) :: iblocking,bloall; Parameter(bloall=200)
    ! Integer(ipr), Dimension(0:bloall,2) :: bloblo,blo123=0,blok1k2=0
    ! Real(pr),    Dimension(0:bloall,2) :: bloqpdif
    ! Integer(ipr), Dimension(3) :: keyblo
    ! Integer(ipr), Dimension(2) :: blocross,blomax,blo123d,blok1k2d,blocanon
    ! Write(lwou) bloblo,blo123,blok1k2,blomax,bloqpdif
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: it,ip,n,i,k,j
    Real(pr) :: p
    Do i=1,n
       k=i; p=bloqpdif(i,it)
       If (i.Lt.n) Then
          Do j=i+1,n
             If (bloqpdif(j,it).Lt.p) Then
                k=j; p=bloqpdif(j,it)
             End If
          End Do
          If (k.Ne.i) Then
             bloqpdif(k,it)=bloqpdif(i,it); bloqpdif(i,it)=p
             ip=bloblo(k,it);  bloblo(k,it)=bloblo(i,it);  bloblo(i,it)=ip
             ip=blo123(k,it);  blo123(k,it)=blo123(i,it);  blo123(i,it)=ip
             ip=blok1k2(k,it); blok1k2(k,it)=blok1k2(i,it); blok1k2(i,it)=ip
          End If
       End If
    End Do
  End Subroutine blosort
  !=======================================================================
  !
  !=======================================================================
  Subroutine tracesln
    !---------------------------------------------------------------------
    !        CALCULATING THE LIPKIN-NOGAMI SUMS IN CANONICAL BASIS
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: iw,it,ib,k1,k2,kkk,k
    Real(pr)    :: AAV,SNtor,SDtor
    Real(pr)    :: S_U1V1,S_U1V3,S_U2V2,S_U3V1,S_U4V4
    Real(pr)    :: U_ACTU,U_ACTU2,U_ACTU3,U_ACTU4
    Real(pr)    :: V_ACTU,V_ACTU2,V_ACTU3,V_ACTU4
    !
    etr=zero
    Do it=itmin,itmax
       S_U1V1=ZERO; S_U1V3=ZERO; S_U2V2=ZERO; S_U3V1=ZERO; S_U4V4=ZERO
       Do ib=1,nb
          k1=ka(ib,it)+1; k2=ka(ib,it)+kd(ib,it)
          If(k1.Le.k2) Then
             kkk=lcanon(ib-1,it)
             Do k=1,id(ib)
                kkk=kkk+1; aav=vk(kkk,it)               ! v^2
                U_ACTU=Sqrt(AAV);       U_ACTU2=U_ACTU*U_ACTU
                U_ACTU3=U_ACTU2*U_ACTU; U_ACTU4=U_ACTU2*U_ACTU2
                V_ACTU=Sqrt(ONE-AAV);   V_ACTU2=V_ACTU*V_ACTU
                V_ACTU3=V_ACTU2*V_ACTU; V_ACTU4=V_ACTU2*V_ACTU2
                S_U1V1=S_U1V1+U_ACTU  * V_ACTU
                S_U1V3=S_U1V3+U_ACTU  * V_ACTU3
                S_U2V2=S_U2V2+U_ACTU2 * V_ACTU2     !Tr r (1-r)
                S_U3V1=S_U3V1+U_ACTU3 * V_ACTU
                S_U4V4=S_U4V4+U_ACTU4 * V_ACTU4     !Tr (1-r)^2 r^2
             End Do
          End If
       End Do !ib
       SNtor=8.0_pr*(S_U3V1*S_U1V3-S_U4V4)
       SDtor=32.0_pr*(S_U2V2*S_U2V2-S_U4V4)
       Geff(it)=del(it)**2/ept(it)
       ala2(it)=-Geff(it)*(SNtor/SDtor)
       If(ala2(it).Ge.10.0_pr) ala2(it)=4.0_pr    ! ala2 goes to hell
       etr(it)=-four*ala2(it)*S_U2V2              ! to total energy
    End Do !it
    etr(3)=etr(1)+etr(2)                          !to total energy
    Do iw=lout,lfile
       Write(iw,'(26x,a,2(1x,f7.3),a,3(1x,f9.3),a,2(1x,f7.3))')  &
            '  f: ala2(n,p)=',ala2,' #eln(n,p,t)=',etr,' #del+ala2=',del+ala2
    End Do
  End Subroutine tracesln
  !=======================================================================
  !
  !=======================================================================
  Subroutine tracesln_qp
    !---------------------------------------------------------------------
    !        CALCULATING THE LIPKIN-NOGAMI TRACES IN QP SPACE
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: iw,nd,ib,i1,i2,n2,ibitnb,it,i1n2nd,i2n2nd,i1i2nd
    Real(pr) :: frit,frit2,ftit
    Real(pr) :: etr2(2),trk(2),trk1(2),SNtor(2),SDtor(2),Sum(2)
    !
    ! initialization
    etr=zero; etr2=zero; trk=zero; trk1=zero
    ! loop over the blocks
    Do ib=1,nb
       nd=id(ib)
       ! Traces for neutrons and protons
       Do i2=1,nd           ! index alpha
          Do i1=i2,nd         ! index beta.ge.alpha
             sum=zero
             Do n2=1,nd
                i1n2nd=Max(i1,n2)+(Min(i1,n2)-1)*nd
                i2n2nd=Max(i2,n2)+(Min(i2,n2)-1)*nd
                Do it=itmin,itmax
                   ibitnb=ib+(it-1)*nbx
                   Sum(it)=Sum(it)+rk(i1n2nd,ibitnb)*rk(i2n2nd,ibitnb)*p14
                End Do !it
             End Do !n2
             i1i2nd=i1+(i2-1)*nd
             Do it=itmin,itmax
                ibitnb=ib+(it-1)*nbx
                frit=rk(i1i2nd,ibitnb)*half
                ftit=ak(i1i2nd,ibitnb)
                frit2=Sum(it)
                If(i1.Eq.i2) Then
                   etr(it)=etr(it)+frit-frit**2                  ! Tr r (1-r)
                   etr2(it)=etr2(it)+(one-two*frit+frit2)*frit2      ! Tr (1-r)^2 r^2
                   trk(it)=trk(it)+frit*ftit                       ! Tr r k
                   trk1(it)=trk1(it)+ftit -frit*ftit                 ! Tr k (1-r)
                Else
                   etr(it)=etr(it) -two*frit**2                     ! Tr r (1-r)
                   etr2(it)=etr2(it)+two*(-two*frit+frit2)*frit2     ! Tr (1-r)^2 r^2
                   trk(it)=trk(it)+two*frit*ftit                   ! Tr r k
                   trk1(it)=trk1(it)-two*ftit*frit                   ! Tr k (1-r)
                End If
             End Do !it
          End Do !i1
       End Do !i2
    End Do !ib
    ! total traces
    Do it=itmin,itmax
       SNtor(it)=8.0_pr*(trk1(it)*trk(it)-etr2(it))
       SDtor(it)=32.0_pr*(etr(it)**2     -etr2(it))
       Geff(it)=del(it)**2/ept(it)
       ala2(it)=-( SNtor(it)/SDtor(it) )*Geff(it)
       If(ala2(it).Ge.10.0_pr) ala2(it)=4.0_pr  ! in case ala2 goes to hell
       etr(it)=-four*ala2(it)*etr(it)           ! to total energy
    End Do
    etr(3)=etr(1)+etr(2)         !to total energy
    Do iw=lout,lfile
       Write(iw,'(26x,a,2(1x,f7.3),a,3(1x,f9.3),a,2(1x,f7.3))')  &
            '  #LN: ala2(n,p)=',ala2,' #eln(n,p,t)=',etr,' #del+ala2=',del+ala2
    End Do
  End Subroutine tracesln_qp
  !=======================================================================
  !
  !=======================================================================
  Subroutine densitln
    !---------------------------------------------------------------------
    ! calculates the densities in r-space at gauss-meshpoints
    ! corrected due to Lipkin-Nogami
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: ih,il,ib,nd,i0,i01,i02,n1,n2,nza,nzb,nra,nrb,nla,  &
         nlb,nsa,nsb,it,ml,ihli,k,k0(2),k00(2),k1,k2
    Real(pr) :: fr(2),vvs,vvc,ssln1,ssln2,ssln3,vks
    Real(pr) :: qla,qlb,qlab,qha,qhb,qhlab,qhab,sro
    !
    k0=0; k00=0
    ro=zero
    ! loop over the blocks
    Do ib=1,nb
       k00=k0
       nd=id(ib); i0=ia(ib)
       Do n2=1,nd
          i02=i0+n2; nzb=nz(i02); nrb=nr(i02);
          nlb=nl(i02); nsb=ns(i02)
          Do n1=1,n2
             i01=i0+n1; nza=nz(i01); nra=nr(i01)
             nla=nl(i01); nsa=ns(i01)
             k0=k00
             Do it=itmin,itmax
                k1=ka(ib,it)+1
                k2=ka(ib,it)+kd(ib,it)
                fr(it)=zero
                If(k1.Le.k2) Then
                   Do k=1,nd
                      k0(it)=k0(it)+1
                      ssln1=ssln(1,it)
                      ssln2=ssln(2,it)
                      ssln3=ssln(3,it)
                      vks=vk(k0(it),it)
                      vvc=vks
                      vvs=Abs(one-vks)
                      If(vvs.Ge.1.0d-40) Then
                         vvs=two*Sqrt(vks*vvs)   !2vu
                         vvc=vks+vvs**2*p14*ssln1*((two*vks-one)*ssln1-ssln2)/ssln3
                      End If
                      fr(it)=fr(it)+two*ddc(n2,k0(it),it)*ddc(n1,k0(it),it)*vvc
                   End Do
                   If (n1.Ne.n2) Then
                      fr(it)=two*fr(it)
                   End If
                End If
             End Do
             !---diagonal in spin
             If (nsa.Eq.nsb) Then
                ml=nla
                Do il=1,ngl
                   qla=ql (nra,ml,il);    qlb=ql (nrb,ml,il)
                   qlab=qla*qlb
                   Do ih=1,ngh
                      ihli=ih+(il-1)*ngh
                      qha=qh (nza,ih); qhb=qh (nzb,ih)
                      qhab=qha*qhb
                      qhlab=qhab*qlab; sro=qhlab
                      ro(ihli,:)=ro(ihli,:)+fr(:)*sro
                   End Do   !ih
                End Do  !il
             End If
          End Do !n2
       End Do !n1
    End Do !ib
    ! set the THO weights
    Do ihli=1,nghl
       ro(ihli,:)=ro(ihli,:)*wdcori(ihli)
    End Do
  End Subroutine densitln
  !=======================================================================
  !
  !=======================================================================
  Subroutine coulom1
    !---------------------------------------------------------------------
    ! Coulomb field (direct part) Vautherin prescription
    ! Ref.: Phys. Rev. C 7, 296 (1973)
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Use EllipticIntegral
    Implicit None
    Integer(ipr), Save :: i,k
    Real(pr) :: zd2,rhl,y1,y2,xx1,xx2,s1,s2,e1,e2,vik,f,r,r1,r4,  &
                rr2,z,z1,zd1,x1,x2,fac1,fac2
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('coulom1',0)
    !
    If(icacou.Eq.0) Then
       !
       icacou=1
       !
       ! For parity-breaking shapes, the Coulomb potential was incorrectly
       ! calculated by assuming the two intervals [0,+\infty[ and ]-infty,0]
       ! were equivalent (see also routine coulom() below). This bug was
       ! corrected in version 200d
       If(Parity) Then
          fac1 = one;  fac2 = one
       Else
          fac1 = zero; fac2 = two
       End If
       !
       f=half*chargee2/pi
       ! See notes in subroutine coulom for explanations about some numerical
       ! factors apparently missing here.
!$OMP PARALLEL DO        &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(nghl,fl,fh,fac1,fac2,wdcor,vc,f) &
!$OMP& PRIVATE(i,r,z,r4,k,r1,z1,rr2,rhl,zd1,y1,xx1,s1,zd2,y2,xx2,s2,vik)
       Do i=1,nghl
          r=fl(i); z=fh(i)
          r4=four*r
          Do k=1,i
             r1=fl(k); z1=fh(k)
             rhl=r4*r1     ! 4 r r'
             rr2=(r+r1)**2 ! (r+r')^2
             ! z>0 part
             zd1=(z-z1)**2 ! (z-z')^2
             y1=zd1+rr2    ! d(r,z) = (r+r')^2 + (z-z')^2
             xx1=rhl/y1    ! 4 r r' / d(r,z)
             s1=Sqrt(y1)   ! sqrt(d(r,z))
             ! z<0 part
             zd2=(z+z1)**2
             y2=zd2+rr2
             xx2=rhl/y2
             s2=Sqrt(y2)
             !
             vik = f*fac2*(s1*CompleteEllipticFunction_2nd(xx1) &
                          +s2*CompleteEllipticFunction_2nd(xx2)*fac1)
             !
             vc(i,k)=vik*wdcor(k)  !wdcor=pi*wh*wl*bz*bp*bp
             vc(k,i)=vik*wdcor(i)  !wdcor=pi*wh*wl*bz*bp*bp
             !
          End Do  !k
       End Do  !i
!$OMP End Parallel Do
    End If
    ! Calculation of the coulomb field (each iteration)
    cou=zero
    Call dgemm('n','n',nghl,1,nghl,1.0_pr,vc,nghl,dro(:,2),nghl,0.0_pr,cou,nghl)
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('coulom1',1)
    !
  End Subroutine coulom1
  !=======================================================================
  !
  !=======================================================================
  Subroutine coulom
    !---------------------------------------------------------------------
    ! Coulomb field (direct part), Gogny prescription
    ! Ref.: Phys. Rev. C 27, 2317 (1983)
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Use bessik
    Implicit None
    Integer(ipr), Save :: i,j,k
    Real(pr), Save :: zd2,y1,y2,xx1,s1,vik,f,r,r1,fac1,fac2,rr2,z,z1,zd1,t,  &
                      bb,r2,r12,rrr,rz1,rz2,rrz1,rrz2,xx,rk1,rip1,rkp1,alpha,&
                      beta,xxx
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('coulom',0)
    !
    If(icacou.Eq.0) Then
       !
       icacou=1
       !
       ! For parity-breaking shapes, the Coulomb potential was incorrectly
       ! calculated by assuming the two intervals [0,+\infty[ and ]-infty,0]
       ! were equivalent (see also below). This bug was corrected in version
       ! 139a
       If(Parity) Then
          fac1 = one;  fac2 = one
       Else
          fac1 = zero; fac2 = two
       End If
       ! Notes:
       !   - Missing factor 2 compared to Eq. (58) CPC paper because the density
       !     ro(:,it) already contains it (see routine DENSIT) due to T-invariance
       !   - Missing factor 1/2 when applying Gauss-Legendre quadrature (from [0,1]
       !     to the proper [-1,1] interval because it will be put back in subroutine
       !     expect() and is cancelled by a factor 2 in the HF field
       !   - For conserved parity, Gauss-Hermite points are all positive, the full
       !     integral over z' is split in z'<0 and z'>0, values of z and z1 below
       !     refer to the absolute values of z' (=-z' if z'<0)
       !
       bb=50.0_pr          ! Length scale L
       beta=2.00_pr
       alpha=one/beta
       f=chargee2/Sqrt(pi) ! e^2/Sqrt(pi)
       !
!$OMP PARALLEL DO        &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(nghl,fl,fh,nleg,xleg,bb,fac1,fac2,wleg,wdcor,vc,f,alpha,beta) &
!$OMP& PRIVATE(i,r,z,k,r1,z1,rrr,rr2,zd1,zd2,rz1,rz2,rrz1,rrz2, &
!$OMP&         xx1,j,xx,y1,s1,t,y2,vik,xxx)
       Do i=1,nghl
          r = fl(i); z = fh(i)
          Do k=1,i
             !
             r1 = fl(k); z1 = fh(k)
             rrr = two*r*r1; rr2 = (r - r1)**2
             ! z>0 part
             zd1 = (z - z1)**2
             rz1 = rr2 + zd1
             ! z<0 part
             zd2 = (z + z1)**2
             rz2 = rr2 + zd2
             ! Gauss-Legendre integration over u from 0 to D
             xx1=zero
             Do j=1,nleg
                xx=(one-xleg(j)**beta)**alpha ! change of variable to 0 <= u <= 1
                xxx=(one-xleg(j)**beta)**(alpha+one)
                y1=(xleg(j)/(bb*xx))**2 ! u^2
                s1=y1*rrr               ! 2 u^2 r r'
                y2=besei0(s1)           ! I0( 2 u^2 r r' ) * exp(-2 u^2 r r')
                xx1=xx1+fac2*wleg(j)*y2*(Exp(-rz1*y1) + fac1*Exp(-rz2*y1)) / xxx
             End Do
             vik=f*xx1/bb
             !
             vc(i,k)=vik*wdcor(k)  !wdcor=pi*wh*wl*bz*bp*bp
             vc(k,i)=vik*wdcor(i)  !wdcor=pi*wh*wl*bz*bp*bp
             !
          End Do  !k
       End Do  !i
!$OMP End Parallel Do
       !
    End If

    ! Calculation of the Coulomb field
    cou=zero
    Call dgemm('n','n',nghl,1,nghl,1.0_pr,vc,nghl,ro(:,2),nghl,0.0_pr,cou,nghl)
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('coulom',1)
    !
  End Subroutine coulom
  !=======================================================================
  !
  !=======================================================================
  Subroutine coulom_test
    !---------------------------------------------------------------------
    ! Coulomb field (direct part), Gogny prescription
    ! Ref.: Phys. Rev. C 27, 2317 (1983)
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Use bessik
    Implicit None
    Integer(ipr), Save :: i,j,k
    Real(pr), Save :: zd2,y1,y2,xx1,s1,vik,f,r,r1,fac1,fac2,rr2,z,z1,zd1,t,  &
                      bb,r2,r12,rrr,rz1,rz2,rrz1,rrz2,xx,rk1,rip1,rkp1,alpha,&
                      beta,xxx,func
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('coulom_test',0)
    !
    !
    ! For parity-breaking shapes, the Coulomb potential was incorrectly
    ! calculated by assuming the two intervals [0,+\infty[ and ]-infty,0]
    ! were equivalent (see also below). This bug was corrected in version
    ! 139a
    If(Parity) Then
       fac1 = one;  fac2 = one
    Else
       fac1 = zero; fac2 = two
    End If
    !
    bb=5.0_pr          ! Length scale L
    beta=2.00_pr
    alpha=one/beta
    !f=chargee2/Sqrt(pi) ! e^2/Sqrt(pi)
    f=one/Sqrt(pi)       ! 1/Sqrt(pi)
    !
    Do j=1,nleg
       ! Gauss-Legendre integration over u from 0 to D
       xx=(one-xleg(j)**beta)**alpha ! change of variable to 0 <= u <= 1
       xxx=(one-xleg(j)**beta)**(alpha+one)
       !
       func=zero
       Do i=1,nghl
          r = fl(i); z = fh(i)
          Do k=1,i
             !
             r1 = fl(k); z1 = fh(k)
             rrr = two*r*r1; rr2 = (r - r1)**2
             ! z>0 part
             zd1 = (z - z1)**2
             rz1 = rr2 + zd1
             ! z<0 part
             zd2 = (z + z1)**2
             rz2 = rr2 + zd2
             y1=(xleg(j)/(bb*xx))**2 ! u^2
             s1=y1*rrr               ! 2 u^2 r r'
             y2=besei0(s1)           ! I0( 2 u^2 r r' ) * exp(-2 u^2 r r')
             xx1=fac2*wleg(j)*y2*(Exp(-rz1*y1) + fac1*Exp(-rz2*y1)) / xxx
             vik=f*xx1/bb
             !
             func=func+vik*wdcor(k)*ro(k,2)*wdcor(i)*ro(i,2)  !wdcor=pi*wh*wl*bz*bp*bp
             !
          End Do ! k
       End Do  ! i
       Write(6,'(2f30.14)') xleg(j),func
       !
    End Do  !j
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('coulom_test',1)
    !
  End Subroutine coulom_test
  !=======================================================================
  !
  !=======================================================================
  Subroutine HartreeDir
    !---------------------------------------------------------------------
    ! Hartree-field (direct part)
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: i,j,k
    Real(pr) :: xx1,vik00,vik01,vik11
    Real(pr) :: r,rr,rrr,r1,r2,rr1,rr2
    Real(pr) :: z,z1,zdm,zdp,rzm,rzp
    Real(pr), Allocatable :: u(:)
    If(icahartree.Eq.0) Then
       icahartree=1
       !
       If(Allocated(u)) Deallocate(u); Allocate(u(nleg))
       u=Cos(HALF*Pi*xleg)
       Do i=1,nghl
          r=fl(i); z=fh(i); rr=r*r
          Do k=1,i
             r1=fl(k); z1=fh(k); rr1=r1*r1;     rr2=two*r*r1; rrr=rr+rr1;
             zdm=(z-z1)**2;      zdp=(z+z1)**2; rzm=rrr+zdm;  rzp=rrr+zdp
             vik00=0.250_pr*Sum(wleg*( &
                  + HartreeV00(Sqrt(rzp-rr2*u)) &
                  + HartreeV00(Sqrt(rzm-rr2*u)) &
                  + HartreeV00(Sqrt(rzp+rr2*u)) &
                  + HartreeV00(Sqrt(rzm+rr2*u))))
             vik01=0.250_pr*Sum(wleg*( &
                  + HartreeV01(Sqrt(rzp-rr2*u)) &
                  + HartreeV01(Sqrt(rzm-rr2*u)) &
                  + HartreeV01(Sqrt(rzp+rr2*u)) &
                  + HartreeV01(Sqrt(rzm+rr2*u))))
             vik11=0.250_pr*Sum(wleg*( &
                  + HartreeV11(Sqrt(rzp-rr2*u)) &
                  + HartreeV11(Sqrt(rzm-rr2*u)) &
                  + HartreeV11(Sqrt(rzp+rr2*u)) &
                  + HartreeV11(Sqrt(rzm+rr2*u))))

             vhart00(i,k)=vik00*wdcor(k)         ! wdcor=pi*wh*wl*bz*bp*bp/fd
             vhart00(k,i)=vik00*wdcor(i)
             vhart01(i,k)=vik01*wdcor(k)
             vhart01(k,i)=vik01*wdcor(i)
             vhart11(i,k)=vik11*wdcor(k)
             vhart11(k,i)=vik11*wdcor(i)
          End Do  !k
       End Do  !i
       Deallocate(u)
    End If
    ! calculation of the Hartree field
    vDHartree=0.0_pr
    Do i=1,nghl
       vDHartree(:,1)=vDHartree(:,1)+vhart00(:,i)*(ro(i,1)+ro(i,2))+vhart01(:,i)*(ro(i,1)-ro(i,2))
       vDHartree(:,2)=vDHartree(:,2)+vhart11(:,i)*(ro(i,1)-ro(i,2))+vhart01(:,i)*(ro(i,1)+ro(i,2))
    End Do
  End Subroutine HartreeDir
  !=======================================================================
  !
  !=======================================================================
  Subroutine optHFBTHO
    !---------------------------------------------------------------------
    ! optimization arrays
    ! NB FI2D_opt(JA,ihil) == Laplacian(r,z) HOwf
    !    FID2D-xlamy2*FID  == Laplacian(r,z,phy) FID
    !    FIU2D-xlapy2*FIU  == Laplacian(r,z,phy) FIU
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: i,ih,il,ib,ibx,nd,nza,nra,nla,nsa
    Integer(ipr) :: ihil,laplus,im,JA,N1,N2,ndnd,n12,n21
    Real(pr)    :: qla,v2,v4,yi,y,y2,qha,qhla,xmi,u,u2,un,up,xxx
    Real(pr)    :: sml2,cnzaa,cnraa,a,b
    Real(pr)    :: FITW1,FITW2,FITW3,FITW4
    Real(pr)    :: fi1r,fi1z,fi2d,QHL1A,QH1LA,vh,vdh,vsh,hbh
    Real(pr)    :: SRFIh,SFIRh,SFIZh,SZFIh,SNABLARh,SNABLAZh
    Real(pr)    :: xlam,xlam2,xlamy,xlamy2,xlap,xlap2,xlapy,xlapy2,XLAMPY
    Real(pr)    :: bpi,bpi2,bzi,bzi2,xh2
    !
    bpi=one/bp; bpi2=bpi*bpi; bzi=one/bz; bzi2=bzi*bzi
    !
    !-----------------------------------------
    ! Allocate the optimization arrays
    !-----------------------------------------
    If(Allocated(QHLA_opt)) Deallocate(QHLA_opt,FI1R_opt,FI1Z_opt,FI2D_opt,y_opt)
    Allocate(QHLA_opt(ntx,nghl),FI1R_opt(ntx,nghl),FI1Z_opt(ntx,nghl),FI2D_opt(ntx,nghl),y_opt(nghl))
    !----------------------------------------------
    ! START BLOCKS
    !----------------------------------------------
    Do ib=1,NB
       ND=ID(ib); IM=ia(ib)
       If(Parity) Then
          LAPLUS=(ib+1)/2 !Yesp
       Else
          LAPLUS=ib       !Nop
       End If
       XLAP=LAPLUS; XLAM=XLAP-ONE; xlap2=xlap*xlap; xlam2=xlam*xlam
       !----------------------------------------------
       ! SUM OVER GAUSS INTEGRATION POINTS
       !----------------------------------------------
       Do IL=1,ngl
          v2=half/xl(il); v4=v2*v2
          Do IH=1,ngh
             ihil=ih+(il-1)*ngh; xh2=xh(ih)**2
             If(iLST1.Eq.0) Then
                ! HO-basis
                yi=Sqrt(xl(il))*bp; y=one/yi; y2=y*y
                xlamy=xlam*y; xlamy2=xlam2*y2; xlapy=xlap*y; xlapy2=xlap2*y2; XLAMPY=XLAMY+XLAPY
             Else
                ! THO-basis
                y=fli(ihil); y2=y*y; xlamy=xlam*y; u=xh(ih); u2=u*u;
                xlamy2=xlam2*y2; xlapy=xlap*y; xlapy2=xlap2*y2; XLAMPY=XLAMY+XLAPY
             End If
             y_opt(ihil)=y
             !----------------------------------------------
             ! SCAN OVER BASIS STATES
             !----------------------------------------------
             Do N1=1,ND
                JA=N1+IM; NLA=NL(JA); NRA=NR(JA); NZA=NZ(JA); NSA=NS(JA)
                SML2=NLA*NLA; CNZAA=NZA+NZA+1; CNRAA=NRA+NRA+NLA+1
                QHA=QH(NZA,IH); QLA=QL(NRA,NLA,IL); QHLA=QHA*QLA
                QHL1A=QHA*QL1(NRA,NLA,IL)*V2; QH1LA=QH1(NZA,IH)*QLA
                If(iLST1.Eq.0) Then
                   ! HO-basis
                   FI1R=(two*Sqrt(xl(il))*bpi)*QHL1A
                   FI1Z=bzi*QH1LA
                   FI2D=((xh2-CNZAA)*bzi2+four*(p14-CNRAA*V2+SML2*V4)*xl(il)*bpi2 )*QHLA
                Else
                   ! THO-basis
                   u=xh(ih); u2=u*u;
                   FI1R=FP4(IHIL)*QHLA+FP5(IHIL)*QH1LA+FP6(IHIL)*QHL1A
                   FI1Z=FP1(IHIL)*QHLA+FP2(IHIL)*QH1LA+FP3(IHIL)*QHL1A
                   FI2D=(FS1(IHIL)*QH1LA*QH1LA+FS2(IHIL)*QHL1A*QHL1A  &
                        +FOUR*FS4(IHIL)*QH1LA*QHL1A  &
                        +TWO*(FS5(IHIL)*QH1LA+FS6(IHIL)*QHL1A)*QHLA  &
                        +((U2-CNZAA)*FS1(IHIL)+(p14-CNRAA*V2+SML2*V4)*FS2(IHIL)  &
                        +FS3(IHIL))*QHLA*QHLA-TWO*(FI1R*FI1R+FI1Z*FI1Z))/(TWO*QHLA)
                End If
                QHLA_opt(JA,ihil)=QHLA; FI2D_opt(JA,ihil)=FI2D; FI1R_opt(JA,ihil)=FI1R; FI1Z_opt(JA,ihil)=FI1Z
             End Do !N1
             !
          End Do !IH
       End Do !IL
    End Do !IB
  End Subroutine optHFBTHO
  !=======================================================================
  !
  !=========================================================================
  Subroutine DENSIT
    !---------------------------------------------------------------------
    ! local densities in coordinate space
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: iw,nsa,nza,nra,nla,k,i,nd,il,ih,ihil,laplus,ii
    Integer(ipr) :: imen,ib,im,it,J,JJ,JA,JN,k0,k1,k2,ibiblo
    Integer(ipr) :: bb,size,ndxmax
    Parameter(ndxmax=(n00max+2)*(n00max+2)/4)
    Real(pr)     :: s,ss,sd,yi,y,y2,sml2,cnzaa,cnraa,u,u2,v2,v4
    Real(pr)     :: anik,pnik,qhla,qh1la,qhl1a,qla,qha,fi1r,fi1z,fi2d,fidd
    Real(pr)     :: xlam,xlam2,xlamy,xlamy2,xlap,xlap2,xlapy,xlapy2,XLAMPY
    Real(pr)     :: TFIU,TFID,TFIUR,TFIDR,TFIUZ,TFIDZ,TFIUD2,TFIDD2
    Real(pr)     :: TPFIU,TPFID,TPFIUR,TPFIDR,TPFIUZ,TPFIDZ,TPFIUD2,TPFIDD2
    Real(pr)     :: PIU,PIUZ,PIUR,PIUD2,PID,PIDZ,PIDR,PIDD2
    Real(pr)     :: TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,TEMP6,TEMP7,TEMP8,TEMP9,TEMP10,TEMP11,TW_T,PW_T,WGT(nghl)
    Real(pr)     :: Takaihil,Troihil,Tdjihil,Ttauihil,Tdroihil,TSRFIihil
    Real(pr)     :: TSFIRihil,TSFIZihil,TSZFIihil,TNABLARIHIL,TNABLAZIHIL
    Real(pr), Pointer :: TAKA(:),TRO(:),TDJ(:),TTAU(:),TDRO(:)
    Real(pr), Pointer :: TSRFI(:),TSFIR(:),TSFIZ(:),TSZFI(:),TNABLAR(:),TNABLAZ(:)
    Real(pr)     :: time1,time2,fk,f1k
    Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:)
    Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
    Real(pr) :: OMPTAKA(nghl,2),OMPTRO(nghl,2),OMPTDJ(nghl,2),OMPTTAU(nghl,2)
    Real(pr) :: OMPTDRO(nghl,2),OMPTSRFI(nghl,2),OMPTSFIR(nghl,2),OMPTSFIZ(nghl,2),OMPTSZFI(nghl,2)
    Real(pr) :: OMPTSZIF(nghl,2),OMPTNABLAR(nghl,2),OMPTNABLAZ(nghl,2)
    !
    Real(pr) :: OMPFIU(ndxmax),OMPFID(ndxmax),OMPFIUR(ndxmax),OMPFIDR(ndxmax),OMPFIUZ(ndxmax)
    Real(pr) :: OMPFIDZ(ndxmax),OMPFIUD2N(ndxmax),OMPFIDD2N(ndxmax)
    !
    Real(pr) :: OMPPFIU(ndxmax),OMPPFID(ndxmax),OMPPFIUR(ndxmax),OMPPFIDR(ndxmax),OMPPFIUZ(ndxmax)
    Real(pr) :: OMPPFIDZ(ndxmax),OMPPFIUD2N(ndxmax),OMPPFIDD2N(ndxmax)
    !
    Real(pr) :: OMPAN(ndxmax*ndxmax),OMPANK(ndxmax*ndxmax)
    Real(pr) :: f_T(ndxmax),f1_T(ndxmax)
    Real(pr) :: dnrm2
    external dnrm2
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('densit',0)
    !
    !-----------------------------------------------
    ! ZERO N & P DENSITIES
    !-----------------------------------------------
    RO=ZERO; TAU=ZERO; DJ=ZERO; DRO=ZERO; AKA=ZERO; SZFI=ZERO; SFIZ=ZERO
    SRFI=ZERO; SFIR=ZERO; NABLAR=ZERO; NABLAZ=ZERO; VARMAS=ZERO
    !
    OMPTAKA   = ZERO; OMPTRO     = ZERO; OMPTDJ     = ZERO; OMPTTAU  = ZERO
    OMPTDRO   = ZERO; OMPTSRFI   = ZERO; OMPTSFIR   = ZERO; OMPTSFIZ = ZERO
    OMPTSZFI  = ZERO; OMPTNABLAR = ZERO; OMPTNABLAZ = ZERO;
    !
    Do bb=0,2*NB-1
       it = bb/NB + 1
       ib = Mod(bb,NB)+1
       !
       ! case of zero particle number, only flush densities
       If((npr_INI(1).Eq.0).And.(it.Eq.1)) Cycle
       If((npr_INI(2).Eq.0).And.(it.Eq.2)) Cycle
       !-----------------------------------------------
       ! SCAN OVER BLOCKS
       !-----------------------------------------------
       ND=ID(ib); IM=ia(ib)
       If(Parity) Then
          LAPLUS=(ib+1)/2 !Yesp
       Else
          LAPLUS=ib       !Nop
       End If
       XLAP=LAPLUS; XLAM=XLAP-ONE; xlap2=xlap*xlap; xlam2=xlam*xlam
       !
       ! blocking
       ibiblo=bloblo(keyblo(it),it)
       K0=0; If(ibiblo.Eq.ib) K0=blo123d(it)
       !
       !----------------------------------------------
       ! PAIRING WINDOW QP WAVE FUNCTIONS
       !----------------------------------------------
       k1=ka(ib,it)+1; k2=ka(ib,it)+kd(ib,it); imen=k2-k1+1
       If(IMEN.Gt.0) Then
          ompan=ZERO; ompank=ZERO; f_T=ZERO; f1_T=ZERO
          J=0
          If(it.Eq.1) then
             Do JJ=1,nd ! basis
                Do K=K1,K2 ! qp
                   J=J+1; I=KpwiN(K)+JJ; ompan(J)=RVqpN(I); ompank(J)=RUqpN(I)
                End Do
             End Do
             J=0
             Do K=K1,K2
                J=J+1;JJ=K !KpwiN(K)
                f_T(J)=one-fn_T(JJ);f1_T(J)=fn_T(JJ)
             End Do
          Else
             Do JJ=1,nd ! basis
                Do K=K1,K2 ! qp
                   J=J+1; I=KpwiP(K)+JJ; ompan(J)=RVqpP(I); ompank(J)=RUqpP(I)
                End Do
             End Do
             J=0
             Do K=K1,K2
                J=J+1;JJ=K !KpwiP(K)
                f_T(J)=one-fp_T(JJ);f1_T(J)=fp_T(JJ)
             End Do
          End If
          !-----------------------------------------------
          ! SCAN OVER GAUSS INTEGRATION POINTS
          !-----------------------------------------------
          Do ihil=1,nghl
             y=y_opt(ihil); xlamy =xlam*y;    xlapy =xlap*y;   XLAMPY=XLAMY+XLAPY
             y2=y*y;        xlamy2=xlam2*y2;  xlapy2=xlap2*y2
             Do K=1,IMEN
                ! V_k components
                OMPFIU(K)    = ZERO; OMPFIUZ(K)   = ZERO; OMPFIUR(K) = ZERO
                OMPFID(K)    = ZERO; OMPFIDZ(K)   = ZERO; OMPFIDR(K) = ZERO
                OMPFIUD2N(K) = ZERO; OMPFIDD2N(K) = ZERO;
                ! U_k components
                OMPPFIU(K)    = ZERO; OMPPFIUZ(K)   = ZERO; OMPPFIUR(K) = ZERO
                OMPPFID(K)    = ZERO; OMPPFIDZ(K)   = ZERO; OMPPFIDR(K) = ZERO
                OMPPFIUD2N(K) = ZERO; OMPPFIDD2N(K) = ZERO;
             End Do
             If(K0.Ne.0) Then
                PIU=ZERO;  PIUZ=ZERO; PIUR=ZERO; PIUD2=ZERO
                PID=ZERO;  PIDZ=ZERO; PIDR=ZERO; PIDD2=ZERO
             End If
             !-----------------------------------------------
             ! SUM OVER BASIS STATES
             !-----------------------------------------------
             JN=0
             Do I=1,ND
                JA=IM+I; NSA=NS(JA); JN=(I-1)*imen
                QHLA=QHLA_opt(JA,ihil); FI2D=FI2D_opt(JA,ihil)
                FI1Z=FI1Z_opt(JA,ihil); FI1R=FI1R_opt(JA,ihil)
                !-----------------------------------------------
                ! QUASIPARTICLE WF IN COORDINATE SPACE
                !-----------------------------------------------
                If (NSA.Gt.0) Then
                   ! SPIN Up
                   Call DAXPY(IMEN,-QHLA,OMPANK(JN+1),1,OMPPFIU,1)
                   ! temperature
                   If(switch_on_temperature) Then
                      Call DAXPY(IMEN,-FI2D,OMPANK(JN+1),1,OMPPFIUD2N,1)
                      Call DAXPY(IMEN,-FI1R,OMPANK(JN+1),1,OMPPFIUR,1)
                      Call DAXPY(IMEN,-FI1Z,OMPANK(JN+1),1,OMPPFIUZ,1)
                   End If
                   Call DAXPY(IMEN, QHLA,OMPAN(JN+1) ,1,OMPFIU,1)
                   Call DAXPY(IMEN, FI2D,OMPAN(JN+1) ,1,OMPFIUD2N,1)
                   Call DAXPY(IMEN, FI1R,OMPAN(JN+1) ,1,OMPFIUR,1)
                   Call DAXPY(IMEN, FI1Z,OMPAN(JN+1) ,1,OMPFIUZ,1)
                   ! blocking
                   If(K0.Ne.0) Then
                      PNIK  = OMPANK(JN+K0)
                      PIU   = PIU   + PNIK*QHLA
                      PIUD2 = PIUD2 + PNIK*FI2D
                      PIUR  = PIUR  + PNIK*FI1R
                      PIUZ  = PIUZ  + PNIK*FI1Z
                   End If
                Else
                 ! SPIN Down
                   Call DAXPY(IMEN,-QHLA,OMPANK(JN+1),1,OMPPFID,1)
                   ! temperature
                   If(switch_on_temperature) Then
                      Call DAXPY(IMEN,-FI2D,OMPANK(JN+1),1,OMPPFIDD2N,1)
                      Call DAXPY(IMEN,-FI1R,OMPANK(JN+1),1,OMPPFIDR,1)
                      Call DAXPY(IMEN,-FI1Z,OMPANK(JN+1),1,OMPPFIDZ,1)
                   End If
                   Call DAXPY(IMEN, QHLA,OMPAN(JN+1) ,1,OMPFID,1)
                   Call DAXPY(IMEN, FI2D,OMPAN(JN+1) ,1,OMPFIDD2N,1)
                   Call DAXPY(IMEN, FI1R,OMPAN(JN+1) ,1,OMPFIDR,1)
                   Call DAXPY(IMEN, FI1Z,OMPAN(JN+1) ,1,OMPFIDZ,1)
                   ! blocking
                   If(K0.Ne.0) Then
                      PNIK  = OMPANK(JN+K0)
                      PID   = PID   + PNIK*QHLA
                      PIDD2 = PIDD2 + PNIK*FI2D
                      PIDR  = PIDR  + PNIK*FI1R
                      PIDZ  = PIDZ  + PNIK*FI1Z
                   End If
                End If
             End Do ! I=1,ND
             !-----------------------------------------------
             ! DENSITIES IN COORDINATE SPACE
             !-----------------------------------------------
             Takaihil=zero;    Troihil=zero;    Tdjihil=zero;   Ttauihil=zero;  Tdroihil=zero
             TSRFIihil=zero;   TSFIRihil=zero;  TSFIZihil=zero; TSZFIihil=zero; TNABLARIHIL=zero; TNABLAZIHIL=zero
             !
             Do K=1,IMEN
                TFIU=OMPFIU(K); TFIUZ=OMPFIUZ(K); TFIUR=OMPFIUR(K); TFIUD2=OMPFIUD2N(K); TPFIU=OMPPFIU(K)
                TFID=OMPFID(K); TFIDZ=OMPFIDZ(K); TFIDR=OMPFIDR(K); TFIDD2=OMPFIDD2N(K); TPFID=OMPPFID(K)
                !
                If(switch_on_temperature) Then
                   !
                   fk=f_T(K); f1k=f1_T(K)
                   !
                   TPFIUZ=OMPPFIUZ(K); TPFIUR=OMPPFIUR(K); TPFIUD2=OMPPFIUD2N(K)
                   TPFIDZ=OMPPFIDZ(K); TPFIDR=OMPPFIDR(K); TPFIDD2=OMPPFIDD2N(K)
                   !
                   TEMP1  = (TPFIU*TFIU+TPFID*TFID)*fk-(TFIU*TPFIU+TFID*TPFID)*f1k
                            TAKAIHIL = TAKAIHIL + TEMP1
                   TEMP2  = (TFIU*TFIU+TFID*TFID)*fk+(TPFIU*TPFIU+TPFID*TPFID)*f1k
                            TROIHIL = TROIHIL + TEMP2
                   TEMP3  = (TFIUR *TFIDZ -TFIDR *TFIUZ +XLAMY*TFIU *(TFIUR -TFIDZ) -XLAPY*TFID *(TFIDR +TFIUZ)) *fk &
                          + (TPFIUR*TPFIDZ-TPFIDR*TPFIUZ+XLAMY*TPFIU*(TPFIUR-TPFIDZ)-XLAPY*TPFID*(TPFIDR+TPFIUZ))*f1k
                            TDJIHIL = TDJIHIL + TEMP3
                   !
                   TW_T=(TFIUR *TFIUR +TFIDR *TFIDR +TFIUZ *TFIUZ +TFIDZ *TFIDZ)*fk&
                       +(TPFIUR*TPFIUR+TPFIDR*TPFIDR+TPFIUZ*TPFIUZ+TPFIDZ*TPFIDZ)*f1k
                   !
                   TEMP4  = (XLAMY2*TFIU *TFIU +XLAPY2*TFID *TFID) *fk &
                          + (XLAMY2*TPFIU*TPFIU+XLAPY2*TPFID*TPFID)*f1k + TW_T
                            TTAUIHIL = TTAUIHIL + TEMP4
                   TEMP5  = (TFIU*TFIUD2+TFID*TFIDD2)*fk + (TPFIU*TPFIUD2+TPFID*TPFIDD2)*f1k + TW_T
                            TDROIHIL = TDROIHIL + TEMP5
                   TEMP6  = (TFIUR*TFID-TFIDR*TFIU)*fk + (TPFIUR*TPFID-TPFIDR*TPFIU)*f1k
                            TSRFIIHIL = TSRFIIHIL + TEMP6
                   TEMP7  = (TFIU*TFID*XLAMPY)*fk + (TPFIU*TPFID*XLAMPY)*f1k
                            TSFIRIHIL = TSFIRIHIL + TEMP7
                   TEMP8  = (XLAMY*TFIU*TFIU-XLAPY*TFID*TFID)*fk + (XLAMY*TPFIU*TPFIU-XLAPY*TPFID*TPFID)*f1k
                            TSFIZIHIL = TSFIZIHIL + TEMP8
                   TEMP9  = (TFIUZ*TFID-TFIDZ*TFIU)*fk + (TPFIUZ*TPFID-TPFIDZ*TPFIU)*f1k
                            TSZFIIHIL = TSZFIIHIL + TEMP9
                   TEMP10 = (TFIUR*TFIU+TFIDR*TFID)*fk + (TPFIUR*TPFIU+TPFIDR*TPFID)*f1k
                            TNABLARIHIL = TNABLARIHIL + TEMP10
                   TEMP11 = (TFIUZ*TFIU+TFIDZ*TFID)*fk + (TPFIUZ*TPFIU+TPFIDZ*TPFID)*f1k
                            TNABLAZIHIL = TNABLAZIHIL + TEMP11
                   !
                Else
                   !
                   TEMP1  = TPFIU*TFIU+TPFID*TFID;                  TAKAIHIL    = TAKAIHIL   + TEMP1
                   TEMP2  = TFIU*TFIU+TFID*TFID;                    TROIHIL     = TROIHIL    + TEMP2
                   TEMP3  = TFIUR*TFIDZ-TFIDR*TFIUZ  &
                           +XLAMY*TFIU*(TFIUR-TFIDZ) &
                           -XLAPY*TFID*(TFIDR+TFIUZ) ;              TDJIHIL     = TDJIHIL    + TEMP3
                   !
                   TW_T=TFIUR*TFIUR+TFIDR*TFIDR+TFIUZ*TFIUZ+TFIDZ*TFIDZ
                   !
                   TEMP4  = XLAMY2*TFIU*TFIU+XLAPY2*TFID*TFID+TW_T; TTAUIHIL    = TTAUIHIL   + TEMP4
                   TEMP5  = TFIU*TFIUD2+TFID*TFIDD2          +TW_T; TDROIHIL    = TDROIHIL   + TEMP5
                   TEMP6  = TFIUR*TFID-TFIDR*TFIU;                  TSRFIIHIL   = TSRFIIHIL  + TEMP6
                   TEMP7  = TFIU*TFID*XLAMPY;                       TSFIRIHIL   = TSFIRIHIL  + TEMP7
                   TEMP8  = XLAMY*TFIU*TFIU-XLAPY*TFID*TFID;        TSFIZIHIL   = TSFIZIHIL  + TEMP8
                   TEMP9  = TFIUZ*TFID-TFIDZ*TFIU;                  TSZFIIHIL   = TSZFIIHIL  + TEMP9
                   TEMP10 = TFIUR*TFIU+TFIDR*TFID;                  TNABLARIHIL = TNABLARIHIL+ TEMP10
                   TEMP11 = TFIUZ*TFIU+TFIDZ*TFID;                  TNABLAZIHIL = TNABLAZIHIL+ TEMP11
                   !
                End If
                !
                If(K.Ne.K0) Cycle
                !
                ! blocking
                TAKAIHIL    = TAKAIHIL    - TEMP1;                 TEMP1  = PIU*PIU+PID*PID
                TROIHIL     = TROIHIL     - HALF*(TEMP2 - TEMP1);  TEMP2  = PIUR*PIDZ-PIDR*PIUZ+XLAMY*PIU*(PIUR-PIDZ) &
                                                                                             -XLAPY*PID*(PIDR+PIUZ)
                !
                PW_T=PIUR*PIUR+PIDR*PIDR+PIUZ*PIUZ+PIDZ*PIDZ
                TDJIHIL     = TDJIHIL     - HALF*(TEMP3 - TEMP2);  TEMP3  = PW_T+XLAMY2*PIU*PIU+XLAPY2*PID*PID
                TTAUIHIL    = TTAUIHIL    - HALF*(TEMP4 - TEMP3);  TEMP4  = PW_T+PIU*PIUD2+PID*PIDD2;
                TDROIHIL    = TDROIHIL    - HALF*(TEMP5 - TEMP4);  TEMP5  = PIUR*PID-PIDR*PIU;
                TSRFIIHIL   = TSRFIIHIL   - HALF*(TEMP6 - TEMP5);  TEMP6  = PIU*PID*XLAMPY;
                TSFIRIHIL   = TSFIRIHIL   - HALF*(TEMP7 - TEMP6);  TEMP7  = XLAMY*PIU*PIU-XLAPY*PID*PID;
                TSFIZIHIL   = TSFIZIHIL   - HALF*(TEMP8 - TEMP7);  TEMP8  = PIUZ*PID-PIDZ*PIU;
                TSZFIIHIL   = TSZFIIHIL   - HALF*(TEMP9 - TEMP8);  TEMP9  = PIUR*PIU+PIDR*PID;
                TNABLARIHIL = TNABLARIHIL - HALF*(TEMP10- TEMP9);  TEMP10 = PIUZ*PIU+PIDZ*PID;
                TNABLAZIHIL = TNABLAZIHIL - HALF*(TEMP11- TEMP10)
             End Do !K
             OMPTaka(ihil,it)    = OMPTaka(ihil,it)    + TAKAIHIL
             OMPTro(ihil,it)     = OMPTro(ihil,it)     + TROIHIL
             OMPTdj(ihil,it)     = OMPTdj(ihil,it)     + TDJIHIL
             OMPTtau(ihil,it)    = OMPTtau(ihil,it)    + TTAUIHIL
             OMPTdro(ihil,it)    = OMPTdro(ihil,it)    + TDROIHIL
             OMPTSRFI(ihil,it)   = OMPTSRFI(ihil,it)   + TSRFIIHIL
             OMPTSFIR(ihil,it)   = OMPTSFIR(ihil,it)   + TSFIRIHIL
             OMPTSFIZ(ihil,it)   = OMPTSFIZ(ihil,it)   + TSFIZIHIL
             OMPTSZFI(ihil,it)   = OMPTSZFI(ihil,it)   + TSZFIIHIL
             OMPTNABLAR(IHIL,IT) = OMPTNABLAR(IHIL,IT) + TNABLARIHIL
             OMPTNABLAZ(IHIL,IT) = OMPTNABLAZ(IHIL,IT) + TNABLAZIHIL
          End Do !ihil
       End If
    End Do !bb
    Do it=1,2
       Do ihil = 1,nghl
          AKA(ihil,it)    = OMPTaka(ihil,it)
          RO(ihil,it)     = OMPTro(ihil,it)
          DJ(ihil,it)     = OMPTdj(ihil,it)
          TAU(ihil,it)    = OMPTtau(ihil,it)
          DRO(ihil,it)    = OMPTdro(ihil,it)
          SRFI(ihil,it)   = OMPTSRFI(ihil,it)
          SFIR(ihil,it)   = OMPTSFIR(ihil,it)
          SFIZ(ihil,it)   = OMPTSFIZ(ihil,it)
          SZFI(ihil,it)   = OMPTSZFI(ihil,it)
          NABLAR(ihil,it) = OMPTNABLAR(ihil,it)
          NABLAZ(ihil,it) = OMPTNABLAZ(ihil,it)
      End Do
    End Do
    Do it = 1,2
       TRO=>ro(:,it);         TTAU=>tau(:,it);       TDJ=>dj(:,it);     TDRO=>dro(:,it)
       TSZFI=>SZFI(:,it);     TSFIZ=>SFIZ(:,it);     TSRFI=>SRFI(:,it); TSFIR=>SFIR(:,it)
       TNABLAR=>NABLAR(:,it); TNABLAZ=>NABLAZ(:,it); TAKA=>aka(:,it)
       s=two*Sum(tro); sd=four*Sum(tdro); drhoi(it)=sd
       !Sumnz(it)=Abs(s-Real(npr(it),Kind=pr))
       Sumnz(it)=Abs(s-tz(it))
       varmas=varmas+s
       !DNFactor(it)=Real(npr(it),Kind=pr)/s
       DNFactor(it)=tz(it)/s
       !----------------------------------------------------
       ! REMOVES INT.WEIGHTS AND MULTIPLIES BY THE JACOBIAN
       !----------------------------------------------------
       !piu=two*Real(npr(it),Kind=pr)/s
       piu=two*tz(it)/s
       WGT=wdcori
       Call dscal(NGHL,piu,WGT,1)
       Tro=Tro*WGT; Ttau=Ttau*WGT; Taka=Half*Taka*WGT;
       TSRFI=TSRFI*WGT; TSFIR=TSFIR*WGT;
       TSFIZ=TSFIZ*WGT; TSZFI=TSZFI*WGT;
       Call dscal(NGHL,two,WGT,1)
       Tdro=Tdro*WGT; Tdj=Tdj*WGT
       TNABLAR=TNABLAR*WGT; TNABLAZ=TNABLAZ*WGT
       !
    End Do !it
    DNFactor(3)=DNFactor(1)+DNFactor(2)
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('densit',1)
    !----------------------------------------------------
    ! COULOMB AND HARTREE FIELDS
    !----------------------------------------------------
    If(nleg.Lt.0) Then
       Call coulom1
    Else
       Call coulom
    End If
    !Call HartreeDir
  End Subroutine DENSIT
  !=======================================================================
  !
  !=======================================================================
  Subroutine field
    !---------------------------------------------------------------------
    ! calculates fields in r-space form axially symmetric densities
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: iw,it,ita,ihli,lambda,icons
    Real(pr) :: ra,ra2,rs,rsa,rsa0,z,rrr
    Real(pr) :: rt,rt1,ds,da,dt,dt1,tts,tta,tt,tt1,djs,dja,djt,djt1
    Real(pr) :: rsa0A,rsa0A1,V0V1,v01a,rns,rps,rsa1,rsa12,rsa10
    Real(pr) :: rsa0An,rsa0An1,rsa0As,rsa0As1
    Real(pr) :: RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1,DJ_0,DJ_1
    Real(pr) :: SZFIN,SFIZN,SRFIN,SFIRN,SZFIP,SFIZP,SRFIP,SFIRP
    Real(pr) :: SZFI_0,SFIZ_0,SRFI_0,SFIR_0,SZFI_1,SFIZ_1,SRFI_1,SFIR_1
    Real(pr) :: SNABLARN,SNABLAZN,SNABLARP,SNABLAZP
    Real(pr) :: SNABLAR_0,SNABLAZ_0,SNABLAR_1,SNABLAZ_1
    Real(pr) :: J2_0,J2_1
    Real(pr) :: cx,x
    Real(pr), Dimension(0:8) :: Qval
    Real(pr),Dimension(2) :: pUr,pUt,pUNr,pUNz,pUDr,pUDj,pUFIZ,pUZFI,pUFIR,pURFI
    Real(pr),Dimension(2) :: tUr,tUt,tUNr,tUNz,tUDr,tUDj,tUFIZ,tUZFI,tUFIR,tURFI
    Real(pr), Save :: ALAMBDA=0.0_pr,AEPSI=1.0_pr,CSPR=1.0_pr
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('field',0)
    !
    ! fields
    Do ihli=1,nghl
       !
       RHO_0 =ro(ihli,1)+ro(ihli,2)     ; RHO_1 =ro(ihli,1)-ro(ihli,2)
       TAU_0 =tau(ihli,1)+tau(ihli,2)   ; TAU_1 =tau(ihli,1)-tau(ihli,2)
       DRHO_0=dro(ihli,1)+dro(ihli,2)   ; DRHO_1=dro(ihli,1)-dro(ihli,2)
       DJ_0  =dj(ihli,1)+dj(ihli,2)     ; DJ_1  =dj(ihli,1)-dj(ihli,2)
       SFIZ_0=SFIZ(ihli,1)+SFIZ(ihli,2) ; SFIZ_1=SFIZ(ihli,1)-SFIZ(ihli,2)
       SFIR_0=SFIR(ihli,1)+SFIR(ihli,2) ; SFIR_1=SFIR(ihli,1)-SFIR(ihli,2)
       SZFI_0=SZFI(ihli,1)+SZFI(ihli,2) ; SZFI_1=SZFI(ihli,1)-SZFI(ihli,2)
       SRFI_0=SRFI(ihli,1)+SRFI(ihli,2) ; SRFI_1=SRFI(ihli,1)-SRFI(ihli,2)
       SNABLAR_0=NABLAR(ihli,1)+NABLAR(ihli,2)
       SNABLAR_1=NABLAR(ihli,1)-NABLAR(ihli,2)
       SNABLAZ_0=NABLAZ(ihli,1)+NABLAZ(ihli,2)
       SNABLAZ_1=NABLAZ(ihli,1)-NABLAZ(ihli,2)
       !
       J2_0=SFIZ_0**2+SFIR_0**2+SZFI_0**2+SRFI_0**2
       J2_1=SFIZ_1**2+SFIR_1**2+SZFI_1**2+SRFI_1**2
       !
       tUr=zero ; tUDr=zero ; tUNr=zero ; tUNz=zero
       tUt=zero ; tUDj=zero ; tUFIZ=zero ; tUZFI=zero
       tUFIR=zero ; tURFI=zero ;
       !
       Call calculate_U_parameters(RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1, &
            (SNABLAR_0**2+SNABLAZ_0**2),(SNABLAR_1**2+SNABLAZ_1**2) )
       !
       ! FUNCTIONAL
       ! E=E+(hb0*(TAU_0+TAU_1)*HALF+hb0*(TAU_0-TAU_1)*HALF  &                         ! tau
       !+Urhotau(0,0)*RHO_0*TAU_0+Urhotau(1,0)*RHO_1*TAU_1  &                         ! rho tau
       !+Urhotau(2,0)*RHO_0*TAU_1+Urhotau(3,0)*RHO_1*TAU_0  &
       !+Urhorho(0,0)*RHO_0**2+Urhorho(1,0)*RHO_1**2  &                               ! rho^2
       !+(Urhorho(2,0)+Urhorho(3,0))*RHO_0*RHO_1  &
       !+UrhoDrho(0,0)*RHO_0*DRHO_0+UrhoDrho(1,0)*RHO_1*DRHO_1  &                     ! rho Delta rho
       !+UrhoDrho(2,0)*RHO_0*DRHO_1+UrhoDrho(3,0)*RHO_1*DRHO_0  &
       !+Unablarho(0,0)*(SNABLAR_0*SNABLAR_0+SNABLAZ_0*SNABLAZ_0)  &                  ! (nabla rho)^2
       !+Unablarho(1,0)*(SNABLAR_1*SNABLAR_1+SNABLAZ_1*SNABLAZ_1)  &
       !+(Unablarho(3,0)+Unablarho(2,0))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1)  &
       !+UrhonablaJ(0,0)*RHO_0*DJ_0+UrhonablaJ(1,0)*RHO_1*DJ_1  &                     ! rho nabla J
       !+UrhonablaJ(2,0)*RHO_0*DJ_1+UrhonablaJ(3,0)*RHO_1*DJ_0  &
       !+UJnablarho(0,0)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))  &     ! J nabla rho
       !+UJnablarho(1,0)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))  &
       !+UJnablarho(2,0)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))  &
       !+UJnablarho(3,0)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1))  &
       !+UJJ(0,0)*J2_0+UJJ(1,0)*J2_1  &                                               ! JJ
       !+(UJJ(3,0)+UJJ(2,0))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)
       !
       ! tUr(1)=dE/d RHO_0;       tUr(2)=dE/d RHO_1
       ! tUt(1)=dE/d TAU_0;       tUt(2)=dE/d TAU_1
       ! tUDr(1)=dE/d DeltaRHO_0; tUDr(2)=dE/d DeltaRHO_1
       ! and so on ...
       !
       !TEST
       !Write(*,'(4(2x,g26.10))') UrhoDrho(0,0)-CrDr(0),CrDr(0),UrhoDrho(1,1)-CrDr(1),CrDr(1); pause
       ! Contributions in the case 'u' depends on RHO_0
       tUr(1)=tUr(1)+two*Urhorho(0,0)*RHO_0+Urhorho(0,1)*RHO_0*RHO_0+Urhorho(1,1)*RHO_1*RHO_1  &  !! rho^2
            +(Urhorho(3,0)+Urhorho(2,0))*RHO_1+(Urhorho(3,1)+Urhorho(2,1))*RHO_0*RHO_1
       tUr(2)=tUr(2)+two*Urhorho(1,0)*RHO_1+Urhorho(0,2)*RHO_0*RHO_0+Urhorho(1,2)*RHO_1*RHO_1  &
            +(Urhorho(3,0)+Urhorho(2,0))*RHO_0+(Urhorho(3,2)+Urhorho(2,2))*RHO_0*RHO_1
       tUr(1)=tUr(1)+vDHartree(ihli,1)
       tUr(2)=tUr(2)+vDHartree(ihli,2)
       !
       tUr(1)=tUr(1)+Urhotau(0,0)*TAU_0+Urhotau(0,1)*TAU_0*RHO_0+Urhotau(1,1)*TAU_1*RHO_1  &  !! rho tau
            +Urhotau(2,0)*TAU_1+Urhotau(2,1)*RHO_0*TAU_1+Urhotau(3,1)*RHO_1*TAU_0
       tUt(1)=tUt(1)+Urhotau(0,0)*RHO_0+Urhotau(3,0)*RHO_1
       tUr(2)=tUr(2)+Urhotau(1,0)*TAU_1+Urhotau(1,2)*TAU_1*RHO_1+Urhotau(0,2)*TAU_0*RHO_0  &
            +Urhotau(3,0)*TAU_0+Urhotau(3,2)*RHO_1*TAU_0+Urhotau(2,2)*RHO_0*TAU_1
       tUt(2)=tUt(2)+Urhotau(1,0)*RHO_1+Urhotau(2,0)*RHO_0
       !
       tUr(1)=tUr(1)+UrhoDrho(0,0)*DRHO_0+UrhoDrho(0,1)*RHO_0*DRHO_0+UrhoDrho(1,1)*RHO_1*DRHO_1  &  !! rho Delta rho
            +UrhoDrho(2,0)*DRHO_1+UrhoDrho(2,1)*RHO_0*DRHO_1+UrhoDrho(3,1)*RHO_1*DRHO_0
       tUDr(1)=tUDr(1)+UrhoDrho(0,0)*RHO_0+UrhoDrho(3,0)*RHO_1
       tUr(2)=tUr(2)+UrhoDrho(1,0)*DRHO_1+UrhoDrho(1,2)*RHO_1*DRHO_1+UrhoDrho(0,2)*RHO_0*DRHO_0  &
            +UrhoDrho(3,0)*DRHO_0+UrhoDrho(3,2)*RHO_1*DRHO_0+UrhoDrho(2,2)*RHO_0*DRHO_1
       tUDr(2)=tUDr(2)+UrhoDrho(1,0)*RHO_1+UrhoDrho(2,0)*RHO_0
       !
       tUr(1)=tUr(1)+Unablarho(0,1)*(SNABLAR_0**2+SNABLAZ_0**2)+Unablarho(1,1)*(SNABLAR_1**2+SNABLAZ_1**2)  &  !! (nabla rho)^2
            +(Unablarho(2,1)+Unablarho(3,1))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1)
       tUNr(1)=tUNr(1)+two*Unablarho(0,0)*SNABLAR_0+(Unablarho(2,0)+Unablarho(3,0))*SNABLAR_1
       tUNz(1)=tUNz(1)+two*Unablarho(0,0)*SNABLAZ_0+(Unablarho(2,0)+Unablarho(3,0))*SNABLAZ_1
       tUr(2)=tUr(2)+Unablarho(0,2)*(SNABLAR_0**2+SNABLAZ_0**2)+Unablarho(1,2)*(SNABLAR_1**2  &
            +SNABLAZ_1**2)+(Unablarho(2,2)+Unablarho(3,2))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1)
       tUNr(2)=tUNr(2)+two*Unablarho(1,0)*SNABLAR_1+(Unablarho(2,0)+Unablarho(3,0))*SNABLAR_0
       tUNz(2)=tUNz(2)+two*Unablarho(1,0)*SNABLAZ_1+(Unablarho(2,0)+Unablarho(3,0))*SNABLAZ_0
       !
       tUr(1)=tUr(1)+UrhonablaJ(0,0)*DJ_0+UrhonablaJ(0,1)*DJ_0*RHO_0+UrhonablaJ(1,1)*DJ_1*RHO_1  &  !! rho nabla J
            +UrhonablaJ(2,0)*DJ_1+UrhonablaJ(2,1)*RHO_0*DJ_1+UrhonablaJ(3,1)*RHO_1*DJ_0
       tUDj(1)=tUDj(1)+UrhonablaJ(0,0)*RHO_0+UrhonablaJ(3,0)*RHO_1
       tUr(2)=tUr(2)+UrhonablaJ(1,0)*DJ_1+UrhonablaJ(1,2)*DJ_1*RHO_1+UrhonablaJ(0,2)*DJ_0*RHO_0  &
            +UrhonablaJ(3,0)*DJ_0+UrhonablaJ(3,2)*RHO_1*DJ_0+UrhonablaJ(2,2)*RHO_0*DJ_1
       tUDj(2)=tUDj(2)+UrhonablaJ(1,0)*RHO_1+UrhonablaJ(2,0)*RHO_0
       !
       tUr(1)=tUr(1)+UJnablarho(0,1)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))      !! J nabla rho
       tUr(1)=tUr(1)+UJnablarho(1,1)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))
       tUr(1)=tUr(1)+UJnablarho(2,1)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))
       tUr(1)=tUr(1)+UJnablarho(3,1)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1))

       tUr(2)=tUr(2)+UJnablarho(0,2)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))
       tUr(2)=tUr(2)+UJnablarho(1,2)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))
       tUr(2)=tUr(2)+UJnablarho(2,2)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))
       tUr(2)=tUr(2)+UJnablarho(3,2)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1))

       tUNr(1)=tUNr(1)+UJnablarho(0,0)*(SFIZ_0-SZFI_0)
       tUNr(2)=tUNr(2)+UJnablarho(1,0)*(SFIZ_1-SZFI_1)
       tUNz(1)=tUNz(1) -UJnablarho(0,0)*(SFIR_0-SRFI_0)
       tUNz(2)=tUNz(2) -UJnablarho(1,0)*(SFIR_1-SRFI_1)

       tUFIZ(1)=tUFIZ(1)+UJnablarho(0,0)*SNABLAR_0*half
       tUFIZ(2)=tUFIZ(2)+UJnablarho(1,0)*SNABLAR_1*half
       tUZFI(1)=tUZFI(1)-UJnablarho(0,0)*SNABLAR_0*half
       tUZFI(2)=tUZFI(2)-UJnablarho(1,0)*SNABLAR_1*half
       tURFI(1)=tURFI(1)+UJnablarho(0,0)*SNABLAZ_0*half
       tURFI(2)=tURFI(2)+UJnablarho(1,0)*SNABLAZ_1*half
       tUFIR(1)=tUFIR(1)-UJnablarho(0,0)*SNABLAZ_0*half
       tUFIR(2)=tUFIR(2)-UJnablarho(1,0)*SNABLAZ_1*half
       !
       !! J.J (Mario: not tested for N2LO)
       tUr(1)=tUr(1)+UJJ(0,1)*J2_0+UJJ(1,1)*J2_1  &
            +(UJJ(3,1)+UJJ(2,1))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)
       tUr(2)=tUr(2)+UJJ(0,2)*J2_0+UJJ(1,2)*J2_1  &
            +(UJJ(3,2)+UJJ(2,2))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)
       tUFIZ(1)=tUFIZ(1)+UJJ(0,0)*SFIZ_0+half*(UJJ(3,0)+UJJ(2,0))*SFIZ_1
       tUFIR(1)=tUFIR(1)+UJJ(0,0)*SFIR_0+half*(UJJ(3,0)+UJJ(2,0))*SFIR_1
       tUZFI(1)=tUZFI(1)+UJJ(0,0)*SZFI_0+half*(UJJ(3,0)+UJJ(2,0))*SZFI_1
       tURFI(1)=tURFI(1)+UJJ(0,0)*SRFI_0+half*(UJJ(3,0)+UJJ(2,0))*SRFI_1
       tUFIZ(2)=tUFIZ(2)+UJJ(1,0)*SFIZ_1+half*(UJJ(3,0)+UJJ(2,0))*SFIZ_0
       tUFIR(2)=tUFIR(2)+UJJ(1,0)*SFIR_1+half*(UJJ(3,0)+UJJ(2,0))*SFIR_0
       tUZFI(2)=tUZFI(2)+UJJ(1,0)*SZFI_1+half*(UJJ(3,0)+UJJ(2,0))*SZFI_0
       tURFI(2)=tURFI(2)+UJJ(1,0)*SRFI_1+half*(UJJ(3,0)+UJJ(2,0))*SRFI_0
       !
       tUr(1)=tUr(1)+UFnonstdr(0)                                        !! other amplitudes
       tUr(2)=tUr(2)+UFnonstdr(1)
       !
       !!  External Field
       !!
       !!tUr(1)=tUr(1)+Vexternal(0,zero,fl(ihli),fh(ihli))
       !!tUr(2)=tUr(2)+Vexternal(1,zero,fl(ihli),fh(ihli))
       !
       ! Contributions in the case 'u' depends on TAU_0
       !
       tUt(1)=tUt(1)+Urhotau(0,6)*RHO_0*TAU_0  &
            +Urhotau(1,6)*RHO_1*TAU_1+Urhotau(2,6)*RHO_0*TAU_1  &
            +Urhotau(3,6)*RHO_1*TAU_0+Urhorho(0,6)*RHO_0**2  &
            +Urhorho(1,6)*RHO_1**2+(Urhorho(2,6)+Urhorho(3,6))*RHO_0*RHO_1  &
            +UrhoDrho(0,6)*RHO_0*DRHO_0+UrhoDrho(1,6)*RHO_1*DRHO_1  &
            +UrhoDrho(2,6)*RHO_0*DRHO_1+UrhoDrho(3,6)*RHO_1*DRHO_0  &
            +Unablarho(0,6)*(SNABLAR_0*SNABLAR_0+SNABLAZ_0*SNABLAZ_0)  &
            +Unablarho(1,6)*(SNABLAR_1*SNABLAR_1+SNABLAZ_1*SNABLAZ_1)  &
            +(Unablarho(2,6)+Unablarho(3,6))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1)  &
            +UrhonablaJ(0,6)*RHO_0*DJ_0+UrhonablaJ(1,6)*RHO_1*DJ_1  &
            +UrhonablaJ(2,6)*RHO_0*DJ_1+UrhonablaJ(3,6)*RHO_1*DJ_0  &
            +UJnablarho(0,6)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))  &
            +UJnablarho(1,6)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))  &
            +UJnablarho(2,6)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))  &
            +UJnablarho(3,6)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1))
       tUt(1)=tUt(1)+UJJ(0,6)*J2_0+UJJ(1,6)*J2_1  &
            +(UJJ(2,6)+UJJ(3,6))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)
       !
       ! Contributions in the case 'u' depends on DeltaRHO_0
       !
       tUDr(1)=tUDr(1)+Urhotau(0,7)*RHO_0*TAU_0  &
            +Urhotau(1,7)*RHO_1*TAU_1+Urhotau(2,7)*RHO_0*TAU_1  &
            +Urhotau(3,7)*RHO_1*TAU_0+Urhorho(0,7)*RHO_0**2  &
            +Urhorho(1,7)*RHO_1**2+(Urhorho(2,7)+Urhorho(3,7))*RHO_0*RHO_1  &
            +UrhoDrho(0,7)*RHO_0*DRHO_0+UrhoDrho(1,7)*RHO_1*DRHO_1  &
            +UrhoDrho(2,7)*RHO_0*DRHO_1+UrhoDrho(3,7)*RHO_1*DRHO_0  &
            +Unablarho(0,7)*(SNABLAR_0*SNABLAR_0+SNABLAZ_0*SNABLAZ_0)  &
            +Unablarho(1,7)*(SNABLAR_1*SNABLAR_1+SNABLAZ_1*SNABLAZ_1)  &
            +(Unablarho(2,7)+Unablarho(3,7))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1)  &
            +UrhonablaJ(0,7)*RHO_0*DJ_0+UrhonablaJ(1,7)*RHO_1*DJ_1  &
            +UrhonablaJ(2,7)*RHO_0*DJ_1+UrhonablaJ(3,7)*RHO_1*DJ_0  &
            +UJnablarho(0,7)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))  &
            +UJnablarho(1,7)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))  &
            +UJnablarho(2,7)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))  &
            +UJnablarho(3,7)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1))
       tUDr(1)=tUDr(1)+UJJ(0,7)*J2_0+UJJ(1,7)*J2_1  &
            +(UJJ(2,7)+UJJ(3,7))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)
       !
       ! proton-neutron representation
       pUr(1)  =tUr(1)+tUr(2);            pUr(2)  =tUr(1)  -tUr(2)
       pUt(1)  =tUt(1)+tUt(2)+hb0*facECM; pUt(2)  =tUt(1)  -tUt(2)+hb0*facECM
       pUDr(1) =tUDr(1)+tUDr(2);          pUDr(2) =tUDr(1) -tUDr(2)
       pUNr(1) =tUNr(1)+tUNr(2);          pUNr(2) =tUNr(1) -tUNr(2)
       pUNz(1) =tUNz(1)+tUNz(2);          pUNz(2) =tUNz(1) -tUNz(2)
       pUDj(1) =tUDj(1)+tUDj(2);          pUDj(2) =tUDj(1) -tUDj(2)
       pUFIZ(1)=tUFIZ(1)+tUFIZ(2);        pUFIZ(2)=tUFIZ(1)-tUFIZ(2)
       pUZFI(1)=tUZFI(1)+tUZFI(2);        pUZFI(2)=tUZFI(1)-tUZFI(2)
       pUFIR(1)=tUFIR(1)+tUFIR(2);        pUFIR(2)=tUFIR(1)-tUFIR(2)
       pURFI(1)=tURFI(1)+tURFI(2);        pURFI(2)=tURFI(1)-tURFI(2)
       !
       Do it=itmin,itmax   !! loop over n  & p
          ita=3-it
          ! constraining potential
          If (numberCons.Gt.0) Then
              z=fh(ihli); rrr=fl(ihli)**2
              Call moments_valueMesh(z,rrr,Qval)
              do icons=1,numberCons
                 lambda=multLambda(icons); pUr(it)= pUr(it) - multLag(lambda)*Qval(lambda)
              end do
          End If
          ! coulomb
          If(it.Eq.2) Then
             If(icou.Ge.1) pUr(it)=pUr(it)+cou(ihli)
             If(icou.Eq.2) pUr(it)=pUr(it)+CExPar*coex*ro(ihli,it)**p13
          End If
          ! pairing contribution to rearrangement term
          If(use_TMR_pairing.eq.0) then
             pUr(it)=pUr(it)-CpV0(it-1) *CpV1(it-1) /rho_c*aka(ihli,it)**2 &
                            -CpV0(ita-1)*CpV1(ita-1)/rho_c*aka(ihli,ita)**2
          endif
          ! pairing contribution to delta dv(ihli,it)
          rsa0=(ro(ihli,it)+ro(ihli,ita))/rho_c
          If(it.Eq.1) Then
             dvn(ihli)=(CpV0(it-1)*(ONE-rsa0*CpV1(it-1)))*aka(ihli,it)
          Else
             dvp(ihli)=(CpV0(it-1)*(ONE-rsa0*CpV1(it-1)))*aka(ihli,it)
          End If
       End Do !it
       !
       vn(ihli)=pUr(1)       ; vp(ihli)=pUr(2)        !* RHO_ij
       vhbn(ihli)=pUt(1)     ; vhbp(ihli)=pUt(2)      !* TAU_ij
       vrn(ihli)=pUNr(1)     ; vrp(ihli)=pUNr(2)      !* NABLAr RHO__ij
       vzn(ihli)=pUNz(1)     ; vzp(ihli)=pUNz(2)      !* NABLAz RHO__ij
       vdn(ihli)=pUDr(1)     ; vdp(ihli)=pUDr(2)      !* DELTA RHO_ij
       vsn(ihli)=pUDj(1)     ; vsp(ihli)=pUDj(2)      !* NABLA . J__ij
       vSFIZn(ihli)=pUFIZ(1) ; vSFIZp(ihli)=pUFIZ(2)  !* JFIZ_ij
       vSZFIn(ihli)=pUZFI(1) ; vSZFIp(ihli)=pUZFI(2)  !* JZFI_ij
       vSFIRn(ihli)=pUFIR(1) ; vSFIRp(ihli)=pUFIR(2)  !* JFIR_ij
       vSRFIn(ihli)=pURFI(1) ; vSRFIp(ihli)=pURFI(2)  !* JRFI_ij
       !
    End Do !ihli
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('field',1)
    !
  End Subroutine field
  !===============================================================================================
  !
  !=======================================================================
  Subroutine gamdel
    !---------------------------------------------------------------------
    ! ph- and pp- matrices in configurational space
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: i,ih,il,ib,ibx,nd,nd2,nza,nra,nla,nsa,nsb,nsab,icons,lambda
    Integer(ipr) :: ihil,laplus,im,JA,N1,N2,ndnd,n12,n21
    Integer(ipr) :: i1,i2,i3
    Real(pr)     :: qla,yi,y,y2,qha,qhla,xmi,u2,un,up,xxx
    Real(pr)     :: sml2,cnzaa,cnraa,SSU,SSD
    Real(pr)     :: FITW1,FITW2,FITW3,FITW4
    Real(pr)     :: fi1r,fi1z,fi2d,QHL1A,QH1LA
    Real(pr)     :: vh,vdh,vsh,hbh,vsum
    Real(pr)     :: SRFIh,SFIRh,SFIZh,SZFIh,SNABLARh,SNABLAZh
    Real(pr)     :: xlam,xlam2,xlamy,xlamy2,xlap,xlap2,xlapy,xlapy2,XLAMPY
    Real(pr)     :: FIUN1,FIDN1,FIURN1,FIDRN1,FIUZN1,FIDZN1,FIUD2N1,FIDD2N1
    Real(pr)     :: FIUN2,FIDN2,FIURN2,FIDRN2,FIUZN2,FIDZN2,FIUD2N2,FIDD2N2
    Real(pr)     :: FIUN12,FIDN12,FIURN12,FIDRN12,FIUZN12,FIDZN12
    Real(pr)     :: vnhl,vrnhl,vznhl,vdnhl,vsnhl,vhbnhl,vSRFInhl,vSFIRnhl
    Real(pr)     :: vSFIZnhl,vSZFInhl,vphl,vrphl,vzphl,vdphl,vsphl,vhbphl
    Real(pr)     :: vSRFIphl,vSFIRphl,vSFIZphl,vSZFIphl,dvnhl,dvphl
    Integer(ipr) :: ibro
    Integer(ipr) :: ndxmax
    Parameter(ndxmax=(n00max+2)*(n00max+2)/4)
    Real(pr) :: OMPFIU(ndxmax),OMPFID(ndxmax),OMPFIUR(ndxmax),OMPFIDR(ndxmax),OMPFIUZ(ndxmax), &
                OMPFIDZ(ndxmax),OMPFIUD2N(ndxmax),OMPFIDD2N(ndxmax)
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('gamdel',0)
    !
    !----------------------------------------------
    ! START BLOCKS
    !----------------------------------------------
    brout=zero; ibro=0
    If (.Not. Allocated(allibro)) Then
       Allocate(allibro(1:NB))
       allibro(1)=0
       Do ib=2,NB
          allibro(ib) = allibro(ib-1) + (ID(ib-1)*(ID(ib-1)+1)/2)
       End Do
    End If
  !$OMP PARALLEL DO        &
  !$OMP& DEFAULT(NONE)     &
  !$OMP& SCHEDULE(DYNAMIC) &
  !$OMP& SHARED(NB,ID,IA,NBX,NS,nghl, &
  !$OMP&        NHHDIM2,NHHDIM3,NHHDIM4,allibro, &
  !$OMP&        vSRFIn,vSFIRn,vSFIZn,vSZFIn, &
  !$OMP&        vSRFIp,vSFIRp,vSFIZp,vSZFIp, &
  !$OMP&        vn,vrn,vzn,vdn,vsn,vhbn,dvn, &
  !$OMP&        vp,vrp,vzp,vdp,vsp,vhbp,dvp, &
  !$OMP&        QHLA_opt,FI1R_opt, FI1Z_opt, FI2D_opt, y_opt, &
  !$OMP&        nhhdim,kindhfb,ALA2,RK,brout,Parity) &
  !$OMP& PRIVATE(I,ND,IB,IM,IBX,LAPLUS,XLAM,XLAP,XLAM2,IL,IH,IHIL,Y,Y2, &
  !$OMP&         XLAMY,XLAMY2,XLAP2,XLAPY,XLAPY2,XLAMPY,N1,JA,NSA,SSU,SSD, &
  !$OMP&         vnhl,vrnhl,vznhl,vdnhl,vsnhl,vhbnhl,dvnhl, &
  !$OMP&         vphl,vrphl,vzphl,vdphl,vsphl,vhbphl,dvphl, &
  !$OMP&         vSRFInhl,vSFIRnhl,vSFIZnhl,vSZFInhl,&
  !$OMP&         vSRFIphl,vSFIRphl,vSFIZphl,vSZFIphl,&
  !$OMP&         FI2D,i1,i2,i3,NSB,NSAB,SNABLARh, SNABLAZh,FI1R,FI1Z, &
  !$OMP&         FIUD2N1,FIDD2N1,FIUD2N2,FIDD2N2,FITW3,FITW4,&
  !$OMP&         OMPFIUD2N,OMPFIDD2N,OMPFIU,OMPFIUR,OMPFIUZ,OMPFID,OMPFIDR,OMPFIDZ, &
  !$OMP&         FIUN1,FIDN1,FIURN1,FIDRN1,FIUZN1,FIDZN1,N2,FIUN2,FIDN2,FIURN2,  &
  !$OMP&         FIDRN2,FIUZN2,FIDZN2,FIUN12,FIDN12,FIURN12,FIDRN12,FIUZN12,FIDZN12,VH,&
  !$OMP&         HBH,VDH,VSH,SRFIH,SFIRH,SFIZH,SZFIH,UN,UP,N12,QHLA)
    Do ib=1,NB
       ND=ID(ib); IM=ia(ib); ibx=ib+nbx
       If(Parity) Then
          LAPLUS=(ib+1)/2 !Yesp
       Else
          LAPLUS=ib       !Nop
       End If
       XLAP=LAPLUS; XLAM=XLAP-ONE; xlap2=xlap*xlap; xlam2=xlam*xlam
       !----------------------------------------------
       ! SUM OVER GAUSS INTEGRATION POINTS
       !----------------------------------------------
       Do ihil=1,nghl
          y=y_opt(ihil); xlamy=xlam*y;     xlapy=xlap*y;   XLAMPY=XLAMY+XLAPY
          y2=y*y;        xlamy2=xlam2*y2;  xlapy2=xlap2*y2
          !
          vnhl=vn(ihil);         vrnhl=vrn(ihil);       vznhl=vzn(ihil);       vdnhl=vdn(ihil)
          vsnhl=vsn(ihil);       vhbnhl=vhbn(ihil);     vSRFInhl=vSRFIn(IHIL); vSFIRnhl=vSFIRn(IHIL)
          vSFIZnhl=vSFIZn(IHIL); vSZFInhl=vSZFIn(IHIL); vphl=vp(ihil);         vrphl=vrp(ihil)
          vzphl=vzp(ihil);       vdphl=vdp(ihil);       vsphl=vsp(ihil);       vhbphl=vhbp(ihil)
          vSRFIphl=vSRFIp(IHIL); vSFIRphl=vSFIRp(IHIL); vSFIZphl=vSFIZp(IHIL); vSZFIphl=vSZFIp(IHIL)
          dvnhl=dvn(ihil);       dvphl=dvp(ihil)
          !
          Do N1=1,ND
             JA=IM+N1;               NSA=NS(JA);             SSU=Max(NSA,0);         SSD=Max(-NSA,0)
             QHLA=QHLA_opt(JA,ihil); FI1R=FI1R_opt(JA,ihil); FI1Z=FI1Z_opt(JA,ihil); FI2D=FI2D_opt(JA,ihil)
             OMPFIU(N1)=QHLA*SSU;    OMPFIUR(N1)=fi1r*SSU
             OMPFIUZ(N1)=fi1z*SSU;   OMPFIUD2N(N1)=(FI2D-XLAMY2*QHLA)*SSU
             OMPFID(N1)=QHLA*SSD;    OMPFIDR(N1)=fi1r*SSD
             OMPFIDZ(N1)=fi1z*SSD;   OMPFIDD2N(N1)=(FI2D-XLAPY2*QHLA)*SSD
          End Do
          !
          I=allibro(ib)
          Do N1=1,ND
             JA=IM+N1;               NSA=NS(JA)
             FIUN1=OMPFIU(N1);       FIURN1=OMPFIUR(N1);
             FIUZN1=OMPFIUZ(N1);     FIUD2N1=OMPFIUD2N(N1)
             FIDN1=OMPFID(N1);       FIDRN1=OMPFIDR(N1);
             FIDZN1=OMPFIDZ(N1);     FIDD2N1=OMPFIDD2N(N1)
             Do N2=1,N1
                I=I+1; i1=i+nhhdim; i2=i+nhhdim2; i3=i+nhhdim3; NSB=NS(N2+IM); NSAB=NSA+NSB
                If (NSAB.Ne.0) Then
                    If (NSB.Gt.0) Then                                    !spin:UpUp
                      FIUN2    = OMPFIU(N2);    FIURN2 = OMPFIUR(N2)
                      FIUD2N2  = OMPFIUD2N(N2); FIUZN2 = OMPFIUZ(N2)
                      vh       = FIUN1*FIUN2
                      hbh      = vh*XLAMY2+FIURN1*FIURN2+FIUZN1*FIUZN2
                      vdh      = hbh+hbh+FIUN1*FIUD2N2+FIUN2*FIUD2N1
                      SNABLARh = FIURN1*FIUN2+FIURN2*FIUN1
                      SNABLAZh = FIUZN1*FIUN2+FIUZN2*FIUN1
                      vsh      = SNABLARh*XLAMY
                      SFIZh    = (vh+vh)*XLAMY ! =SFIZh (v103)
                   Else                                                  !spin:DoDo
                      FIDN2    = OMPFID(N2);  FIDRN2  = OMPFIDR(N2);
                      FIDZN2   = OMPFIDZ(N2); FIDD2N2 = OMPFIDD2N(N2)
                      vh       = FIDN1*FIDN2
                      hbh      = vh*XLAPY2+FIDRN1*FIDRN2+FIDZN1*FIDZN2
                      vdh      = hbh+hbh+FIDN1*FIDD2N2+FIDN2*FIDD2N1;
                      SNABLARh = FIDRN1*FIDN2+FIDRN2*FIDN1
                      SNABLAZh = FIDZN1*FIDN2+FIDZN2*FIDN1
                      vsh      =-SNABLARh*XLAPY
                      SFIZh    =-(vh+vh)*XLAPY ! =SFIZh (v103)
                   End If
                   brout(i )=brout(i )+vSFIZnhl*SFIZh+vh*vnhl+SNABLARh*vrnhl+SNABLAZh*vznhl+vdh*vdnhl+vsh*vsnhl+hbh*vhbnhl
                   brout(i1)=brout(i1)+vSFIZphl*SFIZh+vh*vphl+SNABLARh*vrphl+SNABLAZh*vzphl+vdh*vdphl+vsh*vsphl+hbh*vhbphl
                   brout(i2)=brout(i2)+vh*dvnhl
                   brout(i3)=brout(i3)+vh*dvphl
                Else
                   If (NSB.Gt.0) Then                                                                !spin:DoUp
                      !vh=ZERO; hbh=ZERO; vdh=ZERO; SNABLARh=ZERO; SNABLAZh=ZERO; SFIZh=ZERO
                      FIUN2   = OMPFIU(N2);    FIURN2 = OMPFIUR(N2);
                      FIUD2N2 = OMPFIUD2N(N2); FIUZN2 = OMPFIUZ(N2)
                      FITW3   =-FIDZN1*FIUN2; FITW4=FIUZN2*FIDN1
                      vsh     =-FIDRN1*FIUZN2+FIURN2*FIDZN1+FITW3*XLAMY-FITW4*XLAPY
                      SRFIh   =-FIDRN1*FIUN2+FIURN2*FIDN1
                      SFIRh   = FIDN1*FIUN2*XLAMPY
                      SZFIh   = FITW3+FITW4
                   Else                                                                             !spin:UpDo
                      !vh=ZERO; hbh=ZERO; vdh=ZERO; SNABLARh=ZERO; SNABLAZh=ZERO; SFIZh=ZERO
                      FIDN2  = OMPFID(N2);   FIDRN2  = OMPFIDR(N2);
                      FIDZN2 = OMPFIDZ(N2);  FIDD2N2 = OMPFIDD2N(N2)
                      FITW3  =-FIDZN2*FIUN1; FITW4=FIUZN1*FIDN2
                      vsh    = FIURN1*FIDZN2-FIDRN2*FIUZN1-FITW4*XLAPY+FITW3*XLAMY ! -vsh (v103)
                      SRFIh  = FIURN1*FIDN2-FIDRN2*FIUN1 !=SRFIh (v103)
                      SFIRh  = FIUN1*FIDN2*XLAMPY !=SFIRh(v103)
                      SZFIh  = FITW3+FITW4 !=SZFIh(v103)
                   End If
                   brout(i )=brout(i )+vsh*vsnhl+vSRFInhl*SRFIh+vSFIRnhl*SFIRh+vSZFInhl*SZFIh
                   brout(i1)=brout(i1)+vsh*vsphl+vSRFIphl*SRFIh+vSFIRphl*SFIRh+vSZFIphl*SZFIh
                End If
                !----------------------------------------------
                ! LN PH PART
                !----------------------------------------------
                If(kindhfb.Lt.0) Then
                   If(ihil.Eq.1) Then
                      un=zero; up=zero;
                      If(N1.Eq.N2) Then
                         un=-ala2(1); up=-ala2(2)
                      End If
                      n12=N1+(N2-1)*ND
                      brout(i )=brout(i )+two*(ala2(1)*rk(n12,ib )+un)
                      brout(i1)=brout(i1)+two*(ala2(2)*rk(n12,ibx)+up)
                   End If
                End If
             End Do !N2
          End Do !N1
       End Do !ihil
    End Do !IB
  !$OMP End Parallel Do
    If (IDEBUG.Eq.1) Call get_CPU_time('gamdel',1)
    !
  ! Lagrange parameters for the constraints
    Do lambda=1,lambdaMax
       brout(nhhdim4+lambda)=multLag(lambda)
    End Do
    !----------------------------------------------
    ! BROYDEN/LINEAR MIXING
    !----------------------------------------------
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('broyden',0)
    !
    Call broyden_min(nhhdim4+lambdaMax,brout,brin,alphamix,si,iiter,nbroyden,bbroyden)
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('broyden',1)
    !
    Do lambda=1,lambdaMax
       multLag(lambda)=brin(nhhdim4+lambda)
    End Do
    !
  End Subroutine gamdel
  !=======================================================================
  !
  !======================================
  ! lib PnProjected specIfics Start >>>>>
  !======================================
  !=======================================================================
  Subroutine expectpj(lpr)
    !---------------------------------------------------------------------
    ! calculates expectation values (tz is the particle number)
    ! optimized for half gauge-angle points
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Use UNEDF
    Implicit None
    Logical :: lpr,part_is_ready
    Integer(ipr) :: i,j,it,ihli,iw,iw1=901,iw2=902,icons,lambda
    Complex(pr) :: SZFIN,SFIZN,SRFIN,SFIRN,SZFIP,SFIZP,SRFIP,SFIRP
    Complex(pr) :: SFIZ_0,SFIR_0,SZFI_0,SRFI_0,SFIZ_1,SFIR_1,SZFI_1,SRFI_1
    Complex(pr) :: RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1,DJ_0,DJ_1
    Complex(pr) :: SNABLAR_0,SNABLAZ_0,SNABLAR_1,SNABLAZ_1
    Complex(pr) :: cekt(3),cdel(2),cept(3),cetot,etens,cq2pj(ilpj,ilpj)
    Complex(pr) :: cxn(2),crms(3),cq2(3),cq4(3),xnpj(2),rmspj(3),q2pj(3),q4pj(3)
    Complex(pr) :: evolpj,esurpj,ecdipj,ecexpj,ecoupj,ept1pj,ept2pj,epotpj              &
         ,eki1pj,eki2pj,ekinpj,etotpj,espopj,epa1pj,epa2pj,epirpj,ede1pj,ede2pj,etenspj &
         ,eva,ev3,ev5,es5,eso,ecodi,ecoex,rn,rp,rnp1,rnp2,rt,rt2,tnt,tpt,tt                 &
         ,dn,dp,dt,akn,akp,akn2,akp2,adn,adp,evol,esurf,ecoul,pijk,row,cx,dd1n,dd1p     &
         ,rt1,tt1,dt1,djn,djp,djt,djt1
    Complex(pr) :: rsa,rsa0,rsa0A,rps,rns,rsa1,rsa10,rsa12,rsa0An,rsa0As
    Real(pr) :: whl,x,xn(3),q4(3),def(3),bet2(3),het4(3),r212,r222,rc,z,zz,rrr,p2,p3,p4
    Real(pr)    :: rdelta(2),repair(3),rekin(3),revolpj,resurpj,respopj,recdipj,recexpj,retenspj
    !
    Call densitpj ! calculates complex densities and the direct coulomb field
    !
    evolpj = zero; esurpj = zero; ecdipj = zero; ecexpj = zero; espopj  = zero;
    ept1pj = zero; ept2pj = zero; eki1pj = zero; eki2pj = zero; epj     = zero;
    etotpj = zero; epa1pj = zero; epa2pj = zero; ede1pj = zero; ede2pj  = zero;
    xnpj   = zero; rmspj  = zero; q2pj   = zero; q4pj   = zero; etenspj = zero; cq2pj = zero
    !
    Do i=1,ilpj
       Do j=1,ilpj
          pijk = pjk(i,1)*pjk(j,2)
          !
          cekt = zero; cept  = zero; cdel  = zero;
          cxn  = zero; crms  = zero; cq2   = zero; cq4   = zero;
          eva  = zero; ev3   = zero; ev5   = zero; es5   = zero;
          eso  = zero; ecodi = zero; ecoex = zero; etens = zero;
          !
          Do ihli = 1,nghl
             ! real
             whl = wdcor(ihli)
             z   = fh(ihli); zz = z*z; rrr = zz + fl(ihli)**2
             p2  = p32*zz   - half*rrr    !3/2 z*z-1/2 (z*z+r*r)=1/2(2 z*z-r*2)=1/2 Q
             p3  = p53*z*p2 - p23*rrr*z
             p4  = p74*z*p3 - p34*rrr*p2
             ! complex
             rn  = ropj(ihli,i,1); rp  = ropj(ihli,j,2); rnp2 = rn**2 + rp**2; rnp1=rn - rp
             ! ig - particle number, rms and deformations
             row = whl*rn; cxn(1)=cxn(1)+row; crms(1)=crms(1)+row*rrr; cq2(1)=cq2(1)+row*p2; cq4(1)=cq4(1)+row*p4
             row = whl*rp; cxn(2)=cxn(2)+row; crms(2)=crms(2)+row*rrr; cq2(2)=cq2(2)+row*p2; cq4(2)=cq4(2)+row*p4
             ! ig - energy contributions
             rt   = rn + rp;     rt2  = rt*rt
             tnt  = taupj(ihli,i,1); tpt  = taupj(ihli,j,2); tt   = tnt + tpt
             dn   = dropj(ihli,i,1); dp   = dropj(ihli,j,2); dt   = dn + dp
             akn  = akapj(ihli,i,1); akp  = akapj(ihli,j,2)
             akn2 = akn*akn;         akp2 = akp*akp
             adn  = akn*rn;          adp  = akp*rp
             ! ig-Pairing energy and delta
             rsa0=(rt/rho_c)
             dd1n=CpV0(0)*(ONE-rsa0*CpV1(0))*whl
             dd1p=CpV0(1)*(ONE-rsa0*CpV1(1))*whl
             !
             cept(1) = cept(1) + dd1n*akn2; cept(2) = cept(2) + dd1p*akp2
             cdel(1) = cdel(1) - dd1n*adn;  cdel(2) = cdel(2) - dd1p*adp
             !
             x       = hb0*whl
             cekt(1) = cekt(1) + x*tnt; cekt(2) = cekt(2) + x*tpt       !kinetic
             ev3     = ev3 + (tv1*rt2 - tv2*rnp2)*whl                   !volume
             eva     = eva + (tv3*rt2-tv4*rnp2)*rt**sigma*whl
             ev5     = ev5 + (tv5*rt*tt + tv6*(rn*tnt + rp*tpt))*whl
             es5     = es5 + (ts1*rt*dt + ts2*(rn*dn + rp*dp))*whl        !surface
             eso     = eso + (CrdJ(0)*rt+CrdJ(1)*rnp1)*djpj(ihli,i,1)*whl !spin-orbit
             eso     = eso + (CrdJ(0)*rt-CrdJ(1)*rnp1)*djpj(ihli,j,2)*whl !spin-orbit
             If(icou.Ge.1) ecodi = ecodi + half*coupj(ihli,j)*rp*whl      !Coul.dir
             If(icou.Eq.2) ecoex = ecoex - cex*rp**t4o3*whl               !Coul.exc
             If(use_j2terms) Then
                SFIZN=SFIZpj(IHLI,i,1); SFIRN=SFIRpj(IHLI,i,1); SZFIN=SZFIpj(IHLI,i,1); SRFIN=SRFIpj(IHLI,i,1)
                SFIZP=SFIZpj(IHLI,j,2); SFIRP=SFIRpj(IHLI,j,2); SZFIP=SZFIpj(IHLI,j,2); SRFIP=SRFIpj(IHLI,j,2)
                ETENS=ETENS+whl*(TA7*(SZFIN**2+SFIZN**2+SRFIN**2+SFIRN**2+SZFIP**2+SFIZP**2+SRFIP**2+SFIRP**2)&
                     +TA8*(SZFIN*SZFIP+SFIZN*SFIZP+SRFIN*SRFIP+SFIRN*SFIRP))
             End If
          End Do !ihli
          !
          evol     = ev3 + eva + ev5; esurf = es5; ecoul = ecodi  + ecoex*CExPar
          cekt(3)  = cekt(1) + cekt(2); cept(3)  = cept(1) + cept(2)
          cetot    = cekt(3) + evol + esurf + eso + ecoul + cept(3)+ ETENS
          cdel(1)  = cdel(1)/tz(1); cdel(2)  = cdel(2)/tz(2)
          !------------------------------------------------
          ! half-projected energies required for the matrix elements
          !------------------------------------------------
          epj(i,1) = epj(i,1) + cetot*pjk(j,2)
          epj(j,2) = epj(j,2) + cetot*pjk(i,1)
          !------------------------------------------------
          ! for constraint contributions to half-projected energies
          !------------------------------------------------
          If (icstr.Ne.0) cq2pj(i,j)=two*(cq2(1)+cq2(2))
          !------------------------------------------------
          ! projected energies
          !------------------------------------------------
          evolpj = evolpj + pijk*evol;    esurpj = esurpj + pijk*esurf;   espopj = espopj + pijk*eso
          epa1pj = epa1pj + pijk*cept(1); epa2pj = epa2pj + pijk*cept(2); epirpj = epa1pj + epa2pj
          ede1pj = ede1pj + pijk*cdel(1); ede2pj = ede2pj + pijk*cdel(2)
          ecdipj = ecdipj + pijk*ecodi;   ecexpj = ecexpj + pijk*ecoex;   ecoupj = ecdipj + ecexpj
          ept1pj = ept1pj + pijk*cept(1); ept2pj = ept2pj + pijk*cept(2); epotpj = ept1pj + ept2pj
          eki1pj = eki1pj + pijk*cekt(1); eki2pj = eki2pj + pijk*cekt(2); ekinpj = eki1pj + eki2pj
          !
          etotpj = etotpj + pijk*cetot
          etenspj= etenspj+ pijk*etens
          !------------------------------------------------
          ! unprojected hfb total energy and constraint
          !------------------------------------------------
          If(i.Eq.1.And.j.Eq.1) Then
             rehfbcan=Real(cetot,Kind=pr)
          End If
          !
          ! projected particle numbers, rms, deformations
          If(j.Eq.1) Then
             xnpj(1)  = xnpj(1)  + pjk(i,1)*cxn(1)
             rmspj(1) = rmspj(1) + pjk(i,1)*crms(1)
             q2pj(1)  = q2pj(1)  + pjk(i,1)*cq2(1)
             q4pj(1)  = q4pj(1)  + pjk(i,1)*cq4(1)
          End If
          If(i.Eq.1) Then
             xnpj(2)  = xnpj(2)  + pjk(j,2)*cxn(2)
             rmspj(2) = rmspj(2) + pjk(j,2)*crms(2)
             q2pj(2)  = q2pj(2)  + pjk(j,2)*cq2(2)
             q4pj(2)  = q4pj(2)  + pjk(j,2)*cq4(2)
          End If
          !
       End Do !j
    End Do !i
    !
    ! Real quantities to the end
    !
    !------------------------------------------------
    ! Energies
    !------------------------------------------------
    rdelta(1) = Real(ede1pj,Kind=pr); rdelta(2) = Real(ede2pj,Kind=pr); retotpj   = Real(etotpj,Kind=pr);
    repair(1) = Real(epa1pj,Kind=pr); repair(2) = Real(epa2pj,Kind=pr); repair(3) = Real(epirpj,Kind=pr)
    rekin(1)  = Real(eki1pj,Kind=pr); rekin(2)  = Real(eki2pj,Kind=pr); rekin(3)  = Real(ekinpj,Kind=pr)
    revolpj   = Real(evolpj,Kind=pr); resurpj   = Real(esurpj,Kind=pr); respopj   = Real(espopj,Kind=pr);
    recdipj   = Real(ecdipj,Kind=pr); recexpj   = Real(ecexpj,Kind=pr); retenspj  = Real(etenspj,Kind=pr);
    depnp = retotpj - rehfbcan  !correlation energy due to projection
    !------------------------------------------------
    ! expectation values of multipole moments
    !------------------------------------------------
    Call moments_computeValue()
    !------------------------------------------------
    ! rms and deformations
    !------------------------------------------------
    Do it=itmin,itmax
       xn(it) = Real(xnpj(it),Kind=pr)
       rms(it)= Sqrt(Real(rmspj(it),Kind=pr)/xn(it))
       q2(it) = two*Real(q2pj(it),Kind=pr)    !Qnp=<2r^2P_2(teta)>=<2z^2-x^2-y^2>
       q4(it) = ffdef4*Real(q4pj(it),Kind=pr) !Hn=<8r^4P_4(teta)>=<8z^4-24z^2(x^2+y^2)+3(x^2+y^2)^2>
       def(it)= Sqrt(pi/5.0_pr)*q2(it)/(rms(it)**2*xn(it))
    End Do
    r212    = rms(1)**2; r222 = rms(2)**2
    rms(3)  = Sqrt((xn(1)*r212+xn(2)*r222)/amas)
    q2(3)   = q2(1) + q2(2)  ! quadrupole moment
    q4(3)   = q4(1) + q4(2)  ! hexadecapole moment
    def(3)  = Sqrt(pi/5.0_pr)*q2(3)/(rms(3)**2*amas) !deformation
    !------------------------------------------------
    ! other definitions of the same quantities
    !------------------------------------------------
    bet2(1) = ffdef6*q2(1)/(xn(1)*r02) !beta_n=Qn*Sqrt(5Pi)/(3N x^2)
    bet2(2) = ffdef6*q2(2)/(xn(2)*r02) !x=r0=1.2A^(1/3)
    bet2(3) = ffdef6*q2(3)/(amas*r02)
    het4(1) = ffdef7*q4(1)/(xn(1)*r04)
    het4(2) = ffdef7*q4(2)/(xn(2)*r04)
    het4(3) = ffdef7*q4(3)/(amas*r04)
    xn(3)   = xn(1) + xn(2)
    bet = def(3)
    !------------------------------------------------
    !  constraint constants and contributions to half-projected energies
    !------------------------------------------------
    If(icstr.Ne.0) Then
       cx=0.0_pr
       If (numberCons.Gt.0) Then
           Do icons=1,numberCons
              lambda=multLambda(icons)
              cx = cx - multLag(lambda)*(qmoment(lambda,3)-multRequested(lambda))
           End Do
       End If
       !ty20=Sqrt(5.0_pr/pi)*hom/b0**2/two
       !cx=cqad*(cdef-bet)*ty20;
       Do i=1,ilpj
          Do j=1,ilpj
             epj(i,1) = epj(i,1) + cx*cq2pj(i,j)*pjk(j,2)
             epj(j,2) = epj(j,2) + cx*cq2pj(i,j)*pjk(i,1)
          End Do
       End Do
    End If
    !
    If (lpr) Then
       rc=Sqrt(r222+0.640_pr)
       ! transitions to barn,barn^2,barn^4
       Do i=1,3
          q2(i)=q2(i)/100.0_pr; q4(i)=q4(i)/10000.0_pr
       End Do
       !
       ! STORE to projected buffer 'eresj'
       ! ieresj=50 from module definitions
       ! ' si ','JININ'
       eresj(1)=si; eresj(2)=inin;
       ! ' A','   N ','   Z '
       eresj(3)=npr(1)+npr(2); eresj(4)=npr(1); eresj(5)=npr(2);
       ! ' Jln ',' Jlp '
       eresj(6)=alast(1); eresj(7)=alast(2);
       ! ,'JEtot','Jbett','Jbetn','Jbetp',' JQt ',' JQn ',' JQp '  &
       eresj(8)=retotpj; eresj(9)=def(3); eresj(10)=def(1); eresj(11)=def(2);
       eresj(12)=q2(3); eresj(13)=q2(1); eresj(14)=q2(2);
       ! ' JpEn',' JpEp',' JpDn',' JpDp',' JAsn',' JAsp'  &
       eresj(15)=repair(1); eresj(16)=repair(2);
       eresj(17)=rdelta(1); eresj(18)=rdelta(2); eresj(19)=ass(1); eresj(20)=ass(2);
       ! ,' Jrt ',' Jrn ',' Jrp ',' Jrc ',' Jht ',' Jhn ',' Jhp '  &
       eresj(21)=rms(3); eresj(22)=rms(1); eresj(23)=rms(2); eresj(24)=rc;
       eresj(25)=het4(3); eresj(26)=het4(1); eresj(27)=het4(2);
       ! ,' Jqht',' Jqhn',' Jqhp'  &
       eresj(28)=q4(3); eresj(29)=q4(1); eresj(30)=q4(2);
       ! ,' JKINt',' JKINn','JKINp',' JSO ','JCDIR',' JCEX','JDisn','JDisp'  &
       eresj(31)=rekin(3); eresj(32)=rekin(1); eresj(33)=rekin(2); eresj(34)=respopj;
       eresj(35)=recdipj; eresj(36)=recexpj; eresj(37)=Dispersion(1); eresj(38)=Dispersion(2);
       ! ,'JV2Mn','JV2Mp','JILST','JKIND','  JL '  &
       eresj(39)=v2min(1); eresj(40)=v2min(2)
       eresj(41)=iLST; eresj(42)=kindhfb; eresj(43)=iLpj;
       !  ,'JECMPAV','JECMPAV','JECMPAV'
       eresj(44)=ECMPAV(3); eresj(45)=ECMPAV(1); eresj(46)=ECMPAV(2);
       ! 'JA','JN',JZ'
       eresj(47)=Nint(xn(3)); eresj(48)=Nint(xn(1)); eresj(49)=Nint(xn(2));
       ! 'iter'
       eresj(50)=iiter
       ! nucleus with wrong asymptotic
       If(iasswrong(3).Ne.0) eresj(21)=-eresj(21)
       !
       ! WRITE to screen 'lout' and tape akzout.dat 'lfile'
       Do iw=lout,lfile
          Write(iw,*)
          Write(iw,'(a,9x,a,/)')            '  NB! From expectpj (PNP PAV RESULTS)'
          Write(iw,*)
          If(iLST1.Ne.0)  &
               Write(iw,'(a,6f15.6)') '  hfb decay const. ass ',ass
          Write(iw,'(a,8f15.6)') '  pairing: CpV0,CpV1,pwi... ',CpV0,CpV1,pwi
          Write(iw,'(a,a,a,i3)') '  forces:   ',skyrme,',  Gauge points:',ilpj
          If(keyblo(1).Ne.0)  &
               Write(iw,'(a,i4,a,f10.3)')  '  Blocked neutron block    ',  &
               bloblo(keyblo(1),1)
          If(keyblo(2).Ne.0)  &
               Write(iw,'(a,i4,a,f10.3)')  '  Blocked proton  block    ',  &
               bloblo(keyblo(2),2)
          Write(iw,*)
          Write(iw,'(/,28x,a,8x,a,9x,a)') ' neutrons ','protons','total'
          Write(iw,'(a,6f15.6)') '  Requested part.numbs.',tz,Sum(tz)
          Write(iw,'(a,6f15.6)') '  Projected part.numbs.',xn
          Write(iw,'(a,3f15.6)') '  Dispersion dN2 ......',Dispersion
          Write(iw,'(a,6f15.6)') '  b0, bz, bp ..........',b0,bz,bp
          Write(iw,*)
          Write(iw,'(a,6f15.6)') '  lambda (ala) ........',ala
          Write(iw,'(a,6f15.6)') '  Lambda (alast) ......',alast
          Write(iw,'(a,6f15.6)') '  delta(n,p) ..........',rdelta
          Write(iw,'(a,6f15.6)') '  pairing energy ......',repair
          Write(iw,*)
          Write(iw,'(a,6f15.6)') '  rms-radius ..........',rms
          Write(iw,'(a,15x,2f15.6)') '  charge-radius, r0 ...',rc,r00
          Write(iw,'(a,6f15.6)') '  deformation beta2 ...',def
          Write(iw,'(a,6f15.6)') '  quadrupole moment[b] ',q2
          Write(iw,'(a,6f15.6)') '  hexadecapole moment .',q4
          Write(iw,*)
          Write(iw,'(a,6f15.6)')     '  kinetic energy ......',rekin
          Write(iw,'(a,6f15.6)')     '  cmc-diagonal part ...',rekin/hb0*hbzero-rekin
          Write(iw,'(a,6f15.6)')     '  cmc-PAV .............',ECMPAV
          Write(iw,*)
          Write(iw,'(a,30x,6f15.6)') '  volume energy .......',revolpj
          Write(iw,'(a,30x,6f15.6)') '  surface energy ......',resurpj
          Write(iw,'(a,30x,6f15.6)') '  spin-orbit energy ...',respopj
          Write(iw,'(a,30x,6f15.6)') '  coulomb direct ......',recdipj
          Write(iw,'(a,30x,6f15.6)') '  coulomb exchange ....',recexpj
          Write(iw,'(a,30x,6f15.6)') '  tensor energy .......',retenspj
          Write(iw,*)
          Write(iw,'(a,30x,f15.6)')  '  Energy: ehfb(qp) ....',ehfb
          Write(iw,'(a,30x,f15.6)')  '  Energy: ehfb(can,pj).',rehfbcan
          Write(iw,'(a,30x,f15.6)')  '  ehfb(qp)-ehfb(can,pj)',ehfb-rehfbcan
          Write(iw,'(a,30x,f15.6)')  '  Epj-ehfb(can,pj) ....',depnp
          Write(iw,'(a,30x,6f15.6)') '  Energy: Epj=E(PAV) ..',retotpj
          Write(iw,*)
       End Do
       !
       ! APPEND the results to file 'thodef.dat'
       ! ieres=ieresu+ieresl+ieresj+ierebl from module definitions
       If(iappend.Ne.0) Then
          ierest=0
          ! charge buffers
          Do i=1,ieresj       !charge projected buffer
             ierest=ierest+1
             eres(ierest)=eresj(i)
          End Do
          Do i=1,ieresu       !charge unprojected buffer
             ierest=ierest+1
             eres(ierest)=eresu(i)
          End Do
          Do i=1,ieresl       !charge LN  buffer
             ierest=ierest+1
             eres(ierest)=eresl(i)
          End Do
          Do i=1,ieresbl      !charge Blocking buffer
             ierest=ierest+1
             eres(ierest)=eresbl(i)
          End Do
          If(ierest.Ne.ieres) Then
             ierror_flag=ierror_flag+1
             ierror_info(ierror_flag)='STOP: In expectpj: ierest wrong'
             Return
          End If
          If(Print_Screen) Then
             ! recording results
100        Continue                        ! complications are due to eagle_ornl
             If(iLST1.Le.0) Then
                Open (unit=iw2,file='hodef.dat',err=100,iostat=i,position='append')
                Write(iw2,'(3(1x,a,1x),160(1x,f14.6))') nucname,ereslbl,eres(1:ierest)
                Close(iw2)
             Else
                If (iasswrong(3).Eq.0) Then
                   Open (unit=iw1,file='thodef.dat',err=100,iostat=i,position='append')
                   Write(iw1,'(3(1x,a,1x),160(1x,f14.6))') nucname,ereslbl,eres(1:ierest)
                   Close(iw1)
                End If
             End If
          End If
       End If
    End If
    !
  End Subroutine expectpj
  !=======================================================================
  !
  !=======================================================================
  Subroutine densitpj
    !---------------------------------------------------------------------
    ! calculate local densities
    ! calculate local densities in mixed canonical (\for rho C, Y) and
    ! qp (for\tilde{\rho}) representation therefore E(can) \equiv E(qp)
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Use UNEDF
    Implicit None
    !
    Complex(pr) :: tpfiu1,tpfid1,v2ig,dig,sumsum
    Complex(pr), Allocatable :: ank1(:,:),pfiun1(:,:),pfidn1(:,:)
    Complex(pr), Allocatable :: pakapj(:),propj(:), pdjpj(:), ptaupj(:),pdropj(:)
    Complex(pr), Allocatable :: pszfipj(:),psfizpj(:),psrfipj(:),psfirpj(:)
    Complex(pr), Pointer:: ppjk(:),pcpj(:,:),prpj(:,:),pypj(:,:)
    Real(pr)    :: f,s,sd,su,sud,y,y2,sml2,cnzaa,cnraa,u,v2,tauin,xxx,yyy
    Real(pr)    :: aav,anik,anik2,qhla,qh1la,qhl1a,qla,qha,fi1r,fi1z,fi2d
    Real(pr)    :: xlam,xlam2,xlamy,xlamy2,xlap,xlap2,xlapy,xlapy2,xlampy
    Real(pr)    :: tfiu,tfid,tfiur,tfidr,tfiuz,tfidz,tfiud2,tfidd2,tpfiu2,tpfid2,TW_T
    Real(pr), Allocatable :: an2(:),ank2(:),pfiun2(:),pfidn2(:)
    Integer(ipr)    :: iw,nsa,nza,nra,nla,k,i,nd,il,ih,ihil,laplus,kkymu,n12
    Integer(ipr)    :: imen,ib,m,im,ig,it,j,jj,ja,jn,ILIHLI,k1,k2,kkk,kky,mu,kkkmu
    Integer(ipr)    :: k0(2),ky(2),kk(nqx),kyk(nqx)
    !
    Allocate(ank1(nqx,ilpj),pfiun1(ndx,ilpj),pfidn1(ndx,ilpj))
    Allocate(pfiun2(ndx),pfidn2(ndx),an2(nqx),ank2(nqx))
    Allocate(pakapj(ilnghl),propj(ilnghl), pdjpj(ilnghl), ptaupj(ilnghl),pdropj(ilnghl),  &
             pSZFIpj(ilnghl),pSFIZpj(ilnghl),pSRFIpj(ilnghl),pSFIRpj(ilnghl))
    !
    ! Projection grid points
    ! keypj=max(1,keypj); ilpj=keypj;  ilpj2=ilpj**2 !all
    ! when a value two*pi is used the results are precisely the same
    ! but the accuracy for even L is slow with increasing L.
    ! when 'pi' is used it gives regular and better convergence
    ! with respect to both, odd and even, L.
    ! Write(*,'(2x,a,i2,a,f12.8,a,f12.8)') 'point ig= ',i,' phi= ',yyy,' pi/2= ',pi/two
    xxx = pi/Real(ilpj,Kind=pr) ! equivalent to xxx = two*pi/Real(ilpj)
    Do i=1,ilpj
       yyy          = Real(i-1,Kind=pr)*xxx
       phypj(i)     = onei*yyy
       sinphy(i)    = onei*Sin(yyy)
       exp1iphy(i)  = Exp(onei*yyy)
       exp1iphym(i) = Exp(-onei*yyy)
       exp2iphy(i)  = Exp(two*onei*yyy)
       exp2iphym(i) = Exp(-two*onei*yyy)
    End Do
    !
    ! initialize parameters
    varmas = zero
    !
    Do it=itmin,itmax
       !
       ! zero for densities
       Do J=1,ilnghl
          pakapj(J)=zero; propj(J)=zero; pdjpj(J)=zero; ptaupj(J)=zero; pdropj(J)=zero;
       End Do
       Do J=1,ilnghl
          pszfipj(J)=zero; psfizpj(J)=zero; psrfipj(J)=zero; psfirpj(J)=zero;
       End Do
       !
       ! it-pointers
       prpj => rpj(:,:,it);  pcpj => cpj(:,:,it);
       pypj => ypj(:,:,it);  ppjk => pjk(:,it);
       !
       ! null for all pointers
       pypj=zero; prpj=zero; pcpj=zero; ppjk=one;
       !
       ! particle-init (kkk-even: 2 x number of pairs)
       kkk=npr(it); If(kkk.Ne.2*(kkk/2)) kkk=npr(it)-1
       ppjk(1:ilpj)=exp1iphym(1:ilpj)**kkk
       !
       ! start blocks
       k0(it)=0; ky(it)=0
       Do ib=1,nb
          nd=id(ib); im=ia(ib)
          If(Parity) Then
             LAPLUS=(ib+1)/2 !Yesp
          Else
             LAPLUS=ib       !Nop
          End If
          xlap=laplus; xlap2=xlap*xlap; xlam=xlap-one; xlam2=xlam*xlam
          !
          ! charge block can quantities
          m=ib+(it-1)*nbx; k1=ka(ib,it)+1; k2=ka(ib,it)+kd(ib,it); imen=0
          If(k1.Le.k2) Then
             ! below the pwi cut-off
             imen = nd
             !lcanon(ib,it)=lc
             Do k = 1,nd
                k0(it) = k0(it) + 1; kk(k)  = k0(it); kkk = k0(it)
                ky(it) = ky(it) + 1; kyk(k) = ky(it); kky = ky(it)
                aav    = vk(kkk,it)                                         ! v^2
                Do ig=1,ilpj
                   v2ig = exp2iphy(ig)*aav                                  ! gauged v^2
                   dig  = one - aav + v2ig                                  ! denominator
                   If(kkk.Ne.blocanon(it)) Then
                      ppjk(ig)  = ppjk(ig)*dig                              ! y(ig,it) <<<<<
                   End If
                   prpj(kkk,ig) = v2ig/dig                                  ! rho(mu,ig,it)
                   pcpj(kky,ig) = exp2iphy(ig)/dig                          ! c(mu,ig,it)
                   pypj(kky,ig) = exp1iphy(ig)/dig*onei*Sin(phypj(ig)/onei) ! sinphy(ig) !Y(mu,ig,it)
                End Do
             End Do
             ! At this point density (and related) are strictly equivalent in qp- and can-representation
             ! (up to 10^-14). Pairing density is not so strict (up to 10^-5) due to uv from v^2 but
             ! pairing density is taken directly in qp representation so both representations
             ! qp and can are strictly exact (up to 10^-14).
             j=0
             Do jj = 1,nd
                Do k = 1,nd
                   j=j+1; n12 = jj+(k-1)*nd;
                   an2(j)  = ddc(jj,kk(k),it)
                   ank2(j) = ak(n12,m)                 ! half \tilde{\rho} in q.p. basis
                   Do ig=1,ilpj
                      ank1(j,ig) = zero
                   End Do
                   Do mu=1,nd                          ! for half e^(-i\phy)*C(\phy)*\tilde{\rho} in q.p. basis
                      kkkmu = kk(mu); kkymu=kyk(mu)
                      Do ig=1,ilpj                     ! e^(-i\phy)*C in q.p. basis
                         ank1(j,ig) = ank1(j,ig) + ddc(jj,kkkmu,it)*ddc(k,kkkmu,it)*pcpj(kkymu,ig)*exp1iphym(ig)
                      End Do
                   End Do
                End Do
             End Do
          Else
             ! above the pwi cut-off (NB! Attention)
             ! here imem=0 and the contribution does
             ! not enter the densities but the Hamiltonian matrix
             ! used only in VAP regime
             ky(it)=ky(it)+1; kky = ky(it)
             Do ig=1,ilpj
                pcpj(kky,ig) = exp2iphy(ig)
                pypj(kky,ig) = exp1iphy(ig)*sinphy(ig)
             End Do
          End If
          !
          ! calculate the densities only below the PWI cutoff
          If (imen.Gt.0) Then
             ! gauss integration points
             Do il=1,ngl
                v2 = half/xl(il)
                Do ih=1,ngh
                   ihil = ih + (il-1)*ngh; ilihli=(ihil-1)*ilpj
                   !u = xh(ih); y = fli(ihil); y2=y*y
                   u = xh(ih); y = y_opt(ihil); y2=y*y
                   xlamy=xlam*y; xlamy2=xlam2*y2;
                   xlapy=xlap*y; xlapy2=xlap2*y2;
                   xlampy=xlamy+xlapy
                   !
                   ! initialize spin up/down funct
                   Do k=1,nd
                      fiu(k)=zero; fiuz(k)=zero; fiur(k)=zero; fiud2n(k)=zero; pfiun2(k)=zero;
                      fid(k)=zero; fidz(k)=zero; fidr(k)=zero; fidd2n(k)=zero; pfidn2(k)=zero;
                      Do ig=1,ilpj
                         pfiun1(k,ig)=zero; pfidn1(k,ig)=zero
                      End Do
                   End Do
                   !
                   ! scan over basis states
                   jn=0
                   Do i=1,nd
                      ja = i+im; nla = nl(ja); nra = nr(ja); nza = nz(ja); nsa = ns(ja);
                      sml2  = nla*nla; cnzaa = nza+nza+1; cnraa = nra+nra+nla+1
                      QHLA=QHLA_opt(JA,ihil); FI2D=FI2D_opt(JA,ihil)
                      FI1Z=FI1Z_opt(JA,ihil); FI1R=FI1R_opt(JA,ihil)

                      !qha   = qh(nza,ih); qla = ql(nra,nla,il); qhla = qha*qla
                      !qhl1a = qha*ql1(nra,nla,il)*v2; qh1la = qh1(nza,ih)*qla
                      !fi1z  = fp1(ihil)*qhla+fp2(ihil)*qh1la+fp3(ihil)*qhl1a
                      !fi1r  = fp4(ihil)*qhla+fp5(ihil)*qh1la+fp6(ihil)*qhl1a
                      !fi2d  = (fs1(ihil)*qh1la**2 + four*fs4(ihil)*qh1la*qhl1a        &
                      !      +  fs2(ihil)*qhl1a**2 + two*(fs5(ihil)*qh1la              &
                      !      +  fs6(ihil)*qhl1a)*qhla + ((u*u - cnzaa)*fs1(ihil)       &
                      !      +  (p14-cnraa*v2+sml2*v2*v2)*fs2(ihil)+fs3(ihil))*qhla**2 &
                      !      -  two*(fi1r**2+fi1z**2))/(two*qhla)
                      !
                      ! wave function(spin:up,down; grad:r,z,d2)
                      If (nsa.Gt.0) Then
                         Do k=1,nd
                            jn = jn+1; anik = an2(jn); anik2 = ank2(jn)
                            Do ig=1,ilpj
                               pfiun1(k,ig) = pfiun1(k,ig) + ank1(jn,ig)*qhla
                            End Do
                            pfiun2(k) = pfiun2(k) + anik2*qhla
                            fiu(k)    = fiu(k)    + anik*qhla
                            fiur(k)   = fiur(k)   + anik*fi1r
                            fiuz(k)   = fiuz(k)   + anik*fi1z
                            fiud2n(k) = fiud2n(k) + anik*fi2d
                            !
                         End Do
                      Else
                         Do k=1,nd
                            jn = jn+1; anik = an2(jn); anik2 = ank2(jn)
                            Do ig=1,ilpj
                               pfidn1(k,ig) = pfidn1(k,ig) + ank1(jn,ig)*qhla
                            End Do
                            pfidn2(k) = pfidn2(k) + anik2*qhla
                            fid(k)    = fid(k)    + anik*qhla
                            fidr(k)   = fidr(k)   + anik*fi1r
                            fidz(k)   = fidz(k)   + anik*fi1z
                            fidd2n(k) = fidd2n(k) + anik*fi2d
                            !
                         End Do
                      End If
                   End Do ! i
                   !
                   ! calculate densities
                   Do k=1,nd
                      kkk =kk(k)
                      tfiu=fiu(k); tfiuz=fiuz(k); tfiur=fiur(k); tfiud2=fiud2n(k); tpfiu2=pfiun2(k);
                      tfid=fid(k); tfidz=fidz(k); tfidr=fidr(k); tfidd2=fidd2n(k); tpfid2=pfidn2(k);
                      Do ig=1,ilpj
                         I=ig+ilihli; v2ig=prpj(kkk,ig); tpfiu1=pfiun1(k,ig); tpfid1=pfidn1(k,ig)
                         !
                         pakapj(I)  = pakapj(I)  +  (tpfiu1*tpfiu2+tpfid1*tpfid2)
                         propj(I)   = propj(I)   +  (tfiu**2+tfid**2)*v2ig
                         pdjpj(I)   = pdjpj(I)   +  (tfiur*tfidz-tfidr*tfiuz+xlamy*tfiu*(tfiur-tfidz) &
                                                 -   xlapy*tfid*(tfidr+tfiuz))*v2ig
                         TW_T=(tfiur**2+tfidr**2+tfiuz**2+tfidz**2)
                         tauin      = (xlamy2*tfiu**2+xlapy2*tfid**2+TW_T)
                         ptaupj(I)  = ptaupj(I)  +   tauin*v2ig
                         pdropj(I)  = pdropj(I)  +  (TW_T + tfiu*tfiud2 + tfid*tfidd2)*v2ig
                         psrfipj(I) = psrfipj(I) + (tfiur*tfid - tfidr*tfiu)*v2ig
                         psfirpj(I) = psfirpj(I) + (tfiu*tfid*xlampy)*v2ig
                         psfizpj(I) = psfizpj(I) + (xlamy*tfiu**2 - xlapy*tfid**2)*v2ig
                         pszfipj(I) = pszfipj(I) + (tfiuz*tfid - tfidz*tfiu)*v2ig
                         !
                      End Do !ig
                   End Do !k
                End Do !ih
             End Do !il
          End If
       End Do !ib
       !
       ! normalized pjk
       sumsum = Sum(ppjk(1:ilpj)); ppjk(1:ilpj) = ppjk(1:ilpj)/sumsum
       !
       ! Y minus second term of Y
       Do k=1,ky(it)
          sumsum = Sum(ppjk(1:ilpj)*pypj(k,1:ilpj))
          pypj(k,1:ilpj) = pypj(k,1:ilpj) - sumsum
       End Do
       !
       ! norm of the projected/unprojected density
       s = zero; sd = zero; su = zero; sud = zero;
       Do ihil=1,nghl
          ilihli=(ihil-1)*ilpj
          Do ig=1,ilpj
             I=ig+ilihli
             s=s+Real(two*propj(I)*ppjk(ig)); sd=sd+Real(four*pdropj(I)*ppjk(ig))
          End Do
          I=1+ilihli
          su=su+Real(two*propj(I)); sud=sud+Real(four*pdropj(I))
       End Do
       !
       ! print unprojected normalization
       Do iw=lout,lfile
          Write(iw,'(2(a,2(2x,D15.8)),(a,D15.8),a,i3)') &
               '   pj/unpj  s= ',s,su,'   pj/unpj sd= ',sd,sud,' ala1= ',ala1(it),' inner= ',inner(it)
       End Do
       varmas = varmas + su
       varmasNZ(it) = su; pjmassNZ(it) = s
       !
       s = Real(npr(it),Kind=pr)/s; dnfactor(it) = s; drhoi(it) = sd
       !
       Do ihil = 1,nghl
          ilihli=(ihil-1)*ilpj
          ! wdcor moves out the int.weight and multiply by the jacobian
          f = two*wdcori(ihil)
          Do ig=1,ilpj
             I=ig+ilihli
             ropj (ihil,ig,it)  = f*propj (I)
             taupj(ihil,ig,it)  = f*ptaupj(I)
             dropj(ihil,ig,it)  = f*pdropj(I)*two
             djpj (ihil,ig,it)  = f*pdjpj (I)*two
             akapj(ihil,ig,it)  = f*pakapj(I)*half
             SRFIpj(ihil,ig,it) = f*psrfipj(I)
             SFIRpj(ihil,ig,it) = f*psfirpj(I)
             SFIZpj(ihil,ig,it) = f*psfizpj(I)
             SZFIpj(ihil,ig,it) = f*pszfipj(I)
          End Do !ig
       End Do !ihil
       !
    End Do !it
    !
    dnfactor(3)=dnfactor(1)+dnfactor(2)
    !
    Deallocate(ank1,pfiun1,pfidn1)
    Deallocate(pakapj,propj, pdjpj, ptaupj,pdropj,pSZFIpj,pSFIZpj,pSRFIpj,pSFIRpj)
    !
    Call coulompj !complex coulomb fields
    !
  End Subroutine densitpj
  !=======================================================================
  !
  !=======================================================================
  Subroutine coulompj
    !---------------------------------------------------------------------
    ! Coulom-field (direct part)
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: i,j,k
    Real(pr) :: zd2,y1,y2,xx1,s1,vik,f,r,r1,rr2,z,z1,zd1,t
    Real(pr) :: bb,r2,r12,rrr,rz1,rz2,rrz1,rrz2,xx
    Real(pr) :: bb1=3.5156229_pr,g1=0.39894228_pr,g7=0.02635537_pr,  &
                bb2=3.0899424_pr,g2=0.01328592_pr,g8=0.01647633_pr,  &
                bb3=1.2067492_pr,g3=0.00225319_pr,g9=0.00392377_pr,  &
                bb4=0.2659732_pr,g4=0.00157565_pr,bbxx=3.750_pr,     &
                bb5=0.0360768_pr,g5=0.00916281_pr,                   &
                bb6=0.0045813_pr,g6=0.02057706_pr
    If(icacoupj.Eq.0) Then
       icacoupj = 1; bb = Max(bp,bz)**4; f = chargee2/Sqrt(pi);  ! f=e^2/Sqrt(pi)
       Do i = 1,nghl
          r = fl(i); z = fh(i); r2 = r*r
          Do k = 1,i
             r1   = fl(k);       z1   = fh(k);      r12 = r1*r1
             rrr  = two*r*r1;    rr2  = (r - r1)**2
             zd1  = (z - z1)**2; zd2  = (z + z1)**2
             rz1  = r2+r12+zd1;  rz2  = r2+r12+zd2
             rrz1 = rr2+zd1;     rrz2 = rr2+zd2
             !
             xx1=zero
             Do j=1,nleg
                xx=Sqrt(one-xleg(j)**2); y1=(xleg(j)/(bb*xx))**2; s1=y1*rrr
                If(s1.Le.bbxx) Then
                   t=(s1/bbxx)**2; y2=one+t*(bb1+t*(bb2+t*(bb3+t*(bb4+t*(bb5+t*bb6)))))
                   y2=y2*(Exp(-rz1*y1)+Exp(-rz2*y1))
                Else
                   t=(bbxx/s1); y2=g1+t*(g2+t*(g3+t*(-g4+t*(g5+t*(-g6+t*(g7+t*(-g8+t*g9)))))))
                   y2=y2/Sqrt(s1)*(Exp(-rrz1*y1)+Exp(-rrz2*y1))
                End If
                xx1 = xx1 + wleg(j)*y2/(bb*xx**3)
             End Do
             vik=f*xx1; vc(i,k)=vik*wdcor(k); vc(k,i)=vik*wdcor(i)  !wdcor=pi*wh*wl*bz*bp*bp/fd
          End Do  !k
       End Do  !i
    End If
    ! calculation of the coulomb field
    coupj = zero
    Do i = 1,nghl
       Do k=1,ilpj
          coupj(:,k) = coupj(:,k) + vc(:,i)*ropj(i,k,2)
       End Do
    End Do
  End Subroutine coulompj
  !=======================================================================
  !
  !=======================================================================
  Subroutine broyden_min(N,vout,vin,alpha,si,iter,M,bbroyden)
    !---------------------------------------------------------------------
    ! Modified Broyden's method: D.D.Johnson, PRB 38, 12807 (1988)
    ! Adopted from: (C) 2001 PWSCF group
    ! Input :
    !  N      dimension of arrays vin,vout
    !  vin    outpu at previous iteration
    !  vout   output at current iteration
    !  alpha  mixing factor (0 < alpha <= 1)
    !  iter   current iteration number
    !  M      number of iterations in Broyden history
    !  M=0    Linear mixing
    ! Output:
    !  si     MaxVal(|vout-vin|)
    !  vin    Broyden/Linear mixing result
    !  vout   vout-vin
    !  bbroyden='B' Broyden mixing, curvature>0
    !  bbroyden='L' Linear mixing,  curvature<0
    !---------------------------------------------------------------------
    Use HFBTHO_utilities, Only: pr,ipr
    Use HFBTHO, Only: ierror_flag,ierror_info
    Implicit None
    Integer(ipr),     Intent(In)    :: N,iter,M
    Real(pr),        Intent(In)     :: alpha
    Real(pr),        Intent(Out)    :: si
    Character(1),      Intent(Out)  :: bbroyden
    Real(pr),        Intent(InOut)  :: vout(N),vin(N)
    Integer(ipr)                    :: i,j,iter_used,ipos,inext,info
    Integer(ipr), Allocatable, Save :: iwork(:)
    Real(pr),    Allocatable, Save  :: beta(:,:),work(:)
    Real(pr),    Allocatable, Save  :: df(:,:),dv(:,:),curv(:)
    Real(pr),                 Save  :: w0
    Real(pr)                        :: DDOT,DNRM2,normi,gamma,curvature,sf
    !
    sf=-1.0_pr; Call DAXPY(N,sf,vin,1,vout,1)
    si=Maxval(Abs(vout))
    ! Linear mixing
    If(M.Eq.0.Or.iter.Eq.0) Then
       bbroyden='L'; Call DAXPY(N,alpha,vout,1,vin,1)
       !If(iter.Eq.0) Write(6,*) '  Linear mixing (alpha) : ',alpha
       Return
    End If
    ! Broyden mixing
    iter_used=Min(iter-1,M)
    ipos=iter-1-((iter-2)/M)*M
    inext=iter-((iter-1)/M)*M
    If (iter.Eq.1) Then
       w0=0.010_pr
       If(Allocated(df)) Deallocate(curv,df,dv,beta,work,iwork)
       Allocate(curv(N),df(N,M),dv(N,M),beta(M,M),work(M),iwork(M))
    Else
       df(:,ipos)=vout(:)-df(:,ipos); dv(:,ipos)=vin(:)-dv(:,ipos)
       Normi=1.0_pr/Sqrt((DNRM2(N,df(1,ipos),1))**2)
       Call DSCAL(N,Normi,df(1,ipos),1)
       Call DSCAL(N,Normi,dv(1,ipos),1)
    End If
    Do i=1,iter_used
       Do j=i+1,iter_used
          beta(i,j)=DDOT(N,df(1, j),1,df(1,i),1)
       End Do
       beta(i,i)=1.0_pr+w0*w0
    End Do
#if(SWITCH_ESSL==0)
      Call DSYTRF('U', iter_used, beta, M, iwork, work, M, info)
#else
      Call DPOTRF('U', iter_used, beta, M, info)
#endif
    If(info.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='STOP: In Broyden: info at DSYTRF '
       Return
    End If
#if(SWITCH_ESSL==0)
      Call DSYTRI('U', iter_used, beta, M, iwork, work, info)
#else
      Call DPOTRI('U', iter_used, beta, M, info)
#endif
    If(info.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='STOP: In Broyden: info at DSYTRI '
       Return
    End If
    Do i=1,iter_used
       Do j=i+1,iter_used
          beta(j,i)=beta(i,j)
       End Do
       work(i)=DDOT(N,df(1,i),1,vout,1)
    End Do
    curv=alpha*vout
    Do i=1,iter_used
       gamma=0.0_pr
       Do j=1,iter_used
          gamma=gamma+beta(j,i)*work(j)
       End Do
       curv=curv-gamma*(dv(:,i)+alpha*df(:,i))
    End Do
    Call DCOPY(N,vout,1,df(1,inext),1)
    Call DCOPY(N,vin,1,dv(1,inext),1)
    curvature=DDOT(N,vout,1,curv,1)
    If(curvature.Gt.-1.0_pr) Then
       bbroyden='B'; sf=+1.0_pr; Call DAXPY(N,sf,curv,1,vin,1)
    Else
       bbroyden='L'; sf=alpha*0.50_pr; Call DAXPY(N,sf,vout,1,vin,1)
    End If
  End Subroutine broyden_min
  !=======================================================================
  !
  !=======================================================================
  Subroutine expect(lpr)
    !---------------------------------------------------------------------
    ! calculates expectation values (xn is the particle number)
    ! at lpr=.true. also calculates PAV corrections
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Logical :: lpr
    Integer(ipr) :: i,it,ihli,iw
    Real(pr) :: ekt(3),xn(3),q4(3),def(3),bet2(3),het4(3),econst
    Real(pr) :: z,zz,rrr,p2,p3,p4,row,r212,r222,rc
    Real(pr) :: eso,ecodi,ecoex,rn,rp,rnp1,rnp2,rt,whl,tnt,tpt,tt
    Real(pr) :: dn,dp,dt,akn,akp,akn2,akp2,adn,adp,evol,esurf,ecoul
    Real(pr) :: etens,dd1n,dd1p,rt1,tt1,dt1,djn,djp,djt,djt1
    Real(pr) :: RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1,DJ_0,DJ_1,J2_0,J2_1
    Real(pr) :: SZFIN,SFIZN,SRFIN,SFIRN,SZFIP,SFIZP,SRFIP,SFIRP
    Real(pr) :: SZFI_0,SFIZ_0,SRFI_0,SFIR_0,SZFI_1,SFIZ_1,SRFI_1,SFIR_1
    Real(pr) :: SNABLARN,SNABLAZN,SNABLARP,SNABLAZP
    Real(pr) :: SNABLAR_0,SNABLAZ_0,SNABLAR_1,SNABLAZ_1
    Real(pr) :: xn1,xn2,rms1,rms2,q21,q22,q41,q42,EKIN_N,EKIN_P,ept1,ept2,del1,del2
    Real(pr) :: rsa,rsa0,rsa0A,rps,rns,rsa1,rsa10,rsa12,rsa0An,rsa0As
    Real(pr) :: ESURF_rho_DELTA_rho,ESURF_NABLA_rho_NABLA_rho,ESO_rho_NABLA_J,ESO_NABLA_rho_J
    Real(pr) :: EVOL_rho_tau,EVOL_rho_rho,EExtra,E_HARTREE_DIR,tempE_Crho0,tempREARR
    Real(pr) :: E_EXT_FIELD
    Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:)
    Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
    !------------------------------------------------
    ! Part called during iterations (lpr=F)
    !------------------------------------------------
    !
    Call DENSIT
    If(ierror_flag.Ne.0) Return
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('expect',0)
    !
    !------------------------------------------------
    ! zero energy variables
    !------------------------------------------------
    EKIN_N=zero;              EKIN_P=zero;
    EVOL_rho_tau=zero;        EVOL_rho_rho=zero;
    ESURF_rho_DELTA_rho=zero; ESURF_NABLA_rho_NABLA_rho=zero;
    ESO_rho_NABLA_J=zero;     ESO_NABLA_rho_J=zero; E_HARTREE_DIR=zero
    ept1=zero; ept2=zero;     del1=zero; del2=zero;
    ecodi=zero; ecoex=zero; etens=zero
    EExtra=zero ; E_EXT_FIELD = zero ;
    xn1=zero; xn2=zero; rms1=zero; rms2=zero
    q21=zero; q22=zero; q41=zero; q42=zero
    tempE_Crho0=zero; tempREARR=zero
    DEROT=zero; SQUJ=zero; CRAN=zero; ERIGHFB=zero
    !------------------------------------------------
    ! zero optimization variables
    !------------------------------------------------
    If(DO_FITT) Then
       efit_0=zero; efitV0=zero; dfitV0=zero
       efit_rhorho=zero; efit_rhorhoD=zero;
       efit_rhotau=zero; efit_rhoDrho=zero;
       efit_rhonablaJ=zero; efit_JJ=zero;
    End If
    !------------------------------------------------
    ! Integration in coordinate space
    !------------------------------------------------
    Do ihli=1,nghl
       whl=wdcor(ihli)
       !------------------------------------------------
       ! np-representation
       !------------------------------------------------
       rn=ro(ihli,1);      rp=ro(ihli,2); rnp2=rn**2+rp**2; rnp1=rn-rp
       tnt=tau(ihli,1);    tpt=tau(ihli,2);
       dn=dro(ihli,1);     dp=dro(ihli,2);
       djn=dj(ihli,1);     djp=dj(ihli,2);
       akn=aka(ihli,1);    akp=aka(ihli,2)
       akn2=akn*akn;       akp2=akp*akp
       adn=akn*rn;         adp=akp*rp
       SFIZN=SFIZ(IHLI,1); SFIZP=SFIZ(IHLI,2);
       SFIRN=SFIR(IHLI,1); SFIRP=SFIR(IHLI,2);
       SZFIN=SZFI(IHLI,1); SZFIP=SZFI(IHLI,2);
       SRFIN=SRFI(IHLI,1); SRFIP=SRFI(IHLI,2);
       SNABLARN=NABLAR(IHLI,1); SNABLARP=NABLAR(IHLI,2);
       SNABLAZN=NABLAZ(IHLI,1); SNABLAZP=NABLAZ(IHLI,2);
       !------------------------------------------------
       ! t-representation
       !------------------------------------------------
       RHO_0=rn+rp;        RHO_1=rn-rp;
       TAU_0=tnt+tpt;      TAU_1=tnt-tpt;
       DRHO_0=dn+dp;       DRHO_1=dn-dp;
       DJ_0=djn+djp;       DJ_1=djn-djp;
       SFIZ_0=SFIZN+SFIZP; SFIZ_1=SFIZN-SFIZP;
       SFIR_0=SFIRN+SFIRP; SFIR_1=SFIRN-SFIRP;
       SZFI_0=SZFIN+SZFIP; SZFI_1=SZFIN-SZFIP;
       SRFI_0=SRFIN+SRFIP; SRFI_1=SRFIN-SRFIP;
       SNABLAR_0=SNABLARN+SNABLARP; SNABLAR_1=SNABLARN-SNABLARP;
       SNABLAZ_0=SNABLAZN+SNABLAZP; SNABLAZ_1=SNABLAZN-SNABLAZP;
       J2_0=SFIZ_0**2+SFIR_0**2+SZFI_0**2+SRFI_0**2
       J2_1=SFIZ_1**2+SFIR_1**2+SZFI_1**2+SRFI_1**2
       !
       Call calculate_U_parameters(RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1,  &
            (SNABLAR_0**2+SNABLAZ_0**2),(SNABLAR_1**2+SNABLAZ_1**2))
       !------------------------------------------------
       ! rms and deformations
       !------------------------------------------------
       z=fh(ihli); zz=z*z; rrr=zz+fl(ihli)**2
       p2=p32*zz-half*rrr; p3=p53*z*p2-p23*rrr*z; p4=p74*z*p3-p34*rrr*p2
       row=whl*rn; xn1=xn1+row; rms1=rms1+row*rrr; q21=q21+row*p2; q41=q41+row*p4
       row=whl*rp; xn2=xn2+row; rms2=rms2+row*rrr; q22=q22+row*p2; q42=q42+row*p4
       !------------------------------------------------
       ! PH energies
       !------------------------------------------------
       EKIN_N=EKIN_N+hb0*(TAU_0+TAU_1)*HALF*whl*facECM                          ! kinetic, n
       EKIN_P=EKIN_P+hb0*(TAU_0-TAU_1)*HALF*whl*facECM                          ! kinetic, p
       EVOL_rho_tau=EVOL_rho_tau+(Urhotau(0,0)*RHO_0*TAU_0  &                   ! volume rho tau
            +Urhotau(1,0)*RHO_1*TAU_1+Urhotau(2,0)*RHO_0*TAU_1  &
            +Urhotau(3,0)*RHO_1*TAU_0 )*whl
       EVOL_rho_rho=EVOL_rho_rho+(Urhorho(0,0)*RHO_0**2  &                      ! volume density dependent
            +Urhorho(1,0)*RHO_1**2+(Urhorho(3,0)+Urhorho(2,0))*RHO_0*RHO_1)*whl
       ESURF_rho_DELTA_rho =ESURF_rho_DELTA_rho+(UrhoDrho(0,0)*RHO_0*DRHO_0  &  ! surface: rho delta rho
            +UrhoDrho(1,0)*RHO_1*DRHO_1+UrhoDrho(2,0)*RHO_0*DRHO_1  &
            +UrhoDrho(3,0)*RHO_1*DRHO_0 )*whl
       ESURF_NABLA_rho_NABLA_rho=ESURF_NABLA_rho_NABLA_rho  &                   ! surface: (nabla rho)**2
            +(Unablarho(0,0)*(SNABLAR_0*SNABLAR_0+SNABLAZ_0*SNABLAZ_0)  &
            +Unablarho(1,0)*(SNABLAR_1*SNABLAR_1+SNABLAZ_1*SNABLAZ_1)  &
            +(Unablarho(3,0)+Unablarho(2,0))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1) )*whl
       ESO_rho_NABLA_J=ESO_rho_NABLA_J+(UrhonablaJ(0,0)*RHO_0*DJ_0  &           ! spin-orbit rho Nabla . J
            +UrhonablaJ(1,0)*RHO_1*DJ_1+UrhonablaJ(2,0)*RHO_0*DJ_1  &
            +UrhonablaJ(3,0)*RHO_1*DJ_0 )*whl
       ESO_NABLA_rho_J=ESO_NABLA_rho_J  &
            +(UJnablarho(0,0)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))  &  ! spin-orbit J . Nabla rho
            +UJnablarho(1,0)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))  &
            +UJnablarho(2,0)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))  &
            +UJnablarho(3,0)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1)) )*whl
       ETENS=ETENS+(UJJ(0,0)*J2_0+UJJ(1,0)*J2_1  &                              ! tensor J^2
            +(UJJ(3,0)+UJJ(2,0))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1) )*whl
       EExtra=EExtra+(UEnonstdr(0)+UEnonstdr(1))*whl                            ! extra field if needed
       E_EXT_FIELD=E_EXT_FIELD + ( Vexternal(0,zero,fl(ihli),z)*RHO_0 &         ! external field
            +Vexternal(1,zero,fl(ihli),z)*RHO_1 )*whl
       !------------------------------------------------
       ! Coulomb & Hartree
       !------------------------------------------------
       If (icou.Ge.1) ecodi=ecodi+half*cou(ihli)*rp*whl
       If (icou.Eq.2) ecoex=ecoex-CExPar*cex*rp**p43*whl
       E_HARTREE_DIR=E_HARTREE_DIR +half*vDHartree(ihli,1)*RHO_0*whl+half*vDHartree(ihli,2)*RHO_1*whl
       ! just for printing
       tempE_Crho0=tempE_Crho0+RHO_0**2*whl
       tempREARR=tempREARR+(Cdrho(0)*RHO_0**2+Cdrho(1)*RHO_1**2)*RHO_0**sigma*whl
       !------------------------------------------------
       ! pairing energy and delta
       !------------------------------------------------
       rsa0=(RHO_0/rho_c)
       dd1n=CpV0(0)*(ONE-rsa0*CpV1(0))*whl
       dd1p=CpV0(1)*(ONE-rsa0*CpV1(1))*whl
       ept1=ept1+dd1n*akn2; del1=del1-dd1n*adn
       ept2=ept2+dd1p*akp2; del2=del2-dd1p*adp
       !------------------------------------------------
       ! optimization quantities
       !------------------------------------------------
       If(DO_FITT) Then
          efitV0(0)=efitV0(0)+(ONE-rsa0*CpV1(0))*akn2*whl
          efitV0(1)=efitV0(1)+(ONE-rsa0*CpV1(1))*akp2*whl
          dfitV0(0)=dfitV0(0)-(ONE-rsa0*CpV1(0))*adn*whl
          dfitV0(1)=dfitV0(1)-(ONE-rsa0*CpV1(1))*adp*whl
          !
          efit_rhotau(0)=efit_rhotau(0)+RHO_0*TAU_0*whl              ! rho tau
          efit_rhotau(1)=efit_rhotau(1)+RHO_1*TAU_1*whl              ! rho tau
          efit_rhorho(0)=efit_rhorho(0)+RHO_0**2*whl                 ! rho^2
          efit_rhorho(1)=efit_rhorho(1)+RHO_1**2*whl                 ! rho^2
          efit_rhorhoD(0)=efit_rhorhoD(0)+RHO_0**sigma*RHO_0**2*whl  ! rho^2
          efit_rhorhoD(1)=efit_rhorhoD(1)+RHO_0**sigma*RHO_1**2*whl  ! rho^2
          efit_rhoDrho(0)=efit_rhoDrho(0)+RHO_0*DRHO_0*whl           ! rho Delta rho
          efit_rhoDrho(1)=efit_rhoDrho(1)+RHO_1*DRHO_1*whl           ! rho Delta rho
          efit_rhonablaJ(0)=efit_rhonablaJ(0)+RHO_0*DJ_0*whl         ! rho nabla J J
          efit_rhonablaJ(1)=efit_rhonablaJ(1)+RHO_1*DJ_1*whl         ! rho nabla J J
          efit_JJ(0)=efit_JJ(0)+J2_0*whl                             ! J.J
          efit_JJ(1)=efit_JJ(1)+J2_1*whl                             ! J.J
       End If
    End Do !ihli
    !------------------------------------------------
    ! after the integration
    !------------------------------------------------
    xn(1)=xn1;         xn(2)=xn2;         xn(3)=xn1+xn2;
    rms(1)=rms1;       rms(2)=rms2
    q2(1)=q21;         q2(2)=q22;
    q4(1)=q41;         q4(2)=q42
    ekt(1)=EKIN_N;     ekt(2)=EKIN_P;     ekt(3)=ekt(1)+ekt(2)
    ept(1)=ept1;       ept(2)=ept2;       ept(3)=ept(1)+ept(2)
    del(1)=del1/xn(1); del(2)=del2/xn(2);
    !
    EVOL=EVOL_rho_tau+EVOL_rho_rho+E_HARTREE_DIR
    esurf=ESURF_rho_DELTA_rho+ESURF_NABLA_rho_NABLA_rho
    ESO=ESO_rho_NABLA_J+ESO_NABLA_rho_J
    ecoul=ecodi+ecoex
    etot=ekt(3)+evol+esurf+eso+ecoul+ept(3)+ETENS+EExtra+E_EXT_FIELD
    ehfb=etot
    entropy(3)=entropy(1)+entropy(2)
    !------------------------------------------------
    ! rms and deformations
    !------------------------------------------------
    Do it=itmin,itmax
       rms(it)=Sqrt(rms(it)/xn(it))
       q2(it)=two*q2(it)       !Qnp=<2r^2P_2(teta)>=<2z^2-x^2-y^2>
       q4(it)=ffdef4*q4(it)    !Hn=8r^4P_4(teta)=8z^4-24z^2(x^2+y^2)+3(x^2+y^2)^2
       def(it)=Sqrt(pi/5.0_pr)*q2(it)/(rms(it)**2*xn(it))
    End Do
    r212=rms(1)**2; r222=rms(2)**2
    rms(3)=Sqrt((xn(1)*r212+xn(2)*r222)/amas)
    q2(3)=q2(1)+q2(2)          ! quadrupole moment
    q4(3)=q4(1)+q4(2)          ! hexadecapole moment
    def(3)=Sqrt(pi/5.0_pr)*q2(3)/(rms(3)**2*amas) !deformation
    bet=def(3)
    !bet=ffdef6*q2(3)/(amas*r02)  ! bet=Q2*Sqrt(5 Pi)/(3A x^2);  x=r0 A^(1/3)
    !------------------------------------------------
    ! Lipkin-Nogami energy
    !------------------------------------------------
    If(kindhfb.Lt.0) Then
       Call tracesln
       If(ierror_flag.Ne.0) Return
       etot=etot+etr(3)
    End If
    !------------------------------------------------
    ! optimization quantities
    !------------------------------------------------
    If(DO_FITT) Then
       efV_0=0.0_pr
       If(kindhfb.Lt.0) Then
          efV_0(0)=ala2(1)
          efV_0(1)=ala2(2)
       End If
       dfitV0(0)=dfitV0(0)/xn(1)
       dfitV0(1)=dfitV0(1)/xn(2)
       efit_0=etot-efitV0(0)*CpV0(0)-efitV0(1)*CpV0(1)  &
            -efit_rhotau(0)*Ctau(0)-efit_rhotau(1)*Ctau(1)  &
            -efit_rhorho(0)*Crho(0)-efit_rhorho(1)*Crho(1)  &
            -efit_rhorhoD(0)*Cdrho(0)-efit_rhorhoD(1)*Cdrho(1)  &
            -efit_rhoDrho(0)*CrDr(0)-efit_rhoDrho(1)*CrDr(1)  &
            -efit_rhonablaJ(0)*CrdJ(0)-efit_rhonablaJ(1)*CrdJ(1)  &
            -efit_JJ(0)*CJ(0)-efit_JJ(1)*CJ(1)
    End If
    !------------------------------------------------
    ! expectation values of multipole moments
    !------------------------------------------------
    Call moments_computeValue()
    !------------------------------------------------
    ! debug
    !------------------------------------------------
    If(Print_Screen.And.IDEBUG.Gt.10) Then
       Write(*,'(4(a12,g13.6))')  &
            ' Tn=     ',ekt(1),           ' Tp=     ',ekt(2), &
            ' EPn=    ',ept(1),           ' EPp=    ',ept(2),  &
            ' EVOL=   ',EVOL,             ' Esurf=  ',esurf,  &
            ' NrNr=   ',ESURF_NABLA_rho_NABLA_rho,' rDr=    ',ESURF_rho_DELTA_rho,  &
            ' Etens=  ',ETENS,            ' Eso=   ',eso,  &
            ' rNJ=    ',ESO_rho_NABLA_J,  ' NrJ=   ',ESO_NABLA_rho_J,  &
            ' ECd=    ',ecodi,            ' ECex=  ',ecoex, &
            ' EHd=    ',E_HARTREE_DIR,    ' Ir0^2= ',tempE_Crho0, &
            ' Eextra= ',EExtra,           ' Ext.Fl= ',E_EXT_FIELD, &
            ' Etot=  ',etot
       If(DO_FITT) Then
          Write(*,'(4(a12,g13.6))')
          Write(*,'(4(a12,g13.6))')  &
               ' efrr0= ',efit_rhorho(0),     ' efrr1= ',efit_rhorho(1), &
               ' efrrD0=   ',efit_rhorhoD(0),          ' efrr1D=  ',efit_rhorhoD(1),  &
               ' efrt0= ',efit_rhotau(0),     ' efrt1= ',efit_rhotau(1), &
               ' efrDr0=   ',efit_rhoDrho(0),          ' efrDr1=  ',efit_rhoDrho(1),  &
               ' efrDj0=',efit_rhonablaJ(0),  ' efrDj1=',efit_rhonablaJ(1), &
               ' efjj0=    ',efit_JJ(0),               ' efjj1=   ',efit_JJ(1),  &
               ' efV0_0=',efitV0(0),          ' efV0_1=',efitV0(1), &
               ' dfV0_0=   ',dfitV0(0),                ' dfV0_1=  ',dfitV0(1),  &
               ' efV0=  ',efV_0(0),           ' efV_1= ',efV_0(1), &
               ' ef0=      ',efit_0,                   ' etot=    ',etot
       End If
    End If
    !------------------------------------------------
    ! Part called at the very end only (lpr=T)
    !------------------------------------------------
    If (lpr) Then
       !------------------------------------------------
       ! other definitions of deformations  (ffdef6=Sqrt(5.0_pr*pi)/3.0_pr)
       !------------------------------------------------
       bet2(1)=ffdef6*q2(1)/(xn(1)*r02) ! beta_n=Qn*Sqrt(5 Pi)/(3N x^2)
       bet2(2)=ffdef6*q2(2)/(xn(2)*r02) ! x=r0 A^(1/3)
       bet2(3)=ffdef6*q2(3)/(amas*r02)
       het4(1)=ffdef7*q4(1)/(xn(1)*r04)
       het4(2)=ffdef7*q4(2)/(xn(2)*r04)
       het4(3)=ffdef7*q4(3)/(amas*r04)
       rc=Sqrt(r222+0.640_pr)
       ! transitions to barn,barn^2,barn^4
       Do i=1,3
          q2(i)=q2(i)/100.0_pr; q4(i)=q4(i)/10000.0_pr
       End Do
       !------------------------------------------------
       ! STORE to unprojected buffer 'eresu'
       !------------------------------------------------
       ! ieresu=50 from module definitions
       ! ,'UEtot','Ubett','Ubetn','Ubetp',' UQt ',' UQn ',' UQp '  &
       eresu(1)=etot; eresu(2)=def(3); eresu(3)=def(1); eresu(4)=def(2);
       eresu(5)=q2(3); eresu(6)=q2(1); eresu(7)=q2(2);
       ! ,' Uln ',' Ulp ',' UpEn',' UpEp',' UpDn',' UpDp',' UAsn',' UAsp'  &
       eresu(8)=alast(1); eresu(9)=alast(2); eresu(10)=ept(1); eresu(11)=ept(2);
       eresu(12)=del(1); eresu(13)=del(2); eresu(14)=ass(1); eresu(15)=ass(2);
       ! ,' Urt ',' Urn ',' Urp ',' Urc ',' Uht ',' Uhn ',' Uhp '  &
       eresu(16)=rms(3); eresu(17)=rms(1); eresu(18)=rms(2); eresu(19)=rc;
       eresu(20)=het4(3); eresu(21)=het4(1); eresu(22)=het4(2);
       ! ,' Uqht',' Uqhn',' Uqhp'  &
       eresu(23)=q4(3); eresu(24)=q4(1); eresu(25)=q4(2);
       ! ,'UKINT','UKINN','UKINP',' USO ','UCDIR',' UCEX','UDisn','UDisp'  &
       eresu(26)=ekt(3); eresu(27)=ekt(1); eresu(28)=ekt(2); eresu(29)=eso;
       eresu(30)=ecodi; eresu(31)=ecoex; eresu(32)=Dispersion(1); eresu(33)=Dispersion(2);
       ! ,'UV2Mn','UV2Mp'
       eresu(34)=v2min(1); eresu(35)=v2min(2);
       !  ,'UECMT','UECMN','UECMP'
       eresu(36)=ECMHFB(3); eresu(37)=ECMHFB(1); eresu(38)=ECMHFB(2);
       !  ,'UROTT','UROTN','UROTP'
       eresu(39)=DEROT(3); eresu(40)=DEROT(1); eresu(41)=DEROT(2);
       !  ,'USQUJT','USQUJN','USQUJP'
       eresu(42)=SQUJ(3); eresu(43)=SQUJ(1); eresu(44)=SQUJ(2);
       !  ,'UCRANT','UCRANN','UCRANP'
       eresu(45)=CRAN(3); eresu(46)=CRAN(1); eresu(47)=CRAN(2);
       !  ,'UERIGT','UERIGN','UERIGP'
       eresu(48)=ERIGHFB(3); eresu(49)=ERIGHFB(1); eresu(50)=ERIGHFB(2);
       !
       ! nucleus with wrong assymptotic
       If(iasswrong(3).Ne.0) eresu(16)=-eresu(16)
       !------------------------------------------------
       ! WRITE UNPROJECTED OUTPUT
       !------------------------------------------------
       Do iw=lout,lfile
          Write(iw,*)
          Write(iw,'(a,9x,a)')            '  NB! From expect (UNPROJECTED RESULTS)'
          Write(iw,*)
          If(iLST1.Ne.0)  &
               Write(iw,'(a,3f15.6)') '  hfb decay const. ass ',ass
          Write(iw,'(a,5f15.6)') '  pairing: CpV0,CpV1,...    ',CpV0,CpV1
          Write(iw,'(a,a)')      '  forces:   ',skyrme
          If(keyblo(1).Ne.0)  &
               Write(iw,'(a,i4,a,f10.3)')  '  Blocked neutron block    ',  &
               bloblo(keyblo(1),1)
          If(keyblo(2).Ne.0)  &
               Write(iw,'(a,i4,a,f10.3)')  '  Blocked proton  block    ',  &
               bloblo(keyblo(2),2)
          Write(iw,*)
          Write(iw,'(/,28x,a,8x,a,9x,a)') ' neutrons ','protons','total'
          Write(iw,'(a,6f15.6)')          '  Requested part.numbs.',tz,Sum(tz)
          Write(iw,'(a,6f15.6)')          '  UnPj(av) part.numbs .',xn
          Write(iw,'(a,3f15.6)')          '  b0, bz, bp ..........',b0,bz,bp
          Write(iw,*)
          Write(iw,'(a,3f15.6)') '  lambda (ala) ........',ala
          Write(iw,'(a,3f15.6)') '  Lambda (alast) ......',alast
          Write(iw,'(a,3f15.6)') '  delta(n,p), pwi .....',del,pwi
          Write(iw,'(a,3f15.6)') '  pairing energy ......',ept
          If(kindhfb.Lt.0) Then
             Write(iw,'(a,3f15.6)') '  LN lambda_2 ... ala2 ',ala2
             Write(iw,'(a,3f15.6)') '  LN energies .........',etr
             Write(iw,'(a,3f15.6)') '  delta(n,p)+ala2 .....',del+ala2
             Write(iw,'(a,3f15.6)') '  Geff(n,p) ...........',Geff
          End If
          Write(iw,*)
          Write(iw,'(a,3f15.6)') '  rms-radius ..........',rms
          Write(iw,'(a,15x,2f15.6)') '  charge-radius, r0 ...',rc,r00
          Write(iw,'(a,3f15.6)') '  deformation beta2....',def
          Write(iw,'(a,3f15.6)') '  dipole moment[fm] ...',(qmoment(1,it),it=1,3)
          Write(iw,'(a,3f15.6)') '  quadrupole moment[b] ',(qmoment(2,it),it=1,3)
          Write(iw,'(a,3f15.6)') '  octupole moment .....',(qmoment(3,it),it=1,3)
          Write(iw,'(a,3f15.6)') '  hexadecapole moment .',(qmoment(4,it),it=1,3)
          Write(iw,'(a,3f15.6)') '  q5 ..................',(qmoment(5,it),it=1,3)
          Write(iw,'(a,3f15.6)') '  q6 ..................',(qmoment(6,it),it=1,3)
          Write(iw,'(a,3f15.6)') '  q7 ..................',(qmoment(7,it),it=1,3)
          Write(iw,'(a,3f15.6)') '  q8 ..................',(qmoment(8,it),it=1,3)
          Write(iw,*)
          Write(iw,'(a,3f15.6)')    '  kinetic energy ......',ekt
          Write(iw,'(a,30x,f15.6)') '  volume energy .......',evol
          Write(iw,'(a,30x,f15.6)') '        rho_tau .......',EVOL_rho_tau
          Write(iw,'(a,30x,f15.6)') '        rho_rho .......',EVOL_rho_rho
          Write(iw,'(a,30x,f15.6)') '  surface energy ......',esurf
          Write(iw,'(a,30x,f15.6)') '   rho_DELTA_rho ......',ESURF_rho_DELTA_rho
          Write(iw,'(a,30x,f15.6)') '   (NABLA_rho)^2 ......',ESURF_NABLA_rho_NABLA_rho
          Write(iw,'(a,30x,f15.6)') '  spin-orbit energy ...',eso
          Write(iw,'(a,30x,f15.6)') '        rho_NABLA_J ...',ESO_rho_NABLA_J
          Write(iw,'(a,30x,f15.6)') '        NABLA_rho_J ...',ESO_NABLA_rho_J
          Write(iw,'(a,30x,f15.6)') '  coulomb energy ......',ecodi+ecoex
          Write(iw,'(a,30x,f15.6)') '          direct ......',ecodi
          Write(iw,'(a,30x,f15.6)') '          exchange ....',ecoex
          Write(iw,'(a,30x,f15.6)') '  tensor energy .......',etens
          Write(iw,'(a,30x,f15.6)') '  direct Hartree E  ...',E_HARTREE_DIR
          Write(iw,'(a,30x,f15.6)') '  Extra E .............',EEXTRA
          Write(iw,'(a,30x,f15.6)') '  External field E ....',E_EXT_FIELD
          Write(iw,'(a,3f15.6)')    '  Entropy .............',entropy
          Write(iw,*)
          Write(iw,'(a,30x,f15.6)')    '  tEnergy: ehfb (qp)...',ehfb
          If(kindhfb.Lt.0) Then
             Write(iw,'(a,30x,f15.6)') '  tEnergy: ehfb(qp)+LN ',etot
          End If
          Write(iw,*)
          Write(iw,'(a,6f15.6)')    '  Calculated but not added corrections '
          Write(iw,'(a,6f15.6)')    '===================================='
          Write(iw,'(a,6f15.6)')    '  cmc-diagonal part ...',ekt/hb0*hbzero-ekt
          Write(iw,'(a,6f15.6)')    '  cmc-hfb .............',ECMHFB
          Write(iw,'(a,6f15.6)')    '  cranking rot corr ...',DEROT
          Write(iw,*)
          Write(iw,'(a,6f15.6)')    '  SQUJ ................',SQUJ
          Write(iw,'(a,6f15.6)')    '  CRAN x 4 ............',4.0_pr*CRAN
          Write(iw,'(a,6f15.6)')    '  Rigit Body ..........',ERIGHFB
          Write(iw,'(a,6f15.6)')
       End Do
       !------------------------------------------------
       ! START corrected Lipkin-Nogami characteristics
       !------------------------------------------------
       If(kindhfb.Lt.0) Then
          Call densitln       !density LN corrections
          If(ierror_flag.Ne.0) Return
          Do it=itmin,itmax
             xn(it)=zero
             rms(it)=zero; q2(it)=zero; q4(it)=zero
          End Do
          !
          Do ihli=1,nghl
             whl=wdcor(ihli)
             rn=ro(ihli,1); rp=ro(ihli,2); rnp2=rn**2+rp**2
             ! rms and deformations
             z=fh(ihli); zz=z*z; rrr=zz+fl(ihli)**2
             p2=p32*zz   -half*rrr
             p3=p53*z*p2 -p23*rrr*z
             p4=p74*z*p3 -p34*rrr*p2
             row=whl*rn
             xn(1)=xn(1)+row
             rms(1)=rms(1)+row*rrr
             q2(1)=q2(1)+row*p2
             q4(1)=q4(1)+row*p4
             row=whl*rp
             xn(2)=xn(2)+row
             rms(2)=rms(2)+row*rrr
             q2(2)=q2(2)+row*p2
             q4(2)=q4(2)+row*p4
          End Do !ihli
          !------------------------------------------------
          ! rms and deformations
          !------------------------------------------------
          Do it=itmin,itmax
             rms(it)=Sqrt(rms(it)/xn(it))
             q2(it)=two*q2(it)     !Qnp=<2r^2P_2(teta)>=<2z^2-x^2-y^2>
             q4(it)=ffdef4*q4(it)  !Hn=<8r^4P_4(teta)>=<8z^4-24z^2(x^2+y^2)+3(x^2+y^2)^2>
             def(it)=Sqrt(pi/5.0_pr)*q2(it)/(rms(it)**2*xn(it))
          End Do
          r212=rms(1)**2; r222=rms(2)**2
          rms(3)=Sqrt((xn(1)*r212+xn(2)*r222)/amas)
          q2(3)=q2(1)+q2(2)  ! quadrupole moment
          q4(3)=q4(1)+q4(2)  ! hexadecapole moment
          def(3)=Sqrt(pi/5.0_pr)*q2(3)/(rms(3)**2*amas) !deformation
          ! other definitions of the same quantitsies
          bet2(1)=ffdef6*q2(1)/(xn(1)*r02) !beta_n=Qn*Sqrt(5Pi)/(3N x^2)
          bet2(2)=ffdef6*q2(2)/(xn(2)*r02) !x=r0=1.2A^(1/3)
          bet2(3)=ffdef6*q2(3)/(amas*r02)
          het4(1)=ffdef7*q4(1)/(xn(1)*r04)
          het4(2)=ffdef7*q4(2)/(xn(2)*r04)
          het4(3)=ffdef7*q4(3)/(amas*r04)
          xn(3)=xn(1)+xn(2)
          bet=def(3)
          rc=Sqrt(r222+0.640_pr)
          ! transitions to barn,barn^2,barn^4
          Do i=1,3
             q2(i)=q2(i)/100.0_pr; q4(i)=q4(i)/10000.0_pr
          End Do
          !------------------------------------------------
          ! STORE to unprojected LN buffer 'eresl'
          !------------------------------------------------
          ! ieresl=20 from module definitions
          ! ,' EHFBLN',' EHFB',' LNEt','LNbet','LNben','LNbep',' LNQt',' LNQn',' LNQp'  &
          eresl(1)=etot; eresl(2)=etot-etr(3);
          eresl(3)=def(3); eresl(4)=def(1); eresl(5)=def(2)
          eresl(6)=q2(3); eresl(7)=q2(1); eresl(8)=q2(2);
          ! ,'LNpEn','LNpEp','LNpDn','LNpDp',' LNrt',' LNrn',' LNrC'  &
          eresl(9)=ept(1); eresl(10)=ept(2); eresl(11)=del(1)+ala2(1); eresl(12)=del(2)+ala2(2);
          eresl(13)=rms(3); eresl(14)=rms(1); eresl(15)=rms(2); eresl(16)=rc;
          ! ,' LNam2n',' LNam2p',' LNe2n',' LNe2p'
          eresl(17)=ala2(1); eresl(18)=ala2(2); eresl(19)=etr(1); eresl(20)=etr(2)
          !------------------------------------------------
          ! WRITE UNPROJECTED LN OUTPUT
          !------------------------------------------------
          Do iw=lout,lfile
             Write(iw,'(a,3f15.6)')
             Write(iw,'(a,3f15.6)') '  With Lipkin-Nogami Corrections'
             Write(iw,'(a,3f15.6)') '================================'
             Write(iw,'(a,3f15.6)') '  rms-radius ..........',rms
             Write(iw,'(a,15x,2f15.6)') '  charge-radius, r0 ...',rc,r00
             Write(iw,'(a,3f15.6)') '  deformation beta ....',def
             Write(iw,'(a,3f15.6)') '  quadrupole moment[b] ',q2
             Write(iw,'(a,3f15.6)') '  hexadecapole moment .',q4
             Write(iw,'(a,3f15.6)') '================================'
             Write(iw,'(a,3f15.6)')
          End Do
       End If
       !------------------------------------------------
       ! WRITE all blocking candidates
       !------------------------------------------------
       If(keyblo(3).Eq.0) Then
          Do iw=lout,lfile
             Write(iw,*)
             Do it=itmin,itmax
                If(it.Eq.1) Then
                   EqpPo=>REqpN; VqpPo=>RVqpN; UqpPo=>RUqpN; KpwiPo=>KpwiN; KqpPo=>KqpN
                Else
                   EqpPo=>REqpP; VqpPo=>RVqpP; UqpPo=>RUqpP; KpwiPo=>KpwiP; KqpPo=>KqpP
                End If
                !
                Write(iw,*) ' ',' Blocking candidates are:'
                Write(iw,*) '  ',protn(it),' eqpmin=',eqpmin(it),' pwiblo=',pwiblo
                Do i=1,blomax(it)
                   Write(iw,'(a,i4,a,i4,a,i4,2x,i4,3(a,1x,f12.8,1x),a)') '    num=',i,  &
                        ' block=',bloblo(i,it),  &
                        ' state=',blo123(i,it),blok1k2(i,it),  &
                        ' Eqp=',EqpPo(KqpPo(blok1k2(i,it))),  &
                        ' (1-2N)E=',(one-two*uk(blok1k2(i,it),it))*EqpPo(KqpPo(blok1k2(i,it))),  &
                        ' Ovlp=',vkmax(blok1k2(i,it),it),  &
                        tb(numax(blok1k2(i,it),it))
                End Do
                Write(iw,*)
             End Do
             Write(iw,*)
          End Do
       End If
       !
       !------------------------------------------------
       ! PAV
       !------------------------------------------------
       ! Projecting on different nucleus
       If(iproj.Ne.0) Then
          npr(1)=Int(npr1pj); npr(2)=Int(npr2pj)
          !tz(1)=Real(npr(1),Kind=pr); tz(2)=Real(npr(2),Kind=pr)
          Call expectpj(.True.)
       End If
    End If
    !
    If (IDEBUG.Eq.1) Call get_CPU_time('expect',1)
    !
  End Subroutine expect
  !=======================================================================
  !
  !===============================================================================================
  Subroutine Constraint_or_not(inin_INI0,inin0,icstr0)
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr), Intent(in) :: inin_INI0
    Integer(ipr), Intent(inout) :: inin0,icstr0
    If(SUM(lambda_active).Gt.0) Then
       icstr0=1; inin0=inin_INI0
    Else
       icstr0=0; inin0=inin_INI0
    End If
  End Subroutine Constraint_or_not
  !===============================================================================================
  !
  !===============================================================================================
  Subroutine moments_setUnits()
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: lambda
    Real(pr) :: sqr4pi
    !------------------------------------------------------
    ! Defines standard units for multipole moments
    !------------------------------------------------------
    If (.Not.Allocated(q_units)) Allocate(q_units(0:lambdaMax)); q_units = one

    sqr4pi=Sqrt(pp16*Atan(one))

    q_units(0)=+sqr4pi
    q_units(1)=+sqr4pi/Sqrt(three)
    q_units(2)=+sqr4pi/Sqrt(five)*two

    Do lambda=0,lambdaMax
       q_units(lambda)=q_units(lambda) / ten**lambda
    End Do

    Return
  End Subroutine moments_setUnits
  !===============================================================================================
  !
  !===============================================================================================
  Subroutine moments_computeValue()
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: lambda,ihli
    Real(pr), Dimension(0:8) :: Qval
    Real(pr) :: sqr4pi,z,z2,z3,z4,z5,z6,z7,z8,rrr,rrr4,rrr6
    Real(pr) :: rown,rowp,whl,rn,rp
    !------------------------------------------------------
    ! Expectation value of multipole moments
    !------------------------------------------------------
    sqr4pi=one/Sqrt(pp16*Atan(one))
    !
    qmoment=zero; Qval=zero
    !
    Do ihli=1,nghl
       !
       whl=wdcor(ihli)
       rn=ro(ihli,1); rp=ro(ihli,2)
       rown=whl*rn; rowp=whl*rp;
       z=fh(ihli); rrr=fl(ihli)**2
       !
       Call moments_valueMesh(z,rrr,Qval)
       !
       Do lambda=0,lambdaMax
          qmoment(lambda,1)=qmoment(lambda,1)+rown*Qval(lambda)
          qmoment(lambda,2)=qmoment(lambda,2)+rowp*Qval(lambda)
       End Do
       !
    End Do
    !
    Do lambda=0,lambdaMax
       qmoment(lambda,3)=qmoment(lambda,1)+qmoment(lambda,2)
    End Do

    Return
  End Subroutine moments_computeValue
  !===============================================================================================
  !
  !===============================================================================================
  !
  Subroutine moments_valueMesh(z,rrr,Qval)
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: lambda
    Real(pr), Dimension(0:8) :: Qval
    Real(pr) :: sqr4pi,z,z2,z3,z4,z5,z6,z7,z8,rrr,rrr4,rrr6
    !------------------------------------------------------
    ! Expectation value of multipole moments
    !------------------------------------------------------
    !
    sqr4pi=one/Sqrt(pp16*Atan(one))
    !
    z2=z*z; z3=z2*z; z4=z3*z; z5=z4*z; z6=z5*z; z7=z6*z; z8=z7*z
    rrr4=rrr*rrr; rrr6=rrr4*rrr
    !
    Qval(0) =               sqr4pi
    Qval(1) = Sqrt(three)  *sqr4pi          * z
    Qval(2) = Sqrt(five)   *sqr4pi*half     * (two*z2-        rrr)
    Qval(3) = Sqrt(seven)  *sqr4pi*half     * (two*z3-three*z*rrr)
    Qval(4) = Sqrt(nine)   *sqr4pi*p18      * (eight*z4-24.0_pr*z2*rrr    + three *rrr4)
    Qval(5) = Sqrt(11.0_pr)*sqr4pi*p18      * (eight*z5-   pp40*z3*rrr    +pp15*z *rrr4)
    Qval(6) = Sqrt(13.0_pr)*sqr4pi/pp16     * (pp16*z6-120.0_pr*z4*rrr+ 90.0_pr*z2*rrr4-five     *rrr6)
    Qval(7) = Sqrt(15.0_pr)*sqr4pi/pp16     * (pp16*z7-168.0_pr*z5*rrr+210.0_pr*z3*rrr4-35.0_pr*z*rrr6)
    Qval(8) = Sqrt(17.0_pr)*sqr4pi/128.0_pr * (128.0_pr*z8-1792.0_pr*z6*rrr +3360.0_pr*z4*rrr4 &
                                                          -1120.0_pr*z2*rrr6+  35.0_pr*rrr4*rrr4)
    !
    If(Parity) Then
       Qval(1)=zero; Qval(3)=zero;Qval(5)=zero; Qval(7)=zero
    End If
    !
    Do lambda=0,lambdaMax
       Qval(lambda)=Qval(lambda)*q_units(lambda)
    End Do
    !
    Return
  End Subroutine moments_valueMesh
  !=======================================================================
  !
  !=======================================================================
  Subroutine moments_computeField(lambda,ib)
    !---------------------------------------------------------------------
    ! calculates fields in r-space form axially symmetric densities
    !---------------------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: lambda
    Integer(ipr) :: i,ih,il,ib,nd,nd2,ihli,ihil,im,n1,n2
    Integer(ipr) :: ja,jb,nsa,nsb,nsab,ssu,ssd
    Real(pr)    :: qhla,vh,fiun1,fiun2,fidn1,fidn2,vnhl
    Real(pr), Allocatable :: Vmom(:)
    Real(pr), Dimension(0:8) :: Qval
    Real(pr) :: z,rrr
    Integer(ipr) :: ndxmax
    Parameter(ndxmax=(n00max+2)*(n00max+2)/4)
    Real(pr) :: OMPFIU(ndxmax),OMPFID(ndxmax)
    !
    Allocate(Vmom(1:nghl))
    !
    Qval=zero
    !
    ! Compute moment lambda on integration mesh
    Do ihli = 1,nghl
       z=fh(ihli);rrr=fl(ihli)**2
       Call moments_valueMesh(z,rrr,Qval)
       Vmom(ihli)=Qval(lambda)
    End Do !ihli
    !
    ! Form matrix of the multipole constraint lambda in HO basis
    nd=id(ib); nd2=nd*nd; im=ia(ib)
    ! sum over gauss integration points
    Do ihil=1,nghl
       vnhl=Vmom(ihil)
       ! scan over basis states
       Do n1=1,nd
          ja=n1+im; nsa=NS(ja); ssu=Max(nsa,0); ssd=Max(-nsa,0)
          QHLA=QHLA_opt(ja,ihil)
          OMPFIU(N1)=QHLA*ssu
          OMPFID(N1)=QHLA*ssd
       End Do
       i=0
       Do n1=1,nd
          ja=n1+im; nsa=NS(ja)
          fiun1=OMPFIU(N1); fidn1=OMPFID(N1)
          do n2=1,n1
             i=i+1; nsb=NS(n2+im); nsab=nsa+nsb; vh=0.0_pr
             If (nsab.Ne.0) Then
                 If (nsb.Gt.0) Then
                     fiun2 = OMPFIU(N2)
                     vh    = fiun1*fiun2
                 Else
                     fidn2 = OMPFID(N2)
                     vh    = fidn1*fidn2
                 End If
                 multMatElems(i)=multMatElems(i)+vh*vnhl
             End If
          End Do !n2
       End Do !n1
    End Do !ihil
    Deallocate(Vmom)
    !
    Return
  End Subroutine moments_computeField
  !===============================================================================================
  !
  !===============================================================================================
  Subroutine getLagrange(ite)
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Character(Len=1) :: trans
    Integer(ipr) :: ite,icons,lambda,icount,it,i,j,l,ierror
    Integer(ipr) :: ib,nd,nd2,nhfb,i0,m,k1,k2,n1,n2,nd1,k,kk,ll
    Integer(ipr) :: i_uvN,i_uvP,incx,incy
    Integer(ipr), allocatable :: ipivot(:),iftN(:),iftP(:)
    Real(pr) :: minu,hla,t_term,temp_k,temp_l,result,brakev,epsilo
    Real(pr), allocatable :: EqpN(:),EqpP(:)
    Real(pr), allocatable :: vecold(:),qmultt(:),veclam(:),veccns(:)
    Real(pr), allocatable :: cnsorg(:,:),cnsmat(:,:),cnsvec(:)
    Real(pr), allocatable :: fn12pl(:,:,:),fp12pl(:,:,:)
    Real(pr), allocatable :: fn11pl(:,:,:),fp11pl(:,:,:),fn11mi(:,:,:),fp11mi(:,:,:)
    Real(pr), allocatable :: doubln(:,:),doublp(:,:),dsum_n(:,:),dsum_p(:,:)
    Real(pr), allocatable :: workcn(:),dblmul(:,:),Umatr(:,:),Vmatr(:,:)
    !---------------------------------------------------------------------
    ! This routine updates the  lagrange  multipliers of the multi-
    ! dimensional linear constraints  based on the variation of the
    ! generalized density matrix and the rpa matrix at the cranking
    ! approximation.
    !
    ! references:  (1) phys. rev. c21, 1568 (1980)
    !              (2) phys. rev. c80, 054313 (2009)
    !---------------------------------------------------------------------
    !
    minu=-one
    epsilo=1.E-14
    !
    ! initializing the multipole moment template array
    Allocate(qmultt(0:lambdaMax));qmultt=zero
    Do lambda=0,lambdaMax
       qmultt(lambda)=qmoment(lambda,3)
    End Do
    !
    ! constructing the vector of the deviations of the current constraint
    ! from the requested values
    !
    Allocate(vecold(1:numberCons));vecold=zero
    Allocate(cnsvec(1:numberCons));cnsvec=zero
    Allocate(veclam(1:numberCons));veclam=zero
    !
    Do icons=1,numberCons
       lambda=multLambda(icons)
       cnsvec(icons)=multRequested(lambda)-qmultt(lambda)
       If (nbroyden.lt.1) Then
           vecold(icons)=multLag(lambda)
       Else
           vecold(icons)=brin(nhhdim4+lambda)
       End If
       veclam(icons)=vecold(icons)
    End Do
    !
    ! proceeding to determine the matrix of the constraint operators
    ! in the q.p. basis
    ! loop over the K blocks
    Allocate(cnsmat(numberCons,numberCons));cnsmat=zero
    Allocate(cnsorg(numberCons,numberCons));cnsorg=zero
    !
    i_uvN=0 ! new index referring to all q.p. vectors
    i_uvP=0 ! new index referring to all q.p. vectors
    !
    Do ib = 1,nb
       !
       !------------------------------------------------------
       ! matrix of the constraint in q.p. basis
       !------------------------------------------------------
       !
       !------------------------------------------------------
       ! neutron sector
       !------------------------------------------------------
       !
       it=1
       !
       nd=id(ib); nd2=nd*nd; nhfb=nd+nd; i0=ia(ib); m=ib+(it-1)*nbx
       !
       If(kd(ib,it).Gt.0) Then
          Allocate(doubln(nd,kd(ib,it))); doubln=zero
          Allocate(fn12pl(kd(ib,it),kd(ib,it),numberCons)); fn12pl=zero
          Allocate(Umatr(nd,kd(ib,it))); Umatr=zero
          Allocate(Vmatr(nd,kd(ib,it))); Vmatr=zero
          Allocate(EqpN(kd(ib,it))); EqpN=zero
          Allocate(ifTN(kd(ib,it))); ifTN=1
          ! temperature
          If(switch_on_temperature) Then
             Allocate(fn11pl(kd(ib,it),kd(ib,it),numberCons)); fn11pl=zero
          End If
          !
          ! U and V for this block (v. 101)
          Do k=1,kd(ib,it)
             ifTN(k)=ka(ib,it)+k; kk=KqpN(ka(ib,it)+k); EqpN(k)=REqpN(kk)
             Do n1=1,nd
                i_uvN=i_uvN+1
                Vmatr(n1,k)=RVqpN(i_uvN)
                Umatr(n1,k)=RUqpN(i_uvN)
             End Do
          End Do
          !
          Do icons=1,numberCons
             !
             Allocate(multMatElems(1:nd2)); multMatElems=zero
             !
             lambda=multLambda(icons); Call moments_computeField(lambda,ib)
             !
             ! matrix of the constraints in HO basis (size nd x nd)
             Allocate(dblmul(nd,nd));dblmul=zero
             j=0
             Do n1=1,nd
                Do n2=1,n1
                   j=j+1;hla=multMatElems(j)
                   dblmul(n1,n2)=hla;dblmul(n2,n1)=hla
                End Do
             End Do
             !
             ! matrix of the constraint operator in the qp basis. due to
             ! the q.p. cut-off the actual size of the q.p. basis is not
             ! the same as the s.p. (ho) basis, and it is  not the  same
             ! for protons and neutrons.  the formulas implemented below
             ! differ from the 2 references for 3 reasons:
             !  - different  phase convention for the bogoliubov  matrix
             !  - block structure of the bogoliubov matrix in hfodd
             !  - storage in a() and b() arrays correspond to complex
             !    conjugate of the actual matrices
             !
             ! second term: v^{+} f^{*} u^{*} = v^{T} f u
             Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Umatr,nd,zero,doubln,nd)
             Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Vmatr,nd,doubln,nd,zero,fn12pl(1,1,icons),kd(ib,it))
             !
             ! first term:  u^{+} f v^{*} = u^{T} f v
             Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Vmatr,nd,zero,doubln,nd)
             Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Umatr,nd,doubln,nd,minu,fn12pl(1,1,icons),kd(ib,it))
             !
             ! temperature - computing \tilde{f}^{11}
             If(switch_on_temperature) Then
               !
               ! second term: v^{+} f^{*} v = v^{T} f v
               Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Vmatr,nd,zero,doubln,nd)
               Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Vmatr,nd,doubln,nd,zero,fn11pl(1,1,icons),kd(ib,it))
               ! first term:  u^{+} f u = u^{T} f u
               Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Umatr,nd,zero,doubln,nd)
               Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Umatr,nd,doubln,nd,minu,fn11pl(1,1,icons),kd(ib,it))
               !
             End If
             !
             Deallocate(multMatElems)
             Deallocate(dblmul)
             !
          End Do ! end icons (neutrons)
          !
          Deallocate(doubln,Umatr,Vmatr)
       End If
       !
       !------------------------------------------------------
       ! Proton sector
       !------------------------------------------------------
       !
       it=2
       !
       nd=id(ib); nd2=nd*nd; nhfb=nd+nd; i0=ia(ib); m=ib+(it-1)*nbx
       !
       If(kd(ib,it).Gt.0) Then
          Allocate(doublp(nd,kd(ib,it))); doublp=zero
          Allocate(fp12pl(kd(ib,it),kd(ib,it),numberCons)); fp12pl=zero
          Allocate(Umatr(nd,kd(ib,it))); Umatr=zero
          Allocate(Vmatr(nd,kd(ib,it))); Vmatr=zero
          Allocate(EqpP(kd(ib,it))); EqpP=zero
          Allocate(ifTP(kd(ib,it))); ifTP=1
          ! temperature
          If(switch_on_temperature) Then
             Allocate(fp11pl(kd(ib,it),kd(ib,it),numberCons)); fp11pl=zero
          End If
          !
          ! U and V for this block
          Do k=1,kd(ib,it)
             ifTP(k)=ka(ib,it)+k; kk=KqpP(ka(ib,it)+k); EqpP(k)=REqpP(kk)
             Do n1=1,nd
                i_uvP=i_uvP+1
                Vmatr(n1,k)=RVqpP(i_uvP)
                Umatr(n1,k)=RUqpP(i_uvP)
             End Do
          End Do
          !
          Do icons=1,numberCons
             !
             Allocate(multMatElems(1:nd2)); multMatElems=zero
             !
             lambda=multLambda(icons); Call moments_computeField(lambda,ib)
             !
             ! matrix of the constraints in HO basis (size nd x nd)
             Allocate(dblmul(nd,nd));dblmul=zero
             j=0
             Do n1=1,nd
                Do n2=1,n1
                   j=j+1;hla=multMatElems(j)
                   dblmul(n1,n2)=hla;dblmul(n2,n1)=hla
                End Do
             End Do
             !
             ! matrix of the constraint operator in the qp basis. due to
             ! the q.p. cut-off the actual size of the q.p. basis is not
             ! the same as the s.p. (ho) basis, and it is  not the  same
             ! for protons and neutrons.  the formulas implemented below
             ! differ from the 2 references for 3 reasons:
             !  - different  phase convention for the bogoliubov  matrix
             !  - block structure of the bogoliubov matrix in hfodd
             !  - storage in a() and b() arrays correspond to complex
             !    conjugate of the actual matrices
             !
             ! second term: v^{+} f^{*} u^{*} = v^{t} f u
             Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Umatr,nd,zero,doublp,nd)
             Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Vmatr,nd,doublp,nd,zero,fp12pl(1,1,icons),kd(ib,it))
             !
             ! first term:  u^{+} f v^{*} = u^{t} f v
             Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Vmatr,nd,zero,doublp,nd)
             Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Umatr,nd,doublp,nd,minu,fp12pl(1,1,icons),kd(ib,it))
             !
             ! temperature - computing \tilde{f}^{11}
             If(switch_on_temperature) Then
               !
               ! second term: v f^{*} v = v^{T} f v
               Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Vmatr,nd,zero,doublp,nd)
               Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Vmatr,nd,doublp,nd,zero,fp11pl(1,1,icons),kd(ib,it))
               ! first term:  u f u = u^{T} f u
               Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Umatr,nd,zero,doublp,nd)
               Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Umatr,nd,doublp,nd,minu,fp11pl(1,1,icons),kd(ib,it))
               !
             End If
             !
             Deallocate(dblmul)
             Deallocate(multMatElems)
             !
          End Do ! end icons (protons)
          !
          Deallocate(doublp,Umatr,Vmatr)
       End If
       !
       !------------------------------------------------------
       ! constraint correlation matrix
       !------------------------------------------------------
       !
       !------------------------------------------------------
       ! neutron sector
       !------------------------------------------------------
       !
       it=1
       !
       If(kd(ib,it).Gt.0) Then
          Allocate(doubln(kd(ib,it),kd(ib,it))); doubln=zero
          Allocate(dsum_n(kd(ib,it),kd(ib,it))); dsum_n=zero
          !
          Do i=1,numberCons
             Do j=1,numberCons
                !
                ! temperature
                If((.Not.switch_on_temperature)) Then
                   !
                   Do l=1,kd(ib,it)
                      Do k=1,kd(ib,it)
                         If(Abs(EqpN(k)+EqpN(l)).Gt.Epsilo) Then
                            doubln(k,l)=fn12pl(k,l,i)/(EqpN(k)+EqpN(l))
                         Else
                            doubln(k,l)=zero
                         End If
                      End Do
                   End do
                   !
                   Call dgemm('t','n',kd(ib,it),kd(ib,it),kd(ib,it),one,fn12pl(1,1,j),kd(ib,it),&
                                                         doubln,kd(ib,it),zero,dsum_n,kd(ib,it))
                Else
                   !
                   ! term corresponding to f^12
                   Do l=1,kd(ib,it)
                      Do k=1,kd(ib,it)
                         kk=iftN(k);ll=iftN(l)
                         temp_k=fn_T(kk)
                         temp_l=fn_T(ll)
                         If(Abs(EqpN(k)+EqpN(l)).Gt.Epsilo) Then
                            doubln(k,l)=fn12pl(k,l,i)*(one+temp_k+temp_l)/(EqpN(k)+EqpN(l))
                         Else
                            doubln(k,l)=zero
                         End If
                      End Do
                   End Do
                   !
                   Call dgemm('t','n',kd(ib,it),kd(ib,it),kd(ib,it),one,fn12pl(1,1,j),kd(ib,it),&
                                                         doubln,kd(ib,it),zero,dsum_n,kd(ib,it))
                   !
                   ! first term: positive simplex
                   Do l=1,kd(ib,it)
                      Do k=1,kd(ib,it)
                         kk=iftN(k);ll=iftN(l)
                         temp_k=fn_T(kk)
                         temp_l=fn_T(ll)
                         If(k.ne.l.And.(Abs(EqpN(k)-EqpN(l)).Gt.Epsilo)) Then
                            t_term=-(temp_k-temp_l)/(EqpN(k)-EqpN(l))
                         Else
                            t_term=-temp_k*(temp_k-one)/temper
                         End If
                         doubln(k,l)=half*t_term*fn11pl(k,l,i)
                      End Do
                   End Do
                   !
                   Call dgemm('t','n',kd(ib,it),kd(ib,it),kd(ib,it),one,fn11pl(1,1,j),kd(ib,it),&
                                                          doubln,kd(ib,it),one,dsum_n,kd(ib,it))
                   !
                End If
                !
                ! taking the trace of the resulting matrix
                !
                result=zero
                Do l=1,kd(ib,it)
                   result=result+dsum_n(l,l)
                End Do
                !
                cnsmat(i,j)=cnsmat(i,j)+0.5*result
                !
             End Do ! end of loop over j constraint
          End Do ! end of loop over i constraint
          !
          Deallocate(doubln,dsum_n,fn12pl,EqpN,ifTN)
          If(switch_on_temperature) Deallocate(fn11pl)
       End If
       !
       !------------------------------------------------------
       ! proton sector
       !------------------------------------------------------
       !
       it=2
       !
       If(kd(ib,it).Gt.0) Then
          Allocate(doublp(kd(ib,it),kd(ib,it))); doublp=zero
          Allocate(dsum_p(kd(ib,it),kd(ib,it))); dsum_p=zero
          !
          Do i=1,numberCons
             Do j=1,numberCons
                !
                ! temperature
                If((.Not.switch_on_temperature)) Then
                   !
                   Do l=1,kd(ib,it)
                      Do k=1,kd(ib,it)
                         If(Abs(EqpP(k)+EqpP(l)).Gt.Epsilo) Then
                            doublp(k,l)=fp12pl(k,l,i)/(EqpP(k)+EqpP(l))
                         Else
                            doublp(k,l)=zero
                         End If
                      End Do
                   End do
                   !
                   Call dgemm('t','n',kd(ib,it),kd(ib,it),kd(ib,it),one,fp12pl(1,1,j),kd(ib,it),&
                                                         doublp,kd(ib,it),zero,dsum_p,kd(ib,it))
                Else
                   !
                   ! term corresponding to f^12
                   Do l=1,kd(ib,it)
                      Do k=1,kd(ib,it)
                         kk=iftP(k);ll=iftP(l)
                         temp_k=fp_T(kk)
                         temp_l=fp_T(ll)
                         If(Abs(EqpP(k)+EqpP(l)).Gt.Epsilo) Then
                            doublp(k,l)=fp12pl(k,l,i)*(one+temp_k+temp_l)/(EqpP(k)+EqpP(l))
                         Else
                            doublp(k,l)=zero
                         End If
                      End Do
                   End Do
                   !
                   Call dgemm('t','n',kd(ib,it),kd(ib,it),kd(ib,it),one,fp12pl(1,1,j),kd(ib,it),&
                                                         doublp,kd(ib,it),zero,dsum_p,kd(ib,it))
                   !
                   ! first term: positive simplex
                   Do l=1,kd(ib,it)
                      Do k=1,kd(ib,it)
                         kk=iftP(k);ll=iftP(l)
                         temp_k=fp_T(kk)
                         temp_l=fp_T(ll)
                         If(k.ne.l.And.(Abs(EqpP(k)-EqpP(l)).Gt.Epsilo)) Then
                            t_term=-(temp_k-temp_l)/(EqpP(k)-EqpP(l))
                         Else
                            t_term=-temp_k*(temp_k-one)/temper
                         End If
                         doublp(k,l)=half*t_term*fp11pl(k,l,i)
                      End Do
                   End Do
                   !
                   Call dgemm('t','n',kd(ib,it),kd(ib,it),kd(ib,it),one,fp11pl(1,1,j),kd(ib,it),&
                                                          doublp,kd(ib,it),one,dsum_p,kd(ib,it))
                   !
                End If
                !
                ! taking the trace of the resulting matrix
                !
                result=zero
                Do l=1,kd(ib,it)
                   result=result+dsum_p(l,l)
                End Do
                !
                cnsmat(i,j)=cnsmat(i,j)+0.5*result
                !
             End Do ! end of loop over j constraint
          End Do ! end of loop over i constraint
          !
          Deallocate(doublp,dsum_p,fp12pl,EqpP,ifTP)
          If(switch_on_temperature) Deallocate(fp11pl)
       End If
       !
    End Do ! end of loop over blocks ib
    !
    ! computing the inverse of the correlation matrix
    cnsorg=cnsmat
    !
    ierror=0
    Allocate(ipivot(numberCons))
    Call dgetrf(numberCons,numberCons,cnsmat,numberCons,ipivot,ierror)
    !
    ierror=0
    Allocate(workcn(numberCons))
    Call dgetri(numberCons,cnsmat,numberCons,ipivot,workcn,numberCons,ierror)
    Deallocate(ipivot)
    !
    ! constructing the vector of variations of the linear constraints
    trans='N'; incx=1; incy=1
    Call dgemv(trans,numberCons,numberCons,one,cnsmat,numberCons,cnsvec,incx,zero,workcn,incy)
    !
    ! updating the linear constraint vector (mixing has to be done simultaneously).
    If (ite.Eq.0) Then
        brakev=zero
    Else
        brakev=xmix
    End If
    !
    Allocate(veccns(numberCons))
    Do i=1,numberCons
       veccns(i)=veclam(i)+workcn(i)
       lambda=multLambda(i)
       If(nbroyden.lt.1) Then
          multLag(lambda)=brakev*vecold(i)+(1.0-brakev)*veccns(i)
       Else
          multLag(lambda)=veccns(i)
          brout(nhhdim4+lambda)=multLag(lambda)
       End If
    End Do
    !
    Deallocate(veccns,vecold,workcn)
    Deallocate(cnsmat,cnsorg)
    Deallocate(qmultt,cnsvec,veclam)
    !
    Return
  End Subroutine getLagrange
  !===============================================================================================
  !
  !===============================================================================================
  Subroutine requested_blocked_level(ib,it)
    !------------------------------------------------------
    ! Search for the requested state to block
    !------------------------------------------------------
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr), Intent(in)  :: ib,it
    Integer(ipr) :: nd,im,k,ndk,na2,nad2,iqn,k0,LAPLUS,OMEGA,n1,n2,n3
    Real(pr)     :: s1,s2,UUk,VVk
    k0=0
    If(nkblo(it,2).Eq.0)      Return
    If(Parity) Then
       LAPLUS=(ib+1)/2 !Yesp
    Else
       LAPLUS=ib       !Nop
    End If
    OMEGA=2*LAPLUS-1
    If(nkblo(it,1).Ne.OMEGA) Return
    nd=ID(ib); im=ia(ib);
    Do k=1,nd
       ndk=k+nd; s1=zero
       Do na2=1,nd
          nad2=na2+nd
          UUk=allhfb(ib)%arr(na2,ndk)
          VVk=allhfb(ib)%arr(nad2,ndk)
          s2=Max(s1,Abs(UUk),Abs(VVk))
          If(s2.Gt.s1) Then
             s1=s2
             iqn=na2+im  ! the position in [123] numbering
          End If
       End Do
       ! quantum numbers: Omega,P[n1,n2,n3]=>OMEGA,tpar(npar(iqn))[nz(iqn)+2*nr(iqn)+nl(iqn),nz(iqn),nl(iqn)]
       If(nkblo(it,2).Ne.tpar(npar(iqn))) Cycle
       n3=nl(iqn);          If(nkblo(it,5).Ne.n3) Cycle
       n2=nz(iqn);          If(nkblo(it,4).Ne.n2) Cycle
       n1=n2+2*nr(iqn)+n3;  If(nkblo(it,3).Ne.n1) Cycle
       k0=iqn
       keyblo(it)=1
       bloblo(keyblo(it),it)=ib
       blo123(keyblo(it),it)=k
       Exit
    End Do
  End Subroutine requested_blocked_level
  !===================================================================================================================================
  !#END HFBTHO_SOLVER
  !===================================================================================================================================

End Module hfodd_hfbtho
!===================================================================================================================================
!#END hfodd_hfbtho
!===================================================================================================================================

