#----------------------------------------------------------------------!
#                                                                      !
#                Makefile for HFODD                                    !
#                                                                      !
# Authors: N. Schunck, LLNL                                            !
#                                                                      !
# This compiles HFODD (version > 249o)                                 !
#                                                                      !
#                          First included in official release v249a    !
#                                                                      !
#----------------------------------------------------------------------!

#======================================================================#
#                                                                      #
#              HFODD Suite - Makefile (version >= 249o)                #
#                                                                      #
#  NOTA BENE: When running this Makefile, make sure:                   #
#             - all files are in UNIX format (this includes all        #
#               sources .f and .f90, but also this Makefile)           #
#             - all indentations after the 'Nothing beyond this line   #
#               should be changed, in principle' header are proper     #
#               tabulations, not spaces...                             #
#                                                                      #
#  Compiler options are preset for the following compilers: Portland   #
#  Intel, GNU, Pathscale, Lahey, Cray and IBM (on Blue Gene Q).        #
#                                                                      #
#  Makefile Options                                                    #
#  ================                                                    #
#                                                                      #
#   - COMPILER .....: defines the type of compiler. Choices are:       #
#                     PGI, IFORT, GFORTRAN, PATHSCALE, LAHEY, CRAY,    #
#                     IBM                                              #
#                                                                      #
#   - FORTRAN_MPI ..: gives the executable name of the compiler        #
#                                                                      #
#   - DEBUG ........: turns on all optimizations on (FALSE) or off     #
#                     (TRUE). In the latter case, debugging options    #
#                     are set.                                         #
#                                                                      #
#   - PORTABLE .....: chooses to build your own version of BLAS        #
#                     and LAPACK libraries (=TRUE), based on the       #
#                     sources which must be in the same directory.     #
#                     If FALSE, the compiler will try to link to       #
#                     system libraries contained in the variable       #
#                     LINEAR_ALGEBRA.                                  #
#                                                                      #
#   - LAPACK_LIBRARY: Linking to optimized BLAS and LAPACK             #
#                     libraries. The user is responsible to            #
#                     set up this library on his/her system.           #
#                     Note that ACML (like the competitor MKL          #
#                     from Intel) has threaded BLAS: together          #
#                     with USE_OPENMP=1, this can significantly        #
#                     accelerate the execution.                        #
#                                                                      #
#   - STATIC_LIB ...: TRUE or FALSE, allows static linking             #
#                                                                      #
#   - LARGE_MEMORY .: if TRUE, compiler options are set to enable      #
#                     arrays taking more than 2GB in memory, if        #
#                     FALSE, use of default compiler options           #
#                                                                      #
#  Pre-processor Options: compatibility options                        #
#  =====================                                               #
#                                                                      #
#   - SWITCH_PORT ..: if 1, bold memory aliases in DENSHF, COPDEN,     #
#                     SYMDEN, COPPAI and SYMPAI are de-activated       #
#                     and replaced by safer array assignments. Also    #
#                     a few additional variables are properly ini-     #
#                     tialized. These initializations may result       #
#                     in a noticeable slowdown of the code.            #
#                                                                      #
#   - SWITCH_ESSL ..: if 1, the code calls routines from the IBM ESSL  #
#                     library, if 0 it calls standard LAPACK routines. #
#                     Note that BLAS routines are used throughout the  #
#                     code.                                            #
#                                                                      #
#   - SWITCH_DIAG ..: chooses the method of diagonalization:           #
#                       - if 1, uses ZHPEV                             #
#                       - if 2, uses ZHPEVX                            #
#                       - if 3, uses ZHEEVR (usually the fastest)      #
#                       - if 3, uses ZHEEVD                            #
#                                                                      #
#   - SWITCH_QUAD ..: if 1, promotes a few REAL and COMPLEX  to        #
#                     SWITCH_QUAD precision. Not all compilers         #
#                     accept this, GNU, Lahey and Intel Fortran do.    #
#                                                                      #
#   - SWITCH_VECT ..: if 1, loops are vectorized in INTCEN and         #
#                     INTCOU, otherwise, loops are optimized for       #
#                     scalar architecture.                             #
#                                                                      #
#   - SWITCH_COULEX : if 1, the code will compute the exact Coulomb    #
#                     exchange interaction energy between fragments    #
#                     (using the routine COUENE and the expansion of   #
#                     1/r as a sum of Gaussians). Costly.              #
#                                                                      #
#  Pre-processor Options: conditional compilation                      #
#  =====================                                               #
#                                                                      #
#   - USE_OPENMP ...: chooses to switch on multi-threading (if 1)      #
#                     or not (of 0). For this option to actually       #
#                     have an effect, the user must have set the       #
#                     OMP_NUM_THREADS Unix environment variable to     #
#                     some >1 value (typically 4).                     #
#                                                                      #
#   - USE_MPI ......: chooses to run the code in parallel mode         #
#                     (if 1) or in single-core version (if 0).         #
#                     For this option to actually have an effect,      #
#                     the user must have the Message Passing           #
#                     Interface (MPI) installed and running on         #
#                     his/her system. Programs are usually run         #
#                     using the wrappers 'mpirun' or 'mpiexec'.        #
#                                                                      #
#   - USE_MANYCORES : if 1, the code will use several MPI tasks to     #
#                     solve the HFB equations. Currently, only the     #
#                     routines DENSHF and GAUOPT have been parallized. #
#                     Note that filenames will take the generic form   #
#                     proc_[XXX]_hfodd_[XXXXXX].out (and similarly     #
#                     for the other files), with XXX the current MPI   #
#                     task and XXXXXX the HFB task.                    #
#                                                                      #
#   - USE_SCALAPACK : Parallel diagonalization of the HF or HFB matrix #
#                     for simplex-breaking case only. Can only be      #
#                     active if USE_MPI=1 and the user must have the   #
#                     ScaLAPACK library.                               #
#                                                                      #
#   - USE_FITS      : HFODD is used as a subroutine and fit            #
#                     input and output variables are handled by        #
#                     subroutines FITINP and FITOUT, respectively.     #
#                                                                      #
#======================================================================#

SHELL := /bin/bash

# Version numbers
VERSION_MAIN       = 273y
VERSION_SIZES      = 2
VERSION_MODULES    = 17b
VERSION_HFBTHO     = 200j
VERSION_INTERFACE  = 4
VERSION_FUNCTIONAL = 3
VERSION_MPIIO      = 5b
VERSION_MPIMANAGER = 4b
VERSION_SHELL      = 4b
VERSION_SL         = 3
VERSION_FISSION    = 7b
VERSION_PAIRS      = 1
VERSION_PNP        = 6
VERSION_FITS       = 15
VERSION_LIPCORR    = 7d

# Defines version number of BLAS/LAPACK library on LLNL or ANL computing facilities
#  MKL (LLNL) ...........: mkl-11.3.2
#  ESSL (Vulcan, LLNL) ..: 5.1
VERSION_LIBRARY = mkl-11.3.2

# Makefile Options
#  - Compiler types and names (note: on LC, ifort version 17 seems bugged)
#     LC .............: mpiifort-16.0.258, ifort-16.0.258
#     BG/Q ...........: mpixlf90_r, xlf90_r
#     NERSC ..........: ftn
#     OLCF ...........: ftn
#     local ..........: gfortran
COMPILER       = GFORTRAN
FORTRAN_MPI    = gfortran
DEBUG          = FALSE
PORTABLE       = FALSE
LAPACK_LIBRARY = TRUE
STATIC_LIB     = FALSE
LARGE_MEMORY   = FALSE

# Preprocessor options
SHOW_SOURCE   = FALSE
#  - compatibility options
#       * SWITCH_ESSL=1 may require SWITCH_DIAG=2
#       * USE_SCALAPACK=1 may require SWITCH_DIAG=4
SWITCH_PORT   = 1
SWITCH_ESSL   = 0
SWITCH_DIAG   = 3
SWITCH_QUAD   = 0
SWITCH_VECT   = 1
SWITCH_COULEX = 0
#  - conditional compilation: multi-treading, MPI, multi-core (several MPI
#    tasks per HFB calculation), Scalapack for diagonalization (only for
#    simplex-breaking case).
USE_OPENMP    = 0
USE_MPI       = 0
USE_MANYCORES = 0
USE_SCALAPACK = 0
USE_FITS      = 0
# If multi-core or Scalapack is required, we must define at compilation
# the characteristics of the processor grid. These options will be
# passed to the program using the C preprocessor
ifeq ($(USE_MANYCORES),1)
      M_GRID = 2
      N_GRID = 2
else
      M_GRID =
      N_GRID =
endif

# Names and paths to BLAS and LAPACK Libraries
#  * BG/Q ARCHITECTURES (ESSL)
#      ALCF (OpenMP) ..: -L/soft/libraries/alcf/current/xl/LAPACK/lib -llapack \
#                        -L/soft/libraries/alcf/current/xl/BLAS/lib -lblas
#      ALCF (ESSL) ....: -L/soft/libraries/alcf/current/xl/LAPACK/lib -llapack \
#                        -L/soft/libraries/essl/current/lib64 -lesslsmpbg \
#                        -L/soft/compilers/ibmcmp-may2014/xlf/bg/14.1/bglib64 \
#                        -lxlf90_r -lxlfmath -lxlopt -lxl -Wl,-E
#      Vulcan (ESSL) ..: -L/usr/local/tools/essl/$(VERSION_LIBRARY)/lib -lesslsmpbg
#  * LINUX CLUSTER(MKL OR ACML)
#      LC (Cab) .......: -L/usr/local/tools/$(VERSION_LIBRARY)/lib -I/usr/local/tools/$(VERSION_LIBRARY)/include \
#                        -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
#      LC .............: -L/usr/tce/packages/mkl/$(VERSION_LIBRARY)/lib -I/usr/tce/packages/mkl/$(VERSION_LIBRARY)/include \
#                        -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
#      LC (OpenMP) ....: -L/usr/tce/packages/mkl/$(VERSION_LIBRARY)/lib -I/usr/tce/packages/mkl/$(VERSION_LIBRARY)/include \
#                        -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread#  LC .............: -L/usr/local/tools/$(VERSION_LIBRARY)/lib -I/usr/local/tools/$(VERSION_LIBRARY)/include \
#      NERSC ..........: -mkl (after: module swap PrgEnv-pgi PrgEnv-intel)
#      OLCF (OpenMP)...: -L$(ACML_DIR)/ifort64_mp/lib -lacml_mp -liomp5 -lifcoremt_pic -limf -lirc -lsvml
#      OLCF ...........: -L$(ACML_DIR)/ifort64/lib -lacml -lifcoremt_pic -limf -lirc -lsvml
#  * LINUX DESKTOP
#      local ..........: -L$(HOME)/local/lib -llapack_LINUX -lblas_LINUX
ifeq ($(LAPACK_LIBRARY),TRUE)
      LINEAR_ALGEBRA = -L$(HOME)/local -llapack -lblas 
endif
ifeq ($(PORTABLE),TRUE)
      LINEAR_ALGEBRA =
endif
ifeq ($(LAPACK_LIBRARY),TRUE)
      PORTABLE = FALSE
endif

# Provide names and path to local implementation of ScaLAPAC
ifeq ($(USE_MPI),1)
      SCALAPACK_DIR = $(HOME)/local
      SCALAPACK     = -L$(SCALAPACK_DIR) -lscalapack
endif

#======================================================================#
#  Nothing beyond this line should be changed, in principle            #
#======================================================================#

# Name of the main executable and object file
HFODD_EXE     = hf$(VERSION_MAIN)
HFODD_SOURCE  = hf$(VERSION_MAIN).f
HFODD_OBJ     = hf$(VERSION_MAIN).o

# Object files
HFODD_SIZES_OBJ      = hfodd_sizes_$(VERSION_SIZES).o
HFODD_MODULES_OBJ    = hfodd_modules_$(VERSION_MODULES).o
HFODD_SHELL_OBJ      = hfodd_shell_$(VERSION_SHELL).o
HFODD_HFBTHO_OBJ     = hfodd_hfbtho_$(VERSION_HFBTHO).o
HFODD_INTERFACE_OBJ  = hfodd_interface_$(VERSION_INTERFACE).o
HFODD_FUNCTIONAL_OBJ = hfodd_functional_$(VERSION_FUNCTIONAL).o
HFODD_PAIRS_OBJ      = hfodd_pairs_$(VERSION_PAIRS).o
HFODD_FISSION_OBJ    = hfodd_fission_$(VERSION_FISSION).o
HFODD_PNP_OBJ        = hfodd_pnp_$(VERSION_PNP).o
HFODD_LIPCORR_OBJ    = hfodd_lipcorr_$(VERSION_LIPCORR).o
ifeq ($(USE_FITS),1)
      HFODD_FITS_OBJ = hfodd_fits_$(VERSION_FITS).o
else
      HFODD_FITS_OBJ =
endif
LINPACK_OBJ          = linpack.110112.o
ifeq ($(USE_MPI),1)
      HFODD_MPIIO_OBJ      = hfodd_mpiio_$(VERSION_MPIIO).o
      HFODD_MPIMANAGER_OBJ = hfodd_mpimanager_$(VERSION_MPIMANAGER).o
      ifeq ($(USE_SCALAPACK),1)
            HFODD_SL_OBJ   = hfodd_SLsiz_$(VERSION_SL).o
      endif
endif

# The user can build his/her own versions of BLAS and LAPACK by passing
# the option PORTABLE=TRUE to the compiler. Performance is likely to be
# poor compared to system libraries
ifeq ($(PORTABLE),TRUE)
      LAPACK_OBJ = lapack.110112.o
      BLAS_OBJ   = blas.110112.o
else
      LAPACK_OBJ =
      BLAS_OBJ   =
endif

# Defining compiler options for: PORTLAND FORTRAN COMPILER (pgf90)
ifeq ($(COMPILER),PGI)

      PRECISION    = -Mr8 -Mr8intrinsics=float
      FORMAT_F77   = -fpic -Mfixed
      FORMAT_F90   = -fpic
      STATIC       =
      PREPROCESSOR = -Mpreprocess -DUSE_OPENMP=$(USE_OPENMP) -DUSE_MPI=$(USE_MPI) -DUSE_MANYCORES=$(USE_MANYCORES) \
                     -DSWITCH_QUAD=$(SWITCH_QUAD) -DSWITCH_PORT=$(SWITCH_PORT) \
                     -DSWITCH_DIAG=$(SWITCH_DIAG) -DSWITCH_VECT=$(SWITCH_VECT) -DUSE_FITS=$(USE_FITS) \
                     -DUSE_SCALAPACK=$(USE_SCALAPACK) -DM_GRID=$(M_GRID) -DN_GRID=$(N_GRID) \
                     -DSWITCH_COULEX=$(SWITCH_COULEX) \
                     -DSWITCH_ESSL=$(SWITCH_ESSL)

      ifeq ($(STATIC_LIB),TRUE)
            STATIC = -Bstatic
      endif

      ifeq ($(DEBUG),FALSE)
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -fast
      else
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -g -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -Minfo=all
      endif

      ifeq ($(USE_OPENMP),1)
            OPTIONS_ = $(OPTIONS_FC) -mp
      else
            OPTIONS_ = $(OPTIONS_FC)
      endif

      ifeq ($(LARGE_MEMORY),TRUE)
            OPTIONS = $(OPTIONS_) -mcmodel=medium
      else
            OPTIONS = $(OPTIONS_)
      endif

endif

# Defining compiler options for: IFORT FORTRAN COMPILER (ifort)
ifeq ($(COMPILER),IFORT)

      PRECISION    = -i4 -r8
      FORMAT_F77   = -fPIC -fixed -80
      FORMAT_F90   = -fPIC -free -extend_source
      STATIC       =
      PREPROCESSOR = -cpp -DUSE_OPENMP=$(USE_OPENMP) -DUSE_MPI=$(USE_MPI) -DUSE_MANYCORES=$(USE_MANYCORES) \
                     -DSWITCH_QUAD=$(SWITCH_QUAD) -DSWITCH_PORT=$(SWITCH_PORT) \
                     -DSWITCH_DIAG=$(SWITCH_DIAG) -DSWITCH_VECT=$(SWITCH_VECT) -DUSE_FITS=$(USE_FITS) \
                     -DUSE_SCALAPACK=$(USE_SCALAPACK) -DM_GRID=$(M_GRID) -DN_GRID=$(N_GRID) \
                     -DSWITCH_COULEX=$(SWITCH_COULEX) \
                     -DSWITCH_ESSL=$(SWITCH_ESSL)

      ifeq ($(STATIC_LIB),TRUE)
            STATIC = -static
      endif

      ifeq ($(DEBUG),FALSE)
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -O2
      else
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -g -check all -warn nodeclarations -warn nounused
      endif

      ifeq ($(USE_OPENMP),1)
            OPTIONS_ = $(OPTIONS_FC) -qopenmp
      else
            OPTIONS_ = $(OPTIONS_FC)
      endif

      ifeq ($(LARGE_MEMORY),TRUE)
            OPTIONS = $(OPTIONS_) -mcmodel=medium -shared-intel
      else
            OPTIONS = $(OPTIONS_)
      endif

endif

# Defining compiler options for: GNU FORTRAN COMPILER (gfortran)
ifeq ($(COMPILER),GFORTRAN)

      PRECISION    = -fdefault-real-8 -fdefault-double-8
      FORMAT_F77   = -fpic -ffixed-form
      FORMAT_F90   = -fpic -ffree-form -ffree-line-length-none
      STATIC       =
      PREPROCESSOR = -cpp -DUSE_OPENMP=$(USE_OPENMP) -DUSE_MPI=$(USE_MPI) -DUSE_MANYCORES=$(USE_MANYCORES) \
                     -DSWITCH_QUAD=$(SWITCH_QUAD) -DSWITCH_PORT=$(SWITCH_PORT) \
                     -DSWITCH_DIAG=$(SWITCH_DIAG) -DSWITCH_VECT=$(SWITCH_VECT) -DUSE_FITS=$(USE_FITS) \
                     -DUSE_SCALAPACK=$(USE_SCALAPACK) -DM_GRID=$(M_GRID) -DN_GRID=$(N_GRID) \
                     -DSWITCH_COULEX=$(SWITCH_COULEX) \
                     -DSWITCH_ESSL=$(SWITCH_ESSL)

      ifeq ($(STATIC_LIB),TRUE)
              STATIC = -static
      endif

      ifeq ($(DEBUG),FALSE)
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -O2
      else
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -g -O0 -Wall \
                         -Warray-bounds -Wunderflow -Warray-temporaries \
                         -Wcharacter-truncation -Wtabs -Wintrinsic-shadow -Walign-commons -frange-check \
                         -Wconversion -Wuninitialized -pedantic \
                         -finit-real=nan \
                         -ftrapv
      endif

      ifeq ($(USE_OPENMP),1)
            OPTIONS_ = $(OPTIONS_FC) -fopenmp
      else
            OPTIONS_ = $(OPTIONS_FC)
      endif

      ifeq ($(LARGE_MEMORY),TRUE)
            OPTIONS = $(OPTIONS_) -mcmodel=medium
      else
            OPTIONS = $(OPTIONS_)
      endif

endif

# Defining compiler options for: PATHSCALE FORTRAN COMPILER (pathf90)
ifeq ($(COMPILER),PATHSCALE)

      PRECISION    = -r8
      FORMAT_F77   = -m64 -fixedform
      FORMAT_F90   = -m64 -freeform
      STATIC       =
      PREPROCESSOR = -cpp -DUSE_OPENMP=$(USE_OPENMP) -DUSE_MPI=$(USE_MPI) -DUSE_MANYCORES=$(USE_MANYCORES) \
                     -DSWITCH_QUAD=$(SWITCH_QUAD) -DSWITCH_PORT=$(SWITCH_PORT) \
                     -DSWITCH_DIAG=$(SWITCH_DIAG) -DSWITCH_VECT=$(SWITCH_VECT) -DUSE_FITS=$(USE_FITS) \
                     -DUSE_SCALAPACK=$(USE_SCALAPACK) -DM_GRID=$(M_GRID) -DN_GRID=$(N_GRID) \
                     -DSWITCH_COULEX=$(SWITCH_COULEX) \
                     -DSWITCH_ESSL=$(SWITCH_ESSL)

      ifeq ($(STATIC_LIB),TRUE)
            STATIC = -static
      endif

      ifeq ($(DEBUG),FALSE)
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -woff -O3
      else
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -fullwarn -C
      endif

      ifeq ($(USE_OPENMP),1)
            OPTIONS = $(OPTIONS_FC) -mp
      else
            OPTIONS = $(OPTIONS_FC)
      endif

endif

# Defining compiler options for: LAHEY FORTRAN COMPILER (lf95)
ifeq ($(COMPILER),LAHEY)

      PRECISION    = --dbl
      FORMAT_F77   = --fix
      FORMAT_F90   =
      STATIC       =
      PREPROCESSOR = -Cpp -DUSE_OPENMP=$(USE_OPENMP) -DUSE_MPI=$(USE_MPI) -DUSE_MANYCORES=$(USE_MANYCORES) \
                     -DSWITCH_QUAD=$(SWITCH_QUAD) -DSWITCH_PORT=$(SWITCH_PORT) \
                     -DSWITCH_DIAG=$(SWITCH_DIAG) -DSWITCH_VECT=$(SWITCH_VECT) -DUSE_FITS=$(USE_FITS) \
                     -DUSE_SCALAPACK=$(USE_SCALAPACK) -DM_GRID=$(M_GRID) -DN_GRID=$(N_GRID) \
                     -DSWITCH_COULEX=$(SWITCH_COULEX) \
                     -DSWITCH_ESSL=$(SWITCH_ESSL)

      ifeq ($(STATIC_LIB),TRUE)
            STATIC = --staticlink
      endif

      ifeq ($(DEBUG),FALSE)
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -O3
      else
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) --chk --chkglobal
      endif

      ifeq ($(USE_OPENMP),1)
            OPTIONS = $(OPTIONS_FC) --openmp
      else
            OPTIONS = $(OPTIONS_FC)
      endif

endif

# Defining compiler options for: CRAY FORTRAN COMPILER
ifeq ($(COMPILER),CRAY)

      PRECISION    = -s real64 -d p
      FORMAT_F77   = -f fixed
      FORMAT_F90   = -f free
      STATIC       =
      PREPROCESSOR = -e Z -DUSE_OPENMP=$(USE_OPENMP) -DUSE_MPI=$(USE_MPI) -DUSE_MANYCORES=$(USE_MANYCORES) \
                     -DSWITCH_QUAD=$(SWITCH_QUAD) -DSWITCH_PORT=$(SWITCH_PORT) \
                     -DSWITCH_DIAG=$(SWITCH_DIAG) -DSWITCH_VECT=$(SWITCH_VECT) -DUSE_FITS=$(USE_FITS) \
                     -DUSE_SCALAPACK=$(USE_SCALAPACK) -DM_GRID=$(M_GRID) -DN_GRID=$(N_GRID) \
                     -DSWITCH_COULEX=$(SWITCH_COULEX) \
                     -DSWITCH_ESSL=$(SWITCH_ESSL)

      ifeq ($(STATIC_LIB),TRUE)
            STATIC = --staticlink
      endif

      ifeq ($(DEBUG),FALSE)
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -O3
      else
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -e c -e D
      endif

      ifeq ($(USE_OPENMP),1)
            OPTIONS = $(OPTIONS_FC)
      else
            OPTIONS = $(OPTIONS_FC) -h noomp
      endif

endif

# Defining compiler options for: IBM FORTRAN COMPILER (xlf90)
ifeq ($(COMPILER),IBM)

      PRECISION    = -qautodbl=dbl4
      FORMAT_F77   = -qstrict -qfixed -qsuffix=cpp=f
      FORMAT_F90   = -qstrict -qfree=f90 -qsuffix=cpp=f90
      STATIC       =
      PREPROCESSOR = '-WF,-DUSE_OPENMP=$(USE_OPENMP),-DUSE_MPI=$(USE_MPI),-DUSE_MANYCORES=$(USE_MANYCORES)'\
                     '-WF,-DSWITCH_QUAD=$(SWITCH_QUAD),-DSWITCH_PORT=$(SWITCH_PORT)'\
                     '-WF,-DSWITCH_DIAG=$(SWITCH_DIAG),-DSWITCH_VECT=$(SWITCH_VECT),-DUSE_FITS=$(USE_FITS)'\
                     '-WF,-DUSE_SCALAPACK=$(USE_SCALAPACK),-DM_GRID=$(M_GRID),-DN_GRID=$(N_GRID)'\
                     '-WF,-DSWITCH_COULEX=$(SWITCH_COULEX)'\
                     '-WF,-DSWITCH_ESSL=$(SWITCH_ESSL)'

      ifeq ($(DEBUG),FALSE)
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -O3 -qhot
      else
            OPTIONS_FC = $(PREPROCESSOR) $(STATIC) -C -qflttrap
      endif

      ifeq ($(USE_OPENMP),1)
            OPTIONS = $(OPTIONS_FC) -qsmp=omp
      else
            OPTIONS = $(OPTIONS_FC)
      endif

endif

#=========================#
# Beginning of the action #
#=========================#

ifeq ($(SHOW_SOURCE),TRUE)
	build = $(HFODD_EXE) preprocessed
else
	build = $(HFODD_EXE)
endif

# Build HFODD
all: $(build)

preprocessed :
	$(FORTRAN_MPI) $(FORMAT_F90) $(PRECISION) $(OPTIONS) -E -P hfodd_hfbtho_$(VERSION_HFBTHO).f90 -c > hfodd_hfbtho_$(VERSION_HFBTHO)_PREPROCESSED.f90 $^
	$(FORTRAN_MPI) $(FORMAT_F90) $(PRECISION) $(OPTIONS) -E -P hfodd_mpimanager_$(VERSION_MPIMANAGER).f90 -c > hfodd_mpimanager_$(VERSION_MPIMANAGER)_PREPROCESSED.f90 $^
	$(FORTRAN_MPI) $(FORMAT_F90) $(PRECISION) $(OPTIONS) -E -P hfodd_fission_$(VERSION_FISSION).f90 -c > hfodd_fission_$(VERSION_FISSION)_PREPROCESSED.f90 $^
	$(FORTRAN_MPI) $(FORMAT_F77) $(PRECISION) $(OPTIONS) -E -P hf$(VERSION_MAIN).f -c > hf$(VERSION_MAIN)_PREPROCESSED.f $^

# Portable linear algebra (1/2): calculation of the determinant
$(LINPACK_OBJ) : linpack.110112.f
	$(FORTRAN_MPI) $(FORMAT_F77) $(OPTIONS) -c $<

# Portable linear algebra (2/2): set and activated only if PORTABLE is set to TRUE
$(LAPACK_OBJ) : lapack.110112.f
	$(FORTRAN_MPI) $(FORMAT_F77) $(OPTIONS) -c $<
$(BLAS_OBJ) : blas.110112.f
	$(FORTRAN_MPI) $(FORMAT_F77) $(OPTIONS) -c $<

# Extra modules (1/2): HFBTHO package and interface, MPI framework, Scalapack
$(HFODD_HFBTHO_OBJ) : hfodd_hfbtho_$(VERSION_HFBTHO).f90 $(LAPACK_OBJ) $(BLAS_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(PRECISION) $(OPTIONS) -c $^
$(HFODD_INTERFACE_OBJ) : hfodd_interface_$(VERSION_INTERFACE).f90
	$(FORTRAN_MPI) $(FORMAT_F90) $(PRECISION) $(OPTIONS) -c $<
$(HFODD_FUNCTIONAL_OBJ) : hfodd_functional_$(VERSION_FUNCTIONAL).f90 $(HFODD_HFBTHO_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(PRECISION) $(OPTIONS) -c $^
$(HFODD_SL_OBJ) : hfodd_SLsiz_$(VERSION_SL).f90
	$(FORTRAN_MPI) $(FORMAT_F90) $(PRECISION) $(OPTIONS) -c $<
$(HFODD_MPIIO_OBJ) : hfodd_mpiio_$(VERSION_MPIIO).f90
	$(FORTRAN_MPI) $(FORMAT_F90) $(PRECISION) $(OPTIONS) -c $<
$(HFODD_MPIMANAGER_OBJ) : hfodd_mpimanager_$(VERSION_MPIMANAGER).f90 $(HFODD_MPIIO_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(PRECISION) $(OPTIONS) -c $^

# Extra modules (2/2): modules having dependencies on size of HFODD arrays
$(HFODD_SIZES_OBJ) : hfodd_sizes_$(VERSION_SIZES).f90
	$(FORTRAN_MPI) $(FORMAT_F90) $(PRECISION) $(OPTIONS) -c $<
$(HFODD_MODULES_OBJ) : hfodd_modules_$(VERSION_MODULES).f
	$(FORTRAN_MPI) $(FORMAT_F77) $(PRECISION) $(OPTIONS) -c $<
$(HFODD_SHELL_OBJ) : hfodd_shell_$(VERSION_SHELL).f $(HFODD_SIZES_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F77) $(PRECISION) $(OPTIONS) -c $^
$(HFODD_PAIRS_OBJ) : hfodd_pairs_$(VERSION_PAIRS).f90
	$(FORTRAN_MPI) $(FORMAT_F90) $(PRECISION) $(OPTIONS) -c $<
$(HFODD_FISSION_OBJ) : hfodd_fission_$(VERSION_FISSION).f90 $(HFODD_SIZES_OBJ) $(HFODD_MODULES_OBJ) $(HFODD_PAIRS_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(PRECISION) $(OPTIONS) -c $^
$(HFODD_PNP_OBJ) : hfodd_pnp_$(VERSION_PNP).f90  $(HFODD_SIZES_OBJ) $(HFODD_MODULES_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(PRECISION) $(OPTIONS) -c $^
$(HFODD_LIPCORR_OBJ) : hfodd_lipcorr_$(VERSION_LIPCORR).f90  $(HFODD_SIZES_OBJ) $(HFODD_MODULES_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(PRECISION) $(OPTIONS) -c $^
$(HFODD_FITS_OBJ) : hfodd_fits_$(VERSION_FITS).f90 $(HFODD_SIZES_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F90) $(PRECISION) $(OPTIONS) -c $^

# HFODD (1/3): object file
$(HFODD_OBJ) : $(HFODD_SOURCE) $(HFODD_MODULES_OBJ) $(HFODD_SIZES_OBJ) \
               $(HFODD_SL_OBJ) $(HFODD_MPIIO_OBJ) $(HFODD_MPIMANAGER_OBJ) \
               $(HFODD_FUNCTIONAL_OBJ) $(HFODD_HFBTHO_OBJ) $(HFODD_INTERFACE_OBJ) \
               $(HFODD_SHELL_OBJ) $(HFODD_FISSION_OBJ) $(HFODD_PAIRS_OBJ) \
               $(HFODD_PNP_OBJ) $(HFODD_LIPCORR_OBJ) $(HFODD_FITS_OBJ)
	$(FORTRAN_MPI) $(FORMAT_F77) $(PRECISION) $(OPTIONS) -c $^

# HFODD (2/3): all object files from Fortran modules
GLOBAL_OBJECTS = $(HFODD_SIZES_OBJ) $(HFODD_MODULES_OBJ) \
                 $(LINPACK_OBJ) $(LAPACK_OBJ) $(BLAS_OBJ) \
                 $(HFODD_SL_OBJ) $(HFODD_MPIIO_OBJ) $(HFODD_MPIMANAGER_OBJ) \
                 $(HFODD_FUNCTIONAL_OBJ) $(HFODD_HFBTHO_OBJ) $(HFODD_INTERFACE_OBJ) \
                 $(HFODD_SHELL_OBJ) $(HFODD_FISSION_OBJ) $(HFODD_PAIRS_OBJ) \
                 $(HFODD_PNP_OBJ) $(HFODD_LIPCORR_OBJ) $(HFODD_FITS_OBJ)

# HFODD (3/3): executable
$(HFODD_EXE) : $(HFODD_OBJ) $(GLOBAL_OBJECTS)
	$(FORTRAN_MPI) $(FORMAT_F77) $(PRECISION) $(OPTIONS) -o $@ $^ $(SCALAPACK) $(LINEAR_ALGEBRA)

# Cleaning
clean ::
	-rm -f *.o *.i *.lst *.oo *.ipo *.ipa *.mod *_PREPROCESSED*
