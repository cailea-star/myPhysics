
SHELL := /bin/bash
COMPILE_DIR = $(shell pwd)

# Check if compile directory exists, returns error code and exit otherwise
EXISTS = $(firstword $(wildcard $(COMPILE_DIR)) )
$(info $$COMPILE_DIR is [${COMPILE_DIR}] )
ifeq ($(EXISTS), )
      $(error ${COMPILE_DIR} directory does not exist)
endif

# Basic options
#   - COMPUTER.............: current choice is DESKTOP, QUARTZ (works for ruby), LASSEN;  used to automatically set LINEAR_ALGEBRA below
#   - COMPILER ............: type of compiler: GFORTRAN, IFORT, CRAY or IBM; used to automatically set compiler options
#   - FORTRAN_MPI .........: actual name of the Fortran compiler. Standard choices are gfortran, mpif90, ifort, xlf-gpu, etc.
COMPUTER    = DESKTOP
COMPILER    = GFORTRAN
FORTRAN_MPI = gfortran
# Build options
#   - DEBUG ...............: activates debugging options for each compiler
#   - DOCUMENTATION .......: builds doxygen documentation
#   - VERBOSE .............: enables verbose compilation
#                              0: short display
#                              1: full command displayed
DEBUG         = 0
DOCUMENTATION = 0
VERBOSE       = 1
#  HPC options (integer numbers, usually 0: inactive, 1: active unless specified otherwise)
#   - USE_OPENMP ..........: enables OpenMP multithreading
#   - USE_MPI .............: enables MPI distributed parallelism
#                              1: Multi-core mode, single HFBTHO calculation on a subcommunicator; used in standalone
#                                 applications or ANL optimization
#                              2: Parallel mass table, driplines or PES mode
USE_OPENMP = 1
USE_MPI    = 0
# Physics options
#   - DRIP_LINES ..........: 1, calculates a mass table from dripline to dripline; 0, inactive
#   - DO_MASSTABLE ........: 1, calculates a section of the mass table; 0, inactive
#   - DO_PES ..............: 1, calculates a potential energy surface for a set of nuclei; 0, inactive
#   - READ_FUNCTIONAL .....: 1, reads the parameters of the functional from a file
#                            0, energy functional is defined in the code based on the value of the input keyword 'functional'
DRIP_LINES      = 0
DO_MASSTABLE    = 0
DO_PES          = 0
READ_FUNCTIONAL = 0
# Expert options
#   - USE_ESSL ............: uses IBM ESSL library for LAPACK and BLAS (internally switches to different routine calls)
#   - GOGNY_SYMMETRIES.....: 1, for production runs; assumes several symmetries in the finite range matrix elements.
#                            0, for debugging the finite range matrix elements; makes no assumption of any symmetry in the finite range matrix elements.
#                               Requires compiler compatible with 2008 standard of up to 15-dimensional arrays.
#   - GOGNY_HYPER..........: 1, uses hypergeometric function to calculate finite range matrix elements; accurate for big and small basis size
#                            0, uses the direct Gogny transformation to calculate the finite range matrix elements; accurate for small basis size ONLY (N<20)
#   - USE_QRPA ............: 1, produces output for the QRPA-pnFAM code of UNC Chapel Hill; 0, inactive
USE_ESSL         = 0
GOGNY_SYMMETRIES = 1
GOGNY_HYPER      = 1
USE_QRPA         = 0

# Names and paths to BLAS and LAPACK Libraries and f2py command
#   - LC CLUSTER (MKL)
#       * Ruby, OpenMP ..: -L/usr/tce/packages/mkl/$(VERSION_LIBRARY)/lib -I/usr/tce/packages/mkl/$(VERSION_LIBRARY)/include \
#                          -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_avx2  -liomp
#       * Ruby, serial ..: -L/usr/tce/packages/mkl/$(VERSION_LIBRARY)/lib -I/usr/tce/packages/mkl/$(VERSION_LIBRARY)/include \
#                         -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
#   - LINUX DESKTOP
#       * local ...........: -L$(HOME)/local/lib -llapack_LINUX -lblas_LINUX
ifeq ($(COMPUTER),DESKTOP)
      F2PY = f2py3
      LINEAR_ALGEBRA = -L$(HOME)/local -llapack -L$(HOME)/local -lblas
endif
ifeq ($(COMPUTER),QUARTZ)
      F2PY = $(HOME)/local/python-3.7.2/bin/f2py3 --fcompiler=intelem
      VERSION_LIBRARY = mkl-2020.0
      LINEAR_ALGEBRA = -L/usr/tce/packages/mkl/$(VERSION_LIBRARY)/lib -I/usr/tce/packages/mkl/$(VERSION_LIBRARY)/include \
                       -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
endif
ifeq ($(COMPUTER),LASSEN)
      F2PY = f2py3
      LINEAR_ALGEBRA = -I$(ESSLHEADERDIR) -L$(ESSLLIBDIR64) -R$(ESSLLIBDIR64) -lesslsmp -llapackforessl
endif

#======================================================================#
#  Nothing beyond this line should be changed, in principle            #
#======================================================================#

# Consistency checks between MPI option, dripline, mass table or PES modes
ifeq ($(USE_MPI),2)
     n = 0
     ifeq ($(DO_MASSTABLE),0)
           n += 0
     endif
     ifeq ($(DRIP_LINES),0)
           n += 0
     endif
     ifeq ($(DO_PES),0)
           n += 0
     endif
     ifeq ($(n),0 0 0 0)
           $(error DRIP_LINES, DO_MASSTABLE or DO_PES mode must be set to 1)
     endif
endif
ifeq ($(DRIP_LINES),1)
     ifeq ($(DO_MASSTABLE),1)
           $(error DRIP_LINES and DO_MASSTABLE modes are incompatible)
     endif
     ifeq ($(DO_PES),1)
           $(error DRIP_LINES and DO_PES modes are incompatible)
     endif
     ifneq ($(USE_MPI),2)
           $(error DRIP_LINES mode requires USE_MPI=2)
     endif
endif
ifeq ($(DO_MASSTABLE),1)
     ifeq ($(DRIP_LINES),1)
           $(error DO_MASSTABLE and DRIP_LINES modes are incompatible)
     endif
     ifeq ($(DO_PES),1)
           $(error DO_MASSTABLE and DO_PES modes are incompatible)
     endif
     ifneq ($(USE_MPI),2)
           $(error DO_MASSTABLE mode requires USE_MPI=2)
     endif
endif
ifeq ($(DO_PES),1)
     ifeq ($(DO_MASSTABLE),1)
           $(error DO_PES and DO_MASSTABLE modes are incompatible)
     endif
     ifeq ($(DRIP_LINES),1)
           $(error DO_PES and DRIP_LINES modes are incompatible)
     endif
     ifneq ($(USE_MPI),2)
           $(error DO_PES mode requires USE_MPI=2)
     endif
endif

ifeq ($(DOCUMENTATION),1)
      DOXYGEN = doc
endif

# Names
HFBTHO_EXE    = hfbtho_main
HFBTHO_SOURCE = hfbtho_main.f90
HFBTHO_OBJ    = hfbtho_main.o

# Defining compiler options for: IFORT FORTRAN COMPILER (ifort)
ifeq ($(COMPILER),IFORT)
      FORMAT_F90   = -fPIC -free -extend_source
      PREPROCESSOR = -fpp -DUSE_OPENMP=$(USE_OPENMP) \
                          -DUSE_MPI=$(USE_MPI)  \
                          -DDO_PES=$(DO_PES) \
                          -DDO_MASSTABLE=$(DO_MASSTABLE) \
                          -DDRIP_LINES=$(DRIP_LINES) \
                          -DREAD_FUNCTIONAL=$(READ_FUNCTIONAL) \
                          -DGOGNY_SYMMETRIES=$(GOGNY_SYMMETRIES) \
                          -DGOGNY_HYPER=$(GOGNY_HYPER) \
                          -DUSE_QRPA=$(USE_QRPA) \
                          -DUSE_ESSL=$(USE_ESSL) $(FFLAGS)
      ifeq ($(DEBUG),0)
            OPTIONS = $(PREPROCESSOR) -O3
      else
            OPTIONS = $(PREPROCESSOR) -g -traceback -check bounds,pointers,shape,uninit -init=snan,arrays \
                     -warn nodeclarations -warn nounused
      endif
      ifeq ($(USE_OPENMP),1)
            OPTIONS += -qopenmp
      endif
      FORMAT_F2PY  = $(FORMAT_F90)
      OPTIONS_F2PY = $(OPTIONS)
endif

# Defining compiler options for: GNU FORTRAN COMPILER (gfortran)
ifeq ($(COMPILER),GFORTRAN)
      FORMAT_F90   = -fPIC -ffree-form -ffree-line-length-none
      PREPROCESSOR = -cpp -DUSE_OPENMP=$(USE_OPENMP) \
                          -DUSE_MPI=$(USE_MPI)  \
                          -DDO_PES=$(DO_PES) \
                          -DDO_MASSTABLE=$(DO_MASSTABLE) \
                          -DDRIP_LINES=$(DRIP_LINES) \
                          -DREAD_FUNCTIONAL=$(READ_FUNCTIONAL) \
                          -DGOGNY_SYMMETRIES=$(GOGNY_SYMMETRIES) \
                          -DGOGNY_HYPER=$(GOGNY_HYPER) \
                          -DUSE_QRPA=$(USE_QRPA) \
                          -DUSE_ESSL=$(USE_ESSL) $(FFLAGS)
      ifeq ($(DEBUG),0)
            OPTIONS = $(PREPROCESSOR) -O3
      else
            OPTIONS = $(PREPROCESSOR) -g -O0 -Wall -Warray-bounds -Wunderflow \
                     -Wcharacter-truncation -Wtabs -Wintrinsic-shadow -Walign-commons -frange-check \
                     -fbounds-check -Wconversion -Wuninitialized -pedantic \
                     -finit-real=nan -ffpe-trap=invalid,zero
      endif
      ifeq ($(USE_OPENMP),1)
            OPTIONS += -fopenmp
      endif
      FORMAT_F2PY  = $(FORMAT_F90)
      OPTIONS_F2PY = $(OPTIONS)
endif

# Defining compiler options for: CRAY COMPILER
ifeq ($(COMPILER),CRAY)
      FORMAT_F90   = -f free
      PREPROCESSOR = -e Z -DUSE_OPENMP=$(USE_OPENMP) \
                          -DUSE_MPI=$(USE_MPI)  \
                          -DDO_PES=$(DO_PES) \
                          -DDO_MASSTABLE=$(DO_MASSTABLE) \
                          -DDRIP_LINES=$(DRIP_LINES) \
                          -DREAD_FUNCTIONAL=$(READ_FUNCTIONAL) \
                          -DGOGNY_SYMMETRIES=$(GOGNY_SYMMETRIES) \
                          -DGOGNY_HYPER=$(GOGNY_HYPER) \
                          -DUSE_QRPA=$(USE_QRPA) \
                          -DUSE_ESSL=$(USE_ESSL) $(FFLAGS)
      ifeq ($(DEBUG),0)
            OPTIONS = $(PREPROCESSOR) -O3
      else
            OPTIONS = $(PREPROCESSOR) -e c -e D
      endif
      ifeq ($(USE_OPENMP),0)
            OPTIONS += -h noomp
      endif
      FORMAT_F2PY  = $(FORMAT_F90)
      OPTIONS_F2PY = $(OPTIONS)
endif

# Defining compiler options for: CRAY COMPILER
ifeq ($(COMPILER),IBM)
      FORMAT_F90   = -qpic -qstrict -qfree=f90 -qsuffix=cpp=f90 -qnoescape
      PREPROCESSOR = '-WF,-DUSE_OPENMP=$(USE_OPENMP)' \
                     '-WF,-DUSE_MPI=$(USE_MPI)' \
                     '-WF,-DDO_PES=$(DO_PES)' \
                     '-WF,-DDO_MASSTABLE=$(DO_MASSTABLE)' \
                     '-WF,-DDRIP_LINES=$(DRIP_LINES)' \
                     '-WF,-DREAD_FUNCTIONAL=$(READ_FUNCTIONAL)' \
                     '-WF,-DGOGNY_SYMMETRIES=$(GOGNY_SYMMETRIES)' \
                     '-WF,-DGOGNY_HYPER=$(GOGNY_HYPER)' \
                     '-WF,-DUSE_QRPA=$(USE_QRPA)' \
                     '-WF,-DUSE_ESSL=$(USE_ESSL)' $(FFLAGS)
      ifeq ($(DEBUG),0)
            OPTIONS = $(PREPROCESSOR) -O2 -qhot
      else
            OPTIONS = $(PREPROCESSOR) -g -C -qfloat=nans -qflttrap=invalid:nanq:overflow:underflow:zerodivide \
                     -qstackprotect=all -qsigtrap -qcheck=all -qfullpath
      endif
      ifeq ($(USE_OPENMP),1)
            OPTIONS += -qsmp=omp
      endif
      FORMAT_F2PY  = $(FORMAT_F90)
      OPTIONS_F2PY = $(OPTIONS)
endif

#=========================#
# Beginning of the action #
#=========================#

HFBTHO_DIR = $(COMPILE_DIR)/src

# Doxygen documentation
DOC_DIR = $(HFBTHO_DIR)/doc

# Export HFBTHO names
export HFBTHO_EXE
export HFBTHO_SOURCE
export HFBTHO_OBJ
export HFBTHO_DIR
export DOC_DIR
export USE_MPI
export READ_FUNCTIONAL

# Export main compiler options
export COMPILER
export FORTRAN_MPI
export F2PY
export FORMAT_F90
export FORMAT_F2PY
export OPTIONS
export OPTIONS_F2PY
export VERBOSE
export LINEAR_ALGEBRA

# Options are: anl, $(HFBTHO_EXE), $(DOXYGEN)
all: $(HFBTHO_EXE) $(DOXYGEN)

# Doxygen documentation
$(DOXYGEN): $(HFBTHO_EXE)
	( $(MAKE) -C $(DOC_DIR); cd .. )

# HFBTHO executable
$(HFBTHO_EXE):
	( cd $(HFBTHO_DIR); $(MAKE); cd .. )

# Cleaning
clean :
	( rm -f *.o *.oo *.ipo *.mod; cd $(HFBTHO_DIR); $(MAKE) clean; cd ..; $(MAKE) clean -C $(DOC_DIR) )

mrproper : clean
	$(MAKE) mrproper -C $(HFBTHO_DIR)
