# Compilation and Execution #

###  Setting the Makefile ###

HFBTHO requires a Fortran compiler, a modern version of the GNU Makefile
utility (>= make-4.0), and an implementation of the standard Open Source
libraries BLAS and LAPACK (http://www.netlib.org/lapack). Before building the
code, you **must** adapt the Makefile to your system. This includes looking at
5 categories of options described below.

#### a) Essential Options ####

These options define the compiler.

  - ***COMPUTER***: This variable mostly helps with setting up the path to the
    BLAS and LAPACK libraries. Most users should set this variable to
    ***DESKTOP***  and refer to the corresponding part of the Makefile to
    adjust the variable ***LINEAR_ALGEBRA***.
    Possible values: ***DESKTOP***, ***QUARTZ***, ***LASSEN***;

  - ***COMPILER***: This variable simply activates the relevant list of
    command-line flags for the requested compiler;
    Possible values: ***IFORT***, ***GFORTRAN***, ***CRAY***, ***IBM***.

  - ***FORTRAN_MPI***: This variables specifies the name of the executable for
    the compiler. For example, for the GNU compiler, it could be
    ***FORTRAN_MPI = gfortran***. When compiling in parallel mode using a
    standard MPI installation, it could be ***FORTRAN_MPI = mpif90***, etc.

#### b) Build Options ####

These options are used to generate additional information about the code or the
build process.

  - ***DEBUG***: Depending on the type of compiler set by the variable
    ***COMPILER***, this variable changes the preset compiler flags to activate
    various debugging options.
    Possible values: ***0*** (optimization), ***1*** (debug);

  - ***DOCUMENTATION***: This variable is used to build the DoxyGen
	documentation.
    Possible values: ***0*** (no documentation), ***1*** (documentation);

  - ***VERBOSE***: This variable simply toggles on/off the output of the
    compiler.
    Possible values: ***0*** (condensed compilation), ***1*** (verbose
    compilation).

#### c) HPC Options ####

These options activate either shared-memory parallelism (OpenMP) or
distributed-memory parallelism (MPI).

  - ***USE_OPENMP***: This variable is used to activate multithreading
    throughout the code. This activates explicitly-parallelized sections of the
    code and can accelerate execution significantly.
    Possible values: ***0*** (no multithreading), ***1*** (multithreading);

  - ***USE_MPI***: This variable activates distributed-memory parallelism
    through the MPI library. If ***USE_MPI > 0***, the code must be compiled
    with an MPI-enabled compiler, e.g., a wrapper such as mpif90 for example.
    There are two possible use-case of MPI: (i) AMP, PNP and AMP\&PNP
    calculations distribute the calculation of rotation and/or gauge angles
    across available MPI processes. This case is activated by setting
    ***USE_MPI=1***, (ii) large-scale mass tables, drip-lines or PES
    calculations launch a series of different HFBTHO calculations across
    available MPI processes. This option is activated by setting
    ***USE_MPI=2***.
    Possible values: ***0*** (no MPI), ***1*** (MPI for a single HFBTHO),
    ***2*** (for parallel mass tables, driplines or PES calculations).

#### d) Physics Options ####

These options provides additional capabilities built on top of the HFBTHO
solver.

  - ***DRIP_LINES***: Calculates a mass table from dripline to dripline.
    Possible values: ***0*** (inactive), ***1*** (active);

  - ***DO_MASSTABLE***: Calculates a section of the mass table.
    Possible values:  ***0*** (inactive), ***1*** (active);

  - ***DO_PES***: Calculates a potential energy surface for a set of nuclei.
    Possible values: ***0*** (inactive), ***1*** (active);

  - ***READ_FUNCTIONAL***: Reads the parameters of the functional from a file.
    By default, the energy functional is defined in the code based on the value
    of the input keyword 'functional'.
    Possible values: ***0***, ***1***.

#### e) Expert Options ####

These options are mostly for expert users.

  - ***USE_ESSL***: On IBM machines, it is usually more efficient to link to
    the ESSL libraries rather than standard BLAS and LAPACK. Setting this flag
    to 1 will change the routines used to perform some of the linear algebra
    operations.
    Possible values: ***0***, ***1***;

  - ***GOGNY_SYMMETRIES***: Assumes that two-body matrix elements of
    finite-range interactions have specific symmetries. This option is only for
    development purposes.
    Possible values: ***0***, ***1***;

  - ***GOGNY_HYPER***: Use hypergeometric functions to calculate finite-range
    matrix elements which is accurate for either large and small bases. This
    option is only for development purposes.
    Possible values: ***0***, ***1***;

  - ***USE_QRPA***: Saves the HFB solution into a file readable by the pnFAM
    QRPA code developed at UNC Chapel Hill, USA.
    Possible values: ***0***, ***1***.



### Default Configuration ###

The Makefile included with the code is preset for the GNU compiler, assumes
'gfortran' is the compiler executable and activates optimization (=deactivates
debugging). It also builds the code by linking to existing BLAS and LAPACK
libraries, that are assumed to be available under the names libblas.a and
liblapack.a and be located in $HOME/local. Refer to the instructions given in
the Makefile itself for further guidance.



### Building the Code ###

Once the Makefile has been set up, simply type

    make

to build the code.



### Running the code ###

The call sequence of HFBTHO has been simplified between version 3.00 and
version 4.00. In serial mode (*USE_MPI = 0*), the code HFBTHO is now run by

     ./hfbtho_main [input file]

where *[input file]* is the name of the input file containing the Namelists. If
this argument is not present, the code tries to open the file named
*hfbtho_NAMELIST.dat* to set the user-defined data. If this file is not
present and no input file is passed as argument, an error message will be
returned.

Three files contain the results of an HFBTHO calculation:

   - *hfbtho.out*: ASCII file containing a selection of results, namely some
     of the input data, the self-consistent iterations and various observables
     at convergence;

   - *thoout.dat*: ASCII file containing an extended selection of results. In
     addition to all the data printed in *hfbtho.out*, it also prints the
     list of quasiparticle energies and equivalent single-particle energies,
     for both protons and neutrons;

   - *hfbtho_output.hel*: binary file used to restart iterations. The
     structure of the data has been updated between versions 3.00 and 4.00.
     Backward compatibility has been maintained, i.e., version 4.00 of HFBTHO
     can read binary files generated by verison 3.00 of the code. The opposite
     is not true.

In parallel mode, the code is run with

    mpirun -np 4 ./hfbtho_main

When running the code on personal laptops and desktops with OpenMP
multi-threading (*USE_OPENMP = 1*), it is often necessary to set the following
environment variables (examples below are given for the bash shell),

     export OMP_STACKSIZE=64M
     ulimit -s unlimited

This will ensure there is enough stack memory for each thread. The value of
*OMP_STACKSIZE* may have to be increased even further for large bases.
**Without such instructions, the code may crash with a 'segmentation fault'
error message.**


