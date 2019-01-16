CAESAR

----------------------------------------
Installation
----------------------------------------
To install Caesar, run:

  cmake path
  make

This should be run from the directory where cmake should be built.
'path' should be the path to the CAESAR/src/ directory.
The executable 'caesar' will be placed in a 'bin' directory.

The compiler can be specified using the flag -DCMAKE_Fortran_COMPILER, e.g.

  cmake path -DCMAKE_Fortran_COMPILER:PATH=gfortran
  cmake path -DCMAKE_Fortran_COMPILER:PATH=ifort
  cmake path -DCMAKE_Fortran_COMPILER:PATH=nagfor
  cmake path -DCMAKE_Fortran_COMPILER:PATH=pgfortran

Caesar has been tested with the following compiler versions:

gfortran 5.4.0, 7.3.0
nagfor 6.2 (release 6207)

There is a known gfortran bug affecting shared counters. To test for this bug, call 'caesar check_counter'. If this bug is present, the behaviour can be corrected by setting the CMake variable CORRECT_COUNTER_BUG (e.g. by using the command line arguments -DCORRECT_COUNTER_BUG:LOGICAL=true). If not fixed, this bug will likely manifest as a "Too many open files" error.

----------------------------------------
Spglib
----------------------------------------
By default, Caesar requires the library 'spglib'. The spglib 'lib' directory must be on LIB, and the spglib 'include/spglib' directory must be on PATH in order for compilation to succeed.

Caesar requires the C API of spglib, which can be obtained as:

   git clone https://github.com/atztogo/spglib.git
   cd spglib
   mkdir _build
   cd _build
   cmake ..
   make
   mkdir _install
   make DESTDIR=_install install

spglib/_build/_install/lib should then be placed on LIB, and spglib/_build/_install/include should then be places on PATH.

----------------------------------------
BLAS/LAPACK
----------------------------------------
By default Caesar will use CMake's find_package(LAPACK) to search for a BLAS/LAPACK distribution.

If multiple distributions are present, the specific distribution can be selected using the CMake variables BLAS_DIR, BLA_VENDOR and LAPACK_DIR.

For build systems where BLAS/LAPACK are bundled with the compiler and do not need finding separately (e.g. when using the Cray Compiler Environment) the BLAS/LAPACK finder can be disabled by setting the CMake variable FIND_LAPACK to false, e.g. with the CMake command line argument

   -DFIND_LAPACK:LOGICAL=false

----------------------------------------
ARPACK
----------------------------------------
By default Caesar will look for libarpack.a (or libarpack.so etc.) on LIB.

If multiple distributions are present, the first found on LIB will be used.

----------------------------------------
Compiling without Spglib and BLAS/LAPACK
----------------------------------------
It is possible to suppress the requirement for any of Spglib, BLAS/LAPACK and ARPACK, by setting LINK_TO_SPGLIB, LINK_TO_LAPACK or LINK_TO_ARPACK respectively to false, e.g. with the CMake command line arguments

   -DLINK_TO_SPGLIB:LOGICAL=false
   -DLINK_TO_LAPACK:LOGICAL=false
   -DLINK_TO_ARPACK:LOGICAL=false

N.B. ARPACK requires BLAS/LAPACK, so if LINK_TO_LAPACK is set to false then LINK_TO_ARPACK must also be set to false.

This will disable symmetry finding, linear algebra and the Lanczos algorithm respectively. Symmetry finding (and thus spglib) is required for setup_harmonic, linear algebra (and thus BLAS/LAPACK) is required for most caesar modes not prfixed with run_ or plot_, and the Lanczos algorithm (and thus ARPACK) is required for calculate_states and calculate_anharmonic_observables.

The 'run' modes, run_harmonic and run_anharmonic, do not require spglib or BLAS/LAPACK, so not linking to either can be useful if compiling Caesar on a cluster where they are not available, with the intention of running setup and processing steps elsewhere.

----------------------------------------
Quip
----------------------------------------
Caesar can be linked against Quip by setting LINK_TO_QUIP, e.g. with the CMake command line argument

  -DLINK_TO_QUIP:LOGICAL=true

The directory containing libquip.so and libatoms.so (or the equivalent .a etc. files) must be on LIB and PATH. These libraries must have been compiled using the same compiler used to compile Caesar.

If multiple distributions are present, the first found on LIB will be used.

Quip requires BLAS/LAPACK, so PATH_TO_QUIP should not be specified if linking against BLAS/LAPACK has been suppressed.

N.B. Quip requires atomic numbers, which are read from .cell files as the label on the atomic species, so e.g. titanium must be written as Ti:22.

----------------------------------------
Running
----------------------------------------
For help, call

  caesar -h

----------------------------------------
Inputs
----------------------------------------
Caesar can be given arguments on the command line, from a file or interactively, as explained in the helptext provided by calling caesar -h.

Caesar also requires an input file, and at present only CASTEP .cell files are supported. This .cell file must include a lattice block, an atomic positions block, and a species mass block. This .cell file should be in the working directory specified when Caesar is called.

Arguments have the same names and syntax regardless of how they are supplied, so e.g. a 3x3x3 q-point grid can be specified on the command line as:
  caesar setup_harmonic --q-point_grid 3 3 3
or in a file, by calling:
  caesar setup_harmonic -f filename
and adding the following line to the file:
  q-point_grid 3 3 3
or interactively, as:
  caesar setup_harmonic -i
and inputing "3 3 3" when prompted.

In Caesar input files:
   - blank lines are ignored.
   - '!' is a comment character, and everything after this up to the end of the line will be ignored.
   - The keyword is taken to be the first non-whitespace string on each line, ending at the next whitespace. Everything else on the line, up to any comment characters, is taken to be the argument of that keyword.
     e.g. the line 'key word a b c' will parse as the keyword 'key' followed by the argument 'word a b c'

An error will be thrown if any unrecognised keywords are found, or if a keyword appears more than once from the same source (command line / file). Command line options override file options where both are present.
A file containing the used values of each option will be produced at the start of each Caesar run. This file can be used as the input settings if the calculation is repeated.
An example input file can be found in doc/input_files.

In order to calculate energies, Caesar requires a run script. This script will be called repeatedly from the working directory, and will be passed a number of arguments on the command line each time it is called.
The return code of the run script will be echoed to the output of Caesar after each run.
An example run script can be found in doc/input_files.

----------------------------------------
Documentation
----------------------------------------
html documentation can be generated using either Doxygen or Ford. Both should be run from the doc directory. Doxygen documentation can be generated by running
  doxygen caesar.Doxyfile
and Ford documentation can be generated by running
  ford caesar.md
documentation will be generated in doc/doxygen or doc/ford respectively, and can be viewed using any html reader (e.g. a web browser).

----------------------------------------
Visualisation
----------------------------------------
Caesar uses python and the matplotlib plotting library. This can be located using the python_path keyword when running plot functionality.

The harmonic calculation produces a castep .phonon file, which can be viewed using JMol or similar.
