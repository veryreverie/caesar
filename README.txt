Caesar; a utility for calculating the vibrational free energy of periodic crystals.
Copyright (C) 2021 Mark Johnson

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

----------------------------------------
Purpose
----------------------------------------
Caesar calculates the vibrational free energy of periodic crystals, using the vibrational self-consistent field approximation.

A description of what Caesar does and what it is for can be found at <https://e-cam.readthedocs.io/en/latest/Electronic-Structure-Modules/modules/caesar/readme.html>.

----------------------------------------
Installation
----------------------------------------
To install Caesar, run:

  cmake path_to_src
  make

This should be run from the directory where caesar should be built.
'path_to_src' should be the path to the caesar/src/ directory.
The executable 'caesar' will be placed in a 'bin' directory inside the build directory.

The compiler can be specified using the flag -DCMAKE_Fortran_COMPILER, e.g.

  cmake path -DCMAKE_Fortran_COMPILER:PATH=gfortran
  cmake path -DCMAKE_Fortran_COMPILER:PATH=ifort
  cmake path -DCMAKE_Fortran_COMPILER:PATH=nagfor
  cmake path -DCMAKE_Fortran_COMPILER:PATH=pgfortran

Caesar has been tested with the following compiler versions:

gfortran 10.1.0

There is a known gfortran bug affecting shared counters. To test for this bug, call 'caesar check_counter' or run the unit tests. If this bug is present, the behaviour can be corrected by setting the CMake variable CORRECT_COUNTER_BUG (e.g. by using the CMake command line argument -DCORRECT_COUNTER_BUG:LOGICAL=true). If not fixed, this bug will likely manifest as a "Too many open files" error.

Gfortran also fails to handle finalisation / deallocation correctly, which may lead to catastrophic memory leaks when running the longer parts of Caesar.

----------------------------------------
Spglib
----------------------------------------
By default, Caesar requires the library 'spglib'. The spglib 'lib' directory must be on LIB, and the spglib 'include/spglib' directory must be on PATH in order for compilation to succeed.

Caesar requires the C API of spglib, which can be obtained as:

   git clone https://github.com/spglib/spglib.git
   cd spglib
   mkdir _build
   cd _build
   cmake .. -DCMAKE_INSTALL_PREFIX=""
   make
   mkdir _install
   make DESTDIR=_install install

spglib/_build/_install/lib should then be placed on LIB, and spglib/_build/_install/include should then be placed on PATH.

----------------------------------------
BLAS/LAPACK
----------------------------------------
By default Caesar will use CMake's find_package(LAPACK) to search for a BLAS/LAPACK distribution.

If multiple distributions are present, the specific distribution can be selected using the CMake variables BLAS_DIR, BLA_VENDOR and LAPACK_DIR.

For build systems where BLAS/LAPACK are bundled with the compiler and do not need finding separately (e.g. when using the Cray Compiler Environment) the BLAS/LAPACK finder can be disabled by setting the CMake variable FIND_LAPACK to false, e.g. with the CMake command line argument

   -DFIND_LAPACK:LOGICAL=false

----------------------------------------
Compiling without Spglib and BLAS/LAPACK
----------------------------------------
It is possible to suppress the dependency on Spglib or BLAS/LAPACK, by setting LINK_TO_SPGLIB or LINK_TO_LAPACK respectively to false, e.g. with the CMake command line arguments

   -DLINK_TO_SPGLIB:LOGICAL=false
   -DLINK_TO_LAPACK:LOGICAL=false

Disabling spglib will disable symmetry finding, which is required for setup_harmonic.

Disabling BLAS/LAPACK will disable linear algebra, which is required for most Caesar modes which are not prefixed with run_ or plot_.

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

For inputs which are file paths, e.g. run_script, relative paths will be converted into absolute paths from the working directory, unless the argument is given in an input file, in which case relative paths will be converted into absolute paths from the directory in which the input file is found.

----------------------------------------
Documentation
----------------------------------------
html documentation can be generated using Ford <https://github.com/Fortran-FOSS-Programmers/ford>. This should be run from the doc directory, as

  ford caesar.md

documentation will be generated in the doc/ford directory, and can be viewed using any html reader (e.g. a web browser).

N.B. If Ford returns an error regarding `contains` statements in `submodule procedures` then Ford pull request #321 is needed <https://github.com/Fortran-FOSS-Programmers/ford/pull/321>.

----------------------------------------
Visualisation
----------------------------------------
Caesar uses python and the matplotlib plotting library. This can be located using the python_path keyword when running plot functionality.

The harmonic calculation produces a castep .phonon file, which can be viewed using JMol or similar.

----------------------------------------
Unit tests
----------------------------------------
Unit tests are compiled by default. This requires version 4.2.0 of the pFUnit library <https://github.com/Goddard-Fortran-Ecosystem/pFUnit>. The /bin directory of pFUnit (and its subsidiaries GFTL, GFTL_SHARED and FARGPARSE if present) must be on PATH. To suppress compilation of tests (and remove the dependency on pFUnit), set ENABLE_TESTS to false when running CMake, e.g. using

  -DENABLE_TESTS:LOGICAL=false

Unit tests can be run by calling

  ctest

from the build directory where CMake was run.
