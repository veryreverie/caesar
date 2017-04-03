# This is a placeholder compile script, to be replaced by
# CMake or Make in due course

# clean build directories if requested
if [ "$1" = "clean" ]; then
  rm -r bin
  rm -r mod
  rm -r obj
  exit
fi

# Make directories which may not have been copied by git
directories=(bin doc mod obj src)
for directory in ${directories[@]}; do
	if [ ! -d $directory ]; then
		mkdir $directory
	fi
done

# define directories
bdir=bin
ddir=doc
mdir=mod
odir=obj
sdir=src

# define compiler and flags
cf90=gfortran
cflags=""
cflags="$cflags -g"               # Turn on debugging.
cflags="$cflags -O0"              # Turn off optimisation.
#cflags="$cflags -std=f95"         # Force standards compliance.
#cflags="$cflags -std=f2003"       # Force standards compliance.
cflags="$cflags -W -Wall -Wextra" # Turn on compiler warnings.
cflags="$cflags -pedantic"        # Turn on pedantic mode.
cflags="$cflags -J${mdir}/"       # Set .mod files to exist in mod/.
cflags="$cflags -fmax-errors=1"   # Make compilation stop on first error.
cflags="$cflags -fcheck=all"      # Turn on run-time checks.
cflags="$cflags -Wno-array-temporaries" # Turn off temp array warning.
#cflags="$cflags -ffpe-trap=invalid,zero,overflow,underflow,denormal"

cc=gcc

# Copy python script.
cp $sdir/quadratic/*.py $bdir

# list programs
# programs should be added so they are to the right of their dependencies
programs=(constants string file utils linear_algebra algebra group supercell kpoints moller_plesset structure dft_input_file dft_output_file structure_to_dft calculate_symmetry bands displacement_patterns)

harmonic_programs=(calculate_symmetry_group unique_directions construct_supercell min_images generate_supercells lte hartree_to_eV setup_harmonic run_harmonic lte_harmonic)

quadratic_programs=(mapping calculate_anharmonic calculate_gap quadratic_spline vscf_1d setup_quadratic run_quadratic anharmonic bs_quadratic)

testing_programs=(atom_mapping test_copy_harmonic test_lte test_copy_quadratic)

# Compile c system call.
$cc -c $sdir/system.c -W -Wall -Wextra -pedantic -std=c99 -o$odir/system.o

# Compile compiler-specific module.
if [ "$cf90" = "gfortran" ]; then
  $cf90 -c $sdir/gfortran_specific.f90 $cflags -o$odir/compiler_specific.o
fi

cflags="$cflags -std=f2003"       # Force standards compliance.

# compile objects
for program in ${programs[@]}; do
  $cf90 -c $sdir/$program.f90 $cflags -o$odir/$program.o
done

for program in ${harmonic_programs[@]}; do
  $cf90 -c $sdir/harmonic/$program.f90 $cflags -o$odir/$program.o
done

for program in ${quadratic_programs[@]}; do
  $cf90 -c $sdir/quadratic/$program.f90 $cflags -o$odir/$program.o
done

for program in ${testing_programs[@]}; do
  $cf90 -c $sdir/testing/$program.f90 $cflags -o$odir/$program.o
done

program=caesar
$cf90 -c $sdir/$program.f90 $cflags -o$odir/$program.o

# list objects for linking
objs=()

objs="$odir/system.o $objs"
objs="$odir/compiler_specific.o $objs"

for program in ${programs[@]}; do
  objs="$odir/$program.o $objs"
done

for program in ${harmonic_programs[@]}; do
  objs="$odir/$program.o $objs"
done

for program in ${quadratic_programs[@]}; do
  objs="$odir/$program.o $objs"
done

for program in ${testing_programs[@]}; do
  objs="$odir/$program.o $objs"
done

# compile main caesar program
program=caesar
$cf90 $odir/$program.o $cflags -o$bdir/$program $objs -lblas -llapack
