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
cflags="$cflags -g"             # turn on debugging
cflags="$cflags -O0"            # turn off optimisation
#cflags="$cflags -std=f95"       # force standards compliance
cflags="$cflags -W -Wall"       # turn on compiler warnings
cflags="$cflags -pedantic"      # turn on pedantic mode
cflags="$cflags -J${mdir}/"     # set .mod files to exist in mod/
cflags="$cflags -fmax-errors=1" # make compilation stop on first error

# copy shell and python scripts
cp $sdir/quadratic/*.py $bdir
cp $sdir/caesar $bdir
chmod u+x $bdir/caesar

# list programs
# programs should be added so they are to the right of their dependencies
programs=(constants string file utils linear_algebra algebra supercell rand_no_gen moller_plesset structure dft_output_file structure_to_dft calculate_symmetry_helper bands displacement_patterns)

harmonic_programs=(group calculate_symmetry_group unique_directions construct_supercell min_images fourier_interpolation_symmetry fourier_interpolation generate_supercells lte hartree_to_eV setup_harmonic lte_harmonic)

quadratic_programs=(mapping calculate_anharmonic calculate_gap quadratic_spline vscf_1d anharmonic bs_quadratic setup_quadratic)

testing_programs=(test_copy_harmonic)

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
$cf90 $odir/$program.o $cflags -o$bdir/.$program $objs -lblas -llapack
