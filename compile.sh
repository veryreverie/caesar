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
cflags="$cflags -W -Wall"       # turn on compiler warnings
cflags="$cflags -J${mdir}/"     # set .mod files to exist in mod/
cflags="$cflags -fmax-errors=1" # make compilation stop on first error

# copy shell and python scripts
cp $sdir/harmonic/*.sh $bdir/
cp $sdir/quadratic/*.sh $bdir/
cp $sdir/quadratic/*.py $bdir/

# list programs
# programs should be added so they are to the right of their dependencies
programs=(constants utils linear_algebra rand_no_gen file process moller_plesset)

harmonic_programs=(combine_forces compare_kpoints construct_finite_displacement construct_matrix_force_cnsts construct_supercell equilibrium_frac min_images symmetry fourier_interpolation generate_kgrid generate_supercell_kpoint_mesh_qe generate_supercells lte)

quadratic_programs=(mapping band_folding calculate_anharmonic calculate_bs calculate_gap generate_amplitudes generate_quadratic_configurations generate_sc_path quadratic_spline vscf_1d anharmonic)

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

# compile caesar sub-programs
#for program in ${harmonic_programs[@]}; do
#  if [ ! "$program" = "combine_forces" ]; then
#    $cf90 $odir/$program.o $objs $cflags -o $bdir/$program -lblas -llapack
#  fi
#done
#
#for program in ${quadratic_programs[@]}; do
#  $cf90 $odir/$program.o $objs $cflags -o $bdir/$program -lblas -llapack
#done

# compile main caesar program
program=caesar
$cf90 $odir/$program.o $cflags -o$bdir/$program $objs -lblas -llapack
