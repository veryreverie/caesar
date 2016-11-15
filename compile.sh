# This is a placeholder compile script, to be replaced by
# CMake or Make in due course

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

# define compiler
cf90=gfortran

# copy shell and python scripts
cp $sdir/harmonic/*.sh $bdir/
cp $sdir/quadratic/*.sh $bdir/
cp $sdir/quadratic/*.py $bdir/

# list programs
harmonic_programs=(combine_forces compare_kpoints construct_finite_displacement construct_matrix_force_cnsts construct_supercell convert_forces_from_Rybohr_to_eVang equilibrium_frac fourier_interpolation generate_kgrid generate_supercell_kpoint_mesh_qe generate_supercells lte lte_lower)

quadratic_programs=(band_folding calculate_anharmonic calculate_bs calculate_gap generate_amplitudes generate_quadratic_configurations generate_sc_path quadratic_spline vscf_1d)

programs=(utils)

# list objects for linking
objs=()
for program in ${programs[@]}; do
  objs+=($odir/$program.o)
done

# compile objects
for program in ${harmonic_programs[@]}; do
  $cf90 -c $sdir/harmonic/$program.f90 -J$mdir/ -o$odir/$program.o
done

for program in ${quadratic_programs[@]}; do
  $cf90 -c $sdir/quadratic/$program.f90 -J$mdir/ -o$odir/$program.o
done

for program in ${programs[@]}; do
  $cf90 -c $sdir/$program.f90 -J$mdir/ -o$odir/$program.o
done

program=caesar
$cf90 -c $sdir/$program.f90 -J$mdir/ -o$odir/$program.o

# link objects
for program in ${harmonic_programs[@]}; do
	$cf90 $odir/$program.o -J$mdir/ -o $bdir/$program $objs -lblas -llapack
done

for program in ${quadratic_programs[@]}; do
	$cf90 $odir/$program.o -J$mdir/ -o $bdir/$program $objs -lblas -llapack
done

# compile main caesar program
program=caesar
$cf90 $odir/$program.o -J$mdir/ -o$bdir/$program $objs -lblas -llapack
