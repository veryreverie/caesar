# This is a placeholder compile script, to be replaced by
# CMake or Make in due course

# --------------------
# harmonic
# --------------------
# define directories
src=src/harmonic/programs/
mod=mod/
obj=obj/
bin=bin/

# copy shell scripts
cp src/harmonic/scripts/*.sh bin/

# list programs
programs=(combine_forces compare_kpoints construct_finite_displacement construct_matrix_force_cnsts construct_supercell convert_forces_from_Rybohr_to_eVang equilibrium_frac fourier_interpolation generate_kgrid generate_supercell_kpoint_mesh_qe generate_supercells lte lte_lower)

# compile programs
for program in ${programs[@]}; do
	gfortran $src$program -J$mod -o $bin$program -lblas -llapack
done

# --------------------
# quadratic
# --------------------
# change src directory
src=src/quadratic/programs/

# copy shell scripts
cp src/quadratic/scripts/*.sh bin/

# list programs
programs=(band_folding calculate_anharmonic calculate_bs calculate_gap generate_amplitudes generate_quadratic_configurations generate_sc_path quadratic_spline vscf_1d)

# compile programs
for program in ${programs[@]}; do
	gfortran $src$program -J$mod -o $bin$program -lblas -llapack
done
