#!/bin/bash

# Script to transform structure file to qe in file

#######################################
####        MAIN   PROGRAM         ####
#######################################

# Read in 'structure.dat'
lattice_line=$(awk -v IGNORECASE=1 '$1~/Lattice/{print NR}' structure.dat) 
atoms_line=$(awk -v IGNORECASE=1 '$1~/Atoms/{print NR}' structure.dat)
symmetry_line=$(awk -v IGNORECASE=1 '$1~/Symmetry/{print NR}' structure.dat)
end_line=$(awk -v IGNORECASE=1 '$1~/End/{print NR}' structure.dat)

no_atoms=$(( $symmetry_line-1-$atoms_line ))

echo 'nat='$no_atoms > structure.in
cat top.in structure.in > structure_temp.in
mv structure_temp.in structure.in
echo '/&end' >> structure.in
cat structure.in pseudo.in > structure_temp.in
mv structure_temp.in structure.in
echo 'CELL_PARAMETERS bohr' >> structure.in
awk -v awk_lattice_line=$lattice_line -v awk_atoms_line=$atoms_line 'NR==(awk_lattice_line+1), NR==(awk_atoms_line-1) {print $1 " " $2 " " $3} ' structure.dat >> structure.in
echo 'ATOMIC_POSITIONS bohr' >> structure.in
awk -v awk_atoms_line=$atoms_line -v awk_symmetry_line=$symmetry_line 'NR==(awk_atoms_line+1), NR==(awk_symmetry_line-1) {print $1 " " $3 " " $4 " " $5} ' structure.dat >> structure.in
cat structure.in kpoints.in >> structure_temp.in
mv structure_temp.in structure.in
