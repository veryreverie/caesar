#!/bin/bash

# Script to transform structure file to qe in file
# $1=structure.dat, $2=pseudo.in, $3=kpoints.in, $4=structure.in

# Read in 'structure.dat'
lattice_line=$(awk -v IGNORECASE=1 '$1~/Lattice/{print NR}' $1) 
atoms_line=$(awk -v IGNORECASE=1 '$1~/Atoms/{print NR}' $1)
symmetry_line=$(awk -v IGNORECASE=1 '$1~/Symmetry/{print NR}' $1)
end_line=$(awk -v IGNORECASE=1 '$1~/End/{print NR}' $1)

no_atoms=$(( $symmetry_line-1-$atoms_line ))

echo 'nat='$no_atoms >> $4
echo '/&end' >> $4
cat $2 >> $4
echo 'CELL_PARAMETERS bohr' >> $4
awk "NR==$(($lattice_line + 1)),NR==$(($atoms_line - 1)) " \
  '{print $1 " " $2 " " $3}' $1 >> $4
echo 'ATOMIC_POSITIONS bohr' >> $4
awk "NR==$(($atoms_line + 1)),NR==$(($symmetry_line - 1)) " \
  '{print $1 " " $3 " " $4 " " $5}' $1 >> $4
cat $3 >> $4
