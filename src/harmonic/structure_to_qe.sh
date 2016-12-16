#!/bin/bash

# Script to transform structure file to qe in file

#######################################
####        MAIN   PROGRAM         ####
#######################################

# Read in 'structure.dat'
lattice_line=$(awk -v IGNORECASE=1 '$1~/Lattice/{print NR}' $1/structure.dat) 
atoms_line=$(awk -v IGNORECASE=1 '$1~/Atoms/{print NR}' $1/structure.dat)
symmetry_line=$(awk -v IGNORECASE=1 '$1~/Symmetry/{print NR}' $1/structure.dat)
end_line=$(awk -v IGNORECASE=1 '$1~/End/{print NR}' $1/structure.dat)

no_atoms=$(( $symmetry_line-1-$atoms_line ))

cat $1/top.in > $1/structure.in
echo 'nat='$no_atoms >> $1/structure.in
echo '/&end' >> $1/structure.in
cat $1/pseudo.in >> $1/structure.in
echo 'CELL_PARAMETERS bohr' >> $1/structure.in
awk "NR==$(($lattice_line + 1)),NR==$(($atoms_line - 1)) " \
  '{print $1 " " $2 " " $3}' $1/structure.dat >> $1/structure.in
echo 'ATOMIC_POSITIONS bohr' >> $1/structure.in
awk "NR==$(($atoms_line + 1)),NR==$(($symmetry_line - 1)) " \
  '{print $1 " " $3 " " $4 " " $5}' $1/structure.dat $1/structure.in
cat $1/kpoints.in >> $1/structure.in
