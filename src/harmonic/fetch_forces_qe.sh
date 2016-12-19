#!/bin/bash

# Script to fetch forces from qe output 
 
seedname=$(awk '{print}' seedname.txt)
forces_line=$(awk -v IGNORECASE=1 '/Forces acting/{print NR}' $seedname.out)
atoms_line=$(awk -v IGNORECASE=1 '/Atoms/{print NR}' structure.dat)
symmetry_line=$(awk -v IGNORECASE=1 '/Symmetry/{print NR}' structure.dat)
no_atoms=$(( $(( $symmetry_line-($atoms_line+1))) | bc ))
atom=$(awk '{print $1}' displacement.dat)
disp=$(awk '{print $2}' displacement.dat)

for (( i=1; i<=$no_atoms; i++ )) do
  forces=($(awk "NR==$(($forces_line + 1 + $i)) " \
    '{print $7 " " $8 " " $9}' $seedname.out))
  j=1
  for force in ${forces[@]}; do
    force_eV_ang=$(echo "scale=20; $force*13.605698066/0.529117" | bc)
    echo $atom $disp $i $j $force_eV_ang >> forces.out
    j=$(($j + 1))
  done
done 
