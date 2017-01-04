#!/bin/bash

# Script to fetch forces from qe output 
# $1=seedname.out, $2=structure.dat, $3=atom, $4=disp, $5=forces.out
 
forces_line=$(awk -v IGNORECASE=1 '/Forces acting/{print NR}' $1)
atoms_line=$(awk -v IGNORECASE=1 '/Atoms/{print NR}' $2)
symmetry_line=$(awk -v IGNORECASE=1 '/Symmetry/{print NR}' $2)
no_atoms=$(( $(( $symmetry_line-($atoms_line+1))) | bc ))
atom=$3
disp=$4

for (( i=1; i<=$no_atoms; i++ )) do
  forces=($(awk "NR==$(($forces_line + 1 + $i)) " \
    '{print $7 " " $8 " " $9}' $1))
  j=1
  for force in ${forces[@]}; do
    force_eV_ang=$(echo "scale=20; $force*13.605698066/0.529117" | bc)
    echo $atom $disp $i $j $force_eV_ang >> $5
    j=$(($j + 1))
  done
done 
