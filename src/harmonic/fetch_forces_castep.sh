#!/bin/bash

# Script to fetch forces from castep output 
 
seedname=$(awk '{print}' seedname.txt)
forces_line=$(awk -v IGNORECASE=1 '/Forces/{print NR}' $seedname.castep)
atoms_line=$(awk -v IGNORECASE=1 '/Atoms/{print NR}' structure.dat)
symmetry_line=$(awk -v IGNORECASE=1 '/Symmetry/{print NR}' structure.dat)
no_atoms=$(( $(( $symmetry_line-($atoms_line+1))) | bc ))
atom=$(awk '{print $1}' displacement.dat)
disp=$(awk '{print $2}' displacement.dat)

awk "NR==$(($forces_line + 6)),NR==$(($forces_line + 6 + $no_atoms)) " \
  '{print $4 " " $5 " " $6}' $seedname.castep > forces.dat

for (( i=1; i<=$no_atoms; i++ )) do

  for (( j=1; j<=3; j++ )) do
    awk "NR==$i {print \$$j}" forces.dat | \
      awk "{print $atom ' ' $disp ' ' $i ' ' $j ' ' \$0}" >> forces_global.dat
  done

done 

mv forces_global.dat forces.dat
