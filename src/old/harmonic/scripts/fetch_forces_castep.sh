#!/bin/bash

# Script to fetch forces from castep output 
 
seedname=$(awk '{print}' seedname.txt)
forces_line=$(awk -v IGNORECASE=1 '/Forces/{print NR}' $seedname.castep)
atoms_line=$(awk -v IGNORECASE=1 '/Atoms/{print NR}' structure.dat)
symmetry_line=$(awk -v IGNORECASE=1 '/Symmetry/{print NR}' structure.dat)
no_atoms=$(( $(( $symmetry_line-($atoms_line+1))) | bc ))
atom=$(awk '{print $1}' displacement.dat)
disp=$(awk '{print $2}' displacement.dat)

awk -v awk_forces_line=$forces_line -v awk_no_atoms=$no_atoms 'NR==awk_forces_line+6,NR==awk_forces_line+6+awk_no_atoms {print $4 " " $5 " " $6}' $seedname.castep > forces.dat

for (( i=1; i<=$no_atoms; i++ )) do

  awk -v awk_line=$i 'NR==awk_line {print $1}' forces.dat > forces_temp.dat
  awk -v awk_atom=$atom -v awk_disp=$disp -v awk_atom2=$i '{print awk_atom " " awk_disp " " awk_atom2 " " 1 " " $0}' forces_temp.dat > forces_temp2.dat
  cat forces_temp2.dat >> forces_global.dat
  awk -v awk_line=$i 'NR==awk_line {print $2}' forces.dat > forces_temp.dat
  awk -v awk_atom=$atom -v awk_disp=$disp -v awk_atom2=$i '{print awk_atom " " awk_disp " " awk_atom2 " " 2 " " $0}' forces_temp.dat > forces_temp2.dat
  cat forces_temp2.dat >> forces_global.dat
  awk -v awk_line=$i 'NR==awk_line {print $3}' forces.dat > forces_temp.dat
  awk -v awk_atom=$atom -v awk_disp=$disp -v awk_atom2=$i '{print awk_atom " " awk_disp " " awk_atom2 " " 3 " " $0}' forces_temp.dat > forces_temp2.dat
  cat forces_temp2.dat >> forces_global.dat

done 

rm forces_temp.dat forces_temp2.dat
mv forces_global.dat forces.dat

