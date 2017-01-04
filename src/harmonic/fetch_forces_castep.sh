#!/bin/bash

# Script to fetch forces from castep output 
# $1=seedname.castep, $2=structure.dat, $3=atom, $4=disp, $5=output file
 
forces_line=$(awk -v IGNORECASE=1 '/Forces/{print NR}' $1)
atoms_line=$(awk -v IGNORECASE=1 '/Atoms/{print NR}' $2)
symmetry_line=$(awk -v IGNORECASE=1 '/Symmetry/{print NR}' $2)
no_atoms=$(( $(( $symmetry_line-($atoms_line+1))) | bc ))
atom=$3
disp=$4

awk "NR==$(($forces_line + 6)),NR==$(($forces_line + 6 + $no_atoms)) " \
  '{print $4 " " $5 " " $6}' $1 > $5

f=fetch_forces_castep_tempfile
for (( i=1; i<=$no_atoms; i++ )) do
  for (( j=1; j<=3; j++ )) do
    awk "NR==$i {print \$$j}" $5 | \
      awk "{print $atom ' ' $disp ' ' $i ' ' $j ' ' \$0}" >> $f
  done
done 

mv $f $5
