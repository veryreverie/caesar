#!/bin/bash

# Script to fetch forces from castep output 
# $1=seedname.castep, $2=structure.dat, $3=atom, $4=disp, $5=output file

castep_file=$1
structure_file=$2
atom=$3
disp=$4
ofile=$5

tempfile=fetch_forces_castep_tempfile

forces_line=$(awk -v IGNORECASE=1 '/Forces/{print NR}' $castep_file)
atoms_line=$(awk -v IGNORECASE=1 '/Atoms/{print NR}' $structure_file)
symmetry_line=$(awk -v IGNORECASE=1 '/Symmetry/{print NR}' $structure_file)
no_atoms=$(( $(( $symmetry_line-($atoms_line+1))) | bc ))

awk "NR==$(($forces_line + 6)),NR==$(($forces_line + 6 + $no_atoms)) " \
  '{print $4 " " $5 " " $6}' $castep_file > $tempfile

for (( i=1; i<=$no_atoms; i++ )) do
  for (( j=1; j<=3; j++ )) do
    awk "NR==$i {print \$$j}" $tempfile | \
      awk "{print $atom ' ' $disp ' ' $i ' ' $j ' ' \$0}" >> $ofile
  done
done 

rm $tempfile
