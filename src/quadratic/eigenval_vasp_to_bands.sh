#!/bin/bash

# Script to transform vasp EIGENVAL file to bands file
# $1=input EIGENVAL file, $2=output directory

# Read in EIGENVAL
no_kpoints=$(awk 'NR==6 {print $2}' $1)
no_bands=$(awk 'NR==6 {print $3}' $1)

for i in `seq 1 $no_kpoints`; do
  first_line=$(( 6+($no_bands+2)*($i-1)+3 ))
  last_line=$(( 6+($no_bands+2)*($i-1)+2+$no_bands ))
  for j in `seq 1 $no_bands`; do
    awk "NR==$first_line,NR==$last_line {print \$2}" $1 > $2/kpoint.${i}.dat
  done
done

