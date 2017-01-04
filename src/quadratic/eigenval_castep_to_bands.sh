#!/bin/bash

# Script to transform castep bands file to bands file
# $1=input bands file, $2=output directory

# Read in bands file
no_kpoints=$(awk 'NR==1 {print $4}' $1)
no_bands=$(awk 'NR==4 {print $4}' $1)

for i in `seq 1 $no_kpoints`; do
  first_line=$(( 9+($no_bands+2)*($i-1)+3 ))
  last_line=$(( 9+($no_bands+2)*($i-1)+2+$no_bands ))
  for j in `seq 1 $no_bands`; do
    awk "NR==$first_line,NR==$last_line {print \$1}" $1 > $2/kpoint.$i.dat
  done
done

