#!/bin/bash

# Script to transform vasp EIGENVAL file to bands file

#######################################
####        MAIN   PROGRAM         ####
#######################################

# Read in EIGENVAL
no_kpoints=$(awk 'NR==6 {print $2}' EIGENVAL)
no_bands=$(awk 'NR==6 {print $3}' EIGENVAL)

for i in `seq 1 $no_kpoints`; do
  first_line=$(( 6+($no_bands+2)*($i-1)+3 ))
  last_line=$(( 6+($no_bands+2)*($i-1)+2+$no_bands ))
  for j in `seq 1 $no_bands`; do
    awk -v awk_first_line=$first_line -v awk_last_line=$last_line 'NR==awk_first_line,NR==awk_last_line {print $2}' EIGENVAL > kpoint.${i}.dat
  done
done

