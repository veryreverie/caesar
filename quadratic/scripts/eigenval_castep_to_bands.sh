#!/bin/bash

# Script to transform vasp EIGENVAL file to bands file

#######################################
####        MAIN   PROGRAM         ####
#######################################

seedname=$(awk '{print $1}' seedname.txt)

# Read in EIGENVAL
no_kpoints=$(awk 'NR==1 {print $4}' $seedname.bands)
no_bands=$(awk 'NR==4 {print $4}' $seedname.bands)

for i in `seq 1 $no_kpoints`; do
  first_line=$(( 9+($no_bands+2)*($i-1)+3 ))
  last_line=$(( 9+($no_bands+2)*($i-1)+2+$no_bands ))
  for j in `seq 1 $no_bands`; do
    awk -v awk_first_line=$first_line -v awk_last_line=$last_line 'NR==awk_first_line,NR==awk_last_line {print $1}' $seedname.bands > kpoint.${i}.dat
  done
done

