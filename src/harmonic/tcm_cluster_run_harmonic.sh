#!/bin/bash

# Script to run hamonic in TCM cluster

echo "What is the first supercell to run?"
read first_sc
echo "What is the last supercell to run?"
read last_sc
echo "How many cores per run?"
read num_cores

# Read code and seedname
code=$( awk '{print}' code.txt )
seedname=$( awk '{print}' seedname.txt )

# Loop over supercells
for i in `seq $first_sc $last_sc`;
do
  sdir=Supercell_$i 

  # Loop over force constants
  while read fline ; do
    line=($fline)
    atom=${line[0]}
    disp=${line[1]}
    
    paths=(positive negative)
    for path in ${paths[@]}; do
      ddir=$sdir/atom.$atom.disp.$disp/$path
      if [ "$code" = "castep" ]; then
        cp $code/$seedname.param $ddir
      fi
      caesar rundft $code $ddir $num_cores $seedname
    done

  done < $sdir/force_constants.dat 

done
echo "Done."
