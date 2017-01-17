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
      dir=$sdir/atom.$atom.disp.$disp/$path
      if [ "$code" = "castep" ]; then
        caesar rundft $code $dir $num_cores
        caesar fetch_forces          \
               $code                 \
               $dir/$seedname.castep \
               $atom                 \
               $disp                 \
               $dir/forces.dat
      elif [ "$code" = "qe" ]; then
        caesar rundft $code $dir $num_cores $seedname
        caesar fetch_forces       \
               $code              \
               $dir/$seedname.out \
               $atom              \
               $disp              \
               $dir/forces.dat
      fi
    done

  done < $sdir/force_constants.dat 

done
echo "Done."
