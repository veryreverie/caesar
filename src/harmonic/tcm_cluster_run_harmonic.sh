#!/bin/bash

# Script to run hamonic in TCM cluster

echo "What code do you want to use (castep,vasp,qe)?"
read code

if [ "$code" = "castep" ]; then
elif [ "$code" = "vasp" ]; then
  echo "Error! vasp is not currently supported."
  exit 1
elif [ "$code" = "qe"]; then
else
  echo "Error! The code $code is not supported."
  echo "Please choose one of: castep vasp qe"
  exit 1
fi

if [ ! -d "$code" ];then
  echo "Error! The directory '$code' does not exist." 
  exit 1
fi 

echo "What is the first supercell to run?"
read first_sc
echo "What is the last supercell to run?"
read last_sc
echo "How many cores per run?"
read num_cores

# Get seedname
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
      caesar rundft $code $dir $num_cores
      caesar fetch_forces          \
             $code                 \
             $dir/$seedname.castep \
             $atom                 \
             $disp                 \
             $dir/forces.dat
    done

  done < $sdir/force_constants.dat 

done
echo "Done."
