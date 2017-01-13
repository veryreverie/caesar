#!/bin/bash

# Script to run hamonic in TCM cluster

echo "What code do you want to use (castep,vasp,qe)?"
read code

if [ "$code" = "castep" ];then

  if [ ! -d "castep" ];then
    echo "Error! The directory 'castep' does not exist." 
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
        caesar rundft castep $dir $num_cores
        caesar fetch_forces_castep   \
               $dir/$seedname.castep \
               $dir/structure.dat    \
               $atom                 \
               $disp                 \
               $dir/forces.dat
      done
  
    done < $sdir/force_constants.dat 

  done
  echo "Done."

elif [ "$code" = "vasp" ]; then

  if [ ! -d "vasp" ];then
    echo "Error! The directory 'vasp' does not exist." 
    exit 1
  fi 

elif [ "$code" = "qe" ]; then

  if [ ! -d "qe" ];then
    echo "Error! The directory 'qe' does not exist." 
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
    while read fline; do
      line=($fline)
      atom=${line[0]}
      disp=${line[1]}
      
      paths=(positive negative)
      for path in ${paths[@]}; do
        dir=$sdir/atom.$atom.disp.$disp/path
        caesar rundft qe $dir $num_cores $seedname
        caesar fetch_forces_qe      \
               $dir/$seedname.out   \
               $atom                \
               $disp                \
               $dir/forces.out
      done

    done < $sdir/force_constants.dat

  done
  echo "Done."

else 

  echo "Error! This code is not supported."

fi
