#!/bin/bash

file_type=$1
directory=$2
num_cores=$3
seedname=$4

if [ "$file_type" = "castep" ]; then
  cp $seedname.param $directory
fi

cd $directory

if [ "$file_type" = "castep" ]; then
  mpirun -n $num_cores castep.mpi $seedname
elif [ "$file_type" = "vasp" ]; then
  echo "Error! vasp run script not yet written."
  exit
elif [ "$file_type" = "qe" ]; then
  mpirun -n $num_cores pw.x -i $seedname.in > $seedname.out
fi

cd -
