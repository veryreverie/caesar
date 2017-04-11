#!/bin/bash

dft_code=$1
directory=$2
num_cores=$3
seedname=$4

if [ "$dft_code" = "castep" ]; then
  cp $seedname.param $directory
fi

cd $directory

if [ "$dft_code" = "castep" ]; then
  mpirun -n $num_cores castep.mpi $seedname
elif [ "$dft_code" = "vasp" ]; then
  echo "Error! vasp run script not yet written."
  exit
elif [ "$dft_code" = "qe" ]; then
  mpirun -n $num_cores pw.x -i $seedname.in > $seedname.out
fi

cd -
