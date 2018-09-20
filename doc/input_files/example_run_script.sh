#!/bin/bash

file_type=$1
directory=$2
no_cores=$3
no_nodes=$4
seedname=$5

# There are as many additional arguments as are passed to run_script_data.
run_script_data_1=$6
run_script_data_2=$7
# ...

cores_per_node=$((no_cores / no_nodes))

if [ "$file_type" = "castep" ]; then
  cp $seedname.param $directory
fi

cd $directory

if [ "$file_type" = "castep" ]; then
  mpirun -n $cores_per_node castep.mpi $seedname
elif [ "$file_type" = "vasp" ]; then
  echo "Error! vasp run script not yet written."
  exit
elif [ "$file_type" = "qe" ]; then
  mpirun -n $cores_per_node pw.x -i $seedname.in > $seedname.out
fi

cd -
