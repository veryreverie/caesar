#!/bin/bash
#!
#! Example run script, to call DFT from inside castep.

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

# Copy across Castep .param file or Quantum Espresso .USP files.
# This copies from the directory where Caesar is called (or that specified by
#    the -d flag) to the working directory where the electronic structure
#    calculation will be run.
if [ "$file_type" = "castep" ]; then
  cp $seedname.param $directory
elif [ "$file_type" = "quantum_espresso" ]; then
  for i in *.UPF; do
    cp "$i" "$directory/$i" || break
  done
fi

# cd into the directory where the electronic structure calculation will be run.
cd $directory

# Run the electronic structure calculation.
if [ "$file_type" = "castep" ]; then
  mpirun -ppn $cores_per_node -np $no_cores castep.mpi $seedname
elif [ "$file_type" = "quantum_espresso" ]; then
  mpirun -ppn $cores_per_node -np $no_cores pw.x -i $seedname.in > $seedname.out
fi

# Return to the starting directory.
cd -
