#!/bin/bash

# Script to run quadratic in TCM cluster

harmonic_path=$( awk '{print $1}' harmonic_path.dat)

sampling_amplitude=$( awk 'NR==1 {print $1}' mapping.dat)
sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)
no_sampling_points=$(( $sampling_point_final-$sampling_point_init ))

atoms_line=$(awk -v IGNORECASE=1 '$1~/Atoms/{print NR}' \
   $harmonic_path/structure.dat)
symmetry_line=$(awk -v IGNORECASE=1 '$1~/Lattice/{print NR}' \
   $harmonic_path/structure.dat)
no_atoms=$($symmetry_line-$atoms_line-1)
no_modes=$(( $no_atoms*3 ))

echo "What is the first supercell to run?"
read first_sc
echo "What is the last supercell to run?"
read last_sc
echo "How many cores per run?"
read num_cores

# Read code and seedname
code=$( awk '{print}' code.txt )
seedname=$( awk '{print}' seedname.txt )
seedname_nscf=$seedname.nscf

# Loop over supercells
for i in `seq $first_sc $last_sc`; do
  sdir=Supercell_$i
  no_atoms_sc=$(awk 'NR==1 {print $1}' \
     $harmonic_path/$sdir/super_equilibrium.dat)

  # Run static first
  static_dir=$sdir/static
  if [ "$code" = "castep" ]; then
    caesar rundft $code $static_dir $num_cores
  elif [ "$code" = "qe" ]; then
    caesar rundft $code $static_dir $num_cores $seedname
    if [ -e "$static_dir/$seedname_nscf.in" ]; then
      caesar rundft $code $static_dir $num_cores $seedname_nscf
    fi
  fi
    
  # loop over kpoints
  while read fline ; do
    line=($fline)
    big_point=${line[0]}
    kdir=$sdir/kpoint.$big_point/configurations
    for j in `seq 1 $no_modes`; do
      for k in `seq $sampling_point_init $sampling_point_final`; do
        mdir=$kdir/mode.$j.$k
        if [ -e "$mdir/structure.dat" ]; then
          if [ "$code" = "castep" ]; then
            caesar rundft $code $mdir $num_cores
          elif [ "$code" = "qe" ]; then
            caesar rundft $code $mdir $num_cores $seedname
            if [ -e "$mdir/$seedname_nscf.in" ]; then
              caesar rundft $code $mdir $num_cores $seedname_nscf
            fi
          fi
        fi
      done # loop over sampling points per mode
    done # loop over modes
  done < $harmonic_path/$sdir/list.dat
done
echo "Done."
