#!/bin/bash

# Script to run quadratic in TCM cluster

# Read in user inputs
echo "What is the first supercell to run?"
read first_sc
echo "What is the last supercell to run?"
read last_sc
echo "How many cores per run?"
read num_cores

# Read in previous user inputs
dft_code=$( awk 'NR==1 {print}' user_input.txt )
seedname=$( awk 'NR==2 {print}' user_input.txt )
seedname_nscf=$seedname.nscf
harmonic_path=$( awk 'NR==3 {print}' user_input.txt)

# Read in mapping file
sampling_amplitude=$( awk 'NR==1 {print $1}' mapping.dat)
sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)
no_sampling_points=$(( $sampling_point_final-$sampling_point_init ))

# Read in structure file
atoms_line=$(awk -v IGNORECASE=1 '$1~/Atoms/{print NR}' \
   $harmonic_path/structure.dat)
symmetry_line=$(awk -v IGNORECASE=1 '$1~/Lattice/{print NR}' \
   $harmonic_path/structure.dat)
no_atoms=$($symmetry_line-$atoms_line-1)
no_modes=$(( $no_atoms*3 ))

# Loop over supercells
for sc_id in `seq $first_sc $last_sc`; do
  sdir=Supercell_$sc_id

  # Run static dft
  static_dir=$sdir/static
  if [ "$dft_code" = "castep" ]; then
    cp $dft_code/$seedname.param $static_dir
    caesar rundft $dft_code $static_dir $num_cores $seedname
  elif [ "$dft_code" = "qe" ]; then
    caesar rundft $dft_code $static_dir $num_cores $seedname
    if [ -e "$static_dir/$seedname_nscf.in" ]; then
      caesar rundft $dft_code $static_dir $num_cores $seedname_nscf
    fi
  fi
done

# loop over kpoints
while read fline ; do
  line=($fline)
  kpoint=${line[0]}
  sc_id=${line[2]}
  kdir=kpoint.$kpoint
  
  # Skip kpoints in other supercells
  if ["$sc_id" -ge "$first_sc"] && ["$sc_id" -le "$last_sc"]; then
    # Loop over modes
    for j in `seq 1 $no_modes`; do
      # Loop over sampling points
      for k in `seq $sampling_point_init $sampling_point_final`; do
        mdir=$kdir/mode.$j.$k
        # Check that dft needs to be run
        if [ -e "$mdir/structure.dat" ]; then
          # Run dft
          if [ "$dft_code" = "castep" ]; then
            cp $dft_code/$seedname.param $mdir
            caesar rundft $dft_code $mdir $num_cores $seedname
          elif [ "$dft_code" = "qe" ]; then
            caesar rundft $dft_code $mdir $num_cores $seedname
            if [ -e "$mdir/$seedname_nscf.in" ]; then
              caesar rundft $dft_code $mdir $num_cores $seedname_nscf
            fi
          fi
        fi
      done
    done
  fi
done < $harmonic_path/list.dat
echo "Done."
