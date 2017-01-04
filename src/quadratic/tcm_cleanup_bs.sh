#!/bin/bash

# Script to clean-up quadratic at Rutgers


no_sc=$(ls -1d Supercell_* | wc -l)
seedname=$( awk '{print $1}' seedname.txt )
sampling_amplitude=$( awk 'NR==1 {print $1}' mapping.dat)
sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)
no_sampling_points=$(( $sampling_point_final-$sampling_point_init ))

# Loop over supercells
kpoint_counter=1
for i in `seq 1 $no_sc`;
do
  sdir=Supercell_$i
  if [ ! -f "$sdir/acoustic.dat" ]; then

    no_atoms=$( awk 'NR==1 {print $1}' $sdir/equilibrium.dat )
    no_modes=$(( $no_atoms*3 ))

    # Static calculation
    static_dir=$sdir/static
    caesar eigenval_castep_to_bands $static_dir/$seedname.bands $static_dir
    rm $static_dir/*.orbitals

    while read fline ; do
      line=($fline)
      big_point=${line[0]}
      kdir=$sdir/kpoint/$big_point/configurations
      cd $kdir
      
      for j in `seq 1 $no_modes`; do
        for k in `seq $sampling_point_init $sampling_point_final`; do
          mdir=$kdir/mode.$j.$k
          if [ -d "$mdir" ];then
            if [ -e "$mdir/$seedname.castep" ]; then
              caesar eigenval_castep_to_bands $mdir/$seedname.bands $mdir
              rm $mdir/*.orbitals
            fi
          fi
        done # loop over sampling points per mode
      done # loop over modes
    done < $sdir/list.dat
  fi
done  # Loop over supercells
