#!/bin/bash

# Script to clean-up quadratic at Rutgers

harmonic_path=$( awk '{print $1}' harmonic_path.dat)

no_sc=$( awk '{print $1}' $harmonic_path/no_sc.dat)
seedname=$( awk '{print $1}' seedname.txt )
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

# Loop over supercells
for i in `seq 1 $no_sc`;
do
  sdir=Supercell_$i
  echo "Supercell" $i

  if [ ! -f "$sdir/acoustic.dat" ]; then

    # Static calculation
    f=$sdir/static/$seedname.castep
    line_number=$(awk '/Final energy/{x=NR}END{print x}' $f)
    energy=$(awk "NR==$line_number {print \$5}" $f)
    echo $energy > $sdir/static/energy.dat

    while read fline ; do
      line=($fline)
      big_point=${line[0]}
      kdir=kpoint.$big_point/configurations

      for j in `seq 1 $no_modes`; do
        for k in `seq $sampling_point_init $sampling_point_final`; do
          mdir=$kdir/mode.$j.$k
          if [ -d "$mdir" ];then 
            f=$mdir/seedname.castep
            if [ -e "$f" ]; then
              line_number=$(awk '/Final energy/{x=NR}END{print x}' $f)
              energy=$(awk "NR==$line_number {print \$5}" $f)
              echo $energy > $mdir/energy.dat
            fi
          fi
        done # loop over sampling points per mode
      done # loop over modes
    done < $harmonic_path/$sdir/list.dat
  fi
done  # Loop over supercells
