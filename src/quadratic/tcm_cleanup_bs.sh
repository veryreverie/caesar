#!/bin/bash

# Script to clean-up quadratic at Rutgers


no_sc=$(ls -1d Supercell_* | wc -l)
seedname=$( awk '{print $1}' Supercell_1/seedname.txt )
sampling_amplitude=$( awk 'NR==1 {print $1}' mapping.dat)
sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)
no_sampling_points=$(( $sampling_point_final-$sampling_point_init ))

# Loop over supercells
kpoint_counter=1
for i in `seq 1 $no_sc`;
do

  cd Supercell_$i
  
  if [ ! -f "acoustic.dat" ]; then

  no_atoms=$( awk 'NR==1 {print $1}' equilibrium.dat )
  no_modes=$(( $no_atoms*3 ))

  # Static calculation
  cd static
    cp ../seedname.txt .
    caesar eigenval_castep_to_bands
    rm seedname.txt *.orbitals
  cd ../

    while read line ; do

      big_point=$(echo ${line} | awk '{print $1}')
      cd kpoint.$big_point/configurations

      for j in `seq 1 $no_modes`; do
        for k in `seq $sampling_point_init $sampling_point_final`; do
          if [ -d "mode.${j}.${k}" ];then
            cd mode.${j}.${k}
              if [ -e "$seedname.castep" ]; then
                cp ../../../seedname.txt .
                caesar eigenval_castep_to_bands
                rm seedname.txt *.orbitals
              fi
            cd ../
          fi
        done # loop over sampling points per mode
      done # loop over modes

      cd ../../  # cd kpoint.$big_point/configurations

    done < list.dat


  fi
  cd ../  # cd Supercell_$i

done  # Loop over supercells
