#!/bin/bash

# Script to calculate anharmonic 1-dimensional correction

#######################################
####        MAIN   PROGRAM         ####
#######################################


no_sc=$(ls -1d Supercell_* | wc -l)
seedname=$( awk '{print $1}' Supercell_1/seedname.txt )
sampling_amplitude=$( awk 'NR==1 {print $1}' mapping.dat)
sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)
no_sampling_points=$(( $sampling_point_final-$sampling_point_init ))


mkdir anharmonic
cp mapping.dat anharmonic

# Loop over supercells
for i in `seq 1 $no_sc`;do

  cd Supercell_$i

    if [ "$i" -eq 1 ]; then
      cp list.dat ../anharmonic
    else
      cat ../anharmonic/list.dat list.dat > tomeu.dat
      mv tomeu.dat ../anharmonic/list.dat
    fi

    echo "Working with Supercell" $i
    if [ ! -f "acoustic.dat" ]; then

    no_atoms=$( awk 'NR==1 {print $1}' equilibrium.dat )
    no_atoms_sc=$(awk 'NR==1 {print}' super_equilibrium.dat)
    no_modes=$(( $no_atoms*3 ))
    no_cells=$(( no_atoms_sc/no_atoms ))


    while read line ; do

      big_point=$(echo ${line} | awk '{print $1}')
      echo $no_cells > size.${big_point}.dat
      mv size.${big_point}.dat ../anharmonic
      cd kpoint.$big_point/configurations

      for j in `seq 1 $no_modes`; do
        if [ -d "mode.${j}.$sampling_point_init" ];then
          for k in `seq $sampling_point_init $sampling_point_final`; do
            if [ -d "mode.${j}.${k}" ];then
              cd mode.${j}.${k}
                if [ -e "$seedname.castep" ]; then
                  if [ "$k" -eq "$sampling_point_init" ]; then
                    cp energy.dat ../
                  else
                    cat ../energy.dat energy.dat > tomeu.dat
                    mv tomeu.dat ../energy.dat
                  fi
                else
                  cat ../energy.dat ../../../static/energy.dat > tomeu.dat
                  mv tomeu.dat ../energy.dat
                fi
              cd ../
            elif [ "$k" -eq "0" ]; then
              cat energy.dat ../../static/energy.dat > tomeu.dat  
              mv tomeu.dat energy.dat
            fi
          done # loop over sampling points per mode
          mv energy.dat energy.${big_point}.${j}.dat
        fi
      done # loop over modes
    
      mv energy.*.dat ../../../anharmonic/

      cd ../../  # cd kpoint.$big_point/configurations
      cp kpoint.${big_point}/frequency.*.dat ../anharmonic/

    done < list.dat

    else 
      cp acoustic.dat ../anharmonic
    fi

  cd ../  # cd Supercell_$i

done

# Calculate anharmonic 1-dimensional correction
cp ibz.dat anharmonic
cd anharmonic

integration_points=5000

while read line ; do

  big_point=$(echo ${line} | awk '{print $1}')
  echo $big_point > kpoint.dat

  if [ "$big_point" -eq 1 ] && [ -f "acoustic.dat" ]; then
 
    acoustic=1

  else

  for j in `seq 1 $no_modes`; do

    # Generate amplitudes
    echo $j > mode.dat
    if [ -e "energy.${big_point}.${j}.dat" ];then
      cp energy.${big_point}.${j}.dat working_energy.dat
      cp frequency.${big_point}.${j}.dat working_frequency.dat
      cp size.${big_point}.dat working_size.dat
      caesar generate_amplitudes
      mv amplitude_energy.dat energy.${big_point}.${j}.dat 
      rm working_energy.dat working_frequency.dat working_size.dat

      # Fit splines
      echo $no_sampling_points $integration_points > fit_input.dat
      cp energy.${big_point}.${j}.dat fit_energy.dat
      caesar quadratic_spline
      mv indep_pot.dat interp_energy.${big_point}.${j}.dat
      rm fit_input.dat fit_energy.dat 

      # Calculate 1-d anharmonic energy 
      cp interp_energy.${big_point}.${j}.dat interp_energy.dat
      cp frequency.${big_point}.${j}.dat frequency.dat
      max_amplitude=$( awk 'NR==1 {print $1}' interp_energy.dat)
      #max_amplitude=$(( $max_amplitude*(-1) ))
      echo $max_amplitude > max_amplitude.dat
      echo $integration_points > integration_points.dat
      caesar vscf_1d
      mv eigenvals.dat eigenvals.${big_point}.${j}.dat
      mv eigenvecs.dat eigenvecs.${big_point}.${j}.dat
      mv anh_pot.dat anh_pot.${big_point}.${j}.dat
      rm interp_energy.dat frequency.dat max_amplitude.dat
    fi
    
  done # loop over modes

  fi

done < list.dat


no_kpoints=$( wc -l < ibz.dat  )
echo $no_kpoints $no_modes > input.dat

echo $no_kpoints $no_modes > input.dat
caesar calculate_anharmonic
rm input.dat




