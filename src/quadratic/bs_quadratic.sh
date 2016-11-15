#!/bin/bash

# Script to calculate quadratic band gap correction

#######################################
####        MAIN   PROGRAM         ####
#######################################

# Define some useful variables
echo "What is the k-point of interest?"
read kpoint
echo "What is the band number of interest (for the primitive cell)?"
echo "(In case of band degeneracy, the highest band is required)"
read band
#echo "Do you want the VBM (1) or the CBM (2)?"
#read band_type
echo "What is the band degeneracy?"
read band_degeneracy
#echo "What is the equilibrium band value?"
#read band_energy

no_sc=$(ls -1d Supercell_* | wc -l)

sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)

mkdir bs
cp mapping.dat bs

# Obtain relevant band for each supercell
cd Supercell_1
  cd static
    band_energy=$( awk -v awk_band=$band 'NR==awk_band {print}' kpoint.${kpoint}.dat )
  cd ../
cd ../

# Loop over static supercells
for i in `seq 1 $no_sc`;do

  cd Supercell_$i
    cd static
      no_bands=$( wc -l < kpoint.${kpoint}.dat  )
      echo $band_energy $kpoint $no_bands > input.dat
      caesar band_folding
      rm input.dat
    cd ../
  cd ../

done # loop over static supercells


# Loop over supercells
kpoint_counter=1
for i in `seq 1 $no_sc`;do

  cd Supercell_$i

    band_ref=$( awk '{print}' static/band_number.dat )

    # Collect quadratic results
    no_kpoints=$(ls -1d kpoint.* | wc -l)
    no_atoms=$(awk 'NR==1 {print}' equilibrium.dat)
    no_modes=$(( $no_atoms*3 ))
    no_atoms_sc=$(awk 'NR==1 {print}' super_equilibrium.dat)
    sc_size=$(( $(( $no_atoms_sc/$no_atoms )) | bc ))

    # Loop over k-points
    for j in `seq $kpoint_counter $(( $kpoint_counter+($no_kpoints-1) ))`; do
      cd kpoint.${j}
      cp frequency.*.dat ../../bs/
      cd configurations

        for k in `seq 1 $no_modes`; do
          for l in `seq $sampling_point_init $sampling_point_final`; do
            if [ -e "mode.${k}.${l}" ] && [ "$l" -ne "0" ]; then
              cd mode.${k}.${l}
                band_first_line=$(( $band_ref-($band_degeneracy-1) )) 
                band_last_line=$(( $band_ref ))
                awk -v awk_first_line=$band_first_line -v awk_last_line=$band_last_line 'NR==awk_first_line,NR==awk_last_line' kpoint.$kpoint.dat > bs.${j}.${k}.${l}.dat
                mv bs.${j}.${k}.${l}.dat ../../../../bs/
              cd ../
            fi
          done
        done

      cd ../../
    done # Loop over k-points
    kpoint_counter=$(( $kpoint_counter+$no_kpoints ))
  cd ../

done

# Calculatee renormalised band
cp ibz.dat bs
cd bs
  no_kpoints=$( wc -l < ibz.dat  )
  echo $no_kpoints $no_modes $band_degeneracy > input.dat
  for i in `seq 1 $no_kpoints`; do
    for j in `seq 1 $no_modes`; do
      if [ -e "bs.${i}.${j}.${sampling_point_init}.dat" ];then
        cat bs.${i}.${j}.${sampling_point_init}.dat bs.${i}.${j}.${sampling_point_final}.dat > bs.${i}.${j}.dat            
        echo $band_energy >> bs.${i}.${j}.dat
        rm bs.${i}.${j}.${sampling_point_init}.dat bs.${i}.${j}.${sampling_point_final}.dat
      fi
    done 
  done 
  caesar calculate_bs
  rm input.dat
cd ../

