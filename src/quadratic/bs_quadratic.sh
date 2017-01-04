#!/bin/bash

# Script to calculate quadratic band gap correction

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

sampling_point_amplitude=$( awk 'NR==1 {print $1}' mapping.dat)
sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)

mkdir bs

# Obtain relevant band for each supercell
band_energy=$(awk "NR==$band {print}" Supercell_1/static/kpoint.$kpoint.dat)

# Loop over static supercells
for i in `seq 1 $no_sc`;do
  static_dir=Supercell_$i/static
  f=$static_dir/kpoint.$kpoint.dat
  no_bands=$( wc -l < $f  )
  caesar band_folding       \
         $f                 \
         $band_energy       \
         $no_bands          \
         $static_dir/band_number.dat
done # loop over static supercells


# Loop over supercells
kpoint_counter=1
for i in `seq 1 $no_sc`;do
  sdir=Supercell_$i
  
  band_ref=$( awk '{print}' $sdir/static/band_number.dat )
  
  # Collect quadratic results
  no_kpoints=$(ls -1d $sdir/kpoint.* | wc -l)
  no_atoms=$(awk 'NR==1 {print}' $sdir/equilibrium.dat)
  no_modes=$(( $no_atoms*3 ))
  no_atoms_sc=$(awk 'NR==1 {print}' $sdir/super_equilibrium.dat)
  sc_size=$(( $(( $no_atoms_sc/$no_atoms )) | bc ))
  
  # Loop over k-points
  for j in `seq $kpoint_counter $(( $kpoint_counter+($no_kpoints-1) ))`; do
    cp $sdir/kpoint.$j/frequency.*.dat bs
    
    kdir=$sdir/kpoint.$j/configurations
    for k in `seq 1 $no_modes`; do
      for l in `seq $sampling_point_init $sampling_point_final`; do
        mdir=$kdir/mode.$k.$l
        if [ -e "$mdir" ] && [ "$l" -ne "0" ]; then
          first_line=$(( $band_ref-($band_degeneracy-1) )) 
          last_line=$(( $band_ref ))
          awk "NR==$first_line,NR==$last_line" $mdir/kpoint.$kpoint.dat \
            > bs/bs.$j.$k.$l.dat
        fi
      done
    done
  done # Loop over k-points
  kpoint_counter=$(( $kpoint_counter+$no_kpoints ))

done

# Calculatee renormalised band
no_kpoints=$( wc -l < ibz.dat  )
for i in `seq 1 $no_kpoints`; do
  for j in `seq 1 $no_modes`; do
    if [ -e "bs/bs.$i.$j.$sampling_point_init.dat" ];then
      cat bs/bs.$i.$j.$sampling_point_init.dat > bs/bs/$i.$j.dat
      cat bs/bs.$i.$j.$sampling_point_final.dat >> bs/bs.$i.$j.dat            
      echo $band_energy >> bs/bs.$i.$j.dat
      rm bs/bs.$i.$j.$sampling_point_init.dat
      rm bs/bs.$i.$j.$sampling_point_final.dat
    fi
  done 
done

caesar calculate_bs               \
       $no_kpoints                \
       $no_modes                  \
       $band_degeneracy           \
       $sampling_point_amplitude  \
       ibz.dat                    \
       bs                         \
       bs/band_gap_correction.dat \
       bs/bg_correction_kp.dat
