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

harmonic_path=$( awk '{print $1}' harmonic_path.dat)

no_sc=$( awk '{print $1}' no_sc.dat)

sampling_point_amplitude=$( awk 'NR==1 {print $1}' mapping.dat)
sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)

atoms_line=$(awk -v IGNORECASE=1 '$1~/Atoms/{print NR}' \
   $harmonic_path/structure.dat)
symmetry_line=$(awk -v IGNORECASE=1 '$1~/Lattice/{print NR}' \
   $harmonic_path/structure.dat)
no_atoms=$($symmetry_line-$atoms_line-1)
no_modes=$(( $no_atoms*3 ))

seedname=$( awk '{print $1}' seedname.txt )

kpoint_counter=1
mkdir bs

# Loop over supercells
for sc_id in `seq 1 $no_sc`; do
  sdir=Supercell_$sc_id
  if [ ! -f "$sdir/acoustic.dat" ]; then
    # Static calculation
    static_dir=$sdir/static
    caesar eigenval_castep_to_bands $static_dir/$seedname.bands $static_dir
  fi
done

# Loop over kpoints
while read fline ; do
  line=($fline)
  kpoint2=${line[0]}
  sc_id=${line[2]}
  sdir=Supercell_$sc_id
  
  if [ ! -f "$sdir/acoustic.dat" ]; then
    for j in `seq 1 $no_modes`; do
      for k in `seq $sampling_point_init $sampling_point_final`; do
        mdir=$sdir/kpoint.$kpoint2/configurations/mode.$j.$k
        if [ -e "$mdir/$seedname.castep" ]; then
          caesar eigenval_castep_to_bands $mdir/$seedname.bands $mdir
        fi
      done # loop over sampling points per mode
    done # loop over modes
  fi
done < $harmonic_path/list.dat

# Loop over supercells
for sc_id in `seq 1 $no_sc`; do
  sdir=Supercell_$sc_id
  
  # Obtain relevant band for each supercell
  band_energy=$(awk "NR==$band {print}" $sdir/static/kpoint.$kpoint.dat)
  
  caesar band_folding                   \
         $static_dir/kpoint.$kpoint.dat \
         $band_energy                   \
         $static_dir/band_number.dat
  
  band_ref=$( awk '{print}' $sdir/static/band_number.dat )
  
  # Collect quadratic results
  no_kpoints=$(ls -1d $sdir/kpoint.* | wc -l)
  
  # Loop over k-points
  for j in `seq $kpoint_counter $(( $kpoint_counter+($no_kpoints-1) ))`; do
    kdir=kpoint.$j/configurations
    for k in `seq 1 $no_modes`; do
      cp kpoint.$j/frequency.$k.dat bs
      for l in `seq $sampling_point_init $sampling_point_final`; do
        mdir=$kdir/mode.$k.$l
        if [ -e "$mdir/$seedname.castep" ] && [ "$l" -ne "0" ]; then
          first_line=$(( $band_ref-($band_degeneracy-1) )) 
          last_line=$(( $band_ref ))
          if [ "$l" = "$sampling_point_init" ]; then
            awk "NR==$first_line,NR==$last_line" $mdir/kpoint.$kpoint.dat \
              > kpoint.$j/bs.$k.dat
          elif [ "$l" = "$sampling_point_final" ]; then
            awk "NR==$first_line,NR==$last_line" $mdir/kpoint.$kpoint.dat \
              >> kpoint.$j/bs.$k.dat
            echo $band_energy >> kpoint.$j/bs.$k.dat
          fi
        fi
      done
    done
  done
  kpoint_counter=$(( $kpoint_counter+$no_kpoints ))
done

caesar calculate_bs               \
       $no_modes                  \
       $band_degeneracy           \
       $sampling_point_amplitude  \
       $harmonic_path/ibz.dat     \
       bs                         \
       bs/band_gap_correction.dat \
       bs/bg_correction_kp.dat
