#!/bin/bash

# Script to set-up a generic quadratic calculation 

# Read in harmonic path
echo "What is the path to the harmonic directory?"
read harmonic_path

# Read in number of supercells
no_sc=$( awk '{print $1}' no_sc.dat)

# Define some useful variables
sampling_amplitude=$( awk 'NR==1 {print $1}' mapping.dat)
sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)
no_sampling_points=$(( $sampling_point_final-$sampling_point_init ))
temperature=$( awk 'NR==3 {print $1}' mapping.dat)

atoms_line=$(awk -v IGNORECASE=1 '$1~/Atoms/{print NR}' \
   $harmonic_path/structure.dat)
symmetry_line=$(awk -v IGNORECASE=1 '$1~/Lattice/{print NR}' \
   $harmonic_path/structure.dat)
no_atoms=$($symmetry_line-$atoms_line-1)
no_modes=$(( $no_atoms*3 ))

# Write harmonic_path.dat
echo $harmonic_path > harmonic_path.dat

# Loop over supercells
for i in `seq 1 $no_sc`; do
  sdir=Supercell_$i
  mkdir $sdir
  
  no_atoms_sc=$( awk 'NR==1 {print $1}' \
     $harmonic_path/$sdir/super_equilibrium.dat )
  
  # Generate configurations

  # Generate static calculation
  mkdir $sdir/static
  
  echo Lattice                                   > $sdir/static/structure.dat
  cat $harmonic_path/$sdir/super_lattice.dat    >> $sdir/static/structure.dat
  echo Atoms                                    >> $sdir/static/structure.dat
  awk "NR==2,NR==$(( $no_atoms_sc+1 )) {print}" \
     $harmonic_path/$sdir/super_equilibrium.dat >> $sdir/static/structure.dat
  echo Symmetry                                 >> $sdir/static/structure.dat
  echo End                                      >> $sdir/static/structure.dat

  # loop over kpoints
  while read fline ; do
    line=$(fline)
    big_point=${line[0]}
    small_point=${line[1]}
    
    first_line=$(( (4+$no_atoms_sc)*($no_modes*($small_point-1))+1 ))
    last_line=$(( $first_line+(4+$no_atoms_sc)*$no_modes ))
    
    kdir=kpoint.$big_point
    mkdir $kdir
    awk "NR==$first_line,NR==$last_line {print}" \
       $harmonic_path/$sdir/lte/disp_patterns.dat > $kdir/disp_patterns.dat
    mkdir $kdir/configurations
    for j in `seq 1 $no_modes`; do
      frequency_line=$(( (4+$no_atoms_sc)*($j-1)+1 ))
      
      # write frequency in a.u. to frequency.$big_point.$j.dat
      frequency=$(awk "NR=$frequency_line {print\$3} $kdir/disp_patterns.dat")
      python -c "print($frequency*27.211396132)" > \
        $kdir/frequency.$big_point.$j.dat
      
      # write kpoint.$big_point/configurations/mode.j.k/structure.dat
      for k in `seq $sampling_point_init $sampling_point_final`; do
        mdir=$kdir/mode.$j.$k
        mkdir $mdir
        caesar generate_quadratic_configurations          \
               $sampling_amplitude                        \
               $k                                         \
               $no_sampling_points                        \
               $frequency                                 \
               $frequency_line                            \
               $harmonic_path/$sdir/super_equilibrium.dat \
               $harmonic_path/$sdir/super_lattice.dat     \
               $kdir/disp_patterns.dat                    \
               $mdir/structure.dat
      done # loop over sampling points per mode
    done # loop over modes
  done < $harmonic_path/$sdir/list.dat
done # Loop over supercells
