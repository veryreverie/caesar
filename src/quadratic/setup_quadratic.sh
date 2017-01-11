#!/bin/bash

# Script to set-up a generic quadratic calculation 

# Read in harmonic path
echo "What is the path to the harmonic directory?"
read harmonic_path

# Read in number of supercells
no_sc=$(ls -1d $harmonic_path/Supercell_* | wc -l)
if [[ $no_sc == 0 ]] then
  echo 'Error! No Supercell_* directories exist.'
  exit
fi

# Check relevant files exist
required_files=()
required_files+=("mapping.dat")
required_files+=("$harmonic_path/ibz.dat")
for i in `seq 1 $no_sc`; do
  sdir=Supercell_$i
  required_files+=("$harmonic_path/$sdir/list.dat")
  required_files+=("$harmonic_path/$sdir/supercell.dat")
  required_files+=("$harmonic_path/$sdir/lte/disp_patterns.dat")
done

for required_file in $required_files; do
  if [ ! -e $required_file ]; then
    echo "Error! the file $required_file does not exist."
    exit
  fi
done

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

cp $harmonic_path/ibz.dat .

# Loop over supercells
for i in `seq 1 $no_sc`; do

  sdir=Supercell_$i
  mkdir $sdir
  
  # copy files into Supercell_$i
  cp $harmonic_path/$sdir/list.dat $sdir
  cp $harmonic_path/$sdir/supercell.dat $sdir
  cp $harmonic_path/$sdir/super_lattice.dat $sdir
  cp $harmonic_path/$sdir/super_equilibrium.dat $sdir
  if [ -e "$harmonic_path/$sdir/KPOINTS.${i}" ]; then
    cp $harmonic_path/$sdir/KPOINTS.${i} $sdir
  fi
  cp $harmonic_path/$sdir/lte/disp_patterns.dat $sdir
  
  no_atoms_sc=$( awk 'NR==1 {print $1}' $sdir/super_equilibrium.dat )
  no_sc_atoms=$( awk 'NR==1 {print}' $sdir/super_equilibrium.dat )
  
  # Generate configurations

  # Generate static calculation
  mkdir $sdir/static
  cp $sdir/super_equilibrium.dat $sdir/super_lattice.dat $sdir/static
  f=$sdir/static/structure.dat
  echo Lattice > $f
  cat $sdir/super_lattice.dat >> $f
  echo Atoms >> $f
  awk "NR==2,NR==$(( $no_sc_atoms+1 )) {print}" $sdir/super_equilibrium.dat>>$f
  echo Symmetry >> $f
  echo End >> $f

  while read fline ; do # list.dat
    line=$(fline)
    big_point=${line[0]}
    small_point=${line[1]}
    kdir=kpoint.$big_point
    mkdir $kdir
    cp mapping.dat $sdir/super_equilibrium.dat $sdir/super_lattice.dat $kdir
    first_line=$(( (4+$no_atoms_sc)*($no_modes*($small_point-1))+1 ))
    last_line=$(( $first_line+(4+$no_atoms_sc)*$no_modes ))
    awk "NR==$first_line,NR==$last_line {print}" disp_patterns.dat > \
      $kdir/disp_patterns.dat
    mkdir $kdir/configurations
    for j in `seq 1 $no_modes`; do
      frequency_line=$(( (4+$no_atoms_sc)*($j-1)+1 ))
      
      # write frequency in a.u. to frequency.$big_point.$j.dat
      frequency=$(awk "NR=$frequency_line {print\$3} $kdir/disp_patterns.dat")
      python -c "print($frequency*27.211396132)" > \
        $kdir/frequency.$big_point.$j.dat
      
      # write structure.$j.$k.dat
      for k in `seq $sampling_point_init $sampling_point_final`; do
        caesar generate_quadratic_configurations \
               $sampling_amplitude               \
               $k                                \
               $no_sampling_points               \
               $frequency                        \
               $frequency_line                   \
               $kdir/configuration.dat           \
               $kdir/super_equilibrium.dat       \
               $kdir/disp_patterns.dat           \
               $kdir/positions.dat
        if [ -e "$kdir/positions.dat" ];then # if k.ne.0 and |frequency|>tol
          f=$kdir/configurations/structure.$j.$k.dat
          echo Lattice > $f
          cat $kdir/super_lattice.dat >> $f
          echo Atoms >> $f
          cat $kdir/positions.dat >> $f
          echo Symmetry >> $f
          echo End >> $f
          rm $kdir/positions.dat
        fi
      done # loop over sampling points per mode
    done # loop over modes
  done < $sdir/list.dat
done # Loop over supercells
