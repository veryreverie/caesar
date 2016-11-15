#!/bin/bash

# Script to set-up a generic quadratic calculation 

# Functions
function write_lattice
{
cat <<EOF
Lattice
EOF
}

function write_atoms
{
cat <<EOF
Atoms
EOF
}

function write_symmetry
{
cat <<EOF
Symmetry
EOF
}

function write_end
{
cat <<EOF
End
EOF
}




#######################################
####        MAIN   PROGRAM         ####
#######################################

# Define some useful variables
echo "What is the path to the harmonic directory?"
read harmonic_path

no_sc=$(ls -1d $harmonic_path/Supercell_* | wc -l)

sampling_amplitude=$( awk 'NR==1 {print $1}' mapping.dat)
sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)
no_sampling_points=$(( $sampling_point_final-$sampling_point_init ))
temperature=$( awk 'NR==3 {print $1}' mapping.dat)

cp $harmonic_path/ibz.dat .

# Loop over supercells
for (( i=1; i<=$no_sc; i++ )) do

  mkdir Supercell_$i
  if [ -e "$harmonic_path/Supercell_$i/list.dat" ]; then
    cp $harmonic_path/Supercell_$i/list.dat Supercell_$i
  else
    echo 'Error! The file list.dat does not exist for supercell' $i
    exit
  fi
  if [ -e "$harmonic_path/Supercell_$i/supercell.dat" ]; then
    cp $harmonic_path/Supercell_$i/supercell.dat Supercell_$i
  else
    echo 'Error! The file supercell.dat does not exist for supercell' $i
    exit
  fi
  cp $harmonic_path/Supercell_$i/lattice.dat Supercell_$i
  cp $harmonic_path/Supercell_$i/equilibrium.dat Supercell_$i
  no_atoms=$( awk 'NR==1 {print $1}' Supercell_$i/equilibrium.dat )
  no_modes=$(( $no_atoms*3 ))
  cp $harmonic_path/Supercell_$i/super_lattice.dat Supercell_$i
  cp $harmonic_path/Supercell_$i/super_equilibrium.dat Supercell_$i
  cp $harmonic_path/Supercell_$i/super_equilibrium.dat Supercell_$i
  if [ -e "$harmonic_path/Supercell_$i/KPOINTS.${i}" ]; then
    cp $harmonic_path/Supercell_$i/KPOINTS.${i} Supercell_$i
  fi
  no_atoms_sc=$( awk 'NR==1 {print $1}' Supercell_$i/super_equilibrium.dat )
  if [ -e "$harmonic_path/Supercell_$i/lte/disp_patterns.dat" ]; then
    cp $harmonic_path/Supercell_$i/lte/disp_patterns.dat Supercell_$i
  else
    echo 'Error! The file disp_patterns.dat do not exist for supercell' $i
    exit
  fi
  
  # Generate configurations
  cd Supercell_$i

  # Generate static calculation
  mkdir static
  cp super_equilibrium.dat super_lattice.dat static
  cd static
    write_lattice > lattice.txt
    write_atoms > atoms.txt
    write_symmetry > symmetry.txt
    write_end > end.txt
    no_sc_atoms=$( awk 'NR==1 {print}' super_equilibrium.dat )
    awk -v awk_last_line=$(( $no_sc_atoms+1 )) 'NR==2,NR==awk_last_line {print}' super_equilibrium.dat > positions.dat
    cat lattice.txt super_lattice.dat atoms.txt positions.dat symmetry.txt end.txt > structure.dat
    rm positions.dat lattice.txt atoms.txt symmetry.txt end.txt 
  cd ../ # static


  while read line ; do
    big_point=$(echo ${line} | awk '{print $1}')
    small_point=$(echo ${line} | awk '{print $2}')
    mkdir kpoint.$big_point
    cp ../mapping.dat super_equilibrium.dat super_lattice.dat kpoint.$big_point
    disp_patt_first_line=$(( (4+$no_atoms_sc)*($no_modes*($small_point-1))+1 ))
    disp_patt_last_line=$(( $disp_patt_first_line+(4+$no_atoms_sc)*$no_modes ))
    awk -v awk_disp_patt_first_line=$disp_patt_first_line -v awk_disp_patt_last_line=$disp_patt_last_line 'NR==awk_disp_patt_first_line,NR==awk_disp_patt_last_line {print}' disp_patterns.dat > disp_patterns_temp.dat
    mv disp_patterns_temp.dat kpoint.$big_point/disp_patterns.dat
    cd kpoint.$big_point
      write_lattice > lattice.txt
      write_atoms > atoms.txt
      write_symmetry > symmetry.txt
      write_end > end.txt
      mkdir configurations
      for j in `seq 1 $no_modes`; do
        disp_patt_first_line=$(( (4+$no_atoms_sc)*($j-1)+1 ))
        disp_patt_last_line=$(( (4+$no_atoms_sc)*($j-1)+1+(4+$no_atoms_sc) ))
        awk -v awk_disp_patt_first_line=$disp_patt_first_line -v awk_disp_patt_last_line=$disp_patt_last_line 'NR==awk_disp_patt_first_line,NR==awk_disp_patt_last_line {print}' disp_patterns.dat > disp_patterns_temp.dat
        for k in `seq $sampling_point_init $sampling_point_final`; do
          echo $sampling_amplitude $k $no_sampling_points > configuration.dat
          caesar generate_quadratic_configurations 
          mv frequency.dat frequency.${big_point}.${j}.dat
          if [ -e 'positions.dat' ];then 
            cat lattice.txt super_lattice.dat atoms.txt positions.dat symmetry.txt end.txt > configurations/structure.${j}.${k}.dat    
            rm positions.dat
          fi
          rm configuration.dat 
        done # loop over sampling points per mode
        rm disp_patterns_temp.dat
      done # loop over modes
      rm lattice.txt atoms.txt symmetry.txt end.txt
    cd ../

  done < list.dat
 
  cd ../

done # Loop over supercells

