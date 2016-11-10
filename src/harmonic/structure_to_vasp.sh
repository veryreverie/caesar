#!/bin/bash

# Script to transform structure file to cell file

#######################################
####        MAIN   PROGRAM         ####
#######################################

# Read in 'structure.dat'
lattice_line=$(awk -v IGNORECASE=1 '$1~/Lattice/{print NR}' structure.dat) 
atoms_line=$(awk -v IGNORECASE=1 '$1~/Atoms/{print NR}' structure.dat)
symmetry_line=$(awk -v IGNORECASE=1 '$1~/Symmetry/{print NR}' structure.dat)
end_line=$(awk -v IGNORECASE=1 '$1~/End/{print NR}' structure.dat)

echo 'Structure' > structure.POSCAR
echo '0.52917721092' >> structure.POSCAR
awk -v awk_lattice_line=$lattice_line -v awk_atoms_line=$atoms_line 'NR==(awk_lattice_line+1), NR==(awk_atoms_line-1) {print $1 " " $2 " " $3} ' structure.dat >> structure.POSCAR

number_of_species=0
species_pre="X"
declare -A species_counter
declare -A species_type
for (( i=$atoms_line+1; i<=$symmetry_line-1; i++)) do
  species=$(awk -v awk_line=$i 'NR==awk_line {print $1}' structure.dat)
  if [ "$species" !=  "$species_pre" ]; then
    number_of_species=$(( $number_of_species+1 ))
    species_pre=$species
    species_counter[${number_of_species}]=0
    species_type[${number_of_species}]=$species
  fi
  species_counter[${number_of_species}]=$(( ${species_counter[${number_of_species}]}+1 ))
done

echo ${species_type[*]} >> structure.POSCAR
echo ${species_counter[*]} >> structure.POSCAR
echo 'Cartesian' >> structure.POSCAR
awk -v awk_atoms_line=$atoms_line -v awk_symmetry_line=$symmetry_line 'NR==(awk_atoms_line+1), NR==(awk_symmetry_line-1) {print $3 " " $4 " " $5} ' structure.dat >> structure.POSCAR

