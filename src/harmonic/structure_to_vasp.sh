#!/bin/bash

# Script to transform structure file to cell file

#######################################
####        MAIN   PROGRAM         ####
#######################################

# Read in 'structure.dat'
lattice_line=$(awk -v IGNORECASE=1 '$1~/Lattice/{print NR}' $1/structure.dat) 
atoms_line=$(awk -v IGNORECASE=1 '$1~/Atoms/{print NR}' $1/structure.dat)
symmetry_line=$(awk -v IGNORECASE=1 '$1~/Symmetry/{print NR}' $1/structure.dat)
end_line=$(awk -v IGNORECASE=1 '$1~/End/{print NR}' $1/structure.dat)

echo 'Structure' > $1/structure.POSCAR
echo '0.52917721092' >> $1/structure.POSCAR
awk "NR==$(($lattice_line + 1)),NR==$(($atoms_line - 1)) " \
  '{print $1 " " $2 " " $3}' $1/structure.dat >> $1/structure.POSCAR

number_of_species=0
species_pre="X"
declare -A species_counter
declare -A species_type
for (( i=$atoms_line+1; i<=$symmetry_line-1; i++)) do
  species=$(awk "NR==$i {print \$1} $1/structure.dat")
  if [ "$species" !=  "$species_pre" ]; then
    number_of_species=$(( $number_of_species+1 ))
    species_pre=$species
    species_counter[${number_of_species}]=0
    species_type[${number_of_species}]=$species
  fi
  species_counter[${number_of_species}]=$(( ${species_counter[${number_of_species}]}+1 ))
done

echo ${species_type[*]} >> $1/structure.POSCAR
echo ${species_counter[*]} >> $1/structure.POSCAR
echo 'Cartesian' >> $1/structure.POSCAR
awk "NR==$(($atoms_line + 1)),NR==$(($symmetry_line - 1)) " \
  '{print $3 " " $4 " " $5}' $1/structure.dat >> $1/structure.POSCAR
