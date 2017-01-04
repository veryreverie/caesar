#!/bin/bash

# Script to transform structure file to cell file
# $1=structure.dat $2=structure.POSCAR

# Read in 'structure.dat'
lattice_line=$(awk -v IGNORECASE=1 '$1~/Lattice/{print NR}' $1) 
atoms_line=$(awk -v IGNORECASE=1 '$1~/Atoms/{print NR}' $1)
symmetry_line=$(awk -v IGNORECASE=1 '$1~/Symmetry/{print NR}' $1)
end_line=$(awk -v IGNORECASE=1 '$1~/End/{print NR}' $1)

echo 'Structure' > $2
echo '0.52917721092' >> $2
awk "NR==$(($lattice_line + 1)),NR==$(($atoms_line - 1)) " \
  '{print $1 " " $2 " " $3}' $1 >> $2

number_of_species=0
species_pre="X"
declare -A species_counter
declare -A species_type
for (( i=$atoms_line+1; i<=$symmetry_line-1; i++)) do
  species=$(awk "NR==$i {print \$1} $1")
  if [ "$species" !=  "$species_pre" ]; then
    number_of_species=$(( $number_of_species+1 ))
    species_pre=$species
    species_counter[${number_of_species}]=0
    species_type[${number_of_species}]=$species
  fi
  species_counter[${number_of_species}]=$(( ${species_counter[${number_of_species}]}+1 ))
done

echo ${species_type[*]} >> $2
echo ${species_counter[*]} >> $2
echo 'Cartesian' >> $2
awk "NR==$(($atoms_line + 1)),NR==$(($symmetry_line - 1)) " \
  '{print $3 " " $4 " " $5}' $1 >> $2
