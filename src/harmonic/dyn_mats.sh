#!/bin/bash

TOTAL=$(ls -1d Supercell_* | wc -l)

for i in $(seq 1 ${TOTAL}) ; do

 cd Supercell_${i}/

 cp lte/gvectors_frac.dat .

 if [[ ! -e "gvectors_frac.dat" ]] ; then
  echo "gvectors_frac.dat file does not exist."
  exit
 fi

 rm -f list.dat

 caesar compare_kpoints

 if [[ ! -e "list.dat" ]] ; then
  echo "Unable to generate list.dat file."
  exit
 fi

 rm -f gvectors_frac.dat

 if [[ ! -d "lte" ]] ; then
  echo "lte directory does not exist."
  exit
 fi

 cp list.dat lte/.

 cd lte/

 caesar lte > lte2.out
 
 cp atoms_in_primitive_cell.dat ../../lte/atoms_in_primitive_cell.${i}.dat


 while read LIST ; do
  BIG_POINT=$(echo ${LIST} | awk '{print $1}')
  SMALL_POINT=$(echo ${LIST} | awk '{print $2}')
  cp dyn_mat.${SMALL_POINT}.dat ../../lte/dyn_mat.${BIG_POINT}.dat
 done < list.dat

 rm -f list.dat

 cd ../../

done
