#!/bin/bash

no_sc=$(ls -1d Supercell_* | wc -l)

for i in `seq 1 $no_sc` ; do
  sdir=Supercell_$i
  
  caesar compare_kpoints             \
         $sdir/kpoints.dat           \
         $sdir/lte/gvectors_frac.dat \
         $sdir/list.dat
  
  cd $sdir/lte/
    caesar lte 0.00001 0.00001 0.01 > lte2.out
  cd -
  
  cp $sdir/lte/atoms_in_primitive_cell.dat lte/atoms_in_primitive_cell.$i.dat
  
  while read fline ; do
    line=($fline)
    big_point=${line[0]}
    small_point=${line[1]}
    cp $sdir/lte/dyn_mat.$small_point.dat lte/dyn_mat.$big_point.dat
  done < $sdir/list.dat

done
