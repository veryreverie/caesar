#!/bin/bash

no_sc=$(ls -1d Supercell_* | wc -l)

for i in `seq 1 $no_sc` ; do
  sdir=Supercell_$i
  
  caesar compare_kpoints             \
         $sdir/kpoints.dat           \
         $sdir/lte/gvectors_frac.dat \
         $sdir/list.dat
  
  caesar lte                                      \
         0.00001                                  \
         0.00001                                  \
         0.01                                     \
         $sdir/lte/lte.dat                        \
         $sdir/lte/freq_dos.dat                   \
         $sdir/lte/tdependence1.dat               \
         $sdir/lte/tdependence2.dat               \
         $sdir/lte/dispersion_curve.dat           \
         $sdir/lte/kpairs.dat                     \
         $sdir/lte/freq_grids.dat                 \
         $sdir/lte/disp_patterns.dat              \
         $sdir/lte/kdisp_patterns.dat             \
         $sdir/lte/pol_vec.dat                    \
         $sdir/lte/gvectors.dat                   \
         $sdir/lte/gvectors_frac.dat              \
         $sdir/lte/error.txt                      \
         $sdir/lte/dyn_mat.                       \
         $sdir/lte/atoms_in_primitive_cell.$i.dat \
         > $sdir/lte/lte2.out
  
  while read fline ; do
    line=($fline)
    big_point=${line[0]}
    small_point=${line[1]}
    cp $sdir/lte/dyn_mat.$small_point.dat lte/dyn_mat.$big_point.dat
  done < $sdir/list.dat

done
