#!/bin/bash

# Script to set-up a generic harmonic calculation to use with LTE

#######################################
####        MAIN   PROGRAM         ####
#######################################

# Read in 'structure.dat'
lattice_line=$(awk -v IGNORECASE=1 '$1~/Lattice/{print NR}' structure.dat) 
atoms_line=$(awk -v IGNORECASE=1 '$1~/Atoms/{print NR}' structure.dat)
symmetry_line=$(awk -v IGNORECASE=1 '$1~/Symmetry/{print NR}' structure.dat)
end_line=$(awk -v IGNORECASE=1 '$1~/End/{print NR}' structure.dat)

# break structure.dat into lattice.dat, equilibrium.dat and symmetry.dat
awk "NR==$(($lattice_line + 1)), NR==$(($atoms_line - 1)) "\
  '{print $1 " " $2 " " $3} ' structure.dat > lattice.dat

echo $(($symmetry_line-($atoms_line+1))) | bc > equilibrium.dat
awk "NR==$(($atoms_line + 1)), NR==$(($symmetry_line - 1)) "\
  '{print $1 " " $2 " " $3 " " $4 " " $5}' structure.dat >> equilibrium.dat

echo $(($end_line-($symmetry_line+1)))/4 | bc > symmetry.dat
awk "NR==$(($symmetry_line + 1)), NR==$(($end_line - 1)) "\
  '{print $1 " " $2 " " $3}' structure.dat >> symmetry.dat

# Generate IBZ
caesar generate_kgrid \
       symmetry.dat   \
       grid.dat       \
       lattice.dat    \
       ibz.dat        \
       rotated_gvectors.dat

# Generate non-diagonal supercells
caesar generate_supercells

CELL_COUNT=0
KPOINT_COUNT=0
OLD_CELL=0
NEW_CELL=0

while read fline ; do
  line=($fline)
  
  KPOINT_COUNT=$(( ${KPOINT_COUNT} + 1 ))
  NEW_CELL=${line[3]}
  
  if (( ${NEW_CELL} != ${OLD_CELL} )) ; then
  
    CELL_COUNT=$(( ${CELL_COUNT} + 1 ))
    OLD_CELL=${NEW_CELL}
    sdir=Supercell_${CELL_COUNT}
    cp equilibrium.dat lattice.dat $sdir
    
    caesar construct_supercell     \
           $sdir/lattice.dat       \
           $sdir/equilibrium.dat   \
           $sdir/supercell.dat     \
           $sdir/super_lattice.dat \
           $sdir/super_equilibrium.dat
    
  fi
  
  KPOINT_1=${line[0]}
  KPOINT_2=${line[1]}
  KPOINT_3=${line[2]}
  
  echo ${KPOINT_COUNT} ${KPOINT_1} ${KPOINT_2} ${KPOINT_3} >> $sdir/kpoints.dat
  
done < kpoint_to_supercell.dat

echo $CELL_COUNT > no_sc.dat

# Generate Force Constants 
for (( i=1; i<=$CELL_COUNT; i++ ))do

  sdir=Supercell_$i

  # Construct symmetry operations
  
  echo Lattice > $sdir/structure.dat
  cat $sdir/super_lattice.dat >> $sdir/structure.dat
  echo Atoms >> $sdir/structure.dat
  sed '1d' $sdir/super_equilibrium.dat >> $sdir/structure.dat #all but 1st line
  echo Symmetry >> $sdir/structure.dat
  echo End >> $sdir/structure.dat
  
  caesar structure_to_castep  \
         $sdir/structure.dat  \
         $sdir/sc_bs_path.dat \
         $sdir/structure.cell
  cellsym --symmetry $sdir/structure.cell > $sdir/symmetry.dat
  symmetry_start_line=$(awk -v IGNORECASE=1 '/%block SYMMETRY_OPS/{print NR}' $sdir/symmetry.dat)
  symmetry_end_line=$(awk -v IGNORECASE=1 '/%endblock SYMMETRY_OPS/{print NR}' $sdir/symmetry.dat)
  awk "NR==$(($symmetry_start_line + 1)),NR==$(($symmetry_end_line - 1)) "\
    '{print $1 " " $2 " " $3}' $sdir/symmetry.dat > $sdir/symmetry_temp.dat
  echo $(( $symmetry_end_line-$symmetry_start_line ))/5 | bc > $sdir/symmetry.dat
  cat $sdir/symmetry_temp.dat | awk NF >> $sdir/symmetry.dat
  rm $sdir/symmetry_temp.dat
  
  # Evaluate needed force constants
  caesar construct_matrix_force_cnsts \
         $sdir/symmetry.dat           \
         $sdir/lattice.dat            \
         $sdir/supercell.dat          \
         $sdir/super_equilibrium.dat  \
         $sdir/force_constants.dat
  while read fline ; do
    # TODO : disp.dat and {displacement.dat} are the same files
    line=($fline)
    atom=${line[0]}
    disp=${line[1]}
    ddir=$sdir/atom.$atom.disp.$disp
    mkdir $ddir
    mkdir $ddir/positive $ddir/negative
    echo $fline > $ddir/disp.dat
    echo $fline > $ddir/displacement.dat
    echo $fline > $ddir/positive/displacement.dat
    echo $fline > $ddir/negative/displacement.dat
    cp $sdir/super_lattice.dat $sdir/super_equilibrium.dat $ddir
    caesar construct_finite_displacement \
           $atom                         \
           $disp                         \
           $ddir/super_lattice.dat       \
           $ddir/super_equilibrium.dat   \
           $ddir/positive/structure.dat  \
           $ddir/negative/structure.dat

  done < $sdir/force_constants.dat
done
