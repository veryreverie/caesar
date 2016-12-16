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

echo '%block lattice_cart' > $1/structure.cell
echo 'bohr' >> $1/structure.cell
awk "NR==$(($lattice_line + 1)),NR==$(($atoms_line - 1)) " \
  '{print $1 " " $2 " " $3}' $1/structure.dat >> $1/structure.cell
echo '%endblock lattice_cart' >> $1/structure.cell
echo '' >> $1/structure.cell
echo '%block positions_abs' >> $1/structure.cell
echo 'bohr' >> $1/structure.cell
awk "NR==$(($atoms_line + 1)),NR==$(($symmetry_line - 1)) " \
  '{print $1 " " $3 " " $4 " " $5}' $1/structure.dat >> $1/structure.cell
echo '%endblock positions_abs' >> $1/structure.cell
echo '' >> $1/structure.cell
if [ -f "$1/bottom.cell" ]; then
  cat $1/bottom.cell >> $1/structure.cell
fi
if [ -f "$1/sc_bs_path.dat" ]; then
  echo '' >> $1/structure.cell
  echo '%block bs_kpoint_path' >> $1/structure.cell
  cat $1/sc_bs_path.dat >> $1/structure.cell
  echo '%endblock bs_kpoint_path' >> $1/structure.cell
  echo '' >> $1/structure.cell
  echo 'bs_kpoint_path_spacing = 10.0' >> $1/structure.cell
fi

