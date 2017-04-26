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

echo '%block lattice_cart' > structure.cell
echo 'bohr' >> structure.cell
awk -v awk_lattice_line=$lattice_line -v awk_atoms_line=$atoms_line 'NR==(awk_lattice_line+1), NR==(awk_atoms_line-1) {print $1 " " $2 " " $3} ' structure.dat >> structure.cell
echo '%endblock lattice_cart' >> structure.cell
echo '' >> structure.cell
echo '%block positions_abs' >> structure.cell
echo 'bohr' >> structure.cell
awk -v awk_atoms_line=$atoms_line -v awk_symmetry_line=$symmetry_line 'NR==(awk_atoms_line+1), NR==(awk_symmetry_line-1) {print $1 " " $3 " " $4 " " $5} ' structure.dat >> structure.cell
echo '%endblock positions_abs' >> structure.cell
echo '' >> structure.cell
if [ -f 'bottom.cell' ]; then
  cat structure.cell bottom.cell > structure_temp.cell
  mv structure_temp.cell structure.cell
fi
if [ -f 'sc_bs_path.dat' ]; then
  echo '' >> structure.cell
  echo '%block bs_kpoint_path' >> structure.cell
  cat structure.cell sc_bs_path.dat > structure_temp.cell
  mv structure_temp.cell structure.cell
  echo '%endblock bs_kpoint_path' >> structure.cell
  echo '' >> structure.cell
  echo 'bs_kpoint_path_spacing = 10.0' >> structure.cell
fi

