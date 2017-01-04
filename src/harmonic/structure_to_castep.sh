#!/bin/bash

# Script to transform structure file to cell file
# $1=structure.dat, $2=sc_bs_path.dat, $3=structure.cell

tempfile=structure_to_castep.temp

# move an existing structure.cell to tempfile, to allow it to be prepended
if [ -f $3 ]
  mv $3 $tempfile
fi

# Read in 'structure.dat'
lattice_line=$(awk -v IGNORECASE=1 '$1~/Lattice/{print NR}' $1)
atoms_line=$(awk -v IGNORECASE=1 '$1~/Atoms/{print NR}' $1)
symmetry_line=$(awk -v IGNORECASE=1 '$1~/Symmetry/{print NR}' $1)
end_line=$(awk -v IGNORECASE=1 '$1~/End/{print NR}' $1)

echo '%block lattice_cart' > $3
echo 'bohr' >> $3
awk "NR==$(($lattice_line + 1)),NR==$(($atoms_line - 1)) " \
  '{print $1 " " $2 " " $3}' $1 >> $3
echo '%endblock lattice_cart' >> $3
echo '' >> $3
echo '%block positions_abs' >> $3
echo 'bohr' >> $3
awk "NR==$(($atoms_line + 1)),NR==$(($symmetry_line - 1)) " \
  '{print $1 " " $3 " " $4 " " $5}' $1 >> $3
echo '%endblock positions_abs' >> $3
echo '' >> $3

if [ -f "$tempfile" ]; then
  cat $tempfile >> $3
  rm $tempfile
fi

if [ -f "$2" ]; then
  echo '' >> $3
  echo '%block bs_kpoint_path' >> $3
  cat $2 >> $3
  echo '%endblock bs_kpoint_path' >> $3
  echo '' >> $3
  echo 'bs_kpoint_path_spacing = 10.0' >> $3
fi
