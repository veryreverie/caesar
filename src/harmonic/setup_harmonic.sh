#!/bin/bash

# Script to set-up a generic harmonic calculation to use with LTE

# Functions
function write_lattice
{
cat <<EOF
Lattice
EOF
}

function write_atoms
{
cat <<EOF
Atoms
EOF
}

function write_symmetry
{
cat <<EOF
Symmetry
EOF
}

function write_end
{
cat <<EOF
End
EOF
}




#######################################
####        MAIN   PROGRAM         ####
#######################################

# Read in 'structure.dat'
lattice_line=$(awk -v IGNORECASE=1 '$1~/Lattice/{print NR}' structure.dat) 
atoms_line=$(awk -v IGNORECASE=1 '$1~/Atoms/{print NR}' structure.dat)
symmetry_line=$(awk -v IGNORECASE=1 '$1~/Symmetry/{print NR}' structure.dat)
end_line=$(awk -v IGNORECASE=1 '$1~/End/{print NR}' structure.dat)

awk -v awk_lattice_line=$lattice_line -v awk_atoms_line=$atoms_line 'NR==(awk_lattice_line+1), NR==(awk_atoms_line-1) {print $1 " " $2 " " $3} ' structure.dat > lattice.dat
echo $(($symmetry_line-($atoms_line+1))) | bc > equilibrium.dat
awk -v awk_atoms_line=$atoms_line -v awk_symmetry_line=$symmetry_line 'NR==(awk_atoms_line+1), NR==(awk_symmetry_line-1) {print $1 " " $2 " " $3 " " $4 " " $5} ' structure.dat >> equilibrium.dat
echo $(($end_line-($symmetry_line+1)))/4 | bc > symmetry.dat
awk -v awk_symmetry_line=$symmetry_line -v awk_end_line=$end_line 'NR==(awk_symmetry_line+1), NR==(awk_end_line-1) {print $1 " " $2 " " $3} ' structure.dat >> symmetry.dat

# Generate IBZ
caesar generate_kgrid

# Generate non-diagonal supercells
caesar generate_supercells

CELL_COUNT=0
KPOINT_COUNT=0
OLD_CELL=0
NEW_CELL=0

while read LINE ; do

 KPOINT_COUNT=$(( ${KPOINT_COUNT} + 1 ))
 NEW_CELL=$(echo ${LINE} | awk '{print $4}')
 
 if (( ${NEW_CELL} != ${OLD_CELL} )) ; then

  CELL_COUNT=$(( ${CELL_COUNT} + 1 ))
  OLD_CELL=${NEW_CELL}
  SUPERCELL_DIR=Supercell_${CELL_COUNT}
  mkdir ${SUPERCELL_DIR}
  mv supercell.${CELL_COUNT}.dat ${SUPERCELL_DIR}/supercell.dat
  mv size.${CELL_COUNT}.dat ${SUPERCELL_DIR}/size.dat
  cp equilibrium.dat lattice.dat ${SUPERCELL_DIR}

  cd ${SUPERCELL_DIR}
  caesar construct_supercell
  cd ../

 fi

 KPOINT_1=$(echo ${LINE} | awk '{print $1}')
 KPOINT_2=$(echo ${LINE} | awk '{print $2}')
 KPOINT_3=$(echo ${LINE} | awk '{print $3}')

 echo ${KPOINT_COUNT} ${KPOINT_1} ${KPOINT_2} ${KPOINT_3} >> ${SUPERCELL_DIR}/kpoints.dat

done < kpoint_to_supercell.dat

echo $CELL_COUNT > no_sc.dat

# Generate Force Constants 
for (( i=1; i<=$CELL_COUNT; i++ ))do

  cd Supercell_$i

  # Construct symmetry operations
  write_lattice > lattice.txt
  write_atoms > atoms.txt
  write_symmetry > symmetry.txt
  write_end > end.txt
  sed '1d' super_equilibrium.dat > tmp.dat 
  cat lattice.txt super_lattice.dat atoms.txt tmp.dat symmetry.txt end.txt > structure.dat
  rm lattice.txt atoms.txt tmp.dat symmetry.txt end.txt
  caesar structure_to_castep
  cellsym --symmetry structure.cell > symmetry.dat
  symmetry_start_line=$(awk -v IGNORECASE=1 '/%block SYMMETRY_OPS/{print NR}' symmetry.dat)
  symmetry_end_line=$(awk -v IGNORECASE=1 '/%endblock SYMMETRY_OPS/{print NR}' symmetry.dat)
  awk -v awk_symmetry_start_line=$symmetry_start_line -v awk_symmetry_end_line=$symmetry_end_line 'NR==(awk_symmetry_start_line+1), NR==(awk_symmetry_end_line-1) {print $1 " " $2 " " $3} ' symmetry.dat >> symmetry_temp.dat
  echo $(( $symmetry_end_line-$symmetry_start_line ))/5 | bc > symm_header.txt
  cat symmetry_temp.dat | awk NF > symmetry.dat
  mv symmetry.dat symmetry_temp.dat
  cat symm_header.txt symmetry_temp.dat > symmetry.dat
  rm symmetry_temp.dat symm_header.txt
  
  # Evaluate needed force constants
  caesar construct_matrix_force_cnsts 
  while read LINE ; do

    echo $LINE > disp.dat
    atom=$(awk '{print $1}' disp.dat) 
    disp=$(awk '{print $2}' disp.dat) 
    mkdir atom.${atom}.disp.${disp}
    mv disp.dat atom.${atom}.disp.${disp}
    cp super_lattice.dat super_equilibrium.dat atom.${atom}.disp.${disp}
    cd atom.${atom}.disp.${disp}
    echo $atom $disp > displacement.dat
    caesar construct_finite_displacement
    mkdir positive negative
    mv positive.dat positive/structure.dat
    cp displacement.dat positive
    mv negative.dat negative/structure.dat
    cp displacement.dat negative
    cd ../

  done < force_constants.dat
  
  cd ../

done
