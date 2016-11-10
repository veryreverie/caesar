#!/bin/bash

# Script to construct and execute LTE


# Functions

function write_lte_bottom
{
cat <<EOF
 Program function: (1) calculate thermal energy; (2) calculate dispersion
  4
 Temperature (K)
  0.0
 Number of lines in k-space to plot
  4
 Points in journey through k-space (in terms of rec. prim. LVs)
0.000000 0.000000 0.000000  # GM
0.500000 0.500000 0.500000  # T
0.000000 0.500000 0.500000  # FB
0.000000 0.000000 0.000000  # GM
0.000000 0.500000 0.000000  # L
 Number of points per line
  400
EOF
}

function write_lte_path
{
cat <<EOF
0.000000 0.000000 0.000000  # GM
0.500000 0.500000 0.500000  # T
0.000000 0.500000 0.500000  # FB
0.000000 0.000000 0.000000  # GM
0.000000 0.500000 0.000000  # L
EOF
}


# MAIN PROGRAM

  echo "What temperature (K)?"
  read temperature

  no_sc=$(awk '{print}' no_sc.dat )
  
  # Loop over 
  for (( i=1; i<=$no_sc; i++ )) do
 
    cd Supercell_$i
  
    # Prepare force constants
    force_const=0
    while read LINE ; do

      force_const=$(( $force_const + 1 ))
      echo $LINE > disp.dat
      atom=$(awk '{print $1}' disp.dat)
      disp=$(awk '{print $2}' disp.dat)
      atoms_line=$(awk -v IGNORECASE=1 '/Atoms/{print NR}' structure.dat)
      symmetry_line=$(awk -v IGNORECASE=1 '/Symmetry/{print NR}' structure.dat)
      no_atoms=$(( $(( $symmetry_line-($atoms_line+1))) | bc ))
      
      cd atom.${atom}.disp.${disp}
      cp positive/forces.dat positive.dat
      cp negative/forces.dat negative.dat
      combine_forces
      cd ../

    done < force_constants.dat

    # Construct LTE
    echo "Primitive cell lattice vectors (rows, in a.u.)" > lte.dat
    cat lte.dat lattice.dat > lte_temp.dat
    mv lte_temp.dat lte.dat
    echo "Supercell lattice vectors (rows, in a.u.)" >> lte.dat
    cat lte.dat super_lattice.dat > lte_temp.dat
    mv lte_temp.dat lte.dat
    echo " Number of atoms in supercell" >> lte.dat
    echo $no_atoms >> lte.dat
    echo " Species ; mass (a.u.) ; position of atom in supercell (in terms of SC LVs)" >> lte.dat
    equilibrium_frac
    awk '{if (NR!=1) {print}}' super_equilibrium_frac.dat > atoms.dat
    cat lte.dat atoms.dat > lte_temp.dat
    rm atoms.dat
    mv lte_temp.dat lte.dat
    echo "Number of point symmetry operations" >> lte.dat
    awk '{if (NR==1) {print}}' symmetry.dat > symm.dat
    cat lte.dat symm.dat > lte_temp.dat
    rm symm.dat
    mv lte_temp.dat lte.dat
    echo "Rotation matrices (3 rows) and translations (1 row, fraction of SC LV)" >> lte.dat
    awk '{if (NR!=1) {print}}' symmetry.dat > symm.dat
    cat lte.dat symm.dat > lte_temp.dat
    rm symm.dat
    mv lte_temp.dat lte.dat
    echo "Number of force constants" >> lte.dat
    echo $(( $force_const*$no_atoms*3 )) >> lte.dat
    echo "Atom 1 ; Cartes'n direction ; Atom 2 ; C. direction ; force constant (a.u.)" >> lte.dat
    while read LINE ; do
      echo $LINE > disp.dat
      atom=$(awk '{print $1}' disp.dat)
      disp=$(awk '{print $2}' disp.dat)
      cat lte.dat atom.${atom}.disp.${disp}/forces.dat > lte_temp.dat
      mv lte_temp.dat lte.dat
    done < force_constants.dat
    write_lte_bottom > lte_bottom.dat
    cat lte.dat lte_bottom.dat > lte_temp.dat
    mv lte_temp.dat lte.dat
    rm lte_bottom.dat

    # Execute LTE
    mkdir lte
    mv lte.dat lte
    cd lte
    lte_lower > lte.out
    #lte > lte.out
    if [ -e "error.txt" ];then
      echo "There is an error in lte: check 'error.txt' file."
      exit 1
    fi
    cd ../

    cd ../

  done  # Loop over supercells


  # Collect relevant data
  mkdir lte
  cp lattice.dat equilibrium.dat symmetry.dat grid.dat ibz.dat kpoint_to_supercell.dat lte
  dyn_mats.sh
  cd lte
  write_lte_path > path.dat
  echo $temperature > temperature.dat
  fourier_interpolation > fourier_interpolation.out
  cd ../
 

