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
    sdir=Supercell_$i
  
    # Prepare force constants
    force_const=0
    while read fline ; do
      force_const=$(( $force_const + 1 ))
      
      line=($fline)
      disp=${line[0]}
      atom=${line[1]}
      echo $fline > $sdir/disp.dat
      atoms_line=$(awk -v IGNORECASE=1 '/Atoms/{print NR}' $sdir/structure.dat)
      symmetry_line=$(awk -v IGNORECASE=1 '/Symmetry/{print NR}' $sdir/structure.dat)
      no_atoms=$(( $(( $symmetry_line-($atoms_line+1))) | bc ))
      dir=$sdir/atom.$atom.disp.$disp
      caesar combine_forces                  \
             $dir/super_equilibrium_file.dat \
             $dir/positive/forces.dat        \
             $dir/negative/forces.dat        \
             $dir/forces.dat
      cp $dir/positive/forces.dat $dir/positive.dat
      cp $dir/negative/forces.dat $dir/negative.dat

    done < $sdir/force_constants.dat

    caesar equilibrium_frac            \
           $sdir/super_equilibrium.dat \
           $sdir/super_lattice.dat     \
           $sdir/super_equiibrium_frac.dat
    
    # Construct LTE
    mkdir $sdir/lte
    f=$sdir/lte/lte.dat
    
    echo "Primitive cell lattice vectors (rows, in a.u.)" > $f
    cat $sdir/lattice.dat >> $f
    echo "Supercell lattice vectors (rows, in a.u.)" >> $f
    cat $sdir/super_lattice.dat >> $f
    echo " Number of atoms in supercell" >> $f
    echo $no_atoms >> $f
    echo " Species ; mass (a.u.) ; position of atom in supercell (in terms of SC LVs)" >> $f
    awk '{if (NR!=1) {print}}' $sdir/super_equilibrium_frac.dat >> $f
    echo "Number of point symmetry operations" >> $f
    awk '{if (NR==1) {print}}' $sdir/symmetry.dat >> $f
    echo "Rotation matrices (3 rows) and translations (1 row, fraction of SC LV)" >> $f
    awk '{if (NR!=1) {print}}' $sdir/symmetry.dat >> $f
    echo "Number of force constants" >> $f
    echo $(( $force_const*$no_atoms*3 )) >> $f
    echo "Atom 1 ; Cartes'n direction ; Atom 2 ; C. direction ; force constant (a.u.)" >> $f
    while read fline ; do
      line=($fline)
      atom=${line[0]}
      disp=${line[1]}
      echo $fline > $sdir/disp.dat
      cat atom.${atom}.disp.${disp}/forces.dat >> $f
    done < $sdir/force_constants.dat
    write_lte_bottom >> $f

    # Execute LTE
    cd $sdir/lte
      caesar lte 0.00000001 0.001 0.000001 > lte.out
      if [ -e "error.txt" ];then
        echo "There is an error in lte: check 'error.txt' file."
        exit 1
      fi
    cd -
  done  # Loop over supercells


  # Collect relevant data
  mkdir lte
  cp lattice.dat equilibrium.dat symmetry.dat grid.dat ibz.dat kpoint_to_supercell.dat lte
  caesar dyn_mats
  write_lte_path > lte/path.dat
  echo $temperature > lte/temperature.dat
  cd lte
    caesar fourier_interpolation > fourier_interpolation.out
  cd ../
