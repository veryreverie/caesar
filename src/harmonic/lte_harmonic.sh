#!/bin/bash

# Script to construct and execute LTE

echo "What temperature (K)?"
read temperature

no_sc=$(awk '{print}' no_sc.dat )

lattice_line=$(awk -v IGNORECASE=1 '$1~/Lattice/{print NR}' structure.dat) 
atoms_line=$(awk -v IGNORECASE=1 '/Atoms/{print NR}' structure.dat)
symmetry_line=$(awk -v IGNORECASE=1 '/Symmetry/{print NR}' structure.dat)
no_atoms=$(( $(( $symmetry_line-($atoms_line+1))) | bc ))

code=$(awk '{print}' code.txt)

mkdir lte

# Loop over supercells
for (( i=1; i<=$no_sc; i++ )) do
  sdir=Supercell_$i
  
  super_lattice_line=$(awk -v IGNORECASE=1 '$1~/Lattice/{print NR}' $sdir/structure.dat) 
  super_atoms_line=$(awk -v IGNORECASE=1 '/Atoms/{print NR}' $sdir/structure.dat)
  super_symmetry_line=$(awk -v IGNORECASE=1 '/Symmetry/{print NR}' $sdir/structure.dat)
  super_end_line=$(awk -v IGNORECASE=1 '/End/{print NR}' $sdir/structure.dat)
  super_no_atoms=$(($super_symmetry_line-$super_atoms_line-1))
  super_no_symmetries=$(($(($super_end_line-$super_symmetry_line-1))/5 | bc))

  # Prepare force constants
  force_const=0
  while read fline ; do
    force_const=$(( $force_const + 1 ))
    
    line=($fline)
    disp=${line[0]}
    atom=${line[1]}
    ddir=$sdir/atom.$atom.disp.$disp
    
    paths=(positive negative)
    for path in ${paths[@]}; do
      if [ "$code" = "castep" ]; then
        caesar fetch_forces                 \
               $code                        \
               $ddir/$path/$seedname.castep \
               $atom                        \
               $disp                        \
               $ddir/$path.dat
      elif [ "$code" = "qe" ]; then
        caesar fetch_forces              \
               $code                     \
               $ddir/$path/$seedname.out \
               $atom                     \
               $disp                     \
               $ddir/$path.dat
      fi
    done
    
    caesar combine_forces      \
           $sdir/structure.dat \
           $ddir/positive.dat  \
           $ddir/negative.dat  \
           $ddir/forces.dat

  done < $sdir/force_constants.dat

  caesar equilibrium_frac    \
         $sdir/structure.dat \
         $sdir/super_equiibrium_frac.dat
  
  # Construct LTE
  mkdir $sdir/lte
  f=$sdir/lte/lte.dat
  
  echo "Primitive cell lattice vectors (rows, in a.u.)"                   > $f
  awk "NR==$(($lattice_line + 1)), NR==$(($atoms_line - 1)) "\
     '{print $1 " " $2 " " $3} ' structure.dat                           >> $f
  echo "Supercell lattice vectors (rows, in a.u.)"                       >> $f
  awk "NR==$(($super-lattice_line + 1)), NR==$(($super_atoms_line - 1)) "\
     '{print $1 " " $2 " " $3} ' $sdir/structure.dat                     >> $f
  echo " Number of atoms in supercell"                                   >> $f
  echo $super_no_atoms                                                   >> $f
  echo " Species ; mass (a.u.) ; position of atom in supercell (in terms of SC LVs)" >> $f
  awk '{if (NR!=1) {print}}' $sdir/super_equilibrium_frac.dat            >> $f
  echo "Number of point symmetry operations"                             >> $f
  echo $super_no_symmetries                                              >> $f
  echo "Rotation matrices (3 rows) and translations (1 row, fraction of SC LV)" >> $f
  awk "NR==$(($super_symmetry_line + 1)), NR==$(($super_end_line - 1)) "\
     '{print $1 " " $2 " " $3} ' $sdir/structure.dat                     >> $f
  echo "Number of force constants"                                       >> $f
  echo $(( $force_const*$super_no_atoms*3 ))                             >> $f
  echo "Atom 1 ; Cartes'n direction ; Atom 2 ; C. direction ; force constant (a.u.)" >> $f
  while read fline ; do
    line=($fline)
    atom=${line[0]}
    disp=${line[1]}
    cat atom.${atom}.disp.${disp}/forces.dat                             >> $f
  done < $sdir/force_constants.dat
  echo " Program function: (1) calculate thermal energy; (2) calculate dispersion" >> $f
  echo "  4"                                                             >> $f
  echo " Temperature (K)"                                                >> $f
  echo "  0.0"                                                           >> $f
  echo " Number of lines in k-space to plot"                             >> $f
  echo "  4"                                                             >> $f
  echo " Points in journey through k-space (in terms of rec. prim. LVs)" >> $f
  echo "0.000000 0.000000 0.000000  # GM"                                >> $f
  echo "0.500000 0.500000 0.500000  # T"                                 >> $f
  echo "0.000000 0.500000 0.500000  # FB"                                >> $f
  echo "0.000000 0.000000 0.000000  # GM"                                >> $f
  echo "0.000000 0.500000 0.000000  # L"                                 >> $f
  echo " Number of points per line"                                      >> $f
  echo "  400"                                                           >> $f
  
  # Execute LTE
  # Executes with mode 4.
  # Reads lte.dat
  # Writes all other files
  # Writes dyn_mat.*.dat for * in [1,number of gvectors]
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
  if [ -e "error.txt" ];then
    echo "There is an error in lte: check 'error.txt' file."
    exit 1
  fi
  
  # Collect relevant data
  caesar compare_kpoints             \
         $i                          \
         ibz.dat                     \
         $sdir/lte/gvectors_frac.dat \
         list.dat
done
  
while read fline ; do
  line=($fline)
  kpoint=${line[0]}
  gvector=${line[1]}
  sc_id=${line[2]}
  cp Supercell_$sc_id/lte/dyn_mat.$gvector.dat lte/dyn_mat.$kpoint.dat
done < list.dat

echo "0.000000 0.000000 0.000000  # GM"  > lte/path.dat
echo "0.500000 0.500000 0.500000  # T"  >> lte/path.dat
echo "0.000000 0.500000 0.500000  # FB" >> lte/path.dat
echo "0.000000 0.000000 0.000000  # GM" >> lte/path.dat
echo "0.000000 0.500000 0.000000  # L"  >> lte/path.dat

echo $temperature > lte/temperature.dat

caesar fourier_interpolation           \
       structure.dat                   \
       lte/phonon_dispersion_curve.dat \
       lte/high_symmetry_points.dat    \
       lte/temperature.dat             \
       lte/free_energy.dat             \
       lte/freq_dos.dat                \
       grid.dat                        \
       ibz.dat                         \
       lte/atoms_in_primitive_cell.    \
       lte/dyn_mat.                    \
       lte/path.dat                    \
       > lte/fourier_interpolation.out
