#!/bin/bash

# Script to convert a generic harmonic calculation to:
# CASTEP
# VASP 
  
no_sc=$(awk '{print}' no_sc.dat )

echo "What code do you want to use (castep,vasp,qe)?"
read code

if [ ! "$code"="castep" ] && [ ! "$code"="vasp" ] && [ ! "$code"="qe" ]; then
  echo "Error! The code $code is not supported."
  echo "Please choose one of: castep vasp qe"
  exit 1
fi

if [ ! -d "$code" ]; then
  echo "Error! The directory '$code' does not exist." 
  exit 1
fi 

if [ "$code" = "castep" ] || [ "$code" = "qe" ]; then
  echo "What is the $code seedname?"
  read seedname
  if [ ! -f "$code/$seedname.param" ];then
    echo "Error! The $code input file (.param or .in) does not exist." 
    exit 1
  fi
  echo $seedname > seedname.txt
fi


if [ "$code" = "castep" ];then
  # Loop over 
  for (( i=1; i<=$no_sc; i++ )) do
    sdir=Supercell_$i
    echo "Converting supercell" $i
  
    while read fline ; do
      line=($fline)
      atom=${line[0]}
      disp=${line[1]}
      
      paths=(positive negative)
      for path in ${paths[@]}; do
        ddir=$sdir/atom.$atom.disp.$disp/$path
        cp $code/* $dir
        caesar structure_to_dft    \
               $code               \
               $ddir/structure.dat \
               $ddir/seedname.cell
      done
    done < $sdir/force_constants.dat
  done
  echo "Done."
elif [ "$code" = "vasp" ]; then
  for (( i=1; i<=$no_sc; i++ )) do
    sdir=Supercell_$i
    echo "Converting supercell" $i

    while read fline ; do
      line=($fline)
      atom=${line[0]}
      disp=${line[1]}
      
      paths=(positive negative)
      for path in ${paths[@]}; do
        ddir=$sdir/atom.$atom.disp.$disp/$path
        cp $code/* $ddir
        caesar structure_to_dft $code $ddir/structure.dat $ddir/POSCAR
      done
    done < $sdir/force_constants.dat

  done
  echo "Done."
elif [ "$code" = "qe" ]; then
  for (( i=1; i<=$no_sc; i++ )) do
    sdir=Supercell_$i
    # Generate supercell k-point mesh
    cp $code/kpoints.in $sdir
    caesar generate_supercell_kpoint_mesh_qe \
           $sdir/kpoints.in                  \
           structure.dat                     \
           $sdir/super_lattice.dat           \
           $sdir/sc_kpoints.dat
    header=$(awk 'NR==1,NR==1 {print}' $sdir/kpoints.in)
    echo $header > $sdir/kpoints.in
    cat $sdir/sc_kpoints.dat >> $sdir/kpoints.in

    echo "Converting supercell" $i

    while read fline ; do
      line=($fline)
      atom=${line[0]}
      disp=${line[1]}
      
      paths=(positive negative)
      for path in ${paths[@]}; do
        ddir=$sdir/atom.$atom.disp.$disp/$path
        cp $code/* $ddir
        caesar structure_to_qe    \
               $ddir/structure.dat \
               $ddir/pseudo.in     \
               $sdir/kpoints.in   \
               $ddir/seedname.in
        if [ -f "$ddir/pseudo.in" ]; then
          rm $ddir/pseudo.in 
        fi
        
      done
      
    done < $sdir/force_constants.dat

  done
  echo "Done."

else 

  echo "Error! This code is not supported."

fi
