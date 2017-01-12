#!/bin/bash

# Script to convert a generic harmonic calculation to:
# CASTEP
# VASP 
  
no_sc=$(awk '{print}' no_sc.dat )

echo "What code do you want to use (castep,vasp,qe)?"
read code

if [ "$code" = "castep" ];then

  if [ ! -d "castep" ];then
    echo "Error! The directory 'castep' does not exist." 
    exit 1
  fi 
 
  echo "What is the castep seedname?"
  read seedname
  
  if [ ! -f "castep/$seedname.param" ];then
    echo "Error! The castep 'param' file does not exist." 
    exit 1
  fi
  
  echo $seedname > seedname.txt
  
  
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
        dir=$sdir/atom.$atom.disp.$disp/$path
        cp castep/* $dir
        caesar structure_to_castep \
               $dir/structure.dat  \
               dummy_argument      \
               $dir/seedname.cell
      done
    done < $sdir/force_constants.dat
      
  done
  echo "Done."

elif [ "$code" = "vasp" ]; then

  if [ ! -d "vasp" ];then
    echo "Error! The directory 'vasp' does not exist." 
    exit 1
  fi 

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
        dir=$sdir/atom.$atom.disp.$disp/$path
        cp vasp/* $dir
        caesar structure_to_vasp $dir/structure.dat $dir/POSCAR
      done
    done < $sdir/force_constants.dat

  done
  echo "Done."

elif [ "$code" = "qe" ]; then

  if [ ! -d "qe" ];then
    echo "Error! The directory 'qe' does not exist." 
    exit 1
  fi
 
  echo "What is the qe seedname?"
  read seedname

  if [ ! -f "qe/$seedname.in" ];then
    echo "Error! The qe 'in' file does not exist." 
    exit 1
  fi

  echo $seedname > seedname.txt

  # Loop over 
  for (( i=1; i<=$no_sc; i++ )) do
    
    sdir=Supercell_$i

    # Generate supercell k-point mesh
    cp qe/kpoints.in $sdir
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
        dir=$sdir/atom.$atom.disp.$disp/$path
        
        cp qe/* $dir
        caesar structure_to_qe    \
               $dir/structure.dat \
               $dir/pseudo.in     \
               $sdir/kpoints.in   \
               $dir/seedname.in
        if [ -f "$dir/pseudo.in" ]; then
          rm $dir/pseudo.in 
        fi
        
      done
      
    done < $sdir/force_constants.dat

  done
  echo "Done."

else 

  echo "Error! This code is not supported."

fi
