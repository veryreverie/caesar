#!/bin/bash

# Script to convert a generic harmonic calculation to:
# CASTEP
# VASP 

echo "What code do you want to use (castep,vasp,qe)?"
read code

if [ "$code" = "castep" ];then

  if [ ! -d "castep" ];then
    echo "Error! The directory 'castep' does not exist." 
    exit 1
  fi 
 
  echo "What is the castep seedname?"
  read seedname
  
  cd castep
  if [ ! -f "$seedname.param" ];then
    echo "Error! The castep 'param' file does not exist." 
    exit 1
  fi
  cd ../
  
  no_sc=$(awk '{print}' no_sc.dat )
  
  # Loop over 
  for (( i=1; i<=$no_sc; i++ )) do
 
    cd Supercell_$i
    echo $seedname > seedname.txt
  
 
    echo "Converting supercell" $i
  
    while read LINE ; do
      echo $LINE > disp.dat
      atom=$(awk '{print $1}' disp.dat)
      disp=$(awk '{print $2}' disp.dat)
      cd atom.${atom}.disp.${disp}
      cd positive
      cp ../../../castep/* .
      cp ../../seedname.txt .
      if [ -f "$seedname.cell" ]; then
        mv $seedname.cell bottom.cell
      fi
      caesar structure_to_castep
      mv structure.cell $seedname.cell
      if [ -f 'bottom.cell' ]; then
        rm bottom.cell
      fi
      cd ../
      cd negative
      cp ../../../castep/* .
      cp ../../seedname.txt .
      if [ -f "$seedname.cell" ]; then
        mv $seedname.cell bottom.cell
      fi
      caesar structure_to_castep
      mv structure.cell $seedname.cell
      if [ -f 'bottom.cell' ]; then
        rm bottom.cell
      fi
      cd ../
      cd ../
    done < force_constants.dat
    
    cd ../

  done
  echo "Done."

elif [ "$code" = "vasp" ]; then

  if [ ! -d "vasp" ];then
    echo "Error! The directory 'vasp' does not exist." 
    exit 1
  fi 

  no_sc=$(awk '{print}' no_sc.dat )

  # Loop over 
  for (( i=1; i<=$no_sc; i++ )) do

    cd Supercell_$i

    echo "Converting supercell" $i

    while read LINE ; do
      echo $LINE > disp.dat
      atom=$(awk '{print $1}' disp.dat)
      disp=$(awk '{print $2}' disp.dat)
      cd atom.${atom}.disp.${disp}
      cd positive
      cp ../../../vasp/* .
      caesar structure_to_vasp
      mv structure.POSCAR POSCAR
      cd ../
      cd negative
      cp ../../../vasp/* .
      caesar structure_to_vasp
      mv structure.POSCAR POSCAR
      cd ../
      cd ../
    done < force_constants.dat

    cd ../

  done
  echo "Done."

elif [ "$code" = "qe" ]; then

  if [ ! -d "qe" ];then
    echo "Error! The directory 'qe' does not exist." 
    exit 1
  fi
 
  echo "What is the qe seedname?"
  read seedname

  cd qe
  if [ ! -f "$seedname.in" ];then
    echo "Error! The qe 'in' file does not exist." 
    exit 1
  fi
  cd ../

  no_sc=$(awk '{print}' no_sc.dat )

  # Loop over 
  for (( i=1; i<=$no_sc; i++ )) do

    cd Supercell_$i
    echo $seedname > seedname.txt

    # Generate supercell k-point mesh
    cp ../qe/kpoints.in .
    caesar generate_supercell_kpoint_mesh_qe
    awk 'NR==1,NR==1 {print}' kpoints.in > kpoints.in.temp
    cat kpoints.in.temp sc_kpoints.dat > kpoints.in.temp2
    mv kpoints.in.temp2 kpoints.in
 

    echo "Converting supercell" $i

    while read LINE ; do
      echo $LINE > disp.dat
      atom=$(awk '{print $1}' disp.dat)
      disp=$(awk '{print $2}' disp.dat)
      cd atom.${atom}.disp.${disp}
      cd positive
      cp ../../../qe/* .
      cp ../../seedname.txt .
      cp ../../kpoints.in .
      if [ -f "$seedname.in" ]; then
        mv $seedname.in top.in
      fi
      caesar structure_to_qe
      mv structure.in $seedname.in
      if [ -f 'top.in' ]; then
        rm top.in
      fi
      if [ -f 'kpoints.in' ]; then
        rm kpoints.in
      fi
      if [ -f 'pseudo.in' ]; then
        rm pseudo.in 
      fi
      cd ../
      cd negative
      cp ../../../qe/* .
      cp ../../seedname.txt .
      cp ../../kpoints.in .
      if [ -f "$seedname.in" ]; then
        mv $seedname.in top.in
      fi
      caesar structure_to_qe
      mv structure.in $seedname.in
      if [ -f 'top.in' ]; then
        rm top.in
      fi
      if [ -f 'kpoints.in' ]; then
        rm kpoints.in 
      fi
      if [ -f 'pseudo.in' ]; then
        rm pseudo.in 
      fi
      cd ../
      cd ../
    done < force_constants.dat

    cd ../

  done
  echo "Done."






else 

  echo "Error! This code is not supported."

fi
