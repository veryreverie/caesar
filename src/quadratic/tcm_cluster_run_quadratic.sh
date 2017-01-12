#!/bin/bash

# Script to run quadratic in TCM cluster

harmonic_path=$( awk '{print $1}' harmonic_path.dat)

sampling_amplitude=$( awk 'NR==1 {print $1}' mapping.dat)
sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)
no_sampling_points=$(( $sampling_point_final-$sampling_point_init ))

atoms_line=$(awk -v IGNORECASE=1 '$1~/Atoms/{print NR}' \
   $harmonic_path/structure.dat)
symmetry_line=$(awk -v IGNORECASE=1 '$1~/Lattice/{print NR}' \
   $harmonic_path/structure.dat)
no_atoms=$($symmetry_line-$atoms_line-1)
no_modes=$(( $no_atoms*3 ))

echo "What code do you want to use (castep,vasp,qe)?"
read code

if [ "$code" = "castep" ];then

  if [ ! -d "castep" ];then
    echo "Error! The directory 'castep' does not exist." 
    exit 1
  fi 

  echo "What is the castep seedname?"
  read seedname
  echo "What is the first supercell to run?"
  read first_sc
  echo "What is the last supercell to run?"
  read last_sc
  echo "How many cores per run?"
  read num_cores


  # Loop over supercells
  for i in `seq $first_sc $last_sc`;
  do
    sdir=Supercell_$i

    no_atoms_sc=$( awk 'NR==1 {print $1}' \
       $harmonic_path/$sdir/super_equilibrium.dat )

    # Run static first
    caesar rundft castep $sdir $num_cores

    while read fline ; do
      line=($fline)
      big_point=${line[0]}
      kdir=$sdir/kpoint.$big_point/configurations
      
      for j in `seq 1 $no_modes`; do
        for k in `seq $sampling_point_init $sampling_point_final`; do
          if [ -e "$kdir/structure.$j.$k.dat" ]; then
            mdir=$kdir/mode.$j.$k
            mkdir $mdir
            cp $kdir/$seedname.param $mdir
            mv $kdir/$seedname.$j.$k.cell $mdir/$seedname.cell
            caesar rundft castep $mdir $num_cores
          fi
        done # loop over sampling points per mode
      done # loop over modes
    done < $harmonic_path/$sdir/list.dat

  done
  echo "Done."

elif [ "$code" = "vasp" ]; then

  if [ ! -d "vasp" ];then
    echo "Error! The directory 'vasp' does not exist." 
    exit 1
  fi 

elif [ "$code" = "qe" ]; then

  if [ ! -d "qe" ];then
    echo "Error! The directory 'qe' does not exist." 
    exit 1
  fi 

  echo "What is the qe seedname?"
  read seedname
  seedname_nscf=${seedname}.nscf
  echo "What is the first supercell to run?"
  read first_sc
  echo "What is the last supercell to run?"
  read last_sc 
  echo "How many cores per run?"
  read num_cores

  # Loop over supercells
  for i in `seq $first_sc $last_sc`;
  do
    sdir=Supercell_$i

    no_atoms_sc=$( awk 'NR==1 {print $1}' $sdir/super_equilibrium.dat )

    # Run static first
    caesar rundft qe $sdir/static $num_cores $seedname
    if [ -e "$sdir/$seedname_nscf.in" ]; then
      caesar rundft qe $sdir/static $num_cores $seedname_nscf
    fi
    
    while read line ; do
      line=($fline)
      big_point=${line[0]}
      echo "Big point" $big_point
      kdir=$sdir/kpoint.$big_point/configurations
      for j in `seq 1 $no_modes`; do
        for k in `seq $sampling_point_init $sampling_point_final`; do
          if [ -e "$kdir/structure.$j.$k.dat" ]; then
            mdir=$kdir/mode/$j/$k
            mkdir $mdir
            mv $kdir/$seedname.$j.$k.in $mdir/$seedname.in
            if [ -e "$kdir/$seedname_nscf.$j.$k.in" ];then
              mv $kdir/$seedname_nscf.$j.$k.in $mdir
            fi
            cp $kdir/*UPF $mdir
            caesar rundft qe $mdir $num_cores $seedname
            if [ -e "$mdir/$seedname_nscf.in" ]; then
              caesar rundft qe $mdir $num_cores $seedname_nscf
            fi
          fi
        done # loop over sampling points per mode
      done # loop over modes
    done < $harmonic_path/$sdir/list.dat
  done
  echo "Done."
  
else 

  echo "Error! This code is not supported."

fi
