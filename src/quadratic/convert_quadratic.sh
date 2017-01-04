#!/bin/bash

# Script to convert a generic quadratic calculation to:
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
  echo $seedname > seedname.txt
  
  if [ ! -f "castep/$seedname.param" ];then
    echo "Error! The castep 'param' file does not exist." 
    exit 1
  fi

  
  # Get primitive cell path 
  bs_calculation=0
  if [ -e "castep/path.dat" ];then
    bs_calculation=1
    bs_path_init=$(awk -v IGNORECASE=1 '/block bs_kpoint_path/{print NR; exit}' castep/path.dat )
    bs_path_final=$(awk -v IGNORECASE=1 '/endblock bs_kpoint_path/{print NR}' castep/path.dat )
    echo $(( $bs_path_final-$bs_path_init-1 )) > castep/bs_path.dat
    awk "NR=$(($bs_path_init + 1)),NR==$(($bs_path_final - 1)) " \
      '{print $1 " " $2 " " $3}' castep/path.dat >> castep/bs_path.dat
  fi
  
  no_sc=$(ls -1d Supercell_* | wc -l)
  sampling_amplitude=$( awk 'NR==1 {print $1}' mapping.dat)
  sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
  sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)
  no_sampling_points=$(( $sampling_point_final-$sampling_point_init ))
  
  # Loop over supercells
  for (( i=1; i<=$no_sc; i++ )) do
    sdir=Supercell_$i
    no_atoms=$( awk 'NR==1 {print $1}' $sdir/equilibrium.dat )
    no_modes=$(( $no_atoms*3 ))
    no_atoms_sc=$( awk 'NR==1 {print $1}' $sdir/super_equilibrium.dat )
  
    echo "Converting supercell" $i
    
    static_dir=$sdir/static
    cp castep/* $static_dir
    if [ "$bs_calculateion" -eq "1" ]; then
      cp $sdir/supercell.dat $static_dir
      caesar generate_sc_path          \
             $static_dir/supercell.dat \
             $static_dir/bs_path.dat   \
             $static_dir/sc_bs_path.dat
    fi
    caesar structure_to_castep        \
           $static_dir/structure.dat  \
           $static_dir/sc_bs_path.dat \
           $static_dir/$seedname.cell
    
    while read fline ; do
      line=($fline)
      big_point=${line[0]}
      kdir=$sdir/kpoint.$big_point/configurations
      cp castep/* $kdir
      if [ "$bs_calculateion" -eq "1" ]; then
        cp $static_dir/sc_bs_path.dat $kdir
      fi
      for j in `seq 1 $no_modes`; do
        for k in `seq $sampling_point_init $sampling_point_final`; do
          if [-e "$kdir/structure.$j.$k.dat" ]; then
            cp $kdir/structure.$j.$k.dat $kdir/structure.dat
            cp $kdir/$seedname.cell $kdir/$seedname.$j.$k.cell
            caesar structure_to_castep  \
                   $kdir/structure.dat  \
                   $kdir/sc_bs_path.dat \
                   $kdir/$seedname.$j.$k.cell
          fi
        done
      done
      
      rm $kdir/$seedname.cell
    done < $sdir/list.dat
  done
  echo "Done."

elif [ "$code" = "vasp" ]; then

  if [ ! -d "vasp" ];then
    echo "Error! The directory 'vasp' does not exist." 
    exit 1
  fi 

  sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
  sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)

  
  # Loop over supercells
  no_sc=$(ls -1d Supercell_* | wc -l)
  kpoint_counter=1
  for (( i=1; i<=$no_sc; i++ )) do
    sdir=Supercell_$i

    echo "Converting supercell" $i

    no_atoms=$( awk 'NR==1 {print $1}' $sdir/equilibrium.dat )
    no_modes=$(( $no_atoms*3 ))

    static_dir=$sdir/static
    cp vasp/INCAR* $static_dir
    cp vasp/POTCAR $static_dir
    cp vasp/KPOINTS.band $static_dir
    cp $sdir/KPOINTS.$i $static_dir/KPOINTS
    cp vasp/mpi_submit_static.sh $static_dir
    caesar structure_to_vasp         \
           $static_dir/structure.dat \
           $static_dir/POSCAR
    
    # Loop over k-points
    no_kpoints=$(ls -1d kpoint.* | wc -l)
    for (( j=$kpoint_counter; j<=$(( $kpoint_counter+($no_kpoints-1) )); j++ )) do
      kdir=$sdir/kpoint.$j/configurations
      cp $sdir/kpoint.$j/mapping.dat $sdir/equilibrium.dat $kdir
      cp vasp/INCAR* $kdir
      cp vasp/POTCAR $kdir
      cp vasp/KPOINTS.band $kdir
      cp vasp/mpi_submit_quadratic.sh $kdir
      cp $sdir/KPOINTS.$i $kdir/KPOINTS
      
      # Loop over number of modes
      for (( k=1; k<=$no_modes; k++ )) do
        for l in `seq $sampling_point_init $sampling_point_final`; do
          if [ -e "$kdir/structure.$k.$l.dat" ]; then
            cp $kdir/structure.$k.$l.dat $kdir/structure.dat
            caesar structure_to_vasp $kdir/structure.dat $kdir/POSCAR.$k.$l
          fi
        done
      done
    done # Loop over k-points
    
    kpoint_counter=$(( $kpoint_counter+$no_kpoints ))

  done # Loop over supercells

elif [ "$code" = "qe" ];then

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
  seedname_nscf=${seedname}.nscf
  
  echo $seedname > seedname.txt
  echo $seedname_nscf > seedname.nscf.txt

  no_sc=$(ls -1d Supercell_* | wc -l)
  sampling_amplitude=$( awk 'NR==1 {print $1}' mapping.dat)
  sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
  sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)
  no_sampling_points=$(( $sampling_point_final-$sampling_point_init ))

  # Loop over 
  for (( i=1; i<=$no_sc; i++ )) do
    sdir=Supercell_$i
    
    no_atoms=$( awk 'NR==1 {print $1}' $sdir/equilibrium.dat )
    no_modes=$(( $no_atoms*3 ))
    no_atoms_sc=$( awk 'NR==1 {print $1}' $sdir/super_equilibrium.dat )

    echo "Converting supercell" $i

    # Generate supercell k-point mesh
    cp qe/kpoints.in $sdir
    cp qe/kpoints.nscf.in $sdir
    caesar generate_supercell_kpoint_mesh_qe \
           $sdir/kpoints.in                  \
           $sdir/lattice.dat                 \
           $sdir/super_lattice.dat           \
           $sdir/sc_kpoints.dat
    header=$(awk 'NR==1,NR==1 {print}' $sdir/kpoints.in)
    echo $header > $sdir/kpoints.in
    cat $sdir/sc_kpoints.dat >> $sdir/kpoints.in
    header=$(awk 'NR==1,NR==1 {print}' $sdir/kpoints.nscf.in)
    echo $header > $sdir/kpoints.nscf.in
    cat $sdir/sc_kpoints.nscf.dat >> $sdir/kpoints.nscf.in

    static_dir=$sdir/static
    cp qe/* $static_dir
    caesar structure_to_qe           \
           $static_dir/structure.dat \
           $static_dir/pseudo.in     \
           $sdir/kpoints.in          \
           $static_dir/$seedname.in
    caesar structure_to_qe           \
           $static_dir/structure.dat \
           $static_dir/pseudo.in     \
           $sdir/kpoints.nscf.in     \
           $static_dir/$seedname_nscf.in
    if [ -f "$static_dir/pseudo.in" ]; then
      rm $static_dir/pseudo.in
    fi
    
    while read fline ; do
      line=($fline)
      big_point=${line[0]}
      kdir=$sdir/kpoint.$big_point/configurations
      
      cp qe/* $kdir
      
      cp $sdir/kpoints.in $kdir
      for j in `seq 1 $no_modes`; do
        for k in `seq $sampling_point_init $sampling_point_final`; do
          if [ -e "$kdir/structure.${j}.${k}.dat" ];then
            cp $kdir/structure.${j}.${k}.dat $kdir/structure.dat
            if [ -f "$kdir/$seedname.in" ]; then
              cp $kdir/$seedname.in $kdir/$seedname.$j.$k.in
            fi
            caesar structure_to_qe     \
                   $kdir/structure.dat \
                   $kdir/pseudo.in     \
                   $kdir/kpoints.in    \
                   $kdir/$seedname.$j.$k.in
          fi # structure exists
        done # loop over sampling points per mode
      done # loop over modes
      rm $kdir/$seedname.in
      
      cp $sdir/kpoints.nscf.in $kdir/kpoints.in
      if [ -f "$kdir/$seedname_nscf.in" ]; then
        for j in `seq 1 $no_modes`; do
          for k in `seq $sampling_point_init $sampling_point_final`; do
            if [ -e "$kdir/structure.${j}.${k}.dat" ];then
              cp $kdir/structure.${j}.${k}.dat $kdir/structure.dat
              cp $kdir/$seedname_nscf.in $kdir/$seedname_nscf.$j.$k.in
              caesar structure_to_qe     \
                     $kdir/structure.dat \
                     $kdir/pseudo.in     \
                     $kdir/kpoints.in    \
                     $kdir/$seedname_nscf.$j.$k.in
            fi # structure exists
          done # loop over sampling points per mode
        done # loop over modes
        rm $kdir/$seedname_nscf.in
      fi

    done < $sdir/list.dat

  done
  echo "Done."

else 

  echo "Error! This code is not supported."

fi
