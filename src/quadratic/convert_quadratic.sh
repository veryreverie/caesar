#!/bin/bash

# Script to convert a generic quadratic calculation to:
# CASTEP
# VASP 

harmonic_path=$( awk '{print $1}' harmonic_path.dat)

no_sc=$( awk '{print $1}' $harmonic_path/no_sc.dat)
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
  echo $seedname > seedname.txt
  
  if [ ! -f "castep/$seedname.param" ];then
    echo "Error! The castep 'param' file does not exist." 
    exit 1
  fi
  
  # Get primitive cell path 
  if [ -e "castep/path.dat" ];then
    bs_path_init=$(awk -v IGNORECASE=1 '/block bs_kpoint_path/{print NR; exit}' castep/path.dat )
    bs_path_final=$(awk -v IGNORECASE=1 '/endblock bs_kpoint_path/{print NR}' castep/path.dat )
    echo $(( $bs_path_final-$bs_path_init-1 )) > castep/bs_path.dat
    awk "NR=$(($bs_path_init + 1)),NR==$(($bs_path_final - 1)) " \
      '{print $1 " " $2 " " $3}' castep/path.dat >> castep/bs_path.dat
  fi
  
  
  # Loop over supercells
  for (( i=1; i<=$no_sc; i++ )) do
    sdir=Supercell_$i
    no_atoms_sc=$( awk 'NR==1 {print $1}' \
       $harmonic_path/$sdir/super_equilibrium.dat )
  
    echo "Converting supercell" $i
    
    static_dir=$sdir/static
    cp castep/* $static_dir
    if [ -e "castep/path.dat" ];then
      caesar generate_sc_path                   \
             $harmonic_path/$sdir/supercell.dat \
             $static_dir/bs_path.dat            \
             $static_dir/sc_bs_path.dat
      caesar structure_to_dft           \
             castep                     \
             $static_dir/structure.dat  \
             $static_dir/sc_bs_path.dat \
             $static_dir/$seedname.cell
    else
      caesar structure_to_dft           \
             castep                     \
             $static_dir/structure.dat  \
             $static_dir/$seedname.cell
    fi
    
    while read fline ; do
      line=($fline)
      big_point=${line[0]}
      kdir=$sdir/kpoint.$big_point/configurations
      cp castep/* $kdir
      if [ -e "castep/path.dat" ];then
        cp $static_dir/sc_bs_path.dat $kdir
      fi
      for j in `seq 1 $no_modes`; do
        for k in `seq $sampling_point_init $sampling_point_final`; do
          if [-e "$kdir/structure.$j.$k.dat" ]; then
            cp $kdir/structure.$j.$k.dat $kdir/structure.dat
            cp $kdir/$seedname.cell $kdir/$seedname.$j.$k.cell
            caesar structure_to_dft     \
                   castep               \
                   $kdir/structure.dat  \
                   $kdir/sc_bs_path.dat \
                   $kdir/$seedname.$j.$k.cell
          fi
        done
      done
      
      rm $kdir/$seedname.cell
    done < $harmonic_path/$sdir/list.dat
  done
  echo "Done."

elif [ "$code" = "vasp" ]; then

  if [ ! -d "vasp" ];then
    echo "Error! The directory 'vasp' does not exist." 
    exit 1
  fi 
  
  # Loop over supercells
  kpoint_counter=1
  for i in $(seq 1 $no_sc); do
    sdir=Supercell_$i

    echo "Converting supercell" $i

    static_dir=$sdir/static
    cp vasp/INCAR* $static_dir
    cp vasp/POTCAR $static_dir
    cp vasp/KPOINTS.band $static_dir
    cp $harmonic_path/$sdir/KPOINTS.$i $static_dir/KPOINTS
    cp vasp/mpi_submit_static.sh $static_dir
    caesar structure_to_dft          \
           vasp                      \
           $static_dir/structure.dat \
           $static_dir/POSCAR
    
    # Loop over k-points
    no_kpoints=$(ls -1d kpoint.* | wc -l)
    for j in $(seq $kpoint_counter $(( $kpoint_counter+($no_kpoints-1) ))); do
      kdir=$sdir/kpoint.$j/configurations
      cp vasp/INCAR* $kdir
      cp vasp/POTCAR $kdir
      cp vasp/KPOINTS.band $kdir
      cp vasp/mpi_submit_quadratic.sh $kdir
      cp $harmonic_path/$sdir/KPOINTS.$i $kdir/KPOINTS
      
      # Loop over number of modes
      for (( k=1; k<=$no_modes; k++ )) do
        for l in `seq $sampling_point_init $sampling_point_final`; do
          if [ -e "$kdir/structure.$k.$l.dat" ]; then
            cp $kdir/structure.$k.$l.dat $kdir/structure.dat
            caesar structure_to_dft vasp $kdir/structure.dat $kdir/POSCAR.$k.$l
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

  # Loop over 
  for (( i=1; i<=$no_sc; i++ )) do
    sdir=Supercell_$i
    
    no_atoms_sc=$( awk 'NR==1 {print $1}' \
       $harmonic_path/$sdir/super_equilibrium.dat )

    echo "Converting supercell" $i

    # Generate supercell k-point mesh
    caesar generate_supercell_kpoint_mesh_qe      \
           qe/kpoints.in                          \
           $harmonic_path/structure.dat           \
           $harmonic_path/$sdir/super_lattice.dat \
           $sdir/kpoints.in
    
    caesar generate_supercell_kpoint_mesh_qe      \
           qe/kpoints.nscf.in                     \
           $harmonic_path/structure.dat           \
           $harmonic_path/$sdir/super_lattice.dat \
           $sdir/kpoints.nscf.in

    static_dir=$sdir/static
    cp qe/* $static_dir
    caesar structure_to_dft          \
           $code                     \
           $static_dir/structure.dat \
           $static_dir/pseudo.in     \
           $sdir/kpoints.in          \
           $static_dir/$seedname.in
    caesar structure_to_dft          \
           $code                     \
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
            caesar structure_to_dft    \
                   $code               \
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
              caesar structure_to_dft    \
                     $code               \
                     $kdir/structure.dat \
                     $kdir/pseudo.in     \
                     $kdir/kpoints.in    \
                     $kdir/$seedname_nscf.$j.$k.in
            fi # structure exists
          done # loop over sampling points per mode
        done # loop over modes
        rm $kdir/$seedname_nscf.in
      fi

    done < $harmonic_path/$sdir/list.dat

  done
  echo "Done."

else 

  echo "Error! This code is not supported."

fi
