#!/bin/bash

# Script to set-up a generic quadratic calculation 
# Converts calculation to castep, vasp or qe

# Read in harmonic path
echo "What is the path to the harmonic directory?"
read harmonic_path

echo "What code do you want to use (castep,vasp,qe)?"
read code

if [ "$code" = "vasp" ]; then
  echo "Error! vasp is not currently supported."
  exit 1
elif [ ! "$code" = "castep"] && [ ! "$code" = "qe"]; then
  echo "Error! The code $code is not supported."
  echo "Please choose one of: castep vasp qe"
  exit 1
fi

if [ ! -d "$code" ]; then
  echo "Error! The directory '$code' does not exist." 
  exit 1
fi 

if [ "$code" = "castep" ] || [ "$code" = "qe" ]; then
  if [ "$code" = "castep" ] && [ ! -f "$code/$seedname.param" ];then
    echo "Error! The $code input file (.param) does not exist." 
    exit 1
  else if [ "$code" = "qe" ] && [ ! -f "$code/$seedname.in" ];then
    echo "Error! The $code input file (.in) does not exist." 
    exit 1
  fi
  
  echo "What is the $code seedname?"
  read seedname
  echo $seedname > seedname.txt
  seedname_nscf=$seedname.nscf
fi

# Write harmonic_path.dat
echo $harmonic_path > harmonic_path.dat

# Read in number of supercells
no_sc=$( awk '{print $1}' no_sc.dat)

# Define some useful variables
sampling_amplitude=$( awk 'NR==1 {print $1}' mapping.dat)
sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)
no_sampling_points=$(( $sampling_point_final-$sampling_point_init ))
temperature=$( awk 'NR==3 {print $1}' mapping.dat)

atoms_line=$(awk -v IGNORECASE=1 '$1~/Atoms/{print NR}' \
   $harmonic_path/structure.dat)
symmetry_line=$(awk -v IGNORECASE=1 '$1~/Lattice/{print NR}' \
   $harmonic_path/structure.dat)
no_atoms=$($symmetry_line-$atoms_line-1)
no_modes=$(( $no_atoms*3 ))

# ----------------------------------------------------------------------
# Set up generic calculation
# ----------------------------------------------------------------------
# Loop over supercells
for i in `seq 1 $no_sc`; do
  sdir=Supercell_$i
  mkdir $sdir
  
  no_atoms_sc=$( awk 'NR==1 {print $1}' \
     $harmonic_path/$sdir/super_equilibrium.dat )
  
  # Generate configurations

  # Generate static calculation
  mkdir $sdir/static
  
  echo Lattice                                   > $sdir/static/structure.dat
  cat $harmonic_path/$sdir/super_lattice.dat    >> $sdir/static/structure.dat
  echo Atoms                                    >> $sdir/static/structure.dat
  awk "NR==2,NR==$(( $no_atoms_sc+1 )) {print}" \
     $harmonic_path/$sdir/super_equilibrium.dat >> $sdir/static/structure.dat
  echo Symmetry                                 >> $sdir/static/structure.dat
  echo End                                      >> $sdir/static/structure.dat

  # loop over kpoints
  while read fline ; do
    line=$(fline)
    big_point=${line[0]}
    small_point=${line[1]}
    
    first_line=$(( (4+$no_atoms_sc)*($no_modes*($small_point-1))+1 ))
    last_line=$(( $first_line+(4+$no_atoms_sc)*$no_modes ))
    
    kdir=kpoint.$big_point
    mkdir $kdir
    awk "NR==$first_line,NR==$last_line {print}" \
       $harmonic_path/$sdir/lte/disp_patterns.dat > $kdir/disp_patterns.dat
    mkdir $kdir/configurations
    for j in `seq 1 $no_modes`; do
      frequency_line=$(( (4+$no_atoms_sc)*($j-1)+1 ))
      
      # write frequency in a.u. to frequency.$big_point.$j.dat
      frequency=$(awk "NR=$frequency_line {print\$3} $kdir/disp_patterns.dat")
      python -c "print($frequency*27.211396132)" > $kdir/frequency.$j.dat
      
      # write kpoint.$big_point/configurations/mode.j.k/structure.dat
      for k in `seq $sampling_point_init $sampling_point_final`; do
        mdir=$kdir/mode.$j.$k
        mkdir $mdir
        caesar generate_quadratic_configurations          \
               $sampling_amplitude                        \
               $k                                         \
               $no_sampling_points                        \
               $frequency                                 \
               $frequency_line                            \
               $harmonic_path/$sdir/super_equilibrium.dat \
               $harmonic_path/$sdir/super_lattice.dat     \
               $kdir/disp_patterns.dat                    \
               $mdir/structure.dat
      done # loop over sampling points per mode
    done # loop over modes
  done < $harmonic_path/$sdir/list.dat
done # Loop over supercells

# ----------------------------------------------------------------------
# Convert to specific code
# ----------------------------------------------------------------------
if [ "$code" = "castep" ] && [ -e "$code/path.dat" ]; then
  # Get primitive cell path 
  bs_path_init=$(awk -v IGNORECASE=1 '/block bs_kpoint_path/{print NR; exit}' $code/path.dat )
  bs_path_final=$(awk -v IGNORECASE=1 '/endblock bs_kpoint_path/{print NR}' $code/path.dat )
  echo $(( $bs_path_final-$bs_path_init-1 )) > $code/bs_path.dat
  awk "NR=$(($bs_path_init + 1)),NR==$(($bs_path_final - 1)) " \
    '{print $1 " " $2 " " $3}' $code/path.dat >> $code/bs_path.dat
fi
  
# Loop over supercells
for (( i=1; i<=$no_sc; i++ )) do
  sdir=Supercell_$i
  no_atoms_sc=$( awk 'NR==1 {print $1}' \
     $harmonic_path/$sdir/super_equilibrium.dat )
  echo "Converting supercell" $i
  
  # Run static calculation
  static_dir=$sdir/static
  if [ "$code" = "castep" ];then
    cp $code/* $static_dir
    if [ -e "$code/path.dat" ];then
      caesar generate_sc_path                   \
             $harmonic_path/$sdir/supercell.dat \
             $static_dir/bs_path.dat            \
             $static_dir/sc_bs_path.dat
      caesar structure_to_dft           \
             $code                      \
             $static_dir/structure.dat  \
             $static_dir/sc_bs_path.dat \
             $static_dir/$seedname.cell
    else
      caesar structure_to_dft           \
             $code                      \
             $static_dir/structure.dat  \
             $static_dir/$seedname.cell
    fi
  elif [ "$code" = "vasp" ]; then
    cp $code/INCAR* $static_dir
    cp $code/POTCAR $static_dir
    cp $code/KPOINTS.band $static_dir
    cp $code/mpi_submit_static.sh $static_dir
    cp $harmonic_path/$sdir/KPOINTS.$i $static_dir/KPOINTS
    caesar structure_to_dft          \
           $code                     \
           $static_dir/structure.dat \
           $static_dir/POSCAR
  elif [ "$code" = "qe" ]; then
    # Generate supercell k-point mesh
    caesar generate_supercell_kpoint_mesh_qe      \
           $code/kpoints.in                       \
           $harmonic_path/structure.dat           \
           $harmonic_path/$sdir/super_lattice.dat \
           $sdir/kpoints.in
    
    caesar generate_supercell_kpoint_mesh_qe      \
           $code/kpoints.nscf.in                  \
           $harmonic_path/structure.dat           \
           $harmonic_path/$sdir/super_lattice.dat \
           $sdir/kpoints.nscf.in

    cp $code/* $static_dir
    caesar structure_to_dft          \
           $code                     \
           $static_dir/structure.dat \
           $code/pseudo.in           \
           $sdir/kpoints.in          \
           $static_dir/$seedname.in
    caesar structure_to_dft          \
           $code                     \
           $static_dir/structure.dat \
           $code/pseudo.in           \
           $sdir/kpoints.nscf.in     \
           $static_dir/$seedname_nscf.in
    if [ -f "$static_dir/pseudo.in" ]; then
      rm $static_dir/pseudo.in
    fi
  fi
    
  # Loop over kpoints
  while read fline ; do
    line=($fline)
    big_point=${line[0]}
    kdir=$sdir/kpoint.$big_point/configurations
    
    for j in `seq 1 $no_modes`; do
      for k in `seq $sampling_point_init $sampling_point_final`; do
        mdir=$kdir/mode.$j.$k
        if [-e "$mdir/structure.dat" ]; then
          if [ "$code" = "castep" ];then
            cp $code/seedname.param $mdir
            cp $code/$seedname.cell $mdir
            caesar structure_to_dft           \
                   $code                      \
                   $mdir/structure.dat        \
                   $static_dir/sc_bs_path.dat \
                   $mdir/$seedname.cell
          elif [ "$code" = "vasp" ]; then
            cp $code/INCAR* $mdir
            cp $code/POTCAR $mdir
            cp $code/KPOINTS.band $mdir
            cp $code/mpi_submit_quadratic.sh $mdir
            cp $harmonic_path/$sdir/KPOINTS.$i $mdir/KPOINTS
            caesar structure_to_dft    \
                   $code               \
                   $mdir/structure.dat \
                   $mdir/POSCAR.$j.$k
          elif [ "$code" = "qe" ];then
            cp $code/*UPF $mdir
            cp $code/$seedname.in $mdir
            caesar structure_to_dft    \
                   $code               \
                   $mdir/structure.dat \
                   $code/pseudo.in     \
                   $code/kpoints.in    \
                   $mdir/$seedname.in
            if [ -f "$code/$seedname_nscf.in" ]; then
              cp $code/$seedname_nscf.in $mdir
              caesar structure_to_dft      \
                     $code                 \
                     $mdir/structure.dat   \
                     $code/pseudo.in       \
                     $code/kpoints.nscf.in \
                     $mdir/$seedname_nscf.in
            fi
          fi
        fi # structure exists
      done # loop over sampling points per mode
    done # loop over modes
  done < $harmonic_path/$sdir/list.dat
done
