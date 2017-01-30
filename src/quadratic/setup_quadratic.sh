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
for sc_id in `seq 1 $no_sc`; do
  sdir=Supercell_$sc_id
  mkdir $sdir
  mkdir $sdir/static
done

# loop over kpoints
while read fline ; do
  line=($fline)
  kpoint=${line[0]}
  gvector=${line[1]}
  sc_id=${line[2]}
  
  sdir=Supercell_$sc_id
  
  atoms_line_sc=$(awk -v IGNORECASE=1 '$1~/Atoms/{print NR}' \
     $harmonic_path/$sdir/structure.dat)
  symmetry_line_sc=$(awk -v IGNORECASE=1 '$1~/Lattice/{print NR}' \
     $harmonic_path/$sdir/structure.dat)
  no_atoms_sc=$($symmetry_line_sc-$atoms_line_sc-1)
  
  kdir=kpoint.$kpoint
  mkdir $kdir
  
  for j in `seq 1 $no_modes`; do
    # write kpoint.$kpoint/mode.j.k/structure.dat
    for k in `seq $sampling_point_init $sampling_point_final`; do
      mdir=$kdir/mode.$j.$k
      mkdir $mdir
      caesar generate_quadratic_configurations          \
             $sampling_amplitude                        \
             $gvector                                   \
             $j                                         \
             $k                                         \
             $no_sampling_points                        \
             $harmonic_path/structure.dat               \
             $harmonic_path/$sdir/structure.dat         \
             $harmonic_path/$sdir/lte/disp_patterns.dat \
             $mdir/structure.dat
    done # loop over sampling points per mode
  done # loop over modes
done < $harmonic_path/list.dat

# ----------------------------------------------------------------------
# Convert to specific code
# ----------------------------------------------------------------------

# Loop over supercells
for sc_id in `seq 1 $no_sc`; do
  sdir=Supercell_$sc_id
  # Run static calculation
  static_dir=$sdir/static
  if [ "$code" = "castep" ];then
    caesar structure_to_dft                   \
           $code                              \
           $harmonic_path/$sdir/structure.dat \
           $code/$seedname.cell               \
           $harmonic_path/$sdir/supercell.dat \
           $code/path.dat                     \
           $static_dir/$seedname.cell
  elif [ "$code" = "vasp" ]; then
    caesar structure_to_dft                    \
           $code                               \
           $harmonic_path/$sdir/structure.dat  \
           $static_dir/POSCAR
  elif [ "$code" = "qe" ]; then
    caesar structure_to_dft                   \
           $code                              \
           $harmonic_path/$sdir/structure.dat \
           $code/$seedname.in                 \
           $code/pseudo.in                    \
           $code/kpoints.in                   \
           $harmonic_path/structure.dat       \
           $static_dir/$seedname.in
    caesar structure_to_dft                   \
           $code                              \
           $harmonic_path/$sdir/structure.dat \
           $code/$seedname_nscf.in            \
           $code/pseudo.in                    \
           $code/kpoints.nscf.in              \
           $harmonic_path/structure.dat       \
           $static_dir/$seedname_nscf.in
  fi
done

# Loop over kpoints
while read fline ; do
  line=($fline)
  kpoint=${line[0]}
  kdir=kpoint.$kpoint
  for j in `seq 1 $no_modes`; do
    for k in `seq $sampling_point_init $sampling_point_final`; do
      mdir=$kdir/mode.$j.$k
      if [-e "$mdir/structure.dat" ]; then
        if [ "$code" = "castep" ];then
          caesar structure_to_dft                   \
                 $code                              \
                 $mdir/structure.dat                \
                 $code/$seedname.cell               \
                 $harmonic_path/$sdir/supercell.dat \
                 $code/path.dat                     \
                 $mdir/$seedname.cell
        elif [ "$code" = "vasp" ]; then
          caesar structure_to_dft    \
                 $code               \
                 $mdir/structure.dat \
                 $mdir/POSCAR.$j.$k
        elif [ "$code" = "qe" ];then
          caesar structure_to_dft             \
                 $code                        \
                 $mdir/structure.dat          \
                 $code/$seedname.in           \
                 $code/pseudo.in              \
                 $code/kpoints.in             \
                 $harmonic_path/structure.dat \
                 $mdir/$seedname.in
          if [ -f "$code/$seedname_nscf.in" ]; then
            caesar structure_to_dft             \
                   $code                        \
                   $mdir/structure.dat          \
                   $code/$seedname.in           \
                   $code/pseudo.in              \
                   $code/kpoints.nscf.in        \
                   $harmonic_path/structure.dat \
                   $mdir/$seedname_nscf.in
          fi
        fi
      fi # structure exists
    done # loop over sampling points per mode
  done # loop over modes
done < $harmonic_path/list.dat
