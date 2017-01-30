#!/bin/bash

# ======================================================================
# Script to set-up a generic harmonic calculation to use with LTE
# Also converts the generic harmonic calculation to castep vasp or qe
# ======================================================================

# Get settings from user
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
fi

echo $code > code.txt

# ----------------------------------------------------------------------
# Generate generic calculation
# ----------------------------------------------------------------------

# Add symmetries to structure.dat
caesar calculate_symmetry structure.dat

# Generate IBZ
# Reads Caesar input files. Writes ibz.dat and rotated_gvectors.dat
caesar generate_kgrid \
       structure.dat  \
       grid.dat       \
       ibz.dat        \
       rotated_gvectors.dat

# Generate non-diagonal supercells
# Makes Supercell_* directories
# Adds to Supercell_* directories:
#   supercell.dat
# Reads and writes ibz.dat
#   Adds sc_id to each kpoint
#   Orders kpoints by supercell
caesar generate_supercells     \
       structure.dat           \
       grid.dat                \
       ibz.dat                 \
       no_sc.dat               \
       Supercell_

no_sc=$(awk '{print}' no_sc.dat)

# Generate Force Constants 
for (( i=1; i<=$no_sc; i++ ))do
  sdir=Supercell_$i
  
  caesar construct_supercell     \
         structure.dat           \
         $sdir/supercell.dat     \
         $sdir/structure.dat
  
  # Add symmetries to structure.dat files
  caesar calculate_symmetry $sdir/structure.dat
  
  # Evaluate needed force constants
  caesar construct_matrix_force_cnsts \
         structure.dat                \
         $sdir/supercell.dat          \
         $sdir/structure.dat          \
         $sdir/force_constants.dat
  
  while read fline ; do
    line=($fline)
    atom=${line[0]}
    disp=${line[1]}
    ddir=$sdir/atom.$atom.disp.$disp
    mkdir $ddir
    mkdir $ddir/positive $ddir/negative
    caesar construct_finite_displacement \
           $atom                         \
           $disp                         \
           $sdir/structure.dat           \
           $ddir/positive/structure.dat  \
           $ddir/negative/structure.dat

  done < $sdir/force_constants.dat
  
  # ----------------------------------------------------------------------
  # Convert generic calculation to specific code
  # ----------------------------------------------------------------------
  
  echo "Converting supercell" $i

  while read fline ; do
    line=($fline)
    atom=${line[0]}
    disp=${line[1]}
    
    paths=(positive negative)
    for path in ${paths[@]}; do
      ddir=$sdir/atom.$atom.disp.$disp/$path
      if [ "$code" = "castep" ]; then
        caesar structure_to_dft     \
               $code                \
               $ddir/structure.dat  \
               $code/$seedname.cell \
               $ddir/$seedname.cell
      elif [ "$code" = "vasp" ]; then
        caesar structure_to_dft    \
               $code               \
               $ddir/structure.dat \
               $ddir/POSCAR
      elif [ "$code" = "qe" ]; then
        caesar structure_to_dft    \
               $code               \
               $ddir/structure.dat \
               $code/$seedname.in  \
               $code/pseudo.in     \
               $code/kpoints.in    \
               structure.dat       \
               $ddir/$seedname.in
      fi
    done
  done < $sdir/force_constants.dat
done
