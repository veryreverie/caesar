#!/bin/bash

# Takes a structure.dat file with no symmetry, and appends symmetry operators

# $1=structure.dat

# Input:

# Lattice
#  ~ ~ ~
#  ~ ~ ~
#  ~ ~ ~
# Atoms
#  ~ ~ ~
#   ...
#  ~ ~ ~
# End

# Output:

# Lattice
#  ~ ~ ~
#  ~ ~ ~
#  ~ ~ ~
# Atoms
#  ~ ~ ~
#   ...
#  ~ ~ ~
# Symmetry
#  ~ ~ ~
#   ...
#  ~ ~ ~
# End

caesar structure_to_dft \
       castep           \
       $1               \
       dummy_argument   \
       calculate_symmetry.cell

cellsym --symmetry calculate_symmetry.cell > calculate_symmetry.dat

caesar calculate_symmetry_helper \
       calculate_symmetry.dat    \
       $1
       
rm calculate_symmetry.cell
rm calculate_symmetry.dat
