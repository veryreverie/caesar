#!/bin/bash
echo "Are you sure (y/n)?"
read answer

if [ $answer == 'y' ]; then
  rm -r Supercell_* ibz.dat kpoint_to_supercell.dat no_sc.dat rotated_gvectors.dat
fi
