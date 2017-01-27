#!/bin/bash
echo "Are you sure (y/n)?"
read answer

if [ $answer == 'y' ]; then
  rm -r Supercell_* ibz.dat no_sc.dat rotated_gvectors.dat
fi
