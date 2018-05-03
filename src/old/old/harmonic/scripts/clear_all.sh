echo "Are you sure (y/n)?"
read answer

if [ $answer == 'y' ]; then

  rm -r Supercell_* equilibrium.dat ibz.dat kpoint_to_supercell.dat lattice.dat no_sc.dat rotated_gvectors.dat symmetry.dat

fi
