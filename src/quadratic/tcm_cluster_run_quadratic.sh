#!/bin/bash

# Script to run quadratic in TCM cluster

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

  sampling_amplitude=$( awk 'NR==1 {print $1}' mapping.dat)
  sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
  sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)
  no_sampling_points=$(( $sampling_point_final-$sampling_point_init ))


  # Loop over supercells
  for i in `seq $first_sc $last_sc`;
  do
 
    cd Supercell_$i

    no_atoms=$( awk 'NR==1 {print $1}' equilibrium.dat )
    no_modes=$(( $no_atoms*3 ))
    no_atoms_sc=$( awk 'NR==1 {print $1}' super_equilibrium.dat )


    # Run static first
    cd static
      rundft nnodes $num_cores
      rm *.castep_bin *.cst_esp *.usp machine_file *.bib *.orbitals
    cd ../

    while read fline ; do
      line=($fline)
      big_point=${line[0]}
      cd kpoint.$big_point/configurations

      for j in `seq 1 $no_modes`; do
        for k in `seq $sampling_point_init $sampling_point_final`; do
          if [ -e "structure.${j}.${k}.dat" ]; then
            mkdir mode.${j}.${k}
            mv $seedname.${j}.${k}.cell mode.${j}.${k}
            cp $seedname.param mode.${j}.${k}
            cd mode.${j}.${k}
              mv $seedname.${j}.${k}.cell $seedname.cell
              rundft nnodes $num_cores 
              rm *.castep_bin *.cst_esp *.usp machine_file *.bib *.orbitals
            cd ../
          fi
        done # loop over sampling points per mode
      done # loop over modes

      cd ../../  # cd kpoint.$big_point/configurations

    done < list.dat

    cd ../  # cd Supercell_$i

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

  sampling_amplitude=$( awk 'NR==1 {print $1}' mapping.dat)
  sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
  sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)
  no_sampling_points=$(( $sampling_point_final-$sampling_point_init ))

  # Loop over supercells
  for i in `seq $first_sc $last_sc`;
  do

    cd Supercell_$i

    no_atoms=$( awk 'NR==1 {print $1}' equilibrium.dat )
    no_modes=$(( $no_atoms*3 ))
    no_atoms_sc=$( awk 'NR==1 {print $1}' super_equilibrium.dat )


    # Run static first
    cd static
      mpirun -np $num_cores /rscratch/bm418/espresso-5.1.1/bin/pw.x -i $seedname.in > $seedname.out
      if [ -e "$seedname_nscf.in" ];then
        mv $seedname_nscf.in $seedname.in 
        mpirun -np $num_cores /rscratch/bm418/espresso-5.1.1/bin/pw.x -i $seedname.in > $seedname_nscf.out
      fi
      rm -r $seedname.save
    cd ../

    while read -u 3 line ; do

      big_point=$(echo ${line} | awk '{print $1}')
      echo "Big point" $big_point
      cd kpoint.$big_point/configurations

      echo "In big point" $big_point

      for j in `seq 1 $no_modes`; do
        for k in `seq $sampling_point_init $sampling_point_final`; do
          if [ -e "structure.${j}.${k}.dat" ]; then
            mkdir mode.${j}.${k}
            mv $seedname.${j}.${k}.in mode.${j}.${k}
            if [ -e "${seedname_nscf}.${j}.${k}.in" ];then
              mv $seedname_nscf.${j}.${k}.in mode.${j}.${k}
            fi
            cp *UPF mode.${j}.${k}
            cd mode.${j}.${k}
              mv $seedname.${j}.${k}.in $seedname.in
              mpirun -np $num_cores /rscratch/bm418/espresso-5.1.1/bin/pw.x -i $seedname.in > $seedname.out
              if [ -e "$seedname_nscf.in" ];then
                mv $seedname_nscf.${j}.${k}.in $seedname.in 
                mpirun -np $num_cores /rscratch/bm418/espresso-5.1.1/bin/pw.x -i $seedname.in > $seedname_nscf.out
              fi
              rm -r $seedname.save
            cd ../
          fi
        done # loop over sampling points per mode
      done # loop over modes

      cd ../../  # cd kpoint.$big_point/configurations

    done 3< list.dat


    cd ../  # cd Supercell_$i

  done
  echo "Done."






else 

  echo "Error! This code is not supported."

fi
