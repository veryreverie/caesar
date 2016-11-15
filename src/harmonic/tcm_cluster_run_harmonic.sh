#!/bin/bash

# Script to run hamonic in TCM cluster

echo "What code do you want to use (castep,vasp,qe)?"
read code

if [ "$code" = "castep" ];then

  if [ ! -d "castep" ];then
    echo "Error! The directory 'castep' does not exist." 
    exit 1
  fi 

  echo "What is the first supercell to run?"
  read first_sc
  echo "What is the last supercell to run?"
  read last_sc
  echo "How many cores per run?"
  read num_cores

  # Loop over supercells
  for i in `seq $first_sc $last_sc`;
  do
 
    cd Supercell_$i

    lines=$( awk 'END {print NR}' force_constants.dat )

    # Loop over force constants
    for i in `seq 1 $lines`;
    do

      awk -v awk_line=$i 'NR==awk_line {print}' force_constants.dat > disp.dat
      atom=$(awk '{print $1}' disp.dat)
      disp=$(awk '{print $2}' disp.dat)

  
      cd atom.${atom}.disp.${disp}

      cd positive
      rundft nnodes $num_cores
      rm *.castep_bin *.cst_esp *.usp machine_file *.bands *.bib
      caesar fetch_forces_castep
      cd ../ 

      cd negative
      rundft nnodes $num_cores
      rm *.castep_bin *.cst_esp *.usp machine_file *.bands *.bib
      caesar fetch_forces_castep
      cd ../ 
      
      cd ../

    done 
    
    cd ../

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

  echo "What is the first supercell to run?"
  read first_sc
  echo "What is the last supercell to run?"
  read last_sc
  echo "How many cores per run?"
  read num_cores

  # Get seedname
  cd Supercell_$first_sc
    seedname=$( awk '{print}' seedname.txt )
  cd ../

  # Loop over supercells
  for i in `seq $first_sc $last_sc`;
  do

    cd Supercell_$i

    lines=$( awk 'END {print NR}' force_constants.dat )

    # Loop over force constants
    for i in `seq 1 $lines`;
    do

      awk -v awk_line=$i 'NR==awk_line {print}' force_constants.dat > disp.dat
      atom=$(awk '{print $1}' disp.dat)
      disp=$(awk '{print $2}' disp.dat)


      cd atom.${atom}.disp.${disp}

      cd positive
      pwd
      echo $seedname
      mpirun -np $num_cores /rscratch/bm418/espresso-5.1.1/bin/pw.x -i $seedname.in > $seedname.out
      rm -r $seedname.save
      caesar fetch_forces_qe
      cd ../

      cd negative
      pwd
      mpirun -np $num_cores /rscratch/bm418/espresso-5.1.1/bin/pw.x -i $seedname.in > $seedname.out
      rm -r $seedname.save
      caesar fetch_forces_qe
      cd ../

      cd ../

    done

    cd ../

  done
  echo "Done."

else 

  echo "Error! This code is not supported."

fi
