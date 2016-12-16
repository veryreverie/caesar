#!/bin/bash

# Script to convert a generic quadratic calculation to:
# CASTEP
# VASP 

echo "What code do you want to use (castep,vasp,qe)?"
read code

if [ "$code" = "castep" ];then

  if [ ! -d "castep" ];then
    echo "Error! The directory 'castep' does not exist." 
    exit 1
  fi 
 
  echo "What is the castep seedname?"
  read seedname
  
  cd castep
  if [ ! -f "$seedname.param" ];then
    echo "Error! The castep 'param' file does not exist." 
    exit 1
  fi
  cd ../

  
  # Get primitive cell path 
  bs_calculation=0
  cd castep
    if [ -e "path.dat" ];then
      bs_calculation=1
      bs_path_init=$(awk -v IGNORECASE=1 '/block bs_kpoint_path/{print NR; exit}' path.dat )
      bs_path_final=$(awk -v IGNORECASE=1 '/endblock bs_kpoint_path/{print NR}' path.dat )
      echo $(( $bs_path_final-$bs_path_init-1 )) > bs_path.dat
      awk -v awk_bs_path_init=$bs_path_init -v awk_bs_path_final=$bs_path_final 'NR==(awk_bs_path_init+1), NR==(awk_bs_path_final-1) {print $1 " " $2 " " $3 }' path.dat >> bs_path.dat
    fi
  cd ../
  
  no_sc=$(ls -1d Supercell_* | wc -l)
  sampling_amplitude=$( awk 'NR==1 {print $1}' mapping.dat)
  sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
  sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)
  no_sampling_points=$(( $sampling_point_final-$sampling_point_init ))
  
  # Loop over 
  for (( i=1; i<=$no_sc; i++ )) do
 
    cd Supercell_$i
    echo $seedname > seedname.txt
    no_atoms=$( awk 'NR==1 {print $1}' equilibrium.dat )
    no_modes=$(( $no_atoms*3 ))
    no_atoms_sc=$( awk 'NR==1 {print $1}' super_equilibrium.dat )
  
    echo "Converting supercell" $i
 
    cd static
      cp ../../castep/* .
      if [ "$bs_calculation" -eq "1" ];then
        cp ../supercell.dat .
        caesar generate_sc_path
      fi
      mv ${seedname}.cell bottom.cell 
      caesar structure_to_castep .
      mv structure.cell ${seedname}.cell
      rm bottom.cell
    cd ../  

    while read line ; do
      big_point=$(echo ${line} | awk '{print $1}')
      cd kpoint.$big_point/configurations
      cp ../../../castep/* .
      mv ${seedname}.cell bottom.cell
      if [ "$bs_calculation" -eq "1" ];then
        cp ../../static/sc_bs_path.dat .
      fi
      for j in `seq 1 $no_modes`; do
        for k in `seq $sampling_point_init $sampling_point_final`; do
          if [ -e "structure.${j}.${k}.dat" ];then 
            cp structure.${j}.${k}.dat structure.dat
            caesar structure_to_castep .
            mv structure.cell ${seedname}.${j}.${k}.cell
          fi # structure exists
        done # loop over sampling points per mode
      done # loop over modes

      rm bottom.cell
 
      cd ../../

    done < list.dat

    cd ../

  done
  echo "Done."

elif [ "$code" = "vasp" ]; then

  if [ ! -d "vasp" ];then
    echo "Error! The directory 'vasp' does not exist." 
    exit 1
  fi 

  sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
  sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)

  
  # Loop over supercells
  no_sc=$(ls -1d Supercell_* | wc -l)
  kpoint_counter=1
  for (( i=1; i<=$no_sc; i++ )) do

    cd Supercell_$i

    echo "Converting supercell" $i

    no_atoms=$( awk 'NR==1 {print $1}' equilibrium.dat )
    no_modes=$(( $no_atoms*3 ))

    cd static
      cp ../../vasp/INCAR* .
      cp ../../vasp/POTCAR .
      cp ../../vasp/KPOINTS.band .
      cp ../KPOINTS.${i} KPOINTS
      cp ../../vasp/mpi_submit_static.sh .
      caesar structure_to_vasp .
      mv structure.POSCAR POSCAR
    cd ../
   
    # Loop over k-points
    no_kpoints=$(ls -1d kpoint.* | wc -l)
    for (( j=$kpoint_counter; j<=$(( $kpoint_counter+($no_kpoints-1) )); j++ )) do

      cd kpoint.$j
      cd configurations
      cp ../mapping.dat ../../equilibrium.dat .
      cp ../../../vasp/INCAR* .
      cp ../../../vasp/POTCAR .
      cp ../../../vasp/KPOINTS.band .
      cp ../../../vasp/mpi_submit_quadratic.sh .
      cp ../../KPOINTS.${i} KPOINTS
     
      # Loop over number of modes
      for (( k=1; k<=$no_modes; k++ )) do
 
        for l in `seq $sampling_point_init $sampling_point_final`; do

          if [ -e "structure.${k}.${l}.dat" ];then
            cp structure.${k}.${l}.dat structure.dat
            caesar structure_to_vasp .
            mv structure.POSCAR POSCAR.${k}.${l}
          fi
        
        done # Loop over sampling points

      done # Loop over modes
      
      cd ../
      cd ../

    done # Loop over k-points
    kpoint_counter=$(( $kpoint_counter+$no_kpoints ))
    
    cd ../

  done # Loop over supercells

elif [ "$code" = "qe" ];then

  if [ ! -d "qe" ];then
    echo "Error! The directory 'qe' does not exist." 
    exit 1
  fi

  echo "What is the qe seedname?" 
  read seedname

  cd qe
  if [ ! -f "$seedname.in" ];then
    echo "Error! The qe 'in' file does not exist." 
    exit 1
  fi
  cd ../
  seedname_nscf=${seedname}.nscf

  no_sc=$(ls -1d Supercell_* | wc -l)
  sampling_amplitude=$( awk 'NR==1 {print $1}' mapping.dat)
  sampling_point_init=$( awk 'NR==2 {print $1}' mapping.dat)
  sampling_point_final=$( awk 'NR==2 {print $2}' mapping.dat)
  no_sampling_points=$(( $sampling_point_final-$sampling_point_init ))

  # Loop over 
  for (( i=1; i<=$no_sc; i++ )) do

    cd Supercell_$i
    echo $seedname > seedname.txt
    echo $seedname_nscf > seedname.nscf.txt
    no_atoms=$( awk 'NR==1 {print $1}' equilibrium.dat )
    no_modes=$(( $no_atoms*3 ))
    no_atoms_sc=$( awk 'NR==1 {print $1}' super_equilibrium.dat )

    echo "Converting supercell" $i

    # Generate supercell k-point mesh
    cp ../qe/kpoints.in .
    cp ../qe/kpoints.nscf.in .
    caesar generate_supercell_kpoint_mesh_qe \
           kpoints.in                        \
           lattice.dat                       \
           super_lattice.dat                 \
           sc_kpoints.dat
    caesar generate_supercell_kpoint_mesh_qe_nscf # TODO : bug!
    awk 'NR==1,NR==1 {print}' kpoints.in > kpoints.in.temp
    cat kpoints.in.temp sc_kpoints.dat > kpoints.in.temp2
    mv kpoints.in.temp2 kpoints.in
    awk 'NR==1,NR==1 {print}' kpoints.nscf.in > kpoints.nscf.in.temp
    cat kpoints.nscf.in.temp sc_kpoints.nscf.dat > kpoints.nscf.in.temp2
    mv kpoints.nscf.in.temp2 kpoints.nscf.in

    cd static
      cp ../../qe/* .
      cp ../kpoints.in .
      if [ -f "$seedname.in" ]; then
        mv $seedname.in top.in
      fi
      caesar structure_to_qe .
      mv structure.in $seedname.in
      cp ../kpoints.nscf.in kpoints.in
      if [ -f "$seedname_nscf.in" ]; then
        mv $seedname_nscf.in top.in
      fi
      caesar structure_to_qe .
      mv structure.in $seedname_nscf.in
      if [ -f 'top.in' ]; then
        rm top.in
      fi
      if [ -f 'kpoints.in' ]; then
        rm kpoints.in
      fi
      if [ -f 'pseudo.in' ]; then
        rm pseudo.in
      fi
    cd ../

    while read line ; do
      big_point=$(echo ${line} | awk '{print $1}')
      cd kpoint.$big_point/configurations
      cp ../../../qe/* .
      cp ../../kpoints.in .
      if [ -f "$seedname.in" ]; then
        mv $seedname.in top.in
      fi
      for j in `seq 1 $no_modes`; do
        for k in `seq $sampling_point_init $sampling_point_final`; do
          if [ -e "structure.${j}.${k}.dat" ];then
            cp structure.${j}.${k}.dat structure.dat
            caesar structure_to_qe .
            mv structure.in ${seedname}.${j}.${k}.in
          fi # structure exists
        done # loop over sampling points per mode
      done # loop over modes
      cp ../../kpoints.nscf.in kpoints.in
      if [ -f "$seedname_nscf.in" ]; then
        mv $seedname_nscf.in top.in
        for j in `seq 1 $no_modes`; do
          for k in `seq $sampling_point_init $sampling_point_final`; do
            if [ -e "structure.${j}.${k}.dat" ];then
              cp structure.${j}.${k}.dat structure.dat
              caesar structure_to_qe .
              mv structure.in ${seedname_nscf}.${j}.${k}.in
            fi # structure exists
          done # loop over sampling points per mode
        done # loop over modes
      fi

      cd ../../

    done < list.dat

    cd ../

  done
  echo "Done."


     




else 

  echo "Error! This code is not supported."

fi
