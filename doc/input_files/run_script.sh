dft_code=$1
directory=$2
num_cores=$3
seedname=$4

cd $directory

if [ "$dft_code" = "castep" ]; then
  cp $dft_code/$seedname.param .
  rundft nnodes $num_cores seedname $seedname
  rm *.castep_bin *.cst_esp *.usp machine_file *.bib *orbitals
elif [ "$dft_code" = "vasp" ]; then
  # vasp run script not yet written
  echo "Error! vasp run script not yet written."
  exit
elif [ "$dft_code" = "qe" ]; then
  mpirun -np $num_cores /rscratch/bm418/espresso-5.1.1/bin/pw.x \
     -i $seedname.in > $seedname.out
  rm -r $seedname.save
fi

cd -
