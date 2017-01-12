! A thin wrapper for running dft
module rundft_module
  implicit none
contains

subroutine rundft(args)
  use string_module
  implicit none
  
  character(100), intent(in) :: args(:)
  
  type(String) :: codename
  type(String) :: directory
  type(String) :: num_cores
  
  ! qe specific inputs
  type(String) :: seedname
  
  ! read inputs
  codename = args(1)
  directory = args(2)
  num_cores = args(3)
  
  call system('cd '//directory)
  
  if (codename=='castep') then
    call system('rundft nnodes '//num_cores)
    call system('rm *.castep_bin *.cst_esp *.usp machine_file *.bib *orbitals')
  elseif (codename=='vasp') then
  elseif (codename=='qe') then
    seedname = args(4)
    call system('mpirun -np '//num_cores//&
       & ' /rscratch/bm418/espresso-5.1.1/bin/pw.x -i '//&
       & seedname//'.in > '//seedname//'.out')
    call system('rm -r '//seedname//'.save')
  endif
  
  call system('cd -')
end subroutine
end module
