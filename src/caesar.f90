program caesar
  ! use utility modules
  use utils,       only : command_line_args
  use file_module, only : open_write_file, open_read_file
  
  ! use class modules
  use string_module
  
  ! use common modules
  use structure_to_dft_module
  use calculate_symmetry_helper_module
  
  ! use harmonic modules
  use setup_harmonic_module
  use lte_harmonic_module
  
  ! use quadratic modules
  use setup_quadratic_module
  use anharmonic_module
  use bs_quadratic_module
  
  ! use misc modules
  use calculate_gap_module
  use hartree_to_eV_module
  
  implicit none
  
  ! Command line arguments
  type(String), allocatable :: args(:)
  
  ! read in command line arguments
  args = command_line_args()
  
  if (size(args) == 1) then
    write(*,*) 'No arguments given. For help, call caesar -h'
    stop
  endif
  
  if (args(2) == '-h' .or. args(2) == '--help') then
    write(*,*) 'caesar [-h] [option]'
    write(*,*) ''
    write(*,*) '-h :'
    write(*,*) '  Displays this help text'
    write(*,*) ''
    write(*,*) 'option : utilities :'
    write(*,*) '  hartree_to_eV :'
    write(*,*) '    Provides a Hartree to eV calculator'
    write(*,*) ''
    write(*,*) 'option : harmonic calculations :'
    write(*,*) '  setup_harmonic :'
    write(*,*) '    Sets up calculation'
    write(*,*) '    Converts calculation to specific DFT code'
    write(*,*) '    DFT code choices are castep, vasp and qe'
    write(*,*) '  run_harmonic :'
    write(*,*) '    Runs calculation on the TCM cluster'
    write(*,*) '    Should be called after setup_harmonic'
    write(*,*) '  lte_harmonic :'
    write(*,*) '    Runs harmonic calculations'
    write(*,*) '    Should be called after run_harmonic'
    write(*,*) '  clear_all :'
    write(*,*) '    Deletes all temporary files and folders'
    write(*,*) ''
    write(*,*) 'option : quadratic calculations :'
    write(*,*) '  setup_quadratic :'
    write(*,*) '    Sets up quadratic calculation for use with a DFT code'
    write(*,*) '    DFT code choices are castep, vasp and qe'
    write(*,*) '    Should be called after lte_harmonic'
    write(*,*) '  run_quadratic :'
    write(*,*) '    Runs calculation on the TCM cluster'
    write(*,*) '    Should be called after setup_quadratic'
    write(*,*) '  anharmonic :'
    write(*,*) '    Runs anharmonic calculations'
    write(*,*) '    Should be called after run_quadratic'
    write(*,*) '  bs_quadratic :'
    write(*,*) '    Runs band structure calculations'
    write(*,*) '    Should be called after run_quadratic'
    write(*,*) '  get_kpoints :'
    write(*,*) '    [Help text pending]'
    write(*,*) '  calculate_gap'
    write(*,*) '    [Help text pending]'
  ! Wrappers for top-level Fortran
  elseif (args(2) == 'setup_harmonic') then
    call setup_harmonic(args(1))
  elseif (args(2) == 'lte_harmonic') then
    call lte_harmonic()
  elseif (args(2) == 'setup_quadratic') then
    call setup_quadratic()
  elseif (args(2) == 'anharmonic') then
    call anharmonic()
  elseif (args(2) == 'bs_quadratic') then
    call bs_quadratic()
  ! Wrappers for subsidiary Fortran 
  elseif (args(2) == 'calculate_gap') then
    call calculate_gap()
  elseif (args(2) == 'hartree_to_eV') then
    call hartree_to_eV()
  elseif (args(2) == 'structure_to_dft') then
    call structure_to_dft( dft_code=args(3),              &
                         & structure_sc_filename=args(4), &
                         & output_filename=args(5))
  elseif (args(2) == 'calculate_symmetry_helper') then
    call calculate_symmetry_helper(args(3),args(4))
  ! unrecognised argument
  else
    write(*,*) char('Unrecognised argument : '//args(2))
  endif
end program
