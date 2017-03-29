program caesar
  ! use utility modules
  use utils, only : command_line_args
  
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
  
  ! use testing modules
  use test_copy_harmonic_module
  use test_lte_module
  use test_copy_quadratic_module
  
  ! use misc modules
  use calculate_gap_module
  use hartree_to_eV_module
  use err_module
  
  implicit none
  
  ! Command line arguments.
  type(String), allocatable :: args(:)
  type(String)              :: source_dir
  type(String)              :: cwd
  integer                   :: terminal_width
  type(String)              :: mode
  
  ! Working directory.
  type(String) :: wd
  
  ! read in command line arguments.
  args = command_line_args()
  
  if (size(args) < 3) then
    call print_line('No arguments given. For help, call caesar -h')
    call err()
  endif
  
  source_dir = args(1)
  cwd = args(2)
  terminal_width = int(args(3))
  mode = args(4)
  
  ! Set terminal width for formatting purposes.
  STRING_MODULE_TERMINAL_WIDTH = terminal_width
  
  ! For now, the working directory is always the current working directory.
  ! TODO: add options to change this.
  wd = cwd
  
  if (mode == '-h' .or. mode == '--help') then
    call print_line('caesar [-h] [option]')
    call print_line('')
    call print_line('-h :')
    call print_line('  Displays this help text')
    call print_line('')
    call print_line('option : utilities :')
    call print_line('  hartree_to_eV :')
    call print_line('    Provides a Hartree to eV calculator')
    call print_line('')
    call print_line('option : harmonic calculations :')
    call print_line('  setup_harmonic :')
    call print_line('    Sets up calculation')
    call print_line('    Converts calculation to specific DFT code')
    call print_line('    DFT code choices are castep, vasp and qe')
    call print_line('  run_harmonic :')
    call print_line('    Runs calculation on the TCM cluster')
    call print_line('    Should be called after setup_harmonic')
    call print_line('  lte_harmonic :')
    call print_line('    Runs harmonic calculations')
    call print_line('    Should be called after run_harmonic')
    call print_line('  clear_all :')
    call print_line('    Deletes all temporary files and folders')
    call print_line('')
    call print_line('option : quadratic calculations :')
    call print_line('  setup_quadratic :')
    call print_line('    Sets up quadratic calculation for use with a DFT code')
    call print_line('    DFT code choices are castep, vasp and qe')
    call print_line('    Should be called after lte_harmonic')
    call print_line('  run_quadratic :')
    call print_line('    Runs calculation on the TCM cluster')
    call print_line('    Should be called after setup_quadratic')
    call print_line('  anharmonic :')
    call print_line('    Runs anharmonic calculations')
    call print_line('    Should be called after run_quadratic')
    call print_line('  bs_quadratic :')
    call print_line('    Runs band structure calculations')
    call print_line('    Should be called after run_quadratic')
    call print_line('  get_kpoints :')
    call print_line('    [Help text pending]')
    call print_line('  calculate_gap')
    call print_line('    [Help text pending]')
  elseif (mode == 'test') then
    call print_line('And in the case of the extremely long line, &
       &which does rather insist on going well beyond the end of &
       &the terminal, it can often be helpful to insert a tactical &
       &line-break or two. Of course, these must not be allowed to break up &
       &individual words if at all possible, lest the resulting &
       &output take on a rather disjointed appearance. It should also be &
       &noted that in some cases you just have to give up. &
       &supercalifragilisticexpialidocious-floccinaucinihilipilification-antidisestablishmentarianism.&
       & (Probably in cases like that.)')
  ! Wrappers for top-level Fortran.
  elseif (mode == 'setup_harmonic') then
    call setup_harmonic(wd,source_dir)
  elseif (mode == 'lte_harmonic') then
    call lte_harmonic(wd)
  elseif (mode == 'setup_quadratic') then
    call setup_quadratic(wd,cwd)
  elseif (mode == 'anharmonic') then
    call anharmonic(wd)
  elseif (mode == 'bs_quadratic') then
    call bs_quadratic(wd)
  ! Wrappers for subsidiary Fortran 
  elseif (mode == 'calculate_gap') then
    call calculate_gap()
  elseif (mode == 'hartree_to_eV') then
    call hartree_to_eV()
  elseif (mode == 'structure_to_dft') then
    call structure_to_dft( dft_code              = args(5), &
                         & structure_sc_filename = args(6), &
                         & output_filename       = args(7))
  elseif (mode == 'calculate_symmetry_helper') then
    call calculate_symmetry_helper(args(5),args(6))
  ! Wrappers for testing modules
  elseif (mode == 'test_copy_harmonic') then
    call test_copy_harmonic(wd,cwd)
  elseif (mode == 'test_lte') then
    call test_lte(wd,cwd)
  elseif (mode == 'test_copy_quadratic') then
    call test_copy_quadratic(wd,cwd)
  ! unrecognised argument
  else
    call print_line(char('Unrecognised argument : '//mode))
  endif
end program
