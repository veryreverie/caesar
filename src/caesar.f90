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
  
  ! use misc modules
  use calculate_gap_module
  use hartree_to_eV_module
  
  use constants, only : dp
  use linear_algebra
  
  implicit none
  
  ! Command line arguments
  type(String), allocatable :: args(:)
  
  integer  :: i
  real(dp) :: a(3,3)
  real(dp) :: b(3,3)
  real(dp) :: c(3,3)
  real(dp) :: d(3,3)
  
  ! read in command line arguments
  args = command_line_args()
  
  if (size(args) == 1) then
    call print_line('No arguments given. For help, call caesar -h')
    stop
  endif
  
  if (args(2) == '-h' .or. args(2) == '--help') then
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
  ! Test
  elseif (args(2) == 'test') then
    a(1,:) = (/ 1,0,0 /)
    a(2,:) = (/ 0,1,0 /)
    a(3,:) = (/ 0,0,1 /)
    
    b(1,:) = (/ 1,0,0 /)
    b(2,:) = (/ 0,1,0 /)
    b(3,:) = (/ 0,0,1 /)
    
    c = invert(a)
    d = matmul(b,c)
    
    call print_line('')
    do i=1,3
      call print_line(a(i,:))
    enddo
    
    call print_line('')
    call print_line(determinant(a))
    
    call print_line('')
    do i=1,3
      call print_line(c(i,:))
    enddo
    
    call print_line('')
    do i=1,3
      call print_line(d(i,:))
    enddo
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
  ! Wrappers for testing modules
  elseif (args(2) == 'test_copy_harmonic') then
    call test_copy_harmonic()
  ! unrecognised argument
  else
    call print_line(char('Unrecognised argument : '//args(2)))
  endif
end program
