! ======================================================================
! Runs unit tests on setup_harmonic.
! ======================================================================
module setup_harmonic_test_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

subroutine setup_harmonic_test(arguments)
  use utils_module, only : format_path, make_dft_output_filename, mkdir
  use structure_module
  use unique_directions_module
  use dft_input_file_module
  use dft_output_file_module
  use group_module
  use atom_mapping_module
  use dictionary_module
  use setup_harmonic_module
  use structure_test_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Previous settings.
  type(String)              :: dft_code
  type(String)              :: seedname
  
  type(String), allocatable :: no_supercells_file(:)
  integer                   :: no_supercells
  
  ! Directory and file names.
  type(String) :: new_wd
  type(String) :: old_wd
  type(String) :: dft_input_filename
  type(String) :: sdir_old
  type(String) :: sdir_new
  type(String) :: old_dirs(2)
  type(String) :: new_dirs(2)
  type(String) :: dft_output_filename
  
  ! Structure information.
  type(StructureData) :: structure
  type(StructureData) :: structure_new
  type(StructureData) :: structure_old
  type(Group)         :: new_to_old
  
  ! Grid information.
  integer :: grid_file
  integer :: grid(3)
  
  ! Unique direction information.
  type(UniqueDirections) :: unique_directions
  integer                :: atom
  character(1)           :: direction
  
  ! DFT output files.
  type(DftOutputFile) :: dft_output
  integer             :: copied_output
  
  ! Temporary variables.
  integer :: i,j,k,l,m
  
  ! ----------------------------------------------------------------------
  ! Read in settings.
  ! ----------------------------------------------------------------------
  new_wd = item(arguments, 'working_directory')
  dft_code = item(arguments, 'dft_code')
  seedname = item(arguments, 'seedname')
  grid = int(split(item(arguments, 'q-point_grid')))
  
  ! ----------------------------------------------------------------------
  ! Run new setup_harmonic.
  ! ----------------------------------------------------------------------
  call setup_harmonic(arguments)
  
  ! Read number of supercells.
  no_supercells_file = read_lines(new_wd//'/no_sc.dat')
  no_supercells = int(no_supercells_file(1))
  
  ! Read in structure file.
  structure = read_structure_file(new_wd//'/structure.dat')
  
  ! ----------------------------------------------------------------------
  ! Run old setup_harmonic.
  ! ----------------------------------------------------------------------
  
  ! Make a directory for running the old code in.
  old_wd = new_wd//'/old_caesar'
  call mkdir(old_wd)
  call mkdir(old_wd//'/'//dft_code)
  
  ! Copy DFT input files.
  dft_input_filename = make_dft_input_filename(dft_code,seedname)
  call reduce_dft_input_file( dft_code,                        &
                            & new_wd//'/'//dft_input_filename, &
                            & old_wd//'/'//dft_code//'/'//dft_input_filename)
  
  ! Write structure.dat and grid.dat.
  call write_old_structure_file(structure, old_wd//'/structure.dat')
  
  grid_file = open_write_file(old_wd//'/grid.dat')
  call print_line(grid_file, grid)
  close(grid_file)
  
  ! Call setup_harmonic.
  call system_call('cd '//old_wd//'; setup_harmonic.sh')
  
  ! ----------------------------------------------------------------------
  ! Compare results.
  ! ----------------------------------------------------------------------
  do i=1,no_supercells
    call print_line('')
    call print_line('Supercell '//i//':')
    sdir_old = old_wd//'/Supercell_'//i
    sdir_new = new_wd//'/Supercell_'//i
    
    ! Read in both structure files.
    structure_old = read_structure_file(sdir_old//'/structure.dat')
    structure_new = read_structure_file(sdir_new//'/structure.dat')
    
    new_to_old = atom_mapping(structure_new,structure_old)
    
    call print_line('Atom positions correct.')
    call print_line('Atoms in new calculation map to those in &
       &old calculation as:')
    call print_line(new_to_old%operation)
    
!    unique_directions = &
!       & read_unique_directions_file(sdir_new//'/unique_directions.dat')
!    do j=1,size(unique_directions)
!      atom = unique_directions%atoms(j)
!      direction = unique_directions%directions_char(j)
!      
!      old_dirs = [ sdir_old//'/atom.'//atom//'.disp.'//k//'/positive', &
!                  & sdir_old//'/atom.'//atom//'.disp.'//k//'/negative'  ]
!      
!      new_dirs = [ sdir_new//'/atom.'//atom//'.+d'//direction, &
!                 & sdir_new//'/atom.'//atom//'.-d'//direction  ]
!      do l=1,2
!        ! Read in the old castep file.
!        dft_output_filename = make_dft_output_filename(dft_code,seedname)
!        dft_output_filename = old_dirs(l)//'/'//dft_output_filename
!        dft_output = read_dft_output_file(dft_code,dft_output_filename)
!      
!        ! Make a fake castep output file.
!        copied_output = open_write_file( &
!           & new_dirs(l)//'/'//seedname//'.castep')
!        call print_line( copied_output,&
!                       & 'final energy, E = '//dft_output%energy)
!        call print_line(copied_output,'*********************** forces')
!        do m=1,5
!          call print_line(copied_output,'')
!        enddo
!        do m=1,structure_old%no_atoms
!          call print_line(copied_output, &
!             & ': '//structure_old%species(operate(new_to_old,m))//' : '// &
!             & dft_output%forces(:,operate(new_to_old,m)))
!        enddo
!        call print_line(copied_output,'')
!        call print_line(copied_output,'*************************************&
!                        &*****************')
!        close(copied_output)
!      enddo
!    enddo
  enddo
end subroutine
end module
