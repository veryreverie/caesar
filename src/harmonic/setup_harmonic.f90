! ======================================================================
! The first stage of Caesar.
! Generates supercells, and prepares harmonic DFT calculations.
! ======================================================================
module setup_harmonic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function setup_harmonic_keywords() result(keywords)
  use help_module
  implicit none
  
  type(KeywordData) :: keywords(3)
  
  keywords = [                                                                &
  & make_keyword('dft_code', 'castep', 'dft_code is the DFT code used to &
     &calculate energies. Settings are: castep vasp qe.'),                    &
  & make_keyword('seedname', NO_ARGUMENT, 'seedname is the DFT seedname from &
     &which file names are constructed.'),                                    &
  & make_keyword('q-point_grid', NO_ARGUMENT, 'q-point_grid is the number of &
     &q-points in each direction in a Monkhorst-Pack grid. This should be &
     &specified as three integers separated by spaces.')                      ]
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine setup_harmonic(arguments)
  use utils_module, only : mkdir
  use structure_module
  use supercell_module
  use group_module
  use qpoints_module
  use dictionary_module
  use dft_input_file_module
  use generate_supercells_module
  use construct_supercell_module
  use unique_directions_module
  use calculate_symmetry_group_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! User input variables
  type(String) :: wd
  type(String) :: dft_code
  type(String) :: seedname
  
  ! File input data
  type(StructureData) :: structure
  integer             :: grid(3)
  
  ! Supercell data
  type(GeneratedSupercells)        :: qpoints_and_supercells
  integer                          :: no_supercells
  type(SupercellData), allocatable :: supercells(:)
  type(StructureData)              :: structure_sc
  
  ! Symmetry group data
  type(Group), allocatable :: symmetry_group(:)
  
  ! Directories
  type(String) :: sdir
  type(String) :: paths(2)
  
  ! Perturbation direction information.
  type(UniqueDirections) :: unique_directions
  integer                :: atom
  integer                :: direction_int
  character(1)           :: direction_char
  
  ! Temporary variables.
  integer        :: i,j,k
  
  type(String) :: dft_input_filename
  
  ! File units
  integer :: no_supercells_file
  
  ! ----------------------------------------------------------------------
  ! Get settings from user, and check them.
  ! ----------------------------------------------------------------------
  wd = item(arguments, 'working_directory')
  dft_code = item(arguments, 'dft_code')
  seedname = item(arguments, 'seedname')
  grid = int(split(item(arguments, 'q-point_grid')))
  
  ! Check dft code is supported
  if (dft_code=='vasp') then
    call print_line('')
    call print_line('Error: vasp is not currently supported.')
    stop
  elseif (dft_code/='castep' .and. dft_code/='qe') then
    call print_line('')
    call print_line('Error: The code '//dft_code//' is not supported.')
    call print_line('Please choose one of: castep vasp qe.')
    stop
  endif
  
  ! Check dft input files exist
  dft_input_filename = make_dft_input_filename(dft_code,seedname)
  dft_input_filename = wd//'/'//dft_input_filename
  
  if (.not. file_exists(dft_input_filename)) then
    call print_line('Error: The input file '//dft_input_filename// &
       &' does not exist.')
    stop
  endif
  
  ! ----------------------------------------------------------------------
  ! Read in input files.
  ! ----------------------------------------------------------------------
  structure = dft_input_file_to_StructureData(dft_code,dft_input_filename)
  call write_structure_file(structure,wd//'/structure.dat')
  
  ! ----------------------------------------------------------------------
  ! Generate supercells.
  ! ----------------------------------------------------------------------
  ! Generate IBZ and non-diagonal supercells
  qpoints_and_supercells = generate_supercells(structure,grid)
  supercells = qpoints_and_supercells%supercells
  
  ! Write q-point data to file.
  call write_structure_file( qpoints_and_supercells%structure_grid, &
                           & wd//'/structure_grid.dat')
  call write_qpoints_file( qpoints_and_supercells%qpoints_ibz, &
                         & wd//'/qpoints_ibz.dat')
  
  ! Write no_supercells to file
  no_supercells = size(supercells)
  no_supercells_file = open_write_file(wd//'/no_sc.dat')
  call print_line(no_supercells_file,no_supercells)
  close(no_supercells_file)
  
  ! ----------------------------------------------------------------------
  ! Generate supercell structures.
  ! ----------------------------------------------------------------------
  do i=1,no_supercells
    sdir=wd//'/Supercell_'//i
    
    call mkdir(sdir)
    
    ! Generate supercell structure.
    structure_sc = construct_supercell(structure, supercells(i))
    call write_structure_file(structure_sc, sdir//'/structure.dat')
    
    ! ----------------------------------------------------------------------
    ! Calculate symmetry group.
    ! ----------------------------------------------------------------------
    symmetry_group = calculate_symmetry_group(structure_sc)
    call write_group_file(symmetry_group,sdir//'/symmetry_group.dat')
    
    ! ----------------------------------------------------------------------
    ! Generate force constants
    ! ----------------------------------------------------------------------
    ! Calculate which forces need calculating
    unique_directions = calculate_unique_directions( structure_sc, &
                                                   & symmetry_group)
    call write_unique_directions_file( unique_directions, &
                                     & sdir//'/unique_directions.dat')
    
    ! ----------------------------------------------------------------------
    ! Write DFT code input files.
    ! ----------------------------------------------------------------------
    do j=1,size(unique_directions)
      atom = unique_directions%atoms(j)
      direction_int = unique_directions%directions_int(j)
      direction_char = unique_directions%directions_char(j)
      
      paths = [ sdir//'/atom.'//atom//'.+d'//direction_char, &
              & sdir//'/atom.'//atom//'.-d'//direction_char  ]
      
      ! Make harmonic run directories.
      do k=1,2
        call mkdir(paths(k))
        
        ! Move relevant atom.
        if(k==1) then
          structure_sc%atoms(direction_int,atom) =      &
             &   structure_sc%atoms(direction_int,atom) &
             & + 0.01_dp
        elseif(k==2) then
          structure_sc%atoms(direction_int,atom) =      &
             &   structure_sc%atoms(direction_int,atom) &
             & - 0.02_dp
        endif
          
        ! Write dft input files.
        dft_input_filename = make_dft_input_filename(dft_code,seedname)
        call StructureData_to_dft_input_file( &
           & dft_code,                        &
           & structure_sc,                    &
           & wd//'/'//dft_input_filename,     &
           & paths(k)//'/'//dft_input_filename)
      enddo
      
      ! Reset moved atom.
      structure_sc%atoms(direction_int,atom) =      &
         &   structure_sc%atoms(direction_int,atom) &
         & + 0.01_dp
    enddo
  enddo
end subroutine
end module
