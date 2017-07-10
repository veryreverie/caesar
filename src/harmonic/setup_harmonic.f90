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
  use group_module
  use qpoints_module
  use dictionary_module
  use dft_input_file_module
  use generate_supercells_module
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
  type(StructureData)              :: supercell
  
  ! Symmetry group data
  type(Group), allocatable :: atom_symmetry_group(:)
  type(Group), allocatable :: operator_symmetry_group(:)
  integer,     allocatable :: symmetry_orders(:)
  integer,     allocatable :: irreducible_symmetries(:)
  type(RealMatrix), allocatable :: cartesian_rotations(:)
  
  ! Directories
  type(String) :: sdir
  type(String) :: paths(2)
  
  ! Perturbation direction information.
  type(UniqueDirections) :: unique_directions
  integer                :: atom
  character(1)           :: direction_char
  type(RealVector)       :: displacement
  
  ! Temporary variables.
  integer        :: i,j,k
  
  type(String) :: dft_input_filename
  
  ! File units
  integer :: no_supercells_file
  
  ! --------------------------------------------------
  ! Get settings from user, and check them.
  ! --------------------------------------------------
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
  
  ! --------------------------------------------------
  ! Read in input files.
  ! --------------------------------------------------
  structure = dft_input_file_to_StructureData(dft_code,dft_input_filename)
  call write_structure_file(structure,wd//'/structure.dat')
  
  ! --------------------------------------------------
  ! Generate supercells.
  ! --------------------------------------------------
  ! Generate IBZ and non-diagonal supercells
  qpoints_and_supercells = generate_supercells(structure,grid)
  
  ! Write q-point data to file.
  call write_structure_file( qpoints_and_supercells%structure_grid, &
                           & wd//'/structure_grid.dat')
  call write_qpoints_file( qpoints_and_supercells%qpoints_ibz, &
                         & wd//'/qpoints_ibz.dat')
  
  ! Write no_supercells to file
  no_supercells = size(qpoints_and_supercells%supercells)
  no_supercells_file = open_write_file(wd//'/no_supercells.dat')
  call print_line(no_supercells_file,no_supercells)
  close(no_supercells_file)
  
  ! Loop over supercells.
  do i=1,no_supercells
    sdir=wd//'/Supercell_'//i
    
    call mkdir(sdir)
    
    ! Write out structure files.
    supercell = qpoints_and_supercells%supercells(i)
    call write_structure_file(supercell, sdir//'/structure.dat')
    
    ! Calculate symmetry groups.
    atom_symmetry_group = calculate_atom_symmetry_group(supercell)
    operator_symmetry_group = calculate_operator_symmetry_group(supercell, &
       & atom_symmetry_group)
    !symmetry_orders = calculate_symmetry_orders(supercell, &
    !   & operator_symmetry_group)
    !irreducible_symmetries = calculate_irreducible_symmetries( &
    !   & supercell, operator_symmetry_group, symmetry_orders)
    cartesian_rotations = calculate_cartesian_rotations(supercell)
    !call print_line('')
    !call print_line(irreducible_symmetries)
    !do j=1,size(irreducible_symmetries)
    !  call print_line(cartesian_rotations(irreducible_symmetries(j)))
    !  call print_line(symmetry_orders(irreducible_symmetries(j)))
    !  call print_line('')
    !enddo
    !call print_line('')
    !stop
    call write_group_file(atom_symmetry_group,sdir//'/atom_symmetry_group.dat')
    
    ! Calculate which forces need calculating.
    unique_directions = calculate_unique_directions( supercell, &
                                                   & atom_symmetry_group)
    call write_unique_directions_file( unique_directions, &
                                     & sdir//'/unique_directions.dat')
    
    ! --------------------------------------------------
    ! Write DFT code input files.
    ! --------------------------------------------------
    do j=1,size(unique_directions)
      atom = unique_directions%atoms(j)
      direction_char = unique_directions%directions_char(j)
      if (direction_char=='x') then
        displacement = [1.0_dp,0.0_dp,0.0_dp]
      elseif (direction_char=='y') then
        displacement = [0.0_dp,1.0_dp,0.0_dp]
      else
        displacement = [0.0_dp,0.0_dp,1.0_dp]
      endif
      
      paths = [ sdir//'/atom.'//atom//'.+d'//direction_char, &
              & sdir//'/atom.'//atom//'.-d'//direction_char  ]
      
      ! Make harmonic run directories.
      do k=1,2
        call mkdir(paths(k))
        
        ! Move relevant atom.
        if(k==1) then
          supercell%atoms(atom) = supercell%atoms(atom) + displacement
        elseif(k==2) then
          supercell%atoms(atom) = supercell%atoms(atom) - 2*displacement
        endif
          
        ! Write dft input files.
        dft_input_filename = make_dft_input_filename(dft_code,seedname)
        call StructureData_to_dft_input_file( &
           & dft_code,                        &
           & supercell,                    &
           & wd//'/'//dft_input_filename,     &
           & paths(k)//'/'//dft_input_filename)
      enddo
      
      ! Reset moved atom.
      supercell%atoms(atom) = supercell%atoms(atom) + displacement
    enddo
  enddo
end subroutine
end module
