module setup_harmonic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ======================================================================
! Program to set up a harmonic calculation to use with LTE.
! ======================================================================
subroutine setup_harmonic(wd)
  use constants_module, only : directions
  use utils_module,     only : mkdir, make_dft_input_filename
  use structure_module
  use supercell_module
  use group_module
  use kpoints_module
  
  use dft_input_file_module
  use structure_to_dft_module
  use generate_supercells_module
  use construct_supercell_module
  use unique_directions_module
  use calculate_symmetry_group_module
  use calculate_symmetry_module
  implicit none
  
  ! Working directory.
  type(String), intent(in) :: wd
  
  ! User input variables
  type(String) :: dft_code
  type(String) :: seedname
  
  ! File input data
  type(StructureData) :: structure
  integer             :: grid(3)
  
  ! Supercell data
  type(GeneratedSupercells)        :: kpoints_and_supercells
  integer                          :: no_supercells
  type(SupercellData), allocatable :: supercells(:)
  type(StructureData)              :: structure_sc
  
  ! Symmetry group data
  type(Group), allocatable :: symmetry_group(:)
  
  ! Directories
  type(String) :: sdir
  type(String) :: paths(2)
  
  ! Force constants data
  type(UniqueDirections) :: unique_directions
  integer :: atom
  character(1) :: direction
  
  ! Temporary variables
  integer        :: i,j,k,l
  
  type(String), allocatable :: user_input_file_in(:)
  
  type(String) :: dft_input_filename
  
  ! File units
  integer :: user_input_file_out
  integer :: no_supercells_file
  
  ! ----------------------------------------------------------------------
  ! Get settings from user.
  ! ----------------------------------------------------------------------
  
  if (file_exists(wd//'/user_input.txt')) then
    ! Get settings from file.
    user_input_file_in = read_lines(wd//'/user_input.txt')
    dft_code = user_input_file_in(1)
    seedname = user_input_file_in(2)
    grid = int(split(user_input_file_in(3)))
  else
    ! Get dft code from the command line.
    call print_line('')
    call print_line('What dft code do you want to use (castep,vasp,qe)?')
    dft_code = read_line_from_user()
    
    ! Get seedname from the command line.
    call print_line('')
    call print_line('What is the '//dft_code//' seedname?')
    seedname = read_line_from_user()
    
    ! Get the K-point grid from the command line.
    call print_line('')
    call print_line('What is the x dimension of the K-point grid?')
    grid(1) = int(read_line_from_user())
    call print_line('')
    call print_line('What is the y dimension of the K-point grid?')
    grid(2) = int(read_line_from_user())
    call print_line('')
    call print_line('What is the z dimension of the K-point grid?')
    grid(3) = int(read_line_from_user())
    
    call print_line('')
  endif
  
  ! Check dft code is supported
  if (dft_code=='vasp') then
    call print_line('Error! vasp is not currently supported.')
    call err()
  elseif (dft_code/='castep' .and. dft_code/='qe') then
    call print_line('Error! The code '//dft_code//' is not supported.')
    call print_line('Please choose one of: castep vasp qe.')
    call err()
  endif
  
  ! Check dft input files exist
  dft_input_filename = make_dft_input_filename(dft_code,seedname)
  dft_input_filename = wd//'/'//dft_input_filename
  
  if (.not. file_exists(dft_input_filename)) then
    call print_line('Error! The input file '//dft_input_filename// &
       &' does not exist.')
    call err()
  endif
  
  ! ----------------------------------------------------------------------
  ! Read in input files.
  ! ----------------------------------------------------------------------
  !structure = read_structure_file(wd//'/structure.dat')
  structure = dft_input_file_to_structure(dft_code,dft_input_filename)
  
  ! ----------------------------------------------------------------------
  ! Write user settings to file
  ! ----------------------------------------------------------------------
  if (.not. file_exists(wd//'/user_input.txt')) then
    user_input_file_out = open_write_file(wd//'/user_input.txt')
    call print_line(user_input_file_out, dft_code)
    call print_line(user_input_file_out, seedname)
    call print_line(user_input_file_out, grid)
    close(user_input_file_out)
  endif
  
  ! ----------------------------------------------------------------------
  ! Add symmetries to structure.dat if not already present.
  ! ----------------------------------------------------------------------
  if (structure%no_symmetries == 0) then
    call calculate_symmetry(structure)
  endif
  
  call write_structure_file(structure,wd//'/structure.dat')
  
  ! ----------------------------------------------------------------------
  ! Generate supercells.
  ! ----------------------------------------------------------------------
  ! Generate IBZ and non-diagonal supercells
  kpoints_and_supercells = generate_supercells(structure,grid)
  supercells = kpoints_and_supercells%supercells
  
  ! Write K-point data to file.
  call write_kpoints_grid_file( kpoints_and_supercells%kpoints_grid, &
                              & wd//'/kpoints_grid.dat')
  call write_kpoints_ibz_file( kpoints_and_supercells%kpoints_ibz, &
                             & wd//'/kpoints_ibz.dat')
  
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
    
    ! Add symmetries to supercell structure.dat
    call calculate_symmetry(structure_sc)
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
    ! Write DFT code input files
    ! ----------------------------------------------------------------------
    do k=1,size(unique_directions)
      atom = unique_directions%unique_atoms(k)
      do j=1,3
        ! Skip directions which are covered by symmetry.
        if (j==2 .and. unique_directions%xy_symmetry(k)>0) then
          cycle
        endif
        if (j==3 .and. ( unique_directions%xz_symmetry(k)>0 .or. &
                       & unique_directions%yz_symmetry(k)>0)) then
          cycle
        endif
        
        direction = directions(j)
        paths = (/ sdir//'/atom.'//atom//'.+d'//direction, &
                 & sdir//'/atom.'//atom//'.-d'//direction /)
        
        ! Make harmonic run directories
        do l=1,2
          call mkdir(paths(l))
          
          ! Move relevant atom
          if(l==1) then
            structure_sc%atoms(j,atom) = structure_sc%atoms(j,atom) + 0.01_dp
          elseif(l==2) then
            structure_sc%atoms(j,atom) = structure_sc%atoms(j,atom) - 0.02_dp
          endif
            
          ! Write dft input files
          dft_input_filename = make_dft_input_filename(dft_code,seedname)
          call structure_to_dft( dft_code, &
                               & structure_sc, &
                               & wd//'/'//dft_input_filename, &
                               & paths(l)//'/'//dft_input_filename)
        enddo
        
        ! Reset moved atom
        structure_sc%atoms(j,atom) = structure_sc%atoms(j,atom) + 0.01_dp
      enddo
    enddo
  enddo
end subroutine
end module
