module setup_harmonic_module

contains

! ======================================================================
! Program to set up a generic harmonic calculation to use with LTE
! Also converts the generic calculateion to castep vasp or qe
! ======================================================================
subroutine setup_harmonic()
  use string_module
  use file_module
  use structure_module
  
  use structure_to_dft_module
  use generate_kgrid_module
  use generate_supercells_module
  use construct_supercell_module
  use construct_matrix_force_cnsts_module
  use construct_finite_displacement_module
  implicit none
  
  ! User input variables
  type(String) :: dft_code
  type(String) :: seedname
  
  ! File input data
  type(StructureData) :: structure
  integer             :: grid(3)
  
  ! Supercell data
  integer             :: no_supercells
  type(StructureData) :: structure_sc
  integer             :: supercell(3,3)
  
  ! Directories
  type(String) :: sdir
  type(String) :: ddir
  type(String) :: paths(2)
  
  ! Force constants data
  integer :: atom
  integer :: disp
  
  ! Temporary variables
  integer        :: i,j,k
  character(100) :: line
  type(String)   :: filename
  integer        :: file_length
  
  ! File units
  integer :: dft_code_file
  integer :: seedname_file
  integer :: grid_file
  integer :: no_supercells_file
  integer :: supercell_file
  integer :: force_constants_file
  
  ! ----------------------------------------------------------------------
  ! Get settings from user
  ! ----------------------------------------------------------------------
  
  ! Get dft code
  write(*,"(a)") "What dft code do you want to use (castep,vasp,qe)?"
  read(*,"(a)") line
  dft_code = line
  
  ! Check dft code is supported
  if (dft_code=="vasp") then
    write(*,"(a)") "Error! vasp is not currently supported."
    stop
  elseif (dft_code/="castep" .and. dft_code/="qe") then
    write(*,"(a)") "Error! The code "//char(dft_code)//" is not supported."
    write(*,"(a)") "Please choose one of: castep vap qe."
    stop
  endif
  
  ! Get seedname
  if (dft_code=="castep" .or. dft_code=="vasp") then
    write(*,"(a)") "What is the "//char(dft_code)//" seedname?"
    read(*,*) line
    seedname = line
  endif
  
  ! Check dft input files exist
  if (dft_code=="castep") then
    filename = dft_code//'/'//seedname//'.param'
  elseif (dft_code=="qe") then
    filename = dft_code//'/'//seedname//'.in'
  endif
  
  if (.not. file_exists(filename)) then
    write(*,"(a)") "Error! The input file "//char(filename)//" does not exist."
    stop
  endif
  
  ! ----------------------------------------------------------------------
  ! Write user settings to file
  ! ----------------------------------------------------------------------
  
  dft_code_file = open_write_file('code.txt')
  write(dft_code_file,"(a)") char(dft_code)
  close(dft_code_file)
  
  seedname_file = open_write_file('seedname.txt')
  write(seedname_file,"(a)") char(seedname)
  close(seedname_file)
  
  ! ----------------------------------------------------------------------
  ! Generate generic calculation
  ! ----------------------------------------------------------------------
  
  ! Add symmetries to structure.dat
  call system('caesar calculate_symmetry structure.dat')
  
  ! Read in input files
  structure = read_structure_file('structure.dat')
  
  grid_file = open_read_file('grid.dat')
  read(grid_file,*) grid
  close(grid_file)
  
  ! Generate IBZ
  ! Reads Caesar input files. Writes ibz.dat and rotated_gvectors.dat
  call generate_kgrid( structure,      &
                     & grid,           &
                     & str('ibz.dat'), &
                     & str('rotated_gvectors.dat'))
  
  ! Generate non-diagonal supercells
  ! Makes Supercell_* directories
  ! Adds to Supercell_* directories:
  !   supercell.dat
  ! Reads and writes ibz.dat
  !   Adds sc_id to each kpoint
  !   Orders kpoints by supercell
  call generate_supercells( structure,        &
                          & grid,             &
                          & str('ibz.dat'),   &
                          & str('no_sc.dat'), &
                          & str('Supercell_'))
  
  no_supercells_file = open_read_file('no_sc.dat')
  read(no_supercells_file,*) no_supercells
  close(no_supercells_file)
  
  ! ----------------------------------------------------------------------
  ! Generate force constants
  ! ----------------------------------------------------------------------
  do i=1,no_supercells
    sdir=str('Supercell_')//i
    
    ! Read in supercell data
    supercell_file = open_read_file(sdir//'/supercell.dat')
    do j=1,3
      read(supercell_file,*) supercell(j,:)
    enddo
    close(supercell_file)
    
    ! Make supercell structure.dat file
    call construct_supercell(structure, supercell, sdir//'/structure.dat')
    
    ! Add symmetries to supercell structure.dat
    call system('caesar calculate_symmetry '//sdir//'/structure.dat')
    
    ! Read in supercell structure data
    structure_sc = read_structure_file(sdir//'/structure.dat')
    
    call construct_matrix_force_cnsts( structure,    &
                                     & supercell,    &
                                     & structure_sc, &
                                     & sdir//'/force_constants.dat')
    
    file_length = count_lines(sdir//'/force_constants.dat')
    force_constants_file = open_read_file(sdir//'/force_constants.dat')
    do j=1,file_length
      read(force_constants_file,*) atom,disp
      ddir = sdir//'/atom.'//atom//'.disp.'//disp
      paths = (/ ddir//'/positive', ddir//'/negative' /)
      call system('mkdir '//ddir)
      do k=1,2
        call system('mkdir '//paths(k))
      enddo
      
      ! Write structures with atoms moved by +- 0.01
      structure_sc%atoms(disp,atom) = structure_sc%atoms(disp,atom) + 0.01_dp
      call write_structure_file(structure_sc,paths(1)//'/structure.dat')
      structure_sc%atoms(disp,atom) = structure_sc%atoms(disp,atom) - 0.02_dp
      call write_structure_file(structure_sc,paths(2)//'/structure.dat')
      structure_sc%atoms(disp,atom) = structure_sc%atoms(disp,atom) + 0.01_dp
      
      ! ----------------------------------------------------------------------
      ! Convert generic calculation to specific dft code
      ! ----------------------------------------------------------------------
      
      do k=1,2
        
        ! Move relevant atom
        if(k==1) then
          structure_sc%atoms(disp,atom) = structure_sc%atoms(disp,atom) &
                                      & + 0.01_dp
        elseif(k==2) then
          structure_sc%atoms(disp,atom) = structure_sc%atoms(disp,atom) &
                                      & - 0.02_dp
        endif
          
        ! Write dft input files
        if (dft_code=="castep") then
          call structure_to_dft(                                   &
             & dft_code        = dft_code,                         &
             & structure_sc    = structure_sc,                     &
             & input_filename  = dft_code//'/'//seedname//'.cell', &
             & output_filename = paths(k)//'/'//seedname//'.cell')
        elseif (dft_code=="vasp") then
          call structure_to_dft(               &
             & dft_code        = dft_code,     &
             & structure_sc    = structure_sc, &
             & output_filename = paths(k)//'/POSCAR')
        elseif (dft_code=="qe") then
          call structure_to_dft(                                  &
             & dft_code         = dft_code,                       &
             & structure_sc     = structure_sc,                   &
             & input_filename   = dft_code//'/'//seedname//'.in', &
             & pseudo_filename  = dft_code//'/pseudo.in',         &
             & kpoints_filename = dft_code//'/kpoints.in',        &
             & structure        = structure,                      &
             & output_filename  = paths(k)//'/'//seedname//'.in')
        endif
      enddo
      
      ! Reset structure
      structure_sc%atoms(disp,atom) = structure_sc%atoms(disp,atom) + 0.01_dp
    enddo
    close(force_constants_file)
  enddo
end subroutine
end module
