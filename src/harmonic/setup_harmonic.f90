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
  use generate_supercells_module
  use construct_supercell_module
  use calculate_force_constants_module
  implicit none
  
  ! User input variables
  type(String) :: dft_code
  type(String) :: seedname
  
  ! File input data
  type(StructureData) :: structure
  integer             :: grid(3)
  
  ! Supercell data
  integer              :: no_supercells
  integer, allocatable :: supercells(:,:,:)
  type(StructureData)  :: structure_sc
  
  ! Directories
  type(String) :: sdir
  type(String) :: ddir
  type(String) :: paths(2)
  
  ! Force constants data
  integer, allocatable :: force_constants(:,:)
  integer :: atom
  integer :: disp
  
  ! Temporary variables
  integer        :: i,j,k
  character(100) :: line
  type(String)   :: filename
  
  ! File units
  integer :: dft_code_file
  integer :: seedname_file
  integer :: grid_file
  integer :: supercells_file
  integer :: no_supercells_file
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
  ! Generate supercells
  ! ----------------------------------------------------------------------
  
  ! Add symmetries to structure.dat
  call system('caesar calculate_symmetry structure.dat')
  
  ! Read in input files
  structure = read_structure_file('structure.dat')
  
  grid_file = open_read_file('grid.dat')
  read(grid_file,*) grid
  close(grid_file)
  
  ! Generate IBZ and non-diagonal supercells
  call generate_supercells(structure,grid,str('ibz.dat'),str('supercells.dat'))
  
  ! Read in supercell data
  no_supercells = count_lines('supercells.dat')/4
  allocate(supercells(3,3,no_supercells))
  supercells_file = open_read_file('supercells.dat')
  do i=1,no_supercells
    read(supercells_file,*)
    do j=1,3
      read(supercells_file,*) supercells(j,:,i)
    enddo
  enddo
  close(supercells_file)
  
  no_supercells_file = open_write_file('no_sc.dat')
  write(no_supercells_file,*) no_supercells
  close(no_supercells_file)
  
  ! ----------------------------------------------------------------------
  ! Generate force constants
  ! ----------------------------------------------------------------------
  do i=1,no_supercells
    sdir=str('Supercell_')//i
    
    call system('mkdir '//sdir)
    
    ! Make supercell structure.dat file
    structure_sc = construct_supercell(structure, supercells(:,:,i))
    
    ! Add symmetries to supercell structure.dat
    call write_structure_file(structure_sc, sdir//'/structure.dat')
    call system('caesar calculate_symmetry '//sdir//'/structure.dat')
    structure_sc = read_structure_file(sdir//'/structure.dat')
    
    force_constants = calculate_force_constants( structure,         &
                                               & supercells(:,:,i), &
                                               & structure_sc)
    
    force_constants_file = open_write_file(sdir//'/force_constants.dat')
    do j=1,size(force_constants,2)
      write(force_constants_file,*) force_constants(:,j)
    enddo
    close(force_constants_file)
    
    do j=1,size(force_constants,2)
      atom = force_constants(1,j)
      disp = force_constants(2,j)
      ddir = sdir//'/atom.'//atom//'.disp.'//disp
      paths = (/ ddir//'/positive', ddir//'/negative' /)
      call system('mkdir '//ddir)
      do k=1,2
        call system('mkdir '//paths(k))
      enddo
      
      ! ----------------------------------------------------------------------
      ! Write DFT code input files
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
      
      ! Reset moved atom
      structure_sc%atoms(disp,atom) = structure_sc%atoms(disp,atom) + 0.01_dp
    enddo
  enddo
end subroutine
end module
