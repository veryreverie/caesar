module setup_quadratic_module

contains

subroutine setup_quadratic(wd,cwd)
  use constants, only : dp, thermal
  use utils,     only : mkdir, format_directory
  use string_module
  use file_module
  use mapping_module
  use structure_module
  use displacement_patterns_module
  use structure_to_dft_module
  use supercell_module
  use err_module
  implicit none
  
  ! Working directories.
  type(String), intent(in) :: wd
  type(String), intent(in) :: cwd
  
  ! Parameters.
  real(dp), parameter :: frequency_tol   = 1.0e-5_dp!frequency<=tol if acoustic
  real(dp), parameter :: temperature     = 0.0_dp
  real(dp), parameter :: temperature_tol = 1.0e-6_dp
  real(dp), parameter :: thermal_energy  = temperature/thermal
  
  ! User inputs
  type(String) :: dft_code      ! The dft code name (castep,vasp,qe)
  type(String) :: seedname      ! The dft input file seedname
  type(String) :: seedname_nscf ! Only needed for qe
  type(String) :: harmonic_path ! The path to the harmonic directory
  
  ! File contents
  type(MappingData)    :: mapping
  
  ! Harmonic file contents
  type(StructureData) :: structure
  integer             :: no_sc
  
  ! Harmonic supercell file contents
  type(StructureData), allocatable :: structure_scs(:)
  
  ! Mode frequencies and displacement patterns.
  type(String), allocatable :: frequency_file(:)
  type(String), allocatable :: prefactors_file(:)
  type(String), allocatable :: displacements_file(:)
  real(dp),     allocatable :: frequencies(:)
  real(dp),     allocatable :: prefactors(:,:)
  real(dp),     allocatable :: displacements(:,:,:)
  
  ! K-point data
  integer              :: no_kpoints
  integer, allocatable :: sc_ids(:)
  integer, allocatable :: gvectors(:)
  
  ! Working variables
  type(StructureData) :: structure_sc
  real(dp)            :: occupation
  real(dp)            :: quad_amplitude
  real(dp)            :: amplitude
  real(dp)            :: disp(3)
  
  ! Temporary variables
  integer        :: i, j, k, l
  type(String)   :: filename
  type(String)   :: sdir
  type(String)   :: mdir
  type(String), allocatable :: line(:)
  
  ! Files
  type(String), allocatable :: ibz_file(:)
  type(String), allocatable :: no_sc_file(:)
  integer :: user_input_file
  
  ! ------------------------------------------------------------
  ! Get user inputs
  ! ------------------------------------------------------------
  ! Get code name
  call print_line("What code do you want to use (castep,vasp,qe)?")
  dft_code = read_line_from_user()
  
  ! Check code is supported
  if (dft_code=="vasp") then
    call print_line("Error! vasp is not currently supported.")
    call err()
  elseif (dft_code/="castep" .and. dft_code/="qe") then
    call print_line("Error! The code "//dft_code//" is not supported.")
    call print_line("Please choose one of: castep vap qe.")
    call err()
  endif
  
  ! Get seedname
  if (dft_code=="castep" .or. dft_code=="vasp") then
    call print_line("What is the "//dft_code//" seedname?")
    seedname = read_line_from_user()
    seedname_nscf = seedname//'.nscf' ! only needed for qe
  endif
  
  ! Check dft input files exist
  if (dft_code=="castep") then
    filename = wd//'/'//dft_code//'/'//seedname//'.param'
  elseif (dft_code=="qe") then
    filename = wd//'/'//dft_code//'/'//seedname//'.in'
  endif
  
  if (.not. file_exists(filename)) then
    call print_line("Error! The input file "//filename//" does not exist.")
    call err()
  endif
  
  call print_line("What is the path to the harmonic directory?")
  harmonic_path = format_directory(read_line_from_user(),cwd)
  
  ! ------------------------------------------------------------
  ! Read in mapping file
  ! ------------------------------------------------------------
  mapping = read_mapping_file(wd//'/mapping.dat')
  
  ! ------------------------------------------------------------
  ! Read in data from harmonic calculation
  ! ------------------------------------------------------------
  ! Read in structure file
  structure = read_structure_file(harmonic_path//'/structure.dat')
  
  ! Read in number of supercells
  no_sc_file = read_lines(harmonic_path//'/no_sc.dat')
  no_sc = int(no_sc_file(1))
  
  ! Read in supercell structures.
  allocate(structure_scs(no_sc))
  do i=1,no_sc
    filename = harmonic_path//'/Supercell_'//i//'/structure.dat'
    structure_scs(i) = read_structure_file(filename)
  enddo
  
  ! Read in kpoint data
  ibz_file = read_lines(harmonic_path//'/ibz.dat')
  no_kpoints = size(ibz_file)
  allocate(sc_ids(no_kpoints))
  allocate(gvectors(no_kpoints))
  do i=1,no_kpoints
    line = split(ibz_file(i))
    sc_ids(i) = int(line(5))
    gvectors(i) = int(line(6))
  enddo
  
  ! ------------------------------------------------------------
  ! Make directories
  ! ------------------------------------------------------------
  do i=1,no_sc
    call mkdir(wd//'/Supercell_'//i)
  enddo
  
  do i=1,no_kpoints
    call mkdir(wd//'/kpoint_'//i)
  enddo
  
  ! ------------------------------------------------------------
  ! Set up static calculations
  ! ------------------------------------------------------------
  do i=1,no_sc
    sdir = wd//'/Supercell_'//i
    if (dft_code=="castep") then
      call structure_to_dft(                                               &
         & dft_code           = dft_code,                                  &
         & structure_sc       = structure_scs(i),                          &
         & input_filename     = wd//'/'//dft_code//'/'//seedname//'.cell', &
         & path_filename      = wd//'/'//dft_code//'/path.dat',            &
         & structure          = structure,                                 &
         & output_filename    = sdir//'/'//seedname//'.cell')
    elseif (dft_code=="vasp") then
      call structure_to_dft(                   &
         & dft_code        = dft_code,         &
         & structure_sc    = structure_scs(i), &
         & output_filename = sdir//'/POSCAR')
    elseif (dft_code=="qe") then
      call structure_to_dft(                                           &
         & dft_code         = dft_code,                                &
         & structure_sc     = structure_scs(i),                        &
         & input_filename   = wd//'/'//dft_code//'/'//seedname//'.in', &
         & pseudo_filename  = wd//'/'//dft_code//'/pseudo.in',         &
         & kpoints_filename = wd//'/'//dft_code//'/kpoints.in',        &
         & structure        = structure,                               &
         & output_filename  = sdir//'/'//seedname//'.in')
      call structure_to_dft(                                                &
         & dft_code         = dft_code,                                     &
         & structure_sc     = structure_scs(i),                             &
         & input_filename   = wd//'/'//dft_code//'/'//seedname_nscf//'.in', &
         & pseudo_filename  = wd//'/'//dft_code//'/pseudo.in',              &
         & kpoints_filename = wd//'/'//dft_code//'/kpoints.nscf.in',        &
         & structure        = structure,                                    &
         & output_filename  = sdir//'/'//seedname_nscf//'.in')
    endif
  enddo
  
  ! ------------------------------------------------------------
  ! Generate quadratic configurations
  ! ------------------------------------------------------------
  do i=1,no_kpoints
    structure_sc = structure_scs(sc_ids(i))
    
    ! Read in frequencies, prefactors and displacements.
    frequency_file = read_lines( &
       & harmonic_path//'/kpoint_'//i//'/frequencies.dat')
    allocate(frequencies(structure%no_modes))
    do j=1,structure%no_modes
      frequencies(j) = dble(frequency_file(j))
    enddo
    
    prefactors_file = read_lines( &
       & harmonic_path//'/kpoint_'//i//'/prefactors.dat')
    allocate(prefactors(structure_sc%no_atoms,structure%no_modes))
    do j=1,structure%no_modes
      do k=1,structure_sc%no_atoms
        prefactors(k,j) = &
           & dble(prefactors_file((j-1)*(structure_sc%no_atoms+2)+k+1))
      enddo
    enddo
    
    displacements_file = read_lines( &
       & harmonic_path//'/kpoint_'//i//'/displacements.dat')
    allocate(displacements(3,structure_sc%no_atoms,structure%no_modes))
    do j=1,structure%no_modes
      do k=1,structure_sc%no_atoms
        line = split(displacements_file((j-1)*(structure_sc%no_atoms+2)+k+1))
        displacements(:,k,j) = dble(line)
      enddo
    enddo
    
    do j=1,structure%no_modes
      ! Skip acoustic modes.
      if (frequencies(j)<=frequency_tol) then
        cycle
      endif
      
      call mkdir(wd//'/kpoint_'//i//'/mode_'//j)
      
      do k=mapping%first,mapping%last
        ! Skip equilibrium configurations
        if (k==0) then
          cycle
        endif
    
        call mkdir(wd//'/kpoint_'//i//'/mode_'//j//'/amplitude_'//k)
        
        ! Calculate quad_amplitude
        ! The normal mode amplitudes to sample are calculated
        if (temperature<temperature_tol) then
          ! Here quad_amplitude is really the normal mode amplitude squared,
          !   with units of 1/(E_h)
          occupation = 0.d0
        else
          occupation = 1.d0/(dexp(frequencies(j)/thermal_energy)-1.d0)
        endif
        quad_amplitude = dsqrt((occupation+0.5d0)/frequencies(j))
        
        ! Calculate amplitude
        amplitude = mapping%max*k*quad_amplitude*2.d0/mapping%count
        
        ! Calculate new positions
        do l=1,structure_sc%no_atoms
          disp = amplitude * displacements(:,l,j) * prefactors(l,j)
          structure_sc%atoms(:,l) = structure_sc%atoms(:,l) + disp
        enddo
        
        ! Write dft input files
        sdir = wd//'/Supercell_'//sc_ids(i)
        mdir = wd//'/kpoint_'//i//'/mode_'//j//'/amplitude_'//k
        
        if (dft_code=="castep") then
          call structure_to_dft(                                              &
             & dft_code          = dft_code,                                  &
             & structure_sc      = structure_sc,                              &
             & input_filename    = wd//'/'//dft_code//'/'//seedname//'.cell', &
             & path_filename     = wd//'/'//dft_code//'/path.dat',            &
             & structure         = structure,                                 &
             & output_filename   = mdir//'/'//seedname//'.cell')
        elseif (dft_code=="vasp") then
          call structure_to_dft(               &
             & dft_code        = dft_code,     &
             & structure_sc    = structure_sc, &
             & output_filename = mdir//'/POSCAR.'//j//'.'//k)
        elseif (dft_code=="qe") then
          call structure_to_dft(                                           &
             & dft_code         = dft_code,                                &
             & structure_sc     = structure_sc,                            &
             & input_filename   = wd//'/'//dft_code//'/'//seedname//'.in', &
             & pseudo_filename  = wd//'/'//dft_code//'/pseudo.in',         &
             & kpoints_filename = wd//'/'//dft_code//'/kpoints.in',        &
             & structure        = structure,                               &
             & output_filename  = mdir//'/'//seedname//'.in')
          if (file_exists(dft_code//'/'//seedname_nscf//'.in')) then
            call structure_to_dft(                                         &
               & dft_code         = dft_code,                              &
               & structure_sc     = structure_sc,                          &
               & input_filename   = wd//'/'//dft_code//'/'//               &
               &                    seedname_nscf//'.in',                  &
               & pseudo_filename  = wd//'/'//dft_code//'/pseudo.in',       &
               & kpoints_filename = wd//'/'//dft_code//'/kpoints.nscf.in', &
               & structure        = structure,                             &
               & output_filename  = mdir//'/'//seedname_nscf//'.in')
          endif
        endif
      enddo
    enddo
    
    deallocate(frequencies)
    deallocate(prefactors)
    deallocate(displacements)
  enddo
  
  ! ------------------------------------------------------------
  ! Write user inputs to file
  ! ------------------------------------------------------------
  user_input_file = open_write_file(wd//'/user_input.txt')
  call print_line(user_input_file,dft_code)
  call print_line(user_input_file,seedname)
  call print_line(user_input_file,harmonic_path)
  close(user_input_file)
  
end subroutine
end module
