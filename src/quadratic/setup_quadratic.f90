module setup_quadratic_module

contains

subroutine setup_quadratic()
  use constants, only : dp, thermal
  use string_module
  use file_module
  use mapping_module
  use structure_module
  use displacement_patterns_module
  use structure_to_dft_module
  use supercell_module
  implicit none
  
  ! Parameters
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
  type(SupercellData), allocatable :: supercells(:)
  type(StructureData), allocatable :: structure_scs(:)
  type(DispPatterns),  allocatable :: disp_patterns(:)
  
  ! List data
  integer              :: no_kpoints
  integer, allocatable :: sc_ids(:)
  integer, allocatable :: gvectors(:)
  
  ! Working variables
  real(dp)            :: frequency
  type(StructureData) :: structure_mode
  real(dp)            :: occupation
  real(dp)            :: quad_amplitude
  real(dp)            :: amplitude
  real(dp)            :: disp(3)
  
  ! Temporary variables
  integer        :: i, j, k, l
  type(String)   :: filename
  type(String)   :: sdir
  type(String)   :: mdir
  character(100) :: dump
  
  ! File units
  integer :: user_input_file
  integer :: no_sc_file
  integer :: list_file
  
  ! ------------------------------------------------------------
  ! Get user inputs
  ! ------------------------------------------------------------
  ! Get code name
  call print_line("What code do you want to use (castep,vasp,qe)?")
  dft_code = read_line_from_user()
  
  ! Check code is supported
  if (dft_code=="vasp") then
    call print_line("Error! vasp is not currently supported.")
    stop
  elseif (dft_code/="castep" .and. dft_code/="qe") then
    call print_line("Error! The code "//dft_code//" is not supported.")
    call print_line("Please choose one of: castep vap qe.")
    stop
  endif
  
  ! Get seedname
  if (dft_code=="castep" .or. dft_code=="vasp") then
    call print_line("What is the "//dft_code//" seedname?")
    seedname = read_line_from_user()
    seedname_nscf = seedname//'.nscf' ! only needed for qe
  endif
  
  ! Check dft input files exist
  if (dft_code=="castep") then
    filename = dft_code//'/'//seedname//'.param'
  elseif (dft_code=="qe") then
    filename = dft_code//'/'//seedname//'.in'
  endif
  
  if (.not. file_exists(filename)) then
    call print_line("Error! The input file "//filename//" does not exist.")
    stop
  endif
  
  call print_line("What is the path to the harmonic directory?")
  harmonic_path = read_line_from_user()
  
  ! ------------------------------------------------------------
  ! Read in mapping file
  ! ------------------------------------------------------------
  mapping = read_mapping_file('mapping.dat')
  
  ! ------------------------------------------------------------
  ! Read in data from harmonic calculation
  ! ------------------------------------------------------------
  ! Read in structure file
  structure = read_structure_file(harmonic_path//'/structure.dat')
  
  ! Read in number of supercells
  no_sc_file = open_read_file(harmonic_path//'/no_sc.dat')
  read(no_sc_file,*) no_sc
  close(no_sc_file)
  
  ! Read in kpoint data
  no_kpoints = count_lines(harmonic_path//'/list.dat')
  allocate(sc_ids(no_kpoints))
  allocate(gvectors(no_kpoints))
  list_file = open_read_file(harmonic_path//'/list.dat')
  do i=1,no_kpoints
    read(list_file,*) dump,dump,dump,dump, sc_ids(i), gvectors(i)
  enddo
  close(list_file)
  
  ! ------------------------------------------------------------
  ! Read in supercell structure and disp_patterns files
  ! ------------------------------------------------------------
  supercells = read_supercells_file(str('supercells.dat'))
  
  allocate(structure_scs(no_sc))
  allocate(disp_patterns(i))
  do i=1,no_sc
    filename = harmonic_path//'/Structure_'//i//'/structure.dat'
    structure_scs(i) = read_structure_file(filename)
    
    filename = harmonic_path//'/Structure_'//i//'/lte/disp_patterns.dat'
    disp_patterns(i) = read_disp_patterns_file(filename, structure%no_modes)
  enddo
  
  ! ------------------------------------------------------------
  ! Make directories
  ! ------------------------------------------------------------
  do i=1,no_sc
    call system('mkdir Supercell_'//i)
    call system('mkdir Supercell_'//i//'/static')
  enddo
  
  do i=1,no_kpoints
    call system('mkdir kpoint.'//i)
    do j=1,structure%no_modes
      do k=mapping%first,mapping%last
        call system('mkdir kpoint.'//i//'/mode.'//j//'.'//k)
      enddo
    enddo
  enddo
  
  ! ------------------------------------------------------------
  ! Set up static calculations
  ! ------------------------------------------------------------
  do i=1,no_sc
    sdir = 'Supercell_'//i
    if (dft_code=="caesar") then
      call structure_to_dft(                                      &
         & dft_code           = dft_code,                         &
         & structure_sc       = structure_scs(i),                 &
         & input_filename     = dft_code//'/'//seedname//'.cell', &
         & path_filename      = dft_code//'/path.dat',            &
         & structure          = structure,                        &
         & output_filename    = sdir//'/static/'//seedname//'.cell')
    elseif (dft_code=="vasp") then
      call structure_to_dft(                   &
         & dft_code        = dft_code,         &
         & structure_sc    = structure_scs(i), &
         & output_filename = sdir//'/POSCAR')
    elseif (dft_code=="qe") then
      call structure_to_dft(                                  &
         & dft_code         = dft_code,                       &
         & structure_sc     = structure_scs(i),               &
         & input_filename   = dft_code//'/'//seedname//'.in', &
         & pseudo_filename  = dft_code//'/pseudo.in',         &
         & kpoints_filename = dft_code//'/kpoints.in',        &
         & structure        = structure,                      &
         & output_filename  = sdir//'/static/'//seedname//'.in')
      call structure_to_dft(                                       &
         & dft_code         = dft_code,                            &
         & structure_sc     = structure_scs(i),                    &
         & input_filename   = dft_code//'/'//seedname_nscf//'.in', &
         & pseudo_filename  = dft_code//'/pseudo.in',              &
         & kpoints_filename = dft_code//'/kpoints.nscf.in',        &
         & structure        = structure,                           &
         & output_filename  = sdir//'/static/'//seedname_nscf//'.in')
    endif
  enddo
  
  ! ------------------------------------------------------------
  ! Generate quadratic configurations
  ! ------------------------------------------------------------
  do i=1,no_kpoints
    do j=1,structure%no_modes
      do k=mapping%first,mapping%last
        
        frequency = disp_patterns(sc_ids(i))%frequencies(j,gvectors(i))
        
        ! Skip acoustic modes and equilibrium configurations
        if (frequency<=frequency_tol .or. k==0) cycle
        
        ! Copy structure data
        structure_mode = structure_scs(sc_ids(i))
        
        ! Calculate quad_amplitude
        ! The normal mode amplitudes to sample are calculated
        if (temperature<temperature_tol) then
          ! Here quad_amplitude is really the normal mode amplitude squared,
          !   with units of 1/(E_h)
          occupation = 0.d0
        else
          occupation = 1.d0/(dexp(frequency/thermal_energy)-1.d0)
        endif
        quad_amplitude = dsqrt((occupation+0.5d0)/frequency)
        
        ! Calculate amplitude
        amplitude = mapping%max*k*quad_amplitude*2.d0/mapping%count
        
        ! Calculate new positions
        do l=1,structure_mode%no_atoms
          disp = amplitude                                                   &
             & * disp_patterns(sc_ids(i))%disp_patterns(1:3,l,j,gvectors(i)) &
             & * disp_patterns(sc_ids(i))%disp_patterns(4:6,l,j,gvectors(i))
          structure_mode%atoms(:,l) = structure_mode%atoms(:,l) + disp
        enddo
        
        ! Write dft input files
        sdir = 'Supercell_'//sc_ids(i)
        mdir = 'kpoint.'//i//'/mode.'//j//'.'//k
        
        if (dft_code=="castep") then
          call structure_to_dft(                                     &
             & dft_code          = dft_code,                         &
             & structure_sc      = structure_mode,                   &
             & input_filename    = dft_code//'/'//seedname//'.cell', &
             & path_filename     = dft_code//'/path.dat',            &
             & structure         = structure,                        &
             & output_filename   = mdir//'/'//seedname//'.cell')
        elseif (dft_code=="vasp") then
          call structure_to_dft(                 &
             & dft_code        = dft_code,       &
             & structure_sc    = structure_mode, &
             & output_filename = mdir//'/POSCAR.'//j//'.'//k)
        elseif (dft_code=="qe") then
          call structure_to_dft(                                  &
             & dft_code         = dft_code,                       &
             & structure_sc     = structure_mode,                 &
             & input_filename   = dft_code//'/'//seedname//'.in', &
             & pseudo_filename  = dft_code//'/pseudo.in',         &
             & kpoints_filename = dft_code//'/kpoints.in',        &
             & structure        = structure,                      &
             & output_filename  = mdir//'/'//seedname//'.in')
          if (file_exists(dft_code//'/'//seedname_nscf//'.in')) then
            call structure_to_dft(                                      &
               & dft_code         = dft_code,                            &
               & structure_sc     = structure_mode,                      &
               & input_filename   = dft_code//'/'//seedname_nscf//'.in', &
               & pseudo_filename  = dft_code//'/pseudo.in',              &
               & kpoints_filename = dft_code//'/kpoints.nscf.in',        &
               & structure        = structure,                           &
               & output_filename  = mdir//'/'//seedname_nscf//'.in')
          endif
        endif
      enddo
    enddo
  enddo
  
  ! ------------------------------------------------------------
  ! Write user inputs to file
  ! ------------------------------------------------------------
  user_input_file = open_write_file('user_input.txt')
  call print_line(user_input_file,dft_code)
  call print_line(user_input_file,seedname)
  call print_line(user_input_file,harmonic_path)
  close(user_input_file)
  
end subroutine
end module
