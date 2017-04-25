! ======================================================================
! The first stage of Caesar's anharmonic process.
! Generates configurations along normal modes, and prepares anharmonic
!    DFT calculations.
! ======================================================================
module setup_quadratic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function setup_quadratic_keywords() result(keywords)
  use help_module
  implicit none
  
  type(KeywordData) :: keywords(3)
  
  keywords = [                                                                &
  & make_keyword('dft_code', 'castep', 'dft_code is the DFT code used to &
     &calculate energies. Settings are: castep vasp qe.'),                    &
  & make_keyword('seedname', no_argument, 'seedname is the DFT seedname from &
     &which file names are constructed.'),                                    &
  & make_keyword('harmonic_path', '.', 'harmonic_path is the path to the &
     &directory where harmonic calculations were run.', is_path=.true.)       ]
end function

! ----------------------------------------------------------------------
! The main program.
! ----------------------------------------------------------------------
subroutine setup_quadratic(arguments)
  use constants_module, only : kb_in_au
  use utils_module,     only : mkdir
  use mapping_module
  use structure_module
  use displacement_patterns_module
  use dft_input_file_module
  use supercell_module
  use kpoints_module
  use dictionary_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! Parameters.
  real(dp), parameter :: frequency_tol   = 1.0e-5_dp!frequency<=tol if acoustic
  real(dp), parameter :: temperature     = 0.0_dp
  real(dp), parameter :: temperature_tol = 1.0e-6_dp
  real(dp), parameter :: thermal_energy  = temperature*kb_in_au
  
  ! User inputs
  type(String) :: dft_code      ! The dft code name (castep,vasp,qe).
  type(String) :: seedname      ! The dft input file seedname.
  type(String) :: harmonic_path ! The path to the harmonic directory.
  
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
  type(KpointData), allocatable :: kpoints(:)
  
  ! Working variables
  type(StructureData) :: structure_sc
  real(dp)            :: occupation
  real(dp)            :: quad_amplitude
  real(dp)            :: amplitude
  real(dp)            :: disp(3)
  
  ! Temporary variables
  integer        :: i, j, k, l
  type(String)   :: dft_input_filename
  type(String)   :: sdir
  type(String)   :: mdir
  type(String), allocatable :: line(:)
  
  ! Files
  type(String), allocatable :: no_sc_file(:)
  
  ! ------------------------------------------------------------
  ! Get settings from user.
  ! ------------------------------------------------------------
  wd = item(arguments, 'working_directory')
  dft_code = item(arguments, 'dft_code')
  seedname = item(arguments, 'seedname')
  harmonic_path = item(arguments, 'harmonic_path')
  
  ! Check code is supported
  if (dft_code=="vasp") then
    call print_line('')
    call print_line("Error: vasp is not currently supported.")
    stop
  elseif (dft_code/="castep" .and. dft_code/="qe") then
    call print_line('')
    call print_line("Error: The code "//dft_code//" is not supported.")
    call print_line("Please choose one of: castep vap qe.")
    stop
  endif
  
  ! Check dft input files exist
  dft_input_filename = make_dft_input_filename(dft_code,seedname)
  dft_input_filename = wd//'/'//dft_code//'/'//dft_input_filename
  
  if (.not. file_exists(dft_input_filename)) then
    call print_line("Error: The input file "//dft_input_filename// &
       & " does not exist.")
    stop
  endif
  
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
    structure_scs(i) = read_structure_file( &
       & harmonic_path//'/Supercell_'//i//'/structure.dat')
  enddo
  
  ! Read in kpoint data.
  kpoints = read_kpoints_file(harmonic_path//'/kpoints_ibz.dat')
  
  ! ------------------------------------------------------------
  ! Make directories
  ! ------------------------------------------------------------
  do i=1,no_sc
    call mkdir(wd//'/Supercell_'//i)
  enddo
  
  do i=1,size(kpoints)
    call mkdir(wd//'/kpoint_'//i)
  enddo
  
  ! ------------------------------------------------------------
  ! Set up static calculations
  ! ------------------------------------------------------------
  do i=1,no_sc
    sdir = wd//'/Supercell_'//i
    dft_input_filename = make_dft_input_filename(dft_code,seedname)
    call StructureData_to_dft_input_file(            &
       & dft_code,                                   &
       & structure_scs(i),                           &
       & wd//'/'//dft_code//'/'//dft_input_filename, &
       & sdir//'/'//dft_input_filename)
  enddo
  
  ! ------------------------------------------------------------
  ! Generate quadratic configurations
  ! ------------------------------------------------------------
  do i=1,size(kpoints)
    structure_sc = structure_scs(kpoints(i)%sc_id)
    
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
          occupation = 1.d0/(exp(frequencies(j)/thermal_energy)-1.d0)
        endif
        quad_amplitude = sqrt((occupation+0.5d0)/frequencies(j))
        
        ! Calculate amplitude
        amplitude = mapping%max*k*quad_amplitude*2.d0/mapping%count
        
        ! Calculate new positions
        do l=1,structure_sc%no_atoms
          disp = amplitude * displacements(:,l,j) * prefactors(l,j)
          structure_sc%atoms(:,l) = structure_sc%atoms(:,l) + disp
        enddo
        
        ! Write dft input files
        mdir = wd//'/kpoint_'//i//'/mode_'//j//'/amplitude_'//k
        
        dft_input_filename = make_dft_input_filename(dft_code,seedname)
        call StructureData_to_dft_input_file( &
           & dft_code,                        &
           & structure_sc,                    &
           & wd//'/'//dft_input_filename,     &
           & mdir//'/'//dft_input_filename)
      enddo
    enddo
    
    deallocate(frequencies)
    deallocate(prefactors)
    deallocate(displacements)
  enddo
  
end subroutine
end module
