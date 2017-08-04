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
  
  type(KeywordData) :: keywords(4)
  
  keywords = [                                                                &
  & make_keyword( 'harmonic_path',                                            &
  &               'harmonic_path is the path to the directory where harmonic &
  &calculations were run.',                                                   &
  &               default_value='.',                                          &
  &               is_path=.true.),                                            &
  & make_keyword( 'temperature',                                              &
  &               'temperature is the temperature, in Kelvin, at which the &
  &simulation is run.',                                                       &
  &               default_value='0'),                                         &
  & make_keyword( 'no_samples',                                               &
  &               'no_samples is the number of single-point calculations &
  &which will be run along each axis.'),                                      &
  & make_keyword( 'displacement',                                             &
  &               'displacement is the maximum total distance in bohr over &
  &which any mode is displaced. At finite temperatures, this will be reduced &
  &by a thermal factor.') ]
end function

! ----------------------------------------------------------------------
! The main program.
! ----------------------------------------------------------------------
subroutine setup_quadratic(arguments)
  use constants_module, only : kb_in_au, pi
  use utils_module,     only : mkdir
  use structure_module
  use dft_input_file_module
  use qpoints_module
  use dictionary_module
  use normal_mode_module
  use linear_algebra_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! Parameters.
  real(dp), parameter :: frequency_tol   = 1.0e-5_dp!frequency<=tol if acoustic
  real(dp), parameter :: temperature_tol = 1.0e-6_dp
  
  ! User inputs
  type(String) :: harmonic_path
  real(dp)     :: temperature
  integer      :: no_samples
  real(dp)     :: displacement
  
  ! Previous user inputs.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: dft_code
  type(String)     :: seedname
  
  ! Harmonic file contents
  type(StructureData) :: structure
  integer             :: no_supercells
  
  ! Harmonic supercell file contents
  type(StructureData), allocatable :: supercells(:)
  
  ! Normal modes.
  type(NormalMode) :: mode
  
  ! q-point data
  type(QpointData), allocatable :: qpoints_ibz(:)
  
  ! Working variables
  real(dp)            :: thermal_energy
  type(StructureData) :: supercell
  real(dp)            :: occupation
  real(dp)            :: quad_amplitude
  real(dp)            :: amplitude
  type(RealVector)    :: disp
  
  ! Supercell-to-primitive indexes.
  integer  :: rvec
  integer  :: prim
  real(dp) :: qr
  
  ! Temporary variables
  integer        :: i, j, k, l, ialloc
  type(String)   :: dft_input_filename
  type(String)   :: sdir
  type(String)   :: mdir
  
  ! Files
  type(String), allocatable :: no_supercells_file(:)
  
  ! --------------------------------------------------
  ! Get settings from user.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  dft_code = arguments%value('dft_code')
  harmonic_path = arguments%value('harmonic_path')
  temperature = dble(arguments%value('temperature'))
  no_samples = int(arguments%value('no_samples'))
  displacement = dble(arguments%value('displacement'))
  
  thermal_energy = temperature*kb_in_au
  
  ! --------------------------------------------------
  ! Read in previous settings.
  ! --------------------------------------------------
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  seedname = setup_harmonic_arguments%value('seedname')
  dft_code = setup_harmonic_arguments%value('dft_code')
  
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
  dft_input_filename = wd//'/'//dft_input_filename
  
  if (.not. file_exists(dft_input_filename)) then
    call print_line("Error: The input file "//dft_input_filename// &
       & " does not exist.")
    stop
  endif
  
  ! Check temperature is valid.
  if (temperature < 0.0_dp) then
    call print_line('Error: temperature must be positive.')
    stop
  endif
  
  ! ------------------------------------------------------------
  ! Read in data from harmonic calculation
  ! ------------------------------------------------------------
  ! Read in structure file
  structure = read_structure_file(harmonic_path//'/structure.dat')
  
  ! Read in number of supercells
  no_supercells_file = read_lines(harmonic_path//'/no_supercells.dat')
  no_supercells = int(no_supercells_file(1))
  
  ! Read in supercell structures.
  allocate(supercells(no_supercells), stat=ialloc); call err(ialloc)
  do i=1,no_supercells
    supercells(i) = read_structure_file( &
       & harmonic_path//'/Supercell_'//i//'/structure.dat')
  enddo
  
  ! Read in qpoint data.
  qpoints_ibz = read_qpoints_file(harmonic_path//'/qpoints_ibz.dat')
  
  ! ------------------------------------------------------------
  ! Make directories
  ! ------------------------------------------------------------
  do i=1,no_supercells
    call mkdir(wd//'/Supercell_'//i)
  enddo
  
  do i=1,size(qpoints_ibz)
    call mkdir(wd//'/qpoint_'//i)
  enddo
  
  ! ------------------------------------------------------------
  ! Set up static calculations
  ! ------------------------------------------------------------
  do i=1,no_supercells
    sdir = wd//'/Supercell_'//i
    dft_input_filename = make_dft_input_filename(dft_code,seedname)
    call StructureData_to_dft_input_file(            &
       & dft_code,                                   &
       & supercells(i),                           &
       & wd//'/'//dft_input_filename, &
       & sdir//'/'//dft_input_filename)
  enddo
  
  ! ------------------------------------------------------------
  ! Generate quadratic configurations
  ! ------------------------------------------------------------
  do i=1,size(qpoints_ibz)
    supercell = supercells(qpoints_ibz(i)%sc_id)
    
    do j=1,structure%no_modes
      ! Read in mode.
      mode = read_normal_mode_file( &
         & harmonic_path//'/qpoint_'//i//'/mode_'//j//'.dat')
      
      ! Skip acoustic modes.
      if (mode%frequency <= frequency_tol) then
        cycle
      endif
      
      call mkdir(wd//'/qpoint_'//i//'/mode_'//j)
      
      do k=-no_samples,no_samples
        ! Skip equilibrium configurations
        if (k==0) then
          cycle
        endif
    
        call mkdir(wd//'/qpoint_'//i//'/mode_'//j//'/amplitude_'//k)
        
        ! Calculate quad_amplitude
        ! The normal mode amplitudes to sample are calculated
        if (temperature<temperature_tol) then
          ! Here quad_amplitude is really the normal mode amplitude squared,
          !   with units of 1/(E_h)
          occupation = 0.d0
        else
          occupation = 1.d0/(exp(mode%frequency/thermal_energy)-1.d0)
        endif
        quad_amplitude = sqrt((occupation+0.5d0)/mode%frequency)
        
        ! Calculate amplitude
        amplitude = quad_amplitude*displacement*k/real(no_samples,dp)
        
        ! Calculate new positions
        do l=1,supercell%no_atoms
          rvec = supercell%atom_to_rvec(l)
          prim = supercell%atom_to_prim(l)
          qr = qpoints_ibz(i)%qpoint*supercell%rvectors(rvec)*2*pi
          disp = amplitude * real( mode%displacements(prim) &
                               & * cmplx(cos(qr),sin(qr),dp))
          supercell%atoms(l) = supercell%atoms(l) + disp
        enddo
        
        ! Write dft input files
        mdir = wd//'/qpoint_'//i//'/mode_'//j//'/amplitude_'//k
        
        dft_input_filename = make_dft_input_filename(dft_code,seedname)
        call StructureData_to_dft_input_file( &
           & dft_code,                        &
           & supercell,                    &
           & wd//'/'//dft_input_filename,     &
           & mdir//'/'//dft_input_filename)
      enddo
    enddo
  enddo
end subroutine
end module
