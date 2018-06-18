! ======================================================================
! Calculates anharmonic states using the anharmonic potential and VSCF.
! ======================================================================
module calculate_states_module
  use common_module
  
  use anharmonic_common_module
  use polynomial_module
  
  use setup_harmonic_module
  
  use setup_anharmonic_module
  implicit none
  
  private
  
  public :: calculate_states
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function calculate_states() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'calculate_states'
  output%description = 'Runs VSCF on the anharmonic potential to calculate &
     &states. Should be run after calculate_potential.'
  output%keywords = [                                                         &
  & KeywordData( 'no_single_mode_samples',                                    &
  &              'no_single_mode_samples is the number of points (either side &
  &of zero) along each mode at which the anharmonic potential will be sampled &
  &when determining the effective frequency with which the harmonic basis &
  &along that mode will be constructed.',                                     &
  &              default_value='100'),                                        &
  & KeywordData( 'validate_potential',                                        &
  &              'validate_potential specifies that the anharmonic potential &
  &should be verified against fresh electronic structure calculations when &
  &sampling. Depending on the electronic structure method, this is likely to &
  &be very computationally intensive.',                                       &
  &              default_value='false') ]
  output%main_subroutine => calculate_states_subroutine
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_states_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! Input arguments.
  integer :: no_single_mode_samples
  logical :: validate_potential
  
  ! Arguments to setup_anharmonic.
  type(Dictionary) :: setup_anharmonic_arguments
  type(String)     :: harmonic_path
  type(String)     :: potential_representation
  real(dp)         :: maximum_displacement
  real(dp)         :: frequency_of_max_displacement
  
  ! Arguments to setup_harmonic.
  type(Dictionary) :: setup_harmonic_arguments
  real(dp)         :: symmetry_precision
  
  ! Maximum displacement in mass-weighted co-ordinates.
  real(dp) :: maximum_weighted_displacement
  
  ! Primitive structure.
  type(StructureData) :: structure
  
  ! Large anharmonic supercell and its q-points.
  type(StructureData)           :: anharmonic_supercell
  type(QpointData), allocatable :: qpoints(:)
  
  ! Normal modes.
  type(ComplexMode), allocatable :: complex_modes(:)
  type(RealMode),    allocatable :: real_modes(:)
  
  ! Anharmonic potential.
  type(PotentialPointer) :: potential
  
  ! Variables for generating effective mode frequencies.
  real(dp),                 allocatable :: displacements(:)
  real(dp),                 allocatable :: scaled_displacements(:)
  type(EffectiveFrequency), allocatable :: effective_frequencies(:)
  integer,                  allocatable :: qpoint_modes(:)
  
  ! Files and directories.
  type(IFile)  :: qpoints_file
  type(IFile)  :: complex_modes_file
  type(IFile)  :: real_modes_file
  type(IFile)  :: potential_file
  type(String) :: qpoint_dir
  type(OFile)  :: effective_frequencies_file
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  ! --------------------------------------------------
  ! Read in inputs and previously calculated data.
  ! --------------------------------------------------
  
  wd = arguments%value('working_directory')
  
  no_single_mode_samples = int(arguments%value('no_single_mode_samples'))
  validate_potential = lgcl(arguments%value('validate_potential'))
  
  ! Read in setup_anharmonic arguments.
  setup_anharmonic_arguments = Dictionary(setup_anharmonic())
  call setup_anharmonic_arguments%read_file( &
     & wd//'/setup_anharmonic.used_settings')
  harmonic_path = setup_anharmonic_arguments%value('harmonic_path')
  potential_representation = &
     & setup_anharmonic_arguments%value('potential_representation')
  maximum_displacement = &
     & dble(setup_anharmonic_arguments%value('maximum_displacement'))
  frequency_of_max_displacement = &
     & dble(setup_anharmonic_arguments%value('frequency_of_max_displacement'))
  
  ! Read in setup_harmonic arguments.
  setup_harmonic_arguments = Dictionary(setup_harmonic())
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  symmetry_precision = &
     & dble(setup_harmonic_arguments%value('symmetry_precision'))
  
  ! Read in structure.
  structure = read_structure_file( harmonic_path//'/structure.dat', &
                                 & symmetry_precision)
  
  ! Read in large anharmonic supercell and its q-points.
  anharmonic_supercell = read_structure_file( &
           & wd//'/anharmonic_supercell.dat', &
           & symmetry_precision,              &
           & calculate_symmetry=.false.)
  
  qpoints_file = IFile(wd//'/qpoints.dat')
  qpoints = QpointData(qpoints_file%sections())
  
  ! Read in normal modes.
  complex_modes_file = IFile(wd//'/complex_modes.dat')
  complex_modes = ComplexMode(complex_modes_file%sections())
  
  real_modes_file = IFile(wd//'/real_modes.dat')
  real_modes = RealMode(real_modes_file%sections())
  
  ! Read in anharmonic potential.
  potential_file = IFile(wd//'/potential.dat')
  if (potential_representation=='polynomial') then
    potential = PolynomialPotential(StringArray(potential_file%lines()))
  else
    call print_line( ERROR//': Unrecognised potential representation : '// &
                   & potential_representation)
    call err()
  endif
  
  ! --------------------------------------------------
  ! Re-calculate maximum_weighted_displacement.
  ! --------------------------------------------------
  maximum_weighted_displacement = maximum_displacement &
                              & * sqrt(minval(structure%atoms%mass()))
  
  ! --------------------------------------------------
  ! Calculate effective harmonic potential, from which initial harmonic
  !    states are constructed.
  ! --------------------------------------------------
  ! Calculate displacements before scaling by 1/sqrt(frequency).
  displacements =                                               &
     &   [(j,j=-no_single_mode_samples,no_single_mode_samples)] &
     & * maximum_displacement                                   &
     & / no_single_mode_samples
  
  allocate( effective_frequencies(size(complex_modes)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(complex_modes)
    ! Scale displacement by 1/sqrt(frequency).
    scaled_displacements = displacements                          &
                       & * sqrt( frequency_of_max_displacement    &
                       &       / max( complex_modes(i)%frequency, &
                       &              frequency_of_max_displacement))
    effective_frequencies(i) = EffectiveFrequency( scaled_displacements, &
                                                 & complex_modes(i),     &
                                                 & real_modes,           &
                                                 & potential)
  enddo
  
  ! --------------------------------------------------
  ! Write effective frequencies to file, q-point by q-point.
  ! --------------------------------------------------
  do i=1,size(qpoints)
    qpoint_dir = wd//'/qpoint_'//left_pad(i,str(size(qpoints)))
    call mkdir(qpoint_dir)
    
    effective_frequencies_file = &
       & OFile(qpoint_dir//'/effective_frequencies.dat')
    qpoint_modes = filter(complex_modes%qpoint_id==qpoints(i)%id)
    call effective_frequencies_file%print_lines( &
          & effective_frequencies(qpoint_modes), &
          & separating_line='')
  enddo
  
end subroutine
end module
