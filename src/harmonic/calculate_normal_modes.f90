! ======================================================================
! The third stage of Caesar.
! Uses the forces calculated previously to calculate harmonic normal modes.
! ======================================================================
! Calculates the set of phonon normal modes at each supercell G-vector.
! Calculates the real part of the non-mass-reduced polarisation vector, which
!    is the pattern of displacement corresponding to the normal mode.
module calculate_normal_modes_module
  use common_module
  
  use construct_supercell_hessian_module
  use calculate_dynamical_matrices_module
  implicit none
  
  private
  
  public :: startup_calculate_normal_modes
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_calculate_normal_modes()
  implicit none
  
  type(CaesarMode) :: mode
  
  mode%mode_name = 'calculate_normal_modes'
  mode%description = 'Finds harmonic normal modes. Should be called &
     &after run_harmonic.'
  mode%keywords = [                                                           &
     & KeywordData( 'acoustic_sum_rule',                                      &
     &              'acoustic_sum_rule specifies where the acoustic sum rule &
     &is applied. The options are "off", "forces", "matrices" and "both".',   &
     &              default_value='both'),                                    &
     & KeywordData( 'loto_direction',                                         &
     &              'loto_direction specifies the direction (in reciprocal &
     &co-ordinates from which the gamma point is approached when calculating &
     &LO/TO corrections. See setup_harmonic help for details. loto_direction &
     &may only be specified here if it was not specified in setup_harmonic, &
     &and if every symmetry of the system leaves it invariant, i.e. q.S=q for &
     &all S, where S is the symmetry matrix and q is loto_direction. See &
     &structure.dat for the list of symmetries.',                             &
     &              is_optional = .true.)                                     ]
  mode%main_subroutine => calculate_normal_modes_subroutine
  
  call add_mode(mode)
end subroutine

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_normal_modes_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Input arguments.
  type(String)         :: acoustic_sum_rule
  logical              :: acoustic_sum_rule_forces
  logical              :: acoustic_sum_rule_matrices
  type(FractionVector) :: loto_direction
  logical              :: loto_direction_set
  
  ! Arguments to setup_harmonic.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: seedname
  
  ! Initial data calculated by setup_harmonic.
  integer                          :: no_supercells
  type(StructureData)              :: structure
  type(StructureData), allocatable :: supercells(:)
  
  ! Electronic structure calculation reader.
  type(CalculationReader) :: calculation_reader
  
  ! q-point data.
  type(StructureData)           :: large_supercell
  type(QpointData), allocatable :: qpoints(:)
  
  ! Hesians, dynamical matrices and normal modes.
  type(CartesianHessian), allocatable :: supercell_hessians(:)
  type(MatrixAndModes),   allocatable :: matrices_and_modes(:)
  type(DynamicalMatrix),  allocatable :: dynamical_matrices(:)
  type(ComplexMode),      allocatable :: complex_modes(:,:)
  type(CartesianHessian)              :: full_hessian
  
  ! Files and directories.
  type(IFile)               :: structure_file
  type(IFile)               :: large_supercell_file
  type(IFile)               :: qpoints_file
  type(IFile)               :: no_supercells_file
  type(IFile)               :: supercell_file
  type(String), allocatable :: supercell_directories(:)
  type(String)              :: qpoint_dir
  type(OFile)               :: dynamical_matrix_file
  type(OFile)               :: complex_modes_file
  type(OFile)               :: hessian_logfile
  type(OFile)               :: dynamical_matrix_logfile
  type(OFile)               :: castep_phonon_file
  type(OFile)               :: qe_force_constants_file
  
  ! Temporary variables.
  integer :: i,ialloc
  
  ! --------------------------------------------------
  ! Read in arguments from user.
  ! --------------------------------------------------
  acoustic_sum_rule = arguments%value('acoustic_sum_rule')
  if (.not. any(acoustic_sum_rule==tokens('off forces matrices both'))) then
    call print_line(ERROR//': acoustic_sum_rule has been set to an unexpected &
       &value. Please set acoustic_sum_rule to "off", "forces", "matrices" or &
       &"both".')
    call quit()
  endif
  acoustic_sum_rule_forces = any(acoustic_sum_rule==tokens('forces both'))
  acoustic_sum_rule_matrices = any(acoustic_sum_rule==tokens('matrices both'))
  
  ! --------------------------------------------------
  ! Read in previous arguments.
  ! --------------------------------------------------
  setup_harmonic_arguments = Dictionary(CaesarMode('setup_harmonic'))
  call setup_harmonic_arguments%read_file('setup_harmonic.used_settings')
  seedname = setup_harmonic_arguments%value('seedname')
  
  structure_file = IFile('structure.dat')
  structure = StructureData(structure_file%lines())
  
  large_supercell_file = IFile('large_supercell.dat')
  large_supercell = StructureData(large_supercell_file%lines())
  
  qpoints_file = IFile('qpoints.dat')
  qpoints = QpointData(qpoints_file%sections())
  
  no_supercells_file = IFile('no_supercells.dat')
  no_supercells = int(no_supercells_file%line(1))
  
  allocate( supercell_directories(no_supercells), &
          & supercells(no_supercells),            &
          & stat=ialloc); call err(ialloc)
  do i=1,no_supercells
    supercell_directories(i) = 'Supercell_'//left_pad(i,str(no_supercells))
    supercell_file = IFile(supercell_directories(i)//'/structure.dat')
    supercells(i) = StructureData(supercell_file%lines())
  enddo
  
  ! --------------------------------------------------
  ! Initialise LO/TO splitting if necessary.
  ! --------------------------------------------------
  loto_direction_set = .false.
  if (setup_harmonic_arguments%is_set('loto_direction')) then
    loto_direction = FractionVector(                      &
       & setup_harmonic_arguments%value('loto_direction') )
    loto_direction_set = .true.
  endif
  
  if (arguments%is_set('loto_direction')) then
    if (loto_direction_set) then
      call print_line(ERROR//': loto_direction may not be specified here, &
         &since it was already specified in setup_harmonic.')
      call quit()
    endif
    loto_direction = FractionVector(                      &
       & setup_harmonic_arguments%value('loto_direction') )
    loto_direction_set = .true.
    if (any(loto_breaks_symmetry( structure%symmetries%tensor, &
                                & loto_direction               ))) then
      call print_line(ERROR//': loto_direction has been specified in a &
         &direction which breaks symmetry. To specify this direction, please &
         &set loto_direction when running setup_harmonic.')
      call quit()
    endif
  endif
  
  ! --------------------------------------------------
  ! Initialise calculation reader.
  ! --------------------------------------------------
  if (loto_direction_set) then
    calculation_reader = CalculationReader(loto_direction)
  else
    calculation_reader = CalculationReader()
  endif
  
  ! --------------------------------------------------
  ! Calculate supercell Hessians, dynamical matrices, normal modes
  !    and the full Hessian.
  ! --------------------------------------------------
  ! Calculate the Hessian matrix corresponding to each non-diagonal supercell.
  hessian_logfile = OFile('hessian_log.dat')
  supercell_hessians = [(                                        &
     & construct_supercell_hessian( supercells(i),               &
     &                              supercell_directories(i),    &
     &                              calculation_reader,          &
     &                              acoustic_sum_rule_forces,    &
     &                              hessian_logfile           ), &
     & i=1,                                                      &
     & size(supercells)                                          )]
  
  ! Calculate the dynamical matrix and normal modes at each q-point.
  dynamical_matrix_logfile = OFile('dynamical_matrix_log.dat')
  matrices_and_modes = calculate_dynamical_matrices( structure,               &
                                                   & supercells,              &
                                                   & supercell_hessians,      &
                                                   & qpoints,                 &
                                                   & dynamical_matrix_logfile )
  dynamical_matrices = [( matrices_and_modes(i)%matrix, &
                        & i=1,                          &
                        & size(matrices_and_modes)      )]
  
  ! Copy normal modes from dynamical matrices into a single array.
  allocate( complex_modes(structure%no_modes,size(qpoints)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(qpoints)
    complex_modes(:,i) = matrices_and_modes(i)%modes
  enddo
  
  ! Construct the Hessian for the full harmonic supercell.
  full_hessian = reconstruct_hessian( large_supercell,    &
                                    & qpoints,            &
                                    & dynamical_matrices, &
                                    & hessian_logfile     )
  
  ! --------------------------------------------------
  ! Write out dynamical matrices, normal modes and the full Hessian.
  ! --------------------------------------------------
  do i=1,size(qpoints)
    qpoint_dir = 'qpoint_'//left_pad(i,str(size(qpoints)))
    call mkdir(qpoint_dir)
    
    ! Write out dynamical matrix.
    dynamical_matrix_file = OFile(qpoint_dir//'/dynamical_matrix.dat')
    call dynamical_matrix_file%print_lines(dynamical_matrices(i))
    
    ! Write out normal modes.
    complex_modes_file = OFile(qpoint_dir//'/complex_modes.dat')
    call complex_modes_file%print_lines( complex_modes(:,i), &
                                       & separating_line=''  )
  enddo
  
  ! Write out Castep .phonon file. This contains all normal modes.
  castep_phonon_file = OFile(make_castep_phonon_filename(seedname))
  call write_castep_phonon_file( castep_phonon_file, &
                               & complex_modes,      &
                               & qpoints,            &
                               & structure           )
  
  ! Write out QE force constants file. This contains the full Hessian.
  qe_force_constants_file = OFile(make_qe_force_constants_filename(seedname))
  call write_qe_force_constants_file( qe_force_constants_file, &
                                    & full_hessian,            &
                                    & structure,               &
                                    & large_supercell          )
end subroutine
end module
