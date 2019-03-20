! ======================================================================
! The third stage of Caesar.
! Uses the forces calculated previously to calculate harmonic normal modes.
! ======================================================================
! Calculates the set of phonon normal modes at each supercell G-vector.
! Calculates the real part of the non-mass-reduced polarisation vector, which
!    is the pattern of displacement corresponding to the normal mode.
module calculate_normal_modes_module
  use common_module
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
  
  ! Hessian matrix data.
  type(UniqueDirection),     allocatable :: unique_directions(:)
  type(ElectronicStructure), allocatable :: electronic_structure(:)
  type(CartesianHessian),    allocatable :: hessian(:)
  
  ! q-point data.
  type(StructureData)           :: large_supercell
  type(QpointData), allocatable :: qpoints(:)
  
  ! Lte output data.
  logical, allocatable :: modes_calculated(:)
  
  ! Normal modes and their symmetries.
  type(QpointData) :: transformed_qpoint
  
  type(DynamicalMatrix), allocatable :: dynamical_matrices(:)
  type(DynamicalMatrix)              :: transformed_matrix
  
  integer :: subspace_id
  
  type(ComplexMode), allocatable :: complex_modes(:,:)
  
  type(CartesianHessian) :: supercell_hessian
  
  ! Files and directories.
  type(IFile)  :: structure_file
  type(IFile)  :: large_supercell_file
  type(IFile)  :: qpoints_file
  type(IFile)  :: no_supercells_file
  type(IFile)  :: supercell_file
  type(IFile)  :: unique_directions_file
  type(String) :: supercell_dir
  type(String) :: calculation_dir
  type(String) :: qpoint_dir
  type(OFile)  :: dynamical_matrix_file
  type(OFile)  :: complex_modes_file
  type(OFile)  :: force_logfile
  type(OFile)  :: qpoint_logfile
  type(OFile)  :: castep_phonon_file
  type(OFile)  :: qe_force_constants_file
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Read in arguments from user.
  ! --------------------------------------------------
  acoustic_sum_rule = arguments%value('acoustic_sum_rule')
  if (acoustic_sum_rule=='off') then
    acoustic_sum_rule_forces = .false.
    acoustic_sum_rule_matrices = .false.
  elseif (acoustic_sum_rule=='forces') then
    acoustic_sum_rule_forces = .true.
    acoustic_sum_rule_matrices = .false.
  elseif (acoustic_sum_rule=='matrices') then
    acoustic_sum_rule_forces = .false.
    acoustic_sum_rule_matrices = .true.
  elseif (acoustic_sum_rule=='both') then
    acoustic_sum_rule_forces = .true.
    acoustic_sum_rule_matrices = .true.
  else
    call print_line(ERROR//': acoustic_sum_rule has been set to an unexpected &
       &value. Please set acoustic_sum_rule to "off", "forces", "matrices" or &
       &"both".')
    call quit()
  endif
  
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
  ! Calculate the Hessian matrix corresponding to each supercell.
  ! --------------------------------------------------
  allocate( supercells(no_supercells), &
          & hessian(no_supercells),    &
          & stat=ialloc); call err(ialloc)
  force_logfile = OFile('hessian_log.dat')
  do i=1,no_supercells
    supercell_dir = 'Supercell_'//left_pad(i,str(no_supercells))
    
    ! Read in supercell structure data.
    supercell_file = IFile(supercell_dir//'/structure.dat')
    supercells(i) = StructureData(supercell_file%lines())
    
    ! Read in symmetry group and unique atoms.
    unique_directions_file = IFile(supercell_dir//'/unique_directions.dat')
    unique_directions = UniqueDirection(unique_directions_file%sections())
    
    ! Read in electronic structure.
    allocate( electronic_structure(size(unique_directions)), &
            & stat=ialloc); call err(ialloc)
    do j=1,size(unique_directions)
      calculation_dir = supercell_dir//'/atom.'                            // &
                      & left_pad( unique_directions(j)%atom_id,               &
                      &           str(maxval(unique_directions%atom_id)) ) // &
                      & '.'//unique_directions(j)%direction
      electronic_structure(j) = calculation_reader%read_calculation( &
         & calculation_dir,                                          &
         & CartesianDisplacement(unique_directions(j),supercells(i)) )
      if ( size(electronic_structure(j)%forces%vectors) /= &
         & supercells(i)%no_atoms                          ) then
        call print_line( ERROR//': Wrong number of forces in '//      &
                       & calculation_dir//'/electronic_structure.dat' )
        call err()
      endif
    enddo
    
    ! Calculate Hessian from the set of forces at each displacement.
    hessian(i) = CartesianHessian( supercells(i),            &
                                 & unique_directions,        &
                                 & electronic_structure,     &
                                 & acoustic_sum_rule_forces, &
                                 & force_logfile             )
    deallocate( unique_directions,    &
              & electronic_structure, &
              & stat=ialloc); call err(ialloc)
  enddo
  
  ! --------------------------------------------------
  ! Calculate the dynamical matrix and normal modes at each q-points.
  ! --------------------------------------------------
  allocate( modes_calculated(size(qpoints)),   &
          & dynamical_matrices(size(qpoints)), &
          & stat=ialloc); call err(ialloc)
  qpoint_logfile = OFile('dynamical_matrix_log.dat')
  modes_calculated = .false.
  
  subspace_id = 1
  iter : do while (.not. all(modes_calculated))
    ! ------------------------------
    ! First, check if there is an un-calculated q-point whose conjugate has
    !    already been calculated.
    ! ------------------------------
    do i=1,size(qpoints)
      if (modes_calculated(i)) then
        cycle
      endif
      
      j = first(qpoints%id==qpoints(i)%paired_qpoint_id)
      
      if (modes_calculated(j)) then
        call qpoint_logfile%print_line('Constructing dynamical matrix and &
           &normal modes at q-point '//i//' as the conjugate of those at &
           &q-point '//j)
        dynamical_matrices(i) = conjg(dynamical_matrices(j))
        modes_calculated(i) = .true.
        cycle iter
      endif
    enddo
    
    ! ------------------------------
    ! Second, check if there is a q-point which is related by symmetry to
    !    an already calculated q-point.
    ! ------------------------------
    do i=1,size(qpoints)
      if (modes_calculated(i)) then
        cycle
      endif
      
      do j=1,size(qpoints)
        if (.not. modes_calculated(j)) then
          cycle
        endif
        
        do k=1,size(structure%symmetries)
          transformed_qpoint = structure%symmetries(k) * qpoints(j)
          if (transformed_qpoint == qpoints(i)) then
            call qpoint_logfile%print_line('Constructing dynamical matrix and &
               &normal modes at q-point '//i//' using symmetry from those at &
               &q-point '//j)
            dynamical_matrices(i) = transform_modes( dynamical_matrices(j),   &
                                                   & structure%symmetries(k), &
                                                   & qpoints(j),              &
                                                   & qpoints(i))
            modes_calculated(i) = .true.
            cycle iter
          endif
        enddo
      enddo
    enddo
    
    ! ------------------------------
    ! Finally check if there is a q-point which is a G-vector of any of the
    !    calculated supercells.
    ! ------------------------------
    do i=1,size(qpoints)
      if (modes_calculated(i)) then
        cycle
      endif
      
      do j=1,no_supercells
        ! Check if q-point i is a G-vector of supercell j.
        if (is_int(supercells(j)%supercell * qpoints(i)%qpoint)) then
          call qpoint_logfile%print_line('Constructing dynamical matrix and &
             &normal modes at q-point '//i//' directly from calculated force &
             &constants.')
          dynamical_matrices(i) = DynamicalMatrix( qpoints(i),    &
                                                 & supercells,    &
                                                 & hessian,       &
                                                 & structure,     &
                                                 & subspace_id,   &
                                                 & qpoint_logfile )
          subspace_id =                                                  &
             &   maxval(dynamical_matrices(i)%complex_modes%subspace_id) &
             & + 1
          modes_calculated(i) = .true.
          cycle iter
        endif
      enddo
    enddo
    
    ! ------------------------------
    ! If no q-point can be constructed, throw an error.
    ! ------------------------------
    call print_line(ERROR//': Insufficient data to construct dynamical &
       &matrices at all q-points.')
    call err()
  enddo iter
  
  ! --------------------------------------------------
  ! Check all dynamical matrices.
  ! --------------------------------------------------
  ! Run basic checks on each matrix in turn.
  do i=1,size(qpoints)
    call dynamical_matrices(i)%check(structure, qpoint_logfile)
  enddo
  
  ! Check that dynamical matrices at q-points q and q' s.t. qS=q'
  !    correctly transform onto one another.
  do i=1,size(structure%symmetries)
    do j=1,size(qpoints)
      do k=1,size(qpoints)
        transformed_qpoint = structure%symmetries(i) * qpoints(j)
        if (transformed_qpoint == qpoints(k)) then
          if (qpoints(j)%paired_qpoint_id/=k) then
            cycle
          endif
          transformed_matrix = transform_modes( dynamical_matrices(j),   &
                                              & structure%symmetries(i), &
                                              & qpoints(j),              &
                                              & qpoints(k))
          call qpoint_logfile%print_line('Comparing symmetrically &
             &equivalent dynamical matrices.')
          call compare_dynamical_matrices( dynamical_matrices(k), &
                                         & transformed_matrix,    &
                                         & qpoint_logfile)
        endif
      enddo
    enddo
  enddo
  
  ! --------------------------------------------------
  ! Write out output.
  ! --------------------------------------------------
  allocate( complex_modes(structure%no_modes,size(qpoints)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(qpoints)
    qpoint_dir = 'qpoint_'//left_pad(i,str(size(qpoints)))
    call mkdir(qpoint_dir)
    
    ! Write out dynamical matrix.
    dynamical_matrix_file = OFile(qpoint_dir//'/dynamical_matrix.dat')
    call dynamical_matrix_file%print_lines(dynamical_matrices(i))
    
    ! Write out normal modes.
    complex_modes(:,i) = dynamical_matrices(i)%complex_modes
    complex_modes_file = OFile(qpoint_dir//'/complex_modes.dat')
    call complex_modes_file%print_lines( dynamical_matrices(i)%complex_modes, &
                                       & separating_line='')
  enddo
  
  ! Write out Castep .phonon file.
  castep_phonon_file = OFile(make_castep_phonon_filename(seedname))
  call write_castep_phonon_file( castep_phonon_file, &
                               & complex_modes,      &
                               & qpoints,            &
                               & structure           )
  
  ! Construct supercell Hessian, and write out QE force constants file.
  supercell_hessian = reconstruct_hessian( large_supercell,    &
                                         & qpoints,            &
                                         & dynamical_matrices, &
                                         & force_logfile       )
  qe_force_constants_file = OFile(make_qe_force_constants_filename(seedname))
  call write_qe_force_constants_file( qe_force_constants_file, &
                                    & supercell_hessian,       &
                                    & structure,               &
                                    & large_supercell          )
end subroutine
end module
