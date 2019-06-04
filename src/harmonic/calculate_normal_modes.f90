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
  
  ! Result type for calculate_dynamical_matrices.
  type :: MatrixAndModes
    type(DynamicalMatrix)          :: matrix
    type(ComplexMode), allocatable :: modes(:)
  end type
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
  supercell_hessians = construct_supercell_hessian( supercells,               &
                                                  & supercell_directories,    &
                                                  & calculation_reader,       &
                                                  & acoustic_sum_rule_forces, &
                                                  & hessian_logfile           )
  
  ! Calculate the dynamical matrix and normal modes at each q-point.
  dynamical_matrix_logfile = OFile('dynamical_matrix_log.dat')
  matrices_and_modes = calculate_dynamical_matrices( structure,               &
                                                   & supercells,              &
                                                   & supercell_hessians,      &
                                                   & qpoints,                 &
                                                   & dynamical_matrix_logfile )
  dynamical_matrices = matrices_and_modes%matrix
  
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
    call complex_modes_file%print_lines( matrices_and_modes(i)%modes, &
                                       & separating_line=''           )
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

! Read in the calculated forces, and use them to construct the Hessian.
impure elemental function construct_supercell_hessian(supercell,directory, &
   & calculation_reader,acoustic_sum_rule_forces,logfile) result(output)
  implicit none
  
  type(StructureData),     intent(in)    :: supercell
  type(String),            intent(in)    :: directory
  type(CalculationReader), intent(inout) :: calculation_reader
  logical,                 intent(in)    :: acoustic_sum_rule_forces
  type(OFile),             intent(inout) :: logfile
  type(CartesianHessian)                 :: output
  
  type(UniqueDirection),     allocatable :: unique_directions(:)
  type(ElectronicStructure), allocatable :: electronic_structure(:)
  
  ! Files and directories.
  type(IFile)  :: unique_directions_file
  type(String) :: calculation_directory
  
  ! Temporary variables.
  integer :: i,ialloc
  
  ! Read in symmetry group and unique atoms.
  unique_directions_file = IFile(directory//'/unique_directions.dat')
  unique_directions = UniqueDirection(unique_directions_file%sections())
  
  ! Read in electronic structure.
  allocate( electronic_structure(size(unique_directions)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(unique_directions)
    calculation_directory =                                    &
       & directory//'/atom.'                                // &
       & left_pad( unique_directions(i)%atom_id,               &
       &           str(maxval(unique_directions%atom_id)) ) // &
       & '.'//unique_directions(i)%direction
    
    electronic_structure(i) = calculation_reader%read_calculation( &
           & calculation_directory,                                &
           & CartesianDisplacement(unique_directions(i),supercell) )
    
    if (    size(electronic_structure(i)%forces%vectors) &
       & /= supercell%no_atoms                           ) then
      call print_line( ERROR//': Wrong number of forces in '//            &
                     & calculation_directory//'/electronic_structure.dat' )
      call err()
    endif
  enddo
  
  ! Calculate Hessian from the set of forces at each displacement.
  output = CartesianHessian( supercell,                &
                           & unique_directions,        &
                           & electronic_structure,     &
                           & acoustic_sum_rule_forces, &
                           & logfile                   )
end function

! Calculate dynamical matrices at each q-point from the Hessians of each
!    non-diagonal supercell.
function calculate_dynamical_matrices(structure,supercells, &
   & supercell_hessians,qpoints,logfile) result(output)
  implicit none
  
  type(StructureData),    intent(in)    :: structure
  type(StructureData),    intent(in)    :: supercells(:)
  type(CartesianHessian), intent(in)    :: supercell_hessians(:)
  type(QpointData),       intent(in)    :: qpoints(:)
  type(OFile),            intent(inout) :: logfile
  type(MatrixAndModes), allocatable     :: output(:)
  
  logical, allocatable :: modes_calculated(:)
  
  integer :: subspace_id
  
  type(QpointData)      :: transformed_qpoint
  type(DynamicalMatrix) :: transformed_matrix
  
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Calculate the dynamical matrix and normal modes at each q-point.
  ! --------------------------------------------------
  allocate( modes_calculated(size(qpoints)), &
          & output(size(qpoints)),           &
          & stat=ialloc); call err(ialloc)
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
        call logfile%print_line('Constructing dynamical matrix and &
           &normal modes at q-point '//i//' as the conjugate of those at &
           &q-point '//j)
        
        ! The dynamical matrix at G-q is the complex conjugate of that at q.
        ! N.B. the complex conjugate, not the Hermitian conjugate.
        output(i)%matrix = DynamicalMatrix(                     &
           & matrices      = conjg(output(j)%matrix%matrices()) )
        
        ! The displacements of the normal modes at G-q are
        !    the complex conjugates of those at q.
        output(i)%modes = conjg(output(j)%modes)
        
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
            call logfile%print_line('Constructing dynamical matrix and &
               &normal modes at q-point '//i//' using symmetry from those at &
               &q-point '//j)
            output(i)%matrix = transform_modes( output(j)%matrix,        &
                                              & structure%symmetries(k), &
                                              & qpoints(j),              &
                                              & qpoints(i)               )
            output(i)%modes = transform( output(j)%modes,         &
                                       & structure%symmetries(k), &
                                       & qpoints(j),              &
                                       & qpoints(i)               )
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
      
      do j=1,size(supercells)
        ! Check if q-point i is a G-vector of supercell j.
        if (is_int(supercells(j)%supercell * qpoints(i)%qpoint)) then
          call logfile%print_line('Constructing dynamical matrix and &
             &normal modes at q-point '//i//' directly from calculated force &
             &constants.')
          output(i) = construct_dynamical_matrix( qpoints(i),         &
                                                & supercells,         &
                                                & supercell_hessians, &
                                                & structure,          &
                                                & subspace_id,        &
                                                & logfile             )
          subspace_id = maxval(output(i)%modes%subspace_id) + 1
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
    call output(i)%matrix%check(structure, logfile)
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
          transformed_matrix = transform_modes( output(j)%matrix,        &
                                              & structure%symmetries(i), &
                                              & qpoints(j),              &
                                              & qpoints(k)               )
          call logfile%print_line('Comparing symmetrically &
             &equivalent dynamical matrices.')
          call compare_dynamical_matrices( output(k)%matrix,   &
                                         & transformed_matrix, &
                                         & logfile             )
        endif
      enddo
    enddo
  enddo
end function

! --------------------------------------------------
! Construct the dynamical matrix at the specified q-point from the given
!    Hessian.
! Considers all supercell, S, and symmetry, R, combinations such that
!    S.R.q is a vector of integers, i.e. the q-point is a G-vector of the
!    supercell.
! --------------------------------------------------
function construct_dynamical_matrix(qpoint,supercells,hessian, &
   & structure,subspace_id,logfile) result(output)
  implicit none
  
  type(QpointData),       intent(in)    :: qpoint
  type(StructureData),    intent(in)    :: supercells(:)
  type(CartesianHessian), intent(in)    :: hessian(:)
  type(StructureData),    intent(in)    :: structure
  integer,                intent(in)    :: subspace_id
  type(OFile),            intent(inout) :: logfile
  type(MatrixAndModes)                  :: output
  
  type(QpointData) :: q_prime
  
  logical, allocatable :: is_copy(:,:)
  integer, allocatable :: sym_id(:)
  integer, allocatable :: sup_id(:)
  integer              :: no_copies
  integer              :: i,j,k,ialloc
  
  type(DynamicalMatrix)            :: dynamical_matrix
  type(ComplexMatrix), allocatable :: copy_matrices(:,:,:)
  type(ComplexMatrix), allocatable :: matrices(:,:)
  type(ComplexMode),   allocatable :: complex_modes(:)
  
  logical                          :: check_matrices
  type(ComplexMatrix), allocatable :: average(:,:)
  real(dp),            allocatable :: averages(:,:)
  real(dp),            allocatable :: differences(:,:)
  real(dp)                         :: l2_error
  
  ! Check that the supercells and Hessians correspond to one another.
  if (size(supercells)/=size(hessian)) then
    call print_line(CODE_ERROR//': Force constants and supercells do not &
       &match.')
    call err()
  endif
  
  ! Identify how many parings of supercell and symmetry
  !    allow q to be simulated.
  allocate( is_copy(size(supercells),size(structure%symmetries)), &
          & stat=ialloc); call err(ialloc)
  is_copy = .false.
  do i=1,size(structure%symmetries)
    ! q' is defined such that symmetry i maps qp onto q.
    ! N.B. q-points transform as q->R^-T.q, so q' s.t. q'->q is given by
    !    q' = R^T.q = q.R.
    q_prime = structure%symmetries(i)%inverse_transform(qpoint)
    do j=1,size(supercells)
      if (is_int(supercells(j)%supercell*q_prime%qpoint)) then
        is_copy(j,i) = .true.
      endif
    enddo
  enddo
  no_copies = count(is_copy)
  
  ! Construct a copy of the dynamical matrix from each pair of supercell and
  !    symmetry.
  allocate( copy_matrices(structure%no_atoms,structure%no_atoms,no_copies), &
          & sym_id(no_copies),                                              &
          & sup_id(no_copies),                                              &
          & stat=ialloc); call err(ialloc)
  k = 0
  do i=1,size(structure%symmetries)
    q_prime = structure%symmetries(i)%inverse_transform(qpoint)
    do j=1,size(supercells)
      if (is_copy(j,i)) then
        k = k+1
        dynamical_matrix = DynamicalMatrix( dblevec(q_prime%qpoint), &
                                          & supercells(j),           &
                                          & hessian(j)               )
        dynamical_matrix = transform_modes( dynamical_matrix,        &
                                          & structure%symmetries(i), &
                                          & q_prime,                 &
                                          & qpoint                   )
        copy_matrices(:,:,k) = dynamical_matrix%matrices()
        sym_id(k) = i
        sup_id(k) = j
      endif
    enddo
  enddo
  
  ! Work out if the output matrices should be equivalent.
  ! If the calculation is coming from more than one supercell, or from a
  !    supercell which is larger than necessary to simulate the q-point,
  !    then the separate contributions to the dynamical matrix will differ.
  check_matrices = .true.
  if (count(any(is_copy,2))>1) then
    check_matrices = .false.
  endif
  do i=1,size(supercells)
    if (any(is_copy(i,:))) then
      if (supercells(i)%sc_size/=qpoint%min_sc_size()) then
        check_matrices = .false.
      endif
    endif
  enddo
  
  if (check_matrices) then
    allocate( average(structure%no_atoms,structure%no_atoms),     &
            & averages(structure%no_atoms,structure%no_atoms),    &
            & differences(structure%no_atoms,structure%no_atoms), &
            & stat=ialloc); call err(ialloc)
    average = cmplxmat(zeroes(3,3))
    averages = 0
    differences = 0
    do i=1,no_copies
      average = average + copy_matrices(:,:,i)/no_copies
      do j=1,i-1
        averages = averages &
               & + sum_squares((copy_matrices(:,:,j)+copy_matrices(:,:,i))/2)
        differences = differences &
                  & + sum_squares(copy_matrices(:,:,j)-copy_matrices(:,:,i))
      enddo
    enddo
    l2_error = sqrt(sum(differences)/sum(averages))
    call logfile%print_line('Fractional L2 difference between equivalent &
       &dynamical matrices: '//l2_error)
    if (l2_error>1e-10_dp) then
      call print_line(WARNING//': Symmetrically equivalent dynamical &
         &matrices differ. Please check log files.')
      call print_line('Fractional L2 error: '//l2_error)
    endif
  endif
  
  ! Average over the copies.
  allocate( matrices(structure%no_atoms,structure%no_atoms), &
          & stat=ialloc); call err(ialloc)
  matrices = cmplxmat(zeroes(3,3))
  do i=1,no_copies
    matrices = matrices + copy_matrices(:,:,i)/no_copies
  enddo
  
  dynamical_matrix = DynamicalMatrix(matrices)
  
  ! Diagonalise the dynamical matrix, to obtain the normal mode
  !    co-ordinates (eigenvectors) and harmonic frequencies (eigenvalues).
  complex_modes = ComplexMode( dynamical_matrix,                    &
                             & structure,                           &
                             & modes_real=qpoint%is_paired_qpoint() )
  
  ! Process the modes to calculate ids.
  complex_modes = process_modes(complex_modes, structure, qpoint, subspace_id)
  
  output = MatrixAndModes(dynamical_matrix, complex_modes)
end function

! ----------------------------------------------------------------------
! Compares two dynamical matrices.
! ----------------------------------------------------------------------
subroutine compare_dynamical_matrices(a,b,logfile)
  implicit none
  
  type(DynamicalMatrix), intent(in)    :: a
  type(DynamicalMatrix), intent(in)    :: b
  type(OFile),           intent(inout) :: logfile
  
  integer :: no_atoms
  
  type(ComplexMatrix), allocatable :: mats_a(:,:)
  type(ComplexMatrix), allocatable :: mats_b(:,:)
  
  type(ComplexMatrix) :: mat_a
  type(ComplexMatrix) :: mat_b
  
  real(dp) :: average
  real(dp) :: difference
  
  integer :: i,j
  
  mats_a = a%matrices()
  mats_b = b%matrices()
  
  no_atoms = size(mats_a,1)
  if (size(mats_b,1)/=no_atoms) then
    call print_line(CODE_ERROR//': dynamical matrices a and b have different &
       &sizes.')
    call err()
  endif
  
  average = 0.0_dp
  difference = 0.0_dp
  do i=1,no_atoms
    do j=1,no_atoms
      mat_a = mats_a(j,i)
      mat_b = mats_b(j,i)
      average = average + sum_squares((mat_a+mat_b)/2)
      difference = difference + sum_squares(mat_a-mat_b)
    enddo
  enddo
  
  if (average>1e-30_dp) then
    call logfile%print_line('Fractional L2 difference between dynamical &
       &matrices: '//sqrt(difference/average))
    if (sqrt(difference/average)>1e-10_dp) then
      call print_line(WARNING//': Dynamical matrices differ. Please check &
         &log files.')
      call print_line('Fractional L2 difference between dynamical &
         &matrices: '//sqrt(difference/average))
    endif
  else
    call logfile%print_line('Dynamical matrices too small to compare &
       &fractional differences.')
  endif
end subroutine

! ----------------------------------------------------------------------
! Transform a dynamical matrix and set of normal modes
!    from one q-point to another.
! ----------------------------------------------------------------------
! Construct data at q_new from data at q_old, where
!    R . q_old = q_new.
! N.B. the provided q-point should be q_new not q_old.
! The symmetry, S, maps equilibrium position ri to rj+R, and q-point q to q'.
! The +R needs to be corrected for, by multiplying the phase by exp(-iq'.R).
function transform_modes(input,symmetry,qpoint_from,qpoint_to) result(output)
  implicit none
  
  type(DynamicalMatrix),  intent(in)    :: input
  type(SymmetryOperator), intent(in)    :: symmetry
  type(QpointData),       intent(in)    :: qpoint_from
  type(QpointData),       intent(in)    :: qpoint_to
  type(DynamicalMatrix)                 :: output
  
  output = DynamicalMatrix(                                            &
     & matrices      = transform_dynamical_matrix( input%matrices(),   &
     &                                             symmetry,           &
     &                                             qpoint_from,        &
     &                                             qpoint_to         ) )
end function

function transform_dynamical_matrix(input,symmetry,qpoint_from,qpoint_to) &
   & result(output)
  implicit none
  
  type(ComplexMatrix),    intent(in) :: input(:,:)
  type(SymmetryOperator), intent(in) :: symmetry
  type(QpointData),       intent(in) :: qpoint_from
  type(QpointData),       intent(in) :: qpoint_to
  type(ComplexMatrix), allocatable   :: output(:,:)
  
  integer :: no_atoms
  
  type(FractionVector) :: q
  type(IntVector)      :: r1
  type(IntVector)      :: r2
  
  integer :: atom_1
  integer :: atom_1p
  integer :: atom_2
  integer :: atom_2p
  
  integer :: ialloc
  
  no_atoms = size(input,1)
  q = qpoint_to%qpoint
  
  ! Check that the symmetry transforms the q-point as expected.
  if (symmetry * qpoint_from /= qpoint_to) then
    call print_line(CODE_ERROR//': Symmetry does not transform q-points as &
       &expected.')
    call err()
  endif
  
  ! Check that the number of atoms is consistent.
  if (size(input,2)/=no_atoms) then
    call err()
  endif
  
  ! Transform dynamical matrix.
  allocate(output(no_atoms,no_atoms), stat=ialloc); call err(ialloc)
  do atom_1=1,no_atoms
    atom_1p = symmetry%atom_group * atom_1
    r1 = symmetry%rvectors(atom_1)
    do atom_2=1,no_atoms
      atom_2p = symmetry%atom_group * atom_2
      r2 = symmetry%rvectors(atom_2)
      
      output(atom_2p,atom_1p) =      &
         & symmetry%cartesian_tensor &
         & * exp_2pii(-q*r2)         &
         & * input(atom_2,atom_1)    &
         & * exp_2pii(q*r1)          &
         & * transpose(symmetry%cartesian_tensor)
    enddo
  enddo
end function
end module
