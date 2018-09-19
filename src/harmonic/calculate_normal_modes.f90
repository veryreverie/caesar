! ======================================================================
! The third stage of Caesar.
! Uses the forces calculated previously to calculate harmonic normal modes.
! ======================================================================
! Calculates the set of phonon normal modes at each supercell G-vector.
! Calculates the real part of the non-mass-reduced polarisation vector, which
!    is the pattern of displacement corresponding to the normal mode.
module calculate_normal_modes_module
  use common_module
  
  use setup_harmonic_module
  implicit none
  
  private
  
  public :: calculate_normal_modes
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function calculate_normal_modes() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'calculate_normal_modes'
  output%description = 'Finds harmonic normal modes. Should be called &
     &after run_harmonic.'
  output%keywords = [                                                         &
     & KeywordData( 'acoustic_sum_rule',                                      &
     &              'acoustic_sum_rule specifies where the acoustic sum rule  &
     &is applied. The options are "off", "forces", "matrices" and "both".',   &
     &              default_value='both'),                                    &
     & KeywordData( 'degenerate_energy',                                      &
     &              'degenerate_energy is the minimum energy difference &
     &between states before they are considered degenerate. This should be &
     &given in Hartrees.',                                                    &
     &              default_value='1e-5')                                     ]
  output%main_subroutine => calculate_normal_modes_subroutine
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_normal_modes_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! Input arguments.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: acoustic_sum_rule
  logical          :: acoustic_sum_rule_forces
  logical          :: acoustic_sum_rule_matrices
  real(dp)         :: degenerate_energy
  
  ! Arguments to setup_harmonic.
  type(String) :: seedname
  
  ! Initial data calculated by setup_harmonic.
  integer                          :: no_supercells
  type(StructureData)              :: structure
  type(StructureData), allocatable :: supercells(:)
  
  ! Electronic structure calculation reader.
  type(CalculationReader) :: calculation_reader
  
  ! Force constant data.
  type(UniqueDirection), allocatable :: unique_directions(:)
  type(ForceConstants),  allocatable :: force_constants(:)
  
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
  
  ! Files and directories.
  type(IFile)  :: structure_file
  type(IFile)  :: large_supercell_file
  type(IFile)  :: qpoints_file
  type(IFile)  :: no_supercells_file
  type(IFile)  :: supercell_file
  type(IFile)  :: unique_directions_file
  type(String) :: relative_supercell_dir
  type(String) :: supercell_dir
  type(String) :: qpoint_dir
  type(OFile)  :: dynamical_matrix_file
  type(OFile)  :: complex_modes_file
  type(OFile)  :: force_logfile
  type(OFile)  :: qpoint_logfile
  type(OFile)  :: phonon_file
  
  ! Temporary variables.
  integer      :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Read in arguments from user.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  
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
    stop
  endif
  
  degenerate_energy = dble(arguments%value('degenerate_energy'))
  
  ! --------------------------------------------------
  ! Read in previous arguments.
  ! --------------------------------------------------
  setup_harmonic_arguments = Dictionary(setup_harmonic())
  call setup_harmonic_arguments%read_file(wd//'/setup_harmonic.used_settings')
  seedname = setup_harmonic_arguments%value('seedname')
  
  structure_file = IFile(wd//'/structure.dat')
  structure = StructureData(structure_file%lines())
  
  large_supercell_file = IFile(wd//'/large_supercell.dat')
  large_supercell = StructureData(large_supercell_file%lines())
  
  qpoints_file = IFile(wd//'/qpoints.dat')
  qpoints = QpointData(qpoints_file%sections())
  
  no_supercells_file = IFile(wd//'/no_supercells.dat')
  no_supercells = int(no_supercells_file%line(1))
  
  ! --------------------------------------------------
  ! Initialise calculation reader.
  ! --------------------------------------------------
  calculation_reader = CalculationReader(wd)
  
  ! --------------------------------------------------
  ! Calculate the matrix of force constants corresponding to each supercell.
  ! --------------------------------------------------
  allocate( supercells(no_supercells),      &
          & force_constants(no_supercells), &
          & stat=ialloc); call err(ialloc)
  force_logfile = OFile(wd//'/force_constants_log.dat')
  do i=1,no_supercells
    relative_supercell_dir = 'Supercell_'//left_pad(i,str(no_supercells))
    supercell_dir = wd//'/'//relative_supercell_dir
    
    ! Read in supercell structure data.
    supercell_file = IFile(supercell_dir//'/structure.dat')
    supercells(i) = StructureData(supercell_file%lines())
    
    ! Read in symmetry group and unique atoms.
    unique_directions_file = IFile(supercell_dir//'/unique_directions.dat')
    unique_directions = UniqueDirection(unique_directions_file%sections())
    
    ! Calculate force constants.
    force_constants(i) = ForceConstants( supercells(i),            &
                                       & unique_directions,        &
                                       & relative_supercell_dir,   &
                                       & acoustic_sum_rule_forces, &
                                       & calculation_reader,       &
                                       & force_logfile)
    deallocate(unique_directions, stat=ialloc); call err(ialloc)
  enddo
  
  ! --------------------------------------------------
  ! Calculate the dynamical matrix and normal modes at each q-points.
  ! --------------------------------------------------
  allocate( modes_calculated(size(qpoints)),   &
          & dynamical_matrices(size(qpoints)), &
          & stat=ialloc); call err(ialloc)
  qpoint_logfile = OFile(wd//'/dynamical_matrix_log.dat')
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
          dynamical_matrices(i) = DynamicalMatrix( qpoints(i),        &
                                                 & supercells,        &
                                                 & force_constants,   &
                                                 & structure,         &
                                                 & degenerate_energy, &
                                                 & subspace_id,       &
                                                 & qpoint_logfile)
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
    qpoint_dir = wd//'/qpoint_'//left_pad(i,str(size(qpoints)))
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
  phonon_file = OFile(wd//'/'//make_phonon_filename(seedname))
  call write_phonon_file(phonon_file,complex_modes,qpoints,structure)
end subroutine
end module
