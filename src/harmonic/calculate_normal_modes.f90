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
  use harmonic_properties_module
  use unique_directions_module
  use dynamical_matrix_module
  use force_constants_module
  implicit none
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function calculate_normal_modes_keywords() result(keywords)
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                                &
     & KeywordData( 'acoustic_sum_rule',                                      &
     &              'acoustic_sum_rule specifies where the acoustic sum rule  &
     &is applied. The options are "off", "forces", "matrices" and "both".',   &
     &              default_value='both'),                                    &
     & KeywordData( 'degenerate_energy',                                      &
     &              'degenerate_energy is the minimum energy difference &
     &between states before they are considered degenerate. This should be &
     &given in Hartrees.',                                                    &
     &              default_value='1e-5'),                                    &
     & KeywordData( 'calculation_type',                                       &
     &              'calculation_type specifies whether electronic structure &
     &calculations have been run using a user-defined script and &
     &run_harmonic, or should be run through Quip. Settings are: "script" and &
     &"quip".',                                                               &
     &              default_value='script') ]
end function

function calculate_normal_modes_mode() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'calculate_normal_modes'
  output%description = 'Finds harmonic normal modes. Should be called &
     &after run_harmonic.'
  output%keywords = calculate_normal_modes_keywords()
  output%main_subroutine => calculate_normal_modes
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_normal_modes(arguments)
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
  type(String)     :: calculation_type
  
  ! Setup data.
  integer                          :: no_supercells
  type(String)                     :: file_type
  type(String)                     :: seedname
  type(StructureData)              :: structure
  type(StructureData), allocatable :: supercells(:)
  real(dp)                         :: symmetry_precision
  
  ! Force constant data.
  type(UniqueDirection), allocatable :: unique_directions(:)
  type(ForceConstants),  allocatable :: force_constants(:)
  
  ! q-point data.
  type(StructureData)           :: large_supercell
  type(QpointData), allocatable :: qpoints(:)
  
  ! Lte output data.
  logical, allocatable :: modes_calculated(:)
  integer              :: mode
  
  ! Normal modes and their symmetries.
  type(QpointData) :: rotated_qpoint
  
  type(DynamicalMatrix), allocatable :: dynamical_matrices(:)
  type(DynamicalMatrix)              :: rotated_matrix
  
  integer :: mode_id
  integer :: degeneracy_id
  
  complex(dp), allocatable :: pair_overlap(:)
  integer,     allocatable :: degenerate_ids(:)
  complex(dp)              :: phase
  integer                  :: paired_pos
  
  type(ComplexMode), allocatable :: complex_modes(:)
  
  ! Files.
  type(IFile)                    :: no_supercells_file
  type(IFile)                    :: qpoint_file
  type(IFile)                    :: unique_directions_file
  type(OFile)                    :: dynamical_matrix_file
  type(OFile)                    :: complex_modes_file
  type(OFile)                    :: force_logfile
  type(OFile)                    :: qpoint_logfile
  type(StringArray), allocatable :: file_sections(:)
  
  ! Temporary variables.
  integer      :: i,j,k,ialloc
  type(String) :: sdir,qdir
  
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
  
  calculation_type = arguments%value('calculation_type')
  
  ! --------------------------------------------------
  ! Read in previous arguments.
  ! --------------------------------------------------
  setup_harmonic_arguments = Dictionary(setup_harmonic_keywords())
  call setup_harmonic_arguments%read_file(wd//'/setup_harmonic.used_settings')
  file_type = setup_harmonic_arguments%value('file_type')
  seedname = setup_harmonic_arguments%value('seedname')
  symmetry_precision = &
     & dble(setup_harmonic_arguments%value('symmetry_precision'))
  
  no_supercells_file = IFile(wd//'/no_supercells.dat')
  no_supercells = int(no_supercells_file%line(1))
  
  structure = read_structure_file(wd//'/structure.dat',symmetry_precision)
  
  large_supercell = read_structure_file( wd//'/large_supercell.dat', &
                                       & symmetry_precision,         &
                                       & calculate_symmetry=.false.)
  
  qpoint_file = IFile(wd//'/qpoints.dat')
  file_sections = split_into_sections(qpoint_file%lines())
  allocate(qpoints(size(file_sections)), stat=ialloc); call err(ialloc)
  do i=1,size(qpoints)
    qpoints(i) = file_sections(i)
  enddo
  
  ! --------------------------------------------------
  ! Calculate the matrix of force constants corresponding to each supercell.
  ! --------------------------------------------------
  allocate( supercells(no_supercells),      &
          & force_constants(no_supercells), &
          & stat=ialloc); call err(ialloc)
  force_logfile = OFile(wd//'/force_constants_log.dat')
  do i=1,no_supercells
    sdir = wd//'/Supercell_'//left_pad(i,str(no_supercells))
    
    ! Read in supercell structure data.
    supercells(i) = read_structure_file( sdir//'/structure.dat', &
                                       & symmetry_precision)
    
    ! Read in symmetry group and unique atoms.
    unique_directions_file = IFile(sdir//'/unique_directions.dat')
    file_sections = split_into_sections(unique_directions_file%lines())
    allocate( unique_directions(size(file_sections)), &
            & stat=ialloc); call err(ialloc)
    do j=1,size(unique_directions)
      unique_directions(j) = file_sections(j)
    enddo
    
    ! Calculate force constants.
    force_constants(i) = ForceConstants( supercells(i),            &
                                       & unique_directions,        &
                                       & wd,                       &
                                       & sdir,                     &
                                       & file_type,                &
                                       & seedname,                 &
                                       & acoustic_sum_rule_forces, &
                                       & symmetry_precision,       &
                                       & calculation_type,         &
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
  
  degeneracy_id = 1
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
          rotated_qpoint = structure%symmetries(k) * qpoints(j)
          if (rotated_qpoint == qpoints(i)) then
            call qpoint_logfile%print_line('Constructing dynamical matrix and &
               &normal modes at q-point '//i//' using symmetry from those at &
               &q-point '//j)
            dynamical_matrices(i) = rotate_modes( dynamical_matrices(j),   &
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
                                                 & degeneracy_id,     &
                                                 & qpoint_logfile)
          degeneracy_id =                                                  &
             &   maxval(dynamical_matrices(i)%complex_modes%degeneracy_id) &
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
  ! Set mode q-point ids.
  ! --------------------------------------------------
  do i=1,size(dynamical_matrices)
    dynamical_matrices(i)%complex_modes%qpoint_id = qpoints(i)%id
  enddo
  
  ! --------------------------------------------------
  ! Assign ids to all modes, and pair up modes.
  ! --------------------------------------------------
  ! Assign ids
  mode_id = 0
  do i=1,size(dynamical_matrices)
    do j=1,size(dynamical_matrices(i)%complex_modes)
      mode_id = mode_id+1
      dynamical_matrices(i)%complex_modes(j)%id = mode_id
    enddo
  enddo
  
  ! Pair up modes.
  do i=1,size(dynamical_matrices)
    j = first(qpoints%id==qpoints(i)%paired_qpoint_id)
    
    if (j/=i) then
      ! If -q /= q, the modes are set up so that the conjugate of a
      !    mode at q is simply the mode in the same position at -q.
      dynamical_matrices(i)%complex_modes%paired_id = &
         & dynamical_matrices(j)%complex_modes%id
    else
      ! If q = -q, then it is necessary to identify conjugate modes.
      do k=1,size(dynamical_matrices(i)%complex_modes)
        ! Calculate the overlap of the conjugate of mode k with the other
        !    modes. N.B. since the dot product would normaly involve a conjg,
        !    the two conjugates cancel.
        degenerate_ids =                                                  &
           & filter( dynamical_matrices(i)%complex_modes%degeneracy_id == &
           &         dynamical_matrices(i)%complex_modes(k)%degeneracy_id)
        pair_overlap = dynamical_matrices(i)%complex_modes(degenerate_ids) &
                   & * dynamical_matrices(i)%complex_modes(k)
        ! conjg(mode(k)) should equal one mode (down to a phase change),
        !    and have no overlap with any other mode.
        ! Check that this is true, and identify the one mode.
        if (count(abs(pair_overlap)<1e-2)/=size(pair_overlap)-1) then
          call print_line(ERROR//': Modes at qpoint '//i//' do not transform &
             &as expected under parity.')
          call print_line(pair_overlap)
          call err()
        endif
        
        ! Get the position of the paired mode in the list of degenerate modes.
        paired_pos = first(abs(abs(pair_overlap)-1)<1e-2)
        ! Find and normalise the phase of the mode.
        ! Only relevant for when this mode is its own pair, and must be made
        !    real.
        phase = sqrt(pair_overlap(paired_pos))
        phase = phase / abs(phase)
        ! Convert from position in degenerate modes to position in all modes.
        paired_pos = degenerate_ids(paired_pos)
        
        ! Pair up modes.
        if (paired_pos==k) then
          ! The mode is paired to itself. Set paired_id, and divide by the
          !    phase so that the mode is real.
          dynamical_matrices(i)%complex_modes(k)%paired_id = &
             & dynamical_matrices(i)%complex_modes(k)%id
          dynamical_matrices(i)%complex_modes(k)%primitive_displacements =    &
             & cmplxvec(real(                                                 &
             & dynamical_matrices(i)%complex_modes(k)%primitive_displacements &
             & / phase))
        elseif (paired_pos>k) then
          ! The mode is paired to another.
          dynamical_matrices(i)%complex_modes(k)%paired_id = &
             & dynamical_matrices(i)%complex_modes(paired_pos)%id
          dynamical_matrices(i)%complex_modes(paired_pos) = &
             & conjg(dynamical_matrices(i)%complex_modes(k))
        elseif (paired_pos<k) then
          ! The mode has already been handled, when its pair came up.
          cycle
        endif
      enddo
    endif
  enddo
  
  ! --------------------------------------------------
  ! Check all dynamical matrices.
  ! --------------------------------------------------
  ! Run basic checks on each matrix in turn.
  do i=1,size(qpoints)
    call dynamical_matrices(i)%check(structure, qpoint_logfile)
  enddo
  
  ! Check that dynamical matrices at q-points q and q' s.t. qS=q'
  !    correctly rotate onto one another.
  do i=1,size(structure%symmetries)
    do j=1,size(qpoints)
      do k=1,size(qpoints)
        rotated_qpoint = structure%symmetries(i) * qpoints(j)
        if (rotated_qpoint == qpoints(k)) then
          if (qpoints(j)%paired_qpoint_id/=k) then
            cycle
          endif
          rotated_matrix = rotate_modes( dynamical_matrices(j),   &
                                       & structure%symmetries(i), &
                                       & qpoints(j),              &
                                       & qpoints(k))
          call qpoint_logfile%print_line('Comparing symmetrically &
             &equivalent dynamical matrices.')
          call compare_dynamical_matrices( dynamical_matrices(k), &
                                         & rotated_matrix,        &
                                         & qpoint_logfile)
        endif
      enddo
    enddo
  enddo
  
  ! --------------------------------------------------
  ! Write out output.
  ! --------------------------------------------------
  do i=1,size(qpoints)
    qdir = wd//'/qpoint_'//left_pad(i,str(size(qpoints)))
    call mkdir(qdir)
    
    ! Write out dynamical matrix.
    dynamical_matrix_file = OFile(qdir//'/dynamical_matrix.dat')
    call dynamical_matrix_file%print_lines(dynamical_matrices(i))
    
    ! Write out normal modes.
    complex_modes = dynamical_matrices(i)%complex_modes
    complex_modes_file = OFile(qdir//'/complex_modes.dat')
    call complex_modes_file%print_lines(complex_modes,separating_line='')
  enddo
end subroutine
end module
