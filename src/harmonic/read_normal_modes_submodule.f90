submodule (caesar_read_normal_modes_module) caesar_read_normal_modes_submodule
  use caesar_harmonic_module
contains

module procedure startup_read_normal_modes
  type(CaesarMode) :: mode
  
  mode%mode_name = 'read_normal_modes'
  mode%description = 'Reads in normal modes from a separate calculation, e.g. &
     &a Hessian calculation or DFPT dynamical matrix calculation. Acts as a &
     &replacement for setup_harmonic, run_harmonic and &
     &calculate_normal_modes, for situations where harmonic phonons can be &
     &calculated directly rather than requiring a non-diagonal supercell &
     &calculation.'
  mode%keywords = [                                                           &
     & KeywordData( 'file_type',                                              &
     &              'file_type is the file type which will be used for &
     &single-point energy calculations. Settings are: "castep", &
     &"quantum_espresso", "caesar" and "xyz".',                               &
     &              default_value='castep'),                                  &
     & KeywordData( 'seedname',                                               &
     &              'seedname is the seedname from which file names are &
     &constructed.'),                                                         &
     & KeywordData( 'q-point_grid',                                           &
     &              'q-point_grid is the number of q-points in each direction &
     &in a Monkhorst-Pack grid. This should be specified as three integers &
     &separated by spaces.'),                                                 &
     & KeywordData( 'symmetry_precision',                                     &
     &              'symmetry_precision is the tolerance up to which a&
     &symmetry is accepted. In order for a symmetry to be accepted, it must &
     &transform the position of every atom to within symmetry_precision of an &
     &atom of the same element. There must be no other atom within &
     &symmetry_precision of this point. symmetry_precision should be given in &
     &Bohr. Ideally, symmetry_precision should be much smaller than the &
     &minimum inter-atomic distance, and much larger than the geometry &
     &optimisation tolerance.',                                               &
     &              default_value='0.1'),                                     &
     & KeywordData( 'snap_to_symmetry',                                       &
     &              'snap_to_symmetry specifies whether or not to enforce &
     &exact symmetry on the structure before running calculations. If &
     &snap_to_symmetry is set to false, and the symmetry error is large, this &
     &may lead to errors and potential crashes in later calculations.',       &
     &              default_value='true'),                                    &
     & KeywordData( 'loto_direction',                                         &
     &              'loto_direction specifies the direction (in reciprocal &
     &co-ordinates from which the gamma point is approached when calculating &
     &LO/TO corrections. loto_direction should be specified as three &
     &fractions, e.g. 1/2 1/2 1/2. This does not have to be normalised. If &
     &loto_direction is not set, then LO/TO corrections are not calculated. &
     &N.B. LO/TO corrections can currently only be calculated from Castep &
     &DFTB calculations, and only if the q-point grid is 1x1x1. If &
     &loto_direction is specified here, it will be used for the whole &
     &calculation and may not be re-specified in calculate_normal_modes or &
     &calculate_potential.',                                                  &
     &              is_optional = .true.),                                    &
     & KeywordData( 'acoustic_sum_rule',                                      &
     &              'acoustic_sum_rule specifies where the acoustic sum rule &
     &is applied. The options are "off", "forces", "matrices" and "both".',   &
     &              default_value='both')                                     ]
  mode%main_subroutine => read_normal_modes_subroutine
  
  call add_mode(mode)
end procedure

module procedure read_normal_modes_subroutine
  ! Input arguments.
  type(String)         :: file_type
  type(String)         :: seedname
  integer              :: grid(3)
  real(dp)             :: symmetry_precision
  logical              :: snap_to_symmetry
  type(FractionVector) :: loto_direction
  type(String)         :: acoustic_sum_rule
  logical              :: acoustic_sum_rule_forces
  logical              :: acoustic_sum_rule_matrices
  
  ! Dictionary for reproducing settings from setup_harmonic,
  !    for use by later modes.
  type(Dictionary) :: setup_harmonic_arguments
  
  ! Variables normally calculated by setup_harmonic.
  type(HarmonicData)            :: harmonic_data
  type(StructureData)           :: structure
  type(StructureData)           :: large_supercell
  type(QpointData), allocatable :: qpoints(:)
  
  ! Dynamical matrices, normal modes and the Hessian.
  type(DynamicalMatrix), allocatable :: dynamical_matrices(:)
  type(ComplexMode),     allocatable :: complex_modes(:,:)
  type(CartesianHessian)             :: full_hessian
  
  ! Files and directories.
  type(OFile)  :: hessian_logfile
  type(OFile)  :: structure_file
  type(OFile)  :: large_supercell_file
  type(OFile)  :: qpoints_file
  type(String) :: qpoint_dir
  type(OFile)  :: dynamical_matrix_file
  type(OFile)  :: complex_modes_file
  
  type(OFile) :: fc_file
  
  ! Temporary variables.
  integer :: i
  
  ! Read in arguments from user.
  file_type = arguments%value('file_type')
  seedname = arguments%value('seedname')
  grid = int(split_line(arguments%value('q-point_grid')))
  symmetry_precision = dble(arguments%value('symmetry_precision'))
  snap_to_symmetry = lgcl(arguments%value('snap_to_symmetry'))
  if (arguments%is_set('loto_direction')) then
    loto_direction = FractionVector(arguments%value('loto_direction'))
    if (any(grid/=1)) then
      call print_line(ERROR//': Unable to calculate LO/TO splitting if &
         &q-point grid is not 1x1x1.')
      call quit()
    endif
  endif
  acoustic_sum_rule = arguments%value('acoustic_sum_rule')
  if (.not. any(acoustic_sum_rule==tokens('off forces matrices both'))) then
    call print_line(ERROR//': acoustic_sum_rule has been set to an unexpected &
       &value. Please set acoustic_sum_rule to "off", "forces", "matrices" or &
       &"both".')
    call quit()
  endif
  acoustic_sum_rule_forces = any(acoustic_sum_rule==tokens('forces both'))
  acoustic_sum_rule_matrices = any(acoustic_sum_rule==tokens('matrices both'))
  
  ! Calculate harmonic data. This includes:
  !    - The primitive cell structure.
  !    - The large supercell structure.
  !    - The q-points.
  if (arguments%is_set('loto_direction')) then
    harmonic_data = HarmonicData( file_type,          &
                                & seedname,           &
                                & grid,               &
                                & symmetry_precision, &
                                & snap_to_symmetry,   &
                                & loto_direction      )
  else
    harmonic_data = HarmonicData( file_type,          &
                                & seedname,           &
                                & grid,               &
                                & symmetry_precision, &
                                & snap_to_symmetry    )
  endif
  
  structure = harmonic_data%structure
  large_supercell = harmonic_data%large_supercell
  qpoints = harmonic_data%qpoints
  
  ! --------------------------------------------------
  ! Read in Hessian / dynamical matrix / normal modes.
  ! --------------------------------------------------
  ! These are converted to dynamical matrices, and then the large Hessian and
  !    normal modes are re-calculated to ensure symmetry properties are
  !    correct.
  call print_line('Reading in Hessian / dynamical matrices from file.')
  if (file_type=='castep') then
    dynamical_matrices = read_dynamical_matrices_castep( seedname,  &
                                                       & structure, &
                                                       & qpoints    )
  elseif (file_type=='quantum_espresso') then
    dynamical_matrices = read_dynamical_matrices_qe( seedname,        &
                                                   & structure,       &
                                                   & large_supercell, &
                                                   & qpoints          )
  else
    call print_line('File types other than "castep" and "quantum_espresso" &
       &are not currently supported by read_normal_modes.')
    call quit()
  endif
  
  ! Calculate normal modes.
  call print_line('Calculating normal modes')
  complex_modes = calculate_normal_modes( structure,         &
                                        & qpoints,           &
                                        & dynamical_matrices )
  
  ! Calculate the Hessian for the full supercell.
  call print_line('Calculating full Hessian.')
  hessian_logfile = OFile('hessian_log.dat')
  full_hessian = reconstruct_hessian( large_supercell,    &
                                    & qpoints,            &
                                    & dynamical_matrices, &
                                    & hessian_logfile     )
  
  ! --------------------------------------------------
  ! Write everything to file.
  ! --------------------------------------------------
  setup_harmonic_arguments = Dictionary(CaesarMode('setup_harmonic'))
  call setup_harmonic_arguments%set(arguments)
  call setup_harmonic_arguments%write_file('setup_harmonic.used_settings')
  
  structure_file = OFile('structure.dat')
  call structure_file%print_lines(structure)
  
  large_supercell_file = OFile('large_supercell.dat')
  call large_supercell_file%print_lines(large_supercell)
  
  qpoints_file = OFile('qpoints.dat')
  call qpoints_file%print_lines(qpoints,separating_line='')
  
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
  
  fc_file = OFile('tmp.fc')
  call write_qe_force_constants_file( fc_file,        &
                                    & full_hessian,   &
                                    & structure,      &
                                    & large_supercell )
end procedure

module procedure read_dynamical_matrices_castep
  type(ComplexMode), allocatable :: complex_modes(:,:)
  
  type(IFile) :: phonon_file
  
  integer :: i
  
  phonon_file = IFile(make_castep_phonon_filename(seedname))
  
  complex_modes = read_castep_phonon_file(phonon_file,structure,qpoints)
  
  output = [(DynamicalMatrix(complex_modes(:,i)), i=1, size(qpoints))]
end procedure

module procedure read_dynamical_matrices_qe
  type(CartesianHessian) :: hessian
  
  integer :: i
  
  hessian = read_qe_force_constants_file( directory = str('.'),       &
                                        & seedname  = seedname,       &
                                        & supercell = large_supercell )
  
  output = [( DynamicalMatrix( dblevec(qpoints(i)%qpoint),    &
            &                  large_supercell,               &
            &                  hessian                     ), &
            & i=1,                                            &
            & size(qpoints)                                   )]
end procedure

module procedure calculate_normal_modes
  logical, allocatable :: modes_calculated(:)
  
  integer :: subspace_id
  
  integer :: i,j,k,ialloc
  
  allocate( output(structure%no_modes,size(qpoints)), &
          & stat=ialloc); call err(ialloc)
  modes_calculated = [(.false.,i=1,size(qpoints))]
  subspace_id = 1
  iter : do while(.not. all(modes_calculated))
    ! First, check if there is an un-calculated q-point whose conjugate has
    !    already been calculated.
    do i=1,size(qpoints)
      if (modes_calculated(i)) then
        cycle
      endif
      
      j = first(qpoints%id==qpoints(i)%paired_qpoint_id)
      
      if (modes_calculated(j)) then
        output(:,i) = conjg(output(:,j))
        modes_calculated(i) = .true.
        cycle iter
      endif
    enddo
    
    ! Second, check if there is a q-point which is related by symmetry to an
    !    already calculated q-point.
    do i=1,size(qpoints)
      if (modes_calculated(i)) then
        cycle
      endif
      
      do j=1,size(qpoints)
        if (.not. modes_calculated(j)) then
          cycle
        endif
        
        do k=1,size(structure%symmetries)
          if (structure%symmetries(k)*qpoints(j)==qpoints(i)) then
            output(:,i) = transform( output(:,j),      &
                                   & structure%symmetries(k), &
                                   & qpoints(j),              &
                                   & qpoints(i)               )
            modes_calculated(i) = .true.
            cycle iter
          endif
        enddo
      enddo
    enddo
    
    ! Finally, construct the modes from a dynamical matrix.
    i = first(.not. modes_calculated)
    output(:,i) = ComplexMode( dynamical_matrices(i),                   &
                             & structure,                               &
                             & modes_real=qpoints(i)%is_paired_qpoint() )
    output(:,i) = process_modes(output(:,i),structure,qpoints(i),subspace_id)
    subspace_id = maxval(output(:,i)%subspace_id) + 1
    modes_calculated(i) = .true.
  enddo iter
end procedure
end submodule
