! ======================================================================
! The third stage of Caesar.
! Uses the forces calculated previously to calculate harmonic normal modes.
! ======================================================================
! Calculates the set of phonon normal modes at each supercell G-vector.
! Calculates the real part of the non-mass-reduced polarisation vector, which
!    is the pattern of displacement corresponding to the normal mode.
module calculate_normal_modes_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function calculate_normal_modes_keywords() result(keywords)
  use keyword_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                               &
     & KeywordData( 'acoustic_sum_rule',                                     &
     &              'acoustic_sum_rule specifies where the acoustic sum rule &
     &is applied. The options are "off", "forces", "matrices" and "both".',  &
     &              default_value='both'),                                   &
     & KeywordData( 'degenerate_energy',                                     &
     &              'degenerate_energy is the minimum energy difference &
     &between states before they are considered degenerate. This should be &
     &given in Bohr.',                                                       &
     &              default_value='1e-5') ]
end function

function calculate_normal_modes_mode() result(output)
  use caesar_modes_module
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
  use utils_module, only : mkdir
  use ifile_module
  use linear_algebra_module
  use setup_harmonic_module
  use structure_module
  use output_file_module
  use lte_module
  use unique_directions_module
  use group_module
  use qpoints_module
  use dictionary_module
  use normal_mode_module
  use dynamical_matrix_module
  use force_constants_module
  use atom_module
  use ofile_module
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
  
  ! No. supercells file.
  type(IFile) :: no_supercells_file
  
  ! Setup data.
  integer                          :: no_supercells
  type(String)                     :: file_type
  type(String)                     :: seedname
  type(StructureData)              :: structure
  type(StructureData), allocatable :: supercells(:)
  real(dp)                         :: harmonic_displacement
  
  ! Force constant data.
  type(UniqueDirection), allocatable :: unique_directions(:)
  type(ForceConstants),  allocatable :: force_constants(:)
  
  ! q-point data.
  type(StructureData)           :: large_supercell
  type(QpointData), allocatable :: qpoints(:)
  
  ! Lte output data.
  logical, allocatable :: modes_calculated(:)
  integer              :: mode
  type(String)         :: mode_string
  
  ! Normal modes and their symmetries.
  type(FractionVector) :: rotated_qpoint
  
  type(DynamicalMatrix), allocatable :: dynamical_matrices(:)
  type(DynamicalMatrix)              :: rotated_matrix
  
  ! Logfiles.
  type(OFile) :: force_logfile
  type(OFile) :: qpoint_logfile
  
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
  
  ! --------------------------------------------------
  ! Read in previous arguments.
  ! --------------------------------------------------
  setup_harmonic_arguments = setup_harmonic_keywords()
  call setup_harmonic_arguments%read_file(wd//'/setup_harmonic.used_settings')
  file_type = setup_harmonic_arguments%value('file_type')
  seedname = setup_harmonic_arguments%value('seedname')
  harmonic_displacement = dble(setup_harmonic_arguments%value( &
     & 'harmonic_displacement'))
  
  no_supercells_file = wd//'/no_supercells.dat'
  no_supercells = int(no_supercells_file%line(1))
  
  structure = read_structure_file(wd//'/structure.dat')
  
  large_supercell = read_structure_file(wd//'/large_supercell.dat')
  
  qpoints = read_qpoints_file(wd//'/qpoints.dat')
  
  ! --------------------------------------------------
  ! Calculate the matrix of force constants corresponding to each supercell.
  ! --------------------------------------------------
  allocate( supercells(no_supercells),      &
          & force_constants(no_supercells), &
          & stat=ialloc); call err(ialloc)
  force_logfile = wd//'/force_constants_log.dat'
  do i=1,no_supercells
    sdir = wd//'/Supercell_'//left_pad(i,str(no_supercells))
    
    ! Read in supercell structure data.
    supercells(i) = read_structure_file(sdir//'/structure.dat')
    
    ! Read in symmetry group and unique atoms.
    unique_directions = read_unique_directions_file( &
       & sdir//'/unique_directions.dat')
    
    ! Calculate force constants.
    force_constants(i) = ForceConstants( supercells(i),            &
                                       & unique_directions,        &
                                       & sdir,                     &
                                       & file_type,                &
                                       & seedname,                 &
                                       & acoustic_sum_rule_forces, &
                                       & force_logfile)
  enddo
  
  ! --------------------------------------------------
  ! Calculate the dynamical matrix and normal modes at each q-points.
  ! --------------------------------------------------
  allocate( modes_calculated(size(qpoints)),   &
          & dynamical_matrices(size(qpoints)), &
          & stat=ialloc); call err(ialloc)
  qpoint_logfile = wd//'/dynamical_matrix_log.dat'
  modes_calculated = .false.
  
  iter : do while (.not. all(modes_calculated))
    ! ------------------------------
    ! First, check if there is an un-calculated q-point whose conjugate has
    !    already been calculated.
    ! ------------------------------
    do i=1,size(qpoints)
      if (modes_calculated(i)) then
        cycle
      endif
      
      if (modes_calculated(qpoints(i)%paired_qpoint)) then
        call qpoint_logfile%print_line('Constructing dynamical matrix and &
           &normal modes at q-point '//i//' as the conjugate of those at &
           &q-point '//j)
        dynamical_matrices(i) = conjg( &
           & dynamical_matrices(qpoints(i)%paired_qpoint))
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
          rotated_qpoint = structure%symmetries(k)%recip_rotation &
                       & * qpoints(j)%qpoint
          if (rotated_qpoint == qpoints(i)%qpoint) then
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
                                                 & qpoint_logfile)
          
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
  !    correctly rotate onto one another.
  do i=1,size(structure%symmetries)
    do j=1,size(qpoints)
      do k=1,size(qpoints)
        rotated_qpoint = structure%symmetries(i)%recip_rotation &
                     & * qpoints(j)%qpoint
        if (rotated_qpoint == qpoints(k)%qpoint) then
          if (qpoints(j)%paired_qpoint/=k) then
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
    call write_dynamical_matrix_file( dynamical_matrices(i), &
                                    & qdir//'/dynamical_matrix.dat')
    
    ! Write out normal modes.
    do mode=1,structure%no_modes
      mode_string = left_pad(mode,str(structure%no_modes))
      call dynamical_matrices(i)%complex_modes(mode)%write_file( &
         & qdir//'/mode_'//mode_string//'.dat')
    enddo
  enddo
end subroutine
end module
