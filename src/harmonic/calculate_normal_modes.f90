! ======================================================================
! The third stage of Caesar.
! Uses the forces calculated previously to calculate harmonic normal modes.
! ======================================================================
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
  
  keywords = [KeywordData::]
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
  use lift_degeneracies_module
  use atom_module
  use ofile_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! Input arguments.
  type(Dictionary) :: setup_harmonic_arguments
  
  ! No. supercells file.
  type(IFile) :: no_supercells_file
  
  ! Setup data.
  integer             :: no_supercells
  type(String)        :: file_type
  type(String)        :: seedname
  type(StructureData) :: structure
  type(StructureData) :: supercell
  real(dp)            :: harmonic_displacement
  
  ! Force constant data.
  type(UniqueDirection), allocatable :: unique_directions(:)
  type(ForceConstants)               :: force_constants
  
  ! q-point data.
  type(StructureData)           :: large_supercell
  type(QpointData), allocatable :: qpoints(:)
  type(QpointData)              :: qpoint
  
  type(IntVector) :: gvector
  
  ! Lte output data.
  logical,                      allocatable :: modes_calculated(:)
  integer                                   :: mode
  type(String)                              :: mode_string
  
  ! Normal modes and their symmetries.
  logical,             allocatable :: translational(:)
  
  type(AtomData) :: atom
  type(ComplexVector) :: prim_disp
  
  type(DynamicalMatrix), allocatable :: dynamical_matrices(:)
  
  logical :: at_gamma
  logical :: paired_gvec
  
  ! Logfile.
  type(OFile) :: logfile
  
  ! Temporary variables.
  integer      :: i,j,k,ialloc
  type(String) :: sdir,qdir
  
  ! --------------------------------------------------
  ! Read in arguments from user.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  
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
  
  ! Make output directories.
  do i=1,size(qpoints)
    qdir = wd//'/qpoint_'//left_pad(i,str(size(qpoints)))
    call mkdir(qdir)
  enddo
  
  ! --------------------------------------------------
  ! Loop over supercells, calculating dynamical matrices and normal modes.
  ! Transfer information to the relevant q-points.
  ! --------------------------------------------------
  allocate( modes_calculated(size(qpoints)),   &
          & dynamical_matrices(size(qpoints)), &
          & stat=ialloc); call err(ialloc)
  modes_calculated = .false.
  do i=1,no_supercells
    sdir = wd//'/Supercell_'//left_pad(i,str(no_supercells))
    
    ! Read in supercell structure data.
    supercell = read_structure_file(sdir//'/structure.dat')
    
    ! Read in symmetry group and unique atoms.
    unique_directions = read_unique_directions_file( &
       & sdir//'/unique_directions.dat')
    
    ! Calculate force constants.
    logfile = sdir//'/force_constants_log.dat'
    force_constants = ForceConstants( supercell,         &
                                    & unique_directions, &
                                    & sdir,              &
                                    & file_type,         &
                                    & seedname,          &
                                    & logfile)
    
    ! Pick out Supercell G-vectors corresponding to required q-points.
    do_j : do j=1,size(qpoints)
      qdir = wd//'/qpoint_'//left_pad(j,str(size(qpoints)))
      qpoint = qpoints(j)
      
      if (.not. qpoint%to_simulate) then
        cycle
      elseif (qpoint%supercell_matrix/=supercell%supercell) then
        cycle
      endif
      
      do k=1,supercell%sc_size
        gvector = supercell%gvectors(k)
        
        if (qpoint%supercell_gvector/=gvector) then
          cycle
        endif
        
        ! Construct the dynamical matrix at the q-point.
        at_gamma = all(int(gvector)==0)
        paired_gvec = supercell%paired_gvec(k)==k
        logfile = qdir//'/dynamical_matrix_log.dat'
        call logfile%print_line('Constructing dynamical matrix and normal &
           &modes through direct calculation using supercell '//i)
        dynamical_matrices(j) = DynamicalMatrix( qpoint,          &
                                               & at_gamma,        &
                                               & paired_gvec,     &
                                               & supercell,       &
                                               & force_constants, &
                                               & logfile)
        
        ! Lift degeneracies.
        !dynamical_matices(j)%complex_modes = lift_degeneracies( &
        !                 & dynamical_matrices(j)%complex_modes, &
        !                 & structure)
        
        modes_calculated(j) = .true.
        cycle do_j
      enddo
          
      call print_line(CODE_ERROR//': No matching G-vector found for q-point &
         &'//j//'.')
      call err()
    enddo do_j
  enddo
  
  ! --------------------------------------------------
  ! Calculate information at remaining q-points from symmetry.
  ! --------------------------------------------------
  do_i : do i=1,size(qpoints)
    qdir = wd//'/qpoint_'//left_pad(i,str(size(qpoints)))
    if (.not. modes_calculated(i)) then
      logfile = qdir//'/dynamical_matrix_log.dat'
      do j=1,i-1
        ! Check if q-point j + q-point i is a G-vector.
        if (qpoints(i)%paired_qpoint==j) then
          call logfile%print_line('Constructing dynamical matrix and normal &
             &modes as the conjugate of those at q-point '//j)
          dynamical_matrices(i) = conjg(dynamical_matrices(j))
          modes_calculated(i) = .true.
          cycle do_i
        endif
        
        ! Check if q-point j is symmetrically related to q-point i.
        do k=1,size(structure%symmetries)
          if (structure%symmetries(k)%rotation * qpoints(i)%scaled_qpoint &
                                            & == qpoints(j)%scaled_qpoint) then
            call logfile%print_line('Constructing dynamical matrix and normal &
               &modes using symmetry from those at q-point '//j)
            dynamical_matrices(i) = rotate_modes( dynamical_matrices(j),   &
                                                & structure%symmetries(k), &
                                                & structure,               &
                                                & qpoints(j),              &
                                                & logfile)
            modes_calculated(i) = .true.
            cycle do_i
          endif
        enddo
      enddo
      
      call print_line(CODE_ERROR//': No rotationally equivalent q-point &
         &found for q-point '//i//'.')
      call err()
    endif
  enddo do_i
  
  ! --------------------------------------------------
  ! Write out output.
  ! --------------------------------------------------
  do i=1,size(qpoints)
    qdir = wd//'/qpoint_'//left_pad(i,str(size(qpoints)))
    
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

! ----------------------------------------------------------------------
! Calculates the set of phonon normal modes at each supercell G-vector.
! Returns the real part of the non-mass-reduced polarisation vector, which
!    is the pattern of displacement corresponding to the normal mode.
! ----------------------------------------------------------------------
!function evaluate_normal_modes(supercell,force_constants) &
!   & result(output)
!  use structure_module
!  use min_images_module
!  use atom_module
!  use force_constants_module
!  implicit none
!  
!  type(StructureData),  intent(in)   :: supercell
!  type(ForceConstants), intent(in)   :: force_constants
!  type(DynamicalMatrix), allocatable :: output(:)
!  
!  ! Working variables.
!  type(MinImages),     allocatable :: min_images(:,:,:)
!  type(RealVector)                 :: qpoint
!  type(ComplexVector), allocatable :: pol_vec(:,:)
!  logical,             allocatable :: translational(:)
!  
!  type(AtomData) :: atom
!  type(ComplexVector) :: prim_disp
!  
!  ! Indices.
!  integer :: mode
!  integer :: gvector,gvector_p
!  
!  logical :: at_gamma
!  logical :: paired_gvec
!  
!  ! Temporary variables.
!  integer :: i,j,ialloc
!  
!  ! Allocate output and translational modes label.
!  allocate(output(supercell%sc_size), stat=ialloc); call err(ialloc)
!  do gvector=1,supercell%sc_size
!    allocate( output(gvector)%complex_modes(supercell%no_modes_prim), &
!            & stat=ialloc); call err(ialloc)
!  enddo
!  
!  ! Calculate minimum R-vector separations between atoms.
!  
!  ! Calculate dynamical matrices, and their eigenstuff.
!  do gvector=1,supercell%sc_size
!    gvector_p = supercell%paired_gvec(gvector)
!    
!    if (gvector_p<gvector) then
!      cycle
!    endif
!    
!    ! Construct the dynamical matrix at G.
!    qpoint = dble( supercell%recip_supercell    &
!         &       * supercell%gvectors(gvector)) &
!         & / supercell%sc_size
!    
!    at_gamma = all(int(supercell%gvectors(gvector))==0)
!    paired_gvec = gvector==gvector_p
!    output(gvector) = DynamicalMatrix( qpoint,          &
!                                     & at_gamma,        &
!                                     & paired_gvec,     &
!                                     & supercell,       &
!                                     & force_constants, &
!                                     & min_images)
!    
!    if (.not. paired_gvec) then
!      output(gvector_p) = conjg(output(gvector))
!    endif
!  enddo
!end function
end module
