! ======================================================================
! Maps the VSCF potential along each normal mode.
! ======================================================================
module map_vscf_modes_module
  use common_module
  
  use states_module
  use anharmonic_common_module
  use potentials_module
  
  use setup_anharmonic_module
  use mode_map_module
  implicit none
  
  private
  
  public :: map_vscf_modes
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function map_vscf_modes() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'map_vscf_modes'
  output%description = 'Maps the VSCF potential along normal modes. &
     &Should be run after calculate_states.'
  output%keywords = [                                                         &
     & KeywordData( 'no_single_mode_samples',                                 &
     &              'no_single_mode_samples is the number of points (either &
     &side of zero) along each mode at which the VSCF potential will be &
     &sampled when determining the effective frequency with which the &
     &harmonic basis along that mode will be constructed.')                   ]
  output%main_subroutine => map_vscf_modes_subroutine
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine map_vscf_modes_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! Input arguments.
  integer :: no_single_mode_samples
  
  ! Anharmonic data.
  type(AnharmonicData)                  :: anharmonic_data
  type(ComplexMode),        allocatable :: complex_modes(:)
  type(RealMode),           allocatable :: real_modes(:)
  type(DegenerateSubspace), allocatable :: subspaces(:)
  real(dp)                              :: frequency_of_max_displacement
  
  ! Single-subspace potentials, bases and ground states.
  type(PotentialPointer), allocatable :: subspace_potentials(:)
  type(SubspaceBasis),    allocatable :: basis(:)
  type(VscfGroundState),  allocatable :: ground_states(:)
  
  ! Variables for calculating displacements.
  real(dp),      allocatable :: displacements(:)
  real(dp),      allocatable :: scaled_displacements(:)
  type(RealMode)             :: mode
  type(ModeMap), allocatable :: mode_maps(:)
  
  ! Files and directories.
  type(IFile)  :: anharmonic_data_file
  type(IFile)  :: subspace_potentials_file
  type(IFile)  :: basis_file
  type(IFile)  :: ground_states_file
  type(String) :: subspace_dir
  type(OFile)  :: mode_maps_file
  
  integer :: i,j,ialloc
  
  ! Read in arguments.
  wd = arguments%value('working_directory')
  no_single_mode_samples = int(arguments%value('no_single_mode_samples'))
  
  ! Read in anharmonic data.
  anharmonic_data_file = IFile(wd//'/anharmonic_data.dat')
  anharmonic_data = AnharmonicData(anharmonic_data_file%lines())
  
  complex_modes = anharmonic_data%complex_modes
  real_modes = anharmonic_data%real_modes
  subspaces = anharmonic_data%degenerate_subspaces
  frequency_of_max_displacement = anharmonic_data%frequency_of_max_displacement
  
  ! Read in single-subspace potentials, bases and ground states.
  subspace_potentials_file = IFile(wd//'/subspace_potentials.dat')
  subspace_potentials = PotentialPointer(                                &
     & subspace_potentials_file%sections(separating_line=repeat('=',70)) )
  
  basis_file = IFile(wd//'/basis.dat')
  basis = SubspaceBasis(basis_file%sections())
  
  ground_states_file = IFile(wd//'/ground_state.dat')
  ground_states = VscfGroundState(ground_states_file%sections())
  
  ! --------------------------------------------------
  ! Calculate the value of the VSCF potential along each mode.
  ! --------------------------------------------------
  displacements = [(i,i=-no_single_mode_samples,no_single_mode_samples)] &
              & * anharmonic_data%maximum_weighted_displacement          &
              & / no_single_mode_samples
  do i=1,size(subspaces)
    ! Scale displacement by 1/sqrt(frequency).
    scaled_displacements = displacements                              &
                       & * sqrt( frequency_of_max_displacement        &
                       &       / max( basis(i)%frequency,             &
                       &              frequency_of_max_displacement ) )
    
    
    allocate(mode_maps(size(subspaces(i))), stat=ialloc); call err(ialloc)
    do j=1,size(subspaces(i))
      mode = real_modes(first(real_modes%id==subspaces(i)%mode_ids(j)))
      mode_maps(j) = ModeMap( scaled_displacements,  &
                            & mode,                  &
                            & subspace_potentials(i) )
    enddo
    
    ! Write out displacements.
    subspace_dir = wd//'/subspace_'// &
       & left_pad(subspaces(i)%id,str(maxval(subspaces%id,1)))
    call mkdir(subspace_dir)
    mode_maps_file = OFile(subspace_dir//'/vscf_mode_maps.dat')
    call mode_maps_file%print_line(                    &
       & 'Harmonic frequencies: '//subspaces%frequency )
    call mode_maps_file%print_line('')
    call mode_maps_file%print_lines(mode_maps, separating_line='')
    
    deallocate(mode_maps, stat=ialloc); call err(ialloc)
  enddo
end subroutine
end module
