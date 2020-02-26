! ======================================================================
! Maps the VSCF potential along each normal mode.
! ======================================================================
module map_vscf_modes_module
  use common_module
  
  use states_module
  use anharmonic_common_module
  use potentials_module
  
  use mode_map_module
  implicit none
  
  private
  
  public :: startup_map_vscf_modes
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_map_vscf_modes()
  implicit none
  
  type(CaesarMode) :: mode
  
  mode%mode_name = 'map_vscf_modes'
  mode%description = 'Maps the VSCF potential along normal modes. &
     &Should be run after calculate_vscf_potential.'
  mode%keywords = [                                                           &
     & KeywordData( 'no_single_mode_samples',                                 &
     &              'no_single_mode_samples is the number of points (either &
     &side of zero) along each mode at which the VSCF potential will be &
     &sampled when determining the effective frequency with which the &
     &harmonic basis along that mode will be constructed.')                   ]
  mode%main_subroutine => map_vscf_modes_subroutine
  
  call add_mode(mode)
end subroutine

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine map_vscf_modes_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Input arguments.
  integer :: no_single_mode_samples
  
  ! Anharmonic data.
  type(AnharmonicData)                  :: anharmonic_data
  type(QpointData),         allocatable :: qpoints(:)
  type(ComplexMode),        allocatable :: complex_modes(:)
  type(RealMode),           allocatable :: real_modes(:)
  type(DegenerateSubspace), allocatable :: subspaces(:)
  real(dp)                              :: frequency_of_max_displacement
  
  ! Single-subspace potentials.
  type(PotentialPointer), allocatable :: subspace_potentials(:)
  
  ! Variables for calculating displacements.
  real(dp),      allocatable  :: unscaled_mode_displacements(:)
  real(dp),      allocatable  :: mode_displacements(:)
  real(dp),      allocatable  :: l2_cartesian_displacements(:)
  type(RealMode)              :: mode
  type(ModeMap), allocatable  :: mode_maps(:)
  type(RealModeDisplacement)  :: real_mode_displacement
  type(CartesianDisplacement) :: cartesian_displacement
  
  ! Files and directories.
  type(IFile)  :: anharmonic_data_file
  type(IFile)  :: subspace_potentials_file
  type(String) :: subspace_dir
  type(OFile)  :: mode_maps_file
  
  integer :: i,j,k,ialloc
  
  ! Read in arguments.
  no_single_mode_samples = int(arguments%value('no_single_mode_samples'))
  
  ! Read in anharmonic data.
  anharmonic_data_file = IFile('anharmonic_data.dat')
  anharmonic_data = AnharmonicData(anharmonic_data_file%lines())
  
  qpoints = anharmonic_data%qpoints
  complex_modes = anharmonic_data%complex_modes
  real_modes = anharmonic_data%real_modes
  subspaces = anharmonic_data%degenerate_subspaces
  frequency_of_max_displacement = anharmonic_data%frequency_of_max_displacement
  
  ! Read in single-subspace potentials.
  subspace_potentials_file = IFile('subspace_potentials.dat')
  subspace_potentials = PotentialPointer(                                &
     & subspace_potentials_file%sections(separating_line=repeat('=',70)) )
  
  ! --------------------------------------------------
  ! Calculate the value of the VSCF potential along each mode.
  ! --------------------------------------------------
  unscaled_mode_displacements =                               &
     & [(i,i=-no_single_mode_samples,no_single_mode_samples)] &
     & * anharmonic_data%maximum_weighted_displacement        &
     & / no_single_mode_samples
  allocate( l2_cartesian_displacements(size(unscaled_mode_displacements)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(subspaces)
    ! Scale displacement by 1/sqrt(frequency).
    mode_displacements = unscaled_mode_displacements                &
                     & * sqrt( frequency_of_max_displacement        &
                     &       / max( subspaces(i)%frequency,         &
                     &              frequency_of_max_displacement ) )
    
    
    allocate(mode_maps(size(subspaces(i))), stat=ialloc); call err(ialloc)
    do j=1,size(subspaces(i))
      mode = real_modes(first(real_modes%id==subspaces(i)%mode_ids(j)))
      do k=1,size(l2_cartesian_displacements)
        real_mode_displacement = RealModeDisplacement( &
                             & [real_modes(i)],        &
                             & [mode_displacements(j)] )
        cartesian_displacement = CartesianDisplacement( &
                & real_mode_displacement,               &
                & anharmonic_data%anharmonic_supercell, &
                & real_modes,                           &
                & qpoints                               )
        l2_cartesian_displacements(j) =                  &
           & sqrt(sum( cartesian_displacement%vectors    &
           &         * cartesian_displacement%vectors )) &
           & / anharmonic_data%anharmonic_supercell%sc_size
      enddo
      l2_cartesian_displacements(:no_single_mode_samples) = &
         & -l2_cartesian_displacements(:no_single_mode_samples)

      mode_maps(j) = ModeMap( mode_displacements,               &
                            & l2_cartesian_displacements,       &
                            & mode,                             &
                            & subspace_potentials(i),           &
                            & anharmonic_data = anharmonic_data )
    enddo
    
    ! Write out displacements.
    subspace_dir = 'subspace_'// &
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
