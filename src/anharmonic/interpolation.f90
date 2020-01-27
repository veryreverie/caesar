! ======================================================================
! Interpolate a potential or stress to a fine grid of q-points,
!    and calculate the thermodynamic variables.
! ======================================================================
! N.B. THIS MODULE IS NOT SIZE_CONSISTENT. FOR NOW IT IS COMMENTED OUT.
module interpolation_module
!  use common_module
!  
!  use states_module
!  
!  use anharmonic_common_module
!  implicit none
!  
!  private
!  
!  public :: InterpolatedThermodynamics
!  public :: calculate_interpolated_stress
!  public :: calculate_interpolated_dispersion
!  public :: calculate_interpolated_thermodynamics
!  
!  type, extends(NoDefaultConstructor) :: InterpolatedThermodynamics
!    type(PhononDos)         :: vscha_dos
!    type(ThermodynamicData) :: vscha_thermodynamics
!  end type
!  
!  interface InterpolatedThermodynamics
!    module procedure new_InterpolatedThermodynamics
!  end interface
!contains
!
!function new_InterpolatedThermodynamics(vscha_dos,vscha_thermodynamics) &
!   & result(this)
!  implicit none
!  
!  type(PhononDos),         intent(in) :: vscha_dos
!  type(ThermodynamicData), intent(in) :: vscha_thermodynamics
!  type(InterpolatedThermodynamics)    :: this
!  
!  this%vscha_dos = vscha_dos
!  this%vscha_thermodynamics = vscha_thermodynamics
!end function
!
!function calculate_interpolated_stress(stress,degenerate_frequency, &
!   & fine_qpoints,thermal_energy,min_frequency,harmonic_supercell,  &
!   & harmonic_hessian,harmonic_min_images,subspaces,subspace_bases, &
!   & basis_states,anharmonic_min_images,anharmonic_data) result(output)
!  implicit none
!  
!  class(StressData),        intent(in) :: stress
!  real(dp),                 intent(in) :: degenerate_frequency
!  type(RealVector),         intent(in) :: fine_qpoints(:)
!  real(dp),                 intent(in) :: thermal_energy
!  real(dp),                 intent(in) :: min_frequency
!  type(StructureData),      intent(in) :: harmonic_supercell
!  type(CartesianHessian),   intent(in) :: harmonic_hessian
!  type(MinImages),          intent(in) :: harmonic_min_images(:,:)
!  type(DegenerateSubspace), intent(in) :: subspaces(:)
!  class(SubspaceBasis),     intent(in) :: subspace_bases(:)
!  class(BasisStates),       intent(in) :: basis_states(:)
!  type(MinImages),          intent(in) :: anharmonic_min_images(:,:)
!  type(AnharmonicData),     intent(in) :: anharmonic_data
!  type(RealMatrix)                     :: output
!  
!  ! Variables at a given q-point.
!  type(DynamicalMatrix)                 :: dynamical_matrix
!  type(ComplexMode),        allocatable :: qpoint_modes(:)
!  type(ComplexMode),        allocatable :: modes(:)
!  type(DegenerateSubspace), allocatable :: fine_subspaces(:)
!  
!  ! Variables in a given subspace.
!  type(ComplexMode), allocatable :: subspace_modes(:)
!  type(StressPrefactors)         :: stress_prefactors
!  type(StressPointer)            :: fine_stress
!  type(RealMatrix)               :: potential_stress
!  type(ThermodynamicData)        :: thermodynamic_data
!  
!  ! Temporary variables
!  integer :: i,j,k
!  
!  output = dblemat(zeroes(3,3))
!  do i=1,size(fine_qpoints)
!    ! Construct normal modes from harmonic potential.
!    ! Include modes from both q and -q.
!    dynamical_matrix = DynamicalMatrix( fine_qpoints(i),    &
!                                      & harmonic_supercell, &
!                                      & harmonic_hessian,   &
!                                      & harmonic_min_images )
!    
!    qpoint_modes = ComplexMode(dynamical_matrix, harmonic_supercell)
!    qpoint_modes%id = [(j,j=1,size(qpoint_modes))]
!    qpoint_modes%paired_id = [(j+size(qpoint_modes),j=1,size(qpoint_modes))]
!    
!    modes = [( [qpoint_modes(j),conjg(qpoint_modes(j))], &
!             & j=1,                                      &
!             & size(qpoint_modes)                        )]
!    
!    ! Split modes into subspaces.
!    fine_subspaces = generate_fine_subspaces( qpoint_modes,        &
!                                            & degenerate_frequency )
!    
!    do j=1,size(fine_subspaces)
!      if (fine_subspaces(j)%frequency>min_frequency) then
!        ! Generate the modes in the subspace.
!        subspace_modes = fine_subspaces(j)%modes(modes)
!        
!        ! Calculate stress prefactors.
!        stress_prefactors = StressPrefactors(fine_subspaces(j), subspace_modes)
!        
!        ! Interpolate stress.
!        fine_stress = stress%interpolate( fine_qpoints(i),       &
!                                        & fine_subspaces(j),     &
!                                        & subspace_modes,        &
!                                        & anharmonic_min_images, &
!                                        & anharmonic_data        )
!        
!        ! Take the harmonic expectation of the stress.
!        potential_stress = fine_stress%harmonic_expectation(   &
!                       &        fine_subspaces(j)%frequency,   &
!                       &        thermal_energy,                &
!                       &        anharmonic_data              ) &
!                       & / (2*size(subspace_modes))
!        
!        ! Calculate total stress, including kinetic stress.
!        thermodynamic_data =                                             &
!           & ThermodynamicData( thermal_energy,                          &
!           &                    fine_subspaces(j)%frequency,             &
!           &                    stress_prefactors%average_prefactor(),   &
!           &                    potential_stress,                        &
!           &                    anharmonic_data%structure%volume       ) &
!           & * size(subspace_modes)
!        
!        output = output + thermodynamic_data%stress
!      endif
!    enddo
!    
!    if (modulo(i,size(fine_qpoints)/10)==0) then
!      call print_line('Stress: '//i//' of '//size(fine_qpoints)// &
!         & ' q-points sampled.')
!    endif
!  enddo
!  
!  ! Normalise to be per primitive cell.
!  ! N.B. at each q-point, q and -q are both considered, hence the factor of 2.
!  output = output / (2*size(fine_qpoints))
!  
!  ! Add in the static-lattice stress.
!  ! TODO: Correct VSCF form for when stress is not harmonic.
!  !output = output + stress%undisplaced_stress() &
!  !              & / anharmonic_data%anharmonic_supercell%sc_size
!end function
!
!function calculate_interpolated_dispersion(potential,degenerate_frequency,  &
!   & path_string,no_path_points,thermal_energy,min_frequency,               &
!   & harmonic_supercell,harmonic_hessian,harmonic_min_images,subspaces,     &
!   & subspace_bases,subspace_states,anharmonic_min_images,convergence_data, &
!   & anharmonic_data) result(output)
!  implicit none
!  
!  class(PotentialData),     intent(in)    :: potential
!  real(dp),                 intent(in)    :: degenerate_frequency
!  type(String),             intent(in)    :: path_string
!  integer,                  intent(in)    :: no_path_points
!  real(dp),                 intent(in)    :: thermal_energy
!  real(dp),                 intent(in)    :: min_frequency
!  type(StructureData),      intent(in)    :: harmonic_supercell
!  type(CartesianHessian),   intent(in)    :: harmonic_hessian
!  type(MinImages),          intent(in)    :: harmonic_min_images(:,:)
!  type(DegenerateSubspace), intent(in)    :: subspaces(:)
!  class(SubspaceBasis),     intent(in)    :: subspace_bases(:)
!  class(BasisStates),       intent(inout) :: subspace_states(:)
!  type(MinImages),          intent(in)    :: anharmonic_min_images(:,:)
!  type(ConvergenceData),    intent(in)    :: convergence_data
!  type(AnharmonicData),     intent(in)    :: anharmonic_data
!  type(PhononDispersion)                  :: output
!  
!  ! Variables for calculating dispersion curve.
!  type(QpointPath)                         :: path
!  type(PathQpointAndFraction), allocatable :: path_qpoints(:)
!  
!  ! Variables at a given q-point.
!  type(RealVector)                      :: qpoint
!  type(DynamicalMatrix)                 :: dynamical_matrix
!  type(ComplexMode),        allocatable :: qpoint_modes(:)
!  type(ComplexMode),        allocatable :: modes(:)
!  type(DegenerateSubspace), allocatable :: fine_subspaces(:)
!  
!  ! Variables in a given subspace.
!  type(ComplexMode),        allocatable :: subspace_modes(:)
!  type(PotentialPointer)                :: fine_potential
!  type(DegenerateSubspace), allocatable :: split_subspaces(:)
!  real(dp)                              :: starting_frequency
!  type(HarmonicBasis)                   :: vscha_basis
!  type(HarmonicStates)                  :: vscha_states
!  real(dp)                              :: vscha_frequency
!  
!  ! Output variables.
!  real(dp),              allocatable :: frequencies(:)
!  type(PathFrequencies), allocatable :: path_frequencies(:)
!  
!  integer :: i,j,k,ialloc
!  
!  path = QpointPath(path_string, no_path_points)
!  path_qpoints = path%path_qpoints()
!  
!  allocate( frequencies(3*anharmonic_data%structure%no_atoms), &
!          & path_frequencies(size(path_qpoints)),              &
!          & stat=ialloc); call err(ialloc)
!  do i=1,size(path_qpoints)
!    qpoint = path_qpoints(i)%qpoint
!    
!    ! Construct normal modes from harmonic potential.
!    ! Include modes from both q and -q.
!    dynamical_matrix = DynamicalMatrix( qpoint,             &
!                                      & harmonic_supercell, &
!                                      & harmonic_hessian,   &
!                                      & harmonic_min_images )
!    
!    qpoint_modes = ComplexMode(dynamical_matrix, harmonic_supercell)
!    qpoint_modes%id = [(j,j=1,size(qpoint_modes))]
!    qpoint_modes%paired_id = [(j+size(qpoint_modes),j=1,size(qpoint_modes))]
!    
!    modes = [( [qpoint_modes(j),conjg(qpoint_modes(j))], &
!             & j=1,                                      &
!             & size(qpoint_modes)                        )]
!    
!    ! Split modes into subspaces by harmonic frequency.
!    fine_subspaces = generate_fine_subspaces( modes,               &
!                                            & degenerate_frequency )
!    k = 0
!    do j=1,size(fine_subspaces)
!      ! Generate the modes in the subspace.
!      subspace_modes = fine_subspaces(j)%modes(modes)
!      
!      ! Interpolate the potential.
!      fine_potential = potential%interpolate( &
!                     & qpoint,                &
!                     & fine_subspaces(j),     &
!                     & subspace_modes,        &
!                     & anharmonic_min_images, &
!                     & thermal_energy,        &
!                     & subspaces,             &
!                     & subspace_bases,        &
!                     & subspace_states,       &
!                     & anharmonic_data        )
!      
!      ! Fit VSCHA frequency.
!      starting_frequency = max( fine_subspaces(j)%frequency,                  &
!                              & anharmonic_data%frequency_of_max_displacement )
!      vscha_basis = HarmonicBasis( fine_subspaces(j)%id, &
!                                 & starting_frequency    )
!      vscha_states = HarmonicStates(                          &
!         & vscha_basis%calculate_states( fine_subspaces(j),   &
!         &                               fine_potential,      &
!         &                               thermal_energy,      &
!         &                               convergence_data,    &
!         &                               anharmonic_data    ) )
!      vscha_frequency = vscha_states%frequency
!      
!      ! Record frequencies.
!      frequencies(k+1:k+size(subspace_modes)/2) = vscha_frequency
!      k = k + size(subspace_modes)/2
!    enddo
!    
!    path_frequencies(i) = PathFrequencies(              &
!       & path_fraction = path_qpoints(i)%path_fraction, &
!       & qpoint        = qpoint,                        &
!       & frequencies   = frequencies,                   &
!       & dynamical_matrix = dynamical_matrix            )
!    
!    if (modulo(i,size(path_qpoints)/10)==0) then
!      call print_line('VSCF Dispersion: '//i//' of '// &
!         & size(path_qpoints)//' q-points sampled.')
!    endif
!  enddo
!  
!  output = PhononDispersion(path, path_frequencies)
!end function
!
!function calculate_interpolated_thermodynamics(potential,stress,            &
!   & degenerate_frequency,harmonic_dos,thermal_energy,min_frequency,        &
!   & harmonic_supercell,harmonic_hessian,harmonic_min_images,subspaces,     &
!   & subspace_bases,subspace_states,anharmonic_min_images,convergence_data, &
!   & anharmonic_data) result(output)
!  implicit none
!  
!  class(PotentialData),     intent(in)           :: potential
!  class(StressData),        intent(in), optional :: stress
!  real(dp),                 intent(in)           :: degenerate_frequency
!  type(PhononDos),          intent(in)           :: harmonic_dos
!  real(dp),                 intent(in)           :: thermal_energy
!  real(dp),                 intent(in)           :: min_frequency
!  type(StructureData),      intent(in)           :: harmonic_supercell
!  type(CartesianHessian),   intent(in)           :: harmonic_hessian
!  type(MinImages),          intent(in)           :: harmonic_min_images(:,:)
!  type(DegenerateSubspace), intent(in)           :: subspaces(:)
!  class(SubspaceBasis),     intent(in)           :: subspace_bases(:)
!  class(BasisStates),       intent(inout)        :: subspace_states(:)
!  type(MinImages),          intent(in)           :: anharmonic_min_images(:,:)
!  type(ConvergenceData),    intent(in)           :: convergence_data
!  type(AnharmonicData),     intent(in)           :: anharmonic_data
!  type(InterpolatedThermodynamics)               :: output
!  
!  ! q-points for fine interpolation.
!  type(RealVector), allocatable :: fine_qpoints(:)
!  
!  ! Parameters extracted from harmonic DOS.
!  integer  :: no_bins
!  real(dp) :: min_freq
!  real(dp) :: max_freq
!  real(dp) :: bin_width
!  
!  ! Variables at a given q-point.
!  type(DynamicalMatrix)                 :: dynamical_matrix
!  type(ComplexMode),        allocatable :: qpoint_modes(:)
!  type(ComplexMode),        allocatable :: modes(:)
!  type(DegenerateSubspace), allocatable :: fine_subspaces(:)
!  
!  ! Variables in a given subspace.
!  type(ComplexMode),        allocatable :: subspace_modes(:)
!  type(PotentialPointer)                :: fine_potential
!  type(DegenerateSubspace), allocatable :: split_subspaces(:)
!  real(dp)                              :: starting_frequency
!  type(HarmonicBasis)                   :: vscha_basis
!  type(HarmonicStates)                  :: vscha_states
!  real(dp)                              :: vscha_frequency
!  
!  ! Variables for working with frequency bins.
!  integer :: bin
!  integer :: extra_bins
!  
!  ! Output variables.
!  real(dp),            allocatable :: freq_dos(:)
!  type(ThermodynamicData)          :: thermodynamics
!  type(SampledQpoint), allocatable :: sampled_qpoints(:)
!  type(PdosBin),       allocatable :: pdos(:)
!  type(PhononDos)                  :: phonon_dos
!  
!  ! Temporary variables.
!  integer :: i,j,k
!  
!  fine_qpoints = harmonic_dos%qpoints%qpoint
!  
!  no_bins = size(harmonic_dos%pdos)
!  min_freq = harmonic_dos%pdos(1)%min_frequency
!  max_freq = harmonic_dos%pdos(no_bins)%max_frequency
!  bin_width = (max_freq-min_freq) / no_bins
!  freq_dos = [(0.0_dp, j=1, no_bins)]
!  
!  thermodynamics = ThermodynamicData(thermal_energy, 0.0_dp, 0.0_dp, 0.0_dp)
!  
!  do i=1,size(fine_qpoints)
!    ! Construct normal modes from harmonic potential.
!    ! Include modes from both q and -q.
!    dynamical_matrix = DynamicalMatrix( fine_qpoints(i),    &
!                                      & harmonic_supercell, &
!                                      & harmonic_hessian,   &
!                                      & harmonic_min_images )
!    
!    qpoint_modes = ComplexMode(dynamical_matrix, harmonic_supercell)
!    qpoint_modes%id = [(j,j=1,size(qpoint_modes))]
!    qpoint_modes%paired_id = [(j+size(qpoint_modes),j=1,size(qpoint_modes))]
!    
!    modes = [( [qpoint_modes(j),conjg(qpoint_modes(j))], &
!             & j=1,                                      &
!             & size(qpoint_modes)                        )]
!    
!    ! Split modes into subspaces by harmonic frequency.
!    fine_subspaces = generate_fine_subspaces( modes,               &
!                                            & degenerate_frequency )
!    do j=1,size(fine_subspaces)
!      ! Generate the modes in the subspace.
!      subspace_modes = fine_subspaces(j)%modes(modes)
!      
!      ! Interpolate the potential.
!      fine_potential = potential%interpolate( &
!                     & fine_qpoints(i),       &
!                     & fine_subspaces(j),     &
!                     & subspace_modes,        &
!                     & anharmonic_min_images, &
!                     & thermal_energy,        &
!                     & subspaces,             &
!                     & subspace_bases,        &
!                     & subspace_states,       &
!                     & anharmonic_data        )
!      
!      ! Fit VSCHA frequency.
!      starting_frequency = max( fine_subspaces(j)%frequency,                  &
!                              & anharmonic_data%frequency_of_max_displacement )
!      vscha_basis = HarmonicBasis( subspace_id    = fine_subspaces(j)%id, &
!                                 & frequency      = starting_frequency,   &
!                                 & supercell_size = size(fine_qpoints)    )
!      vscha_states = HarmonicStates(                          &
!         & vscha_basis%calculate_states( fine_subspaces(j),   &
!         &                               fine_potential,      &
!         &                               thermal_energy,      &
!         &                               convergence_data,    &
!         &                               anharmonic_data    ) )
!      vscha_frequency = vscha_states%frequency
!      
!      bin = ceiling( (vscha_frequency-min_freq) / bin_width)
!      if (bin<1) then
!        extra_bins = 1-bin
!        min_freq = min_freq - extra_bins*bin_width
!        freq_dos = [(0.0_dp,k=1,extra_bins), freq_dos]
!        no_bins = no_bins + extra_bins
!        bin = bin + extra_bins
!      elseif (bin>no_bins) then
!        extra_bins = bin-no_bins
!        max_freq = max_freq + extra_bins*bin_width
!        freq_dos = [freq_dos, (0.0_dp,k=1,extra_bins)]
!        no_bins = no_bins + extra_bins
!      endif
!      freq_dos(bin) = freq_dos(bin) + size(subspace_modes)/2
!      
!      thermodynamics = thermodynamics                                     &
!                   & + ThermodynamicData(thermal_energy, vscha_frequency) &
!                   & * size(subspace_modes)/2
!    enddo
!    
!    if (modulo(i,size(fine_qpoints)/10)==0) then
!      call print_line('VSCF Thermodynamics: '//i//' of '// &
!         & size(fine_qpoints)//' q-points sampled.')
!    endif
!  enddo
!  
!  thermodynamics = thermodynamics / size(fine_qpoints)
!  
!  sampled_qpoints = SampledQpoint(fine_qpoints, 0)
!  
!  pdos = [( PdosBin( min_frequency = min_freq+(i-1)*bin_width,    &
!          &          max_frequency = min_freq+i*bin_width,        &
!          &          occupation    = freq_dos(i)               ), &
!          & i=1,                                                  &
!          & no_bins                                               )]
!  
!  phonon_dos = PhononDos(sampled_qpoints, pdos, [thermodynamics])
!  
!  ! Construct output.
!  output = InterpolatedThermodynamics(       &
!     & vscha_dos = phonon_dos,               &
!     & vscha_thermodynamics = thermodynamics )
!end function
!
!function generate_fine_subspaces(modes,degenerate_frequency) result(output)
!  implicit none
!  
!  type(ComplexMode), intent(in)         :: modes(:)
!  real(dp),          intent(in)         :: degenerate_frequency
!  type(DegenerateSubspace), allocatable :: output(:)
!  
!  type(ComplexMode), allocatable :: qpoint_modes(:)
!  integer,           allocatable :: splits(:)
!  integer                        :: no_splits
!  
!  integer :: i,ialloc
!  
!  qpoint_modes = modes(filter(modes%id<=modes%paired_id))
!  
!  allocate(splits(size(qpoint_modes)+1), stat=ialloc); call err(ialloc)
!  no_splits = 1
!  splits(no_splits) = 1
!  do i=2,size(qpoint_modes)
!    if ( abs(qpoint_modes(i)%frequency-qpoint_modes(i-1)%frequency) &
!     & >= degenerate_frequency                                      ) then
!      no_splits = no_splits+1
!      splits(no_splits) = i
!    endif
!  enddo
!  no_splits = no_splits+1
!  splits(no_splits) = size(qpoint_modes)+1
!  
!  allocate(output(no_splits-1), stat=ialloc); call err(ialloc)
!  do i=1,no_splits-1
!    associate(subspace_modes=>qpoint_modes(splits(i):splits(i+1)-1))
!      output(i) = DegenerateSubspace(                                  &
!         & id         = i,                                             &
!         & frequency  = subspace_modes(1)%frequency,                   &
!         & mode_ids   = [subspace_modes%id, subspace_modes%paired_id], &
!         & paired_ids = [subspace_modes%paired_id, subspace_modes%id]  )
!    end associate
!  enddo
!end function
end module
