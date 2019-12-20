! ======================================================================
! Interpolate a potential or stress to a fine grid of q-points,
!    and calculate the thermodynamic variables.
! ======================================================================
module interpolation_module
  use common_module
  use anharmonic_common_module
  implicit none
  
  private
  
  public :: InterpolatedThermodynamics
  public :: calculate_interpolated_stress
  public :: calculate_interpolated_thermodynamics
  
  type, extends(NoDefaultConstructor) :: InterpolatedThermodynamics
  end type
contains

function calculate_interpolated_stress(stress,degenerate_frequency, &
   & fine_qpoints,thermal_energy,min_frequency,harmonic_supercell,  &
   & harmonic_hessian,harmonic_min_images,subspaces,subspace_bases, &
   & basis_states,anharmonic_min_images,anharmonic_data) result(output)
  implicit none
  
  class(StressData),        intent(in) :: stress
  real(dp),                 intent(in) :: degenerate_frequency
  type(RealVector),         intent(in) :: fine_qpoints(:)
  real(dp),                 intent(in) :: thermal_energy
  real(dp),                 intent(in) :: min_frequency
  type(StructureData),      intent(in) :: harmonic_supercell
  type(CartesianHessian),   intent(in) :: harmonic_hessian
  type(MinImages),          intent(in) :: harmonic_min_images(:,:)
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  class(SubspaceBasis),     intent(in) :: subspace_bases(:)
  class(BasisStates),       intent(in) :: basis_states(:)
  type(MinImages),          intent(in) :: anharmonic_min_images(:,:)
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(RealMatrix)                     :: output
  
  ! Variables at a given q-point.
  type(DynamicalMatrix)                 :: dynamical_matrix
  type(ComplexMode),        allocatable :: qpoint_modes(:)
  type(ComplexMode),        allocatable :: modes(:)
  type(DegenerateSubspace), allocatable :: fine_subspaces(:)
  
  ! Variables in a given subspace.
  type(ComplexMode), allocatable :: subspace_modes(:)
  type(StressPrefactors)         :: stress_prefactors
  type(StressPointer)            :: fine_stress
  type(RealMatrix)               :: potential_stress
  type(ThermodynamicData)        :: thermodynamic_data
  
  ! Temporary variables
  integer :: i,j,k
  
  output = dblemat(zeroes(3,3))
  do i=1,size(fine_qpoints)
    ! Construct normal modes from harmonic potential.
    ! Include modes from both q and -q.
    dynamical_matrix = DynamicalMatrix( fine_qpoints(i),    &
                                      & harmonic_supercell, &
                                      & harmonic_hessian,   &
                                      & harmonic_min_images )
    
    qpoint_modes = ComplexMode(dynamical_matrix, harmonic_supercell)
    qpoint_modes%id = [(j,j=1,size(qpoint_modes))]
    qpoint_modes%paired_id = [(j+size(qpoint_modes),j=1,size(qpoint_modes))]
    
    modes = [( [qpoint_modes(j),conjg(qpoint_modes(j))], &
             & j=1,                                      &
             & size(qpoint_modes)                        )]
    
    ! Split modes into subspaces.
    fine_subspaces = generate_fine_subspaces( qpoint_modes,        &
                                            & degenerate_frequency )
    
    do j=1,size(fine_subspaces)
      if (fine_subspaces(j)%frequency>min_frequency) then
        ! Generate the modes in the subspace.
        subspace_modes = fine_subspaces(j)%modes(modes)
        
        ! Calculate stress prefactors.
        stress_prefactors = StressPrefactors(fine_subspaces(j), subspace_modes)
        
        ! Interpolate stress.
        fine_stress = stress%interpolate( fine_qpoints(i),       &
                                        & fine_subspaces(j),     &
                                        & subspace_modes,        &
                                        & anharmonic_min_images, &
                                        & anharmonic_data        )
        
        ! Take the harmonic expectation of the stress.
        potential_stress = fine_stress%harmonic_expectation(   &
                       &        fine_subspaces(j)%frequency,   &
                       &        thermal_energy,                &
                       &        anharmonic_data              ) &
                       & / (2*size(subspace_modes))
        
        ! Calculate total stress, including kinetic stress.
        thermodynamic_data =                                             &
           & ThermodynamicData( thermal_energy,                          &
           &                    fine_subspaces(j)%frequency,             &
           &                    stress_prefactors%average_prefactor(),   &
           &                    potential_stress,                        &
           &                    anharmonic_data%structure%volume       ) &
           & * size(subspace_modes)
        
        output = output + thermodynamic_data%stress
      endif
    enddo
    
    if (modulo(i,size(fine_qpoints)/10)==0) then
      call print_line('Stress: '//i//' of '//size(fine_qpoints)// &
         & ' q-points sampled.')
    endif
  enddo
  
  ! Normalise to be per primitive cell.
  ! N.B. at each q-point, q and -q are both considered, hence the factor of 2.
  output = output / (2*size(fine_qpoints))
  
  ! Add in the static-lattice stress.
  ! TODO: Correct VSCF form for when stress is not harmonic.
  !output = output + stress%undisplaced_stress() &
  !              & / anharmonic_data%anharmonic_supercell%sc_size
end function

function calculate_interpolated_thermodynamics(potential,stress,           &
   & degenerate_frequency,fine_qpoints,thermal_energy,min_frequency,       &
   & harmonic_supercell,harmonic_hessian,harmonic_min_images,subspaces,    &
   & subspace_bases,subspace_states,anharmonic_min_images,anharmonic_data) &
   & result(output)
  implicit none
  
  class(PotentialData),     intent(in)           :: potential
  class(StressData),        intent(in), optional :: stress
  real(dp),                 intent(in)           :: degenerate_frequency
  type(RealVector),         intent(in)           :: fine_qpoints(:)
  real(dp),                 intent(in)           :: thermal_energy
  real(dp),                 intent(in)           :: min_frequency
  type(StructureData),      intent(in)           :: harmonic_supercell
  type(CartesianHessian),   intent(in)           :: harmonic_hessian
  type(MinImages),          intent(in)           :: harmonic_min_images(:,:)
  type(DegenerateSubspace), intent(in)           :: subspaces(:)
  class(SubspaceBasis),     intent(in)           :: subspace_bases(:)
  class(BasisStates),       intent(inout)        :: subspace_states(:)
  type(MinImages),          intent(in)           :: anharmonic_min_images(:,:)
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(InterpolatedThermodynamics)               :: output
  
  ! Variables at a given q-point.
  type(DynamicalMatrix)                 :: dynamical_matrix
  type(ComplexMode),        allocatable :: qpoint_modes(:)
  type(ComplexMode),        allocatable :: modes(:)
  type(DegenerateSubspace), allocatable :: fine_subspaces(:)
  
  ! Variables in a given subspace.
  type(ComplexMode),        allocatable :: subspace_modes(:)
  type(PotentialPointer)                :: fine_potential
  type(DegenerateSubspace), allocatable :: split_subspaces(:)
  
  ! Temporary variables.
  integer :: i,j
  
  do i=1,size(fine_qpoints)
    ! Construct normal modes from harmonic potential.
    ! Include modes from both q and -q.
    dynamical_matrix = DynamicalMatrix( fine_qpoints(i),    &
                                      & harmonic_supercell, &
                                      & harmonic_hessian,   &
                                      & harmonic_min_images )
    
    qpoint_modes = ComplexMode(dynamical_matrix, harmonic_supercell)
    qpoint_modes%id = [(j,j=1,size(qpoint_modes))]
    qpoint_modes%paired_id = [(j+size(qpoint_modes),j=1,size(qpoint_modes))]
    
    modes = [( [qpoint_modes(j),conjg(qpoint_modes(j))], &
             & j=1,                                      &
             & size(qpoint_modes)                        )]
    
    ! Split modes into subspaces by harmonic frequency.
    fine_subspaces = generate_fine_subspaces( modes,               &
                                            & degenerate_frequency )
    do j=1,size(fine_subspaces)
      ! Generate the modes in the subspace.
      subspace_modes = fine_subspaces(j)%modes(modes)
      
      ! Interpolate the potential.
      fine_potential = potential%interpolate( &
                     & fine_qpoints(i),       &
                     & fine_subspaces(j),     &
                     & subspace_modes,        &
                     & anharmonic_min_images, &
                     & thermal_energy,        &
                     & subspaces,             &
                     & subspace_bases,        &
                     & subspace_states,       &
                     & anharmonic_data        )
      
      ! Split the subspace into individual-mode subspaces.
      split_subspaces = generate_fine_subspaces( subspace_modes,             &
                                               & degenerate_frequency=0.0_dp )
    enddo
    
    if (modulo(i,size(fine_qpoints)/10)==0) then
      call print_line('Stress: '//i//' of '//size(fine_qpoints)// &
         & ' q-points sampled.')
    endif
  enddo
end function

function generate_fine_subspaces(modes,degenerate_frequency) result(output)
  implicit none
  
  type(ComplexMode), intent(in)         :: modes(:)
  real(dp),          intent(in)         :: degenerate_frequency
  type(DegenerateSubspace), allocatable :: output(:)
  
  type(ComplexMode), allocatable :: qpoint_modes(:)
  integer,           allocatable :: splits(:)
  integer                        :: no_splits
  
  integer :: i,ialloc
  
  qpoint_modes = modes(filter(modes%id<=modes%paired_id))
  
  allocate(splits(size(qpoint_modes)+1), stat=ialloc); call err(ialloc)
  no_splits = 1
  splits(no_splits) = 1
  do i=2,size(qpoint_modes)
    if ( abs(qpoint_modes(i)%frequency-qpoint_modes(i-1)%frequency) &
     & >= degenerate_frequency                                      ) then
      no_splits = no_splits+1
      splits(no_splits) = i
    endif
  enddo
  no_splits = no_splits+1
  splits(no_splits) = size(qpoint_modes)+1
  
  allocate(output(no_splits-1), stat=ialloc); call err(ialloc)
  do i=1,no_splits-1
    associate(subspace_modes=>qpoint_modes(splits(i):splits(i+1)-1))
      output(i) = DegenerateSubspace(                                  &
         & id         = i,                                             &
         & frequency  = subspace_modes(1)%frequency,                   &
         & mode_ids   = [subspace_modes%id, subspace_modes%paired_id], &
         & paired_ids = [subspace_modes%paired_id, subspace_modes%id]  )
    end associate
  enddo
end function
end module
