submodule (caesar_sample_result_module) caesar_sample_result_submodule
  use caesar_polynomial_module
contains

module procedure new_SampleResult
  this%energy = energy
  this%force  = force
  if (present(stress)) then
    this%stress_ = stress
  endif
end procedure

module procedure has_stress_SampleResult
  output = allocated(this%stress_)
end procedure

module procedure stress_SampleResult
  if (this%has_stress()) then
    output = this%stress_
  else
    call print_line(ERROR//': Sample result does not contain stress.')
    call err()
  endif
end procedure

module procedure new_SampleResult_calculation
  ! Output variables.
  real(dp)                      :: energy
  type(RealModeForce)           :: force
  type(RealMatrix), allocatable :: stress
  
  ! Normalise the energy to be per anharmonic supercell.
  energy = (calculation%energy() / supercell%sc_size) &
       & * anharmonic_data%anharmonic_supercell%sc_size
  
  ! Transform the forces into normal mode co-ordinates.
  force = RealModeForce( calculation%forces(),   &
      &                  supercell,              &
      &                  real_modes,             &
      &                  qpoints               ) &
      & * real(anharmonic_data%anharmonic_supercell%sc_size,dp)
  
  ! Construct output.
  ! Make the stress extensive, and normalised to be per anharmonic supercell.
  if (calculation%has_stress()) then
    stress = calculation%stress() &
         & * anharmonic_data%anharmonic_supercell%sc_size
  endif
  
  this = SampleResult(energy,force,stress)
end procedure

module procedure new_SampleResult_calculations
  type(SampleResult), allocatable :: results(:)
  
  ! Output variables.
  real(dp)            :: energy
  type(RealModeForce) :: force
  type(RealMatrix)    :: stress
  
  ! Temporary variables.
  integer :: i,ialloc
  
  ! Normalise energies and change co-ordinates of forces.
  allocate(results(size(calculations)), stat=ialloc); call err(ialloc)
  do i=1,size(calculations)
    results(i) = SampleResult( calculations(i), &
                             & supercell,       &
                             & real_modes,      &
                             & qpoints,         &
                             & anharmonic_data  )
    
    ! Reverse the VSCF R-vector transformation.
    results(i)%force = vscf_rvectors(i)%inverse_transform( results(i)%force,  &
                                                         & real_modes,        &
                                                         & qpoints            )
  enddo
  
  ! Average over calculations.
  energy = sum(results%energy) / size(results)
  force = sum(results%force) / real(size(results),dp)
  if (all(results%has_stress())) then
    stress = sum(results%stress()) / size(results)
    this = SampleResult(energy, force, stress)
  else
    this = SampleResult(energy, force)
  endif
end procedure

module procedure construct_sample_vector
  real(dp)              :: energy
  type(RealModeForce)   :: forces
  real(dp), allocatable :: weight
  
  integer :: dims
  
  integer :: i,ialloc
  
  dims = 1+size(modes)
  
  allocate( output(size(sampling_points)*dims), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(sampling_points)
    energy = sample_results(i)%energy
    forces = sample_results(i)%force
    if (present(sample_weights)) then
      weight = sample_weights(i)
    endif
    
    ! Subtract the energy and forces from the existing potential.
    if (present(potential)) then
      energy = energy - potential%energy(sampling_points(i))
      forces = forces - potential%force(sampling_points(i))
    endif
    
    ! Construct the vector.
    output((i-1)*dims+1:i*dims) = make_sample_vector( energy,             &
                                                    & forces,             &
                                                    & modes,              &
                                                    & energy_force_ratio, &
                                                    & weight              )
  enddo
end procedure
end submodule
