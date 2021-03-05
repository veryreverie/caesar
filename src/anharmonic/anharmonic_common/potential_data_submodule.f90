submodule (caesar_potential_data_module) caesar_potential_data_submodule
  use caesar_anharmonic_common_module
contains

module procedure new_PotentialPointer
  integer :: ialloc
  
  select type(potential); type is(PotentialPointer)
    this = potential
  class default
    this%representation_ = potential%representation()
    allocate( this%potential_, source=potential, &
            & stat=ialloc); call err(ialloc)
  end select
end procedure

module procedure check_PotentialPointer
  if (.not. allocated(this%potential_)) then
    call print_line(CODE_ERROR//': Trying to use a PotentialPointer before &
       &it has been allocated.')
    call err()
  endif
end procedure

module procedure representation_PotentialPointer
  output = 'pointer'
end procedure

module procedure potential_PotentialPointer
  call this%check()
  
  output = this%potential_
end procedure

module procedure generate_sampling_points_PotentialPointer
  call this%check()
  
  call this%potential_%generate_sampling_points( anharmonic_data,       &
                                               & use_forces,            &
                                               & energy_to_force_ratio, &
                                               & use_hessians,          &
                                               & calculate_stress,      &
                                               & sampling_points_dir,   &
                                               & calculation_writer,    &
                                               & logfile                )
end procedure

module procedure generate_potential_PotentialPointer
  call this%check()
  
  call this%potential_%generate_potential( anharmonic_data,              &
                                         & weighted_energy_force_ratio,  &
                                         & sampling_points_dir,          &
                                         & calculation_reader,           &
                                         & logfile                       )
end procedure

module procedure generate_stress_PotentialPointer
  call this%check()
  
  output = this%potential_%generate_stress( anharmonic_data,           &
                                          & sampling_points_dir,       &
                                          & stress_expansion_order,    &
                                          & stress_subspace_coupling,  &
                                          & vscf_basis_functions_only, &
                                          & calculation_reader,        &
                                          & logfile                    )
end procedure

module procedure zero_energy_PotentialPointer
  call this%check()
  
  call this%potential_%zero_energy()
end procedure

module procedure add_constant_PotentialPointer
  call this%check()
  
  call this%potential_%add_constant(input)
end procedure

module procedure optimise_subspace_potential_PotentialPointer
  call this%check()
  
  call this%potential_%optimise_subspace_potential( subspace,               &
                                                  & subspace_basis,         &
                                                  & old_subspace_potential, &
                                                  & anharmonic_data         )
end procedure

module procedure energy_RealModeDisplacement_PotentialPointer
  call this%check()
  
  output = this%potential_%energy(displacement)
end procedure

module procedure energy_ComplexModeDisplacement_PotentialPointer
  call this%check()
  
  output = this%potential_%energy(displacement)
end procedure

module procedure force_RealModeDisplacement_PotentialPointer
  call this%check()
  
  output = this%potential_%force(displacement)
end procedure

module procedure force_ComplexModeDisplacement_PotentialPointer
  call this%check()
  
  output = this%potential_%force(displacement)
end procedure

module procedure braket_SubspaceBraKet_PotentialPointer
  call this%check()
  
  call this%potential_%braket(braket,whole_subspace,anharmonic_data)
end procedure

module procedure braket_BasisState_PotentialPointer
  call this%check()
  
  call this%potential_%braket( bra,            &
                             & ket,            &
                             & subspace,       &
                             & subspace_basis, &
                             & whole_subspace, &
                             & anharmonic_data )
end procedure

module procedure braket_BasisStates_PotentialPointer
  call this%check()
  
  call this%potential_%braket( states,         &
                             & subspace,       &
                             & subspace_basis, &
                             & whole_subspace, &
                             & anharmonic_data )
end procedure

module procedure harmonic_expectation_PotentialPointer
  call this%check()
  
  output = this%potential_%harmonic_expectation( frequency,      &
                                               & thermal_energy, &
                                               & supercell_size, &
                                               & anharmonic_data )
end procedure

module procedure potential_energy_SubspaceBraKet_PotentialPointer
  call this%check()
  
  output = this%potential_%potential_energy(braket, anharmonic_data)
end procedure

module procedure potential_energy_BasisState_PotentialPointer
  call this%check()
  
  output = this%potential_%potential_energy( bra,            &
                                           & ket,            &
                                           & subspace,       &
                                           & subspace_basis, &
                                           & anharmonic_data )
end procedure

module procedure coefficients_PotentialPointer
  call this%check()
  
  output = this%potential_%coefficients(anharmonic_data)
end procedure

module procedure set_coefficients_PotentialPointer
  call this%check()
  
  call this%potential_%set_coefficients( coefficients,   &
                                       & anharmonic_data )
end procedure

module procedure all_basis_functions_PotentialPointer
  call this%check()
  
  output = this%potential_%all_basis_functions(anharmonic_data)
end procedure

module procedure variable_basis_functions_PotentialPointer
  call this%check()
  
  output = this%potential_%variable_basis_functions(anharmonic_data)
end procedure

module procedure can_be_interpolated_PotentialPointer
  call this%check()
  
  output = this%potential_%can_be_interpolated()
end procedure

module procedure interpolate_potential_PotentialPointer
  call this%check()
  
  call this%potential_%interpolate_potential( anharmonic_min_images,         &
                                            & potential,                     &
                                            & anharmonic_data,               &
                                            & interpolated_anharmonic_data,  &
                                            & difference_dynamical_matrices, &
                                            & logfile                        )
end procedure

module procedure calculate_dynamical_matrices_PotentialPointer
  call this%check()
  
  output = this%potential_%calculate_dynamical_matrices( qpoints,             &
                                                       & thermal_energy,      &
                                                       & subspaces,           &
                                                       & subspace_bases,      &
                                                       & subspace_states,     &
                                                       & anharmonic_data      )
end procedure

module procedure energy_correction_PotentialPointer
  call this%check()
  
  output = this%potential_%energy_correction( subspaces,       &
                                            & subspace_bases,  &
                                            & subspace_states, &
                                            & anharmonic_data  )
end procedure

module procedure can_be_interpolated_PotentialData
  output = .false.
end procedure

module procedure interpolate_potential_PotentialData
  ! This should be gated behind can_be_interpolated.
  call print_line(CODE_ERROR//': interpolate_potential not implemented for &
     &this potential.')
  call err()
end procedure

module procedure calculate_dynamical_matrices_PotentialData
  ! This should be gated behind can_be_interpolated.
  call print_line(CODE_ERROR//': calculate_dynamical_matrices not &
     &implemented for this potential.')
  call err()
end procedure

module procedure energy_correction_PotentialData
  ! This should be gated behind can_be_interpolated.
  call print_line(CODE_ERROR//': energy_correction not implemented for this &
     &potential.')
  call err()
end procedure

module procedure optimise_subspace_potential_PotentialData
  ! By default this doesn't do anything.
end procedure

module procedure read_PotentialPointer
  type(PotentialPointer), allocatable :: types(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(PotentialPointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    types = types_PotentialData()
    i = first([(                                               &
       & types(i)%potential_%representation()==representation, &
       & i=1,                                                  &
       & size(types)                                           )])
    this = types(i)
    
    call this%potential_%read(input(2:))
  class default
    call err()
  end select
end procedure

module procedure write_PotentialPointer
  select type(this); type is(PotentialPointer)
    output = [ 'Potential representation: '//this%representation_, &
             & str(this%potential_)                                ]
  end select
end procedure

module procedure new_PotentialPointer_Strings
  call this%read(input)
end procedure

module procedure new_PotentialPointer_StringArray
  this = PotentialPointer(str(input))
end procedure
end submodule
