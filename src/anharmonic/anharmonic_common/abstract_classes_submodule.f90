submodule (caesar_abstract_classes_module) caesar_abstract_classes_submodule
  use caesar_anharmonic_common_module
contains

module procedure new_SubspaceBasisPointer
  integer :: ialloc
  
  select type(basis); type is(SubspaceBasisPointer)
    this = basis
  class default
    this%representation_ = basis%representation()
    allocate( this%basis_, source=basis, &
            & stat=ialloc); call err(ialloc)
    this%frequency = this%basis_%frequency
  end select
end procedure

module procedure check_SubspaceBasisPointer
  if (.not. allocated(this%basis_)) then
    call print_line(CODE_ERROR//': Trying to use a SubspaceBasisPointer &
       &before it has been allocated.')
    call err()
  endif
end procedure

module procedure representation_SubspaceBasisPointer
  output = 'pointer'
end procedure

module procedure initial_states_SubspaceBasisPointer
  call this%check()
  
  output = this%basis_%initial_states( subspace,       &
                                     & thermal_energy, &
                                     & anharmonic_data )
end procedure

module procedure calculate_states_SubspaceBasisPointer
  call this%check()
  
  output = this%basis_%calculate_states( subspace,                &
                                       & subspace_potential,      &
                                       & thermal_energy,          &
                                       & state_energy_cutoff,     &
                                       & convergence_data,        &
                                       & anharmonic_data          )
end procedure

module procedure process_subspace_potential_SubspaceBasisPointer
  call this%check()
  
  call this%basis_%process_subspace_potential( potential,      &
                                             & states,         &
                                             & subspace,       &
                                             & anharmonic_data )
end procedure

module procedure process_subspace_stress_SubspaceBasisPointer
  call this%check()
  
  call this%basis_%process_subspace_stress( stress,         &
                                          & states,         &
                                          & subspace,       &
                                          & anharmonic_data )
end procedure

module procedure select_symmetries_SubspaceBasisPointer
  call this%check()
  
  output = this%basis_%select_symmetries(symmetries, anharmonic_data)
end procedure

module procedure mode_ids_SubspaceBasisPointer
  call this%check()
  
  output = this%basis_%mode_ids(subspace,anharmonic_data)
end procedure

module procedure paired_mode_ids_SubspaceBasisPointer
  call this%check()
  
  output = this%basis_%paired_mode_ids(subspace,anharmonic_data)
end procedure

module procedure inner_product_SubspaceBasisPointer
  call this%check()
  
  output = this%basis_%inner_product( bra,            &
                                    & ket,            &
                                    & subspace,       &
                                    & anharmonic_data )
end procedure

module procedure integrate_BasisState_SubspaceBasisPointer
  call this%check()
  
  output = this%basis_%integrate( bra,            &
                                & monomial,       &
                                & ket,            &
                                & subspace,       &
                                & anharmonic_data )
end procedure

module procedure free_energy_gradient_SubspaceBasisPointer
  call this%check()
  
  output = this%basis_%free_energy_gradient( subspace_potential,  &
                                           & basis_functions,     &
                                           & subspace,            &
                                           & states,              &
                                           & thermal_energy,      &
                                           & state_energy_cutoff, &
                                           & anharmonic_data      )
end procedure

module procedure kinetic_energy_SubspaceBasisPointer
  call this%check()
  
  output = this%basis_%kinetic_energy( bra,            &
                                     & ket,            &
                                     & subspace,       &
                                     & anharmonic_data )
end procedure

module procedure harmonic_potential_energy_SubspaceBasisPointer
  call this%check()
  
  output = this%basis_%harmonic_potential_energy( bra,            &
                                                & ket,            &
                                                & subspace,       &
                                                & anharmonic_data )
end procedure

module procedure kinetic_stress_SubspaceBasisPointer
  call this%check()
  
  output = this%basis_%kinetic_stress( bra,               &
                                     & ket,               &
                                     & subspace,          &
                                     & stress_prefactors, &
                                     & anharmonic_data    )
end procedure

module procedure thermodynamic_data_SubspaceBasisPointer
  call this%check()
  
  output = this%basis_%thermodynamic_data( thermal_energy,     &
                                         & states,             &
                                         & subspace,           &
                                         & subspace_potential, &
                                         & subspace_stress,    &
                                         & stress_prefactors,  &
                                         & anharmonic_data     )
end procedure

module procedure wavefunctions_SubspaceBasisPointer
  call this%check()
  
  output = this%basis_%wavefunctions( states,         &
                                    & subspace,       &
                                    & anharmonic_data )
end procedure

module procedure integrate_BasisStates_SubspaceBasisPointer
  call this%check()
  
  output = this%basis_%integrate( states,         &
                                & monomial,       &
                                & subspace,       &
                                & anharmonic_data )
end procedure

module procedure read_SubspaceBasisPointer
  type(SubspaceBasisPointer), allocatable :: types(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(SubspaceBasisPointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    types = types_SubspaceBasis()
    i = first([(                                           &
       & types(i)%basis_%representation()==representation, &
       & i=1,                                              &
       & size(types)                                       )])
    
    call this%basis_%read(input(2:))
  class default
    call err()
  end select
end procedure

module procedure write_SubspaceBasisPointer
  select type(this); type is(SubspaceBasisPointer)
    output = [ 'SubspaceBasis representation: '//this%representation_, &
             & str(this%basis_)                                        ]
  end select
end procedure

module procedure new_SubspaceBasisPointer_Strings
  call this%read(input)
end procedure

module procedure new_SubspaceBasisPointer_StringArray
  this = SubspaceBasisPointer(str(input))
end procedure

module procedure new_PotentialBasePointer
  integer :: ialloc
  
  select type(potential); type is(PotentialBasePointer)
    this = potential
  class default
    this%representation_ = potential%representation()
    allocate( this%potential_, source=potential, &
            & stat=ialloc); call err(ialloc)
  end select
end procedure

module procedure check_PotentialBasePointer
  if (.not. allocated(this%potential_)) then
    call print_line(CODE_ERROR//': Trying to use a PotentialBasePointer &
       &before it has been allocated.')
    call err()
  endif
end procedure

module procedure representation_PotentialBasePointer
  output = 'pointer'
end procedure

module procedure potential_PotentialBasePointer
  call this%check()
  
  output = this%potential_
end procedure

module procedure zero_energy_PotentialBasePointer
  call this%check()
  
  call this%potential_%zero_energy()
end procedure

module procedure add_constant_PotentialBasePointer
  call this%check()
  
  call this%potential_%add_constant(input)
end procedure

module procedure energy_RealModeDisplacement_PotentialBasePointer
  call this%check()
  
  output = this%potential_%energy(displacement)
end procedure

module procedure energy_ComplexModeDisplacement_PotentialBasePointer
  call this%check()
  
  output = this%potential_%energy(displacement)
end procedure

module procedure force_RealModeDisplacement_PotentialBasePointer
  call this%check()
  
  output = this%potential_%force(displacement)
end procedure

module procedure force_ComplexModeDisplacement_PotentialBasePointer
  call this%check()
  
  output = this%potential_%force(displacement)
end procedure

module procedure braket_SubspaceBraKet_PotentialBasePointer
  call this%check()
  
  call this%potential_%braket(braket,whole_subspace,anharmonic_data)
end procedure

module procedure braket_BasisState_PotentialBasePointer
  call this%check()
  
  call this%potential_%braket( bra,            &
                             & ket,            &
                             & subspace,       &
                             & subspace_basis, &
                             & whole_subspace, &
                             & anharmonic_data )
end procedure

module procedure braket_BasisStates_PotentialBasePointer
  call this%check()
  
  call this%potential_%braket( states,         &
                             & subspace,       &
                             & subspace_basis, &
                             & whole_subspace, &
                             & anharmonic_data )
end procedure

module procedure harmonic_expectation_PotentialBasePointer
  call this%check()
  
  output = this%potential_%harmonic_expectation( frequency,      &
                                               & thermal_energy, &
                                               & supercell_size, &
                                               & anharmonic_data )
end procedure

module procedure potential_energy_SubspaceBraKet_PotentialBasePointer
  call this%check()
  
  output = this%potential_%potential_energy(braket, anharmonic_data)
end procedure

module procedure potential_energy_BasisState_PotentialBasePointer
  call this%check()
  
  output = this%potential_%potential_energy( bra,            &
                                           & ket,            &
                                           & subspace,       &
                                           & subspace_basis, &
                                           & anharmonic_data )
end procedure

module procedure read_PotentialBasePointer
  type(PotentialBasePointer), allocatable :: types(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(PotentialBasePointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    types = types_PotentialBase()
    i = first([(                                               &
       & types(i)%potential_%representation()==representation, &
       & i=1,                                                  &
       & size(types)                                           )])
    
    call this%potential_%read(input(2:))
  class default
    call err()
  end select
end procedure

module procedure write_PotentialBasePointer
  select type(this); type is(PotentialBasePointer)
    output = [ 'Potential representation: '//this%representation_, &
             & str(this%potential_)                                ]
  end select
end procedure

module procedure new_PotentialBasePointer_Strings
  call this%read(input)
end procedure

module procedure new_PotentialBasePointer_StringArray
  this = PotentialBasePointer(str(input))
end procedure

module procedure new_StressBasePointer
  integer :: ialloc
  
  select type(stress); type is(StressBasePointer)
    this = stress
  class default
    this%representation_ = stress%representation()
    allocate( this%stress_, source=stress, &
            & stat=ialloc); call err(ialloc)
  end select
end procedure

module procedure check_StressBasePointer
  if (.not. allocated(this%stress_)) then
    call print_line(CODE_ERROR//': Trying to use a StressBasePointer before &
       &it has been allocated.')
    call err()
  endif
end procedure

module procedure representation_StressBasePointer
  output = 'pointer'
end procedure

module procedure zero_stress_StressBasePointer
  call this%check()
  
  call this%stress_%zero_stress()
end procedure

module procedure add_constant_StressBasePointer
  call this%check()
  
  call this%stress_%add_constant(input)
end procedure

module procedure stress_RealModeDisplacement_StressBasePointer
  call this%check()
  
  output = this%stress_%stress(displacement)
end procedure

module procedure stress_ComplexModeDisplacement_StressBasePointer
  call this%check()
  
  output = this%stress_%stress(displacement)
end procedure

module procedure braket_SubspaceBraKet_StressBasePointer
  call this%check()
  
  call this%stress_%braket(braket,whole_subspace,anharmonic_data)
end procedure

module procedure braket_BasisState_StressBasePointer
  call this%check()
  
  call this%stress_%braket( bra,            &
                          & ket,            &
                          & subspace,       &
                          & subspace_basis, &
                          & whole_subspace, &
                          & anharmonic_data )
end procedure

module procedure braket_BasisStates_StressBasePointer
  call this%check()
  
  call this%stress_%braket( states,         &
                          & subspace,       &
                          & subspace_basis, &
                          & whole_subspace, &
                          & anharmonic_data )
end procedure

module procedure harmonic_expectation_StressBasePointer
  call this%check()
  
  output = this%stress_%harmonic_expectation( frequency,      &
                                            & thermal_energy, &
                                            & supercell_size, &
                                            & anharmonic_data )
end procedure

module procedure potential_stress_SubspaceBraKet_StressBasePointer
  call this%check()
  
  output = this%stress_%potential_stress(braket, anharmonic_data)
end procedure

module procedure potential_stress_BasisState_StressBasePointer
  call this%check()
  
  output = this%stress_%potential_stress( bra,            &
                                        & ket,            &
                                        & subspace,       &
                                        & subspace_basis, &
                                        & anharmonic_data )
end procedure

module procedure read_StressBasePointer
  type(StressBasePointer), allocatable :: types(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(StressBasePointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    types = types_StressBase()
    i = first([(                                            &
       & types(i)%stress_%representation()==representation, &
       & i=1,                                               &
       & size(types)                                        )])
    
    call this%stress_%read(input(2:))
  class default
    call err()
  end select
end procedure

module procedure write_StressBasePointer
  select type(this); type is(StressBasePointer)
    output = [ 'Stress representation: '//this%representation_, &
             & str(this%stress_)                                ]
  end select
end procedure

module procedure new_StressBasePointer_Strings
  call this%read(input)
end procedure

module procedure new_StressBasePointer_StringArray
  this = StressBasePointer(str(input))
end procedure

module procedure process_subspace_potential_SubspaceBasis
  ! This does nothing by default.
end procedure

module procedure process_subspace_stress_SubspaceBasis
  ! This does nothing by default.
end procedure

module procedure select_symmetries_SubspaceBasis
  ! By default, this just returns the whole list.
  output = symmetries
end procedure

module procedure undisplaced_energy
  type(RealSingleDisplacement) :: zero_displacement(0)
  
  output = this%energy(RealModeDisplacement(zero_displacement))
end procedure

module procedure potential_energy_SubspaceBraKet_PotentialBase
  type(PotentialBasePointer), allocatable :: integrated_potential
  
  integrated_potential = PotentialBasePointer(this)
  call integrated_potential%braket( braket          = braket,         &
                                  & anharmonic_data = anharmonic_data )
  output = integrated_potential%undisplaced_energy()
end procedure

module procedure potential_energy_BasisState_PotentialBase
  type(PotentialBasePointer), allocatable :: integrated_potential
  
  integrated_potential = PotentialBasePointer(this)
  call integrated_potential%braket( bra             = bra,            &
                                  & ket             = ket,            &
                                  & subspace        = subspace,       &
                                  & subspace_basis  = subspace_basis, &
                                  & anharmonic_data = anharmonic_data )
  output = integrated_potential%undisplaced_energy()
end procedure

module procedure undisplaced_stress
  type(RealSingleDisplacement) :: zero_displacement(0)
  
  output = this%stress(RealModeDisplacement(zero_displacement))
end procedure

module procedure potential_stress_SubspaceBraKet_StressBase
  type(StressBasePointer), allocatable :: integrated_stress
  
  integrated_stress = StressBasePointer(this)
  call integrated_stress%braket( braket          = braket,         &
                               & anharmonic_data = anharmonic_data )
  output = integrated_stress%undisplaced_stress()
end procedure

module procedure potential_stress_BasisState_StressBase
  type(StressBasePointer), allocatable :: integrated_stress
  
  integrated_stress = StressBasePointer(this)
  call integrated_stress%braket( bra             = bra,            &
                               & ket             = ket,            &
                               & subspace        = subspace,       &
                               & subspace_basis  = subspace_basis, &
                               & anharmonic_data = anharmonic_data )
  output = integrated_stress%undisplaced_stress()
end procedure
end submodule
