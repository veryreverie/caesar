submodule (caesar_stress_data_module) caesar_stress_data_submodule
  use caesar_anharmonic_common_module
  
  ! An array of all types which extend StressData.
  ! This array will be filled in by startup routines.
  type(StressPointer), allocatable :: TYPES_StressData(:)
contains

module procedure startup_StressData
  type(StressBasePointer) :: base
  
  integer :: i
  
  ! TODO: this doesn't work. It calls the startup() method for
  !    StressBasePointer only.
  base = StressBasePointer(this)
  call base%startup()
  
  if (.not.allocated(TYPES_StressData)) then
    TYPES_StressData = [StressPointer(this)]
  elseif (.not.any([(                                      &
     & this%representation()                               &
     &    == TYPES_StressData(i)%stress_%representation(), &
     & i=1,                                                &
     & size(TYPES_StressData)                              )])) then
    TYPES_StressData = [TYPES_StressData, StressPointer(this)]
  endif
end procedure

module procedure new_StressPointer
  integer :: ialloc
  
  select type(stress); type is(StressPointer)
    this = stress
  class default
    this%representation_ = stress%representation()
    allocate( this%stress_, source=stress, &
            & stat=ialloc); call err(ialloc)
  end select
end procedure

module procedure check_StressPointer
  if (.not. allocated(this%stress_)) then
    call print_line(CODE_ERROR//': Trying to use a StressPointer before &
       &it has been allocated.')
    call err()
  endif
end procedure

module procedure representation_StressPointer
  output = 'pointer'
end procedure

module procedure zero_stress_StressPointer
  call this%check()
  
  call this%stress_%zero_stress()
end procedure

module procedure add_constant_StressPointer
  call this%check()
  
  call this%stress_%add_constant(input)
end procedure

module procedure stress_RealModeDisplacement_StressPointer
  call this%check()
  
  output = this%stress_%stress(displacement)
end procedure

module procedure stress_ComplexModeDisplacement_StressPointer
  call this%check()
  
  output = this%stress_%stress(displacement)
end procedure

module procedure braket_SubspaceBraKet_StressPointer
  call this%check()
  
  call this%stress_%braket(braket,whole_subspace,anharmonic_data)
end procedure

module procedure braket_BasisState_StressPointer
  call this%check()
  
  call this%stress_%braket( bra,            &
                          & ket,            &
                          & subspace,       &
                          & subspace_basis, &
                          & whole_subspace, &
                          & anharmonic_data )
end procedure

module procedure braket_BasisStates_StressPointer
  call this%check()
  
  call this%stress_%braket( states,         &
                          & subspace,       &
                          & subspace_basis, &
                          & whole_subspace, &
                          & anharmonic_data )
end procedure

module procedure harmonic_expectation_StressPointer
  call this%check()
  
  output = this%stress_%harmonic_expectation( frequency,      &
                                            & thermal_energy, &
                                            & supercell_size, &
                                            & anharmonic_data )
end procedure

module procedure potential_stress_SubspaceBraKet_StressPointer
  call this%check()
  
  output = this%stress_%potential_stress(braket, anharmonic_data)
end procedure

module procedure potential_stress_BasisState_StressPointer
  call this%check()
  
  output = this%stress_%potential_stress( bra,            &
                                        & ket,            &
                                        & subspace,       &
                                        & subspace_basis, &
                                        & anharmonic_data )
end procedure

module procedure can_be_interpolated_StressPointer
  call this%check()
  
  output = this%stress_%can_be_interpolated()
end procedure

module procedure calculate_dynamical_matrices_StressPointer
  call this%check()
  
  output = this%stress_%calculate_dynamical_matrices( qpoints,             &
                                                    & thermal_energy,      &
                                                    & subspaces,           &
                                                    & subspace_bases,      &
                                                    & subspace_states,     &
                                                    & anharmonic_data      )
end procedure

module procedure stress_correction_StressPointer
  call this%check()
  
  output = this%stress_%stress_correction( subspaces,       &
                                         & subspace_bases,  &
                                         & subspace_states, &
                                         & anharmonic_data  )
end procedure

module procedure can_be_interpolated_StressData
  output = .false.
end procedure

module procedure calculate_dynamical_matrices_StressData
  ! This should be gated behind can_be_interpolated.
  call print_line(CODE_ERROR//': calculate_dynamical_matrices not &
     &implemented for this stress.')
  call err()
end procedure

module procedure stress_correction_StressData
  ! This should be gated behind can_be_interpolated.
  call print_line(CODE_ERROR//': stress_correction not implemented for this &
     &potential.')
  call err()
end procedure

module procedure read_StressPointer
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(StressPointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    i = first([(                                                       &
       & TYPES_StressData(i)%stress_%representation()==representation, &
       & i=1,                                                          &
       & size(TYPES_StressData)                                        )])
    
    ! Read the input into the element of the correct type,
    !    and copy that element into the output.
    call TYPES_StressData(i)%stress_%read(input(2:))
    this = StressPointer(TYPES_StressData(i))
  class default
    call err()
  end select
end procedure

module procedure write_StressPointer
  select type(this); type is(StressPointer)
    output = [ 'Stress representation: '//this%representation_, &
             & str(this%stress_)                                ]
  end select
end procedure

module procedure new_StressPointer_Strings
  call this%read(input)
end procedure

module procedure new_StressPointer_StringArray
  this = StressPointer(str(input))
end procedure
end submodule
