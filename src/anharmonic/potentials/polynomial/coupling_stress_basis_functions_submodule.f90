submodule (caesar_coupling_stress_basis_functions_module) caesar_coupling_stress_basis_functions_submodule
  use caesar_polynomial_module
contains

module procedure new_CouplingStressBasisFunctions_empty
  integer :: ialloc
  
  this%coupling = coupling
  allocate(this%basis_functions_(0), stat=ialloc); call err(ialloc)
end procedure

module procedure new_CouplingStressBasisFunctions
  this%coupling = coupling
  this%basis_functions_ = basis_functions
end procedure

module procedure representation_CouplingStressBasisFunctions
  output = 'Polynomial stress basis functions'
end procedure

module procedure size_CouplingStressBasisFunctions
  output = size(this%basis_functions_)
end procedure

module procedure basis_functions_CouplingStressBasisFunctions
  output = this%basis_functions_
end procedure

module procedure stress_RealModeDisplacement_CouplingStressBasisFunctions
  if (size(this%basis_functions_)>0) then
    output = sum(this%basis_functions_%stress(displacement))
  else
    output = dblemat(zeroes(3,3))
  endif
end procedure

module procedure stress_ComplexModeDisplacement_CouplingStressBasisFunctions
  if (size(this%basis_functions_)>0) then
    output = sum(this%basis_functions_%stress(displacement))
  else
    output = cmplxmat(zeroes(3,3))
  endif
end procedure

module procedure braket_SubspaceBraKet_CouplingStressBasisFunctions
  integer :: i
  
  ! Check if the subspace is in this basis function's coupling.
  i = first(this%coupling%ids()==braket%subspace_id, default=0)
  if (i/=0) then
    ! If whole_subspace is .true., remove the subspace from the coupling.
    if (set_default(whole_subspace,.true.)) then
      call this%coupling%remove_subspace(i)
    endif
    
    ! Integrate across the basis function, and simplify it.
    do i=1,size(this)
      call this%basis_functions_(i)%braket( braket,                           &
                                          & anharmonic_data = anharmonic_data )
    enddo
    
    call this%basis_functions_%simplify()
  endif
end procedure

module procedure braket_BasisState_CouplingStressBasisFunctions
  integer :: i
  
  ! Check if the subspace is in this basis function's coupling.
  i = first(this%coupling%ids()==bra%subspace_id, default=0)
  if (i/=0) then
    ! If whole_subspace is .true., remove the subspace from the coupling.
    if (set_default(whole_subspace,.true.)) then
      call this%coupling%remove_subspace(i)
    endif
  
    ! Integrate across the basis function, and simplify it.
    do i=1,size(this)
      call this%basis_functions_(i)%braket( bra,                              &
                                          & ket,                              &
                                          & subspace,                         &
                                          & subspace_basis,                   &
                                          & anharmonic_data = anharmonic_data )
    enddo
    
    call this%basis_functions_%simplify()
  endif
end procedure

module procedure braket_BasisStates_CouplingStressBasisFunctions
  integer :: i,j
  
  ! Check if the subspace is in this basis function's coupling.
  i = first(this%coupling%ids()==states%subspace_id, default=0)
  if (i/=0) then
    ! If whole_subspace is .true., remove the subspace from the coupling.
    if (set_default(whole_subspace,.true.)) then
      call this%coupling%remove_subspace(i)
    endif
    
    ! Integrate across the basis function, and simplify it.
    do j=1,size(this)
      call this%basis_functions_(j)%braket( states,                           &
                                          & subspace,                         &
                                          & subspace_basis,                   &
                                          & anharmonic_data = anharmonic_data )
    enddo
    
    call this%basis_functions_%simplify()
  endif
end procedure

module procedure harmonic_expectation_CouplingStressBasisFunctions
  integer :: i
  
  if (size(this%basis_functions_)==0) then
    output = dblemat(zeroes(3,3))
  else
    output = sum([( this%basis_functions_(i)%harmonic_expectation(    &
                  &                                frequency,         &
                  &                                thermal_energy,    &
                  &                                supercell_size,    &
                  &                                anharmonic_data ), &
                  & i=1,                                              &
                  & size(this%basis_functions_)                       )])
  endif
end procedure

module procedure undisplaced_stress_CouplingStressBasisFunctions
  type(RealSingleDisplacement) :: zero_displacement(0)
  
  output = this%stress(RealModeDisplacement(zero_displacement))
end procedure

module procedure zero_stress_CouplingStressBasisFunctions
  call print_line(CODE_ERROR//': zero_energy() cannot be called for type &
     &CouplingStressBasisFunctions.')
  call err()
end procedure

module procedure add_constant_CouplingStressBasisFunctions
  call print_line(CODE_ERROR//': add_constant() cannot be called for type &
     &CouplingStressBasisFunctions.')
  call err()
end procedure

module procedure append_CouplingStressBasisFunctions
  if (this%coupling/=that%coupling) then
    call print_line(CODE_ERROR//': Appending incompatible basis functions.')
    call err()
  endif
  
  this%basis_functions_ = [this%basis_functions_, that%basis_functions_]
end procedure

module procedure generate_stress_basis_functions_SubspaceCoupling
  type(SubspaceCombination), allocatable :: subspace_combinations(:)
  
  type(StressBasisFunction), allocatable :: coupling_basis_functions(:)
  type(StressBasisFunction), allocatable :: combination_basis_functions(:)
  
  integer :: i,ialloc
  
  ! Generate the set of subspace combinations corresponding to the subspace
  !    coupling.
  ! e.g. the coupling [1,2] might have combinations [1,2], [1,1,2] and [1,2,2].
  subspace_combinations = generate_subspace_combinations( &
                 & coupling,                              &
                 & minimum_power = 1,                     &
                 & maximum_power = stress_expansion_order )
  
  ! Loop over the subspace combinations corresponding to the coupling.
  allocate(coupling_basis_functions(0), stat=ialloc); call err(ialloc)
  do i=1,size(subspace_combinations)
    combination_basis_functions = generate_stress_basis_functions( &
                                      & subspace_combinations(i),  &
                                      & maximum_coupling_order,    &
                                      & structure,                 &
                                      & complex_modes,             &
                                      & qpoints,                   &
                                      & subspaces,                 &
                                      & degenerate_symmetries,     &
                                      & vscf_basis_functions_only, &
                                      & logfile                    )
    coupling_basis_functions = [ coupling_basis_functions,    &
                               & combination_basis_functions  ]
  enddo
  
  output = CouplingStressBasisFunctions( coupling,                &
                                       & coupling_basis_functions )
end procedure

module procedure fit_coefficients_CouplingStressBasisFunctions
  type(RealMatrix) :: stress_tensor
  
  real(dp), allocatable :: a(:,:)
  real(dp), allocatable :: b(:)
  
  real(dp), allocatable :: coefficients(:)
  
  integer :: i,j,ialloc
  
  ! Check inputs are consistent.
  if (size(sampling_points)/=size(sample_results)) then
    call print_line(CODE_ERROR//': The number of sampling points does not &
       &match the number of results.')
    call err()
  endif
  
  ! Calculate the stress due to each basis function at each sampling point.
  allocate( a(9*size(sampling_points), size(this%basis_functions_)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this%basis_functions_)
    do j=1,size(sampling_points)
      stress_tensor = this%basis_functions_(i)%stress(sampling_points(j))
      a((j-1)*9+1:j*9, i) = reshape(dble(stress_tensor), [9])
    enddo
  enddo
  
  ! Calculate the stress sampled at each sampling point.
  allocate(b(9*size(sampling_points)), stat=ialloc); call err(ialloc)
  do i=1,size(sampling_points)
    stress_tensor = sample_results(i)%stress()
    
    ! Subtract the existing stress.
    if (present(stress)) then
      stress_tensor = stress_tensor - stress%stress(sampling_points(i))
    endif
    
    b((i-1)*9+1:i*9) = reshape(dble(stress_tensor), [9])
  enddo
  
  ! Run linear least squares to get the basis function coefficients.
  ! This finds x s.t. (a.x-b)^2 is minimised.
  coefficients = dble(linear_least_squares(a, b))
  
  this%basis_functions_ = coefficients * this%basis_functions_
end procedure

module procedure interpolate_coefficients_CouplingStressBasisFunctions
  output = sum(this%basis_functions_%interpolate_coefficients( monomial,    &
                                                             & interpolator ))
end procedure

module procedure calculate_dynamical_matrices_CouplingStressBasisFunctions
  integer, allocatable :: subspaces_in_coupling(:)
  
  integer :: i
  
  subspaces_in_coupling = filter(subspaces%id .in. this%coupling%ids())
  
  output = [( StressDynamicalMatrix(anharmonic_data%structure%no_atoms), &
            & i=1,                                                       &
            & size(qpoints)                                              )]
  do i=1,size(this%basis_functions_)
    output = output                                                 &
         & + this%basis_functions_(i)%calculate_dynamical_matrices( &
         &                                   qpoints,               &
         &                                   thermal_energy,        &
         &                                   subspaces,             &
         &                                   subspace_bases,        &
         &                                   subspace_states,       &
         &                                   subspaces_in_coupling, &
         &                                   anharmonic_data        )
  enddo
end procedure

module procedure stress_correction_CouplingStressBasisFunctions
  integer :: i
  
  output = dblemat(zeroes(3,3))
  do i=1,size(this%basis_functions_)
    output = output + this%basis_functions_(i)%stress_correction( &
                                               & subspaces,       &
                                               & subspace_bases,  &
                                               & subspace_states, &
                                               & anharmonic_data  )
  enddo
end procedure

module procedure read_CouplingStressBasisFunctions
  type(String), allocatable :: line(:)
  
  type(SubspaceCoupling)                 :: coupling
  type(StressBasisFunction), allocatable :: basis_functions(:)
  
  select type(this); type is(CouplingStressBasisFunctions)
    line = split_line(input(1))
    coupling = SubspaceCoupling(join(line(3:)))
    
    basis_functions = StressBasisFunction(split_into_sections( &
                              & input(3:),                     &
                              & separating_line=repeat('-',50) ))
    
    this = CouplingStressBasisFunctions( coupling        = coupling,       &
                                       & basis_functions = basis_functions )
  class default
    call err()
  end select
end procedure

module procedure write_CouplingStressBasisFunctions
  select type(this); type is(CouplingStressBasisFunctions)
    output = [ 'Subspace Coupling: '//this%coupling,                      &
             & str(repeat('-',50)),                                       &
             & str(this%basis_functions_, separating_line=repeat('-',50)) ]
  class default
    call err()
  end select
end procedure

module procedure new_CouplingStressBasisFunctions_Strings
  call this%read(input)
end procedure

module procedure new_CouplingStressBasisFunctions_StringArray
  this = CouplingStressBasisFunctions(str(input))
end procedure
end submodule
