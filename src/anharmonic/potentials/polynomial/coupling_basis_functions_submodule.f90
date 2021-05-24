submodule (caesar_coupling_basis_functions_module) caesar_coupling_basis_functions_submodule
  use caesar_polynomial_module
contains

module procedure new_CouplingBasisFunctions_empty
  integer :: ialloc
  
  this%coupling = coupling
  allocate(this%basis_functions_(0), stat=ialloc); call err(ialloc)
end procedure

module procedure new_CouplingBasisFunctions
  this%coupling         = coupling
  this%basis_functions_ = basis_functions
end procedure

module procedure representation_CouplingBasisFunctions
  output = 'Polynomial basis functions'
end procedure

module procedure size_CouplingBasisFunctions
  output = size(this%basis_functions_)
end procedure

module procedure basis_functions_CouplingBasisFunctions
  output = this%basis_functions_
end procedure

module procedure optimise_CouplingBasisFunctions
  if (size(this%coupling%ids())/=1) then
    call print_line(CODE_ERROR//': Calling optimise_subspace_potential &
       &on a potential with coupled subspaces.')
    call err()
  elseif (any(this%coupling%ids()/=subspace%id)) then
    call print_line(CODE_ERROR//': Calling optimise_subspace_potential &
       &with the wrong subspace.')
    call err()
  endif
  
  ! Remove constant terms and split basis functions by power.
  this%basis_functions_ = optimise( this%basis_functions_,  &
                                  & subspace,               &
                                  & subspace_basis,         &
                                  & old_subspace_potential, &
                                  & anharmonic_data         )
end procedure

module procedure energy_RealModeDisplacement_CouplingBasisFunctions
  output = sum(this%basis_functions_%energy(displacement))
end procedure

module procedure energy_ComplexModeDisplacement_CouplingBasisFunctions
  output = sum(this%basis_functions_%energy(displacement))
end procedure

module procedure force_RealModeDisplacement_CouplingBasisFunctions
  if (size(this%basis_functions_)>0) then
    output = sum(this%basis_functions_%force(displacement))
  else
    output = RealModeForce(RealSingleForce(displacement%vectors%id,0.0_dp))
  endif
end procedure

module procedure force_ComplexModeDisplacement_CouplingBasisFunctions
  if (size(this%basis_functions_)>0) then
    output = sum(this%basis_functions_%force(displacement))
  else
    output = ComplexModeForce(ComplexSingleForce( displacement%vectors%id, &
                                                & (0.0_dp,0.0_dp)          ))
  endif
end procedure

module procedure braket_SubspaceBraKet_CouplingBasisFunctions
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

module procedure braket_BasisState_CouplingBasisFunctions
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

module procedure braket_BasisStates_CouplingBasisFunctions
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

module procedure harmonic_expectation_CouplingBasisFunctions
  integer :: i
  
  output = sum([( this%basis_functions_(i)%harmonic_expectation(    &
                &                                frequency,         &
                &                                thermal_energy,    &
                &                                supercell_size,    &
                &                                anharmonic_data ), &
                & i=1,                                              &
                & size(this%basis_functions_)                       )])
end procedure

module procedure potential_energy_SubspaceBraKet_CouplingBasisFunctions
  integer :: i
  
  output = sum([( this%basis_functions_(i)%potential_energy(    &
                &                            braket,            &
                &                            anharmonic_data ), &
                & i=1,                                          &
                & size(this%basis_functions_)                   )])
end procedure

module procedure potential_energy_BasisState_CouplingBasisFunctions
  integer :: i
  
  output = sum([( this%basis_functions_(i)%potential_energy(    &
                &                            bra,               &
                &                            ket,               &
                &                            subspace,          &
                &                            subspace_basis,    &
                &                            anharmonic_data ), &
                & i=1,                                          &
                & size(this%basis_functions_)                   )])
end procedure

module procedure generate_basis_functions_SubspaceCoupling
  type(SubspaceCombination), allocatable :: subspace_combinations(:)
  
  type(BasisFunctions), allocatable :: basis_functions(:)
  
  type(BasisFunction), allocatable :: coupling_basis_functions(:)
  
  integer :: i,ialloc
  
  ! Generate the set of subspace combinations corresponding to the subspace
  !    coupling.
  ! e.g. the coupling [1,2] might have combinations [1,2], [1,1,2] and [1,2,2].
  subspace_combinations = generate_subspace_combinations( &
              & coupling,                                 &
              & minimum_power = 2,                        &
              & maximum_power = potential_expansion_order )
    
  ! Loop over the subspace combinations corresponding to the coupling.
  allocate( basis_functions(size(subspace_combinations)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(subspace_combinations)
    ! Generate all basis functions for the subspace combination.
    basis_functions(i)%basis_functions = generate_basis_functions( &
                                      & subspace_combinations(i),  &
                                      & maximum_coupling_order,    &
                                      & structure,                 &
                                      & complex_modes,             &
                                      & qpoints,                   &
                                      & subspaces,                 &
                                      & degenerate_symmetries,     &
                                      & vscf_basis_functions_only, &
                                      & qpoint_symmetry_groups,    &
                                      & subspace_qpoint_stars,     &
                                      & logfile                    )
  enddo
  
  ! Concatenate the terms from each subspace combination together.
  coupling_basis_functions = [( basis_functions(i)%basis_functions, &
                              & i=1,                                &
                              & size(subspace_combinations)         )]
  
  output = CouplingBasisFunctions(coupling, coupling_basis_functions)
end procedure

module procedure undisplaced_energy_CouplingBasisFunctions
  type(RealSingleDisplacement) :: zero_displacement(0)
  
  output = this%energy(RealModeDisplacement(zero_displacement))
end procedure

module procedure zero_energy_CouplingBasisFunctions
  call print_line(CODE_ERROR//': zero_energy() cannot be called for type &
     &CouplingBasisFunctions.')
  call err()
end procedure

module procedure add_constant_CouplingBasisFunctions
  call print_line(CODE_ERROR//': add_constant() cannot be called for type &
     &CouplingBasisFunctions.')
  call err()
end procedure

module procedure fit_coefficients_CouplingBasisFunctions
  integer               :: dimensions
  real(dp), allocatable :: energy_differences(:)
  real(dp), allocatable :: sample_weights(:)
  real(dp), allocatable :: a(:,:)
  real(dp), allocatable :: b(:)
  
  real(dp), allocatable :: coefficients(:)
  
  if (size(this%basis_functions_)==0) then
    return
  endif
  
  ! Check inputs are consistent.
  if (size(sampling_points)/=size(sample_results)) then
    call print_line(CODE_ERROR//': The number of sampling points does not &
       &match the number of results.')
    call err()
  endif
  
  ! Each calculation yields size(modes) forces and one energy.
  dimensions = 1+size(modes)
  
  ! Calculate the weights to give to each sample.
  energy_differences = sample_results%energy-minval(sample_results%energy)
  sample_weights = min( 0.00003_dp/(0.00003_dp+energy_differences), &
                      & 0.01_dp)
  
  ! Calculate the energies and forces due to each basis function at each
  !    sampling point.
  a = construct_sample_matrix( this%basis_functions_, &
                             & sampling_points,       &
                             & modes,                 &
                             & energy_force_ratio,    &
                             & sample_weights         )
  
  ! Calculate the energies and forces sampled at each sampling point.
  b = construct_sample_vector( sampling_points,    &
                             & sample_results,     &
                             & potential,          &
                             & modes,              &
                             & energy_force_ratio, &
                             & sample_weights      )
  
  ! Run linear least squares to get the basis function coefficients.
  ! This finds x s.t. (a.x-b)^2 is minimised.
  coefficients = dble(linear_least_squares(a, b))
  
  this%basis_functions_ = coefficients * this%basis_functions_
end procedure

module procedure no_coefficients_CouplingBasisFunctions
  output = count(this%basis_functions_%power()<=maximum_power)
end procedure

module procedure coefficients_CouplingBasisFunctions
  output = this%basis_functions_(                           &
     & filter(this%basis_functions_%power()<=maximum_power) )%coefficient()
end procedure

module procedure set_coefficients_CouplingBasisFunctions
  integer :: i,j
  
  j = 0
  do i=1,size(this%basis_functions_)
    if (this%basis_functions_(i)%power()<=maximum_power) then
      j = j+1
      call this%basis_functions_(i)%set_coefficient(coefficients(j))
    endif
  enddo
  
  if (j/=size(coefficients)) then
    call print_line(CODE_ERROR//': Incorrect number of coefficients.')
    call err()
  endif
end procedure

module procedure zero_coefficients_CouplingBasisFunctions
  integer :: i
  
  do i=1,size(this%basis_functions_)
    call this%basis_functions_(i)%set_coefficient(0.0_dp)
  enddo
end procedure

module procedure all_basis_functions_CouplingBasisFunctions
  type(BasisFunction), allocatable :: basis_functions(:)
  
  basis_functions = this%basis_functions_
  call basis_functions%set_coefficient(1.0_dp)
  
  output = PotentialBasePointer(basis_functions)
end procedure

module procedure variable_basis_functions_CouplingBasisFunctions
  type(BasisFunction), allocatable :: basis_functions(:)
  
  basis_functions = this%basis_functions_(                  &
     & filter(this%basis_functions_%power()<=maximum_power) )
  call basis_functions%set_coefficient(1.0_dp)
  
  output = PotentialBasePointer(basis_functions)
end procedure

module procedure append_CouplingBasisFunctions
  type(BasisFunction), allocatable :: temp(:)
  
  if (this%coupling/=that%coupling) then
    call print_line(CODE_ERROR//': Appending incompatible basis functions.')
    call err()
  endif
  
  ! WORKAROUND: combining the following two lines makes ifort crash,
  !    as of version 19.0.4.227.
  temp = [this%basis_functions_, that%basis_functions_]
  this%basis_functions_ = temp
end procedure

module procedure interpolate_coefficient_CouplingBasisFunctions
  output = sum(this%basis_functions_%interpolate_coefficient( monomial,    &
                                                            & interpolator ))
end procedure

module procedure add_interpolated_contribution_CouplingBasisFunctions
  integer :: i,j
  
  do i=1,size(this%basis_functions_)
    do j=1,size(basis_function%basis_functions_)
      call this%basis_functions_(i)%add_interpolated_contribution( &
                             & basis_function%basis_functions_(j), &
                             & interpolator                        )
    enddo
  enddo
end procedure

module procedure add_harmonic_contribution_CouplingBasisFunctions
  integer :: i
  
  if (size(this%coupling%ids())==1) then
    do i=1,size(this%basis_functions_)
      call this%basis_functions_(i)%add_harmonic_contribution( &
                                         & dynamical_matrices, &
                                         & anharmonic_data     )
    enddo
  endif
end procedure

module procedure calculate_dynamical_matrices_CouplingBasisFunctions
  integer, allocatable :: subspaces_in_coupling(:)
  
  integer :: i
  
  subspaces_in_coupling = filter(subspaces%id .in. this%coupling%ids())
  
  output = [( DynamicalMatrix(anharmonic_data%structure%no_atoms), &
            & i=1,                                                 &
            & size(qpoints)                                        )]
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

module procedure energy_correction_CouplingBasisFunctions
  integer :: i
  
  output = 0
  do i=1,size(this%basis_functions_)
    output = output + this%basis_functions_(i)%energy_correction( &
                                               & subspaces,       &
                                               & subspace_bases,  &
                                               & subspace_states, &
                                               & anharmonic_data  )
  enddo
end procedure

module procedure read_CouplingBasisFunctions
  type(String), allocatable :: line(:)
  
  type(SubspaceCoupling)           :: coupling
  type(BasisFunction), allocatable :: basis_functions(:)
  
  select type(this); type is(CouplingBasisFunctions)
    line = split_line(input(1))
    
    if (lower_case(line(1))=='basis') then
      call print_line(ERROR//': Basis functions file appears to be in old &
         &format. Please run update_basis_functions.')
      call quit()
    endif
    
    coupling = SubspaceCoupling(join(line(3:)))
    
    basis_functions = BasisFunction(split_into_sections(input(3:)))
    
    this = CouplingBasisFunctions( coupling        = coupling,       &
                                 & basis_functions = basis_functions )
  class default
    call err()
  end select
end procedure

module procedure write_CouplingBasisFunctions
  select type(this); type is(CouplingBasisFunctions)
    output = [ 'Subspace Coupling: '//this%coupling,          &
             & str(''),                                       &
             & str(this%basis_functions_, separating_line='') ]
  class default
    call err()
  end select
end procedure

module procedure new_CouplingBasisFunctions_Strings
  call this%read(input)
end procedure

module procedure new_CouplingBasisFunctions_StringArray
  this = CouplingBasisFunctions(str(input))
end procedure
end submodule
