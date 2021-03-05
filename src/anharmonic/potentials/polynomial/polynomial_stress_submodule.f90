submodule (caesar_polynomial_stress_module) caesar_polynomial_stress_submodule
  use caesar_polynomial_module
contains

module procedure new_PolynomialStress
  this%stress_expansion_order_ = stress_expansion_order
  this%reference_stress_       = reference_stress
  this%basis_functions_        = basis_functions
end procedure

module procedure representation_PolynomialStress
  output = 'polynomial'
end procedure

module procedure zero_stress_PolynomialStress
  this%reference_stress_ = this%reference_stress_ - this%undisplaced_stress()
end procedure

module procedure add_constant_PolynomialStress
  this%reference_stress_ = this%reference_stress_ + input
end procedure

module procedure stress_RealModeDisplacement_PolynomialStress
  if (size(this%basis_functions_)>0) then
    output = this%reference_stress_ &
         & + sum(this%basis_functions_%stress(displacement))
  else
    output = this%reference_stress_
  endif
end procedure

module procedure stress_ComplexModeDisplacement_PolynomialStress
  if (size(this%basis_functions_)>0) then
    output = cmplxmat(this%reference_stress_) &
         & + sum(this%basis_functions_%stress(displacement))
  else
    output = cmplxmat(this%reference_stress_)
  endif
end procedure

module procedure braket_SubspaceBraKet_PolynomialStress
  integer :: i
  
  ! Integrate the reference stress (N.B. <i|e|j> = e<i|j> if e is a scalar.).
  this%reference_stress_ = this%reference_stress_ &
                       & * braket%inner_product(anharmonic_data)
  
  ! Integrate each basis function between the bra and the ket.
  do i=1,size(this%basis_functions_)
    call this%basis_functions_(i)%braket( braket,         &
                                        & whole_subspace, &
                                        & anharmonic_data )
  enddo
  
  ! Simplify the stress.
  if (set_default(whole_subspace,.true.)) then
    call this%simplify()
  endif
end procedure

module procedure braket_BasisState_PolynomialStress
  integer :: i
  
  ! Integrate the reference stress (N.B. <i|e|j> = e<i|j> if e is a scalar.).
  this%reference_stress_ = this%reference_stress_                        &
                       & * subspace_basis%inner_product( bra,            &
                       &                                 ket,            &
                       &                                 subspace,       &
                       &                                 anharmonic_data )
  
  ! Integrate each basis function between the bra and the ket.
  do i=1,size(this%basis_functions_)
    call this%basis_functions_(i)%braket( bra,            &
                                        & ket,            &
                                        & subspace,       &
                                        & subspace_basis, &
                                        & whole_subspace, &
                                        & anharmonic_data )
  enddo
  
  ! Simplify the stress.
  if (set_default(whole_subspace,.true.)) then
    call this%simplify()
  endif
end procedure

module procedure braket_BasisStates_PolynomialStress
  integer :: i
  
  ! Integrate each basis function between the bra and the ket.
  do i=1,size(this%basis_functions_)
    call this%basis_functions_(i)%braket( states,         &
                                        & subspace,       &
                                        & subspace_basis, &
                                        & whole_subspace, &
                                        & anharmonic_data )
  enddo
  
  ! Simplify the stress.
  if (set_default(whole_subspace,.true.)) then
    call this%simplify()
  endif
end procedure

module procedure simplify_PolynomialStress
  logical, allocatable :: to_remove(:)
  
  integer :: i,j
  
  to_remove = [(.false., i=1, size(this%basis_functions_))]
  
  do i=1,size(this%basis_functions_)
    if (size(this%basis_functions_(i)%coupling)==0) then
      this%reference_stress_ =      &
         &   this%reference_stress_ &
         & + this%basis_functions_(i)%undisplaced_stress()
      to_remove(i) = .true.
    else
      do j=1,i-1
        if (    this%basis_functions_(i)%coupling &
           & == this%basis_functions_(j)%coupling ) then
          if (to_remove(j)) then
            call err()
          endif
          call this%basis_functions_(j)%append( &
                     & this%basis_functions_(i) )
          to_remove(i) = .true.
        endif
      enddo
    endif
  enddo
  
  ! Remove constant and duplicate terms.
  this%basis_functions_ = this%basis_functions_(filter(.not.to_remove))
end procedure

module procedure harmonic_expectation_PolynomialStress
  integer :: i
  
  output = this%reference_stress_                                   &
       & + sum([( this%basis_functions_(i)%harmonic_expectation(    &
       &                                         frequency,         &
       &                                         thermal_energy,    &
       &                                         supercell_size,    &
       &                                         anharmonic_data ), &
       &          i=1,                                              &
       &          size(this%basis_functions_)                       )])
end procedure

module procedure can_be_interpolated_PolynomialStress
  output = .true.
end procedure

module procedure generate_stress_elements
  type(ComplexMonomial)   :: no_monomials(0)
  type(ComplexPolynomial) :: empty_polynomial
  type(ComplexPolynomial) :: elements(3,3)
  
  integer :: js(6) = [1,2,3,1,1,2]
  integer :: ks(6) = [1,2,3,2,3,3]
  
  integer :: i,j,k,ialloc
  
  empty_polynomial = ComplexPolynomial(terms=no_monomials)
  elements = empty_polynomial
  
  allocate(output(6), stat=ialloc); call err(ialloc)
  do i=1,6
    j = js(i)
    k = ks(i)
    elements(k,j) = ComplexPolynomial(monomials*coefficients%element(k,j))
    elements(j,k) = elements(k,j)
    output(i) = StressBasisFunction( elements    = elements, &
                                   & coefficient = 1.0_dp    )
    elements(k,j) = empty_polynomial
    elements(j,k) = empty_polynomial
  enddo
end procedure

module procedure interpolate_coefficients_PolynomialStress
  output = sum(this%basis_functions_%interpolate_coefficients( monomial,    &
                                                             & interpolator ))
end procedure

module procedure calculate_dynamical_matrices_PolynomialStress
  integer :: i
  
  output = [( StressDynamicalMatrix(anharmonic_data%structure%no_atoms), &
            & i=1,                                                       &
            & size(qpoints)                                              )]
  
  do i=1,size(this%basis_functions_)
    output = output                                                 &
         & + this%basis_functions_(i)%calculate_dynamical_matrices( &
         &                                         qpoints,         &
         &                                         thermal_energy,  &
         &                                         subspaces,       &
         &                                         subspace_bases,  &
         &                                         subspace_states, &
         &                                         anharmonic_data  )
  enddo
end procedure

module procedure stress_correction_PolynomialStress
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

module procedure expansion_order_PolynomialStress
  output = this%stress_expansion_order_
end procedure

module procedure read_PolynomialStress
  integer                                         :: expansion_order
  type(RealMatrix)                                :: reference_stress
  type(CouplingStressBasisFunctions), allocatable :: basis_functions(:)
  
  select type(this); type is(PolynomialStress)
    expansion_order = int(token(input(1),3))
    
    reference_stress = RealMatrix(input(3:5))
    
    basis_functions = CouplingStressBasisFunctions(split_into_sections( &
                                       & input(8:),                     &
                                       & separating_line=repeat('=',50) ))
    
    this = PolynomialStress( expansion_order,  &
                           & reference_stress, &
                           & basis_functions   )
  class default
    call err()
  end select
end procedure

module procedure write_PolynomialStress
  select type(this); type is(PolynomialStress)
    output = [ 'Expansion order: '//this%stress_expansion_order_,         &
             & str('Reference stress:'),                                  &
             & str(this%reference_stress_),                               &
             & str('Basis functions:'),                                   &
             & str(''),                                                   &
             & str(this%basis_functions_, separating_line=repeat('=',50)) ]
  class default
    call err()
  end select
end procedure

module procedure new_PolynomialStress_Strings
  call this%read(input)
end procedure

module procedure new_PolynomialStress_StringArray
  this = PolynomialStress(str(input))
end procedure
end submodule
