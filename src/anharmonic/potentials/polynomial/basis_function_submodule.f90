submodule (caesar_basis_function_module) caesar_basis_function_submodule
  use caesar_polynomial_module
contains

module procedure new_BasisFunction
  this%complex_representation_ = complex_representation
  if (present(coefficient)) then
    this%coefficient_ = coefficient
  else
    this%coefficient_ = 1.0_dp
  endif
end procedure

module procedure new_BasisFunction_PotentialBase
  select type(input); type is(PotentialBasePointer)
    this = new_BasisFunction_PotentialBase(input%potential())
  type is(PotentialPointer)
    this = new_BasisFunction_PotentialBase(input%potential())
  type is(BasisFunction)
    this = input
  class default
    call err()
  end select
end procedure

module procedure representation_BasisFunction
  output = 'Polynomial basis function'
end procedure

module procedure complex_representation
  output = this%coefficient_*this%complex_representation_
end procedure

module procedure generate_basis_functions_SubspaceCombination
  type(QpointStarProduct),     allocatable :: qpoint_star_products(:)
  type(CombinationQpointStar), allocatable :: combination_qpoint_stars(:)
  
  type(ComplexMonomial), allocatable :: complex_monomials(:)
  
  integer :: i,j
  
  call print_line('Generating basis functions in subspace combination '// &
     &subspace_combination)
  
  if (sum(subspace_combination%powers())<2) then
    call print_line(CODE_ERROR//': Trying to generate basis functions with &
       &power less than 2.')
    call err()
  endif
  
  output = [BasisFunction::]
  
  ! Generate the products of the single-subspace q-point stars.
  qpoint_star_products = generate_qpoint_star_products( &
                                & subspace_combination, &
                                & subspace_qpoint_stars )
  do i=1,size(qpoint_star_products)
    ! From the products of the single-subspace q-point stars,
    !    generate the multi-subspace q-point stars.
    combination_qpoint_stars = generate_combination_qpoint_stars( &
                                       & qpoint_star_products(i), &
                                       & qpoint_symmetry_groups,  &
                                       & .true.,                  &
                                       & qpoints                  )
    do j=1,size(combination_qpoint_stars)
      ! For each q-point star, generate the complex monomials associated with
      !    that star, and generate the symmetry-invariant basis functions
      !    from these monomials.
      complex_monomials = combination_qpoint_stars(j)%complex_monomials( &
                                                         & complex_modes )
      
      output = [ output,                                                &
               & monomials_to_basis_functions( complex_monomials,       &
               &                               structure,               &
               &                               complex_modes,           &
               &                               qpoints,                 &
               &                               degenerate_symmetries,   &
               &                               logfile                ) ]
    enddo
  enddo
  
  !! Generate the complex monomials corresponding to the subspace combination,
  !!    with coefficients such that symmetries are unitary.
  !complex_monomials = subspace_combination%complex_monomials( &
  !     & maximum_coupling_order,                              &
  !     & subspaces,                                           &
  !     & complex_modes,                                       &
  !     & qpoints,                                             &
  !     & conserve_momentum=.true.,                            &
  !     & conserve_subspace_momentum=vscf_basis_functions_only )
  !
  !output = monomials_to_basis_functions( complex_monomials,     &
  !                                     & structure,             &
  !                                     & complex_modes,         &
  !                                     & qpoints,               &
  !                                     & degenerate_symmetries, &
  !                                     & logfile                )
end procedure

! Takes an array of complex monomials, and generates the basis functions
!    containing these monomials.
module function monomials_to_basis_functions(complex_monomials,structure, &
   & complex_modes,qpoints,degenerate_symmetries,logfile) result(output)
  type(ComplexMonomial),    intent(in)              :: complex_monomials(:)
  type(StructureData),      intent(in)              :: structure
  type(ComplexMode),        intent(in)              :: complex_modes(:)
  type(QpointData),         intent(in)              :: qpoints(:)
  type(DegenerateSymmetry), intent(in)              :: degenerate_symmetries(:)
  type(OFile),              intent(inout), optional :: logfile
  type(BasisFunction), allocatable                  :: output(:)
  
  type(BasisConversion) :: basis_conversion
  
  ! Symmetry data.
  type(ComplexMatrix) :: symmetry
  integer             :: order
  type(ComplexMatrix) :: projection
  type(RealMatrix)    :: basis_projection
  
  ! Polynomial coefficients, in both bases.
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  
  ! Variables for constructing the output.
  integer, allocatable :: unique_terms(:)
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  ! Construct projection matrix, which has allowed basis functions as
  !    eigenvectors with eigenvalue 1, and sends all other functions to 0.
  projection = cmplxmat(make_identity_matrix(size(complex_monomials)))
  do i=1,size(degenerate_symmetries)
    order = structure%symmetries(                                             &
       & first(structure%symmetries%id==degenerate_symmetries(i)%symmetry_id) &
       & )%symmetry_order()
    
    ! Constuct symmetry in complex monomial co-ordinates.
    symmetry = degenerate_symmetries(i)%calculate_symmetry( &
                              & complex_monomials,          &
                              & complex_modes,              &
                              & include_coefficients=.true. )
    call check_unitary(symmetry, 'symmetry in monomial basis', logfile)
    
    ! Construct the projection matrix for this symmetry,
    !    and multiply the total projection matrix by this.
    projection = projection * projection_matrix(symmetry, order)
  enddo
  call check_hermitian( projection,                   &
                      & 'monomial projection matrix', &
                      & logfile,                      &
                      & ignore_threshold=1e-10_dp     )
  
  ! Find the indices of the conjugates, such that
  !    complex_monomials(conjugates(i)) == conjg(complex_monomials(i)).
  basis_conversion = BasisConversion(complex_monomials)
  
  ! Transform the projection matrix to basis functions rather than monomials.
  basis_projection = basis_conversion%matrix_to_basis( projection, &
                                                     & lhs=.true., &
                                                     & rhs=.true.  )
  call check_symmetric( basis_projection,          &
                      & 'basis projection matrix', &
                      & logfile,                   &
                      & ignore_threshold=1e-10_dp  )
  
  ! Diagonalise the projection matrix,
  !    check its eigenvalues are either 0 or 1,
  !    and select only the eigenvectors with eigenvalue 1.
  estuff = diagonalise_symmetric(basis_projection)
  if (any(abs(estuff%eval-1)>1e-2_dp .and. abs(estuff%eval)>1e-2_dp)) then
    call print_line(ERROR//': Projection matrix has eigenvalues which are &
       &neither 0 nor 1.')
    call print_line('Eigenvalues:')
    call print_line(estuff%eval)
    call err()
  endif
  estuff = estuff(filter(abs(estuff%eval-1)<1e-2_dp))
  
  ! Take linear combinations of basis functions such that each basis function
  !    contains at least one term which is in no other basis function.
  allocate( unique_terms(size(estuff)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(estuff)
    ! Identify the largest term in basis function i.
    unique_terms(i) = maxloc(abs(estuff(i)%evec), 1)
    
    ! Subtract a multiple of basis function i from all other basis functions,
    !    such that the coefficient of unique_term_id(i) in all other basis
    !    functions is zero.
    do j=1,size(estuff)
      if (j/=i) then
        estuff(j)%evec = estuff(j)%evec                  &
                     & - estuff(i)%evec                  &
                     & * estuff(j)%evec(unique_terms(i)) &
                     & / estuff(i)%evec(unique_terms(i))
      endif
    enddo
  enddo
  
  allocate(output(size(estuff)), stat=ialloc); call err(ialloc)
  do i=1,size(estuff)
    ! Construct a basis function from the eigenvectors.
    output(i) = BasisFunction(                                       &
       & basis_conversion%vector_from_basis(estuff(i)%evec, 1e-4_dp) )
  enddo
end function

module procedure projection_matrix
  type(ComplexMatrix) :: identity
  
  integer :: i
  
  if (order<1) then
    call print_line(CODE_ERROR//': symmetry order may not be < 1.')
    call err()
  elseif (order>6) then
    call print_line(CODE_ERROR//': symmetry order may not be > 6.')
    call err()
  endif
  
  identity = cmplxmat(make_identity_matrix(size(input,1)))
  
  output = identity
  do i=2,order
    output = input*output + identity
  enddo
  output = output/order
end procedure

module procedure simplify_BasisFunction
  call this%complex_representation_%simplify()
end procedure

module procedure optimise_BasisFunctions
  type(ComplexMode), allocatable :: subspace_modes(:)
  
  integer :: order
  
  type(ComplexPolynomial), allocatable :: complex_polynomials(:)
  
  type(ComplexMonomial) :: no_complex_monomials(0)
  
  type(PotentialBasePointer), allocatable :: basis_functions_(:)
  type(BasisFunction),        allocatable :: basis_functions(:)
  type(ComplexPolynomial),    allocatable :: basis_polynomials(:)
  
  type(SubspaceQpointStars), allocatable :: subspace_qpoint_stars(:)
  
  type(ComplexPolynomial), allocatable :: order_polynomials(:)
  
  integer :: i,j,ialloc
  
  subspace_modes = anharmonic_data%complex_modes(                     &
     & filter(anharmonic_data%complex_modes%subspace_id==subspace%id) )
  
  order = anharmonic_data%potential_expansion_order
  
  ! Collate the input into an array of polynomials, one for each power.
  complex_polynomials = [( ComplexPolynomial(no_complex_monomials), &
                         & i=1,                                     &
                         & order                                    )]
  do i=1,size(input)
    do j=1,size(input(i)%complex_representation_%terms)
      associate (term => input(i)%complex_representation_%terms(j))
        if (term%total_power()==0) then
          ! The terms with power=0 are dropped as they have already been
          !    accounted for in the reference energy.
          cycle
        elseif (.not. is_int(term%wavevector( &
                    & subspace_modes,         &
                    & anharmonic_data%qpoints ))) then
          ! The terms with (sum q/=0) are dropped as they are not invariant
          !    under translational symmetry.
          cycle
        endif
        complex_polynomials(term%total_power()) =      &
           &   complex_polynomials(term%total_power()) &
           & + term
      end associate
    enddo
  enddo
  
  if (present(old_subspace_potential)) then
    basis_functions_ = old_subspace_potential%all_basis_functions( &
                                                 & anharmonic_data )
    basis_functions = [( BasisFunction(basis_functions_(i)), &
                       & i=1,                                &
                       & size(basis_functions_)              )]
    
    basis_polynomials = basis_functions%complex_representation_
  endif
  
  if (.not. allocated(basis_polynomials)) then
    subspace_qpoint_stars = generate_subspace_qpoint_stars(                  &
       & subspaces              = [subspace],                                &
       & modes                  = anharmonic_data%complex_modes,             &
       & qpoints                = anharmonic_data%qpoints,                   &
       & qpoint_symmetry_groups = anharmonic_data%qpoint_symmetry_groups,    &
       & max_power              = anharmonic_data%potential_expansion_order, &
       & conserve_momentum      = .true.                                     )
  endif
  
  allocate(output(0), stat=ialloc); call err(ialloc)
  do i=1,order
    if (allocated(basis_polynomials)) then
      order_polynomials = basis_polynomials(                             &
         & filter([( basis_polynomials(j)%terms(1)%total_power()==i,     &
         &           j=1,                                                &
         &           size(basis_polynomials)                         )]) )
    else
      order_polynomials = construct_basis_polynomials( &
                           & subspace,                 &
                           & i,                        &
                           & subspace_basis,           &
                           & subspace_qpoint_stars(1), &
                           & subspace_modes,           &
                           & anharmonic_data           )
    endif
    
    if (size(order_polynomials)>0) then
      output = [ output,                                              &
               & fit_basis_functions( complex_polynomials(i)%terms,   &
               &                      order_polynomials             ) ]
    endif
  enddo
end procedure

! Construct the basis functions at a given order.
module function construct_basis_polynomials(subspace,order,               &
   & subspace_basis,subspace_qpoint_stars,subspace_modes,anharmonic_data) &
   & result(output) 
  type(DegenerateSubspace),  intent(in) :: subspace
  integer,                   intent(in) :: order
  class(SubspaceBasis),      intent(in) :: subspace_basis
  type(SubspaceQpointStars), intent(in) :: subspace_qpoint_stars
  type(ComplexMode),         intent(in) :: subspace_modes(:)
  type(AnharmonicData),      intent(in) :: anharmonic_data
  type(ComplexPolynomial), allocatable  :: output(:)
  
  type(DegenerateSymmetry), allocatable :: symmetries(:)
  
  type(QpointStar), allocatable :: qpoint_stars(:)
  
  type(ComplexMonomial), allocatable :: complex_monomials(:)
  
  type(BasisFunction), allocatable :: basis_functions(:)
  
  real(dp) :: energy_scale
  
  integer :: i
  
  symmetries = subspace_basis%select_symmetries( &
        & anharmonic_data%degenerate_symmetries, &
        & anharmonic_data                        )
  
  output = [ComplexPolynomial::]
  
  qpoint_stars = subspace_qpoint_stars%powers(order+1)%stars
  do i=1,size(qpoint_stars)
    complex_monomials = qpoint_stars(i)%complex_monomials( &
                                          & subspace_modes )
    
    basis_functions = monomials_to_basis_functions( &
                       & complex_monomials,         &
                       & anharmonic_data%structure, &
                       & subspace_modes,            &
                       & anharmonic_data%qpoints,   &
                       & symmetries                 )
    
    output = [ output,                                         &
             & ( basis_functions(i)%complex_representation_,   &
             &   i=1,                                          &
             &   size(basis_functions)                       ) ]
  enddo
  
  !subspace_combination = SubspaceCombination( ids    = [subspace%id], &
  !                                          & powers = [order]        )
  !complex_monomials = subspace_combination%complex_monomials( &
  !                & anharmonic_data%maximum_coupling_order,   &
  !                & [subspace],                               &
  !                & anharmonic_data%complex_modes,            &
  !                & anharmonic_data%qpoints,                  &
  !                & conserve_momentum=.true.,                 &
  !                & conserve_subspace_momentum=.true.         )
  !
  !basis_functions = monomials_to_basis_functions( &
  !               & complex_monomials,             &
  !               & anharmonic_data%structure,     &
  !               & anharmonic_data%complex_modes, &
  !               & anharmonic_data%qpoints,       &
  !               & symmetries                     )
  !
  !output = [( basis_functions(i)%complex_representation_, &
  !          & i=1,                                        &
  !          & size(basis_functions)                       )]
  
  ! A term |u|^(2n) scales like (2Nw)^{-n}.
  ! w must be capped so that it is not considered to be too small, otherwise
  !    some terms become excessively weighted.
  energy_scale = ( 2                                              &
            &    * anharmonic_data%anharmonic_supercell%sc_size   &
            &    * max(subspace_basis%frequency,1e-4_dp)        ) &
            & ** (0.5_dp*order)
  
  ! Multiply the basis functions by energy_scale, E,
  !    so that the basis functions are dimensionless,
  !    and the coefficients are real and have dimensions of energy.
  do i=1,size(output)
    output(i)%terms%coefficient = &
       & output(i)%terms%coefficient * energy_scale
  enddo
end function

module procedure fit_basis_functions
  logical, allocatable :: basis_functions_present(:)
  
  complex(dp), allocatable :: complex_conversion(:,:)
  type(RealMatrix)         :: conversion
  type(RealMatrix)         :: inverse_conversion
  
  type(BasisConversion) :: basis_conversion
  
  integer :: i,j,k,ialloc
  
  ! Find the indices of the conjugates, such that
  !    monomials(conjugates(i)) == conjg(monomials(i)).
  basis_conversion = BasisConversion(monomials)
  
  basis_functions_present = [(.false., i=1, size(basis_functions))]
  
  ! Construct the mapping from the basis functions to the monomials.
  allocate( complex_conversion(size(monomials),size(basis_functions)), &
          & stat=ialloc); call err(ialloc)
  complex_conversion = cmplx(0.0_dp,0.0_dp,dp)
  do i=1,size(basis_functions)
    do j=1,size(basis_functions(i)%terms)
      k = first_equivalent( monomials,                   &
                          & basis_functions(i)%terms(j), &
                          & compare_complex_monomials,   &
                          & default=0                    )
      if (k/=0) then
        basis_functions_present(i) = .true.
        complex_conversion(k,i) = basis_functions(i)%terms(j)%coefficient
      endif
    enddo
  enddo
  
  if (.not. any(basis_functions_present)) then
    allocate(output(0), stat=ialloc); call err(ialloc)
    return
  endif
  
  complex_conversion = complex_conversion(:,filter(basis_functions_present))
  conversion = basis_conversion%matrix_to_basis( mat(complex_conversion), &
                                               & lhs=.true.,              &
                                               & rhs=.false.              )
  
  ! Invert the conversion from basis functions to paired monomials,
  !    to give the conversion from paired monomials to basis functions.
  ! Uses X^-1 = (X^T.X)^-1 . X^T.
  inverse_conversion = invert(transpose(conversion)*conversion) &
                   & * transpose(conversion)
  
  output = [(BasisFunction(basis_functions(i)), i=1, size(basis_functions))]
  output = output(filter(basis_functions_present))
  output%coefficient_ = dble( inverse_conversion                       &
                          & * vec(basis_conversion%vector_to_basis(    &
                          &                   monomials%coefficient )) )
end procedure

module procedure energy_RealModeDisplacement_BasisFunction
  output = real( this%coefficient_                                 &
             & * this%complex_representation_%energy(displacement) )
end procedure

module procedure energy_ComplexModeDisplacement_BasisFunction
  output = this%coefficient_ &
       & * this%complex_representation_%energy(displacement)
end procedure

module procedure force_RealModeDisplacement_BasisFunction
  output = this%coefficient_ &
       & * this%complex_representation_%force(displacement)
end procedure

module procedure force_ComplexModeDisplacement_BasisFunction
  output = this%coefficient_ &
       & * this%complex_representation_%force(displacement)
end procedure

module procedure braket_SubspaceBraKet_BasisFunction
  ! Perform integration in complex co-ordinates.
  call integrate( this%complex_representation_%terms, &
                & braket,                             &
                & anharmonic_data                     )
end procedure

module procedure braket_BasisState_BasisFunction
  ! Perform integration in complex co-ordinates.
  call integrate( this%complex_representation_%terms, &
                & bra,                                &
                & ket,                                &
                & subspace,                           &
                & subspace_basis,                     &
                & anharmonic_data                     )
end procedure

module procedure braket_BasisStates_BasisFunction
  integer :: i
  
  ! Perform integration in complex co-ordinates.
  do i=1,size(this%complex_representation_%terms) 
    call integrate( this%complex_representation_%terms(i), &
                  & states,                                &
                  & subspace,                              &
                  & subspace_basis,                        &
                  & anharmonic_data                        )
  enddo
end procedure

module procedure harmonic_expectation_BasisFunction
  output = this%coefficient_                                                  &
       & * this%complex_representation_%harmonic_expectation( frequency,      &
       &                                                      thermal_energy, &
       &                                                      supercell_size  )
end procedure

module procedure potential_energy_SubspaceBraKet
  output = this%coefficient_                                         &
       & * real(integrate_to_constant( this%complex_representation_, &
       &                               braket,                       &
       &                               anharmonic_data               ))
end procedure

module procedure potential_energy_BasisState
  output = this%coefficient_                                         &
       & * real(integrate_to_constant( this%complex_representation_, &
       &                               bra,                          &
       &                               ket,                          &
       &                               subspace,                     &
       &                               subspace_basis,               &
       &                               anharmonic_data               ))
end procedure

module procedure undisplaced_energy_BasisFunction
  type(RealSingleDisplacement) :: zero_displacement(0)
  
  output = this%energy(RealModeDisplacement(zero_displacement))
end procedure

module procedure zero_energy_BasisFunction
  call print_line(CODE_ERROR//': zero_energy() cannot be called for type &
     &BasisFunction.')
  call err()
end procedure

module procedure add_constant_BasisFunction
  call print_line(CODE_ERROR//': add_constant() cannot be called for type &
     &BasisFunction.')
  call err()
end procedure

module procedure power_BasisFunction
  output = this%complex_representation_%terms(1)%total_power()
end procedure

module procedure coefficient_BasisFunction
  output = this%coefficient_
end procedure

module procedure set_coefficient_BasisFunction
  this%coefficient_ = coefficient
end procedure

module procedure multiply_BasisFunction_real
  output = this
  output%coefficient_ = output%coefficient_ * that
end procedure

module procedure multiply_real_BasisFunction
  output = that
  output%coefficient_ = this * output%coefficient_
end procedure

module procedure divide_BasisFunction_real
  output = this
  output%coefficient_ = output%coefficient_ / that
end procedure

module procedure add_BasisFunction_BasisFunction
  output = BasisFunction(                                          &
     & this%complex_representation()+that%complex_representation() )
end procedure

module procedure negative_BasisFunction
  output = this
  output%coefficient_ = -output%coefficient_
end procedure

module procedure subtract_BasisFunction_BasisFunction
  output = BasisFunction(                                          &
     & this%complex_representation()-that%complex_representation() )
end procedure

module procedure interpolate_coefficient_BasisFunction
  output = interpolator%overlap(monomial, this%complex_representation_) &
       & * this%coefficient_
end procedure

module procedure add_interpolated_contribution_BasisFunction
  type(ComplexMonomial) :: monomial
  
  complex(dp) :: coefficient
  
  ! Locate the monomial in the basis function with the largest coefficient
  monomial = this%complex_representation_%terms(                      &
     & maxloc(abs(this%complex_representation_%terms%coefficient), 1) )
  
  ! Calculate the coefficient of this monomial when the input basis function
  !    is interpolated.
  coefficient = basis_function%interpolate_coefficient(monomial, interpolator)
  
  ! Divide by the monomial's coefficient to get the contribution to this basis
  !    function from the input basis function.
  this%coefficient_ = this%coefficient_ &
                  & + real(coefficient/monomial%coefficient)
end procedure

module procedure add_harmonic_contribution_BasisFunction
  type(ComplexMode) :: mode
  
  real(dp) :: frequency
  
  integer :: i
  
  if (size(this%complex_representation_%terms)==0) then
    return
  endif
  
  associate(monomial=>this%complex_representation_%terms(1))
    if (monomial%total_power()==2) then
      mode = anharmonic_data%complex_modes(                        &
         & first(anharmonic_data%complex_modes%id==monomial%id(1)) )
      i = first(anharmonic_data%qpoints%id==mode%qpoint_id)
      
      frequency = dynamical_matrices(i)%expectation(mode)
      
      this%coefficient_ = this%coefficient_                    &
                      & + real(frequency/monomial%coefficient) &
                      & * anharmonic_data%anharmonic_supercell%sc_size
    endif
  end associate
end procedure

module procedure calculate_dynamical_matrices_BasisFunction
  output = calculate_dynamical_matrices( this%complex_representation_,   &
       &                                 qpoints,                        &
       &                                 thermal_energy,                 &
       &                                 subspaces,                      &
       &                                 subspace_bases,                 &
       &                                 subspace_states,                &
       &                                 subspaces_in_coupling,          &
       &                                 anharmonic_data               ) &
       & * this%coefficient_
end procedure

module procedure energy_correction_BasisFunction
  output = calculate_correction( this%complex_representation_,   &
       &                         subspaces,                      &
       &                         subspace_bases,                 &
       &                         subspace_states,                &
       &                         anharmonic_data               ) &
       & * this%coefficient_
end procedure

module procedure terms_BasisFunction
  output = this%complex_representation_%terms
end procedure

module procedure read_BasisFunction
  integer :: partition_line
  
  type(ComplexPolynomial) :: complex_representation
  
  integer :: i
  
  select type(this); type is(BasisFunction)
    
    ! Locate the line between real terms and complex terms.
    ! This is included for legacy reasons.
    partition_line = 1
    do i=2,size(input)
      if (size(split_line(input(i)))>1) then
        partition_line = i
        exit
      endif
    enddo
    
    complex_representation = ComplexPolynomial(          &
       & join(input(partition_line+1:), delimiter=' + ') )
    
    this = BasisFunction(complex_representation = complex_representation)
  class default
    call err()
  end select
end procedure

module procedure write_BasisFunction
  type(ComplexPolynomial) :: complex_representation
  
  select type(this); type is(BasisFunction)
    complex_representation = this%complex_representation()
    output = [ str('Basis function in complex co-ordinates:'), &
             & str(complex_representation%terms)               ]
  class default
    call err()
  end select
end procedure

module procedure new_BasisFunction_Strings
  call this%read(input)
end procedure

module procedure new_BasisFunction_StringArray
  this = BasisFunction(str(input))
end procedure
end submodule
