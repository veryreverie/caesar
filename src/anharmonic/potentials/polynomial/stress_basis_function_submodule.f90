submodule (caesar_stress_basis_function_module) caesar_stress_basis_function_submodule
  use caesar_polynomial_module
contains

module procedure new_StressBasisFunction
  this%elements_ = elements
  if (present(coefficient)) then
    this%coefficient_ = coefficient
  else
    this%coefficient_ = 1
  endif
end procedure

module procedure representation_StressBasisFunction
  output = 'Polynomial stress basis function'
end procedure

module procedure generate_stress_basis_functions_SubspaceMonomial
  ! Complex monomials, and the indices of their conjugates, such that
  !    complex_monomials(conjugates(i)) == conjg(complex_monomials(i)).
  type(ComplexMonomial), allocatable :: complex_monomials(:)
  
  ! Each symmetry, acting on scalar basis functions, the stress tensor,
  !    and the combined basis function.
  type(ComplexMatrix)   :: complex_scalar_symmetry
  real(dp), allocatable :: real_scalar_symmetry(:,:)
  real(dp)              :: tensor_symmetry(6,6)
  real(dp), allocatable :: symmetry(:,:)
  type(RealMatrix)      :: projection
  
  ! Polynomial coefficients, in both bases.
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  
  ! Variables for constructing the output.
  type(ComplexPolynomial) :: elements(3,3)
  type(ComplexPolynomial) :: basis_function
  
  ! Mappings between 3x3 indices and 6 indices.
  integer :: x(6)
  integer :: y(6)
  
  ! Temporary variables.
  integer  :: i,j,k,l,m,n,o,ialloc
  real(dp) :: tensor(3,3)
  
  type(BasisConversion) :: basis_conversion
  
  ! Initialise index mappings.
  x = [1,2,3,1,1,2]
  y = [1,2,3,2,3,3]
  
  ! Generate the complex monomials corresponding to the subspace monomial,
  !    with coefficients such that symmetries are unitary.
  complex_monomials = generate_complex_monomials(            &
      & subspace_monomial,                                   &
      & maximum_coupling_order,                              &
      & subspaces,                                           &
      & complex_modes,                                       &
      & qpoints,                                             &
      & conserve_momentum=.true.,                            &
      & conserve_subspace_momentum=vscf_basis_functions_only )
  
  if (size(complex_monomials)==0) then
    allocate(output(0), stat=ialloc); call err(ialloc)
    return
  endif
  
  ! Locate the conjugate of each complex monomial.
  basis_conversion = BasisConversion(complex_monomials)
  
  ! Construct projection matrix, which has allowed basis functions as
  !    eigenvectors with eigenvalue 1, and sends all other functions to 0.
  projection = dblemat(make_identity_matrix(size(complex_monomials)*6))
  allocate( symmetry(size(complex_monomials)*6,size(complex_monomials)*6), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(degenerate_symmetries)
    ! Construct the symmetry acting on the complex monomials.
    complex_scalar_symmetry = degenerate_symmetries(i)%calculate_symmetry( &
                                             & complex_monomials,          &
                                             & complex_modes,              &
                                             & include_coefficients=.true. )
    call check_unitary( complex_scalar_symmetry,              &
                      & 'symmetry in complex monomial basis', &
                      & logfile                               )
    
    ! Transform the symmetry into real co-ordinates.
    real_scalar_symmetry = dble(basis_conversion%matrix_to_basis( &
                                       & complex_scalar_symmetry, &
                                       & lhs=.true.,              &
                                       & rhs=.true.               ))
    
    ! Construct the symmetry acting on the tensor components.
    do j=1,6
      tensor = 0.0_dp
      if (x(j)==y(j)) then
        tensor(x(j),y(j)) = 1.0_dp
      else
        tensor(x(j),y(j)) = 1/sqrt(2.0_dp)
        tensor(y(j),x(j)) = 1/sqrt(2.0_dp)
      endif
      tensor = dble( structure%symmetries(i)%cartesian_tensor         &
                 & * mat(tensor)                                      &
                 & * invert(structure%symmetries(i)%cartesian_tensor) )
      do k=1,6
        if (x(k)==y(k)) then
          tensor_symmetry(k,j) = tensor(x(k),y(k))
        else
          tensor_symmetry(k,j) = (tensor(x(k),y(k)) + tensor(y(k),x(k))) &
                             & / sqrt(2.0_dp)
        endif
      enddo
    enddo
    call check_orthogonal( mat(tensor_symmetry),                     &
                         & 'symmetry in basis of tensor components', &
                         & logfile                                   )
    
    ! Construct the full symmetry, transforming both scalar and tensor
    !    compontents.
    l = 0
    do j=1,6
      do k=1,size(complex_monomials)
        l = l+1
        
        o = 0
        do m=1,6
          do n=1,size(complex_monomials)
            o = o+1
            
            symmetry(o,l) = tensor_symmetry(m,j) * real_scalar_symmetry(n,k)
          enddo
        enddo
      enddo
    enddo
    
    ! Construct the projection matrix for this symmetry,
    !    and multiply the total projection matrix by this.
    projection = projection                                                  &
             & * projection_matrix( mat(symmetry),                           &
             &                      structure%symmetries(i)%symmetry_order() )
  enddo
  call check_symmetric( projection,               &
                      & 'projection matrix',      &
                      & logfile,                  &
                      & ignore_threshold=1e-10_dp )
  
  ! Diagonalise the projection matrix,
  !    check its eigenvalues are either 0 or 1,
  !    and select only the eigenvectors with eigenvalue 1.
  estuff = diagonalise_symmetric(projection)
  if (any(abs(estuff%eval-1)>1e-2_dp .and. abs(estuff%eval)>1e-2_dp)) then
    call print_line(ERROR//': Projection matrix has eigenvalues which are &
       &neither 0 nor 1.')
    call print_line('Eigenvalues:')
    call print_line(estuff%eval)
    call err()
  endif
  estuff = estuff(filter(abs(estuff%eval-1)<1e-2_dp))
  
  ! Construct basis functions from coefficients.
  allocate(output(size(estuff)), stat=ialloc); call err(ialloc)
  do i=1,size(estuff)
    do j=1,6
      basis_function = basis_conversion%vector_from_basis(    &
         & estuff(i)%evec( size(complex_monomials)*(j-1)+1    &
         &               : size(complex_monomials)*j       ), &
         & 1e-4_dp                                            )
      call basis_function%simplify()
      
      if (x(j)==y(j)) then
        elements(x(j),y(j)) = basis_function
      else
        elements(x(j),y(j)) = basis_function / sqrt(2.0_dp)
        elements(y(j),x(j)) = elements(x(j),y(j))
      endif
    enddo
    
    output(i) = StressBasisFunction(elements)
  enddo
end procedure

module procedure projection_matrix
  type(RealMatrix) :: identity
  
  integer :: i
  
  if (order<1) then
    call print_line(CODE_ERROR//': symmetry order may not be < 1.')
    call err()
  elseif (order>6) then
    call print_line(CODE_ERROR//': symmetry order may not be > 6.')
    call err()
  endif
  
  identity = dblemat(make_identity_matrix(size(input,1)))
  
  output = identity
  do i=2,order
    output = input*output + identity
  enddo
  output = output/order
end procedure

module procedure simplify_StressBasisFunction
  call this%elements_%simplify()
end procedure

module procedure stress_RealModeDisplacement_StressBasisFunction
  real(dp) :: elements(3,3)
  
  integer :: i,j
  
  do i=1,3
    do j=1,3
      elements(j,i) = real(this%elements_(j,i)%energy(displacement))
    enddo
  enddo
  
  output = this%coefficient_ * mat(elements)
end procedure

module procedure stress_ComplexModeDisplacement_StressBasisFunction
  complex(dp) :: elements(3,3)
  
  integer :: i,j
  
  do i=1,3
    do j=1,3
      elements(j,i) = this%elements_(j,i)%energy(displacement)
    enddo
  enddo
  
  output = this%coefficient_ * mat(elements)
end procedure

module procedure braket_SubspaceBraKet_StressBasisFunction
  integer :: i,j
  
  do i=1,3
    do j=1,3
      call integrate( this%elements_(j,i)%terms, &
                    & braket,                    &
                    & anharmonic_data            )
    enddo
  enddo
end procedure

module procedure braket_BasisState_StressBasisFunction
  integer :: i,j
  
  do i=1,3
    do j=1,3
      call integrate( this%elements_(j,i)%terms, &
                    & bra,                       &
                    & ket,                       &
                    & subspace,                  &
                    & subspace_basis,            &
                    & anharmonic_data            )
    enddo
  enddo
end procedure

module procedure braket_BasisStates_StressBasisFunction
  integer :: i,j,k
  
  do i=1,3
    do j=1,3
      do k=1,size(this%elements_(j,i)%terms)
        call integrate( this%elements_(j,i)%terms(k), &
                      & states,                       &
                      & subspace,                     &
                      & subspace_basis,               &
                      & anharmonic_data               )
      enddo
    enddo
  enddo
end procedure

module procedure harmonic_expectation_StressBasisFunction
  real(dp) :: elements(3,3)
  
  integer :: i,j
  
  do i=1,3
    do j=1,3
      elements(j,i) = this%elements_(j,i)%harmonic_expectation( &
                                              & frequency,      &
                                              & thermal_energy, &
                                              & supercell_size  )
    enddo
  enddo
  
  output = this%coefficient_ * mat(elements)
end procedure

module procedure undisplaced_stress_StressBasisFunction
  type(RealSingleDisplacement) :: zero_displacement(0)
  
  output = this%stress(RealModeDisplacement(zero_displacement))
end procedure

module procedure zero_stress_StressBasisFunction
  call print_line(CODE_ERROR//': zero_energy() cannot be called for type &
     &StressBasisFunction.')
  call err()
end procedure

module procedure add_constant_StressBasisFunction
  call print_line(CODE_ERROR//': add_constant() cannot be called for type &
     &StressBasisFunction.')
  call err()
end procedure

module procedure multiply_StressBasisFunction_real
  output = this
  output%coefficient_ = output%coefficient_ * that
end procedure

module procedure multiply_real_StressBasisFunction
  output = that
  output%coefficient_ = this * output%coefficient_
end procedure

module procedure divide_StressBasisFunction_real
  output = this
  output%coefficient_ = output%coefficient_ / that
end procedure

module procedure interpolate_coefficients_StressBasisFunction
  complex(dp) :: elements(3,3)
  
  integer :: i,j
  
  do i=1,3
    do j=i,3
      elements(j,i) = interpolator%overlap(monomial, this%elements_(j,i)) &
                  & * this%coefficient_
    enddo
  enddo
  
  elements(1,2) = elements(2,1)
  elements(1,3) = elements(3,1)
  elements(2,3) = elements(3,2)
  
  output = mat(elements)
end procedure

module procedure calculate_dynamical_matrices_StressBasisFunction
  type(DynamicalMatrix), allocatable :: matrices(:)
  
  integer :: i,j,k
  
  output = [( StressDynamicalMatrix(anharmonic_data%structure%no_atoms), &
            & i=1,                                                       &
            & size(qpoints)                                              )]
  
  do i=1,3
    do j=1,3
      matrices = calculate_dynamical_matrices( this%elements_(j,i),   &
                                             & qpoints,               &
                                             & thermal_energy,        &
                                             & subspaces,             &
                                             & subspace_bases,        &
                                             & subspace_states,       &
                                             & subspaces_in_coupling, &
                                             & anharmonic_data        )
      do k=1,size(output)
        output(k)%elements(j,i) = matrices(k) * this%coefficient_
      enddo
    enddo
  enddo
end procedure

module procedure stress_correction_StressBasisFunction
  real(dp) :: elements(3,3)
  
  integer :: i,j
  
  do i=1,3
    do j=1,3
      elements(j,i) = calculate_correction( this%elements_(j,i), &
           &                                subspaces,           &
           &                                subspace_bases,      &
           &                                subspace_states,     &
           &                                anharmonic_data )    &
           & * this%coefficient_
    enddo
  enddo
  
  output = mat(elements)
end procedure

module procedure read_StressBasisFunction
  type(ComplexMonomial) :: no_monomials(0)
  
  type(ComplexPolynomial) :: elements(3,3)
  
  type(StringArray), allocatable :: sections(:)
  
  type(String) :: element
  
  character(1), parameter :: directions(3) = ['x', 'y', 'z']
  
  integer :: i,j,k
  
  select type(this); type is(StressBasisFunction)
    elements = ComplexPolynomial(no_monomials)
    
    sections = split_into_sections(input)
    do k=1,size(sections)
      element = token(sections(k)%strings(1),3)
      i = first(directions==char(slice(element,1,1)))
      j = first(directions==char(slice(element,2,2)))
      elements(i,j) = ComplexPolynomial(join( sections(k)%strings(2:), &
                                            & delimiter=' + '          ))
    enddo
    
    this = StressBasisFunction(elements)
  class default
    call err()
  end select
end procedure

module procedure write_StressBasisFunction
  character(1), parameter :: directions(3) = ['x', 'y', 'z']
  
  logical :: first_element
  
  integer :: i,j,ialloc
  
  select type(this); type is(StressBasisFunction)
    allocate(output(0), stat=ialloc); call err(ialloc)
    first_element = .true.
    do i=1,3
      do j=1,3
        if (size(this%elements_(i,j))>0) then
          if (.not. first_element) then
            output = [output, str('')]
          endif
          output = [                                                        &
             & output,                                                      &
             & str('Stress component '//directions(i)//directions(j)//':'), &
             & str(this%elements_(i,j)*this%coefficient_)                   ]
          first_element = .false.
        endif
      enddo
    enddo
  class default
    call err()
  end select
end procedure

module procedure new_StressBasisFunction_Strings
  call this%read(input)
end procedure

module procedure new_StressBasisFunction_StringArray
  this = StressBasisFunction(str(input))
end procedure
end submodule
