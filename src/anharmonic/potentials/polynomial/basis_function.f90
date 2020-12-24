! ======================================================================
! Basis functions for generating the Born-Oppenheimer surface mapping.
! ======================================================================
! N.B. basis functions are initially generated with both a real and complex
!    representation.
! The real representation is used for fitting, and the complex representation
!    is used for calculating overlap integrals.
module basis_function_module
  use common_module
  
  use states_module
  use anharmonic_common_module
  
  use polynomial_interpolator_module
  use polynomial_dynamical_matrices_module
  use polynomial_symmetry_module
  implicit none
  
  private
  
  public :: BasisFunction
  public :: BasisFunctions
  public :: generate_basis_functions
  public :: optimise
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  
  type, extends(PotentialBase) :: BasisFunction
    ! The basis function in complex co-ordinates.
    type(ComplexPolynomial), private :: complex_representation_
    
    ! A leading coefficient.
    real(dp), private :: coefficient_
  contains
    procedure, public, nopass :: representation => &
                               & representation_BasisFunction
    
    procedure, public :: complex_representation
    
    procedure, public :: undisplaced_energy => undisplaced_energy_BasisFunction
    
    procedure, public :: zero_energy => zero_energy_BasisFunction
    procedure, public :: add_constant => add_constant_BasisFunction
    
    procedure, public :: energy_RealModeDisplacement => &
                       & energy_RealModeDisplacement_BasisFunction
    procedure, public :: energy_ComplexModeDisplacement => &
                       & energy_ComplexModeDisplacement_BasisFunction
    
    procedure, public :: force_RealModeDisplacement => &
                       & force_RealModeDisplacement_BasisFunction
    procedure, public :: force_ComplexModeDisplacement => &
                       & force_ComplexModeDisplacement_BasisFunction
    
    procedure, public :: braket_SubspaceBraKet => &
                       & braket_SubspaceBraKet_BasisFunction
    procedure, public :: braket_BasisState => &
                       & braket_BasisState_BasisFunction
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_BasisFunction
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_BasisFunction
    
    procedure, public :: potential_energy_SubspaceBraKet
    procedure, public :: potential_energy_BasisState
    
    procedure, public :: simplify => simplify_BasisFunction
    
    procedure, public :: power => power_BasisFunction
    procedure, public :: coefficient => coefficient_BasisFunction
    procedure, public :: set_coefficient => set_coefficient_BasisFunction
    
    procedure, public :: interpolate_coefficient => &
                       & interpolate_coefficient_BasisFunction
    procedure, public :: add_interpolated_contribution => &
                       & add_interpolated_contribution_BasisFunction
    procedure, public :: add_harmonic_contribution => &
                       & add_harmonic_contribution_BasisFunction
    procedure, public :: calculate_dynamical_matrices => &
                       & calculate_dynamical_matrices_BasisFunction
    procedure, public :: energy_correction => &
                       & energy_correction_BasisFunction
    
    procedure, public :: terms => terms_BasisFunction
    
    procedure, public :: read  => read_BasisFunction
    procedure, public :: write => write_BasisFunction
  end type
  
  ! An array of type BasisFunction.
  type :: BasisFunctions
    type(BasisFunction), allocatable :: basis_functions(:)
  end type
  
  interface BasisFunction
    module procedure new_BasisFunction
    module procedure new_BasisFunction_PotentialBase
    module procedure new_BasisFunction_Strings
    module procedure new_BasisFunction_StringArray
  end interface
  
  interface generate_basis_functions
    module procedure generate_basis_functions_SubspaceMonomial
    module procedure generate_basis_functions_ComplexMonomials
  end interface
  
  interface optimise
    module procedure optimise_BasisFunctions
  end interface
  
  interface operator(*)
    module procedure multiply_BasisFunction_real
    module procedure multiply_real_BasisFunction
  end interface
  
  interface operator(/)
    module procedure divide_BasisFunction_real
  end interface
  
  interface operator(+)
    module procedure add_BasisFunction_BasisFunction
  end interface
  
  interface operator(-)
    module procedure negative_BasisFunction
    
    module procedure subtract_BasisFunction_BasisFunction
  end interface
contains

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
impure elemental function new_BasisFunction(complex_representation, &
   & coefficient) result(this)
  implicit none
  
  type(ComplexPolynomial), intent(in)           :: complex_representation
  real(dp),                intent(in), optional :: coefficient
  type(BasisFunction)                           :: this
  
  this%complex_representation_ = complex_representation
  if (present(coefficient)) then
    this%coefficient_ = coefficient
  else
    this%coefficient_ = 1.0_dp
  endif
end function

recursive function new_BasisFunction_PotentialBase(input) result(this)
  implicit none
  
  class(PotentialBase), intent(in) :: input
  type(BasisFunction)              :: this
  
  select type(input); type is(PotentialBasePointer)
    this = new_BasisFunction_PotentialBase(input%potential())
  type is(PotentialPointer)
    this = new_BasisFunction_PotentialBase(input%potential())
  type is(BasisFunction)
    this = input
  class default
    call err()
  end select
end function

! Type representation.
impure elemental function representation_BasisFunction() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'Polynomial basis function'
end function

! ----------------------------------------------------------------------
! Getters.
! ----------------------------------------------------------------------
impure elemental function complex_representation(this) result(output)
  implicit none
  
  class(BasisFunction), intent(in) :: this
  type(ComplexPolynomial)          :: output
  
  output = this%coefficient_*this%complex_representation_
end function

! ----------------------------------------------------------------------
! Generate basis functions.
! ----------------------------------------------------------------------
function generate_basis_functions_SubspaceMonomial(subspace_monomial, &
   & structure,complex_modes,qpoints,subspaces,degenerate_symmetries, &
   & vscf_basis_functions_only,logfile) result(output) 
  implicit none
  
  type(SubspaceMonomial),   intent(in)    :: subspace_monomial
  type(StructureData),      intent(in)    :: structure
  type(ComplexMode),        intent(in)    :: complex_modes(:)
  type(QpointData),         intent(in)    :: qpoints(:)
  type(DegenerateSubspace), intent(in)    :: subspaces(:)
  type(DegenerateSymmetry), intent(in)    :: degenerate_symmetries(:)
  logical,                  intent(in)    :: vscf_basis_functions_only
  type(OFile),              intent(inout) :: logfile
  type(BasisFunctions)                    :: output
  
  type(ComplexMonomial), allocatable :: complex_monomials(:)
  
  if (sum(subspace_monomial%powers)<2) then
    call print_line(CODE_ERROR//': Trying to generate basis functions with &
       &power less than 2.')
    call err()
  endif
  
  ! Generate the complex monomials corresponding to the subspace monomial,
  !    with coefficients such that symmetries are unitary.
  complex_monomials = generate_complex_monomials(            &
      & subspace_monomial,                                   &
      & subspaces,                                           &
      & complex_modes,                                       &
      & qpoints,                                             &
      & conserve_momentum=.true.,                            &
      & conserve_subspace_momentum=vscf_basis_functions_only )
  
  output = generate_basis_functions( complex_monomials,     &
                                   & structure,             &
                                   & complex_modes,         &
                                   & qpoints,               &
                                   & degenerate_symmetries, &
                                   & logfile                )
end function

function generate_basis_functions_ComplexMonomials(complex_monomials, &
   & structure,complex_modes,qpoints,degenerate_symmetries,logfile)   &
   & result(output) 
  implicit none
  
  type(ComplexMonomial),    intent(in)              :: complex_monomials(:)
  type(StructureData),      intent(in)              :: structure
  type(ComplexMode),        intent(in)              :: complex_modes(:)
  type(QpointData),         intent(in)              :: qpoints(:)
  type(DegenerateSymmetry), intent(in)              :: degenerate_symmetries(:)
  type(OFile),              intent(inout), optional :: logfile
  type(BasisFunctions)                              :: output
  
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
  
  allocate( output%basis_functions(size(estuff)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(estuff)
    ! Construct a basis function from the eigenvectors.
    output%basis_functions(i) = BasisFunction(                       &
       & basis_conversion%vector_from_basis(estuff(i)%evec, 1e-4_dp) )
  enddo
end function

! Given a unitary matrix U, s.t. U^n=I, returns the matrix
!    H = sum(j=0,n-1) U^j / n
! H is Hermitian, and is a projection matrix which projects out the
!    eigenvectors of U with eigenvalue 1.
! Uses Horner's rule to calculate the sum with minimal matrix multiplication.
function projection_matrix(input,order) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  integer,             intent(in) :: order
  type(ComplexMatrix)             :: output
  
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
end function

! ----------------------------------------------------------------------
! Simplify the basis function.
! ----------------------------------------------------------------------
impure elemental subroutine simplify_BasisFunction(this)
  implicit none
  
  class(BasisFunction), intent(inout) :: this
  
  call this%complex_representation_%simplify()
end subroutine

function optimise_BasisFunctions(input,subspace,subspace_basis, &
   & old_subspace_potential,anharmonic_data) result(output)
  implicit none
  
  type(BasisFunction),      intent(in)           :: input(:)
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  class(PotentialData),     intent(in), optional :: old_subspace_potential
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(BasisFunction), allocatable               :: output(:)
  
  integer :: order
  
  type(SubspaceMonomial)               :: subspace_monomial
  type(ComplexMonomial),   allocatable :: complex_monomials(:)
  type(ComplexPolynomial), allocatable :: complex_polynomials(:)
  
  type(ComplexMonomial) :: no_complex_monomials(0)
  
  type(PotentialBasePointer), allocatable :: basis_functions_(:)
  type(BasisFunction),        allocatable :: basis_functions(:)
  type(ComplexPolynomial),    allocatable :: basis_polynomials(:)
  
  type(ComplexPolynomial), allocatable :: order_polynomials(:)
  
  integer :: i,j,ialloc
  
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
             & anharmonic_data%complex_modes, &
             & anharmonic_data%qpoints        ))) then
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
  
  allocate(output(0), stat=ialloc); call err(ialloc)
  do i=1,order
    if (allocated(basis_polynomials)) then
      order_polynomials = basis_polynomials(                             &
         & filter([( basis_polynomials(j)%terms(1)%total_power()==i,     &
         &           j=1,                                                &
         &           size(basis_polynomials)                         )]) )
    else
      subspace_monomial = SubspaceMonomial( ids    = [subspace%id], &
                                          & powers = [i]            )
      complex_monomials = generate_complex_monomials( &
          & subspace_monomial,                        &
          & [subspace],                               &
          & anharmonic_data%complex_modes,            &
          & anharmonic_data%qpoints,                  &
          & conserve_momentum=.true.,                 &
          & conserve_subspace_momentum=.true.         )
      order_polynomials = construct_basis_polynomials( &
                                  & complex_monomials, &
                                  & subspace,          &
                                  & subspace_basis,    &
                                  & anharmonic_data    )
    endif
    
    if (size(order_polynomials)>0) then
      output = [ output,                                              &
               & fit_basis_functions( complex_polynomials(i)%terms,   &
               &                      order_polynomials             ) ]
    endif
  enddo
end function

! Construct the basis functions at a given order.
function construct_basis_polynomials(monomials,subspace,subspace_basis, &
   & anharmonic_data) result(output)
  implicit none
  
  type(ComplexMonomial),    intent(in) :: monomials(:)
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(ComplexPolynomial), allocatable :: output(:)
  
  integer :: order
  
  type(BasisFunctions) :: basis_functions
  
  type(DegenerateSymmetry), allocatable :: symmetries(:)
  
  real(dp) :: energy_scale
  
  integer :: i
  
  symmetries = subspace_basis%select_symmetries( &
        & anharmonic_data%degenerate_symmetries, &
        & anharmonic_data                        )
  
  basis_functions = generate_basis_functions( &
             & monomials,                     &
             & anharmonic_data%structure,     &
             & anharmonic_data%complex_modes, &
             & anharmonic_data%qpoints,       &
             & symmetries                     )
  
  output = [( basis_functions%basis_functions(i)%complex_representation_, &
            & i=1,                                                        &
            & size(basis_functions%basis_functions)                       )]
  
  
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

! Construct the basis functions at a given order.
function fit_basis_functions(monomials,basis_functions) result(output) 
  implicit none
  
  type(ComplexMonomial),    intent(in) :: monomials(:)
  type(ComplexPolynomial),  intent(in) :: basis_functions(:)
  type(BasisFunction), allocatable     :: output(:)
  
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
                          &                 monomials%coefficient )) )
end function

! ----------------------------------------------------------------------
! Evaluate the energy and forces due to the basis function.
! ----------------------------------------------------------------------
impure elemental function energy_RealModeDisplacement_BasisFunction(this, &
   & displacement) result(output)
  implicit none
  
  class(BasisFunction),       intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  real(dp)                               :: output
  
  output = real( this%coefficient_                                 &
             & * this%complex_representation_%energy(displacement) )
end function

impure elemental function energy_ComplexModeDisplacement_BasisFunction(this, &
   & displacement) result(output)
  implicit none
  
  class(BasisFunction),          intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  complex(dp)                               :: output
  
  output = this%coefficient_ &
       & * this%complex_representation_%energy(displacement)
end function

impure elemental function force_RealModeDisplacement_BasisFunction(this, &
   & displacement) result(output)
  implicit none
  
  class(BasisFunction),       intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  type(RealModeForce)                    :: output
  
  output = this%coefficient_ &
       & * this%complex_representation_%force(displacement)
end function

impure elemental function force_ComplexModeDisplacement_BasisFunction(this, &
   & displacement) result(output)
  implicit none
  
  class(BasisFunction),          intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexModeForce)                    :: output
  
  output = this%coefficient_ &
       & * this%complex_representation_%force(displacement)
end function

! ----------------------------------------------------------------------
! Integrate the basis function between two states.
! ----------------------------------------------------------------------
impure elemental subroutine braket_SubspaceBraKet_BasisFunction(this,braket, &
   & whole_subspace,anharmonic_data)
  implicit none
  
  class(BasisFunction),  intent(inout)        :: this
  class(SubspaceBraKet), intent(in)           :: braket
  logical,               intent(in), optional :: whole_subspace
  type(AnharmonicData),  intent(in)           :: anharmonic_data
  
  ! Perform integration in complex co-ordinates.
  call integrate( this%complex_representation_%terms, &
                & braket,                             &
                & anharmonic_data                     )
end subroutine

impure elemental subroutine braket_BasisState_BasisFunction(this,bra,ket, &
   & subspace,subspace_basis,whole_subspace,anharmonic_data)
  implicit none
  
  class(BasisFunction),     intent(inout)        :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  logical,                  intent(in), optional :: whole_subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  ! Perform integration in complex co-ordinates.
  call integrate( this%complex_representation_%terms, &
                & bra,                                &
                & ket,                                &
                & subspace,                           &
                & subspace_basis,                     &
                & anharmonic_data                     )
end subroutine

impure elemental subroutine braket_BasisStates_BasisFunction(this,states, &
   & subspace,subspace_basis,whole_subspace,anharmonic_data)
  implicit none
  
  class(BasisFunction),     intent(inout)        :: this
  class(BasisStates),       intent(inout)        :: states
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  logical,                  intent(in), optional :: whole_subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  integer :: i
  
  ! Perform integration in complex co-ordinates.
  do i=1,size(this%complex_representation_%terms) 
    call integrate( this%complex_representation_%terms(i), &
                  & states,                                &
                  & subspace,                              &
                  & subspace_basis,                        &
                  & anharmonic_data                        )
  enddo
end subroutine

! ----------------------------------------------------------------------
! Returns the thermal expectation of the basis function.
! ----------------------------------------------------------------------
impure elemental function harmonic_expectation_BasisFunction(this,frequency, &
   & thermal_energy,supercell_size,anharmonic_data) result(output)
  implicit none
  
  class(BasisFunction), intent(in) :: this
  real(dp),             intent(in) :: frequency
  real(dp),             intent(in) :: thermal_energy
  integer,              intent(in) :: supercell_size
  type(AnharmonicData), intent(in) :: anharmonic_data
  real(dp)                         :: output
  
  output = this%coefficient_                                                  &
       & * this%complex_representation_%harmonic_expectation( frequency,      &
       &                                                      thermal_energy, &
       &                                                      supercell_size  )
end function

function potential_energy_SubspaceBraKet(this,braket,anharmonic_data) &
   & result(output) 
  implicit none
  
  class(BasisFunction),  intent(in) :: this
  class(SubspaceBraKet), intent(in) :: braket
  type(AnharmonicData),  intent(in) :: anharmonic_data
  real(dp)                          :: output
  
  output = this%coefficient_                                         &
       & * real(integrate_to_constant( this%complex_representation_, &
       &                               braket,                       &
       &                               anharmonic_data               ))
end function

function potential_energy_BasisState(this,bra,ket,subspace,subspace_basis, &
   & anharmonic_data) result(output) 
  implicit none
  
  class(BasisFunction),     intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  output = this%coefficient_                                         &
       & * real(integrate_to_constant( this%complex_representation_, &
       &                               bra,                          &
       &                               ket,                          &
       &                               subspace,                     &
       &                               subspace_basis,               &
       &                               anharmonic_data               ))
end function

! ----------------------------------------------------------------------
! Returns the energy at zero displacement.
! ----------------------------------------------------------------------
impure elemental function undisplaced_energy_BasisFunction(this) result(output)
  implicit none
  
  class(BasisFunction), intent(in) :: this
  real(dp)                         :: output
  
  type(RealSingleDisplacement) :: zero_displacement(0)
  
  output = this%energy(RealModeDisplacement(zero_displacement))
end function

impure elemental subroutine zero_energy_BasisFunction(this)
  implicit none
  
  class(BasisFunction), intent(inout) :: this
  
  call print_line(CODE_ERROR//': zero_energy() cannot be called for type &
     &BasisFunction.')
  call err()
end subroutine

impure elemental subroutine add_constant_BasisFunction(this,input) 
  implicit none
  
  class(BasisFunction), intent(inout) :: this
  real(dp),             intent(in)    :: input
  
  call print_line(CODE_ERROR//': add_constant() cannot be called for type &
     &BasisFunction.')
  call err()
end subroutine

! ----------------------------------------------------------------------
! Converts the basis function to and from a single coefficient.
! These methods allow for simpler linear algebra with basis functions.
! ----------------------------------------------------------------------
impure elemental function power_BasisFunction(this) result(output)
  implicit none
  
  class(BasisFunction), intent(in) :: this
  integer                          :: output
  
  output = this%complex_representation_%terms(1)%total_power()
end function

impure elemental function coefficient_BasisFunction(this) result(output)
  implicit none
  
  class(BasisFunction), intent(in) :: this
  real(dp)                         :: output
  
  output = this%coefficient_
end function

impure elemental subroutine set_coefficient_BasisFunction(this,coefficient)
  implicit none
  
  class(BasisFunction), intent(inout) :: this
  real(dp),             intent(in)    :: coefficient
  
  this%coefficient_ = coefficient
end subroutine

! ----------------------------------------------------------------------
! Arithmetic.
! ----------------------------------------------------------------------
impure elemental function multiply_BasisFunction_real(this,that) result(output)
  implicit none
  
  type(BasisFunction), intent(in) :: this
  real(dp),            intent(in) :: that
  type(BasisFunction)             :: output
  
  output = this
  output%coefficient_ = output%coefficient_ * that
end function

impure elemental function multiply_real_BasisFunction(this,that) result(output)
  implicit none
  
  real(dp),            intent(in) :: this
  type(BasisFunction), intent(in) :: that
  type(BasisFunction)             :: output
  
  output = that
  output%coefficient_ = this * output%coefficient_
end function

impure elemental function divide_BasisFunction_real(this,that) result(output)
  implicit none
  
  type(BasisFunction), intent(in) :: this
  real(dp),            intent(in) :: that
  type(BasisFunction)             :: output
  
  output = this
  output%coefficient_ = output%coefficient_ / that
end function

impure elemental function add_BasisFunction_BasisFunction(this,that) &
   & result(output)
  implicit none
  
  type(BasisFunction), intent(in) :: this
  type(BasisFunction), intent(in) :: that
  type(BasisFunction)             :: output
  
  output = BasisFunction(                                          &
     & this%complex_representation()+that%complex_representation() )
end function

impure elemental function negative_BasisFunction(this) result(output)
  implicit none
  
  type(BasisFunction), intent(in) :: this
  type(BasisFunction)             :: output
  
  output = this
  output%coefficient_ = -output%coefficient_
end function

impure elemental function subtract_BasisFunction_BasisFunction(this,that) &
   & result(output)
  implicit none
  
  type(BasisFunction), intent(in) :: this
  type(BasisFunction), intent(in) :: that
  type(BasisFunction)             :: output
  
  output = BasisFunction(                                          &
     & this%complex_representation()-that%complex_representation() )
end function

! Calculate the contribution to a given monomial from the interpolation of
!    this basis function.
impure elemental function interpolate_coefficient_BasisFunction(this, &
   & monomial,interpolator) result(output)
  implicit none
  
  class(BasisFunction),         intent(in) :: this
  type(ComplexMonomial),        intent(in) :: monomial
  type(PolynomialInterpolator), intent(in) :: interpolator
  complex(dp)                              :: output
  
  output = interpolator%overlap(monomial, this%complex_representation_) &
       & * this%coefficient_
end function

! Calculate the contribution to this basis function from
!    another basis function, and add this to this basis function's coefficient.
subroutine add_interpolated_contribution_BasisFunction(this,basis_function, &
   & interpolator) 
  implicit none
  
  class(BasisFunction),         intent(inout) :: this
  type(BasisFunction),          intent(in)    :: basis_function
  type(PolynomialInterpolator), intent(in)    :: interpolator
  
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
end subroutine

! Calculate the contribution to this basis function from a set of harmonic
!    dynamical matrices.
subroutine add_harmonic_contribution_BasisFunction(this,dynamical_matrices, &
   & anharmonic_data)
  implicit none
  
  class(BasisFunction),  intent(inout) :: this
  type(DynamicalMatrix), intent(in)    :: dynamical_matrices(:)
  type(AnharmonicData),  intent(in)    :: anharmonic_data
  
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
end subroutine

! Calculate this basis function's contribution to the effective dynamical
!    matrix from which the potential can be interpolated in the large-supercell
!    limit.
function calculate_dynamical_matrices_BasisFunction(this,qpoints, &
   & thermal_energy,subspaces,subspace_bases,subspace_states,     &
   & subspaces_in_coupling,anharmonic_data) result(output) 
  implicit none
  
  class(BasisFunction),     intent(in)    :: this
  type(QpointData),         intent(in)    :: qpoints(:)
  real(dp),                 intent(in)    :: thermal_energy
  type(DegenerateSubspace), intent(in)    :: subspaces(:)
  class(SubspaceBasis),     intent(in)    :: subspace_bases(:)
  class(BasisStates),       intent(inout) :: subspace_states(:)
  integer,                  intent(in)    :: subspaces_in_coupling(:)
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  type(DynamicalMatrix), allocatable      :: output(:)
  
  output = calculate_dynamical_matrices( this%complex_representation_,   &
       &                                 qpoints,                        &
       &                                 thermal_energy,                 &
       &                                 subspaces,                      &
       &                                 subspace_bases,                 &
       &                                 subspace_states,                &
       &                                 subspaces_in_coupling,          &
       &                                 anharmonic_data               ) &
       & * this%coefficient_
end function

! Calculate the correction due to double counting
!    for the interpolated potential.
function energy_correction_BasisFunction(this,subspaces,subspace_bases, &
   & subspace_states,anharmonic_data) result(output) 
  implicit none
  
  class(BasisFunction),     intent(in)    :: this
  type(DegenerateSubspace), intent(in)    :: subspaces(:)
  class(SubspaceBasis),     intent(in)    :: subspace_bases(:)
  class(BasisStates),       intent(inout) :: subspace_states(:)
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  real(dp)                                :: output
  
  output = calculate_correction( this%complex_representation_,   &
       &                         subspaces,                      &
       &                         subspace_bases,                 &
       &                         subspace_states,                &
       &                         anharmonic_data               ) &
       & * this%coefficient_
end function

! Return the monomial terms of this basis function.
function terms_BasisFunction(this) result(output)
  implicit none
  
  class(BasisFunction), intent(in)   :: this
  type(ComplexMonomial), allocatable :: output(:)
  
  output = this%complex_representation_%terms
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_BasisFunction(this,input)
  implicit none
  
  class(BasisFunction), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
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
end subroutine

function write_BasisFunction(this) result(output)
  implicit none
  
  class(BasisFunction), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  type(ComplexPolynomial) :: complex_representation
  
  select type(this); type is(BasisFunction)
    complex_representation = this%complex_representation()
    output = [ str('Basis function in complex co-ordinates:'), &
             & str(complex_representation%terms)               ]
  class default
    call err()
  end select
end function

function new_BasisFunction_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(BasisFunction)      :: this
  
  call this%read(input)
end function

impure elemental function new_BasisFunction_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(BasisFunction)           :: this
  
  this = BasisFunction(str(input))
end function
end module
