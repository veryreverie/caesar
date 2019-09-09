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
  implicit none
  
  private
  
  public :: BasisFunction
  public :: BasisFunctionsAndMonomials
  public :: generate_basis_functions
  public :: finalise
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  
  type, extends(Stringsable) :: BasisFunction
    ! The basis function in complex co-ordinates.
    type(ComplexPolynomial), private :: complex_representation_
    
    ! A leading coefficient.
    real(dp), private :: coefficient_
  contains
    procedure, public :: complex_representation
    
    procedure, public :: simplify => simplify_BasisFunction
    
    generic,   public  :: energy =>                                  &
                        & energy_RealModeDisplacement_BasisFunction, &
                        & energy_ComplexModeDisplacement_BasisFunction
    procedure, private :: energy_RealModeDisplacement_BasisFunction
    procedure, private :: energy_ComplexModeDisplacement_BasisFunction
    generic,   public  :: force =>                                  &
                        & force_RealModeDisplacement_BasisFunction, &
                        & force_ComplexModeDisplacement_BasisFunction
    procedure, private :: force_RealModeDisplacement_BasisFunction
    procedure, private :: force_ComplexModeDisplacement_BasisFunction
    
    generic,   public :: braket =>             &
                       & braket_SubspaceState, &
                       & braket_BasisState,    &
                       & braket_BasisStates
    procedure, public :: braket_SubspaceState => &
                       & braket_SubspaceState_BasisFunction
    procedure, public :: braket_BasisState => &
                       & braket_BasisState_BasisFunction
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_BasisFunction
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_BasisFunction
    
    procedure, public :: undisplaced_energy => undisplaced_energy_BasisFunction
    
    procedure, public  :: coefficient => coefficient_BasisFunction
    procedure, public  :: set_coefficient => set_coefficient_BasisFunction
    
    procedure, public :: read  => read_BasisFunction
    procedure, public :: write => write_BasisFunction
  end type
  
  ! Return type for generate_basis_functions
  type :: BasisFunctionsAndMonomials
    type(BasisFunction), allocatable :: basis_functions(:)
    type(RealMonomial),  allocatable :: unique_terms(:)
  end type
  
  interface BasisFunction
    module procedure new_BasisFunction
    module procedure new_BasisFunction_Strings
    module procedure new_BasisFunction_StringArray
  end interface
  
  interface generate_basis_functions
    module procedure generate_basis_functions_SubspaceMonomial
  end interface
  
  interface finalise
    module procedure finalise_BasisFunctions
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
impure elemental function new_BasisFunction(complex_representation) &
   & result(this)
  implicit none
  
  type(ComplexPolynomial), intent(in) :: complex_representation
  type(BasisFunction)                 :: this
  
  this%complex_representation_ = complex_representation
  this%coefficient_            = 1.0_dp
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
   & structure,complex_modes,real_modes,qpoints,subspaces,            &
   & degenerate_symmetries,vscf_basis_functions_only,logfile) result(output)
  implicit none
  
  type(SubspaceMonomial),   intent(in)    :: subspace_monomial
  type(StructureData),      intent(in)    :: structure
  type(ComplexMode),        intent(in)    :: complex_modes(:)
  type(RealMode),           intent(in)    :: real_modes(:)
  type(QpointData),         intent(in)    :: qpoints(:)
  type(DegenerateSubspace), intent(in)    :: subspaces(:)
  type(DegenerateSymmetry), intent(in)    :: degenerate_symmetries(:)
  logical,                  intent(in)    :: vscf_basis_functions_only
  type(OFile),              intent(inout) :: logfile
  type(BasisFunctionsAndMonomials)        :: output
  
  ! Monomials, in complex and real representations.
  type(ComplexMonomial), allocatable :: complex_monomials(:)
  type(RealMonomial),    allocatable :: real_monomials(:)
  
  ! The conversion from the complex monomial basis to the real monomial basis,
  !    and vice-versa.
  type(ComplexMatrix) :: complex_to_real_conversion
  type(ComplexMatrix) :: real_to_complex_conversion
  
  ! Symmetry data.
  type(ComplexMatrix) :: symmetry
  type(ComplexMatrix) :: projection
  
  ! Polynomial coefficients, in both bases.
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  real(dp),                  allocatable :: real_coefficients(:)
  complex(dp),               allocatable :: complex_coefficients(:)
  
  ! Variables for constructing the output.
  type(RealPolynomial),    allocatable :: real_representations(:)
  type(ComplexPolynomial), allocatable :: complex_representations(:)
  
  type(BasisFunction), allocatable :: basis_functions(:)
  type(RealMonomial),  allocatable :: unique_terms(:)
  
  integer            :: unique_term_id
  integer            :: matching_term_location
  type(RealMonomial) :: matching_term
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  if (sum(subspace_monomial%powers)<2) then
    call print_line(CODE_ERROR//': Trying to generate basis functions with &
       &power less than 2.')
    call err()
  endif
  
  ! Generate the real monomials and complex monomials corresponding to the
  !    subspace monomial, with coefficients such that symmetries are unitary.
  real_monomials = generate_real_monomials( subspace_monomial, &
                                          & subspaces,         &
                                          & real_modes,        &
                                          & qpoints            )
  
  complex_monomials = generate_complex_monomials(            &
      & subspace_monomial,                                   &
      & subspaces,                                           &
      & complex_modes,                                       &
      & qpoints,                                             &
      & conserve_momentum=.true.,                            &
      & conserve_subspace_momentum=vscf_basis_functions_only )
  
  if (size(complex_monomials)==0) then
    allocate( basis_functions(0), &
            & unique_terms(0),    &
            & stat=ialloc); call err(ialloc)
     output = BasisFunctionsAndMonomials(basis_functions, unique_terms)
    return
  endif
  
  ! Identify the mappings between complex monomials and real monomials.
  complex_to_real_conversion = basis_conversion_matrix( &
                          & real_monomials,             &
                          & complex_monomials,          &
                          & include_coefficients=.true. )
  
  real_to_complex_conversion = basis_conversion_matrix( &
                          & complex_monomials,          &
                          & real_monomials,             &
                          & include_coefficients=.true. )
  
  ! Construct projection matrix, which has allowed basis functions as
  !    eigenvectors with eigenvalue 1, and sends all other functions to 0.
  projection = cmplxmat(make_identity_matrix(size(complex_monomials)))
  do i=1,size(degenerate_symmetries)
    ! Constuct symmetry in complex monomial co-ordinates.
    symmetry = degenerate_symmetries(i)%calculate_symmetry( &
                              & complex_monomials,          &
                              & complex_modes,              &
                              & include_coefficients=.true. )
    call check_unitary(symmetry, 'symmetry in monomial basis', logfile)
    
    ! Construct the projection matrix for this symmetry,
    !    and multiply the total projection matrix by this.
    projection = projection                                                  &
             & * projection_matrix( symmetry,                                &
             &                      structure%symmetries(i)%symmetry_order() )
  enddo
  call check_hermitian( projection,               &
                      & 'projection matrix',      &
                      & logfile,                  &
                      & ignore_threshold=1e-10_dp )
  
  ! Transform the projection matrix into real co-ordinates,
  !    and check that it is real and symmetric.
  projection = complex_to_real_conversion &
           & * projection                 &
           & * real_to_complex_conversion
  call check_real(projection,'projection matrix',logfile)
  call check_symmetric( real(projection),         &
                      & 'projection_matrix',      &
                      & logfile,                  &
                      & ignore_threshold=1e-10_dp )
  
  ! Diagonalise the projection matrix,
  !    check its eigenvalues are either 0 or 1,
  !    and select only the eigenvectors with eigenvalue 1.
  estuff = diagonalise_symmetric(real(projection))
  if (any(abs(estuff%eval-1)>1e-2_dp .and. abs(estuff%eval)>1e-2_dp)) then
    call print_line(ERROR//': Projection matrix has eigenvalues which are &
       &neither 0 nor 1.')
    call print_line('Eigenvalues:')
    call print_line(estuff%eval)
    call err()
  endif
  estuff = estuff(filter(abs(estuff%eval-1)<1e-2_dp))
  
  ! Construct basis functions from the coefficients.
  allocate( real_representations(size(estuff)),    &
          & complex_representations(size(estuff)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(estuff)
    real_coefficients = estuff(i)%evec
    real_representations(i) = RealPolynomial( real_coefficients &
                                          & * real_monomials    )
    call real_representations(i)%simplify()
    
    complex_coefficients = cmplx( real_to_complex_conversion &
                              & * vec(real_coefficients)     )
    complex_representations(i) = ComplexPolynomial( complex_coefficients &
                                                & * complex_monomials    )
    call complex_representations(i)%simplify()
  enddo
  
  ! Take linear combinations of basis functions such that each basis function
  !    contains at least one term which is in no other basis function.
  allocate( unique_terms(size(estuff)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(estuff)
    ! Identify the largest term in basis function i.
    unique_term_id = maxloc(abs(real_representations(i)%terms%coefficient), 1)
    unique_terms(i) = real_representations(i)%terms(unique_term_id)
    
    ! Subtract a multiple of basis function i from all other basis functions,
    !    such that the coefficient of unique_term_id(i) in all other basis
    !    functions is zero.
    do j=1,size(basis_functions)
      if (j/=i) then
        matching_term_location = first_equivalent( &
                  & real_representations(j)%terms, &
                  & unique_terms(i),               &
                  & compare_real_monomials,        &
                  & default=0                      )
        if (matching_term_location/=0) then
          matching_term = real_representations(j)%terms(matching_term_location)
          real_representations(j) = real_representations(j)   &
                                & - real_representations(i)   &
                                & * matching_term%coefficient &
                                & / unique_terms(i)%coefficient
          complex_representations(j) = complex_representations(j)   &
                                   & - complex_representations(i)   &
                                   & * matching_term%coefficient &
                                   & / unique_terms(i)%coefficient
        endif
      endif
    enddo
  enddo
  
  basis_functions = BasisFunction(complex_representations)
  
  output = BasisFunctionsAndMonomials(basis_functions, unique_terms)
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

function finalise_BasisFunctions(input,subspace,subspace_basis, &
   & anharmonic_data) result(output)
  implicit none
  
  type(BasisFunction),      intent(in) :: input(:)
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(BasisFunction), allocatable     :: output(:)
  
  integer :: order
  
  type(ComplexPolynomial), allocatable :: complex_polynomials(:)
  
  type(ComplexPolynomial) :: complex_representation
  
  type(RealMonomial)    :: no_real_monomials(0)
  type(ComplexMonomial) :: no_complex_monomials(0)
  
  type(ComplexMonomial) :: conjugate
  
  logical, allocatable :: mode_found(:)
  
  complex(dp) :: coefficient
  
  real(dp) :: energy_scale
  
  integer :: i,j,i2,j2,k,ialloc
  
  order = anharmonic_data%potential_expansion_order
  
  ! Collate the input into an array of polynomials, one for each power.
  ! N.B. the terms with power=0 are dropped as they have already been
  !    accounted for in the reference energy.
  complex_polynomials = [( ComplexPolynomial(no_complex_monomials), &
                         & i=1,                                     &
                         & order                                    )]
  do i=1,size(input)
    do j=1,size(input(i)%complex_representation_%terms)
      associate (term => input(i)%complex_representation_%terms(j))
        if (term%total_power()>0) then
          complex_polynomials(term%total_power()) =      &
             &   complex_polynomials(term%total_power()) &
             & + term
        endif
      end associate
    enddo
  enddo
  
  ! Split the polynomials into basis functions, one for each monomial in
  !    each complex polynomial.
  ! If the mode is real, it becomes its own basis function.
  ! If the mode is complex, it and its conjugate are paired into a u+u*
  !    and a i(u-u*) basis function, with a real coefficient.
  allocate( output(sum( [(size(complex_polynomials(i)), i=1, order)] )), &
          & stat=ialloc); call err(ialloc)
  k = 0
  do i=1,order
    energy_scale = subspace_basis%frequency**(0.5_dp*i)
    
    mode_found = [(.false., j=1,size(complex_polynomials(i)%terms))]
    do j=1,size(complex_polynomials(i)%terms)
      if (mode_found(j)) then
        cycle
      endif
      
      conjugate = conjg(complex_polynomials(i)%terms(j))
      
      j2 = 0
      do i2=j,size(complex_polynomials(i)%terms)
        if (compare_complex_monomials(         &
           & complex_polynomials(i)%terms(i2), &
           & conjugate                         )) then
          j2 = i2
          exit
        endif
      enddo
      if (j2==0) then
        call err()
      endif
      
      mode_found(j) = .true.
      mode_found(j2) = .true.
      
      if (j2==j) then
        coefficient = complex_polynomials(i)%terms(j)%coefficient / energy_scale
        complex_representation = ComplexPolynomial( &
                & [complex_polynomials(i)%terms(j)] )
        complex_representation%terms%coefficient = energy_scale
        k = k+1
        output(k) = BasisFunction(complex_representation)
        output(k)%coefficient_ = real(coefficient)
      else
        coefficient = ( complex_polynomials(i)%terms(j)%coefficient    &
                    & + complex_polynomials(i)%terms(j2)%coefficient ) &
                    & / energy_scale
        complex_representation = ComplexPolynomial(    &
                & [ complex_polynomials(i)%terms(j),   &
                &   complex_polynomials(i)%terms(j2) ] )
        complex_representation%terms%coefficient = 0.5_dp*energy_scale
        k = k+1
        output(k) = BasisFunction(complex_representation)
        output(k)%coefficient_ = real(coefficient)
        
        coefficient = ( complex_polynomials(i)%terms(j)%coefficient    &
                    & - complex_polynomials(i)%terms(j2)%coefficient ) &
                    & / cmplx(0.0_dp,0.5_dp*energy_scale,dp)
        complex_representation%terms%coefficient = [1,-1]                  &
                                               & * cmplx(0.0_dp,0.5_dp,dp) &
                                               & * energy_scale
        k = k+1
        output(k) = BasisFunction(complex_representation)
        output(k)%coefficient_ = real(coefficient)
      endif
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Evaluate the energy and forces due to the basis function.
! ----------------------------------------------------------------------
impure elemental function energy_RealModeDisplacement_BasisFunction(this, &
   & displacement) result(output)
  implicit none
  
  class(BasisFunction),        intent(in) :: this
  class(RealModeDisplacement), intent(in) :: displacement
  real(dp)                                :: output
  
  output = this%coefficient_ &
       & * this%complex_representation_%energy(displacement)
end function

impure elemental function energy_ComplexModeDisplacement_BasisFunction(this, &
   & displacement) result(output)
  implicit none
  
  class(BasisFunction),           intent(in) :: this
  class(ComplexModeDisplacement), intent(in) :: displacement
  complex(dp)                                :: output
  
  output = this%coefficient_ &
       & * this%complex_representation_%energy(displacement)
end function

impure elemental function force_RealModeDisplacement_BasisFunction(this, &
   & displacement) result(output)
  implicit none
  
  class(BasisFunction),        intent(in) :: this
  class(RealModeDisplacement), intent(in) :: displacement
  type(RealModeForce)                     :: output
  
  output = this%coefficient_ &
       & * this%complex_representation_%force(displacement)
end function

impure elemental function force_ComplexModeDisplacement_BasisFunction(this, &
   & displacement) result(output)
  implicit none
  
  class(BasisFunction),           intent(in) :: this
  class(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexModeForce)                     :: output
  
  output = this%coefficient_ &
       & * this%complex_representation_%force(displacement)
end function

! ----------------------------------------------------------------------
! Integrate the basis function between two states.
! ----------------------------------------------------------------------
impure elemental subroutine braket_SubspaceState_BasisFunction(this,bra,ket, &
   & anharmonic_data)
  implicit none
  
  class(BasisFunction), intent(inout)        :: this
  class(SubspaceState), intent(in)           :: bra
  class(SubspaceState), intent(in), optional :: ket
  type(AnharmonicData), intent(in)           :: anharmonic_data
  
  ! Perform integration in complex co-ordinates.
  this%complex_representation_%terms = integrate( &
            & bra,                                &
            & this%complex_representation_%terms, &
            & ket,                                &
            & anharmonic_data                     )
end subroutine

impure elemental subroutine braket_BasisState_BasisFunction(this,bra,ket, &
   & subspace,subspace_basis,anharmonic_data)
  implicit none
  
  class(BasisFunction),     intent(inout)        :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  ! Perform integration in complex co-ordinates.
  this%complex_representation_%terms = integrate( &
            & bra,                                &
            & this%complex_representation_%terms, &
            & ket,                                &
            & subspace,                           &
            & subspace_basis,                     &
            & anharmonic_data                     )
end subroutine

impure elemental subroutine braket_BasisStates_BasisFunction(this,states, &
   & thermal_energy,subspace,subspace_basis,anharmonic_data)
  implicit none
  
  class(BasisFunction),     intent(inout) :: this
  class(BasisStates),       intent(in)    :: states
  real(dp),                 intent(in)    :: thermal_energy
  type(DegenerateSubspace), intent(in)    :: subspace
  class(SubspaceBasis),     intent(in)    :: subspace_basis
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  
  ! Perform integration in complex co-ordinates.
  this%complex_representation_%terms = integrate( &
            & states,                             &
            & thermal_energy,                     &
            & this%complex_representation_%terms, &
            & subspace,                           &
            & subspace_basis,                     &
            & anharmonic_data                     )
end subroutine

! ----------------------------------------------------------------------
! Returns the thermal expectation of the basis function.
! ----------------------------------------------------------------------
impure elemental function harmonic_expectation_BasisFunction(this,frequency, &
   & thermal_energy,anharmonic_data) result(output)
  implicit none
  
  class(BasisFunction), intent(in) :: this
  real(dp),             intent(in) :: frequency
  real(dp),             intent(in) :: thermal_energy
  type(AnharmonicData), intent(in) :: anharmonic_data
  real(dp)                         :: output
  
  output = this%coefficient_                                  &
       & * this%complex_representation_%harmonic_expectation( &
       &                 frequency,                           &
       &                 thermal_energy,                      &
       &                 anharmonic_data%anharmonic_supercell )
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

! ----------------------------------------------------------------------
! Converts the basis function to and from a single coefficient.
! The coefficient is the sum of the absolute values of the
!    real representation coefficients, with the sign of the first non-zero
!    real representation coefficient.
! These methods allow for simpler linear algebra with basis functions.
! ----------------------------------------------------------------------
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
