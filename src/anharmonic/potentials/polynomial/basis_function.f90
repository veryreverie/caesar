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
  public :: BasisFunctionsAndUniqueTerms
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
  type :: BasisFunctionsAndUniqueTerms
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
  type(BasisFunctionsAndUniqueTerms)      :: output
  
  ! Monomials, in complex and real representations.
  type(ComplexMonomial), allocatable :: complex_monomials(:)
  type(RealMonomial),    allocatable :: real_monomials(:)
  
  integer, allocatable :: conjugates(:)
  
  ! The conversion from the complex monomial basis to the real monomial basis,
  !    and vice-versa.
  type(ComplexMatrix) :: complex_to_real_conversion
  type(ComplexMatrix) :: real_to_complex_conversion
  
  ! Symmetry data.
  type(ComplexMatrix) :: symmetry
  type(ComplexMatrix) :: projection
  type(RealMatrix)    :: basis_projection
  
  ! Polynomial coefficients, in both bases.
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  real(dp),                  allocatable :: real_coefficients(:)
  complex(dp),               allocatable :: complex_coefficients(:)
  
  ! Variables for constructing the output.
  type(RealPolynomial),    allocatable :: real_representations(:)
  type(ComplexPolynomial), allocatable :: complex_representations(:)
  
  type(BasisFunction), allocatable :: basis_functions(:)
  integer,             allocatable :: unique_terms(:)
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  
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
  call check_hermitian( projection,                   &
                      & 'monomial projection matrix', &
                      & logfile,                      &
                      & ignore_threshold=1e-10_dp     )
  
  ! Pair up each complex monomial with its complex conjugate.
  conjugates = find_monomial_conjugates(complex_monomials)
  basis_functions = pair_complex_monomials(complex_monomials, conjugates)
  
  ! Transform the projection matrix to basis functions rather than monomials.
  basis_projection = make_basis_projection(projection, conjugates)
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
          & output%unique_terms(size(estuff)),    &
          & stat=ialloc); call err(ialloc)
  do i=1,size(estuff)
    ! Construct a basis function from the eigenvectors.
    output%basis_functions(i) = generate_basis_function( estuff(i)%evec,    &
                                                       & complex_monomials, &
                                                       & conjugates         )
    
    ! Construct a real monomial which stands in for the basis function.
    output%unique_terms(i) = generate_unique_term( unique_terms(i),   &
                                                 & complex_monomials, &
                                                 & conjugates         )
  enddo
end function

! Convert the projection matrix from the basis of complex monomials to
!    the basis of basis functions.
function make_basis_projection(projection,conjugates) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: projection
  integer,             intent(in) :: conjugates(:)
  type(RealMatrix)                :: output
  
  real(dp), allocatable :: matrix(:,:)
  
  integer :: i,j,ialloc
  
  ! Construct the output matrix.
  ! N.B. each time conjugates(x)<x, x and x' are the other way around,
  !    e.g. when conjugates(j)<j, Pij=element(i,conjugate(j)) and
  !    Pij'=element(i,j).
  allocate( matrix(size(projection,1),size(projection,2)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(projection,1)
    do j=1,size(projection,2)
      if (conjugates(i)==i) then
        if (conjugates(j)==j) then
          ! ui = Pij uj
          matrix(i,j) = real(projection%element(i,j))
        elseif (conjugates(j)>j) then
          ! ui = (Pij+Pij') (uj+uj')/2 + ...
          matrix(i,j) = real( projection%element(i,j)             &
                          & + projection%element(i,conjugates(j)) )
        else
          ! ... + i(Pij-Pij') (uj-uj')/2i
          matrix(i,j) = -aimag( projection%element(i,conjugates(j)) &
                            & - projection%element(i,j)             )
        endif
      elseif (conjugates(i)>i) then
        if (conjugates(j)==j) then
          ! (ui+ui')/2 = (Pij+Pi'j)/2 uj
          matrix(i,j) = real( projection%element(i,j)             &
                          & + projection%element(conjugates(i),j) )/2
        elseif (conjugates(j)>j) then
          ! (ui+ui')/2 = (Pij+Pij'+Pi'j+Pi'j')/2 (uj+uj')/2 + ...
          matrix(i,j) = real(                                    &
             &   projection%element(i,j)                         &
             & + projection%element(i,conjugates(j))             &
             & + projection%element(conjugates(i),j)             &
             & + projection%element(conjugates(i),conjugates(j)) ) / 2
        else
          ! ... + i(Pij-Pij'+Pi'j-Pi'j')/2 (uj-uj')/2i
          matrix(i,j) = -aimag(                                  &
             &   projection%element(i,conjugates(j))             &
             & - projection%element(i,j)                         &
             & + projection%element(conjugates(i),conjugates(j)) &
             & - projection%element(conjugates(i),j)             ) / 2
        endif
      else
        if (conjugates(j)==j) then
          ! (ui-ui')/2i = (Pij-Pi'j)/2i uj
          matrix(i,j) = aimag( projection%element(conjugates(i),j) &
                           & - projection%element(i,j)             )/2
        elseif (conjugates(j)>j) then
          ! (ui-ui')/2i = (Pij+Pij'-Pi'j-Pi'j')/2i (uj+uj')/2 + ...
          matrix(i,j) = aimag(                                   &
             &   projection%element(conjugates(i),j)             &
             & + projection%element(conjugates(i),conjugates(j)) &
             & - projection%element(i,j)                         &
             & - projection%element(i,conjugates(j))             ) / 2
        else
          ! ... +  (Pij-Pij'-Pi'j+Pi'j')/2 (uj-uj')/2i
          matrix(i,j) = real(                                    &
             &   projection%element(conjugates(i),conjugates(j)) &
             & - projection%element(conjugates(i),j)             &
             & - projection%element(i,conjugates(j))             &
             & + projection%element(i,j)                         ) / 2
        endif
      endif
    enddo
  enddo
  
  output = mat(matrix)
end function

! Generate a basis function from an eigenvector.
function generate_basis_function(evec,monomials,conjugates) result(output)
  implicit none
  
  real(dp),              intent(in) :: evec(:)
  type(ComplexMonomial), intent(in) :: monomials(:)
  integer,               intent(in) :: conjugates(:)
  type(BasisFunction)               :: output
  
  real(dp) :: max_coeff
  
  integer, allocatable :: included_terms(:)
  
  type(ComplexMonomial), allocatable :: terms(:)
  
  integer :: i,j,k,ialloc
  
  max_coeff = maxval(abs(evec))
  
  included_terms = filter( abs(evec)>max_coeff/1e4_dp             &
                    & .or. abs(evec(conjugates))>max_coeff/1e4_dp )
  
  allocate(terms(size(included_terms)), stat=ialloc); call err(ialloc)
  do i=1,size(terms)
    j = included_terms(i)
    k = conjugates(j)
    terms(i) = monomials(j)
    if (j==k) then
      ! The monomial is its own pair.
      terms(i)%coefficient = evec(j)*monomials(j)%coefficient
    elseif (j<k) then
      ! The monomial looks like u+ = uc + i*us.
      terms(i)%coefficient =                                          &
         &                           evec(j)*monomials(j)%coefficient &
         & + cmplx(0.0_dp,1.0_dp,dp)*evec(k)*monomials(k)%coefficient
    else
      ! The monomial looks like u- = uc - i*us.
      terms(i)%coefficient =                                          &
         &                           evec(k)*monomials(k)%coefficient &
         & - cmplx(0.0_dp,1.0_dp,dp)*evec(j)*monomials(j)%coefficient
    endif
  enddo
  
  output = BasisFunction(ComplexPolynomial(terms))
end function

! Takes a complex polynomial, and generates a monomial which can be turned
!    into a sampling point.
! It doesn't matter overly much how this is done, as long as the mapping
!    is unique and the resulting sampling points are distinct under symmetry.
! Arbitrarily:
!    -  u+^n * u-^m + u+^m * u-^n      ->   uc^max(m,n) * us^min(m,n)
!    - (u+^n * u-^m - u+^m * u-^n)/i   ->   uc^min(m,n) * us^max(m,n)
function generate_unique_term(id,monomials,conjugates) result(output)
  implicit none
  
  integer,               intent(in) :: id
  type(ComplexMonomial), intent(in) :: monomials(:)
  integer,               intent(in) :: conjugates(:)
  type(RealMonomial)                :: output
  
  associate (term => monomials(id))
    if (conjugates(id)==id) then
      ! The term is its own pair,
      !    so the real monomial is the same as the complex monomial.
      output = RealMonomial(                                           &
         & coefficient = 1.0_dp,                                       &
         & modes = RealUnivariate( id           = term%ids(),          &
         &                         paired_id    = term%paired_ids(),   &
         &                         power        = term%powers(),       &
         &                         paired_power = term%paired_powers() ))
    elseif (conjugates(id)>id) then
      ! The term is a (x+x*)/2 sum of a complex term and its conjugate.
      ! Construct an output with larger powers of uc than us.
      output = RealMonomial(                              &
         & coefficient = 1.0_dp,                          &
         & modes = RealUnivariate(                        &
         &    id           = term%ids(),                  &
         &    paired_id    = term%paired_ids(),           &
         &    power        = max( term%powers(),          &
         &                        term%paired_powers() ), &
         &    paired_power = min( term%powers(),          &
         &                        term%paired_powers() )  ))
    else
      ! The term is a (x-x*)/2i sum of a complex term and its conjugate.
      ! Construct an output with larger powers of us than uc.
      output = RealMonomial(                              &
         & coefficient = 1.0_dp,                          &
         & modes = RealUnivariate(                        &
         &    id           = term%ids(),                  &
         &    paired_id    = term%paired_ids(),           &
         &    power        = min( term%powers(),          &
         &                        term%paired_powers() ), &
         &    paired_power = max( term%powers(),          &
         &                        term%paired_powers() )  ))
    endif
  end associate
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
  integer, allocatable  :: conjugates(:)
  
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
  
  allocate( output(sum( [(size(complex_polynomials(i)), i=1, order)] )), &
          & stat=ialloc); call err(ialloc)
  k = 0
  do i=1,order
    energy_scale = subspace_basis%frequency**(0.5_dp*i)
    
    associate (terms => output(k+1:k+size(complex_polynomials(i)%terms)))
      conjugates = find_monomial_conjugates(complex_polynomials(i)%terms)
      terms = pair_complex_monomials(complex_polynomials(i)%terms, conjugates)
      ! Multiply the monomial coefficients by energy_scale, and divide
      !    the basis function coefficient by energy_scale.
      ! This makes the monomials dimensionless and gives the basis function
      !    coefficient units of energy.
      terms%coefficient_ = terms%coefficient_ / energy_scale
      do j=1,size(terms)
        terms(j)%complex_representation_%terms%coefficient = &
           & terms(j)%complex_representation_%terms%coefficient * energy_scale
      enddo
    end associate
    
    k = k+size(complex_polynomials(i)%terms)
  enddo
end function

! Find the id of the monomial conjugates,
!    s.t. conjg(input(i)) = input(output(i)).
function find_monomial_conjugates(input) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: input(:)
  integer, allocatable              :: output(:)
  
  logical, allocatable :: conjugate_found(:)
  
  type(ComplexMonomial) :: conjugate
  
  integer :: i,j,ialloc
  
  output = [(0, i=1, size(input))]
  conjugate_found = [(.false., i=1, size(input))]
  do i=1,size(input)
    if (conjugate_found(i)) then
      cycle
    endif
    
    conjugate = conjg(input(i))
    do j=1,size(input)
      if (compare_complex_monomials(input(j),conjugate)) then
        output(i) = j
        output(j) = i
        conjugate_found(i) = .true.
        conjugate_found(j) = .true.
        exit
      endif
    enddo
    
    if (output(i)==0) then
      call err()
    endif
  enddo
end function

! Recombine the complex monomials into real basis functions.
! If the mode is real, it becomes its own basis function.
! If the mode is complex, it and its conjugate are paired into a u+u*
!    and a i(u-u*) basis function, both with real coefficients.
function pair_complex_monomials(input,conjugates) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: input(:)
  integer,               intent(in) :: conjugates(:)
  type(BasisFunction), allocatable  :: output(:)
  
  type(ComplexPolynomial) :: representation
  
  integer :: i,j,ialloc
  
  allocate(output(size(input)), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    j = conjugates(i)
    if (j==i) then
      ! If i==j the monomial is its own conjugate.
      ! The basis function is simply the monomial.
      ! The coefficient of the monomial is moved to that of the basis function.
      representation = ComplexPolynomial([input(i)])
      representation%terms%coefficient = 1.0_dp
      output(i) = BasisFunction(representation)
      output(i)%coefficient_ = real(input(i)%coefficient)
    elseif (j>i) then
      ! If i/=j then monomial i and monomial j are conjugates.
      ! If i>j, construct (i+j)/2 basis function.
      representation = ComplexPolynomial([input(i), input(j) ] )
      representation%terms%coefficient = [0.5_dp, 0.5_dp]
      output(i) = BasisFunction(representation)
      output(i)%coefficient_ = real( input(i)%coefficient &
                                 & + input(j)%coefficient )
      
    else
      ! If i/=j then monomial i and monomial j are conjugates.
      ! If i>j, construct (j-i)/2i basis function.
      representation = ComplexPolynomial([input(j), input(i) ] )
      representation%terms%coefficient = [1,-1] * cmplx(0.0_dp,0.5_dp,dp)
      output(i) = BasisFunction(representation)
      output(i)%coefficient_ = real( ( input(j)%coefficient    &
                                 &   - input(i)%coefficient )  &
                                 &   / cmplx(0.0_dp,1.0_dp,dp) )
    endif
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
  class(BasisStates),       intent(inout) :: states
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
