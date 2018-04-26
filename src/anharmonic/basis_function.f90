! ======================================================================
! Generates basis functions.
! ======================================================================
module basis_function_module
  use common_module
  
  use degeneracy_module
  use subspace_monomial_module
  use mode_monomial_module
  use degenerate_symmetry_module
  implicit none
  
  private
  
  public :: BasisFunction
  public :: generate_basis_functions
  
  type, extends(Printable) :: BasisFunction
    ! The basis function in real co-ordinates.
    type(RealPolynomial) :: real_representation
    
    ! The basis function in complex co-ordinates.
    type(ComplexPolynomial) :: complex_representation
    
    ! The term which is non-zero in this basis function but
    !    zero in every other basis function.
    type(RealMonomial)      :: unique_term
  contains
    procedure, public :: to_String => to_String_BasisFunction
  end type
  
  interface generate_basis_functions
    module procedure generate_basis_functions_SubspaceMonomial
  end interface
  
  interface BasisFunction
    module procedure new_BasisFunction_Strings
  end interface
contains

function generate_basis_functions_SubspaceMonomial(coupling,structure, &
   & complex_modes,real_modes,qpoints,subspaces,degenerate_symmetries, &
   & vscf_basis_functions_only,logfile) result(output)
  implicit none
  
  type(SubspaceMonomial),   intent(in)    :: coupling
  type(StructureData),      intent(in)    :: structure
  type(ComplexMode),        intent(in)    :: complex_modes(:)
  type(RealMode),           intent(in)    :: real_modes(:)
  type(QpointData),         intent(in)    :: qpoints(:)
  type(DegenerateModes),    intent(in)    :: subspaces(:)
  type(DegenerateSymmetry), intent(in)    :: degenerate_symmetries(:)
  logical,                  intent(in)    :: vscf_basis_functions_only
  type(OFile),              intent(inout) :: logfile
  type(BasisFunction), allocatable        :: output(:)
  
  ! Coupling data.
  type(ModeMonomial), allocatable :: mode_monomials(:)
  
  ! Symmetry data.
  type(ComplexMatrix) :: symmetry
  type(ComplexMatrix) :: projection
  
  ! Monomials, in complex and real representations.
  type(ComplexMonomial), allocatable :: complex_monomials(:)
  type(ComplexMonomial), allocatable :: unique_complex_monomials(:)
  type(RealMonomial),    allocatable :: real_monomials(:)
  type(RealMonomial),    allocatable :: unique_real_monomials(:)
  
  ! Variables for converting from the basis of monomials corresponding to
  !    mode_monomials into one with no duplicates.
  ! e.g. modes [3,3,7] and [3,7,3] both become u3^2.u7^1.
  integer,  allocatable :: equal_monomials(:)
  real(dp), allocatable :: all_to_unique(:,:)
  
  ! The conversion from complex co-ordinates to real co-ordinates.
  type(ComplexMatrix) :: complex_to_real_conversion
  
  ! Polynomial coefficients, in both bases.
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  real(dp),                  allocatable :: real_coefficients(:)
  complex(dp),               allocatable :: complex_coefficients(:)
  
  ! A list of which term is unique to which basis function.
  integer, allocatable :: unique_term(:)
  
  integer :: i,j,ialloc
  
  if (size(coupling)<2) then
    call print_line(CODE_ERROR//': Trying to generate basis functions with &
       &order less than 2.')
    call err()
  endif
  
  ! Generate every allowed mode coupling within the subspace coupling.
  mode_monomials = generate_mode_coupling( coupling,      &
                                         & subspaces,     &
                                         & complex_modes, &
                                         & qpoints)
  
  ! Convert coupled modes into real monomials with coefficient 1.
  allocate(real_monomials(size(mode_monomials)), stat=ialloc); call err(ialloc)
  do i=1,size(mode_monomials)
    real_monomials(i) = construct_real_monomial(mode_monomials(i),real_modes)
  enddo
  
  ! Filter the mode couplings, to leave only those which conserve momentum.
  if (vscf_basis_functions_only) then
    mode_monomials = mode_monomials(filter(mode_monomials%conserves_vscf))
  else
    mode_monomials = mode_monomials(filter(mode_monomials%conserves_momentum))
  endif
  if (size(mode_monomials)==0) then
    output = [BasisFunction::]
    return
  endif
  
  ! Convert coupled modes into complex monomials with coefficient 1.
  allocate( complex_monomials(size(mode_monomials)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(mode_monomials)
    complex_monomials(i) = construct_complex_monomial( mode_monomials(i), &
                                                     & complex_modes)
  enddo
  
  ! Identify the unique monomials. (Those with all modes the same).
  unique_complex_monomials = complex_monomials(set( complex_monomials, &
                                                  & compare_monomial_modes))
  unique_real_monomials = real_monomials(set( real_monomials, &
                                            & compare_monomial_modes))
  
  ! Set the coefficient of each unique monomial, and construct the mapping
  !    from complex_monomials to unique_complex_monomials.
  ! In order for the symmetry operators to be unitary in both bases,
  !    it is necessary to preserve the L2 norm. As such, the coefficient
  !    of a unique_monomial representing n monomials is sqrt(n).
  allocate( all_to_unique( size(unique_complex_monomials), &
          &                size(complex_monomials)),       &
          & stat=ialloc); call err(ialloc)
  all_to_unique = 0
  do i=1,size(unique_complex_monomials)
    equal_monomials = filter( complex_monomials,      &
                            & compare_monomial_modes, &
                            & unique_complex_monomials(i))
    unique_complex_monomials(i)%coefficient = &
       & sqrt(real(size(equal_monomials),dp))
    all_to_unique(i,equal_monomials) = &
       & 1/real(unique_complex_monomials(i)%coefficient)
  enddo
  
  do i=1,size(unique_real_monomials)
    equal_monomials = filter( real_monomials,         &
                            & compare_monomial_modes, &
                            & unique_real_monomials(i))
    unique_real_monomials(i)%coefficient = sqrt(real(size(equal_monomials),dp))
  enddo
  
  ! Identify the mapping from complex monomials to real monomials,
  !    and check that this conversion is orthonormal (unitary, but m>=n).
  complex_to_real_conversion = conversion_matrix( unique_real_monomials, &
                                                & unique_complex_monomials)
  call check_orthonormal( complex_to_real_conversion,          &
                        & 'complex to real conversion matrix', &
                        & logfile)
  
  ! Construct projection matrix, which has allowed basis functions as
  !    eigenvectors with eigenvalue 1, and sends all other functions to 0.
  projection = cmplxmat(make_identity_matrix(size(unique_complex_monomials)))
  do i=1,size(degenerate_symmetries)
    ! Constuct symmetry in coupled mode co-ordinates,
    !    and transform it into unique-monomial co-ordinates.
    symmetry = mat(all_to_unique) &
           & * degenerate_symmetries(i)%calculate_symmetry(mode_monomials) &
           & * mat(transpose(all_to_unique))
    call check_unitary(symmetry,'coupled symmetry matrix',logfile)
    
    ! Construct the projection matrix for this symmetry,
    !    and multiply the total projection matrix by this.
    projection = projection                   &
             & * projection_matrix( symmetry, &
             &                      structure%symmetries(i)%symmetry_order())
  enddo
  call check_hermitian( projection,          &
                      & 'projection_matrix', &
                      & logfile,             &
                      & ignore_threshold=1e-10_dp)
  
  ! Transform the projection matrix into real co-ordinates,
  !    and check that it is real and symmetric.
  projection = complex_to_real_conversion &
           & * projection                 &
           & * hermitian(complex_to_real_conversion)
  call check_real(projection,'projection_matrix',logfile)
  call check_symmetric( real(projection),    &
                      & 'projection_matrix', &
                      & logfile,             &
                      & ignore_threshold=1e-10_dp)
  
  ! Diagonalise the projection matrix,
  !    check its eigenvalues are either 0 or 1,
  !    and select only the eigenvectors with eigenvalue 1.
  estuff = diagonalise_symmetric(real(projection))
  if (any(abs(estuff%eval-1)>1e-2_dp .and. abs(estuff%eval)>1e-2_dp)) then
    call print_line(ERROR//': Projection matrix has eigenvalues which are &
       &neither 0 nor 1.')
    call err()
  endif
  estuff = estuff(filter(abs(estuff%eval-1)<1e-2_dp))
  
  ! Take linear combinations of basis functions such that each basis function
  !    contains at least term which is in no other basis function.
  allocate(unique_term(size(estuff)), stat=ialloc); call err(ialloc)
  do i=1,size(estuff)
    ! Identify the largest term in basis function i.
    unique_term(i) = maxloc(abs(estuff(i)%evec),1)
    
    ! Subtract a multiple of basis function i from all other basis functions,
    !    such that the coefficient of unique_term(i) in all other basis
    !    functions is zero.
    do j=1,size(estuff)
      if (j/=i) then
        estuff(j)%evec = estuff(j)%evec                 &
                     & - estuff(i)%evec                 &
                     & * estuff(j)%evec(unique_term(i)) &
                     & / estuff(i)%evec(unique_term(i))
      endif
    enddo
  enddo
  
  ! Construct basis functions from the coefficients.
  allocate(output(size(estuff)), stat=ialloc); call err(ialloc)
  do i=1,size(estuff)
    real_coefficients = estuff(i)%evec
    complex_coefficients = cmplx( hermitian(complex_to_real_conversion) &
                              & * vec(real_coefficients))
    output(i)%real_representation = &
       & RealPolynomial(real_coefficients * unique_real_monomials)
    output(i)%complex_representation = &
       & ComplexPolynomial(complex_coefficients * unique_complex_monomials)
    output(i)%unique_term = output(i)%real_representation%terms(unique_term(i))
  enddo
contains
  ! Lambda for comparing monomials.
  function compare_monomial_modes(this,that) result(output)
    implicit none
    
    class(*), intent(in) :: this
    class(*), intent(in) :: that
    logical              :: output
    
    select type(this)
    type is(ComplexMonomial)
      select type(that)
      type is(ComplexMonomial)
        if (size(this%modes)/=size(that%modes)) then
          output = .false.
        else
          output = all( this%modes%id==that%modes%id .and. &
                      & this%modes%power==that%modes%power)
        endif
      class default
        call err()
      end select
    type is(RealMonomial)
      select type(that)
      type is(RealMonomial)
        if (size(this%modes)/=size(that%modes)) then
          output = .false.
        else
          output = all( this%modes%id==that%modes%id .and. &
                      & this%modes%power==that%modes%power)
        endif
      class default
        call err()
      end select
    end select
  end function
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
! I/O.
! ----------------------------------------------------------------------
function to_String_BasisFunction(this) result(output)
  implicit none
  
  class(BasisFunction), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  integer :: i,j,no_lines,ialloc
  
  no_lines = 4                              &
         & + size(this%real_representation) &
         & + size(this%complex_representation)
  
  allocate(output(no_lines), stat=ialloc); call err(ialloc)
  output(1) = 'Real monomial unique to this basis function:'
  output(2) = this%unique_term
  output(3) = 'Basis function in real co-ordinates:'
  j = 3
  do i=1,size(this%real_representation)
    j = j+1
    output(j) = this%real_representation%terms(i)
  enddo
  j = j+1
  output(j) = 'Basis function in complex co-ordinates:'
  do i=1,size(this%complex_representation)
    j = j+1
    output(j) = this%complex_representation%terms(i)
  enddo
end function

function new_BasisFunction_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(BasisFunction)      :: this
  
  integer :: partition_line
  
  integer :: i,ialloc
  
  if (size(input)<4) then
    call print_line(ERROR//': Basis function file shorter than expected.')
    call err()
  endif
  
  this%unique_term = RealMonomial(input(2))
  
  ! Locate the line between real terms and complex terms.
  do i=4,size(input)
    if (size(split(input(i)))>1) then
      partition_line = i
      exit
    endif
  enddo
  
  this%real_representation = RealPolynomial(input(3:partition_line-1))
  this%complex_representation = ComplexPolynomial(input(partition_line+1:))
end function
end module
