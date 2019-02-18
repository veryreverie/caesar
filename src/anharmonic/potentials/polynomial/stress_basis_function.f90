! ======================================================================
! Basis functions for generating the stress tensor mapping.
! ======================================================================
module stress_basis_function_module
  use common_module
  
  use anharmonic_common_module
  use basis_function_module
  implicit none
  
  private
  
  public :: StressBasisFunction
  public :: generate_stress_basis_functions
  
  type, extends(Stringsable) :: StressBasisFunction
    type(BasisFunction) :: elements(3,3)
  contains
    ! I/O.
    procedure, public :: read  => read_StressBasisFunction
    procedure, public :: write => write_StressBasisFunction
  end type
  
  interface StressBasisFunction
    module procedure new_StressBasisFunction
    module procedure new_StressBasisFunction_Strings
    module procedure new_StressBasisFunction_StringArray
  end interface
  
  interface generate_stress_basis_functions
    module procedure generate_stress_basis_functions_SubspaceMonomial
  end interface
contains

! Constructor.
function new_StressBasisFunction(elements) result(this)
  implicit none
  
  type(BasisFunction), intent(in) :: elements(3,3)
  type(StressBasisFunction)       :: this
  
  this%elements = elements
end function

! Generate basis functions.
function generate_stress_basis_functions_SubspaceMonomial(subspace_monomial, &
   & structure,complex_modes,real_modes,qpoints,subspaces,                   &
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
  type(StressBasisFunction), allocatable  :: output(:)
  
  ! Monomials, in complex and real representations.
  type(ComplexMonomial), allocatable :: complex_monomials(:)
  type(RealMonomial),    allocatable :: real_monomials(:)
  
  ! The conversion from the complex monomial basis to the real monomial basis,
  !    and vice-versa.
  type(ComplexMatrix) :: complex_to_real_conversion
  type(ComplexMatrix) :: real_to_complex_conversion
  
  ! Each symmetry, acting on scalar basis functions, the stress tensor,
  !    and the combined basis function.
  type(ComplexMatrix)   :: complex_scalar_symmetry
  real(dp), allocatable :: real_scalar_symmetry(:,:)
  real(dp)              :: tensor_symmetry(9,9)
  real(dp), allocatable :: symmetry(:,:)
  type(RealMatrix)      :: projection
  
  ! Polynomial coefficients, in both bases.
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  real(dp),                  allocatable :: real_coefficients(:)
  complex(dp),               allocatable :: complex_coefficients(:)
  
  ! Variables for constructing the output.
  type(RealPolynomial)    :: real_representation
  type(ComplexPolynomial) :: complex_representation
  type(BasisFunction)     :: elements(3,3)
  
  ! Mappings between 3x3 indices and 9 indices.
  integer :: x(9)
  integer :: y(9)
  integer :: z(3,3)
  
  ! Temporary variables.
  integer  :: i,j,k,l,m,n,o,ialloc
  real(dp) :: tensor(3,3)
  
  ! Initialise index mappings.
  x = [1,1,1,2,2,2,3,3,3]
  y = [1,2,3,1,2,3,1,2,3]
  z = reshape([1,4,7,2,5,8,3,6,9], [3,3])
  
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
    output = [StressBasisFunction::]
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
  projection = dblemat(make_identity_matrix(size(real_monomials)*9))
  allocate( symmetry(size(real_monomials)*9,size(real_monomials)*9), &
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
    
    ! Transform the symmetry into real co-ordinates,
    !    and check that it is real.
    complex_scalar_symmetry = complex_to_real_conversion &
                          & * complex_scalar_symmetry    &
                          & * real_to_complex_conversion
    call check_real( complex_scalar_symmetry,           &
                   & 'symmetry in real monomial basis', &
                   & logfile                            )
    real_scalar_symmetry = dble(real(complex_scalar_symmetry))
    
    ! Construct the symmetry acting on the tensor components.
    do j=1,9
      tensor = 0.0_dp
      tensor(x(j),y(j)) = 1.0_dp
      tensor = dble( structure%symmetries(i)%cartesian_tensor         &
                 & * mat(tensor)                                      &
                 & * invert(structure%symmetries(i)%cartesian_tensor) )
      do k=1,9
        tensor_symmetry(k,j) = tensor(x(k),y(k))
      enddo
    enddo
    call check_orthogonal( mat(tensor_symmetry),                     &
                         & 'symmetry in basis of tensor components', &
                         & logfile                                   )
    
    ! Construct the full symmetry, transforming both scalar and tensor
    !    compontents.
    l = 0
    do j=1,9
      do k=1,size(real_monomials)
        l = l+1
        
        o = 0
        do m=1,9
          do n=1,size(real_monomials)
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
    do j=1,9
      real_coefficients = estuff(i)%evec( size(real_monomials)*(j-1)+1 &
                                      & : size(real_monomials)*j       )
      complex_coefficients = cmplx( real_to_complex_conversion &
                                & * vec(real_coefficients)     )
      
      real_representation = RealPolynomial( real_coefficients &
                                        & * real_monomials    )
      complex_representation = ComplexPolynomial( complex_coefficients &
                                              & * complex_monomials    )
      
      elements(x(j),y(j)) = BasisFunction( real_representation,   &
                                         & complex_representation )
      call elements(x(j),y(j))%simplify()
    enddo
    
    output(i) = StressBasisFunction(elements)
  enddo
end function

! Given an orthogonal matrix U, s.t. U^n=I, returns the matrix
!    H = sum(j=0,n-1) U^j / n
! H is symmetric, and is a projection matrix which projects out the
!    eigenvectors of U with eigenvalue 1.
! Uses Horner's rule to calculate the sum with minimal matrix multiplication.
function projection_matrix(input,order) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  integer,          intent(in) :: order
  type(RealMatrix)             :: output
  
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
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_StressBasisFunction(this,input)
  implicit none
  
  class(StressBasisFunction), intent(out) :: this
  type(String),               intent(in)  :: input(:)
  
  type(BasisFunction) :: elements(3,3)
  
  type(StringArray), allocatable :: sections(:)
  
  integer :: i,j,k
  
  select type(this); type is(StressBasisFunction)
    sections = split_into_sections(input)
    k = 0
    do i=1,3
      do j=1,3
        k = k+1
        elements(i,j) = BasisFunction(sections(k)%strings(2:))
      enddo
    enddo
    
    this = StressBasisFunction(elements)
  class default
    call err()
  end select
end subroutine

function write_StressBasisFunction(this) result(output)
  implicit none
  
  class(StressBasisFunction), intent(in) :: this
  type(String), allocatable              :: output(:)
  
  integer :: i,j
  
  select type(this); type is(StressBasisFunction)
    output = [String::]
    do i=1,3
      do j=1,3
        if (i/=1 .or. j/=1) then
          output = [output, str('')]
        endif
        output = [ output,                       &
                 & 'Element ('//i//' '//j//'):', &
                 & str(this%elements(i,j))       ]
      enddo
    enddo
  class default
    call err()
  end select
end function

function new_StressBasisFunction_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)  :: input(:)
  type(StressBasisFunction) :: this
  
  call this%read(input)
end function

impure elemental function new_StressBasisFunction_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(StressBasisFunction)     :: this
  
  this = StressBasisFunction(str(input))
end function
end module
