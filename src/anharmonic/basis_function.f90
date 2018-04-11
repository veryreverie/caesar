! ======================================================================
! Generates basis functions.
! ======================================================================
module basis_function_module
  use common_module
  
  use degeneracy_module
  use coupling_module
  use degenerate_symmetry_module
  implicit none
  
  private
  
  public :: BasisFunction
  public :: generate_basis_functions
  
  type :: BasisFunction
    type(RealPolynomial)    :: real_representation
    type(ComplexPolynomial) :: complex_representation
  end type
contains

function generate_basis_functions(coupling,structure,complex_modes, &
   & real_modes,qpoints,subspaces,degenerate_symmetries,            &
   & vscf_basis_functions_only,logfile) result(output)
  implicit none
  
  type(CoupledSubspaces),   intent(in)    :: coupling
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
  type(DegenerateModes), allocatable :: coupled_subspaces(:)
  type(CoupledModes),    allocatable :: coupled_modes(:)
  
  ! Symmetry data.
  type(ComplexMatrix) :: symmetry
  integer             :: symmetry_order
  type(ComplexMatrix) :: projection
  
  ! Monomials, in complex and real representations.
  type(ComplexMonomial), allocatable :: monomials(:)
  type(ComplexMonomial), allocatable :: unique_monomials(:)
  type(RealMonomial),    allocatable :: real_monomials(:)
  type(RealMonomial),    allocatable :: unique_real_monomials(:)
  
  ! Variables for converting from the basis of monomials corresponding to
  !    coupled_modes into one with no duplicates.
  ! e.g. modes [3,3,7] and [3,7,3] both become u3^2.u7^1.
  integer,  allocatable :: equal_monomials(:)
  real(dp), allocatable :: all_to_unique(:,:)
  
  ! The conversion from complex co-ordinates to real co-ordinates.
  type(ComplexMatrix) :: complex_to_real_conversion
  
  ! Polynomial coefficients, in both bases.
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  real(dp),                  allocatable :: real_coefficients(:)
  complex(dp),               allocatable :: complex_coefficients(:)
  
  integer :: i,ialloc
  
  if (size(coupling)<2) then
    call print_line(CODE_ERROR//': Trying to generate basis functions with &
       &order less than 2.')
    call err()
  endif
  
  ! List the subspaces which are coupled together in this coupling.
  coupled_subspaces = coupling%coupled_subspaces(subspaces)
  
  ! Generate every allowed mode coupling within the coupled subspaces.
  coupled_modes = generate_mode_coupling( coupled_subspaces, &
                                        & complex_modes,     &
                                        & qpoints)
  
  ! Convert coupled modes into RealMonomials with coefficient 1.
  allocate(real_monomials(size(coupled_modes)), stat=ialloc); call err(ialloc)
  do i=1,size(coupled_modes)
    real_monomials(i) = construct_real_monomial(coupled_modes(i),real_modes)
  enddo
  
  if (vscf_basis_functions_only) then
    coupled_modes = coupled_modes(filter(coupled_modes%conserves_vscf))
  else
    coupled_modes = coupled_modes(filter(coupled_modes%conserves_momentum))
  endif
  
  if (size(coupled_modes)==0) then
    output = [BasisFunction::]
    return
  endif
  
  ! Convert coupled modes into Monomials with coefficient 1.
  allocate(monomials(size(coupled_modes)), stat=ialloc); call err(ialloc)
  do i=1,size(coupled_modes)
    monomials(i) = construct_complex_monomial(coupled_modes(i),complex_modes)
  enddo
  
  ! Identify the unique monomials. (Those with all modes the same).
  unique_monomials = monomials(set(monomials,compare_monomial_modes))
  unique_real_monomials = real_monomials(set(real_monomials,compare_monomial_modes))
  
  ! Set the coefficient of each unique monomial, and construct the mapping
  !    from monomials to unique_monomials.
  ! In order for the symmetry operators to be unitary in both bases,
  !    it is necessary to preserve the L2 norm. As such, the coefficient
  !    of a unique_monomial representing n monomials is sqrt(n).
  allocate( all_to_unique(size(unique_monomials),size(monomials)), &
          & stat=ialloc); call err(ialloc)
  all_to_unique = 0
  do i=1,size(unique_monomials)
    equal_monomials = filter( monomials,              &
                            & compare_monomial_modes, &
                            & unique_monomials(i))
    unique_monomials(i)%coefficient = sqrt(real(size(equal_monomials),dp))
    all_to_unique(i,equal_monomials) = 1/real(unique_monomials(i)%coefficient)
  enddo
  
  do i=1,size(unique_real_monomials)
    equal_monomials = filter( real_monomials,         &
                            & compare_monomial_modes, &
                            & unique_real_monomials(i))
    unique_real_monomials(i)%coefficient = sqrt(real(size(equal_monomials),dp))
  enddo
  
  ! Identify the mapping from complex monomials to real monomials.
  complex_to_real_conversion = conversion_matrix( unique_real_monomials, &
                                                & unique_monomials)
  
  ! Check that the complex to real conversion is unitary.
  call check_orthonormal(complex_to_real_conversion,logfile)
  
  ! Construct projection matrix, which has allowed basis functions as
  !    eigenvectors with eigenvalue 1, and sends all other functions to 0.
  projection = cmplxmat(make_identity_matrix(size(unique_monomials)))
  do i=1,size(degenerate_symmetries)
    ! Constuct symmetry in coupled mode co-ordinates.
    symmetry = degenerate_symmetries(i)%calculate_symmetry(coupled_modes)
    
    ! Convert symmetry into unique-monomial co-ordinates.
    symmetry = mat(all_to_unique)*symmetry*mat(transpose(all_to_unique))
    
    ! Construct the projection matrix for this symmetry,
    !    and multiply the total projection matrix by this.
    symmetry_order = structure%symmetries(i)%symmetry_order()
    projection = projection * projection_matrix(symmetry,symmetry_order)
  enddo
  
  ! Transform the projection matrix into real co-ordinates.
  projection = complex_to_real_conversion &
           & * projection                 &
           & * hermitian(complex_to_real_conversion)
  
  if (sum_squares(aimag(projection))>1e-10_dp) then
    call print_line(WARNING//': Projection matrix has imaginary components.')
    call print_line( 'L2 norm of imaginary components: '// &
                   & sum_squares(aimag(projection)))
  endif
  
  ! Diagonalise the projection matrix, and check its eigenvalues.
  estuff = diagonalise_symmetric(real(projection))
  if (any(abs(estuff%eval-1)>1e-2_dp .and. abs(estuff%eval)>1e-2_dp)) then
    call err()
  endif
  
  ! Select only those eigenvectors with eigenvalue 1, and construct basis
  !    functions from them.
  estuff = estuff(filter(abs(estuff%eval-1)<1e-2_dp))
  allocate(output(size(estuff)), stat=ialloc); call err(ialloc)
  do i=1,size(estuff)
    real_coefficients = estuff(i)%evec
    complex_coefficients = cmplx( hermitian(complex_to_real_conversion) &
                              & * vec(real_coefficients))
    output(i)%real_representation = RealPolynomial( real_coefficients &
                                                & * unique_real_monomials)
    output(i)%complex_representation = ComplexPolynomial( complex_coefficients&
                                                      & * unique_monomials)
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
! H is real and symmetric, and is a projection matrix which projects out the
!    eigenvectors of U with eigenvalue 1.
! Uses Horner's rule to calculate the sum with minimal multiplication.
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
end module
