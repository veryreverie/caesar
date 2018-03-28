! ======================================================================
! Generates basis functions.
! ======================================================================
module basis_function_module
  use common_module
  
  use coupling_module
  use univariate_module
  use monomial_module
  use polynomial_module
  use degeneracy_module
  use degenerate_symmetry_module
  implicit none
  
  private
  
  public :: generate_basis_functions
contains

function generate_basis_functions(coupling,structure,normal_modes,qpoints, &
   & subspaces,degenerate_symmetries,vscf_basis_functions_only) result(output)
  implicit none
  
  type(CoupledSubspaces),   intent(in) :: coupling
  type(StructureData),      intent(in) :: structure
  type(ComplexMode),        intent(in) :: normal_modes(:)
  type(QpointData),         intent(in) :: qpoints(:)
  type(DegenerateModes),    intent(in) :: subspaces(:)
  type(DegenerateSymmetry), intent(in) :: degenerate_symmetries(:)
  logical,                  intent(in) :: vscf_basis_functions_only
  type(Polynomial), allocatable        :: output(:)
  
  type(DegenerateModes), allocatable :: coupled_subspaces(:)
  
  type(CoupledModes), allocatable :: coupled_modes(:)
  
  type(ComplexMatrix) :: symmetry
  type(ComplexMatrix) :: projection
  
  type(HermitianEigenstuff), allocatable :: estuff(:)
  
  integer :: i,ialloc
  
  integer :: symmetry_order
  
  type(Monomial), allocatable :: monomials(:)
  type(Monomial), allocatable :: unique_monomials(:)
  
  integer, allocatable :: equal_monomials(:)
  
  real(dp), allocatable :: all_to_unique(:,:)
  integer,  allocatable :: degeneracy(:)
  
  if (size(coupling)<2) then
    call print_line(CODE_ERROR//': Trying to generate basis functions with &
       &order less than 2.')
    call err()
  endif
  
  ! List the subspaces which are coupled together in this coupling.
  coupled_subspaces = coupling%coupled_subspaces(subspaces)
  
  ! Generate every allowed mode coupling within the coupled subspaces.
  coupled_modes = generate_mode_coupling( coupled_subspaces, &
                                        & normal_modes,      &
                                        & qpoints,           &
                                        & vscf_basis_functions_only)
  
  if (size(coupled_modes)==0) then
    output = [Polynomial::]
    return
  endif
  
  ! Convert coupled modes into Monomials with coefficient 1.
  allocate(monomials(size(coupled_modes)), stat=ialloc); call err(ialloc)
  do i=1,size(coupled_modes)
    monomials(i) = Monomial(coupled_modes(i))
  enddo
  
  ! Identify the unique monomials. (Those with all modes the same).
  unique_monomials = monomials(set(monomials,compare_monomial_modes))
  allocate( all_to_unique(size(unique_monomials),size(monomials)), &
          & degeneracy(size(unique_monomials)),                    &
          & stat=ialloc); call err(ialloc)
  
  ! Construct the mapping from all-monomial co-ordinates
  !    to unique-monomial co-ordinates.
  all_to_unique = 0
  degeneracy = 0
  do i=1,size(unique_monomials)
    equal_monomials = filter( monomials,              &
                            & compare_monomial_modes, &
                            & unique_monomials(i))
    degeneracy(i) = size(equal_monomials)
    all_to_unique(i,equal_monomials) = 1/sqrt(real(degeneracy(i),dp))
  enddo
  
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
  
  ! Diagonalise the projection matrix, and check its eigenvalues.
  estuff = diagonalise_hermitian(projection)
  if (any(abs(estuff%eval-1)>1e-2_dp .and. abs(estuff%eval)>1e-2_dp)) then
    call err()
  endif
  
  ! Select only those eigenvectors with eigenvalue 1, and construct basis
  !    functions from them.
  estuff = estuff(filter(abs(estuff%eval-1)<1e-2_dp))
  allocate(output(size(estuff)), stat=ialloc); call err(ialloc)
  do i=1,size(estuff)
    output(i) = Polynomial(estuff(i)%evec * unique_monomials)
  enddo
contains
  ! Lambda for comparing monomials.
  function compare_monomial_modes(this,that) result(output)
    implicit none
    
    class(*), intent(in) :: this
    class(*), intent(in) :: that
    logical              :: output
    
    select type(this); type is(Monomial)
      select type(that); type is(Monomial)
        if (size(this%modes)/=size(that%modes)) then
          output = .false.
        else
          output = all(this%modes==that%modes)
        endif
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
