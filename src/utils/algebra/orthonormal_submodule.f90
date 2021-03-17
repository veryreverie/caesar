submodule (caesar_orthonormal_module) caesar_orthonormal_submodule
  use caesar_algebra_module
contains

module procedure orthonormal_basis_RealVectors
  ! Working variables for diagonalisation.
  type(RealMatrix)                       :: sum_outer_products
  type(SymmetricEigenstuff), allocatable :: eigenstuff(:)
  
  integer :: vector_size
  integer :: i,ialloc
  
  ! Return early if nothing to process.
  if (size(input)==0) then
    allocate(output(0), stat=ialloc); call err(ialloc)
    return
  endif
  
  ! Check vector dimensionalities are consistent, and record dimensionality.
  vector_size = size(input(1))
  do i=2,size(input)
    if (size(input(i))/=vector_size) then
      call print_line(CODE_ERROR//': Trying to orhonormalise vectors of &
         &inconsistent lengths.')
      call err()
    endif
  enddo
  
  ! Construct the sum of the outer products of the input vectors.
  sum_outer_products = dblemat(zeroes(vector_size,vector_size))
  do i=1,size(input)
    sum_outer_products = sum_outer_products &
                     & + outer_product(input(i),input(i))
  enddo
  eigenstuff = diagonalise_symmetric(sum_outer_products)
  
  ! Select only the vectors {u} with sum(u.x)>longest_invalid.
  eigenstuff = eigenstuff(filter(eigenstuff%eval>longest_invalid))
  
  ! Check no vector has sum(u.x)<shortest_valid.
  if (any(eigenstuff%eval<shortest_valid)) then
    call print_line(ERROR//': orthonormalise produced a basis &
       &vector with a projection between longest_invalid and &
       &shortest_valid.')
    call err()
  endif
  
  ! Transfer basis vectors to output.
  allocate(output(size(eigenstuff)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    output(i) = vec(eigenstuff(i)%evec)
  enddo
end procedure

module procedure orthonormal_basis_ComplexVectors
  ! Working variables for diagonalisation.
  type(ComplexMatrix)                    :: sum_outer_products
  type(HermitianEigenstuff), allocatable :: eigenstuff(:)
  
  integer :: vector_size
  integer :: i,ialloc
  
  ! Return early if nothing to process.
  if (size(input)==0) then
    allocate(output(0), stat=ialloc); call err(ialloc)
    return
  endif
  
  ! Check vector dimensionalities are consistent, and record dimensionality.
  vector_size = size(input(1))
  do i=2,size(input)
    if (size(input(i))/=vector_size) then
      call print_line(CODE_ERROR//': Trying to orhonormalise vectors of &
         &inconsistent lengths.')
      call err()
    endif
  enddo
  
  ! Construct the sum of the outer products of the input vectors.
  sum_outer_products = cmplxmat(zeroes(vector_size,vector_size))
  do i=1,size(input)
    sum_outer_products = sum_outer_products &
                     & + outer_product(input(i),conjg(input(i)))
  enddo
  eigenstuff = diagonalise_hermitian(sum_outer_products)
  
  ! Select only the vectors {u} with sum(u.x)>longest_invalid.
  eigenstuff = eigenstuff(filter(eigenstuff%eval>longest_invalid))
  
  ! Check no vector has sum(u.x)<shortest_valid.
  if (any(eigenstuff%eval<shortest_valid)) then
    call print_line(ERROR//': orthonormalise produced a basis &
       &vector with a projection between longest_invalid and &
       &shortest_valid.')
    call err()
  endif
  
  ! Transfer basis vectors to output.
  allocate(output(size(eigenstuff)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    output(i) = vec(eigenstuff(i)%evec)
  enddo
end procedure

module procedure intersection_basis_ComplexVectors
  ! Working variables for QR decomposition.
  complex(dp), allocatable     :: matrix(:,:)
  type(ComplexQRDecomposition) :: qr
  
  ! The basis vectors, and their projections onto the input vectors.
  type(ComplexVector), allocatable :: left_vectors(:)
  real(dp),            allocatable :: left_projections(:)
  type(ComplexVector), allocatable :: right_vectors(:)
  real(dp),            allocatable :: right_projections(:)
  
  integer :: vector_length
  integer :: union_size
  integer :: i,ialloc
  
  ! Return early if nothing to process.
  if (size(a)==0 .or. size(b)==0) then
    allocate(output(0), stat=ialloc); call err(ialloc)
    return
  endif
  
  ! Check vector dimensionalities are consistent, and record said lengths.
  vector_length = size(a(1))
  do i=2,size(a)
    if (size(a(i))/=vector_length) then
      call print_line(CODE_ERROR//': Trying to orhonormalise vectors of &
         &inconsistent lengths.')
      call err()
    endif
  enddo
  do i=1,size(b)
    if (size(b(i))/=vector_length) then
      call print_line(CODE_ERROR//': Trying to orhonormalise vectors of &
         &inconsistent lengths.')
      call err()
    endif
  enddo
  
  ! Arrange input vectors into matrix and perform QR decomposition.
  ! The Zassenhaus matrix is given by Z = (A A)
  !                                       (B 0)
  ! Where the rows of A and B are or vectors of a and b respectively.
  !
  ! After the QR factorisation of Z, R is given by R = (U *)
  !                                                    (0 I)
  ! Where:
  !  The rows of U are orthogonal vectors spanning the union of a and b.
  !  The rows of I are orthogonal vectors spanning the intersection of a and b.
  allocate( matrix(size(a)+size(b),2*vector_length), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(a)
    matrix(i,                :vector_length) = cmplx(a(i))
    matrix(i, vector_length+1:             ) = cmplx(a(i))
  enddo
  do i=1,size(b)
    matrix(size(a)+i,                :vector_length) = cmplx(b(i))
    matrix(size(a)+i, vector_length+1:             ) = cmplx(0.0_dp,0.0_dp,dp)
  enddo
  qr = qr_decomposition(matrix)
  
  ! Extract the vectors in [U,0] and [*,I], and calculate their lengths.
  allocate( left_vectors(size(a)+size(b)),      &
          & left_projections(size(a)+size(b)),  &
          & right_vectors(size(a)+size(b)),     &
          & right_projections(size(a)+size(b)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(a)+size(b)
    left_vectors(i) = vec(qr%r(i,:vector_length))
    left_projections(i) = l2_norm(left_vectors(i))
    
    right_vectors(i) = vec(qr%r(i,vector_length+1:))
    right_projections(i) = l2_norm(right_vectors(i))
  enddo
  
  ! Calculate the size of U by finding the first row of 0.
  union_size = first(left_projections<0.1_dp,default=size(left_projections)+1)&
           & - 1
  
  ! Check that the U is larger than I, that the rows of 0 are zero,
  !    and that the rows of I are not small.
  if (union_size<size(matrix,1)-union_size) then
    call print_line(CODE_ERROR//': the union is smaller than the &
       &intersection. Check that a and b are orthonormal.')
    call err()
  elseif (any(left_projections(union_size+1:)>1e-10_dp)) then
    call print_line(CODE_ERROR//': non-zero vector below union in reduced &
       &matrix. Check a and b are orthonormal.')
    call err()
  elseif (any(right_projections(union_size+1:)<0.1_dp)) then
    call print_line(CODE_ERROR//': small vector in intersection. Check a and &
       &b are orthonormal.')
    call err()
  endif
  
  ! Construct the output, the normalised rows of I.
  output = right_vectors(union_size+1:) / right_projections(union_size+1:)
end procedure
end submodule