submodule (caesar_fraction_algebra_module) caesar_fraction_algebra_submodule
  use caesar_algebra_module
contains

module procedure vec_IntFractions
  output%contents_ = input
end procedure

module procedure mat_IntFractions
  output%contents_ = input
end procedure

module procedure mat_IntFractions_shape
  type(IntFraction), allocatable :: contents(:,:)
  
  !output%contents_ = transpose(reshape(input, [m,n]))
  ! WORKAROUND to avoid internal compiler error in ifort 19.1.0.166.
  contents = reshape(input, [m,n])
  output = mat(transpose(contents))
end procedure

module procedure frac_FractionVector
  if (allocated(input%contents_)) then
    output = input%contents_
  else
    call print_line(CODE_ERROR//': Trying to use the contents of a vector &
       &before it has been allocated.')
    call err()
  endif
end procedure

module procedure frac_FractionMatrix
  if (allocated(input%contents_)) then
    output = input%contents_
  else
    call print_line(CODE_ERROR//': Trying to use the contents of a matrix &
       &before it has been allocated.')
    call err()
  endif
end procedure

module procedure fracvec_IntVector
  output = vec(frac(int(input)))
end procedure

module procedure fracmat_IntMatrix
  output = mat(frac(int(input)))
end procedure

module procedure intvec_FractionVector
  output = vec(int(frac(input)))
end procedure

module procedure intmat_FractionMatrix
  output = mat(int(frac(input)))
end procedure

module procedure dblevec_FractionVector
  output = vec(dble(frac(input)))
end procedure

module procedure dblemat_FractionMatrix
  output = mat(dble(frac(input)))
end procedure

module procedure size_FractionVector
  output = size(frac(this))
end procedure

module procedure size_FractionMatrix
  output = size(frac(this), dim)
end procedure

module procedure is_int_FractionVector
  output = all(is_int(frac(this)))
end procedure

module procedure is_int_FractionMatrix
  output = all(is_int(frac(this)))
end procedure

module procedure equality_FractionVector_FractionVector
  output = all(frac(this)==frac(that))
end procedure

module procedure equality_FractionVector_IntVector
  output = all(frac(this)==int(that))
end procedure

module procedure equality_IntVector_FractionVector
  output = all(int(this)==frac(that))
end procedure

module procedure equality_FractionMatrix_FractionMatrix
  output = all(frac(this)==frac(that))
end procedure

module procedure equality_FractionMatrix_IntMatrix
  output = all(frac(this)==int(that))
end procedure

module procedure equality_IntMatrix_FractionMatrix
  output = all(int(this)==frac(that))
end procedure

module procedure non_equality_FractionVector_FractionVector
  output = .not. this==that
end procedure

module procedure non_equality_FractionVector_IntVector
  output = .not. this==that
end procedure

module procedure non_equality_IntVector_FractionVector
  output = .not. this==that
end procedure

module procedure non_equality_FractionMatrix_FractionMatrix
  output = .not. this==that
end procedure

module procedure non_equality_FractionMatrix_IntMatrix
  output = .not. this==that
end procedure

module procedure non_equality_IntMatrix_FractionMatrix
  output = .not. this==that
end procedure

module procedure transpose_FractionMatrix
  output = mat(transpose(frac(this)))
end procedure

module procedure invert_integers
  integer :: det
  integer :: i,j
  
  ! Check that the input is a 3x3 matrix.
  if (size(input,1)/=3 .or. size(input,2)/=3) then
    call print_line(CODE_ERROR//': Trying to invert matrix which is not 3x3.')
    call err()
  endif
  
  det = determinant(input)
  if (det==0) then
    call print_line(ERROR//': Trying to invert a Matrix with determinant=0.')
    call err()
  endif
 
  output(1,1) = frac(input(2,2)*input(3,3)-input(3,2)*input(2,3))
  output(1,2) = frac(input(3,2)*input(1,3)-input(1,2)*input(3,3))
  output(1,3) = frac(input(1,2)*input(2,3)-input(2,2)*input(1,3))
  output(2,1) = frac(input(2,3)*input(3,1)-input(3,3)*input(2,1))
  output(2,2) = frac(input(3,3)*input(1,1)-input(1,3)*input(3,1))
  output(2,3) = frac(input(1,3)*input(2,1)-input(2,3)*input(1,1))
  output(3,1) = frac(input(2,1)*input(3,2)-input(3,1)*input(2,2))
  output(3,2) = frac(input(3,1)*input(1,2)-input(1,1)*input(3,2))
  output(3,3) = frac(input(1,1)*input(2,2)-input(2,1)*input(1,2))
  
  do i=1,3
    do j=1,3
      output(j,i) = output(j,i)/det
    enddo
  enddo
end procedure

module procedure invert_IntMatrix
  output = mat(invert(int(input)))
end procedure

module procedure add_FractionVector_FractionVector
  output = vec(frac(this) + frac(that))
end procedure

module procedure add_FractionVector_IntVector
  output = vec(int(this) + frac(that))
end procedure

module procedure add_IntVector_FractionVector
  output = vec(frac(this) + int(that))
end procedure

module procedure add_FractionMatrix_FractionMatrix
  output = mat(frac(this) + frac(that))
end procedure

module procedure add_FractionMatrix_IntMatrix
  output = mat(int(this) + frac(that))
end procedure

module procedure add_IntMatrix_FractionMatrix
  output = mat(frac(this) + int(that))
end procedure

module procedure negative_FractionVector
  output = vec(-frac(this))
end procedure

module procedure negative_FractionMatrix
  output = mat(-frac(this))
end procedure

module procedure subtract_FractionVector_FractionVector
  output = vec(frac(this) - frac(that))
end procedure

module procedure subtract_FractionVector_IntVector
  output = vec(frac(this) - int(that))
end procedure

module procedure subtract_IntVector_FractionVector
  output = vec(int(this) - frac(that))
end procedure

module procedure subtract_FractionMatrix_FractionMatrix
  output = mat(frac(this) - frac(that))
end procedure

module procedure subtract_FractionMatrix_IntMatrix
  output = mat(frac(this) - int(that))
end procedure

module procedure subtract_IntMatrix_FractionMatrix
  output = mat(int(this) - frac(that))
end procedure

module procedure multiply_FractionVector_integer
  output = vec(frac(this) * that)
end procedure

module procedure multiply_integer_FractionVector
  output = vec(this * frac(that))
end procedure

module procedure multiply_FractionVector_IntFraction
  output = vec(frac(this) * that)
end procedure

module procedure multiply_IntFraction_FractionVector
  output = vec(this * frac(that))
end procedure

module procedure multiply_FractionMatrix_integer
  output = mat(frac(this) * that)
end procedure

module procedure multiply_integer_FractionMatrix
  output = mat(this * frac(that))
end procedure

module procedure multiply_FractionMatrix_IntFraction
  output = mat(frac(this) * that)
end procedure

module procedure multiply_IntFraction_FractionMatrix
  output = mat(this * frac(that))
end procedure

module procedure dot_FractionVector_FractionVector
  type(IntFraction), allocatable :: a(:)
  type(IntFraction), allocatable :: b(:)
  
  integer :: i
  
  if (size(this)/=size(that)) then
    call print_line( CODE_ERROR//': Dot product of vectors of different &
       &lengths.')
    call err()
  endif
  
  a = frac(this)
  b = frac(that)
  output = frac(0)
  do i=1,size(this)
    output = output + a(i)*b(i)
  enddo
end procedure

module procedure dot_FractionVector_IntVector
  type(IntFraction), allocatable :: a(:)
  integer,           allocatable :: b(:)
  
  integer :: i
  
  if (size(this)/=size(that)) then
    call print_line( CODE_ERROR//': Dot product of vectors of different &
       &lengths.')
    call err()
  endif
  
  a = frac(this)
  b = int(that)
  output = frac(0)
  do i=1,size(this)
    output = output + a(i)*b(i)
  enddo
end procedure

module procedure dot_IntVector_FractionVector
  integer,           allocatable :: a(:)
  type(IntFraction), allocatable :: b(:)
  
  integer :: i
  
  if (size(this)/=size(that)) then
    call print_line( CODE_ERROR//': Dot product of vectors of different &
       &lengths.')
    call err()
  endif
  
  a = int(this)
  b = frac(that)
  output = frac(0)
  do i=1,size(this)
    output = output + a(i)*b(i)
  enddo
end procedure

module procedure dot_FractionMatrix_FractionVector
  type(IntFraction), allocatable :: a(:,:)
  type(IntFraction), allocatable :: b(:)
  type(IntFraction), allocatable :: contents(:)
  
  integer :: i,ialloc
  
  if (size(this,2)/=size(that)) then
    call print_line(CODE_ERROR//': Dot product of matrix and vector of &
       &incompatible dimensions.')
    call err()
  endif
  
  a = frac(this)
  b = frac(that)
  allocate(contents(size(this,1)), stat=ialloc); call err(ialloc)
  contents = frac(0)
  do i=1,size(that)
    contents = contents + a(:,i)*b(i)
  enddo
  output = vec(contents)
end procedure

module procedure dot_FractionMatrix_IntVector
  type(IntFraction), allocatable :: a(:,:)
  integer,           allocatable :: b(:)
  type(IntFraction), allocatable :: contents(:)
  
  integer :: i,ialloc
  
  if (size(this,2)/=size(that)) then
    call print_line(CODE_ERROR//': Dot product of matrix and vector of &
       &incompatible dimensions.')
    call err()
  endif
  
  a = frac(this)
  b = int(that)
  allocate(contents(size(this,1)), stat=ialloc); call err(ialloc)
  contents = frac(0)
  do i=1,size(that)
    contents = contents + a(:,i)*b(i)
  enddo
  output = vec(contents)
end procedure

module procedure dot_IntMatrix_FractionVector
  integer,           allocatable :: a(:,:)
  type(IntFraction), allocatable :: b(:)
  type(IntFraction), allocatable :: contents(:)
  
  integer :: i,ialloc
  
  if (size(this,2)/=size(that)) then
    call print_line(CODE_ERROR//': Dot product of matrix and vector of &
       &incompatible dimensions.')
    call err()
  endif
  
  a = int(this)
  b = frac(that)
  allocate(contents(size(this,1)), stat=ialloc); call err(ialloc)
  contents = frac(0)
  do i=1,size(that)
    contents = contents + a(:,i)*b(i)
  enddo
  output = vec(contents)
end procedure

module procedure dot_FractionVector_FractionMatrix
  type(IntFraction), allocatable :: a(:)
  type(IntFraction), allocatable :: b(:,:)
  type(IntFraction), allocatable :: contents(:)
  
  integer :: i,ialloc
  
  if (size(this)/=size(that,1)) then
    call print_line(CODE_ERROR//': Dot product of vector and matrix of &
       &incompatible dimensions.')
    call err()
  endif
  
  a = frac(this)
  b = frac(that)
  allocate(contents(size(that,2)), stat=ialloc); call err(ialloc)
  contents = frac(0)
  do i=1,size(this)
    contents = contents + a(i)*b(i,:)
  enddo
  output = vec(contents)
end procedure

module procedure dot_FractionVector_IntMatrix
  type(IntFraction), allocatable :: a(:)
  integer,           allocatable :: b(:,:)
  type(IntFraction), allocatable :: contents(:)
  
  integer :: i,ialloc
  
  if (size(this)/=size(that,1)) then
    call print_line(CODE_ERROR//': Dot product of vector and matrix of &
       &incompatible dimensions.')
    call err()
  endif
  
  a = frac(this)
  b = int(that)
  allocate(contents(size(that,2)), stat=ialloc); call err(ialloc)
  contents = frac(0)
  do i=1,size(this)
    contents = contents + a(i)*b(i,:)
  enddo
  output = vec(contents)
end procedure

module procedure dot_IntVector_FractionMatrix
  integer,           allocatable :: a(:)
  type(IntFraction), allocatable :: b(:,:)
  type(IntFraction), allocatable :: contents(:)
  
  integer :: i,ialloc
  
  if (size(this)/=size(that,1)) then
    call print_line(CODE_ERROR//': Dot product of vector and matrix of &
       &incompatible dimensions.')
    call err()
  endif
  
  a = int(this)
  b = frac(that)
  allocate(contents(size(that,2)), stat=ialloc); call err(ialloc)
  contents = frac(0)
  do i=1,size(this)
    contents = contents + a(i)*b(i,:)
  enddo
  output = vec(contents)
end procedure

module procedure dot_FractionMatrix_FractionMatrix
  type(IntFraction), allocatable :: a(:,:)
  type(IntFraction), allocatable :: b(:,:)
  type(IntFraction), allocatable :: contents(:,:)
  
  integer :: i,j,k,ialloc
  
  if (size(this,2)/=size(that,1)) then
    call print_line(CODE_ERROR//': Dot product of two matrices of &
       &incompatible dimensions.')
    call err()
  endif
  
  a = frac(this)
  b = frac(that)
  allocate(contents(size(this,1),size(that,2)), stat=ialloc); call err(ialloc)
  contents = frac(0)
  do i=1,size(that,2)
    do j=1,size(that,1)
      do k=1,size(this,1)
        contents(k,i) = contents(k,i) + a(k,j)*b(j,i)
      enddo
    enddo
  enddo
  output = mat(contents)
end procedure

module procedure dot_FractionMatrix_IntMatrix
  type(IntFraction), allocatable :: a(:,:)
  integer,           allocatable :: b(:,:)
  type(IntFraction), allocatable :: contents(:,:)
  
  integer :: i,j,k,ialloc
  
  if (size(this,2)/=size(that,1)) then
    call print_line(CODE_ERROR//': Dot product of two matrices of &
       &incompatible dimensions.')
    call err()
  endif
  
  a = frac(this)
  b = int(that)
  allocate(contents(size(this,1),size(that,2)), stat=ialloc); call err(ialloc)
  contents = frac(0)
  do i=1,size(that,2)
    do j=1,size(that,1)
      do k=1,size(this,1)
        contents(k,i) = contents(k,i) + a(k,j)*b(j,i)
      enddo
    enddo
  enddo
  output = mat(contents)
end procedure

module procedure dot_IntMatrix_FractionMatrix
  integer,           allocatable :: a(:,:)
  type(IntFraction), allocatable :: b(:,:)
  type(IntFraction), allocatable :: contents(:,:)
  
  integer :: i,j,k,ialloc
  
  if (size(this,2)/=size(that,1)) then
    call print_line(CODE_ERROR//': Dot product of two matrices of &
       &incompatible dimensions.')
    call err()
  endif
  
  a = int(this)
  b = frac(that)
  allocate(contents(size(this,1),size(that,2)), stat=ialloc); call err(ialloc)
  contents = frac(0)
  do i=1,size(that,2)
    do j=1,size(that,1)
      do k=1,size(this,1)
        contents(k,i) = contents(k,i) + a(k,j)*b(j,i)
      enddo
    enddo
  enddo
  output = mat(contents)
end procedure

module procedure divide_FractionVector_integer
  output = vec(frac(this) / that)
end procedure

module procedure divide_FractionVector_IntFraction
  output = vec(frac(this) / that)
end procedure

module procedure divide_FractionMatrix_integer
  output = mat(frac(this) / that)
end procedure

module procedure divide_FractionMatrix_IntFraction
  output = mat(frac(this) / that)
end procedure

module procedure sum_FractionVector
  integer :: i
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': Trying to take the sum of an empty array.')
    call err()
  endif
  
  output = input(1)
  do i=2,size(input)
    output = output+input(i)
  enddo
end procedure

module procedure sum_FractionMatrix
  integer :: i
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': Trying to take the sum of an empty array.')
    call err()
  endif
  
  output = input(1)
  do i=2,size(input)
    output = output+input(i)
  enddo
end procedure

module procedure exp_2pii_IntFraction
  output = exp_2pii(dble(input))
end procedure

module procedure cos_2pi_IntFraction
  output = cos_2pi(dble(input))
end procedure

module procedure sin_2pi_IntFraction
  output = sin_2pi(dble(input))
end procedure

module procedure read_FractionVector
  select type(this); type is(FractionVector)
    this = vec(frac(split_line(input)))
  end select
end procedure

module procedure write_FractionVector
  select type(this); type is(FractionVector)
    output = join(frac(this))
  end select
end procedure

module procedure new_FractionVector_String
  call this%read(input)
end procedure

module procedure read_FractionMatrix
  type(IntFraction), allocatable :: line(:)
  type(IntFraction), allocatable :: contents(:,:)
  
  integer :: i,ialloc
  
  select type(this); type is(FractionMatrix)
    if (size(input)==0) then
      allocate(contents(0,0), stat=ialloc); call err(ialloc)
    else
      line = frac(split_line(input(1)))
      allocate( contents(size(input),size(line)), &
              & stat=ialloc); call err(ialloc)
      contents(1,:) = line
      do i=2,size(input)
        line = frac(split_line(input(i)))
        if (size(line)/=size(contents,2)) then
          call print_line(ERROR//': Reading matrix: rows of different &
             &lengths.')
          call err()
        endif
        contents(i,:) = line
      enddo
    endif
    
    this = mat(contents)
  class default
    call err()
  end select
end procedure

module procedure write_FractionMatrix
  type(IntFraction), allocatable :: contents(:,:)
  type(IntFraction), allocatable :: row(:)
  
  integer :: i,ialloc
  
  select type(this); type is(FractionMatrix)
    contents = frac(this)
    allocate(output(size(this,1)), stat=ialloc); call err(ialloc)
    do i=1,size(this,1)
      row = contents(i,:)
      output(i) = join(row)
    enddo
  class default
    call err()
  end select
end procedure

module procedure new_FractionMatrix_Strings
  call this%read(input)
end procedure

module procedure new_FractionMatrix_StringArray
  this = FractionMatrix(str(input))
end procedure
end submodule
