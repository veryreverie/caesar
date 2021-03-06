submodule (caesar_polynomial_symmetry_module) caesar_polynomial_symmetry_submodule
  use caesar_polynomial_module
contains

module procedure new_BasisConversion
  logical, allocatable :: conjugate_found(:)
  
  type(ComplexMonomial) :: conjugate
  
  integer :: i,j
  
  this%monomials_ = monomials
  
  this%mapping_ = [(0, i=1, size(monomials))]
  conjugate_found = [(.false., i=1, size(monomials))]
  do i=1,size(monomials)
    if (conjugate_found(i)) then
      cycle
    endif
    
    conjugate = conjg(monomials(i))
    j = first_equivalent(monomials, conjugate, compare_complex_monomials)
    this%mapping_(i) = j
    this%mapping_(j) = i
    conjugate_found(i) = .true.
    conjugate_found(j) = .true.
  enddo
end procedure

module procedure vector_to_basis
  integer, allocatable :: conjugates(:)
  
  integer, allocatable :: a(:)
  integer, allocatable :: b(:)
  integer, allocatable :: c(:)
  
  integer :: i,ialloc
  
  if (size(input)/=size(this%mapping_)) then
    call print_line(CODE_ERROR//': Inconsistent array lengths.')
    call err()
  endif
  
  allocate(output(size(input)), stat=ialloc); call err(ialloc)
  
  conjugates = this%mapping_
  
  a = filter([(conjugates(i)==i, i=1, size(conjugates))])
  b = filter([(conjugates(i)>i,  i=1, size(conjugates))])
  c = conjugates(b)
  
  ! ui=ui', uj=uj'
  ! ui = Pij uj
  output(a) = real(input(a))
  
  ! ui/=ui', uj=uj'
  ! (ui+ui')/sqrt(2) = (Pij+Pi'j)/sqrt(2) uj
  output(b) = real(input(b)+input(c)) / sqrt(2.0_dp)
  ! (ui-ui')/sqrt(2)i = (Pij-Pi'j)/sqrt(2)i uj
  output(c) = aimag(input(b)-input(c)) / sqrt(2.0_dp)
end procedure

module procedure vector_from_basis
  real(dp) :: max_coeff
  
  integer, allocatable :: included_terms(:)
  
  type(ComplexMonomial), allocatable :: terms(:)
  
  integer :: i,j,k,ialloc
  
  max_coeff = maxval(abs(vector))
  
  included_terms = filter( abs(vector)>=max_coeff*threshold                &
                    & .or. abs(vector(this%mapping_))>=max_coeff*threshold )
  
  allocate(terms(size(included_terms)), stat=ialloc); call err(ialloc)
  do i=1,size(terms)
    j = included_terms(i)
    k = this%mapping_(j)
    terms(i) = this%monomials_(j)
    if (j==k) then
      ! The monomial is its own pair.
      terms(i)%coefficient = vector(j)
    elseif (j<k) then
      ! The monomial looks like u+ = uc + i*us.
      terms(i)%coefficient = cmplx(vector(j),vector(k),dp)
    else
      ! The monomial looks like u- = uc - i*us.
      terms(i)%coefficient = cmplx(vector(k),-vector(j),dp)
    endif
  enddo
  
  output = ComplexPolynomial(terms)
end procedure

module procedure matrix_to_basis
  if (lhs .and. rhs) then
    output = mat(project_matrix_lhs_rhs(cmplx(matrix), this%mapping_))
  elseif (lhs) then
    output = mat(project_matrix_lhs(cmplx(matrix), this%mapping_))
  elseif (rhs) then
    output = mat(project_matrix_rhs(cmplx(matrix), this%mapping_))
  else
    call err()
  endif
end procedure

module procedure find_monomial_conjugates
  logical, allocatable :: conjugate_found(:)
  
  type(ComplexMonomial) :: conjugate
  
  integer :: i,j
  
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
end procedure

module procedure project_matrix_lhs_rhs
  integer, allocatable :: a(:)
  integer, allocatable :: b(:)
  integer, allocatable :: c(:)
  
  integer :: i,ialloc
  
  ! Construct the output matrix.
  allocate( output(size(input,1),size(input,2)), &
          & stat=ialloc); call err(ialloc)
  
  a = filter([(conjugates(i)==i, i=1, size(conjugates))])
  b = filter([(conjugates(i)>i,  i=1, size(conjugates))])
  c = conjugates(b)
  
  ! ui=ui', uj=uj'
  ! ui = Pij uj
  output(a,a) = real(input(a,a))
  
  ! ui=ui', uj/=uj'
  ! ui = (Pij+Pij')/sqrt(2) (uj+uj')/sqrt(2) + ...
  output(a,b) = real(input(a,b)+input(a,c)) / sqrt(2.0_dp)
  ! ... + i(Pij-Pij')/sqrt(2) (uj-uj')/sqrt(2)i
  output(a,c) = -aimag(input(a,b)-input(a,c)) / sqrt(2.0_dp)
  
  ! ui/=ui', uj=uj'
  ! (ui+ui')/sqrt(2) = (Pij+Pi'j)/sqrt(2) uj
  output(b,a) = real(input(b,a)+input(c,a)) / sqrt(2.0_dp)
  ! (ui-ui')/sqrt(2)i = (Pij-Pi'j)/sqrt(2)i uj
  output(c,a) = aimag(input(b,a)-input(c,a)) / sqrt(2.0_dp)
  
  ! ui/=ui', uj/=uj'
  ! (ui+ui')/sqrt(2) = (Pij+Pij'+Pi'j+Pi'j')/2 (uj+uj')/sqrt(2) + ...
  output(b,b) = real(input(b,b)+input(b,c)+input(c,b)+input(c,c)) / 2
  ! ... + i(Pij-Pij'+Pi'j-Pi'j')/2 (uj-uj')/sqrt(2)i
  output(b,c) = -aimag(input(b,b)-input(b,c)+input(c,b)-input(c,c)) / 2
  ! (ui-ui')/sqrt(2)i = (Pij+Pij'-Pi'j-Pi'j')/2i (uj+uj')/sqrt(2) + ...
  output(c,b) = aimag(input(b,b)+input(b,c)-input(c,b)-input(c,c)) / 2
  ! ... +  (Pij-Pij'-Pi'j+Pi'j')/2 (uj-uj')/sqrt(2)i
  output(c,c) = real(input(b,b)-input(b,c)-input(c,b)+input(c,c)) / 2
end procedure

module procedure project_matrix_lhs
  integer, allocatable :: a(:)
  integer, allocatable :: b(:)
  integer, allocatable :: c(:)
  
  integer :: i,ialloc
  
  ! Construct the output matrix.
  allocate( output(size(input,1),size(input,2)), &
          & stat=ialloc); call err(ialloc)
  
  a = filter([(conjugates(i)==i, i=1, size(conjugates))])
  b = filter([(conjugates(i)>i,  i=1, size(conjugates))])
  c = conjugates(b)
  
  ! ui=ui', uj=uj'
  ! ui = Pij uj
  output(a,:) = real(input(a,:))
  
  ! ui/=ui', uj=uj'
  ! (ui+ui')/sqrt(2) = (Pij+Pi'j)/sqrt(2) uj
  output(b,:) = real(input(b,:)+input(c,:)) / sqrt(2.0_dp)
  ! (ui-ui')/sqrt(2)i = (Pij-Pi'j)/sqrt(2)i uj
  output(c,:) = aimag(input(b,:)-input(c,:)) / sqrt(2.0_dp)
end procedure

module procedure project_matrix_rhs
  integer, allocatable :: a(:)
  integer, allocatable :: b(:)
  integer, allocatable :: c(:)
  
  integer :: i,ialloc
  
  ! Construct the output matrix.
  allocate( output(size(input,1),size(input,2)), &
          & stat=ialloc); call err(ialloc)
  
  a = filter([(conjugates(i)==i, i=1, size(conjugates))])
  b = filter([(conjugates(i)>i,  i=1, size(conjugates))])
  c = conjugates(b)
  
  ! ui=ui', uj=uj'
  ! ui = Pij uj
  output(:,a) = real(input(:,a))
  
  ! ui=ui', uj/=uj'
  ! ui = (Pij+Pij')/sqrt(2) (uj+uj')/sqrt(2) + ...
  output(:,b) = real(input(:,b)+input(:,c)) / sqrt(2.0_dp)
  ! ... + i(Pij-Pij')/sqrt(2) (uj-uj')/sqrt(2)i
  output(:,c) = -aimag(input(:,b)-input(:,c)) / sqrt(2.0_dp)
end procedure
end submodule
