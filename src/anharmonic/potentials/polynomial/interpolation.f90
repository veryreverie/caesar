! ======================================================================
! Interpolate a polynomial to a single q-point and its pair,
!    integrating across all other q-points
! ======================================================================
module interpolation_module
  use common_module
  use anharmonic_common_module
  implicit none
  
  ! TODO: public/private
  !private
  
  type, extends(NoDefaultConstructor) :: SplitMonomial
    type(ComplexMonomial), allocatable :: head(:)
    type(ComplexMonomial), allocatable :: tail(:)
  end type
  
  interface SplitMonomial
    module procedure new_SplitMonomial_ComplexMonomial
  end interface
  
  interface size
    module procedure size_SplitMonomial
  end interface
contains

impure elemental function new_SplitMonomial_ComplexMonomial(monomial, &
   & anharmonic_data) result(this)
  implicit none
  
  type(ComplexMonomial), intent(in) :: monomial
  type(AnharmonicData),  intent(in) :: anharmonic_data
  type(SplitMonomial)               :: this
  
  integer :: total_size
  
  type(ComplexUnivariate), allocatable :: head(:)
  type(ComplexUnivariate), allocatable :: tail(:)
  
  type(ComplexMonomial) :: head_monomial
  
  type(ComplexUnivariate), allocatable :: univariates(:)
  
  integer, allocatable :: sizes(:)
  
  integer :: power,paired_power
  
  integer :: i,j,k,l,ialloc
  
  allocate( head(size(monomial)),  &
          & tail(size(monomial)),  &
          & sizes(size(monomial)), &
          & stat=ialloc); call err(ialloc)
  univariates = monomial%modes()
  do i=1,size(sizes)
    if (univariates(i)%id==univariates(i)%paired_id) then
      sizes(i) = univariates(i)%power+1
    else
      sizes(i) = (univariates(i)%power+1)*(univariates(i)%paired_power+1)
    endif
  enddo
  
  total_size = product(sizes)
  
  allocate( this%head(total_size), &
          & this%tail(total_size), &
          & stat=ialloc); call err(ialloc)
  this%head = monomial
  this%tail = monomial
  head = univariates
  tail = univariates
  l = 0
  do i=1,size(this%head)
    do j=1,size(univariates)
      k = modulo((i-1)/product(sizes(j+1:)),sizes(j))
      if (univariates(j)%id==univariates(j)%paired_id) then
        head(j)%power = k
        tail(j)%power = univariates(j)%power - head(j)%power
      else
        head(j)%power = modulo(k,univariates(j)%power+1)
        head(j)%paired_power = k/(univariates(j)%power+1)
      endif
    enddo
    head_monomial = ComplexMonomial( coefficient = monomial%coefficient, &
                                   & modes       = head                  )
    if (is_int(head_monomial%wavevector( anharmonic_data%complex_modes, &
                                       & anharmonic_data%qpoints        ))) then
      ! TODO: G-vector per subspace.
      l = l+1
      this%head(l) = ComplexMonomial( coefficient = monomial%coefficient, &
                                    & modes       = head                  )
      this%tail(l) = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                                    & modes       = tail                     )
    endif
  enddo
  this%head = this%head(:l)
  this%tail = this%tail(:l)
end function

function size_SplitMonomial(this) result(output)
  implicit none
  
  type(SplitMonomial), intent(in) :: this
  integer                         :: output
  
  output = size(this%head)
end function

function interpolate_monomial(monomial,subspaces,subspace_bases, &
   & subspace_states,modes,anharmonic_data) result(output)
  implicit none
  
  type(ComplexMonomial),    intent(in) :: monomial
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  class(SubspaceBasis),     intent(in) :: subspace_bases(:)
  class(BasisStates),       intent(in) :: subspace_states(:)
  type(ComplexMode),        intent(in) :: modes(:)
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(ComplexPolynomial)              :: output
end function

function interpolate(input) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: input
  type(ComplexPolynomial)           :: output
end function

! Generates all permutations of a list of (distinguishable) integers.
! N.B. generates duplicates if any elements of the input are the duplicated.
subroutine heaps_algorithm(input)
  implicit none
  
  integer, intent(in) :: input(:)
  
  integer, allocatable :: state(:)
  
  integer, allocatable :: indices(:)
  
  integer :: swap(2)
  
  integer :: i
  
  
  state = [(1,i=1,size(input))]
  
  indices = [(i,i=1,size(input))]
  
  call print_line(input(indices))
  
  i = 1
  do while (i<=size(input))
    if (state(i)<i) then
      if (modulo(i,2)==1) then
        swap = [1,i]
      else
        swap = [i,state(i)]
      endif
      indices(swap) = indices(swap([2,1]))
      call print_line(input(indices))
      state(i) = state(i)+1
      i = 1
    else
      state(i) = 1
      i = i+1
    endif
  enddo
end subroutine

! Generates all permutations of a list of integers, which may include repeats.
subroutine no_repeats(input)
  implicit none
  
  integer, intent(in) :: input(:)
  
  integer, allocatable :: indices(:)
  
  integer :: i,j
  
  integer :: n
  
  indices = sort(input)
  call print_line(input(indices))
  n = 1
  
  do
    ! Find the last element which has larger elements after it.
    ! Label this element j.
    j = 0
    do i=size(input)-1,1,-1
      if (input(indices(i))<input(indices(i+1))) then
        j = i
        exit
      endif
    enddo
    
    ! If there is no element j then input(indices) is now
    !    sorted in reverse order and the algorithm is done.
    if (j==0) then
      exit
    endif
    
    ! Find the last element larger than element j.
    ! Label this element i.
    ! Swap elements i and j.
    do i=size(input),j+1,-1
      if (input(indices(i))>input(indices(j))) then
        indices([i,j]) = indices([j,i])
        exit
      endif
    enddo
    
    ! Reverse the order of all elements after element j.
    indices(j+1:) = indices(size(indices):j+1:-1)
    
    !call print_line(input(indices))
    n = n+1
  enddo
  
  call print_line('n='//n)
end subroutine

! Generates all permutations of a against b.
subroutine perm(a,b)
  implicit none
  
  integer, intent(in) :: a(:)
  integer, intent(in) :: b(:)
  
  integer, allocatable :: a_id(:)
  integer, allocatable :: b_id(:)
  
  integer, allocatable :: bins(:,:)
  integer              :: no_bins
  
  integer, allocatable :: a_id_remainder(:)
  integer, allocatable :: sort_key(:)
  
  integer :: i,j,k,l,ialloc
  
  integer :: n1,n2
  
  a_id = sort(a)
  b_id = sort(b)
  
  call print_line(a(a_id))
  call print_line(b(b_id))
  call print_line('')
  
  ! Split b into bins.
  allocate(bins(2,size(b)), stat=ialloc); call err(ialloc)
  bins(1,1) = 1
  no_bins = 1
  do i=2,size(b)
    if (b(b_id(i))/=b(b_id(i-1))) then
      bins(2,no_bins) = i-1
      no_bins = no_bins+1
      bins(1,no_bins) = i
    endif
  enddo
  bins(2,no_bins) = size(b)
  bins = bins(:,:no_bins)
  
  ! Check if all elements in b are indistinguishable.
  if (no_bins==1) then
    return
  endif
  
  ! Keeps each bin of a in sorted (increasing) order.
  
  do
    ! Find the last bin which contains an element of a which is smaller than
    !    an element of a in a later bin.
    ! Label this bin j.
    do i=no_bins-1,1,-1
      if (a(a_id(bins(1,i)))<a(a_id(bins(2,i+1)))) then
        j = i
        exit
      elseif (i==1) then
        ! If there is no element j then the a(a_id) is now
        !    sorted in reverse order and the algorithm is done.
        return
      endif
    enddo
    
    ! Find the last bin containing an element of a which is larger than
    !    the smallest element in bin j.
    ! Label this bin i.
    ! Find the largest element in bin j which is smaller than any element in i.
    ! Label this element l.
    ! Find the smallest element in bin i which is larger than element l.
    ! Label this element k.
    do i=no_bins,j+1,-1
      l = last(a(a_id(bins(1,j):bins(2,j)))<a(a_id(bins(2,i))), default=0) &
      & + bins(1,j)-1
      if (l/=bins(1,j)-1) then
        k = first(a(a_id(bins(1,i):bins(2,i)))>a(a_id(l)), default=0) &
        & + bins(1,i)-1
        exit
      endif
    enddo
    
    ! Swap elements k and l.
    a_id([k,l]) = a_id([l,k])
    
    ! The elements after element l in bin j must be the smallest elements
    !    anywhere after l which are larger than element l.
    ! The elements in further bins should be sorted in ascending order.
    
    ! This is done by first reversing the order of elements within each bin
    !    after bin j, resulting in a(a_id(bins(1,j+1):)) being in descending
    !    order.
    ! Then elements a(a_id(bins(1,j+1):)) are reversed.
    ! Then element k is re-found.
    ! Then the elements after l in bin j, and the elements immediately after k
    !    are merged and split so that the smallest elements go
    !    after l in bin j, and the remainder go after k.
    
    do n1=j+1,no_bins
      a_id(bins(1,n1):bins(2,n1)) = a_id(bins(2,n1):bins(1,n1):-1)
    enddo
    
    k = bins(2,i)-(k-bins(1,i))
    
    a_id(bins(1,j+1):) = a_id(size(a_id):bins(1,j+1):-1)
    
    k = size(a_id)-(k-bins(1,j+1))
    
    n1 = bins(2,j)-l
    n2 = min(n1,size(a_id)-k)
    a_id_remainder = [a_id(l+1:l+n1), a_id(k+1:k+n2)]
    sort_key = sort(a(a_id_remainder))
    a_id(l+1:l+n1) = a_id_remainder(sort_key(:n1))
    a_id(k+1:k+n2) = a_id_remainder(sort_key(n1+1:))
    
    call print_line(a(a_id))
    call print_line(b(b_id))
    call print_line('')
  enddo
end subroutine
end module
