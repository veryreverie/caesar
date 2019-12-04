! ======================================================================
! Generates and stores the permutation of one list against another.
! ======================================================================
! Takes two lists of integers, a and b, which may or may not contain
!    duplicate elements.
! Generates all unique permutations of a against b.
! See permutation_example for an example of how to use this class.
module permutation_module
  use common_module
  implicit none
  
  private
  
  public :: PermutationData
  
  type, extends(NoDefaultConstructor) :: PermutationData
    integer, allocatable, private :: a_(:)
    integer, allocatable, private :: b_(:)
    integer, allocatable, private :: bins_(:,:)
    integer, allocatable, private :: a_permutation_(:)
    integer, allocatable, private :: b_permutation_(:)
    logical,              private :: all_permutations_done_
  contains
    procedure, public :: next_permutation => &
                       & next_permutation_PermutationData
    procedure, public :: a => &
                       & a_PermutationData
    procedure, public :: b => &
                       & b_PermutationData
    procedure, public :: all_permutations_done => &
                       & all_permutations_done_PermutationData
  end type
  
  interface PermutationData
    module procedure new_PermutationData
  end interface
contains

function new_PermutationData(a,b) result(this)
  implicit none
  
  integer, intent(in)   :: a(:)
  integer, intent(in)   :: b(:)
  type(PermutationData) :: this
  
  integer :: no_bins
  
  integer :: i,ialloc
  
  this%a_ = a
  this%b_ = b
  
  this%a_permutation_ = sort(a)
  this%b_permutation_ = sort(b)
  
  this%all_permutations_done_ = .false.
  
  ! Split b into bins.
  allocate(this%bins_(2,size(b)), stat=ialloc); call err(ialloc)
  this%bins_(1,1) = 1
  no_bins = 1
  do i=2,size(b)
    if (b(this%b_permutation_(i))/=b(this%b_permutation_(i-1))) then
      this%bins_(2,no_bins) = i-1
      no_bins = no_bins+1
      this%bins_(1,no_bins) = i
    endif
  enddo
  this%bins_(2,no_bins) = size(b)
  this%bins_ = this%bins_(:,:no_bins)
end function

! Generates the next permutation of a against b.
! Keeps each bin of a sorted in ascending order.
! Does not permute b.
subroutine next_permutation_PermutationData(this)
  implicit none
  
  class(PermutationData), intent(inout) :: this
  
  integer, allocatable :: permutation_remainder(:)
  integer, allocatable :: sort_key(:)
  
  integer :: i,j,k,l,ialloc
  
  integer :: n1,n2
  
  if (this%all_permutations_done_) then
    call err()
  endif
  
  associate(bins=>this%bins_, a=>this%a_, permutation=>this%a_permutation_)
    ! Find the last bin which contains an element of a which is smaller than
    !    an element of a in a later bin.
    ! Label this bin j.
    do i=size(bins,2)-1,1,-1
      if (a(permutation(bins(1,i)))<a(permutation(bins(2,i+1)))) then
        j = i
        exit
      elseif (i==1) then
        ! If there is no element j then the a(permutation) is now
        !    sorted in reverse order and the algorithm is done.
        this%all_permutations_done_ = .true.
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
    do i=size(bins,2),j+1,-1
      l = last( a(permutation(bins(1,j):bins(2,j)))   &
      &       < a(permutation(bins(2,i))),            &
      &         default=0                           ) &
      & + bins(1,j)-1
      if (l/=bins(1,j)-1) then
        k = first( a(permutation(bins(1,i):bins(2,i)))   &
        &        > a(permutation(l))                   ) &
        & + bins(1,i)-1
        exit
      endif
    enddo
    
    ! Swap elements k and l.
    permutation([k,l]) = permutation([l,k])
    
    ! The elements after element l in bin j must be the smallest elements
    !    anywhere after l which are larger than element l.
    ! The elements in further bins should be sorted in ascending order.
    
    ! This is done by first reversing the order of elements within each bin
    !    after bin j, resulting in a(permutation(bins(1,j+1):)) being in
    !    descending order.
    ! Then elements a(permutation(bins(1,j+1):)) are reversed.
    ! Then element k is re-found.
    ! Then the elements after l in bin j, and the elements immediately after k
    !    are merged and split so that the smallest elements go
    !    after l in bin j, and the remainder go after k.
    
    do n1=j+1,size(bins,2)
      permutation(bins(1,n1):bins(2,n1)) = &
         & permutation(bins(2,n1):bins(1,n1):-1)
    enddo
    
    k = bins(2,i)-(k-bins(1,i))
    
    permutation(bins(1,j+1):) = permutation(size(permutation):bins(1,j+1):-1)
    
    k = size(permutation)-(k-bins(1,j+1))
    
    n1 = bins(2,j)-l
    n2 = min(n1,size(permutation)-k)
    permutation_remainder = [permutation(l+1:l+n1), permutation(k+1:k+n2)]
    sort_key = sort(a(permutation_remainder))
    permutation(l+1:l+n1) = permutation_remainder(sort_key(:n1))
    permutation(k+1:k+n2) = permutation_remainder(sort_key(n1+1:))
  end associate
end subroutine

function a_PermutationData(this) result(output)
  implicit none
  
  class(PermutationData), intent(in) :: this
  integer, allocatable               :: output(:)
  
  if (this%all_permutations_done_) then
    call err()
  endif
  
  output = this%a_(this%a_permutation_)
end function

function b_PermutationData(this) result(output)
  implicit none
  
  class(PermutationData), intent(in) :: this
  integer, allocatable               :: output(:)
  
  if (this%all_permutations_done_) then
    call err()
  endif
  
  output = this%b_(this%b_permutation_)
end function

function all_permutations_done_PermutationData(this) result(output)
  implicit none
  
  class(PermutationData), intent(in) :: this
  logical                            :: output
  
  output = this%all_permutations_done_
end function

! --------------------------------------------------
! Usage example.
! --------------------------------------------------
subroutine permutation_example()
  implicit none
  
  integer, allocatable :: a(:)
  integer, allocatable :: b(:)
  
  type(PermutationData) :: permutation
  
  a = [1,2,3,3,4]
  b = [1,2,2,3,3]
  
  permutation = PermutationData(a,b)
  
  do
    call print_line('')
    call print_line(permutation%a())
    call print_line(permutation%b())
    
    call permutation%next_permutation()
    
    if (permutation%all_permutations_done()) then
      exit
    endif
  enddo
end subroutine

! --------------------------------------------------
! Related algorithms (not currently in use).
! --------------------------------------------------
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
  
  indices = sort(input)
  call print_line(input(indices))
  
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
    call print_line(input(indices))
  enddo
end subroutine
end module
