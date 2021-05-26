submodule (caesar_expectation_cache_module) caesar_expectation_cache_submodule
  use caesar_anharmonic_common_module
contains

module procedure new_ExpectationCache
  integer :: ialloc
  
  allocate( this%monomials_(0),    &
          & this%expectations_(0), &
          & stat=ialloc); call err(ialloc)
  this%length_ = 0
  this%last_monomial_index_ = 0
end procedure

module procedure is_in_cache
  output = cache_index(this, monomial)>0
end procedure

module procedure expectation
  integer :: index
  
  index = cache_index(this, monomial)
  
  if (index>0) then
    output = this%expectations_(index)
  else
    call print_line(ERROR//': Monomial not in cache.')
    call err()
  endif
end procedure

module procedure cache
  integer :: index
  
  type(SparseMonomial), allocatable :: temp_monomials(:)
  complex(dp),          allocatable :: temp_expectations(:)
  
  integer :: ialloc
  
  index = -cache_index(this, monomial)
  
  if (index<0) then
    call print_line(ERROR//': Monomial already in cache.')
    call err()
  endif
  
  this%length_ = this%length_+1
  if (this%length_>size(this%monomials_)) then
    temp_monomials = this%monomials_
    temp_expectations = this%expectations_
    deallocate( this%monomials_,    &
              & this%expectations_, &
              & stat=ialloc); call err(ialloc)
    allocate( this%monomials_(2*this%length_),    &
            & this%expectations_(2*this%length_), &
            & stat=ialloc); call err(ialloc)
    this%monomials_(:index-1) = temp_monomials(:index-1)
    this%monomials_(index) = monomial
    this%monomials_(index+1:this%length_) = &
       & temp_monomials(index:this%length_-1)
    this%expectations_(:index-1) = temp_expectations(:index-1)
    this%expectations_(index) = expectation
    this%expectations_(index+1:this%length_) = &
       & temp_expectations(index:this%length_-1)
  else
    this%monomials_(index+1:this%length_) = &
       & this%monomials_(index:this%length_-1)
    this%monomials_(index) = monomial
    this%expectations_(index+1:this%length_) = &
       & this%expectations_(index:this%length_-1)
    this%expectations_(index) = expectation
  endif
end procedure

! If the `monomial` is in the cache, returns the index `i` such that
!    cache%monomials_(i)==monomial.
! Otherwise, returns the negative of the index `i` such that
!    cache%monomials_(i-1)<monomial<cache%monomials_(i).
! Also caches the input monomial and output index for use in future calls.
function cache_index(cache,monomial) result(output)
  type(ExpectationCache), intent(inout) :: cache
  type(SparseMonomial),   intent(in)    :: monomial
  integer                               :: output
  
  ! Lower and upper bounds for the bisection search.
  integer :: lower
  integer :: upper
  integer :: centre
  
  ! The result of calling compare_monomials.
  integer :: comparison
  
  if (.not. allocated(cache%monomials_)) then
    call print_line(ERROR//': Cache not initialised.')
    call err()
  endif
  
  ! Handle the case of an empty cache.
  ! Also ensures that last_monomial_ is not checked before it is set.
  if (size(cache%monomials_)==0) then
    output = -1
    cache%last_monomial_ = monomial
    cache%last_monomial_index_ = output
    return
  endif
  
  ! Check the previously cached results.
  if (compare_monomials(monomial, cache%last_monomial_)==0) then
    output = cache%last_monomial_index_
    return
  endif
  
  ! Run a bisection search to find the monomial.
  upper = cache%length_+1
  lower = 0
  do
    centre = lower+(upper-lower)/2
    comparison = compare_monomials(monomial, cache%monomials_(centre))
    
    if (comparison==1) then
      lower = centre
    elseif (comparison==-1) then
      upper = centre
    else
      output = centre
      exit
    endif
    
    if (upper-lower==1) then
      output = -upper
      exit
    endif
  enddo
  
  ! Cache the results for future calls.
  cache%last_monomial_ = monomial
  cache%last_monomial_index_ = output
end function

! Compare two monomials. Returns:
!   -1 if lhs<rhs
!    0 if lhs=rhs
!    1 if lhs>rhs
! The ordering is arbitrary.
function compare_monomials(lhs,rhs) result(output) 
  type(SparseMonomial), intent(in) :: lhs
  type(SparseMonomial), intent(in) :: rhs
  integer                          :: output
  
  integer :: i
  
  if (size(lhs%modes)<size(rhs%modes)) then
    output = -1
  elseif (size(lhs%modes)>size(rhs%modes)) then
    output = 1
  else
    do i=1,size(lhs%modes)
      if (lhs%modes(i)%id<rhs%modes(i)%id) then
        output = -1
        return
      elseif (lhs%modes(i)%id>rhs%modes(i)%id) then
        output = 1
        return
      elseif (lhs%modes(i)%power<rhs%modes(i)%power) then
        output = -1
        return
      elseif (lhs%modes(i)%power>rhs%modes(i)%power) then
        output = 1
        return
      elseif (lhs%modes(i)%paired_power<rhs%modes(i)%paired_power) then
        output = -1
        return
      elseif (lhs%modes(i)%paired_power>rhs%modes(i)%paired_power) then
        output = 1
        return
      endif
    enddo
    
    output = 0
  endif
end function
end submodule
