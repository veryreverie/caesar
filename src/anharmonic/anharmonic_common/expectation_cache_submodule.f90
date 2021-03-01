submodule (caesar_expectation_cache_module) caesar_expectation_cache_submodule
  use caesar_anharmonic_common_module
contains

module procedure new_ExpectationCache
  integer :: ialloc
  
  allocate( this%monomials(0),    &
          & this%expectations(0), &
          & stat=ialloc); call err(ialloc)
  this%length_ = 0
end procedure

module procedure cached_location
  integer :: i
  
  output = 0
  
  if (.not. allocated(this%monomials)) then
    return
  endif
  
  do i=1,this%length_
    if (compare_monomials(this%monomials(i),monomial)) then
      output = i
      return
    endif
  enddo
end procedure

module procedure cached_expectation
  output = this%expectations(location)
end procedure

module procedure cache
  type(SparseMonomial), allocatable :: temp_monomials(:)
  complex(dp),          allocatable :: temp_expectations(:)
  
  integer :: i,ialloc
  
  if (.not. allocated(this%monomials)) then
    this%monomials = [monomial]
    this%expectations = [expectation]
    this%length_ = 1
  else
    do i=1,this%length_
      if (compare_monomials(this%monomials(i),monomial)) then
        call err()
      endif
    enddo
    
    this%length_ = this%length_ + 1
    if (this%length_<=size(this%monomials)) then
      this%monomials(this%length_) = monomial
      this%expectations(this%length_) = expectation
    else
      temp_monomials = this%monomials
      temp_expectations = this%expectations
      deallocate( this%monomials,    &
                & this%expectations, &
                & stat=ialloc); call err(ialloc)
      allocate( this%monomials(this%length_*2),    &
              & this%expectations(this%length_*2), &
              & stat=ialloc); call err(ialloc)
      this%monomials(:this%length_) = [temp_monomials, monomial]
      this%expectations(:this%length_) = [temp_expectations, expectation]
    endif
  endif
end procedure

module procedure compare_monomials
  if (size(lhs%modes)/=size(rhs%modes)) then
    output = .false.
  elseif (any(lhs%modes%id/=rhs%modes%id)) then
    output = .false.
  elseif (any(lhs%modes%power/=rhs%modes%power)) then
    output = .false.
  elseif (any(lhs%modes%paired_power/=rhs%modes%paired_power)) then
    output = .false.
  else
    output = .true.
  endif
end procedure
end submodule
