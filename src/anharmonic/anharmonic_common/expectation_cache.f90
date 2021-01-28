! ======================================================================
! Caches the results of monomial integration,
!    to allow for faster repeat calculations.
! Uses C++ std::Vector-like storage of results.
! ======================================================================
module caesar_expectation_cache_module
  use caesar_common_module
  
  use caesar_sparse_monomial_module
  implicit none
  
  private
  
  public :: ExpectationCache
  
  type, extends(NoDefaultConstructor) :: ExpectationCache
    type(SparseMonomial), allocatable :: monomials(:)
    complex(dp),          allocatable :: expectations(:)
    
    integer, private :: length_
  contains
    procedure, public :: cached_location
    procedure, public :: cached_expectation
    procedure, public :: cache
  end type
  
  interface ExpectationCache
    module procedure new_ExpectationCache
  end interface
contains

! Constructor.
function new_ExpectationCache() result(this)
  implicit none
  
  type(ExpectationCache) :: this
  
  integer :: ialloc
  
  allocate( this%monomials(0),    &
          & this%expectations(0), &
          & stat=ialloc); call err(ialloc)
  this%length_ = 0
end function

! Get the location in where a expectation has been cached.
! Returns 0 if the expectation has not been cached.
function cached_location(this,monomial) result(output)
  implicit none
  
  class(ExpectationCache), intent(in) :: this
  type(SparseMonomial),    intent(in) :: monomial
  integer                             :: output
  
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
end function

! Get a cached expectation, given its location.
function cached_expectation(this,location) result(output)
  implicit none
  
  class(ExpectationCache), intent(in) :: this
  integer,                 intent(in) :: location
  complex(dp)                         :: output
  
  output = this%expectations(location)
end function

! Add an expectation to the cache.
subroutine cache(this,monomial,expectation)
  implicit none
  
  class(ExpectationCache), intent(inout) :: this
  type(SparseMonomial),    intent(in)    :: monomial
  complex(dp),             intent(in)    :: expectation
  
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
end subroutine

! Compare two monomials.
function compare_monomials(lhs,rhs) result(output)
  implicit none
  
  type(SparseMonomial), intent(in) :: lhs
  type(SparseMonomial), intent(in) :: rhs
  logical                          :: output
  
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
end function
end module
