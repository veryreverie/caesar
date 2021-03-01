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
    ! Constructor.
    module function new_ExpectationCache() result(this) 
      type(ExpectationCache) :: this
    end function
  end interface
  
  interface
    ! Get the location in where a expectation has been cached.
    ! Returns 0 if the expectation has not been cached.
    module function cached_location(this,monomial) result(output) 
      class(ExpectationCache), intent(in) :: this
      type(SparseMonomial),    intent(in) :: monomial
      integer                             :: output
    end function
  end interface
  
  interface
    ! Get a cached expectation, given its location.
    module function cached_expectation(this,location) result(output) 
      class(ExpectationCache), intent(in) :: this
      integer,                 intent(in) :: location
      complex(dp)                         :: output
    end function
  end interface
  
  interface
    ! Add an expectation to the cache.
    module subroutine cache(this,monomial,expectation) 
      class(ExpectationCache), intent(inout) :: this
      type(SparseMonomial),    intent(in)    :: monomial
      complex(dp),             intent(in)    :: expectation
    end subroutine
  end interface
  
  interface
    ! Compare two monomials.
    module function compare_monomials(lhs,rhs) result(output) 
      type(SparseMonomial), intent(in) :: lhs
      type(SparseMonomial), intent(in) :: rhs
      logical                          :: output
    end function
  end interface
end module
