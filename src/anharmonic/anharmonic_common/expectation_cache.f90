!> Provides the [[ExpectationCache(type)]] class, and related methods.
module caesar_expectation_cache_module
  use caesar_common_module
  
  use caesar_sparse_monomial_module
  implicit none
  
  private
  
  public :: ExpectationCache
  
  !> Caches ("memoises") the results of monomial integration,
  !>    to allow for faster repeated calculations.
  !> Also caches the location in the cache of the last procedure call,
  !>    so multiple procedure calls with the same monomial argument will run
  !>    faster.
  type, extends(NoDefaultConstructor) :: ExpectationCache
    !> The list of cached monomials.
    !> Stored in ascending order (by an internal sort function).
    !> Uses the C++ std::Vector memory model (doubling in size on overflow).
    type(SparseMonomial), allocatable, private :: monomials_(:)
    !> `expectations_(i)` is the expectation of `monomials_(i)`.
    complex(dp),          allocatable, private :: expectations_(:)
    !> The number of entries of monomials_ and expectations_ which are
    !>    actually used.
    integer, private :: length_
    
    !> The monomial input to the last call.
    type(SparseMonomial), private :: last_monomial_
    !> The index in `monomials_` of `last_monomial_`.
    !> (or where this monomial should go when cached.)
    integer,              private :: last_monomial_index_
  contains
    procedure, public :: is_in_cache
    procedure, public :: expectation
    procedure, public :: cache
  end type
  
  interface ExpectationCache
    ! Constructor for objects of type [[ExpectationCache(type)]].
    module function new_ExpectationCache() result(this) 
      type(ExpectationCache) :: this
    end function
  end interface
  
  interface
    !> Returns whether or not `monomial` is in the cache.
    !> Also caches the location for future use.
    module function is_in_cache(this,monomial) result(output) 
      class(ExpectationCache), intent(inout) :: this
      type(SparseMonomial),    intent(in)    :: monomial
      logical                                :: output
    end function
  end interface
  
  interface
    !> Returns the previously-cached expectation of a monomial.
    !> Throws an error if the monomial is not in the cache.
    module function expectation(this,monomial) result(output) 
      class(ExpectationCache), intent(inout) :: this
      type(SparseMonomial),    intent(in)    :: monomial
      complex(dp)                            :: output
    end function
  end interface
  
  interface
    !> Cache the expectation of a monomial.
    !> Throws an error if the monomial is in the cache.
    module subroutine cache(this,monomial,expectation) 
      class(ExpectationCache), intent(inout) :: this
      type(SparseMonomial),    intent(in)    :: monomial
      complex(dp),             intent(in)    :: expectation
    end subroutine
  end interface
end module
