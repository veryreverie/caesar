!> A number of algorithms to do with logic and sorting.
module caesar_logic_module
  use caesar_foundations_module
  implicit none
  
  private
  
  public :: LogicalLambda
  public :: ComparisonLambda
  public :: OperationLambda
  public :: lazy_and
  public :: lazy_or
  public :: first
  public :: last
  public :: first_equivalent
  public :: last_equivalent
  public :: operate
  public :: map
  public :: count
  public :: filter
  public :: locate
  public :: is_sorted
  public :: sort
  public :: set
  public :: set_default
  public :: operator(.in.)
  
  interface
    !> A procedure argument which takes one argument and returns a logical.
    !> Intended for evaluation functions.
    function LogicalLambda(input) result(output)
      class(*), intent(in) :: input
      logical              :: output
    end function
    
    !> A procedure argument which takes two arguments and returns a logical.
    !> Intended for comparison functions.
    function ComparisonLambda(this,that) result(output)
      class(*), intent(in) :: this
      class(*), intent(in) :: that
      logical              :: output
    end function
    
    !> A procedure argument which takes one mutable argument.
    !> Intended for procedures which modify their argument.
    subroutine OperationLambda(this)
      class(*), intent(inout) :: this
    end subroutine
  end interface
  
  interface lazy_and
    !> Works as `.and.`, but the second argument may be absent if the first
    !>    argument is false.
    impure elemental module function lazy_and(this,that) result(output) 
      logical, intent(in)           :: this
      logical, intent(in), optional :: that
      logical                       :: output
    end function
  end interface
  
  interface lazy_or
    !> Works as `.or.`, but the second argument may be absent if the first
    !>    argument is true.
    impure elemental module function lazy_or(this,that) result(output) 
      logical, intent(in)           :: this
      logical, intent(in), optional :: that
      logical                       :: output
    end function
  end interface
  
  interface first
    !> Returns the index of the first `true` in `input`.
    recursive module function first_logicals(input,mask,default,sorted) &
       & result(output) 
      logical, intent(in)           :: input(:)
      !> If `mask` is present, only elements where `mask` is `true` are
      !>    considered.
      !> `mask` must have the same size as `input`.
      logical, intent(in), optional :: mask(:)
      !> `default` is returned if the array contains no `true` values.
      !> If `default` is not present, then an error is thrown in this case.
      integer, intent(in), optional :: default
      !> If `sorted`, then it is assumed that `input` is sorted (i.e. every
      !>    `true` appears after every `false`). The algorithm will then run
      !>    faster for long lists (`O(log(N))` rather `O(N)`).
      !> Defaults to `false`.
      !> WARNING: If `sorted` is true and the list is not sorted then the
      !>    function will fail silently.
      logical, intent(in), optional :: sorted
      integer                       :: output
    end function

    !> Returns the index `i` of the first element in `input` for which
    !>    `lambda(input(i))` returns `true`.
    module function first_LogicalLambda(input,lambda,mask,default,sorted) &
       & result(output) 
      class(*), intent(in)           :: input(:)
      procedure(LogicalLambda)       :: lambda
      !> If `mask` is present, only elements where `mask` is `true` are
      !>    considered.
      !> `mask` must have the same size as `input`.
      logical,  intent(in), optional :: mask(:)
      !> `default` is returned if no elements in `input` evaluate to `true`.
      !> If `default` is not present, then an error is thrown in this case.
      integer,  intent(in), optional :: default
      !> If `sorted`, then it is assumed that `input` is sorted (i.e. every
      !>    `true` appears after every `false`). The algorithm will then run
      !>    faster for long lists (`O(log(N))` rather `O(N)`).
      !> Defaults to `false`.
      !> WARNING: If `sorted` is true and the list is not sorted then the
      !>    function will fail silently.
      logical,  intent(in), optional :: sorted
      integer                        :: output
    end function
  end interface
  
  interface last
    !> Returns the index of the last `true` in `input`.
    recursive module function last_logicals(input,mask,default,sorted) &
       & result(output) 
      logical, intent(in)           :: input(:)
      !> If `mask` is present, only elements where `mask` is `true` are
      !>    considered.
      !> `mask` must have the same size as `input`.
      logical, intent(in), optional :: mask(:)
      !> `default` is returned if the array contains no `true` values.
      !> If `default` is not present, then an error is thrown in this case.
      integer, intent(in), optional :: default
      !> If `sorted`, then it is assumed that `input` is sorted (i.e. every
      !>    `true` appears before every `false`). The algorithm will then run
      !>    faster for long lists (`O(log(N))` rather `O(N)`).
      !> Defaults to `false`.
      !> WARNING: If `sorted` is true and the list is not sorted then the
      !>    function will fail silently.
      logical, intent(in), optional :: sorted
      integer                       :: output
    end function

    !> Returns the index of the last element in `input` for which
    !>    `lambda(input)` returns `true`.
    module function last_LogicalLambda(input,lambda,mask,default,sorted) &
       & result(output) 
      class(*), intent(in)           :: input(:)
      procedure(LogicalLambda)       :: lambda
      !> If `mask` is present, only elements where `mask` is `true` are
      !>    considered.
      !> `mask` must have the same size as `input`.
      logical,  intent(in), optional :: mask(:)
      !> `default` is returned if no elements in `input` evaluate to `true`.
      !> If `default` is not present, then an error is thrown in this case.
      integer,  intent(in), optional :: default
      !> If `sorted`, then it is assumed that `input` is sorted (i.e. every
      !>    `true` appears before every `false`). The algorithm will then run
      !>    faster for long lists (`O(log(N))` rather `O(N)`).
      !> Defaults to `false`.
      !> WARNING: If `sorted` is true and the list is not sorted then the
      !>    function will fail silently.
      logical,  intent(in), optional :: sorted
      integer                        :: output
    end function
  end interface
  
  interface first_equivalent
    !> Returns the index `i` of the first element in `input` for which
    !>    `input(i)==comparison`.
    module function first_equivalent_integers(input,comparison,mask,default, &
       & sorted) result(output) 
      integer, intent(in)           :: input(:)
      integer, intent(in)           :: comparison
      !> If `mask` is present, only elements where `mask` is `true` are
      !>    considered.
      !> `mask` must have the same size as `input`.
      logical, intent(in), optional :: mask(:)
      !> `default` is returned if no elements in `input` are
      !>    equal to `comparison`.
      !> If `default` is not present, then an error is thrown in this case.
      integer, intent(in), optional :: default
      !> If `sorted`, then it is assumed that `input` is sorted in ascending
      !>    order. The algorithm will then run faster for long lists
      !>    (`O(log(N))` rather `O(N)`).
      !> Defaults to `false`.
      !> WARNING: If `sorted` is true and the list is not sorted then the
      !>    function will fail silently.
      logical, intent(in), optional :: sorted
      integer                       :: output
    end function

    !> Returns the index `i` of the first element in `input` for which
    !>    `equality(input(i), comparison)` returns `true`.
    module function first_equivalent_ComparisonLambda(input,comparison, &
       & equality,greater_than_or_equal,mask,default,sorted) result(output) 
      class(*), intent(in)                  :: input(:)
      class(*), intent(in)                  :: comparison
      procedure(ComparisonLambda)           :: equality
      !> Only required if `sorted` is true.
      procedure(ComparisonLambda), optional :: greater_than_or_equal
      !> If `mask` is present, only elements where `mask` is `true` are
      !>    considered.
      !> `mask` must have the same size as `input`.
      logical,  intent(in),        optional :: mask(:)
      !> `default` is returned if no elements in `input` evaluate as
      !>    equal to `comparison`.
      !> If `default` is not present, then an error is thrown in this case.
      integer,  intent(in),        optional :: default
      !> If `sorted`, then it is assumed that `input` is sorted in ascending
      !>    order as defined by `greater_than_or_equal`.
      !> The algorithm will then run faster for long lists
      !>    (`O(log(N))` rather `O(N)`).
      !> Defaults to `false`.
      !> WARNING: If `sorted` is true and the list is not sorted then the
      !>    function will fail silently.
      logical,  intent(in),        optional :: sorted
      integer                               :: output
    end function
  end interface
  
  interface last_equivalent
    !> Returns the index `i` of the last element in `input` for which
    !>    `input(i)==comparison`.
    module function last_equivalent_integers(input,comparison,mask,default, &
       & sorted) result(output) 
      integer, intent(in)           :: input(:)
      integer, intent(in)           :: comparison
      !> If `mask` is present, only elements where `mask` is `true` are
      !>    considered.
      !> `mask` must have the same size as `input`.
      logical, intent(in), optional :: mask(:)
      !> `default` is returned if no elements in `input` are
      !>    equal to `comparison`.
      !> If `default` is not present, then an error is thrown in this case.
      integer, intent(in), optional :: default
      !> If `sorted`, then it is assumed that `input` is sorted in descending
      !>    order. The algorithm will then run faster for long lists
      !>    (`O(log(N))` rather `O(N)`).
      !> Defaults to `false`.
      !> WARNING: If `sorted` is true and the list is not sorted then the
      !>    function will fail silently.
      logical, intent(in), optional :: sorted
      integer                       :: output
    end function

    !> Returns the index `i` of the last element in `input` for which
    !>    `equality(input(i), comparison)` returns `true`.
    module function last_equivalent_ComparisonLambda(input,comparison, &
       & equality,greater_than_or_equal,mask,default,sorted) result(output) 
      class(*), intent(in)                  :: input(:)
      class(*), intent(in)                  :: comparison
      procedure(ComparisonLambda)           :: equality
      !> Only required if `sorted` is true.
      procedure(ComparisonLambda), optional :: greater_than_or_equal
      !> If `mask` is present, only elements where `mask` is `true` are
      !>    considered.
      !> `mask` must have the same size as `input`.
      logical,  intent(in),        optional :: mask(:)
      !> `default` is returned if no elements in `input` evaluate as
      !>    equal to `comparison`.
      !> If `default` is not present, then an error is thrown in this case.
      integer,  intent(in),        optional :: default
      !> If `sorted`, then it is assumed that `input` is sorted in ascending
      !>    order as defined by `greater_than_or_equal`.
      !> The algorithm will then run faster for long lists
      !>    (`O(log(N))` rather `O(N)`).
      !> Defaults to `false`.
      !> WARNING: If `sorted` is true and the list is not sorted then the
      !>    function will fail silently.
      logical,  intent(in),        optional :: sorted
      integer                               :: output
    end function
  end interface
  
  interface operate
    !> Operates with `lambda` on every element of `input`.
    !> e.g. if `input` is `[2,3,7]` then after `operate(list,add_one)`,
    !>    `input` will be `[3,4,8]`.
    module subroutine operate_OperationLambda(input,lambda) 
      class(*), intent(inout)    :: input(:)
      procedure(OperationLambda) :: lambda
    end subroutine
  end interface
  
  interface map
    !> Sets `output(i)` to `lambda(input(i))`.
    module function map_LogicalLambda(input,lambda) result(output) 
      class(*), intent(in)     :: input(:)
      procedure(LogicalLambda) :: lambda
      logical, allocatable     :: output(:)
    end function

    !> Sets `output(i)` to `lambda(input(i),comparison)`
    module function map_ComparisonLambda(input,lambda,comparison) &
       & result(output)
      class(*), intent(in)        :: input(:)
      procedure(ComparisonLambda) :: lambda
      class(*), intent(in)        :: comparison
      logical, allocatable        :: output(:)
    end function
  end interface
  
  interface count
    !> Returns the total number of elements of `input` for which
    !>    `lambda(input(i))` returns `true`.
    module function count_LogicalLambda(input,lambda) result(output) 
      class(*), intent(in)     :: input(:)
      procedure(LogicalLambda) :: lambda
      integer                  :: output
    end function

    !> Returns the total number of elements of `input` for which
    !>    `lambda(input(i),comparison)` returns `true`.
    module function count_ComparisonLambda(input,lambda,comparison) &
       & result(output) 
      class(*), intent(in)        :: input(:)
      procedure(ComparisonLambda) :: lambda
      class(*), intent(in)        :: comparison
      integer                     :: output
    end function
  end interface
  
  interface filter
    !> Returns an array containing the indices `{i}` for which `input(i)`
    !>    is `true`.
    module function filter_logicals(input,mask) result(output) 
      logical, intent(in)           :: input(:)
      !> If `mask` is present, each index `i` will only be included in `output`
      !>    if `mask(i)` is `true`.
      !> `mask` must have the same size as `input`.
      logical, intent(in), optional :: mask(:)
      integer, allocatable          :: output(:)
    end function

    !> Returns an array containing the indices `{i}` for which `mask(input(i))`
    !>    returns `true`.
    module function filter_LogicalLambda(input,lambda,mask) result(output) 
      class(*), intent(in)           :: input(:)
      procedure(LogicalLambda)       :: lambda
      !> If `mask` is present, each index `i` will only be included in `output`
      !>    if `mask(i)` is `true`.
      !> `mask` must have the same size as `input`.
      logical,  intent(in), optional :: mask(:)
      integer, allocatable           :: output(:)
    end function

    !> Returns an array containing the indices `{i}` for which
    !>    `mask(input(i),comparison)` returns `true`.
    module function filter_ComparisonLambda(input,lambda,comparison,mask) &
       & result(output) 
      class(*), intent(in)           :: input(:)
      procedure(ComparisonLambda)    :: lambda
      class(*), intent(in)           :: comparison
      !> If `mask` is present, each index `i` will only be included in `output`
      !>    if `mask(i)` is `true`.
      !> `mask` must have the same size as `input`.
      logical,  intent(in), optional :: mask(:)
      integer, allocatable           :: output(:)
    end function
  end interface
  
  interface locate
    !> Returns the index `i` of the first `input` element such that
    !>    `lambda(input(i),input(j))` does not return `false` for any indices
    !>    $j\neq i$.
    module function locate_ComparisonLambda(input,lambda,mask) result(output) 
      class(*), intent(in)          :: input(:)
      procedure(ComparisonLambda)   :: lambda
      !> If `mask` is present then only elements for which `mask` is `true` are
      !>    considered.
      !> `mask` must be the same size as `input`.
      logical, intent(in), optional :: mask(:)
      integer                       :: output
    end function
  end interface
  
  interface is_sorted
    !> Returns whether or not `input` is sorted in ascending order.
    module function is_sorted_integers(input) result(output) 
      integer, intent(in) :: input(:)
      logical             :: output
    end function

    !> Returns whether or not `input` is sorted in ascending order.
    module function is_sorted_reals(input) result(output) 
      real(dp), intent(in) :: input(:)
      logical              :: output
    end function

    !> Returns whether or not `input` is sorted according to `lambda`.
    !> `input(i)` and `input(i+1)` are in sorted order if
    !>     `lambda(input(i),input(i+1))` returns `true`, or if
    !>     `lambda(input(i+1),input(i))` returns `false`.
    module function is_sorted_ComparisonLambda(input,lambda) result(output) 
      class(*), intent(in)        :: input(:)
      procedure(ComparisonLambda) :: lambda
      logical                     :: output
    end function
  end interface
  
  interface sort
    !> Returns an array containing the indices `{i}` of `input` in order of
    !>    ascending `input(i)`.
    !> An array `input` can be sorted by calling `input(sort(input))`.
    module function sort_integers(input) result(output) 
      integer, intent(in)  :: input(:)
      integer, allocatable :: output(:)
    end function

    !> Returns an array containing the indices `{i}` of `input` in order of
    !>    ascending `input(i)`.
    !> An array `input` can be sorted by calling `input(sort(input))`.
    module function sort_reals(input) result(output) 
      real(dp), intent(in) :: input(:)
      integer, allocatable :: output(:)
    end function

    !> Returns an array containing the indices `{i}` of `input` in an order
    !>    defined by `lambda`.
    !> An array `input` can be sorted by calling `input(sort(input))`.
    !> The sorted list will satisfy `lambda(input(i),input(i+1))=true`.
    module function sort_ComparisonLambda(input,lambda) result(output) 
      class(*), intent(in)        :: input(:)
      procedure(ComparisonLambda) :: lambda
      integer, allocatable        :: output(:)
    end function
  end interface
  
  interface set
    !> Finds a set in the mathematical sense, i.e. an array which contains
    !>    every entry in `input`, but only one copy of each duplicate entry.
    !> Returns the indices corresponding to the first copy of each entry.
    !> The set of `input` can be found by calling `input(set(input))`.
    !> e.g. if `input` is `[8,8,7,8,2,7]` then `set(input)` returns `[1,3,5]`,
    !>    and `input(set(input))` returns `[8,7,2]`.
    module function set_integers(input,mask) result(output) 
      integer, intent(in)           :: input(:)
      !> If `mask` is present then only elements for which `mask` is `true` are
      !>    considered.
      !> `mask` must be the same size as `input`.
      logical, intent(in), optional :: mask(:)
      integer, allocatable          :: output(:)
    end function

    !> Finds a set in the mathematical sense, i.e. an array which contains
    !>    every entry in `input`, but only one copy of each entry.
    !> Returns the indices corresponding to the first copy of each entry.
    !> The set of `input` can be found by calling `input(set(input,lambda))`.
    !> Entries `i` and `j` are considered duplicates if
    !>    `lambda(input(i),input(j))` returns `true`.
    !> e.g. if `list` is `['A','A','D','A','F','D']` and `lambda` is string
    !>    comparison, then `set(input,lambda)` returns `[1,3,5]`,
    !>    and `input(set(input,lambda))` returns `['A','D','F']`.
    module function set_ComparisonLambda(input,lambda,mask) result(output) 
      class(*), intent(in)          :: input(:)
      procedure(ComparisonLambda)   :: lambda
      !> If `mask` is present then only elements for which `mask` is `true` are
      !>    considered.
      !> `mask` must be the same size as `input`.
      logical, intent(in), optional :: mask(:)
      integer, allocatable          :: output(:)
    end function
  end interface
  
  interface set_default
    !> Takes an `optional_argument`, and its `default_value`.
    !> Returns the `optional_argument` if present, or the `default_value`
    !>    if not.
    impure elemental module function set_default_logical(optional_argument, &
       & default_value) result(output) 
      logical, intent(in), optional :: optional_argument
      logical, intent(in)           :: default_value
      logical                       :: output
    end function
    
    !> Takes an `optional_argument`, and its `default_value`.
    !> Returns the `optional_argument` if present, or the `default_value`
    !>    if not.
    impure elemental module function set_default_integer(optional_argument, &
       & default_value) result(output) 
      integer, intent(in), optional :: optional_argument
      integer, intent(in)           :: default_value
      integer                       :: output
    end function
  end interface
  
  interface operator(.in.)
    !> Given an integer `lhs` and an array of integers `rhs`,
    !>    `lhs.in.rhs` returns `true` if `lhs` appears in `rhs`,
    !>    and `false` otherwise.
    module function element_in_list(lhs,rhs) result(output) 
      integer, intent(in) :: lhs
      integer, intent(in) :: rhs(:)
      logical             :: output
    end function

    !> Given two arrays of integers, `lhs` and `rhs`,
    !>    `lhs.in.rhs` returns an array `output` with the same size as `lhs`,
    !>    where `output(i)` is `true` if `lhs(i)` appears in `rhs`,
    !>    and `false` otherwise.
    module function elements_in_list(lhs,rhs) result(output) 
      integer, intent(in)  :: lhs(:)
      integer, intent(in)  :: rhs(:)
      logical, allocatable :: output(:)
    end function
  end interface
end module
