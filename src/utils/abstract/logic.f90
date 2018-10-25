! ======================================================================
! A number of algorithms to do with logic and sorting.
! ======================================================================
! See the example module for some use cases.
module logic_module
  use precision_module
  use io_basic_module
  implicit none
  
  private
  
  ! Functions and subroutines.
  public :: get
  public :: lazy_and
  public :: lazy_or
  public :: first
  !public :: last
  public :: first_equivalent
  !public :: last_equivalent
  public :: operate
  public :: map
  public :: count
  public :: filter
  public :: locate
  public :: is_sorted
  public :: sort
  public :: set
  
  ! Lambdas.
  public :: LogicalLambda
  public :: ComparisonLambda
  public :: OperationLambda
  
  interface first
    module procedure first_logicals
    module procedure first_LogicalLambda
  end interface
  
  !interface last
  !  module procedure last_logicals
  !  module procedure last_LogicalLambda
  !end interface
  
  interface first_equivalent
    module procedure first_equivalent_integers
    module procedure first_equivalent_ComparisonLambda
  end interface
  
  !interface last_equivalent
  !  module procedure last_equivalent_integers
  !  module procedure last_equivalent_ComparisonLambda
  !end interface
  
  interface operate
    module procedure operate_OperationLambda
  end interface
  
  interface map
    module procedure map_LogicalLambda
    module procedure map_ComparisonLambda
  end interface
  
  interface count
    module procedure count_LogicalLambda
    module procedure count_ComparisonLambda
  end interface
  
  interface filter
    module procedure filter_logicals
    module procedure filter_LogicalLambda
    module procedure filter_ComparisonLambda
  end interface
  
  interface locate
    module procedure locate_ComparisonLambda
  end interface
  
  interface is_sorted
    module procedure is_sorted_integers
    module procedure is_sorted_ComparisonLambda
  end interface
  
  interface sort
    module procedure sort_integers
    module procedure sort_ComparisonLambda
  end interface
  
  interface set
    module procedure set_integers
    module procedure set_ComparisonLambda
  end interface
  
  interface
    function LogicalLambda(input) result(output)
      implicit none
      
      class(*), intent(in) :: input
      logical              :: output
    end function
    
    function ComparisonLambda(this,that) result(output)
      implicit none
      
      class(*), intent(in) :: this
      class(*), intent(in) :: that
      logical              :: output
    end function
    
    subroutine OperationLambda(this)
      implicit none
      
      class(*), intent(inout) :: this
    end subroutine
  end interface
contains

! ----------------------------------------------------------------------
! Takes an optional logical input, and a default value.
! Returns the input if present, or the default otherwise.
! ----------------------------------------------------------------------
impure elemental function get(input,default) result(output)
  implicit none
  
  logical, intent(in), optional :: input
  logical, intent(in)           :: default
  logical                       :: output
  
  if (present(input)) then
    output = input
  else
    output = default
  endif
end function

! ----------------------------------------------------------------------
! Lazy logical operators. .lazyand. and .lazyor. work as .and. and .or.
!    respectively, except that they are guaranteed not to look at the second
!    argument if not necessary.
! i.e. (.true. .lazyor. x) will return .true. and (.false. .lazyand. x) will
!    return .false. even if present(x)==.false..
! ----------------------------------------------------------------------
impure elemental function lazy_and(this,that) result(output)
  implicit none
  
  logical, intent(in)           :: this
  logical, intent(in), optional :: that
  logical                       :: output
  
  if (.not. this) then
    output = .false.
  elseif (present(that)) then
    output = that
  else
    call print_line(ERROR//': First argument to lazy_and is .true. and &
       &second argument is not present.')
    call err()
  endif
end function

impure elemental function lazy_or(this,that) result(output)
  implicit none
  
  logical, intent(in)           :: this
  logical, intent(in), optional :: that
  logical                       :: output
  
  if (this) then
    output = .true.
  elseif (present(that)) then
    output = that
  else
    call print_line(ERROR//': First argument to lazy_or is .false. and &
       &second argument is not present.')
    call err()
  endif
end function

! ----------------------------------------------------------------------
! The _logicals variants return the id of the first/last .true. in a 
!    list of logicals. e.g. first([f,f,t,f]) returns 3.
! The _LogicalLambda variants do the same, but for the first/last element for
!    which the lambda retuns .true..
!    e.g. first([1,1,2,1],equals_2) returns 3.
! The _ComparisonLambda variants do the same but for the first/last element
!    which compares to .true. with the given comparison value.
! If default is given, then this value is returned if all(input)==false.
! If default is not given, then an error is thrown if all(input)==false.
! If mask is present, only returns values where mask is .true..
! If sorted is .true., then the input to 'first' is assumed to be an array of
!    .false. followed by an array of .true.. For 'last', the input is assumed
!    to be in the reverse order. This allows for a bisection algorithm.
! If sorted is .true. and mask is present, then only the elements of input
!    where mask is .true. need be sorted.
! WARNING: since 'sorted' is an optimisation, if sorted is set to true when
!    the list is not sorted, it may fail silently. There is no way to test
!    for this without incurring an O(n) cost.
! ----------------------------------------------------------------------

! Details of the bisection method:
!
! Denoting known values of false and true as        (1)      f         t
!    F and T respectively, unknown unmasked              [...F??XXX?XX?T...]
!    values as ?, and masked values as X,
! The input to each step is  shown in (1),          (2)      f j  i    t
!    with the last known false and first known           [...F??XXX?XX?T...]
!    true marked f and t respectively.                       |<-->|<-->|
! The unknown range is bisected, to give the
!    centre-point (marked i in (2)), then the       (3a)          f    t
!    last un-masked value before this                    [...F?FXXX?XX?T...]
!    (marked j in (2)) is evaluated.
! If j is false, all values <= i are either         (3b)     f t
!    false or masked, so f is moved to i (3a).           [...F?TXXX?XX?T...]
! If j is true, all values >= j are either
!    true or masked, so t is moved to j (3b).       (3c)          f    t
! If j = f, then f is moved to i (3c).                   [...FXXXXX?XX?T...]

! N.B. first and last are called a lot, so they are somewhat optimised.
!    This is why they are not as compact as they could otherwise be.

function first_logicals(input,mask,default,sorted) result(output)
  implicit none
  
  logical, intent(in)           :: input(:)
  logical, intent(in), optional :: mask(:)
  integer, intent(in), optional :: default
  logical, intent(in), optional :: sorted
  integer                       :: output
  
  logical :: input_is_sorted
  
  integer :: first_true
  integer :: last_false
  
  integer :: i,j
  
  ! Check that mask is consistent with input.
  if (present(mask)) then
    if (size(input)/=size(mask)) then
      call print_line(CODE_ERROR//': Input and Mask of different sizes.')
      call err()
    endif
  endif
  
  ! Default sorted to false.
  if (present(sorted)) then
    input_is_sorted = sorted
  else
    input_is_sorted = .false.
  endif
  
  ! Find the first true.
  last_false = 0
  first_true = size(input)+1
  if (input_is_sorted) then
    ! If the input is sorted, use a bisection method to find the first true.
    if (present(mask)) then
      do while(first_true-last_false>1)
        ! Bisect the input between the last known false and first known true.
        i = (last_false+first_true)/2
        ! Find the last unmasked value before i.
        do j=i,last_false,-1
          if (mask(j)) then
            exit
          endif
        enddo
        ! Evaluate input(j), and update the bounds accordingly.
        if (j==last_false) then
          last_false = i
        else
          if (input(j)) then
            first_true = j
          else
            last_false = i
          endif
        endif
      enddo
    else
      do while(first_true-last_false>1)
        ! Bisect the input between the last known false and first known true.
        i = (last_false+first_true)/2
        ! Evaluate input(i), and update the bounds accordingly.
        if (input(i)) then
          first_true = i
        else
          last_false = i
        endif
      enddo
    endif
  else
    ! If the input is not sorted, simply iterate through the input.
    if (present(mask)) then
      do i=1,size(input)
        if (mask(i)) then
          if (input(i)) then
            first_true = i
            exit
          endif
        endif
      enddo
    else
      do i=1,size(input)
        if (input(i)) then
          first_true = i
          exit
        endif
      enddo
    endif
  endif
  
  ! Set the output to the first true, or the default value if present.
  ! Throw an error if no true is found and no default is set.
  if (first_true<=size(input)) then
    output = first_true
  elseif (present(default)) then
    output = default
  else
    call print_line(ERROR//': value not found.')
    call err()
  endif
end function

function first_LogicalLambda(input,lambda,mask,default,sorted) result(output)
  implicit none
  
  class(*), intent(in)           :: input(:)
  procedure(LogicalLambda)       :: lambda
  logical,  intent(in), optional :: mask(:)
  integer,  intent(in), optional :: default
  logical,  intent(in), optional :: sorted
  integer                        :: output
  
  logical :: input_is_sorted
  
  integer :: first_true
  integer :: last_false
  
  integer :: i,j
  
  ! Check that mask is consistent with input.
  if (present(mask)) then
    if (size(mask)/=size(input)) then
      call print_line(CODE_ERROR//': Input and Mask of different sizes.')
      call err()
    endif
  endif
  
  ! Default sorted to false.
  if (present(sorted)) then
    input_is_sorted = sorted
  else
    input_is_sorted = .false.
  endif
  
  ! Find the first true.
  last_false = 0
  first_true = size(input)+1
  if (input_is_sorted) then
    ! If the input is sorted, use a bisection method to find the first true.
    if (present(mask)) then
      do while(first_true-last_false>1)
        ! Bisect the input between the last known false and first known true.
        i = (last_false+first_true)/2
        ! Find the last unmasked value before i.
        do j=i,last_false,-1
          if (mask(j)) then
            exit
          endif
        enddo
        ! Evaluate input(j), and update the bounds accordingly.
        if (j==last_false) then
          last_false = i
        else
          if (lambda(input(j))) then
            first_true = j
          else
            last_false = i
          endif
        endif
      enddo
    else
      do while(first_true-last_false>1)
        ! Bisect the input between the last known false and first known true.
        i = (last_false+first_true)/2
        ! Evaluate input(i), and update the bounds accordingly.
        if (lambda(input(i))) then
          first_true = i
        else
          last_false = i
        endif
      enddo
    endif
  else
    ! If the input is not sorted, simply iterate through the input.
    if (present(mask)) then
      do i=1,size(input)
        if (mask(i)) then
          if (lambda(input(i))) then
            first_true = i
            exit
          endif
        endif
      enddo
    else
      do i=1,size(input)
        if (lambda(input(i))) then
          first_true = i
          exit
        endif
      enddo
    endif
  endif
  
  ! Set the output to the first true, or the default value if present.
  ! Throw an error if no true is found and no default is set.
  if (first_true<=size(input)) then
    output = first_true
  elseif (present(default)) then
    output = default
  else
    call print_line(ERROR//': value not found.')
    call err()
  endif
end function

!function last_logicals(input,mask,default,sorted) result(output)
!  implicit none
!  
!  logical, intent(in)           :: input(:)
!  logical, intent(in), optional :: mask(:)
!  integer, intent(in), optional :: default
!  logical, intent(in), optional :: sorted
!  integer                       :: output
!  
!  logical :: input_is_sorted
!  
!  integer :: last_true
!  integer :: first_false
!  
!  integer :: i,j
!  
!  ! Check that mask is consistent with input.
!  if (present(mask)) then
!    if (size(input)/=size(mask)) then
!      call print_line(CODE_ERROR//': Input and Mask of different sizes.')
!      call err()
!    endif
!  endif
!  
!  ! Default sorted to false.
!  if (present(sorted)) then
!    input_is_sorted = sorted
!  else
!    input_is_sorted = .false.
!  endif
!  
!  ! Find the last true.
!  last_true = 0
!  first_false = size(input)+1
!  if (input_is_sorted) then
!    ! If the input is sorted, use a bisection method to find the first true.
!    if (present(mask)) then
!      do while(first_false-last_true>1)
!        ! Bisect the input between the last known false and first known true.
!        i = (last_true+first_false+1)/2
!        ! Find the first unmasked value after i.
!        j = first(mask(i:first_false-1), default=0)
!        ! Evaluate input(j), and update the bounds accordingly.
!        if (j==0) then
!          first_false = i
!        else
!          j = i-1+j
!          if (input(j)) then
!            last_true = j
!          else
!            first_false = i
!          endif
!        endif
!      enddo
!    else
!      do while(first_false-last_true>1)
!        ! Bisect the input between the last known false and first known true.
!        i = (last_true+first_false)/2
!        ! Evaluate input(i), and update the bounds accordingly.
!        if (input(i)) then
!          last_true = i
!        else
!          first_false = i
!        endif
!      enddo
!    endif
!  else
!    ! If the input is not sorted, simply iterate through the input.
!    if (present(mask)) then
!      do i=size(input),1,-1
!        if (mask(i)) then
!          if (input(i)) then
!            last_true = i
!            exit
!          endif
!        endif
!      enddo
!    else
!      do i=size(input),1,-1
!        if (input(i)) then
!          last_true = i
!          exit
!        endif
!      enddo
!    endif
!  endif
!  
!  ! Set the output to the last true, or the default value if present.
!  ! Throw an error if no true is found and no default is set.
!  if (last_true>0) then
!    output = last_true
!  elseif (present(default)) then
!    output = default
!  else
!    call print_line(ERROR//': value not found.')
!    call err()
!  endif
!end function
!
!function last_LogicalLambda(input,lambda,mask,default,sorted) result(output)
!  implicit none
!  
!  class(*), intent(in)           :: input(:)
!  procedure(LogicalLambda)       :: lambda
!  logical,  intent(in), optional :: mask(:)
!  integer,  intent(in), optional :: default
!  logical,  intent(in), optional :: sorted
!  integer                        :: output
!  
!  logical :: input_is_sorted
!  
!  integer :: first_false
!  integer :: last_true
!  
!  integer :: i,j
!  
!  ! Check that mask is consistent with input.
!  if (present(mask)) then
!    if (size(mask)/=size(input)) then
!      call print_line(CODE_ERROR//': Input and Mask of different sizes.')
!      call err()
!    endif
!  endif
!  
!  ! Default sorted to false.
!  if (present(sorted)) then
!    input_is_sorted = sorted
!  else
!    input_is_sorted = .false.
!  endif
!  
!  last_true = 0
!  first_false = size(input)+1
!  if (input_is_sorted) then
!    ! If the input is sorted, use a bisection method to find the first true.
!    if (present(mask)) then
!      do while(first_false-last_true>1)
!        ! Bisect the input between the last known false and first known true.
!        i = (last_true+first_false+1)/2
!        ! Find the first unmasked value after i.
!        j = first(mask(i:first_false-1), default=0)
!        ! Evaluate input(j), and update the bounds accordingly.
!        if (j==0) then
!          first_false = i
!        else
!          j = i-1+j
!          if (lambda(input(j))) then
!            last_true = j
!          else
!            first_false = i
!          endif
!        endif
!      enddo
!    else
!      do while(first_false-last_true>1)
!        ! Bisect the input between the last known false and first known true.
!        i = (last_true+first_false)/2
!        ! Evaluate input(i), and update the bounds accordingly.
!        if (lambda(input(i))) then
!          last_true = i
!        else
!          first_false = i
!        endif
!      enddo
!    endif
!  else
!    ! If the input is not sorted, simply iterate through the input.
!    if (present(mask)) then
!      do i=size(input),1,-1
!        if (mask(i)) then
!          if (lambda(input(i))) then
!            last_true = i
!            exit
!          endif
!        endif
!      enddo
!    else
!      do i=size(input),1,-1
!        if (lambda(input(i))) then
!          last_true = i
!          exit
!        endif
!      enddo
!    endif
!  endif
!  
!  ! Set the output to the last true, or the default value if present.
!  ! Throw an error if no true is found and no default is set.
!  if (last_true>0) then
!    output = last_true
!  elseif (present(default)) then
!    output = default
!  else
!    call print_line(ERROR//': value not found.')
!    call err()
!  endif
!end function

! ----------------------------------------------------------------------
! first_equivalent(input,i) is equivalent to first(input==i),
!    but is faster than first() if sorted=true.
!
! If sorted=true then:
!    - For the _integers variant, the list must be in ascending order.
!
!    - For the _ComparisonLambda variant, inequality must be passed,
!      and the inequality must evaluate to .false. for the first part of the
!      list, and .true. for the second (or vice versa for last_equivalent).
! ----------------------------------------------------------------------
function first_equivalent_integers(input,comparison,mask,default,sorted) &
   & result(output)
  implicit none
  
  integer, intent(in)           :: input(:)
  integer, intent(in)           :: comparison
  logical, intent(in), optional :: mask(:)
  integer, intent(in), optional :: default
  logical, intent(in), optional :: sorted
  integer                       :: output
  
  logical :: input_is_sorted
  
  integer :: first_as_large
  integer :: last_smaller
  
  integer :: i,j
  
  ! Check that mask is consistent with input.
  if (present(mask)) then
    if (size(input)/=size(mask)) then
      call print_line(CODE_ERROR//': Input and Mask of different sizes.')
      call err()
    endif
  endif
  
  ! Default sorted to false.
  if (present(sorted)) then
    input_is_sorted = sorted
  else
    input_is_sorted = .false.
  endif
  
  ! Find the first integer which is at least as large as the comparison.
  last_smaller = 0
  first_as_large = size(input)+1
  if (input_is_sorted) then
    ! If the input is sorted, use a bisection method.
    if (present(mask)) then
      do while(first_as_large-last_smaller>1)
        ! Bisect the input between the last known value smaller than the
        !    comparison and first known value as large as the comparison.
        i = (last_smaller+first_as_large)/2
        ! Find the last unmasked value before i.
        do j=i,last_smaller,-1
          if (mask(j)) then
            exit
          endif
        enddo
        ! Evaluate input(j)>=comparison, and update the bounds accordingly.
        if (j==last_smaller) then
          last_smaller = i
        else
          if (input(j)>=comparison) then
            first_as_large = j
          else
            last_smaller = i
          endif
        endif
      enddo
    else
      do while(first_as_large-last_smaller>1)
        ! Bisect the input between the last known value smaller than the
        !    comparison and first known value as large as the comparison.
        i = (last_smaller+first_as_large)/2
        ! Evaluate input(i)>=comparison, and update the bounds accordingly.
        if (input(i)>=comparison) then
          first_as_large = i
        else
          last_smaller = i
        endif
      enddo
    endif
  else
    ! If the input is not sorted, simply iterate through the input.
    if (present(mask)) then
      do i=1,size(input)
        if (mask(i)) then
          if (input(i)>=comparison) then
            first_as_large = i
            exit
          endif
        endif
      enddo
    else
      do i=1,size(input)
        if (input(i)>=comparison) then
          first_as_large = i
          exit
        endif
      enddo
    endif
  endif
  
  ! Check whether the first value as large as the comparison is equal to the
  !    comparison. If it is, return its position.
  if (first_as_large<=size(input)) then
    if (input(first_as_large)==comparison) then
      output = first_as_large
      return
    endif
  endif
  
  ! If the comparison is not present, return the default if one is given,
  !    or throw an error if not.
  if (present(default)) then
    output = default
  else
    call print_line(ERROR//': value not found.')
    call err()
  endif
end function

function first_equivalent_ComparisonLambda(input,comparison,equality, &
   & greater_than_or_equal,mask,default,sorted) result(output)
  implicit none
  
  class(*), intent(in)                  :: input(:)
  class(*), intent(in)                  :: comparison
  procedure(ComparisonLambda)           :: equality
  procedure(ComparisonLambda), optional :: greater_than_or_equal
  logical,  intent(in),        optional :: mask(:)
  integer,  intent(in),        optional :: default
  logical,  intent(in),        optional :: sorted
  integer                               :: output
  
  logical :: input_is_sorted
  
  integer :: first_as_large
  integer :: last_smaller
  
  integer :: i,j
  
  ! Check that mask is consistent with input.
  if (present(mask)) then
    if (size(input)/=size(mask)) then
      call print_line(CODE_ERROR//': Input and Mask of different sizes.')
      call err()
    endif
  endif
  
  ! Default sorted to false.
  if (present(sorted)) then
    input_is_sorted = sorted
  else
    input_is_sorted = .false.
  endif
  
  ! Check that if sorted is true that inequality has been given.
  if (input_is_sorted .and. .not.present(greater_than_or_equal)) then
    call print_line(CODE_ERROR//': sorted is true, but no >= function has &
       &been given')
    call err()
  endif
  
  ! Find the first integer which is at least as large as the comparison.
  last_smaller = 0
  first_as_large = size(input)+1
  if (input_is_sorted) then
    ! If the input is sorted, use a bisection method.
    if (present(mask)) then
      do while(first_as_large-last_smaller>1)
        ! Bisect the input between the last known value smaller than the
        !    comparison and first known value as large as the comparison.
        i = (last_smaller+first_as_large)/2
        ! Find the last unmasked value before i.
        do j=i,last_smaller,-1
          if (mask(j)) then
            exit
          endif
        enddo
        ! Evaluate input(j)>=comparison, and update the bounds accordingly.
        if (j==last_smaller) then
          last_smaller = i
        else
          if (greater_than_or_equal(input(j),comparison)) then
            first_as_large = j
          else
            last_smaller = i
          endif
        endif
      enddo
    else
      do while(first_as_large-last_smaller>1)
        ! Bisect the input between the last known value smaller than the
        !    comparison and first known value as large as the comparison.
        i = (last_smaller+first_as_large)/2
        ! Evaluate input(i)>=comparison, and update the bounds accordingly.
        if (greater_than_or_equal(input(i),comparison)) then
          first_as_large = i
        else
          last_smaller = i
        endif
      enddo
    endif
  else
    ! If the input is not sorted, simply iterate through the input.
    if (present(greater_than_or_equal) .and. present(mask)) then
      do i=1,size(input)
        if (mask(i)) then
          if (greater_than_or_equal(input(i),comparison)) then
            first_as_large = i
            exit
          endif
        endif
      enddo
    elseif (present(greater_than_or_equal)) then
      do i=1,size(input)
        if (greater_than_or_equal(input(i),comparison)) then
          first_as_large = i
          exit
        endif
      enddo
    elseif (present(mask)) then
      do i=1,size(input)
        if (mask(i)) then
          if (equality(input(i),comparison)) then
            first_as_large = i
            exit
          endif
        endif
      enddo
    else
      do i=1,size(input)
        if (equality(input(i),comparison)) then
          first_as_large = i
          exit
        endif
      enddo
    endif
  endif
  
  ! Check whether the first value as large as the comparison is equal to the
  !    comparison. If it is, return its position.
  if (first_as_large<=size(input)) then
    if (equality(input(first_as_large),comparison)) then
      output = first_as_large
      return
    endif
  endif
  
  ! If the comparison is not present, return the default if one is given,
  !    or throw an error if not.
  if (present(default)) then
    output = default
  else
    call print_line(ERROR//': value not found.')
    call err()
  endif
end function

!function last_equivalent_integers(input,comparison,mask,default,sorted) &
!   & result(output)
!  implicit none
!  
!  integer, intent(in)           :: input(:)
!  integer, intent(in)           :: comparison
!  logical, intent(in), optional :: mask(:)
!  integer, intent(in), optional :: default
!  logical, intent(in), optional :: sorted
!  integer                       :: output
!  
!  integer :: i
!  
!  i = last(input,less_than_or_equal,mask,default,sorted)
!  if (present(default)) then
!    if (i==default) then
!      output = default
!    elseif (input(i)/=comparison) then
!      output = default
!    else
!      output = i
!    endif
!  else
!    if (input(i)==comparison) then
!      output = i
!    else
!      call print_line(ERROR//': value not found.')
!      call err()
!    endif
!  endif
!contains
!  ! Lambda for comparing the input to the comparison.
!  ! less_than_or_equal_comparison(input) is equivalent to
!  !    input<=comparison.
!  ! Captures comparison.
!  function less_than_or_equal(input) result(output)
!    implicit none
!    
!    class(*), intent(in) :: input
!    logical              :: output
!    
!    select type(input); type is(integer)
!      output = input<=comparison
!    end select
!  end function
!end function
!
!function last_equivalent_ComparisonLambda(input,comparison,equality, &
!   & inequality,mask,default,sorted) result(output)
!  implicit none
!  
!  class(*), intent(in)                  :: input(:)
!  class(*), intent(in)                  :: comparison
!  procedure(ComparisonLambda)           :: equality
!  procedure(ComparisonLambda), optional :: inequality
!  logical,  intent(in),        optional :: mask(:)
!  integer,  intent(in),        optional :: default
!  logical,  intent(in),        optional :: sorted
!  integer                               :: output
!  
!  integer :: i
!  
!  if (present(sorted) .and. .not. present(inequality)) then
!    if (.not. sorted) then
!      call print_line(CODE_ERROR//': sorted=true, but inequality not &
!         &present.')
!      call err()
!    endif
!  endif
!  
!  if (present(inequality)) then
!    i = last(input,inequality_lambda,mask,default,sorted)
!    if (present(default)) then
!      if (i==default) then
!        output = default
!      elseif (.not.equality(input(i),comparison)) then
!        output = default
!      else
!        output = i
!      endif
!    else
!      if (equality(input(i),comparison)) then
!        output = i
!      else
!        call print_line(ERROR//': value not found.')
!        call err()
!      endif
!    endif
!  else
!    output = last(input,equality_lambda,mask,default)
!  endif
!contains
!  function inequality_lambda(input) result(output)
!    implicit none
!    
!    class(*), intent(in) :: input
!    logical              :: output
!    
!    output = inequality(input,comparison)
!  end function
!  
!  function equality_lambda(input) result(output)
!    implicit none
!    
!    class(*), intent(in) :: input
!    logical              :: output
!    
!    output = equality(input,comparison)
!  end function
!end function

! ----------------------------------------------------------------------
! Operates with a lambda on every element on a list.
! e.g. if list=[2,3,7] then after operate(list,add_one) list=[3,4,8].
! ----------------------------------------------------------------------
subroutine operate_OperationLambda(input,lambda)
  class(*), intent(inout)    :: input(:)
  procedure(OperationLambda) :: lambda
  
  integer :: i
  
  do i=1,size(input)
    call lambda(input(i))
  enddo
end subroutine

! ----------------------------------------------------------------------
! Applies a logical lambda to a list.
! e.g. map([2,3,7,1], less_than_3) returns [T,F,F,T].
!
! A list of elements for which the lambda returns true can be returned using
!    the pack intrinsic.
! e.g. if list=[2,3,7,1], pack(list,map(list,less_than_3)) returns [2,1].
! ----------------------------------------------------------------------
function map_LogicalLambda(input,lambda) result(output)
  implicit none
  
  class(*), intent(in)     :: input(:)
  procedure(LogicalLambda) :: lambda
  logical, allocatable     :: output(:)
  
  integer :: i
  
  output = [( lambda(input(i)), i=1, size(input) )]
end function

function map_ComparisonLambda(input,lambda,comparison) result(output)
  implicit none
  
  class(*), intent(in)        :: input(:)
  procedure(ComparisonLambda) :: lambda
  class(*), intent(in)        :: comparison
  logical, allocatable        :: output(:)
  
  output = map(input,logical_lambda)
contains
  function logical_lambda(input) result(output)
    implicit none
    
    class(*), intent(in) :: input
    logical              :: output
    
    output = lambda(input,comparison)
  end function
end function

! ----------------------------------------------------------------------
! Works as the intrinsic count function, returning the number of true
!    elements in the input array.
! The _LogicalLambda variant counts the number of array elements for
!    which the lambda returns true.
! The _ComparisonLambda variant counts the number of array elements which
!    compare to true with the given comparison element.
! ----------------------------------------------------------------------
function count_LogicalLambda(input,lambda) result(output)
  implicit none
  
  class(*), intent(in)     :: input(:)
  procedure(LogicalLambda) :: lambda
  integer                  :: output
  
  output = count(map(input,lambda))
end function

function count_ComparisonLambda(input,lambda,comparison) result(output)
  implicit none
  
  class(*), intent(in)        :: input(:)
  procedure(ComparisonLambda) :: lambda
  class(*), intent(in)        :: comparison
  integer                     :: output
  
  output = count(map(input,lambda,comparison))
end function

! ----------------------------------------------------------------------
! The _logicals variant returns the indices of all elements where the
!    value is .true..
! The _lambda variant is equivalent, but uses a lambda to determine whether
!    each input value is true or false.
! e.g. filter([1,8,2,7],is_even) returns [2,3].
!
! To get the elements rather than their ids, the list should be indexed
!    by the output of filter.
! e.g. if list=[1,8,2,7] then list(filter(list,is_even)) returns [8,2].
!
! list(filter(list,lambda)) is equivalent to pack(list,map(list,lamba)).
! ----------------------------------------------------------------------
function filter_logicals(input,mask) result(output)
  implicit none
  
  logical, intent(in)           :: input(:)
  logical, intent(in), optional :: mask(:)
  integer, allocatable          :: output(:)
  
  integer :: i
  
  if (present(mask)) then
    if (size(mask)/=size(input)) then
      call print_line(CODE_ERROR//': Input and Mask of different sizes.')
      call err()
    endif
  endif
  
  ! Make an array [1,2,3,...,size(input)].
  output = [(i,i=1,size(input))]
  if (present(mask)) then
    ! Return only those elements where input=true and mask=true.
    output = pack(output, mask=input.and.mask)
  else
    ! Return only those elements where input=true.
    output = pack(output, mask=input)
  endif
end function

function filter_LogicalLambda(input,lambda,mask) result(output)
  implicit none
  
  class(*), intent(in)           :: input(:)
  procedure(LogicalLambda)       :: lambda
  logical,  intent(in), optional :: mask(:)
  integer, allocatable           :: output(:)
  
  output = filter(map(input,lambda),mask=mask)
end function

function filter_ComparisonLambda(input,lambda,comparison,mask) result(output)
  implicit none
  
  class(*), intent(in)           :: input(:)
  procedure(ComparisonLambda)    :: lambda
  class(*), intent(in)           :: comparison
  logical,  intent(in), optional :: mask(:)
  integer, allocatable           :: output(:)
  
  output = filter(input,logical_lambda,mask)
contains
  function logical_lambda(input) result(output)
    implicit none
    
    class(*), intent(in) :: input
    logical              :: output
    
    output = lambda(input,comparison)
  end function
end function

! ----------------------------------------------------------------------
! Finds the value in a list which compares well with all other values.
! e.g. if the comparison is less than (<), finds the minimum value.
!    similarly, if lambda= greater than (>), finds the maximum value.
! Takes a lambda to compare values.
! e.g. locate([7,1,5], less_than) returns 2.
! If mask is present, only considers values for which mask is true.
! ----------------------------------------------------------------------
! N.B. assumes that the input list can be totally ordered using the comparison.
function locate_ComparisonLambda(input,lambda,mask) result(output)
  implicit none
  
  class(*), intent(in)          :: input(:)
  procedure(ComparisonLambda)   :: lambda
  logical, intent(in), optional :: mask(:)
  integer                       :: output
  
  integer :: first_unmasked
  
  integer :: i
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': Calling locate on an empty list.')
    call err()
  endif
  
  if (present(mask)) then
    if (size(mask)/=size(input)) then
      call print_line(CODE_ERROR//': Input and Mask of different sizes.')
      call err()
    endif
  endif
  
  if (present(mask)) then
    first_unmasked = first(mask)
    if (first_unmasked==0) then
      call print_line(CODE_ERROR//': Calling locate on entirely masked list.')
    endif
  else
    first_unmasked = 1
  endif
  
  output = first_unmasked
  do i=first_unmasked+1,size(input)
    if (present(mask)) then
      if (lambda(input(i),input(output)) .and. mask(i)) then
        output = i
      endif
    else
      if (lambda(input(i),input(output))) then
        output = i
      endif
    endif
  enddo
end function

! ----------------------------------------------------------------------
! The _integer variant returns .true. if the list is sorted in ascending order,
!    and .false. otherwise.
!
! The _lambda variant does the same but with a list of arbitrary type, using
!    the provided lambda for comparison. The function returns true if
!    lambda(input(i),input(i+1)) = true and lambda(input(i+1),input(i)) = false
!    for every element.
! The lambda should be < or <= for checking the list is in ascending order,
!    or > or >= for checking the list is in descending order.
! ----------------------------------------------------------------------
function is_sorted_integers(input) result(output)
  implicit none
  
  integer, intent(in) :: input(:)
  logical             :: output
  
  if (size(input)<2) then
    output = .true.
  else
    output = all(input(:size(input)-1) <= input(2:))
  endif
end function

function is_sorted_ComparisonLambda(input,lambda) result(output)
  implicit none
  
  class(*), intent(in)        :: input(:)
  procedure(ComparisonLambda) :: lambda
  logical                     :: output
  
  integer :: i
  
  if (size(input)<2) then
    output = .true.
  else
    do i=1,size(input)-1
      if (       lambda(input(i+1), input(i)  ) .and. &
         & .not. lambda(input(i),   input(i+1)) ) then
        output = .false.
        return
      endif
    enddo
  endif
end function

! ----------------------------------------------------------------------
! The _integer variant sorts a list of integers in ascending order,
!    and returns the ids of the sorted elements. 
! e.g. sort([1,3,2,1]) returns [1,4,3,2].
!
! The _lambda variant takes a list of items, and a lambda to compare
!    elements of the list.
! The order of the output depends on the lambda. e.g. if the lambda is < then
!    the list will be sorted into ascending order.
! e.g. sort(['a','c','b'],first_alphabetically) returns [1,3,2].
!
! In order to get the sorted list, the output of sort should be used to
!    index the original list.
! e.g. if list = [1,3,2,1] then list(sort(list)) will return [1,1,2,3].
! e.g. if list = ['c','a','c','b'] then list(sort(list,first_alphabetically))
!    will return ['a','b','c','c'].
!
! Both variants check if the list is sorted before sorting, since
!    checking the list is an O(n) operation, whereas sorting the list is
!    an O(n^2) operation.
! ----------------------------------------------------------------------
function sort_integers(input) result(output)
  implicit none
  
  integer, intent(in)  :: input(:)
  integer, allocatable :: output(:)
  
  logical, allocatable :: sorted(:)
  
  integer :: i,j,ialloc
  
  if (is_sorted(input)) then
    output = [(i,i=1,size(input))]
  else
    allocate( sorted(size(input)), &
            & output(size(input)), &
            & stat=ialloc); call err(ialloc)
    sorted = .false.
    
    do i=1,size(input)
      j = minloc(input,1,mask=.not. sorted)
      sorted(j) = .true.
      output(i) = j
    enddo
  endif
end function

function sort_ComparisonLambda(input,lambda) result(output)
  implicit none
  
  class(*), intent(in)        :: input(:)
  procedure(ComparisonLambda) :: lambda
  integer, allocatable        :: output(:)
  
  logical, allocatable :: sorted(:)
  
  integer :: i,j,ialloc
  
  if (is_sorted(input,lambda)) then
    output = [(i,i=1,size(input))]
  else
    allocate( sorted(size(input)), &
            & output(size(input)), &
            & stat=ialloc); call err(ialloc)
    sorted = .false.
    
    do i=1,size(input)
      j = locate(input,lambda,mask=.not. sorted)
      sorted(j) = .true.
      output(i) = j
    enddo
  endif
end function

! ----------------------------------------------------------------------
! The _integers variant returns the ids at which the first copy of each
!    integer can be found.
! e.g. set([1,1,3,1,2,3]) returns [1,3,5].
!
! To get the deduplicated list, the result of this function should be used
!    to index the original list.
! e.g. if list=[1,1,3,1,2,3] then list(set(list)) returns [1,3,2].
!
! The _lambda variant is the same but using the provided comparison lambda
!    instead of integer comparison.
!
! If mask is provided, then only elements where mask=.true. will be considered.
! ----------------------------------------------------------------------

function set_integers(input,mask) result(output)
  implicit none
  
  integer, intent(in)           :: input(:)
  logical, intent(in), optional :: mask(:)
  integer, allocatable          :: output(:)
  
  integer :: no_unique_elements
  
  integer :: i,ialloc
  
  if (present(mask)) then
    if (size(mask)/=size(input)) then
      call print_line(CODE_ERROR//': Input and Mask of different sizes.')
      call err()
    endif
  endif
  
  allocate(output(size(input)), stat=ialloc); call err(ialloc)
  no_unique_elements = 0
  do i=1,size(input)
    if (present(mask)) then
      if (.not. mask(i)) then
        cycle
      endif
    endif
    if (.not. any(input(output(:no_unique_elements))==input(i))) then
      no_unique_elements = no_unique_elements+1
      output(no_unique_elements) = i
    endif
  enddo
  
  output = output(:no_unique_elements)
end function

function set_ComparisonLambda(input,lambda,mask) result(output)
  implicit none
  
  class(*), intent(in)          :: input(:)
  procedure(ComparisonLambda)   :: lambda
  logical, intent(in), optional :: mask(:)
  integer, allocatable          :: output(:)
  
  integer :: no_unique_elements
  
  integer :: i,j,ialloc
  
  if (present(mask)) then
    if (size(mask)/=size(input)) then
      call print_line(CODE_ERROR//': Input and Mask of different sizes.')
      call err()
    endif
  endif
  
  allocate(output(size(input)), stat=ialloc); call err(ialloc)
  no_unique_elements = 0
  do_i : do i=1,size(input)
    if (present(mask)) then
      if (.not. mask(i)) then
        cycle
      endif
    endif
    
    do j=1,no_unique_elements
      if (lambda(input(i),input(output(j)))) then
        cycle do_i
      endif
    enddo
    
    no_unique_elements = no_unique_elements + 1
    output(no_unique_elements) = i
  enddo do_i
  
  output = output(:no_unique_elements)
end function
end module

! ======================================================================
! Example use cases of procedures in logic_module.
! ======================================================================
module logic_example_module
  use io_basic_module
  
  use logic_module
  implicit none
contains

! --------------------------------------------------
! A subroutine to demonstrate logic procedures.
! --------------------------------------------------
subroutine logic_example()
  implicit none
  
  logical, allocatable      :: logicals(:)
  integer, allocatable      :: integers(:)
  type(String), allocatable :: strings(:)
  
  integer :: val
  
  logicals = [ .false., .true., .true., .false., .true., .false. ]
  
  ! Print out the array of logicals.
  call print_line('')
  call print_line('logicals should   = F T T F T F')
  call print_line('logicals actually = '//logicals)
  
  ! first(logicals)=2 because logicals(2) is the first true element.
  call print_line('')
  call print_line('first(logicals) should   = 2')
  call print_line('first(logicals) actually = '//first(logicals))
  
  integers = [ 0, 1, 1, 0, 0 ]
  
  ! Print out the array of integers.
  call print_line('')
  call print_line('integers should   =  0  1  1  0  0')
  call print_line('integers actually = '//integers)
  
  val = 0
  ! first(integers,equals)=1 when val=0, because integers(2) is
  !    the first element which equals val.
  call print_line('')
  call print_line('first(integers==0) should   = 1')
  call print_line('first(integers==0) actually = '//first(integers,equals))
  
  val = 1
  ! first(integers,equals)=2 when val=1, because integers(2) is
  !    the first element which equals val.
  call print_line('')
  call print_line('first(integers==1) should   = 2')
  call print_line('first(integers==1) actually = '//first(integers,equals))
  
  ! map(integers,equals)=[F,T,T,F,F] when val=1.
  call print_line('')
  call print_line('map(integers==1) should   = F T T F F')
  call print_line('map(integers==1) actually = '// &
                 & map(integers,equals))
  
  ! filter(integers,equals)=[2,3] when val=1.
  call print_line('')
  call print_line('filter(integers==1) should   =  2  3')
  call print_line('filter(integers==1) actually = '// &
                 & filter(integers,equals))
  
  ! integers(filter(integers,equals)) = [1,1] when val=1.
  call print_line('')
  call print_line('integers(filter(integers==1)) should   =  1  1')
  call print_line('integers(filter(integers==1)) actually = '// &
                 & integers(filter(integers,equals)))
  
  ! pack(integers,map(integers,equals)) = [1,1] when val=1.
  call print_line('')
  call print_line('pack(integers,map(integers==1)) should   =  1  1')
  call print_line('pack(integers,map(integers==1)) actually = '// &
                 & pack(integers,map(integers,equals)))
  
  ! sort(integers) = [1,4,5,2,3] because elements 1, 4 and 5 are zero, so are
  !    sorted to the beginning, and elements 2 and 3 are 1.
  call print_line('')
  call print_line('sort(integers) should             =  1  4  5  2  3')
  call print_line('sort(integers) actually           = '//sort(integers))
  
  call print_line('')
  call print_line('integers(sort(integers)) should   =  0  0  0  1  1')
  call print_line('integers(sort(integers)) actually = '// &
                 & integers(sort(integers)))
  
  strings = [ str('3'), str('2'), str('5') ]
  
  ! Print out the array of strings.
  call print_line('')
  call print_line('strings should   = 3 2 5')
  call print_line('strings actually = '//join(strings))
  
  ! sort(strings,string_to_int) = [2,1,3] because the lambda string_to_int
  !    converts the string array ['3','2','5'] to the int array [3,2,5],
  !    and this is then sorted.
  ! strings(2)=2 < strings(1)=3 < strings(3)=5, so the output is [2,1,3].
  call print_line('')
  call print_line('sort(strings,less_than_strings) should   =  2  1  3')
  call print_line( 'sort(strings,less_than_strings) actually = '// &
                 & sort(strings, less_than_strings))
  
  integers = [ 0, 1, 1, 0, 0 ]
  call print_line('')
  call print_line('integers should   =  0  1  1  0  0')
  call print_line('integers actually = '//integers)
  
  ! set(integers) = [1,2].
  ! integers(set(integers)) = [0,1].
  call print_line('')
  call print_line('set(integers) should   =  1  2')
  call print_line('set(integers) actually = '//set(integers))
  call print_line('integers(set(integers)) should   =  0  1')
  call print_line( 'integers(set(integers)) actually = '// &
                 & integers(set(integers)))
  
  strings = [ str('3'), str('1'), str('3'), str('2') ]
  call print_line('')
  call print_line('strings should   = 3 1 3 2')
  call print_line('strings actually = '//join(strings))
  
  ! set(strings) = [1,2,4].
  ! strings(set(strings)) = ['3','1','2'].
  call print_line('')
  call print_line('set(strings,equals_strings) should   =  1  2  4')
  call print_line( 'set(strings,equals_strings) actually = '// &
                 & set(strings,equals_strings))
  call print_line( 'strings(set(strings,equals_strings)) should   &
     &= 3 1 2')
  call print_line( 'strings(set(strings,equals_strings)) actually &
     &= '// join(strings(set(strings,equals_strings))))
contains
! Lambda functions for passing to logic functions.
  
  ! Takes an integer, and returns whether or not it equals val,
  !    where val is defined in the enclosing function.
  ! N.B. val is "captured" in the sense of lambda functions.
  function equals(input) result(output)
    implicit none
    
    class(*), intent(in) :: input
    logical              :: output
    
    select type(input)
      type is(integer)
        output = input==val
      class default
        call err()
    end select
  end function

  ! Takes two strings, and compares them (<) as integers.
  function less_than_strings(this,that) result(output)
    implicit none
    
    class(*), intent(in) :: this
    class(*), intent(in) :: that
    logical              :: output
    
    select type(this)
      type is(String)
        
        select type(that)
          type is(String)
            output = int(this) < int(that)
          class default
            call err()
        end select
        
      class default
        call err()
    end select
  end function
  
  ! Takes two strings, and compares them (==) as integers.
  function equals_strings(this,that) result(output)
    implicit none
    
    class(*), intent(in) :: this
    class(*), intent(in) :: that
    logical              :: output
    
    select type(this)
      type is(String)
        
        select type(that)
          type is(String)
            output = int(this) == int(that)
          class default
            call err()
        end select
        
      class default
        call err()
    end select
  end function
end subroutine
end module
