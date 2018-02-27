! ======================================================================
! A number of algorithms to do with logic and sorting.
! ======================================================================
! See the example module for some use cases.
module logic_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  private
  
  public :: first
  public :: last
  public :: evaluate
  public :: filter
  public :: locate
  public :: sort
  
  interface first
    module procedure first_logicals
    module procedure first_lambda
  end interface
  
  interface last
    module procedure last_logicals
    module procedure last_lambda
  end interface
  
  interface map
    module procedure map_lambda
  end interface
  
  interface evaluate
    module procedure evaluate_lambda
  end interface
  
  interface filter
    module procedure filter_logicals
    module procedure filter_lambda
  end interface
  
  interface locate
    module procedure locate_lambda
  end interface
  
  interface sort
    module procedure sort_integers
    module procedure sort_lambda
  end interface
  
  interface
    function logical_lambda(input) result(output)
      implicit none
      
      class(*), intent(in) :: input
      logical              :: output
    end function
    
    function compare_lambda(this,that) result(output)
      implicit none
      
      class(*), intent(in) :: this
      class(*), intent(in) :: that
      logical              :: output
    end function
    
    subroutine operate_lambda(this)
      implicit none
      
      class(*), intent(inout) :: this
    end subroutine
  end interface
contains

! ----------------------------------------------------------------------
! The _logicals variants return the id of the first/last .true. in a 
!    list of logicals. e.g. first([f,f,t,f]) returns 3.
! The _lambda variants do the same, but for the lambda applied to the input,
!    e.g. first([1,1,2,1],equals_2) returns 3.
! In both cases, returns 0 if all values are .false..
! If mask is present, only returns values where mask is .true..
! ----------------------------------------------------------------------

function first_logicals(input,mask) result(output)
  implicit none
  
  logical, intent(in)           :: input(:)
  logical, intent(in), optional :: mask(:)
  integer                       :: output
  
  integer :: i
  
  if (present(mask)) then
    if (size(input)/=size(mask)) then
      call print_line(CODE_ERROR//': Input and Mask of different sizes.')
      call err()
    endif
  endif
  
  do i=1,size(input)
    if (present(mask)) then
      if (input(i).and.mask(i)) then
        output = i
        return
      endif
    else
      if (input(i)) then
        output = i
        return
      endif
    endif
  enddo
  
  output = 0
end function

function first_lambda(input,lambda,mask) result(output)
  implicit none
  
  class(*), intent(in)           :: input(:)
  procedure(logical_lambda)      :: lambda
  logical,  intent(in), optional :: mask(:)
  integer                        :: output
  
  integer :: i
  
  if (present(mask)) then
    if (size(mask)/=size(input)) then
      call print_line(CODE_ERROR//': Input and Mask of different sizes.')
      call err()
    endif
  endif
  
  output = 0
  
  do i=1,size(input)
    if (present(mask)) then
      if (lambda(input(i)) .and. mask(i)) then
        output = i
        exit
      endif
    else
      if (lambda(input(i))) then
        output = i
        exit
      endif
    endif
  enddo
end function

function last_logicals(input,mask) result(output)
  implicit none
  
  logical, intent(in)           :: input(:)
  logical, intent(in), optional :: mask(:)
  integer                       :: output
  
  integer :: i
  
  if (present(mask)) then
    if (size(input)/=size(mask)) then
      call print_line(CODE_ERROR//': Input and Mask of different sizes.')
      call err()
    endif
  endif
  
  do i=size(input),1,-1
    if (present(mask)) then
      if (input(i).and.mask(i)) then
        output = i
        return
      endif
    else
      if (input(i)) then
        output = i
        return
      endif
    endif
  enddo
  
  output = 0
end function

function last_lambda(input,lambda,mask) result(output)
  implicit none
  
  class(*), intent(in)           :: input(:)
  procedure(logical_lambda)      :: lambda
  logical,  intent(in), optional :: mask(:)
  integer                        :: output
  
  integer :: i
  
  if (present(mask)) then
    if (size(mask)/=size(input)) then
      call print_line(CODE_ERROR//': Input and Mask of different sizes.')
      call err()
    endif
  endif
  
  output = 0
  
  do i=size(input),1,-1
    if (present(mask)) then
      if (lambda(input(i)) .and. mask(i)) then
        output = i
        exit
      endif
    else
      if (lambda(input(i))) then
        output = i
        exit
      endif
    endif
  enddo
end function

! ----------------------------------------------------------------------
! Operates with a lambda on every element on a list.
! e.g. if list=[2,3,7] then after map(list,add_one) list=[3,4,8].
! ----------------------------------------------------------------------
subroutine map_lambda(input,lambda)
  class(*), intent(inout)   :: input(:)
  procedure(operate_lambda) :: lambda
  
  integer :: i
  
  do i=1,size(input)
    call lambda(input(i))
  enddo
end subroutine

! ----------------------------------------------------------------------
! Applies a logical lambda to a list.
! e.g. evaluate([2,3,7,1], less_than_3) returns [T,F,F,T].
!
! A list of elements for which the lambda returns true can be returned using
!    the pack intrinsic.
! e.g. if list=[2,3,7,1], pack(list,evaluate(list,less_than_3)) returns [2,1].
! ----------------------------------------------------------------------
function evaluate_lambda(input,lambda) result(output)
  implicit none
  
  class(*), intent(in)      :: input(:)
  procedure(logical_lambda) :: lambda
  logical, allocatable      :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(input)), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    output(i) = lambda(input(i))
  enddo
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
! list(filter(list,lambda)) is equivalent to pack(list,evaluate(list,lamba)).
! ----------------------------------------------------------------------
function filter_logicals(input) result(output)
  implicit none
  
  logical, intent(in)  :: input(:)
  integer, allocatable :: output(:)
  
  integer :: i
  
  ! Make an array [1,2,3,...,size(input)].
  output = [(i,i=1,size(input))]
  ! Return only those elements where input=true.
  output = pack(output, mask=input)
end function

function filter_lambda(input,lambda) result(output)
  implicit none
  
  class(*), intent(in)      :: input(:)
  procedure(logical_lambda) :: lambda
  integer, allocatable      :: output(:)
  
  output = filter(evaluate(input,lambda))
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
function locate_lambda(input,lambda,mask) result(output)
  implicit none
  
  class(*), intent(in)          :: input(:)
  procedure(compare_lambda)     :: lambda
  logical, intent(in), optional :: mask(:)
  integer                       :: output
  
  integer :: first_unmasked
  
  integer :: i
  
  if (size(input)==0) then
    call print_line(CODE_ERROR//': Calling locate on an empty list.')
    call err()
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
! ----------------------------------------------------------------------
function sort_integers(input) result(output)
  implicit none
  
  integer, intent(in)  :: input(:)
  integer, allocatable :: output(:)
  
  logical, allocatable :: sorted(:)
  
  integer :: i,j,ialloc
  
  allocate( sorted(size(input)), &
          & output(size(input)), &
          & stat=ialloc); call err(ialloc)
  sorted = .false.
  
  do i=1,size(input)
    j = minloc(input,1,mask=.not. sorted)
    sorted(j) = .true.
    output(i) = j
  enddo
end function

function sort_lambda(input,lambda) result(output)
  implicit none
  
  class(*), intent(in)      :: input(:)
  procedure(compare_lambda) :: lambda
  integer, allocatable      :: output(:)
  
  logical, allocatable :: sorted(:)
  
  integer :: i,j,ialloc
  
  allocate( sorted(size(input)), &
          & output(size(input)), &
          & stat=ialloc); call err(ialloc)
  sorted = .false.
  
  do i=1,size(input)
    j = locate(input,lambda,mask=.not. sorted)
    sorted(j) = .true.
    output(i) = j
  enddo
end function
end module

! ======================================================================
! Example use cases of procedures in logic_module.
! ======================================================================
module logic_example_module
  use string_module
  use io_module
  
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
  
  ! last(logicals)=5 because logicals(5) is the last true element.
  call print_line('')
  call print_line('last(logicals) should   = 5')
  call print_line('last(logicals) actually = '//last(logicals))
  
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
  
  ! last(integers,equals)=5 when val=0, because integers(2) is
  !    the last element which equals val.
  call print_line('')
  call print_line('last(integers==0) should   = 5')
  call print_line('last(integers==0) actually = '//last(integers,equals))
  
  val = 1
  ! first(integers,equals)=2 when val=1, because integers(2) is
  !    the first element which equals val.
  call print_line('')
  call print_line('first(integers==1) should   = 2')
  call print_line('first(integers==1) actually = '//first(integers,equals))
  
  ! last(integers,equals)=3 when val=1, because integers(2) is
  !    the last element which equals val.
  call print_line('')
  call print_line('last(integers==1) should   = 3')
  call print_line('last(integers==1) actually = '//last(integers,equals))
  
  ! evaluate(integers,equals)=[F,T,T,F,F] when val=1.
  call print_line('')
  call print_line('evaluate(integers==1) should   = F T T F F')
  call print_line('evaluate(integers==1) actually = '// &
                 & evaluate(integers,equals))
  
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
  
  ! pack(integers,evaluate(integers,equals)) = [1,1] when val=1.
  call print_line('')
  call print_line('pack(integers,evaluate(integers==1)) should   =  1  1')
  call print_line('pack(integers,evaluate(integers==1)) actually = '// &
                 & pack(integers,evaluate(integers,equals)))
  
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
  call print_line('sort(strings,string_to_int) should   =  2  1  3')
  call print_line( 'sort(strings,string_to_int) actually = '// &
                 & sort(strings, less_than_strings))
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

  ! Takes two strings, and compares them as integers.
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
end subroutine
end module