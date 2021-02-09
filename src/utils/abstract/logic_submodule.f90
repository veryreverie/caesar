submodule (caesar_logic_module) caesar_logic_submodule
contains
module procedure lazy_and
  if (.not. this) then
    output = .false.
  elseif (present(that)) then
    output = that
  else
    call print_line(ERROR//': First argument to lazy_and is .true. and &
       &second argument is not present.')
    call err()
  endif
end procedure

module procedure lazy_or
  if (this) then
    output = .true.
  elseif (present(that)) then
    output = that
  else
    call print_line(ERROR//': First argument to lazy_or is .false. and &
       &second argument is not present.')
    call err()
  endif
end procedure

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

module procedure first_logicals
  integer :: lower_bound_false
  integer :: lower_bound_masked
  integer :: upper_bound_masked
  integer :: upper_bound_true
  
  integer :: i,j
  
  ! Check that mask is consistent with input.
  if (present(mask)) then
    if (size(input)/=size(mask)) then
      call print_line(CODE_ERROR//': Input and Mask of different sizes.')
      call err()
    endif
  endif
  
  ! Find the first true.
  lower_bound_false = 0
  lower_bound_masked = 0
  upper_bound_masked = size(input)+1
  upper_bound_true = size(input)+1
  if (set_default(sorted, .false.)) then
    ! If the input is sorted, use a bisection method to find the first true.
    if (present(mask)) then
      do while(upper_bound_masked-lower_bound_masked>1)
        ! Bisect the input between the last known false and first known true.
        i = (lower_bound_masked+upper_bound_masked)/2
        ! Find the last unmasked value in [lower_bound_masked+1:i].
        ! If no such values exist,
        !    find the first unmasked value in [i+1:upper_bound_masked-1]
        ! If no such values exist, exit the loop.
        j = lower_bound_masked + last(mask(lower_bound_masked+1:i), default=0)
        if (j==lower_bound_masked) then
          j = i + first(mask(i+1:upper_bound_masked-1), default=0)
          if (j==i) then
            exit
          endif
        endif
        
        ! Evaluate input(j), and update the bounds accordingly.
        if (input(j)) then
          upper_bound_true = j
          if (j<=i) then
            upper_bound_masked = j
          else
            exit
          endif
        else
          lower_bound_false = i
          lower_bound_masked = max(i,j)
        endif
      enddo
    else
      do while(upper_bound_true-lower_bound_false>1)
        ! Bisect the input between the last known false and first known true.
        i = (lower_bound_false+upper_bound_true)/2
        ! Evaluate input(i), and update the bounds accordingly.
        if (input(i)) then
          upper_bound_true = i
        else
          lower_bound_false = i
        endif
      enddo
    endif
  else
    ! If the input is not sorted, simply iterate through the input.
    if (present(mask)) then
      do i=1,size(input)
        if (mask(i)) then
          if (input(i)) then
            upper_bound_true = i
            exit
          endif
        endif
      enddo
    else
      do i=1,size(input)
        if (input(i)) then
          upper_bound_true = i
          exit
        endif
      enddo
    endif
  endif
  
  ! Set the output to the first true, or the default value if present.
  ! Throw an error if no true is found and no default is set.
  if (upper_bound_true<=size(input)) then
    output = upper_bound_true
  elseif (present(default)) then
    output = default
  else
    call print_line(ERROR//': value not found.')
    call err()
  endif
end procedure

! WORKAROUND: requires explicit interface because of
!    polymorphic procedure argument (gfortran 10 bug).
function first_LogicalLambda(input,lambda,mask,default,sorted) &
   & result(output) 
  class(*), intent(in)           :: input(:)
  procedure(LogicalLambda)       :: lambda
  logical,  intent(in), optional :: mask(:)
  integer,  intent(in), optional :: default
  logical,  intent(in), optional :: sorted
  integer                        :: output
  
  integer :: lower_bound_false
  integer :: lower_bound_masked
  integer :: upper_bound_masked
  integer :: upper_bound_true
  
  integer :: i,j
  
  ! Check that mask is consistent with input.
  if (present(mask)) then
    if (size(mask)/=size(input)) then
      call print_line(CODE_ERROR//': Input and Mask of different sizes.')
      call err()
    endif
  endif
  
  lower_bound_false = 0
  lower_bound_masked = 0
  upper_bound_masked = size(input)+1
  upper_bound_true = size(input)+1
  if (set_default(sorted, .false.)) then
    ! If the input is sorted, use a bisection method to find the first true.
    if (present(mask)) then
      do while(upper_bound_masked-lower_bound_masked>1)
        ! Bisect the input between the last known false and first known true.
        i = (lower_bound_masked+upper_bound_masked)/2
        ! Find the last unmasked value in [lower_bound_masked+1:i].
        ! If no such values exist,
        !    find the first unmasked value in [i+1:upper_bound_masked-1]
        ! If no such values exist, exit the loop.
        j = lower_bound_masked + last(mask(lower_bound_masked+1:i), default=0)
        if (j==lower_bound_masked) then
          j = i + first(mask(i+1:upper_bound_masked-1), default=0)
          if (j==i) then
            exit
          endif
        endif
        
        ! Evaluate lambda(input(j)), and update the bounds accordingly.
        if (lambda(input(j))) then
          upper_bound_true = j
          if (j<=i) then
            upper_bound_masked = j
          else
            exit
          endif
        else
          lower_bound_false = i
          lower_bound_masked = max(i,j)
        endif
      enddo
    else
      do while(upper_bound_true-lower_bound_false>1)
        ! Bisect the input between the last known false and first known true.
        i = (lower_bound_false+upper_bound_true)/2
        ! Evaluate input(i), and update the bounds accordingly.
        if (lambda(input(i))) then
          upper_bound_true = i
        else
          lower_bound_false = i
        endif
      enddo
    endif
  else
    ! If the input is not sorted, simply iterate through the input.
    if (present(mask)) then
      do i=1,size(input)
        if (mask(i)) then
          if (lambda(input(i))) then
            upper_bound_true = i
            exit
          endif
        endif
      enddo
    else
      do i=1,size(input)
        if (lambda(input(i))) then
          upper_bound_true = i
          exit
        endif
      enddo
    endif
  endif
  
  ! Set the output to the first true, or the default value if present.
  ! Throw an error if no true is found and no default is set.
  if (upper_bound_true<=size(input)) then
    output = upper_bound_true
  elseif (present(default)) then
    output = default
  else
    call print_line(ERROR//': value not found.')
    call err()
  endif
end function

module procedure last_logicals
  if (present(mask)) then
    output = first( input   = input(size(input):1:-1), &
                  & mask    = mask(size(mask):1:-1),   &
                  & default = 0,                       &
                  & sorted  = sorted                   )
  else
    output = first( input   = input(size(input):1:-1), &
                  & default = 0,                       &
                  & sorted  = sorted                   )
  endif
  
  if (output/=0) then
    output = size(input)+1-output
  elseif (present(default)) then
    output = default
  else
    call print_line(ERROR//': value not found.')
    call err()
  endif
end procedure

module procedure last_LogicalLambda
  if (present(mask)) then
    output = first( input   = input(size(input):1:-1), &
                  & lambda  = lambda,                  &
                  & mask    = mask(size(mask):1:-1),   &
                  & default = 0,                       &
                  & sorted  = sorted                   )
  else
    output = first( input   = input(size(input):1:-1), &
                  & lambda  = lambda,                  &
                  & default = 0,                       &
                  & sorted  = sorted                   )
  endif
  
  if (output/=0) then
    output = size(input)+1-output
  elseif (present(default)) then
    output = default
  else
    call print_line(ERROR//': value not found.')
    call err()
  endif
end procedure

module procedure first_equivalent_integers
  integer :: lower_bound_smaller
  integer :: lower_bound_masked
  integer :: upper_bound_masked
  integer :: upper_bound_larger_or_equal
  
  integer :: i,j
  
  ! Check that mask is consistent with input.
  if (present(mask)) then
    if (size(input)/=size(mask)) then
      call print_line(CODE_ERROR//': Input and Mask of different sizes.')
      call err()
    endif
  endif
  
  ! Find the first larger_or_equal.
  lower_bound_smaller = 0
  lower_bound_masked = 0
  upper_bound_masked = size(input)+1
  upper_bound_larger_or_equal = size(input)+1
  if (set_default(sorted, .false.)) then
    ! If the input is sorted, use a bisection method to find the first
    !    value in input which is >= comparison.
    if (present(mask)) then
      do while(upper_bound_masked-lower_bound_masked>1)
        ! Bisect the input between the last known smaller value
        !    and first known larger_or_equal value.
        i = (lower_bound_masked+upper_bound_masked)/2
        ! Find the last unmasked value in [lower_bound_masked+1:i].
        ! If no such values exist,
        !    find the first unmasked value in [i+1:upper_bound_masked-1]
        ! If no such values exist, exit the loop.
        j = lower_bound_masked + last(mask(lower_bound_masked+1:i), default=0)
        if (j==lower_bound_masked) then
          j = i + first(mask(i+1:upper_bound_masked-1), default=0)
          if (j==i) then
            exit
          endif
        endif
        
        ! Compare input(j) to comparison, and update the bounds accordingly.
        if (input(j)>=comparison) then
          upper_bound_larger_or_equal = j
          if (j<=i) then
            upper_bound_masked = j
          else
            exit
          endif
        else
          lower_bound_smaller = i
          lower_bound_masked = max(i,j)
        endif
      enddo
    else
      do while(upper_bound_larger_or_equal-lower_bound_smaller>1)
        ! Bisect the input between the last known smaller value
        !    and the first known larger_or_equal value.
        i = (lower_bound_smaller+upper_bound_larger_or_equal)/2
        ! Compare input(i) to comparison, and update the bounds accordingly.
        if (input(i)>=comparison) then
          upper_bound_larger_or_equal = i
        else
          lower_bound_smaller = i
        endif
      enddo
    endif
  else
    ! If the input is not sorted, simply iterate through the input.
    if (present(mask)) then
      do i=1,size(input)
        if (mask(i)) then
          if (input(i)==comparison) then
            output = i
            return
          endif
        endif
      enddo
    else
      do i=1,size(input)
        if (input(i)==comparison) then
          output = i
          return
        endif
      enddo
    endif
  endif
  
  ! Check if the first larger_or_equal value is equivalent to comparison.
  if (upper_bound_larger_or_equal<=size(input)) then
    if (input(upper_bound_larger_or_equal)==comparison) then
      output = upper_bound_larger_or_equal
      return
    endif
  endif
  
  ! No equivalent value has been found. Return the default value if set,
  !    or throw an error if not.
  if (present(default)) then
    output = default
  else
    call print_line(ERROR//': value not found.')
    call err()
  endif
end procedure

! WORKAROUND: requires explicit interface because of
!    polymorphic procedure argument (gfortran 10 bug).
function first_equivalent_ComparisonLambda(input,comparison,equality, &
   & greater_than_or_equal,mask,default,sorted) result(output) 
  class(*), intent(in)                  :: input(:)
  class(*), intent(in)                  :: comparison
  procedure(ComparisonLambda)           :: equality
  procedure(ComparisonLambda), optional :: greater_than_or_equal
  logical,  intent(in),        optional :: mask(:)
  integer,  intent(in),        optional :: default
  logical,  intent(in),        optional :: sorted
  integer                               :: output
  
  integer :: lower_bound_smaller
  integer :: lower_bound_masked
  integer :: upper_bound_masked
  integer :: upper_bound_larger_or_equal
  
  integer :: i,j
  
  ! Check that mask is consistent with input.
  if (present(mask)) then
    if (size(input)/=size(mask)) then
      call print_line(CODE_ERROR//': Input and Mask of different sizes.')
      call err()
    endif
  endif
  
  ! Check that greater_than_or_equal is present if sorted is true.
  if (set_default(sorted,.false.)) then
    if (.not. present(greater_than_or_equal)) then
      call print_line(CODE_ERROR//': sorted is true, but &
         &greater_than_or_equal is not present.')
      call err()
    endif
  endif
  
  ! Find the first larger_or_equal.
  lower_bound_smaller = 0
  lower_bound_masked = 0
  upper_bound_masked = size(input)+1
  upper_bound_larger_or_equal = size(input)+1
  if (set_default(sorted, .false.)) then
    ! If the input is sorted, use a bisection method to find the first
    !    value in input which is >= comparison.
    if (present(mask)) then
      do while(upper_bound_masked-lower_bound_masked>1)
        ! Bisect the input between the last known smaller value
        !    and first known larger_or_equal value.
        i = (lower_bound_masked+upper_bound_masked)/2
        ! Find the last unmasked value in [lower_bound_masked+1:i].
        ! If no such values exist,
        !    find the first unmasked value in [i+1:upper_bound_masked-1]
        ! If no such values exist, exit the loop.
        j = lower_bound_masked + last(mask(lower_bound_masked+1:i), default=0)
        if (j==lower_bound_masked) then
          j = i + first(mask(i+1:upper_bound_masked-1), default=0)
          if (j==i) then
            exit
          endif
        endif
        
        ! Compare input(j) to comparison, and update the bounds accordingly.
        if (greater_than_or_equal(input(j),comparison)) then
          upper_bound_larger_or_equal = j
          if (j<=i) then
            upper_bound_masked = j
          else
            exit
          endif
        else
          lower_bound_smaller = i
          lower_bound_masked = max(i,j)
        endif
      enddo
    else
      do while(upper_bound_larger_or_equal-lower_bound_smaller>1)
        ! Bisect the input between the last known smaller value
        !    and the first known larger_or_equal value.
        i = (lower_bound_smaller+upper_bound_larger_or_equal)/2
        ! Compare input(i) to comparison, and update the bounds accordingly.
        if (greater_than_or_equal(input(i),comparison)) then
          upper_bound_larger_or_equal = i
        else
          lower_bound_smaller = i
        endif
      enddo
    endif
  else
    ! If the input is not sorted, simply iterate through the input.
    if (present(mask)) then
      do i=1,size(input)
        if (mask(i)) then
          if (equality(input(i),comparison)) then
            output = i
            return
          endif
        endif
      enddo
    else
      do i=1,size(input)
        if (equality(input(i),comparison)) then
          output = i
          return
        endif
      enddo
    endif
  endif
  
  ! Check if the first larger_or_equal value is equivalent to comparison.
  if (upper_bound_larger_or_equal<=size(input)) then
    if (equality(input(upper_bound_larger_or_equal),comparison)) then
      output = upper_bound_larger_or_equal
      return
    endif
  endif
  
  ! No equivalent value has been found. Return the default value if set,
  !    or throw an error if not.
  if (present(default)) then
    output = default
  else
    call print_line(ERROR//': value not found.')
    call err()
  endif
end function

module procedure last_equivalent_integers
  if (present(mask)) then
    output = first_equivalent( input      = input(size(input):1:-1), &
                             & comparison = comparison,              &
                             & mask       = mask(size(mask):1:-1),   &
                             & default    = 0,                       &
                             & sorted     = sorted                   )
  else
    output = first_equivalent( input      = input(size(input):1:-1), &
                             & comparison = comparison,              &
                             & default    = 0,                       &
                             & sorted     = sorted                   )
  endif
  
  if (output/=0) then
    output = size(input)+1-output
  elseif (present(default)) then
    output = default
  else
    call print_line(ERROR//': value not found.')
    call err()
  endif
end procedure

function last_equivalent_ComparisonLambda(input,comparison,equality, &
   & greater_than_or_equal,mask,default,sorted) result(output) 
  class(*), intent(in)                  :: input(:)
  class(*), intent(in)                  :: comparison
  procedure(ComparisonLambda)           :: equality
  procedure(ComparisonLambda), optional :: greater_than_or_equal
  logical,  intent(in),        optional :: mask(:)
  integer,  intent(in),        optional :: default
  logical,  intent(in),        optional :: sorted
  integer                               :: output
  
  if (present(mask)) then
    output = first_equivalent(                            &
       & input                 = input(size(input):1:-1), &
       & comparison            = comparison,              &
       & equality              = equality,                &
       & greater_than_or_equal = greater_than_or_equal,   &
       & mask                  = mask(size(mask):1:-1),   &
       & default               = 0,                       &
       & sorted                = sorted                   )
  else
    output = first_equivalent(                            &
       & input                 = input(size(input):1:-1), &
       & comparison            = comparison,              &
       & equality              = equality,                &
       & greater_than_or_equal = greater_than_or_equal,   &
       & default               = 0,                       &
       & sorted                = sorted                   )
  endif
  
  if (output/=0) then
    output = size(input)+1-output
  elseif (present(default)) then
    output = default
  else
    call print_line(ERROR//': value not found.')
    call err()
  endif
end function

! WORKAROUND: requires explicit interface because of
!    polymorphic procedure argument (gfortran 10 bug).
subroutine operate_OperationLambda(input,lambda) 
  class(*), intent(inout)    :: input(:)
  procedure(OperationLambda) :: lambda
  
  integer :: i
  
  do i=1,size(input)
    call lambda(input(i))
  enddo
end subroutine

! WORKAROUND: requires explicit interface because of
!    polymorphic procedure argument (gfortran 10 bug).
function map_LogicalLambda(input,lambda) result(output) 
  class(*), intent(in)     :: input(:)
  procedure(LogicalLambda) :: lambda
  logical, allocatable     :: output(:)
  
  integer :: i
  
  output = [( lambda(input(i)), i=1, size(input) )]
end function

! WORKAROUND: requires explicit interface because of
!    polymorphic procedure argument (gfortran 10 bug).
function map_ComparisonLambda(input,lambda,comparison) result(output)
  class(*), intent(in)        :: input(:)
  procedure(ComparisonLambda) :: lambda
  class(*), intent(in)        :: comparison
  logical, allocatable        :: output(:)
  
  output = map(input,logical_lambda)
contains
  module function logical_lambda(input) result(output) 
    class(*), intent(in) :: input
    logical              :: output
    
    output = lambda(input,comparison)
  end function
end function

module procedure count_LogicalLambda
  output = count(map(input,lambda))
end procedure

module procedure count_ComparisonLambda
  output = count(map(input,lambda,comparison))
end procedure

module procedure filter_logicals
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
end procedure

module procedure filter_LogicalLambda
  output = filter(map(input,lambda),mask=mask)
end procedure

! WORKAROUND: requires explicit interface because of
!    polymorphic procedure argument (gfortran 10 bug).
function filter_ComparisonLambda(input,lambda,comparison,mask) &
   & result(output) 
  class(*), intent(in)           :: input(:)
  procedure(ComparisonLambda)    :: lambda
  class(*), intent(in)           :: comparison
  logical,  intent(in), optional :: mask(:)
  integer, allocatable           :: output(:)
  
  output = filter(input,logical_lambda,mask)
contains
  module function logical_lambda(input) result(output) 
    class(*), intent(in) :: input
    logical              :: output
    
    output = lambda(input,comparison)
  end function
end function

! WORKAROUND: requires explicit interface because of
!    polymorphic procedure argument (gfortran 10 bug).
function locate_ComparisonLambda(input,lambda,mask) result(output) 
  class(*), intent(in)          :: input(:)
  procedure(ComparisonLambda)   :: lambda
  logical, intent(in), optional :: mask(:)
  integer                       :: output
  
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
    output = first(mask, default=0)
    if (output==0) then
      call print_line(CODE_ERROR//': Calling locate on entirely masked list.')
      call err()
    endif
    do i=output+1,size(input)
      if (mask(i)) then
        if ( lambda(input(i),input(output)) .and. &
           & .not. lambda(input(output),input(i)) ) then
          output = i
        endif
      endif
    enddo
  else
    output = 1
    do i=output+1,size(input)
      if ( lambda(input(i),input(output)) .and. &
         & .not. lambda(input(output),input(i)) ) then
        output = i
      endif
    enddo
  endif
end function

module procedure is_sorted_integers
  if (size(input)<2) then
    output = .true.
  else
    output = all(input(:size(input)-1) <= input(2:))
  endif
end procedure

module procedure is_sorted_reals
  if (size(input)<2) then
    output = .true.
  else
    output = all(input(:size(input)-1) <= input(2:))
  endif
end procedure

! WORKAROUND: requires explicit interface because of
!    polymorphic procedure argument (gfortran 10 bug).
function is_sorted_ComparisonLambda(input,lambda) result(output) 
  class(*), intent(in)        :: input(:)
  procedure(ComparisonLambda) :: lambda
  logical                     :: output
  
  integer :: i
  
  output = .true.
  
  do i=1,size(input)-1
    if (       lambda(input(i+1), input(i)  ) .and. &
       & .not. lambda(input(i),   input(i+1)) ) then
      output = .false.
      return
    endif
  enddo
end function

module procedure sort_integers
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
end procedure

module procedure sort_reals
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
end procedure

module procedure sort_ComparisonLambda
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
end procedure

module procedure set_integers
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
end procedure

! WORKAROUND: requires explicit interface because of
!    polymorphic procedure argument (gfortran 10 bug).
function set_ComparisonLambda(input,lambda,mask) result(output) 
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

module procedure set_default
  if (present(optional_argument)) then
    output = optional_argument
  else
    output = default_value
  endif
end procedure

module procedure element_in_list
  output = any(lhs==rhs)
end procedure

module procedure elements_in_list
  integer :: i
  
  output = [(any(lhs(i)==rhs), i=1, size(lhs))]
end procedure
end submodule
