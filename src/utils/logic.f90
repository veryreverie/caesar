! ======================================================================
! A number of algorithms to do with logic and sorting.
! ======================================================================
module logic_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
contains

! ----------------------------------------------------------------------
! Returns the id of the first or last .true. in a logical list.
! In both cases, returns 0 if all values are .false..
! If mask is present, only returns values where mask is .true..
! ----------------------------------------------------------------------
function first(input,mask) result(output)
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
    do i=1,size(input)
      if (input(i).and.mask(i)) then
        output = i
        return
      endif
    enddo
  else
    do i=1,size(input)
      if (input(i)) then
        output = i
        return
      endif
    enddo
  endif
  
  output = 0
end function

function last(input,mask) result(output)
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
    do i=size(input),1,-1
      if (input(i).and.mask(i)) then
        output = i
        return
      endif
    enddo
  else
    do i=size(input),1,-1
      if (input(i)) then
        output = i
        return
      endif
    enddo
  endif
  
  output = 0
end function

! ----------------------------------------------------------------------
! Sorts a list of integers, and returns the ids of the sorted elements.
! e.g. sort([1,3,2,1]) returns [1,4,3,2].
! ----------------------------------------------------------------------
function sort(input) result(output)
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
end module
