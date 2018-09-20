! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use common_module
  implicit none
  
  private
  
  public :: test
  
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function test() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'test'
  output%description = 'Runs temporary code for testing purposes.'
  output%keywords = [KeywordData::]
  output%main_subroutine => test_subroutine
  output%suppress_from_helptext = .true.
end function

! ----------------------------------------------------------------------
! Main function.
! ----------------------------------------------------------------------
subroutine test_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  integer :: no_elements
  
  integer, allocatable :: mins_in(:)
  integer, allocatable :: maxs_in(:)
  integer, allocatable :: mins_out(:)
  integer, allocatable :: maxs_out(:)
  
  integer :: split
  
  integer, allocatable :: sums(:)
  
  integer :: nops
  
  integer :: i,j,k
  
  do k=1,64
    no_elements = k
    
    nops = 0
    
    mins_in = [1]
    maxs_in = [no_elements]
    sums    = [(0,i=1,no_elements)]
    do while (any(mins_in/=maxs_in)) 
      mins_out = [integer::]
      maxs_out = [integer::]
      
      do i=1,size(mins_in)
        if (mins_in(i)==maxs_in(i)) then
          mins_out = [mins_out, mins_in(i)]
          maxs_out = [maxs_out, maxs_in(i)]
        elseif (mins_in(i)/=maxs_in(i)) then
          split = mins_in(i) + (maxs_in(i)-mins_in(i)+1)/2
          
          sums(split) = sums(mins_in(i))
          do j=split,maxs_in(i)
            sums(mins_in(i)) = sums(mins_in(i)) + 10**(j-1)
            nops = nops + 1
          enddo
          do j=mins_in(i),split-1
            sums(split) = sums(split) + 10**(j-1)
            nops = nops + 1
          enddo
          
          mins_out = [mins_out, mins_in(i), split     ]
          maxs_out = [maxs_out, split-1   , maxs_in(i)]
        endif
      enddo
      
      mins_in = mins_out
      maxs_in = maxs_out
    enddo
    
    call print_line('')
    call print_line(no_elements)
    call print_line(sums)
    call print_line('N ops: '//nops)
  enddo
end subroutine
end module
