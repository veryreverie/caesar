! ======================================================================
! Parsing and storage of which modes are coupled with which.
! ======================================================================
module coupling_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  ! A list of mode ids of modes which are coupled.
  type :: CoupledModes
    integer, allocatable :: modes(:)
  end type
  
  ! A list of ids for working out combinations of modes.
  type, private :: IdData
    integer, allocatable :: ids(:)
  end type
  
  ! A list of IdData.
  type, private :: IdDatas
    type(IdData), allocatable :: ids(:)
  end type
contains

! ----------------------------------------------------------------------
! Takes a list of couplings, and appends all subsidiary couplings.
! ----------------------------------------------------------------------
! e.g. ['1 2 3'] becomes ['1 2 3', '1 2', '1 3', '2 3', '1', '2', '3']
!
! Algorithmic information:
! It is convenient to include coupling with zero elements, ''.
! The coupling with zero elements only produces itself (which is neglected)
! The nth coupling produces all of the couplings from the previous couplings,
!    both with and without n.
! ['']      -> ['']                                          -> []
! ['1']     -> ['']               + ['1']                    -> [ '1']
! ['1 2']   -> ['' '1']           + ['2' '1 2']              -> ['1' '2' '1 2']
! ['1 2 3'] -> ['' '1' '2' '1 2'] + ['3' '1 3' '2 3' '1 2 3']-> ...
! The '' coupling is ignored.
function calculate_all_coupling(input) result(output)
  implicit none
  
  type(CoupledModes), intent(in)  :: input(:)
  type(CoupledModes), allocatable :: output(:)
  
  integer :: no_couplings
  integer :: max_no_coupled
  
  integer,       allocatable :: sizes(:)
  type(IdDatas), allocatable :: ids(:)
  
  integer :: i,j,k,l,ialloc
  integer :: s
  
  ! ------------------------------
  ! Calculate the largest single coupling (e.g. '1 4 7' is size 3).
  ! ------------------------------
  max_no_coupled = 0
  do i=1,size(input)
    max_no_coupled = max(max_no_coupled, size(input(i)%modes))
  enddo
  
  ! ------------------------------
  ! Calculate the number of individual terms for a given set of coupled modes.
  ! e.g. '1 2' produces '1', '2' and '1 2', and is of size 3.
  ! ------------------------------
  allocate(sizes(max_no_coupled), stat=ialloc); call err(ialloc)
  do i=1,max_no_coupled
    sizes(i) = 2**i-1
  enddo
  
  ! ------------------------------
  ! Calculate ids.
  ! ids = [[[1]], [[1],[2],[1,2]], [[1],[2],[1,2],[3],[1,3],[2,3],[1,2,3]] ...]
  ! ------------------------------
  
  ! Allocate space for ids.
  allocate(ids(max_no_coupled), stat=ialloc); call err(ialloc)
  do i=1,size(ids)
    allocate(ids(i)%ids(sizes(i)), stat=ialloc); call err(ialloc)
  enddo
  
  ! Base case: single mode. ids(1) = [[1]]
  ids(1)%ids(1)%ids = [1]
  
  ! Further cases : ids(i) = [ids(i-1), [i], ids(i-1)//i]
  do i=2,size(ids)
    ! ids(i)
    do j=1,size(ids(i-1)%ids)
      ids(i)%ids(j)%ids = ids(i-1)%ids(j)%ids
    enddo
    ! [i]
    ids(i)%ids(sizes(i-1)+1)%ids = [i]
    ! ids(i)//i  e.g. [[1], [2]] -> [[1,i], [2,i]]
    do j=1,size(ids(i-1)%ids)
      k = sizes(i-1)+1+j
      s = size(ids(i-1)%ids(j)%ids)
      allocate(ids(i)%ids(k)%ids(s+1), stat=ialloc); call err(ialloc)
      ids(i)%ids(k)%ids(:s) = ids(i-1)%ids(j)%ids
      ids(i)%ids(k)%ids(s+1) = i
    enddo
  enddo
  
  ! ------------------------------
  ! Calculate the total number of couplings,
  ! ------------------------------
  no_couplings = 0
  do i=1,size(input)
    no_couplings = no_couplings + sizes(size(input(i)%modes))
  enddo
  
  ! ------------------------------
  ! Calculate output.
  ! ------------------------------
  allocate(output(no_couplings), stat=ialloc); call err(ialloc)
  k = 0
  do i=1,size(input)
    s = size(input(i)%modes)
    do j=l+1,l+sizes(s)
      allocate( output(j)%modes(size(ids(s)%ids(j)%ids)), &
              & stat=ialloc); call err(ialloc)
      do k=1,size(output(j)%modes)
        output(j)%modes(k) = input(i)%modes( ids(s)%ids(j)%ids(k) )
      enddo
    enddo
    l = l + sizes(s)
  enddo
end function
end module
