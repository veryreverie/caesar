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

subroutine write_coupling_file(this, filename)
  implicit none
  
  type(CoupledModes), intent(in) :: this(:)
  type(String),       intent(in) :: filename
  
  integer :: coupling_file
  integer :: i
  
  coupling_file = open_write_file(filename)
  call print_line(coupling_file,'! Couplings between modes.')
  do i=1,size(this)
    call print_line(coupling_file,this(i)%modes)
  enddo
  close(coupling_file)
end subroutine

function read_coupling_file(filename) result(this)
  implicit none
  
  type(String), intent(in)        :: filename
  type(CoupledModes), allocatable :: this(:)
  
  type(String), allocatable :: coupling_file(:)
  integer                   :: i,ialloc
  
  coupling_file = read_lines(filename)
  allocate(this(size(coupling_file)-1), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    this(i)%modes = int(split(coupling_file(i+1)))
  enddo
end function

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
! These ids are then used as indices for the modes, so e.g.
! ['1 3 4'] -> ['1' '3' '1 3' '4' '1 4' '3 4' '1 3 4']
! Then duplicates are removed and missing modes added, so e.g.
! ['1 3 4'], ['1 3'] -> ['1' '3' '1 3' '4' '1 4' '3 4' '1 3 4' '1' '3' '1 3']
!                    -> ['1' '3' '1 3' '4' '1 4' '3 4' '1 3 4' '2']
function calculate_all_coupling(input, no_modes) result(output)
  implicit none
  
  type(CoupledModes), intent(in)  :: input(:)
  integer,            intent(in)  :: no_modes
  type(CoupledModes), allocatable :: output(:)
  
  integer :: max_no_coupled
  integer :: no_couplings
  
  integer,       allocatable :: sizes(:)
  type(IdDatas), allocatable :: ids(:)
  
  type(CoupledModes), allocatable :: couplings(:)
  
  logical, allocatable :: single_mode_present(:)
  logical, allocatable :: duplicate(:)
  
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
  ! Calculate all couplings.
  ! ------------------------------
  allocate(couplings(no_couplings), stat=ialloc); call err(ialloc)
  k = 0
  do i=1,size(input)
    s = size(input(i)%modes)
    do j=l+1,l+sizes(s)
      allocate( couplings(j)%modes(size(ids(s)%ids(j)%ids)), &
              & stat=ialloc); call err(ialloc)
      do k=1,size(couplings(j)%modes)
        couplings(j)%modes(k) = input(i)%modes( ids(s)%ids(j)%ids(k) )
      enddo
    enddo
    l = l + sizes(s)
  enddo
  
  ! ------------------------------
  ! Remove duplicates, and add in uncoupled modes.
  ! ------------------------------
  
  ! Identify missing modes.
  allocate(single_mode_present(no_modes), stat=ialloc); call err(ialloc)
  single_mode_present = .false.
  do i=1,size(couplings)
    do j=1,size(couplings(i)%modes)
      single_mode_present(couplings(i)%modes(j)) = .true.
    enddo
  enddo
  
  ! Identify duplicate modes.
  allocate(duplicate(size(couplings)), stat=ialloc); call err(ialloc)
  duplicate = .false.
  do i=1,size(couplings)
    do_j : do j=1,i-1
      if (size(couplings(i)%modes)==size(couplings(j)%modes)) then
        do_k : do k=1,size(couplings(i)%modes)
          do l=1,size(couplings(j)%modes)
            if (couplings(i)%modes(k)==couplings(j)%modes(l)) then
              ! Mode k is present in both couplings i and j.
              cycle do_k
            endif
          enddo
          
          ! Mode k is not present in j. Couplings i and j are not duplicates.
          cycle do_j
        enddo do_k
        
        ! Couplings i and j are duplicates.
        duplicate(i) = .true.
      endif
    enddo do_j
  enddo
  
  ! Construct output.
  allocate(output( size(couplings) &
               & + count(.not. single_mode_present) &
               & - count(duplicate) ))
  j = 0
  do i=1,size(couplings)
    if (.not. duplicate(i)) then
      j = j + 1
      output(j) = couplings(i)
    endif
  enddo
  do i=1,no_modes
    if (.not. single_mode_present(i)) then
      j = j + 1
      output(j)%modes = [i]
    endif
  enddo
end function
end module
