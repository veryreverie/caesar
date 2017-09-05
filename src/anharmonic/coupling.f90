! ======================================================================
! Parsing and storage of which modes are coupled with which.
! ======================================================================
module coupling_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  ! A list of mode ids of modes which are coupled.
  type :: CoupledModes
    ! The ids of the modes which are coupled together.
    integer, allocatable :: modes(:)
    ! The ids of the couplings which are subsidiary to this one, e.g.
    !    the coupling [1,2] is subsidiary to [1,2,3]
    integer, allocatable :: subsidiaries(:)
  end type
  
  interface size
    module procedure size_CoupledModes
  end interface
contains

function size_CoupledModes(this) result(output)
  implicit none
  
  type(CoupledModes), intent(in) :: this
  integer                        :: output
  
  output = size(this%modes)
end function

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
! e.g. [[3 5 7]] becomes [[], [3], [5], [3,5], [7], [3,7], [5,7], [3,5,7]]
!
! Algorithmic information:
! The coupling with zero elements only produces itself.
! The nth coupling produces all of the couplings from the previous couplings,
!    both with and without n.
! []      -> []
! [1]     -> [], [1]
! [1,2]   -> [], [1], [2], [1,2]
! [1,2,3] -> [], [1], [2], [3], [1,2], [1,3], [2,3], [1,2,3]
! These ids are then used as indices for the modes, so e.g.
! [1,3,4] -> [], [1], [3], [4], [1,3], [1,4], [3,4], [1,3,4]
! Then duplicates are removed and missing modes added, so e.g.
! [1,3,4], [1,3] -> [], [1], [2], [3], [4], [1,3], [1,4], [3,4], [1,3,4]
function calculate_all_coupling(input, no_modes) result(output)
  use integer_arrays_module
  implicit none
  
  type(CoupledModes), intent(in)  :: input(:)
  integer,            intent(in)  :: no_modes
  type(CoupledModes), allocatable :: output(:)
  
  integer :: max_no_coupled
  integer :: no_couplings
  
  integer,          allocatable :: sizes(:)
  type(IntArray2D), allocatable :: ids(:)
  
  type(CoupledModes), allocatable :: couplings(:)
  integer,            allocatable :: couplings_sizes(:)
  
  logical, allocatable :: single_mode_present(:)
  logical, allocatable :: duplicate(:)
  
  integer :: i,j,k,l,k2,ialloc
  integer :: s
  
  ! ------------------------------
  ! Check that all couplings are in ascending order and within [1,no_modes].
  ! ------------------------------
  do i=1,size(input)
    do j=1,size(input(i))
      if (input(i)%modes(j)<1 .or. input(i)%modes(j)>no_modes) then
        call print_line('Error: mode '//j//' of coupling '//i//', '// &
           & input(i)%modes//' is outside of the expected range.')
        stop
      endif
      if (j>1) then
        if (input(i)%modes(j)<=input(i)%modes(j-1)) then
          call print_line('Error: coupling '//i//', '//input(i)%modes// &
             & ' is not in ascending order.')
          stop
        endif
      endif
    enddo
  enddo
  
  ! ------------------------------
  ! Calculate the largest single coupling (e.g. [1,4,7] is size 3).
  ! ------------------------------
  max_no_coupled = 0
  do i=1,size(input)
    max_no_coupled = max(max_no_coupled, size(input(i)))
  enddo
  
  ! ------------------------------
  ! Calculate the number of individual terms for a given set of coupled modes.
  ! e.g. [1,2] produces [1], [2] and [1,2], and is of size 3.
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
  
  ! Base case: single mode. ids(1) = [[1]]
  if (size(ids)>0) then
    ids(1) = [array([1])]
  endif
  
  ! Further cases : ids(i) = [ids(i-1), [i], ids(i-1)//i]
  do i=2,size(ids)
    ids(i) = ids(i-1) // [array([i])] // ids(i-1)
    do j=sizes(i-1)+1,size(ids(i))
      ids(i)%i(j) = ids(i)%i(j) // [i]
    enddo
  enddo
  
  ! ------------------------------
  ! Calculate the total number of couplings.
  ! ------------------------------
  no_couplings = 0
  do i=1,size(input)
    no_couplings = no_couplings + sizes(size(input(i)))
  enddo
  
  ! ------------------------------
  ! Calculate all couplings.
  ! ------------------------------
  allocate(couplings(no_couplings), stat=ialloc); call err(ialloc)
  l = 0
  do i=1,size(input)
    s = size(input(i))
    do j=l+1,l+sizes(s)
      allocate( couplings(j)%modes(size(ids(s)%i(j))), &
              & stat=ialloc); call err(ialloc)
      do k=1,size(couplings(j))
        couplings(j)%modes(k) = input(i)%modes( ids(s)%i(j)%i(k) )
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
    do j=1,size(couplings(i))
      single_mode_present(couplings(i)%modes(j)) = .true.
    enddo
  enddo
  
  ! Identify duplicate modes.
  allocate(duplicate(size(couplings)), stat=ialloc); call err(ialloc)
  duplicate = .false.
  do i=1,size(couplings)
    do j=1,i-1
      if (size(couplings(i))==size(couplings(j))) then
        if (all(couplings(i)%modes==couplings(j)%modes)) then
          duplicate(i) = .true.
        endif
      endif
    enddo
  enddo
  
  ! Construct output.
  ! Couplings are sorted by size order.
  allocate( output( size(couplings)                  &
          &       + count(.not. single_mode_present) &
          &       - count(duplicate)                 &
          &       + 1),                              &
          & stat=ialloc); call err(ialloc)
  
  ! Add the blank coupling.
  output(1)%modes=[integer::]
  
  ! Add in single modes which have not been specified as part of couplings.
  j = 1
  do i=1,size(single_mode_present)
    if (.not. single_mode_present(i)) then
      j = j + 1
      output(j)%modes = [i]
    endif
  enddo
  
  ! Add in all other couplings, in order of size.
  allocate(couplings_sizes(size(couplings)), stat=ialloc); call err(ialloc)
  
  do i=1,size(couplings)
    if (duplicate(i)) then
      couplings_sizes(i) = -1
    else
      couplings_sizes(i) = size(couplings(i))
    endif
  enddo
  
  do i=1,size(couplings)-count(duplicate)
    k = minloc(couplings_sizes, 1, mask=(couplings_sizes/=1))
    j = j + 1
    output(j) = couplings(k)
    couplings_sizes(k) = -1
  enddo
  
  ! ------------------------------
  ! Calculate subsidiary couplings, e.g. coupling [1,2] has
  !    subsidiary couplings [], [1] and [2]
  ! ------------------------------
  do i=1,size(output)
    if (size(output(i)) <= 1) then
      output(i)%subsidiaries = [integer::]
    else
      allocate( output(i)%subsidiaries(sizes(size(output(i)))), &
              & stat=ialloc); call err(ialloc)
      j = 0
      do_k : do k=1,i-1
        if (size(output(k)) < size(output(i))) then
          do_k2 : do k2=1,size(output(k))
            if (.not. any(output(k)%modes(k2)==output(i)%modes)) then
              ! The mode k2 does not appear in coupling i.
              !    -> k is not a subsidiary of i.
              cycle do_k
            endif
          enddo do_k2
          
          ! Coupling k is a subsidiary coupling to coupling i.
          j = j + 1
          output(i)%subsidiaries(j) = k
        endif
      enddo do_k
    endif
  enddo
end function
end module
