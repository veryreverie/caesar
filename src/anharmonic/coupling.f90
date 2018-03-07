! ======================================================================
! Parsing and storage of which degenerate subspaces are coupled with which.
! ======================================================================
module coupling_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use stringable_module
  implicit none
  
  private
  
  public :: CoupledSubspaces
  public :: size
  public :: write_coupling_file
  public :: read_coupling_file
  public :: generate_subspace_coupling
  
  ! A list of ids of degenerate subspaces which are coupled.
  type, extends(Stringable) :: CoupledSubspaces
    integer, allocatable :: ids(:)
  contains
    ! Returns the coupled subspaces.
    procedure, public :: coupled_subspaces
    
    ! Check if this is subsidiary to the given coupling, e.g. the couplings
    !    [],[1],[2],[3],[1,2],[1,3] and [2,3] are subsidiaries of [1,2,3].
    procedure, public :: is_subsidiary_of
    
    ! I/O.
    procedure, public :: str => str_CoupledSubspaces
  end type
  
  interface size
    module procedure size_CoupledSubspaces
  end interface
contains

function size_CoupledSubspaces(this) result(output)
  implicit none
  
  type(CoupledSubspaces), intent(in) :: this
  integer                            :: output
  
  output = size(this%ids)
end function

function coupled_subspaces(this,subspaces) result(output)
  use logic_module
  use degeneracy_module
  implicit none
  
  class(CoupledSubspaces), intent(in) :: this
  type(DegenerateModes),   intent(in) :: subspaces(:)
  type(DegenerateModes), allocatable  :: output(:)
  
  integer :: i,j,ialloc
  
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    j = first(subspaces%id==this%ids(i))
    if (j==0) then
      call print_line(CODE_ERROR//': Unable to locate subspace.')
      call err()
    endif
    output(i) = subspaces(j)
  enddo
end function

function is_subsidiary_of(this,that) result(output)
  implicit none
  
  Class(CoupledSubspaces), intent(in) :: this
  type(CoupledSubspaces),  intent(in) :: that
  logical                             :: output
  
  integer :: i,j
  
  output = .true.
  
  ! Loop over the ids in this, checking if each is in that.
  do_i : do i=1,size(this)
    do j=1,size(that)
      if (this%ids(i)==that%ids(j)) then
        cycle do_i
      endif
    enddo
    
    ! Will only be reached if a mode in this is not present in that.
    output = .false.
  enddo do_i
end function

! ----------------------------------------------------------------------
! Generates all couplings between a given set of degeneracy ids at a given
!    potential expansion order and coupling order.
! ----------------------------------------------------------------------
function generate_subspace_coupling(degenerate_modes, &
   & potential_expansion_order,coupling_order) result(output)
  use degeneracy_module
  implicit none
  
  type(DegenerateModes), intent(in)   :: degenerate_modes(:)
  integer,               intent(in)   :: potential_expansion_order
  integer,               intent(in)   :: coupling_order
  type(CoupledSubspaces), allocatable :: output(:)
  
  integer :: i
  
  if (potential_expansion_order<2) then
    call print_line(ERROR//': potential_expansion_order must be at least 2.')
    stop
  elseif (coupling_order<1) then
    call print_line(ERROR//': coupling_order must be at least 1.')
    stop
  elseif (coupling_order>potential_expansion_order) then
    call print_line(ERROR//': coupling_order must be less than or equal to &
       &potential_expansion_order.')
    stop
  endif
  
  output = [CoupledSubspaces::]
  
  do i=2,potential_expansion_order
    output = [output, coupling_generator(degenerate_modes, i, coupling_order)]
  enddo
end function

! Helper function for generate_subspace coupling.
! Recursively generates all the coupling at a specific
!    potential_expansion_order.
!
! coupling_in is a recursive parameter. Each time coupling_generator is called,
!    it adds one id to the given coupling, and calls itself again.
! e.g. coupling_generator(coupling_in=[1]) will call
!     coupling_generator(coupling_in=[1,i]) with i in degeneracy_ids.
!
! Only returns couplings in ascending order, and with at most coupling_order
!    distinct elements.
recursive function coupling_generator(degenerate_modes, &
   & potential_expansion_order,coupling_order,coupling_in) result(output)
  use logic_module
  use degeneracy_module
  implicit none
  
  type(DegenerateModes),  intent(in)           :: degenerate_modes(:)
  integer,                intent(in)           :: potential_expansion_order
  integer,                intent(in)           :: coupling_order
  type(CoupledSubspaces), intent(in), optional :: coupling_in
  type(CoupledSubspaces), allocatable          :: output(:)
  
  type(CoupledSubspaces) :: coupling
  
  integer :: i,ialloc
  
  if (.not. present(coupling_in)) then
    output = [CoupledSubspaces::]
    do i=1,size(degenerate_modes)
      coupling = CoupledSubspaces([degenerate_modes(i)%id])
      output = [output, coupling_generator( degenerate_modes(i:),      &
                                          & potential_expansion_order, &
                                          & coupling_order,            &
                                          & coupling )                   ]
    enddo
  elseif (size(set(coupling_in%ids))==coupling_order) then
    allocate( coupling%ids(potential_expansion_order), &
            & stat=ialloc); call err(ialloc)
    coupling%ids(:size(coupling_in)) = coupling_in%ids
    do i=size(coupling_in)+1,potential_expansion_order
      coupling%ids(i) = coupling%ids(i-1)
    enddo
    output = [coupling]
  elseif (size(coupling_in)==potential_expansion_order) then
    output = [coupling_in]
  else
    output = [CoupledSubspaces::]
    do i=1,size(degenerate_modes)
      coupling = CoupledSubspaces([coupling_in%ids, degenerate_modes(i)%id])
      output = [output, coupling_generator( degenerate_modes(i:),      &
                                          & potential_expansion_order, &
                                          & coupling_order,            &
                                          & coupling )                   ]
    enddo
  endif
end function

! ----------------------------------------------------------------------
! File I/O.
! ----------------------------------------------------------------------
subroutine write_coupling_file(this, filename)
  use ofile_module
  implicit none
  
  type(CoupledSubspaces), intent(in) :: this(:)
  type(String),           intent(in) :: filename
  
  type(OFile) :: coupling_file
  
  integer :: i
  
  coupling_file = filename
  call coupling_file%print_line('! Couplings between degenerate subspaces.')
  do i=1,size(this)
    call coupling_file%print_line(this(i))
  enddo
end subroutine

function read_coupling_file(filename) result(this)
  use ifile_module
  implicit none
  
  type(String), intent(in)            :: filename
  type(CoupledSubspaces), allocatable :: this(:)
  
  type(IFile) :: coupling_file
  integer     :: i,ialloc
  
  coupling_file = filename
  allocate(this(size(coupling_file)-1), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    this(i)%ids = int(split(coupling_file%line(i+1)))
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
function OLDFUNCTION_calculate_all_coupling(input, modes) result(output)
  use integer_arrays_module
  use normal_mode_module
  implicit none
  
  type(CoupledSubspaces), intent(in)  :: input(:)
  type(ComplexMode),      intent(in)  :: modes(:)
  type(CoupledSubspaces), allocatable :: output(:)
  
  integer :: no_modes
  
  integer :: max_no_coupled
  integer :: no_couplings
  
  integer,          allocatable :: sizes(:)
  type(IntArray2D), allocatable :: ids(:)
  
  type(CoupledSubspaces), allocatable :: couplings(:)
  integer,            allocatable :: couplings_sizes(:)
  
  logical, allocatable :: mode_unaccounted_for(:)
  logical, allocatable :: duplicate(:)
  
  integer :: i,j,k,l,ialloc
  integer :: s
  
  no_modes = size(modes)
  
  ! ------------------------------
  ! Check that no couplings include translational modes.
  ! ------------------------------
  do i=1,size(input)
    do j=1,size(input(i))
      if (modes(input(i)%ids(j))%translational_mode) then
        call print_line('Error: the translational mode '//input(i)%ids(j)// &
           & 'has been included in coupling '//i//' at the gamma-point.')
        stop
      endif
    enddo
  enddo
  
  ! ------------------------------
  ! Check that all couplings are in ascending order and within [1,no_modes].
  ! ------------------------------
  do i=1,size(input)
    do j=1,size(input(i))
      if (input(i)%ids(j)<1 .or. input(i)%ids(j)>no_modes) then
        call print_line('Error: mode '//j//' of coupling '//i//', '// &
           & input(i)%ids//' is outside of the expected range.')
        stop
      endif
      if (j>1) then
        if (input(i)%ids(j)<=input(i)%ids(j-1)) then
          call print_line('Error: coupling '//i//', '//input(i)%ids// &
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
      allocate( couplings(j)%ids(size(ids(s)%i(j))), &
              & stat=ialloc); call err(ialloc)
      do k=1,size(couplings(j))
        couplings(j)%ids(k) = input(i)%ids( ids(s)%i(j)%i(k) )
      enddo
    enddo
    l = l + sizes(s)
  enddo
  
  ! ------------------------------
  ! Remove duplicates, and add in uncoupled modes.
  ! ------------------------------
  
  ! Identify missing modes.
  allocate(mode_unaccounted_for(no_modes), stat=ialloc); call err(ialloc)
  mode_unaccounted_for = .true.
  do i=1,size(couplings)
    do j=1,size(couplings(i))
      mode_unaccounted_for(couplings(i)%ids(j)) = .false.
    enddo
  enddo
  
  ! Mark translational modes as not missing.
  do i=1,no_modes
    if (modes(i)%translational_mode) then
      mode_unaccounted_for(i) = .false.
    endif
  enddo
  
  ! Identify duplicate modes.
  allocate(duplicate(size(couplings)), stat=ialloc); call err(ialloc)
  duplicate = .false.
  do i=1,size(couplings)
    do j=1,i-1
      if (size(couplings(i))==size(couplings(j))) then
        if (all(couplings(i)%ids==couplings(j)%ids)) then
          duplicate(i) = .true.
        endif
      endif
    enddo
  enddo
  
  ! Construct output.
  ! Couplings are sorted by size order.
  allocate( output( size(couplings)             &
          &       + count(mode_unaccounted_for) &
          &       - count(duplicate)            &
          &       + 1),                         &
          & stat=ialloc); call err(ialloc)
  
  ! Add the blank coupling.
  output(1)%ids=[integer::]
  
  ! Add in single modes which have not been specified as part of couplings.
  j = 1
  do i=1,size(mode_unaccounted_for)
    if (mode_unaccounted_for(i)) then
      j = j + 1
      output(j)%ids = [i]
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
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
recursive function str_CoupledSubspaces(this) result(output)
  implicit none
  
  class(CoupledSubspaces), intent(in) :: this
  type(String)                        :: output
  
  output = join(this%ids)
end function
end module
