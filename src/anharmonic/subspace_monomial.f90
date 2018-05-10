! ======================================================================
! A "subspace" is a set of degenerate modes.
! A subspace monomial is a product of subspaces.
! ======================================================================
module subspace_monomial_module
  use common_module
  
  use degeneracy_module
  use coupled_subspaces_module
  implicit none
  
  private
  
  public :: SubspaceMonomial
  public :: size
  public :: write_coupling_file
  public :: read_coupling_file
  public :: generate_subspace_coupling
  public :: generate_subspace_monomials
  public :: operator(//)
  public :: operator(==)
  public :: operator(/=)
  
  type, extends(Stringable) :: SubspaceMonomial
    ! The ids of the degenerate subspaces in the monomial.
    ! e.g. if ids=[1,2,2,5] then the monomial is (s1).(s2)^2.(s5), where s1 is
    !    the subspace with id 1.
    integer, allocatable :: ids(:)
  contains
    ! Returns the coupled subspaces.
    procedure, public :: coupled_subspaces
    
    ! Check if this is subsidiary to the given coupling, e.g. the couplings
    !    [],[1],[2],[3],[1,2],[1,3] and [2,3] are subsidiaries of [1,2,3].
    procedure, public :: is_subsidiary_of
    
    ! I/O.
    procedure, public :: read  => read_SubspaceMonomial
    procedure, public :: write => write_SubspaceMonomial
  end type
  
  interface SubspaceMonomial
    module procedure new_SubspaceMonomial
  end interface
  
  interface operator(//)
    module procedure concatenate_SubspaceMonomial_integer
  end interface
  
  interface size
    module procedure size_SubspaceMonomial
  end interface
  
  interface operator(==)
    module procedure equality_SubspaceMonomial_SubspaceMonomial
  end interface
  
  interface operator(/=)
    module procedure non_equality_SubspaceMonomial_SubspaceMonomial
  end interface
contains

! ----------------------------------------------------------------------
! Basic functionality.
!    - Constructor.
!    - Concatenation with subspace id.
!    - size() function.
!    - equality and non-equality with other SubspaceMonomials.
! ----------------------------------------------------------------------
function new_SubspaceMonomial(ids) result(this)
  implicit none
  
  integer, intent(in), optional :: ids(:)
  type(SubspaceMonomial)        :: this
  
  if (present(ids)) then
    this%ids = ids
  else
    this%ids = [integer::]
  endif
end function

function concatenate_SubspaceMonomial_integer(this,id) result(output)
  implicit none
  
  type(SubspaceMonomial), intent(in) :: this
  integer,                intent(in) :: id
  type(SubspaceMonomial)             :: output
  
  output = SubspaceMonomial([this%ids, id])
end function

function size_SubspaceMonomial(this) result(output)
  implicit none
  
  type(SubspaceMonomial), intent(in) :: this
  integer                            :: output
  
  output = size(this%ids)
end function

impure elemental function equality_SubspaceMonomial_SubspaceMonomial(this, &
   & that) result(output)
  implicit none
  
  type(SubspaceMonomial), intent(in) :: this
  type(SubspaceMonomial), intent(in) :: that
  logical                            :: output
  
  output = all(this%ids==that%ids)
end function

impure elemental function non_equality_SubspaceMonomial_SubspaceMonomial( &
   & this,that) result(output)
  implicit none
  
  type(SubspaceMonomial), intent(in) :: this
  type(SubspaceMonomial), intent(in) :: that
  logical                            :: output
  
  output = .not. this==that
end function

! ----------------------------------------------------------------------
! Use stored IDs to return the actual objects the IDs represent.
! ----------------------------------------------------------------------
function coupled_subspaces(this,subspaces) result(output)
  implicit none
  
  class(SubspaceMonomial), intent(in) :: this
  type(DegenerateModes),   intent(in) :: subspaces(:)
  type(DegenerateModes), allocatable  :: output(:)
  
  integer :: i,j,ialloc
  
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    j = first(subspaces%id==this%ids(i))
    output(i) = subspaces(j)
  enddo
end function

! ----------------------------------------------------------------------
! Check if a coupling is a subsidiary of another coupling.
! ----------------------------------------------------------------------
function is_subsidiary_of(this,that) result(output)
  implicit none
  
  Class(SubspaceMonomial), intent(in) :: this
  type(SubspaceMonomial),  intent(in) :: that
  logical                             :: output
  
  integer :: i,j
  
  ! Loop over the ids in this, checking if each is in that.
  do i=1,size(this)
    j = first(that%ids==this%ids(i), default=0)
    if (j==0) then
      output = .false.
      return
    endif
  enddo
  
  output = .true.
end function

! ----------------------------------------------------------------------
! Generates all coupling monomials corresponding to a given subspace coupling,
!    up to the given potential expansion order.
! ----------------------------------------------------------------------
! e.g. given the subspace coupling with ids [3,5] and potential order 3,
!    will return coupling monomials [3,3,5] and [3,5,5].
function generate_subspace_monomials(coupled_subspaces, &
   & potential_expansion_order) result(output)
  implicit none
  
  type(CoupledSubspaces), intent(in)  :: coupled_subspaces
  integer,                intent(in)  :: potential_expansion_order
  type(SubspaceMonomial), allocatable :: output(:)
  
  if (potential_expansion_order<min(2,size(coupled_subspaces))) then
    call print_line(ERROR//': potential_expansion_order must be at least 2, &
       &and at least as large as maximum_coupling_order.')
    stop
  endif
  
  output = generate_subspace_monomials_helper( SubspaceMonomial(), &
                                             & coupled_subspaces,  &
                                             & potential_expansion_order)
end function

! Helper function for generate_subspace_monomials.
recursive function generate_subspace_monomials_helper(monomial_in, &
   & coupled_subspaces,potential_expansion_order) result(output)
  implicit none
  
  type(SubspaceMonomial), intent(in)  :: monomial_in
  type(CoupledSubspaces), intent(in)  :: coupled_subspaces
  integer,                intent(in)  :: potential_expansion_order
  type(SubspaceMonomial), allocatable :: output(:)
  
  integer                :: first_subspace
  type(CoupledSubspaces) :: remaining_subspaces
  type(SubspaceMonomial) :: monomial
  
  output = [SubspaceMonomial::]
  
  ! This function appends a number of copies of the first subspace in
  !    coupled_subspaces, then calls itself to append the next subspace etc.
  
  ! If there are no more subspaces to append, return the monomial.
  if (size(coupled_subspaces)==0) then
    if (size(monomial_in)>=2) then
      output = [monomial_in]
    endif
    return
  endif
  
  ! Generate the CoupledSubspaces containing the subspaces which are not
  !    handled by this call of the function, that is to say all but the first.
  ! Split the coupled subspaces into the first subspace, which will be handled
  !    by this call of the function, and all the rest, which will be handled
  !    by a recursive call.
  first_subspace      = coupled_subspaces%ids(1)
  remaining_subspaces = CoupledSubspaces(coupled_subspaces%ids(2:))
  
  ! Append copies of the first subspace as many times as still leaves space
  !    for at least one copy of every remaining subspace, and then call this
  !    function to handle the next subspace along.
  monomial = monomial_in
  do while(size(monomial)+size(coupled_subspaces)<=potential_expansion_order)
    monomial = monomial // first_subspace
    output = [ output,                             &
           &   generate_subspace_monomials_helper( &
           &            monomial,                  &
           &            remaining_subspaces,       &
           &            potential_expansion_order) &
           & ]
  enddo
end function

! ----------------------------------------------------------------------
! Generates all couplings between a given set of degeneracy ids at a given
!    potential expansion order and coupling order.
! ----------------------------------------------------------------------
function generate_subspace_coupling(degenerate_modes, &
   & potential_expansion_order,coupling_order) result(output)
  implicit none
  
  type(DegenerateModes), intent(in)   :: degenerate_modes(:)
  integer,               intent(in)   :: potential_expansion_order
  integer,               intent(in)   :: coupling_order
  type(SubspaceMonomial), allocatable :: output(:)
  
  type(SubspaceMonomial), allocatable :: old_couplings(:)
  type(SubspaceMonomial), allocatable :: new_couplings(:)
  
  integer :: i,ialloc
  
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
  
  output = [SubspaceMonomial::]
  
  do i=2,potential_expansion_order
    ! N.B. the contents of this loop is equivalent to the line:
    !    output = [output, coupling_generator(...)]
    ! But that line triggers an apparent Nagfort bug.
    old_couplings = output
    new_couplings = coupling_generator(degenerate_modes, i, coupling_order)
    
    deallocate(output, stat=ialloc); call err(ialloc)
    allocate( output(size(old_couplings)+size(new_couplings)), &
            & stat=ialloc); call err(ialloc)
    output(                      : size(old_couplings)) = old_couplings
    output(size(old_couplings)+1 :                    ) = new_couplings
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
  implicit none
  
  type(DegenerateModes),  intent(in)           :: degenerate_modes(:)
  integer,                intent(in)           :: potential_expansion_order
  integer,                intent(in)           :: coupling_order
  type(SubspaceMonomial), intent(in), optional :: coupling_in
  type(SubspaceMonomial), allocatable          :: output(:)
  
  type(SubspaceMonomial) :: coupling
  type(SubspaceMonomial) :: coupling_out
  
  type(SubspaceMonomial), allocatable :: old_couplings(:)
  type(SubspaceMonomial), allocatable :: new_couplings(:)
  
  integer :: i
  
  if (present(coupling_in)) then
    coupling = coupling_in
  else
    coupling = SubspaceMonomial([integer::])
  endif
  
  if (size(coupling)==potential_expansion_order) then
    ! The input coupling already contains as many terms as required. Return it.
    output = [coupling]
  elseif (size(set(coupling%ids))==coupling_order) then
    ! The input coupling already contains coupling_order distinct subspaces.
    ! Append the last added subspace until potential_expansion_order.
    ! e.g. if coupling_order=2 and potential_expansion_order=6 then
    !    coupling_in=[3,3,4] -> coupling_out=[3,3,4,4,4,4].
    coupling_out = coupling
    do i=1,potential_expansion_order-size(coupling)
      coupling_out = coupling_out//degenerate_modes(1)%id
    enddo
    output = [coupling_out]
  else
    ! Recursively call coupling_generator after appending each subspace id.
    output = [SubspaceMonomial::]
    do i=1,size(degenerate_modes)
      coupling_out = coupling//degenerate_modes(i)%id
      
      ! N.B. this is equivalent to output = [output, coupling_generator(...)]
      ! But again, Nagfort appears to have a compiler bug at that line.
      old_couplings = output
      new_couplings = coupling_generator( degenerate_modes(i:),      &
                                        & potential_expansion_order, &
                                        & coupling_order,            &
                                        & coupling_out)
      output = [old_couplings, new_couplings]
    enddo
  endif
end function

! ----------------------------------------------------------------------
! File I/O.
! ----------------------------------------------------------------------
subroutine write_coupling_file(this, filename)
  implicit none
  
  type(SubspaceMonomial), intent(in) :: this(:)
  type(String),           intent(in) :: filename
  
  type(OFile) :: coupling_file
  
  integer :: i
  
  coupling_file = OFile(filename)
  call coupling_file%print_line('! Couplings between degenerate subspaces.')
  do i=1,size(this)
    call coupling_file%print_line(this(i))
  enddo
end subroutine

function read_coupling_file(filename) result(this)
  implicit none
  
  type(String), intent(in)            :: filename
  type(SubspaceMonomial), allocatable :: this(:)
  
  type(IFile) :: coupling_file
  integer     :: i,ialloc
  
  coupling_file = IFile(filename)
  allocate(this(size(coupling_file)-1), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    this(i)%ids = int(split_line(coupling_file%line(i+1)))
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
  implicit none
  
  type(SubspaceMonomial), intent(in)  :: input(:)
  type(ComplexMode),      intent(in)  :: modes(:)
  type(SubspaceMonomial), allocatable :: output(:)
  
  integer :: no_modes
  
  integer :: max_no_coupled
  integer :: no_couplings
  
  integer,          allocatable :: sizes(:)
  type(IntArray2D), allocatable :: ids(:)
  
  type(SubspaceMonomial), allocatable :: couplings(:)
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
subroutine read_SubspaceMonomial(this,input)
  implicit none
  
  class(SubspaceMonomial), intent(out) :: this
  type(String),            intent(in)  :: input
  
  select type(this); type is(SubspaceMonomial)
    this = SubspaceMonomial(int(split_line(input)))
  end select
end subroutine

function write_SubspaceMonomial(this) result(output)
  implicit none
  
  class(SubspaceMonomial), intent(in) :: this
  type(String)                        :: output
  
  select type(this); type is(SubspaceMonomial)
    output = join(this%ids)
  end select
end function
end module
