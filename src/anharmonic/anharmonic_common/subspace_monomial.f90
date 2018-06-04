! ======================================================================
! A "subspace" is a set of degenerate modes.
! A subspace monomial is a product of subspaces.
! ======================================================================
module subspace_monomial_module
  use common_module
  
  use degeneracy_module
  use subspace_coupling_module
  implicit none
  
  private
  
  public :: SubspaceMonomial
  public :: size
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
    module procedure new_SubspaceMonomial_DegenerateSubspaces
    module procedure new_SubspaceMonomial_String
  end interface
  
  interface operator(//)
    module procedure concatenate_SubspaceMonomial_DegenerateSubspace
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
!    - Constructors.
!    - Concatenation with DegenerateSubspace.
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

function new_SubspaceMonomial_DegenerateSubspaces(subspaces) result(this)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  type(SubspaceMonomial)               :: this
  
  this%ids = subspaces%id
end function

function concatenate_SubspaceMonomial_DegenerateSubspace(this,subspace) &
   & result(output)
  implicit none
  
  type(SubspaceMonomial),   intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  type(SubspaceMonomial)               :: output
  
  output = SubspaceMonomial([this%ids, subspace%id])
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
  
  class(SubspaceMonomial),  intent(in)  :: this
  type(DegenerateSubspace), intent(in)  :: subspaces(:)
  type(DegenerateSubspace), allocatable :: output(:)
  
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
function generate_subspace_monomials(subspace_coupling,subspaces, &
   & potential_expansion_order) result(output)
  implicit none
  
  type(SubspaceCoupling),   intent(in) :: subspace_coupling
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  integer,                  intent(in) :: potential_expansion_order
  type(SubspaceMonomial), allocatable  :: output(:)
  
  type(DegenerateSubspace), allocatable :: coupled_subspaces(:)
  
  ! Check inputs.
  if (size(subspace_coupling)==0) then
    call print_line(CODE_ERROR//': Empty subspace coupling.')
    call err()
  elseif (potential_expansion_order<min(2,size(subspace_coupling))) then
    call print_line(ERROR//': potential_expansion_order must be at least 2, &
       &and at least as large as maximum_coupling_order.')
    stop
  endif
  
  ! Retrieve coupled subspaces from subspace coupling.
  coupled_subspaces = subspace_coupling%coupled_subspaces(subspaces)
  
  ! Call recursive helper function to generate monomials.
  output = generate_subspace_monomials_helper( coupled_subspaces, &
                                             & potential_expansion_order)
end function

! Helper function for generate_subspace_monomials.
! This function appends a number of copies of the first subspace in
!    coupled_subspaces, then calls itself to append the next subspace etc.
! The optional argument monomial_in is a recursive argument which should only
!    be provided by recursive calls.
recursive function generate_subspace_monomials_helper(coupled_subspaces, &
   & potential_expansion_order,monomial_in) result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in)           :: coupled_subspaces(:)
  integer,                  intent(in)           :: potential_expansion_order
  type(SubspaceMonomial),   intent(in), optional :: monomial_in
  type(SubspaceMonomial), allocatable            :: output(:)
  
  type(DegenerateSubspace)              :: first_subspace
  type(DegenerateSubspace), allocatable :: remaining_subspaces(:)
  
  type(SubspaceMonomial) :: monomial
  
  ! If there are no more subspaces to append, return the monomial,
  !    but only if its size is at least 2.
  ! Size 1 monomials are ignored because they correspond to linear terms in the
  !    potential, which are zero because the structure is geometry optimised.
  if (size(coupled_subspaces)==0) then
    if (.not. present(monomial_in)) then
      call print_line(CODE_ERROR//': Empty subspace coupling.')
      call err()
    elseif (size(monomial_in)<2) then
      output = [SubspaceMonomial::]
    elseif (size(monomial_in)>=2) then
      output = [monomial_in]
    endif
    
    return
  endif
  
  ! Split the coupled subspaces into the first subspace, which will be handled
  !    by this call of the function, and all the rest, which will be handled
  !    by a recursive call.
  first_subspace      = coupled_subspaces(1)
  remaining_subspaces = coupled_subspaces(2:)
  
  ! Append copies of the first subspace as many times as still leaves space
  !    for at least one copy of every remaining subspace, and then call this
  !    function to handle the next subspace along.
  if (present(monomial_in)) then
    monomial = monomial_in
  else
    monomial = SubspaceMonomial()
  endif
  
  output = [SubspaceMonomial::]
  do while(size(monomial)+size(remaining_subspaces)<=potential_expansion_order)
    monomial = monomial // first_subspace
    output = [ output,                                                        &
           &   generate_subspace_monomials_helper( remaining_subspaces,       &
           &                                       potential_expansion_order, &
           &                                       monomial)                  &
           & ]
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

impure elemental function new_SubspaceMonomial_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(SubspaceMonomial)   :: this
  
  this = input
end function
end module
