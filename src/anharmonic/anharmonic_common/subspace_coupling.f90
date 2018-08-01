! ======================================================================
! A "subspace" is a set of degenerate modes.
! "Coupled subspaces" are a set of subspaces which are coupled to one another.
! Coupling here means that there are terms in the Hamiltonian which contain
!    products of one subspace with another.
! ======================================================================
module subspace_coupling_module
  use common_module
  implicit none
  
  private
  
  public :: SubspaceCoupling
  public :: generate_coupled_subspaces
  public :: size
  
  type, extends(Stringable) :: SubspaceCoupling
    ! The ids of the degenerate subspaces which are coupled together.
    ! ids is assumed to always be sorted in ascending order.
    integer, allocatable :: ids(:)
  contains
    ! Returns the subspaces in this coupling.
    procedure, public :: coupled_subspaces
    ! I/O.
    procedure, public :: read  => read_SubspaceCoupling
    procedure, public :: write => write_SubspaceCoupling
  end type
  
  interface SubspaceCoupling
    module procedure new_SubspaceCoupling
    module procedure new_SubspaceCoupling_String
  end interface
  
  interface operator(//)
    module procedure concatenate_SubspaceCoupling_DegenerateSubspace
  end interface
  
  interface size
    module procedure size_SubspaceCoupling
  end interface
contains

! ----------------------------------------------------------------------
! Basic functionality:
!    - Constructor.
!    - Concatenation with a DegenerateSubspace subspace.
!    - size() function.
! ----------------------------------------------------------------------
function new_SubspaceCoupling(ids) result(output)
  implicit none
  
  integer, intent(in), optional :: ids(:)
  type(SubspaceCoupling)        :: output
  
  if (present(ids)) then
    output%ids = ids
  else
    output%ids = [integer::]
  endif
end function

function concatenate_SubspaceCoupling_DegenerateSubspace(input,subspace) &
   & result(output)
  implicit none
  
  type(SubspaceCoupling),   intent(in) :: input
  type(DegenerateSubspace), intent(in) :: subspace
  type(SubspaceCoupling)               :: output
  
  output%ids = [input%ids,subspace%id]
end function

function size_SubspaceCoupling(this) result(output)
  implicit none
  
  type(SubspaceCoupling), intent(in) :: this
  integer                            :: output
  
  output = size(this%ids)
end function

! ----------------------------------------------------------------------
! Generates all sets of coupled subspaces up to a given order.
! ----------------------------------------------------------------------
! Calls its helper function, which recursively calls itself.
function generate_coupled_subspaces(subspaces,maximum_coupling_order) &
   & result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  integer,                  intent(in) :: maximum_coupling_order
  type(SubspaceCoupling), allocatable  :: output(:)
  
  integer :: coupling_order
  
  ! Check input.
  if (maximum_coupling_order<1) then
    call print_line(ERROR//': maximum_coupling_order must be at least 1.')
    stop
  endif
  
  ! Call the helper function once for each coupling order.
  ! Each call returns an array of results, and these arrays are concatenated
  !    together.
  output = [SubspaceCoupling::]
  do coupling_order=1,maximum_coupling_order
    output = [ output,                                                &
           &   generate_coupled_subspaces_helper( SubspaceCoupling(), &
           &                                      subspaces,          &
           &                                      coupling_order)     &
           & ]
  enddo
end function

! Recursive helper function for generate_coupled_subspaces.
recursive function generate_coupled_subspaces_helper(coupling_in, &
   & subspaces,coupling_order) result(output)
  implicit none
  
  type(SubspaceCoupling),   intent(in) :: coupling_in
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  integer,                  intent(in) :: coupling_order
  type(SubspaceCoupling), allocatable  :: output(:)
  
  type(SubspaceCoupling) :: coupling
  
  integer :: i
  
  output = [SubspaceCoupling::]
  do i=1,size(subspaces)
    coupling = coupling_in//subspaces(i)
    if (size(coupling)==coupling_order) then
      output = [output, coupling]
    else
      output = [ output,                                             &
             &   generate_coupled_subspaces_helper( coupling,        &
             &                                      subspaces(i+1:), &
             &                                      coupling_order)  &
             & ]
    endif
  enddo
end function

! ----------------------------------------------------------------------
! Returns the coupled subspaces in this subspace coupling.
! ----------------------------------------------------------------------
function coupled_subspaces(this,subspaces) result(output)
  implicit none
  
  class(SubspaceCoupling),  intent(in)  :: this
  type(DegenerateSubspace), intent(in)  :: subspaces(:)
  type(DegenerateSubspace), allocatable :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    output(i) = subspaces(first(subspaces%id==this%ids(i)))
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_SubspaceCoupling(this,input)
  implicit none
  
  class(SubspaceCoupling), intent(out) :: this
  type(String),            intent(in)  :: input
  
  select type(this); type is(SubspaceCoupling)
    this = SubspaceCoupling(int(split_line(input)))
  end select
end subroutine

function write_SubspaceCoupling(this) result(output)
  implicit none
  
  class(SubspaceCoupling), intent(in) :: this
  type(String)                        :: output
  
  select type(this); type is(SubspaceCoupling)
    output = join(this%ids)
  end select
end function

impure elemental function new_SubspaceCoupling_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(SubspaceCoupling)   :: this
  
  call this%read(input)
end function
end module
