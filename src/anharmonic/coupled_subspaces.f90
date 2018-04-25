! ======================================================================
! A "subspace" is a set of degenerate modes.
! "Coupled subspaces" are a set of subspaces which are coupled to one another.
! Coupling here means that there are terms in the Hamiltonian which contain
!    products of one subspace with another.
! ======================================================================
module coupled_subspaces_module
  use common_module
  
  use degeneracy_module
  implicit none
  
  private
  
  public :: CoupledSubspaces
  public :: generate_coupled_subspaces
  public :: size
  
  type, extends(Stringable) :: CoupledSubspaces
    ! The ids of the degenerate subspaces which are coupled together.
    integer, allocatable :: ids(:)
  contains
    ! I/O.
    procedure, public :: to_String => to_String_CoupledSubspaces
  end type
  
  interface CoupledSubspaces
    module procedure new_CoupledSubspaces
  end interface
  
  interface operator(//)
    module procedure concatenate_CoupledSubspaces_DegenerateModes
  end interface
  
  interface size
    module procedure size_CoupledSubspaces
  end interface
contains

! ----------------------------------------------------------------------
! Basic functionality:
!    - Constructor.
!    - Concatenation with a DegenerateModes subspace.
!    - size() function.
! ----------------------------------------------------------------------
function new_CoupledSubspaces() result(output)
  implicit none
  
  type(CoupledSubspaces) :: output
  
  output%ids = [integer::]
end function

function concatenate_CoupledSubspaces_DegenerateModes(input,subspace) &
   & result(output)
  implicit none
  
  type(CoupledSubspaces), intent(in) :: input
  type(DegenerateModes),  intent(in) :: subspace
  type(CoupledSubspaces)             :: output
  
  output%ids = [input%ids,subspace%id]
end function

function size_CoupledSubspaces(this) result(output)
  implicit none
  
  type(CoupledSubspaces), intent(in) :: this
  integer                            :: output
  
  output = size(this%ids)
end function

! ----------------------------------------------------------------------
! Generates all sets of coupled subspaces up to a given order.
! ----------------------------------------------------------------------
! Calls its helper function, which recursively calls itself.
function generate_coupled_subspaces(degenerate_modes,maximum_coupling_order) &
   & result(output)
  implicit none
  
  type(DegenerateModes), intent(in)   :: degenerate_modes(:)
  integer,               intent(in)   :: maximum_coupling_order
  type(CoupledSubspaces), allocatable :: output(:)
  
  integer :: coupling_order
  
  ! Check input.
  if (maximum_coupling_order<1) then
    call print_line(ERROR//': maximum_coupling_order must be at least 1.')
    stop
  endif
  
  ! Call the helper function once for each coupling order.
  ! Each call returns an array of results, and these arrays are concatenated
  !    together.
  output = [CoupledSubspaces::]
  do coupling_order=1,maximum_coupling_order
    output = [ output,                            &
           &   generate_coupled_subspaces_helper( &
           &                  CoupledSubspaces(), &
           &                  degenerate_modes,   &
           &                  coupling_order)     &
           & ]
  enddo
end function

! Recursive helper function for generate_coupled_subspaces.
recursive function generate_coupled_subspaces_helper(coupling_in, &
   & degenerate_modes,coupling_order) result(output)
  implicit none
  
  type(CoupledSubspaces), intent(in)  :: coupling_in
  type(DegenerateModes),  intent(in)  :: degenerate_modes(:)
  integer,                intent(in)  :: coupling_order
  type(CoupledSubspaces), allocatable :: output(:)
  
  type(CoupledSubspaces) :: coupling
  
  integer :: i
  
  output = [CoupledSubspaces::]
  do i=1,size(degenerate_modes)
    coupling = coupling_in//degenerate_modes(i)
    if (size(coupling)==coupling_order) then
      output = [output, coupling]
    else
      output = [ output,                            &
             &   generate_coupled_subspaces_helper( &
             &              coupling,               &
             &              degenerate_modes(i+1:), &
             &              coupling_order)         &
             & ]
    endif
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
recursive function to_String_CoupledSubspaces(this) result(output)
  implicit none
  
  class(CoupledSubspaces), intent(in) :: this
  type(String)                        :: output
  
  output = join(this%ids)
end function
end module
