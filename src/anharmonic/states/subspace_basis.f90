! ======================================================================
! A basis of states which spans a subspace.
! ======================================================================
module subspace_basis_module
  use common_module
  
  use subspace_state_module
  implicit none
  
  private
  
  public :: SubspaceBasis
  public :: size
  public :: generate_subspace_basis
  
  type, extends(Stringsable) :: SubspaceBasis
    integer                          :: subspace_id
    real(dp)                         :: frequency
    type(SubspaceState), allocatable :: states(:)
  contains
    procedure, public :: read  => read_SubspaceBasis
    procedure, public :: write => write_SubspaceBasis
  end type
  
  interface SubspaceBasis
    module procedure new_SubspaceBasis
    module procedure new_SubspaceBasis_Strings
    module procedure new_SubspaceBasis_StringArray
  end interface
  
  interface size
    module procedure size_SubspaceBasis
  end interface
contains

! Constructor and size function.
function new_SubspaceBasis(subspace_id,frequency,states) result(this)
  implicit none
  
  integer,             intent(in) :: subspace_id
  real(dp),            intent(in) :: frequency
  type(SubspaceState), intent(in) :: states(:)
  type(SubspaceBasis)             :: this
  
  if (any(states%subspace_id/=subspace_id)) then
    call print_line(CODE_ERROR//": Subspaces don't match.")
    call err()
  endif
  
  this%subspace_id = subspace_id
  this%frequency   = frequency
  this%states      = states
end function

function size_SubspaceBasis(this) result(output)
  implicit none
  
  type(SubspaceBasis), intent(in) :: this
  integer                         :: output
  
  output = size(this%states)
end function

! Generates states up to a given power.
function generate_subspace_basis(subspace,frequency,modes,maximum_power) &
   & result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency
  type(ComplexMode),        intent(in) :: modes(:)
  integer,                  intent(in) :: maximum_power
  type(SubspaceBasis)                  :: output
  
  type(SubspaceState), allocatable :: states(:)
  
  states = generate_subspace_states( subspace,     &
                                   & frequency,    &
                                   & modes,        &
                                   & maximum_power )
  output = SubspaceBasis(subspace%id, frequency, states)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_SubspaceBasis(this,input)
  implicit none
  
  class(SubspaceBasis), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  integer                            :: subspace_id
  real(dp)                           :: frequency
  type(ComplexMonomial)              :: state
  type(ComplexMonomial), allocatable :: monomials(:)
  type(SubspaceState),   allocatable :: states(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i
  
  select type(this); type is(SubspaceBasis)
    line = split_line(input(1))
    subspace_id = int(line(2))
    
    line = split_line(input(2))
    frequency = dble(line(2))
    
    monomials = [( ComplexMonomial(slice(input(i),1,len(input(i))-3)), &
                 & i=4,                                                &
                 & size(input)                                         )]
    states = [( SubspaceState(subspace_id,frequency,monomials(i)), &
              & i=1,                                               &
              & size(monomials)                                    )]
    
    this = SubspaceBasis(subspace_id, frequency, states)
  end select
end subroutine

function write_SubspaceBasis(this) result(output)
  implicit none
  
  class(SubspaceBasis), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  type(String), allocatable :: states(:)
  
  integer :: i
  
  select type(this); type is(SubspaceBasis)
    states = [( str(this%states(i)%state)//'|0>', i=1, size(this) )]
    output = [ 'Subspace '//this%subspace_id, &
             & 'Frequency '//this%frequency,  &
             & str('States'),                 &
             & states                         ]
  end select
end function

function new_SubspaceBasis_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(SubspaceBasis)      :: this
  
  call this%read(input)
end function

impure elemental function new_SubspaceBasis_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SubspaceBasis)           :: this
  
  this = SubspaceBasis(str(input))
end function
end module
