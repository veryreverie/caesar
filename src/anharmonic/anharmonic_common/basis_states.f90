! ======================================================================
! A set of basis states, defined in terms of a SubspaceBasis.
! ======================================================================
module basis_states_module
  use common_module
  implicit none
  
  private
  
  public :: BasisStates
  public :: BasisStatesPointer
  
  type, abstract, extends(Stringsable) :: BasisStates
    integer :: subspace_id
  contains
    procedure(representation_BasisStates), public, deferred, nopass :: &
       & representation
    procedure, public :: startup => startup_BasisStates
  end type
  
  type, extends(BasisStates) :: BasisStatesPointer
    type(String),                    private :: representation_
    class(BasisStates), allocatable, private :: states_
  contains
    procedure, private :: check => check_BasisStatesPointer
    
    procedure, public, nopass :: representation => &
                               & representation_BasisStatesPointer
    
    procedure, public :: states => states_BasisStatesPointer
    
    ! I/O.
    procedure, public :: read  => read_BasisStatesPointer
    procedure, public :: write => write_BasisStatesPointer
  end type
  
  ! An array of all types which extend BasisStates.
  ! This array will be filled in by startup routines.
  type(BasisStatesPointer), allocatable :: TYPES_BasisStates(:)
  
  abstract interface
    impure elemental function representation_BasisStates() result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
  end interface
  
  interface BasisStatesPointer
    module procedure new_BasisStatesPointer
    module procedure new_BasisStatesPointer_Strings
    module procedure new_BasisStatesPointer_StringArray
  end interface
contains

subroutine startup_BasisStates(this)
  implicit none
  
  class(BasisStates), intent(in) :: this
  
  integer :: i
  
  if (.not.allocated(TYPES_BasisStates)) then
    TYPES_BasisStates = [BasisStatesPointer(this)]
  elseif (.not.any([(                                       &
     & this%representation()                                &
     &    == TYPES_BasisStates(i)%states_%representation(), &
     & i=1,                                                 &
     & size(TYPES_BasisStates)                              )])) then
    TYPES_BasisStates = [TYPES_BasisStates, BasisStatesPointer(this)]
  endif
end subroutine

! ----------------------------------------------------------------------
! BasisStatesPointer methods.
! ----------------------------------------------------------------------
! Construct a BasisStatesPointer from any type which extends BasisStates.
impure elemental function new_BasisStatesPointer(states) result(this)
  implicit none
  
  class(BasisStates), intent(in) :: states
  type(BasisStatesPointer)       :: this
  
  integer :: ialloc
  
  select type(states); type is(BasisStatesPointer)
    this = states
  class default
    this%representation_ = states%representation()
    allocate( this%states_, source=states, &
            & stat=ialloc); call err(ialloc)
    this%subspace_id = this%states_%subspace_id
  end select
end function

! Checks that the pointer has been allocated before it is used.
subroutine check_BasisStatesPointer(this)
  implicit none
  
  class(BasisStatesPointer), intent(in) :: this
  
  if (.not. allocated(this%states_)) then
    call print_line(CODE_ERROR//': Trying to use a BasisStatesPointer &
       &before it has been allocated.')
    call err()
  endif
end subroutine

! Type representation.
impure elemental function representation_BasisStatesPointer() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'pointer'
end function

! BasisStates methods.
function states_BasisStatesPointer(this) result(output)
  implicit none
  
  class(BasisStatesPointer), intent(in) :: this
  class(BasisStates), allocatable       :: output
  
  call this%check()
  
  output = this%states_
end function

! I/O.
subroutine read_BasisStatesPointer(this,input)
  implicit none
  
  class(BasisStatesPointer), intent(out) :: this
  type(String),              intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(BasisStatesPointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    i = first([(                                                        &
       & TYPES_BasisStates(i)%states_%representation()==representation, &
       & i=1,                                                           &
       & size(TYPES_BasisStates)                                        )])
    
    ! Read the input into the element of the correct type,
    !    and copy that element into the output.
    call TYPES_BasisStates(i)%states_%read(input(2:))
    this = BasisStatesPointer(TYPES_BasisStates(i))
  class default
    call err()
  end select
end subroutine

function write_BasisStatesPointer(this) result(output)
  implicit none
  
  class(BasisStatesPointer), intent(in) :: this
  type(String), allocatable             :: output(:)
  
  select type(this); type is(BasisStatesPointer)
    output = [ 'BasisStates representation: '//this%representation_, &
             & str(this%states_)                                     ]
  end select
end function

function new_BasisStatesPointer_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(BasisStatesPointer) :: this
  
  call this%read(input)
end function

impure elemental function new_BasisStatesPointer_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(BasisStatesPointer)      :: this
  
  this = BasisStatesPointer(str(input))
end function
end module
