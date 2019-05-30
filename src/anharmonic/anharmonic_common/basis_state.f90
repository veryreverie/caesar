! ======================================================================
! ======================================================================
module basis_state_module
  use common_module
  implicit none
  
  private
  
  public :: BasisState
  public :: BasisStatePointer
  
  type, abstract, extends(Stringsable) :: BasisState
    integer :: subspace_id
  contains
    procedure(representation_BasisState), public, deferred, nopass :: &
       & representation
    procedure, public :: startup => startup_BasisState
  end type
  
  type, extends(BasisState) :: BasisStatePointer
    type(String),                   private :: representation_
    class(BasisState), allocatable, private :: state_
  contains
    procedure, private :: check => check_BasisStatePointer
    
    procedure, public, nopass :: representation => &
                               & representation_BasisStatePointer
    
    procedure, public :: state => state_BasisStatePointer
    
    ! I/O.
    procedure, public :: read  => read_BasisStatePointer
    procedure, public :: write => write_BasisStatePointer
  end type
  
  ! An array of all types which extend BasisState.
  ! This array will be filled in by startup routines.
  type(BasisStatePointer), allocatable :: TYPES_BasisState(:)
  
  abstract interface
    impure elemental function representation_BasisState() result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
  end interface
  
  interface BasisStatePointer
    module procedure new_BasisStatePointer
    module procedure new_BasisStatePointer_Strings
    module procedure new_BasisStatePointer_StringArray
  end interface
contains

subroutine startup_BasisState(this)
  implicit none
  
  class(BasisState), intent(in) :: this
  
  integer :: i
  
  if (.not.allocated(TYPES_BasisState)) then
    TYPES_BasisState = [BasisStatePointer(this)]
  elseif (.not.any([(                                     &
     & this%representation()                              &
     &    == TYPES_BasisState(i)%state_%representation(), &
     & i=1,                                               &
     & size(TYPES_BasisState)                             )])) then
    TYPES_BasisState = [TYPES_BasisState, BasisStatePointer(this)]
  endif
end subroutine

! ----------------------------------------------------------------------
! BasisStatePointer methods.
! ----------------------------------------------------------------------
! Construct a BasisStatePointer from any type which extends BasisState.
impure elemental function new_BasisStatePointer(state) result(this)
  implicit none
  
  class(BasisState), intent(in) :: state
  type(BasisStatePointer)       :: this
  
  integer :: ialloc
  
  select type(state); type is(BasisStatePointer)
    this = state
  class default
    this%representation_ = state%representation()
    allocate( this%state_, source=state, &
            & stat=ialloc); call err(ialloc)
  end select
end function

! Checks that the pointer has been allocated before it is used.
subroutine check_BasisStatePointer(this)
  implicit none
  
  class(BasisStatePointer), intent(in) :: this
  
  if (.not. allocated(this%state_)) then
    call print_line(CODE_ERROR//': Trying to use a BasisStatePointer &
       &before it has been allocated.')
    call err()
  endif
end subroutine

! Type representation.
impure elemental function representation_BasisStatePointer() &
   & result(output)
  implicit none
  
  type(String) :: output
  
  output = 'pointer'
end function

! BasisState methods.
function state_BasisStatePointer(this) result(output)
  implicit none
  
  class(BasisStatePointer), intent(in) :: this
  class(BasisState), allocatable       :: output
  
  output = this%state_
end function

! I/O.
subroutine read_BasisStatePointer(this,input)
  implicit none
  
  class(BasisStatePointer), intent(out) :: this
  type(String),             intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(BasisStatePointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    i = first([(                                                      &
       & TYPES_BasisState(i)%state_%representation()==representation, &
       & i=1,                                                         &
       & size(TYPES_BasisState)                                       )])
    
    ! Read the input into the element of the correct type,
    !    and copy that element into the output.
    call TYPES_BasisState(i)%state_%read(input(2:))
    this = BasisStatePointer(TYPES_BasisState(i))
  class default
    call err()
  end select
end subroutine

function write_BasisStatePointer(this) result(output)
  implicit none
  
  class(BasisStatePointer), intent(in) :: this
  type(String), allocatable            :: output(:)
  
  select type(this); type is(BasisStatePointer)
    output = [ 'BasisState representation: '//this%representation_, &
             & str(this%state_)                                     ]
  end select
end function

function new_BasisStatePointer_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(BasisStatePointer)  :: this
  
  call this%read(input)
end function

impure elemental function new_BasisStatePointer_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(BasisStatePointer)       :: this
  
  this = BasisStatePointer(str(input))
end function
end module
