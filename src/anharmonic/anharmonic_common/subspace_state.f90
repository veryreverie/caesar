! ======================================================================
! The SubspaceState abstract class, defining states which do not need to
!    reference a basis.
! ======================================================================
! N.B. the Basis and State methods take other Basis and State class() arguments
!    with the TARGET attribute, so that they can be cast to their concrete type
!    pointers.
! The %state_pointer() and %states_pointer() methods should
!    not be called on objects which do not have the TARGET attribute,
!    as this will silently lead to undefined behaviour.
module caesar_subspace_state_module
  use caesar_common_module
  
  use caesar_stress_prefactors_module
  use caesar_anharmonic_data_module
  use caesar_sparse_monomial_module
  implicit none
  
  private
  
  public :: SubspaceState
  public :: SubspaceStatePointer
  
  type, abstract, extends(Stringsable) :: SubspaceState
  contains
    procedure(representation_SubspaceState), public, deferred, nopass :: &
       & representation
    procedure, public :: startup => startup_SubspaceState
    
    ! Return the id of the modes across which the state is defined.
    procedure(mode_ids_SubspaceState), public, deferred :: mode_ids
    procedure(paired_mode_ids_SubspaceState), public, deferred :: &
       & paired_mode_ids
    
    ! Returns the occupation of the state.
    procedure(occupation_SubspaceState), public, deferred :: occupation
    
    ! Returns the wavevector of the state.
    procedure(wavevector_SubspaceState), public, deferred :: wavevector
  end type
  
  type, extends(SubspaceState) :: SubspaceStatePointer
    type(String),                      private :: representation_
    ! N.B. state_ should never be modified.
    ! It is public for performance reasons only.
    class(SubspaceState), allocatable, public :: state_
  contains
    procedure, private :: check => check_SubspaceStatePointer
    
    procedure, public, nopass :: representation => &
                               & representation_SubspaceStatePointer
    
    procedure, public :: state         => state_SubspaceStatePointer
    procedure, public :: state_pointer => state_pointer_SubspaceStatePointer
    
    procedure, public :: mode_ids => &
                       & mode_ids_SubspaceStatePointer
    procedure, public :: paired_mode_ids => &
                       & paired_mode_ids_SubspaceStatePointer
    
    procedure, public :: occupation => occupation_SubspaceStatePointer
    
    procedure, public :: wavevector => wavevector_SubspaceStatePointer
    
    ! I/O.
    procedure, public :: read  => read_SubspaceStatePointer
    procedure, public :: write => write_SubspaceStatePointer
  end type
  
  ! An array of all types which extend SubspaceState.
  ! This array will be filled in by startup routines.
  type(SubspaceStatePointer), allocatable :: TYPES_SubspaceState(:)
  
  ! Abstract interface for SubspaceState functionality.
  abstract interface
    impure elemental function representation_SubspaceState() result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
    
    function mode_ids_SubspaceState(this) result(output)
      import SubspaceState
      implicit none
      
      class(SubspaceState), intent(in) :: this
      integer, allocatable             :: output(:)
    end function
    
    function paired_mode_ids_SubspaceState(this) result(output)
      import SubspaceState
      implicit none
      
      class(SubspaceState), intent(in) :: this
      integer, allocatable             :: output(:)
    end function
    
    impure elemental function occupation_SubspaceState(this) result(output)
      import SubspaceState
      implicit none
      
      class(SubspaceState), intent(in) :: this
      integer                          :: output
    end function
    
    function wavevector_SubspaceState(this,modes,qpoints) &
       & result(output)
      import SubspaceState
      import ComplexMode
      import QpointData
      import FractionVector
      implicit none
      
      class(SubspaceState), intent(in) :: this
      type(ComplexMode),    intent(in) :: modes(:)
      type(QpointData),     intent(in) :: qpoints(:)
      type(FractionVector)             :: output
    end function
  end interface
  
  interface SubspaceStatePointer
    module procedure new_SubspaceStatePointer
    module procedure new_SubspaceStatePointer_Strings
    module procedure new_SubspaceStatePointer_StringArray
  end interface
contains

! Startup method.
subroutine startup_SubspaceState(this)
  implicit none
  
  class(SubspaceState), intent(in) :: this
  
  integer :: i
  
  if (.not. allocated(TYPES_SubspaceState)) then
    TYPES_SubspaceState = [SubspaceStatePointer(this)]
  elseif (.not. any([(                                    &
     &    this%representation()                           &
     & == TYPES_SubspaceState(i)%state_%representation(), &
     & i=1,                                               &
     & size(TYPES_SubspaceState)                          )])) then
    TYPES_SubspaceState = [TYPES_SubspaceState, SubspaceStatePointer(this)]
  endif
end subroutine

! --------------------------------------------------
! SubspaceStatePointer methods.
! --------------------------------------------------
! Construct a SubspaceStatePointer from any type which extends SubspaceState.
impure elemental function new_SubspaceStatePointer(state) result(this)
  implicit none
  
  class(SubspaceState), intent(in) :: state
  type(SubspaceStatePointer)       :: this
  
  integer :: ialloc
  
  select type(state); type is(SubspaceStatePointer)
    this = state
  class default
    this%representation_ = state%representation()
    allocate( this%state_, source=state, &
            & stat=ialloc); call err(ialloc)
  end select
end function

! Checks that the pointer has been allocated before it is used.
subroutine check_SubspaceStatePointer(this)
  implicit none
  
  class(SubspaceStatePointer), intent(in) :: this
  
  if (.not. allocated(this%state_)) then
    call print_line(CODE_ERROR//': Trying to use a SubspaceStatePointer &
       &before it has been allocated.')
    call err()
  endif
end subroutine

! Type representation.
impure elemental function representation_SubspaceStatePointer() &
   & result(output)
  implicit none
  
  type(String) :: output
  
  output = 'pointer'
end function

! SubspaceState methods.
function state_SubspaceStatePointer(this) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in) :: this
  class(SubspaceState), allocatable       :: output
  
  output = this%state_
end function

function state_pointer_SubspaceStatePointer(this) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in), target :: this
  class(SubspaceState), pointer                   :: output
  
  output => this%state_
end function

! --------------------------------------------------
! SubspaceStatePointer wrappers for SubspaceState methods.
! --------------------------------------------------
function mode_ids_SubspaceStatePointer(this) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in) :: this
  integer, allocatable                    :: output(:)
  
  call this%check()
  
  output = this%state_%mode_ids()
end function

function paired_mode_ids_SubspaceStatePointer(this) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in) :: this
  integer, allocatable                    :: output(:)
  
  call this%check()
  
  output = this%state_%paired_mode_ids()
end function

impure elemental function occupation_SubspaceStatePointer(this) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in) :: this
  integer                                 :: output
  
  call this%check()
  
  output = this%state_%occupation()
end function

function wavevector_SubspaceStatePointer(this,modes,qpoints) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in) :: this
  type(ComplexMode),           intent(in) :: modes(:)
  type(QpointData),            intent(in) :: qpoints(:)
  type(FractionVector)                    :: output
  
  call this%check()
  
  output = this%state_%wavevector(modes,qpoints)
end function

! --------------------------------------------------
! SubspaceStatePointer I/O.
! --------------------------------------------------
subroutine read_SubspaceStatePointer(this,input)
  implicit none
  
  class(SubspaceStatePointer), intent(out) :: this
  type(String),                intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(SubspaceStatePointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    i = first([(                                                         &
       & TYPES_SubspaceState(i)%state_%representation()==representation, &
       & i=1,                                                            &
       & size(TYPES_SubspaceState)                                       )])
    
    ! Read the input into the element of the correct type,
    !    and copy that element into the output.
    call TYPES_SubspaceState(i)%state_%read(input(2:))
    this = SubspaceStatePointer(TYPES_SubspaceState(i))
  class default
    call err()
  end select
end subroutine

function write_SubspaceStatePointer(this) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in) :: this
  type(String), allocatable               :: output(:)
  
  select type(this); type is(SubspaceStatePointer)
    output = [ 'SubspaceState representation: '//this%representation_, &
             & str(this%state_)                                        ]
  end select
end function

function new_SubspaceStatePointer_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)   :: input(:)
  type(SubspaceStatePointer) :: this
  
  call this%read(input)
end function

impure elemental function new_SubspaceStatePointer_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SubspaceStatePointer)    :: this
  
  this = SubspaceStatePointer(str(input))
end function
end module
