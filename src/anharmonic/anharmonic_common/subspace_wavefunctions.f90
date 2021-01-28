! ======================================================================
! An abstract class to hold a printable representation of the wavefunctions
!    spanning a subspace.
! ======================================================================
module caesar_subspace_wavefunctions_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: SubspaceWavefunctions
  public :: SubspaceWavefunctionsPointer
  
  type, abstract, extends(Stringsable) :: SubspaceWavefunctions
  contains
    procedure(representation_SubspaceWavefunctions), public, deferred, &
       & nopass :: representation
    procedure, public :: startup => startup_SubspaceWavefunctions
  end type
  
  type, extends(SubspaceWavefunctions) :: SubspaceWavefunctionsPointer
    type(String),                              private :: representation_
    class(SubspaceWavefunctions), allocatable, private :: wavefunctions_
  contains
    procedure, public :: check => check_SubspaceWavefunctionsPointer
    
    procedure, public, nopass :: representation => &
                               & representation_SubspaceWavefunctionsPointer
    
    ! I/O.
    procedure, public :: read  => read_SubspaceWavefunctionsPointer
    procedure, public :: write => write_SubspaceWavefunctionsPointer
  end type
  
  ! An array of all types which extend SubspaceWavefunctions.
  ! This array will be filled in by startup routines.
  type(SubspaceWavefunctionsPointer), allocatable :: &
     & TYPES_SubspaceWavefunctions(:)
  
  ! Abstract SubspaceWavfunction routines.
  abstract interface
    impure elemental function representation_SubspaceWavefunctions() &
       & result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
  end interface
  
  ! SubspaceWavefunctionPointer routines.
  interface SubspaceWavefunctionsPointer
    module procedure new_SubspaceWavefunctionsPointer
    module procedure new_SubspaceWavefunctionsPointer_Strings
    module procedure new_SubspaceWavefunctionsPointer_StringArray
  end interface
contains

! The startup method.
subroutine startup_SubspaceWavefunctions(this)
  implicit none
  
  class(SubspaceWavefunctions), intent(in) :: this
  
  integer :: i
  
  if (.not. allocated(TYPES_SubspaceWavefunctions)) then
    TYPES_SubspaceWavefunctions = [SubspaceWavefunctionsPointer(this)]
  elseif (.not. any([(                          &
     &    this%representation()                 &
     & == TYPES_SubspaceWavefunctions(i         &
     &       )%wavefunctions_%representation(), &
     & i=1,                                     &
     & size(TYPES_SubspaceWavefunctions)        )])) then
    TYPES_SubspaceWavefunctions = [ TYPES_SubspaceWavefunctions,       &
                                  & SubspaceWavefunctionsPointer(this) ]
  endif
end subroutine

! Construct a SubspaceWavefunctionsPointer from any type which extends
!    SubspaceWavefunctions.
impure elemental function new_SubspaceWavefunctionsPointer(wavefunctions) &
   & result(this)
  implicit none
  
  class(SubspaceWavefunctions), intent(in) :: wavefunctions
  type(SubspaceWavefunctionsPointer)       :: this
  
  integer :: ialloc
  
  select type(wavefunctions); type is (SubspaceWavefunctionsPointer)
    this = wavefunctions
  class default
    this%representation_ = wavefunctions%representation()
    allocate( this%wavefunctions_, source=wavefunctions, &
            & stat=ialloc); call err(ialloc)
  end select
end function

! Checks that the pointer has been allocated before it is used.
subroutine check_SubspaceWavefunctionsPointer(this)
  implicit none
  
  class(SubspaceWavefunctionsPointer), intent(in) :: this
  
  if (.not. allocated(this%wavefunctions_)) then
    call print_line(CODE_ERROR//': Trying to use a &
       &SubspaceWavefunctionsPointer before it has been allocated.')
    call err()
  endif
end subroutine

! Type representation.
impure elemental function representation_SubspaceWavefunctionsPointer() &
   & result(output)
  implicit none
  
  type(String) :: output
  
  output = 'pointer'
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_SubspaceWavefunctionsPointer(this,input)
  implicit none
  
  class(SubspaceWavefunctionsPointer), intent(out) :: this
  type(String),                        intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(SubspaceWavefunctionsPointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    i = first([( TYPES_SubspaceWavefunctions(i                         &
               &    )%wavefunctions_%representation()==representation, &
               & i=1,                                                  &
               & size(TYPES_SubspaceWavefunctions)                     )])
    
    ! Read the input into the element of the correct type,
    !    and copy that element into the output.
    call TYPES_SubspaceWavefunctions(i)%wavefunctions_%read(input(2:))
    this = SubspaceWavefunctionsPointer(TYPES_SubspaceWavefunctions(i))
  class default
    call err()
  end select
end subroutine

function write_SubspaceWavefunctionsPointer(this) result(output)
  implicit none
  
  class(SubspaceWavefunctionsPointer), intent(in) :: this
  type(String), allocatable                       :: output(:)
  
  select type(this); type is(SubspaceWavefunctionsPointer)
    output = [ 'Wavefunction representation: '//this%representation_, &
             & str(this%wavefunctions_)                               ]
  end select
end function

function new_SubspaceWavefunctionsPointer_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)           :: input(:)
  type(SubspaceWavefunctionsPointer) :: this
  
  call this%read(input)
end function

impure elemental function new_SubspaceWavefunctionsPointer_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in)      :: input
  type(SubspaceWavefunctionsPointer) :: this
  
  this = SubspaceWavefunctionsPointer(str(input))
end function
end module
