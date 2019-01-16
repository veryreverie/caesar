! ======================================================================
! A wrapped polymorphic pointer to a wavefunction.
! ======================================================================
! Wraps all of SubspaceWavefunctions' methods,
!    calling them on the pointed-to wavefunction.
module subspace_wavefunctions_pointer_module
  use common_module
  
  use subspace_wavefunctions_module
  use full_subspace_wavefunctions_module
  implicit none
  
  private
  
  public :: SubspaceWavefunctionsPointer
  
  type, extends(SubspaceWavefunctions) :: SubspaceWavefunctionsPointer
    type(String)                              :: representation
    class(SubspaceWavefunctions), allocatable :: wavefunctions
  contains
    procedure, public :: check => check_SubspaceWavefunctionsPointer
    
    ! I/O.
    procedure, public :: read  => read_SubspaceWavefunctionsPointer
    procedure, public :: write => write_SubspaceWavefunctionsPointer
  end type
  
  interface SubspaceWavefunctionsPointer
    module procedure new_SubspaceWavefunctionsPointer
    module procedure new_SubspaceWavefunctionsPointer_Strings
    module procedure new_SubspaceWavefunctionsPointer_StringArray
  end interface
contains

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
  type is(FullSubspaceWavefunctions)
    this%representation = 'full'
    allocate( this%wavefunctions, source=wavefunctions, &
            & stat=ialloc); call err(ialloc)
  class default
    call err()
  end select
end function

! Checks that the pointer has been allocated before it is used.
subroutine check_SubspaceWavefunctionsPointer(this)
  implicit none
  
  class(SubspaceWavefunctionsPointer), intent(in) :: this
  
  if (.not. allocated(this%wavefunctions)) then
    call print_line(CODE_ERROR//': Trying to use a &
       &SubspaceWavefunctionsPointer before it has been allocated.')
    call err()
  endif
end subroutine

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_SubspaceWavefunctionsPointer(this,input)
  implicit none
  
  class(SubspaceWavefunctionsPointer), intent(out) :: this
  type(String),                        intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  select type(this); type is(SubspaceWavefunctionsPointer)
    line = split_line(input(1))
    representation = line(3)
    if (representation=='vscf') then
      this = SubspaceWavefunctionsPointer(FullSubspaceWavefunctions(input(2:)))
    else
      call print_line( 'Unrecognised wavefunction representation: '// &
                     & representation)
    endif
  class default
    call err()
  end select
end subroutine

function write_SubspaceWavefunctionsPointer(this) result(output)
  implicit none
  
  class(SubspaceWavefunctionsPointer), intent(in) :: this
  type(String), allocatable                       :: output(:)
  
  select type(this); type is(SubspaceWavefunctionsPointer)
    output = [ 'Wavefunction representation: '//this%representation, &
             & str(this%wavefunctions)                               ]
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
