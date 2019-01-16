! ======================================================================
! A wrapped polymorphic pointer to a subspace basis.
! ======================================================================
! Wraps all of SubspaceBasis's methods,
!    calling them on the pointed-to states.
module subspace_basis_pointer_module
  use common_module
  
  use anharmonic_common_module
  
  use full_subspace_basis_and_states_module
  implicit none
  
  private
  
  public :: SubspaceBasisPointer
  
  type, extends(SubspaceBasis) :: SubspaceBasisPointer
    type(String),                      private :: representation_
    class(SubspaceBasis), allocatable, private :: basis_
  contains
    procedure, public :: check => check_SubspaceBasisPointer
    
    procedure, public :: initial_states => initial_states_SubspaceBasisPointer
    procedure, public :: calculate_states => &
                       & calculate_states_SubspaceBasisPointer
    
    ! I/O.
    procedure, public :: read  => read_SubspaceBasisPointer
    procedure, public :: write => write_SubspaceBasisPointer
  end type
  
  interface SubspaceBasisPointer
    module procedure new_SubspaceBasisPointer
    module procedure new_SubspaceBasisPointer_Strings
    module procedure new_SubspaceBasisPointer_StringArray
  end interface
contains

! Construct a SubspaceBasisPointer from any type which extends SubspaceBasis.
impure elemental function new_SubspaceBasisPointer(basis) result(this)
  implicit none
  
  class(SubspaceBasis), intent(in) :: basis
  type(SubspaceBasisPointer)       :: this
  
  integer :: ialloc
  
  select type(basis); type is(SubspaceBasisPointer)
    this = basis
  type is(FullSubspaceBasis)
    this%representation_ = 'full'
    allocate( this%basis_, source=basis, &
            & stat=ialloc); call err(ialloc)
  class default
    call err()
  end select
end function

! Checks that the pointer has been allocated before it is used.
subroutine check_SubspaceBasisPointer(this)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  
  if (.not. allocated(this%basis_)) then
    call print_line(CODE_ERROR//': Trying to use a SubspaceBasisPointer &
       &before it has been allocated.')
    call err()
  endif
end subroutine

! ----------------------------------------------------------------------
! SubspaceBasis methods.
! ----------------------------------------------------------------------
impure elemental function initial_states_SubspaceBasisPointer(this,subspace, &
   & anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  type(DegenerateSubspace),    intent(in) :: subspace
  type(AnharmonicData),        intent(in) :: anharmonic_data
  class(SubspaceStates), allocatable      :: output
  
  call this%check()
  
  output = this%basis_%initial_states(subspace, anharmonic_data)
end function

impure elemental function calculate_states_SubspaceBasisPointer(this, &
   & subspace,subspace_potential,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  type(DegenerateSubspace),    intent(in) :: subspace
  class(PotentialData),        intent(in) :: subspace_potential
  type(AnharmonicData),        intent(in) :: anharmonic_data
  class(SubspaceStates), allocatable      :: output
  
  call this%check()
  
  output = this%basis_%calculate_states( subspace,           &
                                       & subspace_potential, &
                                       & anharmonic_data     )
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_SubspaceBasisPointer(this,input)
  implicit none
  
  class(SubspaceBasisPointer), intent(out) :: this
  type(String),                intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  select type(this); type is(SubspaceBasisPointer)
    line = split_line(input(1))
    representation = line(3)
    if (representation=='full') then
      this = SubspaceBasisPointer(FullSubspaceBasis(input(2:)))
    else
      call print_line( 'Unrecognised basis representation: '// &
                     & representation                          )
    endif
  class default
    call err()
  end select
end subroutine

function write_SubspaceBasisPointer(this) result(output)
  implicit none
  
  class(SubspaceBasisPointer), intent(in) :: this
  type(String), allocatable               :: output(:)
  
  select type(this); type is(SubspaceBasisPointer)
    output = [ 'SubspaceBasis representation: '//this%representation_, &
             & str(this%basis_)                                        ]
  end select
end function

function new_SubspaceBasisPointer_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)   :: input(:)
  type(SubspaceBasisPointer) :: this
  
  call this%read(input)
end function

impure elemental function new_SubspaceBasisPointer_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SubspaceBasisPointer)    :: this
  
  this = SubspaceBasisPointer(str(input))
end function
end module
