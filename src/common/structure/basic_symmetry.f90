! ======================================================================
! A minimal representation of the SymmetryOperator class.
! ======================================================================
module basic_symmetry_submodule
  use utils_module
  implicit none
  
  private
  
  public :: BasicSymmetry
  
  type, extends(Stringsable) :: BasicSymmetry
    integer                      :: id
    type(IntMatrix)              :: tensor
    type(RealVector)             :: translation
    type(Group)                  :: atom_group
    type(IntVector), allocatable :: rvectors(:)
  contains
    procedure, public :: read  => read_BasicSymmetry
    procedure, public :: write => write_BasicSymmetry
  end type
  
  interface BasicSymmetry
    module procedure new_BasicSymmetry
    module procedure new_BasicSymmetry_Strings
    module procedure new_BasicSymmetry_StringArray
  end interface
contains

! Constructor.
function new_BasicSymmetry(id,tensor,translation,atom_group,rvectors) &
   & result(this)
  implicit none
  
  integer,          intent(in) :: id
  type(IntMatrix),  intent(in) :: tensor
  type(RealVector), intent(in) :: translation
  type(Group),      intent(in) :: atom_group
  type(IntVector),  intent(in) :: rvectors(:)
  type(BasicSymmetry)          :: this
  
  this%id          = id
  this%tensor      = tensor
  this%translation = translation
  this%atom_group  = atom_group
  this%rvectors    = rvectors
end function

! I/O.
subroutine read_BasicSymmetry(this,input)
  implicit none
  
  class(BasicSymmetry), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  integer                      :: id
  type(IntMatrix)              :: tensor
  type(RealVector)             :: translation
  integer,         allocatable :: atom_group(:)
  type(IntVector), allocatable :: rvectors(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(BasicSymmetry)
    line = split_line(input(1))
    id = int(line(2))
    tensor = IntMatrix(input(3:5))
    translation = RealVector(input(7))
    
    allocate( atom_group(size(input)-8), &
            & rvectors(size(input)-8),   &
            & stat=ialloc); call err(ialloc)
    do i=1,size(atom_group)
      line = split_line(input(8+i))
      atom_group(i) = int(line(5))
      rvectors(i) = vec(int(line(8:10)))
    enddo
    
    this = BasicSymmetry(id,tensor,translation,Group(atom_group),rvectors)
  class default
    call err()
  end select
end subroutine

function write_BasicSymmetry(this) result(output)
  implicit none
  
  class(BasicSymmetry), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  type(String), allocatable :: atom_strings(:)
  
  integer :: i,ialloc
  
  allocate( atom_strings(size(this%rvectors)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this%rvectors)
    atom_strings(i) = 'Atom '//i//' -> Atom '//this%atom_group*i// &
                    & ' + R-vector '//this%rvectors(i)
  enddo
  
  select type(this); type is(BasicSymmetry)
    output = [ 'Operation '//this%id,                          &
             & str('Tensor in fractional co-ordinates:'),      &
             & str(this%tensor),                               &
             & str('Translation in fractional co-ordinates:'), &
             & str(this%translation),                          &
             & str('Effect on atoms:'),                        &
             & atom_strings                                    ]
  class default
    call err()
  end select
end function

function new_BasicSymmetry_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(BasicSymmetry)      :: this
  
  call this%read(input)
end function

impure elemental function new_BasicSymmetry_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(BasicSymmetry)           :: this
  
  this = BasicSymmetry(str(input))
end function
end module
