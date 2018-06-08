! ======================================================================
! A vector in complex mode co-ordinates.
! ======================================================================
module complex_mode_vector_submodule
  use utils_module
  
  use structure_module
  
  use complex_mode_submodule
  use complex_single_mode_vector_submodule
  implicit none
  
  private
  
  public :: ComplexModeVector
  public :: size
  
  type, extends(Stringsable) :: ComplexModeVector
    type(ComplexSingleModeVector), allocatable :: vectors(:)
  contains 
    procedure, public :: modes   => modes_ComplexModeVector
    procedure, public :: qpoints => qpoints_ComplexModeVector
    
    procedure, public :: read  => read_ComplexModeVector
    procedure, public :: write => write_ComplexModeVector
  end type
  
  interface ComplexModeVector
    module procedure new_ComplexModeVector
    module procedure new_ComplexModeVector_StringArray
  end interface
  
  interface size
    module procedure size_ComplexModeVector
  end interface
contains

! Constructor.
function new_ComplexModeVector(vectors) result(this)
  implicit none
  
  type(ComplexSingleModeVector), intent(in) :: vectors(:)
  type(ComplexModeVector)                   :: this
  
  this%vectors = vectors
end function

! Return the number of modes along which the vector has vectors.
function size_ComplexModeVector(input) result(output)
  implicit none
  
  class(ComplexModeVector), intent(in) :: input
  integer                              :: output
  
  output = size(input%vectors)
end function

! Returns a list of the modes at which the vector is non-zero.
function modes_ComplexModeVector(this,modes) result(output)
  implicit none
  
  class(ComplexModeVector), intent(in) :: this
  type(ComplexMode),        intent(in) :: modes(:)
  type(ComplexMode), allocatable       :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    output(i) = modes(first(modes%id==this%vectors(i)%id))
  enddo
end function

! Returns a list of the q-points at which the vector is non-zero.
function qpoints_ComplexModeVector(this,modes,qpoints) result(output)
  implicit none
  
  class(ComplexModeVector), intent(in) :: this
  type(ComplexMode),        intent(in) :: modes(:)
  type(QpointData),         intent(in) :: qpoints(:)
  type(QpointData), allocatable        :: output(:)
  
  type(ComplexMode), allocatable :: non_zero_modes(:)
  
  integer, allocatable :: qpoint_ids(:)
  
  integer :: i,j,ialloc
  
  ! List the q-point IDs of the modes in the vector.
  non_zero_modes = this%modes(modes)
  qpoint_ids = non_zero_modes%qpoint_id
  
  ! De-duplicate the q-point IDs.
  qpoint_ids = qpoint_ids(set(qpoint_ids))
  
  ! List the q-points matching the IDs.
  allocate(output(size(qpoint_ids)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    output(i) = qpoints(first(qpoints%id==qpoint_ids(i)))
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_ComplexModeVector(this,input)
  implicit none
  
  class(ComplexModeVector), intent(out) :: this
  type(String),             intent(in)  :: input(:)
  
  select type(this); type is(ComplexModeVector)
    this = ComplexModeVector(ComplexSingleModeVector(input))
  end select
end subroutine

function write_ComplexModeVector(this) result(output)
  implicit none
  
  class(ComplexModeVector), intent(in) :: this
  type(String), allocatable            :: output(:)
  
  select type(this); type is(ComplexModeVector)
    output = str(this%vectors)
  end select
end function

impure elemental function new_ComplexModeVector_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(ComplexModeVector)       :: this
  
  this = input
end function
end module
