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
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
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
  
  interface operator(*)
    module procedure multiply_real_ComplexModeVector
    module procedure multiply_ComplexModeVector_real
    module procedure multiply_complex_ComplexModeVector
    module procedure multiply_ComplexModeVector_complex
  end interface
  
  interface operator(/)
    module procedure divide_ComplexModeVector_complex
  end interface
  
  interface operator(+)
    module procedure add_ComplexModeVector_ComplexModeVector
  end interface
  
  interface sum
    module procedure sum_ComplexModeVectors
  end interface
  
  interface operator(-)
    module procedure negative_ComplexModeVector
    module procedure subtract_ComplexModeVector_ComplexModeVector
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
  
  integer :: i,ialloc
  
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
! Arithmetic.
! ----------------------------------------------------------------------
impure elemental function multiply_real_ComplexModeVector(this,that) &
   & result(output)
  implicit none
  
  real(dp),                intent(in) :: this
  type(ComplexModeVector), intent(in) :: that
  type(ComplexModeVector)             :: output
  
  output = ComplexModeVector(this*that%vectors)
end function

impure elemental function multiply_ComplexModeVector_real(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeVector), intent(in) :: this
  real(dp),                intent(in) :: that
  type(ComplexModeVector)             :: output
  
  output = ComplexModeVector(this%vectors*that)
end function

impure elemental function multiply_complex_ComplexModeVector(this,that) &
   & result(output)
  implicit none
  
  complex(dp),             intent(in) :: this
  type(ComplexModeVector), intent(in) :: that
  type(ComplexModeVector)             :: output
  
  output = ComplexModeVector(this*that%vectors)
end function

impure elemental function multiply_ComplexModeVector_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeVector), intent(in) :: this
  complex(dp),             intent(in) :: that
  type(ComplexModeVector)             :: output
  
  output = ComplexModeVector(this%vectors*that)
end function

impure elemental function divide_ComplexModeVector_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeVector), intent(in) :: this
  complex(dp),             intent(in) :: that
  type(ComplexModeVector)             :: output
  
  output = ComplexModeVector(this%vectors/that)
end function

impure elemental function add_ComplexModeVector_ComplexModeVector(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeVector), intent(in) :: this
  type(ComplexModeVector), intent(in) :: that
  type(ComplexModeVector)             :: output
  
  integer :: i,j
  
  output = this
  do i=1,size(that)
    j = first(this%vectors%id==that%vectors(i)%id, default=0)
    if (j==0) then
      output%vectors = [output%vectors, that%vectors(i)]
    else
      output%vectors(j) = output%vectors(j) + that%vectors(i)
    endif
  enddo
end function

function sum_ComplexModeVectors(this) result(output)
  implicit none
  
  type(ComplexModeVector), intent(in) :: this(:)
  type(ComplexModeVector)             :: output
  
  integer :: i
  
  if (size(this)==0) then
    call print_line(ERROR//': Trying to sum an empty list.')
    call err()
  endif
  
  output = this(1)
  do i=2,size(this)
    output = output + this(i)
  enddo
end function

impure elemental function negative_ComplexModeVector(this) result(output)
  implicit none
  
  type(ComplexModeVector), intent(in) :: this
  type(ComplexModeVector)             :: output
  
  output = ComplexModeVector(-this%vectors)
end function

impure elemental function subtract_ComplexModeVector_ComplexModeVector(this, &
   & that) result(output)
  implicit none
  
  type(ComplexModeVector), intent(in) :: this
  type(ComplexModeVector), intent(in) :: that
  type(ComplexModeVector)             :: output
  
  integer :: i,j
  
  output = this
  do i=1,size(that)
    j = first(this%vectors%id==that%vectors(i)%id, default=0)
    if (j==0) then
      output%vectors = [output%vectors, -that%vectors(i)]
    else
      output%vectors(j) = output%vectors(j) - that%vectors(i)
    endif
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
