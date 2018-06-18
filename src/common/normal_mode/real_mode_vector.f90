! ======================================================================
! A vector in real mode co-ordinates.
! ======================================================================
module real_mode_vector_submodule
  use utils_module
  
  use structure_module
  
  use cartesian_vector_submodule
  use real_mode_submodule
  use real_single_mode_vector_submodule
  implicit none
  
  private
  
  public :: RealModeVector
  public :: size
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
  type, extends(Stringsable) :: RealModeVector
    type(RealSingleModeVector), allocatable :: vectors(:)
  contains
    procedure, public :: modes   => modes_RealModeVector
    procedure, public :: qpoints => qpoints_RealModeVector
    
    procedure, public :: cartesian_vector => &
       & cartesian_vector_RealModeVector
    
    procedure, public :: read  => read_RealModeVector
    procedure, public :: write => write_RealModeVector
  end type
  
  interface RealModeVector
    module procedure new_RealModeVector
    module procedure new_RealModeVector_CartesianVector
    module procedure new_RealModeVector_StringArray
  end interface
  
  interface size
    module procedure size_RealModeVector
  end interface
  
  interface operator(*)
    module procedure multiply_real_RealModeVector
    module procedure multiply_RealModeVector_real
  end interface
  
  interface operator(/)
    module procedure divide_RealModeVector_real
  end interface
  
  interface operator(+)
    module procedure add_RealModeVector_RealModeVector
  end interface
  
  interface sum
    module procedure sum_RealModeVectors
  end interface
  
  interface operator(-)
    module procedure negative_RealModeVector
    module procedure subtract_RealModeVector_RealModeVector
  end interface
contains

! Constructor.
function new_RealModeVector(vectors) result(this)
  implicit none
  
  type(RealSingleModeVector), intent(in) :: vectors(:)
  type(RealModeVector)                   :: this
  
  this%vectors = vectors
end function

! Return the number of modes along which the vector has vectors.
function size_RealModeVector(input) result(output)
  implicit none
  
  class(RealModeVector), intent(in) :: input
  integer                           :: output
  
  output = size(input%vectors)
end function

! Returns a list of the modes at which the vector is non-zero.
function modes_RealModeVector(this,modes) result(output)
  implicit none
  
  class(RealModeVector), intent(in) :: this
  type(RealMode),        intent(in) :: modes(:)
  type(RealMode), allocatable       :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    output(i) = modes(first(modes%id==this%vectors(i)%id))
  enddo
end function

! Returns a list of the q-points at which the vector is non-zero.
function qpoints_RealModeVector(this,modes,qpoints) result(output)
  implicit none
  
  class(RealModeVector), intent(in) :: this
  type(RealMode),        intent(in) :: modes(:)
  type(QpointData),      intent(in) :: qpoints(:)
  type(QpointData), allocatable     :: output(:)
  
  type(RealMode), allocatable :: non_zero_modes(:)
  
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
impure elemental function multiply_real_RealModeVector(this,that) &
   & result(output)
  implicit none
  
  real(dp),             intent(in) :: this
  type(RealModeVector), intent(in) :: that
  type(RealModeVector)             :: output
  
  output = RealModeVector(this*that%vectors)
end function

impure elemental function multiply_RealModeVector_real(this,that) &
   & result(output)
  implicit none
  
  type(RealModeVector), intent(in) :: this
  real(dp),             intent(in) :: that
  type(RealModeVector)             :: output
  
  output = RealModeVector(this%vectors*that)
end function

impure elemental function divide_RealModeVector_real(this,that) &
   & result(output)
  implicit none
  
  type(RealModeVector), intent(in) :: this
  real(dp),             intent(in) :: that
  type(RealModeVector)             :: output
  
  output = RealModeVector(this%vectors/that)
end function

impure elemental function add_RealModeVector_RealModeVector(this,that) &
   & result(output)
  implicit none
  
  type(RealModeVector), intent(in) :: this
  type(RealModeVector), intent(in) :: that
  type(RealModeVector)             :: output
  
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

function sum_RealModeVectors(this) result(output)
  implicit none
  
  type(RealModeVector), intent(in) :: this(:)
  type(RealModeVector)             :: output
  
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

impure elemental function negative_RealModeVector(this) result(output)
  implicit none
  
  type(RealModeVector), intent(in) :: this
  type(RealModeVector)             :: output
  
  output = RealModeVector(-this%vectors)
end function

impure elemental function subtract_RealModeVector_RealModeVector(this,that) &
   & result(output)
  implicit none
  
  type(RealModeVector), intent(in) :: this
  type(RealModeVector), intent(in) :: that
  type(RealModeVector)             :: output
  
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
! Conversions between CartesianVector and RealModeVector.
! ----------------------------------------------------------------------
! Returns the vector in cartesian co-ordinates.
function cartesian_vector_RealModeVector(this,structure, &
   & modes,qpoints) result(output)
  implicit none
  
  class(RealModeVector), intent(in) :: this
  type(StructureData),   intent(in) :: structure
  type(RealMode),        intent(in) :: modes(:)
  type(QpointData),      intent(in) :: qpoints(:)
  type(CartesianVector)             :: output
  
  type(RealMode)   :: mode
  type(QpointData) :: qpoint
  
  type(CartesianVector), allocatable :: vectors(:)
  
  integer :: i,ialloc
  
  ! Calculate the cartesian vector due to each mode.
  allocate(vectors(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    ! Find the mode and q-point associated with vector i.
    mode = modes(first(modes%id==this%vectors(i)%id))
    qpoint = qpoints(first(qpoints%id==mode%qpoint_id))
    
    ! Calculate the vector from mode i.
    vectors(i) =                                      &
       & this%vectors(i)%cartesian_vector( mode,      &
       &                                   structure, &
       &                                   qpoint)
  enddo
  
  ! Add together the contributions from all modes.
  output = sum(vectors)
end function

! Converts a CartesianVector to a RealModeVector.
function new_RealModeVector_CartesianVector(vector, &
   & structure,modes,qpoints) result(this)
  implicit none
  
  class(CartesianVector), intent(in) :: vector
  type(StructureData),    intent(in) :: structure
  type(RealMode),         intent(in) :: modes(:)
  type(QpointData),       intent(in) :: qpoints(:)
  type(RealModeVector)               :: this
  
  type(RealSingleModeVector), allocatable :: vectors(:)
  type(QpointData)                        :: qpoint
  
  integer :: i,ialloc
  
  allocate(vectors(size(modes)), stat=ialloc); call err(ialloc)
  do i=1,size(modes)
    qpoint = qpoints(first(qpoints%id==modes(i)%qpoint_id))
    vectors(i) = RealSingleModeVector( modes(i),  &
                                     & vector,    &
                                     & structure, &
                                     & qpoint)
  enddo
  
  this = RealModeVector(vectors)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_RealModeVector(this,input)
  implicit none
  
  class(RealModeVector), intent(out) :: this
  type(String),          intent(in)  :: input(:)
  
  select type(this); type is(RealModeVector)
    this = RealModeVector(RealSingleModeVector(input))
  end select
end subroutine

function write_RealModeVector(this) result(output)
  implicit none
  
  class(RealModeVector), intent(in) :: this
  type(String), allocatable         :: output(:)
  
  select type(this); type is(RealModeVector)
    output = str(this%vectors)
  end select
end function

impure elemental function new_RealModeVector_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(RealModeVector)          :: this
  
  this = input
end function
end module
