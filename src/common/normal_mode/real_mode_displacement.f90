! ======================================================================
! A displacement in real mode co-ordinates.
! ======================================================================
module real_mode_displacement_module
  use utils_module
  
  use structure_module
  
  use mass_weighted_displacement_module
  use cartesian_displacement_module
  use real_mode_module
  use real_single_mode_displacement_module
  implicit none
  
  private
  
  public :: RealModeDisplacement
  public :: size
  public :: MassWeightedDisplacement
  public :: CartesianDisplacement
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
  type, extends(Stringsable) :: RealModeDisplacement
    type(RealSingleDisplacement), allocatable :: vectors(:)
  contains
    ! The component of the displacement along a given mode.
    generic,   public  :: displacement =>  &
                        & displacement_id, &
                        & displacement_mode
    procedure, private :: displacement_id
    procedure, private :: displacement_mode
    
    ! I/O.
    procedure, public :: read  => read_RealModeDisplacement
    procedure, public :: write => write_RealModeDisplacement
  end type
  
  interface RealModeDisplacement
    module procedure new_RealModeDisplacement
    module procedure new_RealModeDisplacement_RealModes
    module procedure new_RealModeDisplacement_MassWeightedDisplacement
    module procedure new_RealModeDisplacement_CartesianDisplacement
    module procedure new_RealModeDisplacement_Strings
    module procedure new_RealModeDisplacement_StringArray
  end interface
  
  interface size
    module procedure size_RealModeDisplacement
  end interface
  
  interface MassWeightedDisplacement
    module procedure new_MassWeightedDisplacement_RealModeDisplacement
  end interface
  
  interface CartesianDisplacement
    module procedure new_CartesianDisplacement_RealModeDisplacement
  end interface
  
  interface operator(*)
    module procedure multiply_real_RealModeDisplacement
    module procedure multiply_RealModeDisplacement_real
  end interface
  
  interface operator(/)
    module procedure divide_RealModeDisplacement_real
  end interface
  
  interface operator(+)
    module procedure add_RealModeDisplacement_RealModeDisplacement
  end interface
  
  interface sum
    module procedure sum_RealModeDisplacements
  end interface
  
  interface operator(-)
    module procedure negative_RealModeDisplacement
    module procedure subtract_RealModeDisplacement_RealModeDisplacement
  end interface
contains

! Constructors and size() function.
function new_RealModeDisplacement(displacements) result(this)
  implicit none
  
  type(RealSingleDisplacement), intent(in) :: displacements(:)
  type(RealModeDisplacement)               :: this
  
  this%vectors = displacements
end function

function new_RealModeDisplacement_RealModes(modes,displacements) result(this)
  implicit none
  
  type(RealMode), intent(in) :: modes(:)
  real(dp),       intent(in) :: displacements(:)
  type(RealModeDisplacement) :: this
  
  this = RealModeDisplacement(RealSingleDisplacement(modes,displacements))
end function

function size_RealModeDisplacement(this) result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in) :: this
  integer                                :: output
  
  output = size(this%vectors)
end function

! ----------------------------------------------------------------------
! Arithmetic.
! ----------------------------------------------------------------------
impure elemental function multiply_real_RealModeDisplacement(this,that) &
   & result(output)
  implicit none
  
  real(dp),                   intent(in) :: this
  type(RealModeDisplacement), intent(in) :: that
  type(RealModeDisplacement)             :: output
  
  output = RealModeDisplacement(this*that%vectors)
end function

impure elemental function multiply_RealModeDisplacement_real(this,that) &
   & result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in) :: this
  real(dp),                   intent(in) :: that
  type(RealModeDisplacement)             :: output
  
  output = RealModeDisplacement(this%vectors*that)
end function

impure elemental function divide_RealModeDisplacement_real(this,that) &
   & result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in) :: this
  real(dp),                   intent(in) :: that
  type(RealModeDisplacement)             :: output
  
  output = RealModeDisplacement(this%vectors/that)
end function

impure elemental function add_RealModeDisplacement_RealModeDisplacement(this, &
   & that) result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in) :: this
  type(RealModeDisplacement), intent(in) :: that
  type(RealModeDisplacement)             :: output
  
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

function sum_RealModeDisplacements(this) result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in) :: this(:)
  type(RealModeDisplacement)             :: output
  
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

impure elemental function negative_RealModeDisplacement(this) result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in) :: this
  type(RealModeDisplacement)             :: output
  
  output = RealModeDisplacement(-this%vectors)
end function

impure elemental function subtract_RealModeDisplacement_RealModeDisplacement( &
   & this,that) result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in) :: this
  type(RealModeDisplacement), intent(in) :: that
  type(RealModeDisplacement)             :: output
  
  output = this + (-that)
end function

! ----------------------------------------------------------------------
! Conversions to and from cartesian and mass-weighted co-ordinates.
! ----------------------------------------------------------------------
! Returns the displacement in mass-weighted co-ordinates.
function new_MassWeightedDisplacement_RealModeDisplacement(this,structure, &
   & modes,qpoints) result(output)
  implicit none
  
  class(RealModeDisplacement), intent(in) :: this
  type(StructureData),         intent(in) :: structure
  type(RealMode),              intent(in) :: modes(:)
  type(QpointData),            intent(in) :: qpoints(:)
  type(MassWeightedDisplacement)          :: output
  
  type(RealMode),   allocatable :: selected_modes(:)
  type(QpointData), allocatable :: selected_qpoints(:)
  
  if (size(this)==0) then
    output = MassWeightedDisplacement(structure)
  else
    selected_modes = select_modes(this%vectors, modes)
    selected_qpoints = select_qpoints(selected_modes, qpoints)
    output = sum(MassWeightedDisplacement( this%vectors,    &
                                         & selected_modes,  &
                                         & structure,       &
                                         & selected_qpoints ))
  endif
end function

! Returns the displacement in cartesian co-ordinates.
function new_CartesianDisplacement_RealModeDisplacement(this,structure, &
   & modes,qpoints) result(output)
  implicit none
  
  class(RealModeDisplacement), intent(in) :: this
  type(StructureData),         intent(in) :: structure
  type(RealMode),              intent(in) :: modes(:)
  type(QpointData),            intent(in) :: qpoints(:)
  type(CartesianDisplacement)             :: output
  
  type(RealMode),   allocatable :: selected_modes(:)
  type(QpointData), allocatable :: selected_qpoints(:)
  
  if (size(this)==0) then
    output = CartesianDisplacement(structure)
  else
    selected_modes = select_modes(this%vectors, modes)
    selected_qpoints = select_qpoints(selected_modes, qpoints)
    output = sum(CartesianDisplacement( this%vectors,    &
                                      & selected_modes,  &
                                      & structure,       &
                                      & selected_qpoints ))
  endif
end function

! Converts a MassWeightedDisplacement to a RealModeDisplacement.
function new_RealModeDisplacement_MassWeightedDisplacement(displacement, &
   & structure,modes,qpoints) result(this)
  implicit none
  
  type(MassWeightedDisplacement), intent(in) :: displacement
  type(StructureData),            intent(in) :: structure
  type(RealMode),                 intent(in) :: modes(:)
  type(QpointData),               intent(in) :: qpoints(:)
  type(RealModeDisplacement)                 :: this
  
  type(QpointData), allocatable :: selected_qpoints(:)
  
  selected_qpoints = select_qpoints(modes, qpoints)
  
  this = RealModeDisplacement(RealSingleDisplacement( modes,           &
                                                    & displacement,    &
                                                    & structure,       &
                                                    & selected_qpoints ))
end function

! Converts a CartesianDisplacement to a RealModeDisplacement.
function new_RealModeDisplacement_CartesianDisplacement(displacement, &
   & structure,modes,qpoints) result(this)
  implicit none
  
  type(CartesianDisplacement), intent(in) :: displacement
  type(StructureData),         intent(in) :: structure
  type(RealMode),              intent(in) :: modes(:)
  type(QpointData),            intent(in) :: qpoints(:)
  type(RealModeDisplacement)              :: this
  
  type(QpointData), allocatable :: selected_qpoints(:)
  
  selected_qpoints = select_qpoints(modes, qpoints)
  
  this = RealModeDisplacement(RealSingleDisplacement( modes,           &
                                                    & displacement,    &
                                                    & structure,       &
                                                    & selected_qpoints ))
end function

! ----------------------------------------------------------------------
! The displacement along a given mode.
! ----------------------------------------------------------------------
impure elemental function displacement_id(this,id) result(output)
  implicit none
  
  class(RealModeDisplacement), intent(in) :: this
  integer,                     intent(in) :: id
  real(dp)                                :: output
  
  integer :: i
  
  i = first(this%vectors%id==id, default=0)
  if (i==0) then
    output = 0.0_dp
  else
    output = this%vectors(i)%magnitude
  endif
end function

impure elemental function displacement_mode(this,mode) result(output)
  implicit none
  
  class(RealModeDisplacement), intent(in) :: this
  type(RealMode),              intent(in) :: mode
  real(dp)                                :: output
  
  output = this%displacement(mode%id)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_RealModeDisplacement(this,input)
  implicit none
  
  class(RealModeDisplacement), intent(out) :: this
  type(String),                intent(in)  :: input(:)
  
  select type(this); type is(RealModeDisplacement)
    this = RealModeDisplacement(RealSingleDisplacement(input))
  class default
    call err()
  end select
end subroutine

function write_RealModeDisplacement(this) result(output)
  implicit none
  
  class(RealModeDisplacement), intent(in) :: this
  type(String), allocatable               :: output(:)
  
  select type(this); type is(RealModeDisplacement)
    output = str(this%vectors)
  class default
    call err()
  end select
end function

function new_RealModeDisplacement_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)   :: input(:)
  type(RealModeDisplacement) :: this
  
  call this%read(input)
end function

impure elemental function new_RealModeDisplacement_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(RealModeDisplacement)    :: this
  
  this = RealModeDisplacement(str(input))
end function
end module
