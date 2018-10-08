! ======================================================================
! A force in real mode co-ordinates.
! ======================================================================
module real_mode_force_module
  use utils_module
  
  use structure_module
  
  use cartesian_force_module
  use mass_weighted_force_module
  use real_mode_module
  use real_single_mode_force_module
  implicit none
  
  private
  
  public :: RealModeForce
  public :: size
  public :: MassWeightedForce
  public :: CartesianForce
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
  type, extends(Stringsable) :: RealModeForce
    type(RealSingleForce), allocatable :: vectors(:)
  contains
    ! The component of the force along a given mode.
    generic,   public  :: force =>  &
                        & force_id, &
                        & force_mode
    procedure, private :: force_id
    procedure, private :: force_mode
    
    ! I/O.
    procedure, public :: read  => read_RealModeForce
    procedure, public :: write => write_RealModeForce
  end type
  
  interface RealModeForce
    module procedure new_RealModeForce
    module procedure new_RealModeForce_RealModes
    module procedure new_RealModeForce_MassWeightedForce
    module procedure new_RealModeForce_CartesianForce
    module procedure new_RealModeForce_Strings
    module procedure new_RealModeForce_StringArray
  end interface
  
  interface size
    module procedure size_RealModeForce
  end interface
  
  interface MassWeightedForce
    module procedure new_MassWeightedForce_RealModeForce
  end interface
  
  interface CartesianForce
    module procedure new_CartesianForce_RealModeForce
  end interface
  
  interface operator(*)
    module procedure multiply_real_RealModeForce
    module procedure multiply_RealModeForce_real
  end interface
  
  interface operator(/)
    module procedure divide_RealModeForce_real
  end interface
  
  interface operator(+)
    module procedure add_RealModeForce_RealModeForce
  end interface
  
  interface sum
    module procedure sum_RealModeForces
  end interface
  
  interface operator(-)
    module procedure negative_RealModeForce
    module procedure subtract_RealModeForce_RealModeForce
  end interface
contains

! Constructors and size() function.
function new_RealModeForce(forces) result(this)
  implicit none
  
  type(RealSingleForce), intent(in) :: forces(:)
  type(RealModeForce)               :: this
  
  this%vectors = forces
end function

function new_RealModeForce_RealModes(modes,forces) result(this)
  implicit none
  
  type(RealMode), intent(in) :: modes(:)
  real(dp),       intent(in) :: forces(:)
  type(RealModeForce)        :: this
  
  this = RealModeForce(RealSingleForce(modes,forces))
end function

function size_RealModeForce(this) result(output)
  implicit none
  
  type(RealModeForce), intent(in) :: this
  integer                         :: output
  
  output = size(this%vectors)
end function

! ----------------------------------------------------------------------
! Arithmetic.
! ----------------------------------------------------------------------
impure elemental function multiply_real_RealModeForce(this,that) &
   & result(output)
  implicit none
  
  real(dp),            intent(in) :: this
  type(RealModeForce), intent(in) :: that
  type(RealModeForce)             :: output
  
  output = RealModeForce(this*that%vectors)
end function

impure elemental function multiply_RealModeForce_real(this,that) &
   & result(output)
  implicit none
  
  type(RealModeForce), intent(in) :: this
  real(dp),            intent(in) :: that
  type(RealModeForce)             :: output
  
  output = RealModeForce(this%vectors*that)
end function

impure elemental function divide_RealModeForce_real(this,that) &
   & result(output)
  implicit none
  
  type(RealModeForce), intent(in) :: this
  real(dp),            intent(in) :: that
  type(RealModeForce)             :: output
  
  output = RealModeForce(this%vectors/that)
end function

impure elemental function add_RealModeForce_RealModeForce(this, &
   & that) result(output)
  implicit none
  
  type(RealModeForce), intent(in) :: this
  type(RealModeForce), intent(in) :: that
  type(RealModeForce)             :: output
  
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

function sum_RealModeForces(this) result(output)
  implicit none
  
  type(RealModeForce), intent(in) :: this(:)
  type(RealModeForce)             :: output
  
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

impure elemental function negative_RealModeForce(this) result(output)
  implicit none
  
  type(RealModeForce), intent(in) :: this
  type(RealModeForce)             :: output
  
  output = RealModeForce(-this%vectors)
end function

impure elemental function subtract_RealModeForce_RealModeForce( &
   & this,that) result(output)
  implicit none
  
  type(RealModeForce), intent(in) :: this
  type(RealModeForce), intent(in) :: that
  type(RealModeForce)             :: output
  
  output = this + (-that)
end function

! ----------------------------------------------------------------------
! Conversions to and from cartesian and mass-weighted co-ordinates.
! ----------------------------------------------------------------------
! Returns the force in mass-weighted co-ordinates.
function new_MassWeightedForce_RealModeForce(this,structure, &
   & modes,qpoints) result(output)
  implicit none
  
  class(RealModeForce), intent(in) :: this
  type(StructureData),  intent(in) :: structure
  type(RealMode),       intent(in) :: modes(:)
  type(QpointData),     intent(in) :: qpoints(:)
  type(MassWeightedForce)          :: output
  
  type(RealMode),   allocatable :: selected_modes(:)
  type(QpointData), allocatable :: selected_qpoints(:)
  
  selected_modes = select_modes(this%vectors, modes)
  selected_qpoints = select_qpoints(selected_modes, qpoints)
  output = sum(MassWeightedForce( this%vectors,    &
                                & selected_modes,  &
                                & structure,       &
                                & selected_qpoints ))
end function

! Returns the force in cartesian co-ordinates.
function new_CartesianForce_RealModeForce(this,structure, &
   & modes,qpoints) result(output)
  implicit none
  
  class(RealModeForce), intent(in) :: this
  type(StructureData),  intent(in) :: structure
  type(RealMode),       intent(in) :: modes(:)
  type(QpointData),     intent(in) :: qpoints(:)
  type(CartesianForce)             :: output
  
  type(RealMode),   allocatable :: selected_modes(:)
  type(QpointData), allocatable :: selected_qpoints(:)
  
  selected_modes = select_modes(this%vectors, modes)
  selected_qpoints = select_qpoints(selected_modes, qpoints)
  output = sum(CartesianForce( this%vectors,    &
                             & selected_modes,  &
                             & structure,       &
                             & selected_qpoints ))
end function

! Converts a MassWeightedForce to a RealModeForce.
function new_RealModeForce_MassWeightedForce(force, &
   & structure,modes,qpoints) result(this)
  implicit none
  
  type(MassWeightedForce), intent(in) :: force
  type(StructureData),     intent(in) :: structure
  type(RealMode),          intent(in) :: modes(:)
  type(QpointData),        intent(in) :: qpoints(:)
  type(RealModeForce)                 :: this
  
  type(QpointData), allocatable :: selected_qpoints(:)
  
  selected_qpoints = select_qpoints(modes, qpoints)
  
  this = RealModeForce(RealSingleForce( modes,           &
                                      & force,           &
                                      & structure,       &
                                      & selected_qpoints ))
end function

! Converts a CartesianForce to a RealModeForce.
function new_RealModeForce_CartesianForce(force, &
   & structure,modes,qpoints) result(this)
  implicit none
  
  type(CartesianForce), intent(in) :: force
  type(StructureData),  intent(in) :: structure
  type(RealMode),       intent(in) :: modes(:)
  type(QpointData),     intent(in) :: qpoints(:)
  type(RealModeForce)              :: this
  
  type(QpointData), allocatable :: selected_qpoints(:)
  
  selected_qpoints = select_qpoints(modes, qpoints)
  
  this = RealModeForce(RealSingleForce( modes,           &
                                      & force,           &
                                      & structure,       &
                                      & selected_qpoints ))
end function

! ----------------------------------------------------------------------
! The displacement along a given mode.
! ----------------------------------------------------------------------
impure elemental function force_id(this,id) result(output)
  implicit none
  
  class(RealModeForce), intent(in) :: this
  integer,              intent(in) :: id
  real(dp)                         :: output
  
  integer :: i
  
  i = first(this%vectors%id==id, default=0)
  if (i==0) then
    output = 0.0_dp
  else
    output = this%vectors(i)%magnitude
  endif
end function

impure elemental function force_mode(this,mode) result(output)
  implicit none
  
  class(RealModeForce), intent(in) :: this
  type(RealMode),       intent(in) :: mode
  real(dp)                         :: output
  
  output = this%force(mode%id)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_RealModeForce(this,input)
  implicit none
  
  class(RealModeForce), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  select type(this); type is(RealModeForce)
    this = RealModeForce(RealSingleForce(input))
  class default
    call err()
  end select
end subroutine

function write_RealModeForce(this) result(output)
  implicit none
  
  class(RealModeForce), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  select type(this); type is(RealModeForce)
    output = str(this%vectors)
  class default
    call err()
  end select
end function

function new_RealModeForce_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(RealModeForce)      :: this
  
  call this%read(input)
end function

impure elemental function new_RealModeForce_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(RealModeForce)           :: this
  
  this = RealModeForce(str(input))
end function
end module
