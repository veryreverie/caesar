! ======================================================================
! A force along a single real mode, in mass-weighted co-ordinates.
! ======================================================================
module real_single_mode_force_module
  use utils_module
  
  use structure_module
  
  use cartesian_force_module
  use mass_weighted_force_module
  use real_mode_module
  implicit none
  
  private
  
  public :: RealSingleForce
  public :: MassWeightedForce
  public :: CartesianForce
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  public :: select_mode
  public :: select_modes
  
  type, extends(Stringable) :: RealSingleForce
    ! The id of the mode.
    integer :: id
    
    ! The magnitude of vector along the mode, in mass-weighted co-ordinates.
    real(dp) :: magnitude
  contains
    ! I/O.
    procedure, public :: read  => read_RealSingleForce
    procedure, public :: write => write_RealSingleForce
  end type
  
  interface RealSingleForce
    module procedure new_RealSingleForce
    module procedure new_RealSingleForce_RealMode
    module procedure new_RealSingleForce_MassWeightedForce
    module procedure new_RealSingleForce_CartesianForce
    module procedure new_RealSingleForce_String
  end interface
  
  interface MassWeightedForce
    module procedure new_MassWeightedForce_RealSingleForce
  end interface
  
  interface CartesianForce
    module procedure new_CartesianForce_RealSingleForce
  end interface
  
  interface operator(*)
    module procedure multiply_real_RealSingleForce
    module procedure multiply_RealSingleForce_real
  end interface
  
  interface operator(/)
    module procedure divide_RealSingleForce_real
  end interface
  
  interface operator(+)
    module procedure add_RealSingleForce_RealSingleForce
  end interface
  
  interface operator(-)
    module procedure negative_RealSingleForce
    module procedure subtract_RealSingleForce_RealSingleForce
  end interface
  
  interface select_mode
    module procedure select_mode_RealSingleForce
  end interface
  
  interface select_modes
    module procedure select_modes_RealSingleForces
  end interface
contains

! Constructors.
impure elemental function new_RealSingleForce(id,magnitude) result(this)
  implicit none
  
  integer,  intent(in)  :: id
  real(dp), intent(in)  :: magnitude
  type(RealSingleForce) :: this
  
  this%id        = id
  this%magnitude = magnitude
end function

impure elemental function new_RealSingleForce_RealMode(mode,magnitude) &
   & result(this)
  implicit none
  
  type(RealMode), intent(in) :: mode
  real(dp),       intent(in) :: magnitude
  type(RealSingleForce)      :: this
  
  this = RealSingleForce(mode%id, magnitude)
end function

! ----------------------------------------------------------------------
! Arithmetic.
! ----------------------------------------------------------------------
impure elemental function multiply_real_RealSingleForce(this,that) &
   & result(output)
  implicit none
  
  real(dp),              intent(in) :: this
  type(RealSingleForce), intent(in) :: that
  type(RealSingleForce)             :: output
  
  output = RealSingleForce( id        = that%id,            &
                          & magnitude = this*that%magnitude )
end function

impure elemental function multiply_RealSingleForce_real(this,that) &
   & result(output)
  implicit none
  
  type(RealSingleForce), intent(in) :: this
  real(dp),              intent(in) :: that
  type(RealSingleForce)             :: output
  
  output = RealSingleForce( id        = this%id,            &
                          & magnitude = this%magnitude*that )
end function

impure elemental function divide_RealSingleForce_real(this,that) &
   & result(output)
  implicit none
  
  type(RealSingleForce), intent(in) :: this
  real(dp),              intent(in) :: that
  type(RealSingleForce)             :: output
  
  output = RealSingleForce( id        = this%id,            &
                          & magnitude = this%magnitude/that )
end function

impure elemental function add_RealSingleForce_RealSingleForce(this,that) &
   & result(output)
  implicit none
  
  type(RealSingleForce), intent(in) :: this
  type(RealSingleForce), intent(in) :: that
  type(RealSingleForce)             :: output
  
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to add vectors along different modes.')
    call err()
  endif
  
  output = RealSingleForce( id        = this%id,                      &
                          & magnitude = this%magnitude+that%magnitude )
end function

impure elemental function negative_RealSingleForce(this) result(output)
  implicit none
  
  type(RealSingleForce), intent(in) :: this
  type(RealSingleForce)             :: output
  
  output = RealSingleForce(id=this%id, magnitude=-this%magnitude)
end function

impure elemental function subtract_RealSingleForce_RealSingleForce(this,that) &
   & result(output)
  implicit none
  
  type(RealSingleForce), intent(in) :: this
  type(RealSingleForce), intent(in) :: that
  type(RealSingleForce)             :: output
  
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to subtract vectors along different &
       &modes.')
    call err()
  endif
  
  output = RealSingleForce( id        = this%id,                      &
                          & magnitude = this%magnitude-that%magnitude )
end function

! ----------------------------------------------------------------------
! Conversions to and from cartesian and mass-weighted co-ordinates.
! ----------------------------------------------------------------------
! Constructs the MassWeightedForce corresponding to this vector.
impure elemental function new_MassWeightedForce_RealSingleForce(this,mode, &
   & structure,qpoint) result(output)
  implicit none
  
  class(RealSingleForce), intent(in) :: this
  type(RealMode),         intent(in) :: mode
  type(StructureData),    intent(in) :: structure
  type(QpointData),       intent(in) :: qpoint
  type(MassWeightedForce)            :: output
  
  if (mode%id/=this%id) then
    call print_line(CODE_ERROR//': Mode and vector incompatible.')
    call err()
  elseif ( qpoint%id/=mode%qpoint_id_plus .and. &
         & qpoint%id/=mode%qpoint_id_minus      ) then
    call print_line(CODE_ERROR//': Mode and q-point incompatible.')
    call err()
  endif
  
  output = this%magnitude * MassWeightedForce(mode,structure,qpoint)
end function

! Constructs the CartesianForce corresponding to this vector.
impure elemental function new_CartesianForce_RealSingleForce(this,mode, &
   & structure,qpoint) result(output)
  implicit none
  
  class(RealSingleForce), intent(in) :: this
  type(RealMode),         intent(in) :: mode
  type(StructureData),    intent(in) :: structure
  type(QpointData),       intent(in) :: qpoint
  type(CartesianForce)               :: output
  
  if (mode%id/=this%id) then
    call print_line(CODE_ERROR//': Mode and vector incompatible.')
    call err()
  elseif ( qpoint%id/=mode%qpoint_id_plus .and. &
         & qpoint%id/=mode%qpoint_id_minus      ) then
    call print_line(CODE_ERROR//': Mode and q-point incompatible.')
    call err()
  endif
  
  output = CartesianForce( MassWeightedForce(this,mode,structure,qpoint), &
                         & structure                                      )
end function

! Constructs the vector corresponding to the component of a mass-weighted
!    vector along this mode.
impure elemental function new_RealSingleForce_MassWeightedForce(mode,vector, &
   & structure,qpoint) result(this)
  implicit none
  
  type(RealMode),          intent(in) :: mode
  type(MassWeightedForce), intent(in) :: vector
  type(StructureData),     intent(in) :: structure
  type(QpointData),        intent(in) :: qpoint
  type(RealSingleForce)               :: this
  
  type(MassWeightedForce) :: mode_vector
  
  real(dp) :: magnitude
  
  ! Modes are orthonormal in mass-reduced co-ordinates.
  ! If M is the mass-weighting matrix, M_ab = 1/sqrt(m_a*m_b),
  !    where m_a and m_b are the masses of atoms a and b respectively, then
  !
  ! u_i and u_j are the cartesian representations of modes i and j.
  ! r is the cartesian vector.
  !
  ! u_i.M.u_j = 0 if i/=j
  !           = n if i=j, where n is the number of primitive cells.
  !
  ! r = sum_j[ a_j*u_j ]
  ! => u_i.M.r = sum_j[ a_j*u_i.M.u_j ] = a_i*u_i.M.u_i = a_i*n
  ! => a_j = u_i.M.r / n
  
  mode_vector = MassWeightedForce(mode,structure,qpoint)
  
  magnitude = sum( vector%vectors       &
          &      * mode_vector%vectors) &
          & / structure%sc_size
  
  this = RealSingleForce(id=mode%id, magnitude=magnitude)
end function

! Constructs the vector corresponding to the component of a cartesian
!    vector along this mode.
impure elemental function new_RealSingleForce_CartesianForce(mode,vector, &
   & structure,qpoint) result(this)
  implicit none
  
  type(RealMode),       intent(in) :: mode
  type(CartesianForce), intent(in) :: vector
  type(StructureData),  intent(in) :: structure
  type(QpointData),     intent(in) :: qpoint
  type(RealSingleForce)            :: this
  
  this = RealSingleForce( mode,                                &
                        & MassWeightedForce(vector,structure), &
                        & structure,                           &
                        & qpoint                               )
end function

! ----------------------------------------------------------------------
! Select modes corresponding to a given force or forces.
! ----------------------------------------------------------------------
function select_mode_RealSingleForce(force,modes) result(output)
  implicit none
  
  type(RealSingleForce), intent(in) :: force
  type(RealMode),        intent(in) :: modes(:)
  type(RealMode)                    :: output
  
  output = modes(first(modes%id==force%id))
end function

function select_modes_RealSingleForces(forces,modes) result(output)
  implicit none
  
  type(RealSingleForce), intent(in) :: forces(:)
  type(RealMode),        intent(in) :: modes(:)
  type(RealMode), allocatable       :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(forces)), stat=ialloc); call err(ialloc)
  do i=1,size(forces)
    output(i) = select_mode(forces(i), modes)
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_RealSingleForce(this,input)
  implicit none
  
  class(RealSingleForce), intent(out) :: this
  type(String),           intent(in)  :: input
  
  type(String), allocatable :: split_string(:)
  integer                   :: id
  real(dp)                  :: magnitude
  
  select type(this); type is(RealSingleForce)
    split_string = split_line(input)
    if (size(split_string)/=3) then
      call print_line(ERROR//': unable to parse real single mode vector &
         &from string: '//input)
      call err()
    endif
    
    ! If e.g. id=3 and power=2.1 then split_string = ["u3","=","2.1"]
    ! The 'u' needs stripping off the first element to give the id.
    id = int(slice(split_string(1),2,len(split_string(1))))
    magnitude = dble(split_string(3))
    
    this = RealSingleForce(id,magnitude)
  end select
end subroutine

function write_RealSingleForce(this) result(output)
  implicit none
  
  class(RealSingleForce), intent(in) :: this
  type(String)                       :: output
  
  select type(this); type is(RealSingleForce)
    output = 'u'//this%id//' = '//this%magnitude
  end select
end function

impure elemental function new_RealSingleForce_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(RealSingleForce)    :: this
  
  call this%read(input)
end function
end module
