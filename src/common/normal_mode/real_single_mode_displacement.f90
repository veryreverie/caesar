! ======================================================================
! A displacement along a single real mode, in mass-weighted co-ordinates.
! ======================================================================
module caesar_real_single_mode_displacement_module
  use caesar_utils_module
  
  use caesar_structure_module
  
  use caesar_cartesian_displacement_module
  use caesar_mass_weighted_displacement_module
  use caesar_real_mode_module
  implicit none
  
  private
  
  public :: RealSingleDisplacement
  public :: MassWeightedDisplacement
  public :: CartesianDisplacement
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  public :: select_mode
  public :: select_modes
  
  type, extends(Stringable) :: RealSingleDisplacement
    ! The id of the mode.
    integer :: id
    
    ! The magnitude of vector along the mode, in mass-weighted co-ordinates.
    real(dp) :: magnitude
  contains
    ! I/O.
    procedure, public :: read  => read_RealSingleDisplacement
    procedure, public :: write => write_RealSingleDisplacement
  end type
  
  interface RealSingleDisplacement
    module procedure new_RealSingleDisplacement
    module procedure new_RealSingleDisplacement_RealMode
    module procedure new_RealSingleDisplacement_MassWeightedDisplacement
    module procedure new_RealSingleDisplacement_CartesianDisplacement
    module procedure new_RealSingleDisplacement_String
  end interface
  
  interface MassWeightedDisplacement
    module procedure new_MassWeightedDisplacement_RealSingleDisplacement
  end interface
  
  interface CartesianDisplacement
    module procedure new_CartesianDisplacement_RealSingleDisplacement
  end interface
  
  interface operator(*)
    module procedure multiply_real_RealSingleDisplacement
    module procedure multiply_RealSingleDisplacement_real
  end interface
  
  interface operator(/)
    module procedure divide_RealSingleDisplacement_real
  end interface
  
  interface operator(+)
    module procedure add_RealSingleDisplacement_RealSingleDisplacement
  end interface
  
  interface operator(-)
    module procedure negative_RealSingleDisplacement
    module procedure subtract_RealSingleDisplacement_RealSingleDisplacement
  end interface
  
  interface select_mode
    module procedure select_mode_RealSingleDisplacement
  end interface
  
  interface select_modes
    module procedure select_modes_RealSingleDisplacements
  end interface
contains

! Constructors.
impure elemental function new_RealSingleDisplacement(id,magnitude) result(this)
  implicit none
  
  integer,  intent(in)         :: id
  real(dp), intent(in)         :: magnitude
  type(RealSingleDisplacement) :: this
  
  this%id        = id
  this%magnitude = magnitude
end function

impure elemental function new_RealSingleDisplacement_RealMode(mode,magnitude) &
   & result(this)
  implicit none
  
  type(RealMode), intent(in)   :: mode
  real(dp),       intent(in)   :: magnitude
  type(RealSingleDisplacement) :: this
  
  this = RealSingleDisplacement(mode%id, magnitude)
end function

! ----------------------------------------------------------------------
! Arithmetic.
! ----------------------------------------------------------------------
impure elemental function multiply_real_RealSingleDisplacement(this,that) &
   & result(output)
  implicit none
  
  real(dp),                     intent(in) :: this
  type(RealSingleDisplacement), intent(in) :: that
  type(RealSingleDisplacement)             :: output
  
  output = RealSingleDisplacement( id        = that%id,            &
                                 & magnitude = this*that%magnitude )
end function

impure elemental function multiply_RealSingleDisplacement_real(this,that) &
   & result(output)
  implicit none
  
  type(RealSingleDisplacement), intent(in) :: this
  real(dp),                     intent(in) :: that
  type(RealSingleDisplacement)             :: output
  
  output = RealSingleDisplacement( id        = this%id,            &
                                 & magnitude = this%magnitude*that )
end function

impure elemental function divide_RealSingleDisplacement_real(this,that) &
   & result(output)
  implicit none
  
  type(RealSingleDisplacement), intent(in) :: this
  real(dp),                     intent(in) :: that
  type(RealSingleDisplacement)             :: output
  
  output = RealSingleDisplacement( id        = this%id,            &
                                 & magnitude = this%magnitude/that )
end function

impure elemental function add_RealSingleDisplacement_RealSingleDisplacement( &
   & this,that)  result(output)
  implicit none
  
  type(RealSingleDisplacement), intent(in) :: this
  type(RealSingleDisplacement), intent(in) :: that
  type(RealSingleDisplacement)             :: output
  
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to add vectors along different modes.')
    call err()
  endif
  
  output = RealSingleDisplacement( id        = this%id,                      &
                                 & magnitude = this%magnitude+that%magnitude )
end function

impure elemental function negative_RealSingleDisplacement(this) result(output)
  implicit none
  
  type(RealSingleDisplacement), intent(in) :: this
  type(RealSingleDisplacement)             :: output
  
  output = RealSingleDisplacement(id=this%id, magnitude=-this%magnitude)
end function

impure elemental function                                              &
   & subtract_RealSingleDisplacement_RealSingleDisplacement(this,that) &
   & result(output)
  implicit none
  
  type(RealSingleDisplacement), intent(in) :: this
  type(RealSingleDisplacement), intent(in) :: that
  type(RealSingleDisplacement)             :: output
  
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to subtract vectors along different &
       &modes.')
    call err()
  endif
  
  output = RealSingleDisplacement( id        = this%id,                      &
                                 & magnitude = this%magnitude-that%magnitude )
end function

! ----------------------------------------------------------------------
! Conversions to and from cartesian and mass-weighted co-ordinates.
! ----------------------------------------------------------------------
! Constructs the MassWeightedDisplacement corresponding to this vector.
impure elemental function                                                     &
   & new_MassWeightedDisplacement_RealSingleDisplacement(this,mode,structure, &
   & qpoint) result(output)
  implicit none
  
  class(RealSingleDisplacement), intent(in) :: this
  type(RealMode),                intent(in) :: mode
  type(StructureData),           intent(in) :: structure
  type(QpointData),              intent(in) :: qpoint
  type(MassWeightedDisplacement)            :: output
  
  if (mode%id/=this%id) then
    call print_line(CODE_ERROR//': Mode and vector incompatible.')
    call err()
  elseif ( qpoint%id/=mode%qpoint_id_plus .and. &
         & qpoint%id/=mode%qpoint_id_minus      ) then
    call print_line(CODE_ERROR//': Mode and q-point incompatible.')
    call err()
  endif
  
  output = this%magnitude * MassWeightedDisplacement(mode,structure,qpoint)
end function

! Constructs the CartesianDisplacement corresponding to this vector.
impure elemental function new_CartesianDisplacement_RealSingleDisplacement( &
   & this,mode,structure,qpoint) result(output)
  implicit none
  
  class(RealSingleDisplacement), intent(in) :: this
  type(RealMode),                intent(in) :: mode
  type(StructureData),           intent(in) :: structure
  type(QpointData),              intent(in) :: qpoint
  type(CartesianDisplacement)               :: output
  
  if (mode%id/=this%id) then
    call print_line(CODE_ERROR//': Mode and vector incompatible.')
    call err()
  elseif ( qpoint%id/=mode%qpoint_id_plus .and. &
         & qpoint%id/=mode%qpoint_id_minus      ) then
    call print_line(CODE_ERROR//': Mode and q-point incompatible.')
    call err()
  endif
  
  output = CartesianDisplacement(                            &
     & MassWeightedDisplacement(this,mode,structure,qpoint), &
     & structure                                             )
end function

! Constructs the vector corresponding to the component of a mass-weighted
!    vector along this mode.
impure elemental function                                             &
   & new_RealSingleDisplacement_MassWeightedDisplacement(mode,vector, &
   & structure,qpoint) result(this)
  implicit none
  
  type(RealMode),                 intent(in) :: mode
  type(MassWeightedDisplacement), intent(in) :: vector
  type(StructureData),            intent(in) :: structure
  type(QpointData),               intent(in) :: qpoint
  type(RealSingleDisplacement)               :: this
  
  type(MassWeightedDisplacement) :: mode_vector
  
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
  
  mode_vector = MassWeightedDisplacement(mode,structure,qpoint)
  
  magnitude = sum( vector%vectors       &
          &      * mode_vector%vectors) &
          & / structure%sc_size
  
  this = RealSingleDisplacement(id=mode%id, magnitude=magnitude)
end function

! Constructs the vector corresponding to the component of a cartesian
!    vector along this mode.
impure elemental function new_RealSingleDisplacement_CartesianDisplacement( &
   & mode,vector,structure,qpoint) result(this)
  implicit none
  
  type(RealMode),              intent(in) :: mode
  type(CartesianDisplacement), intent(in) :: vector
  type(StructureData),         intent(in) :: structure
  type(QpointData),            intent(in) :: qpoint
  type(RealSingleDisplacement)            :: this
  
  this = RealSingleDisplacement( mode,                                       &
                               & MassWeightedDisplacement(vector,structure), &
                               & structure,                                  &
                               & qpoint                                      )
end function

! ----------------------------------------------------------------------
! Select modes corresponding to a given displacement or displacements.
! ----------------------------------------------------------------------
function select_mode_RealSingleDisplacement(displacement,modes) result(output)
  implicit none
  
  type(RealSingleDisplacement), intent(in) :: displacement
  type(RealMode),        intent(in) :: modes(:)
  type(RealMode)                    :: output
  
  output = modes(first(modes%id==displacement%id))
end function

function select_modes_RealSingleDisplacements(displacements,modes) &
   & result(output)
  implicit none
  
  type(RealSingleDisplacement), intent(in) :: displacements(:)
  type(RealMode),               intent(in) :: modes(:)
  type(RealMode), allocatable              :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(displacements)), stat=ialloc); call err(ialloc)
  do i=1,size(displacements)
    output(i) = select_mode(displacements(i), modes)
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_RealSingleDisplacement(this,input)
  implicit none
  
  class(RealSingleDisplacement), intent(out) :: this
  type(String),                  intent(in)  :: input
  
  type(String), allocatable :: split_string(:)
  integer                   :: id
  real(dp)                  :: magnitude
  
  select type(this); type is(RealSingleDisplacement)
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
    
    this = RealSingleDisplacement(id,magnitude)
  end select
end subroutine

function write_RealSingleDisplacement(this) result(output)
  implicit none
  
  class(RealSingleDisplacement), intent(in) :: this
  type(String)                              :: output
  
  select type(this); type is(RealSingleDisplacement)
    output = 'u'//this%id//' = '//this%magnitude
  end select
end function

impure elemental function new_RealSingleDisplacement_String(input) &
   & result(this)
  implicit none
  
  type(String), intent(in)     :: input
  type(RealSingleDisplacement) :: this
  
  call this%read(input)
end function
end module
