! ======================================================================
! A displacement along a single real mode.
! ======================================================================
module real_single_mode_displacement_submodule
  use utils_module
  
  use structure_module
  
  use cartesian_displacement_submodule
  use real_mode_submodule
  implicit none
  
  private
  
  public :: RealSingleModeDisplacement
  
  ! The displacement along a single complex mode.
  type, extends(Stringable) :: RealSingleModeDisplacement
    ! The id of the mode.
    integer :: id
    
    ! The displacement along the mode.
    real(dp) :: displacement
  contains
    ! Convert to cartesian co-ordinates.
    procedure, public :: cartesian_displacement => &
       & cartesian_displacement_RealSingleModeDisplacement
    
    ! I/O.
    procedure, public :: read  => read_RealSingleModeDisplacement
    procedure, public :: write => write_RealSingleModeDisplacement
  end type
  
  interface RealSingleModeDisplacement
    module procedure new_RealSingleModeDisplacement
    module procedure new_RealSingleModeDisplacement_CartesianDisplacement
    module procedure new_RealSingleModeDisplacement_String
  end interface
contains

! Constructor.
function new_RealSingleModeDisplacement(id,displacement) result(this)
  implicit none
  
  integer,  intent(in)             :: id
  real(dp), intent(in)             :: displacement
  type(RealSingleModeDisplacement) :: this
  
  this%id           = id
  this%displacement = displacement
end function

! ----------------------------------------------------------------------
! Conversions to and from cartesian co-ordinates.
! ----------------------------------------------------------------------
! Constructs the CartesianDisplacement corresponding to this displacement.
function cartesian_displacement_RealSingleModeDisplacement(this,real_mode, &
   & structure,qpoint) result(output)
  implicit none
  
  class(RealSingleModeDisplacement), intent(in) :: this
  type(RealMode),                    intent(in) :: real_mode
  type(StructureData),               intent(in) :: structure
  type(QpointData),                  intent(in) :: qpoint
  type(CartesianDisplacement)                   :: output
  
  if (real_mode%id/=this%id) then
    call print_line(CODE_ERROR//': Mode and displacement incompatible.')
    call err()
  endif
  
  output = this%displacement &
       & * real_mode%cartesian_displacement(structure,qpoint)
end function

! Constructs the displacement corresponding to the component of a cartesian
!    displacement along this mode.
function new_RealSingleModeDisplacement_CartesianDisplacement(mode, &
   & displacement,structure,qpoint) result(this)
  implicit none
  
  type(RealMode),              intent(in) :: mode
  type(CartesianDisplacement), intent(in) :: displacement
  type(StructureData),         intent(in) :: structure
  type(QpointData),            intent(in) :: qpoint
  type(RealSingleModeDisplacement)        :: this
  
  type(CartesianDisplacement) :: mode_displacement
  
  real(dp) :: numerator
  real(dp) :: denominator
  
  integer :: i
  
  ! Modes are orthogonal in mass-reduced co-ordinates,
  !    but normalised in cartesian co-ordinates.
  ! If M is the mass-weighting matrix, M_ab = 1/sqrt(m_a*m_b),
  !    where m_a and m_b are the masses of atoms a and b respectively, then
  !
  ! u_i and u_j are the cartesian representations of modes i and j.
  ! r is the cartesian displacement.
  !
  ! u_i.u_i=1   for all i (Modes are normalised in cartesian co-ordinates.)
  ! u_i.M.u_j=0 for i/=j  (Modes are orthogonal in mass-reduced co-ordinates.
  !
  ! r = sum_j a_j*u_j
  ! => u_i.M.r = sum_j a_j*u_i.M.u_j = a_i*u_i.M.u_i
  ! => a_j = u_i.M.r / u_i.M.u_i
  
  mode_displacement = mode%cartesian_displacement(structure,qpoint)
  
  numerator = 0
  denominator = 0
  do i=1,structure%no_atoms
    numerator = numerator                          &
            & + displacement%displacements(i)      &
            & * mode_displacement%displacements(i) &
            & / structure%atoms(i)%mass()
    numerator = numerator &
            & + mode_displacement%displacements(i) &
            & * mode_displacement%displacements(i) &
            & / structure%atoms(i)%mass()
  enddo
  
  this = RealSingleModeDisplacement( id=mode%id, &
                                   & displacement=numerator/denominator)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_RealSingleModeDisplacement(this,input)
  implicit none
  
  class(RealSingleModeDisplacement), intent(out) :: this
  type(String),                      intent(in)  :: input
  
  type(String), allocatable :: split_string(:)
  integer                   :: id
  real(dp)                  :: displacement
  
  select type(this); type is(RealSingleModeDisplacement)
    split_string = split_line(input)
    if (size(split_string)/=3) then
      call print_line(ERROR//': unable to parse real single mode displacement &
         &from string: '//input)
      call err()
    endif
    
    ! If e.g. id=3 and power=2.1 then split_string = ["u3","=","2.1"]
    ! The 'u' needs stripping off the first element to give the id.
    id = int(slice(split_string(1),2,len(split_string(1))))
    displacement = dble(split_string(3))
    
    this = RealSingleModeDisplacement(id,displacement)
  end select
end subroutine

function write_RealSingleModeDisplacement(this) result(output)
  implicit none
  
  class(RealSingleModeDisplacement), intent(in) :: this
  type(String)                                  :: output
  
  select type(this); type is(RealSingleModeDisplacement)
    output = 'u'//this%id//' = '//this%displacement
  end select
end function

impure elemental function new_RealSingleModeDisplacement_String(input) &
   & result(this)
  implicit none
  
  type(String), intent(in)         :: input
  type(RealSingleModeDisplacement) :: this
  
  this = input
end function
end module
