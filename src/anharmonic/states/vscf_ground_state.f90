! ======================================================================
! The VSCF ground state of a given subspace,
!    stored in terms of coefficients of SubspaceBasis states.
! ======================================================================
module vscf_ground_state_module
  use common_module
  
  use subspace_basis_module
  use sum_state_module
  implicit none
  
  private
  
  public :: VscfGroundState
  public :: SumState
  public :: initial_ground_state
  
  type, extends(Stringsable) :: VscfGroundState
    integer               :: subspace_id
    real(dp)              :: energy
    type(FractionVector)  :: wavevector
    real(dp), allocatable :: coefficients(:)
  contains
    procedure, public :: read  => read_VscfGroundState
    procedure, public :: write => write_VscfGroundState
  end type
  
  interface VscfGroundState
    module procedure new_VscfGroundState
    module procedure new_VscfGroundState_Strings
    module procedure new_VscfGroundState_StringArray
  end interface
  
  interface SumState
    module procedure new_SumState_VscfGroundState
  end interface
contains

! Constructor.
function new_VscfGroundState(subspace_id,energy,wavevector,coefficients) &
   & result(this)
  implicit none
  
  integer,              intent(in) :: subspace_id
  real(dp),             intent(in) :: energy
  type(FractionVector), intent(in) :: wavevector
  real(dp),             intent(in) :: coefficients(:)
  type(VscfGroundState)            :: this
  
  this%subspace_id = subspace_id
  this%energy = energy
  this%wavevector = wavevector
  this%coefficients = coefficients
end function

! Construct a SumState from a VscfGroundState.
function new_SumState_VscfGroundState(state,basis) result(output)
  implicit none
  
  type(VscfGroundState), intent(in) :: state
  type(SubspaceBasis),   intent(in) :: basis
  type(SumState)                    :: output
  
  real(dp), allocatable :: coefficients(:)
  
  integer :: i
  
  i = first(basis%wavevectors%wavevector==state%wavevector)
  
  coefficients = dble( basis%wavevectors(i)%basis_to_states &
                   & * vec(state%coefficients)              )
  
  output = SumState( basis%subspace_id,           &
                   & basis%wavevectors(i)%states, &
                   & coefficients                 )
end function

! Generate initial guess. This is simply the basis state |0>, i.e. the
!    ground state of the effective harmonic potential from which the basis
!    states were generated.
impure elemental function initial_ground_state(basis) result(output)
  implicit none
  
  type(SubspaceBasis), intent(in) :: basis
  type(VscfGroundState)           :: output
  
  integer,  allocatable :: powers(:)
  real(dp), allocatable :: basis_to_states(:,:)
  real(dp), allocatable :: coefficients(:)
  
  integer :: i,j
  
  ! Find the wavevector [0,0,0].
  i = first(is_int(basis%wavevectors%wavevector))
  
  ! Find the state |0> at that wavevector.
  j = first(basis%wavevectors(i)%states%total_power()==0)
  
  ! Extract the coefficients of the basis functions which refer to |0>.
  basis_to_states = dble(basis%wavevectors(i)%basis_to_states)
  coefficients = basis_to_states(j,:)
  
  ! Construct output. The initial energy is just a placeholder.
  output = VscfGroundState(                            &
     & subspace_id  = basis%subspace_id,               &
     & energy       = 0.0_dp,                          &
     & wavevector   = basis%wavevectors(i)%wavevector, &
     & coefficients = coefficients)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_VscfGroundState(this,input)
  implicit none
  
  class(VscfGroundState), intent(out) :: this
  type(String),           intent(in)  :: input(:)
  
  integer               :: subspace_id
  real(dp)              :: energy
  type(FractionVector)  :: wavevector
  real(dp), allocatable :: coefficients(:)
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(VscfGroundState)
    line = split_line(input(1))
    subspace_id = int(line(2))
    
    line = split_line(input(2))
    energy = dble(line(2))
    
    line = split_line(input(3))
    wavevector = FractionVector(join(line(2:4)))
    
    line = split_line(input(4))
    coefficients = dble(line(:))
    
    this = VscfGroundState(subspace_id,energy,wavevector,coefficients)
  class default
    call err()
  end select
end subroutine

function write_VscfGroundState(this) result(output)
  implicit none
  
  class(VscfGroundState), intent(in) :: this
  type(String), allocatable          :: output(:)
  
  select type(this); type is(VscfGroundState)
    output = [ 'Subspace '//this%subspace_id,     &
             & 'Energy   '//this%energy,          &
             & 'Wavevector '//this%wavevector,    &
             & 'Coefficients '//this%coefficients ]
  class default
    call err()
  end select
end function

function new_VscfGroundState_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(VscfGroundState)    :: this
  
  call this%read(input)
end function

impure elemental function new_VscfGroundState_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(VscfGroundState)         :: this
  
  this = VscfGroundState(str(input))
end function
end module
