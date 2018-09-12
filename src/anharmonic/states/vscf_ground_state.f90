! ======================================================================
! The VSCF ground state of a given subspace,
!    stored in terms of coefficients of SubspaceBasis states.
! ======================================================================
module vscf_ground_state_module
  use common_module
  
  use subspace_basis_module
  use polynomial_state_module
  implicit none
  
  private
  
  public :: VscfGroundState
  public :: PolynomialState
  public :: initial_ground_state
  
  type, extends(Stringsable) :: VscfGroundState
    integer               :: subspace_id
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
  
  interface PolynomialState
    module procedure new_PolynomialState_VscfGroundState
  end interface
contains

! Constructor.
function new_VscfGroundState(subspace_id,wavevector,coefficients) &
   & result(this)
  implicit none
  
  integer,              intent(in) :: subspace_id
  type(FractionVector), intent(in) :: wavevector
  real(dp),             intent(in) :: coefficients(:)
  type(VscfGroundState)            :: this
  
  this%subspace_id  = subspace_id
  this%wavevector   = wavevector
  this%coefficients = coefficients
end function

! Construct a PolynomialState from a VscfGroundState.
impure elemental function new_PolynomialState_VscfGroundState(state,basis) &
   & result(output)
  implicit none
  
  type(VscfGroundState), intent(in) :: state
  type(SubspaceBasis),   intent(in) :: basis
  type(PolynomialState)             :: output
  
  real(dp), allocatable :: coefficients(:)
  
  integer :: i
  
  i = first(basis%wavevectors%wavevector==state%wavevector)
  
  coefficients = &
     & basis%wavevectors(i)%coefficients_basis_to_states(state%coefficients)
  
  output = PolynomialState( basis%subspace_id,           &
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
  
  type(SubspaceWavevectorBasis) :: wavevector_basis
  
  real(dp), allocatable :: coefficients(:)
  
  integer :: i
  
  ! Find the wavevector [0,0,0].
  wavevector_basis = basis%wavevectors(            &
     & first(is_int(basis%wavevectors%wavevector)) )
  
  ! Construct the coefficient vector in the basis of monomial states.
  ! All coefficients are zero, except for the coefficient of |0>, which is one.
  coefficients = [( 0.0_dp, i=1, size(wavevector_basis) )]
  coefficients(first(wavevector_basis%states%total_power()==0)) = 1
  
  ! Convert the coefficients into the orthonormal basis.
  coefficients = wavevector_basis%coefficients_states_to_basis(coefficients)
  
  ! Construct output.
  output = VscfGroundState(                        &
     & subspace_id  = basis%subspace_id,           &
     & wavevector   = wavevector_basis%wavevector, &
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
  type(FractionVector)  :: wavevector
  real(dp), allocatable :: coefficients(:)
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(VscfGroundState)
    line = split_line(input(1))
    subspace_id = int(line(2))
    
    line = split_line(input(2))
    wavevector = FractionVector(join(line(2:4)))
    
    line = split_line(input(3))
    coefficients = dble(line(2:))
    
    this = VscfGroundState(subspace_id,wavevector,coefficients)
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
