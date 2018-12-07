! ======================================================================
! A VSCF state in a given subspace,
!    stored in terms of coefficients of harmonic SubspaceBasis states.
! ======================================================================
module vscf_state_module
  use common_module
  
  use wavevector_basis_module
  use subspace_basis_module
  use polynomial_state_module
  implicit none
  
  private
  
  public :: VscfState
  public :: PolynomialState
  public :: initial_ground_state
  
  type, extends(Stringsable) :: VscfState
    integer               :: subspace_id
    type(FractionVector)  :: wavevector
    integer               :: degeneracy
    real(dp)              :: energy
    real(dp), allocatable :: coefficients(:)
  contains
    procedure, public :: read  => read_VscfState
    procedure, public :: write => write_VscfState
  end type
  
  interface VscfState
    module procedure new_VscfState
    module procedure new_VscfState_Strings
    module procedure new_VscfState_StringArray
  end interface
  
  interface PolynomialState
    module procedure new_PolynomialState_VscfState
  end interface
contains

! Constructor.
function new_VscfState(subspace_id,wavevector,degeneracy,energy,coefficients) &
   & result(this)
  implicit none
  
  integer,              intent(in) :: subspace_id
  type(FractionVector), intent(in) :: wavevector
  integer,              intent(in) :: degeneracy
  real(dp),             intent(in) :: energy
  real(dp),             intent(in) :: coefficients(:)
  type(VscfState)                  :: this
  
  this%subspace_id  = subspace_id
  this%wavevector   = wavevector
  this%degeneracy   = degeneracy
  this%energy       = energy
  this%coefficients = coefficients
end function

! Construct a PolynomialState from a VscfState.
impure elemental function new_PolynomialState_VscfState(state,basis) &
   & result(output)
  implicit none
  
  type(VscfState),     intent(in) :: state
  type(SubspaceBasis), intent(in) :: basis
  type(PolynomialState)           :: output
  
  real(dp), allocatable :: coefficients(:)
  
  integer :: i
  
  i = first(basis%wavevectors%wavevector==state%wavevector)
  
  coefficients = &
     & basis%wavevectors(i)%coefficients_basis_to_states(state%coefficients)
  
  output = PolynomialState( basis%subspace_id,                    &
                          & basis%wavevectors(i)%monomial_states, &
                          & coefficients                          )
end function

! Generate initial guess. This is simply the basis state |0>, i.e. the
!    ground state of the effective harmonic potential from which the basis
!    states were generated.
impure elemental function initial_ground_state(basis) result(output)
  implicit none
  
  type(SubspaceBasis), intent(in) :: basis
  type(VscfState)                 :: output
  
  type(WavevectorBasis) :: wavevector_basis
  
  real(dp), allocatable :: coefficients(:)
  
  integer :: i
  
  ! Find the wavevector [0,0,0].
  wavevector_basis = basis%wavevectors(            &
     & first(is_int(basis%wavevectors%wavevector)) )
  
  ! Construct the coefficient vector in the basis of monomial states.
  ! All coefficients are zero, except for the coefficient of |0>, which is one.
  coefficients = [( 0.0_dp, i=1, size(wavevector_basis) )]
  coefficients(                                                      &
     & first(wavevector_basis%harmonic_states%total_occupation()==0) ) = 1
  
  ! Convert the coefficients into the orthonormal basis.
  coefficients = wavevector_basis%coefficients_states_to_basis(coefficients)
  
  ! Construct output.
  output = VscfState(                              &
     & subspace_id  = basis%subspace_id,           &
     & wavevector   = wavevector_basis%wavevector, &
     & degeneracy   = wavevector_basis%degeneracy, &
     & energy       = 0.0_dp,                      &
     & coefficients = coefficients                 )
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_VscfState(this,input)
  implicit none
  
  class(VscfState), intent(out) :: this
  type(String),     intent(in)  :: input(:)
  
  integer               :: subspace_id
  type(FractionVector)  :: wavevector
  integer               :: degeneracy
  real(dp)              :: energy
  real(dp), allocatable :: coefficients(:)
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(VscfState)
    line = split_line(input(1))
    subspace_id = int(line(3))
    
    line = split_line(input(2))
    wavevector = FractionVector(join(line(3:5)))
    
    line = split_line(input(3))
    degeneracy = int(line(3))
    
    line = split_line(input(4))
    energy = dble(line(3))
    
    line = split_line(input(5))
    coefficients = dble(line(3:))
    
    this = VscfState(subspace_id,wavevector,degeneracy,energy,coefficients)
  class default
    call err()
  end select
end subroutine

function write_VscfState(this) result(output)
  implicit none
  
  class(VscfState), intent(in) :: this
  type(String), allocatable    :: output(:)
  
  select type(this); type is(VscfState)
    output = [ 'Subspace     : '//this%subspace_id, &
             & 'Wavevector   : '//this%wavevector,  &
             & 'Degeneracy   : '//this%degeneracy,  &
             & 'Energy       : '//this%energy,      &
             & 'Coefficients : '//this%coefficients ]
  class default
    call err()
  end select
end function

function new_VscfState_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(VscfState)          :: this
  
  call this%read(input)
end function

impure elemental function new_VscfState_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(VscfState)               :: this
  
  this = VscfState(str(input))
end function
end module
