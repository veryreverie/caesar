! ======================================================================
! A sum of monomial states and coefficients.
! ======================================================================
module polynomial_state_module
  use common_module
  
  use subspace_state_module
  use monomial_state_module
  implicit none
  
  private
  
  public :: PolynomialState
  public :: size
  public :: braket_PolynomialState
  public :: kinetic_energy_PolynomialState
  public :: harmonic_potential_energy_PolynomialState
  
  type, extends(SubspaceState) :: PolynomialState
    type(MonomialState), allocatable :: states(:)
    real(dp),            allocatable :: coefficients(:)
  contains
    procedure, public :: read  => read_PolynomialState
    procedure, public :: write => write_PolynomialState
  end type
  
  interface PolynomialState
    module procedure new_PolynomialState
    module procedure new_PolynomialState_Strings
    module procedure new_PolynomialState_StringArray
  end interface
  
  interface size
    module procedure size_PolynomialState
  end interface
  
  interface braket_PolynomialState
    module procedure braket_PolynomialStates
    module procedure braket_PolynomialStates_ComplexMonomial
  end interface
contains

! Constructor and size function
function new_PolynomialState(subspace_id,states,coefficients) result(this)
  implicit none
  
  integer                          :: subspace_id
  type(MonomialState), allocatable :: states(:)
  real(dp),            allocatable :: coefficients(:)
  type(PolynomialState)            :: this
  
  if (size(states)==0) then
    call print_line(CODE_ERROR//': No states.')
    call err()
  elseif (size(states)/=size(coefficients)) then
    call print_line(CODE_ERROR//': States and coefficients do not match.')
    call err()
  endif
  
  this%subspace_id  = subspace_id
  this%states       = states
  this%coefficients = coefficients
end function

function size_PolynomialState(this) result(output)
  implicit none
  
  type(PolynomialState), intent(in) :: this
  integer                           :: output
  
  output = size(this%states)
end function

! ----------------------------------------------------------------------
! Integrals of the form <bra|ket> and <bra|potential|ket>.
! ----------------------------------------------------------------------
impure elemental function braket_PolynomialStates(bra,ket) result(output)
  implicit none
  
  type(PolynomialState), intent(in) :: bra
  type(PolynomialState), intent(in) :: ket
  real(dp)                          :: output
  
  integer :: i,j
  
  output = 0.0_dp
  do i=1,size(bra)
    do j=1,size(ket)
      output = output                                            &
           & + braket_MonomialState(bra%states(i),ket%states(j)) &
           & * bra%coefficients(i)                               &
           & * ket%coefficients(j)
    enddo
  enddo
end function

impure elemental function braket_PolynomialStates_ComplexMonomial(bra,ket, &
   & monomial,subspace,supercell) result(output)
  implicit none
  
  type(PolynomialState),    intent(in) :: bra
  type(PolynomialState),    intent(in) :: ket
  type(ComplexMonomial),    intent(in) :: monomial
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  type(ComplexMonomial)                :: output
  
  type(ComplexMonomial) :: integrated_monomial
  
  integer :: i,j
  
  do i=1,size(bra)
    do j=1,size(ket)
      integrated_monomial = braket_MonomialState( bra%states(i),   &
                        &                         ket%states(j),   &
                        &                         monomial,        &
                        &                         subspace,        &
                        &                         supercell      ) &
                        & * bra%coefficients(i)                    &
                        & * ket%coefficients(j)
      if (i==1 .and. j==1) then
        output = integrated_monomial
      else
        output%coefficient = output%coefficient &
                         & + integrated_monomial%coefficient
      endif
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Evaluates <bra|T|ket>, where T is the kinetic energy operator.
! Gives the result per primitive cell.
! ----------------------------------------------------------------------
function kinetic_energy_PolynomialState(bra,ket,subspace,supercell) &
   & result(output)
  implicit none
  
  type(PolynomialState),    intent(in) :: bra
  type(PolynomialState),    intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  real(dp)                             :: output
  
  integer :: i,j
  
  output = 0.0_dp
  do i=1,size(bra)
    do j=1,size(ket)
      output = output                                         &
           & + kinetic_energy_MonomialState( bra%states(i),   &
           &                                 ket%states(j),   &
           &                                 subspace,        &
           &                                 supercell      ) &
           & * bra%coefficients(i)                            &
           & * ket%coefficients(j)
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Evaluates <bra|V|ket>, where V is the harmonic potential energy operator.
! Gives the result per primitive cell.
! ----------------------------------------------------------------------
function harmonic_potential_energy_PolynomialState(bra,ket,subspace, &
   & supercell) result(output)
  implicit none
  
  type(PolynomialState),    intent(in) :: bra
  type(PolynomialState),    intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  real(dp)                             :: output
  
  integer :: i,j
  
  output = 0.0_dp
  do i=1,size(bra)
    do j=1,size(ket)
      output = output                                                    &
           & + harmonic_potential_energy_MonomialState( bra%states(i),   &
           &                                            ket%states(j),   &
           &                                            subspace,        &
           &                                            supercell      ) &
           & * bra%coefficients(i)                                       &
           & * ket%coefficients(j)
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_PolynomialState(this,input)
  implicit none
  
  class(PolynomialState), intent(out) :: this
  type(String),           intent(in)  :: input(:)
  
  integer               :: subspace_id
  real(dp)              :: frequency
  type(ComplexMonomial) :: state
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(PolynomialState)
    ! TODO
    call err()
  class default
    call err()
  end select
end subroutine

function write_PolynomialState(this) result(output)
  implicit none
  
  class(PolynomialState), intent(in) :: this
  type(String), allocatable          :: output(:)
  
  type(String) :: state_string
  
  select type(this); type is(PolynomialState)
    ! TODO
    call err()
  class default
    call err()
  end select
end function

function new_PolynomialState_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(PolynomialState)    :: this
  
  call this%read(input)
end function

impure elemental function new_PolynomialState_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(PolynomialState)         :: this
  
  this = PolynomialState(str(input))
end function
end module
