! ======================================================================
! A basis of harmonic eigenstates along a single real normal mode.
! ======================================================================
module single_harmonic_basis_module
  use common_module
  
  use single_mode_state_module
  implicit none
  
  private
  
  public :: SingleHarmonicBasis
  
  type, extends(Stringsable) :: SingleHarmonicBasis
    type(SingleModeState), allocatable :: states_(:)
  contains
    procedure, public :: state
    ! I/O.
    procedure, public :: read  => read_SingleHarmonicBasis
    procedure, public :: write => write_SingleHarmonicBasis
  end type
  
  interface SingleHarmonicBasis
    module procedure new_SingleHarmonicBasis
    module procedure new_SingleHarmonicBasis_Strings
    module procedure new_SingleHarmonicBasis_StringArray
  end interface
contains

! Constructor.
function new_SingleHarmonicBasis(states) result(this)
  implicit none
  
  type(SingleModeState), intent(in) :: states(:)
  type(SingleHarmonicBasis)         :: this
  
  this%states_ = states
end function

! this%state(i) returns |i>.
function state(this,i) result(output)
  implicit none
  
  class(SingleHarmonicBasis), intent(in) :: this
  integer,                    intent(in) :: i
  type(SingleModeState)                  :: output
  
  if (i<0 .or. i>=size(this%states_)) then
    call print_line(ERROR//': |'//i//'> is not a valid state.')
    call err()
  endif
  
  output = this%states_(i+1)
end function

! ----------------------------------------------------------------------
! Generates the harmonic basis functions along a specific mode.
! ----------------------------------------------------------------------
! Uses the recurrence relation:
!   |n> = sqrt(2*freq/n) u |n-1> - sqrt((n-1)/n) |n-2>
! N.B. basis_functions(i) = |i-1> because |0> is a state.
function generate_harmonic_basis(mode,cutoff,frequency) &
   & result(output)
  implicit none
  
  type(RealMode), intent(in)           :: mode
  integer,        intent(in)           :: cutoff
  real(dp),       intent(in), optional :: frequency
  type(SingleModeState), allocatable   :: output(:)
  
  real(dp)              :: freq
  real(dp)              :: normalisation
  real(dp), allocatable :: coefficients(:,:)
  
  integer :: i,ialloc
  
  if (present(frequency)) then
    freq = frequency
  else
    freq = mode%frequency
  endif
  
  normalisation = (freq/PI)**0.25_dp
  
  allocate(coefficients(cutoff+1,cutoff+1), stat=ialloc); call err(ialloc)
  coefficients = 0.0_dp
  
  ! Calculate |0> coefficients.
  if (cutoff >= 0) then
    coefficients(1,1) = normalisation
  endif
  
  ! Calculate |1> coefficients.
  if (cutoff >= 1) then
    ! coefficients(1,2)=0.
    coefficients(2,2) = normalisation * sqrt(2*freq)
  endif
  
  ! Calculate |2> to |cutoff> coefficients.
  do i=3,size(output)
    ! Calculate |i-1> from |i-2> and |i-3>.
    ! |i-1> = sqrt(2w/(i-1))u|i-2> - sqrt((i-2)/(i-1))|i-3>.
    ! N.B. coefficients(1:i,i) = |i-1>.
    
    ! sqrt(2w/(i-1)) u |i-2>.
    coefficients(2:i,i) = coefficients(2:i,i)     &
                      & + sqrt(2*frequency/(i-1)) &
                      & * coefficients(:i-1,i-1)
    
    ! -sqrt((i-2)/(i-1)) |i-3>.
    coefficients(:i-2,i) = coefficients(:i-2,i)   &
                       & - sqrt((i-2.0_dp)/(i-1)) &
                       & * coefficients(:i-2,i-2)
  enddo
  
  ! Construct output from coefficients.
  allocate(output(cutoff+1), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    output(i) = SingleModeState(mode, frequency, coefficients(:i,i))
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_SingleHarmonicBasis(this,input)
  implicit none
  
  class(SingleHarmonicBasis), intent(out) :: this
  type(String),               intent(in)  :: input(:)
  
  type(SingleModeState), allocatable :: states(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(SingleHarmonicBasis)
    allocate(states(size(input)), stat=ialloc); call err(ialloc)
    do i=1,size(input)
      line = split_line(input(i), delimiter='>')
      states(i) = SingleModeState(line(2))
    enddo
    
    this = SingleHarmonicBasis(states)
  class default
    call err()
  end select
end subroutine

function write_SingleHarmonicBasis(this) result(output)
  implicit none
  
  class(SingleHarmonicBasis), intent(in) :: this
  type(String), allocatable              :: output(:)
  
  integer :: i
  
  select type(this); type is(SingleHarmonicBasis)
    output = str(this%states_)
    
    do i=1,size(output)
      output(i) = '|'//i-1//'> '//output(i)
    enddo
  class default
    call err()
  end select
end function

function new_SingleHarmonicBasis_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)  :: input(:)
  type(SingleHarmonicBasis) :: this
  
  call this%read(input)
end function

impure elemental function new_SingleHarmonicBasis_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SingleHarmonicBasis)     :: this
  
  this = SingleHarmonicBasis(str(input))
end function
end module
