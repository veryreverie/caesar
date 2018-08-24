! ======================================================================
! A basis of harmonic eigenstates along a single complex normal mode.
! ======================================================================
module harmonic_mode_basis_module
  use common_module
  
  use mode_state_module
  implicit none
  
  private
  
  public :: HarmonicModeBasis
  
  type, extends(Stringsable) :: HarmonicModeBasis
    integer                      :: mode_id
    real(dp)                     :: frequency
    type(ModeState), allocatable :: states_(:)
  contains
    procedure, public :: state
    ! I/O.
    procedure, public :: read  => read_HarmonicModeBasis
    procedure, public :: write => write_HarmonicModeBasis
  end type
  
  interface HarmonicModeBasis
    module procedure new_HarmonicModeBasis
    module procedure new_HarmonicModeBasis_Strings
    module procedure new_HarmonicModeBasis_StringArray
  end interface
contains

! Constructor.
function new_HarmonicModeBasis(mode_id,frequency,states) result(this)
  implicit none
  
  integer,         intent(in) :: mode_id
  real(dp),        intent(in) :: frequency
  type(ModeState), intent(in) :: states(:)
  type(HarmonicModeBasis)     :: this
  
  this%mode_id   = mode_id
  this%frequency = frequency
  this%states_   = states
end function

! this%state(i) returns |i>.
function state(this,i) result(output)
  implicit none
  
  class(HarmonicModeBasis), intent(in) :: this
  integer,                  intent(in) :: i
  type(ModeState)                      :: output
  
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
  
  type(ComplexMode), intent(in)           :: mode
  integer,           intent(in)           :: cutoff
  real(dp),          intent(in), optional :: frequency
  type(HarmonicModeBasis)                 :: output
  
  real(dp)                     :: freq
  real(dp),        allocatable :: coefficients(:,:)
  type(ModeState), allocatable :: states(:)
  
  integer :: i,ialloc
  
  if (present(frequency)) then
    freq = frequency
  else
    freq = mode%frequency
  endif
  
  ! Calculate coefficients.
  ! N.B. the normalisation factor of (m*freq/pi)^(1/4) is already included
  !    in the definition of the ModeState.
  allocate(coefficients(cutoff+1,cutoff+1), stat=ialloc); call err(ialloc)
  coefficients = 0.0_dp
  
  ! Calculate |0> coefficients.
  if (cutoff >= 0) then
    coefficients(1,1) = 1.0_dp
  endif
  
  ! Calculate |1> coefficients.
  if (cutoff >= 1) then
    ! coefficients(1,2)=0.
    coefficients(2,2) = sqrt(2*freq)
  endif
  
  ! Calculate |2> to |cutoff> coefficients.
  do i=3,cutoff+1
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
  
  ! Construct basis states from coefficients.
  allocate(states(cutoff+1), stat=ialloc); call err(ialloc)
  do i=1,size(states)
    states(i) = ModeState(mode, freq, coefficients(:i,i))
  enddo
  
  ! Construct output from basis states.
  output = HarmonicModeBasis( mode_id   = mode%id, &
                            & frequency = freq,    &
                            & states    = states   )
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_HarmonicModeBasis(this,input)
  implicit none
  
  class(HarmonicModeBasis), intent(out) :: this
  type(String),             intent(in)  :: input(:)
  
  integer                      :: mode_id
  real(dp)                     :: frequency
  type(ModeState), allocatable :: states(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(HarmonicModeBasis)
    line = split_line(input(1))
    mode_id = int(line(4))
    
    line = split_line(input(2))
    frequency = dble(line(3))
    
    allocate(states(size(input)-2), stat=ialloc); call err(ialloc)
    do i=1,size(input)-2
      line = split_line(input(i+2), delimiter='>')
      states(i) = ModeState(line(2))
    enddo
    
    this = HarmonicModeBasis(mode_id, frequency, states)
  class default
    call err()
  end select
end subroutine

function write_HarmonicModeBasis(this) result(output)
  implicit none
  
  class(HarmonicModeBasis), intent(in) :: this
  type(String), allocatable            :: output(:)
  
  integer :: i
  
  select type(this); type is(HarmonicModeBasis)
    output = [ 'Mode ID   : '//this%mode_id,   &
             & 'Frequency : '//this%frequency, &
             & str(this%states_)               ]
    
    do i=3,size(output)
      output(i) = '|'//i-1//'> '//output(i)
    enddo
  class default
    call err()
  end select
end function

function new_HarmonicModeBasis_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(HarmonicModeBasis)  :: this
  
  call this%read(input)
end function

impure elemental function new_HarmonicModeBasis_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(HarmonicModeBasis)       :: this
  
  this = HarmonicModeBasis(str(input))
end function
end module
