! ======================================================================
! Provides a random number generator.
! ======================================================================
!    - Generates real numbers between 0 and 1.
!    - Can be given a consistent seed.
module random_module
  use precision_module
  use abstract_module
  use io_module
  implicit none
  
  private
  
  public :: RandomReal
  
  type, extends(NoDefaultConstructor) :: RandomReal
    logical, private :: initialised_ = .false.
    integer, private :: random_seed_
  contains
    procedure, public :: get_seed
    procedure, public :: random_number  => random_number_RandomReal
    procedure, public :: random_numbers => random_numbers_RandomReal
  end type
  
  interface RandomReal
    module procedure new_RandomReal
  end interface
contains

! ----------------------------------------------------------------------
! Constructor.
! Initialises the random number generator based on the current time.
! ----------------------------------------------------------------------
function new_RandomReal(seed) result(this)
  implicit none
  
  integer, intent(in), optional :: seed
  type(RandomReal)              :: this
  
  integer              :: random_size
  integer, allocatable :: random_seeds(:)
  character(10)        :: time
  
  integer :: i,ialloc
  
  ! Set the seed, either to the given seed if present, or to the
  !    current time in milliseconds.
  if (present(seed)) then
    this%random_seed_ = seed
  else
    ! Get the time in seconds (to three decimal places).
    call date_and_time(time=time)
    ! Convert to an integer in milliseconds.
    this%random_seed_ = nint(1000 * dble(time))
  endif
  
  ! Initialise random number generator from the seed.
  call random_seed(size=random_size)
  allocate(random_seeds(random_size), stat=ialloc); call err(ialloc)
  random_seeds = this%random_seed_
  do i=1,size(random_seeds)
    random_seeds(i) = random_seeds(i) + i
  enddo
  call random_seed(put=random_seeds)
  
  ! Set this%initialised_ to .true.
  this%initialised_ = .true.
end function

! ----------------------------------------------------------------------
! Return the seed used to initialise the generator.
! ----------------------------------------------------------------------
function get_seed(this) result(output)
  implicit none
  
  class(RandomReal), intent(in) :: this
  integer                       :: output
  
  output = this%random_seed_
end function

! ----------------------------------------------------------------------
! Return a random number, or an array of random numbers.
! ----------------------------------------------------------------------
function random_number_RandomReal(this) result(output)
  implicit none
  
  class(RandomReal), intent(in) :: this
  real(dp)                      :: output
  
  if (.not. this%initialised_) then
    call print_line(CODE_ERROR//': Calling random number generator before it &
       &has been initialised.')
    call err()
  endif
  
  call random_number(output)
end function

function random_numbers_RandomReal(this,no_numbers) result(output)
  implicit none
  
  class(RandomReal), intent(in) :: this
  integer,           intent(in) :: no_numbers
  real(dp), allocatable         :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(no_numbers), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    output(i) = this%random_number()
  enddo
end function
end module
