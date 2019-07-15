! ======================================================================
! Performs a Pulay self-consistency solver.
! ======================================================================
! See the example function at the bottom of this file for how to use.
module pulay_module
  use precision_module
  use abstract_module
  use io_module
  use linear_algebra_module
  implicit none
  
  private
  
  public :: PulaySolver
  public :: pulay_solver_example
  
  type, extends(NoDefaultConstructor) :: PulaySolver
    ! Starting parameters.
    integer,  private :: pre_pulay_iterations_
    real(dp), private :: pre_pulay_damping_
    integer,  private :: max_pulay_iterations_
    
    ! The current set of x and f(x).
    type(RealVector), allocatable, private :: inputs_(:)
    type(RealVector), allocatable, private :: outputs_(:)
    
  contains
    procedure, public :: get_input
    procedure, public :: set_output
  end type
  
  interface PulaySolver
    module procedure new_PulaySolver
  end interface
contains

function new_PulaySolver(pre_pulay_iterations,pre_pulay_damping, &
   & max_pulay_iterations,initial_input) result(this)
  implicit none
  
  integer,  intent(in) :: pre_pulay_iterations
  real(dp), intent(in) :: pre_pulay_damping
  integer,  intent(in) :: max_pulay_iterations
  real(dp), intent(in) :: initial_input(:)
  type(PulaySolver)    :: this
  
  integer :: ialloc
  
  this%pre_pulay_iterations_ = pre_pulay_iterations
  this%pre_pulay_damping_ = pre_pulay_damping
  this%max_pulay_iterations_ = max_pulay_iterations
  
  this%inputs_ = [vec(initial_input)]
  allocate(this%outputs_(0), stat=ialloc); call err(ialloc)
end function

function get_input(this) result(output)
  implicit none
  
  class(PulaySolver), intent(in) :: this
  real(dp), allocatable          :: output(:)
  
  output = dble(this%inputs_(size(this%inputs_)))
end function

subroutine set_output(this,input)
  implicit none
  
  class(PulaySolver), intent(inout) :: this
  real(dp),           intent(in)    :: input(:)
  
  integer :: first_pulay_step
  integer :: last_step
  
  type(RealVector) :: new_input
  
  this%outputs_ = [this%outputs_, vec(input)]
  last_step = size(this%outputs_)
  if (last_step<this%pre_pulay_iterations_) then
    new_input = this%pre_pulay_damping_     * this%inputs_(last_step) &
            & + (1-this%pre_pulay_damping_) * this%outputs_(last_step)
  else
    first_pulay_step = max(1, last_step-this%max_pulay_iterations_+1)
    new_input = pulay( this%inputs_(first_pulay_step:last_step), &
                     & this%outputs_(first_pulay_step:last_step) )
  endif
  this%inputs_ = [this%inputs_, new_input]
end subroutine

! ----------------------------------------------------------------------
! An example demonstrating the use of the PulaySolver class.
! ----------------------------------------------------------------------
function example_function(input) result(output)
  implicit none
  
  real(dp), intent(in)  :: input(:)
  real(dp), allocatable :: output(:)
  
  type(IntMatrix) :: matrix
  
  ! The only eivenvector is (1,0,0).
  matrix = mat( [ 1,  0,  0,    &
              &   0,  0,  1,    &
              &   0, -1,  0  ], &
              & 3, 3            )
  
  output = dble(matrix*vec(input))
end function

! Finds x such that example_function(x)=x.
subroutine pulay_solver_example()
  implicit none
  
  integer               :: pre_pulay_iterations
  real(dp)              :: pre_pulay_damping
  integer               :: max_pulay_iterations
  real(dp), allocatable :: initial_input(:)
  
  type(PulaySolver) :: solver
  
  real(dp), allocatable :: input(:)
  real(dp), allocatable :: output(:)
  
  integer :: i
  
  pre_pulay_iterations = 2
  pre_pulay_damping    = 0.9_dp
  max_pulay_iterations = 20
  
  initial_input = [0.5_dp, 0.3_dp, 0.7_dp]
  
  ! Initialise solver.
  solver = PulaySolver( pre_pulay_iterations, &
                      & pre_pulay_damping,    &
                      & max_pulay_iterations, &
                      & initial_input         )
  
  i = 0
  do
    ! Get x from the solver.
    input = solver%get_input()
    
    ! Calculate f(x).
    output = example_function(input)
    
    ! Feed f(x) into the solver.
    call solver%set_output(output)
    
    ! Print progress.
    i = i+1
    call print_line('Step '//left_pad(i,'aaaaaa')//', &
                     &x= '//input//', &
                     &f= '//output)
    
    ! Check for convergence.
    if (l2_norm(vec(output-input))<1e-10_dp) then
      call print_line('Converged')
      exit
    endif
  enddo
end subroutine
end module
