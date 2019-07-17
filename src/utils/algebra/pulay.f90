! ======================================================================
! Performs a Pulay self-consistency solver, to find f(x)=x.
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
    logical,  private :: bound_at_zero_
    
    ! The current set of x and f(x).
    type(RealVector), allocatable, private :: xs_(:)
    type(RealVector), allocatable, private :: fs_(:)
    
  contains
    ! Get x_i, from which f(x_i) should be calculated.
    procedure, public :: get_x
    
    ! Set f(x_i).
    ! This also calculates x_{i+1}.
    procedure, public :: set_f
  end type
  
  interface PulaySolver
    module procedure new_PulaySolver
  end interface
contains

function new_PulaySolver(pre_pulay_iterations,pre_pulay_damping, &
   & max_pulay_iterations,initial_x,bound_at_zero) result(this)
  implicit none
  
  integer,  intent(in)           :: pre_pulay_iterations
  real(dp), intent(in)           :: pre_pulay_damping
  integer,  intent(in)           :: max_pulay_iterations
  real(dp), intent(in)           :: initial_x(:)
  logical,  intent(in), optional :: bound_at_zero
  type(PulaySolver)              :: this
  
  integer :: ialloc
  
  this%pre_pulay_iterations_ = pre_pulay_iterations
  this%pre_pulay_damping_ = pre_pulay_damping
  this%max_pulay_iterations_ = max_pulay_iterations
  if (present(bound_at_zero)) then
    this%bound_at_zero_ = bound_at_zero
  else
    this%bound_at_zero_ = .false.
  endif
  
  this%xs_ = [vec(initial_x)]
  allocate(this%fs_(0), stat=ialloc); call err(ialloc)
end function

! Get x_i
function get_x(this) result(output)
  implicit none
  
  class(PulaySolver), intent(in) :: this
  real(dp), allocatable          :: output(:)
  
  output = dble(this%xs_(size(this%xs_)))
end function

! Set f(x_i). Also calculates x_{i+1}.
subroutine set_f(this,f)
  implicit none
  
  class(PulaySolver), intent(inout) :: this
  real(dp),           intent(in)    :: f(:)
  
  integer :: first_pulay_step
  integer :: last_step
  
  type(RealVector) :: x
  
  ! Record f(x_i).
  this%fs_ = [this%fs_, vec(f)]
  
  ! Construct x_{i+1}.
  ! This uses a damped iterative scheme if i<pre_pulay_iterations, and a Pulay
  !    scheme otherwise.
  ! The Pulay scheme may fail, in which case i is reset to 1, and the scheme
  !    re-starts using the latest values of x and f(x).
  last_step = size(this%fs_)
  if (last_step<this%pre_pulay_iterations_) then
    ! The damped iterative scheme.
    ! x_{i+1} = ax_i + bf(x_i), where a=pre_pulay_damping and b=1-a.
    x = this%pre_pulay_damping_     * this%xs_(last_step) &
    & + (1-this%pre_pulay_damping_) * this%fs_(last_step)
  else
    ! The Pulay scheme.
    first_pulay_step = max(1, last_step-this%max_pulay_iterations_+1)
    x = pulay( this%xs_(first_pulay_step:last_step), &
                     & this%fs_(first_pulay_step:last_step) )
    
    ! Check if the Pulay scheme failed.
    if (size(x)==0) then
      call print_line(WARNING//': Pulay scheme failed. &
         &Restarting from damped iterative scheme.')
      this%fs_ = this%fs_(size(this%fs_):)
      this%xs_ = this%xs_(size(this%xs_):)
      last_step = size(this%fs_)
      x = this%pre_pulay_damping_     * this%xs_(last_step) &
      & + (1-this%pre_pulay_damping_) * this%fs_(last_step)
    endif
  endif
  
  ! Bound the input at zero if required.
  ! N.B. this is done by limiting the input to be at least half the previous
  !    output, to avoid sharp changes.
  if (this%bound_at_zero_) then
    x = vec(max(dble(x), f/2.0_dp))
  endif
  
  ! Set the next input.
  this%xs_ = [this%xs_, x]
end subroutine

! A Pulay scheme to find self-consistent solution to f(x)=x.
! N.B. The routine will fail if the error matrix has less than full rank.
! In this case, the output will be an empty array.
function pulay(xs,fs) result(output)
  implicit none
  
  type(RealVector), intent(in) :: xs(:)
  type(RealVector), intent(in) :: fs(:)
  type(RealVector)             :: output
  
  type(Realvector), allocatable :: errors(:)
  real(dp),         allocatable :: error_matrix(:,:)
  real(dp),         allocatable :: lagrange_vector(:)
  real(dp),         allocatable :: coefficients(:)
  
  integer :: m,n
  
  integer :: rank_info
  
  integer :: i,j,ialloc
  
  ! Each input vector x_i produces an output vector f(x_i) = x_i + e_i.
  ! The aim is to construct x_{n+1} such that e_{n+1}=0.
  ! If x_{n+1} = sum_{i=1}^n a_i x_i, it can be approximated that
  !    e_{n+1} = sum_{i=1}^n a_i e_i.
  ! Minimising e_{n+1} subject to sum_{i=1}^n a_i = 1 then becomes a linear
  !    least squares problem, with equations:
  ! sum_{i=1}^n e_i e_j a_j + l = 0
  ! sum_{i=1}^n a_i             = 1
  ! Or equivalently:
  ! ( e1.e1, e1.e2, ... , e1.en,  1  ) (a_1)   ( 0 )
  ! ( e2.e1, e2.e2, ... , e2.en,  1  ) (a_2)   ( 0 )
  ! (  ... ,  ... , ... ,  ... , ... ) (...) = (...)
  ! ( en.e1, en.e2, ... , en.en,  1  ) (a_n)   ( 0 )
  ! (   1  ,   1  ,  1  ,   1  ,  0  ) ( l )   ( 1 )
  
  if (size(xs)/=size(fs)) then
    call print_line(ERROR//': Input and output vectors do not match.')
    call err()
  endif
  
  m = size(xs(1))
  n = size(xs)
  
  if (any([(size(xs(i)),i=1,n)]/=m)) then
    call print_line(ERROR//': Input vectors inconsistent.')
  elseif (any([(size(fs(i)),i=1,n)]/=m)) then
    call print_line(ERROR//': Output vectors inconsistent.')
  endif
  
  ! Construct errors, e_i = f(x_i)-x_i.
  errors = fs - xs
  
  ! Construct error matrix.
  ! ( e1.e1, e1.e2, ... , e1.en,  1 )
  ! ( e2.e1, e2.e2, ... , e2.en,  1 )
  ! (  ... ,  ... , ... ,  ... , ...)
  ! ( en.e1, en.e2, ... , en.en,  1 )
  ! (   1  ,   1  , ... ,   1  ,  0 )
  allocate(error_matrix(n+1,n+1), stat=ialloc); call err(ialloc)
  do i=1,n
    do j=1,n
      error_matrix(j,i) = errors(j) * errors(i)
    enddo
  enddo
  error_matrix(:n , n+1) = 1
  error_matrix(n+1, :n ) = 1
  error_matrix(n+1, n+1) = 0
  
  ! Construct lagrange vector, equal to (0, 0, ..., 0, 1).
  lagrange_vector = [(0.0_dp,i=1,n), 1.0_dp]
  
  ! Perform least-squares optimisation to get coefficients, {a_i}.
  coefficients = dble(linear_least_squares( error_matrix,                &
                                          & lagrange_vector,             &
                                          & return_empty_on_error=.true. ))
  
  ! Return an empty vector if the linear least squares failed.
  if (size(coefficients)==0) then
    output = vec([real(dp)::])
    return
  endif
  
  ! Construct output frequencies as x_n = sum_{i=1}^n a_i x_i.
  output = sum(coefficients(:n)*xs)
end function

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
  real(dp), allocatable :: initial_x(:)
  
  type(PulaySolver) :: solver
  
  real(dp), allocatable :: x(:)
  real(dp), allocatable :: f(:)
  
  integer :: i
  
  pre_pulay_iterations = 2
  pre_pulay_damping    = 0.9_dp
  max_pulay_iterations = 20
  
  initial_x = [0.5_dp, 0.3_dp, 0.7_dp]
  
  ! Initialise solver.
  solver = PulaySolver( pre_pulay_iterations, &
                      & pre_pulay_damping,    &
                      & max_pulay_iterations, &
                      & initial_x             )
  
  i = 0
  do
    ! Get x from the solver.
    x = solver%get_x()
    
    ! Calculate f(x).
    f = example_function(x)
    
    ! Feed f(x) into the solver.
    call solver%set_f(f)
    
    ! Print progress.
    i = i+1
    call print_line('Step '//left_pad(i,'aaaaaa')//', &
                     &x= '//x//', &
                     &f= '//f)
    
    ! Check for convergence.
    if (l2_norm(vec(f-x))<1e-10_dp) then
      call print_line('Converged')
      exit
    endif
  enddo
end subroutine
end module
