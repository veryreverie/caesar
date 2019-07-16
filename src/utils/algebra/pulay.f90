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
  
  ! Set the output.
  this%outputs_ = [this%outputs_, vec(input)]
  
  ! Construct the next input, using either a damped iterative or Pulay scheme.
  last_step = size(this%outputs_)
  if (last_step<this%pre_pulay_iterations_) then
    new_input = this%pre_pulay_damping_     * this%inputs_(last_step) &
            & + (1-this%pre_pulay_damping_) * this%outputs_(last_step)
  else
    first_pulay_step = max(1, last_step-this%max_pulay_iterations_+1)
    new_input = pulay( this%inputs_(first_pulay_step:last_step), &
                     & this%outputs_(first_pulay_step:last_step) )
    
    ! Check if the Pulay scheme failed. If so, restart from a damped iterative
    !    scheme but using the latest values.
    if (size(new_input)==0) then
      call print_line(WARNING//': Pulay scheme failed. &
         &Restarting from damped iterative scheme.')
      this%outputs_ = this%outputs_(size(this%outputs_):)
      this%inputs_ = this%inputs_(size(this%inputs_):)
      last_step = size(this%outputs_)
      new_input = this%pre_pulay_damping_     * this%inputs_(last_step) &
              & + (1-this%pre_pulay_damping_) * this%outputs_(last_step)
    endif
  endif
  
  ! Set the next input.
  this%inputs_ = [this%inputs_, new_input]
end subroutine

! A Pulay scheme to find self-consistent solution to f(x)=x.
! N.B. The routine will fail if the error matrix has less than full rank.
! In this case, the output will be an empty array.
function pulay(input_vectors,output_vectors) result(output)
  implicit none
  
  type(RealVector), intent(in) :: input_vectors(:)
  type(RealVector), intent(in) :: output_vectors(:)
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
  
  if (size(input_vectors)/=size(output_vectors)) then
    call print_line(ERROR//': Input and output vectors do not match.')
    call err()
  endif
  
  m = size(input_vectors(1))
  n = size(input_vectors)
  
  if (any([(size(input_vectors(i)),i=1,n)]/=m)) then
    call print_line(ERROR//': Input vectors inconsistent.')
  elseif (any([(size(output_vectors(i)),i=1,n)]/=m)) then
    call print_line(ERROR//': Output vectors inconsistent.')
  endif
  
  ! Construct errors, e_i = f(x_i)-x_i.
  errors = output_vectors - input_vectors
  
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
  output = sum(coefficients(:n)*input_vectors)
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
