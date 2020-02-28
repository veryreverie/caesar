! ======================================================================
! Performs a Pulay self-consistency solver, to minimise F(x)
!    subject to f(x)=x.
! ======================================================================
! See the example function at the bottom of this file for how to use.
! Since the Pulay scheme is invariant under the ordering of its inputs,
!    only the last max_pulay_iterations are stored, and the arrays are
!    continually overwritten in a sequential 1,2,...,n,1,2,...,n,1,... manner.
! This is acheived by tracking i_, and taking modulo(i_,max_pulay_iterations)
!    to identify where each iteration should be stored.
module pulay_module
  use precision_module
  use abstract_module
  use io_module
  use linear_algebra_module
  use hermitian_eigenstuff_module
  use random_module
  implicit none
  
  private
  
  public :: ConvergenceData
  public :: PulaySolver
  public :: pulay_solver_example
  
  type, extends(NoDefaultConstructor) :: ConvergenceData
    integer  :: pre_pulay_iterations
    real(dp) :: pre_pulay_damping
    integer  :: max_pulay_iterations
    real(dp) :: iterative_mixing
    real(dp) :: iterative_damping
    real(dp) :: pulay_noise
    real(dp) :: energy_convergence
    integer  :: no_converged_calculations
  end type
  
  interface ConvergenceData
    module procedure new_ConvergenceData
  end interface
  
  type, extends(NoDefaultConstructor) :: PulaySolver
    ! Starting parameters.
    type(ConvergenceData), private :: data_
    type(RandomReal),      private :: random_generator_
    logical,               private :: bound_at_zero_
    
    ! The current set of x, f(x) and the free energy F(x).
    type(RealVector), allocatable, private :: xs_(:)
    type(RealVector), allocatable, private :: fs_(:)
    real(dp),         allocatable, private :: free_energies_(:)
    
    integer, private :: i_
  contains
    ! Get x_i, from which f(x_i) and F(x_i) should be calculated.
    procedure, public :: get_x
    
    ! Set f(x_i) and F(x_i).
    ! This also calculates x_{i+1}.
    procedure, public :: set_f
    
    ! The Pulay scheme itself.
    procedure, private :: pulay_
    
    ! Errors for printing progress.
    procedure, public :: self_consistency_error
    procedure, public :: coefficient_change
    procedure, public :: free_energy_change
    
    ! Check if the Pulay scheme has converged.
    procedure, public :: converged
    
    procedure, private :: i
  end type
  
  interface PulaySolver
    module procedure new_PulaySolver
  end interface
contains

impure elemental function new_ConvergenceData(pre_pulay_iterations, &
   & pre_pulay_damping,max_pulay_iterations,iterative_mixing,       &
   & iterative_damping,pulay_noise,energy_convergence,              &
   & no_converged_calculations) result(this) 
  implicit none
  
  integer,  intent(in)  :: pre_pulay_iterations
  real(dp), intent(in)  :: pre_pulay_damping
  integer,  intent(in)  :: max_pulay_iterations
  real(dp), intent(in)  :: iterative_mixing
  real(dp), intent(in)  :: iterative_damping
  real(dp), intent(in)  :: pulay_noise
  real(dp), intent(in)  :: energy_convergence
  integer,  intent(in)  :: no_converged_calculations
  type(ConvergenceData) :: this
  
  if (pre_pulay_iterations<2) then
    call print_line(ERROR//': pre_pulay_iterations must be at least 2.')
    call quit()
  elseif (max_pulay_iterations<2) then
    call print_line(ERROR//': max_pulay_iterations must be at least 2.')
    call quit()
  endif
  
  this%pre_pulay_iterations      = pre_pulay_iterations
  this%pre_pulay_damping         = pre_pulay_damping
  this%max_pulay_iterations      = max_pulay_iterations
  this%iterative_mixing          = iterative_mixing
  this%iterative_damping         = iterative_damping
  this%pulay_noise               = pulay_noise
  this%energy_convergence        = energy_convergence
  this%no_converged_calculations = no_converged_calculations
end function

function new_PulaySolver(convergence_data,random_generator,initial_x, &
   & bound_at_zero) result(this)
  implicit none
  
  type(ConvergenceData), intent(in)           :: convergence_data
  type(RandomReal),      intent(in)           :: random_generator
  real(dp),              intent(in)           :: initial_x(:)
  logical,               intent(in), optional :: bound_at_zero
  type(PulaySolver)                           :: this
  
  integer :: array_size
  
  integer :: ialloc
  
  this%data_             = convergence_data
  this%random_generator_ = random_generator
  if (present(bound_at_zero)) then
    this%bound_at_zero_ = bound_at_zero
  else
    this%bound_at_zero_ = .false.
  endif
  
  array_size = max( this%data_%max_pulay_iterations,     &
                  & this%data_%no_converged_calculations )
  
  allocate( this%xs_(array_size),            &
          & this%fs_(array_size),            &
          & this%free_energies_(array_size), &
          & stat=ialloc); call err(ialloc)
  this%xs_(1) = vec(initial_x)
  
  this%i_ = 1
end function

! Private index function.
function i(this,offset) result(output)
  implicit none
  
  class(PulaySolver), intent(in)           :: this
  integer,            intent(in), optional :: offset
  integer                                  :: output
  
  if (present(offset)) then
    if (offset>0) then
      call print_line(CODE_ERROR//': offset > 0.')
      call err()
    elseif (this%i_+offset<1) then
      call print_line(CODE_ERROR//': i+offset <= 0.')
      call err()
    endif
  endif
  
  if (present(offset)) then
    output = modulo( this%i_+offset-1,                 &
         &           this%data_%max_pulay_iterations ) &
         & + 1
  else
    output = modulo( this%i_-1,                        &
         &           this%data_%max_pulay_iterations ) &
         & + 1
  endif
end function

! Get x_i
function get_x(this) result(output)
  implicit none
  
  class(PulaySolver), intent(in) :: this
  real(dp), allocatable          :: output(:)
  
  output = dble(this%xs_(this%i()))
end function

! Set f(x_i) and F(x_i). Also calculates x_{i+1}.
subroutine set_f(this,f,free_energy)
  implicit none
  
  class(PulaySolver), intent(inout) :: this
  real(dp),           intent(in)    :: f(:)
  real(dp),           intent(in)    :: free_energy
  
  integer              :: first_iteration
  integer              :: last_iteration
  integer, allocatable :: iterations(:)
  
  integer :: i
  
  ! Record f(x_i) and F(x_i).
  this%fs_(this%i()) = vec(f)
  this%free_energies_(this%i()) = free_energy
  
  ! Update i.
  this%i_ = this%i_+1
  
  ! Construct x_{i+1}.
  ! This uses a damped iterative scheme if i<pre_pulay_iterations, and a Pulay
  !    scheme otherwise.
  ! The Pulay scheme may fail, in which case i is reset to 1, and the scheme
  !    re-starts using the latest values of x and f(x).
  if (this%i_<=this%data_%pre_pulay_iterations) then
    ! The damped iterative scheme.
    ! x_{i+1} = ax_i + bf(x_i), where a=pre_pulay_damping and b=1-a.
    this%xs_(this%i()) = this%data_%pre_pulay_damping     &
                     & * this%xs_(this%i(offset=-1))      &
                     & + (1-this%data_%pre_pulay_damping) &
                     & * this%fs_(this%i(offset=-1))
  else
    ! The Pulay scheme.
    
    ! Calculate the positions of the first and last iteration to pass to the
    !    Pulay scheme.
    if (this%i_<=this%data_%max_pulay_iterations) then
      ! If i_ is less than max_pulay_iterations, only the first i_-1 iterations
      !    are used, as these are all that are available.
      first_iteration = 1
      last_iteration  = this%i_-1
    else
      ! Otherwise, the last max_pulay_iterations are used.
      first_iteration = this%i(-this%data_%max_pulay_iterations)
      last_iteration = this%i(-1)
    endif
    
    ! Construct the array [first_iteration,...,last_iteration], taking into
    !    account possible wrapping around the end of this%xs_ etc.
    if (first_iteration<=last_iteration) then
      iterations = [(i,i=first_iteration,last_iteration)]
    else
      iterations = [ [(i,i=first_iteration,size(this%xs_))], &
                   & [(i,i=1,last_iteration)]                ]
    endif
    
    ! Calculate the next iteration using the Pulay scheme.
    this%xs_(this%i()) = this%pulay_(iterations)
  endif
  
  ! Bound the input at zero if required.
  ! N.B. this is done by limiting the input to be at least half the previous
  !    output, to avoid sharp changes.
  if (this%bound_at_zero_) then
    this%xs_(this%i()) = vec(max(dble(this%xs_(this%i())), f/2.0_dp))
  endif
end subroutine

! A Pulay scheme to find self-consistent solution to f(x)=x.
! If the Pulay scheme is ill-conditioned (i.e. there is a subspace of
!    self-consistent solutions) then a gradient-descent scheme will be run
!    within the self-consistent subspace.
function pulay_(this,iterations) result(output)
  implicit none
  
  class(PulaySolver), intent(in) :: this
  integer,            intent(in) :: iterations(:)
  type(RealVector)               :: output
  
  type(RealVector), allocatable :: errors(:)
  real(dp),         allocatable :: error_matrix(:,:)
  real(dp),         allocatable :: coefficients(:)
  
  integer :: min_guess
  integer :: max_guess
  
  type(RealVector) :: projection
  real(dp)         :: projection_norm
  
  real(dp) :: maximum_prefactor
  real(dp) :: prefactor
  
  type(SymmetricEigenstuff), allocatable :: unsorted_estuff(:)
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  
  integer :: first,last
  
  integer :: m,n
  
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
  
  m = size(this%xs_(iterations(1)))
  n = size(iterations)
  
  if (any([(size(this%xs_(iterations(i))),i=1,n)]/=m)) then
    call print_line(ERROR//': Input vectors inconsistent.')
  elseif (any([(size(this%fs_(iterations(i))),i=1,n)]/=m)) then
    call print_line(ERROR//': Output vectors inconsistent.')
  endif
  
  ! Construct errors, e_i = f(x_i)-x_i.
  errors = this%fs_(iterations) - this%xs_(iterations)
  
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
  
  min_guess = minloc([(error_matrix(i,i), i=1, n)], 1)
  max_guess = maxloc([(error_matrix(i,i), i=1, n)], 1)
  
  ! Check for over-convergence.
  if (error_matrix(max_guess,max_guess)<1e-300_dp) then
    output = this%xs_(iterations(min_guess))
    return
  endif
  
  ! Scale error_matrix(:n,:n) such that the problem is well-conditioned.
  ! This has no effect on the result beyond numerical error.
  if (max_guess>0.01_dp*this%data_%energy_convergence**2) then
    error_matrix(:n,:n) = error_matrix(:n,:n) &
                      & / error_matrix(max_guess,max_guess)
  endif
  
  ! Diagonalise the error matrix.
  ! The i'th eigenvector is (c_i1,c_i2,...,c_in,w_i), where:
  !    - c_ij are coefficients,
  !    - w_i is a weight.
  ! If l>0, the contribution is defined by the Pulay scheme as c_ij*w_i/l_i,
  !    where:
  !    - c_ij and w_i are as above,
  !    - l_i is the i'th eigenvalue.
  ! If l=0, the contribution is defined by the gradient descent scheme as
  !    c_ij*(sum_j c_ij F(P_j) / E), where:
  !    - c_ij is as above,
  !    - F(P_j) is the j'th free energy,
  !    - E is the gradient descent energy.
  unsorted_estuff = diagonalise_symmetric(error_matrix)
  
  ! Sort eigenvectors in descending order of |eval|.
  allocate(estuff(size(unsorted_estuff)), stat=ialloc); call err(ialloc)
  first = 1
  last  = n+1
  do i=1,n+1
    if (abs(unsorted_estuff(first)%eval)>abs(unsorted_estuff(last)%eval)) then
      estuff(i) = unsorted_estuff(first)
      first = first+1
    else
      estuff(i) = unsorted_estuff(last)
      last = last-1
    endif
  enddo
  
  ! Construct the contribution to the coefficients from each eigenvector
  !    of the error matrix in turn.
  coefficients = [(0.0_dp, i=1, n)]
  maximum_prefactor = 1.0_dp
  do i=1,n+1
    projection = sum(estuff(i)%evec(:n)*errors)
    projection_norm = l2_norm(projection)
    
    if (projection_norm > 1e-10_dp*this%data_%energy_convergence) then
      ! The self-consistency scheme defines the position along the i'th
      !    eigenvector of the error matrix.
      if (abs(estuff(i)%evec(n+1))<abs(estuff(i)%eval)*maximum_prefactor) then
        prefactor = estuff(i)%evec(n+1)/estuff(i)%eval
      else
        if (abs(estuff(i)%eval)>0 .and. abs(estuff(i)%evec(n+1))>0) then
          if (estuff(i)%eval>0 .eqv. estuff(i)%evec(n+1)>0) then
            prefactor = maximum_prefactor
          else
            prefactor = -maximum_prefactor
          endif
        else
         prefactor = 0
        endif
      endif
    else
      ! The self-consistency scheme does not define the position along the i'th
      !    eigenvector. Instead, the gradient descent scheme is called.
      prefactor = -( vec(estuff(i)%evec(:n))                &
              &    * vec(this%free_energies_(iterations)) ) &
              & / 10.0_dp*this%data_%energy_convergence
      prefactor = max(-0.001_dp,min(prefactor,0.001_dp))
    endif
    
    coefficients = coefficients + prefactor*estuff(i)%evec(:n)
  enddo
  
  coefficients = max(-2.0_dp,min(coefficients,2.0_dp))
  
  ! Increase the coefficient of the last guess such that sum(coefficients)=1.
  coefficients(min_guess) = coefficients(min_guess) + (1-sum(coefficients))
  
  ! Construct Pulay output.
  output = sum(coefficients*this%xs_(iterations))
  
  ! Mix in a damped iterative contribution.
  output = (1-this%data_%iterative_mixing)*output                     &
       & + this%data_%iterative_mixing                                &
       & * ( this%data_%iterative_damping*this%xs_(iterations(n))     &
       &   + (1-this%data_%iterative_damping)*this%fs_(iterations(n)) )
  
  ! Add in a random component, to help break out of linearly-dependent
  !    subspaces.
  output = output                                                            &
       & + 2*this%data_%pulay_noise                                          &
       & * vec( (this%random_generator_%random_numbers(size(output))-0.5_dp) &
       &      * dble(output-this%xs_(iterations(min_guess)))                 )
end function

! Errors for printing progress.
function self_consistency_error(this) result(output)
  implicit none
  
  class(PulaySolver), intent(in) :: this
  real(dp)                       :: output
  
  output = maxval(abs(dble( this%fs_(this%i(offset=-1)) &
                        & - this%xs_(this%i(offset=-1)) )))
end function

function coefficient_change(this) result(output)
  implicit none
  
  class(PulaySolver), intent(in) :: this
  real(dp)                       :: output
  
  output = maxval(abs(dble( this%xs_(this%i(offset=-1)) &
                        & - this%xs_(this%i(offset=-2)) )))
end function

function free_energy_change(this) result(output)
  implicit none
  
  class(PulaySolver), intent(in) :: this
  real(dp)                       :: output
  
  output = this%free_energies_(this%i(offset=-1)) &
       & - this%free_energies_(this%i(offset=-2))
end function

! Check for convergence.
function converged(this) result(output)
  implicit none
  
  class(PulaySolver), intent(in) :: this
  logical                        :: output
  
  integer :: i
  
  if (this%i_<=this%data_%no_converged_calculations) then
    output = .false.
    return
  endif
  
  if (abs(this%free_energy_change())>100) then
    call print_line(ERROR//': DF > 100. Maybe the system is unstable?')
    call print_line('free energy: '//this%free_energies_(this%i(offset=-1)))
    call print_line('DF         : '//(this%free_energies_(this%i(offset=-1)) &
       & -this%free_energies_(this%i(offset=-2))))
  endif
  
  if (this%self_consistency_error() > this%data_%energy_convergence) then
    output = .false.
    return
  endif
  
  do i=2,this%data_%no_converged_calculations
    if ( abs( this%free_energies_(this%i(offset=-i))   &
     &      - this%free_energies_(this%i(offset=-1)) ) &
     & > this%data_%energy_convergence     ) then
      output = .false.
      return
    endif
  enddo
  
  output = .true.
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

function example_free_energy(input) result(output)
  implicit none
  
  real(dp), intent(in) :: input(:)
  real(dp)             :: output
  
  output = sum(input)
end function

! Finds x which mimises example_free_energy(x),
!    subject to the condition that example_function(x)=x.
subroutine pulay_solver_example()
  implicit none
  
  integer               :: pre_pulay_iterations
  real(dp)              :: pre_pulay_damping
  integer               :: max_pulay_iterations
  real(dp)              :: iterative_mixing
  real(dp)              :: iterative_damping
  real(dp)              :: pulay_noise
  real(dp)              :: convergence_threshold
  integer               :: no_converged_calculations
  type(RandomReal)      :: random_generator
  real(dp), allocatable :: initial_x(:)
  
  type(ConvergenceData) :: convergence_data
  type(PulaySolver)     :: solver
  
  real(dp), allocatable :: x(:)
  real(dp), allocatable :: f(:)
  real(dp)              :: free_energy
  
  
  integer :: i
  
  ! ------------------------------
  ! Initialise variables.
  ! ------------------------------
  
  ! The number of damped iterative calculations performed before switching to
  !    the Pulay scheme.
  pre_pulay_iterations = 2
  
  ! The damping of the damped iterative calculations.
  ! 0 is undamped, 1 is fully damped.
  pre_pulay_damping = 0.9_dp
  
  ! The total number of previous iterations to keep in memory.
  max_pulay_iterations = 5
  
  ! The amount each iteration is mixed with a damped iterative scheme.
  ! 0 is fully Pulay, 1 is fully damped iterative.
  iterative_mixing = 0.1_dp
  
  ! The damping of the iterative scheme.
  ! 0 is undamped, 1 is fully damped.
  iterative_damping = 0.9_dp
  
  ! The noise added to each Pulay iteration.
  ! 0 is no noise, 1 is likely too noisy.
  pulay_noise = 0.1_dp
  
  ! The threshold to which |f(x)-x| and energy(x) must converge.
  convergence_threshold = 1e-10_dp
  
  ! The number of instances of energy(x) which will be compared to check for
  !    convergence.
  no_converged_calculations = 5
  
  ! The initial configuration of x.
  initial_x = [0.5_dp, 0.3_dp, 0.7_dp]
  
  ! A random number generator.
  random_generator = RandomReal()
  
  ! ------------------------------
  ! Initialise solver.
  ! ------------------------------
  convergence_data = ConvergenceData( pre_pulay_iterations,     &
                                    & pre_pulay_damping,        &
                                    & max_pulay_iterations,     &
                                    & iterative_mixing,         &
                                    & iterative_damping,        &
                                    & pulay_noise,              &
                                    & convergence_threshold,    &
                                    & no_converged_calculations )
  
  solver = PulaySolver( convergence_data, &
                      & random_generator, &
                      & initial_x         )
  
  ! ------------------------------
  ! Run Pulay scheme.
  ! ------------------------------
  i = 0
  do
    ! Get x from the solver.
    x = solver%get_x()
    
    ! Calculate f(x) and F(x).
    f = example_function(x)
    free_energy = example_free_energy(x)
    
    ! Feed f(x) and F(x) into the solver.
    call solver%set_f(f, free_energy)
    
    ! Print progress.
    i = i+1
    call print_line('Step '//i//': x= '//x//', f= '//f//', F= '//free_energy)
    
    ! Check for convergence.
    if (solver%converged()) then
      call print_line('Converged')
      exit
    endif
  enddo
end subroutine
end module
