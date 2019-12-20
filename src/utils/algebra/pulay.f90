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
    real(dp) :: energy_convergence
    integer  :: no_converged_calculations
  end type
  
  interface ConvergenceData
    module procedure new_ConvergenceData
  end interface
  
  type, extends(NoDefaultConstructor) :: PulaySolver
    ! Starting parameters.
    integer,          private :: pre_pulay_iterations_
    real(dp),         private :: pre_pulay_damping_
    integer,          private :: max_pulay_iterations_
    real(dp),         private :: convergence_threshold_
    integer,          private :: no_converged_calculations_
    type(RandomReal), private :: random_generator_
    logical,          private :: bound_at_zero_
    
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
    
    ! Errors for printing progress.
    procedure, public :: self_consistency_error
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
   & pre_pulay_damping,max_pulay_iterations,energy_convergence,     &
   & no_converged_calculations) result(this)
  implicit none
  
  integer,  intent(in)  :: pre_pulay_iterations
  real(dp), intent(in)  :: pre_pulay_damping
  integer,  intent(in)  :: max_pulay_iterations
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
  
  integer :: ialloc
  
  this%pre_pulay_iterations_ = convergence_data%pre_pulay_iterations
  this%pre_pulay_damping_ = convergence_data%pre_pulay_damping
  this%max_pulay_iterations_ = convergence_data%max_pulay_iterations
  this%convergence_threshold_ = convergence_data%energy_convergence
  this%no_converged_calculations_ = convergence_data%no_converged_calculations
  this%random_generator_ = random_generator
  if (present(bound_at_zero)) then
    this%bound_at_zero_ = bound_at_zero
  else
    this%bound_at_zero_ = .false.
  endif
  
  allocate( this%xs_(this%max_pulay_iterations_),            &
          & this%fs_(this%max_pulay_iterations_),            &
          & this%free_energies_(this%max_pulay_iterations_), &
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
    output = modulo(this%i_+offset-1,this%max_pulay_iterations_)+1
  else
    output = modulo(this%i_-1,this%max_pulay_iterations_)+1
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
  if (this%i_<=this%pre_pulay_iterations_) then
    ! The damped iterative scheme.
    ! x_{i+1} = ax_i + bf(x_i), where a=pre_pulay_damping and b=1-a.
    this%xs_(this%i()) = this%pre_pulay_damping_     &
                     & * this%xs_(this%i(offset=-1)) &
                     & + (1-this%pre_pulay_damping_) &
                     & * this%fs_(this%i(offset=-1))
  else
    ! The Pulay scheme.
    if (this%i_<=this%max_pulay_iterations_) then
      ! If i_ is less than max_pulay_iterations, only the first i_-1 iterations
      !    are used, as these are all that are available.
      this%xs_(this%i()) = pulay( this%xs_(:this%i_-1),            &
                                & this%fs_(:this%i_-1),            &
                                & this%free_energies_(:this%i_-1), &
                                & this%convergence_threshold_,     &
                                & this%i(-1),                      &
                                & this%random_generator_           )
    else
      ! Otherwise, the last max_pulay_iterations are used.
      ! N.B. the pulay scheme is invariant under re-ordering of iterations,
      !    so it does not matter that the latest entries may be in the middle
      !    of the arrays.
      this%xs_(this%i()) = pulay( this%xs_,                    &
                                & this%fs_,                    &
                                & this%free_energies_,         &
                                & this%convergence_threshold_, &
                                & this%i(-1),                  &
                                & this%random_generator_       )
    endif
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
function pulay(xs,fs,free_energies,convergence_threshold,last_guess, &
   & random_generator) result(output)
  implicit none
  
  type(RealVector), intent(in) :: xs(:)
  type(RealVector), intent(in) :: fs(:)
  real(dp),         intent(in) :: free_energies(:)
  real(dp),         intent(in) :: convergence_threshold
  integer,          intent(in) :: last_guess
  type(RandomReal), intent(in) :: random_generator
  type(RealVector)             :: output
  
  type(RealVector), allocatable :: errors(:)
  real(dp),         allocatable :: error_matrix(:,:)
  real(dp),         allocatable :: lagrange_vector(:)
  real(dp),         allocatable :: coefficients(:)
  
  type(RealVector) :: x_projection
  real(dp)         :: x_projection_norm
  
  type(RealVector) :: projection
  real(dp)         :: projection_norm
  
  type(RealVector), allocatable :: previous_vectors(:)
  
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
    if (abs(estuff(first)%eval)>abs(estuff(last)%eval)) then
      estuff(i) = unsorted_estuff(first)
      first = first+1
    else
      estuff(i) = unsorted_estuff(last)
      last = last-1
    endif
  enddo
  
  coefficients = [(0.0_dp, i=1, n)]
  allocate(previous_vectors(0), stat=ialloc); call err(ialloc)
  do i=1,n+1
    
    x_projection = sum(estuff(i)%evec(:n)*xs)
    do j=1,size(previous_vectors)
      x_projection = x_projection &
                 & - (x_projection*previous_vectors(j))*previous_vectors(j)
    enddo
    x_projection_norm = l2_norm(x_projection)
    
    projection = sum(estuff(i)%evec(:n)*errors)
    projection_norm = l2_norm(projection)
    
    if (x_projection_norm < 1e-10_dp) then
      ! The vector has no independent projection in the input space.
      ! It exists only because the number of input vectors is greater than
      !    the number of degrees of freedom.
      prefactor = 0
    elseif (projection_norm > 1e-10_dp) then
      ! The self-consistency scheme defines the position along the i'th
      !    eigenvector of the error matrix.
      prefactor = estuff(i)%evec(n+1)/estuff(i)%eval
      previous_vectors = [previous_vectors, x_projection/x_projection_norm]
      prefactor = max(-0.1_dp,min(prefactor,0.1_dp))
    else
      ! The self-consistency scheme does not define the position along the i'th
      !    eigenvector. Instead, the gradient descent scheme is called.
      prefactor = -(vec(estuff(i)%evec(:n))*vec(free_energies)) &
              & / 10.0_dp*convergence_threshold
      previous_vectors = [previous_vectors, x_projection/x_projection_norm]
      prefactor = max(-0.001_dp,min(prefactor,0.001_dp))
    endif
    
    coefficients = coefficients + prefactor*estuff(i)%evec(:n)
  enddo
  
  coefficients = max(-2.0_dp,min(coefficients,2.0_dp))
  
  ! Construct output frequencies as x_n = sum_{i=1}^n a_i x_i,
  !    plus a small component of a_i f_i.
  !output = 0.8_dp*xs(last_guess) + 0.2_dp*fs(last_guess)
  output = sum(coefficients*(0.8_dp*xs+0.2_dp*fs)) &
       & + (1-sum(coefficients))*(0.8_dp*xs(last_guess)+0.2_dp*fs(last_guess))
  
  ! Add in a small random component, to help break out of linearly-dependent
  !    subspaces.
  output = output                                                      &
       & + 2e-1_dp                                                     &
       & * vec( (random_generator%random_numbers(size(output))-0.5_dp) &
       &      * dble(output-xs(last_guess))                            )
end function

! Errors for printing progress.
function self_consistency_error(this) result(output)
  implicit none
  
  class(PulaySolver), intent(in) :: this
  real(dp)                       :: output
  
  output = maxval(abs(dble( this%fs_(this%i(offset=-1)) &
                        & - this%xs_(this%i(offset=-1)) )))
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
  
  if (this%i_<=this%no_converged_calculations_) then
    output = .false.
    return
  endif
  
  if (abs(this%free_energy_change())>100) then
    call print_line(ERROR//': DF > 100. Maybe the system is unstable?')
    call print_line('free energy: '//this%free_energies_(this%i(offset=-1)))
    call print_line('DF         : '//(this%free_energies_(this%i(offset=-1)) &
       & -this%free_energies_(this%i(offset=-2))))
  endif
  
  if (this%self_consistency_error()>this%convergence_threshold_) then
    output = .false.
    return
  endif
  
  do i=2,this%no_converged_calculations_
    if ( abs( this%free_energies_(this%i(offset=-i))   &
     &      - this%free_energies_(this%i(offset=-1)) ) &
     & > this%convergence_threshold_                   ) then
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
  pre_pulay_iterations    = 2
  
  ! The damping of the damped iterative calculations.
  ! 0 is undamped, 1 is fully damped.
  pre_pulay_damping       = 0.9_dp
  
  ! The total number of previous iterations to keep in memory.
  max_pulay_iterations    = 20
  
  ! The threshold to which |f(x)-x| and energy(x) must converge.
  convergence_threshold     = 1e-10_dp
  
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
