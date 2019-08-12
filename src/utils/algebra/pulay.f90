! ======================================================================
! Performs a Pulay self-consistency solver, to minimise F(x)
!    subject to f(x)=x.
! ======================================================================
! See the example function at the bottom of this file for how to use.
module pulay_module
  use precision_module
  use abstract_module
  use io_module
  use linear_algebra_module
  use hermitian_eigenstuff_module
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
    
    ! Check if the Pulay scheme has converged.
    procedure, public :: converged
    
    procedure, private :: i
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
  
  if (pre_pulay_iterations<2) then
    call print_line(ERROR//': pre_pulay_iterations must be at least 2.')
    call quit()
  elseif (max_pulay_iterations<2) then
    call print_line(ERROR//': max_pulay_iterations must be at least 2.')
    call quit()
  endif
  
  this%pre_pulay_iterations_ = pre_pulay_iterations
  this%pre_pulay_damping_ = pre_pulay_damping
  this%max_pulay_iterations_ = max_pulay_iterations
  if (present(bound_at_zero)) then
    this%bound_at_zero_ = bound_at_zero
  else
    this%bound_at_zero_ = .false.
  endif
  
  allocate( this%xs_(max_pulay_iterations),            &
          & this%fs_(max_pulay_iterations),            &
          & this%free_energies_(max_pulay_iterations), &
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
      this%xs_(this%i()) = pulay(this%xs_(:this%i_-1), this%fs_(:this%i_-1))
    else
      this%xs_(this%i()) = pulay(this%xs_, this%fs_)
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
  
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  
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
  ! Each eigenvector has length n+1, where the first n components are
  !    coefficients, and the final compontent is a weight.
  ! The contribution of each eigenvector to the output is c*w/l,
  !    where c is the vector of coefficients,
  !    w is the weight,
  !    and l is the eigenvalue.
  estuff = diagonalise_symmetric(error_matrix)
  
  coefficients = [(0.0_dp, i=1, n)]
  do i=1,n+1
    if (abs(estuff(i)%evec(n+1))<abs(estuff(i)%eval)*100) then
      coefficients = coefficients        &
                 & + estuff(i)%evec(:n)  &
                 & * estuff(i)%evec(n+1) &
                 & / estuff(i)%eval
    endif
  enddo
  
  ! Construct output frequencies as x_n = sum_{i=1}^n a_i x_i.
  output = sum(coefficients(:n)*(0.99_dp*xs+0.01_dp*fs))
end function

! Check for convergence.
function converged(this,energy_convergence,no_converged_calculations) &
   & result(output)
  implicit none
  
  class(PulaySolver), intent(in) :: this
  real(dp),           intent(in) :: energy_convergence
  integer,            intent(in) :: no_converged_calculations
  logical                        :: output
  
  integer :: i
  
  if (this%i_<=no_converged_calculations) then
    output = .false.
    return
  endif
  
  do i=2,no_converged_calculations
    if ( abs( this%free_energies_(this%i(offset=-i))   &
     &      - this%free_energies_(this%i(offset=-1)) ) &
     & > energy_convergence                            ) then
      output = .false.
      return
    endif
  enddo
  
  if (any( abs(dble( this%fs_(this%i(offset=-1))    &
       &           - this%xs_(this%i(offset=-1)) )) &
       & > energy_convergence                       )) then
    output = .false.
    return
  endif
  
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
  real(dp)              :: free_energy
  
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
    
    ! Calculate f(x) and F(x).
    f = example_function(x)
    free_energy = example_free_energy(x)
    
    ! Feed f(x) into the solver.
    call solver%set_f(f, free_energy)
    
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
