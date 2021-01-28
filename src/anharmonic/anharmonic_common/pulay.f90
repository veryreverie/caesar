! ======================================================================
! Performs a Pulay self-consistency solver, to minimise F(x)
!    subject to f(x)=x.
! ======================================================================
! See the example function at the bottom of this file for how to use.
module caesar_pulay_module
  use caesar_precision_module
  use caesar_abstract_module
  use caesar_io_module
  use caesar_random_module
  
  use caesar_algebra_utils_module
  use caesar_linear_algebra_module
  use caesar_hermitian_eigenstuff_module
  use caesar_qr_decomposition_module
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
    real(dp),         allocatable, private :: error_norms_(:)
    real(dp),         allocatable, private :: free_energies_(:)
    type(RealVector), allocatable, private :: gradients_(:)
    logical,          allocatable, private :: gradients_set_(:)
    
    integer, allocatable, private :: iterations_(:)
    
    integer, private :: last_best_guess_
    integer, private :: best_guess_
    integer, private :: iterations_since_best_guess_
    
    integer, private :: i_
  contains
    ! Get x_i, from which f(x_i) and F(x_i) should be calculated.
    procedure, public :: calculate_x
    procedure, public :: get_x
    
    ! Set f(x_i) and F(x_i).
    procedure, public :: set_f
    
    ! Set the free energy gradient, d(F(x_i))/d(x_i).
    procedure, public :: gradient_requested
    procedure, public :: set_gradient
    
    ! The Pulay scheme itself.
    procedure, private :: mixed_pulay_
    procedure, private :: pulay_
    
    ! Errors for printing progress.
    procedure, private :: self_consistency_error
    procedure, private :: coefficient_change
    procedure, private :: free_energy_change
    
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
          & this%error_norms_(array_size),   &
          & this%free_energies_(array_size), &
          & this%gradients_(array_size),     &
          & this%gradients_set_(array_size), &
          & stat=ialloc); call err(ialloc)
  this%xs_(1) = vec(initial_x)
  this%gradients_set_ = .false.
  
  this%i_ = 0
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
    elseif (size(this%iterations_)+offset<1) then
      call print_line(CODE_ERROR//': offset before start of iterations.')
      call err()
    endif
  endif
  
  if (present(offset)) then
    output = this%iterations_(size(this%iterations_)+offset)
  else
    output = this%iterations_(size(this%iterations_))
  endif
end function

! Get x_i
function get_x(this) result(output)
  implicit none
  
  class(PulaySolver), intent(in) :: this
  real(dp), allocatable          :: output(:)
  
  output = dble(this%xs_(this%i()))
end function

! Set f(x_i) and F(x_i).
subroutine set_f(this,f,free_energy,print_progress)
  implicit none
  
  class(PulaySolver), intent(inout)        :: this
  real(dp),           intent(in)           :: f(:)
  real(dp),           intent(in)           :: free_energy
  logical,            intent(in), optional :: print_progress
  
  integer :: i
  
  i = this%i()
  
  ! Record f(x_i) and F(x_i).
  this%fs_(i) = vec(f)
  this%error_norms_(i) = l2_norm(this%fs_(i)-this%xs_(i))
  this%free_energies_(i) = free_energy
  
  if (this%i_==1) then
    this%last_best_guess_ = i
    this%best_guess_ = i
    this%iterations_since_best_guess_ = 0
  elseif (this%error_norms_(i)<this%error_norms_(this%best_guess_)) then
    this%last_best_guess_ = this%best_guess_
    this%best_guess_ = i
    this%iterations_since_best_guess_ = 0
  else
    this%iterations_since_best_guess_ = this%iterations_since_best_guess_ + 1
  endif
  
  if (set_default(print_progress, .false.) .and. this%i_>1) then
      call print_line('Self-consistency step '//this%i_//'.')
      call print_line( 'Self-consistency error : '   // &
                     & this%self_consistency_error() // &
                     & ' (Ha)'                          )
      call print_line( 'Change in coefficients : ' // &
                     & this%coefficient_change()   // &
                     & ' (Ha)'                        )
      call print_line( 'Free energy            : '   // &
                     & this%free_energies_(this%i()) // &
                     & ' (Ha)'                          )
      call print_line( 'Change in free energy  : ' // &
                     & this%free_energy_change()   // &
                     & ' (Ha)'                        )
  endif
end subroutine

! Check if d(F(x_i))/d(x_i) is requested.
! This is simply an optimisation; calculating the gradient is expensive,
!    so if it is not needed it should not be calculated.
function gradient_requested(this) result(output)
  implicit none
  
  class(PulaySolver), intent(in) :: this
  logical                        :: output
  
  ! Only request the gradient if the latest iteration has the smallest
  !    error norm.
  output = this%best_guess_==this%i()
end function

! Set d(F(x_i))/d(x_i).
subroutine set_gradient(this,gradient)
  implicit none
  
  class(PulaySolver), intent(inout) :: this
  real(dp),           intent(in)    :: gradient(:)
  
  this%gradients_(this%i()) = vec(gradient)
  this%gradients_set_(this%i()) = .true.
end subroutine

! Calculates x_{i+1}.
subroutine calculate_x(this)
  implicit none
  
  class(PulaySolver), intent(inout) :: this
  
  type(RealVector) :: x
  
  integer :: iteration
  
  integer :: i
  
  ! If this is the first iteration, nothing needs doing.
  if (this%i_==0) then
    this%i_ = 1
    this%iterations_ = [1]
    return
  endif
  
  ! Construct x_{i+1}.
  ! This uses a damped iterative scheme if i<pre_pulay_iterations, and a Pulay
  !    scheme otherwise.
  ! The Pulay scheme may fail, in which case i is reset to 1, and the scheme
  !    re-starts using the latest values of x and f(x).
  if (this%i_<=this%data_%pre_pulay_iterations) then
    ! The damped iterative scheme.
    ! x_{i+1} = ax_i + bf(x_i), where a=pre_pulay_damping and b=1-a.
    x = this%data_%pre_pulay_damping     &
    & * this%xs_(this%i())               &
    & + (1-this%data_%pre_pulay_damping) &
    & * this%fs_(this%i())
  else
    ! Calculate the next iteration using the Pulay scheme.
    x = this%mixed_pulay_()
  endif
  
  ! Bound the input at zero if required.
  ! N.B. this is done by limiting the input to be at least half the previous
  !    output, to avoid sharp changes.
  if (this%bound_at_zero_) then
    x = vec(max(dble(x), dble(this%fs_(this%i()))/2.0_dp))
  endif
  
  ! Update i.
  this%i_ = this%i_+1
  
  ! Update iterations_.
  if (this%i_<=size(this%xs_)) then
    ! If there are less iterations than max_pulay_iterations,
    !    append the iteration to the list.
    this%iterations_ = [this%iterations_, this%i_]
  else
    ! If max_pulay_iterations has been reached, find the iteration
    !    with the largest error, and replace it.
    ! The last two iterations are never replaced, to prevent the Pulay scheme
    !    calculating f(x) for the same x values repeatedly.
    i = maxloc( this%error_norms_(                                &
              &    this%iterations_(:size(this%iterations_)-2) ), &
              & 1)
    iteration = this%iterations_(i)
    this%iterations_(i:size(this%iterations_)-1) = &
       & this%iterations_(i+1:size(this%iterations_))
    this%iterations_(size(this%iterations_)) = iteration
  endif
  
  this%xs_(this%i()) = x
  this%gradients_set_(this%i()) = .false.
end subroutine

function mixed_pulay_(this) result(output)
  implicit none
  
  class(PulaySolver), intent(in) :: this
  type(RealVector)               :: output
  
  real(dp) :: min_error
  real(dp) :: max_error
  
  integer :: min_guess
  integer :: max_guess
  
  type(RealVector) :: best_guess
  type(RealVector) :: pulay
  logical          :: downhill_calculated
  type(RealVector) :: downhill
  type(RealVector) :: random
  
  ! Check for over-convergence.
  if (any( this%error_norms_(this%iterations_)     &
       & < 1e-10_dp*this%data_%energy_convergence  )) then
    output = this%xs_(this%iterations_(size(this%iterations_)))
    return
  endif
  
  associate(errors=>this%error_norms_(this%iterations_))
    min_guess = minloc(errors, 1)
    min_error = errors(min_guess)
    max_guess = maxloc(errors, 1)
    max_error = errors(max_guess)
  end associate
  
  best_guess = this%xs_(this%best_guess_)
  pulay      = (1-this%data_%iterative_damping)*(this%pulay_()-best_guess)
  if (this%gradients_set_(this%best_guess_)) then
    associate(gradient=>this%gradients_(this%best_guess_))
      downhill_calculated = .true.
      downhill = - vec(max(-min_error,min(dble(gradient),min_error))) &
             & * 1e-3_dp
    end associate
  else
    downhill_calculated = .false.
    downhill = best_guess*0
  endif
  random = 2*this%data_%pulay_noise                                         &
       & * 10**(3*this%random_generator_%random_number())                   &
       & * vec( ( this%random_generator_%random_numbers(size(best_guess))   &
       &        - 0.5_dp                                                  ) &
       &      * min_error*abs(dble(best_guess))                             )
  
  if ( this%iterations_since_best_guess_==0 .or.                  &
     & this%random_generator_%random_number(0.0_dp,1.0_dp)<0.2_dp ) then
    call print_line(colour('pulay','red'))
    output = best_guess + pulay
  elseif ( this%iterations_since_best_guess_<=1 .and.                      &
         & this%random_generator_%random_number(0.0_dp,1.0_dp)<0.2_dp .or. &
         & ( downhill_calculated                  .and.                    &
         &   this%random_generator_%random_number(0.0_dp,1.0_dp)<0.8_dp ) ) then
    call print_line(colour('downhill','blue'))
    output = best_guess + downhill
  elseif (       modulo(this%iterations_since_best_guess_,2)==1 &
         & .and. downhill_calculated                            ) then
    call print_line(colour('downhill','magenta'))
    output = best_guess &
         & + downhill   &
         & * this%random_generator_%random_number(0.01_dp,1._dp,.true.)
  else
    call print_line(colour('mixed','green'))
    output = best_guess                                                 &
         & + pulay                                                      &
         & * this%random_generator_%random_number(0.01_dp,1._dp,.true.) &
         & + downhill                                                   &
         & * this%random_generator_%random_number(0.01_dp,1._dp,.true.) &
         & + random                                                     &
         & * this%random_generator_%random_number(0.01_dp,1._dp,.true.)
  endif
end function

! A Pulay scheme to find self-consistent solution to f(x)=x.
! If the Pulay scheme is ill-conditioned (i.e. there is a subspace of
!    self-consistent solutions) then a gradient-descent scheme will be run
!    within the self-consistent subspace.
function pulay_(this) result(output)
  implicit none
  
  class(PulaySolver), intent(in) :: this
  type(RealVector)               :: output
  
  type(RealVector), allocatable :: errors(:)
  type(RealVector), allocatable :: preconditioned_errors(:)
  real(dp),         allocatable :: error_matrix(:,:)
  
  real(dp) :: min_error
  real(dp) :: max_error
  
  integer :: min_guess
  integer :: max_guess
  
  real(dp) :: scaling
  
  real(dp) :: regularisation
  
  type(RealVector) :: x_projection
  real(dp)         :: x_projection_norm
  type(RealVector) :: f_projection
  type(RealVector) :: error_projection
  real(dp)         :: error_projection_norm
  
  type(RealVector) :: x_solution
  type(RealVector) :: f_solution
  
  real(dp) :: maximum_projection
  
  type(SymmetricEigenstuff), allocatable :: unsorted_estuff(:)
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  
  integer :: first_,last_
  
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
  
  m = size(this%xs_(this%iterations_(1)))
  n = size(this%iterations_)
  
  if (any([(size(this%xs_(this%iterations_(i))),i=1,n)]/=m)) then
    call print_line(ERROR//': Input vectors inconsistent.')
  elseif (any([(size(this%fs_(this%iterations_(i))),i=1,n)]/=m)) then
    call print_line(ERROR//': Output vectors inconsistent.')
  endif
  
  associate(errors=>this%error_norms_(this%iterations_))
    min_guess = minloc(errors, 1)
    min_error = errors(min_guess)
    max_guess = maxloc(errors, 1)
    max_error = errors(max_guess)
  end associate
  
  ! Construct errors, e_i = f(x_i)-x_i.
  errors = this%fs_(this%iterations_) - this%xs_(this%iterations_)
  
  ! Precondition errors.
  preconditioned_errors = precondition_errors( this%xs_(this%iterations_), &
                                             & errors                )
  
  ! Construct error matrix.
  ! ( e1.e1, e1.e2, ... , e1.en,  1 )
  ! ( e2.e1, e2.e2, ... , e2.en,  1 )
  ! (  ... ,  ... , ... ,  ... , ...)
  ! ( en.e1, en.e2, ... , en.en,  1 )
  ! (   1  ,   1  , ... ,   1  ,  0 )
  allocate(error_matrix(n+1,n+1), stat=ialloc); call err(ialloc)
  do i=1,n
    do j=1,n
      error_matrix(j,i) = preconditioned_errors(j) * preconditioned_errors(i)
    enddo
  enddo
  error_matrix(:n , n+1) = 1
  error_matrix(n+1, :n ) = 1
  error_matrix(n+1, n+1) = 0
  
  ! Regularise the error matrix.
  regularisation = max( min_error**2,                    &
                      & this%data_%energy_convergence**2 ) / 1e2_dp
  do i=1,n
    error_matrix(i,i) = error_matrix(i,i)+regularisation
  enddo
  
  ! Scale error_matrix(:n,:n) such that the problem is well-conditioned.
  ! This has no effect on the result beyond numerical error.
  if (max_guess>0.01_dp*this%data_%energy_convergence**2) then
    scaling = 1/max_error**2
  else
    scaling = 1.0_dp
  endif
  if (max_guess>0.01_dp*this%data_%energy_convergence**2) then
    error_matrix(:n,:n) = error_matrix(:n,:n) * scaling
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
  first_ = 1
  last_  = n+1
  do i=1,n+1
    if ( abs(unsorted_estuff(first_)%eval) &
     & > abs(unsorted_estuff(last_)%eval)  ) then
      estuff(i) = unsorted_estuff(first_)
      first_ = first_+1
    else
      estuff(i) = unsorted_estuff(last_)
      last_ = last_-1
    endif
  enddo
  
  ! Construct the contribution to the coefficients from each eigenvector
  !    of the error matrix in turn.
  maximum_projection = 1e4_dp                                        &
                   & * max(min_error, this%data_%energy_convergence) &
                   & * size(this%iterations_)
  x_solution = this%xs_(this%iterations_(min_guess))
  f_solution = this%fs_(this%iterations_(min_guess))
  do i=1,n+1
    x_projection = sum( estuff(i)%evec(:n)                  &
                      & * ( this%xs_(this%iterations_)              &
                      &   - this%xs_(this%iterations_(min_guess)) ) )
    x_projection_norm = l2_norm(x_projection)
    
    error_projection = sum(estuff(i)%evec(:n)*errors)
    error_projection_norm = l2_norm(error_projection)
    
    f_projection = sum( estuff(i)%evec(:n)                  &
                    & * ( this%fs_(this%iterations_)              &
                    &   - this%fs_(this%iterations_(min_guess)) ) )
    
    if (error_projection_norm > 1e-10_dp*x_projection_norm) then
      ! The self-consistency scheme defines the position along the i'th
      !    eigenvector of the error matrix.
      x_solution = x_solution &
               & + x_projection * estuff(i)%evec(n+1) / estuff(i)%eval
      f_solution = f_solution &
               & + f_projection * estuff(i)%evec(n+1) / estuff(i)%eval
    else
      ! The self-consistency scheme does not define the position along the i'th
      !    eigenvector. Instead, the gradient descent scheme is called.
      x_solution = x_solution                               &
               & - x_projection                             &
               & * ( vec(estuff(i)%evec(:n))                &
               &   * vec(this%free_energies_(this%iterations_)) ) &
               & / (10.0_dp*this%data_%energy_convergence)
    endif
  enddo
  
  output = (1-this%data_%iterative_mixing)*x_solution &
       & + this%data_%iterative_mixing*f_solution
end function

! The Pulay scheme assumes that the error vectors e and input vectors x are
!    linearly related by some tensor a such that e=a.x', where x' is the
!    extended vector x'=[x,l] for some Lagrange parameter l.
! This subroutine performs a linear fit for the tensor a in the subspace of
!    x for which it is fully specified, and then replaces e vectors with the
!    corresponding a.x' vectors.
! If the x' vectors are linearly independent,
!    this will not change the error vectors.
! N.B. a weighted fit is used, favouring small values of |e|, so the fit
!    focusses on the region close to the solution e=0.
function precondition_errors(inputs,errors) result(output)
  implicit none
  
  type(RealVector), intent(in)  :: inputs(:)
  type(RealVector), intent(in)  :: errors(:)
  type(RealVector), allocatable :: output(:)
  
  real(dp), allocatable :: error_norms(:)
  integer,  allocatable :: sort_key(:)
  real(dp), allocatable :: weights(:)
  
  type(RealVector), allocatable :: sorted_inputs(:)
  
  real(dp) :: offset
  
  real(dp), allocatable :: extended_input_matrix(:,:)
  
  type(RealQRDecomposition) :: qr
  
  integer                       :: reduced_dimension
  type(RealVector), allocatable :: extended_inputs(:)
  
  type(RealMatrix) :: xx
  type(RealMatrix) :: ex
  
  real(dp) :: regularisation
  
  type(RealMatrix) :: mapping
  
  integer :: i,ialloc
  
  ! Sort the inputs in ascending order of |error|.
  ! Also subtract the input corresponding to the minimum error
  !    from every input.
  error_norms = l2_norm(errors)
  sort_key = sort(error_norms)
  sorted_inputs = inputs(sort_key)-inputs(sort_key(1))
  
  ! Construct the fitting weights.
  weights = [ max((error_norms(sort_key(1))/error_norms(sort_key))**2.0_dp, &
            & 1e-100_dp)                                                    ]
  
  ! Construct a matrix where each column is an extended input vector,
  !    and the columns are sorted in ascending order of |e|.
  offset = maxval(l2_norm(sorted_inputs)*sqrt(weights))
  allocate( extended_input_matrix(size(inputs(1))+1,size(inputs)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(inputs)
    extended_input_matrix(:,i) = [dble(sorted_inputs(i)), offset]
  enddo
  
  ! Use a QR decomposition to express the input vectors in a minimal basis.
  ! Extract the (sorted) extended inputs in the reduced basis from the
  !    QR decomposition.
  qr = qr_decomposition(extended_input_matrix)
  reduced_dimension = min(size(inputs), size(inputs(1))+1)
  extended_inputs = [(vec(qr%r(:reduced_dimension,i)), i=1, size(inputs))]
  
  ! Construct sum_i w_i x_i^x_i and sum_i w_i e_i^x_i,
  !    where ^ is the outer product.
  xx = sum(weights * outer_product(extended_inputs,extended_inputs))
  ex = sum(weights * outer_product(errors(sort_key),extended_inputs))
  
  ! Regularise xx.
  ! Ideally, the regularisation would be a small fraction of the offset,
  !    but this is sometimes too small for numerical stability.
  regularisation = max(                                             &
     & offset**2/1e4_dp,                                            &
     & maxval([(abs(qr%r(i,i)),i=1,reduced_dimension)],1)**2/1e9_dp )
  
  mapping = ex*regularised_invert(xx, regularisation)
  
  ! Construct the linearised errors (unsorting the inputs in the process).
  allocate(output(size(inputs)), stat=ialloc); call err(ialloc)
  do i=1,size(inputs)
    output(sort_key(i)) = mapping * extended_inputs(i)
  enddo
end function

! Adaptively regularises the matrix and inverts it.
! Replaces the eigenvalues l_i with
!    l_i + r*e^(-l_i/r),
!    where r is the regularisation.
! This regularises the small eigenvalues without overly affecting the large
!    eigenvalues.
function regularised_invert(input,regularisation) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  real(dp),         intent(in) :: regularisation
  type(RealMatrix)             :: output
  
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  
  integer :: i
  
  estuff = diagonalise_symmetric(input)
  
  output = dblemat(zeroes(size(input,1),size(input,2)))
  
  do i=1,size(estuff)
    output = output                                                 &
         & + outer_product(vec(estuff(i)%evec),vec(estuff(i)%evec)) &
         & / ( regularisation*exp(-estuff(i)%eval/regularisation)   &
         &   + estuff(i)%eval                                       )
  enddo
end function

! Errors for printing progress.
function self_consistency_error(this) result(output)
  implicit none
  
  class(PulaySolver), intent(in) :: this
  real(dp)                       :: output
  
  output = maxval(abs(dble( this%fs_(this%i()) &
                        & - this%xs_(this%i()) )))
end function

function coefficient_change(this) result(output)
  implicit none
  
  class(PulaySolver), intent(in) :: this
  real(dp)                       :: output
  
  integer :: i,j
  
  i = this%i()
  if (this%best_guess_/=i) then
    j = this%best_guess_
  else
    j = this%last_best_guess_
  endif
  
  output = maxval(abs(dble(this%xs_(i) - this%xs_(j))))
end function

function free_energy_change(this) result(output)
  implicit none
  
  class(PulaySolver), intent(in) :: this
  real(dp)                       :: output
  
  integer :: i,j
  
  i = this%i()
  if (this%best_guess_/=i) then
    j = this%best_guess_
  else
    j = this%last_best_guess_
  endif
  
  output = this%free_energies_(i) &
       & - this%free_energies_(j)
end function

! Check for convergence.
function converged(this) result(output)
  implicit none
  
  class(PulaySolver), intent(in) :: this
  logical                        :: output
  
  integer :: i
  
  if (size(this%xs_(this%i()))==0) then
    output = .true.
    return
  endif
  
  if (this%i_<=this%data_%no_converged_calculations) then
    output = .false.
    return
  endif
  
  if (abs(this%free_energy_change())>100) then
    call print_line(WARNING//': Change in free energy > 100 (Ha). &
       &Maybe the system is unstable?')
    call print_line('free energy: '//this%free_energies_(this%i(-1)))
    call print_line('DF         : '//(this%free_energies_(this%i(-1)) &
       & -this%free_energies_(this%i(-2))))
    if (abs(this%free_energy_change())>1e4) then
      call err()
    endif
  endif
  
  if (this%self_consistency_error() > this%data_%energy_convergence) then
    output = .false.
    return
  endif
  
  do i=2,this%data_%no_converged_calculations
    if ( abs( this%free_energies_(this%i(-i))   &
     &      - this%free_energies_(this%i(-1)) ) &
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
    ! Calculate and then return x from the solver.
    call solver%calculate_x()
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
