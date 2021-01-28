! ======================================================================
! Performs a bounded Newton-Raphson minimisation.
! ======================================================================
! See the example function at the bottom of this file for how to use.
module caesar_newton_raphson_module
  use caesar_precision_module
  use caesar_abstract_module
  use caesar_io_module
  implicit none
  
  private
  
  public :: NewtonRaphson
  public :: newton_raphson_example
  
  type, extends(NoDefaultConstructor) :: NewtonRaphson
    ! Starting parameters.
    real(dp), private :: finite_displacement_
    real(dp), private :: convergence_threshold_
    logical,  private :: has_lower_bound_
    real(dp), private :: lower_bound_
    logical,  private :: has_upper_bound_
    real(dp), private :: upper_bound_
    
    ! The current set of x and f(x).
    real(dp), allocatable, private :: inputs_(:)
    real(dp), allocatable, private :: outputs_(:)
    
    ! The guesses so far.
    logical,               private :: converged_
    real(dp), allocatable, private :: guesses_(:)
  contains
    procedure, public :: get_inputs
    procedure, public :: set_outputs
    procedure, public :: converged
    procedure, public :: solution
  end type
  
  interface NewtonRaphson
    module procedure new_NewtonRaphson
  end interface
contains

! Construct and initialise the solver.
function new_NewtonRaphson(starting_value,finite_displacement, &
   & convergence_threshold,lower_bound,upper_bound) result(this)
  implicit none
  
  real(dp), intent(in)           :: starting_value
  real(dp), intent(in)           :: finite_displacement
  real(dp), intent(in)           :: convergence_threshold
  real(dp), intent(in), optional :: lower_bound
  real(dp), intent(in), optional :: upper_bound
  type(NewtonRaphson)            :: this
  
  ! Check inputs.
  if (finite_displacement<=0) then
    call print_line(ERROR//': Newton Raphson solver must have a finite &
       &displacement greater than zero.')
    call err()
  endif
  
  if (convergence_threshold<=0) then
    call print_line(ERROR//': Newton Raphson solver must have a convergence &
       &threshold greater than zero.')
    call err()
  endif
  
  if (present(lower_bound)) then
    if (starting_value<lower_bound) then
      call print_line(ERROR//': starting value must be >= lower bound.')
      call err()
    endif
  endif
  
  if (present(upper_bound)) then
    if (starting_value>upper_bound) then
      call print_line(ERROR//': starting value must be <= upper bound.')
      call err()
    endif
  endif
  
  if (present(lower_bound) .and. present(upper_bound)) then
    if (upper_bound-lower_bound<=3*finite_displacement) then
      call print_line(ERROR//': upper_bound-lower_bound must be larger than &
         &3*finite_displacement.')
      call err()
    endif
    
    if (upper_bound-lower_bound<=3*convergence_threshold) then
      call print_line(ERROR//': upper_bound-lower_bound must be larger than &
         &3*convergence_threshold.')
      call err()
    endif
  endif
  
  ! Set starting parameters.
  this%finite_displacement_   = finite_displacement
  this%convergence_threshold_ = convergence_threshold
  if (present(lower_bound)) then
    this%has_lower_bound_ = .true.
    this%lower_bound_ = lower_bound
  else
    this%has_lower_bound_ = .false.
  endif
  if (present(upper_bound)) then
    this%has_upper_bound_ = .true.
    this%upper_bound_ = upper_bound
  else
    this%has_upper_bound_ = .false.
  endif
  
  ! Generate initial inputs.
  this%converged_ = .false.
  this%guesses_ = [starting_value]
  
  if (present(lower_bound)) then
    if (starting_value-finite_displacement<lower_bound) then
      this%inputs_ = [ starting_value,                        &
                     & starting_value + finite_displacement,  &
                     & starting_value + 2*finite_displacement ]
      return
    endif
  endif
  
  if (present(upper_bound)) then
    if (starting_value+finite_displacement>upper_bound) then
      this%inputs_ = [ starting_value - 2*finite_displacement, &
                     & starting_value - finite_displacement,   &
                     & starting_value                          ]
      return
    endif
  endif
  
  this%inputs_ = [ starting_value - finite_displacement, &
                 & starting_value,                       &
                 & starting_value + finite_displacement  ]
end function

! Get the next inputs, x, for which f(x) should be calculated.
function get_inputs(this) result(output)
  implicit none
  
  class(NewtonRaphson), intent(in) :: this
  real(dp), allocatable            :: output(:)
  
  output = this%inputs_
end function

! Set the next outputs, f(x).
! Runs the bounded Newton-Raphson scheme to calculate the next inputs.
subroutine set_outputs(this,input)
  implicit none
  
  class(NewtonRaphson), intent(inout) :: this
  real(dp),             intent(in)    :: input(:)
  
  real(dp) ::xm, x, xp
  real(dp) ::fm, f, fp
  
  real(dp) :: df_dx
  real(dp) :: d2f_dx2
  real(dp) :: lower
  real(dp) :: upper
  
  real(dp) :: guess
  
  if (size(input)/=size(this%inputs_)) then
    call print_line(ERROR//': Given outputs do not match inputs.')
    call err()
  endif
  
  this%outputs_ = input
  
  ! ------------------------------
  ! Calculate the next guess.
  ! ------------------------------
  ! Record x-d, x, x_d, f(x-d), f(x) and f(x+d).
  xm = this%inputs_(1)
  x  = this%inputs_(2)
  xp = this%inputs_(3)
  
  fm = input(1)
  f  = input(2)
  fp = input(3)
  
  ! Calculate the derivatives df/dx and d2f/dx2 at x.
  ! df/dx = [f(x+d)-f(x-d)]/(2d)
  df_dx = (fp-fm)/(2*this%finite_displacement_)
  ! d2f/dx2 = (fm+fp-2f)/d**2
  d2f_dx2 = (fm+fp-2*f)/this%finite_displacement_**2
  
  ! Calculate bounds on the next guess.
  ! If the upper and lower bounds exist, the next guess must lie at most
  !    2/3 of the way from the current guess to either bound.
  ! If a bound does not exist, then the guess shoud be at most twice as far
  !    from the other bound as the last guess, unless that distance is smaller
  !    than finite_displacement.
  if (this%has_lower_bound_ .and. this%has_upper_bound_) then
    lower = (2*this%lower_bound_+x)/3
    upper = (2*this%upper_bound_+x)/3
  elseif (this%has_lower_bound_) then
    lower = (2*this%lower_bound_+x)/3
    upper = x + 2*max( x-this%lower_bound_,      &
                     & this%finite_displacement_ )
  elseif (this%has_upper_bound_) then
    lower = x - 2*max( this%upper_bound_-x,      &
                     & this%finite_displacement_ )
    upper = (2*this%upper_bound_+x)/3
  else
    lower = x - 2*this%finite_displacement_
    upper = x + 2*this%finite_displacement_
  endif
  
  ! Calculate the next guess.
  ! By Newton-Raphson, the next guess g is calculated in terms of the guess x
  !    as
  !          df/dx          (fp-fm)/(2d)
  ! g = x - ------- = x - ---------------
  !         d2f/dx2       (fm+fp-2f)/d**2
  if (input(2)<input(1) .and. input(2)<input(3)) then
    ! f(x-d) and f(x+d) are both greater than f(x).
    ! The next guess lies between the two.
    this%has_lower_bound_ = .true.
    this%lower_bound_ = xm
    this%has_upper_bound_ = .true.
    this%upper_bound_ = xp
    
    guess = x - df_dx / d2f_dx2
  elseif (input(1)<input(3)) then
    ! f(x+d) is greater than f(x-d). The next guess is less than x.
    this%has_upper_bound_ = .true.
    this%upper_bound_ = x
    if ((.not. this%has_lower_bound_) .or. abs(df_dx)>d2f_dx2*(x-lower)) then
      ! The guess is the lower bound.
      guess = lower
    else
      ! The guess is greater than the lower bound.
      guess = x - df_dx/d2f_dx2
    endif
  elseif (input(1)>input(3)) then
    ! f(x+d) is less than f(x-d). The next guess is greater than x.
    this%has_lower_bound_ = .true.
    this%lower_bound_ = x
    if ((.not. this%has_upper_bound_) .or. abs(df_dx)>d2f_dx2*(upper-x)) then
      ! The guess is the upper bound.
      guess = upper
    else
      ! The guess is less than the upper bound.
      guess = x - df_dx/d2f_dx2
    endif
  else
    ! f(x-d) = f(x+d) <= f(x).
    ! Increase the finite displacement, and try again.
    this%finite_displacement_ = this%finite_displacement_ * 3
    if (this%has_lower_bound_ .and. this%has_upper_bound_) then
      guess = x
      if ( this%finite_displacement_               &
       & > (this%upper_bound_-this%lower_bound_)/3 ) then
        ! The whole region between the bounds is flat.
        this%lower_bound_ = x
        this%upper_bound_ = x
      endif
    elseif (this%has_lower_bound_) then
      ! There is currently no upper bound. Attempt to find one.
      guess = upper
    elseif (this%has_upper_bound_) then
      ! There is currently no lower bound. Attempt to find one.
      guess = lower
    else
      ! There are no bounds. The guess cannot be updated.
      guess = x
    endif
  endif
  
  ! ------------------------------
  ! Update the list of guesses, and whether or not the solver has converged.
  ! ------------------------------
  if (this%has_upper_bound_ .and. this%has_lower_bound_) then
    this%converged_ = (this%upper_bound_-this%lower_bound_) &
                  & < this%convergence_threshold_
  else
    this%converged_ = .false.
  endif
  
  this%guesses_ = [this%guesses_, guess]
  
  ! ------------------------------
  ! Calculate the next inputs from the guess, updating finite_displacement
  !    if need be.
  ! ------------------------------
  if (this%has_lower_bound_ .and. this%has_upper_bound_) then
    this%finite_displacement_ = min( this%finite_displacement_,              &
                                   & (this%upper_bound_-this%lower_bound_)/3 )
  endif
  
  if (this%has_lower_bound_) then
    if (guess-this%finite_displacement_<this%lower_bound_) then
      this%inputs_ = [ guess,                              &
                     & guess + this%finite_displacement_,  &
                     & guess + 2*this%finite_displacement_ ]
      if (this%inputs_(3)>x) then
        this%finite_displacement_ = this%finite_displacement_ / 2
      endif
      return
    endif
  endif
  
  if (this%has_upper_bound_) then
    if (guess+this%finite_displacement_>this%upper_bound_) then
      this%inputs_ = [ guess - 2*this%finite_displacement_, &
                     & guess - this%finite_displacement_,   &
                     & guess                                ]
      if (this%inputs_(1)<x) then
        this%finite_displacement_ = this%finite_displacement_ / 2
      endif
      return
    endif
  endif
  
  this%inputs_ = [ guess - this%finite_displacement_, &
                 & guess,                             &
                 & guess + this%finite_displacement_  ]
end subroutine

! Checks whether or not the solver has converged.
function converged(this) result(output)
  implicit none
  
  class(NewtonRaphson), intent(in) :: this
  logical                          :: output
  
  output = this%converged_
end function

! Returns the solution.
function solution(this) result(output)
  implicit none
  
  class(NewtonRaphson), intent(in) :: this
  real(dp)                         :: output
  
  output = this%guesses_(size(this%guesses_))
end function

! ----------------------------------------------------------------------
! An example demonstrating the use of the NewtonRaphson class.
! ----------------------------------------------------------------------
impure elemental function example_potential(input) result(output)
  implicit none
  
  real(dp), intent(in) :: input
  real(dp)             :: output
  
  ! A quartic polynomial with a minimum at input=2, and no other extrema.
  ! The large constant makes the naive solver ill-conditioned.
  output = 1e4 + (input-2)**2 * (input**2+4)
end function

! Finds x such that example_potential(x) is minimised.
subroutine newton_raphson_example()
  implicit none
  
  real(dp) :: starting_value
  real(dp) :: finite_displacement
  real(dp) :: convergence_threshold
  real(dp) :: lower_bound
  
  type(NewtonRaphson) :: solver
  
  real(dp), allocatable :: inputs(:)
  real(dp), allocatable :: outputs(:)
  
  integer :: i
  
  starting_value        = 0
  finite_displacement   = 1e-8
  convergence_threshold = 1e-8
  lower_bound           = 0
  
  ! Initialise solver.
  solver = NewtonRaphson( starting_value,           &
                        & finite_displacement,      &
                        & convergence_threshold,    &
                        & lower_bound = lower_bound )
  i = 0
  do
    ! Get the set of points, x, at which to calculate f(x).
    inputs = solver%get_inputs()
    
    ! Calculate f(x).
    outputs = example_potential(inputs)
    
    ! Hand f(x) back to the solver.
    call solver%set_outputs(outputs)
    
    ! Print progress.
    i = i+1
    call print_line('Step '//left_pad(i,'aaaaaa')//', &
                     &x= '//solver%solution()//', &
                     &f= '//example_potential(solver%solution()))
    
    ! Check for convergence.
    if (solver%converged()) then
      call print_line('Converged')
      exit
    endif
  enddo
end subroutine
end module
