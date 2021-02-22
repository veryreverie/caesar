! ======================================================================
! Performs a bounded Newton-Raphson minimisation.
! ======================================================================
! See the example function at the bottom of this file for how to use.
module caesar_newton_raphson_module
  use caesar_foundations_module
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
    ! Construct and initialise the solver.
    module function new_NewtonRaphson(starting_value,finite_displacement, &
       & convergence_threshold,lower_bound,upper_bound) result(this) 
      real(dp), intent(in)           :: starting_value
      real(dp), intent(in)           :: finite_displacement
      real(dp), intent(in)           :: convergence_threshold
      real(dp), intent(in), optional :: lower_bound
      real(dp), intent(in), optional :: upper_bound
      type(NewtonRaphson)            :: this
    end function
  end interface
  
  interface
    ! Get the next inputs, x, for which f(x) should be calculated.
    module function get_inputs(this) result(output) 
      class(NewtonRaphson), intent(in) :: this
      real(dp), allocatable            :: output(:)
    end function
  end interface
  
  interface
    ! Set the next outputs, f(x).
    ! Runs the bounded Newton-Raphson scheme to calculate the next inputs.
    module subroutine set_outputs(this,input) 
      class(NewtonRaphson), intent(inout) :: this
      real(dp),             intent(in)    :: input(:)
    end subroutine
  end interface
  
  interface
    ! Checks whether or not the solver has converged.
    module function converged(this) result(output) 
      class(NewtonRaphson), intent(in) :: this
      logical                          :: output
    end function
  end interface
  
  interface
    ! Returns the solution.
    module function solution(this) result(output) 
      class(NewtonRaphson), intent(in) :: this
      real(dp)                         :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! An example demonstrating the use of the NewtonRaphson class.
    ! ----------------------------------------------------------------------
    impure elemental module function example_potential(input) result(output) 
      real(dp), intent(in) :: input
      real(dp)             :: output
    end function
  end interface
  
  interface
    ! Finds x such that example_potential(x) is minimised.
    module subroutine newton_raphson_example() 
    end subroutine
  end interface
end module
