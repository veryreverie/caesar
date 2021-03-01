! ======================================================================
! Performs a Pulay self-consistency solver, to minimise F(x)
!    subject to f(x)=x.
! ======================================================================
! See the example function at the bottom of this file for how to use.
module caesar_pulay_module
  use caesar_foundations_module
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
  
  interface ConvergenceData
    impure elemental module function new_ConvergenceData(pre_pulay_iterations,pre_pulay_damping,max_pulay_iterations,iterative_mixing,iterative_damping,pulay_noise,energy_convergence,no_converged_calculations) result(this) 
      integer,  intent(in)  :: pre_pulay_iterations
      real(dp), intent(in)  :: pre_pulay_damping
      integer,  intent(in)  :: max_pulay_iterations
      real(dp), intent(in)  :: iterative_mixing
      real(dp), intent(in)  :: iterative_damping
      real(dp), intent(in)  :: pulay_noise
      real(dp), intent(in)  :: energy_convergence
      integer,  intent(in)  :: no_converged_calculations
      type(ConvergenceData) :: this
    end function
  end interface
  
  interface PulaySolver
    module function new_PulaySolver(convergence_data,random_generator, &
       & initial_x,bound_at_zero) result(this) 
      type(ConvergenceData), intent(in)           :: convergence_data
      type(RandomReal),      intent(in)           :: random_generator
      real(dp),              intent(in)           :: initial_x(:)
      logical,               intent(in), optional :: bound_at_zero
      type(PulaySolver)                           :: this
    end function
  end interface
  
  interface
    ! Private index module function.
    module function i(this,offset) result(output) 
      class(PulaySolver), intent(in)           :: this
      integer,            intent(in), optional :: offset
      integer                                  :: output
    end function
  end interface
  
  interface
    ! Get x_i
    module function get_x(this) result(output) 
      class(PulaySolver), intent(in) :: this
      real(dp), allocatable          :: output(:)
    end function
  end interface
  
  interface
    ! Set f(x_i) and F(x_i).
    module subroutine set_f(this,f,free_energy,print_progress) 
      class(PulaySolver), intent(inout)        :: this
      real(dp),           intent(in)           :: f(:)
      real(dp),           intent(in)           :: free_energy
      logical,            intent(in), optional :: print_progress
    end subroutine
  end interface
  
  interface
    ! Check if d(F(x_i))/d(x_i) is requested.
    ! This is simply an optimisation; calculating the gradient is expensive,
    !    so if it is not needed it should not be calculated.
    module function gradient_requested(this) result(output) 
      class(PulaySolver), intent(in) :: this
      logical                        :: output
    end function
  end interface
  
  interface
    ! Set d(F(x_i))/d(x_i).
    module subroutine set_gradient(this,gradient) 
      class(PulaySolver), intent(inout) :: this
      real(dp),           intent(in)    :: gradient(:)
    end subroutine
  end interface
  
  interface
    ! Calculates x_{i+1}.
    module subroutine calculate_x(this) 
      class(PulaySolver), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    module function mixed_pulay_(this) result(output) 
      class(PulaySolver), intent(in) :: this
      type(RealVector)               :: output
    end function
  end interface
  
  interface
    ! A Pulay scheme to find self-consistent solution to f(x)=x.
    ! If the Pulay scheme is ill-conditioned (i.e. there is a subspace of
    !    self-consistent solutions) then a gradient-descent scheme will be run
    !    within the self-consistent subspace.
    module function pulay_(this) result(output) 
      class(PulaySolver), intent(in) :: this
      type(RealVector)               :: output
    end function
  end interface
  
  interface
    ! The Pulay scheme assumes that the error vectors e and input vectors x are
    !    linearly related by some tensor a such that e=a.x', where x' is the
    !    extended vector x'=[x,l] for some Lagrange parameter l.
    ! This function performs a linear fit for the tensor a in the subspace of
    !    x for which it is fully specified, and then replaces e vectors with the
    !    corresponding a.x' vectors.
    ! If the x' vectors are linearly independent,
    !    this will not change the error vectors.
    ! N.B. a weighted fit is used, favouring small values of |e|, so the fit
    !    focusses on the region close to the solution e=0.
    module function precondition_errors(inputs,errors) result(output) 
      type(RealVector), intent(in)  :: inputs(:)
      type(RealVector), intent(in)  :: errors(:)
      type(RealVector), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! Adaptively regularises the matrix and inverts it.
    ! Replaces the eigenvalues l_i with
    !    l_i + r*e^(-l_i/r),
    !    where r is the regularisation.
    ! This regularises the small eigenvalues without overly affecting the large
    !    eigenvalues.
    module function regularised_invert(input,regularisation) result(output) 
      type(RealMatrix), intent(in) :: input
      real(dp),         intent(in) :: regularisation
      type(RealMatrix)             :: output
    end function
  end interface
  
  interface
    ! Errors for printing progress.
    module function self_consistency_error(this) result(output) 
      class(PulaySolver), intent(in) :: this
      real(dp)                       :: output
    end function
  end interface
  
  interface
    module function coefficient_change(this) result(output) 
      class(PulaySolver), intent(in) :: this
      real(dp)                       :: output
    end function
  end interface
  
  interface
    module function free_energy_change(this) result(output) 
      class(PulaySolver), intent(in) :: this
      real(dp)                       :: output
    end function
  end interface
  
  interface
    ! Check for convergence.
    module function converged(this) result(output) 
      class(PulaySolver), intent(in) :: this
      logical                        :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! An example demonstrating the use of the PulaySolver class.
    ! ----------------------------------------------------------------------
    module function example_function(input) result(output) 
      real(dp), intent(in)  :: input(:)
      real(dp), allocatable :: output(:)
    end function
  end interface
  
  interface
    module function example_free_energy(input) result(output) 
      real(dp), intent(in) :: input(:)
      real(dp)             :: output
    end function
  end interface
  
  interface
    ! Finds x which mimises example_free_energy(x),
    !    subject to the condition that example_function(x) =x.
    module subroutine pulay_solver_example() 
    end subroutine
  end interface
end module
