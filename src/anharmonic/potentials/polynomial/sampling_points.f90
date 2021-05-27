! ======================================================================
! Sampling points for sampling a given set of basis functions.
! ======================================================================
module caesar_sampling_points_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_basis_function_module
  implicit none
  
  private
  
  public :: SamplingPoints
  public :: generate_sampling_points
  public :: size
  public :: construct_sample_matrix
  public :: make_sample_vector
  
  type, extends(Stringsable) :: SamplingPoints
    type(RealModeDisplacement), allocatable :: points(:)
  contains
    procedure, public :: read  => read_SamplingPoints
    procedure, public :: write => write_SamplingPoints
  end type
  
  interface SamplingPoints
    ! Constructor and size() function.
    module function new_SamplingPoints(points) result(this) 
      type(RealModeDisplacement), intent(in) :: points(:)
      type(SamplingPoints)                   :: this
    end function
  end interface
  
  interface size
    module function size_SamplingPoints(this) result(output) 
      type(SamplingPoints), intent(in) :: this
      integer                          :: output
    end function
  end interface
  
  interface
    ! Construct the matrix necessary for L2 fitting of basis function coefficients.
    ! The elements are the energies and forces at a given sampling point due to
    !    a given basis function.
    module function construct_sample_matrix(basis_functions,sampling_points, &
       & modes,energy_force_ratio,sample_weights) result(output) 
      type(BasisFunction),        intent(in)           :: basis_functions(:) 
      type(RealModeDisplacement), intent(in)           :: sampling_points(:)
      type(RealMode),             intent(in)           :: modes(:)
      real(dp),                   intent(in)           :: energy_force_ratio
      real(dp),                   intent(in), optional :: sample_weights(:)
      real(dp), allocatable                            :: output(:,:)
    end function
  end interface
  
  interface
    ! Convert an energy and force into a single vector which can be inserted
    !    into a matrix for passing to LAPACK.
    module function make_sample_vector(energy,force,modes, &
       & energy_force_ratio,weight) result(output) 
      real(dp),            intent(in)           :: energy
      type(RealModeForce), intent(in)           :: force
      type(RealMode),      intent(in)           :: modes(:)
      real(dp),            intent(in)           :: energy_force_ratio
      real(dp),            intent(in), optional :: weight
      real(dp), allocatable                     :: output(:)
    end function
  end interface
  
  interface
    ! Generates a set of sampling points for sampling a set of basis functions.
    module function generate_sampling_points(subspace_coupling,       &
       & basis_functions,potential_expansion_order,                   &
       & maximum_weighted_displacement,frequency_of_max_displacement, &
       & real_modes,energy_to_force_ratio) result(output) 
      type(SubspaceCoupling), intent(in) :: subspace_coupling
      type(BasisFunction),    intent(in) :: basis_functions(:) 
      integer,                intent(in) :: potential_expansion_order
      real(dp),               intent(in) :: maximum_weighted_displacement
      real(dp),               intent(in) :: frequency_of_max_displacement
      type(RealMode),         intent(in) :: real_modes(:)
      real(dp),               intent(in) :: energy_to_force_ratio
      type(SamplingPoints)               :: output
    end function
  end interface
  
  interface
    module function generate_basis_points(basis_function,         &
       & potential_expansion_order,maximum_weighted_displacement, &
       & frequency_of_max_displacement,modes) result(output) 
      type(BasisFunction), intent(in)         :: basis_function
      integer,             intent(in)         :: potential_expansion_order
      real(dp),            intent(in)         :: maximum_weighted_displacement
      real(dp),            intent(in)         :: frequency_of_max_displacement
      type(RealMode),      intent(in)         :: modes(:)
      type(RealModeDisplacement), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_SamplingPoints(this,input) 
      class(SamplingPoints), intent(out) :: this
      type(String),          intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_SamplingPoints(this) result(output) 
      class(SamplingPoints), intent(in) :: this
      type(String), allocatable         :: output(:)
    end function
  end interface
  
  interface SamplingPoints
    module function new_SamplingPoints_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(SamplingPoints)     :: this
    end function
  
    impure elemental module function new_SamplingPoints_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(SamplingPoints)          :: this
    end function
  end interface
end module
