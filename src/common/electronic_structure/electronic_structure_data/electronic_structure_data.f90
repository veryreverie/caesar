! ======================================================================
! A class for holding calculated energy and forces.
! ======================================================================
module caesar_electronic_structure_data_module
  use caesar_utils_module
  
  use caesar_structure_module
  use caesar_normal_mode_module
  
  use caesar_linear_response_module
  implicit none
  
  private
  
  public :: LinearResponse
  
  public :: ElectronicStructure
  
  type, extends(Stringsable) :: ElectronicStructure
    real(dp),                            private :: energy_
    type(CartesianForce),   allocatable, private :: forces_
    type(CartesianHessian), allocatable, private :: hessian_
    type(RealMatrix),       allocatable, private :: stress_
    type(LinearResponse),   allocatable, private :: linear_response_
  contains
    procedure, public :: energy => energy_ElectronicStructure
    procedure, public :: has_forces => has_forces_ElectronicStructure
    procedure, public :: forces => forces_ElectronicStructure
    procedure, public :: has_hessian => has_hessian_ElectronicStructure
    procedure, public :: hessian => hessian_ElectronicStructure
    procedure, public :: has_stress => has_stress_ElectronicStructure
    procedure, public :: stress => stress_ElectronicStructure
    procedure, public :: has_linear_response => &
                       & has_linear_response_ElectronicStructure
    procedure, public :: linear_response => &
                       & linear_response_ElectronicStructure
    ! I/O.
    procedure, public :: read  => read_ElectronicStructure
    procedure, public :: write => write_ElectronicStructure
  end type
  
  interface ElectronicStructure
    ! Constructor.
    module function new_ElectronicStructure(energy,forces,hessian,stress, &
       & linear_response) result(this) 
      real(dp),               intent(in)           :: energy
      type(CartesianForce),   intent(in), optional :: forces
      type(CartesianHessian), intent(in), optional :: hessian
      type(RealMatrix),       intent(in), optional :: stress
      type(LinearResponse),   intent(in), optional :: linear_response
      type(ElectronicStructure)                    :: this
    end function
  end interface
  
  interface
    ! Getters.
    impure elemental module function energy_ElectronicStructure(this) &
       & result(output) 
      class(ElectronicStructure), intent(in) :: this
      real(dp)                               :: output
    end function
  end interface
  
  interface
    impure elemental module function has_forces_ElectronicStructure(this) &
       & result(output) 
      class(ElectronicStructure), intent(in) :: this
      logical                                :: output
    end function
  end interface
  
  interface
    impure elemental module function forces_ElectronicStructure(this) &
       & result(output) 
      class(ElectronicStructure), intent(in) :: this
      type(CartesianForce)                   :: output
    end function
  end interface
  
  interface
    impure elemental module function has_hessian_ElectronicStructure(this) &
       & result(output) 
      class(ElectronicStructure), intent(in) :: this
      logical                                :: output
    end function
  end interface
  
  interface
    impure elemental module function hessian_ElectronicStructure(this) &
       & result(output) 
      class(ElectronicStructure), intent(in) :: this
      type(CartesianHessian)                 :: output
    end function
  end interface
  
  interface
    impure elemental module function has_stress_ElectronicStructure(this) &
       & result(output) 
      class(ElectronicStructure), intent(in) :: this
      logical                                :: output
    end function
  end interface
  
  interface
    impure elemental module function stress_ElectronicStructure(this) &
       & result(output) 
      class(ElectronicStructure), intent(in) :: this
      type(RealMatrix)                       :: output
    end function
  end interface
  
  interface
    impure elemental module function has_linear_response_ElectronicStructure(this) result(output) 
      class(ElectronicStructure), intent(in) :: this
      logical                                :: output
    end function
  end interface
  
  interface
    impure elemental module function linear_response_ElectronicStructure(this) result(output) 
      class(ElectronicStructure), intent(in) :: this
      type(LinearResponse)                   :: output
    end function
  end interface
  
  interface
    ! I/O.
    module subroutine read_ElectronicStructure(this,input) 
      class(ElectronicStructure), intent(out) :: this
      type(String),               intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_ElectronicStructure(this) result(output) 
      class(ElectronicStructure), intent(in) :: this
      type(String), allocatable              :: output(:)
    end function
  end interface
  
  interface ElectronicStructure
    module function new_ElectronicStructure_Strings(input) result(this) 
      type(String), intent(in)  :: input(:)
      type(ElectronicStructure) :: this
    end function
  
    impure elemental module function new_ElectronicStructure_StringArray(input) result(this) 
      type(StringArray), intent(in) :: input
      type(ElectronicStructure)     :: this
    end function
  end interface
end module
