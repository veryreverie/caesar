! ======================================================================
! A force in real mode co-ordinates.
! ======================================================================
module caesar_real_mode_force_module
  use caesar_utils_module
  
  use caesar_structure_module
  
  use caesar_cartesian_force_module
  use caesar_mass_weighted_force_module
  use caesar_real_mode_module
  use caesar_real_single_mode_force_module
  implicit none
  
  private
  
  public :: RealModeForce
  public :: size
  public :: MassWeightedForce
  public :: CartesianForce
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
  type, extends(Stringsable) :: RealModeForce
    type(RealSingleForce), allocatable :: vectors(:)
  contains
    ! The component of the force along a given mode.
    generic,   public  :: force =>  &
                        & force_id, &
                        & force_mode
    procedure, private :: force_id
    procedure, private :: force_mode
    
    ! I/O.
    procedure, public :: read  => read_RealModeForce
    procedure, public :: write => write_RealModeForce
  end type
  
  interface RealModeForce
    ! Constructors and size() module function.
    module function new_RealModeForce(forces) result(this) 
      type(RealSingleForce), intent(in) :: forces(:)
      type(RealModeForce)               :: this
    end function
  
    module function new_RealModeForce_RealModes(modes,forces) result(this) 
      type(RealMode), intent(in) :: modes(:)
      real(dp),       intent(in) :: forces(:)
      type(RealModeForce)        :: this
    end function
  end interface
  
  interface size
    module function size_RealModeForce(this) result(output) 
      type(RealModeForce), intent(in) :: this
      integer                         :: output
    end function
  end interface
  
  interface operator(*)
    ! ----------------------------------------------------------------------
    ! Arithmetic.
    ! ----------------------------------------------------------------------
    impure elemental module function multiply_real_RealModeForce(this,that) &
       & result(output) 
      real(dp),            intent(in) :: this
      type(RealModeForce), intent(in) :: that
      type(RealModeForce)             :: output
    end function
  
    impure elemental module function multiply_RealModeForce_real(this,that) &
       & result(output) 
      type(RealModeForce), intent(in) :: this
      real(dp),            intent(in) :: that
      type(RealModeForce)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_RealModeForce_real(this,that) &
       & result(output) 
      type(RealModeForce), intent(in) :: this
      real(dp),            intent(in) :: that
      type(RealModeForce)             :: output
    end function
  end interface
  
  interface operator(+)
    impure elemental module function add_RealModeForce_RealModeForce(this, &
       & that) result(output) 
      type(RealModeForce), intent(in) :: this
      type(RealModeForce), intent(in) :: that
      type(RealModeForce)             :: output
    end function
  end interface
  
  interface sum
    module function sum_RealModeForces(this) result(output) 
      type(RealModeForce), intent(in) :: this(:)
      type(RealModeForce)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function negative_RealModeForce(this) &
       & result(output) 
      type(RealModeForce), intent(in) :: this
      type(RealModeForce)             :: output
    end function
  
    impure elemental module function subtract_RealModeForce_RealModeForce(   this,that) result(output) 
      type(RealModeForce), intent(in) :: this
      type(RealModeForce), intent(in) :: that
      type(RealModeForce)             :: output
    end function
  end interface
  
  interface MassWeightedForce
    ! ----------------------------------------------------------------------
    ! Conversions to and from cartesian and mass-weighted co-ordinates.
    ! ----------------------------------------------------------------------
    ! Returns the force in mass-weighted co-ordinates.
    module function new_MassWeightedForce_RealModeForce(this,structure, &
       & modes,qpoints) result(output) 
      class(RealModeForce), intent(in) :: this
      type(StructureData),  intent(in) :: structure
      type(RealMode),       intent(in) :: modes(:)
      type(QpointData),     intent(in) :: qpoints(:)
      type(MassWeightedForce)          :: output
    end function
  end interface
  
  interface CartesianForce
    ! Returns the force in cartesian co-ordinates.
    module function new_CartesianForce_RealModeForce(this,structure,modes, &
       & qpoints) result(output) 
      class(RealModeForce), intent(in) :: this
      type(StructureData),  intent(in) :: structure
      type(RealMode),       intent(in) :: modes(:)
      type(QpointData),     intent(in) :: qpoints(:)
      type(CartesianForce)             :: output
    end function
  end interface
  
  interface RealModeForce
    ! Converts a MassWeightedForce to a RealModeForce.
    module function new_RealModeForce_MassWeightedForce(force,structure, &
       & modes,qpoints) result(this) 
      type(MassWeightedForce), intent(in) :: force
      type(StructureData),     intent(in) :: structure
      type(RealMode),          intent(in) :: modes(:)
      type(QpointData),        intent(in) :: qpoints(:)
      type(RealModeForce)                 :: this
    end function
  
    ! Converts a CartesianForce to a RealModeForce.
    module function new_RealModeForce_CartesianForce(force,structure,modes, &
       & qpoints) result(this) 
      type(CartesianForce), intent(in) :: force
      type(StructureData),  intent(in) :: structure
      type(RealMode),       intent(in) :: modes(:)
      type(QpointData),     intent(in) :: qpoints(:)
      type(RealModeForce)              :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! The displacement along a given mode.
    ! ----------------------------------------------------------------------
    impure elemental module function force_id(this,id) result(output) 
      class(RealModeForce), intent(in) :: this
      integer,              intent(in) :: id
      real(dp)                         :: output
    end function
  end interface
  
  interface
    impure elemental module function force_mode(this,mode) result(output) 
      class(RealModeForce), intent(in) :: this
      type(RealMode),       intent(in) :: mode
      real(dp)                         :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_RealModeForce(this,input) 
      class(RealModeForce), intent(out) :: this
      type(String),         intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_RealModeForce(this) result(output) 
      class(RealModeForce), intent(in) :: this
      type(String), allocatable        :: output(:)
    end function
  end interface
  
  interface RealModeForce
    module function new_RealModeForce_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(RealModeForce)      :: this
    end function
  
    impure elemental module function new_RealModeForce_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(RealModeForce)           :: this
    end function
  end interface
end module
