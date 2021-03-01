! ======================================================================
! A force along a single real mode, in mass-weighted co-ordinates.
! ======================================================================
module caesar_real_single_mode_force_module
  use caesar_utils_module
  
  use caesar_structure_module
  
  use caesar_cartesian_force_module
  use caesar_mass_weighted_force_module
  use caesar_real_mode_module
  implicit none
  
  private
  
  public :: RealSingleForce
  public :: MassWeightedForce
  public :: CartesianForce
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  public :: select_mode
  public :: select_modes
  
  type, extends(Stringable) :: RealSingleForce
    ! The id of the mode.
    integer :: id
    
    ! The magnitude of vector along the mode, in mass-weighted co-ordinates.
    real(dp) :: magnitude
  contains
    ! I/O.
    procedure, public :: read  => read_RealSingleForce
    procedure, public :: write => write_RealSingleForce
  end type
  
  interface RealSingleForce
    ! Constructors.
    impure elemental module function new_RealSingleForce(id,magnitude) &
       & result(this) 
      integer,  intent(in)  :: id
      real(dp), intent(in)  :: magnitude
      type(RealSingleForce) :: this
    end function
  
    impure elemental module function new_RealSingleForce_RealMode(mode, &
       & magnitude) result(this) 
      type(RealMode), intent(in) :: mode
      real(dp),       intent(in) :: magnitude
      type(RealSingleForce)      :: this
    end function
  end interface
  
  interface operator(*)
    ! ----------------------------------------------------------------------
    ! Arithmetic.
    ! ----------------------------------------------------------------------
    impure elemental module function multiply_real_RealSingleForce(this,that) &
       & result(output) 
      real(dp),              intent(in) :: this
      type(RealSingleForce), intent(in) :: that
      type(RealSingleForce)             :: output
    end function
  
    impure elemental module function multiply_RealSingleForce_real(this,that) &
       & result(output) 
      type(RealSingleForce), intent(in) :: this
      real(dp),              intent(in) :: that
      type(RealSingleForce)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_RealSingleForce_real(this,that) &
       & result(output) 
      type(RealSingleForce), intent(in) :: this
      real(dp),              intent(in) :: that
      type(RealSingleForce)             :: output
    end function
  end interface
  
  interface operator(+)
    impure elemental module function add_RealSingleForce_RealSingleForce(this,that) result(output) 
      type(RealSingleForce), intent(in) :: this
      type(RealSingleForce), intent(in) :: that
      type(RealSingleForce)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function negative_RealSingleForce(this) &
       & result(output) 
      type(RealSingleForce), intent(in) :: this
      type(RealSingleForce)             :: output
    end function
  
    impure elemental module function subtract_RealSingleForce_RealSingleForce(this,that) result(output) 
      type(RealSingleForce), intent(in) :: this
      type(RealSingleForce), intent(in) :: that
      type(RealSingleForce)             :: output
    end function
  end interface
  
  interface MassWeightedForce
    ! ----------------------------------------------------------------------
    ! Conversions to and from cartesian and mass-weighted co-ordinates.
    ! ----------------------------------------------------------------------
    ! Constructs the MassWeightedForce corresponding to this vector.
    impure elemental module function new_MassWeightedForce_RealSingleForce(this,mode,structure,qpoint) result(output) 
      class(RealSingleForce), intent(in) :: this
      type(RealMode),         intent(in) :: mode
      type(StructureData),    intent(in) :: structure
      type(QpointData),       intent(in) :: qpoint
      type(MassWeightedForce)            :: output
    end function
  end interface
  
  interface CartesianForce
    ! Constructs the CartesianForce corresponding to this vector.
    impure elemental module function new_CartesianForce_RealSingleForce(this,mode,structure,qpoint) result(output) 
      class(RealSingleForce), intent(in) :: this
      type(RealMode),         intent(in) :: mode
      type(StructureData),    intent(in) :: structure
      type(QpointData),       intent(in) :: qpoint
      type(CartesianForce)               :: output
    end function
  end interface
  
  interface RealSingleForce
    ! Constructs the vector corresponding to the component of a mass-weighted
    !    vector along this mode.
    impure elemental module function new_RealSingleForce_MassWeightedForce(mode,vector,structure,qpoint) result(this) 
      type(RealMode),          intent(in) :: mode
      type(MassWeightedForce), intent(in) :: vector
      type(StructureData),     intent(in) :: structure
      type(QpointData),        intent(in) :: qpoint
      type(RealSingleForce)               :: this
    end function
  end interface
  
  interface RealSingleForce
    ! Constructs the vector corresponding to the component of a cartesian
    !    vector along this mode.
    impure elemental module function new_RealSingleForce_CartesianForce(mode,vector,structure,qpoint) result(this) 
      type(RealMode),       intent(in) :: mode
      type(CartesianForce), intent(in) :: vector
      type(StructureData),  intent(in) :: structure
      type(QpointData),     intent(in) :: qpoint
      type(RealSingleForce)            :: this
    end function
  end interface
  
  interface select_mode
    ! ----------------------------------------------------------------------
    ! Select modes corresponding to a given force or forces.
    ! ----------------------------------------------------------------------
    module function select_mode_RealSingleForce(force,modes) result(output) 
      type(RealSingleForce), intent(in) :: force
      type(RealMode),        intent(in) :: modes(:)
      type(RealMode)                    :: output
    end function
  end interface
  
  interface select_modes
    module function select_modes_RealSingleForces(forces,modes) result(output) 
      type(RealSingleForce), intent(in) :: forces(:)
      type(RealMode),        intent(in) :: modes(:)
      type(RealMode), allocatable       :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_RealSingleForce(this,input) 
      class(RealSingleForce), intent(out) :: this
      type(String),           intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_RealSingleForce(this) result(output) 
      class(RealSingleForce), intent(in) :: this
      type(String)                       :: output
    end function
  end interface
  
  interface RealSingleForce
    impure elemental module function new_RealSingleForce_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(RealSingleForce)    :: this
    end function
  end interface
end module
