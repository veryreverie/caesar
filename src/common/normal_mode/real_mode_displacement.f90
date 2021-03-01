! ======================================================================
! A displacement in real mode co-ordinates.
! ======================================================================
module caesar_real_mode_displacement_module
  use caesar_utils_module
  
  use caesar_structure_module
  
  use caesar_mass_weighted_displacement_module
  use caesar_cartesian_displacement_module
  use caesar_real_mode_module
  use caesar_real_single_mode_displacement_module
  implicit none
  
  private
  
  public :: RealModeDisplacement
  public :: size
  public :: MassWeightedDisplacement
  public :: CartesianDisplacement
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
  type, extends(Stringsable) :: RealModeDisplacement
    type(RealSingleDisplacement), allocatable :: vectors(:)
  contains
    ! The component of the displacement along a given mode.
    generic,   public  :: displacement =>  &
                        & displacement_id, &
                        & displacement_mode
    procedure, private :: displacement_id
    procedure, private :: displacement_mode
    
    ! I/O.
    procedure, public :: read  => read_RealModeDisplacement
    procedure, public :: write => write_RealModeDisplacement
  end type
  
  interface RealModeDisplacement
    ! Constructors and size() module function.
    module function new_RealModeDisplacement(displacements) result(this) 
      type(RealSingleDisplacement), intent(in) :: displacements(:)
      type(RealModeDisplacement)               :: this
    end function
  
    module function new_RealModeDisplacement_RealModes(modes,displacements) &
       & result(this) 
      type(RealMode), intent(in) :: modes(:)
      real(dp),       intent(in) :: displacements(:)
      type(RealModeDisplacement) :: this
    end function
  end interface
  
  interface size
    module function size_RealModeDisplacement(this) result(output) 
      type(RealModeDisplacement), intent(in) :: this
      integer                                :: output
    end function
  end interface
  
  interface operator(*)
    ! ----------------------------------------------------------------------
    ! Arithmetic.
    ! ----------------------------------------------------------------------
    impure elemental module function multiply_real_RealModeDisplacement(this,that) result(output) 
      real(dp),                   intent(in) :: this
      type(RealModeDisplacement), intent(in) :: that
      type(RealModeDisplacement)             :: output
    end function
  
    impure elemental module function multiply_RealModeDisplacement_real(this,that) result(output) 
      type(RealModeDisplacement), intent(in) :: this
      real(dp),                   intent(in) :: that
      type(RealModeDisplacement)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_RealModeDisplacement_real(this, &
       & that) result(output) 
      type(RealModeDisplacement), intent(in) :: this
      real(dp),                   intent(in) :: that
      type(RealModeDisplacement)             :: output
    end function
  end interface
  
  interface operator(+)
    impure elemental module function add_RealModeDisplacement_RealModeDisplacement(this,that) result(output) 
      type(RealModeDisplacement), intent(in) :: this
      type(RealModeDisplacement), intent(in) :: that
      type(RealModeDisplacement)             :: output
    end function
  end interface
  
  interface sum
    module function sum_RealModeDisplacements(this) result(output) 
      type(RealModeDisplacement), intent(in) :: this(:)
      type(RealModeDisplacement)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function negative_RealModeDisplacement(this) &
       & result(output) 
      type(RealModeDisplacement), intent(in) :: this
      type(RealModeDisplacement)             :: output
    end function
  
    impure elemental module function subtract_RealModeDisplacement_RealModeDisplacement(   this,that) result(output) 
      type(RealModeDisplacement), intent(in) :: this
      type(RealModeDisplacement), intent(in) :: that
      type(RealModeDisplacement)             :: output
    end function
  end interface
  
  interface MassWeightedDisplacement
    ! ----------------------------------------------------------------------
    ! Conversions to and from cartesian and mass-weighted co-ordinates.
    ! ----------------------------------------------------------------------
    ! Returns the displacement in mass-weighted co-ordinates.
    module function new_MassWeightedDisplacement_RealModeDisplacement(this, &
       & structure,modes,qpoints) result(output) 
      class(RealModeDisplacement), intent(in) :: this
      type(StructureData),         intent(in) :: structure
      type(RealMode),              intent(in) :: modes(:)
      type(QpointData),            intent(in) :: qpoints(:)
      type(MassWeightedDisplacement)          :: output
    end function
  end interface
  
  interface CartesianDisplacement
    ! Returns the displacement in cartesian co-ordinates.
    module function new_CartesianDisplacement_RealModeDisplacement(this, &
       & structure,modes,qpoints) result(output) 
      class(RealModeDisplacement), intent(in) :: this
      type(StructureData),         intent(in) :: structure
      type(RealMode),              intent(in) :: modes(:)
      type(QpointData),            intent(in) :: qpoints(:)
      type(CartesianDisplacement)             :: output
    end function
  end interface
  
  interface RealModeDisplacement
    ! Converts a MassWeightedDisplacement to a RealModeDisplacement.
    module function new_RealModeDisplacement_MassWeightedDisplacement(displacement,structure,modes,qpoints) result(this) 
      type(MassWeightedDisplacement), intent(in) :: displacement
      type(StructureData),            intent(in) :: structure
      type(RealMode),                 intent(in) :: modes(:)
      type(QpointData),               intent(in) :: qpoints(:)
      type(RealModeDisplacement)                 :: this
    end function
  
    ! Converts a CartesianDisplacement to a RealModeDisplacement.
    module function new_RealModeDisplacement_CartesianDisplacement(displacement,structure,modes,qpoints) result(this) 
      type(CartesianDisplacement), intent(in) :: displacement
      type(StructureData),         intent(in) :: structure
      type(RealMode),              intent(in) :: modes(:)
      type(QpointData),            intent(in) :: qpoints(:)
      type(RealModeDisplacement)              :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! The displacement along a given mode.
    ! ----------------------------------------------------------------------
    impure elemental module function displacement_id(this,id) result(output) 
      class(RealModeDisplacement), intent(in) :: this
      integer,                     intent(in) :: id
      real(dp)                                :: output
    end function
  end interface
  
  interface
    impure elemental module function displacement_mode(this,mode) &
       & result(output) 
      class(RealModeDisplacement), intent(in) :: this
      type(RealMode),              intent(in) :: mode
      real(dp)                                :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_RealModeDisplacement(this,input) 
      class(RealModeDisplacement), intent(out) :: this
      type(String),                intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_RealModeDisplacement(this) result(output) 
      class(RealModeDisplacement), intent(in) :: this
      type(String), allocatable               :: output(:)
    end function
  end interface
  
  interface RealModeDisplacement
    module function new_RealModeDisplacement_Strings(input) result(this) 
      type(String), intent(in)   :: input(:)
      type(RealModeDisplacement) :: this
    end function
  
    impure elemental module function new_RealModeDisplacement_StringArray(input) result(this) 
      type(StringArray), intent(in) :: input
      type(RealModeDisplacement)    :: this
    end function
  end interface
end module
