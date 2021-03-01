! ======================================================================
! A displacement along a single real mode, in mass-weighted co-ordinates.
! ======================================================================
module caesar_real_single_mode_displacement_module
  use caesar_utils_module
  
  use caesar_structure_module
  
  use caesar_cartesian_displacement_module
  use caesar_mass_weighted_displacement_module
  use caesar_real_mode_module
  implicit none
  
  private
  
  public :: RealSingleDisplacement
  public :: MassWeightedDisplacement
  public :: CartesianDisplacement
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  public :: select_mode
  public :: select_modes
  
  type, extends(Stringable) :: RealSingleDisplacement
    ! The id of the mode.
    integer :: id
    
    ! The magnitude of vector along the mode, in mass-weighted co-ordinates.
    real(dp) :: magnitude
  contains
    ! I/O.
    procedure, public :: read  => read_RealSingleDisplacement
    procedure, public :: write => write_RealSingleDisplacement
  end type
  
  interface RealSingleDisplacement
    ! Constructors.
    impure elemental module function new_RealSingleDisplacement(id,magnitude) &
       & result(this) 
      integer,  intent(in)         :: id
      real(dp), intent(in)         :: magnitude
      type(RealSingleDisplacement) :: this
    end function
  
    impure elemental module function new_RealSingleDisplacement_RealMode(mode,magnitude) result(this) 
      type(RealMode), intent(in)   :: mode
      real(dp),       intent(in)   :: magnitude
      type(RealSingleDisplacement) :: this
    end function
  end interface
  
  interface operator(*)
    ! ----------------------------------------------------------------------
    ! Arithmetic.
    ! ----------------------------------------------------------------------
    impure elemental module function multiply_real_RealSingleDisplacement(this,that) result(output) 
      real(dp),                     intent(in) :: this
      type(RealSingleDisplacement), intent(in) :: that
      type(RealSingleDisplacement)             :: output
    end function
  
    impure elemental module function multiply_RealSingleDisplacement_real(this,that) result(output) 
      type(RealSingleDisplacement), intent(in) :: this
      real(dp),                     intent(in) :: that
      type(RealSingleDisplacement)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_RealSingleDisplacement_real(this,that) result(output) 
      type(RealSingleDisplacement), intent(in) :: this
      real(dp),                     intent(in) :: that
      type(RealSingleDisplacement)             :: output
    end function
  end interface
  
  interface operator(+)
    impure elemental module function add_RealSingleDisplacement_RealSingleDisplacement(   this,that) result(output) 
      type(RealSingleDisplacement), intent(in) :: this
      type(RealSingleDisplacement), intent(in) :: that
      type(RealSingleDisplacement)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function negative_RealSingleDisplacement(this) &
       & result(output) 
      type(RealSingleDisplacement), intent(in) :: this
      type(RealSingleDisplacement)             :: output
    end function
  
    impure elemental module function                                                subtract_RealSingleDisplacement_RealSingleDisplacement(this,that) result(output) 
      type(RealSingleDisplacement), intent(in) :: this
      type(RealSingleDisplacement), intent(in) :: that
      type(RealSingleDisplacement)             :: output
    end function
  end interface
  
  interface MassWeightedDisplacement
    ! ----------------------------------------------------------------------
    ! Conversions to and from cartesian and mass-weighted co-ordinates.
    ! ----------------------------------------------------------------------
    ! Constructs the MassWeightedDisplacement corresponding to this vector.
    impure elemental module function                                                       new_MassWeightedDisplacement_RealSingleDisplacement(this,mode,structure,qpoint) result(output) 
      class(RealSingleDisplacement), intent(in) :: this
      type(RealMode),                intent(in) :: mode
      type(StructureData),           intent(in) :: structure
      type(QpointData),              intent(in) :: qpoint
      type(MassWeightedDisplacement)            :: output
    end function
  end interface
  
  interface CartesianDisplacement
    ! Constructs the CartesianDisplacement corresponding to this vector.
    impure elemental module function new_CartesianDisplacement_RealSingleDisplacement(   this,mode,structure,qpoint) result(output) 
      class(RealSingleDisplacement), intent(in) :: this
      type(RealMode),                intent(in) :: mode
      type(StructureData),           intent(in) :: structure
      type(QpointData),              intent(in) :: qpoint
      type(CartesianDisplacement)               :: output
    end function
  end interface
  
  interface RealSingleDisplacement
    ! Constructs the vector corresponding to the component of a mass-weighted
    !    vector along this mode.
    impure elemental module function                                               new_RealSingleDisplacement_MassWeightedDisplacement(mode,vector,structure,qpoint) result(this) 
      type(RealMode),                 intent(in) :: mode
      type(MassWeightedDisplacement), intent(in) :: vector
      type(StructureData),            intent(in) :: structure
      type(QpointData),               intent(in) :: qpoint
      type(RealSingleDisplacement)               :: this
    end function
  
    ! Constructs the vector corresponding to the component of a cartesian
    !    vector along this mode.
    impure elemental module function new_RealSingleDisplacement_CartesianDisplacement(   mode,vector,structure,qpoint) result(this) 
      type(RealMode),              intent(in) :: mode
      type(CartesianDisplacement), intent(in) :: vector
      type(StructureData),         intent(in) :: structure
      type(QpointData),            intent(in) :: qpoint
      type(RealSingleDisplacement)            :: this
    end function
  end interface
  
  interface select_mode
    ! ----------------------------------------------------------------------
    ! Select modes corresponding to a given displacement or displacements.
    ! ----------------------------------------------------------------------
    module function select_mode_RealSingleDisplacement(displacement,modes) &
       & result(output) 
      type(RealSingleDisplacement), intent(in) :: displacement
      type(RealMode),        intent(in) :: modes(:)
      type(RealMode)                    :: output
    end function
  end interface
  
  interface select_modes
    module function select_modes_RealSingleDisplacements(displacements,modes) &
       & result(output) 
      type(RealSingleDisplacement), intent(in) :: displacements(:)
      type(RealMode),               intent(in) :: modes(:)
      type(RealMode), allocatable              :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_RealSingleDisplacement(this,input) 
      class(RealSingleDisplacement), intent(out) :: this
      type(String),                  intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_RealSingleDisplacement(this) result(output) 
      class(RealSingleDisplacement), intent(in) :: this
      type(String)                              :: output
    end function
  end interface
  
  interface RealSingleDisplacement
    impure elemental module function new_RealSingleDisplacement_String(input) &
       & result(this) 
      type(String), intent(in)     :: input
      type(RealSingleDisplacement) :: this
    end function
  end interface
end module
