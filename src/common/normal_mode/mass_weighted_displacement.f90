! ======================================================================
! A displacement in mass-weighted cartesian co-ordinates.
! ======================================================================
module caesar_mass_weighted_displacement_module
  use caesar_utils_module
  
  use caesar_structure_module
  
  use caesar_cartesian_displacement_module
  implicit none
  
  private
  
  public :: MassWeightedDisplacement
  public :: size
  public :: CartesianDisplacement
  public :: displace_structure
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
  type, extends(Stringsable) :: MassWeightedDisplacement
    type(RealVector), allocatable :: vectors(:)
  contains
    procedure, public :: read  => read_MassWeightedDisplacement
    procedure, public :: write => write_MassWeightedDisplacement
  end type
  
  interface MassWeightedDisplacement
    ! ----------------------------------------------------------------------
    ! Constructor and size() module function.
    ! ----------------------------------------------------------------------
    module function new_MassWeightedDisplacement(displacements) result(this) 
      type(RealVector), intent(in)   :: displacements(:)
      type(MassWeightedDisplacement) :: this
    end function
  end interface
  
  interface size
    module function size_MassWeightedDisplacement(this) result(output) 
      type(MassWeightedDisplacement), intent(in) :: this
      integer                                    :: output
    end function
  end interface
  
  interface MassWeightedDisplacement
    ! ----------------------------------------------------------------------
    ! Construct a zero displacement.
    ! ----------------------------------------------------------------------
    impure elemental module function new_MassWeightedDisplacement_zero(structure) result(this) 
      type(StructureData), intent(in) :: structure
      type(MassWeightedDisplacement)  :: this
    end function
  
    ! ----------------------------------------------------------------------
    ! Conversion to and from non-mass-weighted co-ordinates.
    ! ----------------------------------------------------------------------
    module function new_MassWeightedDisplacement_CartesianDisplacement(input,structure) result(output) 
      type(CartesianDisplacement), intent(in) :: input
      type(StructureData),         intent(in) :: structure
      type(MassWeightedDisplacement)          :: output
    end function
  end interface
  
  interface CartesianDisplacement
    module function new_CartesianDisplacement_MassWeightedDisplacement(input,structure) result(output) 
      type(MassWeightedDisplacement), intent(in) :: input
      type(StructureData),            intent(in) :: structure
      type(CartesianDisplacement)                :: output
    end function
  end interface
  
  interface displace_structure
    ! ----------------------------------------------------------------------
    ! Construct the structure which is displaced from the input structure.
    ! ----------------------------------------------------------------------
    module function displace_structure_MassWeightedDisplacement(structure, &
       & displacement) result(output) 
      type(StructureData),            intent(in) :: structure
      type(MassWeightedDisplacement), intent(in) :: displacement
      type(StructureData)                        :: output
    end function
  end interface
  
  interface operator(*)
    ! ----------------------------------------------------------------------
    ! Algebra.
    ! ----------------------------------------------------------------------
    impure elemental module function multiply_real_MassWeightedDisplacement(this,that) result(output) 
      real(dp),                       intent(in) :: this
      type(MassWeightedDisplacement), intent(in) :: that
      type(MassWeightedDisplacement)             :: output
    end function
  
    impure elemental module function multiply_MassWeightedDisplacement_real(this,that) result(output) 
      type(MassWeightedDisplacement), intent(in) :: this
      real(dp),                       intent(in) :: that
      type(MassWeightedDisplacement)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_MassWeightedDisplacement_real(this,that) result(output) 
      type(MassWeightedDisplacement), intent(in) :: this
      real(dp),                       intent(in) :: that
      type(MassWeightedDisplacement)             :: output
    end function
  end interface
  
  interface operator(+)
    impure elemental module function                                               add_MassWeightedDisplacement_MassWeightedDisplacement(this,that) result(output) 
      type(MassWeightedDisplacement), intent(in) :: this
      type(MassWeightedDisplacement), intent(in) :: that
      type(MassWeightedDisplacement)             :: output
    end function
  end interface
  
  interface sum
    module function sum_MassWeightedDisplacements(input) result(output) 
      type(MassWeightedDisplacement), intent(in) :: input(:)
      type(MassWeightedDisplacement)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function negative_MassWeightedDisplacement(this) &
       & result(output) 
      type(MassWeightedDisplacement), intent(in) :: this
      type(MassWeightedDisplacement)             :: output
    end function
  
    impure elemental module function                                                    subtract_MassWeightedDisplacement_MassWeightedDisplacement(this,that) result(output) 
      type(MassWeightedDisplacement), intent(in) :: this
      type(MassWeightedDisplacement), intent(in) :: that
      type(MassWeightedDisplacement)             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_MassWeightedDisplacement(this,input) 
      class(MassWeightedDisplacement), intent(out) :: this
      type(String),                    intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_MassWeightedDisplacement(this) result(output) 
      class(MassWeightedDisplacement), intent(in) :: this
      type(String), allocatable                   :: output(:)
    end function
  end interface
  
  interface MassWeightedDisplacement
    module function new_MassWeightedDisplacement_Strings(input) result(this) 
      type(String), intent(in)       :: input(:)
      type(MassWeightedDisplacement) :: this
    end function
  end interface
  
  interface MassWeightedDisplacement
    impure elemental module function new_MassWeightedDisplacement_StringArray(input) result(this) 
      type(StringArray), intent(in)  :: input
      type(MassWeightedDisplacement) :: this
    end function
  end interface
end module
