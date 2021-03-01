! ======================================================================
! A displacement in cartesian co-ordinates.
! ======================================================================
module caesar_cartesian_displacement_module
  use caesar_utils_module
  
  use caesar_structure_module
  implicit none
  
  private
  
  public :: CartesianDisplacement
  public :: size
  public :: displace_structure
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
  type, extends(Stringsable) :: CartesianDisplacement
    type(RealVector), allocatable :: vectors(:)
  contains
    procedure, public :: read  => read_CartesianDisplacement
    procedure, public :: write => write_CartesianDisplacement
  end type
  
  interface CartesianDisplacement
    ! ----------------------------------------------------------------------
    ! Constructor and size() module function.
    ! ----------------------------------------------------------------------
    module function new_CartesianDisplacement(displacements) result(this) 
      type(RealVector), intent(in) :: displacements(:)
      type(CartesianDisplacement)  :: this
    end function
  end interface
  
  interface size
    module function size_CartesianDisplacement(this) result(output) 
      class(CartesianDisplacement), intent(in) :: this
      integer                                  :: output
    end function
  end interface
  
  interface CartesianDisplacement
    ! ----------------------------------------------------------------------
    ! Construct a zero displacement.
    ! ----------------------------------------------------------------------
    impure elemental module function new_CartesianDisplacement_zero(structure) result(this) 
      type(StructureData), intent(in) :: structure
      type(CartesianDisplacement)     :: this
    end function
  end interface
  
  interface displace_structure
    ! ----------------------------------------------------------------------
    ! Construct the structure which is displaced from the input structure.
    ! ----------------------------------------------------------------------
    module function displace_structure_CartesianDisplacement(structure, &
       & displacement) result(output) 
      type(StructureData),         intent(in) :: structure
      type(CartesianDisplacement), intent(in) :: displacement
      type(StructureData)                     :: output
    end function
  end interface
  
  interface operator(*)
    ! ----------------------------------------------------------------------
    ! Algebra.
    ! ----------------------------------------------------------------------
    impure elemental module function multiply_real_CartesianDisplacement(this,that) result(output) 
      real(dp),                    intent(in) :: this
      type(CartesianDisplacement), intent(in) :: that
      type(CartesianDisplacement)             :: output
    end function
  
    impure elemental module function multiply_CartesianDisplacement_real(this,that) result(output) 
      type(CartesianDisplacement), intent(in) :: this
      real(dp),                    intent(in) :: that
      type(CartesianDisplacement)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_CartesianDisplacement_real(this, &
       & that) result(output) 
      type(CartesianDisplacement), intent(in) :: this
      real(dp),                    intent(in) :: that
      type(CartesianDisplacement)             :: output
    end function
  end interface
  
  interface operator(+)
    impure elemental module function add_CartesianDisplacement_CartesianDisplacement(   this,that) result(output) 
      type(CartesianDisplacement), intent(in) :: this
      type(CartesianDisplacement), intent(in) :: that
      type(CartesianDisplacement)             :: output
    end function
  end interface
  
  interface sum
    module function sum_CartesianDisplacements(input) result(output) 
      type(CartesianDisplacement), intent(in) :: input(:)
      type(CartesianDisplacement)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function negative_CartesianDisplacement(this) &
       & result(output) 
      type(CartesianDisplacement), intent(in) :: this
      type(CartesianDisplacement)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function                                              subtract_CartesianDisplacement_CartesianDisplacement(this,that) result(output) 
      type(CartesianDisplacement), intent(in) :: this
      type(CartesianDisplacement), intent(in) :: that
      type(CartesianDisplacement)             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_CartesianDisplacement(this,input) 
      class(CartesianDisplacement), intent(out) :: this
      type(String),                 intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_CartesianDisplacement(this) result(output) 
      class(CartesianDisplacement), intent(in) :: this
      type(String), allocatable                :: output(:)
    end function
  end interface
  
  interface CartesianDisplacement
    module function new_CartesianDisplacement_Strings(input) result(this) 
      type(String), intent(in)    :: input(:)
      type(CartesianDisplacement) :: this
    end function
  end interface
  
  interface CartesianDisplacement
    impure elemental module function new_CartesianDisplacement_StringArray(input) result(this) 
      type(StringArray), intent(in) :: input
      type(CartesianDisplacement)   :: this
    end function
  end interface
end module
