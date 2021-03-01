! ======================================================================
! A displacement in complex mode co-ordinates.
! ======================================================================
module caesar_complex_mode_displacement_module
  use caesar_utils_module
  
  use caesar_structure_module
  
  use caesar_complex_mode_module
  use caesar_complex_single_mode_displacement_module
  implicit none
  
  private
  
  public :: ComplexModeDisplacement
  public :: size
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
  type, extends(Stringsable) :: ComplexModeDisplacement
    type(ComplexSingleDisplacement), allocatable :: vectors(:)
  contains
    ! The component of the displacement along a given mode.
    generic,   public  :: displacement =>  &
                        & displacement_id, &
                        & displacement_mode
    procedure, private :: displacement_id
    procedure, private :: displacement_mode
    
    ! I/O.
    procedure, public :: read  => read_ComplexModeDisplacement
    procedure, public :: write => write_ComplexModeDisplacement
  end type
  
  interface ComplexModeDisplacement
    ! Constructors and size() module function.
    module function new_ComplexModeDisplacement(displacements) result(this) 
      type(ComplexSingleDisplacement), intent(in) :: displacements(:)
      type(ComplexModeDisplacement)               :: this
    end function
  
    module function new_ComplexModeDisplacement_ComplexModes(modes, &
       & displacements) result(this) 
      type(ComplexMode), intent(in) :: modes(:)
      complex(dp),       intent(in) :: displacements(:)
      type(ComplexModeDisplacement) :: this
    end function
  end interface
  
  interface size
    module function size_ComplexModeDisplacement(this) result(output) 
      type(ComplexModeDisplacement), intent(in) :: this
      integer                                   :: output
    end function
  end interface
  
  interface operator(*)
    ! ----------------------------------------------------------------------
    ! Arithmetic.
    ! ----------------------------------------------------------------------
    impure elemental module function multiply_real_ComplexModeDisplacement(this,that) result(output) 
      real(dp),                      intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: that
      type(ComplexModeDisplacement)             :: output
    end function
  
    impure elemental module function multiply_ComplexModeDisplacement_real(this,that) result(output) 
      type(ComplexModeDisplacement), intent(in) :: this
      real(dp),                      intent(in) :: that
      type(ComplexModeDisplacement)             :: output
    end function
  
    impure elemental module function multiply_complex_ComplexModeDisplacement(this,that) result(output) 
      complex(dp),                   intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: that
      type(ComplexModeDisplacement)             :: output
    end function
  
    impure elemental module function multiply_ComplexModeDisplacement_complex(this,that) result(output) 
      type(ComplexModeDisplacement), intent(in) :: this
      complex(dp),                   intent(in) :: that
      type(ComplexModeDisplacement)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_ComplexModeDisplacement_complex(this,that) result(output) 
      type(ComplexModeDisplacement), intent(in) :: this
      complex(dp),                   intent(in) :: that
      type(ComplexModeDisplacement)             :: output
    end function
  end interface
  
  interface operator(+)
    impure elemental module function add_ComplexModeDisplacement_ComplexModeDisplacement(  this,that) result(output) 
      type(ComplexModeDisplacement), intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: that
      type(ComplexModeDisplacement)             :: output
    end function
  end interface
  
  interface sum
    module function sum_ComplexModeDisplacements(this) result(output) 
      type(ComplexModeDisplacement), intent(in) :: this(:)
      type(ComplexModeDisplacement)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function negative_ComplexModeDisplacement(this) &
       & result(output) 
      type(ComplexModeDisplacement), intent(in) :: this
      type(ComplexModeDisplacement)             :: output
    end function
  
    impure elemental module function                                                  subtract_ComplexModeDisplacement_ComplexModeDisplacement(this,that) result(output) 
      type(ComplexModeDisplacement), intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: that
      type(ComplexModeDisplacement)             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! The displacement along a given mode.
    ! ----------------------------------------------------------------------
    impure elemental module function displacement_id(this,id) result(output) 
      class(ComplexModeDisplacement), intent(in) :: this
      integer,                        intent(in) :: id
      complex(dp)                                :: output
    end function
  end interface
  
  interface
    impure elemental module function displacement_mode(this,mode) &
       & result(output) 
      class(ComplexModeDisplacement), intent(in) :: this
      type(ComplexMode),              intent(in) :: mode
      complex(dp)                                :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_ComplexModeDisplacement(this,input) 
      class(ComplexModeDisplacement), intent(out) :: this
      type(String),                   intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_ComplexModeDisplacement(this) result(output) 
      class(ComplexModeDisplacement), intent(in) :: this
      type(String), allocatable                  :: output(:)
    end function
  end interface
  
  interface ComplexModeDisplacement
    module function new_ComplexModeDisplacement_Strings(input) result(this) 
      type(String), intent(in)      :: input(:)
      type(ComplexModeDisplacement) :: this
    end function
  
    impure elemental module function new_ComplexModeDisplacement_StringArray(input) result(this) 
      type(StringArray), intent(in) :: input
      type(ComplexModeDisplacement) :: this
    end function
  end interface
end module
