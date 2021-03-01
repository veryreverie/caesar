! ======================================================================
! A force in complex mode co-ordinates.
! ======================================================================
module caesar_complex_mode_force_module
  use caesar_utils_module
  
  use caesar_structure_module
  
  use caesar_complex_mode_module
  use caesar_complex_single_mode_force_module
  implicit none
  
  private
  
  public :: ComplexModeForce
  public :: size
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
  type, extends(Stringsable) :: ComplexModeForce
    type(ComplexSingleForce), allocatable :: vectors(:)
  contains
    ! The component of the force along a given mode.
    generic,   public  :: force =>  &
                        & force_id, &
                        & force_mode
    procedure, private :: force_id
    procedure, private :: force_mode
    
    ! I/O.
    procedure, public :: read  => read_ComplexModeForce
    procedure, public :: write => write_ComplexModeForce
  end type
  
  interface ComplexModeForce
    ! Constructors and size() module function.
    module function new_ComplexModeForce(forces) result(this) 
      type(ComplexSingleForce), intent(in) :: forces(:)
      type(ComplexModeForce)               :: this
    end function
  
    module function new_ComplexModeForce_ComplexModes(modes,forces) &
       & result(this) 
      type(ComplexMode), intent(in) :: modes(:)
      complex(dp),       intent(in) :: forces(:)
      type(ComplexModeForce)        :: this
    end function
  end interface
  
  interface size
    module function size_ComplexModeForce(this) result(output) 
      type(ComplexModeForce), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface operator(*)
    ! ----------------------------------------------------------------------
    ! Arithmetic.
    ! ----------------------------------------------------------------------
    impure elemental module function multiply_real_ComplexModeForce(this, &
       & that) result(output) 
      real(dp),               intent(in) :: this
      type(ComplexModeForce), intent(in) :: that
      type(ComplexModeForce)             :: output
    end function
  
    impure elemental module function multiply_ComplexModeForce_real(this, &
       & that) result(output) 
      type(ComplexModeForce), intent(in) :: this
      real(dp),               intent(in) :: that
      type(ComplexModeForce)             :: output
    end function
  
    impure elemental module function multiply_complex_ComplexModeForce(this, &
       & that) result(output) 
      complex(dp),            intent(in) :: this
      type(ComplexModeForce), intent(in) :: that
      type(ComplexModeForce)             :: output
    end function
  
    impure elemental module function multiply_ComplexModeForce_complex(this, &
       & that) result(output) 
      type(ComplexModeForce), intent(in) :: this
      complex(dp),            intent(in) :: that
      type(ComplexModeForce)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_ComplexModeForce_complex(this, &
       & that) result(output) 
      type(ComplexModeForce), intent(in) :: this
      complex(dp),            intent(in) :: that
      type(ComplexModeForce)             :: output
    end function
  end interface
  
  interface operator(+)
    impure elemental module function add_ComplexModeForce_ComplexModeForce(this,that) result(output) 
      type(ComplexModeForce), intent(in) :: this
      type(ComplexModeForce), intent(in) :: that
      type(ComplexModeForce)             :: output
    end function
  end interface
  
  interface sum
    module function sum_ComplexModeForces(this) result(output) 
      type(ComplexModeForce), intent(in) :: this(:)
      type(ComplexModeForce)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function negative_ComplexModeForce(this) &
       & result(output) 
      type(ComplexModeForce), intent(in) :: this
      type(ComplexModeForce)             :: output
    end function
  
    impure elemental module function subtract_ComplexModeForce_ComplexModeForce(   this,that) result(output) 
      type(ComplexModeForce), intent(in) :: this
      type(ComplexModeForce), intent(in) :: that
      type(ComplexModeForce)             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! The displacement along a given mode.
    ! ----------------------------------------------------------------------
    impure elemental module function force_id(this,id) result(output) 
      class(ComplexModeForce), intent(in) :: this
      integer,                 intent(in) :: id
      complex(dp)                         :: output
    end function
  end interface
  
  interface
    impure elemental module function force_mode(this,mode) result(output) 
      class(ComplexModeForce), intent(in) :: this
      type(ComplexMode),       intent(in) :: mode
      complex(dp)                         :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_ComplexModeForce(this,input) 
      class(ComplexModeForce), intent(out) :: this
      type(String),            intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_ComplexModeForce(this) result(output) 
      class(ComplexModeForce), intent(in) :: this
      type(String), allocatable           :: output(:)
    end function
  end interface
  
  interface ComplexModeForce
    module function new_ComplexModeForce_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(ComplexModeForce)   :: this
    end function
  
    impure elemental module function new_ComplexModeForce_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(ComplexModeForce)        :: this
    end function
  end interface
end module
