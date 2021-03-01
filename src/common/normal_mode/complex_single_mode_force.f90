! ======================================================================
! A force along a single complex mode.
! ======================================================================
module caesar_complex_single_mode_force_module
  use caesar_utils_module
  
  use caesar_complex_mode_module
  implicit none
  
  private
  
  public :: ComplexSingleForce
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  public :: select_mode
  public :: select_modes
  
  type, extends(Stringable) :: ComplexSingleForce
    ! The id of the mode.
    integer :: id
    
    ! The magnitude of the force along the mode.
    complex(dp) :: magnitude
  contains
    procedure, public :: read  => read_ComplexSingleForce
    procedure, public :: write => write_ComplexSingleForce
  end type
  
  interface ComplexSingleForce
    ! Constructors.
    impure elemental module function new_ComplexSingleForce(id,magnitude) &
       & result(this) 
      integer,     intent(in)  :: id
      complex(dp), intent(in)  :: magnitude
      type(ComplexSingleForce) :: this
    end function
  
    impure elemental module function new_ComplexSingleForce_ComplexMode(mode,magnitude) result(this) 
      type(ComplexMode), intent(in) :: mode
      complex(dp),       intent(in) :: magnitude
      type(ComplexSingleForce)      :: this
    end function
  end interface
  
  interface operator(*)
    ! Arithmetic.
    impure elemental module function multiply_real_ComplexSingleForce(this, &
       & that) result(output) 
      real(dp),                 intent(in) :: this
      type(ComplexSingleForce), intent(in) :: that
      type(ComplexSingleForce)             :: output
    end function
  
    impure elemental module function multiply_ComplexSingleForce_real(this, &
       & that) result(output) 
      type(ComplexSingleForce), intent(in) :: this
      real(dp),                 intent(in) :: that
      type(ComplexSingleForce)             :: output
    end function
  
    impure elemental module function multiply_complex_ComplexSingleForce(this,that) result(output) 
      complex(dp),              intent(in) :: this
      type(ComplexSingleForce), intent(in) :: that
      type(ComplexSingleForce)             :: output
    end function
  
    impure elemental module function multiply_ComplexSingleForce_complex(this,that) result(output) 
      type(ComplexSingleForce), intent(in) :: this
      complex(dp),              intent(in) :: that
      type(ComplexSingleForce)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_ComplexSingleForce_complex(this, &
       & that) result(output) 
      type(ComplexSingleForce), intent(in) :: this
      complex(dp),              intent(in) :: that
      type(ComplexSingleForce)             :: output
    end function
  end interface
  
  interface operator(+)
    impure elemental module function add_ComplexSingleForce_ComplexSingleForce(   this,that) result(output) 
      type(ComplexSingleForce), intent(in) :: this
      type(ComplexSingleForce), intent(in) :: that
      type(ComplexSingleForce)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function negative_ComplexSingleForce(this) &
       & result(output) 
      type(ComplexSingleForce), intent(in) :: this
      type(ComplexSingleForce)             :: output
    end function
  
    impure elemental module function subtract_ComplexSingleForce_ComplexSingleForce(   this,that) result(output) 
      type(ComplexSingleForce), intent(in) :: this
      type(ComplexSingleForce), intent(in) :: that
      type(ComplexSingleForce)             :: output
    end function
  end interface
  
  interface select_mode
    ! Select modes corresponding to a given force or forces.
    module function select_mode_ComplexSingleForce(force,modes) result(output) 
      type(ComplexSingleForce), intent(in) :: force
      type(ComplexMode),        intent(in) :: modes(:)
      type(ComplexMode)                    :: output
    end function
  end interface
  
  interface select_modes
    module function select_modes_ComplexSingleForces(forces,modes) &
       & result(output) 
      type(ComplexSingleForce), intent(in) :: forces(:)
      type(ComplexMode),        intent(in) :: modes(:)
      type(ComplexMode), allocatable       :: output(:)
    end function
  end interface
  
  interface
    ! I/O.
    module subroutine read_ComplexSingleForce(this,input) 
      class(ComplexSingleForce), intent(out) :: this
      type(String),              intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_ComplexSingleForce(this) result(output) 
      class(ComplexSingleForce), intent(in) :: this
      type(String)                          :: output
    end function
  end interface
  
  interface ComplexSingleForce
    impure elemental module function new_ComplexSingleForce_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(ComplexSingleForce) :: this
    end function
  end interface
end module
