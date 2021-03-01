! ======================================================================
! A displacement along a single complex mode.
! ======================================================================
module caesar_complex_single_mode_displacement_module
  use caesar_utils_module
  
  use caesar_complex_mode_module
  implicit none
  
  private
  
  public :: ComplexSingleDisplacement
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  public :: select_mode
  public :: select_modes
  
  type, extends(Stringable) :: ComplexSingleDisplacement
    ! The id of the mode.
    integer :: id
    
    ! The magnitude of the displacement along the mode.
    complex(dp) :: magnitude
  contains
    procedure, public :: read  => read_ComplexSingleDisplacement
    procedure, public :: write => write_ComplexSingleDisplacement
  end type
  
  interface ComplexSingleDisplacement
    ! Constructors.
    impure elemental module function new_ComplexSingleDisplacement(id, &
       & magnitude) result(this) 
      integer,     intent(in)         :: id
      complex(dp), intent(in)         :: magnitude
      type(ComplexSingleDisplacement) :: this
    end function
  
    impure elemental module function new_ComplexSingleDisplacement_ComplexMode(mode,magnitude) result(this) 
      type(ComplexMode), intent(in)   :: mode
      complex(dp),       intent(in)   :: magnitude
      type(ComplexSingleDisplacement) :: this
    end function
  end interface
  
  interface operator(*)
    ! Arithmetic.
    impure elemental module function multiply_real_ComplexSingleDisplacement(this,that) result(output) 
      real(dp),                        intent(in) :: this
      type(ComplexSingleDisplacement), intent(in) :: that
      type(ComplexSingleDisplacement)             :: output
    end function
  
    impure elemental module function multiply_ComplexSingleDisplacement_real(this,that) result(output) 
      type(ComplexSingleDisplacement), intent(in) :: this
      real(dp),                        intent(in) :: that
      type(ComplexSingleDisplacement)             :: output
    end function
  
    impure elemental module function multiply_complex_ComplexSingleDisplacement(this,that) result(output) 
      complex(dp),                     intent(in) :: this
      type(ComplexSingleDisplacement), intent(in) :: that
      type(ComplexSingleDisplacement)             :: output
    end function
  
    impure elemental module function multiply_ComplexSingleDisplacement_complex(this,that) result(output) 
      type(ComplexSingleDisplacement), intent(in) :: this
      complex(dp),                     intent(in) :: that
      type(ComplexSingleDisplacement)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_ComplexSingleDisplacement_complex(this,that) result(output) 
      type(ComplexSingleDisplacement), intent(in) :: this
      complex(dp),                     intent(in) :: that
      type(ComplexSingleDisplacement)             :: output
    end function
  end interface
  
  interface operator(+)
    impure elemental module function                                                 add_ComplexSingleDisplacement_ComplexSingleDisplacement(this,that) result(output) 
      type(ComplexSingleDisplacement), intent(in) :: this
      type(ComplexSingleDisplacement), intent(in) :: that
      type(ComplexSingleDisplacement)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function negative_ComplexSingleDisplacement(this) &
       & result(output) 
      type(ComplexSingleDisplacement), intent(in) :: this
      type(ComplexSingleDisplacement)             :: output
    end function
  
    impure elemental module function                                                      subtract_ComplexSingleDisplacement_ComplexSingleDisplacement(this,that) result(output) 
      type(ComplexSingleDisplacement), intent(in) :: this
      type(ComplexSingleDisplacement), intent(in) :: that
      type(ComplexSingleDisplacement)             :: output
    end function
  end interface
  
  interface select_mode
    ! Select modes corresponding to a given force or forces.
    module function select_mode_ComplexSingleDisplacement(displacement,modes) &
       & result(output) 
      type(ComplexSingleDisplacement), intent(in) :: displacement
      type(ComplexMode),               intent(in) :: modes(:)
      type(ComplexMode)                           :: output
    end function
  end interface
  
  interface select_modes
    module function select_modes_ComplexSingleDisplacements(displacements, &
       & modes) result(output) 
      type(ComplexSingleDisplacement), intent(in) :: displacements(:)
      type(ComplexMode),               intent(in) :: modes(:)
      type(ComplexMode), allocatable              :: output(:)
    end function
  end interface
  
  interface
    ! I/O.
    module subroutine read_ComplexSingleDisplacement(this,input) 
      class(ComplexSingleDisplacement), intent(out) :: this
      type(String),                     intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_ComplexSingleDisplacement(this) result(output) 
      class(ComplexSingleDisplacement), intent(in) :: this
      type(String)                                 :: output
    end function
  end interface
  
  interface ComplexSingleDisplacement
    impure elemental module function new_ComplexSingleDisplacement_String(input) result(this) 
      type(String), intent(in)        :: input
      type(ComplexSingleDisplacement) :: this
    end function
  end interface
end module
