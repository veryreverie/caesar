!> Provides the [[PhaseData(type)]] class and related methods.
module caesar_phase_module
  use caesar_foundations_module
  use caesar_io_module
  
  use caesar_fraction_module
  use caesar_mathematical_constants_module
  implicit none
  
  private
  
  public :: PhaseData
  public :: cmplx
  public :: operator(==)
  public :: operator(/=)
  public :: calculate_phase
  
  !> A complex phase, of the form `exp(2*pi*i*theta)`.
  type, extends(Stringable) :: PhaseData
    type(IntFraction), private :: theta_
  contains
    procedure, public :: theta => theta_PhaseData
    
    procedure, public :: read  => read_PhaseData
    procedure, public :: write => write_PhaseData
  end type
  
  interface PhaseData
    ! Constructor for [[PhaseData(type)]].
    module function new_PhaseData(theta) result(this) 
      type(IntFraction), intent(in) :: theta
      type(PhaseData)               :: this
    end function
  end interface
  
  interface
    module function theta_PhaseData(this) result(output)
      class(PhaseData), intent(in) :: this
      type(IntFraction)            :: output
    end function
  end interface
  
  interface cmplx
    ! Conversion to complex(dp).
    module function cmplx_PhaseData(this) result(output) 
      type(PhaseData) :: this
      complex(dp)     :: output
    end function
  end interface
  
  interface
    !> Calculates the phase of a complex number.
    !> Stores the phase as an exact fraction, with the given denominator
    !>    (or a simplified form of this fraction if possible).
    module function calculate_phase(input,denominator) result(output) 
      complex(dp), intent(in) :: input
      integer,     intent(in) :: denominator
      type(PhaseData)         :: output
    end function
  end interface
  
  interface operator(==)
    ! ----------------------------------------------------------------------
    ! Comparison.
    ! ----------------------------------------------------------------------
    impure elemental module function equality_PhaseData_PhaseData(this,that) &
       & result(output) 
      type(PhaseData), intent(in) :: this
      type(PhaseData), intent(in) :: that
      logical                     :: output
    end function
  end interface
  
  interface operator(/=)
    impure elemental module function non_equality_PhaseData_PhaseData(this, &
       & that) result(output) 
      type(PhaseData), intent(in) :: this
      type(PhaseData), intent(in) :: that
      logical                     :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_PhaseData(this,input) 
      class(PhaseData), intent(out) :: this
      type(String),     intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_PhaseData(this) result(output) 
      class(PhaseData), intent(in) :: this
      type(String)                 :: output
    end function
  end interface
  
  interface PhaseData
    impure elemental module function new_PhaseData_String(input) result(this) 
      type(String), intent(in) :: input
      type(PhaseData)          :: this
    end function
  end interface
end module
