! ======================================================================
! A complex phase, of the form exp(2*pi*i*frac), where frac is an integer
!    fraction.
! ======================================================================
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
  
  type, extends(Stringable) :: PhaseData
    type(IntFraction) :: fraction
  contains
    procedure, public :: read  => read_PhaseData
    procedure, public :: write => write_PhaseData
  end type
  
  interface PhaseData
    ! Constructor.
    module function new_PhaseData(input) result(this) 
      type(IntFraction), intent(in) :: input
      type(PhaseData)               :: this
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
    ! Finds the exact phase of a complex number.
    module function calculate_phase(input,denom) result(output) 
      complex(dp), intent(in) :: input
      integer,     intent(in) :: denom
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
