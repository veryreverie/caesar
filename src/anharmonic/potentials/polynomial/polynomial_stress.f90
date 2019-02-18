! ======================================================================
! A polynomial representation of the stress.
! ======================================================================
module polynomial_stress_module
  use common_module
  
  use states_module
  use anharmonic_common_module
  
  use coupling_stress_basis_functions_module
  implicit none
  
  private
  
  public :: startup_polynomial_stress
  
  public :: PolynomialStress
  
  type, extends(StressData) :: PolynomialStress
    type(RealMatrix), private :: reference_stress_
    type(CouplingStressBasisFunctions), allocatable, private :: &
       & basis_functions_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_PolynomialStress
    ! I/O.
    procedure, public :: read  => read_PolynomialStress
    procedure, public :: write => write_PolynomialStress
  end type
  
  interface PolynomialStress
    module procedure new_PolynomialStress
    module procedure new_PolynomialStress_Strings
    module procedure new_PolynomialStress_StringArray
  end interface
contains

! Startup procedure.
subroutine startup_polynomial_stress()
  implicit none
  
  type(PolynomialStress) :: stress
  
  call stress%startup()
end subroutine

! Constructor.
function new_PolynomialStress(reference_stress,basis_functions) result(this)
  implicit none
  
  type(RealMatrix),                   intent(in) :: reference_stress
  type(CouplingStressBasisFunctions), intent(in) :: basis_functions(:)
  type(PolynomialStress)                         :: this
  
  this%reference_stress_ = reference_stress
  this%basis_functions_  = basis_functions
end function

! Type representation.
impure elemental function representation_PolynomialStress() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'polynomial'
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_PolynomialStress(this,input)
  implicit none
  
  class(PolynomialStress), intent(out) :: this
  type(String),               intent(in)  :: input(:)
  
  type(RealMatrix)                                :: reference_stress
  type(CouplingStressBasisFunctions), allocatable :: basis_functions(:)
  
  select type(this); type is(PolynomialStress)
    reference_stress = RealMatrix(input(2:4))
    
    basis_functions = CouplingStressBasisFunctions(split_into_sections( &
                                       & input(7:),                     &
                                       & separating_line=repeat('=',50) ))
    
    this = PolynomialStress( reference_stress, &
                           & basis_functions   )
  class default
    call err()
  end select
end subroutine

function write_PolynomialStress(this) result(output)
  implicit none
  
  class(PolynomialStress), intent(in) :: this
  type(String), allocatable              :: output(:)
  
  select type(this); type is(PolynomialStress)
    output = [ str('Reference stress:'),                                  &
             & str(this%reference_stress_),                               &
             & str('Basis functions:'),                                   &
             & str(''),                                                   &
             & str(this%basis_functions_, separating_line=repeat('=',50)) ]
  class default
    call err()
  end select
end function

function new_PolynomialStress_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(PolynomialStress)   :: this
  
  call this%read(input)
end function

impure elemental function new_PolynomialStress_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(PolynomialStress)        :: this
  
  this = PolynomialStress(str(input))
end function
end module
