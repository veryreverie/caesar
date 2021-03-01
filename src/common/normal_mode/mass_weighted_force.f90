! ======================================================================
! A force in mass-weighted cartesian co-ordinates.
! ======================================================================
module caesar_mass_weighted_force_module
  use caesar_utils_module
  
  use caesar_structure_module
  
  use caesar_cartesian_force_module
  implicit none
  
  private
  
  public :: MassWeightedForce
  public :: size
  public :: CartesianForce
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  public :: sum
  
  type, extends(Stringsable) :: MassWeightedForce
    type(RealVector), allocatable :: vectors(:)
  contains
    procedure, public :: read  => read_MassWeightedForce
    procedure, public :: write => write_MassWeightedForce
  end type
  
  interface MassWeightedForce
    ! ----------------------------------------------------------------------
    ! Constructor and size() module function.
    ! ----------------------------------------------------------------------
    module function new_MassWeightedForce(forces) result(this) 
      type(RealVector), intent(in) :: forces(:)
      type(MassWeightedForce)      :: this
    end function
  end interface
  
  interface size
    module function size_MassWeightedForce(this) result(output) 
      type(MassWeightedForce), intent(in) :: this
      integer                             :: output
    end function
  end interface
  
  interface MassWeightedForce
    ! ----------------------------------------------------------------------
    ! Construct a zero force.
    ! ----------------------------------------------------------------------
    impure elemental module function new_MassWeightedForce_zero(structure) &
       & result(this) 
      type(StructureData), intent(in) :: structure
      type(MassWeightedForce)         :: this
    end function
  
    ! ----------------------------------------------------------------------
    ! Conversion to and from non-mass-weighted co-ordinates.
    ! ----------------------------------------------------------------------
    module function new_MassWeightedForce_CartesianForce(input,structure) &
       & result(output) 
      type(CartesianForce), intent(in) :: input
      type(StructureData),  intent(in) :: structure
      type(MassWeightedForce)          :: output
    end function
  end interface
  
  interface CartesianForce
    module function new_CartesianForce_MassWeightedForce(input,structure) &
       & result(output) 
      type(MassWeightedForce), intent(in) :: input
      type(StructureData),     intent(in) :: structure
      type(CartesianForce)                :: output
    end function
  end interface
  
  interface operator(*)
    ! ----------------------------------------------------------------------
    ! Algebra.
    ! ----------------------------------------------------------------------
    impure elemental module function multiply_real_MassWeightedForce(this, &
       & that) result(output) 
      real(dp),                intent(in) :: this
      type(MassWeightedForce), intent(in) :: that
      type(MassWeightedForce)             :: output
    end function
  
    impure elemental module function multiply_MassWeightedForce_real(this, &
       & that) result(output) 
      type(MassWeightedForce), intent(in) :: this
      real(dp),                intent(in) :: that
      type(MassWeightedForce)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_MassWeightedForce_real(this,that) &
       & result(output) 
      type(MassWeightedForce), intent(in) :: this
      real(dp),                intent(in) :: that
      type(MassWeightedForce)             :: output
    end function
  end interface
  
  interface operator(+)
    impure elemental module function add_MassWeightedForce_MassWeightedForce(this,that) result(output) 
      type(MassWeightedForce), intent(in) :: this
      type(MassWeightedForce), intent(in) :: that
      type(MassWeightedForce)             :: output
    end function
  end interface
  
  interface sum
    module function sum_MassWeightedForces(input) result(output) 
      type(MassWeightedForce), intent(in) :: input(:)
      type(MassWeightedForce)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function negative_MassWeightedForce(this) &
       & result(output) 
      type(MassWeightedForce), intent(in) :: this
      type(MassWeightedForce)             :: output
    end function
  
    impure elemental module function subtract_MassWeightedForce_MassWeightedForce(this,that) result(output) 
      type(MassWeightedForce), intent(in) :: this
      type(MassWeightedForce), intent(in) :: that
      type(MassWeightedForce)             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_MassWeightedForce(this,input) 
      class(MassWeightedForce), intent(out) :: this
      type(String),             intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_MassWeightedForce(this) result(output) 
      class(MassWeightedForce), intent(in) :: this
      type(String), allocatable            :: output(:)
    end function
  end interface
  
  interface MassWeightedForce
    module function new_MassWeightedForce_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(MassWeightedForce)  :: this
    end function
  
    impure elemental module function new_MassWeightedForce_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(MassWeightedForce)       :: this
    end function
  end interface
end module
