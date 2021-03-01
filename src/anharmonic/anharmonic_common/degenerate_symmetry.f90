! ======================================================================
! Symmetries in normal mode co-ordinates.
! ======================================================================
module caesar_degenerate_symmetry_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: DegenerateSymmetry
  
  ! The mode u_i, and the action of the symmetry on u_i, S.u_i, given in terms
  !    of other modes {u_j}.
  ! S.u_i = sum_j a_j u_j
  ! symmetric_mode_ids gives the set of j for which a_j is non-zero,
  ! and symmetric_mode_coefficients are the corresponding a_j.
  type, extends(Stringable) :: SingleModeSymmetry
    integer                  :: mode_id
    integer,     allocatable :: symmetric_mode_ids(:)
    complex(dp), allocatable :: symmetric_mode_coefficients(:)
  contains
    procedure, public :: read  => read_SingleModeSymmetry
    procedure, public :: write => write_SingleModeSymmetry
  end type
  
  ! A mode-by-mode list of SingleModeSymmetry.
  type, extends(Stringsable) :: DegenerateSymmetry
    integer                                        :: symmetry_id
    type(SingleModeSymmetry), allocatable, private :: symmetries_(:)
  contains
    procedure, public :: calculate_symmetry
    
    procedure, public :: transform_monomial
    
    ! I/O.
    procedure, public :: read  => read_DegenerateSymmetry
    procedure, public :: write => write_DegenerateSymmetry
  end type
  
  interface SingleModeSymmetry
    module function new_SingleModeSymmetry(mode_id,symmetric_mode_ids, &
       & symmetric_mode_coefficients) result(this) 
      integer,     intent(in)  :: mode_id
      integer,     intent(in)  :: symmetric_mode_ids(:)
      complex(dp), intent(in)  :: symmetric_mode_coefficients(:)
      type(SingleModeSymmetry) :: this
    end function
  end interface
  
  interface DegenerateSymmetry
    module function new_DegenerateSymmetry(symmetry,subspaces,modes,qpoints) &
       & result(this) 
      type(SymmetryOperator),   intent(in)    :: symmetry
      type(DegenerateSubspace), intent(in)    :: subspaces(:)
      type(ComplexMode),        intent(in)    :: modes(:)
      type(QpointData),         intent(in)    :: qpoints(:)
      type(DegenerateSymmetry)                :: this
    end function
  end interface
  
  interface
    ! Calculate the symmetry in a basis of complex monomials.
    ! If include_symmetry is .false. then the symmetry will not include the
    !    coefficients of the basis monomials.
    ! If include_symmetry is .true., then the coefficients will be included.
    module function calculate_symmetry(this,input,modes,include_coefficients) &
       & result(output) 
      class(DegenerateSymmetry), intent(in) :: this
      type(ComplexMonomial),     intent(in) :: input(:)
      type(ComplexMode),         intent(in) :: modes(:)
      logical,                   intent(in) :: include_coefficients
      type(ComplexMatrix)                   :: output
    end function
  end interface
  
  interface
    module function transform_monomial(this,input,modes) result(output) 
      class(DegenerateSymmetry), intent(in) :: this
      type(ComplexMonomial),     intent(in) :: input
      type(ComplexMode),         intent(in) :: modes(:)
      type(ComplexPolynomial)               :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_SingleModeSymmetry(this,input) 
      class(SingleModeSymmetry), intent(out) :: this
      type(String),              intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_SingleModeSymmetry(this) result(output) 
      class(SingleModeSymmetry), intent(in) :: this
      type(String)                          :: output
    end function
  end interface
  
  interface SingleModeSymmetry
    impure elemental module function new_SingleModeSymmetry_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(SingleModeSymmetry) :: this
    end function
  end interface
  
  interface
    module subroutine read_DegenerateSymmetry(this,input) 
      class(DegenerateSymmetry), intent(out) :: this
      type(String),              intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_DegenerateSymmetry(this) result(output) 
      class(DegenerateSymmetry), intent(in) :: this
      type(String), allocatable             :: output(:)
    end function
  end interface
  
  interface DegenerateSymmetry
    module function new_DegenerateSymmetry_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(DegenerateSymmetry) :: this
    end function
  
    impure elemental module function new_DegenerateSymmetry_StringArray(input) result(this) 
      type(StringArray), intent(in) :: input
      type(DegenerateSymmetry)      :: this
    end function
  end interface
end module
