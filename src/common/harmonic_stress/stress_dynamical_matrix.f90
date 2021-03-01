! ======================================================================
! The harmonic approximation to the stress, in normal-mode co-ordinates.
! ======================================================================
module caesar_stress_dynamical_matrix_module
  use caesar_utils_module
  
  use caesar_structure_module
  use caesar_dynamical_matrices_module
  implicit none
  
  private
  
  public :: StressDynamicalMatrix
  public :: conjg
  
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  
  type, extends(Stringsable) :: StressDynamicalMatrix
    type(DynamicalMatrix), allocatable :: elements(:,:)
  contains
    procedure, public :: expectation => &
                       & expectation_StressDynamicalMatrix
    ! I/O.
    procedure, public :: read  => read_StressDynamicalMatrix
    procedure, public :: write => write_StressDynamicalMatrix
  end type
  
  interface StressDynamicalMatrix
    ! Constructors.
    module function new_StressDynamicalMatrix(elements) result(this) 
      type(DynamicalMatrix), intent(in) :: elements(:,:)
      type(StressDynamicalMatrix)       :: this
    end function
  
    module function new_StressDynamicalMatrix_zeroes(no_atoms) result(this) 
      integer, intent(in)         :: no_atoms
      type(StressDynamicalMatrix) :: this
    end function
  end interface
  
  interface
    ! The expectation of the stress w/r/t a given mode.
    impure elemental module function expectation_StressDynamicalMatrix(this, &
       & mode) result(output) 
      class(StressDynamicalMatrix), intent(in) :: this
      type(ComplexMode),            intent(in) :: mode
      type(RealMatrix)                         :: output
    end function
  end interface
  
  interface operator(+)
    ! Algebra.
    impure elemental module function add_StressDynamicalMatrix_StressDynamicalMatrix(   this,that) result(output) 
      type(StressDynamicalMatrix), intent(in) :: this
      type(StressDynamicalMatrix), intent(in) :: that
      type(StressDynamicalMatrix)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function negative_StressDynamicalMatrix(this) &
       & result(output) 
      type(StressDynamicalMatrix), intent(in) :: this
      type(StressDynamicalMatrix)             :: output
    end function
  
    impure elemental module function                                              subtract_StressDynamicalMatrix_StressDynamicalMatrix(this,that) result(output) 
      type(StressDynamicalMatrix), intent(in) :: this
      type(StressDynamicalMatrix), intent(in) :: that
      type(StressDynamicalMatrix)             :: output
    end function
  end interface
  
  interface operator(*)
    impure elemental module function multiply_StressDynamicalMatrix_real(this,that) result(output) 
      type(StressDynamicalMatrix), intent(in) :: this
      real(dp),                    intent(in) :: that
      type(StressDynamicalMatrix)             :: output
    end function
  
    impure elemental module function multiply_real_StressDynamicalMatrix(this,that) result(output) 
      real(dp),                    intent(in) :: this
      type(StressDynamicalMatrix), intent(in) :: that
      type(StressDynamicalMatrix)             :: output
    end function
  
    impure elemental module function multiply_StressDynamicalMatrix_complex(this,that) result(output) 
      type(StressDynamicalMatrix), intent(in) :: this
      complex(dp),                 intent(in) :: that
      type(StressDynamicalMatrix)             :: output
    end function
  
    impure elemental module function multiply_complex_StressDynamicalMatrix(this,that) result(output) 
      complex(dp),                 intent(in) :: this
      type(StressDynamicalMatrix), intent(in) :: that
      type(StressDynamicalMatrix)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_StressDynamicalMatrix_real(this, &
       & that) result(output) 
      type(StressDynamicalMatrix), intent(in) :: this
      real(dp),                    intent(in) :: that
      type(StressDynamicalMatrix)             :: output
    end function
  
    impure elemental module function divide_StressDynamicalMatrix_complex(this,that) result(output) 
      type(StressDynamicalMatrix), intent(in) :: this
      complex(dp),                 intent(in) :: that
      type(StressDynamicalMatrix)             :: output
    end function
  end interface
  
  interface conjg
    impure elemental module function conjg_StressDynamicalMatrix(input) &
       & result(output) 
      type(StressDynamicalMatrix), intent(in) :: input
      type(StressDynamicalMatrix)             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_StressDynamicalMatrix(this,input) 
      class(StressDynamicalMatrix), intent(out) :: this
      type(String),                 intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_StressDynamicalMatrix(this) result(output) 
      class(StressDynamicalMatrix), intent(in) :: this
      type(String), allocatable                :: output(:)
    end function
  end interface
  
  interface StressDynamicalMatrix
    module function new_StressDynamicalMatrix_Strings(input) result(this) 
      type(String), intent(in)    :: input(:)
      type(StressDynamicalMatrix) :: this
    end function
  
    impure elemental module function new_StressDynamicalMatrix_StringArray(input) result(this) 
      type(StringArray), intent(in) :: input
      type(StressDynamicalMatrix)   :: this
    end function
  end interface
end module
