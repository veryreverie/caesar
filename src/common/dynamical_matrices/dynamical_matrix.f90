! ======================================================================
! The forces between atoms at a given q-point.
! ======================================================================
module caesar_dynamical_matrix_module
  use caesar_utils_module
  
  use caesar_structure_module
  use caesar_normal_mode_module
  
  use caesar_min_images_module
  implicit none
  
  private
  
  public :: DynamicalMatrix
  public :: reconstruct_hessian
  public :: ComplexMode
  public :: conjg
  
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  
  type, extends(Stringsable) :: DynamicalMatrix
    type(ComplexMatrix), allocatable, private :: elements_(:,:)
  contains
    procedure, public :: check
    
    procedure, public :: elements => elements_DynamicalMatrix
    
    procedure, public :: expectation => expectation_DynamicalMatrix
    
    ! I/O.
    procedure, public :: read  => read_DynamicalMatrix
    procedure, public :: write => write_DynamicalMatrix
  end type
  
  interface DynamicalMatrix
    ! Constructors.
    module function new_DynamicalMatrix(elements) result(this) 
      type(ComplexMatrix), intent(in) :: elements(:,:)
      type(DynamicalMatrix)           :: this
    end function
  
    module function new_DynamicalMatrix_zeroes(no_atoms) result(this) 
      integer, intent(in)   :: no_atoms
      type(DynamicalMatrix) :: this
    end function
  end interface
  
  interface
    ! Getter for elements.
    module function elements_DynamicalMatrix(this) result(output) 
      class(DynamicalMatrix), intent(in) :: this
      type(ComplexMatrix), allocatable   :: output(:,:)
    end function
  end interface
  
  interface
    ! The expectation of the dynamical matrix w/r/t a given mode.
    impure elemental module function expectation_DynamicalMatrix(this,mode) &
       & result(output) 
      class(DynamicalMatrix), intent(in) :: this
      type(ComplexMode),      intent(in) :: mode
      real(dp)                           :: output
    end function
  end interface
  
  interface DynamicalMatrix
    ! ----------------------------------------------------------------------
    ! Constructs the dynamical matrix, which is the matrix of force constants in
    !    q-point co-ordinates,
    !    given the Hessian, which is the matrix of force constants in
    !    cartesian supercell co-ordinates.
    ! ----------------------------------------------------------------------
    
    ! --------------------------------------------------
    ! Construct a dynamical matrix at a q-point which is not commensurate
    !    with the supercell.
    ! This is only an approximation, using a minimum-image convention.
    ! --------------------------------------------------
    module function new_DynamicalMatrix_interpolated(q,supercell,hessian, &
       & min_images) result(this) 
      type(RealVector),       intent(in)           :: q
      type(StructureData),    intent(in)           :: supercell
      type(CartesianHessian), intent(in)           :: hessian
      type(MinImages),        intent(in), optional :: min_images(:,:)
      type(DynamicalMatrix)                        :: this
    end function
  
    module function new_DynamicalMatrix_qpoint(qpoint,supercell,hessian, &
       & min_images) result(this) 
      type(QpointData),       intent(in)           :: qpoint
      type(StructureData),    intent(in)           :: supercell
      type(CartesianHessian), intent(in)           :: hessian
      type(MinImages),        intent(in), optional :: min_images(:,:)
      type(DynamicalMatrix)                        :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Check a dynamical matrix.
    ! ----------------------------------------------------------------------
    ! Always checks that the matrix is Hermitian.
    ! If check_eigenstuff is .true., also checks that the normal modes match
    !    the dynamical matrix.
    ! check_eigenstuff defaults to .true..
    ! Structure may be any supercell.
    module subroutine check(this,structure,logfile) 
      class(DynamicalMatrix), intent(in)    :: this
      type(StructureData),    intent(in)    :: structure
      type(OFile),            intent(inout) :: logfile
    end subroutine
  end interface
  
  interface ComplexMode
    ! ----------------------------------------------------------------------
    ! Calculates complex modes by diagonalising a dynamical matrix.
    ! ----------------------------------------------------------------------
    ! N.B. Structure may be any supercell.
    
    module function new_ComplexMode_DynamicalMatrix(dynamical_matrix, &
       & structure,modes_real) result(output) 
      type(DynamicalMatrix), intent(in)           :: dynamical_matrix
      type(StructureData),   intent(in)           :: structure
      logical,               intent(in), optional :: modes_real
      type(ComplexMode), allocatable              :: output(:)
    end function
  end interface
  
  interface calculate_modes
    ! Calculate modes at a q-point where 2q/=G.
    module function calculate_modes_complex(elements,structure) result(output) 
      type(ComplexMatrix), intent(in) :: elements(:,:)
      type(StructureData), intent(in) :: structure
      type(ComplexMode), allocatable  :: output(:)
    end function
  
    ! Calculate modes at a q-point where 2q=G.
    module function calculate_modes_real(elements,structure) result(output) 
      type(RealMatrix),    intent(in) :: elements(:,:)
      type(StructureData), intent(in) :: structure
      type(ComplexMode), allocatable  :: output(:)
    end function
  end interface
  
  interface DynamicalMatrix
    ! ----------------------------------------------------------------------
    ! Construct a DynamicalMatrix from the ComplexModes at a given q-point.
    ! ----------------------------------------------------------------------
    module function new_DynamicalMatrix_ComplexModes(modes,frequencies) &
       & result(this) 
      type(ComplexMode), intent(in)           :: modes(:)
      real(dp),          intent(in), optional :: frequencies(:)
      type(DynamicalMatrix)                   :: this
    end function
  
    ! ----------------------------------------------------------------------
    ! Construct the contribution to a DynamicalMatrix two ComplexModes.
    ! If they are at the same q-point, the conjugate of the second mode is taken.
    ! If they are at paired q-points, no conjugate is taken.
    ! If they are neither at the same q-point nor at paired q-points,
    !    an error is thrown.
    ! ----------------------------------------------------------------------
    module function new_DynamicalMatrix_ComplexMode_ComplexMode(mode1,mode2, &
       & coefficient) result(this) 
      type(ComplexMode), intent(in)           :: mode1
      type(ComplexMode), intent(in)           :: mode2
      real(dp),          intent(in), optional :: coefficient
      type(DynamicalMatrix)                   :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Construct the Hessian for the large supercell, given the
    !    dynamical matrices at each q-point.
    ! ----------------------------------------------------------------------
    module function reconstruct_hessian(large_supercell,qpoints, &
       & dynamical_matrices,logfile) result(output) 
      type(StructureData),   intent(in)              :: large_supercell
      type(QpointData),      intent(in)              :: qpoints(:)
      type(DynamicalMatrix), intent(in)              :: dynamical_matrices(:)
      type(OFile),           intent(inout), optional :: logfile
      type(CartesianHessian)                         :: output
    end function
  end interface
  
  interface operator(+)
    ! ----------------------------------------------------------------------
    ! Algebra involving dynamical matrices.
    ! ----------------------------------------------------------------------
    impure elemental module function add_DynamicalMatrix_DynamicalMatrix(this,that) result(output) 
      type(DynamicalMatrix), intent(in) :: this
      type(DynamicalMatrix), intent(in) :: that
      type(DynamicalMatrix)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function negative_DynamicalMatrix(this) &
       & result(output) 
      type(DynamicalMatrix), intent(in) :: this
      type(DynamicalMatrix)             :: output
    end function
  
    impure elemental module function subtract_DynamicalMatrix_DynamicalMatrix(this,that) result(output) 
      type(DynamicalMatrix), intent(in) :: this
      type(DynamicalMatrix), intent(in) :: that
      type(DynamicalMatrix)             :: output
    end function
  end interface
  
  interface operator(*)
    impure elemental module function multiply_DynamicalMatrix_real(this,that) &
       & result(output) 
      type(DynamicalMatrix), intent(in) :: this
      real(dp),              intent(in) :: that
      type(DynamicalMatrix)             :: output
    end function
  
    impure elemental module function multiply_real_DynamicalMatrix(this,that) &
       & result(output) 
      real(dp),              intent(in) :: this
      type(DynamicalMatrix), intent(in) :: that
      type(DynamicalMatrix)             :: output
    end function
  
    impure elemental module function multiply_DynamicalMatrix_complex(this, &
       & that) result(output) 
      type(DynamicalMatrix), intent(in) :: this
      complex(dp),           intent(in) :: that
      type(DynamicalMatrix)             :: output
    end function
  
    impure elemental module function multiply_complex_DynamicalMatrix(this, &
       & that) result(output) 
      complex(dp),           intent(in) :: this
      type(DynamicalMatrix), intent(in) :: that
      type(DynamicalMatrix)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_DynamicalMatrix_real(this,that) &
       & result(output) 
      type(DynamicalMatrix), intent(in) :: this
      real(dp),              intent(in) :: that
      type(DynamicalMatrix)             :: output
    end function
  
    impure elemental module function divide_DynamicalMatrix_complex(this, &
       & that) result(output) 
      type(DynamicalMatrix), intent(in) :: this
      complex(dp),           intent(in) :: that
      type(DynamicalMatrix)             :: output
    end function
  end interface
  
  interface conjg
    impure elemental module function conjg_DynamicalMatrix(input) &
       & result(output) 
      type(DynamicalMatrix), intent(in) :: input
      type(DynamicalMatrix)             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_DynamicalMatrix(this,input) 
      class(DynamicalMatrix), intent(out) :: this
      type(String),           intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_DynamicalMatrix(this) result(output) 
      class(DynamicalMatrix), intent(in) :: this
      type(String), allocatable          :: output(:)
    end function
  end interface
  
  interface DynamicalMatrix
    module function new_DynamicalMatrix_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(DynamicalMatrix)    :: this
    end function
  
    impure elemental module function new_DynamicalMatrix_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(DynamicalMatrix)         :: this
    end function
  end interface
end module
