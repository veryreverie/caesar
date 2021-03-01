! ======================================================================
! The second derivatives of the potential, in cartesian co-ordinates.
! ======================================================================
module caesar_cartesian_hessian_module
  use caesar_utils_module
  
  use caesar_structure_module
  implicit none
  
  private
  
  public :: CartesianHessian
  
  ! The Hessian is the matrix of force constants, F, such that F.x=f,
  !    where x and f are the cartesian displacement and cartesian force
  !    respectively.
  type, extends(Stringsable) :: CartesianHessian
    type(RealMatrix), allocatable, private :: elements_(:,:)
  contains
    procedure, public :: elements
    
    ! I/O.
    procedure, public :: read  => read_CartesianHessian
    procedure, public :: write => write_CartesianHessian
  end type
  
  interface CartesianHessian
    ! ----------------------------------------------------------------------
    ! Constructs and checks Hessian from given matrices.
    ! ----------------------------------------------------------------------
    ! If check_symmetry is true, then it will be assumed that the Hessian is taken
    !    at the un-displaced co-ordinates.
    module function new_CartesianHessian_elements(supercell,elements, &
       & check_symmetry,logfile) result(this) 
      type(StructureData), intent(in)              :: supercell
      type(RealMatrix),    intent(in)              :: elements(:,:)
      logical,             intent(in)              :: check_symmetry
      type(OFile),         intent(inout), optional :: logfile
      type(CartesianHessian)                       :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Return the force constants between two atoms.
    ! ----------------------------------------------------------------------
    module function elements(this,a,b) result(output) 
      class(CartesianHessian), intent(in) :: this
      type(AtomData),          intent(in) :: a
      type(AtomData),          intent(in) :: b
      type(RealMatrix)                    :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_CartesianHessian(this,input) 
      class(CartesianHessian), intent(out) :: this
      type(String),            intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_CartesianHessian(this) result(output) 
      class(CartesianHessian), intent(in) :: this
      type(String), allocatable           :: output(:)
    end function
  end interface
  
  interface CartesianHessian
    module function new_CartesianHessian_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(CartesianHessian)   :: this
    end function
  end interface
  
  interface CartesianHessian
    impure elemental module function new_CartesianHessian_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(CartesianHessian)        :: this
    end function
  end interface
end module
