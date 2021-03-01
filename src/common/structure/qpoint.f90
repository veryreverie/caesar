! ======================================================================
! q-points of the primitive cell.
! ======================================================================
module caesar_qpoint_module
  use caesar_utils_module
  implicit none
  
  private
  
  public :: QpointData
  public :: operator(-)
  public :: operator(==)
  public :: operator(/=)
  
  type, extends(Stringsable) :: QpointData
    ! The q-point in fractional primitive reciprocal space co-ordinates.
    type(FractionVector) :: qpoint
    
    ! The id of this q-point, q, and the q-point q' s.t. q+q' is
    !    a G-vector of the primitive cell.
    integer :: id
    integer :: paired_qpoint_id
  contains
    ! Whether or not q is a G-vector of the primitive cell.
    procedure, public :: is_gvector
    
    ! Whether or not 2*q is a G-vector of the primitive cell.
    procedure, public :: is_paired_qpoint
    
    ! The smallest value of |S| s.t. S.q is a vector of integers,
    !    i.e. such that this q-point is a G-vector of the supercell S.
    procedure, public :: min_sc_size
    
    ! Translates the q-point by a G-vector if necessary, to ensure that all
    !    elements are in [-1/2,1/2), i.e. that the q-point is in the primitive
    !    reciprocal cell.
    procedure, public :: translate_to_primitive
    
    ! I/O.
    procedure, public :: read  => read_QpointData
    procedure, public :: write => write_QpointData
  end type
  
  interface QpointData
    ! ----------------------------------------------------------------------
    ! Constructor.
    ! ----------------------------------------------------------------------
    impure elemental module function new_QpointData(qpoint,id, &
       & paired_qpoint_id) result(this) 
      type(FractionVector), intent(in) :: qpoint
      integer,              intent(in) :: id
      integer,              intent(in) :: paired_qpoint_id
      type(QpointData)                 :: this
    end function
  end interface
  
  interface operator(-)
    ! ----------------------------------------------------------------------
    ! Transforms a q-point to its pair.
    ! ----------------------------------------------------------------------
    impure elemental module function pair_QpointData(this) result(output) 
      type(QpointData), intent(in) :: this
      type(QpointData)             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns whether or not the q-point is a G-vector of the primitive cell.
    ! ----------------------------------------------------------------------
    impure elemental module function is_gvector(this) result(output) 
      class(QpointData), intent(in) :: this
      logical                       :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns whether or not 2*q is a G-vector of the primitive cell.
    ! ----------------------------------------------------------------------
    impure elemental module function is_paired_qpoint(this) result(output) 
      class(QpointData), intent(in) :: this
      logical                       :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the smallest integer, s, s.t. there exists the matrix S
    !    with determinant |S|=s s.t. S.q is a vector of integers.
    ! i.e. The supercell with supercell matrix S has this q-point as a G-vector.
    ! ----------------------------------------------------------------------
    module function min_sc_size(this) result(output) 
      class(QpointData), intent(in) :: this
      integer                       :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Translates the q-point by a G-vector if necessary, to ensure that all
    !    elements are in [-1/2,1/2), i.e. that the q-point is in the primitive
    !    reciprocal cell.
    ! ----------------------------------------------------------------------
    module subroutine translate_to_primitive(this) 
      class(QpointData), intent(inout) :: this
    end subroutine
  end interface
  
  interface operator(==)
    ! ----------------------------------------------------------------------
    ! Comparison of q-points.
    ! ----------------------------------------------------------------------
    ! Accounts for translations by G-vectors.
    impure elemental module function equality_QpointData(this,that) &
       & result(output) 
      type(QpointData), intent(in) :: this
      type(QpointData), intent(in) :: that
      logical                      :: output
    end function
  end interface
  
  interface operator(/=)
    impure elemental module function non_equality_QpointData(this,that) &
       & result(output) 
      type(QpointData), intent(in) :: this
      type(QpointData), intent(in) :: that
      logical                      :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_QpointData(this,input) 
      class(QpointData), intent(out) :: this
      type(String),      intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_QpointData(this) result(output) 
      class(QpointData), intent(in) :: this
      type(String), allocatable     :: output(:)
    end function
  end interface
  
  interface QpointData
    module function new_QpointData_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(QpointData)         :: this
    end function
  
    impure elemental module function new_QpointData_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(QpointData)              :: this
    end function
  end interface
end module
