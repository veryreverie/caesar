! ======================================================================
! q-points of the primitive cell.
! ======================================================================
module qpoint_module
  use utils_module
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
  
  ! Constructor.
  interface QpointData
    module procedure new_QpointData
    module procedure new_QpointData_Strings
    module procedure new_QpointData_StringArray
  end interface
  
  ! Transformation from a q-point to its pair (-q, modulo G-vectors).
  interface operator(-)
    module procedure pair_QpointData
  end interface
  
  ! Comparison of q-points.
  interface operator(==)
    module procedure equality_QpointData
  end interface
  
  interface operator(/=)
    module procedure non_equality_QpointData
  end interface
contains

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
impure elemental function new_QpointData(qpoint,id,paired_qpoint_id) &
   & result(this)
  implicit none
  
  type(FractionVector), intent(in) :: qpoint
  integer,              intent(in) :: id
  integer,              intent(in) :: paired_qpoint_id
  type(QpointData)                 :: this
  
  this%qpoint           = qpoint
  this%id               = id
  this%paired_qpoint_id = paired_qpoint_id
end function

! ----------------------------------------------------------------------
! Transforms a q-point to its pair.
! ----------------------------------------------------------------------
impure elemental function pair_QpointData(this) result(output)
  implicit none
  
  type(QpointData), intent(in) :: this
  type(QpointData)             :: output
  
  output = QpointData( qpoint           = -this%qpoint,          &
                     & id               = this%paired_qpoint_id, &
                     & paired_qpoint_id = this%id                )
  call translate_to_primitive(output)
end function

! ----------------------------------------------------------------------
! Returns whether or not the q-point is a G-vector of the primitive cell.
! ----------------------------------------------------------------------
impure elemental function is_gvector(this) result(output)
  implicit none
  
  class(QpointData), intent(in) :: this
  logical                       :: output
  
  output = is_int(this%qpoint)
end function

! ----------------------------------------------------------------------
! Returns whether or not 2*q is a G-vector of the primitive cell.
! ----------------------------------------------------------------------
impure elemental function is_paired_qpoint(this) result(output)
  implicit none
  
  class(QpointData), intent(in) :: this
  logical                       :: output
  
  output = this%id==this%paired_qpoint_id
end function

! ----------------------------------------------------------------------
! Returns the smallest integer, s, s.t. there exists the matrix S
!    with determinant |S|=s s.t. S.q is a vector of integers.
! i.e. The supercell with supercell matrix S has this q-point as a G-vector.
! ----------------------------------------------------------------------
function min_sc_size(this) result(output)
  implicit none
  
  class(QpointData), intent(in) :: this
  integer                       :: output
  
  type(IntFraction) :: q(3)
  
  q = frac(this%qpoint)
  output = lcm( q(1)%denominator(), &
              & q(2)%denominator(), &
              & q(3)%denominator()  )
end function

! ----------------------------------------------------------------------
! Translates the q-point by a G-vector if necessary, to ensure that all
!    elements are in [-1/2,1/2), i.e. that the q-point is in the primitive
!    reciprocal cell.
! ----------------------------------------------------------------------
subroutine translate_to_primitive(this)
  implicit none
  
  class(QpointData), intent(inout) :: this
  
  type(IntFraction) :: qpoint(3)
  integer           :: i
  
  ! Convert the q-point to an array of fractions.
  qpoint = frac(this%qpoint)
  ! Translate all elements to [0,1).
  qpoint = modulo(qpoint,1)
  ! Translate all elements in [1/2,1) to [-1/2,0).
  do i=1,3
    if (qpoint(i)>IntFraction(1,2)) then
      qpoint(i) = qpoint(i) - 1
    endif
  enddo
  ! Convert back to a fraction vector.
  this%qpoint = vec(qpoint)
end subroutine

! ----------------------------------------------------------------------
! Comparison of q-points.
! ----------------------------------------------------------------------
! Accounts for translations by G-vectors.
impure elemental function equality_QpointData(this,that) result(output)
  implicit none
  
  type(QpointData), intent(in) :: this
  type(QpointData), intent(in) :: that
  logical                      :: output
  
  output = is_int(this%qpoint-that%qpoint)
end function

impure elemental function non_equality_QpointData(this,that) result(output)
  implicit none
  
  type(QpointData), intent(in) :: this
  type(QpointData), intent(in) :: that
  logical                      :: output
  
  output = .not. this==that
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_QpointData(this,input)
  implicit none
  
  class(QpointData), intent(out) :: this
  type(String),      intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(FractionVector) :: qpoint
  integer              :: id
  integer              :: paired_qpoint_id
  
  select type(this); type is(QpointData)
    if (size(input)/=3) then
      call print_line(ERROR//': Unable to parse q-point from string input:')
      call print_lines(input)
      call err()
    endif
    
    line = split_line(input(1))
    id = int(line(2))
    
    line = split_line(input(2))
    qpoint = frac(line(3:5))
    
    line = split_line(input(3))
    paired_qpoint_id = int(line(11))
    
    this = QpointData(qpoint,id,paired_qpoint_id)
  class default
    call err()
  end select
end subroutine

function write_QpointData(this) result(output)
  implicit none
  
  class(QpointData), intent(in) :: this
  type(String), allocatable     :: output(:)
  
  select type(this); type is(QpointData)
    output = [ 'q-point '//this%id,                                 &
             & 'q = '//this%qpoint,                                 &
             & "The ID of q' s.t. q+q' is a primitive G-vector: "// &
             &    this%paired_qpoint_id                             ]
  class default
    call err()
  end select
end function

function new_QpointData_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(QpointData)         :: this
  
  call this%read(input)
end function

impure elemental function new_QpointData_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(QpointData)              :: this
  
  this = QpointData(str(input))
end function
end module
