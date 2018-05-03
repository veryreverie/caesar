! ======================================================================
! q-points of the primitive cell.
! ======================================================================
module qpoint_submodule
  use utils_module
  implicit none
  
  private
  
  public :: QpointData
  public :: operator(==)
  public :: operator(/=)
  public :: write_qpoints_file
  public :: read_qpoints_file
  
  type :: QpointData
    ! The q-point in fractional primitive reciprocal space co-ordinates.
    type(FractionVector) :: qpoint
    
    ! Whether or not 2*q is a G-vector of the primitive cell.
    logical :: is_paired_qpoint
    
    ! The id of q-point q' such that q+q' is a G-vector of the primitive cell.
    integer :: paired_qpoint
  contains
    ! The smallest value of |S| s.t. S.q is a vector of integers,
    !    i.e. such that this q-point is a G-vector of the supercell S.
    procedure, public :: min_sc_size
    
    ! Translates the q-point by a G-vector if necessary, to ensure that all
    !    elements are in [-1/2,1/2), i.e. that the q-point is in the primitive
    !    reciprocal cell.
    procedure, public :: translate_to_primitive
  end type
  
  ! Comparison of q-points.
  interface operator(==)
    module procedure equality_QpointData
  end interface
  
  interface operator(/=)
    module procedure non_equality_QpointData
  end interface
contains

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
! I/O.
! ----------------------------------------------------------------------
subroutine write_qpoints_file(this,filename)
  implicit none
  
  type(QpointData), intent(in), allocatable :: this(:)
  type(String),     intent(in)              :: filename
  
  type(OFile) :: qpoints_file
  
  integer :: i
  
  qpoints_file = OFile(filename)
  do i=1,size(this)
    call qpoints_file%print_line( 'q-point '//i//', q=(qx/nx, qy/ny, qz/nz):')
    call qpoints_file%print_line( this(i)%qpoint)
    call qpoints_file%print_line( 'Is 2q equal to a primitive G-vector?:')
    call qpoints_file%print_line( this(i)%is_paired_qpoint)
    call qpoints_file%print_line( "The ID of q' s.t. q+q' is a primitive &
       &G-vector")
    call qpoints_file%print_line( this(i)%paired_qpoint)
    call qpoints_file%print_line( '')
  enddo
end subroutine

function read_qpoints_file(filename) result(this)
  implicit none
  
  type(String), intent(in)      :: filename
  type(QpointData), allocatable :: this(:)
  
  type(IFile) :: qpoints_file
  
  integer :: no_qpoints
  integer :: qpoint_line
  integer :: i,ialloc
  
  qpoints_file = IFile(filename)
  no_qpoints = size(qpoints_file)/7
  allocate(this(no_qpoints), stat=ialloc); call err(ialloc)
  do i=1,no_qpoints
    qpoint_line = (i-1)*7
    
    this(i)%qpoint           = frac( qpoints_file%split_line(qpoint_line+2) )
    this(i)%is_paired_qpoint = lgcl( qpoints_file%line(      qpoint_line+4) )
    this(i)%paired_qpoint    = int(  qpoints_file%line(      qpoint_line+6) )
  enddo
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
end module
