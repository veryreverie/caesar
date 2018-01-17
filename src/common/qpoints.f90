module qpoints_module
  use constants_module, only : dp
  use string_module
  use io_module
  use linear_algebra_module
  implicit none
  
  ! The large supercell has a supercell matrix Z, given by
  !     | Nx, 0 , 0 |
  ! Z = | 0 , Ny, 0 |
  !     | 0 , 0 , Nz|
  ! Then the q-point is defined by three integers, qx,qy and qz,
  !    where (-Nx)/2 < qx <= Nx/2 etc. (floor division intentional).
  !    - qpoint = (qx/Nx, qy/Ny, qz/Nz)
  !    - scaled_qpoint = Nx * Ny * Nz * qpoint = (qxNyNz, qyNzNx, qzNxNy)
  !    - gvector = Z*qpoint = (qx, qy, qz)
  type QpointData
    ! The q-point in fractional primitive reciprocal space co-ordinates.
    type(RealVector) :: qpoint
    
    ! The q-point in scaled fractional primitive reciprocal co-ordinates.
    type(IntVector) :: scaled_qpoint
    
    ! The gvector of the large supercell corresponding to this q-point.
    type(IntVector) :: gvector
    
    ! Whether or not 2q=G.
    logical :: is_paired_qpoint
    
    ! The id of q-point q' such that q+2q=G.
    integer :: paired_qpoint
    
    ! Whether or not this q-point is to be simulated.
    ! Only one q-point out of all symmetrically equivalent points is chosen.
    ! Only one q-point in each q+q'=G pair is chosen.
    ! At the harmonic level, all other q-points are constructed from
    !    the chosen q-points.
    logical :: to_simulate
    
    ! The supercell matrix, S, with smallest |S| such that
    !    all elements of S.q are integers.
    type(IntMatrix) :: supercell_matrix
    
    ! S.q, where S is supercell_matrix.
    type(IntVector) :: supercell_gvector
  end type
contains

subroutine write_qpoints_file(this,filename)
  use ofile_module
  implicit none
  
  type(QpointData), intent(in), allocatable :: this(:)
  type(String),     intent(in)              :: filename
  
  type(OFile) :: qpoints_file
  
  integer :: i
  
  qpoints_file = filename
  do i=1,size(this)
    call qpoints_file%print_line( 'q-point '//i//', q=(qx/nx, qy/ny, qz/nz):')
    call qpoints_file%print_line( this(i)%qpoint)
    call qpoints_file%print_line( 'Scaled q-point, q*nx*ny*nz:')
    call qpoints_file%print_line( this(i)%scaled_qpoint)
    call qpoints_file%print_line( 'G-vector of large supercell, G=(qx,qy,qz):')
    call qpoints_file%print_line( this(i)%gvector)
    call qpoints_file%print_line( 'Is 2q equal to a primitive G-vector?:')
    call qpoints_file%print_line( this(i)%is_paired_qpoint)
    call qpoints_file%print_line( "The ID of q' s.t. q+q' is a primitive G-vector")
    call qpoints_file%print_line( this(i)%paired_qpoint)
    call qpoints_file%print_line( 'Is this q-point to be simulated?:')
    call qpoints_file%print_line( this(i)%to_simulate)
    call qpoints_file%print_line( 'Supercell Matrix, S, s.t. S.q is integer.:')
    call qpoints_file%print_line( this(i)%supercell_matrix)
    call qpoints_file%print_line( 'Supercell G-vector, G=S.q:')
    call qpoints_file%print_line( this(i)%supercell_gvector)
    call qpoints_file%print_line( '')
  enddo
end subroutine

function read_qpoints_file(filename) result(this)
  use ifile_module
  implicit none
  
  type(String), intent(in)      :: filename
  type(QpointData), allocatable :: this(:)
  
  type(IFile)               :: qpoints_file
  
  integer :: no_qpoints
  integer :: qpoint_line
  
  integer :: supercell_matrix(3,3)
  
  integer :: i,j,ialloc
  
  qpoints_file = filename
  no_qpoints = size(qpoints_file)/19
  
  allocate(this(no_qpoints), stat=ialloc); call err(ialloc)
  
  do i=1,no_qpoints
    qpoint_line = (i-1)*19
    
    this(i)%qpoint            = dble( qpoints_file%split_line(qpoint_line+2 ) )
    this(i)%scaled_qpoint     = int(  qpoints_file%split_line(qpoint_line+4 ) )
    this(i)%gvector           = int(  qpoints_file%split_line(qpoint_line+6 ) )
    this(i)%is_paired_qpoint  = lgcl( qpoints_file%line(      qpoint_line+8 ) )
    this(i)%paired_qpoint     = int(  qpoints_file%line(      qpoint_line+10) )
    this(i)%to_simulate       = lgcl( qpoints_file%line(      qpoint_line+12) )
    do j=1,3
      supercell_matrix(j,:) = int(qpoints_file%split_line(qpoint_line+13+j))
    enddo
    this(i)%supercell_matrix  = supercell_matrix
    this(i)%supercell_gvector = int(  qpoints_file%split_line(qpoint_line+18) )
  enddo
end function
end module
