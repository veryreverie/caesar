module qpoints_module
  use constants_module, only : dp
  use string_module
  use io_module
  use linear_algebra_module
  implicit none
  
  type QpointData
    ! The q-point in fractional primitive reciprocal space co-ordinates.
    type(RealVector)     :: qpoint
    ! The id of the supercell in which this q-point is a G-vector, and the
    !    id of that G-vector in said supercell.
    integer              :: sc_id
    integer              :: gvector_id
    ! The G-vectors in the large supercell which transform onto this q-point,
    !    by symmetry, and the ids of said symmetries.
    integer, allocatable :: gvectors(:)
    integer, allocatable :: symmetry_ids(:)
  end type
  
  interface QpointData
    module procedure new_QpointData
  end interface
contains

function new_QpointData(multiplicity) result(this)
  implicit none
  
  integer, intent(in) :: multiplicity
  type(QpointData)    :: this
  
  integer :: ialloc
  
  allocate( this%gvectors(multiplicity),  &
          & this%symmetry_ids(multiplicity), &
          & stat=ialloc); call err(ialloc)
end function

subroutine write_qpoints_file(this,filename)
  use ofile_module
  implicit none
  
  type(QpointData), intent(in), allocatable :: this(:)
  type(String),     intent(in)              :: filename
  
  type(OFile) :: qpoints_file
  
  integer :: i
  
  qpoints_file = filename
  do i=1,size(this)
    call qpoints_file%print_line( 'q-point:')
    call qpoints_file%print_line( this(i)%qpoint)
    call qpoints_file%print_line( 'Corresponding supercell             : '// &
                                & this(i)%sc_id)
    call qpoints_file%print_line( 'Corresponding G-vector in supercell : '// &
                                & this(i)%gvector_id)
    call qpoints_file%print_line( 'Matching G-vectors in grid          : '// &
                                & this(i)%gvectors)
    call qpoints_file%print_line( 'ID of symmetries to grid G-vector   : '// &
                                & this(i)%symmetry_ids)
    call qpoints_file%print_line( '')
  enddo
end subroutine

function read_qpoints_file(filename) result(this)
  use ifile_module
  implicit none
  
  type(String), intent(in)      :: filename
  type(QpointData), allocatable :: this(:)
  
  type(IFile)               :: qpoints_file
  type(String), allocatable :: line(:)
  
  integer :: no_qpoints
  integer :: ialloc
  integer :: i
  
  qpoints_file = filename
  no_qpoints = size(qpoints_file)/7
  
  allocate(this(no_qpoints), stat=ialloc); call err(ialloc)
  
  do i=1,no_qpoints
    line = split(qpoints_file%line((i-1)*7+2))
    this(i)%qpoint = dble(line)
    
    line = split(qpoints_file%line((i-1)*7+3))
    this(i)%sc_id = int(line(4))
    
    line = split(qpoints_file%line((i-1)*7+4))
    this(i)%gvector_id = int(line(6))
    
    line = split(qpoints_file%line((i-1)*7+5))
    this(i)%gvectors = int(line(6:))
    
    line = split(qpoints_file%line((i-1)*7+6))
    this(i)%symmetry_ids = int(line(8:))
  enddo
end function
end module
