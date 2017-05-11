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
    ! The G-vectors in the large supercell which rotate onto this q-point,
    !    and the ids of said rotations.
    integer, allocatable :: gvectors(:)
    integer, allocatable :: rotations(:)
  end type
  
  interface new
    module procedure new_QpointData
  end interface
contains

subroutine new_QpointData(this,multiplicity)
  implicit none
  
  type(QpointData), intent(out) :: this
  integer,          intent(in)  :: multiplicity
  
  integer :: ialloc
  
  allocate( this%gvectors(multiplicity),  &
          & this%rotations(multiplicity), &
          & stat=ialloc); call err(ialloc)
end subroutine

subroutine write_qpoints_file(this,filename)
  implicit none
  
  type(QpointData), intent(in), allocatable :: this(:)
  type(String),     intent(in)              :: filename
  
  integer :: qpoints_file
  integer :: i
  
  qpoints_file = open_write_file(filename)
  do i=1,size(this)
    call print_line(qpoints_file, 'q-point:')
    call print_line(qpoints_file, this(i)%qpoint)
    call print_line(qpoints_file, &
       & 'Corresponding supercell             : '//this(i)%sc_id)
    call print_line(qpoints_file, &
       & 'Corresponding G-vector in supercell : '//this(i)%gvector_id)
    call print_line(qpoints_file, &
       & 'Matching G-vectors in grid       : '//this(i)%gvectors)
    call print_line(qpoints_file, &
       & 'ID of rotations to grid G-vector : '//this(i)%rotations)
    call print_line(qpoints_file, '')
  enddo
  close(qpoints_file)
end subroutine

function read_qpoints_file(filename) result(this)
  implicit none
  
  type(String), intent(in)      :: filename
  type(QpointData), allocatable :: this(:)
  
  type(String), allocatable :: qpoints_file(:)
  type(String), allocatable :: line(:)
  integer :: no_qpoints
  integer :: ialloc
  integer :: i
  
  qpoints_file = read_lines(filename)
  no_qpoints = size(qpoints_file)/7
  
  allocate(this(no_qpoints), stat=ialloc); call err(ialloc)
  
  do i=1,no_qpoints
    line = split(qpoints_file((i-1)*7+2))
    this(i)%qpoint = dble(line)
    
    line = split(qpoints_file((i-1)*7+3))
    this(i)%sc_id = int(line(4))
    
    line = split(qpoints_file((i-1)*7+4))
    this(i)%gvector_id = int(line(6))
    
    line = split(qpoints_file((i-1)*7+5))
    allocate(this(i)%gvectors(size(line)-5), stat=ialloc); call err(ialloc)
    this(i)%gvectors = int(line(6:))
    
    line = split(qpoints_file((i-1)*7+6))
    allocate(this(i)%rotations(size(line)-7), stat=ialloc); call err(ialloc)
    this(i)%rotations = int(line(8:))
  enddo
end function
end module
