module kpoints_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  type KpointData
    ! The k-point in fractional co-ords.
    real(dp)             :: kpoint(3)
    ! The id of the supercell used to simulate this k-point, and the
    !    G-vector in said supercell which matches.
    integer              :: sc_id
    integer              :: gvector_id
    ! The G-vectors in the large supercell which rotate onto this k-point,
    !    and the ids of said rotations.
    integer, allocatable :: gvectors(:)
    integer, allocatable :: rotations(:)
  end type
  
  interface new
    module procedure new_KpointData
  end interface
contains

subroutine new_KpointData(this,multiplicity)
  implicit none
  
  type(KpointData), intent(out) :: this
  integer,          intent(in)  :: multiplicity
  
  integer :: ialloc
  
  allocate( this%gvectors(multiplicity),  &
          & this%rotations(multiplicity), &
          & stat=ialloc); call err(ialloc)
end subroutine

subroutine write_kpoints_file(this,filename)
  implicit none
  
  type(KpointData), intent(in), allocatable :: this(:)
  type(String),     intent(in)              :: filename
  
  integer :: kpoints_file
  integer :: i
  
  kpoints_file = open_write_file(filename)
  do i=1,size(this)
    call print_line(kpoints_file, 'k-point:')
    call print_line(kpoints_file, this(i)%kpoint)
    call print_line(kpoints_file, &
       & 'Corresponding supercell             : '//this(i)%sc_id)
    call print_line(kpoints_file, &
       & 'Corresponding G-vector in supercell : '//this(i)%gvector_id)
    call print_line(kpoints_file, &
       & 'Matching G-vectors in grid       : '//this(i)%gvectors)
    call print_line(kpoints_file, &
       & 'ID of rotations to grid G-vector : '//this(i)%rotations)
    call print_line(kpoints_file, '')
  enddo
  close(kpoints_file)
end subroutine

function read_kpoints_file(filename) result(this)
  implicit none
  
  type(String), intent(in)      :: filename
  type(KpointData), allocatable :: this(:)
  
  type(String), allocatable :: kpoints_file(:)
  type(String), allocatable :: line(:)
  integer :: no_kpoints
  integer :: ialloc
  integer :: i
  
  kpoints_file = read_lines(filename)
  no_kpoints = size(kpoints_file)/7
  
  allocate(this(no_kpoints), stat=ialloc); call err(ialloc)
  
  do i=1,no_kpoints
    line = split(kpoints_file((i-1)*7+2))
    this(i)%kpoint = dble(line)
    
    line = split(kpoints_file((i-1)*7+3))
    this(i)%sc_id = int(line(4))
    
    line = split(kpoints_file((i-1)*7+4))
    this(i)%sc_id = int(line(6))
    
    line = split(kpoints_file((i-1)*7+5))
    allocate(this(i)%gvectors(size(line)-5), stat=ialloc); call err(ialloc)
    this(i)%gvectors = int(line(6:))
    
    line = split(kpoints_file((i-1)*7+6))
    allocate(this(i)%rotations(size(line)-7), stat=ialloc); call err(ialloc)
    this(i)%rotations = int(line(8:))
  enddo
end function
end module
