module kpoints_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  ! N.B. K-points are stored in scaled fractional reciprocal primitive lattice
  !    co-ordinates.
  ! To get K-points in true reciprocal space co-ordinates:
  !    k_true = matmul(recip_lattice,kpoint) / grid
  
  ! ----------------------------------------------------------------------
  ! K-points in a grid across the entire reciprocal unit cell.
  ! ----------------------------------------------------------------------
  ! ibz_ids(i) = the K-point in the IBZ corresponding to K-point i
  !    in the grid, via symmetry operation symmetry_ids(i).
  type KpointsGrid
    integer, allocatable :: kpoints(:,:)
    integer, allocatable :: ibz_ids(:)
    integer, allocatable :: symmetry_ids(:)
  end type
  
  ! ----------------------------------------------------------------------
  ! K-points in the Initial Brillouin Zone (IBZ)
  ! ----------------------------------------------------------------------
  ! multiplicity(i) = the number of K-points in the grid which map onto
  !    K-point i in the IBZ.
  ! sc_ids(i) = the id of the Supercell which has K-point i as a G-vector.
  ! gvector_ids(i) = the id of said G-vector in Supercell sc_ids(i).
  type KpointsIbz
    integer, allocatable :: kpoints(:,:)
    integer, allocatable :: multiplicity(:)
    integer, allocatable :: sc_ids(:)
    integer, allocatable :: gvector_ids(:)
  end type
  
  interface new
    module procedure new_KpointsGrid
    module procedure new_KpointsIbz
  end interface
  
  interface drop
    module procedure drop_KpointsGrid
    module procedure drop_KpointsIbz
  end interface
  
  interface size
    module procedure size_KpointsGrid
    module procedure size_KpointsIbz
  end interface
contains

subroutine new_KpointsGrid(this,no_kpoints)
  implicit none
  
  type(KpointsGrid), intent(out) :: this
  integer,           intent(in)  :: no_kpoints
  
  integer :: ialloc
  
  allocate( this%kpoints(3,no_kpoints),    &
          & this%ibz_ids(no_kpoints),      &
          & this%symmetry_ids(no_kpoints), &
          & stat=ialloc); call err(ialloc)
end subroutine

subroutine new_KpointsIbz(this,no_kpoints)
  implicit none
  
  type(KpointsIbz), intent(out) :: this
  integer,          intent(in)  :: no_kpoints
  
  integer :: ialloc
  
  allocate( this%kpoints(3,no_kpoints),    &
          & this%multiplicity(no_kpoints), &
          & this%sc_ids(no_kpoints),       &
          & this%gvector_ids(no_kpoints),  &
          & stat=ialloc); call err(ialloc)
end subroutine

subroutine drop_KpointsGrid(this)
  implicit none
  
  type(KpointsGrid), intent(inout) :: this
  
  integer :: ialloc
  
  deallocate( this%kpoints,      &
            & this%ibz_ids,      &
            & this%symmetry_ids, &
            & stat=ialloc); call err(ialloc)
end subroutine

subroutine drop_KpointsIbz(this)
  implicit none
  
  type(KpointsIbz), intent(inout) :: this
  
  integer :: ialloc
  
  deallocate( this%kpoints,      &
            & this%multiplicity, &
            & this%sc_ids,       &
            & this%gvector_ids,  &
            & stat=ialloc); call err(ialloc)
end subroutine

function size_KpointsGrid(this) result(output)
  implicit none
  
  type(KpointsGrid), intent(in) :: this
  integer                       :: output
  
  output = size(this%kpoints,2)
end function

function size_KpointsIbz(this) result(output)
  implicit none
  
  type(KpointsIbz), intent(in) :: this
  integer                      :: output
  
  output = size(this%kpoints,2)
end function

function read_kpoints_grid_file(filename) result(this)
  implicit none
  
  type(String), intent(in) :: filename
  type(KpointsGrid)        :: this
  
  integer                   :: i
  type(String), allocatable :: kpoints_grid_file(:)
  type(String), allocatable :: line(:)
  
  kpoints_grid_file = read_lines(filename)
  
  call new(this, size(kpoints_grid_file)-1)
  
  do i=2,size(kpoints_grid_file)
    line = split(kpoints_grid_file(i))
    this%kpoints(:,i) = int(line(1:3))
    this%ibz_ids(i) = int(line(4))
    this%symmetry_ids(i) = int(line(5))
  enddo
end function

function read_kpoints_ibz_file(filename) result(this)
  implicit none
  
  type(String), intent(in) :: filename
  type(KpointsIbz)         :: this
  
  integer                   :: i
  type(String), allocatable :: kpoints_ibz_file(:)
  type(String), allocatable :: line(:)
  
  kpoints_ibz_file = read_lines(filename)
  
  call new(this, size(kpoints_ibz_file)-1)
  
  do i=2,size(kpoints_ibz_file)
    line = split(kpoints_ibz_file(i))
    this%kpoints(:,i) = int(line(1:3))
    this%multiplicity(i) = int(line(4))
    this%sc_ids(i) = int(line(5))
    this%gvector_ids(i) = int(line(6))
  enddo
end function

subroutine write_kpoints_grid_file(this, filename)
  implicit none
  
  type(KpointsGrid), intent(in) :: this
  type(String),      intent(in) :: filename
  
  integer :: i
  integer :: kpoints_grid_file
  
  kpoints_grid_file = open_write_file(filename)
  call print_line(kpoints_grid_file, &
     & 'K-point (x,y,z); Equivalent K-point in IBZ; Rotation to IBZ')
  do i=1,size(this)
    call print_line(kpoints_grid_file, this%kpoints(:,i) //' '// &
                                     & this%ibz_ids(i)   //' '// &
                                     & this%symmetry_ids(i))
  enddo
  close(kpoints_grid_file)
end subroutine

subroutine write_kpoints_ibz_file(this, filename)
  implicit none
  
  type(KpointsIbz), intent(in) :: this
  type(String),     intent(in) :: filename
  
  integer :: i
  integer :: kpoints_ibz_file
  
  kpoints_ibz_file = open_write_file(filename)
  call print_line(kpoints_ibz_file, &
     & 'K-point (x,y,z); No. equivalent K-points; Supercell; G-vector')
  do i=1,size(this)
    call print_line(kpoints_ibz_file, this%kpoints(:,i)    //' '// &
                                    & this%multiplicity(i) //' '// &
                                    & this%sc_ids(i)       //' '// &
                                    & this%gvector_ids(i))
  enddo
  close(kpoints_ibz_file)
end subroutine
end module
