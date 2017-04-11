module supercell_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  ! A class to hold supercell data.
  ! n.b. recip_supercell and gvectors are stored in a representation where
  ! The primitive IBZ is a cube with side length sc_size.
  type SupercellData
    integer              :: sc_size ! The no. of primitive cells in supercell.
    integer              :: supercell(3,3) ! The supercell lattice vectors.
    integer              :: recip_supercell(3,3) ! The supercell recip. latt.
    integer, allocatable :: gvectors(:,:)
  end type
  
  interface new
    module procedure new_SupercellData
  end interface
  
contains

subroutine new_SupercellData(this,sc_size)
  implicit none
  
  type(SupercellData), intent(out) :: this
  integer,             intent(in)  :: sc_size
  
  this%sc_size = sc_size
  allocate(this%gvectors(3,sc_size))
end subroutine

! ----------------------------------------------------------------------
! Reads a supercells file to an array of type SupercellData
! ----------------------------------------------------------------------
! The supercells file should look like:
!
! Supercell ~
! Size: ~
! Supercell matrix:
!  ~ ~ ~ 
!  ~ ~ ~
!  ~ ~ ~
! G-vectors:
!  ~ ~ ~
!   ...
!  ~ ~ ~
!
! This will be repeated once for each supercell.
function read_supercells_file(filename) result(supercells)
  use linear_algebra_module, only : invert_int
  implicit none
  
  type(String),        intent(in)  :: filename
  type(SupercellData), allocatable :: supercells(:)
  
  integer                   :: i,j
  type(String), allocatable :: supercells_file(:)
  type(String), allocatable :: line(:)
  integer                   :: supercell_id
  logical                   :: reading_supercell
  logical                   :: reading_gvectors
  
  supercells_file = read_lines(filename)
  
  ! Read in total number of supercells.
  supercell_id = 0
  do i=1,size(supercells_file)
    line = split(supercells_file(i))
    if (size(line) == 2) then
      if (line(1) == 'Supercell' .and. line(2) /= 'matrix:') then
        supercell_id = int(line(2))
      endif
    endif
  enddo
  
  ! Allocate supercells.
  allocate(supercells(supercell_id))
  
  ! Read in data.
  supercell_id = 0
  reading_supercell = .False.
  reading_gvectors = .False.
  do i=1,size(supercells_file)
    line = split(supercells_file(i))
    if (size(line) == 2) then
      if (line(1) == 'Supercell' .and. line(2) /= 'matrix:') then
        supercell_id = int(line(2))
      elseif (line(1) == 'Size:') then
        call new(supercells(supercell_id), int(line(2)))
      elseif (line(1) == 'Supercell' .and. line(2) == 'matrix:') then
        j = 1
        reading_supercell = .true.
        reading_gvectors = .false.
      endif
    elseif (size(line) == 1) then
      if (line(1) == 'G-vectors:') then
        j = 1
        reading_supercell = .false.
        reading_gvectors = .true.
      endif
    elseif (size(line) == 3) then
      if (reading_supercell) then
        supercells(supercell_id)%supercell(j,:) = int(line)
        j = j + 1
      elseif (reading_gvectors) then
        supercells(supercell_id)%gvectors(:,j) = int(line)
        j = j + 1
      endif
    endif
  enddo
  
  ! Calculate reciprocal supercells.
  do i=1,size(supercells)
    supercells(i)%recip_supercell = &
       & invert_int(transpose(supercells(i)%supercell))
  enddo
end function

! ----------------------------------------------------------------------
! As above, but writes file instead of reading it.
! ----------------------------------------------------------------------
subroutine write_supercells_file(supercells, filename)
  implicit none
  
  type(SupercellData), intent(in) :: supercells(:)
  type(String),        intent(in) :: filename
  
  integer :: supercells_file
  integer :: i,j
  
  supercells_file = open_write_file(filename)
  do i=1,size(supercells)
    call print_line(supercells_file, 'Supercell '//i)
    call print_line(supercells_file, 'Size: '//supercells(i)%sc_size)
    call print_line(supercells_file, 'Supercell matrix:')
    do j=1,3
      call print_line(supercells_file, supercells(i)%supercell(j,:))
    enddo
    call print_line(supercells_file, 'G-vectors:')
    do j=1,supercells(i)%sc_size
      call print_line(supercells_file, supercells(i)%gvectors(:,j))
    enddo
    call print_line(supercells_file,'')
  enddo
  close(supercells_file)
end subroutine
end module
