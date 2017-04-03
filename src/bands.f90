module bands_module
  use constants, only : dp
  
  type BandsData
    integer               :: no_kpoints
    integer               :: no_bands
    real(dp), allocatable :: bands(:,:)
  end type
  
  interface new
    module procedure new_BandsData
  end interface
  
  interface drop
    module procedure drop_BandsData
  end interface
  
  interface read_castep_bands_file
    module procedure read_castep_bands_file_character
    module procedure read_castep_bands_file_String
  end interface
  
  interface read_vasp_bands_file
    module procedure read_vasp_bands_file_character
    module procedure read_vasp_bands_file_String
  end interface
contains

subroutine new_BandsData(this,no_kpoints,no_bands)
  implicit none
  
  type(BandsData), intent(out) :: this
  integer,         intent(in)  :: no_kpoints
  integer,         intent(in)  :: no_bands
  
  this%no_kpoints = no_kpoints
  this%no_bands = no_bands
  
  allocate(this%bands(no_bands,no_kpoints))
end subroutine

subroutine drop_BandsData(this)
  implicit none
  
  type(BandsData), intent(inout) :: this
  
  deallocate(this%bands)
end subroutine

! ----------------------------------------------------------------------
! Reads a .bands file into BandsData
! ----------------------------------------------------------------------
! A .bands file looks like:
!
! Number of k-points             ~
! Number of spin components      ~
! Number of electrons            ~
! Number of eigenvalues          ~
! Fermi energy (in atomic units) ~
! Unit cell vectors
!  ~ ~ ~
!  ~ ~ ~
!  ~ ~ ~
! K-point 1 ~ ~ ~ ~
! Spin component 1
!  ~
! ...
!  ~
! K-point 2 ~ ~ ~ ~
! Spin component 1
!  ~
! ...
!  ~
! ...
function read_castep_bands_file_character(filename) result(this)
  use string_module
  use file_module
  implicit none
  
  character(*), intent(in) :: filename
  type(BandsData)          :: this
  
  ! kpoint and band counts
  integer :: no_kpoints
  integer :: no_bands
  
  ! File contents
  type(String), allocatable :: bands_file(:)
  type(String), allocatable :: line(:)
  
  ! Temporary variables
  integer        :: i, j
  
  ! Read in bands file
  bands_file = read_lines(filename)
  
  ! Get no_kpoints and no_bands
  line = split(bands_file(1))
  no_kpoints = int(line(4))
  
  line = split(bands_file(4))
  no_bands = int(line(4))
  
  call new(this,no_kpoints,no_bands)
  
  do i=1,no_kpoints
    do j=1,no_bands
      this%bands(j,i) = dble(bands_file(11+(no_bands+2)*(i-1)+j))
    enddo
  enddo
end function

function read_castep_bands_file_String(filename) result(this)
  use string_module
  implicit none
  
  type(String), intent(in) :: filename
  type(BandsData)          :: this
  
  this = read_castep_bands_file(char(filename))
end function

function read_vasp_bands_file_character(filename) result(this)
  use string_module
  use file_module
  implicit none
  
  character(*), intent(in) :: filename
  type(BandsData)          :: this
  
  ! kpoint and band counts
  integer :: no_kpoints
  integer :: no_bands
  
  ! File contents
  type(String), allocatable :: bands_file(:)
  type(String), allocatable :: line(:)
  
  ! Temporary variables
  integer        :: i, j
  
  ! Read in bands file
  bands_file = read_lines(filename)
  
  ! Get no_kpoints and no_bands
  line = split(bands_file(6))
  no_kpoints = int(line(2))
  no_bands = int(line(3))
  
  call new(this,no_kpoints,no_bands)
  
  do i=1,no_kpoints
    do j=1,no_bands
      call print_line("Vasp bands file parser needs updating")
      call err()
      this%bands(j,i) = dble(bands_file(11+(no_bands+2)*(i-1)+j))
    enddo
  enddo
end function

function read_vasp_bands_file_String(filename) result(this)
  use string_module
  implicit none
  
  type(String), intent(in) :: filename
  type(BandsData)          :: this
  
  this = read_vasp_bands_file(char(filename))
end function
  
end module
