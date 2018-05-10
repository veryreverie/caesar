! ======================================================================
! Electronic band structure information.
! ======================================================================
! N.B. this module is entirely untested.
module bands_module
  use utils_module
  implicit none
  
  private
  
  public :: BandsData
  public :: read_castep_bands_file
  public :: read_vasp_bands_file
  
  type BandsData
    integer               :: no_qpoints
    integer               :: no_bands
    real(dp), allocatable :: bands(:,:)
  end type
  
  interface BandsData
    module procedure new_BandsData
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

function new_BandsData(no_qpoints,no_bands) result(this)
  implicit none
  
  integer, intent(in) :: no_qpoints
  integer, intent(in) :: no_bands
  type(BandsData)     :: this
  
  integer :: ialloc
  
  this%no_qpoints = no_qpoints
  this%no_bands = no_bands
  
  allocate(this%bands(no_bands,no_qpoints), stat=ialloc); call err(ialloc)
end function

! ----------------------------------------------------------------------
! Reads a .bands file into BandsData
! ----------------------------------------------------------------------
! A .bands file looks like:
!
! Number of q-points             ~
! Number of spin components      ~
! Number of electrons            ~
! Number of eigenvalues          ~
! Fermi energy (in atomic units) ~
! Unit cell vectors
!  ~ ~ ~
!  ~ ~ ~
!  ~ ~ ~
! q-point 1 ~ ~ ~ ~
! Spin component 1
!  ~
! ...
!  ~
! q-point 2 ~ ~ ~ ~
! Spin component 1
!  ~
! ...
!  ~
! ...
function read_castep_bands_file_character(filename) result(this)
  implicit none
  
  character(*), intent(in) :: filename
  type(BandsData)          :: this
  
  ! qpoint and band counts
  integer :: no_qpoints
  integer :: no_bands
  
  ! File contents
  type(IFile)               :: bands_file
  type(String), allocatable :: line(:)
  
  ! Temporary variables
  integer        :: i, j
  
  ! Read in bands file
  bands_file = IFile(filename)
  
  ! Get no_qpoints and no_bands
  line = split_line(bands_file%line(1))
  no_qpoints = int(line(4))
  
  line = split_line(bands_file%line(4))
  no_bands = int(line(4))
  
  this = BandsData(no_qpoints,no_bands)
  
  do i=1,no_qpoints
    do j=1,no_bands
      this%bands(j,i) = dble(bands_file%line(11+(no_bands+2)*(i-1)+j))
    enddo
  enddo
end function

function read_castep_bands_file_String(filename) result(this)
  implicit none
  
  type(String), intent(in) :: filename
  type(BandsData)          :: this
  
  this = read_castep_bands_file(char(filename))
end function

function read_vasp_bands_file_character(filename) result(this)
  implicit none
  
  character(*), intent(in) :: filename
  type(BandsData)          :: this
  
  ! qpoint and band counts
  integer :: no_qpoints
  integer :: no_bands
  
  ! File contents
  type(IFile)               :: bands_file
  type(String), allocatable :: line(:)
  
  ! Temporary variables
  integer        :: i, j
  
  ! Read in bands file
  bands_file = IFile(filename)
  
  ! Get no_qpoints and no_bands
  line = split_line(bands_file%line(6))
  no_qpoints = int(line(2))
  no_bands = int(line(3))
  
  this = BandsData(no_qpoints,no_bands)
  
  do i=1,no_qpoints
    do j=1,no_bands
      call print_line("Vasp bands file parser needs updating")
      call err()
      this%bands(j,i) = dble(bands_file%line(11+(no_bands+2)*(i-1)+j))
    enddo
  enddo
end function

function read_vasp_bands_file_String(filename) result(this)
  implicit none
  
  type(String), intent(in) :: filename
  type(BandsData)          :: this
  
  this = read_vasp_bands_file(char(filename))
end function
  
end module
