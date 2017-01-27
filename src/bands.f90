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
  use file_module
  implicit none
  
  character(*), intent(in) :: filename
  type(BandsData)          :: this
  
  ! kpoint and band data
  integer :: no_kpoints
  integer :: no_bands
  integer :: kpoint
  integer :: band
  
  ! Line numbers
  integer :: bands_file_length
  
  ! File unit
  integer :: bands_file
  
  ! Temporary variables
  integer        :: i
  character(100) :: line
  character(100) :: dump
  
  bands_file_length = count_lines(filename)
  bands_file = open_read_file(filename)
  
  ! Get no_kpoints and no_bands
  do i=1,bands_file_length
    read(bands_file,"(a)") line
    if (i==1) then
      read(line,*) dump,dump,dump,no_kpoints
    elseif (i==4) then
      read(line,*) dump,dump,dump,no_bands
    endif
  enddo
  
  call new(this,no_kpoints,no_bands)
  
  rewind(bands_file)
  
  do i=1,bands_file_length
    read(bands_file,"(a)") line
    kpoint = (i-10)/(no_bands+2)+1
    band = modulo(i-10,no_bands+2)-1
    if (kpoint > 0 .and. band > 0) then
      read(line,*) this%bands(band,kpoint)
    endif
  enddo
  
  close(bands_file)
end function

function read_castep_bands_file_String(filename) result(this)
  use string_module
  implicit none
  
  type(String), intent(in) :: filename
  type(BandsData)          :: this
  
  this = read_castep_bands_file(char(filename))
end function

function read_vasp_bands_file_character(filename) result(this)
  use file_module
  implicit none
  
  character(*), intent(in) :: filename
  type(BandsData)          :: this
  
  ! kpoint and band data
  integer :: no_kpoints
  integer :: no_bands
  integer :: kpoint
  integer :: band
  
  ! Line numbers
  integer :: bands_file_length
  
  ! File unit
  integer :: bands_file
  
  ! Temporary variables
  integer        :: i
  character(100) :: line
  character(100) :: dump
  
  bands_file_length = count_lines(filename)
  bands_file = open_read_file(filename)
  
  ! Get no_kpoints and no_bands
  do i=1,bands_file_length
    read(bands_file,"(a)") line
    if (i==6) then
      read(line,*) dump,no_kpoints,no_bands
    endif
  enddo
  
  call new(this,no_kpoints,no_bands)
  
  rewind(bands_file)
  
  do i=1,bands_file_length
    read(bands_file,"(a)") line
    kpoint = (i-7)/(no_bands+2)+1
    band = modulo(i-7,no_bands+2)-1
    if (kpoint > 0 .and. band > 0) then
      read(line,*) dump, this%bands(band,kpoint)
    endif
  enddo
  
  close(bands_file)
end function

function read_vasp_bands_file_String(filename) result(this)
  use string_module
  implicit none
  
  type(String), intent(in) :: filename
  type(BandsData)          :: this
  
  this = read_vasp_bands_file(char(filename))
end function
  
end module
