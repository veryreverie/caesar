module band_folding_module
  implicit none
contains

! returns the id of the band closest to band_ref
subroutine band_folding(args)
  use constants, only : dp
  use utils,     only : i2s
  use file_module
  use string_module
  implicit none
  
  type(String), intent(in) :: args(:)
  
  ! Working variables
  integer :: i,iref
  real,allocatable :: bands(:)
  real(dp) :: band_ref
  integer :: no_bands
  
  ! File names
  type(String) :: input_file
  type(String) :: output_file
  
  ! File units
  integer :: ifile
  integer :: ofile
  
  input_file = args(1)
  band_ref = dble(args(2))
  output_file = args(3)
  
  no_bands = count_lines(input_file)
  allocate(bands(no_bands))
  
  ifile = open_read_file(input_file)
  do i=1, no_bands
    read(ifile,*) bands(i)
  enddo
  close(ifile)
  
  iref=1
  do i=1,no_bands
    if(abs(bands(i)-band_ref)<=(abs(bands(iref)-band_ref)+0.0002))then
      iref=i
    endif
  enddo ! i
  
  ofile = open_write_file(output_file)
  write(ofile,*) iref
  close(ofile)
end subroutine
end module
