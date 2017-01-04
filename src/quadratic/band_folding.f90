module band_folding_module
  implicit none
contains

! returns the id of the band closest to band_ref
subroutine band_folding(args)
  use utils,   only : i2s
  use file_io, only : open_read_file, open_write_file
  implicit none
  
  character(32), intent(in) :: args(:)
  
  ! Working variables
  integer :: i,iref
  real,allocatable :: bands(:)
  real :: band_ref
  integer :: no_bands
  character(32) :: input_file, output_file
  
  ! file units
  integer :: ifile, ofile
  
  input_file = args(1)
  read(args(2),*) band_ref
  read(args(3),*) no_bands
  output_file = args(4)
  
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
