module band_folding_module
  implicit none
contains

! returns the id of the band closest to band_ref
subroutine band_folding()
  use utils, only : i2s
  implicit none
  
  ! Working variables
  integer :: i,iref
  ! Input variables
  real,allocatable :: bands(:)
  real :: band_ref
  integer :: kpoint,no_bands

  open(1,file='input.dat')
  read(1,*)band_ref,kpoint,no_bands
  close(1)
  
  allocate(bands(no_bands))

  open(1,file='kpoint.'//trim(i2s(kpoint))//'.dat')
  do i=1,no_bands
    read(1,*)bands(i)
  enddo ! i
  close(1)

  iref=1
  do i=1,no_bands
    if(abs(bands(i)-band_ref)<=(abs(bands(iref)-band_ref)+0.0002))then
      iref=i
    endif
  enddo ! i
  
  open(1,file='band_number.dat')
  write(1,*)iref
  close(1)
end subroutine
end module
