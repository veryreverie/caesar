program band_folding
  implicit none
  ! Working variables
  integer :: i,iref
  character(80) :: filename,c_kp
  ! Input variables
  real,allocatable :: bands(:)
  real :: band_ref
  integer :: kpoint,no_bands

  open(1,file='input.dat')
  read(1,*)band_ref,kpoint,no_bands
  close(1)
  allocate(bands(no_bands))

  write(c_kp,*)kpoint; c_kp=adjustl(c_kp)
  filename=trim('kpoint.')//trim(c_kp)//('.dat')
  open(1,file=filename)
  do i=1,no_bands
    read(1,*)bands(i)
  enddo ! i
  close(1)

  iref=1
  do i=1,no_bands
    !write(*,*)iref,i,bands(i),bands(iref),band_ref
    if(abs(bands(i)-band_ref)<=(abs(bands(iref)-band_ref)+0.0002))then
      !write(*,*)'in:',iref,i,bands(i),bands(iref),band_ref
      iref=i
    endif
  enddo ! i
  open(1,file='band_number.dat')
  write(1,*)iref
  !write(*,*)'The final answer is:',iref
  close(1)

end program band_folding
