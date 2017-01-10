module compare_kpoints_module
  implicit none
contains

subroutine compare_kpoints(filenames)
  use constants, only : dp
  use file_io,   only : open_read_file, open_write_file
  use utils,     only : count_lines
  implicit none
  
  character(100), intent(in) :: filenames(:)
  
  real(dp),parameter :: tol=1.d-10
  integer :: i,j,size,label,list
  real(dp) :: kpoint(3),gvec_frac(3)
  integer :: no_lines
  
  ! file units
  integer :: kpoints_file
  integer :: gvectors_frac_file
  integer :: list_file
  
  kpoints_file = open_read_file(filenames(1))
  gvectors_frac_file = open_read_file(filenames(2))
  list_file = open_write_file(filenames(3))

  read(gvectors_frac_file,*) size
  
  no_lines = count_lines(kpoints_file)
  
  do i=1,size
    read(gvectors_frac_file,*) label, gvec_frac(1:3)
    gvec_frac(1:3)=modulo(0.5d0+gvec_frac(1:3)+tol,1.d0)-0.5d0-tol
    do j=1,no_lines
      read(kpoints_file,*)list,kpoint(1:3)
      kpoint(1:3)=modulo(0.5d0+kpoint(1:3)+tol,1.d0)-0.5d0-tol
      if (all(dabs(gvec_frac(1:3)-kpoint(1:3))<tol)) then
        write(list_file,*) list, label
      endif
    enddo
    rewind(kpoints_file)
  enddo

  close(kpoints_file)
  close(gvectors_frac_file)
  close(list_file)

  end subroutine
end module
