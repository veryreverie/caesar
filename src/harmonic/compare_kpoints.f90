module compare_kpoints_module
  implicit none
contains

! Reads kpoints and gvectors, and writes a list of all instances when the two 
! are equal
subroutine compare_kpoints(filenames)
  use constants, only : dp
  use file_io,   only : open_read_file, open_write_file, count_lines
  implicit none
  
  character(100), intent(in) :: filenames(:)
  
  real(dp), parameter :: tol = 1.d-10
  
  ! kpoint variables
  integer               :: no_kpoints
  integer,  allocatable :: list(:)
  real(dp), allocatable :: kpoint(:,:)
  
  ! gvector variables
  integer               :: no_gvectors
  integer,  allocatable :: label(:)
  real(dp), allocatable :: gvec_frac(:,:)
  
  ! loop indices
  integer :: i,j
  
  ! file units
  integer :: kpoints_file
  integer :: gvectors_frac_file
  integer :: list_file
  
  ! read kpoints_file
  kpoints_file = open_read_file(filenames(1))
  no_kpoints = count_lines(kpoints_file)
  allocate(list(no_kpoints))
  allocate(kpoint(3,no_kpoints))
  do i=1,no_kpoints
    read(kpoints_file,*) list(i), kpoint(:,i)
  enddo
  close(kpoints_file)
  
  do i=1,no_kpoints
    do j=1,3
      if (kpoint(j,i)>0.5d0+tol) then
        kpoint(j,i) = kpoint(j,i)-1.d0
      endif
    enddo
  enddo
  
  ! read gvectors_frac_file
  gvectors_frac_file = open_read_file(filenames(2))
  read(gvectors_frac_file,*) no_gvectors
  allocate(label(no_gvectors))
  allocate(gvec_frac(3,no_gvectors))
  do i=1,no_gvectors
    read(gvectors_frac_file,*) label(i), gvec_frac(:,i)
  enddo
  close(gvectors_frac_file)
  
  do i=1,no_gvectors
    do j=1,3
      if (gvec_frac(j,i)>0.5d0+tol) then
        gvec_frac(j,i) = gvec_frac(j,i)-1.d0
      endif
    enddo
  enddo
  
  list_file = open_write_file(filenames(3))
  do i=1,no_gvectors
    do j=1,no_kpoints
      if (all(abs(gvec_frac(:,i)-kpoint(:,j))<tol)) then
        write(list_file,*) list(j), label(i)
      endif
    enddo
  enddo
  close(list_file)
  
  deallocate(list)
  deallocate(kpoint)
  deallocate(label)
  deallocate(gvec_frac)

  end subroutine
end module
