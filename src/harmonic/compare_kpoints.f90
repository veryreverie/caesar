module compare_kpoints_module
  implicit none
contains

! Reads kpoints and gvectors, and writes a list of all instances when the two 
! are equal
subroutine compare_kpoints(filenames)
  use constants, only : dp
  use file_io,   only : open_read_file, open_write_file, count_lines
  use string_module
  implicit none
  
  type(String), intent(in) :: filenames(:)
  
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
  
  ! file names
  type(String) :: kpoints_filename
  type(String) :: gvectors_frac_filename
  type(String) :: list_filename
  
  ! file units
  integer :: kpoints_file
  integer :: gvectors_frac_file
  integer :: list_file
  
  kpoints_filename = filenames(1)
  gvectors_frac_filename = filenames(2)
  list_filename = filenames(3)
  
  ! read kpoints_file
  no_kpoints = count_lines(kpoints_filename)
  allocate(list(no_kpoints))
  allocate(kpoint(3,no_kpoints))
  kpoints_file = open_read_file(kpoints_filename)
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
  gvectors_frac_file = open_read_file(gvectors_frac_filename)
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
  
  list_file = open_write_file(list_filename)
  do i=1,no_gvectors
    do j=1,no_kpoints
      if (all(abs(gvec_frac(:,i)-kpoint(:,j))<tol)) then
        write(list_file,*) list(j), label(i)
      endif
    enddo
  enddo
  close(list_file)
  end subroutine
end module
