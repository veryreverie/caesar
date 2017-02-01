module compare_kpoints_module
  implicit none
contains

! Reads kpoints and gvectors, and writes a list of all instances when the two 
! are equal
subroutine compare_kpoints(args)
  use constants, only : dp
  use utils,     only : reduce_interval
  use file_module
  use string_module
  implicit none
  
  type(String), intent(in) :: args(:)
  
  real(dp), parameter :: tol = 1.d-10
  
  ! Input variables
  integer :: sc_id ! supercell id
  
  ! k-point variables
  integer               :: no_kpoints
  
  real(dp), allocatable :: kpoints(:,:)
  integer,  allocatable :: sc_ids(:)
  
  ! g-vector variables
  integer               :: no_gvectors
  integer,  allocatable :: gvector_ids(:)
  real(dp), allocatable :: gvec_frac(:,:)
  
  ! File names
  type(String) :: ibz_filename
  type(String) :: gvectors_frac_filename
  type(String) :: list_filename
  
  ! File units
  integer :: ibz_file
  integer :: gvectors_frac_file
  integer :: list_file
  
  ! Temporary variables
  integer        :: i,j
  character(100) :: dump
  
  sc_id = int(args(1))
  ibz_filename = args(2)
  gvectors_frac_filename = args(3)
  list_filename = args(4)
  
  ! Read kpoints from ibz.dat
  no_kpoints = count_lines(ibz_filename)
  allocate(kpoints(3,no_kpoints))
  allocate(sc_ids(no_kpoints))
  ibz_file = open_read_file(ibz_filename)
  do i=1,no_kpoints
    read(ibz_file,*) kpoints(:,i),dump,sc_ids(i)
  enddo
  
  ! read gvectors_frac_file
  gvectors_frac_file = open_read_file(gvectors_frac_filename)
  read(gvectors_frac_file,*) no_gvectors
  allocate(gvector_ids(no_gvectors))
  allocate(gvec_frac(3,no_gvectors))
  do i=1,no_gvectors
    read(gvectors_frac_file,*) gvector_ids(i), gvec_frac(:,i)
  enddo
  close(gvectors_frac_file)
  
  list_file = open_append_file(list_filename)
  do i=1,no_gvectors
    do j=1,no_kpoints
      if (sc_ids(j)/=sc_id) cycle ! skip kpoints not in this unit cell
      if (all(abs(gvec_frac(:,i)-kpoints(:,j))<tol)) then
        write(list_file,*) j, gvector_ids(i), sc_id
      endif
    enddo
  enddo
  close(list_file)
end subroutine
end module
