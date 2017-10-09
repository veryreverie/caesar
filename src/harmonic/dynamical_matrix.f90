! ======================================================================
! The forces between atoms at a given q-point.
! ======================================================================
module dynamical_matrix_module
  use constants_module, only : dp
  use string_module
  use io_module
  use linear_algebra_module
  
  type :: DynamicalMatrix
    type(ComplexMatrix), allocatable :: matrices(:,:)
  end type
contains

subroutine write_dynamical_matrix_file(dynamical_matrix,filename)
  use ofile_module
  implicit none
  
  type(DynamicalMatrix), intent(in) :: dynamical_matrix
  type(String),          intent(in) :: filename
  
  type(OFile) :: matrix_file
  
  integer :: no_atoms
  
  integer :: i,j
  
  no_atoms = size(dynamical_matrix%matrices,1)
  if (size(dynamical_matrix%matrices,2)/=no_atoms) then
    call err()
  endif
  
  matrix_file = filename
  do i=1,no_atoms
    do j=1,no_atoms
      call matrix_file%print_line('Atoms: '//i//' and '//j//'.')
      call matrix_file%print_line(dynamical_matrix%matrices(j,i))
      call matrix_file%print_line('')
    enddo
  enddo
end subroutine

function read_dynamical_matrix_file(filename) result(dynamical_matrix)
  use utils_module, only : int_sqrt
  use ifile_module
  implicit none
  
  type(String), intent(in) :: filename
  type(DynamicalMatrix)    :: dynamical_matrix
  
  type(IFile) :: matrix_file
  
  integer :: no_atoms
  
  ! Temporary variables
  integer                   :: i,j,k,ialloc
  type(String), allocatable :: line(:)
  complex(dp)               :: matrix(3,3)
  
  matrix_file = filename
  
  no_atoms = int_sqrt(size(matrix_file)/5)
  
  allocate( dynamical_matrix%matrices(no_atoms,no_atoms), &
          & stat=ialloc); call err(ialloc)
  
  do i=1,no_atoms
    do j=1,no_atoms
      do k=1,3
        line = split(matrix_file%line(5*(no_atoms*(i-1)+(j-1))+1+k))
        matrix(k,:) = cmplx(line)
      enddo
      dynamical_matrix%matrices(j,i) = matrix
    enddo
  enddo
end function
end module
