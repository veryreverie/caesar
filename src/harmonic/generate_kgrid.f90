module generate_kgrid_module
    implicit none
contains

subroutine generate_kgrid(structure,grid,ibz_filename, &
   & rotated_gvectors_filename)
  use constants,      only : dp
  use linear_algebra, only : inv_33
  use utils,          only : reduce_interval
  use file_module
  use structure_module
  use string_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  integer,             intent(in) :: grid(3)
  type(String),        intent(in) :: ibz_filename
  type(String),        intent(in) :: rotated_gvectors_filename
  
  ! Parameters
  real(dp),parameter :: tol=1.d-8
  
  ! File units
  integer :: ibz_file
  integer :: rotated_gvectors_file
  
  integer               :: i_symm
  integer               :: no_gvectors
  integer               :: i_vec
  integer               :: j_vec
  real(dp), allocatable :: gvecs_cart(:,:)
  real(dp), allocatable :: gvecs_frac(:,:)
  real(dp), allocatable :: ibz(:,:)
  real(dp)              :: temp_frac(3)
  real(dp)              :: rvec(3)
  integer,  allocatable :: rot_operation(:)
  integer,  allocatable :: multiplicity(:)
  integer               :: counter
  
  ! Temporary variables
  integer :: i,j,k
  
  
  ! Generate G-vectors 
  no_gvectors = grid(1)*grid(2)*grid(3)
  allocate(gvecs_cart(3,no_gvectors))
  allocate(ibz(3,no_gvectors))
  allocate(gvecs_frac(3,no_gvectors))
  allocate(rot_operation(no_gvectors))
  allocate(multiplicity(no_gvectors))
  counter=0
  do i=0,grid(1)-1
    do j=0,grid(2)-1
      do k=0,grid(3)-1
        counter=counter+1
        gvecs_frac(:,counter) = dble((/i,j,k/)) / grid
      enddo
    enddo
  enddo
  
  gvecs_cart = matmul(structure%recip_lattice, gvecs_frac)
  
  ! Rotate all G-vectors to the IBZ
  ibz=0.d0
  multiplicity=0
  do_i : do i_vec=1,no_gvectors
    if(i_vec==1)then
      ibz(:,i_vec) = gvecs_frac(:,i_vec)
      rot_operation(i_vec)=1
      multiplicity(i_vec)=multiplicity(i_vec)+1
    endif
    
    do j_vec=1,i_vec-1
      do i_symm=1,structure%no_symmetries
        rvec = matmul( structure%rotation_matrices(:,:,i_symm), &
                     & gvecs_cart(:,i_vec))
        temp_frac = reduce_interval(matmul(structure%lattice,rvec),tol)
        
        if(all(abs(temp_frac(:)-gvecs_frac(:,j_vec))<tol))then
          ibz(:,i_vec) = gvecs_frac(:,j_vec)
          rot_operation(i_vec)=i_symm
          multiplicity(j_vec)=multiplicity(j_vec)+1
          cycle do_i
        endif
      enddo
      
      if(j_vec==i_vec-1)then
        ibz(:,i_vec)=gvecs_frac(:,i_vec)
        rot_operation(i_vec)=1
        multiplicity(i_vec)=multiplicity(i_vec)+1
      endif
    enddo
  enddo do_i
  
  ibz_file = open_write_file(ibz_filename)
  rotated_gvectors_file = open_write_file(rotated_gvectors_filename)
  do i_vec=1,no_gvectors
    if(rot_operation(i_vec)==1)then
      write(ibz_file,*) ibz(:,i_vec), multiplicity(i_vec)
    endif
    
    write(rotated_gvectors_file,*) ibz(:,i_vec), rot_operation(i_vec)
  enddo
  close(ibz_file)
  close(rotated_gvectors_file)
end subroutine
end module
