module generate_kgrid_module
    implicit none
contains

subroutine generate_kgrid(filenames)
  use constants,        only : dp
  use linear_algebra,   only : inv_33
  use file_io,          only : open_read_file, open_write_file
  use structure_module, only : StructureData, read_structure_file, drop
  use string_module
  implicit none
  
  type(String), intent(in) :: filenames(:)
  
  type(String) :: structure_filename
  type(String) :: grid_filename
  type(String) :: ibz_filename
  type(String) :: rotated_gvectors_filename
  
  integer :: grid_file
  integer :: ibz_file
  integer :: rotated_gvectors_file
  
  real(dp),parameter :: tol=1.d-8
  
  type(StructureData) :: structure
  
  integer               :: i_symm
  integer               :: no_gvectors
  integer               :: grid1
  integer               :: grid2
  integer               :: grid3
  integer               :: i_frac
  integer               :: j_frac
  integer               :: k_frac
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
  
  ! Read filenames from input
  structure_filename = filenames(1)
  grid_filename = filenames(2)
  ibz_filename = filenames(3)
  rotated_gvectors_filename = filenames(4)
  
  ! Read in structure file
  structure = read_structure_file(structure_filename)
  
  ! Generate G-vectors 
  grid_file = open_read_file(grid_filename)
  read(grid_file,*)grid1,grid2,grid3
  close(grid_file)
  
  no_gvectors=grid1*grid2*grid3
  allocate(gvecs_cart(3,no_gvectors))
  allocate(ibz(3,no_gvectors))
  allocate(gvecs_frac(3,no_gvectors))
  allocate(rot_operation(no_gvectors))
  allocate(multiplicity(no_gvectors))
  counter=0
  do i_frac=0,grid1-1
    do j_frac=0,grid2-1
      do k_frac=0,grid3-1
        counter=counter+1
        gvecs_frac(1,counter)=dble(i_frac)/dble(grid1)
        gvecs_frac(2,counter)=dble(j_frac)/dble(grid2)
        gvecs_frac(3,counter)=dble(k_frac)/dble(grid3)
      enddo ! k
    enddo ! j
  enddo ! i
  
  do i_frac=1,no_gvectors
    do j_frac=1,3
      if (gvecs_frac(j_frac,i_frac)>0.5d0+tol) then
        gvecs_frac(j_frac,i_frac) = gvecs_frac(j_frac,i_frac)-1.d0
      endif
    enddo
  enddo
  
  do i_vec=1,no_gvectors
    gvecs_cart(:,i_vec) = matmul(gvecs_frac(:,i_vec),structure%recip_lattice)
  enddo
  
  ! Rotate all G-vectors to the IBZ
  ibz=0.d0
  multiplicity=0
  do_i : do i_vec=1,no_gvectors
    if(i_vec==1)then
      ibz(:,i_vec) = gvecs_frac(:,i_vec)
      rot_operation(i_vec)=1
      multiplicity(i_vec)=multiplicity(i_vec)+1
    endif ! i_vec
    do j_vec=1,i_vec-1
      do i_symm=1,structure%no_symmetries
        rvec = matmul( structure%rotation_matrices(:,:,i_symm), &
                     & gvecs_cart(:,i_vec))
        temp_frac = matmul(structure%lattice,rvec)
        temp_frac(1:3)=modulo(temp_frac(1:3)+0.5d0+tol,1.d0)-0.5d0-tol
        if(all(abs(temp_frac(1:3)-gvecs_frac(1:3,j_vec))<tol))then
          ibz(:,i_vec) = gvecs_frac(:,j_vec)
          rot_operation(i_vec)=i_symm
          multiplicity(j_vec)=multiplicity(j_vec)+1
          cycle do_i
        endif ! tol
      enddo ! i_symm
      if(j_vec==i_vec-1)then
        ibz(:,i_vec)=gvecs_frac(:,i_vec)
        rot_operation(i_vec)=1
        multiplicity(i_vec)=multiplicity(i_vec)+1
      endif
    enddo ! j_vec
  enddo do_i
  
  ibz_file = open_write_file(ibz_filename)
  rotated_gvectors_file = open_write_file(rotated_gvectors_filename)
  do i_vec=1,no_gvectors
    if(rot_operation(i_vec)==1)then
      write(ibz_file,*)ibz(:,i_vec),multiplicity(i_vec)
    endif ! rot_operation
    write(rotated_gvectors_file,*)ibz(:,i_vec),rot_operation(i_vec)
  enddo ! i_vec
  close(ibz_file)
  close(rotated_gvectors_file)
  
  call drop(structure)
end subroutine
end module
