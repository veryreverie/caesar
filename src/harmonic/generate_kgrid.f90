! Main program
module generate_kgrid_module
  implicit none
contains

subroutine generate_kgrid()
 USE linear_algebra,ONLY : inv_33

 IMPLICIT NONE

 INTEGER,PARAMETER :: dp=KIND(1.d0)
 REAL(dp),PARAMETER :: tol=1.d-8
 INTEGER :: i_symm,no_symm,no_gvectors,grid1,grid2,grid3,i_frac,j_frac,k_frac,&
   &i_vec,j_vec
 REAL(dp),ALLOCATABLE :: gvecs_cart(:,:),gvecs_frac(:,:),rot_matrix(:,:,:)
 REAL(dp),ALLOCATABLE :: ibz(:,:)
 REAL(dp) :: prim(3,3),recip(3,3),temp_frac(3),rvec(3)
 INTEGER,ALLOCATABLE :: rot_operation(:),multiplicity(:)
 INTEGER :: counter,multiplicity_total

! Read in rotation matrices for the point group
 open(1,file='symmetry.dat')
 read(1,*)no_symm
 allocate(rot_matrix(3,3,no_symm))
 do i_symm=1,no_symm
  read(1,*)rot_matrix(1,1:3,i_symm)
  read(1,*)rot_matrix(2,1:3,i_symm)
  read(1,*)rot_matrix(3,1:3,i_symm)
  read(1,*)
 enddo ! i_symm
 close(1)

! Generate G-vectors 
 open(1,file='grid.dat')
 read(1,*)grid1,grid2,grid3
 close(1)
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

 if(counter/=no_gvectors)then
  write(*,*)'Problem'
  stop
 endif ! counter

 do i_vec=1,no_gvectors
  gvecs_frac(1:3,i_vec)=modulo(gvecs_frac(1:3,i_vec)+0.5d0+tol,1.d0)-0.5d0-tol
 enddo ! i_vec 

! Read in primitive cell vectors
 open(1,file='lattice.dat')
 read(1,*)prim(1,1:3)
 read(1,*)prim(2,1:3)
 read(1,*)prim(3,1:3)
 close(1)
 call inv_33(prim,recip)
 recip=transpose(recip)
 do i_vec=1,no_gvectors
  gvecs_cart(1,i_vec)=dot_product(gvecs_frac(1:3,i_vec),recip(1:3,1))
  gvecs_cart(2,i_vec)=dot_product(gvecs_frac(1:3,i_vec),recip(1:3,2))
  gvecs_cart(3,i_vec)=dot_product(gvecs_frac(1:3,i_vec),recip(1:3,3))
 enddo ! i_vec

! Rotate all G-vectors to the IBZ
 ibz=0.d0
 multiplicity=0
 do i_vec=1,no_gvectors
  counter=0
  if(i_vec==1)then
   ibz(1:3,i_vec)=gvecs_frac(1:3,i_vec)
   rot_operation(i_vec)=1
   multiplicity(i_vec)=multiplicity(i_vec)+1
  endif ! i_vec
  do j_vec=1,i_vec-1
   do i_symm=1,no_symm
    rvec(1)=dot_product(rot_matrix(1,1:3,i_symm),gvecs_cart(1:3,i_vec))
    rvec(2)=dot_product(rot_matrix(2,1:3,i_symm),gvecs_cart(1:3,i_vec))
    rvec(3)=dot_product(rot_matrix(3,1:3,i_symm),gvecs_cart(1:3,i_vec))
    temp_frac(1)=dot_product(prim(1,1:3),rvec(1:3))
    temp_frac(2)=dot_product(prim(2,1:3),rvec(1:3))
    temp_frac(3)=dot_product(prim(3,1:3),rvec(1:3))
    temp_frac(1:3)=modulo(temp_frac(1:3)+0.5d0+tol,1.d0)-0.5d0-tol
    if(all(abs(temp_frac(1:3)-gvecs_frac(1:3,j_vec))<tol))then
     ibz(1:3,i_vec)=gvecs_frac(1:3,j_vec)
     rot_operation(i_vec)=i_symm
     multiplicity(j_vec)=multiplicity(j_vec)+1
     counter=1
     exit
    endif ! tol
   enddo ! i_symm
   if(j_vec==i_vec-1.and.counter==0)then
    ibz(1:3,i_vec)=gvecs_frac(1:3,i_vec)
    rot_operation(i_vec)=1
    multiplicity(i_vec)=multiplicity(i_vec)+1
   endif
   if(counter==1)exit
  enddo ! j_vec
 enddo ! i_vec

 open(3,file='ibz.dat')
 open(4,file='rotated_gvectors.dat')
 multiplicity_total=0
 do i_vec=1,no_gvectors
  if(rot_operation(i_vec)==1)then
   write(3,*)ibz(1:3,i_vec),multiplicity(i_vec)
   multiplicity_total=multiplicity_total+multiplicity(i_vec)
  endif ! rot_operation
  write(4,*)ibz(1:3,i_vec),rot_operation(i_vec)
 enddo ! i_vec
 close(3)
 close(4)
end subroutine
end module
