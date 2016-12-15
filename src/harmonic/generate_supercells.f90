module generate_supercells_module
  USE constants
  IMPLICIT NONE
CONTAINS

! ----------------------------------------------------------------------
! Calculate the greatest common divisor of two positive integers using
! Euclid's algorithm.
! ----------------------------------------------------------------------
pure function gcd(int_1,int_2) result(output)
  implicit none
  
  integer, intent(in) :: int_1
  integer, intent(in) :: int_2
  integer             :: output
  
  integer :: a
  integer :: b
  integer :: temp
  
  if (a>=b) then
    a = int_1
    b = int_2
  else
    a = int_2
    b = int_1
  endif
  
  do
    temp=mod(a,b)
    if(temp==0)exit
    a=b
    b=temp
  enddo
  
  output=b

end function


 LOGICAL FUNCTION reduce_vec(vecs)
!------------------------------------------------------------------------------!
! Given three linearly independent input vectors a, b and c, construct the     !
! following linear combinations: a+b-c, a-b+c, -a+b+c, a+b+c and  check if any !
! of the four new vectors is shorter than any of a, b or c. If so, replace the !
! longest of a, b and c with the new (shorter) vector. The resulting three     !
! vectors are also linearly independent.                                       !
!------------------------------------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),INTENT(inout) :: vecs(3,3)
 REAL(dp),PARAMETER :: tol_zero=1.d-7
 INTEGER :: longest,i
 REAL(dp) :: newvecs(4,3),maxlen,nlen

! Determine which of the three input vectors is the longest.
 maxlen=0
 DO i=1,3
  nlen=vecs(i,1)**2+vecs(i,2)**2+vecs(i,3)**2
! Test nlen>maxlen within some tolerance to avoid floating point problems.
  IF(nlen-maxlen>tol_zero*maxlen)THEN
   maxlen=nlen
   longest=i
  ENDIF
 ENDDO ! i

! Construct the four linear combinations
 newvecs(1,1:3)=vecs(1,1:3)+vecs(2,1:3)-vecs(3,1:3)
 newvecs(2,1:3)=vecs(1,1:3)-vecs(2,1:3)+vecs(3,1:3)
 newvecs(3,1:3)=-vecs(1,1:3)+vecs(2,1:3)+vecs(3,1:3)
 newvecs(4,1:3)=vecs(1,1:3)+vecs(2,1:3)+vecs(3,1:3)

! Check if any of the four new vectors is shorter than longest of
! input vectors
 reduce_vec=.FALSE.
 DO i=1,4
  nlen=newvecs(i,1)**2+newvecs(i,2)**2+newvecs(i,3)**2
! Test nlen<maxlen within some tolerance to avoid floating point problems.
  IF(nlen-maxlen<-tol_zero*maxlen)THEN
   vecs(longest,1:3)=newvecs(i,1:3)
   reduce_vec=.TRUE.
   EXIT
  ENDIF
 ENDDO ! i

 END FUNCTION reduce_vec


 SUBROUTINE supercells_generator(num_pcells,num_hnf,hnf)
!-------------------------------------------------------------------------!
! Generate all unique supercells that contain a given number of primitive !
! unit cells. See 'Hart and Forcade, Phys. Rev. B 77, 224115 (2008)' for  !
! details of the algorithm.                                               !
!-------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: num_pcells
 INTEGER,INTENT(out) :: num_hnf
 INTEGER,POINTER :: hnf(:,:,:)
 INTEGER :: a,b,c,d,e,f,ialloc,count_hnf,quotient

 count_hnf=0

 do a=1,num_pcells 
  if(.not.mod(num_pcells,a)==0)cycle
  quotient=num_pcells/a
  do c=1,quotient  
   if(.not.mod(quotient,c)==0)cycle
   f=quotient/c
   count_hnf=count_hnf+c*f**2
  enddo ! c
 enddo ! a

 num_hnf=count_hnf
 count_hnf=0

 allocate(hnf(3,3,num_hnf),stat=ialloc)
 if(ialloc/=0)then
  write(*,*)'Problem allocating hnf array in supercells_generator.'
  stop
 endif

 hnf(1:3,1:3,1:num_hnf)=0

 do a=1,num_pcells 
  if(.not.mod(num_pcells,a)==0)cycle
  quotient=num_pcells/a
  do c=1,quotient  
   if(.not.mod(quotient,c)==0)cycle
   f=quotient/c
   do b=0,c-1
    do d=0,f-1
     do e=0,f-1
      count_hnf=count_hnf+1
      hnf(1,1,count_hnf)=a
      hnf(1,2,count_hnf)=b
      hnf(2,2,count_hnf)=c
      hnf(1,3,count_hnf)=d
      hnf(2,3,count_hnf)=e
      hnf(3,3,count_hnf)=f
     enddo ! e
    enddo ! d
   enddo ! b
  enddo ! c
 enddo ! a

 if(count_hnf/=num_hnf)then
  write(*,*)'Did not generate all HNF matrices.'
  stop
 endif 
 
 END SUBROUTINE supercells_generator


 SUBROUTINE minkowski_reduce(vecs)
!-----------------------------------------------------------------------------!
! Given n vectors a(i) that form a basis for a lattice L in n dimensions, the ! 
! a(i) are said to be Minkowski-reduced if the following conditions are met:  !
!                                                                             !
! - a(1) is the shortest non-zero vector in L                                 !
! - for i>1, a(i) is the shortest possible vector in L such that a(i)>=a(i-1) !
!   and the set of vectors a(1) to a(i) are linearly independent              !
!                                                                             !
! In other words the a(i) are the shortest possible basis vectors for L. This !
! routine, given a set of input vectors a'(i) that are possibly not           !
! Minkowski-reduced, returns the vectors a(i) that are.                       !
!-----------------------------------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),INTENT(inout) :: vecs(3,3)
 INTEGER :: i
 REAL(dp) :: tempvec(3,3)
 LOGICAL :: changed

 iter: DO
  tempvec=vecs
  DO i=1,3
! First check linear combinations involving two vectors.
   vecs(i,1:3)=0
   changed=reduce_vec(vecs)
   vecs(i,1:3)=tempvec(i,1:3)
   IF(changed)CYCLE iter
  ENDDO ! i
! Then check linear combinations involving all three.
  IF(reduce_vec(vecs))CYCLE
  EXIT
 ENDDO iter

 END SUBROUTINE minkowski_reduce

! ----------------------------------------
! GENERATE_SUPERCELLS
! ----------------------------------------
subroutine generate_supercells()
  use utils
  use linear_algebra, only : inv_33
  use file_io,        only : open_read_file, open_write_file
  IMPLICIT NONE
  
  REAL(dp),PARAMETER :: tol=1.d-10
  INTEGER,ALLOCATABLE :: multiplicity(:),int_kpoints(:,:),numerator(:,:),&
    &denominator(:,:),super_size(:),label(:)
  INTEGER :: i,j,k,grid(1:3),ialloc,ierr,istat,lcm,num_kpoints,count,&
    &s11,s12,s13,s22,s23,s33,quotient,hnf(3,3),size_count
  REAL(dp),ALLOCATABLE :: kpoints(:,:)
  REAL(dp) :: prim_latt_vecs(3,3),rec_latt_vecs(3,3),temp_latt_vecs(3,3),&
    &temp_scell(3,3),prim(3)
  LOGICAL,ALLOCATABLE :: found_kpoint(:)
  
  ! file units
  integer :: lattice_file   ! lattice.dat
  integer :: grid_file      ! grid.dat
  integer :: ibz_file       ! ibz.dat
  integer :: k_t_s_file     ! kpoint_to_supercell.dat
  integer :: supercell_file ! Supercell_*/supercell.dat
  integer :: size_file      ! Supercell_*/size.dat
  
  ! Get the primitice cell lattice vectors
  lattice_file = open_read_file('lattice.dat')
  read(lattice_file,*)prim_latt_vecs(1,1:3)
  read(lattice_file,*)prim_latt_vecs(2,1:3)
  read(lattice_file,*)prim_latt_vecs(3,1:3)
  close(lattice_file)
  
  ! Get the dimensions of the k-point grid
  grid_file = open_read_file('grid.dat')
  read(grid_file,*)grid(1:3)
  close(grid_file)
  
  call inv_33(prim_latt_vecs,rec_latt_vecs)
  rec_latt_vecs=transpose(rec_latt_vecs)

  ! Get the number of k-points in the ibz.dat file
  ibz_file = open_read_file('ibz.dat')
  num_kpoints = count_lines(ibz_file)
  
  ! Allocate arrays
  allocate( kpoints(3,num_kpoints),     &
          & multiplicity(num_kpoints),  &
          & int_kpoints(3,num_kpoints), &
          & numerator(3,num_kpoints),   &
          & denominator(3,num_kpoints), &
          & super_size(num_kpoints),    &
          & found_kpoint(num_kpoints),  &
          & label(num_kpoints),         &
          & stat=ialloc)
  if(ialloc/=0)then
    write(*,*)'Problem allocating arrays.'
    stop
  endif
  
  ! Read ibz.dat file
  do i=1,num_kpoints
    read(ibz_file,*) kpoints(1:3,i),multiplicity(i)
  enddo
  close(ibz_file)
  
  ! Express k-points as fractions
  do i=1,num_kpoints
    do j=1,3
      int_kpoints(j,i) = nint(kpoints(j,i)*grid(j))
      if (dabs(dble(int_kpoints(j,i))/dble(grid(j))-kpoints(j,i)) >= tol) then
        write(*,*) 'Unable to find fractional representation of k-point.'
        stop
      endif
    enddo
  enddo
  
  numerator(1:3,1:num_kpoints)=0
  denominator(1:3,1:num_kpoints)=1
  
  ! Reduce fractions
  do i=1,num_kpoints
    do j=1,3
      if(int_kpoints(j,i)/=0)then
        numerator(j,i)=int_kpoints(j,i)/gcd(abs(int_kpoints(j,i)),grid(j))
        denominator(j,i)=grid(j)/gcd(abs(int_kpoints(j,i)),grid(j))
      endif ! int_kpoints
    enddo ! j
    lcm = denominator(2,i)*denominator(3,i) &
      & / gcd(denominator(2,i),denominator(3,i))
    lcm = denominator(1,i)*lcm &
      & / gcd(denominator(1,i),lcm)
    super_size(i)=lcm
  enddo ! i
  
  found_kpoint(1:num_kpoints)=.false.
  label(1:num_kpoints)=0
  count=0
  
  k_t_s_file = open_write_file('kpoint_to_supercell.dat')
  ibz_file = open_write_file('ibz.dat')
  
 do size_count=1,maxval(super_size(1:num_kpoints))
  do i=1,num_kpoints
   if(super_size(i)/=size_count)cycle
   if(found_kpoint(i))cycle
   do s11=1,super_size(i)
    if(.not.mod(super_size(i),s11)==0)cycle
    quotient=super_size(i)/s11
    do s22=1,quotient
     if(.not.mod(quotient,s22)==0)cycle
     s33=quotient/s22
     do s12=0,s22-1
      do s13=0,s33-1
       do s23=0,s33-1
        hnf(1:3,1:3)=0
        hnf(1,1)=s11 ; hnf(1,2)=s12 ; hnf(1,3)=s13
        hnf(2,2)=s22 ; hnf(2,3)=s23
        hnf(3,3)=s33
        temp_scell(1:3,1:3)=dble(hnf(1:3,1:3))
        do k=1,3
         prim(k)=sum(temp_scell(k,1:3)*kpoints(1:3,i))
        enddo ! k
        if(all(abs(prim(1:3)-dble(nint(prim(1:3))))<tol))then
         count=count+1
         found_kpoint(i)=.true.
         label(i)=count
         write(k_t_s_file,*)kpoints(1:3,i),label(i)
         write(ibz_file,*)kpoints(1:3,i),multiplicity(i)
         do j=i+1,num_kpoints
          if(found_kpoint(j))cycle
          if(super_size(j)/=super_size(i))cycle
          do k=1,3
           prim(k)=sum(temp_scell(k,1:3)*kpoints(1:3,j))
          enddo ! k
          if(all(abs(prim(1:3)-dble(nint(prim(1:3))))<tol))then
           found_kpoint(j)=.true.
           label(j)=count
           write(k_t_s_file,*)kpoints(1:3,j),label(j)
           write(ibz_file,*)kpoints(1:3,j),multiplicity(j)
          endif ! tol
         enddo ! j
         do k=1,3
          do j=1,3
           temp_latt_vecs(k,j)=sum(dble(hnf(k,1:3))*prim_latt_vecs(1:3,j))
          enddo ! j
         enddo ! k
         call minkowski_reduce(temp_latt_vecs)
         do k=1,3
          do j=1,3
           hnf(k,j)=nint(sum(temp_latt_vecs(k,1:3)*rec_latt_vecs(j,1:3)))  
          enddo ! j
         enddo ! k
         
         call system('mkdir Supercell_'//trim(i2s(count)))
         supercell_file = open_write_file('Supercell_'//trim(i2s(count))//&
          &'/supercell.dat')
         write(supercell_file,*)hnf(1,1:3)
         write(supercell_file,*)hnf(2,1:3)
         write(supercell_file,*)hnf(3,1:3)
         close(supercell_file)
         size_file = open_write_file('Supercell_'//trim(i2s(count))//&
          &'/size.dat')
         write(size_file,*)super_size(i)
         close(size_file)
         
        endif ! tol
        if(found_kpoint(i))exit
       enddo ! s23
       if(found_kpoint(i))exit
      enddo ! s13
      if(found_kpoint(i))exit
     enddo ! s12
     if(found_kpoint(i))exit
    enddo ! s22
    if(found_kpoint(i))exit
   enddo ! s11
  enddo ! i
 enddo ! size_count

  close(ibz_file)
  close(k_t_s_file)

  if(any(.not.found_kpoint(1:num_kpoints)))then
    write(*,*)'Unable to allocate each k-point to a supercell matrix.'
    stop
  endif ! found_kpoint
end subroutine
end module
