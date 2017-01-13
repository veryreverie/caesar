module generate_supercells_module
  use constants
  implicit none

interface lcm
  module procedure lcm_2 ! lowest common multiple of two positive integers
  module procedure lcm_3 ! lowest common multiple of three positive integers
end interface

contains

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
  
  if (int_1>=int_2) then
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

! ----------------------------------------------------------------------
! Calculate the lowest common multiple of two positive integers
! lcm(a,b) = a*b/gcd(a,b)
! ----------------------------------------------------------------------
pure function lcm_2(int_1,int_2) result(output)
  implicit none
  
  integer, intent(in) :: int_1
  integer, intent(in) :: int_2
  integer             :: output
  
  output = int_1*int_2/gcd(int_1,int_2)
end function

! ----------------------------------------------------------------------
! Calculate the lowest common multiple of three positive integers
! lcm(a,b,c) = lcm(a,lcm(b,c))
! ----------------------------------------------------------------------
pure function lcm_3(int_1,int_2,int_3) result(output)
  implicit none
  
  integer, intent(in) :: int_1
  integer, intent(in) :: int_2
  integer, intent(in) :: int_3
  integer             :: output
  
  integer :: lcm_2
  
  lcm_2 = lcm(int_2,int_3)
  output = lcm(int_1,lcm_2)
end function

!------------------------------------------------------------------------------
! Given three linearly independent input vectors a, b and c, construct the
! following linear combinations: a+b-c, a-b+c, -a+b+c, a+b+c and  check if any
! of the four new vectors is shorter than any of a, b or c. If so, replace the
! longest of a, b and c with the new (shorter) vector. The resulting three
! vectors are also linearly independent.
!------------------------------------------------------------------------------
logical function reduce_vec(vecs)
  implicit none
  
  real(dp), intent(inout) :: vecs(3,3)
  
  integer  :: longest
  integer  :: i
  real(dp) :: newvecs(4,3)
  real(dp) :: maxlen
  real(dp) :: nlen
  
  ! Determine which of the three input vectors is the longest.
  maxlen=0
  do i=1,3
    nlen = norm2(vecs(i,:))
    if (nlen>maxlen) then
      maxlen=nlen
      longest=i
    endif
  enddo
  
  ! Construct the four linear combinations
  newvecs(1,:) =  vecs(1,:)+vecs(2,:)-vecs(3,:)
  newvecs(2,:) =  vecs(1,:)-vecs(2,:)+vecs(3,:)
  newvecs(3,:) = -vecs(1,:)+vecs(2,:)+vecs(3,:)
  newvecs(4,:) =  vecs(1,:)+vecs(2,:)+vecs(3,:)
  
  ! Check if any of the four new vectors is shorter than longest of
  ! input vectors
  reduce_vec=.false.
  do i=1,4
    nlen = norm2(newvecs(i,:))
    if(nlen<maxlen)then
      vecs(longest,:) = newvecs(i,:)
      reduce_vec = .true.
      exit
    endif
  enddo
end function


!-------------------------------------------------------------------------!
! Generate all unique supercells that contain a given number of primitive !
! unit cells. See 'Hart and Forcade, Phys. Rev. B 77, 224115 (2008)' for  !
! details of the algorithm.                                               !
!-------------------------------------------------------------------------!
subroutine supercells_generator(num_pcells,num_hnf,hnf)
 implicit none
 
 integer,intent(in) :: num_pcells
 integer,intent(out) :: num_hnf
 
 integer,pointer :: hnf(:,:,:)
 integer :: a,b,c,d,e,f,ialloc,count_hnf,quotient

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
end subroutine


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
subroutine minkowski_reduce(vecs)
  implicit none
  
  real(dp),intent(inout) :: vecs(3,3)
  
  integer  :: i
  real(dp) :: tempvec(3)
  logical  :: changed
  
  iter: do
    ! First check linear combinations involving two vectors.
    do i=1,3
      tempvec = vecs(i,:)
      vecs(i,:) = 0
      changed = reduce_vec(vecs)
      vecs(i,:) = tempvec
      if(changed)cycle iter
    enddo
    
    ! Then check linear combinations involving all three.
    if(reduce_vec(vecs))cycle
    
    exit
  enddo iter
end subroutine minkowski_reduce

subroutine generate_supercells(args)
  use constants,        only : dp
  use utils,            only : i2s
  use linear_algebra,   only : inv_33
  use file_io,          only : open_read_file, open_write_file, count_lines
  
  use string_module
  use structure_module
  implicit none
  
  ! inputs
  type(String), intent(in) :: args(:)
  
  ! parameters
  real(dp),parameter :: tol=1.d-10
  
  type(StructureData) :: structure
  
  integer,allocatable :: multiplicity(:),int_kpoints(:,:),numerator(:,:),&
    &denominator(:,:),super_size(:),label(:)
  integer :: i,j,k,grid(1:3),ialloc,num_kpoints,count,&
    &s11,s12,s13,s22,s23,s33,quotient,hnf(3,3),size_count
  real(dp),allocatable :: kpoints(:,:)
  real(dp) :: temp_latt_vecs(3,3)
  real(dp) :: temp_scell(3,3)
  real(dp) :: prim(3)
  logical,allocatable :: found_kpoint(:)
  character(100) :: sdir ! Supercell directory
  
  ! file names
  character(100) :: structure_filename
  character(100) :: grid_filename
  character(100) :: ibz_filename
  character(100) :: kpoint_to_supercell_filename
  character(100) :: supercell_directory_root
  
  ! file units
  integer :: grid_file
  integer :: ibz_file
  integer :: k_t_s_file
  integer :: supercell_file
  integer :: size_file
  
  ! read filenames from input
  structure_filename = args(1)
  grid_filename = args(2)
  ibz_filename = args(3)
  kpoint_to_supercell_filename = args(4)
  supercell_directory_root = args(5)
  
  ! Read the structure file
  structure = read_structure_file(structure_filename)
  
  ! Get the dimensions of the k-point grid
  grid_file = open_read_file(grid_filename)
  read(grid_file,*)grid(1:3)
  close(grid_file)

  ! Get the number of k-points in the ibz.dat file
  ibz_file = open_read_file(ibz_filename)
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
  
  numerator = 0
  denominator = 1
  
  ! Reduce fractions
  do i=1,num_kpoints
    do j=1,3
      if(int_kpoints(j,i)/=0)then
        numerator(j,i)=int_kpoints(j,i)/gcd(abs(int_kpoints(j,i)),grid(j))
        denominator(j,i)=grid(j)/gcd(abs(int_kpoints(j,i)),grid(j))
      endif ! int_kpoints
    enddo ! j
    super_size(i)=lcm(denominator(1,i),denominator(2,i),denominator(3,i))
  enddo ! i
  
  found_kpoint = .false.
  label = 0
  count = 0
  
  k_t_s_file = open_write_file(kpoint_to_supercell_filename)
  ibz_file = open_write_file(args(3))
  
  do size_count=1,maxval(super_size(1:num_kpoints))
    do_i : do i=1,num_kpoints
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
                      temp_latt_vecs(k,j) = sum( dble(hnf(k,:)) &
                                             & * structure%lattice(:,j))
                    enddo ! j
                  enddo ! k
                  call minkowski_reduce(temp_latt_vecs)
                  do k=1,3
                    do j=1,3
                      hnf(k,j) = nint(sum( temp_latt_vecs(k,:) &
                                       & * structure%recip_lattice(j,:)))
                    enddo ! j
                  enddo ! k
                  
                  ! sdir="Supercell_*
                  sdir=trim(supercell_directory_root)//trim(i2s(count))
                  call system('mkdir '//trim(sdir))
                  
                  supercell_file =open_write_file(trim(sdir)//'/supercell.dat')
                  write(supercell_file,*)hnf(1,1:3)
                  write(supercell_file,*)hnf(2,1:3)
                  write(supercell_file,*)hnf(3,1:3)
                  close(supercell_file)
                  
                  size_file =open_write_file(trim(sdir)//'/size.dat')
                  write(size_file,*)super_size(i)
                  close(size_file)
                  
                endif ! tol
                if (found_kpoint(i)) cycle do_i
              enddo ! s23
            enddo ! s13
          enddo ! s12
        enddo ! s22
      enddo ! s11
    enddo do_i
  enddo ! size_count
  
  close(ibz_file)
  close(k_t_s_file)

  if(any(.not.found_kpoint(1:num_kpoints)))then
    write(*,*)'Unable to allocate each k-point to a supercell matrix.'
    stop
  endif ! found_kpoint
  
  call drop(structure)
end subroutine
end module
