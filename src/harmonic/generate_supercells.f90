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

subroutine generate_supercells(structure,grid,ibz_filename,no_sc_filename, &
   & supercell_directory_root)
  use constants,        only : dp
  use utils,            only : i2s, reduce_interval
  use linear_algebra,   only : inv_33
  use file_module,      only : open_read_file, open_write_file, count_lines
  
  use string_module
  use structure_module
  implicit none
  
  ! inputs
  type(StructureData), intent(in) :: structure
  integer,             intent(in) :: grid(:)
  type(String),        intent(in) :: ibz_filename
  type(String),        intent(in) :: no_sc_filename
  type(String),        intent(in) :: supercell_directory_root
  
  ! parameters
  real(dp),parameter :: tol=1.d-10
  
  ! Working variables
  integer,  allocatable :: multiplicity(:)
  integer,  allocatable :: int_kpoints(:,:)
  integer,  allocatable :: numerator(:,:)
  integer,  allocatable :: denominator(:,:)
  integer,  allocatable :: super_size(:)
  integer,  allocatable :: label(:)
  integer               :: ialloc
  integer               :: no_gvectors
  integer               :: sc_id
  integer               :: s11
  integer               :: s12
  integer               :: s13
  integer               :: s22
  integer               :: s23
  integer               :: s33
  integer               :: quotient
  integer               :: hnf(3,3)
  real(dp), allocatable :: kpoints(:,:)
  real(dp)              :: temp_latt_vecs(3,3)
  real(dp)              :: prim(3)
  logical,  allocatable :: found_kpoint(:)
  
  ! generate kgrid variables
  real(dp), parameter   :: tol_g = 1.d-10
  integer               :: i_symm
  integer               :: i_vec
  integer               :: j_vec
  real(dp), allocatable :: gvecs_cart(:,:)
  real(dp), allocatable :: gvecs_frac(:,:)
  real(dp)              :: temp_frac(3)
  real(dp)              :: rvec(3)
  integer,  allocatable :: rot_operation(:)
  integer               :: counter
  
  ! Supercell directory
  type(String) :: sdir
  
  ! Temporary variables
  integer :: i,j,k
  
  ! file units
  integer :: ibz_file
  integer :: no_sc_file
  integer :: supercell_file
  
  ! Generate G-vectors 
  no_gvectors = product(grid)
  
  ! Allocate arrays
  allocate( gvecs_cart(3,no_gvectors),  &
          & kpoints(3,no_gvectors),     &
          & gvecs_frac(3,no_gvectors),  &
          & rot_operation(no_gvectors), &
          & multiplicity(no_gvectors),  &
          & int_kpoints(3,no_gvectors), &
          & numerator(3,no_gvectors),   &
          & denominator(3,no_gvectors), &
          & super_size(no_gvectors),    &
          & found_kpoint(no_gvectors),  &
          & label(no_gvectors),         &
          & stat=ialloc)
  if(ialloc/=0)then
    write(*,*)'Problem allocating arrays.'
    stop
  endif
  
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
  kpoints=0.d0
  multiplicity=0
  do_i_vec : do i_vec=1,no_gvectors
    if(i_vec==1)then
      kpoints(:,i_vec) = gvecs_frac(:,i_vec)
      rot_operation(i_vec)=1
      multiplicity(i_vec)=multiplicity(i_vec)+1
    endif
    
    do j_vec=1,i_vec-1
      do i_symm=1,structure%no_symmetries
        rvec = matmul( structure%rotation_matrices(:,:,i_symm), &
                     & gvecs_cart(:,i_vec))
        temp_frac = reduce_interval(matmul(structure%lattice,rvec),tol_g)
        
        if(all(abs(temp_frac(:)-gvecs_frac(:,j_vec))<tol_g))then
          kpoints(:,i_vec) = gvecs_frac(:,j_vec)
          rot_operation(i_vec)=i_symm
          multiplicity(j_vec)=multiplicity(j_vec)+1
          cycle do_i_vec
        endif
      enddo
      
      if(j_vec==i_vec-1)then
        kpoints(:,i_vec)=gvecs_frac(:,i_vec)
        rot_operation(i_vec)=1
        multiplicity(i_vec)=multiplicity(i_vec)+1
      endif
    enddo
  enddo do_i_vec
  
  ! Generate supercells
  
  ! Express k-points as fractions
  do i=1,no_gvectors
    int_kpoints(:,i) = nint(grid*kpoints(:,i))
    do j=1,3
      if (dabs(dble(int_kpoints(j,i))/grid(j)-kpoints(j,i)) >= tol) then
        write(*,*) 'Unable to find fractional representation of k-point.'
        stop
      endif
    enddo
  enddo
  
  numerator = 0
  denominator = 1
  
  ! Reduce fractions
  do i=1,no_gvectors
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
  sc_id = 0
  
  ibz_file = open_write_file(ibz_filename)
  
  hnf = 0
  
  do_i : do i=1,no_gvectors
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
              hnf(1,:) = (/ s11, s12, s13 /)
              hnf(2,:) = (/ s22, s23, 0   /)
              hnf(3,:) = (/ s33, 0  , 0   /)
              prim = matmul(hnf,kpoints(:,i))
              if(all(abs(prim(:)-nint(prim(:)))<tol))then
                sc_id=sc_id+1
                
                ! sdir="Supercell_*
                sdir=supercell_directory_root//sc_id
                call system('mkdir '//sdir)
                
                found_kpoint(i)=.true.
                label(i)=sc_id
                
                write(ibz_file,*) kpoints(:,i),    &
                                & multiplicity(i), &
                                & sc_id
                
                do j=i+1,no_gvectors
                  if(found_kpoint(j))cycle
                  if(super_size(j)/=super_size(i))cycle
                  prim = matmul(hnf,kpoints(:,j))
                  if(all(abs(prim(:)-nint(prim(:)))<tol))then
                    found_kpoint(j)=.true.
                    label(j)=sc_id
                    
                    write(ibz_file,*) kpoints(:,j),    &
                                    & multiplicity(j), &
                                    & sc_id
                  endif
                enddo
                
                temp_latt_vecs = matmul(hnf,structure%lattice)
                
                call minkowski_reduce(temp_latt_vecs)
                
                hnf = nint(matmul( temp_latt_vecs, &
                                 & transpose(structure%recip_lattice)))
                
                supercell_file = open_write_file(sdir//'/supercell.dat')
                do j=1,3
                  write(supercell_file,*) hnf(j,:)
                enddo
                close(supercell_file)
                
              endif ! tol
              if (found_kpoint(i)) cycle do_i
            enddo ! s23
          enddo ! s13
        enddo ! s12
      enddo ! s22
    enddo ! s11
  enddo do_i
  
  close(ibz_file)
  
  no_sc_file = open_write_file(no_sc_filename)
  write(no_sc_file,*) sc_id
  close(no_sc_file)

  if(any(.not.found_kpoint(1:no_gvectors)))then
    write(*,*)'Unable to allocate each k-point to a supercell matrix.'
    stop
  endif ! found_kpoint
end subroutine
end module
