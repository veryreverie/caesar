module generate_supercells_module
  use supercell_module, only : SupercellData
  implicit none
  
  type GeneratedSupercells
    type(SupercellData), allocatable :: supercells(:)
    integer,             allocatable :: kpoints(:,:)
    integer,             allocatable :: multiplicity(:)
    integer,             allocatable :: sc_ids(:)
    integer,             allocatable :: gvector_ids(:)
  end type
  
  interface new
    module procedure new_GeneratedSupercells
  end interface
  
  interface drop
    module procedure drop_GeneratedSupercells
  end interface
  
  interface lcm
    module procedure lcm_2 ! lowest common multiple of two positive integers
    module procedure lcm_3 ! lowest common multiple of three positive integers
  end interface

contains

subroutine new_GeneratedSupercells(this,no_supercells,no_kpoints)
  implicit none
  
  type(GeneratedSupercells), intent(out) :: this
  integer,                   intent(in)  :: no_supercells
  integer,                   intent(in)  :: no_kpoints
  
  allocate(this%supercells(no_supercells))
  allocate(this%kpoints(3,no_kpoints))
  allocate(this%multiplicity(no_kpoints))
  allocate(this%sc_ids(no_kpoints))
  allocate(this%gvector_ids(no_kpoints))
end subroutine

subroutine drop_GeneratedSupercells(this)
  implicit none
  
  type(GeneratedSupercells), intent(inout) :: this
  
  deallocate(this%supercells)
  deallocate(this%kpoints)
  deallocate(this%multiplicity)
  deallocate(this%sc_ids)
  deallocate(this%gvector_ids)
end subroutine

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

! ----------------------------------------------------------------------
! Given three linearly independent input vectors a, b and c, construct the
! following linear combinations: a+b-c, a-b+c, -a+b+c, a+b+c and  check if any
! of the four new vectors is shorter than any of a, b or c. If so, replace the
! longest of a, b and c with the new (shorter) vector. The resulting three
! vectors are also linearly independent.
! ----------------------------------------------------------------------
logical function reduce_vec(vecs)
  use constants, only : dp
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

! ----------------------------------------------------------------------
! Generate all unique supercells that contain a given number of primitive
! unit cells. See 'Hart and Forcade, Phys. Rev. B 77, 224115 (2008)' for
! details of the algorithm.
! ----------------------------------------------------------------------
subroutine supercells_generator(num_pcells,num_hnf,hnf)
  use err_module
  use string_module
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
    call print_line('Problem allocating hnf array in supercells_generator.')
    call err()
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
    call print_line('Did not generate all HNF matrices.')
    call err()
  endif 
end subroutine


! ----------------------------------------------------------------------
! Given n vectors a(i) that form a basis for a lattice L in n dimensions, the
! a(i) are said to be Minkowski-reduced if the following conditions are met:
!
! - a(1) is the shortest non-zero vector in L
! - for i>1, a(i) is the shortest possible vector in L such that a(i)>=a(i-1)
!   and the set of vectors a(1) to a(i) are linearly independent
!
! In other words the a(i) are the shortest possible basis vectors for L. This
! routine, given a set of input vectors a'(i) that are possibly not
! Minkowski-reduced, returns the vectors a(i) that are.
! ----------------------------------------------------------------------
subroutine minkowski_reduce(vecs)
  use constants, only : dp
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
      if (changed) then
        cycle iter
      endif
    enddo
    
    ! Then check linear combinations involving all three.
    if (reduce_vec(vecs)) then
      cycle
    endif
    
    exit
  enddo iter
end subroutine minkowski_reduce

function generate_supercells(structure,grid) result(output)
  use constants,      only : dp
  use utils,          only : reduce_to_ibz
  use linear_algebra, only : invert_int
  use string_module
  use file_module
  use structure_module
  use supercell_module
  use err_module
  implicit none
  
  ! inputs
  type(StructureData), intent(in) :: structure
  integer,             intent(in) :: grid(:)
  type(GeneratedSupercells)       :: output
  
  ! K-point variables
  integer, allocatable :: multiplicity(:)
  integer, allocatable :: kpoints(:,:)
  integer, allocatable :: sc_size(:)
  integer, allocatable :: sc_ids(:)
  integer, allocatable :: gvector_ids(:)
  
  ! Supercell variables
  type(SupercellData), allocatable :: supercells(:)
  integer               :: hnf(3,3) ! Supercell in lattice space in normal form
  real(dp)              :: supercell_cart(3,3) ! The supercell in real space
  
  ! G-vector and K-point variables
  integer :: gvector_id
  integer :: delta_gvec(3)
  integer :: gvector(3)
  integer :: gvector_cart(3)
  integer :: kpoint_cart(3)
  
  ! Working variables
  integer,  allocatable :: denominator(:,:)
  integer               :: ialloc
  integer               :: grid_size
  integer               :: sc_num
  integer               :: s11
  integer               :: s12
  integer               :: s13
  integer               :: s22
  integer               :: s23
  integer               :: s33
  integer               :: quotient
  
  ! generate kgrid variables
  integer               :: i_symm
  integer               :: i_vec
  integer               :: j_vec
  integer,  allocatable :: kpoints_grid(:,:)
  integer,  allocatable :: rot_operation(:)
  integer               :: counter
  
  real(dp) :: rotation(3,3)
  integer  :: rvec(3)
  
  ! Temporary variables
  integer :: i,j,k,l,m
  
  ! Generate G-vectors 
  grid_size = product(grid)
  
  ! Allocate arrays
  allocate( rot_operation(grid_size), &
          & multiplicity(grid_size),  &
          & kpoints(3,grid_size),     &
          & denominator(3,grid_size), &
          & sc_size(grid_size),       &
          & sc_ids(grid_size),        &
          & gvector_ids(grid_size),   &
          & supercells(grid_size),    & ! Length of worst case.
          & stat=ialloc)
  call err(ialloc)
  
  ! ----------------------------------------------------------------------
  ! Generate all k-points in the grid.
  ! ----------------------------------------------------------------------
  allocate(kpoints_grid(3,grid_size))
  counter=0
  ! Loop over integer vectors, such that k-point/grid in [-0.5,0.5)
  ! n.b. integer floor division intentional
  
  ! This is the method used in the old code. It is in place for testing.
  do i=0,grid(1)-1
    do j=0,grid(2)-1
      do k=0,grid(3)-1
      counter = counter + 1
      kpoints_grid(:,counter) = (/i,j,k/)
      enddo
    enddo
  enddo
  
  do i=1,counter
    do j=1,3
      if (kpoints_grid(j,i) > (grid(j)-1)/2) then
        kpoints_grid(j,i) = kpoints_grid(j,i) - grid(j)
      endif
    enddo
  enddo
  
  ! This method gives the k-points in a more sensible order.
  !do i=-grid(1)/2,(grid(1)-1)/2
  !  do j=-grid(2)/2,(grid(2)-1)/2
  !    do k=-grid(3)/2,(grid(3)-1)/2
  !      counter=counter+1
  !      kpoints_grid(:,counter) = (/i,j,k/)
  !    enddo
  !  enddo
  !enddo
  
  ! ----------------------------------------------------------------------
  ! Find equivalent k-points by rotating all k-points into the IBZ.
  ! ----------------------------------------------------------------------
  multiplicity=1
  do_i_vec : do i_vec=1,grid_size
    
    ! Check if an equivalent k-point has already been found.
    do j_vec=1,i_vec-1
      do i_symm=1,structure%no_symmetries
        ! Work out rotation in fractional lattice coordinates.
        rotation = matmul(matmul( structure%lattice, &
                                & structure%rotation_matrices(:,:,i_symm)), &
                                & transpose(structure%recip_lattice))
        
        ! Rotate the k-point, and map it back to the IBZ.
        rvec = reduce_to_ibz(nint(matmul(rotation,kpoints_grid(:,i_vec))),grid)
        
        ! If the rotated k-point = k-point(j)
        if (all(rvec == kpoints_grid(:,j_vec))) then
          kpoints(:,i_vec) = kpoints_grid(:,j_vec)
          rot_operation(i_vec)=i_symm
          multiplicity(j_vec)=multiplicity(j_vec)+1
          cycle do_i_vec
        endif
      enddo
    enddo
    
    ! If kpoint not already found, assign it to the new gvec
    kpoints(:,i_vec) = kpoints_grid(:,i_vec)
    rot_operation(i_vec)=1
  enddo do_i_vec
  
  ! ----------------------------------------------------------------------
  ! Find supercells which match each k-point.
  ! ----------------------------------------------------------------------
  
  ! Calculate minimum supercell size for each k-point.
  denominator = 1
  do i=1,grid_size
    do j=1,3
      if(kpoints(j,i)/=0)then
        denominator(j,i)=grid(j)/gcd(abs(kpoints(j,i)),grid(j))
      endif
    enddo
    sc_size(i)=lcm(denominator(1,i),denominator(2,i),denominator(3,i))
  enddo
  
  sc_ids = 0
  sc_num = 0
  ! Loop over k-points in ascending order of sc_size.
  ! Find a supercell matching each k-point.
  do_i : do while (any(sc_ids==0))
    i = minloc(sc_size,1,sc_ids==0)
    do s11=1,sc_size(i)
      if(.not.mod(sc_size(i),s11)==0)cycle
      quotient=sc_size(i)/s11
      do s22=1,quotient
        if(.not.mod(quotient,s22)==0)cycle
        s33=quotient/s22
        do s12=0,s22-1
          do s13=0,s33-1
            do s23=0,s33-1
              hnf(1,:) = (/ s11, s12, s13 /)
              hnf(2,:) = (/ 0  , s22, s23 /)
              hnf(3,:) = (/ 0  , 0  , s33 /)
              if (all(modulo(matmul(hnf,kpoints(:,i)),grid)==0)) then
                ! --------------------------------------------------
                ! Add the new supercell.
                ! --------------------------------------------------
                sc_num=sc_num+1
                
                ! Reduce the supercell to minimise its real-space volume.
                supercell_cart = matmul(hnf,structure%lattice)
                call minkowski_reduce(supercell_cart)
                
                ! Allocate space for the supercell,
                !    and calculate the supercell matrix.
                call new(supercells(sc_num),sc_size(i))
                supercells(sc_num)%supercell = nint(matmul( supercell_cart, &
                   & transpose(structure%recip_lattice)))
                supercells(sc_num)%recip_supercell = invert_int(transpose( &
                   & supercells(sc_num)%supercell))
                
                ! --------------------------------------------------
                ! Assign matching k-points to the supercell.
                ! --------------------------------------------------
                sc_ids(i)=sc_num
                
                do j=i+1,grid_size
                  if(sc_ids(j) /= 0)cycle
                  if(sc_size(j)/=sc_size(i))cycle
                  if (all(modulo(matmul(hnf,kpoints(:,j)),grid)==0)) then
                    sc_ids(j)=sc_num
                  endif
                enddo
                
                cycle do_i
              endif
            enddo ! s23
          enddo ! s13
        enddo ! s12
      enddo ! s22
    enddo ! s11
  enddo do_i
  
  ! Check that all kpoints have been assigned to supercells.
  if(any(sc_ids==0))then
    call print_line('Unable to allocate each k-point to a supercell matrix.')
    call err()
  endif
  
  ! ----------------------------------------------------------------------
  ! Calculate the supercell G-vectors in the IBZ of the primitive cell.
  ! ----------------------------------------------------------------------
  ! N.B. recip_supercell transforms G-vectors into a representation where
  !   the IBZ of the prim. cell is a cartesian cube with side length sc_size.
  ! In this space, supercell BZs are parallelapipeds of volume sc_size*sc_size,
  !    with integer recip. latt. vectors.
  do i=1,sc_num
    ! Loop over all possible G-vectors in the primitive cell.
    gvector_id = 1
    do_j : do j=0,supercells(i)%sc_size-1
      do k=0,supercells(i)%sc_size-1
        do_l : do l=0,supercells(i)%sc_size-1
          gvector = (/l,k,j/)
          
          ! Check if gvector has already been found.
          do m=1,gvector_id-1
            delta_gvec = matmul( transpose(supercells(i)%recip_supercell), &
                               & gvector - supercells(i)%gvectors(:,m))
            if (all(modulo(delta_gvec,supercells(i)%sc_size)==0)) then
              cycle do_l
            endif
          enddo
          
          supercells(i)%gvectors(:,gvector_id) = gvector
          gvector_id = gvector_id + 1
          
          if (gvector_id > supercells(i)%sc_size) then
            exit do_l
          endif
        enddo do_l
      enddo
    enddo do_j
    
    ! Map G-vectors to [-0.5,0.5).
    do j=1,supercells(i)%sc_size
      gvector_cart = matmul( transpose(supercells(i)%recip_supercell), &
                           & supercells(i)%gvectors(:,j))
      gvector_cart = modulo(gvector_cart,supercells(i)%sc_size)
      do k=1,3
        if (gvector_cart(k) > (supercells(i)%sc_size-1)/2) then
          gvector_cart(k) = gvector_cart(k) - supercells(i)%sc_size
        endif
      enddo
      supercells(i)%gvectors(:,j) = matmul( supercells(i)%supercell, &
                                &           gvector_cart)            &
                                & / supercells(i)%sc_size
    enddo
    
    ! Check that the correct number of gvectors have been found.
    if (gvector_id<=supercells(i)%sc_size) then
      call print_line('Error: Wrong number of G-vectors found.')
      call print_line('Supercell size: '//supercells(i)%sc_size)
      call print_line('No. G-vectors found: '//gvector_id-1)
      call err()
    endif
  enddo
  
  ! ----------------------------------------------------------------------
  ! Match K-points and G-vectors.
  ! ----------------------------------------------------------------------
  ! To convert to real-space cartesian representation:
  !   kpoint_cart = kpoint/grid
  !   gvec_cart   = matmul(transpose(recip_supercell),gvec)/sc_size
  ! These are paired when kpoint_cart = gvec_cart, down to unit translations.
  ! To keep calculation in integer representation,
  !   grid and sc_size are multiplied through.
  gvector_ids = 0
  do i=1,grid_size
    ! Calculate k-point in scaled real space.
    kpoint_cart = kpoints(:,i)*sc_size(i)
    do j=1,sc_size(i)
      ! Calculate g-vector in scaled real space.
      gvector_cart = matmul( transpose(supercells(sc_ids(i))%recip_supercell),&
                 &           supercells(sc_ids(i))%gvectors(:,j))             &
                 & * grid
      ! Check if the k-point and g-vector are equivalent down to translations.
      if (all(modulo(kpoint_cart-gvector_cart,sc_size(i)*grid)==0)) then
        gvector_ids(i) = j
        exit
      endif
    enddo
    
    ! Check that the corresponding gvector has been found.
    if (gvector_ids(i) == 0) then
      call print_line("Error: could not locate G-vector.")
      call err()
    endif
  enddo
  
  ! ----------------------------------------------------------------------
  ! Generate output.
  ! ----------------------------------------------------------------------
  call new(output,sc_num,grid_size)
  output%supercells = supercells(1:sc_num)
  output%kpoints = kpoints
  output%multiplicity = multiplicity
  output%sc_ids = sc_ids
  output%gvector_ids = gvector_ids
end function
end module
