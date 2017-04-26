module generate_supercells_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use supercell_module
  use structure_module
  use qpoints_module
  implicit none
  
  type GeneratedSupercells
    type(StructureData)              :: structure_grid
    type(QpointData),    allocatable :: qpoints_ibz(:)
    type(SupercellData), allocatable :: supercells(:)
  end type
  
  interface new
    module procedure new_GeneratedSupercells
  end interface
  
  interface lcm
    module procedure lcm_2 ! lowest common multiple of two positive integers
    module procedure lcm_3 ! lowest common multiple of three positive integers
  end interface

contains

! ----------------------------------------------------------------------
! Allocates all arrays.
! ----------------------------------------------------------------------
subroutine new_GeneratedSupercells(this,no_qpoints_ibz,no_supercells)
  implicit none
  
  type(GeneratedSupercells), intent(out) :: this
  integer,                   intent(in)  :: no_qpoints_ibz
  integer,                   intent(in)  :: no_supercells
  
  integer :: ialloc
  
  allocate( this%qpoints_ibz(no_qpoints_ibz), &
          & this%supercells(no_supercells),   &
          & stat=ialloc); call err(ialloc)
end subroutine

! ----------------------------------------------------------------------
! Calculate the greatest common divisor of two positive integers using
! Euclid's algorithm.
! ----------------------------------------------------------------------
elemental function gcd(int_1,int_2) result(output)
  implicit none
  
  integer, intent(in) :: int_1
  integer, intent(in) :: int_2
  integer             :: output
  
  integer :: a
  integer :: b
  integer :: temp
  
  a = max(int_1,int_2)
  b = min(int_1,int_2)
  
  if (b==0) then
    output=a
  
  else
    do
      temp=modulo(a,b)
      if(temp==0)exit
      a=b
      b=temp
    enddo
    output=b
  endif
end function

! ----------------------------------------------------------------------
! Calculate the lowest common multiple of two positive integers
! lcm(a,b) = a*b/gcd(a,b)
! ----------------------------------------------------------------------
elemental function lcm_2(int_1,int_2) result(output)
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
elemental function lcm_3(int_1,int_2,int_3) result(output)
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
  use utils_module, only : l2_norm
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
    nlen = l2_norm(vecs(i,:))
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
    nlen = l2_norm(newvecs(i,:))
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
  implicit none
  
  integer,intent(in) :: num_pcells
  integer,intent(out) :: num_hnf
  
  integer,pointer :: hnf(:,:,:)
  integer :: a,b,c,d,e,f,ialloc,count_hnf,quotient
  
  count_hnf=0
  
  do a=1,num_pcells 
    if(.not.modulo(num_pcells,a)==0)cycle
    quotient=num_pcells/a
    do c=1,quotient  
      if(.not.modulo(quotient,c)==0)cycle
      f=quotient/c
      count_hnf=count_hnf+c*f**2
    enddo ! c
  enddo ! a
  
  num_hnf=count_hnf
  count_hnf=0
  
  allocate(hnf(3,3,num_hnf), stat=ialloc); call err(ialloc)
  hnf = 0
  do a=1,num_pcells 
    if(.not.modulo(num_pcells,a)==0)cycle
    quotient=num_pcells/a
    do c=1,quotient  
      if(.not.modulo(quotient,c)==0)cycle
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
  use linear_algebra_module, only : invert_int
  use structure_module
  use supercell_module
  use construct_supercell_module
  implicit none
  
  ! Inputs.
  type(StructureData), intent(in) :: structure
  integer,             intent(in) :: grid(:)
  type(GeneratedSupercells)       :: output
  
  ! Grid supercell structure.
  type(SupercellData) :: supercell_grid
  type(StructureData) :: structure_grid
  
  ! Grid q-point variables.
  integer, allocatable :: grid_to_ibz(:)
  integer, allocatable :: ibz_to_grid(:)
  integer, allocatable :: rotation_ids(:)
  integer              :: counter
  
  ! IBZ q-point variables.
  integer              :: no_qpoints_ibz
  integer, allocatable :: multiplicity(:)
  integer, allocatable :: qpoints_grid(:,:)
  integer, allocatable :: qpoints_ibz(:,:)
  integer, allocatable :: sc_size(:)
  integer, allocatable :: sc_ids(:)
  integer, allocatable :: gvector_ids(:)
  integer              :: gvector_id
  integer              :: delta_gvec(3)
  integer              :: gvector(3)
  integer              :: qpoint(3)
  
  ! Supercell variables
  type(SupercellData), allocatable :: supercells(:)
  integer                          :: hnf(3,3) ! Frac. co-ords, normal form.
  real(dp)                         :: supercell_cart(3,3) ! Real space co-ords.
  
  ! Working variables
  integer :: denominator(3)
  integer :: ialloc
  integer :: sc_num
  integer :: s11
  integer :: s12
  integer :: s13
  integer :: s22
  integer :: s23
  integer :: s33
  integer :: quotient
  
  integer :: rot_qpoint(3)
  
  ! Temporary variables
  integer :: i,j,k,l,m
  
  ! ----------------------------------------------------------------------
  ! Construct a supercell corresponding to the entire q-point grid.
  ! ----------------------------------------------------------------------
  call new(supercell_grid, product(grid))
  
  supercell_grid%supercell(1,:) = [ grid(1), 0      , 0       ]
  supercell_grid%supercell(2,:) = [ 0      , grid(2), 0       ]
  supercell_grid%supercell(3,:) = [ 0      , 0      , grid(3) ]
  
  counter = 0
  do i=-grid(1)/2,(grid(1)-1)/2
    do j=-grid(2)/2,(grid(2)-1)/2
      do k=-grid(3)/2,(grid(3)-1)/2
        counter=counter+1
        supercell_grid%gvectors(:,counter) = [ i, j, k ]
      enddo
    enddo
  enddo
  
  structure_grid = construct_supercell(structure, supercell_grid, .false.)
  
  ! ----------------------------------------------------------------------
  ! Transform G-vectors to equivalent q-points.
  ! ----------------------------------------------------------------------
  ! N.B. all q-points stored in scaled fractional reciprocal lattice 
  !    co-ordinates, such that the reciprocal primitive lattice vectors are
  !    cartesian vectors of length product(grid).
  allocate( qpoints_grid(3,structure_grid%sc_size), &
          & stat=ialloc); call err(ialloc)
  qpoints_grid = matmul( structure_grid%recip_supercell, &
                       & structure_grid%gvectors)
  
  ! ----------------------------------------------------------------------
  ! Find equivalent q-points by rotating all q-points into the IBZ.
  ! ----------------------------------------------------------------------
  allocate( grid_to_ibz(structure_grid%sc_size),  &
          & ibz_to_grid(structure_grid%sc_size),  &
          & rotation_ids(structure_grid%sc_size), &
          & stat=ialloc); call err(ialloc)
  no_qpoints_ibz = 0
  do_i1 : do i=1,structure_grid%sc_size
    ! Check if an equivalent-by-symmetry q-point has already been found.
    do j=1,no_qpoints_ibz
      do k=1,structure%no_symmetries
        
        ! Rotate the q-point.
        rot_qpoint = matmul(structure%rotations(:,:,k), qpoints_grid(:,i))
        
        ! If the rotated q-point = q-point(j), modulo recip. latt. vecs.
        if (all(modulo( rot_qpoint - qpoints_grid(:,ibz_to_grid(j)), &
                      & structure_grid%sc_size) == 0)) then
          grid_to_ibz(i) = j
          rotation_ids(i) = k
          cycle do_i1
        endif
      enddo
    enddo
    
    ! If qpoint not already found, assign it to the new gvec
    no_qpoints_ibz = no_qpoints_ibz+1
    grid_to_ibz(i) = no_qpoints_ibz
    ibz_to_grid(no_qpoints_ibz) = i
    rotation_ids(i) = 1
  enddo do_i1
  
  allocate( qpoints_ibz(3,no_qpoints_ibz), &
          & multiplicity(no_qpoints_ibz),  &
          & stat=ialloc); call err(ialloc)
  multiplicity=0
  do i=1,structure_grid%sc_size
    if (multiplicity(grid_to_ibz(i))==0) then
      qpoints_ibz(:,grid_to_ibz(i)) = qpoints_grid(:,i)
    endif
    multiplicity(grid_to_ibz(i)) = multiplicity(grid_to_ibz(i)) + 1
  enddo
  
  ! ----------------------------------------------------------------------
  ! Find supercells which match each q-point in the IBZ.
  ! ----------------------------------------------------------------------
  ! Calculate minimum supercell size for each q-point.
  allocate(sc_size(no_qpoints_ibz), stat=ialloc); call err(ialloc)
  do i=1,no_qpoints_ibz
    denominator = structure_grid%sc_size &
              & / gcd( abs(qpoints_ibz(:,i)), structure_grid%sc_size)
    sc_size(i)=lcm(denominator(1),denominator(2),denominator(3))
  enddo
  
  ! Loop over q-points in ascending order of sc_size.
  ! Find a supercell matching each q-point.
  allocate( sc_ids(no_qpoints_ibz),     &
          & supercells(no_qpoints_ibz), & ! Length of worst case.
          & stat=ialloc); call err(ialloc)
  sc_ids = 0
  sc_num = 0
  do_i2 : do while (any(sc_ids==0))
    i = minloc(sc_size,1,sc_ids==0)
    do s11=1,sc_size(i)
      if(.not.modulo(sc_size(i),s11)==0)cycle
      quotient=sc_size(i)/s11
      do s22=1,quotient
        if(.not.modulo(quotient,s22)==0)cycle
        s33=quotient/s22
        do s12=0,s22-1
          do s13=0,s33-1
            do s23=0,s33-1
              hnf(1,:) = [ s11, s12, s13 ]
              hnf(2,:) = [ 0  , s22, s23 ]
              hnf(3,:) = [ 0  , 0  , s33 ]
              if (all(modulo( matmul(hnf, qpoints_ibz(:,i)), &
                            & structure_grid%sc_size)==0)) then
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
                ! Assign matching q-points to the supercell.
                ! --------------------------------------------------
                sc_ids(i)=sc_num
                
                do j=i+1,no_qpoints_ibz
                  if(sc_ids(j) /= 0)cycle
                  if(sc_size(j)/=sc_size(i))cycle
                  if (all(modulo( matmul(hnf, qpoints_ibz(:,j)), &
                                & structure_grid%sc_size)==0)) then
                    sc_ids(j)=sc_num
                  endif
                enddo
                
                cycle do_i2
              endif
            enddo ! s23
          enddo ! s13
        enddo ! s12
      enddo ! s22
    enddo ! s11
  enddo do_i2
  
  ! Check that all qpoints have been assigned to supercells.
  if(any(sc_ids==0))then
    call print_line('Unable to allocate each q-point to a supercell matrix.')
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
          gvector = [ l, k, j ]
          
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
    ! Translate G-vectors by lattice vectors such that the corresponding
    !    q-points are in [-0.5,0.5).
    do j=1,supercells(i)%sc_size
      qpoint = matmul( transpose(supercells(i)%recip_supercell), &
                           & supercells(i)%gvectors(:,j))
      qpoint = modulo(qpoint,supercells(i)%sc_size)
      do k=1,3
        if (qpoint(k) > (supercells(i)%sc_size-1)/2) then
          qpoint(k) = qpoint(k) - supercells(i)%sc_size
        endif
      enddo
      supercells(i)%gvectors(:,j) = matmul(supercells(i)%supercell, qpoint) &
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
  ! Match q-points and G-vectors.
  ! ----------------------------------------------------------------------
  ! To convert to fractional primitive cell co-ordinates:
  !   qpoint_prim = qpoint/grid
  !   gvec_cart   = matmul(transpose(recip_supercell),gvec)/sc_size
  ! These are paired when qpoint_prim = gvec_cart, down to unit translations.
  ! To keep calculation in integer representation,
  !   both sc_sizes are multiplied through.
  allocate(gvector_ids(no_qpoints_ibz), stat=ialloc); call err(ialloc)
  gvector_ids = 0
  do i=1,no_qpoints_ibz
    do j=1,sc_size(i)
      ! Calculate q-point corresponding to g-vector.
      qpoint = matmul( transpose(supercells(sc_ids(i))%recip_supercell), &
                     & supercells(sc_ids(i))%gvectors(:,j))
      
      ! Check if the q-points match, modulo reciprocal lattice vectors.
      if (all(modulo(   qpoints_ibz(:,i) * sc_size(i) &
                    & - qpoint           * structure_grid%sc_size, &
                    & sc_size(i)*structure_grid%sc_size)==0)) then
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
  call new(output,no_qpoints_ibz,sc_num)
  
  output%structure_grid = structure_grid
  
  do i=1,no_qpoints_ibz
    call new(output%qpoints_ibz(i),multiplicity(i))
    output%qpoints_ibz(i)%qpoint = qpoints_ibz(:,i) &
                               & / dble(structure_grid%sc_size)
    output%qpoints_ibz(i)%sc_id = sc_ids(i)
    output%qpoints_ibz(i)%gvector_id = gvector_ids(i)
    
    k=1
    do j=1,structure_grid%sc_size
      if (grid_to_ibz(j)==i) then
        output%qpoints_ibz(i)%gvectors(k) = j
        output%qpoints_ibz(i)%rotations(k) = rotation_ids(j)
        k = k+1
      endif
    enddo
  enddo
  
  output%supercells = supercells(1:sc_num)
end function
end module
