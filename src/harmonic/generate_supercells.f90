! ======================================================================
! Generates the supercells needed to simulate phonons at all q-points.
! ======================================================================
module generate_supercells_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use structure_module
  use qpoints_module
  implicit none
  
  type GeneratedSupercells
    type(QpointData),    allocatable :: qpoints_ibz(:)
    type(StructureData), allocatable :: supercells(:)
  end type
  
  interface GeneratedSupercells
    module procedure new_GeneratedSupercells
  end interface
  
  interface lcm
    module procedure lcm_2 ! lowest common multiple of two positive integers
    module procedure lcm_3 ! lowest common multiple of three positive integers
  end interface

contains

! ----------------------------------------------------------------------
! GeneratedSupercells allocator.
! ----------------------------------------------------------------------
function new_GeneratedSupercells(no_qpoints_ibz,no_supercells) result(this)
  implicit none
  
  integer, intent(in)       :: no_qpoints_ibz
  integer, intent(in)       :: no_supercells
  type(GeneratedSupercells) :: this
  
  integer :: ialloc
  
  allocate( this%qpoints_ibz(no_qpoints_ibz), &
          & this%supercells(no_supercells),   &
          & stat=ialloc); call err(ialloc)
end function

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
function minkowski_reduce(input) result(output)
  use linear_algebra_module
  implicit none
  
  type(RealMatrix) :: input
  type(RealMatrix) :: output
  
  integer  :: i
  real(dp) :: tempmat(3,3)
  real(dp) :: tempvec(3)
  logical  :: changed
  
  tempmat = dble(input)
  
  iter: do
    ! First check linear combinations involving two vectors.
    do i=1,3
      tempvec = tempmat(i,:)
      tempmat(i,:) = 0
      changed = reduce_vec(tempmat)
      tempmat(i,:) = tempvec
      
      if (changed) then
        cycle iter
      endif
    enddo
    
    ! Then check linear combinations involving all three.
    if (reduce_vec(tempmat)) then
      cycle
    endif
    
    exit
  enddo iter
  
  output = tempmat
end function

! ----------------------------------------------------------------------
! Find a supercell matrix, S, s.t. S.q is a vector of integers.
! ----------------------------------------------------------------------
! Returns answer in Hermite Normal Form.
function find_supercell_matrix(qpoint,sc_size,scaling) result(output)
  use linear_algebra_module
  implicit none
  
  type(IntVector), intent(in) :: qpoint
  integer,         intent(in) :: sc_size
  integer,         intent(in) :: scaling
  type(IntMatrix)             :: output  ! S
  
  integer :: s11,s12,s13
  integer ::     s22,s23
  integer ::         s33
  
  ! Loop over all matrices in Hermite Normal Form.
  ! s11*s22*s33 = det(S) = sc_size.
  ! 0 <= s12 < s22
  ! 0 <= s13 < s33
  ! 0 <= s23 < s33
  do s11=1,sc_size
    if (modulo(sc_size,s11)==0) then
      do s22=1,sc_size/s11
        if (modulo(sc_size,s11*s22)==0) then
          s33=sc_size/(s11*s22)
          do s12=0,s22-1
            do s13=0,s33-1
              do s23=0,s33-1
                ! Construct S.
                output = mat([ s11, s12, s13, &
                             & 0  , s22, s23, &
                             & 0  , 0  , s33  ], 3,3)
                
                ! Check if S.q is a vector of integers (scaled by scaling).
                if (all(modulo( int(output*qpoint), scaling) == 0)) then
                  return
                endif
              enddo
            enddo
          enddo
        endif
      enddo
    endif
  enddo
  
  ! Throw an error if no supercell could be found.
  call print_line('Error: no supercell matrix found for q-point '//qpoint)
  call err()
end function

! ----------------------------------------------------------------------
! Constructs a supercell which has the given q-point as a G-vector.
! ----------------------------------------------------------------------
function find_supercell(qpoint,sc_size,structure,scaling,symmetry_precision) &
   & result(output)
  use linear_algebra_module
  use construct_supercell_module
  implicit none
  
  type(IntVector),     intent(in) :: qpoint
  integer,             intent(in) :: sc_size
  type(StructureData), intent(in) :: structure
  integer,             intent(in) :: scaling
  real(dp),            intent(in) :: symmetry_precision
  type(StructureData)             :: output
  
  type(IntMatrix)  :: supercell_frac
  type(RealMatrix) :: supercell_cart
  
  ! Find a supercell matrix S, s.t. S.q is a vector of integers.
  supercell_frac = find_supercell_matrix(qpoint, sc_size, scaling)
  
  ! Minkowski reduce the supercell in cartesian co-ordinates.
  ! Converts to cartesian, reduces, and converts back again.
  supercell_cart = supercell_frac*structure%lattice
  supercell_cart = minkowski_reduce(supercell_cart)
  supercell_frac = nint(dble( supercell_cart &
                          & * transpose(structure%recip_lattice)))
  
  ! Generate the supercell from its supercell matrix.
  output = construct_supercell( structure,      &
                              & supercell_frac, &
                              & symmetry_precision=symmetry_precision)
end function

! ----------------------------------------------------------------------
! Checks supercells and q-points.
! ----------------------------------------------------------------------
! Throws an error if there is a problem.
subroutine check_supercells(input,structure,large_supercell)
  implicit none
  
  type(GeneratedSupercells), intent(in) :: input
  type(StructureData),       intent(in) :: structure
  type(StructureData),       intent(in) :: large_supercell
  
  ! Working variables.
  type(IntVector)               :: gvector_large_supercell
  type(QpointData)              :: qpoint
  type(RealVector)              :: gvector_qpoint_real
  type(IntVector)               :: gvector_qpoint
  type(IntMatrix)               :: rotation
  type(StructureData)           :: supercell
  type(IntVector)               :: gvector_supercell
  real(dp)                      :: fractional_position(3)
  integer                       :: atom_1
  integer                       :: atom_2
  type(RealVector)              :: fractional_difference
  type(RealMatrix), allocatable :: rotations(:)
  type(RealMatrix)              :: rotation_identity
  
  ! Temporary variables.
  integer :: i,j,k,l
  
  ! --------------------------------------------------
  ! Check all G-vectors in the M-P grid have been matched with IBZ q-points.
  ! --------------------------------------------------
  do_i : do i=1,large_supercell%sc_size
    gvector_large_supercell = large_supercell%gvectors(i)
    
    ! Find the q-point which this G-vector maps onto.
    do j=1,size(input%qpoints_ibz)
      qpoint = input%qpoints_ibz(j)
      
      do k=1,size(qpoint%gvectors)
        if (qpoint%gvectors(k) == i) then
          
          ! Check that the q-point corresponds to a G-vector.
          gvector_qpoint_real = large_supercell%supercell * qpoint%qpoint
          gvector_qpoint = nint(dble(gvector_qpoint_real))
          if (l2_norm(gvector_qpoint_real-gvector_qpoint) > 1.0e-10_dp) then
            call print_line('Code Error: q-point '//j//' does not lie on the &
               &specified Monkhorst-Pack grid.')
            call err()
          endif
          
          ! Check that the q-point corresponds to the correct G-vector.
          rotation = structure%symmetries(qpoint%symmetry_ids(k))%rotation
          if (rotation*gvector_large_supercell /= gvector_qpoint) then
            call print_line('Code Error: G-vector '//i//' has been &
               &incorrectly matched to a q-point and supercell.')
            call err()
          endif
          
          ! Check the supercell to which the q-point has been assigned.
          supercell = input%supercells(qpoint%sc_id)
          gvector_supercell = supercell%gvectors(qpoint%gvector_id)
          if ( transpose(large_supercell%recip_supercell) &
           & * gvector_qpoint                             &
           & * supercell%sc_size                          &
           & /=                                           &
           &   transpose(supercell%recip_supercell)       &
           & * gvector_supercell                          &
           & * large_supercell%sc_size) then
            call print_line('Code Error: The G-vectors corresponding to &
               &q-point '//j//' do not match.')
            call err()
          endif
          
          cycle do_i
        endif
      enddo
    enddo
    
    call print_line('Code Error: G-vector '//i//' has not been matched with &
       &a suitable q-point and supercell.')
    call err()
  enddo do_i
  
  ! --------------------------------------------------
  ! Check that all supercells have been constructed correctly.
  ! --------------------------------------------------
  do i=1,size(input%supercells)
    supercell = input%supercells(i)
    
    do j=1,supercell%no_atoms
      ! Check that no atoms are on top of one another.
      do k=1,j-1
        if (l2_norm( supercell%atoms(j)%cartesian_position() &
                 & - supercell%atoms(k)%cartesian_position() )<0.5_dp) then
          call print_line('Code Error: atoms '//k//' and '//j//' are within &
             &0.5 Bohr of one another in supercell '//i//'.')
          call err()
        endif
      enddo
      
      ! Check that all atoms are within the supercell's primitive cell.
      fractional_position = dble(supercell%atoms(j)%fractional_position())
      if ( any(fractional_position <  -1.0e-10_dp) .or. &
         & any(fractional_position > 1+1.0e-10_dp)) then
        call print_line('Code Error: atom '//j//' has been placed outside &
           &the primitive cell of supercell '//i//'.')
        call err()
      endif
    enddo
    
    ! Check that all atoms related by G-vectors are correctly placed.
    do j=1,structure%no_atoms
      do k=1,supercell%sc_size
        do l=1,k-1
          atom_1 = supercell%rvec_and_prim_to_atom(j,k)
          atom_2 = supercell%rvec_and_prim_to_atom(j,l)
          fractional_difference =                              &
             &   supercell%atoms(atom_1)%fractional_position() &
             & - supercell%atoms(atom_2)%fractional_position()
          if ( l2_norm( fractional_difference &
           &          - vec(nint(dble(fractional_difference)))) &
           & > 1.0e-10_dp) then
            call print_line('Code Error: atoms '//atom_1//' and '//atom_2// &
               & ' in supercell '//i//' are not related by a G-vector')
            call err()
          endif
        enddo
      enddo
    enddo
    
    ! Check supercell rotations
    rotations = supercell%calculate_cartesian_rotations()
    
    do j=1,size(supercell%symmetries)
      ! Check that R.R^T is the identity.
      rotation_identity = rotations(j) * transpose(rotations(j))
      if (any(abs(dble( rotation_identity &
                    & - make_identity_matrix(3)))>1.0e-10_dp)) then
        call print_line('Code Error: rotation '//j//' in supercell '//i// &
           & ' is not a rotation or an improper rotation.')
        call err()
      endif
    enddo
  enddo
end subroutine

! ----------------------------------------------------------------------
! The main program.
! ----------------------------------------------------------------------
function generate_supercells(structure,large_supercell,symmetry_precision) &
   & result(output)
  use linear_algebra_module
  use structure_module
  use construct_supercell_module
  use group_module
  implicit none
  
  ! Inputs.
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: large_supercell
  real(dp),            intent(in) :: symmetry_precision
  type(GeneratedSupercells)       :: output
  
  ! Grid q-point variables.
  integer, allocatable :: grid_to_ibz(:)
  integer, allocatable :: ibz_to_grid(:)
  integer, allocatable :: rotation_to_ibz(:)
  
  ! IBZ q-point variables.
  integer                      :: no_qpoints_ibz
  integer,         allocatable :: multiplicity(:)
  type(IntVector), allocatable :: qpoints_grid(:)
  type(IntVector), allocatable :: qpoints_ibz(:)
  integer,         allocatable :: sc_size(:)
  integer,         allocatable :: sc_ids(:)
  integer,         allocatable :: gvector_ids(:)
  type(IntVector)              :: qpoint
  
  ! Supercell variables
  type(StructureData), allocatable :: supercells(:)
  
  ! Working variables
  integer :: denominator(3)
  integer :: sc_num
  
  type(IntVector) :: rot_qpoint
  
  ! Temporary variables
  integer :: i,j,k,ialloc
  integer :: symmetry
  
  ! --------------------------------------------------
  ! Construct q-points from G-vectors of large supercell.
  ! --------------------------------------------------
  ! N.B. all q-points stored in scaled fractional reciprocal lattice 
  !    co-ordinates, such that the reciprocal primitive lattice vectors are
  !    cartesian vectors of length det(large_supercell%supercell).
  ! This is to keep everything in an integer representation.
  allocate( qpoints_grid(large_supercell%sc_size), &
          & stat=ialloc); call err(ialloc)
  do i=1,large_supercell%sc_size
    qpoints_grid(i) = transpose(large_supercell%recip_supercell) &
                  & * large_supercell%gvectors(i)
  enddo
  
  ! --------------------------------------------------
  ! Find equivalent q-points by rotating all q-points into the IBZ.
  ! --------------------------------------------------
  allocate( grid_to_ibz(large_supercell%sc_size),     &
          & ibz_to_grid(large_supercell%sc_size),     &
          & rotation_to_ibz(large_supercell%sc_size), &
          & stat=ialloc); call err(ialloc)
  no_qpoints_ibz = 0
  do_i1 : do i=1,large_supercell%sc_size
    ! Check if an equivalent-by-symmetry q-point has already been found.
    do j=1,no_qpoints_ibz
      do k=1,size(structure%symmetries)
        
        ! Rotate the q-point.
        rot_qpoint = structure%symmetries(k)%rotation * qpoints_grid(i)
        
        ! If the rotated q-point = q-point(j).
        if (rot_qpoint-qpoints_grid(ibz_to_grid(j))==vec([0,0,0])) then
          grid_to_ibz(i) = j
          rotation_to_ibz(i) = k
          cycle do_i1
        endif
      enddo
    enddo
    
    ! If qpoint not already found, assign it to the new G-vector.
    no_qpoints_ibz = no_qpoints_ibz+1
    grid_to_ibz(i) = no_qpoints_ibz
    ibz_to_grid(no_qpoints_ibz) = i
    rotation_to_ibz(i) = 1
  enddo do_i1
  
  allocate( qpoints_ibz(no_qpoints_ibz),  &
          & multiplicity(no_qpoints_ibz), &
          & stat=ialloc); call err(ialloc)
  multiplicity=0
  do i=1,large_supercell%sc_size
    if (multiplicity(grid_to_ibz(i))==0) then
      qpoints_ibz(grid_to_ibz(i)) = qpoints_grid(i)
    endif
    multiplicity(grid_to_ibz(i)) = multiplicity(grid_to_ibz(i)) + 1
  enddo
  
  ! ----------------------------------------------------------------------
  ! For each q-point in the IBZ, find the smallest possible supercell
  !    such that that q-point is a G-vector of the supercell.
  ! ----------------------------------------------------------------------
  
  ! Calculate minimum supercell size for each q-point.
  allocate(sc_size(no_qpoints_ibz), stat=ialloc); call err(ialloc)
  do i=1,no_qpoints_ibz
    denominator = large_supercell%sc_size &
              & / gcd( abs(int(qpoints_ibz(i))), large_supercell%sc_size)
    sc_size(i)=lcm(denominator(1),denominator(2),denominator(3))
  enddo
  
  ! Loop over q-points in ascending order of sc_size.
  ! Find a supercell matching each q-point.
  allocate( sc_ids(no_qpoints_ibz),     &
          & supercells(no_qpoints_ibz), & ! Length of worst case.
          & stat=ialloc); call err(ialloc)
  sc_ids = 0
  sc_num = 0
  do while (any(sc_ids==0))
    i = minloc(sc_size,1,sc_ids==0)
    
    sc_num = sc_num+1
    sc_ids(i) = sc_num
    
    ! Construct a supercell for which this q-point is a G-vector.
    supercells(sc_num) = find_supercell( qpoints_ibz(i),          &
                                       & sc_size(i),              &
                                       & structure,               &
                                       & large_supercell%sc_size, &
                                       & symmetry_precision)
    
    ! Find all other q-points which are G-vectors of the supercell.
    do j=1,large_supercell%sc_size
      ! Ignore q-points which cannot match.
      if (sc_ids(grid_to_ibz(j)) /= 0) then
        cycle
      elseif (sc_size(grid_to_ibz(j)) /= sc_size(i)) then
        cycle
      
      ! Check for matching q-points.
      elseif (all(modulo( int( supercells(sc_num)%supercell &
                        &    * qpoints_grid(j)),            &
                        & large_supercell%sc_size) == 0)) then
        ! Assign the q-point to the supercell.
        sc_ids(grid_to_ibz(j)) = sc_num
        
        ! Re-arrange ibz to grid mapping so that simulated q-point
        !    is in ibz.
        qpoints_ibz(grid_to_ibz(j)) = qpoints_grid(j)
        ibz_to_grid(ibz_to_grid(j)) = j
        do k=1,large_supercell%sc_size
          if (grid_to_ibz(k) == grid_to_ibz(j)) then
            symmetry = structure%symmetries(rotation_to_ibz(j))%inverse
            rotation_to_ibz(k) =                                 &
               &   structure%symmetries(symmetry)%operator_group &
               & * rotation_to_ibz(k)
          endif
        enddo
        
      endif
    enddo
  enddo
  
  ! Check that all q-points have been assigned to supercells.
  if(any(sc_ids==0))then
    call print_line('Unable to allocate every q-point to a supercell matrix.')
    call err()
  endif
  
  ! --------------------------------------------------
  ! Match q-points and G-vectors.
  ! --------------------------------------------------
  allocate(gvector_ids(no_qpoints_ibz), stat=ialloc); call err(ialloc)
  gvector_ids = 0
  do i=1,no_qpoints_ibz
    do j=1,sc_size(i)
      ! Calculate q-point corresponding to G-vector in scaled co-ordinates.
      qpoint = transpose(supercells(sc_ids(i))%recip_supercell) &
           & * supercells(sc_ids(i))%gvectors(j)
      
      ! Check if the q-points match, modulo reciprocal lattice vectors,
      !    taking account of the different co-ordinate scaling.
      if (all(modulo( int(  qpoints_ibz(i) * sc_size(i)               &
                    &     - qpoint         * large_supercell%sc_size), &
                    & sc_size(i)*large_supercell%sc_size) == 0)) then
        gvector_ids(i) = j
        exit
      endif
    enddo
    
    ! Check that the corresponding G-vector has been found.
    if (gvector_ids(i) == 0) then
      call print_line('Error: could not locate q-point '//i//':')
      call print_line(qpoints_ibz(i)/real(large_supercell%sc_size,dp))
      call print_line('G-vectors :')
      do j=1,sc_size(i)
        call print_line(supercells(sc_ids(i))%gvectors(j))
      enddo
      call err()
    endif
  enddo
  
  ! --------------------------------------------------
  ! Generate output.
  ! --------------------------------------------------
  output = GeneratedSupercells(no_qpoints_ibz,sc_num)
  
  do i=1,no_qpoints_ibz
    output%qpoints_ibz(i) = QpointData(multiplicity(i))
    output%qpoints_ibz(i)%qpoint = qpoints_ibz(i) &
                               & / dble(large_supercell%sc_size)
    output%qpoints_ibz(i)%sc_id = sc_ids(i)
    output%qpoints_ibz(i)%gvector_id = gvector_ids(i)
    
    k=1
    do j=1,large_supercell%sc_size
      if (grid_to_ibz(j)==i) then
        output%qpoints_ibz(i)%gvectors(k) = j
        output%qpoints_ibz(i)%symmetry_ids(k) = rotation_to_ibz(j)
        k = k+1
      endif
    enddo
  enddo
  
  do i=1,sc_num
    output%supercells(i) = supercells(i)
  enddo
  
  ! --------------------------------------------------
  ! Check output.
  ! --------------------------------------------------
  call check_supercells(output,structure,large_supercell)
end function
end module
