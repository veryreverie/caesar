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
  
  private
  
  public :: generate_qpoints
  public :: generate_supercells
contains
! ----------------------------------------------------------------------
! Given three linearly independent input vectors a, b and c, construct the
! following linear combinations: a+b-c, a-b+c, -a+b+c, a+b+c and  check if any
! of the four new vectors is shorter than any of a, b or c. If so, replace the
! longest of a, b and c with the new (shorter) vector. The resulting three
! vectors are also linearly independent.
! ----------------------------------------------------------------------
function reduce_vec(vecs) result(output)
  use utils_module, only : l2_norm
  implicit none
  
  real(dp), intent(inout) :: vecs(3,3)
  logical                 :: output
  
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
  output=.false.
  do i=1,4
    nlen = l2_norm(newvecs(i,:))
    if(nlen<maxlen)then
      vecs(longest,:) = newvecs(i,:)
      output = .true.
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
function find_hnf_supercell_matrix(qpoint,scaling) result(output)
  use utils_module, only : gcd
  use linear_algebra_module
  implicit none
  
  type(QpointData), intent(in) :: qpoint
  integer,          intent(in) :: scaling
  type(IntMatrix)              :: output ! S.
  
  integer :: scaled_qpoint(3)
  
  ! The elements of the scaled q-point.
  integer :: qx
  integer :: qy
  integer :: qz
  
  ! The determinant of S.
  integer :: sc_size
  
  ! The elements of S.
  integer :: s11,s12,s13
  integer ::     s22,s23
  integer ::         s33
  
  ! Split the q-point into qx,qy,qz
  scaled_qpoint = int(qpoint%scaled_qpoint)
  qx = scaled_qpoint(1)
  qy = scaled_qpoint(2)
  qz = scaled_qpoint(3)
  
  ! Calculate the determinant of the output matrix.
  sc_size = scaling / gcd(qx,qy,qz,scaling)
  
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
                if (all(modulo( int(output*qpoint%scaled_qpoint), &
                              & scaling) == 0)) then
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
  call print_line(CODE_ERROR//': no supercell matrix found for q-point '// &
     & qpoint%qpoint)
  call err()
end function

function find_reduced_supercell_matrix(qpoint,scaling,structure) &
   & result(output)
  implicit none
  
  type(QpointData),    intent(in) :: qpoint
  integer,             intent(in) :: scaling
  type(StructureData), intent(in) :: structure
  type(IntMatrix)                 :: output
  
  type(IntMatrix)  :: supercell_hnf
  type(RealMatrix) :: supercell_cart
  
  ! Find a supercell matrix S, s.t. S.q is a vector of integers.
  supercell_hnf = find_hnf_supercell_matrix(qpoint, scaling)
  
  ! Minkowski reduce the supercell in cartesian co-ordinates.
  ! Converts to cartesian, reduces, and converts back again.
  supercell_cart = supercell_hnf*structure%lattice
  supercell_cart = minkowski_reduce(supercell_cart)
  output         = nint(dble( supercell_cart &
                          & * transpose(structure%recip_lattice)))
  
  ! Check that S still transforms q to an integer vector.
  if (any(modulo( int(output*qpoint%scaled_qpoint), scaling) /= 0)) then
    call print_line(CODE_ERROR//': Error Minkowski reducing supercell.')
    call err()
  endif
end function

! ----------------------------------------------------------------------
! Checks supercells.
! ----------------------------------------------------------------------
! Throws an error if there is a problem.
subroutine check_supercells(supercells,structure)
  implicit none
  
  type(StructureData), intent(in) :: supercells(:)
  type(StructureData), intent(in) :: structure
  
  ! Working variables.
  type(StructureData)           :: supercell
  integer                       :: atom_1
  integer                       :: atom_2
  type(RealVector)              :: fractional_difference
  type(RealMatrix)              :: rotation_identity
  
  ! Temporary variables.
  integer :: i,j,k
  
  ! --------------------------------------------------
  ! Check that all supercells have been constructed correctly.
  ! --------------------------------------------------
  do i=1,size(supercells)
    supercell = supercells(i)
    
    ! Check that no atoms are on top of one another.
    do j=1,supercell%no_atoms
      do k=1,j-1
        if (l2_norm( supercell%atoms(j)%cartesian_position() &
                 & - supercell%atoms(k)%cartesian_position() )<0.5_dp) then
          call print_line(CODE_ERROR//': atoms '//k//' and '//j//' are within &
             &0.5 Bohr of one another in supercell '//i//'.')
          call err()
        endif
      enddo
    enddo
    
    ! Check that all atoms are related by an R-vector to their copies in the
    !    primitive cell.
    do atom_1=1,supercell%no_atoms
      atom_2 = supercell%atoms(atom_1)%prim_id()
      fractional_difference = structure%recip_lattice       &
         & * ( supercell%atoms(atom_1)%cartesian_position() &
         &   - supercell%atoms(atom_2)%cartesian_position())
      if ( l2_norm( fractional_difference             &
       &          - vec(nint(fractional_difference))) &
       & > 1.0e-10_dp) then
        call print_line(CODE_ERROR//': atoms '//atom_1//' and '//atom_2// &
           & ' in supercell '//i//' are not related by an R-vector')
        call err()
      endif
    enddo
    
    ! Check supercell rotations
    do j=1,size(supercell%symmetries)
      ! Check that R.R^T is the identity.
      rotation_identity = supercell%symmetries(j)%cartesian_rotation &
                      & * transpose(supercell%symmetries(j)%cartesian_rotation)
      if (any(abs(dble( rotation_identity &
                    & - make_identity_matrix(3)))>1.0e-10_dp)) then
        call print_line(ERROR//': rotation '//j//' in supercell '//i// &
           & ' is not a rotation or an improper rotation.')
        call print_line('Rotation matrix:')
        call print_line(supercell%symmetries(j)%cartesian_rotation)
        call print_line('')
        call print_line('This may be caused by small deviations from &
           &symmetry. Try adjusting the lattice vectors and atomic positions &
           &so that symmetry is exact.')
        call err()
      endif
    enddo
  enddo
end subroutine

! ----------------------------------------------------------------------
! Generates q-points.
! ----------------------------------------------------------------------
function generate_qpoints(structure,large_supercell) result(output)
  use linear_algebra_module
  use structure_module
  use group_module
  implicit none
  
  ! Inputs.
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: large_supercell
  type(QpointData), allocatable   :: output(:)
  
  ! IBZ q-point variables.
  integer,         allocatable :: paired_qpoints(:)
  
  ! Working variables
  type(IntVector) :: rotated_scaled_qpoint
  type(IntVector) :: scaled_supercell_gvector
  
  ! Temporary variables
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Construct q-points from G-vectors of large supercell.
  ! --------------------------------------------------
  allocate(output(large_supercell%sc_size), stat=ialloc); call err(ialloc)
  do i=1,large_supercell%sc_size
    output(i)%gvector = large_supercell%gvectors(i)
    output(i)%scaled_qpoint = transpose(large_supercell%recip_supercell) &
                          & * output(i)%gvector
    output(i)%qpoint = frac(output(i)%scaled_qpoint) / large_supercell%sc_size
    output(i)%to_simulate = .true.
  enddo
  
  ! --------------------------------------------------
  ! Find which q-points can be rotated onto other q-points.
  ! --------------------------------------------------
  ! N.B. q-points are in fractional reciprocal co-ordinates, so they
  !    transform by the transpose of the fractional rotation matrix.
  do_i : do i=1,size(output)
    do j=1,size(structure%symmetries)
      rotated_scaled_qpoint = transpose(structure%symmetries(j)%rotation) &
                          & * output(i)%scaled_qpoint
      do k=1,i-1
        if (rotated_scaled_qpoint==output(k)%scaled_qpoint) then
          output(i)%to_simulate = .false.
          cycle do_i
        endif
      enddo
    enddo
  enddo do_i
  
  ! --------------------------------------------------
  ! Find paired q-points.
  ! --------------------------------------------------
  ! qpoint + paired_qpoint = G, for a primitive-cell G-vector.
  allocate( paired_qpoints(large_supercell%sc_size), &
          & stat=ialloc); call err(ialloc)
  paired_qpoints = 0
  do i=1,size(output)
    do j=1,size(output)
      if (all(modulo( int(output(i)%scaled_qpoint+output(j)%scaled_qpoint), &
                    & large_supercell%sc_size)==0)) then
        if (paired_qpoints(i)==0 .and. paired_qpoints(j)==0) then
          paired_qpoints(i) = j
          paired_qpoints(j) = i
        else
          if (paired_qpoints(i)/=j .or. paired_qpoints(j)/=i) then
            call print_line(CODE_ERROR//': error pairing q-points.')
          endif
        endif
      endif
    enddo
  enddo
  
  if (any(paired_qpoints==0)) then
    call print_line(CODE_ERROR//': q-points were not succesfully paired up.')
    call err()
  endif
  
  ! --------------------------------------------------
  ! Use paired q-point information to decide which q-points to simulate.
  ! --------------------------------------------------
  do i=1,size(output)
    if (paired_qpoints(i)==i) then
      output(i)%is_paired_qpoint = .true.
    else
      output(i)%is_paired_qpoint = .false.
    endif
    
    output(i)%paired_qpoint = paired_qpoints(i)
    
    if (paired_qpoints(i)<i) then
      output(i)%to_simulate = .false.
    endif
  enddo
  
  ! --------------------------------------------------
  ! Find the smallest supercell matrix for each q-point.
  ! --------------------------------------------------
  do i=1,size(output)
    output(i)%supercell_matrix = find_reduced_supercell_matrix( &
                                     & output(i),               &
                                     & large_supercell%sc_size, &
                                     & structure)
    
    scaled_supercell_gvector = output(i)%supercell_matrix &
                           & * output(i)%scaled_qpoint
    
    if (any(modulo( int(scaled_supercell_gvector), &
                  & large_supercell%sc_size)/=0)) then
      call print_line(CODE_ERROR//': q-point and supercell do not match.')
      call err()
    endif
    
    output(i)%supercell_gvector = scaled_supercell_gvector &
                              & / large_supercell%sc_size
  enddo
end function

! ----------------------------------------------------------------------
! Generates supercells.
! ----------------------------------------------------------------------
function generate_supercells(structure,qpoints,symmetry_precision) &
   & result(output)
  use linear_algebra_module
  use structure_module
  use group_module
  implicit none
  
  ! Inputs.
  type(StructureData), intent(in)  :: structure
  type(QpointData),    intent(in)  :: qpoints(:)
  real(dp),            intent(in)  :: symmetry_precision
  type(StructureData), allocatable :: output(:)
  
  ! Supercell variables.
  integer :: no_supercells
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  allocate(output(size(qpoints)), stat=ialloc); call err(ialloc)
  
  no_supercells = 0
  do_i : do i=1,size(qpoints)
    ! Only make supercells for those q-points which are to be simulated.
    if (qpoints(i)%to_simulate) then
      
      ! Check if a supercell suitable for this q-point already exists.
      do j=1,i-1
        if (qpoints(j)%supercell_matrix==qpoints(i)%supercell_matrix) then
          cycle do_i
        endif
      enddo
      
      ! Construct a supercell for this q-point.
      no_supercells = no_supercells+1
      output(no_supercells) = construct_supercell( &
                    & structure,                   &
                    & qpoints(i)%supercell_matrix, &
                    & symmetry_precision=symmetry_precision)
    endif
  enddo do_i
  
  output = output(:no_supercells)
  
  call check_supercells(output,structure)
end function
end module
