! ======================================================================
! Generates the supercells needed to simulate phonons at all q-points.
! ======================================================================
module generate_supercells_module
  use common_module
  implicit none
  
  private
  
  public :: generate_supercells
contains

! ----------------------------------------------------------------------
! Takes three linearly independent (LI) vectors, a, b and c.
! The longest vector is replaced with the shortest linear combination
!    of the three s.t. the three output vectors are still LI.
! The linear combinations considered are a+b+c, a+b-c, a-b+c and -a+b+c.
! ----------------------------------------------------------------------
function reduce_once(input,metric) result(output)
  implicit none
  
  type(IntVector),  intent(in) :: input(3)
  type(RealMatrix), intent(in) :: metric
  type(IntVector)              :: output(3)
  
  type(IntVector) :: combinations(4)
  real(dp)        :: vector_lengths(3)
  real(dp)        :: combination_lengths(4)
  
  integer :: i,j
  
  ! Construct the linear combinations a+b+c, a+b-c, b+c-a and c+a-b.
  combinations = [ input(1)+input(2)+input(3), &
                 & input(1)+input(2)-input(3), &
                 & input(2)+input(3)-input(1), &
                 & input(3)+input(1)-input(2)  ]
  
  ! Find the lengths of the input vectors and linear combinations.
  ! Lengths are in terms of cartesian distances, and stored as
  !    length squared, since taking the square root is unnecessary.
  do i=1,3
    vector_lengths(i) = input(i) * metric * input(i)
  enddo
  
  do i=1,4
    combination_lengths(i) = combinations(i) * metric * combinations(i)
  enddo
  
  ! Identify the longest input vector, and the shortest linear combination.
  i = maxloc(vector_lengths, 1)
  j = minloc(combination_lengths, 1)
  
  ! Replace said input vector with said linear combination,
  !    if it is shorter.
  output = input
  if (combination_lengths(j)<vector_lengths(i)) then
    output(i) = combinations(j)
  endif
end function

! ----------------------------------------------------------------------
! Takes a lattice matrix, L, whose rows are three vectors, a, b and c.
! Returns a matrix whose rows are linear combinations of a, b and c s.t.:
!    - |A|/=0, i.e. the rows or the output are linearly independent.
!    - Each vector is shorter than any linear combination of those vectors.
! ----------------------------------------------------------------------
function reduce(input,structure) result(output)
  implicit none
  
  type(IntMatrix),     intent(in) :: input
  type(StructureData), intent(in) :: structure
  type(IntMatrix)                 :: output
  
  type(RealMatrix) :: metric
  
  type(IntVector) :: vectors(3)
  type(IntVector) :: reduced(3)
  
  integer :: matrix(3,3)
  
  integer :: i
  
  ! Construct the metric, M, such that v.M.v is the length squared of v
  !    in cartesian co-ordinates.
  ! Since v in cartesian co-ordinates is L^T . v, M=LL^T.
  metric = structure%lattice * transpose(structure%lattice)
  
  ! Construct the three supercell lattice vectors in fractional co-ordinates.
  matrix = int(input)
  do i=1,3
    vectors(i) = matrix(i,:)
  enddo
  
  ! Reduce the three vectors.
  iter : do
    ! Check linear combinations of two vectors.
    do i=1,3
      reduced = vectors
      reduced(i) = zeroes(3)
      reduced = reduce_once(vectors, metric)
      reduced(i) = vectors(i)
      if (any(reduced/=vectors)) then
        vectors = reduced
        cycle iter
      endif
    enddo
    
    ! Check linear combinations of all three vectors.
    reduced = reduce_once(vectors, metric)
    if (all(reduced==vectors)) then
      exit
    endif
    vectors = reduced
  enddo iter
  
  ! Construct output.
  do i=1,3
    matrix(i,:) = int(vectors(i))
  enddo
  output = matrix
end function

! ----------------------------------------------------------------------
! Find a supercell matrix, S, s.t. S.q is a vector of integers.
! ----------------------------------------------------------------------
! S is found s.t. |S| is as small as possible (whilst being >0).
! Returns answer in Hermite Normal Form.
function find_hnf_supercell_matrix(qpoint) result(output)
  implicit none
  
  type(QpointData), intent(in) :: qpoint
  type(IntMatrix)              :: output ! S.
  
  ! The determinant of S.
  integer :: sc_size
  
  ! The elements of S.
  integer :: s11,s12,s13
  integer ::     s22,s23
  integer ::         s33
  
  ! Calculate the determinant of the output matrix.
  sc_size = qpoint%min_sc_size()
  
  ! Loop over all matrices in Hermite Normal Form, defined as
  !    s11*s22*s33 = |S| = sc_size.
  !    s21=s31=s32=0
  !    0 <= s12 < s22
  !    0 <= s13 < s33
  !    0 <= s23 < s33
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
                
                ! Check if S.q is a vector of integers.
                if (is_int(output*qpoint%qpoint)) then
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

! ----------------------------------------------------------------------
! Construct the smallest matrix S s.t. S.q is a G-vector.
! Smallest means:
!    - |S| is as small as possible.
!    - The vectors defined in fractional co-ordinates by the rows of S
!         are as short as possible in cartesian co-ordinates.
! ----------------------------------------------------------------------
function find_reduced_supercell_matrix(qpoint,structure) result(output)
  implicit none
  
  type(QpointData),    intent(in) :: qpoint
  type(StructureData), intent(in) :: structure
  type(IntMatrix)                 :: output
  
  type(IntMatrix)  :: supercell_hnf
  
  ! Find the matrix S, s.t. S.q is a vector of integers.
  ! S is found in Hermite Normal Form.
  supercell_hnf = find_hnf_supercell_matrix(qpoint)
  
  ! Reduce the matrix, s.t. the vectors defined by its rows are as short as
  !    possible when transformed into cartesian co-ordinates.
  output = reduce(supercell_hnf, structure)
  
  ! Check that S still transforms q to an integer vector.
  if (.not. is_int(output*qpoint%qpoint)) then
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
  type(StructureData)   :: supercell
  integer               :: atom_1
  integer               :: atom_2
  type(RealVector)      :: fractional_difference
  type(RealMatrix)      :: rotation_identity
  real(dp), allocatable :: symmetry_errors(:)
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Check that all supercells have been constructed correctly.
  ! --------------------------------------------------
  do i=1,size(supercells)
    supercell = supercells(i)
    
    ! Check that no atoms are on top of one another.
    do j=1,supercell%no_atoms
      do k=1,j-1
        if (l2_norm( supercell%atoms(j)%cartesian_position() &
                 & - supercell%atoms(k)%cartesian_position() )<0.1_dp) then
          call print_line(WARNING//': atoms '//k//' and '//j//' are within &
             &0.1 Bohr of one another in supercell '//i//'.')
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
    allocate( symmetry_errors(size(supercell%symmetries)), &
            & stat=ialloc); call err(ialloc)
    symmetry_errors = 0
    do j=1,size(supercell%symmetries)
      ! Check that R.R^T is the identity.
      rotation_identity = supercell%symmetries(j)%cartesian_rotation &
                      & * transpose(supercell%symmetries(j)%cartesian_rotation)
      symmetry_errors(j) = sqrt(sum_squares( rotation_identity &
                                         & - make_identity_matrix(3)))
    enddo
    if (any(symmetry_errors>1e-10_dp)) then
      call print_line('')
      call print_line(WARNING//': At least one symmetry in supercell '//i// &
         & ' is not an orthogonal transformation.')
      call print_line('Largest L2 deviation in R^T.R=I: '// &
         & maxval(symmetry_errors))
      call print_line('Please check symmetries in &
         &Supercell_'//i//'/structure.dat')
      call print_line('This may be caused by small deviations from &
         &symmetry. Try adjusting the lattice vectors and atomic positions &
         &so that symmetry is exact.')
    endif
    deallocate(symmetry_errors, stat=ialloc); call err(ialloc)
  enddo
end subroutine

! ----------------------------------------------------------------------
! Generates supercells.
! ----------------------------------------------------------------------
function generate_supercells(structure,qpoints,symmetry_precision) &
   & result(output)
  implicit none
  
  ! Inputs.
  type(StructureData), intent(in)  :: structure
  type(QpointData),    intent(in)  :: qpoints(:)
  real(dp),            intent(in)  :: symmetry_precision
  type(StructureData), allocatable :: output(:)
  
  ! q-point variables.
  logical, allocatable :: accounted_for(:)
  integer, allocatable :: min_sc_sizes(:)
  type(QpointData)     :: rotated_qpoint
  
  ! Supercell variables.
  integer         :: no_supercells
  type(IntMatrix) :: supercell_matrix
  
  ! Temporary variables.
  integer :: i,j,k,l,ialloc
  
  allocate( accounted_for(size(qpoints)), &
          & min_sc_sizes(size(qpoints)),  &
          & output(size(qpoints)),        &
          & stat=ialloc); call err(ialloc)
  
  ! Find the minimum size of the supercell required to simulate each q-point.
  do i=1,size(qpoints)
    min_sc_sizes(i) = qpoints(i)%min_sc_size()
  enddo
  
  ! Find the minimal set of supercells required to simulate every q-point.
  accounted_for = .false.
  no_supercells = 0
  do while (any(.not. accounted_for))
    ! Find the q-point with the largest min supercell size
    !    which is not yet accounted for.
    i = maxloc(min_sc_sizes, dim=1, mask=.not. accounted_for)
    
    ! Create a supercell for simulating this q-point.
    no_supercells = no_supercells+1
    supercell_matrix = find_reduced_supercell_matrix(qpoints(i), structure)
    output(no_supercells) = construct_supercell( &
                             & structure,        &
                             & supercell_matrix, &
                             & symmetry_precision=symmetry_precision)
    
    ! Find all q-points which can be simulated using this supercell.
    do j=1,size(qpoints)
      if (.not. accounted_for(j)) then
        if (is_int(output(no_supercells)%supercell * qpoints(j)%qpoint)) then
          ! Mark the q-point and its paired q-point as accounted for.
          accounted_for(j) = .true.
          accounted_for(qpoints(j)%paired_qpoint) = .true.
          
          ! Find any symmetrically equivalent q-points,
          !    and mark them as accounted for.
          do k=1,size(structure%symmetries)
            rotated_qpoint = structure%symmetries(k) * qpoints(j)
            do l=1,size(qpoints)
              if (qpoints(l) == rotated_qpoint) then
                ! Mark the q-point and its paired q-point as accounted for.
                accounted_for(l) = .true.
                accounted_for(qpoints(l)%paired_qpoint) = .true.
              endif
            enddo
          enddo
        endif
      endif
    enddo
    
    ! Check that this supercell does simulate this q-point.
    if (.not. accounted_for(i)) then
      call print_line(CODE_ERROR//': Generated supercell does not match &
         &q-point.')
      call err()
    endif
  enddo
  
  output = output(:no_supercells)
  
  call check_supercells(output,structure)
end function
end module
