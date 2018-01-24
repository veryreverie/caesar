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
  
  public :: generate_supercells
contains

! ----------------------------------------------------------------------
! Given three linearly independent input vectors a, b and c, construct the
! following linear combinations: a+b-c, a-b+c, -a+b+c, a+b+c and  check if any
! of the four new vectors is shorter than any of a, b or c. If so, replace the
! longest of a, b and c with the new (shorter) vector. The resulting three
! vectors are also linearly independent.
! ----------------------------------------------------------------------
function reduce(input,metric) result(output)
  use linear_algebra_module
  implicit none
  
  type(IntVector),  intent(in) :: input(3)
  type(RealMatrix), intent(in) :: metric
  type(IntVector)              :: output(3)
  
  type(IntVector) :: vectors(13)
  
  real(dp) :: lengths(13)
  logical  :: chosen(13)
  
  integer  :: i,j
  
  ! Construct the relevant linear combinations of input vectors.
  vectors(1:3) =  input
  vectors(4)   =  input(1) + input(2)
  vectors(5)   =  input(2) + input(3)
  vectors(6)   =  input(3) + input(1)
  vectors(7)   =  input(1) - input(2)
  vectors(8)   =  input(2) - input(3)
  vectors(9)   =  input(3) - input(1)
  vectors(10)  =  input(1) + input(2) + input(3)
  vectors(11)  =  input(1) + input(2) - input(3)
  vectors(12)  =  input(2) + input(3) - input(1)
  vectors(13)  =  input(3) + input(1) - input(2)
  
  ! Calculate the lengths of the vectors in cartesian co-ordinates.
  ! N.B. all lengths are stored as length**2, to avoid square roots.
  do i=1,13
    lengths(i) = input(i) * metric * input(i)
  enddo
  
  ! Select the shortest three vectors.
  chosen = .false.
  do i=1,3
    j=minloc(lengths, dim=1, mask=.not. chosen)
    output(i) = vectors(j)
    chosen(i) = .true.
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
function minkowski_reduce(input,structure) result(output)
  use utils_module, only : sum_squares
  use linear_algebra_module
  implicit none
  
  type(IntMatrix),     intent(in) :: input
  type(StructureData), intent(in) :: structure
  type(IntMatrix)                 :: output
  
  type(RealMatrix) :: metric
  
  type(IntVector) :: vectors(3)
  type(IntVector) :: reduced(3)
  
  integer :: matrix(3,3)
  
  integer :: i
  
  ! Construct the metric, M, such that v.M.v is the length of v
  !    in cartesian co-ordinates.
  metric = structure%lattice * transpose(structure%lattice)
  
  ! Construct the three supercell lattice vectors in fractional co-ordinates.
  matrix = int(input)
  do i=1,3
    vectors(i) = matrix(i,:)
  enddo
  
  ! Minkowski reduce the three vectors,
  !   i.e. find the shortest linear combinations of the three.
  do
    reduced = reduce(vectors, metric)
    if (all(reduced==vectors)) then
      exit
    endif
    vectors = reduced
  enddo
  
  ! Construct output.
  do i=1,3
    matrix(i,:) = int(vectors(i))
  enddo
  output = matrix
end function

! ----------------------------------------------------------------------
! Find a supercell matrix, S, s.t. S.q is a vector of integers.
! ----------------------------------------------------------------------
! Returns answer in Hermite Normal Form.
function find_hnf_supercell_matrix(qpoint) result(output)
  use linear_algebra_module
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

function find_reduced_supercell_matrix(qpoint,structure) result(output)
  implicit none
  
  type(QpointData),    intent(in) :: qpoint
  type(StructureData), intent(in) :: structure
  type(IntMatrix)                 :: output
  
  type(IntMatrix)  :: supercell_hnf
  
  ! Find a supercell matrix S, s.t. S.q is a vector of integers.
  supercell_hnf = find_hnf_supercell_matrix(qpoint)
  
  ! Minkowski reduce the supercell in cartesian co-ordinates.
  ! Converts to cartesian, reduces, and converts back again.
  output = minkowski_reduce(supercell_hnf, structure)
  
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
  
  ! q-point variables.
  logical, allocatable :: accounted_for(:)
  integer, allocatable :: min_sc_sizes(:)
  type(FractionVector) :: rotated_qpoint
  
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
            rotated_qpoint = transpose(structure%symmetries(k)%rotation) &
                         & * qpoints(j)%qpoint
            do l=1,size(qpoints)
              if (qpoints(l)%qpoint == rotated_qpoint) then
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
