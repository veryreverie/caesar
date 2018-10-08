! ======================================================================
! Methods for generating supercells.
! ======================================================================
module supercell_module
  use utils_module
  
  use qpoint_module
  use basic_structure_module
  use structure_data_module
  implicit none
  
  private
  
  public :: construct_supercell
  public :: construct_supercell_matrix
  
  interface construct_supercell_matrix
    module procedure construct_supercell_matrix_qpoint
    module procedure construct_supercell_matrix_qpoints
  end interface
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
! Find a supercell matrix, S, s.t. S.q is a vector of integers for all q.
! ----------------------------------------------------------------------
! S is found s.t. |S| is as small as possible (whilst being >0).
! Returns answer in Hermite Normal Form.
function find_hnf_supercell_matrix(qpoints) result(output)
  implicit none
  
  type(QpointData), intent(in) :: qpoints(:)
  type(IntMatrix)              :: output
  
  type(IntFraction), allocatable :: qpoint_components(:,:)
  
  ! The elements of S (s21=s31=s32=0).
  integer :: s11, s12, s13
  integer ::      s22, s23
  integer ::           s33
  
  integer :: output_matrix(3,3)
  
  integer :: s11_max,s22_max,s33_max
  
  integer :: i,ialloc
  
  output_matrix = 0
  
  allocate(qpoint_components(3,size(qpoints)), stat=ialloc); call err(ialloc)
  do i=1,size(qpoints)
    qpoint_components(:,i) = frac(qpoints(i)%qpoint)
  enddo
  
  ! Calculate an upper bound on the diagonal elements.
  ! (This is the case where s12=s13=s23=0).
  s11_max = lcm(qpoint_components(1,:)%denominator())
  s22_max = lcm(qpoint_components(2,:)%denominator())
  s33_max = lcm(qpoint_components(3,:)%denominator())
  
  ! s33 is simply s33_max, since there are no matrix elements to the right
  !    of s33.
  s33 = s33_max
  output_matrix(3,:) = [0, 0, s33]
  
  ! Loop over s22 and s23 to find the smallest s22.
  do_s22 : do s22=1,s22_max
    do s23=0,s33-1
      if (all(is_int( s22*qpoint_components(2,:) &
                  & + s23*qpoint_components(3,:) ))) then
        output_matrix(2,:) = [0, s22, s23]
        exit do_s22
      endif
    enddo
  enddo do_s22
  
  if (output_matrix(2,2)==0) then
    call print_line(ERROR//': Could not find supercell matrix.')
    call err()
  endif
  
  ! Loop over s11, s12 and s13 to find the smallest s11.
  do_s11 : do s11=1,s11_max
    do s12=0,s22-1
      do s13=0,s33-1
        if (all(is_int( s11*qpoint_components(1,:) &
                    & + s12*qpoint_components(2,:) &
                    & + s13*qpoint_components(3,:) ))) then
          output_matrix(1,:) = [s11, s12, s13]
          exit do_s11
        endif
      enddo
    enddo
  enddo do_s11
  
  if (output_matrix(1,1)==0) then
    call print_line(ERROR//': Could not find supercell matrix.')
    call err()
  endif
  
  ! Transfer output_matrix to output.
  output = output_matrix
  
  ! Check that S correctly transforms all q-points to integer vectors.
  do i=1,size(qpoints)
    if (.not. is_int(output*qpoints(i)%qpoint)) then
      call print_line(CODE_ERROR//': Supercell matrix does not correctly &
         &transform q-points.')
      call err()
    endif
  enddo
end function

! ----------------------------------------------------------------------
! Construct the smallest matrix S s.t. S.q is a G-vector for all q.
! Smallest means:
!    - |S| is as small as possible.
!    - The vectors defined in fractional co-ordinates by the rows of S
!         are as short as possible in cartesian co-ordinates.
! ----------------------------------------------------------------------
function construct_supercell_matrix_qpoints(qpoints,structure) result(output)
  implicit none
  
  type(QpointData),    intent(in) :: qpoints(:)
  type(StructureData), intent(in) :: structure
  type(IntMatrix)                 :: output
  
  type(IntMatrix)  :: supercell_hnf
  
  integer :: i
  
  ! Find the matrix S, s.t. S.q is a vector of integers.
  ! S is found in Hermite Normal Form.
  supercell_hnf = find_hnf_supercell_matrix(qpoints)
  
  ! Reduce the matrix, s.t. the vectors defined by its rows are as short as
  !    possible when transformed into cartesian co-ordinates.
  output = reduce(supercell_hnf, structure)
  
  ! Check that S still transforms q to an integer vector.
  do i=1,size(qpoints)
    if (.not. is_int(output*qpoints(i)%qpoint)) then
      call print_line(CODE_ERROR//': Error Minkowski reducing supercell.')
      call err()
    endif
  enddo
end function

function construct_supercell_matrix_qpoint(qpoint,structure) result(output)
  implicit none
  
  type(QpointData),    intent(in) :: qpoint
  type(StructureData), intent(in) :: structure
  type(IntMatrix)                 :: output
  
  output = construct_supercell_matrix([qpoint], structure)
end function

! ----------------------------------------------------------------------
! Checks supercells.
! ----------------------------------------------------------------------
! Throws an error if there is a problem.
impure elemental subroutine check_supercell(supercell,structure)
  implicit none
  
  type(StructureData), intent(in) :: supercell
  type(StructureData), intent(in) :: structure
  
  ! Working variables.
  integer               :: atom_1
  integer               :: atom_2
  type(RealVector)      :: fractional_difference
  type(RealMatrix)      :: tensor_identity
  real(dp), allocatable :: symmetry_errors(:)
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  ! Check that no atoms are on top of one another.
  do i=1,supercell%no_atoms
    do j=1,i-1
      if (l2_norm( supercell%atoms(i)%cartesian_position() &
               & - supercell%atoms(j)%cartesian_position() )<0.1_dp) then
        call print_line(WARNING//': atoms '//j//' and '//i//' are within &
           &0.1 Bohr of one another.')
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
         & ' in supercell are not related by an R-vector')
      call err()
    endif
  enddo
  
  ! Check supercell transformations.
  allocate( symmetry_errors(size(supercell%symmetries)), &
          & stat=ialloc); call err(ialloc)
  symmetry_errors = 0
  do i=1,size(supercell%symmetries)
    ! Check that R.R^T is the identity.
    tensor_identity = supercell%symmetries(i)%cartesian_tensor &
                  & * transpose(supercell%symmetries(i)%cartesian_tensor)
    symmetry_errors(i) = sqrt(sum_squares( tensor_identity         &
                                       & - make_identity_matrix(3) ))
  enddo
  if (any(symmetry_errors>1e-10_dp)) then
    call print_line('')
    call print_line(WARNING//': At least one symmetry in supercell '// &
       & ' is not an orthogonal transformation.')
    call print_line('Largest L2 deviation in T^T.T=I: '// &
       & maxval(symmetry_errors))
    call print_line('This may be caused by small deviations from &
       &symmetry. Try adjusting the lattice vectors and atomic positions &
       &so that symmetry is exact.')
  endif
end subroutine

! ----------------------------------------------------------------------
! Calculates the set of vectors which are not related to one another by
!    lattice vectors.
! ----------------------------------------------------------------------
! Either calculates the R-vectors of the primitive cell which are unique
!    in the supercell, or the G-vectors of the reciprocal supercell which are
!    unique in the reciprocal primitive cell.
function calculate_unique_vectors(lattice,centre_on_origin) result(output)
  implicit none
  
  ! Lattice vectors are the rows of lattice.
  type(IntMatrix), intent(in)  :: lattice
  logical,         intent(in)  :: centre_on_origin
  type(IntVector), allocatable :: output(:)
  
  integer         :: lattice_size
  integer         :: no_vectors
  type(IntVector) :: frac_vec
  type(IntVector) :: prim_vec
  integer         :: i,j,k,ialloc
  
  if (size(lattice,1)/=3 .or. size(lattice,2)/=3) then
    call print_line('Error: lattice is not 3x3.')
    call err()
  endif
  
  lattice_size = abs(determinant(lattice))
  
  allocate(output(lattice_size), stat=ialloc); call err(ialloc)
  
  no_vectors = 0
  do i=0,lattice_size-1
    do j=0,lattice_size-1
      do k=0,lattice_size-1
        
        ! Construct vectors in scaled fractional primitive co-ordinates.
        ! (scaled by lattice_size, to preserve integer representation).
        if (centre_on_origin) then
          frac_vec = [i-lattice_size/2,j-lattice_size/2,k-lattice_size/2]
        else
          frac_vec = [i,j,k]
        endif
        
        ! Transform to scaled fractional supercell co-ordinates.
        prim_vec = transpose(lattice)*frac_vec
        
        ! Check if the scaled co-ordinate is scaling*(integer co-ordinate).
        if (all(modulo(int(prim_vec),lattice_size)==0)) then
          no_vectors = no_vectors+1
          output(no_vectors) = prim_vec / lattice_size
        endif
      enddo
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Makes a supercell from a given primitive cell and supercell matrix.
! ----------------------------------------------------------------------
function construct_supercell(structure,supercell_matrix) result(output)
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(IntMatrix),     intent(in) :: supercell_matrix
  type(StructureData)             :: output
  
  ! R-vector and G-vector information.
  type(IntVector), allocatable :: rvectors(:)
  type(IntVector), allocatable :: gvectors(:)
  
  ! Atomic information.
  type(String),     allocatable :: species(:)
  type(String),     allocatable :: species2(:)
  real(dp),         allocatable :: masses(:)
  real(dp),         allocatable :: masses2(:)
  type(RealVector), allocatable :: positions(:,:)
  type(RealVector), allocatable :: positions2(:)
  integer,          allocatable :: atom_rvector_ids(:)
  integer,          allocatable :: atom_prim_ids(:)
  
  ! Temporary variables
  integer :: sc_size
  
  integer :: atom,prim,rvec,ialloc
  
  ! Generate supercell and lattice data.
  sc_size = abs(determinant(supercell_matrix))
  
  ! Construct R-vectors and G-vectors.
  rvectors = calculate_unique_vectors(supercell_matrix, .false.)
  gvectors = calculate_unique_vectors(transpose(supercell_matrix), .true.)
  
  ! Check that rvectors(1) is the Gamma-point.
  if (rvectors(1)/=vec([0,0,0])) then
    call print_line(CODE_ERROR//': The first R-vector is not the Gamma point.')
    call err()
  endif
  
  ! Construct atomic information.
  allocate( species(structure%no_atoms),                  &
          & species2(structure%no_atoms*sc_size),         &
          & masses(structure%no_atoms),                   &
          & masses2(structure%no_atoms*sc_size),          &
          & positions(structure%no_atoms,sc_size),        &
          & positions2(structure%no_atoms*sc_size),       &
          & atom_rvector_ids(structure%no_atoms*sc_size), &
          & atom_prim_ids(structure%no_atoms*sc_size),    &
          & stat=ialloc); call err(ialloc)
  do prim=1,structure%no_atoms
    species(prim) = structure%atoms(prim)%species()
    masses(prim) = structure%atoms(prim)%mass()
    do rvec=1,sc_size
      positions(prim,rvec) = transpose(structure%lattice)                  &
                         & * ( structure%atoms(prim)%fractional_position() &
                         &   + rvectors(rvec))
    enddo
  enddo
  
  atom = 0
  do rvec=1,sc_size
    do prim=1,structure%no_atoms
      atom = atom+1
      species2(atom) = species(prim)
      masses2(atom) = masses(prim)
      positions2(atom) = positions(prim,rvec)
      atom_rvector_ids(atom) = rvec
      atom_prim_ids(atom) = prim
    enddo
  enddo
  
  ! Construct output.
  output = StructureData(                                                     &
     & basic_structure = BasicStructure( supercell_matrix*structure%lattice,  &
     &                                   species2,                            &
     &                                   masses2,                             &
     &                                   positions2),                         &
     & basic_supercell = BasicSupercell(                                      &
     &         supercell_matrix,                                              &
     &         calculate_unique_vectors(supercell_matrix, .false.),           &
     &         calculate_unique_vectors(transpose(supercell_matrix), .true.), &
     &         atom_rvector_ids,                                              &
     &         atom_prim_ids))
  
  ! Check output.
  call check_supercell(output, structure)
end function
end module
