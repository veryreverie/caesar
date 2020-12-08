! ======================================================================
! Calculate dynamical matrices at each q-point from the Hessians of the
!    non-diagonal supercells.
! ======================================================================
module calculate_dynamical_matrices_module
  use utils_module
  
  use structure_module
  use normal_mode_module
  
  use dynamical_matrix_module
  implicit none
  
  private
  
  public :: MatrixAndModes
  public :: calculate_dynamical_matrices
  
  ! Result type for calculate_dynamical_matrices.
  type :: MatrixAndModes
    type(DynamicalMatrix)          :: matrix
    type(ComplexMode), allocatable :: modes(:)
  end type
  
  interface calculate_dynamical_matrices
    module procedure calculate_dynamical_matrices_hessians
  end interface
contains

! Calculate dynamical matrices at each q-point from the Hessians of each
!    non-diagonal supercell.
function calculate_dynamical_matrices_hessians(structure,supercells, &
   & supercell_hessians,qpoints,logfile) result(output)
  implicit none
  
  type(StructureData),    intent(in)    :: structure
  type(StructureData),    intent(in)    :: supercells(:)
  type(CartesianHessian), intent(in)    :: supercell_hessians(:)
  type(QpointData),       intent(in)    :: qpoints(:)
  type(OFile),            intent(inout) :: logfile
  type(MatrixAndModes), allocatable     :: output(:)
  
  integer, allocatable :: qpoint_symmetries(:,:)
  
  logical, allocatable :: modes_calculated(:)
  
  integer :: subspace_id
  
  type(QpointData) :: transformed_qpoint
  
  integer :: i,j,k,ialloc
  
  ! Calculate the symmetry mapping between q-points.
  ! If qpoint_symmetries(i,j)==k, then symmetry k maps qpoint i to qpoint j.
  ! If qpoint_symmetries(i,j)==0, then there is no symmetry mapping i to j.
  allocate( qpoint_symmetries(size(qpoints), size(qpoints)), &
          & stat=ialloc); call err(ialloc)
  qpoint_symmetries = 0
  do i=1,size(qpoints)
    do k=1,size(structure%symmetries)
      transformed_qpoint = structure%symmetries(k) * qpoints(i)
      j = first(qpoints==transformed_qpoint, default=0)
      if (j/=0) then
        if (qpoint_symmetries(i,j)==0) then
          qpoint_symmetries(i,j) = k
        endif
      endif
    enddo
  enddo
  
  ! --------------------------------------------------
  ! Calculate the dynamical matrix and normal modes at each q-point.
  ! --------------------------------------------------
  allocate(output(size(qpoints)), stat=ialloc); call err(ialloc)
  modes_calculated = [(.false.,i=1,size(qpoints))]
  subspace_id = 1
  iter : do while (.not. all(modes_calculated))
    ! ------------------------------
    ! First, check if there is an un-calculated q-point whose conjugate has
    !    already been calculated.
    ! ------------------------------
    do i=1,size(qpoints)
      if (modes_calculated(i)) then
        cycle
      endif
      
      j = first(qpoints%id==qpoints(i)%paired_qpoint_id)
      
      if (modes_calculated(j)) then
        call logfile%print_line('Constructing dynamical matrix and &
           &normal modes at q-point '//i//' as the conjugate of those at &
           &q-point '//j)
        
        ! The dynamical matrix at G-q is the complex conjugate of that at q.
        ! N.B. the complex conjugate, not the Hermitian conjugate.
        output(i)%matrix = DynamicalMatrix(                &
           & elements = conjg(output(j)%matrix%elements()) )
        
        ! The displacements of the normal modes at G-q are
        !    the complex conjugates of those at q.
        output(i)%modes = conjg(output(j)%modes)
        
        modes_calculated(i) = .true.
        call print_line(count(modes_calculated)//' of '//size(qpoints)//' &
           &dynamical matrices calculated.')
        cycle iter
      endif
    enddo
    
    ! ------------------------------
    ! Second, check if there is a q-point which is related by symmetry to
    !    an already calculated q-point.
    ! ------------------------------
    do i=1,size(qpoints)
      if (modes_calculated(i)) then
        cycle
      endif
      
      j = first(qpoint_symmetries(:,i)/=0, mask=modes_calculated, default=0)
      if (j/=0) then
        k = qpoint_symmetries(j,i)
        call logfile%print_line('Constructing dynamical matrix and &
           &normal modes at q-point '//i//' using symmetry from those at &
           &q-point '//j)
        output(i)%matrix = transform_modes( output(j)%matrix,        &
                                          & structure%symmetries(k), &
                                          & qpoints(j),              &
                                          & qpoints(i)               )
        output(i)%modes = transform( output(j)%modes,         &
                                   & structure%symmetries(k), &
                                   & qpoints(j),              &
                                   & qpoints(i)               )
        modes_calculated(i) = .true.
        call print_line(count(modes_calculated)//' of '//size(qpoints)//' &
           &dynamical matrices calculated.')
        cycle iter
      endif
    enddo
    
    ! ------------------------------
    ! Finally, check if there is a q-point which is a G-vector of any of the
    !    calculated supercells.
    ! ------------------------------
    do i=1,size(qpoints)
      if (modes_calculated(i)) then
        cycle
      endif
      
      do j=1,size(supercells)
        ! Check if q-point i is a G-vector of supercell j.
        if (is_int(supercells(j)%supercell * qpoints(i)%qpoint)) then
          call logfile%print_line('Constructing dynamical matrix and &
             &normal modes at q-point '//i//' directly from calculated force &
             &constants.')
          output(i) = construct_dynamical_matrix( qpoints(i),         &
                                                & supercells,         &
                                                & supercell_hessians, &
                                                & structure,          &
                                                & subspace_id,        &
                                                & logfile             )
          subspace_id = maxval(output(i)%modes%subspace_id) + 1
          modes_calculated(i) = .true.
          call print_line(count(modes_calculated)//' of '//size(qpoints)//' &
             &dynamical matrices calculated.')
          cycle iter
        endif
      enddo
    enddo
    
    ! ------------------------------
    ! If no q-point can be constructed, throw an error.
    ! ------------------------------
    call print_line(ERROR//': Insufficient data to construct dynamical &
       &matrices at all q-points.')
    call err()
  enddo iter
  
  ! --------------------------------------------------
  ! Check all dynamical matrices.
  ! --------------------------------------------------
  ! Run basic checks on each matrix in turn.
  call print_line('Checking dynamical matrices.')
  do i=1,size(qpoints)
    call output(i)%matrix%check(structure, logfile)
  enddo
end function

! --------------------------------------------------
! Construct the dynamical matrix at the specified q-point from the given
!    Hessian.
! Considers all supercell, S, and symmetry, R, combinations such that
!    S.R.q is a vector of integers, i.e. the q-point is a G-vector of the
!    supercell.
! --------------------------------------------------
function construct_dynamical_matrix(qpoint,supercells,hessian, &
   & structure,subspace_id,logfile) result(output)
  implicit none
  
  type(QpointData),       intent(in)    :: qpoint
  type(StructureData),    intent(in)    :: supercells(:)
  type(CartesianHessian), intent(in)    :: hessian(:)
  type(StructureData),    intent(in)    :: structure
  integer,                intent(in)    :: subspace_id
  type(OFile),            intent(inout) :: logfile
  type(MatrixAndModes)                  :: output
  
  type(QpointData) :: q_prime
  
  logical, allocatable :: is_copy(:,:)
  integer, allocatable :: sym_id(:)
  integer, allocatable :: sup_id(:)
  integer              :: no_copies
  integer              :: i,j,k,ialloc
  
  type(DynamicalMatrix)            :: dynamical_matrix
  type(ComplexMatrix), allocatable :: copy_matrices(:,:,:)
  type(ComplexMatrix), allocatable :: matrices(:,:)
  type(ComplexMode),   allocatable :: complex_modes(:)
  
  logical                          :: check_matrices
  type(ComplexMatrix), allocatable :: average(:,:)
  real(dp),            allocatable :: averages(:,:)
  real(dp),            allocatable :: differences(:,:)
  real(dp)                         :: l2_error
  
  ! Check that the supercells and Hessians correspond to one another.
  if (size(supercells)/=size(hessian)) then
    call print_line(CODE_ERROR//': Force constants and supercells do not &
       &match.')
    call err()
  endif
  
  ! Identify how many parings of supercell and symmetry
  !    allow q to be simulated.
  allocate( is_copy(size(supercells),size(structure%symmetries)), &
          & stat=ialloc); call err(ialloc)
  is_copy = .false.
  do i=1,size(structure%symmetries)
    ! q' is defined such that symmetry i maps qp onto q.
    ! N.B. q-points transform as q->R^-T.q, so q' s.t. q'->q is given by
    !    q' = R^T.q = q.R.
    q_prime = structure%symmetries(i)%inverse_transform(qpoint)
    do j=1,size(supercells)
      if (is_int(supercells(j)%supercell*q_prime%qpoint)) then
        is_copy(j,i) = .true.
      endif
    enddo
  enddo
  no_copies = count(is_copy)
  
  ! Construct a copy of the dynamical matrix from each pair of supercell and
  !    symmetry.
  allocate( copy_matrices(structure%no_atoms,structure%no_atoms,no_copies), &
          & sym_id(no_copies),                                              &
          & sup_id(no_copies),                                              &
          & stat=ialloc); call err(ialloc)
  k = 0
  do i=1,size(structure%symmetries)
    q_prime = structure%symmetries(i)%inverse_transform(qpoint)
    do j=1,size(supercells)
      if (is_copy(j,i)) then
        k = k+1
        dynamical_matrix = DynamicalMatrix( dblevec(q_prime%qpoint), &
                                          & supercells(j),           &
                                          & hessian(j)               )
        dynamical_matrix = transform_modes( dynamical_matrix,        &
                                          & structure%symmetries(i), &
                                          & q_prime,                 &
                                          & qpoint                   )
        copy_matrices(:,:,k) = dynamical_matrix%elements()
        sym_id(k) = i
        sup_id(k) = j
      endif
    enddo
  enddo
  
  ! Work out if the output matrices should be equivalent.
  ! If the calculation is coming from more than one supercell, or from a
  !    supercell which is larger than necessary to simulate the q-point,
  !    then the separate contributions to the dynamical matrix will differ.
  check_matrices = .true.
  if (count(any(is_copy,2))>1) then
    check_matrices = .false.
  endif
  do i=1,size(supercells)
    if (any(is_copy(i,:))) then
      if (supercells(i)%sc_size/=qpoint%min_sc_size()) then
        check_matrices = .false.
      endif
    endif
  enddo
  
  if (check_matrices) then
    allocate( average(structure%no_atoms,structure%no_atoms),     &
            & averages(structure%no_atoms,structure%no_atoms),    &
            & differences(structure%no_atoms,structure%no_atoms), &
            & stat=ialloc); call err(ialloc)
    average = cmplxmat(zeroes(3,3))
    averages = 0
    differences = 0
    do i=1,no_copies
      average = average + copy_matrices(:,:,i)/no_copies
      do j=1,i-1
        averages = averages &
               & + sum_squares((copy_matrices(:,:,j)+copy_matrices(:,:,i))/2)
        differences = differences &
                  & + sum_squares(copy_matrices(:,:,j)-copy_matrices(:,:,i))
      enddo
    enddo
    l2_error = sqrt(sum(differences)/sum(averages))
    call logfile%print_line('Fractional L2 difference between equivalent &
       &dynamical matrices: '//l2_error)
    if (l2_error>1e-10_dp) then
      call print_line(WARNING//': Symmetrically equivalent dynamical &
         &matrices differ. Please check log files.')
      call print_line('Fractional L2 error: '//l2_error)
    endif
  endif
  
  ! Average over the copies.
  allocate( matrices(structure%no_atoms,structure%no_atoms), &
          & stat=ialloc); call err(ialloc)
  matrices = cmplxmat(zeroes(3,3))
  do i=1,no_copies
    matrices = matrices + copy_matrices(:,:,i)/no_copies
  enddo
  
  dynamical_matrix = DynamicalMatrix(matrices)
  
  ! Diagonalise the dynamical matrix, to obtain the normal mode
  !    co-ordinates (eigenvectors) and harmonic frequencies (eigenvalues).
  complex_modes = ComplexMode( dynamical_matrix,                    &
                             & structure,                           &
                             & modes_real=qpoint%is_paired_qpoint() )
  
  ! Process the modes to calculate ids.
  complex_modes = process_modes(complex_modes, structure, qpoint, subspace_id)
  
  output = MatrixAndModes(dynamical_matrix, complex_modes)
end function

! ----------------------------------------------------------------------
! Compares two dynamical matrices.
! ----------------------------------------------------------------------
subroutine compare_dynamical_matrices(a,b,logfile)
  implicit none
  
  type(DynamicalMatrix), intent(in)    :: a
  type(DynamicalMatrix), intent(in)    :: b
  type(OFile),           intent(inout) :: logfile
  
  integer :: no_atoms
  
  type(ComplexMatrix), allocatable :: mats_a(:,:)
  type(ComplexMatrix), allocatable :: mats_b(:,:)
  
  type(ComplexMatrix) :: mat_a
  type(ComplexMatrix) :: mat_b
  
  real(dp) :: average
  real(dp) :: difference
  
  integer :: i,j
  
  mats_a = a%elements()
  mats_b = b%elements()
  
  no_atoms = size(mats_a,1)
  if (size(mats_b,1)/=no_atoms) then
    call print_line(CODE_ERROR//': dynamical matrices a and b have different &
       &sizes.')
    call err()
  endif
  
  average = 0.0_dp
  difference = 0.0_dp
  do i=1,no_atoms
    do j=1,no_atoms
      mat_a = mats_a(j,i)
      mat_b = mats_b(j,i)
      average = average + sum_squares((mat_a+mat_b)/2)
      difference = difference + sum_squares(mat_a-mat_b)
    enddo
  enddo
  
  if (average>1e-30_dp) then
    call logfile%print_line('Fractional L2 difference between dynamical &
       &matrices: '//sqrt(difference/average))
    if (sqrt(difference/average)>1e-10_dp) then
      call print_line(WARNING//': Dynamical matrices differ. Please check &
         &log files.')
      call print_line('Fractional L2 difference between dynamical &
         &matrices: '//sqrt(difference/average))
    endif
  else
    call logfile%print_line('Dynamical matrices too small to compare &
       &fractional differences.')
  endif
end subroutine

! ----------------------------------------------------------------------
! Transform a dynamical matrix and set of normal modes
!    from one q-point to another.
! ----------------------------------------------------------------------
! Construct data at q_new from data at q_old, where
!    R . q_old = q_new.
! N.B. the provided q-point should be q_new not q_old.
! The symmetry, S, maps equilibrium position ri to rj+R, and q-point q to q'.
! The +R needs to be corrected for, by multiplying the phase by exp(-iq'.R).
function transform_modes(input,symmetry,qpoint_from,qpoint_to) result(output)
  implicit none
  
  type(DynamicalMatrix),  intent(in)    :: input
  type(SymmetryOperator), intent(in)    :: symmetry
  type(QpointData),       intent(in)    :: qpoint_from
  type(QpointData),       intent(in)    :: qpoint_to
  type(DynamicalMatrix)                 :: output
  
  output = DynamicalMatrix(                                       &
     & elements = transform_dynamical_matrix( input%elements(),   &
     &                                        symmetry,           &
     &                                        qpoint_from,        &
     &                                        qpoint_to         ) )
end function

function transform_dynamical_matrix(input,symmetry,qpoint_from,qpoint_to) &
   & result(output)
  implicit none
  
  type(ComplexMatrix),    intent(in) :: input(:,:)
  type(SymmetryOperator), intent(in) :: symmetry
  type(QpointData),       intent(in) :: qpoint_from
  type(QpointData),       intent(in) :: qpoint_to
  type(ComplexMatrix), allocatable   :: output(:,:)
  
  integer :: no_atoms
  
  type(FractionVector) :: q
  type(IntVector)      :: r1
  type(IntVector)      :: r2
  
  integer :: atom_1
  integer :: atom_1p
  integer :: atom_2
  integer :: atom_2p
  
  integer :: ialloc
  
  no_atoms = size(input,1)
  q = qpoint_to%qpoint
  
  ! Check that the symmetry transforms the q-point as expected.
  if (symmetry * qpoint_from /= qpoint_to) then
    call print_line(CODE_ERROR//': Symmetry does not transform q-points as &
       &expected.')
    call err()
  endif
  
  ! Check that the number of atoms is consistent.
  if (size(input,2)/=no_atoms) then
    call err()
  endif
  
  ! Transform dynamical matrix.
  allocate(output(no_atoms,no_atoms), stat=ialloc); call err(ialloc)
  do atom_1=1,no_atoms
    atom_1p = symmetry%atom_group * atom_1
    r1 = symmetry%rvectors(atom_1)
    do atom_2=1,no_atoms
      atom_2p = symmetry%atom_group * atom_2
      r2 = symmetry%rvectors(atom_2)
      
      output(atom_2p,atom_1p) =      &
         & symmetry%cartesian_tensor &
         & * exp_2pii(-q*r2)         &
         & * input(atom_2,atom_1)    &
         & * exp_2pii(q*r1)          &
         & * transpose(symmetry%cartesian_tensor)
    enddo
  enddo
end function
end module
