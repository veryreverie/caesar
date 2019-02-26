! ======================================================================
! The forces between atoms at a given q-point.
! ======================================================================
module dynamical_matrix_module
  use utils_module
  
  use structure_module
  use normal_mode_module
  
  use min_images_module
  use cartesian_hessian_module
  implicit none
  
  private
  
  public :: DynamicalMatrix
  public :: compare_dynamical_matrices
  public :: transform_modes
  public :: conjg
  public :: reconstruct_hessian
  
  type, extends(Stringsable) :: DynamicalMatrix
    type(ComplexMatrix), allocatable, private :: matrices_(:,:)
    type(ComplexMode),   allocatable          :: complex_modes(:)
  contains
    procedure, public :: check
    
    ! I/O.
    procedure, public :: read  => read_DynamicalMatrix
    procedure, public :: write => write_DynamicalMatrix
  end type
  
  interface DynamicalMatrix
    module procedure new_DynamicalMatrix_calculated
    module procedure new_DynamicalMatrix_interpolated
    module procedure new_DynamicalMatrix_ComplexModes
    module procedure new_DynamicalMatrix_Strings
    module procedure new_DynamicalMatrix_StringArray
  end interface
  
  interface conjg
    module procedure conjg_DynamicalMatrix
  end interface
  
  interface ComplexMode
    module procedure new_ComplexMode_interpolated_complex
    module procedure new_ComplexMode_interpolated_real
    module procedure new_ComplexMode_calculated
  end interface
  
  ! Types used for splitting degenerate modes. The type SplitModes is for
  !    modes split by commuting symmetries, and the type AntiSplitModes is for
  !    modes split by anti-commuting symmetries.
  type, extends(NoDefaultConstructor) :: SplitModes
    type(ComplexMode), allocatable :: modes(:)
    integer,           allocatable :: phases(:)
  end type
  
  interface SplitModes
    module procedure new_SplitModes
    module procedure new_SplitModes_ComplexModes
  end interface
  
  type, extends(NoDefaultConstructor) :: AntiSplitModes
    type(ComplexMode), allocatable :: modes(:)
  end type
  
  interface AntiSplitModes
    module procedure new_AntiSplitModes
    module procedure new_AntiSplitModes_ComplexModes
  end interface
contains

! Constructors.
function new_SplitModes(modes,phases) result(this)
  implicit none
  
  type(ComplexMode), intent(in) :: modes(:)
  integer,           intent(in) :: phases(:)
  type(SplitModes)              :: this
  
  this%modes = modes
  this%phases = phases
end function

function new_AntiSplitModes(modes) result(this)
  implicit none
  
  type(ComplexMode), intent(in) :: modes(:)
  type(AntiSplitModes)          :: this
  
  this%modes = modes
end function

! ----------------------------------------------------------------------
! Constructs the dynamical matrix, which is the matrix of force constants in
!    q-point co-ordinates,
!    given the Hessian, which is the matrix of force constants in
!    cartesian supercell co-ordinates.
! ----------------------------------------------------------------------

! --------------------------------------------------
! Construct the dynamical matrix at the specified q-point from the given
!    Hessian.
! Considers all supercell, S, and symmetry, R, combinations such that
!    S.R.q is a vector of integers, i.e. the q-point is a G-vector of the
!    supercell.
! --------------------------------------------------
function new_DynamicalMatrix_calculated(qpoint,supercells,hessian, &
   & structure,subspace_id,logfile) result(this)
  implicit none
  
  type(QpointData),       intent(in)    :: qpoint
  type(StructureData),    intent(in)    :: supercells(:)
  type(CartesianHessian), intent(in)    :: hessian(:)
  type(StructureData),    intent(in)    :: structure
  integer,                intent(in)    :: subspace_id
  type(OFile),            intent(inout) :: logfile
  type(DynamicalMatrix)                 :: this
  
  type(QpointData) :: q_prime
  
  logical, allocatable :: is_copy(:,:)
  integer, allocatable :: sym_id(:)
  integer, allocatable :: sup_id(:)
  integer              :: no_copies
  integer              :: i,j,k,ialloc
  
  type(ComplexMatrix), allocatable :: matrix(:,:)
  type(ComplexMatrix), allocatable :: matrices(:,:,:)
  
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
  allocate( matrices(structure%no_atoms,structure%no_atoms,no_copies), &
          & sym_id(no_copies),                                         &
          & sup_id(no_copies),                                         &
          & stat=ialloc); call err(ialloc)
  k = 0
  do i=1,size(structure%symmetries)
    q_prime = structure%symmetries(i)%inverse_transform(qpoint)
    do j=1,size(supercells)
      if (is_copy(j,i)) then
        k = k+1
        matrix = calculate_dynamical_matrix( dblevec(q_prime%qpoint), &
                                           & supercells(j),           &
                                           & hessian(j)               )
        matrices(:,:,k) = transform_dynamical_matrix( &
                           & matrix,                  &
                           & structure%symmetries(i), &
                           & q_prime,                 &
                           & qpoint                   )
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
      average = average + matrices(:,:,i)/no_copies
      do j=1,i-1
        averages = averages + sum_squares((matrices(:,:,j)+matrices(:,:,i))/2)
        differences = differences &
                  & + sum_squares(matrices(:,:,j)-matrices(:,:,i))
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
  allocate( this%matrices_(structure%no_atoms,structure%no_atoms), &
          & stat=ialloc); call err(ialloc)
  this%matrices_ = cmplxmat(zeroes(3,3))
  do i=1,no_copies
    this%matrices_ = this%matrices_ + matrices(:,:,i)/no_copies
  enddo
  
  ! Diagonalise the dynamical matrix, to obtain the normal mode
  !    co-ordinates (eigenvectors) and harmonic frequencies (eigenvalues).
  this%complex_modes = ComplexMode( this%matrices_,    &
                                  & structure,         &
                                  & qpoint,            &
                                  & subspace_id        )
end function

! --------------------------------------------------
! Construct a dynamical matrix at a q-point which is not commensurate
!    with the supercell.
! This is only an approximation, using a minimum-image convention.
! --------------------------------------------------
function new_DynamicalMatrix_interpolated(q,supercell,hessian,min_images) &
   & result(this)
  implicit none
  
  type(RealVector),       intent(in) :: q
  type(StructureData),    intent(in) :: supercell
  type(CartesianHessian), intent(in) :: hessian
  type(MinImages),        intent(in) :: min_images(:,:)
  type(DynamicalMatrix)              :: this
  
  ! Evaluate the dynamical matrix.
  this%matrices_ = calculate_dynamical_matrix( q,         &
                                             & supercell, &
                                             & hessian,   &
                                             & min_images )
  
  
  ! Diagonalise the dynamical matrix, to obtain the normal mode
  !    co-ordinates (eigenvectors) and harmonic frequencies (eigenvalues).
  this%complex_modes = ComplexMode(this%matrices_, supercell)
end function

! ----------------------------------------------------------------------
! Construct the dynamical matrix itself.
! ----------------------------------------------------------------------
function calculate_dynamical_matrix(q,supercell,hessian,min_images) &
   & result(output)
  implicit none
  
  type(RealVector),       intent(in)           :: q
  type(StructureData),    intent(in)           :: supercell
  type(CartesianHessian), intent(in)           :: hessian
  type(MinImages),        intent(in), optional :: min_images(:,:)
  type(ComplexMatrix), allocatable             :: output(:,:)
  
  type(AtomData)               :: atom_1
  type(AtomData)               :: atom_2
  type(IntVector)              :: rvector
  type(IntVector), allocatable :: rvectors(:)
  
  integer :: i,j,k,ialloc
  
  ! Allocate and zero output.
  allocate( output( supercell%no_atoms_prim,  &
          &         supercell%no_atoms_prim), &
          & stat=ialloc); call err(ialloc)
  output = cmplxmat(zeroes(3,3))
  
  ! Add up contributions to the dynamical matrix from
  !    the force constants between each pair of atoms.
  do i=1,supercell%no_atoms
    atom_1 = supercell%atoms(i)
    
    ! Atom 2 will always be at R=0, so the R-vector from atom 2 to atom 1
    !    is simply that of atom 1.
    rvector = supercell%rvectors(atom_1%rvec_id())
    do j=1,supercell%no_atoms_prim
      atom_2 = supercell%atoms(j)
      
      if (present(min_images)) then
        rvectors = min_images(atom_2%id(),atom_1%id())%image_rvectors
        ! Check that all min image R-vectors actually point from atom 2 to
        !    a copy of atom 1.
        do k=1,size(rvectors)
          if (.not. is_int( supercell%recip_supercell &
                        & * (rvectors(k)-rvector))) then
            call print_line(CODE_ERROR// &
               & ': Problem with minimum image R-vectors.')
            call err()
          endif
        enddo
      else
        rvectors = [rvector]
      endif
      
      output(atom_2%prim_id(),atom_1%prim_id()) =    &
         & output(atom_2%prim_id(),atom_1%prim_id()) &
         & + hessian%elements(atom_2,atom_1)         &
         & * sum(exp_2pii(q*rvectors))               &
         & / (supercell%sc_size*size(rvectors))
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Construct the dynamical matrix and normal modes at the q-point q'=G-q,
!    given the dynamical matrix and normal modes at the q-point q.
! ----------------------------------------------------------------------
function conjg_DynamicalMatrix(input) result(output)
  implicit none
  
  type(DynamicalMatrix), intent(in) :: input
  type(DynamicalMatrix)             :: output
  
  ! The dynamical matrix at G-q is the complex conjugate of that at q.
  ! N.B. the complex conjugate, not the Hermitian conjugate.
  output%matrices_ = conjg(input%matrices_)
  
  ! The displacements of the normal modes at G-q are the complex conjugates
  !    of those at q.
  output%complex_modes = conjg(input%complex_modes)
end function

! ----------------------------------------------------------------------
! Check a dynamical matrix.
! ----------------------------------------------------------------------
! Always checks that the matrix is Hermitian.
! If check_eigenstuff is .true., also checks that the normal modes match
!    the dynamical matrix.
! check_eigenstuff defaults to .true..
! Structure may be any supercell.
subroutine check(this,structure,logfile,check_eigenstuff)
  implicit none
  
  class(DynamicalMatrix), intent(in)           :: this
  type(StructureData),    intent(in)           :: structure
  type(OFile),            intent(inout)        :: logfile
  logical,                intent(in), optional :: check_eigenstuff
  
  type(ComplexMatrix)            :: matrix
  type(ComplexMatrix)            :: hermitian_matrix
  type(ComplexMode), allocatable :: modes(:)
  real(dp)                       :: freq_1
  real(dp)                       :: freq_2
  type(ComplexVector)            :: prim_vec_1
  type(ComplexVector)            :: prim_vec_2
  real(dp)                       :: average
  real(dp)                       :: difference
  logical                        :: check_estuff
  
  integer :: no_atoms
  integer :: i,j
  
  if (present(check_eigenstuff)) then
    check_estuff = check_eigenstuff
  else
    check_estuff = .true.
  endif
  
  no_atoms = size(this%matrices_,1)
  
  ! Check the dynamical matrix is Hermitian.
  average = 0.0_dp
  difference = 0.0_dp
  do i=1,no_atoms
    do j=1,i
      matrix = this%matrices_(j,i)
      hermitian_matrix = hermitian(this%matrices_(i,j))
      
      average = average + sum_squares((matrix+hermitian_matrix)/2.0_dp)
      difference = difference + sum_squares(matrix-hermitian_matrix)
    enddo
  enddo
  call logfile%print_line(                                        &
     & 'Fractional L2 error in hermicity of dynamical matrix :'// &
     & sqrt(difference/average))
  if (sqrt(difference/average)>1.0e-10_dp) then
    call print_line(WARNING//': Dynamical matrix is not hermitian. Please &
       &check log files.')
    call print_line('Fractional L2 error in hermicity: '// &
       & sqrt(difference/average))
  endif
  
  ! Check that dynamical matrix and normal modes match.
  if (check_estuff) then
    modes = ComplexMode(this%matrices_, structure)
    
    ! Check that eigenfrequencies match.
    average = 0.0_dp
    difference = 0.0_dp
    do i=1,structure%no_modes_prim
      ! Ignore translational and degenerate modes.
      if (this%complex_modes(i)%translational_mode) then
        cycle
      elseif (count( this%complex_modes%subspace_id &
                & == this%complex_modes(i)%subspace_id)>1) then
        cycle
      endif
      freq_1 = this%complex_modes(i)%frequency
      freq_2 = modes(i)%frequency
      average = average + ((freq_1+freq_2)/2)**2
      difference = difference + (freq_1-freq_2)**2
    enddo
    call logfile%print_line(                           &
       & 'Fractional L2 error in eigenfrequencies: '// &
       & sqrt(difference/average))
    if (sqrt(difference/average) > 1e-10_dp) then
      call print_line(WARNING//': Eigenfrequencies do not match. Please &
         &check log files.')
      call print_line('Fractional L2 error: '//sqrt(difference/average))
    endif
    
    ! Check that the primitive displacements match.
    ! N.B. the global phase of the displacements is not necessarily consistent.
    average = 0.0_dp
    difference = 0.0_dp
    do i=1,structure%no_modes_prim
      ! Ignore translational and degenerate modes.
      if (this%complex_modes(i)%translational_mode) then
        cycle
      elseif (count( this%complex_modes%subspace_id &
                & == this%complex_modes(i)%subspace_id)>1) then
        cycle
      endif
      
      do j=1,structure%no_atoms_prim
        prim_vec_1 = this%complex_modes(i)%unit_vector(j)
        prim_vec_2 = modes(i)%unit_vector(j)
        ! Ignore phases.
        prim_vec_1 = cmplxvec(vec(abs(cmplx(prim_vec_1))))
        prim_vec_2 = cmplxvec(vec(abs(cmplx(prim_vec_2))))
        average = average + sum_squares((prim_vec_1+prim_vec_2)/2)
        difference = difference + sum_squares(prim_vec_1-prim_vec_2)
      enddo
    enddo
    call logfile%print_line(                                        &
       & 'Fractional L2 error in transformed primitive vectors: '// &
       & sqrt(difference/average)                                   )
    if (sqrt(difference/average) > 1e-10_dp) then
      call print_line( WARNING//': Error in primitive vectors. &
                     &Please check log files.'                 )
      call print_line('Fractional L2 error: '//sqrt(difference/average))
    endif
  endif
end subroutine

! ----------------------------------------------------------------------
! Calculates complex modes by diagonalising a dynamical matrix.
! ----------------------------------------------------------------------
! N.B. Structure may be any supercell.

! Calculate modes at an arbitrary q-point.
function new_ComplexMode_interpolated_complex(matrices,structure) &
   & result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: matrices(:,:)
  type(StructureData), intent(in) :: structure
  type(ComplexMode), allocatable  :: output(:)
  
  complex(dp),               allocatable :: dyn_mat(:,:)
  type(HermitianEigenstuff), allocatable :: estuff(:)
  
  integer :: i,j,ialloc
  
  ! Convert (3x3Matrix) x no_atoms x no_atoms to no_modes x no_modes
  allocate( dyn_mat(structure%no_modes_prim,structure%no_modes_prim), &
          & stat=ialloc); call err(ialloc)
  do i=1,structure%no_atoms_prim
    do j=1,structure%no_atoms_prim
      dyn_mat(3*j-2:3*j, 3*i-2:3*i) = cmplx(matrices(j,i))
    enddo
  enddo
  
  ! Diagonalise dynamical matrix.
  estuff = diagonalise_hermitian(dyn_mat)
  
  ! Calculate normal modes.
  ! Eigenvalues are reversed so that the complex modes are in ascending order.
  output = ComplexMode(estuff(size(estuff):1:-1), structure)
end function

! Calculate modes at an arbitrary q-point.
function new_ComplexMode_interpolated_real(matrices,structure) &
   & result(output)
  implicit none
  
  type(RealMatrix),    intent(in) :: matrices(:,:)
  type(StructureData), intent(in) :: structure
  type(ComplexMode), allocatable  :: output(:)
  
  real(dp),                  allocatable :: dyn_mat(:,:)
  type(SymmetricEigenstuff), allocatable :: real_estuff(:)
  type(HermitianEigenstuff), allocatable :: estuff(:)
  
  integer :: i,j,ialloc
  
  ! Convert (3x3Matrix) x no_atoms x no_atoms to no_modes x no_modes
  allocate( dyn_mat(structure%no_modes_prim,structure%no_modes_prim), &
          & stat=ialloc); call err(ialloc)
  do i=1,structure%no_atoms_prim
    do j=1,structure%no_atoms_prim
      dyn_mat(3*j-2:3*j, 3*i-2:3*i) = dble(matrices(j,i))
    enddo
  enddo
  
  ! Diagonalise dynamical matrix.
  real_estuff = diagonalise_symmetric(dyn_mat)
  estuff = [(                                                         &
     & HermitianEigenstuff( real_estuff(i)%eval,                      &
     &                      [cmplx(real_estuff(i)%evec,0.0_dp,dp)] ), &
     & i=1,                                                           &
     & size(real_estuff)                                              )]
  
  ! Calculate normal modes.
  ! Eigenvalues are reversed so that the complex modes are in ascending order.
  output = ComplexMode(estuff(size(estuff):1:-1), structure)
end function

! Calculate modes for one of the calculated q-points.
! Will choose correct basis in degenerate spaces using symmetry.
function new_ComplexMode_calculated(matrices,structure,qpoint, &
   & subspace_id) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in)    :: matrices(:,:)
  type(StructureData), intent(in)    :: structure
  type(QpointData),    intent(in)    :: qpoint
  integer,             intent(in)    :: subspace_id
  type(ComplexMode), allocatable     :: output(:)
  
  integer, allocatable :: subspace_ids(:)
  integer, allocatable :: subspace_id_set(:)
  
  real(dp) :: energy_difference
  
  ! Symmetry data.
  type(SymmetryOperator), allocatable :: symmetries(:)
  
  integer, allocatable :: states(:)
  
  complex(dp) :: symmetry(2,2)
  
  complex(dp), allocatable :: overlap(:,:)
  
  integer :: i,j,k,ialloc
  
  ! Calculate normal modes as if at an arbitrary q-point.
  if (qpoint%is_paired_qpoint()) then
    output = ComplexMode(real(matrices), structure)
  else
    output = ComplexMode(matrices, structure)
  endif
  
  ! Set q-point ids and mode ids.
  output%qpoint_id = qpoint%id
  output%paired_qpoint_id = qpoint%paired_qpoint_id
  
  do i=1,size(output)
    output(i)%id        = (qpoint%id-1)*structure%no_modes_prim + i
    output(i)%paired_id = (qpoint%paired_qpoint_id-1)*structure%no_modes_prim &
                      & + i
  enddo
  
  ! Identify purely translational modes (at the gamma-point only).
  output%translational_mode = .false.
  if (qpoint%is_gvector()) then
    ! Find the three modes with the minimum abs(frequency).
    do i=1,3
      j = minloc( abs(output%frequency),              &
                & dim=1,                              &
                & mask=.not.output%translational_mode )
      output(j)%translational_mode = .true.
    enddo
  endif
  
  ! Identify the symmetries which map the q-point to itself.
  symmetries = structure%symmetries(               &
     & filter(structure%symmetries*qpoint==qpoint) )
  
  ! Assign subspace ids, which are equal if two states are degenerate,
  !    and different if they are not.
  subspace_ids = [(i,i=1,size(output))]
  do i=1,size(output)
    do j=1,size(output)
      if (subspace_ids(j)==subspace_ids(i)) then
        continue
      elseif ( qpoint%is_paired_qpoint() .and.                            &
             & abs(sum(output(i)%unit_vector*output(j)%unit_vector))>1e-4 &
             & ) then
        ! If conjg(u1).u2 is non-zero then u1 and u2 are degenerate.
        ! This is only possible at 2q=G.
        subspace_ids(filter(subspace_ids==subspace_ids(j))) = subspace_ids(i)
      else
        ! If u1.T.u2 is non-zero then u1 and u2 are degenerate.
        do k=1,size(symmetries)
          symmetry = cmplx(calculate_symmetry_in_normal_coordinates(   &
                                           & [output(i), output(j)],   &
                                           & qpoint,                   &
                                           & symmetries(k)           ) )
          if (abs(symmetry(1,2))>1e-4 .or. abs(symmetry(2,1))>1e-4) then
            subspace_ids(filter(subspace_ids==subspace_ids(j))) = &
               & subspace_ids(i)
            exit
          endif
        enddo
      endif
    enddo
  enddo
  
  subspace_id_set = subspace_ids(set(subspace_ids))
  i = subspace_id
  do j=1,size(subspace_id_set)
    output(filter(subspace_ids==subspace_id_set(j)))%subspace_id = i
    i = i+1
  enddo
  
  ! Loop over subspaces, choosing the correct basis using symmetry operators.
  ! The correct basis has <p|H|q>=0 for all H which are invariant under
  !    the symmetries of the system.
  do i=subspace_id,output(size(output))%subspace_id
    ! Find the set of states with degeneracy id i.
    states = filter(output%subspace_id==i)
    
    ! Set the frequencies of each degenerate state to
    !    the average of their frequencies.
    output(states)%frequency = sum(output(states)%frequency) / size(states)
    
    ! Choose basis using symmetry operators.
    if (size(states)>1) then
      if (qpoint%is_paired_qpoint()) then
        output(states) = choose_basis_real( output(states), &
                                          & structure,      &
                                          & symmetries,     &
                                          & qpoint          )
      else
        output(states) = choose_basis_complex( output(states), &
                                             & structure,      &
                                             & symmetries,     &
                                             & qpoint          )
      endif
      
      ! Check that the chosen basis is orthonormal.
      allocate( overlap(size(states),size(states)), &
              & stat=ialloc); call err(ialloc)
      do j=1,size(states)
        do k=1,size(states)
          overlap(k,j) = sum( conjg(output(states(k))%unit_vector) &
                          & * output(states(j))%unit_vector        )
        enddo
      enddo
      call check_identity(abs(mat(overlap)), 'Overlap matrix')
      deallocate(overlap, stat=ialloc); call err(ialloc)
    endif
  enddo
end function

! --------------------------------------------------
! Chooses the basis for the degenerate subspace using symmetry operators.
! --------------------------------------------------
! Symmetries must all take the q-point to itself.
function choose_basis_complex(input,structure,symmetries,qpoint) &
   & result(output)
  implicit none
  
  type(ComplexMode),      intent(in)    :: input(:)
  type(StructureData),    intent(in)    :: structure
  type(SymmetryOperator), intent(in)    :: symmetries(:)
  type(QpointData),       intent(in)    :: qpoint
  type(ComplexMode), allocatable        :: output(:)
  
  integer, allocatable :: ids(:)
  integer, allocatable :: id_set(:)
  integer, allocatable :: phases(:)
  integer, allocatable :: phases_set(:)
  integer              :: max_id
  type(SplitModes)     :: split_modes
  logical              :: symmetry_used
  
  type(SymmetryOperator), allocatable :: used_symmetries(:)
  type(IntArray1D),       allocatable :: used_phases(:)
  type(AntiSplitModes)                :: anti_split_modes
  
  type(SymmetryOperator), allocatable :: anticommuting_symmetries(:)
  
  integer,           allocatable :: new_ids(:)
  integer,           allocatable :: new_id_set(:)
  logical,           allocatable :: successes(:)
  type(ComplexMode), allocatable :: new_output(:)
  
  integer :: i,j,k,l,m,ialloc
  
  if (size(input)==1) then
    call print_line(CODE_ERROR//': Trying to lift the degeneracy of only one &
       &state.')
    call err()
  endif
  
  output = input
  used_symmetries = [SymmetryOperator::]
  used_phases = [IntArray1D::]
  ids = [(0,i=1,size(output))]
  phases = [(0,i=1,size(output))]
  max_id = 0
  
  ! Find a basis which diagonalises a maximal set of commuting symmetries.
  ! If for every pair of modes u1 and u2 there is a symmetry T for which
  !    u1 and u2 are eigenvectors with different eigenvalues, then
  !    <u1|H|u2>=0 for all u1 and u2, and so the basis is well-defined.
  do i=1,size(symmetries)
    if (all(operators_commute(symmetries(i),used_symmetries,qpoint))) then
      symmetry_used = .false.
      
      ! Loop over positive and negative superpositions of the symmetry.
      do j=1,2
        ! Loop over each distinct id.
        id_set = ids(set(ids))
        do k=1,size(id_set)
          if (count(ids==id_set(k))>1) then
            split_modes = SplitModes( output(filter(ids==id_set(k))), &
                                    & symmetries(i),                  &
                                    & qpoint,                         &
                                    & positive_superposition = j==1   )
            phases(filter(ids==id_set(k))) = split_modes%phases
            phases_set = split_modes%phases(set(split_modes%phases))
            if (size(phases_set)>1) then
              output(filter(ids==id_set(k))) = split_modes%modes
              do l=1,size(phases_set)
                max_id = max_id + 1
                ids(filter(ids==id_set(k).and.phases==phases_set(l))) = max_id
              enddo
              symmetry_used = .true.
            endif
          endif
        enddo
      enddo
      
      if (symmetry_used) then
        used_symmetries = [used_symmetries, symmetries(i)]
        used_phases = [used_phases, IntArray1D(phases)]
      endif
   endif
  enddo
  
  if (size(set(ids))/=size(ids)) then
    call print_line(ERROR//': Unable to lift degeneracies using symmetry.')
    call print_line('q-point '//qpoint%id//': '//qpoint%qpoint)
    stop 1
  endif
end function

function choose_basis_real(input,structure,symmetries,qpoint) &
   & result(output)
  implicit none
  
  type(ComplexMode),      intent(in)    :: input(:)
  type(StructureData),    intent(in)    :: structure
  type(SymmetryOperator), intent(in)    :: symmetries(:)
  type(QpointData),       intent(in)    :: qpoint
  type(ComplexMode), allocatable        :: output(:)
  
  integer, allocatable :: ids(:)
  integer, allocatable :: id_set(:)
  integer, allocatable :: phases(:)
  integer, allocatable :: phases_set(:)
  integer              :: max_id
  type(SplitModes)     :: split_modes
  logical              :: symmetry_used
  
  type(SymmetryOperator), allocatable :: used_symmetries(:)
  type(IntArray1D),       allocatable :: used_phases(:)
  type(AntiSplitModes)                :: anti_split_modes
  
  type(SymmetryOperator), allocatable :: anticommuting_symmetries(:)
  
  type(SymmetryOperator), allocatable :: antisymmetric_symmetries(:)
  
  integer,           allocatable :: new_ids(:)
  integer,           allocatable :: new_id_set(:)
  logical,           allocatable :: successes(:)
  type(ComplexMode), allocatable :: new_output(:)
  
  integer :: i,j,k,l,m,ialloc
  
  if (size(input)==1) then
    call print_line(CODE_ERROR//': Trying to lift the degeneracy of only one &
       &state.')
    call err()
  endif
  
  output = input
  used_symmetries = [SymmetryOperator::]
  used_phases = [IntArray1D::]
  ids = [(0,i=1,size(output))]
  phases = [(0,i=1,size(output))]
  max_id = 0
  
  ! Attempt to find a basis using only commuting symmetric symmetries.
  ! If for every pair of modes u1 and u2 there is a symmetry T for which
  !    u1 and u2 are eigenvectors with different eigenvalues, then
  !    <u1|H|u2>=0 for all u1 and u2, and so the basis is well-defined.
  do i=1,size(symmetries)
    if (all(superposed_operators_commute( symmetries(i),   &
                                        & used_symmetries, &
                                        & qpoint           ))) then
      symmetry_used = .false.
      
      ! Loop over each distinct id.
      id_set = ids(set(ids))
      do j=1,size(id_set)
        if (count(ids==id_set(j))>1) then
          split_modes = SplitModes( output(filter(ids==id_set(j))), &
                                  & symmetries(i),                  &
                                  & qpoint,                         &
                                  & positive_superposition = .true. )
          phases(filter(ids==id_set(j))) = split_modes%phases
          phases_set = split_modes%phases(set(split_modes%phases))
          if (size(phases_set)>1) then
            output(filter(ids==id_set(j))) = split_modes%modes
            do k=1,size(phases_set)
              max_id = max_id + 1
              ids(filter(ids==id_set(j).and.phases==phases_set(k))) = max_id
            enddo
            symmetry_used = .true.
          endif
        endif
      enddo
      
      if (symmetry_used) then
        used_symmetries = [used_symmetries, symmetries(i)]
        used_phases = [used_phases, IntArray1D(phases)]
      endif
    endif
  enddo
  
  ! If the basis can't be fully defined using symmetric symmetries,
  !    lift the remaining ambiguity using anti-symmetric symmetries.
  ! If for every pair of modes u1 and u2 there is a symmetry T for which
  !    T.u1=u2 and T.u2=-u1 then <u1|H|u2>=-(<u1|H|u2>)*. Since the states
  !    and Hamiltonian are real at 2q=G, <u1|H|u2>=0.
  id_set = ids(set(ids))
  if (size(id_set)/=size(ids)) then
    ! List the symmetries which commute with the used symmetries.
    antisymmetric_symmetries = symmetries(filter([(                      &
       & all(operators_commute(symmetries(i), used_symmetries, qpoint)), &
       & i=1,                                                            &
       & size(symmetries)                                                )]))
    
    ! Remove the symmetries with order 2 or less,
    !    since for these symmetries S^T+S, so S-S^T=0.
    antisymmetric_symmetries = antisymmetric_symmetries(           &
       & filter(antisymmetric_symmetries%symmetry_order(qpoint)>2) )
    
    ! Only take one from each commuting set.
    antisymmetric_symmetries = antisymmetric_symmetries(filter([(         &
        & .not.any(operators_commute( antisymmetric_symmetries(:i-1),     &
        &                             antisymmetric_symmetries(i),        &
        &                             qpoint                          )), &
        & i=1,                                                            &
        & size(antisymmetric_symmetries)                                  )]))
    
    ! Use the anti-symmetric symmetries to choose the correct basis,
    !    id by id.
    do i=1,size(id_set)
      if (count(ids==id_set(i))>1) then
        anti_split_modes = AntiSplitModes( output(filter(ids==id_set(i))), &
                                         & antisymmetric_symmetries,       &
                                         & qpoint                          )
        output(filter(ids==id_set(i))) = anti_split_modes%modes
      endif
    enddo
  endif
  
  ! If 2q=G then every mode is real, and is its own pair under inversion.
  output%paired_id = output%id
end function

function new_SplitModes_ComplexModes(input,symmetry,qpoint, &
   & positive_superposition) result(this)
  implicit none
  
  type(ComplexMode),      intent(in)    :: input(:)
  type(SymmetryOperator), intent(in)    :: symmetry
  type(QpointData),       intent(in)    :: qpoint
  logical,                intent(in)    :: positive_superposition
  type(SplitModes)                      :: this
  
  integer                                :: order
  type(ComplexMatrix)                    :: symmetry_matrix
  type(SymmetricEigenstuff), allocatable :: real_estuff(:)
  type(HermitianEigenstuff), allocatable :: estuff(:)
  real(dp),                  allocatable :: phases_real(:)
  integer,                   allocatable :: phases_int(:)
  type(ComplexMode),         allocatable :: modes(:)
  
  integer :: i,j,k,ialloc
  
  if (qpoint%is_paired_qpoint() .and. .not. positive_superposition) then
    call err()
  endif
  
  ! Calculate the order of the symmetry, n s.t. U^n=I.
  order = symmetry%symmetry_order(qpoint)
  
  ! Construct the symmetry, U, in normal mode co-ordinates.
  symmetry_matrix = calculate_symmetry_in_normal_coordinates( input,   &
                                                            & qpoint,  &
                                                            & symmetry )
  
  ! Instead of directly calculating the eigenstuff of the unitary symmetry
  !    matrices {U}, it is more stable to calculate the eigenstuff of the
  !    Hermitian matrices {C=(U+U^T)/2} and {S=(U-U^T)/2i}.
  ! The eigenvalues of U are e^(2*pi*i*j/n), so the eigenvalues of C and S are
  !    cos(2*pi*j/n) and sin(2*pi*j/n) respectively.
  if (positive_superposition) then
    symmetry_matrix = (symmetry_matrix + hermitian(symmetry_matrix)) &
                  & / 2.0_dp
  else
    symmetry_matrix = (symmetry_matrix - hermitian(symmetry_matrix)) &
                  & / cmplx(0.0_dp,2.0_dp,dp)
  endif
  
  ! Diagonalise the Hermitian symmetry,
  !    and convert the eigenvalues into phases.
  if (qpoint%is_paired_qpoint()) then
    real_estuff = diagonalise_symmetric(real(symmetry_matrix))
    estuff = [(                                                            &
       & HermitianEigenstuff( eval=real_estuff(i)%eval,                    &
       &                      evec=cmplx(real_estuff(i)%evec,0.0_dp,dp) ), &
       & i=1,                                                              &
       & size(real_estuff)                                                 )]
  else
    estuff = diagonalise_hermitian(symmetry_matrix)
  endif
  
  ! Correct for numerical errors taking the eigenvalue outside the range
  !    [-1,1].
  if (any(abs(estuff%eval)>1.01_dp)) then
    call print_line(ERROR//': Symmetry eigenvalue outside the range [-1,1].')
    call err()
  endif
  estuff%eval = max(-1.0_dp, min(estuff%eval, 1.0_dp))
  
  ! Convert eigenvalues into phases.
  ! If the eigenvalue is cos(2 pi j/order) then the phase is j.
  ! N.B. because sin(x)=a has two solutions, 
  if (positive_superposition) then
    phases_real = acos(estuff%eval)*order/(2.0_dp*PI)
  else
    phases_real = asin(estuff%eval)*order/(2.0_dp*PI)
  endif
  allocate(phases_int(size(phases_real)), stat=ialloc); call err(ialloc)
  phases_int = nint(phases_real)
  
  ! sin(x)=sin(pi-x). Normally this will not cause problems, since the
  !    actual value of x is unimportant, and it only matters that different
  !    eigenvalues are distinguishable (which they are).
  ! If order is odd, then only one of x and pi-x will give x=2pij/order
  !    with j as an integer.
  if (modulo(order,2)==1 .and. .not. positive_superposition) then
    do i=1,size(phases_real)
      if (abs(0.5_dp-abs(phases_int(i)-phases_real(i)))<0.1_dp) then
        phases_real(i) = (order/2.0_dp)-phases_real(i)
        phases_int(i) = nint(phases_real(i))
      endif
    enddo
  endif
  
  if (any(abs(phases_int-phases_real)>0.1_dp)) then
    call print_line(ERROR//': Symmetry with non-integer phase eigenvalue.')
    call err()
  endif
  phases_int = modulo(phases_int, order)
  
  ! If the symmetry lifts degeneracy (has multiple phases), then transform the
  !    input vectors into the symmetry's eigenbasis.
  modes = input
  if (any(phases_int/=phases_int(1))) then
    do i=1,size(input)
      do j=1,size(modes(i)%unit_vector)
        modes(i)%unit_vector(j) = cmplxvec(zeroes(3))
        do k=1,size(estuff(i)%evec)
          modes(i)%unit_vector(j) = modes(i)%unit_vector(j) &
                                & + estuff(i)%evec(k)*input(k)%unit_vector(j)
        enddo
      enddo
    enddo
  endif
  
  this = SplitModes(modes, phases_int)
end function

function new_AntiSplitModes_ComplexModes(input,symmetries,qpoint) result(this)
  implicit none
  
  type(ComplexMode),      intent(in) :: input(:)
  type(SymmetryOperator), intent(in) :: symmetries(:)
  type(QpointData),       intent(in) :: qpoint
  type(AntiSplitModes)               :: this
  
  real(dp), allocatable :: symmetry_matrix(:,:)
  
  type(ComplexMode), allocatable :: modes(:)
  
  integer :: i,j
  
  if (size(symmetries)/=size(input)-1) then
    call print_line(CODE_ERROR//': The number of symmetries is unexpected.')
    call err()
  endif
  
  modes = input
  ! The first mode is the first input mode.
  ! The other modes are each a symmetry acting on the first mode.
  do i=1,size(symmetries)
    symmetry_matrix = dble(real(calculate_symmetry_in_normal_coordinates( &
                                                          & input,        &
                                                          & qpoint,       &
                                                          & symmetries(i) )))
    symmetry_matrix = (symmetry_matrix - transpose(symmetry_matrix))/2
    
    modes(i+1)%unit_vector = cmplxvec(zeroes(3))
    do j=1,size(input)
      modes(i+1)%unit_vector = modes(i+1)%unit_vector &
                           & + symmetry_matrix(1,j)   &
                           & * input(j)%unit_vector
    enddo
    modes(i+1)%unit_vector = modes(i+1)%unit_vector &
                         & / l2_norm(symmetry_matrix(1,:))
  enddo
  
  this = AntiSplitModes(modes)
end function

! ----------------------------------------------------------------------
! Construct a DynamicalMatrix from the ComplexModes at a given q-point.
! ----------------------------------------------------------------------
function new_DynamicalMatrix_ComplexModes(modes,frequencies) result(this)
  implicit none
  
  type(ComplexMode), intent(in)           :: modes(:)
  real(dp),          intent(in), optional :: frequencies(:)
  type(DynamicalMatrix)                   :: this
  
  integer :: no_atoms
  
  integer :: i,j,k,ialloc
  
  if (size(modes)==0) then
    no_atoms = 1
  else
    no_atoms = size(modes(1)%unit_vector)
  endif
  
  if (size(modes)/=size(frequencies)) then
    call print_line(CODE_ERROR//': modes and frequencies do not match.')
    call err()
  elseif (size(modes)/=3*no_atoms .and. size(modes)/=3*(no_atoms-1)) then
    call print_line(CODE_ERROR//': modes and no_atoms do not match.')
  endif
  
  do i=1,size(modes)
    if (size(modes(i)%unit_vector)/=no_atoms) then
      call print_line(CODE_ERROR//': inconsistent no_atoms.')
      call err()
    endif
  enddo
  
  this%complex_modes = modes
  if (present(frequencies)) then
    do i=1,size(frequencies)
      this%complex_modes(i)%frequency = frequencies(i)
      if (frequencies(i)>=0) then
        this%complex_modes(i)%spring_constant = frequencies(i)**2
      else
        this%complex_modes(i)%spring_constant = -frequencies(i)**2
      endif
    enddo
  endif
  
  ! A dynamical matrix is D = sum_i e_i u_i^(u_i)*.
  ! The eigenvalue e_i is the negative of the mode's spring constant.
  allocate(this%matrices_(no_atoms,no_atoms), stat=ialloc); call err(ialloc)
  this%matrices_ = cmplxmat(zeroes(3,3))
  do i=1,size(modes)
    do j=1,no_atoms
      do k=1,no_atoms
        this%matrices_(k,j) = &
           &   this%matrices_(k,j)                                        &
           & - this%complex_modes(i)%spring_constant                      &
           & * outer_product( this%complex_modes(i)%unit_vector(k),       &
           &                  conjg(this%complex_modes(i)%unit_vector(j)) )
      enddo
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Construct the Hessian for the large supercell, given the
!    dynamical matrices at each q-point.
! ----------------------------------------------------------------------
function reconstruct_hessian(large_supercell,qpoints,dynamical_matrices, &
   & logfile) result(output)
  implicit none
  
  type(StructureData),   intent(in)    :: large_supercell
  type(QpointData),      intent(in)    :: qpoints(:)
  type(DynamicalMatrix), intent(in)    :: dynamical_matrices(:)
  type(OFile),           intent(inout) :: logfile
  type(CartesianHessian)               :: output
  
  type(RealMatrix), allocatable :: hessian(:,:)
  
  type(AtomData) :: atom_1
  type(AtomData) :: atom_2
  
  type(IntVector)  :: r
  type(RealVector) :: q
  
  integer :: i,j,k,ialloc
  
  allocate( hessian( large_supercell%no_atoms,  &
          &          large_supercell%no_atoms), &
          & stat=ialloc); call err(ialloc)
  hessian = dblemat(zeroes(3,3))
  
  ! Loop across q-points, summing up the contribution from
  !    the dynamical matrix at each.
  do i=1,large_supercell%sc_size
    q = dblevec(qpoints(i)%qpoint)
    do j=1,large_supercell%no_atoms
      atom_1 = large_supercell%atoms(j)
      do k=1,large_supercell%no_atoms
        atom_2 = large_supercell%atoms(k)
        
        r = large_supercell%rvectors(atom_2%rvec_id()) &
        & - large_supercell%rvectors(atom_1%rvec_id())
        
        ! Add in the contribution to the Hessian.
        hessian(atom_1%id(),atom_2%id()) =                               &
           &   hessian(atom_1%id(),atom_2%id())                          &
           & + real( dynamical_matrices(i)%matrices_( atom_1%prim_id(),  &
           &                                          atom_2%prim_id() ) &
           &       * exp_2pii(-q*r)                                      )
      enddo
    enddo
  enddo
  
  output = CartesianHessian(large_supercell, hessian, logfile)
end function

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
  
  ! Transform dynamical matrix.
  output%matrices_ = transform_dynamical_matrix( input%matrices_, &
                                               & symmetry,        &
                                               & qpoint_from,     &
                                               & qpoint_to)
  
  ! Transform normal modes.
  output%complex_modes = transform( input%complex_modes, &
                                  & symmetry,            &
                                  & qpoint_from,         &
                                  & qpoint_to)
end function

! ----------------------------------------------------------------------
! Transform a dynamical matrix from one q-point to another.
! ----------------------------------------------------------------------
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

! ----------------------------------------------------------------------
! Compares two dynamical matrices.
! ----------------------------------------------------------------------
subroutine compare_dynamical_matrices(a,b,logfile)
  implicit none
  
  type(DynamicalMatrix), intent(in)    :: a
  type(DynamicalMatrix), intent(in)    :: b
  type(OFile),           intent(inout) :: logfile
  
  integer :: no_atoms
  
  type(ComplexMatrix) :: mat_a
  type(ComplexMatrix) :: mat_b
  
  real(dp) :: average
  real(dp) :: difference
  
  integer :: i,j
  
  no_atoms = size(a%matrices_,1)
  if (size(b%matrices_,1)/=no_atoms) then
    call print_line(CODE_ERROR//': dynamical matrices a and b have different &
       &sizes.')
    call err()
  endif
  
  average = 0.0_dp
  difference = 0.0_dp
  do i=1,no_atoms
    do j=1,no_atoms
      mat_a = a%matrices_(j,i)
      mat_b = b%matrices_(j,i)
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
! I/O.
! ----------------------------------------------------------------------
subroutine read_DynamicalMatrix(this,input)
  implicit none
  
  class(DynamicalMatrix), intent(out) :: this
  type(String),           intent(in)  :: input(:)
  
  type(StringArray), allocatable :: elements(:)
  integer                        :: no_atoms
  
  integer :: i,j,k,ialloc
  
  select type(this); type is(DynamicalMatrix)
    elements = split_into_sections(input)
    no_atoms = int_sqrt(size(elements))
    allocate(this%matrices_(no_atoms,no_atoms), stat=ialloc); call err(ialloc)
    k = 0
    do i=1,no_atoms
      do j=1,no_atoms
        k = k+1
        this%matrices_(j,i) = ComplexMatrix(elements(k)%strings(2:4))
      enddo
    enddo
  class default
    call err()
  end select
end subroutine

function write_DynamicalMatrix(this) result(output)
  implicit none
  
  class(DynamicalMatrix), intent(in) :: this
  type(String), allocatable          :: output(:)
  
  integer :: no_atoms
  
  type(String) :: matrix_strings(3)
  
  integer :: i,j,k,ialloc
  
  select type(this); type is(DynamicalMatrix)
    no_atoms = size(this%matrices_,1)
    if (size(this%matrices_,2)/=no_atoms) then
      call err()
    endif
    
    allocate(output(5*no_atoms*no_atoms), stat=ialloc); call err(ialloc)
    k = 0
    do i=1,no_atoms
      do j=1,no_atoms
        k = k+1
        matrix_strings = str(this%matrices_(j,i))
        output(5*k-4) = 'Atoms: ('//j//' '//i//')'
        output(5*k-3) = matrix_strings(1)
        output(5*k-2) = matrix_strings(2)
        output(5*k-1) = matrix_strings(3)
        output(5*k)   = ''
      enddo
    enddo
  class default
    call err()
  end select
end function

function new_DynamicalMatrix_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(DynamicalMatrix)    :: this
  
  call this%read(input)
end function

impure elemental function new_DynamicalMatrix_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(DynamicalMatrix)         :: this
  
  this = DynamicalMatrix(str(input))
end function
end module
