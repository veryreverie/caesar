! ======================================================================
! The forces between atoms at a given q-point.
! ======================================================================
module dynamical_matrix_module
  use common_module
  
  use min_images_module
  use force_constants_module
  use calculate_modes_module
  implicit none
  
  private
  
  public :: DynamicalMatrix
  public :: compare_dynamical_matrices
  public :: rotate_modes
  public :: conjg
  public :: reconstruct_force_constants
  
  type, extends(Stringsable) :: DynamicalMatrix
    type(ComplexMatrix), allocatable, private :: matrices_(:,:)
    type(ComplexMode),   allocatable :: complex_modes(:)
  contains
    procedure, public :: check
    procedure, public :: frequencies
    
    procedure, public :: read  => read_DynamicalMatrix
    procedure, public :: write => write_DynamicalMatrix
  end type
  
  interface DynamicalMatrix
    module procedure new_DynamicalMatrix_calculated
    module procedure new_DynamicalMatrix_interpolated
    module procedure new_DynamicalMatrix_StringArray
  end interface
  
  interface conjg
    module procedure conjg_DynamicalMatrix
  end interface
  
  interface print_dyn_mat
    module procedure print_dyn_mat_matrices
    module procedure print_dyn_mat_DynamicalMatrix
  end interface
contains

! ----------------------------------------------------------------------
! Constructs the matrix of force constants in q-point co-ordinates,
!    given the matrix of force constants in supercell co-ordinates.
! ----------------------------------------------------------------------

! --------------------------------------------------
! Construct the dynamical matrix at the specified q-point from the given
!    force constants.
! Considers all supercell, S, and symmetry, R, combinations such that
!    S.R.q is a vector of integers, i.e. the q-point is a G-vector of the
!    supercell.
! --------------------------------------------------
function new_DynamicalMatrix_calculated(qpoint,supercells,force_constants, &
   & structure,degenerate_energy,subspace_id,logfile) result(this)
  implicit none
  
  type(QpointData),     intent(in)    :: qpoint
  type(StructureData),  intent(in)    :: supercells(:)
  type(ForceConstants), intent(in)    :: force_constants(:)
  type(StructureData),  intent(in)    :: structure
  real(dp),             intent(in)    :: degenerate_energy
  integer,              intent(in)    :: subspace_id
  type(OFile),          intent(inout) :: logfile
  type(DynamicalMatrix)               :: this
  
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
  
  ! Check that the supercells and force constants correspond to one another.
  if (size(supercells)/=size(force_constants)) then
    call print_line(CODE_ERROR//': Force constants and supercells do not &
       &match.')
    call err()
  endif
  
  ! Identify how many parings of supercell and rotation
  !    allow q to be simulated.
  allocate( is_copy(size(supercells),size(structure%symmetries)), &
          & stat=ialloc); call err(ialloc)
  is_copy = .false.
  do i=1,size(structure%symmetries)
    ! q' is defined such that rotation i maps qp onto q.
    ! N.B. q-points rotate as q->R^-T.q, so q' s.t. q'->q is given by
    !    q' = R^T.q = q.R.
    q_prime = structure%inverse_symmetries(i) * qpoint
    do j=1,size(supercells)
      if (is_int(supercells(j)%supercell*q_prime%qpoint)) then
        is_copy(j,i) = .true.
      endif
    enddo
  enddo
  no_copies = count(is_copy)
  
  ! Construct a copy of the dynamical matrix from each pair of supercell and
  !    rotation.
  allocate( matrices(structure%no_atoms,structure%no_atoms,no_copies), &
          & sym_id(no_copies), &
          & sup_id(no_copies), &
          & stat=ialloc); call err(ialloc)
  k = 0
  do i=1,size(structure%symmetries)
    q_prime = structure%inverse_symmetries(i) * qpoint
    do j=1,size(supercells)
      if (is_copy(j,i)) then
        k = k+1
        matrix = calculate_dynamical_matrix( dblevec(q_prime%qpoint), &
                                           & supercells(j),           &
                                           & force_constants(j))
        matrices(:,:,k) = rotate_dynamical_matrix( matrix,                  &
                                                 & structure%symmetries(i), &
                                                 & q_prime,                 &
                                                 & qpoint)
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
        differences = differences + sum_squares(matrices(:,:,j)-matrices(:,:,i))
      enddo
    enddo
    l2_error = sqrt(sum(differences)/sum(averages))
    call logfile%print_line('Fractional L2 difference between equivalent &
       &dynamical matrices: '//l2_error)
    if (l2_error>1e-10_dp) then
      call print_line(WARNING//': Symmetrically equivalent dynamical &
         &matrices differ. Please check log files.')
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
  this%complex_modes = calculate_modes( this%matrices_,    &
                                      & structure,         &
                                      & qpoint,            &
                                      & degenerate_energy, &
                                      & subspace_id,       &
                                      & logfile)
end function

! --------------------------------------------------
! Construct a dynamical matrix at a q-point which is not commensurate
!    with the supercell.
! This is only an approximation, using a minimum-image convention.
! --------------------------------------------------
function new_DynamicalMatrix_interpolated(q,supercell,force_constants, &
   & min_images) result(this)
  implicit none
  
  type(RealVector),     intent(in) :: q
  type(StructureData),  intent(in) :: supercell
  type(ForceConstants), intent(in) :: force_constants
  type(MinImages),      intent(in) :: min_images(:,:)
  type(DynamicalMatrix)            :: this
  
  ! Evaluate the dynamical matrix.
  this%matrices_ = calculate_dynamical_matrix( q,               &
                                             & supercell,       &
                                             & force_constants, &
                                             & min_images)
  
  
  ! Diagonalise the dynamical matrix, to obtain the normal mode
  !    co-ordinates (eigenvectors) and harmonic frequencies (eigenvalues).
  this%complex_modes = calculate_modes(this%matrices_, supercell)
end function

! ----------------------------------------------------------------------
! Construct the dynamical matrix itself.
! ----------------------------------------------------------------------
function calculate_dynamical_matrix(q,supercell,force_constants,min_images) &
   & result(output)
  implicit none
  
  type(RealVector),     intent(in)           :: q
  type(StructureData),  intent(in)           :: supercell
  type(ForceConstants), intent(in)           :: force_constants
  type(MinImages),      intent(in), optional :: min_images(:,:)
  type(ComplexMatrix), allocatable           :: output(:,:)
  
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
      
      output(atom_2%prim_id(),atom_1%prim_id()) =      &
         & output(atom_2%prim_id(),atom_1%prim_id())   &
         & + force_constants%constants(atom_2,atom_1)  &
         & * sum(exp_2pii(-q*rvectors))                &
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
! If check_eigenstuff, then degenerate_energy must be present.
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
  endif
  
  ! Check that dynamical matrix and normal modes match.
  if (check_estuff) then
    modes = calculate_modes(this%matrices_, structure)
    
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
    call logfile%print_line(                                    &
       & 'Fractional L2 error in rotated primitive vectors: '// &
       & sqrt(difference/average))
    if (sqrt(difference/average) > 1e-10_dp) then
      call print_line(WARNING//': Error in primitive vectors. &
         &Please check log files.')
      call print_line(difference//' / '//average)
    endif
  endif
end subroutine

! ----------------------------------------------------------------------
! Construct the matrix of force constants for the large supercell, given the
!    dynamical matrices at each q-point.
! ----------------------------------------------------------------------
function reconstruct_force_constants(large_supercell,qpoints, &
   & dynamical_matrices,logfile) result(output)
  implicit none
  
  type(StructureData),   intent(in)    :: large_supercell
  type(QpointData),      intent(in)    :: qpoints(:)
  type(DynamicalMatrix), intent(in)    :: dynamical_matrices(:)
  type(OFile),           intent(inout) :: logfile
  type(ForceConstants)                 :: output
  
  type(RealMatrix), allocatable :: force_constants(:,:)
  
  type(AtomData) :: atom_1
  type(AtomData) :: atom_2
  
  type(IntVector)  :: r
  type(RealVector) :: q
  
  integer :: i,j,k,ialloc
  
  allocate( force_constants( large_supercell%no_atoms,  &
          &                  large_supercell%no_atoms), &
          & stat=ialloc); call err(ialloc)
  force_constants = dblemat(zeroes(3,3))
  
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
        
        ! Add in the contribution to the force constant matrix.
        force_constants(atom_1%id(),atom_2%id()) =                     &
           &   force_constants(atom_1%id(),atom_2%id())                &
           & + real( dynamical_matrices(i)%matrices_( atom_1%prim_id(), &
           &                                         atom_2%prim_id()) &
           &       * exp_2pii(q*r))
      enddo
    enddo
  enddo
  
  output = ForceConstants(large_supercell, force_constants, logfile)
end function

! ----------------------------------------------------------------------
! Rotate a dynamical matrix and set of normal modes onto a new q-point.
! ----------------------------------------------------------------------
! Construct data at q_new from data at q_old, where
!    R . q_old = q_new.
! N.B. the provided q-point should be q_new not q_old.
! The symmetry, S, maps equilibrium position ri to rj+R, and q-point q to q'.
! The +R needs to be corrected for, by multiplying the phase by exp(-iq'.R).
function rotate_modes(input,symmetry,qpoint_from,qpoint_to) result(output)
  implicit none
  
  type(DynamicalMatrix),  intent(in)    :: input
  type(SymmetryOperator), intent(in)    :: symmetry
  type(QpointData),       intent(in)    :: qpoint_from
  type(QpointData),       intent(in)    :: qpoint_to
  type(DynamicalMatrix)                 :: output
  
  ! Rotate dynamical matrix.
  output%matrices_ = rotate_dynamical_matrix( input%matrices_, &
                                            & symmetry,        &
                                            & qpoint_from,     &
                                            & qpoint_to)
  
  ! Rotate normal modes.
  output%complex_modes = transform( input%complex_modes, &
                                  & symmetry,            &
                                  & qpoint_from,         &
                                  & qpoint_to)
end function

! ----------------------------------------------------------------------
! Rotate a dynamical matrix from one q-point to another.
! ----------------------------------------------------------------------
function rotate_dynamical_matrix(input,symmetry,qpoint_from,qpoint_to) &
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
  
  ! Check that the symmetry rotates the q-point as expected.
  if (symmetry * qpoint_from /= qpoint_to) then
    call print_line(CODE_ERROR//': Symmetry does not transform q-points as &
       &expected.')
    call err()
  endif
  
  ! Check that the number of atoms is consistent.
  if (size(input,2)/=no_atoms) then
    call err()
  endif
  
  ! Rotate dynamical matrix.
  allocate(output(no_atoms,no_atoms), stat=ialloc); call err(ialloc)
  do atom_1=1,no_atoms
    atom_1p = symmetry%atom_group * atom_1
    r1 = symmetry%rvector(atom_1)
    do atom_2=1,no_atoms
      atom_2p = symmetry%atom_group * atom_2
      r2 = symmetry%rvector(atom_2)
      
      output(atom_2p,atom_1p) =        &
         & symmetry%cartesian_rotation &
         & * exp_2pii(q*r2)            &
         & * input(atom_2,atom_1)      &
         & * exp_2pii(-q*r1)           &
         & * transpose(symmetry%cartesian_rotation)
    enddo
  enddo
end function

subroutine print_dyn_mat_matrices(input,colours)
  implicit none
  
  type(ComplexMatrix), intent(in)           :: input(:,:)
  type(String),        intent(in), optional :: colours(:,:,:,:)
  
  integer :: no_atoms
  integer :: i,j,k,l
  
  complex(dp) :: thingy(3,3)
  character(50) :: thing
  
  type(String) :: line
  
  no_atoms = size(input,1)
  
  do i=1,no_atoms,no_atoms-1
    do j=1,3
      line = ''
      do k=1,no_atoms,no_atoms-1
        do l=1,3
          thingy = cmplx(input(i,k))
          write(thing,fmt='(f7.3,sp,f6.3,"i")') thingy(j,l)*1e6_dp
          if (present(colours)) then
            line = line//colour(trim(thing),colours(i,j,k,l))
          else
            line = line//trim(thing)
          endif
        enddo
      enddo
      call print_line(line)
    enddo
  enddo
end subroutine

subroutine print_dyn_mat_DynamicalMatrix(input,colours)
  implicit none
  
  type(DynamicalMatrix), intent(in)           :: input
  type(String),          intent(in), optional :: colours(:,:,:,:)
  
  call print_dyn_mat(input%matrices_,colours)
end subroutine

! ----------------------------------------------------------------------
! Returns mode frequencies as a single array.
! ----------------------------------------------------------------------
function frequencies(this) result(output)
  implicit none
  
  class(DynamicalMatrix), intent(in) :: this
  real(dp), allocatable              :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(this%complex_modes)), stat=ialloc); call err(ialloc)
  do i=1,size(this%complex_modes)
    output(i) = this%complex_modes(i)%frequency
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
  complex(dp) :: ma(3,3)
  complex(dp) :: mb(3,3)
  
  real(dp) :: average
  real(dp) :: difference
  
  real(dp)                  :: l2_error
  type(String), allocatable :: colours(:,:,:,:)
  
  integer :: i,j,k,l,ialloc
  
  no_atoms = size(a%matrices_,1)
  if (size(b%matrices_,1)/=no_atoms) then
    call print_line(CODE_ERROR//': dynamical matrices a and b have different &
       &sizes.')
    call err()
  endif
  
  allocate(colours(no_atoms,3,no_atoms,3), stat=ialloc); call err(ialloc)
  
  average = 0.0_dp
  difference = 0.0_dp
  do i=1,no_atoms
    do j=1,no_atoms
      mat_a = a%matrices_(j,i)
      mat_b = b%matrices_(j,i)
      average = average + sum_squares((mat_a+mat_b)/2)
      difference = difference + sum_squares(mat_a-mat_b)
      
      ma = cmplx(mat_a)
      mb = cmplx(mat_b)
      do k=1,3
        do l=1,3
          if (abs(ma(l,k)+mb(l,k))<1e-10_dp) then
            ! If the absolute value is small, ignore the error.
            colours(j,l,i,k) = 'green'
            cycle
          endif
          l2_error = abs(ma(l,k)-mb(l,k))/abs((ma(l,k)+mb(l,k))/2)
          if (l2_error>1e-2_dp) then
            colours(j,l,i,k) = 'red'
          elseif (l2_error>1e-4_dp) then
            colours(j,l,i,k) = 'yellow'
          elseif (l2_error>1e-6_dp) then
            colours(j,l,i,k) = 'magenta'
          elseif (l2_error>1e-8_dp) then
            colours(j,l,i,k) = 'blue'
          elseif (l2_error>1e-10_dp) then
            colours(j,l,i,k) = 'cyan'
          else
            colours(j,l,i,k) = 'green'
          endif
        enddo
      enddo
    enddo
  enddo
  call logfile%print_line('Fractional L2 difference between dynamical &
     &matrices: '//sqrt(difference/average))
  if (sqrt(difference/average)>1e-10_dp) then
    call print_line(WARNING//': Dynamical matrices differ. Please check &
       &log files.')
    call print_line(difference//' /'//average)
    call print_line('First matrix:')
    call print_dyn_mat(a,colours)
    call print_line('Second matrix:')
    call print_dyn_mat(b,colours)
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
        this%matrices_(j,i) = elements(k)%strings(2:4)
      enddo
    enddo
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
  end select
end function

impure elemental function new_DynamicalMatrix_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(DynamicalMatrix)         :: this
  
  this = input
end function
end module
