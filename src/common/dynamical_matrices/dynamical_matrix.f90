! ======================================================================
! The forces between atoms at a given q-point.
! ======================================================================
module dynamical_matrix_module
  use utils_module
  
  use structure_module
  use normal_mode_module
  
  use min_images_module
  use force_constants_module
  implicit none
  
  private
  
  public :: DynamicalMatrix
  public :: compare_dynamical_matrices
  public :: transform_modes
  public :: conjg
  public :: reconstruct_force_constants
  
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
    module procedure new_ComplexMode_interpolated
    module procedure new_ComplexMode_calculated
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
          & sym_id(no_copies), &
          & sup_id(no_copies), &
          & stat=ialloc); call err(ialloc)
  k = 0
  do i=1,size(structure%symmetries)
    q_prime = structure%symmetries(i)%inverse_transform(qpoint)
    do j=1,size(supercells)
      if (is_copy(j,i)) then
        k = k+1
        matrix = calculate_dynamical_matrix( dblevec(q_prime%qpoint), &
                                           & supercells(j),           &
                                           & force_constants(j))
        matrices(:,:,k) = transform_dynamical_matrix( &
                           & matrix,                  &
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
  this%complex_modes = ComplexMode(this%matrices_, supercell)
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
         & * sum(exp_2pii(q*rvectors))                 &
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
      call print_line(difference//' / '//average)
    endif
  endif
end subroutine

! ----------------------------------------------------------------------
! Calculates complex modes by diagonalising a dynamical matrix.
! ----------------------------------------------------------------------
! N.B. Structure may be any supercell.

! Calculate modes for a q-point other than one of the calculated q-points.
function new_ComplexMode_interpolated(matrices,structure) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: matrices(:,:)
  type(StructureData), intent(in) :: structure
  type(ComplexMode), allocatable  :: output(:)
  
  complex(dp), allocatable :: dyn_mat(:,:)
  
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

! Calculate modes for one of the calculated q-points.
! Will lift degeneracies using symmetries.
function new_ComplexMode_calculated(matrices,structure,qpoint, &
   &degenerate_energy,subspace_id,logfile) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in)    :: matrices(:,:)
  type(StructureData), intent(in)    :: structure
  type(QpointData),    intent(in)    :: qpoint
  real(dp),            intent(in)    :: degenerate_energy
  integer,             intent(in)    :: subspace_id
  type(OFile),         intent(inout) :: logfile
  type(ComplexMode), allocatable     :: output(:)
  
  real(dp) :: energy_difference
  
  ! Symmetry data.
  integer, allocatable :: symmetry_ids(:)
  
  integer, allocatable :: states(:)
  
  ! Error checking variables.
  type(ComplexMatrix) :: symmetry
  
  ! Variables for pairing up modes.
  complex(dp), allocatable :: pair_overlap(:)
  complex(dp)              :: phase
  
  integer :: i,j,k,ialloc
  
  ! Calculate normal modes as if at an arbitrary q-point.
  output = ComplexMode(matrices,structure)
  
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
  if (is_gvector(qpoint)) then
    ! Find the three modes with the minimum abs(frequency).
    do i=1,3
      j = minloc( abs(output%frequency), &
                & dim=1,                 &
                & mask=.not.output%translational_mode)
      output(j)%translational_mode = .true.
    enddo
  endif
  
  ! Identify the symmetries which map the q-point to itself.
  symmetry_ids = filter(structure%symmetries, leaves_q_invariant)
  
  ! Assign degeneracy ids, which are equal if two states are degenerate,
  !    and different if they are not.
  output(1)%subspace_id = subspace_id
  do i=2,size(output)
    energy_difference = abs(output(i)%frequency-output(i-1)%frequency)
    if (energy_difference<degenerate_energy) then
      output(i)%subspace_id = output(i-1)%subspace_id
      if (energy_difference>degenerate_energy/2) then
        call print_line(WARNING//': Degenerate energies within a factor of &
           &two of degenerate_energy at q-point '//qpoint%qpoint//'.')
      endif
    else
      output(i)%subspace_id = output(i-1)%subspace_id + 1
      if (energy_difference<degenerate_energy*2) then
        call print_line(WARNING//': Non-degenerate energies within a factor &
           &of two of degenerate_energy at q-point '//qpoint%qpoint//'.')
      endif
    endif
  enddo
  
  ! Loop over degeneracy ids, checking each degenerate subspace, and lifting
  !    degeneracy using symmetry operators.
  do i=subspace_id,output(size(output))%subspace_id
    ! Find the set of states with degeneracy id i.
    states = filter(output%subspace_id==i)
    
    ! Check that degenerate states are consistent, i.e. that if two states
    !    are both degenerate with a third state that they are also degenerate
    !    with one another.
    if ( maxval(output(states)%frequency) - minval(output(states)%frequency) &
     & > degenerate_energy ) then
      call print_line(ERROR//': Modes inconsistently degenerate. Please try &
         &adjusting degenerate_energy.')
      call err()
    endif
    
    ! Set the frequencies of each degenerate state to
    !    the average of their frequencies.
    output(states)%frequency = sum(output(states)%frequency) / size(states)
    
    ! Check that all symmetries map the degenerate subspace onto itself.
    do j=1,size(symmetry_ids)
      symmetry = calculate_symmetry_in_normal_coordinates( &
                  & output(states),                        &
                  & qpoint,                                &
                  & structure%symmetries(symmetry_ids(j)), &
                  & logfile)
      call check_unitary(symmetry,'degenerate modes symmetry',logfile)
    enddo
    
    ! Lift each degeneracy using symmetry operators.
    if (size(states)>1) then
      output(states) = lift_degeneracies( output(states), &
                                        & structure,      &
                                        & symmetry_ids,   &
                                        & qpoint,         &
                                        & logfile)
    endif
    
    ! If this q-point is its own pair, identify and pair up conjugate modes.
    if (qpoint%id==qpoint%paired_qpoint_id) then
      allocate(pair_overlap(size(states)), stat=ialloc); call err(ialloc)
      do j=1,size(states)
        ! Calculate the overlap of the conjugate of mode j with the other
        !    degenerate modes.
        ! N.B. since the dot product would normaly involve a conjugate,
        !    the two conjugates cancel.
        do k=1,size(states)
          pair_overlap(k) = sum( output(states(j))%unit_vector &
                             & * output(states(k))%unit_vector )
        enddo
        
        ! conjg(mode(k)) should be equal to a mode (down to a phase change),
        !    and have no overlap with any other mode.
        ! Check that this is true, and identify the one mode.
        if (count(abs(pair_overlap)<1e-2)/=size(pair_overlap)-1) then
          call print_line(ERROR//': Modes at qpoint '//i//' do not transform &
             &as expected under parity.')
          call print_line(pair_overlap)
          call err()
        endif
        
        ! Get the position of the paired mode in both
        !    the list of degenerate modes and the list of modes.
        k = first(abs(abs(pair_overlap)-1)<1e-2)
        
        if (k==j) then
          ! This state is its own conjugate; rotate the phase so that
          !    the vector is real, and set the imaginary part to zero.
          phase = sqrt(pair_overlap(j))
          phase = phase / abs(phase)
          output(states(j))%unit_vector = &
             & cmplxvec(real(output(states(j))%unit_vector/phase))
        elseif (k<j) then
          ! Set the paired state to be the conjugate of this state.
          output(states(j))%paired_id = output(states(k))%id
          output(states(k))%paired_id = output(states(j))%id
          output(states(j))%unit_vector = conjg(output(states(k))%unit_vector)
        endif
      enddo
      deallocate(pair_overlap, stat=ialloc); call err(ialloc)
    endif
  enddo
contains
  ! Lambda for identifying if a symmetry leaves the q-point invariant.
  ! i.e. if R*q=q, modulo G-vector translations.
  ! Captures qpoint.
  function leaves_q_invariant(input) result(output)
    implicit none
    
    class(*), intent(in) :: input
    logical              :: output
    
    select type(input); type is(SymmetryOperator)
      output = input*qpoint == qpoint
    end select
  end function
end function

! --------------------------------------------------
! Recursively lifts degeneracies using symmetries.
! --------------------------------------------------
! Input must be a list of degenerate modes.
! Symmetries must all take the q-point to itself.
recursive function lift_degeneracies(input,structure,symmetry_ids, &
   & qpoint,logfile) result(output)
  implicit none
  
  type(ComplexMode),      intent(in)    :: input(:)
  type(StructureData),    intent(in)    :: structure
  integer,                intent(in)    :: symmetry_ids(:)
  type(QpointData),       intent(in)    :: qpoint
  type(OFile),            intent(inout) :: logfile
  type(ComplexMode), allocatable        :: output(:)
  
  ! Symmetry information.
  integer :: order
  
  type(SymmetryOperator) :: first_symmetry
  type(ComplexMatrix) :: symmetry
  
  type(HermitianEigenstuff), allocatable :: estuff(:)
  
  real(dp), allocatable :: phases_real(:)
  integer,  allocatable :: phases_int(:)
  
  type(SymmetryOperator), allocatable :: symmetries(:)
  integer, allocatable :: commuting_symmetry_ids(:)
  
  logical :: symmetry_is_sin
  
  integer :: i,j,k,ialloc
  
  ! All q-point data and eigenvalues will be unchanged.
  ! Copy over all data, and only change that which changes.
  output = input
  
  if (size(input)==1) then
    call print_line(CODE_ERROR//': Trying to lift the degeneracy of only one &
       &state.')
    call err()
  endif
  
  if (size(symmetry_ids)==0) then
    call print_line(ERROR//': Unable to lift degeneracies using symmetry. &
       &Please try reducing degenerate_energy.')
    stop
  endif

  symmetries = structure%symmetries(symmetry_ids)
  first_symmetry = symmetries(1)
  
  ! Construct the first symmetry in normal mode co-ordinates.
  symmetry = calculate_symmetry_in_normal_coordinates( input,          &
                                                     & qpoint,         &
                                                     & first_symmetry, &
                                                     & logfile)
  
  ! Instead of directly calculating the eigenstuff of the unitary symmetry
  !    matrices {U}, it is more stable to calculate the eigenstuff of the
  !    Hermitian matrices {C=(U+U^T)/2} and {S=(U-U^T)/2i}.
  ! The eigenvalues of U are e^(2*pi*i*j/n), so the eigenvalues of C and S are
  !    cos(2*pi*j/n) and sin(2*pi*j/n) respectively.
  ! Arbitrarily, if this matrix's conjugate has a later id
  !    (or is the same matrix) then C is considered, else S is considered.
  symmetry_is_sin = first_symmetry%inverse_symmetry_id < first_symmetry%id
  if (symmetry_is_sin) then
    symmetry = (symmetry - hermitian(symmetry))/cmplx(0.0_dp,2.0_dp,dp)
  else
    symmetry = (symmetry + hermitian(symmetry))/2.0_dp
  endif
  
  ! Calculate the order of the first symmetry, n s.t. S^n=I.
  order = first_symmetry%symmetry_order(qpoint)
  
  ! Diagonalise the first symmetry, and construct diagonalised displacements.
  ! Only transform displacements if this symmetry lifts degeneracy.
  estuff = diagonalise_hermitian(symmetry)
  
  ! Work out the phases of the eigenvalues (2*j in e^(2*pi*i*j/order)).
  allocate( phases_real(size(input)), &
          & phases_int(size(input)),  &
          & stat=ialloc); call err(ialloc)
  do i=1,size(input)
    ! Correct for numerical errors taking the eigenvalue outside the range
    !    [-1,1].
    if (abs(estuff(i)%eval)>1.01_dp) then
      call print_line(ERROR//': Symmetry eigenvalue outside the range [-1,1].')
      call err()
    endif
    
    if (estuff(i)%eval>1.0_dp) then
      estuff(i)%eval = 1.0_dp
    elseif (estuff(i)%eval<-1.0_dp) then
      estuff(i)%eval = -1.0_dp
    endif
    
    if (symmetry_is_sin) then
      phases_real(i) = asin(estuff(i)%eval)*order/PI
    else
      phases_real(i) = acos(estuff(i)%eval)*order/PI
    endif
    phases_int(i) = nint(phases_real(i))
    if (abs(phases_int(i)-phases_real(i))>0.1_dp) then
      call print_line(ERROR//': Symmetry with non-integer phase eigenvalue.')
      call err()
    endif
    phases_int(i) = modulo(phases_int(i),2*order)
  enddo
  
  ! If the symmetry lifts degeneracy (has multiple phases), then transform the
  !    input vectors into the symmetry's eigenbasis.
  if (any(phases_int/=phases_int(1))) then
    do i=1,size(input)
      do j=1,size(output(i)%unit_vector)
        output(i)%unit_vector(j) = cmplxvec(zeroes(3))
        do k=1,size(estuff(i)%evec)
          output(i)%unit_vector(j) = output(i)%unit_vector(j) &
                                 & + estuff(i)%evec(k)*input(k)%unit_vector(j)
        enddo
      enddo
    enddo
  endif
  
  ! Lift remaining degeneracies using remaining symmetries.
  i = 1
  do while(i<=size(input))
    ! The range i:j is the set of modes degenerate with mode i under the
    !    first symmetry.
    j = last(phases_int==phases_int(i))
    
    if (j<i) then
      call err()
    endif
    
    if (j>i) then
      ! Select only the symmetry operators which commute with the
      !    first symmetry.
      commuting_symmetry_ids = symmetry_ids(filter( symmetries, &
                                                  & commutes_with_first))
      
      ! Lift further degeneracies using symmetries which commute with the
      !    first symmetry, not including the first symmetry.
      output(i:j) = lift_degeneracies( output(i:j),                &
                                     & structure,                  &
                                     & commuting_symmetry_ids(2:), &
                                     & qpoint,                     &
                                     & logfile)
    endif
    i = j+1
  enddo
contains
  ! Lambda for determining whether or not a symmetry commutes with the first
  !    symmetry.
  ! Captures:
  !    - first_symmetry
  !    - qpoint
  function commutes_with_first(input) result(output)
    implicit none
    
    class(*), intent(in) :: input
    logical              :: output
    
    select type(input); type is(SymmetryOperator)
      output = operators_commute(input,first_symmetry,qpoint)
    end select
  end function
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
        force_constants(atom_1%id(),atom_2%id()) =                       &
           &   force_constants(atom_1%id(),atom_2%id())                  &
           & + real( dynamical_matrices(i)%matrices_( atom_1%prim_id(),  &
           &                                          atom_2%prim_id() ) &
           &       * exp_2pii(-q*r)                                    )
      enddo
    enddo
  enddo
  
  output = ForceConstants(large_supercell, force_constants, logfile)
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
