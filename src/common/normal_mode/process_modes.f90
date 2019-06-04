! ======================================================================
! Take the modes at a given q-point, and process them w/r/t the q-point,
!    structure and symmetry data.
! ======================================================================
! Fills in:
!    - id
!    - paired_id
!    - qpoint_id
!    - paired_qpoint_id
!    - subspace_id
! In general, a the modes in a given subspace will be returned as a linear
!    combination of one another.
module process_modes_module
  use utils_module
  
  use structure_module
  
  use complex_mode_module
  use complex_mode_symmetry_module
  implicit none
  
  private
  
  public :: process_modes
  
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

! Process modes.
function process_modes(input,structure,qpoint,subspace_id) result(output)
  implicit none
  
  type(ComplexMode),   intent(in) :: input(:)
  type(StructureData), intent(in) :: structure
  type(QpointData),    intent(in) :: qpoint
  integer,             intent(in) :: subspace_id
  type(ComplexMode), allocatable  :: output(:)
  
  integer, allocatable :: subspace_ids(:)
  integer, allocatable :: subspace_id_set(:)
  
  real(dp) :: energy_difference
  
  ! Symmetry data.
  type(SymmetryOperator), allocatable :: symmetries(:)
  
  integer, allocatable :: states(:)
  
  complex(dp) :: symmetry(2,2)
  
  complex(dp), allocatable :: overlap(:,:)
  
  integer :: i,j,k,ialloc
  
  ! Copy the input to the output.
  output = input
  
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
             & abs(sum(output(i)%unit_vector*output(j)%unit_vector))>1e-3 &
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
          if (abs(symmetry(1,2))>1e-3 .or. abs(symmetry(2,1))>1e-3) then
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
    call print_line(WARNING//': Unable to lift degeneracies using symmetry.')
    call print_line('q-point '//qpoint%id//': '//qpoint%qpoint)
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
end module
