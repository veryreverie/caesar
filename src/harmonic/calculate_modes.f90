! ======================================================================
! Routines to calculate normal modes from dynamical matrices.
! ======================================================================
module calculate_modes_module
  use common_module
  implicit none
  
  private
  
  public :: calculate_modes
  
  interface calculate_modes
    module procedure calculate_modes_interpolated
    module procedure calculate_modes_calculated
  end interface
contains

! ----------------------------------------------------------------------
! Calculates complex modes by diagonalising a dynamical matrix.
! ----------------------------------------------------------------------
! N.B. Structure may be any supercell.

! Calculate modes for a q-point other than one of the calculated q-points.
function calculate_modes_interpolated(matrices,structure) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: matrices(:,:)
  type(StructureData), intent(in) :: structure
  type(ComplexMode), allocatable  :: output(:)
  
  complex(dp), allocatable :: dyn_mat(:,:)
  
  type(HermitianEigenstuff), allocatable :: estuff(:)
  
  integer :: i,j,k,ialloc
  
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
  
  ! Calculate normal mode frequencies and displacements.
  !          V = sum_i[ 0.5 * freq_i**2 * u_i**2]
  ! -> F = -2V = sum_i[ - freq_i**2 * u_i**2 ]
  allocate( output(structure%no_modes_prim),   &
          & stat=ialloc); call err(ialloc)
  do i=1,structure%no_modes_prim
    
    ! The eigenvalues are in descending order, but the normal modes should
    !    be in ascending order of frequency. i->k reverses the order.
    k = structure%no_modes_prim - i + 1
    
    ! Initialise id and paired_id to 0.
    output(i)%id        = 0
    output(i)%paired_id = 0
    
    if (estuff(k)%eval>=0.0_dp) then
      ! Unstable mode.
      output(i)%frequency = - sqrt(estuff(k)%eval)
    else
      ! Stable mode.
      output(i)%frequency = sqrt(- estuff(k)%eval)
    endif
    
    output(i)%soft_mode = output(i)%frequency < -1.0e-6_dp
    output(i)%translational_mode = .false.
    
    ! Calculate displacements in the primitive cell,
    !    which are the non-mass-reduced eigenvectors of the dynamical matrix.
    allocate( output(i)%primitive_displacements(structure%no_atoms_prim), &
            & stat=ialloc); call err(ialloc)
    do j=1,structure%no_atoms_prim
      output(i)%primitive_displacements(j) = estuff(k)%evec(3*j-2:3*j) &
                                         & * sqrt(structure%atoms(j)%mass())
      
    enddo
    
    ! Re-normalise modes, now in non-mass-reduced co-ordinates.
    output(i)%primitive_displacements = output(i)%primitive_displacements &
                                    & / l2_norm(output(i))
  enddo
end function

! Calculate modes for one of the calculated q-points.
! Will lift degeneracies using symmetries.
function calculate_modes_calculated(matrices,structure,qpoint, &
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
  real(dp)            :: check
  
  integer :: i,j
  
  ! Calculate normal modes as if at an arbitrary q-point.
  output = calculate_modes(matrices,structure)
  
  ! Identify purely translational modes (at the gamma-point only).
  output%translational_mode = .false.
  if (is_int(qpoint%qpoint)) then
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
  
  integer :: no_modes
  
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
  symmetry_is_sin = structure%symmetry_is_sin(symmetry_ids(1))
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
  
  ! If the symmetry lifts degeneracy (has multiple phases), then rotate the
  !    input displacements into the symmetry's eigenbasis.
  if (any(phases_int/=phases_int(1))) then
    do i=1,size(input)
      do j=1,size(output(i)%primitive_displacements)
        output(i)%primitive_displacements(j) = cmplxvec(zeroes(3))
        do k=1,size(estuff(i)%evec)
          output(i)%primitive_displacements(j) =    &
             & output(i)%primitive_displacements(j) &
             & + estuff(i)%evec(k) * input(k)%primitive_displacements(j)
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
  ! Captures first_symmetry.
  function commutes_with_first(input) result(output)
    implicit none
    
    class(*), intent(in) :: input
    logical              :: output
    
    select type(input); type is(SymmetryOperator)
      output = operators_commute(input,first_symmetry)
    end select
  end function
end function
end module
