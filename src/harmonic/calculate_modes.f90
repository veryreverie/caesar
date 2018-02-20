! ======================================================================
! Routines to calculate normal modes from dynamical matrices.
! ======================================================================
module calculate_modes_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use normal_mode_module
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
  use structure_module
  use atom_module
  use eigenstuff_module
  implicit none
  
  type(ComplexMatrix), intent(in) :: matrices(:,:)
  type(StructureData), intent(in) :: structure
  type(ComplexMode), allocatable  :: output(:)
  
  complex(dp), allocatable :: dyn_mat(:,:)
  
  type(HermitianEigenstuff), allocatable :: estuff(:)
  
  type(ComplexVector) :: displacement
  
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
   &degenerate_energy,degeneracy_id,logfile) result(output)
  use utils_module, only : sum_squares
  use structure_module
  use qpoints_module
  use atom_module
  use eigenstuff_module
  use ofile_module
  use symmetry_module
  use logic_module
  implicit none
  
  type(ComplexMatrix), intent(in)    :: matrices(:,:)
  type(StructureData), intent(in)    :: structure
  type(QpointData),    intent(in)    :: qpoint
  real(dp),            intent(in)    :: degenerate_energy
  integer,             intent(in)    :: degeneracy_id
  type(OFile),         intent(inout) :: logfile
  type(ComplexMode), allocatable     :: output(:)
  
  real(dp) :: energy_difference
  
  ! Symmetry data.
  type(SymmetryOperator), allocatable :: symmetries(:)
  
  integer, allocatable :: states(:)
  
  ! Error checking variables.
  type(ComplexMatrix) :: symmetry
  real(dp)            :: check
  
  integer :: i,j,k,ialloc
  
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
  symmetries = structure%symmetries( filter( structure%symmetries, &
                                           & leaves_q_invariant) )
  
  ! Assign degeneracy ids, which are equal if two states are degenerate,
  !    and different if they are not.
  output(1)%degeneracy_id = degeneracy_id
  do i=2,size(output)
    energy_difference = abs(output(i)%frequency-output(i-1)%frequency)
    if (energy_difference<degenerate_energy) then
      output(i)%degeneracy_id = output(i-1)%degeneracy_id
      if (energy_difference>degenerate_energy/2) then
        call print_line(WARNING//': Degenerate energies within a factor of &
           &two of degenerate_energy at q-point '//qpoint%qpoint//'.')
      endif
    else
      output(i)%degeneracy_id = output(i-1)%degeneracy_id + 1
      if (energy_difference<degenerate_energy*2) then
        call print_line(WARNING//': Non-degenerate energies within a factor &
           &of two of degenerate_energy at q-point '//qpoint%qpoint//'.')
      endif
    endif
  enddo
  
  ! Loop over degeneracy ids, checking each degenerate subspace, and lifting
  !    degeneracy using symmetry operators.
  do i=degeneracy_id,output(size(output))%degeneracy_id
    ! Find the set of states with degeneracy id i.
    states = filter(output%degeneracy_id==i)
    
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
    do j=1,size(symmetries)
      symmetry = calculate_symmetry_in_normal_coordinates( qpoint,         &
                                                         & output(states), &
                                                         & qpoint,         &
                                                         & output(states), &
                                                         & symmetries(j),  &
                                                         & logfile)
      check = sqrt(sum_squares( symmetry*hermitian(symmetry) &
                            & - cmplxmat(make_identity_matrix(size(states)))))
      if (check>1e-2_dp) then
        call print_line(ERROR//': Rotation between degenerate modes at &
           &q-point '//qpoint%qpoint//' not unitary. Please adjust &
           &degenerate_energy so that it accurately captures degeneracies.')
        stop
      endif
    enddo
    
    ! Lift each degeneracy using symmetry operators.
    if (size(states)>1) then
      output(states) = lift_degeneracies( output(states), &
                                        & symmetries,     &
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
recursive function lift_degeneracies(input,symmetries,qpoint,logfile) &
   & result(output)
  use symmetry_module
  use qpoints_module
  use eigenstuff_module
  use phase_module
  use logic_module
  use ofile_module
  implicit none
  
  type(ComplexMode),      intent(in)    :: input(:)
  type(SymmetryOperator), intent(in)    :: symmetries(:)
  type(QpointData),       intent(in)    :: qpoint
  type(OFile),            intent(inout) :: logfile
  type(ComplexMode), allocatable        :: output(:)
  
  integer :: no_modes
  
  ! Symmetry information.
  integer :: sym_id
  integer :: order
  
  type(SymmetryOperator) :: first_symmetry
  type(ComplexMatrix) :: symmetry
  
  type(UnitaryEigenstuff), allocatable :: estuff(:)
  
  type(SymmetryOperator), allocatable :: commuting_symmetries(:)
  
  ! Error checking.
  real(dp) :: check
  
  integer :: i,j,ialloc
  
  ! All q-point data and eigenvalues will be unchanged.
  ! Copy over all data, and only change that which changes.
  output = input
  
  ! Count the number of modes.
  no_modes = size(input)
  
  if (no_modes==1) then
    call print_line(CODE_ERROR//': Trying to lift the degeneracy of only one &
       &state.')
    call err()
  endif
  
  if (size(symmetries)==0) then
    call print_line(ERROR//': Unable to lift degeneracies using symmetry. &
       &Please try reducing degenerate_energy.')
    stop
  endif

  first_symmetry = symmetries(1)
  
  ! Construct the first symmetry in normal mode co-ordinates.
  symmetry = calculate_symmetry_in_normal_coordinates( qpoint,         &
                                                     & input,          &
                                                     & qpoint,         &
                                                     & input,          &
                                                     & first_symmetry, &
                                                     & logfile)
  
  ! Calculate the order of the first symmetry, n s.t. S^n=I.
  order = calculate_symmetry_order(first_symmetry, qpoint)
  
  ! Diagonalise the first symmetry, and construct diagonalised displacements.
  ! Only transform displacements if this symmetry lifts degeneracy.
  estuff = diagonalise_unitary(symmetry,order)
  if (estuff(1)%eval/=estuff(no_modes)%eval) then
    do i=1,no_modes
      output(i)%primitive_displacements = cmplxvec(zeroes(3))
      do j=1,no_modes
        output(i)%primitive_displacements = output(i)%primitive_displacements &
                                        & + estuff(i)%evec(j)                 &
                                        & * input(j)%primitive_displacements
      enddo
      
      check = abs(l2_norm(output(i))-1)
      call logfile%print_line('Error in mode normalisation: '// &
         & check)
      if (check>1e-10_dp) then
        call print_line(WARNING//': Error in mode normalisation. Please check &
           &log files.')
      endif
    enddo
  endif
  
  ! Lift remaining degeneracies using remaining symmetries.
  i = 1
  do while(i<=no_modes)
    ! The range i:j is the set of modes degenerate with mode i under the
    !    first symmetry.
    j = last(estuff%eval==estuff(i)%eval)
    
    if (j<i) then
      call err()
    endif
    
    if (j>i) then
      ! Select only the symmetry operators which commute with the
      !    first symmetry.
      commuting_symmetries = symmetries(filter(symmetries,commutes_with_first))
      
      ! Lift further degeneracies using symmetries which commute with the
      !    first symmetry, not including the first symmetry.
      output(i:j) = lift_degeneracies( output(i:j),              &
                                     & commuting_symmetries(2:), &
                                     & qpoint,                   &
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

! --------------------------------------------------
! Calculates a symmetry in normal mode co-ordinates.
! --------------------------------------------------
! Takes q1, {u1}, q2, {u2} and S. Outputs {u2.S.u1}.
function calculate_symmetry_in_normal_coordinates(qpoint_from,modes_from, &
   & qpoint_to,modes_to,symmetry,logfile) result(output)
  use utils_module, only : sum_squares
  use qpoints_module
  use symmetry_module
  use ofile_module
  implicit none
  
  type(QpointData),       intent(in)    :: qpoint_from
  type(ComplexMode),      intent(in)    :: modes_from(:)
  type(QpointData),       intent(in)    :: qpoint_to
  type(ComplexMode),      intent(in)    :: modes_to(:)
  type(SymmetryOperator), intent(in)    :: symmetry
  type(OFile),            intent(inout) :: logfile
  type(ComplexMatrix)                   :: output
  
  integer :: no_modes
  
  type(ComplexMode), allocatable :: rotated_modes(:)
  complex(dp),       allocatable :: dot_products(:,:)
  
  real(dp) :: check
  
  integer :: i,j,ialloc
  
  no_modes = size(modes_from)
  if (no_modes/=size(modes_to)) then
    call print_line(CODE_ERROR//': Inconsistent number of modes.')
    call err()
  endif
  
  ! Construct the symmetry in normal mode co-ordinates.
  rotated_modes = rotate_complex_modes( modes_from,  &
                                      & symmetry,    &
                                      & qpoint_from, &
                                      & qpoint_to)
  
  ! Construct the overlap matrix, u2.S.u1.
  allocate( dot_products(no_modes,no_modes), &
          & stat=ialloc); call err(ialloc)
  do i=1,no_modes
    do j=1,no_modes
      dot_products(j,i) = modes_to(j) * rotated_modes(i)
    enddo
  enddo
  
  output = dot_products
  
  ! Check that the symmetry is unitary.
  check = sqrt(sum_squares( output*hermitian(output) &
                        & - cmplxmat(make_identity_matrix(no_modes))))
  call logfile%print_line('Error in unitarity of rotation: ' &
     & //check)
  if (check>1e-10_dp) then
    call print_line(WARNING//': Rotation between degenerate modes not &
       &unitary. Please try adjusting degenerate_energy. Please check log &
       &files.')
  endif
end function

! ----------------------------------------------------------------------
! Calculates the order of a symmetry operation at a give q-point.
! The order is the smallest integer n>0 s.t. S^n=I, where I is the identity.
! ----------------------------------------------------------------------
! N.B. This calculation assumes that the symmetry changes relative phases.
!    If this is not the case, then n will be too large.
function calculate_symmetry_order(symmetry,qpoint) result(output)
  use utils_module, only : lcm
  use symmetry_module
  use qpoints_module
  implicit none
  
  type(SymmetryOperator), intent(in) :: symmetry
  type(QpointData),       intent(in) :: qpoint
  integer                            :: output
  
  type(IntMatrix) :: identity
  type(IntMatrix) :: rotation ! R^n.
  
  integer :: i
  
  identity = make_identity_matrix(3)
  rotation = identity
  
  ! Calculate n s.t. R^n=I.
  output = 0
  do i=1,6
    rotation = symmetry%rotation * rotation
    if (rotation==identity) then
      output = i
      exit
    endif
  enddo
  
  if (output==0) then
    call print_line(CODE_ERROR//': Unable to find order of symmetry.')
    call err()
  endif
  
  ! Assume that the symmetry changes relative phases.
  output = lcm(output,qpoint%min_sc_size())
end function
end module
