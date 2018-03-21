module normal_mode_symmetry_module
  use utils_module
  
  use structure_module
  use normal_mode_module
  implicit none
  
  private
  
  public :: rotate_complex_modes
  public :: calculate_symmetry_in_normal_coordinates
  
  interface calculate_symmetry_in_normal_coordinates
    module procedure calculate_symmetry_in_normal_coordinates_qpoint
    module procedure calculate_symmetry_in_normal_coordinates_qpoints
  end interface
contains

! ----------------------------------------------------------------------
! Rotates a complex mode.
! ----------------------------------------------------------------------
! Checks that the symmetry correctly maps qpoint_from to qpoint_to.
impure elemental function rotate_complex_modes(input,symmetry,qpoint_from, &
   & qpoint_to) result(output)
  implicit none
  
  type(ComplexMode),      intent(in) :: input
  type(SymmetryOperator), intent(in) :: symmetry
  type(QpointData),       intent(in) :: qpoint_from
  type(QpointData),       intent(in) :: qpoint_to
  type(ComplexMode)                  :: output
  
  integer :: atom_from
  integer :: atom_to
  
  type(IntVector) :: r
  
  ! Check that the symmetry rotates the q-point as expected.
  if (symmetry * qpoint_from /= qpoint_to) then
    call print_line(CODE_ERROR//': Symmetry does not transform q-points as &
       &expected.')
    call err()
  endif
  
  ! Allocate output, and transfer across all data.
  ! (Displacements need rotating, but everything else stays the same.)
  output = input
  do atom_from=1,size(input%primitive_displacements)
    atom_to = symmetry%prim_atom_group * atom_from
    r = symmetry%prim_rvector(atom_from)
    output%primitive_displacements(atom_to) =       &
       &   symmetry%cartesian_rotation              &
       & * input%primitive_displacements(atom_from) &
       & * exp_2pii(qpoint_to%qpoint*r)
  enddo
end function

! ----------------------------------------------------------------------
! Calculates a symmetry in normal mode co-ordinates.
! ----------------------------------------------------------------------

! Takes q1, {u1}, q2, {u2} and S. Outputs {u2.S.u1}.
function calculate_symmetry_in_normal_coordinates_qpoint(modes,qpoint, &
   & symmetry,logfile) result(output)
  implicit none
  
  type(ComplexMode),      intent(in)    :: modes(:)
  type(QpointData),       intent(in)    :: qpoint
  type(SymmetryOperator), intent(in)    :: symmetry
  type(OFile),            intent(inout) :: logfile
  type(ComplexMatrix)                   :: output
  
  type(ComplexMode), allocatable :: rotated_modes(:)
  complex(dp),       allocatable :: dot_products(:,:)
  
  real(dp) :: check
  
  integer :: i,j,ialloc
  
  ! Calculate the rotated modes, S.u1.
  rotated_modes = rotate_complex_modes( modes,    &
                                      & symmetry, &
                                      & qpoint,   &
                                      & qpoint)
  
  ! Construct the overlap matrix, u2.S.u1.
  allocate( dot_products(size(modes),size(modes)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(modes)
    do j=1,size(modes)
      dot_products(j,i) = modes(j) * rotated_modes(i)
    enddo
  enddo
  
  output = dot_products
  
  ! Check that the symmetry is unitary.
  check = sqrt(sum_squares( output*hermitian(output) &
                        & - cmplxmat(make_identity_matrix(size(modes)))))
  if (check>1e-14_dp) then
    ! This check happens many times, so the 1e-14 limit on logging prevents
    !    clutter.
    call logfile%print_line('Error in unitarity of rotation: ' &
       & //check)
  endif
  if (check>1e-10_dp) then
    call print_line(WARNING//': Rotation between degenerate modes not &
       &unitary. Please try adjusting degenerate_energy. Please check log &
       &files.')
  endif
end function

! Takes {q1}, {u1}, {q2}, {u2} and S. Outputs {u2.S.u1}.
function calculate_symmetry_in_normal_coordinates_qpoints(modes,qpoints, &
   & symmetry,logfile) result(output)
  implicit none
  
  type(ComplexMode),      intent(in)    :: modes(:)
  type(QpointData),       intent(in)    :: qpoints(:)
  type(SymmetryOperator), intent(in)    :: symmetry
  type(OFile),            intent(inout) :: logfile
  type(ComplexMatrix)                   :: output
  
  type(QpointData),  allocatable :: rotated_qpoints(:)
  type(ComplexMode), allocatable :: rotated_modes(:)
  complex(dp),       allocatable :: dot_products(:,:)
  
  real(dp) :: check
  
  integer :: i,j,ialloc
  
  ! Check input sizes are consistent.
  if (size(modes)/=size(qpoints)) then
    call print_line(CODE_ERROR//': modes and q-points do not match.')
    call err()
  endif
  
  ! Rotate q-points.
  allocate(rotated_qpoints(size(qpoints)), stat=ialloc); call err(ialloc)
  do i=1,size(qpoints)
    rotated_qpoints(i) = symmetry * qpoints(i)
  enddo
  
  ! Calculate all rotated modes, S.u1.
  rotated_modes = rotate_complex_modes( modes,    &
                                      & symmetry, &
                                      & qpoints,  &
                                      & rotated_qpoints)
  
  ! Construct the overlap matrix, u2.S.u1.
  allocate( dot_products(size(modes),size(modes)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(modes)
    do j=1,size(modes)
      if (qpoints(j)==rotated_qpoints(i)) then
        dot_products(j,i) = modes(j) * rotated_modes(i)
      else
        dot_products(j,i) = cmplx(0.0_dp,0.0_dp,dp)
      endif
    enddo
  enddo
  
  output = dot_products
  
  ! Check that the symmetry is unitary.
  check = sqrt(sum_squares( output*hermitian(output) &
                        & - cmplxmat(make_identity_matrix(size(modes)))))
  if (check>1e-14_dp) then
    ! This check happens many times, so the 1e-14 limit on logging prevents
    !    clutter.
    call logfile%print_line('Error in unitarity of rotation: ' &
       & //check)
  endif
  if (check>1e-10_dp) then
    call print_line(WARNING//': Rotation between degenerate modes not &
       &unitary. Please try adjusting degenerate_energy. Please check log &
       &files.')
  endif
end function
end module
