! ======================================================================
! Methods for converting symmetries into complex normal mode co-ordinates.
! ======================================================================
module complex_mode_symmetry_submodule
  use utils_module
  
  use structure_module
  
  use complex_mode_submodule
  implicit none
  
  private
  
  public :: calculate_symmetry_in_normal_coordinates
  
  interface calculate_symmetry_in_normal_coordinates
    module procedure calculate_symmetry_in_normal_coordinates_qpoint
    module procedure calculate_symmetry_in_normal_coordinates_qpoints
  end interface
contains
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
  
  integer :: i,j,ialloc
  
  ! Calculate the rotated modes, S.u1.
  rotated_modes = transform( modes,    &
                           & symmetry, &
                           & qpoint,   &
                           & qpoint)
  
  ! Construct the overlap matrix, u2.S.u1.
  allocate( dot_products(size(modes),size(modes)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(modes)
    do j=1,size(modes)
      dot_products(j,i) = sum( conjg(modes(j)%mass_weighted_vector) &
                           & * rotated_modes(i)%mass_weighted_vector)
    enddo
  enddo
  
  output = dot_products
  
  ! Check that the symmetry is unitary.
  call check_unitary(output,'symmetry in normal co-ordinates',logfile)
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
  rotated_modes = transform( modes,    &
                           & symmetry, &
                           & qpoints,  &
                           & rotated_qpoints)
  
  ! Construct the overlap matrix, u2.S.u1.
  allocate( dot_products(size(modes),size(modes)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(modes)
    do j=1,size(modes)
      if (qpoints(j)==rotated_qpoints(i)) then
        dot_products(j,i) = sum( conjg(modes(j)%mass_weighted_vector) &
                             & * rotated_modes(i)%mass_weighted_vector)
      else
        dot_products(j,i) = cmplx(0.0_dp,0.0_dp,dp)
      endif
    enddo
  enddo
  
  output = dot_products
  
  ! Check that the symmetry is unitary.
  call check_unitary(output,'symmetry in normal co-ordinates',logfile)
end function
end module
