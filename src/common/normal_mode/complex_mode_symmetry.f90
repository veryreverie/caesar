! ======================================================================
! Methods for converting symmetries into complex normal mode co-ordinates.
! ======================================================================
module complex_mode_symmetry_module
  use utils_module
  
  use structure_module
  
  use complex_mode_module
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
   & symmetry) result(output)
  implicit none
  
  type(ComplexMode),      intent(in)    :: modes(:)
  type(QpointData),       intent(in)    :: qpoint
  type(SymmetryOperator), intent(in)    :: symmetry
  type(ComplexMatrix)                   :: output
  
  type(ComplexMode), allocatable :: transformed_modes(:)
  complex(dp),       allocatable :: dot_products(:,:)
  
  integer :: i,j,ialloc
  
  ! Calculate the transformed modes, S.u1.
  transformed_modes = transform( modes,    &
                               & symmetry, &
                               & qpoint,   &
                               & qpoint    )
  
  ! Construct the overlap matrix, u2.S.u1.
  allocate( dot_products(size(modes),size(modes)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(modes)
    do j=1,size(modes)
      dot_products(j,i) = sum( conjg(modes(j)%unit_vector)      &
                           & * transformed_modes(i)%unit_vector )
    enddo
  enddo
  
  output = dot_products
end function

! Takes {q1}, {u1}, {q2}, {u2} and S. Outputs {u2.S.u1}.
function calculate_symmetry_in_normal_coordinates_qpoints(modes,qpoints, &
   & symmetry) result(output)
  implicit none
  
  type(ComplexMode),      intent(in)    :: modes(:)
  type(QpointData),       intent(in)    :: qpoints(:)
  type(SymmetryOperator), intent(in)    :: symmetry
  type(ComplexMatrix)                   :: output
  
  type(QpointData),  allocatable :: transformed_qpoints(:)
  type(ComplexMode), allocatable :: transformed_modes(:)
  complex(dp),       allocatable :: dot_products(:,:)
  
  integer :: i,j,ialloc
  
  ! Check input sizes are consistent.
  if (size(modes)/=size(qpoints)) then
    call print_line(CODE_ERROR//': modes and q-points do not match.')
    call err()
  endif
  
  ! Transform q-points.
  allocate(transformed_qpoints(size(qpoints)), stat=ialloc); call err(ialloc)
  do i=1,size(qpoints)
    transformed_qpoints(i) = symmetry * qpoints(i)
  enddo
  
  ! Calculate all transformed modes, S.u1.
  transformed_modes = transform( modes,              &
                               & symmetry,           &
                               & qpoints,            &
                               & transformed_qpoints )
  
  ! Construct the overlap matrix, u2.S.u1.
  allocate( dot_products(size(modes),size(modes)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(modes)
    do j=1,size(modes)
      if (qpoints(j)==transformed_qpoints(i)) then
        dot_products(j,i) = sum( conjg(modes(j)%unit_vector)      &
                             & * transformed_modes(i)%unit_vector )
      else
        dot_products(j,i) = cmplx(0.0_dp,0.0_dp,dp)
      endif
    enddo
  enddo
  
  output = dot_products
end function
end module
