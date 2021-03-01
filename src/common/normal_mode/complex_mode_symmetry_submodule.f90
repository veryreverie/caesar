submodule (caesar_complex_mode_symmetry_module) caesar_complex_mode_symmetry_submodule
  use caesar_normal_mode_module
contains

module procedure calculate_symmetry_in_normal_coordinates
  type(QpointData),  allocatable :: transformed_qpoints(:)
  type(ComplexMode), allocatable :: transformed_modes(:)
  complex(dp),       allocatable :: dot_products(:,:)
  
  integer :: i,j,ialloc
  
  ! Check input sizes are consistent.
  if (size(modes)/=size(qpoints)) then
    call print_line(CODE_ERROR//': modes and q-points do not match.')
    call err()
  elseif (any(modes%qpoint_id/=qpoints%id)) then
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
  
  output = mat(dot_products)
end procedure
end submodule
