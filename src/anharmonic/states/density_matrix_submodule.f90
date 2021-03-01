submodule (caesar_density_matrix_module) caesar_density_matrix_submodule
  use caesar_states_module
contains

module procedure new_DensityMatrix
  integer               :: no_states
  real(dp), allocatable :: coefficients(:,:)
  
  type(IntArray1d), allocatable :: selected_couplings(:)
  
  integer :: no_elements
  
  integer :: i,j,k,ialloc
  
  output%wavevector = wavevector
  
  if (size(states)==0) then
    return
  endif
  
  no_states = size(states(1)%coefficients)
  
  ! Convert state coefficients to a matrix, weighted by the state weights.
  ! This matrix is constructed as the transpose of the eigenvector matrix,
  !    to allow for fast access when calculating matrix elements.
  allocate(coefficients(size(states),no_states), stat=ialloc); call err(ialloc)
  do i=1,size(states)
    coefficients(i,:) = states(i)%coefficients * sqrt(weights(i))
  enddo
  
  output%state_ids = states(1)%state_ids
  selected_couplings = selected_states_couplings(couplings, output%state_ids)
  no_elements = sum([( size(selected_couplings(i)), &
                     & i=1,                         &
                     & size(selected_couplings)     )])
  allocate( output%bra_ids(no_elements), &
          & output%ket_ids(no_elements), &
          & output%values(no_elements),  &
          & stat=ialloc); call err(ialloc)
  k = 0
  do i=1,size(selected_couplings)
    do j=1,size(selected_couplings(i))
      k = k+1
      output%bra_ids(k) = output%state_ids(i)
      output%ket_ids(k) = output%state_ids(selected_couplings(i)%i(j))
      output%values(k) = dot_product(                  &
         & coefficients(:,selected_couplings(i)%i(j)), &
         & coefficients(:,i)                           )
    enddo
  enddo
end procedure
end submodule
