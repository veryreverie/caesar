submodule (caesar_harmonic_stress_conversion_module) caesar_harmonic_stress_conversion_submodule
  use caesar_harmonic_stress_module
contains

module procedure new_StressDynamicalMatrix_interpolated
  type(DynamicalMatrix), allocatable :: elements(:,:)
  
  integer :: i,j,ialloc
  
  allocate(elements(3,3), stat=ialloc); call err(ialloc)
  do i=1,3
    do j=1,3
      elements(j,i) = DynamicalMatrix( q,                     &
                                     & supercell,             &
                                     & hessian%elements(j,i), &
                                     & min_images             )
    enddo
  enddo
  
  output = StressDynamicalMatrix(elements)
end procedure

module procedure reconstruct_stress_hessian
  type(CartesianHessian), allocatable :: elements(:,:)
  
  integer :: i,j,k,ialloc
  
  allocate(elements(3,3), stat=ialloc); call err(ialloc)
  do i=1,3
    do j=1,3
      elements(j,i) = reconstruct_hessian(                              &
         & large_supercell,                                             &
         & qpoints,                                                     &
         & [(dynamical_matrices(k)%elements(j,i), k=1, size(qpoints))], &
         & logfile                                                      )
    enddo
  enddo
  
  output = StressHessian(elements)
end procedure
end submodule
