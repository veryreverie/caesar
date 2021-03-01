submodule (caesar_interpolation_module) caesar_interpolation_submodule
  use caesar_anharmonic_module
contains

module procedure interpolate_dynamical_matrices
  type(CartesianHessian) :: coarse_hessian
  
  integer :: i
   
  coarse_hessian = reconstruct_hessian( coarse_supercell,         &
                                      & coarse_qpoints,           &
                                      & coarse_dynamical_matrices )
  
  output = [(                                &
     & DynamicalMatrix( fine_qpoints(i),     &
     &                  coarse_supercell,    &
     &                  coarse_hessian,      &
     &                  coarse_min_images ), &
     & i=1,                                  &
     & size(fine_qpoints)                    )]
end procedure

module procedure calculate_difference_dynamical_matrices
  type(DynamicalMatrix), allocatable :: coarse_dynamical_matrices(:)
  type(DynamicalMatrix), allocatable :: interpolated_dynamical_matrices(:)
  
  ! Interpolate the dynamical matrices onto the coarse grid.
  coarse_dynamical_matrices = interpolate_dynamical_matrices( &
                                   & fine_dynamical_matrices, &
                                   & fine_qpoints,            &
                                   & fine_supercell,          &
                                   & fine_min_images,         &
                                   & coarse_qpoints           )
  
  ! Interpolate the dynamical matrices back onto the fine grid.
  interpolated_dynamical_matrices = interpolate_dynamical_matrices( &
                                       & coarse_dynamical_matrices, &
                                       & coarse_qpoints,            &
                                       & coarse_supercell,          &
                                       & coarse_min_images,         &
                                       & fine_qpoints               )
  
  ! The output is the difference between the original dynamical matrices,
  !    and the contribution to these dynamical matrices which has survived
  !    interpolation.
  output = fine_dynamical_matrices - interpolated_dynamical_matrices
end procedure

module procedure calculate_interpolated_hessian
  type(DynamicalMatrix), allocatable :: dynamical_matrices(:)
  
  
  ! Construct the output: the interpolated anharmonic hessian,
  !    plus the difference between the harmonic hessian and the interpolated
  !    coarse-grained harmonic hessian.
  dynamical_matrices =                                                    &
     &   interpolate_dynamical_matrices( anharmonic_dynamical_matrices,   &
     &                                   anharmonic_qpoints,              &
     &                                   anharmonic_supercell,            &
     &                                   anharmonic_min_images,           &
     &                                   harmonic_qpoints               ) &
     & + calculate_difference_dynamical_matrices(                         &
     &               harmonic_dynamical_matrices,                         &
     &               harmonic_qpoints,                                    &
     &               harmonic_supercell,                                  &
     &               harmonic_min_images,                                 &
     &               anharmonic_qpoints,                                  &
     &               anharmonic_supercell,                                &
     &               anharmonic_min_images        )
  
  output = reconstruct_hessian( harmonic_supercell, &
                              & harmonic_qpoints,   &
                              & dynamical_matrices  )
end procedure

module procedure calculate_interpolated_stress_hessian
  type(DynamicalMatrix),  allocatable :: harmonic_matrices(:)
  type(DynamicalMatrix),  allocatable :: anharmonic_matrices(:)
  type(CartesianHessian), allocatable :: elements(:,:)
  
  integer :: i,j,k,ialloc
  
  allocate(elements(3,3), stat=ialloc); call err(ialloc)
  do i=1,3
    do j=1,3
      harmonic_matrices = [(                             &
         & harmonic_dynamical_matrices(k)%elements(j,i), &
         & k=1,                                          &
         & size(harmonic_dynamical_matrices)             )]
      
      anharmonic_matrices = [(                             &
         & anharmonic_dynamical_matrices(k)%elements(j,i), &
         & k=1,                                            &
         & size(anharmonic_dynamical_matrices)             )]
      
      elements(j,i) = calculate_interpolated_hessian( &
                             & anharmonic_supercell,  &
                             & anharmonic_qpoints,    &
                             & anharmonic_matrices,   &
                             & anharmonic_min_images, &
                             & harmonic_supercell,    &
                             & harmonic_qpoints,      &
                             & harmonic_matrices,     &
                             & harmonic_min_images    )
    enddo
  enddo
  
  output = StressHessian(elements)
end procedure
end submodule
