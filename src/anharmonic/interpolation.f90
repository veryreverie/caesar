! ======================================================================
! Interpolate an anharmonic dynamical matrix to the same q-point density
!    as a given harmonic dynamical matrix.
! ======================================================================
module interpolation_module
  use common_module
  
  use anharmonic_common_module
  
  private
  
  public :: calculate_interpolated_hessian
  public :: calculate_interpolated_stress_hessian
contains

! Takes anharmonic dynamical matrices on a coarse q-point grid,
!    and harmonic dynamical matrices on a fine grid,
!    and returns a fine-grid Hessian containing both.
function calculate_interpolated_hessian(anharmonic_supercell,                &
   & anharmonic_qpoints,anharmonic_dynamical_matrices,anharmonic_min_images, &
   & harmonic_supercell,harmonic_qpoints,harmonic_dynamical_matrices,        &
   & harmonic_min_images) result(output) 
  implicit none
  
  type(StructureData),   intent(in) :: anharmonic_supercell
  type(QpointData),      intent(in) :: anharmonic_qpoints(:)
  type(DynamicalMatrix), intent(in) :: anharmonic_dynamical_matrices(:)
  type(MinImages),       intent(in) :: anharmonic_min_images(:,:)
  type(StructureData),   intent(in) :: harmonic_supercell
  type(QpointData),      intent(in) :: harmonic_qpoints(:)
  type(DynamicalMatrix), intent(in) :: harmonic_dynamical_matrices(:)
  type(MinImages),       intent(in) :: harmonic_min_images(:,:)
  type(CartesianHessian)            :: output
  
  type(CartesianHessian) :: fine_harmonic_hessian
  type(CartesianHessian) :: coarse_harmonic_hessian
  type(CartesianHessian) :: coarse_anharmonic_hessian
  
  type(DynamicalMatrix), allocatable :: coarse_harmonic_dynamical_matrices(:)
  type(DynamicalMatrix), allocatable :: interpolated_harmonic_dynamical_matrices(:)
  type(DynamicalMatrix), allocatable :: interpolated_anharmonic_dynamical_matrices(:)
  type(DynamicalMatrix), allocatable :: output_dynamical_matrices(:)
  
  integer :: i
  
  ! Interpolate the harmonic hessian onto the coarse grid and then back to the
  !    fine grid, to remove the fine-grid-only components.
  fine_harmonic_hessian = reconstruct_hessian( harmonic_supercell,         &
                                             & harmonic_qpoints,           &
                                             & harmonic_dynamical_matrices )
  
  coarse_harmonic_dynamical_matrices = [(         &
     & DynamicalMatrix( anharmonic_qpoints(i),    &
     &                  harmonic_supercell,       &
     &                  fine_harmonic_hessian,    &
     &                  harmonic_min_images    ), &
     & i=1,                                       &
     & size(anharmonic_qpoints)                   )]
  
  coarse_harmonic_hessian = reconstruct_hessian( &
            & anharmonic_supercell,              &
            & anharmonic_qpoints,                &
            & coarse_harmonic_dynamical_matrices )
  
  interpolated_harmonic_dynamical_matrices = [(     &
     & DynamicalMatrix( harmonic_qpoints(i),        &
     &                  anharmonic_supercell,       &
     &                  coarse_harmonic_hessian,    &
     &                  anharmonic_min_images    ), &
     & i=1,                                         &
     & size(harmonic_qpoints)                       )]
  
  ! Interpolate the anharmonic Hessian onto the fine grid.
  coarse_anharmonic_hessian = reconstruct_hessian( &
                   & anharmonic_supercell,         &
                   & anharmonic_qpoints,           &
                   & anharmonic_dynamical_matrices )
  
  interpolated_anharmonic_dynamical_matrices = [(     &
     & DynamicalMatrix( harmonic_qpoints(i),          &
     &                  anharmonic_supercell,         &
     &                  coarse_anharmonic_hessian,    &
     &                  anharmonic_min_images      ), &
     & i=1,                                           &
     & size(harmonic_qpoints)                         )]
  
  ! Construct the output: the interpolated anharmonic hessian,
  !    plus the difference between the harmonic hessian and the interpolated
  !    coarse-grained harmonic hessian.
  output_dynamical_matrices = interpolated_anharmonic_dynamical_matrices &
                          & + harmonic_dynamical_matrices                &
                          & - interpolated_harmonic_dynamical_matrices
  
  output = reconstruct_hessian( harmonic_supercell,       &
                              & harmonic_qpoints,         &
                              & output_dynamical_matrices )
end function

function calculate_interpolated_stress_hessian(anharmonic_supercell,         &
   & anharmonic_qpoints,anharmonic_dynamical_matrices,anharmonic_min_images, &
   & harmonic_supercell,harmonic_qpoints,harmonic_dynamical_matrices,        &
   & harmonic_min_images) result(output) 
  implicit none
  
  type(StructureData),         intent(in) :: anharmonic_supercell
  type(QpointData),            intent(in) :: anharmonic_qpoints(:)
  type(StressDynamicalMatrix), intent(in) :: anharmonic_dynamical_matrices(:)
  type(MinImages),             intent(in) :: anharmonic_min_images(:,:)
  type(StructureData),         intent(in) :: harmonic_supercell
  type(QpointData),            intent(in) :: harmonic_qpoints(:)
  type(StressDynamicalMatrix), intent(in) :: harmonic_dynamical_matrices(:)
  type(MinImages),             intent(in) :: harmonic_min_images(:,:)
  type(StressHessian)                     :: output
  
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
end function
end module
