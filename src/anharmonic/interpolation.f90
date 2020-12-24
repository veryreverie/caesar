! ======================================================================
! Interpolate an anharmonic dynamical matrix to the same q-point density
!    as a given harmonic dynamical matrix.
! ======================================================================
module interpolation_module
  use common_module
  
  use anharmonic_common_module
  
  private
  
  public :: interpolate_dynamical_matrices
  public :: calculate_difference_dynamical_matrices
  public :: calculate_interpolated_hessian
  public :: calculate_interpolated_stress_hessian
contains

! Takes a set of dynamical matrices defined on a coarse q-point grid,
!    and interpolates to give the dynamical matrices on a fine q-point grid.
function interpolate_dynamical_matrices(coarse_dynamical_matrices,   &
   & coarse_qpoints,coarse_supercell,coarse_min_images,fine_qpoints) &
   & result(output) 
   
  type(DynamicalMatrix), intent(in)  :: coarse_dynamical_matrices(:)
  type(QpointData),      intent(in)  :: coarse_qpoints(:)
  type(StructureData),   intent(in)  :: coarse_supercell
  type(MinImages),       intent(in)  :: coarse_min_images(:,:)
  type(QpointData),      intent(in)  :: fine_qpoints(:)
  type(DynamicalMatrix), allocatable :: output(:)
  
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
end function

! When dynamical matrices are interpolated from a fine q-point grid
!    onto a coarse q-point grid, some information may be lost.
! This function takes two inputs:
!    - a set of dynamical matrices defined on a fine q-point grid.
!    - a coarse q-point grid.
! The function separates the dynamical matrices into the sum of two parts:
!    - the part which is preserved when the dynamical matrices are interpolated
!         onto the coarse q-point grid and back again.
!    - the part which is lost when the dynamical matrices are interpolated
!         onto the coarse q-point grid and back again.
! The function then returns only the contribution to the dynamical matrices
!    which is lost under interpolation.
function calculate_difference_dynamical_matrices(fine_dynamical_matrices, &
   & fine_qpoints,fine_supercell,fine_min_images,coarse_qpoints,          &
   & coarse_supercell,coarse_min_images) result(output) 
  implicit none
  
  type(DynamicalMatrix), intent(in)  :: fine_dynamical_matrices(:)
  type(QpointData),      intent(in)  :: fine_qpoints(:)
  type(StructureData),   intent(in)  :: fine_supercell
  type(MinImages),       intent(in)  :: fine_min_images(:,:)
  type(QpointData),      intent(in)  :: coarse_qpoints(:)
  type(StructureData),   intent(in)  :: coarse_supercell
  type(MinImages),       intent(in)  :: coarse_min_images(:,:)
  type(DynamicalMatrix), allocatable :: output(:)
  
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
end function

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
