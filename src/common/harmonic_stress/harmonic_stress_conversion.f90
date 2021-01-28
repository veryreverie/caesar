! ======================================================================
! Conversions between cartesian and normal-mode representations of
!    harmonic stress.
! ======================================================================
module caesar_harmonic_stress_conversion_module
  use caesar_utils_module
  use caesar_structure_module
  use caesar_normal_mode_module
  use caesar_dynamical_matrices_module
  
  use caesar_stress_hessian_module
  use caesar_stress_dynamical_matrix_module
  implicit none
  
  private
  
  public :: StressDynamicalMatrix
  public :: reconstruct_stress_hessian
  
  interface StressDynamicalMatrix
    module procedure new_StressDynamicalMatrix_interpolated
  end interface
contains

! Construct a stress dynamical matrix at an arbitrary q-point.
function new_StressDynamicalMatrix_interpolated(q,supercell,hessian, &
   & min_images) result(output) 
  implicit none
  
  type(RealVector),    intent(in)           :: q
  type(StructureData), intent(in)           :: supercell
  type(StressHessian), intent(in)           :: hessian
  type(MinImages),     intent(in), optional :: min_images(:,:)
  type(StressDynamicalMatrix)               :: output
  
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
end function

function reconstruct_stress_hessian(large_supercell,qpoints, &
   & dynamical_matrices,logfile) result(output) 
  implicit none
  
  type(StructureData),         intent(in)    :: large_supercell
  type(QpointData),            intent(in)    :: qpoints(:)
  type(StressDynamicalMatrix), intent(in)    :: dynamical_matrices(:)
  type(OFile),                 intent(inout) :: logfile
  type(StressHessian)                        :: output
  
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
end function
end module
