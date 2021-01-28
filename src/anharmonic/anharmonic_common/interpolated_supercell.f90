! ======================================================================
! Generates a supercell with a given q-point grid,
!    along with the corresponding q-points and normal modes.
! ======================================================================
module caesar_interpolated_supercell_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: InterpolatedSupercell
  
  type, extends(NoDefaultConstructor) :: InterpolatedSupercell
    type(StructureData)                :: supercell
    type(QpointData),      allocatable :: qpoints(:)
    type(DynamicalMatrix), allocatable :: dynamical_matrices(:)
    type(ComplexMode),     allocatable :: complex_modes(:)
    type(RealMode),        allocatable :: real_modes(:)
  end type
  
  interface InterpolatedSupercell
    module procedure new_InterpolatedSupercell
    module procedure new_InterpolatedSupercell_interpolated
  end interface
contains

function new_InterpolatedSupercell(supercell,qpoints,dynamical_matrices, &
   & complex_modes,real_modes) result(this) 
  implicit none
  
  type(StructureData),   intent(in) :: supercell
  type(QpointData),      intent(in) :: qpoints(:)
  type(DynamicalMatrix), intent(in) :: dynamical_matrices(:)
  type(ComplexMode),     intent(in) :: complex_modes(:)
  type(RealMode),        intent(in) :: real_modes(:)
  type(InterpolatedSupercell)       :: this
  
  this%supercell          = supercell
  this%qpoints            = qpoints
  this%dynamical_matrices = dynamical_matrices
  this%complex_modes      = complex_modes
  this%real_modes         = real_modes
end function

function new_InterpolatedSupercell_interpolated(qpoint_grid,structure, &
   & harmonic_supercell,harmonic_qpoints,harmonic_dynamical_matrices,  &
   & harmonic_complex_modes,logfile) result(this) 
  implicit none
  
  integer,               intent(in)    :: qpoint_grid(:)
  type(StructureData),   intent(in)    :: structure
  type(StructureData),   intent(in)    :: harmonic_supercell
  type(QpointData),      intent(in)    :: harmonic_qpoints(:)
  type(DynamicalMatrix), intent(in)    :: harmonic_dynamical_matrices(:)
  type(ComplexMode),     intent(in)    :: harmonic_complex_modes(:,:)
  type(OFile),           intent(inout) :: logfile
  type(InterpolatedSupercell)          :: this
  
  type(IntMatrix)                    :: supercell_matrix
  type(StructureData)                :: supercell
  type(QpointData),      allocatable :: qpoints(:)
  type(CartesianHessian)             :: hessian
  type(MinImages),       allocatable :: min_images(:,:)
  type(DynamicalMatrix), allocatable :: dynamical_matrices(:)
  type(MatrixAndModes),  allocatable :: matrices_and_modes(:)
  type(ComplexMode),     allocatable :: complex_modes(:)
  type(RealMode),        allocatable :: real_modes(:)
  
  integer :: i,ialloc
  
  ! Calculate the supercell matrix, supercell and q-points
  !    for the anharmonic supercell.
  supercell_matrix =                                           &
     & mat([ qpoint_grid(1), 0             , 0            ,    &
     &       0             , qpoint_grid(2), 0            ,    &
     &       0             , 0             , qpoint_grid(3) ], &
     & 3,3)
  supercell = construct_supercell( structure,       &
                                 & supercell_matrix )
  qpoints = generate_qpoints(supercell)
  
  ! Construct the Hessian and minimum images of the harmonic supercell.
  hessian = reconstruct_hessian( harmonic_supercell,          &
                               & harmonic_qpoints,            &
                               & harmonic_dynamical_matrices, &
                               & logfile                      )
  
  min_images = calculate_min_images(harmonic_supercell)
  
  ! Construct the dynamical matrices and Hessian corresponding to
  !    the anharmonic q-points.
  allocate( dynamical_matrices(size(qpoints)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(qpoints)
    dynamical_matrices(i) = DynamicalMatrix( dblevec(qpoints(i)%qpoint), &
                                           & harmonic_supercell,         &
                                           & hessian,                    &
                                           & min_images                  )
  enddo
  
  hessian = reconstruct_hessian( supercell, &
                               & qpoints, &
                               & dynamical_matrices, &
                               & logfile)
  
  ! Re-calculate the dynamical matrices and calculate the normal modes
  !    corresponding to the anahrmonic q-points, using the symmetry relations.
  matrices_and_modes = calculate_dynamical_matrices( structure,   &
                                                   & [supercell], &
                                                   & [hessian],   &
                                                   & qpoints,     &
                                                   & logfile      )
  
  dynamical_matrices = [( matrices_and_modes(i)%matrix, &
                        & i=1,                          &
                        & size(matrices_and_modes)      )]
  
  allocate(complex_modes(0), stat=ialloc); call err(ialloc)
  do i=1,size(qpoints)
    complex_modes = [complex_modes, matrices_and_modes(i)%modes]
  enddo
  complex_modes = complex_modes(filter(.not.complex_modes%translational_mode))
  
  real_modes = complex_to_real(complex_modes)
  
  ! Construct the output.
  this = InterpolatedSupercell( supercell,          &
                              & qpoints,            &
                              & dynamical_matrices, &
                              & complex_modes,      &
                              & real_modes          )
end function
end module
