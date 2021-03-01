submodule (caesar_interpolated_supercell_module) caesar_interpolated_supercell_submodule
  use caesar_anharmonic_common_module
contains

module procedure new_InterpolatedSupercell
  this%supercell          = supercell
  this%qpoints            = qpoints
  this%dynamical_matrices = dynamical_matrices
  this%complex_modes      = complex_modes
  this%real_modes         = real_modes
end procedure

module procedure new_InterpolatedSupercell_interpolated
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
end procedure
end submodule
