! ======================================================================
! Generates a supercell with a given q-point grid,
!    along with the corresponding q-points and normal modes.
! ======================================================================
module interpolated_supercell_module
  use common_module
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

function new_InterpolatedSupercell_interpolated(qpoint_grid,structure,    &
   & harmonic_qpoints,harmonic_dynamical_matrices,harmonic_complex_modes) &
   & result(this) 
  implicit none
  
  integer,               intent(in) :: qpoint_grid(:)
  type(StructureData),   intent(in) :: structure
  type(QpointData),      intent(in) :: harmonic_qpoints(:)
  type(DynamicalMatrix), intent(in) :: harmonic_dynamical_matrices(:)
  type(ComplexMode),     intent(in) :: harmonic_complex_modes(:,:)
  type(InterpolatedSupercell)       :: this
  
  type(IntMatrix) :: supercell_matrix
  
  type(StructureData)                :: supercell
  type(QpointData),      allocatable :: qpoints(:)
  type(DynamicalMatrix), allocatable :: dynamical_matrices(:)
  type(ComplexMode),     allocatable :: complex_modes(:)
  type(RealMode),        allocatable :: real_modes(:)
  
  integer :: i,j,ialloc
  
  supercell_matrix =                                           &
     & mat([ qpoint_grid(1), 0             , 0            ,    &
     &       0             , qpoint_grid(2), 0            ,    &
     &       0             , 0             , qpoint_grid(3) ], &
     & 3,3)
  supercell = construct_supercell( structure,       &
                                 & supercell_matrix )
  qpoints = generate_qpoints(supercell)
  
  allocate( dynamical_matrices(size(qpoints)), &
          & complex_modes(0),                  &
          & real_modes(0),                     &
          & stat=ialloc); call err(ialloc)
  do i=1,size(qpoints)
    ! TODO: interpolate rather than finding.
    j = first(harmonic_qpoints==qpoints(i), default=0)
    if (j==0) then
      call print_line(ERROR//': interpolated q-point '//qpoints(i)%qpoint// &
         &' is not also a harmonic q-point.')
      call quit()
    endif
    
    dynamical_matrices(i) = harmonic_dynamical_matrices(j)
    complex_modes = [complex_modes, harmonic_complex_modes(:,j)]
  enddo
  
  complex_modes = complex_modes(filter(.not.complex_modes%translational_mode))
  real_modes = complex_to_real(complex_modes)
  
  this = InterpolatedSupercell( supercell,          &
                              & qpoints,            &
                              & dynamical_matrices, &
                              & complex_modes,      &
                              & real_modes          )
end function
end module
