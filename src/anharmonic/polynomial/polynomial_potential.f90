! ======================================================================
! A polynomial representation of a potential.
! ======================================================================
module polynomial_potential_module
  use common_module
  
  use anharmonic_common_module
  
  use basis_function_module
  use basis_functions_module
  use vscf_rvectors_module
  use sampling_points_module
  implicit none
  
  private
  
  public :: PolynomialPotential
  
  type, extends(PotentialData) :: PolynomialPotential
    integer                           :: potential_expansion_order
    type(BasisFunctions), allocatable :: basis_functions(:)
    type(SamplingPoints), allocatable :: sampling_points(:)
  contains
    procedure, public :: generate_sampling_points => &
       & generate_sampling_points_PolynomialPotential
  end type
  
  interface PolynomialPotential
    module procedure new_PolynomialPotential
  end interface
contains

! Constructor.
function new_PolynomialPotential(potential_expansion_order) result(this)
  implicit none
  
  integer, intent(in)       :: potential_expansion_order
  type(PolynomialPotential) :: this
  
  this%potential_expansion_order = potential_expansion_order
end function

! Generate sampling points.
subroutine generate_sampling_points_PolynomialPotential(this,             &
   & sampling_points_dir,structure,symmetry_precision,complex_modes,      &
   & real_modes,qpoints,degenerate_subspaces,degenerate_symmetries,       &
   & coupled_subspaces,vscf_basis_functions_only,                         &
   & maximum_weighted_displacement,frequency_of_max_displacement,logfile, &
   & write_lambda)
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  type(String),               intent(in)    :: sampling_points_dir
  type(StructureData),        intent(in)    :: structure
  real(dp),                   intent(in)    :: symmetry_precision
  type(ComplexMode),          intent(in)    :: complex_modes(:)
  type(RealMode),             intent(in)    :: real_modes(:)
  type(QpointData),           intent(in)    :: qpoints(:)
  type(DegenerateModes),      intent(in)    :: degenerate_subspaces(:)
  type(DegenerateSymmetry),   intent(in)    :: degenerate_symmetries(:)
  type(CoupledSubspaces),     intent(in)    :: coupled_subspaces(:)
  logical,                    intent(in)    :: vscf_basis_functions_only
  real(dp),                   intent(in)    :: maximum_weighted_displacement
  real(dp),                   intent(in)    :: frequency_of_max_displacement
  type(OFile),                intent(inout) :: logfile
  procedure(WriteLambda)                    :: write_lambda
  
  ! Variables used when generating sampling points.
  type(SubspaceMonomial), allocatable :: subspace_monomials(:)
  
  ! Supercell variables.
  type(QpointData), allocatable :: sampling_point_qpoints(:)
  type(IntMatrix)               :: supercell_matrix
  type(StructureData)           :: supercell
  
  ! Displacement data.
  type(VscfRvectors), allocatable :: vscf_rvectors(:)
  type(CartesianDisplacement)     :: displacement
  type(StructureData)             :: displaced_structure
  
  ! Directories and files.
  type(String) :: coupling_dir
  type(OFile)  :: basis_function_file
  type(OFile)  :: sampling_points_file
  type(String) :: sampling_dir
  type(OFile)  :: vscf_rvectors_file
  type(String) :: vscf_rvectors_dir
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Generate basis functions and sampling points.
  ! --------------------------------------------------
  
  ! Loop over subspace couplings, generating basis functions and sampling
  !    points for each.
  allocate( this%basis_functions(size(coupled_subspaces)), &
          & this%sampling_points(size(coupled_subspaces)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(coupled_subspaces)
    ! Generate the set of subspace monomials corresponding to the subspace
    !    coupling.
    ! e.g. the coupling [1,2] might have monomials [1,2], [1,1,2] and [1,2,2].
    subspace_monomials = generate_subspace_monomials( &
                              & coupled_subspaces(i), &
                              & this%potential_expansion_order)
    
    
    ! Generate basis functions at each coupling.
    this%basis_functions(i) = generate_basis_functions( &
                           & subspace_monomials,        &
                           & structure,                 &
                           & complex_modes,             &
                           & real_modes,                &
                           & qpoints,                   &
                           & degenerate_subspaces,      &
                           & degenerate_symmetries,     &
                           & vscf_basis_functions_only, &
                           & logfile)
    
    ! Generate a set of sampling points from which the basis functions can
    !    be constructed.
    this%sampling_points(i) = generate_sampling_points( &
                   & this%basis_functions(i)%functions, &
                   & this%potential_expansion_order,    &
                   & maximum_weighted_displacement,     &
                   & frequency_of_max_displacement,     &
                   & real_modes)
  enddo
  
  ! --------------------------------------------------
  ! Write out basis functions and sampling points.
  ! --------------------------------------------------
  do i=1,size(this%sampling_points)
    ! Make a directory for each coupling.
    coupling_dir = sampling_points_dir// &
                 & '/coupling_'//left_pad(i,str(size(this%sampling_points)))
    call mkdir(coupling_dir)
    
    ! Write basis functions to file.
    basis_function_file = OFile(coupling_dir//'/basis_functions.dat')
    call basis_function_file%print_lines(this%basis_functions(i))
    
    ! Write sampling points to file.
    sampling_points_file = OFile(coupling_dir//'/sampling_points.dat')
    call sampling_points_file%print_lines(this%sampling_points(i))
    
    do j=1,size(this%sampling_points(i))
      ! Make a directory for each sampling point.
      sampling_dir = coupling_dir// &
         & '/sampling_point_'//left_pad(j,str(size(this%sampling_points(i))))
      call mkdir(sampling_dir)
      
      ! Construct a supercell for each sampling point.
      sampling_point_qpoints = &
         & this%sampling_points(i)%points(j)%qpoints(real_modes,qpoints)
      supercell_matrix = construct_supercell_matrix( sampling_point_qpoints, &
                                                   & structure)
      supercell = construct_supercell( structure,          &
                                     & supercell_matrix,   &
                                     & symmetry_precision, &
                                     & calculate_symmetry=.false.)
      
      ! Write out the supercell.
      call write_structure_file(supercell, sampling_dir//'/structure.dat')
      
      ! Construct VSCF R-vectors.
      vscf_rvectors = construct_vscf_rvectors( &
          & this%sampling_points(i)%points(j), &
          & supercell,                         &
          & real_modes,                        &
          & qpoints)
      vscf_rvectors_file = OFile(sampling_dir//'/vscf_rvectors.dat')
      call vscf_rvectors_file%print_lines(vscf_rvectors,separating_line='')
      
      do k=1,size(vscf_rvectors)
        ! Construct displaced structure.
        displacement =                                                 &
           & this%sampling_points(i)%points(j)%cartesian_displacement( &
           &                             supercell,                    &
           &                             real_modes,                   &
           &                             qpoints,                      &
           &                             vscf_rvectors(k)%rvectors(real_modes))
        displaced_structure = displace_structure(supercell,displacement)
        
        ! Create directory and structure files for displaced structure.
        vscf_rvectors_dir = sampling_dir// &
           & '/vscf_rvectors_'//left_pad(k,str(size(vscf_rvectors)))
        call write_lambda(displaced_structure,vscf_rvectors_dir)
      enddo
    enddo
  enddo
end subroutine

end module
