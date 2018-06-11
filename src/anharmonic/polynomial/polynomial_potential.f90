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
    procedure, public :: generate_potential => &
       & generate_potential_PolynomialPotential
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
subroutine generate_sampling_points_PolynomialPotential(this,inputs, &
   & sampling_points_dir,logfile,write_lambda)
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  type(AnharmonicData),       intent(in)    :: inputs
  type(String),               intent(in)    :: sampling_points_dir
  type(OFile),                intent(inout) :: logfile
  procedure(WriteLambda)                    :: write_lambda
  
  ! Variables used when generating sampling points.
  type(SubspaceMonomial), allocatable :: subspace_monomials(:)
  
  ! Supercell variables.
  type(QpointData), allocatable :: sampling_point_qpoints(:)
  type(IntMatrix)               :: supercell_matrix
  type(StructureData)           :: supercell
  
  ! Displacement data.
  type(RealModeDisplacement)      :: sampling_point
  type(RealModeDisplacement)      :: transformed_sampling_point
  type(VscfRvectors), allocatable :: vscf_rvectors(:)
  type(CartesianDisplacement)     :: displacement
  type(StructureData)             :: displaced_structure
  
  ! Directories and files.
  type(String) :: coupling_dir
  type(OFile)  :: basis_functions_file
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
  allocate( this%basis_functions(size(inputs%subspace_couplings)), &
          & this%sampling_points(size(inputs%subspace_couplings)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(inputs%subspace_couplings)
    ! Generate the set of subspace monomials corresponding to the subspace
    !    coupling.
    ! e.g. the coupling [1,2] might have monomials [1,2], [1,1,2] and [1,2,2].
    subspace_monomials = generate_subspace_monomials( &
                      & inputs%subspace_couplings(i), &
                      & inputs%degenerate_subspaces,  &
                      & this%potential_expansion_order)
    
    ! Generate basis functions at each coupling.
    this%basis_functions(i) = generate_basis_functions( &
                    & subspace_monomials,               &
                    & inputs%structure,                 &
                    & inputs%complex_modes,             &
                    & inputs%real_modes,                &
                    & inputs%qpoints,                   &
                    & inputs%degenerate_subspaces,      &
                    & inputs%degenerate_symmetries,     &
                    & inputs%vscf_basis_functions_only, &
                    & logfile)
    
    ! Generate a set of sampling points from which the basis functions can
    !    be constructed.
    this%sampling_points(i) = generate_sampling_points( &
                & this%basis_functions(i)%functions,    &
                & this%potential_expansion_order,       &
                & inputs%maximum_weighted_displacement, &
                & inputs%frequency_of_max_displacement, &
                & inputs%real_modes)
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
    basis_functions_file = OFile(coupling_dir//'/basis_functions.dat')
    call basis_functions_file%print_lines(this%basis_functions(i))
    
    ! Write sampling points to file.
    sampling_points_file = OFile(coupling_dir//'/sampling_points.dat')
    call sampling_points_file%print_lines(this%sampling_points(i))
    
    do j=1,size(this%sampling_points(i))
      ! Select the sampling point for clarity.
      sampling_point = this%sampling_points(i)%points(j)
      
      ! Make a directory for each sampling point.
      sampling_dir = coupling_dir// &
         & '/sampling_point_'//left_pad(j,str(size(this%sampling_points(i))))
      call mkdir(sampling_dir)
      
      ! Construct a supercell for each sampling point.
      sampling_point_qpoints = sampling_point%qpoints( inputs%real_modes, &
                                                     & inputs%qpoints)
      supercell_matrix = construct_supercell_matrix( sampling_point_qpoints, &
                                                   & inputs%structure)
      supercell = construct_supercell( inputs%structure,          &
                                     & supercell_matrix,          &
                                     & inputs%symmetry_precision, &
                                     & calculate_symmetry=.false.)
      
      ! Write out the supercell.
      call write_structure_file(supercell, sampling_dir//'/structure.dat')
      
      ! Construct VSCF R-vectors.
      vscf_rvectors = construct_vscf_rvectors( sampling_point,    &
                                             & supercell,         &
                                             & inputs%real_modes, &
                                             & inputs%qpoints)
      vscf_rvectors_file = OFile(sampling_dir//'/vscf_rvectors.dat')
      call vscf_rvectors_file%print_lines(vscf_rvectors,separating_line='')
      
      do k=1,size(vscf_rvectors)
        ! Transform the sampling point by the VSCF R-vector.
        transformed_sampling_point = vscf_rvectors(k)%transform( &
                                            & sampling_point,    &
                                            & inputs%real_modes, &
                                            & inputs%qpoints)
        
        ! Construct displaced structure.
        displacement = transformed_sampling_point%cartesian_displacement( &
                                                     & supercell,         &
                                                     & inputs%real_modes, &
                                                     & inputs%qpoints)
        displaced_structure = displace_structure(supercell,displacement)
        
        ! Create directory and structure files for displaced structure.
        vscf_rvectors_dir = sampling_dir// &
           & '/vscf_rvectors_'//left_pad(k,str(size(vscf_rvectors)))
        call write_lambda(displaced_structure,vscf_rvectors_dir)
      enddo
    enddo
  enddo
end subroutine

! Generate potential.
subroutine generate_potential_PolynomialPotential(this,inputs, &
   & sampling_points_dir,logfile,read_lambda)
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  type(AnharmonicData),       intent(in)    :: inputs
  type(String),               intent(in)    :: sampling_points_dir
  type(OFile),                intent(inout) :: logfile
  procedure(ReadLambda)                     :: read_lambda
  
  type(StructureData)             :: supercell
  type(VscfRvectors), allocatable :: vscf_rvectors(:)
  type(ElectronicStructure)       :: electronic_structure
  
  type(ElectronicStructure), allocatable :: calculations(:)
  type(ElectronicStructure), allocatable :: sampled_points(:)
  
  real(dp)                         :: energy
  type(RealModeForce), allocatable :: forces(:)
  type(RealModeForce)              :: force
  
  ! Directories and files.
  type(String) :: coupling_dir
  type(IFile)  :: basis_functions_file
  type(IFile)  :: sampling_points_file
  type(String) :: sampling_dir
  type(IFile)  :: vscf_rvectors_file
  type(String) :: vscf_rvectors_dir
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  ! Read in basis functions and sampling points.
  allocate( this%basis_functions(size(inputs%subspace_couplings)), &
          & this%sampling_points(size(inputs%subspace_couplings)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this%sampling_points)
    coupling_dir = sampling_points_dir// &
       & '/coupling_'//left_pad(i,str(size(this%sampling_points)))
    
    ! Read in basis functions and sampling points.
    basis_functions_file = IFile(coupling_dir//'/basis_functions.dat')
    this%basis_functions(i) = basis_functions_file%lines()
    
    sampling_points_file = IFile(coupling_dir//'/sampling_points.dat')
    this%sampling_points(i) = sampling_points_file%lines()
    
    allocate( sampled_points(size(this%sampling_points(i))), &
            & stat=ialloc); call err(ialloc)
    do j=1,size(this%sampling_points(i))
      sampling_dir = coupling_dir// &
         & '/sampling_point_'//left_pad(j,str(size(this%sampling_points(i))))
      
      ! Read in supercell and VSCF R-vectors.
      supercell = read_structure_file( sampling_dir//'/structure.dat', &
                                     & inputs%symmetry_precision,      &
                                     & calculate_symmetry=.false.)
      
      vscf_rvectors_file = IFile(sampling_dir//'/vscf_rvectors.dat')
      vscf_rvectors = VscfRvectors(vscf_rvectors_file%sections())
      
      allocate( calculations(size(vscf_rvectors)), &
              & stat=ialloc); call err(ialloc)
      do k=1,size(vscf_rvectors)
        vscf_rvectors_dir = sampling_dir// &
           & '/vscf_rvectors_'//left_pad(k,str(size(vscf_rvectors)))
        calculations(k) = read_lambda(vscf_rvectors_dir)
      enddo
      
      ! Average the energy over all VSCF R-vectors, and normalise to be per
      !    primitive cell.
      energy = sum(calculations%energy) &
           & / (size(calculations) * supercell%sc_size)
      
      ! Transform forces into normal mode co-ordinates, reverse the
      !    VSCF R-vector transformation, and find the average.
      allocate(forces(size(vscf_rvectors)), stat=ialloc); call err(ialloc)
      do k=1,size(vscf_rvectors)
        forces(k) = RealModeForce( calculations(k)%forces, &
                                 & supercell,              &
                                 & inputs%real_modes,      &
                                 & inputs%qpoints)
        forces(k) = vscf_rvectors(k)%inverse_transform( forces(k),         &
                                                      & inputs%real_modes, &
                                                      & inputs%qpoints)
      enddo
      force = sum(forces) / real(size(calculations),dp)
      
      ! sampled_points(j) = average(calculations,vscf_rvectors)
      deallocate(calculations, stat=ialloc); call err(ialloc)
    enddo
    deallocate(sampled_points, stat=ialloc); call err(ialloc)
  enddo
  
  ! TODO
end subroutine

!! Average over VSCF R-vectors.
!function average(calculations,vscf_rvectors) result(output)
!  implicit none
!  
!  type(ElectronicStructure), intent(in) :: calculations(:)
!  type(VscfRvectors),        intent(in) :: vscf_rvectors(:)
!  type(ElectronicStructure)             :: output
!  
!  type(ElectronicStructure), allocatable :: untransformed_calculations(:)
!  
!  real(dp)             :: energy
!  type(CartesianForce) :: forces(:)
!  
!  integer :: i,ialloc
!  
!  if (size(calculations)/=size(vscf_rvectors)) then
!    call err()
!  endif
!  
!  allocate( untransformed_calculations(size(calculations)), &
!          & stat=ialloc); call err(ialloc)
!  do i=1,size(calculations)
!    untransformed_calculations = vscf_rvectors(i)%inverse_transform( &
!  enddo
!end function
end module
