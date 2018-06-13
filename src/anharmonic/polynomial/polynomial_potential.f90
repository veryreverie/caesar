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
  use sample_result_module
  use sample_results_module
  use fit_coefficients_module
  implicit none
  
  private
  
  public :: PolynomialPotential
  
  type, extends(PotentialData) :: PolynomialPotential
    integer,                          private :: potential_expansion_order
    real(dp),                         private :: energy_baseline
    type(BasisFunction), allocatable, private :: basis_functions(:)
    real(dp),            allocatable, private :: coefficients(:)
  contains
    procedure, public :: generate_sampling_points => &
       & generate_sampling_points_PolynomialPotential
    procedure, public :: generate_potential => &
       & generate_potential_PolynomialPotential
    
    procedure, public :: energy => energy_PolynomialPotential
    procedure, public :: force  => force_PolynomialPotential
    
    procedure, public :: read  => read_PolynomialPotential
    procedure, public :: write => write_PolynomialPotential
  end type
  
  interface PolynomialPotential
    module procedure new_PolynomialPotential
    module procedure new_PolynomialPotential_StringArray
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
  type(BasisFunctions),   allocatable :: basis_functions(:)
  type(SamplingPoints),   allocatable :: sampling_points(:)
  
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
  allocate( basis_functions(size(inputs%subspace_couplings)), &
          & sampling_points(size(inputs%subspace_couplings)), &
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
    basis_functions(i) = generate_basis_functions( &
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
    sampling_points(i) = generate_sampling_points( &
           & basis_functions(i)%functions,         &
           & this%potential_expansion_order,       &
           & inputs%maximum_weighted_displacement, &
           & inputs%frequency_of_max_displacement, &
           & inputs%real_modes)
  enddo
  
  ! --------------------------------------------------
  ! Write out basis functions and sampling points.
  ! --------------------------------------------------
  do i=1,size(sampling_points)
    ! Make a directory for each coupling.
    coupling_dir = sampling_points_dir// &
                 & '/coupling_'//left_pad(i,str(size(sampling_points)))
    call mkdir(coupling_dir)
    
    ! Write basis functions to file.
    basis_functions_file = OFile(coupling_dir//'/basis_functions.dat')
    call basis_functions_file%print_lines(basis_functions(i))
    
    ! Write sampling points to file.
    sampling_points_file = OFile(coupling_dir//'/sampling_points.dat')
    call sampling_points_file%print_lines(sampling_points(i))
    
    do j=1,size(sampling_points(i))
      ! Select the sampling point for clarity.
      sampling_point = sampling_points(i)%points(j)
      
      ! Make a directory for each sampling point.
      sampling_dir = coupling_dir// &
         & '/sampling_point_'//left_pad(j,str(size(sampling_points(i))))
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
   & weighted_energy_force_ratio,sampling_points_dir,logfile,read_lambda)
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  type(AnharmonicData),       intent(in)    :: inputs
  real(dp),                   intent(in)    :: weighted_energy_force_ratio
  type(String),               intent(in)    :: sampling_points_dir
  type(OFile),                intent(inout) :: logfile
  procedure(ReadLambda)                     :: read_lambda
  
  ! Variables for processing electronic structure.
  type(StructureData)                    :: supercell
  type(VscfRvectors),        allocatable :: vscf_rvectors(:)
  type(ElectronicStructure), allocatable :: calculations(:)
  
  ! Previously calculated basis functions.
  type(BasisFunctions), allocatable :: basis_functions(:)
  
  ! Electronic structure results.
  type(SamplingPoints), allocatable :: sampling_points(:)
  type(SampleResults),  allocatable :: sample_results(:)
  
  ! Variables for generating coefficients.
  logical,                    allocatable :: uncoupled(:)
  type(BasisFunction),        allocatable :: uncoupled_basis_functions(:)
  type(RealModeDisplacement), allocatable :: uncoupled_sampling_points(:)
  type(SampleResult),         allocatable :: uncoupled_sample_results(:)
  real(dp),                   allocatable :: coefficients(:)
  
  ! Directories and files.
  type(String) :: coupling_dir
  type(IFile)  :: basis_functions_file
  type(IFile)  :: sampling_points_file
  type(String) :: sampling_dir
  type(IFile)  :: vscf_rvectors_file
  type(String) :: vscf_rvectors_dir
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Read in basis functions and sampling points.
  ! --------------------------------------------------
  allocate( basis_functions(size(inputs%subspace_couplings)), &
          & sampling_points(size(inputs%subspace_couplings)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(sampling_points)
    coupling_dir = sampling_points_dir// &
       & '/coupling_'//left_pad(i,str(size(sampling_points)))
    
    ! Read in basis functions and sampling points.
    basis_functions_file = IFile(coupling_dir//'/basis_functions.dat')
    basis_functions(i) = basis_functions_file%lines()
    
    sampling_points_file = IFile(coupling_dir//'/sampling_points.dat')
    sampling_points(i) = sampling_points_file%lines()
  enddo
  
  ! --------------------------------------------------
  ! Read in energies and forces.
  ! --------------------------------------------------
  allocate( sample_results(size(inputs%subspace_couplings)),        &
          & stat=ialloc); call err(ialloc)
  do i=1,size(sampling_points)
    coupling_dir = sampling_points_dir// &
       & '/coupling_'//left_pad(i,str(size(sampling_points)))
    allocate( sample_results(i)%results(size(sampling_points(i))), &
            & stat=ialloc); call err(ialloc)
    do j=1,size(sampling_points(i))
      sampling_dir = coupling_dir// &
         & '/sampling_point_'//left_pad(j,str(size(sampling_points(i))))
      
      ! Read in supercell and VSCF R-vectors.
      supercell = read_structure_file( sampling_dir//'/structure.dat', &
                                     & inputs%symmetry_precision,      &
                                     & calculate_symmetry=.false.)
      
      vscf_rvectors_file = IFile(sampling_dir//'/vscf_rvectors.dat')
      vscf_rvectors = VscfRvectors(vscf_rvectors_file%sections())
      
      ! Read in electronic structure calculations.
      allocate( calculations(size(vscf_rvectors)), &
              & stat=ialloc); call err(ialloc)
      do k=1,size(vscf_rvectors)
        vscf_rvectors_dir = sampling_dir// &
           & '/vscf_rvectors_'//left_pad(k,str(size(vscf_rvectors)))
        calculations(k) = read_lambda(vscf_rvectors_dir)
      enddo
      
      ! Average electronic structure across VSCF R-vectors, and convert
      !    to correct normalisation and real mode co-ordinates.
      sample_results(i)%results(j) = SampleResult( vscf_rvectors,     &
                                                 & calculations,      &
                                                 & supercell,         &
                                                 & inputs%real_modes, &
                                                 & inputs%qpoints)
      
      deallocate(calculations, stat=ialloc); call err(ialloc)
    enddo
  enddo
  
  ! --------------------------------------------------
  ! Calculate basis function coefficients.
  ! --------------------------------------------------
  
  ! Identify which subspace couplings only contain one subspace.
  allocate( uncoupled(size(inputs%subspace_couplings)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(inputs%subspace_couplings)
    if (size(inputs%subspace_couplings(i))==1) then
      uncoupled(i) = .true.
    else
      uncoupled(i) = .false.
    endif
  enddo
  
  ! Calculate the coefficients of all basis functions which do not involve
  !    subspace coupling at once. This takes O(no_modes^3) time.
  uncoupled_basis_functions = [BasisFunction::]
  uncoupled_sampling_points = [RealModeDisplacement::]
  uncoupled_sample_results  = [SampleResult::]
  do i=1,size(inputs%subspace_couplings)
    if (uncoupled(i)) then
      uncoupled_basis_functions = [ uncoupled_basis_functions, &
                                  & basis_functions(i)%functions ]
      uncoupled_sampling_points = [ uncoupled_sampling_points, &
                                  & sampling_points(i)%points ]
      uncoupled_sample_results = [ uncoupled_sample_results, &
                                 & sample_results(i)%results ]
    endif
  enddo
  
  coefficients = fit_coefficients( uncoupled_basis_functions, &
                                 & uncoupled_sampling_points, &
                                 & uncoupled_sample_results,  &
                                 & inputs%real_modes,         &
                                 & weighted_energy_force_ratio)
  
  this%basis_functions = uncoupled_basis_functions
  this%coefficients    = coefficients
  
  ! TODO: calculate energy baseline
  
  ! Calculate the coefficients of all basis functions involving subspace
  !    coupling. These are calculated on a coupling-by-coupling basis.
  do i=1,size(inputs%subspace_couplings)
    if (uncoupled(i)) then
      cycle
    endif
    
    coefficients = fit_coefficients( basis_functions(i)%functions, &
                                   & sampling_points(i)%points,    &
                                   & sample_results(i)%results,    &
                                   & inputs%real_modes,            &
                                   & weighted_energy_force_ratio,  &
                                   & this)
    
    this%basis_functions = [this%basis_functions, basis_functions(i)%functions]
    this%coefficients    = [this%coefficients, coefficients]
  enddo
end subroutine

! Calculate the energy at a given displacement.
impure elemental function energy_PolynomialPotential(this,displacement) &
   & result(output)
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  real(dp)                               :: output
  
  output = this%energy_baseline   &
       & + sum( this%coefficients &
       &      * this%basis_functions%evaluate(displacement))
end function

! Calculate the force at a given displacement.
impure elemental function force_PolynomialPotential(this,displacement) &
   & result(output)
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  type(RealModeForce)                    :: output
  
  output = RealModeForce(sum( this%coefficients &
                          & * this%basis_functions%derivative(displacement)))
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_PolynomialPotential(this,input)
  implicit none
  
  class(PolynomialPotential), intent(out) :: this
  type(String),               intent(in)  :: input(:)
  
  select type(this); type is(PolynomialPotential)
  end select
end subroutine

function write_PolynomialPotential(this) result(output)
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  type(String), allocatable              :: output(:)
  
  select type(this); type is(PolynomialPotential)
  end select
end function

impure elemental function new_PolynomialPotential_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(PolynomialPotential)     :: this
  
  this = input
end function
end module
