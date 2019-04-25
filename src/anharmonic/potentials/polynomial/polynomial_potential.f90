! ======================================================================
! A polynomial representation of a potential.
! ======================================================================
module polynomial_potential_module
  use common_module
  
  use states_module
  use anharmonic_common_module
  
  use basis_function_module
  use coupling_basis_functions_module
  use stress_basis_function_module
  use coupling_stress_basis_functions_module
  use vscf_rvectors_module
  use sampling_points_module
  use sample_result_module
  use sample_results_module
  use fit_coefficients_module
  use fit_stress_coefficients_module
  use polynomial_stress_module
  implicit none
  
  private
  
  public :: startup_polynomial_potential
  
  public :: PolynomialPotential
  
  type, extends(PotentialData) :: PolynomialPotential
    integer,  private :: potential_expansion_order
    real(dp), private :: reference_energy
    type(CouplingBasisFunctions), allocatable, private :: basis_functions_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_PolynomialPotential
    
    procedure, public :: generate_sampling_points => &
       & generate_sampling_points_PolynomialPotential
    procedure, public :: generate_potential => &
       & generate_potential_PolynomialPotential
    procedure, public :: generate_stress => &
       & generate_stress_PolynomialPotential
    
    procedure, public :: zero_energy => zero_energy_PolynomialPotential
    
    procedure, public :: energy_RealModeDisplacement => &
                       & energy_RealModeDisplacement_PolynomialPotential
    procedure, public :: energy_ComplexModeDisplacement => &
                       & energy_ComplexModeDisplacement_PolynomialPotential
    procedure, public :: force_RealModeDisplacement => &
                       & force_RealModeDisplacement_PolynomialPotential
    procedure, public :: force_ComplexModeDisplacement => &
                       & force_ComplexModeDisplacement_PolynomialPotential
    
    procedure, public :: braket_state  => braket_state_PolynomialPotential
    procedure, public :: braket_states => braket_states_PolynomialPotential
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_PolynomialPotential
    
    procedure, public :: iterate_damped => iterate_damped_PolynomialPotential
    procedure, public :: iterate_pulay => iterate_pulay_PolynomialPotential
    
    ! I/O.
    procedure, public :: read  => read_PolynomialPotential
    procedure, public :: write => write_PolynomialPotential
  end type
  
  interface PolynomialPotential
    module procedure new_PolynomialPotential
    module procedure new_PolynomialPotential_PotentialData
    module procedure new_PolynomialPotential_BasisFunctions
    module procedure new_PolynomialPotential_Strings
    module procedure new_PolynomialPotential_StringArray
  end interface
contains

! Startup procedure.
subroutine startup_polynomial_potential()
  implicit none
  
  type(PolynomialPotential) :: potential
  
  call potential%startup()
end subroutine

! Constructors.
function new_PolynomialPotential(potential_expansion_order) result(this)
  implicit none
  
  integer, intent(in)       :: potential_expansion_order
  type(PolynomialPotential) :: this
  
  this%potential_expansion_order = potential_expansion_order
end function

recursive function new_PolynomialPotential_PotentialData(input) result(this)
  implicit none
  
  class(PotentialData), intent(in) :: input
  type(PolynomialPotential)        :: this
  
  select type(input); type is(PolynomialPotential)
    this = input
  type is(PotentialPointer)
    this = PolynomialPotential(input%potential())
  class default
    call err()
  end select
end function

function new_PolynomialPotential_BasisFunctions(potential_expansion_order, &
   & reference_energy,basis_functions) result(this)
  implicit none
  
  integer,                      intent(in) :: potential_expansion_order
  real(dp),                     intent(in) :: reference_energy
  type(CouplingBasisFunctions), intent(in) :: basis_functions(:)
  type(PolynomialPotential)                :: this
  
  this%potential_expansion_order = potential_expansion_order
  this%reference_energy          = reference_energy
  this%basis_functions_          = basis_functions
end function

! Type representation.
impure elemental function representation_PolynomialPotential() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'polynomial'
end function

! Generate sampling points.
subroutine generate_sampling_points_PolynomialPotential(this, &
   & anharmonic_data,sampling_points_dir,calculation_writer,logfile)
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  type(AnharmonicData),       intent(in)    :: anharmonic_data
  type(String),               intent(in)    :: sampling_points_dir
  type(CalculationWriter),    intent(inout) :: calculation_writer
  type(OFile),                intent(inout) :: logfile
  
  ! Variables used when generating sampling points.
  type(BasisFunctionsAndSamplingPoints)     :: basis_functions_and_points
  type(CouplingBasisFunctions), allocatable :: basis_functions(:)
  type(SamplingPoints),         allocatable :: sampling_points(:)
  
  ! Supercell variables.
  type(RealMode),   allocatable :: sampling_point_modes(:)
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
  type(String) :: equilibrium_dir
  type(String) :: coupling_dir
  type(OFile)  :: basis_functions_file
  type(OFile)  :: sampling_points_file
  type(String) :: sampling_dir
  type(OFile)  :: supercell_file
  type(OFile)  :: vscf_rvectors_file
  type(String) :: vscf_rvectors_dir
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Generate basis functions and sampling points.
  ! --------------------------------------------------
  
  ! Loop over subspace couplings, generating basis functions and sampling
  !    points for each.
  allocate( basis_functions(size(anharmonic_data%subspace_couplings)), &
          & sampling_points(size(anharmonic_data%subspace_couplings)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(anharmonic_data%subspace_couplings)
    ! Generate basis functions at each coupling,
    !    and the sampling points from which the basis function coefficients
    !    can be constructed.
    basis_functions_and_points = generate_basis_functions( &
          & anharmonic_data%subspace_couplings(i),         &
          & this%potential_expansion_order,                &
          & anharmonic_data%structure,                     &
          & anharmonic_data%complex_modes,                 &
          & anharmonic_data%real_modes,                    &
          & anharmonic_data%qpoints,                       &
          & anharmonic_data%degenerate_subspaces,          &
          & anharmonic_data%degenerate_symmetries,         &
          & anharmonic_data%vscf_basis_functions_only,     &
          & anharmonic_data%maximum_weighted_displacement, &
          & anharmonic_data%frequency_of_max_displacement, &
          & logfile                                        )
    basis_functions(i) = basis_functions_and_points%basis_functions
    sampling_points(i) = basis_functions_and_points%sampling_points
  enddo
  
  ! --------------------------------------------------
  ! Write out basis functions and sampling points.
  ! --------------------------------------------------
  
  ! Write out sampling point at equilibrium.
  equilibrium_dir = sampling_points_dir//'/equilibrium'
  call calculation_writer%write_calculation( anharmonic_data%structure, &
                                           & equilibrium_dir            )
  
  ! Write out all other sampling points.
  do i=1,size(sampling_points)
    ! Make a directory for each coupling.
    coupling_dir = sampling_points_dir// &
       & '/coupling_'//left_pad(i,str(size(sampling_points)))
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
      sampling_point_modes = select_modes( sampling_point%vectors,    &
                                         & anharmonic_data%real_modes )
      sampling_point_qpoints = select_qpoints( sampling_point_modes,   &
                                             & anharmonic_data%qpoints )
      supercell_matrix = construct_supercell_matrix( &
                         & sampling_point_qpoints,   &
                         & anharmonic_data%structure )
      supercell = construct_supercell( anharmonic_data%structure, &
                                     & supercell_matrix           )
      
      ! Write out the supercell.
      supercell_file = OFile(sampling_dir//'/structure.dat')
      call supercell_file%print_lines(supercell)
      
      ! Construct VSCF R-vectors.
      vscf_rvectors = construct_vscf_rvectors( sampling_point,             &
                                             & supercell,                  &
                                             & anharmonic_data%real_modes, &
                                             & anharmonic_data%qpoints     )
      vscf_rvectors_file = OFile(sampling_dir//'/vscf_rvectors.dat')
      call vscf_rvectors_file%print_lines(vscf_rvectors,separating_line='')
      
      do k=1,size(vscf_rvectors)
        ! Transform the sampling point by the VSCF R-vector.
        transformed_sampling_point = vscf_rvectors(k)%transform( &
                                   & sampling_point,             &
                                   & anharmonic_data%real_modes, &
                                   & anharmonic_data%qpoints     )
        
        ! Construct displaced structure.
        displacement = CartesianDisplacement( transformed_sampling_point, &
                                            & supercell,                  &
                                            & anharmonic_data%real_modes, &
                                            & anharmonic_data%qpoints     )
        displaced_structure = displace_structure(supercell,displacement)
        
        ! Create directory and structure files for displaced structure.
        vscf_rvectors_dir = sampling_dir// &
           & '/vscf_rvectors_'//left_pad(k,str(size(vscf_rvectors)))
        call calculation_writer%write_calculation( displaced_structure, &
                                                 & vscf_rvectors_dir    )
      enddo
    enddo
  enddo
end subroutine

! Generate potential.
subroutine generate_potential_PolynomialPotential(this,anharmonic_data,  &
   & weighted_energy_force_ratio,sampling_points_dir,calculation_reader, &
   & logfile)
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  type(AnharmonicData),       intent(in)    :: anharmonic_data
  real(dp),                   intent(in)    :: weighted_energy_force_ratio
  type(String),               intent(in)    :: sampling_points_dir
  type(CalculationReader),    intent(inout) :: calculation_reader
  type(OFile),                intent(inout) :: logfile
  
  ! Basis functions and sampling points.
  type(CouplingBasisFunctions), allocatable :: basis_functions(:)
  type(SamplingPoints),         allocatable :: sampling_points(:)
  type(SampleResults),          allocatable :: sample_results(:)
  
  ! Basis functions and sampling points corresponding to the un-displaced
  !    structure and the constant basis function.
  type(BasisFunction)         :: constant_basis_function
  type(RealModeDisplacement)  :: equilibrium_sampling_point
  type(SampleResult)          :: equilibrium_sample_result
  
  ! Basis function coefficients.
  real(dp), allocatable :: coefficients(:)
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Read in basis functions and sampling points.
  ! --------------------------------------------------
  basis_functions = read_basis_functions( anharmonic_data%subspace_couplings, &
                                        & sampling_points_dir                 )
  constant_basis_function = generate_constant_basis_function()
  
  sampling_points = read_sampling_points( anharmonic_data%subspace_couplings, &
                                        & sampling_points_dir                 )
  equilibrium_sampling_point = generate_equilibrium_sampling_point()
  
  ! --------------------------------------------------
  ! Read in electronic structure data.
  ! --------------------------------------------------
  sample_results = read_sample_results( anharmonic_data%subspace_couplings, &
                                      & sampling_points,                    &
                                      & sampling_points_dir,                &
                                      & anharmonic_data,                    &
                                      & calculation_reader                  )
  equilibrium_sample_result = read_equilibrium_sample_result( &
                                & equilibrium_sampling_point, &
                                & sampling_points_dir,        &
                                & anharmonic_data,            &
                                & calculation_reader          )
  
  ! --------------------------------------------------
  ! Fit basis function coefficients.
  ! --------------------------------------------------
  
  this%reference_energy = equilibrium_sample_result%energy
  this%basis_functions_ = [CouplingBasisFunctions::]
  
  do i=1,size(anharmonic_data%subspace_couplings)
    coefficients = fit_coefficients( basis_functions(i)%basis_functions(), &
                                   & sampling_points(i)%points,            &
                                   & sample_results(i)%results,            &
                                   & anharmonic_data%real_modes,           &
                                   & weighted_energy_force_ratio,          &
                                   & this                                  )
    
    this%basis_functions_ = [                                        &
       & this%basis_functions_,                                      &
       & CouplingBasisFunctions(                                     &
       &    coupling = basis_functions(i)%coupling,                  &
       &    basis_functions = coefficients                           &
       &                    * basis_functions(i)%basis_functions() ) ]
  enddo
end subroutine

! Generate stress.
function generate_stress_PolynomialPotential(this,anharmonic_data,        &
   & sampling_points_dir,stress_expansion_order,stress_subspace_coupling, &
   & vscf_basis_functions_only,calculation_reader,logfile) result(output)
  implicit none
  
  class(PolynomialPotential), intent(in)    :: this
  type(AnharmonicData),       intent(in)    :: anharmonic_data
  type(String),               intent(in)    :: sampling_points_dir
  integer,                    intent(in)    :: stress_expansion_order
  type(SubspaceCoupling),     intent(in)    :: stress_subspace_coupling(:)
  logical,                    intent(in)    :: vscf_basis_functions_only
  type(CalculationReader),    intent(inout) :: calculation_reader
  type(OFile),                intent(inout) :: logfile
  type(StressPointer)                       :: output
  
  ! Stress basis functions.
  type(CouplingStressBasisFunctions), allocatable :: basis_functions(:)
  
  ! Electronic structure results.
  type(SamplingPoints), allocatable :: sampling_points(:)
  type(SampleResults),  allocatable :: sample_results(:)
  type(SamplingPoints), allocatable :: stress_sampling_points(:)
  type(SampleResults),  allocatable :: stress_sample_results(:)
  
  ! Electronic structure results corresponding to the un-displaced structure.
  type(RealModeDisplacement) :: equilibrium_sampling_point
  type(SampleResult)         :: equilibrium_sample_result
  
  ! The stress.
  type(PolynomialStress) :: stress
  
  ! Basis function coefficients.
  real(dp), allocatable :: coefficients(:)
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  ! Generate basis functions.
  allocate( basis_functions(size(stress_subspace_coupling)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(basis_functions)
    basis_functions(i) = generate_stress_basis_functions( &
                 & stress_subspace_coupling(i),           &
                 & stress_expansion_order,                &
                 & anharmonic_data%structure,             &
                 & anharmonic_data%complex_modes,         &
                 & anharmonic_data%real_modes,            &
                 & anharmonic_data%qpoints,               &
                 & anharmonic_data%degenerate_subspaces,  &
                 & anharmonic_data%degenerate_symmetries, &
                 & vscf_basis_functions_only,             &
                 & logfile                                )
  enddo
  
  ! Read in all sampling points and sample results.
  sampling_points = read_sampling_points( anharmonic_data%subspace_couplings, &
                                        & sampling_points_dir                 )
  equilibrium_sampling_point = generate_equilibrium_sampling_point()
  
  sample_results = read_sample_results( anharmonic_data%subspace_couplings, &
                                      & sampling_points,                    &
                                      & sampling_points_dir,                &
                                      & anharmonic_data,                    &
                                      & calculation_reader                  )
  equilibrium_sample_result = read_equilibrium_sample_result( &
                                & equilibrium_sampling_point, &
                                & sampling_points_dir,        &
                                & anharmonic_data,            &
                                & calculation_reader          )
  
  ! Set the un-displaced stress tensor.
  stress = PolynomialStress( equilibrium_sample_result%stress(), &
                           & [CouplingStressBasisFunctions::]    )
  
  ! Fit basis functions.
  do i=1,size(basis_functions)
    j = first(stress_subspace_coupling==basis_functions(i)%coupling)
    
    coefficients = fit_stress_coefficients(    &
       & basis_functions(i)%basis_functions(), &
       & sampling_points(j)%points,            &
       & sample_results(j)%results,            &
       & stress                                )
    
    basis_functions(i) = CouplingStressBasisFunctions(     &
       & basis_functions(i)%coupling,                      &
       & coefficients*basis_functions(i)%basis_functions() )
  enddo
  
  ! Assemble output.
  stress = PolynomialStress( equilibrium_sample_result%stress(), &
                           & basis_functions                     )
  output = StressPointer(stress)
end function

! Helper functions for generate_potential and generate_stress.

! Reads sampling points.
function read_sampling_points(couplings,sampling_points_directory) &
   & result(output)
  implicit none
  
  type(SubspaceCoupling), intent(in) :: couplings(:)
  type(String),           intent(in) :: sampling_points_directory
  type(SamplingPoints), allocatable  :: output(:)
  
  type(String) :: coupling_directory
  type(IFile)  :: sampling_points_file
  
  integer :: i,ialloc
  
  allocate(output(size(couplings)), stat=ialloc); call err(ialloc)
  do i=1,size(couplings)
    coupling_directory = sampling_points_directory// &
       & '/coupling_'//left_pad(i,str(size(couplings)))
    
    sampling_points_file = IFile(coupling_directory//'/sampling_points.dat')
    output(i) = SamplingPoints(sampling_points_file%lines())
  enddo
end function

! Generates the sampling point for the equilibrium structure.
function generate_equilibrium_sampling_point() result(output)
  implicit none
  
  type(RealModeDisplacement) :: output
  
  output = RealModeDisplacement([RealSingleDisplacement::])
end function

! Reads basis functions.
function read_basis_functions(couplings,sampling_points_directory) &
   & result(output)
  implicit none
  
  type(SubspaceCoupling), intent(in)        :: couplings(:)
  type(String),           intent(in)        :: sampling_points_directory
  type(CouplingBasisFunctions), allocatable :: output(:)
  
  type(String) :: coupling_directory
  type(IFile)  :: basis_functions_file
  
  integer :: i,ialloc
  
  allocate(output(size(couplings)), stat=ialloc); call err(ialloc)
  do i=1,size(couplings)
    coupling_directory = sampling_points_directory// &
       & '/coupling_'//left_pad(i,str(size(couplings)))
    
    basis_functions_file = IFile(coupling_directory//'/basis_functions.dat')
    output(i) = CouplingBasisFunctions(basis_functions_file%lines())
  enddo
end function

! Generates the basis function which is just a constant.
function generate_constant_basis_function() result(output)
  implicit none
  
  type(BasisFunction) :: output
  
  type(RealMonomial)    :: constant_real_monomial
  type(ComplexMonomial) :: constant_complex_monomial
  
  constant_real_monomial = RealMonomial( &
      & coefficient = 1.0_dp,            &
      & modes       = [RealUnivariate::] )
  constant_complex_monomial = ComplexMonomial( &
         & coefficient = (1.0_dp,0.0_dp),      &
         & modes       = [ComplexUnivariate::] )
  output = BasisFunction(                                                     &
     & real_representation = RealPolynomial([constant_real_monomial]),        &
     & complex_representation = ComplexPolynomial([constant_complex_monomial]))
end function

! Reads sample results.
function read_sample_results(couplings,sampling_points,            &
   & sampling_points_directory,anharmonic_data,calculation_reader) &
   & result(output)
  implicit none
  
  type(SubspaceCoupling),  intent(in)    :: couplings(:)
  type(SamplingPoints),    intent(in)    :: sampling_points(:)
  type(String),            intent(in)    :: sampling_points_directory
  type(AnharmonicData),    intent(in)    :: anharmonic_data
  type(CalculationReader), intent(inout) :: calculation_reader
  type(SampleResults), allocatable       :: output(:)
  
  type(StructureData)                    :: supercell
  type(CartesianDisplacement)            :: displacement
  type(VscfRvectors),        allocatable :: vscf_rvectors(:)
  type(ElectronicStructure), allocatable :: calculations(:)
  
  ! Files and directories.
  type(String) :: coupling_directory
  type(String) :: sampling_directory
  type(IFile)  :: supercell_file
  type(IFile)  :: vscf_rvectors_file
  type(String) :: vscf_rvectors_directory
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  allocate(output(size(couplings)), stat=ialloc); call err(ialloc)
  do i=1,size(couplings)
    coupling_directory = sampling_points_directory// &
       & '/coupling_'//left_pad(i,str(size(sampling_points)))
    
    allocate( output(i)%results(size(sampling_points(i))), &
            & stat=ialloc); call err(ialloc)
    do j=1,size(sampling_points(i))
      sampling_directory = coupling_directory// &
         & '/sampling_point_'//left_pad(j,str(size(sampling_points(i))))
      
      ! Read in supercell and VSCF R-vectors.
      supercell_file = IFile(sampling_directory//'/structure.dat')
      supercell = StructureData(supercell_file%lines())
      
      displacement = CartesianDisplacement( sampling_points(i)%points(j), &
                                          & supercell,                    &
                                          & anharmonic_data%real_modes,   &
                                          & anharmonic_data%qpoints       )
      
      vscf_rvectors_file = IFile(sampling_directory//'/vscf_rvectors.dat')
      vscf_rvectors = VscfRvectors(vscf_rvectors_file%sections())
      
      ! Read in electronic structure calculations.
      allocate( calculations(size(vscf_rvectors)), &
              & stat=ialloc); call err(ialloc)
      do k=1,size(vscf_rvectors)
        vscf_rvectors_directory = sampling_directory// &
           & '/vscf_rvectors_'//left_pad(k,str(size(vscf_rvectors)))
        calculations(k) = calculation_reader%read_calculation( &
                                    & vscf_rvectors_directory, &
                                    & displacement             )
      enddo
      
      ! Average electronic structure across VSCF R-vectors, and convert
      !    to correct normalisation and real mode co-ordinates.
      output(i)%results(j) = SampleResult( vscf_rvectors,              &
                                         & calculations,               &
                                         & supercell,                  &
                                         & anharmonic_data%real_modes, &
                                         & anharmonic_data%qpoints     )
      
      deallocate(calculations, stat=ialloc); call err(ialloc)
    enddo
  enddo
end function

! Reads the sample result corresponding to equilibrium.
function read_equilibrium_sample_result(equilibrium_sampling_point, &
   & sampling_points_directory,anharmonic_data,calculation_reader)  &
   & result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in)    :: equilibrium_sampling_point
  type(String),               intent(in)    :: sampling_points_directory
  type(AnharmonicData),       intent(in)    :: anharmonic_data
  type(CalculationReader),    intent(inout) :: calculation_reader
  type(SampleResult)                        :: output
  
  type(CartesianDisplacement) :: displacement
  type(ElectronicStructure)   :: electronic_structure
  
  ! The vector of length 0 as a CartesianDisplacement.
  displacement = CartesianDisplacement( &
          & equilibrium_sampling_point, &
          & anharmonic_data%structure,  &
          & anharmonic_data%real_modes, &
          & anharmonic_data%qpoints     )
  
  electronic_structure = calculation_reader%read_calculation( &
                 & sampling_points_directory//'/equilibrium', &
                 & displacement                               )
  
  output = SampleResult( electronic_structure,       &
                       & anharmonic_data%structure,  &
                       & anharmonic_data%real_modes, &
                       & anharmonic_data%qpoints     )
end function

! Set the undisplaced energy to zero.
impure elemental subroutine zero_energy_PolynomialPotential(this)
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  
  this%reference_energy = this%reference_energy - this%undisplaced_energy()
end subroutine

! Calculate the energy at a given displacement.
impure elemental function energy_RealModeDisplacement_PolynomialPotential( &
   & this,displacement) result(output)
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  real(dp)                               :: output
  
  output = this%reference_energy &
       & + sum(this%basis_functions_%energy(displacement))
end function

impure elemental function energy_ComplexModeDisplacement_PolynomialPotential( &
   & this,displacement) result(output)
  implicit none
  
  class(PolynomialPotential),    intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  complex(dp)                               :: output
  
  output = this%reference_energy &
       & + sum(this%basis_functions_%energy(displacement))
end function

! Calculate the force at a given displacement.
impure elemental function force_RealModeDisplacement_PolynomialPotential( &
   & this,displacement) result(output)
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  type(RealModeForce)                    :: output
  
  if (size(this%basis_functions_)>0) then
    output = sum(this%basis_functions_%force(displacement))
  else
    output = RealModeForce(RealSingleForce(displacement%vectors%id,0.0_dp))
  endif
end function

impure elemental function force_ComplexModeDisplacement_PolynomialPotential( &
   & this,displacement) result(output)
  implicit none
  
  class(PolynomialPotential),    intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexModeForce)                    :: output
  
  if (size(this%basis_functions_)>0) then
    output = sum(this%basis_functions_%force(displacement))
  else
    output = ComplexModeForce(ComplexSingleForce( displacement%vectors%id, &
                                                & (0.0_dp,0.0_dp)          ))
  endif
end function

! Integrate the potential between two states.
subroutine braket_state_PolynomialPotential(this,bra,ket,subspace, &
   & subspace_basis,anharmonic_data,qpoint)
  implicit none
  
  class(PolynomialPotential), intent(inout)        :: this
  class(SubspaceState),       intent(in)           :: bra
  class(SubspaceState),       intent(in), optional :: ket
  type(DegenerateSubspace),   intent(in)           :: subspace
  class(SubspaceBasis),       intent(in)           :: subspace_basis
  type(AnharmonicData),       intent(in)           :: anharmonic_data
  type(QpointData),           intent(in), optional :: qpoint
  
  logical, allocatable :: to_remove(:)
  
  integer :: i,j,k
  
  ! Integrate the reference energy (N.B. <i|e|j> = e<i|j> if e is a scalar.).
  this%reference_energy = this%reference_energy          &
                      & * inner_product( bra,            &
                      &                  ket,            &
                      &                  subspace,       &
                      &                  subspace_basis, &
                      &                  anharmonic_data )
  
  ! Integrate each basis function between the bra and the ket.
  do i=1,size(this%basis_functions_)
    call this%basis_functions_(i)%braket( bra,             &
                                        & ket,             &
                                        & subspace,        &
                                        & subspace_basis,  &
                                        & anharmonic_data, &
                                        & qpoint           )
  enddo
  
  ! Check if the basis function is now a constant.
  ! If so, add the constant energy to the potential's reference energy,
  !    and flag the term for removal.
  ! Then check if a coupling is now the same as a previous coupling.
  ! If so, combine the two and flag the duplicate term for removal.
  to_remove = [(.false., i=1, size(this%basis_functions_))]
  do i=1,size(this%basis_functions_)
    if (size(this%basis_functions_(i)%coupling)==0) then
      this%reference_energy =      &
         &   this%reference_energy &
         & + this%basis_functions_(i)%undisplaced_energy()
      to_remove(i) = .true.
    else
      do j=1,i-1
        if (    this%basis_functions_(i)%coupling &
           & == this%basis_functions_(j)%coupling ) then
          if (to_remove(j)) then
            call err()
          endif
          call this%basis_functions_(j)%append( &
                     & this%basis_functions_(i) )
          to_remove(i) = .true.
        endif
      enddo
    endif
  enddo
  
  ! Remove constant and duplicate terms.
  this%basis_functions_ = this%basis_functions_(filter(.not.to_remove) )
end subroutine

subroutine braket_states_PolynomialPotential(this,states,subspace, &
   & subspace_basis,anharmonic_data)
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  class(SubspaceStates),      intent(in)    :: states
  type(DegenerateSubspace),   intent(in)    :: subspace
  class(SubspaceBasis),       intent(in)    :: subspace_basis
  type(AnharmonicData),       intent(in)    :: anharmonic_data
  
  logical, allocatable :: to_remove(:)
  
  integer :: i,j,k
  
  ! Integrate each basis function between the bra and the ket.
  do i=1,size(this%basis_functions_)
    call this%basis_functions_(i)%braket( states,         &
                                        & subspace,       &
                                        & subspace_basis, &
                                        & anharmonic_data )
  enddo
  
  ! Check if the basis function is now a constant.
  ! If so, add the constant energy to the potential's reference energy,
  !    and flag the term for removal.
  ! Then check if a coupling is now the same as a previous coupling.
  ! If so, combine the two and flag the duplicate term for removal.
  to_remove = [(.false., i=1, size(this%basis_functions_))]
  do i=1,size(this%basis_functions_)
    if (size(this%basis_functions_(i)%coupling)==0) then
      this%reference_energy =      &
         &   this%reference_energy &
         & + this%basis_functions_(i)%undisplaced_energy()
      to_remove(i) = .true.
    else
      do j=1,i-1
        if (    this%basis_functions_(i)%coupling &
           & == this%basis_functions_(j)%coupling ) then
          if (to_remove(j)) then
            call err()
          endif
          call this%basis_functions_(j)%append( &
                     & this%basis_functions_(i) )
          to_remove(i) = .true.
        endif
      enddo
    endif
  enddo
  
  ! Remove constant and duplicate terms.
  this%basis_functions_ = this%basis_functions_(filter(.not.to_remove) )
end subroutine

! Calculate the thermal expectation of the potential, <V>, for a set of
!    harmonic states.
function harmonic_expectation_PolynomialPotential(this,frequency, &
   & thermal_energy,no_states,subspace,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  real(dp),                   intent(in) :: frequency
  real(dp),                   intent(in) :: thermal_energy
  integer,                    intent(in) :: no_states
  type(DegenerateSubspace),   intent(in) :: subspace
  type(AnharmonicData),       intent(in) :: anharmonic_data
  real(dp)                               :: output
  
  output = this%reference_energy                                           &
       & + sum(this%basis_functions_%harmonic_expectation( frequency,      &
       &                                                   thermal_energy, &
       &                                                   no_states,      &
       &                                                   subspace,       &
       &                                                   anharmonic_data ))
end function

! Generates the next iteration of the potential, either following a damped
!    iterative scheme or a pulay scheme.
impure elemental function iterate_damped_PolynomialPotential(this, &
   & new_potential,damping,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  class(PotentialData),       intent(in) :: new_potential
  real(dp),                   intent(in) :: damping
  type(AnharmonicData),       intent(in) :: anharmonic_data
  type(PotentialPointer)                 :: output
  
  type(PolynomialPotential) :: poly_potential
  
  type(PolynomialPotential) :: potential
  
  integer :: i
  
  poly_potential = PolynomialPotential(new_potential)
  
  if (size(this%basis_functions_)/=size(poly_potential%basis_functions_)) then
    call err()
  endif
  potential = this
  potential%reference_energy = damping*this%reference_energy &
                           & + (1-damping)*poly_potential%reference_energy
  do i=1,size(potential%basis_functions_)
    if (    size(this%basis_functions_(i))           &
       & /= size(poly_potential%basis_functions_(i)) ) then
      call err()
    endif
    call potential%basis_functions_(i)%set_coefficients(                   &
       &   damping * this%basis_functions_(i)%coefficients()               &
       & + (1-damping) * poly_potential%basis_functions_(i)%coefficients() )
  enddo
  
  output = PotentialPointer(potential)
end function

function iterate_pulay_PolynomialPotential(this,input_potentials, &
   & output_potentials,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  type(PotentialPointer),     intent(in) :: input_potentials(:)
  type(PotentialPointer),     intent(in) :: output_potentials(:)
  type(AnharmonicData),       intent(in) :: anharmonic_data
  type(PotentialPointer)                 :: output
  
  type(RealVector), allocatable :: input_coefficients(:)
  type(RealVector), allocatable :: output_coefficients(:)
  real(dp),         allocatable :: pulay_coefficients(:)
  
  type(PolynomialPotential), allocatable :: input_potential
  type(PolynomialPotential), allocatable :: output_potential
  
  type(PolynomialPotential) :: potential
  
  integer :: i,j,ialloc
  
  if (size(input_potentials)/=size(output_potentials)) then
    call err()
  endif
  
  potential = this
  
  ! Convert previous potential iterations to vectors of coefficients.
  allocate( input_coefficients(size(input_potentials)),  &
          & output_coefficients(size(input_potentials)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(input_potentials)
    input_potential = PolynomialPotential(input_potentials(i))
    input_coefficients(i) = vec([(                           &
       & input_potential%basis_functions_(j)%coefficients(), &
       & j=1,                                                &
       & size(input_potential%basis_functions_)              )])
    
    output_potential = PolynomialPotential(output_potentials(i))
    output_coefficients(i) = vec([(                           &
       & output_potential%basis_functions_(j)%coefficients(), &
       & j=1,                                                 &
       & size(output_potential%basis_functions_)              )])
  enddo
  
  ! Use a Pulay scheme to construct the next set of coefficients.
  pulay_coefficients = dble(pulay( input_coefficients, &
                                 & output_coefficients ))
  
  ! Set the output coefficients from the Pulay coefficients.
  j = 0
  do i=1,size(potential%basis_functions_)
    call potential%basis_functions_(i)%set_coefficients(                &
        & pulay_coefficients(j+1:j+size(potential%basis_functions_(i))) )
    j = j + size(potential%basis_functions_(i))
  enddo
  if (j/=size(pulay_coefficients)) then
    call err()
  endif
  
  output = PotentialPointer(potential)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_PolynomialPotential(this,input)
  implicit none
  
  class(PolynomialPotential), intent(out) :: this
  type(String),               intent(in)  :: input(:)
  
  integer                                   :: potential_expansion_order
  real(dp)                                  :: reference_energy
  type(CouplingBasisFunctions), allocatable :: basis_functions(:)
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(PolynomialPotential)
    line = split_line(input(1))
    potential_expansion_order = int(line(3))
    
    line = split_line(input(2))
    reference_energy = dble(line(3))
    
    basis_functions = CouplingBasisFunctions(split_into_sections( &
                                 & input(5:),                     &
                                 & separating_line=repeat('=',50) ))
    
    this = PolynomialPotential( potential_expansion_order, &
                              & reference_energy,          &
                              & basis_functions            )
  class default
    call err()
  end select
end subroutine

function write_PolynomialPotential(this) result(output)
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  type(String), allocatable              :: output(:)
  
  select type(this); type is(PolynomialPotential)
    output = [ 'Expansion order: '//this%potential_expansion_order,       &
             & 'Reference energy: '//this%reference_energy,               &
             & str('Basis functions:'),                                   &
             & str(''),                                                   &
             & str(this%basis_functions_, separating_line=repeat('=',50)) ]
  class default
    call err()
  end select
end function

function new_PolynomialPotential_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)  :: input(:)
  type(PolynomialPotential) :: this
  
  call this%read(input)
end function

impure elemental function new_PolynomialPotential_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(PolynomialPotential)     :: this
  
  this = PolynomialPotential(str(input))
end function
end module
