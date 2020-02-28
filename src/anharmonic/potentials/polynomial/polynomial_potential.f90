! ======================================================================
! A polynomial representation of a potential.
! ======================================================================
module polynomial_potential_module
  use common_module
  
  use states_module
  use anharmonic_common_module
  
  use polynomial_interpolator_module
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
    integer,  private :: potential_expansion_order_
    real(dp), private :: reference_energy_
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
    procedure, public :: add_constant => add_constant_PolynomialPotential
    
    procedure, public :: finalise_subspace_potential => &
                       & finalise_subspace_potential_PolynomialPotential
    
    procedure, public :: energy_RealModeDisplacement => &
                       & energy_RealModeDisplacement_PolynomialPotential
    procedure, public :: energy_ComplexModeDisplacement => &
                       & energy_ComplexModeDisplacement_PolynomialPotential
    procedure, public :: force_RealModeDisplacement => &
                       & force_RealModeDisplacement_PolynomialPotential
    procedure, public :: force_ComplexModeDisplacement => &
                       & force_ComplexModeDisplacement_PolynomialPotential
    
    procedure, public :: braket_SubspaceBraket  => &
                       & braket_SubspaceBraket_PolynomialPotential
    procedure, public :: braket_BasisState  => &
                       & braket_BasisState_PolynomialPotential
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_PolynomialPotential
    
    procedure, public :: simplify => &
                       & simplify_PolynomialPotential
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_PolynomialPotential
    
    procedure, public :: potential_energy_SubspaceBraKet => &
                       & potential_energy_SubspaceBraKet_PolynomialPotential
    procedure, public :: potential_energy_BasisState => &
                       & potential_energy_BasisState_PolynomialPotential
    
    procedure, public :: coefficients => &
                       & coefficients_PolynomialPotential
    procedure, public :: set_coefficients => &
                       & set_coefficients_PolynomialPotential
    
    procedure, public :: can_be_interpolated => &
                       & can_be_interpolated_PolynomialPotential
    procedure, public :: interpolate => &
                       & interpolate_PolynomialPotential
    
    procedure, public :: interpolate_coefficient => &
                       & interpolate_coefficient_PolynomialPotential
    
    procedure, public :: calculate_dynamical_matrices => &
                       & calculate_dynamical_matrices_PolynomialPotential
    
    procedure, public :: energy_correction => &
                       & energy_correction_PolynomialPotential
    
    procedure, public :: expansion_order => expansion_order_PolynomialPotential
    
    ! I/O.
    procedure, public :: read  => read_PolynomialPotential
    procedure, public :: write => write_PolynomialPotential
  end type
  
  interface PolynomialPotential
    module procedure new_PolynomialPotential
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
  
  this%potential_expansion_order_ = potential_expansion_order
end function

function new_PolynomialPotential_BasisFunctions(potential_expansion_order, &
   & reference_energy,basis_functions) result(this)
  implicit none
  
  integer,                      intent(in) :: potential_expansion_order
  real(dp),                     intent(in) :: reference_energy
  type(CouplingBasisFunctions), intent(in) :: basis_functions(:)
  type(PolynomialPotential)                :: this
  
  this%potential_expansion_order_ = potential_expansion_order
  this%reference_energy_          = reference_energy
  this%basis_functions_           = basis_functions
end function

! Type representation.
impure elemental function representation_PolynomialPotential() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'polynomial'
end function

! Generate sampling points.
! N.B. does not look at use_forces, use_hessians or calculate stress.
subroutine generate_sampling_points_PolynomialPotential(this,       &
   & anharmonic_data,use_forces,energy_to_force_ratio,use_hessians, &
   & calculate_stress,sampling_points_dir,calculation_writer,logfile)
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  type(AnharmonicData),       intent(in)    :: anharmonic_data
  logical,                    intent(in)    :: use_forces
  real(dp),                   intent(in)    :: energy_to_force_ratio
  logical,                    intent(in)    :: use_hessians
  logical,                    intent(in)    :: calculate_stress
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
  
  if (.not. use_forces) then
    call print_line(ERROR//': The polynomial potential cannot currently not &
       &use forces.')
    call err()
  elseif (use_hessians) then
    call print_line(ERROR//': The polynomial potential cannot currently &
       &use hessians.')
    call err()
  endif
  
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
    if (modulo(i,max(size(anharmonic_data%subspace_couplings)/100,1))==0) then
      call print_line('Generating sampling points in subspace coupling '// &
         & i//' of '//size(anharmonic_data%subspace_couplings)//'.')
    endif
    basis_functions_and_points = generate_basis_functions( &
          & anharmonic_data%subspace_couplings(i),         &
          & this%potential_expansion_order_,               &
          & anharmonic_data%structure,                     &
          & anharmonic_data%complex_modes,                 &
          & anharmonic_data%real_modes,                    &
          & anharmonic_data%qpoints,                       &
          & anharmonic_data%degenerate_subspaces,          &
          & anharmonic_data%degenerate_symmetries,         &
          & anharmonic_data%vscf_basis_functions_only,     &
          & anharmonic_data%maximum_weighted_displacement, &
          & anharmonic_data%frequency_of_max_displacement, &
          & energy_to_force_ratio,                         &
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
    if (modulo(i,max(size(anharmonic_data%subspace_couplings)/100,1))==0) then
      call print_line('Generating calculations in subspace coupling '// &
         & i//' of '//size(anharmonic_data%subspace_couplings)//'.')
    endif
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
  
  ! Basis functions and sampling points corresponding to the un-displaced
  !    structure and the constant basis function.
  type(BasisFunction)         :: constant_basis_function
  type(RealModeDisplacement)  :: equilibrium_sampling_point
  type(SampleResult)          :: equilibrium_sample_result
  
  ! Basis functions and sampling points at each coupling.
  type(SubspaceCoupling)       :: coupling
  type(String)                 :: coupling_directory
  type(CouplingBasisFunctions) :: basis_functions
  type(SamplingPoints)         :: sampling_points
  type(SampleResults)          :: sample_results
  
  ! Basis function coefficients.
  real(dp), allocatable :: coefficients(:)
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  ! Generate the potential at zero displacement.
  call print_line('Generating potential at zero displacement.')
  constant_basis_function = generate_constant_basis_function()
  equilibrium_sampling_point = generate_equilibrium_sampling_point()
  equilibrium_sample_result = read_equilibrium_sample_result( &
                                & equilibrium_sampling_point, &
                                & sampling_points_dir,        &
                                & anharmonic_data,            &
                                & calculation_reader          )
  this%reference_energy_ = equilibrium_sample_result%energy
  
  allocate(this%basis_functions_(0), stat=ialloc); call err(ialloc)
  
  do i=1,size(anharmonic_data%subspace_couplings)
    call print_line('Fitting potential in subspace coupling '//i//' of ' // &
       &size(anharmonic_data%subspace_couplings)//', containing &
       &subspaces '//anharmonic_data%subspace_couplings(i)%ids//'.')
    
    coupling = anharmonic_data%subspace_couplings(i)
    
    ! Generate directory names.
    coupling_directory = sampling_points_dir//'/coupling_'//left_pad( &
                      & i,                                            &
                      & str(size(anharmonic_data%subspace_couplings)) )
    
    ! Read in basis functions.
    basis_functions = read_basis_functions(coupling_directory)
    
    ! Read in sampling points.
    sampling_points = read_sampling_points(coupling_directory)
    
    ! Read in electronic structure data.
    sample_results = read_sample_results( sampling_points,    &
                                        & coupling_directory, &
                                        & anharmonic_data,    &
                                        & calculation_reader  )
    
    ! Fit basis function coefficients.
    coefficients = fit_coefficients( basis_functions%basis_functions(), &
                                   & sampling_points%points,            &
                                   & sample_results%results,            &
                                   & anharmonic_data%real_modes,        &
                                   & weighted_energy_force_ratio,       &
                                   & this                               )
    
    ! Add fitted basis functions to potential.
    this%basis_functions_ = [                                     &
       & this%basis_functions_,                                   &
       & CouplingBasisFunctions(                                  &
       &    coupling = basis_functions%coupling,                  &
       &    basis_functions = coefficients                        &
       &                    * basis_functions%basis_functions() ) ]
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
  
  ! Electronic structure results corresponding to the un-displaced structure.
  type(CouplingStressBasisFunctions) :: zero_basis_functions(0)
  type(RealModeDisplacement)         :: equilibrium_sampling_point
  type(SampleResult)                 :: equilibrium_sample_result
  
  ! Stress basis functions.
  type(CouplingStressBasisFunctions), allocatable :: basis_functions(:)
  
  ! Electronic structure results.
  type(String)         :: coupling_directory
  type(SamplingPoints) :: sampling_points
  type(SampleResults)  :: sample_results
  
  ! The stress.
  type(PolynomialStress) :: stress
  
  ! Basis function coefficients.
  real(dp), allocatable :: coefficients(:)
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  ! Generate the stress tensor at zero displacement.
  call print_line('Generating stress at zero displacement.')
  equilibrium_sampling_point = generate_equilibrium_sampling_point()
  equilibrium_sample_result = read_equilibrium_sample_result( &
                                & equilibrium_sampling_point, &
                                & sampling_points_dir,        &
                                & anharmonic_data,            &
                                & calculation_reader          )
  stress = PolynomialStress( stress_expansion_order,             &
                           & equilibrium_sample_result%stress(), &
                           & zero_basis_functions                )
  
  ! Generate basis functions.
  allocate( basis_functions(size(stress_subspace_coupling)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(basis_functions)
    call print_line('Fitting stress in subspace coupling '//i//' of ' // &
       &size(stress_subspace_coupling)//', containing &
       &subspaces '//stress_subspace_coupling(i)%ids//'.')
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
    
    j = first(anharmonic_data%subspace_couplings==stress_subspace_coupling(i))
    coupling_directory = sampling_points_dir//'/coupling_'//left_pad( &
                      & j,                                            &
                      & str(size(anharmonic_data%subspace_couplings)) )
    
    ! Read in all sampling points and sample results.
    sampling_points = read_sampling_points(coupling_directory)
    
    sample_results = read_sample_results( sampling_points,    &
                                        & coupling_directory, &
                                        & anharmonic_data,    &
                                        & calculation_reader  )
    
    ! Fit basis functions.
    j = first(stress_subspace_coupling==basis_functions(i)%coupling)
    
    coefficients = fit_stress_coefficients(    &
       & basis_functions(i)%basis_functions(), &
       & sampling_points%points,               &
       & sample_results%results,               &
       & stress                                )
    
    basis_functions(i) = CouplingStressBasisFunctions(     &
       & basis_functions(i)%coupling,                      &
       & coefficients*basis_functions(i)%basis_functions() )
  enddo
  
  ! Assemble output.
  stress = PolynomialStress( stress_expansion_order,             &
                           & equilibrium_sample_result%stress(), &
                           & basis_functions                     )
  output = StressPointer(stress)
end function

! Helper functions for generate_potential and generate_stress.

! Reads sampling points.
impure elemental function read_sampling_points(coupling_directory) &
   & result(output)
  implicit none
  
  type(String), intent(in) :: coupling_directory
  type(SamplingPoints)     :: output
  
  type(IFile) :: sampling_points_file
  
  sampling_points_file = IFile(coupling_directory//'/sampling_points.dat')
  output = SamplingPoints(sampling_points_file%lines())
end function

! Generates the sampling point for the equilibrium structure.
function generate_equilibrium_sampling_point() result(output)
  implicit none
  
  type(RealModeDisplacement) :: output
  
  type(RealSingleDisplacement) :: zero_displacement(0)
  
  output = RealModeDisplacement(zero_displacement)
end function

! Reads basis functions.
impure elemental function read_basis_functions(coupling_directory) &
   & result(output)
  implicit none
  
  type(String), intent(in)     :: coupling_directory
  type(CouplingBasisFunctions) :: output
  
  type(IFile) :: basis_functions_file
  
  basis_functions_file = IFile(coupling_directory//'/basis_functions.dat')
  output = CouplingBasisFunctions(basis_functions_file%lines())
end function

! Generates the basis function which is just a constant.
function generate_constant_basis_function() result(output)
  implicit none
  
  type(BasisFunction) :: output
  
  type(ComplexMonomial) :: constant_complex_monomial
  
  type(ComplexUnivariate) :: zero_complex(0)
  
  constant_complex_monomial = ComplexMonomial( &
         & coefficient = (1.0_dp,0.0_dp),      &
         & modes       = zero_complex          )
  output = BasisFunction(                                                     &
     & complex_representation = ComplexPolynomial([constant_complex_monomial]))
end function

! Reads sample results.
impure elemental function read_sample_results(sampling_points, &
   & coupling_directory,anharmonic_data,calculation_reader) result(output)
  implicit none
  
  type(SamplingPoints),    intent(in)    :: sampling_points
  type(String),            intent(in)    :: coupling_directory
  type(AnharmonicData),    intent(in)    :: anharmonic_data
  type(CalculationReader), intent(inout) :: calculation_reader
  type(SampleResults)                    :: output
  
  type(StructureData)                    :: supercell
  type(CartesianDisplacement)            :: displacement
  type(VscfRvectors),        allocatable :: vscf_rvectors(:)
  type(ElectronicStructure), allocatable :: calculations(:)
  
  ! Files and directories.
  type(String) :: sampling_directory
  type(IFile)  :: supercell_file
  type(IFile)  :: vscf_rvectors_file
  type(String) :: vscf_rvectors_directory
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  allocate( output%results(size(sampling_points)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(sampling_points)
    sampling_directory = coupling_directory// &
       & '/sampling_point_'//left_pad(i,str(size(sampling_points)))
    
    ! Read in supercell and VSCF R-vectors.
    supercell_file = IFile(sampling_directory//'/structure.dat')
    supercell = StructureData(supercell_file%lines())
    
    displacement = CartesianDisplacement( sampling_points%points(i),  &
                                        & supercell,                  &
                                        & anharmonic_data%real_modes, &
                                        & anharmonic_data%qpoints     )
    
    vscf_rvectors_file = IFile(sampling_directory//'/vscf_rvectors.dat')
    vscf_rvectors = VscfRvectors(vscf_rvectors_file%sections())
    
    ! Read in electronic structure calculations.
    allocate( calculations(size(vscf_rvectors)), &
            & stat=ialloc); call err(ialloc)
    do j=1,size(vscf_rvectors)
      vscf_rvectors_directory = sampling_directory// &
         & '/vscf_rvectors_'//left_pad(j,str(size(vscf_rvectors)))
      calculations(j) = calculation_reader%read_calculation( &
                                  & vscf_rvectors_directory, &
                                  & displacement             )
    enddo
    
    ! Average electronic structure across VSCF R-vectors, and convert
    !    to correct normalisation and real mode co-ordinates.
    output%results(i) = SampleResult( vscf_rvectors,              &
                                    & calculations,               &
                                    & supercell,                  &
                                    & anharmonic_data%real_modes, &
                                    & anharmonic_data%qpoints,    &
                                    & anharmonic_data             )
    
    deallocate(calculations, stat=ialloc); call err(ialloc)
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
                       & anharmonic_data%qpoints,    &
                       & anharmonic_data             )
end function

! Set the undisplaced energy to zero.
impure elemental subroutine zero_energy_PolynomialPotential(this)
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  
  this%reference_energy_ = this%reference_energy_ - this%undisplaced_energy()
end subroutine

! Add a constant to the energy.
impure elemental subroutine add_constant_PolynomialPotential(this,input)
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  real(dp),                   intent(in)    :: input
  
  this%reference_energy_ = this%reference_energy_ + input
end subroutine

! Finalise a subspace potential.
! Re-arranges basis functions to remove duplicates and separate all monomials.
impure elemental subroutine finalise_subspace_potential_PolynomialPotential( &
   & this,subspace,subspace_basis,anharmonic_data)
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  type(DegenerateSubspace),   intent(in)    :: subspace
  class(SubspaceBasis),       intent(in)    :: subspace_basis
  type(AnharmonicData),       intent(in)    :: anharmonic_data
  
  if (size(this%basis_functions_)/=1) then
    call print_line(CODE_ERROR//': Calling finalise_subspace_potential &
       &on a potential with more than one coupling.')
    call err()
  endif
  
  this%reference_energy_ = this%reference_energy_ &
                       & + this%basis_functions_(1)%undisplaced_energy()
  
  call this%basis_functions_(1)%finalise( subspace,       &
                                        & subspace_basis, &
                                        & anharmonic_data )
end subroutine

! Calculate the energy at a given displacement.
impure elemental function energy_RealModeDisplacement_PolynomialPotential( &
   & this,displacement) result(output)
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  real(dp)                               :: output
  
  output = this%reference_energy_ &
       & + sum(this%basis_functions_%energy(displacement))
end function

impure elemental function energy_ComplexModeDisplacement_PolynomialPotential( &
   & this,displacement) result(output)
  implicit none
  
  class(PolynomialPotential),    intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  complex(dp)                               :: output
  
  output = this%reference_energy_ &
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
subroutine braket_SubspaceBraKet_PolynomialPotential(this,braket, &
   & whole_subspace,anharmonic_data)
  implicit none
  
  class(PolynomialPotential), intent(inout)        :: this
  class(SubspaceBraKet),      intent(in)           :: braket
  logical,                    intent(in), optional :: whole_subspace
  type(AnharmonicData),       intent(in)           :: anharmonic_data
  
  integer :: i
  
  ! Integrate the reference energy (N.B. <i|e|j> = e<i|j> if e is a scalar.).
  this%reference_energy_ = this%reference_energy_ &
                       & * braket%inner_product(anharmonic_data)
  
  ! Integrate each basis function between the bra and the ket.
  do i=1,size(this%basis_functions_)
    call this%basis_functions_(i)%braket( braket,         &
                                        & whole_subspace, &
                                        & anharmonic_data )
  enddo
  
  ! Simplify the potential.
  if (set_default(whole_subspace,.true.)) then
    call this%simplify()
  endif
end subroutine

! Integrate the potential between two states.
subroutine braket_BasisState_PolynomialPotential(this,bra,ket,subspace, &
   & subspace_basis,whole_subspace,anharmonic_data)
  implicit none
  
  class(PolynomialPotential), intent(inout)        :: this
  class(BasisState),          intent(in)           :: bra
  class(BasisState),          intent(in), optional :: ket
  type(DegenerateSubspace),   intent(in)           :: subspace
  class(SubspaceBasis),       intent(in)           :: subspace_basis
  logical,                    intent(in), optional :: whole_subspace
  type(AnharmonicData),       intent(in)           :: anharmonic_data
  
  integer :: i
  
  ! Integrate the reference energy (N.B. <i|e|j> = e<i|j> if e is a scalar.).
  this%reference_energy_ = this%reference_energy_                         &
                       & * subspace_basis%inner_product( bra,            &
                       &                                 ket,            &
                       &                                 subspace,       &
                       &                                 anharmonic_data )
  
  ! Integrate each basis function between the bra and the ket.
  do i=1,size(this%basis_functions_)
    call this%basis_functions_(i)%braket( bra,            &
                                        & ket,            &
                                        & subspace,       &
                                        & subspace_basis, &
                                        & whole_subspace, &
                                        & anharmonic_data )
  enddo
  
  ! Simplify the potential.
  if (set_default(whole_subspace,.true.)) then
    call this%simplify()
  endif
end subroutine

subroutine braket_BasisStates_PolynomialPotential(this,states,thermal_energy, &
   & subspace,subspace_basis,whole_subspace,anharmonic_data)
  implicit none
  
  class(PolynomialPotential), intent(inout)        :: this
  class(BasisStates),         intent(inout)        :: states
  real(dp),                   intent(in)           :: thermal_energy
  type(DegenerateSubspace),   intent(in)           :: subspace
  class(SubspaceBasis),       intent(in)           :: subspace_basis
  logical,                    intent(in), optional :: whole_subspace
  type(AnharmonicData),       intent(in)           :: anharmonic_data
  
  integer :: i
  
  ! Integrate each basis function between the bra and the ket.
  do i=1,size(this%basis_functions_)
    call this%basis_functions_(i)%braket( states,         &
                                        & thermal_energy, &
                                        & subspace,       &
                                        & subspace_basis, &
                                        & whole_subspace, &
                                        & anharmonic_data )
  enddo
  
  ! Simplify the potential.
  if (set_default(whole_subspace,.true.)) then
    call this%simplify()
  endif
end subroutine

! Identify basis functions which are constant,
!    add the constant energy to the potential's reference energy,
!    and remove the term.
! Then identify basis functions with the same coupling as a previous coupling,
!    combine the two and remove the duplicate term.
subroutine simplify_PolynomialPotential(this)
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  
  logical, allocatable :: to_remove(:)
  
  integer :: i,j
  
  to_remove = [(.false., i=1, size(this%basis_functions_))]
  
  do i=1,size(this%basis_functions_)
    if (size(this%basis_functions_(i)%coupling)==0) then
      this%reference_energy_ =      &
         &   this%reference_energy_ &
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
          exit
        endif
      enddo
    endif
  enddo
  
  ! Remove constant and duplicate terms.
  this%basis_functions_ = this%basis_functions_(filter(.not.to_remove))
end subroutine

! Calculate the thermal expectation of the potential, <V>, for a set of
!    harmonic states.
impure elemental function harmonic_expectation_PolynomialPotential(this, &
   & frequency,thermal_energy,supercell_size,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  real(dp),                   intent(in) :: frequency
  real(dp),                   intent(in) :: thermal_energy
  integer,                    intent(in) :: supercell_size
  type(AnharmonicData),       intent(in) :: anharmonic_data
  real(dp)                               :: output
  
  output = this%reference_energy_                                          &
       & + sum(this%basis_functions_%harmonic_expectation( frequency,      &
       &                                                   thermal_energy, &
       &                                                   supercell_size  ))
end function

recursive function potential_energy_SubspaceBraKet_PolynomialPotential(this, &
   & braket,anharmonic_data) result(output) 
  implicit none
  
  class(PolynomialPotential), intent(in)           :: this
  class(SubspaceBraKet),      intent(in)           :: braket
  type(AnharmonicData),       intent(in)           :: anharmonic_data
  real(dp)                                         :: output
  
  output = this%reference_energy_                                      &
       & * braket%inner_product(anharmonic_data)                       &
       & + sum(this%basis_functions_%potential_energy( braket,         &
       &                                               anharmonic_data ))
end function

recursive function potential_energy_BasisState_PolynomialPotential(this,bra, &
   & ket,subspace,subspace_basis,anharmonic_data) result(output) 
  implicit none
  
  class(PolynomialPotential), intent(in)           :: this
  class(BasisState),          intent(in)           :: bra
  class(BasisState),          intent(in), optional :: ket
  type(DegenerateSubspace),   intent(in)           :: subspace
  class(SubspaceBasis),       intent(in)           :: subspace_basis
  type(AnharmonicData),       intent(in)           :: anharmonic_data
  real(dp)                                         :: output
  
  output = this%reference_energy_                                      &
       & * subspace_basis%inner_product( bra,                          &
       &                                 ket,                          &
       &                                 subspace,                     &
       &                                 anharmonic_data )             &
       & + sum(this%basis_functions_%potential_energy( bra,            &
       &                                               ket,            &
       &                                               subspace,       &
       &                                               subspace_basis, &
       &                                               anharmonic_data ))
end function

function coefficients_PolynomialPotential(this,anharmonic_data) &
   & result(output)
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  type(AnharmonicData),       intent(in) :: anharmonic_data
  real(dp), allocatable                  :: output(:)
  
  integer :: i
  
  output = [( this%basis_functions_(i)%coefficients(), &
            & i=1,                                     &
            & size(this%basis_functions_)              )]
end function

subroutine set_coefficients_PolynomialPotential(this,coefficients, &
   & anharmonic_data)
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  real(dp),                   intent(in)    :: coefficients(:)
  type(AnharmonicData),       intent(in)    :: anharmonic_data
  
  integer :: i,j
  
  if (    size(coefficients)                         &
     & /= sum([( size(this%basis_functions_(i)),     &
     &           i=1,                                &
     &           size(this%basis_functions_)     )]) ) then
    call print_line(CODE_ERROR//': Coefficients do not match basis &
       &functions.')
    call err()
  endif
  
  j = 0
  do i=1,size(this%basis_functions_)
    call this%basis_functions_(i)%set_coefficients(         &
       & coefficients(j+1:j+size(this%basis_functions_(i))) )
    j = j+size(this%basis_functions_(i))
  enddo
end subroutine

! Interpolate the potential.
function can_be_interpolated_PolynomialPotential(this) result(output)
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  logical                                :: output
  
  output = .true.
end function

! TODO: integrate and unintegrate polynomial, to give interpolation under
!    VSCF rather than raw interpolation.
function interpolate_PolynomialPotential(this,qpoint,subspace,subspace_modes, &
   & anharmonic_min_images,thermal_energy,subspaces,subspace_bases,           &
   & subspace_states,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialPotential), intent(in)    :: this
  type(RealVector),           intent(in)    :: qpoint
  type(DegenerateSubspace),   intent(in)    :: subspace
  type(ComplexMode),          intent(in)    :: subspace_modes(:)
  type(MinImages),            intent(in)    :: anharmonic_min_images(:,:)
  real(dp),                   intent(in)    :: thermal_energy
  type(DegenerateSubspace),   intent(in)    :: subspaces(:)
  class(SubspaceBasis),       intent(in)    :: subspace_bases(:)
  class(BasisStates),         intent(inout) :: subspace_states(:)
  type(AnharmonicData),       intent(in)    :: anharmonic_data
  type(PotentialPointer)                    :: output
  
  integer :: expansion_order
  
  type(ComplexMode), allocatable :: modes(:)
  
  type(RealVector), allocatable :: qpoints(:)
  
  type(PolynomialInterpolator) :: interpolator
  
  integer :: power
  
  type(ComplexMonomial), allocatable :: monomials(:)
  
  type(BasisFunction), allocatable :: basis_functions(:)
  type(BasisFunction), allocatable :: new_basis_functions(:)
  type(BasisFunction), allocatable :: old_basis_functions(:)
  
  complex(dp) :: coefficient
  
  type(CouplingBasisFunctions) :: coupling_basis_functions
  
  integer :: i,j,k,l,ialloc
  
  expansion_order = this%expansion_order()
  
  modes = subspace_modes(filter(subspace_modes%id<=subspace_modes%paired_id))
  
  allocate(qpoints(size(subspace_modes)), stat=ialloc); call err(ialloc)
  do i=1,size(qpoints)
    if (subspace_modes(i)%id<=subspace_modes(i)%paired_id) then
      qpoints(i) = qpoint
    else
      qpoints(i) = -qpoint
    endif
  enddo
  
  ! Construct the polynomial interpolator.
  interpolator = PolynomialInterpolator(                &
     & fine_modes      = subspace_modes,                &
     & fine_qpoints    = qpoints,                       &
     & coarse_modes    = anharmonic_data%complex_modes, &
     & min_images      = anharmonic_min_images,         &
     & anharmonic_data = anharmonic_data                )
  
  do power=1,expansion_order/2
    ! Construct all translationally-invariant monomials of the given power
    !    using the given modes.
    monomials = construct_fine_monomials(power, modes)
    
    ! Interpolate the stress to find monomial coefficients,
    !    then convert the monomials into real basis functions.
    ! If power=paired_power, construct the u                 basis function.
    ! If power<paired_power, construct the u+u* and (u-u*)/i basis function.
    ! If power>paired_power, ignore the monomial, as it is included above.
    l = 0
    allocate( new_basis_functions(size(monomials)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(monomials)
      j = first( monomials(i)%powers()<monomials(i)%paired_powers(), &
               & default=size(monomials(i))+1                        )
      k = first( monomials(i)%powers()>monomials(i)%paired_powers(), &
               & default=size(monomials(i))+1                        )
      if (j==k) then
        coefficient = this%interpolate_coefficient( monomials(i), &
                                                  & interpolator  )
        l = l+1
        new_basis_functions(l) = BasisFunction( &
           & ComplexPolynomial([monomials(i)]), &
           & real(coefficient)                  )
      elseif (j<k) then
        coefficient = this%interpolate_coefficient( monomials(i), &
                                                  & interpolator  )
        l = l+1
        new_basis_functions(l) = BasisFunction(                      &
           & ComplexPolynomial([monomials(i), conjg(monomials(i))]), &
           & real(coefficient)                                       )
        l = l+1
        new_basis_functions(l) = BasisFunction(                         &
           & ComplexPolynomial( [monomials(i), -conjg(monomials(i))]    &
           &                  / cmplx(0.0_dp,1.0_dp,dp)              ), &
           & aimag(coefficient)                                         )
      endif
    enddo
    
    ! Append the new basis functions to the basis functions
    !    from previous powers.
    if (power==1) then
      basis_functions = new_basis_functions
    else
      old_basis_functions = basis_functions
      basis_functions = [old_basis_functions, new_basis_functions]
    endif
    deallocate(new_basis_functions, stat=ialloc); call err(ialloc)
  enddo
  
  ! Construct output.
  coupling_basis_functions = CouplingBasisFunctions(           &
     & coupling        = SubspaceCoupling(ids=[subspaces%id]), &
     & basis_functions = basis_functions                       )
  
  ! TODO: calculate reference energy correctly.
  output = PotentialPointer(PolynomialPotential(              &
     & potential_expansion_order = expansion_order,           &
     & reference_energy          = 0.0_dp,                    &
     & basis_functions           = [coupling_basis_functions] ))
end function

! Calculate the contribution to a given monomial from the interpolation of
!    this potential.
impure elemental function interpolate_coefficient_PolynomialPotential(this, &
   & monomial,interpolator) result(output)
  implicit none
  
  class(PolynomialPotential),   intent(in) :: this
  type(ComplexMonomial),        intent(in) :: monomial
  type(PolynomialInterpolator), intent(in) :: interpolator
  complex(dp)                              :: output
  
  output = sum(this%basis_functions_%interpolate_coefficient( monomial,    &
                                                            & interpolator ))
end function

! Calculate the effective dynamical matrices from which the potential can be
!    interpolated in the large-supercell limit.
function calculate_dynamical_matrices_PolynomialPotential(this,qpoints,       &
   & thermal_energy,subspaces,subspace_bases,subspace_states,anharmonic_data) &
   & result(output) 
  implicit none
  
  class(PolynomialPotential), intent(in)    :: this
  type(QpointData),           intent(in)    :: qpoints(:)
  real(dp),                   intent(in)    :: thermal_energy
  type(DegenerateSubspace),   intent(in)    :: subspaces(:)
  class(SubspaceBasis),       intent(in)    :: subspace_bases(:)
  class(BasisStates),         intent(inout) :: subspace_states(:)
  type(AnharmonicData),       intent(in)    :: anharmonic_data
  type(DynamicalMatrix), allocatable        :: output(:)
  
  integer :: i
  
  output = [( DynamicalMatrix(anharmonic_data%structure%no_atoms), &
            & i=1,                                                 &
            & size(qpoints)                                        )]
  
  do i=1,size(this%basis_functions_)
    output = output                                                 &
         & + this%basis_functions_(i)%calculate_dynamical_matrices( &
         &                                         qpoints,         &
         &                                         thermal_energy,  &
         &                                         subspaces,       &
         &                                         subspace_bases,  &
         &                                         subspace_states, &
         &                                         anharmonic_data  )
  enddo
end function

! Calculate the correction due to double counting
!    for the interpolated potential.
function energy_correction_PolynomialPotential(this,subspaces,subspace_bases, &
   & subspace_states,thermal_energy,anharmonic_data) result(output) 
  implicit none
  
  class(PolynomialPotential), intent(in)    :: this
  type(DegenerateSubspace),   intent(in)    :: subspaces(:)
  class(SubspaceBasis),       intent(in)    :: subspace_bases(:)
  class(BasisStates),         intent(inout) :: subspace_states(:)
  real(dp),                   intent(in)    :: thermal_energy
  type(AnharmonicData),       intent(in)    :: anharmonic_data
  real(dp)                                  :: output
  
  integer :: i
  
  output = 0
  do i=1,size(this%basis_functions_)
    output = output + this%basis_functions_(i)%energy_correction( &
                                               & subspaces,       &
                                               & subspace_bases,  &
                                               & subspace_states, &
                                               & thermal_energy,  &
                                               & anharmonic_data  )
  enddo
end function

! Expansion order.
impure elemental function expansion_order_PolynomialPotential(this) &
   & result(output)
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  integer                                :: output
  
  output = this%potential_expansion_order_
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
    output = [ 'Expansion order: '//this%potential_expansion_order_,      &
             & 'Reference energy: '//this%reference_energy_,              &
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
