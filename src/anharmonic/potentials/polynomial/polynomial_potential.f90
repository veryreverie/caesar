! ======================================================================
! A polynomial representation of a potential.
! ======================================================================
module polynomial_potential_module
  use common_module
  
  use states_module
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
    real(dp),                         private :: reference_energy
    type(BasisFunction), allocatable, private :: basis_functions(:)
  contains
    procedure, public :: generate_sampling_points => &
       & generate_sampling_points_PolynomialPotential
    procedure, public :: generate_potential => &
       & generate_potential_PolynomialPotential
    
    procedure, public :: zero_energy => zero_energy_PolynomialPotential
    
    procedure, public :: energy_RealModeDisplacement => &
                       & energy_RealModeDisplacement_PolynomialPotential
    procedure, public :: energy_ComplexModeDisplacement => &
                       & energy_ComplexModeDisplacement_PolynomialPotential
    procedure, public :: force_RealModeDisplacement => &
                       & force_RealModeDisplacement_PolynomialPotential
    procedure, public :: force_ComplexModeDisplacement => &
                       & force_ComplexModeDisplacement_PolynomialPotential
    
    procedure, public :: braket => braket_PolynomialPotential
    
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

! Constructors.
function new_PolynomialPotential(potential_expansion_order) result(this)
  implicit none
  
  integer, intent(in)       :: potential_expansion_order
  type(PolynomialPotential) :: this
  
  this%potential_expansion_order = potential_expansion_order
end function

function new_PolynomialPotential_BasisFunctions(potential_expansion_order, &
   & reference_energy,basis_functions) result(this)
  implicit none
  
  integer,             intent(in) :: potential_expansion_order
  real(dp),            intent(in) :: reference_energy
  type(BasisFunction), intent(in) :: basis_functions(:)
  type(PolynomialPotential)       :: this
  
  this%potential_expansion_order = potential_expansion_order
  this%reference_energy          = reference_energy
  this%basis_functions           = basis_functions
end function

! Generate sampling points.
subroutine generate_sampling_points_PolynomialPotential(this,inputs, &
   & sampling_points_dir,calculation_writer,logfile)
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  type(AnharmonicData),       intent(in)    :: inputs
  type(String),               intent(in)    :: sampling_points_dir
  type(CalculationWriter),    intent(inout) :: calculation_writer
  type(OFile),                intent(inout) :: logfile
  
  ! Variables used when generating sampling points.
  type(SubspaceMonomial), allocatable :: subspace_monomials(:)
  type(BasisFunctions),   allocatable :: basis_functions(:)
  type(SamplingPoints),   allocatable :: sampling_points(:)
  
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
           & basis_functions(i),                   &
           & this%potential_expansion_order,       &
           & inputs%maximum_weighted_displacement, &
           & inputs%frequency_of_max_displacement, &
           & inputs%real_modes)
  enddo
  
  ! --------------------------------------------------
  ! Write out basis functions and sampling points.
  ! --------------------------------------------------
  
  ! Write out sampling point at equilibrium.
  equilibrium_dir = sampling_points_dir//'/equilibrium'
  call calculation_writer%write_calculation( inputs%structure, &
                                           & equilibrium_dir   )
  
  ! Write out all other sampling points.
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
      sampling_point_modes = select_modes( sampling_point%vectors, &
                                         & inputs%real_modes       )
      sampling_point_qpoints = select_qpoints( sampling_point_modes, &
                                             & inputs%qpoints        )
      supercell_matrix = construct_supercell_matrix( sampling_point_qpoints, &
                                                   & inputs%structure        )
      supercell = construct_supercell( inputs%structure, &
                                     & supercell_matrix  )
      
      ! Write out the supercell.
      supercell_file = OFile(sampling_dir//'/structure.dat')
      call supercell_file%print_lines(supercell)
      
      ! Construct VSCF R-vectors.
      vscf_rvectors = construct_vscf_rvectors( sampling_point,    &
                                             & supercell,         &
                                             & inputs%real_modes, &
                                             & inputs%qpoints     )
      vscf_rvectors_file = OFile(sampling_dir//'/vscf_rvectors.dat')
      call vscf_rvectors_file%print_lines(vscf_rvectors,separating_line='')
      
      do k=1,size(vscf_rvectors)
        ! Transform the sampling point by the VSCF R-vector.
        transformed_sampling_point = vscf_rvectors(k)%transform( &
                                            & sampling_point,    &
                                            & inputs%real_modes, &
                                            & inputs%qpoints     )
        
        ! Construct displaced structure.
        displacement = CartesianDisplacement( transformed_sampling_point, &
                                            & supercell,                  &
                                            & inputs%real_modes,          &
                                            & inputs%qpoints              )
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
subroutine generate_potential_PolynomialPotential(this,inputs,           &
   & weighted_energy_force_ratio,sampling_points_dir,calculation_reader, &
   & logfile)
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  type(AnharmonicData),       intent(in)    :: inputs
  real(dp),                   intent(in)    :: weighted_energy_force_ratio
  type(String),               intent(in)    :: sampling_points_dir
  type(CalculationReader),    intent(inout) :: calculation_reader
  type(OFile),                intent(inout) :: logfile
  
  ! Variables for processing electronic structure.
  type(StructureData)                    :: supercell
  type(VscfRvectors),        allocatable :: vscf_rvectors(:)
  type(ElectronicStructure), allocatable :: calculations(:)
  
  ! Previously calculated basis functions.
  type(BasisFunctions), allocatable :: basis_functions(:)
  
  ! Electronic structure results.
  type(SamplingPoints), allocatable :: sampling_points(:)
  type(SampleResults),  allocatable :: sample_results(:)
  
  ! Variables associated with the constant basis function
  !    and the sampling point at zero displacement.
  type(RealMonomial)         :: constant_real_monomial
  type(ComplexMonomial)      :: constant_complex_monomial
  type(BasisFunction)        :: constant_basis_function
  type(RealModeDisplacement) :: equilibrium_sampling_point
  type(String)               :: equilibrium_dir
  type(ElectronicStructure)  :: equilibrium_electronic_structure
  type(SampleResult)         :: equilibrium_sample_result
  
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
  type(IFile)  :: supercell_file
  type(IFile)  :: vscf_rvectors_file
  type(String) :: vscf_rvectors_dir
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Construct a basis function representing a constant reference energy,
  !    and a sampling point representing the unperturbed structure.
  ! --------------------------------------------------
  constant_real_monomial = RealMonomial( coefficient=1.0_dp, &
                                       & modes=[RealUnivariate::])
  constant_complex_monomial = ComplexMonomial( &
        & coefficient=cmplx(1.0_dp,0.0_dp,dp), &
        & modes=[ComplexUnivariate::])
  constant_basis_function = BasisFunction(                             &
     & real_representation = RealPolynomial([constant_real_monomial]), &
     & complex_representation =                                        &
     &    ComplexPolynomial([constant_complex_monomial])               )
  
  equilibrium_sampling_point = RealModeDisplacement([RealSingleDisplacement::])
  
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
    basis_functions(i) = BasisFunctions(basis_functions_file%lines())
    
    sampling_points_file = IFile(coupling_dir//'/sampling_points.dat')
    sampling_points(i) = SamplingPoints(sampling_points_file%lines())
  enddo
  
  ! --------------------------------------------------
  ! Read in energies and forces.
  ! --------------------------------------------------
  equilibrium_dir = sampling_points_dir//'/equilibrium'
  equilibrium_electronic_structure  = calculation_reader%read_calculation( &
                                                         & equilibrium_dir )
  equilibrium_sample_result = SampleResult( equilibrium_electronic_structure, &
                                          & inputs%structure,                 &
                                          & inputs%real_modes,                &
                                          & inputs%qpoints                    )
  
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
      supercell_file = IFile(sampling_dir//'/structure.dat')
      supercell = StructureData(supercell_file%lines())
      
      vscf_rvectors_file = IFile(sampling_dir//'/vscf_rvectors.dat')
      vscf_rvectors = VscfRvectors(vscf_rvectors_file%sections())
      
      ! Read in electronic structure calculations.
      allocate( calculations(size(vscf_rvectors)), &
              & stat=ialloc); call err(ialloc)
      do k=1,size(vscf_rvectors)
        vscf_rvectors_dir = sampling_dir// &
           & '/vscf_rvectors_'//left_pad(k,str(size(vscf_rvectors)))
        calculations(k) = calculation_reader%read_calculation( &
           & vscf_rvectors_dir)
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
  !    subspace coupling at once.
  
  uncoupled_basis_functions = [constant_basis_function]
  uncoupled_sampling_points = [equilibrium_sampling_point]
  uncoupled_sample_results  = [equilibrium_sample_result]
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
  
  this%reference_energy = coefficients(1)
  this%basis_functions  = coefficients(2:) * uncoupled_basis_functions(2:)
  
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
                                   & this                          )
    
    this%basis_functions = [ this%basis_functions,                       &
                           & coefficients * basis_functions(i)%functions ]
  enddo
end subroutine

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
       & + sum(this%basis_functions%energy(displacement))
end function

impure elemental function energy_ComplexModeDisplacement_PolynomialPotential( &
   & this,displacement) result(output)
  implicit none
  
  class(PolynomialPotential),    intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  complex(dp)                               :: output
  
  output = this%reference_energy &
       & + sum(this%basis_functions%energy(displacement))
end function

! Calculate the force at a given displacement.
impure elemental function force_RealModeDisplacement_PolynomialPotential( &
   & this,displacement) result(output)
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  type(RealModeForce)                    :: output
  
  output = sum(this%basis_functions%force(displacement))
end function

impure elemental function force_ComplexModeDisplacement_PolynomialPotential( &
   & this,displacement) result(output)
  implicit none
  
  class(PolynomialPotential),    intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexModeForce)                    :: output
  
  output = sum(this%basis_functions%force(displacement))
end function

! Integrate the potential between two states.
subroutine braket_PolynomialPotential(this,bra,ket,inputs)
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  class(SubspaceState),       intent(in)    :: bra
  class(SubspaceState),       intent(in)    :: ket
  type(AnharmonicData),       intent(in)    :: inputs
  
  logical,             allocatable :: is_constant(:)
  type(BasisFunction), allocatable :: constant_terms(:)
  
  integer :: i
  
  ! Integrate the reference energy and each basis function
  !    between the bra and the ket.
  this%reference_energy = this%reference_energy * braket(bra,ket)
  
  do i=1,size(this%basis_functions)
    call this%basis_functions(i)%braket(bra,ket,inputs)
  enddo
  
  ! Simplify the output.
  call this%basis_functions%simplify()
  
  is_constant = this%basis_functions%is_constant()
  constant_terms = this%basis_functions(filter(is_constant))
  this%basis_functions = this%basis_functions(filter(.not.is_constant))
  this%reference_energy = this%reference_energy &
                      & + sum(constant_terms%undisplaced_energy())
end subroutine

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_PolynomialPotential(this,input)
  implicit none
  
  class(PolynomialPotential), intent(out) :: this
  type(String),               intent(in)  :: input(:)
  
  integer                          :: potential_expansion_order
  real(dp)                         :: reference_energy
  type(BasisFunction), allocatable :: basis_functions(:)
  
  type(String),      allocatable :: line(:)
  type(StringArray), allocatable :: sections(:)
  
  select type(this); type is(PolynomialPotential)
    line = split_line(input(1))
    potential_expansion_order = int(line(3))
    
    line = split_line(input(2))
    reference_energy = dble(line(3))
    
    sections = split_into_sections(input)
    basis_functions = BasisFunction(sections(2:))
    
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
    output = [ 'Expansion order: '//this%potential_expansion_order, &
             & 'Reference energy: '//this%reference_energy,         &
             & str('Basis functions:'),                             &
             & str(''),                                             &
             & str(this%basis_functions, separating_line='')        ]
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
