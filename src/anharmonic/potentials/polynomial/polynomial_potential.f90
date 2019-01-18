! ======================================================================
! A polynomial representation of a potential.
! ======================================================================
module polynomial_potential_module
  use common_module
  
  use states_module
  use anharmonic_common_module
  
  use basis_function_module
  use coupling_basis_functions_module
  use vscf_rvectors_module
  use sampling_points_module
  use sample_result_module
  use sample_results_module
  use fit_coefficients_module
  implicit none
  
  private
  
  public :: startup_polynomial_potential
  
  public :: PolynomialPotential
  
  type, extends(PotentialData) :: PolynomialPotential
    integer,  private :: potential_expansion_order
    real(dp), private :: reference_energy
    type(CouplingBasisFunctions), allocatable, private :: basis_functions(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_PolynomialPotential
    
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
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_PolynomialPotential
    
    procedure, public :: iterate_damped => iterate_damped_PolynomialPotential
    procedure, public :: iterate_pulay => iterate_pulay_PolynomialPotential
    
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
  
  this%potential_expansion_order = potential_expansion_order
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
  this%basis_functions           = basis_functions
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
  
  ! Variables for processing electronic structure.
  type(StructureData)                    :: supercell
  type(VscfRvectors),        allocatable :: vscf_rvectors(:)
  type(ElectronicStructure), allocatable :: calculations(:)
  
  ! Previously calculated basis functions.
  type(CouplingBasisFunctions), allocatable :: basis_functions(:)
  
  ! Electronic structure results.
  type(SamplingPoints), allocatable :: sampling_points(:)
  type(CartesianDisplacement)       :: displacement
  type(SampleResults),  allocatable :: sample_results(:)
  
  ! Variables associated with the constant basis function
  !    and the sampling point at zero displacement.
  type(RealMonomial)          :: constant_real_monomial
  type(ComplexMonomial)       :: constant_complex_monomial
  type(BasisFunction)         :: constant_basis_function
  type(RealModeDisplacement)  :: equilibrium_sampling_point
  type(CartesianDisplacement) :: equilibrium_displacement
  type(ElectronicStructure)   :: equilibrium_electronic_structure
  type(SampleResult)          :: equilibrium_sample_result
  
  ! Variables for generating coefficients.
  logical,                    allocatable :: uncoupled(:)
  type(BasisFunction),        allocatable :: uncoupled_basis_functions(:)
  type(RealModeDisplacement), allocatable :: uncoupled_sampling_points(:)
  type(SampleResult),         allocatable :: uncoupled_sample_results(:)
  real(dp),                   allocatable :: coefficients(:)
  
  ! Directories and files.
  type(String) :: equilibrium_dir
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
  equilibrium_displacement = CartesianDisplacement( &
                      & equilibrium_sampling_point, &
                      & anharmonic_data%structure,  &
                      & anharmonic_data%real_modes, &
                      & anharmonic_data%qpoints     )
  
  ! --------------------------------------------------
  ! Read in basis functions and sampling points.
  ! --------------------------------------------------
  allocate( basis_functions(size(anharmonic_data%subspace_couplings)), &
          & sampling_points(size(anharmonic_data%subspace_couplings)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(sampling_points)
    coupling_dir = sampling_points_dir// &
       & '/coupling_'//left_pad(i,str(size(sampling_points)))
    
    ! Read in basis functions and sampling points.
    basis_functions_file = IFile(coupling_dir//'/basis_functions.dat')
    basis_functions(i) = CouplingBasisFunctions(basis_functions_file%lines())
    
    sampling_points_file = IFile(coupling_dir//'/sampling_points.dat')
    sampling_points(i) = SamplingPoints(sampling_points_file%lines())
  enddo
  
  ! --------------------------------------------------
  ! Read in energies and forces.
  ! --------------------------------------------------
  equilibrium_dir = sampling_points_dir//'/equilibrium'
  equilibrium_electronic_structure  = calculation_reader%read_calculation( &
                                                & equilibrium_dir,         &
                                                & equilibrium_displacement )
  equilibrium_sample_result = SampleResult( equilibrium_electronic_structure, &
                                          & anharmonic_data%structure,        &
                                          & anharmonic_data%real_modes,       &
                                          & anharmonic_data%qpoints           )
  
  allocate( sample_results(size(anharmonic_data%subspace_couplings)), &
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
      
      displacement = CartesianDisplacement( sampling_points(i)%points(j), &
                                          & supercell,                    &
                                          & anharmonic_data%real_modes,   &
                                          & anharmonic_data%qpoints       )
      
      vscf_rvectors_file = IFile(sampling_dir//'/vscf_rvectors.dat')
      vscf_rvectors = VscfRvectors(vscf_rvectors_file%sections())
      
      ! Read in electronic structure calculations.
      allocate( calculations(size(vscf_rvectors)), &
              & stat=ialloc); call err(ialloc)
      do k=1,size(vscf_rvectors)
        vscf_rvectors_dir = sampling_dir// &
           & '/vscf_rvectors_'//left_pad(k,str(size(vscf_rvectors)))
        calculations(k) = calculation_reader%read_calculation( &
                                          & vscf_rvectors_dir, &
                                          & displacement       )
      enddo
      
      ! Average electronic structure across VSCF R-vectors, and convert
      !    to correct normalisation and real mode co-ordinates.
      sample_results(i)%results(j) = SampleResult( &
                     & vscf_rvectors,              &
                     & calculations,               &
                     & supercell,                  &
                     & anharmonic_data%real_modes, &
                     & anharmonic_data%qpoints     )
      
      deallocate(calculations, stat=ialloc); call err(ialloc)
    enddo
  enddo
  
  ! --------------------------------------------------
  ! Calculate basis function coefficients.
  ! --------------------------------------------------
  
  ! Identify which subspace couplings only contain one subspace.
  allocate( uncoupled(size(anharmonic_data%subspace_couplings)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(anharmonic_data%subspace_couplings)
    if (size(anharmonic_data%subspace_couplings(i))==1) then
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
  do i=1,size(anharmonic_data%subspace_couplings)
    if (uncoupled(i)) then
      uncoupled_basis_functions = [ uncoupled_basis_functions,         &
                                  & basis_functions(i)%basis_functions ]
      uncoupled_sampling_points = [ uncoupled_sampling_points, &
                                  & sampling_points(i)%points ]
      uncoupled_sample_results = [ uncoupled_sample_results, &
                                 & sample_results(i)%results ]
    endif
  enddo
  
  coefficients = fit_coefficients( uncoupled_basis_functions,  &
                                 & uncoupled_sampling_points,  &
                                 & uncoupled_sample_results,   &
                                 & anharmonic_data%real_modes, &
                                 & weighted_energy_force_ratio )
  
  this%reference_energy = coefficients(1)
  allocate( this%basis_functions(count(uncoupled)), &
          & stat=ialloc); call err(ialloc)
  j = 0
  k = 1
  do i=1,size(anharmonic_data%subspace_couplings)
    if (uncoupled(i)) then
      j = j+1
      this%basis_functions(j) = CouplingBasisFunctions(                   &
         & coupling        = basis_functions(i)%coupling,                 &
         & basis_functions = coefficients(k+1:k+size(basis_functions(i))) &
         &                 * basis_functions(i)%basis_functions           )
      k = k+size(basis_functions(i))
    endif
  enddo
  if (k/=size(coefficients)) then
    call err()
  endif
  
  ! Calculate the coefficients of all basis functions involving subspace
  !    coupling. These are calculated on a coupling-by-coupling basis.
  do i=1,size(anharmonic_data%subspace_couplings)
    if (uncoupled(i)) then
      cycle
    endif
    
    coefficients = fit_coefficients( basis_functions(i)%basis_functions, &
                                   & sampling_points(i)%points,          &
                                   & sample_results(i)%results,          &
                                   & anharmonic_data%real_modes,         &
                                   & weighted_energy_force_ratio,        &
                                   & this                                )
    
    this%basis_functions = [                                       &
       & this%basis_functions,                                     &
       & CouplingBasisFunctions(                                   &
       &    coupling = basis_functions(i)%coupling,                &
       &    basis_functions = coefficients                         &
       &                    * basis_functions(i)%basis_functions ) ]
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
  
  integer :: i
  
  output = this%reference_energy
  do i=1,size(this%basis_functions)
    output = output &
         & + sum(this%basis_functions(i)%basis_functions%energy(displacement))
  enddo
end function

impure elemental function energy_ComplexModeDisplacement_PolynomialPotential( &
   & this,displacement) result(output)
  implicit none
  
  class(PolynomialPotential),    intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  complex(dp)                               :: output
  
  integer :: i
  
  output = this%reference_energy
  do i=1,size(this%basis_functions)
    output = output &
         & + sum(this%basis_functions(i)%basis_functions%energy(displacement))
  enddo
end function

! Calculate the force at a given displacement.
impure elemental function force_RealModeDisplacement_PolynomialPotential( &
   & this,displacement) result(output)
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  type(RealModeForce)                    :: output
  
  integer :: i
  
  output = sum([(                                                        &
     & sum(this%basis_functions(i)%basis_functions%force(displacement)), &
     & i=1,                                                              &
     & size(this%basis_functions)                                        )])
end function

impure elemental function force_ComplexModeDisplacement_PolynomialPotential( &
   & this,displacement) result(output)
  implicit none
  
  class(PolynomialPotential),    intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexModeForce)                    :: output
  
  integer :: i
  
  output = sum([(                                                        &
     & sum(this%basis_functions(i)%basis_functions%force(displacement)), &
     & i=1,                                                              &
     & size(this%basis_functions)                                        )])
end function

! Integrate the potential between two states.
function braket_PolynomialPotential(this,bra,ket,subspace,anharmonic_data) &
   & result(output)
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  class(SubspaceState),       intent(in) :: bra
  class(SubspaceState),       intent(in) :: ket
  type(DegenerateSubspace),   intent(in) :: subspace
  type(AnharmonicData),       intent(in) :: anharmonic_data
  class(PotentialData), allocatable      :: output
  
  logical, allocatable :: to_remove(:)
  
  integer :: i,j,k,l,ialloc
  
  allocate(output, source=this, stat=ialloc); call err(ialloc)
  select type(output); type is(PolynomialPotential)
    ! Integrate the reference energy (N.B. <i|e|j> = e<i|j> if e is a scalar.).
    output%reference_energy = output%reference_energy * braket(bra,ket)
    
    ! Integrate each basis function between the bra and the ket.
    to_remove = [(.false., i=1, size(output%basis_functions))]
    do i=1,size(output%basis_functions)
      j = first(output%basis_functions(i)%coupling%ids==subspace%id, default=0)
      if (j/=0) then
        do k=1,size(output%basis_functions(i))
          call output%basis_functions(i)%basis_functions(k)%braket( &
                                                  & bra,            &
                                                  & ket,            &
                                                  & anharmonic_data )
        enddo
        
        ! Simplify the output.
        call output%basis_functions(i)%basis_functions%simplify()
        
        ! Update the coupling to remove the integrated subspace.
        output%basis_functions(i)%coupling%ids = [         &
           & output%basis_functions(i)%coupling%ids(:j-1), &
           & output%basis_functions(j)%coupling%ids(j+1:)  ]
        
        ! Check if the basis function is now a constant.
        ! If so, add the constant energy to the potential's reference energy,
        !    and flag the term for removal.
        ! Then check if a coupling is now the same as a previous coupling.
        ! If so, combine the two and flag the duplicate term for removal.
        if (size(output%basis_functions(i)%coupling)==0) then
          output%reference_energy =      &
             &   output%reference_energy &
             & + sum(output%basis_functions(i)%basis_functions%undisplaced_energy())
          to_remove(i) = .true.
        else
          k = first([(                                          &
             & all(    output%basis_functions(i)%coupling%ids   &
             &      == output%basis_functions(l)%coupling%ids), &
             & l=1,                                             &
             & i-1                                              )], default=0)
          if (k/=0) then
            output%basis_functions(k)%basis_functions = [   &
               & output%basis_functions(k)%basis_functions, &
               & output%basis_functions(i)%basis_functions  ]
            to_remove(i) = .true.
          endif
        endif
      endif
    enddo
    
    ! Remove constant terms.
    output%basis_functions = output%basis_functions(filter(.not.to_remove))
  end select
end function

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
  
  integer :: i,j
  
  output = this%reference_energy
  do i=1,size(this%basis_functions)
    do j=1,size(this%basis_functions(i))
      output = output                                                        &
           & + this%basis_functions(i)%basis_functions(j)%harmonic_expectation( frequency,      &
           &                                                 thermal_energy, &
           &                                                 no_states,      &
           &                                                 subspace,       &
           &                                                 anharmonic_data )
    enddo
  enddo
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
  class(PotentialData), allocatable      :: output
  
  integer :: i,ialloc
  
  select type(new_potential); type is(PolynomialPotential)
    if (size(this%basis_functions)/=size(new_potential%basis_functions)) then
      call err()
    endif
    allocate(output, source=this, stat=ialloc); call err(ialloc)
    select type(output); type is(PolynomialPotential)
      output%reference_energy = (1-damping)*this%reference_energy &
                            & + damping*new_potential%reference_energy
      do i=1,size(output%basis_functions)
        if (    size(this%basis_functions(i))          &
           & /= size(new_potential%basis_functions(i)) ) then
          call err()
        endif
        call output%basis_functions(i)%basis_functions%set_coefficient(       &
           &   (1-damping)                                                    &
           & * this%basis_functions(i)%basis_functions%coefficient()          &
           & + damping                                                        &
           & * new_potential%basis_functions(i)%basis_functions%coefficient() )
      enddo
    end select
  class default
    call err()
  end select
end function

function iterate_pulay_PolynomialPotential(this,input_potentials, &
   & output_potentials,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  class(PotentialData),       intent(in) :: input_potentials(:)
  class(PotentialData),       intent(in) :: output_potentials(:)
  type(AnharmonicData),       intent(in) :: anharmonic_data
  class(PotentialData), allocatable      :: output
  
  type(RealVector), allocatable :: input_coefficients(:)
  type(RealVector), allocatable :: output_coefficients(:)
  real(dp),         allocatable :: pulay_coefficients(:)
  
  integer :: i,j,ialloc
  
  if (size(input_potentials)/=size(output_potentials)) then
    call err()
  endif
  
  select type(input_potentials); type is(PolynomialPotential)
    select type(output_potentials); type is(PolynomialPotential)
      allocate(output, source=this, stat=ialloc); call err(ialloc)
      select type(output); type is(PolynomialPotential)
        ! Convert previous potential iterations to vectors of coefficients.
        allocate( input_coefficients(size(input_potentials)),  &
                & output_coefficients(size(input_potentials)), &
                & stat=ialloc); call err(ialloc)
        do i=1,size(input_potentials)
            input_coefficients(i) = vec([(                 &
               & input_potentials(i)%basis_functions(j     &
               &          )%basis_functions%coefficient(), &
               & j=1,                                      &
               & size(input_potentials(i)%basis_functions) )])
            output_coefficients(i) = vec([(                &
               & output_potentials(i)%basis_functions(j    &
               &          )%basis_functions%coefficient(), &
               & j=1,                                      &
               & size(input_potentials(i)%basis_functions) )])
        enddo
        
        ! Use a Pulay scheme to construct the next set of coefficients.
        pulay_coefficients = dble(pulay( input_coefficients, &
                                       & output_coefficients ))
        
        ! Set the output coefficients from the Pulay coefficients.
        j = 0
        do i=1,size(output%basis_functions)
          call output%basis_functions(i)%basis_functions%set_coefficient( &
              & pulay_coefficients(j+1:j+size(output%basis_functions(i))) )
          j = j + size(output%basis_functions(i))
        enddo
        if (j/=size(pulay_coefficients)) then
          call err()
        endif
      end select
    class default
      call err()
    end select
  class default
    call err()
  end select
end function

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PART OF PULAY SCHEME.
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! TODO: convert this to a PotentialData function.
!function pulay_scheme()
!  implicit none
!
!  class(PotentialData),     intent(in) :: potential
!  type(DegenerateSubspace), intent(in) :: subspaces(:)
!  type(FullSubspaceBasis),  intent(in) :: subspace_bases(:)
!  real(dp),                 intent(in) :: energy_convergence
!  integer,                  intent(in) :: no_converged_calculations
!  integer,                  intent(in) :: max_pulay_iterations
!  integer,                  intent(in) :: pre_pulay_iterations
!  real(dp),                 intent(in) :: pre_pulay_damping
!  type(AnharmonicData),     intent(in) :: anharmonic_data
!  type(VscfOutput), allocatable        :: output(:)
!  
!  type(PotentialPointer),  allocatable :: subspace_potentials(:)
!  type(FullSubspaceState), allocatable :: subspace_states(:)
!  real(dp),                allocatable :: energies(:)
!  
!  type(RealVector), allocatable :: old_coefficients(:)
!  type(RealVector), allocatable :: new_coefficients(:)
!  
!  type(FullSubspaceState), allocatable :: input_states(:,:)
!  type(FullSubspaceState), allocatable :: output_states(:,:)
!  
!  type(RealVector) :: next_old_coefficients
!  
!  type(FractionVector), allocatable :: state_wavevectors(:)
!  integer                           :: wavevector_changed
!  
!  integer :: i,j,ialloc
!    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    new_coefficients = [ new_coefficients,                       &
!                       & states_to_coefficients(subspace_states) ]
!    
!    if (i==1) then
!      wavevector_changed = i
!    else
!      if (any(subspace_states%wavevector/=state_wavevectors)) then
!        wavevector_changed = i
!      endif
!    endif
!    
!    state_wavevectors = subspace_states%wavevector
!    
!    ! Use a damped iterative scheme or a Pulay scheme to converge towards the
!    !    self-consistent solution where new coefficients = old coefficients.
!    if (i-wavevector_changed<=pre_pulay_iterations) then
!      next_old_coefficients = (1-pre_pulay_damping) * old_coefficients(i) &
!                          & + pre_pulay_damping * new_coefficients(i)
!    else
!      j = max(2, i-max_pulay_iterations+1)
!      next_old_coefficients = pulay(old_coefficients(j:), new_coefficients(j:))
!    endif
!    
!    ! Normalise coefficients, and append them to the list.
!    next_old_coefficients = normalise_coefficients( next_old_coefficients, &
!                                                  & subspace_states        )
!    old_coefficients = [old_coefficients, next_old_coefficients]
!    
!    ! Use the next iteration of old_coefficients to update the current state.
!    ! N.B. the energy of this state is unknown.
!    subspace_states = coefficients_to_states( next_old_coefficients, &
!                                            & subspace_states        )
!    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!end function
! Concatenates the coefficients of all the states together,
!    to allow the Pulay scheme to be called.
!function states_to_coefficients(states) result(output)
!  implicit none
!  
!  type(FullSubspaceState), intent(in) :: states(:)
!  type(RealVector)                    :: output
!  
!  integer :: i
!  
!  output = vec([( states(i)%coefficients, i=1, size(states) )])
!end function
!
!! Reverses coefficients(), generating a set of states from a vector of
!!    coefficients.
!! Takes an array of states to act as a template for the output states.
!function coefficients_to_states(input,states) result(output)
!  implicit none
!  
!  type(RealVector),        intent(in)  :: input
!  type(FullSubspaceState), intent(in)  :: states(:)
!  type(FullSubspaceState), allocatable :: output(:)
!  
!  real(dp), allocatable :: coefficients(:)
!  real(dp), allocatable :: state_coefficients(:)
!  
!  integer :: i,j,ialloc
!  
!  coefficients = dble(input)
!  
!  j = 0
!  allocate(output(size(states)), stat=ialloc); call err(ialloc)
!  do i=1,size(states)
!    state_coefficients = coefficients(j+1:j+size(states(i)%coefficients))
!    j = j+size(state_coefficients)
!    output(i) = FullSubspaceState( subspace_id  = states(i)%subspace_id, &
!                                 & wavevector   = states(i)%wavevector,  &
!                                 & degeneracy   = states(i)%degeneracy,  &
!                                 & energy       = 0.0_dp,                &
!                                 & coefficients = state_coefficients     )
!  enddo
!end function
!
!! Normalises a set of coefficients, using a set of states to identify which
!!    coefficients correspond to which state.
!function normalise_coefficients(input,states) result(output)
!  implicit none
!  
!  type(RealVector),        intent(in) :: input
!  type(FullSubspaceState), intent(in) :: states(:)
!  type(RealVector)                    :: output
!  
!  real(dp), allocatable :: coefficients(:)
!  real(dp), allocatable :: state_coefficients(:)
!  real(dp), allocatable :: new_coefficients(:)
!  
!  integer :: i,j
!  
!  coefficients = dble(input)
!  new_coefficients = [real::]
!  do i=1,size(states)
!    j = size(new_coefficients)
!    state_coefficients = coefficients(j+1:j+size(states(i)%coefficients))
!    state_coefficients = state_coefficients / l2_norm(state_coefficients)
!    new_coefficients = [new_coefficients, state_coefficients]
!  enddo
!  output = new_coefficients
!end function
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    output = [ 'Expansion order: '//this%potential_expansion_order,      &
             & 'Reference energy: '//this%reference_energy,              &
             & str('Basis functions:'),                                  &
             & str(''),                                                  &
             & str(this%basis_functions, separating_line=repeat('=',50)) ]
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
