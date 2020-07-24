! ======================================================================
! A set of basis functions spanning a given subspace coupling.
! ======================================================================
module coupling_basis_functions_module
  use common_module
  
  use anharmonic_common_module
  
  use polynomial_interpolator_module
  use basis_function_module
  use sampling_points_module
  use sample_result_module
  implicit none
  
  private
  
  public :: CouplingBasisFunctions
  public :: BasisFunctionsAndSamplingPoints
  public :: size
  public :: generate_basis_functions
  
  type, extends(PotentialBase) :: CouplingBasisFunctions
    type(SubspaceCoupling)                    :: coupling
    type(BasisFunction), allocatable, private :: basis_functions_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_CouplingBasisFunctions
    
    procedure, public :: basis_functions => &
                       & basis_functions_CouplingBasisFunctions
    
    procedure, public :: undisplaced_energy => &
                       & undisplaced_energy_CouplingBasisFunctions
    
    procedure, public :: zero_energy => zero_energy_CouplingBasisFunctions
    procedure, public :: add_constant => add_constant_CouplingBasisFunctions
    
    procedure, public :: energy_RealModeDisplacement => &
                       & energy_RealModeDisplacement_CouplingBasisFunctions
    procedure, public :: energy_ComplexModeDisplacement => &
                       & energy_ComplexModeDisplacement_CouplingBasisFunctions
    
    procedure, public :: force_RealModeDisplacement => &
                       & force_RealModeDisplacement_CouplingBasisFunctions
    procedure, public :: force_ComplexModeDisplacement => &
                       & force_ComplexModeDisplacement_CouplingBasisFunctions
    
    procedure, public :: braket_SubspaceBraKet => &
                       & braket_SubspaceBraKet_CouplingBasisFunctions
    procedure, public :: braket_BasisState => &
                       & braket_BasisState_CouplingBasisFunctions
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_CouplingBasisFunctions
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_CouplingBasisFunctions
    
    procedure, public :: potential_energy_SubspaceBraKet => &
                       & potential_energy_SubspaceBraKet_CouplingBasisFunctions
    procedure, public :: potential_energy_BasisState => &
                       & potential_energy_BasisState_CouplingBasisFunctions
    
    procedure, public :: fit_coefficients => &
                       & fit_coefficients_CouplingBasisFunctions
    procedure, public :: no_coefficients => &
                       & no_coefficients_CouplingBasisFunctions
    procedure, public :: coefficients => &
                       & coefficients_CouplingBasisFunctions
    procedure, public :: set_coefficients => &
                       & set_coefficients_CouplingBasisFunctions
    procedure, public :: all_basis_functions => &
                       & all_basis_functions_CouplingBasisFunctions
    procedure, public :: variable_basis_functions => &
                       & variable_basis_functions_CouplingBasisFunctions
    
    procedure, public :: optimise => optimise_CouplingBasisFunctions
    
    procedure, public :: append => append_CouplingBasisFunctions
    
    procedure, public :: interpolate_coefficient => &
                       & interpolate_coefficient_CouplingBasisFunctions
    procedure, public :: calculate_dynamical_matrices => &
                       & calculate_dynamical_matrices_CouplingBasisFunctions
    procedure, public :: energy_correction => &
                       & energy_correction_CouplingBasisFunctions
    
    ! I/O.
    procedure, public :: read  => read_CouplingBasisFunctions
    procedure, public :: write => write_CouplingBasisFunctions
  end type
  
  interface CouplingBasisFunctions
    module procedure new_CouplingBasisFunctions_empty
    module procedure new_CouplingBasisFunctions
    module procedure new_CouplingBasisFunctions_Strings
    module procedure new_CouplingBasisFunctions_StringArray
  end interface
  
  interface size
    module procedure size_CouplingBasisFunctions
  end interface
  
  interface generate_basis_functions
    module procedure generate_basis_functions_SubspaceCoupling
  end interface
  
  type :: BasisFunctionsAndSamplingPoints
    type(CouplingBasisFunctions) :: basis_functions
    type(SamplingPoints)         :: sampling_points
  end type
contains

impure elemental function new_CouplingBasisFunctions_empty(coupling) &
   & result(this) 
  implicit none
  
  type(SubspaceCoupling), intent(in) :: coupling
  type(CouplingBasisFunctions)       :: this
  
  integer :: ialloc
  
  this%coupling = coupling
  allocate(this%basis_functions_(0), stat=ialloc); call err(ialloc)
end function

function new_CouplingBasisFunctions(coupling,basis_functions) result(this)
  implicit none
  
  type(SubspaceCoupling), intent(in) :: coupling
  type(BasisFunction),    intent(in) :: basis_functions(:)
  type(CouplingBasisFunctions)       :: this
  
  this%coupling         = coupling
  this%basis_functions_ = basis_functions
end function

impure elemental function representation_CouplingBasisFunctions() &
   & result(output)
  implicit none
  
  type(String) :: output
  
  output = 'Polynomial basis functions'
end function

function size_CouplingBasisFunctions(this) result(output)
  implicit none
  
  type(CouplingBasisFunctions), intent(in) :: this
  integer                                  :: output
  
  output = size(this%basis_functions_)
end function

function basis_functions_CouplingBasisFunctions(this) result(output)
  implicit none
  
  class(CouplingBasisFunctions), intent(in) :: this
  type(BasisFunction), allocatable          :: output(:)
  
  output = this%basis_functions_
end function

impure elemental subroutine optimise_CouplingBasisFunctions(this,subspace, &
   & subspace_basis,old_subspace_potential,anharmonic_data)
  implicit none
  
  class(CouplingBasisFunctions), intent(inout)        :: this
  type(DegenerateSubspace),      intent(in)           :: subspace
  class(SubspaceBasis),          intent(in)           :: subspace_basis
  class(PotentialData),          intent(in), optional :: old_subspace_potential
  type(AnharmonicData),          intent(in)           :: anharmonic_data
  
  if (size(this%coupling%ids)/=1) then
    call print_line(CODE_ERROR//': Calling optimise_subspace_potential &
       &on a potential with coupled subspaces.')
    call err()
  elseif (this%coupling%ids(1)/=subspace%id) then
    call print_line(CODE_ERROR//': Calling optimise_subspace_potential &
       &with the wrong subspace.')
    call err()
  endif
  
  ! Remove constant terms and split basis functions by power.
  this%basis_functions_ = optimise( this%basis_functions_,  &
                                  & subspace,               &
                                  & subspace_basis,         &
                                  & old_subspace_potential, &
                                  & anharmonic_data         )
end subroutine

impure elemental function energy_RealModeDisplacement_CouplingBasisFunctions( &
   & this,displacement) result(output)
  implicit none
  
  class(CouplingBasisFunctions), intent(in) :: this
  type(RealModeDisplacement),    intent(in) :: displacement
  real(dp)                                  :: output
  
  output = sum(this%basis_functions_%energy(displacement))
end function

impure elemental function                                                     &
   & energy_ComplexModeDisplacement_CouplingBasisFunctions(this,displacement) &
   & result(output)
  implicit none
  
  class(CouplingBasisFunctions), intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  complex(dp)                               :: output
  
  output = sum(this%basis_functions_%energy(displacement))
end function

impure elemental function force_RealModeDisplacement_CouplingBasisFunctions( &
   & this,displacement) result(output)
  implicit none
  
  class(CouplingBasisFunctions), intent(in) :: this
  type(RealModeDisplacement),    intent(in) :: displacement
  type(RealModeForce)                       :: output
  
  if (size(this%basis_functions_)>0) then
    output = sum(this%basis_functions_%force(displacement))
  else
    output = RealModeForce(RealSingleForce(displacement%vectors%id,0.0_dp))
  endif
end function

impure elemental function                                                    &
   & force_ComplexModeDisplacement_CouplingBasisFunctions(this,displacement) &
   & result(output)
  implicit none
  
  class(CouplingBasisFunctions), intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexModeForce)                    :: output
  
  if (size(this%basis_functions_)>0) then
    output = sum(this%basis_functions_%force(displacement))
  else
    output = ComplexModeForce(ComplexSingleForce( displacement%vectors%id, &
                                                & (0.0_dp,0.0_dp)          ))
  endif
end function

impure elemental subroutine braket_SubspaceBraKet_CouplingBasisFunctions( &
   & this,braket,whole_subspace,anharmonic_data)
  implicit none
  
  class(CouplingBasisFunctions), intent(inout)        :: this
  class(SubspaceBraKet),         intent(in)           :: braket
  logical,                       intent(in), optional :: whole_subspace
  type(AnharmonicData),          intent(in)           :: anharmonic_data
  
  integer :: i
  
  ! Check if the subspace is in this basis function's coupling.
  i = first(this%coupling%ids==braket%subspace_id, default=0)
  if (i/=0) then
    ! If whole_subspace is .true., remove the subspace from the coupling.
    if (set_default(whole_subspace,.true.)) then
      this%coupling%ids = [this%coupling%ids(:i-1),this%coupling%ids(i+1:)]
    endif
    
    ! Integrate across the basis function, and simplify it.
    do i=1,size(this)
      call this%basis_functions_(i)%braket( braket,                           &
                                          & anharmonic_data = anharmonic_data )
    enddo
    
    call this%basis_functions_%simplify()
  endif
end subroutine

impure elemental subroutine braket_BasisState_CouplingBasisFunctions(this, &
   & bra,ket,subspace,subspace_basis,whole_subspace,anharmonic_data)
  implicit none
  
  class(CouplingBasisFunctions), intent(inout)        :: this
  class(BasisState),             intent(in)           :: bra
  class(BasisState),             intent(in), optional :: ket
  type(DegenerateSubspace),      intent(in)           :: subspace
  class(SubspaceBasis),          intent(in)           :: subspace_basis
  logical,                       intent(in), optional :: whole_subspace
  type(AnharmonicData),          intent(in)           :: anharmonic_data
  
  integer :: i
  
  ! Check if the subspace is in this basis function's coupling.
  i = first(this%coupling%ids==bra%subspace_id, default=0)
  if (i/=0) then
    ! If whole_subspace is .true., remove the subspace from the coupling.
    if (set_default(whole_subspace,.true.)) then
      this%coupling%ids = [this%coupling%ids(:i-1),this%coupling%ids(i+1:)]
    endif
    
    ! Integrate across the basis function, and simplify it.
    do i=1,size(this)
      call this%basis_functions_(i)%braket( bra,                              &
                                          & ket,                              &
                                          & subspace,                         &
                                          & subspace_basis,                   &
                                          & anharmonic_data = anharmonic_data )
    enddo
    
    call this%basis_functions_%simplify()
  endif
end subroutine

impure elemental subroutine braket_BasisStates_CouplingBasisFunctions(this, &
   & states,subspace,subspace_basis,whole_subspace,anharmonic_data)
  implicit none
  
  class(CouplingBasisFunctions), intent(inout)        :: this
  class(BasisStates),            intent(inout)        :: states
  type(DegenerateSubspace),      intent(in)           :: subspace
  class(SubspaceBasis),          intent(in)           :: subspace_basis
  logical,                       intent(in), optional :: whole_subspace
  type(AnharmonicData),          intent(in)           :: anharmonic_data
  
  integer :: i,j
  
  ! Check if the subspace is in this basis function's coupling.
  i = first(this%coupling%ids==states%subspace_id, default=0)
  if (i/=0) then
    ! If whole_subspace is .true., remove the subspace from the coupling.
    if (set_default(whole_subspace,.true.)) then
      this%coupling%ids = [this%coupling%ids(:i-1),this%coupling%ids(i+1:)]
    endif
    
    ! Integrate across the basis function, and simplify it.
    do j=1,size(this)
      call this%basis_functions_(j)%braket( states,                           &
                                          & subspace,                         &
                                          & subspace_basis,                   &
                                          & anharmonic_data = anharmonic_data )
    enddo
    
    call this%basis_functions_%simplify()
  endif
end subroutine

impure elemental function harmonic_expectation_CouplingBasisFunctions(this, &
   & frequency,thermal_energy,supercell_size,anharmonic_data) result(output)
  implicit none
  
  class(CouplingBasisFunctions), intent(in) :: this
  real(dp),                      intent(in) :: frequency
  real(dp),                      intent(in) :: thermal_energy
  integer,                       intent(in) :: supercell_size
  type(AnharmonicData),          intent(in) :: anharmonic_data
  real(dp)                                  :: output
  
  integer :: i
  
  output = sum([( this%basis_functions_(i)%harmonic_expectation(    &
                &                                frequency,         &
                &                                thermal_energy,    &
                &                                supercell_size,    &
                &                                anharmonic_data ), &
                & i=1,                                              &
                & size(this%basis_functions_)                       )])
end function

function potential_energy_SubspaceBraKet_CouplingBasisFunctions(this,braket, &
   & anharmonic_data) result(output) 
  implicit none
  
  class(CouplingBasisFunctions), intent(in) :: this
  class(SubspaceBraKet),         intent(in) :: braket
  type(AnharmonicData),          intent(in) :: anharmonic_data
  real(dp)                                  :: output
  
  integer :: i
  
  output = sum([( this%basis_functions_(i)%potential_energy(    &
                &                            braket,            &
                &                            anharmonic_data ), &
                & i=1,                                          &
                & size(this%basis_functions_)                   )])
end function

function potential_energy_BasisState_CouplingBasisFunctions(this,bra,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output) 
  implicit none
  
  class(CouplingBasisFunctions), intent(in)           :: this
  class(BasisState),             intent(in)           :: bra
  class(BasisState),             intent(in), optional :: ket
  type(DegenerateSubspace),      intent(in)           :: subspace
  class(SubspaceBasis),          intent(in)           :: subspace_basis
  type(AnharmonicData),          intent(in)           :: anharmonic_data
  real(dp)                                            :: output
  
  integer :: i
  
  output = sum([( this%basis_functions_(i)%potential_energy(    &
                &                            bra,               &
                &                            ket,               &
                &                            subspace,          &
                &                            subspace_basis,    &
                &                            anharmonic_data ), &
                & i=1,                                          &
                & size(this%basis_functions_)                   )])
end function

function generate_basis_functions_SubspaceCoupling(coupling,               &
   & potential_expansion_order,structure,complex_modes,real_modes,qpoints, &
   & subspaces,degenerate_symmetries,vscf_basis_functions_only,            &
   & maximum_weighted_displacement,frequency_of_max_displacement,          &
   & energy_to_force_ratio,logfile) result(output)
  implicit none
  
  type(SubspaceCoupling),   intent(in)    :: coupling
  integer,                  intent(in)    :: potential_expansion_order
  type(StructureData),      intent(in)    :: structure
  type(ComplexMode),        intent(in)    :: complex_modes(:)
  type(RealMode),           intent(in)    :: real_modes(:)
  type(QpointData),         intent(in)    :: qpoints(:)
  type(DegenerateSubspace), intent(in)    :: subspaces(:)
  type(DegenerateSymmetry), intent(in)    :: degenerate_symmetries(:)
  logical,                  intent(in)    :: vscf_basis_functions_only
  real(dp),                 intent(in)    :: maximum_weighted_displacement
  real(dp),                 intent(in)    :: frequency_of_max_displacement
  real(dp),                 intent(in)    :: energy_to_force_ratio
  type(OFile),              intent(inout) :: logfile
  type(BasisFunctionsAndSamplingPoints)   :: output
  
  type(SubspaceMonomial), allocatable :: subspace_monomials(:)
  
  type(BasisFunctions), allocatable :: basis_functions(:)
  
  type(BasisFunction), allocatable :: coupling_basis_functions(:)
  
  integer :: i,ialloc
  
  ! Generate the set of subspace monomials corresponding to the subspace
  !    coupling.
  ! e.g. the coupling [1,2] might have monomials [1,2], [1,1,2] and [1,2,2].
  subspace_monomials = generate_subspace_monomials(        &
     & coupling,                                           &
     & subspaces,                                          &
     & minimum_expansion_order = 2,                        &
     & maximum_expansion_order = potential_expansion_order )
    
  ! Loop over the subspace monomials corresponding to the coupling.
  allocate( basis_functions(size(subspace_monomials)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(subspace_monomials)
    ! Generate all basis functions for the subspace monomial.
    basis_functions(i) = generate_basis_functions( subspace_monomials(i),     &
                                                 & structure,                 &
                                                 & complex_modes,             &
                                                 & qpoints,                   &
                                                 & subspaces,                 &
                                                 & degenerate_symmetries,     &
                                                 & vscf_basis_functions_only, &
                                                 & logfile                    )
  enddo
  
  ! Concatenate the terms from each subspace monomial together.
  coupling_basis_functions = [( basis_functions(i)%basis_functions, &
                              & i=1,                                &
                              & size(subspace_monomials)            )]
  
  output = BasisFunctionsAndSamplingPoints(                                 &
     & basis_functions = CouplingBasisFunctions( coupling,                  &
     &                                           coupling_basis_functions), &
     & sampling_points = generate_sampling_points(                          &
     &                                     basis_functions,                 &
     &                                     potential_expansion_order,       &
     &                                     maximum_weighted_displacement,   &
     &                                     frequency_of_max_displacement,   &
     &                                     real_modes,                      &
     &                                     energy_to_force_ratio          ) )
end function

! Return the energy at zero displacement.
impure elemental function undisplaced_energy_CouplingBasisFunctions(this) &
   & result(output)
  implicit none
  
  class(CouplingBasisFunctions), intent(in) :: this
  real(dp)                                  :: output
  
  type(RealSingleDisplacement) :: zero_displacement(0)
  
  output = this%energy(RealModeDisplacement(zero_displacement))
end function

impure elemental subroutine zero_energy_CouplingBasisFunctions(this)
  implicit none
  
  class(CouplingBasisFunctions), intent(inout) :: this
  
  call print_line(CODE_ERROR//': zero_energy() cannot be called for type &
     &CouplingBasisFunctions.')
  call err()
end subroutine

impure elemental subroutine add_constant_CouplingBasisFunctions(this,input)
  implicit none
  
  class(CouplingBasisFunctions), intent(inout) :: this
  real(dp),                      intent(in)    :: input
  
  call print_line(CODE_ERROR//': add_constant() cannot be called for type &
     &CouplingBasisFunctions.')
  call err()
end subroutine

! Fit basis function coefficients using L2 regression.
subroutine fit_coefficients_CouplingBasisFunctions(this,sampling_points, &
   & sample_results,modes,energy_force_ratio,potential)
  implicit none
  
  class(CouplingBasisFunctions), intent(inout)          :: this
  type(RealModeDisplacement),    intent(in)             :: sampling_points(:)
  type(SampleResult),            intent(in)             :: sample_results(:)
  type(RealMode),                intent(in)             :: modes(:)
  real(dp),                      intent(in)             :: energy_force_ratio
  class(PotentialData),          intent(in),   optional :: potential
  
  integer               :: dimensions
  real(dp), allocatable :: a(:,:)
  real(dp), allocatable :: b(:)
  
  real(dp), allocatable :: coefficients(:)
  
  ! Check inputs are consistent.
  if (size(sampling_points)/=size(sample_results)) then
    call print_line(CODE_ERROR//': The number of sampling points does not &
       &match the number of results.')
    call err()
  endif
  
  ! Each calculation yields size(modes) forces and one energy.
  dimensions = 1+size(modes)
  
  ! Calculate the energies and forces due to each basis function at each
  !    sampling point.
  a = construct_sample_matrix( this%basis_functions_, &
                             & sampling_points,       &
                             & modes,                 &
                             & energy_force_ratio     )
  
  ! Calculate the energies and forces sampled at each sampling point.
  b = construct_sample_vector( sampling_points,   &
                             & sample_results,    &
                             & potential,         &
                             & modes,             &
                             & energy_force_ratio )
  
  ! Run linear least squares to get the basis function coefficients.
  ! This finds x s.t. (a.x-b)^2 is minimised.
  coefficients = dble(linear_least_squares(a, b))
  
  this%basis_functions_ = coefficients * this%basis_functions_
end subroutine

! Get and set basis function coefficients.
impure elemental function no_coefficients_CouplingBasisFunctions(this, &
   & maximum_power) result(output)
  implicit none
  
  class(CouplingBasisFunctions), intent(in) :: this
  integer,                       intent(in) :: maximum_power
  integer                                   :: output
  
  output = count(this%basis_functions_%power()<=maximum_power)
end function

function coefficients_CouplingBasisFunctions(this,maximum_power) result(output)
  implicit none
  
  class(CouplingBasisFunctions), intent(in) :: this
  integer,                       intent(in) :: maximum_power
  real(dp), allocatable                     :: output(:)
  
  output = this%basis_functions_(                           &
     & filter(this%basis_functions_%power()<=maximum_power) )%coefficient()
end function

subroutine set_coefficients_CouplingBasisFunctions(this,coefficients, &
   & maximum_power) 
  implicit none
  
  class(CouplingBasisFunctions), intent(inout) :: this
  integer,                       intent(in)    :: maximum_power
  real(dp),                      intent(in)    :: coefficients(:)
  
  integer :: i,j
  
  j = 0
  do i=1,size(this%basis_functions_)
    if (this%basis_functions_(i)%power()<=maximum_power) then
      j = j+1
      call this%basis_functions_(i)%set_coefficient(coefficients(j))
    endif
  enddo
  
  if (j/=size(coefficients)) then
    call print_line(CODE_ERROR//': Incorrect number of coefficients.')
    call err()
  endif
end subroutine

function all_basis_functions_CouplingBasisFunctions(this) result(output) 
  implicit none
  
  class(CouplingBasisFunctions), intent(in) :: this
  type(PotentialBasePointer), allocatable   :: output(:)
  
  type(BasisFunction), allocatable :: basis_functions(:)
  
  basis_functions = this%basis_functions_
  call basis_functions%set_coefficient(1.0_dp)
  
  output = PotentialBasePointer(basis_functions)
end function

function variable_basis_functions_CouplingBasisFunctions(this,maximum_power) &
   & result(output) 
  implicit none
  
  class(CouplingBasisFunctions), intent(in) :: this
  integer,                       intent(in) :: maximum_power
  type(PotentialBasePointer), allocatable   :: output(:)
  
  type(BasisFunction), allocatable :: basis_functions(:)
  
  basis_functions = this%basis_functions_(                  &
     & filter(this%basis_functions_%power()<=maximum_power) )
  call basis_functions%set_coefficient(1.0_dp)
  
  output = PotentialBasePointer(basis_functions)
end function

! Append another CouplingBasisFunctions to this.
subroutine append_CouplingBasisFunctions(this,that)
  implicit none
  
  class(CouplingBasisFunctions), intent(inout) :: this
  type(CouplingBasisFunctions),  intent(in)    :: that
  
  type(BasisFunction), allocatable :: temp(:)
  
  if (this%coupling/=that%coupling) then
    call print_line(CODE_ERROR//': Appending incompatible basis functions.')
    call err()
  endif
  
  ! WORKAROUND: combining the following two lines makes ifort crash,
  !    as of version 19.0.4.227.
  temp = [this%basis_functions_, that%basis_functions_]
  this%basis_functions_ = temp
end subroutine

! Calculate the contribution to a given monomial from the interpolation of
!    this basis function.
impure elemental function interpolate_coefficient_CouplingBasisFunctions( &
   & this,monomial,interpolator) result(output)
  implicit none
  
  class(CouplingBasisFunctions), intent(in) :: this
  type(ComplexMonomial),         intent(in) :: monomial
  type(PolynomialInterpolator),  intent(in) :: interpolator
  complex(dp)                               :: output
  
  output = sum(this%basis_functions_%interpolate_coefficient( monomial,    &
                                                            & interpolator ))
end function

! Calculate this basis function's contribution to the effective dynamical
!    matrix from which the potential can be interpolated in the large-supercell
!    limit.
function calculate_dynamical_matrices_CouplingBasisFunctions(this,qpoints,    &
   & thermal_energy,subspaces,subspace_bases,subspace_states,anharmonic_data) &
   & result(output) 
  implicit none
  
  class(CouplingBasisFunctions), intent(in)    :: this
  type(QpointData),              intent(in)    :: qpoints(:)
  real(dp),                      intent(in)    :: thermal_energy
  type(DegenerateSubspace),      intent(in)    :: subspaces(:)
  class(SubspaceBasis),          intent(in)    :: subspace_bases(:)
  class(BasisStates),            intent(inout) :: subspace_states(:)
  type(AnharmonicData),          intent(in)    :: anharmonic_data
  type(DynamicalMatrix), allocatable           :: output(:)
  
  integer, allocatable :: subspaces_in_coupling(:)
  
  integer :: i
  
  subspaces_in_coupling = filter(subspaces%id .in. this%coupling%ids)
  
  output = [( DynamicalMatrix(anharmonic_data%structure%no_atoms), &
            & i=1,                                                 &
            & size(qpoints)                                        )]
  do i=1,size(this%basis_functions_)
    output = output                                                 &
         & + this%basis_functions_(i)%calculate_dynamical_matrices( &
         &                                   qpoints,               &
         &                                   thermal_energy,        &
         &                                   subspaces,             &
         &                                   subspace_bases,        &
         &                                   subspace_states,       &
         &                                   subspaces_in_coupling, &
         &                                   anharmonic_data        )
  enddo
end function

! Calculate the correction due to double counting
!    for the interpolated potential.
function energy_correction_CouplingBasisFunctions(this,subspaces, &
   & subspace_bases,subspace_states,anharmonic_data) result(output) 
  implicit none
  
  class(CouplingBasisFunctions), intent(in)    :: this
  type(DegenerateSubspace),      intent(in)    :: subspaces(:)
  class(SubspaceBasis),          intent(in)    :: subspace_bases(:)
  class(BasisStates),            intent(inout) :: subspace_states(:)
  type(AnharmonicData),          intent(in)    :: anharmonic_data
  real(dp)                                     :: output
  
  integer :: i
  
  output = 0
  do i=1,size(this%basis_functions_)
    output = output + this%basis_functions_(i)%energy_correction( &
                                               & subspaces,       &
                                               & subspace_bases,  &
                                               & subspace_states, &
                                               & anharmonic_data  )
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_CouplingBasisFunctions(this,input)
  implicit none
  
  class(CouplingBasisFunctions), intent(out) :: this
  type(String),                  intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(SubspaceCoupling)           :: coupling
  type(BasisFunction), allocatable :: basis_functions(:)
  
  select type(this); type is(CouplingBasisFunctions)
    line = split_line(input(1))
    
    if (lower_case(line(1))=='basis') then
      call print_line(ERROR//': Basis functions file appears to be in old &
         &format. Please run update_basis_functions.')
      call quit()
    endif
    
    coupling = SubspaceCoupling(int(line(3:)))
    
    basis_functions = BasisFunction(split_into_sections(input(3:)))
    
    this = CouplingBasisFunctions( coupling        = coupling,       &
                                 & basis_functions = basis_functions )
  class default
    call err()
  end select
end subroutine

function write_CouplingBasisFunctions(this) result(output)
  implicit none
  
  class(CouplingBasisFunctions), intent(in) :: this
  type(String), allocatable                 :: output(:)
  
  select type(this); type is(CouplingBasisFunctions)
    output = [ 'Subspace Coupling: '//this%coupling,          &
             & str(''),                                       &
             & str(this%basis_functions_, separating_line='') ]
  class default
    call err()
  end select
end function

function new_CouplingBasisFunctions_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)     :: input(:)
  type(CouplingBasisFunctions) :: this
  
  call this%read(input)
end function

impure elemental function new_CouplingBasisFunctions_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(CouplingBasisFunctions)  :: this
  
  this = CouplingBasisFunctions(str(input))
end function
end module
