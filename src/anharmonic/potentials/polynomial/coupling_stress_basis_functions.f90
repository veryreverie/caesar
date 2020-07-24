! ======================================================================
! A set of stress basis functions spanning a given subspace coupling.
! ======================================================================
module coupling_stress_basis_functions_module
  use common_module
  
  use anharmonic_common_module
  
  use polynomial_interpolator_module
  use basis_function_module
  use stress_basis_function_module
  use sampling_points_module
  use sample_result_module
  implicit none
  
  private
  
  public :: CouplingStressBasisFunctions
  public :: size
  public :: generate_stress_basis_functions
  
  type, extends(StressBase) :: CouplingStressBasisFunctions
    type(SubspaceCoupling)                          :: coupling
    type(StressBasisFunction), allocatable, private :: basis_functions_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_CouplingStressBasisFunctions
    
    procedure, public :: basis_functions => &
                       & basis_functions_CouplingStressBasisFunctions
    
    procedure, public :: undisplaced_stress => &
                       & undisplaced_stress_CouplingStressBasisFunctions
    
    procedure, public :: zero_stress => &
                       & zero_stress_CouplingStressBasisFunctions
    procedure, public :: add_constant => &
                       & add_constant_CouplingStressBasisFunctions
    
    procedure, public :: stress_RealModeDisplacement => &
       & stress_RealModeDisplacement_CouplingStressBasisFunctions
    procedure, public :: stress_ComplexModeDisplacement => &
       & stress_ComplexModeDisplacement_CouplingStressBasisFunctions
    
    procedure, public :: braket_SubspaceBraKet => &
                       & braket_SubspaceBraKet_CouplingStressBasisFunctions
    procedure, public :: braket_BasisState => &
                       & braket_BasisState_CouplingStressBasisFunctions
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_CouplingStressBasisFunctions
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_CouplingStressBasisFunctions
    
    procedure, public :: append => append_CouplingStressBasisFunctions
    
    procedure, public :: fit_coefficients => &
                       & fit_coefficients_CouplingStressBasisFunctions
    procedure, public :: interpolate_coefficients => &
                       & interpolate_coefficients_CouplingStressBasisFunctions
    procedure, public :: calculate_dynamical_matrices => &
                    & calculate_dynamical_matrices_CouplingStressBasisFunctions
    procedure, public :: stress_correction => &
                       & stress_correction_CouplingStressBasisFunctions
    
    ! I/O.
    procedure, public :: read  => read_CouplingStressBasisFunctions
    procedure, public :: write => write_CouplingStressBasisFunctions
  end type
  
  interface CouplingStressBasisFunctions
    module procedure new_CouplingStressBasisFunctions_empty
    module procedure new_CouplingStressBasisFunctions
    module procedure new_CouplingStressBasisFunctions_Strings
    module procedure new_CouplingStressBasisFunctions_StringArray
  end interface
  
  interface size
    module procedure size_CouplingStressBasisFunctions
  end interface
  
  interface generate_stress_basis_functions
    module procedure generate_stress_basis_functions_SubspaceCoupling
  end interface
contains

! Constructor and size function.
impure elemental function new_CouplingStressBasisFunctions_empty(coupling) &
   & result(this) 
  implicit none
  
  type(SubspaceCoupling),    intent(in) :: coupling
  type(CouplingStressBasisFunctions)    :: this
  
  integer :: ialloc
  
  this%coupling = coupling
  allocate(this%basis_functions_(0), stat=ialloc); call err(ialloc)
end function

function new_CouplingStressBasisFunctions(coupling,basis_functions) &
   & result(this)
  implicit none
  
  type(SubspaceCoupling),    intent(in) :: coupling
  type(StressBasisFunction), intent(in) :: basis_functions(:)
  type(CouplingStressBasisFunctions)    :: this
  
  this%coupling = coupling
  this%basis_functions_ = basis_functions
end function

impure elemental function representation_CouplingStressBasisFunctions() &
   & result(output)
  implicit none
  
  type(String) :: output
  
  output = 'Polynomial stress basis functions'
end function

function size_CouplingStressBasisFunctions(this) result(output)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(in) :: this
  integer                                         :: output
  
  output = size(this%basis_functions_)
end function

function basis_functions_CouplingStressBasisFunctions(this) result(output)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(in) :: this
  type(StressBasisFunction), allocatable          :: output(:)
  
  output = this%basis_functions_
end function

impure elemental function                                           &
   & stress_RealModeDisplacement_CouplingStressBasisFunctions(this, &
   & displacement) result(output)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(in) :: this
  type(RealModeDisplacement),          intent(in) :: displacement
  type(RealMatrix)                                :: output
  
  if (size(this%basis_functions_)>0) then
    output = sum(this%basis_functions_%stress(displacement))
  else
    output = dblemat(zeroes(3,3))
  endif
end function

impure elemental function                                              &
   & stress_ComplexModeDisplacement_CouplingStressBasisFunctions(this, &
   & displacement) result(output)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(in) :: this
  type(ComplexModeDisplacement),       intent(in) :: displacement
  type(ComplexMatrix)                             :: output
  
  if (size(this%basis_functions_)>0) then
    output = sum(this%basis_functions_%stress(displacement))
  else
    output = cmplxmat(zeroes(3,3))
  endif
end function

impure elemental subroutine                                          &
   & braket_SubspaceBraKet_CouplingStressBasisFunctions(this,braket, &
   & whole_subspace,anharmonic_data)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(inout)        :: this
  class(SubspaceBraKet),               intent(in)           :: braket
  logical,                             intent(in), optional :: whole_subspace
  type(AnharmonicData),                intent(in)           :: anharmonic_data
  
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

impure elemental subroutine braket_BasisState_CouplingStressBasisFunctions( &
   & this,bra,ket,subspace,subspace_basis,whole_subspace,anharmonic_data)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(inout)        :: this
  class(BasisState),                   intent(in)           :: bra
  class(BasisState),                   intent(in), optional :: ket
  type(DegenerateSubspace),            intent(in)           :: subspace
  class(SubspaceBasis),                intent(in)           :: subspace_basis
  logical,                             intent(in), optional :: whole_subspace
  type(AnharmonicData),                intent(in)           :: anharmonic_data
  
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

impure elemental subroutine braket_BasisStates_CouplingStressBasisFunctions( &
   & this,states,subspace,subspace_basis,whole_subspace,anharmonic_data) 
  implicit none
  
  class(CouplingStressBasisFunctions), intent(inout)        :: this
  class(BasisStates),                  intent(inout)        :: states
  type(DegenerateSubspace),            intent(in)           :: subspace
  class(SubspaceBasis),                intent(in)           :: subspace_basis
  logical,                             intent(in), optional :: whole_subspace
  type(AnharmonicData),                intent(in)           :: anharmonic_data
  
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

impure elemental function harmonic_expectation_CouplingStressBasisFunctions( &
   & this,frequency,thermal_energy,supercell_size,anharmonic_data)           &
   & result(output) 
  implicit none
  
  class(CouplingStressBasisFunctions), intent(in) :: this
  real(dp),                            intent(in) :: frequency
  real(dp),                            intent(in) :: thermal_energy
  integer,                             intent(in) :: supercell_size
  type(AnharmonicData),                intent(in) :: anharmonic_data
  type(RealMatrix)                                :: output
  
  integer :: i
  
  if (size(this%basis_functions_)==0) then
    output = dblemat(zeroes(3,3))
  else
    output = sum([( this%basis_functions_(i)%harmonic_expectation(    &
                  &                                frequency,         &
                  &                                thermal_energy,    &
                  &                                supercell_size,    &
                  &                                anharmonic_data ), &
                  & i=1,                                              &
                  & size(this%basis_functions_)                       )])
  endif
end function

! Return the stress at zero displacement.
impure elemental function undisplaced_stress_CouplingStressBasisFunctions( &
   & this) result(output)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(in) :: this
  type(RealMatrix)                                :: output
  
  type(RealSingleDisplacement) :: zero_displacement(0)
  
  output = this%stress(RealModeDisplacement(zero_displacement))
end function

impure elemental subroutine zero_stress_CouplingStressBasisFunctions(this)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(inout) :: this
  
  call print_line(CODE_ERROR//': zero_energy() cannot be called for type &
     &CouplingStressBasisFunctions.')
  call err()
end subroutine

impure elemental subroutine add_constant_CouplingStressBasisFunctions(this, &
   & input) 
  implicit none
  
  class(CouplingStressBasisFunctions), intent(inout) :: this
  type(RealMatrix),                    intent(in)    :: input
  
  call print_line(CODE_ERROR//': add_constant() cannot be called for type &
     &CouplingStressBasisFunctions.')
  call err()
end subroutine

! Append another StressCouplingBasisFunctions to this.
subroutine append_CouplingStressBasisFunctions(this,that)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(inout) :: this
  class(CouplingStressBasisFunctions), intent(in)    :: that
  
  if (this%coupling/=that%coupling) then
    call print_line(CODE_ERROR//': Appending incompatible basis functions.')
    call err()
  endif
  
  this%basis_functions_ = [this%basis_functions_, that%basis_functions_]
end subroutine

! Generate stress basis functions.
function generate_stress_basis_functions_SubspaceCoupling(coupling,    &
   & stress_expansion_order,structure,complex_modes,qpoints,subspaces, &
   & degenerate_symmetries,vscf_basis_functions_only,logfile) result(output) 
  implicit none
  
  type(SubspaceCoupling),   intent(in)    :: coupling
  integer,                  intent(in)    :: stress_expansion_order
  type(StructureData),      intent(in)    :: structure
  type(ComplexMode),        intent(in)    :: complex_modes(:)
  type(QpointData),         intent(in)    :: qpoints(:)
  type(DegenerateSubspace), intent(in)    :: subspaces(:)
  type(DegenerateSymmetry), intent(in)    :: degenerate_symmetries(:)
  logical,                  intent(in)    :: vscf_basis_functions_only
  type(OFile),              intent(inout) :: logfile
  type(CouplingStressBasisFunctions)      :: output
  
  type(SubspaceMonomial), allocatable :: subspace_monomials(:)
  
  type(StressBasisFunction), allocatable :: coupling_basis_functions(:)
  type(StressBasisFunction), allocatable :: monomial_basis_functions(:)
  
  integer :: i,ialloc
  
  ! Generate the set of subspace monomials corresponding to the subspace
  !    coupling.
  ! e.g. the coupling [1,2] might have monomials [1,2], [1,1,2] and [1,2,2].
  subspace_monomials = generate_subspace_monomials(     &
     & coupling,                                        &
     & subspaces,                                       &
     & minimum_expansion_order = 1,                     &
     & maximum_expansion_order = stress_expansion_order )
  
  ! Loop over the subspace monomials corresponding to the coupling.
  allocate(coupling_basis_functions(0), stat=ialloc); call err(ialloc)
  do i=1,size(subspace_monomials)
    monomial_basis_functions = generate_stress_basis_functions( &
                                   & subspace_monomials(i),     &
                                   & structure,                 &
                                   & complex_modes,             &
                                   & qpoints,                   &
                                   & subspaces,                 &
                                   & degenerate_symmetries,     &
                                   & vscf_basis_functions_only, &
                                   & logfile                    )
    coupling_basis_functions = [ coupling_basis_functions, &
                               & monomial_basis_functions  ]
  enddo
  
  output = CouplingStressBasisFunctions( coupling,                &
                                       & coupling_basis_functions )
end function

! Uses L2 regression to calculate the coefficients of a set of stress basis
!    functions.
subroutine fit_coefficients_CouplingStressBasisFunctions(this, &
   & sampling_points,sample_results,stress) 
  implicit none
  
  class(CouplingStressBasisFunctions), intent(inout)        :: this
  type(RealModeDisplacement),          intent(in)           :: sampling_points(:)
  type(SampleResult),                  intent(in)           :: sample_results(:)
  class(StressData),                   intent(in), optional :: stress
  
  type(RealMatrix) :: stress_tensor
  
  real(dp), allocatable :: a(:,:)
  real(dp), allocatable :: b(:)
  
  real(dp), allocatable :: coefficients(:)
  
  integer :: i,j,ialloc
  
  ! Check inputs are consistent.
  if (size(sampling_points)/=size(sample_results)) then
    call print_line(CODE_ERROR//': The number of sampling points does not &
       &match the number of results.')
    call err()
  endif
  
  ! Calculate the stress due to each basis function at each sampling point.
  allocate( a(9*size(sampling_points), size(this%basis_functions_)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this%basis_functions_)
    do j=1,size(sampling_points)
      stress_tensor = this%basis_functions_(i)%stress(sampling_points(j))
      a((j-1)*9+1:j*9, i) = reshape(dble(stress_tensor), [9])
    enddo
  enddo
  
  ! Calculate the stress sampled at each sampling point.
  allocate(b(9*size(sampling_points)), stat=ialloc); call err(ialloc)
  do i=1,size(sampling_points)
    stress_tensor = sample_results(i)%stress()
    
    ! Subtract the existing stress.
    if (present(stress)) then
      stress_tensor = stress_tensor - stress%stress(sampling_points(i))
    endif
    
    b((i-1)*9+1:i*9) = reshape(dble(stress_tensor), [9])
  enddo
  
  ! Run linear least squares to get the basis function coefficients.
  ! This finds x s.t. (a.x-b)^2 is minimised.
  coefficients = dble(linear_least_squares(a, b))
  
  this%basis_functions_ = coefficients * this%basis_functions_
end subroutine

! Calculate the contribution to a given monomial from the interpolation of
!    this basis function.
! The result is given as a cartesian tensor.
impure elemental function                                                 &
   & interpolate_coefficients_CouplingStressBasisFunctions(this,monomial, &
   & interpolator) result(output)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(in) :: this
  type(ComplexMonomial),               intent(in) :: monomial
  type(PolynomialInterpolator),        intent(in) :: interpolator
  type(ComplexMatrix)                             :: output
  
  output = sum(this%basis_functions_%interpolate_coefficients( monomial,    &
                                                             & interpolator ))
end function

! Calculate this basis function's contribution to the effective dynamical
!    matrix from which the potential can be interpolated in the large-supercell
!    limit.
function calculate_dynamical_matrices_CouplingStressBasisFunctions(this, &
   & qpoints,thermal_energy,subspaces,subspace_bases,subspace_states,    &
   & anharmonic_data) result(output) 
  implicit none
  
  class(CouplingStressBasisFunctions), intent(in)    :: this
  type(QpointData),                    intent(in)    :: qpoints(:)
  real(dp),                            intent(in)    :: thermal_energy
  type(DegenerateSubspace),            intent(in)    :: subspaces(:)
  class(SubspaceBasis),                intent(in)    :: subspace_bases(:)
  class(BasisStates),                  intent(inout) :: subspace_states(:)
  type(AnharmonicData),                intent(in)    :: anharmonic_data
  type(StressDynamicalMatrix), allocatable           :: output(:)
  
  integer, allocatable :: subspaces_in_coupling(:)
  
  integer :: i
  
  subspaces_in_coupling = filter(subspaces%id .in. this%coupling%ids)
  
  output = [( StressDynamicalMatrix(anharmonic_data%structure%no_atoms), &
            & i=1,                                                       &
            & size(qpoints)                                              )]
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
!    for the interpolated stress.
function stress_correction_CouplingStressBasisFunctions(this,subspaces, &
   & subspace_bases,subspace_states,anharmonic_data) result(output) 
  implicit none
  
  class(CouplingStressBasisFunctions), intent(in)    :: this
  type(DegenerateSubspace),            intent(in)    :: subspaces(:)
  class(SubspaceBasis),                intent(in)    :: subspace_bases(:)
  class(BasisStates),                  intent(inout) :: subspace_states(:)
  type(AnharmonicData),                intent(in)    :: anharmonic_data
  type(RealMatrix)                                   :: output
  
  integer :: i
  
  output = dblemat(zeroes(3,3))
  do i=1,size(this%basis_functions_)
    output = output + this%basis_functions_(i)%stress_correction( &
                                               & subspaces,       &
                                               & subspace_bases,  &
                                               & subspace_states, &
                                               & anharmonic_data  )
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_CouplingStressBasisFunctions(this,input)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(out) :: this
  type(String),                        intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(SubspaceCoupling)                 :: coupling
  type(StressBasisFunction), allocatable :: basis_functions(:)
  
  select type(this); type is(CouplingStressBasisFunctions)
    line = split_line(input(1))
    coupling = SubspaceCoupling(int(line(3:)))
    
    basis_functions = StressBasisFunction(split_into_sections( &
                              & input(3:),                     &
                              & separating_line=repeat('-',50) ))
    
    this = CouplingStressBasisFunctions( coupling        = coupling,       &
                                       & basis_functions = basis_functions )
  class default
    call err()
  end select
end subroutine

function write_CouplingStressBasisFunctions(this) result(output)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(in) :: this
  type(String), allocatable                       :: output(:)
  
  select type(this); type is(CouplingStressBasisFunctions)
    output = [ 'Subspace Coupling: '//this%coupling,                      &
             & str(repeat('-',50)),                                       &
             & str(this%basis_functions_, separating_line=repeat('-',50)) ]
  class default
    call err()
  end select
end function

function new_CouplingStressBasisFunctions_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)           :: input(:)
  type(CouplingStressBasisFunctions) :: this
  
  call this%read(input)
end function

impure elemental function new_CouplingStressBasisFunctions_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in)       :: input
  type(CouplingStressBasisFunctions)  :: this
  
  this = CouplingStressBasisFunctions(str(input))
end function
end module
