! ======================================================================
! A set of basis functions spanning a given subspace coupling.
! ======================================================================
module coupling_basis_functions_module
  use common_module
  
  use anharmonic_common_module
  
  use basis_function_module
  use sampling_points_module
  implicit none
  
  private
  
  public :: CouplingBasisFunctions
  public :: BasisFunctionsAndSamplingPoints
  public :: size
  public :: generate_basis_functions
  
  type, extends(Stringsable) :: CouplingBasisFunctions
    type(SubspaceCoupling)                    :: coupling
    type(BasisFunction), allocatable, private :: basis_functions_(:)
  contains
    procedure, public :: basis_functions => &
                       & basis_functions_CouplingBasisFunctions
    
    procedure, public :: finalise => finalise_CouplingBasisFunctions
    
    generic,   public  :: energy =>                                           &
                        & energy_RealModeDisplacement_CouplingBasisFunctions, &
                        & energy_ComplexModeDisplacement_CouplingBasisFunctions
    procedure, private :: energy_RealModeDisplacement_CouplingBasisFunctions
    procedure, private :: energy_ComplexModeDisplacement_CouplingBasisFunctions
    generic,   public  :: force =>                                           &
                        & force_RealModeDisplacement_CouplingBasisFunctions, &
                        & force_ComplexModeDisplacement_CouplingBasisFunctions
    procedure, private :: force_RealModeDisplacement_CouplingBasisFunctions
    procedure, private :: force_ComplexModeDisplacement_CouplingBasisFunctions
    
    generic,   public :: braket =>             &
                       & braket_SubspaceState, &
                       & braket_BasisState,    &
                       & braket_BasisStates
    procedure, public :: braket_SubspaceState => &
                       & braket_SubspaceState_CouplingBasisFunctions
    procedure, public :: braket_BasisState => &
                       & braket_BasisState_CouplingBasisFunctions
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_CouplingBasisFunctions
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_CouplingBasisFunctions
    
    procedure, public :: undisplaced_energy => &
                       & undisplaced_energy_CouplingBasisFunctions
    
    procedure, public :: coefficients => &
                       & coefficients_CouplingBasisFunctions
    procedure, public :: set_coefficients => &
                       & set_coefficients_CouplingBasisFunctions
    
    procedure, public :: append => append_CouplingBasisFunctions
    
    ! I/O.
    procedure, public :: read  => read_CouplingBasisFunctions
    procedure, public :: write => write_CouplingBasisFunctions
  end type
  
  interface CouplingBasisFunctions
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

function new_CouplingBasisFunctions(coupling,basis_functions) result(this)
  implicit none
  
  type(SubspaceCoupling), intent(in) :: coupling
  type(BasisFunction),    intent(in) :: basis_functions(:)
  type(CouplingBasisFunctions)       :: this
  
  this%coupling         = coupling
  this%basis_functions_ = basis_functions
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

impure elemental subroutine finalise_CouplingBasisFunctions(this,subspace, &
   & subspace_basis,anharmonic_data)
  implicit none
  
  class(CouplingBasisFunctions), intent(inout) :: this
  type(DegenerateSubspace),      intent(in)    :: subspace
  class(SubspaceBasis),          intent(in)    :: subspace_basis
  type(AnharmonicData),          intent(in)    :: anharmonic_data
  
  if (size(this%coupling%ids)/=1) then
    call print_line(CODE_ERROR//': Calling finalise_subspace_potential &
       &on a potential with coupled subspaces.')
    call err()
  elseif (this%coupling%ids(1)/=subspace%id) then
    call print_line(CODE_ERROR//': Calling finalise_subspace_potential &
       &with the wrong subspace.')
    call err()
  endif
  
  ! Remove constant terms and split basis functions by power.
  this%basis_functions_ = finalise( this%basis_functions_, &
                                  & subspace,              &
                                  & subspace_basis,        &
                                  & anharmonic_data        )
end subroutine

impure elemental function energy_RealModeDisplacement_CouplingBasisFunctions( &
   & this,displacement) result(output)
  implicit none
  
  class(CouplingBasisFunctions), intent(in) :: this
  class(RealModeDisplacement),   intent(in) :: displacement
  real(dp)                                  :: output
  
  output = sum(this%basis_functions_%energy(displacement))
end function

impure elemental function                                                     &
   & energy_ComplexModeDisplacement_CouplingBasisFunctions(this,displacement) &
   & result(output)
  implicit none
  
  class(CouplingBasisFunctions),  intent(in) :: this
  class(ComplexModeDisplacement), intent(in) :: displacement
  complex(dp)                                :: output
  
  output = sum(this%basis_functions_%energy(displacement))
end function

impure elemental function force_RealModeDisplacement_CouplingBasisFunctions( &
   & this,displacement) result(output)
  implicit none
  
  class(CouplingBasisFunctions), intent(in) :: this
  class(RealModeDisplacement),   intent(in) :: displacement
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
  
  class(CouplingBasisFunctions),  intent(in) :: this
  class(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexModeForce)                     :: output
  
  if (size(this%basis_functions_)>0) then
    output = sum(this%basis_functions_%force(displacement))
  else
    output = ComplexModeForce(ComplexSingleForce( displacement%vectors%id, &
                                                & (0.0_dp,0.0_dp)          ))
  endif
end function

impure elemental subroutine braket_SubspaceState_CouplingBasisFunctions(this, &
   & bra,ket,whole_subspace,anharmonic_data)
  implicit none
  
  class(CouplingBasisFunctions), intent(inout)        :: this
  class(SubspaceState),          intent(in)           :: bra
  class(SubspaceState),          intent(in), optional :: ket
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
      call this%basis_functions_(i)%braket(bra, ket, anharmonic_data)
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
      call this%basis_functions_(i)%braket( bra,            &
                                          & ket,            &
                                          & subspace,       &
                                          & subspace_basis, &
                                          & anharmonic_data )
    enddo
    
    call this%basis_functions_%simplify()
  endif
end subroutine

impure elemental subroutine braket_BasisStates_CouplingBasisFunctions(this, &
   & states,thermal_energy,subspace,subspace_basis,whole_subspace,          &
   & anharmonic_data)
  implicit none
  
  class(CouplingBasisFunctions), intent(inout)        :: this
  class(BasisStates),            intent(in)           :: states
  real(dp),                      intent(in)           :: thermal_energy
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
      call this%basis_functions_(j)%braket( states,         &
                                          & thermal_energy, &
                                          & subspace,       &
                                          & subspace_basis, &
                                          & anharmonic_data )
    enddo
    
    call this%basis_functions_%simplify()
  endif
end subroutine

impure elemental function harmonic_expectation_CouplingBasisFunctions(this, &
   & frequency,thermal_energy,anharmonic_data) result(output)
  implicit none
  
  class(CouplingBasisFunctions), intent(in) :: this
  real(dp),                      intent(in) :: frequency
  real(dp),                      intent(in) :: thermal_energy
  type(AnharmonicData),          intent(in) :: anharmonic_data
  real(dp)                                  :: output
  
  output = sum(this%basis_functions_%harmonic_expectation( frequency,      &
                                                         & thermal_energy, &
                                                         & anharmonic_data ))
end function

function generate_basis_functions_SubspaceCoupling(coupling,               &
   & potential_expansion_order,structure,complex_modes,real_modes,qpoints, &
   & subspaces,degenerate_symmetries,vscf_basis_functions_only,            &
   & maximum_weighted_displacement,frequency_of_max_displacement,logfile)  &
   & result(output)
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
  type(OFile),              intent(inout) :: logfile
  type(BasisFunctionsAndSamplingPoints)   :: output
  
  type(SubspaceMonomial), allocatable :: subspace_monomials(:)
  
  type(BasisFunctionsAndUniqueTerms) :: basis_functions
  
  type(BasisFunction), allocatable :: coupling_basis_functions(:)
  type(RealMonomial),  allocatable :: coupling_unique_terms(:)
  
  integer :: i,j,k,ialloc
  
  ! Generate the set of subspace monomials corresponding to the subspace
  !    coupling.
  ! e.g. the coupling [1,2] might have monomials [1,2], [1,1,2] and [1,2,2].
  subspace_monomials = generate_subspace_monomials(        &
     & coupling,                                           &
     & subspaces,                                          &
     & minimum_expansion_order = 2,                        &
     & maximum_expansion_order = potential_expansion_order )
    
  ! Loop over the subspace monomials corresponding to the coupling.
  allocate( coupling_basis_functions(0), &
          & coupling_unique_terms(0),    &
          & stat=ialloc); call err(ialloc)
  do i=1,size(subspace_monomials)
    ! Generate all basis functions for the subspace monomial.
    basis_functions = generate_basis_functions( subspace_monomials(i),     &
                                              & structure,                 &
                                              & complex_modes,             &
                                              & real_modes,                &
                                              & qpoints,                   &
                                              & subspaces,                 &
                                              & degenerate_symmetries,     &
                                              & vscf_basis_functions_only, &
                                              & logfile                    )

    ! Concatenate the terms from each monomial together.
    coupling_basis_functions = [ coupling_basis_functions,       &
                               & basis_functions%basis_functions ]
    coupling_unique_terms = [ coupling_unique_terms,       &
                            & basis_functions%unique_terms ]
  enddo
  
  output = BasisFunctionsAndSamplingPoints(                                 &
     & basis_functions = CouplingBasisFunctions( coupling,                  &
     &                                           coupling_basis_functions), &
     & sampling_points = generate_sampling_points(                          &
     &                                     coupling_unique_terms,           &
     &                                     potential_expansion_order,       &
     &                                     maximum_weighted_displacement,   &
     &                                     frequency_of_max_displacement,   &
     &                                     real_modes                     ) )
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

! Get and set basis function coefficients.
function coefficients_CouplingBasisFunctions(this) result(output)
  implicit none
  
  class(CouplingBasisFunctions), intent(in) :: this
  real(dp), allocatable                     :: output(:)
  
  output = this%basis_functions_%coefficient()
end function

subroutine set_coefficients_CouplingBasisFunctions(this,coefficients)
  implicit none
  
  class(CouplingBasisFunctions), intent(inout) :: this
  real(dp),                      intent(in)    :: coefficients(:)
  
  if (size(this)/=size(coefficients)) then
    call print_line(CODE_ERROR//': Incorrect number of coefficients.')
    call err()
  endif
  
  call this%basis_functions_%set_coefficient(coefficients)
end subroutine

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
