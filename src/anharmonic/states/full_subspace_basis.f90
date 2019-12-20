! ======================================================================
! A basis of states which spans a full subspace.
! ======================================================================
module full_subspace_basis_module
  use common_module
  
  use anharmonic_common_module
  
  use wavevector_state_module
  use wavevector_states_module
  use wavevector_basis_module
  use full_subspace_wavefunctions_module
  use core_shell_thermodynamics_module
  use calculate_weights_module
  implicit none
  
  private
  
  public :: startup_full_subspace_basis
  
  public :: FullSubspaceBasis
  
  public :: initial_ground_state
  
  ! All states spanning the subspace.
  type, extends(SubspaceBasis) :: FullSubspaceBasis
    ! The maximum power of the monomial states.
    ! This is also the maximum occupation of the harmonic basis states.
    integer  :: maximum_power
    ! The expansion order of the potential.
    ! This is also the limit on coupling between states.
    integer :: expansion_order
    ! The ID of the subspace.
    integer  :: subspace_id
    ! The states, wavevector by wavevector.
    ! N.B. this only includes one wavevector from each symmetry-related set.
    type(WavevectorBasis), allocatable :: wavevectors(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_FullSubspaceBasis
    
    ! Print the harmonic ground-state wavefunction of the basis.
    procedure, public :: ground_state_wavefunction
    
    ! Generate the first guess at states.
    procedure, public :: initial_states => initial_states_FullSubspaceBasis
    
    ! Generate the eigenstates of a single-subspace potential.
    procedure, public :: calculate_states => calculate_states_FullSubspaceBasis
    
    ! Return the modes spanned by this basis.
    procedure, public :: mode_ids => mode_ids_FullSubspaceBasis
    procedure, public :: paired_mode_ids => paired_mode_ids_FullSubspaceBasis
    
    ! Procedures involving individual states.
    procedure, public :: inner_product => &
                       & inner_product_FullSubspaceBasis
    procedure, public :: integrate_BasisState => &
                       & integrate_BasisState_FullSubspaceBasis
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_FullSubspaceBasis
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_FullSubspaceBasis
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_FullSubspaceBasis
    procedure, public :: wavefunction => wavefunction_FullSubspaceBasis
    
    ! Procedures involving sets of states.
    procedure, public :: thermodynamic_data => &
                       & thermodynamic_data_FullSubspaceBasis
    procedure, public :: wavefunctions => wavefunctions_FullSubspaceBasis
    procedure, public :: integrate_BasisStates => &
                       & integrate_BasisStates_FullSubspaceBasis
    
    ! I/O.
    procedure, public :: read  => read_FullSubspaceBasis
    procedure, public :: write => write_FullSubspaceBasis
  end type
  
  interface FullSubspaceBasis
    module procedure new_FullSubspaceBasis
    module procedure new_FullSubspaceBasis_SubspaceBasis
    module procedure new_FullSubspaceBasis_subspace
    module procedure new_FullSubspaceBasis_Strings
    module procedure new_FullSubspaceBasis_StringArray
  end interface
contains

! Startup procedure.
subroutine startup_full_subspace_basis()
  implicit none
  
  type(FullSubspaceBasis) :: basis
  
  call basis%startup()
end subroutine

! ----------------------------------------------------------------------
! FullSubspaceBasis methods.
! ----------------------------------------------------------------------
! Constructors.
function new_FullSubspaceBasis(maximum_power,expansion_order,subspace_id, &
   & frequency,wavevectors) result(this)
  implicit none
  
  integer,               intent(in) :: maximum_power
  integer,               intent(in) :: expansion_order
  integer,               intent(in) :: subspace_id
  real(dp),              intent(in) :: frequency
  type(WavevectorBasis), intent(in) :: wavevectors(:)
  type(FullSubspaceBasis)           :: this
  
  this%maximum_power   = maximum_power
  this%expansion_order = expansion_order
  this%subspace_id     = subspace_id
  this%frequency       = frequency
  this%wavevectors     = wavevectors
end function

recursive function new_FullSubspaceBasis_SubspaceBasis(input) result(this)
  implicit none
  
  class(SubspaceBasis), intent(in) :: input
  type(FullSubspaceBasis)          :: this
  
  select type(input); type is (FullSubspaceBasis)
    this = input
  type is(SubspaceBasisPointer)
    ! WORKAROUND: ifort doesn't recognise the interface to this function
    !    from within this function, so the full name is used instead.
    !this = FullSubspaceBasis(input%basis())
    this = new_FullSubspaceBasis_SubspaceBasis(input%basis())
  class default
    call err()
  end select
end function

! Type representation.
impure elemental function representation_FullSubspaceBasis() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'full_subspace'
end function

! Generates states up to a given power, spanning the whole subspace.
function new_FullSubspaceBasis_subspace(subspace,frequency,modes,qpoints, &
   & supercell,maximum_power,potential_expansion_order,symmetries)        &
   & result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency
  type(ComplexMode),        intent(in) :: modes(:)
  type(QpointData),         intent(in) :: qpoints(:)
  type(StructureData),      intent(in) :: supercell
  integer,                  intent(in) :: maximum_power
  integer,                  intent(in) :: potential_expansion_order
  type(SymmetryOperator),   intent(in) :: symmetries(:)
  type(FullSubspaceBasis)              :: output
  
  type(WavevectorBasis), allocatable :: wavevectors(:)
  
  wavevectors = WavevectorBasis( subspace,                  &
                               & frequency,                 &
                               & modes,                     &
                               & qpoints,                   &
                               & maximum_power,             &
                               & potential_expansion_order, &
                               & symmetries                 )
  
  output = FullSubspaceBasis( maximum_power,             &
                            & potential_expansion_order, &
                            & subspace%id,               &
                            & frequency,                 &
                            & wavevectors                )
end function

! Returns the harmonic ground-state wavefunction for the basis.
function ground_state_wavefunction(this,subspace,supercell) result(output)
  implicit none
  
  class(FullSubspaceBasis), intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  type(String)                         :: output
  
  real(dp) :: mass
  
  real(dp)                  :: coefficient
  type(String), allocatable :: terms(:)
  
  integer :: i,ialloc
  
  ! Calculate the (geometric) average mass.
  mass = product(supercell%atoms%mass())
  mass = mass**(1.0_dp/size(supercell%atoms))
  
  ! Calculate the coefficient.
  coefficient = 1
  allocate(terms(0), stat=ialloc); call err(ialloc)
  do i=1,size(subspace)
    if (subspace%mode_ids(i)==subspace%paired_ids(i)) then
      ! |0_i> = sqrt(sqrt(m*w/pi)) exp(- 1/2 N w (u_i)^2 )
      coefficient = coefficient * (mass*this%frequency/PI)**0.25_dp
      terms = [terms, 'u'//subspace%mode_ids(i)//'^2']
    elseif (subspace%mode_ids(i)<subspace%paired_ids(i)) then
      ! |0_i,0_j> = sqrt(2*m*w/pi) exp(- N w |u_i|^2 )
      coefficient = coefficient * (2*mass*this%frequency/PI)**0.5_dp
      terms = [ terms,                                                    &
              & '2*u'//subspace%mode_ids(i)//'*u'//subspace%paired_ids(i) ]
    endif
  enddo
  
  output = coefficient//'*e^('               // &
     & (-supercell%sc_size*this%frequency/2) // &
     & '*('//join(terms,'+')//'))'
end function

! Generate an initial guess at states.
impure elemental function initial_states_FullSubspaceBasis(this,subspace, &
   & thermal_energy,anharmonic_data) result(output)
  implicit none
  
  class(FullSubspaceBasis), intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: thermal_energy
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(BasisStatesPointer)             :: output
  
  type(WavevectorState) :: ground_state
  
  integer :: ialloc
  
  ! Generate the state |0>.
  ground_state = initial_ground_state(this)
  
  ! Generate the set of states {|0>}.
  output = BasisStatesPointer(WavevectorStates( &
                & subspace_id = subspace%id,    &
                & states      = [ground_state], &
                & energies    = [0.0_dp],       &
                & weights     = [1.0_dp]        ))
end function

! Generate initial guess. This is simply the basis state |0>, i.e. the
!    ground state of the effective harmonic potential from which the basis
!    states were generated.
impure elemental function initial_ground_state(basis) result(output)
  implicit none
  
  type(FullSubspaceBasis), intent(in) :: basis
  type(WavevectorState)               :: output
  
  type(WavevectorBasis) :: wavevector_basis
  
  ! Find the basis at wavevector [0,0,0].
  wavevector_basis = basis%wavevectors(            &
     & first(is_int(basis%wavevectors%wavevector)) )
  
  ! Construct the ground state at [0,0,0].
  output = wavevector_basis%initial_ground_state()
end function

! Calculate the eigenstates of a single-subspace potential.
impure elemental function calculate_states_FullSubspaceBasis(this,subspace, &
   & subspace_potential,thermal_energy,convergence_data,anharmonic_data)    &
   & result(output)
  implicit none
  
  class(FullSubspaceBasis), intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  class(PotentialData),     intent(in) :: subspace_potential
  real(dp),                 intent(in) :: thermal_energy
  type(ConvergenceData),    intent(in) :: convergence_data
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(BasisStatesPointer)             :: output
  
  type(WavevectorState)              :: state
  type(WavevectorState), allocatable :: states(:)
  real(dp),              allocatable :: energies(:)
  real(dp),              allocatable :: weights(:)
  
  type(WavevectorStates) :: wavevector_states
  
  integer :: i,ialloc
  
  allocate( states(0),   &
          & energies(0), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this%wavevectors)
    wavevector_states = WavevectorStates(        &
       & this%wavevectors(i)%calculate_states(   &
       &            subspace,                    &
       &            subspace_potential,          &
       &            thermal_energy,              &
       &            convergence_data,            &
       &            anharmonic_data            ) )
    states = [states, wavevector_states%states]
    energies = [energies, wavevector_states%energies]
  enddo
  
  weights = calculate_weights(energies, thermal_energy)
  
  output = BasisStatesPointer(WavevectorStates( subspace%id, &
                                              & states,      &
                                              & energies,    &
                                              & weights      ))
end function

impure elemental function wavefunction_FullSubspaceBasis(this,state, &
   & supercell) result(output)
  implicit none
  
  class(FullSubspaceBasis), intent(in) :: this
  type(WavevectorState),    intent(in) :: state
  type(StructureData),      intent(in) :: supercell
  type(String)                         :: output
  
  ! TODO
  output = ''
end function

function mode_ids_FullSubspaceBasis(this,subspace,anharmonic_data) &
   & result(output)
  implicit none
  
  class(FullSubspaceBasis), intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  integer, allocatable                 :: output(:)
  
  output = this%wavevectors(1)%mode_ids(subspace,anharmonic_data)
end function

function paired_mode_ids_FullSubspaceBasis(this,subspace,anharmonic_data) &
   & result(output)
  implicit none
  
  class(FullSubspaceBasis), intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  integer, allocatable                 :: output(:)
  
  output = this%wavevectors(1)%paired_mode_ids(subspace,anharmonic_data)
end function

impure elemental function inner_product_FullSubspaceBasis(this,bra,ket, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(FullSubspaceBasis), intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(WavevectorState)              :: full_bra
  type(WavevectorState), allocatable :: full_ket
  
  integer :: i
  
  full_bra = WavevectorState(bra)
  if (present(ket)) then
    full_ket = WavevectorState(ket)
  endif
  
  i = first(this%wavevectors%wavevector == full_bra%wavevector)
  
  output = this%wavevectors(i)%inner_product( full_bra,       &
                                            & full_ket,       &
                                            & subspace,       &
                                            & anharmonic_data )
end function

impure elemental function integrate_BasisState_FullSubspaceBasis(this, &
   & bra,monomial,ket,subspace,anharmonic_data) result(output)
  implicit none
  
  class(FullSubspaceBasis), intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  type(SparseMonomial),     intent(in)           :: monomial
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  complex(dp)                                    :: output
  
  type(WavevectorState)              :: full_bra
  type(WavevectorState), allocatable :: full_ket
  
  integer :: i
  
  full_bra = WavevectorState(bra)
  if (present(ket)) then
    full_ket = WavevectorState(ket)
  endif
  
  i = first(this%wavevectors%wavevector == full_bra%wavevector)
  
  output = this%wavevectors(i)%integrate( full_bra,       &
                                        & monomial,       &
                                        & full_ket,       &
                                        & subspace,       &
                                        & anharmonic_data )
end function

impure elemental function kinetic_energy_FullSubspaceBasis(this,bra,ket, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(FullSubspaceBasis), intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(WavevectorState)              :: full_bra
  type(WavevectorState), allocatable :: full_ket
  
  integer :: i
  
  full_bra = WavevectorState(bra)
  if (present(ket)) then
    full_ket = WavevectorState(ket)
  endif
  
  i = first(this%wavevectors%wavevector == full_bra%wavevector)
  
  output = this%wavevectors(i)%kinetic_energy( full_bra,       &
                                             & full_ket,       &
                                             & subspace,       &
                                             & anharmonic_data )
end function

impure elemental function harmonic_potential_energy_FullSubspaceBasis( &
   & this,bra,ket,subspace,anharmonic_data) result(output)
  implicit none
  
  class(FullSubspaceBasis), intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(WavevectorState)              :: full_bra
  type(WavevectorState), allocatable :: full_ket
  
  integer :: i
  
  full_bra = WavevectorState(bra)
  if (present(ket)) then
    full_ket = WavevectorState(ket)
  endif
  
  i = first(this%wavevectors%wavevector == full_bra%wavevector)
  
  output = this%wavevectors(i)%harmonic_potential_energy( full_bra,       &
                                                        & full_ket,       &
                                                        & subspace,       &
                                                        & anharmonic_data )
end function

impure elemental function kinetic_stress_FullSubspaceBasis(this,bra,ket, &
   & subspace,stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(FullSubspaceBasis), intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  type(StressPrefactors),   intent(in)           :: stress_prefactors
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(RealMatrix)                               :: output
  
  type(WavevectorState)              :: full_bra
  type(WavevectorState), allocatable :: full_ket
  
  integer :: i
  
  full_bra = WavevectorState(bra)
  if (present(ket)) then
    full_ket = WavevectorState(ket)
  endif
  
  i = first(this%wavevectors%wavevector == full_bra%wavevector)
  
  output = this%wavevectors(i)%kinetic_stress( full_bra,          &
                                             & full_ket,          &
                                             & subspace,          &
                                             & stress_prefactors, &
                                             & anharmonic_data    )
end function

! Thermodynamic data. Energy, entropy, free energy etc.
impure elemental function thermodynamic_data_FullSubspaceBasis(this,    &
   & thermal_energy,states,subspace,subspace_potential,subspace_stress, &
   & stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(FullSubspaceBasis), intent(in)           :: this
  real(dp),                 intent(in)           :: thermal_energy
  class(BasisStates),       intent(in)           :: states
  type(DegenerateSubspace), intent(in)           :: subspace
  class(PotentialData),     intent(in)           :: subspace_potential
  class(StressData),        intent(in), optional :: subspace_stress
  type(StressPrefactors),   intent(in), optional :: stress_prefactors
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(ThermodynamicData)                        :: output
  
  ! Calculate the thermodynamic properties for the system.
  output = core_shell_thermodynamics( &
              & thermal_energy,       &
              & this%frequency,       &
              & size(subspace),       &
              & subspace,             &
              & this%wavevectors,     &
              & states,               &
              & subspace_potential,   &
              & subspace_stress,      &
              & stress_prefactors,    &
              & anharmonic_data       )
end function

! Wavefunctions.
impure elemental function wavefunctions_FullSubspaceBasis(this,states, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(FullSubspaceBasis),  intent(in) :: this
  class(BasisStates),        intent(in) :: states
  type(DegenerateSubspace),  intent(in) :: subspace
  type(AnharmonicData),      intent(in) :: anharmonic_data
  type(SubspaceWavefunctionsPointer)    :: output
  
  type(WavevectorStates)          :: full_states
  type(String)                    :: ground_state
  type(String), allocatable       :: state_wavefunctions(:)
  type(FullSubspaceWavefunctions) :: wavefunctions
  
  integer :: ialloc
  
  full_states = WavevectorStates(states)
  
  ! Construct the wavefunction of |0>.
  ground_state = this%ground_state_wavefunction(      &
               & subspace,                            &
               & anharmonic_data%anharmonic_supercell )
  
  ! Construct the wavefunctions of each state, in terms of |0>.
  state_wavefunctions = this%wavefunction(  &
     & full_states%states,                  &
     & anharmonic_data%anharmonic_supercell )
  
  ! Pack the wavefunctions into the output,
  !    of concrete type FullSubspaceWavefunctions.
  wavefunctions = FullSubspaceWavefunctions( &
                     & subspace%id,          &
                     & subspace%mode_ids,    &
                     & subspace%paired_ids,  &
                     & ground_state,         &
                     & full_states%energies, &
                     & state_wavefunctions   )
  
  ! Convert the output to abstract class SubspaceWavefunctions.
  output = SubspaceWavefunctionsPointer(wavefunctions)
end function

! Integrate a monomial.
impure elemental function integrate_BasisStates_FullSubspaceBasis(this, &
   & states,thermal_energy,monomial,subspace,anharmonic_data) result(output)
  implicit none
  
  class(FullSubspaceBasis), intent(in) :: this
  class(BasisStates),       intent(in) :: states
  real(dp),                 intent(in) :: thermal_energy
  type(SparseMonomial),     intent(in) :: monomial
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  complex(dp)                          :: output
  
  type(WavevectorStates) :: full_states
  
  ! Convert the states to type WavevectorStates.
  full_states = WavevectorStates(states)
  
  output = integrate( this%wavevectors, &
                    & full_states,      &
                    & thermal_energy,   &
                    & monomial,         &
                    & anharmonic_data   )
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_FullSubspaceBasis(this,input)
  implicit none
  
  class(FullSubspaceBasis), intent(out) :: this
  type(String),             intent(in)  :: input(:)
  
  integer                            :: maximum_power
  integer                            :: expansion_order
  integer                            :: subspace_id
  real(dp)                           :: frequency
  type(WavevectorBasis), allocatable :: wavevectors(:)
  
  integer, allocatable :: starting_lines(:)
  integer, allocatable :: ending_lines(:)
  
  type(String), allocatable :: line(:)
  type(String), allocatable :: wavevector_lines(:)
  
  integer :: i,ialloc
  
  select type(this); type is(FullSubspaceBasis)
    line = split_line(input(1))
    maximum_power = int(line(4))
    
    line = split_line(input(2))
    expansion_order = int(line(4))
    
    line = split_line(input(3))
    subspace_id = int(line(3))
    
    line = split_line(input(4))
    frequency = dble(line(3))
    
    starting_lines = [integer::]
    do i=5,size(input)
      line = split_line(input(i))
      if (size(line)>0) then
        if (line(1)=='Wavevector') then
          starting_lines = [starting_lines, i]
        endif
      endif
    enddo
    
    ending_lines = [starting_lines(2:)-1, size(input)]
    
    allocate(wavevectors(size(starting_lines)), stat=ialloc); call err(ialloc)
    do i=1,size(starting_lines)
      wavevector_lines = [input(:4), input(starting_lines(i):ending_lines(i))]
      wavevectors(i) = WavevectorBasis(wavevector_lines)
    enddo
    
    this = FullSubspaceBasis( maximum_power,   &
                            & expansion_order, &
                            & subspace_id,     &
                            & frequency,       &
                            & wavevectors      )
  class default
    call err()
  end select
end subroutine

function write_FullSubspaceBasis(this) result(output)
  implicit none
  
  class(FullSubspaceBasis), intent(in) :: this
  type(String), allocatable            :: output(:)
  
  type(String), allocatable :: wavevector_strings(:)
  
  integer :: i
  
  select type(this); type is(FullSubspaceBasis)
    output = [ 'Maximum power   : '//this%maximum_power,   &
             & 'Expansion order : '//this%expansion_order, &
             & 'Subspace        : '//this%subspace_id,     &
             & 'Frequency       : '//this%frequency        ]
    do i=1,size(this%wavevectors)
      wavevector_strings = str(this%wavevectors(i))
      output = [ output, wavevector_strings(5:) ]
    enddo
  class default
    call err()
  end select
end function

function new_FullSubspaceBasis_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(FullSubspaceBasis)  :: this
  
  call this%read(input)
end function

impure elemental function new_FullSubspaceBasis_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(FullSubspaceBasis)       :: this
  
  this = FullSubspaceBasis(str(input))
end function
end module
