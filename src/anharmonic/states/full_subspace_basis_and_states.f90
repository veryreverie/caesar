! ======================================================================
! A basis of states which spans a full subspace.
! ======================================================================
module full_subspace_basis_and_states_module
  use common_module
  
  use anharmonic_common_module
  
  use state_helper_module
  use monomial_state_module
  use harmonic_state_module
  use polynomial_state_module
  use braket_module
  use state_conversion_module
  use wavevector_basis_module
  implicit none
  
  private
  
  public :: FullSubspaceBasis
  public :: FullSubspaceState
  public :: FullSubspaceStates
  
  public :: size
  public :: PolynomialState
  public :: initial_ground_state
  
  ! All states spanning the subspace.
  type, extends(SubspaceBasis) :: FullSubspaceBasis
    ! The maximum power of the monomial states.
    ! This is also the maximum occupation of the harmonic basis states.
    integer  :: maximum_power
    ! The expansion order of the potential.
    ! This is also the limit on coupling between states.
    integer :: expansion_order
    ! The ID and frequency of the subspace.
    integer  :: subspace_id
    real(dp) :: frequency
    ! The states, wavevector by wavevector.
    ! N.B. this only includes one wavevector from each symmetry-related set.
    type(WavevectorBasis), allocatable :: wavevectors(:)
  contains
    ! Set the frequency of the basis.
    procedure, public :: set_frequency => set_frequency_FullSubspaceBasis
    
    ! Print the harmonic ground-state wavefunction of the basis.
    procedure, public :: ground_state_wavefunction
    
    ! Generate the first guess at states.
    procedure, public :: initial_states => initial_states_FullSubspaceBasis
    
    ! Generate the eigenstates of a single-subspace potential.
    procedure, public :: calculate_states => calculate_states_FullSubspaceBasis
    
    ! I/O.
    procedure, public :: read  => read_FullSubspaceBasis
    procedure, public :: write => write_FullSubspaceBasis
  end type
  
  interface FullSubspaceBasis
    module procedure new_FullSubspaceBasis
    module procedure new_FullSubspaceBasis_subspace
    module procedure new_FullSubspaceBasis_Strings
    module procedure new_FullSubspaceBasis_StringArray
  end interface
  
  interface size
    module procedure size_FullSubspaceBasis
  end interface
  
  type, extends(Stringsable) :: FullSubspaceState
    integer               :: subspace_id
    type(FractionVector)  :: wavevector
    integer               :: degeneracy
    real(dp)              :: energy
    real(dp), allocatable :: coefficients(:)
  contains
    procedure, public :: wavefunction => wavefunction_FullSubspaceState
    ! I/O.
    procedure, public :: read  => read_FullSubspaceState
    procedure, public :: write => write_FullSubspaceState
  end type
  
  interface FullSubspaceState
    module procedure new_FullSubspaceState
    module procedure new_FullSubspaceState_Strings
    module procedure new_FullSubspaceState_StringArray
  end interface
  
  interface PolynomialState
    module procedure new_PolynomialState_FullSubspaceState
  end interface
  
  type, extends(SubspaceStates) :: FullSubspaceStates
    type(FullSubspaceState), allocatable :: vscf_states(:)
  contains
    procedure, public :: states => states_FullSubspaceStates
    procedure, public :: spectra => spectra_FullSubspaceStates
    procedure, public :: wavefunctions => wavefunctions_FullSubspaceStates
    procedure, public :: integrate => integrate_FullSubspaceStates
  end type
  
  interface FullSubspaceStates
    module procedure new_FullSubspaceStates
  end interface
contains

! ----------------------------------------------------------------------
! FullSubspaceBasis methods.
! ----------------------------------------------------------------------
! Constructors and size functions.
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

function size_FullSubspaceBasis(this) result(output)
  implicit none
  
  type(FullSubspaceBasis), intent(in) :: this
  integer                             :: output
  
  output = size(this%wavevectors)
end function

! Set the frequency of the basis.
impure elemental subroutine set_frequency_FullSubspaceBasis(this,frequency)
  implicit none
  
  class(FullSubspaceBasis), intent(inout) :: this
  real(dp),                 intent(in)    :: frequency
  
  this%frequency = frequency
  call this%wavevectors%set_frequency(frequency)
end subroutine

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
  
  integer :: i
  
  ! Calculate the (geometric) average mass.
  mass = product(supercell%atoms%mass())
  mass = mass**(1.0_dp/size(supercell%atoms))
  
  ! Calculate the coefficient.
  coefficient = 1
  terms = [String::]
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
   & anharmonic_data) result(output)
  implicit none
  
  class(FullSubspaceBasis), intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  class(SubspaceStates), allocatable   :: output
  
  type(FullSubspaceState)  :: ground_state
  type(FullSubspaceStates) :: states
  
  integer :: ialloc
  
  ! Generate the state |0>.
  ground_state = initial_ground_state(this)
  
  ! Generate the set of states {|0>}.
  states = FullSubspaceStates([ground_state])
  
  ! Convert the set of states to the abstract type SubspaceStates.
  allocate(output, source=states, stat=ialloc); call err(ialloc)
end function

! Calculate the eigenstates of a single-subspace potential.
impure elemental function calculate_states_FullSubspaceBasis(this,subspace, &
   & subspace_potential,anharmonic_data) result(output)
  implicit none
  
  class(FullSubspaceBasis), intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  class(PotentialData),     intent(in) :: subspace_potential
  type(AnharmonicData),     intent(in) :: anharmonic_data
  class(SubspaceStates), allocatable   :: output
  
  ! Variables for constructing the basis.
  type(StructureData) :: supercell
  
  ! Variables for calculating and diagonalising the VSCF Hamiltonian.
  integer                                :: no_states
  type(HarmonicState)                    :: bra
  type(HarmonicState)                    :: ket
  real(dp),                  allocatable :: hamiltonian(:,:)
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  
  type(FullSubspaceState)              :: vscf_state
  type(FullSubspaceState), allocatable :: vscf_states(:)
  
  integer :: i,j,k,l,ialloc
  
  supercell = anharmonic_data%anharmonic_supercell
  
  vscf_states = [FullSubspaceState::]
  do i=1,size(this)
    no_states = size(this%wavevectors(i))
    allocate( hamiltonian(no_states,no_states), &
            & stat=ialloc); call err(ialloc)
    hamiltonian = 0
    do j=1,no_states
      bra = this%wavevectors(i)%harmonic_states(j)
      do k=1,size(this%wavevectors(i)%harmonic_couplings(j))
        l = this%wavevectors(i)%harmonic_couplings(j)%id(k)
        ket = this%wavevectors(i)%harmonic_states(l)
        
        hamiltonian(j,l) = kinetic_energy(bra, ket, subspace, supercell) &
                       & + subspace_potential%potential_energy(          &
                       &                        bra,                     &
                       &                        ket,                     &
                       &                        subspace,                &
                       &                        anharmonic_data )
      enddo
    enddo
    
    estuff = diagonalise_symmetric(hamiltonian)
    do j=1,size(estuff)
      vscf_state = FullSubspaceState(                     &
         & subspace_id  = this%subspace_id,               &
         & wavevector   = this%wavevectors(i)%wavevector, &
         & degeneracy   = this%wavevectors(i)%degeneracy, &
         & energy       = estuff(j)%eval,                 &
         & coefficients = estuff(j)%evec                  )
      
      ! Make energies extensive rather than intensive.
      vscf_state%energy = vscf_state%energy * supercell%sc_size
      
      vscf_states = [vscf_states, vscf_state]
    enddo
    deallocate(hamiltonian, stat=ialloc); call err(ialloc)
  enddo
  
  allocate( output, source=FullSubspaceStates(vscf_states), &
          & stat=ialloc); call err(ialloc)
end function

! ----------------------------------------------------------------------
! FullSubspaceState methods.
! ----------------------------------------------------------------------
! Constructor.
function new_FullSubspaceState(subspace_id,wavevector,degeneracy,energy, &
   & coefficients) result(this)
  implicit none
  
  integer,              intent(in) :: subspace_id
  type(FractionVector), intent(in) :: wavevector
  integer,              intent(in) :: degeneracy
  real(dp),             intent(in) :: energy
  real(dp),             intent(in) :: coefficients(:)
  type(FullSubspaceState)          :: this
  
  this%subspace_id  = subspace_id
  this%wavevector   = wavevector
  this%degeneracy   = degeneracy
  this%energy       = energy
  this%coefficients = coefficients
end function

! Construct a PolynomialState from a FullSubspaceState.
impure elemental function new_PolynomialState_FullSubspaceState(state,basis) &
   & result(output)
  implicit none
  
  type(FullSubspaceState), intent(in) :: state
  type(FullSubspaceBasis), intent(in) :: basis
  type(PolynomialState)               :: output
  
  real(dp), allocatable :: coefficients(:)
  
  integer :: i
  
  i = first(basis%wavevectors%wavevector==state%wavevector)
  
  coefficients = &
     & basis%wavevectors(i)%coefficients_basis_to_states(state%coefficients)
  
  output = PolynomialState( basis%subspace_id,                    &
                          & basis%wavevectors(i)%monomial_states, &
                          & coefficients                          )
end function

! Generate initial guess. This is simply the basis state |0>, i.e. the
!    ground state of the effective harmonic potential from which the basis
!    states were generated.
impure elemental function initial_ground_state(basis) result(output)
  implicit none
  
  type(FullSubspaceBasis), intent(in) :: basis
  type(FullSubspaceState)             :: output
  
  type(WavevectorBasis) :: wavevector_basis
  
  real(dp), allocatable :: coefficients(:)
  
  integer :: i
  
  ! Find the wavevector [0,0,0].
  wavevector_basis = basis%wavevectors(            &
     & first(is_int(basis%wavevectors%wavevector)) )
  
  ! Construct the coefficient vector in the basis of monomial states.
  ! All coefficients are zero, except for the coefficient of |0>, which is one.
  coefficients = [( 0.0_dp, i=1, size(wavevector_basis) )]
  coefficients(                                                      &
     & first(wavevector_basis%harmonic_states%total_occupation()==0) ) = 1
  
  ! Convert the coefficients into the orthonormal basis.
  coefficients = wavevector_basis%coefficients_states_to_basis(coefficients)
  
  ! Construct output.
  output = FullSubspaceState(                      &
     & subspace_id  = basis%subspace_id,           &
     & wavevector   = wavevector_basis%wavevector, &
     & degeneracy   = wavevector_basis%degeneracy, &
     & energy       = 0.0_dp,                      &
     & coefficients = coefficients                 )
end function

impure elemental function wavefunction_FullSubspaceState(this,basis, &
   & supercell) result(output)
  implicit none
  
  class(FullSubspaceState), intent(in) :: this
  type(FullSubspaceBasis),  intent(in) :: basis
  type(StructureData),      intent(in) :: supercell
  type(String)                         :: output
  
  type(WavevectorBasis) :: wavevector_basis
  
  real(dp),     allocatable :: coefficients(:)
  type(String), allocatable :: terms(:)
  
  type(String) :: state
  
  integer :: i,ialloc
  
  wavevector_basis = basis%wavevectors(                     &
     & first(basis%wavevectors%wavevector==this%wavevector) )
  
  coefficients = wavevector_basis%coefficients_basis_to_states( &
                                            & this%coefficients )
  allocate(terms(size(coefficients)), stat=ialloc); call err(ialloc)
  do i=1,size(coefficients)
    state = wavevector_basis%monomial_states(i)%wavefunction( &
                                           & basis%frequency, &
                                           & supercell        )
    ! Trim the '|0>' from the state.
    state = slice(state,1,len(state)-3)
    
    terms(i) = state
  enddo
  output = '('//join(terms,' + ')//')|0>'
end function

! ----------------------------------------------------------------------
! FullSubspaceStates methods.
! ----------------------------------------------------------------------
! Constructor.
function new_FullSubspaceStates(vscf_states) result(this)
  implicit none
  
  type(FullSubspaceState), intent(in) :: vscf_states(:)
  type(FullSubspaceStates)            :: this
  
  this%vscf_states = vscf_states
end function

! States.
function states_FullSubspaceStates(this,subspace,subspace_basis, &
   & anharmonic_data) result(output)
  implicit none
  
  class(FullSubspaceStates), intent(in) :: this
  type(DegenerateSubspace),  intent(in) :: subspace
  class(SubspaceBasis),      intent(in) :: subspace_basis
  type(AnharmonicData),      intent(in) :: anharmonic_data
  class(SubspaceState), allocatable     :: output(:)
  
  type(PolynomialState), allocatable :: states(:)
  
  integer :: ialloc
  
  select type(subspace_basis); type is(FullSubspaceBasis)
    states = PolynomialState(this%vscf_states, subspace_basis)
  class default
    call err()
  end select
  
  allocate(output, source=states, stat=ialloc); call err(ialloc)
end function

! Energy spectra.
impure elemental function spectra_FullSubspaceStates(this,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(FullSubspaceStates), intent(in) :: this
  type(DegenerateSubspace),  intent(in) :: subspace
  class(SubspaceBasis),      intent(in) :: subspace_basis
  type(AnharmonicData),      intent(in) :: anharmonic_data
  type(EnergySpectra)                   :: output
  
  output = EnergySpectra([EnergySpectrum( this%vscf_states%energy,    &
                                        & this%vscf_states%degeneracy )])
end function

! Wavefunctions.
impure elemental function wavefunctions_FullSubspaceStates(this,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(FullSubspaceStates), intent(in)     :: this
  type(DegenerateSubspace),  intent(in)     :: subspace
  class(SubspaceBasis),      intent(in)     :: subspace_basis
  type(AnharmonicData),      intent(in)     :: anharmonic_data
  class(SubspaceWavefunctions), allocatable :: output
  
  type(String)                    :: ground_state
  type(String), allocatable       :: state_wavefunctions(:)
  type(FullSubspaceWavefunctions) :: wavefunctions
  
  integer :: ialloc
  
  select type(subspace_basis); type is(FullSubspaceBasis)
    ! Construct the wavefunction of |0>.
    ground_state = subspace_basis%ground_state_wavefunction( &
                      & subspace,                            &
                      & anharmonic_data%anharmonic_supercell )
    
    ! Construct the wavefunctions of each state, in terms of |0>.
    state_wavefunctions = this%vscf_states%wavefunction( &
                  & subspace_basis,                      &
                  & anharmonic_data%anharmonic_supercell )
    
    ! Pack the wavefunctions into the output,
    !    of concrete type FullSubspaceWavefunctions.
    wavefunctions = FullSubspaceWavefunctions( subspace%id,                 &
                                             & subspace%mode_ids,           &
                                             & subspace%paired_ids,         &
                                             & ground_state,                &
                                             & this%vscf_states%energy,     &
                                             & this%vscf_states%degeneracy, &
                                             & state_wavefunctions          )
    
    ! Convert the output to abstract class SubspaceWavefunctions.
    allocate(output, source=wavefunctions, stat=ialloc); call err(ialloc)
  class default
    call print_line(CODE_ERROR//': Wrong kind of basis passed to &
       &FullSubspaceStates. Basis kind must be FullSubspaceBasis.')
    call err()
  end select
end function

! Integrate a potential.
impure elemental function integrate_FullSubspaceStates(this,potential, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(FullSubspaceStates), intent(in) :: this
  class(PotentialData),      intent(in) :: potential
  type(DegenerateSubspace),  intent(in) :: subspace
  class(SubspaceBasis),      intent(in) :: subspace_basis
  type(AnharmonicData),      intent(in) :: anharmonic_data
  class(PotentialData), allocatable     :: output
  
  type(FullSubspaceState) :: ground_state
  type(PolynomialState)   :: polynomial_state
  
  integer :: i,ialloc
  
  select type(subspace_basis); type is(FullSubspaceBasis)
    ! Identify the ground state.
    ground_state = this%vscf_states(minloc(this%vscf_states%energy,1))
    polynomial_state = PolynomialState(ground_state, subspace_basis)
    
    ! Braket the potential between the ground state.
    output = potential%braket( polynomial_state, &
                             & polynomial_state, &
                             & subspace,         &
                             & anharmonic_data   )
  class default
    call err()
  end select
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

subroutine read_FullSubspaceState(this,input)
  implicit none
  
  class(FullSubspaceState), intent(out) :: this
  type(String),             intent(in)  :: input(:)
  
  integer               :: subspace_id
  type(FractionVector)  :: wavevector
  integer               :: degeneracy
  real(dp)              :: energy
  real(dp), allocatable :: coefficients(:)
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(FullSubspaceState)
    line = split_line(input(1))
    subspace_id = int(line(3))
    
    line = split_line(input(2))
    wavevector = FractionVector(join(line(3:5)))
    
    line = split_line(input(3))
    degeneracy = int(line(3))
    
    line = split_line(input(4))
    energy = dble(line(3))
    
    line = split_line(input(5))
    coefficients = dble(line(3:))
    
    this = FullSubspaceState( subspace_id, &
                            & wavevector,  &
                            & degeneracy,  &
                            & energy,      &
                            & coefficients )
  class default
    call err()
  end select
end subroutine

function write_FullSubspaceState(this) result(output)
  implicit none
  
  class(FullSubspaceState), intent(in) :: this
  type(String), allocatable            :: output(:)
  
  select type(this); type is(FullSubspaceState)
    output = [ 'Subspace     : '//this%subspace_id, &
             & 'Wavevector   : '//this%wavevector,  &
             & 'Degeneracy   : '//this%degeneracy,  &
             & 'Energy       : '//this%energy,      &
             & 'Coefficients : '//this%coefficients ]
  class default
    call err()
  end select
end function

function new_FullSubspaceState_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(FullSubspaceState)  :: this
  
  call this%read(input)
end function

impure elemental function new_FullSubspaceState_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(FullSubspaceState)       :: this
  
  this = FullSubspaceState(str(input))
end function
end module