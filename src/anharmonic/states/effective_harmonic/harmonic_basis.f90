! ======================================================================
! A basis of harmonic states.
! ======================================================================
module harmonic_basis_module
  use common_module
  
  use anharmonic_common_module
  
  use harmonic_states_module
  implicit none
  
  private
  
  public :: startup_harmonic_basis
  
  public :: HarmonicBasis
  
  type, extends(SubspaceBasis) :: HarmonicBasis
    integer :: subspace_id
    integer :: supercell_size
  contains
    procedure, public, nopass :: representation => representation_HarmonicBasis
    
    procedure, public :: initial_states => initial_states_HarmonicBasis
    procedure, public :: calculate_states => calculate_states_HarmonicBasis
    procedure, public :: mode_ids => mode_ids_HarmonicBasis
    procedure, public :: paired_mode_ids => paired_mode_ids_HarmonicBasis
    
    ! Procedures involving individual states.
    ! N.B. these are all left blank, as individual harmonic states are
    !    currently treated under a different framework.
    procedure, public :: inner_product => &
                       & inner_product_HarmonicBasis
    procedure, public :: integrate_BasisState => &
                       & integrate_BasisState_HarmonicBasis
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_HarmonicBasis
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_HarmonicBasis
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_HarmonicBasis
    
    ! Procedures involving sets of states.
    procedure, public :: thermodynamic_data => &
                       & thermodynamic_data_HarmonicBasis
    procedure, public :: wavefunctions => &
                       & wavefunctions_HarmonicBasis
    procedure, public :: integrate_BasisStates => &
                       & integrate_BasisStates_HarmonicBasis
    
    ! I/O.
    procedure, public :: read  => read_HarmonicBasis
    procedure, public :: write => write_HarmonicBasis
  end type
  
  interface HarmonicBasis
    module procedure new_HarmonicBasis
    module procedure new_HarmonicBasis_Strings
    module procedure new_HarmonicBasis_StringArray
  end interface
contains

! Startup procedure and type representation.
subroutine startup_harmonic_basis()
  implicit none
  
  type(HarmonicBasis) :: basis
  
  call basis%startup()
end subroutine

impure elemental function representation_HarmonicBasis() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'harmonic'
end function

! Constructor.
impure elemental function new_HarmonicBasis(subspace_id,frequency, &
   & supercell_size) result(this)
  implicit none
  
  integer,  intent(in) :: subspace_id
  real(dp), intent(in) :: frequency
  integer,  intent(in) :: supercell_size
  type(HarmonicBasis)  :: this
  
  this%subspace_id = subspace_id
  this%frequency = frequency
  this%supercell_size = supercell_size
end function

! Calculate states.
impure elemental function initial_states_HarmonicBasis(this,subspace, &
   & thermal_energy,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicBasis),     intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: thermal_energy
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(BasisStatesPointer)             :: output
  
  output = BasisStatesPointer(HarmonicStates( subspace%id,    &
                                            & this%frequency, &
                                            & thermal_energy  ))
end function

impure elemental function calculate_states_HarmonicBasis(this,subspace,      &
   & subspace_potential,thermal_energy,state_energy_cutoff,convergence_data, &
   & anharmonic_data) result(output) 
  implicit none
  
  class(HarmonicBasis),     intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  class(PotentialData),     intent(in) :: subspace_potential
  real(dp),                 intent(in) :: thermal_energy
  real(dp),                 intent(in) :: state_energy_cutoff
  type(ConvergenceData),    intent(in) :: convergence_data
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(BasisStatesPointer)             :: output
  
  real(dp) :: energy_convergence
  
  type(NewtonRaphson) :: solver
  
  real(dp)                :: frequencies(3)
  type(ThermodynamicData) :: observables(3)
  
  real(dp) :: frequency
  
  integer :: i
  
  energy_convergence = convergence_data%energy_convergence
  
  solver = NewtonRaphson(                                  &
     & starting_value        = this%frequency,             &
     & finite_displacement   = 0.01_dp*energy_convergence, &
     & convergence_threshold = energy_convergence,         &
     & lower_bound           = 0.0_dp                      )
  i = 0
  do 
    frequencies = solver%get_inputs()
    
    observables = effective_harmonic_observables( &
         & thermal_energy  = thermal_energy,      &
         & potential       = subspace_potential,  &
         & frequency       = frequencies,         &
         & num_dimensions  = size(subspace),      &
         & supercell_size  = this%supercell_size, &
         & anharmonic_data = anharmonic_data      )
    
    call solver%set_outputs(observables%free_energy)
    
    if (solver%converged()) then
      frequency = solver%solution()
      exit
    endif
    
    i = i+1
    if (modulo(i,1000)==0) then
      call print_line(WARNING//': Newton-Raphson scheme taking a long time to &
         &converge.')
      call print_line('Iteration   : '//i)
      call print_line('Frequency   : '//frequencies)
      call print_line('Free energy : '//observables%free_energy)
    endif
  enddo
  
  output = BasisStatesPointer(HarmonicStates( this%subspace_id, &
                                            & frequency,        &
                                            & thermal_energy    ))
end function

function mode_ids_HarmonicBasis(this,subspace,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicBasis),     intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  integer, allocatable                 :: output(:)
  
  type(ComplexMode), allocatable :: modes(:)
  
  modes = subspace%modes(anharmonic_data%complex_modes)
  
  modes = modes(filter(modes%id<=modes%paired_id))
  
  output = modes%id
end function

function paired_mode_ids_HarmonicBasis(this,subspace,anharmonic_data) &
   & result(output)
  implicit none
  
  class(HarmonicBasis),     intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  integer, allocatable                 :: output(:)
  
  type(ComplexMode), allocatable :: modes(:)
  
  modes = subspace%modes(anharmonic_data%complex_modes)
  
  modes = modes(filter(modes%id<=modes%paired_id))
  
  output = modes%paired_id
end function

! Procedures involving individual states.
! N.B. these are all left blank, as individual harmonic states are
!    currently treated under a different framework.
impure elemental function inner_product_HarmonicBasis(this,bra,ket, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicBasis),     intent(in)                   :: this
  class(BasisState),        intent(in),           target :: bra
  class(BasisState),        intent(in), optional, target :: ket
  type(DegenerateSubspace), intent(in)                   :: subspace
  type(AnharmonicData),     intent(in)                   :: anharmonic_data
  real(dp)                                               :: output
  
  call print_line(CODE_ERROR//': Procedures involving individual states have &
     &not been implemented for HarmonicBasis.')
  call err()
end function

impure elemental function integrate_BasisState_HarmonicBasis(this, &
   & bra,monomial,ket,subspace,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicBasis),     intent(in)                   :: this
  class(BasisState),        intent(in),           target :: bra
  type(SparseMonomial),     intent(in)                   :: monomial
  class(BasisState),        intent(in), optional, target :: ket
  type(DegenerateSubspace), intent(in)                   :: subspace
  type(AnharmonicData),     intent(in)                   :: anharmonic_data
  complex(dp)                                            :: output
  
  call print_line(CODE_ERROR//': Procedures involving individual states have &
     &not been implemented for HarmonicBasis.')
  call err()
end function

impure elemental function kinetic_energy_HarmonicBasis(this,bra,ket, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicBasis),     intent(in)                   :: this
  class(BasisState),        intent(in),           target :: bra
  class(BasisState),        intent(in), optional, target :: ket
  type(DegenerateSubspace), intent(in)                   :: subspace
  type(AnharmonicData),     intent(in)                   :: anharmonic_data
  real(dp)                                               :: output
  
  call print_line(CODE_ERROR//': Procedures involving individual states have &
     &not been implemented for HarmonicBasis.')
  call err()
end function

impure elemental function harmonic_potential_energy_HarmonicBasis( &
   & this,bra,ket,subspace,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicBasis),     intent(in)                   :: this
  class(BasisState),        intent(in),           target :: bra
  class(BasisState),        intent(in), optional, target :: ket
  type(DegenerateSubspace), intent(in)                   :: subspace
  type(AnharmonicData),     intent(in)                   :: anharmonic_data
  real(dp)                                               :: output
  
  call print_line(CODE_ERROR//': Procedures involving individual states have &
     &not been implemented for HarmonicBasis.')
  call err()
end function

impure elemental function kinetic_stress_HarmonicBasis(this,bra,ket, &
   & subspace,stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicBasis),     intent(in)                   :: this
  class(BasisState),        intent(in),           target :: bra
  class(BasisState),        intent(in), optional, target :: ket
  type(DegenerateSubspace), intent(in)                   :: subspace
  type(StressPrefactors),   intent(in)                   :: stress_prefactors
  type(AnharmonicData),     intent(in)                   :: anharmonic_data
  type(RealMatrix)                                       :: output
  
  call print_line(CODE_ERROR//': Procedures involving individual states have &
     &not been implemented for HarmonicBasis.')
  call err()
end function

! Procedures involving sets of states.
impure elemental function thermodynamic_data_HarmonicBasis(this,    &
   & thermal_energy,states,subspace,subspace_potential,subspace_stress, &
   & stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicBasis),     intent(in)                  :: this
  real(dp),                 intent(in)                  :: thermal_energy
  class(BasisStates),       intent(in),          target :: states
  type(DegenerateSubspace), intent(in)                  :: subspace
  class(PotentialData),     intent(in)                  :: subspace_potential
  class(StressData),        intent(in), optional        :: subspace_stress
  type(StressPrefactors),   intent(in), optional        :: stress_prefactors
  type(AnharmonicData),     intent(in)                  :: anharmonic_data
  type(ThermodynamicData)                               :: output
  
  type(HarmonicStates) :: harmonic_states
  
  harmonic_states = HarmonicStates(states)
  
  output = effective_harmonic_observables(            &
     & thermal_energy  = thermal_energy,              &
     & potential       = subspace_potential,          &
     & frequency       = harmonic_states%frequency,   &
     & num_dimensions  = size(subspace),              &
     & supercell_size  = this%supercell_size,         &
     & anharmonic_data = anharmonic_data              )
end function

impure elemental function wavefunctions_HarmonicBasis(this,states, &
   & subspace,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicBasis),      intent(in)         :: this
  class(BasisStates),        intent(in), target :: states
  type(DegenerateSubspace),  intent(in)         :: subspace
  type(AnharmonicData),      intent(in)         :: anharmonic_data
  type(SubspaceWavefunctionsPointer)            :: output
  
  call err()
end function

impure elemental function integrate_BasisStates_HarmonicBasis(this, &
   & states,monomial,subspace,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicBasis),     intent(in)         :: this
  class(BasisStates),       intent(in), target :: states
  type(SparseMonomial),     intent(in)         :: monomial
  type(DegenerateSubspace), intent(in)         :: subspace
  type(AnharmonicData),     intent(in)         :: anharmonic_data
  complex(dp)                                  :: output
  
  type(HarmonicStates) :: harmonic_states
  
  harmonic_states = HarmonicStates(states)
  
  output = product(monomial%modes%harmonic_expectation( &
                      & harmonic_states%frequency,      &
                      & harmonic_states%thermal_energy, &
                      & this%supercell_size             ))
end function

! I/O.
subroutine read_HarmonicBasis(this,input)
  implicit none
  
  class(HarmonicBasis), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  integer  :: subspace_id
  real(dp) :: frequency
  integer  :: supercell_size
  
  select type(this); type is(HarmonicBasis)
    subspace_id = int(token(input(1),3))
    frequency = dble(token(input(2),3))
    supercell_size = int(token(input(3),4))
    
    this = HarmonicBasis(subspace_id, frequency, supercell_size)
  class default
    call err()
  end select
end subroutine

function write_HarmonicBasis(this) result(output)
  implicit none
  
  class(HarmonicBasis), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  select type(this); type is(HarmonicBasis)
    output = [ 'Subspace       : '//this%subspace_id,   &
             & 'Frequency      : '//this%frequency,     &
             & 'Supercell Size : '//this%supercell_size ]
  class default
    call err()
  end select
end function

function new_HarmonicBasis_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(HarmonicBasis)      :: this
  
  call this%read(input)
end function

impure elemental function new_HarmonicBasis_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(HarmonicBasis)           :: this
  
  this = HarmonicBasis(str(input))
end function
end module
