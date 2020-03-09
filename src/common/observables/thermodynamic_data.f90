! ======================================================================
! Calculates thermodynamic properties, either for energy spectra or for
!    uncoupled harmonic oscillators.
! ======================================================================
! N.B. Stress is stored with the same extensivity as U, F, S, H and G.
module thermodynamic_data_module
  use utils_module
  implicit none
  
  private
  
  public :: ThermodynamicData
  public :: calculate_bose_factor
  public :: calculate_state_weight
  
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  
  public :: sum
  public :: min
  
  type, extends(Stringable) :: ThermodynamicData
    ! T, U, F and S.
    real(dp) :: thermal_energy
    real(dp) :: energy
    real(dp) :: free_energy
    real(dp) :: entropy
    
    ! P, V, H and G.
    type(RealMatrix), allocatable :: stress
    real(dp),         allocatable :: primitive_volume
    real(dp),         allocatable :: enthalpy
    real(dp),         allocatable :: gibbs
  contains
    procedure, public :: set_stress => set_stress_ThermodynamicData
    
    procedure, public :: add_energy  => add_energy_ThermodynamicData
    procedure, public :: add_entropy => add_entropy_ThermodynamicData
    procedure, public :: add_stress  => add_stress_ThermodynamicData
    
    ! I/O.
    procedure, public :: read  => read_ThermodynamicData
    procedure, public :: write => write_ThermodynamicData
  end type
  
  interface ThermodynamicData
    module procedure new_ThermodynamicData
    module procedure new_ThermodynamicData_harmonic
    module procedure new_ThermodynamicData_spectrum
    module procedure new_ThermodynamicData_String
  end interface
  
  interface operator(+)
    module procedure add_ThermodynamicData_ThermodynamicData
  end interface
  
  interface operator(-)
    module procedure negative_ThermodynamicData
    module procedure subtract_ThermodynamicData_ThermodynamicData
  end interface
  
  interface operator(*)
    module procedure multiply_ThermodynamicData_integer
    module procedure multiply_integer_ThermodynamicData
    module procedure multiply_ThermodynamicData_real
    module procedure multiply_real_ThermodynamicData
  end interface
  
  interface operator(/)
    module procedure divide_ThermodynamicData_integer
    module procedure divide_ThermodynamicData_real
  end interface
  
  interface sum
    module procedure sum_ThermodynamicData
  end interface
  
  interface min
    module procedure min_ThermodynamicData
  end interface
contains

impure elemental function new_ThermodynamicData(thermal_energy,energy, &
   & free_energy,entropy,stress,primitive_volume,enthalpy,gibbs) result(this)
  implicit none
  
  real(dp),         intent(in)           :: thermal_energy
  real(dp),         intent(in)           :: energy
  real(dp),         intent(in)           :: free_energy
  real(dp),         intent(in)           :: entropy
  type(RealMatrix), intent(in), optional :: stress
  real(dp),         intent(in), optional :: primitive_volume
  real(dp),         intent(in), optional :: enthalpy
  real(dp),         intent(in), optional :: gibbs
  type(ThermodynamicData)                :: this
  
  logical :: arguments_present(4)
  
  arguments_present = [ present(stress),           &
                      & present(primitive_volume), &
                      & present(enthalpy),         &
                      & present(gibbs)             ]
  if (any(arguments_present) .and. .not. all(arguments_present)) then
    call print_line(CODE_ERROR//': If any of P, V, H and G are given they &
       &must all be given.')
    call err()
  endif
  
  this%thermal_energy = thermal_energy
  this%energy = energy
  this%free_energy = free_energy
  this%entropy = entropy
  
  if (any(arguments_present)) then
    this%stress = stress
    this%primitive_volume = primitive_volume
    this%enthalpy = enthalpy
    this%gibbs = gibbs
  endif
end function

function new_ThermodynamicData_spectrum(thermal_energy,energies,stresses, &
   & primitive_volume) result(output)
  implicit none
  
  real(dp),         intent(in)           :: thermal_energy
  real(dp),         intent(in)           :: energies(:)
  type(RealMatrix), intent(in), optional :: stresses(:)
  real(dp),         intent(in), optional :: primitive_volume
  type(ThermodynamicData)                :: output
  
  integer :: ground_state
  
  real(dp), allocatable :: reduced_energies(:)
  
  real(dp), allocatable :: boltzmann_factors(:)
  real(dp)              :: partition_function
  
  real(dp)                      :: energy
  real(dp)                      :: free_energy
  real(dp)                      :: entropy
  type(RealMatrix), allocatable :: stress
  real(dp),         allocatable :: enthalpy
  real(dp),         allocatable :: gibbs
  
  if (present(stresses) .neqv. present(primitive_volume)) then
    call print_line(CODE_ERROR//': Stresses and volume must either both be &
       &present or both be absent.')
    call err()
  endif
  
  ground_state = minloc(energies, 1)
  
  reduced_energies = energies - energies(ground_state)
  
  ! Calculate U, F, S and stress using reduced energies.
  if (maxval(energies)-energies(ground_state) < 1e20_dp*thermal_energy) then
    ! Normal temperature range.
    boltzmann_factors = [exp(-reduced_energies/thermal_energy)]
    partition_function = sum(boltzmann_factors)
    energy = sum(boltzmann_factors*reduced_energies) / partition_function
    free_energy = -thermal_energy * log(partition_function)
    entropy = (energy-free_energy)/thermal_energy
    if (present(stresses)) then
      stress = sum(boltzmann_factors*stresses) / partition_function
    endif
  else
    ! Very low temperature limit.
    energy = 0
    free_energy = 0
    entropy = 0
    if (present(stresses)) then
      stress = stresses(ground_state)
    endif
  endif
  
  ! Add ground state energy back in.
  energy = energy + energies(ground_state)
  free_energy = free_energy + energies(ground_state)
  
  ! Calculate stress-related quantities.
  if (present(stresses)) then
    enthalpy = energy + trace(stress)*primitive_volume/3
    gibbs = free_energy + trace(stress)*primitive_volume/3
  endif
  
  output = ThermodynamicData( thermal_energy,   &
                            & energy,           &
                            & free_energy,      &
                            & entropy,          &
                            & stress,           &
                            & primitive_volume, &
                            & enthalpy,         &
                            & gibbs             )
end function

function new_ThermodynamicData_harmonic(thermal_energy,frequency,         &
   & stress_prefactor,potential_stress,primitive_volume) result(output) 
  implicit none
  
  real(dp),         intent(in)           :: thermal_energy ! T.
  real(dp),         intent(in)           :: frequency      ! w.
  type(RealMatrix), intent(in), optional :: stress_prefactor
  type(RealMatrix), intent(in), optional :: potential_stress
  real(dp),         intent(in), optional :: primitive_volume
  type(ThermodynamicData) :: output
  
  real(dp) :: exp_plus
  real(dp) :: exp_minus
  
  real(dp) :: energy
  real(dp) :: free_energy
  real(dp) :: entropy
  
  type(RealMatrix), allocatable :: stress
  real(dp),         allocatable :: enthalpy
  real(dp),         allocatable :: gibbs
  
  logical, allocatable :: arguments_present(:)
  
  if (thermal_energy<0) then
    call print_line(CODE_ERROR//': Trying to calculate thermodynamic &
       & quantities at a negative thermal energy.')
    call err()
  elseif (frequency<0) then
    call print_line(CODE_ERROR//': Trying to calculate thermodynamic &
       &quantities along a negative frequency mode.')
    call err()
  endif
  
  ! Calculate thermodynamic quantities, using alternative strategies if
  !    in low temperature limit to avoid numerical problems.
  ! U(T) = w/2 + w/(exp(w/T)-1),
  ! F(T) = w/2 + T*ln(1-exp(-w/T))
  ! S(T) = (U(T)-F(T))/T
  !    where w is frequency and
  !          T is thermal energy (kB * temperature).
  if (frequency > 690*thermal_energy) then
    ! Very low-temperature regime. exp(-690)<1e-300, below dp floating point.
    energy      = 0.5_dp*frequency
    free_energy = 0.5_dp*frequency
    entropy     = 0.0_dp
  elseif (frequency > 23*thermal_energy) then
    ! Low-temperature regime. exp(-23)<1e-9,
    !    so O(exp(-2w/T)) < 1e-20 is below numerical error.
    ! U(T) = w/2 + w*exp(-w/T) + O(exp(-2w/T)
    ! F(T) = w/2 - T*exp(-w/T) + O(exp(-2w/T)
    ! S(T) = (w/T - 1)*exp(-w/T) + O(exp(-2w/T))
    exp_minus = exp(-frequency/thermal_energy)
    energy      = 0.5_dp*frequency + exp_minus*frequency
    free_energy = 0.5_dp*frequency - exp_minus*thermal_energy
    entropy     = exp_minus*(frequency/thermal_energy-1)
  else
    ! Usual regieme. Neither in low-temperature nor high-temperature limits.
    ! N.B. low frequencies are ignored, so the high-temperature limit is
    !    never considered.
    exp_plus = exp(frequency/thermal_energy)
    exp_minus = 1/exp_plus
    energy = frequency*(0.5_dp + 1/(exp_plus-1))
    free_energy = frequency*0.5_dp + thermal_energy*log(1-exp_minus)
    entropy = (energy-free_energy)/thermal_energy
  endif
  
  ! Calculate stress-related properties.
  arguments_present = [ present(stress_prefactor), &
                      & present(potential_stress), &
                      & present(primitive_volume)  ]
  if (any(arguments_present) .and. .not. all(arguments_present)) then
    call print_line(CODE_ERROR//': If any stress-related arguments are given, &
       &all must be given.')
    call err()
  endif
  if (all(arguments_present)) then
    ! For harmonic states, the kinetic stress is equal to UI/V,
    !    where I is the stress prefactor
    !    and V is the primitive cell volume.
    stress = potential_stress + energy*stress_prefactor/primitive_volume
    enthalpy = energy + primitive_volume*trace(stress)/3
    gibbs = free_energy + primitive_volume*trace(stress)/3
  endif
  
  output = ThermodynamicData( thermal_energy,   &
                            & energy,           &
                            & free_energy,      &
                            & entropy,          &
                            & stress,           &
                            & primitive_volume, &
                            & enthalpy,         &
                            & gibbs             )
end function

! Set the stress of a ThermodynamicData.
impure elemental subroutine set_stress_ThermodynamicData(this,stress, &
   & primitive_volume) 
  implicit none
  
  class(ThermodynamicData), intent(inout) :: this
  type(RealMatrix),         intent(in)    :: stress
  real(dp),                 intent(in)    :: primitive_volume
  
  this%stress = stress
  this%primitive_volume = primitive_volume
  this%enthalpy = this%energy + primitive_volume*trace(stress)/3
  this%gibbs = this%free_energy + primitive_volume*trace(stress)/3
end subroutine

! Add energy, entropy or stress.
impure elemental subroutine add_energy_ThermodynamicData(this,energy)
  implicit none
  
  class(ThermodynamicData), intent(inout) :: this
  real(dp),                 intent(in)    :: energy
  
  this%energy = this%energy + energy
  this%free_energy = this%free_energy + energy
  this%enthalpy = this%enthalpy + energy
  this%gibbs = this%gibbs + energy
end subroutine

impure elemental subroutine add_entropy_ThermodynamicData(this,entropy)
  implicit none
  
  class(ThermodynamicData), intent(inout) :: this
  real(dp),                 intent(in)    :: entropy
  
  this%entropy = this%entropy + entropy
  this%free_energy = this%free_energy - this%thermal_energy*entropy
  this%gibbs = this%gibbs - this%thermal_energy*entropy
end subroutine

impure elemental subroutine add_stress_ThermodynamicData(this,stress)
  implicit none
  
  class(ThermodynamicData), intent(inout) :: this
  type(RealMatrix),         intent(in)    :: stress
  
  this%stress = this%stress + stress
  this%enthalpy = this%enthalpy + trace(stress)*this%primitive_volume/3
  this%gibbs = this%gibbs + trace(stress)*this%primitive_volume/3
end subroutine

! ----------------------------------------------------------------------
! Algebra with ThermodynamicData.
! N.B. + and - assume that the thermal_energy and volume are the same
!    for both arguments.
! ----------------------------------------------------------------------
impure elemental function add_ThermodynamicData_ThermodynamicData(this,that) &
   & result(output)
  implicit none
  
  type(ThermodynamicData), intent(in) :: this
  type(ThermodynamicData), intent(in) :: that
  type(ThermodynamicData)             :: output
  
  real(dp) :: thermal_energy
  real(dp) :: energy
  real(dp) :: free_energy
  real(dp) :: entropy
  
  type(RealMatrix), allocatable :: stress
  real(dp),         allocatable :: primitive_volume
  real(dp),         allocatable :: enthalpy
  real(dp),         allocatable :: gibbs
  
  thermal_energy = this%thermal_energy
  energy = this%energy + that%energy
  free_energy = this%free_energy + that%free_energy
  entropy = this%entropy + that%entropy
  
  if (allocated(this%stress) .and. allocated(that%stress)) then
    stress = this%stress + that%stress
  endif
  
  if (       allocated(this%primitive_volume) &
     & .and. allocated(that%primitive_volume) ) then
    primitive_volume = this%primitive_volume
  endif
  
  if (allocated(this%enthalpy) .and. allocated(that%enthalpy)) then
    enthalpy = this%enthalpy + that%enthalpy
  endif
  
  if (allocated(this%gibbs) .and. allocated(that%gibbs)) then
    gibbs = this%gibbs + that%gibbs
  endif
  
  output = ThermodynamicData( thermal_energy,   &
                            & energy,           &
                            & free_energy,      &
                            & entropy,          &
                            & stress,           &
                            & primitive_volume, &
                            & enthalpy,         &
                            & gibbs             )
end function

impure elemental function negative_ThermodynamicData(input) result(output)
  implicit none
  
  type(ThermodynamicData), intent(in) :: input
  type(ThermodynamicData)             :: output
  
  real(dp) :: thermal_energy
  real(dp) :: energy
  real(dp) :: free_energy
  real(dp) :: entropy
  
  type(RealMatrix), allocatable :: stress
  real(dp),         allocatable :: primitive_volume
  real(dp),         allocatable :: enthalpy
  real(dp),         allocatable :: gibbs
  
  thermal_energy = input%thermal_energy
  energy = -input%energy
  free_energy = -input%free_energy
  entropy = -input%entropy
  
  if (allocated(input%stress)) then
    stress = -input%stress
  endif
  
  if (allocated(input%primitive_volume)) then
    primitive_volume = input%primitive_volume
  endif
  
  if (allocated(input%enthalpy)) then
    enthalpy = -input%enthalpy
  endif
  
  if (allocated(input%gibbs)) then
    gibbs = -input%gibbs
  endif
  
  output = ThermodynamicData( thermal_energy,   &
                            & energy,           &
                            & free_energy,      &
                            & entropy,          &
                            & stress,           &
                            & primitive_volume, &
                            & enthalpy,         &
                            & gibbs             )
end function

impure elemental function subtract_ThermodynamicData_ThermodynamicData(this, &
   & that) result(output)
  implicit none
  
  type(ThermodynamicData), intent(in) :: this
  type(ThermodynamicData), intent(in) :: that
  type(ThermodynamicData)             :: output
  
  real(dp) :: thermal_energy
  real(dp) :: energy
  real(dp) :: free_energy
  real(dp) :: entropy
  
  type(RealMatrix), allocatable :: stress
  real(dp),         allocatable :: primitive_volume
  real(dp),         allocatable :: enthalpy
  real(dp),         allocatable :: gibbs
  
  thermal_energy = this%thermal_energy
  energy = this%energy - that%energy
  free_energy = this%free_energy - that%free_energy
  entropy = this%entropy - that%entropy
  
  if (allocated(this%stress) .and. allocated(that%stress)) then
    stress = this%stress - that%stress
  endif
  
  if (       allocated(this%primitive_volume) &
     & .and. allocated(that%primitive_volume) ) then
    primitive_volume = this%primitive_volume
  endif
  
  if (allocated(this%enthalpy) .and. allocated(that%enthalpy)) then
    enthalpy = this%enthalpy - that%enthalpy
  endif
  
  if (allocated(this%gibbs) .and. allocated(that%gibbs)) then
    gibbs = this%gibbs - that%gibbs
  endif
  
  output = ThermodynamicData( thermal_energy,   &
                            & energy,           &
                            & free_energy,      &
                            & entropy,          &
                            & stress,           &
                            & primitive_volume, &
                            & enthalpy,         &
                            & gibbs             )
end function

impure elemental function multiply_ThermodynamicData_integer(this,that) &
   & result(output)
  implicit none
  
  type(ThermodynamicData), intent(in) :: this
  integer,                 intent(in) :: that
  type(ThermodynamicData)             :: output
  
  real(dp) :: thermal_energy
  real(dp) :: energy
  real(dp) :: free_energy
  real(dp) :: entropy
  
  type(RealMatrix), allocatable :: stress
  real(dp),         allocatable :: primitive_volume
  real(dp),         allocatable :: enthalpy
  real(dp),         allocatable :: gibbs
  
  thermal_energy = this%thermal_energy
  energy = this%energy*that
  free_energy = this%free_energy*that
  entropy = this%entropy*that
  
  if (allocated(this%stress)) then
    stress = this%stress*that
  endif
  
  if (allocated(this%primitive_volume)) then
    primitive_volume = this%primitive_volume
  endif
  
  if (allocated(this%enthalpy)) then
    enthalpy = this%enthalpy*that
  endif
  
  if (allocated(this%gibbs)) then
    gibbs = this%gibbs*that
  endif
  
  output = ThermodynamicData( thermal_energy,   &
                            & energy,           &
                            & free_energy,      &
                            & entropy,          &
                            & stress,           &
                            & primitive_volume, &
                            & enthalpy,         &
                            & gibbs             )
end function

impure elemental function multiply_integer_ThermodynamicData(this,that) &
   & result(output)
  implicit none
  
  integer,                 intent(in) :: this
  type(ThermodynamicData), intent(in) :: that
  type(ThermodynamicData)             :: output
  
  real(dp) :: thermal_energy
  real(dp) :: energy
  real(dp) :: free_energy
  real(dp) :: entropy
  
  type(RealMatrix), allocatable :: stress
  real(dp),         allocatable :: primitive_volume
  real(dp),         allocatable :: enthalpy
  real(dp),         allocatable :: gibbs
  
  thermal_energy = that%thermal_energy
  energy = that%energy*this
  free_energy = that%free_energy*this
  entropy = that%entropy*this
  
  if (allocated(that%stress)) then
    stress = that%stress*this
  endif
  
  if (allocated(that%primitive_volume)) then
    primitive_volume = that%primitive_volume
  endif
  
  if (allocated(that%enthalpy)) then
    enthalpy = that%enthalpy*this
  endif
  
  if (allocated(that%gibbs)) then
    gibbs = that%gibbs*this
  endif
  
  output = ThermodynamicData( thermal_energy,   &
                            & energy,           &
                            & free_energy,      &
                            & entropy,          &
                            & stress,           &
                            & primitive_volume, &
                            & enthalpy,         &
                            & gibbs             )
end function

impure elemental function multiply_ThermodynamicData_real(this,that) &
   & result(output)
  implicit none
  
  type(ThermodynamicData), intent(in) :: this
  real(dp),                intent(in) :: that
  type(ThermodynamicData)             :: output
  
  real(dp) :: thermal_energy
  real(dp) :: energy
  real(dp) :: free_energy
  real(dp) :: entropy
  
  type(RealMatrix), allocatable :: stress
  real(dp),         allocatable :: primitive_volume
  real(dp),         allocatable :: enthalpy
  real(dp),         allocatable :: gibbs
  
  thermal_energy = this%thermal_energy
  energy = this%energy*that
  free_energy = this%free_energy*that
  entropy = this%entropy*that
  
  if (allocated(this%stress)) then
    stress = this%stress*that
  endif
  
  if (allocated(this%primitive_volume)) then
    primitive_volume = this%primitive_volume
  endif
  
  if (allocated(this%enthalpy)) then
    enthalpy = this%enthalpy*that
  endif
  
  if (allocated(this%gibbs)) then
    gibbs = this%gibbs*that
  endif
  
  output = ThermodynamicData( thermal_energy,   &
                            & energy,           &
                            & free_energy,      &
                            & entropy,          &
                            & stress,           &
                            & primitive_volume, &
                            & enthalpy,         &
                            & gibbs             )
end function

impure elemental function multiply_real_ThermodynamicData(this,that) &
   & result(output)
  implicit none
  
  real(dp),                intent(in) :: this
  type(ThermodynamicData), intent(in) :: that
  type(ThermodynamicData)             :: output
  
  real(dp) :: thermal_energy
  real(dp) :: energy
  real(dp) :: free_energy
  real(dp) :: entropy
  
  type(RealMatrix), allocatable :: stress
  real(dp),         allocatable :: primitive_volume
  real(dp),         allocatable :: enthalpy
  real(dp),         allocatable :: gibbs
  
  thermal_energy = that%thermal_energy
  energy = that%energy*this
  free_energy = that%free_energy*this
  entropy = that%entropy*this
  
  if (allocated(that%stress)) then
    stress = that%stress*this
  endif
  
  if (allocated(that%primitive_volume)) then
    primitive_volume = that%primitive_volume
  endif
  
  if (allocated(that%enthalpy)) then
    enthalpy = that%enthalpy*this
  endif
  
  if (allocated(that%gibbs)) then
    gibbs = that%gibbs*this
  endif
  
  output = ThermodynamicData( thermal_energy,   &
                            & energy,           &
                            & free_energy,      &
                            & entropy,          &
                            & stress,           &
                            & primitive_volume, &
                            & enthalpy,         &
                            & gibbs             )
end function

impure elemental function divide_ThermodynamicData_integer(this,that) &
   & result(output)
  implicit none
  
  type(ThermodynamicData), intent(in) :: this
  integer,                 intent(in) :: that
  type(ThermodynamicData)             :: output
  
  real(dp) :: thermal_energy
  real(dp) :: energy
  real(dp) :: free_energy
  real(dp) :: entropy
  
  type(RealMatrix), allocatable :: stress
  real(dp),         allocatable :: primitive_volume
  real(dp),         allocatable :: enthalpy
  real(dp),         allocatable :: gibbs
  
  thermal_energy = this%thermal_energy
  energy = this%energy/that
  free_energy = this%free_energy/that
  entropy = this%entropy/that
  
  if (allocated(this%stress)) then
    stress = this%stress/that
  endif
  
  if (allocated(this%primitive_volume)) then
    primitive_volume = this%primitive_volume
  endif
  
  if (allocated(this%enthalpy)) then
    enthalpy = this%enthalpy/that
  endif
  
  if (allocated(this%gibbs)) then
    gibbs = this%gibbs/that
  endif
  
  output = ThermodynamicData( thermal_energy,   &
                            & energy,           &
                            & free_energy,      &
                            & entropy,          &
                            & stress,           &
                            & primitive_volume, &
                            & enthalpy,         &
                            & gibbs             )
end function

impure elemental function divide_ThermodynamicData_real(this,that) &
   & result(output)
  implicit none
  
  type(ThermodynamicData), intent(in) :: this
  real(dp),                intent(in) :: that
  type(ThermodynamicData)             :: output
  
  real(dp) :: thermal_energy
  real(dp) :: energy
  real(dp) :: free_energy
  real(dp) :: entropy
  
  type(RealMatrix), allocatable :: stress
  real(dp),         allocatable :: primitive_volume
  real(dp),         allocatable :: enthalpy
  real(dp),         allocatable :: gibbs
  
  thermal_energy = this%thermal_energy
  energy = this%energy/that
  free_energy = this%free_energy/that
  entropy = this%entropy/that
  
  if (allocated(this%stress)) then
    stress = this%stress/that
  endif
  
  if (allocated(this%primitive_volume)) then
    primitive_volume = this%primitive_volume
  endif
  
  if (allocated(this%enthalpy)) then
    enthalpy = this%enthalpy/that
  endif
  
  if (allocated(this%gibbs)) then
    gibbs = this%gibbs/that
  endif
  
  output = ThermodynamicData( thermal_energy,   &
                            & energy,           &
                            & free_energy,      &
                            & entropy,          &
                            & stress,           &
                            & primitive_volume, &
                            & enthalpy,         &
                            & gibbs             )
end function

function sum_ThermodynamicData(input) result(output)
  implicit none
  
  type(ThermodynamicData), intent(in) :: input(:)
  type(ThermodynamicData)             :: output
  
  integer :: i
  
  if (size(input)==0) then
    call print_line(ERROR//': Trying to sum an empty list.')
    call err()
  endif
  
  output = input(1)
  
  do i=2,size(input)
    output = output + input(i)
  enddo
end function

function min_ThermodynamicData(input) result(output)
  implicit none
  
  type(ThermodynamicData), intent(in) :: input(:)
  type(ThermodynamicData)             :: output
  
  output = input(minloc(input%free_energy,1))
end function

! ----------------------------------------------------------------------
! Returns n(T,w) = 1/(e^(w/T)-1), using numerically stable strategies.
! ----------------------------------------------------------------------
impure elemental function calculate_bose_factor(thermal_energy,frequency) &
   & result(output)
  implicit none
  
  real(dp), intent(in) :: thermal_energy
  real(dp), intent(in) :: frequency
  real(dp)             :: output
  
  real(dp) :: exp_minus
  real(dp) :: exp_plus
  
  if (frequency > 690*thermal_energy) then
    ! Very low-temperature regime. exp(-690)<1e-300, below dp floating point.
    output = 0.0_dp
  elseif (frequency > 23*thermal_energy) then
    ! Low-temperature regime. exp(-23)<1e-9,
    !    so O(exp(-2E/T)) < 1e-20 is below numerical error.
    ! n(T,w) = exp(-w/T)/(1-exp(-w/T)) = exp(-w/T)(1+exp(-w/T)+O(exp(-2E/T)))
    exp_minus = exp(-frequency/thermal_energy)
    output = exp_minus*(1+exp_minus)
  elseif (frequency*1e10_dp > thermal_energy) then
    ! Usual regieme. Neither in low-temperature nor high-temperature limits.
    exp_plus = exp(frequency/thermal_energy)
    output = 1/(exp_plus-1)
  else
    ! High-temperature regime. O((w/T)^2) < 1e-20 is below numerical error.
    ! n(T,w) = 1/(1+w/T+0.5(w/T)^2-1+O((w/T)^3)) = T/w-0.5+O(w/T)
    output = thermal_energy/frequency - 0.5_dp
  endif
end function

! ----------------------------------------------------------------------
! Returns P(T,w,n) = (1-e^(-w/T))^d * e^(-nw/T),
!    using numerically stable strategies.
! ----------------------------------------------------------------------
impure elemental function calculate_state_weight(thermal_energy,frequency, &
   & occupation,state_dimension) result(output)
  implicit none
  
  real(dp), intent(in) :: thermal_energy
  real(dp), intent(in) :: frequency
  integer,  intent(in) :: occupation
  integer,  intent(in) :: state_dimension
  real(dp)             :: output
  
  real(dp) :: exp_minus
  
  if (thermal_energy<0) then
    call print_line(CODE_ERROR//': Trying to calculate thermodynamic &
       & quantities at a negative thermal energy.')
    call err()
  elseif (frequency<0) then
    call print_line(CODE_ERROR//': Trying to calculate thermodynamic &
       &quantities along a negative frequency mode.')
    call err()
  elseif (occupation<0) then
    call print_line(CODE_ERROR//': Trying to calculate thermodynamic &
       &quantities using a state with negative occupation.')
    call err()
  endif
  
  if (frequency > 690*thermal_energy) then
    ! Very low-temperature regime. exp(-690)<1e-300, below dp floating point.
    if (occupation==0) then
      output = 1.0_dp
    else
      output = 0.0_dp
    endif
  else
    exp_minus = exp(-frequency/thermal_energy)
    output = (1-exp_minus)**state_dimension * exp_minus**occupation
  endif
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_ThermodynamicData(this,input)
  implicit none
  
  class(ThermodynamicData), intent(out) :: this
  type(String),             intent(in)  :: input
  
  real(dp) :: thermal_energy
  real(dp) :: energy
  real(dp) :: free_energy
  real(dp) :: entropy
  
  type(RealMatrix), allocatable :: stress
  real(dp),         allocatable :: primitive_volume
  real(dp),         allocatable :: enthalpy
  real(dp),         allocatable :: gibbs
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(ThermodynamicData)
    line = split_line(input)
    thermal_energy = dble(line(1))
    energy = dble(line(2))
    free_energy = dble(line(3))
    entropy = dble(line(4))
    
    if (size(line)>4) then
      stress = mat(dble(line(5:13)),3,3)
      primitive_volume = dble(line(14))
      enthalpy = dble(line(15))
      gibbs = dble(line(16))
    endif
    
    this = ThermodynamicData( thermal_energy,   &
                            & energy,           &
                            & free_energy,      &
                            & entropy,          &
                            & stress,           &
                            & primitive_volume, &
                            & enthalpy,         &
                            & gibbs             )
  class default
    call err()
  end select
end subroutine

function write_ThermodynamicData(this) result(output)
  implicit none
  
  class(ThermodynamicData), intent(in) :: this
  type(String)                         :: output
  
  select type(this); type is(ThermodynamicData)
    output = this%thermal_energy //' '// &
           & this%energy         //' '// &
           & this%free_energy    //' '// &
           & this%entropy
    if (allocated(this%stress)) then
      output = output                   //' '// &
             & this%stress%element(1,1) //' '// &
             & this%stress%element(1,2) //' '// &
             & this%stress%element(1,3) //' '// &
             & this%stress%element(2,1) //' '// &
             & this%stress%element(2,2) //' '// &
             & this%stress%element(2,3) //' '// &
             & this%stress%element(3,1) //' '// &
             & this%stress%element(3,2) //' '// &
             & this%stress%element(3,3) //' '// &
             & this%primitive_volume    //' '// &
             & this%enthalpy            //' '// &
             & this%gibbs
    endif
  class default
    call err()
  end select
end function

impure elemental function new_ThermodynamicData_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(ThermodynamicData)  :: this
  
  call this%read(input)
end function
end module
