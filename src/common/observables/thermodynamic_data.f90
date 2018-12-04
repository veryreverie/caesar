! ======================================================================
! Calculates thermodynamic properties, either for energy spectra or for
!    uncoupled harmonic oscillators.
! ======================================================================
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
  
  type, extends(Stringable) :: ThermodynamicData
    real(dp) :: thermal_energy
    real(dp) :: energy
    real(dp) :: free_energy
    real(dp) :: entropy
  contains
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
contains

impure elemental function new_ThermodynamicData(thermal_energy,energy, &
   & free_energy,entropy) result(this)
  implicit none
  
  real(dp), intent(in)    :: thermal_energy
  real(dp), intent(in)    :: energy
  real(dp), intent(in)    :: free_energy
  real(dp), intent(in)    :: entropy
  type(ThermodynamicData) :: this
  
  this%thermal_energy = thermal_energy
  this%energy = energy
  this%free_energy = free_energy
  this%entropy = entropy
end function

function new_ThermodynamicData_spectrum(thermal_energy,energies) result(output)
  implicit none
  
  real(dp), intent(in)    :: thermal_energy
  real(dp), intent(in)    :: energies(:)
  type(ThermodynamicData) :: output
  
  real(dp) :: max_energy
  real(dp) :: min_energy
  
  real(dp), allocatable :: reduced_energies(:)
  
  real(dp), allocatable :: boltzmann_factors(:)
  real(dp)              :: partition_function
  
  real(dp) :: energy
  real(dp) :: free_energy
  real(dp) :: entropy
  
  max_energy = maxval(energies)
  min_energy = minval(energies)
  
  reduced_energies = energies - min_energy
  
  ! Calculate U, F and S using reduced energies.
  if (max_energy-min_energy < 1e20_dp*thermal_energy) then
    ! Normal temperature range.
    boltzmann_factors = [exp(-reduced_energies/thermal_energy)]
    partition_function = sum(boltzmann_factors)
    energy = sum(boltzmann_factors*reduced_energies) / partition_function
    free_energy = -thermal_energy * log(partition_function)
    entropy = (energy-free_energy)/thermal_energy
  else
    ! Very low temperature limit.
    energy = 0
    free_energy = 0
    entropy = 0
  endif
  
  ! Add min_energy back in.
  energy = energy + min_energy
  free_energy = free_energy + min_energy
  
  output = ThermodynamicData(thermal_energy, energy, free_energy, entropy)
end function

function new_ThermodynamicData_harmonic(thermal_energy,frequency) &
   & result(output)
  implicit none
  
  real(dp), intent(in)    :: thermal_energy ! T.
  real(dp), intent(in)    :: frequency      ! w.
  type(ThermodynamicData) :: output
  
  real(dp) :: exp_plus
  real(dp) :: exp_minus
  
  real(dp) :: energy
  real(dp) :: free_energy
  real(dp) :: entropy
  
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
  
  output = ThermodynamicData(thermal_energy, energy, free_energy, entropy)
end function

! ----------------------------------------------------------------------
! Algebra with ThermodynamicData.
! N.B. + and - assume that the thermal_energy is the same for both arguments.
! ----------------------------------------------------------------------
impure elemental function add_ThermodynamicData_ThermodynamicData(this,that) &
   & result(output)
  implicit none
  
  type(ThermodynamicData), intent(in) :: this
  type(ThermodynamicData), intent(in) :: that
  type(ThermodynamicData)             :: output
  
  output = ThermodynamicData( this%thermal_energy,                 &
                            & this%energy      + that%energy,      &
                            & this%free_energy + that%free_energy, &
                            & this%entropy     + that%entropy      )
end function

impure elemental function negative_ThermodynamicData(input) result(output)
  implicit none
  
  type(ThermodynamicData), intent(in) :: input
  type(ThermodynamicData)             :: output
  
  output = ThermodynamicData( input%thermal_energy, &
                            & -input%energy         &
                            & -input%free_energy    &
                            & -input%entropy        )
end function

impure elemental function subtract_ThermodynamicData_ThermodynamicData(this, &
   & that) result(output)
  implicit none
  
  type(ThermodynamicData), intent(in) :: this
  type(ThermodynamicData), intent(in) :: that
  type(ThermodynamicData)             :: output
  
  output = ThermodynamicData( this%thermal_energy,                 &
                            & this%energy      - that%energy,      &
                            & this%free_energy - that%free_energy, &
                            & this%entropy     - that%entropy      )
end function

impure elemental function multiply_ThermodynamicData_integer(this,that) &
   & result(output)
  implicit none
  
  type(ThermodynamicData), intent(in) :: this
  integer,                 intent(in) :: that
  type(ThermodynamicData)             :: output
  
  output = ThermodynamicData( this%thermal_energy,     &
                            & this%energy      * that, &
                            & this%free_energy * that, &
                            & this%entropy     * that  )
end function

impure elemental function multiply_integer_ThermodynamicData(this,that) &
   & result(output)
  implicit none
  
  integer,                 intent(in) :: this
  type(ThermodynamicData), intent(in) :: that
  type(ThermodynamicData)             :: output
  
  output = ThermodynamicData( that%thermal_energy,     &
                            & this * that%energy,      &
                            & this * that%free_energy, &
                            & this * that%entropy      )
end function

impure elemental function multiply_ThermodynamicData_real(this,that) &
   & result(output)
  implicit none
  
  type(ThermodynamicData), intent(in) :: this
  real(dp),                intent(in) :: that
  type(ThermodynamicData)             :: output
  
  output = ThermodynamicData( this%thermal_energy,     &
                            & this%energy      * that, &
                            & this%free_energy * that, &
                            & this%entropy     * that  )
end function

impure elemental function multiply_real_ThermodynamicData(this,that) &
   & result(output)
  implicit none
  
  real(dp),                intent(in) :: this
  type(ThermodynamicData), intent(in) :: that
  type(ThermodynamicData)             :: output
  
  output = ThermodynamicData( that%thermal_energy,     &
                            & this * that%energy,      &
                            & this * that%free_energy, &
                            & this * that%entropy      )
end function

impure elemental function divide_ThermodynamicData_integer(this,that) &
   & result(output)
  implicit none
  
  type(ThermodynamicData), intent(in) :: this
  integer,                 intent(in) :: that
  type(ThermodynamicData)             :: output
  
  output = ThermodynamicData( this%thermal_energy,     &
                            & this%energy      / that, &
                            & this%free_energy / that, &
                            & this%entropy     / that  )
end function

impure elemental function divide_ThermodynamicData_real(this,that) &
   & result(output)
  implicit none
  
  type(ThermodynamicData), intent(in) :: this
  real(dp),                intent(in) :: that
  type(ThermodynamicData)             :: output
  
  output = ThermodynamicData( this%thermal_energy,     &
                            & this%energy      / that, &
                            & this%free_energy / that, &
                            & this%entropy     / that  )
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
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(ThermodynamicData)
    line = split_line(input)
    thermal_energy = dble(line(1))
    energy = dble(line(2))
    free_energy = dble(line(3))
    entropy = dble(line(4))
    
    this = ThermodynamicData(thermal_energy,energy,free_energy,entropy)
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
  class default
    call err()
  end select
end function

function new_ThermodynamicData_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(ThermodynamicData)  :: this
  
  call this%read(input)
end function
end module
