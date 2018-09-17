! ======================================================================
! Calculates thermodynamic properties of an isolated harmonic oscillator.
! ======================================================================
module harmonic_thermodynamics_module
  use utils_module
  implicit none
  
  private
  
  public :: ThermodynamicVariables
  public :: calculate_bose_factor
  public :: calculate_state_weight
  
  type, extends(NoDefaultConstructor) :: ThermodynamicVariables
    real(dp) :: energy
    real(dp) :: free_energy
    real(dp) :: entropy
  end type
  
  interface ThermodynamicVariables
    module procedure new_ThermodynamicVariables
  end interface
contains

impure elemental function new_ThermodynamicVariables(thermal_energy,frequency) result(output)
  implicit none
  
  real(dp), intent(in)         :: thermal_energy ! T.
  real(dp), intent(in)         :: frequency      ! w.
  type(ThermodynamicVariables) :: output
  
  real(dp) :: exp_plus
  real(dp) :: exp_minus
  
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
    output%energy      = 0.5_dp*frequency
    output%free_energy = 0.5_dp*frequency
    output%entropy     = 0.0_dp
  elseif (frequency > 23*thermal_energy) then
    ! Low-temperature regime. exp(-23)<1e-9,
    !    so O(exp(-2w/T)) < 1e-20 is below numerical error.
    ! U(T) = w/2 + w*exp(-w/T) + O(exp(-2w/T)
    ! F(T) = w/2 - T*exp(-w/T) + O(exp(-2w/T)
    ! S(T) = (w/T - 1)*exp(-w/T) + O(exp(-2w/T))
    exp_minus = exp(-frequency/thermal_energy)
    output%energy      = 0.5_dp*frequency + exp_minus*frequency
    output%free_energy = 0.5_dp*frequency - exp_minus*thermal_energy
    output%entropy     = exp_minus*(frequency/thermal_energy-1)
  else
    ! Usual regieme. Neither in low-temperature nor high-temperature limits.
    ! N.B. low frequencies are ignored, so the high-temperature limit is
    !    never considered.
    exp_plus = exp(frequency/thermal_energy)
    exp_minus = 1/exp_plus
    output%energy = frequency*(0.5_dp + 1/(exp_plus-1))
    output%free_energy = frequency*0.5_dp + thermal_energy*log(1-exp_minus)
    output%entropy = (output%energy-output%free_energy)/thermal_energy
  endif
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
  elseif (frequency*1e10 > thermal_energy) then
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
! Returns P(T,w,n) = (1-e^(-w/T))e^(-nw/T),
!    using numerically stable strategies.
! ----------------------------------------------------------------------
impure elemental function calculate_state_weight(thermal_energy,frequency, &
   & occupation) result(output)
  implicit none
  
  real(dp), intent(in) :: thermal_energy
  real(dp), intent(in) :: frequency
  integer,  intent(in) :: occupation
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
    output = (1-exp_minus) * exp_minus**occupation
  endif
end function
end module
