! ======================================================================
! Calculates thermodynamic properties of an isolated harmonic oscillator.
! ======================================================================
module harmonic_thermodynamics_module
  use common_module
  implicit none
  
  type :: ThermodynamicVariables
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
       & quantities along a negative frequency mode.')
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
end module
