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
  public :: add_subsystems
  public :: remove_subsystem
  
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  
  public :: sum
  
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
  
  interface sum
    module procedure sum_ThermodynamicData
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
! Combining two subsystems, or removing one subsystem from another.
! ----------------------------------------------------------------------

! Takes the thermodynamic data from two orthogonal systems in isolation,
!    and calculates the thermodynamic data of the to systems in equilibrium.
! If subsystem 1 has energy u1 and free energy f1 in isolation,
!    and subsystem 2 has energy u2 and free energy f2 in isolation,
!    then a system consisting of both subsystems in equilibrium will have:
!    U = (u1e^{-f1/T} + u2e^{-f2/T}) / (e^{-f1/T} + e^{-f2/T})
!    F = -T ln(e^{-f1/T} + e^{-f2/T})
! N.B. this assumes that the subsystems are orthogonal, i.e. every state in
!    subsystem 1 is orthogonal to every state in subsystem 2.
impure elemental function add_subsystems(subsystem1,subsystem2) result(output)
  implicit none
  
  type(ThermodynamicData), intent(in) :: subsystem1
  type(ThermodynamicData), intent(in) :: subsystem2
  type(ThermodynamicData)             :: output
  
  real(dp) :: thermal_energy
  real(dp) :: energy
  real(dp) :: free_energy
  real(dp) :: entropy
  
  ! The energy, free energy and entropy
  !    of the subsystem with the smaller free energy.
  real(dp) :: u1, f1, s1
  ! The energy, free energy and entropy
  !    of the subsystem with the larger free energy.
  real(dp) :: u2, f2, s2
  
  if (abs(subsystem1%thermal_energy-subsystem2%thermal_energy)>1e-10_dp) then
    call print_line(CODE_ERROR//': Adding subsytems at different &
       &temperatures.')
    call err()
  endif
  
  thermal_energy = subsystem1%thermal_energy
  
  if (subsystem1%free_energy<=subsystem2%free_energy) then
    u1 = subsystem1%energy
    f1 = subsystem1%free_energy
    s1 = subsystem1%entropy
    
    u2 = subsystem2%energy
    f2 = subsystem2%free_energy
    s2 = subsystem2%entropy
  else
    u1 = subsystem2%energy
    f1 = subsystem2%free_energy
    s1 = subsystem2%entropy
    
    u2 = subsystem1%energy
    f2 = subsystem1%free_energy
    s2 = subsystem1%entropy
  endif
  
  if (f2-f1 > thermal_energy*690) then
    ! Very low-temperature regime. exp((f1-f2)/T)<1e-300,
    !    below dp floating point.
    energy = u1
    free_energy = f1
    entropy = s1
  elseif (f2-f1 > thermal_energy*23) then
    ! Low-temperature regime. exp((f1-f2)/T)<1e-9,
    !    so O(exp(2(f1-f2)/T)) < 1e-20 is below numerical error.
    !
    ! U = (u1 + u2e^{(f1-f2)/T}) / (1+e^{(f1-f2)/T})
    !   = u1 + (u2-u1)e^{(f1-f2)/T} + ...
    !
    ! F = f1 - T*ln(1 + e^{(f1-f2)/T})
    !   = f1 - T*e^{(f1-f2)/T} + ...
    !
    ! S = (U-F)/T
    !   = s1 + ((u2-u1)/T+1)e^{(f1-f2)/T} + ...
    energy = u1 + (u2-u1)*exp((f1-f2)/thermal_energy)
    free_energy = f1 - thermal_energy*exp((f1-f2)/thermal_energy)
    entropy = s1 + ((u2-u1)/thermal_energy+1)*exp((f1-f2)/thermal_energy)
  elseif (f2-f1 > thermal_energy*1e-10_dp) then
    ! Usual temperature regime.
    ! U = (u1 + u2e^{(f1-f2)/T}) / (1 + e^{(f1-f2)/T})
    ! F = f1 - T*ln(1 + e^{(f1-f2)/T})
    energy = (u1 + u2*exp((f1-f2)/thermal_energy)) &
         & / (1  + exp((f1-f2)/thermal_energy))
    free_energy = f1 - thermal_energy*log(1+exp((f1-f2)/thermal_energy))
    entropy = (energy-free_energy)/thermal_energy
  else
    ! High-temperature regime. O(((f1-f2)/T)^2) < 1e-20,
    !    is below numerical error.
    ! U = (u1 + u2 + u2(f1-f2)/T + ...) / (2 + (f1-f2)/T + ...)
    !   = (u1+u2)/2 - (u1-u2)(f1-f2)/4T + ...
    !
    ! F = f1 - T*ln(2 + (f1-f2)/T + ...)
    !   = (f1+f2)/2 - T*ln(2) + ...
    !
    ! S = (U-F)/T
    !   = (s1+s2)/2 + ln(2) - (u1-u2)(f1-f2)/(4T^2) + ...
    energy = (u1+u2)/2 - (u1-u2)*(f1-f2)/thermal_energy
    free_energy = (f1+f2)/2 - thermal_energy*log(2.0_dp)
    entropy = (s1+s2)/2 + log(2.0_dp) - (u1-u2)*(f1-f2)/(4*thermal_energy**2)
  endif
  
  output = ThermodynamicData(thermal_energy, energy, free_energy, entropy)
end function

! Takes the thermodynamic data from a system, and from a subsystem of that
!    system, and calculates the thermodynamic data of the complement subsystem
!    in isolation.
!    u2 = (Ue^{-F/T} - u1e^{-f1/T}) / (e^{-F/T} - e^{-f1/T})
!    f2 = -T ln(e^{-F/T} - e^{-f1/T})
! N.B. this assumes that the subsystem is entirely contained within the system.
impure elemental function remove_subsystem(system,subsystem) result(output)
  implicit none
  
  type(ThermodynamicData), intent(in) :: system
  type(ThermodynamicData), intent(in) :: subsystem
  type(ThermodynamicData)             :: output
  
  real(dp) :: thermal_energy
  real(dp) :: energy
  real(dp) :: free_energy
  real(dp) :: entropy
  
  ! The energy, free energy and entropy of the system.
  real(dp) :: u, f, s
  ! The energy, free energy and entropy of the subsystem.
  real(dp) :: u1, f1, s1
  
  if (abs(system%thermal_energy-subsystem%thermal_energy)>1e-10_dp) then
    call print_line(CODE_ERROR//': Removing subsystem from system at &
       &different temperature.')
    call err()
  endif
  
  u = system%energy
  f = system%free_energy
  s = system%entropy
  
  u1 = subsystem%energy
  f1 = subsystem%free_energy
  s1 = subsystem%entropy
  
  if (f1<f) then
    call print_line(CODE_ERROR//': Removing subsystem from system with larger &
       &free energy.')
    call err()
  endif
  
  thermal_energy = system%thermal_energy
  
  if (f1-f > thermal_energy*690) then
    ! Very low-temperature regime. exp((F-f1)/T)<1e-300,
    !    below dp floating point.
    if (f1>f) then
      energy = u
      free_energy = f
      entropy = s
    else
      ! N.B. if f1=F, then the output subsystem is not occupied.
      energy = huge(0.0_dp)
      free_energy = huge(0.0_dp)
      entropy = 0.0_dp
    endif
  elseif (f1-f > thermal_energy*23) then
    ! Low-temperature regime. exp((F-f1)/T)<1e-9,
    !    so O(exp(2(f-f1)/T)) < 1e-20 is below numerical error.
    !
    ! u2 = (U - u1e^{(F-f1)/T}) / (1-e^{(F-f1)/T})
    !    = U + (U-u1)e^{(F-f1)/T} + ...
    !
    ! f2 = F - T*ln(1-e^{(F-f1)/T})
    !    = F + T*e^{(F-f1)/T} + ...
    !
    ! s2 = (u2-f2)/T
    !    = S + ((U-u1)/T-1)e^{(F-f1)/T} + ...
    energy = u + (u-u1)*exp((f-f1)/thermal_energy)
    free_energy = f + thermal_energy*exp((f-f1)/thermal_energy)
    entropy = s1 + ((u-u1)/thermal_energy+1)*exp((f-f1)/thermal_energy)
  elseif (f1-f > thermal_energy*1e-10_dp) then
    ! Usual temperature regime.
    ! u2 = (U - u1e^{(F-f1)/T}) / (1 - e^{(F-f1)/T})
    ! f2 = F - T*ln(1 - e^{(F-f1)/T})
    ! s2 = (u2-f2)/T
    energy = (u - u1*exp((f-f1)/thermal_energy)) &
         & / (1 - exp((f-f1)/thermal_energy))
    free_energy = f - thermal_energy*log(1+exp((f-f1)/thermal_energy))
    entropy = (energy-free_energy)/thermal_energy
  else
    ! High-temperature regime. O(((F-f1)/T)^2) < 1e-20,
    !    is below numerical error.
    ! N.B. This is unlikely to come up, since for (f1-F)/T to be small,
    !    (f2-F)/T must be very large.
    ! N.B. another order in (F-f1)/T is kept, since the leading term in each
    !    case contains a factor of (U-u1)/T, which is likely to be small.
    ! 
    ! u2 = (U - u1 - u1(F-f1)/T - u1((F-f1)/T)^2/2 + ...) 
    !    / (-(F-f1)/T - ((F-f1)/T)^2/2 - ((F-f1)/T)^3/6 + ...)
    !
    !    = ((u1-U)T/(F-f1) + u1 + u1(F-f1)/2T + ...) 
    !    / (1 + (F-f1)/2T + ((F-f1)/T)^2/6 + ...)
    !
    !    = ((u1-U)T/(F-f1) + u1 + u1(F-f1)/2T + ...) 
    !    * (1 - (F-f1)/2T + ((F-f1)/T)^2/12 + ...)
    !
    !    = (U+u1)/2 + (u1-U)T/(F-f1) + (u1-U)(F-f1)/12T + ...
    !
    ! f2 = F - T*ln(-(F-f1)/T - ((F-f1)/T)^2/2 - ((F-f1)/T)^3/6 + ...)
    !    = F - T*ln(-(F-f1)/T) - T*ln(1 + (F-f1)/2T + ((F-f1)/T)^2/6 + ...)
    !    = (F+f1)/2 - T*ln(-(F-f1)/T) + (F-f1)^2/24T + ...
    !
    ! s2 = (u2-f2)/T
    !    = (S+s1)/2 + (u1-U)/(F-f1) + ln(-(F-f1)/T) + (2(u1-U)+(F-f1))(F-f1)/24T^2 + ...
    !    = (S+s1)/2 - 1 + (s1-S)T/(F-f1) + (2(s1-S)-(F-f1)/T)(F-f1)/24T + ...
    energy = (u+u1)/2                     &
         & + (u1-u)*thermal_energy/(f-f1) &
         & + (u1-u)*(f-f1)/(12*thermal_energy)
    free_energy = (f+f1)/2                                  &
              & - thermal_energy*log((f1-f)/thermal_energy) &
              & + (f-f1)**2/(24*thermal_energy)
    entropy = (s+s1)/2                     &
          & - 1                            &
          & + (s1-S)*thermal_energy/(f-f1) &
          & + (2*(s1-s)-(f-f1)/thermal_energy)*(f-f1)/(24*thermal_energy)
  endif
  
  output = ThermodynamicData(thermal_energy, energy, free_energy, entropy)
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

impure elemental function new_ThermodynamicData_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(ThermodynamicData)  :: this
  
  call this%read(input)
end function
end module
