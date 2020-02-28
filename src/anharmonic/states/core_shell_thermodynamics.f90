! ======================================================================
! Calculates the thermodynamic properties of a system where the core
!    states (those below a given kinetic energy cutoff) are explicitly
!    treated under VCI, and the shell states (those outside the cutoff) are
!    implicitly treated using an effective harmonic treatment.
! ======================================================================
! First, the thermodynamic data for several combinations of basis, states
!    and potential are calculated directly:
!
!  - full_harmonic  : full basis, harmonic states, harmonic potential.
!  - core_harmonic  : core basis, harmonic states, harmonic potential.
!  - full_effective : full basis, harmonic states, anharmonic potential.
!  - core_effective : core basis, harmonic states, anharmonic potential.
!  - core_vci       : core basis, vci states,      anharmonic potential.
!
! Then, the thermodynamic data is calculated indirectly for:
!  - shell_effective : shell basis, harmonic states, anharmonic potential.
! The shell basis is the complement to the core basis, i.e. the core basis
!    and the shell basis together form the full basis.
!
! In order to obtain shell_effective, core_effective must be removed from
!    full_effective. However, the two are not in equilibrium with one another,
!    since they are weighted by the harmonic treatment.
! full_harmonic and core_harmonic are used to calculate the weighting of the
!    core basis under the harmonic treatment, and then this weighting is used
!    to remove core_effective from full_effective to leave shell_effective.
!
! Finally, core_vci and shell_effective are equilibrated, to give the output.
module core_shell_thermodynamics_module
  use common_module
  
  use anharmonic_common_module
  
  use wavevector_basis_module
  implicit none
  
  private
  
  public :: core_shell_thermodynamics
contains

function core_shell_thermodynamics(thermal_energy,frequency,supercell_size, &
   & no_modes,subspace,bases,states,potential,stress,stress_prefactors,     &
   & anharmonic_data) result(output) 
  implicit none
  
  real(dp),                 intent(in)           :: thermal_energy
  real(dp),                 intent(in)           :: frequency
  integer,                  intent(in)           :: supercell_size
  integer,                  intent(in)           :: no_modes
  type(DegenerateSubspace), intent(in)           :: subspace
  type(WavevectorBasis),    intent(in)           :: bases(:)
  class(BasisStates),       intent(in)           :: states
  class(PotentialData),     intent(in)           :: potential
  class(StressData),        intent(in), optional :: stress
  type(StressPrefactors),   intent(in), optional :: stress_prefactors
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(ThermodynamicData)                        :: output
  
  type(RealMatrix), allocatable :: stress_prefactor
  
  type(ThermodynamicData) :: full_harmonic
  type(ThermodynamicData) :: core_harmonic
  type(ThermodynamicData) :: full_effective
  type(ThermodynamicData) :: core_effective
  type(ThermodynamicData) :: core_vci
  
  real(dp) :: pc
  real(dp) :: ps
  
  real(dp) :: energy
  real(dp) :: free_energy
  real(dp) :: entropy
  
  type(RealMatrix), allocatable :: stress_tensor
  real(dp),         allocatable :: volume
  real(dp),         allocatable :: enthalpy
  real(dp),         allocatable :: gibbs
  
  type(ThermodynamicData) :: shell_effective
  
  if (present(stress_prefactors)) then
    stress_prefactor = stress_prefactors%average_prefactor()
  endif
  
  ! Calculate the thermodynamic properties of the four combinations of
  !    basis, states and potential described above.
  full_harmonic = harmonic_observables(     &
     & thermal_energy   = thermal_energy,   &
     & stress           = stress,           &
     & stress_prefactor = stress_prefactor, &
     & frequency        = frequency,        &
     & num_dimensions   = no_modes,         &
     & supercell_size   = supercell_size,   &
     & anharmonic_data  = anharmonic_data   )
  
  core_harmonic = core_harmonic_observables( bases,             &
                                           & thermal_energy,    &
                                           & stress,            &
                                           & stress_prefactors, &
                                           & anharmonic_data    )
  
  full_effective = effective_harmonic_observables( &
            & thermal_energy   = thermal_energy,   &
            & potential        = potential,        &
            & stress           = stress,           &
            & stress_prefactor = stress_prefactor, &
            & frequency        = frequency,        &
            & num_dimensions   = no_modes,         &
            & supercell_size   = supercell_size,   &
            & anharmonic_data  = anharmonic_data   )
  
  core_effective = core_effective_harmonic_observables( bases,             &
                                                      & thermal_energy,    &
                                                      & potential,         &
                                                      & stress,            &
                                                      & stress_prefactors, &
                                                      & anharmonic_data    )
  
  core_vci = core_vci_observables( bases,             &
                                 & thermal_energy,    &
                                 & states,            &
                                 & subspace,          &
                                 & potential,         &
                                 & stress,            &
                                 & stress_prefactors, &
                                 & anharmonic_data    )
  
  ! Calculate the thermodynamic variables for the core VCI states equilibrated
  !    with the shell effective harmonic states.
  
  ! Check for the case where the shell is completely unoccupied in the harmonic
  !    case.
  if (core_harmonic%free_energy<=full_harmonic%free_energy) then
    output = min([full_effective, core_effective, core_vci])
    return
  endif
  
  ! Check for the very low temperature case, where the shell can be neglected
  !    entirely.
  ! e^(-690)<1e-300, below dp floating point.
  if ( full_effective%free_energy-core_vci%free_energy &
   & > thermal_energy*690                              ) then
    output = core_vci
    return
  endif
  
  ! First calculate the weighting of the core basis under the harmonic
  !    treatment.
  ! Pc = e^{-(Fc-F)/T}
  ! Ps = 1-Pc
  ! Fc and F are the free energies of the core and full bases respectively,
  !    under the harmonic treatment.
  if ( core_harmonic%free_energy-full_harmonic%free_energy &
   & > thermal_energy*1e-10_dp                             ) then
    ! Normal temperature regime.
    pc = exp( ( full_harmonic%free_energy   &
          &   - core_harmonic%free_energy ) &
          & / thermal_energy                )
    ps = 1-pc
  else
    ! High-temperature regime. O(((F-Fc)/T)^2) < 1e-20,
    !    is below numerical error.
    pc = 1 + (full_harmonic%free_energy-core_harmonic%free_energy) &
         & / thermal_energy
    ps = (core_harmonic%free_energy-full_harmonic%free_energy) &
     & / thermal_energy
  endif
  
  ! Calculate the thermodynamics of the shell basis under the effective
  !    harmonic treatment.
  ! U = PcUc + PsUs
  ! F = PcFc + PsFs +T(Pc*ln(Pc)+Ps*ln(Ps))
  ! S = PcSc + PsSs - (Pc*ln(Pc)+Ps*ln(Ps))
  !
  ! Us = U +  (U-Uc)Pc                          / Ps
  ! Fs = F + ((F-Fc)Pc -T(Pc*ln(Pc)+Ps*ln(Ps))) / Ps
  ! Ss = S + ((S-Sc)Pc + (Pc*ln(Pc)+Ps*ln(Ps))) / Ps
  !
  ! N.B. the effective entropies are the same as the corresponding
  !    harmonic entropies.
  energy = full_effective%energy &
       & + (full_effective%energy-core_effective%energy)*pc/ps
  
  entropy = full_harmonic%entropy                                &
        & + ( (full_harmonic%entropy-core_harmonic%entropy)*pc   &
        &   + pc*log(pc)+ps*log(ps)                            ) &
        & / ps
  
  free_energy = energy - thermal_energy*entropy
  
  !! Check that F in the core effective region is greater than that in the
  !!    whole region.
  !! N.B. this comparison can be swamped by numerical error if the core
  !!    region covers most of the thermally accessible space, which is checked
  !!    by the first comparison.
  !if (pc<0.99_dp) then
  !  if (full_effective%energy < core_effective%energy) then
  !    call print_line(ERROR//': Harmonic U outside of core region is less &
  !       &than that within this region. Please increase the number of &
  !       &states.')
  !    call print_line('Potential V:')
  !    call print_lines(potential)
  !    call print_line('')
  !    call print_line('Full Harmonic  U: '//full_harmonic%energy)
  !    call print_line('Core Harmonic  U: '//core_harmonic%energy)
  !    call print_line('Full Effective U: '//full_effective%energy)
  !    call print_line('Core Effective U: '//core_effective%energy)
  !    call print_line('Core VCI       U: '//core_vci%energy)
  !    call print_line('Temperature    T: '//thermal_energy)
  !    call print_line('Frequency      w: '//frequency)
  !    call err()
  !  endif
  !endif
  !
  !free_energy = full_effective%free_energy                                  &
  !          & + ( (full_effective%free_energy-core_effective%free_energy)   &
  !          &   * pc                                                        &
  !          &   - thermal_energy*(pc*log(pc)+ps*log(ps))                  ) &
  !          & / ps
  !
  !entropy = full_effective%entropy                                 &
  !      & + ( (full_effective%entropy-core_effective%entropy)*pc   &
  !      &   + pc*log(pc)+ps*log(ps)                              ) &
  !      & / ps
  
  if (present(stress)) then
    stress_tensor = full_effective%stress &
                & + (full_effective%stress-core_effective%stress)*pc/ps
    volume = anharmonic_data%structure%volume
    enthalpy = energy + volume*trace(stress_tensor)/3
    gibbs = free_energy + volume*trace(stress_tensor)/3
  endif
  
  shell_effective = ThermodynamicData( thermal_energy, &
                                     & energy,         &
                                     & free_energy,    &
                                     & entropy,        &
                                     & stress_tensor,  &
                                     & volume,         &
                                     & enthalpy,       &
                                     & gibbs           )
  
  ! Equilibrate the vci core with the effective harmonic shell.
  output = add_subsystems(core_vci, shell_effective)
end function

! ----------------------------------------------------------------------
! Combining two subsystems, removing one subsystem from another,
!    or substituting one subsystem for another.
! ----------------------------------------------------------------------

! N.B. remove_subsystem is too numerically unstable to be useful.

! Takes the thermodynamic data from two orthogonal systems in isolation,
!    and calculates the thermodynamic data of the to systems in equilibrium.
! If subsystem 1 has energy u1 and free energy f1 in isolation,
!    and subsystem 2 has energy u2 and free energy f2 in isolation,
!    then a system consisting of both subsystems in equilibrium will have:
!    U = (u1e^{-f1/T} + u2e^{-f2/T}) / (e^{-f1/T} + e^{-f2/T})
!    F = -T ln(e^{-f1/T} + e^{-f2/T})
!    H = (h1e^{-f1/T} + h2e^{-f2/T}) / (e^{-f1/T} + e^{-f2/T})
!    G = F + H - U
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
  
  type(RealMatrix), allocatable :: stress
  real(dp),         allocatable :: volume
  real(dp),         allocatable :: enthalpy
  real(dp),         allocatable :: gibbs
  
  ! The energy, free energy and entropy
  !    of the subsystem with the smaller free energy.
  real(dp) :: u1, f1, s1
  ! The energy, free energy and entropy
  !    of the subsystem with the larger free energy.
  real(dp) :: u2, f2, s2
  
  ! The stress tensors of the subsystems as defined above.
  type(RealMatrix), allocatable :: stress1, stress2
  
  logical :: calculate_stress
  
  ! Check if stress is being calculated.
  if (allocated(subsystem1%stress).neqv.allocated(subsystem2%stress)) then
    call print_line(CODE_ERROR//': Either both subsystems must have stress &
       &or neither must.')
    call err()
  endif
  calculate_stress = allocated(subsystem1%stress)
  
  ! Check that temperatures and volumes are conistent between the subspaces.
  if (abs(subsystem1%thermal_energy-subsystem2%thermal_energy)>1e-10_dp) then
    call print_line(CODE_ERROR//': Adding subsytems at different &
       &temperatures.')
    call err()
  endif
  thermal_energy = subsystem1%thermal_energy
  if (calculate_stress) then
    if (abs(subsystem1%volume-subsystem2%volume)>1e-10_dp) then
      call print_line(CODE_ERROR//': Adding subsytems at different &
         &volumes.')
      call err()
    endif
    volume = subsystem1%volume
  endif
  
  ! Identify the subsystem with the smaller free energy.
  if (subsystem1%free_energy<=subsystem2%free_energy) then
    u1 = subsystem1%energy
    f1 = subsystem1%free_energy
    s1 = subsystem1%entropy
    
    u2 = subsystem2%energy
    f2 = subsystem2%free_energy
    s2 = subsystem2%entropy
    
    if (calculate_stress) then
      stress1 = subsystem1%stress
      stress2 = subsystem2%stress
    endif
  else
    u1 = subsystem2%energy
    f1 = subsystem2%free_energy
    s1 = subsystem2%entropy
    
    u2 = subsystem1%energy
    f2 = subsystem1%free_energy
    s2 = subsystem1%entropy
    
    if (calculate_stress) then
      stress1 = subsystem2%stress
      stress2 = subsystem2%stress
    endif
  endif
  
  ! Calculate U, F, S and stress.
  if (f2-f1 > thermal_energy*690) then
    ! Very low-temperature regime. exp((f1-f2)/T)<1e-300,
    !    below dp floating point.
    energy = u1
    free_energy = f1
    entropy = s1
    if (calculate_stress) then
      stress = stress1
    endif
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
    if (calculate_stress) then
      stress = stress1 + (stress2-stress1)*exp((f1-f2)/thermal_energy)
    endif
  elseif (f2-f1 > thermal_energy*1e-10_dp) then
    ! Usual temperature regime.
    ! U = (u1 + u2e^{(f1-f2)/T}) / (1 + e^{(f1-f2)/T})
    ! F = f1 - T*ln(1 + e^{(f1-f2)/T})
    energy = (u1 + u2*exp((f1-f2)/thermal_energy)) &
         & / (1  +    exp((f1-f2)/thermal_energy))
    free_energy = f1 - thermal_energy*log(1+exp((f1-f2)/thermal_energy))
    entropy = (energy-free_energy)/thermal_energy
    if (calculate_stress) then
      stress = (stress1 + stress2*exp((f1-f2)/thermal_energy)) &
           & / (1       +         exp((f1-f2)/thermal_energy))
    endif
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
    if (calculate_stress) then
      stress = (stress1+stress2)/2 - (stress1-stress2)*(f1-f2)/thermal_energy
    endif
  endif
  
  ! Calculate H and G.
  if (calculate_stress) then
    enthalpy = energy + volume*trace(stress)/3
    gibbs = free_energy + volume*trace(stress)/3
  endif
  
  ! Construct the output.
  output = ThermodynamicData( thermal_energy, &
                            & energy,         &
                            & free_energy,    &
                            & entropy,        &
                            & stress,         &
                            & volume,         &
                            & enthalpy,       &
                            & gibbs           )
end function

!! Takes the thermodynamic data from a system, and from a subsystem of that
!!    system, and calculates the thermodynamic data of the complement subsystem
!!    in isolation.
!!    u2 = (Ue^{-F/T} - u1e^{-f1/T}) / (e^{-F/T} - e^{-f1/T})
!!    f2 = -T ln(e^{-F/T} - e^{-f1/T})
!! N.B. this assumes that the subsystem is entirely contained within the system.
!impure elemental function remove_subsystem(system,subsystem) result(output)
!  implicit none
!  
!  type(ThermodynamicData), intent(in) :: system
!  type(ThermodynamicData), intent(in) :: subsystem
!  type(ThermodynamicData)             :: output
!  
!  real(dp) :: thermal_energy
!  real(dp) :: energy
!  real(dp) :: free_energy
!  real(dp) :: entropy
!  
!  ! The energy, free energy and entropy of the system.
!  real(dp) :: u, f, s
!  ! The energy, free energy and entropy of the subsystem.
!  real(dp) :: u1, f1, s1
!  
!  if (abs(system%thermal_energy-subsystem%thermal_energy)>1e-10_dp) then
!    call print_line(CODE_ERROR//': Removing subsystem from system at &
!       &different temperature.')
!    call err()
!  endif
!  
!  u = system%energy
!  f = system%free_energy
!  s = system%entropy
!  
!  u1 = subsystem%energy
!  f1 = subsystem%free_energy
!  s1 = subsystem%entropy
!  
!  thermal_energy = system%thermal_energy
!  
!  ! Check that F<f1. If T=0, this is likely due to numerical noise.
!  if (f1<=f) then
!    if (thermal_energy<1e-10_dp .or. (f-f1)/f<0.01_dp) then
!      energy = huge(0.0_dp)
!      free_energy = huge(0.0_dp)
!      entropy = 0
!      output = ThermodynamicData(thermal_energy, energy, free_energy, entropy)
!      return
!    else
!      call print_line(CODE_ERROR//': Removing subsystem from system with &
!         &larger free energy.')
!      call err()
!    endif
!  endif
!  
!  if (f1-f > thermal_energy*690) then
!    ! Very low-temperature regime. exp((F-f1)/T)<1e-300,
!    !    below dp floating point.
!    if (f1>f) then
!      energy = u
!      free_energy = f
!      entropy = s
!    else
!      ! N.B. if f1=F, then the output subsystem is not occupied.
!      energy = huge(0.0_dp)
!      free_energy = huge(0.0_dp)
!      entropy = 0.0_dp
!    endif
!  elseif (f1-f > thermal_energy*23) then
!    ! Low-temperature regime. exp((F-f1)/T)<1e-9,
!    !    so O(exp(2(f-f1)/T)) < 1e-20 is below numerical error.
!    !
!    ! u2 = (U - u1e^{(F-f1)/T}) / (1-e^{(F-f1)/T})
!    !    = U + (U-u1)e^{(F-f1)/T} + ...
!    !
!    ! f2 = F - T*ln(1-e^{(F-f1)/T})
!    !    = F + T*e^{(F-f1)/T} + ...
!    !
!    ! s2 = (u2-f2)/T
!    !    = S + ((U-u1)/T-1)e^{(F-f1)/T} + ...
!    energy = u + (u-u1)*exp((f-f1)/thermal_energy)
!    free_energy = f + thermal_energy*exp((f-f1)/thermal_energy)
!    entropy = s1 + ((u-u1)/thermal_energy+1)*exp((f-f1)/thermal_energy)
!  elseif (f1-f > thermal_energy*1e-10_dp) then
!    ! Usual temperature regime.
!    ! u2 = (U - u1e^{(F-f1)/T}) / (1 - e^{(F-f1)/T})
!    ! f2 = F - T*ln(1 - e^{(F-f1)/T})
!    ! s2 = (u2-f2)/T
!    energy = (u - u1*exp((f-f1)/thermal_energy)) &
!         & / (1 - exp((f-f1)/thermal_energy))
!    free_energy = f - thermal_energy*log(1+exp((f-f1)/thermal_energy))
!    entropy = (energy-free_energy)/thermal_energy
!  else
!    ! High-temperature regime. O(((F-f1)/T)^2) < 1e-20,
!    !    is below numerical error.
!    ! N.B. This is unlikely to come up, since for (f1-F)/T to be small,
!    !    (f2-F)/T must be very large.
!    ! N.B. another order in (F-f1)/T is kept, since the leading term in each
!    !    case contains a factor of (U-u1)/T, which is likely to be small.
!    ! 
!    ! u2 = (U - u1 - u1(F-f1)/T - u1((F-f1)/T)^2/2 + ...) 
!    !    / (-(F-f1)/T - ((F-f1)/T)^2/2 - ((F-f1)/T)^3/6 + ...)
!    !
!    !    = ((u1-U)T/(F-f1) + u1 + u1(F-f1)/2T + ...) 
!    !    / (1 + (F-f1)/2T + ((F-f1)/T)^2/6 + ...)
!    !
!    !    = ((u1-U)T/(F-f1) + u1 + u1(F-f1)/2T + ...) 
!    !    * (1 - (F-f1)/2T + ((F-f1)/T)^2/12 + ...)
!    !
!    !    = (U+u1)/2 + (u1-U)T/(F-f1) + (u1-U)(F-f1)/12T + ...
!    !
!    ! f2 = F - T*ln(-(F-f1)/T - ((F-f1)/T)^2/2 - ((F-f1)/T)^3/6 + ...)
!    !    = F - T*ln(-(F-f1)/T) - T*ln(1 + (F-f1)/2T + ((F-f1)/T)^2/6 + ...)
!    !    = (F+f1)/2 - T*ln(-(F-f1)/T) + (F-f1)^2/24T + ...
!    !
!    ! s2 = (u2-f2)/T
!    !    = (S+s1)/2 + (u1-U)/(F-f1) + ln(-(F-f1)/T) + (2(u1-U)+(F-f1))(F-f1)/24T^2 + ...
!    !    = (S+s1)/2 - 1 + (s1-S)T/(F-f1) + (2(s1-S)-(F-f1)/T)(F-f1)/24T + ...
!    energy = (u+u1)/2                     &
!         & + (u1-u)*thermal_energy/(f-f1) &
!         & + (u1-u)*(f-f1)/(12*thermal_energy)
!    free_energy = (f+f1)/2                                  &
!              & - thermal_energy*log((f1-f)/thermal_energy) &
!              & + (f-f1)**2/(24*thermal_energy)
!    entropy = (s+s1)/2                     &
!          & - 1                            &
!          & + (s1-S)*thermal_energy/(f-f1) &
!          & + (2*(s1-s)-(f-f1)/thermal_energy)*(f-f1)/(24*thermal_energy)
!  endif
!  
!  output = ThermodynamicData(thermal_energy, energy, free_energy, entropy)
!end function
!
!! Takes the thermodynamic data from a system, a subsystem of that
!!    system, and a replacement subsystem,
!!    and calculates the thermodynamic data of the system with the first
!!    subsystem replaced by the second subsystem.
!! The system, X, has variables U and F,
!!    the subsystem to be removed, x1, has variables u1 and f1,
!!    the subsystem to be added, x', has variables u' and f',
!!    and the resulting system, X', has variables U' and F'.
!!    U'e^{-F'/T} = Ue^{-F/T} - u1e^{-f1/T} + u'e^{-f'/T}
!!      e^{-F'/T} =  e^{-F/T} -   e^{-f1/T} +   e^{-f'/T}
!! N.B. this assumes that the subsystem is entirely contained within the system,
!!    and that both subsystems cover the same basis of states.
!!
!! Noting that these equations are symmetric under exchange of X and x',
!!    it is convenient to redifine
!!          { (X , x')  if F<f'
!! (x,x2) = { 
!!          { (x', X )  if f'<F
!! Then, denoting the properties of x as u and f,
!!    and the properties of x2 as u2 and f2,
!! U' = u + (    (u2-u)e^{-(f2-f)/T} - (u1-u)e^{-(f1-f)/T})
!!        / (1 +       e^{-(f2-f)/T} -       e^{-(f1-f)/T})
!! F' = f - T*ln(1 +   e^{-(f2-f)/T} -       e^{-(f1-f)/T})
!impure elemental function substitute_subsystems(system,removed_subsystem, &
!   & added_subsystem) result(output)
!  use, intrinsic :: ieee_arithmetic
!  implicit none
!  
!  type(ThermodynamicData), intent(in) :: system
!  type(ThermodynamicData), intent(in) :: removed_subsystem
!  type(ThermodynamicData), intent(in) :: added_subsystem
!  type(ThermodynamicData)             :: output
!  
!  real(dp) :: f   ! min(system%free_energy, added_subsystem%free_energy)
!  real(dp) :: df1 ! f1-f, where f1=removed_subsystem%free_energy.
!  real(dp) :: df2 ! f2-f, where f2 is from whichever system f and f1 aren't.
!  real(dp) :: u   ! The energy of whichever system f is the free energy of.
!  real(dp) :: du1 ! u1-u, where u1=removed_subsystem%free_energy.
!  real(dp) :: du2 ! u2-u, where u2 is from whichever system u and u1 aren't.
!  real(dp) :: s   ! The entropy of whichever system f is the free energy of.
!  real(dp) :: ds1 ! s1-s, where s1=removed_subsystem%entropy.
!  real(dp) :: ds2 ! s2-s, where s2 is from whichever system u and u1 aren't.
!  
!  ! Output variables.
!  real(dp) :: thermal_energy
!  real(dp) :: energy
!  real(dp) :: free_energy
!  real(dp) :: entropy
!  
!  ! Working variables.
!  real(dp) :: exp_1 ! e^{-df1/T}
!  real(dp) :: exp_2 ! e^{-df2/T}
!  
!  ! Check that temperatures are the same.
!  if (   abs(system%thermal_energy-removed_subsystem%thermal_energy) &
!     & > 1e-10_dp                                                    ) then
!    call print_line(CODE_ERROR//': Substituting subsystems at &
!       &different temperatures.')
!    call err()
!  elseif (   abs(system%thermal_energy-added_subsystem%thermal_energy) &
!         & > 1e-10_dp                                                  ) then
!    call print_line(CODE_ERROR//': Substituting subsystems at &
!       &different temperatures.')
!    call err()
!  endif
!  
!  ! Set du1, du2, df1 and df2.
!  if (system%free_energy<added_subsystem%free_energy) then
!    u   = system%energy
!    du1 = removed_subsystem%energy - system%energy
!    du2 = added_subsystem%energy - system%energy
!    f   = system%free_energy
!    df1 = removed_subsystem%free_energy - system%free_energy
!    df2 = added_subsystem%free_energy - system%free_energy
!    s   = system%entropy
!    ds1 = removed_subsystem%entropy - system%entropy
!    ds2 = added_subsystem%entropy - system%entropy
!  else
!    u   = added_subsystem%energy
!    du1 = removed_subsystem%energy - added_subsystem%energy
!    du2 = system%energy - added_subsystem%energy
!    f   = added_subsystem%free_energy
!    df1 = removed_subsystem%free_energy - added_subsystem%free_energy
!    df2 = system%free_energy - added_subsystem%free_energy
!    s   = added_subsystem%entropy
!    ds1 = removed_subsystem%entropy - added_subsystem%entropy
!    ds2 = system%entropy - added_subsystem%entropy
!  endif
!  
!  thermal_energy = system%thermal_energy
!  
!  if (      (df1 > thermal_energy*690 .and. df2 > thermal_energy*690) &
!     & .or. thermal_energy<1e-10_dp                                   ) then
!    ! Very low-temperature regime. Both df1 and df2 are >>>> T,
!    !    or T is effectively 0.
!    !    exp(-df/T)<1e-300, below dp floating point.
!    ! The system will simply be in whichever of the system and the added
!    !    subsystem has the lowest free energy.
!    ! The removed subsystem is unimportant.
!    energy = u
!    free_energy = f
!    entropy = s
!  elseif (df1 > thermal_energy*23 .and. df2 > thermal_energy*23) then
!    ! Low-temperature regime. Both df1 and df2 are >> T.
!    !    exp(-df/T)<1e-9, so O(exp(-2df/T)) < 1e-20 is below numerical error.
!    !
!    ! U' = u +     du2*e^{-df2/T} -     du1*e^{-df1/T}
!    ! F' = f -       T*e^{-df2/T} +       T*e^{-df1/T}
!    ! S' = s + (ds2+1)*e^{-df2/T} - (ds1+1)*e^{-df1/T}
!    exp_1 = exp(-df1/thermal_energy)
!    exp_2 = exp(-df2/thermal_energy)
!    energy = u + du2*exp_2 - du1*exp_1
!    free_energy = f + thermal_energy*(exp_2-exp_1)
!    entropy = s + (ds2+1)*exp_2 - (ds1+1)*exp_1
!  elseif ( df1 > thermal_energy*1e-10_dp &
!      .or. df2 > thermal_energy*1e-10_dp ) then
!    ! Usual temperature regime. df1 and/or df2 are comparable with T.
!    ! U' = u + (  du2*e^{-df2/T}-du1*e^{-df1/T})
!    !        / (1+    e^{-df2/T}-    e^{-df1/T})
!    ! F' = f - T*ln(1+e^{-df2/T}-    e^{-df1/T})
!    ! S' = (U'-F')/T
!    exp_1 = exp(-df1/thermal_energy)
!    exp_2 = exp(-df2/thermal_energy)
!    energy = u + (du2*exp_2-du1*exp_1)/(1+exp_2-exp_1)
!    free_energy = f - thermal_energy*log(1+exp_2-exp_1)
!    entropy = (energy-free_energy)/thermal_energy
!  elseif ( df1 > thermal_energy*1e-20_dp &
!     &.or. df2 > thermal_energy*1e-20_dp ) then
!    ! High-temperature regime. Both df1 and df2 are << T.
!    !    O((df/T)^2) < 1e-20, is below numerical error.
!    ! U' = u + du2 - du1 - du2*df1/T + 2*du1*df1/T - du1*df2/T
!    ! F' = f - df2 + df1
!    ! S' = s + ds2 - ds1 - ds2*df1/T + 2*ds1*df1/T - ds1*df2/T
!    energy = u + du2 - du1 + (-du2*df1+2*du1*df1-du1*df2)/thermal_energy
!    free_energy = f + df2 - df1
!    entropy = s + ds2 - ds1 + (-ds2*df1+2*ds1*df1-ds1*df2)/thermal_energy
!  else
!    ! Very high-temperature regime. Both df1 and df2 are <<<< T.
!    !    O((df/T)) < 1e-20, is below numerical error.
!    ! U' = u + du2 - du1
!    ! F' = f + df2 - df1
!    ! S' = s + ds2 - ds1
!    energy = u + du2 - du1
!    free_energy = f + df2 - df1
!    entropy = s + ds2 - ds1
!  endif
!  
!  output = ThermodynamicData(thermal_energy, energy, free_energy, entropy)
!end function
end module
