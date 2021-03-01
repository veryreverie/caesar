submodule (caesar_core_shell_thermodynamics_module) caesar_core_shell_thermodynamics_submodule
  use caesar_states_module
contains

module procedure core_shell_thermodynamics
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
end procedure

module procedure add_subsystems
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
    if ( abs(subsystem1%primitive_volume-subsystem2%primitive_volume) &
     & > 1e-10_dp                                                     ) then
      call print_line(CODE_ERROR//': Adding subsytems at different &
         &volumes.')
      call err()
    endif
    volume = subsystem1%primitive_volume
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
end procedure
end submodule
