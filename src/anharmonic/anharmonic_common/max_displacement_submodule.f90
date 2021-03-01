submodule (caesar_max_displacement_module) caesar_max_displacement_submodule
  use caesar_anharmonic_common_module
contains

module procedure new_MaxDisplacement
  if (count([ present(maximum_weighted_displacement), &
            & present(frequency_of_max_displacement), &
            & present(max_energy_of_displacement)     ])<2) then
    call print_line(ERROR//': At most one argument to MaxDisplacement can be &
       &reconstructed.')
    call err()
  endif
  
  ! Use E=0.5(wu)^2 to reconstruct a missing value if needed.
  if (present(maximum_weighted_displacement)) then
    this%maximum_weighted_displacement = maximum_weighted_displacement
  else
    this%maximum_weighted_displacement = sqrt(2*max_energy_of_displacement) &
                                     & / frequency_of_max_displacement
  endif
  
  if (present(frequency_of_max_displacement)) then
    this%frequency_of_max_displacement = frequency_of_max_displacement
  else
    this%frequency_of_max_displacement = sqrt(2*max_energy_of_displacement) &
                                     & / maximum_weighted_displacement
  endif
  
  if (present(max_energy_of_displacement)) then
    this%max_energy_of_displacement = max_energy_of_displacement
  else
    this%max_energy_of_displacement = 0.5_dp                          &
                                  & * ( frequency_of_max_displacement &
                                  &   * maximum_weighted_displacement )**2
  endif
end procedure

module procedure new_MaxDisplacement_displacement
  real(dp) :: maximum_weighted_displacement
  
  ! Calculate the maximum mass-weighted displacement from the maximum
  !    displacement. This corresponds to a mode made entirely from the
  !    lightest element moving up to maximum_displacement.
  maximum_weighted_displacement = maximum_displacement &
                              & * sqrt(minval(structure%atoms%mass()))
  
  this = MaxDisplacement( maximum_weighted_displacement, &
                        & frequency_of_max_displacement, &
                        & max_energy_of_displacement     )
end procedure

module procedure read_MaxDisplacement
  real(dp) :: maximum_weighted_displacement
  real(dp) :: frequency_of_max_displacement
  real(dp) :: max_energy_of_displacement
  
  select type(this); type is(MaxDisplacement)
    maximum_weighted_displacement = dble(token(input(1),5))
    frequency_of_max_displacement = dble(token(input(2),6))
    max_energy_of_displacement = dble(token(input(3),6))
    
    this = MaxDisplacement( maximum_weighted_displacement, &
                          & frequency_of_max_displacement, &
                          & max_energy_of_displacement     )
  class default
    call err()
  end select
end procedure

module procedure write_MaxDisplacement
  select type(this); type is(MaxDisplacement)
    output = [ 'Max. weighted displacement     : '// &
             &   this%maximum_weighted_displacement, &
             & 'Frequency of max. displacement : '// &
             &   this%frequency_of_max_displacement, &
             & 'Max. energy of displacement    : '// &
             &       this%max_energy_of_displacement ]
  class default
    call err()
  end select
end procedure

module procedure new_MaxDisplacement_Strings
  call this%read(input)
end procedure

module procedure new_MaxDisplacement_StringArray
  this = MaxDisplacement(str(input))
end procedure
end submodule
