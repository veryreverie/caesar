submodule (caesar_calculate_weights_module) caesar_calculate_weights_submodule
  use caesar_states_module
contains

module procedure calculate_weights
  integer :: min_energy
  
  integer :: i,ialloc
  
  if (size(energies)==1) then
    output = [1.0_dp]
    return
  endif
  
  min_energy = minloc(energies, 1)
  allocate(output(size(energies)), stat=ialloc); call err(ialloc)
  output = 0
  output(min_energy) = 1
  do i=1,size(energies)
    if (thermal_energy>1e-300_dp*(energies(i)-energies(min_energy))) then
      output(i) = exp((energies(min_energy)-energies(i))/thermal_energy)
    endif
  enddo
  
  output = output / sum(output)
end procedure
end submodule
