! ======================================================================
! Calculates w=e^(-E/kT)/Z, using numerically stable strategies.
! ======================================================================
module caesar_calculate_weights_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: calculate_weights
contains

! W = e^((E0-E)/kT) / Z
! Z = sum(e^((E0-E)/kT))
function calculate_weights(energies,thermal_energy) result(output)
  implicit none
  
  real(dp), intent(in)  :: energies(:)
  real(dp), intent(in)  :: thermal_energy
  real(dp), allocatable :: output(:)
  
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
end function
end module
