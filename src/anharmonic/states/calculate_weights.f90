! ======================================================================
! Calculates w=e^(-E/kT)/Z, using numerically stable strategies.
! ======================================================================
module calculate_weights_module
  use common_module
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
  
  min_energy = minloc(energies, 1)
  
  allocate(output(size(energies)), stat=ialloc); call err(ialloc)
  do i=1,size(energies)
    if (i==min_energy) then
      output(i) = 1.0_dp
    elseif (thermal_energy>1e-300_dp*(energies(i)-energies(min_energy))) then
      output(i) = exp((energies(min_energy)-energies(i))/thermal_energy)
    else
      output(i) = 0.0_dp
    endif
  enddo
  
  output = output / sum(output)
end function
end module
