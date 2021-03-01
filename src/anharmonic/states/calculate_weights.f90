! ======================================================================
! Calculates w=e^(-E/kT)/Z, using numerically stable strategies.
! ======================================================================
module caesar_calculate_weights_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: calculate_weights
  
  interface
    ! W = e^((E0-E)/kT) / Z
    ! Z = sum(e^((E0-E)/kT))
    module function calculate_weights(energies,thermal_energy) result(output) 
      real(dp), intent(in)  :: energies(:)
      real(dp), intent(in)  :: thermal_energy
      real(dp), allocatable :: output(:)
    end function
  end interface
end module
