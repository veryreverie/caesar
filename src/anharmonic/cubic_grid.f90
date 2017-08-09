! ======================================================================
! Anharmonic calculations for cubic grids.
! ======================================================================
module cubic_grid_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! Takes a set of basis functions and the B-O landscape on a cubic grid,
!    and returns the coulings (<i|j>) between the basis functions.
function harmonic_couplings_cubic_grid(harmonic_states,coupling, &
   & sampling_points,energy,forces) result(harmonic_couplings)
  use coupling_module
  use sampling_points_module
  use linear_algebra_module
  use harmonic_states_module
  implicit none
  
  type(HarmonicState), intent(in) :: harmonic_states(:,:)
  type(CoupledModes),  intent(in) :: coupling(:)
  type(SamplingPoint), intent(in) :: sampling_points(:)
  real(dp),            intent(in) :: energy(:)
  type(RealVector),    intent(in) :: forces(:,:)
  real(dp), allocatable           :: harmonic_couplings(:,:,:)
end function
end module
