! ======================================================================
! An abstract class to hold a printable representation of the wavefunctions
!    spanning a subspace.
! ======================================================================
module subspace_wavefunctions_module
  use common_module
  implicit none
  
  private
  
  public :: SubspaceWavefunctions
  
  type, abstract, extends(Stringsable) :: SubspaceWavefunctions
  end type
end module
