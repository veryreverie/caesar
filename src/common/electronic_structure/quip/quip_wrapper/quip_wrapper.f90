! ======================================================================
! Interfaces to the QUIP library.
! ======================================================================
! Functionality includes:
!    - Reading and writing .xyz files.
!    - Calculating electronic structure.
module quip_wrapper_module
  ! Use modules from QUIP itself.
  use quip_unified_wrapper_module, only : quip_unified_wrapper,initialise
  use libatoms_module, only : Atoms,read,write
  implicit none
  
  private
  
  public :: quip_unified_wrapper
  public :: Atoms,read,write
  public :: QUIP_LINKED
  
  logical, parameter :: QUIP_LINKED = .true.
end module
