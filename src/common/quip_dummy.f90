! ======================================================================
! A dummy module for when QUIP is not being linked against.
! ======================================================================
module quip_wrapper_module
  use utils_module
  implicit none
  
  private
  
  public :: QuipAtoms
  public :: QuipResult
  public :: call_quip
  public :: read_quip_atoms
  public :: write_quip_atoms
  
  type :: QuipResult
    real(dp)              :: energy
    real(dp), allocatable :: forces(:,:)
    real(dp)              :: virial(3,3)
  end type
  
  type :: QuipAtoms
    real(dp)              :: lattice(3,3)
    integer,  allocatable :: atomic_nos(:)
    real(dp), allocatable :: positions(:,:)
  end type
contains

function call_quip(quip_atoms,quip_filename) result(output)
  implicit none
  
  type(QuipAtoms), intent(in) :: quip_atoms
  character(*),    intent(in) :: quip_filename
  type(QuipResult)            :: output
  
  call print_line(ERROR//': Caesar has not been linked against QUIP.')
  stop
end function

function read_quip_atoms(filename) result(output)
  implicit none
  
  character(*), intent(in)  :: filename
  type(QuipAtoms)           :: output
  
  call print_line(ERROR//': Caesar has not been linked against QUIP.')
  stop
end function

subroutine write_quip_atoms(filename,this)
  implicit none
  
  character(*),    intent(in) :: filename
  type(QuipAtoms), intent(in) :: this
  
  call print_line(ERROR//': Caesar has not been linked against QUIP.')
  stop
end subroutine
end module
