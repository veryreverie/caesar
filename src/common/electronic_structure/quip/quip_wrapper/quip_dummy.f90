! ======================================================================
! A dummy module for when QUIP is not being linked.
! ======================================================================
module caesar_quip_wrapper_module
  use caesar_utils_module
  implicit none
  
  private
  
  public :: quip_unified_wrapper
  public :: Atoms,read,write
  public :: QUIP_LINKED
  
  logical, parameter :: QUIP_LINKED = .false.
  
  type :: Atoms
    real(dp)              :: lattice(3,3)
    integer,  allocatable :: z(:)
    real(dp), allocatable :: pos(:,:)
    real(dp), allocatable :: mass(:)
  end type
  
  interface read
    module procedure read_quip
  end interface
  
  interface write
    module procedure write_quip
  end interface
contains

subroutine quip_unified_wrapper(n,lattice,z,pos,init_args_str,           &
   & init_args_str_len,energy,force,virial,do_energy,do_force,do_virial, &
   & quip_param_file,quip_param_file_len,calc_args_str,calc_args_str_len)
  implicit none
  
  integer,      intent(in) :: n
  real(dp),     intent(in) :: lattice(:,:)
  integer,      intent(in) :: z(:)
  real(dp),     intent(in) :: pos(:,:)
  character(*), intent(in) :: init_args_str
  integer,      intent(in) :: init_args_str_len
  real(dp),     intent(in) :: energy
  real(dp),     intent(in) :: force(:,:)
  real(dp),     intent(in) :: virial(:,:)
  logical,      intent(in) :: do_energy
  logical,      intent(in) :: do_force
  logical,      intent(in) :: do_virial
  character(*), intent(in) :: quip_param_file
  integer,      intent(in) :: quip_param_file_len
  character(*), intent(in) :: calc_args_str
  integer,      intent(in) :: calc_args_str_len
  
  call print_line(ERROR//': unable to call Quip because Quip has not been &
     &linked. Please use the CMake flag LINK_TO_QUIP for Quip support.')
  call err()
end subroutine

subroutine write_quip(input,filename)
  implicit none
  
  type(Atoms),  intent(in) :: input
  character(*), intent(in) :: filename
  
  call print_line(ERROR//': unable to write .xyz file because Quip has not &
     &been linked. Please use the CMake flag LINK_TO_QUIP for Quip support.')
  call err()
end subroutine

subroutine read_quip(input,filename)
  implicit none
  
  type(Atoms),  intent(in) :: input
  character(*), intent(in) :: filename
  
  call print_line(ERROR//': unable to write .xyz file because Quip has not &
     &been linked. Please use the CMake flag LINK_TO_QUIP for Quip support.')
  call err()
end subroutine
end module
