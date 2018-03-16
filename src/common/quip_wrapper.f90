! ======================================================================
! Interfaces to the QUIP library.
! ======================================================================
module quip_wrapper_module
  ! N.B. only uses precision_module, rather than all of utils_module,
  !    in order to avoid namespace collisions with Quip.
  use precision_module
  
  ! Use modules from QUIP itself.
  use quip_unified_wrapper_module, only : quip_unified_wrapper,initialise
  use libatoms_module, only : Atoms,read,write
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
  contains
    procedure, private :: check => check_QuipAtoms
  end type
contains

subroutine check_QuipAtoms(this)
  implicit none
  
  class(QuipAtoms), intent(in) :: this
  
  if (.not. allocated(this%atomic_nos)) then
    write(*,*) 'QuipAtoms error: atomic_nos not allocated.'
    stop
  elseif (.not. allocated(this%positions)) then
    write(*,*) 'QuipAtoms error: positions not allocated.'
    stop
  elseif (size(this%positions,1)/=3) then
    write(*,*) 'QuipAtoms error: size(positions,1)/=3'
    stop
  elseif (size(this%positions,2)/=size(this%atomic_nos)) then
    write(*,*) 'QuipAtoms error: size(positions,2)/=size(atomic_nos)'
    stop
  endif
end subroutine

function call_quip(quip_atoms,quip_filename) result(output)
  implicit none
  
  type(QuipAtoms), intent(in) :: quip_atoms
  character(*),    intent(in) :: quip_filename
  type(QuipResult)            :: output
  
  integer :: ialloc
  
  call quip_atoms%check()
  
  allocate(output%forces(3,size(quip_atoms%atomic_nos)), stat=ialloc)
  if (ialloc/=0) then
    write(*,*) 'Quip wrapper error: error allocating forces.'
    stop
  endif
  
  call quip_unified_wrapper(                              &
     & n                   = size(quip_atoms%atomic_nos), &
     & lattice             = quip_atoms%lattice,          &
     & z                   = quip_atoms%atomic_nos,       &
     & pos                 = quip_atoms%positions,        &
     & init_args_str       = 'IP SI_meam',                &
     & init_args_str_len   = 10,                          &
     & energy              = output%energy,               &
     & force               = output%forces,               &
     & virial              = output%virial,               &
     & do_energy           = .true.,                      &
     & do_force            = .true.,                      &
     & do_virial           = .true.,                      &
     & quip_param_file     = quip_filename,               &
     & quip_param_file_len = len(quip_filename),          &
     & calc_args_str       = '',                          &
     & calc_args_str_len   = 0 )
end function

function read_quip_atoms(filename) result(output)
  implicit none
  
  character(*), intent(in)  :: filename
  type(QuipAtoms)           :: output
  
  type(Atoms) :: quip_atoms
  integer     :: error
  
  call read(quip_atoms, filename, error=error)
  
  output%lattice    = quip_atoms%lattice
  output%atomic_nos = quip_atoms%z
  output%positions  = quip_atoms%pos
end function

subroutine write_quip_atoms(filename,this)
  implicit none
  
  character(*),    intent(in) :: filename
  type(QuipAtoms), intent(in) :: this
  
  type(Atoms) :: quip_atoms
  
  call this%check()
  
  quip_atoms%lattice = this%lattice
  quip_atoms%z       = this%atomic_nos
  quip_atoms%pos     = this%positions
  
  call write(quip_atoms, filename)
end subroutine
end module
