! ======================================================================
! Interfaces to the QUIP library.
! ======================================================================
module quip_wrapper_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use linear_algebra_module
  
  ! Use QUIP modules.
  use quip_unified_wrapper_module, only : quip_unified_wrapper!,initialise
  !use libatoms_module, only : Atoms,read,write
  implicit none
  
  type :: QuipResult
    real(dp)              :: energy
    real(dp), allocatable :: forces(:,:)
    real(dp)              :: virial(3,3)
  end type
contains

function call_quip(lattice,atomic_nos,positions) result(output)
  implicit none
  
  real(dp), intent(in) :: lattice(3,3)
  integer,  intent(in) :: atomic_nos(:)
  real(dp), intent(in) :: positions(:,:)
  type(QuipResult)     :: output
  
  integer :: no_atoms
  
  integer :: i,ialloc
  
  no_atoms = size(atomic_nos)
  if (size(positions,1)/=3) then
    call err()
  elseif (size(positions,2)/=no_atoms) then
    call err()
  endif
  
  allocate( output%forces(3,no_atoms), &
          & stat=ialloc); call err(ialloc)
  
  call quip_unified_wrapper( n                   = no_atoms,      &
                           & lattice             = lattice,       &
                           & z                   = atomic_nos,    &
                           & pos                 = positions,     &
                           & init_args_str       = 'IP SI_meam',  &
                           & init_args_str_len   = 10,            &
                           & energy              = output%energy, &
                           & force               = output%forces, &
                           & virial              = output%virial, &
                           & do_energy           = .true.,        &
                           & do_force            = .true.,        &
                           & do_virial           = .true.,        &
                           & quip_param_file     = 'Ti_MEAM.xml', &
                           & quip_param_file_len = 11,            &
                           & calc_args_str       = '',            &
                           & calc_args_str_len   = 0 )
end function
end module
