module supercell_module
  implicit none
  
  ! holds information about each supercell
  type Supercell
    integer :: no_atoms
    integer :: no_atoms_sc
    integer :: no_modes
    integer :: no_cells
  end type
  
contains

! Makes a Supercell from two files ('equilibrium.dat and super_equilibrium.dat')
function read_supercell(eqm_file_unit, super_eqm_file_unit) result(output)
  implicit none
  
  integer, intent(in) :: eqm_file_unit
  integer, intent(in) :: super_eqm_file_unit
  type(Supercell)     :: output
  
  read(eqm_file_unit,*) output%no_atoms
  read(super_eqm_file_unit,*) output%no_atoms_sc
  
  output%no_modes = output%no_atoms*3
  output%no_cells = output%no_atoms_sc/output%no_atoms
end function

end module
