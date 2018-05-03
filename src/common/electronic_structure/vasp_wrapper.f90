! ======================================================================
! Reads and writes VASP .poscar files.
! ======================================================================

! N.B. This module has never been tested.

module vasp_wrapper_submodule
  use utils_module
  
  use structure_module
  
  use electronic_structure_data_submodule
  implicit none
  
  private
  
  public :: make_input_filename_vasp
  public :: make_output_filename_vasp
  public :: read_input_file_vasp
  public :: write_input_file_vasp
  public :: read_output_file_vasp
contains

function make_input_filename_vasp(seedname) result(output)
  implicit none
  
  type(String), intent(in) :: seedname
  type(String)             :: output
  
  output = seedname//'.poscar'
end function

function make_output_filename_vasp(seedname) result(output)
  implicit none
  
  type(String), intent(in) :: seedname
  type(String)             :: output
  
  output = 'OUTCAR'
end function

function read_input_file_vasp(filename,symmetry_precision, &
   & calculate_symmetry) result(output)
  implicit none
  
  type(String), intent(in)           :: filename
  real(dp),     intent(in)           :: symmetry_precision
  logical,      intent(in), optional :: calculate_symmetry
  type(StructureData)                :: output
  
  call print_line(CODE_ERROR//': VASP not yet supported.')
  call err()
end function

subroutine write_input_file_vasp(structure,poscar_filename)
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(String),        intent(in) :: poscar_filename
  
  ! File units
  type(OFile) :: poscar_file
  
  ! Species data
  type(String)              :: previous_species
  integer                   :: no_species
  type(String), allocatable :: species(:)
  integer,      allocatable :: species_counts(:)
  
  ! Temporary variables
  integer      :: i, ialloc
  type(String) :: line
  
  ! Count the number of species
  no_species = 0
  previous_species=''
  do i=1,structure%no_atoms
    if (structure%atoms(i)%species()/=previous_species) then
      previous_species = structure%atoms(i)%species()
      no_species = no_species+1
    endif
  enddo
  
  ! Generate species lists
  allocate( species(no_species),        &
          & species_counts(no_species), &
          & stat=ialloc); call err(ialloc)
  no_species = 0
  previous_species=''
  species_counts = 0
  do i=1,structure%no_atoms
    if (structure%atoms(i)%species()/=previous_species) then
      previous_species = structure%atoms(i)%species()
      no_species = no_species+1
      species(no_species) = structure%atoms(i)%species()
    endif
    species_counts(no_species) = species_counts(no_species)+1
  enddo
  
  ! Write output file
  poscar_file = OFile(poscar_filename)
  
  call poscar_file%print_line('Structure')
  call poscar_file%print_line(ANGSTROM_PER_BOHR)
  call poscar_file%print_lines(structure%lattice)
  
  line = species(1)
  do i=2,no_species
    line = line//' '//species(i)
  enddo
  call poscar_file%print_line(line)
  
  line = species_counts(1)
  do i=2,no_species
    line = line//' '//species_counts(i)
  enddo
  call poscar_file%print_line(line)
  
  call poscar_file%print_line('Cartesian')
  do i=1,structure%no_atoms
    call poscar_file%print_line(structure%atoms(i)%cartesian_position())
  enddo
end subroutine

function read_output_file_vasp(filename,structure) result(output)
  implicit none
  
  type(String),        intent(in) :: filename
  type(StructureData), intent(in) :: structure
  type(ElectronicStructure)       :: output
  
  call print_line(CODE_ERROR//': VASP not yet supported.')
  call err()
end function
end module
