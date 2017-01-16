! Program to transform structure file to cell file
module structure_to_dft_module
  implicit none
  
  private
  
  public :: structure_to_dft

  interface structure_to_dft
    module procedure structure_to_dft_bs
    module procedure structure_to_dft_no_bs
  end interface

contains

subroutine structure_to_castep_no_bs(structure_filename,cell_filename)
  use string_module
  use structure_module
  use file_io
  implicit none
  
  type(String), intent(in) :: structure_filename
  type(String), intent(in) :: cell_filename
  
  ! The contents of the old cell file
  integer                   :: cell_file_length
  type(String), allocatable :: old_cell_file_contents(:)
  
  ! The contents of the structure file
  type(StructureData) :: structure
  
  ! Temporary variables
  character(100) :: line
  integer        :: i
  
  ! File units
  integer :: cell_file
  
  if (file_exists(cell_filename)) then
    cell_file = open_read_file(cell_filename)
    cell_file_length = count_lines(cell_file)
    allocate(old_cell_file_contents(cell_file_length))
    do i=1,cell_file_length
      read(cell_file,"(a)") line
      old_cell_file_contents(i) = trim(line)
    enddo
    close(cell_file)
  endif
  
  structure = read_structure_file(structure_filename)
  
  cell_file = open_write_file(cell_filename)
  write(cell_file,"(a)") '%block lattice_cart'
  write(cell_file,"(a)") 'bohr'
  do i=1,3
    write(cell_file,*) structure%lattice(i,:)
  enddo
  write(cell_file,"(a)") '%endblock lattice_cart'
  write(cell_file,"(a)") ''
  write(cell_file,"(a)") 'block positions_abs'
  write(cell_file,"(a)") 'bohr'
  do i=1,structure%no_atoms
    write(cell_file,*) structure%species(i),structure%atoms(:,i)
  enddo
  write(cell_file,"(a)") '%endblock positions_abs'
  write(cell_file,"(a)") ''
  
  if (allocated(old_cell_file_contents)) then
    do i=1,size(old_cell_file_contents)
      write(cell_file,"(a)") char(old_cell_file_contents(i))
    enddo
  endif
  
  close(cell_file)
end subroutine
 
subroutine structure_to_castep_bs(structure_filename,sc_bs_path_filename, &
   & cell_filename)
  use constants, only : dp
  use string_module
  use file_io
  implicit none
  
  type(String), intent(in) :: structure_filename
  type(String), intent(in) :: sc_bs_path_filename
  type(String), intent(in) :: cell_filename
  
  ! file units
  integer :: sc_bs_path_file
  integer :: cell_file
  
  ! bs_path data
  integer               :: no_points
  real(dp), allocatable :: sc_bs_path(:,:)
  
  ! temporary variables
  integer :: i
  
  sc_bs_path_file = open_read_file(sc_bs_path_filename)
  no_points = count_lines(sc_bs_path_file)
  allocate(sc_bs_path(no_points,3))
  do i=1,no_points
    read(sc_bs_path_file,*) sc_bs_path(i,:)
  enddo
  close(sc_bs_path_file)
  
  call structure_to_castep_no_bs(structure_filename,cell_filename)
  
  cell_file = open_append_file(cell_filename)
  write(cell_file,"(a)") ''
  write(cell_file,"(a)") '%block_bs_kpoints_path'
  do i=1,no_points
    write(cell_file,*) sc_bs_path(i,:)
  enddo
  write(cell_file,"(a)") '%endblock_bs_kpoint_path'
  write(cell_file,"(a)") ''
  write(cell_file,"(a)") 'bs_kpoints_path_spacing = 10.0'
end subroutine

subroutine structure_to_vasp(structure_filename,poscar_filename)
  use string_module
  use structure_module
  use constants, only : bohr
  use file_io
  implicit none
  
  type(String), intent(in) :: structure_filename
  type(String), intent(in) :: poscar_filename
  
  ! Structure file data
  type(StructureData) :: structure
  
  ! file units
  integer :: poscar_file
  
  ! species data
  character(2)              :: previous_species
  integer                   :: no_species
  character(2), allocatable :: species(:)
  integer,      allocatable :: species_counts(:)
  
  ! Temporary variables
  integer :: i
  
  structure = read_structure_file(structure_filename)
  
  ! count the number of species
  no_species = 0
  previous_species='XX'
  do i=1,structure%no_atoms
    if (structure%species(i)/=previous_species) then
      previous_species = structure%species(i)
      no_species = no_species+1
    endif
  enddo
  
  ! generate species lists
  allocate(species(no_species))
  allocate(species_counts(no_species))
  no_species = 0
  previous_species='XX'
  species_counts = 0
  do i=1,structure%no_atoms
    if (structure%species(i)/=previous_species) then
      previous_species = structure%species(i)
      no_species = no_species+1
      species(no_species) = structure%species(i)
    endif
    species_counts(no_species) = species_counts(no_species)+1
  enddo
  
  poscar_file = open_write_file(poscar_filename)
  write(poscar_file,"(a)") 'Structure'
  write(poscar_file,*) bohr
  do i=1,3
    write(poscar_file,*) structure%lattice(:,i)
  enddo
  do i=1,no_species
    write(poscar_file,"(a)",advance="no") species(i)//" "
  enddo
  write(poscar_file,*)
  do i=1,no_species
    write(poscar_file,"(a)",advance="no") char(str(species_counts(i))//" ")
  enddo
  write(poscar_file,*)
  write(poscar_file,"(a)") 'Cartesian'
  do i=1,structure%no_atoms
    write(poscar_file,*) structure%atoms(:,i)
  enddo
  close(poscar_file)
end subroutine

subroutine structure_to_dft_no_bs(codename,structure_filename,output_filename)
  use string_module
  implicit none
  
  type(String), intent(in) :: codename
  type(String), intent(in) :: structure_filename
  type(String), intent(in) :: output_filename
  
  if (codename=="castep") then
    call structure_to_castep_no_bs(structure_filename,output_filename)
  elseif (codename=="vasp") then
    call structure_to_vasp(structure_filename,output_filename)
  endif
end subroutine

subroutine structure_to_dft_bs(codename,structure_filename, &
   & sc_bs_path_filename,cell_filename)
  use string_module
  implicit none
  
  type(String), intent(in) :: codename
  type(String), intent(in) :: structure_filename
  type(String), intent(in) :: sc_bs_path_filename
  type(String), intent(in) :: cell_filename
  
  if (codename=="castep") then
    call structure_to_castep_bs(structure_filename,sc_bs_path_filename, &
      & cell_filename)
  endif
end subroutine
end module
