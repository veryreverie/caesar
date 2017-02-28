! Program to transform structure file to cell file
module structure_to_dft_module
  implicit none
  
  private
  
  public :: structure_to_dft

  interface structure_to_dft
    ! For use with StructureData type
    module procedure structure_to_dft_StructureData
    ! For use with structure file
    module procedure structure_to_dft_filename
  end interface

contains

subroutine structure_to_castep(structure_sc,input_filename, &
   & path_filename,structure,output_filename)
  use string_module
  use structure_module
  use file_module
  implicit none
  
  type(StructureData), intent(in)           :: structure_sc
  type(String),        intent(in), optional :: input_filename
  type(String),        intent(in), optional :: path_filename
  type(StructureData), intent(in), optional :: structure
  type(String),        intent(in)           :: output_filename
  
  ! The contents of the old cell file
  type(String), allocatable :: input_file_contents(:)
  
  ! Band structure path data
  integer               :: no_points
  real(dp), allocatable :: path(:,:)
  
  ! Line numbers
  integer :: path_start_line
  integer :: path_end_line
  
  ! Temporary variables
  integer        :: i
  
  ! File contents
  type(String), allocatable :: path_file(:)
  type(String), allocatable :: line(:)
  
  ! File units
  integer :: cell_file
  
  ! Whether or not path exists
  logical :: path_wanted
  
  path_wanted = .false.
  if (present(structure) .and. present(path_filename)) then
    if (file_exists(path_filename)) then
      path_wanted = .true.
    endif
  endif
  
  if (path_wanted) then
    ! Parse path file
    path_file = read_lines(path_filename)
    do i=1,size(path_file)
      line = split(lower_case(path_file(i)))
      if (line(1)=="block") then
        path_start_line = i
      elseif (line(1)=="endblock") then
        path_end_line = i
      endif
    enddo
    
    no_points = path_end_line-path_start_line-1
    allocate(path(3,no_points))
    
    do i=1,path_end_line-path_start_line-1
      line = split(path_file(path_start_line+i))
      path(:,i) = dble(line)
    enddo
    
    path = matmul(transpose(structure%recip_lattice),path)
    path = matmul(structure_sc%lattice,path)
  endif
  
  ! Write cell file
  cell_file = open_write_file(output_filename)
  write(cell_file,"(a)") '%block lattice_cart'
  write(cell_file,"(a)") 'bohr'
  do i=1,3
    write(cell_file,*) structure_sc%lattice(i,:)
  enddo
  write(cell_file,"(a)") '%endblock lattice_cart'
  write(cell_file,"(a)") ''
  write(cell_file,"(a)") '%block positions_abs'
  write(cell_file,"(a)") 'bohr'
  do i=1,structure_sc%no_atoms
    write(cell_file,*) structure_sc%species(i),structure_sc%atoms(:,i)
  enddo
  write(cell_file,"(a)") '%endblock positions_abs'
  write(cell_file,"(a)") ''
  
  ! Copy the contents of input file to cell file
  if (present(input_filename)) then
    if (file_exists(input_filename)) then
      input_file_contents = read_lines(input_filename)
      do i=1,size(input_file_contents)
        write(cell_file,"(a)") char(input_file_contents(i))
      enddo
    endif
  endif
  
  if (path_wanted) then
    ! Append band structure data to cell file
    write(cell_file,"(a)") ''
    write(cell_file,"(a)") '%block_bs_kpoints_path'
    do i=1,no_points
      write(cell_file,*) path(:,i)
    enddo
    write(cell_file,"(a)") '%endblock_bs_kpoint_path'
    write(cell_file,"(a)") ''
    write(cell_file,"(a)") 'bs_kpoints_path_spacing = 10.0'
  endif
  
  close(cell_file)
end subroutine

subroutine structure_to_vasp(structure_sc,poscar_filename)
  use string_module
  use structure_module
  use constants, only : bohr
  use file_module
  implicit none
  
  type(StructureData), intent(in) :: structure_sc
  type(String),        intent(in) :: poscar_filename
  
  ! File units
  integer :: poscar_file
  
  ! Species data
  character(2)              :: previous_species
  integer                   :: no_species
  character(2), allocatable :: species(:)
  integer,      allocatable :: species_counts(:)
  
  ! Temporary variables
  integer :: i
  
  ! Count the number of species
  no_species = 0
  previous_species='XX'
  do i=1,structure_sc%no_atoms
    if (structure_sc%species(i)/=previous_species) then
      previous_species = structure_sc%species(i)
      no_species = no_species+1
    endif
  enddo
  
  ! Generate species lists
  allocate(species(no_species))
  allocate(species_counts(no_species))
  no_species = 0
  previous_species='XX'
  species_counts = 0
  do i=1,structure_sc%no_atoms
    if (structure_sc%species(i)/=previous_species) then
      previous_species = structure_sc%species(i)
      no_species = no_species+1
      species(no_species) = structure_sc%species(i)
    endif
    species_counts(no_species) = species_counts(no_species)+1
  enddo
  
  ! Write output file
  poscar_file = open_write_file(poscar_filename)
  write(poscar_file,"(a)") 'Structure'
  write(poscar_file,*) bohr
  do i=1,3
    write(poscar_file,*) structure_sc%lattice(:,i)
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
  do i=1,structure_sc%no_atoms
    write(poscar_file,*) structure_sc%atoms(:,i)
  enddo
  close(poscar_file)
end subroutine

subroutine structure_to_qe(structure_sc,input_filename,pseudo_filename, &
   & kpoints_filename,structure,output_filename)
  use string_module
  use structure_module
  use file_module
  implicit none
  
  type(StructureData), intent(in) :: structure_sc
  type(String),        intent(in) :: input_filename
  type(String),        intent(in) :: pseudo_filename
  type(String),        intent(in) :: kpoints_filename
  type(StructureData), intent(in) :: structure
  type(String),        intent(in) :: output_filename
  
  ! File contents
  type(String), allocatable :: pseudo_contents(:)
  
  ! The contents of the old .in file
  type(String), allocatable :: input_file_contents(:)
  
  ! kpoints.in data
  integer        :: primitive_mesh(3)
  real(dp)       :: sc_distance(3)
  real(dp)       :: distance(3)
  
  ! File contents
  type(String), allocatable :: kpoints_file(:)
  
  ! File units
  integer :: output_file
  
  ! Temporary variables
  integer        :: i
  
  ! Read in pseudo file
  pseudo_contents = read_lines(pseudo_filename)
  
  ! Read in kpoints file
  kpoints_file = read_lines(kpoints_filename)
  primitive_mesh = int(split(kpoints_file(2)))
  
  ! Construct reciprocal primitive lattice
  do i=1,3
    distance(i) = norm2(structure%recip_lattice(i,:))
  enddo
  
  ! Construct reciprocal supercell lattice
  do i=1,3
    sc_distance(i) = norm2(structure_sc%recip_lattice(i,:))
  enddo
  
  output_file = open_write_file(output_filename)
  ! Write input file to output file
  if (file_exists(input_filename)) then
    input_file_contents = read_lines(input_filename)
    do i=1,size(input_file_contents)
      write(output_file,"(a)") char(input_file_contents(i))
    enddo
  endif
  
  ! Write output file
  write(output_file,"(a)") char(str('nat=')//structure_sc%no_atoms)
  write(output_file,"(a)") '/&end'
  do i=1,size(pseudo_contents)
    write(output_file,"(a)") char(pseudo_contents(i))
  enddo
  write(output_file,"(a)") 'CELL_PARAMETERS bohr'
  do i=1,3
    write(output_file,"(a)") structure_sc%lattice(i,:)
  enddo
  write(output_file,"(a)") 'ATOMIC_POSITIONS bohr'
  do i=1,structure_sc%no_atoms
    write(output_file,*) structure_sc%species(i),structure_sc%atoms(:,i)
  enddo
  write(output_file,"(a)") char(kpoints_file(1))
  write(output_file,*) int(primitive_mesh*sc_distance/distance)+1,0,0,0
  close(output_file)
end subroutine

subroutine structure_to_dft_StructureData(dft_code,structure_sc,      &
   & input_filename,path_filename,pseudo_filename, &
   & kpoints_filename,structure,output_filename)
  use string_module
  use structure_module
  implicit none
  
  type(String),        intent(in)           :: dft_code
  type(StructureData), intent(in)           :: structure_sc
  type(String),        intent(in), optional :: input_filename     ! castep & qe
  type(String),        intent(in), optional :: path_filename      ! castep only
  type(String),        intent(in), optional :: pseudo_filename    ! qe only
  type(String),        intent(in), optional :: kpoints_filename   ! qe only
  type(StructureData), intent(in), optional :: structure          ! castep & qe
  type(String),        intent(in)           :: output_filename
  
  if (dft_code=="castep") then
    call structure_to_castep(structure_sc,input_filename, &
       & path_filename,structure,output_filename)
  elseif (dft_code=="vasp") then
    call structure_to_vasp(structure_sc,output_filename)
  elseif (dft_code=="qe") then
    call structure_to_qe(structure_sc,input_filename,pseudo_filename, &
       & kpoints_filename,structure,output_filename)
  endif
end subroutine

subroutine structure_to_dft_filename(dft_code,structure_sc_filename,  &
   & input_filename,path_filename,pseudo_filename, &
   & kpoints_filename,structure_filename,output_filename)
  use string_module
  use structure_module
  use supercell_module
  implicit none
  
  type(String), intent(in)           :: dft_code
  type(String), intent(in)           :: structure_sc_filename
  type(String), intent(in), optional :: input_filename     ! castep and qe only
  type(String), intent(in), optional :: path_filename      ! castep only
  type(String), intent(in), optional :: pseudo_filename    ! qe only
  type(String), intent(in), optional :: kpoints_filename   ! qe only
  type(String), intent(in), optional :: structure_filename ! castep and qe only
  type(String), intent(in)           :: output_filename
  
  type(StructureData) :: structure_sc
  type(StructureData) :: structure
  
  ! n.b. supercell is not used here, so a dummy is provided.
  structure_sc = read_structure_file( structure_sc_filename, &
                                    & identity_supercell())
  
  if (present(structure_filename)) then
    ! n.b. supercell is not used here, so a dummy is provided.
    structure = read_structure_file(structure_filename, identity_supercell())
    call structure_to_dft( dft_code=dft_code, &
                         & structure_sc=structure_sc, &
                         & input_filename=input_filename, &
                         & path_filename=path_filename, &
                         & pseudo_filename=pseudo_filename, &
                         & kpoints_filename=kpoints_filename, &
                         & structure=structure, &
                         & output_filename=output_filename)
  else
    call structure_to_dft( dft_code=dft_code, &
                         & structure_sc=structure_sc, &
                         & input_filename=input_filename, &
                         & path_filename=path_filename, &
                         & pseudo_filename=pseudo_filename, &
                         & kpoints_filename=kpoints_filename, &
                         & output_filename=output_filename)
  endif
end subroutine
end module
