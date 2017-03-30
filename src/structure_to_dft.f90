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

subroutine structure_to_castep(structure_sc,old_cell_filename,new_cell_filename)
  use string_module
  use structure_module
  use file_module
  implicit none
  
  type(StructureData), intent(in)           :: structure_sc
  type(String),        intent(in), optional :: old_cell_filename
  type(String),        intent(in)           :: new_cell_filename
  
  
  ! Band structure path data
  real(dp) :: kpoint(3)
  
  ! Line numbers
  integer :: lattice_start_line
  integer :: lattice_end_line
  integer :: positions_start_line
  integer :: positions_end_line
  integer :: path_start_line
  integer :: path_end_line
  
  ! Temporary variables
  integer        :: i
  type(String), allocatable :: line(:)
  
  ! Old and new cell files.
  type(String), allocatable :: old_cell_file(:)
  integer                   :: new_cell_file
  
  if (present(old_cell_filename)) then
    ! --------------------------------------------------
    ! Parse old cell file.
    ! --------------------------------------------------
    lattice_start_line = 0
    lattice_end_line = 0
    positions_start_line = 0
    positions_end_line = 0
    path_start_line = 0
    path_end_line = 0
    
    ! Get line numbers.
    old_cell_file = read_lines(old_cell_filename)
    do i=1,size(old_cell_file)
      line = split(lower_case(old_cell_file(i)))
      if (size(line) >= 2) then
        if (line(1)=='%block' .and. ( line(2)=='lattice_abc' .or. &
                                    & line(2)=='lattice_cart')) then
          lattice_start_line = i
        elseif (line(1)=='%endblock' .and. ( line(2)=='lattice_abc' .or. &
                                           & line(2)=='lattice_cart')) then
          lattice_end_line = i
        elseif (line(1)=='%block' .and. ( line(2)=='positions_abs' .or. &
                                        & line(2)=='positions_frac')) then
          positions_start_line = i
        elseif (line(1)=='%endblock' .and. ( line(2)=='positions_abs' .or. &
                                           & line(2)=='positions_frac')) then
          positions_end_line = i
        elseif (line(1)=='%block' .and. line(2)=='bs_kpoint_path') then
          path_start_line = i
        elseif ( line(1)=='%endblock' .and. line(2)=='bs_kpoint_path') then
          path_end_line = i
        endif
      endif
    enddo
    
    ! Transform k-points from fractional primitive cell co-ordinates into 
    !    fractional supercell co-ordinates.
    do i=path_start_line+1,path_end_line-1
      line = split(lower_case(old_cell_file(i)))
      kpoint = dble(line(1:3))
      kpoint = matmul(structure_sc%supercell,kpoint)
      old_cell_file(i) = kpoint//' '//join(line(4:))
    enddo
  endif
  
  ! --------------------------------------------------
  ! Write cell file.
  ! --------------------------------------------------
  new_cell_file = open_write_file(new_cell_filename)
  call print_line(new_cell_file,'%block lattice_cart')
  call print_line(new_cell_file,'bohr')
  do i=1,3
    call print_line(new_cell_file,structure_sc%lattice(i,:))
  enddo
  call print_line(new_cell_file,'%endblock lattice_cart')
  call print_line(new_cell_file,'')
  call print_line(new_cell_file,'%block positions_abs')
  call print_line(new_cell_file,'bohr')
  do i=1,structure_sc%no_atoms
    call print_line(new_cell_file, structure_sc%species(i)//' '// &
                             & structure_sc%atoms(:,i))
  enddo
  call print_line(new_cell_file,'%endblock positions_abs')
  call print_line(new_cell_file,'')
  
  ! Copy the contents of old cell file to new cell file.
  if (present(old_cell_filename)) then
    do i=1,size(old_cell_file)
      
      ! Skip lines which are changed between supercells.
      if (lattice_start_line<=i .and. i<=lattice_end_line) then
        cycle
      elseif (positions_start_line<=i .and. i<=positions_start_line) then
        cycle
      endif
      
      call print_line(new_cell_file,char(old_cell_file(i)))
    enddo
  endif
  
  close(new_cell_file)
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
  type(String)              :: previous_species
  integer                   :: no_species
  type(String), allocatable :: species(:)
  integer,      allocatable :: species_counts(:)
  
  ! Temporary variables
  integer      :: i
  type(String) :: line
  
  ! Count the number of species
  no_species = 0
  previous_species=''
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
  previous_species=''
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
  
  call print_line(poscar_file,'Structure')
  call print_line(poscar_file,bohr)
  do i=1,3
    call print_line(poscar_file, structure_sc%lattice(:,i))
  enddo
  
  line = species(1)
  do i=2,no_species
    line = line//' '//species(i)
  enddo
  call print_line(poscar_file, line)
  
  line = species_counts(1)
  do i=2,no_species
    line = line//' '//species_counts(i)
  enddo
  call print_line(poscar_file, line)
  
  call print_line(poscar_file,'Cartesian')
  do i=1,structure_sc%no_atoms
    call print_line(poscar_file, structure_sc%atoms(:,i))
  enddo
  close(poscar_file)
end subroutine

subroutine structure_to_qe(structure_sc,old_qe_in_filename,new_qe_in_filename)
  use string_module
  use structure_module
  use file_module
  use err_module
  implicit none
  
  type(StructureData), intent(in)           :: structure_sc
  type(String),        intent(in), optional :: old_qe_in_filename
  type(String),        intent(in)           :: new_qe_in_filename
  
  ! The new and old qe input files.
  type(String), allocatable :: old_qe_in_file(:)
  integer                   :: new_qe_in_file
  
  ! Temporary variables
  integer                   :: i
  type(String), allocatable :: line(:)
  
  if (present(old_qe_in_filename)) then
    old_qe_in_file = read_lines(old_qe_in_filename)
    
    ! --------------------------------------------------
    ! Transform K-points into supercell co-ordinates.
    ! --------------------------------------------------
    do i=1,size(old_qe_in_file)
      line = split(lower_case(old_qe_in_file(i)))
      if (size(line) >= 1) then
        if (line(1)=='k_points') then
          call print_line('qe K-points not yet supported.')
          call err()
        endif
      endif
    enddo
  endif
  
  ! --------------------------------------------------
  ! Write output file
  ! --------------------------------------------------
  new_qe_in_file = open_write_file(new_qe_in_filename)
  call print_line(new_qe_in_file,'nat='//structure_sc%no_atoms)
  call print_line(new_qe_in_file,'/&end')
  call print_line(new_qe_in_file,'CELL_PARAMETERS bohr')
  do i=1,3
    call print_line(new_qe_in_file, structure_sc%lattice(i,:))
  enddo
  call print_line(new_qe_in_file,'ATOMIC_POSITIONS bohr')
  do i=1,structure_sc%no_atoms
    call print_line(new_qe_in_file, structure_sc%species(i)//' '// &
                                  & structure_sc%atoms(:,i))
  enddo
  
  ! Write old qe in file contents to new qe in file.
  if (present(old_qe_in_filename)) then
    do i=1,size(old_qe_in_file)
      call print_line(new_qe_in_file,old_qe_in_file(i))
    enddo
  endif
  
  close(new_qe_in_file)
end subroutine

subroutine structure_to_dft_StructureData(dft_code,structure_sc, &
   & input_filename,output_filename)
  use string_module
  use structure_module
  use err_module
  implicit none
  
  type(String),        intent(in)           :: dft_code
  type(StructureData), intent(in)           :: structure_sc
  type(String),        intent(in), optional :: input_filename
  type(String),        intent(in)           :: output_filename
  
  if (dft_code=="castep") then
    if (present(input_filename)) then
      call structure_to_castep(structure_sc,input_filename,output_filename)
    else
      call structure_to_castep( structure_sc      = structure_sc, &
                              & new_cell_filename = output_filename)
    endif
  elseif (dft_code=="vasp") then
    if (present(input_filename)) then
      call print_line('Conversion of vasp input files not yet supported.')
      call err()
    else
      call structure_to_vasp(structure_sc,output_filename)
    endif
  elseif (dft_code=="qe") then
    if (present(input_filename)) then
      call structure_to_qe(structure_sc,input_filename,output_filename)
    else
      call structure_to_qe( structure_sc       = structure_sc, &
                          & new_qe_in_filename = output_filename)
    endif
  endif
end subroutine

subroutine structure_to_dft_filename(dft_code,structure_sc_filename, &
   & input_filename,output_filename)
  use string_module
  use structure_module
  use supercell_module
  implicit none
  
  type(String), intent(in)           :: dft_code
  type(String), intent(in)           :: structure_sc_filename
  type(String), intent(in), optional :: input_filename
  type(String), intent(in)           :: output_filename
  
  type(StructureData) :: structure_sc
  
  structure_sc = read_structure_file(structure_sc_filename)
  
  if (present(input_filename)) then
    call structure_to_dft(dft_code,structure_sc,input_filename,output_filename)
  else
    call structure_to_dft( dft_code=dft_code, &
                         & structure_sc=structure_sc, &
                         & output_filename=output_filename)
  endif
end subroutine
end module
