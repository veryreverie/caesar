! ======================================================================
! Converts DFT input files to and from StructureData.
! ======================================================================
module dft_input_file_module
  use constants_module, only : dp
  use string_module
  use io_module
  use linear_algebra_module
  
  private
  
  public :: make_dft_input_filename
  public :: dft_input_file_to_StructureData
  public :: StructureData_to_dft_input_file
  public :: reduce_dft_input_file ! For testing purposes only.
  
  type CastepInputFile
    type(String), allocatable :: lattice_block(:)
    type(String), allocatable :: positions_block(:)
    type(String), allocatable :: masses_block(:)
    type(String), allocatable :: qpoints_block(:)
    type(String), allocatable :: remainder(:)
  end type
contains

! ----------------------------------------------------------------------
! Converts a file seedname into the appropriate dft input filename.
! ----------------------------------------------------------------------
function make_dft_input_filename(dft_code,seedname) result(output)
  implicit none
  
  type(String), intent(in) :: dft_code
  type(String), intent(in) :: seedname
  type(String)             :: output
  
  if (dft_code == 'castep') then
    output = seedname//'.cell'
  elseif (dft_code == 'vasp') then
    output = 'POSCAR'
  elseif (dft_code == 'qe') then
    output = seedname//'.in'
  else
    call print_line('Unrecognised dft code: '//dft_code)
    call err()
  endif
end function


! ----------------------------------------------------------------------
! Helper functions to parse DFT input files.
! ----------------------------------------------------------------------

! Splits a .cell file into blocks.
function parse_castep_input_file(filename) result(output)
  use ifile_module
  implicit none
  
  type(String), intent(in) :: filename
  type(CastepInputFile)    :: output
  
  ! File contents.
  type(IFile) :: cell_file
  
  ! Line numbers.
  integer :: lattice_block_start
  integer :: lattice_block_size
  integer :: positions_block_start
  integer :: positions_block_size
  integer :: masses_block_start
  integer :: masses_block_size
  integer :: qpoints_block_start
  integer :: qpoints_block_size
  integer :: remainder_size
  
  ! Temporary variables.
  type(String), allocatable :: line(:)
  integer                   :: i,j,k,ialloc
  
  cell_file = filename
  
  lattice_block_start = 0
  lattice_block_size = 0
  positions_block_start = 0
  positions_block_size = 0
  masses_block_start = 0
  masses_block_size = 0
  qpoints_block_start = 0
  qpoints_block_size = 0
  remainder_size = 0
  
  ! Work out line numbers.
  do i=1,size(cell_file)
    line = split(lower_case(cell_file%line(i)))
    if (size(line) >= 2) then
      if (line(1)=='%block' .and. ( line(2)=='lattice_cart' .or. &
                                  & line(2)=='lattice_abc')) then
        lattice_block_start = i
      elseif (line(1)=='%endblock' .and. ( line(2)=='lattice_cart' .or. &
                                         & line(2)=='lattice_abc')) then
        lattice_block_size = i-lattice_block_start+1
      elseif (line(1)=='%block' .and. ( line(2)=='positions_abs' .or. &
                                      & line(2)=='positions_frac')) then
        positions_block_start = i
      elseif (line(1)=='%endblock' .and. ( line(2)=='positions_abs' .or. &
                                         & line(2)=='positions_frac')) then
        positions_block_size = i-positions_block_start+1
      elseif (line(1)=='%block' .and. line(2)=='species_mass') then
        masses_block_start = i
      elseif (line(1)=='%endblock' .and. line(2)=='species_mass') then
        masses_block_size = i-masses_block_start+1
      elseif (line(1)=='%block' .and. line(2)=='bs_qpoint_path') then
        qpoints_block_start = i
      elseif (line(1)=='%endblock' .and. line(2)=='bs_qpoint_path') then
        qpoints_block_size = i-qpoints_block_start+1
      endif
    endif
  enddo
  
  remainder_size = size(cell_file)      &
               & - lattice_block_size   &
               & - positions_block_size &
               & - masses_block_size    &
               & - qpoints_block_size
  
  ! Check required blocks exist and all blocks are of reasonable sizes.
  if (lattice_block_start==0) then
    call print_line('Error: lattice block not present in '//filename)
    call err()
  elseif (lattice_block_size <= 2) then
    call print_line('Error: lattice block of unexpected size in '//filename)
    call err()
  elseif (positions_block_start==0) then
    call print_line('Error: positions block not present in '//filename)
    call err()
  elseif (positions_block_size <= 2) then
    call print_line('Error: positions block of unexpected size in '//filename)
    call err()
  elseif (masses_block_start==0) then
    call print_line('Error: species_mass block not present in '//filename)
    call err()
  elseif (masses_block_size <= 2) then
    call print_line('Error: species_mass block of unexpected size in '// &
       & filename)
    call err()
  elseif (qpoints_block_start/=0 .and. qpoints_block_size<=2) then
    call print_line('Error: bs_qpoint_path block of unexpected size in '// &
       & filename)
    call err()
  elseif (qpoints_block_start==0 .and. qpoints_block_size/=0) then
    call print_line('Error: bs_qpoint_path block of unexpected size in '// &
       & filename)
    call err()
  endif
  
  ! Allocate space for the blocks.
  allocate( output%lattice_block(lattice_block_size),     &
          & output%positions_block(positions_block_size), &
          & output%masses_block(masses_block_size),       &
          & output%qpoints_block(qpoints_block_size),     &
          & output%remainder(remainder_size),             &
          & stat=ialloc); call err(ialloc)
  
  ! Copy across file contents.
  k = 0
  do i=1,size(cell_file)
    
    j = i-lattice_block_start+1
    if (j>0 .and. j<=lattice_block_size) then
      output%lattice_block(j) = cell_file%line(i)
      cycle
    endif
    
    j = i-positions_block_start+1
    if (j>0 .and. j<=positions_block_size) then
      output%positions_block(j) = cell_file%line(i)
      cycle
    endif
    
    j = i-masses_block_start+1
    if (j>0 .and. j<=masses_block_size) then
      output%masses_block(j) = cell_file%line(i)
      cycle
    endif
    
    j = i-qpoints_block_start+1
    if (j>0 .and. j<qpoints_block_size) then
      output%qpoints_block(j) = cell_file%line(i)
      cycle
    endif
    
    k = k+1
    output%remainder(k) = cell_file%line(i)
  enddo
end function

! ----------------------------------------------------------------------
! Reading DFT input files to StructureData.
! ----------------------------------------------------------------------
function castep_input_file_to_StructureData(filename, symmetry_precision) &
   & result(output)
  use constants_module, only : pi, angstrom_per_bohr, kg_per_me, kg_per_amu
  use structure_module
  use linear_algebra_module
  implicit none
  
  type(String), intent(in) :: filename
  real(dp),     intent(in) :: symmetry_precision
  type(StructureData)      :: output
  
  type(CastepInputFile) :: cell_file
  
  ! Lattice variables.
  logical  :: lattice_is_cart
  real(dp) :: conversion
  real(dp) :: lattice(3,3)
  real(dp) :: lengths(3)
  real(dp) :: angles(3)
  
  ! Atomic variables.
  logical                   :: positions_are_abs
  type(String), allocatable :: species(:)
  type(RealVector), allocatable :: positions(:)
  integer                   :: no_atoms
  logical,      allocatable :: masses_found(:)
  
  ! Temporary variables.
  type(String), allocatable :: line(:)
  integer                   :: i,j,ialloc
  
  cell_file = parse_castep_input_file(filename)
  
  ! Read data format.
  line = split(lower_case(cell_file%lattice_block(1)))
  lattice_is_cart = line(2)=='lattice_cart'
  
  line = split(lower_case(cell_file%positions_block(1)))
  positions_are_abs = line(2)=='positions_abs'
  
  ! Parse lattice.
  conversion = 1.0_dp / angstrom_per_bohr
  j=0
  do i=2,size(cell_file%lattice_block)-1
    line = split(lower_case(cell_file%lattice_block(i)))
    
    ! Ignore comments.
    if (slice(line(1),1,1)=='!') then
      cycle
    
    ! Read units if present.
    elseif (line(1) == 'bohr') then
      conversion = 1.0_dp
    elseif (line(1) == 'a0') then
      conversion = 1.0_dp
    elseif (line(1) == 'm') then
      conversion = 1e10_dp / angstrom_per_bohr
    elseif (line(1) == 'cm') then
      conversion = 1e8_dp / angstrom_per_bohr
    elseif (line(1) == 'nm') then
      conversion = 1e1_dp / angstrom_per_bohr
    elseif (line(1) == 'ang') then
      conversion = 1.0_dp / angstrom_per_bohr
      
    elseif (lattice_is_cart) then
      j = j+1
      if (j>3) then
        call print_line('Error: too many lattice lines found in '//filename)
        call err()
      endif
      lattice(j,:) = dble(line(1:3))*conversion
    
    else
      j = j+1
      if (j==1) then
        lengths = dble(line(1:3)) * conversion
      elseif (j==2) then
        angles = dble(line(1:3)) * pi/180
      else
        call print_line('Error: too many lattice lines found in '//filename)
        call err()
      endif
    endif
  enddo
  
  if (.not. lattice_is_cart) then
    lattice = 0.0_dp
    lattice(1,1) = lengths(1)
    lattice(2,1) = lengths(2)*cos(angles(3))
    lattice(2,2) = lengths(2)*sin(angles(3))
    lattice(3,1) = lengths(3)*cos(angles(2))
    lattice(3,2) = lengths(3)* ( cos(angles(1))                  &
                           &   - cos(angles(2))*cos(angles(3)) ) &
                           & / sin(angles(3))
    lattice(3,3) = sqrt(lengths(3)**2 - lattice(3,1)**2 - lattice(3,2)**2)
  endif
  
  ! Parse atomic positions.
  conversion = angstrom_per_bohr
  j=0
  allocate( species(size(cell_file%positions_block)-2),     &
          & positions(size(cell_file%positions_block)-2), &
          & stat=ialloc); call err(ialloc)
  
  do i=2,size(cell_file%positions_block)-1
    line = split(lower_case(cell_file%positions_block(i)))
    
    ! Ignore comments.
    if (slice(line(1),1,1)=='!') then
      cycle
    
    ! Read units if present.
    elseif (line(1) == 'bohr' .or. line(1) == 'a0') then
      conversion = 1.0_dp
    elseif (line(1) == 'm') then
      conversion = 1e10_dp / angstrom_per_bohr
    elseif (line(1) == 'cm') then
      conversion = 1e8_dp / angstrom_per_bohr
    elseif (line(1) == 'nm') then
      conversion = 1e1_dp / angstrom_per_bohr
    elseif (line(1) == 'ang') then
      conversion = 1.0_dp / angstrom_per_bohr
    
    elseif (j>size(positions)) then
      call print_line('Error: too many atom lines found in '//filename)
    
    ! Read in atomic positions.
    else
      j = j+1
      line = split(cell_file%positions_block(i)) ! N.B. no lower_case
      species(j) = line(1)
      positions(j) = dble(line(2:4))
    endif
  enddo
  
  no_atoms = j
  
  if (positions_are_abs) then
    do i=1,no_atoms
      positions(i) = positions(i) * conversion
    enddo
  else
    do i=1,no_atoms
      positions(i) = transpose(mat(lattice)) * positions(i)
    enddo
  endif
  
  ! Make output.
  output = StructureData(no_atoms, 0, 1)
  
  output%lattice = lattice
  do i=1,no_atoms
    call output%atoms(i)%set_species(species(i))
    call output%atoms(i)%set_cartesian_position(positions(i))
  enddo
  
  ! Parse masses.
  conversion = kg_per_amu / kg_per_me
  allocate(masses_found(no_atoms), stat=ialloc); call err(ialloc)
  masses_found = .false.
  do i=2,size(cell_file%masses_block)-1
    line = split(lower_case(cell_file%masses_block(i)))
    
    ! Ignore comments.
    if (slice(line(1),1,1)=='!') then
      cycle
    
    ! Read units if present.
    elseif (line(1) == 'me') then
      conversion = 1.0_dp
    elseif (line(1) == 'amu') then
      conversion = kg_per_amu / kg_per_me
    elseif (line(1) == 'kg') then
      conversion = 1.0_dp / kg_per_me
    elseif (line(1) == 'g') then
      conversion = 1e-3_dp / kg_per_me
    
    ! Read in masses.
    else
      line = split(cell_file%masses_block(i)) ! N.B. no lower_case
      do j=1,no_atoms
        if (line(1)==species(j)) then
          masses_found(j) = .true.
          call output%atoms(j)%set_mass(dble(line(2))*conversion)
        endif
      enddo
    endif
  enddo
  
  if (.not. all(masses_found)) then
    call print_line('Error: not all masses specified in '//filename)
    call err()
  endif
  
  output%supercell = make_identity_matrix(3)
  output%rvectors(1) = [0,0,0]
  output%gvectors(1) = [0,0,0]
  
  call output%calculate_derived_quantities()
  call output%calculate_symmetry(symmetry_precision)
end function

function dft_input_file_to_StructureData(dft_code,filename, &
   & symmetry_precision) result(output)
  use structure_module
  implicit none
  
  type(String), intent(in) :: dft_code
  type(String), intent(in) :: filename
  real(dp),     intent(in) :: symmetry_precision
  type(StructureData)      :: output
  
  if (dft_code == 'castep') then
    output = castep_input_file_to_StructureData(filename, symmetry_precision)
  else
    call print_line('Reading '//dft_code//' input file not yet supported.')
    call err()
  endif
end function

! ----------------------------------------------------------------------
! Writing StructureData to DFT input file.
! ----------------------------------------------------------------------
subroutine StructureData_to_castep_input_file(structure,old_cell_filename, &
   & new_cell_filename)
  use ofile_module
  use structure_module
  implicit none
  
  type(StructureData), intent(in)           :: structure
  type(String),        intent(in), optional :: old_cell_filename
  type(String),        intent(in)           :: new_cell_filename
  
  ! Band structure path data.
  type(RealVector) :: qpoint
  
  ! Old and new cell files.
  type(CastepInputFile) :: old_cell_file
  type(OFile)           :: new_cell_file
  
  ! Temporary variables.
  integer        :: i
  type(String), allocatable :: line(:)
  
  ! --------------------------------------------------
  ! Parse old cell file.
  ! --------------------------------------------------
  if (present(old_cell_filename)) then
    old_cell_file = parse_castep_input_file(old_cell_filename)
    
    ! Transform q-points from fractional primitive cell co-ordinates into 
    !    fractional supercell co-ordinates.
    do i=2,size(old_cell_file%qpoints_block)-1
      line = split(old_cell_file%qpoints_block(i))
      qpoint = dble(line(1:3))
      qpoint = transpose(structure%supercell) * qpoint
      old_cell_file%qpoints_block(i) = qpoint//' '//join(line(4:))
    enddo
  endif
  
  ! --------------------------------------------------
  ! Write new cell file.
  ! --------------------------------------------------
  new_cell_file = new_cell_filename
  call new_cell_file%print_line('%block lattice_cart')
  call new_cell_file%print_line('bohr')
  call new_cell_file%print_line(structure%lattice)
  call new_cell_file%print_line('%endblock lattice_cart')
  call new_cell_file%print_line('')
  call new_cell_file%print_line('%block positions_abs')
  call new_cell_file%print_line('bohr')
  do i=1,structure%no_atoms
    call new_cell_file%print_line( structure%atoms(i)%species() //' '// &
                                 & structure%atoms(i)%cartesian_position())
  enddo
  call new_cell_file%print_line('%endblock positions_abs')
  call new_cell_file%print_line('')
  
  ! Copy the contents of old cell file to new cell file.
  if (present(old_cell_filename)) then
    do i=1,size(old_cell_file%qpoints_block)
      call new_cell_file%print_line(old_cell_file%qpoints_block(i))
    enddo
    
    do i=1,size(old_cell_file%remainder)
      call new_cell_file%print_line(old_cell_file%remainder(i))
    enddo
  endif
end subroutine

subroutine StructureData_to_vasp_input_file(structure,poscar_filename)
  use constants_module, only : angstrom_per_bohr
  use ofile_module
  use structure_module
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
  poscar_file = poscar_filename
  
  call poscar_file%print_line('Structure')
  call poscar_file%print_line(angstrom_per_bohr)
  call poscar_file%print_line(structure%lattice)
  
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

subroutine StructureData_to_qe_input_file(structure,old_qe_in_filename, &
   & new_qe_in_filename)
  use ifile_module
  use ofile_module
  use structure_module
  implicit none
  
  type(StructureData), intent(in)           :: structure
  type(String),        intent(in), optional :: old_qe_in_filename
  type(String),        intent(in)           :: new_qe_in_filename
  
  ! The new and old qe input files.
  type(IFile) :: old_qe_in_file
  type(OFile) :: new_qe_in_file
  
  ! Temporary variables
  integer                   :: i
  type(String), allocatable :: line(:)
  
  if (present(old_qe_in_filename)) then
    old_qe_in_file = old_qe_in_filename
    
    ! --------------------------------------------------
    ! Transform q-points into supercell co-ordinates.
    ! --------------------------------------------------
    do i=1,size(old_qe_in_file)
      line = split(lower_case(old_qe_in_file%line(i)))
      if (size(line) >= 1) then
        if (line(1)=='k_points') then
          call print_line('qe q-points not yet supported.')
          call err()
        endif
      endif
    enddo
  endif
  
  ! --------------------------------------------------
  ! Write output file
  ! --------------------------------------------------
  new_qe_in_file = new_qe_in_filename
  call new_qe_in_file%print_line('nat='//structure%no_atoms)
  call new_qe_in_file%print_line('/&end')
  call new_qe_in_file%print_line('CELL_PARAMETERS bohr')
  call new_qe_in_file%print_line(structure%lattice)
  call new_qe_in_file%print_line('ATOMIC_POSITIONS bohr')
  do i=1,structure%no_atoms
    call new_qe_in_file%print_line( structure%atoms(i)%species() //' '// &
                                  & structure%atoms(i)%cartesian_position())
  enddo
  
  ! Write old qe in file contents to new qe in file.
  if (present(old_qe_in_filename)) then
    do i=1,size(old_qe_in_file)
      call new_qe_in_file%print_line(old_qe_in_file%line(i))
    enddo
  endif
end subroutine

subroutine StructureData_to_dft_input_file(dft_code,structure,input_filename, &
   & output_filename)
  use structure_module
  implicit none
  
  type(String),        intent(in)           :: dft_code
  type(StructureData), intent(in)           :: structure
  type(String),        intent(in), optional :: input_filename
  type(String),        intent(in)           :: output_filename
  
  if (dft_code=="castep") then
    if (present(input_filename)) then
      call StructureData_to_castep_input_file(structure,input_filename, &
         & output_filename)
    else
      call StructureData_to_castep_input_file(structure      = structure, &
                                             &new_cell_filename = output_filename)
    endif
  elseif (dft_code=="vasp") then
    if (present(input_filename)) then
      call print_line('Conversion of vasp input files not yet supported.')
      call err()
    else
      call StructureData_to_vasp_input_file(structure,output_filename)
    endif
  elseif (dft_code=="qe") then
    if (present(input_filename)) then
      call StructureData_to_qe_input_file(structure,input_filename, &
         & output_filename)
    else
      call StructureData_to_qe_input_file( structure       = structure, &
                                         & new_qe_in_filename = output_filename)
    endif
  endif
end subroutine

! ----------------------------------------------------------------------
! Converts a DFT input file into the form accepted by the old Caesar code.
! ----------------------------------------------------------------------
subroutine reduce_castep_input_file(input_filename,output_filename)
  use ofile_module
  implicit none
  
  type(String), intent(in) :: input_filename
  type(String), intent(in) :: output_filename
  
  type(CastepInputFile) :: cell_file_in
  type(OFile)           :: cell_file_out
  
  integer :: i
  
  cell_file_in = parse_castep_input_file(input_filename)
  
  cell_file_out = output_filename
  
  do i=1,size(cell_file_in%qpoints_block)
    call cell_file_out%print_line(cell_file_in%qpoints_block(i))
  enddo
  
  do i=1,size(cell_file_in%remainder)
    call cell_file_out%print_line(cell_file_in%remainder(i))
  enddo
end subroutine

subroutine reduce_dft_input_file(dft_code,input_filename,output_filename)
  implicit none
  
  type(String), intent(in) :: dft_code
  type(String), intent(in) :: input_filename
  type(String), intent(in) :: output_filename
  
  if (dft_code=='castep') then
    call reduce_castep_input_file(input_filename,output_filename)
  else
    call print_line('Reducing '//dft_code//' input files not yet supported.')
    call err()
  endif
end subroutine
end module
