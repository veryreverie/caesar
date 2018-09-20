! ======================================================================
! Reads and writes CASTEP .cell files, and reads CASTEP .castep files.
! ======================================================================
module castep_wrapper_submodule
  use utils_module
  
  use structure_module
  use normal_mode_module
  
  use electronic_structure_data_submodule
  implicit none
  
  private
  
  type :: CastepInputFile
    type(String), allocatable :: lattice_block(:)
    type(String), allocatable :: positions_block(:)
    type(String), allocatable :: masses_block(:)
    type(String), allocatable :: kpoints_block(:)
    type(String), allocatable :: remainder(:)
  end type
  
  public :: make_input_filename_castep
  public :: make_output_filename_castep
  public :: read_input_file_castep
  public :: write_input_file_castep
  public :: read_output_file_castep
contains

function make_input_filename_castep(seedname) result(output)
  implicit none
  
  type(String), intent(in) :: seedname
  type(String)             :: output
  
  output = seedname//'.cell'
end function

function make_output_filename_castep(seedname) result(output)
  implicit none
  
  type(String), intent(in) :: seedname
  type(String)             :: output
  
  output = seedname//'.castep'
end function

function read_input_file_castep(filename) result(output)
  implicit none
  
  type(String), intent(in)           :: filename
  type(BasicStructure)               :: output
  
  type(CastepInputFile) :: cell_file
  
  ! Lattice variables.
  logical  :: lattice_is_cart
  real(dp) :: conversion
  real(dp) :: lattice(3,3)
  real(dp) :: lengths(3)
  real(dp) :: angles(3)
  
  ! Atomic variables.
  logical                       :: positions_are_abs
  integer                       :: no_atoms
  type(String),     allocatable :: species(:)
  type(RealVector), allocatable :: positions(:)
  logical,          allocatable :: masses_found(:)
  real(dp),         allocatable :: masses(:)
  
  ! Temporary variables.
  type(String), allocatable :: line(:)
  integer                   :: i,j,ialloc
  
  cell_file = parse_castep_input_file(filename)
  
  ! Read data format.
  line = split_line(lower_case(cell_file%lattice_block(1)))
  lattice_is_cart = line(2)=='lattice_cart'
  
  line = split_line(lower_case(cell_file%positions_block(1)))
  positions_are_abs = line(2)=='positions_abs'
  
  ! Parse lattice.
  conversion = 1.0_dp / ANGSTROM_PER_BOHR
  j=0
  do i=2,size(cell_file%lattice_block)-1
    line = split_line(lower_case(cell_file%lattice_block(i)))
    
    ! Ignore comments.
    if (slice(line(1),1,1)=='!') then
      cycle
    
    ! Read units if present.
    elseif (line(1) == 'bohr') then
      conversion = 1.0_dp
    elseif (line(1) == 'a0') then
      conversion = 1.0_dp
    elseif (line(1) == 'm') then
      conversion = 1e10_dp / ANGSTROM_PER_BOHR
    elseif (line(1) == 'cm') then
      conversion = 1e8_dp / ANGSTROM_PER_BOHR
    elseif (line(1) == 'nm') then
      conversion = 1e1_dp / ANGSTROM_PER_BOHR
    elseif (line(1) == 'ang') then
      conversion = 1.0_dp / ANGSTROM_PER_BOHR
      
    elseif (lattice_is_cart) then
      j = j+1
      if (j>3) then
        call print_line(ERROR//': too many lattice lines found in '//filename)
        call err()
      endif
      lattice(j,:) = dble(line(1:3))*conversion
    
    else
      j = j+1
      if (j==1) then
        lengths = dble(line(1:3)) * conversion
      elseif (j==2) then
        angles = dble(line(1:3)) * PI/180
      else
        call print_line(ERROR//': too many lattice lines found in '//filename)
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
  conversion = 1/ANGSTROM_PER_BOHR
  j=0
  allocate( species(size(cell_file%positions_block)-2),   &
          & positions(size(cell_file%positions_block)-2), &
          & stat=ialloc); call err(ialloc)
  
  do i=2,size(cell_file%positions_block)-1
    line = split_line(lower_case(cell_file%positions_block(i)))
    
    ! Ignore comments.
    if (slice(line(1),1,1)=='!') then
      cycle
    
    ! Read units if present.
    elseif (line(1) == 'bohr' .or. line(1) == 'a0') then
      conversion = 1.0_dp
    elseif (line(1) == 'm') then
      conversion = 1e10_dp / ANGSTROM_PER_BOHR
    elseif (line(1) == 'cm') then
      conversion = 1e8_dp / ANGSTROM_PER_BOHR
    elseif (line(1) == 'nm') then
      conversion = 1e1_dp / ANGSTROM_PER_BOHR
    elseif (line(1) == 'ang') then
      conversion = 1.0_dp / ANGSTROM_PER_BOHR
    
    elseif (j>size(positions)) then
      call print_line(ERROR//': too many atom lines found in '//filename)
    
    ! Read in atomic positions.
    else
      j = j+1
      line = split_line(cell_file%positions_block(i)) ! N.B. no lower_case
      species(j) = line(1)
      positions(j) = dble(line(2:4))
    endif
  enddo
  
  no_atoms = j
  species = species(:no_atoms)
  positions = positions(:no_atoms)
  
  if (positions_are_abs) then
    do i=1,no_atoms
      positions(i) = positions(i) * conversion
    enddo
  else
    do i=1,no_atoms
      positions(i) = transpose(mat(lattice)) * positions(i)
    enddo
  endif
  
  ! Parse masses.
  conversion = KG_PER_AMU / KG_PER_ME
  allocate( masses_found(no_atoms), &
          & masses(no_atoms),       &
          & stat=ialloc); call err(ialloc)
  masses_found = .false.
  do i=2,size(cell_file%masses_block)-1
    line = split_line(lower_case(cell_file%masses_block(i)))
    
    ! Ignore comments.
    if (slice(line(1),1,1)=='!') then
      cycle
    
    ! Read units if present.
    elseif (line(1) == 'me') then
      conversion = 1.0_dp
    elseif (line(1) == 'amu') then
      conversion = KG_PER_AMU / KG_PER_ME
    elseif (line(1) == 'kg') then
      conversion = 1.0_dp / KG_PER_ME
    elseif (line(1) == 'g') then
      conversion = 1e-3_dp / KG_PER_ME
    
    ! Read in masses.
    else
      line = split_line(cell_file%masses_block(i)) ! N.B. no lower_case
      do j=1,no_atoms
        if (line(1)==species(j)) then
          masses_found(j) = .true.
          masses(j) = dble(line(2))*conversion
        endif
      enddo
    endif
  enddo
  
  if (.not. all(masses_found)) then
    call print_line(ERROR//': not all masses specified in '//filename)
    call err()
  endif
  
  ! Make output.
  output =  BasicStructure( mat(lattice), &
                            species,      &
                            masses,       &
                            positions)
end function

subroutine write_input_file_castep(structure,old_cell_filename, &
   & new_cell_filename)
  implicit none
  
  type(StructureData), intent(in)           :: structure
  type(String),        intent(in), optional :: old_cell_filename
  type(String),        intent(in)           :: new_cell_filename
  
  ! Band structure path data.
  type(RealVector) :: kpoint
  
  ! Old and new cell files.
  type(CastepInputFile) :: old_cell_file
  type(OFile)           :: new_cell_file
  
  ! Temporary variables.
  integer                   :: i
  type(String), allocatable :: line(:)
  
  ! --------------------------------------------------
  ! Parse old cell file.
  ! --------------------------------------------------
  if (present(old_cell_filename)) then
    old_cell_file = parse_castep_input_file(old_cell_filename)
    
    ! Transform electronic k-points from fractional primitive cell
    !    co-ordinates into fractional supercell co-ordinates.
    do i=2,size(old_cell_file%kpoints_block)-1
      line = split_line(old_cell_file%kpoints_block(i))
      kpoint = dble(line(1:3))
      kpoint = transpose(dblemat(structure%recip_supercell)) * kpoint
      old_cell_file%kpoints_block(i) = kpoint//' '//join(line(4:))
    enddo
  endif
  
  ! --------------------------------------------------
  ! Write new cell file.
  ! --------------------------------------------------
  new_cell_file = OFile(new_cell_filename)
  call new_cell_file%print_line('%block lattice_cart')
  call new_cell_file%print_line('bohr')
  call new_cell_file%print_lines(structure%lattice)
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
    do i=1,size(old_cell_file%kpoints_block)
      call new_cell_file%print_line(old_cell_file%kpoints_block(i))
    enddo
    
    call new_cell_file%print_line('')
    do i=1,size(old_cell_file%masses_block)
      call new_cell_file%print_line(old_cell_file%masses_block(i))
    enddo
    
    do i=1,size(old_cell_file%remainder)
      call new_cell_file%print_line(old_cell_file%remainder(i))
    enddo
  endif
end subroutine

! ----------------------------------------------------------------------
! Helper function which parse castep .cell files into blocks.
! ----------------------------------------------------------------------
function parse_castep_input_file(filename) result(output)
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
  integer :: kpoints_block_start
  integer :: kpoints_block_size
  integer :: remainder_size
  
  ! Temporary variables.
  type(String), allocatable :: line(:)
  integer                   :: i,j,k,ialloc
  
  if (.not. file_exists(filename)) then
    call print_line(ERROR//': '//filename//' does not exist. Please ensure &
       &that the provided run script produces a .castep file.')
    stop
  endif
  
  cell_file = IFile(filename)
  
  lattice_block_start = 0
  lattice_block_size = 0
  positions_block_start = 0
  positions_block_size = 0
  masses_block_start = 0
  masses_block_size = 0
  kpoints_block_start = 0
  kpoints_block_size = 0
  remainder_size = 0
  
  ! Work out line numbers.
  do i=1,size(cell_file)
    line = split_line(lower_case(cell_file%line(i)))
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
      elseif (line(1)=='%block' .and. line(2)=='bs_kpoint_path') then
        kpoints_block_start = i
      elseif (line(1)=='%endblock' .and. line(2)=='bs_kpoint_path') then
        kpoints_block_size = i-kpoints_block_start+1
      endif
    endif
  enddo
  
  remainder_size = size(cell_file)      &
               & - lattice_block_size   &
               & - positions_block_size &
               & - masses_block_size    &
               & - kpoints_block_size
  
  ! Check required blocks exist and all blocks are of reasonable sizes.
  if (lattice_block_start==0) then
    call print_line(ERROR//': lattice block not present in '//filename)
    call err()
  elseif (lattice_block_size <= 2) then
    call print_line(ERROR//': lattice block of unexpected size in '//filename)
    call err()
  elseif (positions_block_start==0) then
    call print_line(ERROR//': positions block not present in '//filename)
    call err()
  elseif (positions_block_size <= 2) then
    call print_line( ERROR//': positions block of unexpected size in '// &
                   & filename)
    call err()
  elseif (masses_block_start==0) then
    call print_line(ERROR//': species_mass block not present in '//filename)
    call err()
  elseif (masses_block_size <= 2) then
    call print_line(ERROR//': species_mass block of unexpected size in '// &
       & filename)
    call err()
  elseif (kpoints_block_start/=0 .and. kpoints_block_size<=2) then
    call print_line(ERROR//': bs_kpoint_path block of unexpected size in '// &
       & filename)
    call err()
  elseif (kpoints_block_start==0 .and. kpoints_block_size/=0) then
    call print_line(ERROR//': bs_kpoint_path block of unexpected size in '// &
       & filename)
    call err()
  endif
  
  ! Allocate space for the blocks.
  allocate( output%lattice_block(lattice_block_size),     &
          & output%positions_block(positions_block_size), &
          & output%masses_block(masses_block_size),       &
          & output%kpoints_block(kpoints_block_size),     &
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
    
    j = i-kpoints_block_start+1
    if (j>0 .and. j<kpoints_block_size) then
      output%kpoints_block(j) = cell_file%line(i)
      cycle
    endif
    
    k = k+1
    output%remainder(k) = cell_file%line(i)
  enddo
end function

function read_output_file_castep(filename,structure) result(output)
  implicit none
  
  type(String),        intent(in) :: filename
  type(StructureData), intent(in) :: structure
  type(ElectronicStructure)       :: output
  
  ! File contents.
  type(IFile)               :: castep_file
  type(String), allocatable :: line(:)
  logical,      allocatable :: atom_found(:)
  type(String)              :: species
  
  ! Line numbers.
  integer :: energy_line
  integer :: forces_start_line
  integer :: forces_end_line
  
  ! Output variables.
  real(dp)                      :: energy
  type(RealVector), allocatable :: forces(:)
  
  ! Temporary variables.
  integer :: i,j,ialloc
    
  if (.not. file_exists(filename)) then
    call print_line(ERROR//': '//filename//' does not exist. Please ensure &
       &that the provided run script produces a .castep file.')
    stop
  endif
  
  castep_file = IFile(filename)
  
  ! Work out line numbers.
  energy_line = 0
  forces_start_line = 0
  forces_end_line = 0
  do i=1,size(castep_file)
    line = split_line(lower_case(castep_file%line(i)))
    ! Energy.
    if (size(line)>=2) then
      if (line(1)=='final' .and. line(2)=='energy,') then
        energy_line = i
      endif
    endif
    ! Forces.
    if (size(line)>=2) then
      if ( line(1)==repeat('*',len(line(1))) .and. &
         & ( line(2)=='forces'      .or.           &
         &   line(2)=='symmetrised' .or.           &
         &   line(2)=='unconstrained'    )         &
         &                                         ) then
        forces_start_line = i
      endif
    endif
    if (size(line)==1) then
      if ( forces_start_line/=0 .and.        &
         & forces_end_line==0   .and.        &
         & line(1)==repeat('*',len(line(1))) ) then
        forces_end_line = i
      endif
    endif
  enddo
  
  ! Check line numbers.
  if (energy_line==0) then
    call print_line('Error: Energy not found in '//char(filename))
    stop
  endif
  if (forces_start_line==0) then
    call print_line('Error: Start of forces not found in '//char(filename))
    stop
  endif
  if (forces_end_line==0) then
    call print_line('Error: End of forces not found in '//char(filename))
    stop
  elseif (forces_end_line-forces_start_line-7/=structure%no_atoms) then
    call print_line(ERROR//': The number of atoms in the Castep output file &
       &does not match that in the input file.')
    call err()
  endif
  
  ! Allocate arrays.
  allocate( forces(structure%no_atoms),     &
          & atom_found(structure%no_atoms), &
          & stat=ialloc); call err(ialloc)
  atom_found = .false.
  
  ! Read data.
  line = split_line(castep_file%line(energy_line))
  energy = dble(line(5)) / EV_PER_HARTREE
  
  do i=1,structure%no_atoms
    line = split_line(castep_file%line(forces_start_line+5+i))
    species = line(2)
    j = first( structure%atoms%species()==species, &
             & mask    = .not.atom_found,          &
             & default = 0)
    if (j==0) then
      call print_line(ERROR//': The atoms in the Castep output file do not &
         &match those in the input file.')
      call err()
    endif
    atom_found(j) = .true.
    forces(j) = dble(line(4:6)) * ANGSTROM_PER_BOHR / EV_PER_HARTREE
  enddo
  
  if (.not. all(atom_found)) then
    call print_line(ERROR//': The atoms in the Castep output file do not &
       &match those in the input file.')
    call err()
  endif
  
  ! Construct output.
  output = ElectronicStructure(energy,CartesianForce(forces))
end function
end module
