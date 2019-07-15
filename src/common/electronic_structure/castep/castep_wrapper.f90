! ======================================================================
! Reads and writes CASTEP .cell files, and reads CASTEP .castep files.
! ======================================================================
module castep_wrapper_module
  use utils_module
  
  use structure_module
  use normal_mode_module
  
  use electronic_structure_data_module
  implicit none
  
  private
  
  public :: make_input_filename_castep
  public :: make_output_filename_castep
  public :: read_input_file_castep
  public :: write_input_file_castep
  public :: read_output_file_castep
  
  type :: CastepInputFile
    type(String), allocatable :: lattice_block(:)
    type(String), allocatable :: positions_block(:)
    type(String), allocatable :: masses_block(:)
    type(String), allocatable :: kpoints_block(:)
    type(String), allocatable :: remainder(:)
  end type
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
  
  type(String), intent(in) :: filename
  type(BasicStructure)     :: output
  
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
      positions(j) = vec(dble(line(2:4)))
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
      kpoint = vec(dble(line(1:3)))
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
  integer :: lattice_block_end
  integer :: positions_block_start
  integer :: positions_block_end
  integer :: masses_block_start
  integer :: masses_block_end
  integer :: kpoints_block_start
  integer :: kpoints_block_end
  
  ! Output variables.
  type(String), allocatable :: lattice_block(:)
  type(String), allocatable :: positions_block(:)
  type(String), allocatable :: masses_block(:)
  type(String), allocatable :: kpoints_block(:)
  type(String), allocatable :: remainder(:)
  
  ! Temporary variables.
  type(String), allocatable :: lines(:)
  type(String), allocatable :: line(:)
  integer                   :: i,j,k,ialloc
  
  if (.not. file_exists(filename)) then
    call print_line(ERROR//': '//filename//' does not exist. Please ensure &
       &that the provided run script produces a .castep file.')
    call quit()
  endif
  
  cell_file = IFile(filename)
  lines = cell_file%lines()
  
  ! Remove blank lines.
  lines = lines(filter([( len(lines(i))>0 ,i=1, size(lines) )]))
  
  ! Work out line numbers.
  lattice_block_start = 0
  lattice_block_end = 0
  positions_block_start = 0
  positions_block_end = 0
  masses_block_start = 0
  masses_block_end = 0
  kpoints_block_start = 0
  kpoints_block_end = 0
  do i=1,size(lines)
    line = split_line(lower_case(lines(i)))
    if (size(line) >= 2) then
      if (line(1)=='%block' .and. ( line(2)=='lattice_cart' .or. &
                                  & line(2)=='lattice_abc')) then
        lattice_block_start = i
      elseif (line(1)=='%endblock' .and. ( line(2)=='lattice_cart' .or. &
                                         & line(2)=='lattice_abc')) then
        lattice_block_end = i
      elseif (line(1)=='%block' .and. ( line(2)=='positions_abs' .or. &
                                      & line(2)=='positions_frac')) then
        positions_block_start = i
      elseif (line(1)=='%endblock' .and. ( line(2)=='positions_abs' .or. &
                                         & line(2)=='positions_frac')) then
        positions_block_end = i
      elseif (line(1)=='%block' .and. line(2)=='species_mass') then
        masses_block_start = i
      elseif (line(1)=='%endblock' .and. line(2)=='species_mass') then
        masses_block_end = i
      elseif (line(1)=='%block' .and. line(2)=='bs_kpoint_path') then
        kpoints_block_start = i
      elseif (line(1)=='%endblock' .and. line(2)=='bs_kpoint_path') then
        kpoints_block_end = i
      endif
    endif
  enddo
  
  ! Check required blocks exist and all blocks are of reasonable sizes.
  if (lattice_block_start==0) then
    call print_line(ERROR//': lattice block not present in '//filename)
    call err()
  elseif (lattice_block_end-lattice_block_start < 3) then
    call print_line(ERROR//': lattice block of unexpected size in '//filename)
    call err()
  elseif (positions_block_start==0) then
    call print_line(ERROR//': positions block not present in '//filename)
    call err()
  elseif (positions_block_end-positions_block_start < 2) then
    call print_line( ERROR//': positions block of unexpected size in '// &
                   & filename)
    call err()
  elseif (masses_block_start==0) then
    call print_line(ERROR//': species_mass block not present in '//filename)
    call err()
  elseif (masses_block_end-masses_block_start < 2) then
    call print_line(ERROR//': species_mass block of unexpected size in '// &
       & filename)
    call err()
  elseif ( kpoints_block_start/=0 .and.             &
         & kpoints_block_end-kpoints_block_start <2 ) then
    call print_line(ERROR//': bs_kpoint_path block of unexpected size in '// &
       & filename)
    call err()
  endif
  
  ! Copy over blocks.
  lattice_block = lines(lattice_block_start:lattice_block_end)
  positions_block = lines(positions_block_start:positions_block_end)
  masses_block = lines(masses_block_start:masses_block_end)
  if (kpoints_block_start==0) then
    allocate(kpoints_block(0), stat=ialloc); call err(ialloc)
  else
    kpoints_block = lines(kpoints_block_start:kpoints_block_end)
  endif
  
  ! Remove comment lines from blocks.
  lattice_block = lattice_block(filter(                          &
     & [( all(char(slice(lattice_block(i),1,1))/=['!','#',';']), &
     &    i=1,                                                   &
     &    size(lattice_block)                                    )]))
  positions_block = positions_block(filter(                        &
     & [( all(char(slice(positions_block(i),1,1))/=['!','#',';']), &
     &    i=1,                                                     &
     &    size(positions_block)                                    )]))
  masses_block = masses_block(filter(                           &
     & [( all(char(slice(masses_block(i),1,1))/=['!','#',';']), &
     &    i=1,                                                  &
     &    size(masses_block)                                    )]))
  kpoints_block = kpoints_block(filter(                          &
     & [( all(char(slice(kpoints_block(i),1,1))/=['!','#',';']), &
     &    i=1,                                                   &
     &    size(kpoints_block)                                    )]))
  
  ! Copy across remainder of file.
  allocate(remainder(0), stat=ialloc); call err(ialloc)
  do i=1,size(lines)
    if (i>=lattice_block_start .and. i<=lattice_block_end) then
      cycle
    elseif (i>=positions_block_start .and. i<=positions_block_end) then
      cycle
    elseif (i>=masses_block_start .and. i<=masses_block_end) then
      cycle
    elseif (i>=kpoints_block_start .and. i<=kpoints_block_end) then
      cycle
    else
      remainder = [remainder, lines(i)]
    endif
  enddo
  
  output = CastepInputFile( lattice_block   = lattice_block,   &
                          & positions_block = positions_block, &
                          & masses_block    = masses_block,    &
                          & kpoints_block   = kpoints_block,   &
                          & remainder       = remainder        )
end function

function read_output_file_castep(directory,seedname,structure,use_forces, &
   & use_hessians,calculate_stress) result(output)
  implicit none
  
  type(String),        intent(in) :: directory
  type(String),        intent(in) :: seedname
  type(StructureData), intent(in) :: structure
  logical,             intent(in) :: use_forces
  logical,             intent(in) :: use_hessians
  logical,             intent(in) :: calculate_stress
  type(ElectronicStructure)       :: output
  
  ! File contents.
  type(String)              :: filename
  type(IFile)               :: castep_file
  type(String), allocatable :: line(:)
  logical,      allocatable :: atom_found(:)
  type(String)              :: species
  
  ! Line numbers.
  integer :: energy_line
  integer :: forces_start_line
  integer :: forces_end_line
  integer :: stress_line
  integer :: permittivity_line
  integer :: born_charges_line
  
  ! Output variables.
  real(dp)                          :: energy
  type(RealVector),     allocatable :: forces_elements(:)
  type(CartesianForce), allocatable :: forces
  real(dp)                          :: stress_elements(3,3)
  type(RealMatrix),     allocatable :: stress
  real(dp)                          :: permittivity(3,3)
  real(dp)                          :: born_charge(3,3)
  type(RealMatrix),     allocatable :: born_charges(:)
  type(LinearResponse), allocatable :: linear_response
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
    
  filename = directory//'/'//make_output_filename_castep(seedname)
  if (.not. file_exists(filename)) then
    call print_line(ERROR//': '//filename//' does not exist. Please ensure &
       &that the provided run script produces a .castep file.')
    call quit()
  endif
  castep_file = IFile(filename)
  
  ! Work out line numbers.
  energy_line = 0
  forces_start_line = 0
  forces_end_line = 0
  stress_line = 0
  permittivity_line = 0
  born_charges_line = 0
  do i=1,size(castep_file)
    line = split_line(lower_case(castep_file%line(i)))
    
    ! Energy.
    ! N.B. this line likely appears in the file three times for reasons of
    !    finite basis correction. The third time is the final answer.
    if (size(line)>=2) then
      if (line(1)=='final' .and. line(2)=='energy,') then
        energy_line = i
      endif
    endif
    
    ! Forces.
    if (forces_start_line==0) then
      if (size(line)>=2) then
        if ( line(1)==repeat('*',len(line(1))) .and. &
           & ( line(2)=='forces'      .or.           &
           &   line(2)=='symmetrised' .or.           &
           &   line(2)=='constrained' .or.           &
           &   line(2)=='unconstrained'    )         &
           &                                         ) then
          forces_start_line = i
        endif
      endif
    endif
    
    if (forces_start_line/=0 .and. forces_end_line==0) then
      if (size(line)==1) then
        if ( forces_start_line/=0 .and.        &
           & forces_end_line==0   .and.        &
           & line(1)==repeat('*',len(line(1))) ) then
          forces_end_line = i
        endif
      endif
    endif
    
    ! Stress tensor.
    if (stress_line==0) then
      if (size(line)>=4) then
        if ( line(1)==repeat('*',len(line(1))) .and. &
           & line(2)=='stress'                 .and. &
           & line(3)=='tensor'                       ) then
          stress_line = i
        endif
      endif
    endif
    
    ! Linear response (permittivity and Born effective charges).
    if (permittivity_line==0) then
      if (size(line)>=3) then
        if ( line(1)=='optical'       .and. &
           & line(2)=='permittivity'  .and. &
           & line(3)=='(f->infinity)'       ) then
          permittivity_line = i
        endif
      endif
    endif
    
    if (born_charges_line==0) then
      if (size(line)>=3) then
        if ( line(1)=='born'      .and. &
               & line(2)=='effective' .and. &
               & line(3)=='charges'         ) then
          born_charges_line = i
        endif
      endif
    endif
  enddo
  
  ! Check line numbers.
  if (energy_line==0) then
    call print_line('Error: Energy not found in '//char(filename))
    call quit()
  elseif (permittivity_line==0 .neqv. born_charges_line==0) then
    call print_line(ERROR//': DFPT linear response is only partially present &
       &in Castep output file.')
    call err()
  endif
  
  ! Read energy.
  line = split_line(castep_file%line(energy_line))
  energy = dble(line(5)) / EV_PER_HARTREE
  
  ! Read forces and Born effective charges.
  if (use_forces) then
    if (forces_start_line==0) then
      call print_line('Error: Start of forces not found in '//char(filename))
      call quit()
    elseif (forces_end_line==0) then
      call print_line('Error: End of forces not found in '//char(filename))
      call quit()
    elseif (forces_end_line-forces_start_line-7/=structure%no_atoms) then
      call print_line(forces_start_line)
      call print_line(forces_end_line)
      call print_line(forces_end_line-forces_start_line-7//' '// &
                     &structure%no_atoms)
      call print_line(ERROR//': The number of atoms in the Castep output file &
         &does not match that in the input file.')
      call err()
    endif
  endif
  allocate( forces_elements(structure%no_atoms), &
          & born_charges(structure%no_atoms),    &
          & atom_found(structure%no_atoms),      &
          & stat=ialloc); call err(ialloc)
  atom_found = .false.
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
    
    if (use_forces) then
      forces_elements(j) = vec( dble(line(4:6))   &
                            & * ANGSTROM_PER_BOHR &
                            & / EV_PER_HARTREE    )
    endif
    
    if (born_charges_line/=0) then
      do k=1,3
        line = split_line(castep_file%line(born_charges_line+1+(i-1)*3+k))
        if (k==1) then
          if (line(1)/=species) then
            call print_line(ERROR//': The Born charges and Forces in the &
               &Castep output file do not match.')
            call err()
          endif
          born_charge(k,:) = dble(line(3:5))
        else
          born_charge(k,:) = dble(line(1:3))
        endif
      enddo
      born_charges(j) = mat(born_charge)
    endif
  enddo
  
  if (.not. all(atom_found)) then
    call print_line(ERROR//': The atoms in the Castep output file do not &
       &match those in the input file.')
    call err()
  endif
  
  if (use_forces) then
    forces = CartesianForce(forces_elements)
  endif
  
  ! Read Hessian.
  if (use_hessians) then
    ! TODO
    call print_line(ERROR//': Reading Hessians from Castep not yet &
       &implemented.')
    call err()
  endif
  
  ! Read stress.
  if (calculate_stress) then
    if (stress_line==0) then
      call print_line(ERROR//': No stress found in Castep output file.')
      call err()
    endif
    
    do i=1,3
      line = split_line(castep_file%line(stress_line+5+i))
      ! Stress is in GPa = 10^9 * J.m^-3 = 10^-21 * J.A^-3
      ! The a.u. unit is Ha.bohr^-3
      stress_elements (i,:) = dble(line(3:5)) &
                          & * 1e-21_dp*ANGSTROM_PER_BOHR**3/JOULES_PER_HARTREE
    enddo
    
    stress = mat(stress_elements)
  endif
  
  ! Read permittivity.
  if (permittivity_line/=0) then
    do i=1,3
      line = split_line(castep_file%line(permittivity_line+1+i))
      permittivity(i,:) = dble(line(1:3))
    enddo
    
    linear_response = LinearResponse(mat(permittivity), born_charges)
  endif
  
  ! Construct output.
  output = ElectronicStructure( energy          = energy,         &
                              & forces          = forces,         &
                              & stress          = stress,         &
                              & linear_response = linear_response )
end function
end module
