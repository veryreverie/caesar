! ======================================================================
! Reads and writes Quantum Espresso PWscf input and output files.
! ======================================================================
module qe_pw_module
  use utils_module
  
  use structure_module
  use normal_mode_module
  
  use electronic_structure_data_module
  use electronic_structure_common_module
  
  use qe_fc_module
  implicit none
  
  private
  
  public :: make_input_filename_qe
  public :: make_output_filename_qe
  public :: read_input_file_qe
  public :: write_input_file_qe
  public :: read_output_file_qe
  
  public :: QeInputFile
  
  type, extends(Stringsable) :: QeInputFile
    type(String), allocatable :: namelists(:)
    type(String), allocatable :: atomic_species(:)
    type(String), allocatable :: atomic_positions(:)
    type(String), allocatable :: k_points(:)
    type(String), allocatable :: cell_parameters(:)
    type(String), allocatable :: occupations(:)
    type(String), allocatable :: constraints(:)
    type(String), allocatable :: atomic_forces(:)
  contains
    procedure, public :: read  => read_QeInputFile
    procedure, public :: write => write_QeInputFile
  end type
  
  interface QeInputFile
    module procedure new_QeInputFile
    module procedure new_QeInputFile_Strings
    module procedure new_QeInputFile_StringArray
  end interface
contains

function make_input_filename_qe(seedname) result(output)
  implicit none
  
  type(String), intent(in) :: seedname
  type(String)             :: output
  
  output = seedname//'.in'
end function

function make_output_filename_qe(seedname) result(output)
  implicit none
  
  type(String), intent(in) :: seedname
  type(String)             :: output
  
  output = seedname//'.out'
end function

function read_input_file_qe(filename) result(output)
  implicit none
  
  type(String), intent(in) :: filename
  type(BasicStructure)     :: output
  
  ! File data.
  type(IFile)       :: input_file
  type(QeInputFile) :: qe_file
  
  ! Required variables.
  type(RealMatrix)              :: lattice
  type(String),     allocatable :: species(:)
  type(RealVector), allocatable :: positions(:)
  real(dp),         allocatable :: masses(:)
  logical,          allocatable :: masses_set(:)
  
  ! Temporary variables.
  type(String), allocatable :: line(:)
  integer                   :: i,ialloc
  
  ! Read file.
  input_file = IFile(filename)
  qe_file = QeInputFile(input_file%lines())
  
  ! Parse lattice.
  if (.not. allocated(qe_file%cell_parameters)) then
    call print_line(ERROR//': Caesar requires cell_parameters card.')
    call quit()
  endif
  line = split_line(lower_case(qe_file%cell_parameters(1)))
  if (size(line)==1) then
    call print_line(ERROR//': No unit given for cell_parameters card. This is &
       &deprecated behaviour.')
    call quit()
  elseif ( line(2)=='bohr'   .or. &
         & line(2)=='(bohr)' .or. &
         & line(2)=='{bohr}'      ) then
    lattice = RealMatrix(qe_file%cell_parameters(2:4))
  elseif ( line(2)=='angstrom'   .or. &
         & line(2)=='(angstrom)' .or. &
         & line(2)=='{angstrom}'      ) then
    lattice = RealMatrix(qe_file%cell_parameters(2:4)) / ANGSTROM_PER_BOHR
  elseif ( line(2)=='alat'   .or. &
         & line(2)=='(alat)' .or. &
         & line(2)=='{alat}'      ) then
    call print_line(ERROR//': Caesar cannot parse cell_parameters card in &
       &"alat" format.')
    call quit()
  else
    call print_line(ERROR//': Unrecognised cell_parameters card format: '// &
       & line(2))
    call quit()
  endif
  
  ! Parse atomic positions.
  line = split_line(lower_case(qe_file%atomic_positions(1)))
  if (size(line)==1) then
    call print_line(ERROR//': No unit given for atomic_positions card. This &
       &is deprecated behaviour.')
    call quit()
  elseif ( line(2)=='alat'   .or. &
         & line(2)=='(alat)' .or. &
         & line(2)=='{alat}'      ) then
    call print_line(ERROR//': Caesar cannot parse atomic_positions card in &
       &"alat" format.')
    call quit()
  elseif ( line(2)=='bohr'   .or. &
         & line(2)=='(bohr)' .or. &
         & line(2)=='{bohr}'      ) then
    allocate( species(size(qe_file%atomic_positions)-1),   &
            & positions(size(qe_file%atomic_positions)-1), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(species)
      line = split_line(qe_file%atomic_positions(i+1))
      species(i) = line(1)
      positions(i) = vec(dble(line(2:4)))
    enddo
  elseif ( line(2)=='angstrom'   .or. &
         & line(2)=='(angstrom)' .or. &
         & line(2)=='{angstrom}'      ) then
    allocate( species(size(qe_file%atomic_positions)-1),   &
            & positions(size(qe_file%atomic_positions)-1), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(species)
      line = split_line(qe_file%atomic_positions(i+1))
      species(i) = line(1)
      positions(i) = vec(dble(line(2:4)) / ANGSTROM_PER_BOHR)
    enddo
  elseif ( line(2)=='crystal'   .or. &
         & line(2)=='(crystal)' .or. &
         & line(2)=='{crystal}'      ) then
    allocate( species(size(qe_file%atomic_positions)-1),   &
            & positions(size(qe_file%atomic_positions)-1), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(species)
      line = split_line(qe_file%atomic_positions(i+1))
      species(i) = line(1)
      positions(i) = vec(dble(line(2:4)))
      positions(i) = transpose(lattice) * positions(i)
    enddo
  elseif ( line(2)=='crystal_sg'   .or. &
         & line(2)=='(crystal_sg)' .or. &
         & line(2)=='{crystal_sg}'      ) then
    call print_line(ERROR//': Caesar cannot parse atomic_positions card in &
       &crystal_sg format.')
    call quit()
  else
    call print_line(ERROR//': Unrecognised atomic_positions card format: '// &
       & line(2))
    call quit()
  endif
  
  ! Parse masses.
  allocate(masses(size(species)), stat=ialloc); call err(ialloc)
  masses_set = [(.false., i=1, size(masses))]
  do i=2,size(qe_file%atomic_species)
    line = split_line(qe_file%atomic_species(i))
    masses(filter(species==line(1))) = dble(line(2))*KG_PER_AMU/KG_PER_ME
    masses_set(filter(species==line(1))) = .true.
  enddo
  if (.not. all(masses_set)) then
    call print_line(ERROR//': Not all species present in atomic_species card.')
    call quit()
  endif
  
  ! Make output.
  output = BasicStructure( lattice,  &
                         & species,  &
                         & masses,   &
                         & positions )
end function

subroutine write_input_file_qe(structure,old_qe_in_filename,new_qe_in_filename)
  implicit none
  
  type(StructureData), intent(in)           :: structure
  type(String),        intent(in), optional :: old_qe_in_filename
  type(String),        intent(in)           :: new_qe_in_filename
  
  type(IFile)       :: old_qe_in_file
  type(OFile)       :: new_qe_in_file
  type(QeInputFile) :: qe_file
  
  type(KpointGrid) :: kpoint_grid
  real(dp)         :: kpoint_spacing
  
  type(String), allocatable :: line(:)
  
  integer :: i
  
  if (.not. present(old_qe_in_filename)) then
    call print_line(ERROR//': Unable to write QE input file without QE input &
       &file to work from.')
    call err()
  endif
  
  ! Read in old QE file.
  old_qe_in_file = IFile(old_qe_in_filename)
  qe_file = QeInputFile(old_qe_in_file%lines())
  
  ! Change cell_parameters and atomic_positions cards to match new structure.
  qe_file%cell_parameters = [ str('CELL_PARAMETERS bohr'), &
                            & str(structure%lattice)       ]
  
  qe_file%atomic_positions = [                         &
     & str('ATOMIC_POSITIONS bohr'),                   &
     & ( structure%atoms(i)%species() //' '//          &
     &      structure%atoms(i)%cartesian_position(),   &
     &   i=1,                                          &
     &   structure%no_atoms                          ) ]
  
  ! Update 'nat' in namelists.
  do i=1,size(qe_file%namelists)
    if (len(qe_file%namelists(i))>=4) then
      if (slice(qe_file%namelists(i),1,4)=='nat=') then
        qe_file%namelists(i) = 'nat='//structure%no_atoms
      endif
    endif
  enddo
  
  ! Update k-point grid.
  line = split_line(lower_case(qe_file%k_points(1)))
  if ( line(2)/='automatic'   .and. &
     & line(2)/='(automatic)' .and. &
     & line(2)/='{automatic}'       ) then
    call print_line(ERROR//': Unable to generate supercell k-point grid for &
       &Quantum Espresso calculation whose k_points card is not in &
       &"automatic" format.')
    call quit()
  endif
  
  line = split_line(qe_file%k_points(2))
  kpoint_grid = KpointGrid(int(line(1:3)))
  kpoint_spacing = calculate_kpoint_spacing( &
     & kpoint_grid,                          &
     & dble(structure%prim_recip_lattice())  )
  kpoint_grid = calculate_kpoint_grid( &
     & kpoint_spacing,                 &
     & dble(structure%recip_lattice)   )
  
  if (size(line)==3) then
    qe_file%k_points(2) = str(kpoint_grid)
  elseif (size(line)==6) then
    qe_file%k_points(2) = join([str(kpoint_grid), line(4:6)])
  else
    call print_line(ERROR//': k_points card has an unexpected number of &
       &entries.')
    call quit()
  endif
  
  
  ! Write out new QE file.
  new_qe_in_file = OFile(new_qe_in_filename)
  call new_qe_in_file%print_lines(qe_file)
end subroutine

function read_output_file_qe(directory,seedname,structure,use_forces, &
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
  type(IFile)               :: qe_file
  type(String), allocatable :: line(:)
  
  ! Line numbers.
  integer :: species_start_line
  integer :: energy_line
  integer :: forces_start_line
  integer :: stress_start_line
  
  ! Counts.
  integer :: no_species
  integer :: no_forces
  
  ! QE 'type' to species conversion.
  type(String), allocatable :: species(:)
  logical,      allocatable :: forces_found(:)
  
  ! Output variables.
  real(dp)                            :: energy
  type(RealVector),       allocatable :: forces_elements(:)
  type(CartesianForce),   allocatable :: forces
  real(dp)                            :: stress_elements(3,3)
  type(RealMatrix),       allocatable :: stress
  type(CartesianHessian), allocatable :: hessian
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  filename = directory//'/'//make_output_filename_qe(seedname)
  qe_file = IFile(filename)
  
  ! Work out line numbers.
  ! N.B. start lines are one before the start of the block.
  ! end lines are the last lines in the block.
  species_start_line = 0
  no_species = 0
  energy_line = 0
  forces_start_line = 0
  no_forces = 0
  stress_start_line = 0
  
  do i=1,size(qe_file)
    line = split_line(lower_case(qe_file%line(i)))
    
    ! Species.
    if (size(line)>=2) then
      if (line(1)=='atomic' .and. line(2)=='species') then
        species_start_line = i
      endif
    endif
    
    if (species_start_line/=0 .and. no_species==0 .and. size(line)==0) then
      no_species = i-1-species_start_line
    endif
    
    ! Energy.
    if (size(line)>=1) then
      if (line(1)=='!') then
        energy_line=i
      endif
    endif
    
    ! Forces.
    if (size(line)>=2) then
      if (line(1)=='forces' .and. line(2)=='acting') then
        forces_start_line = i+1
      elseif (line(1)=='total' .and. line(2)=='force') then
        no_forces = i-2-forces_start_line
      endif
    endif
    
    ! Stresses
    if (size(line)>=2) then
      if (line(1)=='total' .and. line(2)=='stress') then
        stress_start_line = i
      endif
    endif
  enddo
  
  ! Read energy.
  line = split_line(qe_file%line(energy_line))
  energy = dble(line(5)) / RYDBERG_PER_HARTREE
  
  ! Read forces.
  if (use_forces) then
    if (forces_start_line==0) then
      call print_line(ERROR//': No forces found in Quantum Espresso output &
         &file.')
      call quit()
    elseif (no_forces/=structure%no_atoms) then
      call print_line(ERROR//': The number of atoms in the Quantum Espresso &
         &output file does not match that in the input file.')
      call quit()
    endif
    
    allocate( species(no_species),                 &
            & forces_elements(structure%no_atoms), &
            & stat=ialloc); call err(ialloc)
    forces_found = [(.false.,i=1,structure%no_atoms)]
    do i=1,no_species
      line = split_line(qe_file%line(species_start_line+i))
      species(i) = line(1)
    enddo
    
    do i=1,no_forces
      line = split_line(qe_file%line(forces_start_line+i))
      
      j = first( structure%atoms%species()==species(int(line(4))), &
               & mask=.not.forces_found,                           &
               & default=0                                         )
      
      if (j==0) then
        call print_line(ERROR//': Unable to match species in Quantum Espresso &
           &output file with those in input file.')
        call quit()
      endif
      
      forces_elements(j) = vec(dble(line(7:9)) / RYDBERG_PER_HARTREE)
      forces_found(j) = .true.
    enddo
    
    if (.not.all(forces_found)) then
      call print_line(ERROR//': Unable to match forces in Quantum Espresso &
         &output file with atoms in input file.')
      call quit()
    endif
    
    forces = CartesianForce(forces_elements)
  endif
  
  ! Read Hessian.
  if (use_hessians) then
    ! N.B. only the Hessian for the primitive cell is needed,
    !    so the supercell is just structure.
    hessian = read_qe_force_constants_file( directory = directory, &
                                          & seedname  = seedname,  &
                                          & supercell = structure  )
  endif
  
  ! Read stress.
  if (calculate_stress) then
    do i=1,3
      line = split_line(qe_file%line(stress_start_line+i))
      stress_elements(i,:) = dble(line(1:3)) / RYDBERG_PER_HARTREE
    enddo
    
    stress = mat(stress_elements)
  endif
  
  ! Construct output.
  output = ElectronicStructure( energy  = energy,  &
                              & forces  = forces,  &
                              & hessian = hessian, &
                              & stress  = stress   )
end function

! ----------------------------------------------------------------------
! QeInputFile methods.
! ----------------------------------------------------------------------
function new_QeInputFile(namelists,atomic_species,atomic_positions,k_points, &
   & cell_parameters,occupations,constraints,atomic_forces) result(this)
  implicit none
  
  type(String), intent(in)           :: namelists(:)
  type(String), intent(in)           :: atomic_species(:)
  type(String), intent(in)           :: atomic_positions(:)
  type(String), intent(in)           :: k_points(:)
  type(String), intent(in), optional :: cell_parameters(:)
  type(String), intent(in), optional :: occupations(:)
  type(String), intent(in), optional :: constraints(:)
  type(String), intent(in), optional :: atomic_forces(:)
  type(QeInputFile)                  :: this
  
  this%namelists = namelists
  this%atomic_species = atomic_species
  this%atomic_positions = atomic_positions
  this%k_points = k_points
  if (present(cell_parameters)) then
    this%cell_parameters = cell_parameters
  endif
  if (present(occupations)) then
    this%occupations = occupations
  endif
  if (present(constraints)) then
    this%constraints = constraints
  endif
  if (present(atomic_forces)) then
    this%atomic_forces = atomic_forces
  endif
end function

subroutine read_QeInputFile(this,input)
  implicit none
  
  class(QeInputFile), intent(out) :: this
  type(String),       intent(in)  :: input(:)
  
  type(String), allocatable :: lines(:)
  
  type(String), allocatable :: namelists(:)
  type(String), allocatable :: atomic_species(:)
  type(String), allocatable :: atomic_positions(:)
  type(String), allocatable :: k_points(:)
  type(String), allocatable :: cell_parameters(:)
  type(String), allocatable :: occupations(:)
  type(String), allocatable :: constraints(:)
  type(String), allocatable :: atomic_forces(:)
  
  integer, allocatable :: namelists_end_lines(:)
  
  integer :: namelists_end_line
  integer :: atomic_species_line
  integer :: atomic_positions_line
  integer :: k_points_line
  integer :: cell_parameters_line
  integer :: occupations_line
  integer :: constraints_line
  integer :: atomic_forces_line
  
  integer, allocatable :: end_lines(:)
  integer              :: end_line
  
  type(String), allocatable :: line(:)
  
  integer :: i,j,ialloc
  
  select type(this); type is(QeInputFile)
    ! Remove blank lines and comment lines.
    allocate(lines(size(input)), stat=ialloc); call err(ialloc)
    j = 0
    do i=1,size(input)
      line = tokens(input(i))
      if (size(line)==0) then
        cycle
      elseif (slice(line(1),1,1)=='!' .or. slice(line(1),1,1)=='#') then
        cycle
      endif
      j = j+1
      lines(j) = input(i)
    enddo
    lines = lines(:j)
    
    ! Locate lines.
    namelists_end_lines = filter([(lines(i)=='/',i=1,size(lines))])
    namelists_end_line = namelists_end_lines(size(namelists_end_lines))
    
    atomic_species_line = 0
    atomic_positions_line = 0
    k_points_line = 0
    cell_parameters_line = 0
    occupations_line = 0
    constraints_line = 0
    atomic_forces_line = 0
    
    end_lines = [integer::]
    
    do i=namelists_end_line+1,size(lines)
      line = split_line(lower_case(lines(i)))
      if (size(line)==0) then
        cycle
      endif
      
      if (line(1)=='atomic_species') then
        if (atomic_species_line/=0) then
          call print_line(ERROR//': atomic_species card appears twice.')
          call quit()
        endif
        atomic_species_line = i
        end_lines = [end_lines, i-1]
      elseif (line(1)=='atomic_positions') then
        if (atomic_positions_line/=0) then
          call print_line(ERROR//': atomic_positions card appears twice.')
          call quit()
        endif
        atomic_positions_line = i
        end_lines = [end_lines, i-1]
      elseif (line(1)=='k_points') then
        if (k_points_line/=0) then
          call print_line(ERROR//': k_points card appears twice.')
          call quit()
        endif
        k_points_line = i
        end_lines = [end_lines, i-1]
      elseif (line(1)=='cell_parameters') then
        if (cell_parameters_line/=0) then
          call print_line(ERROR//': cell_parameters card appears twice.')
          call quit()
        endif
        cell_parameters_line = i
        end_lines = [end_lines, i-1]
      elseif (line(1)=='occupations') then
        if (occupations_line/=0) then
          call print_line(ERROR//': occupations card appears twice.')
          call quit()
        endif
        occupations_line = i
        end_lines = [end_lines, i-1]
      elseif (line(1)=='constraints') then
        if (constraints_line/=0) then
          call print_line(ERROR//': constraints card appears twice.')
          call quit()
        endif
        constraints_line = i
        end_lines = [end_lines, i-1]
      elseif (line(1)=='atomic_forces') then
        if (atomic_forces_line/=0) then
          call print_line(ERROR//': atomic_forces card appears twice.')
          call quit()
        endif
        atomic_forces_line = i
        end_lines = [end_lines, i-1]
      endif
    enddo
    
    end_lines = [end_lines, size(lines)]
    
    namelists = lines(:namelists_end_line)
    
    if (atomic_species_line==0) then
      call print_line(ERROR//': atomic_species card not present.')
    else
      end_line = minval(end_lines(filter(end_lines>=atomic_species_line)))
      atomic_species = lines(atomic_species_line:end_line)
    endif
    
    if (atomic_positions_line==0) then
      call print_line(ERROR//': atomic_positions card not present.')
    else
      end_line = minval(end_lines(filter(end_lines>=atomic_positions_line)))
      atomic_positions = lines(atomic_positions_line:end_line)
    endif
    
    if (k_points_line==0) then
      call print_line(ERROR//': k_points card not present.')
    else
      end_line = minval(end_lines(filter(end_lines>=k_points_line)))
      k_points = lines(k_points_line:end_line)
    endif
    
    this = QeInputFile(namelists, atomic_species, atomic_positions, k_points)
    
    if (cell_parameters_line/=0) then
      end_line = minval(end_lines(filter(end_lines>=cell_parameters_line)))
      cell_parameters = lines(cell_parameters_line:end_line)
      this%cell_parameters = cell_parameters
    endif
    
    if (occupations_line/=0) then
      end_line = minval(end_lines(filter(end_lines>=occupations_line)))
      occupations = lines(occupations_line:end_line)
      this%occupations = occupations
    endif
    
    if (constraints_line/=0) then
      end_line = minval(end_lines(filter(end_lines>=constraints_line)))
      constraints = lines(constraints_line:end_line)
      this%constraints = constraints
    endif
    
    if (atomic_forces_line/=0) then
      end_line = minval(end_lines(filter(end_lines>=atomic_forces_line)))
      atomic_forces = lines(atomic_forces_line:end_line)
      this%atomic_forces = atomic_forces
    endif
  class default
    call err()
  end select
end subroutine

function write_QeInputFile(this) result(output)
  implicit none
  
  class(QeInputFile), intent(in) :: this
  type(String), allocatable      :: output(:)
  
  select type(this); type is(QeInputFile)
    output = [ this%namelists,        &
             & this%atomic_species,   &
             & this%atomic_positions, &
             & this%k_points          ]
    if (allocated(this%cell_parameters)) then
      output = [output, this%cell_parameters]
    endif
    if (allocated(this%occupations)) then
      output = [output, this%occupations]
    endif
    if (allocated(this%constraints)) then
      output = [output, this%constraints]
    endif
    if (allocated(this%atomic_forces)) then
      output = [output, this%atomic_forces]
    endif
  class default
    call err()
  end select
end function

function new_QeInputFile_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(QeInputFile)        :: this
  
  call this%read(input)
end function

impure elemental function new_QeInputFile_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(QeInputFile)             :: this
  
  this = QeInputFile(str(input))
end function
end module
