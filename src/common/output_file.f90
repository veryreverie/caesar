! ======================================================================
! A class for holding the information in a castep .castep or qe .out file.
! ======================================================================
module output_file_module
  use utils_module
  
  use physical_constants_module
  use structure_module
  use input_file_module
  use quip_wrapper_module
  implicit none
  
  private
  
  public :: make_output_filename
  public :: OutputFile
  public :: read_output_file
  
  type OutputFile
    integer                       :: no_atoms
    type(String), allocatable     :: species(:)
    real(dp)                      :: energy
    type(RealVector), allocatable :: forces(:)
  end type
  
  interface OutputFile
    module procedure new_OutputFile
  end interface
  
contains

! ----------------------------------------------------------------------
! Converts a file seedname into the appropriate dft input or output filename.
! ----------------------------------------------------------------------
function make_output_filename(file_type,seedname) result(output)
  implicit none
  
  type(String), intent(in) :: file_type
  type(String), intent(in) :: seedname
  type(String)             :: output
  
  if (file_type == 'castep') then
    output = seedname//'.castep'
  elseif (file_type == 'vasp') then
    output = 'OUTCAR'
  elseif (file_type == 'qe') then
    output = seedname//'.out'
  elseif (file_type == 'quip') then
    ! N.B. the '.cell' is intentional. When running with QUIP, no output file
    !    is produced. Instead the input file is read back in, and the energy
    !    etc. are calculated on the fly.
    output = seedname//'.cell'
  else
    call print_line('Unrecognised dft code: '//file_type)
    call err()
  endif
end function

function new_OutputFile(no_atoms) result(this)
  implicit none
  
  integer, intent(in) :: no_atoms
  type(OutputFile)    :: this
  
  integer :: ialloc
  
  this%no_atoms = no_atoms
  allocate( this%species(no_atoms), &
          & this%forces(no_atoms),  &
          & stat=ialloc); call err(ialloc)
end function

function read_castep_output_file(filename,structure) result(output)
  implicit none
  
  type(String),        intent(in) :: filename
  type(StructureData), intent(in) :: structure
  type(OutputFile)                :: output
  
  ! file contents
  type(IFile)                   :: castep_file
  type(String),     allocatable :: line(:)
  type(String)                  :: species
  type(RealVector)              :: force
  logical,          allocatable :: atom_found(:)
  
  ! line numbers
  integer :: energy_line
  integer :: forces_start_line
  integer :: forces_end_line
  
  ! temporary variables
  integer :: i,j,ialloc
  
  castep_file = filename
  
  ! Work out line numbers
  energy_line = 0
  forces_start_line = 0
  forces_end_line = 0
  do i=1,size(castep_file)
    line = split(lower_case(castep_file%line(i)))
    ! energy
    if (size(line)>=2) then
      if (line(1)=='final' .and. line(2)=='energy,') then
        energy_line = i
      endif
    endif
    ! forces
    if (size(line)>=2) then
      if ( (line(1)=='***********************' .and. line(2)=='forces') .or. &
         & (line(1)=='*****************' .and. line(2)=='symmetrised')  .or. &
         & (line(1)=='*******************' .and. line(2)=='unconstrained')   &
         & ) then
        forces_start_line = i
      endif
    endif
    if (size(line)==1) then
      if ( forces_start_line/=0 .and. &
         & forces_end_line==0   .and. &
         & len(line(1))>=5) then
        if (slice(line(1),1,5)=='*****') then
          forces_end_line = i
        endif
      endif
    endif
  enddo
  
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
  endif
  
  ! Allocate output
  output = OutputFile(forces_end_line-forces_start_line-7)
  if (output%no_atoms/=structure%no_atoms) then
    call print_line(ERROR//': The number of atoms in the Castep output file &
       &does not match that in the input file.')
    call err()
  endif
  
  ! Read data
  line = split(castep_file%line(energy_line))
  output%energy = dble(line(5)) / EV_PER_HARTREE
  
  allocate(atom_found(output%no_atoms), stat=ialloc); call err(ialloc)
  atom_found = .false.
  do_i : do i=1,output%no_atoms
    line = split(castep_file%line(forces_start_line+5+i))
    species = line(2)
    force = dble(line(4:6)) * ANGSTROM_PER_BOHR / EV_PER_HARTREE
    do j=1,structure%no_atoms
      if (atom_found(j)) then
        cycle
      elseif (structure%atoms(j)%species()==species) then
        atom_found(j) = .true.
        output%species(j) = species
        output%forces(j) = force
        cycle do_i
      endif
    enddo
    
    call print_line(ERROR//': The atoms in the Castep output file do not match &
       &those in the input file.')
    call err()
  enddo do_i
  if (.not. all(atom_found)) then
    call print_line(ERROR//': The atoms in the Castep output file do not match &
       &those in the input file.')
    call err()
  endif
end function

function read_qe_output_file(filename,structure) result(output)
  implicit none
  
  type(String),        intent(in) :: filename
  type(StructureData), intent(in) :: structure
  type(OutputFile)                :: output
  
  ! file contents
  type(IFile)               :: qe_file
  type(String), allocatable :: line(:)
  
  ! line numbers
  integer :: species_start_line
  integer :: species_end_line
  integer :: energy_line
  integer :: forces_start_line
  integer :: forces_end_line
  
  ! qe 'type' to species conversion
  integer                   :: no_species
  type(String), allocatable :: species(:)
  
  ! temporary variables
  integer        :: i
  
  qe_file = filename
  
  ! Work out line numbers
  species_start_line = 0
  do i=1,size(qe_file)
    line = split(lower_case(qe_file%line(i)))
    ! species
    if (line(1)=='atomic' .and. line(2)=='species') then
      species_start_line = i
    elseif ( species_start_line/=0 .and. &
           & species_end_line==0   .and. &
           & size(line)==0) then
      species_end_line = i
    ! energy
    elseif (line(1)=='!') then
      energy_line=i
    ! forces
    elseif (line(1)=='forces' .and. line(2)=='acting') then
      forces_start_line = i
    elseif (line(1)=='total' .and. line(2)=='force') then
      forces_end_line = i
    endif
  enddo
  
  no_species = species_end_line-species_start_line-1
  allocate(species(no_species))
  
  ! Allocate output
  allocate(output%species(forces_end_line-forces_start_line-3))
  allocate(output%forces(forces_end_line-forces_start_line-3))
  
  ! Read data
  do i=1,species_end_line-species_start_line-1
    line = split(qe_file%line(species_start_line+1))
    species(i) = line(1)
  enddo
  
  line = split(qe_file%line(energy_line))
  output%energy = dble(line(5)) * EV_PER_RYDBERG / EV_PER_HARTREE
  
  do i=1,forces_end_line-forces_start_line-3
    line = split(qe_file%line(forces_start_line+1+i))
    output%species(i) = species(int(line(4)))
    output%forces(i) = dble(line(7:9)) * EV_PER_RYDBERG / EV_PER_HARTREE
    
    if (output%species(i)/=structure%atoms(i)%species()) then
      call print_line(ERROR//': The species in the qe output file do not match &
         &those in the input file.')
      call err()
    endif
  enddo
end function

function run_quip_on_file(filename,structure,dir,seedname,symmetry_precision) &
   & result(output)
  implicit none
  
  type(String),        intent(in) :: filename
  type(StructureData), intent(in) :: structure
  type(String),        intent(in) :: dir
  type(String),        intent(in) :: seedname
  real(dp),            intent(in) :: symmetry_precision
  type(OutputFile)                :: output
  
  type(StructureData) :: displaced_structure
  
  type(QuipAtoms) :: quip_atoms
  
  type(String) :: quip_filename
  
  type(QuipResult) :: quip_result
  
  integer :: i,ialloc
  
  ! Read in displaced structure data.
  displaced_structure = input_file_to_StructureData( &
                               & str('castep'),      &
                               & filename,           &
                               & symmetry_precision, &
                               & calculate_symmetry=.false.)
  
  ! Convert structure information into QuipAtoms and QUIP units (eV/Angstrom).
  quip_atoms%lattice = transpose(dble(displaced_structure%lattice)) &
                   & * ANGSTROM_PER_BOHR
  allocate( quip_atoms%atomic_nos(displaced_structure%no_atoms),  &
          & quip_atoms%positions(3,displaced_structure%no_atoms), &
          & stat=ialloc); call err(ialloc)
  do i=1,displaced_structure%no_atoms
    quip_atoms%atomic_nos(i) = int(displaced_structure%atoms(i)%species())
    quip_atoms%positions(:,i) =                                    &
       &   dble(displaced_structure%atoms(i)%cartesian_position()) &
       & * ANGSTROM_PER_BOHR
  enddo
  
  ! Call QUIP.
  quip_filename = dir//'/'//seedname//'_MEAM.xml'
  quip_result = call_quip(quip_atoms, char(quip_filename))
  
  ! Convert QUIP's output into Caesar units (Bohr/Hartree) and types.
  output = OutputFile(displaced_structure%no_atoms)
  output%energy = quip_result%energy / EV_PER_HARTREE
  do i=1,displaced_structure%no_atoms
    output%species(i) = displaced_structure%atoms(i)%species()
    output%forces(i) = quip_result%forces(:,i) &
                   & * ANGSTROM_PER_BOHR       &
                   & / EV_PER_HARTREE
  enddo
end function

function read_output_file(file_type,filename,structure,dir,seedname, &
   & symmetry_precision) result(output)
  implicit none
  
  type(String),        intent(in) :: file_type
  type(String),        intent(in) :: filename
  type(StructureData), intent(in) :: structure
  type(String),        intent(in) :: dir
  type(String),        intent(in) :: seedname
  real(dp),            intent(in) :: symmetry_precision
  type(OutputFile)                :: output
  
  if (file_type=='castep') then
    output = read_castep_output_file(filename,structure)
  elseif (file_type=='qe') then
    output = read_qe_output_file(filename,structure)
  elseif (file_type=='quip') then
    output = run_quip_on_file( filename,  &
                             & structure, &
                             & dir,       &
                             & seedname,  &
                             & symmetry_precision)
  endif
end function
end module
