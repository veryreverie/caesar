! ======================================================================
! Reads and writes Quantum Espresso .in files.
! ======================================================================

! N.B. This module has never been tested.

module qe_wrapper_submodule
  use utils_module
  
  use structure_module
  use normal_mode_module
  
  use electronic_structure_data_submodule
  implicit none
  
  private
  
  public :: make_input_filename_qe
  public :: make_output_filename_qe
  public :: read_input_file_qe
  public :: write_input_file_qe
  public :: read_output_file_qe
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

function read_input_file_qe(filename,symmetry_precision,calculate_symmetry) &
   & result(output)
  implicit none
  
  type(String), intent(in)           :: filename
  real(dp),     intent(in)           :: symmetry_precision
  logical,      intent(in), optional :: calculate_symmetry
  type(StructureData)                :: output
  
  call print_line(CODE_ERROR//': Quantum Espresso not yet supported.')
  call err()
end function

subroutine write_input_file_qe(structure,old_qe_in_filename,new_qe_in_filename)
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
    old_qe_in_file = IFile(old_qe_in_filename)
    
    ! --------------------------------------------------
    ! Transform q-points into supercell co-ordinates.
    ! --------------------------------------------------
    do i=1,size(old_qe_in_file)
      line = split_line(lower_case(old_qe_in_file%line(i)))
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
  new_qe_in_file = OFile(new_qe_in_filename)
  call new_qe_in_file%print_line('nat='//structure%no_atoms)
  call new_qe_in_file%print_line('/&end')
  call new_qe_in_file%print_line('CELL_PARAMETERS bohr')
  call new_qe_in_file%print_lines(structure%lattice)
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

function read_output_file_qe(filename,structure) result(output)
  implicit none
  
  type(String),        intent(in) :: filename
  type(StructureData), intent(in) :: structure
  type(ElectronicStructure)       :: output
  
  ! File contents.
  type(IFile)               :: qe_file
  type(String), allocatable :: line(:)
  
  ! Line numbers.
  integer :: species_start_line
  integer :: species_end_line
  integer :: energy_line
  integer :: forces_start_line
  integer :: forces_end_line
  
  ! QE 'type' to species conversion.
  integer                   :: no_species
  type(String), allocatable :: species(:)
  
  ! Output variables.
  real(dp)                      :: energy
  type(RealVector), allocatable :: forces(:)
  
  ! Temporary variables.
  integer :: i,ialloc
  
  qe_file = IFile(filename)
  
  ! Work out line numbers.
  species_start_line = 0
  do i=1,size(qe_file)
    line = split_line(lower_case(qe_file%line(i)))
    ! Species.
    if (line(1)=='atomic' .and. line(2)=='species') then
      species_start_line = i
    elseif ( species_start_line/=0 .and. &
           & species_end_line==0   .and. &
           & size(line)==0) then
      species_end_line = i
    ! Energy.
    elseif (line(1)=='!') then
      energy_line=i
    ! Forces.
    elseif (line(1)=='forces' .and. line(2)=='acting') then
      forces_start_line = i
    elseif (line(1)=='total' .and. line(2)=='force') then
      forces_end_line = i
    endif
  enddo
  
  ! Check line numbers.
  if (forces_end_line-forces_start_line-3/=structure%no_atoms) then
    call print_line(ERROR//': The number of atoms in the Castep output file &
       &does not match that in the input file.')
    call err()
  endif
  
  ! Allocate arrays.
  no_species = species_end_line-species_start_line-1
  allocate( species(no_species),        &
          & forces(structure%no_atoms), &
          & stat=ialloc); call err(ialloc)
  
  ! Read data.
  do i=1,species_end_line-species_start_line-1
    line = split_line(qe_file%line(species_start_line+1))
    species(i) = line(1)
  enddo
  
  line = split_line(qe_file%line(energy_line))
  energy = dble(line(5)) * EV_PER_RYDBERG / EV_PER_HARTREE
  
  do i=1,forces_end_line-forces_start_line-3
    line = split_line(qe_file%line(forces_start_line+1+i))
    
    if (species(int(line(4)))/=structure%atoms(i)%species()) then
      call print_line(ERROR//': The species in the qe output file do not match &
         &those in the input file.')
      call err()
    endif
    
    forces(i) = dble(line(7:9)) * EV_PER_RYDBERG / EV_PER_HARTREE
  enddo
  
  ! Construct output.
  output = ElectronicStructure(energy,CartesianForce(forces))
end function
end module
