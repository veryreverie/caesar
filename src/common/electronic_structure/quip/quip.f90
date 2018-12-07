! ======================================================================
! The wrapper for Quip.
! ======================================================================
! N.B. Quip and Caesar have clashing modules, so the two-layer separation
!    of quip_module and quip_wrapper module is necessary.
! quip_wrapper_module should depend on no Caesar modules, and the only
!    exposure quip_module should have to Quip should be through
!    quip_wrapper_module.
module quip_module
  use utils_module
  
  use physical_constants_module
  use basic_structure_module
  use electronic_structure_data_module
  use normal_mode_module
  
  use quip_wrapper_module
  implicit none
  
  private
  
  public :: make_input_filename_xyz
  public :: read_input_file_xyz
  public :: write_input_file_xyz
  public :: run_quip_on_structure
  
  public :: QUIP_LINKED
  
  type :: QuipElectronicStructure
    real(dp)              :: energy
    real(dp), allocatable :: forces(:,:)
    real(dp)              :: virial(3,3)
  end type
  
  interface assignment(=)
    module procedure assign_Atoms_BasicStructure
    module procedure assign_BasicStructure_Atoms
    module procedure assign_ElectronicStructure_QuipElectronicStructure
  end interface
contains

! ----------------------------------------------------------------------
! Functionality involving Quip. Includes:
!    - Reading and writing .xyz files.
!    - Calculating electronic structure.
! ----------------------------------------------------------------------
function make_input_filename_xyz(seedname) result(output)
  implicit none
  
  type(String), intent(in) :: seedname
  type(String)             :: output
  
  output = seedname//'.xyz'
end function

function read_input_file_xyz(filename) result(output)
  implicit none
  
  type(String), intent(in)  :: filename
  type(BasicStructure)      :: output
  
  type(Atoms) :: quip_structure
  
  integer :: ierr
  
  call print_line(ERROR//': Procedure to read .xyz files not yet implemented.')
  call err()
  !call read(quip_structure, char(filename), error=ierr)
  !
  !call print_line('Shape pos: '//shape(quip_structure%pos))
  !call print_line('Shape mass: '//shape(quip_structure%mass))
  !
  !output = quip_structure
end function

subroutine write_input_file_xyz(structure,input_filename,output_filename)
  implicit none
  
  type(BasicStructure), intent(in)           :: structure
  type(String),         intent(in), optional :: input_filename
  type(String),         intent(in)           :: output_filename
  
  type(IFile) :: input_file_in
  type(IFile) :: output_file_in
  type(OFile) :: output_file_out
  
  type(Atoms) :: quip_structure
  
  integer :: i
  
  call print_line(ERROR//': Procedure to write .xyz files not yet &
     &implemented.')
  call err()
  
  !quip_structure = structure
  !call write(quip_structure, char(output_filename))
  !
  !! If an input file is given, replace the second line of the new output file
  !!    with that from the input file.
  !if (present(input_filename)) then
  !  input_file_in = IFile(input_filename)
  !  output_file_in = IFile(output_filename)
  !  output_file_out = OFile(output_filename)
  !  do i=1,size(output_file_in)
  !    if (i==2) then
  !      call output_file_out%print_line(input_file_in%line(i))
  !    else
  !      call output_file_out%print_line(output_file_in%line(i))
  !    endif
  !  enddo
  !endif
end subroutine

function run_quip_on_structure(structure,seedname) result(output)
  implicit none
  
  type(BasicStructure), intent(in) :: structure
  type(String),         intent(in) :: seedname
  type(ElectronicStructure)        :: output
  
  type(String) :: quip_filename
  
  type(Atoms)                   :: quip_structure
  type(QuipElectronicStructure) :: quip_electronic_structure
  
  integer :: ialloc
  
  quip_filename = format_path(seedname//'_MEAM.xml')
  quip_structure = structure
  allocate( quip_electronic_structure%forces(3,size(structure%atoms)), &
          & stat=ialloc); call err(ialloc)
  
  call quip_unified_wrapper(                                   &
     & n                   = size(structure%atoms),            &
     & lattice             = quip_structure%lattice,           &
     & z                   = quip_structure%z,                 &
     & pos                 = quip_structure%pos,               &
     & init_args_str       = 'IP SI_meam',                     &
     & init_args_str_len   = 10,                               &
     & energy              = quip_electronic_structure%energy, &
     & force               = quip_electronic_structure%forces, &
     & virial              = quip_electronic_structure%virial, &
     & do_energy           = .true.,                           &
     & do_force            = .true.,                           &
     & do_virial           = .true.,                           &
     & quip_param_file     = char(quip_filename),              &
     & quip_param_file_len = len(quip_filename),               &
     & calc_args_str       = '',                               &
     & calc_args_str_len   = 0 )
  
  output = quip_electronic_structure
end function

! ----------------------------------------------------------------------
! Conversions between Caesar and Quip types.
! ----------------------------------------------------------------------
subroutine assign_Atoms_BasicStructure(output,input)
  implicit none
  
  type(Atoms),          intent(out) :: output
  type(BasicStructure), intent(in)  :: input
  
  type(String), allocatable :: split_species(:)
  
  integer :: i,ialloc
  
  allocate( output%z(size(input%atoms)),     &
          & output%mass(size(input%atoms)),  &
          & output%pos(3,size(input%atoms)), &
          & stat=ialloc); call err(ialloc)
  
  output%lattice = dble(transpose(input%lattice_matrix)) * ANGSTROM_PER_BOHR
  
  do i=1,size(input%atoms)
    split_species = split_line(input%atoms(i)%species,':')
    output%z(i) = int(split_species(2))
  enddo
  
  output%mass = input%atoms%mass * KG_PER_ME / KG_PER_AMU
  
  do i=1,size(input%atoms)
    output%pos(:,i) = dble(input%atoms(i)%cartesian_position) &
                  & * ANGSTROM_PER_BOHR
  enddo
end subroutine

subroutine assign_BasicStructure_Atoms(output,input)
  implicit none
  
  type(BasicStructure), intent(out) :: output
  type(Atoms),          intent(in)  :: input
  
  type(RealMatrix)              :: lattice
  type(RealVector), allocatable :: positions(:)
  
  integer :: i,ialloc
  
  lattice = transpose(input%lattice) / ANGSTROM_PER_BOHR
  
  allocate(positions(size(input%z)), stat=ialloc); call err(ialloc)
  do i=1,size(positions)
    positions(i) = input%pos(:,i) / ANGSTROM_PER_BOHR
  enddo
  
  output = BasicStructure( lattice,                             &
                         & str(input%z),                        &
                         & input%mass * KG_PER_AMU / KG_PER_ME, &
                         & positions)
end subroutine

subroutine assign_ElectronicStructure_QuipElectronicStructure(output,input)
  implicit none
  
  type(ElectronicStructure),     intent(out) :: output
  type(QuipElectronicStructure), intent(in)  :: input
  
  type(RealVector), allocatable :: forces(:)
  
  integer :: i,ialloc
  
  allocate(forces(size(input%forces,2)), stat=ialloc); call err(ialloc)
  do i=1,size(forces)
    forces(i) = input%forces(:,i) * ANGSTROM_PER_BOHR / EV_PER_HARTREE
  enddo
  
  output = ElectronicStructure( input%energy / EV_PER_HARTREE, &
                              & CartesianForce(forces))
end subroutine
end module