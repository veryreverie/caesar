! ======================================================================
! Interfaces to the QUIP library.
! ======================================================================
module quip_wrapper_submodule
  use utils_module
  
  use structure_module
  
  use electronic_structure_data_submodule
  
  ! Use modules from QUIP itself.
  use quip_unified_wrapper_module, only : quip_unified_wrapper,initialise
  use libatoms_module, only : Atoms,read,write
  implicit none
  
  private
  
  public :: QuipAtoms
  public :: QuipResult
  public :: call_quip
  public :: read_input_file_xyz
  public :: write_input_file_xyz
  public :: run_quip_on_file
  
  type :: QuipResult
    real(dp)              :: energy
    real(dp), allocatable :: forces(:,:)
    real(dp)              :: virial(3,3)
  end type
  
  type :: QuipAtoms
    real(dp)              :: lattice(3,3)
    integer,  allocatable :: atomic_nos(:)
    real(dp), allocatable :: masses(:)
    real(dp), allocatable :: positions(:,:)
  contains
    procedure, private :: check => check_QuipAtoms
  end type
contains

subroutine check_QuipAtoms(this)
  implicit none
  
  class(QuipAtoms), intent(in) :: this
  
  if (.not. allocated(this%atomic_nos)) then
    write(*,*) 'QuipAtoms error: atomic_nos not allocated.'
    stop
  elseif (.not. allocated(this%positions)) then
    write(*,*) 'QuipAtoms error: positions not allocated.'
    stop
  elseif (size(this%positions,1)/=3) then
    write(*,*) 'QuipAtoms error: size(positions,1)/=3'
    stop
  elseif (size(this%positions,2)/=size(this%atomic_nos)) then
    write(*,*) 'QuipAtoms error: size(positions,2)/=size(atomic_nos)'
    stop
  elseif (size(this%atomic_nos)/=size(this%masses)) then
    write(*,*) 'QuipAtoms error: size(atomic_nos)/=size(masses)'
    stop
  endif
end subroutine

function call_quip(quip_atoms,quip_filename) result(output)
  implicit none
  
  type(QuipAtoms), intent(in) :: quip_atoms
  character(*),    intent(in) :: quip_filename
  type(QuipResult)            :: output
  
  integer :: ialloc
  
  call quip_atoms%check()
  
  allocate( output%forces(3,size(quip_atoms%atomic_nos)), &
          & stat=ialloc); call err(ialloc)
  
  call quip_unified_wrapper(                              &
     & n                   = size(quip_atoms%atomic_nos), &
     & lattice             = quip_atoms%lattice,          &
     & z                   = quip_atoms%atomic_nos,       &
     & pos                 = quip_atoms%positions,        &
     & init_args_str       = 'IP SI_meam',                &
     & init_args_str_len   = 10,                          &
     & energy              = output%energy,               &
     & force               = output%forces,               &
     & virial              = output%virial,               &
     & do_energy           = .true.,                      &
     & do_force            = .true.,                      &
     & do_virial           = .true.,                      &
     & quip_param_file     = quip_filename,               &
     & quip_param_file_len = len(quip_filename),          &
     & calc_args_str       = '',                          &
     & calc_args_str_len   = 0 )
end function

function read_input_file_xyz(filename) result(output)
  implicit none
  
  type(String), intent(in)  :: filename
  type(BasicStructure)      :: output
  
  type(Atoms)                   :: quip_atoms
  type(RealVector), allocatable :: positions(:)
  
  integer :: i,ialloc,ierr
  
  call read(quip_atoms, char(filename), error=ierr)
  
  allocate(positions(size(quip_atoms%pos)), stat=ialloc); call err(ialloc)
  do i=1,size(positions)
    positions(i) = quip_atoms%pos(:,i) / ANGSTROM_PER_BOHR
  enddo
  
  output = BasicStructure( mat(quip_atoms%lattice) / ANGSTROM_PER_BOHR, &
                         & str(quip_atoms%z),                           &
                         & quip_atoms%mass,                             &
                         & positions)
end function

subroutine write_input_file_xyz(structure,filename)
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(String),        intent(in) :: filename
  
  type(Atoms) :: quip_atoms
  
  integer :: i,ialloc
  
  quip_atoms%lattice = dble(structure%lattice) * ANGSTROM_PER_BOHR
  quip_atoms%z       = int(structure%atoms%species())
  quip_atoms%mass    = structure%atoms%mass()
  allocate( quip_atoms%pos(3,size(structure%atoms)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(structure%atoms)
  quip_atoms%pos(:,i) = dble(structure%atoms(i)%cartesian_position()) &
                    & * ANGSTROM_PER_BOHR
  enddo
  
  call write(quip_atoms, char(filename))
end subroutine

function run_quip_on_file(displaced_structure,file_type,filename,structure,dir,seedname, &
   & symmetry_precision) result(output)
  implicit none
  
  type(StructureData), intent(in) :: displaced_structure
  type(String),        intent(in) :: file_type
  type(String),        intent(in) :: filename
  type(StructureData), intent(in) :: structure
  type(String),        intent(in) :: dir
  type(String),        intent(in) :: seedname
  real(dp),            intent(in) :: symmetry_precision
  type(ElectronicStructure)       :: output
  
  type(QuipAtoms) :: quip_atoms
  
  type(String) :: quip_filename
  
  type(QuipResult) :: quip_result
  
  ! Output variables.
  real(dp)                      :: energy
  type(RealVector), allocatable :: forces(:)
  
  integer :: i,ialloc
  
  ! Convert structure information into QuipAtoms and QUIP units (eV/Angstrom).
  quip_atoms%lattice = transpose(dble(displaced_structure%lattice)) &
                   & * ANGSTROM_PER_BOHR
  allocate( quip_atoms%atomic_nos(displaced_structure%no_atoms),  &
          & quip_atoms%positions(3,displaced_structure%no_atoms), &
          & stat=ialloc); call err(ialloc)
  quip_atoms%atomic_nos = int(displaced_structure%atoms%species())
  quip_atoms%masses = displaced_structure%atoms%mass()
  do i=1,displaced_structure%no_atoms
    quip_atoms%positions(:,i) =                                    &
       &   dble(displaced_structure%atoms(i)%cartesian_position()) &
       & * ANGSTROM_PER_BOHR
  enddo
  
  ! Call QUIP.
  quip_filename = dir//'/'//seedname//'_MEAM.xml'
  quip_result = call_quip(quip_atoms, char(quip_filename))
  
  ! Convert QUIP's output into Caesar units (Bohr/Hartree) and types.
  energy = quip_result%energy / EV_PER_HARTREE
  allocate(forces(displaced_structure%no_atoms), stat=ialloc); call err(ialloc)
  do i=1,displaced_structure%no_atoms
    forces(i) = quip_result%forces(:,i) &
            & * ANGSTROM_PER_BOHR       &
            & / EV_PER_HARTREE
  enddo
  
  output = ElectronicStructure(energy,forces)
end function
end module
