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
  
  public :: run_quip_on_structure
  
  public :: QUIP_LINKED
  
  type :: QuipElectronicStructure
    real(dp)              :: energy
    real(dp), allocatable :: forces(:,:)
    real(dp), allocatable :: stress(:,:)
  end type
  
  interface assignment(=)
    module procedure assign_Atoms_BasicStructure
    module procedure assign_BasicStructure_Atoms
    module procedure assign_ElectronicStructure_QuipElectronicStructure
  end interface
contains

! ----------------------------------------------------------------------
! Functionality involving Quip. Includes:
!    - Calculating electronic structure.
! ----------------------------------------------------------------------
function run_quip_on_structure(structure,seedname,use_forces,use_hessians, &
   & calculate_stress) result(output)
  implicit none
  
  type(BasicStructure), intent(in) :: structure
  type(String),         intent(in) :: seedname
  logical,              intent(in) :: use_forces
  logical,              intent(in) :: use_hessians
  logical,              intent(in) :: calculate_stress
  type(ElectronicStructure)        :: output
  
  type(String) :: quip_filename
  
  type(Atoms)                   :: quip_structure
  type(QuipElectronicStructure) :: quip_electronic_structure
  
  integer :: ialloc
  
  if (use_hessians) then
    call print_line(ERROR//': Reading Hessians from Quip is not implemented.')
    call quit()
  endif
  
  quip_filename = format_path(seedname//'_MEAM.xml')
  quip_structure = structure
  allocate( quip_electronic_structure%forces(3,size(structure%atoms)), &
          & quip_electronic_structure%stress(3,3),                     &
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
     & virial              = quip_electronic_structure%stress, &
     & do_energy           = .true.,                           &
     & do_force            = use_forces,                       &
     & do_virial           = calculate_stress,                 &
     & quip_param_file     = char(quip_filename),              &
     & quip_param_file_len = len(quip_filename),               &
     & calc_args_str       = '',                               &
     & calc_args_str_len   = 0                                 )
  
  if (.not. use_forces) then
    deallocate(quip_electronic_structure%forces, stat=ialloc); call err(ialloc)
  endif
  
  if (.not. calculate_stress) then
    deallocate(quip_electronic_structure%stress, stat=ialloc); call err(ialloc)
  endif
  
  ! Quip calculates stress*volume rather than stress.
  if (calculate_stress) then
    quip_electronic_structure%stress = quip_electronic_structure%stress &
                                   & / structure%volume()
  endif
  
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
  
  lattice = mat(transpose(input%lattice) / ANGSTROM_PER_BOHR)
  
  allocate(positions(size(input%z)), stat=ialloc); call err(ialloc)
  do i=1,size(positions)
    positions(i) = vec(input%pos(:,i) / ANGSTROM_PER_BOHR)
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
  
  real(dp)                          :: energy
  type(CartesianForce), allocatable :: forces
  type(RealMatrix),     allocatable :: stress
  
  integer :: i,ialloc
  
  energy = input%energy / EV_PER_HARTREE
  
  if (allocated(input%forces)) then
    forces = CartesianForce([(                                    &
       & vec(input%forces(:,i))*ANGSTROM_PER_BOHR/EV_PER_HARTREE, &
       & i=1,                                                     &
       & size(input%forces,2)                                     )])
  endif
  
  if (allocated(input%stress)) then
    stress = mat(input%stress)/EV_PER_HARTREE
  endif
  
  output = ElectronicStructure(energy=energy, forces=forces, stress=stress)
end subroutine
end module
