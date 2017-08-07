! ======================================================================
! Uses DFT results to calculate anharmonic properties.
! ======================================================================
module calculate_anharmonic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function calculate_anharmonic_keywords() result(keywords)
  use help_module
  implicit none
  
  type(KeywordData) :: keywords(1)
  
  keywords = [                                                               &
  & make_keyword( 'no_basis_functions',                                      &
  &               'no_basis_functions is the number of basis functions along &
  &each normal mode.') ]
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_anharmonic(arguments)
  use dictionary_module
  use structure_module
  use qpoints_module
  use coupling_module
  use sampling_points_module
  use linear_algebra_module
  use dft_output_file_module
  use basis_function_module
  use cubic_grid_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Inputs.
  type(String) :: wd
  integer      :: no_basis_functions
  
  ! Previous inputs.
  type(Dictionary) :: setup_harmonic_arguments
  type(Dictionary) :: setup_anharmonic_arguments
  type(String)     :: harmonic_path
  type(String)     :: dft_code
  type(String)     :: seedname
  type(String)     :: grid_type
  
  ! Previously calculated data.
  type(StructureData)              :: structure
  type(String),        allocatable :: no_supercells_file(:)
  integer                          :: no_supercells
  type(StructureData), allocatable :: supercells(:)
  type(StructureData), allocatable :: supercell
  type(QpointData),    allocatable :: qpoints(:)
  type(CoupledModes),  allocatable :: coupling(:)
  type(SamplingPoint), allocatable :: sample_points(:)
  
  ! DFT results (Indexed as sample_points).
  type(String)                  :: dft_output_filename
  type(DftOutputFile)           :: dft_output_file
  real(dp),         allocatable :: energy(:)
  type(RealVector), allocatable :: forces(:,:)
  
  ! Basis functions and coupling between them.
  type(BasisFunction), allocatable :: basis_functions(:)
  real(dp),            allocatable :: harmonic_couplings(:,:)
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  ! Read in inputs.
  wd = arguments%value('working_directory')
  no_basis_functions = int(arguments%value('no_basis_functions'))
  
  ! Read in setup_anharmonic settings.
  call setup_anharmonic_arguments%read_file( &
     & wd//'/setup_anharmonic.used_settings')
  harmonic_path = setup_anharmonic_arguments%value('harmonic_path')
  grid_type = setup_anharmonic_arguments%value('grid_type')
  
  ! Read in setup_harmonic settings.
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  dft_code = setup_harmonic_arguments%value('dft_code')
  seedname = setup_harmonic_arguments%value('seedname')
  
  ! Read in structure and supercells.
  structure = read_structure_file(harmonic_path//'/structure.dat')
  no_supercells_file = read_lines(harmonic_path//'/no_supercells.dat')
  no_supercells = int(no_supercells_file(1))
  allocate(supercells(no_supercells), stat=ialloc); call err(ialloc)
  do i=1,no_supercells
    supercells(i) = read_structure_file( &
       & harmonic_path//'/Supercell_'//i//'/structure.dat')
  enddo
  
  ! Read in q-points.
  qpoints = read_qpoints_file(harmonic_path//'/qpoints_ibz.dat')
  
  ! Construct DFT output filename.
  dft_output_filename = make_dft_output_filename(dft_code,seedname)
  
  ! Loop over q-points, calculating anharmonic correction at each point.
  do i=1,size(qpoints)
    ! Identify the relevant supercell.
    supercell = supercells(qpoints(i)%sc_id)
    
    ! Read in coupling and sampling points.
    coupling = read_coupling_file(wd//'/qpoint_'//i//'/coupling.dat')
    sample_points = read_sampling_points_file( &
       & wd//'/qpoint_'//i//'/sampling_points.dat')
    
    ! Read in DFT outputs.
    allocate( energy(size(sample_points)),                     &
            & forces(supercell%no_atoms,size(sample_points)),  &
            & stat=ialloc); call err(ialloc)
    do j=1,size(sample_points)
      dft_output_file = read_dft_output_file(dft_code, &
         & wd//'/qpoint_'//i//'/sampling_point_'//j//'/'//dft_output_filename)
      energy(j) = dft_output_file%energy
      forces(:,j) = dft_output_file%forces
    enddo
    
    ! calculate harmonic basis functions.
    basis_functions = generate_basis_functions()
    
    ! Calculate coupling between harmonic basis functions.
    if (grid_type=='cubic') then
      harmonic_couplings = harmonic_couplings_cubic_grid(basis_functions, &
         & coupling,sample_points,energy,forces)
    endif
  enddo
end subroutine
end module
