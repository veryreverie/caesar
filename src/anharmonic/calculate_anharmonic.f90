! ======================================================================
! Calculates anharmonic properties, using the results of run_anharmonic.
! ======================================================================
module calculate_anharmonic_module
  use common_module
  
  use setup_harmonic_module
  use calculate_normal_modes_module
  
  use anharmonic_common_module
  use polynomial_module
  
  use setup_anharmonic_module
  implicit none
  
  private
  
  public :: calculate_anharmonic_keywords
  public :: calculate_anharmonic_mode
  public :: calculate_anharmonic
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function calculate_anharmonic_keywords() result(keywords)
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                             &
  & KeywordData( 'energy_to_force_ratio',                                  &
  &              'energy_to_force_ratio is the ratio of how penalised &
  &deviations in energy are compared to deviations in forces when the &
  &potential is being fitted. This should be given in units of (Hartree &
  &per primitive cell) divided by (Hartree per bohr). Due to the use of &
  &mass-weighted co-ordinates, in systems containing different elements &
  &forces along modes with higher contributions from heavier elements will &
  &be weighted less than this.')]
end function

function calculate_anharmonic_mode() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'calculate_anharmonic'
  output%description = 'Uses the results of run_anharmonic to calculate &
     &anharmonic properties. Should be run after run_anharmonic.'
  output%keywords = calculate_anharmonic_keywords()
  output%main_subroutine => calculate_anharmonic
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_anharmonic(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! Input arguments.
  real(dp) :: energy_to_force_ratio
  real(dp) :: weighted_energy_force_ratio
  
  ! Arguments to setup_anharmonic.
  type(Dictionary) :: setup_anharmonic_arguments
  type(String)     :: harmonic_path
  type(String)     :: potential_representation
  integer          :: potential_expansion_order
  logical          :: vscf_basis_functions_only
  real(dp)         :: maximum_displacement
  real(dp)         :: frequency_of_max_displacement
  
  ! Arguments to setup_harmonic.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: file_type
  type(String)     :: seedname
  real(dp)         :: symmetry_precision
  
  ! Arguments to calculate_normal_modes.
  type(Dictionary) :: calculate_normal_modes_arguments
  type(String)     :: calculation_type
  
  ! Primitive structure.
  type(StructureData) :: structure
  
  ! Large anharmonic supercell and its q-points.
  type(StructureData)           :: anharmonic_supercell
  type(QpointData), allocatable :: qpoints(:)
  
  ! Normal modes.
  type(ComplexMode),        allocatable :: complex_modes(:)
  type(RealMode),           allocatable :: real_modes(:)
  
  ! maximum_displacement in mass-weighted co-ordinates.
  real(dp) :: maximum_weighted_displacement
  
  ! Degeneracy and coupling data.
  type(DegenerateSubspace), allocatable :: degenerate_subspaces(:)
  type(DegenerateSymmetry), allocatable :: degenerate_symmetries(:)
  type(SubspaceCoupling),   allocatable :: subspace_coupling(:)
  
  ! Anharmonic data container.
  type(AnharmonicData) :: anharmonic_data
  
  ! Data specific to the chosen representation of the potential.
  type(PotentialPointer) :: potential
  
  ! Files and directories.
  type(OFile)               :: logfile
  type(IFile)               :: qpoints_file
  type(IFile)               :: complex_modes_file
  type(IFile)               :: real_modes_file
  type(IFile)               :: subspaces_file
  type(IFile)               :: subspace_coupling_file
  type(IFile)               :: calculation_directories_file
  type(String), allocatable :: calculation_directories(:)
  type(String)              :: sampling_points_dir
  
  ! Temporary variables.
  integer :: i,ialloc
  
  ! ----------------------------------------------------------------------
  ! Read in inputs and previously calculated data which is
  !    independent of the choice of potential representation.
  ! ----------------------------------------------------------------------
  
  wd = arguments%value('working_directory')
  energy_to_force_ratio = dble(arguments%value('energy_to_force_ratio'))
  
  ! Read in setup_anharmonic arguments.
  setup_anharmonic_arguments = Dictionary(setup_anharmonic_keywords())
  call setup_anharmonic_arguments%read_file( &
     & wd//'/setup_anharmonic.used_settings')
  harmonic_path = setup_anharmonic_arguments%value('harmonic_path')
  potential_representation = &
     & setup_anharmonic_arguments%value('potential_representation')
  potential_expansion_order = &
     & int(setup_anharmonic_arguments%value('potential_expansion_order'))
  vscf_basis_functions_only = &
     & lgcl(setup_anharmonic_arguments%value('vscf_basis_functions_only'))
  maximum_displacement = &
     & dble(setup_anharmonic_arguments%value('maximum_displacement'))
  frequency_of_max_displacement = &
     & dble(setup_anharmonic_arguments%value('frequency_of_max_displacement'))
  
  ! Read in setup_harmonic arguments.
  setup_harmonic_arguments = Dictionary(setup_harmonic_keywords())
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  file_type = setup_harmonic_arguments%value('file_type')
  seedname = setup_harmonic_arguments%value('seedname')
  symmetry_precision = &
     & dble(setup_harmonic_arguments%value('symmetry_precision'))
  
  ! Read in calculate_normal_modes arguments.
  calculate_normal_modes_arguments = &
     & Dictionary(calculate_normal_modes_keywords())
  call calculate_normal_modes_arguments%read_file( &
     & harmonic_path//'/calculate_normal_modes.used_settings')
  calculation_type = calculate_normal_modes_arguments%value('calculation_type')
  
  ! Read in structure.
  structure = read_structure_file( harmonic_path//'/structure.dat', &
                                 & symmetry_precision)
  
  ! Read in large anharmonic supercell and its q-points.
  anharmonic_supercell = read_structure_file( &
           & wd//'/anharmonic_supercell.dat', &
           & symmetry_precision,              &
           & calculate_symmetry=.false.)
  
  qpoints_file = IFile(wd//'/qpoints.dat')
  qpoints = QpointData(qpoints_file%sections())
  
  ! Read in normal modes.
  complex_modes_file = IFile(wd//'/complex_modes.dat')
  complex_modes = ComplexMode(complex_modes_file%sections())
  
  real_modes_file = IFile(wd//'/real_modes.dat')
  real_modes = RealMode(real_modes_file%sections())
  
  ! Read in degenerate subspaces and subspace coupling.
  subspaces_file = IFile(wd//'/degenerate_subspaces.dat')
  degenerate_subspaces = DegenerateSubspace(subspaces_file%sections())
  
  subspace_coupling_file = IFile(wd//'/subspace_coupling.dat')
  subspace_coupling = SubspaceCoupling(subspace_coupling_file%lines())
  
  ! ----------------------------------------------------------------------
  ! Re-calculate symmetries in degenerate subspaces.
  ! ----------------------------------------------------------------------
  ! Open logfile.
  logfile = OFile(wd//'/calculate_anharmonic_logfile.dat')
  
  allocate( degenerate_symmetries(size(structure%symmetries)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(structure%symmetries)
    degenerate_symmetries(i) = DegenerateSymmetry( structure%symmetries(i), &
                                                 & degenerate_subspaces,    &
                                                 & complex_modes,           &
                                                 & qpoints,                 &
                                                 & logfile)
  enddo
  
  ! ----------------------------------------------------------------------
  ! Load anharmonic data into container.
  ! ----------------------------------------------------------------------
  anharmonic_data = AnharmonicData( structure,                     &
                                  & symmetry_precision,            &
                                  & anharmonic_supercell,          &
                                  & qpoints,                       &
                                  & complex_modes,                 &
                                  & real_modes,                    &
                                  & degenerate_subspaces,          &
                                  & degenerate_symmetries,         &
                                  & subspace_coupling,             &
                                  & vscf_basis_functions_only,     &
                                  & maximum_weighted_displacement, &
                                  & frequency_of_max_displacement )
  
  ! ----------------------------------------------------------------------
  ! Read in calculation directories.
  ! ----------------------------------------------------------------------
  calculation_directories_file = IFile(wd//'/calculation_directories.dat')
  calculation_directories = calculation_directories_file%lines()
  
  ! ----------------------------------------------------------------------
  ! Run representation-specific code.
  ! ----------------------------------------------------------------------
  
  ! Calculate weighted energy to force ratio.
  weighted_energy_force_ratio = energy_to_force_ratio &
                            & * sqrt(maxval(structure%atoms%mass()))
  
  ! Initialise potential to the chosen representation
  if (potential_representation=='polynomial') then
    potential = PolynomialPotential(potential_expansion_order)
  else
    call print_line( ERROR//': Unrecognised potential representation: '// &
                   & potential_representation)
    call err()
  endif
  
  ! Call potential-specific function.
  sampling_points_dir = wd//'/sampling_points'
  call potential%generate_potential( anharmonic_data,             &
                                   & weighted_energy_force_ratio, &
                                   & sampling_points_dir,         &
                                   & logfile,                     &
                                   & read_calculation_directory)
contains
  ! Lambda to read electronic structure results from a directory.
  ! Captures:
  !    - calculation_directories
  !    - symmetry_precision
  !    - calculation_type
  !    - file_type
  !    - seedname
  !    - wd
  function read_calculation_directory(directory) result(output)
    implicit none
    
    type(String), intent(in)  :: directory
    type(ElectronicStructure) :: output
    
    type(StructureData) :: structure
    type(String)        :: filename
    
    ! Check that directory is valid.
    if (.not. any(directory==calculation_directories)) then
      call print_line(CODE_ERROR//': Trying to read electronic structure from &
         &a directory where such calculations were not run.')
      call print_line('Directory: '//directory)
      call err()
    endif
    
    ! Read in structure.
    structure = read_structure_file( directory//'/structure.dat', &
                                   & symmetry_precision,          &
                                   & calculate_symmetry=.false.)
    
    ! Construct filename.
    if (calculation_type=='script') then
      filename = make_output_filename(file_type,seedname)
    else
      filename = make_input_filename(file_type,seedname)
    endif
    filename = directory//'/'//filename
    
    ! Read in electronic structure.
    output = read_output_file( file_type,          &
                             & filename,           &
                             & structure,          &
                             & wd,                 &
                             & seedname,           &
                             & symmetry_precision, &
                             & calculation_type)
  end function
end subroutine
end module
