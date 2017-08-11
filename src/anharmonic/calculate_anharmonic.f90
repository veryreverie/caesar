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
  
  type(KeywordData) :: keywords(4)
  
  keywords = [                                                                &
  & make_keyword( 'no_harmonic_states',                                       &
  &               'no_harmonic_states is the number of harmonic eigenstates &
  &in the direction of each normal mode.'),                                   &
  & make_keyword( 'potential_expansion_order',                                &
  &               'potential_expansion_order is the order up to which the &
  &potential is expanded. e.g. a cubic expansion would be order 3.'),         &
  & make_keyword( 'scf_convergence_threshold',                                &
  &               'scf_convergence_threshold is the energy to within which &
  &the VSCF calculation will be converged.'),                                 &
  & make_keyword( 'max_scf_cycles',                                           &
  &               'max_scf_cycles is the maximum number of SCF cycles which &
  &will be carried out as part of the VSCF calculation.') ]
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_anharmonic(arguments)
  use dictionary_module
  use structure_module
  use qpoints_module
  use normal_mode_module
  use coupling_module
  use sampling_points_module
  use linear_algebra_module
  use dft_output_file_module
  use harmonic_states_module
  use potential_basis_module
  use cubic_grid_module
  use scf_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Inputs.
  type(String) :: wd
  integer      :: no_harmonic_states
  integer      :: no_basis_functions
  real(dp)     :: scf_convergence_threshold
  integer      :: max_scf_cycles
  
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
  type(NormalMode),    allocatable :: modes(:)
  type(CoupledModes),  allocatable :: coupling(:)
  type(SamplingPoint), allocatable :: sampling_points(:)
  
  ! DFT results (Indexed as sampling_points).
  type(String)                  :: dft_output_filename
  type(DftOutputFile)           :: dft_output_file
  real(dp),         allocatable :: energy(:)
  type(RealVector), allocatable :: forces(:,:)
  
  ! Harmonic eigenstates coupling between them.
  type(HarmonicState), allocatable :: harmonic_states(:)
  real(dp),            allocatable :: gaussian_integrals(:)
  type(RealMatrix),    allocatable :: harmonic_couplings(:,:)
  
  ! The Born-Oppenheimer potential.
  type(Monomial), allocatable :: potential(:)
  
  ! SCF variables
  integer                       :: scf_step
  type(RealMatrix), allocatable :: eigenstates(:)
  type(ScfResult)               :: scf_result
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  ! Read in inputs.
  wd = arguments%value('working_directory')
  no_harmonic_states = int(arguments%value('no_harmonic_states'))
  no_basis_functions = int(arguments%value('potential_expansion_order'))
  scf_convergence_threshold = &
     & dble(arguments%value('scf_convergence_threshold'))
  max_scf_cycles = int(arguments%value('max_scf_cycles'))
  
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
    
    ! Read in normal modes.
    allocate(modes(structure%no_modes), stat=ialloc); call err(ialloc)
    do j=1,structure%no_modes
      modes(j) = read_normal_mode_file( &
         & harmonic_path//'/qpoint_'//i//'/mode_'//j//'.dat')
    enddo
    
    ! Read in coupling and sampling points.
    coupling = read_coupling_file(wd//'/qpoint_'//i//'/coupling.dat')
    sampling_points = read_sampling_points_file( &
       & wd//'/qpoint_'//i//'/sampling_points.dat')
    
    ! Calculate harmonic states and coupling between them.
    allocate( harmonic_states(no_harmonic_states),                       &
            & harmonic_couplings(no_basis_functions,structure%no_modes), &
            & stat=ialloc); call err(ialloc)
    do j=1,structure%no_modes
      ! calculate harmonic basis functions, {|a>}.
      harmonic_states = generate_harmonic_basis( modes(j)%frequency, &
                                               & no_harmonic_states)
      ! Calculate gaussian integrals, I(n) = integral[u^n e^(-freq*u*u)].
      gaussian_integrals = calculate_gaussian_integrals( modes(j)%frequency, &
                                                       & no_harmonic_states, &
                                                       & no_basis_functions)
      ! Calculate coupling between harmonic states and monomial potentials,
      !    {<a|u^n|b>}.
      harmonic_couplings(:,j) = generate_harmonic_couplings( harmonic_states, &
                                                      & no_basis_functions,   &
                                                      & gaussian_integrals)
    enddo
    
    ! Calculate potential basis functions.
    ! TODO
    
    ! Read in DFT outputs.
    allocate( energy(size(sampling_points)),                     &
            & forces(supercell%no_atoms,size(sampling_points)),  &
            & stat=ialloc); call err(ialloc)
    do j=1,size(sampling_points)
      dft_output_file = read_dft_output_file(dft_code, &
         & wd//'/qpoint_'//i//'/sampling_point_'//j//'/'//dft_output_filename)
      energy(j) = dft_output_file%energy
      forces(:,j) = dft_output_file%forces
    enddo
    
    ! Calculate potential basis function coefficients using regression.
    ! TODO
    
    ! Convert the potential from basis functions into monomials.
    ! TODO
    
    ! Simplify potential, summing across duplicate monomials.
    potential = simplify(potential)
    
    ! Initialise eigenstates.
    ! TODO
    
    ! SCF cycles.
    do scf_step=1,max_scf_cycles
      ! Run SCF calculations.
      scf_result = scf(harmonic_couplings,potential,eigenstates)
      eigenstates = scf_result%eigenstates
      
      ! Check for convergence
      if (scf_result%energy_change<scf_convergence_threshold) then
        exit
      elseif (scf_step==max_scf_cycles) then
        call err()
      endif
    enddo
  enddo
end subroutine
end module
