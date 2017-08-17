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
  use potential_module
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
  type(StructureData)                 :: structure
  type(String),           allocatable :: no_supercells_file(:)
  integer                             :: no_supercells
  type(StructureData),    allocatable :: supercells(:)
  type(StructureData),    allocatable :: supercell
  type(QpointData),       allocatable :: qpoints(:)
  type(NormalMode),       allocatable :: modes(:)
  type(CoupledModes),     allocatable :: coupling(:)
  integer                             :: no_sampling_points
  type(CouplingSampling), allocatable :: sampling(:)
  
  ! DFT results (Indexed as sampling_points).
  type(String)                  :: dft_output_filename
  type(DftOutputFile)           :: dft_output_file
  
  ! Harmonic eigenstates coupling between them.
  type(HarmonicState), allocatable :: harmonic_states(:,:)
  
  ! The Born-Oppenheimer potential.
  type(PolynomialPotential) :: potential
  
  ! SCF variables
  integer                           :: scf_step
  type(RealEigenstuff), allocatable :: vscf_eigenstuff(:)
  real(dp),             allocatable :: energy_change(:,:)
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
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
    
    ! Calculate harmonic eigenstates, {|a>}, along each normal mode.
    allocate( harmonic_states(no_harmonic_states,structure%no_modes),    &
            & stat=ialloc); call err(ialloc)
    do j=1,structure%no_modes
      harmonic_states(:,j) = generate_harmonic_basis( modes(j)%frequency, &
                                                    & no_harmonic_states)
    enddo
    
    ! Read in coupling and sampling points.
    coupling = read_coupling_file(wd//'/qpoint_'//i//'/coupling.dat')
    allocate(sampling(size(coupling)), stat=ialloc); call err(ialloc)
    do j=1,size(coupling)
      sampling(j)%coupling = coupling(j)
      sampling(j)%sampling_points = read_sampling_points_file( &
         & wd//'/qpoint_'//i//'/coupling_'//j//'/sampling_points.dat')
      no_sampling_points = size(sampling(j)%sampling_points)
      allocate( sampling(j)%energy(no_sampling_points),                    &
              & sampling(j)%forces(supercell%no_atoms,no_sampling_points), &
              & stat=ialloc); call err(ialloc)
      
      ! Read in DFT outputs.
      do k=1,no_sampling_points
        dft_output_file = read_dft_output_file(dft_code, wd//                 &
                                                    & '/qpoint_'//i//         &
                                                    & '/coupling_'//j//       &
                                                    & '/sampling_point_'//k// &
                                                    & '/'//dft_output_filename)
        sampling(j)%energy(k) = dft_output_file%energy
        sampling(j)%forces(:,k) = dft_output_file%forces
      enddo
    enddo
    
    ! Calculate the Born-Oppenheimer potential.
    potential = calculate_potential( no_basis_functions, &
                                   & sampling,           &
                                   & modes,              &
                                   & harmonic_states)
    
    ! Initialise eigenstuff to harmonic values.
    allocate( vscf_eigenstuff(structure%no_modes), &
            & stat=ialloc); call err(ialloc)
    do j=1,structure%no_modes
      vscf_eigenstuff(j)%evecs = real(int(identity(no_basis_functions)))
      do k=1,no_basis_functions
        vscf_eigenstuff(j)%evals(k) = (k-0.5_dp)*modes(j)%frequency
      enddo
    enddo
    
    ! Initialise energy change array.
    allocate( energy_change(no_basis_functions, structure%no_modes), &
            & stat=ialloc); call err(ialloc)
    
    ! SCF cycles.
    do scf_step=1,max_scf_cycles
      ! Record old eigenvalues.
      do j=1,structure%no_modes
        energy_change(:,j) = vscf_eigenstuff(j)%evals
      enddo
      
      ! Run SCF calculations.
      vscf_eigenstuff = scf(potential,vscf_eigenstuff)
      
      ! Calculate L2-norm change in eigenvalues.
      do j=1,structure%no_modes
        energy_change(:,j) = l2_norm(vec( vscf_eigenstuff(j)%evals &
                                        & - energy_change(:,j) ))
      enddo
      
      ! Check for convergence
      if (sum(energy_change)<scf_convergence_threshold) then
        exit
      elseif (scf_step==max_scf_cycles) then
        call err()
      endif
    enddo
  enddo
end subroutine
end module
