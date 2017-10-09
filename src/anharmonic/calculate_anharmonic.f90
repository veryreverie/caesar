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
  use keyword_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                                &
  & make_keyword( 'energy_error',                                             &
  &               'energy_error is the expected error in the calculated &
  &energies.'),                                                               &
  & make_keyword( 'force_error',                                              &
  &               'force_error is the expected error in the calculated &
  &forces.'),                                                                 &
  & make_keyword( 'harmonic_states_cutoff',                                   &
  &               'harmonic_states_cutoff is the number of harmonic &
  &eigenstates in the direction of each normal mode.'),                       &
  & make_keyword( 'potential_basis_cutoff',                                &
  &               'potential_basis_cutoff is the order up to which the &
  &potential is expanded. e.g. a cubic expansion would be order 3.'),         &
  & make_keyword( 'scf_convergence_threshold',                                &
  &               'scf_convergence_threshold is the energy to within which &
  &the VSCF calculation will be converged.'),                                 &
  & make_keyword( 'max_scf_cycles',                                           &
  &               'max_scf_cycles is the maximum number of SCF cycles which &
  &will be carried out as part of the VSCF calculation.'),                    &
  & make_keyword( 'perturbation_order',                                       &
  &               'perturbation_order is the order up to which perturbation &
  &theory will be run',                                                       &
  &               is_optional=.true.),                                        &
  & make_keyword( 'perturb_states_to_same_order',                             &
  &               'perturb_states_to_same_order specifies whether or not to &
  &calculate state correction at the same order as energy corrections.',      &
  &               default_value='false') ]
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_anharmonic(arguments)
  use ifile_module
  use setup_harmonic_module
  use setup_anharmonic_module
  use dictionary_module
  use structure_module
  use qpoints_module
  use normal_mode_module
  use coupling_module
  use sampling_points_module
  use linear_algebra_module
  use dft_output_file_module
  use single_mode_states_module
  use harmonic_states_module
  use vscf_states_module
  use product_states_module
  use potential_module
  use scf_module
  use perturbation_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Inputs.
  type(String) :: wd
  real(dp)     :: energy_error
  real(dp)     :: force_error
  integer      :: harmonic_states_cutoff
  integer      :: potential_basis_cutoff
  real(dp)     :: scf_convergence_threshold
  integer      :: max_scf_cycles
  integer      :: energy_perturbation_order
  integer      :: state_perturbation_order
  
  ! Previous inputs.
  type(Dictionary) :: setup_harmonic_arguments
  type(Dictionary) :: setup_anharmonic_arguments
  type(String)     :: harmonic_path
  type(String)     :: dft_code
  type(String)     :: seedname
  type(String)     :: grid_type
  
  ! Previously calculated data.
  type(StructureData)                 :: structure
  type(IFile)                         :: no_supercells_file
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
  
  ! The Born-Oppenheimer potential.
  type(PolynomialPotential) :: potential
  
  ! The number of non-translational modes.
  integer :: no_modes
  
  ! Eigenstates in various representations.
  type(HarmonicStates), allocatable :: harmonic_states(:)
  type(VscfStates),     allocatable :: vscf_states(:)
  type(ProductStates)               :: vscf_product_states
  
  ! Moller-Plesset variables.
  integer                :: no_product_states
  real(dp), allocatable  :: vscf_product_energies(:)
  real(dp), allocatable  :: perturbative_potential(:,:)
  type(RealPerturbation) :: perturbation
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Read in data.
  ! --------------------------------------------------
  
  ! Read in inputs.
  wd = arguments%value('working_directory')
  energy_error = dble(arguments%value('energy_error'))
  force_error = dble(arguments%value('force_error'))
  harmonic_states_cutoff = int(arguments%value('harmonic_states_cutoff'))
  potential_basis_cutoff = int(arguments%value('potential_basis_cutoff'))
  scf_convergence_threshold = &
     & dble(arguments%value('scf_convergence_threshold'))
  max_scf_cycles = int(arguments%value('max_scf_cycles'))
  if (arguments%is_set('perturbation_order')) then
    energy_perturbation_order = int(arguments%value('perturbation_order'))
    if (lgcl(arguments%value('perturb_states_to_same_order'))) then
      state_perturbation_order = energy_perturbation_order
    else
      state_perturbation_order = energy_perturbation_order - 1
    endif
  else
    energy_perturbation_order = 0
    state_perturbation_order = 0
  endif
  
  ! Read in setup_anharmonic settings.
  setup_anharmonic_arguments = setup_anharmonic_keywords()
  call setup_anharmonic_arguments%read_file( &
     & wd//'/setup_anharmonic.used_settings')
  harmonic_path = setup_anharmonic_arguments%value('harmonic_path')
  grid_type = setup_anharmonic_arguments%value('grid_type')
  
  ! Read in setup_harmonic settings.
  setup_harmonic_arguments = setup_harmonic_keywords()
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  dft_code = setup_harmonic_arguments%value('dft_code')
  seedname = setup_harmonic_arguments%value('seedname')
  
  ! Read in structure and supercells.
  structure = read_structure_file(harmonic_path//'/structure.dat')
  
  no_supercells_file = harmonic_path//'/no_supercells.dat'
  no_supercells = int(no_supercells_file%line(1))
  
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
    
    ! --------------------------------------------------
    ! Run calculations.
    ! --------------------------------------------------
    
    ! Calculate the Born-Oppenheimer potential.
    potential = calculate_potential( potential_basis_cutoff, &
                                   & sampling,               &
                                   & energy_error,           &
                                   & force_error,            &
                                   & modes,                  &
                                   & qpoints(i),             &
                                   & supercell)
    
    call print_line('')
    call print_line('Born-Oppenheimer Potential:')
    call print_line(potential)
    
    ! Calculate harmonic eigenstates, {|a>}, along each normal mode.
    no_modes = 0
    do j=1,size(modes)
      if (.not. modes(j)%translational_mode) then
        no_modes = no_modes+1
      endif
    enddo
    
    if (no_modes==0) then
      cycle
    endif
    
    allocate(harmonic_states(no_modes), stat=ialloc); call err(ialloc)
    k = 0
    do j=1,structure%no_modes
      if (.not. modes(j)%translational_mode) then
        k = k+1
        harmonic_states(k) = calculate_harmonic_states( &
                              & k,                      &
                              & modes(j)%frequency,     &
                              & harmonic_states_cutoff, &
                              & potential_basis_cutoff)
      endif
    enddo
    
    ! Run VSCF calculation to find VSCF eigenstates.
    vscf_states = vscf( potential,       &
                      & harmonic_states, &
                      & max_scf_cycles,  &
                      & scf_convergence_threshold)
    
    do j=1,size(vscf_states)
      call print_line('')
      call print_line('Mode '//j//' VSCF states (frequency='// &
         & modes(j)%frequency**2//'):')
      do k=0,vscf_states(j)%cutoff()
        call print_line('|'//k//'>, energy= '// vscf_states(j)%vscf_energy(k))
      enddo
    enddo
    
    ! Construct product states from VSCF basis.
    vscf_product_states = construct_product_states(vscf_states, coupling)
    
    allocate( vscf_product_energies(size(vscf_product_states)), &
            & stat=ialloc); call err(ialloc)
    do j=1,size(vscf_product_states)
      vscf_product_energies(j) = vscf_product_states%vscf_energy(j)
    enddo
    
    call print_line('')
    call print_line('VSCF product states:')
    do j=1,size(vscf_product_states)
      call print_line(vscf_product_states%state_as_ket_string(j)// &
         & ', energy= '//vscf_product_energies(j)-vscf_product_energies(1))
    enddo
    
    ! Calculate (potential - VSCF potential) in VSCF product state basis.
    no_product_states = size(vscf_product_states)
    allocate( perturbative_potential(no_product_states,no_product_states), &
            & stat=ialloc); call err(ialloc)
    
    do j=1,no_product_states
      do k=1,j
        ! Construct <i|V|j> = <i|H-T|j>
        perturbative_potential(k,j) = potential%integrate_to_constant( &
                                  &               vscf_product_states, &
                                  &               j,                   &
                                  &               k)                   &
                                  & - vscf_product_states%kinetic_energy(j,k)
        
        if (k==j) then
          ! Subtract VSCF energies from the diagonal.
          perturbative_potential(k,j) = perturbative_potential(k,j) &
                                    & - vscf_product_energies(j)
        else
          perturbative_potential(j,k) = perturbative_potential(k,j)
        endif
      enddo
    enddo
    
    call print_line('Perturbation:')
    do j=1,no_product_states
      do k=1,no_product_states
        call print_line(vscf_product_states%state_as_bra_string(j)//'V'// &
           & vscf_product_states%state_as_ket_string(k)//' = '//          &
           & perturbative_potential(j,k))
      enddo
    enddo
    
    ! Run Moller-Plesset perturbation theory to construct new basis.
    perturbation = calculate_perturbation( vscf_product_energies,     &
                                         & perturbative_potential,    &
                                         & energy_perturbation_order, &
                                         & state_perturbation_order)
    
    call print_line('')
    call print_line('Moller-Plesset corrections:')
    do j=1,size(perturbation%energy)
      call print_line(perturbation%energy(j))
    enddo
    
    ! Deallocate variables.
    deallocate( modes,                  &
              & sampling,               &
              & harmonic_states,        &
              & vscf_product_energies,  &
              & perturbative_potential, &
              & stat=ialloc); call err(ialloc)
  enddo
end subroutine
end module
