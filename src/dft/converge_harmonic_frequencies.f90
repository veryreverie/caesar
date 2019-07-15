! ======================================================================
! Converges harmonic frequencies and free energies w/r/t
!    DFT cutoff energy, k-point spacing and electronic smearing.
! Runs multiple harmonic caesar calculations with different
!    cutoff energies and k-point spacings.
! ======================================================================
module converge_harmonic_frequencies_module
  use common_module

  implicit none
  
  private
  
  public :: startup_converge_harmonic_frequencies
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext
! ----------------------------------------------------------------------
subroutine startup_converge_harmonic_frequencies()
  implicit none
  
  type(CaesarMode) :: mode
  
  type(CaesarMode) :: setup_harmonic_mode
  type(CaesarMode) :: run_harmonic_mode
  type(CaesarMode) :: calculate_normal_modes_mode
  type(CaesarMode) :: calculate_harmonic_observables_mode
  
  setup_harmonic_mode = CaesarMode('setup_harmonic')
  run_harmonic_mode = CaesarMode('run_harmonic')
  calculate_normal_modes_mode = CaesarMode('calculate_normal_modes')
  calculate_harmonic_observables_mode = &
     & CaesarMode('calculate_harmonic_observables')
  
  mode%mode_name = 'converge_harmonic_frequencies'
  mode%description = 'Converges harmonic frequencies and free energies &
     &w/r/t cutoff energy, k-point spacing and electronic smearing. N.B. only &
     &the value being converged will be changed in each input file.'
  mode%keywords = [                                                           &
     & KeywordData( 'converge_cutoff',                                      &
     &              'converge_cutoff specifies whether or not electronic &
     &cut-off energy will be converged.',                                     &
     &              default_value='true'),                                    &
     & KeywordData( 'minimum_cutoff',                                         &
     &              'minimum_cutoff is the smallest cutoff energy which will &
     &be tested. minimum_cutoff should be given in Hartree.'),                &
     & KeywordData( 'no_cutoffs',                                             &
     &              'no_cutoffs is the number of cutoff energies which will &
     &be sampled.'),                                                          &
     & KeywordData( 'maximum_cutoff',                                         &
     &              'maximum_cutoff is the cutoff energy at which calculation &
     &will be terminated if convergence has not been reached. maximum_cutoff &
     &should be given in Hartree.'),                                          &
     & KeywordData( 'converge_k-point_spacing',                               &
     &              'converge_k-point_spacing specifies whether or not &
     &k-point spacing will be converged.',                                    &
     &              default_value='true'),                                    &
     & KeywordData( 'maximum_k-point_spacing',                                &
     &              'maximum_k-point_spacing is the largest k-point spacing &
     &which will be sampled. maximum_k-point_spacing should be given in &
     &inverse Bohr.'),                                                        &
     & KeywordData( 'no_k-point_spacings',                                    &
     &              'no_k-point_spacings is the approximate number of k-point &
     &spacings which will be sampled. N.B. only k-point spacings &
     &corresponding to integer k-point grids will be sampled.'),              &
     & KeywordData( 'minimum_k-point_spacing',                                &
     &              'minimum_k-point_spacing is the k-point spacing at which &
     &calculation will be terminated if convergence has not been reached. &
     &minimum_k-point_spacing should be given in inverse Bohr.'),             &
     & KeywordData( 'k-point_parity',                                         &
     &              'k-point_parity specifies which k-point grids are &
     &sampled. If k-point_parity is "odd" or "even" then only k-point grids &
     &with an odd or even number of k-points in all directions are sampled. &
     &If k-point_parity is "both" then all k-point grids are samples.',       &
     &              default_value='both'),                                    &
     & KeywordData( 'converge_smearing',                                      &
     &              'converge_smearing specifies whether or not electron &
     &smearing will be converged.',                                           &
     &              default_value='true'),                                    &
     & KeywordData( 'maximum_smearing',                                       &
     &              'maximum_smearing is the largest electronic smearing &
     &which will be sampled. maximum_smearing should be given in Hartree.'),  &
     & KeywordData( 'no_smearings',                                           &
     &              'no_smearings is the number of electronic smearings which &
     &which will be sampled.'),                                               &
     & KeywordData( 'minimum_smearing',                                       &
     &              'minimum_smearing is the electronic smearing at which &
     &calculation will be terminated if convergence has not been reached. &
     &minimum_smearing should be given in Hartree.'),                         &
     & KeywordData( 'converge_energies',                                      &
     &              'converge_energies determines whether convergence of the &
     &free energies is monitored as well as convergence of harmonic &
     &frequencies.',                                                          &
     &              default_value='false'),                                   &
     & KeywordData( 'freq_tolerance',                                         &
     &              'freq_tolerance is the tolerance used in the convergence &
     &of the harmonic frequencies.',                                          &
     &              default_value='1e-8'),                                    &
     & KeywordData( 'energy_tolerance',                                       &
     &              'energy_tolerance is the tolerance used in the &
     &convergence of the vibrational free energies. It is given in Hartree &
     &per primitive unit cell.',                                              &
     &              default_value='1e-6'),                                    &
     & KeywordData( 'convergence_count',                                      &
     &              'convergence_count is the number of consecutive &
     &frequencies and/or energies that must be within the defined tolerance &
     &of each other for convergence to be considered reached',                &
     &              default_value='3')                                        ]
  mode%keywords = [ mode%keywords,                               &
                  & setup_harmonic_mode%keywords,                &
                  & run_harmonic_mode%keywords,                  &
                  & calculate_normal_modes_mode%keywords,        &
                  & calculate_harmonic_observables_mode%keywords ]
  mode%main_subroutine => converge_harmonic_frequencies_subroutine
  
  call mode%remove_keyword('supercells_to_run')
  call mode%remove_keyword('calculations_to_run')
  call mode%remove_keyword('exit_on_error')
  
  call add_mode(mode)
end subroutine

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine converge_harmonic_frequencies_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! User inputs.
  type(String) :: seedname
  type(String) :: file_type
  
  logical      :: converge_cutoff
  real(dp)     :: minimum_cutoff
  integer      :: no_cutoffs
  real(dp)     :: maximum_cutoff
  
  logical      :: converge_spacing
  real(dp)     :: maximum_spacing
  integer      :: no_spacings
  real(dp)     :: minimum_spacing
  type(String) :: kpoint_parity
  
  logical      :: converge_smearing
  real(dp)     :: minimum_smearing
  integer      :: no_smearings
  real(dp)     :: maximum_smearing
  
  integer      :: grid(3)
  logical      :: converge_energies
  logical      :: repeat_calculations
  real(dp)     :: freq_tolerance
  real(dp)     :: energy_tolerance
  integer      :: convergence_count
  integer      :: random_seed
  
  ! Structure.
  type(String)        :: input_filename
  type(StructureData) :: structure

  ! Energy cut-offs and kpoint spacings
  real(dp),         allocatable :: cutoffs(:)
  real(dp),         allocatable :: kpoint_spacings(:)
  type(KpointGrid), allocatable :: kpoint_grids(:)
  logical,          allocatable :: unique_grids(:)
  real(dp),         allocatable :: smearings(:)
  
  ! Variables to converge.
  type(RealVector), allocatable :: cutoff_frequencies(:)
  type(RealVector), allocatable :: cutoff_free_energies(:)
  type(RealVector), allocatable :: kpoint_frequencies(:)
  type(RealVector), allocatable :: kpoint_free_energies(:)
  type(RealVector), allocatable :: smearing_frequencies(:)
  type(RealVector), allocatable :: smearing_free_energies(:)
  
  ! Working variables
  integer :: no_qpoints
  logical :: frequencies_converged
  logical :: energies_converged
  
  ! Random generator for generating seed if not set.
  type(RandomReal) :: random_generator
  
  ! Files and directories.
  type(String) :: directory
  
  ! Temporary variables
  integer :: i,j,ialloc
  
  ! --------------------------------------------------
  ! Get settings from user, and check them.
  ! --------------------------------------------------
  file_type = arguments%value('file_type')
  seedname = arguments%value('seedname')
  
  converge_cutoff = lgcl(arguments%value('converge_cutoff'))
  minimum_cutoff = dble(arguments%value('minimum_cutoff'))
  no_cutoffs = int(arguments%value('no_cutoffs'))
  maximum_cutoff = dble(arguments%value('maximum_cutoff'))
  
  converge_spacing = lgcl(arguments%value('converge_k-point_spacing'))
  maximum_spacing = dble(arguments%value('maximum_k-point_spacing'))
  no_spacings = int(arguments%value('no_k-point_spacings'))
  minimum_spacing = dble(arguments%value('minimum_k-point_spacing'))
  kpoint_parity = arguments%value('k-point_parity')
  
  converge_smearing = lgcl(arguments%value('converge_smearing'))
  maximum_smearing = dble(arguments%value('maximum_smearing'))
  no_smearings = int(arguments%value('no_smearings'))
  minimum_smearing = dble(arguments%value('minimum_smearing'))
  
  grid = int(split_line(arguments%value('q-point_grid')))
  no_qpoints = grid(1)*grid(2)*grid(3)
  converge_energies = lgcl(arguments%value('converge_energies'))
  repeat_calculations = lgcl(arguments%value('repeat_calculations'))
  freq_tolerance = dble(arguments%value('freq_tolerance'))
  energy_tolerance = dble(arguments%value('energy_tolerance'))
  convergence_count = int(arguments%value('convergence_count'))
  
  if (arguments%is_set('random_seed')) then
    random_generator = RandomReal(int(arguments%value('random_seed')))
  else
    random_generator = RandomReal()
  endif
  random_seed = random_generator%get_seed()
  
  if (maximum_cutoff<=minimum_cutoff) then
    call print_line(ERROR//': maximum_cutoff is smaller than minimum_cutoff.')
    call quit()
  elseif (maximum_spacing<=minimum_spacing) then
    call print_line(ERROR//': maximum_k-point_spacing is smaller than &
       &minimum_k-point_spacing.')
    call quit()
  elseif (maximum_smearing<=minimum_smearing) then
    call print_line(ERROR//': maximum_smearing is smaller than &
       &minimum_smearing.')
    call quit()
  elseif (file_type/='castep' .and. file_type/='quantum_espresso') then
    call print_line(ERROR//': castep and quantum_espresso are the only &
       &accepted file types for this mode.')
    call quit()
  elseif (convergence_count<1) then
    call print_line(ERROR//': convergence_count must be at least 1.')
    call quit()
  endif
  
  ! Read in structure.
  input_filename = make_input_filename(file_type, seedname)
  structure = input_file_to_StructureData(file_type, input_filename)
  
  ! Initialise arrays.
  allocate( cutoff_frequencies(0),     &
          & cutoff_free_energies(0),   &
          & kpoint_frequencies(0),     &
          & kpoint_free_energies(0),   &
          & smearing_frequencies(0),   &
          & smearing_free_energies(0), &
          & stat=ialloc); call err(ialloc)
  
  ! Calculate spacings, smearings and cutoffs.
  allocate( cutoffs(no_cutoffs),          &
          & kpoint_spacings(no_spacings), &
          & kpoint_grids(no_spacings),    &
          & smearings(no_smearings),      &
          & stat=ialloc); call err(ialloc)
  do i=1,no_cutoffs
    cutoffs(i) = ( (no_cutoffs-i)*minimum_cutoff   &
             &   + (i-1)*maximum_cutoff          ) &
             & / (no_cutoffs-1)
  enddo
  
  do i=1,no_smearings
    smearings(i) = ( (no_smearings-i)*maximum_smearing &
               &   + (i-1)*minimum_smearing          ) &
               & / (no_smearings-1)
  enddo
  
  do i=1,no_spacings
    ! Evenly space k-point spacings in 1/spacing co-ordinates.
    kpoint_spacings(i) = (no_spacings-1)                   &
                     & / ( (no_spacings-i)/maximum_spacing &
                     &   + (i-1)/minimum_spacing           )
    
    ! Calculate the k-point grid corresponding to each spacing.
    kpoint_grids(i) = calculate_kpoint_grid( kpoint_spacings(i),           &
                                           & dble(structure%recip_lattice) )
    
    ! Increase any grid elements which have the wrong parity.
    if (kpoint_parity=='both') then
      continue
    elseif (kpoint_parity=='odd') then
      do j=1,3
        if (modulo(kpoint_grids(i)%grid(j),2)==0) then
          kpoint_grids(i)%grid(j) = kpoint_grids(i)%grid(j)+1
        endif
      enddo
    elseif (kpoint_parity=='even') then
      do j=1,3
        if (modulo(kpoint_grids(i)%grid(j),2)==1) then
          kpoint_grids(i)%grid(j) = kpoint_grids(i)%grid(j)+1
        endif
      enddo
    else
      call print_line(ERROR//': Unexpected value for k-point_parity.')
      call quit()
    endif
    
    ! Re-calculate spacing to reflect the actual k-point grid.
    kpoint_spacings(i) = calculate_kpoint_spacing( &
                   & kpoint_grids(i),              &
                   & dble(structure%recip_lattice) )
  enddo
  
  ! Remove duplicate k-point grids.
  unique_grids = [ .true.,                                 &
                 & ( kpoint_grids(i)/=kpoint_grids(i-1),   &
                 &   i=2,                                  &
                 &   no_spacings                         ) ]
  kpoint_spacings = kpoint_spacings(filter(unique_grids))
  kpoint_grids = kpoint_grids(filter(unique_grids))
  no_spacings = size(kpoint_spacings)
  
  ! --------------------------------------------------
  ! Run cut-off energy convergence.
  ! --------------------------------------------------
  if (converge_cutoff) then
    ! Initialise convergence checks and observable arrays.
    frequencies_converged = .false.
    energies_converged = .false.
    if (.not. converge_energies) then
     energies_converged = .true.
    endif
    
    ! Loop over cutoffs, running calculations at each.
    do i=1,no_cutoffs
      directory = 'cutoff_'//trim(str(                           &
         & cutoffs(i),                                           &
         & settings=PrintSettings( decimal_places        = 2,    &
         &                         floating_point_format = 'f' ) ))
      call mkdir(directory)
      
      cutoff_frequencies = [                               &
         & cutoff_frequencies,                             &
         & calculate_frequencies( directory  = directory,  &
         &                        cutoff     = cutoffs(i), &
         &                        seedname   = seedname,   &
         &                        no_qpoints = no_qpoints, &
         &                        file_type  = file_type,  &
         &                        arguments  = arguments   ) ]
      
      if (i>convergence_count) then
        frequencies_converged =                                    &
           & all( maximum_difference(                              &
           &         cutoff_frequencies(i),                        &
           &         cutoff_frequencies(i-convergence_count:i-1) ) &
           &    < freq_tolerance                                   )
      endif

      if (converge_energies .and. .not. energies_converged) then
        cutoff_free_energies = [                            &
           & cutoff_free_energies,                          &
           & calculate_free_energies( directory,            &
           &                          repeat_calculations,  &
           &                          random_seed,          &
           &                          arguments             ) ]

        if (i>convergence_count) then
          energies_converged =                                         &
             & all( maximum_difference(                                &
             &         cutoff_free_energies(i),                        &
             &         cutoff_free_energies(i-convergence_count:i-1) ) &
             &    < energy_tolerance                                   )
        endif
      endif
      
      ! Write output file.
      call write_output_file( structure,             &
                            & cutoffs,               &
                            & cutoff_frequencies,    &
                            & cutoff_free_energies,  &
                            & kpoint_spacings,       &
                            & kpoint_frequencies,    &
                            & kpoint_free_energies,  &
                            & smearings,             &
                            & smearing_frequencies,  &
                            & smearing_free_energies )
      
      ! Check for convergence.
      call print_line('Cut-off energy: '//cutoffs(i)//' (Ha)')
      if (frequencies_converged .and. energies_converged) then
        call print_line('Convergence reached.')
        exit
      endif
    enddo
    
    if (.not. (frequencies_converged .or. energies_converged)) then
      call print_line('Convergence not reached.')
    endif
  endif
  
  ! --------------------------------------------------
  ! Run k-points spacing convergence.
  ! --------------------------------------------------
  if (converge_spacing) then
    ! Initialise convergence checks and observable arrays.
    frequencies_converged = .false.
    energies_converged = .false.
    if (.not. converge_energies) then
     energies_converged = .true.
    endif
    
    ! Loop over spacings, running calculations at each.
    do i=1,no_spacings
      directory = 'kpoints_'//join(str(kpoint_grids(i)%grid), delimiter='_')
      call mkdir(directory)
      
      kpoint_frequencies = [                                     &
         & kpoint_frequencies,                                   &
         & calculate_frequencies( directory   = directory,       &
         &                        kpoint_grid = kpoint_grids(i), &
         &                        seedname    = seedname,        &
         &                        no_qpoints  = no_qpoints,      &
         &                        file_type   = file_type,       &
         &                        arguments   = arguments        ) ]
      
      if (i>convergence_count) then
        frequencies_converged =                                    &
           & all( maximum_difference(                              &
           &         kpoint_frequencies(i),                        &
           &         kpoint_frequencies(i-convergence_count:i-1) ) &
           &    < freq_tolerance                                   )
      endif

      if (converge_energies .and. .not. energies_converged) then
        kpoint_free_energies = [                            &
           & kpoint_free_energies,                          &
           & calculate_free_energies( directory,            &
           &                          repeat_calculations,  &
           &                          random_seed,          &
           &                          arguments             ) ]

        if (i>convergence_count) then
          energies_converged =                                         &
             & all( maximum_difference(                                &
             &         kpoint_free_energies(i),                        &
             &         kpoint_free_energies(i-convergence_count:i-1) ) &
             &    < energy_tolerance                                   )
        endif
      endif
      
      ! Write output file.
      call write_output_file( structure,             &
                            & cutoffs,               &
                            & cutoff_frequencies,    &
                            & cutoff_free_energies,  &
                            & kpoint_spacings,       &
                            & kpoint_frequencies,    &
                            & kpoint_free_energies,  &
                            & smearings,             &
                            & smearing_frequencies,  &
                            & smearing_free_energies )
      
      call print_line('K-point spacing: '//kpoint_spacings(i)//' (Bohr^-1)')
      if (frequencies_converged .and. energies_converged) then
        call print_line('Convergence reached.')
        exit
      endif
    enddo
    
    if (.not. (frequencies_converged .or. energies_converged)) then
      call print_line('Convergence not reached.')
    endif
  endif
  
  ! Run electronic smearing convergence.
  if (converge_smearing) then
    ! Initialise convergence checks and observable arrays.
    frequencies_converged = .false.
    energies_converged = .false.
    if (.not. converge_energies) then
     energies_converged = .true.
    endif
    
    ! Loop over cutoffs, running calculations at each.
    do i=1,no_smearings
      directory = 'smearing_'//trim(str(                         &
         & smearings(i),                                         &
         & settings=PrintSettings( decimal_places        = 6,    &
         &                         floating_point_format = 'f' ) ))
      call mkdir(directory)
      
      smearing_frequencies = [                               &
         & smearing_frequencies,                             &
         & calculate_frequencies( directory  = directory,    &
         &                        smearing   = smearings(i), &
         &                        seedname   = seedname,     &
         &                        no_qpoints = no_qpoints,   &
         &                        file_type  = file_type,    &
         &                        arguments  = arguments     ) ]
      
      if (i>convergence_count) then
        frequencies_converged =                                    &
           & all( maximum_difference(                              &
           &         smearing_frequencies(i),                        &
           &         smearing_frequencies(i-convergence_count:i-1) ) &
           &    < freq_tolerance                                   )
      endif

      if (converge_energies .and. .not. energies_converged) then
        smearing_free_energies = [                          &
           & smearing_free_energies,                        &
           & calculate_free_energies( directory,            &
           &                          repeat_calculations,  &
           &                          random_seed,          &
           &                          arguments             ) ]

        if (i>convergence_count) then
          energies_converged =                                         &
             & all( maximum_difference(                                &
             &         smearing_free_energies(i),                        &
             &         smearing_free_energies(i-convergence_count:i-1) ) &
             &    < energy_tolerance                                   )
        endif
      endif
      
      ! Write output file.
      call write_output_file( structure,             &
                            & cutoffs,               &
                            & cutoff_frequencies,    &
                            & cutoff_free_energies,  &
                            & kpoint_spacings,       &
                            & kpoint_frequencies,    &
                            & kpoint_free_energies,  &
                            & smearings,             &
                            & smearing_frequencies,  &
                            & smearing_free_energies )
      
      ! Check for convergence.
      call print_line('Electronic smearing: '//smearings(i)//' (Ha)')
      if (frequencies_converged .and. energies_converged) then
        call print_line('Convergence reached.')
        exit
      endif
    enddo
    
    if (.not. (frequencies_converged .or. energies_converged)) then
      call print_line('Convergence not reached.')
    endif
    
  endif
end subroutine

subroutine write_output_file(structure,cutoffs,cutoff_frequencies, &
   & cutoff_free_energies,kpoint_spacings,kpoint_frequencies,      &
   & kpoint_free_energies,smearings,smearing_frequencies,          &
   & smearing_free_energies)
  implicit none
  
  type(StructureData), intent(in) :: structure
  real(dp),            intent(in) :: cutoffs(:)
  type(RealVector),    intent(in) :: cutoff_frequencies(:)
  type(RealVector),    intent(in) :: cutoff_free_energies(:)
  real(dp),            intent(in) :: kpoint_spacings(:)
  type(RealVector),    intent(in) :: kpoint_frequencies(:)
  type(RealVector),    intent(in) :: kpoint_free_energies(:)
  real(dp),            intent(in) :: smearings(:)
  type(RealVector),    intent(in) :: smearing_frequencies(:)
  type(RealVector),    intent(in) :: smearing_free_energies(:)
  
  type(OFile)  :: output_file
  
  integer :: i
  
  ! Open output file.
  output_file = OFile('convergence.dat')
  
  call output_file%print_line('No. atoms   : '//structure%no_atoms)
  call output_file%print_line('Cell volume : '//structure%volume)
  call output_file%print_line('')
  
  if (size(cutoff_frequencies)>0) then
    call output_file%print_line('Cutoff energy (Ha) | Mode frequencies (Ha)')
    do i=1,size(cutoff_frequencies)
      call output_file%print_line(cutoffs(i)//' '//cutoff_frequencies(i))
    enddo
  endif
  
  if (size(cutoff_free_energies)>0) then
    call output_file%print_line('')
    call output_file%print_line('Cutoff energy (Ha) | Free energies &
       &(Ha per primitive cell)')
    do i=1,size(cutoff_free_energies)
      call output_file%print_line(cutoffs(i)//' '//cutoff_free_energies(i))
    enddo
  endif
  
  if (size(cutoff_frequencies)>0 .and. size(kpoint_frequencies)>0) then
    call output_file%print_line('')
  endif
  
  if (size(kpoint_frequencies)>0) then
    call output_file%print_line('k-point spacing (Bohr^-1) | Mode frequencies &
       &(Ha)')
    do i=1,size(kpoint_frequencies)
      call output_file%print_line( kpoint_spacings(i) //' '// &
                                 & kpoint_frequencies(i)      )
    enddo
  endif
  
  if (size(kpoint_free_energies)>0) then
    call output_file%print_line('')
    call output_file%print_line('k-point spacing (Bohr^-1) | Free energies &
       &(Ha per primitive cell)')
    do i=1,size(kpoint_free_energies)
      call output_file%print_line( kpoint_spacings(i)//' '// &
                                 & kpoint_free_energies(i)   )
    enddo
  endif
  
  if ( ( size(cutoff_frequencies)>0 .or.         &
     &   size(kpoint_frequencies)>0      ) .and. &
     & size(smearing_frequencies)>0              ) then
    call output_file%print_line('')
  endif
  
  if (size(smearing_frequencies)>0) then
    call output_file%print_line('Electronic smearing (Hartree) | Mode &
       &frequencies (Ha)')
    do i=1,size(smearing_frequencies)
      call output_file%print_line( smearings(i)    //' '// &
                                 & smearing_frequencies(i) )
    enddo
  endif
  
  if (size(smearing_free_energies)>0) then
    call output_file%print_line('')
    call output_file%print_line('Electronic smearing (Hartree) | Free &
       &energies (Ha per primitive cell)')
    do i=1,size(smearing_free_energies)
      call output_file%print_line( smearings(i)      //' '// &
                                 & smearing_free_energies(i) )
    enddo
  endif
end subroutine

function calculate_frequencies(directory,cutoff,kpoint_grid,smearing,     &
   & seedname,no_qpoints,file_type,arguments) result(output)
  implicit none
  
  type(String),     intent(in)           :: directory
  real(dp),         intent(in), optional :: cutoff
  type(KpointGrid), intent(in), optional :: kpoint_grid
  real(dp),         intent(in), optional :: smearing
  type(String),     intent(in)           :: seedname
  integer,          intent(in)           :: no_qpoints
  type(String),     intent(in)           :: file_type
  type(Dictionary), intent(in)           :: arguments
  type(RealVector)                       :: output
  
  type(Dictionary) :: setup_harmonic_arguments
  type(Dictionary) :: run_harmonic_arguments
  type(Dictionary) :: calculate_normal_modes_arguments
  
  type(String) :: qpoint_dir
  type(IFile)  :: complex_modes_file
  
  type(ComplexMode), allocatable :: modes(:)
  
  integer :: i,ialloc
  
  ! Check input.
  if (count([ present(cutoff),      &
            & present(kpoint_grid), &
            & present(smearing)     ])/=1) then
    call print_line(CODE_ERROR//': calculate_frequencies must be passed &
       &exactly one of cutoff, kpoint_grid and smearing.')
  endif
  
  ! Write DFT files
  if (file_type=='castep') then
    call write_castep_files(directory, seedname, cutoff, kpoint_grid, smearing)
  elseif (file_type=='quantum_espresso') then
    call write_qe_file(directory, seedname, cutoff, kpoint_grid, smearing)
  endif
  
  ! Call setup_harmonic.
  setup_harmonic_arguments = Dictionary(CaesarMode('setup_harmonic'))
  call setup_harmonic_arguments%set(arguments)
  call setup_harmonic_arguments%set('working_directory', directory)
  call setup_harmonic_arguments%set( 'output_file',                      &
                                   & directory//'/setup_harmonic.output' )
  call call_caesar('setup_harmonic', setup_harmonic_arguments)
  
  ! Call run_harmonic.
  run_harmonic_arguments = Dictionary(CaesarMode('run_harmonic'))
  call run_harmonic_arguments%set(arguments)
  call run_harmonic_arguments%set('working_directory', directory)
  call run_harmonic_arguments%set( 'output_file',                    &
                                 & directory//'/run_harmonic.output' )
  call call_caesar('run_harmonic', run_harmonic_arguments)
  
  ! Call calculate_normal_modes.
  calculate_normal_modes_arguments = &
     & Dictionary(CaesarMode('calculate_normal_modes'))
  call calculate_normal_modes_arguments%set(arguments)
  call calculate_normal_modes_arguments%set('working_directory', directory)
  call calculate_normal_modes_arguments%set(       &
     & 'output_file',                              &
     & directory//'/calculate_normal_modes.output' )
  call call_caesar('calculate_normal_modes', calculate_normal_modes_arguments)
  
  ! Read in normal modes.
  allocate(modes(0), stat=ialloc); call err(ialloc)
  do i=1,no_qpoints
    qpoint_dir = directory//'/qpoint_'//left_pad(i, str(no_qpoints))
    complex_modes_file = IFile(qpoint_dir//'/complex_modes.dat')
    modes = [modes, ComplexMode(complex_modes_file%sections())]
  enddo
  
  output = vec(modes%frequency)
end function

function calculate_free_energies(directory,repeat_calculations,random_seed, &
   & arguments) result(output)
  implicit none
  
  type(String),     intent(in) :: directory
  logical,          intent(in) :: repeat_calculations
  integer,          intent(in) :: random_seed
  type(Dictionary), intent(in) :: arguments
  type(RealVector)             :: output
  
  type(Dictionary) :: calculate_harmonic_observables_arguments
  
  type(String)                         :: file_name
  type(IFile)                          :: thermodynamics_file
  type(String),            allocatable :: lines(:)
  type(ThermodynamicData), allocatable :: thermodynamics(:)

  ! Call calculate_harmonic_observables.
  file_name = directory//'/harmonic_observables/thermodynamic_variables.dat'
  if (repeat_calculations .or. .not. file_exists(file_name)) then
    calculate_harmonic_observables_arguments = &
       & Dictionary(CaesarMode('calculate_harmonic_observables'))
    call calculate_harmonic_observables_arguments%set(arguments)
    call calculate_harmonic_observables_arguments%set( 'working_directory', &
                                                     & directory            )
    call calculate_harmonic_observables_arguments%set(       &
       & 'output_file',                                      &
       & directory//'/calculate_harmonic_observables.output' )
    call calculate_harmonic_observables_arguments%set( 'random_seed',   &
                                                     & str(random_seed) )
    call call_caesar( 'calculate_harmonic_observables',        &
                    & calculate_harmonic_observables_arguments )
  endif
  
  ! Read in thermodynamic data.
  thermodynamics_file = IFile(file_name)
  lines = thermodynamics_file%lines()
  thermodynamics = ThermodynamicData(lines(2:))
  
  output = vec(thermodynamics%free_energy)
end function

subroutine write_castep_files(directory,seedname,cutoff,kpoint_grid,smearing)
  implicit none
  
  type(String),     intent(in)           :: directory
  type(String),     intent(in)           :: seedname
  real(dp),         intent(in), optional :: cutoff
  type(KpointGrid), intent(in), optional :: kpoint_grid
  real(dp),         intent(in), optional :: smearing
  
  type(IFile) :: input_cell_file
  type(OFile) :: output_cell_file
  type(IFile) :: input_param_file
  type(OFile) :: output_param_file
  
  logical :: line_replaced
  
  type(String), allocatable :: lines(:)
  type(String), allocatable :: line(:)
  
  integer :: i
  
  ! --------------------------------------------------
  ! Update cell file.
  ! --------------------------------------------------
  
  ! Read input cell file.
  input_cell_file = IFile(seedname//'.cell')
  lines = input_cell_file%lines()
  
  ! Update k-point spacing.
  if (present(kpoint_grid)) then
    line_replaced = .false.
    do i=1,size(lines)
      line = split_line(lower_case(lines(i)))
      if (size(line)>=1) then
        if ( line(1)=='kpoint_mp_spacing'  .or. &
           & line(1)=='kpoint_mp_grid'     .or. &
           & line(1)=='kpoints_mp_spacing' .or. &
           & line(1)=='kpoints_mp_grid'         ) then
          lines(i) = 'kpoint_mp_grid '//kpoint_grid
          line_replaced = .true.
          exit
        endif
      endif
    enddo
    
    if (.not. line_replaced) then
      lines = [lines, 'kpoint_mp_grid '//kpoint_grid]
    endif
  endif
  
  ! Write output cell file.
  output_cell_file = OFile(directory//'/'//seedname//'.cell')
  call output_cell_file%print_lines(lines)
  
  ! --------------------------------------------------
  ! Update param file.
  ! --------------------------------------------------
  
  ! Read input param file.
  input_param_file = IFile(seedname//'.param')
  lines = input_param_file%lines()
  
  ! Update electronic cut-off energy.
  if (present(cutoff)) then
    line_replaced = .false.
    do i=1,size(lines)
      line = split_line(lower_case(lines(i)))
      if (line(1)=='cut_off_energy') then
        lines(i) = 'cut_off_energy '//cutoff*EV_PER_HARTREE
        line_replaced = .true.
        exit
      endif
    enddo
    
    if (.not. line_replaced) then
      lines = [lines, 'cut_off_energy '//cutoff*EV_PER_HARTREE]
    endif
  endif
  
  ! Update electronic smearing.
  if (present(smearing)) then
    line_replaced = .false.
    do i=1,size(lines)
      line = split_line(lower_case(lines(i)))
      if (line(1)=='smearing_width') then
        lines(i) = 'smearing_width '//smearing/KB_IN_AU//' K'
        line_replaced = .true.
        exit
      endif
    enddo
    
    if (.not. line_replaced) then
      lines = [lines, 'smearing_width '//smearing/KB_IN_AU//' K']
    endif
  endif
  
  ! Write output param file.
  output_param_file = OFile(directory//'/'//seedname//'.param')
  call output_param_file%print_lines(lines)
end subroutine

subroutine write_qe_file(directory,seedname,cutoff,kpoint_grid,smearing)
  implicit none
  
  type(String),     intent(in)           :: directory
  type(String),     intent(in)           :: seedname
  real(dp),         intent(in), optional :: cutoff
  type(KpointGrid), intent(in), optional :: kpoint_grid
  real(dp),         intent(in), optional :: smearing
  
  type(IFile)       :: input_file
  type(QeInputFile) :: qe_input_file
  type(OFile)       :: output_file
  
  integer  :: ecutwfc_line
  integer  :: ecutrho_line
  real(dp) :: ecutwfc
  real(dp) :: ecutrho
  
  integer :: smearing_line
  
  type(String), allocatable :: line(:)
  
  integer :: i
  
  ! Read input file.
  input_file = IFile(seedname//'.in')
  qe_input_file = QeInputFile(input_file%lines())
  
  ! Update k-point grid.
  if (present(kpoint_grid)) then
    line = split_line(qe_input_file%k_points(2))
    if (size(line)==3) then
      qe_input_file%k_points(2) = str(kpoint_grid)
    elseif (size(line)==6) then
      qe_input_file%k_points(2) = join([str(kpoint_grid), line(4:6)])
    else
      call print_line(ERROR//': k_points card has an unexpected number of &
         &entries.')
      call quit()
    endif
  endif
  
  ! Update electronic cutoff.
  if (present(cutoff)) then
    ecutwfc_line = 0
    ecutrho_line = 0
    do i=1,size(qe_input_file%namelists)
      line = split_line(lower_case(qe_input_file%namelists(i)), delimiter='=')
      if (size(line)==2) then
        if (line(1)=='ecutwfc') then
          if (ecutwfc_line/=0) then
            call print_line(ERROR//': ecutwfc appears twice in input file.')
            call quit()
          endif
          ecutwfc_line = i
        elseif (line(1)=='ecutrho') then
          if (ecutrho_line/=0) then
            call print_line(ERROR//': ecutrho appears twice in input file.')
            call quit()
          endif
          ecutrho_line = i
        endif
      endif
    enddo
    
    if (ecutwfc_line==0) then
      call print_line(ERROR//': ecutwfc not present.')
      call err()
    endif
    
    line = split_line(qe_input_file%namelists(ecutwfc_line), delimiter='=')
    ecutwfc = dble(line(2))
    qe_input_file%namelists(ecutwfc_line) = &
       & 'ecutwfc='//cutoff*RYDBERG_PER_HARTREE//','
    
    if (ecutrho_line/=0) then
      line = split_line(qe_input_file%namelists(ecutrho_line), delimiter='=')
      ecutrho = dble(line(2))
      qe_input_file%namelists(ecutrho_line) = 'ecutrho='// &
         & ecutrho*cutoff/ecutwfc*RYDBERG_PER_HARTREE//','
    endif
  endif
  
  ! Update electronic smearing.
  if (present(smearing)) then
    smearing_line = 0
    do i=1,size(qe_input_file%namelists)
      line = split_line(qe_input_file%namelists(i), delimiter='=')
      if (line(1)=='degauss') then
        smearing_line = i
      endif
    enddo
    
    if (smearing_line==0) then
      call print_line(ERROR//': "degauss" not present in .in file.')
      call err()
    else
      qe_input_file%namelists(smearing_line) = 'degauss='// &
         & smearing*RYDBERG_PER_HARTREE//','
    endif
  endif
  
  ! Write output file.
  output_file = OFile(directory//'/'//seedname//'.in')
  call output_file%print_lines(qe_input_file)
end subroutine

! Find the maximum element-wise difference between two vectors.
impure elemental function maximum_difference(this,that) result(output)
  implicit none
  
  type(RealVector), intent(in) :: this
  type(RealVector), intent(in) :: that
  real(dp)                     :: output
  
  if (size(this)/=size(that)) then
    call print_line(CODE_ERROR//': Vectors of different lengths.')
    call err()
  endif
  
  output = maxval(abs(dble(this-that)))
end function
end module
