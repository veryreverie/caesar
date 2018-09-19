! ======================================================================
! Converges the cutoff energy.
! ======================================================================
! Runs multiple single-point calculations with different cutoff energies.
module converge_cutoff_and_kpoints_module
  use common_module
  
  use castep_output_file_module
  implicit none
  
  private
  
  public :: converge_cutoff_and_kpoints
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function converge_cutoff_and_kpoints() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'converge_cutoff_and_kpoints'
  output%description = 'Converges cutoff energy and k-point spacing.'
  output%keywords = [                                                         &
  & KeywordData( 'seedname',                                                 &
  &               'seedname is the CASTEP seedname from which file names are &
  &constructed.'),                                                            &
  & KeywordData( 'run_script',                                               &
  &               'run_script is the path to the script for running CASTEP. &
  &An example run script can be found in doc/input_files.',                   &
  &               is_path=.true.),                                            &
  & KeywordData( 'no_cores',                                                 &
  &               'no_cores is the number of cores on which CASTEP will be &
  &run. This is passed to the specified run script.',                         &
  &               default_value='1'),                                         &
  & KeywordData( 'minimum_cutoff',                                           &
  &               'minimum_cutoff is the smallest cutoff energy which will be &
  &tested. minimum_cutoff must be an integer.',                               &
  &               default_value='300'),                                       &
  & KeywordData( 'cutoff_step',                                              &
  &               'cutoff_step is the step between each cutoff energy which &
  &will be tested. cutoff_step must be an integer.',                          &
  &               default_value='50'),                                        &
  & KeywordData( 'maximum_cutoff',                                           &
  &               'maximum_cutoff is the cutoff energy at which calculation &
  &will be terminated if convergence has not been reached. maximum_cutoff &
  &must be an integer.',                                                      &
  &               default_value='1500'),                                      &
  & KeywordData( 'minimum_kpoints',                                          &
  &               'minimum_kpoints is the smallest number of k-points which &
  &will be tested. This is the average number of k-points in each direction. &
  &minimum_kpoints must be an integer.',                                      &
  &               default_value='1'),                                         &
  & KeywordData( 'kpoints_step',                                             &
  &               'kpoints_step is the step between each number of k-points &
  &which will be tested. kpoints_step must be an integer, and should be even &
  &if consistent sampling of the gamma point is desirable.',                  &
  &               default_value='2'),                                         &
  & KeywordData( 'maximum_kpoints',                                          &
  &               'maximum_kpoints is the number of k-points at which &
  &calculation will be terminated if convergence has not been reached. &
  &maximum_kpoints must be an integer.',                                      &
  &               default_value='30'),                                        &
  & KeywordData( 'energy_convergence_threshold',                             &
  &               'energy_convergence_threshold is the accuracy to which the &
  &energy  must converge in order for the cutoff energy to be accepted.'),    &
  & KeywordData( 'force_convergence_threshold',                              &
  &               'force_convergence_threshold is the accuracy to which the &
  &forces must converge in order for the cutoff energy to be accepted.'),     &
  & KeywordData( 'no_converged_calculations',                                &
  &               'no_converged_calculations is the number of calculations &
  &which must be within both thresholds in order to accept convergence. The &
  &cutoff energy and k-point spacing of the first such calculation will be &
  &accepted.',                                                                &
  &               default_value='10'),                                        &
  & KeywordData( 'symmetry_precision',                                        &
  &              'In order for a symmetry to be accepted, it must transform &
  &the position of every atom to within symmetry_precision of an atom of the &
  &same element. symmetry_precision should be given in Bohr.',                &
  &              default_value='0.1') ]
  output%main_subroutine => converge_cutoff_and_kpoints_subroutine
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine converge_cutoff_and_kpoints_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! User input variables.
  type(String) :: wd
  type(String) :: seedname
  type(String) :: run_script
  integer      :: no_cores
  integer      :: minimum_cutoff
  integer      :: cutoff_step
  integer      :: maximum_cutoff
  integer      :: minimum_kpoints
  integer      :: kpoints_step
  integer      :: maximum_kpoints
  real(dp)     :: energy_convergence_threshold
  real(dp)     :: force_convergence_threshold
  integer      :: no_converged_calculations
  real(dp)     :: symmetry_precision
  
  ! Electronic structure calculation handlers.
  type(CalculationWriter) :: calculation_writer
  type(CalculationRunner) :: calculation_runner
  type(CalculationReader) :: calculation_reader
  
  ! The structure being converged.
  type(StructureData) :: structure
  real(dp)            :: recip_lattice(3,3)
  real(dp)            :: average_reciprocal_length
  
  ! Files.
  type(IFile)            :: cell_file
  type(IFile)            :: param_file
  type(CastepOutputFile) :: castep_file
  type(OFile)            :: progress_file
  type(OFile)            :: cutoff_file
  type(OFile)            :: kpoints_file
  
  ! Energies and forces.
  integer                       :: no_cutoffs
  integer                       :: no_kpoints
  integer,          allocatable :: cutoffs(:)
  real(dp),         allocatable :: kpoint_spacings(:)
  integer,          allocatable :: kpoint_numbers(:)
  integer,          allocatable :: kpoint_grids(:,:)
  real(dp),         allocatable :: cutoffs_energies(:)
  real(dp),         allocatable :: kpoints_energies(:)
  type(RealVector), allocatable :: cutoffs_forces(:,:)
  type(RealVector), allocatable :: kpoints_forces(:,:)
  
  ! Convergence variables.
  real(dp) :: energy_difference
  real(dp) :: force_difference
  
  ! Working variables.
  type(String) :: dir
  integer      :: converged_step
  integer      :: final_step
  logical      :: converged
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  ! --------------------------------------------------
  ! Get settings from user, and check them.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  seedname = arguments%value('seedname')
  run_script = arguments%value('run_script')
  no_cores = int(arguments%value('no_cores'))
  minimum_cutoff = int(arguments%value('minimum_cutoff'))
  cutoff_step = int(arguments%value('cutoff_step'))
  maximum_cutoff = int(arguments%value('maximum_cutoff'))
  minimum_kpoints = int(arguments%value('minimum_kpoints'))
  kpoints_step = int(arguments%value('kpoints_step'))
  maximum_kpoints = int(arguments%value('maximum_kpoints'))
  energy_convergence_threshold = &
     & dble(arguments%value('energy_convergence_threshold'))
  force_convergence_threshold = &
     & dble(arguments%value('force_convergence_threshold'))
  no_converged_calculations = int(arguments%value('no_converged_calculations'))
  symmetry_precision = dble(arguments%value('symmetry_precision'))
  
  if (cutoff_step<10) then
    call print_line('Error: cutoff_step is too small.')
    stop
  elseif (maximum_cutoff<=minimum_cutoff) then
    call print_line('Error: maximum_cutoff is smaller than minimum_cutoff.')
    stop
  elseif (kpoints_step<1) then
    call print_line('Error: kpoints_step is too small.')
    stop
  elseif (maximum_kpoints<=minimum_kpoints) then
    call print_line('Error: maximum_kpoints is smaller than minimum_kpoints.')
    stop
  elseif (energy_convergence_threshold<=0.0_dp) then
    call print_line('Error: energy_convergence_threshold must be positive.')
    stop
  elseif (force_convergence_threshold<=0.0_dp) then
    call print_line('Error: force_convergence_threshold must be positive.')
    stop
  elseif (no_converged_calculations<2) then
    call print_line('Error: no_converged_calculations must be at least 2.')
    stop
  endif
  
  ! --------------------------------------------------
  ! Initialise calculation handlers.
  ! --------------------------------------------------
  calculation_writer = CalculationWriter( working_directory = wd,            &
                                        & file_type         = str('castep'), &
                                        & seedname          = seedname       )
  
  calculation_runner = CalculationRunner( working_directory = wd,            &
                                        & file_type         = str('castep'), &
                                        & seedname          = seedname,      &
                                        & run_script        = run_script,    &
                                        & no_cores          = no_cores,      &
                                        & calculation_type  = str('none')    )
  
  calculation_reader = CalculationReader(wd)
  
  ! --------------------------------------------------
  ! Read .cell and .param files.
  ! --------------------------------------------------
  cell_file = IFile(wd//'/'//seedname//'.cell')
  param_file = IFile(wd//'/'//seedname//'.param')
  
  structure = input_file_to_StructureData( str('castep'),             &
                                         & wd//'/'//seedname//'.cell' )
  
  recip_lattice = dble(structure%recip_lattice)
  average_reciprocal_length = 0.0_dp
  do i=1,3
    average_reciprocal_length = average_reciprocal_length        &
                            & + l2_norm(vec(recip_lattice(i,:))) &
                            & / 3.0_dp
  enddo
  
  ! --------------------------------------------------
  ! Calculate the number of steps in both convergence calculations.
  ! --------------------------------------------------
  no_cutoffs = (maximum_cutoff-minimum_cutoff)/cutoff_step &
           & + no_converged_calculations
  
  no_kpoints = (maximum_kpoints-minimum_kpoints)/kpoints_step &
           & + no_converged_calculations
  
  ! --------------------------------------------------
  ! Allocate arrays.
  ! --------------------------------------------------
  allocate( cutoffs(no_cutoffs),                            &
          & kpoint_spacings(no_kpoints),                    &
          & kpoint_numbers(no_kpoints),                     &
          & kpoint_grids(3,no_kpoints),                     &
          & cutoffs_energies(no_cutoffs),                   &
          & kpoints_energies(no_kpoints),                   &
          & cutoffs_forces(structure%no_atoms, no_cutoffs), &
          & kpoints_forces(structure%no_atoms, no_kpoints), &
          & stat=ialloc); call err(ialloc)
  
  ! --------------------------------------------------
  ! Run cutoff convergence.
  ! --------------------------------------------------
  
  progress_file = OFile(wd//'/convergence_progress.dat')
  call progress_file%print_line( 'Energy cutoff convergence:')
  converged = .false.
  do i=1,no_cutoffs
    cutoffs(i) = minimum_cutoff + (i-1)*cutoff_step
    
    dir = wd//'/cutoff_'//cutoffs(i)
    castep_file = run_castep( cutoffs(i),                                &
                            & average_reciprocal_length/minimum_kpoints, &
                            & wd,                                        &
                            & dir,                                       &
                            & seedname,                                  &
                            & run_script,                                &
                            & no_cores,                                  &
                            & cell_file,                                 &
                            & param_file,                                &
                            & structure,                                 &
                            & symmetry_precision)
    cutoffs_energies(i) = castep_file%energy
    cutoffs_forces(:,i) = castep_file%forces%vectors
    
    if (i==1) then
      call progress_file%print_line( 'k-point grid : '// &
                                   & castep_file%kpoints_mp_grid)
    endif
    
    call progress_file%print_line('')
    call progress_file%print_line( 'Cutoff : '// &
                               & cutoffs(i)//' eV')
    call progress_file%print_line( 'Energy : '// &
                               & castep_file%energy//' eV')
    
    if (i>=no_converged_calculations) then
      converged_step = i-no_converged_calculations+1
      energy_difference = maxval(abs( cutoffs_energies(converged_step:i) &
                                  & - cutoffs_energies(i)))
      force_difference = 0
      do j=converged_step,i
        force_difference = max( force_difference,                   &
                              & maxval(l2_norm( cutoffs_forces(:,j) &
                              &               - cutoffs_forces(:,i))))
      enddo
      
      call progress_file%print_line( '')
      call progress_file%print_line( 'Cutoff               : '// &
                                 & cutoffs(converged_step)//' eV')
      call progress_file%print_line( 'Energy difference    : '// &
                                 & energy_difference//' eV')
      call progress_file%print_line( 'Force difference     : '// &
                                 & force_difference//' eV/A')
      
      if ( energy_difference<energy_convergence_threshold .and. &
         & force_difference<force_convergence_threshold ) then
        final_step = i
        converged = .true.
        exit
      endif
    endif
  enddo
  
  call progress_file%print_line( '')
  
  if (converged) then
    call progress_file%print_line( 'Final energy cutoff : '// &
                               & cutoffs(converged_step)//' eV')
    call progress_file%print_line( 'Final energy        : '// &
                               & cutoffs_energies(converged_step)//' eV')
  else
    call print_line('Error: cutoff energy convergence not acheived.')
    stop
  endif
  
  call progress_file%print_line( '')
  
  ! Write out energy convergence file.
  cutoff_file = OFile(wd//'/cutoff_convergence.dat')
  call cutoff_file%print_line( 'Final Energy Cutoff : '// &
                             & cutoffs(converged_step)//' eV')
  call cutoff_file%print_line( 'Energy Tolerance    : '// &
                             & energy_convergence_threshold//' eV/cell')
  call cutoff_file%print_line( 'Force Tolerance     : '// &
                             & force_convergence_threshold//' eV/A')
  call cutoff_file%print_line( '')
  call cutoff_file%print_line( 'Energy cutoff (eV), &
     &energy difference (eV/cell), max force difference (eV/A)')
  do i=1,final_step
    energy_difference = abs(cutoffs_energies(i)-cutoffs_energies(final_step))
    force_difference = maxval(l2_norm( cutoffs_forces(:,i) &
                                   & - cutoffs_forces(:,final_step)))
    call cutoff_file%print_line( cutoffs(i)        //' '// &
                               & energy_difference //' '// &
                               & force_difference)
  enddo
  
  ! --------------------------------------------------
  ! Run k-point convergence.
  ! --------------------------------------------------
  call progress_file%print_line( 'k-point convergence:')
  call progress_file%print_line( 'Cutoff : '// &
                             & minimum_cutoff//' eV')
  converged = .false.
  do i=1,no_kpoints
    kpoint_spacings(i) = average_reciprocal_length &
                     & / (minimum_kpoints + (i-1)*kpoints_step)
    dir = wd//'/kpoints_'//minimum_kpoints+(i-1)*kpoints_step
    castep_file = run_castep( minimum_cutoff,     &
                            & kpoint_spacings(i), &
                            & wd,                 &
                            & dir,                &
                            & seedname,           &
                            & run_script,         &
                            & no_cores,           &
                            & cell_file,          &
                            & param_file,         &
                            & structure,          &
                            & symmetry_precision)
    kpoint_numbers(i) = castep_file%no_kpoints
    kpoint_grids(:,i) = castep_file%kpoints_mp_grid
    kpoints_energies(i) = castep_file%energy
    kpoints_forces(:,i) = castep_file%forces%vectors
    
    call progress_file%print_line( '')
    call progress_file%print_line( 'k-point spacing : '// &
                               & kpoint_spacings(i)//' 1/bohr')
    call progress_file%print_line( 'no. k-points    : '// &
                               & kpoint_numbers(i))
    call progress_file%print_line( 'k-point grid    : '// &
                               & kpoint_grids(:,i))
    call progress_file%print_line( 'Energy          : '// &
                               & castep_file%energy//' eV')
    
    if (i>=no_converged_calculations) then
      converged_step = i-no_converged_calculations+1
      energy_difference = maxval(abs( kpoints_energies(converged_step:i) &
                                  & - kpoints_energies(i)))
      force_difference = 0
      do j=converged_step,i
        force_difference = max( force_difference,                   &
                              & maxval(l2_norm( kpoints_forces(:,j) &
                              &               - kpoints_forces(:,i))))
      enddo
      
      call progress_file%print_line( '')
      call progress_file%print_line( 'k-point spacing       : '// &
                                 & kpoint_spacings(converged_step)//' 1/bohr')
      call progress_file%print_line( 'no. k-points          : '// &
                                 & kpoint_numbers(converged_step))
      call progress_file%print_line( 'k-point grid          : '// &
                                 & kpoint_grids(:,converged_step))
      call progress_file%print_line( 'Energy difference     : '// &
                                 & energy_difference//' eV')
      call progress_file%print_line( 'Force difference      : '// &
                                 & force_difference//' eV/A')
      
      if ( energy_difference<energy_convergence_threshold .and. &
         & force_difference<force_convergence_threshold ) then
        final_step = i
        converged = .true.
        exit
      endif
    endif
  enddo
  
  if (converged) then
    call progress_file%print_line( '')
    call progress_file%print_line( 'Final k-point spacing : '// &
                               & kpoint_spacings(converged_step)//' 1/A')
    call progress_file%print_line( 'Final no. k-points    : '// &
                               & kpoint_numbers(converged_step))
    call progress_file%print_line( 'Final k-point grid    : '// &
                               & kpoint_grids(:,converged_step))
    call progress_file%print_line( 'Final energy          : '// &
                               & kpoints_energies(converged_step)//' eV')
  else
    call print_line('Error: k-point convergence not acheived.')
    stop
  endif
  
  ! Write out k-points file.
  kpoints_file = OFile(wd//'/kpoints_convergence.dat')
  call kpoints_file%print_line( 'Final k-point spacing : '// &
                              & kpoint_spacings(converged_step)//' 1/A')
  call kpoints_file%print_line( 'Final no. k-points    : '// &
                              & kpoint_numbers(converged_step))
  call kpoints_file%print_line( 'Final k-point grid    : '// &
                              & kpoint_grids(:,converged_step))
  call kpoints_file%print_line( 'Energy Tolerance      : '// &
                              & energy_convergence_threshold//' eV/cell')
  call kpoints_file%print_line( 'Force Tolerance       : '// &
                              & force_convergence_threshold//' eV/A')
  call kpoints_file%print_line( '')
  call kpoints_file%print_line('k-point spacing (1/bohr), no. k-points, &
     &k-point grid, energy difference (eV), max force difference (eV/A)')
  do i=1,final_step
    energy_difference = abs(kpoints_energies(i)-kpoints_energies(final_step))
    force_difference = maxval(l2_norm( kpoints_forces(:,i) &
                                   & - kpoints_forces(:,final_step)))
    call kpoints_file%print_line( kpoint_spacings(i) //' '// &
                                & kpoint_numbers(i)  //' '// &
                                & kpoint_grids(:,i)  //' '// &
                                & energy_difference  //' '// &
                                & force_difference)
  enddo
end subroutine

! ----------------------------------------------------------------------
! Runs CASTEP with a particular k-point spacing and cutoff.
! ----------------------------------------------------------------------
function run_castep(cutoff,kpoint_spacing,wd,dir,seedname,run_script, &
   & no_cores,cell_file,param_file,structure,symmetry_precision) result(output)
  implicit none
  
  integer,             intent(in) :: cutoff
  real(dp),            intent(in) :: kpoint_spacing
  type(String),        intent(in) :: wd
  type(String),        intent(in) :: dir
  type(String),        intent(in) :: seedname
  type(String),        intent(in) :: run_script
  integer,             intent(in) :: no_cores
  type(IFile),         intent(in) :: cell_file
  type(IFile),         intent(in) :: param_file
  type(StructureData), intent(in) :: structure
  real(dp),            intent(in) :: symmetry_precision
  type(CastepOutputFile)          :: output
  
  type(OFile)  :: output_cell_file
  type(OFile)  :: output_param_file
  integer      :: result_code
  
  integer :: i
  
  ! Write CASTEP inputs.
  call mkdir(dir)
  
  output_cell_file = OFile(dir//'/'//seedname//'.cell')
  do i=1,size(cell_file)
    call output_cell_file%print_line(cell_file%line(i))
  enddo
  call output_cell_file%print_line( 'kpoints_mp_spacing : '// &
                                  & kpoint_spacing//' 1/bohr')
  
  output_param_file = OFile(dir//'/'//seedname//'.param')
  do i=1,size(param_file)
    call output_param_file%print_line(param_file%line(i))
  enddo
  call output_param_file%print_line('cut_off_energy = '//cutoff//' eV')
  
  ! Run CASTEP.
  call print_line('')
  call print_line('Running calculation in directory '//dir)
  result_code = system_call( 'cd '//wd//'; '//run_script//' castep '//dir// &
     & ' '//no_cores//' '//seedname)
  call print_line('Result code: '//result_code)
  
  ! Read CASTEP file.
  output = read_castep_output_file( dir//'/'//seedname//'.castep', &
                                  & structure,                     &
                                  & wd,                            &
                                  & seedname                       )
end function
end module
