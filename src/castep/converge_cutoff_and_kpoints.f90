! ======================================================================
! Converges the cutoff energy.
! ======================================================================
! Runs multiple single-point calculations with different cutoff energies.
module converge_cutoff_and_kpoints_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function converge_cutoff_and_kpoints_keywords() result(keywords)
  use keyword_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                                &
  & make_keyword( 'seedname',                                                 &
  &               'seedname is the CASTEP seedname from which file names are &
  &constructed.'),                                                            &
  & make_keyword( 'run_script',                                               &
  &               'run_script is the path to the script for running CASTEP. &
  &An example run script can be found in doc/input_files.',                   &
  &               is_path=.true.),                                            &
  & make_keyword( 'no_cores',                                                 &
  &               'no_cores is the number of cores on which CASTEP will be &
  &run. This is passed to the specified run script.',                         &
  &               default_value='1'),                                         &
  & make_keyword( 'minimum_cutoff',                                           &
  &               'minimum_cutoff is the smallest cutoff energy which will be &
  &tested. minimum_cutoff should be an integer.',                             &
  &               default_value='300'),                                       &
  & make_keyword( 'cutoff_step',                                              &
  &               'cutoff_step is the step between each cutoff energy which &
  &will be tested. cutoff_step should be an integer.',                        &
  &               default_value='50'),                                        &
  & make_keyword( 'maximum_cutoff',                                           &
  &               'maximum_cutoff is the largest cutoff energy which will be &
  &tested. Caesar will terminate with an error if convergence has not been &
  &acheived by the maximum cutoff energy. maximum_cutoff should be an &
  &integer.',                                                                 &
  &               default_value='1500'),                                      &
  & make_keyword( 'maximum_spacing',                                          &
  &               'maximum_spacing is the largest k-point spacing which will &
  &be tested.',                                                    &
  &               default_value='0.1'),                                       &
  & make_keyword( 'spacing_step',                                             &
  &               'spacing_step is the step between each k-point spacing &
  &which will be tested.',                                                    &
  &               default_value='0.005'),                                     &
  & make_keyword( 'minimum_spacing',                                          &
  &               'minimum_spacing is the smallest k-point spacing which will &
  &be tested. Caesar will terminate with an error if convergence has not been &
  &acheived by the minimum k-point spacing.',                                 &
  &               default_value='0.01'),                                      &
  & make_keyword( 'convergence_threshold',                                    &
  &               'convergence_threshold is the energy to which the &
  &calculations must converge in order for the cutoff energy to be &
  &accepted.'),                                                               &
  & make_keyword( 'no_converged_calculations',                                &
  &               'no_converged_calculations is the number of calculations &
  &which must be within convergence_threshold in order to accept the cutoff &
  &energy. The cutoff energy of the first such calculation will be &
  &accepted.',                                                                &
  &               default_value='3'),                                         &
  & make_keyword( 'generate_plots',                                           &
  &               'generate_plots specifies that further calculations should &
  &be run to allow plots of the convergence testing to be made. This will &
  &take extra time.',                                                         &
  &               default_value='true') ]
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine converge_cutoff_and_kpoints(arguments)
  use dictionary_module
  use dft_output_file_module
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
  real(dp)     :: maximum_spacing
  real(dp)     :: spacing_step
  real(dp)     :: minimum_spacing
  real(dp)     :: convergence_threshold
  integer      :: no_converged_calculations
  logical      :: generate_plots
  
  ! Files.
  type(String), allocatable :: cell_file(:)
  type(String), allocatable :: param_file(:)
  integer                   :: output_file
  
  ! Arrays.
  integer,  allocatable :: cutoffs(:)
  real(dp), allocatable :: spacings(:)
  real(dp), allocatable :: energies(:)
  real(dp), allocatable :: finer_cutoff_energies(:)
  real(dp), allocatable :: finer_spacing_energies(:)
  
  ! Working variables.
  integer  :: max_no_steps
  integer  :: step
  integer  :: converged_step
  logical  :: converged
  integer  :: no_cutoffs
  integer  :: no_spacings
  
  ! Output variables.
  integer  :: cutoff
  real(dp) :: spacing
  real(dp) :: energy
  
  ! Temporary variables.
  integer      :: i,ialloc
  
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
  maximum_spacing = dble(arguments%value('maximum_spacing'))
  spacing_step = dble(arguments%value('spacing_step'))
  minimum_spacing = dble(arguments%value('minimum_spacing'))
  convergence_threshold = dble(arguments%value('convergence_threshold'))
  no_converged_calculations = int(arguments%value('no_converged_calculations'))
  generate_plots = lgcl(arguments%value('generate_plots'))
  
  if (maximum_cutoff<=minimum_cutoff) then
    call print_line('Error: maximum_cutoff must be greater than &
       &minimum_cutoff')
    stop
  elseif (cutoff_step<10) then
    call print_line('Error: cutoff_step is too small.')
    stop
  elseif (maximum_spacing<=minimum_spacing) then
    call print_line('Error: maximum_spacing must be greater than &
       &minimum_spacing')
    stop
  elseif (spacing_step<0.0001_dp) then
    call print_line('Error: spacing_step is too small.')
    stop
  elseif (convergence_threshold<=0.0_dp) then
    call print_line('Error: convergence_threshold must be positive.')
    stop
  elseif (no_converged_calculations<2) then
    call print_line('Error: no_converged_calculations must be at least 2.')
    stop
  endif
  
  ! Calculate the maximum number of steps.
  ! The +2 is to allow for the first step at coarsest values, and to prevent
  !    overflow when calculating the final step.
  max_no_steps = (maximum_cutoff-minimum_cutoff)/cutoff_step           &
             & + floor((maximum_spacing-minimum_spacing)/spacing_step) &
             & + 2
  
  ! --------------------------------------------------
  ! Read .cell and .param files.
  ! --------------------------------------------------
  cell_file = read_lines(wd//'/'//seedname//'.cell')
  param_file = read_lines(wd//'/'//seedname//'.param')
  
  ! --------------------------------------------------
  ! Run convergence tests.
  ! --------------------------------------------------
  allocate( cutoffs(max_no_steps),                &
          & spacings(max_no_steps),               &
          & energies(max_no_steps),               &
          & finer_cutoff_energies(max_no_steps),  &
          & finer_spacing_energies(max_no_steps), &
          & stat=ialloc); call err(ialloc)
  converged = .false.
  cutoffs(1) = minimum_cutoff
  spacings(1) = maximum_spacing
  energies(1) = calculate_energy( cutoffs(1),  &
                                & spacings(1), &
                                & wd,          &
                                & seedname,    &
                                & run_script,  &
                                & no_cores,    &
                                & cell_file,   &
                                & param_file)
  do step=1,max_no_steps
    ! Check for convergence.
    if (step>=no_converged_calculations) then
      converged_step = step-no_converged_calculations+1
      if (all( abs( energies(converged_step:step) &
           &      - energies(converged_step))     &
           & < convergence_threshold)) then
        converged = .true.
        exit
      endif
    endif
    
    if (step==max_no_steps) then
      call print_line('Code Error: max_no_steps reached without convergence.')
      call err()
    endif
    
    ! Calculate energy at finer energy cutoff and k-point spacing.
    finer_cutoff_energies(step) = calculate_energy( cutoffs(1)+cutoff_step, &
                                                  & spacings(1),            &
                                                  & wd,                     &
                                                  & seedname,               &
                                                  & run_script,             &
                                                  & no_cores,               &
                                                  & cell_file,              &
                                                  & param_file)
    
    finer_spacing_energies(step) = calculate_energy( cutoffs(1),              &
                                                   & spacings(1)-spacing_step,&
                                                   & wd,                      &
                                                   & seedname,                &
                                                   & run_script,              &
                                                   & no_cores,                &
                                                   & cell_file,               &
                                                   & param_file)
    
    ! Copy whichever change made the largest difference to the next step.
    if ( abs(finer_cutoff_energies(step)-energies(step)) >= &
       & abs(finer_spacing_energies(step)-energies(step)) ) then
      cutoffs(step+1)  = cutoffs(step)+cutoff_step
      spacings(step+1) = spacings(step)
      energies(step+1) = finer_cutoff_energies(step)
    else
      cutoffs(step+1)  = cutoffs(step)
      spacings(step+1) = spacings(step)-spacing_step
      energies(step+1) = finer_spacing_energies(step)
    endif
    
    ! Check that neither cutoff nor spacing has exceeded its bounds.
    if (cutoffs(step+1)>maximum_cutoff) then
      converged_step = step
      call print_line('Error: convergence not acheived within the maximum &
         &energy cutoff.')
      exit
    elseif (spacings(step+1)<minimum_spacing) then
      converged_step = step
      call print_line('Error: convergence not acheived within the minimum &
         &k-point spacing.')
      exit
    endif
  enddo
  
  ! --------------------------------------------------
  ! Write out progress if no convergence.
  ! --------------------------------------------------
  if (.not. converged) then
    output_file = open_write_file(wd//'/convergence_progress.dat')
    call print_line(output_file, 'cutoff | spacing | structure energy')
    do i=1,converged_step
      call print_line( output_file, &
                     & cutoffs(i)//' '//spacings(i)//' '//energies(i))
    enddo
    close(output_file)
    stop
  endif
  
  ! --------------------------------------------------
  ! Generate output.
  ! --------------------------------------------------
  output_file = open_write_file(wd//'/convergence_results.dat')
  call print_line(output_file, 'cut_off_energy = '//cutoffs(converged_step))
  call print_line( output_file, &
                 & 'kpoints_mp_spacing : '//spacings(converged_step))
  call print_line(output_file, 'Energy '//energies(converged_step))
  
  if (.not. generate_plots) then
    close(output_file)
    stop
  endif
  
  ! Generate cutoff plot data.
  no_cutoffs = (cutoffs(converged_step)-cutoffs(1))/cutoff_step &
           & + no_converged_calculations
  
  call print_line(output_file, '')
  call print_line(output_file, 'Energy cutoff convergence:')
  call print_line(output_file, 'cutoff energy | structure energy')
  
  spacing = spacings(converged_step)
  do i=1,no_cutoffs
    cutoff = minimum_cutoff + (i-1)*cutoff_step
    energy = calculate_energy( cutoff,     &
                             & spacing,    &
                             & wd,         &
                             & seedname,   &
                             & run_script, &
                             & no_cores,   &
                             & cell_file,  &
                             & param_file)
    call print_line(output_file, cutoff//' '//energy)
  enddo
  
  ! Generate spacing plot data.
  no_spacings = nint( (spacings(1)-spacings(converged_step))/spacing_step ) &
            & + no_converged_calculations
  
  call print_line(output_file, '')
  call print_line(output_file, 'k-point spacing convergence:')
  call print_line(output_file, 'k-point spacing | structure energy')
  
  cutoff = cutoffs(converged_step)
  do i=1,no_spacings
    spacing = maximum_spacing - (i-1)*spacing_step
    energy = calculate_energy( cutoff,     &
                             & spacing,    &
                             & wd,         &
                             & seedname,   &
                             & run_script, &
                             & no_cores,   &
                             & cell_file,  &
                             & param_file)
    call print_line(output_file, spacing//' '//energy)
  enddo
  
  close(output_file)
end subroutine

function calculate_energy(cutoff,spacing,wd,seedname,run_script,no_cores, &
   & cell_file,param_file) result(output)
  use utils_module, only : mkdir
  use dft_output_file_module
  implicit none
  
  integer,          intent(in) :: cutoff
  real(dp),         intent(in) :: spacing
  type(String),     intent(in) :: wd
  type(String),     intent(in) :: seedname
  type(String),     intent(in) :: run_script
  integer,          intent(in) :: no_cores
  type(String),     intent(in) :: cell_file(:)
  type(String),     intent(in) :: param_file(:)
  real(dp)                     :: output
  
  type(String)        :: dir
  character(6)        :: spacing_char
  integer             :: output_cell_file
  integer             :: output_param_file
  integer             :: result_code
  type(DftOutputFile) :: castep_file
  
  integer :: i
  
  ! Convert spacing to a string for the purposes of naming the directory.
  write(spacing_char,'(f1.4)') spacing
  
  ! Write CASTEP inputs.
  dir = wd//'/cutoff_'//cutoff//'_spacing_'//spacing_char
  call mkdir(dir)
  
  output_cell_file = open_write_file(dir//'/'//seedname//'.cell')
  do i=1,size(cell_file)
    call print_line(output_cell_file, cell_file(i))
  enddo
  call print_line( output_cell_file, &
                 & 'kpoints_mp_spacing : '//spacing)
  close(output_cell_file)
  
  output_param_file = open_write_file(dir//'/'//seedname//'.param')
  do i=1,size(param_file)
    call print_line(output_param_file, param_file(i))
  enddo
  call print_line( output_param_file, &
                 & 'cut_off_energy = '//cutoff//' eV')
  close(output_param_file)
  
  ! Run CASTEP.
  call print_line('')
  call print_line('Running calculation in directory '//dir)
  result_code = system_call( 'cd '//wd//'; '//run_script//' castep '//dir// &
     & ' '//no_cores//' '//seedname)
  call print_line('Result code: '//result_code)
  
  ! Read CASTEP energy.
  castep_file = read_dft_output_file( str('castep'), &
                                    & dir//'/'//seedname//'.castep')
  output = castep_file%energy
end function

end module
