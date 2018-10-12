! ======================================================================
! Converges harmonic frequencies and free energies w/r/t
!    CASTEP cutoff energy and k-point spacing.
! Runs multiple harmonic caesar calculations with different
!    cutoff energies and k-point spacings.
! ======================================================================
module converge_harmonic_frequencies_module
  use common_module

  implicit none
  
  private
  
  public :: converge_harmonic_frequencies
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext
! ----------------------------------------------------------------------
function converge_harmonic_frequencies() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'converge_harmonic_frequencies'
  output%description = 'Converges harmonic frequencies w/r/t cutoff energy &
     &and k-point spacing.'
  output%keywords = [                                                         &
     & KeywordData( 'file_type',                                              &
     &              'file_type is the file type which will be used for &
     &single-point energy calculations. Usual settings are: "castep", &
     &"caesar" and "xyz", but only CASTEP is supported by this convergence &
     &module at present',                                                     &
     &              default_value='castep'),                                  &
     & KeywordData( 'seedname',                                               &
     &               'seedname is the CASTEP seedname from which filenames &
     &are constructed.'),                                                     &
     & KeywordData( 'minimum_cutoff',                                         &
     &              'minimum_cutoff is the smallest cutoff energy which will &
     &be tested. minimum_cutoff must be an integer. This cut-off energy will &
     &be used when converging harmonic frequencies and free energies with &
     &respect to kpoint sampling',                                            &
     &               default_value='300'),                                    &
     & KeywordData( 'cutoff_step',                                            &
     &               'cutoff_step is the step between each cutoff energy &
     &which will be tested. cutoff_step must be an integer.',                 &
     &               default_value='100'),                                    &
     & KeywordData( 'maximum_cutoff',                                         &
     &              'maximum_cutoff is the cutoff energy at which calculation &
     &will be terminated if convergence has not been reached. maximum_cutoff &
     &must be an integer.',                                                   &
     &               default_value='1200'),                                   &
     & KeywordData( 'maximum_kpoint_spacing',                                 &
     &              'maximum_kpoint_spacing is the largest k-point spacing &
     &which will be tested. This is the kpoint spacing that will be used when &
     &converging harmonic frequencies and free energies with respect to &
     &cut-off energy',                                                        &
     &               default_value='0.07'),                                   &
     & KeywordData( 'kpoint_spacing_step',                                    &
     &              'kpoint_spacing_step is the step between each kpoint &
     &spacing which will be tested.',                                         &
     &               default_value='0.01'),                                   &
     & KeywordData( 'minimum_kpoint_spacing',                                 &
     &              'minimum_kpoint_spacing is the kpoint spacing at which &
     &calculation will be terminated if convergence has not been reached.',   &
     &               default_value='0.01'),                                   &
     & KeywordData( 'q-point_grid',                                           &
     &              'q-point_grid is the number of q-points in each direction &
     &in a Monkhorst-Pack grid. This should be specified as three integers &
     &separated by spaces.'),                                                 &
     & KeywordData( 'symmetry_precision',                                     &
     &              'In order for a symmetry to be accepted, it must &
     &transform the position of every atom to within symmetry_precision of an &
     &atom of the same element. symmetry_precision should be given in Bohr.', &
     &              default_value='0.1'),                                     &
     & KeywordData( 'harmonic_displacement',                                  &
     &              'harmonic_displacement is the distance in bohr by which &
     &atoms will be displaced when mapping the harmonic Born-Oppenheimer &
     &surface.',                                                              &
     &              default_value='0.01'),                                    &
     & KeywordData( 'run_script',                                             &
     &              'run_script is the path to the script for running DFT. An &
     &example run script can be found in doc/input_files. This script should &
     &have executable file permissions and should run CASTEP',                &
     &              is_path=.true.),                                          &
     & KeywordData( 'no_cores',                                               &
     &              'no_cores is the number of cores on which the electronic &
     &structure calculation will be run. This is passed to the specified run &
     &script.',                                                               &
     &              default_value='1'),                                       &
     & KeywordData( 'no_nodes',                                               &
     &              'no_nodes is the number of nodes on which the electronic &
     &structure calculation will be run. This is passed to the specified run &
     &script.',                                                               &
     &              default_value='1'),                                      &
     & KeywordData( 'converge_energies',                                      &
     &              'converge_energies dictates whether you want to monitor &
     &the convergence of the free energies as well as the harmonic &
     &frequencies. Options are "yes" and "no".',                              &
     &              default_value='no'),                                    &
     & KeywordData( 'convergence_mode',                                      &
     &              'convergence_mode determines which CASTEP parameters you &
     &wish to vary in your convergence testing.                              &
     &Options are "cutoff", "kpoints" and "both".',                          &
     &              default_value='both'),                                  &
     & KeywordData( 'acoustic_sum_rule',                                      &
     &              'acoustic_sum_rule specifies where the acoustic sum rule  &
     &is applied. The options are "off", "forces", "matrices" and "both".',   &
     &              default_value='both'),                                    &
     & KeywordData( 'min_temperature',                                        &
     &              'min_temperature is the minimum temperature at which &
     &thermodynamic quantities are calculated. min_temperature should be &
     &given in Kelvin.',                                                     &
     &              default_value='0'),                                       &
     & KeywordData( 'max_temperature',                                        &
     &              'max_temperature is the maximum temperature at which &
     &thermodynamic quantities are calculated. min_temperature should be &
     &given in Kelvin.',                                                     &
     &              default_value='500'),                                     &
     & KeywordData( 'no_temperature_steps',                                   &
     &              'no_temperature_steps is the number of temperatures at &
     &which thermodynamic quantities are calculated.',                        &
     &              default_value='6'),                                       &
     & KeywordData( 'min_frequency',                                          &
     &              'min_frequency is the frequency below which modes will be &
     &ignored when calculating thermodynamic quantities. min_frequency should &
     &be given in Hartree.',                                                  &
     &              default_value='1e-8'),                                    &
     & KeywordData( 'path',                                                   &
     &              'path is the path through fractional reciprocal space &
     &which will be mapped by the phonon dispersion curve. The path should be &
     &specified as labels and q-points, separated by commas. The Gamma-point &
     &should be labelled G.',                                                 &
     &              default_value='G 0.0 0.0 0.0, R 0.5 0.5 0.5, &
     &M 0.0 0.5 0.5, G 0.0 0.0 0.0, X 0.0 0.0 0.5'),                          &
     & KeywordData( 'no_dos_samples',                                         &
     &              'no_dos_samples is the number of points in reciprocal &
     &space at which the normal modes are calculated when calculating the &
     &vibrational density of states.',                                        &
     &              default_value='100000'),                                  &
     & KeywordData( 'pseudopotential_file',                                   &
     &              'pseudopotential_file is the name of the pseudopotential &
     &to be used in the CASTEP calculations. If no pseudpotential is needed &
     &(eg. if you are using OTF generation) the value of this should be &
     &"none".',                                                               &
     &              default_value='none'),                                    &
     & KeywordData( 'freq_tolerance',                                         &
     &              'freq_tolerance is the tolerance used in the convergence &
     &of the harmonic frequencies.',                                          &
     &              default_value='1e-8'),                                    &
     & KeywordData( 'energy_tolerance',                                       &
     &              'energy_tolerance is the tolerance used in the &
     &convergence of the vibrational free energies. It is given in Hartree &
     &per unit cell',                                                         &
     &              default_value='1e-6'),                                    &
     & KeywordData( 'convergence_count',                                      &
     &              'convergence_count is the number of consecutive &
     &frequencies and/or energies that must be within the defined tolerance &
     &of each other for convergence to be considered reached',                &
     &              default_value='3')                                        ]
  output%main_subroutine => converge_harmonic_frequencies_subroutine
end function

! ----------------------------------------------------------------------
! Main program
! ----------------------------------------------------------------------
subroutine converge_harmonic_frequencies_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! User inputs
  type(String) :: wd
  type(String) :: seedname
  type(String) :: file_type
  integer      :: minimum_cutoff, cutoff_step, maximum_cutoff
  real(dp)     :: maximum_kpoint_spacing, kpoint_spacing_step, minimum_kpoint_spacing
  integer      :: grid(3)
  real(dp)     :: symmetry_precision
  real(dp)     :: harmonic_displacement
  integer      :: no_cores, no_nodes
  type(String) :: run_script
  type(String) :: converge_energies
  type(String) :: convergence_mode
  type(String) :: acoustic_sum_rule
  real(dp)     :: min_temperature
  real(dp)     :: max_temperature
  integer      :: no_temperature_steps
  real(dp)     :: min_frequency
  type(String) :: path
  integer      :: no_dos_samples
  type(String) :: pseudopotential_file
  real(dp)     :: freq_tolerance, energy_tolerance
  integer      :: convergence_count
  
  ! Files
  type(IFile)            :: cell_file, param_file, supercell_file, pseudopot_file_in
  type(OFile)            :: output_cell_file, output_param_file
  type(OFile)            :: setup_harmonic_inputs, run_harmonic_inputs
  type(OFile)            :: normal_mode_inputs, harm_obs_input, normal_modes, pseudopot_file_out
  type(IFile)            :: complex_mode_file, therm_vars, freqfile1, freqfile2, enfile1, enfile2
  type(OFile)            :: vib_free_en, freq_conv, energy_conv

  ! Energy cut-offs and kpoint spacings
  integer                       :: no_cutoffs, no_kpoint_spacings
  integer,          allocatable :: cutoffs(:), kpoint_spacings(:)
  
  ! Working variables
  type(String) :: dir, qpoint_dir, copyfile, rmfile
  real(dp)     :: max_val
  
  ! Temporary variables
  integer      :: i, j, k, h, ialloc, no_qpoints, result_code, conv_count1, conv_count2
  type(String) :: line
  
  ! --------------------------------------------------
  ! Get settings from user, and check them.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  file_type = arguments%value('file_type')
  seedname = arguments%value('seedname')
  minimum_cutoff = int(arguments%value('minimum_cutoff'))
  cutoff_step = int(arguments%value('cutoff_step'))
  maximum_cutoff = int(arguments%value('maximum_cutoff'))
  maximum_kpoint_spacing = dble(arguments%value('maximum_kpoint_spacing'))
  kpoint_spacing_step = dble(arguments%value('kpoint_spacing_step'))
  minimum_kpoint_spacing = dble(arguments%value('minimum_kpoint_spacing'))
  grid = int(split_line(arguments%value('q-point_grid')))
  symmetry_precision = dble(arguments%value('symmetry_precision'))
  harmonic_displacement = dble(arguments%value('harmonic_displacement'))
  no_cores = int(arguments%value('no_cores'))
  no_nodes = int(arguments%value('no_nodes'))
  run_script = arguments%value('run_script')
  converge_energies = arguments%value('converge_energies')
  convergence_mode = arguments%value('convergence_mode')
  acoustic_sum_rule = arguments%value('acoustic_sum_rule')
  min_temperature = dble(arguments%value('min_temperature'))
  max_temperature = dble(arguments%value('max_temperature'))
  no_temperature_steps = int(arguments%value('no_temperature_steps'))
  min_frequency = dble(arguments%value('min_frequency'))
  path = arguments%value('path')
  no_dos_samples = int(arguments%value('no_dos_samples'))
  pseudopotential_file = arguments%value('pseudopotential_file')
  freq_tolerance = dble(arguments%value('freq_tolerance'))
  energy_tolerance = dble(arguments%value('energy_tolerance'))
  convergence_count = int(arguments%value('convergence_count'))
  
  if (cutoff_step<10) then
    call print_line('Error: cutoff_step is too small.')
    stop
  elseif (maximum_cutoff<=minimum_cutoff) then
    call print_line('Error: maximum_cutoff is smaller than minimum_cutoff.')
    stop
  elseif (maximum_kpoint_spacing<=minimum_kpoint_spacing) then
    call print_line('Error: maximum_kpoints is smaller than minimum_kpoints.')
    stop
  elseif (minimum_kpoint_spacing<=0.004_dp) then
    call print_line('Error: minimum_kpoint_spacing is too small.')
    stop
  elseif (maximum_kpoint_spacing>=0.09_dp) then
    call print_line('Error: maximum_kpoint_spacing is too large.')
    stop
  elseif (kpoint_spacing_step<=0.004_dp) then
    call print_line('Error: kpoint_spacing_step is too small.')
    stop
  elseif (kpoint_spacing_step<=0.004_dp) then
    call print_line('Error: kpoint_spacing_step is too small.')
    stop
  elseif (file_type/='castep') then
    call print_line('Error: CASTEP is the only accepted file type for this mode.')
    stop
  elseif ((converge_energies/='yes').AND.(converge_energies/='no')) then
    call print_line('Error: converge_energies must be "yes" or "no".')
    stop
  elseif ((convergence_mode/='cutoff').AND.(convergence_mode/='kpoints').AND.(convergence_mode/='both')) then
    call print_line('Error: convergence_mode must be "cutoff", "kpoints" or "both".')
    stop
  elseif (max_temperature<=min_temperature) then
    call print_line('Error: max_temperature is smaller than min_temperature.')
    stop
  elseif (convergence_count<1) then
    call print_line('Error: convergence_count must be at least 1.')
    stop
  endif
  
  ! --------------------------------------------------
  ! Read .cell and .param files
  ! --------------------------------------------------
  cell_file = IFile(seedname//'.cell')
  param_file = IFile(seedname//'.param')
  if (pseudopotential_file/='none') then
      pseudopot_file_in = IFile(pseudopotential_file)
  endif
  
  ! ---------------------------------------------------------------------------
  ! Calculate number of steps in convergence calculations and number of qpoints
  ! ---------------------------------------------------------------------------
  no_cutoffs = (maximum_cutoff-minimum_cutoff)/cutoff_step &
           & + 1
  
  no_kpoint_spacings = nint((maximum_kpoint_spacing-minimum_kpoint_spacing)/kpoint_spacing_step) &
           & + 1

  no_qpoints = grid(1)*grid(2)*grid(3)
  
  ! --------------------------------------------------
  ! Allocate arrays
  ! --------------------------------------------------
  allocate( cutoffs(no_cutoffs),                            &
          & kpoint_spacings(no_kpoint_spacings),            &
          & stat=ialloc); call err(ialloc)
  
  ! ---------------------------------------------------------
  ! Create input files for cutoff convergence and run caesar
  ! ---------------------------------------------------------
  freq_conv = OFile('frequency_conv')
  energy_conv = OFile('energy_conv')
  conv_count1 = 0
  conv_count2 = 0

  if ((convergence_mode=='cutoff').OR.(convergence_mode=='both')) then
    call freq_conv%print_line( 'Result for kpoint spacing: '//maximum_kpoint_spacing)
    call energy_conv%print_line( 'Results for kpoint spacing: '//maximum_kpoint_spacing)
    do i=1,no_cutoffs
      cutoffs(i) = minimum_cutoff + (i-1)*cutoff_step
      dir = 'cutoff_'//cutoffs(i)

! Write CASTEP cell and param files
      output_cell_file = write_castep_cell( nint(maximum_kpoint_spacing*1000),   &
                                        & dir, seedname, cell_file)

      output_param_file = write_castep_param( cutoffs(i), dir, seedname, param_file)

      if (pseudopotential_file/='none') then
         pseudopot_file_out = OFile(dir//'/'//pseudopotential_file)
         do j=1,size(pseudopot_file_in)
              call pseudopot_file_out%print_line( pseudopot_file_in%line(j))
         enddo
      endif
         
! Write setup_harmonic input file and call caesar in setup_harmonic mode
      setup_harmonic_inputs = write_setup( file_type, seedname, grid, &
                                       & symmetry_precision, harmonic_displacement, dir)

      call call_caesar('setup_harmonic -f setup_inputs -d '//dir)

! Write run_harmonic input files and call caesar in run_harmonic mode
      supercell_file = IFile(dir//'/no_supercells.dat')
      run_harmonic_inputs = write_runinput( supercell_file, run_script, no_cores, no_nodes, dir)
  
      call call_caesar('run_harmonic -f run_inputs -d '//dir)

      normal_mode_inputs=write_normmode( acoustic_sum_rule, dir)
      call call_caesar('calculate_normal_modes -f normal_mode_inputs -d '//dir)

      normal_modes = OFile(dir//'/normal_mode_freq')

      do j=1,no_qpoints
        if (no_qpoints<10) then
           qpoint_dir=dir//'/qpoint_'//j
        elseif (no_qpoints>=10 .AND. no_qpoints<100) then
           if (j<10) then
              qpoint_dir=dir//'/qpoint_0'//j
           else
              qpoint_dir=dir//'/qpoint_'//j
           endif
        elseif (no_qpoints>=100 .AND. no_qpoints<1000) then
           if (j<10) then
              qpoint_dir=dir//'/qpoint_00'//j
           elseif (j>=10 .AND. j<100) then
              qpoint_dir=dir//'/qpoint_0'//j
           else
              qpoint_dir=dir//'/qpoint_'//j
           endif
        else 
           call print_line('Consider reducing vibrational grid size for convergence testing.')
           stop
        endif

        complex_mode_file = IFile(qpoint_dir//'/complex_modes.dat')
        do k=1,size(complex_mode_file)
           h = findindex_string(complex_mode_file%line(k),str('Mode frequency'))
           if (h>=1) then
              line = complex_mode_file%line(k)
              line = slice(line,30,len(line))
              call normal_modes%print_line(line)
           endif
        enddo
      enddo

      if (converge_energies=='yes') then
        harm_obs_input=write_obs_input( min_temperature, max_temperature, no_temperature_steps, &
                                         & min_frequency, path, no_dos_samples, dir)
        call call_caesar('calculate_harmonic_observables -f harm_obs_input -d '//dir)

        therm_vars = IFile(dir//'/harmonic_observables/thermodynamic_variables.dat')

        vib_free_en = OFile(dir//'/free_energies')

        do k=2,size(therm_vars)
           line = therm_vars%line(k)
           line = slice(line,55,77)
           call vib_free_en%print_line(line)
        enddo
      endif

      if (i>1) then
         copyfile='cp '//dir//'/normal_mode_freq '//dir//'/normal_mode_freq2'
         result_code=system_call(copyfile)
         if (result_code/=0) then
            call print_line(ERROR//': Copying file failed.')
         endif
         freqfile1=IFile(dir//'/normal_mode_freq2')
         copyfile='cp cutoff_'//cutoffs(i-1)//'/normal_mode_freq cutoff_'//cutoffs(i-1)//'/normal_mode_freq2'
         result_code=system_call(copyfile)
         if (result_code/=0) then
            call print_line(ERROR//': Coping file failed.')
         endif
         freqfile2=IFile('cutoff_'//cutoffs(i-1)//'/normal_mode_freq2')
         max_val = maxval(abs(dble(freqfile1%lines())-dble(freqfile2%lines())))
         call freq_conv%print_line( 'Cutoff energy changed from '//cutoffs(i-1)//' to '//cutoffs(i)//'eV')
         call freq_conv%print_line( 'Maximum frequency change: '//max_val)

         if (max_val<freq_tolerance) then
            conv_count1 = conv_count1 + 1
         endif
         call freq_conv%print_line( 'Convergence count: '//conv_count1)
         if (conv_count1>=convergence_count) then
            call freq_conv%print_line( ' ')
            call freq_conv%print_line( 'Desired harmonic phonon frequency convergence reached at '//cutoffs(i)//'eV')
         endif
         call freq_conv%print_line( ' ')

         rmfile='rm '//dir//'/normal_mode_freq2'
         result_code=system_call(rmfile)
         if (result_code/=0) then
            call print_line(ERROR//': Removing file failed.')
         endif
         rmfile='rm cutoff_'//cutoffs(i-1)//'/normal_mode_freq2'
         result_code=system_call(rmfile)
         if (result_code/=0) then
            call print_line(ERROR//': Removing file failed.')
         endif

         if (converge_energies=='yes') then

            copyfile='cp '//dir//'/free_energies '//dir//'/free_energies2'
            result_code=system_call(copyfile)
            if (result_code/=0) then
               call print_line(ERROR//': Copying file failed.')
            endif
            enfile1=IFile(dir//'/free_energies2')
            copyfile='cp cutoff_'//cutoffs(i-1)//'/free_energies cutoff_'//cutoffs(i-1)//'/free_energies2'
            result_code=system_call(copyfile)
            if (result_code/=0) then
               call print_line(ERROR//': Coping file failed.')
            endif
            enfile2=IFile('cutoff_'//cutoffs(i-1)//'/free_energies2')
            max_val = maxval([(                                    &
               & abs(dble(enfile1%line(j))-dble(enfile2%line(j))), &
               & j=1,                                              &
               & size(therm_vars)-1                                )])
            call energy_conv%print_line( 'Cutoff energy changed from '//cutoffs(i-1)//' to '//cutoffs(i)//'eV')
            call energy_conv%print_line( 'Maximum free energy change: '//max_val//' Hartree per cell')

            if (max_val<energy_tolerance) then
               conv_count2 = conv_count2 + 1
            endif
            call energy_conv%print_line( 'Convergence count: '//conv_count2)
            if (conv_count2>=convergence_count)then
               call energy_conv%print_line( ' ')
               call energy_conv%print_line( 'Desired vibrational free energy convergence reached at '//cutoffs(i)//'eV')
            endif
            call energy_conv%print_line( ' ')

            rmfile='rm '//dir//'/free_energies2'
            result_code=system_call(rmfile)
            if (result_code/=0) then
               call print_line(ERROR//': Removing file failed.')
            endif
            rmfile='rm cutoff_'//cutoffs(i-1)//'/free_energies2'
            result_code=system_call(rmfile)
            if (result_code/=0) then
               call print_line(ERROR//': Removing file failed.')
            endif

         endif

      endif

    enddo

    if (conv_count1<convergence_count) then
       call freq_conv%print_line( ' ')
       call freq_conv%print_line( 'Desired convergence of harmonic frequencies wrt cut-off energy not achieved')
       call freq_conv%print_line( ' ')
    endif
    if ((converge_energies=='yes').AND.(conv_count2<convergence_count)) then
       call energy_conv%print_line( ' ')
       call energy_conv%print_line( 'Desired convergence of vibrational free energy wrt cut-off energy not achieved')
       call energy_conv%print_line( ' ')
    endif

  endif
  
  ! -----------------------------------------------------------------
  ! Create input files for kpoint spacing convergence and run caesar
  ! -----------------------------------------------------------------
  
  conv_count1 = 0
  conv_count2 = 0

  if ((convergence_mode=='kpoints').OR.(convergence_mode=='both')) then
    call freq_conv%print_line( ' ')
    call freq_conv%print_line( 'Result for cut-off energy: '//minimum_cutoff)
    call energy_conv%print_line( ' ')
    call energy_conv%print_line( 'Result for cut-off energy: '//minimum_cutoff)
    do i=1,no_kpoint_spacings
      kpoint_spacings(i) = nint((maximum_kpoint_spacing-(i-1)*kpoint_spacing_step)*1000)

      dir = 'kpoints_0.0'//kpoint_spacings(i)

! Write CASTEP cell and param files
      output_cell_file = write_castep_cell( kpoint_spacings(i),   &
                                        & dir, seedname, cell_file)

      output_param_file = write_castep_param( minimum_cutoff, dir, seedname, param_file)

      if (pseudopotential_file/='none') then
         pseudopot_file_out = OFile(dir//'/'//pseudopotential_file)
         do j=1,size(pseudopot_file_in)
              call pseudopot_file_out%print_line( pseudopot_file_in%line(j))
         enddo
      endif

! Write setup_harmonic input file and call caesar in setup_harmonic mode
      setup_harmonic_inputs = write_setup( file_type, seedname, grid, &
                                       & symmetry_precision, harmonic_displacement, dir)

      call call_caesar('setup_harmonic -f setup_inputs -d '//dir)

! Write run_harmonic input files and call caesar in run_harmonic mode
      supercell_file = IFile(dir//'/no_supercells.dat')
      run_harmonic_inputs = write_runinput( supercell_file, run_script, no_cores, no_nodes, dir)

      call call_caesar('run_harmonic -f run_inputs -d '//dir)

      normal_mode_inputs=write_normmode( acoustic_sum_rule, dir)
      call call_caesar('calculate_normal_modes -f normal_mode_inputs -d '//dir)

      normal_modes = OFile(dir//'/normal_mode_freq')

      do j=1,no_qpoints
        if (no_qpoints<10) then
           qpoint_dir=dir//'/qpoint_'//j
        elseif (no_qpoints>=10 .AND. no_qpoints<100) then
           if (j<10) then
              qpoint_dir=dir//'/qpoint_0'//j
           else
              qpoint_dir=dir//'/qpoint_'//j
           endif
        elseif (no_qpoints>=100 .AND. no_qpoints<1000) then
           if (j<10) then
              qpoint_dir=dir//'/qpoint_00'//j
           elseif (j>=10 .AND. j<100) then
              qpoint_dir=dir//'/qpoint_0'//j
           else
              qpoint_dir=dir//'/qpoint_'//j
           endif
        else 
           call print_line('Consider reducing vibrational grid size for convergence testing.')
           stop
        endif

        complex_mode_file = IFile(qpoint_dir//'/complex_modes.dat')
        do k=1,size(complex_mode_file)
           h = findindex_string(complex_mode_file%line(k),str('Mode frequency'))
           if (h>=1) then
             line = complex_mode_file%line(k)
             line = slice(line,30,len(line))
             call normal_modes%print_line(line)
           endif

        enddo
      enddo

      if (converge_energies=='yes') then
        harm_obs_input=write_obs_input( min_temperature, max_temperature, no_temperature_steps, &
                                         & min_frequency, path, no_dos_samples, dir)
        call call_caesar('calculate_harmonic_observables -f harm_obs_input -d '//dir)

        therm_vars = IFile(dir//'/harmonic_observables/thermodynamic_variables.dat')

        vib_free_en = OFile(dir//'/free_energies')

        do k=2,size(therm_vars)
          line = therm_vars%line(k)
          line = slice(line,55,77)
          call vib_free_en%print_line(line)
        enddo
      endif

      if (i>1) then
         copyfile='cp '//dir//'/normal_mode_freq '//dir//'/normal_mode_freq2'
         result_code=system_call(copyfile)
         if (result_code/=0) then
            call print_line(ERROR//': Copying file failed.')
         endif
         freqfile1=IFile(dir//'/normal_mode_freq2')
         copyfile='cp kpoints_0.0'//kpoint_spacings(i-1)//'/normal_mode_freq kpoints_0.0'&
                                              &//kpoint_spacings(i-1)//'/normal_mode_freq2'
         result_code=system_call(copyfile)
         if (result_code/=0) then
            call print_line(ERROR//': Coping file failed.')
         endif
         freqfile2=IFile('kpoints_0.0'//kpoint_spacings(i-1)//'/normal_mode_freq2')
         max_val = maxval(abs(dble(freqfile1%lines())-dble(freqfile2%lines())))
         call freq_conv%print_line( 'Kpoint spacing changing from 0.0'//kpoint_spacings(i-1)//' to 0.0'//kpoint_spacings(i))
         call freq_conv%print_line( 'Maximum frequency change: '//max_val)

         if (max_val<freq_tolerance) then
            conv_count1 = conv_count1 + 1
         endif
         call freq_conv%print_line( 'Convergence count: '//conv_count1)
         if (conv_count1>=convergence_count) then
            call freq_conv%print_line( ' ')
            call freq_conv%print_line( 'Desired harmonic phonon frequency convergence reached at 0.0'&
                                                              &//kpoint_spacings(i)//' kpoint spacing')
         endif
         call freq_conv%print_line( ' ')

         rmfile='rm '//dir//'/normal_mode_freq2'
         result_code=system_call(rmfile)
         if (result_code/=0) then
            call print_line(ERROR//': Removing file failed.')
         endif
         rmfile='rm kpoints_0.0'//kpoint_spacings(i-1)//'/normal_mode_freq2'
         result_code=system_call(rmfile)
         if (result_code/=0) then
            call print_line(ERROR//': Removing file failed.')
         endif

         if (converge_energies=='yes') then
            copyfile='cp '//dir//'/free_energies '//dir//'/free_energies2'
            result_code=system_call(copyfile)
            if (result_code/=0) then
               call print_line(ERROR//': Copying file failed.')
            endif
            enfile1=IFile(dir//'/free_energies2')

            copyfile='cp kpoints_0.0'//kpoint_spacings(i-1)//'/free_energies kpoints_0.0'&
                                               &//kpoint_spacings(i-1)//'/free_energies2'
            result_code=system_call(copyfile)
            if (result_code/=0) then
               call print_line(ERROR//': Coping file failed.')
            endif
            enfile2=IFile('kpoints_0.0'//kpoint_spacings(i-1)//'/free_energies2')
            max_val = maxval([(                                    &
               & abs(dble(enfile1%line(j))-dble(enfile2%line(j))), &
               & j=1,                                              &
               & size(therm_vars)-1                                )])
            call energy_conv%print_line( 'Kpoint spacing changing from 0.0'//kpoint_spacings(i-1)//' to 0.0'//kpoint_spacings(i))
            call energy_conv%print_line( 'Maximum free energy change: '//max_val//' Hartree per cell')

            if (max_val<energy_tolerance) then
               conv_count2 = conv_count2 + 1
            endif
            call energy_conv%print_line( 'Convergence count: '//conv_count2)
            if (conv_count2>=convergence_count) then
               call energy_conv%print_line( ' ')
               call energy_conv%print_line( 'Desired vibrational free energy convergence reached at 0.0'&
                                                                         &//kpoint_spacings(i)//' kpoint spacing')
            endif
            call energy_conv%print_line( ' ')

            rmfile='rm '//dir//'/free_energies2'
            result_code=system_call(rmfile)
            if (result_code/=0) then
               call print_line(ERROR//': Removing file failed.')
            endif
            rmfile='rm kpoints_0.0'//kpoint_spacings(i-1)//'/free_energies2'
            result_code=system_call(rmfile)
            if (result_code/=0) then
               call print_line(ERROR//': Removing file failed.')
            endif
         endif

      endif

    enddo

    if (conv_count1<convergence_count) then
       call freq_conv%print_line( ' ')
       call freq_conv%print_line( 'Desired convergence of harmonic frequencies wrt kpoint sampling not achieved')
       call freq_conv%print_line( ' ')
    endif
    if ((converge_energies=='yes').AND.(conv_count2<convergence_count)) then
       call energy_conv%print_line( ' ')
       call energy_conv%print_line( 'Desired convergence of vibrational free energy wrt kpoint sampling not achieved')
       call energy_conv%print_line( ' ')
    endif

  endif

end subroutine

! ----------------------------------------------------------------------
! FUNCTIONS
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Make directory and write CASTEP cell file with particular k-point
! spacing (replaces original kpoint spacing/grid in cell file if present)
! ----------------------------------------------------------------------
function write_castep_cell(kpoint_spacing,dir,seedname, &
   & cell_file) result(output_cell_file)
  implicit none
  
  integer,             intent(in) :: kpoint_spacing
  type(String),        intent(in) :: dir, seedname
  type(IFile),         intent(in) :: cell_file
  
  type(OFile)  :: output_cell_file
  
  integer :: i, j, k, l, h, replace
  
  call mkdir(dir)
  replace = 0
  output_cell_file = OFile(dir//'/'//seedname//'.cell')

  do i=1,size(cell_file)
    j = findindex_string(cell_file%line(i),str('mp_spacing'))
    k = findindex_string(cell_file%line(i),str('mp_grid'))
    l = findindex_string(cell_file%line(i),str('MP_SPACING'))
    h = findindex_string(cell_file%line(i),str('MP_GRID'))
! Replace if kpoint line, otherwise print original line
    if ((j>=1).OR.(k>=1).OR.(l>=1).OR.(h>=1)) then
       call output_cell_file%print_line( 'kpoints_mp_spacing 0.0'//kpoint_spacing)
       replace = 1
    else
       call output_cell_file%print_line(cell_file%line(i))
    endif
  enddo

! Print kpoint spacing if not present in original file
  if (replace == 0) then
     call output_cell_file%print_line( 'kpoints_mp_spacing 0.0'//kpoint_spacing)
  endif

end function

! ----------------------------------------------------------------------
! Write CASTEP param file with a particular cutoff
! (replaces original cut-off energy in param file if present)
! ----------------------------------------------------------------------
function write_castep_param(cutoff,dir,seedname, &
   & param_file) result(output_param_file)
  implicit none
  
  integer,             intent(in) :: cutoff
  type(String),        intent(in) :: dir, seedname
  type(IFile),         intent(in) :: param_file
  
  type(OFile)  :: output_param_file
  
  integer :: i, j, k, replace

  replace = 0
  output_param_file = OFile(dir//'/'//seedname//'.param')

  do i=1,size(param_file)
    j = findindex_string(param_file%line(i),str('cut_off_energy'))
    k = findindex_string(param_file%line(i),str('CUT_OFF_ENERGY'))
! Replace if cut-off energy line, otherwise print original line
    if ((j>=1).OR.(k>=1)) then
       call output_param_file%print_line( 'cut_off_energy    : '//cutoff)
       replace = 1
    else
       call output_param_file%print_line(param_file%line(i))
    endif
  enddo

! Write energy cut-off if not present in original file
  if (replace == 0) then
     call output_param_file%print_line( 'cut_off_energy    : '//cutoff)
  endif
  
end function

! ----------------------------------------------------------------------
! Write setup_harmonic input file in specified directory
! ----------------------------------------------------------------------
function write_setup(file_type, seedname, grid, symmetry_precision, &
   & harmonic_displacement, dir) result(setup_harmonic_inputs)
  implicit none
  
  integer,             intent(in) :: grid(3)
  type(String),        intent(in) :: dir, seedname, file_type
  real(dp),            intent(in) :: symmetry_precision, harmonic_displacement
  
  type(OFile)  :: setup_harmonic_inputs
  
  setup_harmonic_inputs = OFile(dir//'/setup_inputs')
  call setup_harmonic_inputs%print_line( 'file_type '//file_type)
  call setup_harmonic_inputs%print_line( 'seedname '//seedname)
  call setup_harmonic_inputs%print_line( 'q-point_grid '//grid)
  call setup_harmonic_inputs%print_line( 'symmetry_precision '//symmetry_precision)
  call setup_harmonic_inputs%print_line( 'harmonic_displacement '//harmonic_displacement)
  
end function

! ----------------------------------------------------------------------
! Write run_harmonic input file in specified directory
! ----------------------------------------------------------------------
function write_runinput( supercell_file, run_script, no_cores, no_nodes, dir) result(run_harmonic_inputs)
  implicit none
  
  type(IFile),         intent(in) :: supercell_file
  type(String),        intent(in) :: dir
  type(String),        intent(in) :: run_script
  integer,             intent(in) :: no_cores, no_nodes

  integer :: i

  type(OFile)  :: run_harmonic_inputs
  
  run_harmonic_inputs = OFile(dir//'/run_inputs')
  do i=1,size(supercell_file)
    call run_harmonic_inputs%print_line( 'supercells_to_run 1 '//supercell_file%line(i))
  enddo
  call run_harmonic_inputs%print_line( 'run_script '//run_script)
  call run_harmonic_inputs%print_line( 'no_cores '//no_cores)
  call run_harmonic_inputs%print_line( 'no_nodes '//no_nodes)
  
end function

! ----------------------------------------------------------------------
! Index function for two strings
! ----------------------------------------------------------------------
function findindex_character(string1,string2) result(index_found)
  implicit none
  character(*), intent(in)     :: string1, string2
  integer      :: index_found
  
  index_found=index(string1,string2)
end function

function findindex_string(string1,string2) result(index_found)
  implicit none
  type(String),        intent(in) :: string1, string2
  integer      :: index_found
  
  index_found = findindex_character(char(string1),char(string2))
end function

! ----------------------------------------------------------------------
! Write calculate_normal_modes input file in specified directory
! ----------------------------------------------------------------------
function write_normmode( acoustic_sum_rule, dir) result(normal_mode_inputs)
  implicit none
  type(String),        intent(in) :: dir, acoustic_sum_rule

  type(OFile)  :: normal_mode_inputs
  
  normal_mode_inputs = OFile(dir//'/normal_mode_inputs')

  call normal_mode_inputs%print_line( 'acoustic_sum_rule '//acoustic_sum_rule)
  
end function

! ----------------------------------------------------------------------
! Write calculate_harmonic_observables input file in specified directory
! ----------------------------------------------------------------------
function write_obs_input( min_temperature, max_temperature, no_temperature_steps, &
                         & min_frequency, path, no_dos_samples, dir) result(harm_obs_input)
  implicit none

  real(dp),            intent(in) :: min_temperature
  real(dp),            intent(in) :: max_temperature
  integer,             intent(in) :: no_temperature_steps
  real(dp),            intent(in) :: min_frequency
  type(String),        intent(in) :: path
  integer,             intent(in) :: no_dos_samples
  type(String),        intent(in) :: dir

  type(OFile)  :: harm_obs_input
  
  harm_obs_input = OFile(dir//'/harm_obs_input')

  call harm_obs_input%print_line( 'min_temperature '//min_temperature)
  call harm_obs_input%print_line( 'max_temperature '//max_temperature)
  call harm_obs_input%print_line( 'no_temperature_steps '//no_temperature_steps)
  call harm_obs_input%print_line( 'min_frequency '//min_frequency)
  call harm_obs_input%print_line( 'path '//path)
  call harm_obs_input%print_line( 'no_dos_samples '//no_dos_samples)
  
end function

end module
