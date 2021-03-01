submodule (caesar_converge_qpoint_grid_module) caesar_converge_qpoint_grid_submodule
  use caesar_harmonic_module
contains

module procedure converge_qpoint_grid
  output%mode_name = 'converge_qpoint_grid'
  output%description = 'Monitors convergence of harmonic frequencies and free energies wrt qpoint grid.'
  output%keywords = [                                                         &
  & KeywordData( 'file_type',                                                 &
  &              'file_type is the file type which will be used for &
  &single-point energy calculations. Usual settings are: "castep", "caesar" and &
  &"xyz" - note this mode only works for CASTEP inputs at present',           &
  &              default_value='castep'),                                     &
  & KeywordData( 'seedname',                                                 &
  &               'seedname is the CASTEP seedname from which filenames are &
  &constructed.'),                                                            &
  & KeywordData( 'qgrid_start',                                              &
  &              'qgrid_start is the starting (minimum) number of q-points in each &
  &direction in a Monkhorst-Pack grid. The number of points in each direction &
  &will be the same, so this should be specified as a single integer.'),    &
  & KeywordData( 'qgrid_end',                                              &
  &              'qgrid_end is the final (maximum) number of q-points in each direction &
  &in a Monkhorst-Pack grid. The number of points in each direction will be the &
  &same, so this should be specified as a single integer.'),                  &
  & KeywordData( 'symmetry_precision',                                        &
  &              'In order for a symmetry to be accepted, it must transform &
  &the position of every atom to within symmetry_precision of an atom of the &
  &same element. symmetry_precision should be given in Bohr.',                &
  &              default_value='0.1'),                                        &
  & KeywordData( 'harmonic_displacement',                                     &
  &              'harmonic_displacement is the distance in bohr by which &
  &atoms will be displaced when mapping the harmonic Born-Oppenheimer &
  &surface.',                                                                 &
  &              default_value='0.01'),                                    &
  & KeywordData( 'run_script',                                             &
  &              'run_script is the path to the script for running DFT. An &
  &example run script can be found in doc/input_files. This script should   &
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
  &              'converge_energies determines whether convergence of the free &
  &energies is monitored as well as convergence of harmonic frequencies.     &
  &Options are "yes" and "no".',                                            &
  &              default_value='no'),                                    &
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
  &              default_value='500'),                                       &
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
  &              default_value='100000'),                                &
  & KeywordData( 'pseudopotential_file',                                 &
  &              'pseudopotential_file is the name of the pseudopotential &
  &to be used in the CASTEP calculations. If no pseudpotential is needed &
  &(eg. if you are using OTF generation) the value of this should be "none".', &
  &              default_value='none'),                                  &
  & KeywordData( 'freq_tolerance',                                     &
  &              'freq_tolerance is the tolerance used in the convergence &
  &of the harmonic frequencies.',                                       &
  &              default_value='1e-8'),                                &
  & KeywordData( 'energy_tolerance',                                     &
  &              'energy_tolerance is the tolerance used in the convergence &
  &of the vibrational free energies. It is given in Hartree per unit cell',   &
  &              default_value='1e-6'),                                 &
  & KeywordData( 'convergence_count',                                     &
  &              'convergence_count is the number of consecutive frequencies &
  &and/or energies that must be within the defined tolerance of each other &
  &for convergence to be considered reached',   &
  &              default_value='3')]
  output%main_subroutine => converge_qpoint_grid_subroutine
end procedure

module procedure converge_qpoint_grid_subroutine
  ! User inputs
  type(String) :: wd
  type(String) :: seedname
  type(String) :: file_type
  integer      :: qgrid_start, qgrid_end
  real(dp)     :: symmetry_precision
  real(dp)     :: harmonic_displacement
  integer      :: no_cores, no_nodes
  type(String) :: run_script
  type(String) :: converge_energies
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
  type(IFile)            :: supercell_file
  type(OFile)            :: setup_harmonic_inputs, run_harmonic_inputs
  type(OFile)            :: normal_mode_inputs, harm_obs_input, normal_modes
  type(IFile)            :: complex_mode_file, therm_vars, freqfile1, freqfile2, enfile1, enfile2
  type(OFile)            :: vib_free_en, freq_conv, energy_conv
  
  ! Working variables
  type(String) :: dir, qpoint_dir, grid, copyfile, rmfile
  real(dp)     :: var1, var2, var_diff, max_val
  
  ! Temporary variables
  integer      :: i, j, k, h, no_qpoints, result_code, conv_count1, conv_count2
  integer      :: freqconv_check, enconv_check
  character (len = 100) :: line_to_trim, var1char, var2char
  
  ! --------------------------------------------------
  ! Get settings from user, and check them.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  file_type = arguments%value('file_type')
  seedname = arguments%value('seedname')
  qgrid_start = int(arguments%value('qgrid_start'))
  qgrid_end = int(arguments%value('qgrid_end'))
  symmetry_precision = dble(arguments%value('symmetry_precision'))
  harmonic_displacement = dble(arguments%value('harmonic_displacement'))
  no_cores = int(arguments%value('no_cores'))
  no_nodes = int(arguments%value('no_nodes'))
  run_script = arguments%value('run_script')
  converge_energies = arguments%value('converge_energies')
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
  
  if (file_type/='castep') then
    call print_line('Error: CASTEP is the only accepted file type for this mode.')
    call quit()
  elseif ((converge_energies/='yes').AND.(converge_energies/='no')) then
    call print_line('Error: converge_energies must be "yes" or "no".')
    call quit()
  elseif (convergence_count<1) then
    call print_line('Error: convergence_count must be at least 1.')
    call quit()
  elseif (qgrid_start>=qgrid_end) then
    call print_line('Error: qgrid_start must be less than qgrid_end.')
    call quit()
  endif

  ! ---------------------------------------------------------
  ! Create input files for cutoff convergence and run caesar
  ! ---------------------------------------------------------
  freq_conv = OFile('frequency_conv')
  energy_conv = OFile('energy_conv')
  conv_count1 = 0
  conv_count2 = 0
  freqconv_check = 0
  if (converge_energies=='yes') then
     enconv_check = 0
  else 
     enconv_check = 1
  endif

    do i=qgrid_start,qgrid_end
      if ((freqconv_check==0).OR.(enconv_check==0)) then
        dir = 'grid_'//i
        call mkdir(dir)

        result_code = system_call('cp *.cell '//dir)
        if (result_code/=0) then
           call print_line(ERROR//': Error copying .cell file.')
           call err()
        endif

        result_code = system_call('cp *.param '//dir)
        if (result_code/=0) then
           call print_line(ERROR//': Error copying .param file.')
           call err()
        endif

        if (pseudopotential_file/='none') then
           result_code = system_call('cp '//pseudopotential_file//' '//dir)
           if (result_code/=0) then
              call print_line(ERROR//': Error copying .cell file.')
              call err()
           endif
        endif
         
         ! Write setup_harmonic input file and call caesar in setup_harmonic mode
         grid = i//' '//i//' '//i
         setup_harmonic_inputs = write_setupgrid( file_type, seedname, grid, &
                                       & symmetry_precision, harmonic_displacement, dir)

         call call_caesar('setup_harmonic -f setup_inputs -d '//dir)

         ! Write run_harmonic input files and call caesar in run_harmonic mode
         supercell_file = IFile(dir//'/no_supercells.dat')
         run_harmonic_inputs = write_runinput( supercell_file, run_script, no_cores, no_nodes, dir)
  
         call call_caesar('run_harmonic -f run_inputs -d '//dir)

         normal_mode_inputs=write_normmode( acoustic_sum_rule, dir)
         call call_caesar('calculate_normal_modes -f normal_mode_inputs -d '//dir)

         normal_modes = OFile(dir//'/normal_mode_freq')

         no_qpoints = i**3

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
              call print_line('Very large vibrational grid.')
              call quit()
           endif

           complex_mode_file = IFile(qpoint_dir//'/complex_modes.dat')
           do k=1,size(complex_mode_file)
              h = findindex_string(complex_mode_file%line(k),str('Mode frequency'))
              if (h>=1) then
                 line_to_trim = char(complex_mode_file%line(k))
                 call normal_modes%print_line( trim(line_to_trim(30:)))
              endif
           enddo
         enddo

         if ((converge_energies=='yes').AND.(enconv_check==0)) then
           harm_obs_input=write_obs_input( min_temperature, max_temperature, no_temperature_steps, &
                                         & min_frequency, path, no_dos_samples, dir)
           call call_caesar('calculate_harmonic_observables -f harm_obs_input -d '//dir)

           therm_vars = IFile(dir//'/harmonic_observables/thermodynamic_variables.dat')

           vib_free_en = OFile(dir//'/free_energies')

           do k=2,size(therm_vars)
              line_to_trim = char(therm_vars%line(k))
              call vib_free_en%print_line( trim(line_to_trim(55:77)))
           enddo
         endif

         if (i>qgrid_start) then
            copyfile='cp '//dir//'/normal_mode_freq '//dir//'/normal_mode_freq2'
            result_code=system_call(copyfile)
            if (result_code/=0) then
               call print_line(ERROR//': Copying file failed.')
            endif
            freqfile1=IFile(dir//'/normal_mode_freq2')
            copyfile='cp grid_'//(i-1)//'/normal_mode_freq grid_'//(i-1)//'/normal_mode_freq2'
            result_code=system_call(copyfile)
            if (result_code/=0) then
               call print_line(ERROR//': Coping file failed.')
            endif
            freqfile2=IFile('grid_'//(i-1)//'/normal_mode_freq2')
            max_val=0.0_dp
            do j=1,size(freqfile1)
               var1char=char(freqfile1%line(j))
               var2char=char(freqfile2%line(j))
               read(var1char,*) var1
               read(var2char,*) var2
               var_diff=abs(var1-var2)
               if (var_diff>max_val) then
                  max_val=var_diff
               endif
            enddo

            call freq_conv%print_line( 'qpoint grid changed from '// &
                                     & (i-1)//' '//(i-1)//' '//(i-1)// &
                                     & ' to '//i//' '//i//' '//i)
            call freq_conv%print_line( 'Maximum frequency change: '//max_val)

            if (max_val<freq_tolerance) then
               conv_count1 = conv_count1 + 1
            else 
               conv_count1 = 0
            endif
            call freq_conv%print_line( 'Convergence count: '//conv_count1)
            if (conv_count1>=convergence_count) then
               call freq_conv%print_line( ' ')
               call freq_conv%print_line( 'Desired frequency convergence reached at qpoint grid of '//i//' '//i//' '//i)
            endif
            call freq_conv%print_line( ' ')

            rmfile='rm '//dir//'/normal_mode_freq2'
            result_code=system_call(rmfile)
            if (result_code/=0) then
               call print_line(ERROR//': Removing file failed.')
            endif
            rmfile='rm grid_'//(i-1)//'/normal_mode_freq2'
            result_code=system_call(rmfile)
            if (result_code/=0) then
               call print_line(ERROR//': Removing file failed.')
            endif

            if ((converge_energies=='yes').AND.(enconv_check==0)) then

               copyfile='cp '//dir//'/free_energies '//dir//'/free_energies2'
               result_code=system_call(copyfile)
               if (result_code/=0) then
                  call print_line(ERROR//': Copying file failed.')
               endif
               enfile1=IFile(dir//'/free_energies2')
               copyfile='cp grid_'//(i-1)//'/free_energies grid_'//(i-1)//'/free_energies2'
               result_code=system_call(copyfile)
               if (result_code/=0) then
                  call print_line(ERROR//': Coping file failed.')
               endif
               enfile2=IFile('grid_'//(i-1)//'/free_energies2')
               max_val=0.0_dp
               do j=1,(size(therm_vars)-1)
                  var1char=char(enfile1%line(j))
                  var2char=char(enfile2%line(j))
                  read(var1char,*) var1
                  read(var2char,*) var2
                  var_diff=abs(var1-var2)
                  if (var_diff>max_val) then
                     max_val=var_diff
                  endif
               enddo

               call energy_conv%print_line( 'qpoint grid changed from '//   &
                                          & (i-1)//' '//(i-1)//' '//(i-1)// &
                                          & ' to '//i//' '//i//' '//i       )
               call energy_conv%print_line( 'Maximum free energy change: '// &
                                          & max_val//' Hartree per cell')

               if (max_val<energy_tolerance) then
                  conv_count2 = conv_count2 + 1
               else
                  conv_count2 = 0
               endif
               call energy_conv%print_line( 'Convergence count: '//conv_count2)
               if (conv_count2>=convergence_count)then
                  call energy_conv%print_line( ' ')
                  call energy_conv%print_line( 'Free energy convergence reached at qpoint grid of '//i//' '//i//' '//i)
               endif
               call energy_conv%print_line( ' ')

               rmfile='rm '//dir//'/free_energies2'
               result_code=system_call(rmfile)
               if (result_code/=0) then
                  call print_line(ERROR//': Removing file failed.')
               endif
               rmfile='rm grid_'//(i-1)//'/free_energies2'
               result_code=system_call(rmfile)
               if (result_code/=0) then
                  call print_line(ERROR//': Removing file failed.')
               endif

            endif

         endif

      endif

      if (conv_count1>=convergence_count) then
         freqconv_check = 1
      endif
      if ((converge_energies=='yes').AND.(conv_count2>=convergence_count)) then
         enconv_check = 1
      endif

    enddo

    if (conv_count1<convergence_count) then
       call freq_conv%print_line( ' ')
       call freq_conv%print_line( 'Desired convergence of harmonic frequencies wrt qpoint grid not achieved')
       call freq_conv%print_line( ' ')
    endif
    if ((converge_energies=='yes').AND.(conv_count2<convergence_count)) then
       call energy_conv%print_line( ' ')
       call energy_conv%print_line( 'Desired convergence of vibrational free energy wrt qpoint grid not achieved')
       call energy_conv%print_line( ' ')
    endif

end procedure

module procedure write_setupgrid
  setup_harmonic_inputs = OFile(dir//'/setup_inputs')
  call setup_harmonic_inputs%print_line( 'file_type '//file_type)
  call setup_harmonic_inputs%print_line( 'seedname '//seedname)
  call setup_harmonic_inputs%print_line( 'q-point_grid '//grid)
  call setup_harmonic_inputs%print_line( 'symmetry_precision '//symmetry_precision)
  call setup_harmonic_inputs%print_line( 'harmonic_displacement '//harmonic_displacement)
  
end procedure

module procedure write_runinput
  integer :: i
  
  run_harmonic_inputs = OFile(dir//'/run_inputs')
  do i=1,size(supercell_file)
    call run_harmonic_inputs%print_line( 'supercells_to_run 1 '//supercell_file%line(i))
  enddo
  call run_harmonic_inputs%print_line( 'run_script '//run_script)
  call run_harmonic_inputs%print_line( 'no_cores '//no_cores)
  call run_harmonic_inputs%print_line( 'no_nodes '//no_nodes)
  
end procedure

module procedure write_normmode
  normal_mode_inputs = OFile(dir//'/normal_mode_inputs')

  call normal_mode_inputs%print_line( 'acoustic_sum_rule '//acoustic_sum_rule)
  
end procedure

module procedure write_obs_input
  harm_obs_input = OFile(dir//'/harm_obs_input')

  call harm_obs_input%print_line( 'min_temperature '//min_temperature)
  call harm_obs_input%print_line( 'max_temperature '//max_temperature)
  call harm_obs_input%print_line( 'no_temperature_steps '//no_temperature_steps)
  call harm_obs_input%print_line( 'min_frequency '//min_frequency)
  call harm_obs_input%print_line( 'path '//path)
  call harm_obs_input%print_line( 'no_dos_samples '//no_dos_samples)
  
end procedure

module procedure findindex_character
  index_found=index(string1,string2)
end procedure

module procedure findindex_string
  index_found = findindex_character(char(string1),char(string2))
end procedure
end submodule
