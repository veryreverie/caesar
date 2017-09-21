!  Title : BO_MAPPING
!  Author : B. Monserrat
!  Last updated : 29 November 2011
!  Description : This script calculates CASTEP total energies as a function of normal mode amplitude

!  It requires the following input files:
!    1. caesar.input : general input parameters
!    2. disp_patterns.dat : atomic displacement patterns (calculated with 'LTE')
!    3. equilibrium.dat : atomic equilibrium positions and atomic masses
!    4. symmetry.dat : symmetry operations
!    5. CASTEP : .param and .cell input files
!

! Modified by J.C.A. Prentice, 23 June 2015
! Added finding gradient of BO surface to be used in fitting process
! Added in coupled mapping, 09 Feb 2016

module bo_mappingforcescoupled_module
 use string_module
 use io_module
 use constants_module, only : dp
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function anharmonic_keywords() result(keywords)
  use keyword_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                                &
  & make_keyword( 'seed_name',                                                &
  &               'seed_name is the DFT seedname from which file names are &
  &constructed.'),                                                            &
  & make_keyword( 'num_indep_data',                                           &
  &               'num_indep_data is the number of data points per mode. It &
  &should be an odd integer.',                                                &
  &               default_value='11'),                                        &
  & make_keyword( 'temperature', &
  &               'temperature is the temperature in Kelvin at which &
  &thermodynamic quantities are calculated.'),                                &
  & make_keyword( 'first_mode',                                               &
  &               'first_mode is the first mode to be considered.',           &
  &               default_value='4'),                                         &
  & make_keyword( 'last_mode',                                                &
  &               'last_mode is the last mode to be considered.',             &
  &               default_value='4'),                                         &
  & make_keyword( 'first_amp',                                                &
  &               'first_amp is the first amplitude to be considered.',       &
  &               default_value='1'),                                         &
  & make_keyword( 'last_amp',                                                 &
  &               'last_amp is the last amplitude to be considered.',         &
  &               default_keyword='num_indep_data'),                          &
  & make_keyword( 'mc_sampling',                                              &
  &               'mc_sampling should be specified to turn on Monte-Carlo &
  &sampling.',                                                                &
  &               is_boolean=.true.),                                         &
  & make_keyword( 'mc_data_points',                                           &
  &               'mc_data_points is the number of Monte-Carlo data points. &
  &Should only be set if mc_sampling is set.',                                &
  &               default_value='20'),                                        &
  & make_keyword( 'mc_continuation',                                          &
  &               'mc_continuation is the Monte-Carlo continuation. Should &
  &only be set if mc_sampling is set.',                                       &
  &               default_value='20'),                                        &
  & make_keyword( 'coupled_sampling',                                         &
  &               'coupled_sampling should be specified to turn on coupled &
  &sampling.',                                                                &
  &                is_boolean=.true.),                                        &
  & make_keyword( 'num_2body_data',                                           &
  &               'num_2body_data is the number of two-body data points. &
  &Should only be set if coupled_sampling is set.',                           &
  &               default_value='11'),                                        &
  & make_keyword( 'first_amp_2body',                                          &
  &               'first_amp_2body is the first two-body amplitude &
  &considered. Must be <= num_2body_data.',                                   &
  &               default_value='1'),                                         &
  & make_keyword( 'last_amp_2body',                                           &
  &               'last_amp_2body is the last two-body amplitude considered. &
  &Must be >= first_amp_2body and <= num_2body data.',                        &
  &               default_keyword='num_2body_data'),                          &
  & make_keyword( 'magres',                                                   &
  &               'magres specifies whether or not the DFT calculation is &
  &magres.',                                                                  &
  &               is_boolean=.true.) ]
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine bo_mappingforcescoupled(arguments)
  use structure_module
  use dictionary_module
  use positions_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(StructureData) :: supercell
  
  ! ----------------------------------------------------------------------
  ! Read inputs
  ! ----------------------------------------------------------------------
  wd = arguments%value("working_directory")
  
  supercell = read_structure_file('')

  ! Read in CASTEP seed name (required)
  seed_name = arguments%value("seed_name")

  ! Read in number of data points per mode (odd integer) (default=11)
  num_indep_data = int(arguments%value("num_indep_data"))

  ! Read in temperature (K)
  temperature = dble(arguments%value("temperature"))
  if (temperature==0) then
    num_data = 0
  endif

  ! Read in number of unit cells (required)
  no_unit_cells = supercell%sc_size

  ! Read in number of atoms in unit cell (required)
  basis = supercell%no_atoms_prim

  ! Read in first mode to consider (default=4)
  ! TODO: calculate this automatically
  first_mode = int(arguments%value("first_mode"))

  ! Read in last mode to consider (default=4)
  ! TODO: calculate this automatically
  last_mode = int(arguments%value("last_mode"))

  ! Read in first amplitude to consider (default=1) (must be <= num_indep_data)
  ! TODO: calculate this automatically
  first_amp = int(arguments%value("first_amp"))

  ! Read in last amplitude to consider (default=num_indep_data) (must be <= num_indep_data) (first <= last)
  last_amp = int(arguments%value("last_amp"))

  ! Read in sampling method (default=false)
  mc_sampling = exists(arguments, "mc_sampling")
  if (mc_sampling) then
    ! (default=20)
    mc_data_points = int(arguments%value("mc_data_points"))
    ! (default=0)
    mc_continuation = int(arguments%value("mc_continuation"))
  endif

  ! (default=false)
  coupled_sampling = exists(arguments, "coupled_sampling")
  if (coupled_sampling) then
    ! (default=11)
    num_2body_data = int(arguments%value("num_2body_data"))
    
    ! (default=1) ( <= num_2body_data)
    first_amp_2body = int(arguments%value("first_amp_2body"))
    
    ! (default=num_2body_data) ( <= num_2body_data)
    last_amp_2body = int(arguments%value("last_amp_2body"))
  endif

  ! Whether or not CASTEP task = magres
  magres = exists(arguments, "magres")
  if (magres) then
    call  mkdir(wd//'/'//magres)
  endif

  ! ----------------------------------------------------------------------
  ! Create output directories.
  ! ----------------------------------------------------------------------
  call mkdir(wd//'/'//band_structure)
  call mkdir(wd//'/'//total_output)
  call mkdir(wd//'/'//max_amplitude)
  call mkdir(wd//'/'//gradients)

  ! WORKING VARIABLES
  
  ! Declare variable 'freq_line' that stores the line in the
  !    'disp_patterns.dat' file at which the frequency of the working
  !    normal mode is
  freq_line = 1
  
  ! Declare variable 'disp_line' that stores the line in the
  !    'disp_patterns.dat' file at which the atomic displacement we
  !    are working on is
  disp_line = freq_line + 3
  
  ! Initialise variable 'blank_line' to 0
  blank_line_original = 0

  ! ----------------------------------------------------------------------
  !   MAIN PROGRAM
  ! ----------------------------------------------------------------------
  
  disp_patterns_file = read_lines(wd//'/disp_patterns.dat')
  
  output_file = open_write_file(wd//'/caesar.output')

  ! PAA sampling or MC sampling
  if (.not. mc_sampling) then 
    ! Declare variable 'blank_line_original' that finds the next relevant blank line in 'disp_patterns.dat' file
    do i=1,size(disp_patterns_file)
      if (len(disp_patterns_file(i))==0) then
        blank_line_original = i
        exit
      endif
    enddo

    ! Calculate ground state energy
    if (file_exists(wd//'/gs_energy.dat')) then
      gs_energy_file_in = read_lines(wd//'/gs_energy.dat')
      line = split(gs_energy_file_in(1))
      static_energy = dble(line(2))
    else
      call print_line(output_file, 'Ground state calculation')

      call write_dft_input_file(supercell, wd//'/'//seed_name//'.cell')

      call run_dft

      ! Collect partial output information
      energy = dft_output_file%energy ! corrected for finite basis set
      static_energy = energy
      call print_line(output_file, ' Ground state energy (eV): '//energy)

      ! Organise output files
      mv $seed_name.bands band_structure
      
      gs_energy_file_out = open_write_file(wd//'/gs_energy.dat')
      call print_line(gs_energy_file_out, '0.0 '//static_energy)
      close(gs_energy_file_out)
      
      mv $seed_name.castep total_output

      if (magres) then
        mv $seed_name.magres magres
      endif
    endif

    if (coupled_sampling) then
      do n=first_mode,last_mode
        do_m : do m=first_mode,last_mode
          ! Calculate blank_line, freq_line and disp_line for current normal mode
          blank_line1 = n * blank_line_original
          freq_line1 = blank_line1 - blank_line_original + 1
          disp_line1 = freq_line1 + 3 
          gvec_line1 = freq_line1 + 1
          
          blank_line2 = m * blank_line_original
          freq_line2 = blank_line2 - blank_line_original + 1
          disp_line2 = freq_line2 + 3
          gvec_line2 = freq_line2 + 1

          ! Calculate only non-symmetry related modes
          reference1=$(awk 'NR==$n {print $2}' symmetry.dat)
          reference2=$(awk 'NR==$m {print $2}' symmetry.dat)
          if (n == reference1) then
            if (m > n .or. (m < n .and. m /= reference2)) then
              ! Calculate frequency of current normal mode
              frequency1=$(awk 'NR==$freq_line1 {print $3}' disp_patterns.dat)
              frequency2=$(awk 'NR==$freq_line2 {print $3}' disp_patterns.dat)
              
              ! Write frequencies to the output 'coupled_energy.dat' file
              !echo $frequency1 $frequency2 >> coupled_energy.dat
              
              call print_line(output_file, ' Calculating coupling of modes '//n//' & '//m)
              
              ! Loop over the atomic displacements for a particular normal mode
              if (n == 7 .and. m < 10) then
                i_initial = last_amp_2body + 1
              elseif (n == 7 .and. m == 10) then
                i_initial = 10
              else
                i_initial = first_amp_2body
              endif
              
              do i=i_initial,last_amp_2body
                do j=first_amp_2body,last_amp_2body
                  !Run 'positions2.f90' to evaluate the atomic positions for the current normal mode amplitudes
                  do k=disp_line1,blank_line1-1
                    disp1(k-disp_line1+1) = &
                       & vec(dble(split(disp_pattern_file(k))))
                  enddo
                  do k=disp_line2,blank_line2-1
                    disp2(k-disp_line2+1) = &
                       & vec(dble(split(disp_pattern_file(k))))
                  enddo
                  displaced_position =  positions(supercell, temperature, &
                     & num_2body_data, frequency1, frequency2, i, j)
                  scaling_check = displaced_position%scaling_check
                  scaling1 = displaced_position%scaling1
                  scaling2 = displaced_position%scaling2
                  
                  ! Consider only non-translational modes
                  if (scaling_check) then
                    
                    ! Write CASTEP input file 'seed.cell' with calculated atomic positions
                    write_cell(positions.dat) > $seed_name.cell
                    
                    check_static_energy = (num_2body_data - 1) / 2 + 1
                    if (check_static_energy == i .and. check_static_energy == j) then
                      echo 0.0 0.0 $static_energy 0.0 0.0 >> coupled_energy.dat
                      echo '    - Amplitude ' $i ' ' $j ' energy (eV):' $static_energy >> caesar.output
                      
                    elseif (check_static_energy == i .and. num_2body_data == num_indep_data) then
                      cp band_structure/$seed_name.$m.$j.bands ./$seed_name.bands
                      cp total_output/$seed_name.$m.$j.castep ./$seed_name.castep
                      if (magres) then
                        cp magres/$seed_name.$m.$j.magres ./$seed_name.magres
                      endif
                      
                      ! Collect partial output information
                      line_number=$(awk '/Final energy/{x=NR}END{print x}' $seed_name.castep)
                      energy=$(awk 'NR==$line_number {print $9}' $seed_name.castep)
                      read -a gradient <<< $(cat gradient.dat)   ! JCAP
                      echo $scaling1 $scaling2 $energy ${gradient[0]} ${gradient[1]} >> coupled_energy.dat
                      echo '    - Amplitude ' $i ' ' $j ' energy (eV):' $energy >> caesar.output
                      echo '    - Amplitude ' $i ' ' $j ' forces (a.u.):' ${gradient[0]} ' ' ${gradient[1]} >> caesar.output ! JCAP
                      
                      mv $seed_name.bands band_structure
                      mv $seed_name.castep total_output
                      if (magres) then
                        mv $seed_name.magres magres
                      endif
                    elseif (check_static_energy == j .and. num_2body_data == num_indep_data) then
                      cp band_structure/$seed_name.$n.$i.bands ./$seed_name.bands
                      cp total_output/$seed_name.$n.$i.castep ./$seed_name.castep
                      if (magres) then
                        cp magres/$seed_name.$n.$i.magres ./$seed_name.magres
                      endif
                      
                      ! Collect partial output information
                      line_number=$(awk '/Final energy/{x=NR}END{print x}' $seed_name.castep)
                      energy=$(awk 'NR==$line_number {print $9}' $seed_name.castep)
                      read -a gradient <<< $(cat gradient.dat)   ! JCAP
                      echo $scaling1 $scaling2 $energy ${gradient[0]} ${gradient[1]} >> coupled_energy.dat
                      echo '    - Amplitude ' $i ' ' $j ' energy (eV):' $energy >> caesar.output
                      echo '    - Amplitude ' $i ' ' $j ' forces (a.u.):' ${gradient[0]} ' ' ${gradient[1]} >> caesar.output ! JCAP
                      
                      mv $seed_name.bands band_structure
                      mv $seed_name.castep total_output
                      if (magres) then
                        mv $seed_name.magres magres
                      endif
                    else
                      ! Run CASTEP
                      run_castep
                      
                      ! Clean-up after run to remove large output files
                      call rm machine_file
                      cp input_files/$seed_name.* .
                      
                      ! Collect partial output information
                      line_number=$(awk '/Final energy/{x=NR}END{print x}' $seed_name.castep)
                      energy=$(awk 'NR==$line_number {print $9}' $seed_name.castep)
                      read -a gradient <<< $(cat gradient.dat)   ! JCAP
                      echo $scaling1 $scaling2 $energy ${gradient[0]} ${gradient[1]} >> coupled_energy.dat
                      echo '    - Amplitude ' $i ' ' $j ' energy (eV):' $energy >> caesar.output
                      echo '    - Amplitude ' $i ' ' $j ' forces (a.u.):' ${gradient[0]} ' ' ${gradient[1]} >> caesar.output ! JCAP
                      mv gradient.dat gradients/gradient.$n.$i.$m.$j.dat ! JCAP
                      
                      ! Organise output files
                      mv $seed_name.bands band_structure
                      mv $seed_name.castep total_output
                      if (magres) then
                        mv $seed_name.magres magres
                      endif
                    endif
                  else
                    echo ' Skipping translational modes...' >> caesar.output
                    cycle do_m
                  endif
                enddo
              enddo
            endif
          endif
        enddo do_m
      enddo
    endif
  else

    ! HARMONIC CALCULATION ONLY FOR NOW

    ! Calculate sample points
    echo $mc_data_points $mc_continuation > no_data.dat
    call q_dist ! Program doesn't exist

    ! Main loop for CASTEP calculations
    do i=mc_continuation+1,mc_data_points+mc_continuation
      ! Prepare input file for mc_positions.f90
      awk 'NR==$i { print }' phonon_coord.dat > mc_positions_input.dat
      call mc_positions ! Program doesn't exist.

      ! Write CASTEP input file 'seed.cell' with calculated atomic positions
      write_cell(positions.dat) > $seed_name.cell
      
      ! Run CASTEP
      run_castep
      
      ! Collect CASTEP output data
      dft_output_file = read_dft_output_file(seed_name//'.castep')
      energy = dft_output_file%energy
      read gradient <<< $(cat gradient.dat)   ! JCAP
      echo $scaling $energy $gradient >> neural_energy.dat
      mv gradient.dat gradients/gradient.$i.dat  ! JCAP

      mv $seed_name.castep total_output
      mv $seed_name.bands band_structure
      if (magres) then
        mv $seed_name.magres magres
      endif
    enddo
  endif
end subroutine
end module
