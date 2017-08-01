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
  use help_module
  implicit none
  
  type(KeywordData) :: keywords(15)
  
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
  use dictionary_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! ----------------------------------------------------------------------
  ! Read inputs
  ! ----------------------------------------------------------------------
  wd = arguments%value("working_directory")

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
  no_unit_cells = structure%sc_size

  ! Read in number of atoms in unit cell (required)
  basis = structure%no_atoms_prim

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
  !   CREATE OUTPUT DIRECTORIES
  ! ----------------------------------------------------------------------
  call mkdir(wd//'/'//band_structure)
  call mkdir(wd//'/'//total_output)
  call mkdir(wd//'/'//max_amplitude)
  call mkdir(wd//'/'//gradients)

  ! WORKING VARIABLES
  ! Declare variable 'freq_line' that stores the line in the 'disp_patterns.dat' file at which the frequency of the working normal mode is
  freq_line = 1
  ! Declare variable 'disp_line' that stores the line in the 'disp_patterns.dat' file at which the atomic displacement we are working on is
  disp_line = freq_line + 3
  ! Initialise variable 'blank_line' to 0
  blank_line_original = 0

  ! ----------------------------------------------------------------------
  !   MAIN PROGRAM
  ! ----------------------------------------------------------------------
  noofatoms = basis * no_unit_cells
  
  output_file = open_write_file(wd//'/caesar.output')

  ! PAA sampling or MC sampling
  if (.not. mc_sampling) then 
    ! Declare variable 'blank_line_original' that finds the next relevant blank line in 'disp_patterns.dat' file
    blank_line_original=$(awk -v awk_blank_line_original=$blank_line_original 'BEGIN{x=0} {if(!NF && NR>awk_blank_line_original && x==0){x++; print NR}}' disp_patterns.dat)

    ! Calculate ground state energy

    if (file_exists(wd//'/'//'gs_energy.dat')) then
      static_energy=$(awk 'NR==1 {print $2}' gs_energy.dat)
    else
      call print_line(output_file, 'Ground state calculation')

      call write_dft_input_file(supercell, wd//'/'//seed_name//'.cell')

      call run_dft

      ! Collect partial output information
      energy = dft_output_file%energy ! corrected for finite basis set
      static_energy = energy
      call print_line(output_file, ' Ground state energy (eV): '//energy)

      ! Organise output files
      mv $seed_name.bands $seed_name.gs.bands
      mv $seed_name.gs.bands band_structure
      echo 0.0 $static_energy > gs_energy.dat
      mv $seed_name.castep $seed_name.gs.castep
      mv $seed_name.gs.castep total_output

      if (magres) then
          mv $seed_name.magres $seed_name.gs.magres
          mv $seed_name.gs.magres magres
      endif
    endif

    !! Loop over the different normal modes
    !for (( n=$((first_mode)); n<=$((last_mode)); n++ )); do
    !  num_data=$num_indep_data
    !
    !  ! Calculate blank_line, freq_line and disp_line for current normal mode
    !  blank_line=$(( $n*$blank_line_original ))
    !  freq_line=$(( $blank_line-$blank_line_original+1 ))
    !  disp_line=$(($freq_line + 3)) 
    !  gvec_line=$(($freq_line + 1))
    !
    !  ! Calculate only non-symmetry related modes
    !  reference=$(awk -v awk_line=$n 'NR==awk_line {print $2}' symmetry.dat)
    !  if [ $n -eq $reference ]; then
    !  ! Calculate frequency of current normal mode
    !  frequency=$(awk -v awk_freq_line=$freq_line 'NR==awk_freq_line {print $3}' disp_patterns.dat)
    !
    !  ! Write frequency to the output 'energy.dat' file
    !  if [ ! -f energy.dat ]; then
    !    echo $frequency > energy.dat
    !  else
    !    echo $frequency >> energy.dat
    !  fi
    !
    !  i=$first_amp
    !  i_final=$last_amp
    !
    !  echo '' >> caesar.output
    !  echo ' Calculating mode number ' $n >> caesar.output
    !
    !  ! Loop over the atomic displacements for a particular normal mode
    !  while [ $i -le $i_final ]; do
    !
    !    ! Prepare working file 'disp.dat'
    !    echo $frequency $temperature $num_data $i > disp.dat
    !    awk -v awk_disp_line=$disp_line -v awk_blank_line=$blank_line 'NR==awk_disp_line, NR==(awk_blank_line-1) {print $1 " " $2 " " $3 " " $4}' disp_patterns.dat >> disp.dat
    !
    !    ! Compile and run 'positions.f90' to evaluate the atomic positions for the current normal mode amplitude
    !    call positions >> output.txt
    !
    !    scaling=$(awk 'NR==1 {print $1}' positions.dat)
    !    scaling_check=$(awk 'NR==1 {print $2}' positions.dat)
    !
    !    ! Consider only non-translational modes
    !    if [ $scaling_check -eq 0 ]; then
    !      sed '1d' positions.dat > positions_temp.dat ! Delete first line from file positions.dat
    !
    !      ! Write CASTEP input file 'seed.cell' with calculated atomic positions
    !      mv positions_temp.dat positions.dat
    !      write_cell(positions.dat) > $seed_name.cell
    !
    !      ! check_static_energy=$(( ($num_data-1)/2+1 ))
    !      if [ $check_static_energy -eq $i ]; then
    !        echo 0.0 $static_energy >> energy.dat
    !        echo '    - Amplitude ' $i ' energy (eV):' $static_energy >> caesar.output
    !      else
    !        run_castep
    !
    !        cp input_files/$seed_name.cell .
    !
    !        echo $seed_name > force_extract
    !        echo $noofatoms >> force_extract
    !        ! Collect partial output information
    !        line_number=$(awk '/Total energy corrected for finite basis set/{x=NR}END{print x}' $seed_name.castep)
    !        line_number=$(awk '/Final energy/{x=NR}END{print x}' $seed_name.castep)
    !        energy=$(awk -v awk_line_number=$line_number 'NR==awk_line_number {print $9}' $seed_name.castep)
    !        read gradient <<< $(cat gradient.dat)   ! JCAP
    !        echo $scaling $energy $gradient >> energy.dat
    !        echo '    - Amplitude ' $i ' energy (eV):' $energy >> caesar.output
    !        echo '    - Amplitude ' $i ' force (a.u.):' $gradient >> caesar.output ! JCAP
    !        mv gradient.dat gradients/gradient.$n.$i.dat ! JCAP
    !
    !        ! Organise output files
    !        mv $seed_name.bands $seed_name.$n.$i.bands
    !        mv $seed_name.$n.$i.bands band_structure
    !        mv $seed_name.castep $seed_name.$n.$i.castep
    !        mv $seed_name.$n.$i.castep total_output
    !        if [ "$magres" = "true" ]; then
    !          mv $seed_name.magres $seed_name.$n.$i.magres
    !          mv $seed_name.$n.$i.magres magres
    !        fi
    !      fi
    !
    !      echo "Finished CASTEP run $i."
    !
    !      i=$(( $i+1 ))
    !    else
    !      echo ' Skipping translational modes...' >> caesar.output
    !      i=$(( $i_final+1 ))
    !    fi  ! skip translational modes
    !
    !  done ! i 
    !
    !  echo '' >> caesar.output
    !
    !  ! Format 'energy.dat' output file
    !   echo '' >> energy.dat
    !
    !  ! Organise output files
    !  mv max_amplitude.dat max_amplitude.$n.dat
    !  mv max_amplitude.$n.dat max_amplitude
    !
    !  echo
    !  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    !  echo "Finished normal mode $n."
    !  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    !  echo
    !
    !  fi ! symmetry
    !
    !done ! n

    if (coupled_sampling) then
      do n=first_mode,last_mode
        do m=first_mode,last_mode
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
          reference1=$(awk -v awk_line=$n 'NR==awk_line {print $2}' symmetry.dat)
          reference2=$(awk -v awk_line=$m 'NR==awk_line {print $2}' symmetry.dat)
          if (n == reference1) then
            if (m > n .or. (m < n .and. m /= reference2)) then
              ! Calculate frequency of current normal mode
              frequency1=$(awk -v awk_freq_line=$freq_line1 'NR==awk_freq_line {print $3}' disp_patterns.dat)
              frequency2=$(awk -v awk_freq_line=$freq_line2 'NR==awk_freq_line {print $3}' disp_patterns.dat)
              
              ! Write frequencies to the output 'coupled_energy.dat' file
              !echo $frequency1 $frequency2 >> coupled_energy.dat
              
              i = first_amp_2body
              i_final = last_amp_2body
              j = first_amp_2body
              j_final = last_amp_2body
              
              call print_line(output_file, ' Calculating coupling of modes '//n//' & '//m)
              
              ! Loop over the atomic displacements for a particular normal mode
              if (n == 7 .and. m < 10) then
                i_initial = last_amp_2body + 1
              elseif (n == 7 .and. m == 10) then
                i_initial = 10
              endif
              
              do_i : do i=i_initial,i_final
                j = first_amp_2body
                do j=first_amp_2body,j_final
                  ! Prepare working file 'disp.dat'
                  frequency=$frequency1
                  num_data=$num_2body_data
                  echo $frequency $temperature $num_data $i > disp1.dat
                  awk -v awk_disp_line=$disp_line1 -v awk_blank_line=$blank_line1 'NR==awk_disp_line, NR==(awk_blank_line-1) {print $1 " " $2 " " $3 " " $4}' disp_patterns.dat >> disp1.dat
                  frequency=$frequency2
                  echo $frequency $temperature $num_data $j > disp2.dat
                  awk -v awk_disp_line=$disp_line2 -v awk_blank_line=$blank_line2 'NR==awk_disp_line, NR==(awk_blank_line-1) {print $1 " " $2 " " $3 " " $4}' disp_patterns.dat >> disp2.dat
              
                  !Run 'positions2.f90' to evaluate the atomic positions for the current normal mode amplitudes
                  call positions2 >> output.txt
                  
                  scaling1=$(awk 'NR==1 {print $1}' positions.dat)
                  scaling2=$(awk 'NR==1 {print $2}' positions.dat)
                  scaling_check=$(awk 'NR==1 {print $3}' positions.dat)
                  
                  ! Consider only non-translational modes
                  if (scaling_check == 0) then
                    
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
                      
                      echo $seed_name > force_extract2
                      echo $noofatoms >> force_extract2
                      
                      ! Collect partial output information
                      line_number=$(awk '/Final energy/{x=NR}END{print x}' $seed_name.castep)
                      energy=$(awk -v awk_line_number=$line_number 'NR==awk_line_number {print $9}' $seed_name.castep)
                      read -a gradient <<< $(cat gradient.dat)   ! JCAP
                      echo $scaling1 $scaling2 $energy ${gradient[0]} ${gradient[1]} >> coupled_energy.dat
                      echo '    - Amplitude ' $i ' ' $j ' energy (eV):' $energy >> caesar.output
                      echo '    - Amplitude ' $i ' ' $j ' forces (a.u.):' ${gradient[0]} ' ' ${gradient[1]} >> caesar.output ! JCAP
                      
                      mv $seed_name.bands $seed_name.$n.$i.$m.$j.bands
                      mv $seed_name.$n.$i.$m.$j.bands band_structure
                      mv $seed_name.castep $seed_name.$n.$i.$m.$j.castep
                      mv $seed_name.$n.$i.$m.$j.castep total_output
                      if (magres) then
                        mv $seed_name.magres $seed_name.$n.$i.$m.$j.magres
                        mv $seed_name.$n.$i.$m.$j.magres magres
                      endif
                    elseif (check_static_energy == j .and. num_2body_data == num_indep_data) then
                      cp band_structure/$seed_name.$n.$i.bands ./$seed_name.bands
                      cp total_output/$seed_name.$n.$i.castep ./$seed_name.castep
                      if (magres) then
                        cp magres/$seed_name.$n.$i.magres ./$seed_name.magres
                      endif
                      
                      echo $seed_name > force_extract2
                      echo $noofatoms >> force_extract2
                      
                      ! Collect partial output information
                      line_number=$(awk '/Final energy/{x=NR}END{print x}' $seed_name.castep)
                      energy=$(awk -v awk_line_number=$line_number 'NR==awk_line_number {print $9}' $seed_name.castep)
                      read -a gradient <<< $(cat gradient.dat)   ! JCAP
                      echo $scaling1 $scaling2 $energy ${gradient[0]} ${gradient[1]} >> coupled_energy.dat
                      echo '    - Amplitude ' $i ' ' $j ' energy (eV):' $energy >> caesar.output
                      echo '    - Amplitude ' $i ' ' $j ' forces (a.u.):' ${gradient[0]} ' ' ${gradient[1]} >> caesar.output ! JCAP
                      
                      mv $seed_name.bands $seed_name.$n.$i.$m.$j.bands
                      mv $seed_name.$n.$i.$m.$j.bands band_structure
                      mv $seed_name.castep $seed_name.$n.$i.$m.$j.castep
                      mv $seed_name.$n.$i.$m.$j.castep total_output
                      if (magres) then
                        mv $seed_name.magres $seed_name.$n.$i.$m.$j.magres
                        mv $seed_name.$n.$i.$m.$j.magres magres
                      endif
                    else
                      ! Run CASTEP
                      run_castep
                      
                      ! Clean-up after run to remove large output files
                      call rm machine_file
                      cp input_files/$seed_name.* .
                      
                      echo $seed_name > force_extract
                      echo $noofatoms >> force_extract
                      
                      ! Collect partial output information
                      line_number=$(awk '/Final energy/{x=NR}END{print x}' $seed_name.castep)
                      energy=$(awk -v awk_line_number=$line_number 'NR==awk_line_number {print $9}' $seed_name.castep)
                      read -a gradient <<< $(cat gradient.dat)   ! JCAP
                      echo $scaling1 $scaling2 $energy ${gradient[0]} ${gradient[1]} >> coupled_energy.dat
                      echo '    - Amplitude ' $i ' ' $j ' energy (eV):' $energy >> caesar.output
                      echo '    - Amplitude ' $i ' ' $j ' forces (a.u.):' ${gradient[0]} ' ' ${gradient[1]} >> caesar.output ! JCAP
                      mv gradient.dat gradients/gradient.$n.$i.$m.$j.dat ! JCAP
                      
                      ! Organise output files
                      mv $seed_name.bands $seed_name.$n.$i.$m.$j.bands
                      mv $seed_name.$n.$i.$m.$j.bands band_structure
                      mv $seed_name.castep $seed_name.$n.$i.$m.$j.castep
                      mv $seed_name.$n.$i.$m.$j.castep total_output
                      if (magres) then
                        mv $seed_name.magres $seed_name.$n.$i.$m.$j.magres
                        mv $seed_name.$n.$i.$m.$j.magres magres
                      endif
                    endif
                  else
                    echo ' Skipping translational modes...' >> caesar.output
                    exit do_i
                  endif
                enddo
              enddo do_i
            endif
          endif
        enddo
      enddo
    endif
  else ! PAA sampling or MC sampling

    ! HARMONIC CALCULATION ONLY FOR NOW

    ! Calculate sample points
    echo $mc_data_points $mc_continuation > no_data.dat
    call q_dist ! Program doesn't exist

    ! Main loop for CASTEP calculations
    do i=mc_continuation+1,mc_data_points+mc_continuation
      ! Prepare input file for mc_positions.f90
      awk -v awk_line=$i 'NR==awk_line { print }' phonon_coord.dat > mc_positions_input.dat
      call mc_positions

      ! Write CASTEP input file 'seed.cell' with calculated atomic positions
      write_cell(positions.dat) > $seed_name.cell
      
      ! Run CASTEP
      run_castep
      
      echo $seed_name > force_extract
      echo $noofatoms >> force_extract

      ! Collect CASTEP output data
      energy_line_number=$(awk '/Final energy/{x=NR}END{print x}' $seed_name.castep)
      energy=$(awk -v awk_line_number=$energy_line_number 'NR==awk_line_number {print $5}' $seed_name.castep)
      read gradient <<< $(cat gradient.dat)   ! JCAP
      echo $scaling $energy $gradient >> neural_energy.dat
      mv gradient.dat gradients/gradient.$i.dat  ! JCAP

      mv $seed_name.castep $seed_name.$i.castep
      mv $seed_name.$i.castep total_output
      mv $seed_name.bands $seed_name.$i.bands
      mv $seed_name.$i.bands band_structure
      if (magres) then
        mv $seed_name.magres $seed_name.$i.magres
        mv $seed_name.$i.magres magres
      endif

    enddo

    ! Clean-up
    rm phonon_coord.dat

  endif ! PAA sampling or MC sampling
end subroutine
end module
