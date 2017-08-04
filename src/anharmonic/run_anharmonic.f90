! ======================================================================
! Runs anharmonic DFT calculations.
! ======================================================================
module run_anharmonic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function run_anharmonic_keywords() result(keywords)
  use help_module
  implicit none
  
  type(KeywordData) :: keywords(4)
  
  keywords = [                                                                &
  & make_keyword( 'qpoints_to_run',                                           &
  &               'qpoints_to_run is the first and last q-point to run. &
  &These should be specified as two integers separated by a space.'),         &
  & make_keyword( 'sampling_points_to_run',                                   &
  &               'sampling_points_to_run is the first and last sampling &
  &point to run. Should only be specified if only one q-point is being run. &
  &These should be specified as two integers separated by a space.',          &
  &               is_optional=.true.),                                        &
  & make_keyword( 'no_cores',                                                 &
  &               'no_cores is the number of cores on which DFT will be run. &
  &This is passed to the specified run script.',                              &
  &               default_value='1'),                                         &
  & make_keyword( 'run_script',                                               &
  &               'run_script is the path to the script for running DFT. An &
  &example run script can be found in doc/input_files.',                      &
  &               is_path=.true.) ]
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine run_anharmonic(arguments)
  use dictionary_module
  use qpoints_module
  use sampling_points_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! Inputs
  integer      :: first_qpoint
  integer      :: last_qpoint
  integer      :: first_sampling_point
  integer      :: last_sampling_point
  integer      :: no_cores
  type(String) :: run_script
  
  ! Previous inputs.
  type(Dictionary) :: setup_harmonic_arguments
  type(Dictionary) :: setup_anharmonic_arguments
  type(String)     :: harmonic_path
  type(String)     :: dft_code
  type(String)     :: seedname
  
  ! Previously calculated data.
  type(QpointData),    allocatable :: qpoints(:)
  type(SamplingPoint), allocatable :: sampling_points(:)
  
  ! Working variables.
  logical :: only_run_part_qpoint
  integer :: qpoint
  integer :: sampling_point
  
  ! Temporary variables.
  type(String), allocatable :: line(:)
  type(String)              :: dir
  integer                   :: result_code
  
  ! Read in inputs.
  wd = arguments%value('working_directory')
  
  line = split(arguments%value('qpoints_to_run'))
  first_qpoint = int(line(1))
  last_qpoint = int(line(2))
  
  only_run_part_qpoint = .false.
  if (arguments%is_set('sampling_points_to_run')) then
    if (first_qpoint/=last_qpoint) then
      call err()
    endif
    only_run_part_qpoint = .true.
    line = split(arguments%value('sampling_points_to_run'))
    first_sampling_point = int(line(1))
    last_sampling_point = int(line(2))
  endif
  
  no_cores = int(arguments%value('no_cores'))
  run_script = arguments%value('run_script')
  
  ! Read in setup_anharmonic settings.
  call setup_anharmonic_arguments%read_file( &
     & wd//'/setup_anharmonic.used_settings')
  harmonic_path = setup_anharmonic_arguments%value('harmonic_path')
  
  ! Read in setup_harmonic settings.
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  dft_code = setup_harmonic_arguments%value('dft_code')
  seedname = setup_harmonic_arguments%value('seedname')
  
  ! Read in previously calculated data.
  qpoints = read_qpoints_file(harmonic_path//'/qpoints_ibz.dat')
  
  if (first_qpoint<1) then
    call err()
  elseif (first_qpoint>last_qpoint) then
    call err()
  elseif (last_qpoint>size(qpoints)) then
    call err()
  endif
  
  ! Run calculations.
  do qpoint=first_qpoint,last_qpoint
    sampling_points = read_sampling_points_file( &
       & wd//'/qpoint_'//qpoint//'/sampling_points.dat')
    
    if (only_run_part_qpoint) then
      if (first_sampling_point<1) then
        call err()
      elseif (first_sampling_point>last_sampling_point) then
        call err()
      elseif (last_sampling_point>size(sampling_points)) then
        call err()
      endif
    else
      first_sampling_point = 1
      last_sampling_point = size(sampling_points)
    endif
    
    do sampling_point=first_sampling_point,last_sampling_point
      call print_line('')
      call print_line('Running sampling point '//sampling_point//' at q-point &
         &'//qpoint//'.')
      dir = wd//'/qpoint_'//qpoint//'sampling_point_'//sampling_point
      result_code = system_call( 'cd '//wd//'; '//run_script//' '// &
         & dft_code//' '//dir//' '//no_cores//' '//seedname)
      call print_line('Result code: '//result_code)
    enddo
  enddo
end subroutine
end module
