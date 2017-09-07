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
  use keyword_module
  implicit none
  
  type(KeywordData) :: keywords(8)
  
  keywords = [                                                                &
  & make_keyword( 'first_qpoint',                                             &
  &               'first_qpoint is the id of the first q-point at which &
  &calculations will be run. Details of q-points can be found in &
  &qpoint.dat.',                                                              &
  &               default_value='1'),                                         &
  & make_keyword( 'first_coupling',                                           &
  &               'first_coupling is the id of the first coupling at the &
  &first q-point at which calculations will be run. Details of the couplings &
  &at q-point i can be found in qpoint_i/coupling.dat.',                      &
  &               default_value='1'),                                         &
  & make_keyword( 'first_sampling_point',                                     &
  &               'first_sampling_point is the id of the first sampling point &
  &at the first coupling at the q-point at which calculations will be run. &
  &Details of the couplings at sampling point j at q-point i can be found in &
  &qpoint_i/coupling_j/sampling_points.dat.',                                 &
  &               default_value='1'),                                         &
  & make_keyword( 'last_qpoint',                                              &
  &               'last_qpoint is the id of the last q-point at which &
  &calculations will be run.',                                                &
  &               is_optional=.true.),                                        &
  & make_keyword( 'last_coupling',                                            &
  &               'last_coupling is the id of the last coupling at the &
  &last q-point at which calculations will be run.',                          &
  &               is_optional=.true.),                                        &
  & make_keyword( 'last_sampling_point',                                      &
  &               'last_sampling_point is the id of the last sampling point &
  &at the last coupling at the q-point at which calculations will be run.',   &
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
  use setup_harmonic_module
  use setup_anharmonic_module
  use dictionary_module
  use qpoints_module
  use coupling_module
  use sampling_points_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! Inputs
  integer      :: first_qpoint
  integer      :: first_coupling
  integer      :: first_sampling_point
  integer      :: last_qpoint
  integer      :: last_coupling
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
  type(CoupledModes),  allocatable :: couplings(:)
  type(SamplingPoint), allocatable :: sampling_points(:)
  
  ! Working variables.
  integer :: qpoint
  integer :: coupling,min_coupling,max_coupling
  integer :: sampling_point,min_sampling_point,max_sampling_point
  
  ! Temporary variables.
  type(String) :: qdir,cdir,sdir
  integer      :: result_code
  
  ! Read in inputs.
  wd = arguments%value('working_directory')
  first_qpoint = int(arguments%value('first_qpoint'))
  first_coupling = int(arguments%value('first_coupling'))
  first_sampling_point = int(arguments%value('first_sampling_point'))
  if (arguments%is_set('last_qpoint')) then
    last_qpoint = int(arguments%value('last_qpoint'))
  else
  if (arguments%is_set('last_coupling')) then
    last_coupling = int(arguments%value('last_coupling'))
  endif
  if (arguments%is_set('last_sampling_point')) then
    last_sampling_point = int(arguments%value('last_sampling_point'))
  endif
  no_cores = int(arguments%value('no_cores'))
  run_script = arguments%value('run_script')
  
  call print_line(run_script)
  call err()
  
  ! Check that calculation range is well defined.
  if (first_qpoint<1) then
    call err()
  elseif (first_coupling<1) then
    call err()
  elseif (first_sampling_point<1) then
    call err()
  endif
  
  if (arguments%is_set('last_qpoint')) then
    elseif (first_qpoint>last_qpoint) then
      call err()
    elseif (first_qpoint==last_qpoint) then
      if (arguments%is_set('last_coupling')) then
        if (first_coupling>last_coupling) then
          call err()
        elseif (first_coupling==last_coupling) then
          if (arguments%is_set('last_sampling_point')) then
            if (first_sampling_point>last_sampling_point) then
              call err()
            endif
          endif
        endif
      endif
    endif
  endif
  
  ! Read in setup_anharmonic settings.
  setup_anharmonic_arguments = setup_anharmonic_keywords()
  call setup_anharmonic_arguments%read_file( &
     & wd//'/setup_anharmonic.used_settings')
  harmonic_path = setup_anharmonic_arguments%value('harmonic_path')
  
  ! Read in setup_harmonic settings.
  setup_harmonic_arguments = setup_harmonic_keywords()
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  dft_code = setup_harmonic_arguments%value('dft_code')
  seedname = setup_harmonic_arguments%value('seedname')
  
  ! Run calculations at each q-point.
  qpoints = read_qpoints_file(harmonic_path//'/qpoints_ibz.dat')
  
  if (arguments%is_set('last_qpoint')) then
    if (last_qpoint>size(qpoints)) then
      call err()
    endif
  else
    last_qpoint = size(qpoints)
  endif
  
  do qpoint=first_qpoint,last_qpoint
    qdir = wd//'/qpoint_'//qpoint
    
    ! Run calculations at each coupling at the q-point.
    couplings = read_coupling_file(qdir//'/coupling.dat')
    
    if (qpoint==first_qpoint) then
      min_coupling = first_coupling
    else
      min_coupling = 1
    endif
    
    if (qpoint==last_qpoint) then
      if (arguments%is_set('last_coupling')) then
        if (last_coupling>size(couplings)) then
          call err()
        endif
      else
        last_coupling = size(couplings)
      endif
      max_coupling = last_coupling
    else
      max_coupling = size(couplings)
    endif
    
    do coupling=min_coupling,max_coupling
      cdir = qdir//'/coupling_'//coupling
      
      ! Run calculations at each sampling point at the coupling.
      sampling_points = read_sampling_points_file(cdir//'/sampling_points.dat')
      
      if (qpoint==first_qpoint .and. coupling==first_coupling) then
        min_sampling_point = first_sampling_point
      else
        min_sampling_point = size(sampling_points)
      endif
      
      if (qpoint==last_qpoint .and. coupling==last_coupling) then
        if (arguments%is_set('last_sampling_point')) then
          if (last_sampling_point>size(sampling_points) ) then
            call err()
          endif
        else
          last_sampling_point = size(sampling_points)
        endif
        max_sampling_point = last_sampling_point
      else
        max_sampling_point = size(sampling_points)
      endif
      
      do sampling_point=min_sampling_point,max_sampling_point
        call print_line('')
        call print_line('Running calculations &
           &at sampling point '//sampling_point//' &
           &at coupling '//coupling//' &
           &at q-point '//qpoint//'.')
        sdir = cdir//'/sampling_point_'//sampling_point
        result_code = system_call( 'cd '//wd//'; ' //      &
                                 & run_script      //' '// &
                                 & dft_code        //' '// &
                                 & sdir            //' '// &
                                 & no_cores        //' '// &
                                 & seedname)
        call print_line('Result code: '//result_code)
      enddo
    enddo
  enddo
end subroutine
end module
