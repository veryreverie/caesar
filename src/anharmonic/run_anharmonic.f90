! ======================================================================
! Runs anharmonic DFT calculations, as set up by setup_anharmonic.
! ======================================================================
module run_anharmonic_module
  use common_module
  
  use setup_harmonic_module
  
  use shared_module
  use polynomial_module
  
  use setup_anharmonic_module
  implicit none
  
  private
  
  public :: run_anharmonic_keywords
  public :: run_anharmonic_mode
  public :: run_anharmonic
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function run_anharmonic_keywords() result(keywords)
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                                &
  & KeywordData( 'couplings_to_run',                                          &
  &              'couplings_to_run specifies the couplings at which &
  &electronic structure calculations should be run. Couplings should be &
  &specified as integers linked by hyphens and separated by spaces, e.g. to &
  &run calculations at coupling 5 and coupling 3-7, couplings_to_run should &
  &be "5 3-7"',                                                               &
  &              is_optional=.true.),                                         &
  & KeywordData( 'run_script',                                                &
  &              'run_script is the path to the script for running DFT. An &
  &example run script can be found in doc/input_files.',                      &
  &               is_path=.true.),                                            &
  & KeywordData( 'no_cores',                                                  &
  &              'no_cores is the number of cores on which DFT will be run. &
  &This is passed to the specified run script.',                              &
  &               default_value='1') ]
end function

function run_anharmonic_mode() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'run_anharmonic'
  output%description = 'Runs DFT calculations set up by setup_anharmonic. &
     &Should be run after setup_anharmonic.'
  output%keywords = run_anharmonic_keywords()
  output%main_subroutine => run_anharmonic
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine run_anharmonic(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! Inputs
  type(String), allocatable :: couplings_to_run(:)
  type(String)              :: run_script
  integer                   :: no_cores
  
  ! Previous inputs.
  type(Dictionary) :: setup_harmonic_arguments
  type(Dictionary) :: setup_anharmonic_arguments
  type(String)     :: harmonic_path
  type(String)     :: file_type
  type(String)     :: seedname
  
  ! Coupling variables.
  integer                             :: max_subspace_id
  type(String)                        :: max_subspace_id_string
  type(CoupledSubspaces), allocatable :: coupled_subspaces(:)
  integer,                allocatable :: coupling_ints(:)
  type(String),           allocatable :: coupling_strings(:)
  type(String)                        :: coupling_string
  
  ! Sampling point variables.
  type(SamplingPoints) :: sampling_points
  
  ! VSCF R-vector variables.
  type(VscfRvectors), allocatable :: vscf_rvectors(:)
  
  ! Files and directories.
  type(IFile)                    :: coupled_subspaces_file
  type(String)                   :: coupled_subspaces_dir
  type(IFile)                    :: sampling_points_file
  type(String)                   :: sampling_points_dir
  type(IFile)                    :: vscf_rvectors_file
  type(String)                   :: vscf_rvectors_dir
  type(StringArray), allocatable :: file_sections(:)
  type(String),      allocatable :: directories(:)
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  integer :: result_code
  
  ! Read in inputs.
  wd = arguments%value('working_directory')
  if (arguments%is_set('couplings_to_run')) then
    couplings_to_run = split_line(arguments%value('couplings_to_run'))
  endif
  run_script = arguments%value('run_script')
  no_cores = int(arguments%value('no_cores'))
  
  ! Read in setup_anharmonic settings.
  setup_anharmonic_arguments = Dictionary(setup_anharmonic_keywords())
  call setup_anharmonic_arguments%read_file( &
     & wd//'/setup_anharmonic.used_settings')
  harmonic_path = setup_anharmonic_arguments%value('harmonic_path')
  
  ! Read in setup_harmonic settings.
  setup_harmonic_arguments = Dictionary(setup_harmonic_keywords())
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  file_type = setup_harmonic_arguments%value('file_type')
  seedname = setup_harmonic_arguments%value('seedname')
  
  ! Read in coupled subspaces.
  coupled_subspaces_file = IFile(wd//'/coupling.dat')
  allocate( coupled_subspaces(size(coupled_subspaces_file)), &
          & stat=ialloc); call err(ialloc)
  max_subspace_id = 0
  do i=1,size(coupled_subspaces_file%lines())
    coupled_subspaces(i) = coupled_subspaces_file%line(i)
    max_subspace_id = max(max_subspace_id,maxval(coupled_subspaces(i)%ids,1))
  enddo
  max_subspace_id_string = str(max_subspace_id)
  
  ! Convert input argument 'couplings' into convenient format.
  if (arguments%is_set('couplings_to_run')) then
    allocate( coupling_strings(size(couplings_to_run)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(couplings_to_run)
      coupling_ints = int(split_line(couplings_to_run(i),delimiter='-'))
      coupling_ints = coupling_ints(sort(coupling_ints))
      if (coupling_ints(size(coupling_ints)) > max_subspace_id) then
        call print_line(ERROR//': Input coupling ID larger than maximum.')
        call err()
      endif
      coupling_strings(i) = join( left_pad( coupling_ints,         &
                                &           str(max_subspace_id)), &
                                & delimiter='-')
    enddo
  endif
  
  ! --------------------------------------------------
  ! List directories for running calculations in.
  ! --------------------------------------------------
  directories = [String::]
  
  ! Loop over couplings.
  do i=1,size(coupled_subspaces_file%lines())
    coupling_string = join( left_pad( coupled_subspaces(i)%ids, &
                          &           max_subspace_id_string),  &
                          & delimiter='-')
    
    ! Ignore this coupling if it is not in the input list.
    if (arguments%is_set('couplings_to_run')) then
      if (.not. any(coupling_strings==coupling_string)) then
        cycle
      endif
    endif
    
    coupled_subspaces_dir = wd//'/coupling_'//coupling_string
    
    ! Read in sampling points for this coupling.
    sampling_points_file = IFile(coupled_subspaces_dir//'/sampling_points.dat')
    sampling_points = sampling_points_file%lines()
    
    ! Loop over sampling points.
    do j=1,size(sampling_points)
      sampling_points_dir = coupled_subspaces_dir//'/sampling_point_'// &
                          & left_pad(j,str(size(sampling_points)))
      
      ! Read in VSCF R-vectors.
      vscf_rvectors_file = IFile(sampling_points_dir//'/vscf_rvectors.dat')
      file_sections = split_into_sections(vscf_rvectors_file%lines())
      allocate( vscf_rvectors(size(file_sections)), &
              & stat=ialloc); call err(ialloc)
      do k=1,size(file_sections)
        vscf_rvectors(k) = file_sections(k)
      enddo
      
      ! Loop over VSCF R-vectors.
      do k=1,size(vscf_rvectors)
        vscf_rvectors_dir = sampling_points_dir//'/vscf_rvector_'// &
           & left_pad(k,str(size(vscf_rvectors)))
        
        ! Add directory to directories list.
        directories = [directories, vscf_rvectors_dir]
      enddo
      deallocate(vscf_rvectors,stat=ialloc); call err(ialloc)
    enddo
  enddo
  
  ! --------------------------------------------------
  ! Loop over directories, running calculations in each.
  ! --------------------------------------------------
  do i=1,size(directories)
    ! Run calculation at each R-vector.
    call print_line('')
    call print_line( 'Running calculation in directory '//i// &
                   & ' of '//size(directories)//':')
    call print_line(vscf_rvectors_dir)
    result_code = system_call(  &
       & 'cd '//wd//';' //' '// &
       & run_script     //' '// &
       & file_type      //' '// &
       & directories(i) //' '// &
       & no_cores       //' '// &
       & seedname)
    call print_line('Result code: '//result_code)
  enddo
end subroutine
end module
