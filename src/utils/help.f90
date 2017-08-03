! ======================================================================
! Handles input keywords and helptext.
! ======================================================================
module help_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  private
  
  public :: KeywordData  ! A keywords and its helptext.
  public :: make_keyword ! Bundles a keyword and helptext into a KeywordData.
  public :: help         ! Prints help text.
  
  ! --------------------------------------------------
  ! A key:value object with descriptors. 
  ! --------------------------------------------------
  type KeywordData
    ! Keyword and flag (alternative one-character keyword).
    type(String) :: keyword
    character(1) :: flag
    
    ! Helptext for help calls.
    type(String) :: helptext
    
    ! Defaults.
    type(String) :: default_value
    type(String) :: default_keyword
    
    ! Properties of the value.
    logical      :: is_boolean
    logical      :: is_optional
    logical      :: is_path
    
    ! The value itself.
    logical      :: is_set
    type(String) :: value
  end type
  
  interface help
    module procedure help_default
    module procedure help_specific
    module procedure help_keyword
  end interface
  
contains

! ----------------------------------------------------------------------
! Takes a keyword, helptext and options and returns a KeywordData.
! ----------------------------------------------------------------------
! default_value is the value which will be defaulted to if the keyword is not
!    specified. Defaults to not set.
! default_keyword is the keyword whose value will be copied if the keyword is
!    not specified. Defaults to not set.
! At most one of default_value and default_keyword should be specified.
! If is_boolean is .true., the keyword will not take an argument.
!    Boolean keywords should not have defaults specified.
!    is_boolean defaults to .false..
! is_optional specifies whether or not the keyword is required to be given.
!    is_optional defaults to .false. unless a default is given or is_boolean
!    is .true., where it defaults to and is required to be .true..
! is_path specifies that the value will be a path to a file or directory.
!    Paths are converted to absolute paths (from /) automatically.
!    Should not be specified if is_boolean is .true..
!    Defaults to .false..
! flag specifies the flag by which the keyword can be alternately called.
!    Defaults to ' '.
function make_keyword(keyword,helptext,default_value,default_keyword, &
   & is_boolean,is_optional,is_path,flag) result(this)
  implicit none
  
  character(*), intent(in)           :: keyword
  character(*), intent(in)           :: helptext
  character(*), intent(in), optional :: default_value
  character(*), intent(in), optional :: default_keyword
  logical,      intent(in), optional :: is_boolean
  logical,      intent(in), optional :: is_optional
  logical,      intent(in), optional :: is_path
  character(1), intent(in), optional :: flag
  type(KeywordData)                  :: this
  
  ! Set keyword; flag, if present; and helptext.
  this%keyword = lower_case(keyword)
  
  if (present(flag)) then
    this%flag = flag
  else
    this%flag = ' '
  endif
  this%helptext = helptext
  
  ! Set type of keyword.
  if (present(is_boolean)) then
    this%is_boolean = is_boolean
  else
    this%is_boolean = .false.
  endif
  
  if (present(is_optional)) then
    this%is_optional = is_optional
  else
    if (this%is_boolean) then
      this%is_optional = .true.
    elseif (present(default_value) .or. present(default_keyword)) then
      this%is_optional = .true.
    else
      this%is_optional = .false.
    endif
  endif
  
  if (present(is_path)) then
    this%is_path = is_path
  else
    this%is_path = .false.
  endif
  
  if (this%is_boolean) then
    ! Check for conflicts with the keyword being boolean.
    if (this%is_path) then
      call print_line('Error: the keyword '//keyword//' cannot be boolean and &
         &be a path.')
      call err()
    elseif (.not. this%is_optional) then
      call print_line('Error: the keyword '//keyword//' cannot be boolean and &
         &not be optional.')
      call err()
    elseif (present(default_value)) then
      call print_line('Error: the keyword '//keyword//' cannot be boolean and &
         &have a default value.')
      call err()
    endif
    
    ! Set default value, if present.
    this%default_value = ''
    if (present(default_keyword)) then
      this%default_keyword = lower_case(default_keyword)
    else
      this%default_keyword = ''
    endif
  else
    ! Check for conflicts between defaults.
    if (present(default_value) .and. present(default_keyword)) then
      call print_line('Error: the keyword '//keyword//' cannot have both a &
         &default value and a default keyword set.')
      call err()
    elseif ( (present(default_value) .or. present(default_keyword)) .and. &
           & .not. this%is_optional) then
      call print_line('Error: the keyword '//keyword//' cannot have a default &
         &and not be optional.')
        call err()
    endif
    
    ! Set defaults, if present.
    if (present(default_value)) then
      this%default_value = default_value
      this%default_keyword = ''
    elseif (present(default_keyword)) then
      this%default_value = ''
      this%default_keyword = lower_case(default_keyword)
    else
      this%default_value = ''
      this%default_keyword = ''
    endif
  endif
  
  ! Set value and is_set to unset state.
  this%is_set = .false.
  this%value = ''
end function

! ----------------------------------------------------------------------
! Prints helptext.
! ----------------------------------------------------------------------
  
! Print default helptext.
subroutine help_default()
  call print_line('caesar mode [-h [keyword]] [-i] [-f input_file] &
     &[-d working_directory] [--options]')
  call print_line('')
  call print_line('Flags')
  call print_line('  -h [keyword] | --help [keyword]')
  call print_line('      "caesar -h" displays this help text.')
  call print_line('      "caesar mode -h" displays help text relevant to the &
     &specified mode.')
  call print_line('      "caesar mode -h keyword" displays help text relevant &
     &to the specified keyword, in the context of the specified mode.')
  call print_line('')
  call print_line('  -i | --interactive')
  call print_line('      Runs interactively, prompting the user to review and &
     &set all options.')
  call print_line('')
  call print_line('  -f filename | --input_file filename')
  call print_line('      Reads additional settings from specified file.')
  call print_line('      These should be of the form:')
  call print_line('         keyword1 argument')
  call print_line('         keyword2 argument1 argument2 argument3')
  call print_line('      Keywords are the same as command-line --keywords.')
  call print_line('      The "--" prefix. should not be given.')
  call print_line('      The keywords "filename", "interactive" and "help" &
     &should not be specified in a file.')
  call print_line('')
  call print_line('  -d dirname | --working_directory dirname')
  call print_line('      Specifies the directory where Caesar should work.')
  call print_line('      All files and folders will be created here.')
  call print_line('      This is also where any run scripts will be called.')
  call print_line('')
  call print_line('Harmonic modes')
  call print_line('  setup_harmonic')
  call print_line('      Sets up harmonic calculation.')
  call print_line('      DFT code choices are: castep.')
  call print_line('  run_harmonic')
  call print_line('      Runs harmonic calculation.')
  call print_line('      Should be called after setup_harmonic.')
  call print_line('  lte_harmonic')
  call print_line('      Runs harmonic calculations.')
  call print_line('      Should be called after run_harmonic.')
  call print_line('')
  call print_line('Quadratic modes')
  call print_line('  setup_quadratic')
  call print_line('      Sets up quadratic calculation.')
  call print_line('      DFT code choices are: castep.')
  call print_line('      Should be called after lte_harmonic.')
  call print_line('  run_quadratic')
  call print_line('      Runs quadratic calculation.')
  call print_line('      Should be called after setup_quadratic.')
  call print_line('  anharmonic')
  call print_line('      Runs anharmonic calculations.')
  call print_line('      Should be called after run_quadratic.')
  call print_line('  bs_quadratic')
  call print_line('      Runs band structure calculations.')
  call print_line('      Should be called after run_quadratic.')
  call print_line('')
  call print_line('Utility modes')
  call print_line('  hartree_to_eV')
  call print_line('      Provides a Hartree to eV calculator.')
  call print_line('  get_kpoints')
  call print_line('      [Help text pending]')
  call print_line('  calculate_gap')
  call print_line('      [Help text pending]')
end subroutine

! Prints the helptext for a particular mode or keyword.
subroutine help_specific(keyword,mode,keywords)
  implicit none
  
  type(String),      intent(in) :: keyword
  type(String),      intent(in) :: mode
  type(KeywordData), intent(in) :: keywords(:)
  
  integer :: i
  
  if (keyword=='') then
    do i=1,size(keywords)
      call help(keywords(i))
    enddo
    stop
  else
    do i=1,size(keywords)
      if (keywords(i)%keyword==keyword) then
        call help(keywords(i))
        stop
      endif
    enddo
    
    call print_line('')
    call print_line('Keyword '//keyword//' not recognised. For a list of &
       &keywords associated with mode '//mode//', call:')
    call print_line('  caesar '//mode//' -h')
  endif
end subroutine

! Prints help corresponding to a specific keyword.
subroutine help_keyword(keyword)
  implicit none
  
  type(KeywordData), intent(in) :: keyword
  
  call print_line('')
  call print_line(keyword%keyword)
  call print_line(keyword%helptext)
  if (keyword%is_boolean) then
    call print_line(keyword%keyword//' is either set or unset, and &
       &takes no argument.')
  elseif (keyword%default_keyword/='') then
    call print_line(keyword%keyword//' defaults to the same value as &
       &keyword "'//keyword%default_keyword//'".')
  elseif (keyword%default_value/='') then
    call print_line(keyword%keyword//' has a default value of "'// &
       & keyword%default_value//'".')
  elseif (keyword%is_optional) then
    call print_line(keyword%keyword//' is optional.')
  else
    call print_line(keyword%keyword//' has no default value, and &
       &must be set.')
  endif
end subroutine
end module
