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
  
  type KeywordData
    type(String) :: keyword
    type(String) :: default_value
    type(String) :: helptext
    logical      :: is_path
  end type
  
  interface make_keyword
    module procedure make_keyword_characters
    module procedure make_keyword_Strings
  end interface
  
  interface help
    module procedure help_default
    module procedure help_keyword
  end interface
  
contains

! ----------------------------------------------------------------------
! Takes a keyword, its default value and its helptext and returns KeywordData.
! ----------------------------------------------------------------------
! If default value is NO_ARGUMENT (from io_module), Caesar will insist on
!    a value being given, and will abort if this does not happen.
! If default value is NOT_SET (from io_module), Caesar will not print the
!    keyword to file.
! If is_path is .true., Caesar will convert the path to an absolute path.
function make_keyword_characters(keyword,default_value,helptext,is_path) &
   & result(this)
  implicit none
  
  character(*), intent(in)           :: keyword
  character(*), intent(in)           :: default_value
  character(*), intent(in)           :: helptext
  logical,      intent(in), optional :: is_path
  type(KeywordData)                  :: this
  
  this%keyword = keyword
  this%default_value = default_value
  this%helptext = helptext
  if (present(is_path)) then
    this%is_path = is_path
  else
    this%is_path = .false.
  endif
end function

function make_keyword_Strings(keyword,default_value,helptext,is_path) &
   & result(this)
  implicit none
  
  type(String), intent(in)           :: keyword
  type(String), intent(in)           :: default_value
  type(String), intent(in)           :: helptext
  logical,      intent(in), optional :: is_path
  type(KeywordData)                  :: this
  
  if (present(is_path)) then
    this = make_keyword(char(keyword),char(default_value),char(helptext), &
       & is_path)
  else
    this = make_keyword(char(keyword),char(default_value),char(helptext))
  endif
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

subroutine help_keyword(keyword,mode,keywords)
  implicit none
  
  type(String),                   intent(in) :: keyword
  type(String),                   intent(in) :: mode
  type(KeywordData), allocatable, intent(in) :: keywords(:)
  
  integer :: i
  logical :: success
  
  if (keyword==NO_ARGUMENT) then
    do i=1,size(keywords)
      call print_line('')
      call print_line(keywords(i)%keyword)
      call print_line(keywords(i)%helptext)
      if (keywords(i)%default_value==NOT_SET) then
        call print_line(keywords(i)%keyword//' is optional.')
      elseif (keywords(i)%default_value==NO_ARGUMENT) then
        call print_line(keywords(i)%keyword//' has no default value, and &
           &must be set.')
      else
        call print_line(keywords(i)%keyword//' has a default value of: '// &
           & keywords(i)%default_value)
      endif
    enddo
  else
    success = .false.
    do i=1,size(keywords)
      if (keywords(i)%keyword==keyword) then
        call print_line(keywords(i)%helptext)
        success = .true.
        exit
      endif
    enddo
    
    if (.not. success) then
      call print_line('')
      call print_line('Keyword '//keyword//' not recognised. For a list of &
         &keywords associated with mode '//mode//', call:')
      call print_line('  caesar '//mode//' -h')
    endif
  endif
end subroutine
end module
