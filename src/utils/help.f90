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
    type(String) :: helptext
  end type
  
  interface make_keyword
    module procedure make_keyword_character_character
    module procedure make_keyword_character_String
    module procedure make_keyword_String_character
    module procedure make_keyword_String_String
  end interface
  
contains

! ----------------------------------------------------------------------
! Takes a keyword and its helptext and returns a KeywordData.
! ----------------------------------------------------------------------
function make_keyword_character_character(keyword,helptext) result(this)
  implicit none
  
  character(*), intent(in) :: keyword
  character(*), intent(in) :: helptext
  type(KeywordData)        :: this
  
  this%keyword  = keyword
  this%helptext = helptext
end function

function make_keyword_character_String(keyword,helptext) result(this)
  implicit none
  
  character(*), intent(in) :: keyword
  type(String), intent(in) :: helptext
  type(KeywordData)        :: this
  
  this = make_keyword(keyword, char(helptext))
end function

function make_keyword_String_character(keyword,helptext) result(this)
  implicit none
  
  type(String), intent(in) :: keyword
  character(*), intent(in) :: helptext
  type(KeywordData)        :: this
  
  this = make_keyword(char(keyword), helptext)
end function

function make_keyword_String_String(keyword,helptext) result(this)
  implicit none
  
  type(String), intent(in) :: keyword
  type(String), intent(in) :: helptext
  type(KeywordData)        :: this
  
  this = make_keyword(char(keyword), char(helptext))
end function

! ----------------------------------------------------------------------
! Prints helptext.
! ----------------------------------------------------------------------
subroutine help(keyword,mode,keywords)
  implicit none
  
  type(String),                   intent(in)           :: keyword
  type(String),                   intent(in), optional :: mode
  type(KeywordData), allocatable, intent(in), optional :: keywords(:)
  
  integer :: i
  logical :: success
  
  ! Print default helptext.
  if (keyword=='') then
    call print_line('caesar mode [-h] [-i] [-f input_file] &
       &[-d working_directory] [--options]')
    call print_line('')
    call print_line('Flags')
    call print_line('  -h, --help')
    call print_line('      caesar -h displays this help text.')
    call print_line('      caesar mode -h displays mode-relevant help text.')
    call print_line('      caesar --help arg displays help relevant to arg.')
    call print_line('')
    call print_line('  -i, --interactive')
    call print_line('      Runs interactively.')
    call print_line('')
    call print_line('  -f filename | --input_file filename')
    call print_line('      Reads additional settings from specified file.')
    call print_line('      These should be of the form:')
    call print_line('         keyword1 argument')
    call print_line('         keyword2 argument1 argument2 argument3')
    call print_line('      Keywords are the same as command-line --keywords.')
    call print_line('      The "--" prefix. should not be given.')
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
  elseif (.not. present(keywords)) then
    call print_line('')
    call print_line('For keyword-specific help, please also specify the mode, &
       &e.g.')
    call print_line('')
    call print_line('  caesar setup_harmonic --help dft_code')
  elseif (keyword=='setup_harmonic') then
    do i=1,size(keywords)
      call print_line('')
      call print_line(keywords(i)%keyword)
      call print_line(keywords(i)%helptext)
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
