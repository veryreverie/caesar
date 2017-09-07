! ======================================================================
! Handles input keywords and helptext.
! ======================================================================
module help_module
  use constants_module, only : dp
  use string_module
  use io_module
  use keyword_module
  implicit none
  
  private
  
  ! Prints help text.
  public :: help
  
  interface help
    module procedure help_default
    module procedure help_specific
    module procedure help_keyword
  end interface
  
contains

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
