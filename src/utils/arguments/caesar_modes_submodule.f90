submodule (caesar_caesar_modes_module) caesar_caesar_modes_submodule
  use caesar_arguments_module
  
  ! An array of modes, which will be populated at startup.
  type(CaesarMode), allocatable :: CAESAR_MODES(:)
contains

module procedure new_CaesarMode_character
  type(CaesarModes) :: modes
  
  modes = CaesarModes()
  
  this = modes%mode(mode)
end procedure

module procedure new_CaesarMode_String
  this = CaesarMode(char(mode))
end procedure

module procedure new_CaesarMode_character_character
  output%mode_name = lower_case(mode_name)
  output%description = description
  output%keywords = keywords
  output%main_subroutine => main_subroutine
  
  if (present(suppress_from_helptext)) then
    output%suppress_from_helptext = suppress_from_helptext
  else
    output%suppress_from_helptext = .false.
  endif
  
  if (present(suppress_settings_file)) then
    output%suppress_settings_file = suppress_settings_file
  else
    output%suppress_settings_file = .false.
  endif
end procedure

module procedure new_CaesarMode_character_String
  output = CaesarMode( mode_name,              &
                     & char(description),      &
                     & keywords,               &
                     & main_subroutine,        &
                     & suppress_from_helptext, &
                     & suppress_settings_file)
end procedure

module procedure new_CaesarMode_String_character
  output = CaesarMode( char(mode_name),        &
                     & description,            &
                     & keywords,               &
                     & main_subroutine,        &
                     & suppress_from_helptext, &
                     & suppress_settings_file)
end procedure

module procedure new_CaesarMode_String_String
  output = CaesarMode( char(mode_name),        &
                     & char(description),      &
                     & keywords,               &
                     & main_subroutine,        &
                     & suppress_from_helptext, &
                     & suppress_settings_file)
end procedure

module procedure remove_keyword_String
  integer :: i
  
  i = first(this%keywords%keyword()==keyword, default=0)
  if (i/=0) then
    this%keywords = [this%keywords(:i-1), this%keywords(i+1:)]
  endif
end procedure

module procedure remove_keyword_character
  call this%remove_keyword(str(keyword))
end procedure

module procedure new_Dictionary_CaesarMode
  this = Dictionary(mode%keywords)
end procedure

module procedure print_help_CaesarMode
  call print_line('')
  call print_line(colour(this%mode_name,'cyan'))
  call print_line(this%description)
end procedure

module procedure new_CaesarModes
  if (.not. allocated(CAESAR_MODES)) then
    call print_line(CODE_ERROR//': CAESAR_MODES has not been allocated.')
    call err()
  else
    this%modes_ = CAESAR_MODES
  endif
end procedure

module procedure new_CaesarModes_CaesarMode
  this%modes_ = modes
end procedure

module procedure mode_String
  type(String) :: lower_case_mode_name
  
  integer :: i
  
  lower_case_mode_name = lower_case(mode_name)
  
  i = first(this%modes_%mode_name==lower_case_mode_name, default=0)
  
  if (i==0) then
    call print_line(ERROR//': Unrecognised mode: '//mode_name)
    call print_line('Call '//colour('caesar -h','white')//' for help.')
    call quit()
  endif
  
  output = this%modes_(i)
end procedure

module procedure mode_character
  output = this%mode(str(mode_name))
end procedure

module procedure print_help_CaesarModes
  integer :: i
  
  do i=1,size(this%modes_)
    call this%modes_(i)%print_help()
  enddo
end procedure

module procedure add_mode_CaesarMode
  if (.not. allocated(CAESAR_MODES)) then
    CAESAR_MODES = [mode]
  else
    if (any(CAESAR_MODES%mode_name==mode%mode_name)) then
      call print_line(CODE_ERROR//': Trying to add the same mode twice.')
      call err()
    else
      CAESAR_MODES = [CAESAR_MODES, mode]
    endif
  endif
end procedure
end submodule
