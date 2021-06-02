submodule (caesar_program_mode_module) caesar_program_mode_submodule
  use caesar_arguments_module
contains

module procedure new_ProgramMode_character_character
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

module procedure new_ProgramMode_character_String
  output = ProgramMode( mode_name,              &
                      & char(description),      &
                      & keywords,               &
                      & main_subroutine,        &
                      & suppress_from_helptext, &
                      & suppress_settings_file)
end procedure

module procedure new_ProgramMode_String_character
  output = ProgramMode( char(mode_name),        &
                      & description,            &
                      & keywords,               &
                      & main_subroutine,        &
                      & suppress_from_helptext, &
                      & suppress_settings_file)
end procedure

module procedure new_ProgramMode_String_String
  output = ProgramMode( char(mode_name),        &
                      & char(description),      &
                      & keywords,               &
                      & main_subroutine,        &
                      & suppress_from_helptext, &
                      & suppress_settings_file)
end procedure

module procedure remove_keyword_String
  integer :: i
  
  i = first(this%keywords%keyword()==keyword, default=0)
  if (i==0) then
    call print_line(CODE_ERROR//': Program mode does not contain keyword '// &
       &keyword)
    call err()
  else
    this%keywords = [this%keywords(:i-1), this%keywords(i+1:)]
  endif
end procedure

module procedure remove_keyword_character
  call this%remove_keyword(str(keyword))
end procedure

module procedure new_Dictionary_ProgramMode
  this = Dictionary(mode%keywords)
end procedure

module procedure print_help_ProgramMode
  call print_line('')
  call print_line(colour(this%mode_name,'cyan'))
  call print_line(this%description)
end procedure
end submodule
