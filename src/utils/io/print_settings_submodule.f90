submodule (caesar_print_settings_module) caesar_print_settings_submodule
  use caesar_io_module
  
  !> Defines whether or not a [[PrintSettings(type)]] should be used.
  logical             :: USE_CURRENT_PRINT_SETTINGS = .false.
  !> The [[PrintSettings(type)]] which should be used.
  type(PrintSettings) :: CURRENT_PRINT_SETTINGS
contains
module procedure new_PrintSettings_character
  if (.not. USE_CURRENT_PRINT_SETTINGS) then
    ! Default arguments.
    this%indent = 0
    this%overhang = 3
    this%decimal_places = 17
    this%floating_point_format = 'es'
    this%integer_digits = 1
  else
    ! Currently set settings.
    this = CURRENT_PRINT_SETTINGS
  endif
  
  ! Settings specific to this call.
  if (present(indent)) then
    this%indent = indent
  endif
  
  if (present(overhang)) then
    this%overhang = overhang
  endif
  
  if (present(decimal_places)) then
    this%decimal_places = decimal_places
  endif
  
  if (present(floating_point_format)) then
    this%floating_point_format = floating_point_format
  endif
  
  if (present(integer_digits)) then
    this%integer_digits = integer_digits
  endif
end procedure

module procedure new_PrintSettings_String
  this = PrintSettings( indent,                      &
                      & overhang,                    &
                      & decimal_places,              &
                      & char(floating_point_format), &
                      & integer_digits               )
end procedure

module procedure set_print_settings_PrintSettings
  CURRENT_PRINT_SETTINGS = settings
  USE_CURRENT_PRINT_SETTINGS = .true.
end procedure

module procedure set_print_settings_arguments_character
  call set_print_settings(PrintSettings( indent,                &
                                       & overhang,              &
                                       & decimal_places,        &
                                       & floating_point_format, &
                                       & integer_digits         ))
end procedure

module procedure set_print_settings_arguments_String
  call set_print_settings(PrintSettings( indent,                &
                                       & overhang,              &
                                       & decimal_places,        &
                                       & floating_point_format, &
                                       & integer_digits         ))
end procedure

module procedure unset_print_settings
  USE_CURRENT_PRINT_SETTINGS = .false.
end procedure
end submodule
