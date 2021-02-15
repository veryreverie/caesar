submodule (caesar_file_module) caesar_file_submodule
contains
module procedure open_file
  type(String) :: formatted_filename
  
  integer :: iostat
  logical :: opened
  
  logical :: unit_found
  
  unit_num = 0
  unit_found = .false.
  
  ! Loop over file units until an unused unit is found.
  do while (.not. unit_found)
    unit_num = unit_num + 1
    
    ! Ignore standard units 0, 5 and 6.
    if (unit_num==5 .or. unit_num==6) then
      cycle
    endif
    
    ! Ignore standard units 100, 101 and 102.
    if (unit_num>=100 .and. unit_num<=102) then
      cycle
    endif
    
    ! Check file unit.
    inquire(unit=unit_num,opened=opened,iostat=iostat)
    if (iostat==0 .and. .not. opened) then
      unit_found = .true.
    endif
  enddo
  
  ! Format filename.
  formatted_filename = format_path(filename)
  
  open( unit   = unit_num,                 &
      & file   = char(formatted_filename), &
      & status = status,                   &
      & action = action,                   &
      & access = access,                   &
      & iostat = iostat                    )
  if (iostat /= 0) then
    call print_line('Error opening '//formatted_filename//' file.')
    call err()
  endif
end procedure

module procedure open_read_file
  unit_num = open_file(filename, 'old', 'read', 'sequential')
end procedure

module procedure open_write_file
  unit_num = open_file(filename, 'unknown', 'write', 'sequential')
end procedure
end submodule
