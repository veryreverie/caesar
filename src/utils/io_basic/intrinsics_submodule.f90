submodule (caesar_intrinsics_module) caesar_intrinsics_submodule
contains
module procedure str_logical
  if (input) then
    output = 'T'
  else
    output = 'F'
  endif
end procedure

module procedure str_integer
  integer, parameter       :: max_int_width = 12
  character(max_int_width) :: int_string
  
  write(int_string,"(I0)") input
  output = trim(int_string)
end procedure

module procedure str_real
  type(PrintSettings)       :: print_settings
  integer                   :: real_width
  type(String)              :: exponent_format
  type(String)              :: format_string
  character(:), allocatable :: real_string
  integer                   :: ialloc
  
  if (present(settings)) then
    print_settings = settings
  else
    print_settings = PrintSettings()
  endif
  
  if (print_settings%floating_point_format=='es') then
    ! - One character for the leading '-' or ' '.
    ! - One character for the decimal point.
    ! - Five characters for the trailing 'E+~~~' or 'E-~~~'.
    real_width = print_settings%decimal_places + 8
    exponent_format = 'e3'
  elseif (print_settings%floating_point_format=='f') then
    ! - One character for the leading '-' or ' '.
    real_width = print_settings%decimal_places &
             & + print_settings%integer_digits &
             & + 2
    exponent_format = ''
  else
    call err()
  endif
  
  format_string = '('                                  // &
                & print_settings%floating_point_format // &
                & str(real_width)                      // &
                & '.'                                  // &
                & str(print_settings%decimal_places)   // &
                & exponent_format                      // &
                & ")"
  
  allocate( character(real_width) :: real_string, &
          & stat=ialloc); call err(ialloc)
  
  write(real_string, char(format_string)) input
  
  output = real_string
end procedure

module procedure str_complex
  type(PrintSettings)       :: print_settings
  integer                   :: real_width
  integer                   :: imag_width
  type(String)              :: exponent_format
  type(String)              :: format_string
  character(:), allocatable :: complex_string
  integer                   :: i,j,ialloc
  
  if (present(settings)) then
    print_settings = settings
  else
    print_settings = PrintSettings()
  endif
  
  if (print_settings%floating_point_format=='es') then
    real_width = print_settings%decimal_places + 8
    imag_width = real_width
    exponent_format = 'e3'
  elseif (print_settings%floating_point_format=='f') then
    real_width = print_settings%decimal_places &
             & + print_settings%integer_digits &
             & + 2
    imag_width = real_width
    exponent_format = ''
  else
    call err()
  endif
  
  format_string =                                                    &
     & '('                                                        // &
     & print_settings%floating_point_format//str(real_width) //'.'// &
     & str(print_settings%decimal_places)//exponent_format   //','// &
     & 'sp'                                                  //','// &
     & print_settings%floating_point_format//str(imag_width) //'.'// &
     & str(print_settings%decimal_places)//exponent_format   //','// &
     & '"i"'                                                      // &
     & ')'
  
  allocate( character(real_width+imag_width+1) :: complex_string, &
          & stat=ialloc); call err(ialloc)
  
  write(complex_string, char(format_string)) input
  
  ! The output should be a single token.
  ! However, if fixed point formatting is used the output can look like
  !    '  12.0 +12.0i'. The imaginary part should be padded with zeroes
  !    rather than spaces.
  if (print_settings%floating_point_format=='f') then
    do i=len(complex_string),real_width+1,-1
      if (complex_string(i:i)==' ') then
        complex_string(real_width+1:real_width+1) = complex_string(i+1:i+1)
        do j=real_width+2,i+1
          complex_string(j:j) = '0'
        enddo
      endif
    enddo
  endif
  
  output = complex_string
end procedure

module procedure str_logicals_character
  integer :: i,ialloc
  
  allocate(output(0), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    output = [output, str(input(i))]
    if (i<size(input) .and. present(separating_line)) then
      output = [output, str(separating_line)]
    endif
  enddo
end procedure

module procedure str_logicals_String
  output = str(input,char(separating_line),settings)
end procedure

module procedure str_integers_character
  integer :: i,ialloc
  
  allocate(output(0), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    output = [output, str(input(i))]
    if (i<size(input) .and. present(separating_line)) then
      output = [output, str(separating_line)]
    endif
  enddo
end procedure

module procedure str_integers_String
  output = str(input,char(separating_line),settings)
end procedure

module procedure str_reals_character
  integer :: i,ialloc
  
  allocate(output(0), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    output = [output, str(input(i))]
    if (i<size(input) .and. present(separating_line)) then
      output = [output, str(separating_line)]
    endif
  enddo
end procedure

module procedure str_reals_String
  output = str(input,char(separating_line),settings)
end procedure

module procedure str_complexes_character
  integer :: i,ialloc
  
  allocate(output(0), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    output = [output, str(input(i))]
    if (i<size(input) .and. present(separating_line)) then
      output = [output, str(separating_line)]
    endif
  enddo
end procedure

module procedure str_complexes_String
  output = str(input,char(separating_line),settings)
end procedure

module procedure left_pad_character
  if (input<0) then
    call err()
  endif
  
  output = str(input)
  
  if (len(output)>len(match)) then
    call err()
  endif
  
  if (present(pad_character)) then
    output = repeat(pad_character, len(match)-len(output)) // output
  else
    output = repeat('0', len(match)-len(output)) // output
  endif
end procedure

module procedure left_pad_String
  output = left_pad(input, char(match), pad_character)
end procedure

module procedure pad_int_to_str
  output = str(this,settings)
  if (slice(output,1,1)/='-') then
    output = ' '//output
  endif
end procedure

module procedure lgcl_character
  integer :: ierr
  
  read(this,*,iostat=ierr) output
  if (ierr/=0) then
    call print_line(ERROR//': unable to convert the string "'//this//'" to &
       &a logical.')
    call err()
  endif
end procedure

module procedure lgcl_String
  output = lgcl(char(this))
end procedure

module procedure int_character
  integer :: ierr
  
  read(this,*,iostat=ierr) output
  if (ierr/=0) then
    call print_line(ERROR//': unable to convert the string "'//this//'" to &
       &an integer.')
    call err()
  endif
end procedure

module procedure int_String
  output = int(char(this))
end procedure

module procedure dble_character
  type(String), allocatable :: split_string(:)
  
  integer :: ierr
  
  split_string = split_line(this,'/')
  
  if (size(split_string)==1) then
    ! The string does not contain a '/': treat it as a real.
    read(this,*,iostat=ierr) output
    if (ierr/=0) then
      call print_line(ERROR//': unable to convert the string "'//this//'" to &
         &a real number.')
      call err()
    endif
  elseif (size(split_string)==2) then
    ! The string contains a '/': treat it as a fraction.
    output = dble(char(split_string(1)))/dble(char(split_string(2)))
  else
    ! The string contains multiple '/'s: this is an error.
    call print_line(ERROR//': unable to convert the string "'//this//'" to &
       &a real number.')
    call err()
  endif
end procedure

module procedure dble_String
  output = dble(char(this))
end procedure

module procedure cmplx_character
  logical :: characters_before_sign
  integer :: i,j
  
  ! Find the '+' or '-' which separates the real and imaginary parts.
  characters_before_sign = .false.
  do i=1,len(this)
    ! Ignore an opening '+' or '-'.
    if (.not. characters_before_sign) then
      ! N.B. this comparison is [space] or [tab].
      if (this(i:i)/=' ' .and. this(i:i)/='	') then
        characters_before_sign = .true.
      endif
      cycle
    endif
    
    if (this(i:i)=='+' .or. this(i:i)=='-') then
      ! Ignore the '+' or '-' if it is preceeded by an 'E' or 'e'.
      if (i>1) then
        if (lower_case(this(i-1:i-1))=='e') then
          cycle
        endif
      endif
      
      ! Find the 'i' after the imaginary part.
      j = index(this(i+1:), 'i')
      if (j==0) then
        call print_line(ERROR//': Complex number contains real and imaginary &
           &parts but not "i".')
        call print_line(this)
        call err()
      endif
      
      output = cmplx( dble(this(1:i-1)), dble(this(i:i+j-1)), dp)
      return
    endif
  enddo
  
  ! There is no '+' or '-', so the number must be purely real
  !    or purely imaginary.
  j = index(this, 'i')
  if (j==0) then
    output = dble(this)
  else
    output = cmplx(0.0_dp, dble(this(1:j-1)), dp)
  endif
end procedure

module procedure cmplx_String
  output = cmplx(char(this))
end procedure

module procedure join_logicals
  output = join( str(this,settings), &
               & delimiter           )
end procedure

module procedure join_integers
  logical :: to_pad
  
  if (present(pad_sign)) then
    to_pad = pad_sign
  else
    to_pad = .true.
  endif
  
  if (to_pad) then
    output = join( pad_int_to_str(this, settings), &
                 & delimiter                       )
  else
    output = join( str(this, settings), &
                 & delimiter            )
  endif
end procedure

module procedure join_reals
  output = join( str(this, settings), &
               & delimiter            )
end procedure

module procedure join_complexes
  output = join( str(this,settings), &
               & delimiter           )
end procedure

module procedure concatenate_String_logical
  output = this//str(that)
end procedure

module procedure concatenate_logical_String
  output = str(this)//that
end procedure

module procedure concatenate_String_integer
  output = this//str(that)
end procedure

module procedure concatenate_integer_String
  output = str(this)//that
end procedure

module procedure concatenate_String_real
  output = this//str(that)
end procedure

module procedure concatenate_real_String
  output = str(this)//that
end procedure

module procedure concatenate_String_complex
  output = this//str(that)
end procedure

module procedure concatenate_complex_String
  output = str(this)//that
end procedure

module procedure concatenate_String_logicals
  output = this//join(that)
end procedure

module procedure concatenate_logicals_String
  output = join(this)//that
end procedure

module procedure concatenate_String_integers
  output = this//join(that)
end procedure

module procedure concatenate_integers_String
  output = join(this)//that
end procedure

module procedure concatenate_String_reals
  output = this//join(that)
end procedure

module procedure concatenate_reals_String
  output = join(this)//that
end procedure

module procedure concatenate_String_complexes
  output = this//join(that)
end procedure

module procedure concatenate_complexes_String
  output = join(this)//that
end procedure

module procedure concatenate_character_logical
  output = this//str(that)
end procedure

module procedure concatenate_logical_character
  output = str(this)//that
end procedure

module procedure concatenate_character_integer
  output = this//str(that)
end procedure

module procedure concatenate_integer_character
  output = str(this)//that
end procedure

module procedure concatenate_character_real
  output = this//str(that)
end procedure

module procedure concatenate_real_character
  output = str(this)//that
end procedure

module procedure concatenate_character_complex
  output = this//str(that)
end procedure

module procedure concatenate_complex_character
  output = str(this)//that
end procedure

module procedure concatenate_character_logicals
  output = this//join(that)
end procedure

module procedure concatenate_logicals_character
  output = join(this)//that
end procedure

module procedure concatenate_character_integers
  output = this//join(that)
end procedure

module procedure concatenate_integers_character
  output = join(this)//that
end procedure

module procedure concatenate_character_reals
  output = this//join(that)
end procedure

module procedure concatenate_reals_character
  output = join(this)//that
end procedure

module procedure concatenate_character_complexes
  output = this//join(that)
end procedure

module procedure concatenate_complexes_character
  output = join(this)//that
end procedure

module procedure print_line_logical
  call print_line(str(this,settings))
end procedure

module procedure print_line_integer
  call print_line(str(this,settings))
end procedure

module procedure print_line_real
  call print_line(str(this,settings))
end procedure

module procedure print_line_complex
  call print_line(str(this,settings))
end procedure

module procedure print_line_logicals
  call print_line(join(this, settings=settings))
end procedure

module procedure print_line_integers
  call print_line(join(this, settings=settings))
end procedure

module procedure print_line_reals
  call print_line(join(this, settings=settings))
end procedure

module procedure print_line_complexes
  call print_line(join(this, settings=settings))
end procedure

module procedure print_lines_logicals_character
  call print_lines(str(this,separating_line,settings))
end procedure

module procedure print_lines_logicals_String
  call print_lines(str(this,separating_line,settings))
end procedure

module procedure print_lines_integers_character
  call print_lines(str(this,separating_line,settings))
end procedure

module procedure print_lines_integers_String
  call print_lines(str(this,separating_line,settings))
end procedure

module procedure print_lines_reals_character
  call print_lines(str(this,separating_line,settings))
end procedure

module procedure print_lines_reals_String
  call print_lines(str(this,separating_line,settings))
end procedure

module procedure print_lines_complexes_character
  call print_lines(str(this,separating_line,settings))
end procedure

module procedure print_lines_complexes_String
  call print_lines(str(this,separating_line,settings))
end procedure
end submodule
