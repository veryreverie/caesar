! ======================================================================
! Provides I/O functionality for intrinsic types.
! ======================================================================
module intrinsics_module
  use precision_module
  
  use string_module
  use error_module
  use print_settings_module
  use print_module
  implicit none
  
  private
  
  public :: str
  public :: left_pad
  public :: pad_int_to_str
  public :: join
  public :: lgcl
  public :: int
  public :: dble
  public :: cmplx
  public :: operator(//)
  public :: print_line
  public :: print_lines
  
  interface str
    module procedure str_logical
    module procedure str_integer
    module procedure str_real
    module procedure str_complex
    
    module procedure str_logicals_character
    module procedure str_logicals_String
    module procedure str_integers_character
    module procedure str_integers_String
    module procedure str_reals_character
    module procedure str_reals_String
    module procedure str_complexes_character
    module procedure str_complexes_String
  end interface
  
  interface left_pad
    module procedure left_pad_character
    module procedure left_pad_String
  end interface
  
  interface join
    module procedure join_logicals
    module procedure join_integers
    module procedure join_reals
    module procedure join_complexes
  end interface
  
  interface lgcl
    module procedure lgcl_character
    module procedure lgcl_String
  end interface
  
  interface int
    module procedure int_character
    module procedure int_String
  end interface
  
  interface dble
    module procedure dble_character
    module procedure dble_String
  end interface
  
  interface cmplx
    module procedure cmplx_character
    module procedure cmplx_String
  end interface

  interface operator(//)
    module procedure concatenate_character_logical
    module procedure concatenate_logical_character
    module procedure concatenate_String_logical
    module procedure concatenate_logical_String
    
    module procedure concatenate_character_integer
    module procedure concatenate_integer_character
    module procedure concatenate_String_integer
    module procedure concatenate_integer_String
    
    module procedure concatenate_character_real
    module procedure concatenate_real_character
    module procedure concatenate_String_real
    module procedure concatenate_real_String
    
    module procedure concatenate_String_complex
    module procedure concatenate_complex_String
    module procedure concatenate_character_complex
    module procedure concatenate_complex_character
    
    module procedure concatenate_character_logicals
    module procedure concatenate_logicals_character
    module procedure concatenate_String_logicals
    module procedure concatenate_logicals_String
    
    module procedure concatenate_character_integers
    module procedure concatenate_integers_character
    module procedure concatenate_String_integers
    module procedure concatenate_integers_String
    
    module procedure concatenate_character_reals
    module procedure concatenate_reals_character
    module procedure concatenate_String_reals
    module procedure concatenate_reals_String
    
    module procedure concatenate_character_complexes
    module procedure concatenate_complexes_character
    module procedure concatenate_String_complexes
    module procedure concatenate_complexes_String
  end interface
  
  interface print_line
    module procedure print_line_logical
    module procedure print_line_integer
    module procedure print_line_real
    module procedure print_line_complex
    
    module procedure print_line_logicals
    module procedure print_line_integers
    module procedure print_line_reals
    module procedure print_line_complexes
  end interface
  
  interface print_lines
    module procedure print_lines_logicals_character
    module procedure print_lines_logicals_String
    module procedure print_lines_integers_character
    module procedure print_lines_integers_String
    module procedure print_lines_reals_character
    module procedure print_lines_reals_String
    module procedure print_lines_complexes_character
    module procedure print_lines_complexes_String
  end interface
contains

! ----------------------------------------------------------------------
! The str() function, which converts types to String.
! ----------------------------------------------------------------------

impure elemental function str_logical(input,settings) result(output)
  implicit none
  
  logical, intent(in)                       :: input
  type(PrintSettings), intent(in), optional :: settings
  type(String)                              :: output
  
  if (input) then
    output = 'T'
  else
    output = 'F'
  endif
end function

impure elemental function str_integer(input,settings) result(output)
  implicit none
  
  integer, intent(in)                       :: input
  type(PrintSettings), intent(in), optional :: settings
  type(String)                              :: output
  
  integer, parameter       :: max_int_width = 12
  character(max_int_width) :: int_string
  
  write(int_string,"(I0)") input
  output = trim(int_string)
end function

impure elemental function str_real(input,settings) result(output)
  implicit none
  
  real(dp), intent(in)                      :: input
  type(PrintSettings), intent(in), optional :: settings
  type(String)                              :: output
  
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
end function

impure elemental function str_complex(input,settings) result(output)
  implicit none
  
  complex(dp), intent(in)                   :: input
  type(PrintSettings), intent(in), optional :: settings
  type(String)                              :: output
  
  type(PrintSettings)       :: print_settings
  integer                   :: real_width
  integer                   :: imag_width
  type(String)              :: exponent_format
  type(String)              :: format_string
  character(:), allocatable :: complex_string
  integer                   :: ialloc
  
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
  
  output = complex_string
end function

function str_logicals_character(input,separating_line,settings) result(output)
  implicit none
  
  logical,             intent(in)           :: input(:)
  character(*),        intent(in), optional :: separating_line
  type(PrintSettings), intent(in), optional :: settings
  type(String), allocatable                 :: output(:)
  
  integer :: i
  
  output = [String::]
  do i=1,size(input)
    output = [output, str(input(i))]
    if (i<size(input) .and. present(separating_line)) then
      output = [output, str(separating_line)]
    endif
  enddo
end function

function str_logicals_String(input,separating_line,settings) result(output)
  implicit none
  
  logical,             intent(in)           :: input(:)
  type(String),        intent(in)           :: separating_line
  type(PrintSettings), intent(in), optional :: settings
  type(String), allocatable                 :: output(:)
  
  output = str(input,char(separating_line),settings)
end function

function str_integers_character(input,separating_line,settings) result(output)
  implicit none
  
  integer,             intent(in)           :: input(:)
  character(*),        intent(in), optional :: separating_line
  type(PrintSettings), intent(in), optional :: settings
  type(String), allocatable                 :: output(:)
  
  integer :: i
  
  output = [String::]
  do i=1,size(input)
    output = [output, str(input(i))]
    if (i<size(input) .and. present(separating_line)) then
      output = [output, str(separating_line)]
    endif
  enddo
end function

function str_integers_String(input,separating_line,settings) result(output)
  implicit none
  
  integer,             intent(in)           :: input(:)
  type(String),        intent(in)           :: separating_line
  type(PrintSettings), intent(in), optional :: settings
  type(String), allocatable                 :: output(:)
  
  output = str(input,char(separating_line),settings)
end function

function str_reals_character(input,separating_line,settings) result(output)
  implicit none
  
  real(dp),            intent(in)           :: input(:)
  character(*),        intent(in), optional :: separating_line
  type(PrintSettings), intent(in), optional :: settings
  type(String), allocatable                 :: output(:)
  
  integer :: i
  
  output = [String::]
  do i=1,size(input)
    output = [output, str(input(i))]
    if (i<size(input) .and. present(separating_line)) then
      output = [output, str(separating_line)]
    endif
  enddo
end function

function str_reals_String(input,separating_line,settings) result(output)
  implicit none
  
  real(dp),            intent(in)           :: input(:)
  type(String),        intent(in)           :: separating_line
  type(PrintSettings), intent(in), optional :: settings
  type(String), allocatable                 :: output(:)
  
  output = str(input,char(separating_line),settings)
end function

function str_complexes_character(input,separating_line,settings) result(output)
  implicit none
  
  complex(dp),         intent(in)           :: input(:)
  character(*),        intent(in), optional :: separating_line
  type(PrintSettings), intent(in), optional :: settings
  type(String), allocatable                 :: output(:)
  
  integer :: i
  
  output = [String::]
  do i=1,size(input)
    output = [output, str(input(i))]
    if (i<size(input) .and. present(separating_line)) then
      output = [output, str(separating_line)]
    endif
  enddo
end function

function str_complexes_String(input,separating_line,settings) result(output)
  implicit none
  
  complex(dp),         intent(in)           :: input(:)
  type(String),        intent(in)           :: separating_line
  type(PrintSettings), intent(in), optional :: settings
  type(String), allocatable                 :: output(:)
  
  output = str(input,char(separating_line),settings)
end function

! ----------------------------------------------------------------------
! Converts an integer to String, then left pads that string until it is
!    as long as the given match string.
! Uses pad_character to pad, which defaults to '0'.
! ----------------------------------------------------------------------

impure elemental function left_pad_character(input,match,pad_character) &
   & result(output)
  implicit none
  
  integer,      intent(in)           :: input
  character(*), intent(in)           :: match
  character(1), intent(in), optional :: pad_character
  type(String)                       :: output
  
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
end function

impure elemental function left_pad_String(input,match,pad_character) &
   & result(output)
  implicit none
  
  integer,      intent(in)           :: input
  type(String), intent(in)           :: match
  character(1), intent(in), optional :: pad_character
  type(String)                       :: output
  
  output = left_pad(input, char(match), pad_character)
end function

! --------------------------------------------------
! Converts an integer to a String, and then pads the string so that positive
!    and negative integers end up being the same length.
! --------------------------------------------------
impure elemental function pad_int_to_str(this,settings) result(output)
  implicit none
  
  integer,             intent(in)           :: this
  type(PrintSettings), intent(in), optional :: settings
  type(String)                              :: output
  
  output = str(this,settings)
  if (slice(output,1,1)/='-') then
    output = ' '//output
  endif
end function

! ----------------------------------------------------------------------
! Conversion from String or character(*) to intrinsic type.
! ----------------------------------------------------------------------

function lgcl_character(this) result(output)
  implicit none
  
  character(*), intent(in) :: this
  logical                  :: output
  
  integer :: ierr
  
  read(this,*,iostat=ierr) output
  if (ierr/=0) then
    call print_line(ERROR//': unable to convert the string "'//this//'" to &
       &a logical.')
    call err()
  endif
end function

impure elemental function lgcl_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  logical                  :: output
  
  output = lgcl(char(this))
end function

function int_character(this) result(output)
  implicit none
  
  character(*), intent(in) :: this
  integer                  :: output
  
  integer :: ierr
  
  read(this,*,iostat=ierr) output
  if (ierr/=0) then
    call print_line(ERROR//': unable to convert the string "'//this//'" to &
       &an integer.')
    call err()
  endif
end function

impure elemental function int_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  integer                  :: output
  
  output = int(char(this))
end function

! Can also convert fractions as strings to real(dp), e.g. '4/5' -> 0.8_dp.
recursive function dble_character(this) result(output)
  implicit none
  
  character(*), intent(in) :: this
  real(dp)                 :: output
  
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
end function

impure elemental function dble_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  real(dp)                 :: output
  
  output = dble(char(this))
end function

function cmplx_character(this) result(output)
  implicit none
  
  character(*), intent(in) :: this
  complex(dp)              :: output
  
  integer :: i
  logical :: split_allowed
  
  split_allowed = .false.
  do i=len(this)-1,1,-1
    if (this(i:i)=='E') then
      split_allowed = .true.
    elseif (split_allowed .and. this(i:i)=='+') then
      output = cmplx(  dble(this(1:i-1)),           &
                    &  dble(this(i+1:len(this)-1)), &
                    &  dp)
      exit
    elseif (split_allowed .and. this(i:i)=='-') then
      output = cmplx(  dble(this(1:i-1)),           &
                    & -dble(this(i+1:len(this)-1)), &
                    &  dp)
      exit
    endif
  enddo
end function

impure elemental function cmplx_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  complex(dp)              :: output
  
  output = cmplx(char(this))
end function

! ----------------------------------------------------------------------
! Convert an intrinsic array to an array of Strings,
!    then concatenate them into a single string.
! ----------------------------------------------------------------------

function join_logicals(this,delimiter,settings) result(output)
  implicit none
  
  logical,             intent(in)           :: this(:)
  character(*),        intent(in), optional :: delimiter
  type(PrintSettings), intent(in), optional :: settings
  type(String)                              :: output
  
  output = join( str(this,settings), &
               & delimiter           )
end function

function join_integers(this,delimiter,pad_sign,settings) result(output)
  implicit none
  
  integer,             intent(in)           :: this(:)
  character(*),        intent(in), optional :: delimiter
  logical,             intent(in), optional :: pad_sign
  type(PrintSettings), intent(in), optional :: settings
  type(String)                              :: output
  
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
end function

function join_reals(this,delimiter,settings) result(output)
  implicit none
  
  real(dp),            intent(in)           :: this(:)
  character(*),        intent(in), optional :: delimiter
  type(PrintSettings), intent(in), optional :: settings
  type(String)                              :: output
  
  output = join( str(this, settings), &
               & delimiter            )
end function

function join_complexes(this,delimiter,settings) result(output)
  implicit none
  
  complex(dp),         intent(in)           :: this(:)
  character(*),        intent(in), optional :: delimiter
  type(PrintSettings), intent(in), optional :: settings
  type(String)                              :: output
  
  output = join( str(this,settings), &
               & delimiter           )
end function

! ----------------------------------------------------------------------
! Concatenation of character(*) or String with intrinsic types.
! Always returns a String.
! ----------------------------------------------------------------------

! --------------------------------------------------
! Concatenation of string types and default types.
! --------------------------------------------------

! String = String//logical
function concatenate_String_logical(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  logical,       intent(in) :: that
  type(String)              :: output
  
  output = this//str(that)
end function

! String = logical//String
function concatenate_logical_String(this,that) result(output)
  implicit none
  
  logical,       intent(in) :: this
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = str(this)//that
end function

! String = String//integer
function concatenate_String_integer(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  integer,       intent(in) :: that
  type(String)              :: output
  
  output = this//str(that)
end function

! String = integer//String
function concatenate_integer_String(this,that) result(output)
  implicit none
  
  integer,       intent(in) :: this
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = str(this)//that
end function

! String = String//real(dp)
function concatenate_String_real(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  real(dp),      intent(in) :: that
  type(String)              :: output
  
  output = this//str(that)
end function

! String = real(dp)//String
function concatenate_real_String(this,that) result(output)
  implicit none
  
  real(dp),      intent(in) :: this
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = str(this)//that
end function

! String = String//complex
function concatenate_String_complex(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  complex(dp),   intent(in) :: that
  type(String)              :: output
  
  output = this//str(that)
end function

! String = complex//String
function concatenate_complex_String(this,that) result(output)
  implicit none
  
  complex(dp),   intent(in) :: this
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = str(this)//that
end function

! String = String//logical(:)
function concatenate_String_logicals(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  logical,       intent(in) :: that(:)
  type(String)              :: output
  
  output = this//join(that)
end function

! String = logical(:)//String
function concatenate_logicals_String(this,that) result(output)
  implicit none
  
  logical,       intent(in) :: this(:)
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = join(this)//that
end function

! String = String//integer(:)
function concatenate_String_integers(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  integer,       intent(in) :: that(:)
  type(String)              :: output
  
  output = this//join(that)
end function

! String = integer(:)//String
function concatenate_integers_String(this,that) result(output)
  implicit none
  
  integer,       intent(in) :: this(:)
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = join(this)//that
end function

! String = String//real(dp)(:)
function concatenate_String_reals(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  real(dp),      intent(in) :: that(:)
  type(String)              :: output
  
  output = this//join(that)
end function

! String = real(dp)(:)//String
function concatenate_reals_String(this,that) result(output)
  implicit none
  
  real(dp),      intent(in) :: this(:)
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = join(this)//that
end function

! String = String//complex(:)
function concatenate_String_complexes(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  complex(dp),   intent(in) :: that(:)
  type(String)              :: output
  
  output = this//join(that)
end function

! String = complex(:)//String
function concatenate_complexes_String(this,that) result(output)
  implicit none
  
  complex(dp),   intent(in) :: this(:)
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = join(this)//that
end function

! String = character//logical
function concatenate_character_logical(this,that) result(output)
  implicit none
  
  character(*), intent(in) :: this
  logical,      intent(in) :: that
  type(String)             :: output
  
  output = this//str(that)
end function

! String = logical//character
function concatenate_logical_character(this,that) result(output)
  implicit none
  
  logical,      intent(in) :: this
  character(*), intent(in) :: that
  type(String)             :: output
  
  output = str(this)//that
end function

! String = character//integer
function concatenate_character_integer(this,that) result(output)
  implicit none
  
  character(*), intent(in) :: this
  integer,      intent(in) :: that
  type(String)             :: output
  
  output = this//str(that)
end function

! String = integer//character
function concatenate_integer_character(this,that) result(output)
  implicit none
  
  integer,      intent(in) :: this
  character(*), intent(in) :: that
  type(String)             :: output
  
  output = str(this)//that
end function

! String = character//real(dp)
function concatenate_character_real(this,that) result(output)
  implicit none
  
  character(*), intent(in) :: this
  real(dp),     intent(in) :: that
  type(String)             :: output
  
  output = this//str(that)
end function

! String = real(dp)//character
function concatenate_real_character(this,that) result(output)
  implicit none
  
  real(dp),     intent(in) :: this
  character(*), intent(in) :: that
  type(String)             :: output
  
  output = str(this)//that
end function

! String = character//complex
function concatenate_character_complex(this,that) result(output)
  implicit none
  
  character(*), intent(in) :: this
  complex(dp),  intent(in) :: that
  type(String)             :: output
  
  output = this//str(that)
end function

! String = complex//character
function concatenate_complex_character(this,that) result(output)
  implicit none
  
  complex(dp),  intent(in) :: this
  character(*), intent(in) :: that
  type(String)             :: output
  
  output = str(this)//that
end function

! String = character//logical(:)
function concatenate_character_logicals(this,that) result(output)
  implicit none
  
  character(*), intent(in) :: this
  logical,      intent(in) :: that(:)
  type(String)             :: output
  
  output = this//join(that)
end function

! String = logical(:)//character
function concatenate_logicals_character(this,that) result(output)
  implicit none
  
  logical,      intent(in) :: this(:)
  character(*), intent(in) :: that
  type(String)             :: output
  
  output = join(this)//that
end function

! String = character//integer(:)
function concatenate_character_integers(this,that) result(output)
  implicit none
  
  character(*), intent(in) :: this
  integer,      intent(in) :: that(:)
  type(String)             :: output
  
  output = this//join(that)
end function

! String = integer(:)//character
function concatenate_integers_character(this,that) result(output)
  implicit none
  
  integer,      intent(in) :: this(:)
  character(*), intent(in) :: that
  type(String)             :: output
  
  output = join(this)//that
end function

! String = character//real(dp)(:)
function concatenate_character_reals(this,that) result(output)
  implicit none
  
  character(*), intent(in) :: this
  real(dp),     intent(in) :: that(:)
  type(String)             :: output
  
  output = this//join(that)
end function

! String = real(dp)(:)//character
function concatenate_reals_character(this,that) result(output)
  implicit none
  
  real(dp),     intent(in) :: this(:)
  character(*), intent(in) :: that
  type(String)             :: output
  
  output = join(this)//that
end function

! String = character//complex(:)
function concatenate_character_complexes(this,that) result(output)
  implicit none
  
  character(*), intent(in) :: this
  complex(dp),  intent(in) :: that(:)
  type(String)             :: output
  
  output = this//join(that)
end function

! String = complex(:)//character
function concatenate_complexes_character(this,that) result(output)
  implicit none
  
  complex(dp),  intent(in) :: this(:)
  character(*), intent(in) :: that
  type(String)             :: output
  
  output = join(this)//that
end function

! ----------------------------------------------------------------------
! Provides the print_line and print_lines functions for non-character types.
! ----------------------------------------------------------------------
subroutine print_line_logical(this,settings)
  implicit none
  
  logical,             intent(in)           :: this
  type(PrintSettings), intent(in), optional :: settings
  
  call print_line(str(this,settings))
end subroutine

subroutine print_line_integer(this,settings)
  implicit none
  
  integer,             intent(in)           :: this
  type(PrintSettings), intent(in), optional :: settings
  
  call print_line(str(this,settings))
end subroutine

subroutine print_line_real(this,settings)
  implicit none
  
  real(dp),            intent(in)           :: this
  type(PrintSettings), intent(in), optional :: settings
  
  call print_line(str(this,settings))
end subroutine

subroutine print_line_complex(this,settings)
  implicit none
  
  complex(dp),         intent(in)           :: this
  type(PrintSettings), intent(in), optional :: settings
  
  call print_line(str(this,settings))
end subroutine

subroutine print_line_logicals(this,settings)
  implicit none
  
  logical,             intent(in)           :: this(:)
  type(PrintSettings), intent(in), optional :: settings
  
  call print_line(join(this, settings=settings))
end subroutine

subroutine print_line_integers(this,settings)
  implicit none
  
  integer,             intent(in)           :: this(:)
  type(PrintSettings), intent(in), optional :: settings
  
  call print_line(join(this, settings=settings))
end subroutine

subroutine print_line_reals(this,settings)
  implicit none
  
  real(dp),            intent(in)           :: this(:)
  type(PrintSettings), intent(in), optional :: settings
  
  call print_line(join(this, settings=settings))
end subroutine

subroutine print_line_complexes(this,settings)
  implicit none
  
  complex(dp),         intent(in)           :: this(:)
  type(PrintSettings), intent(in), optional :: settings
  
  call print_line(join(this, settings=settings))
end subroutine

subroutine print_lines_logicals_character(this,separating_line,settings)
  implicit none
  
  logical,             intent(in)           :: this(:)
  character(*),        intent(in), optional :: separating_line
  type(PrintSettings), intent(in), optional :: settings
  
  call print_lines(str(this,separating_line,settings))
end subroutine

subroutine print_lines_logicals_String(this,separating_line,settings)
  implicit none
  
  logical,             intent(in)           :: this(:)
  type(String),        intent(in)           :: separating_line
  type(PrintSettings), intent(in), optional :: settings
  
  call print_lines(str(this,separating_line,settings))
end subroutine

subroutine print_lines_integers_character(this,separating_line,settings)
  implicit none
  
  integer,             intent(in)           :: this(:)
  character(*),        intent(in), optional :: separating_line
  type(PrintSettings), intent(in), optional :: settings
  
  call print_lines(str(this,separating_line,settings))
end subroutine

subroutine print_lines_integers_String(this,separating_line,settings)
  implicit none
  
  integer,             intent(in)           :: this(:)
  type(String),        intent(in)           :: separating_line
  type(PrintSettings), intent(in), optional :: settings
  
  call print_lines(str(this,separating_line,settings))
end subroutine

subroutine print_lines_reals_character(this,separating_line,settings)
  implicit none
  
  real(dp),            intent(in)           :: this(:)
  character(*),        intent(in), optional :: separating_line
  type(PrintSettings), intent(in), optional :: settings
  
  call print_lines(str(this,separating_line,settings))
end subroutine

subroutine print_lines_reals_String(this,separating_line,settings)
  implicit none
  
  real(dp),            intent(in)           :: this(:)
  type(String),        intent(in)           :: separating_line
  type(PrintSettings), intent(in), optional :: settings
  
  call print_lines(str(this,separating_line,settings))
end subroutine

subroutine print_lines_complexes_character(this,separating_line,settings)
  implicit none
  
  complex(dp),         intent(in)           :: this(:)
  character(*),        intent(in), optional :: separating_line
  type(PrintSettings), intent(in), optional :: settings
  
  call print_lines(str(this,separating_line,settings))
end subroutine

subroutine print_lines_complexes_String(this,separating_line,settings)
  implicit none
  
  complex(dp),         intent(in)           :: this(:)
  type(String),        intent(in)           :: separating_line
  type(PrintSettings), intent(in), optional :: settings
  
  call print_lines(str(this,separating_line,settings))
end subroutine
end module
