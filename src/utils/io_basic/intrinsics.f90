! ======================================================================
! Provides I/O functionality for intrinsic types.
! ======================================================================
module intrinsics_submodule
  use precision_module
  
  use string_submodule
  use error_submodule
  use print_submodule
  implicit none
  
  private
  
  public :: assignment(=)
  public :: str
  public :: left_pad
  public :: pad_int_to_str
  public :: join
  public :: lgcl
  public :: int
  public :: dble
  public :: cmplx
  public :: operator(//)
  
  interface assignment(=)
    module procedure assign_String_logical
    module procedure assign_String_integer
    module procedure assign_String_real
    module procedure assign_String_complex
  end interface
  
  interface str
    module procedure str_integer
    module procedure str_real
    module procedure str_logical
    module procedure str_complex
  end interface
  
  interface left_pad
    module procedure left_pad_character
    module procedure left_pad_String
  end interface
  
  interface join
    module procedure join_reals
    module procedure join_integers
    module procedure join_logicals
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
    module procedure concatenate_character_integer
    module procedure concatenate_integer_character
    module procedure concatenate_String_integer
    module procedure concatenate_integer_String
    
    module procedure concatenate_character_real
    module procedure concatenate_real_character
    module procedure concatenate_String_real
    module procedure concatenate_real_String
    
    module procedure concatenate_character_logical
    module procedure concatenate_logical_character
    module procedure concatenate_String_logical
    module procedure concatenate_logical_String
    
    module procedure concatenate_String_complex
    module procedure concatenate_complex_String
    module procedure concatenate_character_complex
    module procedure concatenate_complex_character
    
    module procedure concatenate_character_integers
    module procedure concatenate_integers_character
    module procedure concatenate_String_integers
    module procedure concatenate_integers_String
    
    module procedure concatenate_character_reals
    module procedure concatenate_reals_character
    module procedure concatenate_String_reals
    module procedure concatenate_reals_String
    
    module procedure concatenate_character_logicals
    module procedure concatenate_logicals_character
    module procedure concatenate_String_logicals
    module procedure concatenate_logicals_String
    
    module procedure concatenate_character_complexes
    module procedure concatenate_complexes_character
    module procedure concatenate_String_complexes
    module procedure concatenate_complexes_String
  end interface
contains

! ----------------------------------------------------------------------
! Assignment to String from intrinsic types.
! ----------------------------------------------------------------------

subroutine assign_String_logical(output,input)
  implicit none
  
  type(String), intent(out) :: output
  logical,      intent(in)  :: input
  
  if (input) then
    output = 'T'
  else
    output = 'F'
  endif
end subroutine

subroutine assign_String_integer(output,input)
  implicit none
  
  type(String), intent(out) :: output
  integer,      intent(in)  :: input
  
  integer, parameter   :: int_width = 12
  character(int_width) :: int_string
  
  write(int_string,"(I0)") input
  output = trim(int_string)
end subroutine

subroutine assign_String_real(output,input)
  implicit none
  
  type(String), intent(out) :: output
  real(dp),     intent(in)  :: input
  
  integer, parameter    :: real_width = 25
  integer, parameter    :: decimal_places = 17
  type(String)          :: format_string
  character(real_width) :: real_string
  
  format_string = "(ES"//real_width//'.'//decimal_places//")"
  write(real_string, char(format_string)) input
  
  output = real_string
end subroutine

subroutine assign_String_complex(output,input)
  implicit none
  
  type(String), intent(out) :: output
  complex(dp),  intent(in)  :: input
  
  integer, parameter                 :: real_width = 25
  integer, parameter                 :: imag_width = 24
  integer, parameter                 :: decimal_places = 17
  type(String)                       :: format_string
  character(real_width+imag_width+1) :: complex_string
  
  format_string = '('                                        // &
                & 'ES'//real_width//'.'//decimal_places //','// &
                & 'sp'                                  //','// &
                & 'ES'//imag_width//'.'//decimal_places //','// &
                & '"i"'                                      // &
                & ')'
  write(complex_string,char(format_string)) input
  output = complex_string
end subroutine

! ----------------------------------------------------------------------
! The str() function, which converts types to String.
! a=str(b) is equivalent to a=b.
! ----------------------------------------------------------------------

impure elemental function str_logical(this) result(output)
  implicit none
  
  logical, intent(in) :: this
  type(String)        :: output
  
  output = this
end function

impure elemental function str_integer(this) result(output)
  implicit none
  
  integer, intent(in) :: this
  type(String)        :: output
  
  output = this
end function

impure elemental function str_real(this) result(output)
  implicit none
  
  real(dp), intent(in) :: this
  type(String)         :: output
  
  output = this
end function

impure elemental function str_complex(this) result(output)
  implicit none
  
  complex(dp), intent(in) :: this
  type(String)            :: output
  
  output = this
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
  
  character(1) :: pad
  integer      :: i
  
  if (present(pad_character)) then
    pad = pad_character
  else
    pad = '0'
  endif
  
  if (input<0) then
    call err()
  endif
  
  output = str(input)
  
  if (len(output)>len(match)) then
    call err()
  endif
  
  do i=1,len(match)-len(output)
    output = pad//output
  enddo
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
impure elemental function pad_int_to_str(this) result(output)
  implicit none
  
  integer, intent(in) :: this
  type(String)        :: output
  
  output = this
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

function join_reals(this,delimiter) result(output)
  implicit none
  
  real(dp),     intent(in)           :: this(:)
  character(*), intent(in), optional :: delimiter
  type(String)                       :: output
  
  output = join(str(this), delimiter)
end function

function join_integers(this,delimiter,pad_sign) result(output)
  implicit none
  
  integer,      intent(in)           :: this(:)
  character(*), intent(in), optional :: delimiter
  logical,      intent(in), optional :: pad_sign
  type(String)                       :: output
  
  logical :: to_pad
  
  if (present(pad_sign)) then
    to_pad = pad_sign
  else
    to_pad = .true.
  endif
  
  if (to_pad) then
    output = join(pad_int_to_str(this), delimiter)
  else
    output = join(str(this), delimiter)
  endif
end function

function join_logicals(this,delimiter) result(output)
  implicit none
  
  logical,      intent(in)           :: this(:)
  character(*), intent(in), optional :: delimiter
  type(String)                       :: output
  
  output = join(str(this), delimiter)
end function

function join_complexes(this,delimiter) result(output)
  implicit none
  
  complex(dp),  intent(in)           :: this(:)
  character(*), intent(in), optional :: delimiter
  type(String)                       :: output
  
  output = join(str(this), delimiter)
end function

! ----------------------------------------------------------------------
! Concatenation of character(*) or String with intrinsic types.
! Always returns a String.
! ----------------------------------------------------------------------

! --------------------------------------------------
! Concatenation of string types and default types.
! --------------------------------------------------

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
end module
