!> Provides the `token` and `tokens` functions,
!>    which split strings by delimiters.
module caesar_token_module
  use caesar_error_module
  use caesar_string_module
  use caesar_print_module
  use caesar_intrinsics_module
  implicit none
  
  private
  
  public :: token
  public :: tokens
  
  interface token
    !> Split a string into tokens by the given `delimiter` or `delimiters`,
    !>    and return the token at the given index.
    module function token_String(input,index,delimiter,delimiters) &
       & result(output) 
      type(String), intent(in)           :: input
      integer,      intent(in)           :: index
      character(1), intent(in), optional :: delimiter
      character(1), intent(in), optional :: delimiters(:)
      type(String)                       :: output
    end function

    !> Split a string into tokens by the given `delimiter` or `delimiters`,
    !>    and return the token at the given index.
    module function token_character(input,index,delimiter,delimiters) &
       & result(output) 
      character(*), intent(in)           :: input
      integer,      intent(in)           :: index
      character(1), intent(in), optional :: delimiter
      character(1), intent(in), optional :: delimiters(:)
      type(String)                       :: output
    end function
  end interface
  
  interface tokens
    !> Split a string into a `tokens` array by the given
    !>    `delimiter` or `delimiters`,
    !> Returns `tokens(first:last)`.
    !> `first` defaults to 1, `last` defaults to `len(tokens)`.
    module function tokens_String(input,first,last,delimiter,delimiters) &
       & result(output) 
      type(String), intent(in)           :: input
      integer,      intent(in), optional :: first
      integer,      intent(in), optional :: last
      character(1), intent(in), optional :: delimiter
      character(1), intent(in), optional :: delimiters(:)
      type(String), allocatable          :: output(:)
    end function

    !> Split a string into a `tokens` array by the given
    !>    `delimiter` or `delimiters`,
    !> Returns `tokens(first:last)`.
    !> `first` defaults to 1, `last` defaults to `len(tokens)`.
    module function tokens_character(input,first,last,delimiter,delimiters) &
       & result(output) 
      character(*), intent(in)           :: input
      integer,      intent(in), optional :: first
      integer,      intent(in), optional :: last
      character(1), intent(in), optional :: delimiter
      character(1), intent(in), optional :: delimiters(:)
      type(String), allocatable          :: output(:)
    end function
  end interface
end module
