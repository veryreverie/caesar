!> Provides the [[OFile(type)]] class, which handles output file writing.
module caesar_ofile_module
  use caesar_precision_module
  use caesar_io_basic_module
  use caesar_abstract_module
  
  use caesar_ofile_target_module
  use caesar_string_writeable_module
  use caesar_strings_writeable_module
  implicit none
  
  private
  
  public :: OFile
  public :: assignment(=)
  
  !> Handles output file writing.
  !> Constructing this class constructs an [[OFileTarget(type)]],
  !>    which opens a file for writing.
  !> An [[OFile(type)]] can be copied, creating multiple [[OFile(type)]]
  !>    instances which all share a single [[OFileTarget(type)]] and thus a
  !>    single file.
  !> When the last [[OFile(type)]] which shares a given [[OFileTarget(type)]]
  !>    is finalised, the [[OFileTarget(type)]] is also finalised and the file
  !>    is closed.
  type :: OFile
    !> Handles opening and closing the file, and writing to the file.
    type(OFileTarget),   pointer     :: ofile_target
    !> Counts how many extant [[OFile(type)]] instances share `ofile_target`.
    type(SharedCounter), allocatable :: counter
  contains
    final :: final_OFile
    
    procedure, private :: check_associated
    
    procedure, public :: make_stdout
    
    generic,   public  :: print_line =>               &
                        & print_line_character,       &
                        & print_line_String,          &
                        & print_line_StringWriteable, &
                        & print_line_logical,         &
                        & print_line_integer,         &
                        & print_line_real,            &
                        & print_line_complex,         &
                        & print_line_logicals,        &
                        & print_line_integers,        &
                        & print_line_reals,           &
                        & print_line_complexes
    procedure, private :: print_line_character
    procedure, private :: print_line_String
    procedure, private :: print_line_StringWriteable
    procedure, private :: print_line_logical
    procedure, private :: print_line_integer
    procedure, private :: print_line_real
    procedure, private :: print_line_complex
    procedure, private :: print_line_logicals
    procedure, private :: print_line_integers
    procedure, private :: print_line_reals
    procedure, private :: print_line_complexes
    
    generic,   public  :: print_lines =>                           &
                        & print_lines_Strings_character,           &
                        & print_lines_Strings_String,              &
                        & print_lines_StringWriteables_character,  &
                        & print_lines_StringWriteables_String,     &
                        & print_lines_StringsWriteable,            &
                        & print_lines_StringsWriteables_character, &
                        & print_lines_StringsWriteables_String,    &
                        & print_lines_logicals_character,          &
                        & print_lines_logicals_String,             &
                        & print_lines_integers_character,          &
                        & print_lines_integers_String,             &
                        & print_lines_reals_character,             &
                        & print_lines_reals_String,                &
                        & print_lines_complexes_character,         &
                        & print_lines_complexes_String
    procedure, private :: print_lines_Strings_character
    procedure, private :: print_lines_Strings_String
    procedure, private :: print_lines_StringWriteables_character
    procedure, private :: print_lines_StringWriteables_String
    procedure, private :: print_lines_StringsWriteable
    procedure, private :: print_lines_StringsWriteables_character
    procedure, private :: print_lines_StringsWriteables_String
    procedure, private :: print_lines_logicals_character
    procedure, private :: print_lines_logicals_String
    procedure, private :: print_lines_integers_character
    procedure, private :: print_lines_integers_String
    procedure, private :: print_lines_reals_character
    procedure, private :: print_lines_reals_String
    procedure, private :: print_lines_complexes_character
    procedure, private :: print_lines_complexes_String
  end type
  
  interface OFile
    !> Constructor. Opens the file specified by `filename` for writing.
    module function new_OFile_character(filename) result(this) 
      character(*), intent(in) :: filename
      type(OFile)              :: this
    end function

    !> Constructor. Opens the file specified by `filename` for writing.
    module function new_OFile_String(filename) result(this) 
      type(String), intent(in) :: filename
      type(OFile)              :: this
    end function
  end interface
  
  interface assignment(=)
    !> Copies one instance of [[OFile(type)]] to another, sharing a single
    !>    [[OFileTarget(type)]] between the two.
    module subroutine assign_OFile_OFile(output,input) 
      type(OFile), intent(out) :: output
      type(OFile), intent(in)  :: input
    end subroutine
  end interface
  
  interface
    !> Finalisation. Decrements `this%counter`, and finalises
    !>    `this%ofile_target` if there are no extant [[OFile(type)]] instances
    !>    which share `this%ofile_target`.
    module subroutine final_OFile(this) 
      type(OFile), intent(inout) :: this
    end subroutine

    !> Checks that this instance has been associated with a file.
    module subroutine check_associated(this) 
      class(OFile), intent(in) :: this
    end subroutine

    !> Makes the file the target of
    !>    [[caesar_print_module:print_line(interface)]] and
    !>    [[caesar_print_module:print_lines(interface)]] statements.
    module subroutine make_stdout(this) 
      class(OFile), intent(inout) :: this
    end subroutine

    !> Prints `input` to the file as a single line followed by a line break.
    module subroutine print_line_character(this,input,settings) 
      class(OFile),        intent(inout)        :: this
      character(*),        intent(in)           :: input
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file as a single line followed by a line break.
    module subroutine print_line_String(this,input,settings) 
      class(OFile),        intent(inout)        :: this
      type(String),        intent(in)           :: input
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file as a single line followed by a line break.
    module subroutine print_line_StringWriteable(this,input,settings) 
      class(OFile),           intent(inout)        :: this
      class(StringWriteable), intent(in)           :: input
      !> Controls formatting.
      type(PrintSettings),    intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file as a single line followed by a line break.
    module subroutine print_line_logical(this,input,settings) 
      class(OFile),        intent(inout)        :: this
      logical,             intent(in)           :: input
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file as a single line followed by a line break.
    module subroutine print_line_integer(this,input,settings) 
      class(OFile),        intent(inout)        :: this
      integer,             intent(in)           :: input
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file as a single line followed by a line break.
    module subroutine print_line_real(this,input,settings) 
      class(OFile),        intent(inout)        :: this
      real(dp),            intent(in)           :: input
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file as a single line followed by a line break.
    module subroutine print_line_complex(this,input,settings) 
      class(OFile),        intent(inout)        :: this
      complex(dp),         intent(in)           :: input
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file as a single line followed by a line break.
    module subroutine print_line_logicals(this,input,settings) 
      class(OFile),        intent(inout)        :: this
      logical,             intent(in)           :: input(:)
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file as a single line followed by a line break.
    module subroutine print_line_integers(this,input,settings) 
      class(OFile),        intent(inout)        :: this
      integer,             intent(in)           :: input(:)
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file as a single line followed by a line break.
    module subroutine print_line_reals(this,input,settings) 
      class(OFile),        intent(inout)        :: this
      real(dp),            intent(in)           :: input(:)
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file as a single line followed by a line break.
    module subroutine print_line_complexes(this,input,settings) 
      class(OFile),        intent(inout)        :: this
      complex(dp),         intent(in)           :: input(:)
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file, as multiple lines.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module subroutine print_lines_Strings_character(this,input, &
       & separating_line,settings) 
      class(OFile),        intent(inout)        :: this
      type(String),        intent(in)           :: input(:)
      character(*),        intent(in), optional :: separating_line
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file, as multiple lines.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module subroutine print_lines_Strings_String(this,input,separating_line, &
       & settings) 
      class(OFile),        intent(inout)        :: this
      type(String),        intent(in)           :: input(:)
      type(String),        intent(in)           :: separating_line
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file, as multiple lines.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module subroutine print_lines_StringWriteables_character(this,input, &
       & separating_line,settings) 
      class(OFile),           intent(inout)        :: this
      class(StringWriteable), intent(in)           :: input(:)
      character(*),           intent(in), optional :: separating_line
      !> Controls formatting.
      type(PrintSettings),    intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file, as multiple lines.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module subroutine print_lines_StringWriteables_String(this,input, &
       & separating_line,settings) 
      class(OFile),           intent(inout)        :: this
      class(StringWriteable), intent(in)           :: input(:)
      type(String),           intent(in)           :: separating_line
      !> Controls formatting.
      type(PrintSettings),    intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file, as multiple lines.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module subroutine print_lines_StringsWriteable(this,input,settings) 
      class(OFile),            intent(inout)        :: this
      class(StringsWriteable), intent(in)           :: input
      !> Controls formatting.
      type(PrintSettings),     intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file, as multiple lines.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module subroutine print_lines_StringsWriteables_String(this,input, &
       & separating_line,settings) 
      class(OFile),            intent(inout)        :: this
      class(StringsWriteable), intent(in)           :: input(:)
      type(String),            intent(in), optional :: separating_line
      !> Controls formatting.
      type(PrintSettings),     intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file, as multiple lines.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module subroutine print_lines_StringsWriteables_character(this,input, &
       & separating_line,settings) 
      class(OFile),            intent(inout)        :: this
      class(StringsWriteable), intent(in)           :: input(:)
      character(*),            intent(in)           :: separating_line
      !> Controls formatting.
      type(PrintSettings),     intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file, as multiple lines.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module subroutine print_lines_logicals_character(this,input, &
       & separating_line,settings) 
      class(OFile),        intent(inout)        :: this
      logical,             intent(in)           :: input(:)
      character(*),        intent(in), optional :: separating_line
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file, as multiple lines.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module subroutine print_lines_logicals_String(this,input, &
       & separating_line,settings) 
      class(OFile),        intent(inout)        :: this
      logical,             intent(in)           :: input(:)
      type(String),        intent(in)           :: separating_line
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file, as multiple lines.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module subroutine print_lines_integers_character(this,input, &
       & separating_line,settings) 
      class(OFile),        intent(inout)        :: this
      integer,             intent(in)           :: input(:)
      character(*),        intent(in), optional :: separating_line
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file, as multiple lines.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module subroutine print_lines_integers_String(this,input, &
       & separating_line,settings) 
      class(OFile),        intent(inout)        :: this
      integer,             intent(in)           :: input(:)
      type(String),        intent(in)           :: separating_line
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file, as multiple lines.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module subroutine print_lines_reals_character(this,input, &
       & separating_line,settings) 
      class(OFile),        intent(inout)        :: this
      real(dp),            intent(in)           :: input(:)
      character(*),        intent(in), optional :: separating_line
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file, as multiple lines.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module subroutine print_lines_reals_String(this,input,separating_line, &
       & settings) 
      class(OFile),        intent(inout)        :: this
      real(dp),            intent(in)           :: input(:)
      type(String),        intent(in)           :: separating_line
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file, as multiple lines.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module subroutine print_lines_complexes_character(this,input, &
       & separating_line,settings) 
      class(OFile),        intent(inout)        :: this
      complex(dp),         intent(in)           :: input(:)
      character(*),        intent(in), optional :: separating_line
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints `input` to the file, as multiple lines.
    !> If `separating_line` is given then it will be inserted between
    !>    each line of the output.
    module subroutine print_lines_complexes_String(this,input, &
       & separating_line,settings) 
      class(OFile),        intent(inout)        :: this
      complex(dp),         intent(in)           :: input(:)
      type(String),        intent(in)           :: separating_line
      !> Controls formatting.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine
  end interface
end module
