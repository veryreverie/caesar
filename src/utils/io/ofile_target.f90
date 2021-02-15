!> Provides the [[OFileTarget(type)]] class, which is used by [[OFile(type)]].
module caesar_ofile_target_module
  use caesar_foundations_module
  
  use caesar_string_module
  use caesar_print_settings_module
  use caesar_string_writeable_module
  use caesar_strings_writeable_module
  use caesar_file_module
  implicit none
  
  private
  
  public :: OFileTarget
  
  !> Keeps track of an output file.
  !> The file is opened when the [[OFileTarget(type)]] is constructed,
  !>    and closed when the [[OFileTarget(type)]] is finalised.
  type :: OFileTarget
    !> Whether or not the file is open.
    logical,      private :: open_ = .false.
    !> The file name.
    type(String), private :: filename_
    !> The file unit associated with the file.
    integer,      private :: file_unit_
    !> Tracks whether or not this file has been set as the write location
    !>    for [[caesar_print_module:print_line(interface)]] and
    !>    [[caesar_print_module:print_lines(interface)]] calls.
    logical,      private :: is_stdout_
  contains
    procedure, public :: close => close_OFileTarget
    
    procedure, public :: is_open
    
    procedure, public :: make_stdout
    
    generic,   public  :: print_line => &
                        & print_line_character
    procedure, private :: print_line_character
  end type
  
  interface OFileTarget
    !> Constructor. Opens the file specified by `filename` for writing.
    module function new_OFileTarget(filename) result(this) 
      character(*), intent(in) :: filename
      type(OFileTarget)        :: this
    end function
  end interface
  
  interface
    !> Finalisation. Closes the file, and redirects the write location of
    !>    [[caesar_print_module:print_line(interface)]] and
    !>    [[caesar_print_module:print_lines(interface)]] if relevant.
    module subroutine close_OFileTarget(this) 
      class(OFileTarget), intent(inout) :: this
    end subroutine

    !> Returns whether or not this file is open.
    module function is_open(this) result(output) 
      class(OFileTarget), intent(in) :: this
      logical                        :: output
    end function

    !> Makes this file the write location for
    !>    [[caesar_print_module:print_line(interface)]] and
    !>    [[caesar_print_module:print_lines(interface)]].
    module subroutine make_stdout(this) 
      class(OFileTarget), intent(inout) :: this
    end subroutine

    !> Writes a line to the file.
    module subroutine print_line_character(this,input,settings) 
      class(OFileTarget),  intent(inout)        :: this
      type(PrintSettings), intent(in), optional :: settings
      character(*),        intent(in)           :: input
    end subroutine
  end interface
end module
