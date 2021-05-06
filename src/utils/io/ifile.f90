!> Provides the [[IFile(type)]] class, which handles input file reading.
module caesar_ifile_module
  use caesar_foundations_module
  
  use caesar_string_module
  use caesar_file_module
  use caesar_string_array_module
  implicit none
  
  private
  
  public :: IFile
  public :: size
  
  !> Handles input file reading.
  !> Constructing this class reads a file and stores its contents.
  type :: IFile
    !> The file name.
    type(String), private              :: filename_
    !> The contents of the file, split by line breaks.
    type(String), private, allocatable :: lines_(:)
  contains
    procedure, public :: line
    
    generic,   public  :: lines =>   &
                        & lines_all, &
                        & lines_slice
    procedure, private :: lines_all
    procedure, private :: lines_slice
    
    generic,   public  :: sections =>         &
                        & sections_character, &
                        & sections_String
    procedure, private :: sections_character
    procedure, private :: sections_String
  end type
  
  interface IFile
    !> Read the file with the specified `filename`, and store the contents of
    !>    that file within `this`.
    !> N.B. the file is closed immediately after reading.
    module function new_IFile_character(filename) result(this) 
      character(*), intent(in) :: filename
      type(IFile)              :: this
    end function

    !> Read the file with the specified `filename`, and store the contents of
    !>    that file within `this`.
    !> N.B. the file is closed immediately after reading.
    module function new_IFile_String(filename) result(this) 
      type(String), intent(in) :: filename
      type(IFile)              :: this
    end function
  end interface
  
  interface
    !> Returns the line from the file with the specified `line_number`.
    module function line(this,line_number) result(output) 
      class(IFile), intent(in) :: this
      integer,      intent(in) :: line_number
      type(String)             :: output
    end function

    !> Returns a [[String(type)]] array containing every line in the file.
    module function lines_all(this) result(output) 
      class(IFile), intent(in)  :: this
      type(String), allocatable :: output(:)
    end function

    !> Returns a [[String(type)]] array containing the lines with line numbers
    !>    between `first_line_number` and `last_line_number` inclusive.
    module function lines_slice(this,first_line_number,last_line_number) &
       & result(output) 
      class(IFile), intent(in)  :: this
      integer,      intent(in)  :: first_line_number
      integer,      intent(in)  :: last_line_number
      type(String), allocatable :: output(:)
    end function

    !> Returns the contents of the file as a `[[StringArray(type)]]` array,
    !>    split into sections by `separating_line`.
    !> `separating_line` defaults to an empty string.
    module function sections_String(this,separating_line) result(output) 
      class(IFile), intent(in)           :: this
      type(String), intent(in), optional :: separating_line
      type(StringArray), allocatable     :: output(:)
    end function

    !> Returns the sections of the file, using
    !>    [[split_into_sections(procedure)]].
    module function sections_character(this,separating_line) result(output) 
      class(IFile), intent(in)       :: this
      character(*), intent(in)       :: separating_line
      type(StringArray), allocatable :: output(:)
    end function
  end interface
  
  interface count_lines
    !> Returns the number of lines in the file specified by `filename`.
    module function count_lines(filename) result(output) 
      character(*), intent(in) :: filename
      integer                  :: output
    end function
  end interface
  
  interface size
    !> Returns the number of lines in the file.
    module function size_IFile(this) result(output) 
      type(IFile), intent(in) :: this
      integer                 :: output
    end function
  end interface
end module  
