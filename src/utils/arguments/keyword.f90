! ======================================================================
! Keywords for input arguments.
! ======================================================================
module caesar_keyword_module
  use caesar_foundations_module
  use caesar_abstract_module
  use caesar_io_module
  implicit none
  
  private
  
  public :: KeywordData
  
  ! A keyword for getting input arguments from the user.
  type KeywordData
    ! Keyword.
    type(String), private :: keyword_
    
    ! Flag (alternative one-character keyword).
    ! 0=no flag,1=flag without arguments,2=flag which takes arguments.
    integer,      private :: flag_type_
    character(1), private :: flag_
    
    ! Helptext for help calls.
    type(String), private :: helptext_
    
    ! Properties of the value.
    logical, private :: is_path_
    logical, private :: allowed_in_file_
    logical, private :: can_be_interactive_
    logical, private :: pass_to_python_
    
    ! Defaults.
    ! 0=must be set,1=no default,2=default to value,3=default to keyword.
    integer,      private :: default_type_
    type(String), private :: default_
    
    ! The value itself.
    logical,      private :: is_set_
    type(String), private :: value_
    
    ! Dependencies on other keywords.
    type(String), allocatable, private :: exclusive_with_(:)
    
    ! If the keyword is a path and is set from a file, then relative paths
    !    are taken to be relative to that file.
    ! The file path must then be stored.
    type(String), allocatable, private :: working_directory_
  contains
    ! Flag-related procedures.
    procedure, public :: has_flag            => has_flag_KeywordData
    procedure, public :: flag_takes_argument => flag_takes_argument_KeywordData
    procedure, public :: flag                => flag_KeywordData
    procedure, public :: set_flag            => set_flag_KeywordData
    
    ! Default-related procedures.
    procedure, public :: defaults_to_keyword => &
                       & defaults_to_keyword_KeywordData
    procedure, public :: set_default => set_default_KeywordData
    
    ! Setters.
    procedure, public  :: unset   => unset_KeywordData
    
    generic,   public  :: set =>                     &
                        & set_KeywordData_character, &
                        & set_KeywordData_String
    procedure, private :: set_KeywordData_character
    procedure, private :: set_KeywordData_String
    
    generic,   public  :: append =>                     &
                        & append_KeywordData_character, &
                        & append_KeywordData_String
    procedure, private :: append_KeywordData_character
    procedure, private :: append_KeywordData_String
    
    ! Getters.
    procedure, public :: keyword => keyword_KeywordData
    procedure, public :: is_set  => is_set_KeywordData
    procedure, public :: value   => value_KeywordData
    
    procedure, public :: is_path            => is_path_KeywordData
    procedure, public :: allowed_in_file    => allowed_in_file_KeywordData
    procedure, public :: can_be_interactive => can_be_interactive_KeywordData
    procedure, public :: pass_to_python     => pass_to_python_KeywordData
    
    procedure, public :: exclusive_with => exclusive_with_KeywordData
    
    ! Set interactively from user.
    procedure, public :: set_interactively
    
    ! Processes and checks value.
    procedure, public :: process_and_check
    
    ! Prints helptext.
    procedure, public :: print_help
  end type
  
  interface
    ! ----------------------------------------------------------------------
    ! Flag-related procedures.
    ! ----------------------------------------------------------------------
    module function has_flag_KeywordData(this) result(output) 
      class(KeywordData), intent(in) :: this
      logical                        :: output
    end function
  end interface
  
  interface
    module function flag_takes_argument_KeywordData(this) result(output) 
      class(KeywordData), intent(in) :: this
      logical                        :: output
    end function
  end interface
  
  interface
    module function flag_KeywordData(this) result(output) 
      class(KeywordData), intent(in) :: this
      character(1)                   :: output
    end function
  end interface
  
  interface
    module subroutine set_flag_KeywordData(this,flag,flag_takes_arguments) 
      class(KeywordData), intent(inout) :: this
      character(1),       intent(in)    :: flag
      logical,            intent(in)    :: flag_takes_arguments
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Default-related procedures.
    ! ----------------------------------------------------------------------
    
    ! If the keyword defaults to another keyword, returns that keyword.
    ! Returns '' if the keyword does not default to a keyword.
    module function defaults_to_keyword_KeywordData(this) result(output) 
      class(KeywordData), intent(in) :: this
      type(String)                   :: output
    end function
  end interface
  
  interface
    ! Sets an unset keyword to its default value.
    ! Does nothing if the keyword is set or has no default.
    ! Throws an error if the keyword defaults to another keyword.
    module subroutine set_default_KeywordData(this) 
      class(KeywordData), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Setters.
    ! ----------------------------------------------------------------------
    module subroutine unset_KeywordData(this) 
      class(KeywordData), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    module subroutine set_KeywordData_character(this,value, &
       & only_update_if_unset,working_directory) 
      class(KeywordData), intent(inout)        :: this
      character(*),       intent(in)           :: value
      logical,            intent(in), optional :: only_update_if_unset
      type(String),       intent(in), optional :: working_directory
    end subroutine
  end interface
  
  interface
    module subroutine set_KeywordData_String(this,value, &
       & only_update_if_unset,working_directory) 
      class(KeywordData), intent(inout)        :: this
      type(String),       intent(in)           :: value
      logical,            intent(in), optional :: only_update_if_unset
      type(String),       intent(in), optional :: working_directory
    end subroutine
  end interface
  
  interface
    module subroutine append_KeywordData_character(this,value) 
      class(KeywordData), intent(inout)        :: this
      character(*),       intent(in)           :: value
    end subroutine
  end interface
  
  interface
    module subroutine append_KeywordData_String(this,value) 
      class(KeywordData), intent(inout)        :: this
      type(String),       intent(in)           :: value
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Getters.
    ! ----------------------------------------------------------------------
    impure elemental module function keyword_KeywordData(this) result(output) 
      class(KeywordData), intent(in) :: this
      type(String)                   :: output
    end function
  end interface
  
  interface
    impure elemental module function is_set_KeywordData(this) result(output) 
      class(KeywordData), intent(in) :: this
      logical                        :: output
    end function
  end interface
  
  interface
    impure elemental module function value_KeywordData(this) result(output) 
      class(KeywordData), intent(in) :: this
      type(String)                   :: output
    end function
  end interface
  
  interface
    module function is_path_KeywordData(this) result(output) 
      class(KeywordData), intent(in) :: this
      logical                        :: output
    end function
  end interface
  
  interface
    module function allowed_in_file_KeywordData(this) result(output) 
      class(KeywordData), intent(in) :: this
      logical                        :: output
    end function
  end interface
  
  interface
    module function can_be_interactive_KeywordData(this) result(output) 
      class(KeywordData), intent(in) :: this
      logical                        :: output
    end function
  end interface
  
  interface
    module function pass_to_python_KeywordData(this) result(output) 
      class(KeywordData), intent(in) :: this
      logical                        :: output
    end function
  end interface
  
  interface
    module function exclusive_with_KeywordData(this) result(output) 
      class(KeywordData), intent(in) :: this
      type(String), allocatable      :: output(:)
    end function
  end interface
  
  interface KeywordData
    ! ----------------------------------------------------------------------
    ! Takes a keyword, helptext and options and returns a KeywordData.
    ! ----------------------------------------------------------------------
    ! keyword is the name of the keyword, by which it can be accesed both in the
    !    code and by the user.
    ! helptext is the text which will print for the keyword when help is called
    !    for.
    ! default_value is the value which will be defaulted to if the keyword is not
    !    specified. Defaults to not set.
    ! default_keyword is the keyword whose value will be copied if the keyword is
    !    not specified. Defaults to not set.
    ! N.B. At most one of default_value and default_keyword should be specified.
    ! is_optional specifies whether or not the keyword is required to be given.
    !    is_optional defaults to .false. unless a default or exclusive_with is
    !    given, when it defaults to, and is required to be, .true..
    ! is_path specifies that the value will be a path to a file or directory.
    !    Paths are converted to absolute paths (from /) automatically.
    !    Defaults to .false..
    ! allowed_in_file specifies that the value may be read from or printed to file.
    !    Defaults to .true..
    ! can_be_interactive specifies that the value may be set interactively by the
    !    user.
    !    Defaults to .true.
    ! flag specifies the flag by which the keyword can be alternately called.
    !    Defaults to ' '.
    ! pass_to_python specifies that the keyword should be passed to python.
    ! exclusive_with specifies which other keywords this keyword is mutually
    !    exclusive with. The code guarantees that at most one of each mutually
    !    exclusive set of keywords is set. N.B. it is not guaranteed that at least
    !    one of each set is set; they could all be unset.
    module function new_KeywordData(keyword,helptext,default_value,     &
       & default_keyword,is_optional,is_path,allowed_in_file,           &
       & can_be_interactive,flag_without_arguments,flag_with_arguments, &
       & pass_to_python,exclusive_with) result(this) 
      character(*), intent(in)           :: keyword
      character(*), intent(in)           :: helptext
      character(*), intent(in), optional :: default_value
      character(*), intent(in), optional :: default_keyword
      logical,      intent(in), optional :: is_optional
      logical,      intent(in), optional :: is_path
      logical,      intent(in), optional :: allowed_in_file
      logical,      intent(in), optional :: can_be_interactive
      character(1), intent(in), optional :: flag_without_arguments
      character(1), intent(in), optional :: flag_with_arguments
      logical,      intent(in), optional :: pass_to_python
      type(String), intent(in), optional :: exclusive_with(:)
      type(KeywordData)                  :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Sets the keyword interactively, asking the user for a value.
    ! ----------------------------------------------------------------------
    recursive module subroutine set_interactively(this) 
      class(KeywordData), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Process and check value.
    ! ----------------------------------------------------------------------
    ! Throws and error if a value is required but has not been set.
    ! Turns paths into absolute form.
    module subroutine process_and_check(this) 
      class(KeywordData), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Prints helptext and relevant keyword settings.
    ! ----------------------------------------------------------------------
    module subroutine print_help(this) 
      class(KeywordData), intent(in) :: this
    end subroutine
  end interface
end module
