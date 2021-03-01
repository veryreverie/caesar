! ======================================================================
! A dictionary of keys and values, both of type String.
! ======================================================================
module caesar_dictionary_module
  use caesar_foundations_module
  use caesar_abstract_module
  use caesar_io_module
  
  use caesar_keyword_module
  use caesar_common_keywords_module
  implicit none
  
  private
  
  public :: Dictionary
  public :: size
  public :: call_caesar
  
  ! ------------------------------
  ! A dictionary of keys and values.
  ! ------------------------------
  type, extends(NoDefaultConstructor) :: Dictionary
    ! Accepted keywords.
    type(KeywordData), allocatable, private :: keywords_(:)
  contains
    ! ----------
    ! Getters.
    ! ----------
    ! Private function to return the index of a keyword.
    generic,   private :: index =>                 &
                        & index_Dictionary_String, &
                        & index_Dictionary_character
    procedure, private :: index_Dictionary_String
    procedure, private :: index_Dictionary_character
    
    ! As above, but by flag rather than by keyword.
    generic,   private :: index_by_flag =>                    &
                        & index_by_flag_Dictionary_character, &
                        & index_by_flag_Dictionary_String
    procedure, private :: index_by_flag_Dictionary_character
    procedure, private :: index_by_flag_Dictionary_String
    
    ! Returns the corresponding keyword.
    generic,   private :: flag_to_keyword =>                    &
                        & flag_to_keyword_Dictionary_character, &
                        & flag_to_keyword_Dictionary_String
    procedure, private :: flag_to_keyword_Dictionary_character
    procedure, private :: flag_to_keyword_Dictionary_String
    
    ! Returns whether or not a keyword is set with a value.
    generic,   public  :: is_set =>                    &
                        & is_set_Dictionary_character, &
                        & is_set_Dictionary_String
    procedure, private :: is_set_Dictionary_character
    procedure, private :: is_set_Dictionary_String
    
    ! Returns the value matching a keyword.
    generic,   public  :: value =>                    &
                        & value_Dictionary_character, &
                        & value_Dictionary_String
    procedure, private :: value_Dictionary_character
    procedure, private :: value_Dictionary_String
    
    ! Returns a list of python arguments.
    procedure, public :: python_arguments
    
    ! ----------
    ! Setters.
    ! ----------
    ! Unsets a keyword.
    generic,   public  :: unset =>                    &
                        & unset_Dictionary_character, &
                        & unset_Dictionary_String
    procedure, private :: unset_Dictionary_character
    procedure, private :: unset_Dictionary_String
  
    ! Sets a keyword and sets its value.
    generic,   public  :: set =>                              &
                        & set_Dictionary_character_character, &
                        & set_Dictionary_character_String,    &
                        & set_Dictionary_String_character,    &
                        & set_Dictionary_String_String,       &
                        & set_Dictionary_Dictionary
    procedure, private :: set_Dictionary_character_character
    procedure, private :: set_Dictionary_character_String
    procedure, private :: set_Dictionary_String_character
    procedure, private :: set_Dictionary_String_String
    procedure, private :: set_Dictionary_Dictionary
    
    ! Appends to the value of a keyword.
    ! Returns an error if the keyword has no value set.
    generic,   private :: append_to_value =>                              &
                        & append_to_value_Dictionary_character_character, &
                        & append_to_value_Dictionary_character_String,    &
                        & append_to_value_Dictionary_String_character,    &
                        & append_to_value_Dictionary_String_String
    procedure, private :: append_to_value_Dictionary_character_character
    procedure, private :: append_to_value_Dictionary_character_String
    procedure, private :: append_to_value_Dictionary_String_character
    procedure, private :: append_to_value_Dictionary_String_String
    
    ! ----------
    ! File operations.
    ! ----------
    ! Write to a file.
    generic,   public  :: write_file =>                    &
                        & write_file_Dictionary_character, &
                        & write_file_Dictionary_String
    procedure, private :: write_file_Dictionary_character
    procedure, private :: write_file_Dictionary_String
  
    ! Read keywords from a file.
    generic,   public  :: read_file =>                    &
                        & read_file_Dictionary_character, &
                        & read_file_Dictionary_String
    procedure, private :: read_file_Dictionary_character
    procedure, private :: read_file_Dictionary_String
  
    
    ! ----------
    ! Set defaults, process paths and check non-optional keywords are set.
    ! ----------
    procedure, public  :: set_interactively => set_interactively_Dictionary
    procedure, private :: process_and_check_inputs => &
                        & process_and_check_inputs_Dictionary
  end type
  
  interface Dictionary
    ! ----------------------------------------------------------------------
    ! Constructors.
    ! ----------------------------------------------------------------------
    ! Construct a Dictionary from and array of KeywordData.
    module function new_Dictionary_KeywordDatas(keywords) result(this) 
      type(KeywordData), intent(in) :: keywords(:)
      type(Dictionary)              :: this
    end function
  
    ! Construct a Dictionary from the command line input to Caesar.
    module function new_Dictionary_arguments(args,keywords_in) &
       & result(arguments) 
      type(String),      intent(in) :: args(:)
      type(KeywordData), intent(in) :: keywords_in(:)
      type(Dictionary)              :: arguments
    end function
  end interface
  
  interface size
    ! ----------------------------------------------------------------------
    ! size(Dictionary).
    ! ----------------------------------------------------------------------
    module function size_Dictionary(this) result(output) 
      type(Dictionary), intent(in) :: this
      integer                      :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Get the index where the keyword is stored.
    ! ----------------------------------------------------------------------
    ! Private module function.
    ! Throws an error if the keyword is not found.
    ! If there are duplicate keys, returns the first match.
    impure elemental module function index_Dictionary_String(this,keyword) &
       & result(output) 
      class(Dictionary), intent(in) :: this
      type(String),      intent(in) :: keyword
      integer                       :: output
    end function
  end interface
  
  interface
    module function index_Dictionary_character(this,keyword) result(output) 
      class(Dictionary), intent(in) :: this
      character(*),      intent(in) :: keyword
      integer                       :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! As above, but by flag rather than keyword.
    ! ----------------------------------------------------------------------
    module function index_by_flag_Dictionary_character(this,flag) &
       & result(output) 
      class(Dictionary), intent(in) :: this
      character(1),      intent(in) :: flag
      integer                       :: output
    end function
  end interface
  
  interface
    impure elemental module function index_by_flag_Dictionary_String(this, &
       & flag) result(output) 
      class(Dictionary), intent(in) :: this
      type(String),      intent(in) :: flag
      integer                       :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    !  Returns the keyword corresponding to the given flag.
    ! ----------------------------------------------------------------------
    module function flag_to_keyword_Dictionary_character(this,flag) &
       & result(output) 
      class(Dictionary), intent(in) :: this
      character(1),      intent(in) :: flag
      type(String)                  :: output
    end function
  end interface
  
  interface
    module function flag_to_keyword_Dictionary_String(this,flag) &
       & result(output) 
      class(Dictionary), intent(in) :: this
      type(String),      intent(in) :: flag
      type(String)                  :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Get whether or not a keyword has been set with a value.
    ! ----------------------------------------------------------------------
    module function is_set_Dictionary_character(this,keyword) result(output) 
      class(Dictionary), intent(in) :: this
      character(*),      intent(in) :: keyword
      logical                       :: output
    end function
  end interface
  
  interface
    module function is_set_Dictionary_String(this,keyword) result(output) 
      class(Dictionary), intent(in) :: this
      type(String),      intent(in) :: keyword
      logical                       :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Get the value corresponding to a given key.
    ! ----------------------------------------------------------------------
    ! Throws an error if the key has not been set or is boolean.
    module function value_Dictionary_character(this,keyword) result(output) 
      class(Dictionary), intent(in) :: this
      character(*),      intent(in) :: keyword
      type(String)                  :: output
    end function
  end interface
  
  interface
    module function value_Dictionary_String(this,keyword) result(output) 
      class(Dictionary), intent(in) :: this
      type(String),      intent(in) :: keyword
      type(String)                  :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Return a list of python arguments.
    ! ----------------------------------------------------------------------
    module function python_arguments(this) result(output) 
      class(Dictionary), intent(in) :: this
      type(String), allocatable     :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Unsets a given keyword.
    ! ----------------------------------------------------------------------
    module subroutine unset_Dictionary_character(this,keyword) 
      class(Dictionary), intent(inout) :: this
      character(*),      intent(in)    :: keyword
    end subroutine
  end interface
  
  interface
    module subroutine unset_Dictionary_String(this,keyword) 
      class(Dictionary), intent(inout) :: this
      type(String),      intent(in)    :: keyword
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Sets the value corresponding to a given keyword.
    ! ----------------------------------------------------------------------
    ! If only_update_if_unset is set, the value will not be overwritten if set.
    !    defaults to .false..
    module subroutine set_Dictionary_character_character(this,keyword,value, &
       & only_update_if_unset) 
      class(Dictionary), intent(inout)        :: this
      character(*),      intent(in)           :: keyword
      character(*),      intent(in)           :: value
      logical,           intent(in), optional :: only_update_if_unset
    end subroutine
  end interface
  
  interface
    module subroutine set_Dictionary_character_String(this,keyword,value, &
       & only_update_if_unset) 
      class(Dictionary), intent(inout)        :: this
      character(*),      intent(in)           :: keyword
      type(String),      intent(in)           :: value
      logical,           intent(in), optional :: only_update_if_unset
    end subroutine
  end interface
  
  interface
    module subroutine set_Dictionary_String_character(this,keyword,value, &
       & only_update_if_unset) 
      class(Dictionary), intent(inout)        :: this
      type(String),      intent(in)           :: keyword
      character(*),      intent(in)           :: value
      logical,           intent(in), optional :: only_update_if_unset
    end subroutine
  end interface
  
  interface
    module subroutine set_Dictionary_String_String(this,keyword,value, &
       & only_update_if_unset) 
      class(Dictionary), intent(inout)        :: this
      type(String),      intent(in)           :: keyword
      type(String),      intent(in)           :: value
      logical,           intent(in), optional :: only_update_if_unset
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Sets all values from another dictionary.
    ! ----------------------------------------------------------------------
    module subroutine set_Dictionary_Dictionary(this,that) 
      class(Dictionary), intent(inout) :: this
      type(Dictionary),  intent(in)    :: that
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Appends to the value corresponding to a given keyword.
    ! ----------------------------------------------------------------------
    module subroutine append_to_value_Dictionary_character_character(this, &
       & keyword,value) 
      class(Dictionary), intent(inout) :: this
      character(*),      intent(in)    :: keyword
      character(*),      intent(in)    :: value
    end subroutine
  end interface
  
  interface
    module subroutine append_to_value_Dictionary_character_String(this, &
       & keyword,value) 
      class(Dictionary), intent(inout) :: this
      character(*),      intent(in)    :: keyword
      type(String),      intent(in)    :: value
    end subroutine
  end interface
  
  interface
    module subroutine append_to_value_Dictionary_String_character(this, &
       & keyword,value) 
      class(Dictionary), intent(inout) :: this
      type(String),      intent(in)    :: keyword
      character(*),      intent(in)    :: value
    end subroutine
  end interface
  
  interface
    module subroutine append_to_value_Dictionary_String_String(this,keyword, &
       & value) 
      class(Dictionary), intent(inout) :: this
      type(String),      intent(in)    :: keyword
      type(String),      intent(in)    :: value
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Writes a Dictionary to file.
    ! ----------------------------------------------------------------------
    module subroutine write_file_Dictionary_character(this,filename) 
      class(Dictionary), intent(in) :: this
      character(*),      intent(in) :: filename
    end subroutine
  end interface
  
  interface
    module subroutine write_file_Dictionary_String(this,filename) 
      class(Dictionary), intent(in) :: this
      type(String),      intent(in) :: filename
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Reads Dictionary entries from a file.
    ! ----------------------------------------------------------------------
    ! The dictionary must already have been initialised from a list of keywords
    !    before this is called.
    ! The contents of the file will be checked against the dictionary's keywords.
    ! If only_update_if_unset is set, then only keywords which have not already
    !    been set will be modified. This defaults to .false..
    module subroutine read_file_Dictionary_character(this,filename, &
       & only_update_if_unset) 
      class(Dictionary), intent(inout)        :: this
      character(*),      intent(in)           :: filename
      logical,           intent(in), optional :: only_update_if_unset
    end subroutine
  end interface
  
  interface
    module subroutine read_file_Dictionary_String(this,filename, &
       & only_update_if_unset) 
      class(Dictionary), intent(inout)        :: this
      type(String),      intent(in)           :: filename
      logical,           intent(in), optional :: only_update_if_unset
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Set defaults for keywords which have not been set.
    ! Convert all paths to absolute format (from /).
    ! Check that all non-optional keywords have been set.
    ! ----------------------------------------------------------------------
    module subroutine set_interactively_Dictionary(this) 
      class(Dictionary), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    module subroutine process_and_check_inputs_Dictionary(this) 
      class(Dictionary), intent(inout) :: this
    end subroutine
  end interface
  
  interface call_caesar
    ! ----------------------------------------------------------------------
    ! Calls Caesar with a given dictionary of arguments.
    ! ----------------------------------------------------------------------
    module subroutine call_caesar_String_Dictionary(mode,arguments) 
      type(String),     intent(in) :: mode
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
  
  interface call_caesar
    module subroutine call_caesar_character_Dictionary(mode,arguments) 
      character(*),     intent(in) :: mode
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module
