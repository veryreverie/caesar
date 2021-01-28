!> Provides the [[StringBase(type)]] class.
module caesar_string_base_module
  use caesar_error_module
  implicit none
  
  private
  
  public :: StringBase
  public :: char
  
  !> A basic allocatable string class.
  !> This class exists to provide an interface between the actual contents of
  !>    a String and all of the procedures involving a String.
  !> This ensures that all of the memory management is in one place,
  !>    and that the String class cannot cause segfaults or allocation errors.
  type :: StringBase
    !> The contents of the string.
    character(:), allocatable, private :: contents_
  contains
    procedure, private :: check_
    
    generic,   public  :: assignment(=) => &
                        & assign_StringBase_character_
    procedure, private :: assign_StringBase_character_
    
    generic,   public             :: operator(==) =>                 &
                                   & equality_StringBase_character_, &
                                   & equality_character_StringBase_, &
                                   & equality_StringBase_StringBase_
    procedure, public             :: equality_StringBase_character_
    procedure, public, pass(that) :: equality_character_StringBase_
    procedure, public             :: equality_StringBase_StringBase_
    
    generic,   public             :: operator(/=) =>                     &
                                   & non_equality_StringBase_character_, &
                                   & non_equality_character_StringBase_, &
                                   & non_equality_StringBase_StringBase_
    procedure, public             :: non_equality_StringBase_character_
    procedure, public, pass(that) :: non_equality_character_StringBase_
    procedure, public             :: non_equality_StringBase_StringBase_
  end type
  
  interface
    !> Check that the string has been allocated. Abort otherwise.
    module subroutine check_(this)
      implicit none
      
      class(StringBase), intent(in) :: this
    end subroutine

    !> Create a [[StringBase(type)]] from a `character(*)`.
    !> Effectively a setter.
    module subroutine assign_StringBase_character_(output,input)
      implicit none
      
      class(StringBase), intent(out) :: output
      character(*),      intent(in)  :: input
    end subroutine

    !> Equality comparison between [[StringBase(type)]] and `character(*)`.
    module function equality_StringBase_character_(this,that) result(output)
      implicit none
      
      class(StringBase), intent(in) :: this
      character(*),      intent(in) :: that
      logical                       :: output
    end function
    
    !> Equality comparison between `character(*)` and [[StringBase(type)]].
    module function equality_character_StringBase_(this,that) result(output)
      implicit none
      
      character(*),      intent(in) :: this
      class(StringBase), intent(in) :: that
      logical                       :: output
    end function
    
    !> Equality comparision between [[StringBase(type)]] and
    !>    [[StringBase(type)]].
    impure elemental module function equality_StringBase_StringBase_(this, &
       & that) result(output)
      implicit none
      
      class(StringBase), intent(in) :: this
      class(StringBase), intent(in) :: that
      logical                       :: output
    end function

    !> Non-equality comparison between [[StringBase(type)]] and `character(*)`.
    module function non_equality_StringBase_character_(this,that) &
       & result(output)
      implicit none
      
      class(StringBase), intent(in) :: this
      character(*),      intent(in) :: that
      logical                       :: output
    end function
    
    !> Non-equality comparison between `character(*)` and [[StringBase(type)]].
    module function non_equality_character_StringBase_(this,that) &
       & result(output)
      implicit none
      
      character(*),      intent(in) :: this
      class(StringBase), intent(in) :: that
      logical                       :: output
    end function
    
    !> Non-equality comparison between [[StringBase(type)]] and
    !>    [[StringBase(type)]].
    impure elemental module function non_equality_StringBase_StringBase_( &
       & this,that) result(output)
      implicit none
      
      class(StringBase), intent(in) :: this
      class(StringBase), intent(in) :: that
      logical                       :: output
    end function
  end interface
  
  interface char
    !> Conversion from [[StringBase(type)]] to `character(*)`.
    !> Effectively a getter.
    module function char_StringBase(this) result(output)
      implicit none
      
      class(StringBase), intent(in) :: this
      character(:), allocatable     :: output
    end function
  end interface
end module
