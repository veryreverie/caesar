! ======================================================================
! Provides keywords common to all modes.
! ======================================================================
module caesar_common_keywords_module
  use caesar_random_module
  use caesar_io_module
  
  use caesar_keyword_module
  implicit none
  
  private
  
  public :: common_keywords
  
  interface
    module function common_keywords() result(output) 
      type(KeywordData), allocatable :: output(:)
    end function
  end interface
end module
