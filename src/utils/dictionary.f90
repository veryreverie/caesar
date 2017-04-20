! ======================================================================
! A dictionary of keys and values, both of type String.
! ======================================================================
module dictionary_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  private
  
  public :: Dictionary      ! The Dictionary type.
  public :: size            ! The number of key/value pairs.
  public :: make_dictionary ! Make from (array of keys, array of values).
  public :: operator(//)    ! Concatenate two Dictrionaries.
  public :: item            ! Get the value corresponding to a given key.
  public :: write_dictionary_file ! Write Dictionary to file.
  
  type Dictionary
    type(String), allocatable :: keys(:)
    type(String), allocatable :: values(:)
  end type
  
  interface new
    module procedure new_Dictionary
  end interface
  
  interface size
    module procedure size_Dictionary
  end interface
  
  interface operator(//)
    module procedure concatenate_Dictionary_Dictionary
  end interface
  
  interface item
    module procedure item_character
    module procedure item_String
  end interface
  
contains

! ----------------------------------------------------------------------
! Private allocate(Dictionary) subroutine.
! ----------------------------------------------------------------------
subroutine new_Dictionary(this,no_entries)
  implicit none
  
  type(Dictionary), intent(out) :: this
  integer,          intent(in)  :: no_entries
  
  integer :: ialloc
  
  allocate( this%keys(no_entries),   &
          & this%values(no_entries), &
          & stat=ialloc); call err(ialloc)
end subroutine

! ----------------------------------------------------------------------
! size(Dictionary).
! ----------------------------------------------------------------------
function size_Dictionary(this) result(output)
  implicit none
  
  type(Dictionary), intent(in) :: this
  integer                      :: output
  
  output = size(this%keys)
end function

! ----------------------------------------------------------------------
! Takes an array of keys and an array of values and returns a Dictionary.
! ----------------------------------------------------------------------
function make_dictionary(keys,values) result(output)
  implicit none
  
  type(String), intent(in) :: keys(:)
  type(String), intent(in) :: values(:)
  type(Dictionary)         :: output
  
  if (size(keys) /= size(values)) then
    call print_line('Error: cannor make a dictionary with different numbers &
       &of keys and values.')
    call err()
  endif
  
  call new(output,size(keys))
  output%keys = keys
  output%values = values
end function

! ----------------------------------------------------------------------
! Concatenates two Dictionaries.
! ----------------------------------------------------------------------
function concatenate_Dictionary_Dictionary(dict1,dict2) result(output)
  implicit none
  
  type(Dictionary), intent(in) :: dict1
  type(Dictionary), intent(in) :: dict2
  type(Dictionary)             :: output
  
  call new(output,size(dict1)+size(dict2))
  output%keys(:size(dict1)) = dict1%keys
  output%values(:size(dict1)) = dict1%values
  output%keys(size(dict1)+1:) = dict2%keys
  output%values(size(dict1)+1:) = dict2%values
end function

! ----------------------------------------------------------------------
! Get the value corresponding to a given key.
! ----------------------------------------------------------------------
! Returns '' if the key is not found.
! If there are duplicate keys, returns the first match.
function item_character(this,key) result(output)
  implicit none
  
  type(Dictionary), intent(in) :: this
  character(*),     intent(in) :: key
  type(String)                 :: output
  
  integer :: i
  
  output = ''
  do i=1,size(this)
    if (this%keys(i) == key) then
      output = this%values(i)
      exit
    endif
  enddo
end function

function item_String(this,key) result(output)
  implicit none
  
  type(Dictionary), intent(in) :: this
  type(String),     intent(in) :: key
  type(String)                 :: output
  
  output = item(this,char(key))
end function

! ----------------------------------------------------------------------
! Writes a Dictionary to file.
! ----------------------------------------------------------------------
subroutine write_dictionary_file(this,filename)
  implicit none
  
  type(Dictionary), intent(in) :: this
  type(String),     intent(in) :: filename
  
  integer :: i,j
  integer :: dictionary_file
  
  dictionary_file = open_write_file(filename)
  do_i : do i=1,size(this)
    
    ! Ignore duplicate entries.
    do j=1,i-1
      if (this%keys(i) == this%keys(j)) then
        cycle do_i
      endif
    enddo
    
    call print_line(dictionary_file,this%keys(i)//' '//this%values(i))
  enddo do_i
  close(dictionary_file)
end subroutine

end module
