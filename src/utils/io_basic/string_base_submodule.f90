submodule (caesar_string_base_module) caesar_string_base_submodule
contains
module procedure check_
  if (.not. allocated(this%contents_)) then
    write(*,*) CODE_ERROR//': Trying to use the contents of a string before &
       &it has been allocated.'
    call err()
  endif
end procedure

module procedure assign_StringBase_character
  ! Currently uses naive memory allocation, reallocating all memory for every
  !    operation.
  ! If it is desired that String be sped up, this should probably be changed to
  !    a smarter memory scheme, e.g. C++'s std::vector model.
  output%contents_ = input
end procedure

module procedure equality_StringBase_character
  call this%check_()
  output = this%contents_==that
end procedure

module procedure equality_character_StringBase
  call that%check_()
  output = this==that%contents_
end procedure

module procedure equality_StringBase_StringBase
  call this%check_()
  call that%check_()
  output = this%contents_==that%contents_
end procedure

module procedure non_equality_StringBase_character
  output = .not. this==that
end procedure

module procedure non_equality_character_StringBase
  output = .not. this==that
end procedure

module procedure non_equality_StringBase_StringBase
  output = .not. this==that
end procedure

module procedure char_StringBase
  call this%check_()
  output = this%contents_
end procedure

module procedure char_StringBases
  integer :: i
  
  output = ''
  do i=1,size(this)
    output = output//char(this(i))
    if (i<size(this)) then
      output = output//new_line(' ')
    endif
  enddo
end procedure
end submodule
