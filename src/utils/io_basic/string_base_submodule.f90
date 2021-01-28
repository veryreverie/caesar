submodule (caesar_string_base_module) caesar_string_base_submodule
contains
module procedure check_
  if (.not. allocated(this%contents_)) then
    write(*,*) CODE_ERROR//': Trying to use the contents of a string before &
       &it has been allocated.'
    call err()
  endif
end procedure

module procedure assign_StringBase_character_
  ! Currently uses naive memory allocation, reallocating all memory for every
  !    operation.
  ! If it is desired that String be sped up, this should probably be changed to
  !    a smarter memory scheme, e.g. C++'s std::vector model.
  output%contents_ = input
end procedure

module procedure equality_StringBase_character_
  call this%check_()
  output = this%contents_==that
end procedure

module procedure equality_character_StringBase_
  call that%check_()
  output = this==that%contents_
end procedure

module procedure equality_StringBase_StringBase_
  call this%check_()
  call that%check_()
  output = this%contents_==that%contents_
end procedure

module procedure non_equality_StringBase_character_
  output = .not. this==that
end procedure

module procedure non_equality_character_StringBase_
  output = .not. this==that
end procedure

module procedure non_equality_StringBase_StringBase_
  output = .not. this==that
end procedure

module procedure char_StringBase
  call this%check_()
  output = this%contents_
end procedure
end submodule
