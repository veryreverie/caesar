submodule (caesar_shared_counter_module) caesar_shared_counter_submodule
  use caesar_io_module
contains
module procedure new_SharedCounter
  integer :: ialloc
  
  allocate(this%no_copies_, stat=ialloc); call err(ialloc)
  if (SHARED_COUNTER_BUG) then
    this%no_copies_ = 0
  else
    this%no_copies_ = 1
  endif
end procedure

module procedure assign_SharedCounter_SharedCounter
  output%no_copies_ => input%no_copies_
  output%no_copies_ =  output%no_copies_ + 1
end procedure

module procedure final_SharedCounter
  if (associated(this%no_copies_)) then
    this%no_copies_ = this%no_copies_ - 1
  endif
end procedure

module procedure is_only_copy_SharedCounter
  if (.not. associated(this%no_copies_)) then
    output = .false.
  else
    output = this%no_copies_==1
  endif
end procedure
end submodule
