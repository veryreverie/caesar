submodule (caesar_ofile_target_module) caesar_ofile_target_submodule
  use caesar_io_module
contains
module procedure new_OFileTarget
  this%open_      = .true.
  this%filename_  = filename
  this%file_unit_ = open_write_file(filename)
  this%is_stdout_ = .false.
end procedure

module procedure close_OFileTarget
  integer :: ierr
  
  ! Check file is open.
  if (.not. this%is_open()) then
    call print_line(ERROR//': Trying to close file which is not open.')
    call err()
  endif
  
  ! Close file
  close(this%file_unit_,iostat=ierr)
  if (ierr/=0) then
    call print_line(ERROR//': could not close file. Error code '//ierr)
    call err()
  endif
  
  ! Reset stdout if necessary.
  if (this%is_stdout_) then
    call unset_output_unit()
    this%is_stdout_ = .false.
  endif
end procedure

module procedure is_open
  output = this%open_
end procedure

module procedure make_stdout
  if (.not. this%is_open()) then
    call print_line(CODE_ERROR//': attempted to point stdout to a file which &
       &has either not been opened or has already been closed.')
    call err()
  endif
  
  call set_output_unit(this%file_unit_)
  this%is_stdout_ = .true.
end procedure

module procedure print_line_character
  integer :: ierr
  
  if (.not. this%is_open()) then
    call print_line(CODE_ERROR//': attempted to write to a file which has &
       &either not been opened or has already been closed.')
    call err()
  endif
  
  write(this%file_unit_,'(a)',iostat=ierr) input
  
  if (ierr/=0) then
    call print_line(ERROR//': writing to file failed. Error code '//ierr)
    call err()
  endif
  
  flush(this%file_unit_,iostat=ierr)
  
  if (ierr/=0) then
    call print_line(ERROR//': flushing file failed. Error code '//ierr)
    call err()
  endif
end procedure
end submodule
