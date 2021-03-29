submodule (caesar_ifile_module) caesar_ifile_submodule
  use caesar_io_module
contains
module procedure new_IFile_character
  integer          :: file_length
  integer          :: file_unit
  character(10000) :: line
  
  integer :: i,ierr,ialloc
  
  if (.not. file_exists(filename)) then
    call print_line(ERROR//': file does not exist: '//filename)
    call err()
  endif
  
  this%filename_ = filename
  
  file_length = count_lines(filename)
  
  allocate(this%lines_(file_length),stat=ialloc); call err(ialloc)
  
  file_unit = open_read_file(filename)
  do i=1,file_length
    read(file_unit,'(a)',iostat=ierr) line
    if (ierr/=0) then
      call print_line(ERROR//': failed to read file: '//filename)
      call err()
    endif
    this%lines_(i) = trim(line)
  enddo
  close(file_unit)
end procedure

module procedure new_IFile_String
  this = IFile(char(filename))
end procedure

module procedure size_IFile
  output = size(this%lines_)
end procedure

module procedure line
  output = this%lines_(line_number)
end procedure

module procedure lines_all
  output = this%lines_
end procedure

module procedure lines_slice
  output = this%lines_(first_line_number:last_line_number)
end procedure

module procedure sections_character
  type(String), allocatable :: lines(:)
  
  lines = this%lines()
  
  output = split_into_sections(lines, separating_line)
end procedure

module procedure sections_String
  output = this%sections(char(separating_line))
end procedure

module procedure count_lines
  integer      :: file_unit
  integer      :: iostat
  
  file_unit = open_read_file(filename)
  output = 0
  iostat = 0
  do while (iostat==0)
    read(file_unit, *, iostat=iostat)
    if (iostat==0) then
      output = output+1
    elseif (iostat>0) then
      call print_line('Error counting lines of '//filename)
      call err()
    endif
  enddo
  close(file_unit)
end procedure
end submodule
