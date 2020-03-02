! ======================================================================
! Checks the memory use of the running process.
! ======================================================================
module memory_checker_module
  use precision_module
  implicit none
  
  private
  
  public :: MemoryChecker
  
  type :: MemoryChecker
    integer, private :: pid_
  contains
    procedure, public :: memory => memory_MemoryChecker
  end type
  
  interface MemoryChecker
    module procedure new_MemoryChecker
  end interface
  
  ! C getpid interface.
  interface
    function getpid_c() bind(c) result(output)
      use, intrinsic :: iso_c_binding
      implicit none
      
      integer(kind=c_int) :: output
    end function
  end interface
contains

impure elemental function new_MemoryChecker() result(this)
  implicit none
  
  type(MemoryChecker) :: this
  
  this%pid_ = getpid_c()
end function

! Returns the current memory use in kB.
! N.B. does not use the machinery of IFile so that profiling does not depend
!    on io, and can be used in io_basic.
impure elemental function memory_MemoryChecker(this) result(output)
  implicit none
  
  class(MemoryChecker), intent(in) :: this
  integer                          :: output
  
  logical         :: unit_found
  logical         :: opened
  integer         :: iostat
  integer         :: file_unit
  character(1024) :: filename
  character(1024) :: line
  logical         :: line_found
  
  ! Find an available file unit.
  file_unit = 0
  unit_found = .false.
  do while(.not. unit_found)
    file_unit = file_unit+1
    
    if (any(file_unit==[0,5,6,100,101,102])) then
      cycle
    endif
    
    inquire(unit=file_unit, opened=opened, iostat=iostat)
    if (iostat==0 .and. .not. opened) then
      unit_found = .true.
    endif
  enddo
  
  ! Construct the filename /proc/pid/status from this process' pid.
  write(filename, "(A6,I0,A7)") "/proc/", this%pid_, "/status"
  
  ! Open /proc/pid/status for reading.
  open( unit=file_unit,      &
      & file=filename,       &
      & status='old',        &
      & action='read',       &
      & access='sequential', &
      & iostat=iostat        )
  call mem_err(iostat, filename, 0)
  
  ! Read
  line_found = .false.
  do while(.not. line_found)
    read(file_unit, '(a)', iostat=iostat) line
    call mem_err(iostat, filename, 1)
    if (len(line)>=7) then
      if (line(:7)=='VmSize:') then
        read(line(8:len(line)-2),*,iostat=iostat) output
        call mem_err(iostat, filename, 2)
        line_found = .true.
      endif
    endif
  enddo
  
  ! Close /proc/pid/status.
  close(file_unit)
end function

subroutine mem_err(iostat,filename,flag)
  implicit none
  
  integer,         intent(in) :: iostat
  character(1024), intent(in) :: filename
  integer,         intent(in) :: flag
  
  if (iostat/=0) then
    write(*, '(a,i0)') 'Error reading '//trim(filename)//' for profiling memory.',flag
    stop 1
  endif
end subroutine
end module
