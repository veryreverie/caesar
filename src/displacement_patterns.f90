! ======================================================================
! Parses disp_patterns.dat
! ======================================================================

! disp_patterns.dat is structured as:

! Frequency : ~
!  ~ ~ ~
! Displacement pattern for each atom:
!  ~ ~
!  ...
!  ~ ~
! [blank line]

! This is repeated no_gvectors*no_modes times

module displacement_patterns_module
  use constants, only : dp
  implicit none
  
  type DispPatterns
    real(dp), allocatable :: frequencies(:,:)
    real(dp), allocatable :: disp_patterns(:,:,:,:)
  end type
  
  interface new
    module procedure new_DispPatterns
  end interface
  
  interface drop
    module procedure drop_DispPatterns
  end interface
  
  interface read_disp_patterns_file
    module procedure read_disp_patterns_file_character
    module procedure read_disp_patterns_file_String
  end interface

contains

subroutine new_DispPatterns(this,no_modes,no_gvectors,no_atoms)
  implicit none
  
  type(DispPatterns), intent(out) :: this
  integer,            intent(in)  :: no_modes
  integer,            intent(in)  :: no_gvectors
  integer,            intent(in)  :: no_atoms
  
  allocate(this%frequencies(no_modes,no_gvectors))
  allocate(this%disp_patterns(6,no_atoms,no_modes,no_gvectors))
end subroutine

subroutine drop_DispPatterns(this)
  implicit none
  
  type(DispPatterns), intent(inout) :: this
  
  deallocate(this%frequencies)
  deallocate(this%disp_patterns)
end subroutine

function read_disp_patterns_file_character(filename,no_modes) result(this)
  use constants, only : eV
  use utils,     only : lower_case
  use file_module
  implicit none
  
  character(*), intent(in) :: filename
  integer,      intent(in) :: no_modes
  type(DispPatterns)       :: this
  
  ! Sizes
  integer        :: file_length
  integer        :: lines_per_mode
  integer        :: no_gvectors
  integer        :: no_atoms
  
  ! File unit
  integer        :: disp_patterns_file
  
  ! Temporary variables
  integer        :: i, j, k
  character(100) :: line
  character(100) :: dump
  
  file_length = count_lines(filename)
  disp_patterns_file = open_read_file(filename)
  do i=1,file_length
    read(disp_patterns_file,"(a)") line
    line = lower_case(trim(line))
    if (i/=1 .and. line(1:9)=="frequency") then
      lines_per_mode = i-1
      no_atoms = lines_per_mode-4
      no_gvectors = file_length/(no_modes*lines_per_mode)
      exit
    endif
  enddo
  
  call new(this,no_atoms,no_modes,no_gvectors)
  rewind(disp_patterns_file)
  
  do i=1,no_gvectors
    do j=1,no_modes
      read(disp_patterns_file,*) dump,dump,this%frequencies(j,i)
      read(disp_patterns_file,*)
      read(disp_patterns_file,*)
      do k=1,no_atoms
        read(disp_patterns_file,*) this%disp_patterns(:,k,j,i)
      enddo
      read(disp_patterns_file,*)
    enddo
  enddo
  
  close(disp_patterns_file)
  
  this%frequencies = this%frequencies*eV
end function

function read_disp_patterns_file_String(filename,no_modes) result(this)
  use string_module
  implicit none
  
  type(String), intent(in) :: filename
  integer,      intent(in) :: no_modes
  type(DispPatterns)       :: this
  
  this = read_disp_patterns_file(char(filename),no_modes)
end function
end module
