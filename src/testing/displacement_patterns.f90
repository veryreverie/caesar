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

module caesar_displacement_patterns_module
  use caesar_common_module
  implicit none
  
  type DispPatterns
    real(dp), allocatable :: frequencies(:,:)
    real(dp), allocatable :: disp_patterns(:,:,:,:)
    real(dp), allocatable :: prefactors(:,:,:)
  end type
  
  interface DispPatterns
    module procedure new_DispPatterns
  end interface
  
  interface read_disp_patterns_file
    module procedure read_disp_patterns_file_character
    module procedure read_disp_patterns_file_String
  end interface

contains

function new_DispPatterns(no_gvectors,no_modes,no_atoms) result(this)
  implicit none
  
  integer, intent(in) :: no_modes
  integer, intent(in) :: no_gvectors
  integer, intent(in) :: no_atoms
  type(DispPatterns)  :: this
  
  integer :: ialloc
  
  allocate( this%frequencies(no_modes,no_gvectors),              &
          & this%disp_patterns(3,no_atoms,no_modes,no_gvectors), &
          & this%prefactors(no_atoms,no_modes,no_gvectors),      &
          & stat=ialloc); call err(ialloc)
end function

function read_disp_patterns_file_character(filename,no_modes) result(this)
  implicit none
  
  character(*), intent(in) :: filename
  integer,      intent(in) :: no_modes
  type(DispPatterns)       :: this
  
  ! Sizes
  integer        :: lines_per_mode
  integer        :: no_gvectors
  integer        :: no_atoms
  
  ! File contents
  type(IFile)               :: disp_patterns_file
  type(String), allocatable :: line(:)
  
  ! Temporary variables
  integer :: i,j,k
  integer :: line_no
  
  disp_patterns_file = IFile(filename)
  
  ! Find no_atoms and no_gvectors
  do i=1,size(disp_patterns_file)
    line = split_line(lower_case(disp_patterns_file%line(i)))
    if (size(line)>=1) then
      if (i/=1 .and. line(1)=="frequency") then
        lines_per_mode = i-1
        no_atoms = lines_per_mode-4
        no_gvectors = size(disp_patterns_file)/(no_modes*lines_per_mode)
        exit
      endif
    endif
  enddo
  
  this = DispPatterns(no_gvectors,no_modes,no_atoms)
  
  do i=1,no_gvectors
    do j=1,no_modes
      line_no = ( (i-1)*no_modes + (j-1) )*(no_atoms+4) + 1
      line = split_line(disp_patterns_file%line(line_no))
      this%frequencies(j,i) = dble(line(3))
      do k=1,no_atoms
        line = split_line(disp_patterns_file%line(line_no + 2 + k))
        this%disp_patterns(:,k,j,i) = dble(line(1:3))
        this%prefactors(k,j,i) = dble(line(4))
      enddo
    enddo
  enddo
  
  this%frequencies = this%frequencies
end function

function read_disp_patterns_file_String(filename,no_modes) result(this)
  implicit none
  
  type(String), intent(in) :: filename
  integer,      intent(in) :: no_modes
  type(DispPatterns)       :: this
  
  this = read_disp_patterns_file(char(filename),no_modes)
end function
end module
