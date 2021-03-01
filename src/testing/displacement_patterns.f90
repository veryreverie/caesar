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

  interface
    module function new_DispPatterns(no_gvectors,no_modes,no_atoms) &
       & result(this) 
      integer, intent(in) :: no_modes
      integer, intent(in) :: no_gvectors
      integer, intent(in) :: no_atoms
      type(DispPatterns)  :: this
    end function
  end interface
  
  interface
    module function read_disp_patterns_file_character(filename,no_modes) &
       & result(this) 
      character(*), intent(in) :: filename
      integer,      intent(in) :: no_modes
      type(DispPatterns)       :: this
    end function
  end interface
  
  interface
    module function read_disp_patterns_file_String(filename,no_modes) &
       & result(this) 
      type(String), intent(in) :: filename
      integer,      intent(in) :: no_modes
      type(DispPatterns)       :: this
    end function
  end interface
end module
