! ======================================================================
! Generate the frequency of maximum displacement from
!    the maximum energy of displacement,
!    or vice versa.
! ======================================================================
module caesar_max_displacement_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: MaxDisplacement
  
  type, extends(Stringsable) :: MaxDisplacement
    real(dp) :: maximum_weighted_displacement
    real(dp) :: frequency_of_max_displacement
    real(dp) :: max_energy_of_displacement
  contains
    ! I/O.
    procedure, public :: read  => read_MaxDisplacement
    procedure, public :: write => write_MaxDisplacement
  end type
  
  interface MaxDisplacement
    impure elemental module function new_MaxDisplacement(maximum_weighted_displacement,frequency_of_max_displacement,max_energy_of_displacement) result(this) 
      real(dp), intent(in), optional :: maximum_weighted_displacement
      real(dp), intent(in), optional :: frequency_of_max_displacement
      real(dp), intent(in), optional :: max_energy_of_displacement
      type(MaxDisplacement)          :: this
    end function
  
    impure elemental module function new_MaxDisplacement_displacement(          maximum_displacement,structure,frequency_of_max_displacement,max_energy_of_displacement) result(this) 
      real(dp),            intent(in)           :: maximum_displacement
      type(StructureData), intent(in)           :: structure
      real(dp),            intent(in), optional :: frequency_of_max_displacement
      real(dp),            intent(in), optional :: max_energy_of_displacement
      type(MaxDisplacement)                     :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_MaxDisplacement(this,input) 
      class(MaxDisplacement), intent(out) :: this
      type(String),           intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_MaxDisplacement(this) result(output) 
      class(MaxDisplacement), intent(in) :: this
      type(String), allocatable          :: output(:)
    end function
  end interface
  
  interface MaxDisplacement
    module function new_MaxDisplacement_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(MaxDisplacement)    :: this
    end function
  
    impure elemental module function new_MaxDisplacement_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(MaxDisplacement)         :: this
    end function
  end interface
end module
