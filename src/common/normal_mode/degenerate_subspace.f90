! ======================================================================
! Holds information about a degenerate subspace.
! ======================================================================
module caesar_degenerate_subspace_module
  use caesar_utils_module
  
  use caesar_structure_module
  
  use caesar_complex_mode_module
  use caesar_real_mode_module
  implicit none
  
  private
  
  public :: DegenerateSubspace
  public :: process_degeneracies
  public :: size
  
  type, extends(Stringsable) :: DegenerateSubspace
    ! The id of the degeneracy.
    integer, public :: id
    
    ! The harmonic frequency of the modes in the degenerate subspace.
    real(dp), public :: frequency
    
    ! The IDs of the modes in the degenerate subspace.
    integer, allocatable, public :: mode_ids(:)
    
    ! paired_ids(i) is the ID of the mode paired to mode_ids(i).
    integer, allocatable, public :: paired_ids(:)
  contains
    generic,   public  :: modes =>                               &
                        & modes_DegenerateSubspace_ComplexModes, &
                        & modes_DegenerateSubspace_RealModes
    procedure, private :: modes_DegenerateSubspace_ComplexModes
    procedure, private :: modes_DegenerateSubspace_RealModes
    
    generic,   public  :: qpoints =>                               &
                        & qpoints_DegenerateSubspace_ComplexModes, &
                        & qpoints_DegenerateSubspace_RealModes
    procedure, private :: qpoints_DegenerateSubspace_ComplexModes
    procedure, private :: qpoints_DegenerateSubspace_RealModes
    
    ! I/O.
    procedure, public :: read  => read_DegenerateSubspace
    procedure, public :: write => write_DegenerateSubspace
  end type
  
  interface DegenerateSubspace
    ! ----------------------------------------------------------------------
    ! Basic functionality: constructor and size() module function.
    ! ----------------------------------------------------------------------
    module function new_DegenerateSubspace(id,frequency,mode_ids,paired_ids) &
       & result(this) 
      integer,  intent(in)     :: id
      real(dp), intent(in)     :: frequency
      integer,  intent(in)     :: mode_ids(:)
      integer,  intent(in)     :: paired_ids(:)
      type(DegenerateSubspace) :: this
    end function
  end interface
  
  interface size
    module function size_DegenerateSubspace(input) result(output) 
      type(DegenerateSubspace), intent(in) :: input
      integer                              :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Construct the degenerate subspaces.
    ! ----------------------------------------------------------------------
    module function process_degeneracies(modes) result(output) 
      type(ComplexMode), intent(in)         :: modes(:)
      type(DegenerateSubspace), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the degenerate modes.
    ! ----------------------------------------------------------------------
    module function modes_DegenerateSubspace_ComplexModes(this,modes) &
       & result(output) 
      class(DegenerateSubspace), intent(in) :: this
      type(ComplexMode),         intent(in) :: modes(:)
      type(ComplexMode), allocatable        :: output(:)
    end function
  end interface
  
  interface
    module function modes_DegenerateSubspace_RealModes(this,modes) &
       & result(output) 
      class(DegenerateSubspace), intent(in) :: this
      type(RealMode),            intent(in) :: modes(:)
      type(RealMode), allocatable           :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the q-points corresponding to the degenerate modes.
    ! ----------------------------------------------------------------------
    module function qpoints_DegenerateSubspace_ComplexModes(this,modes, &
       & qpoints) result(output) 
      class(DegenerateSubspace), intent(in) :: this
      type(ComplexMode),         intent(in) :: modes(:)
      type(QpointData),          intent(in) :: qpoints(:)
      type(QpointData), allocatable         :: output(:)
    end function
  end interface
  
  interface
    module function qpoints_DegenerateSubspace_RealModes(this,modes,qpoints) &
       & result(output) 
      class(DegenerateSubspace), intent(in) :: this
      type(RealMode),            intent(in) :: modes(:)
      type(QpointData),          intent(in) :: qpoints(:)
      type(QpointData), allocatable         :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_DegenerateSubspace(this,input) 
      class(DegenerateSubspace), intent(out) :: this
      type(String),              intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_DegenerateSubspace(this) result(output) 
      class(DegenerateSubspace), intent(in) :: this
      type(String), allocatable             :: output(:)
    end function
  end interface
  
  interface DegenerateSubspace
    module function new_DegenerateSubspace_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(DegenerateSubspace) :: this
    end function
  
    impure elemental module function new_DegenerateSubspace_StringArray(input) result(this) 
      type(StringArray), intent(in) :: input
      type(DegenerateSubspace)      :: this
    end function
  end interface
end module
