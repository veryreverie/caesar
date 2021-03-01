! ======================================================================
! Determines which directions are related by symmetry,
!    and thus which atoms should be displaced in which directions 
!    in order to construct the Hessian.
! ======================================================================
! The symmetry operations to consider are outlined here:
! http://www.homepages.ucl.ac.uk/~ucfbdxa/phon/node4.html
module caesar_unique_directions_module
  use caesar_utils_module
  
  use caesar_structure_module
  use caesar_normal_mode_module
  implicit none
  
  private
  
  public :: UniqueDirection
  public :: calculate_unique_directions
  public :: CartesianDisplacement
  
  ! A direction in which an atom needs to be perturbed in order to map out
  !    the harmonic Born-Oppenheimer surface.
  type, extends(Stringsable) :: UniqueDirection
    ! Atom id.
    integer :: atom_id
    
    ! The direction in which the atom is displaced: '+dx', '-dx', ..., '-dz'.
    type(String) :: direction
    
    ! The vector through which the atom is displaced.
    type(RealVector) :: atomic_displacement
  contains
    ! I/O.
    procedure, public :: read  => read_UniqueDirection
    procedure, public :: write => write_UniqueDirection
  end type
  
  interface UniqueDirection
    ! ----------------------------------------------------------------------
    ! Constructor
    ! ----------------------------------------------------------------------
    module function new_UniqueDirection(atom_id,direction, &
       & atomic_displacement) result(output) 
      integer,          intent(in) :: atom_id
      type(String),     intent(in) :: direction
      type(Realvector), intent(in) :: atomic_displacement
      type(UniqueDirection)        :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Works out which atoms need to be perturbed in which directions
    !    in order to map out the harmonic Born-Oppenheimer surface.
    ! ----------------------------------------------------------------------
    module function calculate_unique_directions(structure, &
       & harmonic_displacement) result(output) 
      type(StructureData), intent(in)    :: structure
      real(dp),            intent(in)    :: harmonic_displacement
      type(UniqueDirection), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Check that symmetry can be used to reconstruct all directions from
    !    the unique directions.
    ! ----------------------------------------------------------------------
    module subroutine check_unique_directions(unique_directions,structure,a, &
       & b,c) 
      type(UniqueDirection), intent(in) :: unique_directions(:)
      type(StructureData),   intent(in) :: structure
      type(RealVector),      intent(in) :: a
      type(RealVector),      intent(in) :: b
      type(RealVector),      intent(in) :: c
    end subroutine
  end interface
  
  interface CartesianDisplacement
    ! ----------------------------------------------------------------------
    ! Construct the displacement for all atoms. This is atomic_displacement
    !    for the displaced atom, and zero for all other atoms.
    ! ----------------------------------------------------------------------
    impure elemental module function new_CartesianDisplacement_UniqueDirection(input,structure) result(this) 
      type(UniqueDirection), intent(in) :: input
      type(StructureData),   intent(in) :: structure
      type(CartesianDisplacement)       :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_UniqueDirection(this,input) 
      class(UniqueDirection), intent(out) :: this
      type(String),           intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_UniqueDirection(this) result(output) 
      class(UniqueDirection), intent(in) :: this
      type(String), allocatable          :: output(:)
    end function
  end interface
  
  interface UniqueDirection
    module function new_UniqueDirection_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(UniqueDirection)    :: this
    end function
  end interface
  
  interface UniqueDirection
    impure elemental module function new_UniqueDirection_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(UniqueDirection)         :: this
    end function
  end interface
end module
