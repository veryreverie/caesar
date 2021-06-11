!> Provides the [[MaxDisplacement(type)]] class, and related methods.
module caesar_max_displacement_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: MaxDisplacement
  
  !> Stores details of the maximum displacement from the static-lattice
  !>    configuration at which the potential will be sampled.
  !> Modes with frequency `w < frequency_of_max_displacement` are
  !>    displaced up to `u < maximum_weighted_displacement`.
  !> Modes with frequency `w > frequency_of_max_displacement` are
  !>    displaced up to `1/2 w^2 u^2 < max_energy_of_displacement`.
  type, extends(Stringsable) :: MaxDisplacement
    real(dp) :: maximum_weighted_displacement
    real(dp) :: frequency_of_max_displacement
    real(dp) :: max_energy_of_displacement
  contains
    procedure, public :: max_displacement
    
    procedure, public :: read  => read_MaxDisplacement
    procedure, public :: write => write_MaxDisplacement
  end type
  
  interface MaxDisplacement
    !> Construct a [[MaxDisplacement(type)]] from at least two of
    !>   `maximum_weighted_displacement`, `frequency_of_max_displacement`,
    !>    and `max_energy_of_displacement`.
    impure elemental module function new_MaxDisplacement( &
       & maximum_weighted_displacement,frequency_of_max_displacement, &
       & max_energy_of_displacement) result(this) 
      real(dp), intent(in), optional :: maximum_weighted_displacement
      real(dp), intent(in), optional :: frequency_of_max_displacement
      real(dp), intent(in), optional :: max_energy_of_displacement
      type(MaxDisplacement)          :: this
    end function
  
    !> Construct a [[MaxDisplacement(type)]] from the un-mass-weighted
    !>    `maximum_displacement`, a `structure` (which may be a supercell),
    !>    and at least one of `frequency_of_max_displacement` and
    !>    `max_energy_of_displacement`.
    impure elemental module function new_MaxDisplacement_displacement( &
       & maximum_displacement,structure,frequency_of_max_displacement, &
       & max_energy_of_displacement) result(this) 
      real(dp),            intent(in)           :: maximum_displacement
      type(StructureData), intent(in)           :: structure
      real(dp),            intent(in), optional :: frequency_of_max_displacement
      real(dp),            intent(in), optional :: max_energy_of_displacement
      type(MaxDisplacement)                     :: this
    end function
  end interface
  
  interface
    !> Returns the maximum displacement along a mode, given the `frequency`
    !>    of that mode.
    impure elemental module function max_displacement(this,frequency) &
       & result(output)
      class(MaxDisplacement), intent(in) :: this
      real(dp),               intent(in) :: frequency
      real(dp)                           :: output
    end function
    
    !> Convert a [[String(type)]] array to [[MaxDisplacement(type)]].
    module subroutine read_MaxDisplacement(this,input) 
      class(MaxDisplacement), intent(out) :: this
      type(String),           intent(in)  :: input(:)
    end subroutine
  
    !> Convert a [[MaxDisplacement(type)]] to a [[String(type)]] array.
    module function write_MaxDisplacement(this) result(output) 
      class(MaxDisplacement), intent(in) :: this
      type(String), allocatable          :: output(:)
    end function
  end interface
  
  interface MaxDisplacement
    !> Convert a [[String(type)]] array to [[MaxDisplacement(type)]].
    module function new_MaxDisplacement_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(MaxDisplacement)    :: this
    end function
  
    !> Convert a [[StringArray(type)]] to [[MaxDisplacement(type)]].
    impure elemental module function new_MaxDisplacement_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(MaxDisplacement)         :: this
    end function
  end interface
end module
