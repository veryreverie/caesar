!> Provides the [[KpointGrid(type)]] class, and related methods.
module caesar_kpoint_grid_module
  use caesar_utils_module
  
  use caesar_structure_module
  implicit none
  
  private
  
  public :: KpointGrid
  public :: operator(==)
  public :: operator(/=)
  
  !> A k-point grid from which a diagonal supercell is constructed.
  type, extends(Stringable) :: KpointGrid
    !> The number of k-points along each reciprocal lattice vector.
    !> Each element must be >= 1.
    integer :: grid(3)
  contains
    procedure, public :: kpoint_spacing => kpoint_spacing_KpointGrid
    
    procedure, public :: read  => read_KpointGrid
    procedure, public :: write => write_KpointGrid
  end type
  
  interface KpointGrid
    !> Constructor for objects of type [[KpointGrid(type)]].
    module function new_KpointGrid(grid) result(this) 
      integer, intent(in) :: grid(3)
      type(KpointGrid)    :: this
    end function
  end interface
  
  interface
    !> Returns the k-point spacing corresponding to a [[KpointGrid(type)]].
    !> This is the maximum k-point spacing across the three reciprocal
    !>    lattice vectors.
    impure elemental module function kpoint_spacing_KpointGrid(this, &
       & recip_lattice) result(output)
      class(KpointGrid), intent(in) :: this
      !> The primitive reciprocal lattice matrix.
      type(RealMatrix),  intent(in) :: recip_lattice
      real(dp)                      :: output
    end function
  
    !> Convert a [[String(type)]] to a [[KpointGrid(type)]].
    module subroutine read_KpointGrid(this,input)
      class(KpointGrid), intent(out) :: this
      type(String),      intent(in)  :: input
    end subroutine
  
    !> Convert a [[KpointGrid(type)]] to [[String(type)]].
    module function write_KpointGrid(this) result(output) 
      class(KpointGrid), intent(in) :: this
      type(String)                  :: output
    end function
  end interface
  
  interface KpointGrid
    !> Convert a [[String(type)]] to a [[KpointGrid(type)]].
    impure elemental module function new_KpointGrid_String(input) result(this) 
      type(String), intent(in) :: input
      type(KpointGrid)         :: this
    end function
  end interface
  
  interface KpointGrid
    !> Construct the smallest k-point grid whose spacing is at least as large
    !>    as `kpoint_spacing`.
    !> N.B. calling `KpointGrid(kpoint_grid%spacing())` is likely to give
    !>    a larger k-point grid than the input grid, due to finite precision.
    impure elemental module function new_KpointGrid_spacing(kpoint_spacing, &
       & recip_lattice) result(this)
      real(dp),         intent(in) :: kpoint_spacing
      !> The primitive reciprocal lattice matrix.
      type(RealMatrix), intent(in) :: recip_lattice
      type(KpointGrid)             :: this
    end function
  end interface
  
  interface operator(==)
    !> Equality comparison between two [[KpointGrid(type)]]s.
    impure elemental module function equality_KpointGrid(this,that) &
       & result(output) 
      type(KpointGrid), intent(in) :: this
      type(KpointGrid), intent(in) :: that
      logical                      :: output
    end function
  end interface
  
  interface operator(/=)
    !> Non-equality comparison between two [[KpointGrid(type)]]s.
    impure elemental module function non_equality_KpointGrid(this,that) &
       & result(output) 
      type(KpointGrid), intent(in) :: this
      type(KpointGrid), intent(in) :: that
      logical                      :: output
    end function
  end interface
end module
