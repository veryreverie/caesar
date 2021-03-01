! ======================================================================
! Calculates k-point grids for supercells.
! ======================================================================
module caesar_kpoint_grid_module
  use caesar_utils_module
  
  use caesar_structure_module
  implicit none
  
  private
  
  public :: KpointGrid
  public :: calculate_kpoint_spacing
  public :: calculate_kpoint_grid
  public :: operator(==)
  public :: operator(/=)
  
  type, extends(Stringable) :: KpointGrid
    integer :: grid(3)
  contains
    procedure, public :: read  => read_KpointGrid
    procedure, public :: write => write_KpointGrid
  end type
  
  interface KpointGrid
    module function new_KpointGrid(grid) result(this) 
      integer, intent(in) :: grid(3)
      type(KpointGrid)    :: this
    end function
  end interface
  
  interface
    module function calculate_kpoint_spacing(kpoint_grid,recip_lattice) &
       & result(output) 
      type(KpointGrid), intent(in) :: kpoint_grid
      real(dp),         intent(in) :: recip_lattice(3,3)
      real(dp)                     :: output
    end function
  end interface
  
  interface
    module function calculate_kpoint_grid(kpoint_spacing,recip_lattice) &
       & result(output) 
      real(dp), intent(in) :: kpoint_spacing
      real(dp), intent(in) :: recip_lattice(3,3)
      type(KpointGrid)     :: output
    end function
  end interface
  
  interface operator(==)
    ! ----------------------------------------------------------------------
    ! Comparison.
    ! ----------------------------------------------------------------------
    impure elemental module function equality_KpointGrid(this,that) &
       & result(output) 
      type(KpointGrid), intent(in) :: this
      type(KpointGrid), intent(in) :: that
      logical                      :: output
    end function
  end interface
  
  interface operator(/=)
    impure elemental module function non_equality_KpointGrid(this,that) &
       & result(output) 
      type(KpointGrid), intent(in) :: this
      type(KpointGrid), intent(in) :: that
      logical                      :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_KpointGrid(this,input) 
      class(KpointGrid), intent(out) :: this
      type(String),      intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_KpointGrid(this) result(output) 
      class(KpointGrid), intent(in) :: this
      type(String)                  :: output
    end function
  end interface
  
  interface KpointGrid
    impure elemental module function new_KpointGrid_String(input) result(this) 
      type(String), intent(in) :: input
      type(KpointGrid)         :: this
    end function
  end interface
end module
