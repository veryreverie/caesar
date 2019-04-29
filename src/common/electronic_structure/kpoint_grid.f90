! ======================================================================
! Calculates k-point grids for supercells.
! ======================================================================
module kpoint_grid_module
  use utils_module
  
  use structure_module
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
    module procedure new_KpointGrid
    module procedure new_KpointGrid_String
  end interface
  
  interface operator(==)
    module procedure equality_KpointGrid
  end interface
  
  interface operator(/=)
    module procedure non_equality_KpointGrid
  end interface
contains

function new_KpointGrid(grid) result(this)
  implicit none
  
  integer, intent(in) :: grid(3)
  type(KpointGrid)    :: this
  
  this%grid = grid
end function

function calculate_kpoint_spacing(kpoint_grid,recip_lattice) result(output)
  implicit none
  
  type(KpointGrid), intent(in) :: kpoint_grid
  real(dp),         intent(in) :: recip_lattice(3,3)
  real(dp)                     :: output
  
  type(RealVector) :: a
  type(RealVector) :: b
  type(RealVector) :: c
  
  a = vec(recip_lattice(1,:))
  b = vec(recip_lattice(2,:))
  c = vec(recip_lattice(3,:))
  
  output = minval([ l2_norm(a)/kpoint_grid%grid(1), &
                  & l2_norm(b)/kpoint_grid%grid(2), &
                  & l2_norm(c)/kpoint_grid%grid(3)  ])
end function

function calculate_kpoint_grid(kpoint_spacing,recip_lattice) result(output)
  implicit none
  
  real(dp), intent(in) :: kpoint_spacing
  real(dp), intent(in) :: recip_lattice(3,3)
  type(KpointGrid)     :: output
  
  type(RealVector) :: a
  type(RealVector) :: b
  type(RealVector) :: c
  
  a = vec(recip_lattice(1,:))
  b = vec(recip_lattice(2,:))
  c = vec(recip_lattice(3,:))
  
  ! N.B. the -0.1 is to ensure that calculating a grid from a spacing which
  !    was in turn calculated from a grid will give the same answer.
  ! The max(x,1) is to ensure the above correction does not lead to a k-point
  !    grid with 0 k-points.
  output = KpointGrid(max(                                  &
     & ceiling( [ l2_norm(a)/kpoint_spacing,                &
     &            l2_norm(b)/kpoint_spacing,                &
     &            l2_norm(c)/kpoint_spacing  ] - 0.1_dp  ), &
     & 1                                                    ))
end function

! ----------------------------------------------------------------------
! Comparison.
! ----------------------------------------------------------------------
impure elemental function equality_KpointGrid(this,that) result(output)
  implicit none
  
  type(KpointGrid), intent(in) :: this
  type(KpointGrid), intent(in) :: that
  logical                      :: output
  
  output = all(this%grid==that%grid)
end function

impure elemental function non_equality_KpointGrid(this,that) result(output)
  implicit none
  
  type(KpointGrid), intent(in) :: this
  type(KpointGrid), intent(in) :: that
  logical                      :: output
  
  output = .not. this==that
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_KpointGrid(this,input)
  implicit none
  
  class(KpointGrid), intent(out) :: this
  type(String),      intent(in)  :: input
  
  select type(this); type is(KpointGrid)
    this = KpointGrid(int(split_line(input)))
  end select
end subroutine

function write_KpointGrid(this) result(output)
  implicit none
  
  class(KpointGrid), intent(in) :: this
  type(String)                  :: output
  
  select type(this); type is(KpointGrid)
    output = join(str(this%grid))
  end select
end function

impure elemental function new_KpointGrid_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(KpointGrid)         :: this
  
  call this%read(input)
end function
end module
