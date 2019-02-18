! ======================================================================
! Calculates k-point grids for supercells.
! ======================================================================
module kpoint_grid_module
  use utils_module
  
  use structure_module
  implicit none
  
  private
  
  public :: calculate_kpoint_spacing
  public :: calculate_kpoint_grid
contains

function calculate_kpoint_spacing(kpoint_grid,recip_lattice) result(output)
  implicit none
  
  integer,  intent(in) :: kpoint_grid(3)
  real(dp), intent(in) :: recip_lattice(3,3)
  real(dp)             :: output
  
  type(RealVector) :: a
  type(RealVector) :: b
  type(RealVector) :: c
  
  a = vec(recip_lattice(1,:))
  b = vec(recip_lattice(2,:))
  c = vec(recip_lattice(3,:))
  
  output = minval([ l2_norm(a)/kpoint_grid(1), &
                  & l2_norm(b)/kpoint_grid(2), &
                  & l2_norm(c)/kpoint_grid(3)  ])
end function

function calculate_kpoint_grid(kpoint_spacing,recip_lattice) result(output)
  implicit none
  
  real(dp), intent(in) :: kpoint_spacing
  real(dp), intent(in) :: recip_lattice(3,3)
  integer              :: output(3)
  
  type(RealVector) :: a
  type(RealVector) :: b
  type(RealVector) :: c
  
  a = vec(recip_lattice(1,:))
  b = vec(recip_lattice(2,:))
  c = vec(recip_lattice(3,:))
  
  output = ceiling( [ l2_norm(a)/kpoint_spacing,    &
                  &   l2_norm(b)/kpoint_spacing,    &
                  &   l2_norm(c)/kpoint_spacing  ], &
                  & dp                              )
end function
end module
