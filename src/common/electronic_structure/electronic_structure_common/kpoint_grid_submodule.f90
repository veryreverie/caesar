submodule (caesar_kpoint_grid_module) caesar_kpoint_grid_submodule
  use caesar_electronic_structure_common_module
contains

module procedure new_KpointGrid
  this%grid = grid
end procedure

module procedure calculate_kpoint_spacing
  type(RealVector) :: a
  type(RealVector) :: b
  type(RealVector) :: c
  
  a = vec(recip_lattice(1,:))
  b = vec(recip_lattice(2,:))
  c = vec(recip_lattice(3,:))
  
  output = minval([ l2_norm(a)/kpoint_grid%grid(1), &
                  & l2_norm(b)/kpoint_grid%grid(2), &
                  & l2_norm(c)/kpoint_grid%grid(3)  ])
end procedure

module procedure calculate_kpoint_grid
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
end procedure

module procedure equality_KpointGrid
  output = all(this%grid==that%grid)
end procedure

module procedure non_equality_KpointGrid
  output = .not. this==that
end procedure

module procedure read_KpointGrid
  select type(this); type is(KpointGrid)
    this = KpointGrid(int(split_line(input)))
  end select
end procedure

module procedure write_KpointGrid
  select type(this); type is(KpointGrid)
    output = join(str(this%grid))
  end select
end procedure

module procedure new_KpointGrid_String
  call this%read(input)
end procedure
end submodule
