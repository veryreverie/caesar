submodule (caesar_kpoint_grid_module) caesar_kpoint_grid_submodule
  use caesar_electronic_structure_common_module
contains

module procedure new_KpointGrid
  if (any(grid<1)) then
    call print_line(ERROR//': The elements of a k-point grid must be >= 1.')
    call err()
  endif
  this%grid = grid
end procedure

module procedure kpoint_spacing_KpointGrid
  integer :: i
  
  output = maxval([( l2_norm(recip_lattice%row(i))/this%grid(i), &
                   & i=1,                                        &
                   & 3                                           )])
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

module procedure new_KpointGrid_spacing
  integer :: i
  
  this%grid = [(                                                      &
     & max( ceiling(l2_norm(recip_lattice%row(i))/kpoint_spacing),    &
     &      1                                                      ), &
     & i=1,                                                           &
     & 3                                                              )]
end procedure

module procedure equality_KpointGrid
  output = all(this%grid==that%grid)
end procedure

module procedure non_equality_KpointGrid
  output = .not. this==that
end procedure
end submodule
