submodule (caesar_cartesian_force_module) caesar_cartesian_force_submodule
  use caesar_normal_mode_module
contains

module procedure new_CartesianForce
  this%vectors = forces
end procedure

module procedure size_CartesianForce
  output = size(this%vectors)
end procedure

module procedure new_CartesianForce_zero
  integer :: i
  
  this%vectors = [(dblevec(zeroes(3)), i=1, structure%no_atoms)]
end procedure

module procedure multiply_real_CartesianForce
  output = CartesianForce(this * that%vectors)
end procedure

module procedure multiply_CartesianForce_real
  output = CartesianForce(this%vectors * that)
end procedure

module procedure divide_CartesianForce_real
  output = CartesianForce(this%vectors / that)
end procedure

module procedure add_CartesianForce_CartesianForce
  output = CartesianForce(this%vectors + that%vectors)
end procedure

module procedure sum_CartesianForces
  integer :: i
  
  if (size(input)==0) then
    call print_line(ERROR//': Trying to sum an empty array.')
    call err()
  endif
  
  output = input(1)
  do i=2,size(input)
    output = output + input(i)
  enddo
end procedure

module procedure negative_CartesianForce
  output = CartesianForce(-this%vectors)
end procedure

module procedure subtract_CartesianForce_CartesianForce
  output = CartesianForce(this%vectors - that%vectors)
end procedure

module procedure read_CartesianForce
  select type(this); type is(CartesianForce)
    this = CartesianForce(RealVector(input))
  class default
    call err()
  end select
end procedure

module procedure write_CartesianForce
  select type(this); type is(CartesianForce)
    output = str(this%vectors)
  class default
    call err()
  end select
end procedure

module procedure new_CartesianForce_Strings
  call this%read(input)
end procedure

module procedure new_CartesianForce_StringArray
  this = CartesianForce(str(input))
end procedure
end submodule
