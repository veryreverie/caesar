submodule (caesar_mass_weighted_force_module) caesar_mass_weighted_force_submodule
  use caesar_normal_mode_module
contains

module procedure new_MassWeightedForce
  this%vectors = forces
end procedure

module procedure size_MassWeightedForce
  output = size(this%vectors)
end procedure

module procedure new_MassWeightedForce_zero
  integer :: i
  
  this%vectors = [(dblevec(zeroes(3)), i=1, structure%no_atoms)]
end procedure

module procedure new_MassWeightedForce_CartesianForce
  output = MassWeightedForce(input%vectors/sqrt(structure%atoms%mass()))
end procedure

module procedure new_CartesianForce_MassWeightedForce
  output = CartesianForce(input%vectors*sqrt(structure%atoms%mass()))
end procedure

module procedure multiply_real_MassWeightedForce
  output = MassWeightedForce(this * that%vectors)
end procedure

module procedure multiply_MassWeightedForce_real
  output = MassWeightedForce(this%vectors * that)
end procedure

module procedure divide_MassWeightedForce_real
  output = MassWeightedForce(this%vectors / that)
end procedure

module procedure add_MassWeightedForce_MassWeightedForce
  output = MassWeightedForce(this%vectors + that%vectors)
end procedure

module procedure sum_MassWeightedForces
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

module procedure negative_MassWeightedForce
  output = MassWeightedForce(-this%vectors)
end procedure

module procedure subtract_MassWeightedForce_MassWeightedForce
  output = MassWeightedForce(this%vectors - that%vectors)
end procedure

module procedure read_MassWeightedForce
  select type(this); type is(MassWeightedForce)
    this = MassWeightedForce(RealVector(input))
  class default
    call err()
  end select
end procedure

module procedure write_MassWeightedForce
  select type(this); type is(MassWeightedForce)
    output = str(this%vectors)
  class default
    call err()
  end select
end procedure

module procedure new_MassWeightedForce_Strings
  call this%read(input)
end procedure

module procedure new_MassWeightedForce_StringArray
  this = MassWeightedForce(str(input))
end procedure
end submodule
