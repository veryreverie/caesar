submodule (caesar_mass_weighted_displacement_module) caesar_mass_weighted_displacement_submodule
  use caesar_normal_mode_module
contains

module procedure new_MassWeightedDisplacement
  this%vectors = displacements
end procedure

module procedure size_MassWeightedDisplacement
  output = size(this%vectors)
end procedure

module procedure new_MassWeightedDisplacement_zero
  integer :: i
  
  this%vectors = [(dblevec(zeroes(3)), i=1, structure%no_atoms)]
end procedure

module procedure new_MassWeightedDisplacement_CartesianDisplacement
  output = MassWeightedDisplacement(input%vectors*sqrt(structure%atoms%mass()))
end procedure

module procedure new_CartesianDisplacement_MassWeightedDisplacement
  output = CartesianDisplacement(input%vectors/sqrt(structure%atoms%mass()))
end procedure

module procedure displace_structure_MassWeightedDisplacement
  type(CartesianDisplacement) :: cartesian_displacement
  
  cartesian_displacement = CartesianDisplacement(displacement,structure)
  output = displace_structure(structure,cartesian_displacement)
end procedure

module procedure multiply_real_MassWeightedDisplacement
  output = MassWeightedDisplacement(this * that%vectors)
end procedure

module procedure multiply_MassWeightedDisplacement_real
  output = MassWeightedDisplacement(this%vectors * that)
end procedure

module procedure divide_MassWeightedDisplacement_real
  output = MassWeightedDisplacement(this%vectors / that)
end procedure

module procedure add_MassWeightedDisplacement_MassWeightedDisplacement
  output = MassWeightedDisplacement(this%vectors + that%vectors)
end procedure

module procedure sum_MassWeightedDisplacements
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

module procedure negative_MassWeightedDisplacement
  output = MassWeightedDisplacement(-this%vectors)
end procedure

module procedure subtract_MassWeightedDisplacement_MassWeightedDisplacement
  output = MassWeightedDisplacement(this%vectors - that%vectors)
end procedure

module procedure read_MassWeightedDisplacement
  select type(this); type is(MassWeightedDisplacement)
    this = MassWeightedDisplacement(RealVector(input))
  class default
    call err()
  end select
end procedure

module procedure write_MassWeightedDisplacement
  select type(this); type is(MassWeightedDisplacement)
    output = str(this%vectors)
  class default
    call err()
  end select
end procedure

module procedure new_MassWeightedDisplacement_Strings
  call this%read(input)
end procedure

module procedure new_MassWeightedDisplacement_StringArray
  this = MassWeightedDisplacement(str(input))
end procedure
end submodule
