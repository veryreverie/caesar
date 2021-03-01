submodule (caesar_cartesian_displacement_module) caesar_cartesian_displacement_submodule
  use caesar_normal_mode_module
contains

module procedure new_CartesianDisplacement
  this%vectors = displacements
end procedure

module procedure size_CartesianDisplacement
  output = size(this%vectors)
end procedure

module procedure new_CartesianDisplacement_zero
  integer :: i
  
  this%vectors = [(dblevec(zeroes(3)), i=1, structure%no_atoms)]
end procedure

module procedure displace_structure_CartesianDisplacement
  integer :: i
  
  if (structure%no_atoms/=size(displacement)) then
    call print_line(CODE_ERROR//': Trying to displace a structure by a &
       &displacement which does not match the number of atoms.')
    call err()
  endif
  
  output = structure
  
  do i=1,size(displacement)
    call output%atoms(i)%set_cartesian_position( &
        &   output%atoms(i)%cartesian_position() &
        & + displacement%vectors(i))
  enddo
end procedure

module procedure multiply_real_CartesianDisplacement
  output = CartesianDisplacement(this * that%vectors)
end procedure

module procedure multiply_CartesianDisplacement_real
  output = CartesianDisplacement(this%vectors * that)
end procedure

module procedure divide_CartesianDisplacement_real
  output = CartesianDisplacement(this%vectors / that)
end procedure

module procedure add_CartesianDisplacement_CartesianDisplacement
  output = CartesianDisplacement(this%vectors + that%vectors)
end procedure

module procedure sum_CartesianDisplacements
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

module procedure negative_CartesianDisplacement
  output = CartesianDisplacement(-this%vectors)
end procedure

module procedure subtract_CartesianDisplacement_CartesianDisplacement
  output = CartesianDisplacement(this%vectors - that%vectors)
end procedure

module procedure read_CartesianDisplacement
  select type(this); type is(CartesianDisplacement)
    this = CartesianDisplacement(RealVector(input))
  class default
    call err()
  end select
end procedure

module procedure write_CartesianDisplacement
  select type(this); type is(CartesianDisplacement)
    output = str(this%vectors)
  class default
    call err()
  end select
end procedure

module procedure new_CartesianDisplacement_Strings
  call this%read(input)
end procedure

module procedure new_CartesianDisplacement_StringArray
  this = CartesianDisplacement(str(input))
end procedure
end submodule
