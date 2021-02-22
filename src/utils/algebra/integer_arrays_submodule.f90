submodule (caesar_integer_arrays_module) caesar_integer_arrays_submodule
  use caesar_algebra_module
contains

module procedure new_IntArray1D
  this%i = i
end procedure

module procedure new_IntArray2D
  this%i = i
end procedure

module procedure array_IntArray1D_integers
  output = IntArray1D(input)
end procedure

module procedure array_IntArray2D_IntArray1Ds
  output = IntArray2D(input)
end procedure

module procedure size_IntArray1D
  output = size(input%i)
end procedure

module procedure size_IntArray2D
  output = size(input%i)
end procedure

module procedure concatenate_IntArray1D_integer
  output = array([this%i, that])
end procedure

module procedure concatenate_integer_IntArray1D
  output = array([this, that%i])
end procedure

module procedure concatenate_IntArray1D_integers
  output = array([this%i, that])
end procedure

module procedure concatenate_integers_IntArray1D
  output = array([this, that%i])
end procedure

module procedure concatenate_IntArray1D_IntArray1D
  output = array([this%i, that%i])
end procedure

module procedure concatenate_IntArray2D_IntArray1D
  output = array([this%i, that])
end procedure

module procedure concatenate_IntArray1D_IntArray2D
  output = array([this, that%i])
end procedure

module procedure concatenate_IntArray2D_IntArray1Ds
  output = array([this%i, that])
end procedure

module procedure concatenate_IntArray1Ds_IntArray2D
  output = array([this, that%i])
end procedure

module procedure concatenate_IntArray2D_IntArray2D
  output = array([this%i, that%i])
end procedure

module procedure equality_IntArray1D_IntArray1D
  output = all(this%i==that%i)
end procedure

module procedure equality_IntArray2D_IntArray2D
  output = all(this%i==that%i)
end procedure

module procedure non_equality_IntArray1D_IntArray1D
  output = .not. this==that
end procedure

module procedure non_equality_IntArray2D_IntArray2D
  output = .not. this==that
end procedure

module procedure read_IntArray1D
  select type(this); type is(IntArray1D)
    this = IntArray1D(int(split_line(input)))
  end select
end procedure

module procedure write_IntArray1D
  select type(this); type is(IntArray1D)
    output = join(this%i)
  end select
end procedure

module procedure new_IntArray1D_String
  call this%read(input)
end procedure

module procedure read_IntArray2D
  select type(this); type is(IntArray2D)
    this = IntArray2D(IntArray1D(input))
  class default
    call err()
  end select
end procedure

module procedure write_IntArray2D
  select type(this); type is(IntArray2D)
    output = str(this%i)
  class default
    call err()
  end select
end procedure

module procedure new_IntArray2D_Strings
  call this%read(input)
end procedure

module procedure new_IntArray2D_StringArray
  this = IntArray2D(str(input))
end procedure
end submodule
