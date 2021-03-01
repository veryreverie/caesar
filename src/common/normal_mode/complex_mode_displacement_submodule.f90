submodule (caesar_complex_mode_displacement_module) caesar_complex_mode_displacement_submodule
  use caesar_normal_mode_module
contains

module procedure new_ComplexModeDisplacement
  this%vectors = displacements
end procedure

module procedure new_ComplexModeDisplacement_ComplexModes
  this = ComplexModeDisplacement(ComplexSingleDisplacement( modes,        &
                                                          & displacements ))
end procedure

module procedure size_ComplexModeDisplacement
  output = size(this%vectors)
end procedure

module procedure multiply_real_ComplexModeDisplacement
  output = ComplexModeDisplacement(this*that%vectors)
end procedure

module procedure multiply_ComplexModeDisplacement_real
  output = ComplexModeDisplacement(this%vectors*that)
end procedure

module procedure multiply_complex_ComplexModeDisplacement
  output = ComplexModeDisplacement(this*that%vectors)
end procedure

module procedure multiply_ComplexModeDisplacement_complex
  output = ComplexModeDisplacement(this%vectors*that)
end procedure

module procedure divide_ComplexModeDisplacement_complex
  output = ComplexModeDisplacement(this%vectors/that)
end procedure

module procedure add_ComplexModeDisplacement_ComplexModeDisplacement
  integer :: i,j
  
  output = this
  do i=1,size(that)
    j = first(this%vectors%id==that%vectors(i)%id, default=0)
    if (j==0) then
      output%vectors = [output%vectors, that%vectors(i)]
    else
      output%vectors(j) = output%vectors(j) + that%vectors(i)
    endif
  enddo
end procedure

module procedure sum_ComplexModeDisplacements
  integer :: i
  
  if (size(this)==0) then
    call print_line(ERROR//': Trying to sum an empty list.')
    call err()
  endif
  
  output = this(1)
  do i=2,size(this)
    output = output + this(i)
  enddo
end procedure

module procedure negative_ComplexModeDisplacement
  output = ComplexModeDisplacement(-this%vectors)
end procedure

module procedure subtract_ComplexModeDisplacement_ComplexModeDisplacement
  output = this + (-that)
end procedure

module procedure displacement_id
  integer :: i
  
  i = first(this%vectors%id==id, default=0)
  if (i==0) then
    output = 0.0_dp
  else
    output = this%vectors(i)%magnitude
  endif
end procedure

module procedure displacement_mode
  output = this%displacement(mode%id)
end procedure

module procedure read_ComplexModeDisplacement
  select type(this); type is(ComplexModeDisplacement)
    this = ComplexModeDisplacement(ComplexSingleDisplacement(input))
  class default
    call err()
  end select
end procedure

module procedure write_ComplexModeDisplacement
  select type(this); type is(ComplexModeDisplacement)
    output = str(this%vectors)
  class default
    call err()
  end select
end procedure

module procedure new_ComplexModeDisplacement_Strings
  call this%read(input)
end procedure

module procedure new_ComplexModeDisplacement_StringArray
  this = ComplexModeDisplacement(str(input))
end procedure
end submodule
