submodule (caesar_complex_mode_force_module) caesar_complex_mode_force_submodule
  use caesar_normal_mode_module
contains

module procedure new_ComplexModeForce
  this%vectors = forces
end procedure

module procedure new_ComplexModeForce_ComplexModes
  this = ComplexModeForce(ComplexSingleForce(modes,forces))
end procedure

module procedure size_ComplexModeForce
  output = size(this%vectors)
end procedure

module procedure multiply_real_ComplexModeForce
  output = ComplexModeForce(this*that%vectors)
end procedure

module procedure multiply_ComplexModeForce_real
  output = ComplexModeForce(this%vectors*that)
end procedure

module procedure multiply_complex_ComplexModeForce
  output = ComplexModeForce(this*that%vectors)
end procedure

module procedure multiply_ComplexModeForce_complex
  output = ComplexModeForce(this%vectors*that)
end procedure

module procedure divide_ComplexModeForce_complex
  output = ComplexModeForce(this%vectors/that)
end procedure

module procedure add_ComplexModeForce_ComplexModeForce
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

module procedure sum_ComplexModeForces
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

module procedure negative_ComplexModeForce
  output = ComplexModeForce(-this%vectors)
end procedure

module procedure subtract_ComplexModeForce_ComplexModeForce
  output = this + (-that)
end procedure

module procedure force_id
  integer :: i
  
  i = first(this%vectors%id==id, default=0)
  if (i==0) then
    output = 0.0_dp
  else
    output = this%vectors(i)%magnitude
  endif
end procedure

module procedure force_mode
  output = this%force(mode%id)
end procedure

module procedure read_ComplexModeForce
  select type(this); type is(ComplexModeForce)
    this = ComplexModeForce(ComplexSingleForce(input))
  class default
    call err()
  end select
end procedure

module procedure write_ComplexModeForce
  select type(this); type is(ComplexModeForce)
    output = str(this%vectors)
  class default
    call err()
  end select
end procedure

module procedure new_ComplexModeForce_Strings
  call this%read(input)
end procedure

module procedure new_ComplexModeForce_StringArray
  this = ComplexModeForce(str(input))
end procedure
end submodule
