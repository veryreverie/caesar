submodule (caesar_real_mode_displacement_module) caesar_real_mode_displacement_submodule
  use caesar_normal_mode_module
contains

module procedure new_RealModeDisplacement
  this%vectors = displacements
end procedure

module procedure new_RealModeDisplacement_RealModes
  this = RealModeDisplacement(RealSingleDisplacement(modes,displacements))
end procedure

module procedure size_RealModeDisplacement
  output = size(this%vectors)
end procedure

module procedure multiply_real_RealModeDisplacement
  output = RealModeDisplacement(this*that%vectors)
end procedure

module procedure multiply_RealModeDisplacement_real
  output = RealModeDisplacement(this%vectors*that)
end procedure

module procedure divide_RealModeDisplacement_real
  output = RealModeDisplacement(this%vectors/that)
end procedure

module procedure add_RealModeDisplacement_RealModeDisplacement
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

module procedure sum_RealModeDisplacements
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

module procedure negative_RealModeDisplacement
  output = RealModeDisplacement(-this%vectors)
end procedure

module procedure subtract_RealModeDisplacement_RealModeDisplacement
  output = this + (-that)
end procedure

module procedure new_MassWeightedDisplacement_RealModeDisplacement
  type(RealMode),   allocatable :: selected_modes(:)
  type(QpointData), allocatable :: selected_qpoints(:)
  
  if (size(this)==0) then
    output = MassWeightedDisplacement(structure)
  else
    selected_modes = select_modes(this%vectors, modes)
    selected_qpoints = select_qpoints(selected_modes, qpoints)
    output = sum(MassWeightedDisplacement( this%vectors,    &
                                         & selected_modes,  &
                                         & structure,       &
                                         & selected_qpoints ))
  endif
end procedure

module procedure new_CartesianDisplacement_RealModeDisplacement
  type(RealMode),   allocatable :: selected_modes(:)
  type(QpointData), allocatable :: selected_qpoints(:)
  
  if (size(this)==0) then
    output = CartesianDisplacement(structure)
  else
    selected_modes = select_modes(this%vectors, modes)
    selected_qpoints = select_qpoints(selected_modes, qpoints)
    output = sum(CartesianDisplacement( this%vectors,    &
                                      & selected_modes,  &
                                      & structure,       &
                                      & selected_qpoints ))
  endif
end procedure

module procedure new_RealModeDisplacement_MassWeightedDisplacement
  type(QpointData), allocatable :: selected_qpoints(:)
  
  selected_qpoints = select_qpoints(modes, qpoints)
  
  this = RealModeDisplacement(RealSingleDisplacement( modes,           &
                                                    & displacement,    &
                                                    & structure,       &
                                                    & selected_qpoints ))
end procedure

module procedure new_RealModeDisplacement_CartesianDisplacement
  type(QpointData), allocatable :: selected_qpoints(:)
  
  selected_qpoints = select_qpoints(modes, qpoints)
  
  this = RealModeDisplacement(RealSingleDisplacement( modes,           &
                                                    & displacement,    &
                                                    & structure,       &
                                                    & selected_qpoints ))
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

module procedure read_RealModeDisplacement
  select type(this); type is(RealModeDisplacement)
    this = RealModeDisplacement(RealSingleDisplacement(input))
  class default
    call err()
  end select
end procedure

module procedure write_RealModeDisplacement
  select type(this); type is(RealModeDisplacement)
    output = str(this%vectors)
  class default
    call err()
  end select
end procedure

module procedure new_RealModeDisplacement_Strings
  call this%read(input)
end procedure

module procedure new_RealModeDisplacement_StringArray
  this = RealModeDisplacement(str(input))
end procedure
end submodule
