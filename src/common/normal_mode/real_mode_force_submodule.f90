submodule (caesar_real_mode_force_module) caesar_real_mode_force_submodule
  use caesar_normal_mode_module
contains

module procedure new_RealModeForce
  this%vectors = forces
end procedure

module procedure new_RealModeForce_RealModes
  this = RealModeForce(RealSingleForce(modes,forces))
end procedure

module procedure size_RealModeForce
  output = size(this%vectors)
end procedure

module procedure multiply_real_RealModeForce
  output = RealModeForce(this*that%vectors)
end procedure

module procedure multiply_RealModeForce_real
  output = RealModeForce(this%vectors*that)
end procedure

module procedure divide_RealModeForce_real
  output = RealModeForce(this%vectors/that)
end procedure

module procedure add_RealModeForce_RealModeForce
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

module procedure sum_RealModeForces
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

module procedure negative_RealModeForce
  output = RealModeForce(-this%vectors)
end procedure

module procedure subtract_RealModeForce_RealModeForce
  output = this + (-that)
end procedure

module procedure new_MassWeightedForce_RealModeForce
  type(RealMode),   allocatable :: selected_modes(:)
  type(QpointData), allocatable :: selected_qpoints(:)
  
  selected_modes = select_modes(this%vectors, modes)
  selected_qpoints = select_qpoints(selected_modes, qpoints)
  output = sum(MassWeightedForce( this%vectors,    &
                                & selected_modes,  &
                                & structure,       &
                                & selected_qpoints ))
end procedure

module procedure new_CartesianForce_RealModeForce
  type(RealMode),   allocatable :: selected_modes(:)
  type(QpointData), allocatable :: selected_qpoints(:)
  
  selected_modes = select_modes(this%vectors, modes)
  selected_qpoints = select_qpoints(selected_modes, qpoints)
  output = sum(CartesianForce( this%vectors,    &
                             & selected_modes,  &
                             & structure,       &
                             & selected_qpoints ))
end procedure

module procedure new_RealModeForce_MassWeightedForce
  type(QpointData), allocatable :: selected_qpoints(:)
  
  selected_qpoints = select_qpoints(modes, qpoints)
  
  this = RealModeForce(RealSingleForce( modes,           &
                                      & force,           &
                                      & structure,       &
                                      & selected_qpoints ))
end procedure

module procedure new_RealModeForce_CartesianForce
  type(QpointData), allocatable :: selected_qpoints(:)
  
  selected_qpoints = select_qpoints(modes, qpoints)
  
  this = RealModeForce(RealSingleForce( modes,           &
                                      & force,           &
                                      & structure,       &
                                      & selected_qpoints ))
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

module procedure read_RealModeForce
  select type(this); type is(RealModeForce)
    this = RealModeForce(RealSingleForce(input))
  class default
    call err()
  end select
end procedure

module procedure write_RealModeForce
  select type(this); type is(RealModeForce)
    output = str(this%vectors)
  class default
    call err()
  end select
end procedure

module procedure new_RealModeForce_Strings
  call this%read(input)
end procedure

module procedure new_RealModeForce_StringArray
  this = RealModeForce(str(input))
end procedure
end submodule
