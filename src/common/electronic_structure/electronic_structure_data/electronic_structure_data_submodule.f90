submodule (caesar_electronic_structure_data_module) caesar_electronic_structure_data_submodule
  use caesar_electronic_structure_data_module
contains

module procedure new_ElectronicStructure
  this%energy_ = energy
  
  if (present(forces)) then
    this%forces_ = forces
  endif
  
  if (present(hessian)) then
    this%hessian_ = hessian
  endif
  
  if (present(stress)) then
    this%stress_ = stress
  endif
  
  if (present(linear_response)) then
    this%linear_response_ = linear_response
  endif
end procedure

module procedure energy_ElectronicStructure
  output = this%energy_
end procedure

module procedure has_forces_ElectronicStructure
  output = allocated(this%forces_)
end procedure

module procedure forces_ElectronicStructure
  if (this%has_forces()) then
    output = this%forces_
  else
    call print_line(ERROR//': Sample result does not contain forces.')
    call err()
  endif
end procedure

module procedure has_hessian_ElectronicStructure
  output = allocated(this%hessian_)
end procedure

module procedure hessian_ElectronicStructure
  if (this%has_hessian()) then
    output = this%hessian_
  else
    call print_line(ERROR//': Sample result does not contain the Hessian.')
    call err()
  endif
end procedure

module procedure has_stress_ElectronicStructure
  output = allocated(this%stress_)
end procedure

module procedure stress_ElectronicStructure
  if (this%has_stress()) then
    output = this%stress_
  else
    call print_line(ERROR//': Sample result does not contain stress.')
    call err()
  endif
end procedure

module procedure has_linear_response_ElectronicStructure
  output = allocated(this%linear_response_)
end procedure

module procedure linear_response_ElectronicStructure
  if (this%has_linear_response()) then
    output = this%linear_response_
  else
    call print_line(ERROR//': Sample result does not contain linear_response.')
    call err()
  endif
end procedure

module procedure read_ElectronicStructure
  real(dp)                            :: energy
  type(CartesianForce),   allocatable :: forces
  type(CartesianHessian), allocatable :: hessian
  type(RealMatrix),       allocatable :: stress
  type(LinearResponse),   allocatable :: linear_response
  
  ! Line numbers.
  integer :: forces_start_line
  integer :: forces_end_line
  integer :: hessian_start_line
  integer :: hessian_end_line
  integer :: stress_start_line
  integer :: stress_end_line
  integer :: linear_response_start_line
  integer :: linear_response_end_line
  
  integer :: end_lines(5)
  
  integer :: i
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(ElectronicStructure)
    ! Find the lines which start each section.
    forces_start_line = 0
    hessian_start_line = 0
    stress_start_line = 0
    linear_response_start_line = 0
    do i=1,size(input)
      line = split_line(lower_case(input(i)))
      if (size(line)>=1) then
        if (line(1)=='forces') then
          forces_start_line = i+1
        elseif (line(1)=='hessian') then
          hessian_start_line = i+1
        elseif (line(1)=='stress' .or. line(1)=='virial') then
          stress_start_line = i+1
        elseif (line(1)=='permittivity') then
          linear_response_start_line = i+1
        endif
      endif
    enddo
    
    ! Find the lines which end each section.
    end_lines = [ forces_start_line-2,          &
                & hessian_start_line-2,         &
                & stress_start_line-2,          &
                & linear_response_start_line-2, &
                & size(input)                   ]
    if (forces_start_line/=0) then
      forces_end_line = minval( end_lines,                          &
                              & 1,                                  &
                              & mask = end_lines>=forces_start_line )
    endif
    
    if (hessian_start_line/=0) then
      hessian_end_line = minval( end_lines,                           &
                               & 1,                                   &
                               & mask = end_lines>=hessian_start_line )
    endif
    
    if (stress_start_line/=0) then
      stress_end_line = minval( end_lines,                          &
                              & 1,                                  &
                              & mask = end_lines>=stress_start_line )
    endif
    
    if (linear_response_start_line/=0) then
      linear_response_end_line = minval(                &
         & end_lines,                                   &
         & 1,                                           &
         & mask = end_lines>=linear_response_start_line )
    endif
    
    ! Read file contents.
    energy = dble(input(2))
    if (forces_start_line/=0) then
      forces = CartesianForce(input(forces_start_line:forces_end_line))
    endif
    if (hessian_start_line/=0) then
      hessian = CartesianHessian(input(hessian_start_line:hessian_end_line))
    endif
    if (stress_start_line/=0) then
      stress = RealMatrix(input(stress_start_line:stress_end_line))
    endif
    if (linear_response_start_line/=0) then
      ! N.B. the linear response class expects its header line, hence the -1.
      linear_response = LinearResponse(input(                    &
         & linear_response_start_line-1:linear_response_end_line ))
    endif
    
    ! Construct output.
    this = ElectronicStructure( energy,         &
                              & forces,         &
                              & hessian,        &
                              & stress,         &
                              & linear_response )
  class default
    call err()
  end select
end procedure

module procedure write_ElectronicStructure
  select type(this); type is(ElectronicStructure)
    output = [ str('Energy (Hartree):'), &
             & str(this%energy_)         ]
    
    if (allocated(this%forces_)) then
      output = [ output,                        &
               & str('Forces (Hartree/Bohr):'), &
               & str(this%forces_)              ]
    endif
    
    if (allocated(this%hessian_)) then
      output = [ output,                           &
               & str('Hessian (Hartree/Bohr^2):'), &
               & str(this%hessian_)                ]
    endif
    
    if (allocated(this%stress_)) then
      output = [ output,                          &
               & str('Stress (Hartree/Bohr^3):'), &
               & str(this%stress_)                ]
    endif
    
    if (allocated(this%linear_response_)) then
      output = [ output,                    &
               & str(this%linear_response_) ]
    endif
  class default
    call err()
  end select
end procedure

module procedure new_ElectronicStructure_Strings
  call this%read(input)
end procedure

module procedure new_ElectronicStructure_StringArray
  this = ElectronicStructure(str(input))
end procedure
end submodule
