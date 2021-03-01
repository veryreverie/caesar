submodule (caesar_integrated_module) caesar_integrated_submodule
  use caesar_integration_module
contains

module procedure new_Integrated_scalar
  this%scalar_ = input
end procedure

module procedure new_Integrated_vector
  this%vector_ = input
end procedure

module procedure new_Integrated_tensor
  this%tensor_ = input
end procedure

module procedure scalar_Integrated
  if (allocated(this%scalar_)) then
    output = this%scalar_
  else
    call err()
  endif
end procedure

module procedure vector_Integrated
  if (allocated(this%vector_)) then
    output = this%vector_
  else
    call err()
  endif
end procedure

module procedure tensor_Integrated
  if (allocated(this%tensor_)) then
    output = this%tensor_
  else
    call err()
  endif
end procedure

module procedure add_Integrated_Integrated
  if (allocated(this%scalar_)) then
    output = Integrated(this%scalar_+that%scalar_)
  elseif (allocated(this%vector_)) then
    output = Integrated(this%vector_+that%vector_)
  elseif (allocated(this%tensor_)) then
    output = Integrated(this%tensor_+that%tensor_)
  else
    call err()
  endif
end procedure

module procedure sum_Integrateds
  integer :: i
  
  if (size(input)>0) then
    output = input(1)
    do i=2,size(input)
      output = output + input(i)
    enddo
  else
    call err()
  endif
end procedure
end submodule
