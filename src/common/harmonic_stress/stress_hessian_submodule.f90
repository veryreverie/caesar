submodule (caesar_stress_hessian_module) caesar_stress_hessian_submodule
  use caesar_harmonic_stress_module
contains

module procedure new_StressHessian
  this%elements = elements
end procedure

module procedure read_StressHessian
  type(StringArray), allocatable :: sections(:)
  
  type(CartesianHessian), allocatable :: elements(:,:)
  
  integer :: i
  
  select type(this); type is(StressHessian)
    sections = split_into_sections( input,                         &
                                  & separating_line=repeat('-',50) )
    do i=1,size(sections)
      sections(i)%strings = sections(i)%strings(2:)
    enddo
    elements = transpose(reshape(CartesianHessian(sections), [3,3]))
    this = StressHessian(elements)
  class default
    call err()
  end select
end procedure

module procedure write_StressHessian
  character(1), parameter :: directions(3) = ['x', 'y', 'z']
  
  integer :: i,j
  
  select type(this); type is(StressHessian)
    output = [str(repeat('-',50))]
    do i=1,3
      do j=1,3
        output = [                                                        &
           & output,                                                      &
           & str('Stress component '//directions(i)//directions(j)//':'), &
           & str(this%elements(i,j)),                                     &
           & str(repeat('-',50))                                          ]
      enddo
    enddo
  class default
    call err()
  end select
end procedure

module procedure new_StressHessian_Strings
  call this%read(input)
end procedure

module procedure new_StressHessian_StringArray
  this = StressHessian(str(input))
end procedure
end submodule
