submodule (caesar_stress_dynamical_matrix_module) caesar_stress_dynamical_matrix_submodule
  use caesar_harmonic_stress_module
contains

module procedure new_StressDynamicalMatrix
  this%elements = elements
end procedure

module procedure new_StressDynamicalMatrix_zeroes
  integer :: ialloc
  
  allocate(this%elements(3,3), stat=ialloc); call err(ialloc)
  this%elements = DynamicalMatrix(no_atoms)
end procedure

module procedure expectation_StressDynamicalMatrix
  real(dp), allocatable :: elements(:,:)
  
  integer :: i,j,ialloc
  
  allocate(elements(3,3), stat=ialloc); call err(ialloc)
  elements = 0
  do i=1,3
    do j=1,3
      elements(j,i) = this%elements(j,i)%expectation(mode)
    enddo
  enddo
  
  output = mat(elements)
end procedure

module procedure add_StressDynamicalMatrix_StressDynamicalMatrix
  output = StressDynamicalMatrix(this%elements + that%elements)
end procedure

module procedure negative_StressDynamicalMatrix
  output = StressDynamicalMatrix(-this%elements)
end procedure

module procedure subtract_StressDynamicalMatrix_StressDynamicalMatrix
  output = StressDynamicalMatrix(this%elements - that%elements)
end procedure

module procedure multiply_StressDynamicalMatrix_real
  output = StressDynamicalMatrix(this%elements*that)
end procedure

module procedure multiply_real_StressDynamicalMatrix
  output = StressDynamicalMatrix(this*that%elements)
end procedure

module procedure multiply_StressDynamicalMatrix_complex
  output = StressDynamicalMatrix(this%elements*that)
end procedure

module procedure multiply_complex_StressDynamicalMatrix
  output = StressDynamicalMatrix(this*that%elements)
end procedure

module procedure divide_StressDynamicalMatrix_real
  output = StressDynamicalMatrix(this%elements/that)
end procedure

module procedure divide_StressDynamicalMatrix_complex
  output = StressDynamicalMatrix(this%elements/that)
end procedure

module procedure conjg_StressDynamicalMatrix
  output = StressDynamicalMatrix(conjg(input%elements))
end procedure

module procedure read_StressDynamicalMatrix
  type(StringArray), allocatable :: sections(:)
  
  type(DynamicalMatrix), allocatable :: elements(:,:)
  
  integer :: i
  
  select type(this); type is(StressDynamicalMatrix)
    sections = split_into_sections( input,                         &
                                  & separating_line=repeat('-',50) )
    do i=1,size(sections)
      sections(i)%strings = sections(i)%strings(2:)
    enddo
    elements = transpose(reshape(DynamicalMatrix(sections), [3,3]))
    this = StressDynamicalMatrix(elements)
  class default
    call err()
  end select
end procedure

module procedure write_StressDynamicalMatrix
  character(1), parameter :: directions(3) = ['x', 'y', 'z']
  
  integer :: i,j
  
  select type(this); type is(StressDynamicalMatrix)
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

module procedure new_StressDynamicalMatrix_Strings
  call this%read(input)
end procedure

module procedure new_StressDynamicalMatrix_StringArray
  this = StressDynamicalMatrix(str(input))
end procedure
end submodule
