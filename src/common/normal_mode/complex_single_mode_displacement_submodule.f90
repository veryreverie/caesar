submodule (caesar_complex_single_mode_displacement_module) caesar_complex_single_mode_displacement_submodule
  use caesar_normal_mode_module
contains

module procedure new_ComplexSingleDisplacement
  this%id        = id
  this%magnitude = magnitude
end procedure

module procedure new_ComplexSingleDisplacement_ComplexMode
  this = ComplexSingleDisplacement(id=mode%id, magnitude=magnitude)
end procedure

module procedure multiply_real_ComplexSingleDisplacement
  output = ComplexSingleDisplacement( id        = that%id,            &
                                    & magnitude = this*that%magnitude )
end procedure

module procedure multiply_ComplexSingleDisplacement_real
  output = ComplexSingleDisplacement( id        = this%id,            &
                                    & magnitude = this%magnitude*that )
end procedure

module procedure multiply_complex_ComplexSingleDisplacement
  output = ComplexSingleDisplacement( id        = that%id,            &
                                    & magnitude = this*that%magnitude )
end procedure

module procedure multiply_ComplexSingleDisplacement_complex
  output = ComplexSingleDisplacement( id        = this%id,            &
                                    & magnitude = this%magnitude*that )
end procedure

module procedure divide_ComplexSingleDisplacement_complex
  output = ComplexSingleDisplacement( id        = this%id,            &
                                    & magnitude = this%magnitude/that )
end procedure

module procedure add_ComplexSingleDisplacement_ComplexSingleDisplacement
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to add vectors along different modes.')
    call err()
  endif
  
  output = ComplexSingleDisplacement(            &
     & id        = this%id,                      &
     & magnitude = this%magnitude+that%magnitude )
end procedure

module procedure negative_ComplexSingleDisplacement
  output = ComplexSingleDisplacement(id=this%id, magnitude=-this%magnitude)
end procedure

module procedure subtract_ComplexSingleDisplacement_ComplexSingleDisplacement
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to subtract vectors along different &
       &modes.')
    call err()
  endif
  
  output = ComplexSingleDisplacement(            &
     & id        = this%id,                      &
     & magnitude = this%magnitude-that%magnitude )
end procedure

module procedure select_mode_ComplexSingleDisplacement
  output = modes(first(modes%id==displacement%id))
end procedure

module procedure select_modes_ComplexSingleDisplacements
  integer :: i,ialloc
  
  allocate(output(size(displacements)), stat=ialloc); call err(ialloc)
  do i=1,size(displacements)
    output(i) = select_mode(displacements(i), modes)
  enddo
end procedure

module procedure read_ComplexSingleDisplacement
  type(String), allocatable :: split_string(:)
  integer                   :: id
  complex(dp)               :: magnitude
  
  select type(this); type is(ComplexSingleDisplacement)
    split_string = split_line(input)
    if (size(split_string)/=3) then
      call print_line(ERROR//': unable to parse complex single mode &
         &vector from string: '//input)
      call err()
    endif
    
    ! If e.g. id=3 and power=2.1+1.2i then split_string = ["u3","=","2.1+1.2i"]
    ! The 'u' needs stripping off the first element to give the id.
    id = int(slice(split_string(1),2,len(split_string(1))))
    magnitude = cmplx(split_string(3))
    
    this = ComplexSingleDisplacement(id,magnitude)
  end select
end procedure

module procedure write_ComplexSingleDisplacement
  select type(this); type is(ComplexSingleDisplacement)
    output = 'u'//this%id//' = '//this%magnitude
  end select
end procedure

module procedure new_ComplexSingleDisplacement_String
  call this%read(input)
end procedure
end submodule
