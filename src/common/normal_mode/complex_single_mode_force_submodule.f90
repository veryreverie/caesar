submodule (caesar_complex_single_mode_force_module) caesar_complex_single_mode_force_submodule
  use caesar_normal_mode_module
contains

module procedure new_ComplexSingleForce
  this%id        = id
  this%magnitude = magnitude
end procedure

module procedure new_ComplexSingleForce_ComplexMode
  this = ComplexSingleForce(id=mode%id, magnitude=magnitude)
end procedure

module procedure multiply_real_ComplexSingleForce
  output = ComplexSingleForce( id        = that%id,            &
                             & magnitude = this*that%magnitude )
end procedure

module procedure multiply_ComplexSingleForce_real
  output = ComplexSingleForce( id        = this%id,            &
                             & magnitude = this%magnitude*that )
end procedure

module procedure multiply_complex_ComplexSingleForce
  output = ComplexSingleForce( id        = that%id,            &
                             & magnitude = this*that%magnitude )
end procedure

module procedure multiply_ComplexSingleForce_complex
  output = ComplexSingleForce( id        = this%id,            &
                             & magnitude = this%magnitude*that )
end procedure

module procedure divide_ComplexSingleForce_complex
  output = ComplexSingleForce( id        = this%id,            &
                             & magnitude = this%magnitude/that )
end procedure

module procedure add_ComplexSingleForce_ComplexSingleForce
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to add vectors along different modes.')
    call err()
  endif
  
  output = ComplexSingleForce( id        = this%id,                      &
                             & magnitude = this%magnitude+that%magnitude )
end procedure

module procedure negative_ComplexSingleForce
  output = ComplexSingleForce(id=this%id, magnitude=-this%magnitude)
end procedure

module procedure subtract_ComplexSingleForce_ComplexSingleForce
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to subtract vectors along different &
       &modes.')
    call err()
  endif
  
  output = ComplexSingleForce( id        = this%id,                      &
                             & magnitude = this%magnitude-that%magnitude )
end procedure

module procedure select_mode_ComplexSingleForce
  output = modes(first(modes%id==force%id))
end procedure

module procedure select_modes_ComplexSingleForces
  integer :: i,ialloc
  
  allocate(output(size(forces)), stat=ialloc); call err(ialloc)
  do i=1,size(forces)
    output(i) = select_mode(forces(i), modes)
  enddo
end procedure

module procedure read_ComplexSingleForce
  type(String), allocatable :: split_string(:)
  integer                   :: id
  complex(dp)               :: magnitude
  
  select type(this); type is(ComplexSingleForce)
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
    
    this = ComplexSingleForce(id,magnitude)
  end select
end procedure

module procedure write_ComplexSingleForce
  select type(this); type is(ComplexSingleForce)
    output = 'u'//this%id//' = '//this%magnitude
  end select
end procedure

module procedure new_ComplexSingleForce_String
  call this%read(input)
end procedure
end submodule
