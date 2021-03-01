submodule (caesar_qpoint_module) caesar_qpoint_submodule
  use caesar_structure_module
contains

module procedure new_QpointData
  this%qpoint           = qpoint
  this%id               = id
  this%paired_qpoint_id = paired_qpoint_id
end procedure

module procedure pair_QpointData
  output = QpointData( qpoint           = -this%qpoint,          &
                     & id               = this%paired_qpoint_id, &
                     & paired_qpoint_id = this%id                )
  call translate_to_primitive(output)
end procedure

module procedure is_gvector
  output = is_int(this%qpoint)
end procedure

module procedure is_paired_qpoint
  output = this%id==this%paired_qpoint_id
end procedure

module procedure min_sc_size
  type(IntFraction) :: q(3)
  
  q = frac(this%qpoint)
  output = lcm( q(1)%denominator(), &
              & q(2)%denominator(), &
              & q(3)%denominator()  )
end procedure

module procedure translate_to_primitive
  type(IntFraction) :: qpoint(3)
  integer           :: i
  
  ! Convert the q-point to an array of fractions.
  qpoint = frac(this%qpoint)
  ! Translate all elements to [0,1).
  qpoint = modulo(qpoint,1)
  ! Translate all elements in [1/2,1) to [-1/2,0).
  do i=1,3
    if (qpoint(i)>IntFraction(1,2)) then
      qpoint(i) = qpoint(i) - 1
    endif
  enddo
  ! Convert back to a fraction vector.
  this%qpoint = vec(qpoint)
end procedure

module procedure equality_QpointData
  output = is_int(this%qpoint-that%qpoint)
end procedure

module procedure non_equality_QpointData
  output = .not. this==that
end procedure

module procedure read_QpointData
  type(String), allocatable :: line(:)
  
  type(FractionVector) :: qpoint
  integer              :: id
  integer              :: paired_qpoint_id
  
  select type(this); type is(QpointData)
    if (size(input)/=3) then
      call print_line(ERROR//': Unable to parse q-point from string input:')
      call print_lines(input)
      call err()
    endif
    
    line = split_line(input(1))
    id = int(line(2))
    
    line = split_line(input(2))
    qpoint = vec(frac(line(3:5)))
    
    line = split_line(input(3))
    paired_qpoint_id = int(line(11))
    
    this = QpointData(qpoint,id,paired_qpoint_id)
  class default
    call err()
  end select
end procedure

module procedure write_QpointData
  select type(this); type is(QpointData)
    output = [ 'q-point '//this%id,                                 &
             & 'q = '//this%qpoint,                                 &
             & "The ID of q' s.t. q+q' is a primitive G-vector: "// &
             &    this%paired_qpoint_id                             ]
  class default
    call err()
  end select
end procedure

module procedure new_QpointData_Strings
  call this%read(input)
end procedure

module procedure new_QpointData_StringArray
  this = QpointData(str(input))
end procedure
end submodule
