submodule (caesar_integer_complex_module) caesar_integer_complex_submodule
  use caesar_algebra_module
contains

module procedure new_IntComplex
  this%real_ = real_part
  if (present(imaginary_part)) then
    this%imag_ = imaginary_part
  else
    this%imag_ = 0
  endif
end procedure

module procedure real_IntComplex
  output = this%real_
end procedure

module procedure aimag_IntComplex
  output = this%imag_
end procedure

module procedure conjg_IntComplex
  output = IntComplex(real(this),-aimag(this))
end procedure

module procedure add_IntComplex_IntComplex
  output%real_ = real(this) + real(that)
  output%imag_ = aimag(this) + aimag(that)
end procedure

module procedure add_IntComplex_integer
  output%real_ = real(this) + that
  output%imag_ = aimag(this)
end procedure

module procedure add_integer_IntComplex
  output%real_ = this + real(that)
  output%imag_ =        aimag(that)
end procedure

module procedure negative_IntComplex
  output%real_ = -real(this)
  output%imag_ = -aimag(this)
end procedure

module procedure subtract_IntComplex_IntComplex
  output%real_ = real(this) - real(that)
  output%imag_ = aimag(this) - aimag(that)
end procedure

module procedure subtract_IntComplex_integer
  output%real_ = real(this) - that
  output%imag_ = aimag(this)
end procedure

module procedure subtract_integer_IntComplex
  output%real_ = this - real(that)
  output%imag_ =      - aimag(that)
end procedure

module procedure multiply_IntComplex_IntComplex
  output%real_ = real(this)*real(that) - aimag(this)*aimag(that)
  output%imag_ = real(this)*aimag(that) + aimag(this)*real(that)
end procedure

module procedure multiply_IntComplex_integer
  output%real_ = real(this) * that
  output%imag_ = aimag(this) * that
end procedure

module procedure multiply_integer_IntComplex
  output%real_ = this * real(that)
  output%imag_ = this * aimag(that)
end procedure

module procedure read_IntComplex
  type(String), allocatable :: split_string(:)
  integer                   :: real_part
  integer                   :: imag_part
  
  select type(this); type is(IntComplex)
    split_string = split_line(input,'+')
    if (size(split_string)==2) then
      ! input is of the form "a+bi" or "-a+bi".
      real_part = int(split_string(1))
      imag_part = int(slice(split_string(2),1,len(split_string(2))-1))
    elseif (size(split_string)==1) then
      ! input is of the form "a-bi" or "-a-bi".
      split_string = split_line(input,'-')
      if (size(split_string)==2) then
        if (slice(input,1,1)=='-') then
          real_part = -int(split_string(1))
        else
          real_part = int(split_string(1))
        endif
        imag_part = int(slice(split_string(2),1,len(split_string(2))-1))
      else
        call print_line(ERROR//': Unable to read IntComplex from string: '// &
           & input)
        call err()
      endif
    else
      call print_line(ERROR//': Unable to read IntComplex from string: '// &
         & input)
      call err()
    endif
    
    this = IntComplex(real_part,imag_part)
  end select
end procedure

module procedure write_IntComplex
  select type(this); type is(IntComplex)
    ! N.B. abs() is called in both cases to stop +/-0 from breaking formatting.
    if (aimag(this)>=0) then
      output = real(this)//'+'//abs(aimag(this))//'i'
    else
      output = real(this)//'-'//abs(aimag(this))//'i'
    endif
  end select
end procedure

module procedure new_IntComplex_String
  call this%read(input)
end procedure
end submodule
