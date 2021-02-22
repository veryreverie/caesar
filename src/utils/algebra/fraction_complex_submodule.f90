submodule (caesar_fraction_complex_module) caesar_fraction_complex_submodule
  use caesar_algebra_module
contains

module procedure new_FractionComplex
  this%real_ = real_part
  if (present(imaginary_part)) then
    this%imag_ = imaginary_part
  else
    this%imag_ = frac(0)
  endif
end procedure

module procedure real_FractionComplex
  output = this%real_
end procedure

module procedure aimag_FractionComplex
  output = this%imag_
end procedure

module procedure conjg_FractionComplex
  output = FractionComplex(real(this),-aimag(this))
end procedure

module procedure add_FractionComplex_FractionComplex
  output%real_ = real(this) + real(that)
  output%imag_ = aimag(this) + aimag(that)
end procedure

module procedure add_FractionComplex_IntComplex
  output%real_ = real(this) + real(that)
  output%imag_ = aimag(this) + aimag(that)
end procedure

module procedure add_IntComplex_FractionComplex
  output%real_ = real(this) + real(that)
  output%imag_ = aimag(this) + aimag(that)
end procedure

module procedure add_FractionComplex_IntFraction
  output%real_ = real(this) + that
  output%imag_ = aimag(this)
end procedure

module procedure add_IntFraction_FractionComplex
  output%real_ = this + real(that)
  output%imag_ =        aimag(that)
end procedure

module procedure add_FractionComplex_integer
  output%real_ = real(this) + that
  output%imag_ = aimag(this)
end procedure

module procedure add_integer_FractionComplex
  output%real_ = this + real(that)
  output%imag_ =        aimag(that)
end procedure

module procedure subtract_FractionComplex_FractionComplex
  output%real_ = real(this) - real(that)
  output%imag_ = aimag(this) - aimag(that)
end procedure

module procedure subtract_FractionComplex_IntComplex
  output%real_ = real(this) - real(that)
  output%imag_ = aimag(this) - aimag(that)
end procedure

module procedure subtract_IntComplex_FractionComplex
  output%real_ = real(this) - real(that)
  output%imag_ = aimag(this) - aimag(that)
end procedure

module procedure subtract_FractionComplex_IntFraction
  output%real_ = real(this) - that
  output%imag_ = aimag(this)
end procedure

module procedure subtract_IntFraction_FractionComplex
  output%real_ = this - real(that)
  output%imag_ =      - aimag(that)
end procedure

module procedure subtract_FractionComplex_integer
  output%real_ = real(this) - that
  output%imag_ = aimag(this)
end procedure

module procedure subtract_integer_FractionComplex
  output%real_ = this - real(that)
  output%imag_ =      - aimag(that)
end procedure

module procedure multiply_FractionComplex_FractionComplex
  output%real_ = real(this)*real(that) - aimag(this)*aimag(that)
  output%imag_ = real(this)*aimag(that) + aimag(this)*real(that)
end procedure

module procedure multiply_FractionComplex_IntComplex
  output%real_ = real(this)*real(that) - aimag(this)*aimag(that)
  output%imag_ = real(this)*aimag(that) + aimag(this)*real(that)
end procedure

module procedure multiply_IntComplex_FractionComplex
  output%real_ = real(this)*real(that) - aimag(this)*aimag(that)
  output%imag_ = real(this)*aimag(that) + aimag(this)*real(that)
end procedure

module procedure multiply_FractionComplex_IntFraction
  output%real_ = real(this) * that
  output%imag_ = aimag(this) * that
end procedure

module procedure multiply_IntFraction_FractionComplex
  output%real_ = this * real(that)
  output%imag_ = this * aimag(that)
end procedure

module procedure multiply_FractionComplex_integer
  output%real_ = real(this) * that
  output%imag_ = aimag(this) * that
end procedure

module procedure multiply_integer_FractionComplex
  output%real_ = this * real(that)
  output%imag_ = this * aimag(that)
end procedure

module procedure read_FractionComplex
  type(String), allocatable :: split_string(:)
  type(IntFraction)         :: real_part
  type(IntFraction)         :: imag_part
  
  select type(this); type is(FractionComplex)
    split_string = split_line(input,'+')
    if (size(split_string)==2) then
      ! input is of the form "a+bi" or "-a+bi".
      real_part = frac(split_string(1))
      imag_part = frac(slice(split_string(2),1,len(split_string(2))-1))
    elseif (size(split_string)==1) then
      ! input is of the form "a-bi" or "-a-bi".
      split_string = split_line(input,'-')
      if (size(split_string)==2) then
        if (slice(input,1,1)=='-') then
          real_part = -frac(split_string(1))
        else
          real_part = frac(split_string(1))
        endif
        imag_part = frac(slice(split_string(2),1,len(split_string(2))-1))
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
    
    this = FractionComplex(real_part,imag_part)
  end select
end procedure

module procedure write_FractionComplex
  select type(this); type is(FractionComplex)
    ! N.B. abs() is called in both cases to stop +/-0 from breaking formatting.
    if (aimag(this)>=0) then
      output = real(this)//'+'//abs(aimag(this))//'i'
    else
      output = real(this)//'-'//abs(aimag(this))//'i'
    endif
  end select
end procedure

module procedure new_FractionComplex_String
  call this%read(input)
end procedure
end submodule
