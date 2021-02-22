submodule (caesar_fraction_module) caesar_fraction_submodule
  use caesar_algebra_module
contains

module procedure new_IntFraction
  this%n_ = numerator
  this%d_ = denominator
  call this%simplify()
end procedure

module procedure numerator
  output = this%n_
end procedure

module procedure denominator
  output = this%d_
end procedure

module procedure simplify
  integer :: common_denominator
  
  common_denominator = gcd(this%n_, this%d_)
  this%n_ = this%n_ / common_denominator
  this%d_ = this%d_ / common_denominator
  
  if (this%d_<0) then
    this%n_ = - this%n_
    this%d_ = - this%d_
  endif
end procedure

module procedure int_IntFraction
  output = this%n_/this%d_
end procedure

module procedure dble_IntFraction
  output = real(this%n_,dp) / this%d_
end procedure

module procedure frac_character
  output = frac(str(input))
end procedure

module procedure frac_String
  output = IntFraction(input)
end procedure

module procedure frac_integer
  output = IntFraction(input,1)
end procedure

module procedure frac_integers
  output = IntFraction(numerator,denominator)
end procedure

module procedure equality_IntFraction_IntFraction
  output = this%n_==that%n_ .and. this%d_==that%d_
end procedure

module procedure equality_IntFraction_integer
  output = this%n_==that .and. this%d_==1
end procedure

module procedure equality_integer_IntFraction
  output = this==that%n_ .and. that%d_==1
end procedure

module procedure non_equality_IntFraction_IntFraction
  output = .not. this==that
end procedure

module procedure non_equality_IntFraction_integer
  output = .not. this==that
end procedure

module procedure non_equality_integer_IntFraction
  output = .not. this==that
end procedure

module procedure lt_IntFraction_IntFraction
  output = this%n_*that%d_ < that%n_*this%d_
end procedure

module procedure lt_IntFraction_integer
  output = this%n_ < that*this%d_
end procedure

module procedure lt_integer_IntFraction
  output = this*that%d_ < that%n_
end procedure

module procedure gt_IntFraction_IntFraction
  output = this%n_*that%d_ > that%n_*this%d_
end procedure

module procedure gt_IntFraction_integer
  output = this%n_ > that*this%d_
end procedure

module procedure gt_integer_IntFraction
  output = this*that%d_ > that%n_
end procedure

module procedure le_IntFraction_IntFraction
  output = .not. this > that
end procedure

module procedure le_IntFraction_integer
  output = .not. this > that
end procedure

module procedure le_integer_IntFraction
  output = .not. this > that
end procedure

module procedure ge_IntFraction_IntFraction
  output = .not. this < that
end procedure

module procedure ge_IntFraction_integer
  output = .not. this < that
end procedure

module procedure ge_integer_IntFraction
  output = .not. this < that
end procedure

module procedure add_IntFraction_IntFraction
  output = IntFraction(this%n_*that%d_+that%n_*this%d_, this%d_*that%d_)
end procedure

module procedure add_IntFraction_integer
  output = IntFraction(this%n_+that*this%d_, this%d_)
end procedure

module procedure add_integer_IntFraction
  output = IntFraction(this*that%d_+that%n_, that%d_)
end procedure

module procedure subtract_IntFraction_IntFraction
  output = IntFraction(this%n_*that%d_-that%n_*this%d_, this%d_*that%d_)
end procedure

module procedure subtract_IntFraction_integer
  output = IntFraction(this%n_-that*this%d_, this%d_)
end procedure

module procedure subtract_integer_IntFraction
  output = IntFraction(this*that%d_-that%n_, that%d_)
end procedure

module procedure multiply_IntFraction_IntFraction
  output = IntFraction(this%n_*that%n_, this%d_*that%d_)
end procedure

module procedure multiply_IntFraction_integer
  output = IntFraction(this%n_*that, this%d_)
end procedure

module procedure multiply_integer_IntFraction
  output = IntFraction(this*that%n_, that%d_)
end procedure

module procedure divide_IntFraction_IntFraction
  output = IntFraction(this%n_*that%d_, this%d_*that%n_)
end procedure

module procedure divide_IntFraction_integer
  output = IntFraction(this%n_, this%d_*that)
end procedure

module procedure divide_integer_IntFraction
  output = IntFraction(this*that%d_, that%n_)
end procedure

module procedure is_int_IntFraction
  output = this%d_==1
end procedure

module procedure modulo_IntFraction_integer
  output = IntFraction(modulo(this%n_,this%d_*that), this%d_)
end procedure

module procedure negative_IntFraction
  output = IntFraction(-this%n_, this%d_)
end procedure

module procedure abs_IntFraction
  output = IntFraction(abs(this%n_), this%d_)
end procedure

module procedure read_IntFraction
  type(String), allocatable :: split_string(:)
  
  select type(this); type is(IntFraction)
    split_string = split_line(input, '/')
    if (size(split_string)==1) then
      ! Assume the string is an integer.
      this = frac(int(split_string(1)))
    elseif (size(split_string)==2) then
      ! Assume the string is of the form 'a/b'.
      this = IntFraction(int(split_string(1)),int(split_string(2)))
    else
      call print_line('Error parsing fraction from string: '//input)
      call err()
    endif
  class default
    call err()
  end select
end procedure

module procedure write_IntFraction
  select type(this); type is(IntFraction)
    if (is_int(this)) then
     output = pad_int_to_str(this%n_)
    else
      output = pad_int_to_str(this%n_)//'/'//this%d_
    endif
  class default
    call err()
  end select
end procedure

module procedure new_IntFraction_String
  call this%read(input)
end procedure
end submodule
