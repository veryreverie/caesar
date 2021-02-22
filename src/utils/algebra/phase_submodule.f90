submodule (caesar_phase_module) caesar_phase_submodule
  use caesar_algebra_module
contains

module procedure new_PhaseData
  this%fraction = input
end procedure

module procedure cmplx_PhaseData
  real(dp) :: exponent
  
  exponent = 2*PI*dble(this%fraction)
  output = cmplx(cos(exponent),sin(exponent),dp)
end procedure

module procedure calculate_phase
  real(dp)          :: phase_real
  type(IntFraction) :: phase_frac
  
  phase_real = atan2(aimag(input), real(input)) / (2*PI)
  phase_frac = IntFraction(nint(phase_real*denom),denom)
  if (abs(dble(phase_frac)-phase_real)>0.01_dp) then
    call print_line(ERROR//': Phase incompatible with given denominator.')
    call print_line('Phase: '//phase_real)
    call print_line('Denominator: '//denom)
    call err()
  endif
  output = PhaseData(phase_frac)
end procedure

module procedure equality_PhaseData_PhaseData
  output = this%fraction==that%fraction
end procedure

module procedure non_equality_PhaseData_PhaseData
  output = .not. this==that
end procedure

module procedure read_PhaseData
  select type(this); type is(PhaseData)
    this = PhaseData(frac(slice(input,10,len(input)-1)))
  end select
end procedure

module procedure write_PhaseData
  select type(this); type is(PhaseData)
    output = 'exp(2pii*'//this%fraction//')'
  end select
end procedure

module procedure new_PhaseData_String
  call this%read(input)
end procedure
end submodule
