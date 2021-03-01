submodule (caesar_monomial_state_1d_module) caesar_monomial_state_1d_submodule
  use caesar_states_module
contains

module procedure new_MonomialState1D
  this%id_    = id
  this%power_ = power
end procedure

module procedure id_MonomialState1D
  output = this%id_
end procedure

module procedure power_MonomialState1D
  output = this%power_
end procedure

module procedure total_power_MonomialState1D
  output = this%power_
end procedure

module procedure wavevector_MonomialState1D
  integer :: mode
  integer :: qpoint
  
  mode   = first(modes%id==this%id_)
  qpoint = first(qpoints%id==modes(mode)%qpoint_id)
  output = qpoints(qpoint)%qpoint * this%power_
end procedure

module procedure finite_overlap_MonomialState1D
  integer :: p,q
  
  p = bra%power_
  q = ket%power_
  
  output = modulo(p+q, 2)==0
end procedure

module procedure inner_product_MonomialState1D
  integer :: p,q
  
  p = bra%power_
  q = ket%power_
  
  if (modulo(p+q,2)==1) then
    output = 0
  else
    output = exp( log_odd_factorial((p+q)/2)                         &
              & - 0.5_dp*(log_odd_factorial(p)+log_odd_factorial(q)) )
  endif
end procedure

module procedure braket_MonomialState1D
  integer :: p,q,n
  
  p = bra%power_
  q = ket%power_
  n = potential%power
  
  if (modulo(p+n+q,2)==1) then
    output = 0
  else
    output = exp( log_odd_factorial((p+n+q)/2)                       &
              & - 0.5_dp*(log_odd_factorial(p)+log_odd_factorial(q)) &
              & - 0.5_dp*n*log_2nw                                   )
  endif
end procedure

module procedure first_derivative_MonomialState1D
  integer :: p,q
  
  p = bra%power_
  q = ket%power_
  
  if (modulo(p+q,2)==0) then
    output = 0
  else
    output = ((q-p)/2)                                               &
         & * exp( log_odd_factorial((p+q-1)/2)                       &
         &      - 0.5_dp*(log_odd_factorial(p)+log_odd_factorial(q)) )
  endif
end procedure

module procedure second_derivative_MonomialState1D
  integer :: p,q
  
  p = bra%power_
  q = ket%power_
  
  if (modulo(p+q,2)==1) then
    output = 0
  else
    output = (((p-q)*(p-q)-1)/(2.0_dp*(p+q-1)) - 1)                  &
         & * exp( log_odd_factorial((p+q)/2)                         &
         &      - 0.5_dp*(log_odd_factorial(p)+log_odd_factorial(q)) )
  endif
end procedure

module procedure read_MonomialState1D
  type(String), allocatable :: line(:)
  
  select type(this); type is(MonomialState1D)
    ! input = '|ui^j>', where i and j are integers.
    ! First remove the '|' and '>',
    !    and then split by the '^', to give ['ui^j'] 
    line = split_line(slice(input,2,len(input)-1),delimiter='^')
    
    ! Construct the state. id=i, power=j.
    this = MonomialState1D( id    = int(line(1)), &
                          & power = int(line(2))  )
  class default
    call err()
  end select
end procedure

module procedure write_MonomialState1D
  select type(this); type is(MonomialState1D)
    output = '|u'//this%id_//'^'//this%power_//'>'
  class default
    call err()
  end select
end procedure

module procedure new_MonomialState1D_String
  call this%read(input)
end procedure
end submodule
