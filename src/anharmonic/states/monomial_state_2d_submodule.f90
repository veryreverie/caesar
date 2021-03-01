submodule (caesar_monomial_state_2d_module) caesar_monomial_state_2d_submodule
  use caesar_states_module
contains

module procedure new_MonomialState2D
  if (id==paired_id) then
    call print_line(CODE_ERROR//': Mode is its own pair.')
    call err()
  endif
  
  this%id_           = id
  this%paired_id_    = paired_id
  this%power_        = power
  this%paired_power_ = paired_power
end procedure

module procedure id_MonomialState2D
  output = this%id_
end procedure

module procedure paired_id_MonomialState2D
  output = this%paired_id_
end procedure

module procedure power_MonomialState2D
  output = this%power_
end procedure

module procedure paired_power_MonomialState2D
  output = this%paired_power_
end procedure

module procedure total_power_MonomialState2D
  output = this%power_ + this%paired_power_
end procedure

module procedure wavevector_MonomialState2D
  integer :: mode
  integer :: qpoint
  
  mode   = first(modes%id==this%id_)
  qpoint = first(qpoints%id==modes(mode)%qpoint_id)
  output = qpoints(qpoint)%qpoint * (this%power_-this%paired_power_)
end procedure

module procedure finite_overlap_MonomialState2D
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%power_
  p_j = bra%paired_power_
  q_i = ket%power_
  q_j = ket%paired_power_
  
  output = p_i-p_j-q_i+q_j==0
end procedure

module procedure inner_product_MonomialState2D
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%power_
  p_j = bra%paired_power_
  q_i = ket%power_
  q_j = ket%paired_power_
  
  if (p_i-p_j-q_i+q_j/=0) then
    output = 0.0_dp
  else
    output = exp( log_factorial((p_i+p_j+q_i+q_j)/2)                     &
              & - 0.5_dp*(log_factorial(p_i+p_j)+log_factorial(q_i+q_j)) )
  endif
end procedure

module procedure braket_MonomialState2D
  integer :: p_i,p_j,q_i,q_j,n_i,n_j
  
  p_i = bra%power_
  p_j = bra%paired_power_
  q_i = ket%power_
  q_j = ket%paired_power_
  n_i = potential%power
  n_j = potential%paired_power
  
  if (p_i-p_j-n_i+n_j-q_i+q_j/=0) then
    output = 0
  else
    output = exp( log_factorial((p_i+p_j+n_i+n_j+q_i+q_j)/2)             &
              & - 0.5_dp*(log_factorial(p_i+p_j)+log_factorial(q_i+q_j)) &
              & - 0.5_dp*(n_i+n_j)*log_2nw                               )
  endif
end procedure

module procedure plus_derivative_MonomialState2D
  integer :: p_i,p_j,q_i,q_j
  integer :: p,q
  
  p_i = bra%power_
  p_j = bra%paired_power_
  q_i = ket%power_
  q_j = ket%paired_power_
  
  p = p_i+p_j
  q = q_i+q_j
  
  if (q_i-p_i-q_j+p_j/=1) then
    output = 0
  else
    output = (q_i - sqrt((q+1)/(p+1.0_dp))*(p+q+1)/4.0_dp)   &
         & * exp( log_factorial(p+q-1)                       &
         &      - 0.5_dp*(log_factorial(p)+log_factorial(q)) )
  endif
end procedure

module procedure minus_derivative_MonomialState2D
  integer :: p_i,p_j,q_i,q_j
  integer :: p,q
  
  p_i = bra%power_
  p_j = bra%paired_power_
  q_i = ket%power_
  q_j = ket%paired_power_
  
  p = p_i+p_j
  q = q_i+q_j
  
  if (q_j-p_j-q_i+p_i/=1) then
    output = 0
  else
    output = (q_j - sqrt((q+1)/(p+1.0_dp))*(p+q+1)/4.0_dp)   &
         & * exp( log_factorial(p+q-1)                       &
         &      - 0.5_dp*(log_factorial(p)+log_factorial(q)) )
  endif
end procedure

module procedure second_derivative_MonomialState2D
  integer :: p_i,p_j,q_i,q_j
  integer :: p,q
  
  p_i = bra%power_
  p_j = bra%paired_power_
  q_i = ket%power_
  q_j = ket%paired_power_
  
  p = p_i+p_j
  q = q_i+q_j
  
  if (p_i-p_j-q_i+q_j/=0) then
    output = 0
  else
    output = ((p_i-q_j)*(p_j-q_i)/(2.0_dp*(p+q))-0.25_dp)    &
         & * exp( log_factorial((p+q)/2)                     &
         &      - 0.5_dp*(log_factorial(p)+log_factorial(q)) )
  endif
end procedure

module procedure read_MonomialState2D
  type(String), allocatable :: line(:)
  
  select type(this); type is(MonomialState2D)
    ! input = '|ui^j,uk^l>', where i, j, k and l are integers.
    ! First remove the '|' and '>', and then split by the comma,
    !    to give ['ui^j','uk^l'].
    line = split_line(slice(input,2,len(input)-1),delimiter=',')
    
    ! Split ['ui^j','uk^l'] into ['i','j','k','l'].
    line = [ split_line(slice(line(1),2,len(line(1))),delimiter='^'), &
           & split_line(slice(line(2),2,len(line(2))),delimiter='^')  ]
    
    ! Construct the state. id=i, power=j, paired_id=k, paired_power=l.
    this = MonomialState2D( id           = int(line(1)), &
                          & power        = int(line(2)), &
                          & paired_id    = int(line(3)), &
                          & paired_power = int(line(4))  )
  class default
    call err()
  end select
end procedure

module procedure write_MonomialState2D
  select type(this); type is(MonomialState2D)
    output = '|u'                                     // &
           &  this%id_//'^'//this%power_              // &
           & ',u'                                     // &
           & this%paired_id_//'^'//this%paired_power_ // &
           & '>'
  class default
    call err()
  end select
end procedure

module procedure new_MonomialState2D_String
  call this%read(input)
end procedure
end submodule
