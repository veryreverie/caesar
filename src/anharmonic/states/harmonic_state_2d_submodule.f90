submodule (caesar_harmonic_state_2d_module) caesar_harmonic_state_2d_submodule
  use caesar_states_module
contains

module procedure new_HarmonicState2D
  if (id==paired_id) then
    call print_line(CODE_ERROR//': Mode is its own pair.')
    call err()
  endif
  
  this%id_                = id
  this%paired_id_         = paired_id
  this%occupation_        = occupation
  this%paired_occupation_ = paired_occupation
end procedure

module procedure id_HarmonicState2D
  output = this%id_
end procedure

module procedure paired_id_HarmonicState2D
  output = this%paired_id_
end procedure

module procedure occupation_HarmonicState2D
  output = this%occupation_
end procedure

module procedure paired_occupation_HarmonicState2D
  output = this%paired_occupation_
end procedure

module procedure total_occupation_HarmonicState2D
  output = this%occupation_ + this%paired_occupation_
end procedure

module procedure wavevector_HarmonicState2D
  integer :: mode
  integer :: qpoint
  
  mode   = first(modes%id==this%id_)
  qpoint = first(qpoints%id==modes(mode)%qpoint_id)
  output = qpoints(qpoint)%qpoint * (this%occupation_-this%paired_occupation_)
end procedure

module procedure finite_overlap_HarmonicState2D
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%occupation_
  p_j = bra%paired_occupation_
  q_i = ket%occupation_
  q_j = ket%paired_occupation_
  
  output = p_i==q_i .and. p_j==q_j
end procedure

module procedure inner_product_HarmonicState2D
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%occupation_
  p_j = bra%paired_occupation_
  q_i = ket%occupation_
  q_j = ket%paired_occupation_
  
  if (p_i==q_i .and. p_j==q_j) then
    output = 1
  else
    output = 0
  endif
end procedure

module procedure braket_HarmonicState2D
  output = exp(bra%log_braket( ket,            &
                             & potential,      &
                             & log_2nw,        &
                             & maximum_power,  &
                             & expansion_order ))
end procedure

module procedure log_braket_HarmonicState2D
  real(dp), allocatable, save :: cache(:,:,:,:,:,:)
  logical,  allocatable, save :: cached(:,:,:,:,:,:)
  
  integer :: p_i,p_j,q_i,q_j,n_i,n_j
  
  integer :: k,ialloc
  
  if (.not. allocated(cache)) then
    allocate( cache( 0               :expansion_order,     &
            &        0               :expansion_order,     &
            &        -expansion_order:expansion_order,     &
            &        -expansion_order:expansion_order,     &
            &        0               :maximum_power,       &
            &        0               :maximum_power    ),  &
            & cached( 0               :expansion_order,    &
            &         0               :expansion_order,    &
            &         -expansion_order:expansion_order,    &
            &         -expansion_order:expansion_order,    &
            &         0               :maximum_power,      &
            &         0               :maximum_power    ), &
            & stat=ialloc); call err(ialloc)
    cached = .false.
  endif
  
  p_i = bra%occupation_
  p_j = bra%paired_occupation_
  q_i = ket%occupation_
  q_j = ket%paired_occupation_
  n_i = potential%power
  n_j = potential%paired_power
  
  if (p_i-p_j-n_i+n_j-q_i+q_j/=0 .or. n_i<p_i-q_i .or. n_i<q_j-p_j) then
    output = -1e100_dp
  else
    ! Rather than calculating the result directly, the result is calculated
    !    and cached the first time around, and then the cached result is used
    !    thereafter.
    if (.not.cached(n_i,n_j,p_i-q_i,p_j-q_j,p_i,p_j)) then
      cache(n_i,n_j,p_i-q_i,p_j-q_j,p_i,p_j) =               &
         &   0.5_dp*( log_factorial(p_i)                     &
         &          + log_factorial(q_i)                     &
         &          - log_factorial(p_j)                     &
         &          - log_factorial(q_j) )                   &
         & + log(sum([( exp( log_binomial(n_i,p_i-k)         &
         &                 + log_binomial(n_j,q_i-k)         &
         &                 + log_factorial(n_i+p_j-p_i+k)    &
         &                 - log_factorial(k)             ), &
         &              k=max(0,p_i-n_i,q_i-n_j),            &
         &              min(p_i,q_i)                         )]))
      cached(n_i,n_j,p_i-q_i,p_j-q_j,p_i,p_j) = .true.
    endif
    output = cache(n_i,n_j,p_i-q_i,p_j-q_j,p_i,p_j) &
         & -0.5_dp*(n_i+n_j)*log_2nw
  endif
end procedure

module procedure plus_derivative_HarmonicState2D
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%occupation_
  p_j = bra%paired_occupation_
  q_i = ket%occupation_
  q_j = ket%paired_occupation_
  
  if (p_i-q_i==1 .and. p_j==q_j) then
    output = -sqrt(p_i/4.0_dp)
  elseif (q_j-p_j==1 .and. p_i==q_i) then
    output = sqrt(q_j/4.0_dp)
  else
    output = 0
  endif
end procedure

module procedure minus_derivative_HarmonicState2D
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%occupation_
  p_j = bra%paired_occupation_
  q_i = ket%occupation_
  q_j = ket%paired_occupation_
  
  if (p_j-q_j==1 .and. p_i==q_i) then
    output = -sqrt(p_j/4.0_dp)
  elseif (q_i-p_i==1 .and. p_j==q_j) then
    output = sqrt(q_i/4.0_dp)
  else
    output = 0
  endif
end procedure

module procedure second_derivative_HarmonicState2D
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%occupation_
  p_j = bra%paired_occupation_
  q_i = ket%occupation_
  q_j = ket%paired_occupation_
  
  if (p_i==q_i .and. p_j==q_j) then
    output = -(p_i+p_j+1)/4.0_dp
  elseif (p_i-q_i==1 .and. p_j-q_j==1) then
    output = sqrt(p_i*p_j/16.0_dp)
  elseif (q_i-p_i==1 .and. q_j-p_j==1) then
    output = sqrt(q_i*q_j/16.0_dp)
  else
    output = 0
  endif
end procedure

module procedure read_HarmonicState2D
  type(String), allocatable :: line(:)
  
  select type(this); type is(HarmonicState2D)
    ! input = '|ai^j,ak^l>', where i, j, k and l are integers.
    ! First remove the '|' and '>', and then split by the comma,
    !    to give ['ai^j','ak^l'].
    line = split_line(slice(input,2,len(input)-1),delimiter=',')
    
    ! Split ['ai^j','ak^l'] into ['i','j','k','l'].
    line = [ split_line(slice(line(1),2,len(line(1))),delimiter='^'), &
           & split_line(slice(line(2),2,len(line(2))),delimiter='^')  ]
    
    ! Construct the state. id=i, occupation=j, paired_id=k, paired_occupation=l.
    this = HarmonicState2D( id                = int(line(1)), &
                          & occupation        = int(line(2)), &
                          & paired_id         = int(line(3)), &
                          & paired_occupation = int(line(4))  )
  class default
    call err()
  end select
end procedure

module procedure write_HarmonicState2D
  select type(this); type is(HarmonicState2D)
    output = '|a'                                          // &
           &  this%id_//'^'//this%occupation_              // &
           & ',a'                                          // &
           & this%paired_id_//'^'//this%paired_occupation_ // &
           & '>'
  class default
    call err()
  end select
end procedure

module procedure new_HarmonicState2D_String
  call this%read(input)
end procedure
end submodule
