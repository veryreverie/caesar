submodule (caesar_harmonic_state_1d_module) caesar_harmonic_state_1d_submodule
  use caesar_states_module
contains

module procedure new_HarmonicState1D
  this%id_         = id
  this%occupation_ = occupation
end procedure

module procedure id_HarmonicState1D
  output = this%id_
end procedure

module procedure occupation_HarmonicState1D
  output = this%occupation_
end procedure

module procedure total_occupation_HarmonicState1D
  output = this%occupation_
end procedure

module procedure wavevector_HarmonicState1D
  integer :: mode
  integer :: qpoint
  
  mode   = first(modes%id==this%id_)
  qpoint = first(qpoints%id==modes(mode)%qpoint_id)
  output = qpoints(qpoint)%qpoint * this%occupation_
end procedure

module procedure finite_overlap_HarmonicState1D
  integer :: p,q
  
  p = bra%occupation_
  q = ket%occupation_
  
  output = p==q
end procedure

module procedure inner_product_HarmonicState1D
  integer :: p,q
  
  p = bra%occupation_
  q = ket%occupation_
  
  if (p==q) then
    output = 1
  else
    output = 0
  endif
end procedure

module procedure braket_HarmonicState1D
  output = exp(bra%log_braket( ket,            &
                             & potential,      &
                             & log_2nw,        &
                             & maximum_power,  &
                             & expansion_order ))
end procedure

module procedure log_braket_HarmonicState1D
  real(dp), allocatable, save :: cache(:,:,:)
  logical,  allocatable, save :: cached(:,:,:)
  
  integer :: min_p,max_p,n,d
  
  integer :: k,ialloc
  
  if (.not. allocated(cache)) then
    allocate( cache(0:expansion_order,0:expansion_order,0:maximum_power),  &
            & cached(0:expansion_order,0:expansion_order,0:maximum_power), &
            & stat=ialloc); call err(ialloc)
    cached = .false.
  endif
  
  min_p = min(bra%occupation_, ket%occupation_)
  max_p = max(bra%occupation_, ket%occupation_)
  n = potential%power
  d = max_p-min_p
  
  if (modulo(min_p+max_p+n,2)==1 .or. n<d) then
    output = -1e100_dp
  else
    ! Rather than calculating the result directly, the result is calculated
    !    and cached the first time around, and then the cached result is used
    !    thereafter.
    if (.not. cached(n,d,min_p)) then
      cache(n,d,min_p) = log_factorial(n)                              &
                     & - log_factorial((n+d)/2)                        &
                     & - ((n-d)/2) * log(2.0_dp)                       &
                     & + 0.5_dp*( log_factorial(max_p)                 &
                     &          - log_factorial(min_p) )               &
                     & + log(sum([( exp( k*log(2.0_dp)                 &
                     &                 + log_binomial((n+d)/2, d+k)    &
                     &                 + log_binomial(min_p, k)     ), &
                     &              k=0,                               &
                     &              min(min_p,(n-d)/2)                 )]))
      cached(n,d,min_p) = .true.
    endif
    output = cache(n,d,min_p)-0.5_dp*n*log_2nw
  endif
end procedure

module procedure first_derivative_HarmonicState1D
  integer :: p,q
  
  p = bra%occupation_
  q = ket%occupation_
  
  if (p==q+1) then
    output = -sqrt(p/4.0_dp)
  elseif (q==p+1) then
    output =  sqrt(q/4.0_dp)
  else
    output = 0
  endif
end procedure

module procedure second_derivative_HarmonicState1D
  integer :: p,q
  
  p = bra%occupation_
  q = ket%occupation_
  
  if (p==q) then
    output = -(p+0.5_dp)/2
  elseif (abs(p-q)==2) then
    output = sqrt(max(p,q)*(max(p,q)-1.0_dp))/4
  else
    output = 0
  endif
end procedure

module procedure read_HarmonicState1D
  type(String), allocatable :: line(:)
  
  select type(this); type is(HarmonicState1D)
    ! input = '|ai^j>', where i and j are integers.
    ! First remove the '|' and '>',
    !    and then split by the '^', to give ['ai^j'] 
    line = split_line(slice(input,2,len(input)-1),delimiter='^')
    
    ! Construct the state. id=i, occupation=j.
    this = HarmonicState1D( id         = int(line(1)), &
                          & occupation = int(line(2))  )
  class default
    call err()
  end select
end procedure

module procedure write_HarmonicState1D
  select type(this); type is(HarmonicState1D)
    output = '|a'//this%id_//'^'//this%occupation_//'>'
  class default
    call err()
  end select
end procedure

module procedure new_HarmonicState1D_String
  call this%read(input)
end procedure
end submodule
