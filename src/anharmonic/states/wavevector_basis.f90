! ======================================================================
! A basis of monomial and harmonic states for a given subspace
!    at a given wavevector.
! ======================================================================
module wavevector_basis_module
  use common_module
  
  use monomial_state_module
  use braket_module
  use state_conversion_module
  implicit none
  
  private
  
  public :: WavevectorBasis
  public :: size
  public :: generate_subspace_basis_states
  
  type, extends(NoDefaultConstructor) :: WavevectorBasis
    integer                            :: maximum_power
    type(MonomialState),   allocatable :: monomial_states(:)
    type(StateConversion), allocatable :: states_to_basis(:)
    type(StateConversion), allocatable :: basis_to_states(:)
  end type
  
  interface WavevectorBasis
    module procedure new_WavevectorBasis
  end interface
  
  interface size
    module procedure size_WavevectorBasis
  end interface
contains
function new_WavevectorBasis(maximum_power,monomial_states,states_to_basis, &
   & basis_to_states) result(this)
  implicit none
  
  integer,               intent(in) :: maximum_power
  type(MonomialState),   intent(in) :: monomial_states(:)
  type(StateConversion), intent(in) :: states_to_basis(:)
  type(StateConversion), intent(in) :: basis_to_states(:)
  type(WavevectorBasis)                 :: this
  
  if (size(monomial_states)/=size(states_to_basis)) then
    call print_line(CODE_ERROR//': monomial states and states to basis do not &
       & match.')
    call err()
  elseif (size(monomial_states)/=size(basis_to_states)) then
    call print_line(CODE_ERROR//': monomial states and basis to states do not &
       & match.')
    call err()
  endif
  
  this%maximum_power   = maximum_power
  this%monomial_states = monomial_states
  this%states_to_basis = states_to_basis
  this%basis_to_states = basis_to_states
end function

function size_WavevectorBasis(this) result(output)
  implicit none
  
  type(WavevectorBasis), intent(in) :: this
  integer                       :: output
  
  output = size(this%monomial_states)
end function

! ----------------------------------------------------------------------
! Generates states in a given subspace, up to a given power.
! ----------------------------------------------------------------------
function generate_subspace_basis_states(subspace,frequency,modes, &
   & maximum_power) result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency
  type(ComplexMode),        intent(in) :: modes(:)
  integer,                  intent(in) :: maximum_power
  type(WavevectorBasis)                    :: output
  
  ! The state |0>.
  type(ComplexMonomial) :: constant_monomial
  type(MonomialState)   :: ground_state
  type(StateConversion) :: ground_state_conversion
  
  ! Variables for generating single-mode bases.
  type(ComplexMode), allocatable :: subspace_modes(:)
  type(WavevectorBasis)              :: mode_basis
  
  ! Variables for updating output with each mode.
  integer, allocatable :: no_states(:)
  integer, allocatable :: old_to_new(:)
  
  integer  :: id
  real(dp) :: coefficient
  
  type(MonomialState),   allocatable :: monomial_states(:)
  type(StateConversion), allocatable :: states_to_basis(:)
  type(StateConversion), allocatable :: basis_to_states(:)
  
  integer :: i,j,k,l,m,ialloc
  
  ! Construct the state |0>, to seed the generator.
  constant_monomial = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                                     & modes       = [ComplexUnivariate::]    )
  ground_state = MonomialState( subspace_id = subspace%id,      &
                              & frequency   = frequency,        &
                              & state       = constant_monomial )
  ground_state_conversion = StateConversion( ids         =  [1],     &
                                           & coefficients = [1.0_dp] )
  output = WavevectorBasis( maximum_power   = maximum_power,             &
                      & monomial_states = [ground_state],            &
                      & states_to_basis = [ground_state_conversion], &
                      & basis_to_states = [ground_state_conversion]  )
  
  ! List modes in the subspace, but only including one of each conjugate pair.
  subspace_modes = subspace%modes(modes)
  subspace_modes = subspace_modes(filter(          &
     & subspace_modes%paired_id>=subspace_modes%id ))
  
  do i=1,size(subspace_modes)
    ! Generate the basis along a single mode (or a mode and its pair).
    if (subspace_modes(i)%paired_id==subspace_modes(i)%id) then
      mode_basis = mode_basis_1d( subspace,          &
                                & frequency,         &
                                & subspace_modes(i), &
                                & maximum_power      )
    else
      mode_basis = mode_basis_2d( subspace,          &
                                & frequency,         &
                                & subspace_modes(i), &
                                & maximum_power      )
    endif
    
    ! Count the number of states corresponding to each state in the old output.
    allocate( no_states(size(output)),  &
            & old_to_new(size(output)), &
            & stat=ialloc); call err(ialloc)
    old_to_new(1) = 0
    do j=1,size(no_states)
      no_states(j) = count( mode_basis%monomial_states%total_power() &
                       &  + output%monomial_states(j)%total_power()  &
                       & <= maximum_power                            )
      if (j<size(no_states)) then
        old_to_new(j+1) = old_to_new(j) + no_states(j)
      endif
    enddo
    
    ! Generate states.
    allocate(monomial_states(sum(no_states)), stat=ialloc); call err(ialloc)
    do j=1,size(output)
      do k=1,no_states(j)
        monomial_states(old_to_new(j)+k) = output%monomial_states(j) &
                                       & * mode_basis%monomial_states(k)
      enddo
    enddo
    
    ! Generate states-to-basis conversion.
    allocate( states_to_basis(size(monomial_states)), &
            & stat=ialloc); call err(ialloc)
    do j=1,size(output)
      do k=1,no_states(j)
        states_to_basis(old_to_new(j)+k) = StateConversion()
        do l=1,size(output%states_to_basis(j))
          do m=1,size(mode_basis%states_to_basis(k))
            id = old_to_new(output%states_to_basis(j)%id(l)) &
               + mode_basis%states_to_basis(k)%id(m)
            coefficient = output%states_to_basis(j)%coefficient(l) &
                        * mode_basis%states_to_basis(k)%coefficient(m)
            call states_to_basis(old_to_new(j)+k)%add_element( &
                                   & id          = id,         &
                                   & coefficient = coefficient )
          enddo
        enddo
      enddo
    enddo
    
    ! Generate basis-to-states conversion.
    allocate( basis_to_states(size(monomial_states)), &
            & stat=ialloc); call err(ialloc)
    do j=1,size(output)
      do k=1,no_states(j)
        basis_to_states(old_to_new(j)+k) = StateConversion()
        do l=1,size(output%basis_to_states(j))
          do m=1,size(mode_basis%basis_to_states(k))
            id = old_to_new(output%basis_to_states(j)%id(l)) &
               + mode_basis%basis_to_states(k)%id(m)
            coefficient = output%basis_to_states(j)%coefficient(l) &
                        * mode_basis%basis_to_states(k)%coefficient(m)
            call basis_to_states(old_to_new(j)+k)%add_element( &
                                   & id          = id,         &
                                   & coefficient = coefficient )
          enddo
        enddo
      enddo
    enddo
    
    ! Update output.
    output = WavevectorBasis( maximum_power   = maximum_power,   &
                        & monomial_states = monomial_states, &
                        & states_to_basis = states_to_basis, &
                        & basis_to_states = basis_to_states )
    
    deallocate( no_states,       &
              & old_to_new,      &
              & monomial_states, &
              & states_to_basis, &
              & basis_to_states, &
              & stat=ialloc); call err(ialloc)
  enddo
end function

function mode_basis_1d(subspace,frequency,mode,maximum_power) result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency
  type(ComplexMode),        intent(in) :: mode
  integer,                  intent(in) :: maximum_power
  type(WavevectorBasis)                    :: output
  
  integer, allocatable :: ids(:)
  
  type(ComplexUnivariate), allocatable :: univariates(:)
  type(ComplexMonomial),   allocatable :: monomials(:)
  type(MonomialState),     allocatable :: monomial_states(:)
  type(StateConversion),   allocatable :: states_to_basis(:)
  type(StateConversion),   allocatable :: basis_to_states(:)
  
  integer               :: no_states
  integer               :: power
  integer               :: state
  integer               :: term
  real(dp)              :: coefficient
  type(StateConversion) :: old
  
  integer :: i,j,ialloc
  
  no_states = 1+maximum_power
  
  ! Generate the monomial states.
  ! The state monomial_states(id) is |u^{id-1}>, so ID=j+1.
  allocate(ids(0:maximum_power), stat=ialloc); call err(ialloc)
  do power=0,maximum_power
    ids(power) = power+1
  enddo
  univariates = [(ComplexUnivariate(mode,power), power=0, maximum_power)]
  monomials = [( ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp),    &
               &                  modes       = [univariates(i)]         ), &
               & i=1,                                                       &
               & no_states                                                  )]
  monomials(1)%modes = [ComplexUnivariate::]
  monomial_states = [( MonomialState(subspace%id, frequency, monomials(i)), &
                     & i=1,                                                 &
                     & no_states                                            )]
  
  ! A harmonic basis state |i> is defined in terms of monomial states |j> as
  !    |i> = sum h[i,j]|j>.
  ! h[0,0] = 1.
  ! h[i,j] = sqrt((2j-1)/i)       h[i-1,j-1]
  !        - (j+1)/sqrt(i(2j+1))) h[i-1,j+1]
  allocate(states_to_basis(no_states), stat=ialloc); call err(ialloc)
  states_to_basis(1) = StateConversion(ids=[1], coefficients=[1.0_dp])
  do state=2,no_states
    ! Construct h[i,j] from old = h[i-1,{j}].
    old = states_to_basis(state-1)
    
    states_to_basis(state) = StateConversion()
    i  = state-1
    do term=1,size(old)
      ! Add sqrt((2j-1)/i) h[i-1,j-1]
      j = univariates(old%id(term))%power + 1
      coefficient = sqrt((2*j-1.0_dp)/i) * old%coefficient(term)
      call states_to_basis(state)%add_element( id          = ids(j),     &
                                             & coefficient = coefficient )
      
      ! Subtract (j+1)/sqrt(i(2j+1)) h[i-1,j+1]
      j = univariates(old%id(term))%power - 1
      if (j>=0) then
        coefficient = (-(j+1)/sqrt(i*(2*j+1.0_dp))) * old%coefficient(term)
        call states_to_basis(state)%add_element( id          = ids(j),     &
                                               & coefficient = coefficient )
      endif
    enddo
  enddo
  
  ! A monomial state |i> is defined in terms of harmonic basis states |j> as
  !    |i> = sum g[i,j]|j>, where g is the inverse of h above.
  ! g[0,0] = 1.
  ! g[i,j] = sqrt(j/(2i-1))     g[i-1,j-1]
  !        + sqrt((j+1)/(2i-1)) g[i-1,j+1]
  allocate(basis_to_states(no_states), stat=ialloc); call err(ialloc)
  basis_to_states(1) = StateConversion(ids=[1], coefficients=[1.0_dp])
  do state=2,no_states
    ! Construct g[i,j] from g[i-1,{j}].
    old = basis_to_states(state-1)
    
    basis_to_states(state) = StateConversion()
    i  = state-1
    do term=1,size(old)
      ! Add sqrt(j/(2i-1)) g[i-1,j-1]
      j = univariates(old%id(term))%power + 1
      coefficient = sqrt(j/(2*i-1.0_dp)) * old%coefficient(term)
      call basis_to_states(state)%add_element( id          = ids(j),     &
                                             & coefficient = coefficient )
      
      ! Add sqrt((j+1)/(2i-1)) g[i-1,j+1]
      j = univariates(old%id(term))%power - 1
      if (j>=0) then
        coefficient = sqrt((j+1)/(2*i-1.0_dp)) * old%coefficient(term)
        call basis_to_states(state)%add_element( id          = ids(j),     &
                                               & coefficient = coefficient )
      endif
    enddo
  enddo
  
  ! Construct output.
  output = WavevectorBasis( maximum_power   = maximum_power,   &
                      & monomial_states = monomial_states, &
                      & states_to_basis = states_to_basis, &
                      & basis_to_states = basis_to_states  )
end function

function mode_basis_2d(subspace,frequency,mode,maximum_power) result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency
  type(ComplexMode),        intent(in) :: mode
  integer,                  intent(in) :: maximum_power
  type(WavevectorBasis)                    :: output
  
  integer, allocatable :: ids(:,:) ! N.B. ids is zero-indexed.
  
  type(ComplexUnivariate), allocatable :: univariates(:)
  type(ComplexMonomial),   allocatable :: monomials(:)
  type(MonomialState),     allocatable :: monomial_states(:)
  type(StateConversion),   allocatable :: states_to_basis(:)
  type(StateConversion),   allocatable :: basis_to_states(:)
  
  integer               :: no_states
  integer               :: power
  integer               :: state
  integer               :: term
  integer               :: id
  real(dp)              :: coefficient
  type(StateConversion) :: old
  
  integer :: ip,im,jp,jm
  
  integer :: i,j,ialloc
  
  no_states = ((1+maximum_power)*(2+maximum_power))/2
  
  ! Generate the monomial states.
  allocate( ids(0:maximum_power,0:maximum_power), &
          & univariates(no_states),               &
          & stat=ialloc); call err(ialloc)
  id = 0
  do i=0,maximum_power
    do j=0,i
      ip = j
      im = i-j
      id = id+1
      ! Generate the monomial state |ip,im> = |(u+)^(ip) * (u-)^(im)> at
      !    monomial_states(id), and record the (ip,im)->id mapping.
      ids(ip,im) = id
      univariates(id) = ComplexUnivariate(mode, power=ip, paired_power=im)
    enddo
  enddo
  monomials = [( ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp),    &
               &                  modes       = [univariates(i)]         ), &
               & i=1,                                                       &
               & no_states                                                  )]
  monomials(1)%modes = [ComplexUnivariate::]
  monomial_states = [( MonomialState(subspace%id, frequency, monomials(i)), &
                     & i=1,                                                 &
                     & no_states                                            )]
  
  ! A harmonic basis state |ip,im> is defined in terms of monomial states
  !    |jp,jm> as |ip,im> = sum h[ip,im,jp,jm]|jp,jm>.
  ! h[0,0,0,0} = 1.
  ! h[ip,im,jp,jm] = sqrt((jp+jm)/ip)         h[ip-1, im, jp-1, jm  ]
  !                - (jm+1)/sqrt(ip(jp+jm+1)) h[ip-1, im, jp,   jm+1]
  ! also
  ! h[ip,im,jp,jm] = sqrt((jp+jm)/im)          h[ip, im-1, jp,   jm-1]
  !                - (jp+1)/sqrt(im(jp+jm+1))  h[ip, im-1, jp+1, jm  ]
  allocate(states_to_basis(no_states), stat=ialloc); call err(ialloc)
  states_to_basis(1) = StateConversion(ids=[1], coefficients=[1.0_dp])
  do state=2,no_states
    ip = univariates(state)%power
    im = univariates(state)%paired_power
    
    ! If ip>0, construct h[ip,im,jp,jm] from h[ip-1,im,{jp},{jm}].
    ! Otherwise, construct h[ip,im,jp,jm] from h[ip,im-1,{jp},{jm}].
    
    states_to_basis(state) = StateConversion()
    if (ip>0) then
      old = states_to_basis(ids(ip-1,im))
      do term=1,size(old)
        ! Add sqrt((jp+jm)/ip) h[ip-1,im,jp-1,jm]
        jp = univariates(old%id(term))%power + 1
        jm = univariates(old%id(term))%paired_power
        coefficient = sqrt((1.0_dp*(jp+jm))/ip) &
                  & * old%coefficient(term)
        call states_to_basis(state)%add_element( id          = ids(jp,jm), &
                                               & coefficient = coefficient )
        
        ! Subtract (jm+1)/sqrt(ip(jp+jm+1)) h[ip-1,im,jp,jm+1]
        jp = univariates(old%id(term))%power
        jm = univariates(old%id(term))%paired_power - 1
        if (jm>=0) then
          coefficient = ((jm+1)/sqrt(ip*(jp+jm+1.0_dp))) &
                    & * old%coefficient(term)
          call states_to_basis(state)%add_element( id          = ids(jp,jm), &
                                                 & coefficient = coefficient )
        endif
      enddo
    else
      old = states_to_basis(ids(ip,im-1))
      do term=1,size(old)
        ! Add sqrt((jp+jm)/im) h[ip,im-1,jp-1,jm-1]
        jp = univariates(old%id(term))%power
        jm = univariates(old%id(term))%paired_power + 1
        coefficient = sqrt((1.0_dp*(jp+jm))/im) &
                  & * old%coefficient(term)
        call states_to_basis(state)%add_element( id          = ids(jp,jm), &
                                               & coefficient = coefficient )
        
        ! Subtract (jp+1)/sqrt(im(jp+jm+1)) h[ip,im-1,jp+1,jm]
        jp = univariates(old%id(term))%power - 1
        jm = univariates(old%id(term))%paired_power
        if (jp>=0) then
          coefficient = ((jp+1)/sqrt(im*(jp+jm+1.0_dp))) &
                    & * old%coefficient(term)
          call states_to_basis(state)%add_element( id          = ids(jp,jm), &
                                                 & coefficient = coefficient )
        endif
      enddo
    endif
  enddo
  
  ! A monomial state |ip,im> is defined in terms of harmonic basis states
  !    |jp,jm> as |ip,im> = sum g[ip,im,jp,jm]|jp,jm>,
  !    where g is the inverse of h above.
  ! g[0,0,0,0] = 1.
  ! g[ip,im,jp,jm] = sqrt(jp/(ip+im))    g[ip-1,im,jp-1,jm  ]
  !                + sqrt((jm+1)/(ip+im) g[ip-1,im,jp,  jm+1]
  ! also
  ! g[ip,im,jp,jm] = sqrt(jm/(ip+im))    g[ip-1,im,jp,  jm-1]
  !                + sqrt((jp+1)/(ip+im) g[ip-1,im,jp+1,jm  ]
  allocate(basis_to_states(no_states), stat=ialloc); call err(ialloc)
  basis_to_states(1) = StateConversion(ids=[1], coefficients=[1.0_dp])
  do state=2,no_states
    ip = univariates(state)%power
    im = univariates(state)%paired_power
    
    ! If ip>0, construct g[ip,im,jp,jm] from g[ip-1,im,{jp},{jm}].
    ! Otherwise, construct g[ip,im,jp,jm] from f[ip,im-1,{jp},{jm}].
    
    basis_to_states(state) = StateConversion()
    if (ip>0) then
      old = basis_to_states(ids(ip-1,im))
      do term=1,size(old)
        ! Add sqrt(jp/(ip+im)) g[ip-1,im,jp-1,jm]
        jp = univariates(old%id(term))%power + 1
        jm = univariates(old%id(term))%paired_power
        coefficient = sqrt((1.0_dp*jp)/(ip+im)) * old%coefficient(term)
        call basis_to_states(state)%add_element( id          = ids(jp,jm), &
                                               & coefficient = coefficient )
        
        ! Add sqrt((jm-1)/(ip+im)) g[ip-1,im,jp,jm+1]
        jp = univariates(old%id(term))%power
        jm = univariates(old%id(term))%paired_power - 1
        if (jm>=0) then
          coefficient = sqrt((jm+1.0_dp)/(ip+im)) * old%coefficient(term)
          call basis_to_states(state)%add_element( id          = ids(jp,jm), &
                                                 & coefficient = coefficient )
        endif
      enddo
    else
      old = basis_to_states(ids(ip,im-1))
      do term=1,size(old)
        ! Add sqrt(jm/(ip+im)) g[ip,im-1,jp,jm-1]
        jp = univariates(old%id(term))%power
        jm = univariates(old%id(term))%paired_power + 1
        coefficient = sqrt((1.0_dp*jm)/(ip+im)) * old%coefficient(term)
        call basis_to_states(state)%add_element( id          = ids(jp,jm), &
                                               & coefficient = coefficient )
        
        ! Add sqrt((jp-1)/(ip+im)) g[ip,im-1,jp+1,jm]
        jp = univariates(old%id(term))%power - 1
        jm = univariates(old%id(term))%paired_power
        if (jp>=0) then
          coefficient = sqrt((jp+1.0_dp)/(ip+im)) * old%coefficient(term)
          call basis_to_states(state)%add_element( id          = ids(jp,jm), &
                                                 & coefficient = coefficient )
        endif
      enddo
    endif
  enddo
  
  ! Generate output.
  output = WavevectorBasis( maximum_power   = maximum_power,   &
                      & monomial_states = monomial_states, &
                      & states_to_basis = states_to_basis, &
                      & basis_to_states = basis_to_states  )
end function
end module
