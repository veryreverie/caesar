! ======================================================================
! Conversions between bases. Will eventually become a sparse matrix class.
! ======================================================================
module state_conversion_module
  use common_module
  
  use monomial_state_module
  use braket_module
  implicit none
  
  private
  
  public :: StateConversion
  public :: StatesBasis
  public :: generate_subspace_basis_states
  public :: size
  
  ! Conversions between bases of states.
  type, extends(Stringable) :: StateConversion
    integer,  allocatable :: ids(:)
    real(dp), allocatable :: coefficients(:)
  contains
    ! Add an id and coefficient to the StateConversion.
    procedure, private :: add_element
    ! I/O.
    procedure, public :: read  => read_StateConversion
    procedure, public :: write => write_StateConversion
  end type
  
  interface StateConversion
    module procedure new_StateConversion_null
    module procedure new_StateConversion
    module procedure new_StateConversion_String
  end interface
  
  interface size
    module procedure size_StateConversion
  end interface
  
  ! Output type.
  type, extends(NoDefaultConstructor) :: StatesBasis
    integer                            :: maximum_power
    type(MonomialState),   allocatable :: monomial_states(:)
    type(StateConversion), allocatable :: states_to_basis(:)
    type(StateConversion), allocatable :: basis_to_states(:)
  end type
  
  interface StatesBasis
    module procedure new_StatesBasis
  end interface
  
  interface size
    module procedure size_StatesBasis
  end interface
contains

! Constructors and size functions.
function new_StateConversion_null() result(this)
  implicit none
  
  type(StateConversion) :: this
  
  this%ids = [integer::]
  this%coefficients = [real(dp)::]
end function

function new_StateConversion(ids,coefficients) result(this)
  implicit none
  
  integer,  intent(in)  :: ids(:)
  real(dp), intent(in)  :: coefficients(:)
  type(StateConversion) :: this
  
  if (size(ids)/=size(coefficients)) then
    call print_line(CODE_ERROR//': ids and coefficients do not match.')
    call err()
  endif
  
  this%ids          = ids
  this%coefficients = coefficients
end function

function size_StateConversion(this) result(output)
  implicit none
  
  type(StateConversion), intent(in) :: this
  integer                           :: output
  
  output = size(this%ids)
end function

function new_StatesBasis(maximum_power,monomial_states,states_to_basis, &
   & basis_to_states) result(this)
  implicit none
  
  integer,               intent(in) :: maximum_power
  type(MonomialState),   intent(in) :: monomial_states(:)
  type(StateConversion), intent(in) :: states_to_basis(:)
  type(StateConversion), intent(in) :: basis_to_states(:)
  type(StatesBasis)                 :: this
  
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

function size_StatesBasis(this) result(output)
  implicit none
  
  type(StatesBasis), intent(in) :: this
  integer                       :: output
  
  output = size(this%monomial_states)
end function

! ----------------------------------------------------------------------
! Add and ID and coefficient to a StateConversion.
! ----------------------------------------------------------------------
subroutine add_element(this,id,coefficient)
  implicit none
  
  class(StateConversion), intent(inout) :: this
  integer,                intent(in)    :: id
  real(dp),               intent(in)    :: coefficient
  
  integer :: i
  
  if (.not. (allocated(this%ids) .and. allocated(this%coefficients))) then
    call print_line(CODE_ERROR//': Trying to add and element to a &
       &StateConversion which has not been allocated.')
    call err()
  endif
  
  i = first(this%ids==id, default=0)
  if (i==0) then
    this%ids = [this%ids, id]
    this%coefficients = [this%coefficients, coefficient]
  else
    this%coefficients(i) = this%coefficients(i) + coefficient
  endif
end subroutine

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
  type(StatesBasis)                    :: output
  
  ! The state |0>.
  type(ComplexMonomial) :: constant_monomial
  type(MonomialState)   :: ground_state
  type(StateConversion) :: ground_state_conversion
  
  ! Variables for generating single-mode bases.
  type(ComplexMode), allocatable :: subspace_modes(:)
  type(StatesBasis)              :: mode_basis
  
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
  output = StatesBasis( maximum_power   = maximum_power,             &
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
            id = old_to_new(output%states_to_basis(j)%ids(l)) &
               + mode_basis%states_to_basis(k)%ids(m)
            coefficient = output%states_to_basis(j)%coefficients(l) &
                        * mode_basis%states_to_basis(k)%coefficients(m)
            call states_to_basis(old_to_new(j)+k)%add_element( &
                                   & id          = id,         &
                                   & coefficient = coefficient )
          enddo
        enddo
      enddo
    enddo
    
    ! Generate basis-to-states conversion.
    basis_to_states = calculate_basis_to_states(states_to_basis)
    
    ! Update output.
    output = StatesBasis( maximum_power   = maximum_power,   &
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
  type(StatesBasis)                    :: output
  
  integer, allocatable :: ids(:)
  
  type(ComplexUnivariate), allocatable :: univariates(:)
  type(ComplexMonomial),   allocatable :: monomials(:)
  type(MonomialState),     allocatable :: monomial_states(:)
  type(StateConversion),   allocatable :: states_to_basis(:)
  type(StateConversion),   allocatable :: basis_to_states(:)
  
  integer  :: no_states
  integer  :: power
  integer  :: state
  integer  :: term
  integer  :: id
  real(dp) :: coefficient
  
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
  !    |i> = sum h_{i,j}|j>.
  ! h_{0,0} = 1.
  ! h_{i,j} = sqrt(2j-1)        h_{i-1,j-1}
  !         - (j+1)/sqrt(2j+1)) h_{i-1,j+1}
  allocate(states_to_basis(no_states), stat=ialloc); call err(ialloc)
  states_to_basis(1) = StateConversion(ids=[1], coefficients=[1.0_dp])
  do i=2,no_states
    ! Construct h_{i,j} from {h_{i-1,j}}.
    id = i-1
    states_to_basis(i) = StateConversion()
    do term=1,size(states_to_basis(id))
      ! Add sqrt(2j-1) h_{i-1,j-1}
      j = univariates(states_to_basis(id)%ids(term))%power + 1
      coefficient = sqrt(2*j-1.0_dp) &
                & * states_to_basis(id)%coefficients(term)
      call states_to_basis(i)%add_element( id          = ids(j),     &
                                         & coefficient = coefficient )
      
      ! Subtract (j+1)/sqrt(2j+1) h_{i-1,j+1}.
      j = univariates(states_to_basis(id)%ids(term))%power - 1
      if (j>=0) then
        coefficient = (-(j+1)/sqrt(2*j+1.0_dp)) &
                  & * states_to_basis(id)%coefficients(term)
        call states_to_basis(i)%add_element( id          = ids(j),     &
                                           & coefficient = coefficient )
      endif
    enddo
  enddo
  
  ! Construct basis to states conversion from states to basis conversion.
  basis_to_states = calculate_basis_to_states(states_to_basis)
  
  output = StatesBasis( maximum_power   = maximum_power,   &
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
  type(StatesBasis)                    :: output
  
  integer, allocatable :: ids(:,:) ! N.B. ids is zero-indexed.
  
  type(ComplexUnivariate), allocatable :: univariates(:)
  type(ComplexMonomial),   allocatable :: monomials(:)
  type(MonomialState),     allocatable :: monomial_states(:)
  type(StateConversion),   allocatable :: states_to_basis(:)
  type(StateConversion),   allocatable :: basis_to_states(:)
  
  integer  :: no_states
  integer  :: power
  integer  :: state
  integer  :: term
  integer  :: id
  real(dp) :: coefficient
  
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
  !    |jp,jm> as |ip,im> = sum h_{ip,im,jp,jm}|jp,jm>.
  ! h_{0,0,0,0} = 1.
  ! h_{ip,im,jp,jm} = sqrt(jp+jm)          h_{ip-1, im, jp-1, jm  }
  !                 - (jm+1)/sqrt(jp+jm+1) h_{ip-1, im, jp,   jm+1}
  ! also
  ! h_{ip,im,jp,jm} = sqrt(jp+jm)           h_{ip, im-1, jp,   jm-1}
  !                 - (jp+1)/sqrt(jp+jm+1)  h_{ip, im-1, jp+1, jm  }
  allocate(states_to_basis(no_states), stat=ialloc); call err(ialloc)
  states_to_basis(1) = StateConversion(ids=[1], coefficients=[1.0_dp])
  do state=2,no_states
    ip = univariates(state)%power
    im = univariates(state)%paired_power
    
    ! If ip>0, construct h_{ip,im,jp,jm} from h_{ip-1,im,jp,jm}.
    ! Otherwise, construct h_{ip,im,jp,jm} from h_{ip,im-1,jp,jm}.
    
    states_to_basis(state) = StateConversion()
    if (ip>0) then
      id = ids(ip-1,im)
      do term=1,size(states_to_basis(id))
        ! Add sqrt(jp+jm) h_{ip-1,im,jp-1,jm}
        jp = univariates(states_to_basis(id)%ids(term))%power + 1
        jm = univariates(states_to_basis(id)%ids(term))%paired_power
        coefficient = sqrt(1.0_dp*(jp+jm)) &
                  & * states_to_basis(id)%coefficients(term)
        call states_to_basis(state)%add_element( id          = ids(jp,jm), &
                                               & coefficient = coefficient )
        
        ! Subtract (jm+1)/sqrt(jp+jm+1) h_{ip-1,im,jp,jm+1}
        jp = univariates(states_to_basis(id)%ids(term))%power
        jm = univariates(states_to_basis(id)%ids(term))%paired_power - 1
        if (jm>=0) then
          coefficient = ((jm+1)/sqrt(jp+jm+1.0_dp)) &
                    & * states_to_basis(id)%coefficients(term)
          call states_to_basis(state)%add_element( id          = ids(jp,jm), &
                                                 & coefficient = coefficient )
        endif
      enddo
    else
      id = ids(ip,im-1)
      do term=1,size(states_to_basis(id))
        ! Add sqrt(jp+jm) h_{ip,im-1,jp-1,jm-1}
        jp = univariates(states_to_basis(id)%ids(term))%power
        jm = univariates(states_to_basis(id)%ids(term))%paired_power + 1
        coefficient = sqrt(1.0_dp*(jp+jm)) &
                  & * states_to_basis(id)%coefficients(term)
        call states_to_basis(state)%add_element( id          = ids(jp,jm), &
                                               & coefficient = coefficient )
        
        ! Subtract (jp+1)/sqrt(jp+jm+1) h_{ip,im-1,jp+1,jm}
        jp = univariates(states_to_basis(id)%ids(term))%power - 1
        jm = univariates(states_to_basis(id)%ids(term))%paired_power
        if (jp>=0) then
          coefficient = ((jp+1)/sqrt(jp+jm+1.0_dp)) &
                    & * states_to_basis(id)%coefficients(term)
          call states_to_basis(state)%add_element( id          = ids(jp,jm), &
                                                 & coefficient = coefficient )
        endif
      enddo
    endif
  enddo
  
  ! Construct basis to states conversion from states to basis conversion.
  basis_to_states = calculate_basis_to_states(states_to_basis)
  
  ! Generate output.
  output = StatesBasis( maximum_power   = maximum_power,   &
                      & monomial_states = monomial_states, &
                      & states_to_basis = states_to_basis, &
                      & basis_to_states = basis_to_states  )
end function

! Calculate basis_to_states from states_to_basis.
function calculate_basis_to_states(input) result(output)
  implicit none
  
  type(StateConversion), intent(in)  :: input(:)
  type(StateConversion), allocatable :: output(:)
  
  integer :: i,j
  
  output = [(StateConversion(),i=1,size(input))]
  do i=1,size(input)
    do j=1,size(input(i))
      call output(input(i)%ids(j))%add_element(   &
         & id          = i,                       &
         & coefficient = input(i)%coefficients(j) )
    enddo
  enddo
  
  do i=1,size(output)
    output(i)%coefficients = output(i)%coefficients        &
                         & / ( vec(output(i)%coefficients) &
                         &   * vec(output(i)%coefficients) )
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_StateConversion(this,input)
  implicit none
  
  class(StateConversion), intent(out) :: this
  type(String),           intent(in)  :: input
  
  integer,  allocatable :: ids(:)
  real(dp), allocatable :: coefficients(:)
  
  type(String), allocatable :: states(:)
  type(String), allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(StateConversion)
    states = split_line(input)
    allocate( ids(size(states)),          &
            & coefficients(size(states)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(states)
      line = split_line(states(i), delimiter='|')
      coefficients(i) = dble(line(1))
      ids(i) = int(slice(line(2),1,len(line(2))-1))
    enddo
    
    this = StateConversion(ids, coefficients)
  class default
    call err()
  end select
end subroutine

function write_StateConversion(this) result(output)
  implicit none
  
  class(StateConversion), intent(in) :: this
  type(String)                       :: output
  
  integer :: i
  
  select type(this); type is(StateConversion)
    output = join([( this%coefficients(i)//'|'//this%ids(i)//'>', &
                   & i=1,                                         &
                   & size(this)                                   )])
  class default
    call err()
  end select
end function

function new_StateConversion_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(StateConversion)    :: this
  
  call this%read(input)
end function
end module
