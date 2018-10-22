! ======================================================================
! A state along each complex mode in a degenerate subspace,
!    defined as propotional to a monomial times |0>, where |0> is the
!    harmonic ground state.
! ======================================================================
module monomial_state_module
  use common_module
  
  use subspace_state_module
  implicit none
  
  private
  
  public :: MonomialState
  public :: operator(*)
  public :: generate_subspace_states
  public :: finite_overlap
  public :: braket_MonomialState
  public :: kinetic_energy_MonomialState
  public :: harmonic_potential_energy_MonomialState
  
  ! A state |n_1,n_2,...,n_> = |n_1>|n_2>...|n_>.
  !
  ! If (u_i)*=u_i then
  !    |n_i> = prod_{k=1}^i[ sqrt(2Nw_i/(2k-1)) ] (u_i)^(n_i) |0_i>
  !    |0_i> = sqrt(sqrt(m*w_i/pi)) exp(- 1/2 N w_i (u_i)^2 )
  !
  ! If (u_i)^*=u_j then
  !    |n_i,n_j> = prod_{k=1}^{n_i+n_j}[ sqrt(4Nw_i/k) ] ui^ni uj^nj |0_i,0_j>
  !    |0_i,0_j> = sqrt(2*m*w_i/pi) exp(- N w_i |u_i|^2 )
  ! N.B. w_i=w_j and |u_i|=|u_j|.
  !
  ! In both cases, states are normalised (<n_i|n_i>=1, <n_i,n_j|n_i,n_j>=1),
  !    but are in general not orthogonal (<p_i|q_i>/=0).
  !
  ! w_i is the effective frequency of mode u_i
  !    (in general not the same as the harmonic frequency).
  !
  ! m is the geometric average mass, arising as a result of mass-reduction
  !    of co-ordinates. This cancels in integrals, and so is not stored.
  !
  ! N is the number of primitive cells in the anharmonic supercell.
  !
  ! N.B. the coefficient of state_ is not used. Instead, states are implicitly
  !    normalised such that <state|state>=1. This removes the need to keep
  !    track of the various factors of two, pi, m (the geometric mean of the
  !    atomic masses, arising as a result of mass reduction) and N (the number
  !    of primitive cells in the supercell) , and allows the frequency to be
  !    changed without having to re-calculate coefficients.
  type, extends(SubspaceState) :: MonomialState
    real(dp)                       :: frequency
    type(ComplexMonomial), private :: state_
  contains
    procedure, public :: total_power => total_power_MonomialState
    procedure, public :: wavevector => wavevector_MonomialState
    ! I/O.
    procedure, public :: read  => read_MonomialState
    procedure, public :: write => write_MonomialState
  end type
  
  interface MonomialState
    module procedure new_MonomialState
    module procedure new_MonomialState_Strings
    module procedure new_MonomialState_StringArray
  end interface
  
  interface operator(*)
    module procedure multiply_MonomialState_MonomialState
  end interface
  
  interface finite_overlap
    module procedure finite_overlap_MonomialStates
  end interface
  
  interface braket_MonomialState
    module procedure braket_MonomialStates
    module procedure braket_MonomialStates_ComplexMonomial
  end interface
  
  interface kinetic_energy_MonomialState
    module procedure kinetic_energy_MonomialStates
  end interface
  
  interface harmonic_potential_energy_MonomialState
    module procedure harmonic_potential_energy_MonomialStates
  end interface
contains

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function new_MonomialState(subspace_id,frequency,state) result(this)
  implicit none
  
  integer,               intent(in) :: subspace_id
  real(dp),              intent(in) :: frequency
  type(ComplexMonomial), intent(in) :: state
  type(MonomialState)               :: this
  
  this%subspace_id = subspace_id
  this%frequency   = frequency
  this%state_      = state
end function

! ----------------------------------------------------------------------
! Multiply two monomial states together.
! ----------------------------------------------------------------------
impure elemental function multiply_MonomialState_MonomialState(this,that) &
   & result(output)
  implicit none
  
  type(MonomialState), intent(in) :: this
  type(MonomialState), intent(in) :: that
  type(MonomialState)             :: output
  
  output = this
  output%state_ = output%state_*that%state_
end function

! ----------------------------------------------------------------------
! Generates subspace states up to a given power.
! ----------------------------------------------------------------------
function generate_subspace_states(subspace,frequency,modes,maximum_power) &
   & result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency
  type(ComplexMode),        intent(in) :: modes(:)
  integer,                  intent(in) :: maximum_power
  type(MonomialState), allocatable     :: output(:)
  
  type(ComplexMode), allocatable :: subspace_modes(:)
  type(ComplexMonomial)          :: monomial
  type(MonomialState)            :: state
  
  subspace_modes = subspace%modes(modes)
  monomial = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                            & modes       = [ComplexUnivariate::]    )
  state = MonomialState( subspace_id = subspace%id, &
                       & frequency   = frequency,   &
                       & state       = monomial     )
  output = generate_subspace_states_helper(subspace_modes,maximum_power,state)
end function

recursive function generate_subspace_states_helper(modes,power,state) &
   & result(output)
  implicit none
  
  type(ComplexMode),   intent(in)  :: modes(:)
  integer,             intent(in)  :: power
  type(MonomialState), intent(in)  :: state
  type(MonomialState), allocatable :: output(:)
  
  type(MonomialState) :: output_state
  
  integer :: i
  
  if (size(modes)==0 .or. power==0) then
    output = [state]
  else
    output = [ generate_subspace_states_helper( modes(2:),   &
             &                                  power,       &
             &                                  state      ) ]
    do i=1,power
      output_state = state
      output_state%state_ = output_state%state_ &
                        & * ComplexUnivariate(mode=modes(1), power=i)
      output = [ output,                                         &
               & generate_subspace_states_helper( modes(2:),     &
               &                                  power-i,       &
               &                                  output_state ) ]
    enddo
  endif
end function

! ----------------------------------------------------------------------
! Returns the total power of a given state.
! ----------------------------------------------------------------------
! The total power of the state |product_{q,i} (u_{q,i})^(n_{q,i})> is equal to
!    sum_{q,i}} n_{q,i}.
impure elemental function total_power_MonomialState(this) result(output)
  implicit none
  
  class(MonomialState), intent(in) :: this
  integer                          :: output
  
  output = this%state_%total_power()
end function

! ----------------------------------------------------------------------
! Returns the wavevector of a given state.
! ----------------------------------------------------------------------
! The wavevector of the state |product_{q,i} (u_{q,i})^(n_{q,i})> is equal to
!    sum_{q,i}} n_{q,i}q.
function wavevector_MonomialState(this,modes,qpoints) result(output)
  implicit none
  
  class(MonomialState), intent(in) :: this
  type(ComplexMode),    intent(in) :: modes(:)
  type(QpointData),     intent(in) :: qpoints(:)
  type(QpointData)                 :: output
  
  output = this%state_%wavevector(modes,qpoints)
end function

! ----------------------------------------------------------------------
! Returns whether or not braket(bra,ket) is non-zero.
! ----------------------------------------------------------------------
impure elemental function finite_overlap_MonomialStates(bra,ket) result(output)
  implicit none
  
  type(MonomialState),      intent(in) :: bra
  type(MonomialState),      intent(in) :: ket
  logical                              :: output
  
  logical, allocatable :: bra_mode_integrated(:)
  logical, allocatable :: ket_mode_integrated(:)
  
  type(ComplexUnivariate) :: bra_mode
  type(ComplexUnivariate) :: ket_mode
  
  integer :: id,paired_id
  
  integer :: i,j
  
  ! Check that the bra and the ket cover the same subspace.
  if (bra%subspace_id/=ket%subspace_id) then
    call print_line(ERROR//': bra and ket from different subspaces.')
    call err()
  endif
  
  ! Integrate over modes in subspace.
  bra_mode_integrated = [(.false., i=1, size(bra%state_))]
  ket_mode_integrated = [(.false., i=1, size(ket%state_))]
  
  do while (any(.not.bra_mode_integrated) .or. any(.not.ket_mode_integrated))
    ! Identify an un-integrated mode.
    ! Find its ID and paired ID, and its location in the bra and the ket.
    i = first(.not. bra_mode_integrated, default=0)
    if (i==0) then
      j = first(.not. ket_mode_integrated)
      id = ket%state_%id(j)
      paired_id = ket%state_%paired_id(j)
    else
      id = bra%state_%id(i)
      paired_id = bra%state_%paired_id(i)
      j = first_equivalent(ket%state_%ids(), id, default=0, sorted=.true.)
    endif
    
    ! Multiply the output by the contribution from the mode.
    if (i==0) then
      bra_mode = ComplexUnivariate(id,paired_id,0,0)
    else
      bra_mode = bra%state_%mode(i)
    endif
    
    if (j==0) then
      ket_mode = ComplexUnivariate(id,paired_id,0,0)
    else
      ket_mode = ket%state_%mode(j)
    endif
    
    if (.not. finite_overlap_modes(bra_mode,ket_mode)) then
      output = .false.
      return
    endif
    
    ! Update bra_mode_integrated and ket_mode_integrated.
    if (i/=0) then
      bra_mode_integrated(i) = .true.
    endif
    
    if (j/=0) then
      ket_mode_integrated(j) = .true.
    endif
  enddo
  
  output = .true.
end function

impure elemental function finite_overlap_modes(bra,ket) result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: bra
  type(ComplexUnivariate), intent(in) :: ket
  logical                             :: output
  
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%power
  p_j = bra%paired_power
  q_i = ket%power
  q_j = ket%paired_power
  
  if (bra%id==bra%paired_id) then
    ! <p_i|q_i> = 0                                               if p+q odd.
    !
    !           = prod_{k=1}^{(p_i+q_i)/2} [ 2k-1 ]
    !           / sqrt( prod_{k=1}^{p_i} [ 2k-1 ]
    !                 * prod_{k=1}^{q_i} [ 2k-1 ] )                otherwise.
    !
    if (modulo(p_i+q_i,2)==1) then
      output = .false.
    else
      output = .true.
    endif
  else
    ! <p_i,p_j|q_i,q_j> = 0                          if p_i-p_j-q_i+q_j /= 0.
    !
    !                   = prod_{k=1}^{p_i+p_j+q_i+q_j} [ k ]
    !                   / sqrt( prod_{k=1}^{p_i+p_j} [ k ]
    !                         * prod_{k=1}^{q_i+q_j} [ k ] )       otherwise.
    if (p_i-p_j-q_i+q_j/=0) then
      output = .false.
    else
      output = .true.
    endif
  endif
end function

! ----------------------------------------------------------------------
! Evaluates integrals of the form <bra|ket>.
! ----------------------------------------------------------------------
impure elemental function braket_MonomialStates(bra,ket) result(output)
  implicit none
  
  type(MonomialState), intent(in) :: bra
  type(MonomialState), intent(in) :: ket
  real(dp)                        :: output
  
  logical, allocatable :: bra_mode_integrated(:)
  logical, allocatable :: ket_mode_integrated(:)
  
  type(ComplexUnivariate) :: bra_mode
  type(ComplexUnivariate) :: ket_mode
  
  integer :: id,paired_id
  
  integer :: i,j
  
  ! Check that the bra and the ket cover the same subspace.
  if (bra%subspace_id/=ket%subspace_id) then
    call print_line(ERROR//': bra and ket from different subspaces.')
    call err()
  endif
  
  ! Integrate over modes in subspace.
  bra_mode_integrated = [(.false., i=1, size(bra%state_))]
  ket_mode_integrated = [(.false., i=1, size(ket%state_))]
  
  output = 1.0_dp
  
  do while (any(.not.bra_mode_integrated) .or. any(.not.ket_mode_integrated))
    ! Identify an un-integrated mode.
    ! Find its ID and paired ID, and its location in the bra and the ket.
    i = first(.not. bra_mode_integrated, default=0)
    if (i==0) then
      j = first(.not. ket_mode_integrated)
      id = ket%state_%id(j)
      paired_id = ket%state_%paired_id(j)
    else
      id = bra%state_%id(i)
      paired_id = bra%state_%paired_id(i)
      j = first_equivalent(ket%state_%ids(), id, default=0, sorted=.true.)
    endif
    
    ! Multiply the output by the contribution from the mode.
    if (i==0) then
      bra_mode = ComplexUnivariate(id,paired_id,0,0)
    else
      bra_mode = bra%state_%mode(i)
    endif
    
    if (j==0) then
      ket_mode = ComplexUnivariate(id,paired_id,0,0)
    else
      ket_mode = ket%state_%mode(j)
    endif
    
    output = output * braket_modes(bra_mode, ket_mode)
    
    ! Update bra_mode_integrated and ket_mode_integrated.
    if (i/=0) then
      bra_mode_integrated(i) = .true.
    endif
    
    if (j/=0) then
      ket_mode_integrated(j) = .true.
    endif
  enddo
end function

impure elemental function braket_modes(bra,ket) result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: bra
  type(ComplexUnivariate), intent(in) :: ket
  real(dp)                            :: output
  
  integer :: p_i,p_j,q_i,q_j
  
  integer :: k
  
  p_i = bra%power
  p_j = bra%paired_power
  q_i = ket%power
  q_j = ket%paired_power
  
  if (bra%id==bra%paired_id) then
    ! <p_i|q_i> = 0                                               if p+q odd.
    !
    !           = prod_{k=1}^{(p_i+q_i)/2} [ 2k-1 ]
    !           / sqrt( prod_{k=1}^{p_i} [ 2k-1 ]
    !                 * prod_{k=1}^{q_i} [ 2k-1 ] )                otherwise.
    !
    if (modulo(p_i+q_i,2)==1) then
      output = 0.0_dp
    else
      output = 1.0_dp
      do k=2,(p_i+q_i)/2
        output = output * (2*k-1)
      enddo
      do k=2,p_i
        output = output / sqrt(real(2*k-1,dp))
      enddo
      do k=2,q_i
        output = output / sqrt(real(2*k-1,dp))
      enddo
    endif
  else
    ! <p_i,p_j|q_i,q_j> = 0                          if p_i-p_j-q_i+q_j /= 0.
    !
    !                   = prod_{k=1}^{(p_i+p_j+q_i+q_j)/2} [ k ]
    !                   / sqrt( prod_{k=1}^{p_i+p_j} [ k ]
    !                         * prod_{k=1}^{q_i+q_j} [ k ] )       otherwise.
    if (p_i-p_j-q_i+q_j/=0) then
      output = 0.0_dp
      return
    else
      output = 1.0_dp
      do k=2,(p_i+p_j+q_i+q_j)/2
        output = output * k
      enddo
      do k=2,(p_i+p_j)
        output = output / sqrt(real(k,dp))
      enddo
      do k=2,(q_i+q_j)
        output = output / sqrt(real(k,dp))
      enddo
    endif
  endif
end function

! ----------------------------------------------------------------------
! Evaluates integrals of the form <bra|monomial|ket>.
! ----------------------------------------------------------------------
impure elemental function braket_MonomialStates_ComplexMonomial(bra,ket, &
   & monomial,subspace,supercell) result(output)
  implicit none
  
  type(MonomialState),      intent(in) :: bra
  type(MonomialState),      intent(in) :: ket
  type(ComplexMonomial),    intent(in) :: monomial
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  type(ComplexMonomial)                :: output
  
  real(dp) :: sqrt_two_n_omega
  
  real(dp) :: coefficient
  
  integer :: sum_n
  
  logical, allocatable :: subspace_mode_integrated(:)
  logical, allocatable :: monomial_mode_integrated(:)
  
  type(ComplexUnivariate), allocatable :: monomial_modes(:)
  
  type(ComplexUnivariate) :: bra_mode
  type(ComplexUnivariate) :: ket_mode
  type(ComplexUnivariate) :: monomial_mode
  
  integer :: i_subspace,i_bra,i_ket,i_monomial
  integer :: id,paired_id
  
  type(String), allocatable :: lines(:)
  
  integer :: i
  
  if (bra%subspace_id/=ket%subspace_id) then
    call print_line(ERROR//': bra and ket from different subspaces.')
    call err()
  elseif (bra%subspace_id/=subspace%id) then
    call print_line(ERROR//': bra and subspace do not match.')
    call err()
  endif
  
  ! Integrate over modes in subspace.
  subspace_mode_integrated = [(.false., i=1, size(subspace))]
  monomial_mode_integrated = [(.false., i=1, size(monomial))]
  
  coefficient = 1.0_dp
  do while ((.not. all(subspace_mode_integrated)))
    ! Find the first mode in the subspace which has not yet been integrated.
    i_subspace = first(.not. subspace_mode_integrated)
    id = subspace%mode_ids(i_subspace)
    
    ! Locate this mode in the bra, the ket and the monomial.
    i_bra = first_equivalent( bra%state_%ids(), &
                            & id,               &
                            & default = 0,      &
                            & sorted  = .true.  )
    if (i_bra==0) then
      i_bra = first_equivalent( bra%state_%paired_ids(), &
                              & id,                      &
                              & default = 0,             &
                              & sorted  = .true.         )
    endif
    
    i_ket = first_equivalent( ket%state_%ids(), &
                            & id,               &
                            & default = 0,      &
                            & sorted  = .true.  )
    if (i_ket==0) then
      i_ket = first_equivalent( ket%state_%paired_ids(), &
                              & id,                      &
                              & default = 0,             &
                              & sorted  = .true.         )
    endif
    
    i_monomial = first_equivalent( monomial%ids(),  &
                                 & id,              &
                                 & default = 0,     &
                                 & sorted  = .true. )
    if (i_monomial==0) then
      i_monomial = first_equivalent( monomial%paired_ids(), &
                                   & id,                    &
                                   & default = 0,           &
                                   & sorted  = .true.       )
    endif
    
    ! If the mode does not appear in any, then the contribution is either
    !    <0|(u_i)^0|0> = 1 or <0,0|(u_i)^0(u_j)^0|0,0>=1,
    !    so the mode can be neglected.
    if (i_bra==0 .and. i_ket==0 .and. i_monomial==0) then
      subspace_mode_integrated(i_subspace) = .true.
      cycle
    endif
    
    ! Use mode locations to find p_i, p_j, n_i, n_j, q_i, q_j,
    !    and id and paired_id.
    ! N.B. id<=paired_id, so id may be changed from its value above.
    if (i_bra/=0) then
      id = bra%state_%id(i_bra)
      paired_id = bra%state_%paired_id(i_bra)
    elseif (i_ket/=0) then
      id = ket%state_%id(i_ket)
      paired_id = ket%state_%paired_id(i_ket)
    elseif (i_monomial/=0) then
      id = monomial%id(i_monomial)
      paired_id = monomial%paired_id(i_monomial)
    else
      call err()
    endif
    
    if (i_bra==0) then
      bra_mode = ComplexUnivariate(id,paired_id,0,0)
    else
      bra_mode = bra%state_%mode(i_bra)
    endif
    
    if (i_ket==0) then
      ket_mode = ComplexUnivariate(id,paired_id,0,0)
    else
      ket_mode = ket%state_%mode(i_ket)
    endif
    
    if (i_monomial==0) then
      monomial_mode = ComplexUnivariate(id,paired_id,0,0)
    else
      monomial_mode = monomial%mode(i_monomial)
    endif
    
    coefficient = coefficient &
              & * braket_modes_potential(bra_mode,ket_mode,monomial_mode)
    
    ! Update subspace_mode_integrated and monomial_mode_integrated.
    subspace_mode_integrated(first(subspace%mode_ids==id)) = .true.
    subspace_mode_integrated(first(subspace%mode_ids==paired_id)) = .true.
    
    if (i_monomial/=0) then
      monomial_mode_integrated(i_monomial) = .true.
    endif
  enddo
  
  ! Calculate sqrt(4Nw_i).
  sqrt_two_n_omega = sqrt(2.0_dp * supercell%sc_size * bra%frequency)
  
  monomial_modes = monomial%modes()
  
  sum_n = sum(monomial_modes(filter(monomial_mode_integrated))%total_power())
  
  coefficient = coefficient / sqrt_two_n_omega**sum_n
  
  ! The output is the un-integrated modes of the monomial, multiplied by
  !    the new coefficient.
  output = ComplexMonomial(                                                &
     & coefficient = monomial%coefficient*coefficient,                     &
     & modes       = monomial_modes(filter(.not.monomial_mode_integrated)) )
end function

impure elemental function braket_modes_potential(bra,ket,potential) &
   & result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: bra
  type(ComplexUnivariate), intent(in) :: ket
  type(ComplexUnivariate), intent(in) :: potential
  real(dp)                            :: output
  
  integer :: p_i,p_j,q_i,q_j,n_i,n_j
  
  integer :: k
  
  p_i = bra%power
  p_j = bra%paired_power
  q_i = ket%power
  q_j = ket%paired_power
  n_i = potential%power
  n_j = potential%paired_power
  
  if (bra%id==bra%paired_id) then
    ! <p_i|(u_i)^(n_i)|q_i> = 0                                 if p+n+q odd.
    !
    !                       = 1/sqrt(2Nw_i)^{n_i}
    !                       * prod_{k=1}^{(p_i+n_i+q_i)/2} [ 2k-1 ]
    !                       / sqrt( prod_{k=1}^{p_i} [ 2k-1 ]
    !                             * prod_{k=1}^{q_i} [ 2k-1 ] )    otherwise.
    !
    if (modulo(p_i+n_i+q_i,2)==1) then
      output = 0
    else
      output = 1
      do k=2,(p_i+n_i+q_i)/2
        output = output * (2*k-1)
      enddo
      do k=2,p_i
        output = output / sqrt(real(2*k-1,dp))
      enddo
      do k=2,q_i
        output = output / sqrt(real(2*k-1,dp))
      enddo
    endif
  else
    ! <p_i,p_j|(u_i)^(n_i) (u_j)^(n_j)|q_i,q_j> =
    !    = 0                                 if p_i-p_j-n_i+n_j-q_i+q_j /= 0.
    !
    !    = 1/sqrt(2Nw_i)^{n_i+n_j}
    !    * prod_{k=1}^{(p_i+p_j+n_i+n_j+q_i+q_j)/2} [ k ]
    !    / sqrt( prod_{k=1}^{p_i+p_j} [ k ]
    !          * prod_{k=1}^{q_i+q_j} [ k ] )                      otherwise.
    if (p_i-p_j-n_i+n_j-q_i+q_j/=0) then
      output = 0
    else
      output = 1
      do k=2,(p_i+p_j+n_i+n_j+q_i+q_j)/2
        output = output * k
      enddo
      do k=2,(p_i+p_j)
        output = output / sqrt(real(k,dp))
      enddo
      do k=2,(q_i+q_j)
        output = output / sqrt(real(k,dp))
      enddo
    endif
  endif
end function

! ----------------------------------------------------------------------
! Evaluates <bra|T|ket>, where T is the kinetic energy operator.
! Gives the result per primitive cell.
! ----------------------------------------------------------------------
function kinetic_energy_MonomialStates(bra,ket,subspace,supercell) &
   & result(output)
  implicit none
  
  type(MonomialState),      intent(in) :: bra
  type(MonomialState),      intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  real(dp)                             :: output
  
  real(dp) :: prefactor
  real(dp) :: frequency
  
  type(ComplexUnivariate) :: bra_mode
  type(ComplexUnivariate) :: ket_mode
  
  logical, allocatable :: bra_mode_integrated(:)
  logical, allocatable :: ket_mode_integrated(:)
  logical, allocatable :: subspace_mode_integrated(:)
  
  integer :: i_bra,i_ket
  
  integer :: id,paired_id
  
  integer :: i
  
  ! Check that the bra and the ket cover the same subspace.
  if (bra%subspace_id/=ket%subspace_id) then
    call print_line(ERROR//': bra and ket from different subspaces.')
    call err()
  endif
  
  ! Calculate the prefactor, s.t. <bra|T|ket> = prefactor * <bra|ket>.
  prefactor = 0.0_dp
  
  frequency = bra%frequency
  
  bra_mode_integrated = [(.false., i=1, size(bra%state_))]
  ket_mode_integrated = [(.false., i=1, size(ket%state_))]
  subspace_mode_integrated = [(.false., i=1, size(subspace))]
  
  do while (any(.not. bra_mode_integrated) .or. any(.not. ket_mode_integrated))
    ! Identify an un-integrated mode.
    ! Find i_bra and i_ket, the locations of the mode within the bra and
    !    the ket respectively.
    ! Find p_i,p_j,q_i and q_j.
    i_bra = first(.not. bra_mode_integrated, default=0)
    if (i_bra==0) then
      i_ket = first(.not. ket_mode_integrated)
      id = ket%state_%id(i_ket)
      paired_id = ket%state_%paired_id(i_ket)
    else
      id = bra%state_%id(i_bra)
      paired_id = bra%state_%paired_id(i_bra)
      i_ket = first_equivalent(ket%state_%ids(), id, default=0, sorted=.true.)
    endif
    
    if (i_bra==0) then
      bra_mode = ComplexUnivariate(id,paired_id,0,0)
    else
      bra_mode = bra%state_%mode(i_bra)
    endif
    
    if (i_ket==0) then
      ket_mode = ComplexUnivariate(id,paired_id,0,0)
    else
      ket_mode = ket%state_%mode(i_ket)
    endif
    
    ! Calculate the contribution to the prefactor from the mode.
    prefactor = prefactor &
            & + kinetic_energy_prefactor(bra_mode,ket_mode,frequency)
    
    ! Update bra_mode_integrated ket_mode_integrated and
    !    subspace_mode_integrated.
    if (i_bra/=0) then
      bra_mode_integrated(i_bra) = .true.
    endif
    
    if (i_ket/=0) then
      ket_mode_integrated(i_ket) = .true.
    endif
    
    subspace_mode_integrated(first(subspace%mode_ids==id)) = .true.
    subspace_mode_integrated(first(subspace%mode_ids==paired_id)) = .true.
  enddo
  
  ! Account for modes in subspace which do not appear in bra or ket.
  prefactor = prefactor &
          & + 0.25_dp * frequency * count(.not. subspace_mode_integrated)
  
  ! Calculate output, and normalise by the number of primitive cells.
  output = prefactor * braket_MonomialState(bra,ket) / supercell%sc_size
end function

impure elemental function kinetic_energy_prefactor(bra,ket,frequency) &
   & result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: bra
  type(ComplexUnivariate), intent(in) :: ket
  real(dp),                intent(in) :: frequency
  real(dp)                            :: output
  
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%power
  p_j = bra%paired_power
  q_i = ket%power
  q_j = ket%paired_power
  
  if (bra%id==bra%paired_id) then
    ! <p_i|T|q_i> = 1/2 w_i ( (1-p_i-q_i)/2 - 2p_iq_i/(p_i+q_i+1) ) <p_i|q_i>
    if (modulo(p_i+q_i,2)==0) then
      output = 0.5_dp * frequency &
           & * ( (1-p_i-q_i)/2.0_dp + (2.0_dp*p_i*q_i)/(p_i+q_i-1) )
    else
      output = 0.0_dp
    endif
  else
    ! <p_i,p_j|T|q_i,q_j> = w_i ( 1/2
    !                           + ((p_i-p_j)^2-(p_i-q_i)^2)
    !                           / (p_i+p_j+q_i+q_j)         )
    !                     * <p_i,p_j|q_i,q_j>
    !
    ! N.B. if p_i=p_j=q_i=q_j=0 then the second term is zero.
    if (p_i+q_j==p_j+q_i) then
      output = frequency/2
      if (p_i+p_j+q_i+q_j/=0) then
        output = output                        &
             & + frequency                     &
             & * ((p_i-p_j)**2 - (p_i-q_i)**2) &
             & / (1.0_dp*(p_i+p_j+q_i+q_j))
      endif
    else
      output = 0.0_dp
    endif
  endif
end function

! ----------------------------------------------------------------------
! Evaluates <bra|V|ket>, where V is the harmonic potential energy operator.
! Gives the result per primitive cell.
! ----------------------------------------------------------------------
function harmonic_potential_energy_MonomialStates(bra,ket,subspace,supercell) &
   & result(output)
  implicit none
  
  type(MonomialState),      intent(in) :: bra
  type(MonomialState),      intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  real(dp)                             :: output
  
  real(dp) :: prefactor
  real(dp) :: frequency
  
  type(ComplexUnivariate) :: bra_mode
  type(ComplexUnivariate) :: ket_mode
  
  logical, allocatable :: bra_mode_integrated(:)
  logical, allocatable :: ket_mode_integrated(:)
  logical, allocatable :: subspace_mode_integrated(:)
  
  integer :: i_bra,i_ket
  
  integer :: id,paired_id
  
  integer :: i
  
  ! Check that the bra and the ket cover the same subspace.
  if (bra%subspace_id/=ket%subspace_id) then
    call print_line(ERROR//': bra and ket from different subspaces.')
    call err()
  endif
  
  ! Calculate the prefactor, s.t. <bra|T|ket> = prefactor * <bra|ket>.
  prefactor = 0.0_dp
  
  frequency = bra%frequency
  
  bra_mode_integrated = [(.false., i=1, size(bra%state_))]
  ket_mode_integrated = [(.false., i=1, size(ket%state_))]
  subspace_mode_integrated = [(.false., i=1, size(subspace))]
  
  do while (any(.not. bra_mode_integrated) .or. any(.not. ket_mode_integrated))
    ! Identify an un-integrated mode.
    ! Find i_bra and i_ket, the locations of the mode within the bra and
    !    the ket respectively.
    ! Find p_i,p_j,q_i and q_j.
    i_bra = first(.not. bra_mode_integrated, default=0)
    if (i_bra==0) then
      i_ket = first(.not. ket_mode_integrated)
      id = ket%state_%id(i_ket)
      paired_id = ket%state_%paired_id(i_ket)
    else
      id = bra%state_%id(i_bra)
      paired_id = bra%state_%paired_id(i_bra)
      i_ket = first_equivalent(ket%state_%ids(), id, default=0, sorted=.true.)
    endif
    
    if (i_bra==0) then
      bra_mode = ComplexUnivariate(id,paired_id,0,0)
    else
      bra_mode = bra%state_%mode(i_bra)
    endif
    
    if (i_ket==0) then
      ket_mode = ComplexUnivariate(id,paired_id,0,0)
    else
      ket_mode = ket%state_%mode(i_ket)
    endif
    
    ! Calculate the contribution to the prefactor from the mode.
    prefactor = prefactor                                      &
            & + harmonic_potential_energy_prefactor( bra_mode, &
            &                                        ket_mode, &
            &                                        frequency )
    
    ! Update bra_mode_integrated ket_mode_integrated and
    !    subspace_mode_integrated.
    if (i_bra/=0) then
      bra_mode_integrated(i_bra) = .true.
    endif
    
    if (i_ket/=0) then
      ket_mode_integrated(i_ket) = .true.
    endif
    
    subspace_mode_integrated(first(subspace%mode_ids==id)) = .true.
    subspace_mode_integrated(first(subspace%mode_ids==paired_id)) = .true.
  enddo
  
  ! Account for modes in subspace which do not appear in bra or ket.
  prefactor = prefactor &
          & + 0.25_dp * frequency * count(.not. subspace_mode_integrated)
  
  ! Calculate output, and normalise by the number of primitive cells.
  output = prefactor * braket_MonomialState(bra,ket) / supercell%sc_size
end function

impure elemental function harmonic_potential_energy_prefactor(bra,ket, &
   & frequency) result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: bra
  type(ComplexUnivariate), intent(in) :: ket
  real(dp),                intent(in) :: frequency
  real(dp)                            :: output
  
  integer :: p_i,p_j,q_i,q_j
  
  p_i = bra%power
  p_j = bra%paired_power
  q_i = ket%power
  q_j = ket%paired_power
  
  if (bra%id==bra%paired_id) then
    ! <p_i|V|q_i> = 1/4 w_i ( 1 + p_i + q_i ) <p_i|q_i>
    if (modulo(p_i+q_i,2)==0) then
      output = 0.25_dp * frequency * (1 + p_i + q_i)
    else
      output = 0.0_dp
    endif
  else
    ! <p_i,p_j|V|q_i,q_j> = 1/4 w_i ( 2 + (p_i+p_j+q_i+q_j) ) <p_i,p_j|q_i,q_j>
    if (p_i+q_j==p_j+q_i) then
      output = 0.25_dp * frequency * (2+p_i+p_j+q_i+q_j)
    else
      output = 0.0_dp
    endif
  endif
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_MonomialState(this,input)
  implicit none
  
  class(MonomialState), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  integer               :: subspace_id
  real(dp)              :: frequency
  type(ComplexMonomial) :: state
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(MonomialState)
    line = split_line(input(1))
    subspace_id = int(line(3))
    
    line = split_line(input(2))
    frequency = dble(line(3))
    
    if (input(4)=='|0>') then
      state = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                             & modes       = [ComplexUnivariate::]    )
    else
      state = ComplexMonomial(1.0_dp//'*'//slice(input(4),2,len(input(4))-1))
    endif
    
    this = MonomialState(subspace_id,frequency,state)
  class default
    call err()
  end select
end subroutine

function write_MonomialState(this) result(output)
  implicit none
  
  class(MonomialState), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  type(String) :: state_string
  
  select type(this); type is(MonomialState)
    if (this%state_%total_power()==0) then
      state_string = '|0>'
    else
      state_string = '|'//join(this%state_%modes(),delimiter='*')//'>'
    endif
    output = [ 'Subspace  : '//this%subspace_id, &
             & 'Frequency : '//this%frequency,  &
             & str('State'),                  &
             & state_string                   ]
  class default
    call err()
  end select
end function

function new_MonomialState_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(MonomialState)      :: this
  
  call this%read(input)
end function

impure elemental function new_MonomialState_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(MonomialState)           :: this
  
  this = MonomialState(str(input))
end function
end module
