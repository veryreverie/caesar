! ======================================================================
! A product of states along each complex mode in a degenerate subspace.
! ======================================================================
module subspace_state_module
  use common_module
  implicit none
  
  private
  
  public :: SubspaceState
  public :: generate_subspace_states
  public :: braket
  public :: kinetic_energy
  
  ! A state |n_1,n_2,...,n_> = |n_1>|n_2>...|n_>.
  !
  ! If (u_i)*=u_i then
  !    |n_i> = prod_{k=1}^i[ sqrt(2Nw_i/(2k-1)) ] |0_i>
  !    |0_i> = sqrt(sqrt(m*w_i/pi)) exp(- 1/2 N w_i (u_i)^2 )
  !
  ! If (u_i)^*=u_j then
  !    |n_i,n_j> = prod_{k=1}^{n_i+n_j}[ sqrt(2Nw_i/k) ] |0_i,0_j>
  !    |0_i,0_j> = sqrt(m*w_i/pi) exp(- N w_i |u_i|^2 )
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
  ! --------------------
  ! <bra|ket>:
  ! --------------------
  !
  ! (u_i)^* = u_i:
  ! <p_i|q_i> = 0                                                   if p+q odd.
  !
  !           = prod_{k=1}^{(p_i+q_i)/2} [ 2k-1 ]
  !           / sqrt( prod_{k=1}^{p_i} [ 2k-1 ]
  !                 * prod_{k=1}^{q_i} [ 2k-1 ] )                    otherwise.
  !
  ! (u_i)^* /= u_i:
  ! <p_i,p_j|q_i,q_j> = 0                              if p_i-p_j-q_i+q_j /= 0.
  !
  !                   = prod_{k=1}^{p_i+p_j+q_i+q_j} [ k ]
  !                   / sqrt( prod_{k=1}^{p_i+p_j} [ k ]
  !                         * prod_{k=1}^{q_i+q_j} [ k ] )           otherwise.
  !
  ! --------------------
  ! Potential energy:
  ! --------------------
  !
  ! (u_i)^* = u_i:
  ! <p_i|(u_i)^(n_i)|q_i> = 0                                     if p+n+q odd.
  !
  !                       = 1/sqrt(2Nw_i)^{n_i}
  !                       * prod_{k=1}^{(p_i+n_i+q_i)/2} [ 2k-1 ]
  !                       / sqrt( prod_{k=1}^{p_i} [ 2k-1 ]
  !                             * prod_{k=1}^{q_i} [ 2k-1 ] )        otherwise.
  !
  ! (u_i)^* /= u_i:
  ! <p_i,p_j|(u_i)^(n_i) (u_j)^(n_j)|q_i,q_j> =
  !    = 0                                     if p_i-p_j-n_i+n_j-q_i+q_j /= 0.
  !
  !    = 1/sqrt(2Nw_i)^{n_i+n_j}
  !    * prod_{k=1}^{p_i+p_j+n_i+n_j+q_i+q_j} [ k ]
  !    / sqrt( prod_{k=1}^{p_i+p_j} [ k ]
  !          * prod_{k=1}^{q_i+q_j} [ k ] )                          otherwise.
  !
  ! --------------------
  ! Kinetic energy:
  ! --------------------
  !
  ! (u_i)^* = u_i:
  ! <p_i|T|q_i> = 1/2 w_i  ( (1-p_i-q_i)/2 - 2p_iq_i/(p_i+q_i+1) ) <p_i|q_i>
  !
  ! (u_i)^* /= u_i:
  ! <p_i,p_j|T|q_i,q_j> = w_i ( (p_i+p_j+q_i+q_j+2)/4 
  !                           - 2(p_ip_j+q_iq_j)/(p_i+p_j+q_i+q_j) )
  !                     * <p_i,p_j|q_i,q_j>
  !
  ! N.B. (p_ip_j+q_iq_j)/(p_i+p_j+q_i+q_j) is zero if p_i=p_j=q_i=q_j=0.
  ! The denominator is zero in this case, so care should be taken.
  
  type, extends(Stringsable) :: SubspaceState
    integer               :: subspace_id
    real(dp)              :: frequency
    type(ComplexMonomial) :: state
  contains
    procedure, public :: read  => read_SubspaceState
    procedure, public :: write => write_SubspaceState
  end type
  
  interface SubspaceState
    module procedure new_SubspaceState
    module procedure new_SubspaceState_Strings
    module procedure new_SubspaceState_StringArray
  end interface
  
  interface braket
    module procedure braket_SubspaceStates
    module procedure braket_SubspaceStates_ComplexMonomial
    module procedure braket_SubspaceStates_ComplexPolynomial
  end interface
  
  interface kinetic_energy
    module procedure kinetic_energy_SubspaceStates
  end interface
contains

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function new_SubspaceState(subspace_id,frequency,state) result(this)
  implicit none
  
  integer,               intent(in) :: subspace_id
  real(dp),              intent(in) :: frequency
  type(ComplexMonomial), intent(in) :: state
  type(SubspaceState)               :: this
  
  this%subspace_id = subspace_id
  this%frequency   = frequency
  this%state       = state
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
  type(SubspaceState), allocatable     :: output(:)
  
  type(ComplexMode), allocatable :: subspace_modes(:)
  type(ComplexMonomial)          :: monomial
  type(SubspaceState)            :: state
  
  subspace_modes = subspace%modes(modes)
  monomial = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                            & modes       = [ComplexUnivariate::]    )
  state = SubspaceState( subspace_id = subspace%id, &
                       & frequency   = frequency,   &
                       & state       = monomial     )
  output = generate_subspace_states_helper(subspace_modes,maximum_power,state)
end function

recursive function generate_subspace_states_helper(modes,power,state) &
   & result(output)
  implicit none
  
  type(ComplexMode),   intent(in)  :: modes(:)
  integer,             intent(in)  :: power
  type(SubspaceState), intent(in)  :: state
  type(SubspaceState), allocatable :: output(:)
  
  type(SubspaceState) :: output_state
  
  integer :: i
  
  if (size(modes)==0 .or. power==0) then
    output = [state]
  else
    output = [SubspaceState::]
    output_state = state
    do i=0,power
      output = [ output,                                         &
               & generate_subspace_states_helper( modes(2:),     &
               &                                  power-i,       &
               &                                  output_state ) ]
      output_state%state = output_state%state &
                       & * ComplexUnivariate(mode=modes(1), power=i)
    enddo
  endif
end function

! ----------------------------------------------------------------------
! Evaluates integrals of the form <bra|ket>.
! ----------------------------------------------------------------------
! (u_i)^* = u_i:
! <p_i|q_i> = 0                                                   if p+q odd.
!
!           = prod_{k=1}^{(p_i+q_i)/2} [ 2k-1 ]
!           / sqrt( prod_{k=1}^{p_i} [ 2k-1 ]
!                 * prod_{k=1}^{q_i} [ 2k-1 ] )                    otherwise.
!
! (u_i)^* /= u_i:
! <p_i,p_j|q_i,q_j> = 0                              if p_i-p_j-q_i+q_j /= 0.
!
!                   = prod_{k=1}^{p_i+p_j+q_i+q_j} [ k ]
!                   / sqrt( prod_{k=1}^{p_i+p_j} [ k ]
!                         * prod_{k=1}^{q_i+q_j} [ k ] )           otherwise.
impure elemental function braket_SubspaceStates(bra,ket) result(output)
  implicit none
  
  type(SubspaceState),      intent(in) :: bra
  type(SubspaceState),      intent(in) :: ket
  complex(dp)                          :: output
  
  logical, allocatable :: bra_mode_integrated(:)
  logical, allocatable :: ket_mode_integrated(:)
  
  integer :: i_bra,j_bra,i_ket,j_ket
  integer :: i,j
  integer :: p_i,p_j,q_i,q_j
  
  integer :: k
  
  ! Check that the bra and the ket cover the same subspace.
  if (bra%subspace_id/=ket%subspace_id) then
    call print_line(ERROR//': bra and ket from different subspaces.')
    call err()
  endif
  
  ! Integrate over modes in subspace.
  bra_mode_integrated = [(.false., i=1, size(bra%state))]
  ket_mode_integrated = [(.false., i=1, size(ket%state))]
  
  output = 1.0_dp
  do while (any(.not.bra_mode_integrated) .or. any(.not.ket_mode_integrated))
    ! Identify an un-integrated mode with ID i.
    ! Find j s.t. (u_i)^* = u_j.
    ! Find i_bra,j_bra,i_ket and j_ket, the locations of modes i and j
    !    respectively within the bra and ket respectively.
    ! Find p_i,p_j,q_i and q_j, the powers of modes i and j respectively
    !    within the bra and ket respectively.
    i_bra = first(.not. bra_mode_integrated, default=0)
    if (i_bra==0) then
      j_bra = 0
      p_i = 0
      p_j = 0
      
      i_ket = first(.not. ket_mode_integrated)
      i = ket%state%modes(i_ket)%id
      j = ket%state%modes(i_ket)%paired_id
      q_i = ket%state%modes(i_ket)%power
      if (i==j) then
        j_ket = i_ket
        q_j   = q_i
      else
        j_ket = first(ket%state%modes%id==j, default=0)
        if (j_ket==0) then
          q_j = 0
        else
          q_j = ket%state%modes(j_ket)%power
        endif
      endif
    else
      i = bra%state%modes(i_bra)%id
      j = bra%state%modes(i_bra)%paired_id
      p_i = bra%state%modes(i)%power
      
      i_ket = first(ket%state%modes%id==i, default=0)
      if (i_ket==0) then
        q_i = 0
      else
        q_i = ket%state%modes(i_ket)%power
      endif
      
      if (i==j) then
        j_bra = i_bra
        j_ket = i_ket
        p_j   = p_i
        q_j   = q_i
      else
        j_bra = first(bra%state%modes%id==j, default=0)
        if (j_bra==0) then
          p_j = 0
        else
          p_j = bra%state%modes(j_bra)%power
        endif
        
        j_ket = first(ket%state%modes%id==j, default=0)
        if (j_ket==0) then
          q_j = 0
        else
          q_j = ket%state%modes(j_ket)%power
        endif
      endif
    endif
    
    ! Integrate over mode i and j (or just mod i if i==j).
    if (i==j) then
      ! <p_i|q_i> = 0                                               if p+q odd.
      !
      !           = prod_{k=1}^{(p_i+q_i)/2} [ 2k-1 ]
      !           / sqrt( prod_{k=1}^{p_i} [ 2k-1 ]
      !                 * prod_{k=1}^{q_i} [ 2k-1 ] )                otherwise.
      !
      if (modulo(p_i+q_i,1)==1) then
        output = 0.0_dp
        return
      endif
      
      do k=2,(p_i+q_i)/2
        output = output * (2*k-1)
      enddo
      
      do k=2,p_i
        output = output / sqrt(real(2*k-1,dp))
      enddo
      
      do k=2,q_i
        output = output / sqrt(real(2*k-1,dp))
      enddo
    else
      ! <p_i,p_j|q_i,q_j> = 0                          if p_i-p_j-q_i+q_j /= 0.
      !
      !                   = prod_{k=1}^{p_i+p_j+q_i+q_j} [ k ]
      !                   / sqrt( prod_{k=1}^{p_i+p_j} [ k ]
      !                         * prod_{k=1}^{q_i+q_j} [ k ] )       otherwise.
      if (p_i-p_j-q_i+q_j/=0) then
        output = 0.0_dp
        return
      endif
      
      do k=2,(p_i+p_j+q_i+q_j)
        output = output * k
      enddo
      
      do k=2,(p_i+p_j)
        output = output / sqrt(real(k,dp))
      enddo
      
      do k=2,(q_i+q_j)
        output = output / sqrt(real(k,dp))
      enddo
    endif
    
    ! Update bra_mode_integrated and ket_mode_integrated.
    if (i_bra/=0) then
      bra_mode_integrated(i_bra) = .true.
    endif
    
    if (j_bra/=0) then
      bra_mode_integrated(j_bra) = .true.
    endif
    
    if (i_ket/=0) then
      ket_mode_integrated(i_ket) = .true.
    endif
    
    if (j_ket/=0) then
      ket_mode_integrated(j_ket) = .true.
    endif
  enddo
end function

! ----------------------------------------------------------------------
! Evaluates integrals of the form <bra|monomial|ket>.
! ----------------------------------------------------------------------
! (u_i)^* = u_i:
! <p_i|(u_i)^(n_i)|q_i> = 0                                     if p+n+q odd.
!
!                       = 1/sqrt(2Nw_i)^{n_i}
!                       * prod_{k=1}^{(p_i+n_i+q_i)/2} [ 2k-1 ]
!                       / sqrt( prod_{k=1}^{p_i} [ 2k-1 ]
!                             * prod_{k=1}^{q_i} [ 2k-1 ] )        otherwise.
!
! (u_i)^* /= u_i:
! <p_i,p_j|(u_i)^(n_i) (u_j)^(n_j)|q_i,q_j> =
!    = 0                                     if p_i-p_j-n_i+n_j-q_i+q_j /= 0.
!
!    = 1/sqrt(2Nw_i)^{n_i+n_j}
!    * prod_{k=1}^{p_i+p_j+n_i+n_j+q_i+q_j} [ k ]
!    / sqrt( prod_{k=1}^{p_i+p_j} [ k ]
!          * prod_{k=1}^{q_i+q_j} [ k ] )                          otherwise.
impure elemental function braket_SubspaceStates_ComplexMonomial(bra,ket, &
   & monomial,subspace,supercell) result(output)
  implicit none
  
  type(SubspaceState),      intent(in) :: bra
  type(SubspaceState),      intent(in) :: ket
  type(ComplexMonomial),    intent(in) :: monomial
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  type(ComplexMonomial)                :: output
  
  real(dp) :: sqrt_two_n_omega
  
  real(dp) :: coefficient
  
  logical, allocatable :: subspace_mode_integrated(:)
  logical, allocatable :: monomial_mode_integrated(:)
  
  integer :: i_bra,j_bra,i_ket,j_ket
  integer :: i_subspace,j_subspace,i_monomial,j_monomial
  integer :: i,j
  integer :: p_i,p_j,q_i,q_j,n_i,n_j
  
  integer :: k
  
  if (bra%subspace_id/=ket%subspace_id) then
    call print_line(ERROR//': bra and ket from different subspaces.')
    call err()
  elseif (bra%subspace_id/=subspace%id) then
    call print_line(ERROR//': bra and subspace do not match.')
    call err()
  endif
  
  ! Calculate sqrt(2Nw_i).
  sqrt_two_n_omega = sqrt(2.0_dp * supercell%sc_size * bra%frequency)
  
  ! Integrate over modes in subspace.
  subspace_mode_integrated = [(.false., i=1, size(subspace))]
  monomial_mode_integrated = [(.false., i=1, size(monomial))]
  
  coefficient = 1.0_dp
  do while ((.not. all(subspace_mode_integrated)))
    ! Find the first mode in the subspace which has not yet been integrated.
    i_subspace = first(.not. subspace_mode_integrated)
    i = subspace%mode_ids(i_subspace)
    
    ! Locate this mode (mode i) in the bra, the ket and the monomial.
    i_bra = first(bra%state%modes%id==i, default=0)
    i_ket = first(ket%state%modes%id==i, default=0)
    i_monomial = first(monomial%modes%id==i, default=0)
    
    ! If the mode does not appear in any, then the contribution is either
    !    <0|(u_i)^0|0> = 0, or <0,p_j|(u_i)^0(u_j)^(n_j)|0,q_j>, in which
    !    case the contribution will be handled when mode j is integrated.
    if (i_bra==0 .and. i_ket==0 .and. i_monomial==0) then
      subspace_mode_integrated(i_subspace) = .true.
      cycle
    endif
    
    ! Use mode locations to find p_i, n_i, q_i and j.
    if (i_bra==0) then
      p_i = 0
    else
      p_i = bra%state%modes(i_bra)%power
      j = bra%state%modes(i_bra)%paired_id
    endif
    
    if (i_ket==0) then
      q_i = 0
    else
      q_i = ket%state%modes(i_ket)%power
      j = ket%state%modes(i_ket)%paired_id
    endif
    
    if (i_monomial==0) then
      n_i = 0
    else
      n_i = monomial%modes(i_monomial)%power
      j = monomial%modes(i_monomial)%paired_id
    endif
    
    ! Locate mode j in the subspace, the bra, the ket and the monomial,
    !    and use these locations to find p_j, n_j and q_j.
    if (i==j) then
      j_subspace = i_subspace
      j_bra      = i_bra
      j_ket      = i_ket
      j_monomial = i_monomial
      p_j        = p_i
      n_j        = n_i
      q_j        = q_i
    else
      j_subspace = first(subspace%mode_ids==j)
      j_bra = first(bra%state%modes%id==j, default=0)
      j_ket = first(ket%state%modes%id==j, default=0)
      j_monomial = first(monomial%modes%id==j, default=0)
      
      if (j_bra==0) then
        p_j = 0
      else
        p_j = bra%state%modes(j_bra)%power
      endif
      
      if (j_ket==0) then
        q_j = 0
      else
        q_j = ket%state%modes(j_ket)%power
      endif
      
      if (j_monomial==0) then
        n_j = 0
      else
        n_j = monomial%modes(j_monomial)%power
      endif
    endif
    
    ! Integrate over mode i and j (or just mode i if i==j).
    if (i==j) then
      ! <p_i|(u_i)^(n_i)|q_i> = 0                                 if p+n+q odd.
      !
      !                       = 1/sqrt(2Nw_i)^{n_i}
      !                       * prod_{k=1}^{(p_i+n_i+q_i)/2} [ 2k-1 ]
      !                       / sqrt( prod_{k=1}^{p_i} [ 2k-1 ]
      !                             * prod_{k=1}^{q_i} [ 2k-1 ] )    otherwise.
      !
      if (modulo(p_i+n_i+q_i,1)==1) then
        output = ComplexMonomial( coefficient=cmplx(0.0_dp,0.0_dp,dp), &
                                & modes=[ComplexUnivariate::]          )
        return
      endif
      
      coefficient = coefficient / sqrt_two_n_omega**n_i
      
      do k=2,(p_i+n_i+q_i)/2
        coefficient = coefficient * (2*k-1)
      enddo
      
      do k=2,p_i
        coefficient = coefficient / sqrt(real(2*k-1,dp))
      enddo
      
      do k=2,q_i
        coefficient = coefficient / sqrt(real(2*k-1,dp))
      enddo
    else
      ! <p_i,p_j|(u_i)^(n_i) (u_j)^(n_j)|q_i,q_j> =
      !    = 0                                 if p_i-p_j-n_i+n_j-q_i+q_j /= 0.
      !
      !    = 1/sqrt(2Nw_i)^{n_i+n_j}
      !    * prod_{k=1}^{p_i+p_j+n_i+n_j+q_i+q_j} [ k ]
      !    / sqrt( prod_{k=1}^{p_i+p_j} [ k ]
      !          * prod_{k=1}^{q_i+q_j} [ k ] )                      otherwise.
      if (p_i-p_j-n_i+n_j-q_i+q_j/=0) then
        output = ComplexMonomial( coefficient=cmplx(0.0_dp,0.0_dp,dp), &
                                & modes=[ComplexUnivariate::]          )
        return
      endif
      
      coefficient = coefficient / sqrt_two_n_omega**(n_i+n_j)
      
      do k=2,(p_i+p_j+n_i+n_j+q_i+q_j)/2
        coefficient = coefficient * k
      enddo
      
      do k=2,(p_i+p_j)
        coefficient = coefficient / sqrt(real(k,dp))
      enddo
      
      do k=2,(q_i+q_j)
        coefficient = coefficient / sqrt(real(k,dp))
      enddo
    endif
    
    ! Update subspace_mode_integrated and monomial_mode_integrated.
    subspace_mode_integrated(i_subspace) = .true.
    subspace_mode_integrated(j_subspace) = .true.
    
    if (i_monomial/=0) then
      monomial_mode_integrated(i_monomial) = .true.
    endif
    
    if (j_monomial/=0) then
      monomial_mode_integrated(j_monomial) = .true.
    endif
  enddo
  
  ! The output is the un-integrated modes of the monomial, multiplied by
  !    the new coefficient.
  output = ComplexMonomial(                                                &
     & coefficient = monomial%coefficient*coefficient,                     &
     & modes       = monomial%modes(filter(.not.monomial_mode_integrated)) )
end function

! ----------------------------------------------------------------------
! Evaluates integrals of the form <bra|polynomial|ket>.
! ----------------------------------------------------------------------
impure elemental function braket_SubspaceStates_ComplexPolynomial(bra,ket, &
   & polynomial,subspace,supercell) result(output)
  implicit none
  
  type(SubspaceState),      intent(in) :: bra
  type(SubspaceState),      intent(in) :: ket
  type(ComplexPolynomial),  intent(in) :: polynomial
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  type(ComplexPolynomial)              :: output
  
  output = ComplexPolynomial(braket( bra,              &
                                   & ket,              &
                                   & polynomial%terms, &
                                   & subspace,         &
                                   & supercell         ))
end function

! ----------------------------------------------------------------------
! Evaluates <bra|T|ket>, where T is the kinetic energy operator.
! Gives the result per primitive cell.
! ----------------------------------------------------------------------
! (u_i)^* = u_i:
! <p_i|T|q_i> = 1/2 w_i  ( (1-p_i-q_i)/2 - 2p_iq_i/(p_i+q_i+1) ) <p_i|q_i>
!
! (u_i)^* /= u_i:
! <p_i,p_j|T|q_i,q_j> = w_i ( (p_i+p_j+q_i+q_j+2)/4 
!                           - 2(p_ip_j+q_iq_j)/(p_i+p_j+q_i+q_j) )
!                     * <p_i,p_j|q_i,q_j>
!
! N.B. (p_ip_j+q_iq_j)/(p_i+p_j+q_i+q_j) is zero if p_i=p_j=q_i=q_j=0.
! The denominator is zero in this case, so care should be taken.
function kinetic_energy_SubspaceStates(bra,ket,subspace,supercell) &
   & result(output)
  implicit none
  
  type(SubspaceState),      intent(in) :: bra
  type(SubspaceState),      intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  real(dp)                             :: output
  
  real(dp) :: prefactor
  real(dp) :: frequency
  
  logical, allocatable :: bra_mode_integrated(:)
  logical, allocatable :: ket_mode_integrated(:)
  logical, allocatable :: subspace_mode_integrated(:)
  
  integer :: i_bra,j_bra,i_ket,j_ket
  integer :: p_i,p_j,q_i,q_j
  
  integer :: i,j
  
  ! Check that the bra and the ket cover the same subspace.
  if (bra%subspace_id/=ket%subspace_id) then
    call print_line(ERROR//': bra and ket from different subspaces.')
    call err()
  endif
  
  ! Calculate the prefactor, s.t. <bra|T|ket> = prefactor * <bra|ket>.
  prefactor = 0.0_dp
  
  frequency = bra%frequency
  
  bra_mode_integrated = [(.false., i=1, size(bra%state))]
  ket_mode_integrated = [(.false., i=1, size(ket%state))]
  subspace_mode_integrated = [(.false., i=1, size(subspace))]
  
  do while (any(.not. bra_mode_integrated) .or. any(.not. ket_mode_integrated))
    ! Identify an un-integrated mode with ID i.
    ! Find j s.t. (u_i)^* = u_j.
    ! Find i_bra,j_bra,i_ket and j_ket, the locations of modes i and j
    !    respectively within the bra and ket respectively.
    ! Find p_i,p_j,q_i and q_j, the powers of modes i and j respectively
    !    within the bra and ket respectively.
    i_bra = first(.not. bra_mode_integrated, default=0)
    if (i_bra==0) then
      j_bra = 0
      p_i = 0
      p_j = 0
      
      i_ket = first(.not. ket_mode_integrated)
      i = ket%state%modes(i_ket)%id
      j = ket%state%modes(i_ket)%paired_id
      q_i = ket%state%modes(i_ket)%power
      if (i==j) then
        j_ket = i_ket
        q_j   = q_i
      else
        j_ket = first(ket%state%modes%id==j, default=0)
        if (j_ket==0) then
          q_j = 0
        else
          q_j = ket%state%modes(j_ket)%power
        endif
      endif
    else
      i = bra%state%modes(i_bra)%id
      j = bra%state%modes(i_bra)%paired_id
      p_i = bra%state%modes(i)%power
      
      i_ket = first(ket%state%modes%id==i, default=0)
      if (i_ket==0) then
        q_i = 0
      else
        q_i = ket%state%modes(i_ket)%power
      endif
      
      if (i==j) then
        j_bra = i_bra
        j_ket = i_ket
        p_j   = p_i
        q_j   = q_i
      else
        j_bra = first(bra%state%modes%id==j, default=0)
        if (j_bra==0) then
          p_j = 0
        else
          p_j = bra%state%modes(j_bra)%power
        endif
        
        j_ket = first(ket%state%modes%id==j, default=0)
        if (j_ket==0) then
          q_j = 0
        else
          q_j = ket%state%modes(j_ket)%power
        endif
      endif
    endif
    
    ! Calculate the contribution to the prefactor from modes i and j,
    !    or just from mode i if i=j.
    if (i==j) then
      ! <p_i|T|q_i> = 1/2 w_i ( (1-p_i-q_i)/2 - 2p_iq_i/(p_i+q_i+1) ) <p_i|q_i>
      prefactor = prefactor &
              & + 0.5_dp    &
              & * frequency &
              & * ( (1-p_i-q_i)/2.0_dp - (2.0_dp*p_i*q_i)/(p_i+q_i+1) )
    else
      ! <p_i,p_j|T|q_i,q_j> = w_i ( (p_i+p_j+q_i+q_j+2)/4 
      !                           - 2(p_ip_j+q_iq_j)/(p_i+p_j+q_i+q_j) )
      !                     * <p_i,p_j|q_i,q_j>
      !
      ! N.B. (p_ip_j+q_iq_j)/(p_i+p_j+q_i+q_j) is zero if p_i=p_j=q_i=q_j=0.
      ! The denominator is zero in this case, so care should be taken.
      if (p_i+p_j+q_i+q_j==0) then
        prefactor = prefactor + frequency / 2
      else
        prefactor = prefactor                                    &
                & + frequency                                    &
                & * ( (p_i+p_j+q_i+q_j+2)/4.0_dp                 &
                &   - 2.0_dp*(p_i*p_j+q_i*q_j)/(p_i+p_j+q_i+q_j) )
      endif
    endif
    
    ! Update bra_mode_integrated ket_mode_integrated and
    !    subspace_mode_integrated.
    if (i_bra/=0) then
      bra_mode_integrated(i_bra) = .true.
    endif
    
    if (j_bra/=0) then
      bra_mode_integrated(j_bra) = .true.
    endif
    
    if (i_ket/=0) then
      ket_mode_integrated(i_ket) = .true.
    endif
    
    if (j_ket/=0) then
      ket_mode_integrated(j_ket) = .true.
    endif
    
    subspace_mode_integrated(first(subspace%mode_ids==i)) = .true.
    subspace_mode_integrated(first(subspace%mode_ids==j)) = .true.
  enddo
  
  ! Account for modes in subspace which do not appear in bra or ket.
  prefactor = prefactor &
          & + 0.25_dp * frequency * count(.not. subspace_mode_integrated)
  
  ! Calculate output, and normalise by the number of primitive cells.
  output = prefactor * braket(bra,ket) / supercell%sc_size
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_SubspaceState(this,input)
  implicit none
  
  class(SubspaceState), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  integer               :: subspace_id
  real(dp)              :: frequency
  type(ComplexMonomial) :: state
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(SubspaceState)
    line = split_line(input(1))
    subspace_id = int(line(2))
    
    line = split_line(input(2))
    frequency = dble(line(2))
    
    state = ComplexMonomial(input(4))
    
    this = SubspaceState(subspace_id,frequency,state)
  end select
end subroutine

function write_SubspaceState(this) result(output)
  implicit none
  
  class(SubspaceState), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  select type(this); type is(SubspaceState)
    output = [ 'Subspace '//this%subspace_id, &
             & 'Frequency '//this%frequency,  &
             & str('State'),                  &
             & str(this%state)                ]
  end select
end function

function new_SubspaceState_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(SubspaceState)      :: this
  
  call this%read(input)
end function

impure elemental function new_SubspaceState_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SubspaceState)           :: this
  
  this = SubspaceState(str(input))
end function
end module
