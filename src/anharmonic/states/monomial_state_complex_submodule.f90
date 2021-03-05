submodule (caesar_monomial_state_complex_module) caesar_monomial_state_complex_submodule
  use caesar_states_module
contains

module procedure new_MonomialStateComplex
  this%supercell_size = supercell_size
  this%frequency      = frequency
  this%modes_         = modes
  
  this%log_2nw_ = log(2*this%supercell_size*this%frequency)
end procedure

module procedure new_MonomialStateComplex_SubspaceState
  select type(input); type is(MonomialStateComplex)
    this = input
  type is(SubspaceStatePointer)
    ! WORKAROUND: ifort doesn't recognise the interface to this module procedure
    !    from within this module procedure, so the full name is used instead.
    this = new_MonomialStateComplex_SubspaceState(input%state())
  class default
    call err()
  end select
end procedure

module procedure monomial_state_complex_pointer
  select type(input); type is(MonomialStateComplex)
    this => input
  type is(SubspaceStatePointer)
    this => monomial_state_complex_pointer(input%state_pointer())
  class default
    call err()
  end select
end procedure

module procedure representation_MonomialStateComplex
  output = 'monomial complex'
end procedure

module procedure mode_ids_MonomialStateComplex
  output = this%modes_%id()
end procedure

module procedure paired_mode_ids_MonomialStateComplex
  output = this%modes_%paired_id()
end procedure

module procedure occupation_MonomialStateComplex
  output = sum(this%modes_%total_power())
end procedure

module procedure wavevector_MonomialStateComplex
  integer :: i
  
  output = sum([( this%modes_(i)%wavevector(modes,qpoints), &
                & i=1,                                      &
                & size(this%modes_)                         )])
end procedure

module procedure wavefunction_MonomialStateComplex
  real(dp) :: log_2nw
  real(dp) :: coefficient
  
  integer :: i
  
  ! |n_i,n_j> = sqrt((2Nw)^(n_i+n_j)/(n_i+n_j)!) (u_i)^(n_i) (u_j)^(n_j) |0>,
  
  log_2nw = log(2*supercell%sc_size*frequency)
  coefficient = product(exp(0.5_dp*(                                       &
     &   (this%modes_%power()+this%modes_%paired_power())*log_2nw          &
     & - log_odd_factorial(this%modes_%power()+this%modes_%paired_power()) )))
  
  if (size(this%modes_)==0) then
    output = str(coefficient)//'|0,0>'
  else
    output = coefficient                                          // &
       & '*'                                                      // &
       & join( [( '(u'//this%modes_(i)%id()          //              &
       &           '^'//this%modes_(i)%power()       //              &
       &           'u'//this%modes_(i)%paired_id()   //              &
       &           '^'//this%modes_(i)%paired_power()//')',          &
       &          i=1,                                               &
       &          size(this%modes_)                         )],      &
       &       delimiter='*'                                    ) // &
       & '|0,0>'
  endif
end procedure

module procedure finite_overlap_MonomialStateComplexs
  output = all(bra%modes_%finite_overlap(ket%modes_))
end procedure

module procedure inner_product_MonomialStateComplex
  type(MonomialStateComplex), pointer :: monomial_ket
  
  if (present(ket)) then
    monomial_ket => monomial_state_complex_pointer(ket)
    output = product(this%modes_%inner_product(monomial_ket%modes_))
  else
    ! Modes are normalised, so <p|p>=1.
    output = 1
  endif
end procedure

module procedure integrate_MonomialStateComplex
  type(MonomialStateComplex), pointer :: monomial_ket
  
  if (present(ket)) then
    monomial_ket => monomial_state_complex_pointer(ket)
    output = product(this%modes_%braket( monomial_ket%modes_, &
                                       & monomial%modes,      &
                                       & this%log_2nw_        ))
  else
    output = product(this%modes_%braket( this%modes_,    &
                                       & monomial%modes, &
                                       & this%log_2nw_   ))
  endif
end procedure

module procedure kinetic_energy_MonomialStateComplex
  type(MonomialStateComplex), pointer :: monomial_ket
  
  type(MonomialState2D), allocatable :: bra_modes(:)
  type(MonomialState2D), allocatable :: ket_modes(:)
  real(dp),              allocatable :: overlap(:)
  
  ! |p> = product_i |p_i>
  ! The kinetic energy is given by T = -(1/N) sum_i d^2/d(u_i^2).
  ! <p_i|d^2/d(u_i)^2|q_i> is calculated up to a factor of 2Nw, so
  !    <p|T|q> = -2w * sum_i <p_i|%second_derivative(|q_i>) * (<p'|q'>),
  !    where |p'> is |p> excluding |p_i>, so |p>=|p_i>|p'>.
  
  bra_modes = this%modes_
  
  if (present(ket)) then
    monomial_ket => monomial_state_complex_pointer(ket)
    ket_modes = monomial_ket%modes_
    
    if (all(bra_modes%finite_overlap(ket_modes))) then
      ! All <p_i|q_i>/=0, so all <p'|q'>=<p|q>/<p_i|q_i>,
      !    so <p|T|q> = -2w<p|q>
      !               * sum_i <p_i|%second_derivative(|q_i>) / <p_i|q_i>
      overlap = bra_modes%inner_product(ket_modes)
      output = -2*this%frequency*product(overlap) &
           & * sum(bra_modes%second_derivative(ket_modes)/overlap)
    else
      ! At least one <p_i|q_i>=0, and for monomial states,
      !    if <p_i|q_i>=0 then <p_i|T|q_i>=0 as well.
      output = 0
    endif
  else
    ! |p>=|q>, so <p'|q'>=1, and so
    !    <p|T|q> = -2w*sum_i <p_i|%second_derivative().
    output = -2*this%frequency*sum(bra_modes%second_derivative(bra_modes))
  endif
end procedure

module procedure harmonic_potential_energy_MonomialStateComplex
  type(MonomialStateComplex), pointer :: monomial_ket
  
  type(MonomialState2D),   allocatable :: bra_modes(:)
  type(MonomialState2D),   allocatable :: ket_modes(:)
  type(ComplexUnivariate), allocatable :: harmonic_potential(:)
  real(dp),                allocatable :: overlap(:)
  
  ! |p> = product_i |p_i>
  ! The harmonic potential energy is given by V = Nw^2 sum_i (u_i^2).
  ! <p_i|(u_i)^2|q_i> is calculated up to a factor of 2Nw, so
  !    <p|V|q> = (w/2) * sum_i <p_i|%braket(|q_i>,V_i) * (<p'|q'>),
  !    where |p'> is |p> excluding |p_i>, so |p>=|p_i>|p'>.
  
  bra_modes = this%modes_
  
  harmonic_potential = ComplexUnivariate(    &
     & id           = bra_modes%id(),        &
     & paired_id    = bra_modes%paired_id(), &
     & power        = 1,                     &
     & paired_power = 1                      )
  
  if (present(ket)) then
    monomial_ket => monomial_state_complex_pointer(ket)
    ket_modes = monomial_ket%modes_
    
    if (all(bra_modes%finite_overlap(ket_modes))) then
      ! All <p_i|q_i>/=0, so all <p'|q'>=<p|q>/<p_i|q_i>,
      !    so <p|T|q> = (w^2/2)<p|q>
      !               * sum_i <p_i|%braket(|q_i>,V_i) / <p_i|q_i>
      overlap = bra_modes%inner_product(ket_modes)
      output = (this%frequency**2/2)*product(overlap)       &
           & * sum( bra_modes%braket( ket_modes,            &
           &                          harmonic_potential,   &
           &                          this%log_2nw_       ) &
           &      / overlap                                 )
    else
      ! At least one <p_i|q_i>=0, and for monomial states,
      !    if <p_i|q_i>=0 then <p_i|V_i|q_i>=0 as well.
      output = 0
    endif
  else
    ! |p>=|q>, so <p'|q'>=1, and so
    !    <p|V|q> = (w^2/2)*sum_i <p_i|%braket(|p_i>,V_i).
    output = (this%frequency**2/2)                     &
         & * sum(bra_modes%braket( bra_modes,          &
         &                         harmonic_potential, &
         &                         this%log_2nw_       ))
  endif
end procedure

module procedure kinetic_stress_MonomialStateComplex
  type(MonomialStateComplex), pointer :: monomial_ket
  
  logical,  allocatable :: finite_overlap(:)
  real(dp), allocatable :: overlaps(:)
  
  integer :: i,j
  
  ! |p> = product_i |p_i>
  ! The kinetic stress is given by
  !    S = -(1/NV) sum_i (I_{i,i}d^2/d(u_i^2) + sum_{j/=i}I_{i,j}d^2/du_idu_j).
  ! <p_i|d^2/d(u_i)^2|q_i> and <p_i|d^2/du_idu_j|q_i> are calculated up to
  !    a factor of 2Nw, so
  !    <p|S|q> = -2w * (
  !      sum_i prefactor_{i,i}*<p_i|%second_derivative(|q_i>) * (<p'|q'>)
  !    + sum_{i,j} prefactor_{i,j}*<p_i|%plus_derivative(|q_i>)
  !                               *<p_j|%minus_derivative(|q_j>)*(<p''|q''>) ),
  ! where |p'> is |p> excluding |p_i>, so |p>=|p_i>|p'>,
  ! and |p''> is |p> excluding |p_i> and |p_j>, so |p>=|p_i>|p_j>|p''>.
  
  output = dblemat(zeroes(3,3))
  
  if (present(ket)) then
    monomial_ket => monomial_state_complex_pointer(ket)
    
    finite_overlap = this%modes_%finite_overlap(monomial_ket%modes_)
    if (count(.not.finite_overlap)==0) then
      ! All <p_i|q_i> are finite, so
      !    <p'|q'> = <p|q>/<p_i|q_i>
      !    <p''|q''> = <p|q>/(<p_i|q_i><p_j|q_j>)
      ! For monomial states, if <p_i|q_i>/=0 then
      !    <p_i|d/d(u_i)|q_i>=0 by parity.
      ! S = -2w<p|q>
      !   *  sum_i prefactor_{i,i}*<p_i|%second_derivative(|q_i>)/<p_i|q_i>
      
      overlaps = this%modes_%inner_product(monomial_ket%modes_)
      output = sum( stress_prefactors%prefactor( this%modes_%id(),       &
           &                                     this%modes_%id()  )     &
           &      * this%modes_%second_derivative(monomial_ket%modes_)   &
           &      / overlaps                                           ) &
           & * (-2*this%frequency)*product(overlaps)
    elseif (count(.not.finite_overlap)==2) then
      ! Exactly two <p_i|q_i> are zero. Label these i and j.
      ! <p|S|q> = -2w*
      !   (prefactor_{i,j}+prefactor_{j,i})
      !   * <p_i|%first_derivative(|q_i>)
      !   * <p_j|%first_derivative(|q_j>)
      !   * <p''|q''>
      i = first(.not. finite_overlap)
      j = i + first(.not. finite_overlap(i+1:))
      
      ! Calculate <p|S|q> up to the factor of <p''|q''>.
      output = ( stress_prefactors%prefactor( this%modes_(i)%id(),         &
           &                                  this%modes_(j)%id()  )       &
           &   * this%modes_(i)%plus_derivative(monomial_ket%modes_(i))    &
           &   * this%modes_(j)%minus_derivative(monomial_ket%modes_(j))   &
           &   + stress_prefactors%prefactor( this%modes_(j)%id(),         &
           &                                  this%modes_(i)%id()  )       &
           &   * this%modes_(j)%plus_derivative(monomial_ket%modes_(j))    &
           &   * this%modes_(i)%minus_derivative(monomial_ket%modes_(i)) ) &
           & * (-2*this%frequency)
      
      ! Calculate and include the factor of <p''|q''>.
      overlaps = this%modes_%inner_product(monomial_ket%modes_)
      output = output                     &
           & * product(overlaps(:i-1))    &
           & * product(overlaps(i+1:j-1)) &
           & * product(overlaps(j+1:))
    else
      ! More than two <p_i|q_i>=0, so the whole expression is zero.
      return
    endif
  else
    ! |p>=|q>, so all first derivative expectations are zero.
    ! Also, <p|p>=1, so <p'|p'>=1.
    ! -> <p|S|p> = -2w sum_i prefactor_{i,i}<p_i|%second_derivative(|p_i>)
    output = -2*this%frequency                                   &
         & * sum( stress_prefactors%prefactor( this%modes_%id(), &
         &                                     this%modes_%id()) &
         &      * this%modes_%second_derivative(this%modes_)     )
  endif
  
  ! Divide by the volume.
  output = output / anharmonic_data%structure%volume
end procedure

module procedure change_modes_MonomialStateComplex
  integer, allocatable :: ids(:)
  integer, allocatable :: paired_ids(:)
  integer, allocatable :: powers(:)
  integer, allocatable :: paired_powers(:)
  integer, allocatable :: sort_key(:)
  
  ! Get the ids and powers of the single-mode terms.
  ids = this%modes_%id()
  paired_ids = this%modes_%paired_id()
  powers = this%modes_%power()
  paired_powers = this%modes_%paired_power()
  
  ! Change the ids according to the given group.
  ids = mode_group*ids
  
  ! Sort the modes by id.
  sort_key = sort(ids)
  ids = ids(sort_key)
  paired_ids = paired_ids(sort_key)
  powers = powers(sort_key)
  paired_powers = paired_powers(sort_key)
  
  ! Construct output using the new ids.
  output = MonomialStateComplex(                                        &
     & supercell_size = this%supercell_size,                            &
     & frequency      = this%frequency,                                 &
     & modes          = MonomialState2D( id           = ids,            &
     &                                   paired_id    = paired_ids,     &
     &                                   power        = powers,         &
     &                                   paired_power = paired_powers ) )
end procedure

module procedure read_MonomialStateComplex
  integer                            :: supercell_size
  real(dp)                           :: frequency
  type(MonomialState2D), allocatable :: modes(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i
  
  select type(this); type is(MonomialStateComplex)
    supercell_size = int(token(input(2),4))
    
    frequency = dble(token(input(3),3))
    
    line = split_line(input(5),delimiter='>')
    line = [(line(i)//'>',i=1,size(line))]
    modes = MonomialState2D(line)
    
    this = MonomialStateComplex(supercell_size,frequency,modes)
  class default
    call err()
  end select
end procedure

module procedure write_MonomialStateComplex
  select type(this); type is(MonomialStateComplex)
    output = [ 'Supercell size : '//this%supercell_size, &
             & 'Frequency      : '//this%frequency,      &
             & str('State'),                             &
             & join(str(this%modes_), delimiter='')      ]
  class default
    call err()
  end select
end procedure

module procedure new_MonomialStateComplex_Strings
  call this%read(input)
end procedure

module procedure new_MonomialStateComplex_StringArray
  this = MonomialStateComplex(str(input))
end procedure
end submodule
