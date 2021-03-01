submodule (caesar_harmonic_braket_real_module) caesar_harmonic_braket_real_submodule
  use caesar_states_module
contains

module procedure new_HarmonicBraKetReal
  this%subspace_id      = subspace_id
  this%mode_ids         = mode_ids
  this%paired_mode_ids  = paired_mode_ids
  this%frequency_       = frequency
  this%log_2nw_         = log(2*supercell_size*frequency)
  this%maximum_power_   = maximum_power
  this%expansion_order_ = expansion_order
end procedure

module procedure set_bra_pointer_HarmonicBraKetReal
  this%bra_ => harmonic_state_real_pointer(bra)
end procedure

module procedure set_ket_pointer_HarmonicBraKetReal
  this%ket_ => harmonic_state_real_pointer(ket)
end procedure

module procedure finite_overlap_HarmonicBraKetReal
  if (.not. associated(this%bra_)) then
    call err()
  endif
  
  if (associated(this%ket_)) then
    output = all(this%bra_%modes_%finite_overlap(this%ket_%modes_))
  else
    output = .true.
  endif
end procedure

module procedure inner_product_HarmonicBraKetReal
  if (.not. associated(this%bra_)) then
    call err()
  endif
  
  if (associated(this%ket_)) then
    output = product(this%bra_%modes_%inner_product(this%ket_%modes_))
  else
    ! Modes are normalised, so <bra|bra>=1.
    output = 1
  endif
end procedure

module procedure integrate_HarmonicBraKetReal
  if (.not. associated(this%bra_)) then
    call err()
  endif
  
  if (associated(this%ket_)) then
    output = exp(sum(this%bra_%modes_%log_braket( this%ket_%modes_,     &
                                                & monomial%modes,       &
                                                & this%log_2nw_,        &
                                                & this%maximum_power_,  &
                                                & this%expansion_order_ )))
  else
    output = exp(sum(this%bra_%modes_%log_braket( this%bra_%modes_,     &
                                                & monomial%modes,       &
                                                & this%log_2nw_,        &
                                                & this%maximum_power_,  &
                                                & this%expansion_order_ )))
  endif
end procedure

module procedure kinetic_energy_HarmonicBraKetReal
  logical, allocatable :: modes_overlap(:)
  
  integer :: i
  
  ! |p> = product_i |p_i>
  ! The kinetic energy is given by T = -(1/N) sum_i d^2/d(u_i^2).
  ! <p_i|d^2/d(u_i)^2|q_i> is calculated up to a factor of 2Nw, so
  !    <p|T|q> = -w * sum_i <p_i|%second_derivative(|q_i>) * (<p'|q'>),
  !    where |p'> is |p> excluding |p_i>, so |p>=|p_i>|p'>.
  ! Since harmonic states are orthonormal, <p'|q'> is either 0 or 1.
  
  if (.not. associated(this%bra_)) then
    call err()
  endif
  
  if (associated(this%ket_)) then
    modes_overlap = this%bra_%modes_%finite_overlap(this%ket_%modes_)
    if (all(modes_overlap)) then
      ! All <p_i|q_i>/=0, so |q>=|p>,
      !    so <p|T|q> = -w * sum_i <p_i|%second_derivative(|q_i>).
      output = -this%frequency_ &
           & * sum(this%bra_%modes_%second_derivative(this%ket_%modes_))
    elseif (count(.not. modes_overlap)==1) then
      ! <p_i|q_i>=0, but all other <p_j|q_j>=1.
      ! <p|T|q> = -w<p_i%second_derivative(|q_i>).
      i = first(.not. modes_overlap)
      output = -this%frequency_ &
           & * this%bra_%modes_(i)%second_derivative(this%ket_%modes_(i))
    else
      ! More than one <p_i|q_i>=0, so <p|V|q>=0.
      output = 0
    endif
  else
    ! |p>=|q>, so <p'|q'>=1, and so
    !    <p|T|q> = -w*sum_i <p_i|%second_derivative().
    output = -this%frequency_ &
         & * sum(this%bra_%modes_%second_derivative(this%bra_%modes_))
  endif
end procedure

module procedure harmonic_potential_energy_HarmonicBraKetReal
  type(ComplexUnivariate), allocatable :: harmonic_potential(:)
  
  logical, allocatable :: modes_overlap(:)
  
  integer :: i
  
  ! |p> = product_i |p_i>
  ! The harmonic potential energy is given by V = (Nw^2/2) sum_i (u_i^2).
  ! <p_i|(u_i)^2|q_i> is calculated up to a factor of 2Nw, so
  !    <p|V|q> = (w/4) * sum_i <p_i|%braket(|q_i>,V_i) * (<p'|q'>),
  !    where |p'> is |p> excluding |p_i>, so |p>=|p_i>|p'>.
  ! Since harmonic states are orthonormal, <p'|q'> is either 0 or 1.
  
  if (.not. associated(this%bra_)) then
    call err()
  endif
  
  harmonic_potential = ComplexUnivariate(    &
     & id           = this%bra_%modes_%id(), &
     & paired_id    = this%bra_%modes_%id(), &
     & power        = 2,                     &
     & paired_power = 2                      )
  
  if (associated(this%ket_)) then
    modes_overlap = this%bra_%modes_%finite_overlap(this%ket_%modes_)
    if (all(modes_overlap)) then
      ! All <p_i|q_i>/=0, so all <p'|q'>=<p|q>/<p_i|q_i>,
      !    so <p|T|q> = (w^2/4)
      !               * sum_i <p_i|%braket(|q_i>,V_i)
      output = (this%frequency_**2/4)                             &
           & * sum(this%bra_%modes_%braket( this%ket_%modes_,     &
           &                                harmonic_potential,   &
           &                                this%log_2nw_,        &
           &                                this%maximum_power_,  &
           &                                this%expansion_order_ ))
    elseif (count(.not. modes_overlap)==1) then
      ! <p_i|q_i>=0, but all other <p_j|q_j>=1.
      ! <p|V|q> = (w^2/4)<p_i%second_derivative(|q_i>).
      i = first(.not. modes_overlap)
      output = (this%frequency_**2/4)                             &
           & * this%bra_%modes_(i)%braket( this%ket_%modes_(i),   &
           &                               harmonic_potential(i), &
           &                               this%log_2nw_,         &
           &                               this%maximum_power_,   &
           &                               this%expansion_order_  )
    else
      ! More than one <p_i|q_i>=0, so <p|V|q>=0.
      output = 0
    endif
  else
    ! |p>=|q>, so <p'|q'>=1, and so
    !    <p|V|q> = (w^2/4)*sum_i <p_i|%braket(|p_i>,V_i).
    output = (this%frequency_**2/4)                             &
         & * sum(this%bra_%modes_%braket( this%bra_%modes_,     &
         &                                harmonic_potential,   &
         &                                this%log_2nw_,        &
         &                                this%maximum_power_,  &
         &                                this%expansion_order_ ))
  endif
end procedure

module procedure kinetic_stress_HarmonicBraKetReal
  logical,  allocatable :: modes_overlap(:)
  
  integer :: i,j
  
  ! |p> = product_i |p_i>
  ! The kinetic stress is given by
  !    S = -(1/NV) sum_i (I_{i,i}d^2/d(u_i^2) + sum_{j/=i}I_{i,j}d^2/du_idu_j).
  ! <p_i|d^2/d(u_i)^2|q_i> and <p_i|d^2/du_idu_j|q_i> are calculated up to
  !    a factor of 2Nw, so
  !    <p|S|q> = -2w * (
  !      sum_i prefactor_{i,i}*<p_i|%second_derivative(|q_i>) * (<p'|q'>)
  !    + sum_{i,j} prefactor_{i,j}*<p_i|%first_derivative(|q_i>)
  !                               *<p_j|%first_derivative(|q_j>)*(<p''|q''>) ),
  ! where |p'> is |p> excluding |p_i>, so |p>=|p_i>|p'>,
  ! and |p''> is |p> excluding |p_i> and |p_j>, so |p>=|p_i>|p_j>|p''>.
  ! Since harmonic states are orthonormal, <p'|q'> and <p''|q''> are
  !    either 0 or 1.
  
  if (.not. associated(this%bra_)) then
    call err()
  endif
  
  output = dblemat(zeroes(3,3))
  
  if (associated(this%ket_)) then
    modes_overlap = this%bra_%modes_%finite_overlap(this%ket_%modes_)
    if (all(modes_overlap)) then
      ! All <p_i|q_i> are finite, so |q>=|p>.
      ! -> <p|S|p> = -w sum_i prefactor_{i,i}<p_i|%second_derivative(|p_i>)
      output = -2*this%frequency_                                         &
           & * sum( stress_prefactors%prefactor( this%bra_%modes_%id(),   &
           &                                     this%bra_%modes_%id()  ) &
           &      * this%bra_%modes_%second_derivative(this%bra_%modes_)  )
    elseif (count(.not. modes_overlap)==2) then
      ! Exactly two <p_i|q_i> are zero. Label these i and j.
      ! <p|S|q> = -w*
      !   (prefactor_{i,j}+prefactor_{j,i})
      !   * <p_i|%first_derivative(|q_i>)
      !   * <p_j|%first_derivative(|q_j>)
      !   * <p''|q''>
      i = first(.not. modes_overlap)
      j = i + first(.not. modes_overlap(i+1:))
      
      output = ( stress_prefactors%prefactor( this%bra_%modes_(i)%id(),     &
           &                                  this%bra_%modes_(j)%id()  )   &
           &   + stress_prefactors%prefactor( this%bra_%modes_(j)%id(),     &
           &                                  this%bra_%modes_(i)%id()  ) ) &
           & * this%bra_%modes_(i)%first_derivative(this%ket_%modes_(i)) &
           & * this%bra_%modes_(j)%first_derivative(this%ket_%modes_(j)) &
           & * (-2*this%frequency_)
    else
      ! More than two <p_i|q_i>=0, so the whole expression is zero.
      return
    endif
  else
    ! |p>=|q>, so all first derivative expectations are zero.
    ! -> <p|S|p> = -w sum_i prefactor_{i,i}<p_i|%second_derivative(|p_i>)
    output = -2*this%frequency_                                         &
         & * sum( stress_prefactors%prefactor( this%bra_%modes_%id(),   &
         &                                     this%bra_%modes_%id()  ) &
         &      * this%bra_%modes_%second_derivative(this%bra_%modes_)  )
  endif
  
  ! Divide by the volume.
  output = output / anharmonic_data%structure%volume
end procedure
end submodule
