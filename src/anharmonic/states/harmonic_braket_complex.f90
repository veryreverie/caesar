! ======================================================================
! A pair of HarmonicStateComplex.
! ======================================================================
module harmonic_braket_complex_module
  use common_module
  
  use anharmonic_common_module
  
  use harmonic_state_complex_module
  implicit none
  
  private
  
  public :: HarmonicBraKetComplex
  
  type, extends(SubspaceBraKet) :: HarmonicBraKetComplex
    type(HarmonicStateComplex), pointer :: bra_
    type(HarmonicStateComplex), pointer :: ket_
    
    real(dp), private :: frequency_
    real(dp), private :: log_2nw_
    integer,  private :: maximum_power_
    integer,  private :: expansion_order_
  contains
    procedure, public :: set_bra_pointer => &
                       & set_bra_pointer_HarmonicBraKetComplex
    procedure, public :: set_ket_pointer => &
                       & set_ket_pointer_HarmonicBraKetComplex
    
    procedure, public :: finite_overlap => &
                       & finite_overlap_HarmonicBraKetComplex
    procedure, public :: inner_product => &
                       & inner_product_HarmonicBraKetComplex
    procedure, public :: integrate => &
                       & integrate_HarmonicBraKetComplex
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_HarmonicBraKetComplex
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_HarmonicBraKetComplex
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_HarmonicBraKetComplex
  end type
  
  interface HarmonicBraKetComplex
    module procedure new_HarmonicBraKetComplex
  end interface
contains

function new_HarmonicBraKetComplex(subspace_id,mode_ids,paired_mode_ids, &
   & frequency,supercell_size,maximum_power,expansion_order) result(this) 
  implicit none
  
  integer,  intent(in)        :: subspace_id
  integer,  intent(in)        :: mode_ids(:)
  integer,  intent(in)        :: paired_mode_ids(:)
  real(dp), intent(in)        :: frequency
  integer,  intent(in)        :: supercell_size
  integer,  intent(in)        :: maximum_power
  integer,  intent(in)        :: expansion_order
  type(HarmonicBraKetComplex) :: this
  
  this%subspace_id      = subspace_id
  this%mode_ids         = mode_ids
  this%paired_mode_ids  = paired_mode_ids
  this%frequency_       = frequency
  this%log_2nw_         = log(2*supercell_size*frequency)
  this%maximum_power_   = maximum_power
  this%expansion_order_ = expansion_order
end function

subroutine set_bra_pointer_HarmonicBraKetComplex(this,bra)
  implicit none
  
  class(HarmonicBraKetComplex), intent(inout)      :: this
  class(SubspaceState),         intent(in), target :: bra
  
  this%bra_ => harmonic_state_complex_pointer(bra)
end subroutine

subroutine set_ket_pointer_HarmonicBraKetComplex(this,ket)
  implicit none
  
  class(HarmonicBraKetComplex), intent(inout)      :: this
  class(SubspaceState),         intent(in), target :: ket
  
  this%ket_ => harmonic_state_complex_pointer(ket)
end subroutine

impure elemental function finite_overlap_HarmonicBraKetComplex(this, &
   & anharmonic_data) result(output)
  implicit none
  
  class(HarmonicBraKetComplex), intent(in) :: this
  type(AnharmonicData),         intent(in) :: anharmonic_data
  logical                                  :: output
  
  if (.not. associated(this%bra_)) then
    call err()
  endif
  
  if (associated(this%ket_)) then
    output = .true.
  else
    output = all(this%bra_%modes_%finite_overlap(this%ket_%modes_))
  endif
end function

impure elemental function inner_product_HarmonicBraKetComplex(this, &
   & anharmonic_data) result(output)
  implicit none
  
  class(HarmonicBraKetComplex), intent(in) :: this
  type(AnharmonicData),         intent(in) :: anharmonic_data
  real(dp)                                 :: output
  
  if (.not. associated(this%bra_)) then
    call err()
  endif
  
  if (associated(this%ket_)) then
    output = product(this%bra_%modes_%inner_product(this%ket_%modes_))
  else
    ! Modes are normalised, so <bra|bra>=1.
    output = 1
  endif
end function

impure elemental function integrate_HarmonicBraKetComplex(this,monomial, &
   & anharmonic_data) result(output) 
  implicit none
  
  class(HarmonicBraKetComplex), intent(in) :: this
  type(SparseMonomial),         intent(in) :: monomial
  type(AnharmonicData),         intent(in) :: anharmonic_data
  complex(dp)                              :: output
  
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
end function

impure elemental function kinetic_energy_HarmonicBraKetComplex(this, &
   & anharmonic_data) result(output)
  implicit none
  
  class(HarmonicBraKetComplex), intent(in) :: this
  type(AnharmonicData),         intent(in) :: anharmonic_data
  real(dp)                                 :: output
  
  logical, allocatable :: modes_overlap(:)
  
  integer :: i
  
  ! |p> = product_i |p_i>
  ! The kinetic energy is given by T = -(1/N) sum_i d^2/d(u_i^2).
  ! <p_i|d^2/d(u_i)^2|q_i> is calculated up to a factor of 2Nw, so
  !    <p|T|q> = -2w * sum_i <p_i|%second_derivative(|q_i>) * (<p'|q'>),
  !    where |p'> is |p> excluding |p_i>, so |p>=|p_i>|p'>.
  ! Since harmonic states are orthonormal, <p'|q'> is either 0 or 1.
  
  if (.not. associated(this%bra_)) then
    call err()
  endif
  
  if (associated(this%ket_)) then
    modes_overlap = this%bra_%modes_%finite_overlap(this%ket_%modes_)
    if (all(modes_overlap)) then
      ! All <p_i|q_i>/=0, so |q>=|p>,
      !    so <p|T|q> = -2w * sum_i <p_i|%second_derivative(|q_i>).
      output = -2*this%frequency_ &
           & * sum(this%bra_%modes_%second_derivative(this%ket_%modes_))
    elseif (count(.not. modes_overlap)==1) then
      ! <p_i|q_i>=0, but all other <p_j|q_j>=1.
      ! <p|T|q> = -2w<p_i%second_derivative(|q_i>).
      i = first(.not. modes_overlap)
      output = -2*this%frequency_ &
           & * this%bra_%modes_(i)%second_derivative(this%ket_%modes_(i))
    else
      ! More than one <p_i|q_i>=0, so <p|V|q>=0.
      output = 0
    endif
  else
    ! |p>=|q>, so <p'|q'>=1, and so
    !    <p|T|q> = -2w*sum_i <p_i|%second_derivative().
    output = -2*this%frequency_ &
         & * sum(this%bra_%modes_%second_derivative(this%bra_%modes_))
  endif
end function

impure elemental function harmonic_potential_energy_HarmonicBraKetComplex( &
   & this,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicBraKetComplex), intent(in) :: this
  type(AnharmonicData),         intent(in) :: anharmonic_data
  real(dp)                                 :: output
  
  type(ComplexUnivariate), allocatable :: harmonic_potential(:)
  
  logical, allocatable :: modes_overlap(:)
  
  integer :: i
  
  ! |p> = product_i |p_i>
  ! The harmonic potential energy is given by V = Nw^2 sum_i (u_i^2).
  ! <p_i|(u_i)^2|q_i> is calculated up to a factor of 2Nw, so
  !    <p|V|q> = (w/2) * sum_i <p_i|%braket(|q_i>,V_i) * (<p'|q'>),
  !    where |p'> is |p> excluding |p_i>, so |p>=|p_i>|p'>.
  ! Since harmonic states are orthonormal, <p'|q'> is either 0 or 1.
  
  if (.not. associated(this%bra_)) then
    call err()
  endif
  
  harmonic_potential = ComplexUnivariate(           &
     & id           = this%bra_%modes_%id(),        &
     & paired_id    = this%bra_%modes_%paired_id(), &
     & power        = 1,                            &
     & paired_power = 1                             )
  
  if (associated(this%ket_)) then
    modes_overlap = this%bra_%modes_%finite_overlap(this%ket_%modes_)
    if (all(modes_overlap)) then
      ! All <p_i|q_i>/=0, so |q>=|p>,
      !    so <p|T|q> = (w^2/2)
      !               * sum_i <p_i|%braket(|q_i>,V_i)
      output = (this%frequency_**2/2)                             &
           & * sum(this%bra_%modes_%braket( this%ket_%modes_,     &
           &                                harmonic_potential,   &
           &                                this%log_2nw_,        &
           &                                this%maximum_power_,  &
           &                                this%expansion_order_ ))
    elseif (count(.not. modes_overlap)==1) then
      ! <p_i|q_i>=0, but all other <p_j|q_j>=1.
      ! <p|T|q> = (w^2/2)<p_i%second_derivative(|q_i>).
      i = first(.not. modes_overlap)
      output = (this%frequency_**2/2)                             &
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
    !    <p|V|q> = (w^2/2)*sum_i <p_i|%braket(|p_i>,V_i).
    output = (this%frequency_**2/2)                             &
         & * sum(this%bra_%modes_%braket( this%bra_%modes_,     &
         &                                harmonic_potential,   &
         &                                this%log_2nw_,        &
         &                                this%maximum_power_,  &
         &                                this%expansion_order_ ))
  endif
  
  output = output
end function

impure elemental function kinetic_stress_HarmonicBraKetComplex(this, &
   & stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(HarmonicBraKetComplex), intent(in) :: this
  type(StressPrefactors),       intent(in) :: stress_prefactors
  type(AnharmonicData),         intent(in) :: anharmonic_data
  type(RealMatrix)                         :: output
  
  logical,  allocatable :: modes_overlap(:)
  
  integer :: i,j,ialloc
  
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
  ! Since harmonic states are orthonormal, <p'|q'> and <p''|q''> are
  !    either 0 or 1.
  
  if (.not. associated(this%bra_)) then
    call err()
  endif
  
  output = dblemat(zeroes(3,3))
  
  if (associated(this%ket_)) then
    modes_overlap = this%bra_%modes_%finite_overlap(this%ket_%modes_)
    if (all(modes_overlap)) then
      ! All <p_i|q_i> are finite, |q>=|p>.
      ! -> <p|S|p> = -2w sum_i prefactor_{i,i}<p_i|%second_derivative(|p_i>)
      output = -2*this%frequency_                                         &
           & * sum( stress_prefactors%prefactor( this%bra_%modes_%id(),   &
           &                                     this%bra_%modes_%id()  ) &
           &      * this%bra_%modes_%second_derivative(this%bra_%modes_)  )
    elseif (count(.not. modes_overlap)==2) then
      ! Exactly two <p_i|q_i> are zero. Label these i and j.
      ! <p|S|q> = -2w*
      !   (prefactor_{i,j}+prefactor_{j,i})
      !   * <p_i|%first_derivative(|q_i>)
      !   * <p_j|%first_derivative(|q_j>)
      i = first(.not. modes_overlap)
      j = i + first(.not. modes_overlap(i+1:))
      
      output = ( stress_prefactors%prefactor( this%bra_%modes_(i)%id(),      &
           &                                  this%bra_%modes_(j)%id()  )    &
           &   * this%bra_%modes_(i)%plus_derivative(this%ket_%modes_(i))    &
           &   * this%bra_%modes_(j)%minus_derivative(this%ket_%modes_(j))   &
           &   + stress_prefactors%prefactor( this%bra_%modes_(j)%id(),      &
           &                                  this%bra_%modes_(i)%id()  )    &
           &   * this%bra_%modes_(j)%plus_derivative(this%ket_%modes_(j))    &
           &   * this%bra_%modes_(i)%minus_derivative(this%ket_%modes_(i)) ) &
           & * (-2*this%frequency_)
    else
      ! More than two <p_i|q_i>=0, so the whole expression is zero.
      return
    endif
  else
    ! |p>=|q>, so all first derivative expectations are zero.
    ! Also, <p|p>=1, so <p'|p'>=1.
    ! -> <p|S|p> = -2w sum_i prefactor_{i,i}<p_i|%second_derivative(|p_i>)
    output = -2*this%frequency_                                         &
         & * sum( stress_prefactors%prefactor( this%bra_%modes_%id(),   &
         &                                     this%bra_%modes_%id()  ) &
         &      * this%bra_%modes_%second_derivative(this%bra_%modes_)  )
  endif
  
  ! Divide by the volume.
  output = output / anharmonic_data%structure%volume
end function
end module