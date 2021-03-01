! ======================================================================
! A harmonic state along a mode at q and its -q pair.
!    |n_i,n_j> ~ (a_i^n_i)(a_j^n_j)|0_i,0_j>.
! ======================================================================
! N.B. throughout, f(p) is the odd factorial of p, f(p) = (2p)!/(p! 2^p).
! N.B. for speed, any factors of (2Nw) are neglected here,
!    and should be added elsewhere.
!
! If a mode is not its own conjugate, (u_i)* = u_j,
!    then the double-mode states along modes u_i and u_j are:
! |p_i,p_j> = 1/sqrt(p_i!p_j!) a_i^p_i a_j^p_j |0_i,0_j>,
!    where a_i is the harmonic creation operator along mode i.
! |0_i,0_j> = sqrt(2*m*w/pi) exp(- N w |u_i|^2 )
!
! States are orthonormalised, <p_i,p_j|p_i,p_j>=1 and
!    <p_i,p_j|q_i,q_j>=0 if p_i/=q_i or p_j/=q_j.
!
! w is the effective frequency the pair of modes
!    (in general not the same as the harmonic frequency).
!
! m is the geometric average mass, arising as a result of mass-reduction
!    of co-ordinates. This cancels in integrals, and so is not stored.
!
! N is the number of primitive cells in the anharmonic supercell.
module caesar_harmonic_state_2d_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  implicit none
  
  private
  
  public :: HarmonicState2D
  
  type, extends(Stringable) :: HarmonicState2D
    integer, private :: id_
    integer, private :: paired_id_
    integer, private :: occupation_
    integer, private :: paired_occupation_
  contains
    procedure, public :: id                => id_HarmonicState2D
    procedure, public :: paired_id         => paired_id_HarmonicState2D
    procedure, public :: occupation        => occupation_HarmonicState2D
    procedure, public :: paired_occupation => paired_occupation_HarmonicState2D
    procedure, public :: total_occupation  => total_occupation_HarmonicState2D
    procedure, public :: wavevector        => wavevector_HarmonicState2D
    
    ! ------------------------------
    ! Objects of the form <p_i,p_j|X|q_i,q_j>
    ! ------------------------------
    ! Whether <p_i,p_j|q_i,q_j> is finite.
    procedure, public :: finite_overlap => &
                       & finite_overlap_HarmonicState2D
    ! <p_i,p_j|q_i,q_j>.
    procedure, public :: inner_product => &
                       & inner_product_HarmonicState2D
    ! <p_i,p_j|(u_i)^(n_i)(u_j)^(n_j)|q_i,q_j>.
    procedure, public :: braket => &
                       & braket_HarmonicState2D
    procedure, public :: log_braket => &
                       & log_braket_HarmonicState2D
    ! <p_i,p_j|d/d(u_i)|q_i,q_j>.
    procedure, public :: plus_derivative => &
                       & plus_derivative_HarmonicState2D
    ! <p_i,p_j|d/d(u_j)|q_i,q_j>.
    procedure, public :: minus_derivative => &
                       & minus_derivative_HarmonicState2D
    ! <p_i,p_j|d2/d(u_i)d(u_j)|q_i,q_j>.
    procedure, public :: second_derivative => &
                       & second_derivative_HarmonicState2D
    
    ! ------------------------------
    ! I/O.
    ! ------------------------------
    procedure, public :: read  => read_HarmonicState2D
    procedure, public :: write => write_HarmonicState2D
  end type
  
  interface HarmonicState2D
    ! Constructor.
    impure elemental module function new_HarmonicState2D(id,paired_id, &
       & occupation,paired_occupation) result(this) 
      integer, intent(in)   :: id
      integer, intent(in)   :: paired_id
      integer, intent(in)   :: occupation
      integer, intent(in)   :: paired_occupation
      type(HarmonicState2D) :: this
    end function
  end interface
  
  interface
    ! Getters.
    impure elemental module function id_HarmonicState2D(this) result(output) 
      class(HarmonicState2D), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface
    impure elemental module function paired_id_HarmonicState2D(this) &
       & result(output) 
      class(HarmonicState2D), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface
    impure elemental module function occupation_HarmonicState2D(this) &
       & result(output) 
      class(HarmonicState2D), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface
    impure elemental module function paired_occupation_HarmonicState2D(this) &
       & result(output) 
      class(HarmonicState2D), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface
    ! Returns the total occupation of the state.
    impure elemental module function total_occupation_HarmonicState2D(this) &
       & result(output) 
      class(HarmonicState2D), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface
    ! Returns the wavevector of the state.
    module function wavevector_HarmonicState2D(this,modes,qpoints) &
       & result(output) 
      class(HarmonicState2D), intent(in) :: this
      type(ComplexMode),      intent(in) :: modes(:)
      type(QpointData),       intent(in) :: qpoints(:)
      type(FractionVector)               :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Objects of the form <p|X|q>, for a variety of operators X.
    ! ----------------------------------------------------------------------
    
    ! Whether <p_i,p_j|q_i,q_j> is non-zero.
    !    = true    if p_i=q_i and p_j=q_j.
    !    = false   otherwise.
    impure elemental module function finite_overlap_HarmonicState2D(bra,ket) &
       & result(output) 
      class(HarmonicState2D), intent(in) :: bra
      class(HarmonicState2D), intent(in) :: ket
      logical                            :: output
    end function
  end interface
  
  interface
    ! <p_i,p_j|q_i,q_j>.
    !    = 1   if p_i=q_i and p_j=q_j.
    !    = 0   otherwise.
    impure elemental module function inner_product_HarmonicState2D(bra,ket) &
       & result(output) 
      class(HarmonicState2D), intent(in) :: bra
      class(HarmonicState2D), intent(in) :: ket
      integer                            :: output
    end function
  end interface
  
  interface
    ! <p_i,p_j|(u_i)^(n_i)(u_j)^(n_j)|q_i,q_j>.
    !    = 0                                 if p_i-p_j-n_i+n_j-q_i+q_j /= 0.
    !    = 0                                 if n_i < p_i-q_i.
    !    = 0                                 if n_i < q_j-p_j.
    !    = 0                                 if n_j < p_j-q_j.
    !    = 0                                 if n_j < q_i-p_i.
    !    = 1/sqrt(2Nw)^{n_i+n_j}
    !    * sqrt(p_i!q_i!/(p_j!q_j!))
    !    * sum_{k=max(0,p_i-n_i,q_i-n_j)}^{min(p_i,q_i)}[
    !                                    binom(n_i,p_i-k)
    !                                  * binom(n_j,q_i-k)
    !                                  * (n_i+p_j-p_i+k)!
    !                                  / k!               ] otherwise.
    impure elemental module function braket_HarmonicState2D(bra,ket, &
       & potential,log_2nw,maximum_power,expansion_order) result(output) 
      class(HarmonicState2D),  intent(in) :: bra
      class(HarmonicState2D),  intent(in) :: ket
      type(ComplexUnivariate), intent(in) :: potential
      real(dp),                intent(in) :: log_2nw
      integer,                 intent(in) :: maximum_power
      integer,                 intent(in) :: expansion_order
      real(dp)                            :: output
    end function
  end interface
  
  interface
    impure elemental module function log_braket_HarmonicState2D(bra,ket, &
       & potential,log_2nw,maximum_power,expansion_order) result(output) 
      class(HarmonicState2D),  intent(in) :: bra
      class(HarmonicState2D),  intent(in) :: ket
      type(ComplexUnivariate), intent(in) :: potential
      real(dp),                intent(in) :: log_2nw
      integer,                 intent(in) :: maximum_power
      integer,                 intent(in) :: expansion_order
      real(dp)                            :: output
    end function
  end interface
  
  interface
    ! <p_i,p_j|d/d(u_i)|q_i,q_j>.
    !    = -sqrt(2Nw) * sqrt(p_i)/2   if p_i-q_i=1 and p_j=q_j.
    !    =  sqrt(2Nw) * sqrt(q_j)/2   if q_j-p_j=1 and p_i=q_i.
    !    =  0                         otherwise.
    !
    ! N.B. the factor of sqrt(2Nw) is neglected.
    impure elemental module function plus_derivative_HarmonicState2D(bra,ket) &
       & result(output) 
      class(HarmonicState2D), intent(in) :: bra
      class(HarmonicState2D), intent(in) :: ket
      real(dp)                           :: output
    end function
  end interface
  
  interface
    ! <p_i,p_j|d/d(u_j)|q_i,q_j>.
    !    = -sqrt(2Nw) * sqrt(p_j)/2   if p_j-q_j=1 and p_i=q_i.
    !    =  sqrt(2Nw) * sqrt(q_i)/2   if q_i-p_i=1 and p_j=q_j.
    !    =  0                         otherwise.
    !
    ! N.B. the factor of sqrt(2Nw) is neglected.
    impure elemental module function minus_derivative_HarmonicState2D(bra, &
       & ket) result(output) 
      class(HarmonicState2D), intent(in) :: bra
      class(HarmonicState2D), intent(in) :: ket
      real(dp)                           :: output
    end function
  end interface
  
  interface
    ! <p_i,p_j|d2/d(u_i)d(u_j)|q_i,q_j>.
    !    = -2Nw * (p_i+p_j+1)/4     if p_i=q_i and p_j=q_j.
    !    =  2Nw * sqrt(p_i*p_j)/4   if p_i-q_i=p_j-q_j=1.
    !    =  2Nw * sqrt(q_i*q_j)/4   if q_i-p_i=q_j-p_j=1.
    !    =  0                       otherwise.
    !
    ! N.B. the factor of 2Nw is neglected.
    impure elemental module function second_derivative_HarmonicState2D(bra, &
       & ket) result(output) 
      class(HarmonicState2D), intent(in) :: bra
      class(HarmonicState2D), intent(in) :: ket
      real(dp)                           :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_HarmonicState2D(this,input) 
      class(HarmonicState2D), intent(out) :: this
      type(String),           intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_HarmonicState2D(this) result(output) 
      class(HarmonicState2D), intent(in) :: this
      type(String)                       :: output
    end function
  end interface
  
  interface HarmonicState2D
    impure elemental module function new_HarmonicState2D_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(HarmonicState2D)    :: this
    end function
  end interface
end module
