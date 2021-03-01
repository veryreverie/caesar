! ======================================================================
! A monomial state along a mode at q and its -q pair,
!    |p_i,p_j> ~ (u_i^p_i)(u_j^p_j)|0_i,0_j>.
! ======================================================================
! N.B. throughout, f(p) is the odd factorial of p, f(p) = (2p)!/(p! 2^p).
! N.B. for speed, any factors of (2Nw) are neglected here,
!    and should be added elsewhere.
!
! If a mode is not its own conjugate, (u_i)* = u_j,
!    then the double-mode states along modes u_i and u_j are:
! |p_i,p_j> = sqrt((2Nw)^(p_i+p_j) / (p_i+p_j)!) u_i^p_i u_j^p_j |0_i,0_j>
! |0_i,0_j> = sqrt(2*m*w/pi) exp(- N w |u_i|^2 )
!
! States are normalised <p_i,p_j|p_i,p_j>=1,
!    but are in general not orthogonal <p_i,p_j|q_i,q_j>/=0.
!
! w is the effective frequency the pair of modes
!    (in general not the same as the harmonic frequency).
!
! m is the geometric average mass, arising as a result of mass-reduction
!    of co-ordinates. This cancels in integrals, and so is not stored.
!
! N is the number of primitive cells in the anharmonic supercell.
module caesar_monomial_state_2d_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  implicit none
  
  private
  
  public :: MonomialState2D
  
  type, extends(Stringable) :: MonomialState2D
    integer, private :: id_
    integer, private :: paired_id_
    integer, private :: power_
    integer, private :: paired_power_
  contains
    procedure, public :: id           => id_MonomialState2D
    procedure, public :: paired_id    => paired_id_MonomialState2D
    procedure, public :: power        => power_MonomialState2D
    procedure, public :: paired_power => paired_power_MonomialState2D
    procedure, public :: total_power  => total_power_MonomialState2D
    procedure, public :: wavevector   => wavevector_MonomialState2D
    
    ! ------------------------------
    ! Objects of the form <p_i,p_j|X|q_i,q_j>
    ! ------------------------------
    ! Whether <p_i,p_j|q_i,q_j> is finite.
    procedure, public :: finite_overlap => &
                       & finite_overlap_MonomialState2D
    ! <p_i,p_j|q_i,q_j>.
    procedure, public :: inner_product => &
                       & inner_product_MonomialState2D
    ! <p_i,p_j|(u_i)^(n_i)(u_j)^(n_j)|q_i,q_j>.
    procedure, public :: braket => &
                       & braket_MonomialState2D
    ! <p_i,p_j|d/d(u_i)|q_i,q_j>.
    procedure, public :: plus_derivative => &
                       & plus_derivative_MonomialState2D
    ! <p_i,p_j|d/d(u_j)|q_i,q_j>.
    procedure, public :: minus_derivative => &
                       & minus_derivative_MonomialState2D
    ! <p_i,p_j|d2/d(u_i)d(u_j)|q_i,q_j>.
    procedure, public :: second_derivative => &
                       & second_derivative_MonomialState2D
    
    ! ------------------------------
    ! I/O.
    ! ------------------------------
    procedure, public :: read  => read_MonomialState2D
    procedure, public :: write => write_MonomialState2D
  end type
  
  interface MonomialState2D
    ! Constructor.
    impure elemental module function new_MonomialState2D(id,paired_id,power, &
       & paired_power) result(this) 
      integer, intent(in)   :: id
      integer, intent(in)   :: paired_id
      integer, intent(in)   :: power
      integer, intent(in)   :: paired_power
      type(MonomialState2D) :: this
    end function
  end interface
  
  interface
    ! Getters.
    impure elemental module function id_MonomialState2D(this) result(output) 
      class(MonomialState2D), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface
    impure elemental module function paired_id_MonomialState2D(this) &
       & result(output) 
      class(MonomialState2D), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface
    impure elemental module function power_MonomialState2D(this) &
       & result(output) 
      class(MonomialState2D), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface
    impure elemental module function paired_power_MonomialState2D(this) &
       & result(output) 
      class(MonomialState2D), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface
    ! Returns the total power of the state.
    impure elemental module function total_power_MonomialState2D(this) &
       & result(output) 
      class(MonomialState2D), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface
    ! Returns the wavevector of the state.
    module function wavevector_MonomialState2D(this,modes,qpoints) &
       & result(output) 
      class(MonomialState2D), intent(in) :: this
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
    !    = false    if p_i-p_j-q_i+q_j /= 0.
    !    = true     otherwise.
    impure elemental module function finite_overlap_MonomialState2D(bra,ket) &
       & result(output) 
      class(MonomialState2D), intent(in) :: bra
      class(MonomialState2D), intent(in) :: ket
      logical                            :: output
    end function
  end interface
  
  interface
    ! <p_i,p_j|q_i,q_j>.
    !    = 0                          if p_i-p_j-q_i+q_j /= 0.
    !    = ((p_i+p_j+q_i+q_j)/2)!
    !    / sqrt((p_i+p_j)!(q_i+q_j)!) otherwise.
    impure elemental module function inner_product_MonomialState2D(bra,ket) &
       & result(output) 
      class(MonomialState2D), intent(in) :: bra
      class(MonomialState2D), intent(in) :: ket
      real(dp)                           :: output
    end function
  end interface
  
  interface
    ! <p_i,p_j|(u_i)^(n_i)(u_j)^(n_j)|q_i,q_j>.
    !    = 0                                 if p_i-p_j-n_i+n_j-q_i+q_j /= 0.
    !    = 1/sqrt(2Nw)^{n_i+n_j}
    !    * ((p_i+p_j+n_i+n_j+q_i+q_j)/2)!
    !    / sqrt( (p_i+p_j)! * (q_i+q_j)! )   otherwise.
    impure elemental module function braket_MonomialState2D(bra,ket, &
       & potential,log_2nw) result(output) 
      class(MonomialState2D),  intent(in) :: bra
      class(MonomialState2D),  intent(in) :: ket
      type(ComplexUnivariate), intent(in) :: potential
      real(dp),                intent(in) :: log_2nw
      real(dp)                            :: output
    end function
  end interface
  
  interface
    ! <p_i,p_j|d/d(u_i)|q_i,q_j>.
    !    = 0                                       if q_i-p_i-q_j+p_j /= 1
    !    = sqrt(2Nw)
    !    * (q_i - sqrt((q+1)/(p+1))*(p+q+1)/4)
    !    * ((p+q-1)/2)!/sqrt(p!q!)                 otherwise,
    ! where p=p_i+p_j and q=q_i+q_j
    !
    ! N.B. the factor of sqrt(2Nw) is neglected.
    impure elemental module function plus_derivative_MonomialState2D(bra,ket) &
       & result(output) 
      class(MonomialState2D), intent(in) :: bra
      class(MonomialState2D), intent(in) :: ket
      real(dp)                           :: output
    end function
  end interface
  
  interface
    ! <p_i,p_j|d/d(u_j)|q_i,q_j>.
    !    = 0                                       if q_j-p_j-q_i+p_i /= 1
    !    = sqrt(2Nw)
    !    * (q_j - sqrt((q+1)/(p+1))*(p+q+1)/4)
    !    * ((p+q-1)/2)!/sqrt(p!q!)                 otherwise,
    ! where p=p_i+p_j and q=q_i+q_j
    !
    ! N.B. the factor of sqrt(2Nw) is neglected.
    impure elemental module function minus_derivative_MonomialState2D(bra, &
       & ket) result(output) 
      class(MonomialState2D), intent(in) :: bra
      class(MonomialState2D), intent(in) :: ket
      real(dp)                           :: output
    end function
  end interface
  
  interface
    ! <p_i,p_j|d2/d(u_i)d(u_j)|q_i,q_j>.
    !    = 0  if p_i-p_j-q_i+q_j /= 0.
    !    = 2Nw
    !    * ((p_i-q_j)(p_j-q_i)/(2(p+q))-1/4)
    !    * ((p+q)/2)!/sqrt(p!q!)
    ! where p=p_i+p_j and q=q_i+q_j.
    !
    ! N.B. the factor of 2Nw is neglected.
    impure elemental module function second_derivative_MonomialState2D(bra, &
       & ket) result(output) 
      class(MonomialState2D), intent(in) :: bra
      class(MonomialState2D), intent(in) :: ket
      real(dp)                           :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_MonomialState2D(this,input) 
      class(MonomialState2D), intent(out) :: this
      type(String),           intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_MonomialState2D(this) result(output) 
      class(MonomialState2D), intent(in) :: this
      type(String)                       :: output
    end function
  end interface
  
  interface MonomialState2D
    impure elemental module function new_MonomialState2D_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(MonomialState2D)    :: this
    end function
  end interface
end module
