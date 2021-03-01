! ======================================================================
! A monomial state along a single mode at 2q=G,
!    |p> ~ u^p|0>.
! ======================================================================
! N.B. throughout, f(p) is the odd factorial of p, f(p) = (2p)!/(p! 2^p).
! N.B. for speed, any factors of (2Nw) are neglected here,
!    and should be added elsewhere.
!
! If a mode is its own conjugate, (u)* = u,
!    then the single-mode states along mode u are:
! |p> = sqrt((2Nw)^(p) / f(p)) (u)^(p) |0>
! The wavefunction of |0> is sqrt(sqrt(m*w/pi)) exp(- 1/2 N w (u)^2 )
!
! States are normalised <p|p>=1,
!    but are in general not orthogonal <p|q>/=0.
!
! w is the effective frequency of the mode
!    (in general not the same as the harmonic frequency).
!
! m is the geometric average mass, arising as a result of mass-reduction
!    of co-ordinates. This cancels in integrals, and so is not stored.
!
! N is the number of primitive cells in the anharmonic supercell.
module caesar_monomial_state_1d_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  implicit none
  
  private
  
  public :: MonomialState1D
  
  type, extends(Stringable) :: MonomialState1D
    integer, private :: id_
    integer, private :: power_
  contains
    procedure, public :: id           => id_MonomialState1D
    procedure, public :: power        => power_MonomialState1D
    procedure, public :: total_power  => total_power_MonomialState1D
    procedure, public :: wavevector   => wavevector_MonomialState1D
    
    ! ------------------------------
    ! Objects of the form <p|X|q>
    ! ------------------------------
    ! Whether <p|q> is finite.
    procedure, public :: finite_overlap => &
                       & finite_overlap_MonomialState1D
    ! <p|q>.
    procedure, public :: inner_product => &
                       & inner_product_MonomialState1D
    ! <p|u^n|q>.
    procedure, public :: braket => &
                       & braket_MonomialState1D
    ! <p|d/du|q>.
    procedure, public :: first_derivative => &
                       & first_derivative_MonomialState1D
    ! <p|d2/du2|q>.
    procedure, public :: second_derivative => &
                       & second_derivative_MonomialState1D
    
    ! ------------------------------
    ! I/O.
    ! ------------------------------
    procedure, public :: read  => read_MonomialState1D
    procedure, public :: write => write_MonomialState1D
  end type
  
  interface MonomialState1D
    ! Constructor.
    impure elemental module function new_MonomialState1D(id,power) &
       & result(this) 
      integer, intent(in)   :: id
      integer, intent(in)   :: power
      type(MonomialState1D) :: this
    end function
  end interface
  
  interface
    ! Getters.
    impure elemental module function id_MonomialState1D(this) result(output) 
      class(MonomialState1D), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface
    impure elemental module function power_MonomialState1D(this) &
       & result(output) 
      class(MonomialState1D), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface
    ! Returns the total power of the state.
    impure elemental module function total_power_MonomialState1D(this) &
       & result(output) 
      class(MonomialState1D), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface
    ! Returns the wavevector of the state.
    module function wavevector_MonomialState1D(this,modes,qpoints) &
       & result(output) 
      class(MonomialState1D), intent(in) :: this
      type(ComplexMode),      intent(in) :: modes(:)
      type(QpointData),       intent(in) :: qpoints(:)
      type(FractionVector)               :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Objects of the form <p|X|q>, for a variety of operators X.
    ! ----------------------------------------------------------------------
    
    ! Whether <p|q> is non-zero.
    !    = false  if p+q odd.
    !    = true   otherwise.
    impure elemental module function finite_overlap_MonomialState1D(bra,ket) &
       & result(output) 
      class(MonomialState1D), intent(in) :: bra
      class(MonomialState1D), intent(in) :: ket
      logical                            :: output
    end function
  end interface
  
  interface
    ! <p|q>.
    !    = 0                              if p+q odd.
    !    = f((p+q)/2) / sqrt(f(p)*f(q))   otherwise,
    ! where f(x) is the odd factorial of x, f(x) = prod_{k=1}^x [ 2k-1 ].
    impure elemental module function inner_product_MonomialState1D(bra,ket) &
       & result(output) 
      class(MonomialState1D), intent(in) :: bra
      class(MonomialState1D), intent(in) :: ket
      real(dp)                           :: output
    end function
  end interface
  
  interface
    ! <p|(u)^(n)|q>.
    !     = 0                 if p+n+q odd.
    !     = 1/sqrt(2Nw)^{n}
    !     * f((p+n+q)/2)
    !     / sqrt(f(p)*f(q))   otherwise,
    ! where f(x) is the odd factorial of x, f(x) = prod_{k=1}^x [ 2k-1 ].
    impure elemental module function braket_MonomialState1D(bra,ket, &
       & potential,log_2nw) result(output) 
      class(MonomialState1D),  intent(in) :: bra
      class(MonomialState1D),  intent(in) :: ket
      type(ComplexUnivariate), intent(in) :: potential
      real(dp),                intent(in) :: log_2nw
      real(dp)                            :: output
    end function
  end interface
  
  interface
    ! <p|d/du|q>.
    !    = 0                 if p+q even.
    !    = sqrt(2Nw)
    !    * (q-p)/2
    !    * f((p+q-1)/2)
    !    / sqrt(f(p)*f(q))   otherwise,
    ! where f(x) is the odd factorial of x, f(x) = prod_{k=1}^x [ 2k-1 ].
    !
    ! N.B. the factor of sqrt(2Nw) is neglected.
    impure elemental module function first_derivative_MonomialState1D(bra, &
       & ket) result(output) 
      class(MonomialState1D), intent(in) :: bra
      class(MonomialState1D), intent(in) :: ket
      real(dp)                           :: output
    end function
  end interface
  
  interface
    ! <p|d2/du2|q>.
    !    = 0                              if p+q odd.
    !    = 2Nw
    !    * (((p-q)^2-1)/(2(p+q-1)) - 1)
    !    * f((p+q)/2) / sqrt(f(p)*f(q))   otherwise,
    ! where f(x) is the odd factorial of x, f(x) = prod_{k=1}^x [ 2k-1 ].
    !
    ! N.B. the factor of 2Nw is neglected.
    impure elemental module function second_derivative_MonomialState1D(bra, &
       & ket) result(output) 
      class(MonomialState1D), intent(in) :: bra
      class(MonomialState1D), intent(in) :: ket
      real(dp)                           :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_MonomialState1D(this,input) 
      class(MonomialState1D), intent(out) :: this
      type(String),           intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_MonomialState1D(this) result(output) 
      class(MonomialState1D), intent(in) :: this
      type(String)                       :: output
    end function
  end interface
  
  interface MonomialState1D
    impure elemental module function new_MonomialState1D_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(MonomialState1D)    :: this
    end function
  end interface
end module
