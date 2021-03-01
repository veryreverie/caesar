! ======================================================================
! A harmonic state along a single mode at 2q=G,
!    |p> = a^p|0>.
! ======================================================================
! N.B. throughout, f(p) is the odd factorial of p, f(p) = (2p)!/(p! 2^p).
! N.B. for speed, any factors of (2Nw) are neglected here,
!    and should be added elsewhere.
!
! If a mode is its own conjugate, (u)* = u,
!    then the single-mode states along mode u are:
! |p> = 1/sqrt(p!) (a)^(p) |0>,
!    where a is the harmonic creation operator.
! The wavefunction of |0> is sqrt(sqrt(m*w/pi)) exp(- 1/2 N w (u)^2 ).
!
! States are orthonormalised, <p|p>=1 and <p|q>=0 if p/=q.
!
! w is the effective frequency of the mode
!    (in general not the same as the harmonic frequency).
!
! m is the geometric average mass, arising as a result of mass-reduction
!    of co-ordinates. This cancels in integrals, and so is not stored.
!
! N is the number of primitive cells in the anharmonic supercell.
module caesar_harmonic_state_1d_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  implicit none
  
  private
  
  public :: HarmonicState1D
  
  type, extends(Stringable) :: HarmonicState1D
    integer, private :: id_
    integer, private :: occupation_
  contains
    procedure, public :: id               => id_HarmonicState1D
    procedure, public :: occupation       => occupation_HarmonicState1D
    procedure, public :: total_occupation => total_occupation_HarmonicState1D
    procedure, public :: wavevector       => wavevector_HarmonicState1D
    
    ! ------------------------------
    ! Objects of the form <p|X|q>
    ! ------------------------------
    ! Whether <p|q> is finite.
    procedure, public :: finite_overlap => &
                       & finite_overlap_HarmonicState1D
    ! <p|q>.
    procedure, public :: inner_product => &
                       & inner_product_HarmonicState1D
    ! <p|u^n|q>.
    procedure, public :: braket => &
                       & braket_HarmonicState1D
    procedure, public :: log_braket => &
                       & log_braket_HarmonicState1D
    ! <p|d/du|q>.
    procedure, public :: first_derivative => &
                       & first_derivative_HarmonicState1D
    ! <p|d2/du2|q>.
    procedure, public :: second_derivative => &
                       & second_derivative_HarmonicState1D
    
    ! ------------------------------
    ! I/O.
    ! ------------------------------
    procedure, public :: read  => read_HarmonicState1D
    procedure, public :: write => write_HarmonicState1D
  end type
  
  interface HarmonicState1D
    ! Constructor.
    impure elemental module function new_HarmonicState1D(id,occupation) &
       & result(this) 
      integer, intent(in)   :: id
      integer, intent(in)   :: occupation
      type(HarmonicState1D) :: this
    end function
  end interface
  
  interface
    ! Getters.
    impure elemental module function id_HarmonicState1D(this) result(output) 
      class(HarmonicState1D), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface
    impure elemental module function occupation_HarmonicState1D(this) &
       & result(output) 
      class(HarmonicState1D), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface
    ! Returns the total occupation of the state.
    impure elemental module function total_occupation_HarmonicState1D(this) &
       & result(output) 
      class(HarmonicState1D), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface
    ! Returns the wavevector of the state.
    module function wavevector_HarmonicState1D(this,modes,qpoints) &
       & result(output) 
      class(HarmonicState1D), intent(in) :: this
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
    !    = true   if p=q.
    !    = false  otherwise.
    impure elemental module function finite_overlap_HarmonicState1D(bra,ket) &
       & result(output) 
      class(HarmonicState1D), intent(in) :: bra
      class(HarmonicState1D), intent(in) :: ket
      logical                            :: output
    end function
  end interface
  
  interface
    ! <p|q>.
    !    = 1   if p=q.
    !    = 0   otherwise.
    impure elemental module function inner_product_HarmonicState1D(bra,ket) &
       & result(output) 
      class(HarmonicState1D), intent(in) :: bra
      class(HarmonicState1D), intent(in) :: ket
      integer                            :: output
    end function
  end interface
  
  interface
    ! <p|(u)^(n)|q>.  (Defining d=|p-q|).
    !     = 0                 if p+n+q odd.
    !     = 0                 if n<d.
    !     = 1/sqrt(2Nw)^{n}
    !     * (n!/(((n+d)/2)!2^((n-d)/2))
    !     * sqrt(max(p,q)!/min(p,q)!)
    !     * sum_{k=0}^{min((n-d)/2,p,q)}[ 2^k 
    !                                   * binom((n+d)/2,d+k)
    !                                   * binom(min(p,q),k)  ] otherwise.
    impure elemental module function braket_HarmonicState1D(bra,ket, &
       & potential,log_2nw,maximum_power,expansion_order) result(output) 
      class(HarmonicState1D),  intent(in) :: bra
      class(HarmonicState1D),  intent(in) :: ket
      type(ComplexUnivariate), intent(in) :: potential
      real(dp),                intent(in) :: log_2nw
      integer,                 intent(in) :: maximum_power
      integer,                 intent(in) :: expansion_order
      real(dp)                            :: output
    end function
  end interface
  
  interface
    impure elemental module function log_braket_HarmonicState1D(bra,ket, &
       & potential,log_2nw,maximum_power,expansion_order) result(output) 
      class(HarmonicState1D),  intent(in) :: bra
      class(HarmonicState1D),  intent(in) :: ket
      type(ComplexUnivariate), intent(in) :: potential
      real(dp),                intent(in) :: log_2nw
      integer,                 intent(in) :: maximum_power
      integer,                 intent(in) :: expansion_order
      real(dp)                            :: output
    end function
  end interface
  
  interface
    ! <p|d/du|q>.
    !    = -sqrt(2Nw) * sqrt(p)/2   if p=q+1.
    !    =  sqrt(2Nw) * sqrt(q)/2   if q=p+1.
    !    =  0                       otherwise.
    !
    ! N.B. the factor of sqrt(2Nw) is neglected.
    impure elemental module function first_derivative_HarmonicState1D(bra, &
       & ket) result(output) 
      class(HarmonicState1D), intent(in) :: bra
      class(HarmonicState1D), intent(in) :: ket
      real(dp)                           :: output
    end function
  end interface
  
  interface
    ! <p|d2/du2|q>.
    !    = -2Nw (p+0.5)/2                       if p=q.
    !    =  2Nw sqrt(max(p,q)*(max(p,q)-1))/4   if |p-q|=2.
    !    =  0                                   otherwise.
    !
    ! N.B. the factor of 2Nw is neglected.
    impure elemental module function second_derivative_HarmonicState1D(bra, &
       & ket) result(output) 
      class(HarmonicState1D), intent(in) :: bra
      class(HarmonicState1D), intent(in) :: ket
      real(dp)                           :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_HarmonicState1D(this,input) 
      class(HarmonicState1D), intent(out) :: this
      type(String),           intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_HarmonicState1D(this) result(output) 
      class(HarmonicState1D), intent(in) :: this
      type(String)                       :: output
    end function
  end interface
  
  interface HarmonicState1D
    impure elemental module function new_HarmonicState1D_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(HarmonicState1D)    :: this
    end function
  end interface
end module
