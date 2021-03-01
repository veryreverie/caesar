! ======================================================================
! Calculates and stores 1-D mode overlaps for interpolating polynomials
!    from a coarse q-point grid to a fine q-point grid.
! ======================================================================
! The absolute overlap for atom j and modes u1 and u2 is just <u1|j><j|u2>,
!    i.e. the inner product of the projections of u1 and u2 onto atom j.
! The relative overlap for atom j and modes u1 and u2 is
!    sum_k <u1|k>e^{2 pi i (q1-q2).R_{jk}}<k|u2>, where R_{jk} is the minimum
!    image R-vector between atoms j and k.
! The overlap between two products of modes {u1} and {u2} is
!    sum_atom sum_j [ absolute(u1j,u2j,atom) 
!                   * product_{k/=j} [ relative(u1k,u2k,atom) ] ]
! N.B. the product overlap must be averaged over permutations of {u1} against
!    {u2}.
module caesar_polynomial_interpolator_module
  use caesar_common_module
  use caesar_anharmonic_common_module
  use caesar_permutation_module
  implicit none
  
  private
  
  public :: PolynomialInterpolator
  public :: construct_fine_monomials
  
  type, extends(NoDefaultConstructor) :: PolynomialInterpolator
    complex(dp), allocatable, private :: absolute_(:,:,:)
    complex(dp), allocatable, private :: relative_(:,:,:)
    integer,     allocatable, private :: fine_map_(:)
    integer,     allocatable, private :: coarse_map_(:)
  contains
    procedure, private :: absolute_overlap
    procedure, private :: relative_overlap
    procedure, private :: product_overlap
    
    generic,   public  :: overlap =>                               &
                        & overlap_ComplexMonomial_ComplexMonomial, &
                        & overlap_ComplexMonomial_ComplexPolynomial
    procedure, private :: overlap_ComplexMonomial_ComplexMonomial
    procedure, private :: overlap_ComplexMonomial_ComplexPolynomial
  end type
  
  interface PolynomialInterpolator
    ! Constructor.
    module function new_PolynomialInterpolator(min_images,anharmonic_data, &
       & interpolated_anharmonic_data) result(this) 
      type(MinImages),      intent(in) :: min_images(:,:)
      type(AnharmonicData), intent(in) :: anharmonic_data
      type(AnharmonicData), intent(in) :: interpolated_anharmonic_data
      type(PolynomialInterpolator)     :: this
    end function
  end interface
  
  interface
    ! Calculates <u1|j><j|u2>.
    impure elemental module function absolute_overlap(this,fine_mode, &
       & coarse_mode,atom) result(output) 
      class(PolynomialInterpolator), intent(in) :: this
      integer,                       intent(in) :: fine_mode
      integer,                       intent(in) :: coarse_mode
      integer,                       intent(in) :: atom
      complex(dp)                               :: output
    end function
  end interface
  
  interface
    ! Calculates sum_k <u1|k>e^{2 pi i (q1-q2).R_{jk}}<k|u2>.
    impure elemental module function relative_overlap(this,fine_mode, &
       & coarse_mode,atom) result(output) 
      class(PolynomialInterpolator), intent(in) :: this
      integer,                       intent(in) :: fine_mode
      integer,                       intent(in) :: coarse_mode
      integer,                       intent(in) :: atom
      complex(dp)                               :: output
    end function
  end interface
  
  interface
    ! Calculates <{u1}|{u2}> for a single permutation of {u1} and {u2}.
    module function product_overlap(this,fine_modes,coarse_modes) &
       & result(output) 
      class(PolynomialInterpolator), intent(in) :: this
      integer,                       intent(in) :: fine_modes(:)
      integer,                       intent(in) :: coarse_modes(:)
      complex(dp)                               :: output
    end function
  end interface
  
  interface
    ! Calculates the overlap between a monomial on the fine grid and a monomial
    !    on the coarse grid.
    impure elemental module function overlap_ComplexMonomial_ComplexMonomial(this,fine,coarse) result(output) 
      class(PolynomialInterpolator), intent(in) :: this
      type(ComplexMonomial),         intent(in) :: fine
      type(ComplexMonomial),         intent(in) :: coarse
      complex(dp)                               :: output
    end function
  end interface
  
  interface
    module function to_ids(monomial) result(output) 
      type(ComplexMonomial), intent(in) :: monomial
      integer, allocatable              :: output(:)
    end function
  end interface
  
  interface
    ! Calculates the number of ways a given monomial can be permuted.
    module function log_no_permutations(monomial) result(output) 
      type(ComplexMonomial), intent(in) :: monomial
      real(dp)                          :: output
    end function
  end interface
  
  interface
    ! Calculates the overlap between an monomial on the fine grid and a polynomial
    !    on the coarse grid.
    impure elemental module function overlap_ComplexMonomial_ComplexPolynomial(   this,fine,coarse) result(output) 
      class(PolynomialInterpolator), intent(in) :: this
      type(ComplexMonomial),         intent(in) :: fine
      type(ComplexPolynomial),       intent(in) :: coarse
      complex(dp)                               :: output
    end function
  end interface
  
  interface
    ! Construct the monomials at a q-point on the fine grid.
    module function construct_fine_monomials(power,modes) result(output) 
      integer,           intent(in)      :: power
      type(ComplexMode), intent(in)      :: modes(:)
      type(ComplexMonomial), allocatable :: output(:)
    end function
  end interface
end module
