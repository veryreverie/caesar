! ======================================================================
! Converts objects between complex and real co-ordinates.
! ======================================================================

! Complex modes with ids i and j, with i<j, are paired under the transformation
!    q -> -q if they transform as u_i -> u_j and u_j -> u_i.
! N.B. u_i = u_j*.
! These modes can then be transformed to and from the real modes i and j
!    as follows.
! The complex u_i and u_j will be written u+ and u-,
!    and the real u_i and u_j will be written c and s respectively.
! The choice of which mode is u+/c and which is u-/s is arbitrary,
!    but the choice of the mode with the smaller id is u+/c is used throughout.
! 
! u+ = (c  + i*s) /  sqrt(2)
! u- = (c  - i*s) /  sqrt(2)
! c  = (u+ +  u-) /  sqrt(2)
! s  = (u+ -  u-) / (sqrt(2)*i)
!
! u+ = (a+ib)e^{ 2pi i q.R}
! u- = (a-ib)e^{-2pi i q.R}
! c  = sqrt(2)(a*cos(2pi q.R) - b*sin(2pi q.R))
! s  = sqrt(2)(a*sin(2pi q.R) + b*cos(2pi q.R))
! 
! Under q -> -q : u+ ->   u-
!                 u- ->   u+
!                 c  ->   c
!                 s  -> - s
!
! N.B. if i=j, then the mode is real (2q=G, so e^{+iq.r}=e^{-iq.r}=cos(q.r)).
! In this case, u+ = c, and the mode has no pair (u- = s = 0).

module caesar_real_complex_conversion_module
  use caesar_utils_module
  
  use caesar_complex_mode_module
  use caesar_real_mode_module
  use caesar_complex_polynomial_module
  use caesar_real_polynomial_module
  use caesar_paired_polynomial_module
  implicit none
  
  private
  
  public :: complex_to_real
  public :: real_to_complex
  public :: basis_conversion_matrix
  public :: coefficient_conversion_matrix
  public :: PairedMonomial
  public :: ComplexMonomial
  public :: PairedPolynomial
  public :: ComplexPolynomial
  
  interface complex_to_real
    ! ----------------------------------------------------------------------
    ! Functions for converting between real and complex co-ordinates.
    ! ----------------------------------------------------------------------
    ! Construct a set of real modes from a set of complex modes, and vice versa.
    impure elemental module function complex_to_real_mode(input) &
       & result(output) 
      type(ComplexMode), intent(in) :: input
      type(RealMode)                :: output
    end function
  end interface
  
  interface real_to_complex
    impure elemental module function real_to_complex_mode(input) &
       & result(output) 
      type(RealMode), intent(in) :: input
      type(ComplexMode)          :: output
    end function
  end interface
  
  interface basis_conversion_matrix
    ! Construct the matrix M which converts from one basis set to the other,
    !    s.t. this = M * that.
    ! If include_coefficients is .true. then the coefficients of the monomials are
    !    included in the definition of the basis. If .false., then the coefficients
    !    are not included.
    ! e.g. if the real basis is A*u+ and B*u-, and the complex basis is C*c and
    !    D*s, where c=(u+)+i(u-) and s=(u+)-i(u-), then:
    !    - if .false. then the conversion matrix is (c) = ( 1  i) . (u+)
    !                                               (s)   ( 1 -i)   (u-)
    !    - if .true. then the conversion matrix is (Cc) = ( C/A  Ci/A) . (Au+)
    !                                              (Ds)   ( D/B -Di/B)   (Bu-)
    module function basis_conversion_matrix_ComplexMonomials_RealMonomials(this,that,include_coefficients) result(output) 
      type(ComplexMonomial), intent(in) :: this(:)
      type(RealMonomial),    intent(in) :: that(:)
      logical,               intent(in) :: include_coefficients
      type(ComplexMatrix)               :: output
    end function
  
    module function basis_conversion_matrix_RealMonomials_ComplexMonomials(this,that,include_coefficients) result(output) 
      type(RealMonomial),    intent(in) :: this(:)
      type(ComplexMonomial), intent(in) :: that(:)
      logical,               intent(in) :: include_coefficients
      type(ComplexMatrix)               :: output
    end function
  end interface
  
  interface coefficient_conversion_matrix
    ! ----------------------------------------------------------------------
    ! Conversion matrices for coefficients rather than basis functions.
    ! ----------------------------------------------------------------------
    ! Coefficients transform in the opposite manner to basis functions.
    ! If f(x) = sum_i a_ig_i(x) = sum_j b_jh_j(x) then if
    !    f_i(x) = sum_j C_ij h_j(x)
    !    h_j(x) = sum_i D_ji g_j(x)
    ! where sum_j C_ij D_jk is the identity, then
    !    a_i    = sum_j D_ji b_j
    !    b_j    = sum_j C_ij a_i
    ! So the conversion matrices are D^T and C^T rather than C and D.
    module function coefficient_conversion_matrix_ComplexMonomials_RealMonomials(this,that,include_coefficients) result(output) 
      type(ComplexMonomial), intent(in) :: this(:)
      type(RealMonomial),    intent(in) :: that(:)
      logical,               intent(in) :: include_coefficients
      type(ComplexMatrix)               :: output
    end function
  
    module function coefficient_conversion_matrix_RealMonomials_ComplexMonomials(this,that,include_coefficients) result(output) 
      type(RealMonomial),    intent(in) :: this(:)
      type(ComplexMonomial), intent(in) :: that(:)
      logical,               intent(in) :: include_coefficients
      type(ComplexMatrix)               :: output
    end function
  end interface
  
  interface element
    ! ----------------------------------------------------------------------
    ! Private helper functions.
    ! ----------------------------------------------------------------------
    
    ! Calculate an element for conversion_matrix.
    !
    ! Calculates the coefficient of 'that' in the expansion of 'this'.
    ! i.e. 'this' = sum[ element(this,that) * 'that' ].
    impure elemental module function element_ComplexMonomial_RealMonomial(this,that,include_coefficients) result(output) 
      type(ComplexMonomial), intent(in) :: this
      type(RealMonomial),    intent(in) :: that
      logical,               intent(in) :: include_coefficients
      complex(dp)                       :: output
    end function
  
    ! Calculate an element for conversion_matrix.
    !
    ! Calculates the coefficient of 'that' in the expansion of 'this'.
    ! i.e. 'this' = sum[ element(this,that) * 'that' ].
    impure elemental module function element_RealMonomial_ComplexMonomial(this,that,include_coefficients) result(output) 
      type(RealMonomial),    intent(in) :: this
      type(ComplexMonomial), intent(in) :: that
      logical,               intent(in) :: include_coefficients
      complex(dp)                       :: output
    end function
  end interface
  
  interface PairedMonomial
    ! Convert between ComplexMonomial and PairedMonomial arrays.
    module function new_PairedMonomials_ComplexMonomials(input) result(output) 
      type(ComplexMonomial), intent(in) :: input(:)
      type(PairedMonomial), allocatable :: output(:)
    end function
  end interface
  
  interface ComplexMonomial
    module function new_ComplexMonomials_PairedMonomials(input) result(output) 
      type(PairedMonomial), intent(in)   :: input(:)
      type(ComplexMonomial), allocatable :: output(:)
    end function
  end interface
  
  interface PairedPolynomial
    module function new_PairedPolynomial_ComplexPolynomial(input) &
       & result(output) 
      type(ComplexPolynomial), intent(in) :: input
      type(PairedPolynomial)              :: output
    end function
  end interface
  
  interface ComplexPolynomial
    module function new_ComplexPolynomial_PairedPolynomial(input) &
       & result(output) 
      type(PairedPolynomial), intent(in) :: input
      type(ComplexPolynomial)            :: output
    end function
  end interface
end module
