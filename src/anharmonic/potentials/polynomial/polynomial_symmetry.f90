! ======================================================================
! Procedures for working with symmetry in polynomial co-ordinates.
! ======================================================================
module caesar_polynomial_symmetry_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  implicit none
  
  private
  
  public :: BasisConversion
  
  type, extends(NoDefaultConstructor) :: BasisConversion
    type(ComplexMonomial), allocatable, private :: monomials_(:)
    integer,               allocatable, private :: mapping_(:)
  contains
    procedure, public :: vector_to_basis
    procedure, public :: vector_from_basis
    procedure, public :: matrix_to_basis
  end type
  
  interface BasisConversion
    module function new_BasisConversion(monomials) result(this) 
      type(ComplexMonomial), intent(in) :: monomials(:)
      type(BasisConversion)             :: this
    end function
  end interface
  
  interface
    module function vector_to_basis(this,input) result(output) 
      class(BasisConversion), intent(in) :: this
      complex(dp),            intent(in) :: input(:)
      real(dp), allocatable              :: output(:)
    end function
  end interface
  
  interface
    module function vector_from_basis(this,vector,threshold) result(output) 
      class(BasisConversion), intent(in) :: this
      real(dp),               intent(in) :: vector(:)
      real(dp),               intent(in) :: threshold
      type(ComplexPolynomial)            :: output
    end function
  end interface
  
  interface
    module function matrix_to_basis(this,matrix,lhs,rhs) result(output) 
      class(BasisConversion), intent(in) :: this
      type(ComplexMatrix),    intent(in) :: matrix
      logical,                intent(in) :: lhs
      logical,                intent(in) :: rhs
      type(RealMatrix)                   :: output
    end function
  end interface
  
  interface
    ! Find the id of the monomial conjugates,
    !    s.t. conjg(input(i)) = input(output(i)).
    module function find_monomial_conjugates(input) result(output) 
      type(ComplexMonomial), intent(in) :: input(:)
      integer, allocatable              :: output(:)
    end function
  end interface
  
  interface
    ! Convert the projection matrix from the basis of complex monomials to
    !    the basis of basis functions.
    module function project_matrix_lhs_rhs(input,conjugates) result(output) 
      complex(dp), intent(in) :: input(:,:)
      integer,     intent(in) :: conjugates(:)
      real(dp), allocatable   :: output(:,:)
    end function
  end interface
  
  interface
    module function project_matrix_lhs(input,conjugates) result(output) 
      complex(dp), intent(in) :: input(:,:)
      integer,     intent(in) :: conjugates(:)
      real(dp), allocatable   :: output(:,:)
    end function
  end interface
  
  interface
    module function project_matrix_rhs(input,conjugates) result(output) 
      complex(dp), intent(in) :: input(:,:)
      integer,     intent(in) :: conjugates(:)
      real(dp), allocatable   :: output(:,:)
    end function
  end interface
end module
