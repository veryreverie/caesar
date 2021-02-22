! ======================================================================
! Routines for finding the eigenvalues and eigenvectors of Hermitian,
!    symmetric and general matrices.
! ======================================================================
module caesar_hermitian_eigenstuff_module
  use caesar_foundations_module
  use caesar_io_module
  
  use caesar_lapack_wrapper_module
  use caesar_linear_algebra_module
  use caesar_phase_module
  use caesar_qr_decomposition_module
  implicit none
  
  private
  
  ! Classes.
  public :: SymmetricEigenstuff
  public :: HermitianEigenstuff
  public :: ComplexEigenstuff
  
  ! Diagonalisation routines.
  public :: diagonalise_symmetric
  public :: diagonalise_hermitian
  public :: diagonalise_complex
  
  ! --------------------------------------------------
  ! The eigen(values and vectors) of a matrix.
  ! --------------------------------------------------
  type SymmetricEigenstuff
    real(dp)              :: eval
    real(dp), allocatable :: evec(:)
  end type
  
  type HermitianEigenstuff
    real(dp)                 :: eval
    complex(dp), allocatable :: evec(:)
  end type
  
  type ComplexEigenstuff
    complex(dp)              :: eval
    complex(dp), allocatable :: left_evec(:)
    complex(dp), allocatable :: right_evec(:)
  end type
  
  interface diagonalise_symmetric
    ! --------------------------------------------------
    ! Calculates the eigenvalues and eigenvectors of a real, symmetric matrix.
    ! --------------------------------------------------
    module function diagonalise_symmetric_reals(input,basis) result(output) 
      real(dp),         intent(in)           :: input(:,:)
      type(RealVector), intent(in), optional :: basis(:)
      type(SymmetricEigenstuff), allocatable :: output(:)
    end function
  
    module function diagonalise_symmetric_RealMatrix(input,basis) &
       & result(output) 
      type(RealMatrix), intent(in)           :: input
      type(RealVector), intent(in), optional :: basis(:)
      type(SymmetricEigenstuff), allocatable :: output(:)
    end function
  end interface
  
  interface diagonalise_hermitian
    ! --------------------------------------------------
    ! Calculates the eigenvalues and eigenvectors of a complex, Hermitian matrix.
    ! --------------------------------------------------
    module function diagonalise_hermitian_complexes(input,basis) &
       & result(output) 
      complex(dp),         intent(in)           :: input(:,:)
      type(ComplexVector), intent(in), optional :: basis(:)
      type(HermitianEigenstuff), allocatable    :: output(:)
    end function
  
    module function diagonalise_hermitian_ComplexMatrix(input,basis) &
       & result(output) 
      type(ComplexMatrix), intent(in)           :: input
      type(ComplexVector), intent(in), optional :: basis(:)
      type(HermitianEigenstuff), allocatable    :: output(:)
    end function
  end interface
  
  interface diagonalise_complex
    ! --------------------------------------------------
    ! Calculates the eigenvalues and eigenvectors of
    !    a general square complex matrix.
    ! --------------------------------------------------
    module function diagonalise_complex_complexes(input,basis) result(output) 
      complex(dp),         intent(in)           :: input(:,:)
      type(ComplexVector), intent(in), optional :: basis(:)
      type(ComplexEigenstuff), allocatable      :: output(:)
    end function
  
    module function diagonalise_complex_ComplexMatrix(input,basis) &
       & result(output) 
      type(ComplexMatrix), intent(in)           :: input
      type(ComplexVector), intent(in), optional :: basis(:)
      type(ComplexEigenstuff), allocatable      :: output(:)
    end function
  end interface
end module
