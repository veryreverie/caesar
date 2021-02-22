! ======================================================================
! Routines for finding the eigenvalues and eigenvectors of orthogonal
!    and unitary matrices.
! ======================================================================
! N.B. This module is untested, and is likely unstable.
! Instead of diagonalising the unitary matrix U, it is numerically superior to
!    construct and diagonalise the Hermitian matrices (U+U*) and (U-U*)/i.
module caesar_unitary_eigenstuff_module
  use caesar_foundations_module
  use caesar_abstract_module
  use caesar_io_module
  
  use caesar_mathematical_constants_module
  use caesar_linear_algebra_module
  use caesar_fraction_module
  use caesar_phase_module
  use caesar_qr_decomposition_module
  use caesar_hermitian_eigenstuff_module
  use caesar_orthonormal_module
  implicit none
  
  private
  
  public :: UnitaryEigenstuff
  public :: diagonalise_orthogonal
  public :: diagonalise_unitary
  
  type :: UnitaryEigenstuff
    type(PhaseData)          :: eval
    complex(dp), allocatable :: evec(:)
  end type
  
  interface diagonalise_orthogonal
    ! --------------------------------------------------
    ! Calculates the eigenvalues and eigenvectors of a real orthogonal matrix.
    ! --------------------------------------------------
    ! Sorts outputs in order of the phase of the eigenvalue.
    module function diagonalise_orthogonal_reals(input,order) result(output) 
      real(dp), intent(in)                 :: input(:,:)
      integer,  intent(in)                 :: order
      type(UnitaryEigenstuff), allocatable :: output(:)
    end function
  
    module function diagonalise_orthogonal_RealMatrix(input,order) &
       & result(output) 
      type(RealMatrix), intent(in)         :: input
      integer,          intent(in)         :: order
      type(UnitaryEigenstuff), allocatable :: output(:)
    end function
  end interface
  
  interface diagonalise_unitary
    ! --------------------------------------------------
    ! Calculates the eigenvalues and eigenvectors of a unitary complex matrix.
    ! --------------------------------------------------
    ! Sorts outputs in order of the phase of the eigenvalue.
    module function diagonalise_unitary_complexes(input,order) result(output) 
      complex(dp), intent(in)              :: input(:,:)
      integer,     intent(in)              :: order
      type(UnitaryEigenstuff), allocatable :: output(:)
    end function
  
    module function diagonalise_unitary_ComplexMatrix(input,order) &
       & result(output) 
      type(ComplexMatrix), intent(in)      :: input
      integer,             intent(in)      :: order
      type(UnitaryEigenstuff), allocatable :: output(:)
    end function
  end interface
end module
