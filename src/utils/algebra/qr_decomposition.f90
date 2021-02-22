! ======================================================================
! Routines for finding the QR decomposition of various matrices.
! ======================================================================
module caesar_qr_decomposition_module
  use caesar_foundations_module
  use caesar_abstract_module
  use caesar_io_module
  
  use caesar_lapack_wrapper_module
  use caesar_linear_algebra_module
  implicit none
  
  private
  
  public :: RealQRDecomposition
  public :: ComplexQRDecomposition
  public :: qr_decomposition
  public :: determinant_qr
  
  type :: RealQRDecomposition
    real(dp), allocatable :: q(:,:)
    real(dp), allocatable :: r(:,:)
  end type
  
  type :: ComplexQRDecomposition
    complex(dp), allocatable :: q(:,:)
    complex(dp), allocatable :: r(:,:)
  end type
  
  interface qr_decomposition
    ! --------------------------------------------------
    ! Calculates the QR decomposition of a real matrix.
    ! --------------------------------------------------
    module function qr_decomposition_reals(input) result(output) 
      real(dp), intent(in)      :: input(:,:)
      type(RealQRDecomposition) :: output
    end function
  
    module function qr_decomposition_RealMatrix(input) result(output) 
      type(RealMatrix), intent(in) :: input
      type(RealQRDecomposition)    :: output
    end function
  
    ! --------------------------------------------------
    ! Calculates the QR decomposition of a complex matrix.
    ! --------------------------------------------------
    module function qr_decomposition_complexes(input) result(output) 
      complex(dp), intent(in)      :: input(:,:)
      type(ComplexQRDecomposition) :: output
    end function
  
    module function qr_decomposition_ComplexMatrix(input) result(output) 
      type(ComplexMatrix), intent(in) :: input
      type(ComplexQRDecomposition)    :: output
    end function
  end interface
  
  interface determinant_qr
    ! --------------------------------------------------
    ! Calculate the determinant of a real square matrix using QR decomposition.
    ! --------------------------------------------------
    module function determinant_qr_reals(input) result(output) 
      real(dp), intent(in) :: input(:,:)
      real(dp)             :: output
    end function
  
    module function determinant_qr_RealMatrix(input) result(output) 
      type(RealMatrix), intent(in) :: input
      real(dp)                     :: output
    end function
  end interface
end module
