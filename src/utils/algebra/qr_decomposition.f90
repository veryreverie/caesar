!> Provides routines for finding the QR decomposition of various matrices.
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
  
  !> The QR decomposition of a real matrix.
  type :: RealQRDecomposition
    real(dp), allocatable :: q(:,:)
    real(dp), allocatable :: r(:,:)
  end type
  
  !> The QR decomposition of a complex matrix.
  type :: ComplexQRDecomposition
    complex(dp), allocatable :: q(:,:)
    complex(dp), allocatable :: r(:,:)
  end type
  
  interface qr_decomposition
    !> Calculates the QR decomposition of a real matrix.
    module function qr_decomposition_reals(input) result(output) 
      real(dp), intent(in)      :: input(:,:)
      type(RealQRDecomposition) :: output
    end function
  
    !> Calculates the QR decomposition of a real matrix.
    module function qr_decomposition_RealMatrix(input) result(output) 
      type(RealMatrix), intent(in) :: input
      type(RealQRDecomposition)    :: output
    end function
  
    !> Calculates the QR decomposition of a complex matrix.
    module function qr_decomposition_complexes(input) result(output) 
      complex(dp), intent(in)      :: input(:,:)
      type(ComplexQRDecomposition) :: output
    end function
  
    !> Calculates the QR decomposition of a complex matrix.
    module function qr_decomposition_ComplexMatrix(input) result(output) 
      type(ComplexMatrix), intent(in) :: input
      type(ComplexQRDecomposition)    :: output
    end function
  end interface
end module
