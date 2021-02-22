! ======================================================================
! Provides the Lanczos algorithm, which calculates
!    the lowest n eigenvalues and eigenvectors of a real symmetric matrix.
! ======================================================================
module caesar_lanczos_module
  use caesar_foundations_module
  use caesar_io_module
  
  use caesar_arpack_wrapper_module
  
  use caesar_linear_algebra_module
  use caesar_hermitian_eigenstuff_module
  implicit none
  
  private
  
  public :: lanczos
  
  interface
    module function lanczos(input,no_eigenvalues,eigenvalue_convergence) &
       & result(output) 
      type(RealMatrix), intent(in)           :: input
      integer,          intent(in)           :: no_eigenvalues
      real(dp),         intent(in)           :: eigenvalue_convergence
      type(SymmetricEigenstuff), allocatable :: output(:)
    end function
  end interface
end module
