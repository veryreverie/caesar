! ======================================================================
! Various matrix tests.
! ======================================================================
module caesar_tests_module
  use caesar_foundations_module
  use caesar_io_module
  
  use caesar_linear_algebra_module
  use caesar_algebra_utils_module
  implicit none
  
  private
  
  public :: check_identity
  public :: check_symmetric
  public :: check_hermitian
  public :: check_orthogonal
  public :: check_unitary
  public :: check_orthonormal
  public :: check_real
  public :: check_int
  
  interface check_identity
    ! ----------------------------------------------------------------------
    ! Various checks for matrices.
    ! ----------------------------------------------------------------------
    ! Checks will print a warning message if they fail.
    ! matrix_name will be used to make messages more descriptive.
    ! If a logfile is present then it will always be written to.
    ! warning_threshold is the threshold at which the test is failed.
    ! ignore_threshold is the norm of the matrix below which the test is ignored.
    
    module subroutine check_identity_String(input,matrix_name,logfile, &
       & warning_threshold,ignore_threshold) 
      type(RealMatrix), intent(in)              :: input
      type(String),     intent(in)              :: matrix_name
      type(OFile),      intent(inout), optional :: logfile
      real(dp),         intent(in),    optional :: warning_threshold
      real(dp),         intent(in),    optional :: ignore_threshold
    end subroutine
  
    module subroutine check_identity_character(input,matrix_name,logfile, &
       & warning_threshold,ignore_threshold) 
      type(RealMatrix), intent(in)              :: input
      character(*),     intent(in)              :: matrix_name
      type(OFile),      intent(inout), optional :: logfile
      real(dp),         intent(in),    optional :: warning_threshold
      real(dp),         intent(in),    optional :: ignore_threshold
    end subroutine
  end interface
  
  interface check_symmetric
    module subroutine check_symmetric_String(input,matrix_name,logfile, &
       & warning_threshold,ignore_threshold) 
      type(RealMatrix), intent(in)              :: input
      type(String),     intent(in)              :: matrix_name
      type(OFile),      intent(inout), optional :: logfile
      real(dp),         intent(in),    optional :: warning_threshold
      real(dp),         intent(in),    optional :: ignore_threshold
    end subroutine
  
    module subroutine check_symmetric_character(input,matrix_name,logfile, &
       & warning_threshold,ignore_threshold) 
      type(RealMatrix), intent(in)              :: input
      character(*),     intent(in)              :: matrix_name
      type(OFile),      intent(inout), optional :: logfile
      real(dp),         intent(in),    optional :: warning_threshold
      real(dp),         intent(in),    optional :: ignore_threshold
    end subroutine
  end interface
  
  interface check_hermitian
    module subroutine check_hermitian_String(input,matrix_name,logfile, &
       & warning_threshold,ignore_threshold) 
      type(ComplexMatrix), intent(in)              :: input
      type(String),        intent(in)              :: matrix_name
      type(OFile),         intent(inout), optional :: logfile
      real(dp),            intent(in),    optional :: warning_threshold
      real(dp),            intent(in),    optional :: ignore_threshold
    end subroutine
  
    module subroutine check_hermitian_character(input,matrix_name,logfile, &
       & warning_threshold,ignore_threshold) 
      type(ComplexMatrix), intent(in)              :: input
      character(*),        intent(in)              :: matrix_name
      type(OFile),         intent(inout), optional :: logfile
      real(dp),            intent(in),    optional :: warning_threshold
      real(dp),            intent(in),    optional :: ignore_threshold
    end subroutine
  end interface
  
  interface check_orthogonal
    module subroutine check_orthogonal_String(input,matrix_name,logfile, &
       & warning_threshold) 
      type(RealMatrix), intent(in)              :: input
      type(String),     intent(in)              :: matrix_name
      type(OFile),      intent(inout), optional :: logfile
      real(dp),         intent(in),    optional :: warning_threshold
    end subroutine
  
    module subroutine check_orthogonal_character(input,matrix_name,logfile, &
       & warning_threshold) 
      type(RealMatrix), intent(in)              :: input
      character(*),     intent(in)              :: matrix_name
      type(OFile),      intent(inout), optional :: logfile
      real(dp),         intent(in),    optional :: warning_threshold
    end subroutine
  end interface
  
  interface check_unitary
    module subroutine check_unitary_String(input,matrix_name,logfile, &
       & warning_threshold) 
      type(ComplexMatrix), intent(in)              :: input
      type(String),        intent(in)              :: matrix_name
      type(OFile),         intent(inout), optional :: logfile
      real(dp),            intent(in),    optional :: warning_threshold
    end subroutine
  
    module subroutine check_unitary_character(input,matrix_name,logfile, &
       & warning_threshold) 
      type(ComplexMatrix), intent(in)              :: input
      character(*),        intent(in)              :: matrix_name
      type(OFile),         intent(inout), optional :: logfile
      real(dp),            intent(in),    optional :: warning_threshold
    end subroutine
  end interface
  
  interface check_orthonormal
    module subroutine check_orthonormal_RealMatrix_String(input,matrix_name, &
       & logfile,warning_threshold) 
      type(RealMatrix), intent(in)              :: input
      type(String),     intent(in)              :: matrix_name
      type(OFile),      intent(inout), optional :: logfile
      real(dp),         intent(in),    optional :: warning_threshold
    end subroutine
  
    module subroutine check_orthonormal_RealMatrix_character(input, &
       & matrix_name,logfile,warning_threshold) 
      type(RealMatrix), intent(in)              :: input
      character(*),     intent(in)              :: matrix_name
      type(OFile),      intent(inout), optional :: logfile
      real(dp),         intent(in),    optional :: warning_threshold
    end subroutine
  
    module subroutine check_orthonormal_ComplexMatrix_String(input, &
       & matrix_name,logfile,warning_threshold) 
      type(ComplexMatrix), intent(in)              :: input
      type(String),        intent(in)              :: matrix_name
      type(OFile),         intent(inout), optional :: logfile
      real(dp),            intent(in),    optional :: warning_threshold
    end subroutine
  
    module subroutine check_orthonormal_ComplexMatrix_character(input, &
       & matrix_name,logfile,warning_threshold) 
      type(ComplexMatrix), intent(in)              :: input
      character(*),        intent(in)              :: matrix_name
      type(OFile),         intent(inout), optional :: logfile
      real(dp),            intent(in),    optional :: warning_threshold
    end subroutine
  end interface
  
  interface check_real
    module subroutine check_real_String(input,matrix_name,logfile, &
       & warning_threshold) 
      type(ComplexMatrix), intent(in)              :: input
      type(String),        intent(in)              :: matrix_name
      type(OFile),         intent(inout), optional :: logfile
      real(dp),            intent(in),    optional :: warning_threshold
    end subroutine
  
    module subroutine check_real_character(input,matrix_name,logfile, &
       & warning_threshold) 
      type(ComplexMatrix), intent(in)              :: input
      character(*),        intent(in)              :: matrix_name
      type(OFile),         intent(inout), optional :: logfile
      real(dp),            intent(in),    optional :: warning_threshold
    end subroutine
  end interface
  
  interface check_int
    module subroutine check_int_String(input,matrix_name,logfile, &
       & warning_threshold) 
      type(RealMatrix), intent(in)              :: input
      type(String),     intent(in)              :: matrix_name
      type(OFile),      intent(inout), optional :: logfile
      real(dp),         intent(in),    optional :: warning_threshold
    end subroutine
  
    module subroutine check_int_character(input,matrix_name,logfile, &
       & warning_threshold) 
      type(RealMatrix), intent(in)              :: input
      character(*),     intent(in)              :: matrix_name
      type(OFile),      intent(inout), optional :: logfile
      real(dp),         intent(in),    optional :: warning_threshold
    end subroutine
  end interface
end module
