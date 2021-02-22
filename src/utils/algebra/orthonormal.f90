! ======================================================================
! Routines for finding orthonormal bases for vector subspaces.
! ======================================================================
module caesar_orthonormal_module
  use caesar_foundations_module
  use caesar_abstract_module
  use caesar_io_module
  
  use caesar_linear_algebra_module
  use caesar_hermitian_eigenstuff_module
  use caesar_qr_decomposition_module
  implicit none
  
  private
  
  public :: orthonormal_basis
  public :: intersection_basis
  
  interface orthonormal_basis
    ! --------------------------------------------------
    ! Construct the minimal orthonormal basis which spans the given vectors.
    ! Only returns basis vectors with an L2 projection
    !    onto the input vectors of at least shortest_valid.
    ! If any basis vectors have an L2 projection onto the input vectors of between
    !    longest_invalid and shortest_valid, an error is thrown.
    ! --------------------------------------------------
    module function orthonormal_basis_RealVectors(input,shortest_valid, &
       & longest_invalid) result(output) 
      type(RealVector), intent(in)  :: input(:)
      real(dp),         intent(in)  :: shortest_valid
      real(dp),         intent(in)  :: longest_invalid
      type(RealVector), allocatable :: output(:)
    end function
  
    module function orthonormal_basis_ComplexVectors(input,shortest_valid, &
       & longest_invalid) result(output) 
      type(ComplexVector), intent(in)  :: input(:)
      real(dp),            intent(in)  :: shortest_valid
      real(dp),            intent(in)  :: longest_invalid
      type(ComplexVector), allocatable :: output(:)
    end function
  end interface
  
  interface intersection_basis
    ! ----------------------------------------------------------------------
    ! Uses the Zassenhaus algorithm, based on, QR factorisation, to find
    !    an orthonormal basis for the intersection of two vector subspaces,
    !    given an orthonormal basis for each subspace.
    ! ----------------------------------------------------------------------
    module function intersection_basis_ComplexVectors(a,b) result(output) 
      type(ComplexVector), intent(in)  :: a(:)
      type(ComplexVector), intent(in)  :: b(:)
      type(ComplexVector), allocatable :: output(:)
    end function
  end interface
end module
