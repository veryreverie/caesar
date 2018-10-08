! ======================================================================
! Routines for finding the eigenvalues and eigenvectors of orthogonal
!    and unitary matrices.
! ======================================================================
! N.B. This module is untested, and is likely unstable.
! Instead of diagonalising the unitary matrix U, it is numerically superior to
!    construct and diagonalise the Hermitian matrices (U+U*) and (U-U*)/i.
module unitary_eigenstuff_module
  use precision_module
  use abstract_module
  use io_module
  
  use mathematical_constants_module
  use linear_algebra_module
  use fraction_module
  use phase_module
  use qr_decomposition_module
  use hermitian_eigenstuff_module
  use orthonormal_module
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
    module procedure diagonalise_orthogonal_reals
    module procedure diagonalise_orthogonal_RealMatrix
  end interface
  
  interface diagonalise_unitary
    module procedure diagonalise_unitary_complexes
    module procedure diagonalise_unitary_ComplexMatrix
  end interface
contains

! --------------------------------------------------
! Calculates the eigenvalues and eigenvectors of a real orthogonal matrix.
! --------------------------------------------------
! Sorts outputs in order of the phase of the eigenvalue.
function diagonalise_orthogonal_reals(input,order) result(output)
  implicit none
  
  real(dp), intent(in)                 :: input(:,:)
  integer,  intent(in)                 :: order
  type(UnitaryEigenstuff), allocatable :: output(:)
  
  type(RealMatrix)    :: cos_matrix
  type(ComplexMatrix) :: sin_matrix
  
  type(SymmetricEigenstuff), allocatable :: cos_estuff(:)
  type(HermitianEigenstuff), allocatable :: sin_estuff(:)
  
  real(dp) :: phase
  
  integer, allocatable :: ids(:)
  
  integer :: i,ialloc
  
  integer, allocatable :: cs(:)
  
  type(ComplexVector), allocatable :: degenerate_cos_basis(:)
  
  integer :: c
  integer :: s
  
  if (size(input,1)/=size(input,2)) then
    call print_line(CODE_ERROR//': Trying to diagonalise non-square &
       &orthogonal matrix.')
    call err()
  endif
  
  allocate( cs(size(input,1)),        &
          & output(size(input,1)),    &
          & stat=ialloc); call err(ialloc)
  
  ! M is orthogonal, so M.x=e^(2*pi*i*phase/n).x, where x is an eigenvector
  !    of M, n is the order of M (s.t. M^n=I), and f is an integer.
  
  ! C = (M+M^T)/2 is a symmetric matrix. C.x = cos(2*pi*phase/n).x .
  ! cos(2*pi*c/n) = cos(2*pi*phase/n).
  cos_matrix = (input + transpose(input))/2.0_dp
  
  ! S = (M-m^T)/2i is a Hermitian matrix. S.x = sin(2*pi*phase/n).x .
  ! sin(2*pi*s/n) = sin(2*pi*phase/n).
  sin_matrix = (input - transpose(input))/cmplx(0.0_dp,2.0_dp,dp)
  
  ! Find the eigenstuff of C, and use it to calculate {c}.
  ! N.B. phase may be c or order - c.
  ! There may also be additional degeneracy as a result.
  cos_estuff = diagonalise_symmetric(cos_matrix)
  do i=1,size(cos_estuff)
    phase = acos(cos_estuff(i)%eval)*order/(2*PI)
    c = nint(phase)
    if (abs(phase-c)>0.1_dp) then
      call err()
    endif
    cs(i) = modulo(c,order)
  enddo
  
  ! Use S to lift the degeneracies of C, and to distinguish phase=c
  !                                                     and phase=order-c.
  do c=0,order/2
    ! Get the set of eigenvectors with C eigenvalue cos(2*pi*c/order).
    ids = filter(cs==c)
    if (size(ids)==0) then
      cycle
    endif
    
    ! Construct the matrix S in this basis, and diagonalise it.
    ! N.B. S is Hermitian and purely imaginary.
    allocate(degenerate_cos_basis(size(ids)), stat=ialloc); call err(ialloc)
    do i=1,size(ids)
      degenerate_cos_basis(i) = cmplx(cos_estuff(ids(i))%evec, 0.0_dp, dp)
    enddo
    sin_estuff = diagonalise_hermitian(sin_matrix,degenerate_cos_basis)
    deallocate(degenerate_cos_basis, stat=ialloc); call err(ialloc)
    
    ! Use the eigenbases of C and S to calculate the eigenbasis of the
    !    original matrix.
    do i=1,size(sin_estuff)
      phase = asin(sin_estuff(i)%eval)*order/(2*PI)
      s = nint(phase)
      if (abs(phase-s)>0.1_dp) then
        call err()
      endif
      s = modulo(s,order)
      
      if (s<=order/2) then
        ! 2*pi*s/order <= pi, so phase=c.
        output(ids(i))%eval = PhaseData(frac(c,order))
      else
        ! 2*pi*s/order > pi, so phase=order-c.
        output(ids(i))%eval = PhaseData(frac(order-c,order))
      endif
      
      output(ids(i))%evec = sin_estuff(i)%evec
    enddo
  enddo
end function

function diagonalise_orthogonal_RealMatrix(input,order) result(output)
  implicit none
  
  type(RealMatrix), intent(in)         :: input
  integer,          intent(in)         :: order
  type(UnitaryEigenstuff), allocatable :: output(:)
  
  output = diagonalise_orthogonal(dble(input),order)
end function

! --------------------------------------------------
! Calculates the eigenvalues and eigenvectors of a unitary complex matrix.
! --------------------------------------------------
! Sorts outputs in order of the phase of the eigenvalue.
function diagonalise_unitary_complexes(input,order) result(output)
  implicit none
  
  complex(dp), intent(in)              :: input(:,:)
  integer,     intent(in)              :: order
  type(UnitaryEigenstuff), allocatable :: output(:)
  
  type(ComplexEigenstuff), allocatable :: estuff(:)
  
  integer,             allocatable :: degenerate_modes(:)
  type(ComplexVector), allocatable :: degenerate_evecs(:)
  
  integer :: i,j,ialloc
  
  ! Diagonalise matrix.
  estuff = diagonalise_complex(input)
  
  ! Convert eigenvalues to exact representation.
  allocate(output(size(estuff)), stat=ialloc); call err(ialloc)
  do i=1,size(estuff)
    output(i)%eval = calculate_phase(estuff(i)%eval, order)
    output(i)%evec = estuff(i)%right_evec
  enddo
  
  ! Sort output in ascending order of phase.
  output = output(sort(output,compare_phases))
  
  ! Orthonormalise degenerate eigenvectors.
  do i=1,size(output)
    degenerate_modes = filter(output%eval==output(i)%eval)
    if (degenerate_modes(1)<i) then
      ! This mode has already been processed.
      cycle
    elseif (size(degenerate_modes)==1) then
      ! This mode is not degenerate with any other mode.
      cycle
    endif
    
    allocate( degenerate_evecs(size(degenerate_modes)), &
            & stat=ialloc); call err(ialloc)
    do j=1,size(degenerate_modes)
      degenerate_evecs(j) = vec(output(degenerate_modes(j))%evec)
    enddo
    degenerate_evecs = orthonormal_basis(degenerate_evecs,1e-3_dp,1e-10_dp)
    do j=1,size(degenerate_modes)
      output(degenerate_modes(j))%evec = cmplx(degenerate_evecs(j))
    enddo
    deallocate(degenerate_evecs, stat=ialloc); call err(ialloc)
  enddo
contains
  ! A lambda for ordering phases.
  ! Captures nothing.
  function compare_phases(this,that) result(output)
    implicit none
    
    class(*), intent(in) :: this
    class(*), intent(in) :: that
    logical              :: output
    
    select type(this); type is(UnitaryEigenstuff)
      select type(that); type is(UnitaryEigenstuff)
        output = this%eval%fraction < that%eval%fraction
      end select
    end select
  end function
end function

function diagonalise_unitary_ComplexMatrix(input,order) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in)      :: input
  integer,             intent(in)      :: order
  type(UnitaryEigenstuff), allocatable :: output(:)
  
  output = diagonalise_unitary(cmplx(input), order)
end function
end module
