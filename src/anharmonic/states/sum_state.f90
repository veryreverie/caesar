! ======================================================================
! A sum of states and coefficients.
! ======================================================================
module sum_state_module
  use common_module
  
  use subspace_state_module
  implicit none
  
  private
  
  public :: SumState
  public :: size
  public :: braket
  
  type, extends(NoDefaultConstructor) :: SumState
    integer                          :: subspace_id
    type(SubspaceState), allocatable :: states(:)
    real(dp),            allocatable :: coefficients(:)
  end type
  
  interface SumState
    module procedure new_SumState
  end interface
  
  interface size
    module procedure size_SumState
  end interface
  
  interface braket
    module procedure braket_SumStates
    module procedure braket_SumStates_ComplexMonomial
    module procedure braket_SumStates_ComplexPolynomial
  end interface
contains

! Constructor and size function
function new_SumState(subspace_id,states,coefficients) result(this)
  implicit none
  
  integer                          :: subspace_id
  type(SubspaceState), allocatable :: states(:)
  real(dp),            allocatable :: coefficients(:)
  type(SumState)                   :: this
  
  if (size(states)==0) then
    call print_line(CODE_ERROR//': No states.')
    call err()
  elseif (size(states)/=size(coefficients)) then
    call print_line(CODE_ERROR//': States and coefficients do not match.')
    call err()
  endif
  
  this%subspace_id  = subspace_id
  this%states       = states
  this%coefficients = coefficients
end function

impure elemental function size_SumState(this) result(output)
  implicit none
  
  type(SumState), intent(in) :: this
  integer                    :: output
  
  output = size(this%states)
end function

! ----------------------------------------------------------------------
! Integrals of the form <bra|ket> and <bra|potential|ket>.
! ----------------------------------------------------------------------
impure elemental function braket_SumStates(bra,ket) result(output)
  implicit none
  
  type(SumState), intent(in) :: bra
  type(SumState), intent(in) :: ket
  real(dp)                   :: output
  
  integer :: i,j
  
  output = 0.0_dp
  do i=1,size(bra)
    do j=1,size(ket)
      output = output                              &
           & + braket(bra%states(i),ket%states(j)) &
           & * bra%coefficients(i)                 &
           & * ket%coefficients(j)
    enddo
  enddo
end function

impure elemental function braket_SumStates_ComplexMonomial(bra,ket,monomial, &
   & subspace,supercell) result(output)
  implicit none
  
  type(SumState),           intent(in) :: bra
  type(SumState),           intent(in) :: ket
  type(ComplexMonomial),    intent(in) :: monomial
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  type(ComplexMonomial)                :: output
  
  type(ComplexMonomial) :: integrated_monomial
  
  integer :: i,j
  
  do i=1,size(bra)
    do j=1,size(ket)
      integrated_monomial = braket( bra%states(i),   &
                        &           ket%states(j),   &
                        &           monomial,        &
                        &           subspace,        &
                        &           supercell      ) &
                        & * bra%coefficients(i)      &
                        & * ket%coefficients(j)
      if (i==1 .and. j==1) then
        output = integrated_monomial
      else
        output%coefficient = output%coefficient &
                         & + integrated_monomial%coefficient
      endif
    enddo
  enddo
end function
  
impure elemental function braket_SumStates_ComplexPolynomial(bra,ket, &
   & polynomial,subspace,supercell) result(output)
  implicit none
  
  type(SumState),           intent(in) :: bra
  type(SumState),           intent(in) :: ket
  type(ComplexPolynomial),  intent(in) :: polynomial
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  type(ComplexPolynomial)              :: output
  
  type(ComplexPolynomial) :: integrated_polynomial
  
  integer :: i,j
  
  do i=1,size(bra)
    do j=1,size(ket)
      integrated_polynomial = braket( bra%states(i),   &
                          &           ket%states(j),   &
                          &           polynomial,      &
                          &           subspace,        &
                          &           supercell      ) &
                          & * bra%coefficients(i)      &
                          & * ket%coefficients(j)
      if (i==1 .and. j==1) then
        output = integrated_polynomial
      else
        output%terms%coefficient = output%terms%coefficient &
                               & + integrated_polynomial%terms%coefficient
      endif
    enddo
  enddo
end function
end module
