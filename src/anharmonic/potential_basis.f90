! ======================================================================
! The basis functions from which the Born-Oppenheimer potential is made.
! ======================================================================
module potential_basis_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  ! A monomial, e.g.
  !    C * (u1)**a * (u2)**b * (u4)**d => coef=C, powers=[a,b,0,d]
  type Monomial
    real(dp)             :: coefficient
    integer, allocatable :: powers(:)
  end type
  
  ! A sum of monomials.
  type BasisFunction
    type(Monomial), allocatable :: monomials(:)
  end type
contains

subroutine calculate_basis_functions(coupling,no_basis_functions)
end subroutine

function simplify(potential) result(output)
  implicit none
  
  type(Monomial), intent(in)  :: potential(:)
  type(Monomial), allocatable :: output(:)
  
  integer, allocatable :: unique_id(:)
  
  integer :: i,j,id,ialloc
  
  allocate(unique_id(size(potential)), stat=ialloc); call err(ialloc)
  id = 0
  do_i : do i=1,size(potential)
    do j=1,i-1
      if (all(potential(i)%powers==potential(j)%powers)) then
        unique_id(i) = unique_id(j)
        cycle do_i
      endif
    enddo
    
    id = id + 1
    unique_id(i) = id
  enddo do_i
  
  allocate(output(id), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    output(i)%coefficient = 0
  enddo
  
  do i=1,size(potential)
    id = unique_id(i)
    output(id)%coefficient = output(id)%coefficient &
                         & + potential(i)%coefficient
    output(id)%powers = potential(i)%powers
  enddo
end function

! ----------------------------------------------------------------------
! Takes a potential, a mode id and the harmonic couplings of that mode,
!    and integrates the potential between the average of the mode.
! ----------------------------------------------------------------------
! 
! f( V(u1,u2), 1, {<u1|u^n|u1>}, {|U1>} ) -> V(u2) = average(<U1|V(u1,u2)|U1>)
!
! where |u1> is a harmonic eigenstate, and |U1> is an anharmonic eigenstate.
function integrate_over_mode_average(potential,mode,harmonic_couplings, &
   & eigenstates) result(output)
  use linear_algebra_module
  implicit none
  
  type(Monomial),   intent(in) :: potential(:)
  integer,          intent(in) :: mode
  type(RealMatrix), intent(in) :: harmonic_couplings(:)
  type(RealMatrix), intent(in) :: eigenstates
  type(Monomial), allocatable  :: output(:)
  
  integer :: i
  
  ! Copy input into output.
  output = potential
  
  ! Integrate across each monomial.
  do i=1,size(output)
    output(i)%coefficient = output(i)%coefficient                             &
                        & + trace( eigenstates                                &
                        &        * harmonic_couplings(output(i)%powers(mode)) &
                        &        * transpose(eigenstates) )                   &
                        & / size(eigenstates,1)
    output(i)%powers(mode) = 0
  enddo
  
  ! Simplify across now identical monomials.
  output = simplify(output)
end function

end module
