! ======================================================================
! Interpolate a polynomial to a single q-point and its pair,
!    integrating across all other q-points
! ======================================================================
module polynomial_interpolation_module
  use common_module
  use anharmonic_common_module
  use permutation_module
  use stress_basis_function_module
  use coupling_stress_basis_functions_module
  use polynomial_interpolator_module
  implicit none
  
  ! TODO: public/private
  !private
  
  type, extends(NoDefaultConstructor) :: SplitMonomial
    type(ComplexMonomial), allocatable :: head(:)
    type(ComplexMonomial), allocatable :: tail(:)
  end type
  
  interface SplitMonomial
    module procedure new_SplitMonomial_ComplexMonomial
  end interface
  
  interface size
    module procedure size_SplitMonomial
  end interface
contains

impure elemental function new_SplitMonomial_ComplexMonomial(monomial, &
   & anharmonic_data) result(this)
  implicit none
  
  type(ComplexMonomial), intent(in) :: monomial
  type(AnharmonicData),  intent(in) :: anharmonic_data
  type(SplitMonomial)               :: this
  
  integer :: total_size
  
  type(ComplexUnivariate), allocatable :: head(:)
  type(ComplexUnivariate), allocatable :: tail(:)
  
  type(ComplexMonomial) :: head_monomial
  
  type(ComplexUnivariate), allocatable :: univariates(:)
  
  integer, allocatable :: sizes(:)
  
  integer :: i,j,k,l,ialloc
  
  allocate( head(size(monomial)),  &
          & tail(size(monomial)),  &
          & sizes(size(monomial)), &
          & stat=ialloc); call err(ialloc)
  univariates = monomial%modes()
  do i=1,size(sizes)
    if (univariates(i)%id==univariates(i)%paired_id) then
      sizes(i) = univariates(i)%power+1
    else
      sizes(i) = (univariates(i)%power+1)*(univariates(i)%paired_power+1)
    endif
  enddo
  
  total_size = product(sizes)
  
  allocate( this%head(total_size), &
          & this%tail(total_size), &
          & stat=ialloc); call err(ialloc)
  this%head = monomial
  this%tail = monomial
  head = univariates
  tail = univariates
  l = 0
  do i=1,size(this%head)
    do j=1,size(univariates)
      k = modulo((i-1)/product(sizes(j+1:)),sizes(j))
      if (univariates(j)%id==univariates(j)%paired_id) then
        head(j)%power = k
        tail(j)%power = univariates(j)%power - head(j)%power
      else
        head(j)%power = modulo(k,univariates(j)%power+1)
        head(j)%paired_power = k/(univariates(j)%power+1)
      endif
    enddo
    head_monomial = ComplexMonomial( coefficient = monomial%coefficient, &
                                   & modes       = head                  )
    if (is_int(head_monomial%wavevector( anharmonic_data%complex_modes, &
                                       & anharmonic_data%qpoints        ))) then
      ! TODO: G-vector per subspace.
      l = l+1
      this%head(l) = ComplexMonomial( coefficient = monomial%coefficient, &
                                    & modes       = head                  )
      this%tail(l) = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                                    & modes       = tail                     )
    endif
  enddo
  this%head = this%head(:l)
  this%tail = this%tail(:l)
end function

function size_SplitMonomial(this) result(output)
  implicit none
  
  type(SplitMonomial), intent(in) :: this
  integer                         :: output
  
  output = size(this%head)
end function
end module
