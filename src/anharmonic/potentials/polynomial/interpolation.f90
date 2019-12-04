! ======================================================================
! Interpolate a polynomial to a single q-point and its pair,
!    integrating across all other q-points
! ======================================================================
module interpolation_module
  use common_module
  use anharmonic_common_module
  use permutation_module
  use stress_basis_function_module
  use coupling_stress_basis_functions_module
  use polynomial_potential_module
  use polynomial_stress_module
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
  
  integer :: power,paired_power
  
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

function interpolate_stress(stress,degenerate_frequency,harmonic_dos,        &
   & thermal_energy,harmonic_supercell,harmonic_hessian,harmonic_min_images, &
   & subspaces,subspace_bases,basis_states,anharmonic_min_images,            &
   & anharmonic_data) result(output)
  implicit none
  
  type(PolynomialStress),   intent(in) :: stress
  real(dp),                 intent(in) :: degenerate_frequency
  type(PhononDos),          intent(in) :: harmonic_dos
  real(dp),                 intent(in) :: thermal_energy
  type(StructureData),      intent(in) :: harmonic_supercell
  type(CartesianHessian),   intent(in) :: harmonic_hessian
  type(MinImages),          intent(in) :: harmonic_min_images(:,:)
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  class(SubspaceBasis),     intent(in) :: subspace_bases(:)
  class(BasisStates),       intent(in) :: basis_states(:)
  type(MinImages),          intent(in) :: anharmonic_min_images(:,:)
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(ThermodynamicData)              :: output
  
  type(DynamicalMatrix)          :: dynamical_matrix
  type(ComplexMode), allocatable :: qpoint_modes(:)
  type(ComplexMode), allocatable :: modes(:)
  type(RealVector),  allocatable :: qpoints(:)
  
  type(PolynomialInterpolator) :: interpolator
  
  type(DegenerateSubspace), allocatable :: fine_subspaces(:)
  type(PolynomialStress),   allocatable :: fine_stresses(:)
  
  type(RealMatrix) :: stress_tensor
  
  integer :: i,j
  
  stress_tensor = dblemat(zeroes(3,3))
  do i=1,size(harmonic_dos%qpoints)
    associate(qpoint => harmonic_dos%qpoints(i)%qpoint)
      ! Construct normal modes from harmonic potential.
      ! Include modes from both q and -q.
      dynamical_matrix = DynamicalMatrix( harmonic_dos%qpoints(i)%qpoint, &
                                        & harmonic_supercell,             &
                                        & harmonic_hessian,               &
                                        & harmonic_min_images             )
      
      qpoint_modes = ComplexMode(dynamical_matrix, harmonic_supercell)
      qpoint_modes%id = [(j,j=1,size(qpoint_modes))]
      qpoint_modes%paired_id = [(j+size(qpoint_modes),j=1,size(qpoint_modes))]
      
      modes = [( [qpoint_modes(j),conjg(qpoint_modes(j))], &
               & j=1,                                      &
               & size(qpoint_modes)                        )]
      
      qpoints = [([qpoint,-qpoint], j=1, size(qpoint_modes))]
      
      ! Construct the polynomial interpolator.
      interpolator = PolynomialInterpolator(                &
         & fine_modes      = modes,                         &
         & fine_qpoints    = qpoints,                       &
         & coarse_modes    = anharmonic_data%complex_modes, &
         & min_images      = anharmonic_min_images,         &
         & anharmonic_data = anharmonic_data                )
      
      ! Split modes into subspaces.
      fine_subspaces = generate_fine_subspaces( qpoint_modes,        &
                                              & degenerate_frequency )
      
      ! Construct stress basis functions.
      fine_stresses = [( generate_fine_stress( fine_subspaces(j),          &
                       &                       qpoint_modes,               &
                       &                       stress%expansion_order() ), &
                       & j=1,                                              &
                       & size(subspaces)                                   )]
      
      ! Interpolate stress to the q-point.
      call fine_stresses%add_overlap(stress,interpolator)
      
      ! Take the harmonic expectation of the stress.
      do j=1,size(fine_subspaces)
        if (fine_subspaces(j)%frequency>0) then
          stress_tensor = stress_tensor                          &
                      & + fine_stresses(j)%harmonic_expectation( &
                      &             fine_subspaces(j)%frequency, &
                      &             thermal_energy,              &
                      &             anharmonic_data              )
        endif
      enddo
    end associate
  enddo
  
  ! Normalise to be per primitive cell.
  ! N.B. at each q-point, q and -q are both considered, hence the factor of 2.
  stress_tensor = stress_tensor / (2*size(harmonic_dos%qpoints))
  
  output = harmonic_dos%thermodynamic_data(1)
  output%stress = stress_tensor
  output%enthalpy = output%energy &
                & - trace(stress_tensor)*anharmonic_data%structure%volume/3
  output%gibbs = output%free_energy &
             & - trace(stress_tensor)*anharmonic_data%structure%volume/3
end function

function generate_fine_subspaces(modes,degenerate_frequency) result(output)
  implicit none
  
  type(ComplexMode), intent(in)         :: modes(:)
  real(dp),          intent(in)         :: degenerate_frequency
  type(DegenerateSubspace), allocatable :: output(:)
  
  integer, allocatable :: splits(:)
  integer              :: no_splits
  
  integer :: i,ialloc
  
  allocate(splits(size(modes)+1), stat=ialloc); call err(ialloc)
  no_splits = 1
  splits(no_splits) = 1
  do i=2,size(modes)
    if (abs(modes(i)%frequency-modes(i-1)%frequency)>degenerate_frequency) then
      no_splits = no_splits+1
      splits(no_splits) = i
    endif
  enddo
  no_splits = no_splits+1
  splits(no_splits) = size(modes)+1
  
  allocate(output(no_splits-1), stat=ialloc); call err(ialloc)
  do i=1,no_splits-1
    associate(subspace_modes=>modes(splits(i):splits(i+1)-1))
      output(i) = DegenerateSubspace(                                  &
         & id         = i,                                             &
         & frequency  = subspace_modes(1)%frequency,                   &
         & mode_ids   = [subspace_modes%id, subspace_modes%paired_id], &
         & paired_ids = [subspace_modes%paired_id, subspace_modes%id]  )
    end associate
  enddo
end function

function generate_fine_stress(subspace,modes,expansion_order) result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  type(ComplexMode),        intent(in) :: modes(:)
  integer,                  intent(in) :: expansion_order
  type(PolynomialStress)               :: output
  
  integer :: power
  
  type(ComplexUnivariate) :: no_modes(0)
  
  type(ComplexMonomial), allocatable :: monomials(:)
  type(ComplexMonomial), allocatable :: old(:)
  
  type(ComplexPolynomial) :: polynomial
  
  type(StressBasisFunction), allocatable :: basis_functions(:)
  type(StressBasisFunction), allocatable :: old_basis_functions(:)
  type(StressBasisFunction), allocatable :: new_basis_functions(:)
  
  integer :: i,j,k,l,ialloc
  
  if (modulo(expansion_order,2)/=0) then
    call print_line(ERROR//': stress expansion order must be even.')
    call err()
  elseif (expansion_order<2) then
    call print_line(ERROR//': stress expansion order must be at least 2.')
    call err()
  endif
  
  do power=1,expansion_order/2
    ! Construct the monomial containing no modes.
    old = [ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                          & modes       = no_modes                 )]
    
    ! Loop over each mode in modes.
    ! For each mode, loop over all monomials containing the previous modes,
    !    and for each monomial append the new mode with all valid powers and
    !    paired_powers such that sum(monomial%powers()) and
    !    sum(monomial%paired_powers()) are both between 0 and power inclusive.
    do i=1,size(modes)-1
      monomials = [(                                                   &
         & (                                                           &
         &   ( ComplexMonomial(                                        &
         &        coefficient = cmplx(1.0_dp,0.0_dp,dp),               &
         &        modes       = [ old(l)%modes(),                      &
         &                        ComplexUnivariate(                   &
         &                           id           = modes(i)%id,       &
         &                           paired_id    = modes(i)%id,       &
         &                           power        = j,                 &
         &                           paired_power = k            )] ), &
         &     j=0,                                                    &
         &     power-sum(old(l)%powers()) ),                           &
         &   k=0,                                                      &
         &   power-sum(old(l)%paired_powers()) ),                      &
         & l=1,                                                        &
         & size(old)                                                   )]
      old = monomials
    enddo
    
    ! Add powers and paired_powers of the final mode such that
    !    sum(monomial%powers()) = sum(monomial%paired_powers()) = power.
    monomials = [(                                                        &
       & ComplexMonomial(                                                 &
       &    coefficient = cmplx(1.0_dp,0.0_dp,dp),                        &
       &    modes       = [                                               &
       &       old(i)%modes(),                                            &
       &       ComplexUnivariate(                                         &
       &          id           = modes(size(modes))%id,                   &
       &          paired_id    = modes(size(modes))%paired_id,            &
       &          power        = power-sum(old(i)%powers()),              &
       &          paired_power = power-sum(old(i)%paired_powers()) ) ] ), &
       & i=1,                                                             &
       & size(old)                                                        )]
    
    call monomials%simplify()
    
    ! Convert the monomials into real basis functions.
    ! If power=paired_power, construct the u        basis function.
    ! If power<paired_power, construct the u+u*     basis function.
    ! If power>paired_power, construct the (u-u*)/i basis function
    allocate( new_basis_functions(6*size(monomials)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(monomials)
      j = first( monomials(i)%powers()<monomials(i)%paired_powers(), &
               & default=size(monomials(i))+1                        )
      k = first( monomials(i)%powers()>monomials(i)%paired_powers(), &
               & default=size(monomials(i))+1                        )
      if (j==k) then
        polynomial = ComplexPolynomial(terms=[monomials(i)])
      elseif (j<k) then
        polynomial = ComplexPolynomial(terms=[monomials(i), conjg(monomials)])
        polynomial%terms%coefficient = [1,1]/sqrt(2.0_dp)
      else
        polynomial = ComplexPolynomial(terms=[monomials(i), conjg(monomials)])
        polynomial%terms%coefficient = [1,-1]/cmplx(0.0_dp,sqrt(2.0_dp),dp)
      endif
      
      ! Construct the six independent stress element basis functions.
      new_basis_functions(6*i-5:6*i) = generate_stress_elements(polynomial)
      l = l+6
    enddo
    
    ! Append the new basis functions to the basis functions
    !    from previous powers.
    if (power==1) then
      basis_functions = new_basis_functions
    else
      old_basis_functions = basis_functions
      basis_functions = [old_basis_functions, new_basis_functions]
    endif
  enddo
  
  output = PolynomialStress(                                        &
     & stress_expansion_order = expansion_order,                    &
     & reference_stress       = dblemat(zeroes(3,3)),               &
     & basis_functions        = [CouplingStressBasisFunctions(      &
     &    coupling        = SubspaceCoupling(ids=[subspace%id]),    &
     &    basis_functions = basis_functions                      )] )
end function

! Constructs the six tensor elements from a given polynomial.
function generate_stress_elements(polynomial) result(output)
  implicit none
  
  type(ComplexPolynomial), intent(in)    :: polynomial
  type(StressBasisFunction), allocatable :: output(:)
  
  type(ComplexMonomial)   :: no_monomials(0)
  type(ComplexPolynomial) :: empty_polynomial
  type(ComplexPolynomial) :: elements(3,3)
  
  integer :: i,ialloc
  
  empty_polynomial = ComplexPolynomial(terms=no_monomials)
  elements = empty_polynomial
  
  allocate(output(6), stat=ialloc); call err(ialloc)
  elements(1,1) = polynomial
  output(1) = StressBasisFunction( elements    = elements, &
                                 & coefficient = 0.0_dp    )
  elements(1,1) = empty_polynomial
  
  elements(2,2) = polynomial
  output(2) = StressBasisFunction( elements    = elements, &
                                 & coefficient = 0.0_dp    )
  elements(2,2) = empty_polynomial
  
  elements(3,3) = polynomial
  output(3) = StressBasisFunction( elements    = elements, &
                                 & coefficient = 0.0_dp    )
  elements(3,3) = empty_polynomial
  
  elements(1,2) = polynomial
  elements(2,1) = polynomial
  output(4) = StressBasisFunction( elements    = elements, &
                                 & coefficient = 0.0_dp    )
  elements(1,2) = empty_polynomial
  elements(2,1) = empty_polynomial
  
  elements(2,3) = polynomial
  elements(3,2) = polynomial
  output(5) = StressBasisFunction( elements    = elements, &
                                 & coefficient = 0.0_dp    )
  elements(2,3) = empty_polynomial
  elements(3,2) = empty_polynomial
  
  elements(3,1) = polynomial
  elements(1,3) = polynomial
  output(6) = StressBasisFunction( elements    = elements, &
                                 & coefficient = 0.0_dp    )
end function
end module
