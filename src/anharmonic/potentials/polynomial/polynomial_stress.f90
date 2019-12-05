! ======================================================================
! A polynomial representation of the stress.
! ======================================================================
module polynomial_stress_module
  use common_module
  
  use states_module
  use anharmonic_common_module
  
  use polynomial_interpolator_module
  use coupling_stress_basis_functions_module
  use interpolation_module
  implicit none
  
  private
  
  public :: startup_polynomial_stress
  
  public :: PolynomialStress
  
  type, extends(StressData) :: PolynomialStress
    integer,          private :: stress_expansion_order_
    type(RealMatrix), private :: reference_stress_
    type(CouplingStressBasisFunctions), allocatable, private :: &
       & basis_functions_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_PolynomialStress
    
    procedure, public :: zero_stress => zero_stress_PolynomialStress
    procedure, public :: add_constant => add_constant_PolynomialStress
    
    procedure, public :: stress_RealModeDisplacement => &
                       & stress_RealModeDisplacement_PolynomialStress
    procedure, public :: stress_ComplexModeDisplacement => &
                       & stress_ComplexModeDisplacement_PolynomialStress
    
    procedure, public :: braket_SubspaceState  => &
                       & braket_SubspaceState_PolynomialStress
    procedure, public :: braket_BasisState  => &
                       & braket_BasisState_PolynomialStress
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_PolynomialStress
    
    procedure, public :: simplify => &
                       & simplify_PolynomialStress
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_PolynomialStress
    
    procedure, public :: can_be_interpolated => &
                       & can_be_interpolated_PolynomialStress
    procedure, public :: calculate_interpolated_stress => &
                       & calculate_interpolated_stress_PolynomialStress
    
    procedure, public :: interpolate => interpolate_PolynomialStress
    
    procedure, public :: expansion_order => expansion_order_PolynomialStress
    
    ! I/O.
    procedure, public :: read  => read_PolynomialStress
    procedure, public :: write => write_PolynomialStress
  end type
  
  interface PolynomialStress
    module procedure new_PolynomialStress
    module procedure new_PolynomialStress_Strings
    module procedure new_PolynomialStress_StringArray
  end interface
contains

! Startup procedure.
subroutine startup_polynomial_stress()
  implicit none
  
  type(PolynomialStress) :: stress
  
  call stress%startup()
end subroutine

! Constructor.
function new_PolynomialStress(stress_expansion_order,reference_stress, &
   & basis_functions) result(this)
  implicit none
  
  integer,                            intent(in) :: stress_expansion_order
  type(RealMatrix),                   intent(in) :: reference_stress
  type(CouplingStressBasisFunctions), intent(in) :: basis_functions(:)
  type(PolynomialStress)                         :: this
  
  this%stress_expansion_order_ = stress_expansion_order
  this%reference_stress_       = reference_stress
  this%basis_functions_        = basis_functions
end function

! Type representation.
impure elemental function representation_PolynomialStress() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'polynomial'
end function

! Set the undisplaced stress to zero.
impure elemental subroutine zero_stress_PolynomialStress(this)
  implicit none
  
  class(PolynomialStress), intent(inout) :: this
  
  this%reference_stress_ = this%reference_stress_ - this%undisplaced_stress()
end subroutine

! Add a constant to the stress.
impure elemental subroutine add_constant_PolynomialStress(this,input)
  implicit none
  
  class(PolynomialStress), intent(inout) :: this
  type(RealMatrix),        intent(in)    :: input
  
  this%reference_stress_ = this%reference_stress_ + input
end subroutine

! Calculate the stress at a given displacement.
impure elemental function stress_RealModeDisplacement_PolynomialStress(this, &
   & displacement) result(output)
  implicit none
  
  class(PolynomialStress),    intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  type(RealMatrix)                       :: output
  
  if (size(this%basis_functions_)>0) then
    output = this%reference_stress_ &
         & + sum(this%basis_functions_%stress(displacement))
  else
    output = this%reference_stress_
  endif
end function

impure elemental function stress_ComplexModeDisplacement_PolynomialStress( &
   & this,displacement) result(output)
  implicit none
  
  class(PolynomialStress),       intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexMatrix)                       :: output
  
  if (size(this%basis_functions_)>0) then
    output = cmplxmat(this%reference_stress_) &
         & + sum(this%basis_functions_%stress(displacement))
  else
    output = cmplxmat(this%reference_stress_)
  endif
end function

! Integrate the stress between two states.
subroutine braket_SubspaceState_PolynomialStress(this,bra,ket,whole_subspace, &
   & anharmonic_data)
  implicit none
  
  class(PolynomialStress),  intent(inout)        :: this
  class(SubspaceState),     intent(in)           :: bra
  class(SubspaceState),     intent(in), optional :: ket
  logical,                  intent(in), optional :: whole_subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  integer :: i
  
  ! Integrate the reference stress (N.B. <i|e|j> = e<i|j> if e is a scalar.).
  this%reference_stress_ = this%reference_stress_ &
                       & * bra%inner_product(ket,anharmonic_data)
  
  ! Integrate each basis function between the bra and the ket.
  do i=1,size(this%basis_functions_)
    call this%basis_functions_(i)%braket( bra,            &
                                        & ket,            &
                                        & whole_subspace, &
                                        & anharmonic_data )
  enddo
  
  ! Simplify the stress.
  if (set_default(whole_subspace,.true.)) then
    call this%simplify()
  endif
end subroutine

! Integrate the stress between two states.
subroutine braket_BasisState_PolynomialStress(this,bra,ket,subspace, &
   & subspace_basis,whole_subspace,anharmonic_data)
  implicit none
  
  class(PolynomialStress),  intent(inout)        :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  logical,                  intent(in), optional :: whole_subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  integer :: i
  
  ! Integrate the reference stress (N.B. <i|e|j> = e<i|j> if e is a scalar.).
  this%reference_stress_ = this%reference_stress_                        &
                       & * subspace_basis%inner_product( bra,            &
                       &                                 ket,            &
                       &                                 subspace,       &
                       &                                 anharmonic_data )
  
  ! Integrate each basis function between the bra and the ket.
  do i=1,size(this%basis_functions_)
    call this%basis_functions_(i)%braket( bra,            &
                                        & ket,            &
                                        & subspace,       &
                                        & subspace_basis, &
                                        & whole_subspace, &
                                        & anharmonic_data )
  enddo
  
  ! Simplify the stress.
  if (set_default(whole_subspace,.true.)) then
    call this%simplify()
  endif
end subroutine

subroutine braket_BasisStates_PolynomialStress(this,states,thermal_energy, &
   & subspace,subspace_basis,whole_subspace,anharmonic_data)
  implicit none
  
  class(PolynomialStress),  intent(inout)        :: this
  class(BasisStates),       intent(inout)        :: states
  real(dp),                 intent(in)           :: thermal_energy
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  logical,                  intent(in), optional :: whole_subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  integer :: i
  
  ! Integrate each basis function between the bra and the ket.
  do i=1,size(this%basis_functions_)
    call this%basis_functions_(i)%braket( states,         &
                                        & thermal_energy, &
                                        & subspace,       &
                                        & subspace_basis, &
                                        & whole_subspace, &
                                        & anharmonic_data )
  enddo
  
  ! Simplify the stress.
  if (set_default(whole_subspace,.true.)) then
    call this%simplify()
  endif
end subroutine

! Identify basis functions which are constant,
!    add the constant energy to the potential's reference energy,
!    and remove the term.
! Then identify basis functions with the same coupling as a previous coupling,
!    combine the two and remove the duplicate term.
subroutine simplify_PolynomialStress(this)
  implicit none
  
  class(PolynomialStress), intent(inout) :: this
  
  logical, allocatable :: to_remove(:)
  
  integer :: i,j
  
  to_remove = [(.false., i=1, size(this%basis_functions_))]
  
  do i=1,size(this%basis_functions_)
    if (size(this%basis_functions_(i)%coupling)==0) then
      this%reference_stress_ =      &
         &   this%reference_stress_ &
         & + this%basis_functions_(i)%undisplaced_stress()
      to_remove(i) = .true.
    else
      do j=1,i-1
        if (    this%basis_functions_(i)%coupling &
           & == this%basis_functions_(j)%coupling ) then
          if (to_remove(j)) then
            call err()
          endif
          call this%basis_functions_(j)%append( &
                     & this%basis_functions_(i) )
          to_remove(i) = .true.
        endif
      enddo
    endif
  enddo
  
  ! Remove constant and duplicate terms.
  this%basis_functions_ = this%basis_functions_(filter(.not.to_remove))
end subroutine

! Calculate the thermal expectation of the stress, <stress>, for a set of
!    harmonic states.
impure elemental function harmonic_expectation_PolynomialStress(this, &
   & frequency,thermal_energy,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialStress), intent(in) :: this
  real(dp),                intent(in) :: frequency
  real(dp),                intent(in) :: thermal_energy
  type(AnharmonicData),    intent(in) :: anharmonic_data
  type(RealMatrix)                    :: output
  
  output = this%reference_stress_                                          &
       & + sum(this%basis_functions_%harmonic_expectation( frequency,      &
       &                                                   thermal_energy, &
       &                                                   anharmonic_data ))
end function

! Calculate the contribution to a given monomial from the interpolation of
!    this stress.
! The result is given as a cartesian tensor.
function can_be_interpolated_PolynomialStress(this) result(output)
  implicit none
  
  class(PolynomialStress), intent(in) :: this
  logical                             :: output
  
  output = .true.
end function

function calculate_interpolated_stress_PolynomialStress(this,             &
   & degenerate_frequency,fine_qpoints,thermal_energy,harmonic_supercell, &
   & harmonic_hessian,harmonic_min_images,subspaces,subspace_bases,       &
   & basis_states,anharmonic_min_images,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialStress),  intent(in) :: this
  real(dp),                 intent(in) :: degenerate_frequency
  type(RealVector),         intent(in) :: fine_qpoints(:)
  real(dp),                 intent(in) :: thermal_energy
  type(StructureData),      intent(in) :: harmonic_supercell
  type(CartesianHessian),   intent(in) :: harmonic_hessian
  type(MinImages),          intent(in) :: harmonic_min_images(:,:)
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  class(SubspaceBasis),     intent(in) :: subspace_bases(:)
  class(BasisStates),       intent(in) :: basis_states(:)
  type(MinImages),          intent(in) :: anharmonic_min_images(:,:)
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(RealMatrix)                     :: output
  
  type(DynamicalMatrix)          :: dynamical_matrix
  type(ComplexMode), allocatable :: qpoint_modes(:)
  type(ComplexMode), allocatable :: modes(:)
  type(RealVector),  allocatable :: qpoints(:)
  
  type(PolynomialInterpolator) :: interpolator
  
  type(DegenerateSubspace), allocatable :: fine_subspaces(:)
  type(ComplexMode),        allocatable :: subspace_modes(:)
  type(PolynomialStress)                :: fine_stress
  
  type(RealMatrix) :: stress_tensor
  
  integer :: i,j,k
  
  stress_tensor = dblemat(zeroes(3,3))
  do i=1,size(fine_qpoints)
    ! Construct normal modes from harmonic potential.
    ! Include modes from both q and -q.
    dynamical_matrix = DynamicalMatrix( fine_qpoints(i),    &
                                      & harmonic_supercell, &
                                      & harmonic_hessian,   &
                                      & harmonic_min_images )
    
    qpoint_modes = ComplexMode(dynamical_matrix, harmonic_supercell)
    qpoint_modes%id = [(j,j=1,size(qpoint_modes))]
    qpoint_modes%paired_id = [(j+size(qpoint_modes),j=1,size(qpoint_modes))]
    
    modes = [( [qpoint_modes(j),conjg(qpoint_modes(j))], &
             & j=1,                                      &
             & size(qpoint_modes)                        )]
    
    qpoints = [([fine_qpoints(i),-fine_qpoints(i)], j=1, size(qpoint_modes))]
    
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
    
    do j=1,size(fine_subspaces)
      if (fine_subspaces(j)%frequency>0) then
        ! Interpolate stress.
        subspace_modes = qpoint_modes(filter([(                   &
           & any(qpoint_modes(k)%id==fine_subspaces(j)%mode_ids), &
           & k=1,                                                 &
           & size(subspace_modes)                                 )]))
        fine_stress = generate_fine_stress( fine_subspaces(j), &
                                          & subspace_modes,    &
                                          & this,              &
                                          & interpolator       )
        
        ! Take the harmonic expectation of the stress.
        stress_tensor = stress_tensor                     &
                    & + fine_stress%harmonic_expectation( &
                    &        fine_subspaces(j)%frequency, &
                    &        thermal_energy,              &
                    &        anharmonic_data              )
      endif
    enddo
  enddo
  
  ! Normalise to be per primitive cell.
  ! N.B. at each q-point, q and -q are both considered, hence the factor of 2.
  stress_tensor = stress_tensor / (2*size(fine_qpoints))
  
  output = stress_tensor
end function

! Helper functions for interpolation.
function generate_fine_stress(subspace,modes,stress,interpolator) &
   & result(output)
  implicit none
  
  type(DegenerateSubspace),     intent(in) :: subspace
  type(ComplexMode),            intent(in) :: modes(:)
  type(PolynomialStress),       intent(in) :: stress
  type(PolynomialInterpolator), intent(in) :: interpolator
  type(PolynomialStress)                   :: output
  
  integer :: expansion_order
  integer :: power
  
  type(ComplexUnivariate) :: no_modes(0)
  
  type(ComplexMonomial), allocatable :: monomials(:)
  type(ComplexMonomial), allocatable :: old(:)
  
  type(ComplexMatrix) :: coefficients
  
  type(StressBasisFunction), allocatable :: basis_functions(:)
  type(StressBasisFunction), allocatable :: old_basis_functions(:)
  type(StressBasisFunction), allocatable :: new_basis_functions(:)
  
  integer :: i,j,k,l,ialloc
  
  expansion_order = stress%expansion_order()
  
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
    
    ! Interpolate the stress to find monomial coefficients,
    !  then convert the monomials into real basis functions.
    ! If power=paired_power, construct the u                 basis function.
    ! If power<paired_power, construct the u+u* and (u-u*)/i basis function.
    ! If power>paired_power, ignore the monomial, as it is included above.
    allocate( new_basis_functions(6*size(monomials)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(monomials)
      j = first( monomials(i)%powers()<monomials(i)%paired_powers(), &
               & default=size(monomials(i))+1                        )
      k = first( monomials(i)%powers()>monomials(i)%paired_powers(), &
               & default=size(monomials(i))+1                        )
      if (j==k) then
        coefficients = stress%interpolate(monomials(i), interpolator)
        new_basis_functions(6*l-5:6*l) = generate_stress_elements( &
                                             & [monomials(i)],     &
                                             & real(coefficients)  )
        l = l+6
      elseif (j<k) then
        coefficients = stress%interpolate(monomials(i), interpolator)
        new_basis_functions(6*l-5:6*l) = generate_stress_elements( &
                            & [monomials(i), conjg(monomials(i))], &
                            & real(coefficients)                   )
        l = l+6
        new_basis_functions(6*l-5:6*l) = generate_stress_elements(         &
           & [monomials(i), -conjg(monomials(i))]/cmplx(0.0_dp,1.0_dp,dp), &
           & aimag(coefficients)                                           )
        l = l+6
      endif
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

! Constructs the six tensor elements from an array of monomials and a matrix
!    of coefficients.
function generate_stress_elements(monomials,coefficients) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in)      :: monomials(:)
  type(RealMatrix),      intent(in)      :: coefficients
  type(StressBasisFunction), allocatable :: output(:)
  
  type(ComplexMonomial)   :: no_monomials(0)
  type(ComplexPolynomial) :: empty_polynomial
  type(ComplexPolynomial) :: elements(3,3)
  
  integer :: js(6) = [1,2,3,1,1,2]
  integer :: ks(6) = [1,2,3,2,3,3]
  
  integer :: i,j,k,ialloc
  
  empty_polynomial = ComplexPolynomial(terms=no_monomials)
  elements = empty_polynomial
  
  allocate(output(6), stat=ialloc); call err(ialloc)
  do i=1,6
    j = js(i)
    k = ks(i)
    elements(k,j) = ComplexPolynomial(monomials*coefficients%element(k,j))
    elements(j,k) = elements(k,j)
    output(i) = StressBasisFunction( elements    = elements, &
                                   & coefficient = 1.0_dp    )
    elements(k,j) = empty_polynomial
    elements(j,k) = empty_polynomial
  enddo
end function

impure elemental function interpolate_PolynomialStress(this,monomial, &
   & interpolator) result(output)
  implicit none
  
  class(PolynomialStress),      intent(in) :: this
  type(ComplexMonomial),        intent(in) :: monomial
  type(PolynomialInterpolator), intent(in) :: interpolator
  type(ComplexMatrix)                      :: output
  
  output = sum(this%basis_functions_%interpolate(monomial, interpolator))
end function

! Expansion order.
impure elemental function expansion_order_PolynomialStress(this) result(output)
  implicit none
  
  class(PolynomialStress), intent(in) :: this
  integer                             :: output
  
  output = this%stress_expansion_order_
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_PolynomialStress(this,input)
  implicit none
  
  class(PolynomialStress), intent(out) :: this
  type(String),               intent(in)  :: input(:)
  
  integer                                         :: expansion_order
  type(RealMatrix)                                :: reference_stress
  type(CouplingStressBasisFunctions), allocatable :: basis_functions(:)
  
  select type(this); type is(PolynomialStress)
    expansion_order = int(token(input(1),3))
    
    reference_stress = RealMatrix(input(3:5))
    
    basis_functions = CouplingStressBasisFunctions(split_into_sections( &
                                       & input(8:),                     &
                                       & separating_line=repeat('=',50) ))
    
    this = PolynomialStress( expansion_order,  &
                           & reference_stress, &
                           & basis_functions   )
  class default
    call err()
  end select
end subroutine

function write_PolynomialStress(this) result(output)
  implicit none
  
  class(PolynomialStress), intent(in) :: this
  type(String), allocatable              :: output(:)
  
  select type(this); type is(PolynomialStress)
    output = [ 'Expansion order: '//this%stress_expansion_order_,         &
             & str('Reference stress:'),                                  &
             & str(this%reference_stress_),                               &
             & str('Basis functions:'),                                   &
             & str(''),                                                   &
             & str(this%basis_functions_, separating_line=repeat('=',50)) ]
  class default
    call err()
  end select
end function

function new_PolynomialStress_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(PolynomialStress)   :: this
  
  call this%read(input)
end function

impure elemental function new_PolynomialStress_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(PolynomialStress)        :: this
  
  this = PolynomialStress(str(input))
end function
end module
