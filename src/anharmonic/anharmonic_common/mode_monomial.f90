! ======================================================================
! A product of modes, without specifying real or complex co-ordinates.
! Treats products of modes as non-commutative, e.g. u1*u2 /= u2*u1, in order
!    to calculate symmetries.
! ======================================================================
module mode_monomial_module
  use common_module
  
  use subspace_monomial_module
  use mode_coupling_module
  implicit none
  
  private
  
  public :: ModeMonomial
  public :: size
  public :: generate_mode_monomials
  public :: operator(//)
  public :: operator(==)
  public :: operator(/=)
  public :: ComplexMonomial
  public :: RealMonomial
  public :: ModeCoupling
  
  ! A list of ids of modes which are coupled.
  type :: ModeMonomial
    integer, allocatable :: mode_ids(:)
    ! Whether or not the complex basis functions corresponding to this coupling
    !    are part of the Hamiltonian. i.e. they conserve Bloch momentum.
    logical :: conserves_momentum
    ! Whether or not the complex basis functions corresponding to this coupling
    !    are part of the VSCF Hamiltonian. i.e. they conserve Bloch momentum
    !    within every coupled subspace.
    logical :: conserves_vscf
  end type
  
  interface size
    module procedure size_ModeMonomial
  end interface
  
  interface operator(//)
    module procedure concatenate_ModeMonomial_integer
  end interface
  
  interface operator(==)
    module procedure equality_ModeMonomial_ModeMonomial
  end interface
  
  interface operator(/=)
    module procedure non_equality_ModeMonomial_ModeMonomial
  end interface
  
  ! Construct the ComplexMonomial, RealMonomial or ModeCoupling
  !    corresponding to this ModeMonomial.
  interface ComplexMonomial
    module procedure new_ComplexMonomial_ModeMonomial
  end interface
  
  interface RealMonomial
    module procedure new_RealMonomial_ModeMonomial
  end interface
  
  interface ModeCoupling
    module procedure new_ModeCoupling_ModeMonomial
  end interface
contains

! size() functions.
function size_ModeMonomial(this) result(output)
  implicit none
  
  type(ModeMonomial), intent(in) :: this
  integer                        :: output
  
  output = size(this%mode_ids)
end function

! Append an id to the coupling.
function concatenate_ModeMonomial_integer(this,mode_id) result(output)
  implicit none
  
  type(ModeMonomial), intent(in) :: this
  integer,            intent(in) :: mode_id
  type(ModeMonomial)             :: output
  
  output = ModeMonomial( mode_ids           = [this%mode_ids, mode_id], &
                       & conserves_momentum = this%conserves_momentum,  &
                       & conserves_vscf     = this%conserves_vscf       )
end function

! Compare couplings.
impure elemental function equality_ModeMonomial_ModeMonomial(this,that) &
   & result(output)
  implicit none
  
  type(ModeMonomial), intent(in) :: this
  type(ModeMonomial), intent(in) :: that
  logical                        :: output
  
  output = all(this%mode_ids==that%mode_ids)
end function

impure elemental function non_equality_ModeMonomial_ModeMonomial(this,that) &
   & result(output)
  implicit none
  
  type(ModeMonomial), intent(in) :: this
  type(ModeMonomial), intent(in) :: that
  logical                        :: output
  
  output = .not. this==that
end function

! ----------------------------------------------------------------------
! Recursively generates all sets of mode monomials within a set of coupled
!    subspaces.
! Only returns couplings with sum(q)=0, modulo G-vectors.
! ----------------------------------------------------------------------
function generate_mode_monomials(subspace_monomial,subspaces,modes,qpoints) &
   & result(output)
  implicit none
  
  type(SubspaceMonomial),   intent(in) :: subspace_monomial
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  type(ComplexMode),        intent(in) :: modes(:)
  type(QpointData),         intent(in) :: qpoints(:)
  type(ModeMonomial), allocatable      :: output(:)
  
  type(ModeMonomial)   :: root_monomial
  type(FractionVector) :: root_sum_q
  
  type(DegenerateSubspace), allocatable :: coupled_subspaces(:)
  
  ! List the subspaces which are coupled together in this coupling.
  coupled_subspaces = subspace_monomial%coupled_subspaces(subspaces)
  
  root_monomial = ModeMonomial( mode_ids           = [integer::], &
                              & conserves_momentum = .true.,      &
                              & conserves_vscf     = .true.       )
  root_sum_q = fracvec(zeroes(3))
  
  ! Generate all the sets of coupled modes within the coupled subspaces.
  output = generate_mode_monomials_helper( coupled_subspaces,        &
                                         & subspace_monomial%powers, &
                                         & modes,                    &
                                         & qpoints,                  &
                                         & root_monomial,            &
                                         & root_sum_q                )
end function

recursive function generate_mode_monomials_helper(subspaces,powers, &
   & normal_modes,qpoints,mode_monomial,sum_q) result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  integer,                  intent(in) :: powers(:)
  type(ComplexMode),        intent(in) :: normal_modes(:)
  type(QpointData),         intent(in) :: qpoints(:)
  type(ModeMonomial),       intent(in) :: mode_monomial
  type(FractionVector),     intent(in) :: sum_q
  type(ModeMonomial), allocatable      :: output(:)
  
  integer :: first_subspace
  
  type(QpointData), allocatable :: subspace_qpoints(:)
  
  integer, allocatable :: new_powers(:)
  type(ModeMonomial)   :: new_mode_monomial
  type(FractionVector) :: new_sum_q
  
  logical :: last_mode_in_coupling
  
  integer :: i,j
  
  i = first(powers/=0, default=0)
  
  if (i==0) then
    ! There is nothing else to append. Check that the sum across q-points of
    !    the mode coupling is zero, and return the mode couplings.
    output = [mode_monomial]
    if (.not. is_int(sum_q)) then
      output%conserves_momentum = .false.
      output%conserves_vscf     = .false.
    endif
  else
    ! If vscf_basis_functions_only is true, then only mode couplings which have
    !    sum(q)=0 for all degenerate subspaces are allowed.
    last_mode_in_coupling = .false.
    if (i==size(subspaces)) then
      last_mode_in_coupling = .true.
    elseif (subspaces(i+1)%id/=subspaces(i)%id) then
      last_mode_in_coupling = .true.
    endif
    
    ! Update the array of powers.
    new_powers = powers
    new_powers(i) = new_powers(i)-1
    
    ! Loop over modes in this subspaces, recursively calling this function for
    !    each in turn.
    subspace_qpoints = subspaces(i)%qpoints(normal_modes,qpoints)
    output = [ModeMonomial::]
    do j=1,size(subspaces(i))
      new_mode_monomial = mode_monomial//subspaces(i)%mode_ids(j)
      new_sum_q = sum_q + subspace_qpoints(j)%qpoint
      if (last_mode_in_coupling .and. .not. is_int(new_sum_q)) then
        new_mode_monomial%conserves_vscf = .false.
      endif
      output = [ output,                                             &
             &   generate_mode_monomials_helper( subspaces,          &
             &                                   new_powers,         &
             &                                   normal_modes,       &
             &                                   qpoints,            &
             &                                   new_mode_monomial,  &
             &                                   new_sum_q          )]
    enddo
  endif
end function

! ----------------------------------------------------------------------
! Construct a ComplexMonomial, RealMonomial or ModeCoupling
!    corresponding to this ModeMonomial.
! ----------------------------------------------------------------------
function new_ComplexMonomial_ModeMonomial(this,complex_modes) result(output)
  implicit none
  
  type(ModeMonomial), intent(in) :: this
  type(ComplexMode),  intent(in) :: complex_modes(:)
  type(ComplexMonomial)          :: output
  
  integer, allocatable :: mode_ids(:)
  
  type(ComplexMode) :: mode
  integer           :: power
  
  integer :: i,ialloc
  
  mode_ids = this%mode_ids
  mode_ids = mode_ids(set(mode_ids))
  mode_ids = mode_ids(sort(mode_ids))
  output = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                          & modes       = [ComplexUnivariate::]    )
  do i=1,size(mode_ids)
    mode = complex_modes(first(complex_modes%id==mode_ids(i)))
    power = count(this%mode_ids==mode_ids(i))
    output = output * ComplexUnivariate(mode=mode, power=power)
  enddo
end function

function new_RealMonomial_ModeMonomial(this,real_modes) result(output)
  implicit none
  
  type(ModeMonomial), intent(in) :: this
  type(RealMode),     intent(in) :: real_modes(:)
  type(RealMonomial)             :: output
  
  integer, allocatable :: mode_ids(:)
  
  type(RealMode) :: mode
  
  integer :: id
  integer :: paired_id
  integer :: power
  
  integer :: i,ialloc
  
  mode_ids = this%mode_ids
  mode_ids = mode_ids(set(mode_ids))
  mode_ids = mode_ids(sort(mode_ids))
  output = RealMonomial( coefficient = 1.0_dp,            &
                       & modes       = [RealUnivariate::] )
  do i=1,size(mode_ids)
    mode = real_modes(first(real_modes%id==mode_ids(i)))
    power = count(this%mode_ids==mode_ids(i))
    output = output * RealUnivariate(mode=mode, power=power)
  enddo
end function

! e.g. if this has ids [2,1,2,3,1], the corresponding ModeCoupling
!    has ids [1,2,3].
impure elemental function new_ModeCoupling_ModeMonomial(this) result(output)
  implicit none
  
  class(ModeMonomial), intent(in) :: this
  type(ModeCoupling)              :: output
  
  integer, allocatable :: mode_ids(:)
  
  mode_ids = this%mode_ids
  mode_ids = mode_ids(set(mode_ids))
  mode_ids = mode_ids(sort(mode_ids))
  output = ModeCoupling(mode_ids)
end function
end module
