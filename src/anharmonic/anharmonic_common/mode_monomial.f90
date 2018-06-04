! ======================================================================
! A product of modes, without specifying real or complex co-ordinates.
! Treats products of modes as non-commutative, e.g. u1*u2 /= u2*u1, in order
!    to calculate symmetries.
! ======================================================================
module mode_monomial_module
  use common_module
  
  use degeneracy_module
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
  
  ! A list of ids of modes which are coupled.
  type :: ModeMonomial
    integer, allocatable :: ids(:)
    ! Whether or not the complex basis functions corresponding to this coupling
    !    are part of the Hamiltonian. i.e. they conserve Bloch momentum.
    logical :: conserves_momentum
    ! Whether or not the complex basis functions corresponding to this coupling
    !    are part of the VSCF Hamiltonian. i.e. they conserve Bloch momentum
    !    within every coupled subspace.
    logical :: conserves_vscf
  contains
    ! Construct the ComplexMonomial, RealMonomial or ModeCoupling
    !    corresponding to this ModeMonomial.
    procedure, public :: complex_monomial
    procedure, public :: real_monomial
    procedure, public :: mode_coupling
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
contains

! size() functions.
function size_ModeMonomial(this) result(output)
  implicit none
  
  type(ModeMonomial), intent(in) :: this
  integer                        :: output
  
  output = size(this%ids)
end function

! Append an id to the coupling.
function concatenate_ModeMonomial_integer(this,id) result(output)
  implicit none
  
  type(ModeMonomial), intent(in) :: this
  integer,            intent(in) :: id
  type(ModeMonomial)             :: output
  
  output = ModeMonomial( ids                = [this%ids, id],          &
                       & conserves_momentum = this%conserves_momentum, &
                       & conserves_vscf     = this%conserves_vscf)
end function

! Compare couplings.
impure elemental function equality_ModeMonomial_ModeMonomial(this,that) &
   & result(output)
  implicit none
  
  type(ModeMonomial), intent(in) :: this
  type(ModeMonomial), intent(in) :: that
  logical                        :: output
  
  output = all(this%ids==that%ids)
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
function generate_mode_monomials(coupling,subspaces,normal_modes,qpoints) &
   & result(output)
  implicit none
  
  type(SubspaceMonomial),   intent(in) :: coupling
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  type(ComplexMode),        intent(in) :: normal_modes(:)
  type(QpointData),         intent(in) :: qpoints(:)
  type(ModeMonomial), allocatable     :: output(:)
  
  type(DegenerateSubspace), allocatable :: coupled_subspaces(:)
  
  ! List the subspaces which are coupled together in this coupling.
  coupled_subspaces = coupling%coupled_subspaces(subspaces)
  
  ! Generate all the sets of coupled modes within the coupled subspaces.
  output = generate_mode_monomials_helper( coupled_subspaces, &
                                         & normal_modes,      &
                                         & qpoints)
end function

recursive function generate_mode_monomials_helper(coupled_subspaces, &
   & normal_modes,qpoints,mode_monomial_in,sum_q_in) result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in)        :: coupled_subspaces(:)
  type(ComplexMode),     intent(in)           :: normal_modes(:)
  type(QpointData),      intent(in)           :: qpoints(:)
  type(ModeMonomial),    intent(in), optional :: mode_monomial_in
  type(FractionVector),  intent(in), optional :: sum_q_in
  type(ModeMonomial), allocatable             :: output(:)
  
  type(QpointData), allocatable :: subspace_qpoints(:)
  
  type(ModeMonomial)   :: mode_monomial
  type(FractionVector) :: sum_q
  
  type(ModeMonomial)   :: mode_monomial_out
  type(FractionVector) :: sum_q_out
  
  logical :: last_mode_in_coupling
  
  integer :: i
  
  if (present(mode_monomial_in) .neqv. present(sum_q_in)) then
    call print_line(CODE_ERROR//': generate_mode_monomial must be called with &
       &all optional arguments or none.')
    call err()
  endif
  
  if (present(mode_monomial_in)) then
    mode_monomial = mode_monomial_in
    sum_q = sum_q_in
  else
    mode_monomial = ModeMonomial( ids                = [integer::], &
                                & conserves_momentum = .true.,      &
                                & conserves_vscf     = .true.)
    sum_q = fracvec(zeroes(3))
  endif
  
  if (size(coupled_subspaces)==0) then
    ! There is nothing else to append. Check that the sum across q-points of
    !    the mode coupling is zero, and return the mode couplings.
    if (.not. is_int(sum_q_in)) then
      mode_monomial%conserves_momentum = .false.
      mode_monomial%conserves_vscf     = .false.
    endif
    output = [mode_monomial]
  else
    ! If vscf_basis_functions_only is true, then only mode couplings which have
    !    sum(q)=0 for all degenerate subspaces are allowed.
    last_mode_in_coupling = .false.
    if (size(coupled_subspaces)==1) then
      last_mode_in_coupling = .true.
    elseif (coupled_subspaces(2)%id/=coupled_subspaces(1)%id) then
      last_mode_in_coupling = .true.
    endif
    
    ! Loop over modes in this subspaces, recursively calling this function for
    !    each in turn.
    subspace_qpoints = coupled_subspaces(1)%qpoints(normal_modes,qpoints)
    output = [ModeMonomial::]
    do i=1,size(coupled_subspaces(1))
      mode_monomial_out = mode_monomial//coupled_subspaces(1)%mode_ids(i)
      sum_q_out = sum_q + subspace_qpoints(i)%qpoint
      if (last_mode_in_coupling .and. .not. is_int(sum_q_out)) then
        mode_monomial_out%conserves_vscf = .false.
      endif
      output = [ output,                                                &
             &   generate_mode_monomials_helper( coupled_subspaces(2:), &
             &                                   normal_modes,          &
             &                                   qpoints,               &
             &                                   mode_monomial_out,     &
             &                                   sum_q_out)             &
             & ]
    enddo
  endif
end function

! ----------------------------------------------------------------------
! Construct a ComplexMonomial, RealMonomial or ModeCoupling
!    corresponding to this ModeMonomial.
! ----------------------------------------------------------------------
function complex_monomial(this,complex_modes) result(output)
  implicit none
  
  class(ModeMonomial), intent(in) :: this
  type(ComplexMode),   intent(in) :: complex_modes(:)
  type(ComplexMonomial)           :: output
  
  integer, allocatable :: mode_ids(:)
  
  integer :: id
  integer :: paired_id
  integer :: power
  
  integer :: i,ialloc
  
  mode_ids = this%ids
  mode_ids = mode_ids(set(mode_ids))
  mode_ids = mode_ids(sort(mode_ids))
  output%coefficient = 1
  allocate(output%modes(size(mode_ids)), stat=ialloc); call err(ialloc)
  do i=1,size(mode_ids)
    id = mode_ids(i)
    paired_id = complex_modes(first(complex_modes%id==id))%paired_id
    power = count(this%ids==id)
    output%modes(i) = ComplexUnivariate( id        = id,        &
                                       & paired_id = paired_id, &
                                       & power     = power)
  enddo
end function

function real_monomial(this,real_modes) result(output)
  implicit none
  
  class(ModeMonomial), intent(in) :: this
  type(RealMode),      intent(in) :: real_modes(:)
  type(RealMonomial)              :: output
  
  integer, allocatable :: mode_ids(:)
  
  integer :: id
  integer :: paired_id
  integer :: power
  
  integer :: i,ialloc
  
  mode_ids = this%ids
  mode_ids = mode_ids(set(mode_ids))
  mode_ids = mode_ids(sort(mode_ids))
  output%coefficient = 1
  allocate(output%modes(size(mode_ids)), stat=ialloc); call err(ialloc)
  do i=1,size(mode_ids)
    id = mode_ids(i)
    paired_id = real_modes(first(real_modes%id==id))%paired_id
    power = count(this%ids==id)
    output%modes(i) = RealUnivariate( id        = id,        &
                                    & paired_id = paired_id, &
                                    & power     = power)
  enddo
end function

! e.g. if this has ids [2,1,2,3,1], the corresponding ModeCoupling
!    has ids [1,2,3].
impure elemental function mode_coupling(this) result(output)
  implicit none
  
  class(ModeMonomial), intent(in) :: this
  type(ModeCoupling)              :: output
  
  integer, allocatable :: ids(:)
  
  ids = this%ids
  ids = ids(set(ids))
  ids = ids(sort(ids))
  output = ModeCoupling(ids)
end function
end module
