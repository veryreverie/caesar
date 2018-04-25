! ======================================================================
! Records a single set of modes which are coupled together.
! ======================================================================
module coupled_modes_module
  use common_module
  
  use degeneracy_module
  use subspace_monomial_module
  implicit none
  
  private
  
  public :: CoupledModes
  public :: size
  public :: generate_mode_coupling
  public :: operator(//)
  public :: operator(==)
  public :: operator(/=)
  public :: construct_complex_monomial
  public :: construct_real_monomial
  
  ! A list of ids of modes which are coupled.
  type, extends(Stringable) :: CoupledModes
    integer, allocatable :: ids(:)
    ! Whether or not the complex basis functions corresponding to this coupling
    !    are part of the Hamiltonian. i.e. they conserve Bloch momentum.
    logical :: conserves_momentum
    ! Whether or not the complex basis functions corresponding to this coupling
    !    are part of the VSCF Hamiltonian. i.e. they conserve Bloch momentum
    !    within every coupled subspace.
    logical :: conserves_vscf
  contains
    ! I/O.
    procedure, public :: to_String => to_String_CoupledModes
  end type
  
  interface size
    module procedure size_CoupledModes
  end interface
  
  interface operator(//)
    module procedure concatenate_CoupledModes_integer
  end interface
  
  interface operator(==)
    module procedure equality_CoupledModes_CoupledModes
  end interface
  
  interface operator(/=)
    module procedure non_equality_CoupledModes_CoupledModes
  end interface
contains

! Constructs a complex monomial from a mode coupling.
function construct_complex_monomial(input,complex_modes) result(output)
  implicit none
  
  type(CoupledModes), intent(in) :: input
  type(ComplexMode),  intent(in) :: complex_modes(:)
  type(ComplexMonomial)          :: output
  
  integer, allocatable :: mode_ids(:)
  
  integer :: id
  integer :: paired_id
  integer :: power
  
  integer :: i,ialloc
  
  mode_ids = input%ids
  mode_ids = mode_ids(set(mode_ids))
  mode_ids = mode_ids(sort(mode_ids))
  output%coefficient = 1
  allocate(output%modes(size(mode_ids)), stat=ialloc); call err(ialloc)
  do i=1,size(mode_ids)
    id = mode_ids(i)
    paired_id = complex_modes(first(complex_modes%id==id))%paired_id
    power = count(input%ids==id)
    output%modes(i) = ComplexUnivariate( id        = id,        &
                                       & paired_id = paired_id, &
                                       & power     = power)
  enddo
end function

! Constructs a real monomial from a mode coupling.
function construct_real_monomial(input,real_modes) result(output)
  implicit none
  
  type(CoupledModes), intent(in) :: input
  type(RealMode),     intent(in) :: real_modes(:)
  type(RealMonomial)             :: output
  
  integer, allocatable :: mode_ids(:)
  
  integer :: id
  integer :: paired_id
  integer :: power
  
  integer :: i,ialloc
  
  mode_ids = input%ids
  mode_ids = mode_ids(set(mode_ids))
  mode_ids = mode_ids(sort(mode_ids))
  output%coefficient = 1
  allocate(output%modes(size(mode_ids)), stat=ialloc); call err(ialloc)
  do i=1,size(mode_ids)
    id = mode_ids(i)
    paired_id = real_modes(first(real_modes%id==id))%paired_id
    power = count(input%ids==id)
    output%modes(i) = RealUnivariate( id        = id,        &
                                    & paired_id = paired_id, &
                                    & power     = power)
  enddo
end function

! size() functions.
function size_CoupledModes(this) result(output)
  implicit none
  
  type(CoupledModes), intent(in) :: this
  integer                        :: output
  
  output = size(this%ids)
end function

! Append an id to the coupling.
function concatenate_CoupledModes_integer(this,id) result(output)
  implicit none
  
  type(CoupledModes), intent(in) :: this
  integer,            intent(in) :: id
  type(CoupledModes)             :: output
  
  output = CoupledModes( ids                = [this%ids, id],          &
                       & conserves_momentum = this%conserves_momentum, &
                       & conserves_vscf     = this%conserves_vscf)
end function

! Compare couplings.
impure elemental function equality_CoupledModes_CoupledModes(this,that) &
   & result(output)
  implicit none
  
  type(CoupledModes), intent(in) :: this
  type(CoupledModes), intent(in) :: that
  logical                        :: output
  
  output = all(this%ids==that%ids)
end function

impure elemental function non_equality_CoupledModes_CoupledModes(this,that) &
   & result(output)
  implicit none
  
  type(CoupledModes), intent(in) :: this
  type(CoupledModes), intent(in) :: that
  logical                        :: output
  
  output = .not. this==that
end function

! ----------------------------------------------------------------------
! Recursively generates all sets of coupled modes within a set of coupled
!    subspaces.
! Only returns couplings with sum(q)=0, modulo G-vectors.
! ----------------------------------------------------------------------
function generate_mode_coupling(coupling,subspaces,normal_modes,qpoints) &
   & result(output)
  implicit none
  
  type(SubspaceMonomial), intent(in) :: coupling
  type(DegenerateModes),  intent(in) :: subspaces(:)
  type(ComplexMode),      intent(in) :: normal_modes(:)
  type(QpointData),       intent(in) :: qpoints(:)
  type(CoupledModes), allocatable    :: output(:)
  
  type(DegenerateModes), allocatable :: coupled_subspaces(:)
  
  ! List the subspaces which are coupled together in this coupling.
  coupled_subspaces = coupling%coupled_subspaces(subspaces)
  
  ! Generate all the sets of coupled modes within the coupled subspaces.
  output = generate_mode_coupling_helper( coupled_subspaces, &
                                        & normal_modes,      &
                                        & qpoints)
end function

recursive function generate_mode_coupling_helper(coupled_subspaces, &
   & normal_modes,qpoints,coupled_modes_in,sum_q_in) result(output)
  implicit none
  
  type(DegenerateModes), intent(in)           :: coupled_subspaces(:)
  type(ComplexMode),     intent(in)           :: normal_modes(:)
  type(QpointData),      intent(in)           :: qpoints(:)
  type(CoupledModes),    intent(in), optional :: coupled_modes_in
  type(FractionVector),  intent(in), optional :: sum_q_in
  type(CoupledModes), allocatable             :: output(:)
  
  type(QpointData), allocatable :: subspace_qpoints(:)
  
  type(CoupledModes)   :: coupled_modes
  type(FractionVector) :: sum_q
  
  type(CoupledModes)   :: coupled_modes_out
  type(FractionVector) :: sum_q_out
  
  logical :: last_mode_in_coupling
  
  integer :: i
  
  if (present(coupled_modes_in) .neqv. present(sum_q_in)) then
    call print_line(CODE_ERROR//': generate_coupled_modes must be called with &
       &all optional arguments or none.')
    call err()
  endif
  
  if (present(coupled_modes_in)) then
    coupled_modes = coupled_modes_in
    sum_q = sum_q_in
  else
    coupled_modes = CoupledModes( ids                = [integer::], &
                                & conserves_momentum = .true.,      &
                                & conserves_vscf     = .true.)
    sum_q = fracvec(zeroes(3))
  endif
  
  if (size(coupled_subspaces)==0) then
    ! There is nothing else to append. Check that the sum across q-points of
    !    the mode coupling is zero, and return the mode couplings.
    if (.not. is_int(sum_q_in)) then
      coupled_modes%conserves_momentum = .false.
      coupled_modes%conserves_vscf     = .false.
    endif
    output = [coupled_modes]
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
    subspace_qpoints = coupled_subspaces(1)%qpoints(qpoints)
    output = [CoupledModes::]
    do i=1,size(coupled_subspaces(1))
      coupled_modes_out = coupled_modes//coupled_subspaces(1)%mode_ids(i)
      sum_q_out = sum_q + subspace_qpoints(i)%qpoint
      if (last_mode_in_coupling .and. .not. is_int(sum_q_out)) then
        coupled_modes_out%conserves_vscf = .false.
      endif
      output = [ output,                                               &
             &   generate_mode_coupling_helper( coupled_subspaces(2:), &
             &                                  normal_modes,          &
             &                                  qpoints,               &
             &                                  coupled_modes_out,     &
             &                                  sum_q_out)             &
             & ]
    enddo
  endif
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
recursive function to_String_CoupledModes(this) result(output)
  implicit none
  
  class(CoupledModes), intent(in) :: this
  type(String)                    :: output
  
  output = join(this%ids)
end function
end module
