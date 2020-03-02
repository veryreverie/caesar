! ======================================================================
! A "subspace" is a set of degenerate modes.
! A subspace monomial is a product of subspaces.
! ======================================================================
module subspace_monomial_module
  use common_module
  
  use subspace_coupling_module
  implicit none
  
  private
  
  public :: SubspaceMonomial
  public :: size
  public :: generate_subspace_monomials
  public :: operator(//)
  public :: operator(==)
  public :: operator(/=)
  public :: generate_complex_monomials
  public :: generate_real_monomials
  
  type, extends(Stringable) :: SubspaceMonomial
    ! The ids and powers of the degenerate subspaces in the monomial.
    ! e.g. if ids=[1,2,5] and powers=[1,2,1] then the monomial is
    !    (s1).(s2)^2.(s5), where s1 is the subspace with id 1.
    integer, allocatable :: ids(:)
    integer, allocatable :: powers(:)
  contains
    ! Returns the coupled subspaces.
    procedure, public :: coupled_subspaces
    
    ! Check if this is subsidiary to the given coupling, e.g. the couplings
    !    [],[1],[2],[3],[1,2],[1,3] and [2,3] are subsidiaries of [1,2,3].
    procedure, public :: is_subsidiary_of
    
    ! I/O.
    procedure, public :: read  => read_SubspaceMonomial
    procedure, public :: write => write_SubspaceMonomial
  end type
  
  interface SubspaceMonomial
    module procedure new_SubspaceMonomial
    module procedure new_SubspaceMonomial_ids_powers
    module procedure new_SubspaceMonomial_DegenerateSubspaces
    module procedure new_SubspaceMonomial_String
  end interface
  
  interface operator(//)
    module procedure concatenate_SubspaceMonomial_DegenerateSubspace
  end interface
  
  interface size
    module procedure size_SubspaceMonomial
  end interface
  
  interface operator(==)
    module procedure equality_SubspaceMonomial_SubspaceMonomial
  end interface
  
  interface operator(/=)
    module procedure non_equality_SubspaceMonomial_SubspaceMonomial
  end interface
contains

! ----------------------------------------------------------------------
! Basic functionality.
!    - Constructors.
!    - Concatenation with DegenerateSubspace.
!    - size() function.
!    - equality and non-equality with other SubspaceMonomials.
! ----------------------------------------------------------------------
function new_SubspaceMonomial() result(this)
  implicit none
  
  type(SubspaceMonomial) :: this
  
  integer :: ialloc
  
  allocate( this%ids(0),    &
          & this%powers(0), &
          & stat=ialloc); call err(ialloc)
end function

function new_SubspaceMonomial_ids_powers(ids,powers) result(this)
  implicit none
  
  integer, intent(in)    :: ids(:)
  integer, intent(in)    :: powers(:)
  type(SubspaceMonomial) :: this
  
  if (size(ids)/=size(powers)) then
    call print_line(CODE_ERROR//': IDs and powers do not match.')
    call err()
  endif
end function

function new_SubspaceMonomial_DegenerateSubspaces(subspaces) result(this)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  type(SubspaceMonomial)               :: this
  
  integer :: i
  
  this%ids = subspaces%id
  this%ids = this%ids(sort(this%ids))
  this%ids = this%ids(set(this%ids))
  this%powers = [(count(subspaces%id==this%ids(i)), i=1, size(subspaces))]
end function

function concatenate_SubspaceMonomial_DegenerateSubspace(this,subspace) &
   & result(output)
  implicit none
  
  type(SubspaceMonomial),   intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  type(SubspaceMonomial)               :: output
  
  integer :: i
  
  output = this
  
  i = first(this%ids==subspace%id, default=0)
  if (i==0) then
    output%ids = [output%ids, subspace%id]
    output%powers = [output%powers, 1]
  else
    output%powers(i) = output%powers(i) + 1
  endif
end function

function size_SubspaceMonomial(this) result(output)
  implicit none
  
  type(SubspaceMonomial), intent(in) :: this
  integer                            :: output
  
  output = size(this%ids)
end function

impure elemental function equality_SubspaceMonomial_SubspaceMonomial(this, &
   & that) result(output)
  implicit none
  
  type(SubspaceMonomial), intent(in) :: this
  type(SubspaceMonomial), intent(in) :: that
  logical                            :: output
  
  if (size(this)/=size(that)) then
    output = .false.
  else
    output = all(this%ids==that%ids) .and. all(this%powers==that%powers)
  endif
end function

impure elemental function non_equality_SubspaceMonomial_SubspaceMonomial( &
   & this,that) result(output)
  implicit none
  
  type(SubspaceMonomial), intent(in) :: this
  type(SubspaceMonomial), intent(in) :: that
  logical                            :: output
  
  output = .not. this==that
end function

! ----------------------------------------------------------------------
! Use stored IDs to return the actual objects the IDs represent.
! ----------------------------------------------------------------------
function coupled_subspaces(this,subspaces) result(output)
  implicit none
  
  class(SubspaceMonomial),  intent(in)  :: this
  type(DegenerateSubspace), intent(in)  :: subspaces(:)
  type(DegenerateSubspace), allocatable :: output(:)
  
  integer :: i
  
  output = [( subspaces(first(subspaces%id==this%ids(i))), i=1, size(this) )]
end function

! ----------------------------------------------------------------------
! Check if a coupling is a subsidiary of another coupling.
! ----------------------------------------------------------------------
function is_subsidiary_of(this,that) result(output)
  implicit none
  
  Class(SubspaceMonomial), intent(in) :: this
  type(SubspaceMonomial),  intent(in) :: that
  logical                             :: output
  
  integer :: i,j
  
  ! Loop over the ids in this, checking if each is in that.
  do i=1,size(this)
    j = first(that%ids==this%ids(i), default=0)
    if (j==0) then
      output = .false.
      return
    elseif (this%powers(i)>that%powers(j)) then
      output = .false.
      return
    endif
  enddo
  
  output = .true.
end function

! ----------------------------------------------------------------------
! Generates all coupling monomials corresponding to a given subspace coupling,
!    up to the given potential expansion order.
! ----------------------------------------------------------------------
! e.g. given the subspace coupling with ids [3,5] and expansion order 3,
!    will return coupling monomials [3,3,5] and [3,5,5].
function generate_subspace_monomials(subspace_coupling,subspaces, &
   & minimum_expansion_order,maximum_expansion_order) result(output)
  implicit none
  
  type(SubspaceCoupling),   intent(in) :: subspace_coupling
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  integer,                  intent(in) :: minimum_expansion_order
  integer,                  intent(in) :: maximum_expansion_order
  type(SubspaceMonomial), allocatable  :: output(:)
  
  type(DegenerateSubspace), allocatable :: coupled_subspaces(:)
  
  ! Check inputs.
  if (size(subspace_coupling)==0) then
    call print_line(CODE_ERROR//': Empty subspace coupling.')
    call err()
  elseif (minimum_expansion_order<0) then
    call print_line(ERROR//': minimum_expansion_order must be non-negative.')
    call quit()
  elseif (   maximum_expansion_order                               &
         & < min(minimum_expansion_order, size(subspace_coupling)) ) then
    call print_line(ERROR//': maximum_expansion_order must be at least as &
       &large as minimum_expansion_order and the size of the coupling.')
    call quit()
  endif
  
  ! Retrieve coupled subspaces from subspace coupling.
  coupled_subspaces = subspace_coupling%coupled_subspaces(subspaces)
  
  ! Call recursive helper function to generate monomials.
  output = generate_subspace_monomials_helper( coupled_subspaces,       &
                                             & minimum_expansion_order, &
                                             & maximum_expansion_order  )
end function

! Helper function for generate_subspace_monomials.
! This function appends a number of copies of the first subspace in
!    coupled_subspaces, then calls itself to append the next subspace etc.
! The optional argument monomial_in is a recursive argument which should only
!    be provided by recursive calls.
recursive function generate_subspace_monomials_helper(coupled_subspaces, &
   & minimum_expansion_order,maximum_expansion_order,monomial_in)        &
   & result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in)           :: coupled_subspaces(:)
  integer,                  intent(in)           :: minimum_expansion_order
  integer,                  intent(in)           :: maximum_expansion_order
  type(SubspaceMonomial),   intent(in), optional :: monomial_in
  type(SubspaceMonomial), allocatable            :: output(:)
  
  type(DegenerateSubspace)              :: first_subspace
  type(DegenerateSubspace), allocatable :: remaining_subspaces(:)
  
  type(SubspaceMonomial) :: monomial
  
  integer :: ialloc
  
  ! If there are no more subspaces to append, return the monomial,
  !    but only if its size is at least minimum_expansion_order.
  ! Size 1 monomials are ignored because they correspond to linear terms in the
  !    potential, which are zero because the structure is geometry optimised.
  if (size(coupled_subspaces)==0) then
    if (.not. present(monomial_in)) then
      call print_line(CODE_ERROR//': Empty subspace coupling.')
      call err()
    elseif (sum(monomial_in%powers)<minimum_expansion_order) then
      allocate(output(0), stat=ialloc); call err(ialloc)
    else
      output = [monomial_in]
    endif
    
    return
  endif
  
  ! Split the coupled subspaces into the first subspace, which will be handled
  !    by this call of the function, and all the rest, which will be handled
  !    by a recursive call.
  first_subspace      = coupled_subspaces(1)
  remaining_subspaces = coupled_subspaces(2:)
  
  ! Append copies of the first subspace as many times as still leaves space
  !    for at least one copy of every remaining subspace, and then call this
  !    function to handle the next subspace along.
  if (present(monomial_in)) then
    monomial = monomial_in
  else
    monomial = SubspaceMonomial()
  endif
  
  allocate(output(0), stat=ialloc); call err(ialloc)
  do while(    sum(monomial%powers)+size(coupled_subspaces) &
          & <= maximum_expansion_order                      )
    monomial = monomial // first_subspace
    output = [ output,                                                      &
           &   generate_subspace_monomials_helper( remaining_subspaces,     &
           &                                       minimum_expansion_order, &
           &                                       maximum_expansion_order, &
           &                                       monomial)                &
           & ]
  enddo
end function

! ----------------------------------------------------------------------
! Generate the ComplexMonomials or RealMonomials corresponding to
!    a given SubspaceMonomial.
! ----------------------------------------------------------------------
! The coefficients are chosen such that symmetry operations are unitary
!    in the basis of monomials.
function generate_complex_monomials(this,subspaces,modes,qpoints, &
   & conserve_momentum,conserve_subspace_momentum) result(output)
  implicit none
  
  type(SubspaceMonomial),   intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  type(ComplexMode),        intent(in) :: modes(:)
  type(QpointData),         intent(in) :: qpoints(:)
  logical,                  intent(in) :: conserve_momentum
  logical,                  intent(in) :: conserve_subspace_momentum
  type(ComplexMonomial), allocatable   :: output(:)
  
  type(ComplexUnivariate) :: zero_univariate(0)
  type(ComplexMonomial)   :: root
  
  type(DegenerateSubspace)       :: subspace
  type(ComplexMode), allocatable :: subspace_modes(:)
  
  type(ComplexMonomial), allocatable :: subspace_monomials(:)
  
  integer :: i,j,k,ialloc
  
  if (size(this)==0) then
    allocate(output(0), stat=ialloc); call err(ialloc)
    return
  endif
  
  root = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                        & modes       = zero_univariate          )
  
  do i=1,size(this)
    subspace = subspaces(first(subspaces%id==this%ids(i)))
    subspace_modes = subspace%modes(modes)
    subspace_monomials = generate_subspace_complex_monomials( subspace_modes, &
                                                            & this%powers(i), &
                                                            & root            )
    if (conserve_subspace_momentum) then
      subspace_monomials = subspace_monomials(filter( subspace_monomials, &
                                                    & conserves_momentum  ))
    endif
    
    if (i==1) then
      output = subspace_monomials
    else
      output = [(                                                            &
         & (output(k)*subspace_monomials(j), j=1, size(subspace_monomials)), &
         & k=1,                                                              &
         & size(output)                                                      )]
    endif
  enddo
  
  if (conserve_momentum) then
    output = output(filter(output, conserves_momentum))
  endif
contains
  ! Lambda for checking if a monomial conserves momentum.
  ! Captures:
  !    - modes
  !    - qpoints
  ! N.B. a monomial given by
  !    prod_i (u_{q_i,i})^{n_i}
  ! conserves momentum iff sum_i n_i q_i is a G-vector.
  function conserves_momentum(input) result(output)
    implicit none
    
    class(*), intent(in) :: input
    logical              :: output
    
    select type(input); type is(ComplexMonomial)
      output = input%wavevector(modes,qpoints)==fracvec(zeroes(3))
    end select
  end function
end function

function generate_real_monomials(this,subspaces,modes,qpoints) result(output)
  implicit none
  
  type(SubspaceMonomial),   intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  type(RealMode),           intent(in) :: modes(:)
  type(QpointData),         intent(in) :: qpoints(:)
  type(RealMonomial), allocatable      :: output(:)
  
  type(RealUnivariate) :: zero_univariate(0)
  type(RealMonomial)   :: root
  
  type(DegenerateSubspace)    :: subspace
  type(RealMode), allocatable :: subspace_modes(:)
  
  type(RealMonomial), allocatable :: subspace_monomials(:)
  
  integer :: i,j,k,ialloc
  
  if (size(this)==0) then
    allocate(output(0), stat=ialloc); call err(ialloc)
    return
  endif
  
  root = RealMonomial(coefficient=1.0_dp, modes=zero_univariate)
  
  do i=1,size(this)
    subspace = subspaces(first(subspaces%id==this%ids(i)))
    subspace_modes = subspace%modes(modes)
    subspace_monomials = generate_subspace_real_monomials( subspace_modes, &
                                                         & this%powers(i), &
                                                         & root            )
    if (i==1) then
      output = subspace_monomials
    else
      output = [(                                                            &
         & (output(k)*subspace_monomials(j), j=1, size(subspace_monomials)), &
         & k=1,                                                              &
         & size(output)                                                      )]
    endif
  enddo
end function

! Helper functions for the above. These generate monomials for a single
!    subspace.
recursive function generate_subspace_complex_monomials(modes,power,root) &
   & result(output)
  implicit none
  
  type(ComplexMode),     intent(in)  :: modes(:)
  integer,               intent(in)  :: power
  type(ComplexMonomial), intent(in)  :: root
  type(ComplexMonomial), allocatable :: output(:)
  
  integer :: i
  
  if (size(modes)==0) then
    if (power>0) then
      call print_line(CODE_ERROR//': Monomial requires further modes, but no &
         &modes remain to append.')
      call err()
    endif
  endif
  
  if (power==0) then
    output = [root]
    output%coefficient = sqrt(no_permutations(output(1)))
  elseif (size(modes)==1) then
    output = [root*ComplexUnivariate(mode=modes(1), power=power)]
    output%coefficient = sqrt(no_permutations(output(1)))
  else
    output = [generate_subspace_complex_monomials(modes(2:),power,root)]
    do i=1,power
      output = [ output,                                            &
               & generate_subspace_complex_monomials(               &
               &   modes(2:),                                       &
               &   power-i,                                         &
               &   root*ComplexUnivariate(mode=modes(1), power=i) ) ]
    enddo
  endif
end function

recursive function generate_subspace_real_monomials(modes,power,root) &
   & result(output)
  implicit none
  
  type(RealMode),     intent(in)  :: modes(:)
  integer,            intent(in)  :: power
  type(RealMonomial), intent(in)  :: root
  type(RealMonomial), allocatable :: output(:)
  
  integer :: i
  
  if (size(modes)==0) then
    if (power>0) then
      call print_line(CODE_ERROR//': Monomial requires further modes, but no &
         &modes remain to append.')
      call err()
    endif
  endif
  
  if (power==0) then
    output = [root]
    output%coefficient = sqrt(no_permutations(output(1)))
  elseif (size(modes)==1) then
    output = [root*RealUnivariate(mode=modes(1), power=power)]
    output%coefficient = sqrt(no_permutations(output(1)))
  else
    output = [generate_subspace_real_monomials(modes(2:),power,root)]
    do i=1,power
      output = [ output,                                         &
               & generate_subspace_real_monomials(               &
               &   modes(2:),                                    &
               &   power-i,                                      &
               &   root*RealUnivariate(mode=modes(1), power=i) ) ]
    enddo
  endif
end function

function no_permutations(input) result(output)
  implicit none
  
  class(*), intent(in) :: input
  real(dp)             :: output
  
  integer, allocatable :: powers(:)
  
  integer :: i,ialloc
  
  allocate(powers(0), stat=ialloc); call err(ialloc)
  select type(input); type is (ComplexMonomial)
    do i=1,size(input)
      powers = [powers, input%power(i)]
      if (input%id(i)/=input%paired_id(i)) then
        powers = [powers, input%paired_power(i)]
      endif
    enddo
  type is(RealMonomial)
    do i=1,size(input)
      powers = [powers, input%power(i)]
      if (input%id(i)/=input%paired_id(i)) then
        powers = [powers, input%paired_power(i)]
      endif
    enddo
  class default
    call err()
  end select
  
  output = real_multinomial(sum(powers), powers)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_SubspaceMonomial(this,input)
  implicit none
  
  class(SubspaceMonomial), intent(out) :: this
  type(String),            intent(in)  :: input
  
  integer, allocatable :: ids(:)
  integer, allocatable :: powers(:)
  
  type(String), allocatable :: subspace_strings(:)
  type(String), allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(SubspaceMonomial)
    subspace_strings = split_line(input, delimiter='*')
    allocate( ids(size(subspace_strings)),    &
            & powers(size(subspace_strings)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(subspace_strings)
      line = split_line(subspace_strings(i), delimiter='^')
      ids(i) = int(slice(line(1),2,len(line(1))))
      powers(i) = int(slice(line(2),1,len(line(2))-1))
    enddo
    this = SubspaceMonomial(ids, powers)
  end select
end subroutine

function write_SubspaceMonomial(this) result(output)
  implicit none
  
  class(SubspaceMonomial), intent(in) :: this
  type(String)                        :: output
  
  type(String), allocatable :: subspace_strings(:)
  
  integer :: i
  
  select type(this); type is(SubspaceMonomial)
    subspace_strings = [( '(s'//this%ids(i)//'^'//this%powers(i)//')', &
                        & i=1,                                         &
                        & size(this)                                   )]
    output = join(subspace_strings, delimiter='*')
  end select
end function

impure elemental function new_SubspaceMonomial_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(SubspaceMonomial)   :: this
  
  call this%read(input)
end function
end module
