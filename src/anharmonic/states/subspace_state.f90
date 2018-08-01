! ======================================================================
! A product of states along each complex mode in a degenerate subspace.
! ======================================================================
module subspace_state_module
  use common_module
  
  use degenerate_subspace_module
  use mode_state_module
  implicit none
  
  private
  
  public :: SubspaceState
  public :: size
  public :: braket
  
  type, extends(Stringsable) :: SubspaceState
    integer               :: subspace_id
    real(dp)              :: frequency
    type(ComplexMonomial) :: state
  contains
    procedure, public :: read  => read_SubspaceState
    procedure, public :: write => write_SubspaceState
  end type
  
  interface SubspaceState
    module procedure new_SubspaceState
    module procedure new_SubspaceState_Strings
    module procedure new_SubspaceState_StringArray
  end interface
  
  interface size
    module procedure size_SubspaceState
  end interface
  
  interface braket
    module procedure braket_SubspaceStates
    module procedure braket_SubspaceStates_ComplexMonomial
    module procedure braket_SubspaceStates_ComplexPolynomial
  end interface
contains

! ----------------------------------------------------------------------
! Basic functionality:
!    - constructor.
!    - size() function.
! ----------------------------------------------------------------------
function new_SubspaceState(subspace_id,frequency,state) result(this)
  implicit none
  
  integer,               intent(in) :: subspace_id
  real(dp),              intent(in) :: frequency
  type(ComplexMonomial), intent(in) :: state
  type(SubspaceState)               :: this
  
  this%subspace_id = subspace_id
  this%frequency   = frequency
  this%state       = state
end function

function size_SubspaceState(this) result(output)
  implicit none
  
  type(SubspaceState), intent(in) :: this
  integer                         :: output
  
  output = size(this%state)
end function

! ----------------------------------------------------------------------
! Evaluates integrals of the form <bra|ket> and <bra|monomial|ket>.
! ----------------------------------------------------------------------
function braket_SubspaceStates(bra,ket,subspace,supercell) result(output)
  implicit none
  
  type(SubspaceState),      intent(in) :: bra
  type(SubspaceState),      intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  complex(dp)                          :: output
  
  real(dp) :: denominator
  
  type(ComplexMonomial) :: powers
  
  logical, allocatable :: mode_integrated(:)
  
  integer :: n,m
  
  integer :: i,j,k,ialloc
  
  ! u=u*:
  ! <|u^n|> = 0                                if n odd
  !         = prod_{k=1}^{n/2}[ (2k-1)/(2Nw) ] if n even
  
  ! u/=u*:
  ! <|u^nu^*^m|> = 0 if n/=m
  !              = prod_{k=1}^n[ k/(2Nw) ]
  
  if (bra%subspace_id/=ket%subspace_id) then
    call print_line(ERROR//': bra and ket from different subspaces.')
    call err()
  elseif (bra%subspace_id/=subspace%id) then
    call print_line(ERROR//': bra and subspace do not match.')
    call err()
  endif
  
  denominator = 2*bra%frequency*supercell%sc_size
  
  powers = conjg(bra%state) * ket%state
  
  allocate(mode_integrated(size(powers)), stat=ialloc); call err(ialloc)
  mode_integrated = .false.
  output = powers%coefficient
  do i=1,size(powers)
    if (mode_integrated(i)) then
      cycle
    endif
    
    if (.not. any(powers%modes(i)%id==subspace%mode_ids)) then
      call print_line(ERROR//': <bra|ket> integration contains modes outside &
         &of the given degenerate subspace.')
    endif
    
    if (powers%modes(i)%id==powers%modes(i)%paired_id) then
      n = powers%modes(i)%power
      
      if (modulo(n,2)==1) then
        output = 0.0_dp
        return
      endif
      
      do k=1,n/2
        output = output * (2*k-1)/denominator
      enddo
      
      mode_integrated(i) = .true.
    else
      j = first(powers%modes%id==powers%modes(i)%paired_id)
      
      n = powers%modes(i)%power
      m = powers%modes(j)%power
      
      if (n/=m) then
        output = 0.0_dp
        return
      endif
      
      do k=1,n
        output = output * k/denominator
      enddo
      
      mode_integrated(i) = .true.
      mode_integrated(j) = .true.
    endif
  enddo
  
  if (.not. all(mode_integrated)) then
    call print_line(CODE_ERROR//': Failed to integrate <bra|ket>.')
    call err()
  endif
end function

function braket_SubspaceStates_ComplexMonomial(bra,ket,monomial,subspace, &
   & supercell) result(output)
  implicit none
  
  type(SubspaceState),      intent(in) :: bra
  type(SubspaceState),      intent(in) :: ket
  type(ComplexMonomial),    intent(in) :: monomial
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  type(ComplexMonomial)                :: output
  
  logical,                 allocatable :: mode_in_subspace(:)
  type(ComplexUnivariate), allocatable :: subspace_modes(:)
  type(ComplexUnivariate), allocatable :: non_subspace_modes(:)
  
  type(ComplexMonomial) :: subspace_monomial
  type(SubspaceState)   :: monomial_times_ket
  
  complex(dp) :: coefficient
  
  integer :: i,ialloc
  
  allocate(mode_in_subspace(size(monomial)), stat=ialloc); call err(ialloc)
  mode_in_subspace = .false.
  do i=1,size(monomial)
    if (any(monomial%modes(i)%id==subspace%mode_ids)) then
      mode_in_subspace(i) = .true.
    endif
  enddo
  
  subspace_modes = monomial%modes(filter(mode_in_subspace))
  non_subspace_modes = monomial%modes(filter(.not.mode_in_subspace))
  
  subspace_monomial = ComplexMonomial( coefficient = monomial%coefficient, &
                                     & modes       = subspace_modes )
  
  monomial_times_ket = ket
  monomial_times_ket%state = monomial_times_ket%state * subspace_monomial
  
  coefficient = braket(bra,monomial_times_ket,subspace,supercell)
  
  output = ComplexMonomial( coefficient = coefficient, &
                          & modes = non_subspace_modes )
end function

function braket_SubspaceStates_ComplexPolynomial(bra,ket,polynomial,subspace, &
   & supercell) result(output)
  implicit none
  
  type(SubspaceState),      intent(in) :: bra
  type(SubspaceState),      intent(in) :: ket
  type(ComplexPolynomial),  intent(in) :: polynomial
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  type(ComplexPolynomial)              :: output
  
  integer :: i
  
  output = polynomial
  
  do i=1,size(output)
    output%terms(i) = braket(bra,ket,output%terms(i),subspace,supercell)
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_SubspaceState(this,input)
  implicit none
  
  class(SubspaceState), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  integer               :: subspace_id
  real(dp)              :: frequency
  type(ComplexMonomial) :: state
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(SubspaceState)
    line = split_line(input(1))
    subspace_id = int(line(2))
    
    line = split_line(input(2))
    frequency = dble(line(2))
    
    state = ComplexMonomial(input(4))
    
    this = SubspaceState(subspace_id,frequency,state)
  end select
end subroutine

function write_SubspaceState(this) result(output)
  implicit none
  
  class(SubspaceState), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  select type(this); type is(SubspaceState)
    output = [ 'Subspace '//this%subspace_id, &
             & 'Frequency '//this%frequency,  &
             & str('State'),                  &
             & str(this%state)                ]
  end select
end function

function new_SubspaceState_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(SubspaceState)      :: this
  
  call this%read(input)
end function

impure elemental function new_SubspaceState_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SubspaceState)           :: this
  
  this = SubspaceState(str(input))
end function
end module
