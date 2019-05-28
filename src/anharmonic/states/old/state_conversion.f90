! ======================================================================
! Conversions between bases. Will eventually become a sparse matrix class.
! ======================================================================
module state_conversion_module
  use common_module
  
  use monomial_state_module
  implicit none
  
  private
  
  public :: StateConversion
  public :: size
  
  ! Conversions between bases of states.
  type, extends(Stringable) :: StateConversion
    integer,  allocatable, private :: ids_(:)
    real(dp), allocatable, private :: coefficients_(:)
  contains
    ! Getters for ids and coefficients.
    procedure, public :: id => id_StateConversion
    procedure, public :: ids => ids_StateConversion
    procedure, public :: coefficient => coefficient_StateConversion
    procedure, public :: coefficients => coefficients_StateConversion
    
    ! Add an id and coefficient to the StateConversion.
    procedure, public :: add_element
    
    ! I/O.
    procedure, public :: read  => read_StateConversion
    procedure, public :: write => write_StateConversion
  end type
  
  interface StateConversion
    module procedure new_StateConversion_null
    module procedure new_StateConversion
    module procedure new_StateConversion_String
  end interface
  
  interface size
    module procedure size_StateConversion
  end interface
contains

! Constructors and size function.
function new_StateConversion_null() result(this)
  implicit none
  
  type(StateConversion) :: this
  
  this%ids_ = [integer::]
  this%coefficients_ = [real(dp)::]
end function

function new_StateConversion(ids,coefficients) result(this)
  implicit none
  
  integer,  intent(in)  :: ids(:)
  real(dp), intent(in)  :: coefficients(:)
  type(StateConversion) :: this
  
  integer, allocatable :: sort_key(:)
  
  if (size(ids)/=size(coefficients)) then
    call print_line(CODE_ERROR//': IDs and coefficients do not match.')
    call err()
  endif
  
  sort_key = sort(ids)
  
  this%ids_          = ids(sort_key)
  this%coefficients_ = coefficients(sort_key)
  
  if (size(this)>1) then
    if (any(this%ids_(2:)==this%ids_(:size(this)-1))) then
      call print_line(CODE_ERROR//': An ID has been given twice.')
      call err()
    endif
  endif
end function

function size_StateConversion(this) result(output)
  implicit none
  
  type(StateConversion), intent(in) :: this
  integer                           :: output
  
  output = size(this%ids_)
end function

! Getters.
impure elemental function id_StateConversion(this,index) result(output)
  implicit none
  
  class(StateConversion), intent(in) :: this
  integer,                intent(in) :: index
  integer                            :: output
  
  output = this%ids_(index)
end function

function ids_StateConversion(this) result(output)
  implicit none
  
  class(StateConversion), intent(in) :: this
  integer, allocatable               :: output(:)
  
  output = this%ids_
end function

impure elemental function coefficient_StateConversion(this,index) &
   & result(output)
  implicit none
  
  class(StateConversion), intent(in) :: this
  integer,                intent(in) :: index
  real(dp)                           :: output
  
  output = this%coefficients_(index)
end function

function coefficients_StateConversion(this) result(output)
  implicit none
  
  class(StateConversion), intent(in) :: this
  real(dp), allocatable              :: output(:)
  
  output = this%coefficients_
end function

! ----------------------------------------------------------------------
! Add and ID and coefficient to a StateConversion.
! ----------------------------------------------------------------------
! Assumes that this%ids_ is sorted, and keeps it sorted.
subroutine add_element(this,id,coefficient)
  implicit none
  
  class(StateConversion), intent(inout) :: this
  integer,                intent(in)    :: id
  real(dp),               intent(in)    :: coefficient
  
  integer :: i
  
  if (.not. allocated(this%ids_)) then
    call print_line(CODE_ERROR//': Trying to add an element to a &
       &StateConversion which has not been allocated.')
    call err()
  endif
  
  i = first(this%ids_>=id, default=0, sorted=.true.)
  if (i==0) then
    this%ids_ = [this%ids_, id]
    this%coefficients_ = [this%coefficients_, coefficient]
  else
    if (this%ids_(i)==id) then
      this%coefficients_(i) = this%coefficients_(i) + coefficient
    else
      this%ids_ = [this%ids_(:i-1), id, this%ids_(i:)]
      this%coefficients_ = [ this%coefficients_(:i-1), &
                           & coefficient,              &
                           & this%coefficients_(i:)    ]
    endif
  endif
end subroutine

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_StateConversion(this,input)
  implicit none
  
  class(StateConversion), intent(out) :: this
  type(String),           intent(in)  :: input
  
  integer,  allocatable :: ids(:)
  real(dp), allocatable :: coefficients(:)
  
  type(String), allocatable :: states(:)
  type(String), allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(StateConversion)
    states = split_line(input)
    states = states(filter([(modulo(i,2)==1,i=1,size(states))]))
    allocate( ids(size(states)),          &
            & coefficients(size(states)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(states)
      line = split_line(states(i), delimiter='|')
      coefficients(i) = dble(line(1))
      ids(i) = int(slice(line(2),1,len(line(2))-1))
    enddo
    
    this = StateConversion(ids, coefficients)
  class default
    call err()
  end select
end subroutine

function write_StateConversion(this) result(output)
  implicit none
  
  class(StateConversion), intent(in) :: this
  type(String)                       :: output
  
  integer :: i
  
  select type(this); type is(StateConversion)
    output = join( [( this%coefficients_(i)//'|'//this%ids_(i)//'>',     &
                 &    i=1,                                               &
                 &    size(this)                                     )], &
                 & delimiter=' + '                                       )
  class default
    call err()
  end select
end function

function new_StateConversion_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(StateConversion)    :: this
  
  call this%read(input)
end function
end module
