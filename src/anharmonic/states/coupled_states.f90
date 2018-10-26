! ======================================================================
! For a given state |i>, records which states |j> have non-zero <i|H|j>.
! ======================================================================
module coupled_states_module
  use common_module
  implicit none
  
  private
  
  public :: CoupledStates
  public :: size
  
  type, extends(Stringable) :: CoupledStates
    integer, allocatable, private :: ids_(:)
    integer, allocatable, private :: separations_(:)
  contains
    ! Getters for ids and separations.
    procedure, public :: id => id_CoupledStates
    procedure, public :: ids => ids_CoupledStates
    procedure, public :: separation => separation_CoupledStates
    procedure, public :: separations => separations_CoupledStates
    
    ! Add and id and separation to the CoupledStates.
    procedure, public :: add_coupling
    
    ! I/O.
    procedure, public :: read  => read_CoupledStates
    procedure, public :: write => write_CoupledStates
  end type
  
  interface CoupledStates
    module procedure new_CoupledStates_null
    module procedure new_CoupledStates
    module procedure new_CoupledStates_String
  end interface
  
  interface size
    module procedure size_CoupledStates
  end interface
contains

! Constructors and size function.
function new_CoupledStates_null() result(this)
  implicit none
  
  type(CoupledStates) :: this
  
  this%ids_ = [integer::]
  this%separations_ = [integer::]
end function

function new_CoupledStates(ids,separations) result(this)
  implicit none
  
  integer, intent(in) :: ids(:)
  integer, intent(in) :: separations(:)
  type(CoupledStates) :: this
  
  integer, allocatable :: sort_key(:)
  
  if (size(ids)/=size(separations)) then
    call print_line(CODE_ERROR//': ids and separations do not match.')
    call err()
  endif
  
  sort_key = sort(ids)
  
  this%ids_ = ids(sort_key)
  this%separations_ = separations(sort_key)
  
  if (size(this)>1) then
    if (any(this%ids_(2:)==this%ids_(:size(this)-1))) then
      call print_line(CODE_ERROR//': An ID has been given twice.')
      call err()
    endif
  endif
end function

function size_CoupledStates(this) result(output)
  implicit none
  
  type(CoupledStates), intent(in) :: this
  integer                         :: output
  
  output = size(this%ids_)
end function

! Getters.
impure elemental function id_CoupledStates(this,index) result(output)
  implicit none
  
  class(CoupledStates), intent(in) :: this
  integer,              intent(in) :: index
  integer                          :: output
  
  output = this%ids_(index)
end function

function ids_CoupledStates(this) result(output)
  implicit none
  
  class(CoupledStates), intent(in) :: this
  integer, allocatable             :: output(:)
  
  output = this%ids_
end function

impure elemental function separation_CoupledStates(this,index) &
   & result(output)
  implicit none
  
  class(CoupledStates), intent(in) :: this
  integer,              intent(in) :: index
  real(dp)                         :: output
  
  output = this%separations_(index)
end function

function separations_CoupledStates(this) result(output)
  implicit none
  
  class(CoupledStates), intent(in) :: this
  real(dp), allocatable            :: output(:)
  
  output = this%separations_
end function

! ----------------------------------------------------------------------
! Add a coupling to a CoupledStates.
! ----------------------------------------------------------------------
! Assumes that this%ids_ is sorted, and keeps it sorted.
subroutine add_coupling(this,id,separation)
  implicit none
  
  class(CoupledStates), intent(inout) :: this
  integer,              intent(in)    :: id
  integer,              intent(in)    :: separation
  
  integer :: i
  
  if (.not. allocated(this%ids_)) then
    call print_line(CODE_ERROR//': Trying to add a coupling to a &
       &CoupledStates which has not been allocated.')
    call err()
  endif
  
  i = first(this%ids_>=id, default=0, sorted=.true.)
  if (i==0) then
    this%ids_ = [this%ids_, id]
    this%separations_ = [this%separations_, separation]
  else
    if (this%ids_(i)==id) then
      call print_line(CODE_ERROR//': Trying to add a coupling which has &
         &already been added.')
      call err()
    else
      this%ids_ = [this%ids_(:i-1), id, this%ids_(i:)]
      this%separations_ = [ this%separations_(:i-1), &
                          & separation,              &
                          & this%separations_(i:)    ]
    endif
  endif
end subroutine

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_CoupledStates(this,input)
  implicit none
  
  class(CoupledStates), intent(out) :: this
  type(String),         intent(in)  :: input
  
  integer, allocatable :: ids(:)
  integer, allocatable :: separations(:)
  
  type(String), allocatable :: line(:)
  type(String), allocatable :: token(:)
  
  integer :: i,ialloc
  
  select type(this); type is(CoupledStates)
    line = split_line(input)
    allocate( ids(size(line)),         &
            & separations(size(line)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(line)
      token = split_line(line(i),delimiter=':')
      ids(i) = int(token(1))
      separations(i) = int(token(2))
    enddo
    this = CoupledStates(ids,separations)
  class default
    call err()
  end select
end subroutine

function write_CoupledStates(this) result(output)
  implicit none
  
  class(CoupledStates), intent(in) :: this
  type(String)                     :: output
  
  integer :: i
  
  select type(this); type is(CoupledStates)
    output = join([(this%ids_(i)//':'//this%separations_(i),i=1,size(this))])
  class default
    call err()
  end select
end function

impure elemental function new_CoupledStates_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(CoupledStates)      :: this
  
  call this%read(input)
end function
end module
