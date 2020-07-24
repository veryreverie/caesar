! ======================================================================
! For a given state |i>, records which states |j> have non-zero <i|H|j>.
! ======================================================================
! This class gets appended to millions of times, so it handles its storage
!    like a C++ vector for speed reasons.
module coupled_states_module
  use common_module
  implicit none
  
  private
  
  public :: CoupledStates
  public :: size
  public :: tensor_product_couplings
  public :: selected_states_couplings
  
  type, extends(Stringable) :: CoupledStates
    integer,              private :: length_
    integer, allocatable, private :: ids_(:)
    integer, allocatable, private :: separations_(:)
  contains
    ! Getters for ids and separations.
    procedure, public :: id => id_CoupledStates
    procedure, public :: ids => ids_CoupledStates
    procedure, public :: separation => separation_CoupledStates
    procedure, public :: separations => separations_CoupledStates
    
    ! Add an id and separation to the CoupledStates.
    procedure, public :: add_coupling
    
    ! Re-map the ids.
    procedure, public :: map_ids
    
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
  
  this%length_ = 0
  this%ids_ = [integer::]
  this%separations_ = [integer::]
end function

function new_CoupledStates(ids,separations) result(this)
  implicit none
  
  integer, intent(in) :: ids(:)
  integer, intent(in) :: separations(:)
  type(CoupledStates) :: this
  
  if (size(ids)/=size(separations)) then
    call print_line(CODE_ERROR//': ids and separations do not match.')
    call err()
  endif
  
  this%length_ = size(ids)
  this%ids_ = ids
  this%separations_ = separations
end function

function size_CoupledStates(this) result(output)
  implicit none
  
  type(CoupledStates), intent(in) :: this
  integer                         :: output
  
  output = this%length_
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
  
  output = this%ids_(:this%length_)
end function

impure elemental function separation_CoupledStates(this,index) &
   & result(output)
  implicit none
  
  class(CoupledStates), intent(in) :: this
  integer,              intent(in) :: index
  integer                          :: output
  
  output = this%separations_(index)
end function

function separations_CoupledStates(this) result(output)
  implicit none
  
  class(CoupledStates), intent(in) :: this
  integer, allocatable             :: output(:)
  
  output = this%separations_(:this%length_)
end function

! Add a coupling to a CoupledStates.
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
  
  this%length_ = this%length_ + 1
  if (size(this%ids_)>=this%length_) then
    ! There is space left in the arrays; simply add the coupling.
    this%ids_(this%length_) = id
    this%separations_(this%length_) = separation
  else
    ! There is not sufficient space left in the array.
    ! First double the allocated space, then add the coupling.
    this%ids_ = [this%ids_, id, [(0,i=1,this%length_)]]
    this%separations_ = [this%separations_, separation, [(0,i=1,this%length_)]]
  endif
end subroutine

! Re-map the ids.
subroutine map_ids(this,id_map)
  implicit none
  
  class(CoupledStates), intent(inout) :: this
  integer,              intent(in)    :: id_map(:)
  
  this%ids_(:this%length_) = id_map(this%ids_(:this%length_))
end subroutine

! Take the tensor product of two couplings.
function tensor_product_couplings(unmapped_this,this,no_states,that, &
   & max_separation) result(output) 
  implicit none
  
  type(CoupledStates), intent(in) :: unmapped_this
  type(CoupledStates), intent(in) :: this
  integer,             intent(in) :: no_states(:)
  type(CoupledStates), intent(in) :: that
  integer,             intent(in) :: max_separation
  type(CoupledStates)             :: output
  
  integer :: length
  
  integer :: id
  integer :: separation
  
  integer :: i,j,k,ialloc
  
  length = 0
  do i=1,this%length_
    length = length                                                       &
         & + count( this%separations_(i)+that%separations_(:that%length_) &
         &       <= max_separation                                        &
         &    .and. that%ids_(:that%length_)                              &
         &       <= no_states(unmapped_this%ids_(i))                      )
  enddo
  
  allocate( output%ids_(length),         &
          & output%separations_(length), &
          & stat=ialloc); call err(ialloc)
  output%length_ = length
  
  k = 0
  do i=1,this%length_
    do j=1,that%length_
      if (that%ids_(j)>no_states(unmapped_this%ids_(i))) then
        exit
      endif
      separation = this%separations_(i) + that%separations_(j)
      if (separation>max_separation) then
        cycle
      endif
      id = this%ids_(i) + that%ids_(j)
      k = k+1
      output%ids_(k) = id
      output%separations_(k) = separation
    enddo
  enddo
end function

! Generate state couplings between a selected list of states.
function selected_states_couplings(input,selected_states) result(output)
  implicit none
  
  class(CoupledStates), intent(in) :: input(:)
  integer,              intent(in) :: selected_states(:)
  type(IntArray1D), allocatable    :: output(:)
  
  integer, allocatable :: state_map(:)
  
  integer :: i,ialloc
  
  ! Generate a list such that state_map(selected_states(i))=i,
  !    and state_map(i)=0 if i is not in selected_states.
  !state_map = [( first(selected_states==i, default=0), &
  !             & i=1,                                  &
  !             & size(input)                           )]
  ! WORKAROUND to avoid memory leak in ifort 19.1.0.166
  allocate(state_map(size(input)), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    state_map(i) = first(selected_states==i, default=0)
  enddo
  
  ! Generate the couplings.
  allocate(output(size(selected_states)), stat=ialloc); call err(ialloc)
  do i=1,size(selected_states)
    output(i)%i = state_map(input(selected_states(i))%ids())
    output(i)%i = output(i)%i(filter(output(i)%i/=0))
  enddo
end function

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
  type(String), allocatable :: coupling(:)
  
  integer :: i,ialloc
  
  select type(this); type is(CoupledStates)
    line = split_line(input)
    allocate( ids(size(line)),         &
            & separations(size(line)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(line)
      coupling = split_line(line(i),delimiter=':')
      ids(i) = int(coupling(1))
      separations(i) = int(coupling(2))
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
