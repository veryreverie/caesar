! ======================================================================
! A single R-vector at which a subspace should be sampled.
! ======================================================================
module vscf_rvector_module
  use common_module
  implicit none
  
  private
  
  public :: VscfRvector
  
  type, extends(Stringable) :: VscfRvector
    integer         :: subspace_id
    type(IntVector) :: rvector
  contains
    procedure, public :: read  => read_VscfRvector
    procedure, public :: write => write_VscfRvector
  end type
  
  interface VscfRvector
    module procedure new_VscfRvector
    module procedure new_VscfRvector_String
  end interface
contains

! Constructor.
function new_VscfRvector(subspace_id,rvector) result(this)
  implicit none
  
  integer,         intent(in) :: subspace_id
  type(IntVector), intent(in) :: rvector
  type(VscfRvector)           :: this
  
  this%subspace_id = subspace_id
  this%rvector     = rvector
end function

! I/O.
subroutine read_VscfRvector(this,input)
  implicit none
  
  class(VscfRvector), intent(out) :: this
  type(String),       intent(in)  :: input
  
  integer         :: subspace_id
  type(IntVector) :: rvector
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(VscfRvector)
    line = split_line(input)
    subspace_id = int(line(2))
    rvector = vec(int(line(4:6)))
    
    this = VscfRvector(subspace_id,rvector)
  end select
end subroutine

function write_VscfRvector(this) result(output)
  implicit none
  
  class(VscfRvector), intent(in) :: this
  type(String)                   :: output
  
  select type(this); type is(VscfRvector)
    output = 'Subspace: '//this%subspace_id//' R-vector: '//this%rvector
  end select
end function

impure elemental function new_VscfRvector_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(VscfRvector)        :: this
  
  call this%read(input)
end function
end module
