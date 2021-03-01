! ======================================================================
! A single R-vector at which a subspace should be sampled.
! ======================================================================
module caesar_vscf_rvector_module
  use caesar_common_module
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
    ! Constructor.
    module function new_VscfRvector(subspace_id,rvector) result(this) 
      integer,         intent(in) :: subspace_id
      type(IntVector), intent(in) :: rvector
      type(VscfRvector)           :: this
    end function
  end interface
  
  interface
    ! I/O.
    module subroutine read_VscfRvector(this,input) 
      class(VscfRvector), intent(out) :: this
      type(String),       intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_VscfRvector(this) result(output) 
      class(VscfRvector), intent(in) :: this
      type(String)                   :: output
    end function
  end interface
  
  interface VscfRvector
    impure elemental module function new_VscfRvector_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(VscfRvector)        :: this
    end function
  end interface
end module
