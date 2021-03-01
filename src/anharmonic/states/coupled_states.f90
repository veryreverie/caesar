! ======================================================================
! For a given state |i>, records which states |j> have non-zero <i|H|j>.
! ======================================================================
! This class gets appended to millions of times, so it handles its storage
!    like a C++ vector for speed reasons.
module caesar_coupled_states_module
  use caesar_common_module
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
    ! Constructors and size module function.
    module function new_CoupledStates_null() result(this) 
      type(CoupledStates) :: this
    end function
  
    module function new_CoupledStates(ids,separations) result(this) 
      integer, intent(in) :: ids(:)
      integer, intent(in) :: separations(:)
      type(CoupledStates) :: this
    end function
  end interface
  
  interface size
    module function size_CoupledStates(this) result(output) 
      type(CoupledStates), intent(in) :: this
      integer                         :: output
    end function
  end interface
  
  interface
    ! Getters.
    impure elemental module function id_CoupledStates(this,index) &
       & result(output) 
      class(CoupledStates), intent(in) :: this
      integer,              intent(in) :: index
      integer                          :: output
    end function
  end interface
  
  interface
    module function ids_CoupledStates(this) result(output) 
      class(CoupledStates), intent(in) :: this
      integer, allocatable             :: output(:)
    end function
  end interface
  
  interface
    impure elemental module function separation_CoupledStates(this,index) &
       & result(output) 
      class(CoupledStates), intent(in) :: this
      integer,              intent(in) :: index
      integer                          :: output
    end function
  end interface
  
  interface
    module function separations_CoupledStates(this) result(output) 
      class(CoupledStates), intent(in) :: this
      integer, allocatable             :: output(:)
    end function
  end interface
  
  interface
    ! Add a coupling to a CoupledStates.
    module subroutine add_coupling(this,id,separation) 
      class(CoupledStates), intent(inout) :: this
      integer,              intent(in)    :: id
      integer,              intent(in)    :: separation
    end subroutine
  end interface
  
  interface
    ! Re-map the ids.
    module subroutine map_ids(this,id_map) 
      class(CoupledStates), intent(inout) :: this
      integer,              intent(in)    :: id_map(:)
    end subroutine
  end interface
  
  interface
    ! Take the tensor product of two couplings.
    module function tensor_product_couplings(unmapped_this,this,no_states, &
       & that,max_separation) result(output) 
      type(CoupledStates), intent(in) :: unmapped_this
      type(CoupledStates), intent(in) :: this
      integer,             intent(in) :: no_states(:)
      type(CoupledStates), intent(in) :: that
      integer,             intent(in) :: max_separation
      type(CoupledStates)             :: output
    end function
  end interface
  
  interface
    ! Generate state couplings between a selected list of states.
    module function selected_states_couplings(input,selected_states) &
       & result(output) 
      class(CoupledStates), intent(in) :: input(:)
      integer,              intent(in) :: selected_states(:)
      type(IntArray1D), allocatable    :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_CoupledStates(this,input) 
      class(CoupledStates), intent(out) :: this
      type(String),         intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_CoupledStates(this) result(output) 
      class(CoupledStates), intent(in) :: this
      type(String)                     :: output
    end function
  end interface
  
  interface CoupledStates
    impure elemental module function new_CoupledStates_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(CoupledStates)      :: this
    end function
  end interface
end module
