!> Provides the [[SubspaceCoupling(type)]] class, and related methods.
module caesar_subspace_coupling_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: SubspaceCoupling
  public :: size
  public :: operator(==)
  public :: operator(/=)
  public :: generate_coupled_subspaces
  
  !> A set of [[DegenerateSubspace(type)]]s which are coupled to one another.
  type, extends(Stringable) :: SubspaceCoupling
    !> The `id`s of the degenerate subspaces which are coupled together.
    !> `ids_` is sorted in ascending order.
    integer, allocatable, private :: ids_(:)
  contains
    generic,   public  :: ids => &
                        & ids_SubspaceCoupling, &
                        & ids_SubspaceCoupling_index
    procedure, private :: ids_SubspaceCoupling
    procedure, private :: ids_SubspaceCoupling_index
    
    procedure, public :: subspaces => subspaces_SubspaceCoupling
    
    procedure, public :: remove_subspace
    
    ! I/O.
    procedure, public :: read  => read_SubspaceCoupling
    procedure, public :: write => write_SubspaceCoupling
  end type
  
  interface SubspaceCoupling
    !> Constructor for [[SubspaceCoupling(type)]] objects.
    !> Sorts `ids`, and throws an error if `ids` contains duplicates.
    module function new_SubspaceCoupling(ids) result(output) 
      integer, intent(in)    :: ids(:)
      type(SubspaceCoupling) :: output
    end function
  end interface
  
  interface
    !> A getter for `%ids_`.
    module function ids_SubspaceCoupling(this) result(output)
      class(SubspaceCoupling), intent(in) :: this
      integer, allocatable                :: output(:)
    end function
    
    !> A getter for `%ids_(index)`.
    impure elemental module function ids_SubspaceCoupling_index(this,index) &
       & result(output)
      class(SubspaceCoupling), intent(in) :: this
      integer,                 intent(in) :: index
      integer                             :: output
    end function
  
    !> Returns the [[DegenerateSubspace(type)]]s represented by a
    !>    [[SubspaceCoupling(type)]].
    module function subspaces_SubspaceCoupling(this,subspaces) result(output)
      class(SubspaceCoupling),  intent(in)  :: this
      type(DegenerateSubspace), intent(in)  :: subspaces(:)
      type(DegenerateSubspace), allocatable :: output(:)
    end function
  
    !> Removes the subspace at the given `index` from the coupling.
    module subroutine remove_subspace(this,index)
      class(SubspaceCoupling), intent(inout) :: this
      integer,                 intent(in)    :: index
    end subroutine
  
    !> Convert a [[String(type)]] to a [[SubspaceCoupling(type)]].
    module subroutine read_SubspaceCoupling(this,input) 
      class(SubspaceCoupling), intent(out) :: this
      type(String),            intent(in)  :: input
    end subroutine
  
    !> Convert a [[SubspaceCoupling(type)]] to a [[String(type)]].
    module function write_SubspaceCoupling(this) result(output) 
      class(SubspaceCoupling), intent(in) :: this
      type(String)                        :: output
    end function
  end interface
  
  interface SubspaceCoupling
    !> Convert a [[String(type)]] to a [[SubspaceCoupling(type)]].
    impure elemental module function new_SubspaceCoupling_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(SubspaceCoupling)   :: this
    end function
  end interface
  
  interface size
    !> Returns the number of [[DegenerateSubspace(type)]]s in the
    !>    [[SubspaceCoupling(type)]].
    module function size_SubspaceCoupling(this) result(output) 
      type(SubspaceCoupling), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface operator(==)
    !> Equality between two [[SubspaceCoupling(type)]]s.
    !> Compares `ids`.
    impure elemental module function equality_SubspaceCoupling(this,that) &
       & result(output) 
      type(SubspaceCoupling), intent(in) :: this
      type(SubspaceCoupling), intent(in) :: that
      logical                            :: output
    end function
  end interface
  
  interface operator(/=)
    !> Non-equality between two [[SubspaceCoupling(type)]]s.
    !> Compares `ids`.
    impure elemental module function non_equality_SubspaceCoupling(this,that) &
       & result(output) 
      type(SubspaceCoupling), intent(in) :: this
      type(SubspaceCoupling), intent(in) :: that
      logical                            :: output
    end function
  end interface
  
  interface
    !> Generates the set of all [[SubspaceCombination(type)]]s
    !>    containing subspaces from the list `subspaces`.
    !> The couplings contain between `1` and `max_subspace_coupling`
    !>    separate subspaces.
    !> `max_subspace_coupling` must be at least 1.
    module function generate_coupled_subspaces(subspaces, &
       & max_subspace_coupling) result(output) 
      type(DegenerateSubspace), intent(in) :: subspaces(:)
      integer,                  intent(in) :: max_subspace_coupling
      type(SubspaceCoupling), allocatable  :: output(:)
    end function
  end interface
end module
