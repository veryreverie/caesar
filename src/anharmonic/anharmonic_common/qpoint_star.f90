!> Provides the [[QpointStar(type)]] class, and related methods.
module caesar_qpoint_star_module
  use caesar_common_module
  
  use caesar_qpoint_power_module
  use caesar_qpoint_combination_module
  implicit none
  
  private
  
  public :: QpointStar
  public :: operator(==)
  public :: operator(/=)
  public :: generate_qpoint_stars
  
  !> Records the total power of modes at a combination of q-points.
  type, extends(Stringsable) :: QpointStar
    type(QpointCombination), allocatable, private :: combinations_(:)
  contains
    procedure, public :: combinations => combinations_QpointStar
    
    procedure, public :: total_power => total_power_QpointStar
    
    procedure, public :: read  => read_QpointStar
    procedure, public :: write => write_QpointStar
  end type
  
  interface QpointStar
    !> Construct a [[QpointStar(type)]].
    !> `combinations` must contain at least one element,
    !>    and no duplicate elements.
    !> Sorts `combinations`.
    module function new_QpointStar(combinations) result(this)
      type(QpointCombination), intent(in) :: combinations(:)
      type(QpointStar)                    :: this
    end function
  end interface
  
  interface
    !> Getter for `this%combinations_`.
    module function combinations_QpointStar(this) result(output)
      class(QpointStar), intent(in)        :: this
      type(QpointCombination), allocatable :: output(:)
    end function
  end interface
  
  interface
    !> Returns the total power of modes at each combination in the star.
    !> Assumes that the total power is the same for all combinations.
    impure elemental module function total_power_QpointStar(this) &
       & result(output)
      class(QpointStar), intent(in) :: this
      integer                       :: output
    end function
  end interface
  
  interface operator(==)
    !> Equality comparison between two [[QpointStar(type)]] objects.
    !> Assumes that the two stars are either the same or mutually exclusive,
    !>    so only compares the first element of `combinations`.
    impure elemental module function &
       & equality_QpointStar_QpointStar(this,that) result(output)
      type(QpointStar), intent(in) :: this
      type(QpointStar), intent(in) :: that
      logical                      :: output
    end function
  end interface
  
  interface operator(/=)
    !> Non-equality comparison between two [[QpointStar(type)]] objects.
    !> Assumes that the two stars are either the same or mutually exclusive,
    !>    so only compares the first element of `combinations`.
    impure elemental module function                                 &
       & non_equality_QpointStar_QpointStar(this,that) &
       & result(output)
      type(QpointStar), intent(in) :: this
      type(QpointStar), intent(in) :: that
      logical                      :: output
    end function
  end interface
  
  interface
    !> Convert a [[String(type)]] array to a [[QpointStar(type)]].
    module subroutine read_QpointStar(this,input)
      class(QpointStar), intent(out) :: this
      type(String),      intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    !> Convert a [[QpointStar(type)]] to a [[String(type)]] array.
    module function write_QpointStar(this) result(output)
      class(QpointStar), intent(in) :: this
      type(String), allocatable     :: output(:)
    end function
  end interface
  
  interface QpointStar
    !> Convert a [[String(type)]] array to a [[QpointStar(type)]].
    module function new_QpointStar_Strings(input) &
       & result(this)
      type(String), intent(in) :: input(:)
      type(QpointStar)         :: this
    end function
    
    !> Convert a [[StringArray(type)]] to a [[QpointStar(type)]].
    impure elemental module function new_QpointStar_StringArray(input) &
      & result(this)
      type(StringArray), intent(in) :: input
      type(QpointStar)              :: this
    end function
  end interface
  
  interface
    !> Generates the set of q-point stars with a given `power` from a given set
    !>    of `qpoints`, based on how these q-points transform under symmetry,
    !>    as described by `qpoint_groups`.
    !> If `conserve_momentum` is `true` then only stars of q-point combinations
    !>    which conserve momentum (i.e. sum q = G) are returned.
    !> `conserve_momentum` defaults to `false`.
    module function generate_qpoint_stars(qpoints,qpoint_groups,power, &
       & conserve_momentum) result(output)
      type(QpointData),       intent(in)           :: qpoints(:)
      type(Group),            intent(in)           :: qpoint_groups(:)
      integer,                intent(in)           :: power
      logical,                intent(in), optional :: conserve_momentum
      type(QpointStar), allocatable                :: output(:)
    end function
  end interface
end module
