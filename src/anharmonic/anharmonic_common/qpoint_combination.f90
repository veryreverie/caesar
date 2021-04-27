!> Provides the [[QpointCombination(type)]] class, and related methods.
module caesar_qpoint_combination_module
  use caesar_common_module
  
  use caesar_qpoint_power_module
  implicit none
  
  private
  
  public :: QpointCombination
  public :: operator(==)
  public :: operator(/=)
  public :: operator(<)
  public :: operator(<=)
  public :: operator(>)
  public :: operator(>=)
  public :: conjg
  public :: operator(*)
  public :: generate_qpoint_combinations
  
  !> Records the total power of modes at a combination of q-points.
  type, extends(Stringable) :: QpointCombination
    type(QpointPower), allocatable, private :: qpoints_(:)
  contains
    procedure, public :: qpoints => qpoints_QpointCombination
    
    procedure, public :: total_power => total_power_QpointCombination
    
    procedure, public :: read  => read_QpointCombination
    procedure, public :: write => write_QpointCombination
  end type
  
  interface QpointCombination
    !> Construct a [[QpointCombination(type)]].
    !> Sorts `qpoints` by `id`.
    module function new_QpointCombination(qpoints) result(this)
      type(QpointPower), intent(in) :: qpoints(:)
      type(QpointCombination)       :: this
    end function
  end interface
  
  interface
    !> Getter for `this%qpoints_`.
    module function qpoints_QpointCombination(this) result(output)
      class(QpointCombination), intent(in) :: this
      type(QpointPower), allocatable       :: output(:)
    end function
  end interface
  
  interface
    !> Returns the total power of modes across the q-point combination.
    impure elemental module function total_power_QpointCombination(this) &
       & result(output)
      class(QpointCombination), intent(in) :: this
      integer                              :: output
    end function
  end interface
  
  interface operator(==)
    !> Equality comparison between two [[QpointCombination(type)]] objects.
    impure elemental module function &
       & equality_QpointCombination_QpointCombination(this,that) result(output)
      type(QpointCombination), intent(in) :: this
      type(QpointCombination), intent(in) :: that
      logical                             :: output
    end function
  end interface
  
  interface operator(/=)
    !> Non-equality comparison between two [[QpointCombination(type)]] objects.
    impure elemental module function                                 &
       & non_equality_QpointCombination_QpointCombination(this,that) &
       & result(output)
      type(QpointCombination), intent(in) :: this
      type(QpointCombination), intent(in) :: that
      logical                             :: output
    end function
  end interface
  
  interface operator(<)
    !> Less-than comparison between two [[QpointCombination(type)]] objects.
    !> Sorts first by `total_power`, then by comparing the q-point powers
    !>    in order.
    impure elemental module function lt_QpointCombination_QpointCombination( &
       & this,that) result(output)
      type(QpointCombination), intent(in) :: this
      type(QpointCombination), intent(in) :: that
      logical                             :: output
    end function
  end interface
  
  interface operator(<=)
    !> Less-than-or-equal comparison between two [[QpointCombination(type)]]
    !>    objects.
    !> Sorts first by `total_power`, then by comparing the q-point powers
    !>    in order.
    impure elemental module function le_QpointCombination_QpointCombination( &
       & this,that) result(output)
      type(QpointCombination), intent(in) :: this
      type(QpointCombination), intent(in) :: that
      logical                             :: output
    end function
  end interface
  
  interface operator(>)
    !> Greater-than comparison between two [[QpointCombination(type)]] objects.
    !> Sorts first by `total_power`, then by comparing the q-point powers
    !>    in order.
    impure elemental module function gt_QpointCombination_QpointCombination( &
       & this,that) result(output)
      type(QpointCombination), intent(in) :: this
      type(QpointCombination), intent(in) :: that
      logical                             :: output
    end function
  end interface
  
  interface operator(>=)
    !> Greater-than-or-equal comparison between two [[QpointCombination(type)]]
    !>    objects.
    !> Sorts first by `total_power`, then by comparing the q-point powers
    !>    in order.
    impure elemental module function ge_QpointCombination_QpointCombination( &
       & this,that) result(output)
      type(QpointCombination), intent(in) :: this
      type(QpointCombination), intent(in) :: that
      logical                             :: output
    end function
  end interface
  
  interface
    !> Convert a [[String(type)]] to a [[QpointCombination(type)]].
    module subroutine read_QpointCombination(this,input)
      class(QpointCombination), intent(out) :: this
      type(String),             intent(in)  :: input
    end subroutine
  end interface
  
  interface
    !> Convert a [[QpointCombination(type)]] to a [[String(type)]].
    module function write_QpointCombination(this) result(output)
      class(QpointCombination), intent(in) :: this
      type(String)                         :: output
    end function
  end interface
  
  interface QpointCombination
    !> Convert a [[String(type)]] to a [[QpointCombination(type)]].
    impure elemental module function new_QpointCombination_String(input) &
       & result(this)
      type(String), intent(in) :: input
      type(QpointCombination)  :: this
    end function
  end interface
  
  interface conjg
    !> Returns a [[QpointCombination(type)]] where the conjugate has been taken
    !>    of each [[QpointPower(type)]].
    impure elemental module function conjg_QpointCombination(this) &
       & result(output)
      type(QpointCombination), intent(in) :: this
      type(QpointCombination)             :: output
    end function
  end interface
  
  interface operator(*)
    !> The action of a [[Group(type)]] on the q-point indices of a
    !>    [[QpointCombination(type)]].
    impure elemental module function operate_Group_QpointCombination( &
       & qpoint_group,this) result(output)
      type(Group),             intent(in) :: qpoint_group
      type(QpointCombination), intent(in) :: this
      type(QpointCombination)             :: output
    end function
  end interface
  
  interface
    !> Generates all q-point combinations with a given `power` from a given set
    !>    of `qpoints`.
    !> The combinations are returned in ascending order by '<'.
    !> If `conserve_momentum` is `true` then only q-point combinations
    !>    which conserve momentum (i.e. sum q = G) are returned.
    !> `conserve_momentum` defaults to `false`.
    module function generate_qpoint_combinations(qpoints,power, &
       & conserve_momentum) result(output)
      type(QpointData), intent(in)           :: qpoints(:)
      integer,          intent(in)           :: power
      logical,          intent(in), optional :: conserve_momentum
      type(QpointCombination), allocatable   :: output(:)
    end function
  end interface
end module
