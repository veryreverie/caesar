!> Provides the [[QpointPower(type)]] class, and related methods.
module caesar_qpoint_power_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: QpointPower
  public :: operator(==)
  public :: operator(/=)
  public :: operator(<)
  public :: operator(<=)
  public :: operator(>)
  public :: operator(>=)
  public :: conjg
  public :: operator(*)
  public :: generate_qpoint_powers
  
  !> Records the total power of modes at a q-point and its pair.
  type, extends(Stringable) :: QpointPower
    !> The `id` of the q-point.
    integer, private :: id_
    !> The total power of modes at the q-point.
    integer, private :: power_
    !> The `id` of the paired q-point.
    integer, private :: paired_id_
    !> The total power of modes at the paired q-point.
    integer, private :: paired_power_
  contains
    procedure, public :: id => id_QpointPower
    procedure, public :: power => power_QpointPower
    procedure, public :: paired_id => paired_id_QpointPower
    procedure, public :: paired_power => paired_power_QpointPower
    
    procedure, public :: total_power => total_power_QpointPower
    
    procedure, public :: read  => read_QpointPower
    procedure, public :: write => write_QpointPower
  end type
  
  interface QpointPower
    !> Construct a [[QpointPower(type)]].
    !> If `paired_id` and `paired_power` default to `id` and `power`
    !>    respectively. `paired_id` and `paired_power` must both be present or
    !>    both be absent.
    !> Reverses `id` and `paired_id` (and `power` and `paired_power`)
    !>    if `paired_id<id`.
    module function new_QpointPower(id,power,paired_id,paired_power) &
       & result(this)
      integer, intent(in)           :: id
      integer, intent(in)           :: power
      integer, intent(in), optional :: paired_id
      integer, intent(in), optional :: paired_power
      type(QpointPower)             :: this
    end function
  end interface
  
  interface
    !> Getter for `this%id_`.
    impure elemental module function id_QpointPower(this) result(output)
      class(QpointPower), intent(in) :: this
      integer                        :: output
    end function
    
    !> Getter for `this%power_`.
    impure elemental module function power_QpointPower(this) result(output)
      class(QpointPower), intent(in) :: this
      integer                        :: output
    end function
  
    !> Getter for `this%paired_id_`.
    impure elemental module function paired_id_QpointPower(this) result(output)
      class(QpointPower), intent(in) :: this
      integer                        :: output
    end function
    
    !> Getter for `this%paired_power_`.
    impure elemental module function paired_power_QpointPower(this) &
       & result(output)
      class(QpointPower), intent(in) :: this
      integer                        :: output
    end function
  end interface
  
  interface
    !> Returns the total power of modes at the q-point pair.
    impure elemental module function total_power_QpointPower(this) &
       & result(output)
      class(QpointPower), intent(in) :: this
      integer                        :: output
    end function
  end interface
  
  interface operator(==)
    !> Equality comparison between two [[QpointPower(type)]] objects.
    impure elemental module function equality_QpointPower_QpointPower(this, &
       & that) result(output)
      type(QpointPower), intent(in) :: this
      type(QpointPower), intent(in) :: that
      logical                       :: output
    end function
  end interface
  
  interface operator(/=)
    !> Non-equality comparison between two [[QpointPower(type)]] objects.
    impure elemental module function non_equality_QpointPower_QpointPower( &
       & this,that) result(output)
      type(QpointPower), intent(in) :: this
      type(QpointPower), intent(in) :: that
      logical                       :: output
    end function
  end interface
  
  interface operator(<)
    !> Less-than comparison between two [[QpointPower(type)]] objects.
    !> Sorts first by `id`, then by `-total_power`, then by `-power`.
    impure elemental module function lt_QpointPower_QpointPower(this,that) &
       & result(output)
      type(QpointPower), intent(in) :: this
      type(QpointPower), intent(in) :: that
      logical                       :: output
    end function
  end interface
  
  interface operator(<=)
    !> Less-than-or-equal comparison between two [[QpointPower(type)]] objects.
    !> Sorts first by `id`, then by `-total_power`, then by `-power`.
    impure elemental module function le_QpointPower_QpointPower(this,that) &
       & result(output)
      type(QpointPower), intent(in) :: this
      type(QpointPower), intent(in) :: that
      logical                       :: output
    end function
  end interface
  
  interface operator(>)
    !> Greater-than comparison between two [[QpointPower(type)]] objects.
    !> Sorts first by `id`, then by `-total_power`, then by `-power`.
    impure elemental module function gt_QpointPower_QpointPower(this,that) &
       & result(output)
      type(QpointPower), intent(in) :: this
      type(QpointPower), intent(in) :: that
      logical                       :: output
    end function
  end interface
  
  interface operator(>=)
    !> Greater-than-or-equal comparison between two [[QpointPower(type)]]
    !>    objects.
    !> Sorts first by `id`, then by `-total_power`, then by `-power`.
    impure elemental module function ge_QpointPower_QpointPower(this,that) &
       & result(output)
      type(QpointPower), intent(in) :: this
      type(QpointPower), intent(in) :: that
      logical                       :: output
    end function
  end interface
  
  interface
    !> Convert a [[String(type)]] to a [[QpointPower(type)]].
    module subroutine read_QpointPower(this,input)
      class(QpointPower), intent(out) :: this
      type(String),       intent(in)  :: input
    end subroutine
  end interface
  
  interface
    !> Convert a [[QpointPower(type)]] to a [[String(type)]].
    module function write_QpointPower(this) result(output)
      class(QpointPower), intent(in) :: this
      type(String)                   :: output
    end function
  end interface
  
  interface QpointPower
    !> Convert a [[String(type)]] to a [[QpointPower(type)]].
    impure elemental module function new_QpointPower_String(input) result(this)
      type(String), intent(in) :: input
      type(QpointPower)        :: this
    end function
  end interface
  
  interface conjg
    !> Returns a [[QpointPower(type)]] with reversed `id` and `paired_id`.
    impure elemental module function conjg_QpointPower(this) result(output)
      type(QpointPower), intent(in) :: this
      type(QpointPower)             :: output
    end function
  end interface
  
  interface operator(*)
    !> The action of a [[Group(type)]] on the q-point indices of a
    !>    [[QpointPower(type)]].
    impure elemental module function operate_Group_QpointPower(qpoint_group, &
       & this) result(output)
      type(Group),       intent(in) :: qpoint_group
      type(QpointPower), intent(in) :: this
      type(QpointPower)             :: output
    end function
  end interface
  
  interface
    !> Generates the set of q-point powers at a given `qpoint`, with
    !>    total_power between `1` and `max_power` inclusive.
    !> Returns output sorted in descending order by '<'.
    module function generate_qpoint_powers(qpoint,max_power) result(output)
      type(QpointData), intent(in)   :: qpoint
      integer,          intent(in)   :: max_power
      type(QpointPower), allocatable :: output(:)
    end function
  end interface
end module
