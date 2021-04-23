!> Provides the [[QpointCombination(type)]] class, and related methods.
module caesar_qpoint_combination_module
  use caesar_common_module
  
  use caesar_qpoint_power_module
  implicit none
  
  private
  
  public :: QpointCombination
  public :: operator(==)
  public :: operator(/=)
  
  !> Records the total power of modes at a combination of q-points.
  type, extends(Stringable) :: QpointCombination
    type(QpointPower), allocatable, private :: qpoints_(:)
  contains
    procedure, public :: qpoints => qpoints_QpointCombination
    
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
end module
