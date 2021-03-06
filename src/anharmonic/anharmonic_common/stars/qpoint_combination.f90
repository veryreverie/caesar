!> Provides the [[QpointCombination(type)]] and [[QpointCombinations(type)]]
!>    classes, and related methods.
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
  public :: QpointCombinations
  public :: generate_qpoint_combinations
  
  !> Records the total power of modes at a set of q-points.
  type, extends(Stringable) :: QpointCombination
    type(QpointPower), allocatable, private :: qpoints_(:)
  contains
    procedure, public :: qpoints => qpoints_QpointCombination
    
    procedure, public :: total_power => total_power_QpointCombination
    
    procedure, public :: wavevector => wavevector_QpointCombination
    
    procedure, public :: complex_monomials => &
                       & complex_monomials_QpointCombination
    
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
  
    !> Returns the total power of modes across the q-point combination.
    impure elemental module function total_power_QpointCombination(this) &
       & result(output)
      class(QpointCombination), intent(in) :: this
      integer                              :: output
    end function
    
    !> Returns the sum of the wavevectors of each mode in the combination.
    module function wavevector_QpointCombination(this,qpoints) result(output)
      class(QpointCombination), intent(in) :: this
      type(QpointData),         intent(in) :: qpoints(:)
      type(FractionVector)                 :: output
    end function
    
    !> Returns the [[ComplexMonomial(type)]]s containing the given `modes`
    !>    which match the q-point combination.
    !> The monomials are generated with coefficients such that symmetry
    !>    operations are unitary.
    module function complex_monomials_QpointCombination(this,modes) &
       & result(output)
      class(QpointCombination), intent(in) :: this
      type(ComplexMode),        intent(in) :: modes(:)
      type(ComplexMonomial), allocatable   :: output(:)
    end function
    
    !> Convert a [[String(type)]] to a [[QpointCombination(type)]].
    module subroutine read_QpointCombination(this,input)
      class(QpointCombination), intent(out) :: this
      type(String),             intent(in)  :: input
    end subroutine
  
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
  
  !> An array of [[QpointCombination(type)]]s with a single `total_power`.
  type, extends(Stringsable) :: QpointCombinations
    integer                              :: power
    type(QpointCombination), allocatable :: combinations(:)
  contains
    procedure, public :: read  => read_QpointCombinations
    procedure, public :: write => write_QpointCombinations
  end type
  
  interface QpointCombinations
    !> Constructor for [[QpointCombinations(type)]] objects.
    module function new_QpointCombinations(power,combinations) result(this)
      integer,                 intent(in) :: power
      type(QpointCombination), intent(in) :: combinations(:)
      type(QpointCombinations)            :: this
    end function
  end interface
  
  interface
    !> Convert a [[String(type)]] array to a [[QpointCombinations(type)]].
    module subroutine read_QpointCombinations(this,input)
      class(QpointCombinations), intent(out) :: this
      type(String),              intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    !> Convert a [[QpointCombinations(type)]] to a [[String(type)]] array.
    module function write_QpointCombinations(this) result(output)
      class(QpointCombinations), intent(in) :: this
      type(String), allocatable             :: output(:)
    end function
  end interface
  
  interface QpointCombinations
    !> Convert a [[String(type)]] array to a [[QpointCombinations(type)]].
    module function new_QpointCombinations_Strings(input) &
       & result(this)
      type(String), intent(in) :: input(:)
      type(QpointCombinations) :: this
    end function
    
    !> Convert a [[StringArray(type)]] to a [[QpointCombinations(type)]].
    impure elemental module function new_QpointCombinations_StringArray( &
       & input) result(this)
      type(StringArray), intent(in) :: input
      type(QpointCombinations)      :: this
    end function
  end interface
  
  interface
    !> Generates all q-point combinations with `total_power` up to `max_power`,
    !>    from a given set of `qpoints`.
    !> Returns an array of [[QpointCombinations(type)]] of length
    !>    `max_power+1`, in ascending order of power.
    !> The combinations at each power are returned in ascending order by '<'.
    module function generate_qpoint_combinations(qpoints,max_power, &
       & max_qpoint_coupling,conserve_momentum) result(output)
      type(QpointData), intent(in)           :: qpoints(:)
      integer,          intent(in)           :: max_power
      !> If `max_qpoint_coupling` is given then only q-point combinations
      !>    containing up to `max_qpoint_coupling` distinct q-points
      !>    are returned.
      !> For the purposes of `max_qpoint_coupling`, a q-point and its pair are
      !>    counted as a single q-point.
      integer,          intent(in), optional :: max_qpoint_coupling
      !> If `conserve_momentum` is `true` then only q-point combinations
      !>    which conserve momentum (i.e. sum q = G) are returned.
      !> `conserve_momentum` defaults to `false`.
      logical,          intent(in), optional :: conserve_momentum
      type(QpointCombinations), allocatable  :: output(:)
    end function
  end interface
end module
