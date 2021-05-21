!> Provides the [[CombinationQpointCombination(type)]] class,
!>    and related methods.
module caesar_combination_qpoint_combination_module
  use caesar_common_module
  
  use caesar_subspaces_module
  
  use caesar_qpoint_combination_module
  use caesar_qpoint_star_product_module
  implicit none
  
  private
  
  public :: CombinationQpointCombination
  public :: operator(==)
  public :: operator(/=)
  public :: operator(<)
  public :: operator(<=)
  public :: operator(>)
  public :: operator(>=)
  public :: conjg
  public :: operator(*)
  public :: generate_combination_qpoint_combinations
  
  !> Records the total power of modes at a set of q-points,
  !>    across a subspace combination.
  type, extends(Stringsable) :: CombinationQpointCombination
    type(SubspaceCombination)            :: subspace_combination
    type(QpointCombination), allocatable :: qpoint_combinations(:)
  contains
    procedure, public :: total_power => &
                       & total_power_CombinationQpointCombination
    
    procedure, public :: wavevector => &
                       & wavevector_CombinationQpointCombination
    
    procedure, public :: complex_monomials => &
                       & complex_monomials_CombinationQpointCombination
    
    procedure, public :: read  => read_CombinationQpointCombination
    procedure, public :: write => write_CombinationQpointCombination
  end type
  
  interface CombinationQpointCombination
    !> Construct a [[CombinationQpointCombination(type)]].
    module function new_CombinationQpointCombination(subspace_combination, &
       & qpoint_combinations) result(this)
      type(SubspaceCombination), intent(in) :: subspace_combination
      type(QpointCombination),   intent(in) :: qpoint_combinations(:)
      type(CombinationQpointCombination)    :: this
    end function
  end interface
  
  interface
    !> Returns the total power of modes across the combination.
    impure elemental module function &
       & total_power_CombinationQpointCombination(this) result(output)
      class(CombinationQpointCombination), intent(in) :: this
      integer                                         :: output
    end function
    
    !> Returns the sum of the wavevectors across the combination.
    module function wavevector_CombinationQpointCombination(this,qpoints) &
       & result(output)
      class(CombinationQpointCombination), intent(in) :: this
      type(QpointData),                    intent(in) :: qpoints(:)
      type(FractionVector)                            :: output
    end function
    
    !> Returns the [[ComplexMonomial(type)]]s containing the given `modes`
    !>    which match the subspace and q-point combinations.
    module function complex_monomials_CombinationQpointCombination(this, &
       & modes) result(output)
      class(CombinationQpointCombination), intent(in) :: this
      type(ComplexMode),                   intent(in) :: modes(:)
      type(ComplexMonomial), allocatable              :: output(:)
    end function
    
    !> Convert a [[String(type)]] to a [[CombinationQpointCombination(type)]].
    module subroutine read_CombinationQpointCombination(this,input)
      class(CombinationQpointCombination), intent(out) :: this
      type(String),                        intent(in)  :: input(:)
    end subroutine
  
    !> Convert a [[CombinationQpointCombination(type)]] to a [[String(type)]].
    module function write_CombinationQpointCombination(this) result(output)
      class(CombinationQpointCombination), intent(in) :: this
      type(String), allocatable                       :: output(:)
    end function
  end interface
  
  interface CombinationQpointCombination
    !> Convert a [[String(type)]] array to a
    !>    [[CombinationQpointCombination(type)]].
    module function new_CombinationQpointCombination_Strings(input) &
       & result(this)
      type(String), intent(in)           :: input(:)
      type(CombinationQpointCombination) :: this
    end function
    
    !> Convert a [[StringArray(type)]] to a [[QpointStarProduct(type)]].
    impure elemental module function                   &
       & new_CombinationQpointCombination_StringArray( &
       & input) result(this)
      type(StringArray), intent(in)      :: input
      type(CombinationQpointCombination) :: this
    end function
  end interface
  
  interface operator(==)
    !> Equality comparison between two [[CombinationQpointCombination(type)]]
    !>    objects.
    impure elemental module function                                   &
       & eq_CombinationQpointCombination_CombinationQpointCombination( &
       & this,that) result(output)
      type(CombinationQpointCombination), intent(in) :: this
      type(CombinationQpointCombination), intent(in) :: that
      logical                                       :: output
    end function
  end interface
  
  interface operator(/=)
    !> Non-equality comparison between two
    !>    [[CombinationQpointCombination(type)]] objects.
    impure elemental module function                                   &
       & ne_CombinationQpointCombination_CombinationQpointCombination( &
       & this,that) result(output)
      type(CombinationQpointCombination), intent(in) :: this
      type(CombinationQpointCombination), intent(in) :: that
      logical                                        :: output
    end function
  end interface
  
  interface operator(<)
    !> Less-than comparison between two [[CombinationQpointCombination(type)]]
    !>    objects.
    !> Sorts by `subspace_combination`,
    !>    then by `qpoint_combinations(1)`,
    !>    then by `qpoint_combinations(2)` etc.
    impure elemental module function                                   &
       & lt_CombinationQpointCombination_CombinationQpointCombination( &
       & this,that) result(output)
      type(CombinationQpointCombination), intent(in) :: this
      type(CombinationQpointCombination), intent(in) :: that
      logical                                        :: output
    end function
  end interface
  
  interface operator(<=)
    !> Less-than-or-equal comparison between two
    !>    [[CombinationQpointCombination(type)]] objects.
    !> Sorts by `subspace_combination`,
    !>    then by `qpoint_combinations(1)`,
    !>    then by `qpoint_combinations(2)` etc.
    impure elemental module function                                   &
       & le_CombinationQpointCombination_CombinationQpointCombination( &
       & this,that) result(output)
      type(CombinationQpointCombination), intent(in) :: this
      type(CombinationQpointCombination), intent(in) :: that
      logical                                        :: output
    end function
  end interface
  
  interface operator(>)
    !> Greater-than comparison between two
    !>    [[CombinationQpointCombination(type)]] objects.
    !> Sorts by `subspace_combination`,
    !>    then by `qpoint_combinations(1)`,
    !>    then by `qpoint_combinations(2)` etc.
    impure elemental module function                                   &
       & gt_CombinationQpointCombination_CombinationQpointCombination( &
       & this,that) result(output)
      type(CombinationQpointCombination), intent(in) :: this
      type(CombinationQpointCombination), intent(in) :: that
      logical                                        :: output
    end function
  end interface
  
  interface operator(>=)
    !> Greater-than-or-equal comparison between two
    !>    [[CombinationQpointCombination(type)]] objects.
    !> Sorts by `subspace_combination`,
    !>    then by `qpoint_combinations(1)`,
    !>    then by `qpoint_combinations(2)` etc.
    impure elemental module function &
       & ge_CombinationQpointCombination_CombinationQpointCombination( &
       & this,that) result(output)
      type(CombinationQpointCombination), intent(in) :: this
      type(CombinationQpointCombination), intent(in) :: that
      logical                                        :: output
    end function
  end interface
  
  interface conjg
    !> Returns a [[CombinationQpointCombination(type)]] where the conjugate has
    !>    been taken of each [[QpointCombination(type)]].
    !> N.B. this does not affect the subspace combination.
    impure elemental module function conjg_CombinationQpointCombination(this) &
       & result(output)
      type(CombinationQpointCombination), intent(in) :: this
      type(CombinationQpointCombination)             :: output
    end function
  end interface
  
  interface operator(*)
    !> The action of a [[Group(type)]] on the q-point indices of a
    !>    [[CombinationQpointCombination(type)]].
    impure elemental module function                 &
       & operate_Group_CombinationQpointCombination( &
       & qpoint_group,this) result(output)
      type(Group),                        intent(in) :: qpoint_group
      type(CombinationQpointCombination), intent(in) :: this
      type(CombinationQpointCombination)             :: output
    end function
  end interface
  
  interface
    !> Generates the [[CombinationQpointCombination(type)]]s corresponding to a
    !>    [[QpointStarProduct(type)]].
    module function generate_combination_qpoint_combinations( &
       & qpoint_star_product,conserve_momentum,qpoints) result(output)
      type(QpointStarProduct), intent(in)             :: qpoint_star_product
      !> If `conserve_momentum` is `.true.` then only q-point combinations
      !>    which conserve momentum (i.e. sum q=G) are generated.
      !> `conserve_momentum` defaults to `.false.`.
      logical,                 intent(in), optional   :: conserve_momentum
      !> `qpoints` is only required if `conserve_momentum` is `.true.`.
      !> `qpoints` must contain every q-point in `qpoint_star_product`.
      type(QpointData),        intent(in), optional   :: qpoints(:)
      type(CombinationQpointCombination), allocatable :: output(:)
    end function
  end interface
end module
