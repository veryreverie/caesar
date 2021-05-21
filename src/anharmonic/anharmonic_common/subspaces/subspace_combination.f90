!> Provides the [[SubspaceCombination(type)]] class, and related methods.
module caesar_subspace_combination_module
  use caesar_common_module
  
  use caesar_subspace_coupling_module
  implicit none
  
  private
  
  public :: SubspaceCombination
  public :: size
  public :: operator(==)
  public :: operator(/=)
  public :: operator(<)
  public :: operator(<=)
  public :: operator(>)
  public :: operator(>=)
  public :: generate_subspace_combinations
  
  !> Records the total power of modes in a set of [[DegenerateSubspace]]s.
  type, extends(Stringable) :: SubspaceCombination
    !> The `id`s of the subspaces with non-zero mode powers.
    !> Stored in ascending order.
    integer, allocatable, private :: ids_(:)
    !> `powers_(i)` is the total power of the modes in the subspace with `id`
    !>    equal to `ids_(i)`.
    integer, allocatable, private :: powers_(:)
  contains
    generic,   public  :: ids =>                   &
                        & ids_SubspaceCombination, &
                        & ids_SubspaceCombination_index
    procedure, private :: ids_SubspaceCombination
    procedure, private :: ids_SubspaceCombination_index
    
    generic,   public  :: powers =>                   &
                        & powers_SubspaceCombination, &
                        & powers_SubspaceCombination_index
    procedure, private :: powers_SubspaceCombination
    procedure, private :: powers_SubspaceCombination_index
    
    procedure, public :: subspaces => subspaces_SubspaceCombination
    
    procedure, public :: complex_monomials => &
                       & complex_monomials_SubspaceCombination
    procedure, public :: paired_monomials => &
                       & paired_monomials_SubspaceCombination
    
    procedure, public :: is_subsidiary_of
    
    ! I/O.
    procedure, public :: read  => read_SubspaceCombination
    procedure, public :: write => write_SubspaceCombination
  end type
  
  interface SubspaceCombination
    !> Constructor for objects of type [[SubspaceCombination(type)]].
    module function new_SubspaceCombination(ids,powers) result(this) 
      integer, intent(in)       :: ids(:)
      integer, intent(in)       :: powers(:)
      type(SubspaceCombination) :: this
    end function
  end interface
  
  interface
    !> Getter for `this%_ids`.
    module function ids_SubspaceCombination(this) result(output)
      class(SubspaceCombination), intent(in) :: this
      integer, allocatable                   :: output(:)
    end function
    
    !> Getter for `this%_ids(index)`.
    impure elemental module function ids_SubspaceCombination_index(this, &
       & index) result(output)
      class(SubspaceCombination), intent(in) :: this
      integer,                    intent(in) :: index
      integer                                :: output
    end function
    
    !> Getter for `this%_powers`.
    module function powers_SubspaceCombination(this) result(output)
      class(SubspaceCombination), intent(in) :: this
      integer, allocatable                   :: output(:)
    end function
    
    !> Getter for `this%_powers(index)`.
    impure elemental module function powers_SubspaceCombination_index(this, &
       & index) result(output)
      class(SubspaceCombination), intent(in) :: this
      integer,                    intent(in) :: index
      integer                                :: output
    end function
  
    !> Returns the [[DegenerateSubspace(type)]]s with non-zero power in the
    !>    combination.
    module function subspaces_SubspaceCombination(this,subspaces) &
       & result(output) 
      class(SubspaceCombination), intent(in) :: this
      type(DegenerateSubspace),   intent(in) :: subspaces(:)
      type(DegenerateSubspace), allocatable  :: output(:)
    end function
  
    !> Generate the [[ComplexMonomial(type)]]s corresponding to a given
    !>    SubspaceCombination.
    !> The coefficients are chosen such that symmetry operations are unitary
    !>    in the basis of monomials.
    module function complex_monomials_SubspaceCombination(this,            &
       & maximum_coupling_order,subspaces,modes,qpoints,conserve_momentum, &
       & conserve_subspace_momentum) result(output)
      class(SubspaceCombination), intent(in) :: this
      integer,                    intent(in) :: maximum_coupling_order
      type(DegenerateSubspace),   intent(in) :: subspaces(:)
      type(ComplexMode),          intent(in) :: modes(:)
      type(QpointData),           intent(in) :: qpoints(:)
      logical,                    intent(in) :: conserve_momentum
      logical,                    intent(in) :: conserve_subspace_momentum
      type(ComplexMonomial), allocatable     :: output(:)
    end function
  
    !> Generate the [[PairedMonomial(type)]]s corresponding to a given
    !>    SubspaceCombination.
    !> The coefficients are chosen such that symmetry operations are unitary
    !>    in the basis of monomials.
    module function paired_monomials_SubspaceCombination(this,             &
       & maximum_coupling_order,subspaces,modes,qpoints,conserve_momentum, &
       & conserve_subspace_momentum) result(output)
      class(SubspaceCombination), intent(in) :: this
      integer,                    intent(in) :: maximum_coupling_order
      type(DegenerateSubspace),   intent(in) :: subspaces(:)
      type(ComplexMode),          intent(in) :: modes(:)
      type(QpointData),           intent(in) :: qpoints(:)
      logical,                    intent(in) :: conserve_momentum
      logical,                    intent(in) :: conserve_subspace_momentum
      type(PairedMonomial), allocatable      :: output(:)
    end function
  
    !> Check if this is subsidiary to the given combination,
    !>    e.g. the coupling (s1^2*s1^1) has subsidiary combinations:
    !>    (), (s1^1), (s2^1), (s1^2), (s1^1*s2^1) and (s1^2*s1^1).
    impure elemental module function is_subsidiary_of(this,that) &
       & result(output) 
      Class(SubspaceCombination), intent(in) :: this
      type(SubspaceCombination),  intent(in) :: that
      logical                                :: output
    end function
  
    !> Convert a [[String(type)]] to a [[SubspaceCombination(type)]].
    module subroutine read_SubspaceCombination(this,input) 
      class(SubspaceCombination), intent(out) :: this
      type(String),               intent(in)  :: input
    end subroutine
  
    !> Convert a [[SubspaceCombination(type)]] to a [[String(type)]].
    module function write_SubspaceCombination(this) result(output) 
      class(SubspaceCombination), intent(in) :: this
      type(String)                           :: output
    end function
  end interface
  
  interface SubspaceCombination
    !> Convert a [[String(type)]] to a [[SubspaceCombination(type)]].
    impure elemental module function new_SubspaceCombination_String(input) &
       & result(this) 
      type(String), intent(in)  :: input
      type(SubspaceCombination) :: this
    end function
  end interface
  
  interface size
    !> The number of subspaces with non-zero powers in the combination.
    !> This is equal to `size(this%ids)` and `size(this%powers)`.
    module function size_SubspaceCombination(this) result(output) 
      type(SubspaceCombination), intent(in) :: this
      integer                               :: output
    end function
  end interface
  
  interface operator(==)
    !> Equality comparison between two [[SubspaceCombination(type)]]s.
    !> Compares `ids` and `powers`.
    impure elemental module function                                 &
       & equality_SubspaceCombination_SubspaceCombination(this,that) &
       & result(output) 
      type(SubspaceCombination), intent(in) :: this
      type(SubspaceCombination), intent(in) :: that
      logical                               :: output
    end function
  end interface
  
  interface operator(/=)
    !> Non-equality comparison between two [[SubspaceCombination(type)]]s.
    !> Compares `ids` and `powers`.
    impure elemental module function                                     &
       & non_equality_SubspaceCombination_SubspaceCombination(this,that) &
       & result(output) 
      type(SubspaceCombination), intent(in) :: this
      type(SubspaceCombination), intent(in) :: that
      logical                               :: output
    end function
  end interface
  
  interface operator(<)
    !> Less-than comparison between two [[SubspaceCombination(type)]] objects.
    !> Sorts first by `sum(powers)`,
    !>    then by `ids(1)`, then by `-powers(1)`,
    !>    then by `ids(2)`, then by `-powers(2)` etc.
    impure elemental module function &
       & lt_SubspaceCombination_SubspaceCombination(this,that) result(output)
      type(SubspaceCombination), intent(in) :: this
      type(SubspaceCombination), intent(in) :: that
      logical                               :: output
    end function
  end interface
  
  interface operator(<=)
    !> Less-than-or-equal comparison between two [[SubspaceCombination(type)]]
    !>    objects.
    !> Sorts first by `sum(powers)`,
    !>    then by `ids(1)`, then by `-powers(1)`,
    !>    then by `ids(2)`, then by `-powers(2)` etc.
    impure elemental module function &
       & le_SubspaceCombination_SubspaceCombination(this,that) result(output)
      type(SubspaceCombination), intent(in) :: this
      type(SubspaceCombination), intent(in) :: that
      logical                               :: output
    end function
  end interface
  
  interface operator(>)
    !> Greater-than comparison between two [[SubspaceCombination(type)]]
    !>    objects.
    !> Sorts first by `sum(powers)`,
    !>    then by `ids(1)`, then by `-powers(1)`,
    !>    then by `ids(2)`, then by `-powers(2)` etc.
    impure elemental module function &
       & gt_SubspaceCombination_SubspaceCombination(this,that) result(output)
      type(SubspaceCombination), intent(in) :: this
      type(SubspaceCombination), intent(in) :: that
      logical                               :: output
    end function
  end interface
  
  interface operator(>=)
    !> Greater-than-or-equal comparison between two
    !>    [[SubspaceCombination(type)]] objects.
    !> Sorts first by `sum(powers)`,
    !>    then by `ids(1)`, then by `-powers(1)`,
    !>    then by `ids(2)`, then by `-powers(2)` etc.
    impure elemental module function &
       & ge_SubspaceCombination_SubspaceCombination( &
       & this,that) result(output)
      type(SubspaceCombination), intent(in) :: this
      type(SubspaceCombination), intent(in) :: that
      logical                               :: output
    end function
  end interface
  
  interface
    !> Generates all [[CouplingCombination(type)]]s corresponding to a given
    !>    [[SubspaceCoupling(type)]], with total powers between
    !>    `minimum_power` and `maximum_power`.
    !> e.g. given the subspace coupling (s3*s5) and expansion order 3,
    !>    will return subspace combinations (s3^2)*(s5^1) and (s3^1)*(s5^2).
    module function generate_subspace_combinations(subspace_coupling, &
       & minimum_power,maximum_power) result(output)
      type(SubspaceCoupling),   intent(in)    :: subspace_coupling
      integer,                  intent(in)    :: minimum_power
      integer,                  intent(in)    :: maximum_power
      type(SubspaceCombination), allocatable  :: output(:)
    end function
  end interface
end module
