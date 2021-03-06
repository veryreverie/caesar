!> Provides the [[CombinationQpointStar(type)]] class, and related methods.
module caesar_combination_qpoint_star_module
  use caesar_common_module
  
  use caesar_subspaces_module
  
  use caesar_qpoint_combination_module
  use caesar_qpoint_star_product_module
  use caesar_combination_qpoint_combination_module
  implicit none
  
  private
  
  public :: CombinationQpointStar
  public :: generate_combination_qpoint_stars
  
  !> A star of [[CombinationQpointCombination(type)]]s corresponding to a
  !>    [[SubspaceCombination(type)]].
  type, extends(Stringsable) :: CombinationQpointStar
    type(SubspaceCombination)                       :: subspace_combination
    type(CombinationQpointCombination), allocatable :: qpoint_combinations(:)
  contains
    procedure, public :: complex_monomials => &
                       & complex_monomials_CombinationQpointStar
    
    procedure, public :: read  => read_CombinationQpointStar
    procedure, public :: write => write_CombinationQpointStar
  end type
  
  interface CombinationQpointStar
    !> Constructor for [[CombinationQpointStar(type)]] objects.
    module function new_CombinationQpointStar(qpoint_combinations) result(this)
      type(CombinationQpointCombination), intent(in) :: qpoint_combinations(:)
      type(CombinationQpointStar)                    :: this
    end function
  end interface
  
  interface
    !> Returns the [[ComplexMonomial(type)]]s containing the given `modes`
    !>    which match the subspace and q-point combinations in the star.
    module function complex_monomials_CombinationQpointStar(this,modes) &
       & result(output)
      class(CombinationQpointStar), intent(in) :: this
      type(ComplexMode),            intent(in) :: modes(:)
      type(ComplexMonomial), allocatable       :: output(:)
    end function
    
    !> Convert a [[String(type)]] array to a [[CombinationQpointStar(type)]].
    module subroutine read_CombinationQpointStar(this,input)
      class(CombinationQpointStar), intent(out) :: this
      type(String),                 intent(in)  :: input(:)
    end subroutine
    
    !> Convert a [[CombinationQpointStar(type)]] to a [[String(type)]] array.
    module function write_CombinationQpointStar(this) result(output)
      class(CombinationQpointStar), intent(in) :: this
      type(String), allocatable                :: output(:)
    end function
  end interface
  
  interface CombinationQpointStar
    !> Convert a [[String(type)]] array to a [[CombinationQpointStar(type)]].
    module function new_CombinationQpointStar_Strings(input) &
       & result(this)
      type(String), intent(in)    :: input(:)
      type(CombinationQpointStar) :: this
    end function
    
    !> Convert a [[StringArray(type)]] to a [[CombinationQpointStar(type)]].
    impure elemental module function new_CombinationQpointStar_StringArray( &
       & input) result(this)
      type(StringArray), intent(in) :: input
      type(CombinationQpointStar)   :: this
    end function
  end interface
  
  interface
    !> Generates the set of [[CombinationQpointStar(type)]]s for a given
    !>    [[QpointStarProduct]], based on how the [[QpointCombination(type)]]s
    !>    transform under symmetry, as described by `qpoint_groups`.
    module function generate_combination_qpoint_stars(qpoint_star_product, &
       & qpoint_groups,max_qpoint_coupling,conserve_momentum,qpoints)      &
       & result(output)
      type(QpointStarProduct), intent(in)           :: qpoint_star_product
      type(Group),             intent(in)           :: qpoint_groups(:)
      !> If `max_qpoint_coupling` is given then only q-point combinations
      !>    containing up to `max_qpoint_coupling` distinct q-points
      !>    are returned.
      !> For the purposes of `max_qpoint_coupling`, a q-point and its pair are
      !>    counted as a single q-point.
      integer,                 intent(in), optional   :: max_qpoint_coupling
      !> If `conserve_momentum` is `.true.` then only q-point combinations
      !>    which conserve momentum (i.e. sum q=G) are included in the stars.
      !> `conserve_momentum` defaults to `.false.`.
      logical,                 intent(in), optional :: conserve_momentum
      !> `qpoints` is only required if `conserve_momentum` is `.true.`.
      type(QpointData),        intent(in), optional :: qpoints(:)
      type(CombinationQpointStar), allocatable      :: output(:)
    end function
  end interface
end module
