!> Provides the [[QpointStarProduct(type)]] class, and related methods.
module caesar_qpoint_star_product_module
  use caesar_common_module
  
  use caesar_subspaces_module
  
  use caesar_qpoint_star_module
  use caesar_subspace_qpoint_stars_module
  implicit none
  
  private
  
  public :: QpointStarProduct
  public :: generate_qpoint_star_products
  
  !> A product of [[QpointStar(type)]]s corresponding to a particular
  !>    [[SubspaceCombination]].
  type, extends(Stringsable) :: QpointStarProduct
    type(SubspaceCombination)     :: subspace_combination
    !> `stars(i)` corresponds to the `i`th subspace in the combination.
    type(QpointStar), allocatable :: qpoint_stars(:)
  contains
    procedure, public :: read  => read_QpointStarProduct
    procedure, public :: write => write_QpointStarProduct
  end type
  
  interface QpointStarProduct
    !> Constructor for [[QpointStarProduct(type)]] objects.
    module function new_QpointStarProduct(subspace_combination, &
       & qpoint_stars) result(this)
      type(SubspaceCombination), intent(in) :: subspace_combination
      type(QpointStar),          intent(in) :: qpoint_stars(:)
      type(QpointStarProduct)               :: this
    end function
  end interface
  
  interface
    !> Convert a [[String(type)]] array to a [[QpointStarProduct(type)]].
    module subroutine read_QpointStarProduct(this,input)
      class(QpointStarProduct), intent(out) :: this
      type(String),             intent(in)  :: input(:)
    end subroutine
    
    !> Convert a [[QpointStarProduct(type)]] to a [[String(type)]] array.
    module function write_QpointStarProduct(this) result(output)
      class(QpointStarProduct), intent(in) :: this
      type(String), allocatable            :: output(:)
    end function
  end interface
  
  interface QpointStarProduct
    !> Convert a [[String(type)]] array to a [[QpointStarProduct(type)]].
    module function new_QpointStarProduct_Strings(input) &
       & result(this)
      type(String), intent(in) :: input(:)
      type(QpointStarProduct)  :: this
    end function
    
    !> Convert a [[StringArray(type)]] to a [[QpointStarProduct(type)]].
    impure elemental module function new_QpointStarProduct_StringArray( &
       & input) result(this)
      type(StringArray), intent(in) :: input
      type(QpointStarProduct)       :: this
    end function
  end interface
  
  interface
    !> Generates the set of [[QpointStarProduct(type)]]s for a given
    !>    `subspace_combination`, from a set of `subspace_stars`.
    module function generate_qpoint_star_products(subspace_combination, &
       & subspace_qpoint_stars) result(output)
      type(SubspaceCombination), intent(in) :: subspace_combination
      type(SubspaceQpointStars), intent(in) :: subspace_qpoint_stars(:)
      type(QpointStarProduct), allocatable  :: output(:)
    end function
  end interface
end module
