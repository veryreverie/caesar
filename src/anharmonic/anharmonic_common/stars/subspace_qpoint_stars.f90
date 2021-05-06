!> Provides the [[SubspaceQpointStars(type)]] class, and related methods.
module caesar_subspace_qpoint_stars_module
  use caesar_common_module
  
  use caesar_qpoint_star_module
  implicit none
  
  private
  
  public :: SubspaceQpointStars
  public :: generate_subspace_qpoint_stars
  
  !> The set of [[QpointStars(type)]]s corresponding to a given
  !>    [[DegenerateSubspace(type)]].
  type, extends(Stringsable) :: SubspaceQpointStars
    integer                        :: subspace_id
    type(QpointStars), allocatable :: powers(:)
  contains
    procedure, public :: read  => read_SubspaceQpointStars
    procedure, public :: write => write_SubspaceQpointStars
  end type
  
  interface SubspaceQpointStars
    !> Constructor for [[SubspaceQpointStars(type)]] objects.
    module function new_SubspaceQpointStars(subspace_id,powers) result(this)
      integer,           intent(in) :: subspace_id
      type(QpointStars), intent(in) :: powers(:)
      type(SubspaceQpointStars)     :: this
    end function
  end interface
  
  interface
    !> Convert a [[String(type)]] array to a [[SubspaceQpointStars(type)]].
    module subroutine read_SubspaceQpointStars(this,input)
      class(SubspaceQpointStars), intent(out) :: this
      type(String),               intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    !> Convert a [[SubspaceQpointStars(type)]] to a [[String(type)]] array.
    module function write_SubspaceQpointStars(this) result(output)
      class(SubspaceQpointStars), intent(in) :: this
      type(String), allocatable              :: output(:)
    end function
  end interface
  
  interface SubspaceQpointStars
    !> Convert a [[String(type)]] array to a [[SubspaceQpointStars(type)]].
    module function new_SubspaceQpointStars_Strings(input) &
       & result(this)
      type(String), intent(in)  :: input(:)
      type(SubspaceQpointStars) :: this
    end function
    
    !> Convert a [[StringArray(type)]] to a [[SubspaceQpointStars(type)]].
    impure elemental module function new_SubspaceQpointStars_StringArray( &
       & input) result(this)
      type(StringArray), intent(in) :: input
      type(SubspaceQpointStars)     :: this
    end function
  end interface
  
  interface
    !> Generates the q-point stars corresponding to a set of `subspaces`.
    module function generate_subspace_qpoint_stars(subspaces,modes,qpoints, &
       & symmetries,max_power,conserve_momentum) result(output)
      type(DegenerateSubspace), intent(in)           :: subspaces(:)
      !> A list of complex modes,
      !>    which must include all modes in the subspaces.
      type(ComplexMode),        intent(in)           :: modes(:)
      !> A list of q-points,
      !>    which must include all q-points corresponding to the modes.
      type(QpointData),         intent(in)           :: qpoints(:)
      !> The symmetries of the primitive cell.
      type(SymmetryOperator),   intent(in)           :: symmetries(:)
      !> The q-point stars will be generated with total powers <= `max_power`.
      integer,                  intent(in)           :: max_power
      !> If true, only q-point stars which conserve momentum (i.e. sum q = G)
      !>    will be returned. Defaults to `false`.
      logical,                  intent(in), optional :: conserve_momentum
      type(SubspaceQpointStars), allocatable         :: output(:)
    end function
  end interface
end module
