! ======================================================================
! A path through reciprocal space, for the generation of dispersion curves.
! ======================================================================
module caesar_qpoint_path_module
  use caesar_utils_module
  implicit none
  
  private
  
  public :: QpointPathVertex
  public :: QpointPathEdge
  public :: QpointPath
  public :: PathQpointAndFraction
  
  ! The high-symmetry points along the q-point path.
  type, extends(Stringsable) :: QpointPathVertex
    type(String)     :: label
    type(RealVector) :: qpoint
    real(dp)         :: path_fraction
    integer          :: point_index
  contains
    procedure, public :: path_qpoint => path_qpoint_QpointPathVertex
    ! I/O.
    procedure, public :: read  => read_QpointPathVertex
    procedure, public :: write => write_QpointPathVertex
  end type
  
  ! The edges between the high-symmetry points.
  type, extends(NoDefaultConstructor) :: QpointPathEdge
    integer  :: start_vertex_id
    integer  :: end_vertex_id
    real(dp) :: path_fraction
    integer  :: no_points
  end type
  
  ! The entire path.
  type, extends(NoDefaultConstructor) :: QpointPath
    type(QpointPathVertex), allocatable :: vertices(:)
    type(QpointPathEdge),   allocatable :: edges(:)
  contains
    procedure, public :: path_qpoints => path_qpoints_QpointPath
  end type
  
  ! A return type for a q-point and its fraction along the path.
  type, extends(NoDefaultConstructor) :: PathQpointAndFraction
    type(RealVector) :: qpoint
    real(dp)         :: path_fraction
  end type
  
  interface QpointPathVertex
    ! Element constructors.
    impure elemental module function new_QpointPathVertex(label,qpoint, &
       & path_fraction,point_index) result(this) 
      type(String),     intent(in) :: label
      type(RealVector), intent(in) :: qpoint
      real(dp),         intent(in) :: path_fraction
      integer,          intent(in) :: point_index
      type(QpointPathVertex)       :: this
    end function
  end interface
  
  interface QpointPathEdge
    impure elemental module function new_QpointPathEdge(start_vertex_id, &
       & end_vertex_id,path_fraction,no_points) result(this) 
      integer,  intent(in) :: start_vertex_id
      integer,  intent(in) :: end_vertex_id
      real(dp), intent(in) :: path_fraction
      integer,  intent(in) :: no_points
      type(QpointPathEdge) :: this
    end function
  end interface
  
  interface QpointPath
    module function new_QpointPath(vertices,edges) result(this) 
      type(QpointPathVertex), intent(in) :: vertices(:)
      type(QpointPathEdge),   intent(in) :: edges(:)
      type(QpointPath)                   :: this
    end function
  end interface
  
  interface PathQpointAndFraction
    impure elemental module function new_PathQpointAndFraction(qpoint, &
       & path_fraction) result(this) 
      type(RealVector), intent(in) :: qpoint
      real(dp),         intent(in) :: path_fraction
      type(PathQpointAndFraction)  :: this
    end function
  end interface
  
  interface QpointPath
    ! Constructor from input string.
    ! N.B. the total number of points along the path will only approximately equal
    !    total_path_points, due to rounding errors.
    module function new_QpointPath_string(path,total_path_points) result(this) 
      type(String), intent(in) :: path
      integer,      intent(in) :: total_path_points
      type(QpointPath)         :: this
    end function
  end interface
  
  interface
    ! Generate all of the q-points along the path.
    impure elemental module function path_qpoint_QpointPathVertex(this) &
       & result(output) 
      class(QpointPathVertex), intent(in) :: this
      type(PathQpointAndFraction)         :: output
    end function
  end interface
  
  interface
    module function path_qpoints_QpointPath(this) result(output) 
      class(QpointPath), intent(in)            :: this
      type(PathQpointAndFraction), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! I/O.
    module subroutine read_QpointPathVertex(this,input) 
      class(QpointPathVertex), intent(out) :: this
      type(String),            intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_QpointPathVertex(this) result(output) 
      class(QpointPathVertex), intent(in) :: this
      type(String), allocatable           :: output(:)
    end function
  end interface
  
  interface QpointPathVertex
    module function new_QpointPathVertex_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(QpointPathVertex)   :: this
    end function
  
    impure elemental module function new_QpointPathVertex_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(QpointPathVertex)        :: this
    end function
  end interface
end module
