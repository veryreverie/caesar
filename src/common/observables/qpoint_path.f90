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
  
  interface QpointPathVertex
    module procedure new_QpointPathVertex
    module procedure new_QpointPathVertex_Strings
    module procedure new_QpointPathVertex_StringArray
  end interface
  
  ! The edges between the high-symmetry points.
  type, extends(NoDefaultConstructor) :: QpointPathEdge
    integer  :: start_vertex_id
    integer  :: end_vertex_id
    real(dp) :: path_fraction
    integer  :: no_points
  end type
  
  interface QpointPathEdge
    module procedure new_QpointPathEdge
  end interface
  
  ! The entire path.
  type, extends(NoDefaultConstructor) :: QpointPath
    type(QpointPathVertex), allocatable :: vertices(:)
    type(QpointPathEdge),   allocatable :: edges(:)
  contains
    procedure, public :: path_qpoints => path_qpoints_QpointPath
  end type
  
  interface QpointPath
    module procedure new_QpointPath
    module procedure new_QpointPath_string
  end interface
  
  ! A return type for a q-point and its fraction along the path.
  type, extends(NoDefaultConstructor) :: PathQpointAndFraction
    type(RealVector) :: qpoint
    real(dp)         :: path_fraction
  end type
  
  interface PathQpointAndFraction
    module procedure new_PathQpointAndFraction
  end interface
contains

! Element constructors.
impure elemental function new_QpointPathVertex(label,qpoint,path_fraction, &
   & point_index) result(this)
  implicit none
  
  type(String),     intent(in) :: label
  type(RealVector), intent(in) :: qpoint
  real(dp),         intent(in) :: path_fraction
  integer,          intent(in) :: point_index
  type(QpointPathVertex)       :: this
  
  this%label = label
  this%qpoint = qpoint
  this%path_fraction = path_fraction
  this%point_index   = point_index
end function

impure elemental function new_QpointPathEdge(start_vertex_id,end_vertex_id, &
   & path_fraction,no_points) result(this)
  implicit none
  
  integer,  intent(in) :: start_vertex_id
  integer,  intent(in) :: end_vertex_id
  real(dp), intent(in) :: path_fraction
  integer,  intent(in) :: no_points
  type(QpointPathEdge) :: this
  
  this%start_vertex_id = start_vertex_id
  this%end_vertex_id = end_vertex_id
  this%path_fraction = path_fraction
  this%no_points = no_points
end function

function new_QpointPath(vertices,edges) result(this)
  implicit none
  
  type(QpointPathVertex), intent(in) :: vertices(:)
  type(QpointPathEdge),   intent(in) :: edges(:)
  type(QpointPath)                   :: this
  
  this%vertices = vertices
  this%edges    = edges
end function

impure elemental function new_PathQpointAndFraction(qpoint,path_fraction) &
   & result(this)
  implicit none
  
  type(RealVector), intent(in) :: qpoint
  real(dp),         intent(in) :: path_fraction
  type(PathQpointAndFraction)  :: this
  
  this%qpoint = qpoint
  this%path_fraction = path_fraction
end function

! Constructor from input string.
! N.B. the total number of points along the path will only approximately equal
!    total_path_points, due to rounding errors.
function new_QpointPath_string(path,total_path_points) result(this)
  implicit none
  
  type(String), intent(in) :: path
  integer,      intent(in) :: total_path_points
  type(QpointPath)         :: this
  
  ! Path segments.
  type(String),     allocatable :: segments(:)
  type(String),     allocatable :: segment_vertices(:)
  
  ! Initial edge data.
  integer                       :: no_edges
  integer,          allocatable :: edge_start_vertices(:)
  integer,          allocatable :: edge_end_vertices(:)
  
  ! Initial vertex data.
  integer                       :: no_vertices
  type(String),     allocatable :: vertex_labels(:)
  type(RealVector), allocatable :: vertex_qpoints(:)
  
  ! Path lengths.
  real(dp), allocatable :: edge_lengths(:)
  real(dp), allocatable :: fractional_edge_lengths(:)
  integer,  allocatable :: edge_no_points(:)
  real(dp), allocatable :: vertex_distances(:)
  integer,  allocatable :: point_indices(:)
  
  ! Final forms for vertex and edge data.
  type(QpointPathVertex), allocatable :: vertices(:)
  type(QpointPathEdge),   allocatable :: edges(:)
  
  integer :: i,ialloc
  
  ! Split the path into segments.
  segments = tokens(path, delimiters=['|'])
  
  ! Parse vertices and edges, segment by segment.
  no_edges = 0
  no_vertices = 0
  allocate( edge_start_vertices(0), &
          & edge_end_vertices(0),   &
          & vertex_labels(0),       &
          & vertex_qpoints(0),      &
          & stat=ialloc); call err(ialloc)
  do i=1,size(segments)
    segment_vertices = tokens(segments(i), delimiter=',')
    if (size(segment_vertices)<2) then
      call print_line(ERROR//': There must be at least two high-symmetry &
         &points along each path segment.')
      call err()
    endif
    
    ! Add edges, initially defined by the ids of the vertices at either end
    !    of each edge.
    no_edges = no_edges + size(segment_vertices)-1
    edge_start_vertices = [ edge_start_vertices,           &
                          & [( size(vertex_labels)+i,      &
                          &    i=1,                        &
                          &    size(segment_vertices)-1 )] ]
    edge_end_vertices = [ edge_end_vertices,           &
                        & [( size(vertex_labels)+i,    &
                        &    i=2,                      &
                        &    size(segment_vertices) )] ]
    
    ! Add vertices.
    no_vertices = no_vertices + size(segment_vertices)
    vertex_labels  = [ vertex_labels,                      &
                     & [( token(segment_vertices(i),1),    &
                     &    i=1,                             &
                     &    size(segment_vertices)        )] ]
    vertex_qpoints = [ vertex_qpoints,                                   &
                     & [( vec(dble(tokens(segment_vertices(i),2,4))),    &
                     &    i=1,                                           &
                     &    size(segment_vertices)                      )] ]
  enddo
  
  ! Calculate edge lengths.
  edge_lengths = l2_norm(   vertex_qpoints(edge_end_vertices)   &
                        & - vertex_qpoints(edge_start_vertices) )
  fractional_edge_lengths = edge_lengths / sum(edge_lengths)
  edge_no_points = [nint(total_path_points*fractional_edge_lengths)]
  
  if (any(edge_no_points<2)) then
    call print_line(ERROR//': A dispersion curve edge has less than two &
       &points along its length. Please increase no_path_points.')
    call err()
  endif
  
  ! Calculate the cumulative distance to each vertex.
  vertex_distances = [(0.0_dp, i=1, size(vertex_labels))]
  point_indices = [(i, i=1, size(vertex_labels))]
  do i=1,no_edges
    vertex_distances(edge_end_vertices(i):) =      &
       &   vertex_distances(edge_end_vertices(i):) &
       & + fractional_edge_lengths(i)
    point_indices(edge_end_vertices(i):) =      &
       &   point_indices(edge_end_vertices(i):) &
       & + edge_no_points(i)-1
  enddo
  
  ! Construct output.
  vertices = QpointPathVertex( vertex_labels,    &
                             & vertex_qpoints,   &
                             & vertex_distances, &
                             & point_indices     )
  edges = QpointPathEdge( start_vertex_id = edge_start_vertices,     &
                        & end_vertex_id   = edge_end_vertices,       &
                        & path_fraction   = fractional_edge_lengths, &
                        & no_points       = edge_no_points           )
  this = QpointPath(vertices, edges)
end function

! Generate all of the q-points along the path.
impure elemental function path_qpoint_QpointPathVertex(this) result(output)
  implicit none
  
  class(QpointPathVertex), intent(in) :: this
  type(PathQpointAndFraction)         :: output
  
  output = PathQpointAndFraction(this%qpoint, this%path_fraction)
end function

function path_qpoints_QpointPath(this) result(output)
  implicit none
  
  class(QpointPath), intent(in)            :: this
  type(PathQpointAndFraction), allocatable :: output(:)
  
  type(PathQpointAndFraction) :: start_qpoint
  type(PathQpointAndFraction) :: end_qpoint
  
  integer :: no_points
  
  integer :: i,j,k,ialloc
  
  no_points = this%vertices(size(this%vertices))%point_index
  allocate(output(no_points), stat=ialloc); call err(ialloc)
  k = 0
  do i=1,size(this%edges)
    start_qpoint = this%vertices(this%edges(i)%start_vertex_id)%path_qpoint()
    end_qpoint = this%vertices(this%edges(i)%end_vertex_id)%path_qpoint()
    
    ! Add in the q-point at the start of the edge, if it has not been added as
    !    the q-point at the end of the last edge.
    if (i==1) then
      k = k+1
      output(k) = start_qpoint
    elseif (this%edges(i)%start_vertex_id/=this%edges(i-1)%end_vertex_id) then
      k = k+1
      output(k) = start_qpoint
    endif
    
    ! Add in the q-points along the edge.
    do j=1,this%edges(i)%no_points-1
      k = k+1
      output(k) = PathQpointAndFraction(                   &
         & qpoint        = ( (this%edges(i)%no_points-j)   &
         &                   * start_qpoint%qpoint         &
         &                   + j                           &
         &                   * end_qpoint%qpoint         ) &
         &                 / this%edges(i)%no_points,      &
         & path_fraction = ( (this%edges(i)%no_points-j)   &
         &                   * start_qpoint%path_fraction  &
         &                   + j                           &
         &                   * end_qpoint%path_fraction  ) &
         &                 / this%edges(i)%no_points       )
    enddo
    
    ! Add in the q-point at the end of the edge.
    k = k+1
    output(k) = end_qpoint
  enddo
  
  if (k/=size(output)) then
    call print_line(CODE_ERROR//': Unable to construct q-point path.')
    call err()
  endif
end function

! I/O.
subroutine read_QpointPathVertex(this,input)
  implicit none
  
  class(QpointPathVertex), intent(out) :: this
  type(String),            intent(in)  :: input(:)
  
  type(String)     :: label
  type(RealVector) :: qpoint
  real(dp)         :: path_fraction
  integer          :: point_index
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(QpointPathVertex)
    label = token(line(1),4)
    qpoint = vec(dble(tokens(line(2),3,5)))
    path_fraction = dble(token(line(3),5))
    point_index = int(token(line(4),5))
    this = QpointPathVertex(label,qpoint,path_fraction,point_index)
  class default
    call err()
  end select
end subroutine

function write_QpointPathVertex(this) result(output)
  implicit none
  
  class(QpointPathVertex), intent(in) :: this
  type(String), allocatable           :: output(:)
  
  select type(this); type is(QpointPathVertex)
    output = [ 'q-point label       : '//this%label,         &
             & 'q-point             : '//this%qpoint,        &
             & 'Fraction along path : '//this%path_fraction, &
             & 'Index along path    : '//this%point_index    ]
  class default
    call err()
  end select
end function

function new_QpointPathVertex_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(QpointPathVertex)   :: this
  
  call this%read(input)
end function

impure elemental function new_QpointPathVertex_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(QpointPathVertex)        :: this
  
  this = QpointPathVertex(str(input))
end function
end module
