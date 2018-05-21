! ======================================================================
! The R-vectors at which the displacements in a particular subspace must be
!    sampled in order to construct the VSCF Hamiltonian.
! ======================================================================
module vscf_rvectors_module
  use common_module
  use sampling_points_module
  implicit none
  
  private
  
  public :: VscfRvectors
  public :: size
  public :: construct_vscf_rvectors
  
  type, extends(Stringsable) :: VscfRvectors
    integer                      :: subspace_id
    type(IntVector), allocatable :: rvectors(:)
  contains
    procedure, public :: read  => read_VscfRvectors
    procedure, public :: write => write_VscfRvectors
  end type
  
  interface VscfRvectors
    module procedure new_VscfRvectors
  end interface
  
  interface size
    module procedure size_VscfRvectors
  end interface
contains

! ----------------------------------------------------------------------
! Basic functionality: Constructor and size() function.
! ----------------------------------------------------------------------
function new_VscfRvectors(subspace_id,rvectors) result(this)
  implicit none
  
  integer,         intent(in) :: subspace_id
  type(IntVector), intent(in) :: rvectors(:)
  type(VscfRvectors)          :: this
  
  this%subspace_id = subspace_id
  this%rvectors    = rvectors
end function

function size_VscfRvectors(this) result(output)
  implicit none
  
  type(VscfRvectors), intent(in) :: this
  integer                        :: output
  
  output = size(this%rvectors)
end function

! ----------------------------------------------------------------------
! Construct a set of R-vectors for a given sampling point.
! ----------------------------------------------------------------------
function construct_vscf_rvectors(sampling_point,supercell,real_modes,qpoints) &
   & result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in) :: sampling_point
  type(StructureData),        intent(in) :: supercell
  type(RealMode),             intent(in) :: real_modes(:)
  type(QpointData),           intent(in) :: qpoints(:)
  type(VscfRvectors), allocatable        :: output(:)
  
  ! Modes and q-points.
  type(RealMode),   allocatable :: modes(:)
  type(QpointData), allocatable :: mode_qpoints(:)
  integer,          allocatable :: subspace_ids(:)
  
  ! Modes and q-points within a given subspace.
  integer,          allocatable :: modes_in_subspace(:)
  type(RealMode),   allocatable :: subspace_modes(:)
  type(QpointData), allocatable :: subspace_qpoints(:)
  
  ! Variables for calculating R-vectors for a given subspace.
  integer                        :: no_rvectors
  type(IntVector),   allocatable :: rvectors(:)
  type(IntFraction), allocatable :: q_dot_r(:,:)
  
  ! Variables for removing global R-vector shifts.
  integer, allocatable :: rvector_counts(:)
  
  integer :: i,j,k,ialloc
  
  ! List the modes at which the sampling point has non-zero displacement,
  !    and list the corresponding subspace ids.
  modes = sampling_point%modes(real_modes)
  subspace_ids = modes%subspace_id
  
  ! List q-points corresponding to said modes.
  allocate(mode_qpoints(size(modes)), stat=ialloc); call err(ialloc)
  do i=1,size(modes)
    mode_qpoints(i) = qpoints(first(qpoints%id==modes(i)%qpoint_id))
  enddo
  
  ! De-duplicate the subspace ids.
  subspace_ids = subspace_ids(set(subspace_ids))
  
  ! Identify the R-vectors in the supercell which lead to different
  !    displacements within each subspace.
  allocate(output(size(subspace_ids)), stat=ialloc); call err(ialloc)
  do i=1,size(subspace_ids)
    modes_in_subspace = filter(modes%subspace_id==subspace_ids(i))
    subspace_modes = modes(modes_in_subspace)
    subspace_qpoints = mode_qpoints(modes_in_subspace)
    
    allocate( rvectors(size(supercell%rvectors)),                       &
            & q_dot_r(size(subspace_qpoints),size(supercell%rvectors)), &
            & stat=ialloc); call err(ialloc)
    no_rvectors = 0
    do_j : do j=1,size(supercell%rvectors)
      ! Copy over each R-vector from the supercell,
      !   and calculate q.R over all subspace q-points.
      no_rvectors = no_rvectors+1
      rvectors(no_rvectors) = supercell%rvectors(j)
      do k=1,size(subspace_qpoints)
        q_dot_r(k,no_rvectors) = subspace_qpoints(k)%qpoint &
                             & * rvectors(no_rvectors)
      enddo
      
      ! Ignore the R-vector if every q.R is the same as a previous R-vector.
      do k=1,no_rvectors-1
        if (all(q_dot_r(:,k)==q_dot_r(:,no_rvectors))) then
          no_rvectors = no_rvectors-1
          cycle do_j
        endif
      enddo
    enddo do_j
    
    output(i) = VscfRvectors(subspace_ids(i),rvectors(:no_rvectors))
    deallocate( rvectors, &
              & q_dot_r,  &
              & stat=ialloc); call err(ialloc)
  enddo
  
  ! If every subspace is shifted by the same R-vector,
  !    then nothing has changed.
  ! To avoid this over-sampling, one subspace can be sampled at just one
  !    R-vector (which w.l.g. can be [0,0,0]).
  allocate(rvector_counts(size(output)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    rvector_counts(i) = size(output(i))
  enddo
  output(minloc(rvector_counts,1))%rvectors = [zeroes(3)]
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_VscfRvectors(this,input)
  implicit none
  
  class(VscfRvectors), intent(out) :: this
  type(String),        intent(in)  :: input(:)
  
  integer                      :: subspace_id
  type(IntVector), allocatable :: rvectors(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(VscfRvectors)
    ! Read in subspace id.
    line = split_line(input(1))
    subspace_id = int(line(4))
    
    ! Read in R-vectors.
    allocate(rvectors(size(input)-2), stat=ialloc); call err(ialloc)
    do i=1,size(rvectors)
      rvectors(i) = input(i+2)
    enddo
    
    this = VscfRvectors(subspace_id,rvectors)
  end select
end subroutine

function write_VscfRvectors(this) result(output)
  implicit none
  
  class(VscfRvectors), intent(in) :: this
  type(String), allocatable       :: output(:)
  
  select type(this); type is(VscfRvectors)
    output = [ 'Subspace ID : '//this%subspace_id, &
             & str('R-vectors :'),                 &
             & str(this%rvectors)                             ]
  end select
end function
end module
