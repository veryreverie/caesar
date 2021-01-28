! ======================================================================
! As caesar_vscf_rvector_module, but produces a list of R-vector combinations
!    rather than the R-vectors for each subspace individually.
! ======================================================================
module caesar_vscf_rvectors_module
  use caesar_common_module
  
  use caesar_sampling_points_module
  use caesar_vscf_rvector_module
  implicit none
  
  private
  
  public :: VscfRvectors
  public :: size
  public :: construct_vscf_rvectors
  
  type, extends(Stringsable) :: VscfRvectors
    type(VscfRvector), allocatable :: vscf_rvectors(:)
  contains
    ! List R-vectors for given modes.
    procedure, public :: rvectors => rvectors_VscfRvectors
    
    ! Transform vectors in normal-mode co-ordinates.
    generic,   public  :: transform =>                       &
                        & transform_ComplexModeDisplacement, &
                        & transform_ComplexModeForce,        &
                        & transform_RealModeDisplacement,    &
                        & transform_RealModeForce
    procedure, private :: transform_ComplexModeDisplacement
    procedure, private :: transform_ComplexModeForce
    procedure, private :: transform_RealModeDisplacement
    procedure, private :: transform_RealModeForce
    
    ! The inverse operation to transform.
    generic,   public  :: inverse_transform =>                       &
                        & inverse_transform_ComplexModeDisplacement, &
                        & inverse_transform_ComplexModeForce,        &
                        & inverse_transform_RealModeDisplacement,    &
                        & inverse_transform_RealModeForce
    procedure, private :: inverse_transform_ComplexModeDisplacement
    procedure, private :: inverse_transform_ComplexModeForce
    procedure, private :: inverse_transform_RealModeDisplacement
    procedure, private :: inverse_transform_RealModeForce
    
    ! Helper functions for transform and inverse_transform.
    procedure, private :: transform_complex_magnitudes
    procedure, private :: transform_real_magnitudes
    
    ! I/O.
    procedure, public :: read  => read_VscfRvectors
    procedure, public :: write => write_VscfRvectors
  end type
  
  interface VscfRvectors
    module procedure new_VscfRvectors
    module procedure new_VscfRvectors_Strings
    module procedure new_VscfRvectors_StringArray
  end interface
  
  interface size
    module procedure size_VscfRvectors
  end interface
  
  interface operator(//)
    module procedure concatenate_VscfRvectors_VscfRvector
  end interface
  
  ! Helper type.
  type :: RvectorArray
    integer                      :: subspace_id
    type(IntVector), allocatable :: rvectors(:)
  end type
  
  interface RvectorArray
    module procedure new_RvectorArray
  end interface
  
  interface size
    module procedure size_RvectorArray
  end interface
contains

! ----------------------------------------------------------------------
! Constructor and size() function for VscfRvectors and RvectorArray.
! ----------------------------------------------------------------------
function new_VscfRvectors(vscf_rvectors) result(this)
  implicit none
  
  type(VscfRvector), intent(in), optional :: vscf_rvectors(:)
  type(VscfRvectors)                      :: this
  
  integer :: ialloc
  
  if (present(vscf_rvectors)) then
    this%vscf_rvectors = vscf_rvectors
  else
    allocate(this%vscf_rvectors(0), stat=ialloc); call err(ialloc)
  endif
end function

function size_VscfRvectors(this) result(output)
  implicit none
  
  type(VscfRvectors), intent(in) :: this
  integer                        :: output
  
  output = size(this%vscf_rvectors)
end function

function new_RvectorArray(subspace_id,rvectors) result(this)
  implicit none
  
  integer,         intent(in) :: subspace_id
  type(IntVector), intent(in) :: rvectors(:)
  type(RvectorArray)          :: this
  
  this%subspace_id = subspace_id
  this%rvectors    = rvectors
end function

function size_RvectorArray(this) result(output)
  implicit none
  
  type(RvectorArray), intent(in) :: this
  integer                        :: output
  
  output = size(this%rvectors)
end function

! ----------------------------------------------------------------------
! Concatenation of a VscfRvector to a VscfRvectors.
! ----------------------------------------------------------------------
function concatenate_VscfRvectors_VscfRvector(this,that) result(output)
  implicit none
  
  type(VscfRvectors), intent(in) :: this
  type(VscfRvector),  intent(in) :: that
  type(VscfRvectors)             :: output
  
  output = VscfRvectors([this%vscf_rvectors, that])
end function

! ----------------------------------------------------------------------
! List R-vectors for the given set of modes.
! ----------------------------------------------------------------------
function rvectors_VscfRvectors(this,real_modes) result(output)
  implicit none
  
  class(VscfRvectors), intent(in) :: this
  type(RealMode),      intent(in) :: real_modes(:)
  type(IntVector), allocatable    :: output(:)
  
  integer :: i,j,ialloc
  
  allocate(output(size(real_modes)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    j = first( this%vscf_rvectors%subspace_id==real_modes(i)%subspace_id, &
             & default=0)
    if (j==0) then
      output(i) = zeroes(3)
    else
      output(i) = this%vscf_rvectors(j)%rvector
    endif
  enddo
end function

! ----------------------------------------------------------------------
! Transform a displacement.
! ----------------------------------------------------------------------
function transform_ComplexModeDisplacement(this,displacement,modes,qpoints) &
   & result(output)
  implicit none
  
  class(VscfRvectors),           intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexMode),             intent(in) :: modes(:)
  type(QpointData),              intent(in) :: qpoints(:)
  type(ComplexModeDisplacement)             :: output
  
  type(ComplexMode), allocatable :: output_modes(:)
  complex(dp),       allocatable :: magnitudes(:)
  
  integer :: i
  
  output_modes = select_modes(displacement%vectors, modes)
  magnitudes = displacement%vectors%magnitude
  
  do i=1,size(output_modes)
    if (.not. any(output_modes%id==output_modes(i)%paired_id)) then
      output_modes = [ output_modes,                                     &
                     & modes(first(modes%id==output_modes(i)%paired_id)) ]
      magnitudes = [magnitudes, cmplx(0.0_dp,0.0_dp,dp)]
    endif
  enddo
  
  output = ComplexModeDisplacement(                          &
     & output_modes,                                         &
     & this%transform_complex_magnitudes( output_modes,      &
     &                                    magnitudes,        &
     &                                    qpoints,           &
     &                                    inverse = .false. ))
end function

function transform_ComplexModeForce(this,force,modes,qpoints) &
   & result(output)
  implicit none
  
  class(VscfRvectors),    intent(in) :: this
  type(ComplexModeForce), intent(in) :: force
  type(ComplexMode),      intent(in) :: modes(:)
  type(QpointData),       intent(in) :: qpoints(:)
  type(ComplexModeForce)             :: output
  
  type(ComplexMode), allocatable :: output_modes(:)
  complex(dp),       allocatable :: magnitudes(:)
  
  integer :: i
  
  output_modes = select_modes(force%vectors, modes)
  magnitudes = force%vectors%magnitude
  
  do i=1,size(output_modes)
    if (.not. any(output_modes%id==output_modes(i)%paired_id)) then
      output_modes = [ output_modes,                                     &
                     & modes(first(modes%id==output_modes(i)%paired_id)) ]
      magnitudes = [magnitudes, cmplx(0.0_dp,0.0_dp,dp)]
    endif
  enddo
  
  output = ComplexModeForce(                                 &
     & output_modes,                                         &
     & this%transform_complex_magnitudes( output_modes,      &
     &                                    magnitudes,        &
     &                                    qpoints,           &
     &                                    inverse = .false. ))
end function

function transform_RealModeDisplacement(this,displacement,modes,qpoints) &
   & result(output)
  implicit none
  
  class(VscfRvectors),        intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  type(RealMode),             intent(in) :: modes(:)
  type(QpointData),           intent(in) :: qpoints(:)
  type(RealModeDisplacement)             :: output
  
  type(RealMode), allocatable :: output_modes(:)
  real(dp),       allocatable :: magnitudes(:)
  
  integer :: i
  
  output_modes = select_modes(displacement%vectors, modes)
  magnitudes = displacement%vectors%magnitude
  
  do i=1,size(output_modes)
    if (.not. any(output_modes%id==output_modes(i)%paired_id)) then
      output_modes = [ output_modes,                                     &
                     & modes(first(modes%id==output_modes(i)%paired_id)) ]
      magnitudes = [magnitudes, 0.0_dp]
    endif
  enddo
  
  output = RealModeDisplacement(                          &
     & output_modes,                                      &
     & this%transform_real_magnitudes( output_modes,      &
     &                                 magnitudes,        &
     &                                 qpoints,           &
     &                                 inverse = .false. ))
end function

function transform_RealModeForce(this,force,modes,qpoints) &
   & result(output)
  implicit none
  
  class(VscfRvectors), intent(in) :: this
  type(RealModeForce), intent(in) :: force
  type(RealMode),      intent(in) :: modes(:)
  type(QpointData),    intent(in) :: qpoints(:)
  type(RealModeForce)             :: output
  
  type(RealMode), allocatable :: output_modes(:)
  real(dp),       allocatable :: magnitudes(:)
  
  integer :: i
  
  output_modes = select_modes(force%vectors, modes)
  magnitudes = force%vectors%magnitude
  
  do i=1,size(output_modes)
    if (.not. any(output_modes%id==output_modes(i)%paired_id)) then
      output_modes = [ output_modes,                                     &
                     & modes(first(modes%id==output_modes(i)%paired_id)) ]
      magnitudes = [magnitudes, 0.0_dp]
    endif
  enddo
  
  output = RealModeForce(                                 &
     & output_modes,                                      &
     & this%transform_real_magnitudes( output_modes,      &
     &                                 magnitudes,        &
     &                                 qpoints,           &
     &                                 inverse = .false. ))
end function

function inverse_transform_ComplexModeDisplacement(this,displacement,modes, &
   & qpoints)  result(output)
  implicit none
  
  class(VscfRvectors),           intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexMode),             intent(in) :: modes(:)
  type(QpointData),              intent(in) :: qpoints(:)
  type(ComplexModeDisplacement)             :: output
  
  type(ComplexMode), allocatable :: output_modes(:)
  complex(dp),       allocatable :: magnitudes(:)
  
  integer :: i
  
  output_modes = select_modes(displacement%vectors, modes)
  magnitudes = displacement%vectors%magnitude
  
  do i=1,size(output_modes)
    if (.not. any(output_modes%id==output_modes(i)%paired_id)) then
      output_modes = [ output_modes,                                     &
                     & modes(first(modes%id==output_modes(i)%paired_id)) ]
      magnitudes = [magnitudes, cmplx(0.0_dp,0.0_dp,dp)]
    endif
  enddo
  
  output = ComplexModeDisplacement(                         &
     & output_modes,                                        &
     & this%transform_complex_magnitudes( output_modes,     &
     &                                    magnitudes,       &
     &                                    qpoints,          &
     &                                    inverse = .true. ))
end function

function inverse_transform_ComplexModeForce(this,force,modes,qpoints) &
   & result(output)
  implicit none
  
  class(VscfRvectors),    intent(in) :: this
  type(ComplexModeForce), intent(in) :: force
  type(ComplexMode),      intent(in) :: modes(:)
  type(QpointData),       intent(in) :: qpoints(:)
  type(ComplexModeForce)             :: output
  
  type(ComplexMode), allocatable :: output_modes(:)
  complex(dp),       allocatable :: magnitudes(:)
  
  integer :: i
  
  output_modes = select_modes(force%vectors, modes)
  magnitudes = force%vectors%magnitude
  
  do i=1,size(output_modes)
    if (.not. any(output_modes%id==output_modes(i)%paired_id)) then
      output_modes = [ output_modes,                                     &
                     & modes(first(modes%id==output_modes(i)%paired_id)) ]
      magnitudes = [magnitudes, cmplx(0.0_dp,0.0_dp,dp)]
    endif
  enddo
  
  output = ComplexModeForce(                                &
     & output_modes,                                        &
     & this%transform_complex_magnitudes( output_modes,     &
     &                                    magnitudes,       &
     &                                    qpoints,          &
     &                                    inverse = .true. ))
end function

function inverse_transform_RealModeDisplacement(this,displacement,modes, &
   & qpoints) result(output)
  implicit none
  
  class(VscfRvectors),        intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  type(RealMode),             intent(in) :: modes(:)
  type(QpointData),           intent(in) :: qpoints(:)
  type(RealModeDisplacement)             :: output
  
  type(RealMode), allocatable :: output_modes(:)
  real(dp),       allocatable :: magnitudes(:)
  
  integer :: i
  
  output_modes = select_modes(displacement%vectors, modes)
  magnitudes = displacement%vectors%magnitude
  
  do i=1,size(output_modes)
    if (.not. any(output_modes%id==output_modes(i)%paired_id)) then
      output_modes = [ output_modes,                                     &
                     & modes(first(modes%id==output_modes(i)%paired_id)) ]
      magnitudes = [magnitudes, 0.0_dp]
    endif
  enddo
  
  output = RealModeDisplacement(                         &
     & output_modes,                                     &
     & this%transform_real_magnitudes( output_modes,     &
     &                                 magnitudes,       &
     &                                 qpoints,          &
     &                                 inverse = .true. ))
end function

function inverse_transform_RealModeForce(this,force,modes,qpoints) &
   & result(output)
  implicit none
  
  class(VscfRvectors), intent(in) :: this
  type(RealModeForce), intent(in) :: force
  type(RealMode),      intent(in) :: modes(:)
  type(QpointData),    intent(in) :: qpoints(:)
  type(RealModeForce)             :: output
  
  type(RealMode), allocatable :: output_modes(:)
  real(dp),       allocatable :: magnitudes(:)
  
  integer :: i
  
  output_modes = select_modes(force%vectors, modes)
  magnitudes = force%vectors%magnitude
  
  do i=1,size(output_modes)
    if (.not. any(output_modes%id==output_modes(i)%paired_id)) then
      output_modes = [ output_modes,                                     &
                     & modes(first(modes%id==output_modes(i)%paired_id)) ]
      magnitudes = [magnitudes, 0.0_dp]
    endif
  enddo
  
  output = RealModeForce(                                &
     & output_modes,                                     &
     & this%transform_real_magnitudes( output_modes,     &
     &                                 magnitudes,       &
     &                                 qpoints,          &
     &                                 inverse = .true. ))
end function

function transform_complex_magnitudes(this,modes,magnitudes,qpoints, &
   & inverse) result(output)
  implicit none
  
  class(VscfRvectors), intent(in)  :: this
  type(ComplexMode),   intent(in)  :: modes(:)
  complex(dp),         intent(in)  :: magnitudes(:)
  type(QpointData),    intent(in)  :: qpoints(:)
  logical,             intent(in)  :: inverse
  complex(dp), allocatable         :: output(:)
  
  type(ComplexMode) :: mode
  type(QpointData)  :: qpoint
  type(IntVector)   :: rvector
  
  integer :: i,j,ialloc
  
  allocate(output(size(magnitudes)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    mode = modes(i)
    qpoint = select_qpoint(mode, qpoints)
    
    j = first(this%vscf_rvectors%subspace_id==mode%subspace_id, default=0)
    if (j==0) then
      ! This mode is unaffected by the VSCF R-vector translation.
      output(i) = magnitudes(i)
      cycle
    endif
    
    rvector = this%vscf_rvectors(j)%rvector
    if (inverse) then
      rvector = -rvector
    endif
    
    ! u_q -> u_q * exp{2*pi*i*q.R}
    output(i) = magnitudes(i) &
              * exp_2pii(qpoint%qpoint*rvector)
  enddo
end function

function transform_real_magnitudes(this,modes,magnitudes,qpoints, &
   & inverse) result(output)
  implicit none
  
  class(VscfRvectors),   intent(in) :: this
  type(RealMode),        intent(in) :: modes(:)
  real(dp),              intent(in) :: magnitudes(:)
  type(QpointData),      intent(in) :: qpoints(:)
  logical,               intent(in) :: inverse
  real(dp), allocatable             :: output(:)
  
  type(RealMode)   :: mode
  type(QpointData) :: qpoint
  type(IntVector)  :: rvector
  
  logical, allocatable :: mode_transformed(:)
  
  integer :: i,j,ialloc
  
  allocate( output(size(magnitudes)),           &
          & mode_transformed(size(magnitudes)), &
          & stat=ialloc); call err(ialloc)
  mode_transformed = .false.
  do i=1,size(magnitudes)
    if (mode_transformed(i)) then
      cycle
    endif
    
    mode   = modes(i)
    qpoint = select_qpoint(mode, qpoints)
    
    j = first(this%vscf_rvectors%subspace_id==mode%subspace_id, default=0)
    if (j==0) then
      ! This mode is unaffected by the VSCF R-vector translation.
      output(i) = magnitudes(i)
      mode_transformed(i) = .true.
      cycle
    endif
    
    rvector = this%vscf_rvectors(j)%rvector
    if (inverse) then
      rvector = -rvector
    endif
    
    ! Calculate the transformed vector along vector i,
    !    and along the mode paired to vector i.
    j = first(modes%id==mode%paired_id)
    if (j==i) then
      ! The mode is its own pair, so cos(2*pi*q.R) = +/-1.
      ! u -> cos(2*pi*q.r) * u.
      output(i) = cos_2pi(qpoint%qpoint*rvector) * magnitudes(i)
    else
      if (mode%id<mode%paired_id) then
        ! Displacement i is u_c, vector j is u_s.
        output(i) = cos_2pi(qpoint%qpoint*rvector)*magnitudes(i) &
                & - sin_2pi(qpoint%qpoint*rvector)*magnitudes(j)
        output(j) = cos_2pi(qpoint%qpoint*rvector)*magnitudes(j) &
                & + sin_2pi(qpoint%qpoint*rvector)*magnitudes(i)
      else
        ! Displacement i is u_s, vector j is u_c.
        output(i) = cos_2pi(qpoint%qpoint*rvector)*magnitudes(i) &
                & + sin_2pi(qpoint%qpoint*rvector)*magnitudes(j)
        output(j) = cos_2pi(qpoint%qpoint*rvector)*magnitudes(j) &
                & - sin_2pi(qpoint%qpoint*rvector)*magnitudes(i)
      endif
    endif
    
    mode_transformed(i) = .true.
    mode_transformed(j) = .true.
  enddo
  
  ! Check that all modes have been transformed.
  if (.not. all(mode_transformed)) then
    call print_line(ERROR//': error transforming modes.')
    call err()
  endif
end function

! ----------------------------------------------------------------------
! Construct VSCF R-vectors.
! ----------------------------------------------------------------------
! Constructs the list of all R-vectors for each subspace,
!    and then calls a recursive helper function to generate all permutations.
function construct_vscf_rvectors(sampling_point,supercell,real_modes, &
   & qpoints) result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in) :: sampling_point
  type(StructureData),        intent(in) :: supercell
  type(RealMode),             intent(in) :: real_modes(:)
  type(QpointData),           intent(in) :: qpoints(:)
  type(VscfRvectors), allocatable        :: output(:)
  
  type(RvectorArray), allocatable :: subspace_rvectors(:)
  
  subspace_rvectors = construct_rvector_arrays( sampling_point, &
                                              & supercell,      &
                                              & real_modes,     &
                                              & qpoints)
  output = list_rvector_permutations(subspace_rvectors)
end function

! First helper function for construct_vscf_rvectors.
! Constructs the R-vectors at which each subspace needs to be sampled.
function construct_rvector_arrays(sampling_point,supercell,real_modes, &
   & qpoints) result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in) :: sampling_point
  type(StructureData),        intent(in) :: supercell
  type(RealMode),             intent(in) :: real_modes(:)
  type(QpointData),           intent(in) :: qpoints(:)
  type(RvectorArray), allocatable        :: output(:)
  
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
  modes = select_modes(sampling_point%vectors, real_modes)
  subspace_ids = modes%subspace_id
  
  ! List q-points corresponding to said modes.
  allocate(mode_qpoints(size(modes)), stat=ialloc); call err(ialloc)
  do i=1,size(modes)
    mode_qpoints(i) = qpoints(first(qpoints%id==modes(i)%qpoint_id_plus))
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
    
    output(i) = RvectorArray(subspace_ids(i),rvectors(:no_rvectors))
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

! Second helper function for construct_vscf_rvectors.
! Recursively constructs the permutations of R-vectors across subspaces.
! N.B. rvectors_in should only be included when this function calls itself.
recursive function list_rvector_permutations(rvector_arrays,vscf_rvectors_in) &
   & result(output)
  implicit none
  
  type(RvectorArray), intent(in)           :: rvector_arrays(:)
  type(VscfRvectors), intent(in), optional :: vscf_rvectors_in
  type(VscfRvectors), allocatable          :: output(:)
  
  type(VscfRvectors) :: vscf_rvectors
  type(VscfRvector)  :: vscf_rvector
  
  integer :: i,ialloc
  
  if (present(vscf_rvectors_in)) then
    vscf_rvectors = vscf_rvectors_in
  else
    vscf_rvectors = VscfRvectors()
  endif
  
  allocate(output(0), stat=ialloc); call err(ialloc)
  
  do i=1,size(rvector_arrays(1))
    vscf_rvector = VscfRvector( rvector_arrays(1)%subspace_id, &
                              & rvector_arrays(1)%rvectors(i))
    if (size(rvector_arrays)==1) then
      output = [output, vscf_rvectors//vscf_rvector]
    else
      output = [                                                         &
             &   output,                                                 &
             &   list_rvector_permutations( rvector_arrays(2:),          &
             &                              vscf_rvectors//vscf_rvector) &
             & ]
    endif
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_VscfRvectors(this,input)
  implicit none
  
  class(VscfRvectors), intent(out) :: this
  type(String),        intent(in)  :: input(:)
  
  select type(this); type is(VscfRvectors)
    this = VscfRvectors(VscfRvector(input))
  class default
    call err()
  end select
end subroutine

function write_VscfRvectors(this) result(output)
  implicit none
  
  class(VscfRvectors), intent(in) :: this
  type(String), allocatable       :: output(:)
  
  select type(this); type is(VscfRvectors)
    output = str(this%vscf_rvectors)
  class default
    call err()
  end select
end function

function new_VscfRvectors_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(VscfRvectors)       :: this
  
  call this%read(input)
end function

impure elemental function new_VscfRvectors_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(VscfRvectors)            :: this
  
  this = VscfRvectors(str(input))
end function
end module
