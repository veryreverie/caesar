! ======================================================================
! As complex_mode_module, but in real co-ordinates.
! ======================================================================
module real_mode_submodule
  use utils_module
  
  use structure_module
  
  use cartesian_displacement_submodule
  use cartesian_force_submodule
  use mass_weighted_displacement_submodule
  use mass_weighted_force_submodule
  implicit none
  
  private
  
  public :: RealMode
  public :: MassWeightedDisplacement
  public :: MassWeightedForce
  public :: CartesianDisplacement
  public :: CartesianForce
  public :: select_qpoint
  public :: select_qpoints
  
  ! A normal mode in real co-ordinates.
  type, extends(Stringsable) :: RealMode
    integer :: id
    integer :: paired_id
    
    ! The frequency, and frequency-relevant information.
    real(dp) :: frequency
    real(dp) :: spring_constant    ! k, from w=sqrt(k/m), where w=frequency.
    logical  :: soft_mode          ! True if frequency < -epsilon.
    logical  :: translational_mode ! True if frequency=0 and at gamma.
    
    ! Unit vectors along the normal mode, in cartesian and mass-weighted
    !    co-ordinates.
    ! Normal modes are orthogonal in mass-weighted co-ordinates, so in general
    !    they are not orthogonal in cartesian co-ordinates.
    ! In both cases, the unit vectors are stored as displacements from
    !    equilibrium for each atom in the primitive cell.
    type(RealVector), allocatable :: mass_weighted_vector(:)
    type(RealVector), allocatable :: cartesian_vector(:)
    
    ! The IDs of the q-points at which this mode exists.
    ! N.B. unlike complex modes, real modes are not localised to a single
    !    q-point, but rather a superposition across +/- pair.
    integer :: qpoint_id_plus
    integer :: qpoint_id_minus
    
    ! The ID of the subspace of modes which are degenerate with this mode.
    integer :: subspace_id
  contains
    procedure, private :: construct_vector
    ! I/O.
    procedure, public :: read  => read_RealMode
    procedure, public :: write => write_RealMode
  end type
  
  interface RealMode
    module procedure new_RealMode
    module procedure new_RealMode_StringArray
  end interface
  
  interface MassWeightedDisplacement
    module procedure new_MassWeightedDisplacement_RealMode
  end interface
  
  interface MassWeightedForce
    module procedure new_MassWeightedForce_RealMode
  end interface
  
  interface CartesianDisplacement
    module procedure new_CartesianDisplacement_RealMode
  end interface
  
  interface CartesianForce
    module procedure new_CartesianForce_RealMode
  end interface
  
  interface select_qpoint
    module procedure select_qpoint_RealMode
  end interface
  
  interface select_qpoints
    module procedure select_qpoints_RealModes
  end interface
contains

! ----------------------------------------------------------------------
! Basic constructor.
! ----------------------------------------------------------------------
function new_RealMode(id,paired_id,frequency,spring_constant,soft_mode,       &
   & translational_mode,mass_weighted_vector,cartesian_vector,qpoint_id_plus, &
   & qpoint_id_minus,subspace_id) result(this)
  implicit none
  
  integer,          intent(in) :: id
  integer,          intent(in) :: paired_id
  real(dp),         intent(in) :: frequency
  real(dp),         intent(in) :: spring_constant
  logical,          intent(in) :: soft_mode
  logical,          intent(in) :: translational_mode
  type(RealVector), intent(in) :: mass_weighted_vector(:)
  type(RealVector), intent(in) :: cartesian_vector(:)
  integer,          intent(in) :: qpoint_id_plus
  integer,          intent(in) :: qpoint_id_minus
  integer,          intent(in) :: subspace_id
  type(RealMode)               :: this
  
  this%id                   = id
  this%paired_id            = paired_id
  this%frequency            = frequency
  this%spring_constant      = spring_constant
  this%soft_mode            = soft_mode
  this%translational_mode   = translational_mode
  this%mass_weighted_vector = mass_weighted_vector
  this%cartesian_vector     = cartesian_vector
  this%qpoint_id_plus       = qpoint_id_plus
  this%qpoint_id_minus      = qpoint_id_minus
  this%subspace_id          = subspace_id
end function

! ----------------------------------------------------------------------
! Return the mode in the mass-weighted or cartesian co-ordinates of
!    a given supercell.
! ----------------------------------------------------------------------
impure elemental function new_MassWeightedDisplacement_RealMode(this, &
   & structure,qpoint) result(output)
  implicit none
  
  class(RealMode),     intent(in) :: this
  type(StructureData), intent(in) :: structure
  type(QpointData),    intent(in) :: qpoint
  type(MassWeightedDisplacement)  :: output
  
  output = MassWeightedDisplacement(this%construct_vector( &
                              & this%mass_weighted_vector, &
                              & structure,                 &
                              & qpoint                     ))
end function

impure elemental function new_MassWeightedForce_RealMode(this,structure, &
   & qpoint) result(output)
  implicit none
  
  class(RealMode),     intent(in) :: this
  type(StructureData), intent(in) :: structure
  type(QpointData),    intent(in) :: qpoint
  type(MassWeightedForce)         :: output
  
  output = MassWeightedForce(this%construct_vector( &
                       & this%mass_weighted_vector, &
                       & structure,                 &
                       & qpoint                     ))
end function

impure elemental function new_CartesianDisplacement_RealMode(this,structure, &
   & qpoint) result(output)
  implicit none
  
  class(RealMode),     intent(in) :: this
  type(StructureData), intent(in) :: structure
  type(QpointData),    intent(in) :: qpoint
  type(CartesianDisplacement)     :: output
  
  output = CartesianDisplacement( MassWeightedDisplacement( this,         &
                                &                           structure,    &
                                &                           qpoint     ), &
                                & structure                               )
end function

impure elemental function new_CartesianForce_RealMode(this,structure,qpoint) &
   & result(output)
  implicit none
  
  class(RealMode),     intent(in) :: this
  type(StructureData), intent(in) :: structure
  type(QpointData),    intent(in) :: qpoint
  type(CartesianForce)            :: output
  
  output = CartesianForce( MassWeightedForce( this,         &
                         &                    structure,    &
                         &                    qpoint     ), &
                         & structure                        )
end function

function construct_vector(this,prim_vectors,structure,qpoint) result(output)
  implicit none
  
  class(RealMode),     intent(in) :: this
  type(RealVector),    intent(in) :: prim_vectors(:)
  type(StructureData), intent(in) :: structure
  type(QpointData),    intent(in) :: qpoint
  type(RealVector), allocatable   :: output(:)
  
  type(AtomData)   :: atom
  real(dp)         :: qr
  
  integer :: i,ialloc
  
  if (qpoint%id/=this%qpoint_id_plus.and.qpoint%id/=this%qpoint_id_minus) then
    call print_line(CODE_ERROR//': Mode and q-point incompatible.')
    call err()
  endif
  
  allocate(output(structure%no_atoms), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    atom         = structure%atoms(i)
    
    if (qpoint%id==this%qpoint_id_plus) then
      qr =  dble(qpoint%qpoint*structure%rvectors(atom%rvec_id()))
    else
      qr = -dble(qpoint%qpoint*structure%rvectors(atom%rvec_id()))
    endif
    
    if (this%id<=this%paired_id) then
      ! This is a cos(2 pi q.R) mode.
      output(i) = prim_vectors(atom%prim_id()) * cos(2*PI*qr)
    else
      ! This is a sin(2 pi q.R) mode.
      output(i) = prim_vectors(atom%prim_id()) * sin(2*PI*qr)
    endif
  enddo
end function

! ----------------------------------------------------------------------
! Select q-points corresponding to a given mode or modes.
! ----------------------------------------------------------------------
function select_qpoint_RealMode(mode,qpoints) result(output)
  implicit none
  
  type(RealMode),   intent(in) :: mode
  type(QpointData), intent(in) :: qpoints(:)
  type(QpointData)             :: output
  
  output = qpoints(first(qpoints%id==mode%qpoint_id_plus))
end function

function select_qpoints_RealModes(modes,qpoints) result(output)
  implicit none
  
  type(RealMode),   intent(in)  :: modes(:)
  type(QpointData), intent(in)  :: qpoints(:)
  type(QpointData), allocatable :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(modes)), stat=ialloc); call err(ialloc)
  do i=1,size(modes)
    output(i) = select_qpoint(modes(i), qpoints)
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_RealMode(this,input)
  implicit none
  
  class(RealMode), intent(out) :: this
  type(String),    intent(in)  :: input(:)
  
  integer                       :: id
  integer                       :: paired_id
  real(dp)                      :: frequency
  real(dp)                      :: spring_constant
  logical                       :: soft_mode
  logical                       :: translational_mode
  type(RealVector), allocatable :: mass_weighted_vector(:)
  type(RealVector), allocatable :: cartesian_vector(:)
  integer                       :: qpoint_id_plus
  integer                       :: qpoint_id_minus
  integer                       :: subspace_id
  
  type(String), allocatable :: line(:)
  
  integer :: no_atoms
  
  select type(this); type is(RealMode)
    ! Read the id of this mode.
    line = split_line(input(1))
    id = int(line(4))
    
    ! Read the id of this mode's pair.
    line = split_line(input(2))
    paired_id = int(line(6))
    
    ! Read frequency and spring constant.
    line = split_line(input(3))
    frequency = dble(line(4))
    
    line = split_line(input(4))
    spring_constant = dble(line(4))
    
    ! Read whether or not this mode is soft.
    line = split_line(input(5))
    soft_mode = lgcl(line(5))
    
    ! Read whether or not this mode is purely translational.
    line = split_line(input(6))
    translational_mode = lgcl(line(5))
    
    ! Read the q-qpoint id of the plus mode.
    line = split_line(input(7))
    qpoint_id_plus = int(line(5))
    
    ! Read the q-qpoint id of the minus mode.
    line = split_line(input(8))
    qpoint_id_minus = int(line(5))
    
    ! Read the degeneracy id of this mode.
    line = split_line(input(9))
    subspace_id = int(line(4))
    
    ! Read in the vector associated with the mode.
    no_atoms = (size(input)-11)/2
    mass_weighted_vector = RealVector(input(11:10+no_atoms))
    cartesian_vector = RealVector(input(size(input)-no_atoms+1:))
    
    this = RealMode( id,                   &
                   & paired_id,            &
                   & frequency,            &
                   & spring_constant,      &
                   & soft_mode,            &
                   & translational_mode,   &
                   & mass_weighted_vector, &
                   & cartesian_vector,     &
                   & qpoint_id_plus,       &
                   & qpoint_id_minus,      &
                   & subspace_id           )
  end select
end subroutine

function write_RealMode(this) result(output)
  implicit none
  
  class(RealMode), intent(in) :: this
  type(String), allocatable   :: output(:)
  
  select type(this); type is(RealMode)
    output = [ 'Mode ID                   : '//this%id,                 &
             & 'ID of paired mode         : '//this%paired_id,          &
             & 'Mode frequency            : '//this%frequency,          &
             & 'Spring constant           : '//this%spring_constant,    &
             & 'Mode is soft              : '//this%soft_mode,          &
             & 'Mode purely translational : '//this%translational_mode, &
             & 'Positive q-point id       : '//this%qpoint_id_plus,     &
             & 'Negative q-point id       : '//this%qpoint_id_minus,    &
             & 'Subspace id               : '//this%subspace_id,        &
             & str('Mass-weighted displacements in primitive cell:'),   &
             & str(this%mass_weighted_vector),                          &
             & str('Cartesian displacements in primitive cell:'),       &
             & str(this%cartesian_vector)                               ]
  end select
end function

impure elemental function new_RealMode_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(RealMode)                :: this
  
  this = input
end function
end module
