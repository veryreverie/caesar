! ======================================================================
! Harmonic normal modes, in complex co-ordinates.
! ======================================================================
module complex_mode_submodule
  use utils_module
  
  use structure_module
  implicit none
  
  private
  
  public :: ComplexMode
  public :: operator(*)
  public :: l2_norm
  public :: conjg
  
  ! A normal mode in complex co-ordinates.
  type, extends(Stringsable) :: ComplexMode
    ! An id which is unique to each mode, and the unique id of the equivalent
    !    mode at the q-point q' s.t. q+q' is a G-vector.
    ! If 2q is a G-vector, then paired_id=id.
    integer :: id
    integer :: paired_id
    
    ! The frequency, and frequency-relevant information.
    real(dp) :: frequency
    logical  :: soft_mode          ! True if frequency < -epsilon.
    logical  :: translational_mode ! True if frequency=0 and at gamma.
    
    ! The displacements of atoms in the primitive cell.
    type(ComplexVector), allocatable :: primitive_displacements(:)
    
    ! The ID of the q-point at which this mode exits.
    integer :: qpoint_id
    
    ! The ID of the subspace modes which are degenerate with this mode.
    integer :: subspace_id
  contains
    procedure, public :: read  => read_ComplexMode
    procedure, public :: write => write_ComplexMode
  end type
  
  interface ComplexMode
    module procedure new_ComplexMode
  end interface
  
  interface operator(*)
    module procedure dot_ComplexMode_ComplexMode
  end interface
  
  interface l2_norm
    module procedure l2_norm_ComplexMode
  end interface
  
  interface conjg
    module procedure conjg_ComplexMode
  end interface
contains

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function new_ComplexMode(id,paired_id,frequency,soft_mode,translational_mode, &
   & primitive_displacements,qpoint_id,subspace_id) result(this)
  implicit none
  
  integer,             intent(in) :: id
  integer,             intent(in) :: paired_id
  real(dp),            intent(in) :: frequency
  logical,             intent(in) :: soft_mode
  logical,             intent(in) :: translational_mode
  type(ComplexVector), intent(in) :: primitive_displacements(:)
  integer,             intent(in) :: qpoint_id
  integer,             intent(in) :: subspace_id
  type(ComplexMode)               :: this
  
  this%id                      = id
  this%paired_id               = paired_id
  this%frequency               = frequency
  this%soft_mode               = soft_mode
  this%translational_mode      = translational_mode
  this%primitive_displacements = primitive_displacements
  this%qpoint_id               = qpoint_id
  this%subspace_id             = subspace_id
end function

! ----------------------------------------------------------------------
! Find the dot product between two complex modes.
! ----------------------------------------------------------------------
impure elemental function dot_ComplexMode_ComplexMode(this,that) result(output)
  implicit none
  
  type(ComplexMode), intent(in) :: this
  type(ComplexMode), intent(in) :: that
  complex(dp)                   :: output
  
  output = sum( this%primitive_displacements &
            & * that%primitive_displacements )
end function

impure elemental function l2_norm_ComplexMode(this) result(output)
  implicit none
  
  type(ComplexMode), intent(in) :: this
  real(dp)                      :: output
  
  output = sqrt(real(this*conjg(this)))
end function

! ----------------------------------------------------------------------
! Take the conjugate of the mode.
! This is equivalent to the transformation q -> -q.
! ----------------------------------------------------------------------
impure elemental function conjg_ComplexMode(input) result(output)
  implicit none
  
  type(ComplexMode), intent(in) :: input
  type(ComplexMode)             :: output
  
  output                         = input
  output%id                      = input%paired_id
  output%paired_id               = input%id
  output%primitive_displacements = conjg(input%primitive_displacements)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_ComplexMode(this,input)
  implicit none
  
  class(ComplexMode), intent(out) :: this
  type(String),       intent(in)  :: input(:)
  
  integer                          :: id
  integer                          :: paired_id
  real(dp)                         :: frequency
  logical                          :: soft_mode
  logical                          :: translational_mode
  type(ComplexVector), allocatable :: primitive_displacements(:)
  integer                          :: qpoint_id
  integer                          :: subspace_id
  
  type(String), allocatable :: line(:)
  
  integer :: no_atoms
  
  integer :: i,ialloc
  
  select type(this); type is(ComplexMode)
    ! Read the id of this mode.
    line = split_line(input(1))
    id = int(line(4))
    
    ! Read the id of this mode's pair.
    line = split_line(input(2))
    paired_id = int(line(6))
    
    ! Read frequency.
    line = split_line(input(3))
    frequency = dble(line(4))
    
    ! Read whether or not this mode is soft.
    line = split_line(input(4))
    soft_mode = lgcl(line(5))
    
    ! Read whether or not this mode is purely translational.
    line = split_line(input(5))
    translational_mode = lgcl(line(5))
    
    ! Read the qpoint id of this mode.
    line = split_line(input(6))
    qpoint_id = int(line(4))
    
    ! Read the degeneracy id of this mode.
    line = split_line(input(7))
    subspace_id = int(line(4))
    
    ! Read in the displacement associated with the mode.
    no_atoms = size(input)-8
    allocate( primitive_displacements(no_atoms), &
            & stat=ialloc); call err(ialloc)
    do i=1,no_atoms
      line = split_line(input(8+i))
      primitive_displacements(i) = cmplx(line)
    enddo
    
    this = ComplexMode( id,                      &
                      & paired_id,               &
                      & frequency,               &
                      & soft_mode,               &
                      & translational_mode,      &
                      & primitive_displacements, &
                      & qpoint_id,               &
                      & subspace_id)
  end select
end subroutine

function write_ComplexMode(this) result(output)
  implicit none
  
  class(ComplexMode), intent(in) :: this
  type(String), allocatable      :: output(:)
  
  integer :: i,ialloc
  
  
  select type(this); type is(ComplexMode)
    allocate( output(8+size(this%primitive_displacements)), &
            & stat=ialloc); call err(ialloc)
    output(1) = 'Mode ID                   : '//this%id
    output(2) = 'ID of paired mode         : '//this%paired_id
    output(3) = 'Mode frequency            : '//this%frequency
    output(4) = 'Mode is soft              : '//this%soft_mode
    output(5) = 'Mode purely translational : '//this%translational_mode
    output(6) = 'q-point id                : '//this%qpoint_id
    output(7) = 'Degeneracy id             : '//this%subspace_id
    output(8) = 'Displacements in primitive cell:'
    do i=1,size(this%primitive_displacements)
      output(8+i) = str(this%primitive_displacements(i))
    enddo
  end select
end function
end module
