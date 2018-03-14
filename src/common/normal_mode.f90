! ======================================================================
! Harmonic normal modes.
! ======================================================================
!
! Most calculation happens in normal mode co-ordinates, but setting up
!    DFT calculations requires conversion to cartesian, and the results of
!    DFT have to be converted from cartesian forces.
!
! - Co-ordinates can be converted between real and normal modes via the
!      relevant class constructors.
! - Displacements in real mode co-ordinates can be converted into cartesian
!      co-ordinates using real_mode_to_displacement.
! - Forces can be converted into real mode co-ordinates using
!      force_to_real_mode.
!
! d(R,j) = sqrt(m_j)(r(R,j)-r0(R,j))
! f(R,j) = -dV/dd(R,j)
!
! x(q,j) = 1/N * sum_R [ d(R,j) * exp( i q.R) ]
! d(R,j) =       sum_q [ x(q,j) * exp(-i q.R) ]
!
! N.B. if 2q=G then q.R is a half-integer.
!    => exp(i q.R) = cos(q.R) = +/- 1
!    => x(q,j) is real for all j.
!
! x(q,j) = 1/N       * sum_R [ d(R,j) * cos(q.R) ]           (for 2q=G only)
! c(q,j) = sqrt(2)/N * sum_R [ d(R,j) * cos(q.R) ]           (for 2q/=G)
! s(q,j) = sqrt(2)/N * sum_R [ d(R,j) * sin(q.R) ]           (for 2q/=G)
! d(R,j) =           sum_q [x_q * cos(q,R)]                  (sum over 2q=G)
!        + sqrt(2) * sum_q [c_q * cos(q.R) + s_q * sin(q.R)] (sum over 2q/=G)
!
! ux(q,k) = sum_j [ U(q,j,k) * x(q,j) ]        (for 2q=G only)
! uc(q,k) = sum_j [ U(q,j,k) * c(q,j) ]        (for 2q/=G)
! us(q,k) = sum_j [ U(q,j,k) * s(q,j) ]        (for 2q/=G)
module normal_mode_module
  use utils_module
  
  use qpoints_module
  implicit none
  
  private
  
  public :: ComplexMode
  
  public :: operator(*)
  public :: l2_norm
  
  ! A normal mode in complex co-ordinates.
  type ComplexMode
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
    
    ! An id which is shared between degenerate states, and different otherwise.
    integer :: degeneracy_id
  contains
    procedure, public :: write_file => write_file_ComplexMode
  end type
  
  interface ComplexMode
    module procedure read_file_ComplexMode
  end interface
  
  interface operator(*)
    module procedure dot_ComplexModes
  end interface
  
  interface l2_norm
    module procedure l2_norm_ComplexMode
  end interface
contains

! ----------------------------------------------------------------------
! ComplexMode procedures.
! ----------------------------------------------------------------------
subroutine write_file_ComplexMode(this,filename)
  implicit none
  
  class(ComplexMode), intent(in) :: this
  type(String),       intent(in) :: filename
  
  type(OFile) :: mode_file
  
  integer :: i
  
  mode_file = filename
  call mode_file%print_line('Mode ID                   : '//this%id)
  call mode_file%print_line('ID of paired mode         : '//this%paired_id)
  call mode_file%print_line('Mode frequency            : '//this%frequency)
  call mode_file%print_line('Mode is soft              : '//this%soft_mode)
  call mode_file%print_line('Mode purely translational : '// &
     & this%translational_mode)
  call mode_file%print_line('Degeneracy id             : '//this%degeneracy_id)
  call mode_file%print_line('Displacements in primitive cell:')
  do i=1,size(this%primitive_displacements)
    call mode_file%print_line(this%primitive_displacements(i))
  enddo
end subroutine

function read_file_ComplexMode(filename) result(this)
  implicit none
  
  type(String), intent(in) :: filename
  type(ComplexMode)        :: this
  
  type(IFile)               :: mode_file
  type(String), allocatable :: line(:)
  integer                   :: no_atoms
  
  integer :: i,ialloc
  
  ! Read in mode file.
  mode_file = filename
  
  ! Read the id of this mode.
  line = split(mode_file%line(1))
  this%id = int(line(4))
  
  ! Read the id of this mode's pair.
  line = split(mode_file%line(2))
  this%paired_id = int(line(6))
  
  ! Read frequency.
  line = split(mode_file%line(3))
  this%frequency = dble(line(4))
  
  ! Read whether or not this mode is soft.
  line = split(mode_file%line(4))
  this%soft_mode = lgcl(line(5))
  
  ! Read whether or not this mode is purely translational.
  line = split(mode_file%line(5))
  this%translational_mode = lgcl(line(5))
  
  ! Read the degeneracy id of this mode.
  line = split(mode_file%line(6))
  this%degeneracy_id = int(line(4))
  
  ! Read in the displacement associated with the mode.
  no_atoms = size(mode_file)-7
  allocate( this%primitive_displacements(no_atoms), &
          & stat=ialloc); call err(ialloc)
  do i=1,no_atoms
    line = split(mode_file%line(7+i))
    this%primitive_displacements(i) = cmplx(line)
  enddo
end function

! ----------------------------------------------------------------------
! Find the dot product between two complex modes.
! ----------------------------------------------------------------------
function dot_ComplexModes(this,that) result(output)
  implicit none
  
  type(ComplexMode), intent(in) :: this
  type(ComplexMode), intent(in) :: that
  complex(dp)                   :: output
  
  output = sum( this%primitive_displacements &
            & * conjg(that%primitive_displacements) )
end function

function l2_norm_ComplexMode(this) result(output)
  implicit none
  
  type(ComplexMode), intent(in) :: this
  real(dp)                      :: output
  
  output = sqrt(real(this*this))
end function
end module
