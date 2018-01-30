! ======================================================================
! q-points of the primitive cell.
! ======================================================================
module qpoints_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use linear_algebra_module
  use fraction_algebra_module
  implicit none
  
  type QpointData
    ! The q-point in fractional primitive reciprocal space co-ordinates.
    type(FractionVector) :: qpoint
    
    ! Whether or not 2*q is a G-vector of the primitive cell.
    logical :: is_paired_qpoint
    
    ! The id of q-point q' such that q+q' is a G-vector of the primitive cell.
    integer :: paired_qpoint
  contains
    ! The smallest value of |S| s.t. S.q is a vector of integers,
    !    i.e. such that this q-point is a G-vector of the supercell S.
    procedure, public :: min_sc_size
  end type
contains

! ----------------------------------------------------------------------
! Returns the smallest integer, s, s.t. there exists the matrix S
!    with determinant |S|=s s.t. S.q is a vector of integers.
! i.e. The supercell with supercell matrix S has this q-point as a G-vector.
! ----------------------------------------------------------------------
function min_sc_size(this) result(output)
  use utils_module, only : lcm
  implicit none
  
  class(QpointData), intent(in) :: this
  integer                       :: output
  
  type(IntFraction) :: q(3)
  
  q = frac(this%qpoint)
  output = lcm( q(1)%denominator(), &
              & q(2)%denominator(), &
              & q(3)%denominator()  )
end function

! ----------------------------------------------------------------------
! Generates the set of q-points of the input structure which correspond
!    to G-vectors of the large supercell.
! ----------------------------------------------------------------------
function generate_qpoints(structure,large_supercell) result(output)
  use linear_algebra_module
  use fraction_algebra_module
  use structure_module
  use group_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: large_supercell
  type(QpointData), allocatable   :: output(:)
  
  ! Working variables
  integer, allocatable :: paired_qpoints(:)
  type(FractionVector) :: rotated_qpoint
  
  ! Temporary variables
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Construct q-points from G-vectors of large supercell.
  ! --------------------------------------------------
  allocate(output(large_supercell%sc_size), stat=ialloc); call err(ialloc)
  do i=1,large_supercell%sc_size
    output(i)%qpoint = transpose(large_supercell%recip_supercell) &
                   & * large_supercell%gvectors(i)
  enddo
  
  ! --------------------------------------------------
  ! Find paired q-points.
  ! --------------------------------------------------
  ! qpoint + paired_qpoint = G, for a primitive-cell G-vector.
  allocate( paired_qpoints(large_supercell%sc_size), &
          & stat=ialloc); call err(ialloc)
  paired_qpoints = 0
  do i=1,size(output)
    do j=1,size(output)
      if (is_int(output(i)%qpoint+output(j)%qpoint)) then
        if (paired_qpoints(i)==0 .and. paired_qpoints(j)==0) then
          paired_qpoints(i) = j
          paired_qpoints(j) = i
        else
          if (paired_qpoints(i)/=j .or. paired_qpoints(j)/=i) then
            call print_line(CODE_ERROR//': error pairing q-points.')
          endif
        endif
      endif
    enddo
  enddo
  
  if (any(paired_qpoints==0)) then
    call print_line(CODE_ERROR//': q-points were not succesfully paired up.')
    call err()
  endif
  
  do i=1,size(output)
    if (paired_qpoints(i)==i) then
      output(i)%is_paired_qpoint = .true.
    else
      output(i)%is_paired_qpoint = .false.
    endif
    
    output(i)%paired_qpoint = paired_qpoints(i)
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine write_qpoints_file(this,filename)
  use ofile_module
  implicit none
  
  type(QpointData), intent(in), allocatable :: this(:)
  type(String),     intent(in)              :: filename
  
  type(OFile) :: qpoints_file
  
  integer :: i
  
  qpoints_file = filename
  do i=1,size(this)
    call qpoints_file%print_line( 'q-point '//i//', q=(qx/nx, qy/ny, qz/nz):')
    call qpoints_file%print_line( this(i)%qpoint)
    call qpoints_file%print_line( 'Is 2q equal to a primitive G-vector?:')
    call qpoints_file%print_line( this(i)%is_paired_qpoint)
    call qpoints_file%print_line( "The ID of q' s.t. q+q' is a primitive G-vector")
    call qpoints_file%print_line( this(i)%paired_qpoint)
    call qpoints_file%print_line( '')
  enddo
end subroutine

function read_qpoints_file(filename) result(this)
  use ifile_module
  implicit none
  
  type(String), intent(in)      :: filename
  type(QpointData), allocatable :: this(:)
  
  type(IFile) :: qpoints_file
  
  integer :: no_qpoints
  integer :: qpoint_line
  integer :: i,j,ialloc
  
  qpoints_file = filename
  no_qpoints = size(qpoints_file)/7
  allocate(this(no_qpoints), stat=ialloc); call err(ialloc)
  do i=1,no_qpoints
    qpoint_line = (i-1)*7
    
    this(i)%qpoint           = frac( qpoints_file%split_line(qpoint_line+2) )
    this(i)%is_paired_qpoint = lgcl( qpoints_file%line(      qpoint_line+4) )
    this(i)%paired_qpoint    = int(  qpoints_file%line(      qpoint_line+6) )
  enddo
end function
end module
