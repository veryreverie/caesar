! This program determines which directionss are related by symmetry,
!    and thus which atoms should be displaced in which directions 
!    to construct the matrix of force constants.

! The symmetry operations to consider are nicely outlined here:
! http://www.homepages.ucl.ac.uk/~ucfbdxa/phon/node4.html

module unique_directions_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  type UniqueDirections
    integer, allocatable :: unique_atoms(:)
    integer, allocatable :: xy_symmetry(:)
    integer, allocatable :: xz_symmetry(:)
    integer, allocatable :: yz_symmetry(:)
  end type
  
  interface new
    module procedure new_UniqueDirections
  end interface
  
  interface drop
    module procedure drop_UniqueDirections
  end interface
  
  interface size
    module procedure size_UniqueDirections
  end interface
contains

subroutine new_UniqueDirections(this,no_unique_atoms)
  implicit none
  
  type(UniqueDirections), intent(out) :: this
  integer,                intent(in)  :: no_unique_atoms
  
  allocate(this%unique_atoms(no_unique_atoms))
  allocate(this%xy_symmetry(no_unique_atoms))
  allocate(this%xz_symmetry(no_unique_atoms))
  allocate(this%yz_symmetry(no_unique_atoms))
end subroutine

subroutine drop_UniqueDirections(this)
  implicit none
  
  type(UniqueDirections), intent(inout) :: this
  
  deallocate(this%unique_atoms)
  deallocate(this%xy_symmetry)
  deallocate(this%xz_symmetry)
  deallocate(this%yz_symmetry)
end subroutine

function size_UniqueDirections(this) result(output)
  implicit none
  
  type(UniqueDirections), intent(in) :: this
  integer                            :: output
  
  output = size(this%unique_atoms)
end function

function read_unique_directions_file(filename) result(this)
  implicit none
  
  type(String), intent(in) :: filename
  type(UniqueDirections)   :: this
  
  type(String), allocatable :: contents(:)
  type(String), allocatable :: line(:)
  
  integer :: i
  
  contents = read_lines(filename)
  
  call new(this,size(contents))
  
  do i=1,size(contents)
    line = split(contents(i))
    this%unique_atoms(i) = int(line(1))
    this%xy_symmetry(i) = int(line(2))
    this%xz_symmetry(i) = int(line(3))
    this%yz_symmetry(i) = int(line(4))
  enddo
end function

subroutine write_unique_directions_file(this,filename)
  implicit none
  
  type(UniqueDirections), intent(in) :: this
  type(String),           intent(in) :: filename
  
  integer      :: unique_directions_file
  type(String) :: line
  
  integer :: i
  
  unique_directions_file = open_write_file(filename)
  do i=1,size(this)
    line = this%unique_atoms(i) //' '// &
         & this%xy_symmetry(i)  //' '// &
         & this%xz_symmetry(i)  //' '// &
         & this%yz_symmetry(i)
    call print_line(unique_directions_file,line)
  enddo
  close(unique_directions_file)
end subroutine

function calculate_unique_directions(structure,symmetry_group) result(this)
  use structure_module
  use group_module
  implicit none
  
  ! Inputs
  type(StructureData),              intent(in) :: structure
  type(Group),         allocatable, intent(in) :: symmetry_group(:)
  type(UniqueDirections)                       :: this
  
  ! Parameters
  ! Maximum dot product between two unit vectors before they are no longer
  !    considered linearly independent. (30 degree separation).
  real(dp), parameter :: max_dot = dsqrt(3.0_dp)/2 + 1.0e-5_dp
  
  ! Rotations in cartesian co-ordinates.
  real(dp), allocatable :: rotations_cart(:,:,:)
  
  ! Unique atom variables.
  integer              :: no_unique_atoms
  integer, allocatable :: unique_atoms(:)
  
  ! Symmetry vector variables.
  logical  :: vec1_found
  real(dp) :: vec1(3)
  real(dp) :: vec2(3)
  
  ! Temporary variables
  integer :: i,j,k
  
  ! Identify a minimal set of atoms from which the others can be constructed
  !    by symmetry operations.
  allocate(unique_atoms(structure%no_atoms))
  unique_atoms = -1 ! To make it obvious if an uninitialised value is used.
  no_unique_atoms = 0
  do_i : do i=1,structure%no_atoms
    do j=1,size(symmetry_group)
      do k=1,no_unique_atoms
        if (operate(symmetry_group(j),i)==unique_atoms(k)) then
          cycle do_i
        endif
      enddo
    enddo
    no_unique_atoms = no_unique_atoms + 1
    unique_atoms(no_unique_atoms) = i
  enddo do_i
  
  call new(this,no_unique_atoms)
  this%unique_atoms = unique_atoms(1:no_unique_atoms)
  
  rotations_cart = calculate_cartesian_rotations(structure)
  
  ! Identify which directions are related by symmetry, and record the operators
  !    which relate them
  this%xy_symmetry = 0
  this%xz_symmetry = 0
  this%yz_symmetry = 0
  do i=1,no_unique_atoms
    vec1_found = .false.
    do j=1,size(symmetry_group)
      ! Ignore symmetries which do not map this atom to itself.
      if (operate(symmetry_group(j),i)/=i) then
        cycle
      endif
      
      ! Check it the symmetry maps (1,0,0) onto a linearly independent vector.
      if (abs(rotations_cart(1,1,j)) < max_dot) then
        if (.not. vec1_found) then
          ! Only one lin. indep. vector has been found so far. Check if it
          !    lies closer to (0,1,0) or (0,0,1).
          vec1 = rotations_cart(:,1,j)
          vec1_found = .true.
          if (abs(vec1(2)) > abs(vec1(3))) then
            this%xy_symmetry(i) = j
          else
            this%xz_symmetry(i) = j
          endif
          ! At this point, y->z symmetries are ignored, since if x->y and y->z
          !    symmetries exist, then so do x->z symmetries, and these make
          !    later maths easier.
          this%yz_symmetry(i) = 0
        else
          ! Another lin. indep. vector has already been found. Check if the
          !    new vector is lin. indep. to the old vector.
          vec2 = rotations_cart(:,1,j)
          if (abs(dot_product(vec1,vec2)) < max_dot) then
            ! vec2 is lin. indep. Note the symmetry operation.
            ! n.b. at this stage, which is xy and which is xz is irrelevant.
            if (this%xy_symmetry(i)==0) then
              this%xy_symmetry(i) = j
            else
              this%xz_symmetry(i) = j
            endif
            cycle
          endif
        endif
      endif
      
      ! Check if the symmetry maps (0,1,0) onto a linearly independent vector.
      ! n.b. This will be ignored if symmetries from (1,0,0) are found.
      if (.not. vec1_found) then
        if (abs(rotations_cart(2,2,j)) < max_dot) then
          this%yz_symmetry(i) = j
        endif
      endif
    enddo
  enddo
end function
end module
