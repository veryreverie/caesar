! ======================================================================
! Determines which directions are related by symmetry,
!    and thus which atoms should be displaced in which directions 
!    in order to construct the matrix of force constants.
! ======================================================================

! The symmetry operations to consider are nicely outlined here:
! http://www.homepages.ucl.ac.uk/~ucfbdxa/phon/node4.html

module unique_directions_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  ! All the directions in which atoms need to be perturbed in order to map out
  !    the harmonic Born-Oppenheimer surface.
  type UniqueDirections
    ! Atom id.
    integer,      allocatable :: atoms(:)
    ! Direction, in 1/2/3 form and 'x'/'y'/'z' form. (1='x' etc.)
    integer,      allocatable :: directions_int(:)
    character(1), allocatable :: directions_char(:)
    ! Mode id. ( = 3*(atom id - 1) + direction_int ).
    integer,      allocatable :: modes(:)
  end type
  
  interface UniqueDirections
    module procedure new_UniqueDirections
  end interface
  
  interface size
    module procedure size_UniqueDirections
  end interface
contains

function new_UniqueDirections(no_unique_directions) result(this)
  implicit none
  
  integer, intent(in)    :: no_unique_directions
  type(UniqueDirections) :: this
  
  integer :: ialloc
  
  allocate( this%atoms(no_unique_directions),           &
          & this%directions_int(no_unique_directions),  &
          & this%directions_char(no_unique_directions), &
          & this%modes(no_unique_directions),           &
          & stat=ialloc); call err(ialloc)
end function

function size_UniqueDirections(this) result(output)
  implicit none
  
  type(UniqueDirections), intent(in) :: this
  integer                            :: output
  
  output = size(this%atoms)
end function

! ----------------------------------------------------------------------
! Converts (1,2,3) to ('x','y','z')
! ----------------------------------------------------------------------
elemental function direction_int_to_char(direction_int) result(direction_char)
  implicit none
  
  integer, intent(in) :: direction_int
  character(1)        :: direction_char
  
  if (direction_int==1) then
    direction_char = 'x'
  elseif (direction_int==2) then
    direction_char = 'y'
  elseif (direction_int==3) then
    direction_char = 'z'
  endif
end function

! ----------------------------------------------------------------------
! Converts ('x','y','z') to (1,2,3)
! ----------------------------------------------------------------------
elemental function direction_char_to_int(direction_char) result(direction_int)
  implicit none
  
  character(1), intent(in) :: direction_char
  integer                  :: direction_int
  
  if (direction_char=='x') then
    direction_int = 1
  elseif (direction_char=='y') then
    direction_int = 2
  elseif (direction_char=='z') then
    direction_int = 3
  endif
end function

! ----------------------------------------------------------------------
! Reads a UniqueDirections from a file.
! ----------------------------------------------------------------------
function read_unique_directions_file(filename) result(this)
  use ifile_module
  implicit none
  
  type(String), intent(in) :: filename
  type(UniqueDirections)   :: this
  
  type(IFile) :: unique_directions_file
  
  type(String), allocatable :: line(:)
  integer                   :: i
  
  unique_directions_file = filename
  
  this = UniqueDirections(size(unique_directions_file)-1)
  
  do i=1,size(unique_directions_file)-1
    line = split(unique_directions_file%line(i+1))
    this%atoms(i) = int(line(2))
    this%directions_char(i) = line(4)
    this%directions_int(i) = direction_char_to_int(this%directions_char(i))
    this%modes(i) = (this%atoms(i)-1)*3 + this%directions_int(i)
  enddo
end function

! ----------------------------------------------------------------------
! Writes a UniqueDirections to a file.
! ----------------------------------------------------------------------
subroutine write_unique_directions_file(this,filename)
  use ofile_module
  implicit none
  
  type(UniqueDirections), intent(in) :: this
  type(String),           intent(in) :: filename
  
  type(OFile) :: unique_directions_file
  
  integer :: i
  
  unique_directions_file = filename
  call unique_directions_file%print_line('Atoms to be perturbed in order to &
     &map out the harmonic Born-Oppenheimer surface.')
  do i=1,size(this)
    call unique_directions_file%print_line('atom: '//this%atoms(i)// &
       & ' direction: '//this%directions_char(i))
  enddo
end subroutine

! ----------------------------------------------------------------------
! Works out which atoms need to be perturbed in which directions
!    in order to map out the harmonic Born-Oppenheimer surface.
! ----------------------------------------------------------------------
function calculate_unique_directions(structure,atom_symmetry_group) &
   & result(this)
  use constants_module, only : pi
  use structure_module
  use group_module
  use linear_algebra_module
  implicit none
  
  ! Inputs
  type(StructureData),              intent(in) :: structure
  type(Group),         allocatable, intent(in) :: atom_symmetry_group(:)
  type(UniqueDirections)                       :: this
  
  ! A parameter to determine whether or not two vectors are independent.
  ! All rotations will be at most six-fold, so a.b will be at most cos(2*pi/6) 
  !    if a and b are independent.
  ! The parameter is set at half way between cos(pi/3) and 1, to account for
  !    numerical errors.
  real(dp), parameter :: max_dot_product = (1+cos(pi/3))/2.0_dp
  
  ! Whether or not directions are linearly independent.
  logical, allocatable :: unique_dirs(:,:)
  
  ! Symmetry id.
  integer :: previous_symmetry
  
  ! Unique atom variables.
  integer              :: no_unique_atoms
  integer, allocatable :: unique_atoms(:)
  
  ! Cartesian rotations.
  type(RealMatrix), allocatable :: cartesian_rotations(:)
  real(dp),         allocatable :: rotations_cart(:,:,:)
  
  ! Temporary variables
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Identify a minimal set of atoms from which the others can be constructed
  !    by symmetry operations.
  ! --------------------------------------------------
  allocate(unique_atoms(structure%no_atoms), stat=ialloc); call err(ialloc)
  no_unique_atoms = 0
  do_i : do i=1,structure%no_atoms
    do j=1,size(atom_symmetry_group)
      do k=1,no_unique_atoms
        if (atom_symmetry_group(j)*i == unique_atoms(k)) then
          cycle do_i
        endif
      enddo
    enddo
    no_unique_atoms = no_unique_atoms + 1
    unique_atoms(no_unique_atoms) = i
  enddo do_i
  
  ! --------------------------------------------------
  ! Identify which directions (in cartesian co-ordinates) are
  !    related by symmetry.
  ! --------------------------------------------------
  cartesian_rotations = structure%calculate_cartesian_rotations()
  allocate( rotations_cart(3,3,size(structure%symmetries)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(cartesian_rotations)
    rotations_cart(:,:,i) = dble(cartesian_rotations(i))
  enddo
  
  allocate(unique_dirs(3,no_unique_atoms), stat=ialloc); call err(ialloc)
  unique_dirs = .true.
  do i=1,no_unique_atoms
    previous_symmetry = 0
    do j=1,size(atom_symmetry_group)
      
      ! Ignore symmetries which do not map this atom onto itself.
      if (atom_symmetry_group(j) * unique_atoms(i) /= unique_atoms(i)) then
        cycle
      endif
      
      ! Ignore symmetries only diagonal elements.
      if (all( abs([ rotations_cart(1,1,j),  &
         &           rotations_cart(2,2,j),  &
         &           rotations_cart(3,3,j)]) &
         & > max_dot_product)) then
        cycle
      endif
      
      if (previous_symmetry==0) then
        ! This is the first symmetry found.
        previous_symmetry = j
        
        ! Find which direction is mapped onto.
        if ( abs(rotations_cart(2,1,j)) > abs(rotations_cart(3,1,j)) .and. &
           & abs(rotations_cart(2,1,j)) > abs(rotations_cart(3,2,j)) ) then
          unique_dirs(2,i) = .false.
        else
          unique_dirs(3,i) = .false.
        endif
      else
        ! This is not the first symmetry found.
        ! Check if the two symmetries produce independent vectors.
        if ( abs(dot_product( rotations_cart(1,:,j),                  &
           &                  rotations_cart(1,:,previous_symmetry))) &
           & < max_dot_product) then
          unique_dirs(:,i) = [.true., .false., .false.]
          exit
        endif
      endif
    enddo
  enddo
  
  ! --------------------------------------------------
  ! Create output.
  ! --------------------------------------------------
  this = UniqueDirections(count(unique_dirs))
  
  k = 1
  do i=1,no_unique_atoms
    do j=1,3
      if (.not. unique_dirs(j,i)) then
        cycle
      endif
      
      this%atoms(k) = unique_atoms(i)
      this%directions_int(k) = j
      this%directions_char(k) = direction_int_to_char(j)
      this%modes(k) = (i-1)*3+j
      
      k = k+1
    enddo
  enddo
end function
end module
