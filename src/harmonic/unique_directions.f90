! This program determines which directions are related by symmetry,
!    and thus which atoms should be displaced in which directions 
!    to construct the matrix of force constants.

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
    ! Mode id. ( = 3*(atom id - 1) + direction ).
    integer,      allocatable :: modes(:)
  end type
  
  interface new
    module procedure new_UniqueDirections
  end interface
  
  interface size
    module procedure size_UniqueDirections
  end interface
contains

subroutine new_UniqueDirections(this,no_unique_directions)
  implicit none
  
  type(UniqueDirections), intent(out) :: this
  integer,                intent(in)  :: no_unique_directions
  
  integer :: ialloc
  
  allocate( this%atoms(no_unique_directions),           &
          & this%directions_int(no_unique_directions),  &
          & this%directions_char(no_unique_directions), &
          & this%modes(no_unique_directions),           &
          & stat=ialloc); call err(ialloc)
end subroutine

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
  implicit none
  
  type(String), intent(in) :: filename
  type(UniqueDirections)   :: this
  
  type(String), allocatable :: unique_directions_file(:)
  type(String), allocatable :: line(:)
  
  integer :: i
  
  unique_directions_file = read_lines(filename)
  
  call new(this,size(unique_directions_file)-1)
  
  do i=1,size(unique_directions_file)-1
    line = split(unique_directions_file(i+1))
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
  implicit none
  
  type(UniqueDirections), intent(in) :: this
  type(String),           intent(in) :: filename
  
  integer      :: unique_directions_file
  
  integer :: i
  
  unique_directions_file = open_write_file(filename)
  call print_line(unique_directions_file, 'Atoms to be perturbed in order to &
     &map out the harmonic Born-Oppenheimer surface.')
  do i=1,size(this)
    call print_line(unique_directions_file, 'atom: '//this%atoms(i)// &
       & ' direction: '//this%directions_char(i))
  enddo
  close(unique_directions_file)
end subroutine

! ----------------------------------------------------------------------
! Works out which atoms need to be perturbed in which directions
!    in order to map out the harmonic Born-Oppenheimer surface.
! ----------------------------------------------------------------------
function calculate_unique_directions(structure,symmetry_group) result(this)
  use constants_module, only : identity
  use structure_module
  use group_module
  implicit none
  
  ! Inputs
  type(StructureData),              intent(in) :: structure
  type(Group),         allocatable, intent(in) :: symmetry_group(:)
  type(UniqueDirections)                       :: this
  
  ! Whether or not directions are linearly independent.
  logical, allocatable :: unique_dirs_frac(:,:) ! Fractional co-ordinates.
  logical, allocatable :: unique_dirs_cart(:,:) ! Cartesian co-ordinates.
  
  ! Symmetry id.
  integer :: previous_symmetry
  
  ! Unique atom variables.
  integer              :: no_unique_atoms
  integer, allocatable :: unique_atoms(:)
  
  ! Temporary variables
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Identify a minimal set of atoms from which the others can be constructed
  !    by symmetry operations.
  ! --------------------------------------------------
  allocate(unique_atoms(structure%no_atoms), stat=ialloc); call err(ialloc)
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
  
  ! --------------------------------------------------
  ! Identify which cardinal directions (in fractional co-ordinates)
  !    on the minimal set of atoms are related by symmetry.
  ! --------------------------------------------------
  
  ! Identify related directions in fractional co-ordinates.
  allocate(unique_dirs_frac(3,no_unique_atoms), stat=ialloc); call err(ialloc)
  unique_dirs_frac = .true.
  do i=1,no_unique_atoms
    previous_symmetry = 0
    do j=1,size(symmetry_group)
      ! Ignore symmetries which do not map this atom onto itself.
      if (operate(symmetry_group(j), unique_atoms(i)) /= unique_atoms(i)) then
        cycle
      endif
      
      ! Ignore symmetries which only map directions onto themselves.
      if (all(abs(structure%rotations(:,:,j))==identity)) then
        cycle
      endif
      
      if (previous_symmetry==0) then
        ! This is the first symmetry found.
        previous_symmetry = j
        
        ! Find which direction is mapped onto.
        if ( abs(structure%rotations(2,1,j))>abs(structure%rotations(3,1,j))  &
           & .or.                                                             &
           & abs(structure%rotations(2,1,j))>abs(structure%rotations(3,2,j))) &
           & then
          unique_dirs_frac(2,i) = .false.
        else
          unique_dirs_frac(3,i) = .false.
        endif
      else
        ! This is not the first symmetry found.
        ! Check if the symmetries produce independent vectors.
        if (abs(dot_product( structure%rotations(1,:,j),                 &
             &               structure%rotations(1,:,previous_symmetry)) &
             & ) /= 1) then
          unique_dirs_frac(:,i) = [ .true., .false., .true. ]
          exit
        endif
      endif
    enddo
  enddo
  
  ! --------------------------------------------------
  ! Identify the cardinal directions in cartesian co-ordinates which
  !    most closely match those in fractional co-ordinates.
  ! --------------------------------------------------
  allocate(unique_dirs_cart(3,no_unique_atoms), stat=ialloc); call err(ialloc)
  do i=1,no_unique_atoms
    unique_dirs_cart(:,i) = [ .true.,.true.,.true. ]
    if (all(unique_dirs_frac(:,i) .eqv. .true.)) then
      ! All directions are independent in both co-ordinate systems.
      continue
    elseif (all(unique_dirs_frac(:,i) .eqv. [.true.,.false.,.false.])) then
      ! All directions are dependent in both co-ordinate systems.
      unique_dirs_cart(2:3,i) = .false.
    else
      ! There are two independent directions.
      if (unique_dirs_frac(2,i) .eqv. .false.) then
        if ( abs(structure%lattice(2,1)) > abs(structure%lattice(2,2)) .and. &
           & abs(structure%lattice(2,1)) > abs(structure%lattice(2,3))) then
          unique_dirs_cart(1,i) = .false.
        elseif (abs(structure%lattice(2,2)) > abs(structure%lattice(2,3))) then
          unique_dirs_cart(2,i) = .false.
        else
          unique_dirs_cart(3,i) = .false.
        endif
      else
        if ( abs(structure%lattice(3,1)) > abs(structure%lattice(3,2)) .and. &
           & abs(structure%lattice(3,1)) > abs(structure%lattice(3,3))) then
          unique_dirs_cart(1,i) = .false.
        elseif (abs(structure%lattice(3,2)) > abs(structure%lattice(3,3))) then
          unique_dirs_cart(2,i) = .false.
        else
          unique_dirs_cart(3,i) = .false.
        endif
      endif
    endif
  enddo
  
  ! --------------------------------------------------
  ! Create output.
  ! --------------------------------------------------
  call new(this,count(unique_dirs_cart))
  
  k = 1
  do i=1,no_unique_atoms
    do j=1,3
      if (.not. unique_dirs_cart(j,i)) then
        cycle
      endif
      
      this%atoms(k) = unique_atoms(i)
      this%directions_int(k) = j
      this%directions_char(k) = direction_int_to_char(j)
      this%modes(k) = (i-1)*3+j
      
      k = k+1
    enddo
  enddo
!  
!  
!  
!  
!  
!  
!  
!  rotations_cart = calculate_cartesian_rotations(structure)
!  
!  ! Identify which directions are related by symmetry, and record the operators
!  !    which relate them
!  this%xy_symmetry = 0
!  this%xz_symmetry = 0
!  this%yz_symmetry = 0
!  do i=1,no_unique_atoms
!    vec1_found = .false.
!    do j=1,size(symmetry_group)
!      ! Ignore symmetries which do not map this atom to itself.
!      if (operate(symmetry_group(j),i)/=i) then
!        cycle
!      endif
!      
!      ! Check it the symmetry maps (1,0,0) onto a linearly independent vector.
!      if (abs(rotations_cart(1,1,j)) < max_dot) then
!        if (.not. vec1_found) then
!          ! Only one lin. indep. vector has been found so far. Check if it
!          !    lies closer to (0,1,0) or (0,0,1).
!          vec1 = rotations_cart(:,1,j)
!          vec1_found = .true.
!          if (abs(vec1(2)) > abs(vec1(3))) then
!            this%xy_symmetry(i) = j
!          else
!            this%xz_symmetry(i) = j
!          endif
!          ! At this point, y->z symmetries are ignored, since if x->y and y->z
!          !    symmetries exist, then so do x->z symmetries, and these make
!          !    later maths easier.
!          this%yz_symmetry(i) = 0
!        else
!          ! Another lin. indep. vector has already been found. Check if the
!          !    new vector is lin. indep. to the old vector.
!          vec2 = rotations_cart(:,1,j)
!          if (abs(dot_product(vec1,vec2)) < max_dot) then
!            ! vec2 is lin. indep. Note the symmetry operation.
!            ! n.b. at this stage, which is xy and which is xz is irrelevant.
!            if (this%xy_symmetry(i)==0) then
!              this%xy_symmetry(i) = j
!            else
!              this%xz_symmetry(i) = j
!            endif
!            cycle
!          endif
!        endif
!      endif
!      
!      ! Check if the symmetry maps (0,1,0) onto a linearly independent vector.
!      ! n.b. This will be ignored if symmetries from (1,0,0) are found.
!      if (.not. vec1_found) then
!        if (abs(rotations_cart(2,2,j)) < max_dot) then
!          this%yz_symmetry(i) = j
!        endif
!      endif
!    enddo
!  enddo
end function
end module
