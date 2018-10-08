! ======================================================================
! A minimal representation of the SymmetryOperator class.
! ======================================================================
module basic_symmetry_module
  use utils_module
  
  use atom_module
  use spglib_module
  implicit none
  
  private
  
  public :: BasicSymmetry
  
  type, extends(Stringsable) :: BasicSymmetry
    integer                      :: id
    type(IntMatrix)              :: tensor
    type(RealVector)             :: translation
    type(Group)                  :: atom_group
    type(IntVector), allocatable :: rvectors(:)
  contains
    procedure, public :: read  => read_BasicSymmetry
    procedure, public :: write => write_BasicSymmetry
  end type
  
  interface BasicSymmetry
    module procedure new_BasicSymmetry
    module procedure new_BasicSymmetry_SpglibSymmetries
    module procedure new_BasicSymmetry_Strings
    module procedure new_BasicSymmetry_StringArray
  end interface
contains

! Constructor.
function new_BasicSymmetry(id,tensor,translation,atom_group,rvectors) &
   & result(this)
  implicit none
  
  integer,          intent(in) :: id
  type(IntMatrix),  intent(in) :: tensor
  type(RealVector), intent(in) :: translation
  type(Group),      intent(in) :: atom_group
  type(IntVector),  intent(in) :: rvectors(:)
  type(BasicSymmetry)          :: this
  
  this%id          = id
  this%tensor      = tensor
  this%translation = translation
  this%atom_group  = atom_group
  this%rvectors    = rvectors
end function

! ----------------------------------------------------------------------
! Calculates symmetry data from Spglib symmetries.
! ----------------------------------------------------------------------
function new_BasicSymmetry_SpglibSymmetries(input,atoms,symmetry_precision) &
   & result(output)
  implicit none
  
  type(SpglibSymmetries), intent(in) :: input
  type(AtomData),         intent(in) :: atoms(:)
  real(dp),               intent(in) :: symmetry_precision
  type(BasicSymmetry), allocatable   :: output(:)
  
  ! Variables for calculating atom groups and R-vectors.
  type(Group),      allocatable :: atom_groups(:)
  type(IntVector),  allocatable :: rvectors(:,:)
  type(RealVector)              :: transformed_position
  type(RealVector), allocatable :: kj_displacements(:)
  type(IntVector),  allocatable :: kj_rvectors(:)
  real(dp),         allocatable :: kj_distances(:)
  integer,          allocatable :: atom_group(:)
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Calculate atom group and R-vectors.
  ! --------------------------------------------------
  ! The symmetry S transforms atomic equilibrium positions, {r_i} in
  !    fractional co-ordinates as:
  ! r_i -> r_j + R
  !
  ! atom_group * i = j
  ! rvectors(i) = R
  allocate( atom_groups(size(input)),           &
          & rvectors(size(atoms), size(input)), &
          & kj_displacements(size(atoms)),      &
          & kj_rvectors(size(atoms)),           &
          & kj_distances(size(atoms)),          &
          & atom_group(size(atoms)),            &
          & stat=ialloc); call err(ialloc)
  do i=1,size(input)
    do j=1,size(atoms)
      ! Calculate the position of the transformed atom.
      transformed_position = input%tensors(i)               &
                         & * atoms(j)%fractional_position() &
                         & + input%translations(i)
      
      do k=1,size(atoms)
        ! Calculate the displacement from atom k to the transfomed atom j.
        kj_displacements(k) = transformed_position &
                          & - atoms(k)%fractional_position()
        ! Calculate the nearest R-vector to the displacement.
        kj_rvectors(k) = nint(kj_displacements(k))
        ! Calculate the distance between the displacement and the nearest
        !    R-vector.
        kj_distances(k) = l2_norm(kj_displacements(k)-kj_rvectors(k))
      enddo
      
      ! Check that exactly one k->j distance is small.
      if (count(kj_distances<symmetry_precision)==0) then
        call print_line(ERROR//': A symmetry maps an atom onto space containing no atoms. Try increasing &
                       &symmetry_precision.')
        call err()
      elseif (count(kj_distances<symmetry_precision)>1) then
        call print_line(ERROR//': A symmetry maps an atom onto more than one atom. Try reducing &
                       &symmetry_precision.')
        call err()
      endif
      
      ! Identify the atom which symmetry i maps atom j to.
      k = minloc(kj_distances,1)
      
      if (atoms(k)%species()/=atoms(j)%species()) then
        call print_line(ERROR//': symmetry operation maps between atoms of &
           &different species.')
        call err()
      endif
      
      atom_group(j) = atoms(k)%id()
      rvectors(j,i) = kj_rvectors(k)
    enddo
    atom_groups(i) = Group(atom_group)
  enddo
  
  ! --------------------------------------------------
  ! Construct output.
  ! --------------------------------------------------
  
  allocate(output(size(input)), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    output(i) = BasicSymmetry( id          = i,                     &
                             & tensor      = input%tensors(i),      &
                             & translation = input%translations(i), &
                             & atom_group  = atom_groups(i),        &
                             & rvectors    = rvectors(:,i)          )
  enddo
end function

! I/O.
subroutine read_BasicSymmetry(this,input)
  implicit none
  
  class(BasicSymmetry), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  integer                      :: id
  type(IntMatrix)              :: tensor
  type(RealVector)             :: translation
  integer,         allocatable :: atom_group(:)
  type(IntVector), allocatable :: rvectors(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(BasicSymmetry)
    line = split_line(input(1))
    id = int(line(2))
    tensor = IntMatrix(input(3:5))
    translation = RealVector(input(7))
    
    allocate( atom_group(size(input)-8), &
            & rvectors(size(input)-8),   &
            & stat=ialloc); call err(ialloc)
    do i=1,size(atom_group)
      line = split_line(input(8+i))
      atom_group(i) = int(line(5))
      rvectors(i) = vec(int(line(8:10)))
    enddo
    
    this = BasicSymmetry(id,tensor,translation,Group(atom_group),rvectors)
  class default
    call err()
  end select
end subroutine

function write_BasicSymmetry(this) result(output)
  implicit none
  
  class(BasicSymmetry), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  type(String), allocatable :: atom_strings(:)
  
  integer :: i,ialloc
  
  allocate( atom_strings(size(this%rvectors)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this%rvectors)
    atom_strings(i) = 'Atom '//i//' -> Atom '//this%atom_group*i// &
                    & ' + R-vector '//this%rvectors(i)
  enddo
  
  select type(this); type is(BasicSymmetry)
    output = [ 'Operation '//this%id,                          &
             & str('Tensor in fractional co-ordinates:'),      &
             & str(this%tensor),                               &
             & str('Translation in fractional co-ordinates:'), &
             & str(this%translation),                          &
             & str('Effect on atoms:'),                        &
             & atom_strings                                    ]
  class default
    call err()
  end select
end function

function new_BasicSymmetry_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(BasicSymmetry)      :: this
  
  call this%read(input)
end function

impure elemental function new_BasicSymmetry_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(BasicSymmetry)           :: this
  
  this = BasicSymmetry(str(input))
end function
end module
