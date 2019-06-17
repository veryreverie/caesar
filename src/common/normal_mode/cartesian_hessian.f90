! ======================================================================
! The second derivatives of the potential, in cartesian co-ordinates.
! ======================================================================
module cartesian_hessian_module
  use utils_module
  
  use structure_module
  implicit none
  
  private
  
  public :: CartesianHessian
  
  ! The Hessian is the matrix of force constants, F, such that F.x=f,
  !    where x and f are the cartesian displacement and cartesian force
  !    respectively.
  type, extends(Stringsable) :: CartesianHessian
    type(RealMatrix), allocatable, private :: elements_(:,:)
  contains
    procedure, public :: elements
    
    ! I/O.
    procedure, public :: read  => read_CartesianHessian
    procedure, public :: write => write_CartesianHessian
  end type
  
  interface CartesianHessian
    module procedure new_CartesianHessian
    module procedure new_CartesianHessian_elements
    module procedure new_CartesianHessian_Strings
    module procedure new_CartesianHessian_StringArray
  end interface
contains

! Constructor.
function new_CartesianHessian(elements) result(this)
  implicit none
  
  type(RealMatrix), intent(in) :: elements(:,:)
  type(CartesianHessian)       :: this
  
  this%elements_ = elements
end function

! ----------------------------------------------------------------------
! Constructs and checks Hessian from given matrices.
! ----------------------------------------------------------------------
function new_CartesianHessian_elements(structure,elements,logfile) &
   & result(this)
  implicit none
  
  type(StructureData), intent(in)              :: structure
  type(RealMatrix),    intent(in)              :: elements(:,:)
  type(OFile),         intent(inout), optional :: logfile
  type(CartesianHessian)                       :: this
  
  ! Variables for checking Hessian.
  type(AtomData)   :: atom_i
  type(AtomData)   :: atom_j
  type(AtomData)   :: atom_k
  type(AtomData)   :: atom_i2
  type(AtomData)   :: atom_j2
  integer          :: rvec_i
  integer          :: rvec_j
  integer          :: rvec_k
  integer          :: rvec_ij
  type(RealMatrix) :: matrix
  type(RealMatrix) :: symmetric
  type(RealMatrix) :: copy
  real(dp)         :: average
  real(dp)         :: difference
  
  integer :: i,j,k
  
  ! Copy elements into CartesianHessian.
  this = CartesianHessian(elements)
  
  ! Check that the force constants between atom i and atom j are the same as
  !    those between atom i + R and atom j + R.
  average = 0.0_dp
  difference = 0.0_dp
  do i=1,structure%no_atoms
    atom_i = structure%atoms(i)
    rvec_i = atom_i%rvec_id()
    do j=1,structure%no_atoms
      atom_j = structure%atoms(j)
      rvec_j = atom_j%rvec_id()
      rvec_ij = structure%paired_rvector_group(rvec_i) * rvec_j
      
      do k=1,structure%no_atoms
        atom_k = structure%atoms(k)
        rvec_k = atom_k%rvec_id()
        if (atom_k%prim_id()/=atom_j%prim_id()) then
          cycle
        elseif (rvec_k/=rvec_ij) then
          cycle
        endif
        matrix = this%elements_(atom_i%id(),atom_j%id())
        copy = this%elements_(atom_i%prim_id(),atom_k%id())
        average = average + sum_squares((matrix+copy)/2)
        difference = difference + sum_squares(matrix-copy)
      enddo
    enddo
  enddo
  
  if (present(logfile)) then
    call logfile%print_line('Fractional L2 difference in force constants at &
       &different R-vectors: '//sqrt(difference/average))
  endif
  
  if (sqrt(difference/average)>1.0e-10_dp) then
    call print_line(WARNING//': Hessian does not obey R-vector symmetries. &
       & Please check log files.')
    call err()
  endif
  
  ! Check F(i1,i2) transforms correctly under symmetry operators.
  average = 0
  difference = 0
  do i=1,size(structure%symmetries)
    do j=1,structure%no_atoms_prim
      atom_i = structure%atoms(j)
      atom_i2 = structure%atoms( structure%symmetries(i)%atom_group &
                             & * atom_i%id())
      do k=1,structure%no_atoms
        atom_j = structure%atoms(k)
        atom_j2 = structure%atoms( structure%symmetries(i)%atom_group &
                               & * atom_j%id())
        matrix = this%elements(atom_i2,atom_j2)
        symmetric = structure%symmetries(i)%cartesian_tensor &
                & * this%elements(atom_i,atom_j)             &
                & * transpose(structure%symmetries(i)%cartesian_tensor)
        average = average + sum_squares((matrix+symmetric)/2)
        difference = difference + sum_squares(matrix-symmetric)
      enddo
    enddo
  enddo
  call logfile%print_line(                                  &
     & 'Fractional L2 error in symmetry of F(i1,i2)   : '// &
     & sqrt(difference/average))
  if (sqrt(difference/average) > 1e-10_dp) then
    call print_line(WARNING//': F(i1,i2) is not as symmetric as expected. &
       &Please check log files.')
    call print_line('Fractional L2 error: '//sqrt(difference/average))
  endif
  
  ! Check F(i1,i2)=F(i2,i1).
  average = 0
  difference = 0
  do i=1,structure%no_atoms
    atom_i = structure%atoms(i)
    do j=1,structure%no_atoms
      atom_j = structure%atoms(j)
      matrix = this%elements(atom_j,atom_i)
      symmetric = transpose(this%elements(atom_i,atom_j))
      average = average + sum_squares((matrix+symmetric)/2)
      difference = difference + sum_squares(matrix-symmetric)
    enddo
  enddo
  call logfile%print_line(                                  &
     & 'Fractional L2 error in F(i1,i2)=F(i2,i1)      : '// &
     & sqrt(difference/average))
  if (sqrt(difference/average) > 1e-10_dp) then
    call print_line(WARNING//': F(i1,i2)/=F(i2,i1). Please check log files.')
  endif
end function

! ----------------------------------------------------------------------
! Return the force constants between two atoms.
! ----------------------------------------------------------------------
function elements(this,a,b) result(output)
  implicit none
  
  class(CartesianHessian), intent(in) :: this
  type(AtomData),          intent(in) :: a
  type(AtomData),          intent(in) :: b
  type(RealMatrix)                    :: output
  
  output = this%elements_(a%id(),b%id())
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_CartesianHessian(this,input)
  implicit none
  
  class(CartesianHessian), intent(out) :: this
  type(String),            intent(in)  :: input(:)
  
  type(StringArray), allocatable :: elements(:)
  integer                        :: no_atoms
  
  integer :: i,j,k,ialloc
  
  select type(this); type is(CartesianHessian)
    elements = split_into_sections(input)
    no_atoms = int_sqrt(size(elements))
    allocate(this%elements_(no_atoms,no_atoms), stat=ialloc); call err(ialloc)
    k = 0
    do i=1,no_atoms
      do j=1,no_atoms
        k = k+1
        this%elements_(j,i) = RealMatrix(elements(k)%strings(2:4))
      enddo
    enddo
  class default
    call err()
  end select
end subroutine

function write_CartesianHessian(this) result(output)
  implicit none
  
  class(CartesianHessian), intent(in) :: this
  type(String), allocatable           :: output(:)
  
  integer :: no_atoms
  
  type(String) :: matrix_strings(3)
  
  integer :: i,j,k,ialloc
  
  select type(this); type is(CartesianHessian)
    no_atoms = size(this%elements_,1)
    if (size(this%elements_,2)/=no_atoms) then
      call err()
    endif
    
    allocate(output(5*no_atoms*no_atoms), stat=ialloc); call err(ialloc)
    k = 0
    do i=1,no_atoms
      do j=1,no_atoms
        k = k+1
        matrix_strings = str(this%elements_(j,i))
        output(5*k-4) = 'Atoms: ('//j//' '//i//')'
        output(5*k-3) = matrix_strings(1)
        output(5*k-2) = matrix_strings(2)
        output(5*k-1) = matrix_strings(3)
        output(5*k)   = ''
      enddo
    enddo
  class default
    call err()
  end select
end function

function new_CartesianHessian_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(CartesianHessian)   :: this
  
  call this%read(input)
end function

impure elemental function new_CartesianHessian_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(CartesianHessian)        :: this
  
  this = CartesianHessian(str(input))
end function
end module
