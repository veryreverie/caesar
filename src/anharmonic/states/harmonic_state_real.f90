! ======================================================================
! A state along each complex mode in a degenerate subspace,
!    at a set of q-points for which 2q=G.
! ======================================================================
module harmonic_state_real_module
  use common_module
  
  use anharmonic_common_module
  
  use harmonic_state_1d_module
  implicit none
  
  private
  
  public :: startup_harmonic_state_real
  
  public :: HarmonicStateReal
  
  public :: harmonic_state_real_pointer
  
  public :: prod_real
  
  type, extends(SubspaceState) :: HarmonicStateReal
    type(HarmonicState1D), allocatable :: modes_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_HarmonicStateReal
    
    procedure, public :: mode_ids => mode_ids_HarmonicStateReal
    procedure, public :: paired_mode_ids => paired_mode_ids_HarmonicStateReal
    
    procedure, public :: change_modes => change_modes_HarmonicStateReal
    
    procedure, public :: occupation => occupation_HarmonicStateReal
    
    procedure, public :: wavevector => wavevector_HarmonicStateReal
    
    procedure, public :: wavefunction => wavefunction_HarmonicStateReal
    
    ! I/O.
    procedure, public :: read  => read_HarmonicStateReal
    procedure, public :: write => write_HarmonicStateReal
  end type
  
  interface HarmonicStateReal
    module procedure new_HarmonicStateReal
    module procedure new_HarmonicStateReal_SubspaceState
    module procedure new_HarmonicStateReal_Strings
    module procedure new_HarmonicStateReal_StringArray
  end interface
contains

impure elemental function prod_real(lhs,rhs) result(output)
  implicit none
  
  type(HarmonicStateReal), intent(in) :: lhs
  type(HarmonicStateReal), intent(in) :: rhs
  type(HarmonicStateReal)             :: output
  
  output = HarmonicStateReal([lhs%modes_,rhs%modes_])
end function

! Startup procedure.
subroutine startup_harmonic_state_real()
  implicit none
  
  type(HarmonicStateReal) :: state
  
  call state%startup()
end subroutine

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function new_HarmonicStateReal(modes) result(this) 
  implicit none
  
  type(HarmonicState1D), intent(in) :: modes(:)
  type(HarmonicStateReal)           :: this
  
  this%modes_ = modes
end function

recursive function new_HarmonicStateReal_SubspaceState(input) result(this)
  implicit none
  
  class(SubspaceState), intent(in) :: input
  type(HarmonicStateReal)          :: this
  
  select type(input); type is(HarmonicStateReal)
    this = input
  type is(SubspaceStatePointer)
    ! WORKAROUND: ifort doesn't recognise the interface to this function
    !    from within this function, so the full name is used instead.
    this = new_HarmonicStateReal_SubspaceState(input%state())
  class default
    call err()
  end select
end function

! Cast a class(SubspaceState) to a pointer of type(HarmonicStateReal).
! N.B. this must only be called on inputs with the TARGET attribute.
recursive function harmonic_state_real_pointer(input) result(this)
  implicit none
  
  class(SubspaceState), intent(in), target :: input
  type(HarmonicStateReal), pointer         :: this
  
  select type(input); type is(HarmonicStateReal)
    this => input
  type is(SubspaceStatePointer)
    this => harmonic_state_real_pointer(input%state_pointer())
  class default
    call err()
  end select
end function

! ----------------------------------------------------------------------
! Type representation.
! ----------------------------------------------------------------------
impure elemental function representation_HarmonicStateReal() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'harmonic real'
end function

! ----------------------------------------------------------------------
! Returns the modes spanned by the state.
! ----------------------------------------------------------------------
function mode_ids_HarmonicStateReal(this) result(output)
  implicit none
  
  class(HarmonicStateReal), intent(in) :: this
  integer, allocatable                 :: output(:)
  
  output = this%modes_%id()
end function

function paired_mode_ids_HarmonicStateReal(this) result(output)
  implicit none
  
  class(HarmonicStateReal), intent(in) :: this
  integer, allocatable                 :: output(:)
  
  output = this%modes_%id()
end function

! ----------------------------------------------------------------------
! Returns the total occupation of a given state.
! ----------------------------------------------------------------------
! The total occupation of the state product_{q,i} |n_{q,i}> is equal to
!    sum_{q,i} n_{q,i}.
impure elemental function occupation_HarmonicStateReal(this) result(output)
  implicit none
  
  class(HarmonicStateReal), intent(in) :: this
  integer                              :: output
  
  output = sum(this%modes_%total_occupation())
end function

! ----------------------------------------------------------------------
! Returns the wavevector of a given state.
! ----------------------------------------------------------------------
! The wavevector of the state product_{q,i} |n_{q,i}> is equal to
!    sum_{q,i} n_{q,i}q.
function wavevector_HarmonicStateReal(this,modes,qpoints) result(output)
  implicit none
  
  class(HarmonicStateReal), intent(in) :: this
  type(ComplexMode),        intent(in) :: modes(:)
  type(QpointData),         intent(in) :: qpoints(:)
  type(FractionVector)                 :: output
  
  integer :: i
  
  ! Effectively sum(this%modes_(:)%wavevector(modes,qpoints)).
  output = this%modes_(1)%wavevector(modes,qpoints)
  do i=2,size(this%modes_)
    output = output + this%modes_(i)%wavevector(modes,qpoints)
  enddo
end function

! ----------------------------------------------------------------------
! Returns the wavefunction of the state,
!    with all coefficients accounted for.
! ----------------------------------------------------------------------
impure elemental function wavefunction_HarmonicStateReal(this,frequency, &
   & supercell) result(output)
  implicit none
  
  class(HarmonicStateReal), intent(in) :: this
  real(dp),                 intent(in) :: frequency
  type(StructureData),      intent(in) :: supercell
  type(String)                         :: output
  
  ! TODO
  call err()
end function

! ----------------------------------------------------------------------
! Change the modes of the state by the specified group.
! ----------------------------------------------------------------------
impure elemental function change_modes_HarmonicStateReal(this,mode_group) &
   & result(output)
  implicit none
  
  class(HarmonicStateReal), intent(in) :: this
  type(Group),              intent(in) :: mode_group
  type(HarmonicStateReal)              :: output
  
  integer, allocatable :: ids(:)
  integer, allocatable :: occupations(:)
  integer, allocatable :: sort_key(:)
  
  ! Get the ids and occupations of the single-mode terms.
  ids = this%modes_%id()
  occupations = this%modes_%occupation()
  
  ! Change the ids according to the given group.
  ids = mode_group*ids
  
  ! Sort the modes by id.
  sort_key = sort(ids)
  ids = ids(sort_key)
  occupations = occupations(sort_key)
  
  ! Construct output using the new ids.
  output = HarmonicStateReal(modes = HarmonicState1D(ids,occupations))
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_HarmonicStateReal(this,input)
  implicit none
  
  class(HarmonicStateReal), intent(out) :: this
  type(String),             intent(in)  :: input(:)
  
  type(HarmonicState1D), allocatable :: modes(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i
  
  select type(this); type is(HarmonicStateReal)
    line = split_line(input(1),delimiter='>')
    line = [(line(i)//'>',i=1,size(line))]
    modes = HarmonicState1D(line)
    
    this = HarmonicStateReal(modes)
  class default
    call err()
  end select
end subroutine

function write_HarmonicStateReal(this) result(output)
  implicit none
  
  class(HarmonicStateReal), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  select type(this); type is(HarmonicStateReal)
    output = [join(str(this%modes_), delimiter='')]
  class default
    call err()
  end select
end function

function new_HarmonicStateReal_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(HarmonicStateReal)  :: this
  
  call this%read(input)
end function

impure elemental function new_HarmonicStateReal_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(HarmonicStateReal)       :: this
  
  this = HarmonicStateReal(str(input))
end function
end module
