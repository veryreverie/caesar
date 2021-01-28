! ======================================================================
! A state along each complex mode in a degenerate subspace,
!    at a set of q-points for which 2q/=G.
! ======================================================================
module caesar_harmonic_state_complex_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_harmonic_state_2d_module
  implicit none
  
  private
  
  public :: startup_harmonic_state_complex
  
  public :: HarmonicStateComplex
  
  public :: harmonic_state_complex_pointer
  
  public :: prod_complex
  
  type, extends(SubspaceState) :: HarmonicStateComplex
    type(HarmonicState2D), allocatable :: modes_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_HarmonicStateComplex
    
    procedure, public :: mode_ids => &
                       & mode_ids_HarmonicStateComplex
    procedure, public :: paired_mode_ids => &
                       & paired_mode_ids_HarmonicStateComplex
    
    procedure, public :: occupation => occupation_HarmonicStateComplex
    
    procedure, public :: change_modes => change_modes_HarmonicStateComplex
    
    procedure, public :: wavevector => wavevector_HarmonicStateComplex
    
    procedure, public :: wavefunction => wavefunction_HarmonicStateComplex
    
    ! I/O.
    procedure, public :: read  => read_HarmonicStateComplex
    procedure, public :: write => write_HarmonicStateComplex
  end type
  
  interface HarmonicStateComplex
    module procedure new_HarmonicStateComplex
    module procedure new_HarmonicStateComplex_SubspaceState
    module procedure new_HarmonicStateComplex_Strings
    module procedure new_HarmonicStateComplex_StringArray
  end interface
contains

impure elemental function prod_complex(lhs,rhs) result(output)
  implicit none
  
  type(HarmonicStateComplex), intent(in) :: lhs
  type(HarmonicStateComplex), intent(in) :: rhs
  type(HarmonicStateComplex)             :: output
  
  output = HarmonicStateComplex([lhs%modes_,rhs%modes_])
end function

! Startup procedure.
subroutine startup_harmonic_state_complex()
  implicit none
  
  type(HarmonicStateComplex) :: state
  
  call state%startup()
end subroutine

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function new_HarmonicStateComplex(modes) result(this) 
  implicit none
  
  type(HarmonicState2D), intent(in) :: modes(:)
  type(HarmonicStateComplex)        :: this
  
  this%modes_ = modes
end function

recursive function new_HarmonicStateComplex_SubspaceState(input) result(this)
  implicit none
  
  class(SubspaceState), intent(in) :: input
  type(HarmonicStateComplex)       :: this
  
  select type(input); type is(HarmonicStateComplex)
    this = input
  type is(SubspaceStatePointer)
    ! WORKAROUND: ifort doesn't recognise the interface to this function
    !    from within this function, so the full name is used instead.
    this = new_HarmonicStateComplex_SubspaceState(input%state())
  class default
    call err()
  end select
end function

! Cast a class(SubspaceState) to a pointer of type(HarmonicStateComplex).
! N.B. this must only be called on inputs with the TARGET attribute.
recursive function harmonic_state_complex_pointer(input) result(this)
  implicit none
  
  class(SubspaceState), intent(in), target :: input
  type(HarmonicStateComplex), pointer      :: this
  
  select type(input); type is(HarmonicStateComplex)
    this => input
  type is(SubspaceStatePointer)
    this => harmonic_state_complex_pointer(input%state_pointer())
  class default
    call err()
  end select
end function

! ----------------------------------------------------------------------
! Type representation.
! ----------------------------------------------------------------------
impure elemental function representation_HarmonicStateComplex() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'harmonic complex'
end function

! ----------------------------------------------------------------------
! Returns the modes spanned by the state.
! ----------------------------------------------------------------------
function mode_ids_HarmonicStateComplex(this) result(output)
  implicit none
  
  class(HarmonicStateComplex), intent(in) :: this
  integer, allocatable                    :: output(:)
  
  output = this%modes_%id()
end function

function paired_mode_ids_HarmonicStateComplex(this) result(output)
  implicit none
  
  class(HarmonicStateComplex), intent(in) :: this
  integer, allocatable                    :: output(:)
  
  output = this%modes_%paired_id()
end function

! ----------------------------------------------------------------------
! Returns the total occupation of a given state.
! ----------------------------------------------------------------------
! The total occupation of the state product_{q,i}|n_{q,i}> is equal to
!    sum_{q,i} n_{q,i}.
impure elemental function occupation_HarmonicStateComplex(this) &
   & result(output)
  implicit none
  
  class(HarmonicStateComplex), intent(in) :: this
  integer                                 :: output
  
  output = sum(this%modes_%total_occupation())
end function

! ----------------------------------------------------------------------
! Returns the wavevector of a given state.
! ----------------------------------------------------------------------
! The wavevector of the state product_{q,i} |n_{q,i}> is equal to
!    sum_{q,i} n_{q,i}q.
function wavevector_HarmonicStateComplex(this,modes,qpoints) result(output)
  implicit none
  
  class(HarmonicStateComplex), intent(in) :: this
  type(ComplexMode),           intent(in) :: modes(:)
  type(QpointData),            intent(in) :: qpoints(:)
  type(FractionVector)                    :: output
  
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
impure elemental function wavefunction_HarmonicStateComplex(this,frequency, &
   & supercell) result(output)
  implicit none
  
  class(HarmonicStateComplex), intent(in) :: this
  real(dp),                    intent(in) :: frequency
  type(StructureData),         intent(in) :: supercell
  type(String)                            :: output
  
  ! TODO
  call err()
end function

! ----------------------------------------------------------------------
! Change the modes of the state by the specified group.
! ----------------------------------------------------------------------
impure elemental function change_modes_HarmonicStateComplex(this,mode_group) &
   & result(output)
  implicit none
  
  class(HarmonicStateComplex), intent(in) :: this
  type(Group),                 intent(in) :: mode_group
  type(HarmonicStateComplex)              :: output
  
  integer, allocatable :: ids(:)
  integer, allocatable :: paired_ids(:)
  integer, allocatable :: occupations(:)
  integer, allocatable :: paired_occupations(:)
  integer, allocatable :: sort_key(:)
  
  ! Get the ids and occupations of the single-mode terms.
  ids = this%modes_%id()
  paired_ids = this%modes_%paired_id()
  occupations = this%modes_%occupation()
  paired_occupations = this%modes_%paired_occupation()
  
  ! Change the ids according to the given group.
  ids = mode_group*ids
  
  ! Sort the modes by id.
  sort_key = sort(ids)
  ids = ids(sort_key)
  paired_ids = paired_ids(sort_key)
  occupations = occupations(sort_key)
  paired_occupations = paired_occupations(sort_key)
  
  ! Construct output using the new ids.
  output = HarmonicStateComplex(                                         &
     & modes = HarmonicState2D( id                = ids,                 &
     &                          paired_id         = paired_ids,          &
     &                          occupation        = occupations,         &
     &                          paired_occupation = paired_occupations ) )
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_HarmonicStateComplex(this,input)
  implicit none
  
  class(HarmonicStateComplex), intent(out) :: this
  type(String),                intent(in)  :: input(:)
  
  type(HarmonicState2D), allocatable :: modes(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i
  
  select type(this); type is(HarmonicStateComplex)
    line = split_line(input(1),delimiter='>')
    line = [(line(i)//'>',i=1,size(line))]
    modes = HarmonicState2D(line)
    
    this = HarmonicStateComplex(modes)
  class default
    call err()
  end select
end subroutine

function write_HarmonicStateComplex(this) result(output)
  implicit none
  
  class(HarmonicStateComplex), intent(in) :: this
  type(String), allocatable               :: output(:)
  
  select type(this); type is(HarmonicStateComplex)
    output = [join(str(this%modes_), delimiter='')]
  class default
    call err()
  end select
end function

function new_HarmonicStateComplex_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)   :: input(:)
  type(HarmonicStateComplex) :: this
  
  call this%read(input)
end function

impure elemental function new_HarmonicStateComplex_StringArray(input) &
   & result(this) 
  implicit none
  
  type(StringArray), intent(in) :: input
  type(HarmonicStateComplex)    :: this
  
  this = HarmonicStateComplex(str(input))
end function
end module
