! ======================================================================
! The SubspaceState abstract class, defining states which do not need to
!    reference a basis.
! ======================================================================
module subspace_state_module
  use common_module
  
  use stress_prefactors_module
  use anharmonic_data_module
  use sparse_monomial_module
  implicit none
  
  private
  
  public :: SubspaceState
  public :: SubspaceStatePointer
  
  type, abstract, extends(Stringsable) :: SubspaceState
    integer :: subspace_id
  contains
    procedure(representation_SubspaceState), public, deferred, nopass :: &
       & representation
    procedure, public :: startup => startup_SubspaceState
    
    ! Return the id of the modes across which the state is defined.
    procedure(mode_ids_SubspaceState), public, deferred :: mode_ids
    procedure(paired_mode_ids_SubspaceState), public, deferred :: &
       & paired_mode_ids
    
    ! Returns the occupation of the state.
    procedure(occupation_SubspaceState), public, deferred :: occupation
    
    ! Returns the wavevector of the state.
    procedure(wavevector_SubspaceState), public, deferred :: wavevector
    
    ! If ket is not given, <this|this>, otherwise <this|ket>.
    procedure(inner_product_SubspaceState), public, deferred :: inner_product
    
    ! Integrals of the form <i|V|j>
    ! Either <this|V|this> or <this|V|ket>, where V is a ComplexMonomial.
    procedure(integrate_SubspaceState), public, deferred :: integrate
    
    ! Either <this|T|this> or <this|T|ket>, where T is the kinetic energy.
    procedure(kinetic_energy_SubspaceState), public, deferred :: &
       & kinetic_energy
    
    ! Either <this|V|this> or <this|V|ket>, where V is the harmonic potential
    !    energy.
    procedure(harmonic_potential_energy_SubspaceState), public, deferred :: &
       & harmonic_potential_energy
    
    ! Either <this|stress|this> or <this|stress|ket>, where stress is the
    !    kinetic stress.
    procedure(kinetic_stress_SubspaceState), public, deferred :: &
       & kinetic_stress
  end type
  
  type, extends(SubspaceState) :: SubspaceStatePointer
    type(String),                      private :: representation_
    ! N.B. state_ should never be modified.
    ! It is public for performance reasons only.
    class(SubspaceState), allocatable, public :: state_
  contains
    procedure, private :: check => check_SubspaceStatePointer
    
    procedure, public, nopass :: representation => &
                               & representation_SubspaceStatePointer
    
    procedure, public :: state => state_SubspaceStatePointer
    
    procedure, public :: mode_ids => &
                       & mode_ids_SubspaceStatePointer
    procedure, public :: paired_mode_ids => &
                       & paired_mode_ids_SubspaceStatePointer
    
    procedure, public :: occupation => occupation_SubspaceStatePointer
    
    procedure, public :: wavevector => wavevector_SubspaceStatePointer
    
    procedure, public :: inner_product => &
                       & inner_product_SubspaceStatePointer
    
    procedure, public :: integrate => &
                       & integrate_SubspaceStatePointer
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_SubspaceStatePointer
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_SubspaceStatePointer
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_SubspaceStatePointer
    
    ! I/O.
    procedure, public :: read  => read_SubspaceStatePointer
    procedure, public :: write => write_SubspaceStatePointer
  end type
  
  ! An array of all types which extend SubspaceState.
  ! This array will be filled in by startup routines.
  type(SubspaceStatePointer), allocatable :: TYPES_SubspaceState(:)
  
  ! Abstract interface for SubspaceState functionality.
  abstract interface
    impure elemental function representation_SubspaceState() result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
    
    function mode_ids_SubspaceState(this) result(output)
      import SubspaceState
      implicit none
      
      class(SubspaceState), intent(in) :: this
      integer, allocatable             :: output(:)
    end function
    
    function paired_mode_ids_SubspaceState(this) result(output)
      import SubspaceState
      implicit none
      
      class(SubspaceState), intent(in) :: this
      integer, allocatable             :: output(:)
    end function
    
    impure elemental function occupation_SubspaceState(this) result(output)
      import SubspaceState
      implicit none
      
      class(SubspaceState), intent(in) :: this
      integer                          :: output
    end function
    
    function wavevector_SubspaceState(this,modes,qpoints) &
       & result(output)
      import SubspaceState
      import ComplexMode
      import QpointData
      import FractionVector
      implicit none
      
      class(SubspaceState), intent(in) :: this
      type(ComplexMode),    intent(in) :: modes(:)
      type(QpointData),     intent(in) :: qpoints(:)
      type(FractionVector)             :: output
    end function
    
    impure elemental function inner_product_SubspaceState(this,ket, &
       & anharmonic_data) result(output)
      import SubspaceState
      import AnharmonicData
      import dp
      implicit none
      
      class(SubspaceState), intent(in)           :: this
      class(SubspaceState), intent(in), optional :: ket
      type(AnharmonicData), intent(in)           :: anharmonic_data
      real(dp)                                   :: output
    end function
    
    impure elemental function integrate_SubspaceState(this,monomial,ket, &
       & anharmonic_data) result(output)
      import SubspaceState
      import SparseMonomial
      import AnharmonicData
      import dp
      implicit none
      
      class(SubspaceState), intent(in)           :: this
      type(SparseMonomial), intent(in)           :: monomial
      class(SubspaceState), intent(in), optional :: ket
      type(AnharmonicData), intent(in)           :: anharmonic_data
      complex(dp)                                :: output
    end function
    
    impure elemental function kinetic_energy_SubspaceState(this,ket, &
       & anharmonic_data) result(output)
      import SubspaceState
      import AnharmonicData
      import dp
      implicit none
      
      class(SubspaceState), intent(in)           :: this
      class(SubspaceState), intent(in), optional :: ket
      type(AnharmonicData), intent(in)           :: anharmonic_data
      real(dp)                                   :: output
    end function
    
    impure elemental function harmonic_potential_energy_SubspaceState(this, &
       & ket,anharmonic_data) result(output)
      import SubspaceState
      import AnharmonicData
      import dp
      implicit none
      
      class(SubspaceState), intent(in)           :: this
      class(SubspaceState), intent(in), optional :: ket
      type(AnharmonicData), intent(in)           :: anharmonic_data
      real(dp)                                   :: output
    end function
    
    impure elemental function kinetic_stress_SubspaceState(this,ket, &
       & stress_prefactors,anharmonic_data) result(output)
      import SubspaceState
      import StressPrefactors
      import AnharmonicData
      import RealMatrix
      implicit none
      
      class(SubspaceState),   intent(in)           :: this
      class(SubspaceState),   intent(in), optional :: ket
      type(StressPrefactors), intent(in)           :: stress_prefactors
      type(AnharmonicData),   intent(in)           :: anharmonic_data
      type(RealMatrix)                             :: output
    end function
  end interface
  
  interface SubspaceStatePointer
    module procedure new_SubspaceStatePointer
    module procedure new_SubspaceStatePointer_Strings
    module procedure new_SubspaceStatePointer_StringArray
  end interface
contains

! Startup method.
subroutine startup_SubspaceState(this)
  implicit none
  
  class(SubspaceState), intent(in) :: this
  
  integer :: i
  
  if (.not. allocated(TYPES_SubspaceState)) then
    TYPES_SubspaceState = [SubspaceStatePointer(this)]
  elseif (.not. any([(                                    &
     &    this%representation()                           &
     & == TYPES_SubspaceState(i)%state_%representation(), &
     & i=1,                                               &
     & size(TYPES_SubspaceState)                          )])) then
    TYPES_SubspaceState = [TYPES_SubspaceState, SubspaceStatePointer(this)]
  endif
end subroutine

! --------------------------------------------------
! SubspaceStatePointer methods.
! --------------------------------------------------
! Construct a SubspaceStatePointer from any type which extends SubspaceState.
impure elemental function new_SubspaceStatePointer(state) result(this)
  implicit none
  
  class(SubspaceState), intent(in) :: state
  type(SubspaceStatePointer)       :: this
  
  integer :: ialloc
  
  select type(state); type is(SubspaceStatePointer)
    this = state
  class default
    this%representation_ = state%representation()
    allocate( this%state_, source=state, &
            & stat=ialloc); call err(ialloc)
    this%subspace_id = this%state_%subspace_id
  end select
end function

! Checks that the pointer has been allocated before it is used.
subroutine check_SubspaceStatePointer(this)
  implicit none
  
  class(SubspaceStatePointer), intent(in) :: this
  
  if (.not. allocated(this%state_)) then
    call print_line(CODE_ERROR//': Trying to use a SubspaceStatePointer &
       &before it has been allocated.')
    call err()
  endif
end subroutine

! Type representation.
impure elemental function representation_SubspaceStatePointer() &
   & result(output)
  implicit none
  
  type(String) :: output
  
  output = 'pointer'
end function

! SubspaceState methods.
function state_SubspaceStatePointer(this) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in) :: this
  class(SubspaceState), allocatable       :: output
  
  output = this%state_
end function

! --------------------------------------------------
! SubspaceStatePointer wrappers for SubspaceState methods.
! --------------------------------------------------
function mode_ids_SubspaceStatePointer(this) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in) :: this
  integer, allocatable                    :: output(:)
  
  call this%check()
  
  output = this%state_%mode_ids()
end function

function paired_mode_ids_SubspaceStatePointer(this) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in) :: this
  integer, allocatable                    :: output(:)
  
  call this%check()
  
  output = this%state_%paired_mode_ids()
end function

impure elemental function occupation_SubspaceStatePointer(this) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in) :: this
  integer                                 :: output
  
  call this%check()
  
  output = this%state_%occupation()
end function

function wavevector_SubspaceStatePointer(this,modes,qpoints) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in) :: this
  type(ComplexMode),           intent(in) :: modes(:)
  type(QpointData),            intent(in) :: qpoints(:)
  type(FractionVector)                    :: output
  
  call this%check()
  
  output = this%state_%wavevector(modes,qpoints)
end function

impure elemental function inner_product_SubspaceStatePointer(this,ket, &
   & anharmonic_data) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in)           :: this
  class(SubspaceState),        intent(in), optional :: ket
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  real(dp)                                          :: output
  
  call this%check()
  
  output = this%state_%inner_product(ket, anharmonic_data)
end function

impure elemental function integrate_SubspaceStatePointer(this,monomial,ket, &
   & anharmonic_data) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in)           :: this
  type(SparseMonomial),        intent(in)           :: monomial
  class(SubspaceState),        intent(in), optional :: ket
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  complex(dp)                                       :: output
  
  call this%check()
  
  output = this%state_%integrate(monomial, ket, anharmonic_data)
end function

impure elemental function kinetic_energy_SubspaceStatePointer(this,ket, &
   & anharmonic_data) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in)           :: this
  class(SubspaceState),        intent(in), optional :: ket
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  real(dp)                                          :: output
  
  call this%check()
  
  output = this%state_%kinetic_energy(ket, anharmonic_data)
end function

impure elemental function harmonic_potential_energy_SubspaceStatePointer( &
   & this,ket,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in)           :: this
  class(SubspaceState),        intent(in), optional :: ket
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  real(dp)                                          :: output
  
  call this%check()
  
  output = this%state_%harmonic_potential_energy(ket, anharmonic_data)
end function

impure elemental function kinetic_stress_SubspaceStatePointer(this,ket, &
   & stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in)           :: this
  class(SubspaceState),        intent(in), optional :: ket
  type(StressPrefactors),      intent(in)           :: stress_prefactors
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  type(RealMatrix)                                  :: output
  
  call this%check()
  
  output = this%state_%kinetic_stress(ket, stress_prefactors, anharmonic_data)
end function

! --------------------------------------------------
! SubspaceStatePointer I/O.
! --------------------------------------------------
subroutine read_SubspaceStatePointer(this,input)
  implicit none
  
  class(SubspaceStatePointer), intent(out) :: this
  type(String),                intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(SubspaceStatePointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    i = first([(                                                         &
       & TYPES_SubspaceState(i)%state_%representation()==representation, &
       & i=1,                                                            &
       & size(TYPES_SubspaceState)                                       )])
    
    ! Read the input into the element of the correct type,
    !    and copy that element into the output.
    call TYPES_SubspaceState(i)%state_%read(input(2:))
    this = SubspaceStatePointer(TYPES_SubspaceState(i))
  class default
    call err()
  end select
end subroutine

function write_SubspaceStatePointer(this) result(output)
  implicit none
  
  class(SubspaceStatePointer), intent(in) :: this
  type(String), allocatable               :: output(:)
  
  select type(this); type is(SubspaceStatePointer)
    output = [ 'SubspaceState representation: '//this%representation_, &
             & str(this%state_)                                        ]
  end select
end function

function new_SubspaceStatePointer_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)   :: input(:)
  type(SubspaceStatePointer) :: this
  
  call this%read(input)
end function

impure elemental function new_SubspaceStatePointer_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SubspaceStatePointer)    :: this
  
  this = SubspaceStatePointer(str(input))
end function
end module
