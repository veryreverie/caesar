! ======================================================================
! Defines the StressData type and the StressPointer.
! This type extends StressBase with methods for constructing
!    and interpolating the stress.
! ======================================================================
module caesar_stress_data_module
  use caesar_common_module
  
  use caesar_subspace_coupling_module
  use caesar_anharmonic_data_module
  use caesar_subspace_state_module
  use caesar_subspace_braket_module
  use caesar_basis_state_module
  use caesar_basis_states_module
  use caesar_abstract_classes_module
  implicit none
  
  public :: StressData
  public :: StressPointer
  
  type, abstract, extends(StressBase) :: StressData
  contains
    procedure(representation_StressData), public, deferred, nopass :: &
       & representation
    procedure, public :: startup => startup_StressData
    
    ! Interpolation of the stress.
    procedure, public :: can_be_interpolated => &
                       & can_be_interpolated_StressData
    procedure, public :: calculate_dynamical_matrices => &
                       & calculate_dynamical_matrices_StressData
    procedure, public :: stress_correction => &
                       & stress_correction_StressData
  end type
  
  type, extends(StressData) :: StressPointer
    type(String),                   private :: representation_
    class(StressData), allocatable, private :: stress_
  contains
    procedure, private :: check => check_StressPointer
    
    procedure, public, nopass :: representation => &
                               & representation_StressPointer
    
    procedure, public :: zero_stress => zero_stress_StressPointer
    procedure, public :: add_constant => add_constant_StressPointer
    
    procedure, public :: stress_RealModeDisplacement => &
                       & stress_RealModeDisplacement_StressPointer
    procedure, public :: stress_ComplexModeDisplacement => &
                       & stress_ComplexModeDisplacement_StressPointer
    
    procedure, public :: braket_SubspaceBraKet  => &
                       & braket_SubspaceBraKet_StressPointer
    procedure, public :: braket_BasisState  => &
                       & braket_BasisState_StressPointer
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_StressPointer
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_StressPointer
    
    procedure, public :: potential_stress_SubspaceBraKet => &
                       & potential_stress_SubspaceBraKet_StressPointer
    procedure, public :: potential_stress_BasisState => &
                       & potential_stress_BasisState_StressPointer
    
    procedure, public :: can_be_interpolated => &
                       & can_be_interpolated_StressPointer
    procedure, public :: calculate_dynamical_matrices => &
                       & calculate_dynamical_matrices_StressPointer
    procedure, public :: stress_correction => &
                       & stress_correction_StressPointer
    
    ! I/O.
    procedure, public :: read  => read_StressPointer
    procedure, public :: write => write_StressPointer
  end type
  
  ! An array of all types which extend StressData.
  ! This array will be filled in by startup routines.
  type(StressPointer), allocatable :: TYPES_StressData(:)
  
  interface StressPointer
    module procedure new_StressPointer
    module procedure new_StressPointer_Strings
    module procedure new_StressPointer_StringArray
  end interface
  
  abstract interface
    ! StressData procedures.
    impure elemental function representation_StressData() result(output)
      import String
      implicit none
      
      type(String) :: output
    end function
  end interface
contains

subroutine startup_StressData(this)
  implicit none
  
  class(StressData), intent(in) :: this
  
  type(StressBasePointer) :: base
  
  integer :: i
  
  ! TODO: this doesn't work. It calls the startup() method for
  !    StressBasePointer only.
  base = StressBasePointer(this)
  call base%startup()
  
  if (.not.allocated(TYPES_StressData)) then
    TYPES_StressData = [StressPointer(this)]
  elseif (.not.any([(                                      &
     & this%representation()                               &
     &    == TYPES_StressData(i)%stress_%representation(), &
     & i=1,                                                &
     & size(TYPES_StressData)                              )])) then
    TYPES_StressData = [TYPES_StressData, StressPointer(this)]
  endif
end subroutine

! Construct a StressPointer from any type which extends StressData.
impure elemental function new_StressPointer(stress) result(this)
  implicit none
  
  class(StressData), intent(in) :: stress
  type(StressPointer)           :: this
  
  integer :: ialloc
  
  select type(stress); type is(StressPointer)
    this = stress
  class default
    this%representation_ = stress%representation()
    allocate( this%stress_, source=stress, &
            & stat=ialloc); call err(ialloc)
  end select
end function

! Checks that the pointer has been allocated before it is used.
impure elemental subroutine check_StressPointer(this)
  implicit none
  
  class(StressPointer), intent(in) :: this
  
  if (.not. allocated(this%stress_)) then
    call print_line(CODE_ERROR//': Trying to use a StressPointer before &
       &it has been allocated.')
    call err()
  endif
end subroutine

! Type representation.
impure elemental function representation_StressPointer() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'pointer'
end function

impure elemental subroutine zero_stress_StressPointer(this)
  implicit none
  
  class(StressPointer), intent(inout) :: this
  
  call this%check()
  
  call this%stress_%zero_stress()
end subroutine

impure elemental subroutine add_constant_StressPointer(this,input)
  implicit none
  
  class(StressPointer), intent(inout) :: this
  type(RealMatrix),     intent(in)    :: input
  
  call this%check()
  
  call this%stress_%add_constant(input)
end subroutine

impure elemental function stress_RealModeDisplacement_StressPointer(this, &
   & displacement) result(output)
  implicit none
  
  class(StressPointer),       intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  type(RealMatrix)                       :: output
  
  call this%check()
  
  output = this%stress_%stress(displacement)
end function

impure elemental function stress_ComplexModeDisplacement_StressPointer( &
   & this,displacement) result(output)
  implicit none
  
  class(StressPointer),          intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexMatrix)                       :: output
  
  call this%check()
  
  output = this%stress_%stress(displacement)
end function

impure elemental subroutine braket_SubspaceBraKet_StressPointer(this,braket, &
   & whole_subspace,anharmonic_data) 
  implicit none
  
  class(StressPointer),  intent(inout)        :: this
  class(SubspaceBraKet), intent(in)           :: braket
  logical,               intent(in), optional :: whole_subspace
  type(AnharmonicData),  intent(in)           :: anharmonic_data
  
  call this%check()
  
  call this%stress_%braket(braket,whole_subspace,anharmonic_data)
end subroutine

impure elemental subroutine braket_BasisState_StressPointer(this,bra,ket, &
   & subspace,subspace_basis,whole_subspace,anharmonic_data) 
  implicit none
  
  class(StressPointer),     intent(inout)        :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  logical,                  intent(in), optional :: whole_subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  call this%check()
  
  call this%stress_%braket( bra,            &
                          & ket,            &
                          & subspace,       &
                          & subspace_basis, &
                          & whole_subspace, &
                          & anharmonic_data )
end subroutine

impure elemental subroutine braket_BasisStates_StressPointer(this,states, &
   & subspace,subspace_basis,whole_subspace,anharmonic_data) 
  implicit none
  
  class(StressPointer),     intent(inout)        :: this
  class(BasisStates),       intent(inout)        :: states
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  logical,                  intent(in), optional :: whole_subspace
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  call this%check()
  
  call this%stress_%braket( states,         &
                          & subspace,       &
                          & subspace_basis, &
                          & whole_subspace, &
                          & anharmonic_data )
end subroutine

impure elemental function harmonic_expectation_StressPointer(this,frequency, &
   & thermal_energy,supercell_size,anharmonic_data) result(output)
  implicit none
  
  class(StressPointer), intent(in) :: this
  real(dp),             intent(in) :: frequency
  real(dp),             intent(in) :: thermal_energy
  integer,              intent(in) :: supercell_size
  type(AnharmonicData), intent(in) :: anharmonic_data
  type(RealMatrix)                 :: output
  
  call this%check()
  
  output = this%stress_%harmonic_expectation( frequency,      &
                                            & thermal_energy, &
                                            & supercell_size, &
                                            & anharmonic_data )
end function

recursive function potential_stress_SubspaceBraKet_StressPointer(this,braket, &
   & anharmonic_data) result(output) 
  implicit none
  
  class(StressPointer),  intent(in) :: this
  class(SubspaceBraKet), intent(in) :: braket
  type(AnharmonicData),  intent(in) :: anharmonic_data
  type(RealMatrix)                  :: output
  
  call this%check()
  
  output = this%stress_%potential_stress(braket, anharmonic_data)
end function

recursive function potential_stress_BasisState_StressPointer(this,bra,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output) 
  implicit none
  
  class(StressPointer),     intent(in)           :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(RealMatrix)                               :: output
  
  call this%check()
  
  output = this%stress_%potential_stress( bra,            &
                                        & ket,            &
                                        & subspace,       &
                                        & subspace_basis, &
                                        & anharmonic_data )
end function

function can_be_interpolated_StressPointer(this) result(output)
  implicit none
  
  class(StressPointer), intent(in) :: this
  logical                          :: output
  
  call this%check()
  
  output = this%stress_%can_be_interpolated()
end function

function calculate_dynamical_matrices_StressPointer(this,qpoints,             &
   & thermal_energy,subspaces,subspace_bases,subspace_states,anharmonic_data) &
   & result(output) 
  implicit none
  
  class(StressPointer),     intent(in)     :: this
  type(QpointData),         intent(in)     :: qpoints(:)
  real(dp),                 intent(in)     :: thermal_energy
  type(DegenerateSubspace), intent(in)     :: subspaces(:)
  class(SubspaceBasis),     intent(in)     :: subspace_bases(:)
  class(BasisStates),       intent(inout)  :: subspace_states(:)
  type(AnharmonicData),     intent(in)     :: anharmonic_data
  type(StressDynamicalMatrix), allocatable :: output(:)
  
  call this%check()
  
  output = this%stress_%calculate_dynamical_matrices( qpoints,             &
                                                    & thermal_energy,      &
                                                    & subspaces,           &
                                                    & subspace_bases,      &
                                                    & subspace_states,     &
                                                    & anharmonic_data      )
end function

function stress_correction_StressPointer(this,subspaces,subspace_bases, &
   & subspace_states,anharmonic_data) result(output) 
  implicit none
  
  class(StressPointer),     intent(in)    :: this
  type(DegenerateSubspace), intent(in)    :: subspaces(:)
  class(SubspaceBasis),     intent(in)    :: subspace_bases(:)
  class(BasisStates),       intent(inout) :: subspace_states(:)
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  type(RealMatrix)                        :: output
  
  call this%check()
  
  output = this%stress_%stress_correction( subspaces,       &
                                         & subspace_bases,  &
                                         & subspace_states, &
                                         & anharmonic_data  )
end function

function can_be_interpolated_StressData(this) result(output)
  implicit none
  
  class(StressData), intent(in) :: this
  logical                       :: output
  
  output = .false.
end function

function calculate_dynamical_matrices_StressData(this,qpoints,                &
   & thermal_energy,subspaces,subspace_bases,subspace_states,anharmonic_data) &
   & result(output) 
  implicit none
  
  class(StressData),        intent(in)     :: this
  type(QpointData),         intent(in)     :: qpoints(:)
  real(dp),                 intent(in)     :: thermal_energy
  type(DegenerateSubspace), intent(in)     :: subspaces(:)
  class(SubspaceBasis),     intent(in)     :: subspace_bases(:)
  class(BasisStates),       intent(inout)  :: subspace_states(:)
  type(AnharmonicData),     intent(in)     :: anharmonic_data
  type(StressDynamicalMatrix), allocatable :: output(:)
  
  ! This should be gated behind can_be_interpolated.
  call print_line(CODE_ERROR//': calculate_dynamical_matrices not &
     &implemented for this stress.')
  call err()
end function

function stress_correction_StressData(this,subspaces,subspace_bases, &
   & subspace_states,anharmonic_data) result(output) 
  implicit none
  
  class(StressData),        intent(in)    :: this
  type(DegenerateSubspace), intent(in)    :: subspaces(:)
  class(SubspaceBasis),     intent(in)    :: subspace_bases(:)
  class(BasisStates),       intent(inout) :: subspace_states(:)
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  type(RealMatrix)                        :: output
  
  ! This should be gated behind can_be_interpolated.
  call print_line(CODE_ERROR//': stress_correction not implemented for this &
     &potential.')
  call err()
end function

! I/O.
subroutine read_StressPointer(this,input)
  implicit none
  
  class(StressPointer), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(StressPointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    i = first([(                                                       &
       & TYPES_StressData(i)%stress_%representation()==representation, &
       & i=1,                                                          &
       & size(TYPES_StressData)                                        )])
    
    ! Read the input into the element of the correct type,
    !    and copy that element into the output.
    call TYPES_StressData(i)%stress_%read(input(2:))
    this = StressPointer(TYPES_StressData(i))
  class default
    call err()
  end select
end subroutine

function write_StressPointer(this) result(output)
  implicit none
  
  class(StressPointer), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  select type(this); type is(StressPointer)
    output = [ 'Stress representation: '//this%representation_, &
             & str(this%stress_)                                ]
  end select
end function

function new_StressPointer_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(StressPointer)      :: this
  
  call this%read(input)
end function

impure elemental function new_StressPointer_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(StressPointer)           :: this
  
  this = StressPointer(str(input))
end function
end module
