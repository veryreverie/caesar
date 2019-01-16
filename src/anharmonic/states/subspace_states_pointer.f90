! ======================================================================
! A wrapped polymorphic pointer to a set of subspace states.
! ======================================================================
! Wraps all of SubspaceStates's methods,
!    calling them on the pointed-to states.
module subspace_states_pointer_module
  use common_module
  
  use anharmonic_common_module
  
  use full_subspace_basis_and_states_module
  implicit none
  
  private
  
  public :: SubspaceStatesPointer
  
  type, extends(SubspaceStates) :: SubspaceStatesPointer
    type(String),                       private :: representation_
    class(SubspaceStates), allocatable, private :: states_
  contains
    procedure, public :: check => check_SubspaceStatesPointer
    
    procedure, public :: states => states_SubspaceStatesPointer
    procedure, public :: spectra => spectra_SubspaceStatesPointer
    procedure, public :: wavefunctions => wavefunctions_SubspaceStatesPointer
    procedure, public :: integrate => integrate_SubspaceStatesPointer
  end type
  
  interface SubspaceStatesPointer
    module procedure new_SubspaceStatesPointer
  end interface
contains

! Construct a SubspaceStatesPointer from any type which extends SubspaceStates.
impure elemental function new_SubspaceStatesPointer(states) result(this)
  implicit none
  
  class(SubspaceStates), intent(in) :: states
  type(SubspaceStatesPointer)       :: this
  
  integer :: ialloc
  
  select type(states); type is(SubspaceStatesPointer)
    this = states
  type is(FullSubspaceStates)
    this%representation_ = 'full'
    allocate( this%states_, source=states, &
            & stat=ialloc); call err(ialloc)
  class default
    call err()
  end select
end function

! Checks that the pointer has been allocated before it is used.
subroutine check_SubspaceStatesPointer(this)
  implicit none
  
  class(SubspaceStatesPointer), intent(in) :: this
  
  if (.not. allocated(this%states_)) then
    call print_line(CODE_ERROR//': Trying to use a SubspaceStatesPointer &
       &before it has been allocated.')
    call err()
  endif
end subroutine

! ----------------------------------------------------------------------
! SubspaceStates methods.
! ----------------------------------------------------------------------
function states_SubspaceStatesPointer(this,subspace,subspace_basis, &
   & anharmonic_data) result(output)
  implicit none
  
  class(SubspaceStatesPointer), intent(in) :: this
  type(DegenerateSubspace),     intent(in) :: subspace
  class(SubspaceBasis),         intent(in) :: subspace_basis
  type(AnharmonicData),         intent(in) :: anharmonic_data
  class(SubspaceState), allocatable        :: output(:)
  
  call this%check()
  
  output = this%states_%states(subspace, subspace_basis, anharmonic_data)
end function

impure elemental function spectra_SubspaceStatesPointer(this,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceStatesPointer), intent(in) :: this
  type(DegenerateSubspace),     intent(in) :: subspace
  class(SubspaceBasis),         intent(in) :: subspace_basis
  type(AnharmonicData),         intent(in) :: anharmonic_data
  type(EnergySpectra)                      :: output
  
  call this%check()
  
  output = this%states_%spectra(subspace, subspace_basis, anharmonic_data)
end function

impure elemental function wavefunctions_SubspaceStatesPointer(this,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceStatesPointer), intent(in)  :: this
  type(DegenerateSubspace),     intent(in)  :: subspace
  class(SubspaceBasis),         intent(in)  :: subspace_basis
  type(AnharmonicData),         intent(in)  :: anharmonic_data
  class(SubspaceWavefunctions), allocatable :: output
  
  call this%check()
  
  output = this%states_%wavefunctions( subspace,       &
                                     & subspace_basis, &
                                     & anharmonic_data )
end function

impure elemental function integrate_SubspaceStatesPointer(this,potential, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceStatesPointer), intent(in) :: this
  class(PotentialData),         intent(in) :: potential
  type(DegenerateSubspace),     intent(in) :: subspace
  class(SubspaceBasis),         intent(in) :: subspace_basis
  type(AnharmonicData),         intent(in) :: anharmonic_data
  class(PotentialData), allocatable        :: output
  
  call this%check()
  
  output = this%states_%integrate( potential,      &
                                 & subspace,       &
                                 & subspace_basis, &
                                 & anharmonic_data )
end function
end module
